/* subtree pruning-regrafting (spr) library
 * Peter Cordes <peter@cordes.ca>, Dalhousie University
 * license: GPLv2 or later
 * some algorithms and code in this file come from
 * David Petry's post on sci.math.num-analysis 91/9/23, and
 * Donald Knuth's wonderful work, The Art of Computer Programming
 */

/* routines for the lcg, mainly setup.
 * Intended for use only by other library code
 */

#define _GNU_SOURCE
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define SPR_PRIVATE
#include "spr.h"

// startstate = UINT_MAX means we've just started
// UINT_MAX as a return value means we've looped.
unsigned int lcg(struct lcg *lcgp) {
	// use long just to avoid overflow on multiply
	unsigned long a = lcgp->a, c = lcgp->c, m = lcgp->m;
	unsigned long old = lcgp->state;

	if (lcgp->startstate == UINT_MAX)
		lcgp->startstate = old;
	else if (old == lcgp->startstate)
		return UINT_MAX;
	// will continue to return this until startstate is reset

	lcgp->state = (a*lcgp->state + c) % m;
	return old;
}



/*********** LCG generator for integer sequence covering 0..n ********/

/********** primality testing ******/
/* prime code originally by David Petry:
 *  posted on sci.math.num-analysis  91/9/23 */

/* primeset is a bitmap of all the odd numbers.
 * After running the sieve of Eratosthenes,
 * the only bits left set will be in prime indices.
 * This code should still work with 64bit ints,
 *  but it will waste half the space.  (AMD64 int=32bits, though)
 */

static unsigned int sieved = 0, maxptest = 0;
static unsigned int *primeset = NULL;

#define is_prime_macro(p)	((primeset[(p)>>6]>>(((p)>>1)&31))&1)
#define exclude_macro(p)	( primeset[(p)>>6] &= ((-1) ^ (1<<(((p)>>1)&31))) )

/* can't easily be grown; easiest to re-sieve from 0 every time it has to grow
 * could change inner loop (over n) to start low*low+m (where m makes it a
 * multiple of i)...  too much work.
 */
static void sieve(unsigned int sqrthigh){
	unsigned int n,i;
	unsigned int high = sqrthigh*sqrthigh;

	memset(primeset, -1, sizeof(*primeset) * (1 + (high>>6)));
//	for (i=0 ; i < (high>>6) ; primeset[i++] = -1 );
	for (i=3 ; i <= sqrthigh ; i += 2){
		if (is_prime_macro(i)){
			for (n=i*i ; n <= high ; n+=2*i){
				exclude_macro(n);
			}
		}
	}
}

static void primesetup(int k){
	int sqrsize = (int)ceilf(sqrtf(k));
	if (sqrsize & 0xffff0000){ // Could be relaxed with 64bit int
		fprintf(stderr, "allspr error: lcg: max prime too big: 0x%x\n", sqrsize);
	}else if (sqrsize < 3)
		sqrsize = 7;

	sqrsize = max(21, sqrsize); // for allspr, just go big the first time

	if (sqrsize > sieved){
		primeset = xrealloc(primeset, sizeof(*primeset) *
				    (((sqrsize*sqrsize)>>6) + 1) );
		sieve(sqrsize);
		sieved = sqrsize;
		maxptest = sieved * sieved;
	}
}

static inline int is_prime(int p){
	if (p > maxptest){
		int old = maxptest;
		primesetup(p*2);
		if (spr_debug>=1)
			fprintf(stderr, "spr debug: sieving more primes. old max %u, new %u\n", old, maxptest);
	}
	return p>1 && (p==2 || (p & is_prime_macro(p)));
}

static inline int next_prime (int i)
{
	if (i<=2) return 2;
	i |= 1;			// make i odd
	while (!is_prime(i)) i+=2;

	return i;
}

/******************** Linear Congruential Generator setup ************/
/* to generate all possible SPRs in a pseudo-random order, we generate
 * all the numbers between 0 and the number of possible SPRs once each
 * without repetition using an LCG of the form: x_n+1 = x_n*a + c mod m.
 */
/*
 * Knuth: TAOCP 3.2.1: ex 2: if a and m are relatively prime,
 * the number X_0 will always appear in the period.  (will return to start?)
 *
 *  3.2.1.2: Theorem A: An LCG will have period m iff:
 *   -  c is relatively prime to m
 *   -  b = a-1 is a multiple of p, for every prime p dividing m
 *   -  b is a multiple of 4, if m is a multiple of 4
 *
 *  3.2.1.3: ex 4:  m = 2**e >= 8  ->  maximum potency when a mod 8 = 5.
 *   small multipliers are to be avoided.
 *
 * Numerical recipies suggests c = a prime close to (1/2 - sqrt(3)/6)*m
 */

/* make up some parameters for an LCG that will have the maximum period
 * equal to the range, so every value is generated once.
 * When maxval doesn't have any repeated prime factors, a = m+1,
 * which is the same as a=1.  It's not exactly random, but it does still
 * mix up which SPRs are done.
 *
 * TODO: take advantage of the fact that maxval = floor(sqrt(maxval))*ceil(sqrt(maxval))
 *  could do that, but then the code would be less general-purpose
 * successfully brute-force tested for maxval=1..1000.
 */
void findlcg(struct lcg *lcg_params, int maxval)
{
	unsigned int a, b, c, m = maxval;
	int i;

	primesetup (maxval+maxval/2);

	if (m<=6){ // will be either 6 or 2.  Just loop in order
		b=0;
		c=1;
	}else{
		int divlimit = m;
		b=1; // b must be a multiple of all of m's prime factors
		if (!(m%2)){
			b=2;
			while (divlimit%2 == 0) divlimit /= 2;
		}
		for (i=3 ; i <= divlimit ; i+=2){
			if (is_prime(i) && m%i == 0){
				b *= i;
				while (divlimit%i == 0) divlimit /= i;
			}
		}

		if (!(m%4)){	// if m is a mult of 4, b must be.
			while (b%4) b *= 2;
		}

		/* make sure a isn't too small */
		while (b<sqrtf(m)) b*=7;

		if (b == m) b=0;  // just give up and avoid overflow

// Numerical Recipies says there is "lore" behind this... :)
// TAOCP says it's useless unless the multiplier sucks (section 3.3.3, eq. 40)
// That would be us.
		c = next_prime(max(5, (0.5 - sqrtf(3)/6.0)*m - 2));
		while (m%c == 0) c = next_prime(c+1);
		// Luckily we don't have to test for c>m, because it doesn't
		// happen with any m<100, and there are enough primes later...
/* I've observed that when a == m, (e.g. a=13, c=13, m=72) you often get
 * two consecutive numbers...  Do something to avoid that if it's a problem */
	}

	a = b+1;

	unsigned long long l = (unsigned long long)a * m;
	if (l > ULONG_MAX){
		fprintf(stderr, "spr: chosen Linear Congruential Generator is bogus\n"
		"   x_n+1 = x_n*%u+%u mod %u\n"
		"   a*m > ULONG_MAX, so it would overflow :(\n", a, c, m);
	}

	lcg_params->a = a;
	lcg_params->c = c;
	lcg_params->m = m;
	lcg_params->startstate = UINT_MAX;
	lcg_params->state = rand() % m;
}



#if 0

/*
int cyclic(int s) {
	int t = ((s & 0x1ff)<<22) + (s >> 9) + s;
	if (t<0) return  (t & 0x7fffffff ) + 1;
	return t;
}
*/

#include <time.h>

struct lcg lcg_params;
#define MAGIC 1234
int main (int argc, char *argv[])
{
	int i, maxlcg;
	long n;
	srand( time(NULL) );
	if (argc>=2) maxlcg = atoi(argv[1]);
	else maxlcg = 132; // 11*12

	findlcg( &lcg_params, maxlcg );
	if( argc == 4 ){
		lcg_params.a = atoi(argv[2]);
		lcg_params.c = atoi(argv[3]);
	}
	printf("using LCG x_n+1 = x_n*%u+%u mod %u\n",
	       lcg_params.a, lcg_params.c, lcg_params.m);

	i = 0; n=0;
	while (42){
		i = lcg(i, &lcg_params);
		printf ("%d\n", i);
		++n;
		if (n > maxlcg*2){ puts("big n"); break; }
		if (i==0){ puts("cycled"); break; }
	}
	printf( "lcg: %ld/%d\n", n, i );
	return 0;
}
#endif

void spr_lcg_staticfree(void)
{
	sieved = 0;
	free( primeset );
	primeset = NULL;
}
