/* subtree pruning-regrafting (spr) library
 * Peter Cordes <peter@cordes.ca>, Dalhousie University
 */

/* setup routines that are called when giving the library a new tree */

#define _GNU_SOURCE
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <assert.h>

#include "spr_private.h"
#include "spr.h"

/* format of sprmap (transposed):
01 0122 012333  src
10 2201 333012  dest
 size 2: 2 2
 size 3: 4 6
 size 4: 6 12
space required: sum( i=2..n, 2*(i-1) ) = n*(n-1).
* TODO: src = 0 entries are useless (if the root is at 0). can't SPR the root!
* Putting entries with src=root at the end could make it possible to vary the length
* of the array that is covered.  A length with repeated prime factors could be
* chosen, so the LCG that generates a sequence covering the whole space will be more
* random
*/

int (*sprmap)[2] = NULL;
int sprmapnodes = 0;

/* TODO: replace sprmap with a function call, if there's a fast formula...
 * TODO: this should probably have some locking for re-entrancy, then the whole
 * library would be re-entrant.
 */

static void grow_map( int nnodes )
{
	int i, n;	      // the number we're currently filling in
	int newsize = nnodes*(nnodes-1);

	if ((n=sprmapnodes) < 2) n=2; // start where we left off last time

	if (nnodes < 2 || nnodes <= sprmapnodes) return; 
	sprmap = xrealloc( sprmap, newsize * sizeof(*sprmap) ); // NULL is ok.
	sprmapnodes = nnodes;	// TODO: real locking
 
	for( ; n<=nnodes ; n++ ){
		for( i=0 ; i<n-1 ; i++ ){
			sprmap[(n-1)*(n-2) + i][0] = i;
			sprmap[(n-1)*(n-2) + i][1] = n-1;

			sprmap[(n-1)*(n-2) + i + n-1][0] = n-1;
			sprmap[(n-1)*(n-2) + i + n-1][1] = i;
			// == (n-1)*(n-1) + i
		}
	}
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

unsigned int sieved = 0, maxptest = 0;
unsigned int *primeset = NULL;

#define is_prime_(p)	((primeset[(p)>>6]>>(((p)>>1)&31))&1)
#define exclude_(p)	( primeset[(p)>>6] &= ((-1) ^ (1<<(((p)>>1)&31))) )

/* can't easily be grown; easiest to re-sieve from 0 every time it has to grow
 * could change inner loop (over n) to start low*low+m (where m makes it a
 * multiple of i)...  too much work.
 */
static void sieve(unsigned int sqrthigh){
	unsigned int n,i;
	unsigned int high = sqrthigh*sqrthigh;

	memset( primeset, -1, sizeof(*primeset) * (1 + (high>>6)) );
//	for( i=0 ; i < (high>>6) ; primeset[i++] = -1 );
	for( i=3 ; i <= sqrthigh ; i += 2){
		if( is_prime_(i) ){
			for( n=i*i; n <= high ; n+=2*i ){
				exclude_(n);
			}
		}
	}
}

static inline int next_prime (int i)
{
	if (i<=2) return 2;
	i |= 1;			// make i odd
	while (i < maxptest && !is_prime_(i)) i+=2;
	assert( i < maxptest ); // FIXME: better error handling

	return i;
}

static void primesetup(int k);

static inline int is_prime(int p){
//#ifndef NDEBUG
	if (p > maxptest){ // TODO: better error handling?
		fprintf(stderr,"spr debug: sieving more primes. old max %d",
			maxptest);
		primesetup( p*2 );
		fprintf( stderr, ", new %d\n", maxptest );
	}
//#endif
	return p>1 && (p==2 || (p & is_prime_(p)));
}

static void primesetup(int k){
	int sqrsize = (int)ceilf(sqrtf(k));
	if( sqrsize & 0xffff0000 ){
		fprintf(stderr, "max prime too big: 0x%x\n", sqrsize );
	}else if (sqrsize < 3) 
		sqrsize = 7;

	// for allspr, just go big the first time
	sqrsize = min( 21, sqrsize );

	if (sqrsize > sieved){
		primeset = xrealloc(primeset, sizeof(*primeset) * 
				    (((sqrsize*sqrsize)>>6) + 1) );
		sieve( sqrsize );
		sieved = sqrsize;
		maxptest = sieved * sieved;
	}
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
static void findlcg( struct lcg *lcg_params, int maxval){
	unsigned int a, b, c, m = maxval;
	int i;

	primesetup (maxval+maxval/2);

//	assert( maxval >= 1 );
	if( m<=6 ){ // will be either 6 or 2.  Just loop in order
		b=0;
		c=1;
	}else{
		int divlimit = m;
		b=1; // b must be a multiple of all of m's prime factors
		if ( !(m%2) ){ 
			b=2;
			while (divlimit%2 == 0) divlimit /= 2;
		}
		for (i=3 ; i <= divlimit ; i+=2 ){
			if( is_prime(i) && m%i == 0 ){
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
		c = next_prime( max(5, (0.5 - sqrtf(3)/6.0)*m - 2));
		while ( m%c == 0 ) c=next_prime(c+1);
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
// 0xcafe
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
	while(42){
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


/************* Library API functions: init and free *****/

/* set state stuff from the tree: taxa, nodes, and an array of pointers to
 * the nodes, so we can map integers to nodes. */
static void initspr( struct spr_tree *state, struct spr_node *tree )
{
	struct spr_node *p, *q;
	int n=0, nsize=15;
	struct spr_node **nodelist = xmalloc( nsize*sizeof(*nodelist) );

	state->nodes = state->taxa = 0;

	p=tree; q=NULL;
	/* non-recursive traversal, using parent pointers
	 * to ascend and detect where we came from */
	while (p){  
		if (p->parent == q){ // first time to the node
			if( n >= nsize ){
				nsize *= 2;
				nodelist = xrealloc( nodelist, nsize*sizeof(*nodelist) );
			}
			nodelist[n] = p;
			n++;

			if (!p->left ^ !p->right){ /* consistency check */
				fprintf( stderr, 
  "libspr: invalid tree detected:\n"
  "internal nodes must have left and right subtrees\n"
  "node \"%s\" has one but not the other.\n", p->data->name );
				state->nodes=-1;
				return;
			}

			if (p->left){ //internal node
				q=p; p=p->left;
			}else{	//leaf
				state->taxa++;
				q=p; p=p->parent;
			}
		}else if (p->left == q){  // returning from subtrees
			q=p; p=p->right;
		}else if (p->right == q){
			q=p; p=p->parent;
		}else{
			assert ( 0 && "shouldn't be here!");
		}
	}

	nodelist = xrealloc( nodelist, n*sizeof(*nodelist) );
	state->nodes = n;	// taxa were counted as we went
	state->nodelist = nodelist; // could have moved on realloc
}


/* return malloc()ed library state, or NULL on error */
struct spr_tree *
spr_init( struct spr_node *root, void (*callback)(struct spr_node *) )
{
	int nnodes;
	struct spr_tree *tree;

	if (!root) return NULL;

	tree = xmalloc( sizeof(*tree) );
	tree->root = root;
	initspr( tree, root );
	// TODO: sort nodelist?

	nnodes = tree->nodes;
	if (nnodes < 4) goto out_err;
// TODO: either be re-entrant or forget about it...
	if (nnodes > (volatile int)sprmapnodes) grow_map( nnodes );

 // set up an LCG that will cover the whole SPR space
	findlcg( &tree->lcg, nnodes*(nnodes-1) );
	tree->lcg.startstate = tree->lcg.state = rand() % tree->lcg.m;

	tree->callback = callback;
	tree->unspr_dest = NULL;

	return tree;
 out_err:
	spr_statefree( tree );
	return NULL;
}


/********** free() functions ***************/

/* depth-first traversal, freeing as we go */
void spr_treefree( struct spr_node *p, int freedata )
{
	if (p->data && freedata) free( p->data );
	if (p->left){
		spr_treefree(p->left, freedata);
		spr_treefree(p->right, freedata);
	}
	free( p );
}	

void spr_statefree( struct spr_tree *p )
{
	if( p->nodelist )
		free( p->nodelist );
	free( p );
}

/* free the private resources allocated by the library */
void spr_staticfree( void )
{
	sprmapnodes = 0;
	free( sprmap );
	sprmap = NULL;

	sieved = 0;
	free( primeset );
	primeset = NULL;
}


void spr_libsprtest( struct spr_tree *st )
{
	printf("nodes = %d\n", st->nodes);
	int i, n;
	n=sprmapnodes;
	puts("sprmap:");
	for( i=0 ; i<n*(n-1) ; i++ ){
		printf("%d %d\n", sprmap[i][0], sprmap[i][1] );
	}
}
