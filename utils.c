/* standard error-exit wrappers for memory allocation
 * allspr library copyright Peter Cordes <peter@cordes.ca>
 * license: GPLv2 or later
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#define SPR_PRIVATE // spr.h warns without this or SPR_NODE_DATAPTR_TYPE defined
#include "spr.h"

void *xcalloc (size_t n, size_t s)
{
	void *p = calloc (n, s);
	if (!p){
		perror("allocating memory");
		exit (2);
	}
	return p;
}

void *xmalloc (size_t s)
{
	void *p = malloc (s);
	if (!p){
		perror("allocating memory");
		exit (2);
	}
	return p;
}

void *xrealloc (void *p, size_t n)
{
	p = realloc(p, n);
	if (!p){
		perror("allocating memory");
		exit (2);
	}
	return p;
}
