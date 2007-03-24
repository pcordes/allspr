/* subtree pruning-regrafting (spr) library
 * Peter Cordes <peter@cordes.ca>, Dalhousie University
 * license: GPLv2 or later
 */

/* setup routines that are called when giving the library a new tree */

#define _GNU_SOURCE
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <assert.h>

#define SPR_PRIVATE
#include "spr.h"

/* format of sprmap (transposed):
01 0122 012333  src
10 2201 333012  dest
 size 2: 2 2
 size 3: 4 6
 size 4: 6 12
space required: sum( i=2..n, 2*(i-1) ) = n*(n-1).
* src = root entries are useless, but the root can move from 0.
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

// only called from spr_init.
static void grow_map( int nnodes )
{
	int i, n;	      // the number we're currently filling in
	int newsize = nnodes*(nnodes-1);

	if ((n=sprmapnodes) < 2) n=2; // start where we left off last time

	if (nnodes < 2 || nnodes <= sprmapnodes) return; 
	sprmap = xrealloc( sprmap, newsize * sizeof(*sprmap) ); // realloc(NULL,...) is malloc
	sprmapnodes = nnodes;
 
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

	state->nodelist = xrealloc( nodelist, n*sizeof(*nodelist) );
	state->nodes = n;	// taxa were counted as we went
}


/* return malloc()ed library state, or NULL on error */
struct spr_tree *
spr_init( struct spr_node *root, void (*callback)(struct spr_node *), int dup )
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

 // set up and initialize an LCG that will cover all source/dest pairs
	findlcg( &tree->lcg, nnodes*(nnodes-1) );

	tree->callback = callback;
	tree->unspr_dest = NULL;

	if(dup) tree->dups = NULL;
	else{
		tree->dups = xmalloc(sizeof(*tree->dups));
		tree->dups->next = NULL;
		tree->dups->tree = xmalloc(nnodes * sizeof(*tree->dups->tree));
		spr_copytoarray(tree->dups->tree, tree->root);
	}

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

void spr_statefree( struct spr_tree *tree )
{
	for( struct spr_duplist *d = tree->dups ; d ; d=d->next )
		free(d->tree); // the dup list is block-allocated
	if (tree->nodelist)
		free(tree->nodelist);
	free(tree);
}

/* free the private resources allocated by the library */
void spr_staticfree( void )
{
	sprmapnodes = 0;
	free( sprmap );
	sprmap = NULL;
	spr_lcg_staticfree();
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
