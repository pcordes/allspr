/* subtree pruning-regrafting (spr) library
 * Peter Cordes <peter@cordes.ca>, Dalhousie University
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "spr_private.h"
#include "spr.h"


/***** Randomness functions *****/
// all the setup for this is in init.c
unsigned int lcg(struct lcg *lcgp) {
	unsigned long a = lcgp->a, c = lcgp->c, m = lcgp->m;

	lcgp->state = (a*lcgp->state + c) % m;
	return lcgp->state;
}


/************ Tree functions ***************/


void inorder(const struct spr_node *p, const void (*func)(const struct spr_node *))
{
	if (p->left) inorder(p->left, func);
	func( p );
	if (p->right) inorder(p->right, func);
}

struct spr_node *spr_treesearchbyname( const struct spr_tree *t, const char *s )
{
// FIXME: binary search the nodelist...
	return NULL;
}

struct spr_node *spr_treesearch( const struct spr_tree *t, const struct spr_node *query )
{
// FIXME: search the node list
	return NULL;
}

/* can't count on tree being sorted by name, so search it all */
struct spr_node *spr_searchbyname( const struct spr_node *p, const char *s )
{
	struct spr_node *q;
	if (!p->data)
		puts("node with NULL data");
	else if (0 == strcmp(s, p->data->name))
		return p;

	/* nodes have either two or zero children */
	if (p->left){
		if ((q = spr_searchbyname(p->left,s))) return q;
		else if ((q = spr_searchbyname(p->right,s))) return q;
	}
	return NULL;
}

struct spr_node *spr_search( const struct spr_node *tree, const struct spr_node *query )
{
	struct spr_node *q;
	if (tree == query) return tree;

	if (tree->left){
		if ((q = spr_search(tree->left,query))) return q;
		else if ((q = spr_search(tree->right,query))) return q;
	}
	return NULL;
}

int spr_countnodes( const struct spr_node *p )
{
	if (p->left)
		return 1 + spr_countnodes(p->left)
			+ spr_countnodes(p->right);
	else
		return 1;
}

/* check if ancestor is an ancestor of p (but not vice versa) */
int spr_isancestor( const struct spr_node *ancestor, const struct spr_node *p )
{
	while( p!=NULL ){
		if (p == ancestor) return TRUE;
		p=p->parent;
	}
	return FALSE;
}

/* follow the linked list all the way up */
struct spr_node *spr_findroot( const struct spr_node *p )
{
	while( p->parent != NULL ) p = p->parent;
	return p;
}

/* find whether a node is pointed to by the left (0) or right (1)
 * pointer in its parent */
static int isrightchild( const struct spr_node *p )
{
	return p == p->parent->right;
}

static struct spr_node **meinparent( const struct spr_node *p )
{
	return isrightchild(p)? &(p->parent->right) : &(p->parent->left);
}

/* reattach src (and it's parent node,
 * which would otherwise have to be deleted)
 * to the branch between dest and it's parent.
 * return success/fail
 */
static int dospr( struct spr_node *src, struct spr_node *dest )
{
	// FIXME: use the callback
	struct spr_node *sp = src->parent, *dp = dest->parent;

	// check for impossible SPRs; the wrapper catches non-useful ones
	// (which are useful when undoing a previous SPR).
	
//	if (!sp || spr_isancestor(sp, dest)) return FALSE;
	if (!sp) return FALSE;  // can't SPR the root
	if (spr_isancestor(src, dest)) return FALSE;

	// This can result in dest->parent having two pointers to sp,
	// resulting in getting the mirror image unspr, for example with cox2
	// spr number 68 (int2->cox2_trybb), because isrightchild will be true!
	if (dp) *meinparent(dest) = sp;
	dest->parent = sp;

	if( isrightchild(src) ){  // TODO: sort?
		sp->left->parent = sp->parent;
		if (sp->parent) *meinparent(sp) = sp->left;
		sp->left = dest;
	}else{ // src is a left child
		sp->right->parent = sp->parent;
		if (sp->parent) *meinparent(sp) = sp->right;
		sp->right = dest;
	}

	sp->parent = dp;
	return TRUE;
}


int spr_unspr( struct spr_tree *tree )
{
	return spr( tree, NULL, NULL );
}


/* A wrapper returns the tree to its original topology if needed, then 
 * returns TRUE if the requested SPR was done, else FALSE.  (except when
 * src=dest=NULL, return TRUE if the tree was modified back to the original.)
 * some of the allowed-spr checks are duplicated in dospr().
 */
int spr( struct spr_tree *tree, struct spr_node *src, struct spr_node *dest )
{
	int tmp, unspr_success=FALSE;

	if (tree->unspr_dest){	// back to starting tree
		unspr_success = dospr( tree->unspr_src, tree->unspr_dest );
		printf("  unspr back to: "); //DEBUG
		newickprint( spr_findroot(tree->root) );
		assert( unspr_success );
		tree->unspr_dest = NULL;
	}
	if (!src && !dest && unspr_success) return TRUE;

	if ( !(src->parent) || !(dest->parent) ||
	     src->parent == dest->parent ||  // don't switch siblings
	     src->parent == dest
	     // ancestor(dest,src) is ok, though
		){
		return FALSE;
	}
	/* src or dest == root is not useful because we actually deal
	 * with internal parent nodes, which the root doesn't have. */
	/* FIXME: if( src is too close to dest ) return; more checking */

	tree->unspr_src = src;
	tree->unspr_dest = isrightchild(src) ? 
		src->parent->left : src->parent->right; // sibling of src.

	tmp = dospr( src, dest );
	if (tmp){
		if (tree->root->parent){
			tree->root = spr_findroot( dest );
			fputs("allspr: tree has new root!\n", stderr);
		}
		// DEBUG
		printf( "  did spr %s -> %s\n", src->data->name, dest->data->name);
	}else
		tree->unspr_dest = NULL;

	return tmp;
}


/****************** SPR iteration ******************/

/* return 0 for all done, else a positive SPR number */
int spr_next_spr( struct spr_tree *tree )
{
	int tmp;
	unsigned sprnum;

	do{  // try SPRs until we find a legal one, or loop to beginning
		sprnum = lcg( &tree->lcg );
		tmp = spr(tree, tree->nodelist[ sprmap[sprnum][0] ],
				tree->nodelist[ sprmap[sprnum][1] ]);
	}while( !tmp && sprnum != tree->lcg.startstate );

	if (tree->lcg.startstate == sprnum){
		return 0;
	}else
		return 1 + sprnum;
}
