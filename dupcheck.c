/* subtree pruning-regrafting (spr) library
 * Peter Cordes <peter@cordes.ca>, Dalhousie University
 * license: GPLv2 or later
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define SPR_PRIVATE
#include "spr.h"

// return neighbours (not necessarily leaves).
// only nodes near the root have two
static inline struct spr_node *neighbour1( struct spr_node *p ){
	if (isroot(p->parent))
		return (sibling(p))->left;
	else return sibling(p);
}
static inline struct spr_node *neighbour2( struct spr_node *p ){
	if (isroot(p->parent))
		return (sibling(p))->right;
	else if (isroot(p->parent->parent)){
		return sibling(p->parent);
	}else
		return NULL;
}

void checktree( const struct spr_node *p )
{
	if(p){
		assert((p->left && p->right) || (!p->left && !p->right));
		checktree(p->left);
		checktree(p->right);
	}
}


/* compare two trees.  return TRUE if they have the same topology,
 * ignoring the position of the root. i.e. as if they represent unrooted trees.
 * Actually operate on copies of the trees, because we destructively reduce
 * the two trees until we find a difference, or get down to three nodes.
 * (there is only possible unrooted topology for three nodes).
 *
 * A cherry is an internal node with two leaf nodes as children.  If a cherry
 * exists in tree A, it must exist in tree B if they are topologically the same.
 * Our algorithm is: find a cherry in A, find it in B, and then replace it
 * with just one of its constituent leaves.
 */
static int sametopo( struct spr_node *tree1, struct spr_node *tree2, int nnodes, int ntaxa )
{
	struct spr_node *A, *B, *cherryA, *cherryB, *p, *q;
	int retval = FALSE;
//	checktree(tree1);
//	checktree(tree2);
	A = spr_copytree (tree1);
	B = spr_copytree (tree2);
	
	while (ntaxa>3){
//		checktree(A);
//		checktree(B);
		cherryA = A;
		while(42){ // find a cherry in A
			if (isleaf(cherryA->left)){
				if (isleaf(cherryA->right)) break; // found
				else cherryA = cherryA->right;
			}else cherryA = cherryA->left;
		}

		p = spr_searchbypointer(B, cherryA->left->data);
		assert(p /* trees must share a set of data pointers */ );
		cherryB = p->parent; // only a potential cherry so far

		// find the cherry in B and reduce both trees
		if ((q=neighbour1(p))->data == cherryA->right->data ||
		    ((q=neighbour2(p)) && q->data == cherryA->right->data)){
			// if the cherry in B spans the root, always delete the root and
			// the leaf attached to it, with the remaining tree needing no modification
			if (isroot(p->parent) || isroot(q->parent)){
			// point q at the leaf child of the root, which we are going to delete
				if (isroot(p->parent)){
					q = p;
					cherryA->data = cherryA->right->data;
				}else
					cherryA->data = cherryA->left->data;
				B = sibling(q);
				B->parent = NULL;
				free(q->parent);
				free(q);
			}else{
				// the simple case not involving the root.
				cherryB->data = q->data;
				free(cherryB->left);
				free(cherryB->right);
				cherryB->left = cherryB->right = NULL;
				cherryA->data = cherryA->right->data;
			}
			// make the matching change in tree A
			free(cherryA->left);
			free(cherryA->right);
			cherryA->left = cherryA->right = NULL;
		}else
			goto out_FALSE;

		ntaxa--;
		nnodes-=2;
	}

	retval = TRUE;
 out_FALSE:
	spr_treefree(A, FALSE);
	spr_treefree(B, FALSE);
	return retval;
}

struct spr_node *spr_find_dup( struct spr_tree *tree, struct spr_node *root ){
	struct spr_duplist *p = tree->dups;
	for (p=tree->dups ; p ; p=p->next){
		if (sametopo(root, p->tree, tree->nodes, tree->taxa)) break;
	}
	return p ? p->tree : NULL;
}


// same code as procov's pc2spr
struct spr_node *spr_copytree( const struct spr_node *node )
{
	struct spr_node *p = NULL;
	if(node){
		// point bookkeeping pointer at corresponding node
		p = spr_newnode(spr_copytree(node->left),
				spr_copytree(node->right),NULL, node->data);
		if (p->left ) p->left ->parent = p;
		if (p->right) p->right->parent = p;
	}
	return p;
}

int spr_add_dup( struct spr_tree *tree, struct spr_node *root )
{
	if (!spr_find_dup( tree, root )){
		struct spr_duplist *p = xmalloc(sizeof(*p));
		p->tree = spr_copytree(root);
		p->next = tree->dups;
		tree->dups = p;
		return 1;
	}
	return 0;
}
