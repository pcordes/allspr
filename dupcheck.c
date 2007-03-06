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
 *
 * As an optimization, we operate on trees where all the node structs are known to be
 * in two contiguous arrays.  Instead of a spr_searchybypointer, we linearly search the
 * array.  Instead of free() we set ->data = NULL (so searches don't find deleted nodes).
 */
static int sametopo( struct spr_node *array1, struct spr_node *array2, size_t asize, int ntaxa )
{
	struct spr_node *A, *B, *cherryA, *cherryB, *p, *q;
	int retval = FALSE;

	A = spr_findroot(array1);
	B = spr_findroot(array2);

	while (ntaxa>3){
		cherryA = A;
		while(42){ // find a cherry in A
			if (isleaf(cherryA->left)){
				if (isleaf(cherryA->right)) break; // found
				else cherryA = cherryA->right;
			}else cherryA = cherryA->left;
		}

		//p = spr_searchbypointer(B, cherryA->left->data);
		for( p=array2 ; p<=array2+asize && p->data!=cherryA->left->data ; p++ );
		assert( p && p <= array2+asize /* trees must share a set of data pointers */ );
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
				q->parent->data = q->data = NULL;
				//free(q->parent); free(q);
			}else{
				// the simple case not involving the root.
				cherryB->data = q->data;
				//free(cherryB->left); free(cherryB->right);
				cherryB->left->data = cherryB->right->data = NULL;
				cherryB->left = cherryB->right = NULL;
				cherryA->data = cherryA->right->data;
			}
			// make the matching change in tree A
			//free(cherryA->left); free(cherryA->right);
			// cherryA->left->data = cherryA->right->data = NULL; // not needed
			cherryA->left = cherryA->right = NULL;
		}else
			goto out_FALSE;

		ntaxa--;
	}

	retval = TRUE;
 out_FALSE:
	//spr_treefree(A, FALSE);
	//spr_treefree(B, FALSE);
	return retval;
}

struct spr_node *spr_find_dup( struct spr_tree *tree, struct spr_node *root ){
	struct spr_duplist *p = tree->dups;
	struct spr_node *A, *B, *saveA;
	int n = tree->nodes, tmp;

	assert( tree->nodes == spr_countnodes(root) );

	A = xmalloc(n*sizeof(*A));
	B = xmalloc(n*sizeof(*B));

	tmp = spr_copytoarray(A, root);
	assert( n == tmp /* copytoarray had better copy the right number of nodes */ );

	saveA = xmalloc(n*sizeof(*A));
	memcpy(saveA, A, n*sizeof(*A));
	for (p=tree->dups ; p ; p=p->next){
		memcpy(A, saveA, n*sizeof(*A));
		spr_copytoarray(B, p->tree);
		if (sametopo(A, B, n, tree->taxa)) break;
	}
	free(saveA);
	free(B);
	free(A);
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

size_t spr_copytoarray( struct spr_node *A, const struct spr_node *node )
{
	struct spr_node *p = A, *left, *right;
	if(isleaf(node)){
		p->left = p->right = NULL;
		p->data = node->data;
		return 1;
	}else{
		// depth-first traversal.  leafs are close to the beginning (for search).
		// do left then right subtree, then set up their ancestor.
		p += spr_copytoarray(p, node->left); left = p-1;
		p += spr_copytoarray(p, node->right); right = p-1;
		p->left = left;	p->right = right; p->parent = NULL;
		p->data = node->data;

		p->left->parent = p->right->parent = p;
		return p-A + 1;
	}
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
