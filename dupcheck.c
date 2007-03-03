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

static int sametopo( struct spr_node *p1, struct spr_node *p2 )
{
//	return TRUE;
	return FALSE;
}

struct spr_node *spr_find_dup( struct spr_tree *tree, struct spr_node *root ){
	struct spr_duplist *p = tree->dups;
	for (p=tree->dups ; p ; p=p->next){
		if (sametopo(root, p->tree)) break;
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
