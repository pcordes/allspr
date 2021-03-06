/* subtree pruning-regrafting (spr) library
 * Peter Cordes <peter@cordes.ca>, Dalhousie University
 * license: GPLv2 or later
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#define SPR_PRIVATE
#include "spr.h"


/************ node search functions ***************/

void inorder(const struct spr_node *p, void (*func)(const struct spr_node *))
{
	if (p->left) inorder(p->left, func);
	func( p );
	if (p->right) inorder(p->right, func);
}

struct spr_node *spr_treesearchbyname( struct spr_tree *t, const char *s )
{
	struct spr_node *p, *array = t->nodelist[0];
	for( p=array ; p < array+t->nodes && 0!=strcmp(s,p->data->name) ; p++ );
	return p;
//	return spr_searchbyname(t->root, s);
}

struct spr_node *spr_treesearch( struct spr_tree *t, const struct spr_node *query )
{
	struct spr_node *p, *array = t->nodelist[0];
	for( p=array ; p < array+t->nodes && p != query ; p++ );
	return p;
//	return spr_search(t->root, query);
}

/* can't count on tree being sorted by name, so search it all */
struct spr_node *spr_searchbyname( struct spr_node *p, const char *s )
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

struct spr_node *spr_search( struct spr_node *tree, const struct spr_node *query )
{
	struct spr_node *q;
	if (tree == query) return tree;

	if (tree->left){
		if ((q = spr_search(tree->left,query))) return q;
		else if ((q = spr_search(tree->right,query))) return q;
	}
	return NULL;
}

struct spr_node *spr_searchbypointer( struct spr_node *tree, const void *query )
{
	struct spr_node *p, *q;

	for( p=tree ; p ; p=p->left ){
		if (p->data == query) return p;
		if ((q = spr_searchbypointer(p->right, query))) return q;
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


/******** dospr: the real SPR function at the heart of the library ********/
/* reattach src (and it's parent node, which would otherwise have to be deleted)
 * to the branch between dest and its parent.  This makes src and dest siblings.
 * return success/fail
 * Only SPRs which would actually break the tree are rejected here.  see spr()
 */
static int dospr( struct spr_node *src, struct spr_node *dest )
{
	// FIXME: use the callback
	struct spr_node *sp = src->parent, *dp = dest->parent;

	if (spr_isancestor(src, dest) || // dest inside the subtree being pruned
	    dest == sp)		// src parent goes with src, so can't be dest
		return FALSE;
	assert( src->parent != NULL /* isancestor should have caught src==root */ );

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


/* A wrapper around dospr():
 * return the tree to its original topology if needed.
 * rejects some useless SPRs (e.g. that don't change the topology)
 * save info so unspr can get back to original topology.
 * returns TRUE if dospr() succeeds and the tree is modified.
 *
 * It's probably easy to make a mess if you call this directly while root
 * moving is going on. e.g. tree->lastspr is checked for some things.
 * Root moving was hacked in as an afterthought.
 */
int spr( struct spr_tree *tree, struct spr_node *src, struct spr_node *dest )
{
	int tmp, unspr_success=FALSE;

	if (tree->unspr_dest){	// back to starting tree
		unspr_success = dospr(tree->unspr_src, tree->unspr_dest);
		tree->root = spr_findroot(tree->unspr_dest);
		if (spr_debug>=2){
			fputs("  unspr back to: ", stderr);
			newickprint(tree->root, stderr);
		}
		assert( unspr_success );
		tree->unspr_dest = NULL;
	}
	if (!src && !dest) return unspr_success;

	// We used to exclude dest==root, but it doesn't break unspr or anything.
	// It always has the same (unrooted) topology as two other trees that spr_next_spr finds.
	// (the root node is the "extra" node, for unrooted vs. rooted tree) */
	if ( !src || !dest ||	// protect against silly callers
	     spr_isancestor(src, dest) || // does this really always catch !(src->parent)?
	     src->parent == dest->parent)  // don't switch siblings
		return FALSE;

	tree->unspr_src = src;
	tree->unspr_dest = sibling(src);

	tmp = dospr(src, dest);
	if (tmp){
		if (!isroot(tree->root)){
			tree->root = spr_findroot(dest);
			if (spr_debug>=2) fputs("allspr: tree has new root!\n", stderr);
		}
		if (spr_debug>=1)
			printf("  did spr %s -> %s\n", src->data->name, dest->data->name);
	}else
		tree->unspr_dest = NULL;

	return tmp;
}


/* sort out a tree where nodes' parent pointers are pointing to nodes that are actually their children.
 * borrowed from procov */
static void spr_organize_tree(struct spr_node *from, struct spr_node *nd)
{
	struct spr_node *p;

	if(isleaf(nd)) return;
	if(from==nd->left){  p=nd->parent; nd->parent=nd->left;  nd->left =p; }
	if(from==nd->right){ p=nd->parent; nd->parent=nd->right; nd->right=p; }
	spr_organize_tree(nd, nd->left);
	spr_organize_tree(nd, nd->right);
}

// insert the root along the branch connecting child to its parent
static void placeroot(struct spr_tree *tree, struct spr_node *child)
{
	if(spr_debug>=5){ spr_treedump(tree, stderr); }
	struct spr_node *r = tree->root;
	if(!tree->rootsave1){ // only update undo info if we were at the original tree
		tree->rootsave1 = r->left;
		tree->rootsave2 = r->right;
	}
	r->left ->parent = r->right;  // extract root from old location
	r->right->parent = r->left;

	r->right = child->parent;
	r->left = child;
	*meinparent(child) = child->parent = r;
	// root is in the tree, but branches are pointing strange directions.
	if(spr_debug>=5){ spr_treedump(tree, stderr);	putc('\n', stderr); }
	spr_organize_tree(NULL, r);
}

// return the tree to it's original state
static void unrootmove(struct spr_tree *tree)
{
	struct spr_node *root = tree->root, *s1 = tree->rootsave1, *s2 = tree->rootsave2;
	if(!s1) return;
	if(s1->parent == root && s2->parent == root){ if(spr_debug>=3)fprintf(stderr, "allspr: something weird probably happened, %s\n", __func__); return; }

	spr_unspr(tree);
	if(s1->parent == s2) placeroot(tree, s1);
	else if(s2->parent == s1) placeroot(tree, s2);
	else assert( FALSE /* two children of old root aren't connected */ );
	tree->rootsave1 = tree->rootsave2 = NULL;
}

// decode an SPR number and do it.
int spr_sprnum(struct spr_tree *tree, int coded_sprnum)
{
	int tmp, sprnum;
	if(!coded_sprnum) return FALSE;
	else if(coded_sprnum>0){ // classic rooted-tree SPRs not spanning the root
		sprnum = coded_sprnum-1;
		if(sprnum > tree->lcg.m) return FALSE;
		if(tree->lastspr < 0) unrootmove(tree);
		tmp = spr(tree, tree->nodelist[ sprmap(sprnum, 0) ],
				tree->nodelist[ sprmap(sprnum, 1) ]);
		return tree->lastspr = tmp ? coded_sprnum : 0;
	}else{ // root moving
		sprnum = (-coded_sprnum)-1;
		int rootpos = sprnum / (tree->nodes*tree->nodes);
		if(rootpos >= tree->nodes) return FALSE;

		spr_unspr(tree);
//		unrootmove(tree);
		struct spr_node *r = tree->root, *c = tree->nodelist[rootpos];
//		if(isleaf(c) || r==c) return FALSE;
		if(r==c) return FALSE; //ROOTMOVE ONLY
		if(c->parent != r) placeroot(tree, c); // don't repeat ourselves

		sprnum = sprnum % (tree->nodes*tree->nodes);
		tmp = sprnum % tree->nodes; sprnum /= tree->nodes;
		tmp = spr(tree, tree->nodelist[tmp], tree->nodelist[sprnum]);
		tree->lastspr = coded_sprnum; // root is out of place whether we succeed or not
		return tmp ? coded_sprnum : 0;
	}
}


/****************** SPR iteration ******************/

/* return 0 for all done, else 1+SPR number.  Zero makes a nicer sentinel than
 * UINT_MAX for users of the library, but beware of the offset when debugging.
 *
 * root moving: loop through normal SPRs, then return negative sprnums.
 */
int spr_next_spr( struct spr_tree *tree )
{
	int tmp = FALSE;
	unsigned sprnum;

//	tree->rootmove=1; //ROOTMOVE ONLY
	if(tree->rootmove == 0){
		do{  // try SPRs until we find a legal one, or loop to beginning
			sprnum = lcg( &tree->lcg );
			if(UINT_MAX == sprnum) break;
			tmp = spr_sprnum(tree, sprnum+1);
			if (tmp && tree->dups)
				tmp = spr_add_dup(tree, tree->root);
		}while(!tmp);
		if(UINT_MAX != sprnum) return 1 + sprnum;
	}
#ifdef NO_ROOT_MOVING
	return FALSE;
#else
	do{
		if(tree->rootmove++ > (tree->nodes*tree->nodes*tree->nodes)) return FALSE;
		tmp = spr_sprnum(tree, -(int)(tree->rootmove+1));
		if (tmp && tree->dups)
			tmp = spr_add_dup(tree, tree->root);
	}while(!tmp);
	return tree->lastspr = -(int)(tree->rootmove+1);
#endif
}

/* move to a new tree.  also called from spr_init() */
void spr_apply(struct spr_tree *tree)
{
	tree->unspr_dest = tree->unspr_src = NULL;
	tree->lcg.startstate = UINT_MAX;
	tree->rootmove = tree->lastspr = 0;
	tree->rootsave1 = tree->rootsave2 = NULL;
}

int spr_apply_sprnum(struct spr_tree *tree, int sprnum)
{
	int tmp;
	if((tmp = spr_sprnum(tree, sprnum)))
		spr_apply(tree);
	return tmp;
}
