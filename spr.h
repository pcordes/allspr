/* subtree pruning-regrafting (spr) library
 * Peter Cordes <peter@cordes.ca>, Dalhousie University
 * license: GPLv2 or later
 */

#define ALLSPR_VERSION "1.3"

#ifdef SPR_PRIVATE // intended for internal library use.
// this has to be up here near the beginning of spr.h

// maybe TODO: have procov compile us with -D..., so the allspr code
// builds stand-alone by default.  Or just take out the printing/debugging
// stuff so the ->data pointer can be opaque to allspr.
#  define SPR_PROCOV_DATA
#  ifdef SPR_PROCOV_DATA
#    define nom name
#    include "../procov/mynhmlg.h"
#    define SPR_NODE_DATAPTR_TYPE struct noeud
#  else
#    define SPR_NODE_DATAPTR_TYPE struct spr_nodename
#  endif
#endif // SPR_PRIVATE

// make sure we always have this defined even without SPR_PRIVATE.
#ifndef SPR_NODE_DATAPTR_TYPE
#  define SPR_NODE_DATAPTR_TYPE void
#endif

struct spr_nodename{ char *name; };

#ifndef TRUE
#define TRUE (1)
#define FALSE (0)
#endif

// backward compat for old gcc, from gcc info page.
#if __STDC_VERSION__ < 199901L
# if __GNUC__ >= 2
#  define __func__ __FUNCTION__
# else
#  define __func__ "<unknown>"
# endif
#endif

#ifndef max
#ifdef __GNUC__
#define max(a,b) ({ typeof (a) _a = (a); typeof (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ typeof (a) _a = (a); typeof (b) _b = (b); _a < _b ? _a : _b; })
#else
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif
#endif


/* libspr requires that there be a char * stored at the address pointed to by 
 * data.  This can be the first member of a custom-defined structure,
 * of course.  The char * is used as the name of the node. 
 * If you compile a version of the library with spr_private.h including
 * your own headers, you can have the ->name structure member anywhere.

 * internal nodes in phylogenetic trees only exist with two children,
 * so left!=NULL implies right!=NULL.  The root node has parent == NULL
 */

struct spr_node{
	struct spr_node *left, *right, *parent;
	SPR_NODE_DATAPTR_TYPE *data;
};

struct spr_duplist{
	struct spr_node *tree;
	struct spr_duplist *next;
};

/* a linear congruential generator is used to generate all integers 
 * between 0 and N without repetition, in a pseudo-random order,
 * to determine the order to try SPRs in.
 * startstate = UINT_MAX means we haven't returned any numbers yet.
 */
struct lcg {
	unsigned int state;
	unsigned int a, c, m;
	unsigned int startstate;
};

struct spr_tree{
	struct spr_node *root;
	struct spr_node **nodelist; // not sorted
//	struct spr_node **nodesbyname;
	struct spr_node *unspr_dest;
	struct spr_node *unspr_src;  // could be an index into nodelist
	struct spr_duplist *dups;
	void (*callback)(struct spr_node *);  // not implemented
	struct spr_node *rootsave1, *rootsave2;
	unsigned int rootmove;
	struct lcg lcg;

	int lastspr;
	int nodes;
	int taxa;
};






/********** Public Functions ***************/

// hopefully this doesn't cause too much of a problem,
// since progs that link libspr can just use these
void *xmalloc (size_t s);  // perror and exit on error
void *xcalloc (size_t n, size_t s);
void *xrealloc (void *p, size_t n);


// Library API stuff
/* call this on the root of the tree before using the other functions.
 * callback is a function to be called on nodes that undergo an SPR.
 * the return value is a pointer to a malloc()ed struct holding info about the
 * tree, such as number of nodes, a pointer to the root, etc. */
struct spr_tree *spr_init( struct spr_node *tree, void (*callback)(struct spr_node *), int allow_dups );

void spr_statefree( struct spr_tree *p ); /* use _instead_ of free( p ). 
   * frees just the struct spr_tree and related stuff, not the tree itself */
void spr_staticfree( void ); // free memory internally allocated by the lib
void spr_treefree( struct spr_node *tree, int freenodedata );
/* traverse the tree, calling free() on all the nodes, and optionally on
 * all the .data pointers, too. */


/******** Tree topology ********/
/* do an spr, attaching the subtree rooted at src to the branch between
 * dest and its parent.  src and dest become siblings.
 * return TRUE(success)/FALSE(invalid or useless spr)
 * Will be undone by the next spr call, because unspr info is saved.
 */
int spr( struct spr_tree *tree, struct spr_node *src, struct spr_node *dest );
int spr_sprnum(struct spr_tree *tree, int sprnum);
/* return 0 for all done, else a positive or negative SPR number */
int spr_next_spr( struct spr_tree *tree );
/* return the tree to its original topology */
static inline int spr_unspr(struct spr_tree *tree){ return spr(tree, NULL, NULL); }
// make last SPR permanent spr: don't save unspr info.  preserves duplicate checking list.
// resets the spr_next_spr() iterator.
void spr_apply(struct spr_tree *tree);
int spr_apply_sprnum(struct spr_tree *tree, int sprnum);

/******** Duplicate checking ********/
/* add a tree topology to the dup list (copies the tree).
 * ->data pointers in nodes must be unique
 * return: TRUE if added ok (implies not already present)
 * FALSE if a dup of a tree already there, so not added. */
int spr_add_dup( struct spr_tree *tree, struct spr_node *root );
/* pointer to root of dup tree, or NULL if not a dup. */
struct spr_node *spr_find_dup( struct spr_tree *tree, struct spr_node *root );

/******** IO ********/
char *newick( const struct spr_node *subtree ); // return a malloc()ed string. no bl
#ifdef BUFSIZ // proxy for stdio.h.  skip these if we don't have FILE.
void newickprint(const struct spr_node *subtree, FILE *stream);
void treeprint(const struct spr_node *p, FILE *stream); // in-order traversal
void spr_treedump(const struct spr_tree *t, FILE *stream); // dump t->nodelist with names for all pointers
#endif // stdio


/******** Debugging ********/
int spr_debug;
static inline void spr_setdebug(int level){ spr_debug=level; }
void spr_libsprtest(struct spr_tree *state);

/***** internal functions that might be useful  ****/
static inline struct spr_node *spr_findroot( struct spr_node *p )
{ /* follow the linked list all the way up */
	while( p->parent != NULL ) p = p->parent;
	return p;
}
int spr_countnodes( const struct spr_node *p );
int spr_isancestor( const struct spr_node *ancestor, const struct spr_node *child );

// xmalloc()ed copy of each node, with ->data pointers the same.
struct spr_node *spr_copytree(const struct spr_node *node);
/* Copy a tree to an array, which must be of size >= spr_countnodes(node).
 * Avoids malloc overhead for each node.  Returns # of nodes copied */
size_t spr_copytoarray(struct spr_node *array, const struct spr_node *root);

/* Return a pointer to the found node, or NULL.  Not very fast because it has
 * to traverse the whole tree with a recursive function */
struct spr_node *spr_search( struct spr_node *HAYSTACK, const struct spr_node *NEEDLE);
struct spr_node *spr_searchbyname( struct spr_node *HAYSTACK, const char *NEEDLE );
// find the node that has the same ->data pointer.
struct spr_node *spr_searchbypointer( struct spr_node *HAYSTACK, const void *NEEDLE );

// These search the node list, which is faster than traversing the tree
struct spr_node *spr_treesearchbyname( struct spr_tree *t, const char *s );
struct spr_node *spr_treesearch( struct spr_tree *t, const struct spr_node *query );

static inline struct spr_node *spr_newnode(struct spr_node *left, struct spr_node *right, struct spr_node *parent, void *data)
{
	struct spr_node *p = xmalloc(sizeof(*p));
	p->parent=parent; p->left=left; p->right=right;
	p->data=data;
	return p;
}


#ifdef SPR_PRIVATE // intended for internal library use.  might be useful generally
unsigned int lcg(struct lcg *lcgp);
void findlcg(struct lcg *lcg_params, int maxval);
void spr_lcg_staticfree(void);

int (*sprmap_table)[2];
int sprmapnodes; // number of nodes the map is good for
static inline int sprmap(int sprnum, int pos){ return sprmap_table[sprnum][pos]; }

// node relationship helpers
#define isleaf(p) (!(p)->left)
#define isroot(p) (!(p)->parent)

/* find whether a node is pointed to by the left (0) or right (1)
 * pointer in its parent */
static inline int isrightchild(const struct spr_node *p){
	return p == p->parent->right; }
static inline struct spr_node **meinparent(const struct spr_node *p){
	return (isrightchild(p)? &p->parent->right : &p->parent->left); }
static inline struct spr_node **siblinginparent(const struct spr_node *p){
	return (isrightchild(p)? &p->parent->left  : &p->parent->right);}
#define sibling(p) (*siblinginparent(p))

#endif // SPR_PRIVATE
