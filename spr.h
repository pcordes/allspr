/* subtree pruning-regrafting (spr) library
 * Peter Cordes <peter@cordes.ca>, Dalhousie University
 * license: GPLv2 or later
 */

#define ALLSPR_VERSION "1.1"

#ifdef SPR_PRIVATE // intended for internal library use.
int (*sprmap)[2];
int sprmapnodes;

// this has to be up here near the beginning of spr.h
#define nom name
#include "../procov/mynhmlg.h"

//#define SPR_NODE_DATAPTR_TYPE struct spr_nodename
#define SPR_NODE_DATAPTR_TYPE struct noeud
//#define name nom
#endif // SPR_PRIVATE


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
#define max(a,b) \
       ({ typeof (a) _a = (a); \
           typeof (b) _b = (b); \
         _a > _b ? _a : _b; })

#define min(a,b) \
       ({ typeof (a) _a = (a); \
           typeof (b) _b = (b); \
         _a < _b ? _a : _b; })
#else
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif
#endif


struct spr_nodename{
	char *name;
};

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
#ifdef SPR_NODE_DATAPTR_TYPE
	SPR_NODE_DATAPTR_TYPE *data;
#else
	void *data;
#endif
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

	struct lcg lcg;

	void (*callback)(struct spr_node *);  // not implemented
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


// Tree topology
/* do an spr, attaching the subtree rooted at src 
 * to the branch between dest and its parent. return success(TRUE)/fail(FALSE)
 * Will be undone by the next spr call, because unspr info is saved.
 */
int spr( struct spr_tree *tree, struct spr_node *src, struct spr_node *dest );

/* return 0 for all done, else a positive SPR number */
int spr_next_spr( struct spr_tree *tree );
/* return the tree to its original topology, without doing a new SPR */
int spr_unspr( struct spr_tree *tree );
// do a permanent spr: don't save unspr info.  preserves duplicate checking list.
// resets the spr_next_spr() iterator.
int spr_apply_sprnum(struct spr_tree *tree, int sprnum);
// keep current topology, whatever it is.
void spr_apply(struct spr_tree *tree);

// duplicate checking
/* add a tree topology to the dup list (copies the tree).
 * ->data pointers in nodes must be unique
 * return: TRUE if added ok (implies not already present)
 * FALSE if a dup of a tree already there, so not added. */
int spr_add_dup( struct spr_tree *tree, struct spr_node *root );
/* pointer to root of dup tree, or NULL if not a dup. */
struct spr_node *spr_find_dup( struct spr_tree *tree, struct spr_node *root );

// IO
char *newick( const struct spr_node *subtree ); // return a malloc()ed string. no bl
void newickprint( const struct spr_node *subtree ); // print to stdout
void treeprint(const struct spr_node *p); // in-order traversal printing to stdout


// debugging
void spr_libsprtest( struct spr_tree *state );


/***** internal functions that might be useful  ****/
static inline struct spr_node *spr_findroot( struct spr_node *p )
{ /* follow the linked list all the way up */
	while( p->parent != NULL ) p = p->parent;
	return p;
}
int spr_countnodes( const struct spr_node *p );
int spr_isancestor( const struct spr_node *ancestor, const struct spr_node *child );

// xmalloc()ed copy of each node, with ->data pointers the same.
struct spr_node *spr_copytree( const struct spr_node *node );
// copy a tree to an array of spr_countnodes(node) struct spr_nodes.  avoids malloc overhead for each node
size_t spr_copytoarray( struct spr_node *A, const struct spr_node *node );

/* Return a pointer to the found node, or NULL.  inefficient because the tree
 * isn't sorted on the name or the address of the nodes. */
struct spr_node *spr_search( struct spr_node *HAYSTACK, const struct spr_node *NEEDLE);
struct spr_node *spr_searchbyname( struct spr_node *HAYSTACK, const char *NEEDLE );
// find the node that has the same ->data pointer.
struct spr_node *spr_searchbypointer( struct spr_node *HAYSTACK, const void *NEEDLE );

// These will be faster (once they're written)
struct spr_node *spr_treesearchbyname( struct spr_tree *t, const char *s );
struct spr_node *spr_treesearch( struct spr_tree *t, const struct spr_node *query );


#ifdef __GNUC__
static inline struct spr_node *spr_newnode( struct spr_node *left, struct spr_node *right, struct spr_node *parent, typeof(parent->data) data)
{
  struct spr_node *p;
  p = xmalloc(sizeof(*p));
  p->parent=parent; p->left=left; p->right=right;
  p->data=data;
  return p;
}
#endif



#ifdef SPR_PRIVATE // intended for internal library use.  might be useful generally
unsigned int lcg(struct lcg *lcgp);
void findlcg(struct lcg *lcg_params, int maxval);
void spr_lcg_staticfree(void);


// node relationship helpers
#define isleaf(p) (!(p)->left)
#define isroot(p) (!(p)->parent)

/* find whether a node is pointed to by the left (0) or right (1)
 * pointer in its parent */
static inline int isrightchild( const struct spr_node *p )
{
	return p == p->parent->right;
}

static inline struct spr_node **meinparent( const struct spr_node *p )
{
	return (isrightchild(p)? &p->parent->right : &p->parent->left);
}

static inline struct spr_node **siblinginparent( const struct spr_node *p )
{
	return (isrightchild(p)? &p->parent->left  : &p->parent->right);
}
#define sibling(p) (*siblinginparent(p))


#endif // SPR_PRIVATE
