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


struct spr_nodename{
	char *name;
};

/* libspr requires that there be a char * stored at the address pointed to by 
 * data.  This can be the first member of a custom-defined structure,
 * of course.  The char * is used as the name of the node. 

 * internal nodes in phylogenetic trees only exist with two children,
 * so left!=NULL implies right!=NULL.  The root node has parent == NULL
 */

struct spr_node{
	struct spr_node *parent, *left, *right;
#ifdef SPR_NODE_DATAPTR_TYPE
	SPR_NODE_DATAPTR_TYPE *data;
#else
	void *data;
#endif

/* Union declaration makes things too clumsy :(
	union {
		struct spr_nodename *def;
		void *custom;
	} data;
*/
};

/* a linear congruential generator is used to generate all integers 
 * between 0 and N without repetition, in a pseudo-random order,
 * to determine the order to try SPRs in */
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

	struct lcg lcg;

	void (*callback)(struct spr_node *);
	int nodes;
	int taxa;
};






/********** Public Functions ***************/


// hopefully this doesn't cause too much of a problem,
// since progs that link libspr can just use these
void *xmalloc (size_t s);  // perror and exit on error
void *xcalloc (size_t n, size_t s);
void *xrealloc (void *p, size_t n);

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


// Library API stuff
/* call this on the root of the tree before using the other functions.
 * callback is a function to be called on nodes that undergo an SPR.
 * the return value is a pointer to a malloc()ed struct holding info about the
 * tree, such as number of nodes, a pointer to the root, etc. */
struct spr_tree *spr_init( struct spr_node *tree, void (*callback)(struct spr_node *) );

void spr_statefree( struct spr_tree *p ); /* use _instead_ of free( p ). 
   * frees just the struct spr_tree and related stuff, not the tree itself */
void spr_staticfree( void ); // free memory internally allocated by the lib
void spr_treefree( struct spr_node *tree, int freenodedata );
/* traverse the tree, calling free() on all the nodes, and optionally on
 * all the .data pointers, too. */


// Tree topology
/* do an spr, attaching the subtree rooted at src 
 * to the branch between dest and its parent. return success(TRUE)/fail(FALSE)
 */
int spr( struct spr_tree *tree, struct spr_node *src, struct spr_node *dest );

/* return 0 for all done, else a positive SPR number */
int spr_next_spr( struct spr_tree *tree );

// IO
char *newick( const struct spr_node *subtree ); // return a malloc()ed string. no bl
void newickprint( const struct spr_node *subtree ); // print to stdout
void treeprint(const struct spr_node *p); // in-order traversal printing to stdout


// debugging
void spr_libsprtest( struct spr_tree *state );


/***** internal functions that might be useful  ****/
struct spr_node *spr_findroot( const struct spr_node *p ); // walk parent ptr to root
int spr_countnodes( const struct spr_node *p );
int spr_isancestor( const struct spr_node *ancestor, const struct spr_node *child );

/* Return a pointer to the found node, or NULL.  inefficient because the tree
 * isn't sorted on the name or the address of the nodes. */
struct spr_node *spr_search( const struct spr_node *HAYSTACK, const struct spr_node *NEEDLE);
struct spr_node *spr_searchbyname( const struct spr_node *HAYSTACK, const char *NEEDLE );

// These will be faster (once they're written)
struct spr_node *spr_treesearchbyname(const struct spr_tree *t, const char *s );
struct spr_node *spr_treesearch(const struct spr_tree *t, const struct spr_node *query );


#ifdef __GNUC__
static inline struct spr_node *spr_newnode( struct spr_node *parent,struct spr_node *left, struct spr_node *right, typeof(parent->data) data)
{
  struct spr_node *p;
  p = xmalloc(sizeof(*p));
  p->parent=parent; p->left=left; p->right=right;
  p->data=data;
  return p;
}
#endif
