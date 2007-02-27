/* generate all SPRs of a given tree.
 * 2004/11/23.  Peter Cordes <peter@cordes.ca>
 * license: GPLv2 or later
 *
 * about the name: brontler appeared during a late-night game of Scrabble :)
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>


struct nodedata{
	char *name;
	int dirty;		/* needs likelihood update after SPR? */
	float bl;		/* branch length to parent */
	float (*dna)[4];	/* 4xn matrix */
};

/* request that spr_node.data be declared as a pointer to the struct
 * that we will actually store in it, instead of a void *
 * so we don't have to cast all the time */
//#define SPR_NODE_DATAPTR_TYPE struct nodedata
#include <spr_private.h>
#include <spr.h>


/* a "constructor".  returns a malloc()ed spr_node */
struct spr_node *newnode(char *name)
{
	struct spr_node *p = xcalloc(1, sizeof(*p));
//	struct nodedata *d = xcalloc(1, sizeof(*d));
	typeof(p->data) d = xcalloc(1, sizeof(*d));

	p->left=p->right=p->parent=NULL;
	p->data = d;

	if (name) strcpy (d->name, name);
//	d->name=name;
//	d->dirty=0;
//	d->bl=0;
//	d->dna=NULL;

	return p;
}


/* grammar:
 *   subtree: (subtree,subtree) | taxon
 *   taxon: name | name:bl
 *   name: string not including (, ) or :.
 *   bl: floating point number
 *
 * half-assed simple parser:
 * doesn't quite handle white space everywhere it should.
 * doesn't handle multifurcating trees.
 * doesn't handle names on internal nodes.
 * should be re-written to use yacc/bison and lex/flex.
 *
 * returns root node of a tree 
 * len is the number of characters this subtree was */
struct spr_node *parsenewick( char *str, int *len )
{
	static char nextname[] = { 'A', '\0' };
	int tmp;
//	char *name;
	struct spr_node *node = newnode( NULL );
	typeof (node->data) data = node->data;

	if (*str == '('){  // internal node
//		data->name = strdup(nextname);
		strcpy (data->name, nextname);
		++*nextname;

		node->left = parsenewick(str+1, &tmp);
		node->left->parent = node;
		*len = tmp+1;

		tmp = -1;  // skip white space around the comma
		sscanf(str+*len, " , %n", &tmp);
		if (-1 == tmp)
			fprintf(stderr,"no comma following left subtree: \"%s\"\n", str+*len);
		else
			*len += tmp;

		node->right = parsenewick(str+*len, &tmp);
		node->right->parent = node;
		*len += tmp;

		if (')' != str[*len])
			fprintf(stderr,"no close paren following right subtree: \"%s\"\n", str+*len);
		else
			++*len;

	}else if (strspn(str, " -_.abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")){
		// leaf node
		tmp = 0;  // let scanf handle taxon names
//		if ( 1 > sscanf( str, " %a[^][)(:;, \f\n\r\t\v] %n", &data->name, len))
		if ( 1 > sscanf( str, " %[^][)(:;, \f\n\r\t\v] %n", data->name, len))
			fprintf(stderr, "bad newick element: \"%s\"\n",str);
	}else{ // not a node!?
		fprintf(stderr, "bad newick element, not an open paren or a taxon name: \"%s\"\n",str);
		free( node );
		return NULL;
	}

	float dummy;
	// internal and leaf nodes can have branch lengths
//	if( sscanf(str+*len, " : %f %n", &data->bl, &tmp) )
	if( sscanf(str+*len, " : %f %n", &dummy, &tmp) )
		*len += tmp;

	return node;
}

/* build up a little tree by hand for testing */
void sprtest(void)
{
	struct spr_node *root;
	struct spr_tree *libstate;
#define decnode(x) struct spr_node *x = newnode(#x)
	decnode(A);
	decnode(B);
	decnode(C);
	decnode(D);
	decnode(E);
	decnode(F);
	decnode(G);
// do both directions
#define link(p,rel,c) p->rel = c; c->parent=p;
	link(A,left,B);
	link(A,right,D);

	link(B,left,C);
	link(B,right,G);
	link(D,left,E);
	link(D,right,F);
	// constructor inits other pointers to NULL.

	newickprint(A);
	puts("\n");
	newickprint(D);
	treeprint(D);
	
	if (NULL == (libstate = spr_init( A, NULL, 0 ))){
		fputs("couldn't init libspr\n", stderr );
		exit( 2 );
	}
//	spr_libsprtest( libstate );

	spr( libstate, C, F );
	newickprint( A );
	root = spr_findroot(A);
	newickprint( root );
	spr_statefree( libstate );
	spr_treefree( root, TRUE );
}
char *usage=
"usage: brontler '(((a,b),(d,(e,f))),(c,d))' - do a fixed SPR\n"
"brontler '(a,(c,b))' foo\t- dump all unique SPR\n"
"brontler '(a,(c,b))' a c\t- SPR from a to c.  (internal nodes have capital letter names...)\n"
"debugging output is always printed before doing anything\n";

int main (int argc, char *argv[])
{
	struct spr_tree *tree;
	struct spr_node *root;

	int i, tmp;
	
//	srand( time(NULL) );
	srand( 42 );

	// sprtest();

	if(argc<2){
		fputs(usage,stderr);
		return 1;
	}
	// parse tree from cmdline and print it out
	root = parsenewick( argv[1], &tmp ); // assert (tmp == strlen)...
	assert( root == spr_findroot(root) );
	treeprint( root );
	newickprint( root );

	if (NULL == (tree = spr_init( root, NULL, FALSE ))){
		fputs("couldn't init libspr\n", stderr );
		return 2;
	}


	if (argc != 3){
		puts("doing SPR...");
		if( argc < 2 ){
			tmp = spr( tree, root->left->left, root->right->left );
		}else if (argc == 4){
			tmp = spr( tree, spr_searchbyname(root, argv[2]), 
			     spr_searchbyname(root, argv[3]));
			// This is lame because internal node names are weird
		}
		// newickprint( root );
		puts(tmp ? "SPR succeeded" : "SPR failed");
		treeprint( tree->root );
		newickprint( tree->root );
	}

	if (argc == 3){
		i=0;
		printf ("tree: taxa: %d, nodes: %d, possible SPRs <= %d\n",
			tree->taxa, tree->nodes, tree->lcg.m );
		tree->lcg.startstate=tree->lcg.state=16;
		lcg(&tree->lcg);
		while ( (tmp = spr_next_spr(tree)) ){
//			treeprint( tree->root );
			printf("%d: tree %d: ", ++i, tmp );
			newickprint( tree->root );
		}
	}

	spr_statefree( tree );
	spr_treefree( root, TRUE );
	spr_staticfree();
	return 0;
}
