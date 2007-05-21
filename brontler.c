/* generate all SPRs of a given tree.
 * 2004/11/23.  Peter Cordes <peter@cordes.ca>
 * license: GPLv2 or later
 *
 * about the name: brontler appeared during a late-night game of Scrabble :)
 */

#define _GNU_SOURCE
#define _XOPEN_SOURCE 600
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
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
// actually just use whatever the library was compiled with
// procov's data->name field is an array in the struct, not a char *.
// this necessitates some alternate code (strcpy instead of strdup)... see the #ifs below
#define SPR_PRIVATE
#include <spr.h>

// globals
int debug = 1;

// TODO: option to control printing the starting tree?
const char *usage=
"usage: brontler [options] tree [src dest]\n"
" brontler '(((a,b),(g,(e,f))),(c,d))' dump all unique SPRs\n"
" brontler '(a,(c,b))' a c\t- SPR from a to c.  (internal nodes have capital letter names...)\n"
"options: -h, -V: help and version\n"
"\t-t tree\tread tree from a file instead of the command line\n"
"\t-d number\tdebug/verbosity level (default 0)\n"
"\t-D number\tallspr library debug/verbosity level (default 0)\n"
"\t-m mode\t0: just exhaust SPRs from the starting tree. (default)\n"
"\t  1: exhaust SPRs from the starting tree, then start from the last SPR.\n"
"\t  2: take the SPRed topology as a new start point 50% of the time.\n"
"\t-T n\twhen mode>0, don't start a new tree after n unique topologies. (default 0, unlimited)\n"
"\tboth non-zero modes only stop when no non-duplicate SPRs can be done.\n";

const char *version="brontler v2.0. allspr library version " ALLSPR_VERSION "\n";



/* a "constructor".  returns a malloc()ed spr_node */
static struct spr_node *newnode(char *name)
{
	SPR_NODE_DATAPTR_TYPE *d = xcalloc(1, sizeof(*d));
#ifdef SPR_PROCOV_DATA
	if (name) strcpy (d->name, name);
#else
	d->name=name;
	d->dirty=0; d->bl=0; d->dna=NULL;
#endif // procov

	return spr_newnode( NULL, NULL, NULL, d);
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
	struct spr_node *node = newnode(NULL);
	SPR_NODE_DATAPTR_TYPE *data = node->data;

	if (*str == '('){  // internal node
#ifdef SPR_PROCOV_DATA
		strcpy (data->name, nextname);
#else
		data->name = strdup(nextname);
#endif // procov
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

		// handle unrooted trees by creating a root node here.
		if(',' == str[*len]){
			struct spr_node *newroot = newnode("root");
			newroot->left = node; node->parent = newroot;
			newroot->right = parsenewick(str + ++*len, &tmp);
			*len += tmp;
			newroot->right->parent = newroot;
			node = newroot;
		}

		if(')' == str[*len]) ++*len;
		else fprintf(stderr,"no close paren following right subtree: \"%s\"\n", str+*len);

	}else if (strspn(str, " -_.abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")){
		// leaf node
		tmp = 0;  // let scanf handle taxon names
#ifdef SPR_PROCOV_DATA
		if ( 1 > sscanf( str, " %[^][)(:;, \f\n\r\t\v] %n", data->name, len))
#else
		if ( 1 > sscanf( str, " %a[^][)(:;, \f\n\r\t\v] %n", &data->name, len))
#endif // procov
			fprintf(stderr, "bad newick element: \"%s\"\n",str);
	}else{ // not a node!?
		fprintf(stderr, "bad newick element, not an open paren or a taxon name: \"%s\"\n",str);
		free( node );
		return NULL;
	}

	// internal and leaf nodes can have branch lengths, or bootstrap:branchlen
	*len += strspn(str+*len, "0123456789.: \f\n\r\t\v");

	return node;
}

/* build up a little tree by hand for testing */
static void sprtest(void)
{
	struct spr_node *root;
	struct spr_tree *libstate;
#define decnode(x) struct spr_node *x = newnode(#x)
	decnode(A);	decnode(B);	decnode(C);	decnode(D);
	decnode(E);	decnode(F);	decnode(G);
// do both directions
#define link(p,rel,c) p->rel = c; c->parent=p;
	link(A,left,B);	link(A,right,D);

	link(B,left,C);	link(B,right,G);
	link(D,left,E);	link(D,right,F);
	// constructor inits other pointers to NULL.

	newickprint(A, stdout);
	puts("\n");
	newickprint(D, stdout);
	treeprint(D, stderr);
	
	if (NULL == (libstate = spr_init( A, NULL, 0 ))){
		fputs("couldn't init libspr\n", stderr );
		exit( 2 );
	}
//	spr_libsprtest( libstate );

	spr(libstate, C, F);
	newickprint(A, stdout);
	root = spr_findroot(A);
	newickprint(root, stdout);
	spr_statefree(libstate);
	spr_treefree(root, TRUE);
}

// return a malloc()ed buffer holding the entire contents of the file, nul terminated.
char *readfile(char *file)
{
	int offset = 0, size = 4000, tmp;
	char *buf = xmalloc(size);
	FILE *f = fopen(file, "r");
	if (!f){ perror("brontler: error opening tree file:" ); exit(2); }
	while ((tmp=fread(buf+offset, 1, size, f))  > 0){
		offset += tmp;
		size*=2;
		buf=xrealloc(buf, size*2); // yes, 4 times as big as first iteration
		// after the first time through the loop, read goes from 1/4 to 3/4 of buf
	}
	if (!feof(f) || ferror(f) || fclose(f)){ perror("brontler: error reading tree file:" ); exit(2); }
	buf[offset] = '\0';
	buf = xrealloc(buf, strlen(buf)+1);
	return buf;
}

// This is where the action is:
// enumerate the possible SPRs, one per line with various counters.
// see usage string for meaning of mode.
static int allspr(struct spr_tree *sprtree, int spr_mode, long topolimit)
{
	int treeiter, bestspr, treecount, sprnum;
	int oldtreecount=0, tmp;
	printf ("tree: taxa: %d, nodes: %d, possible SPRs <= %d\n",
		sprtree->taxa, sprtree->nodes, sprtree->lcg.m );
	// tree->lcg.state = 16;

	for(treecount=0, treeiter=1 ; ; treeiter++){
		bestspr = 0;
		while ( (sprnum = spr_next_spr(sprtree)) ){
			++treecount;
			if (debug>=4) spr_treedump(sprtree, stderr);
			if (debug != 3){ // in case you want just #trees/iteration
				printf("%d: tree %d.%d: ", treecount, treeiter, sprnum);
				newickprint(sprtree->root, stdout);
			}
			bestspr = sprnum;
			if (spr_mode==2 && rand()%2) break;
		}

		if (debug>=1){
			printf("tree iteration %d gave %d new trees\n", treeiter, treecount-oldtreecount);
			oldtreecount = treecount;
		}

		if (spr_mode > 0 && (!topolimit || treecount < topolimit) && bestspr){
			tmp = spr_apply_sprnum(sprtree, bestspr);
			assert ( tmp /* spr_apply_sprnum should always succeed */ );
		}else break;
	}
	return TRUE;
}

int main (int argc, char *argv[])
{
	struct spr_tree *sprtree;
	struct spr_node *root, *src, *dest;
	char *treestring = NULL;
	int spr_mode=0, topolimit=0;
	int i, tmp, retval=0;
	
//	srand( time(NULL) );
	srand( 42 );

	opterr = 1; // make getopt print specific error messages for us
	while ((i = getopt (argc, argv, "hVD:d:m:t:T:")) != -1){
	  switch(i){
	  case 'h': puts(usage);   return 0;
	  case 'V': puts(version); return 0;
	  case 'd': debug=atoi(optarg); break;
	  case 'D': spr_setdebug(atoi(optarg)); break;
	  case 'm': spr_mode=atoi(optarg); break;
	  case 't': treestring=readfile(optarg); break;
	  case 'T': topolimit=atoi(optarg); break;
	  case '?':
		  fputs("you need -h (help)\n", stderr);
		  return 1;
	  default:
		  abort ();	// this should never happen
	  }
	}

	if (42 == debug) sprtest();

	if (!treestring){ // take the tree from the command line
		if ((argc - optind) < 1){
			fputs("brontler: error: no tree specified\n", stderr);
			fputs(usage,stderr);
			return 1;
		}else
			treestring = argv[optind++];
	}

	// parse tree and print it out
	root = parsenewick(treestring, &tmp); // assert (tmp == strlen)...
	assert( root == spr_findroot(root) );
	if (debug>=2){
		puts("starting tree:");
		if (debug>=5) treeprint(root, stderr);
		newickprint(root, stdout);
	}

	if (NULL == (sprtree = spr_init(root, NULL, FALSE))){
		fputs("couldn't init libspr\n", stderr);
		return 2;
	}
	if (debug>=6) spr_treedump(sprtree, stderr);

	switch (argc - optind){
	case 0: retval = !allspr(sprtree, spr_mode, topolimit); break;
	case 2:
		src  = spr_treesearchbyname(sprtree, argv[argc-optind]);
		dest = spr_treesearchbyname(sprtree, argv[argc-optind+1]);
		if (debug>=2) puts("doing SPR...");
		tmp = spr(sprtree, src, dest);
		if (debug>=1){
			puts(tmp ? "SPR succeeded" : "SPR failed");
			if (debug>=4) spr_treedump(sprtree, stderr);
		}
		newickprint(sprtree->root, stdout);
		retval = !tmp;
		break;
	default:
		fputs("wrong number of non-option arguments. you need -h (help)\n", stderr);
		return 1;
	}

	spr_statefree(sprtree);
	spr_treefree(root, TRUE);
	spr_staticfree();
	return retval;
}
