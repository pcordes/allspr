/* subtree pruning-regrafting (spr) library
 * Peter Cordes <peter@cordes.ca>, Dalhousie University
 * license: GPLv2 or later
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define SPR_PRIVATE
#include "spr.h"

// half-assed newick parser is in the front-end brontler.c.  A better parser
// will go here...

struct newickprinter {
  char *string;
  int total_len;
};

/* FIXME: this function needs to be able to realloc the string
 * to handle long taxa names.  
 * And while we're at it, optionally include branch lengths */
static int newick_recurse_unsafe( char *s, const struct spr_node *p )
{
	int tmp;
	if (!p->left){ // leaf
		assert( !p->right );
		return sprintf( s, "%s", p->data->name );
	}else{			// internal
	   // single-child internal nodes make no sense in phylogenetic trees
		assert( p->left && p->right );
		*s = '('; tmp=1;
		tmp += newick_recurse_unsafe( s+tmp, p->left );
		s[tmp] = ','; tmp++;
		tmp += newick_recurse_unsafe( s+tmp, p->right );
		s[tmp] = ')'; tmp++;
		return tmp;
	}		
}
/*
struct newickprint{
	char *str;
	size_t space;
};
// not finished
static int newick_recurse(const struct spr_node *p, struct newickprint *str, int pos)
{
	int tmp;
	char *s = str.str+pos;
	if (!p->left){ // leaf
		assert( !p->right );
		return sprintf( str+pos, "%s", p->data->name );
	}else{			// internal
	   // single-child internal nodes make no sense in phylogenetic trees
		assert( p->left && p->right );
		*s = '('; tmp=1;
		tmp += newick_recurse( s+tmp, p->left );
		s[tmp] = ','; tmp++;
		tmp += newick_recurse( s+tmp, p->right );
		s[tmp] = ')'; tmp++;
		return tmp;
	}		
}
*/

/*
      A
   B    D
 ((CG),(EF)) */
/* return a malloc()ed string holding a Newick rep of the tree 
 * (without branch lengths) */
char *newick( const struct spr_node *tree )
{
	int len;
	int maxlen = 100*spr_countnodes(tree);
	char *string = xmalloc( maxlen ); // plenty of space
	len = newick_recurse_unsafe( string, tree );
	string[len++] = ';';
	string[len++] = '\0';
	assert (len <= maxlen);
	string = xrealloc( string, len );
	return string;
}

void newickprint(const struct spr_node *tree, FILE *stream)
{
	char *s = newick(tree);
	fputs(s, stream); putc('\n', stream);
	free(s);
}


void treeprint(const struct spr_node *p, FILE *stream)
{
	if (p->left) treeprint(p->left, stream);
	if (!p->data)
		fputs("node with NULL data\n", stream);
	else{
		fputs(p->data->name, stream); putc('\n', stream);
	}
	if (p->right) treeprint(p->right, stream);
}
