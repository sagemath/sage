
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef ENABLE_LONG_OPTIONS
#include <getopt.h>
#endif

#include "cliquer/cliquer.h"


#define TRYFORHELP  "Try `%s -h' for more information.\n",argv[0]

void printhelp(char *prog);
void read_options(int argc, char **argv);
void print_search(graph_t *g);
boolean record_clique_func(set_t s,graph_t *g,clique_options *opts);
boolean print_clique_func(set_t s,graph_t *g,clique_options *opts);
void print_clique(set_t s,graph_t *g);

// As the global variables remain between two SAGE call, they
// have to be reset each time
void sage_reset_global_variables(void);
// The opt structure has to be initialised in each SAGE function
clique_options * sage_init_clique_opt(void);
// Computes a maximum clique of the graph g and return its size
// The table list contains the ID of the vertices
int sage_clique_max(graph_t *g,int **list);
int sage_all_clique_max(graph_t *g,int **list);
int sage_clique_number(graph_t *g);


