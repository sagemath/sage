/* sage_setup: distribution = sagemath-graphs
 */
#include "cliquer/cliquer.h"

static int maximal;
static int sage_clique_count=0;
static set_t *sage_clique_list;
static int sage_clique_list_size=0;

// As the global variables remain between two Sage calls, they
// have to be reset each time
void sage_reset_global_variables(void){
  maximal=FALSE; // reachable from read_options
  sage_clique_count=0;
  sage_clique_list_size=0;
}

static boolean sage_record_clique_func(set_t s,graph_t *g,clique_options *opts) {
	if (sage_clique_count>=sage_clique_list_size) {
		sage_clique_list=realloc(sage_clique_list,(sage_clique_list_size+512) * 
				    sizeof(set_t));
		sage_clique_list_size+=512;
	}
	sage_clique_list[sage_clique_count]=set_duplicate(s);
	sage_clique_count++;
	return TRUE;
}

static int *(*reorder)(graph_t *, boolean)=reorder_by_default;
static int quiet=0;
// The opt structure has to be initialised in each Sage function
clique_options * sage_init_clique_opt(){
  sage_reset_global_variables();
  clique_options *opts;
  quiet++;
  opts=malloc(sizeof(clique_options));
  if (quiet)
    opts->time_function=NULL;
  else
    opts->time_function=clique_print_time;
  opts->output=stderr;
  opts->reorder_function=reorder;
  opts->reorder_map=NULL;
  opts->user_function=sage_record_clique_func;
  opts->user_data=NULL;
  opts->clique_list=NULL;
  opts->clique_list_length=0;
  return opts;
}

// Computes a maximum clique of the graph g and return its size
// The table list contains the ID of the vertices
int sage_clique_max(graph_t *g,int **list){
  sage_reset_global_variables();
  quiet++;
  set_t s;
  int i,l;
  clique_options *opts = sage_init_clique_opt();
  s=clique_unweighted_find_single(g,/*min_weight*/0,
				  /*max_weight*/0,/*maximal*/TRUE,
				  opts);
  free(opts);

  // Writing the answer into a int [] to be read by Sage
  int size=set_size(s);
  *list=malloc(sizeof(int)*size);
  l=0;
  for (i=0; i<SET_MAX_SIZE(s); i++) {
    if (SET_CONTAINS(s,i)) {
      *((*list)+l)=i;
      l++;
    }
  }
  return size;
}

int sage_all_clique_max(graph_t *g,int **list){
  sage_reset_global_variables();
  quiet++;
  maximal=TRUE;
  int i,j,l;

  clique_options *opts = sage_init_clique_opt();
  clique_unweighted_find_all(g,/*min_weight*/0,/*max_weight*/0,
			     maximal,opts);
  free(opts);

  int size=set_size(sage_clique_list[0]);
  *list=malloc(sizeof(int)*(size+1)*sage_clique_count);
  l=0;

  for (j=0; j<sage_clique_count; j++) {
    for (i=0; i<SET_MAX_SIZE(sage_clique_list[j]); i++) {
      if (SET_CONTAINS(sage_clique_list[j],i)) {
        *((*list)+l)=i;
        l++;
      }
    }
    set_free(sage_clique_list[j]);
    *((*list)+l)=-1;
    l++;
  }
  return (1+size)*sage_clique_count;
}

int sage_find_all_clique(graph_t *g,int **list, int min_size, int max_size){
  sage_reset_global_variables();
  quiet++;
  maximal = FALSE;
  int i, j, l;

  clique_options *opts = sage_init_clique_opt();
  clique_unweighted_find_all(g, min_size, max_size, maximal, opts);
  free(opts);

  int size = 0;
  // All the cliques in the output set will be delimited with -1 value.
  // The size of the output with delimiters has to be computed now.
  for (i = 0; i < sage_clique_count; i++) {
    size += set_size(sage_clique_list[i]) + 1;
  }
  *list = malloc(sizeof(int) * size);
  l = 0;

  for (j = 0; j < sage_clique_count; j++) {
    for (i = 0; i < SET_MAX_SIZE(sage_clique_list[j]); i++) {
      if (SET_CONTAINS(sage_clique_list[j],i)) {
        *((*list) + l) = i;
        l++;
      }
    }
    set_free(sage_clique_list[j]);
    *((*list) +l ) = -1;
    l++;
  }
  return size;
}

int sage_clique_number(graph_t *g){
  sage_reset_global_variables();
  maximal=TRUE;
  clique_options *opts;
  opts=sage_init_clique_opt();
  int n = clique_unweighted_max_weight(g,opts);
  free(opts);
  opts = NULL;
  return n;
}
