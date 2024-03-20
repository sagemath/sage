# sage_setup: distribution = sagemath-graphs
# distutils: libraries = cliquer

cdef extern from "cliquer/graph.h":
    ctypedef struct graph_t

    cdef graph_t * graph_new(int n)
    cdef void graph_print(graph_t *g)
    cdef void graph_free(graph_t *g)
    cdef void GRAPH_ADD_EDGE(graph_t *g, int, int)

cdef extern from "cliquer/cliquer.h":
    struct clique_options

    bint clique_print_time(intlevel, int i, int n, int max, double cputime, double realtime, clique_options *opts)

cdef extern from "cliquer/reorder.h":
    cdef int *reorder_by_greedy_coloring(graph_t *g, bint weighted)
    cdef int *reorder_by_degree(graph_t *g, bint weighted)
