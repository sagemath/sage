from sage.graphs.base.static_sparse_graph cimport short_digraph

cdef maximum_cardinality_search_M_short_digraph(short_digraph sd,
                                                int initial_vertex,
                                                int* alpha,
                                                int* alpha_inv,
                                                list F,
                                                bint* X)
