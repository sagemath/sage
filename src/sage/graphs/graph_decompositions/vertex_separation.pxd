# sage_setup: distribution = sagemath-graphs
from libc.stdint cimport uint8_t
from sage.graphs.graph_decompositions.fast_digraph cimport FastDigraph

cdef list find_order(FastDigraph, uint8_t *, int)
