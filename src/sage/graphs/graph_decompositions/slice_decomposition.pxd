from libcpp.vector cimport vector

from sage.structure.sage_object cimport SageObject
from sage.graphs.base.c_graph cimport CGraph

cdef void extended_lex_BFS(
        CGraph cg, vector[int] &sigma, vector[int] *sigma_inv,
        int initial_v_int, vector[int] *pred, vector[size_t] *xslice_len,
        vector[vector[int]] *lex_label) except *

cdef class SliceDecomposition(SageObject):
    cdef tuple sigma
    cdef dict sigma_inv
    cdef vector[size_t] xslice_len
    cdef dict lex_label
    cdef object _graph_class
    cdef object _underlying_graph
