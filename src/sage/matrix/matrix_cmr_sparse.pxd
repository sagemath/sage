from sage.libs.cmr.cmr cimport CMR_CHRMAT, CMR_GRAPH, CMR_GRAPH_EDGE

from .matrix_sparse cimport Matrix_sparse

cdef class Matrix_cmr_sparse(Matrix_sparse):
    pass

# cdef class Matrix_cmr_double_sparse(Matrix_cmr_sparse):
#     pass

# cdef class Matrix_cmr_int_sparse(Matrix_cmr_sparse):
#     pass

cdef class Matrix_cmr_chr_sparse(Matrix_cmr_sparse):

    cdef CMR_CHRMAT *_mat
    cdef object _root

    cdef _init_from_dict(self, dict d, int nrows, int ncols)

cdef _sage_edge(CMR_GRAPH *graph, CMR_GRAPH_EDGE e)
cdef _sage_graph(CMR_GRAPH *graph)
