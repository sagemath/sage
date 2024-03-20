# sage_setup: distribution = sagemath-graphs
from sage.graphs.generic_graph_pyx cimport GenericGraph_pyx

cdef class TreelengthConnected:
    cdef unsigned int n
    cdef unsigned short * c_distances
    cdef unsigned short ** distances
    cdef unsigned int diameter
    cdef str name
    cdef dict perm_inv
    cdef bint certificate
    cdef bint k_is_defined
    cdef unsigned int k
    cdef GenericGraph_pyx tree  # The final tree decomposition is stored
    cdef unsigned int length
    cdef bint leq_k
    cdef bint _treelength(self, g, k) noexcept
