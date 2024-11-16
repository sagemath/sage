from sage.matroids.matroid cimport Matroid
from sage.graphs.generic_graph_pyx cimport GenericGraph_pyx

cdef class GraphicMatroid(Matroid):
    cdef frozenset _groundset
    cdef readonly GenericGraph_pyx _G
    cdef dict _vertex_map
    cdef dict _groundset_edge_map
    cpdef frozenset groundset(self)
    cpdef int _rank(self, frozenset X) except? -1
    cpdef _vertex_stars(self)
    cpdef _minor(self, contractions, deletions)
    cpdef _has_minor(self, N, bint certificate=*)
    cpdef int _corank(self, frozenset X) noexcept
    cpdef bint _is_circuit(self, frozenset X) noexcept
    cpdef frozenset _closure(self, frozenset X)
    cpdef frozenset _max_independent(self, frozenset X)
    cpdef frozenset _max_coindependent(self, frozenset X)
    cpdef frozenset _circuit(self, frozenset X)
    cpdef frozenset _coclosure(self, frozenset X)
    cpdef bint _is_closed(self, frozenset X) noexcept
    cpdef _is_isomorphic(self, other, certificate=*)
    cpdef _isomorphism(self, other)
    cpdef is_valid(self, certificate=*)
    cpdef bint is_graphic(self) noexcept
    cpdef bint is_regular(self) noexcept
    cpdef graph(self)
    cpdef vertex_map(self)
    cpdef list groundset_to_edges(self, X)
    cpdef _groundset_to_edges(self, X)
    cpdef subgraph_from_set(self, X)
    cpdef _subgraph_from_set(self, X)
    cpdef graphic_extension(self, u, v=*, element=*)
    cpdef graphic_coextension(self, u, v=*, X=*, element=*)
    cpdef twist(self, X)
    cpdef one_sum(self, X, u, v)
    cpdef regular_matroid(self)
    cpdef relabel(self, mapping)
