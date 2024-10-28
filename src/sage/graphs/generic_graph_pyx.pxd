from sage.structure.sage_object cimport SageObject
from sage.graphs.base.dense_graph cimport DenseGraph
from memory_allocator cimport MemoryAllocator

ctypedef int * D_TWO
ctypedef char * D_THREE
ctypedef fused dimension_t:
    D_TWO
    D_THREE

cdef run_spring(int, dimension_t, double*, int*, int, int, bint)

cdef class GenericGraph_pyx(SageObject):
    pass


cdef class SubgraphSearch:
    cdef int ng
    cdef int nh
    cdef (bint) (*is_admissible) (int, int *, int *) noexcept
    cdef DenseGraph g
    cdef DenseGraph h
    cdef int *busy
    cdef int *stack
    cdef int active
    cdef int *vertices
    cdef int **line_h_out
    cdef int **line_h_in
    cdef list g_vertices
    cdef int i
    cdef bint directed
    cdef int * tmp_array
    cdef MemoryAllocator mem
