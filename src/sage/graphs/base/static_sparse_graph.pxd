from cpython.object cimport PyObject
from libc.stdint cimport uint32_t, uint64_t, UINT32_MAX
from sage.data_structures.bitset cimport *

ctypedef unsigned short ushort
ctypedef unsigned int uint

cdef extern from "stdlib.h":
    ctypedef void const_void "const void"
    void qsort(void *base, int nmemb, int size,
               int(*compar)(const_void *, const_void *)) nogil

    void *bsearch(const_void *key, const_void *base, size_t nmemb,
                  size_t size, int(*compar)(const_void *, const_void *)) nogil

cdef extern from "search.h":
    void *lfind(const_void *key, const_void *base, size_t *nmemb,
                  size_t size, int(*compar)(const_void *, const_void *)) nogil

ctypedef struct short_digraph_s:
    uint32_t * edges
    uint32_t ** neighbors
    PyObject * edge_labels
    int m
    int n
    bint sorted_neighbors

ctypedef short_digraph_s short_digraph[1]

cdef int init_short_digraph(short_digraph g, G, edge_labelled=?, vertex_list=?, sort_neighbors=?) except -1
cdef void free_short_digraph(short_digraph g) noexcept
cdef int init_reverse(short_digraph dst, short_digraph src) except -1
cdef int out_degree(short_digraph g, int u) noexcept
cdef uint32_t * has_edge(short_digraph g, int u, int v) noexcept
cdef object edge_label(short_digraph g, uint32_t * edge)
cdef int tarjan_strongly_connected_components_C(short_digraph g, int *scc) noexcept
cdef void strongly_connected_components_digraph_C(short_digraph g, int nscc, int *scc, short_digraph output) noexcept

cdef uint32_t simple_BFS(short_digraph g,
                         uint32_t source,
                         uint32_t *distances,
                         uint32_t *predecessors,
                         uint32_t *waiting_list,
                         bitset_t seen) noexcept
