from cpython cimport PyObject
from libc.stdint cimport uint32_t
from sage.misc.bitset cimport *
ctypedef unsigned short ushort
ctypedef unsigned int uint

cdef extern from "stdlib.h":
    ctypedef void const_void "const void"
    void qsort(void *base, int nmemb, int size,
            int(*compar)(const_void *, const_void *)) nogil

    void *bsearch(const_void *key, const_void *base, size_t nmemb,
                  size_t size, int(*compar)(const_void *, const_void *)) nogil

ctypedef struct short_digraph_s:
   uint32_t *  edges
   uint32_t ** neighbors
   PyObject * edge_labels
   int m
   int n

ctypedef short_digraph_s short_digraph[1]

cdef int init_short_digraph(short_digraph g, G, edge_labelled = ?) except -1
cdef void free_short_digraph(short_digraph g)
cdef int init_reverse(short_digraph dst, short_digraph src) except -1
cdef int out_degree(short_digraph g, int u)
cdef inline uint32_t * has_edge(short_digraph g, int u, int v)
cdef inline object edge_label(short_digraph g, uint32_t * edge)
