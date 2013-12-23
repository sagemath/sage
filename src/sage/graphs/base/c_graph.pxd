#**************************************************************************
#        Copyright (C) 2008-9 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#**************************************************************************

from sage.misc.bitset cimport bitset_t

cdef class CGraph:
    cdef int num_verts
    cdef int num_arcs
    cdef int *in_degrees
    cdef int *out_degrees

    cdef int add_arc_unsafe(self, int, int)
    cdef int has_arc_unsafe(self, int, int)
    cdef int del_arc_unsafe(self, int, int)
    cdef int out_neighbors_unsafe(self, int, int *, int)
    cdef int in_neighbors_unsafe(self, int, int *, int)

    cdef bitset_t active_vertices

    cpdef bint has_vertex(self, int n)
    cpdef check_vertex(self, int n)
    cpdef del_vertex(self, int v)
    cpdef int current_allocation(self)
    cpdef add_arc(self, int u, int v)
    cpdef bint has_arc(self, int u, int v)
    cpdef del_all_arcs(self, int u, int v)
    cdef adjacency_sequence_in(self, int n, int *vertices, int v, int* sequence)
    cdef adjacency_sequence_out(self, int n, int *vertices, int v, int* sequence)
    cpdef list all_arcs(self, int u, int v)
    cpdef list in_neighbors(self, int v)
    cpdef list out_neighbors(self, int u)
    cpdef list verts(self)
    cpdef add_vertices(self, object verts)
    cdef int del_vertex_unsafe(self, int)
    cpdef realloc(self, int)
    cdef int add_vertex_unsafe(self, int)

cdef int get_vertex(object u, dict vertex_ints, dict vertex_labels, CGraph G) except ? -2
cdef object vertex_label(int u_int, dict vertex_ints, dict vertex_labels, CGraph G)
cdef int check_vertex(object u, dict vertex_ints, dict vertex_labels, CGraph G, CGraph G_revx, bint reverse) except ? -1
# TODO: edge functions!


