#**************************************************************************
#        Copyright (C) 2008-9 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#**************************************************************************

from sage.data_structures.bitset cimport bitset_t
from sage.graphs.base.graph_backends cimport GenericGraphBackend
from libc.stdint cimport uint32_t

cdef class CGraph:
    cdef size_t num_verts
    cdef size_t num_arcs
    cdef int *in_degrees
    cdef int *out_degrees
    cdef bitset_t active_vertices

    ###################################
    # Vertex Functions
    ###################################

    cpdef bint has_vertex(self, int n) except -1
    cpdef check_vertex(self, int n)
    cpdef del_vertex(self, int v)
    cpdef int current_allocation(self) noexcept
    cpdef list verts(self)
    cpdef add_vertices(self, verts)
    cdef int del_vertex_unsafe(self, int) except -1
    cpdef realloc(self, int)
    cdef int add_vertex_unsafe(self, int) except -1

    ###################################
    # Edge Functions
    ###################################

    cdef inline int add_arc_unsafe(self, int u, int v) except -1:
        return self.add_arc_label_unsafe(u, v, 0)

    cdef inline int has_arc_unsafe(self, int u, int v) except -1:
        return self.has_arc_label_unsafe(u, v, -1)

    cdef int del_arc_unsafe(self, int, int) except -1

    cpdef add_arc(self, int u, int v)
    cpdef bint has_arc(self, int u, int v) except -1
    cpdef del_all_arcs(self, int u, int v)

    ###################################
    # Labeled Edge Functions
    ###################################

    cdef int add_arc_label_unsafe(self, int, int, int) except -1
    cdef int has_arc_label_unsafe(self, int, int, int) except -1
    cdef int del_arc_label_unsafe(self, int, int, int) except -1
    cdef int arc_label_unsafe(self, int, int) except -1
    cdef int all_arcs_unsafe(self, int, int, int *, int) except -1

    cpdef int arc_label(self, int u, int v) noexcept
    cpdef list all_arcs(self, int u, int v)
    cpdef del_arc_label(self, int u, int v, int l)
    cpdef bint has_arc_label(self, int u, int v, int l) noexcept

    ###################################
    # Neighbor Functions
    ###################################

    cdef int out_neighbors_unsafe(self, int, int *, int) except -2
    cdef int in_neighbors_unsafe(self, int, int *, int) except -2

    cdef inline int _next_neighbor_unsafe(self, int v, int u, bint out, int* l) except -2:
        if out:
            return self.next_out_neighbor_unsafe(v, u, l)
        else:
            return self.next_in_neighbor_unsafe(v, u, l)

    cdef int next_out_neighbor_unsafe(self, int, int, int*) except -2
    cdef int next_in_neighbor_unsafe(self, int, int, int*) except -2
    cdef adjacency_sequence_out(self, int n, int *vertices, int v, int* sequence)
    cdef adjacency_sequence_in(self, int n, int *vertices, int v, int* sequence)
    cpdef list in_neighbors(self, int v)
    cpdef list out_neighbors(self, int u)


cdef class CGraphBackend(GenericGraphBackend):
    cdef int get_vertex(self, u) except ? -2
    cdef int get_vertex_checked(self, u) except ? -2
    cdef vertex_label(self, int u_int)
    cdef int check_labelled_vertex(self, u, bint reverse) except ? -1
    #cdef CGraph _cg  # a child class should declare this accordingly
    cdef bint _directed
    cdef dict vertex_labels
    cdef dict vertex_ints
    cdef dict edge_labels
    cdef bint _loops
    cdef bint _multiple_edges
    cdef CGraph cg(self)
    cpdef add_edge(self, object u, object v, object l, bint directed)
    cpdef del_edge(self, object u, object v, object l, bint directed)
    cdef bint _has_labeled_edge_unsafe(self, int, int, object) except -1
    cdef bint _delete_edge_before_adding(self) noexcept
    cdef int new_edge_label(self, object l) except -1
    cdef int free_edge_label(self, int l_int) except -1
    cdef int _use_edge_iterator_on_subgraph(self, CGraphBackend other, object vertices, const int modus) except -1
    cdef list _all_edge_labels(self, int u, int v, uint32_t* edge=*)
