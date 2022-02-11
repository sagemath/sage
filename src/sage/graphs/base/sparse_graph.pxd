#*****************************************************************************
#       Copyright (C) 2008-2009 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .c_graph cimport CGraph, CGraphBackend
cimport cython

cdef struct SparseGraphLLNode:
    int label
    int number
    SparseGraphLLNode *next

cdef struct SparseGraphBTNode:
    int vertex
    int number
    SparseGraphLLNode *labels
    SparseGraphBTNode *left
    SparseGraphBTNode *right

@cython.final
cdef class SparseGraph(CGraph):
    cdef int hash_length
    cdef int hash_mask
    cdef SparseGraphBTNode **vertices
    cdef SparseGraphBTNode **vertices_rev
    cdef bint _directed
    cpdef bint is_directed(self)

    cdef int _del_arc_unsafe(self, int, int, SparseGraphBTNode **) except -1
    cdef int _add_arc_label_unsafe(self, int, int, int, SparseGraphBTNode **) except -1
    cdef int _del_arc_label_unsafe(self, int, int, int, SparseGraphBTNode **)
    cdef SparseGraphLLNode* arc_labels_unsafe(self, int u, int v)
    cpdef int out_degree(self, int u)
    cpdef int in_degree(self, int u)

    cdef int out_neighbors_BTNode_unsafe(self, int u, SparseGraphBTNode *** p_pointers)
    cdef int in_neighbors_BTNode_unsafe(self, int u, SparseGraphBTNode *** p_pointers)

    cdef inline SparseGraphBTNode* next_out_neighbor_BTNode_unsafe(self, int u, int v):
        """
        Return the next out-neighbor of ``u`` that is greater than ``v``.

        If ``v`` is ``-1`` return the first neighbor of ``u``.

        Return ``NULL`` in case there does not exist such an out-neighbor.
        """
        return self.next_neighbor_BTNode_unsafe(self.vertices, u, v)

    cdef inline SparseGraphBTNode* next_in_neighbor_BTNode_unsafe(self, int v, int u):
        """
        Return the next in-neighbor of ``v`` that is greater than ``u``.

        If ``u`` is ``-1`` return the first neighbor of ``v``.

        Return ``NULL`` in case there does not exist such an in-neighbor.
        """
        return self.next_neighbor_BTNode_unsafe(self.vertices_rev, v, u)

    cdef inline SparseGraphBTNode* next_neighbor_BTNode_unsafe(self, SparseGraphBTNode** vertices, int u, int v)


cdef class SparseGraphBackend(CGraphBackend):
    cdef int edge_labels_max
    cdef list edge_labels_available_ids
    cdef SparseGraph _cg
    cdef inline CGraph cg(self):
        return <CGraph> self._cg
