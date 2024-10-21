# sage_setup: distribution = sagemath-graphs
# ****************************************************************************
#       Copyright (C) 2008-2009 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.graphs.base.c_graph cimport CGraph, CGraphBackend
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
    cpdef bint is_directed(self) noexcept

    cdef int _del_arc_unsafe(self, int, int, SparseGraphBTNode **) except -1
    cdef int _add_arc_label_unsafe(self, int, int, int, SparseGraphBTNode **) except -1
    cdef int _del_arc_label_unsafe(self, int, int, int, SparseGraphBTNode **) noexcept
    cdef SparseGraphLLNode* arc_labels_unsafe(self, int u, int v) noexcept
    cpdef int out_degree(self, int u) noexcept
    cpdef int in_degree(self, int u) noexcept

    cdef inline int _neighbors_unsafe (self, int u, bint out, int *neighbors, int size) except -2

    cdef inline int _neighbors_BTNode_unsafe (self, int u, bint out, SparseGraphBTNode **res, int size) except -2

    cdef inline SparseGraphBTNode* next_out_neighbor_BTNode_unsafe(self, int u, int v) noexcept:
        """
        Return the next out-neighbor of ``u`` that is greater than ``v``.

        If ``v`` is ``-1`` return the first neighbor of ``u``.

        Return ``NULL`` in case there does not exist such an out-neighbor.

        .. WARNING::

            Repeated calls to this function until NULL is returned DOES NOT
            yield a linear time algorithm in the number of neighbors of u.
            To list the neighbors of a vertex in linear time, one should use
            _neighbors_BTNode_unsafe.
        """
        return self.next_neighbor_BTNode_unsafe(self.vertices, u, v)

    cdef inline SparseGraphBTNode* next_in_neighbor_BTNode_unsafe(self, int v, int u) noexcept:
        """
        Return the next in-neighbor of ``v`` that is greater than ``u``.

        If ``u`` is ``-1`` return the first neighbor of ``v``.

        Return ``NULL`` in case there does not exist such an in-neighbor.

        .. WARNING::

            Repeated calls to this function until NULL is returned DOES NOT
            yield a linear time algorithm in the number of neighbors of u.
            To list the neighbors of a vertex in linear time, one should use
            _neighbors_BTNode_unsafe.
        """
        return self.next_neighbor_BTNode_unsafe(self.vertices_rev, v, u)

    cdef inline SparseGraphBTNode* next_neighbor_BTNode_unsafe(self, SparseGraphBTNode** vertices, int u, int v) noexcept


cdef class SparseGraphBackend(CGraphBackend):
    cdef int edge_labels_max
    cdef list edge_labels_available_ids
    cdef SparseGraph _cg
    cdef inline CGraph cg(self):
        return <CGraph> self._cg
