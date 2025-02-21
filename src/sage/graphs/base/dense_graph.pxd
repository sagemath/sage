#*****************************************************************************
#       Copyright (C) 2008-2009 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.graphs.base.c_graph cimport CGraph, CGraphBackend
from sage.data_structures.binary_matrix cimport binary_matrix_t

cdef class DenseGraph(CGraph):
    cdef bint _directed
    cdef binary_matrix_t edges
    cdef inline int _add_arc_unsafe(self, int, int) except -1
    cdef inline int _del_arc_unsafe(self, int u, int v) except -1

cdef int copy_dense_graph(DenseGraph dest, DenseGraph src) except -1

cdef class DenseGraphBackend(CGraphBackend):
    cdef DenseGraph _cg
    cdef inline CGraph cg(self):
        return <CGraph> self._cg
