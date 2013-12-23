
#*******************************************************************************
#        Copyright (C) 2008-9 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

from c_graph cimport CGraph
include 'sage/ext/stdsage.pxi'

cdef class DenseGraph(CGraph):
    # Values inherited from CGraph:
    # cdef int num_verts
    # cdef int num_arcs
    # cdef int *in_degrees
    # cdef int *out_degrees
    # Values specific to DenseGraph:
    cdef int radix_div_shift
    cdef int radix_mod_mask
    cdef int num_longs
    cdef unsigned long *edges

    # Method declarations inherited from CGraph:
    # cdef int add_arc_unsafe(self, int, int)
    # cdef int has_arc_unsafe(self, int, int)
    # cdef int del_arc_unsafe(self, int, int)
    # cdef int out_neighbors_unsafe(self, int, int *, int)
    # cdef int in_neighbors_unsafe(self, int, int *, int)






