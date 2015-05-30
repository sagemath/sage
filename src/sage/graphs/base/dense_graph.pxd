
#*******************************************************************************
#        Copyright (C) 2008-9 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

from c_graph cimport CGraph

cdef class DenseGraph(CGraph):
    cdef int radix_div_shift
    cdef int radix_mod_mask
    cdef int num_longs
    cdef unsigned long *edges
