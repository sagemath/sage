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

cdef class DenseGraph(CGraph):
    cdef int radix_div_shift
    cdef int radix_mod_mask
    cdef size_t num_longs
    cdef unsigned long *edges

cdef class DenseGraphBackend(CGraphBackend):
    cdef DenseGraph _cg
    cdef inline CGraph cg(self):
        return <CGraph> self._cg
