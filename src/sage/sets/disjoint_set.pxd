#*****************************************************************************
#       Copyright (C) 2009 Sebastien Labbe <slabqc at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.perm_gps.partn_ref.data_structures cimport OrbitPartition
from sage.structure.sage_object cimport SageObject

cpdef DisjointSet(arg)

cdef class DisjointSet_class(SageObject):
    cdef OrbitPartition *_nodes
    cpdef cardinality(self)
    cpdef number_of_subsets(self)

cdef class DisjointSet_of_integers(DisjointSet_class):
    cpdef int find(self, int i) noexcept
    cpdef void union(self, int i, int j) noexcept
    cpdef root_to_elements_dict(self)
    cpdef element_to_root_dict(self)
    cpdef to_digraph(self)

cdef class DisjointSet_of_hashables(DisjointSet_class):
    cdef list _int_to_el
    cdef dict _el_to_int
    cpdef find(self, e)
    cpdef void union(self, e, f)
    cpdef root_to_elements_dict(self) noexcept
    cpdef element_to_root_dict(self)
    cpdef to_digraph(self)
