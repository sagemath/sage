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

cdef class DisjointSet_class(SageObject):
    cdef OrbitPartition *_nodes

cdef class DisjointSet_of_integers(DisjointSet_class):
    pass

cdef class DisjointSet_of_hashables(DisjointSet_class):
    cdef list _int_to_el
    cdef dict _el_to_int
    cdef DisjointSet_of_integers _d
