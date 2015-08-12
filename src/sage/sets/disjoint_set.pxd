#*****************************************************************************
#      Copyright (C) 2009 Sebastien Labbe <slabqc at gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/groups/perm_gps/partn_ref/data_structures_pxd.pxi'

from sage.structure.sage_object cimport SageObject

cdef class DisjointSet_class(SageObject):
    cdef OrbitPartition *_nodes

cdef class DisjointSet_of_integers(DisjointSet_class):
    pass

cdef class DisjointSet_of_hashables(DisjointSet_class):
    cdef list _int_to_el
    cdef dict _el_to_int
    cdef DisjointSet_of_integers _d

