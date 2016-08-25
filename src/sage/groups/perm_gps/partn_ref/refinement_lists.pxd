#*****************************************************************************
#      Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#      Copyright (C) 2009 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include 'data_structures_pxd.pxi' # includes bitsets


# name of the three functions to customize
cdef int refine_list(PartitionStack *, void *, int *, int)
cdef int compare_lists(int *, int *, void *, void *, int)
cdef bint all_list_children_are_equivalent(PartitionStack *, void *)
