
#*****************************************************************************
#      Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************


include 'sage/ext/cdefs.pxi'
include 'data_structures_pxd.pxi' # includes bitsets


from sage.rings.integer cimport Integer
from automorphism_group_canonical_label cimport \
    get_aut_gp_and_can_lab, aut_gp_and_can_lab, agcl_work_space, \
    allocate_agcl_output, deallocate_agcl_output, \
    allocate_agcl_work_space, deallocate_agcl_work_space
from double_coset cimport double_coset


cdef class PythonPartitionStack:
    cdef PartitionStack *c_ps

