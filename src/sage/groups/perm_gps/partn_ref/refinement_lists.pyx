# sage_setup: distribution = sagemath-categories
"""
Partition backtrack functions for lists -- a simple example of using partn_ref

EXAMPLES::

    sage: import sage.groups.perm_gps.partn_ref.refinement_lists
"""

#*****************************************************************************
#       Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#       Copyright (C) 2009 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.memory cimport sig_malloc, sig_free

from sage.groups.perm_gps.partn_ref.data_structures cimport *
from sage.groups.perm_gps.partn_ref.double_coset cimport double_coset, int_cmp


def is_isomorphic(self, other):
    r"""
    Return the bijection as a permutation if two lists are isomorphic, return
    ``False`` otherwise.

    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref.refinement_lists import is_isomorphic
        sage: is_isomorphic([0,0,1],[1,0,0])
        [1, 2, 0]
    """
    cdef int i, n = len(self)
    cdef PartitionStack *part
    cdef int *output
    cdef int *ordering
    part = PS_new(n, 1)
    ordering = <int *> sig_malloc((len(self)) * sizeof(int))
    output = <int *> sig_malloc((len(self)) * sizeof(int))
    if part is NULL or ordering is NULL or output is NULL:
        PS_dealloc(part)
        sig_free(ordering)
        sig_free(output)
        raise MemoryError
    for i from 0 <= i < (len(self)):
        ordering[i] = i

    cdef bint isomorphic = double_coset(<void *> self, <void *> other, part, ordering, (len(self)), &all_list_children_are_equivalent, &refine_list, &compare_lists, NULL, NULL, output)

    PS_dealloc(part)
    sig_free(ordering)
    if isomorphic:
        output_py = [output[i] for i from 0 <= i < (len(self))]
    else:
        output_py = False
    sig_free(output)
    return output_py


cdef bint all_list_children_are_equivalent(PartitionStack *PS, void *S) noexcept:
    return 0

cdef int refine_list(PartitionStack *PS, void *S, int *cells_to_refine_by, int ctrb_len) noexcept:
    return 0

cdef int compare_lists(int *gamma_1, int *gamma_2, void *S1, void *S2, int degree) noexcept:
    r"""
    Compare two lists according to the lexicographic order.
    """
    cdef list MS1 = <list> S1
    cdef list MS2 = <list> S2
    cdef int i, j
    for i in range(degree):
        j = int_cmp(MS1[gamma_1[i]], MS2[gamma_2[i]])
        if j != 0:
            return j
    return 0
