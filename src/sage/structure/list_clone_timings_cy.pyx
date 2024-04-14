# sage_setup: distribution = sagemath-objects
"""
Cython Functions for Timing Clone Protocol
"""
#*****************************************************************************
#  Copyright (C) 2009-2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.list_clone cimport ClonableArray


#####################################################################
######                    Timings functions                    ######
#####################################################################
cpdef ClonableArray cy_add1_internal(ClonableArray bla):
    """
    TESTS::

        sage: from sage.structure.list_clone_timings_cy import *
        sage: from sage.structure.list_clone_timings import *
        sage: cy_add1_internal(IncreasingArrays()([1,4,5]))
        [2, 5, 6]
    """
    cdef int i
    cdef list lst
    cdef ClonableArray blo
    blo = bla.__copy__()
    lst = blo._get_list()
    for i in range(len(lst)): lst[i] += 1
    blo.set_immutable()
    blo.check()
    return blo


cpdef ClonableArray cy_add1_immutable(ClonableArray bla):
    """
    TESTS::

        sage: from sage.structure.list_clone_timings import *
        sage: from sage.structure.list_clone_timings_cy import *
        sage: cy_add1_immutable(IncreasingArrays()([1,4,5]))
        [2, 5, 6]
    """
    cdef int i
    cdef list lbla
    lbla = bla[:]
    for i in range(len(lbla)): lbla[i] += 1
    return bla.__class__(bla._parent, lbla)

cpdef ClonableArray cy_add1_mutable(ClonableArray bla):
    """
    TESTS::

        sage: from sage.structure.list_clone_timings import *
        sage: from sage.structure.list_clone_timings_cy import *
        sage: cy_add1_mutable(IncreasingArrays()([1,4,5]))
        [2, 5, 6]
    """
    cdef int i
    cdef ClonableArray blo
    blo = bla.__copy__()
    for i in range(len(blo)):
        blo._setitem(i, blo._getitem(i)+1)
    blo.set_immutable()
    blo.check()
    return blo


cpdef ClonableArray cy_add1_with(ClonableArray bla):
    """
    TESTS::

        sage: from sage.structure.list_clone_timings import *
        sage: from sage.structure.list_clone_timings_cy import *
        sage: cy_add1_with(IncreasingArrays()([1,4,5]))
        [2, 5, 6]
    """
    cdef int i
    cdef ClonableArray blo
    with bla.__copy__() as blo:
        for i in range(len(blo)):
            blo._setitem(i, blo._getitem(i)+1)
    return blo
