# -*- coding: utf-8 -*-
"""
Performance Test for Clone Protocol

see :class:`sage.structure.list_clone.ClonableArray`

EXAMPLES::

    sage: from sage.structure.list_clone_timings import *
    sage: cmd =["",
    ....:     "e.__copy__()",
    ....:     "copy(e)",
    ....:     "e.clone()",
    ....:     "e.__class__(e.parent(), e._get_list())",
    ....:     "e.__class__(e.parent(), e[:])",
    ....:     "e.check()",
    ....:     "",
    ....:     "add1_internal(e)",
    ....:     "add1_immutable(e)",
    ....:     "add1_mutable(e)",
    ....:     "add1_with(e)",
    ....:     "",
    ....:     "cy_add1_internal(e)",
    ....:     "cy_add1_immutable(e)",
    ....:     "cy_add1_mutable(e)",
    ....:     "cy_add1_with(e)"]

Various timings using a Cython class::

    sage: size = 5
    sage: e = IncreasingArrays()(range(size))
    sage: # random
    ....: for p in cmd:
    ....:     print("{0:36} : ".format(p), end=""); timeit(p)
                                         :
    e.__copy__()                         :  625 loops, best of 3: 446 ns per loop
    copy(e)                              :  625 loops, best of 3: 1.94 µs per loop
    e.clone()                            :  625 loops, best of 3: 736 ns per loop
    e.__class__(e.parent(), e._get_list()) :  625 loops, best of 3: 1.34 µs per loop
    e.__class__(e.parent(), e[:])        :  625 loops, best of 3: 1.35 µs per loop
    e.check()                            :  625 loops, best of 3: 342 ns per loop
                                         :
    add1_internal(e)                     :  625 loops, best of 3: 3.53 µs per loop
    add1_immutable(e)                    :  625 loops, best of 3: 3.72 µs per loop
    add1_mutable(e)                      :  625 loops, best of 3: 3.42 µs per loop
    add1_with(e)                         :  625 loops, best of 3: 4.05 µs per loop
                                         :
    cy_add1_internal(e)                  :  625 loops, best of 3: 752 ns per loop
    cy_add1_immutable(e)                 :  625 loops, best of 3: 1.28 µs per loop
    cy_add1_mutable(e)                   :  625 loops, best of 3: 861 ns per loop
    cy_add1_with(e)                      :  625 loops, best of 3: 1.51 µs per loop

Various timings using a Python class::

    sage: e = IncreasingArraysPy()(range(size))
    sage: # random
    ....: for p in cmd: print("{0:36} : ".format(p), end=""); timeit(p)
                                         :
    e.__copy__()                         :  625 loops, best of 3: 869 ns per loop
    copy(e)                              :  625 loops, best of 3: 2.13 µs per loop
    e.clone()                            :  625 loops, best of 3: 1.86 µs per loop
    e.__class__(e.parent(), e._get_list()) :  625 loops, best of 3: 7.52 µs per loop
    e.__class__(e.parent(), e[:])        :  625 loops, best of 3: 7.27 µs per loop
    e.check()                            :  625 loops, best of 3: 4.02 µs per loop
                                         :
    add1_internal(e)                     :  625 loops, best of 3: 9.34 µs per loop
    add1_immutable(e)                    :  625 loops, best of 3: 9.91 µs per loop
    add1_mutable(e)                      :  625 loops, best of 3: 12.6 µs per loop
    add1_with(e)                         :  625 loops, best of 3: 15.9 µs per loop
                                         :
    cy_add1_internal(e)                  :  625 loops, best of 3: 7.13 µs per loop
    cy_add1_immutable(e)                 :  625 loops, best of 3: 6.95 µs per loop
    cy_add1_mutable(e)                   :  625 loops, best of 3: 14.1 µs per loop
    cy_add1_with(e)                      :  625 loops, best of 3: 17.5 µs per loop
"""
#*****************************************************************************
#  Copyright (C) 2009-2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.list_clone import ClonableArray
from sage.structure.list_clone_demo import IncreasingArrays


class IncreasingArraysPy(IncreasingArrays):

    class Element(ClonableArray):
        """
        A small class for testing :class:`ClonableArray`: Increasing Lists

        TESTS::

            sage: from sage.structure.list_clone_timings import IncreasingArraysPy
            sage: TestSuite(IncreasingArraysPy()([1,2,3])).run()
        """

        def check(self):
            """
            Check that ``self`` is increasing.

            EXAMPLES::

                sage: from sage.structure.list_clone_timings import IncreasingArraysPy
                sage: IncreasingArraysPy()([1,2,3]) # indirect doctest
                [1, 2, 3]
                sage: IncreasingArraysPy()([3,2,1]) # indirect doctest
                Traceback (most recent call last):
                ...
                ValueError: Lists is not increasing
            """
            for i in range(len(self)-1):
                if self[i] >= self[i+1]:
                    raise ValueError("Lists is not increasing")


#####################################################################
######                    Timings functions                    ######
#####################################################################
def add1_internal(bla):
    """
    TESTS::

        sage: from sage.structure.list_clone_timings import *
        sage: add1_internal(IncreasingArrays()([1,4,5]))
        [2, 5, 6]
    """
    blo = bla.__copy__()
    lst = blo._get_list()
    for i in range(len(blo)):
        lst[i] += 1
    blo.set_immutable()
    blo.check()
    return blo

def add1_immutable(bla):
    """
    TESTS::

        sage: from sage.structure.list_clone_timings import *
        sage: add1_immutable(IncreasingArrays()([1,4,5]))
        [2, 5, 6]
    """
    lbla = bla[:]
    for i in range(len(lbla)):
        lbla[i] += 1
    return bla.__class__(bla.parent(), lbla)

def add1_mutable(bla):
    """
    TESTS::

        sage: from sage.structure.list_clone_timings import *
        sage: add1_mutable(IncreasingArrays()([1,4,5]))
        [2, 5, 6]
    """
    blo = bla.__copy__()
    for i in range(len(blo)):
        blo[i] += 1
    blo.set_immutable()
    blo.check()
    return blo

def add1_with(bla):
    """
    TESTS::

        sage: from sage.structure.list_clone_timings import *
        sage: add1_with(IncreasingArrays()([1,4,5]))
        [2, 5, 6]
    """
    with bla.clone() as blo:
        for i in range(len(blo)):
            blo[i] += 1
    return blo
