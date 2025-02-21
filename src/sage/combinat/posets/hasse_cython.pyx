# cython: binding=True
r"""
Some fast computations for finite posets
"""
# ****************************************************************************
#       Copyright (C) 2020 Frédéric Chapoton <chapoton@unistra.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.lazy_import import LazyImport
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest

coxeter_matrix_fast = LazyImport('sage.combinat.posets.hasse_cython_flint', 'coxeter_matrix_fast',
                                 deprecation=35741)
moebius_matrix_fast = LazyImport('sage.combinat.posets.hasse_cython_flint', 'moebius_matrix_fast',
                                 deprecation=35741)


class IncreasingChains(RecursivelyEnumeratedSet_forest):
    r"""
    The enumerated set of increasing chains.

    INPUT:

    - ``positions`` -- list of sets of integers describing the poset,
      as given by the lazy attribute ``_leq_storage`` of Hasse diagrams

    - ``element_constructor`` -- used to determine the type of chains,
      for example :class:`list` or :class:`tuple`

    - ``exclude`` -- list of integers that should not belong to the chains

    - ``conversion`` -- (optional) list of elements of the poset

    If ``conversion`` is provided, it is used to convert chain elements
    to elements of this list.

    EXAMPLES::

        sage: from sage.combinat.posets.hasse_cython import IncreasingChains
        sage: D = IncreasingChains([{0,1},{1}], list, []); D
        An enumerated set with a forest structure
        sage: D.cardinality()
        4
        sage: list(D)
        [[], [0], [0, 1], [1]]
    """
    def __init__(self, list positions, element_constructor,
                 list exclude, conversion=None):
        """
        The enumerated set of increasing chains.

        TESTS::

            sage: from sage.combinat.posets.hasse_cython import IncreasingChains
            sage: D = IncreasingChains([{0,1},{1}], list, [])
            sage: TestSuite(D).run(skip='_test_pickling')
        """
        cdef int i
        cdef Py_ssize_t n = len(positions)
        self._n = n
        if exclude is not None:
            self._greater_than = [gti.difference(exclude) for gti in positions]
            self._vertices = [u for u in range(n) if u not in exclude]
        else:
            self._greater_than = positions
            self._vertices = list(range(n))

        self._constructor = element_constructor
        self._conversion = conversion
        if conversion is not None:
            self._from_poset = {elt: i for i, elt in enumerate(conversion)}

        self._roots = (tuple(),)

        RecursivelyEnumeratedSet_forest.__init__(self, algorithm='depth',
                    category=FiniteEnumeratedSets())

    def __contains__(self, tup):
        """
        Membership testing.

        If ``conversion`` was provided, it first converts elements of ``tup``
        to integers.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_cython import IncreasingChains
            sage: D = IncreasingChains([{0,1},{1}], list, [])
            sage: all(x in D for x in D)
            True
            sage: [2] in D
            False
            sage: [1,1] in D
            False

            sage: P = Poset({'a':['b'],'b':[]})
            sage: ['a'] in P.chains()
            True

        TESTS::

            sage: from sage.combinat.posets.hasse_cython import IncreasingChains
            sage: D = IncreasingChains([{0,1},{1}], list, [])
            sage: all(tuple(x) in D for x in D)
            True
        """
        cdef int k
        cdef Py_ssize_t i, x, y
        if not tup:
            return True
        if self._conversion is not None:
            tup = [self._from_poset[elt] for elt in tup]
        if any(not(0 <= i < self._n) for i in tup):
            return False
        y = tup[0]
        for k in range(1, len(tup)):
            x = y
            y = tup[k]
            if x == y or y not in self._greater_than[x]:
                return False
        return True

    def post_process(self, chain):
        """
        Create a chain from the internal object.

        If ``conversion`` was provided, it first converts elements of the
        chain to elements of this list.

        Then the given ``element_constructor`` is applied to the chain.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_cython import IncreasingChains
            sage: D = IncreasingChains([{0,1},{1}], list, [])
            sage: D.post_process((0,1))
            [0, 1]

            sage: P = Poset({'a':['b'],'b':[]})
            sage: list(P.chains())
            [[], ['a'], ['a', 'b'], ['b']]
        """
        cdef Py_ssize_t i
        if self._conversion is not None:
            return self._constructor(self._conversion[i] for i in chain)
        return self._constructor(chain)

    def children(self, chain):
        """
        Return the children of a chain, by adding one largest element.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_cython import IncreasingChains
            sage: D = IncreasingChains([{0,1},{1}], list, [])
            sage: D.children((0,))
            [(0, 1)]

            sage: P = Poset({'a':['b'],'b':[]})
            sage: next(iter(P.chains()))
            []
        """
        cdef Py_ssize_t x, y
        if not chain:
            return [(x,) for x in self._vertices]
        x = chain[-1]
        return [chain + (y,) for y in self._greater_than[x] if x != y]
