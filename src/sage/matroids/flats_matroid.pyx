r"""
Flats matroids

Matroids are characterized by a set of flats, which are sets invariant under
closure. The ``FlatsMatroid`` class implements matroids using this information
as data.

A ``FlatsMatroid`` can be created from another matroid or from a dictionary of
flats. For a full description of allowed inputs, see
:class:`below <sage.matroids.flats_matroid.FlatsMatroid>`. It is
recommended to use the :class:`Matroid() <sage.matroids.matroid.Matroid>`
class for a more flexible way of constructing a ``FlatsMatroid`` and other
classes of matroids. For direct access to the ``FlatsMatroid`` constructor,
run::

    sage: from sage.matroids.flats_matroid import FlatsMatroid

AUTHORS:

- Giorgos Mousa (2024-01-01): initial version
"""

# ****************************************************************************
#       Copyright (C) 2024 Giorgos Mousa <gmousa@proton.me>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport Py_EQ, Py_NE
from sage.structure.richcmp cimport rich_to_bool, richcmp
from sage.matroids.matroid cimport Matroid
from sage.matroids.set_system cimport SetSystem
from sage.combinat.posets.lattices import LatticePoset, FiniteLatticePoset

cdef class FlatsMatroid(Matroid):
    r"""
    INPUT:

    - ``M`` -- matroid (default: ``None``)
    - ``groundset`` -- list (default: ``None``); the groundset of the matroid
    - ``flats`` -- (default: ``None``) the dictionary of the lists of flats
      (indexed by their rank), or the list of all flats, or the lattice of
      flats of the matroid

    .. NOTE::

        For a more flexible means of input, use the ``Matroid()`` function.
    """

    # necessary (__init__, groundset, _rank)

    def __init__(self, M=None, groundset=None, flats=None):
        """
        Initialization of the matroid. See class docstring for full
        documentation.

        TESTS::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.catalog.Fano())
            sage: TestSuite(M).run()
        """
        self._F = {}
        self._L = None
        if M is not None:
            self._groundset = M.groundset()
            for i in range(M.rank() + 1):
                for F in M.flats(i):
                    try:
                        self._F[i].add(F)
                    except KeyError:
                        self._F[i] = set()
                        self._F[i].add(F)
        else:
            self._groundset = frozenset(groundset)
            if isinstance(flats, dict):
                for i in sorted(flats):
                    for F in flats[i]:
                        try:
                            self._F[i].add(frozenset(F))
                        except KeyError:
                            self._F[i] = set()
                            self._F[i].add(frozenset(F))
            else:  # store lattice of flats
                if isinstance(flats, FiniteLatticePoset):
                    self._L = flats
                else:  # assume iterable of flats
                    self._L = LatticePoset(([frozenset(F) for F in flats], lambda x, y: x < y))
                self._matroid_rank = self._L.rank()
                for i in range(self._matroid_rank + 1):
                    self._F[i] = set()
                for x in self._L:
                    self._F[self._L.rank(x)].add(x)
        self._matroid_rank = max(self._F, default=-1)

    cpdef frozenset groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT: :class:`frozenset`

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.Theta(2))
            sage: sorted(M.groundset())
            ['x0', 'x1', 'y0', 'y1']
        """
        return self._groundset

    cpdef int _rank(self, frozenset X) except? -1:
        """
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and ``X`` may be assumed to have
        the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.Theta(3))
            sage: M._rank(frozenset(['x1', 'y0', 'y2']))
            2

        TESTS::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = matroids.catalog.NonDesargues()
            sage: F = FlatsMatroid(M)
            sage: for S in powerset(M.groundset()):
            ....:     assert M.rank(S) == F.rank(S)
        """
        for i in range(self._matroid_rank + 1):
            for f in self._F[i]:
                if f >= X:
                    return i

    # optional

    cpdef full_rank(self):
        r"""
        Return the rank of the matroid.

        The *rank* of the matroid is the size of the largest independent
        subset of the groundset.

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.Theta(6))
            sage: M.full_rank()
            6
        """
        return self._matroid_rank

    cpdef frozenset _closure(self, frozenset X):
        """
        Return the closure of a set.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: :class:`frozenset`

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.catalog.Vamos())
            sage: sorted(M._closure(frozenset(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']
        """
        cdef int i
        for i in range(self._matroid_rank + 1):
            for f in self._F[i]:
                if f >= X:
                    return f

    cpdef bint _is_closed(self, frozenset X) noexcept:
        """
        Test if input is a closed set.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.catalog.Vamos())
            sage: M._is_closed(frozenset(['a', 'b', 'c', 'd']))
            True
            sage: M._is_closed(frozenset(['a', 'b', 'c', 'e']))
            False
        """
        cdef int i
        for i in self._F:
            if X in self._F[i]:
                return True
        return False

    cpdef _is_isomorphic(self, other, certificate=False):
        """
        Test if ``self`` is isomorphic to ``other``.

        INPUT:

        - ``other`` -- matroid
        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: boolean, and, if ``certificate=True``, a dictionary giving the
        isomorphism or ``None``

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = matroids.catalog.NonDesargues()
            sage: N = FlatsMatroid(M)
            sage: N._is_isomorphic(M)
            True
            sage: N._is_isomorphic(matroids.catalog.R9())
            False

        .. NOTE::

            Internal version that does no input checking.
        """
        if certificate:
            return self._is_isomorphic(other), self._isomorphism(other)
        N = FlatsMatroid(other)
        flats_self = [F for i in self._F for F in self._F[i]]
        flats_other = [F for i in N._F for F in N._F[i]]
        SS = SetSystem(self._groundset, flats_self)
        OS = SetSystem(N._groundset, flats_other)
        return SS._isomorphism(OS) is not None

    # representation

    def _repr_(self):
        """
        Return a string representation of the matroid.

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.Uniform(6, 6)); M
            Matroid of rank 6 on 6 elements with 64 flats
        """
        flats_num = sum(1 for i in self._F for F in self._F[i])
        return f'{Matroid._repr_(self)} with {flats_num} flats'

    # comparison

    def __hash__(self):
        r"""
        Return an invariant of the matroid.

        This function is called when matroids are added to a set. It is very
        desirable to override it so it can distinguish matroids on the same
        groundset, which is a very typical use case!

        .. WARNING::

            This method is linked to __richcmp__ (in Cython) and __cmp__ or
            __eq__/__ne__ (in Python). If you override one, you should
            (and in Cython: MUST) override the other!

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.catalog.Vamos())
            sage: N = FlatsMatroid(matroids.catalog.Vamos())
            sage: hash(M) == hash(N)
            True
            sage: O = FlatsMatroid(matroids.catalog.NonVamos())
            sage: hash(M) == hash(O)
            False
        """
        flats = frozenset([F for i in self._F for F in self._F[i]])
        return hash(tuple([self._groundset, flats]))

    def __richcmp__(left, right, int op):
        r"""
        Compare two matroids.

        We take a very restricted view on equality: the objects need to be of
        the exact same type (so no subclassing) and the internal data need to
        be the same. For FlatsMatroids, this means that the groundsets and the
        dictionaries of flats of the two matroids are equal.

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.catalog.Pappus())
            sage: N = FlatsMatroid(matroids.catalog.NonPappus())
            sage: M == N
            False
            sage: N = Matroid(M.bases())
            sage: M == N
            False
        """
        cdef FlatsMatroid lt, rt
        if op not in [Py_EQ, Py_NE]:
            return NotImplemented
        if type(left) is not type(right):
            return NotImplemented
        lt = <FlatsMatroid> left
        rt = <FlatsMatroid> right
        if lt.groundset() != rt.groundset():
            return rich_to_bool(op, 1)
        if lt.full_rank() != rt.full_rank():
            return rich_to_bool(op, 1)
        return richcmp(lt._F, rt._F, op)

    # copying, loading, saving

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle, (version, data))``, where ``unpickle`` is the
        name of a function that, when called with ``(version, data)``,
        produces a matroid isomorphic to ``self``. ``version`` is an integer
        (currently 0) and ``data`` is a tuple ``(E, F, name)`` where ``E`` is
        the groundset, ``F`` is the dictionary of flats, and ``name`` is a
        custom name.

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.catalog.Vamos())
            sage: M == loads(dumps(M))  # indirect doctest
            True
            sage: M.reset_name()
            sage: loads(dumps(M))
            Matroid of rank 4 on 8 elements with 79 flats
        """
        import sage.matroids.unpickling
        data = (self._groundset, self._F, self.get_custom_name())
        version = 0
        return sage.matroids.unpickling.unpickle_flats_matroid, (version, data)

    cpdef relabel(self, mapping):
        r"""
        Return an isomorphic matroid with relabeled groundset.

        The output is obtained by relabeling each element `e` by
        ``mapping[e]``, where ``mapping`` is a given injective map. If
        ``mapping[e]`` is not defined, then the identity map is assumed.

        INPUT:

        - ``mapping`` -- a Python object such that ``mapping[e]`` is the new
          label of `e`

        OUTPUT: matroid

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.catalog.RelaxedNonFano())
            sage: sorted(M.groundset())
            [0, 1, 2, 3, 4, 5, 6]
            sage: N = M.relabel({'g': 'x', 0: 'z'})  # 'g': 'x' is ignored
            sage: from sage.matroids.utilities import cmp_elements_key
            sage: sorted(N.groundset(), key=cmp_elements_key)
            [1, 2, 3, 4, 5, 6, 'z']
            sage: M.is_isomorphic(N)
            True

        TESTS::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.catalog.RelaxedNonFano())
            sage: f = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g'}
            sage: N = M.relabel(f)
            sage: for S in powerset(M.groundset()):
            ....:     assert M.rank(S) == N.rank([f[x] for x in S])
        """
        d = self._relabel_map(mapping)
        E = [d[x] for x in self._groundset]
        F = {}
        for i in self._F:
            F[i] = []
            F[i] += [[d[y] for y in x] for x in self._F[i]]
        M = FlatsMatroid(groundset=E, flats=F)
        return M

    # enumeration

    cpdef SetSystem flats(self, long k):
        r"""
        Return the flats of the matroid of specified rank.

        INPUT:

        - ``k`` -- integer

        OUTPUT: :class:`SetSystem`

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.Uniform(3, 4))
            sage: sorted(M.flats(2), key=str)
            [frozenset({0, 1}),
             frozenset({0, 2}),
             frozenset({0, 3}),
             frozenset({1, 2}),
             frozenset({1, 3}),
             frozenset({2, 3})]
        """
        if k in self._F:
            return SetSystem(self._groundset, self._F[k])
        return SetSystem(self._groundset)

    def flats_iterator(self, k):
        r"""
        Return an iterator over the flats of the matroid of specified rank.

        INPUT:

        - ``k`` -- integer

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.Uniform(3, 4))
            sage: sorted([list(F) for F in M.flats_iterator(2)])
            [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
        """
        if k in self._F:
            for F in self._F[k]:
                yield F

    def lattice_of_flats(self):
        """
        Return the lattice of flats of the matroid.

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.catalog.Fano())
            sage: M.lattice_of_flats()
            Finite lattice containing 16 elements
        """
        if self._L is None:
            flats = [F for i in range(self._matroid_rank + 1) for F in self._F[i]]
            self._L = LatticePoset((flats, lambda x, y: x < y))
        return self._L

    cpdef list whitney_numbers(self):
        r"""
        Return the Whitney numbers of the first kind of the matroid.

        The Whitney numbers of the first kind -- here encoded as a vector
        `(w_0=1, \ldots, w_r)` -- are numbers of alternating sign, where `w_i`
        is the value of the coefficient of the `(r-i)`-th degree term of the
        matroid's characteristic polynomial. Moreover, `|w_i|` is the number of
        `(i-1)`-dimensional faces of the broken circuit complex of the matroid.

        OUTPUT: list of integers

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.catalog.BetsyRoss())
            sage: M.whitney_numbers()
            [1, -11, 35, -25]

        TESTS::

            sage: M = Matroid(flats=[[0], [0, 1]])
            sage: M.whitney_numbers()
            []
        """
        if self.loops():
            return []
        cdef list w = [0] * (self._matroid_rank + 1)
        if self._L is None:
            for S in self.no_broken_circuits_sets_iterator():
                w[len(S)] += 1
            from sage.rings.integer_ring import ZZ
            return [ZZ((-1)**i * abs_w) for (i, abs_w) in enumerate(w)]
        else:
            mu = self._L.moebius_function_matrix()
            for (i, F) in enumerate(self._L.list()):
                w[self._L.rank(F)] += mu[0, i]
            return w

    cpdef list whitney_numbers2(self):
        r"""
        Return the Whitney numbers of the second kind of the matroid.

        The Whitney numbers of the second kind are here encoded as a vector
        `(W_0, ..., W_r)`, where `W_i` is the number of flats of rank `i`, and
        `r` is the rank of the matroid.

        OUTPUT: list of integers

        EXAMPLES::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: M = FlatsMatroid(matroids.catalog.XY13())
            sage: M.whitney_numbers2()
            [1, 13, 78, 250, 394, 191, 1]

        TESTS::

            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: for M in matroids.AllMatroids(4):  # optional - matroid_database
            ....:     assert M.whitney_numbers2() == FlatsMatroid(M).whitney_numbers2()
        """
        cdef list W = []
        cdef int i
        for i in self._F:
            W.append(len(self._F[i]))
        return W

    # verification

    cpdef is_valid(self, certificate=False):
        r"""
        Test if ``self`` obeys the matroid axioms.

        For a matroid defined by its flats, we check the flats axioms.

        If the lattice of flats has already been computed, we instead perform
        the equivalent check of whether it forms a geometric lattice.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: boolean, or (boolean, dictionary)

        EXAMPLES::

            sage: Matroid(flats={0: [[]], 1: [[0], [1]], 2: [[0, 1]]}).is_valid()
            True
            sage: Matroid(flats={0: [''], 1: ['a', 'b'], 2: ['ab']}).is_valid()
            True
            sage: M = Matroid(flats={0: [[]], 1: [[0], [1]]})  # missing groundset
            sage: M.is_valid()
            False
            sage: M = Matroid(flats={0: [''],
            ....:                    1: ['0','1','2','3','4','5','6','7','8','9','a','b','c'],
            ....:                    2: ['45','46','47','4c','56','57','5c','67','6c','7c',
            ....:                        '048','149','24a','34b','059','15a','25b','358',
            ....:                        '06a','16b','268','369','07b','178','279','37a',
            ....:                        '0123c','89abc'],
            ....:                    3: ['0123456789abc']})
            sage: M.is_valid()
            True
            sage: from sage.matroids.flats_matroid import FlatsMatroid
            sage: FlatsMatroid(matroids.catalog.NonVamos()).is_valid()
            True

        If we input a list or a lattice of flats, the method checks whether the
        lattice of flats is geometric::

            sage: Matroid(flats=[[], [0], [1], [0, 1]]).is_valid()
            True
            sage: Matroid(flats=['', 'a', 'b', 'ab']).is_valid()
            True
            sage: M = Matroid(flats=['',  # missing an extension of flat ['5'] by '6'
            ....:                    '0','1','2','3','4','5','6','7','8','9','a','b','c',
            ....:                    '45','46','47','4c','57','5c','67','6c','7c',
            ....:                    '048','149','24a','34b','059','15a','25b','358',
            ....:                    '06a','16b','268','369','07b','178','279','37a',
            ....:                    '0123c','89abc',
            ....:                    '0123456789abc'])
            sage: M.is_valid(certificate=True)
            (False, {'error': 'the lattice of flats is not geometric'})
            sage: Matroid(matroids.catalog.Fano().lattice_of_flats()).is_valid()
            True

        Some invalid lists of flats are recognized before calling ``is_valid``,
        upon the attempted matroid definition::

            sage: M = Matroid(flats=[[], [0], [1]])  # missing groundset
            Traceback (most recent call last):
            ...
            ValueError: not a join-semilattice: no top element
            sage: Matroid(flats=[[0], [1], [0, 1]])  # missing an intersection
            Traceback (most recent call last):
            ...
            ValueError: not a meet-semilattice: no bottom element
            sage: Matroid(flats=[[], [0, 1], [2], [0], [1], [0, 1, 2]])
            Traceback (most recent call last):
            ...
            ValueError: the poset is not ranked

        TESTS::

            sage: Matroid(flats={0: [], 1: [[0], [1]], 2: [[0, 1]]}).is_valid(certificate=True)  # missing an intersection
            (False, {'error': 'flats dictionary has invalid ranks'})
            sage: Matroid(flats={0: [[]], 2: [[0], [1]], 3: [[0, 1]]}).is_valid(certificate=True)  # invalid ranks
            (False, {'error': 'flats dictionary has invalid ranks'})
            sage: Matroid(flats={0: [[]], 1: [[0], [1]], 2: [[0], [0, 1]]}).is_valid(certificate=True)  # duplicates
            (False, {'error': 'flats dictionary has repeated flats'})
            sage: Matroid(flats={0: [[]], 1: [[0], [1], [0, 1]]}).is_valid(certificate=True)
            (False,
             {'error': 'a single element extension of a flat must be a subset of exactly one flat',
              'flat': frozenset()})
            sage: Matroid(flats={0: [[]], 1: [[0, 1], [2]], 2: [[0], [1], [0, 1, 2]]}).is_valid(certificate=True)
            (False,
             {'error': 'the intersection of two flats must be a flat',
              'flat 1': frozenset({0, 1}),
              'flat 2': frozenset({1})})
            sage: M = Matroid(flats={0: [''],  # missing an extension of flat ['5'] by '6'
            ....:                    1: ['0','1','2','3','4','5','6','7','8','9','a','b','c'],
            ....:                    2: ['45','46','47','4c','57','5c','67','6c','7c',
            ....:                        '048','149','24a','34b','059','15a','25b','358',
            ....:                        '06a','16b','268','369','07b','178','279','37a',
            ....:                        '0123c','89abc'],
            ....:                    3: ['0123456789abc']})
            sage: M.is_valid(certificate=True)
            (False,
             {'error': 'a single element extension of a flat must be a subset of exactly one flat',
              'flat': frozenset({'...'})})
            sage: M = Matroid(flats=[[], [0], [1], [0], [0, 1]])  # duplicates are ignored
            sage: M.lattice_of_flats()
            Finite lattice containing 4 elements
            sage: M.is_valid(certificate=True)
            (True, {})
            sage: M = Matroid(flats=['',
            ....:                    '0','1','2','3','4','5','6','7','8','9','a','b','c',
            ....:                    '45','46','47','4c','56','57','5c','67','6c','7c',
            ....:                    '048','149','24a','34b','059','15a','25b','358',
            ....:                    '06a','16b','268','369','07b','178','279','37a',
            ....:                    '0123c','89abc',
            ....:                    '0123456789abc'])
            sage: M.is_valid(certificate=True)
            (True, {})
        """
        if self._L is not None:  # if the lattice of flats is available
            if certificate:
                if not self._is_closed(self._groundset):
                    return False, {"error": "the groundset must be a flat"}
                if not self._L.is_geometric():
                    return False, {"error": "the lattice of flats is not geometric"}
                return True, {}
            return self._is_closed(self._groundset) and self._L.is_geometric()

        cdef int i, j, k
        cdef frozenset F1, F2, I12
        cdef list ranks, cover, flats_lst
        cdef bint flag

        # check flats dictionary for invalid ranks and repeated flats
        ranks = list(self._F)
        if ranks != list(range(len(ranks))):
            return False if not certificate else (False, {"error": "flats dictionary has invalid ranks"})
        flats_lst = [F for i in self._F for F in self._F[i]]
        if len(flats_lst) != len(set(flats_lst)):
            return False if not certificate else (False, {"error": "flats dictionary has repeated flats"})

        # the groundset must be a flat
        if not self._is_closed(self._groundset):
            return False if not certificate else (False, {"error": "the groundset must be a flat"})

        # a single element extension of a flat must be a subset of exactly one flat
        for i in ranks[:-1]:
            for F1 in self._F[i]:
                cover = []
                for F2 in self._F[i+1]:
                    if F2 >= F1:
                        cover.extend(F1 ^ F2)
                if len(cover) != len(F1 ^ self._groundset) or set(cover) != F1 ^ self._groundset:
                    return False if not certificate else (False, {"error": "a single element extension of a flat must be a subset of exactly one flat", "flat": F1})

        # the intersection of two flats must be a flat
        for i in ranks:
            for j in ranks[i:]:
                for F1 in self._F[i]:
                    for F2 in self._F[j]:
                        flag = False
                        I12 = F1 & F2
                        for k in self._F:
                            if k <= i:
                                if I12 in self._F[k]:
                                    flag = True
                                    break
                                if flag:
                                    break
                        if not flag:
                            return False if not certificate else (False, {"error": "the intersection of two flats must be a flat", "flat 1": F1, "flat 2": F2})

        return True if not certificate else (True, {})
