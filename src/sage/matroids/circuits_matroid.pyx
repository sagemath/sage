r"""
Circuits matroids

Matroids are characterized by a list of circuits, which are minimal dependent
sets. The ``CircuitsMatroid`` class implements matroids using this information
as data.

A ``CircuitsMatroid`` can be created from another matroid or from a list of
circuits. For a full description of allowed inputs, see
:class:`below <sage.matroids.circuits_matroid.CircuitsMatroid>`. It is
recommended to use the :func:`Matroid() <sage.matroids.constructor.Matroid>`
function for a more flexible way of constructing a ``CircuitsMatroid`` and
other classes of matroids. For direct access to the ``CircuitsMatroid``
constructor, run::

    sage: from sage.matroids.circuits_matroid import CircuitsMatroid

AUTHORS:

- Giorgos Mousa (2023-12-23): initial version
"""

# ****************************************************************************
#       Copyright (C) 2023 Giorgos Mousa <gmousa@proton.me>
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

cdef class CircuitsMatroid(Matroid):
    r"""
    A matroid defined by its circuits.

    INPUT:

    - ``M`` -- matroid (default: ``None``)
    - ``groundset`` -- list (default: ``None``); the groundset of the matroid
    - ``circuits`` -- list (default: ``None``); the collection of circuits of
      the matroid
    - ``nsc_defined`` -- boolean (default: ``False``); whether the matroid was
      defined by its nonspanning circuits

    .. NOTE::

        For a more flexible means of input, use the ``Matroid()`` function.
    """

    # necessary (__init__, groundset, _rank)

    def __init__(self, M=None, groundset=None, circuits=None, nsc_defined=False):
        """
        Initialization of the matroid. See class docstring for full
        documentation.

        TESTS::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.catalog.Fano())
            sage: TestSuite(M).run()
        """
        if M is not None:
            self._groundset = M.groundset()
            self._C = set(M.circuits())
        else:
            self._groundset = frozenset(groundset)
            self._C = set([frozenset(C) for C in circuits])
        # k-circuits
        self._k_C = {}
        for C in self._C:
            try:
                self._k_C[len(C)].add(C)
            except KeyError:
                self._k_C[len(C)] = set()
                self._k_C[len(C)].add(C)
        self._sorted_C_lens = sorted(self._k_C)
        self._matroid_rank = self.rank(self._groundset)
        self._nsc_defined = nsc_defined

    cpdef frozenset groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT: set

        EXAMPLES::

            sage: M = matroids.Theta(2)
            sage: sorted(M.groundset())
            ['x0', 'x1', 'y0', 'y1']
        """
        return self._groundset

    cpdef int _rank(self, frozenset X) noexcept:
        """
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and ``X`` may be assumed to have
        the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface

        OUTPUT: integer

        EXAMPLES::

            sage: M = matroids.Theta(3)
            sage: M._rank(frozenset(['x1', 'y0', 'y2']))
            2
        """
        return len(self._max_independent(X))

    # optional

    cpdef full_rank(self):
        r"""
        Return the rank of the matroid.

        The *rank* of the matroid is the size of the largest independent
        subset of the groundset.

        OUTPUT: integer

        EXAMPLES::

            sage: M = matroids.Theta(20)
            sage: M.full_rank()
            20
        """
        return self._matroid_rank

    cpdef bint _is_independent(self, frozenset X) noexcept:
        """
        Test if input is independent.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: boolean

        EXAMPLES::

            sage: M = matroids.Theta(4)
            sage: M._is_independent(frozenset(['y0', 'y1', 'y3', 'x2']))
            False
            sage: M._is_independent(frozenset(['y0', 'y2', 'y3', 'x2']))
            True
        """
        cdef int i, l = len(X)
        for i in self._sorted_C_lens:
            if i > l:
                break
            for C in self._k_C[i]:
                if C <= X:
                    return False
        return True

    cpdef frozenset _max_independent(self, frozenset X):
        """
        Compute a maximal independent subset.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: frozenset; a maximal independent subset of ``X``

        EXAMPLES::

            sage: M = matroids.Theta(6)
            sage: len(M._max_independent(M.groundset()))
            6
        """
        cdef set XX = set(X)
        cdef int i
        cdef frozenset C
        while True:
            try:
                C = self._circuit(frozenset(XX))
                e = next(iter(C))
                XX.remove(e)
            except (ValueError, StopIteration):
                return frozenset(XX)

    cpdef frozenset _circuit(self, frozenset X):
        """
        Return a minimal dependent subset.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: frozenset; a circuit contained in ``X``, if it exists.
        Otherwise an error is raised.

        EXAMPLES::

            sage: M = matroids.Theta(4)
            sage: sorted(M._circuit(frozenset(['y0', 'y1', 'y3', 'x2'])))
            ['x2', 'y0', 'y1', 'y3']
            sage: M._circuit(frozenset(['y0', 'y2', 'y3', 'x2']))
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set
        """
        cdef int i, l = len(X)
        for i in self._sorted_C_lens:
            if i > l:
                break
            for C in self._k_C[i]:
                if C <= X:
                    return C
        raise ValueError("no circuit in independent set")

    cpdef frozenset _closure(self, frozenset X):
        """
        Return the closure of a set.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: :class:`frozenset`

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.catalog.Vamos())
            sage: sorted(M._closure(frozenset(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']
        """
        cdef set XX = set(X)
        cdef frozenset S
        cdef int i
        for i in self._sorted_C_lens:
            if i > len(XX) + 1:
                break
            for C in self._k_C[i]:
                S = C - XX
                if len(S) == 1:
                    XX.add(next(iter(S)))
        return frozenset(XX)

    cpdef _is_isomorphic(self, other, certificate=False):
        """
        Test if ``self`` is isomorphic to ``other``.

        INPUT:

        - ``other`` -- matroid
        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: boolean, and, if ``certificate=True``, a dictionary giving the
        isomorphism or ``None``

        EXAMPLES::

            sage: M = matroids.Spike(3)
            sage: from sage.matroids.basis_matroid import BasisMatroid
            sage: N = BasisMatroid(M)
            sage: M.is_isomorphic(N)
            True
            sage: N = matroids.catalog.Vamos()
            sage: M.is_isomorphic(N)
            False

        .. NOTE::

            Internal version that does no input checking.
        """
        if certificate:
            return self._is_isomorphic(other), self._isomorphism(other)
        N = CircuitsMatroid(other)
        S = SetSystem(self._groundset, self._C)
        O = SetSystem(N._groundset, N._C)
        return S._isomorphism(O) is not None

    # representation

    def _repr_(self):
        """
        Return a string representation of the matroid.

        EXAMPLES::

            sage: matroids.Theta(10)
            Theta_10: Matroid of rank 10 on 20 elements with 490 circuits
            sage: matroids.catalog.NonDesargues()
            NonDesargues: Matroid of rank 3 on 10 elements with 9 nonspanning circuits
        """
        if self._nsc_defined:
            return f'{Matroid._repr_(self)} with {len(self.nonspanning_circuits())} nonspanning circuits'
        else:
            return f'{Matroid._repr_(self)} with {len(self._C)} circuits'

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

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.catalog.Vamos())
            sage: N = CircuitsMatroid(matroids.catalog.Vamos())
            sage: hash(M) == hash(N)
            True
            sage: O = CircuitsMatroid(matroids.catalog.NonVamos())
            sage: hash(M) == hash(O)
            False
        """
        return hash(tuple([self._groundset, frozenset(self._C)]))

    def __richcmp__(left, right, int op):
        r"""
        Compare two matroids.

        We take a very restricted view on equality: the objects need to be of
        the exact same type (so no subclassing) and the internal data need to
        be the same. For CircuitsMatroids, this means that the groundsets and
        the sets of circuits of the two matroids are equal.

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.catalog.Pappus())
            sage: N = CircuitsMatroid(matroids.catalog.NonPappus())
            sage: M == N
            False
            sage: N = Matroid(circuits=M.circuits())
            sage: M == N
            True
        """
        cdef CircuitsMatroid lt, rt
        if op not in [Py_EQ, Py_NE]:
            return NotImplemented
        if type(left) is not type(right):
            return NotImplemented
        lt = <CircuitsMatroid> left
        rt = <CircuitsMatroid> right
        if lt.groundset() != rt.groundset():
            return rich_to_bool(op, 1)
        if lt.full_rank() != rt.full_rank():
            return rich_to_bool(op, 1)
        return richcmp(frozenset(lt._C), frozenset(rt._C), op)

    # copying, loading, saving

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle, (version, data))``, where ``unpickle`` is the
        name of a function that, when called with ``(version, data)``,
        produces a matroid isomorphic to ``self``. ``version`` is an integer
        (currently 0) and ``data`` is a tuple ``(E, C, name)`` where ``E`` is
        the groundset, ``C`` is the list of circuits, and ``name`` is a custom
        name.

        EXAMPLES::

            sage: M = matroids.Theta(5)
            sage: M == loads(dumps(M))  # indirect doctest
            True
            sage: M.reset_name()
            sage: loads(dumps(M))
            Matroid of rank 5 on 10 elements with 45 circuits
        """
        import sage.matroids.unpickling
        data = (self._groundset, frozenset(self._C), self.get_custom_name())
        version = 0
        return sage.matroids.unpickling.unpickle_circuits_matroid, (version, data)

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

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.catalog.RelaxedNonFano())
            sage: sorted(M.groundset())
            [0, 1, 2, 3, 4, 5, 6]
            sage: N = M.relabel({'g': 'x', 0: 'z'})  # 'g': 'x' is ignored
            sage: from sage.matroids.utilities import cmp_elements_key
            sage: sorted(N.groundset(), key=cmp_elements_key)
            [1, 2, 3, 4, 5, 6, 'z']
            sage: M.is_isomorphic(N)
            True

        TESTS::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.catalog.RelaxedNonFano())
            sage: f = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g'}
            sage: N = M.relabel(f)
            sage: for S in powerset(M.groundset()):
            ....:     assert M.rank(S) == N.rank([f[x] for x in S])
        """
        d = self._relabel_map(mapping)
        E = [d[x] for x in self._groundset]
        C = []
        for i in self._k_C:
            C += [[d[y] for y in x] for x in self._k_C[i]]
        M = CircuitsMatroid(groundset=E, circuits=C)
        return M

    # enumeration

    def bases_iterator(self):
        r"""
        Return an iterator over the bases of the matroid.

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.Uniform(2, 4))
            sage: it = M.bases_iterator()
            sage: it.__next__()
            frozenset({0, 1})
            sage: sorted(M.bases_iterator(), key=str)
            [frozenset({0, 1}),
             frozenset({0, 2}),
             frozenset({0, 3}),
             frozenset({1, 2}),
             frozenset({1, 3}),
             frozenset({2, 3})]
        """
        from itertools import combinations
        cdef set NB = set(self.nonbases())
        cdef frozenset S
        for St in combinations(self._groundset, self._matroid_rank):
            S = frozenset(St)
            if S not in NB:
                yield S

    cpdef SetSystem independent_sets(self, long k=-1):
        r"""
        Return the independent sets of the matroid.

        INPUT:

        - ``k`` -- integer (optional); if specified, return the size-`k`
          independent sets of the matroid

        OUTPUT: :class:`SetSystem`

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.catalog.Pappus())
            sage: M.independent_sets(4)
            SetSystem of 0 sets over 9 elements
            sage: M.independent_sets(3)
            SetSystem of 75 sets over 9 elements
            sage: frozenset({'a', 'c', 'e'}) in _
            True

        TESTS::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.CompleteGraphic(4))
            sage: len(M.bases())
            16

        .. SEEALSO::

            :meth:`M.bases() <sage.matroids.circuits_matroid.bases>`
        """
        if k == -1:  # all independent sets
            return self._independent_sets()

        # independent k-sets
        from itertools import combinations
        cdef SetSystem I_k = SetSystem(self._groundset)
        cdef set D_k = set(self.dependent_sets(k))
        cdef frozenset S
        for St in combinations(self._groundset, k):
            S = frozenset(St)
            if S not in D_k:
                I_k.append(S)
        return I_k

    cpdef SetSystem dependent_sets(self, long k):
        r"""
        Return the dependent sets of fixed size.

        INPUT:

        - ``k`` -- integer

        OUTPUT: :class:`SetSystem`

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.catalog.Vamos())
            sage: M.dependent_sets(3)
            SetSystem of 0 sets over 8 elements
            sage: sorted([sorted(X) for X in M.dependent_sets(4)])
            [['a', 'b', 'c', 'd'], ['a', 'b', 'e', 'f'], ['a', 'b', 'g', 'h'],
             ['c', 'd', 'e', 'f'], ['e', 'f', 'g', 'h']]

        TESTS::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.Uniform(2, 4))
            sage: len(M.nonbases())
            0
            sage: M = CircuitsMatroid(matroids.CompleteGraphic(6))
            sage: len(M.nonbases())
            1707
        """
        cdef int i
        cdef set D_k = set()
        cdef frozenset S
        for i in range(min(self._k_C), k + 1):
            if i in self._k_C:
                for S in self._k_C[i]:
                    D_k.add(S)
            if i == k:
                break
            for S in D_k.copy():
                D_k.remove(S)
                for e in S ^ self._groundset:
                    D_k.add(S | set([e]))
        return SetSystem(self._groundset, D_k)

    cpdef SetSystem circuits(self, k=None):
        """
        Return the circuits of the matroid.

        INPUT:

        - ``k`` -- integer (optional); the length of the circuits

        OUTPUT: :class:`SetSystem`

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.Uniform(2, 4))
            sage: M.circuits()
            SetSystem of 4 sets over 4 elements
            sage: list(M.circuits(0))
            []
            sage: sorted(M.circuits(3), key=str)
            [frozenset({0, 1, 2}),
             frozenset({0, 1, 3}),
             frozenset({0, 2, 3}),
             frozenset({1, 2, 3})]
        """
        cdef SetSystem C = SetSystem(self._groundset)
        if k is not None:
            if k in self._k_C:
                for c in self._k_C[k]:
                    C.append(c)
        else:
            for i in self._k_C:
                for c in self._k_C[i]:
                    C.append(c)
        return C

    def circuits_iterator(self, k=None):
        """
        Return an iterator over the circuits of the matroid.

        INPUT:

        - ``k`` -- integer (optional); the length of the circuits

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.Uniform(2, 4))
            sage: sum(1 for C in M.circuits_iterator())
            4
            sage: list(M.circuits_iterator(0))
            []
            sage: sorted(M.circuits_iterator(3), key=str)
            [frozenset({0, 1, 2}),
             frozenset({0, 1, 3}),
             frozenset({0, 2, 3}),
             frozenset({1, 2, 3})]
        """
        if k is not None:
            if k in self._k_C:
                for C in self._k_C[k]:
                    yield C
        else:
            for i in self._k_C:
                for C in self._k_C[i]:
                    yield C

    cpdef SetSystem nonspanning_circuits(self):
        """
        Return the nonspanning circuits of the matroid.

        OUTPUT: :class:`SetSystem`

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.Uniform(2, 4))
            sage: M.nonspanning_circuits()
            SetSystem of 0 sets over 4 elements
            sage: M = matroids.Theta(5)
            sage: M.nonspanning_circuits()
            SetSystem of 15 sets over 10 elements
        """
        cdef SetSystem NSC = SetSystem(self._groundset)
        cdef int i
        cdef frozenset S
        for i in self._sorted_C_lens:
            if i > self._matroid_rank:
                break
            for S in self._k_C[i]:
                NSC.append(S)
        return NSC

    def nonspanning_circuits_iterator(self):
        """
        Return an iterator over the nonspanning circuits of the matroid.

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.Uniform(2, 4))
            sage: list(M.nonspanning_circuits_iterator())
            []
        """
        cdef int i
        for i in self._k_C:
            if i <= self._matroid_rank:
                for C in self._k_C[i]:
                    yield C

    cpdef SetSystem no_broken_circuits_facets(self, ordering=None, reduced=False):
        r"""
        Return the no broken circuits (NBC) facets of ``self``.

        INPUT:

        - ``ordering`` -- list (optional); a total ordering of the groundset
        - ``reduced`` -- boolean (default: ``False``)

        OUTPUT: :class:`SetSystem`

        EXAMPLES::

            sage: M = Matroid(circuits=[[0, 1, 2]])
            sage: M.no_broken_circuits_facets(ordering=[1, 0, 2])
            SetSystem of 2 sets over 3 elements
            sage: sorted([sorted(X) for X in _])
            [[0, 1], [1, 2]]
            sage: M.no_broken_circuits_facets(ordering=[1, 0, 2], reduced=True)
            SetSystem of 2 sets over 3 elements
            sage: sorted([sorted(X) for X in _])
            [[0], [2]]
        """
        from itertools import combinations
        from sage.matroids.utilities import cmp_elements_key
        if ordering is None:
            ordering = sorted(self._groundset, key=cmp_elements_key)
        else:
            if frozenset(ordering) != self._groundset:
                raise ValueError("not an ordering of the groundset")

        cdef int i, r = self._matroid_rank
        cdef frozenset min_e = frozenset([ordering[0]])
        cdef frozenset S
        # compute broken circuits (with minimal element added)
        cdef dict BC = {}
        for i in range(r + 1):
            BC[i] = set()
        for C in self._C:
            for e in ordering:
                if e in C:
                    S = C - frozenset([e]) | min_e
                    if len(S) <= r:
                        BC[len(S)].add(S)
                    break

        for i in range(r):
            for S in BC[i].copy():
                BC[i].remove(S)
                for e in self._groundset ^ S:
                    BC[i+1].add(S | set([e]))

        cdef SetSystem B = SetSystem(self._groundset)
        for St in combinations(ordering[1:], self._matroid_rank - 1):
            S = frozenset(St)
            if S | min_e not in BC[r]:
                if not reduced:
                    B.append(S | min_e)
                else:
                    B.append(S)

        return B

    cpdef SetSystem no_broken_circuits_sets(self, ordering=None, reduced=False):
        r"""
        Return the no broken circuits (NBC) sets of ``self``.

        An NBC set is a subset `A` of the groundset under some total
        ordering `<` such that `A` contains no broken circuit.

        INPUT:

        - ``ordering`` -- list (optional); a total ordering of the groundset
        - ``reduced`` -- boolean (default: ``False``)

        OUTPUT: :class:`SetSystem`

        EXAMPLES::

            sage: M = Matroid(circuits=[[1, 2, 3], [3, 4, 5], [1, 2, 4, 5]])
            sage: SimplicialComplex(M.no_broken_circuits_sets())
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 2, 4), (1, 2, 5), (1, 3, 4), (1, 3, 5)}
            sage: SimplicialComplex(M.no_broken_circuits_sets([5, 4, 3, 2, 1]))
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 3, 5), (1, 4, 5), (2, 3, 5), (2, 4, 5)}

        ::

            sage: M = Matroid(circuits=[[1, 2, 3], [1, 4, 5], [2, 3, 4, 5]])
            sage: SimplicialComplex(M.no_broken_circuits_sets([5, 4, 3, 2, 1]))
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 3, 5), (2, 3, 5), (2, 4, 5), (3, 4, 5)}

        TESTS::

            sage: M = Matroid(circuits=[[1, 2, 3], [3, 4, 5], [1, 2, 4, 5]])
            sage: C1 = SimplicialComplex(M.no_broken_circuits_sets())
            sage: from sage.matroids.basis_matroid import BasisMatroid
            sage: M = BasisMatroid(Matroid(circuits=[[1, 2, 3], [3, 4, 5], [1, 2, 4, 5]]))
            sage: C2 = SimplicialComplex(M.no_broken_circuits_sets())
            sage: C1 == C2
            True
        """
        from sage.topology.simplicial_complex import SimplicialComplex
        cdef SetSystem NBC = SetSystem(self._groundset)
        for f in SimplicialComplex(self.no_broken_circuits_facets(ordering, reduced),
                                   maximality_check=False).face_iterator():
            NBC.append(frozenset(f))
        return NBC

    cpdef broken_circuit_complex(self, ordering=None, reduced=False):
        r"""
        Return the broken circuit complex of ``self``.

        The broken circuit complex of a matroid with a total ordering `<`
        on the groundset is obtained from the
        :meth:`NBC sets <no_broken_circuits_sets>` under subset inclusion.

        INPUT:

        - ``ordering`` -- list (optional); a total ordering of the groundset
        - ``reduced`` -- boolean (default: ``False``); whether to return the
          reduced broken circuit complex (the link at the smallest element)

        OUTPUT: a simplicial complex of the NBC sets under inclusion

        EXAMPLES::

            sage: M = Matroid(circuits=[[1, 2, 3], [3, 4, 5], [1, 2, 4, 5]])
            sage: M.broken_circuit_complex()
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 2, 4), (1, 2, 5), (1, 3, 4), (1, 3, 5)}
            sage: M.broken_circuit_complex([5, 4, 3, 2, 1])
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 3, 5), (1, 4, 5), (2, 3, 5), (2, 4, 5)}
            sage: M.broken_circuit_complex([5, 4, 3, 2, 1], reduced=True)
            Simplicial complex with vertex set (1, 2, 3, 4)
             and facets {(1, 3), (1, 4), (2, 3), (2, 4)}

        For a matroid with loops, the broken circuit complex is not defined,
        and the method yields an error::

            sage: M = Matroid(groundset=[0, 1, 2], circuits=[[0]])
            sage: M.broken_circuit_complex()
            Traceback (most recent call last):
            ...
            ValueError: broken circuit complex of matroid with loops is not defined
        """
        from sage.topology.simplicial_complex import SimplicialComplex
        if self.loops():
            raise ValueError("broken circuit complex of matroid with loops is not defined")
        return SimplicialComplex(self.no_broken_circuits_facets(ordering, reduced), maximality_check=False)

    # properties

    cpdef girth(self):
        r"""
        Return the girth of the matroid.

        The girth is the size of the smallest circuit. In case the matroid has
        no circuits the girth is `\infty`.

        EXAMPLES::

            sage: matroids.Theta(10).girth()
            3

        REFERENCES:

        [Oxl2011]_, p. 327.
        """
        from sage.rings.infinity import infinity
        return min(self._k_C, default=infinity)

    cpdef bint is_paving(self) noexcept:
        """
        Return if ``self`` is paving.

        A matroid is paving if each of its circuits has size `r` or `r+1`.

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.catalog.Vamos())
            sage: M.is_paving()
            True
        """
        return self.girth() >= self._matroid_rank

    # verification

    cpdef bint is_valid(self) noexcept:
        r"""
        Test if ``self`` obeys the matroid axioms.

        For a matroid defined by its circuits, we check the circuit axioms.

        OUTPUT: boolean

        EXAMPLES::

            sage: C = [[1, 2, 3], [3, 4, 5], [1, 2, 4, 5]]
            sage: M = Matroid(circuits=C)
            sage: M.is_valid()
            True
            sage: C = [[1, 2], [1, 2, 3], [3, 4, 5], [1, 2, 4, 5]]
            sage: M = Matroid(circuits=C)
            sage: M.is_valid()
            False
            sage: C = [[3, 6], [1, 2, 3], [3, 4, 5], [1, 2, 4, 5]]
            sage: M = Matroid(circuits=C)
            sage: M.is_valid()
            False
            sage: C = [[3, 6], [1, 2, 3], [3, 4, 5], [1, 2, 6], [6, 4, 5], [1, 2, 4, 5]]
            sage: M = Matroid(circuits=C)
            sage: M.is_valid()
            True
            sage: C = [[], [1, 2, 3], [3, 4, 5], [1, 2, 4, 5]]
            sage: M = Matroid(circuits=C)
            sage: M.is_valid()
            False
            sage: C = [[1, 2, 3], [3, 4, 5]]
            sage: M = Matroid(circuits=C)
            sage: M.is_valid()
            False
        """
        from itertools import combinations_with_replacement
        cdef int i, j
        cdef frozenset C1, C2, I12, U12
        for (i, j) in combinations_with_replacement(self._sorted_C_lens, 2):
            # loop through all circuit length pairs (i, j) with i <= j
            for C1 in self._k_C[i]:
                if not C1:  # the empty set can't be a circuit
                    return False
                for C2 in self._k_C[j]:
                    I12 = C1 & C2
                    if not I12:  # C1 and C2 are disjoint; nothing to test
                        continue
                    if len(I12) == len(C1):
                        if len(C1) == len(C2):  # they are the same circuit
                            break
                        # C1 < C2; a circuit can't be a subset of another circuit
                        return False
                    # check circuit elimination axiom
                    U12 = C1 | C2
                    for e in I12:
                        if self._is_independent(U12 - {e}):
                            return False
        return True
