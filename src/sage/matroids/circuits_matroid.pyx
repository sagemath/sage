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

from sage.structure.richcmp cimport rich_to_bool, richcmp
from sage.matroids.matroid cimport Matroid
from sage.matroids.set_system cimport SetSystem
from cpython.object cimport Py_EQ, Py_NE


cdef class CircuitsMatroid(Matroid):
    r"""
    A matroid defined by its circuits.

    INPUT:

    - ``M`` -- a matroid (default: ``None``)
    - ``groundset`` -- a list (default: ``None``); the groundset of the matroid
    - ``circuits`` -- a list (default: ``None``); the collection of circuits of
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
            self._groundset = frozenset(M.groundset())
            self._C = SetSystem(list(M.groundset()), frozenset([frozenset(C) for C in M.circuits()]))
        else:
            self._groundset = frozenset(groundset)
            self._C = SetSystem(list(groundset), frozenset([frozenset(C) for C in circuits]))
        # k-circuits
        self._k_C = {}
        for C in self._C:
            try:
                self._k_C[len(C)] += [C]
            except KeyError:
                self._k_C[len(C)] = []
                self._k_C[len(C)] += [C]
        self._matroid_rank = self.rank(self._groundset)
        self._nsc_defined = nsc_defined

    cpdef groundset(self) noexcept:
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT: a set

        EXAMPLES::

            sage: M = matroids.Theta(2)
            sage: sorted(M.groundset())
            ['x0', 'x1', 'y0', 'y1']
        """
        return self._groundset

    cpdef _rank(self, X) noexcept:
        """
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and ``X`` may be assumed to have
        the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface

        OUTPUT: an integer; the rank of ``X`` in the matroid

        EXAMPLES::

            sage: M = matroids.Theta(3)
            sage: M._rank(['x1', 'y0', 'y2'])
            2
        """
        return len(self._max_independent(X))

    # optional

    cpdef full_rank(self) noexcept:
        r"""
        Return the rank of the matroid.

        The *rank* of the matroid is the size of the largest independent
        subset of the groundset.

        OUTPUT: an integer; the rank of the matroid

        EXAMPLES::

            sage: M = matroids.Theta(20)
            sage: M.full_rank()
            20
        """
        return self._matroid_rank

    cpdef _is_independent(self, F) noexcept:
        """
        Test if input is independent.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: boolean

        EXAMPLES::

            sage: M = matroids.Theta(4)
            sage: M._is_independent(['y0', 'y1', 'y3', 'x2'])
            False
            sage: M._is_independent(['y0', 'y2', 'y3', 'x2'])
            True
        """
        cdef set I = set(F)
        cdef int s = len(F)
        for i in self._k_C:
            if i <= s:
                for C in self._k_C[i]:
                    if C <= I:
                        return False
        return True

    cpdef _max_independent(self, F) noexcept:
        """
        Compute a maximal independent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: a frozenset; a maximal independent subset of ``X``

        EXAMPLES::

            sage: M = matroids.Theta(6)
            sage: len(M._max_independent(M.groundset()))
            6
        """
        cdef set I = set(F)
        for i in self._k_C:
            for C in self._k_C[i]:
                if i <= len(I) and i > 0:
                    if C <= I:
                        e = next(iter(C))
                        I.remove(e)

        return frozenset(I)

    cpdef _circuit(self, F) noexcept:
        """
        Return a minimal dependent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT: a frozenset; a circuit contained in ``X``, if it exists.
        Otherwise an error is raised.

        EXAMPLES::

            sage: M = matroids.Theta(4)
            sage: sorted(M._circuit(['y0', 'y1', 'y3', 'x2']))
            ['x2', 'y0', 'y1', 'y3']
            sage: M._circuit(['y0', 'y2', 'y3', 'x2'])
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set
        """
        cdef set I = set(F)
        for C in self._C:
            if C <= I:
                return C
        raise ValueError("no circuit in independent set")

    cpdef _is_isomorphic(self, other, certificate=False) noexcept:
        """
        Test if ``self`` is isomorphic to ``other``.

        INPUT:

        - ``other`` -- a matroid
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
        return self._C._isomorphism(N._C) is not None

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
            return Matroid._repr_(self) + " with " + str(len(self.nonspanning_circuits())) + " nonspanning circuits"
        else:
            return Matroid._repr_(self) + " with " + str(len(self._C)) + " circuits"

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
        return hash(tuple([self.groundset(), frozenset(self._C)]))

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

    def __copy__(self):
        """
        Create a shallow copy.

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.catalog.Vamos())
            sage: N = copy(M)  # indirect doctest
            sage: M == N
            True
            sage: M.groundset() is N.groundset()
            True
        """
        N = CircuitsMatroid(groundset=[], circuits=[])
        N._groundset = self._groundset
        N._C = self._C
        N._k_C = self._k_C
        N._nsc_defined = self._nsc_defined
        N._matroid_rank = self._matroid_rank
        N.rename(self.get_custom_name())
        return N

    def __deepcopy__(self, memo=None):
        """
        Create a deep copy.

        .. NOTE::

            Since matroids are immutable, a shallow copy normally suffices.

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.catalog.Vamos())
            sage: N = deepcopy(M)  # indirect doctest
            sage: M == N
            True
            sage: M.groundset() is N.groundset()
            False
        """
        if memo is None:
            memo = {}
        from copy import deepcopy
        # Since matroids are immutable, N cannot reference itself in correct code, so no need to worry about the recursion.
        N = CircuitsMatroid(groundset=deepcopy(self._groundset, memo), circuits=deepcopy(frozenset(self._C), memo))
        N.rename(deepcopy(self.get_custom_name(), memo))
        return N

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

    # enumeration

    cpdef bases(self) noexcept:
        r"""
        Return the bases of the matroid.

        OUTPUT: a SetSystem

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.Uniform(2, 4))
            sage: len(M.bases())
            6
        """
        cdef SetSystem B, NSC
        cdef bint flag
        B = SetSystem(list(self.groundset()))
        NSC = self.nonspanning_circuits()
        from itertools import combinations
        for S in combinations(self._groundset, self._matroid_rank):
            flag = True
            S = frozenset(S)
            for C in NSC:
                if C <= S:
                    flag = False
                    break
            if flag:
                B.append(S)
        return B

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
        cdef SetSystem NSC = self.nonspanning_circuits()
        for B in combinations(self._groundset, self._matroid_rank):
            B = frozenset(B)
            if not any(C <= B for C in NSC):
                yield B

    cpdef circuits(self, k=None) noexcept:
        """
        Return the circuits of the matroid.

        INPUT:

        - ``k`` -- an integer (optional); the length of the circuits

        OUTPUT: a SetSystem

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.Uniform(2, 4))
            sage: M.circuits()
            Iterator over a system of subsets
            sage: list(M.circuits(0))
            []
            sage: sorted(M.circuits(3), key=str)
            [frozenset({0, 1, 2}),
             frozenset({0, 1, 3}),
             frozenset({0, 2, 3}),
             frozenset({1, 2, 3})]
        """
        cdef SetSystem C
        C = SetSystem(list(self.groundset()))
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

        - ``k`` -- an integer (optional); the length of the circuits

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

    cpdef nonspanning_circuits(self) noexcept:
        """
        Return the nonspanning circuits of the matroid.

        OUTPUT: a SetSystem

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.Uniform(2, 4))
            sage: M.nonspanning_circuits()
            Iterator over a system of subsets
        """
        cdef list NSC = []
        for i in self._k_C:
            if i <= self.rank():
                NSC.extend(self._k_C[i])
        return SetSystem(list(self.groundset()), NSC)

    def nonspanning_circuits_iterator(self):
        """
        Return an iterator over the nonspanning circuits of the matroid.

        EXAMPLES::

            sage: from sage.matroids.circuits_matroid import CircuitsMatroid
            sage: M = CircuitsMatroid(matroids.Uniform(2, 4))
            sage: list(M.nonspanning_circuits_iterator())
            []
        """
        for i in self._k_C:
            if i <= self.rank():
                for C in self._k_C[i]:
                    yield C

    cpdef no_broken_circuits_sets(self, ordering=None) noexcept:
        r"""
        Return the no broken circuits (NBC) sets of ``self``.

        An NBC set is a subset `A` of the ground set under some total
        ordering `<` such that `A` contains no broken circuit.

        INPUT:

        - ``ordering`` -- a total ordering of the groundset given as a list

        OUTPUT: a list of frozensets

        EXAMPLES::

            sage: M = Matroid(circuits=[[1,2,3], [3,4,5], [1,2,4,5]])
            sage: SimplicialComplex(M.no_broken_circuits_sets())
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 2, 4), (1, 2, 5), (1, 3, 4), (1, 3, 5)}
            sage: SimplicialComplex(M.no_broken_circuits_sets([5,4,3,2,1]))
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 3, 5), (1, 4, 5), (2, 3, 5), (2, 4, 5)}

        ::

            sage: M = Matroid(circuits=[[1,2,3], [1,4,5], [2,3,4,5]])
            sage: SimplicialComplex(M.no_broken_circuits_sets([5,4,3,2,1]))
            Simplicial complex with vertex set (1, 2, 3, 4, 5)
             and facets {(1, 3, 5), (2, 3, 5), (2, 4, 5), (3, 4, 5)}

        TESTS::

            sage: M = Matroid(circuits=[[1,2,3], [3,4,5], [1,2,4,5]])
            sage: C1 = SimplicialComplex(M.no_broken_circuits_sets())
            sage: from sage.matroids.basis_matroid import BasisMatroid
            sage: M = BasisMatroid(Matroid(circuits=[[1,2,3], [3,4,5], [1,2,4,5]]))
            sage: C2 = SimplicialComplex(M.no_broken_circuits_sets())
            sage: C1 == C2
            True
        """
        if ordering is None:
            ordering = sorted(self.groundset(), key=str)
        else:
            if frozenset(ordering) != self.groundset():
                raise ValueError("not an ordering of the groundset")

        # compute broken circuits
        cdef list BC = []
        for C in self._C:
            for e in ordering:
                if e in C:
                    BC.append(C - set([e]))
                    break

        cdef list F = []  # broken circuit complex facets
        for B in self.bases():
            flag = True
            for bc in BC:
                if bc <= B:
                    flag = False
                    break
            if flag:
                F.append(B)

        from sage.topology.simplicial_complex import SimplicialComplex
        return [frozenset(f) for f in SimplicialComplex(F).face_iterator()]

    # properties

    cpdef girth(self) noexcept:
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
        return min(self._k_C, default=float('inf'))

    cpdef is_paving(self) noexcept:
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
        return self.girth() >= self.rank()

    # verification

    cpdef is_valid(self) noexcept:
        r"""
        Test if ``self`` obeys the matroid axioms.

        For a matroid defined by its circuits, we check the circuit axioms.

        OUTPUT: boolean

        EXAMPLES::

            sage: C = [[1, 2, 3], [3, 4, 5], [1, 2, 4, 5]]
            sage: M = Matroid(circuits=C)
            sage: M.is_valid()
            True
            sage: C = [[1,2], [1, 2, 3], [3, 4, 5], [1, 2, 4, 5]]
            sage: M = Matroid(circuits=C)
            sage: M.is_valid()
            False
            sage: C = [[3,6], [1, 2, 3], [3, 4, 5], [1, 2, 4, 5]]
            sage: M = Matroid(circuits=C)
            sage: M.is_valid()
            False
            sage: C = [[3,6], [1, 2, 3], [3, 4, 5], [1, 2, 6], [6, 4, 5], [1, 2, 4, 5]]
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
        cdef int i, j, k, S_len
        cdef frozenset C1, C2, C3, I12, U12
        cdef bint flag
        for (i, j) in combinations_with_replacement(sorted(self._k_C), 2):
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
                    S_len = len(U12) - 1  # the size of S below
                    for e in I12:
                        flag = False
                        S = U12 - {e}
                        for k in self._k_C:
                            if k <= S_len:
                                for C3 in self._k_C[k]:
                                    if C3 <= S:
                                        flag = True
                                        break
                            if flag:
                                break
                        if not flag:
                            return False
        return True
