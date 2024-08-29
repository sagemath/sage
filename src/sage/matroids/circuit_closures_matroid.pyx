r"""
Circuit closures matroids

Matroids are characterized by a list of all tuples `(C, k)`, where `C` is the
closure of a circuit, and `k` the rank of `C`. The CircuitClosuresMatroid
class implements matroids using this information as data.

Construction
============

A ``CircuitClosuresMatroid`` can be created from another matroid or from a
list of circuit-closures. For a full description of allowed inputs, see
:class:`below <sage.matroids.circuit_closures_matroid.CircuitClosuresMatroid>`.
It is recommended to use the
:func:`Matroid() <sage.matroids.constructor.Matroid>` function for a more
flexible construction of a ``CircuitClosuresMatroid``. For direct access to
the ``CircuitClosuresMatroid`` constructor, run::

    sage: from sage.matroids.advanced import *

See also :mod:`sage.matroids.advanced`.

EXAMPLES::

    sage: from sage.matroids.advanced import *
    sage: M1 = CircuitClosuresMatroid(groundset='abcdef',
    ....:                 circuit_closures={2: ['abc', 'ade'], 3: ['abcdef']})
    sage: M2 = Matroid(circuit_closures={2: ['abc', 'ade'], 3: ['abcdef']})
    sage: M3 = Matroid(circuit_closures=[(2, 'abc'),
    ....:                                (3, 'abcdef'), (2, 'ade')])
    sage: M1 == M2
    True
    sage: M1 == M3
    True

Note that the class does not implement custom minor and dual operations::

    sage: from sage.matroids.advanced import *
    sage: M = CircuitClosuresMatroid(groundset='abcdef',
    ....:                 circuit_closures={2: ['abc', 'ade'], 3: ['abcdef']})
    sage: isinstance(M.contract('a'), MinorMatroid)
    True
    sage: isinstance(M.dual(), DualMatroid)
    True

AUTHORS:

- Rudi Pendavingh, Stefan van Zwam (2013-04-01): initial version
"""

# ****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport Py_EQ, Py_NE
from sage.structure.richcmp cimport rich_to_bool, richcmp
from sage.matroids.matroid cimport Matroid
from sage.matroids.set_system cimport SetSystem
from sage.matroids.utilities import setprint_s, cmp_elements_key

cdef class CircuitClosuresMatroid(Matroid):
    r"""
    A general matroid `M` is characterized by its rank `r(M)` and the set of
    pairs

    `(k, \{` closure `(C) : C ` circuit of ` M, r(C)=k\})` for `k=0, .., r(M)-1`

    As each independent set of size `k` is in at most one closure(`C`) of rank
    `k`, and each closure(`C`) of rank `k` contains at least `k + 1`
    independent sets of size `k`, there are at most `\binom{n}{k}/(k + 1)`
    such closures-of-circuits of rank `k`. Each closure(`C`) takes `O(n)` bits
    to store, giving an upper bound of `O(2^n)` on the space complexity of the
    entire matroid.

    A subset `X` of the groundset is independent if and only if

    `| X \cap ` closure `(C) | \leq k` for all circuits `C` of `M` with
    `r(C)=k`.

    So determining whether a set is independent takes time proportional to the
    space complexity of the matroid.

    INPUT:

    - ``M`` -- matroid (default: ``None``)
    - ``groundset`` -- groundset of a matroid (default: ``None``)
    - ``circuit_closures`` -- dictionary (default: ``None``); the collection of
      circuit closures of a matroid presented as a dictionary whose keys are
      ranks, and whose values are sets of circuit closures of the specified rank

    OUTPUT:

    - If the input is a matroid ``M``, return a ``CircuitClosuresMatroid``
      instance representing ``M``.
    - Otherwise, return a ``CircuitClosuresMatroid`` instance based on
      ``groundset`` and ``circuit_closures``.

    .. NOTE::

        For a more flexible means of input, use the ``Matroid()`` function.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = CircuitClosuresMatroid(matroids.catalog.Fano())
        sage: M
        Matroid of rank 3 on 7 elements with circuit-closures
        {2: {{'a', 'b', 'f'}, {'a', 'c', 'e'}, {'a', 'd', 'g'},
             {'b', 'c', 'd'}, {'b', 'e', 'g'}, {'c', 'f', 'g'},
             {'d', 'e', 'f'}}, 3: {{'a', 'b', 'c', 'd', 'e', 'f', 'g'}}}
        sage: M = CircuitClosuresMatroid(groundset='abcdefgh',
        ....:            circuit_closures={3: ['edfg', 'acdg', 'bcfg', 'cefh',
        ....:                 'afgh', 'abce', 'abdf', 'begh', 'bcdh', 'adeh'],
        ....:                              4: ['abcdefgh']})
        sage: M.equals(matroids.catalog.P8())
        True
    """

    # necessary (__init__, groundset, _rank)

    def __init__(self, M=None, groundset=None, circuit_closures=None):
        """
        Initialization of the matroid. See class docstring for full
        documentation.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = CircuitClosuresMatroid(matroids.catalog.Fano())
            sage: M
            Matroid of rank 3 on 7 elements with circuit-closures
            {2: {{'a', 'b', 'f'}, {'a', 'c', 'e'}, {'a', 'd', 'g'},
                 {'b', 'c', 'd'}, {'b', 'e', 'g'}, {'c', 'f', 'g'},
                 {'d', 'e', 'f'}},
             3: {{'a', 'b', 'c', 'd', 'e', 'f', 'g'}}}

            sage: M = CircuitClosuresMatroid(groundset='abcdefgh',
            ....:        circuit_closures={3: ['edfg', 'acdg', 'bcfg', 'cefh',
            ....:             'afgh', 'abce', 'abdf', 'begh', 'bcdh', 'adeh'],
            ....:                          4: ['abcdefgh']})
            sage: M.equals(matroids.catalog.P8())
            True

        TESTS::

            sage: from sage.matroids.advanced import *
            sage: M = CircuitClosuresMatroid(matroids.catalog.Fano())
            sage: TestSuite(M).run()
        """
        if M is not None:
            self._groundset = M.groundset()
            self._circuit_closures = M.circuit_closures()
        else:
            self._groundset = frozenset(groundset)
            self._circuit_closures = {}
            for k in circuit_closures:
                self._circuit_closures[k] = frozenset([frozenset(X) for X in circuit_closures[k]])
        self._matroid_rank = self.rank(self._groundset)

    cpdef frozenset groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT: frozenset

        EXAMPLES::

            sage: M = matroids.catalog.Pappus()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
        """
        return frozenset(self._groundset)

    cpdef int _rank(self, frozenset X):
        """
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and
        ``X`` may be assumed to have the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface

        OUTPUT: the rank of ``X`` in the matroid

        EXAMPLES::

            sage: M = matroids.catalog.NonPappus()
            sage: M._rank(frozenset('abc'))
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

            sage: M = matroids.catalog.Vamos()
            sage: M.full_rank()
            4
            sage: M.dual().full_rank()
            4
        """
        return self._matroid_rank

    cpdef bint _is_independent(self, frozenset F):
        """
        Test if input is independent.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: boolean

        EXAMPLES::

            sage: M = matroids.catalog.Vamos()
            sage: M._is_independent(frozenset(['a', 'b', 'c']))
            True
            sage: M._is_independent(frozenset(['a', 'b', 'c', 'd']))
            False
        """
        for r in sorted(self._circuit_closures):
            if len(F) <= r:
                break
            for C in self._circuit_closures[r]:
                S = F & C
                if len(S) > r:
                    return False
        return True

    cpdef frozenset _max_independent(self, frozenset F):
        """
        Compute a maximal independent subset.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: a maximal independent subset of ``X``

        EXAMPLES::

            sage: M = matroids.catalog.Vamos()
            sage: X = M._max_independent(frozenset(['a', 'c', 'd', 'e', 'f']))
            sage: sorted(X)  # random
            ['a', 'd', 'e', 'f']
            sage: M.is_independent(X)
            True
            sage: all(M.is_dependent(X.union([y])) for y in M.groundset() if y not in X)
            True
        """
        I = set(F)
        for r in sorted(self._circuit_closures.keys()):
            if len(I) == 0:
                break
            for C in self._circuit_closures[r]:
                if len(I) == 0:
                    break
                S = I & C
                while(len(S) > r):
                    I.discard(S.pop())

        return frozenset(I)

    cpdef frozenset _circuit(self, frozenset F):
        """
        Return a minimal dependent subset.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: a circuit contained in ``X``, if it exists; otherwise, an error
        is raised

        EXAMPLES::

            sage: M = matroids.catalog.Vamos()
            sage: sorted(M._circuit(frozenset(['a', 'c', 'd', 'e', 'f'])))
            ['c', 'd', 'e', 'f']
            sage: sorted(M._circuit(frozenset(['a', 'c', 'd'])))
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set
        """
        for r in sorted(self._circuit_closures):
            for C in self._circuit_closures[r]:
                S = set(F & C)
                if len(S) > r:
                    while len(S) > r + 1:
                        S.pop()
                    return frozenset(S)
        raise ValueError("no circuit in independent set")

    cpdef dict circuit_closures(self):
        """
        Return the closures of circuits of the matroid.

        A *circuit closure* is a closed set containing a circuit.

        OUTPUT: dictionary containing the circuit closures of the matroid,
        indexed by their ranks

        .. SEEALSO::

            :meth:`Matroid.circuit() <sage.matroids.matroid.Matroid.circuit>`,
            :meth:`Matroid.closure() <sage.matroids.matroid.Matroid.closure>`

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = CircuitClosuresMatroid(matroids.catalog.Fano())
            sage: CC = M.circuit_closures()
            sage: len(CC[2])
            7
            sage: len(CC[3])
            1
            sage: len(CC[1])
            Traceback (most recent call last):
            ...
            KeyError: 1
            sage: [sorted(X) for X in CC[3]]
            [['a', 'b', 'c', 'd', 'e', 'f', 'g']]
        """
        return self._circuit_closures

    cpdef _is_isomorphic(self, other, certificate=False):
        """
        Test if ``self`` is isomorphic to ``other``.

        Internal version that performs no checks on input.

        INPUT:

        - ``other`` -- matroid
        - ``certificate`` -- boolean (default: ``False``)

        OUTPUT: boolean, and, if ``certificate = True``, a dictionary giving
        the isomorphism or ``None``

        .. NOTE::

            Internal version that does no input checking.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M1 = CircuitClosuresMatroid(matroids.Wheel(3))
            sage: M2 = matroids.CompleteGraphic(4)                                      # needs sage.graphs
            sage: M1._is_isomorphic(M2)                                                 # needs sage.graphs
            True
            sage: M1._is_isomorphic(M2, certificate=True)                               # needs sage.graphs
            (True, {0: 0, 1: 1, 2: 2, 3: 3, 4: 5, 5: 4})
            sage: M1 = CircuitClosuresMatroid(matroids.catalog.Fano())
            sage: M2 = matroids.catalog.NonFano()
            sage: M1._is_isomorphic(M2)
            False
            sage: M1._is_isomorphic(M2, certificate=True)
            (False, None)
        """
        if certificate:
            return self._is_isomorphic(other), self._isomorphism(other)
        N = CircuitClosuresMatroid(other)
        if sorted(self._circuit_closures.keys()) != sorted(N._circuit_closures.keys()):
            return False
        SM = SetSystem(self.groundset())
        for r in self._circuit_closures:
            for C in self._circuit_closures[r]:
                SM.append(C)
        SN = SetSystem(N.groundset())
        for r in N._circuit_closures:
            for C in N._circuit_closures[r]:
                SN.append(C)
        return SM._isomorphism(SN) is not None

    # representation

    def _repr_(self):
        """
        Return a string representation of the matroid.

        EXAMPLES::

            sage: M = matroids.catalog.Vamos()
            sage: print(M._repr_())
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
                 {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'},
                 {'e', 'f', 'g', 'h'}},
             4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        """
        return Matroid._repr_(self) + " with circuit-closures\n" + setprint_s(self._circuit_closures)

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

            sage: M = matroids.catalog.Vamos()
            sage: N = matroids.catalog.Vamos()
            sage: hash(M) == hash(N)
            True
            sage: O = matroids.catalog.NonVamos()
            sage: hash(M) == hash(O)
            False
        """
        return hash(tuple([self.groundset(), tuple([(r, len(self._circuit_closures[r])) for r in sorted(self._circuit_closures.keys())])]))

    def __richcmp__(left, right, int op):
        r"""
        Compare two matroids.

        We take a very restricted view on equality: the objects need to be of
        the exact same type (so no subclassing) and the internal data need to
        be the same. For ``BasisMatroid``s, this means that the groundsets and
        the sets of bases of the two matroids are equal.

        EXAMPLES::

            sage: M = matroids.catalog.Pappus()
            sage: N = matroids.catalog.NonPappus()
            sage: M == N
            False
            sage: N = Matroid(M.bases())
            sage: M == N
            False
        """
        cdef CircuitClosuresMatroid lt, rt
        if op not in [Py_EQ, Py_NE]:
            return NotImplemented
        if type(left) is not type(right):
            return NotImplemented
        lt = <CircuitClosuresMatroid> left
        rt = <CircuitClosuresMatroid> right
        if lt.groundset() != rt.groundset():
            return rich_to_bool(op, 1)
        if lt.full_rank() != rt.full_rank():
            return rich_to_bool(op, 1)
        return richcmp(lt._circuit_closures, rt._circuit_closures, op)

    # copying, loading, saving

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle, (version, data))``, where ``unpickle`` is the
        name of a function that, when called with ``(version, data)``,
        produces a matroid isomorphic to ``self``. ``version`` is an integer
        (currently 0) and ``data`` is a tuple ``(E, CC, name)`` where ``E`` is
        the groundset, ``CC`` is the dictionary of circuit closures, and
        ``name`` is a custom name.

        EXAMPLES::

            sage: M = matroids.catalog.Vamos()
            sage: M == loads(dumps(M))  # indirect doctest
            True
            sage: M.reset_name()
            sage: loads(dumps(M))
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
                 {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'},
                 {'e', 'f', 'g', 'h'}},
             4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        """
        import sage.matroids.unpickling
        data = (self._groundset, self._circuit_closures, self.get_custom_name())
        version = 0
        return sage.matroids.unpickling.unpickle_circuit_closures_matroid, (version, data)

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

            sage: from sage.matroids.circuit_closures_matroid import CircuitClosuresMatroid
            sage: M = CircuitClosuresMatroid(matroids.catalog.RelaxedNonFano())
            sage: sorted(M.groundset())
            [0, 1, 2, 3, 4, 5, 6]
            sage: N = M.relabel({'g': 'x', 0: 'z'})  # 'g': 'x' is ignored
            sage: from sage.matroids.utilities import cmp_elements_key
            sage: sorted(N.groundset(), key=cmp_elements_key)
            [1, 2, 3, 4, 5, 6, 'z']
            sage: M.is_isomorphic(N)
            True

        TESTS::

            sage: from sage.matroids.circuit_closures_matroid import CircuitClosuresMatroid
            sage: M = CircuitClosuresMatroid(matroids.catalog.RelaxedNonFano())
            sage: f = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e', 5: 'f', 6: 'g'}
            sage: N = M.relabel(f)
            sage: for S in powerset(M.groundset()):
            ....:     assert M.rank(S) == N.rank([f[x] for x in S])
        """
        d = self._relabel_map(mapping)
        E = [d[x] for x in self.groundset()]
        CC = {}
        for i in self.circuit_closures():
            CC[i] = [[d[y] for y in x] for x in self._circuit_closures[i]]
        M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
        return M

# todo: customized minor, extend methods.
