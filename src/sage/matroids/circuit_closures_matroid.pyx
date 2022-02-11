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

TESTS::

    sage: from sage.matroids.advanced import *
    sage: M = CircuitClosuresMatroid(matroids.named_matroids.Fano())
    sage: TestSuite(M).run()

Methods
=======
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

from sage.structure.richcmp cimport rich_to_bool, richcmp
from .matroid cimport Matroid
from .set_system cimport SetSystem
from .utilities import setprint_s
from cpython.object cimport Py_EQ, Py_NE


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

    A subset `X` of the ground set is independent if and only if

    `| X \cap ` closure `(C) | \leq k` for all circuits `C` of `M` with
    `r(C)=k`.

    So determining whether a set is independent takes time proportional to the
    space complexity of the matroid.

    INPUT:

    - ``M`` -- (default: ``None``) an arbitrary matroid.
    - ``groundset`` -- (default: ``None``) the groundset of a matroid.
    - ``circuit_closures`` -- (default: ``None``) the collection of circuit
      closures of a matroid, presented as a dictionary whose keys are ranks,
      and whose values are sets of circuit closures of the specified rank.

    OUTPUT:

    - If the input is a matroid ``M``, return a ``CircuitClosuresMatroid``
      instance representing ``M``.
    - Otherwise, return a ``CircuitClosuresMatroid`` instance based on
      ``groundset`` and ``circuit_closures``.

    .. NOTE::

        For a more flexible means of input, use the ``Matroid()`` function.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = CircuitClosuresMatroid(matroids.named_matroids.Fano())
        sage: M
        Matroid of rank 3 on 7 elements with circuit-closures
        {2: {{'a', 'b', 'f'}, {'a', 'c', 'e'}, {'a', 'd', 'g'},
             {'b', 'c', 'd'}, {'b', 'e', 'g'}, {'c', 'f', 'g'},
             {'d', 'e', 'f'}}, 3: {{'a', 'b', 'c', 'd', 'e', 'f', 'g'}}}
        sage: M = CircuitClosuresMatroid(groundset='abcdefgh',
        ....:            circuit_closures={3: ['edfg', 'acdg', 'bcfg', 'cefh',
        ....:                 'afgh', 'abce', 'abdf', 'begh', 'bcdh', 'adeh'],
        ....:                              4: ['abcdefgh']})
        sage: M.equals(matroids.named_matroids.P8())
        True
    """

    # NECESSARY
    def __init__(self, M=None, groundset=None, circuit_closures=None):
        """
        Initialization of the matroid. See class docstring for full
        documentation.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = CircuitClosuresMatroid(matroids.named_matroids.Fano())
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
            sage: M.equals(matroids.named_matroids.P8())
            True
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

    cpdef groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT:

        A set.

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
        """
        return frozenset(self._groundset)

    cpdef _rank(self, X):
        """
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and
        ``X`` may be assumed to have the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface.

        OUTPUT:

        The rank of ``X`` in the matroid.

        EXAMPLES::

            sage: M = matroids.named_matroids.NonPappus()
            sage: M._rank('abc')
            2
        """
        return len(self._max_independent(X))

    # OPTIONAL, OPTIMIZED FOR THIS CLASS
    cpdef full_rank(self):
        r"""
        Return the rank of the matroid.

        The *rank* of the matroid is the size of the largest independent
        subset of the groundset.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M.full_rank()
            4
            sage: M.dual().full_rank()
            4
        """
        return self._matroid_rank

    cpdef _is_independent(self, F):
        """
        Test if input is independent.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: M._is_independent(set(['a', 'b', 'c']))
            True
            sage: M._is_independent(set(['a', 'b', 'c', 'd']))
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

    cpdef _max_independent(self, F):
        """
        Compute a maximal independent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        A maximal independent subset of ``X``.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: X = M._max_independent(set(['a', 'c', 'd', 'e', 'f']))
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

    cpdef _circuit(self, F):
        """
        Return a minimal dependent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        A circuit contained in ``X``, if it exists. Otherwise an error is
        raised.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: sorted(M._circuit(set(['a', 'c', 'd', 'e', 'f'])))
            ['c', 'd', 'e', 'f']
            sage: sorted(M._circuit(set(['a', 'c', 'd'])))
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

    cpdef circuit_closures(self):
        """
        Return the list of closures of circuits of the matroid.

        A *circuit closure* is a closed set containing a circuit.

        OUTPUT:

        A dictionary containing the circuit closures of the matroid, indexed
        by their ranks.

        .. SEEALSO::

            :meth:`Matroid.circuit() <sage.matroids.matroid.Matroid.circuit>`,
            :meth:`Matroid.closure() <sage.matroids.matroid.Matroid.closure>`

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = CircuitClosuresMatroid(matroids.named_matroids.Fano())
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

        - ``other`` -- A matroid,
        - optional parameter ``certificate`` -- Boolean.

        OUTPUT:

        Boolean,
        and, if certificate = True, a dictionary giving the isomorphism or None

        .. NOTE::

            Internal version that does no input checking.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M1 = CircuitClosuresMatroid(matroids.Wheel(3))
            sage: M2 = matroids.CompleteGraphic(4)
            sage: M1._is_isomorphic(M2)
            True
            sage: M1._is_isomorphic(M2, certificate=True)
            (True, {0: 0, 1: 1, 2: 2, 3: 3, 4: 5, 5: 4})
            sage: M1 = CircuitClosuresMatroid(matroids.named_matroids.Fano())
            sage: M2 = matroids.named_matroids.NonFano()
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
        SM = SetSystem(list(self.groundset()))
        for r in self._circuit_closures:
            for C in self._circuit_closures[r]:
                SM.append(C)
        SN = SetSystem(list(N.groundset()))
        for r in N._circuit_closures:
            for C in N._circuit_closures[r]:
                SN.append(C)
        return SM._isomorphism(SN) is not None

    # REPRESENTATION
    def _repr_(self):
        """
        Return a string representation of the matroid.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: print(M._repr_())
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
                 {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'},
                 {'e', 'f', 'g', 'h'}},
             4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        """
        return Matroid._repr_(self) + " with circuit-closures\n" + setprint_s(self._circuit_closures)

    # COMPARISON

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

            sage: M = matroids.named_matroids.Vamos()
            sage: N = matroids.named_matroids.Vamos()
            sage: hash(M) == hash(N)
            True
            sage: O = matroids.named_matroids.NonVamos()
            sage: hash(M) == hash(O)
            False
        """
        return hash(tuple([self.groundset(), tuple([(r, len(self._circuit_closures[r])) for r in sorted(self._circuit_closures.keys())])]))

    def __richcmp__(left, right, int op):
        r"""
        Compare two matroids.

        We take a very restricted view on equality: the objects need to be of
        the exact same type (so no subclassing) and the internal data need to
        be the same. For BasisMatroids, this means that the groundsets and the
        sets of bases of the two matroids are equal.

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus()
            sage: N = matroids.named_matroids.NonPappus()
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

    # COPYING, LOADING, SAVING

    def __copy__(self):
        """
        Create a shallow copy.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: N = copy(M)  # indirect doctest
            sage: M == N
            True
            sage: M.groundset() is N.groundset()
            True

        """
        N = CircuitClosuresMatroid(groundset=[], circuit_closures={})
        N._groundset = self._groundset
        N._circuit_closures = self._circuit_closures
        N._matroid_rank = self._matroid_rank
        if getattr(self, '__custom_name') is not None:  # because of name wrangling, this is not caught by the default copy
            N.rename(getattr(self, '__custom_name'))
        return N

    def __deepcopy__(self, memo={}):
        """
        Create a deep copy.

        .. NOTE::

            Since matroids are immutable, a shallow copy normally suffices.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos()
            sage: N = deepcopy(M)  # indirect doctest
            sage: M == N
            True
            sage: M.groundset() is N.groundset()
            False
        """
        from copy import deepcopy
        # Since matroids are immutable, N cannot reference itself in correct code, so no need to worry about the recursion.
        N = CircuitClosuresMatroid(groundset=deepcopy(self._groundset, memo), circuit_closures=deepcopy(self._circuit_closures, memo))
        if getattr(self, '__custom_name') is not None:  # because of name wrangling, this is not caught by the default deepcopy
            N.rename(deepcopy(getattr(self, '__custom_name'), memo))
        return N

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

            sage: M = matroids.named_matroids.Vamos()
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
        data = (self._groundset, self._circuit_closures, getattr(self, '__custom_name'))
        version = 0
        return sage.matroids.unpickling.unpickle_circuit_closures_matroid, (version, data)

# todo: customized minor, extend methods.
