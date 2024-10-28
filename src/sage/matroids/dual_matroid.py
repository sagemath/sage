r"""
Dual matroids

Theory
======

Let `M` be a matroid with groundset `E`. If `B` is the set of bases of `M`,
then the set `\{E - b : b \in B\}` is the set of bases of another matroid, the
dual of `M`.

EXAMPLES::

    sage: M = matroids.catalog.Fano()
    sage: N = M.dual()
    sage: M.is_basis('abc')
    True
    sage: N.is_basis('defg')
    True
    sage: M.dual().dual() == M
    True

Implementation
==============

The class ``DualMatroid`` wraps around a matroid instance to represent its
dual. Only useful for classes that don't have an explicit construction of the
dual (such as :class:`RankMatroid <sage.matroids.rank_matroid.RankMatroid>` and
:class:`CircuitClosuresMatroid <sage.matroids.circuit_closures_matroid.CircuitClosuresMatroid>`).
It is also used as the default implementation of the method
:meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`.
For direct access to the ``DualMatroid`` constructor, run::

    sage: from sage.matroids.advanced import *

See also :mod:`sage.matroids.advanced`.

AUTHORS:

- Rudi Pendavingh, Michael Welsh, Stefan van Zwam (2013-04-01): initial version
"""

# ****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Michael Welsh <michael@welsh.co.nz>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matroids.matroid import Matroid

class DualMatroid(Matroid):
    r"""
    Dual of a matroid.

    For some matroid representations it can be computationally expensive to
    derive an explicit representation of the dual. This class wraps around any
    matroid to provide an abstract dual. It also serves as the default
    implementation of the dual.

    INPUT:

    - ``matroid`` -- matroid

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.catalog.Vamos()
        sage: Md = DualMatroid(M)  # indirect doctest
        sage: Md.rank('abd') == M.corank('abd')
        True
        sage: Md
        Dual of 'Vamos: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
             {'c', 'd', 'e', 'f'}, {'e', 'f', 'g', 'h'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}'
    """

    def __init__(self, matroid):
        """
        See class definition for documentation.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = matroids.catalog.Vamos()
            sage: Md = DualMatroid(M)  # indirect doctest
            sage: Md.rank('abd') == M.corank('abd')
            True
            sage: Md
            Dual of 'Vamos:
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
                 {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'},
                 {'e', 'f', 'g', 'h'}},
             4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}'

        TESTS::

            sage: from sage.matroids.dual_matroid import DualMatroid
            sage: DualMatroid([])
            Traceback (most recent call last):
            ...
            TypeError: no matroid provided to take the dual of
        """
        if not isinstance(matroid, Matroid):
            raise TypeError("no matroid provided to take the dual of")
        self._matroid = matroid

    def groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT: set

        EXAMPLES::

            sage: M = matroids.catalog.Pappus().dual()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
        """
        return self._matroid.groundset()

    def _rank(self, X):
        r"""
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and
        ``X`` may be assumed to have the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface

        OUTPUT: integer

        EXAMPLES::

            sage: M = matroids.catalog.NonPappus().dual()
            sage: M._rank(frozenset(['a', 'b', 'c']))
            3
        """
        return self._matroid._corank(X)

    def _corank(self, X):
        """
        Return the corank of a set.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``

        OUTPUT: integer

        EXAMPLES::

            sage: M = matroids.catalog.Vamos().dual()
            sage: M._corank(frozenset(['a', 'e', 'g', 'd', 'h']))
            4
        """
        return self._matroid._rank(X)

    def _max_independent(self, X):
        """
        Compute a maximal independent subset.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: a maximal independent subset of ``X``

        EXAMPLES::

            sage: M = matroids.catalog.Vamos().dual()
            sage: X = M._max_independent(frozenset(['a', 'c', 'd', 'e', 'f']))
            sage: sorted(X) # random
            ['a', 'c', 'd', 'e']
            sage: M.is_independent(X)
            True
            sage: all(M.is_dependent(X.union([y])) for y in M.groundset()
            ....:     if y not in X)
            True
        """
        return self._matroid._max_coindependent(X)

    def _circuit(self, X):
        """
        Return a minimal dependent subset.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: a circuit contained in ``X``, if it exists; otherwise, an error is
        raised

        EXAMPLES::

            sage: M = matroids.catalog.Vamos().dual()
            sage: sorted(sage.matroids.matroid.Matroid._circuit(M,
            ....:                             frozenset(['a', 'c', 'd', 'e', 'f'])))
            ['c', 'd', 'e', 'f']
            sage: sorted(sage.matroids.matroid.Matroid._circuit(M,
            ....:                                       frozenset(['a', 'c', 'd'])))
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set.
        """
        return self._matroid._cocircuit(X)

    def _closure(self, X):
        """
        Return the closure of a set.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: the smallest closed set containing ``X``

        EXAMPLES::

            sage: M = matroids.catalog.Vamos().dual()
            sage: sorted(M._closure(frozenset(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']
        """
        return self._matroid._coclosure(X)

    def _max_coindependent(self, X):
        """
        Compute a maximal coindependent subset.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: a maximal coindependent subset of ``X``

        EXAMPLES::

            sage: M = matroids.catalog.Vamos().dual()
            sage: X = M._max_coindependent(frozenset(['a', 'c', 'd', 'e', 'f']))
            sage: sorted(X) # random
            ['a', 'd', 'e', 'f']
            sage: M.is_coindependent(X)
            True
            sage: all(M.is_codependent(X.union([y])) for y in M.groundset()
            ....:     if y not in X)
            True
        """
        return self._matroid._max_independent(X)

    def _coclosure(self, X):
        """
        Return the coclosure of a set.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: the smallest coclosed set containing ``X``

        EXAMPLES::

            sage: M = matroids.catalog.Vamos().dual()
            sage: sorted(M._coclosure(frozenset(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']
        """
        return self._matroid._closure(X)

    def _cocircuit(self, X):
        """
        Return a minimal codependent subset.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``

        OUTPUT: a cocircuit contained in ``X``, if it exists; otherwise, an error is
        raised

        EXAMPLES::

            sage: M = matroids.catalog.Vamos().dual()
            sage: sorted(M._cocircuit(frozenset(['a', 'c', 'd', 'e', 'f'])))
            ['c', 'd', 'e', 'f']
            sage: sorted(M._cocircuit(frozenset(['a', 'c', 'd'])))
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set
        """
        return self._matroid._circuit(X)

    def _minor(self, contractions=None, deletions=None):
        r"""
        Return a minor.

        INPUT:

        - ``contractions`` -- an object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``
        - ``deletions`` -- an object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``

        OUTPUT: a :class:`DualMatroid` representing
        `(``self._matroid`` / ``deletions`` \ ``contractions``)^*`

        .. NOTE::

            This method does NOT do any checks. Besides the assumptions above,
            we assume the following:

            - ``contractions`` is independent;
            - ``deletions`` is coindependent;
            - ``contractions`` and ``deletions`` are disjoint.

        EXAMPLES::

            sage: M = matroids.catalog.Vamos().dual()
            sage: N = M._minor(contractions=set(['a']), deletions=set([]))
            sage: N._minor(contractions=set([]), deletions=set(['b', 'c']))
            Dual of 'M / {'b', 'c'} \ {'a'}, where M is Vamos:
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
                 {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'},
                 {'e', 'f', 'g', 'h'}},
             4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}'
        """
        # Assumption: if self._matroid cannot make a dual, neither can
        # its minor.
        return DualMatroid(self._matroid._minor(contractions=deletions,
                                                deletions=contractions))

    def dual(self):
        r"""
        Return the dual of the matroid.

        Let `M` be a matroid with groundset `E`. If `B` is the set of bases
        of `M`, then the set `\{E - b : b \in B\}` is the set of bases of
        another matroid, the *dual* of `M`. Note that the dual of the dual of
        `M` equals `M`, so if this is the :class:`DualMatroid` instance
        wrapping `M` then the returned matroid is `M`.

        OUTPUT: the dual matroid

        EXAMPLES::

            sage: M = matroids.catalog.Pappus().dual()
            sage: N = M.dual()
            sage: N.rank()
            3
            sage: N
            Pappus: Matroid of rank 3 on 9 elements with 9 nonspanning circuits
        """
        return self._matroid

    # representation

    def _repr_(self):
        """
        Return a string representation of the matroid.

        EXAMPLES::

            sage: M = matroids.catalog.Vamos().dual()
            sage: print(M._repr_())
            Dual of 'Vamos:
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
                 {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'},
                 {'e', 'f', 'g', 'h'}},
             4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}'
        """
        return "Dual of '" + repr(self._matroid) + "'"

    # comparison

    def __hash__(self):
        r"""
        Return an invariant of the matroid.

        This function is called when matroids are added to a set. It is very
        desirable to override it so it can distinguish matroids on the same
        groundset, which is a very typical use case!

        .. WARNING::

            This method is linked to __richcmp__ (in Cython) and __cmp__ or
            __eq__/__ne__ (in Python). If you override one, you should (and in
            Cython: MUST) override the other!

        EXAMPLES::

            sage: M = matroids.catalog.Vamos().dual()
            sage: N = matroids.catalog.Vamos().dual()
            sage: O = matroids.catalog.Vamos()
            sage: hash(M) == hash(N)
            True
            sage: hash(M) == hash(O)
            False
        """
        return hash(("Dual", hash(self._matroid)))

    def __eq__(self, other):
        """
        Compare two matroids.

        INPUT:

        - ``other`` -- matroid

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M1 = matroids.catalog.Fano()
            sage: M2 = CircuitClosuresMatroid(M1.dual())
            sage: M3 = CircuitClosuresMatroid(M1).dual()
            sage: M4 = CircuitClosuresMatroid(groundset='abcdefg',
            ....:   circuit_closures={3: ['abcdefg'], 2: ['beg', 'cdb', 'cfg',
            ....:   'ace', 'fed', 'gad', 'fab']}).dual()
            sage: M1.dual() == M2  # indirect doctest
            False
            sage: M2 == M3
            False
            sage: M3 == M4
            True
        """
        if not isinstance(other, DualMatroid):
            return False
        return self._matroid == other._matroid

    def __ne__(self, other):
        """
        Compare two matroids.

        INPUT:

        - ``other`` -- matroid

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M1 = matroids.catalog.Fano()
            sage: M2 = CircuitClosuresMatroid(M1.dual())
            sage: M3 = CircuitClosuresMatroid(M1).dual()
            sage: M4 = CircuitClosuresMatroid(groundset='abcdefg',
            ....:   circuit_closures={3: ['abcdefg'], 2: ['beg', 'cdb', 'cfg',
            ....:                'ace', 'fed', 'gad', 'fab']}).dual()
            sage: M1.dual() != M2  # indirect doctest
            True
            sage: M2 != M3
            True
            sage: M3 != M4
            False
        """
        return not self == other

    # copying, loading, saving

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        OUTPUT:

        A tuple ``(unpickle, (version, data))``, where ``unpickle`` is the
        name of a function that, when called with ``(version, data)``,
        produces a matroid isomorphic to ``self``. ``version`` is an integer
        (currently 0) and ``data`` is a tuple ``(M, name)`` where ``M`` is
        the internal matroid, and ``name`` is a custom name.

        EXAMPLES::

            sage: M = matroids.catalog.Vamos().dual()
            sage: M == loads(dumps(M))  # indirect doctest
            True
            sage: loads(dumps(M))
            Dual of 'Vamos:
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
                 {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'},
                 {'e', 'f', 'g', 'h'}},
             4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}'
        """
        import sage.matroids.unpickling
        data = (self._matroid, self.get_custom_name())
        version = 0
        return sage.matroids.unpickling.unpickle_dual_matroid, (version, data)

    def relabel(self, mapping):
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

            sage: M = matroids.catalog.K5dual(range(10))
            sage: type(M)
            <class 'sage.matroids.dual_matroid.DualMatroid'>
            sage: sorted(M.groundset())
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: N = M.dual().relabel({0:10})
            sage: sorted(N.groundset())
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            sage: N.is_isomorphic(matroids.catalog.K5())
            True

        TESTS::

            sage: M = matroids.catalog.K5dual(range(10))
            sage: f = {0: 'a', 1: 'b', 2: 'c', 3: 'd', 4: 'e',
            ....:      5: 'f', 6: 'g', 7: 'h', 8: 'i', 9: 'j'}
            sage: N = M.relabel(f)
            sage: for S in powerset(M.groundset()):
            ....:     assert M.rank(S) == N.rank([f[x] for x in S])
        """
        M = self._matroid.relabel(mapping).dual()
        return M

    def is_valid(self):
        """
        Test if ``self`` obeys the matroid axioms.

        For a :class:`DualMatroid`, we check its dual.

        OUTPUT: boolean

        EXAMPLES::

            sage: M = matroids.catalog.K5dual()
            sage: type(M)
            <class 'sage.matroids.dual_matroid.DualMatroid'>
            sage: M.is_valid()
            True
            sage: M = Matroid([[0, 1], [2, 3]])
            sage: M.dual().is_valid()
            False
        """
        return self._matroid.is_valid()
