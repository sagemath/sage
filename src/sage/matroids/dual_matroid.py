r"""
Dual matroids

Theory
======

Let `M` be a matroid with ground set `E`. If `B` is the set of bases of `M`,
then the set `\{E - b : b \in B\}` is the set of bases of another matroid, the
dual of `M`.

EXAMPLES::

    sage: M = matroids.named_matroids.Fano()
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
dual (such as :class:`RankMatroid <sage.matroids.rank_matroid.RankMatroid>`
and
:class:`CircuitClosuresMatroid <sage.matroids.circuit_closures_matroid.CircuitClosuresMatroid>`).
It is also used as default implementation of the method
:meth:`M.dual() <sage.matroids.matroid.Matroid.dual>`.
For direct access to the ``DualMatroid`` constructor, run::

    sage: from sage.matroids.advanced import *

See also :mod:`sage.matroids.advanced`.


AUTHORS:

- Rudi Pendavingh, Michael Welsh, Stefan van Zwam (2013-04-01): initial version

Methods
=======
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

from .matroid import Matroid


class DualMatroid(Matroid):
    r"""
    Dual of a matroid.

    For some matroid representations it can be computationally expensive to
    derive an explicit representation of the dual. This class wraps around any
    matroid to provide an abstract dual. It also serves as default
    implementation.

    INPUT:

    - ``matroid`` - a matroid.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.named_matroids.Vamos()
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
            sage: M = matroids.named_matroids.Vamos()
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
        """
        if not isinstance(matroid, Matroid):
            raise TypeError("no matroid provided to take dual of.")
        self._matroid = matroid

    def groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT:

        A set.

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus().dual()
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

        - ``X`` -- an object with Python's ``frozenset`` interface.

        OUTPUT:

        The rank of ``X`` in the matroid.

        EXAMPLES::

            sage: M = matroids.named_matroids.NonPappus().dual()
            sage: M._rank(['a', 'b', 'c'])
            3

        """
        return self._matroid._corank(X)

    def _corank(self, X):
        """
        Return the corank of a set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``.

        OUTPUT:

        The corank of ``X``.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos().dual()
            sage: M._corank(set(['a', 'e', 'g', 'd', 'h']))
            4
        """
        return self._matroid._rank(X)

    def _max_independent(self, X):
        """
        Compute a maximal independent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        A maximal independent subset of ``X``.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos().dual()
            sage: X = M._max_independent(set(['a', 'c', 'd', 'e', 'f']))
            sage: sorted(X) # random
            ['a', 'c', 'd', 'e']
            sage: M.is_independent(X)
            True
            sage: all(M.is_dependent(X.union([y])) for y in M.groundset() if y not in X)
            True
        """
        return self._matroid._max_coindependent(X)

    def _circuit(self, X):
        """
        Return a minimal dependent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        A circuit contained in ``X``, if it exists. Otherwise an error is
        raised.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos().dual()
            sage: sorted(sage.matroids.matroid.Matroid._circuit(M,
            ....:                             set(['a', 'c', 'd', 'e', 'f'])))
            ['c', 'd', 'e', 'f']
            sage: sorted(sage.matroids.matroid.Matroid._circuit(M,
            ....:                                       set(['a', 'c', 'd'])))
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set.

        """
        return self._matroid._cocircuit(X)

    def _closure(self, X):
        """
        Return the closure of a set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        The smallest closed set containing ``X``.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos().dual()
            sage: sorted(M._closure(set(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']

        """
        return self._matroid._coclosure(X)

    def _max_coindependent(self, X):
        """
        Compute a maximal coindependent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        A maximal coindependent subset of ``X``.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos().dual()
            sage: X = M._max_coindependent(set(['a', 'c', 'd', 'e', 'f']))
            sage: sorted(X) # random
            ['a', 'd', 'e', 'f']
            sage: M.is_coindependent(X)
            True
            sage: all(M.is_codependent(X.union([y])) for y in M.groundset() if y not in X)
            True
        """
        return self._matroid._max_independent(X)

    def _coclosure(self, X):
        """
        Return the coclosure of a set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        The smallest coclosed set containing ``X``.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos().dual()
            sage: sorted(M._coclosure(set(['a', 'b', 'c'])))
            ['a', 'b', 'c', 'd']

        """
        return self._matroid._closure(X)

    def _cocircuit(self, X):
        """
        Return a minimal codependent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        A cocircuit contained in ``X``, if it exists. Otherwise an error is
        raised.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos().dual()
            sage: sorted(M._cocircuit(set(['a', 'c', 'd', 'e', 'f'])))
            ['c', 'd', 'e', 'f']
            sage: sorted(M._cocircuit(set(['a', 'c', 'd'])))
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set
        """
        return self._matroid._circuit(X)

    def _minor(self, contractions=None, deletions=None):
        r"""
        Return a minor.

        INPUT:

        - ``contractions`` -- An object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``.
        - ``deletions`` -- An object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``.

        OUTPUT:

        A ``DualMatroid`` instance representing
        `(``self._matroid`` / ``deletions`` \ ``contractions``)^*`

        .. NOTE::

            This method does NOT do any checks. Besides the assumptions above,
            we assume the following:

            - ``contractions`` is independent
            - ``deletions`` is coindependent
            - ``contractions`` and ``deletions`` are disjoint.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos().dual()
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

        Let `M` be a matroid with ground set `E`. If `B` is the set of bases
        of `M`, then the set `\{E - b : b \in B\}` is the set of bases of
        another matroid, the *dual* of `M`. Note that the dual of the dual of
        `M` equals `M`, so if this is the
        :class:`DualMatroid` instance
        wrapping `M` then the returned matroid is `M`.

        OUTPUT:

        The dual matroid.

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus().dual()
            sage: N = M.dual()
            sage: N.rank()
            3
            sage: N
            Pappus: Matroid of rank 3 on 9 elements with circuit-closures
            {2: {{'a', 'b', 'c'}, {'a', 'e', 'i'}, {'a', 'f', 'h'},
                 {'b', 'd', 'i'}, {'b', 'f', 'g'}, {'c', 'd', 'h'},
                 {'c', 'e', 'g'}, {'d', 'e', 'f'}, {'g', 'h', 'i'}},
             3: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'}}}
        """
        return self._matroid

    # REPRESENTATION
    def _repr_(self):
        """
        Return a string representation of the matroid.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos().dual()
            sage: print(M._repr_())
            Dual of 'Vamos:
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
                 {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'},
                 {'e', 'f', 'g', 'h'}},
             4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}'
        """
        return "Dual of '" + repr(self._matroid) + "'"

    # COMPARISON

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

            sage: M = matroids.named_matroids.Vamos().dual()
            sage: N = matroids.named_matroids.Vamos().dual()
            sage: O = matroids.named_matroids.Vamos()
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

        - ``other`` -- A matroid.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M1 = matroids.named_matroids.Fano()
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

        - ``other`` -- A matroid.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M1 = matroids.named_matroids.Fano()
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
        N = DualMatroid(self._matroid)
        if getattr(self, '__custom_name') is not None:
            # because of name wrangling, this is not caught by the default copy
            N.rename(getattr(self, '__custom_name'))
        return N

    def __deepcopy__(self, memo={}):
        """
        Create a deep copy.

        .. NOTE::

            Since matroids are immutable, a shallow copy normally suffices.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos().dual()
            sage: N = deepcopy(M)  # indirect doctest
            sage: M == N
            True
            sage: M.groundset() is N.groundset()
            False
        """
        from copy import deepcopy
        N = DualMatroid(deepcopy(self._matroid, memo))
        if getattr(self, '__custom_name') is not None:
            # because of name wrangling, this is not caught by the
            # default deepcopy
            N.rename(deepcopy(getattr(self, '__custom_name'), memo))
        return N

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

            sage: M = matroids.named_matroids.Vamos().dual()
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
        data = (self._matroid, getattr(self, '__custom_name'))
        version = 0
        return sage.matroids.unpickling.unpickle_dual_matroid, (version, data)
