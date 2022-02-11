r"""
Minors of matroids

Theory
======

Let `M` be a matroid with ground set `E`. There are two standard ways to
remove an element from `E` so that the result is again a matroid, *deletion*
and *contraction*. Deletion is simply omitting the elements from a set `D`
from `E` and keeping all remaining independent sets. This is denoted ``M \ D``
(this also works in Sage). Contraction is the dual operation:
``M / C == (M.dual() \ C).dual()``.

EXAMPLES::

    sage: M = matroids.named_matroids.Fano()
    sage: M \ ['a', 'c' ] == M.delete(['a', 'c'])
    True
    sage: M / 'a' == M.contract('a')
    True
    sage: M / 'c' \ 'ab' == M.minor(contractions='c', deletions='ab')
    True

If a contraction set is not independent (or a deletion set not coindependent),
this is taken care of::

    sage: M = matroids.named_matroids.Fano()
    sage: M.rank('abf')
    2
    sage: M / 'abf' == M / 'ab' \ 'f'
    True
    sage: M / 'abf' == M / 'af' \ 'b'
    True

.. SEEALSO::

    :meth:`M.delete() <sage.matroids.matroid.Matroid.delete>`,
    :meth:`M.contract() <sage.matroids.matroid.Matroid.contract>`,
    :meth:`M.minor() <sage.matroids.matroid.Matroid.minor>`,

Implementation
==============

The class :class:`MinorMatroid <sage.matroids.minor_matroid.MinorMatroid>`
wraps around a matroid instance to represent a minor. Only useful for classes
that don't have an explicit construction of minors
(such as :class:`RankMatroid <sage.matroids.rank_matroid.RankMatroid>` and
:class:`CircuitClosuresMatroid <sage.matroids.circuit_closures_matroid.CircuitClosuresMatroid>`).
It is also used as default implementation of the minor methods
:meth:`M.minor(C, D) <sage.matroids.matroid.Matroid.minor>`,
:meth:`M.delete(D) <sage.matroids.matroid.Matroid.delete>`,
:meth:`M.contract(C) <sage.matroids.matroid.Matroid.contract>`.
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
from .utilities import setprint_s


class MinorMatroid(Matroid):
    r"""
    Minor of a matroid.

    For some matroid representations, it can be computationally
    expensive to derive an explicit representation of a minor. This
    class wraps around any matroid to provide an abstract minor. It
    also serves as default implementation.

    Return a minor.

    INPUT:

    - ``matroid`` -- a matroid.
    - ``contractions`` -- An object with Python's ``frozenset`` interface
      containing a subset of ``self.groundset()``.
    - ``deletions`` -- An object with Python's ``frozenset`` interface
      containing a subset of ``self.groundset()``.

    OUTPUT:

    A ``MinorMatroid`` instance representing
    ``matroid / contractions \ deletions``.

    .. WARNING::

        This class does NOT do any checks. Besides the assumptions above, we
        assume the following:

        - ``contractions`` is independent
        - ``deletions`` is coindependent
        - ``contractions`` and ``deletions`` are disjoint.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.named_matroids.Vamos()
        sage: N = MinorMatroid(matroid=M, contractions=set(['a']),
        ....:                  deletions=set())
        sage: N._minor(contractions=set(), deletions=set(['b', 'c']))
        M / {'a'} \ {'b', 'c'}, where M is Vamos:
        Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
             {'c', 'd', 'e', 'f'}, {'e', 'f', 'g', 'h'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
    """

    def __init__(self, matroid, contractions=None, deletions=None):
        """
        See class docstring for documentation.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = MinorMatroid(matroids.named_matroids.Fano(),
            ....:                  contractions=set(), deletions=set(['g']))  # indirect doctest
            sage: M.is_isomorphic(matroids.Wheel(3))
            True
        """
        if not isinstance(matroid, Matroid):
            raise TypeError("no matroid provided to take minor of.")
        self._matroid = matroid
        self._contractions = frozenset(contractions)
        self._deletions = frozenset(deletions)
        self._delsize = len(self._deletions)
        self._consize = len(self._contractions)
        self._groundset = matroid.groundset().difference(self._deletions.union(self._contractions))

    def groundset(self):
        """
        Return the groundset of the matroid.

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus().contract(['c'])
            sage: sorted(M.groundset())
            ['a', 'b', 'd', 'e', 'f', 'g', 'h', 'i']
        """
        return self._groundset

    def _rank(self, X):
        """
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and
        ``X`` may be assumed to have the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface.

        OUTPUT:

        The rank of ``X`` in the matroid.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = MinorMatroid(matroids.named_matroids.NonPappus(),
            ....:                  contractions=set(), deletions={'f', 'g'})
            sage: M._rank(frozenset('abc'))
            2
        """
        return self._matroid._rank(self._contractions.union(X)) - self._consize

    def _corank(self, X):
        """
        Return the corank of a set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        The corank of ``X``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = MinorMatroid(matroids.named_matroids.Vamos(),
            ....:                 contractions=set('c'), deletions={'b', 'f'})
            sage: M._corank(set(['a', 'e', 'g', 'd', 'h']))
            2
        """
        return self._matroid._corank(self._deletions.union(X)) - self._delsize

    def _max_independent(self, X):
        """
        Compute a maximal independent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        A maximal independent subset of ``X``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = MinorMatroid(matroids.named_matroids.Vamos(),
            ....:                 contractions=set('c'), deletions={'b', 'f'})
            sage: X = M._max_independent(set(['a', 'd', 'e', 'g']))
            sage: sorted(X) # random
            ['a', 'd', 'e']
            sage: M.is_independent(X)
            True
            sage: all(M.is_dependent(X.union([y])) for y in M.groundset() if y not in X)
            True
        """
        return self._matroid._augment(self._contractions, X)

    def _closure(self, X):
        """
        Return the closure of a set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        The smallest closed set containing ``X``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = MinorMatroid(matroids.named_matroids.Vamos(),
            ....:                 contractions=set('c'), deletions={'b', 'f'})
            sage: sorted(M._closure(set(['a', 'e', 'd'])))
            ['a', 'd', 'e', 'g', 'h']

        """
        return self._matroid._closure(self._contractions.union(X)).difference(self._contractions.union(self._deletions))

    def _max_coindependent(self, X):
        """
        Compute a maximal coindependent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        A maximal coindependent subset of ``X``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = MinorMatroid(matroids.named_matroids.Vamos(),
            ....:                 contractions=set('c'), deletions={'b', 'f'})
            sage: X = M._max_coindependent(set(['a', 'd', 'e', 'g']))
            sage: sorted(X) # random
            ['d', 'g']
            sage: M.is_coindependent(X)
            True
            sage: all(M.is_codependent(X.union([y])) for y in M.groundset() if y not in X)
            True
        """
        return X - self._matroid._augment(self._contractions.union(self._groundset - X), X)

    def _coclosure(self, X):
        """
        Return the coclosure of a set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        The smallest coclosed set containing ``X``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = MinorMatroid(matroids.named_matroids.Vamos(),
            ....:                 contractions=set('c'), deletions={'b', 'f'})
            sage: sorted(M._coclosure(set(['a', 'b', 'c'])))
            ['a', 'd', 'e', 'g', 'h']

        """
        return self._matroid._coclosure(self._deletions.union(X)).difference(self._contractions.union(self._deletions))

    def _minor(self, contractions, deletions):
        r"""
        Return a minor.

        INPUT:

        - ``contractions`` -- An object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``.
        - ``deletions`` -- An object with Python's ``frozenset`` interface
          containing a subset of ``self.groundset()``.

        OUTPUT:

        A ``MinorMatroid`` instance representing
        `(``self._matroid`` / ``deletions`` \ ``contractions``)^*`

        .. NOTE::

            This method does NOT do any checks. Besides the assumptions above, we assume the following:

            - ``contractions`` is independent
            - ``deletions`` is coindependent
            - ``contractions`` and ``deletions`` are disjoint.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = MinorMatroid(matroids.named_matroids.Vamos(), contractions=set('c'), deletions={'b', 'f'})
            sage: N = M._minor(contractions=set(['a']), deletions=set([]))
            sage: N._minor(contractions=set([]), deletions=set(['d']))
            M / {'a', 'c'} \ {'b', 'd', 'f'}, where M is Vamos:
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
                 {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'},
                 {'e', 'f', 'g', 'h'}},
             4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        """
        return MinorMatroid(self._matroid, self._contractions.union(contractions), self._deletions.union(deletions))

    # representation

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
        s = "M"
        if self._contractions:
            s += r" / " + setprint_s(self._contractions, toplevel=True)
        if self._deletions:
            s += r" \ " + setprint_s(self._deletions, toplevel=True)
        s += ", where M is " + repr(self._matroid)
        return s

    # Comparison:

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

            sage: from sage.matroids.advanced import *
            sage: M = MinorMatroid(matroids.named_matroids.Vamos(),
            ....:                 contractions=set('c'), deletions={'b', 'f'})
            sage: N = MinorMatroid(matroids.named_matroids.Vamos(),
            ....:                 deletions={'b', 'f'}, contractions=set('c'))
            sage: O = MinorMatroid(matroids.named_matroids.Vamos(),
            ....:                 contractions={'b', 'f'}, deletions=set('c'))
            sage: hash(M) == hash(N)
            True
            sage: hash(M) == hash(O)
            False
        """
        return hash((self._matroid, self._contractions, self._deletions))

    def __eq__(self, other):
        """
        Compare two matroids.

        INPUT:

        - ``other`` -- A matroid.

        OUTPUT:

        ``True`` if ``self`` and ``other`` have the same underlying matroid,
        same set of contractions, and same set of deletions; ``False``
        otherwise.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = matroids.named_matroids.Fano()
            sage: M1 = MinorMatroid(M, set('ab'), set('f'))
            sage: M2 = MinorMatroid(M, set('af'), set('b'))
            sage: M3 = MinorMatroid(M, set('a'), set('f'))._minor(set('b'), set())
            sage: M1 == M2  # indirect doctest
            False
            sage: M1.equals(M2)
            True
            sage: M1 == M3
            True
        """
        if not isinstance(other, MinorMatroid):
            return False
        return (self._contractions == other._contractions) and (self._deletions == other._deletions) and (self._matroid == other._matroid)

    def __ne__(self, other):
        """
        Compare two matroids.

        INPUT:

        - ``other`` -- A matroid.

        OUTPUT:

        ``False`` if ``self`` and ``other`` have the same underlying matroid,
        same set of contractions, and same set of deletions; ``True``
        otherwise.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = matroids.named_matroids.Fano()
            sage: M1 = MinorMatroid(M, set('ab'), set('f'))
            sage: M2 = MinorMatroid(M, set('af'), set('b'))
            sage: M3 = MinorMatroid(M, set('a'), set('f'))._minor(set('b'), set())
            sage: M1 != M2  # indirect doctest
            True
            sage: M1.equals(M2)
            True
            sage: M1 != M3
            False
        """
        return not self == other

    # Copying, loading, saving:

    def __copy__(self):
        """
        Create a shallow copy.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = MinorMatroid(matroid=matroids.named_matroids.Vamos(),
            ....:                  contractions={'a', 'b'}, deletions={'f'})
            sage: N = copy(M)  # indirect doctest
            sage: M == N
            True
            sage: M._matroid is N._matroid
            True
        """
        N = MinorMatroid(self._matroid, self._contractions, self._deletions)
        if getattr(self, '__custom_name') is not None:  # because of name wrangling, this is not caught by the default copy
            N.rename(getattr(self, '__custom_name'))
        return N

    def __deepcopy__(self, memo={}):
        """
        Create a deep copy.

        .. NOTE::

            Since matroids are immutable, a shallow copy normally suffices.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = MinorMatroid(matroid=matroids.named_matroids.Vamos(),
            ....:                  contractions={'a', 'b'}, deletions={'f'})
            sage: N = deepcopy(M)  # indirect doctest
            sage: M == N
            True
            sage: M._matroid is N._matroid
            False
        """
        from copy import deepcopy
        # Since matroids are immutable, N cannot reference itself in correct code, so no need to worry about the recursion.
        N = MinorMatroid(deepcopy(self._matroid, memo), deepcopy(self._contractions, memo), deepcopy(self._deletions, memo))
        if getattr(self, '__custom_name') is not None:  # because of name wrangling, this is not caught by the default deepcopy
            N.rename(deepcopy(getattr(self, '__custom_name'), memo))
        return N

    def __reduce__(self):
        r"""
        Save the matroid for later reloading.

        EXAMPLES::

            sage: M = matroids.named_matroids.Vamos().minor('abc', 'g')
            sage: M == loads(dumps(M))  # indirect doctest
            True
            sage: loads(dumps(M))
            M / {'a', 'b', 'c'} \ {'g'}, where M is Vamos:
            Matroid of rank 4 on 8 elements with circuit-closures
            {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'},
                 {'a', 'b', 'g', 'h'}, {'c', 'd', 'e', 'f'},
                 {'e', 'f', 'g', 'h'}},
             4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        """
        import sage.matroids.unpickling
        data = (self._matroid, self._contractions, self._deletions, getattr(self, '__custom_name'))
        version = 0
        return sage.matroids.unpickling.unpickle_minor_matroid, (version, data)
