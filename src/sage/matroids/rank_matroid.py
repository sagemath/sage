r"""
Rank function matroids

The easiest way to define arbitrary matroids in Sage might be through the
class ``RankMatroid``. All that is required is a groundset and a function that
computes the rank for each given subset.

Of course, since the rank function is used as black box, matroids so defined
cannot take advantage of any extra structure your class might have, and rely
on default implementations. Besides this, matroids in this class can't be
saved.

Constructions
=============
Any function can be used, but no checks are performed, so be careful.

EXAMPLES::

    sage: def f(X):
    ....:     return min(len(X), 3)
    sage: M = Matroid(groundset=range(6), rank_function=f)
    sage: M.is_valid()
    True
    sage: M.is_isomorphic(matroids.Uniform(3, 6))
    True

    sage: def g(X):
    ....:     if len(X) >= 3:
    ....:         return 1
    ....:     else:
    ....:         return 0
    sage: N = Matroid(groundset='abc', rank_function=g)
    sage: N.is_valid()
    False

See :class:`below <sage.matroids.rank_matroid.RankMatroid>` for more. It is
recommended to use the :func:`Matroid() <sage.matroids.constructor.Matroid>`
function for easy construction of a ``RankMatroid``. For direct access to the
``RankMatroid`` constructor, run::

        sage: from sage.matroids.advanced import *

AUTHORS:

- Rudi Pendavingh, Stefan van Zwam (2013-04-01): initial version

Methods
=======
"""
# ****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.matroids.matroid import Matroid


class RankMatroid(Matroid):
    """
    Matroid specified by its rank function.

    INPUT:

    - ``groundset`` -- the groundset of a matroid
    - ``rank_function`` -- a function mapping subsets of ``groundset`` to
      nonnegative integers

    OUTPUT: a matroid on ``groundset`` whose rank function equals ``rank_function``

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: def f(X):
        ....:     return min(len(X), 3)
        sage: M = RankMatroid(groundset=range(6), rank_function=f)
        sage: M.is_valid()
        True
        sage: M.is_isomorphic(matroids.Uniform(3, 6))
        True
    """
    def __init__(self, groundset, rank_function):
        """
        Initialize the rank matroid.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = RankMatroid(range(6),
            ....:                 rank_function=matroids.Uniform(3, 6).rank)
            sage: M
            Matroid of rank 3 on 6 elements
        """
        self._groundset = frozenset(groundset)
        self._rank_function = rank_function

    def groundset(self):
        r"""
        Return the groundset of ``self``.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = RankMatroid(range(6),
            ....:                 rank_function=matroids.Uniform(3, 6).rank)
            sage: sorted(M.groundset())
            [0, 1, 2, 3, 4, 5]
        """
        return self._groundset

    def _rank(self, X):
        r"""
        Return the rank of set `X`.

        Internal function without any sanity checks. May assume that ``X``
        has Python's ``frozenset`` interface and is a subset of
        self.groundset().

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = RankMatroid(range(6),
            ....:                 rank_function=matroids.Uniform(3, 6).rank)
            sage: M._rank([0, 2, 4, 5])
            3
        """
        return self._rank_function(X)

    # Comparison:

    def __hash__(self):
        r"""
        Return a string invariant of the matroid.

        This function is called when matroids are added to a set. It is very
        desirable to override it so it can distinguish matroids on the same
        groundset, which is a very typical use case!

        .. WARNING::

            This method is linked to ``__richcmp__`` (in Cython) and ``__cmp__``
            or ``__eq__``/``__ne__`` (in Python). If you override one, you
            should (and, in Cython, \emph{must}) override the other!

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: M = Matroid(groundset=range(10),
            ....:             rank_function=lambda X: min(len(X), 4))
            sage: N = Matroid(groundset=range(10),
            ....:             rank_function=lambda X: min(len(X), 4))
            sage: O = Matroid(groundset=range(10),
            ....:             rank_function=lambda X: min(len(X), 3))
            sage: hash(M) == hash(N)
            True
            sage: hash(M) == hash(O)
            False
        """
        return hash((self.groundset(), self.full_rank()))

    def __eq__(self, other):
        """
        Compare two matroids.

        INPUT:

        - ``other`` -- matroid

        OUTPUT:

        ``True`` if ``self`` and ``other have the same groundset and the same
        rank function; ``False`` otherwise.

        .. NOTE::

            Note that rank functions ``f`` and ``g`` are normally deemed equal
            only if ``f is g``. It would be too time-consuming to check all
            their values.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: def f(X):
            ....:     return min(len(X), 3)
            sage: def g(X):
            ....:     return min(len(X), 3)
            sage: M1 = RankMatroid(groundset=range(6), rank_function=f)
            sage: M2 = RankMatroid(groundset=range(6), rank_function=g)
            sage: M3 = RankMatroid(groundset=range(7), rank_function=f)
            sage: M4 = RankMatroid(groundset=range(6), rank_function=f)
            sage: M1 == M2  # indirect doctest
            False
            sage: M1 == M3
            False
            sage: M1 == M4
            True
        """
        if not isinstance(other, RankMatroid):
            return False
        return (self.groundset() == other.groundset()) and (self._rank_function == other._rank_function)

    def __ne__(self, other):
        """
        Compare two matroids.

        INPUT:

        - ``other`` -- matroid

        OUTPUT:

        ``False`` if ``self`` and ``other have the same groundset and the
        same rank function; ``True`` otherwise.

        .. NOTE::

            Rank functions ``f`` and ``g`` are normally deemed equal only if
            ``f is g``. It would be too time-consuming to check all their
            values.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: def f(X):
            ....:     return min(len(X), 3)
            sage: def g(X):
            ....:     return min(len(X), 3)
            sage: M1 = RankMatroid(groundset=range(6), rank_function=f)
            sage: M2 = RankMatroid(groundset=range(6), rank_function=g)
            sage: M3 = RankMatroid(groundset=range(7), rank_function=f)
            sage: M4 = RankMatroid(groundset=range(6), rank_function=f)
            sage: M1 != M2  # indirect doctest
            True
            sage: M1 != M3
            True
            sage: M1 != M4
            False
        """
        return not self == other

    # Copying, loading, saving:

    def __reduce__(self):
        """
        Save the matroid for later reloading.

        .. NOTE::

            Unfortunately, functions cannot be pickled reliably, so this class
            doesn't have load/save support

        EXAMPLES::

            sage: M = Matroid(groundset=range(10),
            ....:             rank_function=lambda X: min(len(X), 4))
            sage: M == loads(dumps(M))  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: unfortunately, functions cannot be saved reliably, so
            this class doesn't have load/save support. Convert to another
            class, such as BasisMatroid, instead.
        """
        raise TypeError("unfortunately, functions cannot be saved reliably, so this class doesn't have load/save support. Convert to another class, such as BasisMatroid, instead.")
