# -*- coding: utf-8 -*-
r"""
Skew Tableaux

AUTHORS:

- Mike Hansen: Initial version
- Travis Scrimshaw, Arthur Lubovsky (2013-02-11):
  Factored out ``CombinatorialClass``
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#       Copyright (C) 2013 Arthur Lubovsky
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

from sage.rings.all import Integer, QQ, ZZ
from sage.arith.all import factorial
from sage.rings.infinity import PlusInfinity
from sage.matrix.all import zero_matrix

from sage.structure.list_clone import ClonableList
from sage.combinat.partition import Partition
from sage.combinat.tableau import (Tableau, Tableaux,
                                   StandardTableau, SemistandardTableau)
from sage.combinat.skew_partition import SkewPartition, SkewPartitions
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.words.words import Words

from sage.misc.persist import register_unpickle_override


class SkewTableau(ClonableList,
                  metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A skew tableau.

    Note that Sage by default uses the English convention for partitions and
    tableaux. To change this, see :meth:`Tableaux.options`.

    EXAMPLES::

         sage: st = SkewTableau([[None, 1],[2,3]]); st
         [[None, 1], [2, 3]]
         sage: st.inner_shape()
         [1]
         sage: st.outer_shape()
         [2, 2]

    The ``expr`` form of a skew tableau consists of the inner partition
    followed by a list of the entries in each row from bottom to top::

        sage: SkewTableau(expr=[[1,1],[[5],[3,4],[1,2]]])
        [[None, 1, 2], [None, 3, 4], [5]]

    The ``chain`` form of a skew tableau consists of a list of
    partitions `\lambda_1,\lambda_2,\ldots,`, such that all cells in
    `\lambda_{i+1}` that are not in `\lambda_i` have entry `i`::

        sage: SkewTableau(chain=[[2], [2, 1], [3, 1], [4, 3, 2, 1]])
        [[None, None, 2, 3], [1, 3, 3], [3, 3], [3]]
    """
    @staticmethod
    def __classcall_private__(cls, st=None, expr=None, chain=None):
        """
        Return the skew tableau object corresponding to ``st``.

        EXAMPLES::

            sage: SkewTableau([[None,1],[2,3]])
            [[None, 1], [2, 3]]
            sage: SkewTableau(expr=[[1,1],[[5],[3,4],[1,2]]])
            [[None, 1, 2], [None, 3, 4], [5]]
        """
        if isinstance(st, cls):
            return st
        if expr is not None:
            return SkewTableaux().from_expr(expr)
        if chain is not None:
            return SkewTableaux().from_chain(chain)

        return SkewTableaux()(st)

    def __init__(self, parent, st):
        """
        TESTS::

            sage: st = SkewTableau([[None, 1],[2,3]])
            sage: st = SkewTableau([[None,1,1],[None,2],[4]])
            sage: TestSuite(st).run()

        A skew tableau is immutable, see :trac:`15862`::

            sage: T = SkewTableau([[None,2],[2]])
            sage: t0 = T[0]
            sage: t0[1] = 3
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
            sage: T[0][1] = 5
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
        """
        try:
            st = [tuple(t) for t in st]
        except TypeError:
            raise TypeError("each element of the skew tableau must be an iterable")

        ClonableList.__init__(self, parent, st)

    def __eq__(self, other):
        r"""
        Check whether ``self`` is equal to ``other``.

        .. TODO::

            This overwrites the equality check of
            :class:`~sage.structure.list_clone.ClonableList`
            in order to circumvent the coercion framework.
            Eventually this should be solved more elegantly,
            for example along the lines of what was done for
            `k`-tableaux.

            For now, two elements are equal if their underlying
            defining lists compare equal.

        INPUT:

        ``other`` -- the element that ``self`` is compared to

        OUTPUT:

        A Boolean.

        TESTS::

            sage: t = SkewTableau([[None,1,2]])
            sage: t == 0
            False
            sage: t == SkewTableaux()([[None,1,2]])
            True

            sage: s = SkewTableau([[1,2]])
            sage: s == 0
            False
            sage: s == Tableau([[1,2]])
            True
        """
        if isinstance(other, (Tableau, SkewTableau)):
            return list(self) == list(other)
        else:
            return list(self) == other

    def __ne__(self, other):
        r"""
        Check whether ``self`` is unequal to ``other``.

        See the documentation of :meth:`__eq__`.

        INPUT:

        ``other`` -- the element that ``self`` is compared to

        OUTPUT:

        A Boolean.

        TESTS::

            sage: t = Tableau([[2,3],[1]])
            sage: t != []
            True
        """
        if isinstance(other, (Tableau, SkewTableau)):
            return list(self) != list(other)
        else:
            return list(self) != other

    def check(self):
        r"""
        Check that ``self`` is a valid skew tableau. This is currently far too
        liberal, and only checks some trivial things.

        EXAMPLES::

            sage: t = SkewTableau([[None,1,1],[2]])
            sage: t.check()

            sage: t = SkewTableau([[None, None, 1], [2, 4], [], [3, 4, 5]])
            Traceback (most recent call last):
            ...
            TypeError: a skew tableau cannot have an empty list for a row

            sage: s = SkewTableau([[1, None, None],[2, None],[3]])
            Traceback (most recent call last):
            ...
            TypeError: not a valid skew tableau
        """
        for row in self:
            if not row:
                raise TypeError("a skew tableau cannot have an empty list for a row")
            inside = False
            for x in row:
                if x is not None:
                    inside = True
                elif inside:
                    raise TypeError('not a valid skew tableau')

    def _repr_(self):
        """
        Return a string representation of ``self``.

        For more on the display options, see
        :obj:`SkewTableaux.options`.

        EXAMPLES::

            sage: SkewTableau([[None,2,3],[None,4],[5]])
            [[None, 2, 3], [None, 4], [5]]
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def _repr_list(self):
        """
        Return a string representation of ``self`` as a list of lists.

        EXAMPLES::

            sage: print(SkewTableau([[None,2,3],[None,4],[5]])._repr_list())
            [[None, 2, 3], [None, 4], [5]]
        """
        return repr(self.to_list())

    # See #18024. CombinatorialObject provided __str__, though ClonableList
    # doesn't. Emulate the old functionality. Possibly remove when
    # CombinatorialObject is removed.
    __str__ = _repr_list

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as a diagram.

        EXAMPLES::

            sage: print(SkewTableau([[None,2,3],[None,4],[5]])._repr_diagram())
              .  2  3
              .  4
              5
        """
        def none_str(x):
            return "  ." if x is None else "%3s" % str(x)
        if self.parent().options('convention') == "French":
            new_rows = ["".join(map(none_str, row)) for row in reversed(self)]
        else:
            new_rows = ["".join(map(none_str, row)) for row in self]
        return '\n'.join(new_rows)

    def _repr_compact(self):
        """
        Return a compact string representation of ``self``.

        EXAMPLES::

            sage: SkewTableau([[None,None,3],[4,5]])._repr_compact()
            '.,.,3/4,5'
            sage: Tableau([])._repr_compact()
            '-'
        """
        if not self:
            return '-'

        def str_rep(x):
            return '%s' % x if x is not None else '.'
        return '/'.join(','.join(str_rep(r) for r in row) for row in self)

    def pp(self):
        """
        Return a pretty print string of the tableau.

        EXAMPLES::

            sage: SkewTableau([[None,2,3],[None,4],[5]]).pp()
              .  2  3
              .  4
              5
        """
        print(self._repr_diagram())

    def _ascii_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        TESTS::

            sage: T1 = SkewTableau([[None,2,3],[None,4],[5]])
            sage: T2 = SkewTableau([[None,None,3],[4,5]])
            sage: ascii_art([T1, T2])
            [   .  2  3            ]
            [   .  4       .  .  3 ]
            [   5      ,   4  5    ]
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self._repr_diagram().splitlines())

    def _unicode_art_(self):
        """
        Return a unicode art representation of ``self``.

        TESTS::

            sage: T = SkewTableau([[None,None,1,1,2],[None,2,3],[None,4],[None,5],[6]])
            sage: unicode_art(T)
                    ┌───┬───┬───┐
                    │ 1 │ 1 │ 2 │
                ┌───┼───┼───┴───┘
                │ 2 │ 3 │
                ├───┼───┘
                │ 4 │
                ├───┤
                │ 5 │
            ┌───┼───┘
            │ 6 │
            └───┘
        """
        from sage.combinat.output import ascii_art_table
        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt(ascii_art_table(self, use_unicode=True).splitlines())

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: latex(SkewTableau([[None,2,3],[None,4],[5]]))
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{2-3}
            &\lr{2}&\lr{3}\\\cline{2-3}
            &\lr{4}\\\cline{1-2}
            \lr{5}\\\cline{1-1}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        return tex_from_array(self)

    def outer_shape(self):
        """
        Return the outer shape of ``self``.

        EXAMPLES::

            sage: SkewTableau([[None,1,2],[None,3],[4]]).outer_shape()
            [3, 2, 1]
        """
        return Partition([len(row) for row in self])

    def inner_shape(self):
        """
        Return the inner shape of ``self``.

        EXAMPLES::

            sage: SkewTableau([[None,1,2],[None,3],[4]]).inner_shape()
            [1, 1]
            sage: SkewTableau([[1,2],[3,4],[7]]).inner_shape()
            []
            sage: SkewTableau([[None,None,None,2,3],[None,1],[None],[2]]).inner_shape()
            [3, 1, 1]
        """
        return Partition([x for x in (row.count(None) for row in self) if x != 0])

    def shape(self):
        r"""
        Return the shape of ``self``.

        EXAMPLES::

            sage: SkewTableau([[None,1,2],[None,3],[4]]).shape()
            [3, 2, 1] / [1, 1]
        """
        return SkewPartition([self.outer_shape(), self.inner_shape()])

    def outer_size(self):
        """
        Return the size of the outer shape of ``self``.

        EXAMPLES::

            sage: SkewTableau([[None, 2, 4], [None, 3], [1]]).outer_size()
            6
            sage: SkewTableau([[None, 2], [1, 3]]).outer_size()
            4
        """
        return self.outer_shape().size()

    def inner_size(self):
        """
        Return the size of the inner shape of ``self``.

        EXAMPLES::

            sage: SkewTableau([[None, 2, 4], [None, 3], [1]]).inner_size()
            2
            sage: SkewTableau([[None, 2], [1, 3]]).inner_size()
            1
        """
        return self.inner_shape().size()

    def size(self):
        """
        Return the number of cells in ``self``.

        EXAMPLES::

            sage: SkewTableau([[None, 2, 4], [None, 3], [1]]).size()
            4
            sage: SkewTableau([[None, 2], [1, 3]]).size()
            3
        """
        one = ZZ.one()
        return sum(one for row in self for x in row if x is not None)

    def conjugate(self):
        """
        Return the conjugate of ``self``.

        EXAMPLES::

            sage: SkewTableau([[None,1],[2,3]]).conjugate()
            [[None, 2], [1, 3]]
        """
        conj_shape = self.outer_shape().conjugate()

        conj = [[None] * row_length for row_length in conj_shape]

        for i in range(len(conj)):
            for j in range(len(conj[i])):
                conj[i][j] = self[j][i]

        return SkewTableau(conj)

    def to_word_by_row(self):
        """
        Return a word obtained from a row reading of ``self``.

        This is the word obtained by concatenating the rows from
        the bottommost one (in English notation) to the topmost one.

        EXAMPLES::

            sage: s = SkewTableau([[None,1],[2,3]])
            sage: s.pp()
              .  1
              2  3
            sage: s.to_word_by_row()
            word: 231
            sage: s = SkewTableau([[None, 2, 4], [None, 3], [1]])
            sage: s.pp()
              .  2  4
              .  3
              1
            sage: s.to_word_by_row()
            word: 1324

        TESTS::

            sage: SkewTableau([[None, None, None], [None]]).to_word_by_row()
            word:
            sage: SkewTableau([]).to_word_by_row()
            word:
        """
        word = [x for row in reversed(self) for x in row if x is not None]
        return Words("positive integers")(word)

    def to_word_by_column(self):
        """
        Return the word obtained from a column reading of the skew
        tableau.

        This is the word obtained by concatenating the columns from
        the rightmost one (in English notation) to the leftmost one.

        EXAMPLES::

            sage: s = SkewTableau([[None,1],[2,3]])
            sage: s.pp()
              .  1
              2  3
            sage: s.to_word_by_column()
            word: 132

        ::

            sage: s = SkewTableau([[None, 2, 4], [None, 3], [1]])
            sage: s.pp()
            .  2  4
            .  3
            1
            sage: s.to_word_by_column()
            word: 4231
        """
        return self.conjugate().to_word_by_row()

    to_word = to_word_by_row

    def to_permutation(self):
        """
        Return a permutation with the entries of ``self`` obtained by reading
        ``self`` row by row, from the bottommost to the topmost row, with
        each row being read from left to right, in English convention.
        See :meth:`to_word_by_row()`.

        EXAMPLES::

            sage: SkewTableau([[None,2],[3,4],[None],[1]]).to_permutation()
            [1, 3, 4, 2]
            sage: SkewTableau([[None,2],[None,4],[1],[3]]).to_permutation()
            [3, 1, 4, 2]
            sage: SkewTableau([[None]]).to_permutation()
            []
        """
        from sage.combinat.permutation import Permutation
        perm = [i for row in reversed(self) for i in row if i is not None]
        return Permutation(perm)

    def weight(self):
        r"""
        Return the weight (aka evaluation) of the tableau ``self``.
        Trailing zeroes are omitted when returning the weight.

        The weight of a skew tableau `T` is the sequence
        `(a_1, a_2, a_3, \ldots )`, where `a_k` is the number of
        entries of `T` equal to `k`. This sequence contains only
        finitely many nonzero entries.

        The weight of a skew tableau `T` is the same as the weight
        of the reading word of `T`, for any reading order.

        :meth:`evaluation` is a synonym for this method.

        EXAMPLES::

            sage: SkewTableau([[1,2],[3,4]]).weight()
            [1, 1, 1, 1]

            sage: SkewTableau([[None,2],[None,4],[None,5],[None]]).weight()
            [0, 1, 0, 1, 1]

            sage: SkewTableau([]).weight()
            []

            sage: SkewTableau([[None,None,None],[None]]).weight()
            []

            sage: SkewTableau([[None,3,4],[None,6,7],[4,8],[5,13],[6],[7]]).weight()
            [0, 0, 1, 2, 1, 2, 2, 1, 0, 0, 0, 0, 1]

        TESTS:

        We check that this agrees with going to the word::

            sage: t = SkewTableau([[None,None,4,7,15],[6,2,16],[2,3,19],[4,5],[7]])
            sage: def by_word(T):
            ....:     ed = T.to_word().evaluation_dict()
            ....:     m = max(ed) + 1
            ....:     return [ed.get(k, 0) for k in range(1, m)]
            sage: by_word(t) == t.weight()
            True
            sage: SST = SemistandardTableaux(shape=[3,1,1])
            sage: all(by_word(t) == SkewTableau(t).weight() for t in SST)
            True
        """
        if (not self) or all(c is None for row in self for c in row):
            return []
        m = max(c for row in self for c in row if c is not None)
        if m is None:
            return []
        res = [0] * m
        for row in self:
            for i in row:
                if (i is not None) and i > 0:
                    res[i - 1] += 1
        return res

    evaluation = weight

    def is_standard(self):
        """
        Return ``True`` if ``self`` is a standard skew tableau and ``False``
        otherwise.

        EXAMPLES::

            sage: SkewTableau([[None, 2], [1, 3]]).is_standard()
            True
            sage: SkewTableau([[None, 2], [2, 4]]).is_standard()
            False
            sage: SkewTableau([[None, 3], [2, 4]]).is_standard()
            False
            sage: SkewTableau([[None, 2], [2, 4]]).is_standard()
            False
        """
        # Check to make sure that it is filled with 1...size
        w = [i for row in self for i in row if i is not None]
        return sorted(w) == list(range(1, len(w) + 1)) and self.is_semistandard()

    def is_semistandard(self):
        """
        Return ``True`` if ``self`` is a semistandard skew tableau and
        ``False`` otherwise.

        EXAMPLES::

            sage: SkewTableau([[None, 2, 2], [1, 3]]).is_semistandard()
            True
            sage: SkewTableau([[None, 2], [2, 4]]).is_semistandard()
            True
            sage: SkewTableau([[None, 3], [2, 4]]).is_semistandard()
            True
            sage: SkewTableau([[None, 2], [1, 2]]).is_semistandard()
            False
            sage: SkewTableau([[None, 2, 3]]).is_semistandard()
            True
            sage: SkewTableau([[None, 3, 2]]).is_semistandard()
            False
            sage: SkewTableau([[None, 2, 3], [1, 4]]).is_semistandard()
            True
            sage: SkewTableau([[None, 2, 3], [1, 2]]).is_semistandard()
            False
            sage: SkewTableau([[None, 2, 3], [None, None, 4]]).is_semistandard()
            False
        """
        if not self:
            return True

        # Is it weakly increasing along the rows?
        for row in self:
            if any(row[c] is not None and row[c] > row[c + 1]
                   for c in range(len(row) - 1)):
                return False

        # Is it strictly increasing down columns?
        for row, next in zip(self, self[1:]):
            if any(row[c] is not None and row[c] >= next[c] for c in range(len(next))):
                return False

        return True

    def to_tableau(self):
        """
        Returns a tableau with the same filling. This only works if the
        inner shape of the skew tableau has size zero.

        EXAMPLES::

            sage: SkewTableau([[1,2],[3,4]]).to_tableau()
            [[1, 2], [3, 4]]
        """

        if self.inner_size() != 0:
            raise ValueError("the inner size of the skew tableau must be 0")
        from sage.combinat.tableau import Tableau
        return Tableau(self[:])

    def restrict(self, n):
        """
        Return the restriction of the (semi)standard skew tableau to all
        the numbers less than or equal to ``n``.

        .. NOTE::

            If only the outer shape of the restriction, rather than
            the whole restriction, is needed, then the faster method
            :meth:`restriction_outer_shape` is preferred. Similarly if
            only the skew shape is needed, use :meth:`restriction_shape`.

        EXAMPLES::

            sage: SkewTableau([[None,1],[2],[3]]).restrict(2)
            [[None, 1], [2]]
            sage: SkewTableau([[None,1],[2],[3]]).restrict(1)
            [[None, 1]]
            sage: SkewTableau([[None,1],[1],[2]]).restrict(1)
            [[None, 1], [1]]
        """
        data = ([y for y in x if y is None or y <= n] for x in self)
        return SkewTableau([z for z in data if z])

    def restriction_outer_shape(self, n):
        """
        Return the outer shape of the restriction of the semistandard skew
        tableau ``self`` to `n`.

        If `T` is a semistandard skew tableau and `n` is a nonnegative
        integer, then the restriction of `T` to `n` is defined as the
        (semistandard) skew tableau obtained by removing all cells filled
        with entries greater than `n` from `T`.

        This method computes merely the outer shape of the restriction.
        For the restriction itself, use :meth:`restrict`.

        EXAMPLES::

            sage: SkewTableau([[None,None],[2,3],[3,4]]).restriction_outer_shape(3)
            [2, 2, 1]
            sage: SkewTableau([[None,2],[None],[4],[5]]).restriction_outer_shape(2)
            [2, 1]
            sage: T = SkewTableau([[None,None,3,5],[None,4,4],[17]])
            sage: T.restriction_outer_shape(0)
            [2, 1]
            sage: T.restriction_outer_shape(2)
            [2, 1]
            sage: T.restriction_outer_shape(3)
            [3, 1]
            sage: T.restriction_outer_shape(4)
            [3, 3]
            sage: T.restriction_outer_shape(19)
            [4, 3, 1]
        """
        from sage.combinat.partition import _Partitions
        one = ZZ.one()
        res = [sum(one for y in row if y is None or y <= n) for row in self]
        return _Partitions(res)

    def restriction_shape(self, n):
        """
        Return the skew shape of the restriction of the semistandard
        skew tableau ``self`` to ``n``.

        If `T` is a semistandard skew tableau and `n` is a nonnegative
        integer, then the restriction of `T` to `n` is defined as the
        (semistandard) skew tableau obtained by removing all cells filled
        with entries greater than `n` from `T`.

        This method computes merely the skew shape of the restriction.
        For the restriction itself, use :meth:`restrict`.

        EXAMPLES::

            sage: SkewTableau([[None,None],[2,3],[3,4]]).restriction_shape(3)
            [2, 2, 1] / [2]
            sage: SkewTableau([[None,2],[None],[4],[5]]).restriction_shape(2)
            [2, 1] / [1, 1]
            sage: T = SkewTableau([[None,None,3,5],[None,4,4],[17]])
            sage: T.restriction_shape(0)
            [2, 1] / [2, 1]
            sage: T.restriction_shape(2)
            [2, 1] / [2, 1]
            sage: T.restriction_shape(3)
            [3, 1] / [2, 1]
            sage: T.restriction_shape(4)
            [3, 3] / [2, 1]
        """
        return SkewPartition([self.restriction_outer_shape(n), self.inner_shape()])

    def to_chain(self, max_entry=None):
        r"""
        Return the chain of partitions corresponding to the (semi)standard
        skew tableau ``self``.

        The optional keyword parameter ``max_entry`` can be used to
        customize the length of the chain. Specifically, if this parameter
        is set to a nonnegative integer ``n``, then the chain is
        constructed from the positions of the letters `1, 2, \ldots, n`
        in the tableau.

        EXAMPLES::

            sage: SkewTableau([[None,1],[2],[3]]).to_chain()
            [[1], [2], [2, 1], [2, 1, 1]]
            sage: SkewTableau([[None,1],[1],[2]]).to_chain()
            [[1], [2, 1], [2, 1, 1]]
            sage: SkewTableau([[None,1],[1],[2]]).to_chain(max_entry=2)
            [[1], [2, 1], [2, 1, 1]]
            sage: SkewTableau([[None,1],[1],[2]]).to_chain(max_entry=3)
            [[1], [2, 1], [2, 1, 1], [2, 1, 1]]
            sage: SkewTableau([[None,1],[1],[2]]).to_chain(max_entry=1)
            [[1], [2, 1]]
            sage: SkewTableau([[None,None,2],[None,3],[None,5]]).to_chain(max_entry=6)
            [[2, 1, 1], [2, 1, 1], [3, 1, 1], [3, 2, 1], [3, 2, 1], [3, 2, 2], [3, 2, 2]]
            sage: SkewTableau([]).to_chain()
            [[]]
            sage: SkewTableau([]).to_chain(max_entry=1)
            [[], []]

        TESTS:

        Check that :meth:`to_chain()` does not skip letters::

            sage: t = SkewTableau([[None, 2, 3], [3]])
            sage: t.to_chain()
            [[1], [1], [2], [3, 1]]

            sage: T = SkewTableau([[None]])
            sage: T.to_chain()
            [[1]]
        """
        if max_entry is None:
            if (not self) or all(c is None for row in self for c in row):
                max_entry = 0
            else:
                max_entry = max(c for row in self for c in row if c is not None)
        return [self.restriction_outer_shape(x) for x in range(max_entry + 1)]

    def slide(self, corner=None, return_vacated=False):
        """
        Apply a jeu-de-taquin slide to ``self`` on the specified inner
        corner and return the resulting tableau.

        If no corner is given, the topmost inner corner is chosen.

        The optional parameter ``return_vacated=True`` causes
        the output to be the pair ``(t, (i, j))`` where ``t`` is the new
        tableau and ``(i, j)`` are the coordinates of the vacated square.

        See [Ful1997]_ p12-13.

        EXAMPLES::

            sage: st = SkewTableau([[None, None, None, None, 2], [None, None, None, None, 6], [None, 2, 4, 4], [2, 3, 6], [5, 5]])
            sage: st.slide((2, 0))
            [[None, None, None, None, 2], [None, None, None, None, 6], [2, 2, 4, 4], [3, 5, 6], [5]]
            sage: st2 = SkewTableau([[None, None, 3], [None, 2, 4], [1, 5]])
            sage: st2.slide((1, 0), True)
            ([[None, None, 3], [1, 2, 4], [5]], (2, 1))

        TESTS::

            sage: st
            [[None, None, None, None, 2], [None, None, None, None, 6],
             [None, 2, 4, 4], [2, 3, 6], [5, 5]]
        """
        new_st = self.to_list()
        inner_corners = self.inner_shape().corners()
        outer_corners = self.outer_shape().corners()
        if corner is not None:
            if tuple(corner) not in inner_corners:
                raise ValueError("corner must be an inner corner")
        else:
            if not inner_corners:
                return self
            else:
                corner = inner_corners[0]

        spotl, spotc = corner
        while (spotl, spotc) not in outer_corners:
            # Check to see if there is nothing to the right
            if spotc == len(new_st[spotl]) - 1:
                # Swap the hole with the cell below
                new_st[spotl][spotc] = new_st[spotl + 1][spotc]
                new_st[spotl + 1][spotc] = None
                spotl += 1
                continue

            # Check to see if there is nothing below
            if spotl == len(new_st) - 1 or len(new_st[spotl + 1]) <= spotc:
                # Swap the hole with the cell to the right
                new_st[spotl][spotc] = new_st[spotl][spotc + 1]
                new_st[spotl][spotc + 1] = None
                spotc += 1
                continue

            # If we get to this stage, we need to compare
            below = new_st[spotl + 1][spotc]
            right = new_st[spotl][spotc + 1]
            if below <= right:
                # Swap with the cell below
                new_st[spotl][spotc] = new_st[spotl + 1][spotc]
                new_st[spotl + 1][spotc] = None
                spotl += 1
                continue

            # Otherwise swap with the cell to the right
            new_st[spotl][spotc] = new_st[spotl][spotc + 1]
            new_st[spotl][spotc + 1] = None
            spotc += 1

        # Clean up to remove the "None" at an outside corner
        # Remove the last row if there is nothing left in it
        new_st[spotl].pop()
        if not new_st[spotl]:
            new_st.pop()

        if return_vacated:
            return (SkewTableau(new_st), (spotl, spotc))
        return SkewTableau(new_st)

    def rectify(self, algorithm=None):
        """
        Return a :class:`StandardTableau`, :class:`SemistandardTableau`,
        or just :class:`Tableau` formed by applying the jeu de taquin
        process to ``self``.

        See page 15 of [Ful1997]_.

        INPUT:

        - ``algorithm`` -- optional: if set to ``'jdt'``, rectifies by jeu de
          taquin; if set to ``'schensted'``, rectifies by Schensted insertion
          of the reading word; otherwise, guesses which will be faster.

        EXAMPLES::

            sage: S = SkewTableau([[None,1],[2,3]])
            sage: S.rectify()
            [[1, 3], [2]]
            sage: T = SkewTableau([[None, None, None, 4],[None,None,1,6],[None,None,5],[2,3]])
            sage: T.rectify()
            [[1, 3, 4, 6], [2, 5]]
            sage: T.rectify(algorithm='jdt')
            [[1, 3, 4, 6], [2, 5]]
            sage: T.rectify(algorithm='schensted')
            [[1, 3, 4, 6], [2, 5]]
            sage: T.rectify(algorithm='spaghetti')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'jdt', 'schensted', or None

        TESTS::

            sage: S
            [[None, 1], [2, 3]]
            sage: T
            [[None, None, None, 4], [None, None, 1, 6], [None, None, 5], [2, 3]]
        """
        mu_size = self.inner_shape().size()

        # Roughly, use jdt with a small inner shape, Schensted with a large one
        if algorithm is None:
            la = self.outer_shape()
            la_size = la.size()
            if mu_size ** 2 < len(la) * (la_size - mu_size):
                algorithm = 'jdt'
            else:
                algorithm = 'schensted'

        if algorithm == 'jdt':
            rect = self
            for i in range(mu_size):
                rect = rect.slide()
        elif algorithm == 'schensted':
            w = [x for row in reversed(self) for x in row if x is not None]
            rect = Tableau([]).insert_word(w)
        else:
            raise ValueError("algorithm must be 'jdt', 'schensted', or None")
        if self in StandardSkewTableaux():
            return StandardTableau(rect[:])
        if self in SemistandardSkewTableaux():
            return SemistandardTableau(rect[:])
        return Tableau(rect)

    def to_list(self):
        r"""
        Return a (mutable) list representation of ``self``.

        EXAMPLES::

            sage: stlist = [[None, None, 3], [None, 1, 3], [2, 2]]
            sage: st = SkewTableau(stlist)
            sage: st.to_list()
            [[None, None, 3], [None, 1, 3], [2, 2]]
            sage: st.to_list() == stlist
            True

        """
        return [list(row) for row in self]

    def shuffle(self, t2):
        r"""
        Shuffle the standard tableaux ``self`` and ``t2``.

        Let ``t1 = self``. The shape of ``t2`` must extend the shape of
        ``t1``, that is, ``self.outer_shape() == t2.inner_shape()``. Then
        this function computes the pair of tableaux ``(t2_new, t1_new)``
        obtained by using jeu de taquin slides to move the boxes of ``t2``
        behind the boxes of ``self``.

        The entries of ``t2_new`` are obtained by performing successive
        inwards jeu de taquin slides on ``t2`` in the order indicated by
        the entries of ``t1``, from largest to smallest. The entries of
        ``t1`` then slide outwards one by one and land in the squares
        vacated successively by ``t2``, forming ``t1_new``.

        .. NOTE::

            Equivalently, the entries of ``t1_new`` are obtained by performing
            outer jeu de taquin slides on ``t1`` in the order indicated by the
            entries of ``t2``, from smallest to largest. In this case the
            entries of ``t2`` slide backwards and fill the squares
            successively vacated by ``t1`` and so form ``t2_new``.
            (This is not how the algorithm is implemented.)

        INPUT:

        - ``self``, ``t2`` -- a pair of standard SkewTableaux with
          ``self.outer_shape() == t2.inner_shape()``

        OUTPUT:

        - ``t2_new, t1_new`` -- a pair of standard :class:`SkewTableaux`
          with ``t2_new.outer_shape() == t1_new.inner_shape()``

        EXAMPLES::

            sage: t1 = SkewTableau([[None, 1, 2], [3, 4]])
            sage: t2 = SkewTableau([[None, None, None, 3], [None, None, 4], [1, 2, 5]])
            sage: (t2_new, t1_new) = t1.shuffle(t2)
            sage: t1_new
            [[None, None, None, 2], [None, None, 1], [None, 3, 4]]
            sage: t2_new
            [[None, 2, 3], [1, 4], [5]]
            sage: t1_new.outer_shape() == t2.outer_shape()
            True
            sage: t2_new.inner_shape() == t1.inner_shape()
            True

        Shuffling is an involution::

            sage: t1 = SkewTableau([[None, 1, 2], [3, 4]])
            sage: t2 = SkewTableau([[None, None, None, 3], [None, None, 4], [1, 2, 5]])
            sage: sh = lambda x,y : x.shuffle(y)
            sage: (t1, t2) == sh(*sh(t1, t2))
            True

        Both tableaux must be standard::

            sage: t1 = SkewTableau([[None, 1, 2], [2, 4]])
            sage: t2 = SkewTableau([[None, None, None, 3], [None, None, 4], [1, 2, 5]])
            sage: t1.shuffle(t2)
            Traceback (most recent call last):
            ...
            ValueError: the tableaux must be standard
            sage: t1 = SkewTableau([[None, 1, 2], [3, 4]])
            sage: t2 = SkewTableau([[None, None, None, 3], [None, None, 4], [1, 2, 6]])
            sage: t1.shuffle(t2)
            Traceback (most recent call last):
            ...
            ValueError: the tableaux must be standard

        The shapes (not just the nonempty cells) must be adjacent::

            sage: t1 = SkewTableau([[None, None, None], [1]])
            sage: t2 = SkewTableau([[None], [None], [1]])
            sage: t1.shuffle(t2)
            Traceback (most recent call last):
            ...
            ValueError: the shapes must be adjacent

        TESTS:

        A corner case, where one tableau has no cells::

            sage: t1 = SkewTableau([[None]])
            sage: t2 = SkewTableau([[None, 1, 2], [3, 4]])
            sage: (t2_new, t1_new) = t1.shuffle(t2)
            sage: t1_new
            [[None, None, None], [None, None]]
            sage: t2_new == t2
            True
            sage: t2_new.shuffle(t1_new) == (t1, t2)
            True
        """
        if self.outer_shape() != t2.inner_shape():
            raise ValueError("the shapes must be adjacent")
        if not self.is_standard() or not t2.is_standard():
            raise ValueError("the tableaux must be standard")

        # start with t2_new = t2, which we will slide backwards
        t2_new = t2

        # make a blank copy of t2 (to fill in iteratively), which will become t1_new
        t1_new = [[None] * len(x) for x in list(t2)]

        # perform reverse slides according to the entries of t1,
        # from largest to smallest
        m = self.size()
        for i in range(m, 0, -1):
            # the square of t1 containing i
            corner = self.cells_containing(i)[0]

            # slide t2_new backwards, record i in the vacated square
            (t2_new, (x, y)) = t2_new.slide(corner, True)
            t1_new[x][y] = i

        t1_new = SkewTableau(t1_new)
        return t2_new, t1_new

    def standardization(self, check=True):
        r"""
        Return the standardization of ``self``, assuming ``self`` is a
        semistandard skew tableau.

        The standardization of a semistandard skew tableau `T` is the standard
        skew tableau `\mathrm{st}(T)` of the same shape as `T` whose
        reversed reading word is the standardization of the reversed reading
        word of `T`.

        The standardization of a word `w` can be formed by replacing all `1`'s
        in `w` by `1, 2, \ldots, k_1` from left to right, all `2`'s in `w` by
        `k_1 + 1, k_1 + 2, \ldots, k_2`, and repeating for all letters that
        appear in `w`.
        See also :meth:`Word.standard_permutation()`.

        INPUT:

        - ``check`` -- (Default: ``True``) Check to make sure ``self`` is
          semistandard. Set to ``False`` to avoid this check.

        EXAMPLES::

            sage: t = SkewTableau([[None,None,3,4,7,19],[None,4,4,8],[None,5,16,17],[None],[2],[3]])
            sage: t.standardization()
            [[None, None, 3, 6, 8, 12], [None, 4, 5, 9], [None, 7, 10, 11], [None], [1], [2]]

        Standard skew tableaux are fixed under standardization::

            sage: p = Partition([4,3,3,2])
            sage: q = Partitions(3).random_element()
            sage: all((t == t.standardization() for t in StandardSkewTableaux([p, q])))
            True

        The reading word of the standardization is the
        standardization of the reading word::

            sage: t = SkewTableau([[None,3,4,4],[None,6,10],[7,7,11],[18]])
            sage: t.to_word().standard_permutation() == t.standardization().to_permutation()
            True

        TESTS:

        Some corner cases::

            sage: t = SkewTableau([[None,None],[None]])
            sage: t.standardization()
            [[None, None], [None]]
            sage: t = SkewTableau([])
            sage: t.standardization()
            []
        """
        if check and not self.is_semistandard():
            raise ValueError("the skew tableau must be semistandard")
        # This should be a SkewStandardTableau
        return StandardSkewTableaux().from_shape_and_word(self.shape(), self.to_word_by_row().standard_permutation())

    def bender_knuth_involution(self, k, rows=None, check=True):
        r"""
        Return the image of ``self`` under the `k`-th Bender--Knuth
        involution, assuming ``self`` is a skew semistandard tableau.

        Let `T` be a tableau, then a *lower free `k` in `T`* means a cell of
        `T` which is filled with the integer `k` and whose direct lower
        neighbor is not filled with the integer `k + 1` (in particular,
        this lower neighbor might not exist at all). Let an *upper free `k + 1`
        in `T`* mean a cell of `T` which is filled with the integer `k + 1`
        and whose direct upper neighbor is not filled with the integer `k`
        (in particular, this neighbor might not exist at all). It is clear
        that for any row `r` of `T`, the lower free `k`'s and the upper
        free `k + 1`'s in `r` together form a contiguous interval or `r`.

        The *`k`-th Bender--Knuth switch at row `i`* changes the entries of
        the cells in this interval in such a way that if it used to have
        `a` entries of `k` and `b` entries of `k + 1`, it will now
        have `b` entries of `k` and `a` entries of `k + 1`. For fixed `k`, the
        `k`-th Bender--Knuth switches for different `i` commute. The
        composition of the `k`-th Bender--Knuth switches for all rows is
        called the *`k`-th Bender--Knuth involution*. This is used to show that
        the Schur functions defined by semistandard (skew) tableaux are
        symmetric functions.

        INPUT:

        - ``k`` -- an integer

        - ``rows`` -- (Default ``None``) When set to ``None``, the method
          computes the `k`-th Bender--Knuth involution as defined above.
          When an iterable, this computes the composition of the `k`-th
          Bender--Knuth switches at row `i` over all `i` in ``rows``. When set
          to an integer `i`, the method computes the `k`-th Bender--Knuth
          switch at row `i`. Note the indexing of the rows starts with `1`.

        - ``check`` -- (Default: ``True``) Check to make sure ``self`` is
          semistandard. Set to ``False`` to avoid this check.

        OUTPUT:

        The image of ``self`` under either the `k`-th Bender--Knuth
        involution, the `k`-th Bender--Knuth switch at a certain row, or
        the composition of such switches, as detailed in the INPUT section.

        EXAMPLES::

            sage: t = SkewTableau([[None,None,None,4,4,5,6,7],[None,2,4,6,7,7,7],
            ....:                  [None,4,5,8,8,9],[None,6,7,10],[None,8,8,11],[None],[4]])
            sage: t
            [[None, None, None, 4, 4, 5, 6, 7], [None, 2, 4, 6, 7, 7, 7],
             [None, 4, 5, 8, 8, 9], [None, 6, 7, 10], [None, 8, 8, 11], [None], [4]]
            sage: t.bender_knuth_involution(1)
            [[None, None, None, 4, 4, 5, 6, 7], [None, 1, 4, 6, 7, 7, 7],
             [None, 4, 5, 8, 8, 9], [None, 6, 7, 10], [None, 8, 8, 11], [None], [4]]
            sage: t.bender_knuth_involution(4)
            [[None, None, None, 4, 5, 5, 6, 7], [None, 2, 4, 6, 7, 7, 7],
             [None, 5, 5, 8, 8, 9], [None, 6, 7, 10], [None, 8, 8, 11], [None], [5]]
            sage: t.bender_knuth_involution(5)
            [[None, None, None, 4, 4, 5, 6, 7], [None, 2, 4, 5, 7, 7, 7],
             [None, 4, 6, 8, 8, 9], [None, 5, 7, 10], [None, 8, 8, 11], [None], [4]]
            sage: t.bender_knuth_involution(6)
            [[None, None, None, 4, 4, 5, 6, 6], [None, 2, 4, 6, 6, 7, 7],
             [None, 4, 5, 8, 8, 9], [None, 6, 7, 10], [None, 8, 8, 11], [None], [4]]
            sage: t.bender_knuth_involution(666) == t
            True
            sage: t.bender_knuth_involution(4, 2) == t
            True
            sage: t.bender_knuth_involution(4, 3)
            [[None, None, None, 4, 4, 5, 6, 7], [None, 2, 4, 6, 7, 7, 7],
             [None, 5, 5, 8, 8, 9], [None, 6, 7, 10], [None, 8, 8, 11], [None], [4]]

        The Bender--Knuth involution is an involution::

            sage: t = SkewTableau([[None,3,4,4],[None,6,10],[7,7,11],[18]])
            sage: all(t.bender_knuth_involution(k).bender_knuth_involution(k)
            ....:     == t for k in range(1,4))
            True

        The same for the single switches::

            sage: all(t.bender_knuth_involution(k, j).bender_knuth_involution(k, j)
            ....:     == t for k in range(1,5) for j in range(1, 5))
            True

        Locality of the Bender--Knuth involutions::

            sage: all(t.bender_knuth_involution(k).bender_knuth_involution(l)
            ....:     == t.bender_knuth_involution(l).bender_knuth_involution(k)
            ....:     for k in range(1,5) for l in range(1,5) if abs(k - l) > 1)
            True

        TESTS::

            sage: t = SkewTableau([])
            sage: t.bender_knuth_involution(3)
            []
            sage: t = SkewTableau([[None,None],[None]])
            sage: t.bender_knuth_involution(3)
            [[None, None], [None]]

        The `(s_1 s_2)^6 = id` identity that holds for Bender--Knuth
        involutions on straight shapes does not generally hold for
        skew shapes::

            sage: p = lambda t, k: t.bender_knuth_involution(k).bender_knuth_involution(k + 1)
            sage: t = SkewTableau([[None,1,2],[2,3]])
            sage: x = t
            sage: for i in range(6): x = p(x, 1)
            sage: x
            [[None, 2, 2], [1, 3]]
            sage: x == t
            False

        AUTHORS:

        - Darij Grinberg (2013-05-14)
        """
        if check and not self.is_semistandard():
            raise ValueError("the skew tableau must be semistandard")
        l = len(self)    # l is the number of rows of self.
        # Sanitizing the rows input so that it always becomes a list of
        # nonnegative integers. We also subtract 1 from these integers
        # because the i-th row of a tableau T is T[i - 1].
        if rows is None:
            rows = list(range(l))
        elif rows in ZZ:
            rows = [rows - 1]
        else:
            rows = [i - 1 for i in rows]
        # Now, rows should be iterable.

        # result_tab is going to be the result tableau (as a list of lists);
        # we will build it up step by step, starting with a deep copy of self.
        result_tab = self.to_list()
        for i in rows:
            if i >= l:
                continue
            # Setup the previous and next rows
            if i == 0:
                prev_row = [None] * len(result_tab[i])
            else:
                prev_row = result_tab[i - 1]
            if i == l - 1:
                next_row = [None] * len(result_tab[i])
            else:
                next_row = result_tab[i + 1] + [None] * (len(result_tab[i]) - len(result_tab[i + 1]))
            a = 0
            b = 0
            sk = None  # The first entry of k
            sk1 = None  # The first entry of k+1
            for j, val in enumerate(result_tab[i]):
                if val == k and next_row[j] != k + 1:
                    if sk is None:
                        sk = j
                    a += 1
                elif val == k + 1 and prev_row[j] != k:
                    if sk1 is None:
                        sk1 = j
                    b += 1
            if sk1 is not None:
                if a > b:
                    for j in range(sk1 - (a - b), sk1):
                        result_tab[i][j] = k + 1
                elif a < b:
                    for j in range(sk1, sk1 + b - a):
                        result_tab[i][j] = k
            elif sk is not None:
                for j in range(sk, sk + a):
                    result_tab[i][j] = k + 1

        return SkewTableau(result_tab)  # This should be a SkewSemistandardTableau

    def to_expr(self):
        """
        The first list in a result corresponds to the inner partition of
        the skew shape. The second list is a list of the rows in the skew
        tableau read from the bottom up.

        Provided for compatibility with MuPAD-Combinat. In MuPAD-Combinat,
        if ``t`` is a skew tableau, then to_expr gives the same result as
        ``expr(t)`` would give in MuPAD-Combinat.

        EXAMPLES::

            sage: SkewTableau([[None,1,1,3],[None,2,2],[1]]).to_expr()
            [[1, 1], [[1], [2, 2], [1, 1, 3]]]
            sage: SkewTableau([]).to_expr()
            [[], []]
        """
        rows = self.filling()
        rows.reverse()
        return [self.inner_shape(), rows]

    def is_ribbon(self):
        r"""
        Return ``True`` if and only if the shape of ``self`` is a
        ribbon, that is, if it has exactly one cell in each of `q`
        consecutive diagonals for some nonnegative integer `q`.

        EXAMPLES::

            sage: S=SkewTableau([[None, None, 1, 2],[None, None, 3],[1, 3, 4]])
            sage: S.pp()
              .  .  1  2
              .  .  3
              1  3  4
            sage: S.is_ribbon()
            True

            sage: S=SkewTableau([[None, 1, 1, 2],[None, 2, 3],[1, 3, 4]])
            sage: S.pp()
              .  1  1  2
              .  2  3
              1  3  4
            sage: S.is_ribbon()
            False

            sage: S=SkewTableau([[None, None, 1, 2],[None, None, 3],[1]])
            sage: S.pp()
              .  .  1  2
              .  .  3
              1
            sage: S.is_ribbon()
            False

            sage: S=SkewTableau([[None, None, None, None],[None, None, 3],[1, 2, 4]])
            sage: S.pp()
              .  .  .  .
              .  .  3
              1  2  4
            sage: S.is_ribbon()
            True

            sage: S=SkewTableau([[None, None, None, None],[None, None, 3],[None, 2, 4]])
            sage: S.pp()
              .  .  .  .
              .  .  3
              .  2  4
            sage: S.is_ribbon()
            True

            sage: S=SkewTableau([[None, None],[None]])
            sage: S.pp()
              .  .
              .
            sage: S.is_ribbon()
            True

        """
        lam = list(self.outer_shape())
        mu = list(self.inner_shape())
        l_out = len(lam)
        l_in = len(mu)
        mu += [0] * (l_out - l_in)

        if l_out == 0:
            return True
        else:
            # Find the least u for which lam[u]>mu[u], if it exists.
            # If it does not exist then u will equal l_out.
            u = 0
            u_test = True
            while u_test:
                if u >= l_out or lam[u] > mu[u]:
                    u_test = False
                else:
                    u += 1

            # Find the least v strictly greater than u for which
            # lam[v] != mu[v-1]+1
            v = u + 1
            v_test = True
            while v_test:
                if v >= l_out or lam[v] != mu[v - 1] + 1:
                    v_test = False
                else:
                    v += 1

            # Check if lam[i]==mu[i] for all i >= v
            for i in range(v, l_out):
                if lam[i] != mu[i]:
                    return False

        return True

    def to_ribbon(self, check_input=True):
        """
        Return ``self`` as a ribbon-shaped tableau
        (:class:`~sage.combinat.ribbon_shaped_tableau.RibbonShapedTableau`),
        provided that the shape of ``self`` is a ribbon.

        INPUT:

        - ``check_input`` -- (default: ``True``) whether or not to check
          that ``self`` indeed has ribbon shape

        EXAMPLES::

            sage: SkewTableau([[None,1],[2,3]]).to_ribbon()
            [[None, 1], [2, 3]]
        """
        if check_input and not self.is_ribbon():
            raise ValueError("self must be a ribbon")
        from sage.combinat.ribbon_shaped_tableau import RibbonShapedTableau
        r = [[i for i in row if i is not None] for row in self]
        return RibbonShapedTableau(r)

    def filling(self):
        """
        Return a list of the non-empty entries in ``self``.

        EXAMPLES::

            sage: t = SkewTableau([[None,1],[2,3]])
            sage: t.filling()
            [[1], [2, 3]]
        """
        return [[i for i in row if i is not None] for row in self]

    def cells_by_content(self, c):
        """
        Return the coordinates of the cells in ``self`` with content ``c``.

        EXAMPLES::

            sage: s = SkewTableau([[None,1,2],[3,4,5],[6]])
            sage: s.cells_by_content(0)
            [(1, 1)]
            sage: s.cells_by_content(1)
            [(0, 1), (1, 2)]
            sage: s.cells_by_content(2)
            [(0, 2)]
            sage: s.cells_by_content(-1)
            [(1, 0)]
            sage: s.cells_by_content(-2)
            [(2, 0)]
        """
        if not self:
            return []

        if c >= 0:
            if c >= len(self[0]):
                return []
            i, j = 0, c
        else:
            c = -c
            if c >= len(self):
                return []
            i, j = c, 0

        res = []
        while True:
            if self[i][j] is not None:
                res.append((i, j))
            i, j = i + 1, j + 1
            if i >= len(self) or j >= len(self[i]):
                break
        return res

    def entries_by_content(self, c):
        """
        Return the entries in ``self`` with content ``c``.

        EXAMPLES::

            sage: s = SkewTableau([[None,1,2],[3,4,5],[6]])
            sage: s.entries_by_content(0)
            [4]
            sage: s.entries_by_content(1)
            [1, 5]
            sage: s.entries_by_content(2)
            [2]
            sage: s.entries_by_content(-1)
            [3]
            sage: s.entries_by_content(-2)
            [6]
        """
        return [self[i][j] for i, j in self.cells_by_content(c)]

    def cells(self):
        """
        Return the cells in ``self``.

        EXAMPLES::

            sage: s = SkewTableau([[None,1,2],[3],[6]])
            sage: s.cells()
            [(0, 1), (0, 2), (1, 0), (2, 0)]
        """
        res = []
        for i in range(len(self)):
            for j in range(len(self[i])):
                if self[i][j] is not None:
                    res.append((i, j))
        return res

    def cells_containing(self, i):
        r"""
        Return the list of cells in which the letter ``i`` appears in the
        tableau ``self``. The list is ordered with cells appearing from
        left to right.

        Cells are given as pairs of coordinates `(a, b)`, where both
        rows and columns are counted from `0` (so `a = 0` means the cell
        lies in the leftmost column of the tableau, etc.).

        EXAMPLES::

            sage: t = SkewTableau([[None,None,3],[None,3,5],[4,5]])
            sage: t.cells_containing(5)
            [(2, 1), (1, 2)]
            sage: t.cells_containing(4)
            [(2, 0)]
            sage: t.cells_containing(2)
            []

            sage: t = SkewTableau([[None,None,None,None],[None,4,5],[None,5,6],[None,9],[None]])
            sage: t.cells_containing(2)
            []
            sage: t.cells_containing(4)
            [(1, 1)]
            sage: t.cells_containing(5)
            [(2, 1), (1, 2)]

            sage: SkewTableau([]).cells_containing(3)
            []

            sage: SkewTableau([[None,None],[None]]).cells_containing(3)
            []
        """
        cell_list = []
        for r in range(len(self) - 1, -1, -1):
            rth_row = self[r]
            for c, val in enumerate(rth_row):
                if val == i:
                    cell_list.append((r, c))
        return cell_list

    def is_k_tableau(self, k):
        r"""
        Checks whether ``self`` is a valid skew weak `k`-tableau.

        EXAMPLES::

            sage: t = SkewTableau([[None,2,3],[2,3],[3]])
            sage: t.is_k_tableau(3)
            True
            sage: t = SkewTableau([[None,1,3],[2,2],[3]])
            sage: t.is_k_tableau(3)
            False
        """
        shapes = self.to_chain()
        kshapes = [la.k_conjugate(k) for la in shapes]
        return all(kshapes[i + 1].contains(kshapes[i])
                   for i in range(len(shapes) - 1))


def _label_skew(list_of_cells, sk):
    """
    Return a filled-in standard skew tableau given an
    ordered list ``list_of_cells`` of the coordinates to fill in
    (as pairs) and an empty shape ``sk``.

    EXAMPLES::

        sage: import sage.combinat.skew_tableau as skew_tableau
        sage: l = [(0, 0), (1, 1), (1, 0), (0, 1)]
        sage: empty = [[None,None],[None,None]]
        sage: skew_tableau._label_skew(l, empty)
        [[1, 4], [3, 2]]
    """
    i = 1
    skew = [list(row) for row in sk]
    for row, column in list_of_cells:
        skew[row][column] = i
        i += 1
    return skew


class SkewTableaux(UniqueRepresentation, Parent):
    r"""
    Class of all skew tableaux.
    """
    def __init__(self, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SkewTableaux()
            sage: TestSuite(S).run()
        """
        if category is None:
            Parent.__init__(self, category=Sets())
        else:
            Parent.__init__(self, category=category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SkewTableaux()
            Skew tableaux
        """
        return "Skew tableaux"

    def _element_constructor_(self, st):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: S = SkewTableaux()
            sage: elt = S([[None,1],[2,3]]); elt
            [[None, 1], [2, 3]]
            sage: elt.parent() is S
            True
        """
        return self.element_class(self, st)

    Element = SkewTableau
    options = Tableaux.options

    def __contains__(self, x):
        """
        Checks if ``x`` is a skew tableau.

        EXAMPLES::

            sage: T = SkewTableau([[None, None, 1], [3], [4]])
            sage: T in SkewTableaux()
            True
            sage: [[None,1],[2,3]] in SkewTableaux()
            True
        """
        if isinstance(x, SkewTableau):
            return True
        try:
            self.element_class(self, x)
        except Exception:
            return False
        return True

    def from_expr(self, expr):
        """
        Return a :class:`SkewTableau` from a MuPAD-Combinat expr for
        a skew tableau.

        The first list in ``expr`` is the inner shape of the skew
        tableau. The second list are the entries in the rows of the skew
        tableau from bottom to top.

        Provided primarily for compatibility with MuPAD-Combinat.

        EXAMPLES::

            sage: SkewTableaux().from_expr([[1,1],[[5],[3,4],[1,2]]])
            [[None, 1, 2], [None, 3, 4], [5]]
        """
        skp = []
        outer = expr[1]
        inner = expr[0] + [0] * (len(outer) - len(expr[0]))

        for i in range(len(outer)):
            skp.append([None] * (inner[i]) + outer[-(i + 1)])

        return self.element_class(self, skp)

    def from_chain(self, chain):
        """
        Return the tableau corresponding to the chain of partitions.

        EXAMPLES::

            sage: SkewTableaux().from_chain([[1,1],[2,1],[3,1],[3,2],[3,3],[3,3,1]])
            [[None, 1, 2], [None, 3, 4], [5]]
        """
        shape = chain[-1]
        T = [[None]*r for r in shape]
        for i in range(1, len(chain)):
            la = chain[i]
            mu = chain[i - 1]
            mu += [0] * (len(la) - len(mu))

            for r in range(len(la)):
                for c in range(mu[r], la[r]):
                    T[r][c] = i

        return self.element_class(self, T)

    def from_shape_and_word(self, shape, word):
        """
        Return the skew tableau corresponding to the skew partition ``shape``
        and the word ``word`` obtained from the row reading.

        EXAMPLES::

            sage: t = SkewTableau([[None, 1, 3], [None, 2], [4]])
            sage: shape = t.shape()
            sage: word  = t.to_word()
            sage: SkewTableaux().from_shape_and_word(shape, word)
            [[None, 1, 3], [None, 2], [4]]
        """
        outer, inner = shape
        st = [[None] * row_length for row_length in outer]
        w_count = 0
        for i in reversed(range(len(inner), len(outer))):
            for j in range(outer[i]):
                st[i][j] = word[w_count]
                w_count += 1
        for i in reversed(range(len(inner))):
            for j in range(outer[i]):
                if j >= inner[i]:
                    st[i][j] = word[w_count]
                    w_count += 1
        return self.element_class(self, st)


class StandardSkewTableaux(SkewTableaux):
    """
    Standard skew tableaux.

    EXAMPLES::

        sage: S = StandardSkewTableaux(); S
        Standard skew tableaux
        sage: S.cardinality()
        +Infinity

    ::

        sage: S = StandardSkewTableaux(2); S
        Standard skew tableaux of size 2
        sage: S.cardinality()
        4

    ::

        sage: StandardSkewTableaux([[3, 2, 1], [1, 1]]).list()
        [[[None, 2, 3], [None, 4], [1]],
         [[None, 1, 2], [None, 3], [4]],
         [[None, 1, 2], [None, 4], [3]],
         [[None, 1, 3], [None, 4], [2]],
         [[None, 1, 4], [None, 3], [2]],
         [[None, 1, 4], [None, 2], [3]],
         [[None, 1, 3], [None, 2], [4]],
         [[None, 2, 4], [None, 3], [1]]]
    """
    @staticmethod
    def __classcall_private__(cls, skp=None):
        """
        Return the class of standard skew tableaux of skew shape ``skp``.

        EXAMPLES::

            sage: SST1 = StandardSkewTableaux([[3, 2, 1], [1, 1]])
            sage: SST2 = StandardSkewTableaux(SkewPartition([[3, 2, 1], [1, 1]]))
            sage: SST1 is SST2
            True
        """
        if skp is None:
            return StandardSkewTableaux_all()
        elif isinstance(skp, (int, Integer)):
            return StandardSkewTableaux_size(skp)
        elif skp in SkewPartitions():
            return StandardSkewTableaux_shape(skp)
        else:
            raise TypeError("Invalid argument")

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[None, 2], [1, 3]] in StandardSkewTableaux()
            True
            sage: [[None, 2], [2, 4]] in StandardSkewTableaux()
            False
            sage: [[None, 3], [2, 4]] in StandardSkewTableaux()
            False
            sage: [[None, 2], [1, 4]] in StandardSkewTableaux()
            False
        """
        if x not in SkewTableaux():
            return False

        return SkewTableau(x).is_standard()


class StandardSkewTableaux_all(StandardSkewTableaux):
    """
    Class of all standard skew tableaux.
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: s = StandardSkewTableaux()
            sage: TestSuite(s).run()
        """
        StandardSkewTableaux.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: StandardSkewTableaux()
            Standard skew tableaux
        """
        return "Standard skew tableaux"

    def __iter__(self):
        """
        Iterate through all standard skew tableaux having
        no empty rows (before nonempty rows) and no empty columns
        (before nonempty columns).

        EXAMPLES::

            sage: it = StandardSkewTableaux().__iter__()
            sage: [next(it) for x in range(10)]
            [[],
             [[1]],
             [[1, 2]], [[1], [2]], [[None, 2], [1]], [[None, 1], [2]],
             [[1, 2, 3]], [[1, 2], [3]], [[1, 3], [2]], [[None, 2, 3], [1]]]
        """
        n = 0
        while True:
            for st in StandardSkewTableaux_size(n):
                yield self.element_class(self, st)
            n += 1


class StandardSkewTableaux_size(StandardSkewTableaux):
    """
    Standard skew tableaux of a fixed size `n`.
    """
    def __init__(self, n):
        """
        EXAMPLES::

            sage: S = StandardSkewTableaux(3)
            sage: TestSuite(S).run()
        """
        self.n = n
        StandardSkewTableaux.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: StandardSkewTableaux(3)
            Standard skew tableaux of size 3
        """
        return "Standard skew tableaux of size %s" % self.n

    def cardinality(self):
        """
        EXAMPLES::

            sage: StandardSkewTableaux(1).cardinality()
            1
            sage: StandardSkewTableaux(2).cardinality()
            4
            sage: StandardSkewTableaux(3).cardinality()
            24
            sage: StandardSkewTableaux(4).cardinality()
            194
        """
        count = 0
        for skp in SkewPartitions(self.n):
            count += StandardSkewTableaux_shape(skp).cardinality()
        return count

    def __iter__(self):
        """
        Iterate through all standard skew tableaux of size `n` having
        no empty rows (before nonempty rows) and no empty columns
        (before nonempty columns). (The last two requirements
        ensure that the iterator terminates after finitely many steps.)

        EXAMPLES::

            sage: StandardSkewTableaux(2).list()
            [[[1, 2]], [[1], [2]], [[None, 2], [1]], [[None, 1], [2]]]

            sage: StandardSkewTableaux(3).list()
            [[[1, 2, 3]],
             [[1, 2], [3]], [[1, 3], [2]],
             [[None, 2, 3], [1]], [[None, 1, 2], [3]], [[None, 1, 3], [2]],
             [[None, 2], [1, 3]], [[None, 1], [2, 3]],
             [[None, None, 2], [1, 3]], [[None, None, 3], [1, 2]],
             [[None, None, 1], [2, 3]],
             [[1], [2], [3]], [[None, 2], [None, 3], [1]],
             [[None, 1], [None, 2], [3]], [[None, 1], [None, 3], [2]],
             [[None, 2], [1], [3]], [[None, 3], [1], [2]],
             [[None, 1], [2], [3]], [[None, None, 3], [None, 2], [1]],
             [[None, None, 2], [None, 3], [1]],
             [[None, None, 1], [None, 3], [2]],
             [[None, None, 1], [None, 2], [3]],
             [[None, None, 2], [None, 1], [3]],
             [[None, None, 3], [None, 1], [2]]]
        """
        for skp in SkewPartitions(self.n):
            for sst in StandardSkewTableaux_shape(skp):
                yield self.element_class(self, sst)


class StandardSkewTableaux_shape(StandardSkewTableaux):
    r"""
    Standard skew tableaux of a fixed skew shape `\lambda / \mu`.
    """
    @staticmethod
    def __classcall_private__(cls, skp):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = StandardSkewTableaux([[3, 2, 1], [1, 1]])
            sage: S2 = StandardSkewTableaux(SkewPartition([[3, 2, 1], [1, 1]]))
            sage: S is S2
            True
        """
        return super(StandardSkewTableaux_shape, cls).__classcall__(cls, SkewPartition(skp))

    def __init__(self, skp):
        """
        TESTS::

            sage: S = StandardSkewTableaux([[3, 2, 1], [1, 1]])
            sage: TestSuite(S).run()
        """
        self.skp = skp
        StandardSkewTableaux.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: StandardSkewTableaux([[3, 2, 1], [1, 1]])
            Standard skew tableaux of shape [3, 2, 1] / [1, 1]
        """
        return "Standard skew tableaux of shape %s" % repr(self.skp)

    def cardinality(self):
        """
        Return the number of standard skew tableaux with shape of the skew
        partition ``skp``. This uses a formula due to Aitken
        (see Cor. 7.16.3 of [Sta-EC2]_).

        EXAMPLES::

            sage: StandardSkewTableaux([[3, 2, 1], [1, 1]]).cardinality()
            8
        """
        outer, inner = self.skp
        m = len(outer)
        n = sum(outer) - sum(inner)
        outer = list(outer)
        inner = list(inner) + [0] * (m - len(inner))
        a = zero_matrix(QQ, m)
        for i in range(m):
            for j in range(m):
                v = outer[i] - inner[j] - i + j
                if v < 0:
                    a[i, j] = 0
                else:
                    a[i, j] = 1 / factorial(v)
        return ZZ(factorial(n) * a.det())

    def __iter__(self):
        """
        An iterator for all the standard skew tableaux whose shape is
        the skew partition ``skp``. The standard skew tableaux are
        ordered lexicographically by the word obtained from their row
        reading.

        EXAMPLES::

            sage: StandardSkewTableaux([[3, 2, 1], [1, 1]]).list()
            [[[None, 2, 3], [None, 4], [1]],
             [[None, 1, 2], [None, 3], [4]],
             [[None, 1, 2], [None, 4], [3]],
             [[None, 1, 3], [None, 4], [2]],
             [[None, 1, 4], [None, 3], [2]],
             [[None, 1, 4], [None, 2], [3]],
             [[None, 1, 3], [None, 2], [4]],
             [[None, 2, 4], [None, 3], [1]]]
        """
        dag = self.skp.to_dag(format="tuple")
        le_list = list(dag.topological_sort_generator())

        empty = [[None] * row_length for row_length in self.skp.outer()]

        for le in le_list:
            yield self.element_class(self, _label_skew(le, empty))


class SemistandardSkewTableaux(SkewTableaux):
    r"""
    Semistandard skew tableaux.

    This class can be initialized with several optional variables:
    the size of the skew tableaux (as a nameless integer variable),
    their shape (as a nameless skew partition variable), their
    weight (:meth:`~sage.combinat.skew_tableau.SkewTableau.weight`,
    as a nameless second variable after either the size or the
    shape) and their maximum entry (as an optional keyword variable
    called ``max_entry``, unless the weight has been specified). If
    neither the weight nor the maximum entry is specified, the
    maximum entry defaults to the size of the tableau.

    Note that "maximum entry" does not literally mean the highest
    entry; instead it is just an upper bound that no entry is
    allowed to surpass.

    EXAMPLES:

    The (infinite) class of all semistandard skew tableaux::

        sage: SemistandardSkewTableaux()
        Semistandard skew tableaux

    The (still infinite) class of all semistandard skew tableaux
    with maximum entry `2`::

        sage: SemistandardSkewTableaux(max_entry=2)
        Semistandard skew tableaux with maximum entry 2

    The class of all semistandard skew tableaux of given size `3`
    and maximum entry `3`::

        sage: SemistandardSkewTableaux(3)
        Semistandard skew tableaux of size 3 and maximum entry 3

    To set a different maximum entry::

        sage: SemistandardSkewTableaux(3, max_entry = 7)
        Semistandard skew tableaux of size 3 and maximum entry 7

    Specifying a shape::

        sage: SemistandardSkewTableaux([[2,1],[]])
        Semistandard skew tableaux of shape [2, 1] / [] and maximum entry 3

    Specifying both a shape and a maximum entry::

        sage: S = SemistandardSkewTableaux([[2,1],[1]], max_entry = 3); S
        Semistandard skew tableaux of shape [2, 1] / [1] and maximum entry 3
        sage: S.list()
        [[[None, 1], [1]],
         [[None, 2], [1]],
         [[None, 1], [2]],
         [[None, 3], [1]],
         [[None, 1], [3]],
         [[None, 2], [2]],
         [[None, 3], [2]],
         [[None, 2], [3]],
         [[None, 3], [3]]]

        sage: for n in range(5):
        ....:     print("{} {}".format(n, len(SemistandardSkewTableaux([[2,2,1],[1]], max_entry = n))))
        0 0
        1 0
        2 1
        3 9
        4 35

    Specifying a shape and a weight::

        sage: SemistandardSkewTableaux([[2,1],[]],[2,1])
        Semistandard skew tableaux of shape [2, 1] / [] and weight [2, 1]

    (the maximum entry is redundant in this case and thus is ignored).

    Specifying a size and a weight::

        sage: SemistandardSkewTableaux(3, [2,1])
        Semistandard skew tableaux of size 3 and weight [2, 1]

    .. WARNING::

        If the shape is not specified, the iterator of this class
        yields only skew tableaux whose shape is reduced, in the
        sense that there are no empty rows before the last nonempty
        row, and there are no empty columns before the last
        nonempty column. (Otherwise it would go on indefinitely.)

    .. WARNING::

        This class acts as a factory. The resulting classes are mainly
        useful for iteration. Do not rely on their containment tests,
        as they are not correct, e. g.::

            sage: SkewTableau([[None]]) in SemistandardSkewTableaux(2)
            True
    """
    @staticmethod
    def __classcall_private__(cls, p=None, mu=None, max_entry=None):
        """
        Return the correct parent based upon the input.

        EXAMPLES::

            sage: SSST1 = SemistandardSkewTableaux([[3, 2, 1], [1, 1]])
            sage: SSST2 = SemistandardSkewTableaux(SkewPartition([[3, 2, 1], [1, 1]]))
            sage: SSST1 is SSST2
            True
        """
        if p is None:
            if mu is None:
                return SemistandardSkewTableaux_all(max_entry)
            raise ValueError("You must specify either a size or a shape")

        if isinstance(p, (int, Integer)):
            if mu is None:
                return SemistandardSkewTableaux_size(p, max_entry)
            else:
                return SemistandardSkewTableaux_size_weight(p, mu)

        if p in SkewPartitions():
            if mu is None:
                return SemistandardSkewTableaux_shape(p, max_entry)
            else:
                return SemistandardSkewTableaux_shape_weight(p, mu)

        raise ValueError("Invalid input")

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[None, 2], [1, 3]] in SemistandardSkewTableaux()
            True
            sage: [[None, 2], [2, 4]] in SemistandardSkewTableaux()
            True
            sage: [[None, 3], [2, 4]] in SemistandardSkewTableaux()
            True
            sage: [[None, 2], [2, 4]] in SemistandardSkewTableaux()
            True
        """
        if x not in SkewTableaux():
            return False

        try:
            x = self.element_class(self, x)
        except Exception:
            return False
        return x.is_semistandard()


class SemistandardSkewTableaux_all(SemistandardSkewTableaux):
    """
    Class of all semistandard skew tableaux, possibly with a given
    maximum entry.
    """
    def __init__(self, max_entry):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SemistandardSkewTableaux()
            sage: TestSuite(S).run()

            sage: S = SemistandardSkewTableaux(3)
            sage: TestSuite(S).run()
        """
        SemistandardSkewTableaux.__init__(self, category=InfiniteEnumeratedSets())
        if max_entry is None:
            self.max_entry = PlusInfinity()
        else:
            self.max_entry = max_entry

    def _repr_(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux()
            Semistandard skew tableaux
        """
        if self.max_entry == PlusInfinity():
            return "Semistandard skew tableaux"
        return "Semistandard skew tableaux with maximum entry {}".format(self.max_entry)

    def __iter__(self):
        """
        Iterate over the elements of ``self``.

        EXAMPLES::

            sage: it = SemistandardSkewTableaux(max_entry = 5).__iter__()
            sage: [next(it) for x in range(12)]
            [[],
             [[1]],
             [[2]],
             [[3]],
             [[4]],
             [[5]],
             [[1, 1]],
             [[1, 2]],
             [[1, 3]],
             [[1, 4]],
             [[1, 5]],
             [[2, 2]]]

        If no max entry is specified, the iteration goes over all
        semistandard skew tableaux of size `n` with max entry `n`,
        for all `n`::

            sage: it = SemistandardSkewTableaux().__iter__()
            sage: [next(it) for x in range(10)]
            [[],
             [[1]],
             [[1, 1]],
             [[1, 2]],
             [[2, 2]],
             [[1], [2]],
             [[None, 1], [1]],
             [[None, 2], [1]],
             [[None, 1], [2]],
             [[None, 2], [2]]]
        """
        if self.max_entry == PlusInfinity():
            # Old behavior, kept here for backwards compatibility.
            # The usefulness of this iterator is questionable.
            n = 0
            while True:
                for ssst in SemistandardSkewTableaux_size(n, n):
                    yield self.element_class(self, ssst)
                n += 1
        else:
            n = 0
            while True:
                for ssst in SemistandardSkewTableaux_size(n, self.max_entry):
                    yield self.element_class(self, ssst)
                n += 1


class SemistandardSkewTableaux_size(SemistandardSkewTableaux):
    """
    Class of all semistandard skew tableaux of a fixed size `n`,
    possibly with a given maximum entry.
    """
    def __init__(self, n, max_entry):
        """
        EXAMPLES::

            sage: S = SemistandardSkewTableaux(3)
            sage: TestSuite(S).run()
        """
        self.n = n
        if max_entry is None:
            self.max_entry = n
        else:
            self.max_entry = max_entry
        SemistandardSkewTableaux.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(3)
            Semistandard skew tableaux of size 3 and maximum entry 3
            sage: SemistandardSkewTableaux(3, max_entry=8)
            Semistandard skew tableaux of size 3 and maximum entry 8
        """
        return "Semistandard skew tableaux of size %s and maximum entry %s" % (repr(self.n), repr(self.max_entry))

    def cardinality(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(2).cardinality()
            8
        """
        count = 0
        for p in SkewPartitions(self.n):
            count += SemistandardSkewTableaux_shape(p, self.max_entry).cardinality()
        return count

    def __iter__(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(2).list()
            [[[1, 1]],
             [[1, 2]],
             [[2, 2]],
             [[1], [2]],
             [[None, 1], [1]],
             [[None, 2], [1]],
             [[None, 1], [2]],
             [[None, 2], [2]]]
        """
        for p in SkewPartitions(self.n):
            for ssst in SemistandardSkewTableaux_shape(p, self.max_entry):
                yield self.element_class(self, ssst)


class SemistandardSkewTableaux_size_weight(SemistandardSkewTableaux):
    r"""
    Class of semistandard tableaux of a fixed size `n` and weight `\mu`.
    """
    @staticmethod
    def __classcall_private__(cls, n, mu):
        """
        Normalize our input to ensure we have a unique representation.

        EXAMPLES::

            sage: S = SemistandardSkewTableaux(3, [2,1])
            sage: S2 = SemistandardSkewTableaux(int(3), (2,1))
            sage: S is S2
            True
        """
        return super(SemistandardSkewTableaux_size_weight, cls).__classcall__(cls, n, tuple(mu))

    def __init__(self, n, mu):
        """
        EXAMPLES::

            sage: S = SemistandardSkewTableaux(3,[2,1])
            sage: TestSuite(S).run()
        """
        self.n = n
        self.mu = mu
        SemistandardSkewTableaux.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(3,[2,1])
            Semistandard skew tableaux of size 3 and weight [2, 1]
        """
        return "Semistandard skew tableaux of size %s and weight %s" % (repr(self.n), list(self.mu))

    def cardinality(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(2,[1,1]).cardinality()
            4
        """
        count = 0
        for p in SkewPartitions(self.n):
            count += SemistandardSkewTableaux_shape_weight(p, self.mu).cardinality()
        return count

    def __iter__(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux(2,[1,1]).list()
            [[[1, 2]], [[1], [2]], [[None, 2], [1]], [[None, 1], [2]]]
        """
        for p in SkewPartitions(self.n):
            for ssst in SemistandardSkewTableaux_shape_weight(p, self.mu):
                yield self.element_class(self, ssst)


class SemistandardSkewTableaux_shape(SemistandardSkewTableaux):
    r"""
    Class of semistandard skew tableaux of a fixed skew shape
    `\lambda / \mu` with a given max entry.

    A semistandard skew tableau with max entry `i` is required to have all
    its entries less or equal to `i`. It is not required to actually
    contain an entry `i`.

    INPUT:

    - ``p`` -- A skew partition

    - ``max_entry`` -- The max entry; defaults to the size of ``p``.

    .. WARNING::

        Input is not checked; please use :class:`SemistandardSkewTableaux` to
        ensure the options are properly parsed.
    """
    @staticmethod
    def __classcall_private__(cls, p, max_entry=None):
        """
        Normalize our input to ensure we have a unique representation.

        EXAMPLES::

            sage: S = SemistandardSkewTableaux([[2,1],[]])
            sage: S2 = SemistandardSkewTableaux(SkewPartition([[2,1],[]]))
            sage: S is S2
            True
        """
        if max_entry is None:
            max_entry = sum(p[0]) - sum(p[1])
        return super(SemistandardSkewTableaux_shape, cls).__classcall__(cls, SkewPartition(p), max_entry)

    def __init__(self, p, max_entry):
        """
        EXAMPLES::

            sage: S = SemistandardSkewTableaux([[2,1],[]])
            sage: S == loads(dumps(S))
            True
            sage: TestSuite(S).run()
        """
        self.p = p
        self.max_entry = max_entry
        SemistandardSkewTableaux.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]])
            Semistandard skew tableaux of shape [2, 1] / [] and maximum entry 3
        """
        return "Semistandard skew tableaux of shape %s and maximum entry %s" % (repr(self.p), repr(self.max_entry))

    def cardinality(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]]).cardinality()
            8
            sage: SemistandardSkewTableaux([[2,1],[]], max_entry=2).cardinality()
            2
        """
        count = 0
        for mu in IntegerVectors(self.p.size(), self.max_entry):
            count += SemistandardSkewTableaux_shape_weight(self.p, mu).cardinality()
        return count

    def __iter__(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]]).list()
            [[[1, 1], [2]],
             [[1, 1], [3]],
             [[1, 2], [2]],
             [[1, 3], [2]],
             [[1, 2], [3]],
             [[1, 3], [3]],
             [[2, 2], [3]],
             [[2, 3], [3]]]
            sage: from sage.combinat.skew_tableau import SemistandardSkewTableaux_shape
            sage: SemistandardSkewTableaux_shape([[2,1],[]], max_entry=2).list()
            [[[1, 1], [2]], [[1, 2], [2]]]
        """
        for mu in IntegerVectors(self.p.size(), self.max_entry):
            for ssst in SemistandardSkewTableaux_shape_weight(self.p, mu):
                yield self.element_class(self, ssst)


class SemistandardSkewTableaux_shape_weight(SemistandardSkewTableaux):
    r"""
    Class of semistandard skew tableaux of a fixed skew shape `\lambda / \nu`
    and weight `\mu`.
    """
    @staticmethod
    def __classcall_private__(cls, p, mu):
        """
        Normalize our input to ensure we have a unique representation.

        EXAMPLES::

            sage: S = SemistandardSkewTableaux([[2,1],[]], [2,1])
            sage: S2 = SemistandardSkewTableaux(SkewPartition([[2,1],[]]), (2,1))
            sage: S is S2
            True
        """
        p = SkewPartition(p)
        mu = tuple(mu)
        return super(SemistandardSkewTableaux_shape_weight, cls).__classcall__(cls, p, mu)

    def __init__(self, p, mu):
        """
        EXAMPLES::

            sage: S = SemistandardSkewTableaux([[2,1],[]],[2,1])
            sage: S == loads(dumps(S))
            True
            sage: TestSuite(S).run()
        """
        self.p = p
        self.mu = mu
        SemistandardSkewTableaux.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]],[2,1])
            Semistandard skew tableaux of shape [2, 1] / [] and weight [2, 1]
        """
        return "Semistandard skew tableaux of shape %s and weight %s" % (repr(self.p), list(self.mu))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: SemistandardSkewTableaux([[2,1],[]],[2,1]).list()
            [[[1, 1], [2]]]
        """
        from .ribbon_tableau import RibbonTableaux_shape_weight_length
        for x in RibbonTableaux_shape_weight_length(self.p, self.mu, 1):
            yield self.element_class(self, x)


class SkewTableau_class(SkewTableau):
    """
    This exists solely for unpickling ``SkewTableau_class`` objects.
    """
    def __setstate__(self, state):
        r"""
        Unpickle old ``SkewTableau_class`` objects.

        TESTS::

            sage: loads(b'x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1+H,*\xc9,\xc9\xcc\xcf\xe3\n\x80\xb1\xe2\x93s\x12\x8b\x8b\xb9\n\x195\x1b\x0b\x99j\x0b\x995BY\xe33\x12\x8b3\nY\xfc\x80\xac\x9c\xcc\xe2\x92B\xd6\xd8B6\r\x88IE\x99y\xe9\xc5z\x99y%\xa9\xe9\xa9E\\\xb9\x89\xd9\xa9\xf10N!{(\xa3qkP!G\x06\x90a\x04dp\x82\x18\x86@\x06Wji\x92\x1e\x00x0.\xb5')
            [3, 2, 1]
            sage: loads(dumps( SkewTableau([[1,1], [3,2,1]]) ))  # indirect doctest
            [[1, 1], [3, 2, 1]]
        """
        self.__class__ = SkewTableau
        self.__init__(SkewTableaux(), state['_list'])


# October 2012: fixing outdated pickles which use the classes being deprecated


register_unpickle_override('sage.combinat.skew_tableau', 'StandardSkewTableaux_n', StandardSkewTableaux_size)
register_unpickle_override('sage.combinat.skew_tableau', 'SemistandardSkewTableaux_n', SemistandardSkewTableaux_size)
register_unpickle_override('sage.combinat.skew_tableau', 'SemistandardSkewTableaux_nmu', SemistandardSkewTableaux_size_weight)
register_unpickle_override('sage.combinat.skew_tableau', 'SemistandardSkewTableaux_p', SemistandardSkewTableaux_shape)
register_unpickle_override('sage.combinat.skew_tableau', 'SemistandardSkewTableaux_pmu', SemistandardSkewTableaux_shape_weight)
# July 2013: But wait, there more!
register_unpickle_override('sage.combinat.skew_tableau', 'StandardSkewTableaux_skewpartition', StandardSkewTableaux_shape)
register_unpickle_override('sage.combinat.skew_tableau', 'SkewTableau_class', SkewTableau_class)
