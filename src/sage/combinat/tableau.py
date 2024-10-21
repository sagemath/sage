# sage_setup: distribution = sagemath-combinat
r"""
Tableaux

AUTHORS:

- Mike Hansen (2007): initial version

- Jason Bandlow (2011): updated to use Parent/Element model, and many
  minor fixes

- Andrew Mathas (2012-13): completed the transition to the parent/element model
  begun by Jason Bandlow

- Travis Scrimshaw (11-22-2012): Added tuple options, changed ``*katabolism*``
  to ``*catabolism*``. Cleaned up documentation.

- Andrew Mathas (2016-08-11): Row standard tableaux added

- Oliver Pechenik (2018): Added increasing tableaux.

This file consists of the following major classes:

Element classes:

* :class:`Tableau`
* :class:`SemistandardTableau`
* :class:`StandardTableau`
* :class:`RowStandardTableau`
* :class:`IncreasingTableau`

Factory classes:

* :class:`Tableaux`
* :class:`SemistandardTableaux`
* :class:`StandardTableaux`
* :class:`RowStandardTableaux`
* :class:`IncreasingTableaux`

Parent classes:

* :class:`Tableaux_all`
* :class:`Tableaux_size`
* :class:`SemistandardTableaux_all` (facade class)
* :class:`SemistandardTableaux_size`
* :class:`SemistandardTableaux_size_inf`
* :class:`SemistandardTableaux_size_weight`
* :class:`SemistandardTableaux_shape`
* :class:`SemistandardTableaux_shape_inf`
* :class:`SemistandardTableaux_shape_weight`
* :class:`StandardTableaux_all` (facade class)
* :class:`StandardTableaux_size`
* :class:`StandardTableaux_shape`
* :class:`IncreasingTableaux_all` (facade class)
* :class:`IncreasingTableaux_size`
* :class:`IncreasingTableaux_size_inf`
* :class:`IncreasingTableaux_size_weight`
* :class:`IncreasingTableaux_shape`
* :class:`IncreasingTableaux_shape_inf`
* :class:`IncreasingTableaux_shape_weight`
* :class:`RowStandardTableaux_all` (facade class)
* :class:`RowStandardTableaux_size`
* :class:`RowStandardTableaux_shape`

For display options, see :meth:`Tableaux.options`.

.. TODO::

    - Move methods that only apply to semistandard tableaux from tableau to
      semistandard tableau

    - Copy/move functionality to skew tableaux

    - Add a class for tableaux of a given shape (eg Tableaux_shape)
"""

# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2011 Jason Bandlow <jbandlow@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from itertools import repeat

import sage.misc.prandom as random

from sage.arith.misc import binomial, factorial, multinomial
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_cat import Sets
from sage.combinat import permutation
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.composition import Composition, Compositions
from sage.combinat.integer_vector import IntegerVectors, integer_vectors_nk_fast_iter
from sage.combinat.subset import powerset
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.misc_c import prod
from sage.misc.lazy_import import lazy_import
from sage.misc.persist import register_unpickle_override
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.infinity import PlusInfinity
from sage.rings.integer import Integer
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.global_options import GlobalOptions
from sage.structure.list_clone import ClonableList
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp, richcmp_method
from sage.structure.unique_representation import UniqueRepresentation

lazy_import('sage.combinat.posets.posets', 'Poset')
lazy_import('sage.groups.perm_gps.permgroup', 'PermutationGroup')
lazy_import('sage.libs.symmetrica', 'all', as_='symmetrica')


@richcmp_method
class Tableau(ClonableList, metaclass=InheritComparisonClasscallMetaclass):
    """
    A class to model a tableau.

    INPUT:

    - ``t`` -- a Tableau, a list of iterables, or an empty list

    OUTPUT: a Tableau object constructed from ``t``

    A tableau is abstractly a mapping from the cells in a partition to
    arbitrary objects (called entries). It is often represented as a
    finite list of nonempty lists (or, more generally an iterator of
    iterables) of weakly decreasing lengths. This list,
    in particular, can be empty, representing the empty tableau.

    Note that Sage uses the English convention for partitions and
    tableaux; the longer rows are displayed on top.

    EXAMPLES::

        sage: t = Tableau([[1,2,3],[4,5]]); t
        [[1, 2, 3], [4, 5]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty printing
        1 2 3
        4 5
        sage: t.is_standard()
        True

        sage: Tableau([['a','c','b'],[[],(2,1)]])
        [['a', 'c', 'b'], [[], (2, 1)]]
        sage: Tableau([]) # The empty tableau
        []

    When using code that will generate a lot of tableaux, it is slightly more
    efficient to construct a Tableau from the appropriate Parent object::

        sage: T = Tableaux()
        sage: T([[1, 2, 3], [4, 5]])
        [[1, 2, 3], [4, 5]]

    .. SEEALSO::

        - :class:`Tableaux`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`

    TESTS::

        sage: Tableau([[1],[2,3]])
        Traceback (most recent call last):
        ...
        ValueError: a tableau must be a list of iterables of weakly decreasing length
        sage: Tableau([1,2,3])
        Traceback (most recent call last):
        ...
        ValueError: a tableau must be a list of iterables
    """
    @staticmethod
    def __classcall_private__(cls, t):
        r"""
        This ensures that a tableau is only ever constructed as an
        ``element_class`` call of an appropriate parent.

        TESTS::

            sage: t = Tableau([[1,1],[1]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Tableaux
            sage: t.category()
            Category of elements of Tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.Tableaux_all_with_category.element_class'>
        """
        if isinstance(t, cls):
            return t

        # We must verify ``t`` is a list of iterables, and also
        # normalize it to be a list of tuples.
        try:
            t = [tuple(_) for _ in t]
        except TypeError:
            raise ValueError("a tableau must be a list of iterables")

        return Tableaux_all().element_class(Tableaux_all(), t)

    def __init__(self, parent, t, check=True):
        r"""
        Initialize a tableau.

        TESTS::

            sage: t = Tableaux()([[1,1],[1]])
            sage: s = Tableaux(3)([[1,1],[1]])
            sage: s == t
            True
            sage: t.parent()
            Tableaux
            sage: s.parent()
            Tableaux of size 3
            sage: r = Tableaux()(s); r.parent()
            Tableaux
            sage: s is t # identical tableaux are distinct objects
            False

        A tableau is shallowly immutable. See :issue:`15862`. The entries
        themselves may be mutable objects, though in that case the
        resulting Tableau should be unhashable. ::

            sage: T = Tableau([[1,2],[2]])
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
        if isinstance(t, Tableau):
            # Since we are (supposed to be) immutable, we can share the underlying data
            ClonableList.__init__(self, parent, t, check=False)
            return

        # Normalize t to be a list of tuples.
        t = [tuple(_) for _ in t]

        ClonableList.__init__(self, parent, t, check=check)
        # This dispatches the input verification to the :meth:`check`
        # method.

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` to ``other``.

        .. TODO::

            This overwrites the comparison check of
            :class:`~sage.structure.list_clone.ClonableList`
            in order to circumvent the coercion framework.
            Eventually this should be solved more elegantly,
            for example along the lines of what was done for
            `k`-tableaux.

            For now, this compares two elements by their underlying
            defining lists.

        INPUT:

        - ``other`` -- the element that ``self`` is compared to

        OUTPUT: boolean

        TESTS::

            sage: t = Tableau([[1,2]])
            sage: t == 0
            False
            sage: t == Tableaux(2)([[1,2]])
            True

            sage: s = Tableau([[2,3],[1]])
            sage: s != []
            True

            sage: t < s
            True
            sage: s < t
            False
            sage: s > t
            True
        """
        if isinstance(other, Tableau):
            return richcmp(list(self), list(other), op)
        else:
            return richcmp(list(self), other, op)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: t = Tableau([[1,1],[2]])
            sage: hash(tuple(t)) == hash(t)
            True
        """
        return hash(tuple(self))

    def check(self):
        r"""
        Check that ``self`` is a valid straight-shape tableau.

        EXAMPLES::

            sage: t = Tableau([[1,1],[2]])
            sage: t.check()

            sage: t = Tableau([[None, None, 1], [2, 4], [3, 4, 5]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: a tableau must be a list of iterables of weakly decreasing length
        """
        # Check that it has partition shape. That's all we require from a
        # general tableau.
        lens = [len(row) for row in self]
        for (a, b) in zip(lens, lens[1:]):
            if a < b:
                raise ValueError("a tableau must be a list of iterables of weakly decreasing length")
        if lens and lens[-1] == 0:
            raise ValueError("a tableau must not have empty rows")

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: Tableaux.options.display="list"
            sage: t
            [[1, 2, 3], [4, 5]]
            sage: Tableaux.options.display="array"
            sage: t
              1  2  3
              4  5
            sage: Tableaux.options.display="compact"; t
            1,2,3/4,5
            sage: Tableaux.options._reset()
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def _repr_list(self):
        """
        Return a string representation of ``self`` as a list.

        EXAMPLES::

            sage: T = Tableau([[1,2,3],[4,5]])
            sage: T._repr_list()
            '[[1, 2, 3], [4, 5]]'
        """
        return repr([list(_) for _ in self])

    # See #18024. CombinatorialObject provided __str__, though ClonableList
    # doesn't. Emulate the old functionality. Possibly remove when
    # CombinatorialObject is removed.
    __str__ = _repr_list

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as an array.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: print(t._repr_diagram())
              1  2  3
              4  5
            sage: Tableaux.options.convention="french"
            sage: print(t._repr_diagram())
              4  5
              1  2  3
            sage: Tableaux.options.convention='russian'
            sage: print(t._repr_diagram())
               5    3
             4    2
                1
            sage: Tableaux.options._reset()

        TESTS:

        Check that :issue:`20768` is fixed::

            sage: T = Tableau([[1523, 1, 2],[1,12341, -2]])
            sage: T.pp()
             1523     1  2
                1 12341 -2
        """
        if not self:
            return "  -"

        # Get the widths of the columns
        str_tab = [[str(data) for data in row] for row in self]
        col_widths = [2]*len(str_tab[0])
        for row in str_tab:
            for i, e in enumerate(row):
                col_widths[i] = max(col_widths[i], len(e))

        if self.parent().options('convention') == "Russian":
            col_width = max(col_widths) + 1
            max_height = max(a + len(val) for a, val in enumerate(str_tab))
            str_list = []
            for i in range(max_height):
                st = ' ' * ((max_height - i - 1) * col_width)
                for a in range(i + 1):
                    b = i - a
                    if len(str_tab[b:]) > 0 and len(str_tab[b][a:]) > 0:
                        st += str_tab[b][a].rjust(col_width, ' ').ljust(col_width * 2 - 1, ' ')
                    else:
                        st += ' ' * (col_width * 2 - 1)
                str_list.append(st)
            import re
            mm = min(len(re.search('^ +', sline)[0]) for sline in str_list) - 1
            str_list = [sline[mm:] for sline in str_list]
            str_list.reverse()
            return '\n'.join(str_list)

        if self.parent().options('convention') == "French":
            str_tab = reversed(str_tab)

        return "\n".join(" "
                         + " ".join("{:>{width}}".format(e, width=col_widths[i])
                                    for i, e in enumerate(row))
                         for row in str_tab)

    def _repr_compact(self):
        """
        Return a compact string representation of ``self``.

        EXAMPLES::

            sage: Tableau([[1,2,3],[4,5]])._repr_compact()
            '1,2,3/4,5'
            sage: Tableau([])._repr_compact()
            '-'
        """
        if not self:
            return '-'
        return '/'.join(','.join('%s' % r for r in row) for row in self)

    def _ascii_art_(self):
        r"""
        TESTS::

            sage: ascii_art(list(StandardTableaux(3)))
            [                              1 ]
            [              1  3    1  2    2 ]
            [   1  2  3,   2   ,   3   ,   3 ]
            sage: Tableaux.options(ascii_art='compact')
            sage: ascii_art(list(StandardTableaux(3)))
            [                        |1| ]
            [          |1|3|  |1|2|  |2| ]
            [ |1|2|3|, |2|  , |3|  , |3| ]
            sage: Tableaux.options(convention='french', ascii_art='table')
            sage: ascii_art(list(StandardTableaux(3)))
            [                                      +---+ ]
            [                                      | 3 | ]
            [                +---+      +---+      +---+ ]
            [                | 2 |      | 3 |      | 2 | ]
            [ +---+---+---+  +---+---+  +---+---+  +---+ ]
            [ | 1 | 2 | 3 |  | 1 | 3 |  | 1 | 2 |  | 1 | ]
            [ +---+---+---+, +---+---+, +---+---+, +---+ ]
            sage: Tableaux.options(ascii_art='repr')
            sage: ascii_art(list(StandardTableaux(3)))
            [                              3 ]
            [              2       3       2 ]
            [   1  2  3,   1  3,   1  2,   1 ]
            sage: Tableaux.options(convention='russian', ascii_art='table')
            sage: ascii_art(list(StandardTableaux(3)))
            [       / \                             / \      ]
            [      / 3 /                           \ 3 \     ]
            [     / \ /     / \ / \     / \ / \     \ / \    ]
            [    / 2 /     \ 2 X 3 /   \ 3 X 2 /     \ 2 \   ]
            [   / \ /       \ / \ /     \ / \ /       \ / \  ]
            [  \ 1 /         \ 1 /       \ 1 /         \ 1 / ]
            [   \ /     ,     \ /   ,     \ /   ,       \ /  ]
            sage: Tableaux.options(ascii_art='repr')
            sage: ascii_art(list(StandardTableaux(3)))
            [      3                           3             ]
            [    2       2    3     3    2        2          ]
            [  1      ,     1    ,     1    ,        1       ]
            sage: Tableaux.options._reset()
        """
        ascii = self.parent().options._dispatch(self, '_ascii_art_', 'ascii_art')
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(ascii.splitlines())

    def _unicode_art_(self):
        r"""
        TESTS::

            sage: unicode_art(Tableau([[1,2,3],[4],[5]]))
            ┌───┬───┬───┐
            │ 1 │ 2 │ 3 │
            ├───┼───┴───┘
            │ 4 │
            ├───┤
            │ 5 │
            └───┘
            sage: unicode_art(Tableau([]))
            <BLANKLINE>
        """
        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt(self._ascii_art_table(use_unicode=True).splitlines())

    _ascii_art_repr = _repr_diagram

    def _ascii_art_table(self, use_unicode=False):
        r"""
        TESTS:

        We check that :issue:`16487` is fixed::

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: print(t._ascii_art_table())
            +---+---+---+
            | 1 | 2 | 3 |
            +---+---+---+
            | 4 | 5 |
            +---+---+
            sage: Tableaux.options.convention="french"
            sage: print(t._ascii_art_table())
            +---+---+
            | 4 | 5 |
            +---+---+---+
            | 1 | 2 | 3 |
            +---+---+---+
            sage: Tableaux.options.convention="russian"
            sage: print(t._ascii_art_table())
                / \ / \
               / 5 X 3 /
              / \ / \ /
             \ 4 X 2 /
              \ / \ /
               \ 1 /
                \ /
            sage: t = Tableau([]); print(t._ascii_art_table())
            <BLANKLINE>

            sage: Tableaux.options.convention="french"
            sage: print(t._ascii_art_table())
            <BLANKLINE>
            sage: Tableaux.options._reset()

            sage: t = Tableau([[1,2,3,10,15],[12,15,17]])
            sage: print(t._ascii_art_table())
            +----+----+----+----+----+
            |  1 |  2 |  3 | 10 | 15 |
            +----+----+----+----+----+
            | 12 | 15 | 17 |
            +----+----+----+

            sage: t = Tableau([[1,2,15,7],[12,5,6],[8,10],[9]])
            sage: Tableaux.options(ascii_art='table')
            sage: ascii_art(t)
            +----+----+----+---+
            |  1 |  2 | 15 | 7 |
            +----+----+----+---+
            | 12 |  5 |  6 |
            +----+----+----+
            |  8 | 10 |
            +----+----+
            |  9 |
            +----+
            sage: Tableaux.options.convention='french'
            sage: ascii_art(t)
            +----+
            |  9 |
            +----+----+
            |  8 | 10 |
            +----+----+----+
            | 12 |  5 |  6 |
            +----+----+----+---+
            |  1 |  2 | 15 | 7 |
            +----+----+----+---+
            sage: Tableaux.options.convention="russian"
            sage: ascii_art(t)
               / \   / \   / \   / \
              /   \ /   \ /   \ /   \
             \  9  X 10  X  6  X  7  /
              \   / \   / \   / \   /
               \ /   \ /   \ /   \ /
                \  8  X  5  X 15  /
                 \   / \   / \   /
                  \ /   \ /   \ /
                   \ 12  X  2  /
                    \   / \   /
                     \ /   \ /
                      \  1  /
                       \   /
                        \ /
            sage: Tableaux.options._reset()

        Unicode version::

            sage: t = Tableau([[1,2,15,7],[12,5],[8,10],[9]])
            sage: print(t._ascii_art_table(use_unicode=True))
            ┌────┬────┬────┬───┐
            │ 1  │ 2  │ 15 │ 7 │
            ├────┼────┼────┴───┘
            │ 12 │ 5  │
            ├────┼────┤
            │ 8  │ 10 │
            ├────┼────┘
            │ 9  │
            └────┘
            sage: Tableaux().options.convention='french'
            sage: t = Tableau([[1,2,15,7],[12,5],[8,10],[9]])
            sage: print(t._ascii_art_table(use_unicode=True))
            ┌────┐
            │ 9  │
            ├────┼────┐
            │ 8  │ 10 │
            ├────┼────┤
            │ 12 │ 5  │
            ├────┼────┼────┬───┐
            │ 1  │ 2  │ 15 │ 7 │
            └────┴────┴────┴───┘
            sage: Tableaux.options.convention="russian"
            sage: print(t._ascii_art_table(use_unicode=True))
               ╱ ╲   ╱ ╲         ╱ ╲
              ╱   ╲ ╱   ╲       ╱   ╲
             ╲  9  ╳ 10  ╲     ╱  7  ╱
              ╲   ╱ ╲   ╱ ╲   ╱ ╲   ╱
               ╲ ╱   ╲ ╱   ╲ ╱   ╲ ╱
                ╲  8  ╳  5  ╳ 15  ╱
                 ╲   ╱ ╲   ╱ ╲   ╱
                  ╲ ╱   ╲ ╱   ╲ ╱
                   ╲ 12  ╳  2  ╱
                    ╲   ╱ ╲   ╱
                     ╲ ╱   ╲ ╱
                      ╲  1  ╱
                       ╲   ╱
                        ╲ ╱
            sage: Tableaux.options._reset()
        """
        from sage.combinat.output import ascii_art_table
        self.parent().options('convention')
        return ascii_art_table(self, use_unicode=use_unicode,
                               convention=self.parent().options('convention'))

    def _ascii_art_compact(self):
        r"""
        TESTS:

        We check that :issue:`16487` is fixed::

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: print(t._ascii_art_compact())
            |1|2|3|
            |4|5|
            sage: Tableaux.options.convention="french"
            sage: print(t._ascii_art_compact())
            |4|5|
            |1|2|3|
            sage: Tableaux.options.convention="russian"
            sage: print(t._ascii_art_compact())
              /5X3/
             \4X2/
              \1/
            sage: Tableaux.options._reset()

            sage: t = Tableau([[1,2,3,10,15],[12,15,17]])
            sage: print(t._ascii_art_compact())
            |1 |2 |3 |10|15|
            |12|15|17|

            sage: t = Tableau([])
            sage: print(t._ascii_art_compact())
            .
        """
        if not self:
            return "."

        if self.parent().options('convention') == "Russian":
            from sage.combinat.output import ascii_art_table_russian
            return ascii_art_table_russian(self, compact=True)
        if self.parent().options('convention') == "English":
            T = self
        else:
            T = reversed(self)

        # Get the widths of the columns
        str_tab = [[str(_) for _ in row] for row in T]
        col_widths = [1]*len(self[0])
        for row in str_tab:
            for i, e in enumerate(row):
                col_widths[i] = max(col_widths[i], len(e))

        return "\n".join("|"
                         + "|".join("{:^{width}}".format(e, width=col_widths[i])
                                    for i, e in enumerate(row))
                         + "|" for row in str_tab)

    def _latex_(self):
        r"""
        Return a LaTeX version of ``self``.

        EXAMPLES::

            sage: t = Tableau([[1,1,2],[2,3],[3]])
            sage: latex(t)    # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{1-3}
            \lr{1}&\lr{1}&\lr{2}\\\cline{1-3}
            \lr{2}&\lr{3}\\\cline{1-2}
            \lr{3}\\\cline{1-1}
            \end{array}$}
            }
            sage: Tableaux.options.convention="french"
            sage: latex(t)    # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[t]{*{3}c}\cline{1-1}
            \lr{3}\\\cline{1-2}
            \lr{2}&\lr{3}\\\cline{1-3}
            \lr{1}&\lr{1}&\lr{2}\\\cline{1-3}
            \end{array}$}
            }
            sage: Tableaux.options.convention="russian"
            sage: latex(t)    # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{\rotatebox{45}{$\begin{array}[t]{*{3}c}\cline{1-1}
            \lr{\rotatebox{-45}{3}}\\\cline{1-2}
            \lr{\rotatebox{-45}{2}}&\lr{\rotatebox{-45}{3}}\\\cline{1-3}
            \lr{\rotatebox{-45}{1}}&\lr{\rotatebox{-45}{1}}&\lr{\rotatebox{-45}{2}}\\\cline{1-3}
            \end{array}$}}
            }

            sage: Tableaux.options._reset()
        """
        return self.parent().options._dispatch(self, '_latex_', 'latex')

    _latex_list = _repr_list

    def _latex_diagram(self):
        r"""
        Return a LaTeX representation of ``self`` as a Young diagram.

        EXAMPLES::

            sage: t = Tableau([[1,1,2],[2,3],[3]])
            sage: print(t._latex_diagram())
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{1-3}
            \lr{1}&\lr{1}&\lr{2}\\\cline{1-3}
            \lr{2}&\lr{3}\\\cline{1-2}
            \lr{3}\\\cline{1-1}
            \end{array}$}
            }
        """
        if len(self) == 0:
            return "{\\emptyset}"
        from sage.combinat.output import tex_from_array
        return tex_from_array(self)

    def __truediv__(self, t):
        """
        Return the skew tableau ``self``/``t``, where ``t`` is a partition
        contained in the shape of ``self``.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[3,4],[5]])
            sage: t/[1,1]
            [[None, 2, 3], [None, 4], [5]]
            sage: t/[3,1]
            [[None, None, None], [None, 4], [5]]
            sage: t/[2,1,1,1]
            Traceback (most recent call last):
            ...
            ValueError: the shape of the tableau must contain the partition
        """
        from sage.combinat.partition import _Partitions
        # if t is a list, convert it to a partition first
        if isinstance(t, list):
            t = _Partitions(t)

        # Check to make sure that tableau shape contains t
        if not self.shape().contains(t):
            raise ValueError("the shape of the tableau must contain the partition")

        st = [list(row) for row in self]    # create deep copy of t

        for i, t_i in enumerate(t):
            st_i = st[i]
            for j in range(t_i):
                st_i[j] = None

        from sage.combinat.skew_tableau import SkewTableau
        return SkewTableau(st)

    def __call__(self, *cell):
        r"""

        INPUT:

        - ``cell`` -- a pair of integers, tuple, or list specifying a cell in
          the tableau

        OUTPUT: the value in the corresponding cell

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: t(1,0)
            4
            sage: t((1,0))
            4
            sage: t(3,3)
            Traceback (most recent call last):
            ...
            IndexError: the cell (3,3) is not contained in [[1, 2, 3], [4, 5]]
        """
        try:
            i, j = cell
        except ValueError:
            i, j = cell[0]

        try:
            return self[i][j]
        except IndexError:
            raise IndexError("the cell (%d,%d) is not contained in %s" % (i, j, repr(self)))

    def level(self):
        """
        Return the level of ``self``, which is always 1.

        This function exists mainly for compatibility with
        :class:`TableauTuple`.

        EXAMPLES::

            sage: Tableau([[1,2,3],[4,5]]).level()
            1
        """
        return 1

    def components(self):
        """
        This function returns a list containing itself. It exists mainly for
        compatibility with :class:`TableauTuple` as it allows constructions
        like the example below.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: for s in t.components(): print(s.to_list())
            [[1, 2, 3], [4, 5]]
        """
        return [self]

    @combinatorial_map(name='shape')
    def shape(self):
        r"""
        Return the shape of a tableau ``self``.

        EXAMPLES::

            sage: Tableau([[1,2,3],[4,5],[6]]).shape()
            [3, 2, 1]
        """
        from sage.combinat.partition import Partition
        return Partition([len(row) for row in self])

    def size(self):
        """
        Return the size of the shape of the tableau ``self``.

        EXAMPLES::

            sage: Tableau([[1, 4, 6], [2, 5], [3]]).size()
            6
            sage: Tableau([[1, 3], [2, 4]]).size()
            4
        """
        return sum([len(row) for row in self])

    def corners(self):
        """
        Return the corners of the tableau ``self``.

        EXAMPLES::

            sage: Tableau([[1, 4, 6], [2, 5], [3]]).corners()
            [(0, 2), (1, 1), (2, 0)]
            sage: Tableau([[1, 3], [2, 4]]).corners()
            [(1, 1)]
        """
        return self.shape().corners()

    @combinatorial_map(order=2, name='conjugate')
    def conjugate(self):
        """
        Return the conjugate of ``self``.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).conjugate()
            [[1, 3], [2, 4]]
            sage: c = StandardTableau([[1,2],[3,4]]).conjugate()
            sage: c.parent()
            Standard tableaux
        """
        if self:
            conj = [[] for _ in repeat(None, len(self[0]))]
            for row in self:
                for j, x in enumerate(row):
                    conj[j].append(x)
        else:
            conj = []

        if isinstance(self, StandardTableau):
            return StandardTableau(conj)
        return Tableau(conj)

    def pp(self):
        """
        Pretty print a string of the tableau.

        EXAMPLES::

            sage: T = Tableau([[1,2,3],[3,4],[5]])
            sage: T.pp()
              1  2  3
              3  4
              5
            sage: Tableaux.options.convention="french"
            sage: T.pp()
              5
              3  4
              1  2  3
            sage: Tableaux.options.convention="russian"
            sage: T.pp()
             5    4    3
                3    2
                   1
            sage: Tableaux.options._reset()
        """
        print(self._repr_diagram())

    def plot(self, descents=False):
        r"""
        Return a plot ``self``.

        If English notation is set then the first row of the tableau is on the
        top:

        .. PLOT::
            :width: 200 px

            t = Tableau([[1,2,3,4],[2,3],[5]])
            Tableaux.options.convention="english"
            sphinx_plot(t.plot())

        If French notation is set, the first row of the tableau is on
        the bottom:

        .. PLOT::
            :width: 200 px

            t = Tableau([[1,2,3,4],[2,3],[5]])
            Tableaux.options.convention="french"
            sphinx_plot(t.plot())
            Tableaux.options.convention="english"

        If Russian notation is set, we tilt the French notation by 45 degrees:

        .. PLOT::
            :width: 200 px

            t = Tableau([[1,2,3,4],[2,3],[5]])
            Tableaux.options.convention="russian"
            sphinx_plot(t.plot())
            Tableaux.options.convention="english"

        INPUT:

        - ``descents`` -- boolean (default: ``False``); if ``True``,
          then the descents are marked in the tableau; only valid if
          ``self`` is a standard tableau

        EXAMPLES::

            sage: t = Tableau([[1,2,4],[3]])
            sage: t.plot()                                                              # needs sage.plot
            Graphics object consisting of 11 graphics primitives
            sage: t.plot(descents=True)                                                 # needs sage.plot
            Graphics object consisting of 12 graphics primitives

            sage: t = Tableau([[2,2,4],[3]])
            sage: t.plot()                                                              # needs sage.plot
            Graphics object consisting of 11 graphics primitives
            sage: t.plot(descents=True)                                                 # needs sage.plot
            Traceback (most recent call last):
            ...
            ValueError: the tableau must be standard for 'descents=True'
        """
        from sage.plot.polygon import polygon
        from sage.plot.line import line
        from sage.plot.text import text

        if descents and not self.is_standard():
            raise ValueError("the tableau must be standard for 'descents=True'")

        # For English we build up to down for French, down to up
        if self.parent().options('convention') == "English":
            m = 1
        else:
            m = -1

        p = self.shape()
        r = p.conjugate()

        # Increases the height
        if self.parent().options('convention') == "Russian":
            pp = p  # h
            rr = r
            h = [-i-1 for i in range(len(p))]
            v = [i+1 for i in range(len(r))]

        else:
            pp = [0] * len(p)
            rr = [0] * len(r)
            h = [0] * len(p)
            v = [0] * len(r)

        G = line([(0, 0), (p[0], pp[0])], axes=False, figsize=1.5)
        for i in range(len(p)):
            G += line([(h[i], m*(-i-1)), (h[i]+p[i], pp[i]+m*(-i-1))])

        G += line([(0, 0), (-rr[0], m*-r[0])])
        for i in range(len(r)):
            G += line([(i+1, v[i]), (i+1-rr[i], v[i]+m*-r[i])])

        if descents:
            t = StandardTableau(self)
            for i in t.standard_descents():
                c = t.cells_containing(i)[0]
                if self.parent().options('convention') == "Russian":
                    G += polygon([(c[1]+1-v[c[0]], m*(-c[1]-c[0])),
                                  (c[1]+2-v[c[0]], m*(-c[1]-c[0]-1)),
                                  (c[1]+1-v[c[0]], m*(-c[1]-c[0]-2)),
                                  (c[1]-v[c[0]], m*(-c[1]-c[0]-1))
                                  ],
                                 rgbcolor=(1, 0, 1)
                                 )
                else:
                    G += polygon([(c[1], m*-c[0]),
                                  (c[1]+1, m*-c[0]),
                                  (c[1]+1, m*(-c[0]-1)),
                                  (c[1], m*(-c[0]-1))
                                  ],
                                 rgbcolor=(1, 0, 1))

        if self.parent().options('convention') == "Russian":
            for c in self.cells():
                G += text(str(self.entry(c)), (c[1]+1-v[c[0]], m*(-c[1]-c[0]-1)))
        else:
            for c in self.cells():
                G += text(str(self.entry(c)), (c[1]+0.5, m*(-c[0]-0.5)))

        return G

    def to_word_by_row(self):
        """
        Return the word obtained from a row reading of the tableau ``self``
        (starting with the lowermost row, reading every row from left
        to right).

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).to_word_by_row()
            word: 3412
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word_by_row()
            word: 325146
        """
        from sage.combinat.words.word import Word
        w = []
        for row in reversed(self):
            w += row
        return Word(w)

    def to_word_by_column(self):
        """
        Return the word obtained from a column reading of the tableau ``self``
        (starting with the leftmost column, reading every column from bottom
        to top).

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).to_word_by_column()
            word: 3142
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word_by_column()
            word: 321546
        """
        from sage.combinat.words.word import Word
        w = []
        for row in self.conjugate():
            w += row[::-1]
        return Word(w)

    def to_word(self):
        """
        An alias for :meth:`to_word_by_row`.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).to_word()
            word: 3412
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word()
            word: 325146
        """
        return self.to_word_by_row()

    def descents(self):
        """
        Return a list of the cells ``(i,j)`` such that
        ``self[i][j] > self[i-1][j]``.

        .. WARNING::

            This is not to be confused with the descents of a standard tableau.

        EXAMPLES::

            sage: Tableau( [[1,4],[2,3]] ).descents()
            [(1, 0)]
            sage: Tableau( [[1,2],[3,4]] ).descents()
            [(1, 0), (1, 1)]
            sage: Tableau( [[1,2,3],[4,5]] ).descents()
            [(1, 0), (1, 1)]
        """
        return [(i, j)
                for i in range(1, len(self))
                for j, selfij in enumerate(self[i])
                if selfij > self[i-1][j]]

    def major_index(self):
        """
        Return the major index of ``self``.

        The major index of a tableau `T` is defined to be the sum of the number
        of descents of ``T`` (defined in :meth:`descents`) with the sum of
        their legs' lengths.

        .. WARNING::

            This is not to be confused with the major index of a
            standard tableau.

        EXAMPLES::

            sage: Tableau( [[1,4],[2,3]] ).major_index()
            1
            sage: Tableau( [[1,2],[3,4]] ).major_index()
            2

        If the major index would be defined in the sense of standard tableaux
        theory, then the following would give 3 for a result::

            sage: Tableau( [[1,2,3],[4,5]] ).major_index()
            2
        """
        descents = self.descents()
        p = self.shape()
        return len(descents) + sum([p.leg_length(*d) for d in descents])

    def inversions(self):
        """
        Return a list of the inversions of ``self``.

        Let `T` be a tableau. An inversion is an attacking pair `(c,d)` of
        the shape of `T` (see
        :meth:`~sage.combinat.partition.Partition.attacking_pairs` for
        a definition of this) such that the entry of `c` in `T` is
        greater than the entry of `d`.

        .. WARNING::

            Do not mistake this for the inversions of a standard tableau.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[2,5]])
            sage: t.inversions()
            [((1, 1), (0, 0))]
            sage: t = Tableau([[1,4,3],[5,2],[2,6],[3]])
            sage: t.inversions()
            [((0, 1), (0, 2)), ((1, 0), (1, 1)), ((1, 1), (0, 0)), ((2, 1), (1, 0))]
        """
        inversions = []
        previous_row = None
        for i, row in enumerate(self):
            for j, entry in enumerate(row):
                # c is in position (i,j)
                # find the d that satisfy condition 1
                inversions.extend(((i, j), (i, k))
                                  for k in range(j + 1, len(row))
                                  if entry > row[k])
                # find the d that satisfy condition 2
                if i == 0:
                    continue
                inversions.extend(((i, j), (i - 1, k))
                                  for k in range(j)
                                  if entry > previous_row[k])
            previous_row = row
        return inversions

    def inversion_number(self):
        """
        Return the inversion number of ``self``.

        The inversion number is defined to be the number of inversions of
        ``self`` minus the sum of the arm lengths of the descents of ``self``
        (see the :meth:`inversions` and :meth:`descents` methods for the
        relevant definitions).

        .. WARNING::

            This has none of the meanings in which the word "inversion"
            is used in the theory of standard tableaux.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[2,5]])
            sage: t.inversion_number()
            0
            sage: t = Tableau([[1,2,4],[3,5]])
            sage: t.inversion_number()
            0
        """
        p = self.shape()
        return len(self.inversions()) - sum(p.arm_length(*cell)
                                            for cell in self.descents())

    def to_sign_matrix(self, max_entry=None):
        r"""
        Return the sign matrix of ``self``.

        A sign matrix is an `m \times n` matrix of 0s, 1s, and -1s such that the
        partial sums of each column is either 0 or 1 and the partial sums of
        each row is nonnegative. [Ava2007]_

        INPUT:

        - ``max_entry`` -- nonnegative integer; the maximum allowable number in
          the tableau (defaults to the largest entry in the tableau if not specified)

        EXAMPLES::

            sage: t = SemistandardTableau([[1,1,1,2,4],[3,3,4],[4,5],[6,6]])
            sage: t.to_sign_matrix(6)                                                   # needs sage.modules
            [ 0  0  0  1  0  0]
            [ 0  1  0 -1  0  0]
            [ 1 -1  0  1  0  0]
            [ 0  0  1 -1  1  1]
            [ 0  0  0  1 -1  0]
            sage: t = Tableau([[1,2,4],[3,5]])
            sage: t.to_sign_matrix(7)                                                   # needs sage.modules
            [ 0  0  0  1  0  0  0]
            [ 0  1  0 -1  1  0  0]
            [ 1 -1  1  0 -1  0  0]
            sage: t = Tableau([(4,5,4,3),(2,1,3)])
            sage: t.to_sign_matrix(5)                                                   # needs sage.modules
            [ 0  0  1  0  0]
            [ 0  0  0  1  0]
            [ 1  0 -1 -1  1]
            [-1  1  0  1 -1]
            sage: s = Tableau([(1,0,-2,4),(3,4,5)])
            sage: s.to_sign_matrix(6)
            Traceback (most recent call last):
            ...
            ValueError: the entries must be nonnegative integers
        """
        from sage.rings.integer_ring import ZZ
        from sage.sets.positive_integers import PositiveIntegers
        PI = PositiveIntegers()
        for row in self:
            if any(c not in PI for c in row):
                raise ValueError("the entries must be nonnegative integers")
        from sage.matrix.matrix_space import MatrixSpace
        if max_entry is None:
            max_entry = max([max(c) for c in self])
        MS = MatrixSpace(ZZ, len(self[0]), max_entry)
        Tconj = self.conjugate()
        conj_len = len(Tconj)
        d = {(conj_len-i-1, elem-1): 1 for i, row in enumerate(Tconj) for elem in row}
        partial_sum_matrix = MS(d)
        from copy import copy
        sign_matrix = copy(MS.zero())
        for j in range(max_entry):
            sign_matrix[0, j] = partial_sum_matrix[0, j]
        for i in range(1, conj_len):
            for j in range(max_entry):
                sign_matrix[i, j] = partial_sum_matrix[i, j] - partial_sum_matrix[i-1, j]
        return sign_matrix

    def schuetzenberger_involution(self, n=None, check=True):
        r"""
        Return the Schuetzenberger involution of the tableau ``self``.

        This method relies on the analogous method on words, which reverts the
        word and then complements all letters within the underlying ordered
        alphabet. If `n` is specified, the underlying alphabet is assumed to
        be `[1, 2, \ldots, n]`. If no alphabet is specified, `n` is the maximal
        letter appearing in ``self``.

        INPUT:

        - ``n`` -- integer specifying the maximal letter in the
          alphabet (optional)
        - ``check`` -- boolean (default: ``True``); check to make sure ``self``
          is semistandard

        OUTPUT: a tableau, the Schuetzenberger involution of ``self``

        EXAMPLES::

            sage: t = Tableau([[1,1,1],[2,2]])
            sage: t.schuetzenberger_involution(3)
            [[2, 2, 3], [3, 3]]

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: t.schuetzenberger_involution()
            [[1, 2, 5], [3, 4]]

            sage: t = Tableau([[1,3,5,7],[2,4,6],[8,9]])
            sage: t.schuetzenberger_involution()
            [[1, 2, 6, 8], [3, 4, 9], [5, 7]]

            sage: t = Tableau([])
            sage: t.schuetzenberger_involution()
            []

            sage: t = StandardTableau([[1,2,3],[4,5]])
            sage: s = t.schuetzenberger_involution()
            sage: s.parent()
            Standard tableaux
        """
        if check and self not in SemistandardTableaux():
            raise ValueError("the tableau must be semistandard")
        w = [i for row in self for i in reversed(row)]
        # ``w`` is now the Semitic reading word of ``self`` (that is,
        # the reverse of the reading word of ``self``).
        if not w:
            return self
        if n is None:
            n = max(w)
        N = n + 1
        wi = [N - i for i in w]
        t = Tableau([[wi[0]]])
        for k in wi[1:]:
            t = t.bump(k)
        if isinstance(self, StandardTableau):
            return StandardTableau(list(t))
        if isinstance(self, SemistandardTableau):
            return SemistandardTableau(list(t))
        return t

    @combinatorial_map(order=2, name='evacuation')
    def evacuation(self, n=None, check=True):
        r"""
        Return the evacuation of the tableau ``self``.

        This is an alias for :meth:`schuetzenberger_involution`.

        This method relies on the analogous method on words, which reverts the
        word and then complements all letters within the underlying ordered
        alphabet. If `n` is specified, the underlying alphabet is assumed to
        be `[1, 2, \ldots, n]`. If no alphabet is specified, `n` is the maximal
        letter appearing in ``self``.

        INPUT:

        - ``n`` -- integer specifying the maximal letter in the
          alphabet (optional)
        - ``check`` -- boolean (default: ``True``); check to make sure ``self``
          is semistandard

        OUTPUT: a tableau, the evacuation of ``self``

        EXAMPLES::

            sage: t = Tableau([[1,1,1],[2,2]])
            sage: t.evacuation(3)
            [[2, 2, 3], [3, 3]]

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: t.evacuation()
            [[1, 2, 5], [3, 4]]

            sage: t = Tableau([[1,3,5,7],[2,4,6],[8,9]])
            sage: t.evacuation()
            [[1, 2, 6, 8], [3, 4, 9], [5, 7]]

            sage: t = Tableau([])
            sage: t.evacuation()
            []

            sage: t = StandardTableau([[1,2,3],[4,5]])
            sage: s = t.evacuation()
            sage: s.parent()
            Standard tableaux
        """
        return self.schuetzenberger_involution(n, check)

    @combinatorial_map(name='standardization')
    def standardization(self, check=True):
        r"""
        Return the standardization of ``self``, assuming ``self`` is a
        semistandard tableau.

        The standardization of a semistandard tableau `T` is the standard
        tableau `\mathrm{st}(T)` of the same shape as `T` whose
        reversed reading word is the standardization of the reversed reading
        word of `T`.

        The standardization of a word `w` can be formed by replacing all `1`'s in
        `w` by `1, 2, \ldots, k_1` from left to right, all `2`'s in `w` by
        `k_1 + 1, k_1 + 2, \ldots, k_2`, and repeating for all letters which
        appear in `w`.
        See also :meth:`Word.standard_permutation()`.

        INPUT:

        - ``check`` -- (default: ``True``) check to make sure ``self`` is
          semistandard; set to ``False`` to avoid this check

        EXAMPLES::

            sage: t = Tableau([[1,3,3,4],[2,4,4],[5,16]])
            sage: t.standardization()
            [[1, 3, 4, 7], [2, 5, 6], [8, 9]]

        Standard tableaux are fixed under standardization::

            sage: all((t == t.standardization() for t in StandardTableaux(6)))
            True
            sage: t = Tableau([])
            sage: t.standardization()
            []

        The reading word of the standardization is the standardization of
        the reading word::

            sage: T = SemistandardTableaux(shape=[5,2,2,1], max_entry=4)
            sage: all(t.to_word().standard_permutation() == t.standardization().reading_word_permutation() for t in T) # long time
            True
        """
        if check and self not in SemistandardTableaux():
            raise ValueError("the tableau must be semistandard")
        T = from_shape_and_word(self.shape(), self.to_word_by_row().standard_permutation())
        return StandardTableaux()(T)

    def bender_knuth_involution(self, k, rows=None, check=True):
        r"""
        Return the image of ``self`` under the `k`-th Bender--Knuth
        involution, assuming ``self`` is a semistandard tableau.

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
        called the *`k`-th Bender-Knuth involution*. This is used to show that
        the Schur functions defined by semistandard tableaux are symmetric
        functions.

        INPUT:

        - ``k`` -- integer

        - ``rows`` -- (default: ``None``) when set to ``None``, the method
          computes the `k`-th Bender--Knuth involution as defined above.
          When an iterable, this computes the composition of the `k`-th
          Bender--Knuth switches at row `i` over all `i` in ``rows``. When set
          to an integer `i`, the method computes the `k`-th Bender--Knuth
          switch at row `i`. Note the indexing of the rows starts with `1`.

        - ``check`` -- (default: ``True``) check to make sure ``self`` is
          semistandard; Set to ``False`` to avoid this check

        OUTPUT:

        The image of ``self`` under either the `k`-th Bender--Knuth
        involution, the `k`-th Bender--Knuth switch at a certain row, or
        the composition of such switches, as detailed in the INPUT section.

        EXAMPLES::

            sage: t = Tableau([[1,1,3,4,4,5,6,7],[2,2,4,6,7,7,7],[3,4,5,8,8,9],[6,6,7,10],[7,8,8,11],[8]])
            sage: t.bender_knuth_involution(1) == t
            True
            sage: t.bender_knuth_involution(2)
            [[1, 1, 2, 4, 4, 5, 6, 7], [2, 3, 4, 6, 7, 7, 7], [3, 4, 5, 8, 8, 9], [6, 6, 7, 10], [7, 8, 8, 11], [8]]
            sage: t.bender_knuth_involution(3)
            [[1, 1, 3, 3, 3, 5, 6, 7], [2, 2, 4, 6, 7, 7, 7], [3, 4, 5, 8, 8, 9], [6, 6, 7, 10], [7, 8, 8, 11], [8]]
            sage: t.bender_knuth_involution(4)
            [[1, 1, 3, 4, 5, 5, 6, 7], [2, 2, 4, 6, 7, 7, 7], [3, 5, 5, 8, 8, 9], [6, 6, 7, 10], [7, 8, 8, 11], [8]]
            sage: t.bender_knuth_involution(5)
            [[1, 1, 3, 4, 4, 5, 6, 7], [2, 2, 4, 5, 7, 7, 7], [3, 4, 6, 8, 8, 9], [5, 5, 7, 10], [7, 8, 8, 11], [8]]
            sage: t.bender_knuth_involution(666) == t
            True
            sage: t.bender_knuth_involution(4, 2) == t
            True
            sage: t.bender_knuth_involution(4, 3)
            [[1, 1, 3, 4, 4, 5, 6, 7], [2, 2, 4, 6, 7, 7, 7], [3, 5, 5, 8, 8, 9], [6, 6, 7, 10], [7, 8, 8, 11], [8]]

        The ``rows`` keyword can be an iterator::

            sage: t.bender_knuth_involution(6, iter([1,2])) == t
            False
            sage: t.bender_knuth_involution(6, iter([3,4])) == t
            True

        The Bender--Knuth involution is an involution::

            sage: T = SemistandardTableaux(shape=[3,1,1], max_entry=4)
            sage: all(t.bender_knuth_involution(k).bender_knuth_involution(k) == t      # needs sage.modules
            ....:     for k in range(1, 5) for t in T)
            True

        The same holds for the single switches::

            sage: all(t.bender_knuth_involution(k, j).bender_knuth_involution(k, j) == t            # needs sage.modules
            ....:     for k in range(1, 5) for j in range(1, 5) for t in T)
            True

        Locality of the Bender--Knuth involutions::

            sage: all(t.bender_knuth_involution(k).bender_knuth_involution(l)           # needs sage.modules
            ....:       == t.bender_knuth_involution(l).bender_knuth_involution(k)
            ....:     for k in range(1, 5) for l in range(1, 5) if abs(k - l) > 1
            ....:     for t in T)
            True

        Berenstein and Kirillov [KB1995]_ have shown that
        `(s_1 s_2)^6 = id` (for tableaux of straight shape)::

            sage: p = lambda t, k: t.bender_knuth_involution(k).bender_knuth_involution(k + 1)
            sage: all(p(p(p(p(p(p(t,1),1),1),1),1),1) == t for t in T)                  # needs sage.modules
            True

        However, `(s_2 s_3)^6 = id` is false::

            sage: p = lambda t, k: t.bender_knuth_involution(k).bender_knuth_involution(k + 1)
            sage: t = Tableau([[1,2,2],[3,4]])
            sage: x = t
            sage: for i in range(6): x = p(x, 2)
            sage: x
            [[1, 2, 3], [2, 4]]
            sage: x == t
            False

        TESTS::

            sage: t = Tableau([])
            sage: t.bender_knuth_involution(3)
            []
        """
        if check and self not in SemistandardTableaux():
            raise ValueError("the tableau must be semistandard")
        from sage.combinat.skew_tableau import SkewTableau
        sk = SkewTableau(self).bender_knuth_involution(k, rows, False)
        return SemistandardTableaux()(list(sk))

    @combinatorial_map(name='reading word permutation')
    def reading_word_permutation(self):
        """
        Return the permutation obtained by reading the entries of the
        standardization of ``self`` row by row, starting with the
        bottommost row (in English notation).

        EXAMPLES::

            sage: StandardTableau([[1,2],[3,4]]).reading_word_permutation()
            [3, 4, 1, 2]

        Check that :issue:`14724` is fixed::

            sage: SemistandardTableau([[1,1]]).reading_word_permutation()
            [1, 2]
        """
        return permutation.Permutation(self.standardization().to_word())

    def entries(self):
        """
        Return the tuple of all entries of ``self``, in the order obtained
        by reading across the rows from top to bottom (in English
        notation).

        EXAMPLES::

            sage: t = Tableau([[1,3], [2]])
            sage: t.entries()
            (1, 3, 2)
        """
        return sum(self, ())

    def entry(self, cell):
        """
        Return the entry of cell ``cell`` in the tableau ``self``. Here,
        ``cell`` should be given as a tuple `(i,j)` of zero-based
        coordinates (so the northwesternmost cell in English notation
        is `(0,0)`).

        EXAMPLES::

            sage: t = Tableau([[1,2],[3,4]])
            sage: t.entry( (0,0) )
            1
            sage: t.entry( (1,1) )
            4
        """
        i, j = cell
        return self[i][j]

    def weight(self):
        r"""
        Return the weight of the tableau ``self``. Trailing zeroes are
        omitted when returning the weight.

        The weight of a tableau `T` is the sequence `(a_1, a_2, a_3, \ldots )`,
        where `a_k` is the number of entries of `T` equal to `k`. This
        sequence contains only finitely many nonzero entries.

        The weight of a tableau `T` is the same as the weight of the
        reading word of `T`, for any reading order.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).weight()
            [1, 1, 1, 1]

            sage: Tableau([]).weight()
            []

            sage: Tableau([[1,3,3,7],[4,2],[2,3]]).weight()
            [1, 2, 3, 1, 0, 0, 1]

        TESTS:

        We check that this agrees with going to the word::

            sage: t = Tableau([[1,3,4,7],[6,2],[2,3]])
            sage: def by_word(T):
            ....:     ed = T.to_word().evaluation_dict()
            ....:     m = max(ed) + 1
            ....:     return [ed.get(k, 0) for k in range(1, m)]
            sage: by_word(t) == t.weight()
            True
            sage: SST = SemistandardTableaux(shape=[3,1,1])
            sage: all(by_word(t) == t.weight() for t in SST)                            # needs sage.modules
            True
        """
        if len(self) == 0:
            return []
        m = max(max(row) for row in self)
        res = [0] * m
        for row in self:
            for i in row:
                if i > 0:
                    res[i - 1] += 1
        return res

    evaluation = weight

    def is_row_strict(self):
        """
        Return ``True`` if ``self`` is a row strict tableau and ``False``
        otherwise.

        A tableau is row strict if the entries in each row are in
        (strictly) increasing order.

        EXAMPLES::

            sage: Tableau([[1, 3], [2, 4]]).is_row_strict()
            True
            sage: Tableau([[1, 2], [2, 4]]).is_row_strict()
            True
            sage: Tableau([[2, 3], [2, 4]]).is_row_strict()
            True
            sage: Tableau([[5, 3], [2, 4]]).is_row_strict()
            False
        """
        return all(row[i] < row[i+1] for row in self for i in range(len(row)-1))

    def is_row_increasing(self, weak=False):
        r"""
        Return ``True`` if the entries in each row are in increasing order,
        and ``False`` otherwise.

        By default, this checks for strictly increasing rows. Set ``weak``
        to ``True`` to test for weakly increasing rows.

        EXAMPLES::

            sage: T = Tableau([[1, 1, 3], [1, 2]])
            sage: T.is_row_increasing(weak=True)
            True
            sage: T.is_row_increasing()
            False
            sage: Tableau([[2, 1]]).is_row_increasing(weak=True)
            False
        """
        if weak:
            def test(a, b):
                return a <= b
        else:
            def test(a, b):
                return a < b
        return all(test(a, b) for row in self for (a, b) in zip(row, row[1:]))

    def is_column_increasing(self, weak=False):
        r"""
        Return ``True`` if the entries in each column are in increasing order,
        and ``False`` otherwise.

        By default, this checks for strictly increasing columns. Set ``weak``
        to ``True`` to test for weakly increasing columns.

        EXAMPLES::

            sage: T = Tableau([[1, 1, 3], [1, 2]])
            sage: T.is_column_increasing(weak=True)
            True
            sage: T.is_column_increasing()
            False
            sage: Tableau([[2], [1]]).is_column_increasing(weak=True)
            False
        """
        if weak:
            def test(a, b):
                return a <= b
        else:
            def test(a, b):
                return a < b

        def tworow(a, b):
            return all(test(a[i], b_i) for i, b_i in enumerate(b))
        return all(tworow(self[r], self[r + 1]) for r in range(len(self) - 1))

    def is_column_strict(self):
        """
        Return ``True`` if ``self`` is a column strict tableau and ``False``
        otherwise.

        A tableau is column strict if the entries in each column are in
        (strictly) increasing order.

        EXAMPLES::

            sage: Tableau([[1, 3], [2, 4]]).is_column_strict()
            True
            sage: Tableau([[1, 2], [2, 4]]).is_column_strict()
            True
            sage: Tableau([[2, 3], [2, 4]]).is_column_strict()
            False
            sage: Tableau([[5, 3], [2, 4]]).is_column_strict()
            False
            sage: Tableau([]).is_column_strict()
            True
            sage: Tableau([[1, 4, 2]]).is_column_strict()
            True
            sage: Tableau([[1, 4, 2], [2, 5]]).is_column_strict()
            True
            sage: Tableau([[1, 4, 2], [2, 3]]).is_column_strict()
            False
        """
        def tworow(a, b):
            return all(a[i] < b_i for i, b_i in enumerate(b))
        return all(tworow(self[r], self[r+1]) for r in range(len(self)-1))

    def is_semistandard(self):
        r"""
        Return ``True`` if ``self`` is a semistandard tableau, and ``False``
        otherwise.

        A tableau is semistandard if its rows weakly increase and its columns
        strictly increase.

        EXAMPLES::

            sage: Tableau([[1,1],[1,2]]).is_semistandard()
            False
            sage: Tableau([[1,2],[1,2]]).is_semistandard()
            False
            sage: Tableau([[1,1],[2,2]]).is_semistandard()
            True
            sage: Tableau([[1,2],[2,3]]).is_semistandard()
            True
            sage: Tableau([[4,1],[3,2]]).is_semistandard()
            False
        """
        return self.is_row_increasing(weak=True) and self.is_column_increasing()

    def is_standard(self):
        """
        Return ``True`` if ``self`` is a standard tableau and ``False``
        otherwise.

        EXAMPLES::

            sage: Tableau([[1, 3], [2, 4]]).is_standard()
            True
            sage: Tableau([[1, 2], [2, 4]]).is_standard()
            False
            sage: Tableau([[2, 3], [2, 4]]).is_standard()
            False
            sage: Tableau([[5, 3], [2, 4]]).is_standard()
            False
        """
        entries = sorted(self.entries())
        return entries == list(range(1, self.size() + 1)) and self.is_row_strict() and self.is_column_strict()

    def is_increasing(self):
        """
        Return ``True`` if ``self`` is an increasing tableau and
        ``False`` otherwise.

        A tableau is increasing if it is both row strict and column strict.

        EXAMPLES::

            sage: Tableau([[1, 3], [2, 4]]).is_increasing()
            True
            sage: Tableau([[1, 2], [2, 4]]).is_increasing()
            True
            sage: Tableau([[2, 3], [2, 4]]).is_increasing()
            False
            sage: Tableau([[5, 3], [2, 4]]).is_increasing()
            False
            sage: Tableau([[1, 2, 3], [2, 3], [3]]).is_increasing()
            True
        """
        return self.is_row_strict() and self.is_column_strict()

    def is_rectangular(self):
        """
        Return ``True`` if the tableau ``self`` is rectangular and
        ``False`` otherwise.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).is_rectangular()
            True
            sage: Tableau([[1,2,3],[4,5],[6]]).is_rectangular()
            False
            sage: Tableau([]).is_rectangular()
            True
        """
        if len(self) == 0:
            return True
        return len(self[-1]) == len(self[0])

    def vertical_flip(self):
        """
        Return the tableau obtained by vertically flipping the tableau
        ``self``.

        This only works for rectangular tableaux.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).vertical_flip()
            [[3, 4], [1, 2]]
        """
        if not self.is_rectangular():
            raise TypeError("the tableau must be rectangular to use vertical_flip()")

        return Tableau(list(reversed(self)))

    def rotate_180(self):
        """
        Return the tableau obtained by rotating ``self`` by `180` degrees.

        This only works for rectangular tableaux.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).rotate_180()
            [[4, 3], [2, 1]]
        """
        if not self.is_rectangular():
            raise TypeError("the tableau must be rectangular to use rotate_180()")

        return Tableau([list(reversed(row)) for row in reversed(self)])

    def cells(self):
        """
        Return a list of the coordinates of the cells of ``self``.

        Coordinates start at `0`, so the northwesternmost cell (in
        English notation) has coordinates `(0, 0)`.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).cells()
            [(0, 0), (0, 1), (1, 0), (1, 1)]
        """
        s = []
        for i, row in enumerate(self):
            s += [(i, j) for j in range(len(row))]
        return s

    def cells_containing(self, i):
        r"""
        Return the list of cells in which the letter `i` appears in the
        tableau ``self``. The list is ordered with cells appearing from
        left to right.

        Cells are given as pairs of coordinates `(a, b)`, where both
        rows and columns are counted from `0` (so `a = 0` means the cell
        lies in the leftmost column of the tableau, etc.).

        EXAMPLES::

            sage: t = Tableau([[1,1,3],[2,3,5],[4,5]])
            sage: t.cells_containing(5)
            [(2, 1), (1, 2)]
            sage: t.cells_containing(4)
            [(2, 0)]
            sage: t.cells_containing(6)
            []

            sage: t = Tableau([[1,1,2,4],[2,4,4],[4]])
            sage: t.cells_containing(4)
            [(2, 0), (1, 1), (1, 2), (0, 3)]

            sage: t = Tableau([[1,1,2,8,9],[2,5,6,11],[3,7,7,13],[4,8,9],[5],[13],[14]])
            sage: t.cells_containing(8)
            [(3, 1), (0, 3)]

            sage: Tableau([]).cells_containing(3)
            []
        """
        cell_list = []
        for r in range(len(self)-1, -1, -1):
            rth_row = self[r]
            for c, val in enumerate(rth_row):
                if val == i:
                    cell_list.append((r, c))
        return cell_list

    def leq(self, secondtab):
        """
        Check whether each entry of ``self`` is less-or-equal to the
        corresponding entry of a further tableau ``secondtab``.

        INPUT:

        - ``secondtab`` -- a tableau of the same shape as ``self``

        EXAMPLES::

            sage: T = Tableau([[1, 2], [3]])
            sage: S = Tableau([[1, 3], [3]])
            sage: G = Tableau([[2, 1], [4]])
            sage: H = Tableau([[1, 2], [4]])
            sage: T.leq(S)
            True
            sage: T.leq(T)
            True
            sage: T.leq(G)
            False
            sage: T.leq(H)
            True
            sage: S.leq(T)
            False
            sage: S.leq(G)
            False
            sage: S.leq(H)
            False
            sage: G.leq(H)
            False
            sage: H.leq(G)
            False

        TESTS::

            sage: StandardTableau(T).leq(S)
            True
            sage: T.leq(SemistandardTableau(S))
            True
        """
        if secondtab not in Tableaux():
            raise TypeError("{} must be a tableau".format(secondtab))
        sh = self.shape()
        if sh != secondtab.shape():
            raise TypeError("the tableaux must be the same shape")
        return all(self[a][b] <= secondtab[a][b]
                   for a in range(len(self))
                   for b in range(len(self[a])))

    def k_weight(self, k):
        r"""
        Return the `k`-weight of ``self``.

        A tableau has `k`-weight `\alpha = (\alpha_1, ..., \alpha_n)`
        if there are exactly `\alpha_i` distinct residues for the
        cells occupied by the letter `i` for each `i`.  The residue
        of a cell in position `(a,b)` is `a-b` modulo `k+1`.

        This definition is the one used in [Ive2012]_ (p. 12).

        EXAMPLES::

            sage: Tableau([[1,2],[2,3]]).k_weight(1)
            [1, 1, 1]
            sage: Tableau([[1,2],[2,3]]).k_weight(2)
            [1, 2, 1]
            sage: t = Tableau([[1,1,1,2,5],[2,3,6],[3],[4]])
            sage: t.k_weight(1)
            [2, 1, 1, 1, 1, 1]
            sage: t.k_weight(2)
            [3, 2, 2, 1, 1, 1]
            sage: t.k_weight(3)
            [3, 1, 2, 1, 1, 1]
            sage: t.k_weight(4)
            [3, 2, 2, 1, 1, 1]
            sage: t.k_weight(5)
            [3, 2, 2, 1, 1, 1]
        """
        res = []
        w = self.weight()
        s = self.cells()

        for l in range(1, len(w)+1):
            new_s = [(i, j) for i, j in s if self[i][j] == l]

            # If there are no elements that meet the condition
            if new_s == []:
                res.append(0)
                continue
            x = set((i-j) % (k+1) for i, j in new_s)
            res.append(len(x))

        return res

    def is_k_tableau(self, k):
        r"""
        Check whether ``self`` is a valid weak `k`-tableau.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[2,3],[3]])
            sage: t.is_k_tableau(3)
            True
            sage: t = Tableau([[1,1,3],[2,2],[3]])
            sage: t.is_k_tableau(3)
            False
        """
        shapes = self.to_chain()
        kshapes = [la.k_conjugate(k) for la in shapes]
        return all(kshapes[i+1].contains(kshapes[i]) for i in range(len(shapes)-1))

    def restrict(self, n):
        """
        Return the restriction of the semistandard tableau ``self``
        to ``n``. If possible, the restricted tableau will have the same
        parent as this tableau.

        If `T` is a semistandard tableau and `n` is a nonnegative integer,
        then the restriction of `T` to `n` is defined as the
        (semistandard) tableau obtained by removing all cells filled with
        entries greater than `n` from `T`.

        .. NOTE::

            If only the shape of the restriction, rather than the whole
            restriction, is needed, then the faster method
            :meth:`restriction_shape` is preferred.

        EXAMPLES::

            sage: Tableau([[1,2],[3],[4]]).restrict(3)
            [[1, 2], [3]]
            sage: StandardTableau([[1,2],[3],[4]]).restrict(2)
            [[1, 2]]
            sage: Tableau([[1,2,3],[2,4,4],[3]]).restrict(0)
            []
            sage: Tableau([[1,2,3],[2,4,4],[3]]).restrict(2)
            [[1, 2], [2]]
            sage: Tableau([[1,2,3],[2,4,4],[3]]).restrict(3)
            [[1, 2, 3], [2], [3]]
            sage: Tableau([[1,2,3],[2,4,4],[3]]).restrict(5)
            [[1, 2, 3], [2, 4, 4], [3]]

        If possible the restricted tableau will belong to the same category as
        the original tableau::

            sage: S = StandardTableau([[1,2,4,7],[3,5],[6]]); S.category()
            Category of elements of Standard tableaux
            sage: S.restrict(4).category()
            Category of elements of Standard tableaux
            sage: SS=StandardTableaux([4,2,1])([[1,2,4,7],[3,5],[6]]); SS.category()
            Category of elements of Standard tableaux of shape [4, 2, 1]
            sage: SS.restrict(4).category()
            Category of elements of Standard tableaux

            sage: Tableau([[1,2],[3],[4]]).restrict(3)
            [[1, 2], [3]]
            sage: Tableau([[1,2],[3],[4]]).restrict(2)
            [[1, 2]]
            sage: SemistandardTableau([[1,1],[2]]).restrict(1)
            [[1, 1]]
            sage: _.category()
            Category of elements of Semistandard tableaux
        """
        res = [[y for y in row if y <= n] for row in self]
        res = [row for row in res if row]
        # attempt to return a tableau of the same type
        try:
            return self.parent()(res)
        except Exception:
            try:
                return self.parent().Element(res)
            except Exception:
                return Tableau(res)

    def restriction_shape(self, n):
        """
        Return the shape of the restriction of the semistandard tableau
        ``self`` to ``n``.

        If `T` is a semistandard tableau and `n` is a nonnegative integer,
        then the restriction of `T` to `n` is defined as the
        (semistandard) tableau obtained by removing all cells filled with
        entries greater than `n` from `T`.

        This method computes merely the shape of the restriction. For
        the restriction itself, use :meth:`restrict`.

        EXAMPLES::

            sage: Tableau([[1,2],[2,3],[3,4]]).restriction_shape(3)
            [2, 2, 1]
            sage: StandardTableau([[1,2],[3],[4],[5]]).restriction_shape(2)
            [2]
            sage: Tableau([[1,3,3,5],[2,4,4],[17]]).restriction_shape(0)
            []
            sage: Tableau([[1,3,3,5],[2,4,4],[17]]).restriction_shape(2)
            [1, 1]
            sage: Tableau([[1,3,3,5],[2,4,4],[17]]).restriction_shape(3)
            [3, 1]
            sage: Tableau([[1,3,3,5],[2,4,4],[17]]).restriction_shape(5)
            [4, 3]

            sage: all( T.restriction_shape(i) == T.restrict(i).shape()
            ....:      for T in StandardTableaux(5) for i in range(1, 5) )
            True
        """
        from sage.combinat.partition import Partition
        res = [len([y for y in row if y <= n]) for row in self]
        return Partition(res)

    def to_chain(self, max_entry=None):
        r"""
        Return the chain of partitions corresponding to the (semi)standard
        tableau ``self``.

        The optional keyword parameter ``max_entry`` can be used to
        customize the length of the chain. Specifically, if this parameter
        is set to a nonnegative integer ``n``, then the chain is
        constructed from the positions of the letters `1, 2, \ldots, n`
        in the tableau.

        EXAMPLES::

            sage: Tableau([[1,2],[3],[4]]).to_chain()
            [[], [1], [2], [2, 1], [2, 1, 1]]
            sage: Tableau([[1,1],[2]]).to_chain()
            [[], [2], [2, 1]]
            sage: Tableau([[1,1],[3]]).to_chain()
            [[], [2], [2], [2, 1]]
            sage: Tableau([]).to_chain()
            [[]]
            sage: Tableau([[1,1],[2],[3]]).to_chain(max_entry=2)
            [[], [2], [2, 1]]
            sage: Tableau([[1,1],[2],[3]]).to_chain(max_entry=3)
            [[], [2], [2, 1], [2, 1, 1]]
            sage: Tableau([[1,1],[2],[3]]).to_chain(max_entry=4)
            [[], [2], [2, 1], [2, 1, 1], [2, 1, 1]]
            sage: Tableau([[1,1,2],[2,3],[4,5]]).to_chain(max_entry=6)
            [[], [2], [3, 1], [3, 2], [3, 2, 1], [3, 2, 2], [3, 2, 2]]
        """
        if max_entry is None:
            if len(self) == 0:
                max_entry = 0
            else:
                max_entry = max(max(row) for row in self)
        return [self.restriction_shape(k) for k in range(max_entry+1)]

    @combinatorial_map(name='to Gelfand-Tsetlin pattern')
    def to_Gelfand_Tsetlin_pattern(self):
        """
        Return the :class:`Gelfand-Tsetlin pattern <GelfandTsetlinPattern>`
        corresponding to ``self`` when semistandard.

        EXAMPLES::

            sage: T = Tableau([[1,2,3],[2,3],[3]])
            sage: G = T.to_Gelfand_Tsetlin_pattern(); G
            [[3, 2, 1], [2, 1], [1]]
            sage: G.to_tableau() == T
            True
            sage: T = Tableau([[1,3],[2]])
            sage: T.to_Gelfand_Tsetlin_pattern()
            [[2, 1, 0], [1, 1], [1]]
        """
        from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPatterns
        return GelfandTsetlinPatterns()(self)

    def anti_restrict(self, n):
        """
        Return the skew tableau formed by removing all of the cells from
        ``self`` that are filled with a number at most `n`.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[4,5]]); t
            [[1, 2, 3], [4, 5]]
            sage: t.anti_restrict(1)
            [[None, 2, 3], [4, 5]]
            sage: t.anti_restrict(2)
            [[None, None, 3], [4, 5]]
            sage: t.anti_restrict(3)
            [[None, None, None], [4, 5]]
            sage: t.anti_restrict(4)
            [[None, None, None], [None, 5]]
            sage: t.anti_restrict(5)
            [[None, None, None], [None, None]]
        """
        t_new = [[None if g <= n else g for g in row] for row in self]
        from sage.combinat.skew_tableau import SkewTableau
        return SkewTableau(t_new)

    def to_list(self):
        """
        Return ``self`` as a list of lists (not tuples!).

        EXAMPLES::

            sage: t = Tableau([[1,2],[3,4]])
            sage: l = t.to_list(); l
            [[1, 2], [3, 4]]
            sage: l[0][0] = 2
            sage: t
            [[1, 2], [3, 4]]
        """
        return [list(row) for row in self]

    def bump(self, x):
        """
        Insert ``x`` into ``self`` using Schensted's row-bumping (or
        row-insertion) algorithm.

        EXAMPLES::

            sage: t = Tableau([[1,2],[3]])
            sage: t.bump(1)
            [[1, 1], [2], [3]]
            sage: t
            [[1, 2], [3]]
            sage: t.bump(2)
            [[1, 2, 2], [3]]
            sage: t.bump(3)
            [[1, 2, 3], [3]]
            sage: t
            [[1, 2], [3]]
            sage: t = Tableau([[1,2,2,3],[2,3,5,5],[4,4,6],[5,6]])
            sage: t.bump(2)
            [[1, 2, 2, 2], [2, 3, 3, 5], [4, 4, 5], [5, 6, 6]]
            sage: t.bump(1)
            [[1, 1, 2, 3], [2, 2, 5, 5], [3, 4, 6], [4, 6], [5]]
        """
        to_insert = x
        new_t = self.to_list()
        for row in new_t:
            i = 0
            # try to insert to_insert into row
            while i < len(row):
                if to_insert < row[i]:
                    t = to_insert
                    to_insert = row[i]
                    row[i] = t
                    break
                i += 1

            # if we haven't already inserted to_insert
            # append it to the end of row
            if i == len(row):
                row.append(to_insert)
                if isinstance(self, SemistandardTableau):
                    return SemistandardTableau(new_t)
                return Tableau(new_t)
        # if we got here, we are at the end of the tableau
        # add to_insert as the last row
        new_t.append([to_insert])
        if isinstance(self, SemistandardTableau):
            return SemistandardTableau(new_t)
        return Tableau(new_t)

    def schensted_insert(self, i, left=False):
        """
        Insert ``i`` into ``self`` using Schensted's row-bumping (or
        row-insertion) algorithm.

        INPUT:

        - ``i`` -- a number to insert
        - ``left`` -- boolean (default: ``False``); if set to
          ``True``, the insertion will be done from the left. That
          is, if one thinks of the algorithm as appending a letter
          to the reading word of ``self``, we append the letter to
          the left instead of the right

        EXAMPLES::

            sage: t = Tableau([[3,5],[7]])
            sage: t.schensted_insert(8)
            [[3, 5, 8], [7]]
            sage: t.schensted_insert(8, left=True)
            [[3, 5], [7], [8]]
        """
        if left:
            return self._left_schensted_insert(i)
        else:
            return self.bump(i)

    def _left_schensted_insert(self, letter):
        """
        EXAMPLES::

            sage: t = Tableau([[3,5],[7]])
            sage: t._left_schensted_insert(8)
            [[3, 5], [7], [8]]
            sage: t._left_schensted_insert(6)
            [[3, 5], [6, 7]]
            sage: t._left_schensted_insert(2)
            [[2, 3, 5], [7]]
        """
        h = len(self)
        if h == 0:
            return Tableau([[letter]])
        h1 = h + 1
        rep = self.to_list()
        rep.reverse()

        width = len(rep[h-1])
        heights = self._heights() + [h1]

        for j in range(1, width+2):
            i = heights[j-1]
            while i != h1 and rep[i-1][j-1] >= letter:
                i += 1
            if i == heights[j-1]:  # add on top of column j
                if j == 1:
                    rep = [[letter]] + rep
                else:
                    rep[i-2].append(letter)
                break
            elif i == h1 and j == width:  # add on right of line i
                if rep[i-2][j-1] < letter:
                    rep[i-2].append(letter)
                else:
                    new_letter = rep[i-2][j-1]
                    rep[i-2][j-1] = letter
                    rep[i-2].append(new_letter)
                break
            else:
                new_letter = rep[i-2][j-1]
                rep[i-2][j-1] = letter
                letter = new_letter

        rep.reverse()
        return Tableau(rep)

    def insert_word(self, w, left=False):
        """
        Insert the word ``w`` into the tableau ``self`` letter by letter
        using Schensted insertion. By default, the word ``w`` is being
        processed from left to right, and the insertion used is row
        insertion. If the optional keyword ``left`` is set to ``True``,
        the word ``w`` is being processed from right to left, and column
        insertion is used instead.

        EXAMPLES::

            sage: t0 = Tableau([])
            sage: w = [1,1,2,3,3,3,3]
            sage: t0.insert_word(w)
            [[1, 1, 2, 3, 3, 3, 3]]
            sage: t0.insert_word(w,left=True)
            [[1, 1, 2, 3, 3, 3, 3]]
            sage: w.reverse()
            sage: t0.insert_word(w)
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t0.insert_word(w,left=True)
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t1 = Tableau([[1,3],[2]])
            sage: t1.insert_word([4,5])
            [[1, 3, 4, 5], [2]]
            sage: t1.insert_word([4,5], left=True)
            [[1, 3], [2, 5], [4]]
        """
        if left:
            w = list(reversed(w))
        res = self
        for i in w:
            res = res.schensted_insert(i, left=left)
        return res

    def reverse_bump(self, loc):
        r"""
        Reverse row bump the entry of ``self`` at the specified
        location ``loc`` (given as a row index or a
        corner ``(r, c)`` of the tableau).

        This is the reverse of Schensted's row-insertion algorithm.
        See Section 1.1, page 8, of Fulton's [Ful1997]_.

        INPUT:

        - ``loc`` -- can be either of the following:

          - The coordinates ``(r, c)`` of the square to reverse-bump
            (which must be a corner of the tableau);
          - The row index ``r`` of this square.

          Note that both ``r`` and ``c`` are `0`-based, i.e., the
          topmost row and the leftmost column are the `0`-th row
          and the `0`-th column.

        OUTPUT: an ordered pair consisting of:

        1. The resulting (smaller) tableau;
        2. The entry bumped out at the end of the process.

        .. SEEALSO::

            :meth:`bump`

        EXAMPLES:

        This is the reverse of Schensted's bump::

            sage: T = Tableau([[1, 1, 2, 2, 4], [2, 3, 3], [3, 4], [4]])
            sage: T.reverse_bump(2)
            ([[1, 1, 2, 3, 4], [2, 3, 4], [3], [4]], 2)
            sage: T == T.reverse_bump(2)[0].bump(2)
            True
            sage: T.reverse_bump((3, 0))
            ([[1, 2, 2, 2, 4], [3, 3, 3], [4, 4]], 1)

        Some errors caused by wrong input::

            sage: T.reverse_bump((3, 1))
            Traceback (most recent call last):
            ...
            ValueError: invalid corner
            sage: T.reverse_bump(4)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: Tableau([[2, 2, 1], [3, 3]]).reverse_bump(0)
            Traceback (most recent call last):
            ...
            ValueError: reverse bumping is only defined for semistandard tableaux

        Some edge cases::

            sage: Tableau([[1]]).reverse_bump(0)
            ([], 1)
            sage: Tableau([[1,1]]).reverse_bump(0)
            ([[1]], 1)
            sage: Tableau([]).reverse_bump(0)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range

        .. NOTE::

            Reverse row bumping is only implemented for tableaux with weakly increasing
            and strictly increasing columns (though the tableau does not need to be an
            instance of class :class:`SemistandardTableau`).
        """
        if not (self.is_semistandard()):
            raise ValueError("reverse bumping is only defined for semistandard tableaux")
        try:
            (r, c) = loc
            if (r, c) not in self.corners():
                raise ValueError("invalid corner")
        except TypeError:
            r = loc
            c = len(self[r]) - 1

        # make a copy of self
        new_t = self.to_list()

        # remove the last entry of row r from the tableau
        to_move = new_t[r].pop()

        # delete the row if it's now empty
        if not new_t[r]:
            new_t.pop()

        from bisect import bisect_left

        for row in reversed(new_t[:r]):
            # Decide where to insert:
            # the bisect_left command returns the greatest index such that
            # every entry to its left is strictly less than to_move
            c = bisect_left(row, to_move, lo=c) - 1

            # swap it with to_move
            row[c], to_move = to_move, row[c]

        if isinstance(self, SemistandardTableau):
            return SemistandardTableau(new_t), to_move
        return Tableau(new_t), to_move

    def bump_multiply(self, other):
        """
        Multiply two tableaux using Schensted's bump.

        This product makes the set of semistandard tableaux into an
        associative monoid. The empty tableau is the unit in this monoid.
        See pp. 11-12 of [Ful1997]_.

        The same product operation is implemented in a different way in
        :meth:`slide_multiply`.

        EXAMPLES::

            sage: t = Tableau([[1,2,2,3],[2,3,5,5],[4,4,6],[5,6]])
            sage: t2 = Tableau([[1,2],[3]])
            sage: t.bump_multiply(t2)
            [[1, 1, 2, 2, 3], [2, 2, 3, 5], [3, 4, 5], [4, 6, 6], [5]]
        """
        if not isinstance(other, Tableau):
            raise TypeError("other must be a Tableau")

        row = len(other)
        product = Tableau([list(a) for a in self])  # create deep copy of self
        while row:
            row -= 1
            for i in other[row]:
                product = product.bump(i)
        return product

    def slide_multiply(self, other):
        """
        Multiply two tableaux using jeu de taquin.

        This product makes the set of semistandard tableaux into an
        associative monoid. The empty tableau is the unit in this monoid.

        See pp. 15 of [Ful1997]_.

        The same product operation is implemented in a different way in
        :meth:`bump_multiply`.

        EXAMPLES::

            sage: t = Tableau([[1,2,2,3],[2,3,5,5],[4,4,6],[5,6]])
            sage: t2 = Tableau([[1,2],[3]])
            sage: t.slide_multiply(t2)
            [[1, 1, 2, 2, 3], [2, 2, 3, 5], [3, 4, 5], [4, 6, 6], [5]]
        """
        if len(self) == 0:
            return other

        l = len(self[0])
        st = [(None,) * l + row for row in other]
        st.extend(self)

        from sage.combinat.skew_tableau import SkewTableau
        return SkewTableau(st).rectify()

    def _slide_up(self, c):
        r"""
        Auxiliary method used for promotion, which removes cell `c` from ``self``,
        slides the letters of ``self`` up using jeu de taquin slides, and
        then fills the empty cell at `(0,0)` with the value `0`.

        TESTS::

            sage: t = Tableau([[1,1,2],[2,3,5],[4,5]])
            sage: t._slide_up((2,1))
            [[0, 1, 2], [1, 3, 5], [2, 4]]

            sage: t._slide_up((1,2))
            [[0, 1, 2], [1, 2, 3], [4, 5]]

            sage: t = Tableau([[1,1,3],[2,3,5],[4,5]])
            sage: t._slide_up((1,2))
            [[0, 1, 1], [2, 3, 3], [4, 5]]
        """
        new_st = self.to_list()
        spotl, spotc = c
        while [spotl, spotc] != [0, 0]:
            # once moving box is in first column, just move letters up
            # (French notation!)
            if spotc == 0:
                new_st[spotl][spotc] = new_st[spotl-1][spotc]
                spotl -= 1
                continue
            # once moving box is in first row, just move letters up
            elif spotl == 0:
                new_st[spotl][spotc] = new_st[spotl][spotc-1]
                spotc -= 1
                continue
            else:
                # If we get to this stage, we need to compare
                below = new_st[spotl-1][spotc]
                left = new_st[spotl][spotc-1]
                if below >= left:
                    # Swap with the cell below
                    new_st[spotl][spotc] = new_st[spotl-1][spotc]
                    spotl -= 1
                    continue
                else:
                    # Swap with the cell to the left
                    new_st[spotl][spotc] = new_st[spotl][spotc-1]
                    spotc -= 1
                    continue
        # set box in position (0,0) to 0
        new_st[0][0] = 0
        return Tableau(new_st)

    def _slide_down(self, c, n):
        r"""
        Auxiliary method used for promotion, which removes cell `c` from ``self``,
        slides the letters of ``self`` down using jeu de taquin slides, and
        then fills the empty cell with the value `n + 2`.

        When the entries of ``self`` are positive integers, and cell `c` is
        filled with `1`, then the position of `c` is irrelevant.

        TESTS::

            sage: t = Tableau([[1,1,2],[2,3,5],[4,5]])
            sage: t._slide_down((0, 0), 8)
            [[1, 2, 5], [2, 3, 10], [4, 5]]

            sage: t._slide_down((0, 1), 8)
            [[1, 2, 5], [2, 3, 10], [4, 5]]

            sage: t = Tableau([[1,1,2,2,2,3],[2,2,4,6,6],[4,4,5,7],[5,8]])
            sage: t._slide_down((0, 1), 9)
            [[1, 2, 2, 2, 2, 3], [2, 4, 4, 6, 6], [4, 5, 7, 11], [5, 8]]
        """
        new_st = self.to_list()
        # new_st is a deep copy of self, so as not to mess around with self.
        new_st_shape = [len(x) for x in self]
        spotl, spotc = c
        # spotl and spotc are the coordinates of the wandering hole.
        # All comments and variable names below refer to French notation.
        while True:
            # "right_neighbor" and "upper_neighbor" refer to neighbors of the
            # hole.
            go_right = None
            if len(new_st_shape) > spotl + 1 and new_st_shape[spotl + 1] >= spotc + 1:
                upper_neighbor = new_st[spotl + 1][spotc]
                go_right = False
            if new_st_shape[spotl] != spotc + 1:
                right_neighbor = new_st[spotl][spotc + 1]
                if go_right is None or upper_neighbor > right_neighbor:
                    go_right = True
            if go_right is True:
                new_st[spotl][spotc] = right_neighbor
                spotc += 1
            elif go_right is False:
                new_st[spotl][spotc] = upper_neighbor
                spotl += 1
            else:
                break
        new_st[spotl][spotc] = n + 2
        return Tableau(new_st)

    def promotion_inverse(self, n):
        r"""
        Return the image of ``self`` under the inverse promotion operator.

        .. WARNING::

            You might know this operator as the promotion operator
            (without "inverse") -- literature does not agree on the
            name.

        The inverse promotion operator, applied to a tableau `t`, does the
        following:

        Iterate over all letters `1` in the tableau `t`, from right to left.
        For each of these letters, do the following:

        - Remove the letter from `t`, thus leaving a hole where it used to be.

        - Apply jeu de taquin to move this hole northeast (in French notation)
          until it reaches the outer boundary of `t`.

        - Fill `n+2` into the hole once jeu de taquin has completed.

        Once this all is done, subtract `1` from each letter in the tableau.
        This is not always well-defined. Restricted to the class of
        semistandard tableaux whose entries are all `\leq n + 1`, this is the
        usual inverse promotion operator defined on this class.

        When ``self`` is a standard tableau of size ``n + 1``, this definition of
        inverse promotion is the map called "promotion" in [Sag2011]_ (p. 23) and
        in [Stan2009]_, and is the inverse of the map called "promotion" in
        [Hai1992]_ (p. 90).

        .. WARNING::

            To my (Darij's) knowledge, the fact that the above "inverse
            promotion operator" really is the inverse of the promotion
            operator :meth:`promotion` for semistandard tableaux has never
            been proven in literature. Corrections are welcome.

        EXAMPLES::

            sage: t = Tableau([[1,2],[3,3]])
            sage: t.promotion_inverse(2)
            [[1, 2], [2, 3]]

            sage: t = Tableau([[1,2],[2,3]])
            sage: t.promotion_inverse(2)
            [[1, 1], [2, 3]]

            sage: t = Tableau([[1,2,5],[3,3,6],[4,7]])
            sage: t.promotion_inverse(8)
            [[1, 2, 4], [2, 5, 9], [3, 6]]

            sage: t = Tableau([])
            sage: t.promotion_inverse(2)
            []

        TESTS:

        We check the equivalence of two definitions of inverse promotion
        on semistandard tableaux::

            sage: ST = SemistandardTableaux(shape=[4,2,1], max_entry=7)
            sage: def bk_promotion_inverse7(st):
            ....:     st2 = st
            ....:     for i in range(1, 7):
            ....:         st2 = st2.bender_knuth_involution(i, check=False)
            ....:     return st2
            sage: all( bk_promotion_inverse7(st) == st.promotion_inverse(6) for st in ST ) # long time
            True
            sage: ST = SemistandardTableaux(shape=[2,2,2], max_entry=7)
            sage: all( bk_promotion_inverse7(st) == st.promotion_inverse(6) for st in ST ) # long time
            True

        A test for :issue:`13203`::

            sage: T = Tableau([[1]])
            sage: type(T.promotion_inverse(2)[0][0])
            <class 'sage.rings.integer.Integer'>
        """
        if self.is_rectangular():
            n = Integer(n)
            if self.size() == 0:
                return self
            s = self.shape()[0]
            l = self.weight()[0]
            word = [i-1 for row in reversed(self) for i in row if i > 1]
            t = Tableau([])
            t = t.insert_word(word)
            t = t.to_list()
            if l < s:
                for i in range(l):
                    t[len(t)-1].append(n+1)
            else:
                t.append([n+1 for i in range(s)])
            return Tableau(t)
        # Now, the non-rectangular case.
        p = self
        for c in reversed(self.cells_containing(1)):
            p = p._slide_down(c, n)
        return Tableau([[i-1 for i in row] for row in p])

    def promotion(self, n):
        r"""
        Return the image of ``self`` under the promotion operator.

        .. WARNING::

            You might know this operator as the inverse promotion
            operator -- literature does not agree on the name. You
            might also be looking for the Lapointe-Lascoux-Morse
            promotion operator (:meth:`promotion_operator`).

        The promotion operator, applied to a tableau `t`, does the following:

        Iterate over all letters `n+1` in the tableau `t`, from left to right.
        For each of these letters, do the following:

        - Remove the letter from `t`, thus leaving a hole where it used to be.

        - Apply jeu de taquin to move this hole southwest (in French notation)
          until it reaches the inner boundary of `t`.

        - Fill `0` into the hole once jeu de taquin has completed.

        Once this all is done, add `1` to each letter in the tableau.
        This is not always well-defined. Restricted to the class of
        semistandard tableaux whose entries are all `\leq n + 1`, this is the
        usual promotion operator defined on this class.

        When ``self`` is a standard tableau of size ``n + 1``, this definition of
        promotion is precisely the one given in [Hai1992]_ (p. 90). It is the
        inverse of the maps called "promotion" in [Sag2011]_ (p. 23) and in [Stan2009]_.

        .. WARNING::

            To my (Darij's) knowledge, the fact that the above promotion
            operator really is the inverse of the "inverse promotion
            operator" :meth:`promotion_inverse` for semistandard tableaux
            has never been proven in literature. Corrections are welcome.

        REFERENCES:

        - [Hai1992]_

        - [Sag2011]_

        EXAMPLES::

            sage: t = Tableau([[1,2],[3,3]])
            sage: t.promotion(2)
            [[1, 1], [2, 3]]

            sage: t = Tableau([[1,1,1],[2,2,3],[3,4,4]])
            sage: t.promotion(3)
            [[1, 1, 2], [2, 2, 3], [3, 4, 4]]

            sage: t = Tableau([[1,2],[2]])
            sage: t.promotion(3)
            [[2, 3], [3]]

            sage: t = Tableau([[1,1,3],[2,2]])
            sage: t.promotion(2)
            [[1, 2, 2], [3, 3]]

            sage: t = Tableau([[1,1,3],[2,3]])
            sage: t.promotion(2)
            [[1, 1, 2], [2, 3]]

            sage: t = Tableau([])
            sage: t.promotion(2)
            []

        TESTS:

        We check the equivalence of two definitions of promotion on
        semistandard tableaux::

            sage: ST = SemistandardTableaux(shape=[3,2,2,1], max_entry=6)
            sage: def bk_promotion6(st):
            ....:     st2 = st
            ....:     for i in range(5, 0, -1):
            ....:         st2 = st2.bender_knuth_involution(i, check=False)
            ....:     return st2
            sage: all( bk_promotion6(st) == st.promotion(5) for st in ST ) # long time
            True
            sage: ST = SemistandardTableaux(shape=[4,4], max_entry=6)
            sage: all( bk_promotion6(st) == st.promotion(5) for st in ST ) # long time
            True

        We also check :meth:`promotion_inverse()` is the inverse
        of :meth:`promotion()`::

            sage: ST = SemistandardTableaux(shape=[3,2,1], max_entry=7)
            sage: all( st.promotion(6).promotion_inverse(6) == st for st in ST ) # long time
            True
        """
        if self.is_rectangular():
            t = self.rotate_180()
            t = [tuple(n+2-i for i in row) for row in t]
            t = Tableau(t).promotion_inverse(n)
            t = [tuple(n+2-i for i in row) for row in t]
            return Tableau(t).rotate_180()
        p = self
        for c in self.cells_containing(n+1):
            p = p._slide_up(c)
        return Tableau([[i+1 for i in row] for row in p])

    def row_stabilizer(self):
        """
        Return the PermutationGroup corresponding to the row stabilizer of
        ``self``.

        This assumes that every integer from `1` to the size of ``self``
        appears exactly once in ``self``.

        EXAMPLES::

            sage: # needs sage.groups
            sage: rs = Tableau([[1,2,3],[4,5]]).row_stabilizer()
            sage: rs.order() == factorial(3)*factorial(2)
            True
            sage: PermutationGroupElement([(1,3,2),(4,5)]) in rs
            True
            sage: PermutationGroupElement([(1,4)]) in rs
            False
            sage: rs = Tableau([[1, 2],[3]]).row_stabilizer()
            sage: PermutationGroupElement([(1,2),(3,)]) in rs
            True
            sage: rs.one().domain()
            [1, 2, 3]
            sage: rs = Tableau([[1],[2],[3]]).row_stabilizer()
            sage: rs.order()
            1
            sage: rs = Tableau([[2,4,5],[1,3]]).row_stabilizer()
            sage: rs.order()
            12
            sage: rs = Tableau([]).row_stabilizer()
            sage: rs.order()
            1
        """
        # Ensure that the permutations involve all elements of the
        # tableau, by including the identity permutation on the set [1..k].
        k = self.size()
        gens = [list(range(1, k + 1))]
        for row in self:
            for j in range(len(row) - 1):
                gens.append((row[j], row[j + 1]))
        return PermutationGroup(gens)

    def column_stabilizer(self):
        """
        Return the PermutationGroup corresponding to the column stabilizer
        of ``self``.

        This assumes that every integer from `1` to the size of ``self``
        appears exactly once in ``self``.

        EXAMPLES::

            sage: # needs sage.groups
            sage: cs = Tableau([[1,2,3],[4,5]]).column_stabilizer()
            sage: cs.order() == factorial(2)*factorial(2)
            True
            sage: PermutationGroupElement([(1,3,2),(4,5)]) in cs
            False
            sage: PermutationGroupElement([(1,4)]) in cs
            True
        """
        # Ensure that the permutations involve all elements of the
        # tableau, by including the identity permutation on the set [1..k].
        k = self.size()
        gens = [list(range(1, k + 1))]
        ell = len(self)
        while ell > 1:
            ell -= 1
            for i, val in enumerate(self[ell]):
                gens.append((val, self[ell-1][i]))
        return PermutationGroup(gens)

    def height(self):
        """
        Return the height of ``self``.

        EXAMPLES::

            sage: Tableau([[1,2,3],[4,5]]).height()
            2
            sage: Tableau([[1,2,3]]).height()
            1
            sage: Tableau([]).height()
            0
        """
        return len(self)

    def _heights(self):
        """
        EXAMPLES::

            sage: Tableau([[1,2,3,4],[5,6],[7],[8]])._heights()
            [1, 3, 4, 4]
            sage: Tableau([])._heights()
            []
            sage: Tableau([[1]])._heights()
            [1]
            sage: Tableau([[1,2]])._heights()
            [1, 1]
            sage: Tableau([[1,2],[3],[4]])._heights()
            [1, 3]
        """
        cor = self.corners()
        if not cor:
            return []
        ncor = len(cor)
        k = len(self)
        cor = [[k - i, j + 1] for i, j in reversed(cor)]

        heights = [1]*(cor[0][1])
        for i in range(1, ncor):
            heights += [cor[i][0]]*(cor[i][1]-cor[i-1][1])

        return heights

    def last_letter_lequal(self, tab2):
        """
        Return ``True`` if ``self`` is less than or equal to ``tab2`` in the last
        letter ordering.

        EXAMPLES::

            sage: st = StandardTableaux([3,2])
            sage: f = lambda b: 1 if b else 0
            sage: matrix([[f(t1.last_letter_lequal(t2)) for t2 in st] for t1 in st])    # needs sage.modules
            [1 1 1 1 1]
            [0 1 1 1 1]
            [0 0 1 1 1]
            [0 0 0 1 1]
            [0 0 0 0 1]
        """
        n = self.size()
        if not isinstance(tab2, Tableau):
            try:
                tab2 = Tableau(tab2)
            except Exception:
                raise TypeError("tab2 must be a standard tableau")

        if tab2.size() != n:
            raise ValueError("tab2 must be the same size as self")

        if self == tab2:
            return True

        for j in range(n, 1, -1):
            self_j_pos = None
            for i in range(len(self)):
                if j in self[i]:
                    self_j_pos = i
                    break

            tab2_j_pos = None
            for i in range(len(tab2)):
                if j in tab2[i]:
                    tab2_j_pos = i
                    break

            if self_j_pos < tab2_j_pos:
                return True
            if tab2_j_pos < self_j_pos:
                return False

    def charge(self):
        r"""
        Return the charge of the reading word of ``self``.  See
        :meth:`~sage.combinat.words.finite_word.FiniteWord_class.charge` for more information.

        EXAMPLES::

            sage: Tableau([[1,1],[2,2],[3]]).charge()
            0
            sage: Tableau([[1,1,3],[2,2]]).charge()
            1
            sage: Tableau([[1,1,2],[2],[3]]).charge()
            1
            sage: Tableau([[1,1,2],[2,3]]).charge()
            2
            sage: Tableau([[1,1,2,3],[2]]).charge()
            2
            sage: Tableau([[1,1,2,2],[3]]).charge()
            3
            sage: Tableau([[1,1,2,2,3]]).charge()
            4
        """
        return self.to_word().charge()

    def cocharge(self):
        r"""
        Return the cocharge of the reading word of ``self``.  See
        :meth:`~sage.combinat.words.finite_word.FiniteWord_class.cocharge` for more information.

        EXAMPLES::

            sage: Tableau([[1,1],[2,2],[3]]).cocharge()
            4
            sage: Tableau([[1,1,3],[2,2]]).cocharge()
            3
            sage: Tableau([[1,1,2],[2],[3]]).cocharge()
            3
            sage: Tableau([[1,1,2],[2,3]]).cocharge()
            2
            sage: Tableau([[1,1,2,3],[2]]).cocharge()
            2
            sage: Tableau([[1,1,2,2],[3]]).cocharge()
            1
            sage: Tableau([[1,1,2,2,3]]).cocharge()
            0
        """
        return self.to_word().cocharge()

    def add_entry(self, cell, m):
        """
        Return the result of setting the entry in cell ``cell`` equal
        to ``m`` in the tableau ``self``.

        This tableau has larger size than ``self`` if ``cell`` does not
        belong to the shape of ``self``; otherwise, the tableau has the
        same shape as ``self`` and has the appropriate entry replaced.

        INPUT:

        - ``cell`` -- a pair of nonnegative integers

        OUTPUT:

        The tableau ``self`` with the entry in cell ``cell`` set to ``m``. This
        entry overwrites an existing entry if ``cell`` already belongs to
        ``self``, or is added to the tableau if ``cell`` is a cocorner of the
        shape ``self``. (Either way, the input is not modified.)

        .. NOTE::

            Both coordinates of ``cell`` are interpreted as starting at `0`.
            So, ``cell == (0, 0)`` corresponds to the northwesternmost cell.

        EXAMPLES::

            sage: s = StandardTableau([[1,2,5],[3,4]]); s.pp()
              1  2  5
              3  4
            sage: t = s.add_entry( (1,2), 6); t.pp()
              1  2  5
              3  4  6
            sage: t.category()
            Category of elements of Standard tableaux
            sage: s.add_entry( (2,0), 6).pp()
              1  2  5
              3  4
              6
            sage: u = s.add_entry( (1,2), 3); u.pp()
              1  2  5
              3  4  3
            sage: u.category()
            Category of elements of Tableaux
            sage: s.add_entry( (2,2),3)
            Traceback (most recent call last):
            ...
            IndexError: (2, 2) is not an addable cell of the tableau
        """
        tab = self.to_list()
        (r, c) = cell
        try:
            tab[r][c] = m   # will work if we are replacing an entry
        except IndexError:
            # Only add a new row if (r,c) is an addable cell (previous code
            # added m to the end of row r independently of the value of c)
            if r >= len(tab):
                if r == len(tab) and c == 0:
                    tab.append([m])
                else:
                    raise IndexError('%s is not an addable cell of the tableau' % ((r, c),))
            else:
                tab_r = tab[r]
                if c == len(tab_r):
                    tab_r.append(m)
                else:
                    raise IndexError('%s is not an addable cell of the tableau' % ((r, c),))

        # attempt to return a tableau of the same type as self
        if tab in self.parent():
            return self.parent()(tab)
        else:
            try:
                return self.parent().Element(tab)
            except Exception:
                return Tableau(tab)

    ##############
    # catabolism #
    ##############

    def catabolism(self):
        """
        Remove the top row of ``self`` and insert it back in using
        column Schensted insertion (starting with the largest letter).

        EXAMPLES::

            sage: Tableau([]).catabolism()
            []
            sage: Tableau([[1,2,3,4,5]]).catabolism()
            [[1, 2, 3, 4, 5]]
            sage: Tableau([[1,1,3,3],[2,3],[3]]).catabolism()
            [[1, 1, 2, 3, 3, 3], [3]]
            sage: Tableau([[1, 1, 2, 3, 3, 3], [3]]).catabolism()
            [[1, 1, 2, 3, 3, 3, 3]]
        """
        h = self.height()
        if h == 0:
            return self
        else:
            # Remove the top row and insert it back in
            return Tableau(self[1:]).insert_word(self[0], left=True)

    def catabolism_sequence(self):
        """
        Perform :meth:`catabolism` on ``self`` until it returns a
        tableau consisting of a single row.

        EXAMPLES::

            sage: t = Tableau([[1,2,3,4,5,6,8],[7,9]])
            sage: t.catabolism_sequence()
            [[[1, 2, 3, 4, 5, 6, 8], [7, 9]],
             [[1, 2, 3, 4, 5, 6, 7, 9], [8]],
             [[1, 2, 3, 4, 5, 6, 7, 8], [9]],
             [[1, 2, 3, 4, 5, 6, 7, 8, 9]]]
            sage: Tableau([]).catabolism_sequence()
            [[]]
        """
        h = self.height()
        res = [self]
        newterm = self
        while h > 1:
            newterm = newterm.catabolism()
            res.append(newterm)
            h = newterm.height()
        return res

    def lambda_catabolism(self, part):
        r"""
        Return the ``part``-catabolism of ``self``, where ``part`` is a
        partition (which can be just given as an array).

        For a partition `\lambda` and a tableau `T`, the
        `\lambda`-catabolism of `T` is defined by performing the following
        steps.

        1. Truncate the parts of `\lambda` so that `\lambda` is contained
           in the shape of `T`.  Let `m` be the length of this partition.

        2. Let `T_a` be the first `m` rows of `T`, and `T_b` be the
           remaining rows.

        3. Let `S_a` be the skew tableau `T_a / \lambda`.

        4. Concatenate the reading words of `S_a` and `T_b`, and insert
           into a tableau.

        EXAMPLES::

            sage: Tableau([[1,1,3],[2,4,5]]).lambda_catabolism([2,1])
            [[3, 5], [4]]
            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: t.lambda_catabolism([])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t.lambda_catabolism([1])
            [[1, 2, 3, 3, 3], [3]]
            sage: t.lambda_catabolism([1,1])
            [[1, 3, 3, 3], [3]]
            sage: t.lambda_catabolism([2,1])
            [[3, 3, 3, 3]]
            sage: t.lambda_catabolism([4,2,1])
            []
            sage: t.lambda_catabolism([5,1])
            [[3, 3]]
            sage: t.lambda_catabolism([4,1])
            [[3, 3]]
        """
        # Reduce the partition if it is too big for the tableau
        part = [min(part[i], len(self[i]))
                for i in range(min(len(self), len(part)))]
        if self.shape() == part:
            return Tableau([])

        m = len(part)

        w1 = list(sum((row for row in reversed(self[m:])), ()))

        w2 = []
        for i, row in enumerate(reversed(self[:m])):
            w2 += row[part[-1 - i]:]

        return Tableau([]).insert_word(w2 + w1)

    def reduced_lambda_catabolism(self, part):
        """
        EXAMPLES::

            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: t.reduced_lambda_catabolism([])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t.reduced_lambda_catabolism([1])
            [[1, 2, 3, 3, 3], [3]]
            sage: t.reduced_lambda_catabolism([1,1])
            [[1, 3, 3, 3], [3]]
            sage: t.reduced_lambda_catabolism([2,1])
            [[3, 3, 3, 3]]
            sage: t.reduced_lambda_catabolism([4,2,1])
            []
            sage: t.reduced_lambda_catabolism([5,1])
            0
            sage: t.reduced_lambda_catabolism([4,1])
            0
        """
        part1 = part

        if self == []:
            return self

        res = self.lambda_catabolism(part)

        if res == []:
            return res

        if res == 0:
            return 0

        a = self[0][0]

        part = [min(part1[i], len(self[i])) for i in range(min(len(part1), len(self)))]
        tt_part = Tableau([[a+i]*part[i] for i in range(len(part))])
        t_part = Tableau([[self[i][j] for j in range(part[i])] for i in range(len(part))])

        if t_part == tt_part:
            return res
        else:
            return 0

    def catabolism_projector(self, parts):
        """
        EXAMPLES::

            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: t.catabolism_projector([[4,2,1]])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t.catabolism_projector([[1]])
            []
            sage: t.catabolism_projector([[2,1],[1]])
            []
            sage: t.catabolism_projector([[1,1],[4,1]])
            [[1, 1, 3, 3], [2, 3], [3]]
        """
        res = self
        for p in parts:
            res = res.reduced_lambda_catabolism(p)
            if res == 0:
                return 0

        if res == []:
            return self
        else:
            return Tableau([])

    def promotion_operator(self, i):
        r"""
        Return a list of semistandard tableaux obtained by the `i`-th
        Lapointe-Lascoux-Morse promotion operator from the
        semistandard tableau ``self``.

        .. WARNING::

            This is not Schuetzenberger's jeu de taquin promotion!
            For the latter, see :meth:`promotion` and
            :meth:`promotion_inverse`.

        This operator is defined by taking the maximum entry `m` of
        `T`, then adding a horizontal `i`-strip to `T` in all possible
        ways, each time filling this strip with `m+1`'s, and finally
        letting the permutation
        `\sigma_1 \sigma_2 \cdots \sigma_m = (2, 3, \ldots, m+1, 1)`
        act on each of the resulting tableaux via the
        Lascoux-Schuetzenberger action
        (:meth:`symmetric_group_action_on_values`). This method
        returns the list of all resulting tableaux. See [LLM2003]_ for
        the purpose of this operator.

        EXAMPLES::

            sage: t = Tableau([[1,2],[3]])
            sage: t.promotion_operator(1)
            [[[1, 2, 4], [3]], [[1, 2], [3, 4]], [[1, 2], [3], [4]]]
            sage: t.promotion_operator(2)
            [[[1, 1, 2, 4], [3]],
             [[1, 1, 4], [2, 3]],
             [[1, 1, 2], [3], [4]],
             [[1, 1], [2, 3], [4]]]
            sage: Tableau([[1]]).promotion_operator(2)
            [[[1, 1, 2]], [[1, 1], [2]]]
            sage: Tableau([[1,1],[2]]).promotion_operator(3)
            [[[1, 1, 1, 2, 3], [2]],
             [[1, 1, 1, 3], [2, 2]],
             [[1, 1, 1, 2], [2], [3]],
             [[1, 1, 1], [2, 2], [3]]]

        The example from [LLM2003]_ p. 12::

            sage: Tableau([[1,1],[2,2]]).promotion_operator(3)
            [[[1, 1, 1, 3, 3], [2, 2]],
             [[1, 1, 1, 3], [2, 2], [3]],
             [[1, 1, 1], [2, 2], [3, 3]]]

        TESTS::

            sage: Tableau([]).promotion_operator(2)
            [[[1, 1]]]
            sage: Tableau([]).promotion_operator(1)
            [[[1]]]
        """
        chain = self.to_chain()
        part = self.shape()
        weight = self.weight()
        perm = permutation.from_reduced_word(range(1, len(weight)+1))
        l = part.add_horizontal_border_strip(i)
        ltab = [from_chain(chain + [next]) for next in l]
        return [x.symmetric_group_action_on_values(perm) for x in ltab]

    ##################################
    # actions on tableaux from words #
    ##################################

    def raise_action_from_words(self, f, *args):
        """
        EXAMPLES::

            sage: from sage.combinat.tableau import symmetric_group_action_on_values
            sage: import functools
            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: f = functools.partial(t.raise_action_from_words, symmetric_group_action_on_values)
            sage: f([1,2,3])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: f([3,2,1])
            [[1, 1, 1, 1], [2, 3], [3]]
            sage: f([1,3,2])
            [[1, 1, 2, 2], [2, 2], [3]]
        """
        w = self.to_word()
        w = f(w, *args)
        return from_shape_and_word(self.shape(), w)

    def symmetric_group_action_on_values(self, perm):
        r"""
        Return the image of the semistandard tableau ``self`` under the
        action of the permutation ``perm`` using the
        Lascoux-Schuetzenberger action of the symmetric group `S_n` on
        the semistandard tableaux with ceiling `n`.

        If `n` is a nonnegative integer, then the
        Lascoux-Schuetzenberger action is a group action of the
        symmetric group `S_n` on the set of semistandard Young tableaux
        with ceiling `n` (that is, with entries taken from the set
        `\{1, 2, \ldots, n\}`). It is defined as follows:

        Let `i \in \{1, 2, \ldots, n-1\}`, and let `T` be a
        semistandard tableau with ceiling `n`. Let `w` be the reading
        word (:meth:`to_word`) of `T`. Replace all letters `i` in `w`
        by closing parentheses, and all letters `i+1` in `w` by
        opening parentheses. Whenever an opening parenthesis stands
        left of a closing parenthesis without there being any
        parentheses in between (it is allowed to have letters
        in-between as long as they are not parentheses), consider these
        two parentheses as matched with each other, and replace them
        back by the letters `i+1` and `i`. Repeat this procedure until
        there are no more opening parentheses standing left of closing
        parentheses. Then, let `a` be the number of opening
        parentheses in the word, and `b` the number of closing
        parentheses (notice that all opening parentheses are right of
        all closing parentheses). Replace the first `a` parentheses
        by the letters `i`, and replace the remaining `b` parentheses
        by the letters `i+1`. Let `w'` be the resulting word. Let
        `T'` be the tableau with the same shape as `T` but with reading
        word `w'`. This tableau `T'` can be shown to be semistandard.
        We define the image of `T` under the action of the simple
        transposition `s_i = (i, i+1) \in S_n` to be this tableau `T'`.
        It can be shown that these actions of the transpositions
        `s_1, s_2, \ldots, s_{n-1}` satisfy the Moore-Coxeter relations
        of `S_n`, and thus this extends to a unique action of the
        symmetric group `S_n` on the set of semistandard tableaux with
        ceiling `n`. This is the Lascoux-Schuetzenberger action.

        This action of the symmetric group `S_n` on the set of all
        semistandard tableaux of given shape `\lambda` with entries
        in `\{ 1, 2, \ldots, n \}` is the one defined in
        [Loth02]_ Theorem 5.6.3. In particular, the action of `s_i`
        is denoted by `\sigma_i` in said source. (Beware of the typo
        in the definition of `\sigma_i`: it should say
        `\sigma_i ( a_i^r a_{i+1}^s ) = a_i^s a_{i+1}^r`, not
        `\sigma_i ( a_i^r a_{i+1}^s ) = a_i^s a_{i+1}^s`.)

        EXAMPLES::

            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: t.symmetric_group_action_on_values([1,2,3])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t.symmetric_group_action_on_values([2,1,3])
            [[1, 2, 3, 3], [2, 3], [3]]
            sage: t.symmetric_group_action_on_values([3,1,2])
            [[1, 2, 2, 2], [2, 3], [3]]
            sage: t.symmetric_group_action_on_values([2,3,1])
            [[1, 1, 1, 1], [2, 2], [3]]
            sage: t.symmetric_group_action_on_values([3,2,1])
            [[1, 1, 1, 1], [2, 3], [3]]
            sage: t.symmetric_group_action_on_values([1,3,2])
            [[1, 1, 2, 2], [2, 2], [3]]

        TESTS::

            sage: t = Tableau([])
            sage: t.symmetric_group_action_on_values([])
            []
        """
        return self.raise_action_from_words(symmetric_group_action_on_values, perm)

    #########
    # atoms #
    #########
    def socle(self):
        """
        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).socle()
            2
            sage: Tableau([[1,2,3,4]]).socle()
            4
        """
        h = self.height()
        if h == 0:
            return 0
        w1row = self[0]
        i = 0
        while i < len(w1row)-1:
            if w1row[i+1] != w1row[i] + 1:
                break
            i += 1
        return i+1

    def atom(self):
        """
        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).atom()
            [2, 2]
            sage: Tableau([[1,2,3],[4,5],[6]]).atom()
            [3, 2, 1]
        """
        ll = [t.socle() for t in self.catabolism_sequence()]
        lres = ll[:]
        for i in range(1, len(ll)):
            lres[i] = ll[i] - ll[i-1]
        return lres

    def symmetric_group_action_on_entries(self, w):
        r"""
        Return the tableau obtained form this tableau by acting by the
        permutation ``w``.

        Let `T` be a standard tableau of size `n`, then the action of
        `w \in S_n` is defined by permuting the entries of `T` (recall they
        are `1, 2, \ldots, n`). In particular, suppose the entry at cell
        `(i, j)` is `a`, then the entry becomes `w(a)`. In general, the
        resulting tableau `wT` may *not* be standard.

        .. NOTE::

            This is different than :meth:`symmetric_group_action_on_values`
            which is defined on semistandard tableaux and is guaranteed to
            return a semistandard tableau.

        INPUT:

        - ``w`` -- a permutation

        EXAMPLES::

            sage: StandardTableau([[1,2,4],[3,5]]).symmetric_group_action_on_entries( Permutation(((4,5))) )
            [[1, 2, 5], [3, 4]]
            sage: _.category()
            Category of elements of Standard tableaux
            sage: StandardTableau([[1,2,4],[3,5]]).symmetric_group_action_on_entries( Permutation(((1,2))) )
            [[2, 1, 4], [3, 5]]
            sage: _.category()
            Category of elements of Tableaux
        """
        w = w + [i + 1 for i in range(len(w), self.size())]  # need to ensure that it belongs to Sym_size
        try:
            return self.parent()([[w[entry-1] for entry in row] for row in self])
        except Exception:
            return Tableau([[w[entry-1] for entry in row] for row in self])

    def is_key_tableau(self):
        r"""
        Return ``True`` if ``self`` is a key tableau or ``False`` otherwise.

        A tableau is a *key tableau* if the set of entries in the `j`-th
        column is a subset of the set of entries in the `(j-1)`-st column.

        REFERENCES:

        - [LS1990]_

        - [Wil2010]_

        EXAMPLES::

            sage: t = Tableau([[1,1,1],[2,3],[3]])
            sage: t.is_key_tableau()
            True

            sage: t = Tableau([[1,1,2],[2,3],[3]])
            sage: t.is_key_tableau()
            False
        """
        T_conj = self.conjugate()
        return all(x in T_conj[i-1] for i in range(1, len(T_conj)) for x in T_conj[i])

    def right_key_tableau(self):
        """
        Return the right key tableau of ``self``.

        The right key tableau of a tableau `T` is a key tableau whose entries
        are weakly greater than the corresponding entries in `T`, and whose column
        reading word is subject to certain conditions. See [LS1990]_ for the full definition.

        ALGORITHM:

        The following algorithm follows [Wil2010]_. Note that if `T` is a key tableau
        then the output of the algorithm is `T`.

        To compute the right key tableau `R` of a tableau `T` we iterate over the columns
        of `T`. Let `T_j` be the `j`-th column of `T` and iterate over the entries
        in `T_j` from bottom to top. Initialize the corresponding entry `k` in `R` to be
        the largest entry in `T_j`. Scan the bottom of each column of `T` to the right of
        `T_j`, updating `k` to be the scanned entry whenever the scanned entry is weakly
        greater than `k`. Update `T_j` and all columns to the right by removing all
        scanned entries.

        .. SEEALSO::

            - :meth:`is_key_tableau()`

        EXAMPLES::

            sage: t = Tableau([[1,2],[2,3]])
            sage: t.right_key_tableau()
            [[2, 2], [3, 3]]
            sage: t = Tableau([[1,1,2,4],[2,3,3],[4],[5]])
            sage: t.right_key_tableau()
            [[2, 2, 2, 4], [3, 4, 4], [4], [5]]

        TESTS:

        We check that if we have a key tableau, we return the same tableau::

            sage: t = Tableau([[1,1,1,2], [2,2,2], [4], [5]])
            sage: t.is_key_tableau()
            True
            sage: t.right_key_tableau() == t
            True

        We check that the empty tableau, which is a key tableau, yields itself::

            sage: Tableau([]).right_key_tableau()
            []
        """
        if not self:
            return self
        cols_list = self.conjugate()
        key = [[] for _ in cols_list]

        for i, col_a in enumerate(cols_list):
            right_cols = cols_list[i+1:]
            for elem in reversed(col_a):
                key_val = elem
                update = []
                for col_b in right_cols:
                    if col_b and key_val <= col_b[-1]:
                        key_val = col_b[-1]
                        update.append(col_b[:-1])
                    else:
                        update.append(col_b)
                key[i].insert(0, key_val)
                right_cols = update
        return Tableau(key).conjugate()

    def left_key_tableau(self):
        """
        Return the left key tableau of ``self``.

        The left key tableau of a tableau `T` is the key tableau whose entries
        are weakly less than the corresponding entries in `T`, and whose column
        reading word is subject to certain conditions. See [LS1990]_ for the full definition.

        ALGORITHM:

        The following algorithm follows [Wil2010]_. Note that if `T` is a key tableau
        then the output of the algorithm is `T`.

        To compute the left key tableau `L` of a tableau `T` we iterate over the columns
        of `T`. Let `T_j` be the `j`-th column of `T` and iterate over the entries
        in `T_j` from bottom to top. Initialize the corresponding entry `k` in `L` as the
        largest entry in `T_j`. Scan the columns to the left of `T_j` and with each column
        update `k` to be the lowest entry in that column which is weakly less than `k`.
        Update `T_j` and all columns to the left by removing all scanned entries.

        .. SEEALSO::

            - :meth:`is_key_tableau()`

        EXAMPLES::

            sage: t = Tableau([[1,2],[2,3]])
            sage: t.left_key_tableau()
            [[1, 1], [2, 2]]
            sage: t = Tableau([[1,1,2,4],[2,3,3],[4],[5]])
            sage: t.left_key_tableau()
            [[1, 1, 1, 2], [2, 2, 2], [4], [5]]

        TESTS:

        We check that if we have a key tableau, we return the same tableau::

            sage: t = Tableau([[1,1,1,2], [2,2,2], [4], [5]])
            sage: t.is_key_tableau()
            True
            sage: t.left_key_tableau() == t
            True

        We check that the empty tableau, which is a key tableau, yields itself::

            sage: Tableau([]).left_key_tableau()
            []
        """
        if not self:
            return self
        cols_list = self.conjugate()
        key = [[] for _ in cols_list]
        key[0] = list(cols_list[0])

        from bisect import bisect_right
        for i, col_a in enumerate(cols_list[1:], 1):
            left_cols = cols_list[:i]
            for elem in reversed(col_a):
                key_val = elem
                update = []
                for col_b in reversed(left_cols):
                    j = bisect_right(col_b, key_val) - 1
                    key_val = col_b[j]
                    update.insert(0, col_b[:j])
                left_cols = update
                key[i].insert(0, key_val)
        return Tableau(key).conjugate()

    #################
    # seg and flush #
    #################
    def _segments(self):
        r"""
        Internal function returning the set of segments of a tableau as
        a dictionary.

        OUTPUT:

        - A dictionary with items of the form ``{(r,k):c}``, where ``r`` is the
          row the ``k``-segment appears and ``c`` is the column the left-most
          box of the ``k``-segment appears.

        EXAMPLES::

            sage: t = Tableau([[1,1,2,3,5],[2,3,5,5],[3,4]])
            sage: sorted(t._segments().items())
            [((0, 2), 2), ((0, 3), 3), ((0, 5), 4), ((1, 3), 1), ((1, 5), 2), ((2, 4), 1)]

            sage: B = crystals.Tableaux("A4", shape=[4,3,2,1])                          # needs sage.modules
            sage: t = B[31].to_tableau()                                                # needs sage.modules
            sage: sorted(t._segments().items())                                         # needs sage.modules
            [((0, 5), 3), ((1, 4), 2), ((2, 4), 1)]
        """
        segments = {}
        for r, row in enumerate(self):
            for c in range(len(row)):
                for j in range(c + 1):
                    if row[j] != r + 1 and (r, row[j]) not in segments:
                        segments[(r, row[j])] = j
        return segments

    def seg(self):
        r"""
        Return the total number of segments in ``self``, as in [Sal2014]_.

        Let `T` be a tableaux.  We define a `k`-*segment* of `T` (in the `i`-th
        row) to be a maximal consecutive sequence of `k`-boxes in the `i`-th
        row for any `i+1 \le k \le r+1`.  Denote the total number of
        `k`-segments in `T` by `\mathrm{seg}(T)`.

        REFERENCES:

        - [Sal2014]_

        EXAMPLES::

            sage: t = Tableau([[1,1,2,3,5],[2,3,5,5],[3,4]])
            sage: t.seg()
            6

            sage: B = crystals.Tableaux("A4", shape=[4,3,2,1])                          # needs sage.modules
            sage: t = B[31].to_tableau()                                                # needs sage.modules
            sage: t.seg()                                                               # needs sage.modules
            3
        """
        return len(self._segments())

    def flush(self):
        r"""
        Return the number of flush segments in ``self``, as in [Sal2014]_.

        Let `1 \le i < k \le r+1` and suppose `\ell` is the smallest integer
        greater than `k` such that there exists an `\ell`-segment in the
        `(i+1)`-st row of `T`.  A `k`-segment in the `i`-th row of `T` is
        called *flush* if the leftmost box in the `k`-segment and the leftmost
        box of the `\ell`-segment are in the same column of `T`.  If, however,
        no such `\ell` exists, then this `k`-segment is said to be *flush* if
        the number of boxes in the `k`-segment is equal to `\theta_i`, where
        `\theta_i = \lambda_i - \lambda_{i+1}` and the shape of `T` is
        `\lambda = (\lambda_1 > \lambda_2 > \cdots > \lambda_r)`.  Denote the
        number of flush `k`-segments in `T` by `\mathrm{flush}(T)`.

        EXAMPLES::

            sage: t = Tableau([[1,1,2,3,5],[2,3,5,5],[3,4]])
            sage: t.flush()
            3

            sage: B = crystals.Tableaux("A4", shape=[4,3,2,1])                          # needs sage.modules
            sage: t = B[32].to_tableau()                                                # needs sage.modules
            sage: t.flush()                                                             # needs sage.modules
            4
        """
        for i in range(len(self)-1):
            if len(self[i]) <= len(self[i+1]):
                raise ValueError('only defined for tableaux with strictly decreasing parts')
        f = 0
        S = self._segments().items()
        for s in S:
            if (s[0][0] != len(self)-1 and s[1] == len(self[s[0][0]+1])
                and self[s[0][0]+1][-1] <= s[0][1]) \
                    or (s[0][0] == len(self)-1 and s[1] == 0):
                f += 1
            else:
                for t in S:
                    if s[0][0]+1 == t[0][0] and s[1] == t[1] and (
                            (s[1] >= 1 and self[s[0][0]+1][s[1]-1] <= self[s[0][0]][s[1]])
                            or (s[1] < 1 and self[s[0][0]+1][s[1]] != s[0][0]+2)):
                        f += 1
        return f

    ##################################
    # contents, residues and degrees #
    ##################################

    def content(self, k, multicharge=[0]):
        """
        Return the content of ``k`` in the standard tableau ``self``.

        The content of `k` is `c - r` if `k` appears in row `r` and
        column `c` of the tableau.

        The ``multicharge`` is a list of length 1 which gives an offset for
        all of the contents. It is included mainly for compatibility with
        :meth:`sage.combinat.tableau_tuple.TableauTuple`.

        EXAMPLES::

            sage: StandardTableau([[1,2],[3,4]]).content(3)
            -1

            sage: StandardTableau([[1,2],[3,4]]).content(6)
            Traceback (most recent call last):
            ...
            ValueError: 6 does not appear in tableau
        """
        for r, row in enumerate(self):
            try:
                return row.index(k) - r + multicharge[0]
            except ValueError:
                pass
        raise ValueError("%d does not appear in tableau" % k)

    def residue(self, k, e, multicharge=(0,)):
        r"""
        Return the residue of the integer ``k`` in the tableau ``self``.

        The *residue* of `k` in a standard tableau is `c - r + m`
        in `\ZZ / e\ZZ`, where `k` appears in row `r` and column `c`
        of the tableau with multicharge `m`.

        INPUT:

        - ``k`` -- integer in `\{1, 2, \ldots, n\}`
        - ``e`` -- integer in `\{0, 2, 3, 4, 5, \ldots\}`
        - ``multicharge`` -- (default: ``[0]``) a list of length 1

        Here `n` is its size of ``self``.

        The ``multicharge`` is a list of length 1 which gives an offset for
        all of the contents. It is included mainly for compatibility with
        :meth:`~sage.combinat.tableau_tuples.TableauTuple.residue`.

        OUTPUT: the residue in `\ZZ / e\ZZ`

        EXAMPLES::

            sage: StandardTableau([[1,2,5],[3,4]]).residue(1,3)
            0
            sage: StandardTableau([[1,2,5],[3,4]]).residue(2,3)
            1
            sage: StandardTableau([[1,2,5],[3,4]]).residue(3,3)
            2
            sage: StandardTableau([[1,2,5],[3,4]]).residue(4,3)
            0
            sage: StandardTableau([[1,2,5],[3,4]]).residue(5,3)
            2
            sage: StandardTableau([[1,2,5],[3,4]]).residue(6,3)
            Traceback (most recent call last):
            ...
            ValueError: 6 does not appear in the tableau
        """
        for r, row in enumerate(self):
            try:
                return IntegerModRing(e)(row.index(k) - r + multicharge[0])
            except ValueError:
                pass
        raise ValueError('%d does not appear in the tableau' % k)

    def residue_sequence(self, e, multicharge=(0,)):
        r"""
        Return the :class:`sage.combinat.tableau_residues.ResidueSequence`
        of the tableau ``self``.

        INPUT:

        - ``e`` -- integer in `\{0, 2, 3, 4, 5, \ldots\}`
        - ``multicharge`` -- (default: ``[0]``) a sequence of integers
          of length 1

        The `multicharge` is a list of length 1 which gives an offset for
        all of the contents. It is included mainly for compatibility with
        :meth:`~sage.combinat.tableau_tuples.StandardTableauTuple.residue`.

        OUTPUT:

        The corresponding residue sequence of the tableau;
        see :class:`ResidueSequence`.

        EXAMPLES::

            sage: StandardTableauTuple([[1,2],[3,4]]).residue_sequence(2)               # needs sage.groups
            2-residue sequence (0,1,1,0) with multicharge (0)
            sage: StandardTableauTuple([[1,2],[3,4]]).residue_sequence(3)               # needs sage.groups
            3-residue sequence (0,1,2,0) with multicharge (0)
            sage: StandardTableauTuple([[1,2],[3,4]]).residue_sequence(4)               # needs sage.groups
            4-residue sequence (0,1,3,0) with multicharge (0)
        """
        res = [0] * self.size()
        for r, row in enumerate(self):
            for c, entry in enumerate(row):
                res[entry-1] = multicharge[0] - r + c
        from sage.combinat.tableau_residues import ResidueSequence
        return ResidueSequence(e, multicharge, res, check=False)

    def degree(self, e, multicharge=(0,)):
        r"""
        Return the Brundan-Kleshchev-Wang [BKW2011]_ degree of ``self``.

        The *degree* is an integer that is defined recursively by successively
        stripping off the number `k`, for `k = n, n-1, \ldots, 1` and at stage
        adding the number of addable cell of the same residue minus the number
        of removable cells of the same residue as `k` and which are below `k`
        in the diagram.

        The degrees of the tableau `T` gives the degree of the homogeneous
        basis element of the graded Specht module that is indexed by `T`.

        INPUT:

        - ``e`` -- the *quantum characteristic*
        - ``multicharge`` -- (default: ``[0]``) the multicharge

        OUTPUT: the degree of the tableau ``self``, which is an integer

        EXAMPLES::

            sage: StandardTableau([[1,2,5],[3,4]]).degree(3)                            # needs sage.groups
            0
            sage: StandardTableau([[1,2,5],[3,4]]).degree(4)                            # needs sage.groups
            1
        """
        n = self.size()
        if n == 0:
            return 0

        deg = self.shape()._initial_degree(e, multicharge)
        res = self.shape().initial_tableau().residue_sequence(e, multicharge)
        for r in self.reduced_row_word():
            if res[r] == res[r+1]:
                deg -= 2
            elif res[r] == res[r+1] + 1 or res[r] == res[r+1] - 1:
                deg += (e == 2 and 2 or 1)
            res = res.swap_residues(r, r+1)
        return deg

    def codegree(self, e, multicharge=(0,)):
        r"""
        Return the Brundan-Kleshchev-Wang [BKW2011]_ codegree of the
        standard tableau ``self``.

        The *codegree* of a tableau is an integer that is defined recursively by
        successively stripping off the number `k`, for `k = n, n-1, \ldots, 1`
        and at stage adding the number of addable cell of the same residue
        minus the number of removable cells of the same residue as `k` and
        are above `k` in the diagram.

        The codegree of the tableau `T` gives the degree of  "dual"
        homogeneous basis element of the Graded Specht module that
        is indexed by `T`.

        INPUT:

        - ``e`` -- the *quantum characteristic*
        - ``multicharge`` -- (default: ``[0]``) the multicharge

        OUTPUT: the codegree of the tableau ``self``, which is an integer

        EXAMPLES::

            sage: StandardTableau([[1,3,5],[2,4]]).codegree(3)                          # needs sage.groups
            0
            sage: StandardTableau([[1,2,5],[3,4]]).codegree(3)                          # needs sage.groups
            1
            sage: StandardTableau([[1,2,5],[3,4]]).codegree(4)                          # needs sage.groups
            0
        """
        if not self:  # the trivial case
            return 0

        conj_shape = self.shape().conjugate()
        codeg = conj_shape._initial_degree(e)
        res = conj_shape.initial_tableau().residue_sequence(e)
        for r in self.reduced_column_word():
            if res[r] == res[r+1]:
                codeg -= 2
            elif res[r] == res[r+1] + 1 or res[r] == res[r+1] - 1:
                codeg += (e == 2 and 2 or 1)
            res = res.swap_residues(r, r+1)
        return codeg

    def first_row_descent(self):
        r"""
        Return the first cell where the tableau ``self`` is not row standard.

        Cells are ordered left to right along the rows and then top to bottom.
        That is, the cell `(r,c)` with `r` and `c` minimal such that the entry
        in position `(r,c)` is bigger than the entry in position `(r, c+1)`.
        If there is no such cell then ``None`` is returned - in this case the
        tableau is row strict.

        OUTPUT:

        The first cell which there is a descent or ``None`` if no such
        cell exists.

        EXAMPLES::

            sage: t = Tableau([[1,3,2],[4]]); t.first_row_descent()
            (0, 1)
            sage: Tableau([[1,2,3],[4]]).first_row_descent() is  None
            True
        """
        for row in range(len(self)):
            for col in range(len(self[row])-1):
                if self[row][col] > self[row][col+1]:
                    return (row, col)
        return None

    def first_column_descent(self):
        r"""
        Return the first cell where ``self`` is not column standard.

        Cells are ordered left to right along the rows and then top to bottom.
        That is, the cell `(r, c)` with `r` and `c` minimal such that
        the entry in position `(r, c)` is bigger than the entry in position
        `(r, c+1)`. If there is no such cell then ``None`` is returned - in
        this case the tableau is column strict.

        OUTPUT:

        The first cell which there is a descent or ``None`` if no such
        cell exists.

        EXAMPLES::

            sage: Tableau([[1,4,5],[2,3]]).first_column_descent()
            (0, 1)
            sage: Tableau([[1,2,3],[4]]).first_column_descent() is None
            True
        """
        for row in range(len(self)-1):
            col = 0
            while col < len(self[row+1]):
                if self[row][col] > self[row+1][col]:
                    return (row, col)
                col += 1
        return None

    def reduced_row_word(self):
        r"""
        Return the lexicographically minimal reduced expression for the
        permutation that maps the :meth:`initial_tableau` to ``self``.

        Ths reduced expression is a minimal length coset representative for the
        corresponding Young subgroup.  In one line notation, the permutation is
        obtained by concatenating the rows of the tableau in order from top to
        bottom.

        EXAMPLES::

            sage: StandardTableau([[1,2,3],[4,5],[6]]).reduced_row_word()
            []
            sage: StandardTableau([[1,2,3],[4,6],[5]]).reduced_row_word()
            [5]
            sage: StandardTableau([[1,2,4],[3,6],[5]]).reduced_row_word()
            [3, 5]
            sage: StandardTableau([[1,2,5],[3,6],[4]]).reduced_row_word()
            [3, 5, 4]
            sage: StandardTableau([[1,2,6],[3,5],[4]]).reduced_row_word()
            [3, 4, 5, 4]
        """
        return permutation.Permutation(list(self.entries())).inverse().reduced_word_lexmin()

    def reduced_column_word(self):
        r"""
        Return the lexicographically minimal reduced expression for the
        permutation that maps the conjugate of the :meth:`initial_tableau`
        to ``self``.

        Ths reduced expression is a minimal length coset representative for
        the corresponding Young subgroup.  In one line notation, the
        permutation is obtained by concatenating the columns of the
        tableau in order from top to bottom.

        EXAMPLES::

            sage: StandardTableau([[1,4,6],[2,5],[3]]).reduced_column_word()
            []
            sage: StandardTableau([[1,4,5],[2,6],[3]]).reduced_column_word()
            [5]
            sage: StandardTableau([[1,3,6],[2,5],[4]]).reduced_column_word()
            [3]
            sage: StandardTableau([[1,3,5],[2,6],[4]]).reduced_column_word()
            [3, 5]
            sage: StandardTableau([[1,2,5],[3,6],[4]]).reduced_column_word()
            [3, 2, 5]
        """
        data = list(self.conjugate().entries())
        return permutation.Permutation(data).inverse().reduced_word_lexmin()

    def hillman_grassl(self):
        r"""
        Return the image of the `\lambda`-array ``self`` under the
        Hillman-Grassl correspondence (as a
        :class:`~sage.combinat.hillman_grassl.WeakReversePlanePartition`).

        This relies on interpreting ``self`` as a `\lambda`-array
        in the sense of :mod:`~sage.combinat.hillman_grassl`.

        Fix a partition `\lambda`
        (see :meth:`~sage.combinat.partition.Partition`).
        We draw all partitions and tableaux in English notation.

        A `\lambda`-*array* will mean a tableau of shape `\lambda` whose
        entries are nonnegative integers. (No conditions on the order of
        these entries are made. Note that `0` is allowed.)

        A *weak reverse plane partition of shape* `\lambda` (short:
        `\lambda`-*rpp*) will mean a `\lambda`-array whose entries weakly
        increase along each row and weakly increase along each column.

        The Hillman-Grassl correspondence `H` is the map that sends a
        `\lambda`-array `M` to a `\lambda`-rpp `H(M)` defined recursively
        as follows:

        * If all entries of `M` are `0`, then `H(M) = M`.

        * Otherwise, let `s` be the index of the leftmost column of `M`
          containing a nonzero entry.
          Let `r` be the index of the bottommost nonzero entry in the
          `s`-th column of `M`. Let `M'` be the `\lambda`-array obtained
          from `M` by subtracting `1` from the `(r, s)`-th entry of `M`.
          Let `Q = (q_{i, j})` be the image `H(M')` (which is already
          defined by recursion).

        * Define a sequence `((i_1, j_1), (i_2, j_2), \ldots,
          (i_n, j_n))` of boxes in the diagram of `\lambda` (actually a
          lattice path made of southward and westward steps) as follows:
          Set `(i_1, j_1) = (r, \lambda_r)` (the rightmost box in the
          `r`-th row of `\lambda`). If `(i_k, j_k)` is defined for some
          `k \geq 1`, then `(i_{k+1}, j_{k+1})` is constructed as follows:
          If `q_{i_k + 1, j_k}` is well-defined and equals `q_{i_k, j_k}`,
          then we set `(i_{k+1}, j_{k+1}) = (i_k + 1, j_k)`.
          Otherwise, if `j_k = s`, then the sequence ends here.
          Otherwise, we set `(i_{k+1}, j_{k+1}) = (i_k, j_k - 1)`.

        * Let `H(M)` be the array obtained from `Q` by adding `1` to
          the `(i_k, j_k)`-th entry of `Q` for each
          `k \in \{1, 2, \ldots, n\}`.

        See [Gans1981]_ (Section 3) for this construction.

        .. SEEALSO::

            :meth:`~sage.combinat.hillman_grassl.hillman_grassl`
            for the Hillman-Grassl correspondence as a standalone
            function.

            :meth:`~sage.combinat.hillman_grassl.WeakReversePlanePartition.hillman_grassl_inverse`
            for the inverse map.

        EXAMPLES::

            sage: a = Tableau([[2, 1, 1], [0, 2, 0], [1, 1]])
            sage: A = a.hillman_grassl(); A
            [[2, 2, 4], [2, 3, 4], [3, 5]]
            sage: A.parent(), a.parent()
            (Weak Reverse Plane Partitions, Tableaux)
        """
        from sage.combinat.hillman_grassl import (hillman_grassl,
                                                  WeakReversePlanePartition)
        return WeakReversePlanePartition(hillman_grassl(list(self)))

    def sulzgruber_correspondence(self):
        r"""
        Return the image of the `\lambda`-array ``self`` under the
        Sulzgruber correspondence (as a
        :class:`~sage.combinat.hillman_grassl.WeakReversePlanePartition`).

        This relies on interpreting ``self`` as a `\lambda`-array
        in the sense of :mod:`~sage.combinat.hillman_grassl`.
        See :mod:`~sage.combinat.hillman_grassl` for definitions
        of the objects involved.

        The Sulzgruber correspondence is the map `\Phi_\lambda`
        from [Sulzgr2017]_ Section 7, and is the map
        `\xi_\lambda^{-1}` from [Pak2002]_ Section 5.
        It is denoted by `\mathcal{RSK}` in [Hopkins2017]_.
        It is the inverse of the Pak correspondence
        (:meth:`pak_correspondence`).
        The following description of the Sulzgruber correspondence
        follows [Hopkins2017]_ (which denotes it by `\mathcal{RSK}`):

        Fix a partition `\lambda`
        (see :meth:`~sage.combinat.partition.Partition`).
        We draw all partitions and tableaux in English notation.

        A `\lambda`-*array* will mean a tableau of shape `\lambda` whose
        entries are nonnegative integers. (No conditions on the order of
        these entries are made. Note that `0` is allowed.)

        A *weak reverse plane partition of shape* `\lambda` (short:
        `\lambda`-*rpp*) will mean a `\lambda`-array whose entries weakly
        increase along each row and weakly increase along each column.

        We shall also use the following notation:
        If `(u, v)` is a cell of `\lambda`, and if `\pi` is a
        `\lambda`-rpp, then:

        * the *lower bound* of `\pi` at `(u, v)` (denoted by
          `\pi_{<(u, v)}`) is defined to be
          `\max \{ \pi_{u-1, v} , \pi_{u, v-1} \}` (where
          `\pi_{0, v}` and `\pi_{u, 0}` are understood to mean `0`).

        * the *upper bound* of `\pi` at `(u, v)` (denoted by
          `\pi_{>(u, v)}`) is defined to be
          `\min \{ \pi_{u+1, v} , \pi_{u, v+1} \}`
          (where `\pi_{i, j}` is understood to mean `+ \infty`
          if `(i, j)` is not in `\lambda`; thus, the upper
          bound at a corner cell is `+ \infty`).

        * *toggling* `\pi` at `(u, v)` means replacing the entry
          `\pi_{u, v}` of `\pi` at `(u, v)` by
          `\pi_{<(u, v)} + \pi_{>(u, v)} - \pi_{u, v}`
          (this is well-defined as long as `(u, v)` is not a
          corner of `\lambda`).

        Note that every `\lambda`-rpp `\pi` and every cell
        `(u, v)` of `\lambda` satisfy
        `\pi_{<(u, v)} \leq \pi_{u, v} \leq \pi_{>(u, v)}`.
        Note that toggling a `\lambda`-rpp (at a cell that is not
        a corner) always results in a `\lambda`-rpp. Also,
        toggling is an involution).

        The Pak correspondence `\xi_\lambda` sends a `\lambda`-rpp `\pi`
        to a `\lambda`-array `\xi_\lambda(\pi)`. It is defined by
        recursion on `\lambda` (that is, we assume that `\xi_\mu` is
        already defined for every partition `\mu` smaller than
        `\lambda`), and its definition proceeds as follows:

        * If `\lambda = \varnothing`, then `\xi_\lambda` is the
          obvious bijection sending the only `\varnothing`-rpp
          to the only `\varnothing`-array.

        * Pick any corner `c = (i, j)` of `\lambda`, and let `\mu`
          be the result of removing this corner `c` from the partition
          `\lambda`.
          (The exact choice of `c` is immaterial.)

        * Let `\pi'` be what remains of `\pi` when the corner cell `c`
          is removed.

        * For each positive integer `k` such that `(i-k, j-k)` is a
          cell of `\lambda`, toggle `\pi'` at `(i-k, j-k)`.
          (All these togglings commute, so the order in which they
          are made is immaterial.)

        * Let `M = \xi_\mu(\pi')`.

        * Extend the `\mu`-array `M` to a `\lambda`-array `M'` by
          adding the cell `c` and writing the number
          `\pi_{i, j} - \pi_{<(i, j)}` into this cell.

        * Set `\xi_\lambda(\pi) = M'`.

        .. SEEALSO::

            :meth:`~sage.combinat.hillman_grassl.sulzgruber_correspondence`
            for the Sulzgruber correspondence as a standalone function.

            :meth:`~sage.combinat.hillman_grassl.WeakReversePlanePartition.pak_correspondence`
            for the inverse map.

        EXAMPLES::

            sage: a = Tableau([[2, 1, 1], [0, 2, 0], [1, 1]])
            sage: A = a.sulzgruber_correspondence(); A
            [[0, 1, 4], [1, 5, 5], [3, 6]]
            sage: A.parent(), a.parent()
            (Weak Reverse Plane Partitions, Tableaux)

            sage: a = Tableau([[1, 3], [0, 1]])
            sage: a.sulzgruber_correspondence()
            [[0, 4], [1, 5]]
        """
        from sage.combinat.hillman_grassl import (sulzgruber_correspondence,
                                                  WeakReversePlanePartition)
        return WeakReversePlanePartition(sulzgruber_correspondence(list(self)))


class SemistandardTableau(Tableau):
    """
    A class to model a semistandard tableau.

    INPUT:

    - ``t`` -- a tableau, a list of iterables, or an empty list

    OUTPUT: a SemistandardTableau object constructed from ``t``

    A semistandard tableau is a tableau whose entries are positive integers,
    which are weakly increasing in rows and strictly increasing down columns.

    EXAMPLES::

        sage: t = SemistandardTableau([[1,2,3],[2,3]]); t
        [[1, 2, 3], [2, 3]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty printing
        1 2 3
        2 3
        sage: t = Tableau([[1,2],[2]])
        sage: s = SemistandardTableau(t); s
        [[1, 2], [2]]
        sage: SemistandardTableau([]) # The empty tableau
        []

    When using code that will generate a lot of tableaux, it is slightly more
    efficient to construct a SemistandardTableau from the appropriate
    :class:`Parent` object::

        sage: SST = SemistandardTableaux()
        sage: SST([[1, 2, 3], [4, 5]])
        [[1, 2, 3], [4, 5]]

    .. SEEALSO::

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`

    TESTS::

        sage: t = Tableaux()([[1,1],[2]])
        sage: s = SemistandardTableaux(3)([[1,1],[2]])
        sage: s == t
        True
        sage: s.parent()
        Semistandard tableaux of size 3 and maximum entry 3
        sage: r = SemistandardTableaux(3)(t); r.parent()
        Semistandard tableaux of size 3 and maximum entry 3
        sage: isinstance(r, Tableau)
        True
        sage: s2 = SemistandardTableaux(3)([(1,1),(2,)])
        sage: s2 == s
        True
        sage: s2.parent()
        Semistandard tableaux of size 3 and maximum entry 3
    """
    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a SemistandardTableau is only ever constructed as an
        element_class call of an appropriate parent.

        TESTS::

            sage: t = SemistandardTableau([[1,1],[2]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Semistandard tableaux
            sage: t.category()
            Category of elements of Semistandard tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.SemistandardTableaux_all_with_category.element_class'>
        """
        if isinstance(t, SemistandardTableau):
            return t
        if t in SemistandardTableaux():
            return SemistandardTableaux_all().element_class(SemistandardTableaux_all(), t)

        # t is not a semistandard tableau so we give an appropriate error message
        if t not in Tableaux():
            raise ValueError('%s is not a tableau' % t)

        for (rix, row) in enumerate(t):
            for (cix, v) in enumerate(row):
                if not isinstance(v, (int, Integer)):
                    raise ValueError("expected entry to be an integer at (row=%s, col=%s)" % (rix, cix))
                if v <= 0:
                    raise ValueError("expected entry to be a positive integer at (row=%s, col=%s). Found (%s)" % (rix, cix, v))

        for (rix, row) in enumerate(t):
            for cix in range(len(row)-1):
                if row[cix] > row[cix+1]:
                    raise ValueError("row (%s) is not weakly increasing between columns (%s, %s)" % (rix, cix, cix+1))

        # If we're still here ``t`` cannot be column strict
        for rix in range(len(t)-1):
            rcur = t[rix]
            rnext = t[rix+1]

            # check that SST is strictly increasing in columns
            # we know that len(rnext) <= len(rcur) as the SST cannot have
            # more columns in the next row than the current row.
            assert len(rnext) <= len(rcur)

            for cix in range(len(rnext)):
                if rnext[cix] <= rcur[cix]:
                    raise ValueError("column (%s) is not strictly increasing between rows (%s, %s)" % (cix, rix, rix+1))

        # we should have found an error by now.
        raise ValueError('we should have found an error by now in tableau %s' % t)

    def check(self):
        """
        Check that ``self`` is a valid semistandard tableau.

        TESTS::

            sage: SemistandardTableau([[1,2,3],[1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: column (0) is not strictly increasing between rows (0, 1)

            sage: SemistandardTableau([[1,2,1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: row (0) is not weakly increasing between columns (1, 2)

            sage: SemistandardTableau([[0,1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: expected entry to be a positive integer at (row=0, col=0). Found (0)
        """
        super().check()

        # Tableau() has checked that t is tableau, so it remains to check that
        # the entries of t are positive integers which are weakly increasing
        # along rows
        from sage.sets.positive_integers import PositiveIntegers
        PI = PositiveIntegers()

        for row in self:
            if any(c not in PI for c in row):
                raise ValueError("the entries of a semistandard tableau must be nonnegative integers")
            if any(row[c] > row[c+1] for c in range(len(row)-1)):
                raise ValueError("the entries in each row of a semistandard tableau must be weakly increasing")

        # and strictly increasing down columns
        if self:
            for row, next in zip(self, self[1:]):
                if not all(row[c] < next[c] for c in range(len(next))):
                    raise ValueError("the entries of each column of a semistandard tableau must be strictly increasing")


class RowStandardTableau(Tableau):
    """
    A class to model a row standard tableau.

    A row standard tableau is a tableau whose entries are
    positive integers from 1 to `m` that increase along rows.

    INPUT:

    - ``t`` -- a Tableau, a list of iterables, or an empty list

    EXAMPLES::

        sage: t = RowStandardTableau([[3,4,5],[1,2]]); t
        [[3, 4, 5], [1, 2]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty printing
        3  4  5
        1  2
        sage: t.is_standard()
        False
        sage: RowStandardTableau([]) # The empty tableau
        []
        sage: RowStandardTableau([[3,4,5],[1,2]]) in StandardTableaux()
        False
        sage: RowStandardTableau([[1,2,5],[3,4]]) in StandardTableaux()
        True

    When using code that will generate a lot of tableaux, it is more
    efficient to construct a :class:`RowStandardTableau` from the
    appropriate :class:`Parent` object::

        sage: ST = RowStandardTableaux()
        sage: ST([[3, 4, 5], [1, 2]])
        [[3, 4, 5], [1, 2]]

    .. SEEALSO::

        - :class:`Tableau`
        - :class:`StandardTableau`
        - :class:`SemistandardTableau`
        - :class:`Tableaux`
        - :class:`StandardTableaux`
        - :class:`RowStandardTableaux`
        - :class:`SemistandardTableaux`

    TESTS::

        sage: t = Tableaux()([[1,2],[3]])
        sage: s = RowStandardTableaux(3)([[1,2],[3]])
        sage: s == t
        True
        sage: u = RowStandardTableaux(3)([[2, 3], [1]])
        sage: s == u
        False
        sage: u.parent()
        Row standard tableaux of size 3
        sage: u = RowStandardTableaux(3)(u)
        sage: u.parent()
        Row standard tableaux of size 3
        sage: isinstance(u, Tableau)
        True
    """
    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a :class:`RowStandardTableau` is only
        constructed as an ``element_class`` call of an appropriate parent.

        TESTS::

            sage: t = RowStandardTableau([[2,3],[1]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Row standard tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.RowStandardTableaux_all_with_category.element_class'>
        """
        if isinstance(t, RowStandardTableau):
            return t

        RST = RowStandardTableaux_all()
        return RST.element_class(RST, t)

    def check(self):
        r"""
        Check that ``self`` is a valid row standard tableau.

        TESTS::

            sage: RowStandardTableau([[1,4,4],[2,3]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the entries in a row standard tableau must increase
             along rows and contain the numbers 1,2,...,n

            sage: RowStandardTableau([[1,3,2]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the entries in a row standard tableau must increase
             along rows and contain the numbers 1,2,...,n
        """
        super().check()
        # We have checked that t is tableau, so it remains to check that
        #   the entries of t are positive integers that increase along rows.
        flatx = sorted(c for row in self for c in row)
        if (flatx != list(range(1, len(flatx)+1))
                or any(row[i] >= row[i+1] for row in self for i in range(len(row)-1))):
            raise ValueError("the entries in a row standard tableau must increase"
                             " along rows and contain the numbers 1,2,...,n")


class StandardTableau(SemistandardTableau):
    """
    A class to model a standard tableau.

    INPUT:

    - ``t`` -- a Tableau, a list of iterables, or an empty list

    A standard tableau is a semistandard tableau whose entries are exactly the
    positive integers from 1 to `n`, where `n` is the size of the tableau.

    EXAMPLES::

        sage: t = StandardTableau([[1,2,3],[4,5]]); t
        [[1, 2, 3], [4, 5]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty printing
        1 2 3
        4 5
        sage: t.is_standard()
        True
        sage: StandardTableau([]) # The empty tableau
        []
        sage: StandardTableau([[1,2,3],[4,5]]) in RowStandardTableaux()
        True

    When using code that will generate a lot of tableaux, it is more
    efficient to construct a StandardTableau from the appropriate
    :class:`Parent` object::

        sage: ST = StandardTableaux()
        sage: ST([[1, 2, 3], [4, 5]])
        [[1, 2, 3], [4, 5]]

    .. SEEALSO::

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`

    TESTS::

        sage: t = Tableaux()([[1,2],[3]])
        sage: s = StandardTableaux(3)([[1,2],[3]])
        sage: s == t
        True
        sage: s.parent()
        Standard tableaux of size 3
        sage: r = StandardTableaux(3)(t); r.parent()
        Standard tableaux of size 3
        sage: isinstance(r, Tableau)
        True
    """
    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a :class:`StandardTableau` is only ever constructed
        as an ``element_class`` call of an appropriate parent.

        TESTS::

            sage: t = StandardTableau([[1,2],[3]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Standard tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.StandardTableaux_all_with_category.element_class'>
        """
        if isinstance(t, StandardTableau):
            return t

        return StandardTableaux_all().element_class(StandardTableaux_all(), t)

    def check(self):
        """
        Check that ``self`` is a standard tableau.

        TESTS::

            sage: StandardTableau([[1,2,3],[4,4]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the entries in a standard tableau must be in bijection with 1,2,...,n

            sage: StandardTableau([[1,3,2]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the entries in each row of a semistandard tableau must be weakly increasing
        """
        super().check()
        # t is semistandard so we only need to check
        # that its entries are in bijection with {1, 2, ..., n}
        flattened_list = [i for row in self for i in row]
        if sorted(flattened_list) != list(range(1, len(flattened_list) + 1)):
            raise ValueError("the entries in a standard tableau must be in bijection with 1,2,...,n")

    def dominates(self, t):
        r"""
        Return ``True`` if ``self`` dominates the tableau ``t``.

        That is, if the shape of the tableau restricted to `k`
        dominates the shape of ``t`` restricted to `k`, for `k = 1, 2,
        \ldots, n`.

        When the two tableaux have the same shape, then this ordering
        coincides with the Bruhat ordering for the corresponding permutations.

        INPUT:

        - ``t`` -- a tableau

        EXAMPLES::

            sage: s = StandardTableau([[1,2,3],[4,5]])
            sage: t = StandardTableau([[1,2],[3,5],[4]])
            sage: s.dominates(t)
            True
            sage: t.dominates(s)
            False
            sage: all(StandardTableau(s).dominates(t) for t in StandardTableaux([3,2]))
            True
            sage: s.dominates([[1,2,3,4,5]])
            False
        """
        t = StandardTableau(t)
        return all(self.restrict(m).shape().dominates(t.restrict(m).shape())
                   for m in range(1, 1 + self.size()))

    def is_standard(self):
        """
        Return ``True`` since ``self`` is a standard tableau.

        EXAMPLES::

            sage: StandardTableau([[1, 3], [2, 4]]).is_standard()
            True
        """
        return True

    def up(self):
        """
        An iterator for all the standard tableaux that can be
        obtained from ``self`` by adding a cell.

        EXAMPLES::

            sage: t = StandardTableau([[1,2]])
            sage: [x for x in t.up()]
            [[[1, 2, 3]], [[1, 2], [3]]]
        """
        # Get a list of all places where we can add a cell
        # to the shape of self

        outside_corners = self.shape().outside_corners()

        n = self.size()

        # Go through and add n+1 to the end of each
        # of the rows
        for row, _ in outside_corners:
            new_t = [list(_) for _ in self]
            if row != len(self):
                new_t[row] += [n+1]
            else:
                new_t.append([n+1])
            yield StandardTableau(new_t)

    def up_list(self):
        """
        Return a list of all the standard tableaux that can be obtained
        from ``self`` by adding a cell.

        EXAMPLES::

            sage: t = StandardTableau([[1,2]])
            sage: t.up_list()
            [[[1, 2, 3]], [[1, 2], [3]]]
        """
        return list(self.up())

    def down(self):
        """
        An iterator for all the standard tableaux that can be obtained
        from ``self`` by removing a cell. Note that this iterates just
        over a single tableau (or nothing if ``self`` is empty).

        EXAMPLES::

            sage: t = StandardTableau([[1,2],[3]])
            sage: [x for x in t.down()]
            [[[1, 2]]]
            sage: t = StandardTableau([])
            sage: [x for x in t.down()]
            []
        """
        if self:
            yield self.restrict(self.size() - 1)

    def down_list(self):
        """
        Return a list of all the standard tableaux that can be obtained
        from ``self`` by removing a cell. Note that this is just a singleton
        list if ``self`` is nonempty, and an empty list otherwise.

        EXAMPLES::

            sage: t = StandardTableau([[1,2],[3]])
            sage: t.down_list()
            [[[1, 2]]]
            sage: t = StandardTableau([])
            sage: t.down_list()
            []
        """
        return list(self.down())

    def standard_descents(self):
        """
        Return a list of the integers `i` such that `i` appears
        strictly further north than `i + 1` in ``self`` (this is not
        to say that `i` and `i + 1` must be in the same column). The
        list is sorted in increasing order.

        EXAMPLES::

            sage: StandardTableau( [[1,3,4],[2,5]] ).standard_descents()
            [1, 4]
            sage: StandardTableau( [[1,2],[3,4]] ).standard_descents()
            [2]
            sage: StandardTableau( [[1,2,5],[3,4],[6,7],[8],[9]] ).standard_descents()
            [2, 5, 7, 8]
            sage: StandardTableau( [] ).standard_descents()
            []
        """
        descents = []
        # whatpart gives the number for which self is a partition
        whatpart = sum(self.shape())
        # now find the descents
        for i in range(1, whatpart):
            # find out what row i and i+1 are in (we're using the
            # standardness of self here)
            for row in self:
                if row.count(i + 1):
                    break
                if row.count(i):
                    descents.append(i)
                    break
        return descents

    def standard_number_of_descents(self):
        """
        Return the number of all integers `i` such that `i` appears
        strictly further north than `i + 1` in ``self`` (this is not
        to say that `i` and `i + 1` must be in the same column). A
        list of these integers can be obtained using the
        :meth:`standard_descents` method.

        EXAMPLES::

            sage: StandardTableau( [[1,2],[3,4],[5]] ).standard_number_of_descents()
            2
            sage: StandardTableau( [] ).standard_number_of_descents()
            0
            sage: tabs = StandardTableaux(5)
            sage: all( t.standard_number_of_descents() == t.schuetzenberger_involution().standard_number_of_descents() for t in tabs )
            True
        """
        return len(self.standard_descents())

    def standard_major_index(self):
        """
        Return the major index of the standard tableau ``self`` in the
        standard meaning of the word. The major index is defined to be
        the sum of the descents of ``self`` (see :meth:`standard_descents`
        for their definition).

        EXAMPLES::

            sage: StandardTableau( [[1,4,5],[2,6],[3]] ).standard_major_index()
            8
            sage: StandardTableau( [[1,2],[3,4]] ).standard_major_index()
            2
            sage: StandardTableau( [[1,2,3],[4,5]] ).standard_major_index()
            3
        """
        return sum(self.standard_descents())

    def promotion_inverse(self, n=None):
        """
        Return the image of ``self`` under the inverse promotion operator.
        The optional variable `m` should be set to the size of ``self`` minus
        `1` for a minimal speedup; otherwise, it defaults to this number.

        The inverse promotion operator, applied to a standard tableau `t`,
        does the following:

        Remove the letter `1` from `t`, thus leaving a hole where it used to be.
        Apply jeu de taquin to move this hole northeast (in French notation)
        until it reaches the outer boundary of `t`. Fill `n + 1` into this hole,
        where `n` is the size of `t`. Finally, subtract `1` from each letter in
        the tableau. This yields a new standard tableau.

        This definition of inverse promotion is the map called "promotion" in
        [Sag2011]_ (p. 23) and in [Stan2009]_, and is the inverse of the map
        called "promotion" in [Hai1992]_ (p. 90).

        See the :meth:`~sage.combinat.tableau.promotion_inverse` method for a
        more general operator.

        EXAMPLES::

            sage: t = StandardTableau([[1,3],[2,4]])
            sage: t.promotion_inverse()
            [[1, 2], [3, 4]]

        We check the equivalence of two definitions of inverse promotion on
        standard tableaux::

            sage: ST = StandardTableaux(7)
            sage: def bk_promotion_inverse7(st):
            ....:     st2 = st
            ....:     for i in range(1, 7):
            ....:         st2 = st2.bender_knuth_involution(i, check=False)
            ....:     return st2
            sage: all( bk_promotion_inverse7(st) == st.promotion_inverse() for st in ST ) # long time
            True
        """
        if n is None:
            n = self.size() - 1
        return StandardTableau(Tableau(self[:]).promotion_inverse(n))

    def promotion(self, n=None):
        r"""
        Return the image of ``self`` under the promotion operator.

        The promotion operator, applied to a standard tableau `t`, does the
        following:

        Remove the letter `n` from `t`, thus leaving a hole where it used to be.
        Apply jeu de taquin to move this hole southwest (in French notation)
        until it reaches the inner boundary of `t`. Fill `0` into the hole once
        jeu de taquin has completed. Finally, add `1` to each letter in the
        tableau. The resulting standard tableau is the image of `t` under the
        promotion operator.

        This definition of promotion is precisely the one given in [Hai1992]_
        (p. 90). It is the inverse of the maps called "promotion" in [Sag2011]_
        (p. 23) and in [Stan2009]_.

        See the :meth:`~sage.combinat.tableau.promotion` method for a
        more general operator.

        EXAMPLES::

            sage: ST = StandardTableaux(7)
            sage: all( st.promotion().promotion_inverse() == st for st in ST ) # long time
            True
            sage: all( st.promotion_inverse().promotion() == st for st in ST ) # long time
            True
            sage: st = StandardTableau([[1,2,5],[3,4]])
            sage: parent(st.promotion())
            Standard tableaux
        """
        if n is None:
            n = self.size() - 1
        return StandardTableau(Tableau(self[:]).promotion(n))


def from_chain(chain):
    """
    Return a semistandard tableau from a chain of partitions.

    EXAMPLES::

        sage: from sage.combinat.tableau import from_chain
        sage: from_chain([[], [2], [2, 1], [3, 2, 1]])
        [[1, 1, 3], [2, 3], [3]]
    """
    res = [[0] * chain[-1][i] for i in range(len(chain[-1]))]
    for i in reversed(range(2, len(chain) + 1)):
        for j in range(len(chain[i - 1])):
            for k in range(chain[i - 1][j]):
                res[j][k] = i - 1
    T = Tableaux_all()
    return T.element_class(T, res)


def from_shape_and_word(shape, w, convention='French'):
    r"""
    Return a tableau from a shape and word.

    INPUT:

    - ``shape`` -- a partition

    - ``w`` -- a word whose length equals that of the partition

    - ``convention`` -- string (default: ``'French'``); can take values
      ``'French'`` or ``'English'``

    OUTPUT:

    A tableau, whose shape is ``shape`` and whose reading word is ``w``.
    If the ``convention`` is specified as ``'French'``, the reading word is to be read
    starting from the top row in French convention (= the bottom row in English
    convention). If the ``convention`` is specified as ``'English'``, the reading word
    is to be read starting with the top row in English convention.

    EXAMPLES::

        sage: from sage.combinat.tableau import from_shape_and_word
        sage: t = Tableau([[1, 3], [2], [4]])
        sage: shape = t.shape(); shape
        [2, 1, 1]
        sage: word = t.to_word(); word
        word: 4213
        sage: from_shape_and_word(shape, word)
        [[1, 3], [2], [4]]
        sage: word = Word(flatten(t))
        sage: from_shape_and_word(shape, word, convention='English')
        [[1, 3], [2], [4]]
    """
    res = []
    j = 0
    if convention == "French":
        shape = reversed(shape)
    for l in shape:
        res.append(tuple(w[j:j+l]))
        j += l
    if convention == "French":
        res.reverse()
    T = Tableaux_all()
    return T.element_class(T, res)


class IncreasingTableau(Tableau):
    """
    A class to model an increasing tableau.

    INPUT:

    - ``t`` -- a tableau, a list of iterables, or an empty list

    An *increasing tableau* is a tableau whose entries are positive
    integers that are strictly increasing across rows and strictly
    increasing down columns.

    EXAMPLES::

        sage: t = IncreasingTableau([[1,2,3],[2,3]]); t
        [[1, 2, 3], [2, 3]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty printing
        1 2 3
        2 3
        sage: t = Tableau([[1,2],[2]])
        sage: s = IncreasingTableau(t); s
        [[1, 2], [2]]
        sage: IncreasingTableau([]) # The empty tableau
        []

    You can also construct an :class:`IncreasingTableau` from the
    appropriate :class:`Parent` object::

        sage: IT = IncreasingTableaux()
        sage: IT([[1, 2, 3], [4, 5]])
        [[1, 2, 3], [4, 5]]

    .. SEEALSO::

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`
        - :class:`IncreasingTableaux`

    TESTS::

        sage: t = Tableaux()([[1,2],[2]])
        sage: s = IncreasingTableaux(3)([[1,2],[2]])
        sage: s == t
        True
        sage: s.parent()
        Increasing tableaux of size 3 and maximum entry 3
        sage: r = IncreasingTableaux(3)(t); r.parent()
        Increasing tableaux of size 3 and maximum entry 3
        sage: isinstance(r, Tableau)
        True
        sage: s2 = IncreasingTableaux(3)([(1,2),(2,)])
        sage: s2 == s
        True
        sage: s2.parent()
        Increasing tableaux of size 3 and maximum entry 3
    """
    @staticmethod
    def __classcall_private__(self, t):
        r"""
        Construct an :class:`IncreasingTableau` from the appropriate parent.

        TESTS::

            sage: t = IncreasingTableau([[1,2],[2]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Increasing tableaux
            sage: t.category()
            Category of elements of Increasing tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.IncreasingTableaux_all_with_category.element_class'>
        """
        if isinstance(t, IncreasingTableau):
            return t
        IT = IncreasingTableaux()
        return IT.element_class(IT, t)  # The check() will raise the appropriate error

    def check(self):
        """
        Check that ``self`` is a valid increasing tableau.

        TESTS::

            sage: IncreasingTableau([[1,2,3],[1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the entries of each column of an increasing tableau must be strictly increasing

            sage: IncreasingTableau([[1,2,2]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the entries in each row of an increasing tableau must be strictly increasing

            sage: IncreasingTableau([[0,1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the entries of an increasing tableau must be nonnegative integers
        """
        if not self:
            # Empty tableau, so trivially an increasing tableau
            return

        super().check()

        # Tableau() has checked that t is tableau, so it remains to check that
        # the entries of t are positive integers which are weakly increasing
        # along rows
        from sage.sets.positive_integers import PositiveIntegers
        PI = PositiveIntegers()

        for row in self:
            if any(c not in PI for c in row):
                raise ValueError("the entries of an increasing tableau"
                                 " must be nonnegative integers")
            if any(row[c] >= row[c+1] for c in range(len(row)-1)):
                raise ValueError("the entries in each row of an increasing"
                                 " tableau must be strictly increasing")

        # and strictly increasing down columns
        for row, next in zip(self, self[1:]):
            if not all(row[c] < next[c] for c in range(len(next))):
                raise ValueError("the entries of each column of an increasing"
                                 " tableau must be strictly increasing")

    def descent_set(self):
        r"""
        Compute the descents of the increasing tableau ``self``
        as defined in [DPS2017]_.

        The number `i` is a *descent* of an increasing tableau if some
        instance of `i + 1` appears in a lower row than some instance
        of `i`.

        .. NOTE::

            This notion is close to the notion of descent for a standard
            tableau but is unrelated to the notion for semistandard tableaux.

        EXAMPLES::

            sage: T = IncreasingTableau([[1,2,4],[3,5,6]])
            sage: T.descent_set()
            [2, 4]
            sage: U = IncreasingTableau([[1,3,4],[2,4,5]])
            sage: U.descent_set()
            [1, 3, 4]
            sage: V = IncreasingTableau([[1,3,4],[3,4,5],[4,5]])
            sage: V.descent_set()
            [3, 4]
        """
        ans = set()
        ell = len(self)
        for r1, row in enumerate(self):
            for val in row:
                for r2 in range(r1+1, ell):
                    if val + 1 in self[r2]:
                        ans.add(val)
        return sorted(ans)

    @combinatorial_map(order=2, name='K-Bender-Knuth involution')
    def K_bender_knuth(self, i):
        r"""
        Return the ``i``-th K-Bender-Knuth operator (as defined in
        [DPS2017]_) applied to ``self``.

        The `i`-th K-Bender-Knuth operator swaps the letters `i` and
        `i + 1` everywhere where doing so would not break increasingness.

        EXAMPLES::

            sage: T = IncreasingTableau([[1,3,4],[2,4,5]])
            sage: T.K_bender_knuth(2)
            [[1, 2, 4], [3, 4, 5]]
            sage: T.K_bender_knuth(3)
            [[1, 3, 4], [2, 4, 5]]
        """
        newtab = [[0] * k for k in self.shape()]
        for r, row in enumerate(self):
            for c, val in enumerate(row):
                if val == i:
                    if (c + 1 < len(row) and row[c+1] == i + 1):
                        newtab[r][c] = i
                    elif (r + 1 < len(self) and c < len(self[r+1]) and self[r+1][c] == i + 1):
                        newtab[r][c] = i
                    else:
                        newtab[r][c] = i + 1
                elif val == i + 1:
                    if c > 0 and row[c-1] == i:
                        newtab[r][c] = i + 1
                    elif r > 0 and self[r-1][c] == i:
                        newtab[r][c] = i + 1
                    else:
                        newtab[r][c] = i
                else:
                    newtab[r][c] = val
        return IncreasingTableau(newtab)

    @combinatorial_map(name='K-promotion')
    def K_promotion(self, ceiling=None):
        """
        Return the K-promotion operator from [Pec2014]_ applied to ``self``.

        EXAMPLES::

            sage: T = IncreasingTableau([[1,3,4],[2,4,5]])
            sage: T.K_promotion()
            [[1, 2, 3], [3, 4, 5]]
            sage: T.K_promotion(6)
            [[1, 2, 3], [3, 4, 6]]
            sage: U = IncreasingTableau([[1,3,4],[3,4,5],[5]])
            sage: U.K_promotion()
            [[2, 3, 4], [3, 4, 5], [4]]
        """
        if ceiling is None:
            ceiling = max(self.entries())
        ans = self
        for i in range(1, ceiling):
            ans = ans.K_bender_knuth(i)
        return ans

    @combinatorial_map(name='K-promotion inverse')
    def K_promotion_inverse(self, ceiling=None):
        """
        Return the inverse of K-promotion operator applied to ``self``.

        EXAMPLES::

            sage: T = IncreasingTableau([[1,3,4],[2,4,5]])
            sage: T.K_promotion_inverse()
            [[1, 2, 4], [3, 4, 5]]
            sage: T.K_promotion_inverse(6)
            [[2, 4, 5], [3, 5, 6]]
            sage: U = IncreasingTableau([[1,3,4],[3,4,5],[5]])
            sage: U.K_promotion_inverse()
            [[1, 2, 4], [2, 4, 5], [4]]

        TESTS::

            sage: V = IncreasingTableau([[1,3,4],[3,4,5],[5,6]])
            sage: V == V.K_promotion().K_promotion_inverse()
            True
            sage: V == V.K_promotion_inverse().K_promotion()
            True
        """
        if ceiling is None:
            ceiling = max(self.entries())
        ans = self
        for i in reversed(range(1, ceiling)):
            ans = ans.K_bender_knuth(i)
        return ans

    @combinatorial_map(order=2, name='K-evacuation')
    def K_evacuation(self, ceiling=None):
        """
        Return the K-evacuation involution from [TY2009]_ to ``self``.

        EXAMPLES::

            sage: T = IncreasingTableau([[1,3,4],[2,4,5]])
            sage: T.K_evacuation()
            [[1, 2, 4], [2, 3, 5]]
            sage: T.K_evacuation(6)
            [[2, 3, 5], [3, 4, 6]]
            sage: U = IncreasingTableau([[1,3,4],[3,4,5],[5]])
            sage: U.K_evacuation()
            [[1, 2, 3], [2, 3, 5], [3]]

        TESTS::

            sage: V = IncreasingTableau([[1,3,4],[3,4,5],[5,6]])
            sage: V == V.K_evacuation().K_evacuation()
            True
            sage: V.K_promotion().K_evacuation() == V.K_evacuation().K_promotion_inverse()
            True
        """
        if ceiling is None:
            ceiling = max(self.entries())
        ans = self
        for j in reversed(range(1, ceiling)):
            for i in range(1, j+1):
                ans = ans.K_bender_knuth(i)
        return ans

    @combinatorial_map(order=2, name='dual K-evacuation')
    def dual_K_evacuation(self, ceiling=None):
        """
        Return the dual K-evacuation involution applied to ``self``.

        EXAMPLES::

            sage: T = IncreasingTableau([[1,3,4],[2,4,5]])
            sage: T.dual_K_evacuation()
            [[1, 2, 4], [2, 3, 5]]
            sage: T.dual_K_evacuation(6)
            [[2, 3, 5], [3, 4, 6]]
            sage: U = IncreasingTableau([[1,3,4],[3,4,5],[5]])
            sage: U.dual_K_evacuation()
            [[1, 2, 3], [2, 3, 5], [3]]

        TESTS::

            sage: V = IncreasingTableau([[1,3,4],[3,4,5],[5,6]])
            sage: V == V.dual_K_evacuation().dual_K_evacuation()
            True
            sage: W = IncreasingTableau([[1,2,4],[2,3,5]])
            sage: W.K_evacuation() == W.dual_K_evacuation()
            True
            sage: X = IncreasingTableau([[1,2,4,7],[3,5,6,8],[5,7,8,10],[7,9,10,11]])
            sage: X.K_evacuation() == X.dual_K_evacuation()
            False
            sage: X.K_promotion().dual_K_evacuation() == X.dual_K_evacuation().K_promotion_inverse()
            True
        """
        if ceiling is None:
            ceiling = max(self.entries())
        ans = self
        for j in range(1, ceiling):
            for i in reversed(range(j, ceiling)):
                ans = ans.K_bender_knuth(i)
        return ans


class Tableaux(UniqueRepresentation, Parent):
    """
    A factory class for the various classes of tableaux.

    INPUT:

    - ``n`` -- (optional) nonnegative integer

    OUTPUT:

    - If ``n`` is specified, the class of tableaux of size ``n``. Otherwise,
      the class of all tableaux.

    A tableau in Sage is a finite list of lists, whose lengths are weakly
    decreasing, or an empty list, representing the empty tableau.  The entries
    of a tableau can be any Sage objects. Because of this, no enumeration
    through the set of Tableaux is possible.

    EXAMPLES::

        sage: T = Tableaux(); T
        Tableaux
        sage: T3 = Tableaux(3); T3
        Tableaux of size 3
        sage: [['a','b']] in T
        True
        sage: [['a','b']] in T3
        False
        sage: t = T3([[1,1,1]]); t
        [[1, 1, 1]]
        sage: t in T
        True
        sage: t.parent()
        Tableaux of size 3
        sage: T([]) # the empty tableau
        []
        sage: T.category()
        Category of sets

    .. SEEALSO::

        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`

    TESTS::

        sage: TestSuite( Tableaux() ).run()
        sage: TestSuite( Tableaux(5) ).run()
        sage: t = Tableaux(3)([[1,2],[3]])
        sage: t.parent()
        Tableaux of size 3
        sage: Tableaux(t)
        Traceback (most recent call last):
        ...
        ValueError: the argument to Tableaux() must be a nonnegative integer
        sage: Tableaux(3)([[1, 1]])
        Traceback (most recent call last):
        ...
        ValueError: [[1, 1]] is not an element of Tableaux of size 3

        sage: t0 = Tableau([[1]])
        sage: t1 = Tableaux()([[1]])
        sage: t2 = Tableaux()(t1)
        sage: t0 == t1 == t2
        True
        sage: t1 in Tableaux()
        True
        sage: t1 in Tableaux(1)
        True
        sage: t1 in Tableaux(2)
        False

        sage: [[1]] in Tableaux()
        True
        sage: [] in Tableaux(0)
        True

    Check that :issue:`14145` has been fixed::

        sage: 1 in Tableaux()
        False
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`Tableaux` for more
        information.

        TESTS::

            sage: Tableaux()
            Tableaux
            sage: Tableaux(3)
            Tableaux of size 3
            sage: Tableaux(n=3)
            Tableaux of size 3
        """
        if args:
            n = args[0]
        elif 'n' in kwargs:
            n = kwargs['n']
        else:
            n = None

        if n is None:
            return Tableaux_all()
        else:
            if not isinstance(n, (int, Integer)) or n < 0:
                raise ValueError("the argument to Tableaux() must be a nonnegative integer")
            return Tableaux_size(n)

    Element = Tableau

    # add options to class
    class options(GlobalOptions):
        r"""
        Set the global options for elements of the tableau, skew_tableau,
        and tableau tuple classes. The defaults are for tableau to be
        displayed as a list, latexed as a Young diagram using the English
        convention.

        @OPTIONS@

        .. NOTE::

            Changing the ``convention`` for tableaux also changes the
            ``convention`` for partitions.

        If no parameters are set, then the function returns a copy of the
        options dictionary.

        EXAMPLES::

            sage: T = Tableau([[1,2,3],[4,5]])
            sage: T
            [[1, 2, 3], [4, 5]]
            sage: Tableaux.options.display="array"
            sage: T
              1  2  3
              4  5
            sage: Tableaux.options.convention="french"
            sage: T
              4  5
              1  2  3

        Changing the ``convention`` for tableaux also changes the ``convention``
        for partitions and vice versa::

            sage: P = Partition([3,3,1])
            sage: print(P.ferrers_diagram())
            *
            ***
            ***
            sage: Partitions.options.convention="english"
            sage: print(P.ferrers_diagram())
            ***
            ***
            *
            sage: T
              1  2  3
              4  5

        The ASCII art can also be changed::

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: ascii_art(t)
              1  2  3
              4  5
            sage: Tableaux.options.ascii_art = "table"
            sage: ascii_art(t)
            +---+---+---+
            | 1 | 2 | 3 |
            +---+---+---+
            | 4 | 5 |
            +---+---+
            sage: Tableaux.options.ascii_art = "compact"
            sage: ascii_art(t)
            |1|2|3|
            |4|5|
            sage: Tableaux.options._reset()
        """
        NAME = 'Tableaux'
        module = 'sage.combinat.tableau'
        display = dict(default='list',
                     description='Controls the way in which tableaux are printed',
                     values=dict(list='print tableaux as lists',
                                 diagram='display as Young diagram (similar to :meth:`~sage.combinat.tableau.Tableau.pp()`',
                                 compact='minimal length string representation'),
                     alias=dict(array='diagram', ferrers_diagram='diagram', young_diagram='diagram'),
                     case_sensitive=False)
        ascii_art = dict(default='repr',
                     description='Controls the ascii art output for tableaux',
                     values=dict(repr='display using the diagram string representation',
                                 table='display as a table',
                                 compact='minimal length ascii art'),
                     case_sensitive=False)
        latex = dict(default='diagram',
                   description='Controls the way in which tableaux are latexed',
                   values=dict(list='as a list', diagram='as a Young diagram'),
                   alias=dict(array='diagram', ferrers_diagram='diagram', young_diagram='diagram'),
                   case_sensitive=False)
        convention = dict(default='English',
                          description='Sets the convention used for displaying tableaux and partitions',
                          values=dict(
                              English='use the English convention',
                              French='use the French convention',
                              Russian='use the Russian convention',
                          ),
                          case_sensitive=False)
        notation = dict(alt_name='convention')

    def _element_constructor_(self, t):
        r"""
        Construct an object from ``t`` as an element of ``self``, if
        possible. This is inherited by all Tableaux, SemistandardTableaux, and
        StandardTableaux classes.

        INPUT:

        - ``t`` -- data which can be interpreted as a tableau

        OUTPUT: the corresponding tableau object

        TESTS::

            sage: T = Tableaux(3)
            sage: T([[1,2,1]]).parent() is T     # indirect doctest
            True
            sage: T( StandardTableaux(3)([[1, 2, 3]])).parent() is T
            True
            sage: T([[1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 2]] is not an element of Tableaux of size 3
        """
        if t not in self:
            raise ValueError("%s is not an element of %s" % (t, self))

        return self.element_class(self, t)

    def __contains__(self, x):
        """
        TESTS::

            sage: T = sage.combinat.tableau.Tableaux()
            sage: [[1,2],[3,4]] in T
            True
            sage: [[1,2],[3]] in T
            True
            sage: [] in T
            True
            sage: [['a','b']] in T
            True
            sage: Tableau([['a']]) in T
            True

            sage: [1,2,3] in T
            False
            sage: [[1],[1,2]] in T
            False

        Check that :issue:`14145` is fixed::

            sage: 1 in sage.combinat.tableau.Tableaux()
            False
        """
        from sage.combinat.partition import _Partitions
        if isinstance(x, Tableau):
            return True
        if isinstance(x, list):
            try:
                for row in x:
                    iter(row)
            except TypeError:
                return False
            # any list of lists of partition shape is a tableau
            return [len(row) for row in x] in _Partitions
        else:
            return False


class Tableaux_all(Tableaux):

    def __init__(self):
        r"""
        Initialize the class of all tableaux.

        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_all()
            sage: TestSuite(T).run()
        """
        super().__init__(category=Sets())

    def _repr_(self):
        """
        TESTS::

            sage: repr(Tableaux())    # indirect doctest
            'Tableaux'
        """
        return "Tableaux"

    def an_element(self):
        r"""
        Return a particular element of the class.

        TESTS::

            sage: T = Tableaux()
            sage: T.an_element()
            [[1, 1], [1]]
        """
        return self.element_class(self, [[1, 1], [1]])


class Tableaux_size(Tableaux):
    """
    Tableaux of a fixed size `n`.
    """

    def __init__(self, n):
        r"""
        Initialize the class of tableaux of size `n`.

        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_size(3)
            sage: TestSuite(T).run()

            sage: T = sage.combinat.tableau.Tableaux_size(0)
            sage: TestSuite(T).run()
        """
        super().__init__(category=Sets())
        self.size = n

    def __contains__(self, x):
        """
        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_size(3)
            sage: [[2,4], [1]] in T
            True

            sage: [[2,4],[1,3]] in T
            False

        Check that :issue:`14145` is fixed::

            sage: 1 in sage.combinat.tableau.Tableaux_size(3)
            False
        """
        return Tableaux.__contains__(self, x) and sum(len(row) for row in x) == self.size

    def _repr_(self):
        """
        TESTS::

            sage: repr(Tableaux(4))    # indirect doctest
            'Tableaux of size 4'
        """
        return "Tableaux of size %s" % self.size

    def an_element(self):
        r"""
        Return a particular element of the class.

        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_size(3)
            sage: T.an_element()
            [[1, 1], [1]]
            sage: T = sage.combinat.tableau.Tableaux_size(0)
            sage: T.an_element()
            []
        """
        if self.size == 0:
            return self.element_class(self, [])

        if self.size == 1:
            return self.element_class(self, [[1]])

        return self.element_class(self, [[1]*(self.size-1), [1]])


##########################
# Semi-standard tableaux #
##########################
class SemistandardTableaux(Tableaux):
    """
    A factory class for the various classes of semistandard tableaux.

    INPUT:

    Keyword arguments:

    - ``size`` -- the size of the tableaux
    - ``shape`` -- the shape of the tableaux
    - ``eval`` -- the weight (also called content or evaluation) of
      the tableaux
    - ``max_entry`` -- a maximum entry for the tableaux.  This can be a
      positive integer or infinity (``oo``). If ``size`` or ``shape`` are
      specified, ``max_entry`` defaults to be ``size`` or the size of
      ``shape``.

    Positional arguments:

    - The first argument is interpreted as either ``size`` or ``shape``
      according to whether it is an integer or a partition
    - The second keyword argument will always be interpreted as ``eval``

    OUTPUT:

    - The appropriate class, after checking basic consistency tests. (For
      example, specifying ``eval`` implies a value for `max_entry`).

    A semistandard tableau is a tableau whose entries are positive integers,
    which are weakly increasing in rows and strictly increasing down columns.
    Note that Sage uses the English convention for partitions and tableaux;
    the longer rows are displayed on top.

    Classes of semistandard tableaux can be iterated over if and only if there
    is some restriction.

    EXAMPLES::

        sage: SST = SemistandardTableaux([2,1]); SST
        Semistandard tableaux of shape [2, 1] and maximum entry 3
        sage: SST.list()                                                                # needs sage.modules
        [[[1, 1], [2]],
         [[1, 1], [3]],
         [[1, 2], [2]],
         [[1, 2], [3]],
         [[1, 3], [2]],
         [[1, 3], [3]],
         [[2, 2], [3]],
         [[2, 3], [3]]]

        sage: SST = SemistandardTableaux(3); SST
        Semistandard tableaux of size 3 and maximum entry 3
        sage: SST.list()                                                                # needs sage.modules
        [[[1, 1, 1]],
         [[1, 1, 2]],
         [[1, 1, 3]],
         [[1, 2, 2]],
         [[1, 2, 3]],
         [[1, 3, 3]],
         [[2, 2, 2]],
         [[2, 2, 3]],
         [[2, 3, 3]],
         [[3, 3, 3]],
         [[1, 1], [2]],
         [[1, 1], [3]],
         [[1, 2], [2]],
         [[1, 2], [3]],
         [[1, 3], [2]],
         [[1, 3], [3]],
         [[2, 2], [3]],
         [[2, 3], [3]],
         [[1], [2], [3]]]

        sage: SST = SemistandardTableaux(3, max_entry=2); SST
        Semistandard tableaux of size 3 and maximum entry 2
        sage: SST.list()                                                                # needs sage.modules
        [[[1, 1, 1]],
         [[1, 1, 2]],
         [[1, 2, 2]],
         [[2, 2, 2]],
         [[1, 1], [2]],
         [[1, 2], [2]]]

        sage: SST = SemistandardTableaux(3, max_entry=oo); SST
        Semistandard tableaux of size 3
        sage: SST[123]                                                                  # needs sage.modules
        [[3, 4], [6]]

        sage: SemistandardTableaux(max_entry=2)[11]                                     # needs sage.modules
        [[1, 1], [2]]

        sage: SemistandardTableaux()[0]                                                 # needs sage.modules
        []

    .. SEEALSO::

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`SemistandardTableaux`
        for more information.

        TESTS::

            sage: SemistandardTableaux()
            Semistandard tableaux
            sage: SemistandardTableaux(3)
            Semistandard tableaux of size 3 and maximum entry 3
            sage: SemistandardTableaux(size=3)
            Semistandard tableaux of size 3 and maximum entry 3
            sage: SemistandardTableaux(0)
            Semistandard tableaux of size 0 and maximum entry 0
            sage: SemistandardTableaux([2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux(shape=[2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux([])
            Semistandard tableaux of shape [] and maximum entry 0
            sage: SemistandardTableaux(eval=[2,1])
            Semistandard tableaux of size 3 and weight [2, 1]
            sage: SemistandardTableaux(max_entry=3)
            Semistandard tableaux with maximum entry 3
            sage: SemistandardTableaux(3, [2,1])
            Semistandard tableaux of size 3 and weight [2, 1]
            sage: SemistandardTableaux(3, shape=[2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux(3, [2,1], shape=[2,1])
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: SemistandardTableaux(3, max_entry=4)
            Semistandard tableaux of size 3 and maximum entry 4
            sage: SemistandardTableaux(3, max_entry=oo)
            Semistandard tableaux of size 3
            sage: SemistandardTableaux([2, 1], max_entry=oo)
            Semistandard tableaux of shape [2, 1]
            sage: SemistandardTableaux([2, 1], [2, 1])
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: mu = Partition([2,1]); SemistandardTableaux(mu, mu)
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: SemistandardTableaux(3, [2, 1], max_entry=2)
            Semistandard tableaux of size 3 and weight [2, 1]

            sage: SemistandardTableaux(3, shape=[2])
            Traceback (most recent call last):
            ...
            ValueError: size and shape are different sizes

            sage: SemistandardTableaux(3, [2])
            Traceback (most recent call last):
            ...
            ValueError: size and eval are different sizes

            sage: SemistandardTableaux([2],[3])
            Traceback (most recent call last):
            ...
            ValueError: shape and eval are different sizes

            sage: SemistandardTableaux(2,[2], max_entry=4)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must match the weight

            sage: SemistandardTableaux(eval=[2], max_entry=oo)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must match the weight

            sage: # needs sage.modules
            sage: SemistandardTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: shape must be a (skew) partition
        """
        from sage.combinat.partition import Partition, _Partitions
        # Process the keyword arguments -- allow for original syntax where
        #   n == size,  p== shape and mu == eval
        n = kwargs.get('n', None)
        size = kwargs.get('size', n)

        p = kwargs.get('p', None)
        shape = kwargs.get('shape', p)

        mu = kwargs.get('eval', None)
        mu = kwargs.get("mu", mu)

        max_entry = kwargs.get('max_entry', None)

        # Process the positional arguments
        if args:
            # The first arg could be either a size or a shape
            if isinstance(args[0], (int, Integer)):
                if size is not None:
                    raise ValueError("size was specified more than once")
                else:
                    size = args[0]
            else:
                if shape is not None:
                    raise ValueError("the shape was specified more than once")
                shape = args[0]  # we check it's a partition later

        if len(args) == 2:
            # The second non-keyword argument is the weight
            if mu is not None:
                raise ValueError("the weight was specified more than once")
            else:
                mu = args[1]

        # Consistency checks
        if size is not None:
            if not isinstance(size, (int, Integer)):
                raise ValueError("size must be an integer")
            if size < 0:
                raise ValueError("size must be nonnegative")

        if shape is not None:
            from sage.combinat.skew_partition import SkewPartitions
            # use in (and not isinstance) below so that lists can be used as
            # shorthand
            if shape in _Partitions:
                shape = Partition(shape)
            elif shape in SkewPartitions():
                from sage.combinat.skew_tableau import SemistandardSkewTableaux
                return SemistandardSkewTableaux(shape, mu)
            else:
                raise ValueError("shape must be a (skew) partition")

        if mu is not None:
            if mu not in Compositions() and mu not in _Partitions:
                raise ValueError("mu must be a composition")
            mu = Composition(mu)

        is_inf = max_entry is PlusInfinity()

        if max_entry is not None:
            if not is_inf and not isinstance(max_entry, (int, Integer)):
                raise ValueError("max_entry must be an integer or PlusInfinity")
            elif max_entry <= 0:
                raise ValueError("max_entry must be positive")

        if (mu is not None) and (max_entry is not None):
            if max_entry != len(mu):
                raise ValueError("the maximum entry must match the weight")

        if (size is not None) and (shape is not None):
            if sum(shape) != size:
                # This could return an empty class instead of an error
                raise ValueError("size and shape are different sizes")

        if (size is not None) and (mu is not None):
            if sum(mu) != size:
                # This could return an empty class instead of an error
                raise ValueError("size and eval are different sizes")

        # Dispatch appropriately
        if (shape is not None) and (mu is not None):
            if sum(shape) != sum(mu):
                # This could return an empty class instead of an error
                raise ValueError("shape and eval are different sizes")
            else:
                return SemistandardTableaux_shape_weight(shape, mu)

        if (shape is not None):
            if is_inf:
                return SemistandardTableaux_shape_inf(shape)
            return SemistandardTableaux_shape(shape, max_entry)

        if (mu is not None):
            return SemistandardTableaux_size_weight(sum(mu), mu)

        if (size is not None):
            if is_inf:
                return SemistandardTableaux_size_inf(size)
            return SemistandardTableaux_size(size, max_entry)

        return SemistandardTableaux_all(max_entry)

    Element = SemistandardTableau

    def __init__(self, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SemistandardTableaux()
            sage: TestSuite(S).run()                                                    # needs sage.modules
        """
        if 'max_entry' in kwds:
            self.max_entry = kwds['max_entry']
            kwds.pop('max_entry')
        else:
            self.max_entry = None
        Tableaux.__init__(self, **kwds)

    def __getitem__(self, r):
        r"""
        The default implementation of ``__getitem__`` for enumerated sets
        does not allow slices so we override it.

        EXAMPLES::

            sage: StandardTableaux([4,3,3,2])[10:20]     # indirect doctest
            [[[1, 3, 9, 12], [2, 5, 10], [4, 6, 11], [7, 8]],
             [[1, 2, 9, 12], [3, 5, 10], [4, 6, 11], [7, 8]],
             [[1, 3, 9, 12], [2, 4, 10], [5, 6, 11], [7, 8]],
             [[1, 2, 9, 12], [3, 4, 10], [5, 6, 11], [7, 8]],
             [[1, 5, 8, 12], [2, 6, 10], [3, 7, 11], [4, 9]],
             [[1, 4, 8, 12], [2, 6, 10], [3, 7, 11], [5, 9]],
             [[1, 3, 8, 12], [2, 6, 10], [4, 7, 11], [5, 9]],
             [[1, 2, 8, 12], [3, 6, 10], [4, 7, 11], [5, 9]],
             [[1, 4, 8, 12], [2, 5, 10], [3, 7, 11], [6, 9]],
             [[1, 3, 8, 12], [2, 5, 10], [4, 7, 11], [6, 9]]]

            sage: SemistandardTableaux(size=2, max_entry=oo)[5]                         # needs sage.modules
            [[2, 3]]

            sage: SemistandardTableaux([2,1], max_entry=oo)[3]                          # needs sage.modules
            [[1, 2], [3]]

            sage: SemistandardTableaux(3, max_entry=2)[0:5]    # indirect doctest       # needs sage.modules
            [[[1, 1, 1]],
            [[1, 1, 2]],
            [[1, 2, 2]],
            [[2, 2, 2]],
            [[1, 1], [2]]]

            sage: SemistandardTableaux([2,2], [2, 1, 1])[0]    # indirect doctest       # needs sage.modules
            [[1, 1], [2, 3]]

            sage: SemistandardTableaux([1,1,1], max_entry=4)[0:4]                       # needs sage.modules
            [[[1], [2], [3]],
             [[1], [2], [4]],
             [[1], [3], [4]],
             [[2], [3], [4]]]

            sage: SemistandardTableaux(3, [2,1])[1]    # indirect doctest               # needs sage.modules
            [[1, 1], [2]]

            sage: StandardTableaux(3)[:]  # indirect doctest                            # needs sage.modules
            [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]

            sage: StandardTableaux([2,2])[1]   # indirect doctest                       # needs sage.modules
            [[1, 2], [3, 4]]

        TESTS::

            sage: SemistandardTableaux()[5]                                             # needs sage.modules
            [[1], [2]]

            sage: SemistandardTableaux(max_entry=2)[5]                                  # needs sage.modules
            [[2, 2]]

            sage: SemistandardTableaux()[:]                                             # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: infinite set

            sage: SemistandardTableaux(size=2, max_entry=oo)[:]                         # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: infinite set
        """
        if isinstance(r, (int, Integer)):
            return self.unrank(r)
        elif isinstance(r, slice):
            start = 0 if r.start is None else r.start
            stop = r.stop
            if stop is None and not self.is_finite():
                raise ValueError('infinite set')
        else:
            raise ValueError('r must be an integer or a slice')
        count = 0
        tabs = []
        for t in self:
            if count == stop:
                break
            if count >= start:
                tabs.append(t)
            count += 1

        # this is to cope with empty slices endpoints like [:6] or [:}
        if count == stop or stop is None:
            return tabs
        raise IndexError('value out of range')

    def __contains__(self, t):
        """
        Return ``True`` if ``t`` can be interpreted as a
        :class:`SemistandardTableau`.

        TESTS::

            sage: T = sage.combinat.tableau.SemistandardTableaux_all()
            sage: [[1,2],[2]] in T
            True
            sage: [] in T
            True
            sage: Tableau([[1]]) in T
            True
            sage: StandardTableau([[1]]) in T
            True

            sage: [[1,2],[1]] in T
            False
            sage: [[1,1],[5]] in T
            True
            sage: [[1,3,2]] in T
            False

        Check that :issue:`14145` is fixed::

            sage: 1 in sage.combinat.tableau.SemistandardTableaux()
            False
        """
        if isinstance(t, SemistandardTableau):
            return (self.max_entry is None or
                    len(t) == 0 or
                    max(max(row) for row in t) <= self.max_entry)
        if not t:
            return True
        if Tableaux.__contains__(self, t):
            for row in t:
                if not all(c > 0 for c in row):
                    return False
                if not all(row[i] <= row[i+1] for i in range(len(row)-1)):
                    return False
            for row, next in zip(t, t[1:]):
                if not all(row[c] < next[c] for c in range(len(next))):
                    return False
            return self.max_entry is None or max(max(row) for row in t) <= self.max_entry
        else:
            return False


class SemistandardTableaux_all(SemistandardTableaux, DisjointUnionEnumeratedSets):
    """
    All semistandard tableaux.
    """

    def __init__(self, max_entry=None):
        r"""
        Initialize the class of all semistandard tableaux.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: T = sage.combinat.tableau.SemistandardTableaux_all()
            sage: TestSuite(T).run()                                                    # needs sage.modules

            sage: T = sage.combinat.tableau.SemistandardTableaux_all(max_entry=3)
            sage: TestSuite(T).run()            # long time                             # needs sage.modules
        """
        if max_entry is not PlusInfinity():
            self.max_entry = max_entry

            def SST_n(n):
                return SemistandardTableaux_size(n, max_entry)
            DisjointUnionEnumeratedSets.__init__(self,
                Family(NonNegativeIntegers(), SST_n),
                facade=True, keepkey=False)
        else:
            self.max_entry = None

    def _repr_(self):
        """
        TESTS::

            sage: SemistandardTableaux()    # indirect doctest
            Semistandard tableaux

            sage: SemistandardTableaux(max_entry=3)
            Semistandard tableaux with maximum entry 3
        """
        if self.max_entry is not None:
            return "Semistandard tableaux with maximum entry %s" % str(self.max_entry)
        return "Semistandard tableaux"

    def list(self):
        """
        TESTS::

            sage: SemistandardTableaux().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class SemistandardTableaux_size_inf(SemistandardTableaux):
    """
    Semistandard tableaux of fixed size `n` with no maximum entry.
    """

    def __init__(self, n):
        r"""
        Initialize the class of semistandard tableaux of size ``n`` with no
        maximum entry.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: T = sage.combinat.tableau.SemistandardTableaux_size_inf(3)
            sage: TestSuite(T).run()                                                    # needs sage.modules
        """
        super().__init__(category=InfiniteEnumeratedSets())
        self.size = n

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux(3, max_entry=oo))    # indirect doctest
            'Semistandard tableaux of size 3'
        """
        return "Semistandard tableaux of size %s" % str(self.size)

    def __contains__(self, t):
        """
        Return ``True`` if ``t`` can be interpreted as an element of this
        class.

        TESTS::

            sage: T = SemistandardTableaux(3, max_entry=oo)
            sage: [[1,2],[5]] in T
            True
            sage: StandardTableau([[1, 2], [3]]) in T
            True

            sage: [] in T
            False
            sage: Tableau([[1]]) in T
            False

        Check that :issue:`14145` is fixed::

            sage: 1 in SemistandardTableaux(3, max_entry=oo)
            False
        """
        return SemistandardTableaux.__contains__(self, t) and sum(map(len, t)) == self.size

    def __iter__(self):
        """
        EXAMPLES::

            sage: sst = SemistandardTableaux(3, max_entry=oo)
            sage: [sst[t] for t in range(5)]                                            # needs sage.modules
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 2, 2]],
             [[2, 2, 2]],
             [[1, 1], [2]]]
            sage: sst[1000]                                                             # needs sage.modules
            [[2, 12], [7]]
            sage: sst[0].parent() is sst                                                # needs sage.modules
            True
        """
        from sage.combinat.partition import Partitions
        # Iterates through with maximum entry as order
        i = 1
        while True:
            for part in Partitions(self.size):
                if i != 1:
                    for k in range(1, self.size+1):
                        for c in integer_vectors_nk_fast_iter(self.size - k, i-1):
                            c.append(k)
                            for sst in SemistandardTableaux_shape_weight(part, Composition(c)):
                                yield self.element_class(self, sst)
                else:
                    for sst in SemistandardTableaux_shape_weight(part, Composition([self.size])):
                        yield self.element_class(self, sst)
            i += 1

    def list(self):
        """
        TESTS::

            sage: SemistandardTableaux(3, max_entry=oo).list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class SemistandardTableaux_shape_inf(SemistandardTableaux):
    """
    Semistandard tableaux of fixed shape `p` and no maximum entry.
    """

    def __init__(self, p):
        r"""
        Initialize the class of semistandard tableaux of shape ``p`` and no
        maximum entry.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: SST = SemistandardTableaux([2,1], max_entry=oo)
            sage: type(SST)
            <class 'sage.combinat.tableau.SemistandardTableaux_shape_inf_with_category'>
            sage: TestSuite(SST).run()                                                  # needs sage.modules
        """
        super().__init__(category=InfiniteEnumeratedSets())
        self.shape = p

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: SST = SemistandardTableaux([2,1], max_entry=oo)
            sage: [[13, 67], [1467]] in SST
            True
            sage: SST = SemistandardTableaux([3,1], max_entry=oo)
            sage: [[13, 67], [1467]] in SST
            False

        Check that :issue:`14145` is fixed::

            sage: SST = SemistandardTableaux([3,1], max_entry=oo)
            sage: 1 in SST
            False
        """
        return SemistandardTableaux.__contains__(self, x) and [len(_) for _ in x] == self.shape

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux([2,1], max_entry=oo))    # indirect doctest
            'Semistandard tableaux of shape [2, 1]'
        """
        return "Semistandard tableaux of shape %s" % str(self.shape)

    def __iter__(self):
        """
        An iterator for the semistandard partitions of shape ``p`` and no
        maximum entry. Iterates through with maximum entry as order.

        EXAMPLES::

            sage: SST = SemistandardTableaux([3, 1], max_entry=oo)
            sage: SST[1000]                                                             # needs sage.modules
            [[1, 1, 10], [6]]
            sage: [ SST[t] for t in range(5) ]                                          # needs sage.modules
            [[[1, 1, 1], [2]],
             [[1, 1, 2], [2]],
             [[1, 2, 2], [2]],
             [[1, 1, 1], [3]],
             [[1, 1, 2], [3]]]
            sage: SST[0].parent() is SST                                                # needs sage.modules
            True
        """
        # Iterates through with maximum entry as order
        i = 1
        n = sum(self.shape)
        while True:
            if i != 1:
                for k in range(1, n+1):
                    for c in integer_vectors_nk_fast_iter(n - k, i-1):
                        c.append(k)
                        for sst in SemistandardTableaux_shape_weight(self.shape, Composition(c)):
                            yield self.element_class(self, sst)
            else:
                for sst in SemistandardTableaux_shape_weight(self.shape, Composition([n])):
                    yield self.element_class(self, sst)
            i += 1


class SemistandardTableaux_size(SemistandardTableaux):
    """
    Semistandard tableaux of fixed size `n`.
    """

    def __init__(self, n, max_entry=None):
        r"""
        Initialize the class of semistandard tableaux of size `n`.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardTableaux`
            to ensure the options are properly parsed.

        TESTS::

            sage: SST = SemistandardTableaux(3); SST
            Semistandard tableaux of size 3 and maximum entry 3
            sage: type(SST)
            <class 'sage.combinat.tableau.SemistandardTableaux_size_with_category'>
            sage: TestSuite(SST).run()                                                  # needs sage.modules

            sage: SST = SemistandardTableaux(3, max_entry=6)
            sage: type(SST)
            <class 'sage.combinat.tableau.SemistandardTableaux_size_with_category'>
            sage: TestSuite(SST).run()                                                  # needs sage.modules
        """

        if max_entry is None:
            max_entry = n
        super().__init__(max_entry=max_entry,
                         category=FiniteEnumeratedSets())
        self.size = n

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux(3))    # indirect doctest
            'Semistandard tableaux of size 3 and maximum entry 3'

            sage: repr(SemistandardTableaux(3, max_entry=6))
            'Semistandard tableaux of size 3 and maximum entry 6'
        """
        return "Semistandard tableaux of size %s and maximum entry %s" % (str(self.size), str(self.max_entry))

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[1,2],[3,3]] in SemistandardTableaux(3)
            False
            sage: [[1,2],[3,3]] in SemistandardTableaux(4)
            True
            sage: [[1,2],[3,3]] in SemistandardTableaux(4, max_entry=2)
            False
            sage: SST = SemistandardTableaux(4)
            sage: all(sst in SST for sst in SST)                                        # needs sage.modules
            True

        Check that :issue:`14145` is fixed::

            sage: SST = SemistandardTableaux(4)
            sage: 1 in SST
            False
        """
        if self.size == 0:
            return x == []

        return (SemistandardTableaux.__contains__(self, x)
            and sum(map(len, x)) == self.size
            and max(max(row) for row in x) <= self.max_entry)

    def random_element(self):
        r"""
        Generate a random :class:`SemistandardTableau` with uniform probability.

        The RSK algorithm gives a bijection between symmetric `k\times k` matrices
        of nonnegative integers that sum to `n` and semistandard tableaux with size `n`
        and maximum entry `k`.

        The number of `k\times k` symmetric matrices of nonnegative integers
        having sum of elements on the diagonal `i` and sum of elements above
        the diagonal `j` is `\binom{k + i - 1}{k - 1}\binom{\binom{k}{2} + j - 1}{\binom{k}{2} - 1}`.
        We first choose the sum of the elements on the diagonal randomly weighted by the
        number of matrices having that trace.  We then create random integer vectors
        of length `k` having that sum and use them to generate a `k\times k` diagonal matrix.
        Then we take a random integer vector of length `\binom{k}{2}` summing to half the
        remainder and distribute it symmetrically to the remainder of the matrix.

        Applying RSK to the random symmetric matrix gives us a pair of identical
        :class:`SemistandardTableau` of which we choose the first.

        EXAMPLES::

            sage: SemistandardTableaux(6).random_element()  # random                    # needs sage.modules
            [[1, 1, 2], [3, 5, 5]]
            sage: SemistandardTableaux(6, max_entry=7).random_element()  # random       # needs sage.modules
            [[2, 4, 4, 6, 6, 6]]
        """
        from sage.rings.integer_ring import ZZ
        from sage.matrix.constructor import diagonal_matrix
        from sage.combinat.rsk import RSK
        kchoose2m1 = self.max_entry * (self.max_entry - 1) // 2 - 1
        km1 = self.max_entry - 1
        weights = [binomial(self.size - i + km1, km1) * binomial((i//2) + kchoose2m1, kchoose2m1)
                   for i in range(0, self.size + 1, 2)]
        randpos = ZZ.random_element(sum(weights))
        tot = weights[0]
        pos = 0
        while randpos >= tot:
            pos += 1
            tot += weights[pos]
        # we now have pos elements over the diagonal and n - 2 * pos on it
        m = diagonal_matrix(list(IntegerVectors(self.size - 2 * pos,
                                                self.max_entry).random_element()))
        above_diagonal = list(IntegerVectors(pos, kchoose2m1 + 1).random_element())
        index = 0
        for i in range(self.max_entry - 1):
            for j in range(i + 1, self.max_entry):
                m[i, j] = above_diagonal[index]
                m[j, i] = above_diagonal[index]
                index += 1
        return RSK(m)[0]

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: SemistandardTableaux(3).cardinality()
            19
            sage: SemistandardTableaux(4).cardinality()
            116
            sage: SemistandardTableaux(4, max_entry=2).cardinality()
            9
            sage: SemistandardTableaux(4, max_entry=10).cardinality()
            4225
            sage: ns = list(range(1, 6))
            sage: ssts = [ SemistandardTableaux(n) for n in ns ]
            sage: all(sst.cardinality() == len(sst.list()) for sst in ssts)             # needs sage.modules
            True
        """
        from sage.combinat.partition import Partitions
        c = 0
        for part in Partitions(self.size):
            c += SemistandardTableaux_shape(part, self.max_entry).cardinality()
        return c

    def __iter__(self):
        """
        EXAMPLES::

            sage: [ t for t in SemistandardTableaux(2) ]                                # needs sage.modules
            [[[1, 1]], [[1, 2]], [[2, 2]], [[1], [2]]]
            sage: [ t for t in SemistandardTableaux(3) ]                                # needs sage.modules
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 1, 3]],
             [[1, 2, 2]],
             [[1, 2, 3]],
             [[1, 3, 3]],
             [[2, 2, 2]],
             [[2, 2, 3]],
             [[2, 3, 3]],
             [[3, 3, 3]],
             [[1, 1], [2]],
             [[1, 1], [3]],
             [[1, 2], [2]],
             [[1, 2], [3]],
             [[1, 3], [2]],
             [[1, 3], [3]],
             [[2, 2], [3]],
             [[2, 3], [3]],
             [[1], [2], [3]]]

            sage: [ t for t in SemistandardTableaux(3, max_entry=2) ]                   # needs sage.modules
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 2, 2]],
             [[2, 2, 2]],
             [[1, 1], [2]],
             [[1, 2], [2]]]

            sage: sst = SemistandardTableaux(3)
            sage: sst[0].parent() is sst                                                # needs sage.modules
            True
        """
        from sage.combinat.partition import Partitions
        for part in Partitions(self.size):
            for sst in SemistandardTableaux_shape(part, self.max_entry):
                yield self.element_class(self, sst)


class SemistandardTableaux_shape(SemistandardTableaux):
    """
    Semistandard tableaux of fixed shape `p` with a given max entry.

    A semistandard tableau with max entry `i` is required to have all
    its entries less or equal to `i`. It is not required to actually
    contain an entry `i`.

    INPUT:

    - ``p`` -- a partition
    - ``max_entry`` -- the max entry; defaults to the size of `p`
    """

    def __init__(self, p, max_entry=None):
        r"""
        Initialize the class of semistandard tableaux of shape `p`, with a
        given ``max_entry``.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: SST = SemistandardTableaux([2,1])
            sage: TestSuite(SST).run()                                                  # needs sage.modules

            sage: SST = SemistandardTableaux([2,1], max_entry=5)
            sage: TestSuite(SST).run()                                                  # needs sage.modules
        """
        if max_entry is None:
            max_entry = sum(p)
        super().__init__(max_entry=max_entry,
                         category=FiniteEnumeratedSets())
        self.shape = p

    def __iter__(self):
        """
        An iterator for the semistandard tableaux of the specified shape
        with the specified max entry.

        EXAMPLES::

            sage: [ t for t in SemistandardTableaux([3]) ]                              # needs sage.modules
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 1, 3]],
             [[1, 2, 2]],
             [[1, 2, 3]],
             [[1, 3, 3]],
             [[2, 2, 2]],
             [[2, 2, 3]],
             [[2, 3, 3]],
             [[3, 3, 3]]]
            sage: [ t for t in SemistandardTableaux([2,1]) ]                            # needs sage.modules
            [[[1, 1], [2]],
             [[1, 1], [3]],
             [[1, 2], [2]],
             [[1, 2], [3]],
             [[1, 3], [2]],
             [[1, 3], [3]],
             [[2, 2], [3]],
             [[2, 3], [3]]]
            sage: [ t for t in SemistandardTableaux([1,1,1]) ]                          # needs sage.modules
            [[[1], [2], [3]]]

            sage: [ t for t in SemistandardTableaux([1,1,1], max_entry=4) ]             # needs sage.modules
            [[[1], [2], [3]],
             [[1], [2], [4]],
             [[1], [3], [4]],
             [[2], [3], [4]]]

            sage: sst = SemistandardTableaux([3])
            sage: sst[0].parent() is sst                                                # needs sage.modules
            True
        """
        for c in integer_vectors_nk_fast_iter(sum(self.shape), self.max_entry):
            for sst in SemistandardTableaux_shape_weight(self.shape, Composition(c)):
                yield self.element_class(self, sst)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: SST = SemistandardTableaux([2,1])
            sage: all(sst in SST for sst in SST)                                        # needs sage.modules
            True
            sage: len([x for x in SemistandardTableaux(3) if x in SST])                 # needs sage.modules
            8
            sage: SST.cardinality()
            8

            sage: SST = SemistandardTableaux([2,1], max_entry=4)
            sage: all(sst in SST for sst in SST)                                        # needs sage.modules
            True
            sage: SST.cardinality()
            20
        """
        return SemistandardTableaux.__contains__(self, x) and [len(_) for _ in x] == self.shape

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux([2,1]))    # indirect doctest
            'Semistandard tableaux of shape [2, 1] and maximum entry 3'

            sage: repr(SemistandardTableaux([2,1], max_entry=5))
            'Semistandard tableaux of shape [2, 1] and maximum entry 5'
        """
        return "Semistandard tableaux of shape %s and maximum entry %s" % (str(self.shape), str(self.max_entry))

    def random_element(self):
        """
        Return a uniformly distributed random tableau of the given ``shape`` and ``max_entry``.

        Uses the algorithm from [Kra1999]_ based on the Novelli-Pak-Stoyanovskii bijection
           http://www.sciencedirect.com/science/article/pii/0012365X9290368P

        EXAMPLES::

            sage: S = SemistandardTableaux([2, 2, 1, 1])
            sage: S.random_element() in S
            True
            sage: S = SemistandardTableaux([2, 2, 1, 1], max_entry=7)
            sage: S.random_element() in S
            True
        """
        from sage.misc.prandom import randint
        with_sentinels = [max(i, j) for i, j in zip([0]+list(self.shape), [k+1 for k in self.shape]+[0])]
        t = [[self.max_entry+1]*i for i in with_sentinels]
        for i, l in enumerate(self.shape):
            for j in range(l):
                content = j - i
                t[i][j] = randint(1 - content, self.max_entry)
        conj = self.shape.conjugate()
        for i in range(len(conj) - 1, -1, -1):
            for j in range(conj[i] - 1, -1, -1):
                row = j
                col = i
                s = t[row][col]
                x = t[row][col + 1]
                y = t[row + 1][col]
                while s > x or s >= y:
                    if x + 1 < y:
                        t[row][col] = x + 1
                        t[row][col + 1] = s
                        col += 1
                    else:
                        t[row][col] = y - 1
                        t[row + 1][col] = s
                        row += 1
                    x = t[row][col + 1]
                    y = t[row + 1][col]
        return SemistandardTableau([l[:c] for l, c in zip(t, self.shape)])

    def cardinality(self, algorithm='hook'):
        r"""
        Return the cardinality of ``self``.

        INPUT:

        - ``algorithm`` -- (default: ``'hook'``) any one of the following:

          - ``'hook'`` -- use Stanley's hook length formula

          - ``'sum'`` -- sum over the compositions of ``max_entry`` the
            number of semistandard tableau with ``shape`` and given
            weight vector

        This is computed using *Stanley's hook length formula*:

        .. MATH::

           f_{\lambda} = \prod_{u\in\lambda} \frac{n+c(u)}{h(u)}.

        where `n` is the ``max_entry``, `c(u)` is the content of `u`,
        and `h(u)` is the hook length of `u`.
        See [Sta-EC2]_ Corollary 7.21.4.

        EXAMPLES::

            sage: SemistandardTableaux([2,1]).cardinality()
            8
            sage: SemistandardTableaux([2,2,1]).cardinality()
            75
            sage: SymmetricFunctions(QQ).schur()([2,2,1]).expand(5)(1,1,1,1,1)  # cross check       # needs sage.modules
            75
            sage: SemistandardTableaux([5]).cardinality()
            126
            sage: SemistandardTableaux([3,2,1]).cardinality()
            896
            sage: SemistandardTableaux([3,2,1], max_entry=7).cardinality()
            2352
            sage: SemistandardTableaux([6,5,4,3,2,1], max_entry=30).cardinality()
            208361017592001331200
            sage: ssts = [SemistandardTableaux(p, max_entry=6) for p in Partitions(5)]
            sage: all(sst.cardinality() == sst.cardinality(algorithm='sum')             # needs sage.modules
            ....:     for sst in ssts)
            True
        """
        if algorithm == 'hook':
            conj = self.shape.conjugate()
            num = Integer(1)
            den = Integer(1)
            for i, l in enumerate(self.shape):
                for j in range(l):
                    num *= self.max_entry + j - i
                    den *= l + conj[j] - i - j - 1
            return Integer(num / den)
        elif algorithm == 'sum':
            c = 0
            for comp in integer_vectors_nk_fast_iter(sum(self.shape), self.max_entry):
                c += SemistandardTableaux_shape_weight(self.shape, Composition(comp)).cardinality()
            return c
        raise ValueError("unknown algorithm {}".format(algorithm))


class SemistandardTableaux_shape_weight(SemistandardTableaux_shape):
    r"""
    Semistandard tableaux of fixed shape `p` and weight `\mu`.
    """

    def __init__(self, p, mu):
        r"""
        Initialize the class of all semistandard tableaux of shape ``p`` and
        weight ``mu``.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: SST = SemistandardTableaux([2,1], [2,1])
            sage: TestSuite(SST).run()                                                  # needs sage.modules
        """
        super().__init__(p, len(mu))
        self.weight = mu

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux([2,1],[2,1]))    # indirect doctest
            'Semistandard tableaux of shape [2, 1] and weight [2, 1]'
        """
        return "Semistandard tableaux of shape %s and weight %s" % (self.shape, self.weight)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: SST = SemistandardTableaux([2,1], [2,1])
            sage: all(sst in SST for sst in SST)                                        # needs sage.modules
            True
            sage: len([x for x in SemistandardTableaux(3) if x in SST])                 # needs sage.modules
            1
            sage: SST.cardinality()                                                     # needs sage.modules
            1
        """
        if x not in SemistandardTableaux_shape(self.shape, self.max_entry):
            return False
        n = sum(self.shape)

        if n == 0 and len(x) == 0:
            return True

        content = {}
        for row in x:
            for i in row:
                content[i] = content.get(i, 0) + 1
        content_list = [0]*int(max(content))

        for key in content:
            content_list[key-1] = content[key]

        if content_list != self.weight:
            return False

        return True

    def cardinality(self):
        """
        Return the number of semistandard tableaux of the given shape and
        weight, as computed by ``kostka_number`` function of ``symmetrica``.

        EXAMPLES::

            sage: # needs sage.modules
            sage: SemistandardTableaux([2,2], [2, 1, 1]).cardinality()
            1
            sage: SemistandardTableaux([2,2,2], [2, 2, 1,1]).cardinality()
            1
            sage: SemistandardTableaux([2,2,2], [2, 2, 2]).cardinality()
            1
            sage: SemistandardTableaux([3,2,1], [2, 2, 2]).cardinality()
            2
        """
        return symmetrica.kostka_number(self.shape, self.weight)

    def __iter__(self):
        """
        TESTS::

            sage: sst = SemistandardTableaux([3,1],[2,1,1])
            sage: [sst[i] for i in range(2)]                                            # needs sage.modules
            [[[1, 1, 2], [3]], [[1, 1, 3], [2]]]
            sage: sst[0].parent() is sst                                                # needs sage.modules
            True
        """
        for t in symmetrica.kostka_tab(self.shape, self.weight):
            yield self.element_class(self, t)

    def list(self):
        """
        Return a list of all semistandard tableaux in ``self`` generated
        by symmetrica.

        EXAMPLES::

            sage: # needs sage.modules
            sage: SemistandardTableaux([2,2], [2, 1, 1]).list()
            [[[1, 1], [2, 3]]]
            sage: SemistandardTableaux([2,2,2], [2, 2, 1,1]).list()
            [[[1, 1], [2, 2], [3, 4]]]
            sage: SemistandardTableaux([2,2,2], [2, 2, 2]).list()
            [[[1, 1], [2, 2], [3, 3]]]
            sage: SemistandardTableaux([3,2,1], [2, 2, 2]).list()
            [[[1, 1, 2], [2, 3], [3]], [[1, 1, 3], [2, 2], [3]]]
        """
        return symmetrica.kostka_tab(self.shape, self.weight)


class SemistandardTableaux_size_weight(SemistandardTableaux):
    r"""
    Semistandard tableaux of fixed size `n` and weight `\mu`.
    """

    def __init__(self, n, mu):
        r"""
        Initialize the class of semistandard tableaux of size ``n`` and
        weight ``mu``.

        .. WARNING::

            Input is not checked; please use :class:`SemistandardTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: SST = SemistandardTableaux(3, [2,1])
            sage: TestSuite(SST).run()                                                  # needs sage.modules
        """
        super().__init__(max_entry=len(mu),
                         category=FiniteEnumeratedSets())
        self.size = n
        self.weight = mu

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux(3, [2,1]))    # indirect doctest
            'Semistandard tableaux of size 3 and weight [2, 1]'
        """
        return "Semistandard tableaux of size %s and weight %s" % (self.size, self.weight)

    def __iter__(self):
        """
        EXAMPLES::

            sage: # needs sage.modules
            sage: [ t for t in SemistandardTableaux(3, [2,1]) ]
            [[[1, 1, 2]], [[1, 1], [2]]]
            sage: [ t for t in SemistandardTableaux(4, [2,2]) ]
            [[[1, 1, 2, 2]], [[1, 1, 2], [2]], [[1, 1], [2, 2]]]
            sage: sst = SemistandardTableaux(4, [2,2])
            sage: sst[0].parent() is sst
            True
        """
        from sage.combinat.partition import Partitions
        for p in Partitions(self.size):
            for sst in SemistandardTableaux_shape_weight(p, self.weight):
                yield self.element_class(self, sst)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: SemistandardTableaux(3, [2,1]).cardinality()                          # needs sage.modules
            2
            sage: SemistandardTableaux(4, [2,2]).cardinality()                          # needs sage.modules
            3
        """
        from sage.combinat.partition import Partitions
        c = 0
        for p in Partitions(self.size):
            c += SemistandardTableaux_shape_weight(p, self.weight).cardinality()
        return c

    def __contains__(self, x):
        """
        TESTS::

            sage: SST = SemistandardTableaux(6, [2,2,2])
            sage: all(sst in SST for sst in SST)                                        # needs sage.modules
            True
            sage: all(sst in SST for sst in SemistandardTableaux([3,2,1],[2,2,2]))      # needs sage.modules
            True
        """
        from sage.combinat.partition import Partition
        return x in SemistandardTableaux_shape_weight(Partition(
            [len(_) for _ in x]), self.weight)


#########################
# Row standard Tableaux #
#########################

class RowStandardTableaux(Tableaux):
    r"""
    A factory for the various classes of row standard tableaux.

    INPUT:

    - either a nonnegative integer (possibly specified with the keyword
      ``n``) or a partition

    OUTPUT: with no argument, the class of all standard tableaux

    - with a nonnegative integer argument, ``n``, the class of all standard
      tableaux of size ``n``

    - with a partition argument, the class of all standard tableaux of that
      shape

    A row standard tableau is a tableau that contains each of the
    entries from `1` to `n` exactly once and is increasing along rows.

    All classes of row standard tableaux are iterable.

    EXAMPLES::

        sage: ST = RowStandardTableaux(3); ST
        Row standard tableaux of size 3
        sage: ST.first()                                                                # needs sage.graphs
        [[1, 2, 3]]
        sage: ST.last()                                                                 # needs sage.graphs sage.modules
        [[3], [1], [2]]
        sage: ST.cardinality()                                                          # needs sage.graphs sage.modules
        10
        sage: ST.list()                                                                 # needs sage.graphs sage.modules
        [[[1, 2, 3]],
         [[2, 3], [1]],
         [[1, 2], [3]],
         [[1, 3], [2]],
         [[3], [2], [1]],
         [[2], [3], [1]],
         [[1], [3], [2]],
         [[1], [2], [3]],
         [[2], [1], [3]],
         [[3], [1], [2]]]

    .. SEEALSO::

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`RowStandardTableau`
        - :class:`StandardSkewTableaux`

    TESTS::

        sage: RowStandardTableaux()([])
        []
        sage: ST = RowStandardTableaux([2,2]); ST
        Row standard tableaux of shape [2, 2]
        sage: ST.first()                                                                # needs sage.graphs
        [[2, 4], [1, 3]]
        sage: ST.last()                                                                 # needs sage.graphs sage.modules
        [[2, 3], [1, 4]]
        sage: ST.cardinality()                                                          # needs sage.graphs sage.modules
        6
        sage: ST.list()                                                                 # needs sage.graphs sage.modules
        [[[2, 4], [1, 3]],
         [[3, 4], [1, 2]],
         [[1, 4], [2, 3]],
         [[1, 3], [2, 4]],
         [[1, 2], [3, 4]],
         [[2, 3], [1, 4]]]
        sage: RowStandardTableau([[3,4,5],[1,2]]).residue_sequence(3).standard_tableaux()
        Standard tableaux with 3-residue sequence (2,0,0,1,2) and multicharge (0)
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`RowStandardTableaux` for
        more information.

        TESTS::

            sage: RowStandardTableaux()
            Row standard tableaux
            sage: RowStandardTableaux(3)
            Row standard tableaux of size 3
            sage: RowStandardTableaux([2,1])
            Row standard tableaux of shape [2, 1]
            sage: RowStandardTableaux(0)
            Row standard tableaux of size 0

            sage: RowStandardTableaux(-1)
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a nonnegative integer or a partition
            sage: RowStandardTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a nonnegative integer or a partition
        """
        from sage.combinat.partition import _Partitions
        from sage.combinat.skew_partition import SkewPartitions

        if args:
            n = args[0]
        elif 'n' in kwargs:
            n = kwargs[n]
        else:
            n = None

        if n is None:
            return RowStandardTableaux_all()

        elif n in _Partitions:
            return RowStandardTableaux_shape(_Partitions(n))

        elif n in SkewPartitions():
            # from sage.combinat.skew_tableau import RowStandardSkewTableaux
            # return RowStandardSkewTableaux(n)
            raise NotImplementedError("row standard skew tableaux not yet implemented")

        if not isinstance(n, (int, Integer)) or n < 0:
            raise ValueError("the argument must be a nonnegative integer or a partition")

        return RowStandardTableaux_size(n)

    Element = RowStandardTableau

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[1,1],[2,3]] in RowStandardTableaux()
            False
            sage: [[1,2],[3,4]] in RowStandardTableaux()
            True
            sage: [[1,3],[2,4]] in RowStandardTableaux()
            True
            sage: [[1,3],[2,5]] in RowStandardTableaux()
            False
            sage: [] in RowStandardTableaux()
            True

        Check that integers are not contained in ``self``
        (see :issue:`14145`)::

            sage: 1 in RowStandardTableaux()
            False
        """
        if isinstance(x, RowStandardTableau):
            return True
        if Tableaux.__contains__(self, x):
            flatx = sorted(c for row in x for c in row)
            return (flatx == list(range(1, len(flatx)+1))
                    and all(row[i] < row[i+1] for row in x for i in range(len(row)-1)))
        return False


class RowStandardTableaux_all(RowStandardTableaux, DisjointUnionEnumeratedSets):
    """
    All row standard tableaux.
    """

    def __init__(self):
        r"""
        Initialize the class of all standard tableaux.

        .. WARNING::

            Input is not checked; please use :class:`RowStandardTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: ST = RowStandardTableaux()
            sage: TestSuite(ST).run()                                                   # needs sage.graphs
        """
        RowStandardTableaux.__init__(self)
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(NonNegativeIntegers(), RowStandardTableaux_size),
                                             facade=True, keepkey=False)

    def _repr_(self):
        """
        TESTS::

            sage: repr(RowStandardTableaux())    # indirect doctest
            'Row standard tableaux'
        """
        return "Row standard tableaux"


class RowStandardTableaux_size(RowStandardTableaux, DisjointUnionEnumeratedSets):
    """
    Row standard tableaux of fixed size `n`.

    EXAMPLES::

        sage: [t for t in RowStandardTableaux(1)]                                       # needs sage.graphs
        [[[1]]]
        sage: [t for t in RowStandardTableaux(2)]                                       # needs sage.graphs
        [[[1, 2]], [[2], [1]], [[1], [2]]]
        sage: list(RowStandardTableaux(3))                                              # needs sage.graphs
        [[[1, 2, 3]],
         [[2, 3], [1]],
         [[1, 2], [3]],
         [[1, 3], [2]],
         [[3], [2], [1]],
         [[2], [3], [1]],
         [[1], [3], [2]],
         [[1], [2], [3]],
         [[2], [1], [3]],
         [[3], [1], [2]]]

    TESTS::

        sage: TestSuite(RowStandardTableaux(4)).run()                                   # needs sage.graphs

        sage: RowStandardTableaux(3).cardinality()                                      # needs sage.libs.flint
        10
        sage: ns = [1,2,3,4,5,6]
        sage: sts = [RowStandardTableaux(n) for n in ns]
        sage: all(st.cardinality() == len(st.list()) for st in sts)                     # needs sage.graphs
        True
        sage: RowStandardTableaux(40).cardinality()  # not tested, too long
        2063837185739279909309355007659204891024472174278
    """

    def __init__(self, n):
        r"""
        Initialize the class of all row standard tableaux of size ``n``.

        .. WARNING::

            Input is not checked; please use :class:`RowStandardTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: TestSuite(RowStandardTableaux(0)).run()                               # needs sage.graphs
            sage: TestSuite(RowStandardTableaux(3)).run()                               # needs sage.graphs
        """
        RowStandardTableaux.__init__(self)
        from sage.combinat.partition import Partitions_n
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(Partitions_n(n), RowStandardTableaux_shape),
                                             facade=True, keepkey=False)
        self._size = Integer(n)

    def _repr_(self):
        """
        TESTS::

            sage: RowStandardTableaux(3)
            Row standard tableaux of size 3
        """
        return "Row standard tableaux of size %s" % self._size

    def __contains__(self, x):
        """
        TESTS::

            sage: ST3 = RowStandardTableaux(3)
            sage: all(st in ST3 for st in ST3)                                          # needs sage.graphs
            True
            sage: ST4 = RowStandardTableaux(4)
            sage: [x for x in ST4 if x in ST3]                                          # needs sage.graphs
            []

        Check that :issue:`14145` is fixed::

            sage: 1 in RowStandardTableaux(4)
            False
        """
        return RowStandardTableaux.__contains__(self, x) and sum(map(len, x)) == self._size

    def an_element(self):
        r"""
        Return a particular element of the class.

        EXAMPLES::

            sage: RowStandardTableaux(4).an_element()
            [[1, 2, 3, 4]]
        """
        if self._size == 0:
            return self.element_class(self, [])
        elif self._size == 1:
            return self.element_class(self, [[1]])
        else:
            return self.element_class(self, [range(1, self._size + 1)])


class RowStandardTableaux_shape(RowStandardTableaux):
    """
    Row Standard tableaux of a fixed shape `p`.
    """

    def __init__(self, p):
        r"""
        Initialize the class of all row standard tableaux of a given shape.

        .. WARNING::

            Input is not checked; please use :class:`RowStandardTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: TestSuite( RowStandardTableaux([2,1,1]) ).run()                       # needs sage.graphs
        """
        super().__init__(category=FiniteEnumeratedSets())
        self.shape = p

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: ST = RowStandardTableaux([2,1,1])
            sage: all(st in ST for st in ST)                                            # needs sage.graphs
            True
            sage: len([x for x in RowStandardTableaux(4) if x in ST])                   # needs sage.graphs
            12
            sage: ST.cardinality()
            12
        """
        return RowStandardTableaux.__contains__(self, x) and [len(_) for _ in x] == self.shape

    def _repr_(self):
        """
        TESTS::

            sage: RowStandardTableaux([2,1,1])
            Row standard tableaux of shape [2, 1, 1]
        """
        return "Row standard tableaux of shape {}".format(self.shape)

    def __iter__(self):
        r"""
        An iterator for the row standard Young tableaux associated to the
        shape `p` of ``self``.

        EXAMPLES::

            sage: [t for t in RowStandardTableaux([2,2])]                               # needs sage.graphs
            [[[2, 4], [1, 3]],
             [[3, 4], [1, 2]],
             [[1, 4], [2, 3]],
             [[1, 3], [2, 4]],
             [[1, 2], [3, 4]],
             [[2, 3], [1, 4]]]
            sage: [t for t in RowStandardTableaux([3,2])]                               # needs sage.graphs
            [[[2, 4, 5], [1, 3]],
             [[3, 4, 5], [1, 2]],
             [[1, 4, 5], [2, 3]],
             [[1, 3, 5], [2, 4]],
             [[1, 2, 5], [3, 4]],
             [[1, 2, 3], [4, 5]],
             [[1, 2, 4], [3, 5]],
             [[1, 3, 4], [2, 5]],
             [[2, 3, 4], [1, 5]],
             [[2, 3, 5], [1, 4]]]
            sage: st = RowStandardTableaux([2,1])
            sage: st[0].parent() is st                                                  # needs sage.graphs
            True
        """
        partial_sums = [sum(self.shape[:i]) for i in range(len(self.shape)+1)]

        # convert self.shape into a poset
        relations = []
        m = 1
        for row in self.shape:
            relations += [(m+i, m+i+1) for i in range(row-1)]
            m += row
        P = Poset((range(1, self.shape.size()+1), relations))
        L = P.linear_extensions()
        # now run through the linear extensions and return the corresponding tableau
        for lin in L:
            linear_tab = list(permutation.Permutation(lin).inverse())
            tab = [linear_tab[partial_sums[i]:partial_sums[i+1]]
                   for i in range(len(self.shape))]
            yield self.element_class(self, tab)

    def cardinality(self):
        r"""
        Return the number of row standard tableaux of this shape.

        This is just the index of the corresponding Young subgroup in
        the full symmetric group.

        EXAMPLES::

            sage: RowStandardTableaux([3,2,1]).cardinality()
            60
            sage: RowStandardTableaux([2,2]).cardinality()
            6
            sage: RowStandardTableaux([5]).cardinality()
            1
            sage: RowStandardTableaux([6,5,5,3]).cardinality()
            1955457504
            sage: RowStandardTableaux([]).cardinality()
            1
        """
        return Integer(multinomial(list(self.shape)))


########################
# Standard Tableaux    #
########################

class StandardTableaux(SemistandardTableaux):
    """
    A factory for the various classes of standard tableaux.

    INPUT:

    - Either a nonnegative integer (possibly specified with the keyword ``n``)
      or a partition.

    OUTPUT: with no argument, the class of all standard tableaux

    - With a nonnegative integer argument, ``n``, the class of all standard
      tableaux of size ``n``

    - With a partition argument, the class of all standard tableaux of that
      shape.

    A standard tableau is a semistandard tableaux which contains each of the
    entries from 1 to ``n`` exactly once.

    All classes of standard tableaux are iterable.

    EXAMPLES::

        sage: ST = StandardTableaux(3); ST
        Standard tableaux of size 3
        sage: ST.first()
        [[1, 2, 3]]
        sage: ST.last()
        [[1], [2], [3]]
        sage: ST.cardinality()
        4
        sage: ST.list()
        [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]

    .. SEEALSO::

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableau`
        - :class:`StandardSkewTableaux`

    TESTS::

        sage: StandardTableaux()([])
        []
        sage: ST = StandardTableaux([2,2]); ST
        Standard tableaux of shape [2, 2]
        sage: ST.first()
        [[1, 3], [2, 4]]
        sage: ST.last()
        [[1, 2], [3, 4]]
        sage: ST.cardinality()
        2
        sage: ST.list()
        [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
        sage: StandardTableau([[1,2,3],[4,5]]).residue_sequence(3).standard_tableaux()  # needs sage.groups
        Standard tableaux with 3-residue sequence (0,1,2,2,0) and multicharge (0)
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`StandardTableaux` for
        more information.

        TESTS::

            sage: StandardTableaux()
            Standard tableaux
            sage: StandardTableaux(3)
            Standard tableaux of size 3
            sage: StandardTableaux([2,1])
            Standard tableaux of shape [2, 1]
            sage: StandardTableaux(0)
            Standard tableaux of size 0
            sage: StandardTableaux(n=3)
            Standard tableaux of size 3

            sage: StandardTableaux(-1)
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a nonnegative integer or a partition
            sage: StandardTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a nonnegative integer or a partition
        """
        from sage.combinat.partition import _Partitions
        from sage.combinat.skew_partition import SkewPartitions

        if args:
            n = args[0]
        elif 'n' in kwargs:
            n = kwargs['n']
        else:
            n = None

        if n is None:
            return StandardTableaux_all()

        elif n in _Partitions:
            return StandardTableaux_shape(_Partitions(n))

        elif n in SkewPartitions():
            from sage.combinat.skew_tableau import StandardSkewTableaux
            return StandardSkewTableaux(n)

        if not isinstance(n, (int, Integer)) or n < 0:
            raise ValueError("the argument must be a nonnegative integer or a partition")

        return StandardTableaux_size(n)

    Element = StandardTableau

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[1,1],[2,3]] in StandardTableaux()
            False
            sage: [[1,2],[3,4]] in StandardTableaux()
            True
            sage: [[1,3],[2,4]] in StandardTableaux()
            True
            sage: [[1,3],[2,5]] in StandardTableaux()
            False
            sage: [] in StandardTableaux()
            True

        Check that :issue:`14145` is fixed::

            sage: 1 in StandardTableaux()
            False
        """
        if isinstance(x, StandardTableau):
            return True
        elif Tableaux.__contains__(self, x):
            flatx = sorted(c for row in x for c in row)
            return flatx == list(range(1, len(flatx)+1)) and (len(x) == 0 or
                     (all(row[i] < row[i+1] for row in x for i in range(len(row)-1)) and
                      all(x[r][c] < x[r+1][c] for r in range(len(x)-1)
                          for c in range(len(x[r+1])))
                      ))
        return False


class StandardTableaux_all(StandardTableaux, DisjointUnionEnumeratedSets):
    """
    All standard tableaux.
    """

    def __init__(self):
        r"""
        Initialize the class of all standard tableaux.

        TESTS::

            sage: ST = StandardTableaux()
            sage: TestSuite(ST).run()
        """
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(NonNegativeIntegers(), StandardTableaux_size),
                                             facade=True, keepkey=False)

    def _repr_(self):
        """
        TESTS::

            sage: repr(StandardTableaux())    # indirect doctest
            'Standard tableaux'
        """
        return "Standard tableaux"


class StandardTableaux_size(StandardTableaux, DisjointUnionEnumeratedSets):
    """
    Standard tableaux of fixed size `n`.

    EXAMPLES::

        sage: [ t for t in StandardTableaux(1) ]
        [[[1]]]
        sage: [ t for t in StandardTableaux(2) ]
        [[[1, 2]], [[1], [2]]]
        sage: [ t for t in StandardTableaux(3) ]
        [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]
        sage: StandardTableaux(4)[:]
        [[[1, 2, 3, 4]],
         [[1, 3, 4], [2]],
         [[1, 2, 4], [3]],
         [[1, 2, 3], [4]],
         [[1, 3], [2, 4]],
         [[1, 2], [3, 4]],
         [[1, 4], [2], [3]],
         [[1, 3], [2], [4]],
         [[1, 2], [3], [4]],
         [[1], [2], [3], [4]]]

    TESTS::

        sage: TestSuite( StandardTableaux(4) ).run()
    """

    def __init__(self, n):
        r"""
        Initialize the class of all standard tableaux of size `n`.

        .. WARNING::

            Input is not checked; please use :class:`StandardTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: TestSuite( StandardTableaux(0) ).run()
            sage: TestSuite( StandardTableaux(3) ).run()
        """
        StandardTableaux.__init__(self)
        from sage.combinat.partition import Partitions_n
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(Partitions_n(n), StandardTableaux_shape),
                                             category=FiniteEnumeratedSets(),
                                             facade=True, keepkey=False)
        self.size = Integer(n)

    def _repr_(self):
        """
        TESTS::

            sage: StandardTableaux(3)
            Standard tableaux of size 3
        """
        return "Standard tableaux of size %s" % self.size

    def __contains__(self, x):
        """
        TESTS::

            sage: ST3 = StandardTableaux(3)
            sage: all(st in ST3 for st in ST3)
            True
            sage: ST4 = StandardTableaux(4)
            sage: [x for x in ST4 if x in ST3]
            []

        Check that :issue:`14145` is fixed::

            sage: 1 in StandardTableaux(4)
            False
        """
        return StandardTableaux.__contains__(self, x) and sum(map(len, x)) == self.size

    def cardinality(self):
        r"""
        Return the number of all standard tableaux of size ``n``.

        The number of standard tableaux of size `n` is equal to the
        number of involutions in the symmetric group `S_n`.
        This is a consequence of the symmetry of the RSK
        correspondence, that if `\sigma \mapsto (P, Q)`, then
        `\sigma^{-1} \mapsto (Q, P)`. For more information, see
        :wikipedia:`Robinson-Schensted-Knuth_correspondence#Symmetry`.

        ALGORITHM:

        The algorithm uses the fact that standard tableaux of size
        ``n`` are in bijection with the involutions of size ``n``,
        (see page 41 in section 4.1 of [Ful1997]_).  For each number of
        fixed points, you count the number of ways to choose those
        fixed points multiplied by the number of perfect matchings on
        the remaining values.

        EXAMPLES::

            sage: StandardTableaux(3).cardinality()
            4
            sage: ns = [1,2,3,4,5,6]
            sage: sts = [StandardTableaux(n) for n in ns]
            sage: all(st.cardinality() == len(st.list()) for st in sts)
            True

        The cardinality can be computed without constructing all elements in
        this set, so this computation is fast (see also :issue:`28273`)::

            sage: StandardTableaux(500).cardinality()
            423107565308608549951551753690...221285999236657443927937253376

        TESTS::

            sage: def cardinality_using_hook_formula(n):
            ....:     c = 0
            ....:     for p in Partitions(n):
            ....:         c += StandardTableaux(p).cardinality()
            ....:     return c
            sage: all(cardinality_using_hook_formula(i) == StandardTableaux(i).cardinality() for i in range(10))
            True
        """
        tableaux_number = self.size % 2  # identity involution
        fixed_point_numbers = list(range(tableaux_number, self.size + 1 - tableaux_number, 2))

        # number of involutions of size "size" (number of ways to
        # choose "fixed_point_number" out of "size" elements *
        # number of involutions without fixed point of size
        # "size" - "fixed_point_number")
        for fixed_point_number in fixed_point_numbers:
            tableaux_number += (self.size.binomial(fixed_point_number) *
                                prod(range(1, self.size - fixed_point_number, 2)))

        return tableaux_number

    def random_element(self):
        r"""
        Return a random ``StandardTableau`` with uniform probability.

        This algorithm uses the fact that the Robinson-Schensted
        correspondence returns a pair of identical standard Young
        tableaux (SYTs) if and only if the permutation was an involution.
        Thus, generating a random SYT is equivalent to generating a
        random involution.

        To generate an involution, we first need to choose its number of
        fixed points `k` (if the size of the involution is even, the
        number of fixed points will be even, and if the size is odd, the
        number of fixed points will be odd). To do this, we choose a
        random integer `r` between 0 and the number `N` of all
        involutions of size `n`. We then decompose the interval
        `\{ 1, 2, \ldots, N \}` into subintervals whose lengths are the
        numbers of involutions of size `n` with respectively `0`, `1`,
        `\ldots`, `\left \lfloor N/2 \right \rfloor` fixed points. The
        interval in which our random integer `r` lies then decides how
        many fixed points our random involution will have. We then
        place those fixed points randomly and then compute a perfect
        matching (an involution without fixed points) on the remaining
        values.

        EXAMPLES::

            sage: StandardTableaux(10).random_element() # random
            [[1, 3, 6], [2, 5, 7], [4, 8], [9], [10]]
            sage: StandardTableaux(0).random_element()
            []
            sage: StandardTableaux(1).random_element()
            [[1]]

        TESTS::

            sage: all(StandardTableaux(10).random_element() in StandardTableaux(10) for i in range(20))
            True
        """
        from sage.misc.prandom import randrange
        from sage.misc.prandom import sample
        from sage.combinat.perfect_matching import PerfectMatchings
        from sage.combinat.permutation import from_cycles
        # We compute the number of involutions of size ``size``.
        involution_index = randrange(0, StandardTableaux(self.size).cardinality())
        # ``involution_index`` is our random integer `r`.
        partial_sum = 0
        fixed_point_number = self.size % 2
        # ``fixed_point_number`` will become `k`.
        while True:
            # We add the number of involutions with ``fixed_point_number``
            # fixed points.
            partial_sum += binomial(self.size, fixed_point_number) * \
                prod(range(1, self.size - fixed_point_number, 2))
            # If the partial sum is greater than the involution index,
            # then the random involution that we want to generate has
            # ``fixed_point_number`` fixed points.
            if partial_sum > involution_index:
                break
            fixed_point_number += 2
        # We generate a subset of size "fixed_point_number" of the set {1,
        # ..., size}.
        fixed_point_positions = set(sample(range(1, self.size + 1), fixed_point_number))
        # We generate a list of tuples which will form the cycle
        # decomposition of our random involution. This list contains
        # singletons (corresponding to the fixed points of the
        # involution) and pairs (forming a perfect matching on the
        # remaining values).
        matching = PerfectMatchings(set(range(1, self.size + 1))
                                    - set(fixed_point_positions)).random_element()
        permutation_cycle_rep = ([(fixed_point,) for fixed_point in fixed_point_positions]
                                 + [tuple(ab) for ab in matching])
        return from_cycles(self.size, permutation_cycle_rep).robinson_schensted()[0]


class StandardTableaux_shape(StandardTableaux):
    """
    Semistandard tableaux of a fixed shape `p`.
    """

    def __init__(self, p):
        r"""
        Initialize the class of all semistandard tableaux of a given shape.

        .. WARNING::

            Input is not checked; please use :class:`StandardTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: TestSuite( StandardTableaux([2,1,1]) ).run()
        """
        super().__init__(category=FiniteEnumeratedSets())
        self.shape = p

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: ST = StandardTableaux([2,1,1])
            sage: all(st in ST for st in ST)
            True
            sage: len([x for x in StandardTableaux(4) if x in ST])
            3
            sage: ST.cardinality()
            3

        Check that :issue:`14145` is fixed::

            sage: 1 in StandardTableaux([2,1,1])
            False
        """
        return StandardTableaux.__contains__(self, x) and [len(_) for _ in x] == self.shape

    def _repr_(self):
        """
        TESTS::

            sage: repr(StandardTableaux([2,1,1]))    # indirect doctest
            'Standard tableaux of shape [2, 1, 1]'
        """
        return "Standard tableaux of shape %s" % str(self.shape)

    def cardinality(self):
        r"""
        Return the number of standard Young tableaux of this shape.

        This method uses the so-called *hook length formula*, a formula
        for the number of Young tableaux associated with a given
        partition. The formula says the following: Let `\lambda` be a
        partition. For each cell `c` of the Young diagram of `\lambda`,
        let the *hook length* of `c` be defined as `1` plus the number of
        cells horizontally to the right of `c` plus the number of cells
        vertically below `c`. The number of standard Young tableaux of
        shape `\lambda` is then `n!` divided by the product of the hook
        lengths of the shape of `\lambda`, where `n = |\lambda|`.

        For example, consider the partition ``[3,2,1]`` of ``6`` with
        Ferrers diagram::

            # # #
            # #
            #

        When we fill in the cells with their respective hook lengths, we
        obtain::

            5 3 1
            3 1
            1

        The hook length formula returns

        .. MATH::

            \frac{6!}{5 \cdot 3 \cdot 1 \cdot 3 \cdot 1 \cdot 1} = 16.

        EXAMPLES::

            sage: StandardTableaux([3,2,1]).cardinality()
            16
            sage: StandardTableaux([2,2]).cardinality()
            2
            sage: StandardTableaux([5]).cardinality()
            1
            sage: StandardTableaux([6,5,5,3]).cardinality()
            6651216
            sage: StandardTableaux([]).cardinality()
            1

        REFERENCES:

        - http://mathworld.wolfram.com/HookLengthFormula.html
        """
        pi = self.shape
        number = factorial(sum(pi))
        for row in pi.hook_lengths():
            for col in row:
                number //= col
        return number

    def __iter__(self):
        r"""
        An iterator for the standard Young tableaux associated to the
        shape `p` of ``self``.

        EXAMPLES::

            sage: [t for t in StandardTableaux([2,2])]
            [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
            sage: [t for t in StandardTableaux([3,2])]
            [[[1, 3, 5], [2, 4]],
             [[1, 2, 5], [3, 4]],
             [[1, 3, 4], [2, 5]],
             [[1, 2, 4], [3, 5]],
             [[1, 2, 3], [4, 5]]]
            sage: st = StandardTableaux([2,1])
            sage: st[0].parent() is st
            True
        """

        pi = self.shape
        # Set the initial tableau by filling it in going down the columns
        tableau = [[None]*n for n in pi]
        size = sum(pi)
        row = 0
        col = 0
        for i in range(size):
            tableau[row][col] = i+1

            # If we can move down, then do it;
            # otherwise, move to the next column over
            if (row + 1 < len(pi) and col < pi[row+1]):
                row += 1
            else:
                row = 0
                col += 1

        yield self.element_class(self, tableau)

        # iterate until we reach the last tableau which is
        # filled with the row indices.
        last_tableau = sum([[row]*l for (row, l) in enumerate(pi)], [])

        # Convert the tableau to "vector format"
        # tableau_vector[i] is the row that number i
        # is in
        tableau_vector = [None]*size
        for row in range(len(pi)):
            for col in range(pi[row]):
                tableau_vector[tableau[row][col]-1] = row

        while tableau_vector != last_tableau:
            # Locate the smallest integer j such that j is not
            # in the lowest corner of the subtableau T_j formed by
            # 1,...,j.  This happens to be first j such that
            # ntableau_vector[j]<tableau_vector[j-1].
            # l will correspond to the shape of T_j
            l = [0]*size
            l[0] = 1
            j = 0
            for i in range(1, size):
                l[tableau_vector[i]] += 1
                if (tableau_vector[i] < tableau_vector[i-1]):
                    j = i
                    break

            # Find the last nonzero row of l and store it in k
            i = size - 1
            while l[i] == 0:
                i -= 1
            k = i

            # Find a new row for the letter j (next lowest corner)
            t = l[1 + tableau_vector[j]]
            i = k
            while l[i] != t:
                i -= 1

            # Move the letter j to row i
            tableau_vector[j] = i
            l[i] -= 1

            # Fill in the columns of T_j using 1,...,j-1 in increasing order
            m = 0
            while m < j:
                r = 0
                while l[r] != 0:
                    tableau_vector[m] = r
                    l[r] -= 1
                    m += 1
                    r += 1

            # Convert the tableau vector back to the regular tableau
            # format
            row_count = [0] * len(pi)
            tableau = [[None] * n for n in pi]

            for i in range(size):
                tableau[tableau_vector[i]][row_count[tableau_vector[i]]] = i+1
                row_count[tableau_vector[i]] += 1

            yield self.element_class(self, tableau)

        return

    def list(self):
        r"""
        Return a list of the standard Young tableaux of the specified shape.

        EXAMPLES::

            sage: StandardTableaux([2,2]).list()
            [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
            sage: StandardTableaux([5]).list()
            [[[1, 2, 3, 4, 5]]]
            sage: StandardTableaux([3,2,1]).list()
            [[[1, 4, 6], [2, 5], [3]],
             [[1, 3, 6], [2, 5], [4]],
             [[1, 2, 6], [3, 5], [4]],
             [[1, 3, 6], [2, 4], [5]],
             [[1, 2, 6], [3, 4], [5]],
             [[1, 4, 5], [2, 6], [3]],
             [[1, 3, 5], [2, 6], [4]],
             [[1, 2, 5], [3, 6], [4]],
             [[1, 3, 4], [2, 6], [5]],
             [[1, 2, 4], [3, 6], [5]],
             [[1, 2, 3], [4, 6], [5]],
             [[1, 3, 5], [2, 4], [6]],
             [[1, 2, 5], [3, 4], [6]],
             [[1, 3, 4], [2, 5], [6]],
             [[1, 2, 4], [3, 5], [6]],
             [[1, 2, 3], [4, 5], [6]]]
        """
        return list(self)

    def random_element(self):
        """
        Return a random standard tableau of the given shape using the
        Greene-Nijenhuis-Wilf Algorithm.

        EXAMPLES::

            sage: t = StandardTableaux([2,2]).random_element()
            sage: t.shape()
            [2, 2]
            sage: StandardTableaux([]).random_element()
            []
        """
        p = self.shape
        t = [[None] * n for n in p]

        # Get the cells in the Young diagram
        cells = []
        for i in range(len(p)):
            for j in range(p[i]):
                cells.append((i, j))

        m = sum(p)
        while m > 0:
            # Choose a cell at random
            cell = random.choice(cells)

            # Find a corner
            inner_corners = p.corners()
            while cell not in inner_corners:
                hooks = []
                for k in range(cell[1] + 1, p[cell[0]]):
                    hooks.append((cell[0], k))
                for k in range(cell[0] + 1, len(p)):
                    if p[k] > cell[1]:
                        hooks.append((k, cell[1]))
                cell = random.choice(hooks)

            # Assign m to cell
            t[cell[0]][cell[1]] = m
            p = p.remove_cell(cell[0])
            cells.remove(cell)
            m -= 1

        return self.element_class(self, t)


##########################
# Symmetric group action #
##########################
def unmatched_places(w, open, close):
    """
    Given a word ``w`` and two letters ``open`` and
    ``close`` to be treated as opening and closing
    parentheses (respectively), return a pair ``(xs, ys)``
    that encodes the positions of the unmatched
    parentheses after the standard parenthesis matching
    procedure is applied to ``w``.

    More precisely, ``xs`` will be the list of all ``i``
    such that ``w[i]`` is an unmatched closing parenthesis,
    while ``ys`` will be the list of all ``i`` such that
    ``w[i]`` is an unmatched opening parenthesis. Both
    lists returned are in increasing order.

    EXAMPLES::

        sage: from sage.combinat.tableau import unmatched_places
        sage: unmatched_places([2,2,2,1,1,1],2,1)
        ([], [])
        sage: unmatched_places([1,1,1,2,2,2],2,1)
        ([0, 1, 2], [3, 4, 5])
        sage: unmatched_places([], 2, 1)
        ([], [])
        sage: unmatched_places([1,2,4,6,2,1,5,3],2,1)
        ([0], [1])
        sage: unmatched_places([2,2,1,2,4,6,2,1,5,3], 2, 1)
        ([], [0, 3])
        sage: unmatched_places([3,1,1,1,2,1,2], 2, 1)
        ([1, 2, 3], [6])
    """
    lw = len(w)
    places_open = []
    places_close = []
    for i in range(lw):
        letter = w[i]
        if letter == open:
            places_open.append(i)
        elif letter == close:
            if places_open == []:
                places_close.append(i)
            else:
                places_open.pop()
    return places_close, places_open


def symmetric_group_action_on_values(word, perm):
    """
    Return the image of the word ``word`` under the
    Lascoux-Schuetzenberger action of the permutation
    ``perm``.

    See :meth:`Tableau.symmetric_group_action_on_values`
    for the definition of the Lascoux-Schuetzenberger
    action on semistandard tableaux. The transformation that
    the reading word of the tableau undergoes in said
    definition is precisely the Lascoux-Schuetzenberger
    action on words.

    EXAMPLES::

        sage: from sage.combinat.tableau import symmetric_group_action_on_values
        sage: symmetric_group_action_on_values([1,1,1],[1,3,2])
        [1, 1, 1]
        sage: symmetric_group_action_on_values([1,1,1],[2,1,3])
        [2, 2, 2]
        sage: symmetric_group_action_on_values([1,2,1],[2,1,3])
        [2, 2, 1]
        sage: symmetric_group_action_on_values([2,2,2],[2,1,3])
        [1, 1, 1]
        sage: symmetric_group_action_on_values([2,1,2],[2,1,3])
        [2, 1, 1]
        sage: symmetric_group_action_on_values([2,2,3,1,1,2,2,3],[1,3,2])
        [2, 3, 3, 1, 1, 2, 3, 3]
        sage: symmetric_group_action_on_values([2,1,1],[2,1])
        [2, 1, 2]
        sage: symmetric_group_action_on_values([2,2,1],[2,1])
        [1, 2, 1]
        sage: symmetric_group_action_on_values([1,2,1],[2,1])
        [2, 2, 1]
    """
    w = list(word)
    ts = permutation.Permutations()(perm).reduced_word()
    for r in reversed(ts):
        l = r + 1
        places_r, places_l = unmatched_places(w, l, r)

        # Now change the number of l's and r's in the new word
        nbl = len(places_l)
        nbr = len(places_r)
        ma = max(nbl, nbr)
        dif = ma - min(nbl, nbr)
        if ma == nbl:
            for i in places_l[:dif]:
                w[i] = r
        else:
            for i in places_r[nbr-dif:]:
                w[i] = l
    return w


class Tableau_class(Tableau):
    """
    This exists solely for unpickling ``Tableau_class`` objects.
    """

    def __setstate__(self, state):
        r"""
        Unpickle old ``Tableau_class`` objects.

        TESTS::

            sage: loads(b'x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1+IL\xcaIM,\xe5\n\x81\xd0\xf1\xc99\x89\xc5\xc5\\\x85\x8c\x9a\x8d\x85L\xb5\x85\xcc\x1a\xa1\xac\xf1\x19\x89\xc5\x19\x85,~@VNfqI!kl![l!;\xc4\x9c\xa2\xcc\xbc\xf4b\xbd\xcc\xbc\x92\xd4\xf4\xd4"\xae\xdc\xc4\xec\xd4x\x18\xa7\x90#\x94\xd1\xb05\xa8\x9031\xb14I\x0f\x00\xf6\xae)7')
            [[1]]
            sage: loads(dumps( Tableau([[1]]) ))
            [[1]]
        """
        self.__class__ = Tableau
        self.__init__(Tableaux(), state['_list'])


##########################
# Increasing tableaux #
##########################

class IncreasingTableaux(Tableaux):
    """
    A factory class for the various classes of increasing tableaux.

    An *increasing tableau* is a tableau whose entries are positive
    integers that are strictly increasing across rows and strictly
    increasing down columns. Note that Sage uses the English convention
    for partitions and tableaux; the longer rows are displayed on top.

    INPUT:

    Keyword arguments:

    - ``size`` -- the size of the tableaux
    - ``shape`` -- the shape of the tableaux
    - ``eval`` -- the weight (also called binary content) of
      the tableaux, where values can be either `0` or `1` with position
      `i` being `1` if and only if `i` can appear in the tableaux
    - ``max_entry`` -- positive integer or infinity (``oo``); the maximum
      entry for the tableaux; if ``size`` or ``shape`` are specified,
      ``max_entry`` defaults to be ``size`` or the size of ``shape``

    Positional arguments:

    - the first argument is interpreted as either ``size`` or ``shape``
      according to whether it is an integer or a partition
    - the second keyword argument will always be interpreted as ``eval``

    .. WARNING::

        The ``eval`` is not the usual notion of ``eval`` or ``weight``,
        where the `i`-th entry counts how many `i`'s appear in the tableau.

    EXAMPLES::

        sage: IT = IncreasingTableaux([2,1]); IT
        Increasing tableaux of shape [2, 1] and maximum entry 3
        sage: IT.list()
        [[[1, 3], [2]], [[1, 2], [3]], [[1, 2], [2]], [[1, 3], [3]], [[2, 3], [3]]]

        sage: IT = IncreasingTableaux(3); IT
        Increasing tableaux of size 3 and maximum entry 3
        sage: IT.list()
        [[[1, 2, 3]],
         [[1, 3], [2]],
         [[1, 2], [3]],
         [[1, 2], [2]],
         [[1, 3], [3]],
         [[2, 3], [3]],
         [[1], [2], [3]]]

        sage: IT = IncreasingTableaux(3, max_entry=2); IT
        Increasing tableaux of size 3 and maximum entry 2
        sage: IT.list()
        [[[1, 2], [2]]]

        sage: IT = IncreasingTableaux(3, max_entry=4); IT
        Increasing tableaux of size 3 and maximum entry 4
        sage: IT.list()
        [[[1, 2, 3]],
         [[1, 2, 4]],
         [[1, 3, 4]],
         [[2, 3, 4]],
         [[1, 3], [2]],
         [[1, 2], [3]],
         [[1, 4], [2]],
         [[1, 2], [4]],
         [[1, 2], [2]],
         [[1, 4], [3]],
         [[1, 3], [4]],
         [[1, 3], [3]],
         [[1, 4], [4]],
         [[2, 4], [3]],
         [[2, 3], [4]],
         [[2, 3], [3]],
         [[2, 4], [4]],
         [[3, 4], [4]],
         [[1], [2], [3]],
         [[1], [2], [4]],
         [[1], [3], [4]],
         [[2], [3], [4]]]

        sage: IT = IncreasingTableaux(3, max_entry=oo); IT
        Increasing tableaux of size 3
        sage: IT[123]
        [[5, 7], [6]]

        sage: IT = IncreasingTableaux(max_entry=2)
        sage: list(IT)
        [[], [[1]], [[2]], [[1, 2]], [[1], [2]]]
        sage: IT[4]
        [[1], [2]]

        sage: IncreasingTableaux()[0]
        []

    .. SEEALSO::

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`
        - :class:`IncreasingTableau`
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`IncreasingTableaux`
        for more information.

        TESTS::

            sage: IncreasingTableaux()
            Increasing tableaux
            sage: IncreasingTableaux(3)
            Increasing tableaux of size 3 and maximum entry 3
            sage: IncreasingTableaux(size=3)
            Increasing tableaux of size 3 and maximum entry 3
            sage: IncreasingTableaux(0)
            Increasing tableaux of size 0 and maximum entry 0
            sage: IncreasingTableaux([2,1])
            Increasing tableaux of shape [2, 1] and maximum entry 3
            sage: IncreasingTableaux(shape=[2,1])
            Increasing tableaux of shape [2, 1] and maximum entry 3
            sage: IncreasingTableaux([])
            Increasing tableaux of shape [] and maximum entry 0
            sage: IncreasingTableaux(eval=(1,0,1))
            Increasing tableaux of size 2 and weight (1, 0, 1)
            sage: IncreasingTableaux(max_entry=3)
            Increasing tableaux with maximum entry 3
            sage: IncreasingTableaux(3, (1,0,1))
            Increasing tableaux of size 3 and weight (1, 0, 1)
            sage: IncreasingTableaux(3, shape=[2,1])
            Increasing tableaux of shape [2, 1] and maximum entry 3
            sage: IncreasingTableaux(3, (1,0,1), shape=[2,1])
            Increasing tableaux of shape [2, 1] and weight (1, 0, 1)
            sage: IncreasingTableaux(3, max_entry=4)
            Increasing tableaux of size 3 and maximum entry 4
            sage: IncreasingTableaux(3, max_entry=oo)
            Increasing tableaux of size 3
            sage: IncreasingTableaux([2, 1], max_entry=oo)
            Increasing tableaux of shape [2, 1]
            sage: IncreasingTableaux([2, 1], (1,0,1))
            Increasing tableaux of shape [2, 1] and weight (1, 0, 1)
            sage: mu = Partition([2,1]); IncreasingTableaux(mu, (1,0,1))
            Increasing tableaux of shape [2, 1] and weight (1, 0, 1)
            sage: IncreasingTableaux(3, (1,0,1), max_entry=3)
            Increasing tableaux of size 3 and weight (1, 0, 1)

            sage: IncreasingTableaux(3, shape=[2])
            Traceback (most recent call last):
            ...
            ValueError: size and shape are different sizes

            sage: IncreasingTableaux(3, (1,0,1,1,1))
            Traceback (most recent call last):
            ...
            ValueError: size is smaller than the number of labels

            sage: IncreasingTableaux([2],(1,0,1,1))
            Traceback (most recent call last):
            ...
            ValueError: number of boxes is smaller than the number of labels

            sage: IncreasingTableaux(2,(1,0,1), max_entry=4)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must match the weight

            sage: IncreasingTableaux(eval=(1,0,1), max_entry=oo)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must match the weight

            sage: IncreasingTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: shape must be a (skew) partition
        """
        from sage.combinat.partition import Partition, _Partitions
        # Process the keyword arguments -- allow for original syntax where
        #   n == size,  p== shape and mu == eval
        n = kwargs.get('n', None)
        size = kwargs.get('size', n)

        p = kwargs.get('p', None)
        shape = kwargs.get('shape', p)

        wt = kwargs.get("eval", None)
        wt = kwargs.get("wt", wt)
        wt = kwargs.get("weight", wt)

        max_entry = kwargs.get('max_entry', None)

        # Process the positional arguments
        if args:
            # The first arg could be either a size or a shape
            if isinstance(args[0], (int, Integer)):
                if size is not None:
                    raise ValueError("size was specified more than once")
                else:
                    size = args[0]
            else:
                if shape is not None:
                    raise ValueError("the shape was specified more than once")
                shape = args[0]  # we check it's a partition later

        if len(args) == 2:
            # The second non-keyword argument is the weight
            if wt is not None:
                raise ValueError("the weight was specified more than once")
            else:
                wt = args[1]

        # Consistency checks
        if size is not None:
            if size not in NonNegativeIntegers():
                raise ValueError("size must be a nonnegative integer")
            size = Integer(size)

        if shape is not None:
            from sage.combinat.skew_partition import SkewPartitions
            # use in (and not isinstance) below so that lists can be used as
            # shorthand
            if shape in _Partitions:
                shape = Partition(shape)
            elif shape in SkewPartitions():
                raise NotImplementedError("skew increasing tableaux are not"
                                          " currently implemented")
                # from sage.combinat.skew_tableau import IncreasingSkewTableaux
                # return IncreasingSkewTableaux(shape, wt)
            else:
                raise ValueError("shape must be a (skew) partition")

        if wt is not None:
            wt = list(wt)
            k = len(wt) - 1
            while k >= 0 and wt[k] == 0:
                k -= 1
            wt = tuple(wt[:k+1])
            if not all(k in [0, 1] for k in wt):
                raise ValueError("wt must be a binary vector")
            if max_entry is not None and max_entry != len(wt):
                raise ValueError("the maximum entry must match the weight")

        is_inf = bool(max_entry == PlusInfinity())

        if max_entry is not None:
            if not is_inf and not isinstance(max_entry, (int, Integer)):
                raise ValueError("max_entry must be an integer or PlusInfinity")
            elif max_entry < 0:
                raise ValueError("max_entry must be nonnegative")

        if size is not None and shape is not None:
            if sum(shape) != size:
                # This could return an empty class instead of an error
                raise ValueError("size and shape are different sizes")

        if size is not None and wt is not None:
            if sum(wt) > size:
                # This could return an empty class instead of an error
                raise ValueError("size is smaller than the number of labels")

        # Dispatch appropriately
        if shape is not None and wt is not None:
            if sum(shape) < sum(wt):
                # This could return an empty class instead of an error
                raise ValueError("number of boxes is smaller than the number of labels")
            else:
                return IncreasingTableaux_shape_weight(shape, wt)

        if shape is not None:
            if is_inf:
                return IncreasingTableaux_shape_inf(shape)
            return IncreasingTableaux_shape(shape, max_entry)

        if wt is not None and size is not None:
            return IncreasingTableaux_size_weight(size, wt)

        if wt is not None:
            return IncreasingTableaux_size_weight(sum(wt), wt)

        if size is not None:
            if is_inf:
                return IncreasingTableaux_size_inf(size)
            return IncreasingTableaux_size(size, max_entry)

        return IncreasingTableaux_all(max_entry)

    Element = IncreasingTableau

    def __init__(self, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = IncreasingTableaux()
            sage: TestSuite(S).run()  # long time
        """
        if 'max_entry' in kwds:
            self.max_entry = kwds['max_entry']
            kwds.pop('max_entry')
        else:
            self.max_entry = None
        Tableaux.__init__(self, **kwds)

    def __getitem__(self, r):
        r"""
        The default implementation of ``__getitem__`` for enumerated sets
        does not allow slices so we override it.

        EXAMPLES::

            sage: IncreasingTableaux([4,3,3,2])[10:20]     # long time
            [[[1, 5, 8, 10], [2, 6, 9], [3, 7, 12], [4, 11]],
             [[1, 5, 8, 10], [2, 6, 9], [3, 7, 11], [4, 12]],
             [[1, 5, 8, 9], [2, 6, 11], [3, 7, 12], [4, 10]],
             [[1, 5, 8, 9], [2, 6, 10], [3, 7, 12], [4, 11]],
             [[1, 5, 8, 9], [2, 6, 10], [3, 7, 11], [4, 12]],
             [[1, 5, 7, 12], [2, 6, 10], [3, 8, 11], [4, 9]],
             [[1, 5, 7, 11], [2, 6, 10], [3, 8, 12], [4, 9]],
             [[1, 5, 7, 10], [2, 6, 11], [3, 8, 12], [4, 9]],
             [[1, 5, 7, 12], [2, 6, 9], [3, 8, 11], [4, 10]],
             [[1, 5, 7, 11], [2, 6, 9], [3, 8, 12], [4, 10]]]

            sage: IncreasingTableaux(size=2, max_entry=oo)[5]
            [[2], [3]]

            sage: IncreasingTableaux([2,1], max_entry=oo)[3]
            [[1, 2], [4]]

            sage: IncreasingTableaux(3, max_entry=4)[0:5]    # indirect doctest
            [[[1, 2, 3]], [[1, 2, 4]], [[1, 3, 4]], [[2, 3, 4]], [[1, 3], [2]]]

            sage: IncreasingTableaux([2,2], (1,0,1,1,1))[0]    # indirect doctest
            [[1, 4], [3, 5]]

            sage: IncreasingTableaux([1,1,1], max_entry=4)[0:4]
            [[[1], [2], [3]], [[1], [2], [4]], [[1], [3], [4]], [[2], [3], [4]]]

            sage: IncreasingTableaux(3, (1,0,1,1))[1]    # indirect doctest
            [[1, 4], [3]]

            sage: IncreasingTableaux(3)[:]  # indirect doctest
            [[[1, 2, 3]],
             [[1, 3], [2]],
             [[1, 2], [3]],
             [[1, 2], [2]],
             [[1, 3], [3]],
             [[2, 3], [3]],
             [[1], [2], [3]]]

            sage: IncreasingTableaux([2,2])[1]   # indirect doctest
            [[1, 2], [3, 4]]

        TESTS::

            sage: IncreasingTableaux()[5]
            [[1, 3], [2]]

            sage: IncreasingTableaux(max_entry=4)[5]
            [[1, 2]]

            sage: IncreasingTableaux()[:]
            Traceback (most recent call last):
            ...
            ValueError: infinite set

            sage: IncreasingTableaux(size=2, max_entry=oo)[:]
            Traceback (most recent call last):
            ...
            ValueError: infinite set
        """
        if isinstance(r, (int, Integer)):
            return self.unrank(r)
        elif isinstance(r, slice):
            start = 0 if r.start is None else r.start
            stop = r.stop
            if stop is None and not self.is_finite():
                raise ValueError('infinite set')
        else:
            raise ValueError('r must be an integer or a slice')
        count = 0
        tabs = []
        for t in self:
            if count == stop:
                break
            if count >= start:
                tabs.append(t)
            count += 1

        # this is to cope with empty slices endpoints like [:6] or [:}
        if count == stop or stop is None:
            return tabs
        raise IndexError('value out of range')

    def __contains__(self, t):
        """
        Return ``True`` if ``t`` can be interpreted as an
        :class:`IncreasingTableau`.

        TESTS::

            sage: T = sage.combinat.tableau.IncreasingTableaux_all()
            sage: [[1,2],[2]] in T
            True
            sage: [] in T
            True
            sage: Tableau([[1]]) in T
            True
            sage: StandardTableau([[1]]) in T
            True

            sage: [[1,2],[1]] in T
            False
            sage: [[1,1],[5]] in T
            False
            sage: [[1,7],[5]] in T
            True
            sage: [[1,3,2]] in T
            False

            sage: None in T
            False
        """
        if not t:
            return isinstance(t, (IncreasingTableau, list))

        if isinstance(t, IncreasingTableau):
            return (self.max_entry is None
                    or max(max(row) for row in t) <= self.max_entry)

        if not Tableaux.__contains__(self, t):
            return False

        for row in t:
            if not all(c > 0 for c in row):
                return False
            if not all(row[i] < row[i+1] for i in range(len(row)-1)):
                return False
        for row, next in zip(t, t[1:]):
            if not all(row[c] < next[c] for c in range(len(next))):
                return False
        return self.max_entry is None or max(max(row) for row in t) <= self.max_entry


class IncreasingTableaux_all(IncreasingTableaux, DisjointUnionEnumeratedSets):
    """
    All increasing tableaux.

    EXAMPLES::

        sage: T = IncreasingTableaux()
        sage: T.cardinality()
        +Infinity

        sage: T = IncreasingTableaux(max_entry=3)
        sage: list(T)
        [[],
         [[1]],
         [[2]],
         [[3]],
         [[1, 2]],
         [[1, 3]],
         [[2, 3]],
         [[1], [2]],
         [[1], [3]],
         [[2], [3]],
         [[1, 2, 3]],
         [[1, 3], [2]],
         [[1, 2], [3]],
         [[1, 2], [2]],
         [[1, 3], [3]],
         [[2, 3], [3]],
         [[1], [2], [3]]]

    TESTS::

        sage: T = IncreasingTableaux(max_entry=0)
        sage: list(T)
        [[]]
    """

    def __init__(self, max_entry=None):
        r"""
        Initialize the class of all increasing tableaux.

        .. WARNING::

            Input is not checked; please use :class:`IncreasingTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: from sage.combinat.tableau import IncreasingTableaux_all
            sage: T = IncreasingTableaux_all()
            sage: TestSuite(T).run()  # long time

            sage: T = IncreasingTableaux_all(max_entry=3)
            sage: TestSuite(T).run()

            sage: T = IncreasingTableaux_all(max_entry=0)
            sage: TestSuite(T).run()
        """

        def SST_n(n):
            return IncreasingTableaux_size(n, max_entry)
        if max_entry is None or max_entry == PlusInfinity():
            self.max_entry = None
            DisjointUnionEnumeratedSets.__init__(self,
                Family(NonNegativeIntegers(), SST_n),
                facade=True, keepkey=False)
        else:
            self.max_entry = max_entry
            DisjointUnionEnumeratedSets.__init__(self,
                Family(list(range(max_entry + 1)), SST_n),
                facade=True, keepkey=False)

    def _repr_(self):
        """
        TESTS::

            sage: IncreasingTableaux()
            Increasing tableaux

            sage: IncreasingTableaux(max_entry=3)
            Increasing tableaux with maximum entry 3
        """
        if self.max_entry is not None:
            return "Increasing tableaux with maximum entry %s" % str(self.max_entry)
        return "Increasing tableaux"


class IncreasingTableaux_size_inf(IncreasingTableaux):
    """
    Increasing tableaux of fixed size `n` with no maximum entry.
    """

    def __init__(self, n):
        r"""
        Initialize the class of increasing tableaux of size `n` with no
        maximum entry.

        .. WARNING::

            Input is not checked; please use :class:`IncreasingTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: from sage.combinat.tableau import IncreasingTableaux_size_inf
            sage: T = IncreasingTableaux_size_inf(3)
            sage: TestSuite(T).run()
        """
        super().__init__(category=InfiniteEnumeratedSets())
        self.size = n

    def _repr_(self):
        """
        TESTS::

            sage: IncreasingTableaux(3, max_entry=oo)
            Increasing tableaux of size 3
        """
        return "Increasing tableaux of size %s" % str(self.size)

    def __contains__(self, t):
        """
        Return ``True`` if ``t`` can be interpreted as an element of ``self``.

        TESTS::

            sage: T = IncreasingTableaux(3, max_entry=oo)
            sage: [[1,2],[5]] in T
            True
            sage: StandardTableau([[1, 2], [3]]) in T
            True

            sage: [] in T
            False
            sage: Tableau([[1]]) in T
            False
        """
        return IncreasingTableaux.__contains__(self, t) and sum(map(len, t)) == self.size

    def __iter__(self):
        """
        EXAMPLES::

            sage: IT = IncreasingTableaux(3, max_entry=oo)
            sage: [IT[t] for t in range(5)]
            [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]], [[1, 2, 4]]]
            sage: IT[1000]
            [[3, 13], [10]]
            sage: IT[0].parent() is IT
            True
        """
        from sage.combinat.partition import Partitions
        # Iterates through with maximum entry as order
        i = 1
        while True:
            for part in Partitions(self.size):
                if i != 1:
                    for k in range(1, self.size+1):
                        for c in integer_vectors_nk_fast_iter(self.size - k, i-1):
                            c.append(k)
                            for sst in IncreasingTableaux_shape_weight(part, tuple(c)):
                                yield self.element_class(self, sst)
                else:
                    for sst in IncreasingTableaux_shape_weight(part, (self.size,)):
                        yield self.element_class(self, sst)
            i += 1


class IncreasingTableaux_shape_inf(IncreasingTableaux):
    """
    Increasing tableaux of fixed shape `p` and no maximum entry.
    """

    def __init__(self, p):
        r"""
        Initialize the class of increasing tableaux of shape `p` and no
        maximum entry.

        .. WARNING::

            Input is not checked; please use :class:`IncreasingTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: IT = IncreasingTableaux([2,1], max_entry=oo)
            sage: type(IT)
            <class 'sage.combinat.tableau.IncreasingTableaux_shape_inf_with_category'>
            sage: TestSuite(IT).run()
        """
        super().__init__(category=InfiniteEnumeratedSets())
        self.shape = p

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: IT = IncreasingTableaux([2,1], max_entry=oo)
            sage: [[13, 67], [1467]] in IT
            True
            sage: IT = IncreasingTableaux([3,1], max_entry=oo)
            sage: [[13, 67], [1467]] in IT
            False
        """
        return IncreasingTableaux.__contains__(self, x) and [len(row) for row in x] == self.shape

    def _repr_(self):
        """
        TESTS::

            sage: IncreasingTableaux([2,1], max_entry=oo)
            Increasing tableaux of shape [2, 1]
        """
        return "Increasing tableaux of shape %s" % str(self.shape)

    def __iter__(self):
        """
        An iterator for the increasing partitions of shape ``p`` and no
        maximum entry. Iterates through with maximum entry as order.

        EXAMPLES::

            sage: IT = IncreasingTableaux([3, 1], max_entry=oo)
            sage: IT[1000]
            [[1, 2, 12], [6]]
            sage: [IT[t] for t in range(5)]
            [[[1, 3, 4], [2]],
             [[1, 2, 4], [3]],
             [[1, 2, 3], [4]],
             [[1, 3, 5], [2]],
             [[1, 2, 5], [3]]]
            sage: IT[0].parent() is IT
            True
        """
        # Iterates through with maximum entry as order
        i = 1
        n = sum(self.shape)
        while True:
            if i != 1:
                for k in range(1, n+1):
                    for c in integer_vectors_nk_fast_iter(n - k, i-1):
                        c.append(k)
                        for sst in IncreasingTableaux_shape_weight(self.shape, tuple(c)):
                            yield self.element_class(self, sst)
            else:
                for sst in IncreasingTableaux_shape_weight(self.shape, (n,)):
                    yield self.element_class(self, sst)
            i += 1


class IncreasingTableaux_size(IncreasingTableaux):
    """
    Increasing tableaux of fixed size `n`.
    """

    def __init__(self, n, max_entry=None):
        r"""
        Initialize the class of increasing tableaux of size `n`.

        .. WARNING::

            Input is not checked; please use :class:`IncreasingTableaux`
            to ensure the options are properly parsed.

        TESTS::

            sage: IT = IncreasingTableaux(3); IT
            Increasing tableaux of size 3 and maximum entry 3
            sage: type(IT)
            <class 'sage.combinat.tableau.IncreasingTableaux_size_with_category'>
            sage: TestSuite(IT).run()

            sage: IT = IncreasingTableaux(3, max_entry=6)
            sage: type(IT)
            <class 'sage.combinat.tableau.IncreasingTableaux_size_with_category'>
            sage: TestSuite(IT).run()
        """
        if max_entry is None:
            max_entry = n
        super().__init__(max_entry=max_entry,
                         category=FiniteEnumeratedSets())
        self.size = n

    def _repr_(self):
        """
        TESTS::

            sage: IncreasingTableaux(3)
            Increasing tableaux of size 3 and maximum entry 3

            sage: IncreasingTableaux(3, max_entry=6)
            Increasing tableaux of size 3 and maximum entry 6
        """
        return "Increasing tableaux of size %s and maximum entry %s" % (str(self.size), str(self.max_entry))

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[1,2],[2,3]] in IncreasingTableaux(3)
            False
            sage: [[1,2],[2,3]] in IncreasingTableaux(4)
            True
            sage: [[1,2],[2,3]] in IncreasingTableaux(4, max_entry=2)
            False
            sage: IT = IncreasingTableaux(4)
            sage: all(it in IT for it in IT)
            True
        """
        if self.size == 0:
            return x == []

        return (IncreasingTableaux.__contains__(self, x)
                and sum(map(len, x)) == self.size
                and max(max(row) for row in x) <= self.max_entry)

    def __iter__(self):
        """
        EXAMPLES::

            sage: [ t for t in IncreasingTableaux(2) ]
            [[[1, 2]], [[1], [2]]]
            sage: [ t for t in IncreasingTableaux(3) ]
            [[[1, 2, 3]],
             [[1, 3], [2]],
             [[1, 2], [3]],
             [[1, 2], [2]],
             [[1, 3], [3]],
             [[2, 3], [3]],
             [[1], [2], [3]]]

            sage: [ t for t in IncreasingTableaux(4, max_entry=3) ]
            [[[1, 2, 3], [2]],
             [[1, 2, 3], [3]],
             [[1, 2], [2, 3]],
             [[1, 2], [2], [3]],
             [[1, 3], [2], [3]]]

            sage: IT = IncreasingTableaux(3)
            sage: IT[0].parent() is IT
            True

            sage: list(IncreasingTableaux(size=4, max_entry=2))
            []
            sage: list(IncreasingTableaux(size=3, max_entry=2))
            [[[1, 2], [2]]]
            sage: list(IncreasingTableaux(0, max_entry=4))
            [[]]
            sage: list(IncreasingTableaux(3, max_entry=0))
            []
        """
        if self.size == 0:
            yield self.element_class(self, [])
            return

        from sage.combinat.partition import Partitions
        for part in Partitions(self.size):
            for sst in IncreasingTableaux_shape(part, self.max_entry):
                yield self.element_class(self, sst)


class IncreasingTableaux_shape(IncreasingTableaux):
    """
    Increasing tableaux of fixed shape `p` with a given max entry.

    An increasing tableau with max entry `i` is required to have all
    its entries less or equal to `i`. It is not required to actually
    contain an entry `i`.

    INPUT:

    - ``p`` -- a partition
    - ``max_entry`` -- the max entry; defaults to the size of ``p``
    """

    def __init__(self, p, max_entry=None):
        r"""
        Initialize the class of increasing tableaux of shape `p`, with a
        given ``max_entry``.

        .. WARNING::

            Input is not checked; please use :class:`IncreasingTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: IT = IncreasingTableaux([2,1])
            sage: TestSuite(IT).run()

            sage: IT = IncreasingTableaux([2,1], max_entry=5)
            sage: TestSuite(IT).run()
        """
        if max_entry is None:
            max_entry = sum(p)
        super().__init__(max_entry=max_entry,
                         category=FiniteEnumeratedSets())
        self.shape = p

    def __iter__(self):
        """
        An iterator for the increasing tableaux of the specified shape
        with the specified max entry.

        EXAMPLES::

            sage: [ t for t in IncreasingTableaux([3]) ]
            [[[1, 2, 3]]]
            sage: [ t for t in IncreasingTableaux([2,1]) ]
            [[[1, 3], [2]], [[1, 2], [3]], [[1, 2], [2]], [[1, 3], [3]], [[2, 3], [3]]]
            sage: [ t for t in IncreasingTableaux([3,1]) ]
            [[[1, 3, 4], [2]],
             [[1, 2, 4], [3]],
             [[1, 2, 3], [4]],
             [[1, 2, 3], [2]],
             [[1, 2, 3], [3]],
             [[1, 2, 4], [2]],
             [[1, 2, 4], [4]],
             [[1, 3, 4], [3]],
             [[1, 3, 4], [4]],
             [[2, 3, 4], [3]],
             [[2, 3, 4], [4]]]

            sage: [ t for t in IncreasingTableaux([3,1], max_entry=3) ]
            [[[1, 2, 3], [2]], [[1, 2, 3], [3]]]

            sage: IT = IncreasingTableaux([3])
            sage: IT[0].parent() is IT
            True

            sage: list(IncreasingTableaux(shape=[4,2], max_entry=2))
            []
            sage: list(IncreasingTableaux(shape=[2,1], max_entry=0))
            []
            sage: list(IncreasingTableaux(shape=[], max_entry=4))
            [[]]
        """
        n = sum(self.shape)
        if n == 0:
            yield self.element_class(self, [])
            return

        list_of_partial_binary_vecs = [[]]
        list_of_binary_vecs = []
        while list_of_partial_binary_vecs:
            active_vec = list_of_partial_binary_vecs.pop()
            if len(active_vec) < self.max_entry:
                list_of_partial_binary_vecs.append(active_vec + [0])
                list_of_partial_binary_vecs.append(active_vec + [1])
            else:
                list_of_binary_vecs.append(tuple(active_vec))
        for wt in list_of_binary_vecs:
            for sst in IncreasingTableaux_shape_weight(self.shape, wt):
                yield self.element_class(self, sst)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: IT = IncreasingTableaux([2,1])
            sage: all(it in IT for it in IT)
            True
            sage: len([x for x in IncreasingTableaux(3) if x in IT])
            5
            sage: IT.cardinality()
            5

            sage: IT = IncreasingTableaux([2,1], max_entry=4)
            sage: all(it in IT for it in IT)
            True
            sage: IT.cardinality()
            14
        """
        return (IncreasingTableaux.__contains__(self, x)
                and [len(row) for row in x] == self.shape)

    def _repr_(self):
        """
        TESTS::

            sage: IncreasingTableaux([2,1])
            Increasing tableaux of shape [2, 1] and maximum entry 3

            sage: IncreasingTableaux([2,1], max_entry=5)
            Increasing tableaux of shape [2, 1] and maximum entry 5
        """
        return "Increasing tableaux of shape %s and maximum entry %s" % (str(self.shape), str(self.max_entry))


class IncreasingTableaux_shape_weight(IncreasingTableaux_shape):
    r"""
    Increasing tableaux of fixed shape `p` and binary weight `wt`.
    """

    def __init__(self, p, wt):
        r"""
        Initialize the class of all increasing tableaux of shape ``p`` and
        weight ``mu``.

        .. WARNING::

            Input is not checked; please use :class:`IncreasingTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: IT = IncreasingTableaux([2,1], (1,0,1))
            sage: TestSuite(IT).run()
        """
        super().__init__(p, len(wt))
        self.weight = wt

    def _repr_(self):
        """
        TESTS::

            sage: IncreasingTableaux([2,1], (1,0,1))
            Increasing tableaux of shape [2, 1] and weight (1, 0, 1)
        """
        return "Increasing tableaux of shape %s and weight %s" % (self.shape, self.weight)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: IT = IncreasingTableaux([2,1], (1,0,1))
            sage: all(it in IT for it in IT)
            True
            sage: len([x for x in IncreasingTableaux(3) if x in IT])
            1
            sage: IT.cardinality()
            1
            sage: None in IT
            False

        Corner case check::

            sage: IT = IncreasingTableaux([], ())
            sage: [] in IT
            True
            sage: None in IT
            False
        """
        if not IncreasingTableaux_shape.__contains__(self, x):
            return False
        n = sum(self.shape)

        if n == 0 and not x:
            return True

        content_list = [0] * int(self.max_entry)
        for row in x:
            for i in row:
                content_list[i-1] = 1

        if tuple(content_list) != self.weight:
            return False

        return True

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: IncreasingTableaux([2,2], (1,0,1,1)).list()
            [[[1, 3], [3, 4]]]
            sage: IncreasingTableaux([2,2,2], (1,0,1,1,0,1,1)).list()
            [[[1, 3], [3, 6], [4, 7]],
             [[1, 3], [3, 4], [6, 7]],
             [[1, 4], [3, 6], [4, 7]],
             [[1, 4], [3, 6], [6, 7]],
             [[1, 3], [4, 6], [6, 7]]]

        TESTS::

            sage: IT = IncreasingTableaux([3,1], (1,0,1,1))
            sage: [IT[i] for i in range(2)]
            [[[1, 3, 4], [3]], [[1, 3, 4], [4]]]
            sage: IT[0].parent() is IT
            True

            sage: list(IncreasingTableaux(shape=[], eval=(0,0,0)))
            [[]]

        We explicitly check a corner case::

            sage: from sage.combinat.tableau import IncreasingTableaux_shape_weight
            sage: list(IncreasingTableaux_shape_weight(Partition([]), (0,1,1)))
            []
        """
        if not self.shape:
            if sum(self.weight) == 0:
                yield IncreasingTableau([])
            return

        tab = Tableau([[0] * k for k in self.shape])
        wt = self.weight
        list_of_partial_inc_tabs = [tab]
        list_of_inc_tabs = []
        while list_of_partial_inc_tabs:
            active_tab = list_of_partial_inc_tabs.pop()
            unfilled_spots = []
            for (r, c) in active_tab.cells():
                if active_tab[r][c] == 0:
                    unfilled_spots.append((r, c))
            if not unfilled_spots:
                top_value = max(active_tab.entries())
                if top_value == len(wt) - wt[::-1].index(1):
                    list_of_inc_tabs.append(self.element_class(self, active_tab))
                continue
            growth_spots = []
            for (r, c) in unfilled_spots:
                if (r-1, c) not in active_tab.cells() or active_tab[r-1][c] != 0:
                    if (r, c-1) not in active_tab.cells() or active_tab[r][c-1] != 0:
                        growth_spots.append((r, c))
            growth_choices = list(powerset(growth_spots))
            top_value = max(active_tab.entries())
            try:
                growth_num = wt[top_value:].index(1) + top_value + 1
            except ValueError:
                continue
            for growth_choice in growth_choices[1:]:
                new_tab = [[0] * k for k in self.shape]
                for (r, c) in active_tab.cells():
                    new_tab[r][c] = active_tab[r][c]
                for (r, c) in growth_choice:
                    new_tab[r][c] = growth_num
                list_of_partial_inc_tabs.append(Tableau(new_tab))
        yield from list_of_inc_tabs


class IncreasingTableaux_size_weight(IncreasingTableaux):
    r"""
    Increasing tableaux of fixed size `n` and weight `wt`.
    """

    def __init__(self, n, wt):
        r"""
        Initialize the class of increasing tableaux of size `n` and
        weight ``wt``.

        .. WARNING::

            Input is not checked; please use :class:`IncreasingTableaux` to
            ensure the options are properly parsed.

        TESTS::

            sage: IT = IncreasingTableaux(3, (1,0,1))
            sage: TestSuite(IT).run()
        """
        super().__init__(max_entry=len(wt),
                         category=FiniteEnumeratedSets())
        self.size = n
        self.weight = wt

    def _repr_(self):
        """
        TESTS::

            sage: IncreasingTableaux(3, (1,0,1))
            Increasing tableaux of size 3 and weight (1, 0, 1)
        """
        return "Increasing tableaux of size %s and weight %s" % (self.size, self.weight)

    def __iter__(self):
        """
        EXAMPLES::

            sage: [ T for T in IncreasingTableaux(3, (1,0,1)) ]
            [[[1, 3], [3]]]
            sage: [ T for T in IncreasingTableaux(4, (1,0,1,1)) ]
            [[[1, 3, 4], [3]],
             [[1, 3, 4], [4]],
             [[1, 3], [3, 4]],
             [[1, 3], [3], [4]],
             [[1, 4], [3], [4]]]
            sage: IT = IncreasingTableaux(4, (1,0,1,1))
            sage: IT[0].parent() is IT
            True
        """
        from sage.combinat.partition import Partitions
        for p in Partitions(self.size):
            for sst in IncreasingTableaux_shape_weight(p, self.weight):
                yield self.element_class(self, sst)

    def __contains__(self, x):
        """
        TESTS::

            sage: IT = IncreasingTableaux(4, (1,0,1,1))
            sage: all(it in IT for it in IT)
            True
            sage: all(it in IT for it in IncreasingTableaux([2,2], (1,0,1,1)))
            True
        """
        from sage.combinat.partition import _Partitions
        shape = [len(row) for row in x]
        if shape not in _Partitions:
            return False
        return x in IncreasingTableaux_shape_weight(_Partitions(shape), self.weight)


# October 2012: fixing outdated pickles which use classed being deprecated
register_unpickle_override('sage.combinat.tableau', 'Tableau_class', Tableau_class)
register_unpickle_override('sage.combinat.tableau', 'Tableaux_n', Tableaux_size)
register_unpickle_override('sage.combinat.tableau', 'StandardTableaux_n', StandardTableaux_size)
register_unpickle_override('sage.combinat.tableau', 'StandardTableaux_partition', StandardTableaux_shape)
register_unpickle_override('sage.combinat.tableau', 'SemistandardTableaux_n', SemistandardTableaux_size)
register_unpickle_override('sage.combinat.tableau', 'SemistandardTableaux_p', SemistandardTableaux_shape)
register_unpickle_override('sage.combinat.tableau', 'SemistandardTableaux_nmu', SemistandardTableaux_size_weight)
register_unpickle_override('sage.combinat.tableau', 'SemistandardTableaux_pmu', SemistandardTableaux_shape_weight)
