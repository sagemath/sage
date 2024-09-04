# sage.doctest: needs sage.combinat sage.modules
"""
Shifted primed tableaux

AUTHORS:

- Kirill Paramonov (2017-08-18): initial implementation
- Chaman Agrawal (2019-08-12): add parameter to allow primed diagonal entry
"""

# ****************************************************************************
#       Copyright (C) 2017 Kirill Paramonov <kbparamonov at ucdavis.edu>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.partition import Partition, Partitions, _Partitions, OrderedPartitions
from sage.combinat.partitions import ZS1_iterator
from sage.combinat.tableau import Tableaux
from sage.combinat.skew_partition import SkewPartition
from sage.combinat.integer_vector import IntegerVectors
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import

from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.sage_object import SageObject

from sage.categories.regular_crystals import RegularCrystals
from sage.categories.regular_supercrystals import RegularSuperCrystals
from sage.categories.sets_cat import Sets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.combination import Combinations

lazy_import('sage.combinat.root_system.cartan_type', 'CartanType')


class ShiftedPrimedTableau(ClonableArray,
        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A shifted primed tableau.

    A primed tableau is a tableau of shifted shape in the alphabet
    `X' = \{1' < 1 < 2' < 2 < \cdots < n' < n\}` such that

    1. the entries are weakly increasing along rows and columns;
    2. a row cannot have two repeated primed elements, and a column
       cannot have two repeated non-primed elements;

    Skew shape of the shifted primed tableaux is specified either
    with an optional argument ``skew`` or with ``None`` entries.

    Primed entries in the main diagonal can be allowed with the optional
    boolean parameter ``primed_diagonal``(default: ``False``).

    EXAMPLES::

        sage: T = ShiftedPrimedTableaux([4,2])
        sage: T([[1,"2'","3'",3],[2,"3'"]])[1]
        (2, 3')
        sage: t = ShiftedPrimedTableau([[1,"2p",2.5,3],[2,2.5]])
        sage: t[1]
        (2, 3')
        sage: ShiftedPrimedTableau([["2p",2,3],["2p","3p"],[2]], skew=[2,1])
        [(None, None, 2', 2, 3), (None, 2', 3'), (2,)]
        sage: ShiftedPrimedTableau([[None,None,"2p"],[None,"2p"]])
        [(None, None, 2'), (None, 2')]
        sage: T = ShiftedPrimedTableaux([4,2], primed_diagonal=True)
        sage: T([[1,"2'","3'",3],["2'","3'"]])[1] # With primed diagonal entry
        (2', 3')

    TESTS::

        sage: t = ShiftedPrimedTableau([[1,2,2.5,3],[2,2.5]])
        Traceback (most recent call last):
        ...
        ValueError: [[1, 2, 2.50000000000000, 3], [2, 2.50000000000000]]
         is not an element of Shifted Primed Tableaux

        sage: ShiftedPrimedTableau([[1,1,2.5],[1.5,2.5]])
        Traceback (most recent call last):
        ...
        ValueError: [[1, 1, 2.50000000000000], [1.50000000000000, 2.50000000000000]]
         is not an element of Shifted Primed Tableaux

        sage: ShiftedPrimedTableau([[1,1,2.5],[1.5,2.5]], primed_diagonal=True)
        [(1, 1, 3'), (2', 3')]
    """
    @staticmethod
    def __classcall_private__(cls, T, skew=None, primed_diagonal=False):
        r"""
        Ensure that a shifted tableau is only ever constructed as an
        ``element_class`` call of an appropriate parent.

        EXAMPLES::

            sage: data = [[1,"2'","2",3],[2,"3'"]]
            sage: t = ShiftedPrimedTableau(data)
            sage: T = ShiftedPrimedTableaux(shape=[4,2],weight=(1,3,2))
            sage: t == T(data)
            True
            sage: S = ShiftedPrimedTableaux(shape=[4,2])
            sage: t == S(data)
            True
            sage: t = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: t.parent()
            Shifted Primed Tableaux skewed by [2, 1]
            sage: s = ShiftedPrimedTableau([[None, None,"2p",2,3],[None,"2p"]])
            sage: s.parent()
            Shifted Primed Tableaux skewed by [2, 1]

        TESTS::

            sage: ShiftedPrimedTableau([])
            []
            sage: ShiftedPrimedTableau([tuple()])
            []
            sage: ShiftedPrimedTableau([], primed_diagonal=True)
            []
            sage: ShiftedPrimedTableau([tuple()], primed_diagonal=True)
            []
        """
        if (isinstance(T, ShiftedPrimedTableau) and T._skew == skew
            and T.parent()._primed_diagonal == primed_diagonal):
            return T

        skew_ = Partition([row.count(None) for row in T])
        if skew_:
            if skew and Partition(skew) != skew_:
                raise ValueError("skew shape does not agree with None entries")
            skew = skew_
        return ShiftedPrimedTableaux(skew=skew, primed_diagonal=primed_diagonal)(T)

    def __init__(self, parent, T, skew=None, check=True, preprocessed=False):
        r"""
        Initialize a shifted tableau.

        TESTS::

            sage: s = ShiftedPrimedTableau([[1,"2'","3'",3], [2,"3'"]])
            sage: t = ShiftedPrimedTableaux([4,2])([[1,"2p","3p",3], [2,"3p"]])
            sage: s == t
            True
            sage: t.parent()
            Shifted Primed Tableaux of shape [4, 2]
            sage: s.parent()
            Shifted Primed Tableaux
            sage: r = ShiftedPrimedTableaux([4, 2])(s); r.parent()
            Shifted Primed Tableaux of shape [4, 2]
            sage: s is t  # identical shifted tableaux are distinct objects
            False

        A shifted primed tableau is deeply immutable as the rows are
        stored as tuples::

            sage: t = ShiftedPrimedTableau([[1,"2p","3p",3],[2,"3p"]])
            sage: t[0][1] = 3
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
        """
        if not preprocessed:
            T = self._preprocess(T, skew=skew)
        self._skew = skew
        ClonableArray.__init__(self, parent, T, check=check)

    @staticmethod
    def _preprocess(T, skew=None):
        """
        Preprocessing list ``T`` to initialize the tableau.

        The output is a list of rows as tuples, with explicit
        ``None`` to indicate the skew shape, and entries being
        ``PrimedEntry`` instances.

        Trailing empty rows are removed.

        TESTS::

            sage: ShiftedPrimedTableau._preprocess([["2'", "3p", 3.5]],
            ....: skew=[1])
            [(None, 2', 3', 4')]
            sage: ShiftedPrimedTableau._preprocess([[None]], skew=[1])
            [(None,)]
            sage: ShiftedPrimedTableau._preprocess([], skew=[2,1])
            [(None, None), (None,)]
            sage: ShiftedPrimedTableau._preprocess([], skew=[])
            []
        """
        if isinstance(T, ShiftedPrimedTableau):
            return T
        # Preprocessing list t for primes and other symbols
        T = [[PrimedEntry(entry) for entry in row if entry is not None]
             for row in T]
        while T and not T[-1]:
            T = T[:-1]
        row_min = min(len(skew), len(T)) if skew else 0
        T_ = [(None,)*skew[i] + tuple(T[i]) for i in range(row_min)]

        if row_min < len(T):
            T_ += [tuple(T[i]) for i in range(row_min, len(T))]
        elif skew:
            T_ += [(None,)*skew[i] for i in range(row_min, len(skew))]
        return T_

    def check(self):
        """
        Check that ``self`` is a valid primed tableau.

        EXAMPLES::

            sage: T = ShiftedPrimedTableaux([4,2])
            sage: t = T([[1,'2p',2,2],[2,'3p']])
            sage: t.check()
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"],[2]],skew=[2,1])
            sage: s.check()
            sage: t = T([['1p','2p',2,2],[2,'3p']])
            Traceback (most recent call last):
            ...
            ValueError: [['1p', '2p', 2, 2], [2, '3p']] is not an element of
             Shifted Primed Tableaux of shape [4, 2]

            sage: T = ShiftedPrimedTableaux([4,2], primed_diagonal=True)
            sage: t = T([['1p','2p',2,2],[2,'3p']])  # primed_diagonal allowed
            sage: t.check()
            sage: t = T([['1p','1p',2,2],[2,'3p']])
            Traceback (most recent call last):
            ...
            ValueError: [['1p', '1p', 2, 2], [2, '3p']] is not an element of
             Shifted Primed Tableaux of shape [4, 2] and maximum entry 6
        """
        if not self.parent()._contains_tableau(self):
            raise ValueError("{} is not an element of Shifted Primed Tableaux".format(self))

    def is_standard(self):
        r"""
        Return ``True`` if the entries of ``self`` are in bijection with
        positive primed integers `1', 1, 2', \ldots, n`.

        EXAMPLES::

            sage: ShiftedPrimedTableau([["1'", 1, "2'"], [2, "3'"]],
            ....:                      primed_diagonal=True).is_standard()
            True
            sage: ShiftedPrimedTableau([["1'", 1, 2], ["2'", "3'"]],
            ....:                      primed_diagonal=True).is_standard()
            True
            sage: ShiftedPrimedTableau([["1'", 1, 1], ["2'", 2]],
            ....:                      primed_diagonal=True).is_standard()
            False
            sage: ShiftedPrimedTableau([[1, "2'"], [2]]).is_standard()
            False
            sage: s = ShiftedPrimedTableau([[None, None,"1p","2p",2],[None,"1"]])
            sage: s.is_standard()
            True
        """
        flattened = set([i for row in self for i in row if i is not None])
        if len(flattened) != sum(len(row) - row.count(None) for row in self):
            return False

        a = PrimedEntry('1p')
        while flattened:
            if a not in flattened:
                return False
            flattened.remove(a)
            a = a.increase_half()
        return True

    def __eq__(self, other):
        """
        Check whether ``self`` is equal to ``other``.

        INPUT:

        - ``other`` -- the element that ``self`` is compared to

        OUTPUT: boolean

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,"2p"]])
            sage: t == ShiftedPrimedTableaux([2])([[1,3/2]])
            True
            sage: s = ShiftedPrimedTableau([["2p",3]], skew=[1])
            sage: s == [[None, "2p", 3]]
            True
        """
        if isinstance(other, ShiftedPrimedTableau):
            return self._skew == other._skew and list(self) == list(other)
        try:
            Tab = ShiftedPrimedTableau(other, primed_diagonal=self.parent()._primed_diagonal)
        except (ValueError, TypeError):
            return False
        return self._skew == Tab._skew and list(self) == list(Tab)

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        INPUT:

        - ``other`` -- the element that ``self`` is compared to

        OUTPUT: boolean

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,"2p"]])
            sage: t != ShiftedPrimedTableaux([2])([[1,1]])
            True
            sage: s = ShiftedPrimedTableau([["2p",3]], skew=[1])
            sage: s != [[None, "2p", 3]]
            False
        """
        return not (self == other)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,"2p"]])
            sage: hash(t) == hash(ShiftedPrimedTableaux([2])([[1,3/2]]))
            True
        """
        return hash((self._skew, tuple(self)))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t
            [(1, 2', 2, 2), (2, 3')]
            sage: ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            [(None, None, 2', 2, 3), (None, 2')]
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def _repr_list(self):
        """
        Return a string representation of ``self`` as a list of tuples.

        EXAMPLES::

            sage: ShiftedPrimedTableau([['2p',3],[2,2]], skew=[2])._repr_list()
            "[(None, None, 2', 3), (2, 2)]"
        """
        return repr(list(self))

    def _repr_tab(self):
        """
        Return a nested list of strings representing the elements.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t._repr_tab()
            [[' 1 ', " 2'", ' 2 ', ' 2 '], [' 2 ', " 3'"]]
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: s._repr_tab()
            [[' . ', ' . ', " 2'", ' 2 ', ' 3 '], [' . ', " 2'"]]
        """
        max_len = len(str(self.max_entry())) + 2
        repr_tab = []
        for row in self:
            repr_row = []
            for entry in row:
                if entry is None:
                    repr_row.append('. '.rjust(max_len))
                elif entry.is_primed():
                    repr_row.append(repr(entry).rjust(max_len))
                elif entry.is_unprimed():
                    repr_row.append(repr(entry).rjust(max_len-1)+" ")
            repr_tab.append(repr_row)
        return repr_tab

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as an array.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: print(t._repr_diagram())
             1  2' 2  2
                2  3'
            sage: t = ShiftedPrimedTableau([[10,'11p',11,11],[11,'12']])
            sage: print(t._repr_diagram())
             10  11' 11  11
                 11  12
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: print(s._repr_diagram())
             .  .  2' 2  3
                .  2'
        """
        max_len = len(str(self.max_entry()))+2
        return "\n".join([" "*max_len*i + "".join(val)
                          for i, val in enumerate(self._repr_tab())])

    _repr_compact = _repr_diagram

    def _ascii_art_(self):
        """
        Return ASCII representation of ``self``.

        EXAMPLES::

            sage: ascii_art(ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']]))
            +---+---+---+---+
            | 1 | 2'| 2 | 2 |
            +---+---+---+---+
                | 2 | 3'|
                +---+---+
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: ascii_art(s)
            +---+---+---+---+---+
            | . | . | 2'| 2 | 3 |
            +---+---+---+---+---+
                | . | 2'|
                +---+---+

        TESTS::

            sage: ascii_art(ShiftedPrimedTableau([]))
            ++
            ++

            sage: ascii_art(ShiftedPrimedTableau([], skew=[1]))
            +---+
            | . |
            +---+
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self._ascii_art_table(unicode=False).splitlines())

    def _unicode_art_(self):
        """
        Return a Unicode representation of ``self``.

        EXAMPLES::

            sage: unicode_art(ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']]))
            ┌───┬───┬───┬───┐
            │ 1 │ 2'│ 2 │ 2 │
            └───┼───┼───┼───┘
                │ 2 │ 3'│
                └───┴───┘
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: unicode_art(s)
            ┌───┬───┬───┬───┬───┐
            │ . │ . │ 2'│ 2 │ 3 │
            └───┼───┼───┼───┴───┘
                │ . │ 2'│
                └───┴───┘

        TESTS::

            sage: unicode_art(ShiftedPrimedTableau([]))
            ┌┐
            └┘
            sage: unicode_art(ShiftedPrimedTableau([], skew=[1]))
            ┌───┐
            │ . │
            └───┘
        """
        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt(self._ascii_art_table(unicode=True).splitlines())

    def _ascii_art_table(self, unicode=False):
        """
        TESTS::

            sage: t = ShiftedPrimedTableau([[1,'2p',2],[2,'3p']])
            sage: print(t._ascii_art_table(unicode=True))
            ┌───┬───┬───┐
            │ 1 │ 2'│ 2 │
            └───┼───┼───┤
                │ 2 │ 3'│
                └───┴───┘
            sage: print(t._ascii_art_table())
            +---+---+---+
            | 1 | 2'| 2 |
            +---+---+---+
                | 2 | 3'|
                +---+---+
            sage: s = ShiftedPrimedTableau([[1,'2p',2, 23],[2,'30p']])
            sage: print(s._ascii_art_table(unicode=True))
            ┌────┬────┬────┬────┐
            │  1 │  2'│  2 │ 23 │
            └────┼────┼────┼────┘
                 │  2 │ 30'│
                 └────┴────┘
            sage: print(s._ascii_art_table(unicode=False))
            +----+----+----+----+
            |  1 |  2'|  2 | 23 |
            +----+----+----+----+
                 |  2 | 30'|
                 +----+----+
            sage: s = ShiftedPrimedTableau([["2p",2,10],["2p"]],skew=[2,1])
            sage: print(s._ascii_art_table(unicode=True))
            ┌────┬────┬────┬────┬────┐
            │  . │  . │  2'│  2 │ 10 │
            └────┼────┼────┼────┴────┘
                 │  . │  2'│
                 └────┴────┘
        """
        if unicode:
            import unicodedata
            v = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL')
            h = unicodedata.lookup('BOX DRAWINGS LIGHT HORIZONTAL')
            dl = unicodedata.lookup('BOX DRAWINGS LIGHT DOWN AND LEFT')
            dr = unicodedata.lookup('BOX DRAWINGS LIGHT DOWN AND RIGHT')
            ul = unicodedata.lookup('BOX DRAWINGS LIGHT UP AND LEFT')
            ur = unicodedata.lookup('BOX DRAWINGS LIGHT UP AND RIGHT')
            vl = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL AND LEFT')
            uh = unicodedata.lookup('BOX DRAWINGS LIGHT UP AND HORIZONTAL')
            dh = unicodedata.lookup('BOX DRAWINGS LIGHT DOWN AND HORIZONTAL')
            vh = unicodedata.lookup(
                'BOX DRAWINGS LIGHT VERTICAL AND HORIZONTAL')
        else:
            v = '|'
            h = '-'
            dl = dr = ul = ur = vl = uh = dh = vh = '+'

        if not self.shape():
            return dr + dl + '\n' + ur + ul

        # Get the widths of the columns
        str_tab = self._repr_tab()
        width = len(str_tab[0][0])
        str_list = [dr + (h*width + dh)*(len(str_tab[0])-1) + h*width + dl]
        for nrow, row in enumerate(str_tab):
            l1 = " " * (width+1) * nrow
            l2 = " " * (width+1) * nrow
            n = len(str_tab[nrow+1]) if nrow+1 < len(str_tab) else -1
            for i, e in enumerate(row):
                if i == 0:
                    l1 += ur + h*width
                elif i <= n+1:
                    l1 += vh + h*width
                else:
                    l1 += uh + h*width
                if unicode:
                    l2 += "{}{:^{width}}".format(v, e, width=width)
                else:
                    l2 += "{}{:^{width}}".format(v, e, width=width)
            if i <= n:
                l1 += vl
            else:
                l1 += ul
            l2 += v
            str_list.append(l2)
            str_list.append(l1)
        return "\n".join(str_list)

    def pp(self):
        """
        Pretty print ``self``.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.pp()
            1  2' 2  2
               2  3'
            sage: t = ShiftedPrimedTableau([[10,'11p',11,11],[11,'12']])
            sage: t.pp()
            10  11' 11  11
                11  12
            sage: s = ShiftedPrimedTableau([['2p',2,3],['2p']],skew=[2,1])
            sage: s.pp()
             .  .  2' 2  3
                .  2'

        TESTS::

            sage: ShiftedPrimedTableau([],skew=[1]).pp()
            .
            sage: ShiftedPrimedTableau([]).pp()
            <BLANKLINE>
        """
        print(self._repr_diagram())

    def _latex_(self):
        r"""
        Return LaTex code for ``self``.

        EXAMPLES::

            sage: T = ShiftedPrimedTableaux([4,2])
            sage: latex(T([[1,"2p",2,"3p"],[2,3]]))
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{1-4}
            \lr{ 1 }&\lr{ 2'}&\lr{ 2 }&\lr{ 3'}\\\cline{1-4}
            &\lr{ 2 }&\lr{ 3 }\\\cline{2-3}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        L = [[None]*i + row for i, row in enumerate(self._repr_tab())]
        return tex_from_array(L)

    def max_entry(self):
        r"""
        Return the minimum unprimed letter `x > y` for all `y` in ``self``.

        EXAMPLES::

            sage: Tab = ShiftedPrimedTableau([(1,1,'2p','3p'),(2,2)])
            sage: Tab.max_entry()
            3

        TESTS::

            sage: Tab = ShiftedPrimedTableau([], skew=[2,1])
            sage: Tab.max_entry()
            0
            sage: Tab = ShiftedPrimedTableau([["1p"]], skew=[2,1])
            sage: Tab.max_entry()
            1
        """
        flat = [entry.unprimed() for row in self
                for entry in row if entry is not None]
        if len(flat) == 0:
            return 0
        return max(flat)

    def shape(self):
        r"""
        Return the shape of the underlying partition of ``self``.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.shape()
            [4, 2]
            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]],skew=[2,1])
            sage: s.shape()
            [5, 2] / [2, 1]
        """
        if self._skew is None:
            return Partition([len(row) for row in self])
        return SkewPartition(([len(row) for row in self], self._skew))

    def restrict(self, n):
        """
        Return the restriction of the shifted tableau to all
        the numbers less than or equal to ``n``.

        .. NOTE::

            If only the outer shape of the restriction, rather than
            the whole restriction, is needed, then the faster method
            :meth:`restriction_outer_shape` is preferred. Similarly if
            only the skew shape is needed, use :meth:`restriction_shape`.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.restrict(2).pp()
            1  2' 2  2
               2

            sage: t.restrict("2p").pp()
            1  2'

            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]], skew=[2,1])
            sage: s.restrict(2).pp()
            .  .  2' 2
               .  2'
            sage: s.restrict(1.5).pp()
            .  .  2'
               .  2'
        """
        t = self[:]
        n = PrimedEntry(n)
        loop = ([y for y in x if y is not None and y <= n] for x in t)
        return ShiftedPrimedTableau([z for z in loop if z], skew=self._skew)

    def restriction_outer_shape(self, n):
        """
        Return the outer shape of the restriction of the shifted
        tableau ``self`` to `n`.

        If `T` is a (skew) shifted tableau and `n` is a half-integer,
        then the restriction of `T` to `n` is defined as the (skew)
        shifted tableau obtained by removing all cells filled with
        entries greater than `n` from `T`.

        This method computes merely the outer shape of the restriction.
        For the restriction itself, use :meth:`restrict`.

        EXAMPLES::

            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]], skew=[2,1])
            sage: s.pp()
            .  .  2' 2  3
               .  2'
            sage: s.restriction_outer_shape(2)
            [4, 2]
            sage: s.restriction_outer_shape("2p")
            [3, 2]
        """
        n = PrimedEntry(n)
        if self._skew is None:
            res = [len([y for y in row if y <= n]) for row in self]
        else:
            res = [len([y for y in row if y is None or y <= n])
                   for i, row in enumerate(self)]

        return Partition(res)

    def restriction_shape(self, n):
        """
        Return the skew shape of the restriction of the skew tableau
        ``self`` to ``n``.

        If `T` is a shifted tableau and `n` is a half-integer, then
        the restriction of `T` to `n` is defined as the
        (skew) shifted tableau obtained by removing all cells
        filled with entries greater than `n` from `T`.

        This method computes merely the skew shape of the restriction.
        For the restriction itself, use :meth:`restrict`.

        EXAMPLES::

            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]], skew=[2,1])
            sage: s.pp()
             .  .  2' 2  3
                .  2'

            sage: s.restriction_shape(2)
            [4, 2] / [2, 1]
        """
        if self._skew is None:
            return Partition(self.restriction_outer_shape(n))
        else:
            return SkewPartition([self.restriction_outer_shape(n), self._skew])

    def to_chain(self):
        """
        Return the chain of partitions corresponding to the (skew)
        shifted tableau ``self``, interlaced by one of the colours
        ``1`` is the added cell is on the diagonal, ``2`` if an
        ordinary entry is added and ``3`` if a primed entry is added.

        EXAMPLES::

            sage: s = ShiftedPrimedTableau([(1, 2, 3.5, 5, 6.5), (3, 5.5)])
            sage: s.pp()
             1  2  4' 5  7'
                3  6'

            sage: s.to_chain()
            [[], 1, [1], 2, [2], 1, [2, 1], 3, [3, 1], 2, [4, 1], 3, [4, 2], 3, [5, 2]]


            sage: s = ShiftedPrimedTableau([(1, 3.5), (2.5,), (6,)], skew=[2,1])
            sage: s.pp()
             .  .  1  4'
                .  3'
                   6

            sage: s.to_chain()
            [[2, 1], 2, [3, 1], 0, [3, 1], 3, [3, 2], 3, [4, 2], 0, [4, 2], 1, [4, 2, 1]]

        TESTS::

            sage: s = ShiftedPrimedTableau([["2p",2,3],["2p"]], skew=[2,1])
            sage: s.pp()
             .  .  2' 2  3
                .  2'
            sage: s.to_chain()
            Traceback (most recent call last):
            ...
            ValueError: can compute a chain of partitions only for
             skew shifted tableaux without repeated entries
        """
        if any(e not in [0, 1] for e in self.weight()):
            raise ValueError("can compute a chain of partitions only for skew"
                             " shifted tableaux without repeated entries")
        entries = sorted(e for row in self for e in row if e is not None)
        if self._skew is None:
            mu = Partition([])
            m = 0
        else:
            mu = self._skew
            m = len(self._skew)
        chain = [mu]
        f = 0
        for e in entries:
            n = e.integer()
            chain.extend([0, mu]*int(n-f-1))
            mu = self.restriction_outer_shape(e)
            if n == e:
                if any(e == row[0] for i, row in enumerate(self) if i >= m or self._skew[i] == 0):
                    chain.append(1)
                else:
                    chain.append(2)
            else:
                chain.append(3)
            chain.append(mu)
            f = n
        return chain

    def weight(self):
        r"""
        Return the weight of ``self``.

        The weight of a shifted primed tableau is defined to be the vector
        with `i`-th component equal to the number of entries `i` and `i'`
        in the tableau.

        EXAMPLES::

           sage: t = ShiftedPrimedTableau([['2p',2,2],[2,'3p']], skew=[1])
           sage: t.weight()
           (0, 4, 1)
        """
        flat = [entry.integer() for row in self
                for entry in row if entry is not None]
        if not flat:
            return ()

        weight = tuple([flat.count(i+1) for i in range(max(flat))])
        return weight


class CrystalElementShiftedPrimedTableau(ShiftedPrimedTableau):
    """
    Class for elements of ``crystals.ShiftedPrimedTableau``.
    """

    def _to_matrix(self):
        """
        Return a 2-dimensional array representation of a shifted tableau.

        EXAMPLES::

            sage: SPT = ShiftedPrimedTableaux([4,2,1])
            sage: t = SPT([[1,'2p',2,2],[2,'3p'],[3]])
            sage: mat = t._to_matrix()
            sage: mat
            [[1, 2', 2, 2], [None, 2, 3', None], [None, None, 3, None]]
        """
        m = len(self[0])
        return [[None]*i + list(row) + [None]*(m-i-len(row))
                for i, row in enumerate(self)]

    def _reading_word_with_positions(self):
        """
        Iterate over the reading word of ``self`` together with positions
        of the corresponding letters in ``self``.

        The reading word of a shifted primed tableau is constructed
        as follows:

        1. List all primed entries in the tableau, column by
           column, in decreasing order within each column, moving
           from the rightmost column to the left, and with all
           the primes removed (i.e. all entries are increased by
           half a unit).

        2. Then list all unprimed entries, row by row, in
           increasing order within each row, moving from the
           bottommost row to the top.

        EXAMPLES::

            sage: SPT = ShiftedPrimedTableaux([4,2])
            sage: t = SPT([[1,'2p',2,2],[2,'3p']])
            sage: list(t._reading_word_with_positions())
            [((1, 2), 3), ((0, 1), 2), ((1, 1), 2), ((0, 0), 1),
             ((0, 2), 2), ((0, 3), 2)]
        """
        mat = self._to_matrix()
        ndim, mdim = len(mat), len(mat[0])
        for j in reversed(range(mdim)):
            for i in range(ndim):
                x = mat[i][j]
                if x is not None and x.is_primed():
                    yield ((i, j), x.integer())
        for i in reversed(range(ndim)):
            for j in range(mdim):
                x = mat[i][j]
                if x is not None and x.is_unprimed():
                    yield ((i, j), x.integer())

    def reading_word(self):
        """
        Return the reading word of ``self``.

        The reading word of a shifted primed tableau is constructed
        as follows:

        1. List all primed entries in the tableau, column by
           column, in decreasing order within each column, moving
           from the rightmost column to the left, and with all
           the primes removed (i.e. all entries are increased by
           half a unit).

        2. Then list all unprimed entries, row by row, in
           increasing order within each row, moving from the
           bottommost row to the top.

        EXAMPLES::

            sage: SPT = ShiftedPrimedTableaux([4,2])
            sage: t = SPT([[1,'2p',2,2],[2,'3p']])
            sage: t.reading_word()
            [3, 2, 2, 1, 2, 2]
        """
        if self._skew is not None:
            raise NotImplementedError('skew tableau must be empty')
        return [tup[1] for tup in self._reading_word_with_positions()]

    def f(self, ind):
        r"""
        Compute the action of the crystal operator `f_i` on a shifted primed
        tableau using cases from the papers [HPS2017]_ and [AO2018]_.

        INPUT:

        - ``ind`` -- element in the index set of the crystal

        OUTPUT: primed tableau or ``None``

        EXAMPLES::

            sage: SPT = ShiftedPrimedTableaux([5,4,2])
            sage: t = SPT([[1,1,1,1,'3p'],[2,2,2,'3p'],[3,3]])
            sage: t.pp()
            1  1  1  1  3'
               2  2  2  3'
                  3  3
            sage: s = t.f(2)
            sage: s is None
            True

            sage: t = SPT([[1,1,1,'2p','3p'],[2,2,3,3],[3,4]])
            sage: t.pp()
            1  1  1  2' 3'
               2  2  3  3
                  3  4
            sage: s = t.f(2)
            sage: s.pp()
            1  1  1  2' 3'
               2  3' 3  3
                  3  4

            sage: SPT = ShiftedPrimedTableaux([2,1])
            sage: t = SPT([[1,1],[2]])
            sage: t.f(-1).pp()
            1  2'
               2
            sage: t.f(1).pp()
            1  2'
               2
            sage: t.f(2).pp()
            1  1
               3

            sage: r = SPT([[1,'2p'],[2]])
            sage: r.f(-1) is None
            True
            sage: r.f(1) is None
            True
            sage: r.f(2).pp()
            1  2'
               3

            sage: r = SPT([[1,1],[3]])
            sage: r.f(-1).pp()
            1  2'
               3
            sage: r.f(1).pp()
            1  2
               3
            sage: r.f(2) is None
            True

            sage: r = SPT([[1,2],[3]])
            sage: r.f(-1).pp()
            2  2
               3
            sage: r.f(1).pp()
            2  2
               3
            sage: r.f(2) is None
            True

            sage: t = SPT([[1,1],[2]])
            sage: t.f(-1).f(2).f(2).f(-1) == t.f(2).f(1).f(-1).f(2)
            True
            sage: t.f(-1).f(2).f(2).f(-1).pp()
            2  3'
               3
            sage: all(t.f(-1).f(2).f(2).f(-1).f(i) is None for i in {-1, 1, 2})
            True

            sage: SPT = ShiftedPrimedTableaux([4])
            sage: t = SPT([[1,1,1,1]])
            sage: t.f(-1).pp()
            1  1  1  2'
            sage: t.f(1).pp()
            1  1  1  2
            sage: t.f(-1).f(-1) is None
            True
            sage: t.f(1).f(-1).pp()
            1  1  2' 2
            sage: t.f(1).f(1).pp()
            1  1  2  2
            sage: t.f(1).f(1).f(-1).pp()
            1  2' 2  2
            sage: t.f(1).f(1).f(1).pp()
            1  2  2  2
            sage: t.f(1).f(1).f(1).f(-1).pp()
            2  2  2  2
            sage: t.f(1).f(1).f(1).f(1).pp()
            2  2  2  2
            sage: t.f(1).f(1).f(1).f(1).f(-1) is None
            True
        """
        T = self._to_matrix()

        # special logic for queer lowering operator f_{-1}
        if ind == -1:
            read_word = [num for num in self._reading_word_with_positions() if num[1] in {1, 2}]

            # f_{-1} acts as zero if tableau contains 2'
            if any(elt == 2 and T[pos[0]][pos[1]].is_primed() for pos, elt in read_word):
                return None

            ones = sorted([pos for pos, elt in read_word if elt == 1], key=lambda x: x[1])

            # f_{-1} acts as zero if tableau contains no entries equal to 1
            if len(ones) == 0:
                return None
            # otherwise, f_{-1} changes last 1 in first row to 2' (if off diagonal) or 2 (if on diagonal)
            else:
                r, c = ones[-1]
                assert r == 0
                T[r][c] = PrimedEntry('2p') if r != c else PrimedEntry(2)
                T = [tuple(elmt for elmt in row if elmt is not None) for row in T]
                return type(self)(self.parent(), T, check=False, preprocessed=True)

        read_word = [num for num in self._reading_word_with_positions()
                     if num[1] == ind or num[1] == ind+1]

        element_to_change = None
        count = 0

        for element in read_word:
            if element[1] == ind+1:
                count += 1
            elif count == 0:
                element_to_change = element
            else:
                count -= 1

        if element_to_change is None:
            return None

        (r, c), elt = element_to_change
        ind_e = PrimedEntry(ind)
        ind_plus_one = ind_e.increase_one()
        ind_plus_half = ind_e.increase_half()

        if T[r][c].is_primed():
            T = [[elmt.increase_half() if elmt is not None else elmt
                  for elmt in row] for row in T]
            T = [list(z) for z in zip(*T)]
            r, c = c, r
        h, l = len(T), len(T[0])

        if (c+1 == l or T[r][c+1] is None or T[r][c+1] >= ind_plus_one):
            (tp_r, tp_c) = (r, c)
            while True:
                if tp_r+1 == h or T[tp_r+1][tp_c] is None or T[tp_r+1][tp_c] > ind_plus_one:
                    break
                if tp_r <= tp_c and T[tp_r+1][tp_r+1] == ind_plus_one:
                    tp_r += 1
                    tp_c = tp_r
                    break
                if ind_plus_half not in T[tp_r+1]:
                    break
                tp_r += 1
                tp_c = T[tp_r].index(ind_plus_half)

            if tp_r == r:
                T[r][c] = T[r][c].increase_one()
            elif tp_r == tp_c:
                T[r][c] = T[r][c].increase_half()
            else:
                T[r][c] = T[r][c].increase_half()
                T[tp_r][tp_c] = T[tp_r][tp_c].increase_half()

        elif T[r][c+1] == ind_plus_half:
            T[r][c+1] = T[r][c+1].increase_half()
            T[r][c] = T[r][c].increase_half()

        if r > c:
            T = [[elmt.decrease_half() if elmt is not None else elmt
                  for elmt in row] for row in T]
            T = [list(z) for z in zip(*T)]

        T = [tuple(elmt for elmt in row if elmt is not None) for row in T]
        return type(self)(self.parent(), T, check=False, preprocessed=True)

    def e(self, ind):
        r"""
        Compute the action of the crystal operator `e_i` on a shifted primed
        tableau using cases from the papers [HPS2017]_ and [AO2018]_.

        INPUT:

        - ``ind`` -- an element in the index set of the crystal

        OUTPUT: primed tableau or ``None``

        EXAMPLES::

            sage: SPT = ShiftedPrimedTableaux([5,4,2])
            sage: t = SPT([[1,1,1,'2p','3p'], [2,'3p',3,3],[3,4]])
            sage: t.pp()
            1  1  1  2' 3'
               2  3' 3  3
                  3  4
            sage: s = t.e(2)
            sage: s.pp()
            1  1  1  2' 3'
               2  2  3  3
                  3  4
            sage: t == s.f(2)
            True

            sage: SPT = ShiftedPrimedTableaux([2,1])
            sage: t = SPT([[2,'3p'],[3]])
            sage: t.e(-1).pp()
            1  3'
               3
            sage: t.e(1).pp()
            1  3'
               3
            sage: t.e(2).pp()
            2  2
               3

            sage: r = SPT([[2, 2],[3]])
            sage: r.e(-1).pp()
            1  2
               3
            sage: r.e(1).pp()
            1  2
               3
            sage: r.e(2) is None
            True

            sage: r = SPT([[1,'3p'],[3]])
            sage: r.e(-1) is None
            True
            sage: r.e(1) is None
            True
            sage: r.e(2).pp()
            1  2'
               3
            sage: r = SPT([[1,'2p'],[3]])
            sage: r.e(-1).pp()
            1  1
               3
            sage: r.e(1) is None
            True
            sage: r.e(2).pp()
            1  2'
               2
            sage: t = SPT([[2,'3p'],[3]])
            sage: t.e(-1).e(2).e(2).e(-1) == t.e(2).e(1).e(1).e(2)
            True
            sage: t.e(-1).e(2).e(2).e(-1).pp()
            1  1
               2
            sage: all(t.e(-1).e(2).e(2).e(-1).e(i) is None for i in {-1, 1, 2})
            True

            sage: SPT = ShiftedPrimedTableaux([4])
            sage: t = SPT([[2,2,2,2]])
            sage: t.e(-1).pp()
            1  2  2  2
            sage: t.e(1).pp()
            1  2  2  2
            sage: t.e(-1).e(-1) is None
            True
            sage: t.e(1).e(1).pp()
            1  1  2  2
        """
        T = self._to_matrix()

        # special logic for queer raising operator e_{-1}
        if ind == -1:
            read_word = [num for num in self._reading_word_with_positions() if num[1] in {1, 2}]

            two_primes = sorted(
                [pos for pos, elt in read_word if elt == 2 and T[pos[0]][pos[1]].is_primed()],
                key=lambda x: x[1]
            )

            # e_{-1} acts as zero if tableau contains no 2' and first diagonal entry is not 2
            if len(two_primes) == 0:
                if T[0][0] != PrimedEntry(2):
                    return None
                # if tableau has no 2' but first diagonal entry is 2, then e_{-1} changes this to 1
                else:
                    T[0][0] = PrimedEntry(1)
            # if tableau has at least one 2', then e_{-1} changes the first 2' to 1
            else:
                r, c = two_primes[0]
                assert r == 0
                T[r][c] = PrimedEntry(1)

            T = [tuple(elmt for elmt in row if elmt is not None) for row in T]
            return type(self)(self.parent(), T, check=False, preprocessed=True)

        read_word = [num for num in self._reading_word_with_positions()
                     if num[1] == ind or num[1] == ind+1]

        element_to_change = None
        count = 0

        for element in reversed(read_word):
            if element[1] == ind:
                count += 1
            elif count == 0:
                element_to_change = element
            else:
                count -= 1

        if element_to_change is None:
            return None
        (r, c), elt = element_to_change

        ind_e = PrimedEntry(ind)
        ind_plus_half = ind_e.increase_half()

        if T[r][c].is_primed():
            T = [[elmt.increase_half() if elmt is not None else elmt
                  for elmt in row] for row in T]
            T = [list(z) for z in zip(*T)]
            r, c = c, r

        if (c == 0 or T[r][c-1] is None or T[r][c-1] <= ind_e):
            (tp_r, tp_c) = (r, c)
            while True:
                if tp_r == 0 or T[tp_r-1][tp_c] is None or T[tp_r-1][tp_c] < ind_e:
                    break
                if ind_plus_half not in T[tp_r-1]:
                    break
                tp_r -= 1
                tp_c = T[tp_r].index(ind_plus_half)

            if tp_r == r:
                T[r][c] = T[r][c].decrease_one()
            elif tp_r == tp_c:
                T[r][c] = T[r][c].decrease_half()
            else:
                T[r][c] = T[r][c].decrease_half()
                T[tp_r][tp_c] = T[tp_r][tp_c].decrease_half()

        elif T[r][c-1] == ind_plus_half:
            T[r][c-1] = T[r][c-1].decrease_half()
            T[r][c] = T[r][c].decrease_half()
        if r > c:
            T = [[elmt.decrease_half() if elmt is not None else elmt
                  for elmt in row] for row in T]
            T = [list(z) for z in zip(*T)]

        T = [tuple(elmt for elmt in row if elmt is not None) for row in T]
        return type(self)(self.parent(), T, check=False, preprocessed=True)

    def is_highest_weight(self, index_set=None):
        r"""
        Return whether ``self`` is a highest weight element of the crystal.

        An element is highest weight if it vanishes under all crystal
        operators `e_i`.

        EXAMPLES::

            sage: SPT = ShiftedPrimedTableaux([5,4,2])
            sage: t = SPT([(1, 1, 1, 1, 1), (2, 2, 2, "3p"), (3, 3)])
            sage: t.is_highest_weight()
            True

            sage: SPT = ShiftedPrimedTableaux([5,4])
            sage: s = SPT([(1, 1, 1, 1, 1), (2, 2, "3p", 3)])
            sage: s.is_highest_weight(index_set=[1])
            True
        """
        read_w = self.reading_word()
        max_entry = max(read_w)
        count = {i: 0 for i in range(max_entry+1)}
        if index_set is None:
            index_set = self.parent().index_set()
        for l in reversed(read_w):
            count[l] += 1
            if l-1 in index_set and l > 1 and count[l] > count[l-1]:
                return False
        return True

    def weight(self):
        r"""
        Return the weight of ``self``.

        The weight of a shifted primed tableau is defined to be the vector
        with `i`-th component equal to the number of entries `i` and `i'`
        in the tableau.

        EXAMPLES::

           sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
           sage: t.weight()
           (1, 4, 1)
        """
        flat = [entry.integer() for row in self for entry in row]
        if flat == []:
            max_ind = 0
        else:
            max_ind = max(flat)
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        return self.parent().weight_lattice_realization()(weight)


class PrimedEntry(SageObject):
    r"""
    The class of entries in shifted primed tableaux.

    An entry in a shifted primed tableau is an element in
    the alphabet `\{1' < 1 < 2' < 2 < \cdots < n' < n\}`.
    The difference between two elements `i` and `i-1` counts as a
    whole unit, whereas the difference between `i` and `i'` counts
    as half a unit.
    Internally, we represent an unprimed element `x` as `2x`
    and the primed elements as the corresponding odd integer
    that respects the total order.

    INPUT:

    - ``entry`` -- half integer or string of an integer
      possibly ending in ``p`` or ``'``
    - ``double`` -- the doubled value
    """

    def __init__(self, entry=None, double=None):
        """
        Normalize the entry.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: PrimedEntry(2)
            2
            sage: PrimedEntry("2p")
            2'
            sage: PrimedEntry("2'")
            2'
            sage: a = PrimedEntry(2.5)
            sage: PrimedEntry(a)
            3'
            sage: PrimedEntry(None)
            Traceback (most recent call last):
            ...
            ValueError: primed entry must not be None
        """
        # store primed numbers as odd, unprimed numbers as even integers
        if isinstance(entry, self.__class__):
            self._entry = entry._entry
            return

        if double is not None:
            self._entry = Integer(double)
            return

        if isinstance(entry, str):
            if (entry[-1] == "'" or entry[-1] == "p") and entry[:-1].isdigit():
                # Check if an element has "'" or "p" at the end
                self._entry = 2 * Integer(entry[:-1]) - 1
            else:
                self._entry = 2 * Integer(entry)
            return

        if entry is None:
            raise ValueError("primed entry must not be None")
        try:
            self._entry = Integer(2*entry)
        except (TypeError, ValueError):
            raise ValueError("primed entries must be half-integers")

    def __hash__(self):
        """
        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry("2p")
            sage: b = PrimedEntry("2'")
            sage: a == b
            True
        """
        return hash(self._entry)

    def __repr__(self):
        """
        Represent ``self`` as primed or unprimed integer.

        TESTS::

            sage: ShiftedPrimedTableau([[1,"2p"]])[0][1]
            2'
        """
        if self.is_unprimed():
            return repr(self._entry // 2)
        else:
            return repr((self._entry+1) // 2) + "'"

    def integer(self):
        """
        Return the corresponding integer `i` for primed entries
        of the form `i` or `i'`.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: b = PrimedEntry("2p").integer()
            sage: b
            2
            sage: b.category()
            Category of elements of Integer Ring
        """
        return (self._entry + 1) // 2

    def __eq__(self, other):
        """
        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry("2p")
            sage: b = PrimedEntry("2'")
            sage: a == b
            True
        """
        try:
            other = PrimedEntry(other)
        except ValueError:
            return False
        return self._entry == other._entry

    def __ne__(self, other):
        """
        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry("1")
            sage: b = PrimedEntry(1)
            sage: a != b
            False
        """
        try:
            other = PrimedEntry(other)
        except ValueError:
            return True
        return self._entry != other._entry

    def __lt__(self, other):
        """
        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry("2p")
            sage: b = PrimedEntry(2)
            sage: a < b
            True
        """
        return self._entry < PrimedEntry(other)._entry

    def __le__(self, other):
        """
        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry(2)
            sage: b = PrimedEntry("3p")
            sage: a <= b
            True
        """
        return self._entry <= PrimedEntry(other)._entry

    def __gt__(self, other):
        """
        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry("2p")
            sage: b = PrimedEntry(2)
            sage: b > a
            True
        """
        return self._entry > PrimedEntry(other)._entry

    def __ge__(self, other):
        """
        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry(2)
            sage: b = PrimedEntry("3p")
            sage: a >= b
            False
        """
        return self._entry >= PrimedEntry(other)._entry

    def is_unprimed(self):
        """
        Check if ``self`` is an unprimed element.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry("2p")
            sage: a.is_unprimed()
            False
        """
        return self._entry % 2 == 0

    def is_primed(self):
        """
        Check if ``self`` is a primed element.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry("3p")
            sage: a.is_primed()
            True
        """
        return self._entry % 2 == 1

    def unprimed(self):
        """
        Unprime ``self`` if it is a primed element.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry("2p")
            sage: a.unprimed()
            2
        """
        if self.is_unprimed():
            return self
        else:
            return PrimedEntry(double=self._entry + 1)

    def primed(self):
        """
        Prime ``self`` if it is an unprimed element.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry(1)
            sage: a.primed()
            1'
        """
        if self.is_unprimed():
            return PrimedEntry(double=self._entry - 1)
        else:
            return self

    def increase_half(self):
        """
        Increase ``self`` by half a unit.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry(1)
            sage: a.increase_half()
            2'
        """
        return PrimedEntry(double=self._entry + 1)

    def decrease_half(self):
        """
        Decrease ``self`` by half a unit.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry(1)
            sage: a.decrease_half()
            1'
        """
        return PrimedEntry(double=self._entry - 1)

    def increase_one(self):
        """
        Increase ``self`` by one unit.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry("2p")
            sage: a.increase_one()
            3'
        """
        return PrimedEntry(double=self._entry + 2)

    def decrease_one(self):
        """
        Decrease ``self`` by one unit.

        TESTS::

            sage: from sage.combinat.shifted_primed_tableau import PrimedEntry
            sage: a = PrimedEntry("2p")
            sage: a.decrease_one()
            1'
        """
        return PrimedEntry(double=self._entry - 2)


class ShiftedPrimedTableaux(UniqueRepresentation, Parent):
    r"""
    Return the combinatorial class of shifted primed tableaux subject
    to the constraints given by the arguments.

    A primed tableau is a tableau of shifted shape on the alphabet
    `X' = \{1' < 1 < 2' < 2 < \cdots < n' < n\}` such that

    1. the entries are weakly increasing along rows and columns

    2. a row cannot have two repeated primed entries, and a column
       cannot have two repeated non-primed entries

    INPUT:

    Valid optional keywords:

    - ``shape`` -- the (outer skew) shape of tableaux

    - ``weight`` -- the weight of tableaux

    - ``max_entry`` -- the maximum entry of tableaux

    - ``skew`` -- the inner skew shape of tableaux

    - ``primed_diagonal`` -- allow primed entries in main diagonal of tableaux

    The weight of a tableau is defined to be the vector with `i`-th
    component equal to the number of entries `i` and `i'` in the tableau.
    The sum of the coordinates in the weight vector must be equal to the
    number of entries in the partition.

    The ``shape`` and ``skew`` must be strictly decreasing partitions.
    The ``primed_diagonal`` is a boolean (default: ``False``).

    EXAMPLES::

        sage: SPT = ShiftedPrimedTableaux(weight=(1,2,2), shape=[3,2]); SPT
        Shifted Primed Tableaux of weight (1, 2, 2) and shape [3, 2]
        sage: SPT.list()
        [[(1, 2, 2), (3, 3)],
         [(1, 2', 3'), (2, 3)],
         [(1, 2', 3'), (2, 3')],
         [(1, 2', 2), (3, 3)]]
        sage: SPT = ShiftedPrimedTableaux(weight=(1,2,2), shape=[3,2],
        ....:                             primed_diagonal=True); SPT
        Shifted Primed Tableaux of weight (1, 2, 2) and shape [3, 2]
        sage: SPT.list()
        [[(1, 2, 2), (3, 3)],
         [(1, 2, 2), (3', 3)],
         [(1, 2', 3'), (2, 3)],
         [(1, 2', 3'), (2, 3')],
         [(1, 2', 3'), (2', 3)],
         [(1, 2', 3'), (2', 3')],
         [(1, 2', 2), (3, 3)],
         [(1, 2', 2), (3', 3)],
         [(1', 2, 2), (3, 3)],
         [(1', 2, 2), (3', 3)],
         [(1', 2', 3'), (2, 3)],
         [(1', 2', 3'), (2, 3')],
         [(1', 2', 3'), (2', 3)],
         [(1', 2', 3'), (2', 3')],
         [(1', 2', 2), (3, 3)],
         [(1', 2', 2), (3', 3)]]
        sage: SPT = ShiftedPrimedTableaux(weight=(1,2)); SPT
        Shifted Primed Tableaux of weight (1, 2)
        sage: list(SPT)
        [[(1, 2, 2)], [(1, 2', 2)], [(1, 2'), (2,)]]
        sage: SPT = ShiftedPrimedTableaux(weight=(1,2), primed_diagonal=True)
        sage: list(SPT)
        [[(1, 2, 2)],
         [(1, 2', 2)],
         [(1', 2, 2)],
         [(1', 2', 2)],
         [(1, 2'), (2,)],
         [(1, 2'), (2',)],
         [(1', 2'), (2,)],
         [(1', 2'), (2',)]]
        sage: SPT = ShiftedPrimedTableaux([3,2], max_entry=2); SPT
        Shifted Primed Tableaux of shape [3, 2] and maximum entry 2
        sage: list(SPT)
        [[(1, 1, 1), (2, 2)], [(1, 1, 2'), (2, 2)]]
        sage: SPT = ShiftedPrimedTableaux([3,2], max_entry=2,
        ....:                             primed_diagonal=True)
        sage: list(SPT)
        [[(1, 1, 1), (2, 2)],
         [(1, 1, 1), (2', 2)],
         [(1', 1, 1), (2, 2)],
         [(1', 1, 1), (2', 2)],
         [(1, 1, 2'), (2, 2)],
         [(1, 1, 2'), (2', 2)],
         [(1', 1, 2'), (2, 2)],
         [(1', 1, 2'), (2', 2)]]

    TESTS::

        sage: [(1,'2p',2,2),(2,'3p')] in ShiftedPrimedTableaux()
        True
        sage: [('1p','2p',2,2),(2,'3p')] in ShiftedPrimedTableaux()
        False
        sage: [(1,1),(2,2)] in ShiftedPrimedTableaux()
        False
        sage: [] in ShiftedPrimedTableaux()
        True
        sage: [('1p','2p',2,2),(2,'3p')] in ShiftedPrimedTableaux(
        ....:                                      primed_diagonal=True)
        True
        sage: [] in ShiftedPrimedTableaux(primed_diagonal=True)
        True

    .. SEEALSO::

        - :class:`ShiftedPrimedTableau`
    """
    Element = ShiftedPrimedTableau
    options = Tableaux.options

    @staticmethod
    def __classcall_private__(cls, shape=None, weight=None, max_entry=None,
                              skew=None, primed_diagonal=False):
        r"""
        Normalize and process input to return the correct parent and
        ensure a unique representation.

        TESTS::

            sage: ShiftedPrimedTableaux([])
            Shifted Primed Tableaux of shape []
            sage: ShiftedPrimedTableaux(3)
            Traceback (most recent call last):
            ...
            ValueError: invalid shape argument
            sage: ShiftedPrimedTableaux(weight=(2,2,2), shape=[3,2])
            Traceback (most recent call last):
            ...
            ValueError: weight and shape are incompatible
            sage: ShiftedPrimedTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: invalid shape argument
            sage: ShiftedPrimedTableaux(weight=(2,2,2), max_entry=2)
            Traceback (most recent call last):
            ...
            ValueError: maximum entry is incompatible with the weight
            sage: ShiftedPrimedTableaux(shape=[4,1],skew=[3,2])
            Traceback (most recent call last):
            ...
            ValueError: skew shape must be inside the given tableau shape

            sage: SPT1 = ShiftedPrimedTableaux(weight=())
            sage: SPT2 = ShiftedPrimedTableaux(weight=(0,0,0))
            sage: SPT1 is SPT2
            True
        """
        if skew is not None:
            try:
                skew = Partition(skew)
            except ValueError:
                raise ValueError('invalid skew argument')
            if not all(skew[i] > skew[i+1] for i in range(len(skew)-1)):
                raise ValueError('skew shape must be a strict partition')

        if weight is not None:
            weight = tuple(weight)

        if shape is not None:
            if isinstance(shape, SkewPartition):
                skew = shape.inner()
                shape = shape.outer()
            try:
                shape = Partition(shape)
            except (ValueError, TypeError):
                raise ValueError('invalid shape argument')

            if not all(shape[i] > shape[i+1] for i in range(len(shape)-1)):
                raise ValueError("shape {} is not a strict partition".format(shape))

            if (skew is not None and not all(skew[i] <= shape[i]
                                             for i in range(len(skew)))):
                raise ValueError('skew shape must be inside the given tableau shape')

        if weight is not None:
            while weight and weight[-1] == 0:
                weight = weight[:-1]

        if max_entry is not None and weight is not None:
            if len(weight) > max_entry:
                raise ValueError("maximum entry is incompatible with the weight")

        if shape is None:
            if weight is None:
                if max_entry is not None:
                    raise ValueError("specify shape or weight argument")
                return ShiftedPrimedTableaux_all(skew=skew, primed_diagonal=primed_diagonal)
            else:
                return ShiftedPrimedTableaux_weight(weight, skew=skew, primed_diagonal=primed_diagonal)
        else:
            if weight is None:
                return ShiftedPrimedTableaux_shape(shape, max_entry=max_entry, skew=skew, primed_diagonal=primed_diagonal)

            if (skew is not None and sum(shape) - sum(skew) != sum(weight)
                    or skew is None and sum(shape) != sum(weight)):
                raise ValueError("weight and shape are incompatible")

            return ShiftedPrimedTableaux_weight_shape(weight, shape, skew=skew, primed_diagonal=primed_diagonal)

    def __init__(self, skew=None, primed_diagonal=False):
        """
        Initialization of the parent class with given skew shape.

        TESTS::

            sage: SPT = ShiftedPrimedTableaux(skew=[1])
            sage: TestSuite(SPT).run()  # known bug
        """
        self._skew = skew
        self._primed_diagonal = primed_diagonal

    def _element_constructor_(self, T):
        """
        Construct an object from ``T`` as an element of shifted primed
        tableaux, if possible.

        INPUT:

        - ``T`` -- data which can be interpreted as a primed tableau

        OUTPUT: the corresponding primed tableau object

        EXAMPLES::

            sage: SPT = ShiftedPrimedTableaux()
            sage: Primed_SPT = ShiftedPrimedTableaux(primed_diagonal=True)

            sage: tab = SPT([[1,1,"2p"]]); tab
            [(1, 1, 2')]
            sage: tab.parent() is SPT
            True
            sage: tab = Primed_SPT([["1p",1,"2p"]]); tab
            [(1', 1, 2')]
            sage: tab.parent() is Primed_SPT
            True

            sage: tab = SPT([[1,1,2],[2,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1, 2], [2, 2]] is not an element of Shifted Primed Tableaux
            sage: SPT([[1,"2p","2p"]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, '2p', '2p']] is not an element of Shifted Primed Tableaux

            sage: SPT = ShiftedPrimedTableaux(skew=[1])
            sage: SPT([["2p",2]])
            [(None, 2', 2)]

            sage: SPT = ShiftedPrimedTableaux(weight=(2,1))
            sage: tab = SPT([[1,1,1.5]]); tab
            [(1, 1, 2')]
            sage: tab.parent() is SPT
            True

            sage: SPT = ShiftedPrimedTableaux([3])
            sage: tab = SPT([[1,1,1.5]]); tab
            [(1, 1, 2')]
            sage: tab.parent() is SPT
            True
            sage: SPT([[1,1]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1]] is not an element of Shifted Primed Tableaux
             of shape [3]

            sage: SPT = ShiftedPrimedTableaux([3], weight=(2,1))
            sage: tab = SPT([[1,1,1.5]]); tab
            [(1, 1, 2')]
            sage: tab.parent() is SPT
            True
            sage: SPT([[1,1]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 1]] is not an element of Shifted Primed Tableaux
             of weight (2, 1) and shape [3]
        """
        try:
            return self.element_class(self, T, skew=self._skew)
        except ValueError:
            raise ValueError("{} is not an element of {}".format(T, self))

    def _contains_tableau(self, T):
        """
        Check if ``self`` contains preprocessed tableau ``T``.

        TESTS::

            sage: Tabs = ShiftedPrimedTableaux()
            sage: tab = ShiftedPrimedTableau._preprocess(
            ....: [[1,"2p","3p","3p"]])
            sage: tab
            [(1, 2', 3', 3')]
            sage: Tabs._contains_tableau(tab)
            False
            sage: tab = ShiftedPrimedTableau._preprocess(
            ....: [["1p","2p","3p",3]])
            sage: Tabs._contains_tableau(tab)
            False
            sage: Tabs = ShiftedPrimedTableaux(primed_diagonal=True)
            sage: Tabs._contains_tableau(tab)
            True
            sage: Tabs = ShiftedPrimedTableaux(skew=[1])
            sage: tab = ShiftedPrimedTableau._preprocess(
            ....: [["2p","3p",3]], skew=[1])
            sage: tab
            [(None, 2', 3', 3)]
            sage: Tabs._contains_tableau(tab)
            True
        """
        if not all(len(T[i]) > len(T[i+1]) for i in range(len(T)-1)):
            return False
        if self._skew is not None:
            skew = self._skew + [0]*(len(T)-len(self._skew))
        else:
            skew = [0] * len(T)
        for i, row in enumerate(T):
            if i > 0:
                if not all(val > T[i-1][j+1]
                           for j, val in enumerate(row)
                           if j+1 >= skew[i-1] and val.is_unprimed()):
                    return False
                if not all(val >= T[i-1][j+1]
                           for j, val in enumerate(row)
                           if j+1 >= skew[i-1] and val.is_primed()):
                    return False
            if not all(row[j] <= row[j+1]
                       for j in range(skew[i], len(row)-1)
                       if row[j].is_unprimed()):
                return False
            if not all(row[j] < row[j+1]
                       for j in range(skew[i], len(row)-1)
                       if row[j].is_primed()):
                return False
        if not (self._primed_diagonal or all(row[0].is_unprimed()
                   for i, row in enumerate(T)
                   if skew[i] == 0)):
            return False
        return True


class ShiftedPrimedTableaux_all(ShiftedPrimedTableaux):
    """
    The class of all shifted primed tableaux.
    """

    def __init__(self, skew=None, primed_diagonal=False):
        """
        Initialize the class of all shifted tableaux.

        TESTS::

            sage: SPT = ShiftedPrimedTableaux()
            sage: [[1,1.5],[2]] in SPT
            True
            sage: [[1,1.5],[1.5]] in SPT
            False
            sage: [[1,1],[1]] in SPT
            False
            sage: [[1,1],[2,2]] in SPT
            False
            sage: TestSuite(SPT).run()  # long time

            sage: Primed_SPT = ShiftedPrimedTableaux(primed_diagonal=True)
            sage: [[0.5,1.5],[2]] in Primed_SPT
            True
            sage: [[0.5,1.5],[1.5]] in Primed_SPT
            True
            sage: [[0.5,1],[1]] in Primed_SPT
            False
            sage: TestSuite(Primed_SPT).run()  # long time
        """
        if skew is None:
            Parent.__init__(self, category=InfiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=Sets().Infinite())
        ShiftedPrimedTableaux.__init__(self, skew=skew, primed_diagonal=primed_diagonal)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedPrimedTableaux()
            Shifted Primed Tableaux
        """
        if self._skew is None:
            return "Shifted Primed Tableaux"
        return "Shifted Primed Tableaux skewed by {}".format(self._skew)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = ShiftedPrimedTableaux()
            sage: Tabs[:5]
            [[], [(1,)], [(2,)], [(1, 2)], [(1, 2')]]
            sage: Tabs = ShiftedPrimedTableaux(primed_diagonal=True)
            sage: Tabs[:5]
            [[], [(1,)], [(1',)], [(2,)], [(2',)]]
        """
        if self._skew is not None:
            raise NotImplementedError('skew tableau must be empty')
        yield self.element_class(self, [])

        max_entry = 1
        while True:
            for size in range(1, max_entry+1):
                for shape in Partitions(size, max_slope=-1):
                    for weight in OrderedPartitions(size+max_entry-1,
                                                    k=max_entry):
                        weight = [weight[i]-1 for i in range(max_entry)]
                        weight[-1] += 1
                        for tab in ShiftedPrimedTableaux(shape=shape,
                                                         weight=weight,
                                                         primed_diagonal=self._primed_diagonal):
                            yield self.element_class(self, tab, check=False,
                                                     preprocessed=True)
            max_entry += 1


class ShiftedPrimedTableaux_shape(ShiftedPrimedTableaux):
    r"""
    Shifted primed tableaux of a fixed shape.

    Shifted primed tableaux admit a type `A_n` classical crystal structure
    with highest weights corresponding to a given shape.

    The list of module generators consists of all elements of the
    crystal with nonincreasing weight entries.

    The crystal is constructed following operations described in [HPS2017]_
    and [AO2018]_.

    The optional ``primed_diagonal`` allows primed entries in the main diagonal
    of all the Shifted primed tableaux of a fixed shape. If the ``max_entry``
    is ``None`` then ``max_entry`` is set to the total number of entries in the
    tableau if ``primed_diagonal`` is ``True``.

    EXAMPLES::

        sage: ShiftedPrimedTableaux([4,3,1], max_entry=4)
        Shifted Primed Tableaux of shape [4, 3, 1] and maximum entry 4
        sage: ShiftedPrimedTableaux([4,3,1], max_entry=4).cardinality()
        384

    We compute some of the crystal structure::

        sage: SPTC = crystals.ShiftedPrimedTableaux([3,2], 3)
        sage: T = SPTC.module_generators[-1]
        sage: T
        [(1, 1, 2'), (2, 3')]
        sage: T.f(2)
        [(1, 1, 3'), (2, 3')]
        sage: len(SPTC.module_generators)
        7
        sage: SPTC[0]
        [(1, 1, 1), (2, 2)]
        sage: SPTC.cardinality()
        24

    We compare this implementation with the `q(n)`-crystal
    on (tensor products) of letters::

        sage: tableau_crystal = crystals.ShiftedPrimedTableaux([4,1], 3)
        sage: tableau_digraph = tableau_crystal.digraph()
        sage: c = crystals.Letters(['Q', 3])
        sage: tensor_crystal = tensor([c]*5)
        sage: u = tensor_crystal(c(1), c(1), c(1), c(2), c(1))
        sage: subcrystal = tensor_crystal.subcrystal(generators=[u],
        ....:                                        index_set=[1,2,-1])
        sage: tensor_digraph = subcrystal.digraph()
        sage: tensor_digraph.is_isomorphic(tableau_digraph, edge_labels=True)
        True

    If we allow primed entries in the main diagonal::

        sage: ShiftedPrimedTableaux([4,3,1], max_entry=4,
        ....:                        primed_diagonal=True)
        Shifted Primed Tableaux of shape [4, 3, 1] and maximum entry 4
        sage: ShiftedPrimedTableaux([4,3,1], max_entry=4,
        ....:                       primed_diagonal=True).cardinality()
        3072
        sage: SPTC = ShiftedPrimedTableaux([3,2], max_entry=3,
        ....:                              primed_diagonal=True)
        sage: T = SPTC[-1]
        sage: T
        [(1', 2', 2), (3', 3)]
        sage: SPTC[0]
        [(1, 1, 1), (2, 2)]
        sage: SPTC.cardinality()
        96
    """
    @staticmethod
    def __classcall_private__(cls, shape, max_entry=None, skew=None, primed_diagonal=False):
        """
        Normalize the attributes for the class.

        TESTS::

            sage: SPT = ShiftedPrimedTableaux(shape=[2,1])
            sage: SPT._shape.parent()
            Partitions

            sage: SPT1 = ShiftedPrimedTableaux(shape=(2,1), max_entry=3)
            sage: SPT2 = ShiftedPrimedTableaux(shape=[2,1], max_entry=3)
            sage: SPT1 is SPT2
            True
        """
        shape = _Partitions(shape)
        return super().__classcall__(cls,
                     shape=shape, max_entry=max_entry, skew=skew, primed_diagonal=primed_diagonal)

    def __init__(self, shape, max_entry=None, skew=None, primed_diagonal=False):
        """
        Initialize the class of shifted primed tableaux of a given shape.

        If ``max_elt`` is specified, a finite set with entries smaller
        or equal to ``max_elt``.

        TESTS::

            sage: SPT = ShiftedPrimedTableaux([4,2,1], max_entry=4)
            sage: TestSuite(SPT).run()  # long time
            sage: SPT.cartan_type()
            ['Q', 4]

            sage: SPT = ShiftedPrimedTableaux([4,2,1])
            sage: SPT.cartan_type()
            ['A', NN]
            sage: SPT = ShiftedPrimedTableaux([4,2,1], max_entry=4,
            ....:                             primed_diagonal=True)
            sage: TestSuite(SPT).run()  # long time
        """
        if not primed_diagonal:
            # Determine the correct category
            if max_entry is None:
                if skew is None:
                    category = RegularCrystals().Infinite()
                    self._cartan_type = CartanType(['A+oo'])
                    self.Element = CrystalElementShiftedPrimedTableau
                else:
                    category = Sets().Infinite()
            else:
                if skew is None:
                    category = RegularSuperCrystals()
                    self._cartan_type = CartanType(['Q', max_entry])
                    self.Element = CrystalElementShiftedPrimedTableau
                else:
                    category = Sets().Finite()
        else:
            if max_entry is None:
                max_entry = sum(shape)

            if skew is None:
                category = FiniteEnumeratedSets()
            else:
                category = Sets().Finite()

        ShiftedPrimedTableaux.__init__(self, skew=skew, primed_diagonal=primed_diagonal)
        Parent.__init__(self, category=category)
        self._max_entry = max_entry
        if skew is None:
            self._shape = Partition(shape)
        else:
            self._shape = SkewPartition((shape, skew))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedPrimedTableaux([3,2,1])
            Shifted Primed Tableaux of shape [3, 2, 1]
            sage: ShiftedPrimedTableaux([3,2,1], primed_diagonal=True)
            Shifted Primed Tableaux of shape [3, 2, 1] and maximum
            entry 6
        """
        base = "Shifted Primed Tableaux of shape " + self._shape._repr_()
        if self._max_entry is not None:
            base += " and maximum entry {}".format(self._max_entry)
        return base

    def _contains_tableau(self, T):
        """
        TESTS::

           sage: t = ShiftedPrimedTableau._preprocess([[1,'2p',2,2],[2,'3p']])
           sage: ShiftedPrimedTableaux([4,2],max_entry=4)._contains_tableau(t)
           True
           sage: s = ShiftedPrimedTableau._preprocess([[1,'2p',2],[2,'3p']])
           sage: ShiftedPrimedTableaux([4,2])._contains_tableau(s)
           False
           sage: t = ShiftedPrimedTableau._preprocess([['1p','2p',2,2],
           ....:                                       ['2p','3p']])
           sage: ShiftedPrimedTableaux([4,2],max_entry=4,
           ....:                      primed_diagonal=True)._contains_tableau(t)
           True
        """
        if not super()._contains_tableau(T):
            return False

        shape = [len(row) for row in T]
        skew = [row.count(None) for row in T]
        if sum(skew) == 0:
            shape = Partition(shape)
        else:
            shape = SkewPartition((shape, skew))
        if self._shape != shape:
            return False

        if self._max_entry is not None:
            flat = [item.integer() for sublist in T for item in sublist]
            if flat == []:
                max_entry = 0
            else:
                max_entry = int(max(flat))
            if max_entry > self._max_entry:
                return False

        return True

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = ShiftedPrimedTableaux(shape=(3,2), max_entry=5)
            sage: Tabs[:4]
            [[(1, 1, 1), (2, 2)],
             [(1, 1, 1), (2, 3)],
             [(1, 1, 1), (2, 3')],
             [(1, 1, 2), (2, 3)]]
            sage: len(list(Tabs))
            376
            sage: Tabs = ShiftedPrimedTableaux(shape=(3,2), primed_diagonal=True)
            sage: Tabs[:4]
            [[(1, 1, 1), (2, 2)],
             [(1, 1, 1), (2', 2)],
             [(1', 1, 1), (2, 2)],
             [(1', 1, 1), (2', 2)]]
            sage: len(list(Tabs))
            1504
        """
        if not self._primed_diagonal:
            for T in super().__iter__():
                yield T
            return

        from sage.combinat.permutation import Permutations
        list_weights = []
        for partition in Partitions(sum(self._shape)):
            if len(partition) <= self._max_entry:
                for c in Combinations(range(self._max_entry), len(partition)):
                    for p in Permutations(partition):
                        weight = [0] * self._max_entry
                        for i,val in enumerate(p):
                            weight[c[i]] = val
                        list_weights.append(weight)
        for weight in list_weights:
            for T in ShiftedPrimedTableaux(weight=tuple(weight), shape=self._shape,
                                           primed_diagonal=self._primed_diagonal):
                yield self.element_class(self, T, preprocessed=True, check=False)

    @lazy_attribute
    def module_generators(self):
        """
        Return the generators of ``self`` as a crystal.

        TESTS::

            sage: SPT = ShiftedPrimedTableaux(shape=[2,1])
            sage: SPT.module_generators
            ([(1, 1), (2,)], [(1, 2), (3,)], [(1, 2'), (3,)])
        """
        if self._skew is not None:
            raise NotImplementedError("only for non-skew shapes")
        list_dw = []
        if self._max_entry is None:
            max_entry = sum(self._shape)
        else:
            max_entry = self._max_entry
        for weight in (Partition(self._shape).dominated_partitions(rows=max_entry)):
            list_dw.extend([self.element_class(self, T, check=False,
                                               preprocessed=True)
                            for T in ShiftedPrimedTableaux(weight=tuple(weight),
                                                           shape=self._shape,
                                                           primed_diagonal=self._primed_diagonal)])
        return tuple(list_dw)

    def shape(self):
        """
        Return the shape of the shifted tableaux ``self``.

        TESTS::

            sage: ShiftedPrimedTableaux([6,4,3,1]).shape()
            [6, 4, 3, 1]
        """
        return self._shape


class ShiftedPrimedTableaux_weight(ShiftedPrimedTableaux):
    """
    Shifted primed tableaux of fixed weight.

    EXAMPLES::

        sage: ShiftedPrimedTableaux(weight=(2,3,1))
        Shifted Primed Tableaux of weight (2, 3, 1)
        sage: ShiftedPrimedTableaux(weight=(2,3,1)).cardinality()
        17
        sage: SPT = ShiftedPrimedTableaux(weight=(2,3,1), primed_diagonal=True)
        sage: SPT.cardinality()
        64
        sage: T = ShiftedPrimedTableaux(weight=(3,2), primed_diagonal=True)
        sage: T[:5]
        [[(1, 1, 1, 2, 2)],
         [(1, 1, 1, 2', 2)],
         [(1', 1, 1, 2, 2)],
         [(1', 1, 1, 2', 2)],
         [(1, 1, 1, 2), (2,)]]
        sage: T.cardinality()
        16
    """

    def __init__(self, weight, skew=None, primed_diagonal=False):
        """
        Initialize the class of shifted primed tableaux of a given weight.

        TESTS::

            sage: TestSuite( ShiftedPrimedTableaux(weight=(3,2,1)) ).run()
            sage: TestSuite( ShiftedPrimedTableaux(weight=(3,2,1),
            ....:                               primed_diagonal=True) ).run()
        """
        ShiftedPrimedTableaux.__init__(self, skew=skew, primed_diagonal=primed_diagonal)
        if skew is None:
            Parent.__init__(self, category=FiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=Sets().Finite())
        self._weight = weight

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedPrimedTableaux(weight=(3,2,1))
            Shifted Primed Tableaux of weight (3, 2, 1)
        """
        if self._skew is None:
            return "Shifted Primed Tableaux of weight {}".format(self._weight)
        return "Shifted Primed Tableaux of weight {} skewed by {}".format(self._weight, self._skew)

    def _contains_tableau(self, T):
        """
        Check if ``self`` contains preprocessed tableau ``T``.

        TESTS::

            sage: t = ShiftedPrimedTableau._preprocess([[1,1.5],[2]])
            sage: ShiftedPrimedTableaux(weight=(1,2))._contains_tableau(t)
            True
            sage: t = ShiftedPrimedTableau._preprocess([[0.5,1.5],[2]])
            sage: ShiftedPrimedTableaux(weight=(1,2),
            ....:                 primed_diagonal=True)._contains_tableau(t)
            True
            sage: s = ShiftedPrimedTableau._preprocess([[1,1.5],[3]])
            sage: ShiftedPrimedTableaux(weight=(1,2))._contains_tableau(s)
            False

            sage: u = ShiftedPrimedTableau._preprocess([])
            sage: ShiftedPrimedTableaux(weight=())._contains_tableau(u)
            True
            sage: ShiftedPrimedTableaux(weight=(1,2))._contains_tableau(u)
            False
        """
        if not super()._contains_tableau(T):
            return False

        flat = [item.integer() for sublist in T for item in sublist]
        if not flat:
            return not self._weight
        max_ind = max(flat)
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        return self._weight == weight

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = ShiftedPrimedTableaux(weight=(2,3))
            sage: Tabs[:4]
            [[(1, 1, 2, 2, 2)],
             [(1, 1, 2', 2, 2)],
             [(1, 1, 2, 2), (2,)],
             [(1, 1, 2', 2), (2,)]]
            sage: len(list(Tabs))
            5
            sage: Tabs = ShiftedPrimedTableaux(weight=(2,3),
            ....:                              primed_diagonal=True)
            sage: Tabs[:4]
            [[(1, 1, 2, 2, 2)],
             [(1, 1, 2', 2, 2)],
             [(1', 1, 2, 2, 2)],
             [(1', 1, 2', 2, 2)]]
            sage: len(list(Tabs))
            16
        """
        for shape_ in ZS1_iterator(sum(self._weight)):
            if all(shape_[i] > shape_[i+1] for i in range(len(shape_)-1)):
                for tab in ShiftedPrimedTableaux(shape=shape_, weight=self._weight,
                                                 skew=self._skew,
                                                 primed_diagonal=self._primed_diagonal
                                                 ):
                    yield self.element_class(self, tab, check=False, preprocessed=True)


class ShiftedPrimedTableaux_weight_shape(ShiftedPrimedTableaux):
    """
    Shifted primed tableaux of the fixed weight and shape.

    EXAMPLES::

        sage: ShiftedPrimedTableaux([4,2,1], weight=(2,3,2))
        Shifted Primed Tableaux of weight (2, 3, 2) and shape [4, 2, 1]
        sage: ShiftedPrimedTableaux([4,2,1], weight=(2,3,2)).cardinality()
        4
        sage: T = ShiftedPrimedTableaux([4,2,1], weight=(2,3,2),
        ....:                           primed_diagonal=True)
        sage: T[:6]
        [[(1, 1, 2, 2), (2, 3'), (3,)],
         [(1, 1, 2, 2), (2, 3'), (3',)],
         [(1, 1, 2, 2), (2', 3'), (3,)],
         [(1, 1, 2, 2), (2', 3'), (3',)],
         [(1, 1, 2', 3), (2, 2), (3,)],
         [(1, 1, 2', 3), (2, 2), (3',)]]
        sage: T.cardinality()
        32
    """

    def __init__(self, weight, shape, skew=None, primed_diagonal=False):
        """
        Initialize the class of shifted primed tableaux of the given weight
        and shape.

        TESTS::

            sage: TestSuite( ShiftedPrimedTableaux([4,2,1],
            ....:                                  weight=(3,2,2)) ).run()
            sage: TestSuite( ShiftedPrimedTableaux([4,2,1],
            ....:                 weight=(3,2,2), primed_diagonal=True) ).run()
        """
        ShiftedPrimedTableaux.__init__(self, skew=skew, primed_diagonal=primed_diagonal)
        if skew is None:
            Parent.__init__(self, category=FiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=Sets().Finite())
        self._weight = weight
        if skew is None:
            self._shape = _Partitions(shape)
        else:
            self._shape = SkewPartition((shape, skew))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedPrimedTableaux([3,2,1], weight=(4,2))
            Shifted Primed Tableaux of weight (4, 2) and shape [3, 2, 1]
        """
        return ("Shifted Primed Tableaux of weight {} and shape {}"
                .format(self._weight, self._shape))

    def _contains_tableau(self, T):
        """
        Check if ``self`` contains preprocessed tableau ``T``.

        TESTS::

            sage: t = ShiftedPrimedTableau._preprocess([[1,1.5],[2]])
            sage: ShiftedPrimedTableaux([2,1],
            ....:                       weight=(1,2))._contains_tableau(t)
            True
            sage: t = ShiftedPrimedTableau._preprocess([[0.5,1.5],[2]])
            sage: ShiftedPrimedTableaux([2,1],
            ....:                       weight=(1,2))._contains_tableau(t)
            False
            sage: t = ShiftedPrimedTableau._preprocess([[0.5,1.5],[2]])
            sage: ShiftedPrimedTableaux([2,1], weight=(1,2),
            ....:                     primed_diagonal=True)._contains_tableau(t)
            True
            sage: ShiftedPrimedTableaux([2,1],
            ....:                           weight=(2,1))._contains_tableau(t)
            False
            sage: s = ShiftedPrimedTableau._preprocess([[1,1.5,2,3],[3]])
            sage: ShiftedPrimedTableaux([3,2],
            ....:                         weight=(1,2,2))._contains_tableau(s)
            False

            sage: u = ShiftedPrimedTableau._preprocess([])
            sage: ShiftedPrimedTableaux([3,2],
            ....:                       weight=(1,2,2))._contains_tableau(u)
            False
            sage: ShiftedPrimedTableaux([], weight=())._contains_tableau(u)
            True
        """
        if not super()._contains_tableau(T):
            return False

        flat = [item.integer() for sublist in T for item in sublist]
        if not flat:
            # It is sufficient only to check this because the weight
            #   and shape must be compatible
            return not self._weight

        max_ind = max(flat)
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        if self._weight != weight:
            return False

        shape = [len(row) for row in T]
        skew = [row.count(None) for row in T]
        if sum(skew) == 0:
            shape = _Partitions(shape)
        else:
            shape = SkewPartition((shape, skew))
        return self._shape == shape

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: Tabs = ShiftedPrimedTableaux([3,2], weight=(1,2,2))
            sage: Tabs[:4]
            [[(1, 2, 2), (3, 3)],
             [(1, 2', 3'), (2, 3)],
             [(1, 2', 3'), (2, 3')],
             [(1, 2', 2), (3, 3)]]
            sage: len(list(Tabs))
            4
            sage: Tabs = ShiftedPrimedTableaux([3,2],weight=(1,2,2),
            ....:                              primed_diagonal=True)
            sage: Tabs[:4]
            [[(1, 2, 2), (3, 3)],
             [(1, 2, 2), (3', 3)],
             [(1, 2', 3'), (2, 3)],
             [(1, 2', 3'), (2, 3')]]
            sage: len(list(Tabs))
            16

        TESTS::

            sage: Tabs = ShiftedPrimedTableaux([3,2], weight=(1,4))
            sage: list(Tabs)
            []
        """
        if self._skew is not None:
            raise NotImplementedError('skew tableau must be empty')

        if not self._shape.dominates(sorted(self._weight, reverse=True)):
            return
        full_shape = self._shape
        sub_tab = []
        tab_list_new = [[]]
        half = ~QQ(2)
        for i, w in enumerate(self._weight):
            tab_list_old = tab_list_new
            tab_list_new = []
            for sub_tab in tab_list_old:
                sub_shape = [len(row) for row in sub_tab]
                for strip in _add_strip(sub_shape, full_shape, w):
                    l = len(strip) // 2
                    new_tab = []
                    new_tab1 = None
                    if len(sub_shape) < len(full_shape):
                        new_tab = [sub_tab[r] + [i+half]*strip[r] + [i+1]*strip[-r-1]
                                   for r in range(l-1)]
                        if strip[l] != 0:
                            if self._primed_diagonal:
                                new_tab1 = new_tab[:]
                                new_tab1.append([i+half] + [i+1] * (strip[l]-1))
                            new_tab.append([i+1] * strip[l])
                    else:
                        new_tab = [sub_tab[r] + [i+half]*strip[r] + [i+1]*strip[-r-1]
                                   for r in range(l)]
                    tab_list_new.append(new_tab)
                    if new_tab1:
                        tab_list_new.append(new_tab1)
        for tab in tab_list_new:
            yield self.element_class(self, tab)


####################
# Helper functions #
####################


def _add_strip(sub_tab, full_tab, length):
    """
    Helper function used in the algorithm to generate all shifted primed
    tableaux of the fixed weight and shape.

    TESTS::

        sage: list(ShiftedPrimedTableaux([3,1], weight=(2,2)))  # indirect doctest
        [[(1, 1, 2), (2,)], [(1, 1, 2'), (2,)]]
        sage: list(ShiftedPrimedTableaux([3,1], weight=(2,2),   # indirect doctest
        ....:                            primed_diagonal=True))
        [[(1, 1, 2), (2,)],
         [(1, 1, 2), (2',)],
         [(1, 1, 2'), (2,)],
         [(1, 1, 2'), (2',)],
         [(1', 1, 2), (2,)],
         [(1', 1, 2), (2',)],
         [(1', 1, 2'), (2,)],
         [(1', 1, 2'), (2',)]]
    """
    if sum(sub_tab) + length > sum(full_tab):
        raise ValueError("strip does not fit")

    if not sub_tab:
        cliff_list = []
    else:
        cliff_list = [int(sub_tab[0] != full_tab[0])]

    for row in range(1, len(sub_tab)):
        if sub_tab[row] == full_tab[row]:
            cliff_list.append(0)
        elif sub_tab[row-1] - 1 == sub_tab[row]:
            cliff_list[-1] += 1
        else:
            cliff_list.append(1)

    if len(sub_tab) < len(full_tab):
        cliff_list.append(0)

    for primes_num in range(min(sum(cliff_list), length) + 1):
        for primed_list in IntegerVectors(n=primes_num, k=len(cliff_list),
                                          outer=cliff_list):
            row = 0
            primed_strip = []
            for i, cliff in enumerate(cliff_list):
                if cliff == 0:
                    row += 1
                    primed_strip.append(0)
                primed_strip.extend([int(primed_list[i] > j)
                                     for j in range(cliff)])
                row += cliff
            plat_list = []

            if sub_tab and len(sub_tab) < len(full_tab):
                plat_list.append(min(sub_tab[-1] + primed_strip[-2] - 1,
                                     full_tab[len(sub_tab)]))
            for row in reversed(range(1, len(sub_tab))):
                plat_list.append(min(sub_tab[row-1]+primed_strip[row-1]-1, full_tab[row])
                                 - sub_tab[row] - primed_strip[row])
            if sub_tab:
                plat_list.append(full_tab[0] - sub_tab[0] - primed_strip[0])
            else:
                plat_list.append(full_tab[0])

            for non_primed_strip in IntegerVectors(n=length-primes_num,
                                                   k=len(plat_list),
                                                   outer=plat_list):
                yield list(primed_strip) + list(non_primed_strip)
