r"""
Quasi-ribbon tableaux

AUTHORS:

- Pranita Urlam (2026): initial version

A *quasi-ribbon tableau* is a filling of a ribbon-shaped skew Young diagram
with positive integers, weakly increasing along rows (left to right) and
satisfying a strict column condition at the joints between consecutive rows.
These tableaux arise as the insertion tableaux in the Krob-Thibon (KT) RSK
correspondence for the hypoplactic monoid.

Internally a quasi-ribbon tableau is stored as a list of weakly increasing
rows ordered **bottom to top**. Consecutive rows must satisfy the
*join condition*: the last entry of the lower row is strictly less than the
first entry of the row above.

REFERENCES:

- [KT1997]_ D. Krob and J.-Y. Thibon, *Noncommutative symmetric functions IV:
  Quantum linear groups and Hecke algebras at q = 0*, J. Algebraic Combin. **6**
  (1997), 339--376.

- [Nov2000]_ J.-C. Novelli, *On the hypoplactic monoid*, Discrete Math. **217**
  (2000), 315--336.

EXAMPLES::

    sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
    sage: QuasiRibbonTableau([])
    []
    sage: QuasiRibbonTableau([[1, 2], [3, 4]])
    [[1, 2], [3, 4]]
    sage: QuasiRibbonTableau([[1, 3], [4, 4, 6]]).shape()
    [2, 3]
    sage: QuasiRibbonTableau([[1, 2], [3, 4]]).reading_word()
    [1, 2, 3, 4]
"""
# ****************************************************************************
#       Copyright (C) 2026 Pranita Urlam
#
#  Distributed under the terms of the GNU General Public License (GPL)
#              https://www.gnu.org/licenses/
# ****************************************************************************

from bisect import bisect_right

from sage.combinat.combinat import CombinatorialElement
from sage.combinat.composition import Composition
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.categories.sets_cat import Sets
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class QuasiRibbonTableau(CombinatorialElement, metaclass=ClasscallMetaclass):
    r"""
    A semistandard quasi-ribbon tableau.

    A *quasi-ribbon tableau* is a list of weakly increasing rows
    `(R_1, R_2, \ldots, R_k)` (indexed bottom to top) of positive integers
    satisfying the *join condition*:

    .. MATH::

        \mathrm{last}(R_i) < \mathrm{first}(R_{i+1})
        \quad \text{for all } 1 \le i < k.

    Geometrically, this is a semistandard filling of a ribbon-shaped skew
    Young diagram: rows are weakly increasing and the two cells that share a
    column at each row-junction satisfy a strict inequality.

    INPUT:

    - ``rows`` -- list of lists of positive integers, ordered bottom to top

    EXAMPLES::

        sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
        sage: QuasiRibbonTableau([])
        []
        sage: QuasiRibbonTableau([[3]])
        [[3]]
        sage: QuasiRibbonTableau([[1, 2], [3, 3, 4]])
        [[1, 2], [3, 3, 4]]

    Two quasi-ribbon tableaux are equal iff they have identical rows::

        sage: QuasiRibbonTableau([[1, 2], [3]]) == QuasiRibbonTableau([[1, 2], [3]])
        True
        sage: QuasiRibbonTableau([[1, 2], [3]]) == QuasiRibbonTableau([[1, 3], [4]])
        False

    TESTS::

        sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
        sage: QuasiRibbonTableau([[2, 1]])
        Traceback (most recent call last):
        ...
        ValueError: each row must be weakly increasing
        sage: QuasiRibbonTableau([[1, 3], [2]])
        Traceback (most recent call last):
        ...
        ValueError: join condition violated: last entry of a row must be strictly less than the first entry of the row above
    """

    @staticmethod
    def __classcall_private__(cls, rows):
        r"""
        Normalise input and route to the parent's element constructor.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: T = QuasiRibbonTableau([[1, 2], [3]])
            sage: T.parent()
            Quasi-ribbon tableaux
        """
        rows = [list(r) for r in rows]
        return QuasiRibbonTableaux_all().element_class(QuasiRibbonTableaux_all(), rows)

    def __init__(self, parent, rows):
        r"""
        Initialise ``self``.

        TESTS::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: T = QuasiRibbonTableau([[1, 2], [3]])
            sage: TestSuite(T).run()
        """
        if not all(isinstance(row, list) for row in rows):
            raise ValueError("rows must be a list of lists")
        for row in rows:
            if any(row[i] > row[i + 1] for i in range(len(row) - 1)):
                raise ValueError("each row must be weakly increasing")
        for i in range(len(rows) - 1):
            if rows[i][-1] >= rows[i + 1][0]:
                raise ValueError(
                    "join condition violated: last entry of a row must be "
                    "strictly less than the first entry of the row above"
                )
        CombinatorialElement.__init__(self, parent, rows)

    def _repr_(self):
        r"""
        Return a string representation.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: QuasiRibbonTableau([[1, 2], [3, 3]])
            [[1, 2], [3, 3]]
        """
        return repr(self._list)

    def rows(self):
        r"""
        Return the list of rows ordered bottom to top.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: QuasiRibbonTableau([[1, 2], [3]]).rows()
            [[1, 2], [3]]
        """
        return list(self._list)

    def shape(self):
        r"""
        Return the shape of ``self`` as a :class:`~sage.combinat.composition.Composition`.

        The `i`-th part of the shape is the length of the `i`-th row (bottom to top).

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: QuasiRibbonTableau([]).shape()
            []
            sage: QuasiRibbonTableau([[1, 2], [3, 3, 4]]).shape()
            [2, 3]
            sage: QuasiRibbonTableau([[5], [6], [7]]).shape()
            [1, 1, 1]
        """
        if not self._list:
            return Composition([])
        return Composition([len(row) for row in self._list])

    def size(self):
        r"""
        Return the number of cells (total number of entries).

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: QuasiRibbonTableau([]).size()
            0
            sage: QuasiRibbonTableau([[1, 2], [3, 3, 4]]).size()
            5
        """
        return sum(len(row) for row in self._list)

    def reading_word(self):
        r"""
        Return the row-reading word of ``self``.

        The reading word is obtained by reading each row left to right,
        starting from the bottom row and ending at the top row.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: QuasiRibbonTableau([]).reading_word()
            []
            sage: QuasiRibbonTableau([[1, 2], [3, 4]]).reading_word()
            [1, 2, 3, 4]
            sage: QuasiRibbonTableau([[1], [2], [3]]).reading_word()
            [1, 2, 3]
        """
        return [x for row in self._list for x in row]

    def insert(self, x):
        r"""
        Return the quasi-ribbon tableau obtained by KT row-inserting letter ``x``.

        The Krob-Thibon row insertion of ``x`` into ``self`` proceeds as
        follows.  Let `(R_1, R_2, \ldots, R_k)` be the rows of ``self``
        (bottom to top).

        1. If ``self`` is empty, return a one-cell ribbon `[[x]]`.

        2. Find the leftmost position `j` in `R_1` with `R_1[j] > x`
           (equivalently `j = \mathrm{bisect\_right}(R_1, x)`).

           - If `j = 0` (all of `R_1` exceeds `x`): prepend a new bottom
             row `[x]`.

           - If `j = |R_1|` (`x \geq` all of `R_1`):

             - If there is no upper row, or `x <` first entry of `R_2`:
               append `x` to `R_1` (join condition is maintained).

             - Otherwise, leave `R_1` unchanged and recursively insert
               `x` into the upper rows.

           - Otherwise (`0 < j < |R_1|`): replace `R_1[j]` by `x`,
             then recursively insert the bumped value `R_1[j]` into the
             upper rows (or start a new top row if there are none).

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: T = QuasiRibbonTableau([])
            sage: T.insert(3)
            [[3]]
            sage: T.insert(3).insert(1)
            [[1], [3]]
            sage: T.insert(3).insert(1).insert(3)
            [[1], [3, 3]]
            sage: T.insert(3).insert(1).insert(3).insert(2)
            [[1, 2], [3, 3]]

        Inserting an increasing sequence gives a single row::

            sage: T = QuasiRibbonTableau([])
            sage: for x in [1, 2, 3, 4]:
            ....:     T = T.insert(x)
            sage: T
            [[1, 2, 3, 4]]

        Inserting a decreasing sequence gives a single column::

            sage: T = QuasiRibbonTableau([])
            sage: for x in [4, 3, 2, 1]:
            ....:     T = T.insert(x)
            sage: T
            [[1], [2], [3], [4]]
        """
        new_rows = _kt_insert(self._list, x)
        return QuasiRibbonTableau(new_rows)

    def check(self):
        r"""
        Check that ``self`` is a valid quasi-ribbon tableau.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: QuasiRibbonTableau([[1, 2], [3]]).check()

        TESTS::

            sage: T = QuasiRibbonTableau([[1, 2], [3]])
            sage: T.check()
        """
        for row in self._list:
            if any(row[i] > row[i + 1] for i in range(len(row) - 1)):
                raise ValueError("each row must be weakly increasing")
        for i in range(len(self._list) - 1):
            if self._list[i][-1] >= self._list[i + 1][0]:
                raise ValueError("join condition violated")


def _kt_insert(rows, x):
    r"""
    Core KT row-insertion step.

    Insert the letter ``x`` into the list of rows ``rows`` (ordered bottom to
    top), returning a new list of rows.  This is a pure Python helper called
    by :meth:`QuasiRibbonTableau.insert`.

    EXAMPLES::

        sage: from sage.combinat.quasi_ribbon_tableau import _kt_insert
        sage: _kt_insert([], 5)
        [[5]]
        sage: _kt_insert([[1, 3]], 2)
        [[1, 2], [3]]
        sage: _kt_insert([[1, 2], [4]], 3)
        [[1, 2, 3], [4]]
        sage: _kt_insert([[1, 2], [3]], 3)
        [[1, 2], [3, 3]]
    """
    if not rows:
        return [[x]]

    R = list(rows[0])
    upper = rows[1:]

    j = bisect_right(R, x)  # leftmost position where R[j] > x

    if j == 0:
        # x is smaller than every entry in R: start a new bottom row.
        return [[x]] + list(rows)
    if j == len(R):
        # x >= all entries in R.
        if not upper or x < upper[0][0]:
            # Safe to extend R (join condition maintained).
            return [R + [x]] + list(upper)
        # x would violate the join condition with the row above;
        # leave R intact and insert into the upper rows.
        return [R] + _kt_insert(upper, x)
    # 0 < j < len(R): bump R[j] out of R.
    bumped = R[j]
    new_R = R[:j] + [x] + R[j + 1:]
    if upper:
        return [new_R] + _kt_insert(upper, bumped)
    return [new_R, [bumped]]


# ---------------------------------------------------------------------------
# Parent class
# ---------------------------------------------------------------------------

class QuasiRibbonTableaux_all(UniqueRepresentation, Parent):
    r"""
    The parent of all quasi-ribbon tableaux.

    EXAMPLES::

        sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
        sage: QRT = QuasiRibbonTableaux()
        sage: QRT
        Quasi-ribbon tableaux
        sage: QuasiRibbonTableaux()([[1, 2], [3]])
        [[1, 2], [3]]
    """

    def __init__(self):
        r"""
        Initialise ``self``.

        TESTS::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: QRT = QuasiRibbonTableaux()
            sage: TestSuite(QRT).run()
        """
        Parent.__init__(self, category=Sets())

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: QuasiRibbonTableaux()
            Quasi-ribbon tableaux
        """
        return "Quasi-ribbon tableaux"

    def _element_constructor_(self, rows):
        r"""
        Construct a :class:`QuasiRibbonTableau` from ``rows``.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: QuasiRibbonTableaux()([[1, 2], [3]])
            [[1, 2], [3]]
        """
        return self.element_class(self, [list(r) for r in rows])

    Element = QuasiRibbonTableau


def QuasiRibbonTableaux():
    r"""
    Return the parent of all quasi-ribbon tableaux.

    EXAMPLES::

        sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
        sage: QuasiRibbonTableaux()
        Quasi-ribbon tableaux
    """
    return QuasiRibbonTableaux_all()
