r"""
Quasi-ribbon tableaux

AUTHORS:

- Daniel Chen, Lisa Johnston, Junbok Lee, Evuilynn Nguyen, Heather Ross, Anne Schilling, Chenchen Zhao (2026): initial version

This file implements quasi-ribbon tableaux and finite families of
quasi-ribbon tableaux. A quasi-ribbon tableau is entered as a list of rows
from top to bottom. The entries in each row are weakly increasing left to right, the
entries in each column are strictly increasing top to bottom, and ``None`` entries are used
to record the shifted positions in the ribbon shape.

The main functionality includes constructing quasi-ribbon tableaux,
computing their associated compositions, computing column-reading words,
generating finite families of fixed composition shape, and performing
Krob--Thibon insertion of letters and words. For references, see [KT1997]_ and [Nov2000]_.
"""

# ****************************************************************************
#       Copyright (C) 2026 Daniel Chen, Lisa Johnston, Junbok Lee, Evuilynn Nguyen, Heather Ross, Anne Schilling, Chenchen Zhao
#
#  Distributed under the terms of the GNU General Public License (GPL)
#              https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.sets_cat import Sets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.skew_tableau import SkewTableau, SkewTableaux
from sage.combinat.composition import Composition, Compositions
from sage.rings.integer import Integer
from sage.sets.positive_integers import PositiveIntegers
from sage.combinat.words.words import Words
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.functions.other import binomial


class QuasiRibbonTableau(SkewTableau):
    r"""
    A quasi-ribbon tableau.

    A quasi-ribbon shape is obtained
    from a composition by drawing rows of cells whose lengths are the parts of
    the composition, and then shifting the rows so that each lower row starts
    under the last cell of the row above. For the purposes of this class, a quasi-ribbon
    tableau is a filling of a quasi-ribbon shape with positive integers which:

    1) has at least one entry in row `1`;
    2) weakly increases from left to right in each row; and
    3) strictly increases from top to bottom in each column.

    A quasi-ribbon tableau is given by a list of the rows from top to bottom.

    EXAMPLES::

        sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
        sage: Q = QuasiRibbonTableau([[1, 2, 3],
        ....:                         [None, None, 4, 5]])
        sage: Q
        [[1, 2, 3], [None, None, 4, 5]]
        sage: Q.pp()
        1  2  3
              4  5
        sage: Q.to_composition()
        [3, 2]

    The entries labeled by ``None`` correspond to shifted positions before the
    first actual entry in a row.  Using ``None`` is optional;
    the entries will be shifted accordingly::

        sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
        sage: Q = QuasiRibbonTableau([[1, 2, 3], [4, 5]])
        sage: Q.pp()
        1  2  3
              4  5

    TESTS::

        sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
        sage: Q = QuasiRibbonTableau([[1], [2, 3], [4, 4]])
        sage: Q.to_word()
        word: 12344
        sage: QuasiRibbonTableau([[1,2],[3,4]]).evaluation()
        [1, 1, 1, 1]
        sage: hash(tuple(Q)) == hash(Q)
        True
        sage: TestSuite(Q).run()
    """
    @staticmethod
    def __classcall_private__(cls, rows):
        r"""
        Construct a quasi-ribbon tableau from ``rows``.

        Construction is delegated to :class:`QuasiRibbonTableaux`, whose
        element constructor normalizes and validates the input.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: QuasiRibbonTableau([[1, 2], [None, 3]])
            [[1, 2], [None, 3]]
            sage: QuasiRibbonTableau([[1, 2, 3], [4, 5]])
            [[1, 2, 3], [None, None, 4, 5]]
            sage: QuasiRibbonTableau([[1, 2, 3], [None, 4, 5]])
            [[1, 2, 3], [None, None, 4, 5]]
            sage: QuasiRibbonTableau([[1, 2, 3], [4, 5], [6]])
            [[1, 2, 3], [None, None, 4, 5], [None, None, None, 6]]
        """
        return QuasiRibbonTableaux()(rows)

    def __init__(self, parent, rows):
        r"""
        Initialize ``self`` from normalized quasi-ribbon rows.

        The rows are given from top to bottom, contain no ``None`` entries,
        and are assumed to satisfy the quasi-ribbon tableau conditions.
        Validation is performed by
        :meth:`QuasiRibbonTableaux._element_constructor_` during ordinary
        construction.

        This method adds the shifts determined by the quasi-ribbon shape and
        reverses the rows before initializing :class:`SkewTableau`.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([[1, 2, 3],
            ....:                         [None, None, 4, 5]])
            sage: Q.rows()
            [[1, 2, 3], [None, None, 4, 5]]
        """
        self._filling = rows

        shifted_rows = []
        shift = 0
        for row in self._filling:
            shifted_rows.append([None] * shift + list(row))
            shift += len(row) - 1

        SkewTableau.__init__(self, parent, shifted_rows[::-1])

    def rows(self):
        """
        Return the quasi-ribbon rows of ``self``.

        This returns the user-facing row-list representation, including
        ``None`` entries.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([[1, 2, 3],
            ....:                         [None, None, 4, 5]])
            sage: Q.rows()
            [[1, 2, 3], [None, None, 4, 5]]
        """
        return [list(row) for row in self[::-1]]

    def to_composition(self):
        """
        Return the composition shape of ``self``.

        The shape counts the actual entries in each row, ignoring ``None``.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([[1, 2, 3],
            ....:                         [None, None, 4, 5]])
            sage: Q.to_composition()
            [3, 2]
        """
        # Each row length in the composition is the number of non-None entries
        # in that row.
        return Composition(
            [sum(entry is not None for entry in row) for row in self.rows()])

    def width(self):
        """
        Return the maximum row length of ``self``, including ``None`` entries.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([[1, 2, 3],
            ....:                         [None, None, 4, 5]])
            sage: Q.width()
            4
        """
        # The empty tableau has width 0.
        if not self.rows():
            return 0

        # We include None entries because they represent shifted positions.
        return max(len(row) for row in self.rows())

    def to_word_by_column(self):
        """
        Return the reading word of ``self``.

        The reading convention is column by column from left to right.
        Within each column, read from bottom to top.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([[1, 2, 3],
            ....:                         [None, None, 4, 5]])
            sage: Q.to_word_by_column()
            word: 12435

            sage: T = QuasiRibbonTableau([[1], [2, 2], [None, 3, 3], [None, None, 4]])
            sage: T.to_word_by_column()
            word: 213243
        """
        rows = self.rows()
        word = []

        for col in range(self.width()):
            for row in rows[::-1]:
                if col < len(row) and row[col] is not None:
                    word.append(row[col])

        return Words(PositiveIntegers())(word)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        This displays the user-facing quasiribbon rows, not the internally
        flipped skew-tableau rows.

        EXAMPLES::
            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([[1,2,2],[3]])
            sage: Q
            [[1, 2, 2], [None, None, 3]]
        """
        # Keep the raw row-list output for debugging because it shows None
        # entries explicitly.
        return repr(self.rows())

    def pp(self):
        """
        Return a pretty print of ``self``.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([[1, 2, 3],
            ....:                         [None, None, 4, 5]])
            sage: Q.pp()
            1  2  3
                  4  5
        """
        lines = []
        for row in self.rows():
            # Print None entries as blank spaces rather than the word "None".
            line = "  ".join(" " if entry is None else str(entry) for entry in row)
            lines.append(line)
        print("\n".join(lines))

    def _entries_with_positions(self):
        """
        Return all entries of ``self`` with their positions in the cleaned rows.

        Each item is a tuple (entry, row_index, col_index).

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([[1, 2, 3], [4, 5]])
            sage: Q._entries_with_positions()
            [(1, 0, 0), (2, 0, 1), (3, 0, 2), (4, 1, 0), (5, 1, 1)]
        """
        entries = []
        for row_index, row in enumerate(self._filling):
            for col_index, entry in enumerate(row):
                if entry is not None:
                    entries.append((entry, row_index, col_index))
        return entries

    def _rightmost_bottommost_leq(self, a):
        """
        Return rightmost, bottommost entry of ``self`` less than or equal to ``a``
        as ``(entry, row_index, col_index)`` or ``None`` if no such entry exists.

        We compare positions by column first, then row. So "rightmost" means
        largest column index. If there is a tie, "bottommost" means largest
        row index.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([[1, 2, 3], [4, 5]])
            sage: Q._rightmost_bottommost_leq(3)
            (3, 0, 2)
            sage: Q = QuasiRibbonTableau([[1, 2, 2], [3, 3]])
            sage: Q._rightmost_bottommost_leq(3)
            (3, 1, 1)
        """
        candidates = []
        for entry, row_index, col_index in self._entries_with_positions():
            if entry <= a:
                candidates.append((entry, row_index, col_index))
        if not candidates:
            return None
        # Pick the candidate with largest row index, then largest column index.
        return max(candidates, key=lambda item: (item[1], item[2]))

    def split(self, row_index, col_index):
        """
        Return two quasi-ribbon tableaux split from ``self``.

        The first quasi-ribbon tableau includes
        all entries above or weakly left of the (``row_index``, ``col_index``)-position
        and the second includes all entries below or strictly right of that
        position. ``col_index`` indicates the position among integers in a row, ignoring
        the ``None`` entries; both it and ``row_index`` count starting from 0.

        When ``row_index`` exceeds the number of rows, the first quasi-ribbon
        tableau is ``self`` and the second is empty. When ``column_index`` exceeds the
        number of integers in ``self.rows()[row_index]``, the first quasi-ribbon
        contains the entirety of ``self.rows()[row_index]`` and the second starts
        at the row below.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([[1, 2, 3], [4, 5]])
            sage: Q1, Q2 = Q.split(0, 0)
            sage: Q1.rows()
            [[1]]
            sage: Q2.rows()
            [[2, 3], [None, 4, 5]]
            sage: Q3, Q4 = Q.split(1, 0)
            sage: Q3.rows()
            [[1, 2, 3], [None, None, 4]]
            sage: Q4.rows()
            [[5]]
        """
        rows = self._filling

        if len(rows) <= row_index:
            parent = self.parent()
            return self, parent.element_class(parent, [])

        top_rows = rows[:row_index] + [rows[row_index][:col_index + 1]]
        tail = rows[row_index][col_index + 1:]
        if tail:
            bottom_rows = [tail] + rows[row_index + 1:]
        else:
            bottom_rows = rows[row_index + 1:]

        parent = self.parent()
        element_class = parent.element_class
        return (element_class(parent, top_rows),
                element_class(parent, bottom_rows))

    @staticmethod
    def _insert_letter(rows, a):
        r"""
        Insert a letter into a list of rows.

        This returns the result as a list of rows rather than a
        quasi-ribbon tableau.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([])
            sage: Q._insert_letter(Q._filling, 3)
            [[3]]

            sage: Q = QuasiRibbonTableau([[4, 5]])
            sage: Q._insert_letter(Q._filling, 3)
            [[3], [4, 5]]

            sage: Q = QuasiRibbonTableau([[1, 2, 3], [4, 5]])
            sage: Q._insert_letter(Q._filling, 4)
            [[1, 2, 3], [4, 4], [5]]
        """
        if not rows:
            return [[a]]

        site = None
        for row_index, row in enumerate(rows):
            for col_index, entry in enumerate(row):
                if entry <= a:
                    if site is None or (row_index, col_index) > (site[0], site[1]):
                        site = (row_index, col_index)

        if site is None:
            new_rows = [[a]] + rows
        else:
            row_index, col_index = site

            new_rows = rows[:row_index]
            new_rows.append(rows[row_index][:col_index + 1] + [a])

            tail = rows[row_index][col_index + 1:]
            if tail:
                new_rows.append(tail)

            new_rows.extend(rows[row_index + 1:])

        return new_rows

    def insert_letter(self, a):
        """
        Insert one letter a into a quasi-ribbon tableau ``self``.

        Case 1:
        If no entry of ``self`` is less or equal to ``a``, add a as a new top row.

        Case 2:
        Otherwise, find `x`, the rightmost and bottommost entry of ``self`` that is less or equal to `a`.
        Put a immediately to the right of `x`. Then split off the entries that
        were originally to the right of `x` and move them below.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([])
            sage: Q.insert_letter(3)
            [[3]]

            sage: Q = QuasiRibbonTableau([[4, 5]])
            sage: Q.insert_letter(3)
            [[3], [4, 5]]

            sage: Q = QuasiRibbonTableau([[1, 2, 3], [4, 5]])
            sage: Q.insert_letter(4)
            [[1, 2, 3], [None, None, 4, 4], [None, None, None, 5]]
        """
        if not isinstance(a, (int, Integer)) or a <= 0:
            raise ValueError("the inserted letter must be a positive integer")

        parent = self.parent()
        element_class = parent.element_class
        return element_class(parent, self._insert_letter(self._filling, a))

    def insert_word(self, word):
        """
        Insert the letters of a ``word`` one at a time from left to right into ``self``.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([])
            sage: Q.insert_word([3])
            [[3]]

            sage: Q = QuasiRibbonTableau([[4, 5]])
            sage: Q.insert_word([3,3,1,2])
            [[1, 2], [None, 3, 3], [None, None, 4, 5]]

            sage: Q = QuasiRibbonTableau([[1], [2,3]])
            sage: Q.insert_word([])
            [[1], [2, 3]]
            sage: Q.insert_word([4,1,2])
            [[1, 1], [None, 2, 2], [None, None, 3, 4]]
        """
        rows = self._filling
        for a in word:
            if not isinstance(a, (int, Integer)) or a <= 0:
                raise ValueError("the inserted letters must be positive integers")
            rows = self._insert_letter(rows, a)

        parent = self.parent()
        element_class = parent.element_class
        return element_class(parent, rows)

    def _test_quasi_ribbon(self, **options):
        r"""
        Test that ``self`` satisfies the quasi-ribbon tableau conditions.

        More precisely, test that:

        - every row contains at least one non-``None`` entry;
        - the entries in each row are weakly increasing;
        - the entries in each column are strictly increasing; and
        - each row begins one position before the preceding row ends.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableau
            sage: Q = QuasiRibbonTableau([[1, 2, 3],
            ....:                         [None, None, 4, 5]])
            sage: Q._test_quasi_ribbon()

        The test is also run as part of the element test suite::

            sage: TestSuite(Q).run()
        """
        tester = self._tester(**options)
        rows = self.rows()

        expected_shift = 0
        for row_index, row in enumerate(rows):
            shift = 0
            while shift < len(row) and row[shift] is None:
                shift += 1

            entries = row[shift:]

            tester.assertTrue(
                entries,
                f"row {row_index} contains no entries",
            )
            tester.assertTrue(
                all(entry is not None for entry in entries),
                f"row {row_index} contains None after its first entry",
            )
            tester.assertEqual(
                shift,
                expected_shift,
                f"row {row_index} begins in column {shift}, "
                f"but should begin in column {expected_shift}",
            )
            tester.assertTrue(
                all(entries[i] <= entries[i + 1]
                    for i in range(len(entries) - 1)),
                f"row {row_index} is not weakly increasing",
            )

            expected_shift += len(entries) - 1

        width = max((len(row) for row in rows), default=0)
        for col in range(width):
            entries = [
                row[col]
                for row in rows
                if col < len(row) and row[col] is not None
            ]
            tester.assertTrue(
                all(entries[i] < entries[i + 1]
                    for i in range(len(entries) - 1)),
                f"column {col} is not strictly increasing",
            )


class QuasiRibbonTableaux(SkewTableaux):
    r"""
    The set of quasi-ribbon tableaux.

    INPUT:

    - ``shape`` -- (optional) the composition shape of the rows
    - ``max_entry`` -- (optional) the largest allowed entry for finite generation

    EXAMPLES::

        sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
        sage: QRT = QuasiRibbonTableaux()
        sage: QRT
        Quasi-ribbon tableaux
        sage: QRT = QuasiRibbonTableaux(shape=[2, 1])
        sage: QRT
        Quasi-ribbon tableaux of shape [2, 1]
        sage: QRT = QuasiRibbonTableaux(shape=[2, 1], max_entry=3)
        sage: list(QRT)
        [[[1, 1], [None, 2]],
        [[1, 1], [None, 3]],
        [[1, 2], [None, 3]],
        [[2, 2], [None, 3]]]
    """
    @staticmethod
    def __classcall_private__(cls, shape=None, max_entry=None, size=None, category=None):
        """
        Normalize input before constructing the parent object.

        TESTS::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: QRT = QuasiRibbonTableaux()
            sage: QRT
            Quasi-ribbon tableaux
            sage: QRT = QuasiRibbonTableaux(shape=[2, 1])
            sage: QRT
            Quasi-ribbon tableaux of shape [2, 1]
            sage: QuasiRibbonTableaux((3, 1))
            Quasi-ribbon tableaux of shape [3, 1]
            sage: QuasiRibbonTableaux(4)
            Quasi-ribbon tableaux of size 4
            sage: QRT = QuasiRibbonTableaux(shape=[2, 1], max_entry=3)
            sage: QRT
            Quasi-ribbon tableaux of shape [2, 1] with entries at most 3
            sage: QuasiRibbonTableaux(4, max_entry=3)
            Quasi-ribbon tableaux of size 4 with entries at most 3
            sage: TestSuite(QRT).run()
        """
        if shape is not None and size is None:
            if shape in ZZ:
                size = Integer(shape)
                shape = None
            else:
                shape = tuple(shape)
        # Do not allow both shape and size
        if shape is not None and size is not None:
            raise ValueError("specify either shape or size, but not both")
        # Normalize shape
        if shape is not None:
            shape = tuple(shape)
            if any(s <= 0 for s in shape):
                raise ValueError("shape must be a composition with positive parts")
        # Normalize size
        if size is not None:
            size = Integer(size)
            if size < 0:
                raise ValueError("size must be nonnegative")
        # Normalize max_entry
        if max_entry is not None:
            max_entry = Integer(max_entry)
            if max_entry < 0:
                raise ValueError("max_entry must be nonnegative")

        return super().__classcall__(cls, shape=shape,
                                     max_entry=max_entry,
                                     size=size,
                                     category=category)

    def __init__(self, shape=None, max_entry=None, size=None, category=None):
        """
        Initialize ``self``.

        INPUT:

        - ``shape`` -- optional composition shape.
        - ``max_entry`` -- optional largest allowed entry for finite generation.
        - ``category`` -- optional Sage category.

        TESTS::

            sage: QuasiRibbonTableaux(2)
            Quasi-ribbon tableaux of size 2
            sage: QuasiRibbonTableaux([1,2])
            Quasi-ribbon tableaux of shape [1, 2]
            sage: QuasiRibbonTableaux(max_entry=2)
            Quasi-ribbon tableaux with entries at most 2
        """
        # If no category is specified, use the category of sets.
        if category is None:
            if max_entry is not None and (shape is not None or size is not None):
                category = FiniteEnumeratedSets()
            else:
                category = Sets().Infinite()

        # Initialize the inherited SkewTableaux parent class.
        SkewTableaux.__init__(self, category=category)

        # If shape is None, this represents quasiribbon tableaux in general.
        # If shape is given, this represents a fixed-shape family.
        self._shape = None if shape is None else Composition(shape)

        # max_entry is only used when generating finite examples.
        # Without this cutoff, there are infinitely many possible fillings.
        self._max_entry = max_entry
        self._size = size

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: QuasiRibbonTableaux([1,3,1])
            Quasi-ribbon tableaux of shape [1, 3, 1]
            sage: QuasiRibbonTableaux(2)
            Quasi-ribbon tableaux of size 2
        """
        if self._shape is not None:
            if self._max_entry is None:
                return "Quasi-ribbon tableaux of shape {}".format(self._shape)
            return "Quasi-ribbon tableaux of shape {} with entries at most {}".format(
                self._shape, self._max_entry)

        if self._size is not None:
            if self._max_entry is None:
                return "Quasi-ribbon tableaux of size {}".format(self._size)
            return "Quasi-ribbon tableaux of size {} with entries at most {}".format(
                self._size, self._max_entry)

        if self._max_entry is None:
            return "Quasi-ribbon tableaux"

        return "Quasi-ribbon tableaux with entries at most {}".format(self._max_entry)

    Element = QuasiRibbonTableau

    def _element_constructor_(self, rows):
        r"""
        Construct an element of ``self`` from ``rows``.

        The rows are given from top to bottom. Entries equal to ``None`` are
        removed, and all remaining entries are coerced to integers.

        The resulting rows must be nonempty, consist of positive integers,
        be weakly increasing from left to right, and be strictly increasing
        from top to bottom in each overlapping column.

        Internal methods that produce rows already known to be normalized and
        valid may bypass this method by calling
        ``self.element_class(self, rows)`` directly.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: QRT = QuasiRibbonTableaux()
            sage: QRT([[1, 2, 3], [4, 5]])
            [[1, 2, 3], [None, None, 4, 5]]
            sage: QRT([[1, 2, 3], [None, 4, 5]])
            [[1, 2, 3], [None, None, 4, 5]]

        TESTS::

            sage: QRT([1, 2])
            Traceback (most recent call last):
            ...
            TypeError: entries must be positive integers or None
            sage: QRT([[1, -2]])
            Traceback (most recent call last):
            ...
            TypeError: entries must be positive integers or None
            sage: QRT([[None]])
            Traceback (most recent call last):
            ...
            TypeError: a skew tableau cannot have an empty list for a row
            sage: QRT([[2, 1]])
            Traceback (most recent call last):
            ...
            ValueError: rows must be weakly increasing
            sage: QRT([[1, 3], [2]])
            Traceback (most recent call last):
            ...
            ValueError: columns must be strictly increasing
        """
        try:
            clean_rows = [
                [ZZ(entry) for entry in row if entry is not None]
                for row in rows
            ]
        except (TypeError, ValueError):
            raise TypeError("entries must be positive integers or None")

        if any(not row for row in clean_rows):
            raise TypeError("a skew tableau cannot have an empty list for a row")

        if any(entry <= 0 for row in clean_rows for entry in row):
            raise TypeError("entries must be positive integers or None")

        for row in clean_rows:
            if any(row[i] > row[i + 1] for i in range(len(row) - 1)):
                raise ValueError("rows must be weakly increasing")

        for upper, lower in zip(clean_rows, clean_rows[1:]):
            if upper[-1] >= lower[0]:
                raise ValueError("columns must be strictly increasing")

        return self.element_class(self, clean_rows)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: QuasiRibbonTableaux().an_element()
            [[1]]
            sage: QuasiRibbonTableaux(shape=[3, 2]).an_element()
            [[1, 1, 1], [None, None, 2, 2]]
            sage: QuasiRibbonTableaux(size=4).an_element()
            [[1, 1, 1, 1]]
        """
        if self._shape is not None:
            shape = self._shape
        elif self._size is not None:
            if self._size == 0:
                return self([])
            shape = (self._size,)
        else:
            return self([[1]])
        rows = []
        for i, row_length in enumerate(shape):
            rows.append([i + 1] * row_length)

        return self(rows)

    def shape(self):
        """
        Return the fixed composition shape of ``self``.

        Raise an error if this family does not have a fixed shape.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: QRT = QuasiRibbonTableaux(shape=[3, 2])
            sage: QRT.shape()
            [3, 2]
            sage: QRT = QuasiRibbonTableaux(3)
            sage: QRT.shape()
            Traceback (most recent call last):
            ...
            ValueError: this family does not have a fixed shape
        """
        if self._shape is None:
            raise ValueError("this family does not have a fixed shape")
        return self._shape

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: list(QuasiRibbonTableaux(shape=[2, 1], max_entry=2))
            [[[1, 1], [None, 2]]]
            sage: list(QuasiRibbonTableaux(size=2, max_entry=2))
            [[[1, 1]], [[1, 2]], [[2, 2]], [[1], [2]]]
        """
        if self._max_entry is None:
            raise NotImplementedError("iteration requires max_entry to be specified")

        if self._shape is not None:
            S = [self._shape]
        elif self._size is not None:
            S = Compositions(self._size, max_length=self._max_entry)
        else:
            raise NotImplementedError("iteration requires either shape or size to be specified")

        from itertools import accumulate, combinations_with_replacement

        for shape in S:
            shape = tuple(shape)
            r = len(shape)
            n = sum(shape)

            M = self._max_entry - (r - 1)
            if M <= 0:
                continue

            descents = set(accumulate(shape[:-1]))

            offset = [0] * n
            k = 0
            for i in range(1, n):
                if i in descents:
                    k += 1
                offset[i] = k

            for b in combinations_with_replacement(range(1, M + 1), n):
                word = [x + off for x, off in zip(b, offset)]

                rows = []
                pos = 0
                for s in shape:
                    rows.append(word[pos:pos+s])
                    pos += s

                yield self.element_class(self, rows)

    def __contains__(self, x):
        r"""
        Return whether ``x`` is an element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: QRT = QuasiRibbonTableaux()
            sage: [[1, 2], [3]] in QRT
            True

            sage: QRT = QuasiRibbonTableaux(shape=[2, 1])
            sage: [[1, 2], [3]] in QRT
            True
            sage: [[1, 2, 3]] in QRT
            False

            sage: QRT = QuasiRibbonTableaux(shape=[2, 1], max_entry=2)
            sage: [[1, 1], [2]] in QRT
            True
            sage: [[1, 1], [3]] in QRT
            False
        """
        if isinstance(x, QuasiRibbonTableau):
            Q = x
        else:
            try:
                Q = QuasiRibbonTableau(x)
            except (TypeError, ValueError):
                return False

        if (self._shape is not None
                and Q.to_composition() != self._shape):
            return False

        if self._max_entry is not None:
            for row in Q._filling:
                for entry in row:
                    if entry > self._max_entry:
                        return False

        return True

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: QuasiRibbonTableaux().cardinality()
            +Infinity
            sage: QuasiRibbonTableaux(shape=[2, 1]).cardinality()
            +Infinity
            sage: QuasiRibbonTableaux(shape=[2, 1], max_entry=2).cardinality()
            1
            sage: QuasiRibbonTableaux(1, max_entry=2).cardinality()
            2
        """
        if self._max_entry is None or (self._shape is None and self._size is None):
            return infinity
        if self._shape is None:
            n = self._size
            return sum(binomial(n - 1, r - 1) * binomial(self._max_entry + n - r, n)
                       for r in range(1, n+1))

        n = self._shape.size()
        return binomial(self._max_entry + n - len(self._shape), n)

    def insert_word(self, word):
        """
        Insert the letters of ``word`` one at a time into ``self``.

        EXAMPLES::

            sage: from sage.combinat.quasi_ribbon_tableau import QuasiRibbonTableaux
            sage: H = QuasiRibbonTableaux()
            sage: H.insert_word([])
            []
            sage: H.insert_word([3, 4])
            [[3, 4]]
            sage: H.insert_word([4, 3])
            [[3], [4]]
            sage: H.insert_word([1, 3, 2, 4, 4])
            [[1, 2], [None, 3, 4, 4]]
        """
        Q = QuasiRibbonTableau([]).insert_word(word)
        return Q
