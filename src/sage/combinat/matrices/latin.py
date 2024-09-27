# sage.doctest: needs sage.combinat sage.groups sage.modules
r"""
Latin Squares

A *latin square* of order `n` is an `n \times n` array such that
each symbol `s \in \{ 0, 1, \dots, n-1\}` appears precisely once in each
row, and precisely once in each column. A *partial latin square* of
order `n` is an `n \times n` array such that
each symbol `s \in \{ 0, 1, \dots, n-1\}` appears at most once in each
row, and at most once in each column. Empty cells are denoted by `-1`.
A latin square `L` is a
*completion* of a partial latin square `P` if `P \subseteq L`. If
`P` completes to just `L` then `P` *has unique completion*.

A *latin bitrade* `(T_1,\, T_2)` is a pair of partial
latin squares such that:

#. `\{ (i,\,j) \mid (i,\,j,\,k) \in T_1 \text{ for some symbol }k \}
   = \{ (i,\,j) \mid (i,\,j,\,k') \in T_2 \text{ for some symbol }k' \};`

#. for each `(i,\,j,\,k) \in T_1` and `(i,\,j,\,k') \in T_2`,
   `k \neq k'`;

#. the symbols appearing in row `i` of `T_1` are the same as those of
   row `i` of `T_2`; the symbols appearing in column `j` of `T_1` are
   the same as those of column `j` of `T_2`.


Intuitively speaking, a bitrade gives the difference between two latin
squares, so if `(T_1,\, T_2)` is a bitrade
for the pair of latin squares `(L_1,\, L_2)`, then
`L1 = (L2 \setminus T_1) \cup T_2`
and
`L2 = (L1 \setminus T_2) \cup T_1`.

This file contains


#. LatinSquare class definition;

#. some named latin squares (back circulant, forward circulant, abelian
   `2`-group);

#. methods :meth:`is_partial_latin_square` and :meth:`is_latin_square` to test
   if a :class:`LatinSquare` object satisfies the definition of a latin square
   or partial latin square, respectively;

#. tests for completion and unique completion (these use the C++
   implementation of Knuth's dancing links algorithm to solve the
   problem as a instance of `0-1` matrix exact cover);

#. functions for calculating the `\tau_i` representation of a bitrade
   and the genus of the associated hypermap embedding;

#. Markov chain of Jacobson and Matthews (1996) for generating latin
   squares uniformly at random (provides a generator interface);

#. a few examples of `\tau_i` representations of bitrades constructed
   from the action of a group on itself by right multiplication,
   functions for converting to a pair of :class:`LatinSquare` objects.

EXAMPLES::

    sage: from sage.combinat.matrices.latin import *
    sage: B = back_circulant(5)
    sage: B
    [0 1 2 3 4]
    [1 2 3 4 0]
    [2 3 4 0 1]
    [3 4 0 1 2]
    [4 0 1 2 3]
    sage: B.is_latin_square()
    True
    sage: B[0, 1] = 0
    sage: B.is_latin_square()
    False

    sage: (a, b, c, G) = alternating_group_bitrade_generators(1)
    sage: (T1, T2) = bitrade_from_group(a, b, c, G)
    sage: T1
    [ 0 -1  3  1]
    [-1  1  0  2]
    [ 1  3  2 -1]
    [ 2  0 -1  3]
    sage: T2
    [ 1 -1  0  3]
    [-1  0  2  1]
    [ 2  1  3 -1]
    [ 0  3 -1  2]
    sage: T1.nr_filled_cells()
    12
    sage: genus(T1, T2)
    1

.. TODO::

    #. Latin squares with symbols from a ring instead of the integers
       `\{ 0, 1, \dots, n-1 \}`.

    #. Isotopism testing of latin squares and bitrades via graph
       isomorphism (nauty?).

    #. Combinatorial constructions for bitrades.


AUTHORS:

- Carlo Hamalainen (2008-03-23): initial version


TESTS::

    sage: L = elementary_abelian_2group(3)
    sage: L == loads(dumps(L))
    True
"""
# ****************************************************************************
#       Copyright (C) 2008 Carlo Hamalainen <carlo.hamalainen@gmail.com>,
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

from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.groups.perm_gps.constructor import PermutationGroupElement as PermutationConstructor
from sage.libs.gap.libgap import libgap
from sage.libs.gap.element import GapElement
from sage.combinat.permutation import Permutation
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.arith.misc import is_prime
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.misc.flatten import flatten

from .dlxcpp import DLXCPP
from functools import reduce


class LatinSquare:
    def __init__(self, *args):
        """
        Latin squares.

        This class implements a latin square of order n with rows and
        columns indexed by the set 0, 1, ..., n-1 and symbols from the same
        set. The underlying latin square is a matrix(ZZ, n, n). If L is a
        latin square, then the cell at row r, column c is empty if and only
        if L[r, c] < 0. In this way we allow partial latin squares and can
        speak of completions to latin squares, etc.

        There are two ways to declare a latin square:

        Empty latin square of order n::

            sage: n = 3
            sage: L = LatinSquare(n)
            sage: L
            [-1 -1 -1]
            [-1 -1 -1]
            [-1 -1 -1]

        Latin square from a matrix::

            sage: M = matrix(ZZ, [[0, 1], [2, 3]])
            sage: LatinSquare(M)
            [0 1]
            [2 3]
        """
        if len(args) == 1 and isinstance(args[0], (Integer, int)):
            self.square = matrix(ZZ, args[0], args[0])
            self.clear_cells()
        elif len(args) == 2 and all(isinstance(a, (Integer, int))
                                    for a in args):
            self.square = matrix(ZZ, args[0], args[1])
            self.clear_cells()
        elif len(args) == 1 and isinstance(args[0], Matrix_integer_dense):
            self.square = args[0]
        else:
            raise TypeError("bad input for latin square")

    def dumps(self):
        """
        Since the latin square class does not hold any other private
        variables we just call dumps on self.square:

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(2) == loads(dumps(back_circulant(2)))
            True
        """
        from sage.misc.persist import dumps
        return dumps(self.square)

    def __str__(self):
        """
        The string representation of a latin square is the same as the
        underlying matrix.

        EXAMPLES::

            sage: print(LatinSquare(matrix(ZZ, [[0, 1], [2, 3]])).__str__())
            [0 1]
            [2 3]
        """
        return str(self.square)

    def __repr__(self):
        """
        The representation of a latin square is the same as the underlying
        matrix.

        EXAMPLES::

            sage: print(LatinSquare(matrix(ZZ, [[0, 1], [2, 3]])).__repr__())
            [0 1]
            [2 3]
        """
        return repr(self.square)

    def __getitem__(self, rc):
        """
        If L is a LatinSquare then this method allows us to evaluate L[r,
        c].

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(3)
            sage: B[1, 1]
            2
        """

        r = rc[0]
        c = rc[1]

        return self.square[r, c]

    def __setitem__(self, rc, val):
        """
        If L is a LatinSquare then this method allows us to set L[r, c].

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(3)
            sage: B[1, 1] = 10
            sage: B[1, 1]
            10
        """

        r = rc[0]
        c = rc[1]

        self.square[r, c] = val

    def set_immutable(self):
        """
        A latin square is immutable if the underlying matrix is immutable.

        EXAMPLES::

            sage: L = LatinSquare(matrix(ZZ, [[0, 1], [2, 3]]))
            sage: L.set_immutable()
            sage: {L : 0}   # this would fail without set_immutable()
            {[0 1]
            [2 3]: 0}
        """

        self.square.set_immutable()

    def __hash__(self):
        """
        The hash of a latin square is precisely the hash of the underlying
        matrix.

        EXAMPLES::

            sage: L = LatinSquare(matrix(ZZ, [[0, 1], [2, 3]]))
            sage: L.set_immutable()
            sage: L.__hash__()
            1677951251422179082  # 64-bit
            -479138038           # 32-bit
        """
        return hash(self.square)

    def __eq__(self, Q):
        """
        Two latin squares are equal if the underlying matrices are equal.

        EXAMPLES::

            sage: A = LatinSquare(matrix(ZZ, [[0, 1], [2, 3]]))
            sage: B = LatinSquare(matrix(ZZ, [[0, 4], [2, 3]]))
            sage: A == B
            False
            sage: B[0, 1] = 1
            sage: A == B
            True
        """

        return self.square == Q.square

    def __copy__(self):
        """
        To copy a latin square we must copy the underlying matrix.

        EXAMPLES::

            sage: A = LatinSquare(matrix(ZZ, [[0, 1], [2, 3]]))
            sage: B = copy(A)
            sage: B
            [0 1]
            [2 3]
        """
        C = LatinSquare(self.square.nrows(), self.square.ncols())
        from copy import copy
        C.square = copy(self.square)
        return C

    def clear_cells(self):
        """
        Mark every cell in ``self`` as being empty.

        EXAMPLES::

            sage: A = LatinSquare(matrix(ZZ, [[0, 1], [2, 3]]))
            sage: A.clear_cells()
            sage: A
            [-1 -1]
            [-1 -1]
        """
        for r in range(self.square.nrows()):
            for c in range(self.square.ncols()):
                self.square[r, c] = -1

    def nrows(self):
        """
        Number of rows in the latin square.

        EXAMPLES::

            sage: LatinSquare(3).nrows()
            3
        """
        return self.square.nrows()

    def ncols(self):
        """
        Number of columns in the latin square.

        EXAMPLES::

            sage: LatinSquare(3).ncols()
            3
        """
        return self.square.ncols()

    def row(self, x):
        """
        Return row x of the latin square.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(3).row(0)
            (0, 1, 2)
        """
        return self.square.row(x)

    def column(self, x):
        """
        Return column x of the latin square.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(3).column(0)
            (0, 1, 2)
        """
        return self.square.column(x)

    def list(self):
        """
        Convert the latin square into a list, in a row-wise manner.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(3).list()
            [0, 1, 2, 1, 2, 0, 2, 0, 1]
        """
        return self.square.list()

    def nr_filled_cells(self):
        """
        Return the number of filled cells (i.e. cells with a positive
        value) in the partial latin square ``self``.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: LatinSquare(matrix([[0, -1], [-1, 0]])).nr_filled_cells()
            2
        """
        s = 0
        for r in range(self.nrows()):
            for c in range(self.ncols()):
                if self[r, c] >= 0:
                    s += 1
        return s

    def actual_row_col_sym_sizes(self):
        """
        Bitrades sometimes end up in partial latin squares with unused
        rows, columns, or symbols. This function works out the actual
        number of used rows, columns, and symbols.

        .. warning::

           We assume that the unused rows/columns occur in the lower
           right of self, and that the used symbols are in the range
           {0, 1, ..., m} (no holes in that list).

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(3)
            sage: B[0,2] = B[1,2] = B[2,2] = -1
            sage: B[0,0] = B[2,1] = -1
            sage: B
            [-1  1 -1]
            [ 1  2 -1]
            [ 2 -1 -1]
            sage: B.actual_row_col_sym_sizes()
            (3, 2, 2)
        """
        row_max = self.nrows()
        col_max = self.ncols()
        sym_max = self.nr_distinct_symbols()

        while self.is_empty_row(row_max-1):
            row_max -= 1
        while self.is_empty_column(col_max-1):
            col_max -= 1

        return row_max, col_max, sym_max

    def is_empty_column(self, c):
        """
        Check if column c of the partial latin square ``self`` is empty.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: L = back_circulant(4)
            sage: L.is_empty_column(0)
            False
            sage: L[0,0] = L[1,0] = L[2,0] = L[3,0] = -1
            sage: L.is_empty_column(0)
            True
        """
        return list(set(self.column(c))) == [-1]

    def is_empty_row(self, r):
        """
        Check if row r of the partial latin square ``self`` is empty.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: L = back_circulant(4)
            sage: L.is_empty_row(0)
            False
            sage: L[0,0] = L[0,1] = L[0,2] = L[0,3] = -1
            sage: L.is_empty_row(0)
            True
        """
        return list(set(self.row(r))) == [-1]

    def nr_distinct_symbols(self):
        """
        Return the number of distinct symbols in the partial latin square
        ``self``.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(5).nr_distinct_symbols()
            5
            sage: L = LatinSquare(10)
            sage: L.nr_distinct_symbols()
            0
            sage: L[0, 0] = 0
            sage: L[0, 1] = 1
            sage: L.nr_distinct_symbols()
            2
        """
        symbols = set(flatten([list(x) for x in list(self.square)]))
        return sum(1 for x in symbols if x >= 0)

    def apply_isotopism(self, row_perm, col_perm, sym_perm):
        """
        An isotopism is a permutation of the rows, columns, and symbols of
        a partial latin square ``self``. Use isotopism() to convert a tuple
        (indexed from 0) to a Permutation object.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(5)
            sage: B
            [0 1 2 3 4]
            [1 2 3 4 0]
            [2 3 4 0 1]
            [3 4 0 1 2]
            [4 0 1 2 3]
            sage: alpha = isotopism((0,1,2,3,4))
            sage: beta  = isotopism((1,0,2,3,4))
            sage: gamma = isotopism((2,1,0,3,4))
            sage: B.apply_isotopism(alpha, beta, gamma)
            [3 4 2 0 1]
            [0 2 3 1 4]
            [1 3 0 4 2]
            [4 0 1 2 3]
            [2 1 4 3 0]
        """
        Q = LatinSquare(self.nrows(), self.ncols())

        for r in range(self.nrows()):
            for c in range(self.ncols()):
                try:
                    if self[r, c] < 0:
                        s2 = -1
                    else:
                        s2 = sym_perm[self[r, c]] - 1
                except IndexError:
                    s2 = self[r, c]  # we must be leaving the symbol fixed?

                Q[row_perm[r]-1, col_perm[c]-1] = s2

        return Q

    def filled_cells_map(self):
        """
        Number the filled cells of ``self`` with integers from {1, 2, 3, ...}.

        INPUT:

        - ``self`` -- partial latin square ``self`` (empty cells
          have negative values)

        OUTPUT:

        A dictionary ``cells_map`` where ``cells_map[(i,j)] = m`` means that
        ``(i,j)`` is the ``m``-th filled cell in ``P``,
        while ``cells_map[m] = (i,j)``.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: (a, b, c, G) = alternating_group_bitrade_generators(1)
            sage: (T1, T2) = bitrade_from_group(a, b, c, G)
            sage: D = T1.filled_cells_map()
            sage: {i: v for i,v in D.items() if i in ZZ}
            {1: (0, 0),
             2: (0, 2),
             3: (0, 3),
             4: (1, 1),
             5: (1, 2),
             6: (1, 3),
             7: (2, 0),
             8: (2, 1),
             9: (2, 2),
             10: (3, 0),
             11: (3, 1),
             12: (3, 3)}
            sage: {i: v for i,v in D.items() if i not in ZZ}
            {(0, 0): 1,
             (0, 2): 2,
             (0, 3): 3,
             (1, 1): 4,
             (1, 2): 5,
             (1, 3): 6,
             (2, 0): 7,
             (2, 1): 8,
             (2, 2): 9,
             (3, 0): 10,
             (3, 1): 11,
             (3, 3): 12}
        """
        cells_map = {}
        k = 1

        for r in range(self.nrows()):
            for c in range(self.ncols()):
                e = self[r, c]

                if e < 0:
                    continue

                cells_map[(r, c)] = k
                cells_map[k] = (r, c)

                k += 1

        return cells_map

    def top_left_empty_cell(self):
        """
        Return the least ``[r, c]`` such that ``self[r, c]`` is an empty cell.
        If all cells are filled then we return ``None``.

        INPUT:

        - ``self`` -- LatinSquare

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(5)
            sage: B[3, 4] = -1
            sage: B.top_left_empty_cell()
            [3, 4]
        """
        for r in range(self.nrows()):
            for c in range(self.ncols()):
                if self[r, c] < 0:
                    return [r, c]

        return None

    def is_partial_latin_square(self):
        """
        ``self`` is a partial latin square if it is an n by n matrix, and each
        symbol in [0, 1, ..., n-1] appears at most once in each row, and at
        most once in each column.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: LatinSquare(4).is_partial_latin_square()
            True
            sage: back_circulant(3).gcs().is_partial_latin_square()
            True
            sage: back_circulant(6).is_partial_latin_square()
            True
        """

        assert self.nrows() == self.ncols()

        n = self.nrows()

        for r in range(n):
            vals_in_row = {}

            for c in range(n):
                e = self[r, c]

                if e < 0:
                    continue

                # Entry out of range 0, 1, ..., n-1:
                if e >= n:
                    return False

                # Entry has already appeared in this row:
                if e in vals_in_row:
                    return False

                vals_in_row[e] = True

        for c in range(n):
            vals_in_col = {}

            for r in range(n):
                e = self[r, c]

                if e < 0:
                    continue

                # Entry out of range 0, 1, ..., n-1:
                if e >= n:
                    return False

                # Entry has already appeared in this column:
                if e in vals_in_col:
                    return False

                vals_in_col[e] = True

        return True

    def is_latin_square(self):
        """
        ``self`` is a latin square if it is an n by n matrix, and each symbol
        in [0, 1, ..., n-1] appears exactly once in each row, and exactly
        once in each column.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: elementary_abelian_2group(4).is_latin_square()
            True

        ::

            sage: forward_circulant(7).is_latin_square()
            True
        """
        # We do not allow latin rectangles:
        if self.nrows() != self.ncols():
            return False

        # Every cell must be filled:
        if any(x < 0 for x in self.list()):
            return False

        # By necessity self must be a partial latin square:
        if not self.is_partial_latin_square():
            return False

        return True

    def permissable_values(self, r, c):
        """
        Find all values that do not appear in row r and column c of the
        latin square ``self``. If ``self[r, c]`` is filled then we return the
        empty list.

        INPUT:

        - ``self`` -- LatinSquare

        - ``r`` -- integer; row of the latin square

        - ``c`` -- integer; column of the latin square

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: L = back_circulant(5)
            sage: L[0, 0] = -1
            sage: L.permissable_values(0, 0)
            [0]
        """

        if self[r, c] >= 0:
            return []

        assert self.nrows() == self.ncols()

        n = self.nrows()

        vals = {}
        for e in range(n):
            vals[e] = True

        for i in range(n):
            if self[i, c] >= 0:
                del vals[self[i, c]]

        for j in range(n):
            if self[r, j] >= 0:
                try:
                    del vals[self[r, j]]
                except KeyError:
                    # We may have already removed a symbol
                    # in the previous for-loop.
                    pass

        return list(vals)

    def random_empty_cell(self):
        """
        Find an empty cell of self, uniformly at random.

        INPUT:

        - ``self`` -- LatinSquare

        OUTPUT:

        - ``[r, c]`` -- cell such that ``self[r, c]`` is empty, or returns
          ``None`` if ``self`` is a (full) latin square

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: P = back_circulant(2)
            sage: P[1,1] = -1
            sage: P.random_empty_cell()
            [1, 1]
        """

        cells = {}

        for r in range(self.nrows()):
            for c in range(self.ncols()):
                if self[r, c] < 0:
                    cells[(r, c)] = True

        cells = list(cells)

        if not cells:
            return None

        rc = cells[ZZ.random_element(len(cells))]

        return [rc[0], rc[1]]

    def is_uniquely_completable(self):
        """
        Return ``True`` if the partial latin square ``self`` has exactly one
        completion to a latin square. This is just a wrapper for the
        current best-known algorithm, Dancing Links by Knuth. See
        dancing_links.spyx

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(4).gcs().is_uniquely_completable()
            True

        ::

            sage: G = elementary_abelian_2group(3).gcs()
            sage: G.is_uniquely_completable()
            True

        ::

            sage: G[0, 0] = -1
            sage: G.is_uniquely_completable()
            False
        """

        return self.dlxcpp_has_unique_completion()

    def is_completable(self):
        """
        Return ``True`` if the partial latin square can be completed to a
        latin square.

        EXAMPLES:

        The following partial latin square has no completion because there
        is nowhere that we can place the symbol 0 in the third row::

            sage: B = LatinSquare(3)

        ::

            sage: B[0, 0] = 0
            sage: B[1, 1] = 0
            sage: B[2, 2] = 1

        ::

            sage: B
            [ 0 -1 -1]
            [-1  0 -1]
            [-1 -1  1]

        ::

            sage: B.is_completable()
            False

        ::

            sage: B[2, 2] = 0
            sage: B.is_completable()
            True
        """
        return bool(dlxcpp_find_completions(self, nr_to_find=1))

    def gcs(self):
        """
        A greedy critical set of a latin square ``self`` is found by
        successively removing elements in a row-wise (bottom-up) manner,
        checking for unique completion at each step.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: A = elementary_abelian_2group(3)
            sage: G = A.gcs()
            sage: A
            [0 1 2 3 4 5 6 7]
            [1 0 3 2 5 4 7 6]
            [2 3 0 1 6 7 4 5]
            [3 2 1 0 7 6 5 4]
            [4 5 6 7 0 1 2 3]
            [5 4 7 6 1 0 3 2]
            [6 7 4 5 2 3 0 1]
            [7 6 5 4 3 2 1 0]
            sage: G
            [ 0  1  2  3  4  5  6 -1]
            [ 1  0  3  2  5  4 -1 -1]
            [ 2  3  0  1  6 -1  4 -1]
            [ 3  2  1  0 -1 -1 -1 -1]
            [ 4  5  6 -1  0  1  2 -1]
            [ 5  4 -1 -1  1  0 -1 -1]
            [ 6 -1  4 -1  2 -1  0 -1]
            [-1 -1 -1 -1 -1 -1 -1 -1]
        """

        n = self.nrows()

        from copy import copy
        G = copy(self)

        for r in range(n-1, -1, -1):
            for c in range(n-1, -1, -1):
                e = G[r, c]
                G[r, c] = -1

                if not G.dlxcpp_has_unique_completion():
                    G[r, c] = e

        return G

    def dlxcpp_has_unique_completion(self):
        """
        Check if the partial latin square ``self`` of order n can be embedded
        in precisely one latin square of order n.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(2).dlxcpp_has_unique_completion()
            True
            sage: P = LatinSquare(2)
            sage: P.dlxcpp_has_unique_completion()
            False
            sage: P[0, 0] = 0
            sage: P.dlxcpp_has_unique_completion()
            True
        """
        return len(dlxcpp_find_completions(self, nr_to_find=2)) == 1

    def vals_in_row(self, r):
        """
        Return a dictionary with key e if and only if row r of ``self`` has
        the symbol e.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(3)
            sage: B[0, 0] = -1
            sage: back_circulant(3).vals_in_row(0)
            {0: True, 1: True, 2: True}
        """

        n = self.ncols()
        vals_in_row = {}

        for c in range(n):
            e = self[r, c]
            if e >= 0:
                vals_in_row[e] = True

        return vals_in_row

    def vals_in_col(self, c):
        """
        Return a dictionary with key e if and only if column c of ``self`` has
        the symbol e.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(3)
            sage: B[0, 0] = -1
            sage: back_circulant(3).vals_in_col(0)
            {0: True, 1: True, 2: True}
        """
        n = self.nrows()
        vals_in_col = {}

        for r in range(n):
            e = self[r, c]
            if e >= 0:
                vals_in_col[e] = True

        return vals_in_col

    def latex(self):
        r"""
        Return LaTeX code for the latin square.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: print(back_circulant(3).latex())
            \begin{array}{|c|c|c|}\hline 0 & 1 & 2\\\hline 1 & 2 & 0\\\hline 2 & 0 & 1\\\hline\end{array}
        """

        a = ""
        a += r"\begin{array}{" + self.ncols()*"|c" + "|}"
        for r in range(self.nrows()):
            a += r"\hline "
            for c in range(self.ncols()):
                s = self[r, c]
                if s < 0:
                    a += "~"
                else:
                    a += str(s)

                if c < self.ncols()-1:
                    a += " & "
                else:
                    a += "\\\\"
        a += r"\hline"
        a += r"\end{array}"
        return a

    def disjoint_mate_dlxcpp_rows_and_map(self, allow_subtrade):
        """
        Internal function for find_disjoint_mates.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(4)
            sage: B.disjoint_mate_dlxcpp_rows_and_map(allow_subtrade = True)
            ([[0, 16, 32],
              [1, 17, 32],
              [2, 18, 32],
              [3, 19, 32],
              [4, 16, 33],
              [5, 17, 33],
              [6, 18, 33],
              [7, 19, 33],
              [8, 16, 34],
              [9, 17, 34],
              [10, 18, 34],
              [11, 19, 34],
              [12, 16, 35],
              [13, 17, 35],
              [14, 18, 35],
              [15, 19, 35],
              [0, 20, 36],
              [1, 21, 36],
              [2, 22, 36],
              [3, 23, 36],
              [4, 20, 37],
              [5, 21, 37],
              [6, 22, 37],
              [7, 23, 37],
              [8, 20, 38],
              [9, 21, 38],
              [10, 22, 38],
              [11, 23, 38],
              [12, 20, 39],
              [13, 21, 39],
              [14, 22, 39],
              [15, 23, 39],
              [0, 24, 40],
              [1, 25, 40],
              [2, 26, 40],
              [3, 27, 40],
              [4, 24, 41],
              [5, 25, 41],
              [6, 26, 41],
              [7, 27, 41],
              [8, 24, 42],
              [9, 25, 42],
              [10, 26, 42],
              [11, 27, 42],
              [12, 24, 43],
              [13, 25, 43],
              [14, 26, 43],
              [15, 27, 43],
              [0, 28, 44],
              [1, 29, 44],
              [2, 30, 44],
              [3, 31, 44],
              [4, 28, 45],
              [5, 29, 45],
              [6, 30, 45],
              [7, 31, 45],
              [8, 28, 46],
              [9, 29, 46],
              [10, 30, 46],
              [11, 31, 46],
              [12, 28, 47],
              [13, 29, 47],
              [14, 30, 47],
              [15, 31, 47]],
             {(0, 16, 32): (0, 0, 0),
              (0, 20, 36): (1, 0, 0),
              (0, 24, 40): (2, 0, 0),
              (0, 28, 44): (3, 0, 0),
              (1, 17, 32): (0, 0, 1),
              (1, 21, 36): (1, 0, 1),
              (1, 25, 40): (2, 0, 1),
              (1, 29, 44): (3, 0, 1),
              (2, 18, 32): (0, 0, 2),
              (2, 22, 36): (1, 0, 2),
              (2, 26, 40): (2, 0, 2),
              (2, 30, 44): (3, 0, 2),
              (3, 19, 32): (0, 0, 3),
              (3, 23, 36): (1, 0, 3),
              (3, 27, 40): (2, 0, 3),
              (3, 31, 44): (3, 0, 3),
              (4, 16, 33): (0, 1, 0),
              (4, 20, 37): (1, 1, 0),
              (4, 24, 41): (2, 1, 0),
              (4, 28, 45): (3, 1, 0),
              (5, 17, 33): (0, 1, 1),
              (5, 21, 37): (1, 1, 1),
              (5, 25, 41): (2, 1, 1),
              (5, 29, 45): (3, 1, 1),
              (6, 18, 33): (0, 1, 2),
              (6, 22, 37): (1, 1, 2),
              (6, 26, 41): (2, 1, 2),
              (6, 30, 45): (3, 1, 2),
              (7, 19, 33): (0, 1, 3),
              (7, 23, 37): (1, 1, 3),
              (7, 27, 41): (2, 1, 3),
              (7, 31, 45): (3, 1, 3),
              (8, 16, 34): (0, 2, 0),
              (8, 20, 38): (1, 2, 0),
              (8, 24, 42): (2, 2, 0),
              (8, 28, 46): (3, 2, 0),
              (9, 17, 34): (0, 2, 1),
              (9, 21, 38): (1, 2, 1),
              (9, 25, 42): (2, 2, 1),
              (9, 29, 46): (3, 2, 1),
              (10, 18, 34): (0, 2, 2),
              (10, 22, 38): (1, 2, 2),
              (10, 26, 42): (2, 2, 2),
              (10, 30, 46): (3, 2, 2),
              (11, 19, 34): (0, 2, 3),
              (11, 23, 38): (1, 2, 3),
              (11, 27, 42): (2, 2, 3),
              (11, 31, 46): (3, 2, 3),
              (12, 16, 35): (0, 3, 0),
              (12, 20, 39): (1, 3, 0),
              (12, 24, 43): (2, 3, 0),
              (12, 28, 47): (3, 3, 0),
              (13, 17, 35): (0, 3, 1),
              (13, 21, 39): (1, 3, 1),
              (13, 25, 43): (2, 3, 1),
              (13, 29, 47): (3, 3, 1),
              (14, 18, 35): (0, 3, 2),
              (14, 22, 39): (1, 3, 2),
              (14, 26, 43): (2, 3, 2),
              (14, 30, 47): (3, 3, 2),
              (15, 19, 35): (0, 3, 3),
              (15, 23, 39): (1, 3, 3),
              (15, 27, 43): (2, 3, 3),
              (15, 31, 47): (3, 3, 3)})
        """

        assert self.nrows() == self.ncols()

        n = self.nrows()

        # We will need 3n^2 columns in total:
        #
        # n^2 for the xCy columns
        # n^2 for the xRy columns
        # n^2 for the xy columns

        dlx_rows = []
        cmap = {}

        max_column_nr = -1

        for r in range(n):
            valsrow = self.vals_in_row(r)

            for c in range(n):
                valscol = self.vals_in_col(c)

                # If this is an empty cell of self then we do nothing.
                if self[r, c] < 0:
                    continue

                for e in sorted(set(list(valsrow) + list(valscol))):
                    # These should be constants
                    c_OFFSET = e + c*n
                    r_OFFSET = e + r*n + n*n
                    xy_OFFSET = 2*n*n + r*n + c

                    cmap[(c_OFFSET, r_OFFSET, xy_OFFSET)] = (r, c, e)

                    # The disjoint mate has to be disjoint.
                    if (not allow_subtrade) and self[r, c] == e:
                        continue

                    # The permissible symbols must come from this row/column.
                    if e not in valsrow:
                        continue
                    if e not in valscol:
                        continue

                    dlx_rows.append([c_OFFSET, r_OFFSET, xy_OFFSET])

                    max_column_nr = max(max_column_nr, max(c_OFFSET, r_OFFSET, xy_OFFSET))

        # We will have missed some columns. We
        # have to add 'dummy' rows so that the C++ DLX solver will find
        # a solution.
        used_columns = flatten(dlx_rows)
        for i in range(max_column_nr + 1):
            if i not in used_columns:
                dlx_rows.append([i])

        return dlx_rows, cmap

    def find_disjoint_mates(self, nr_to_find=None, allow_subtrade=False):
        r"""
        .. warning::

            If allow_subtrade is ``True`` then we may return a partial
            latin square that is *not* disjoint to ``self``. In that case,
            use bitrade(P, Q) to get an actual bitrade.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(4)
            sage: g = B.find_disjoint_mates(allow_subtrade = True)
            sage: B1 = next(g)
            sage: B0, B1 = bitrade(B, B1)
            sage: assert is_bitrade(B0, B1)
            sage: print(B0)
            [-1  1  2 -1]
            [-1  2 -1  0]
            [-1 -1 -1 -1]
            [-1  0  1  2]
            sage: print(B1)
            [-1  2  1 -1]
            [-1  0 -1  2]
            [-1 -1 -1 -1]
            [-1  1  2  0]
        """

        assert self.nrows() == self.ncols()

        dlx_rows, cmap = self.disjoint_mate_dlxcpp_rows_and_map(allow_subtrade)

        nr_found = 0

        for x in DLXCPP(dlx_rows):
            nr_found += 1

            from copy import deepcopy
            Q = deepcopy(self)

            for y in x:
                if len(dlx_rows[y]) == 1:
                    continue  # dummy row
                (r, c, e) = cmap[tuple(dlx_rows[y])]
                Q[r, c] = e

            yield Q

            if nr_to_find is not None and nr_found >= nr_to_find:
                return

    def contained_in(self, Q):
        r"""
        Return ``True`` if ``self`` is a subset of `Q`.

        EXAMPLES::

            sage: from sage.combinat.matrices.latin import *
            sage: P = elementary_abelian_2group(2)
            sage: P[0, 0] = -1
            sage: P.contained_in(elementary_abelian_2group(2))
            True
            sage: back_circulant(4).contained_in(elementary_abelian_2group(2))
            False
        """
        for r in range(self.nrows()):
            for c in range(self.ncols()):
                if self[r, c] >= 0 and Q[r, c] < 0:
                    return False
                if self[r, c] >= 0 and self[r, c] != Q[r, c]:
                    return False
        return True


def genus(T1, T2):
    """
    Return the genus of hypermap embedding associated with the bitrade
    (T1, T2).

    Informally, we compute the [tau_1, tau_2, tau_3]
    permutation representation of the bitrade. Each cycle of tau_1,
    tau_2, and tau_3 gives a rotation scheme for a black, white, and
    star vertex (respectively). The genus then comes from Euler's
    formula.

    For more details see Carlo Hamalainen: *Partitioning
    3-homogeneous latin bitrades*. To appear in Geometriae Dedicata,
    available at :arxiv:`0710.0938`

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: (a, b, c, G) = alternating_group_bitrade_generators(1)
        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: genus(T1, T2)
        1
        sage: (a, b, c, G) = pq_group_bitrade_generators(3, 7)
        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: genus(T1, T2)
        3
    """
    cells_map, t1, t2, t3 = tau123(T1, T2)
    return (len(t1.to_cycles()) + len(t2.to_cycles()) + len(t3.to_cycles()) - T1.nr_filled_cells() - 2) // (-2)


def tau123(T1, T2):
    r"""
    Compute the tau_i representation for a bitrade (T1, T2).

    See the
    functions tau1, tau2, and tau3 for the mathematical definitions.

    OUTPUT:

    - (cells_map, t1, t2, t3)

    where cells_map is a map to/from the filled cells of T1, and t1,
    t2, t3 are the tau1, tau2, tau3 permutations.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: (a, b, c, G) = pq_group_bitrade_generators(3, 7)
        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: T1
        [ 0  1  3 -1 -1 -1 -1]
        [ 1  2  4 -1 -1 -1 -1]
        [ 2  3  5 -1 -1 -1 -1]
        [ 3  4  6 -1 -1 -1 -1]
        [ 4  5  0 -1 -1 -1 -1]
        [ 5  6  1 -1 -1 -1 -1]
        [ 6  0  2 -1 -1 -1 -1]
        sage: T2
        [ 1  3  0 -1 -1 -1 -1]
        [ 2  4  1 -1 -1 -1 -1]
        [ 3  5  2 -1 -1 -1 -1]
        [ 4  6  3 -1 -1 -1 -1]
        [ 5  0  4 -1 -1 -1 -1]
        [ 6  1  5 -1 -1 -1 -1]
        [ 0  2  6 -1 -1 -1 -1]
        sage: (cells_map, t1, t2, t3) = tau123(T1, T2)
        sage: D = cells_map
        sage: {i: v for i,v in D.items() if i in ZZ}
        {1: (0, 0),
         2: (0, 1),
         3: (0, 2),
         4: (1, 0),
         5: (1, 1),
         6: (1, 2),
         7: (2, 0),
         8: (2, 1),
         9: (2, 2),
         10: (3, 0),
         11: (3, 1),
         12: (3, 2),
         13: (4, 0),
         14: (4, 1),
         15: (4, 2),
         16: (5, 0),
         17: (5, 1),
         18: (5, 2),
         19: (6, 0),
         20: (6, 1),
         21: (6, 2)}
        sage: {i: v for i,v in D.items() if i not in ZZ}
        {(0, 0): 1,
         (0, 1): 2,
         (0, 2): 3,
         (1, 0): 4,
         (1, 1): 5,
         (1, 2): 6,
         (2, 0): 7,
         (2, 1): 8,
         (2, 2): 9,
         (3, 0): 10,
         (3, 1): 11,
         (3, 2): 12,
         (4, 0): 13,
         (4, 1): 14,
         (4, 2): 15,
         (5, 0): 16,
         (5, 1): 17,
         (5, 2): 18,
         (6, 0): 19,
         (6, 1): 20,
         (6, 2): 21}
        sage: cells_map_as_square(cells_map, max(T1.nrows(), T1.ncols()))
        [ 1  2  3 -1 -1 -1 -1]
        [ 4  5  6 -1 -1 -1 -1]
        [ 7  8  9 -1 -1 -1 -1]
        [10 11 12 -1 -1 -1 -1]
        [13 14 15 -1 -1 -1 -1]
        [16 17 18 -1 -1 -1 -1]
        [19 20 21 -1 -1 -1 -1]
        sage: t1
        [3, 1, 2, 6, 4, 5, 9, 7, 8, 12, 10, 11, 15, 13, 14, 18, 16, 17, 21, 19, 20]
        sage: t2
        [4, 8, 15, 7, 11, 18, 10, 14, 21, 13, 17, 3, 16, 20, 6, 19, 2, 9, 1, 5, 12]
        sage: t3
        [20, 18, 10, 2, 21, 13, 5, 3, 16, 8, 6, 19, 11, 9, 1, 14, 12, 4, 17, 15, 7]

    ::

        sage: t1.to_cycles()
        [(1, 3, 2), (4, 6, 5), (7, 9, 8), (10, 12, 11), (13, 15, 14), (16, 18, 17), (19, 21, 20)]
        sage: t2.to_cycles()
        [(1, 4, 7, 10, 13, 16, 19), (2, 8, 14, 20, 5, 11, 17), (3, 15, 6, 18, 9, 21, 12)]
        sage: t3.to_cycles()
        [(1, 20, 15), (2, 18, 4), (3, 10, 8), (5, 21, 7), (6, 13, 11), (9, 16, 14), (12, 19, 17)]

    The product t1\*t2\*t3 is the identity, i.e. it fixes every point::

        sage: len((t1*t2*t3).fixed_points()) == T1.nr_filled_cells()
        True
    """
    assert is_bitrade(T1, T2)

    cells_map = T1.filled_cells_map()

    t1 = tau1(T1, T2, cells_map)
    t2 = tau2(T1, T2, cells_map)
    t3 = tau3(T1, T2, cells_map)

    return (cells_map, t1, t2, t3)


def isotopism(p):
    """
    Return a Permutation object that represents an isotopism (for rows,
    columns or symbols of a partial latin square).

    Technically, all this function does is take as input a
    representation of a permutation of `0,...,n-1` and return a
    :class:`Permutation` object defined on `1,...,n`.

    For a definition of isotopism, see the :wikipedia:`wikipedia section on
    isotopism <Latin_square#Equivalence_classes_of_Latin_squares>`.

    INPUT:

    According to the type of input (see examples below):

    - an integer `n` -- the function returns the identity on `1,...,n`

    - a string representing a permutation in disjoint cycles notation,
      e.g. `(0,1,2)(3,4,5)` -- the corresponding permutation is returned,
      shifted by 1 to act on `1,...,n`

    - list/tuple of tuples -- assumes disjoint cycle notation, see previous
      entry

    - a list of integers -- the function adds `1` to each member of the
      list, and returns the corresponding permutation

    - a :class:`PermutationGroupElement` ``p`` -- returns a permutation
      describing ``p`` **without** any shift

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: isotopism(5) # identity on 5 points
        [1, 2, 3, 4, 5]

    ::

        sage: G = PermutationGroup(['(1,2,3)(4,5)'])
        sage: g = G.gen(0)
        sage: isotopism(g)
        [2, 3, 1, 5, 4]

    ::

        sage: isotopism([0,3,2,1]) # 0 goes to 0, 1 goes to 3, etc.
        [1, 4, 3, 2]

    ::

        sage: isotopism( (0,1,2) ) # single cycle, presented as a tuple
        [2, 3, 1]

    ::

        sage: x = isotopism( ((0,1,2), (3,4)) ) # tuple of cycles
        sage: x
        [2, 3, 1, 5, 4]
        sage: x.to_cycles()
        [(1, 2, 3), (4, 5)]
    """

    # Identity isotopism on p points:
    if isinstance(p, (Integer, int)):
        return Permutation(range(1, p + 1))

    if isinstance(p, PermutationGroupElement):
        # fixme Ask the Sage mailing list about the tuple/list issue!
        return Permutation(list(p.tuple()))

    if isinstance(p, list):
        # We expect a list like [0,3,2,1] which means
        # that 0 goes to 0, 1 goes to 3, etc.
        return Permutation([x+1 for x in p])

    if isinstance(p, tuple):
        # We have a single cycle:
        if isinstance(p[0], Integer):
            return Permutation(tuple(x+1 for x in p))

        # We have a tuple of cycles:
        if isinstance(p[0], tuple):
            x = isotopism(p[0])

            for i in range(1, len(p)):
                x = x._left_to_right_multiply_on_left(isotopism(p[i]))

            return x

    # Not sure what we got!
    raise TypeError("unable to convert {!r} to isotopism".format(p))


def cells_map_as_square(cells_map, n):
    """
    Return a LatinSquare with cells numbered from 1, 2, ... to given
    the dictionary cells_map.

    .. NOTE::

       The value n should be the maximum of the number of rows and
       columns of the original partial latin square

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: (a, b, c, G) = alternating_group_bitrade_generators(1)
        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: T1
        [ 0 -1  3  1]
        [-1  1  0  2]
        [ 1  3  2 -1]
        [ 2  0 -1  3]

    There are 12 filled cells in T::

        sage: cells_map_as_square(T1.filled_cells_map(), max(T1.nrows(), T1.ncols()))
        [ 1 -1  2  3]
        [-1  4  5  6]
        [ 7  8  9 -1]
        [10 11 -1 12]
    """

    assert n > 1

    L = LatinSquare(n, n)

    for r in range(n):
        for c in range(n):
            try:
                L[r, c] = cells_map[(r, c)]
            except KeyError:
                # There is no cell (r,c) so skip it
                L[r, c] = -1

    return L


def beta1(rce, T1, T2):
    """
    Find the unique (x, c, e) in T2 such that (r, c, e) is in T1.

    INPUT:

    - ``rce`` -- tuple (or list) (r, c, e) in T1

    - ``T1``, ``T2`` -- latin bitrade

    OUTPUT: (x, c, e) in T2

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: is_bitrade(T1, T2)
        True
        sage: beta1([0, 0, 0], T1, T2)
        (1, 0, 0)
    """
    r = rce[0]
    c = rce[1]
    e = rce[2]

    assert T1[r, c] == e
    assert e >= 0

    for x in range(T1.nrows()):
        if T2[x, c] == e:
            return (x, c, e)

    raise ValueError


def beta2(rce, T1, T2):
    """
    Find the unique (r, x, e) in T2 such that (r, c, e) is in T1.

    INPUT:

    - ``rce`` -- tuple (or list) (r, c, e) in T1

    - ``T1``, ``T2`` -- latin bitrade

    OUTPUT: (r, x, e) in T2

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: is_bitrade(T1, T2)
        True
        sage: beta2([0, 0, 0], T1, T2)
        (0, 1, 0)
    """
    r = rce[0]
    c = rce[1]
    e = rce[2]

    assert T1[r, c] == e
    assert e >= 0

    for x in range(T1.ncols()):
        if T2[r, x] == e:
            return (r, x, e)

    raise ValueError


def beta3(rce, T1, T2):
    """
    Find the unique (r, c, x) in T2 such that (r, c, e) is in T1.

    INPUT:

    - ``rce`` -- tuple (or list) (r, c, e) in T1

    - ``T1, T2`` -- latin bitrade

    OUTPUT: (r, c, x) in T2.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: is_bitrade(T1, T2)
        True
        sage: beta3([0, 0, 0], T1, T2)
        (0, 0, 4)
    """
    r = rce[0]
    c = rce[1]
    e = rce[2]

    assert T1[r, c] == e
    assert e >= 0

    # fixme this range could be iffy if we
    # work with latin bitrade rectangles...
    for x in range(T1.nrows()):
        if T2[r, c] == x:
            return (r, c, x)

    raise ValueError


def tau1(T1, T2, cells_map):
    r"""
    The definition of `\tau_1` is

    .. MATH::

       \tau_1 : T1 \rightarrow T1 \\
       \tau_1 = \beta_2^{-1} \beta_3

    where the composition is left to right and `\beta_i : T2 \rightarrow T1`
    changes just the `i`-th coordinate of a triple.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: is_bitrade(T1, T2)
        True
        sage: (cells_map, t1, t2, t3) = tau123(T1, T2)
        sage: t1 = tau1(T1, T2, cells_map)
        sage: t1
        [2, 3, 4, 5, 1, 7, 8, 9, 10, 6, 12, 13, 14, 15, 11, 17, 18, 19, 20, 16, 22, 23, 24, 25, 21]
        sage: t1.to_cycles()
        [(1, 2, 3, 4, 5), (6, 7, 8, 9, 10), (11, 12, 13, 14, 15), (16, 17, 18, 19, 20), (21, 22, 23, 24, 25)]
    """

    # The cells_map has both directions, i.e. integer to
    # cell and cell to integer, so the size of T1 is
    # just half of len(cells_map).
    x = (int(len(cells_map)/2) + 1) * [-1]

    for r in range(T1.nrows()):
        for c in range(T1.ncols()):
            e = T1[r, c]

            if e < 0:
                continue

            (r2, c2, e2) = beta2((r, c, e), T1, T2)
            (r3, c3, e3) = beta3((r2, c2, e2), T2, T1)

            x[cells_map[(r, c)]] = cells_map[(r3, c3)]

    x.pop(0)
    # remove the head of the list since we
    # have permutations on 1..(something).

    return Permutation(x)


def tau2(T1, T2, cells_map):
    r"""
    The definition of `\tau_2` is

    .. MATH::

       \tau_2 : T1 \rightarrow T1 \\
       \tau_2 = \beta_3^{-1} \beta_1

    where the composition is left to right and `\beta_i : T2 \rightarrow T1`
    changes just the `i`-th coordinate of a triple.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: is_bitrade(T1, T2)
        True
        sage: (cells_map, t1, t2, t3) = tau123(T1, T2)
        sage: t2 = tau2(T1, T2, cells_map)
        sage: t2
        [21, 22, 23, 24, 25, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        sage: t2.to_cycles()
        [(1, 21, 16, 11, 6), (2, 22, 17, 12, 7), (3, 23, 18, 13, 8), (4, 24, 19, 14, 9), (5, 25, 20, 15, 10)]
    """

    # The cells_map has both directions, i.e. integer to
    # cell and cell to integer, so the size of T1 is
    # just half of len(cells_map).
    x = (int(len(cells_map)/2) + 1) * [-1]

    for r in range(T1.nrows()):
        for c in range(T1.ncols()):
            e = T1[r, c]

            if e < 0:
                continue

            (r2, c2, e2) = beta3((r, c, e), T1, T2)
            (r3, c3, e3) = beta1((r2, c2, e2), T2, T1)

            x[cells_map[(r, c)]] = cells_map[(r3, c3)]

    x.pop(0)
    # remove the head of the list since we
    # have permutations on 1..(something).

    return Permutation(x)


def tau3(T1, T2, cells_map):
    r"""
    The definition of `\tau_3` is

    .. MATH::

       \tau_3 : T1 \rightarrow T1 \\
       \tau_3 = \beta_1^{-1} \beta_2

    where the composition is left to right and `\beta_i : T2 \rightarrow T1`
    changes just the `i`-th coordinate of a triple.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: is_bitrade(T1, T2)
        True
        sage: (cells_map, t1, t2, t3) = tau123(T1, T2)
        sage: t3 = tau3(T1, T2, cells_map)
        sage: t3
        [10, 6, 7, 8, 9, 15, 11, 12, 13, 14, 20, 16, 17, 18, 19, 25, 21, 22, 23, 24, 5, 1, 2, 3, 4]
        sage: t3.to_cycles()
        [(1, 10, 14, 18, 22), (2, 6, 15, 19, 23), (3, 7, 11, 20, 24), (4, 8, 12, 16, 25), (5, 9, 13, 17, 21)]
    """

    # The cells_map has both directions, i.e. integer to
    # cell and cell to integer, so the size of T1 is
    # just half of len(cells_map).
    x = (int(len(cells_map)/2) + 1) * [-1]

    for r in range(T1.nrows()):
        for c in range(T1.ncols()):
            e = T1[r, c]

            if e < 0:
                continue

            (r2, c2, e2) = beta1((r, c, e), T1, T2)
            (r3, c3, e3) = beta2((r2, c2, e2), T2, T1)

            x[cells_map[(r, c)]] = cells_map[(r3, c3)]

    x.pop(0)
    # remove the head of the list since we
    # have permutations on 1..(something).

    return Permutation(x)


def back_circulant(n):
    """
    The back-circulant latin square of order n is the Cayley table for
    (Z_n, +), the integers under addition modulo n.

    INPUT:

    - ``n`` -- integer; order of the latin square

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: back_circulant(5)
        [0 1 2 3 4]
        [1 2 3 4 0]
        [2 3 4 0 1]
        [3 4 0 1 2]
        [4 0 1 2 3]
    """
    assert n >= 1

    L = LatinSquare(n, n)

    for r in range(n):
        for c in range(n):
            L[r, c] = (r + c) % n

    return L


def forward_circulant(n):
    """
    The forward-circulant latin square of order n is the Cayley table
    for the operation r + c = (n-c+r) mod n.

    INPUT:

    - ``n`` -- integer; order of the latin square

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: forward_circulant(5)
        [0 4 3 2 1]
        [1 0 4 3 2]
        [2 1 0 4 3]
        [3 2 1 0 4]
        [4 3 2 1 0]
    """
    assert n >= 1

    L = LatinSquare(n, n)

    for r in range(n):
        for c in range(n):
            L[r, c] = (n-c+r) % n

    return L


def direct_product(L1, L2, L3, L4):
    """
    The 'direct product' of four latin squares L1, L2, L3, L4 of order
    n is the latin square of order 2n consisting of

    ::

        -----------
        | L1 | L2 |
        -----------
        | L3 | L4 |
        -----------

    where the subsquares L2 and L3 have entries offset by n.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: direct_product(back_circulant(4), back_circulant(4), elementary_abelian_2group(2), elementary_abelian_2group(2))
        [0 1 2 3 4 5 6 7]
        [1 2 3 0 5 6 7 4]
        [2 3 0 1 6 7 4 5]
        [3 0 1 2 7 4 5 6]
        [4 5 6 7 0 1 2 3]
        [5 4 7 6 1 0 3 2]
        [6 7 4 5 2 3 0 1]
        [7 6 5 4 3 2 1 0]
    """
    assert L1.nrows() == L2.nrows() == L3.nrows() == L4.nrows()
    assert L1.ncols() == L2.ncols() == L3.ncols() == L4.ncols()
    assert L1.nrows() == L1.ncols()

    n = L1.nrows()

    D = LatinSquare(2*n, 2*n)

    for r in range(n):
        for c in range(n):
            D[r, c] = L1[r, c]
            D[r, c+n] = L2[r, c] + n
            D[r+n, c] = L3[r, c] + n
            D[r+n, c+n] = L4[r, c]

    return D


def elementary_abelian_2group(s):
    """
    Return the latin square based on the Cayley table for the
    elementary abelian 2-group of order 2s.

    INPUT:

    - ``s`` -- integer; order of the latin square will be 2s

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: elementary_abelian_2group(3)
        [0 1 2 3 4 5 6 7]
        [1 0 3 2 5 4 7 6]
        [2 3 0 1 6 7 4 5]
        [3 2 1 0 7 6 5 4]
        [4 5 6 7 0 1 2 3]
        [5 4 7 6 1 0 3 2]
        [6 7 4 5 2 3 0 1]
        [7 6 5 4 3 2 1 0]
    """
    assert s > 0

    if s == 1:
        L = LatinSquare(2, 2)
        L[0, 0] = 0
        L[0, 1] = 1
        L[1, 0] = 1
        L[1, 1] = 0

        return L
    else:
        L_prev = elementary_abelian_2group(s-1)
        L = LatinSquare(2**s, 2**s)

        offset = L.nrows() // 2

        for r in range(L_prev.nrows()):
            for c in range(L_prev.ncols()):
                L[r, c] = L_prev[r, c]
                L[r+offset, c] = L_prev[r, c] + offset
                L[r, c+offset] = L_prev[r, c] + offset
                L[r+offset, c+offset] = L_prev[r, c]
    return L


def coin():
    """
    Simulate a fair coin (returns ``True`` or ``False``) using
    ZZ.random_element(2).

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import coin
        sage: x = coin()
        sage: x == 0 or x == 1
        True
    """
    return ZZ.random_element(2) == 0


def next_conjugate(L):
    """
    Permute L[r, c] = e to the conjugate L[c, e] = r.

    We assume that L is an n by n matrix and has values in the range 0,
    1, ..., n-1.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: L = back_circulant(6)
        sage: L
        [0 1 2 3 4 5]
        [1 2 3 4 5 0]
        [2 3 4 5 0 1]
        [3 4 5 0 1 2]
        [4 5 0 1 2 3]
        [5 0 1 2 3 4]
        sage: next_conjugate(L)
        [0 1 2 3 4 5]
        [5 0 1 2 3 4]
        [4 5 0 1 2 3]
        [3 4 5 0 1 2]
        [2 3 4 5 0 1]
        [1 2 3 4 5 0]
        sage: L == next_conjugate(next_conjugate(next_conjugate(L)))
        True
    """
    assert L.nrows() == L.ncols()

    n = L.nrows()

    C = LatinSquare(n, n)

    for r in range(n):
        for c in range(n):
            e = L[r, c]
            assert e >= 0 and e < n

            C[c, e] = r

    return C


def row_containing_sym(L, c, x):
    """
    Given an improper latin square L with L[r1, c] = L[r2, c] = x,
    return r1 or r2 with equal probability. This is an internal
    function and should only be used in LatinSquare_generator().

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: L = matrix([(0, 1, 0, 3), (3, 0, 2, 1), (1, 0, 3, 2), (2, 3, 1, 0)])
        sage: L
        [0 1 0 3]
        [3 0 2 1]
        [1 0 3 2]
        [2 3 1 0]
        sage: c = row_containing_sym(L, 1, 0)
        sage: c == 1 or c == 2
        True
    """
    r1 = -1
    r2 = -1

    for r in range(L.nrows()):
        if r1 >= 0 and r2 >= 0:
            break

        if L[r, c] == x and r1 < 0:
            r1 = r
            continue

        if L[r, c] == x and r2 < 0:
            r2 = r
            break

    assert r1 >= 0 and r2 >= 0

    return r1 if coin() else r2


def column_containing_sym(L, r, x):
    """
    Given an improper latin square L with L[r, c1] = L[r, c2] = x,
    return c1 or c2 with equal probability. This is an internal
    function and should only be used in LatinSquare_generator().

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: L = matrix([(1, 0, 2, 3), (0, 2, 3, 0), (2, 3, 0, 1), (3, 0, 1, 2)])
        sage: L
        [1 0 2 3]
        [0 2 3 0]
        [2 3 0 1]
        [3 0 1 2]
        sage: c = column_containing_sym(L, 1, 0)
        sage: c == 0 or c == 3
        True
    """

    c1 = -1
    c2 = -1

    for c in range(L.ncols()):
        if c1 >= 0 and c2 >= 0:
            break

        if L[r, c] == x and c1 < 0:
            c1 = c
            continue

        if L[r, c] == x and c2 < 0:
            c2 = c
            break

    assert c1 >= 0 and c2 >= 0

    return c1 if coin() else c2


def LatinSquare_generator(L_start, check_assertions=False):
    """
    Generator for a sequence of uniformly distributed latin squares,
    given L_start as the initial latin square.

    This code implements
    the Markov chain algorithm of Jacobson and Matthews (1996), see
    below for the BibTex entry. This generator will never throw the
    ``StopIteration`` exception, so it provides an infinite sequence of
    latin squares.

    EXAMPLES:

    Use the back circulant latin square of order 4 as the initial
    square and print the next two latin squares given by the Markov
    chain::

        sage: from sage.combinat.matrices.latin import *
        sage: g = LatinSquare_generator(back_circulant(4))
        sage: next(g).is_latin_square()
        True

    REFERENCES:

    .. [JacMat96] Mark T. Jacobson and Peter Matthews, "Generating uniformly
       distributed random Latin squares", Journal of Combinatorial Designs,
       4 (1996)
    """
    if check_assertions:
        assert L_start.is_latin_square()

    n = L_start.nrows()

    r1 = r2 = c1 = c2 = x = y = z = -1
    proper = True

    from copy import copy
    L = copy(L_start)

    L_cer = LatinSquare(n, n)
    L_erc = LatinSquare(n, n)

    while True:
        if proper:
            if check_assertions:
                assert L.is_latin_square()

            #################################
            # Update the other two conjugates
            for r in range(n):
                for c in range(n):
                    e = L[r, c]

                    L_cer[c, e] = r
                    L_erc[e, r] = c
            #################################

            yield L

            r1 = ZZ.random_element(n)
            c1 = ZZ.random_element(n)
            x = L[r1, c1]

            y = x
            while y == x:
                y = ZZ.random_element(n)

            # Now find y in row r1 and column c1.
            if check_assertions:
                r2 = 0
                c2 = 0
                while L[r1, c2] != y:
                    c2 += 1
                while L[r2, c1] != y:
                    r2 += 1

                assert L_erc[y, r1] == c2
                assert L_cer[c1, y] == r2

            c2 = L_erc[y, r1]
            r2 = L_cer[c1, y]

            if check_assertions:
                assert L[r1, c2] == y
            if check_assertions:
                assert L[r2, c1] == y

            L[r1, c1] = y
            L[r1, c2] = x
            L[r2, c1] = x

            # Now deal with the unknown point.
            # We want to form z + (y - x)
            z = L[r2, c2]

            if z == x:
                L[r2, c2] = y
            else:
                # z and y have positive coefficients
                # x is the improper term with a negative coefficient
                proper = False
        else:  # improper square,
            # L[r2, c2] = y + z - x
            # y and z are proper while x is the
            # improper symbol in the cell L[r2, c2].

            r1 = row_containing_sym(L, c2, x)
            c1 = column_containing_sym(L, r2, x)

            if check_assertions:
                assert L[r1, c2] == x
            if check_assertions:
                assert L[r2, c1] == x

            # choose one of the proper symbols
            # uniformly at random (we will use whatever
            # lands in variable y).
            if coin():
                y, z = z, y

            # Add/subtract the symbolic difference (y - x)
            L[r2, c2] = z
            L[r1, c2] = y
            L[r2, c1] = y

            if L[r1, c1] == y:
                L[r1, c1] = x
                proper = True
            else:  # got another improper square
                z = L[r1, c1]
                x, y = y, x
                r2 = r1
                c2 = c1

                # Now we have L[r2, c2] = z+y-x as
                # usual
                proper = False  # for emphasis


def group_to_LatinSquare(G):
    """
    Construct a latin square on the symbols [0, 1, ..., n-1] for a
    group with an n by n Cayley table.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import group_to_LatinSquare

        sage: group_to_LatinSquare(DihedralGroup(2))
        [0 1 2 3]
        [1 0 3 2]
        [2 3 0 1]
        [3 2 1 0]

    ::

        sage: G = libgap.Group(PermutationGroupElement((1,2,3)))
        sage: group_to_LatinSquare(G)
        [0 1 2]
        [1 2 0]
        [2 0 1]
    """
    if isinstance(G, GapElement):
        rows = (list(x) for x in list(libgap.MultiplicationTable(G)))
        new_rows = []

        for x in rows:
            new_rows.append([int(xx) - 1 for xx in x])

        return matrix(new_rows)

    # Otherwise we must have some kind of Sage permutation group object,
    # such as sage.groups.perm_gps.permgroup.PermutationGroup_generic
    # or maybe sage.groups.perm_gps.permgroup_named.

    T = G.cayley_table()
    return matrix(ZZ, T.table())


def alternating_group_bitrade_generators(m):
    r"""
    Construct generators a, b, c for the alternating group on 3m+1
    points, such that a\*b\*c = 1.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: a, b, c, G = alternating_group_bitrade_generators(1)
        sage: (a, b, c, G)
        ((1,2,3), (1,4,2), (2,4,3), Permutation Group with generators [(1,2,3), (1,4,2)])
        sage: a*b*c
        ()

    ::

        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: T1
        [ 0 -1  3  1]
        [-1  1  0  2]
        [ 1  3  2 -1]
        [ 2  0 -1  3]
        sage: T2
        [ 1 -1  0  3]
        [-1  0  2  1]
        [ 2  1  3 -1]
        [ 0  3 -1  2]
    """
    assert m >= 1

    a = tuple(range(1, 2*m+1 + 1))

    b = tuple(range(m + 1, 0, -1)) + tuple(range(2*m+2, 3*m+1 + 1))

    a = PermutationConstructor(a)
    b = PermutationConstructor(b)
    c = PermutationConstructor((a*b)**(-1))

    G = PermutationGroup([a, b])

    return (a, b, c, G)


def pq_group_bitrade_generators(p, q):
    """
    Generators for a group of order pq where p and q are primes such
    that (q % p) == 1.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: pq_group_bitrade_generators(3,7)
        ((2,3,5)(4,7,6), (1,2,3,4,5,6,7), (1,4,2)(3,5,6), Permutation Group with generators [(2,3,5)(4,7,6), (1,2,3,4,5,6,7)])
    """
    assert is_prime(p)
    assert is_prime(q)
    assert (q % p) == 1

    # beta is a primitive root of the
    # congruence x^p = 1 mod q
    F = FiniteField(q)
    fgen = F.multiplicative_generator()
    beta = fgen**((q-1)/p)

    assert beta != 1
    assert (beta**p % q) == 1

    Q = tuple(range(1, q+1))

    P = []
    seenValues = {}
    for i in range(2, q):
        if i in seenValues:
            continue

        cycle = []
        for k in range(p):
            x = (1 + (i-1)*beta**k) % q
            if x == 0:
                x = q

            seenValues[x] = True
            cycle.append(x)
        P.append(tuple(map(Integer, cycle)))

    G = PermutationGroup([P, Q])
    assert G.order() == p*q
    assert not G.is_abelian()

    a = PermutationConstructor(P)
    b = PermutationConstructor(Q)
    c = PermutationConstructor((a*b)**(-1))

    return (a, b, c, PermutationGroup([P, Q]))


def p3_group_bitrade_generators(p):
    """
    Generators for a group of order p3 where p is a prime.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: p3_group_bitrade_generators(3)
        ((2,6,7)(3,8,9),
         (1,2,3)(4,7,8)(5,6,9),
         (1,9,2)(3,7,4)(5,8,6),
         Permutation Group with generators [(2,6,7)(3,8,9), (1,2,3)(4,7,8)(5,6,9)])
    """
    assert is_prime(p)

    F = libgap.FreeGroup(3)
    a, b, c = F.GeneratorsOfGroup()

    rels = []
    rels.append(a**p)
    rels.append(b**p)
    rels.append(c**p)
    rels.append(a*b*((b*a*c)**(-1)))
    rels.append(c*a*((a*c)**(-1)))
    rels.append(c*b*((b*c)**(-1)))

    G = F.FactorGroupFpGroupByRels(rels)
    u, v, _ = G.GeneratorsOfGroup()

    iso = libgap.IsomorphismPermGroup(G)

    x = PermutationConstructor(libgap.Image(iso, u))
    y = PermutationConstructor(libgap.Image(iso, v))

    return (x, y, (x*y)**(-1), PermutationGroup([x, y]))


def check_bitrade_generators(a, b, c):
    r"""
    Three group elements a, b, c will generate a bitrade if a\*b\*c = 1
    and the subgroups a, b, c intersect (pairwise) in just the
    identity.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: a, b, c, G = p3_group_bitrade_generators(3)
        sage: check_bitrade_generators(a, b, c)
        True
        sage: check_bitrade_generators(a, b, libgap(gap('()')))
        False
    """
    A = PermutationGroup([a])
    B = PermutationGroup([b])
    C = PermutationGroup([c])

    if a*b != c**(-1):
        return False

    X = libgap.Intersection(libgap.Intersection(A, B), C)
    return X.Size() == 1


def is_bitrade(T1, T2):
    """
    Combinatorially, a pair (T1, T2) of partial latin squares is a
    bitrade if they are disjoint, have the same shape, and have row and
    column balance. For definitions of each of these terms see the
    relevant function in this file.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: is_bitrade(T1, T2)
        True
    """
    return (is_disjoint(T1, T2) and is_same_shape(T1, T2) and
            is_row_and_col_balanced(T1, T2))


def is_primary_bitrade(a, b, c, G):
    """
    A bitrade generated from elements a, b, c is primary if a, b, c =
    G.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: (a, b, c, G) = p3_group_bitrade_generators(5)
        sage: is_primary_bitrade(a, b, c, G)
        True
    """
    return G == PermutationGroup([a, b, c])


def tau_to_bitrade(t1, t2, t3):
    """
    Given permutations t1, t2, t3 that represent a latin bitrade,
    convert them to an explicit latin bitrade (T1, T2). The result is
    unique up to isotopism.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: _, t1, t2, t3 = tau123(T1, T2)
        sage: U1, U2 = tau_to_bitrade(t1, t2, t3)
        sage: assert is_bitrade(U1, U2)
        sage: U1
        [0 1 2 3 4]
        [1 2 3 4 0]
        [2 3 4 0 1]
        [3 4 0 1 2]
        [4 0 1 2 3]
        sage: U2
        [4 0 1 2 3]
        [0 1 2 3 4]
        [1 2 3 4 0]
        [2 3 4 0 1]
        [3 4 0 1 2]
    """
    c1 = t1.to_cycles()
    c2 = t2.to_cycles()
    c3 = t3.to_cycles()

    pt_to_cycle1 = {}
    pt_to_cycle2 = {}
    pt_to_cycle3 = {}

    for i in range(len(c1)):
        for j in range(len(c1[i])):
            pt_to_cycle1[c1[i][j]] = i

    for i in range(len(c2)):
        for j in range(len(c2[i])):
            pt_to_cycle2[c2[i][j]] = i

    for i in range(len(c3)):
        for j in range(len(c3[i])):
            pt_to_cycle3[c3[i][j]] = i

    n = max(len(c1), len(c2), len(c3))

    T1 = LatinSquare(n)
    T2 = LatinSquare(n)

    for r in range(len(c1)):
        for c in range(len(c2)):
            for s in range(len(c3)):
                nr_common = len(reduce(set.intersection,
                                       [set(c1[r]), set(c2[c]), set(c3[s])]))
                assert nr_common in [0, 1]

                if nr_common == 1:
                    T1[r, c] = s

    for cycle in c1:
        for pt1 in cycle:
            pt2 = t1[pt1 - 1]
            pt3 = t2[pt2 - 1]
            assert t3[pt3 - 1] == pt1

            r = pt_to_cycle1[pt1]
            c = pt_to_cycle2[pt2]
            s = pt_to_cycle3[pt3]

            T2[r, c] = s

    return T1, T2


def bitrade_from_group(a, b, c, G):
    """
    Given group elements a, b, c in G such that abc = 1 and the
    subgroups a, b, c intersect (pairwise) only in the identity,
    construct a bitrade (T1, T2) where rows, columns, and symbols
    correspond to cosets of a, b, and c, respectively.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: a, b, c, G = alternating_group_bitrade_generators(1)
        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: T1
        [ 0 -1  3  1]
        [-1  1  0  2]
        [ 1  3  2 -1]
        [ 2  0 -1  3]
        sage: T2
        [ 1 -1  0  3]
        [-1  0  2  1]
        [ 2  1  3 -1]
        [ 0  3 -1  2]
    """
    hom = libgap.ActionHomomorphism(G, libgap.RightCosets(G, libgap.TrivialSubgroup(G)), libgap.OnRight)

    t1 = libgap.Image(hom, a)
    t2 = libgap.Image(hom, b)
    t3 = libgap.Image(hom, c)

    t1 = Permutation(str(t1).replace('\n', ''))
    t2 = Permutation(str(t2).replace('\n', ''))
    t3 = Permutation(str(t3).replace('\n', ''))

    return tau_to_bitrade(t1, t2, t3)


def is_disjoint(T1, T2):
    """
    The partial latin squares T1 and T2 are disjoint if T1[r, c] !=
    T2[r, c] or T1[r, c] == T2[r, c] == -1 for each cell [r, c].

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import is_disjoint, back_circulant, isotopism
        sage: is_disjoint(back_circulant(2), back_circulant(2))
        False

    ::

        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: is_disjoint(T1, T2)
        True
    """
    for i in range(T1.nrows()):
        for j in range(T1.ncols()):
            if T1[i, j] < 0 and T2[i, j] < 0:
                continue
            if T1[i, j] == T2[i, j]:
                return False

    return True


def is_same_shape(T1, T2):
    """
    Two partial latin squares T1, T2 have the same shape if T1[r, c] =
    0 if and only if T2[r, c] = 0.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: is_same_shape(elementary_abelian_2group(2), back_circulant(4))
        True
        sage: is_same_shape(LatinSquare(5), LatinSquare(5))
        True
        sage: is_same_shape(forward_circulant(5), LatinSquare(5))
        False
    """
    for i in range(T1.nrows()):
        for j in range(T1.ncols()):
            if T1[i, j] < 0 and T2[i, j] < 0:
                continue
            if T1[i, j] >= 0 and T2[i, j] >= 0:
                continue

            return False

    return True


def is_row_and_col_balanced(T1, T2):
    """
    Partial latin squares T1 and T2 are balanced if the symbols
    appearing in row r of T1 are the same as the symbols appearing in
    row r of T2, for each r, and if the same condition holds on
    columns.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: T1 = matrix([[0,1,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1]])
        sage: T2 = matrix([[0,1,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1]])
        sage: is_row_and_col_balanced(T1, T2)
        True
        sage: T2 = matrix([[0,3,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1]])
        sage: is_row_and_col_balanced(T1, T2)
        False
    """
    for r in range(T1.nrows()):
        val1 = set(x for x in T1.row(r) if x >= 0)
        val2 = set(x for x in T2.row(r) if x >= 0)
        if val1 != val2:
            return False

    for c in range(T1.ncols()):
        val1 = set(x for x in T1.column(c) if x >= 0)
        val2 = set(x for x in T2.column(c) if x >= 0)
        if val1 != val2:
            return False

    return True


def dlxcpp_rows_and_map(P):
    """
    Internal function for ``dlxcpp_find_completions``. Given a partial
    latin square P we construct a list of rows of a 0-1 matrix M such
    that an exact cover of M corresponds to a completion of P to a
    latin square.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: dlxcpp_rows_and_map(LatinSquare(2))
        ([[0, 4, 8],
          [1, 5, 8],
          [2, 4, 9],
          [3, 5, 9],
          [0, 6, 10],
          [1, 7, 10],
          [2, 6, 11],
          [3, 7, 11]],
         {(0, 4, 8): (0, 0, 0),
          (0, 6, 10): (1, 0, 0),
          (1, 5, 8): (0, 0, 1),
          (1, 7, 10): (1, 0, 1),
          (2, 4, 9): (0, 1, 0),
          (2, 6, 11): (1, 1, 0),
          (3, 5, 9): (0, 1, 1),
          (3, 7, 11): (1, 1, 1)})
    """
    assert P.nrows() == P.ncols()

    n = P.nrows()

    # We will need 3n^2 columns in total:
    #
    # n^2 for the xCy columns
    # n^2 for the xRy columns
    # n^2 for the xy columns

    dlx_rows = []
    cmap = {}

    for r in range(n):
        valsrow = P.vals_in_row(r)

        for c in range(n):
            valscol = P.vals_in_col(c)

            for e in range(n):
                # These should be constants
                c_OFFSET = e + c*n
                r_OFFSET = e + r*n + n*n
                xy_OFFSET = 2*n*n + r*n + c

                cmap[(c_OFFSET, r_OFFSET, xy_OFFSET)] = (r, c, e)

                # We only want the correct value to pop in here
                if P[r, c] >= 0 and P[r, c] != e:
                    continue
                if P[r, c] < 0 and e in valsrow:
                    continue
                if P[r, c] < 0 and e in valscol:
                    continue

                dlx_rows.append([c_OFFSET, r_OFFSET, xy_OFFSET])

    return dlx_rows, cmap


def dlxcpp_find_completions(P, nr_to_find=None):
    """
    Return a list of all latin squares L of the same order as P such
    that P is contained in L. The optional parameter nr_to_find
    limits the number of latin squares that are found.

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: dlxcpp_find_completions(LatinSquare(2))
        [[0 1]
        [1 0], [1 0]
        [0 1]]

    ::

        sage: dlxcpp_find_completions(LatinSquare(2), 1)
        [[0 1]
        [1 0]]
    """
    assert P.nrows() == P.ncols()

    dlx_rows, cmap = dlxcpp_rows_and_map(P)

    SOLUTIONS = {}
    for x in DLXCPP(dlx_rows):
        x.sort()
        SOLUTIONS[tuple(x)] = True

        if nr_to_find is not None and len(SOLUTIONS) >= nr_to_find:
            break

    comps = []

    for i in SOLUTIONS:
        soln = list(i)

        from copy import deepcopy
        Q = deepcopy(P)

        for x in soln:
            (r, c, e) = cmap[tuple(dlx_rows[x])]

            if Q[r, c] >= 0:
                assert Q[r, c] == e
            else:
                Q[r, c] = e

        comps.append(Q)

    return comps


def bitrade(T1, T2):
    r"""
    Form the bitrade (Q1, Q2) from (T1, T2) by setting empty the cells
    (r, c) such that T1[r, c] == T2[r, c].

    EXAMPLES::

        sage: from sage.combinat.matrices.latin import *
        sage: B1 = back_circulant(5)
        sage: alpha = isotopism((0,1,2,3,4))
        sage: beta  = isotopism((1,0,2,3,4))
        sage: gamma = isotopism((2,1,0,3,4))
        sage: B2 = B1.apply_isotopism(alpha, beta, gamma)
        sage: T1, T2 = bitrade(B1, B2)
        sage: T1
        [ 0  1 -1  3  4]
        [ 1 -1 -1  4  0]
        [ 2 -1  4  0  1]
        [ 3  4  0  1  2]
        [ 4  0  1  2  3]
        sage: T2
        [ 3  4 -1  0  1]
        [ 0 -1 -1  1  4]
        [ 1 -1  0  4  2]
        [ 4  0  1  2  3]
        [ 2  1  4  3  0]
    """
    assert T1.nrows() == T1.ncols()
    assert T2.nrows() == T2.ncols()
    assert T1.nrows() == T2.nrows()

    n = T1.nrows()

    from copy import copy
    Q1 = copy(T1)
    Q2 = copy(T2)

    for r in range(n):
        for c in range(n):
            if T1[r, c] == T2[r, c]:
                Q1[r, c] = -1
                Q2[r, c] = -1

    return Q1, Q2
