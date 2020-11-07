r"""
Sudoku Puzzles

This module provides algorithms to solve Sudoku puzzles, plus tools
for inputting, converting and displaying various ways of writing a
puzzle or its solution(s).  Primarily this is accomplished with the
:class:`sage.games.sudoku.Sudoku` class, though the legacy top-level
:func:`sage.games.sudoku.sudoku` function is also available.

AUTHORS:

- Tom Boothby (2008/05/02): Exact Cover, Dancing Links algorithm
- Robert Beezer (2009/05/29): Backtracking algorithm, Sudoku class
"""
######################################################################
#       Copyright (C) 2009, Robert A. Beezer <beezer@ups.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
######################################################################

from sage.structure.sage_object import SageObject


def sudoku(m):
    r"""
    Solves Sudoku puzzles described by matrices.

    INPUT:

    - ``m`` - a square Sage matrix over `\ZZ`, where zeros are blank entries

    OUTPUT:

    A Sage matrix over `\ZZ` containing the first solution found,
    otherwise ``None``.

    This function matches the behavior of the prior Sudoku solver
    and is included only to replicate that behavior.  It could be
    safely deprecated, since all of its functionality is included in the :class:`~sage.games.sudoku.Sudoku` class.

    EXAMPLES:

    An example that was used in previous doctests. ::

        sage: A = matrix(ZZ,9,[5,0,0, 0,8,0, 0,4,9, 0,0,0, 5,0,0, 0,3,0, 0,6,7, 3,0,0, 0,0,1, 1,5,0, 0,0,0, 0,0,0,  0,0,0, 2,0,8, 0,0,0, 0,0,0, 0,0,0, 0,1,8, 7,0,0, 0,0,4, 1,5,0, 0,3,0, 0,0,2, 0,0,0, 4,9,0, 0,5,0, 0,0,3])
        sage: A
        [5 0 0 0 8 0 0 4 9]
        [0 0 0 5 0 0 0 3 0]
        [0 6 7 3 0 0 0 0 1]
        [1 5 0 0 0 0 0 0 0]
        [0 0 0 2 0 8 0 0 0]
        [0 0 0 0 0 0 0 1 8]
        [7 0 0 0 0 4 1 5 0]
        [0 3 0 0 0 2 0 0 0]
        [4 9 0 0 5 0 0 0 3]
        sage: sudoku(A)
        [5 1 3 6 8 7 2 4 9]
        [8 4 9 5 2 1 6 3 7]
        [2 6 7 3 4 9 5 8 1]
        [1 5 8 4 6 3 9 7 2]
        [9 7 4 2 1 8 3 6 5]
        [3 2 6 7 9 5 4 1 8]
        [7 8 2 9 3 4 1 5 6]
        [6 3 5 1 7 2 8 9 4]
        [4 9 1 8 5 6 7 2 3]

    Using inputs that are possible with the
    :class:`~sage.games.sudoku.Sudoku` class,
    other than a matrix, will cause an error. ::

        sage: sudoku('.4..32....14..3.')
        Traceback (most recent call last):
        ...
        ValueError: sudoku function expects puzzle to be a matrix, perhaps use the Sudoku class
    """
    from sage.structure.element import is_Matrix

    if not is_Matrix(m):
        raise ValueError('sudoku function expects puzzle to be a matrix, perhaps use the Sudoku class')
    solution = next(Sudoku(m).solve(algorithm='dlx'))
    return (solution.to_matrix() if solution else None)


class Sudoku(SageObject):
    r"""
    An object representing a Sudoku puzzle.  Primarily the purpose is to
    solve the puzzle, but conversions between formats are also provided.

    INPUT:

    - puzzle -- the first argument can take one of three forms
        * list - a Python list with elements of the puzzle in row-major order,
          where a blank entry is a zero
        * matrix - a square Sage matrix over `\ZZ`
        * string - a string where each character is an entry of
          the puzzle. For two-digit entries, a = 10, b = 11, etc.
    - verify_input -- default = ``True``, use ``False`` if you know the input is valid

    EXAMPLES::

        sage: a = Sudoku('5...8..49...5...3..673....115..........2.8..........187....415..3...2...49..5...3')
        sage: print(a)
        +-----+-----+-----+
        |5    |  8  |  4 9|
        |     |5    |  3  |
        |  6 7|3    |    1|
        +-----+-----+-----+
        |1 5  |     |     |
        |     |2   8|     |
        |     |     |  1 8|
        +-----+-----+-----+
        |7    |    4|1 5  |
        |  3  |    2|     |
        |4 9  |  5  |    3|
        +-----+-----+-----+
        sage: print(next(a.solve()))
        +-----+-----+-----+
        |5 1 3|6 8 7|2 4 9|
        |8 4 9|5 2 1|6 3 7|
        |2 6 7|3 4 9|5 8 1|
        +-----+-----+-----+
        |1 5 8|4 6 3|9 7 2|
        |9 7 4|2 1 8|3 6 5|
        |3 2 6|7 9 5|4 1 8|
        +-----+-----+-----+
        |7 8 2|9 3 4|1 5 6|
        |6 3 5|1 7 2|8 9 4|
        |4 9 1|8 5 6|7 2 3|
        +-----+-----+-----+
    """
    def __init__(self, puzzle, verify_input = True):
        r"""
        Initialize a Sudoku puzzle, determine its size, sanity-check the inputs.

        TESTS::

            sage: d = Sudoku('1.......2.9.4...5...6...7...5.9.3.......7.......85..4.7.....6...3...9.8...2.....1')
            sage: d == loads(dumps(d))
            True

        A lame attempt to construct a puzzle from a single integer::

            sage: Sudoku(8)
            Traceback (most recent call last):
            ...
            ValueError: Sudoku puzzle must be specified as a matrix, list or string

        An attempt to construct a puzzle from a non-square matrix::

            sage: Sudoku(matrix(2,range(6)))
            Traceback (most recent call last):
            ...
            ValueError: Sudoku puzzle must be a square matrix

        An attempt to construct a puzzle from a string of an impossible length (9)::

            sage: Sudoku('.........')
            Traceback (most recent call last):
            ...
            ValueError: Sudoku puzzle dimension of 3 must be a perfect square

        An attempt to construct a `4\times 4` puzzle from a string with a bad entry (5)::

            sage: Sudoku('.1.2.......5....')
            Traceback (most recent call last):
            ...
            ValueError: Sudoku puzzle has an invalid entry
        """
        from math import sqrt
        from sage.structure.element import is_Matrix

        if isinstance(puzzle, list):
            puzzle_size = int(round(sqrt(len(puzzle))))
            self.puzzle = tuple(puzzle)
        elif is_Matrix(puzzle):
            puzzle_size = puzzle.ncols()
            if verify_input and not(puzzle.is_square()):
                raise ValueError('Sudoku puzzle must be a square matrix')
            self.puzzle = tuple([int(x) for x in puzzle.list()])
        elif isinstance(puzzle, str):
            puzzle_size = int(round(sqrt(len(puzzle))))
            puzzle_numeric = []
            for char in puzzle:
                if char.isdigit():
                    puzzle_numeric.append(int(char))
                elif char == '.':
                    puzzle_numeric.append(0)
                else:
                    puzzle_numeric.append(ord(char.upper()) - ord('A')+10)
            self.puzzle = tuple(puzzle_numeric)
        else:
            raise ValueError('Sudoku puzzle must be specified as a matrix, list or string')
        self.n = int(sqrt(puzzle_size))
        if verify_input:
            if self.n**4 != len(self.puzzle):
                raise ValueError('Sudoku puzzle dimension of %s must be a perfect square' % puzzle_size)
            for x in self.puzzle:
                if (x < 0) or (x > self.n*self.n):
                    raise ValueError('Sudoku puzzle has an invalid entry')

    def __eq__(self, other):
        r"""
        Compares two Sudoku puzzles, based on the underlying
        representation of the puzzles as tuples.

        EXAMPLES::

            sage: a = Sudoku('.4..32....14..3.')
            sage: b = Sudoku('8..6..9.5.............2.31...7318.6.24.....73...........279.1..5...8..36..3......')
            sage: c = Sudoku('1..6..9.5.............2.31...7318.6.24.....73...........279.1..5...8..36..3......')
            sage: d = Sudoku('81.6..9.5.............2.31...7318.6.24.....73...........279.1..5...8..36..3......')

            sage: a == b
            False
            sage: b == b
            True
            sage: b == c
            False
            sage: b == d
            False
        """
        return self.puzzle == tuple(other.to_list())

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: a = Sudoku('.4..32....14..3.')
            sage: hash(a) == hash(a.puzzle)
            True
        """
        return hash(self.puzzle)

    def __ne__(self, other):
        """
        Check that ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: a = Sudoku('.4..32....14..3.')
            sage: b = Sudoku('8..6..9.5.............2.31...7318.6.24.....73...........279.1..5...8..36..3......')
            sage: c = Sudoku('1..6..9.5.............2.31...7318.6.24.....73...........279.1..5...8..36..3......')
            sage: d = Sudoku('81.6..9.5.............2.31...7318.6.24.....73...........279.1..5...8..36..3......')

            sage: a != b
            True
            sage: b != b
            False
            sage: b != c
            True
            sage: b != d
            True
        """
        return not (self == other)

    def _repr_(self):
        r"""
        Return a concise description of a Sudoku puzzle using a string representation.

        See the docstring for :func:`to_ascii` for more information on the format.

        EXAMPLES::

            sage: s = Sudoku('.4..32....14..3.')
            sage: s._repr_()
            '+---+---+\n|  4|   |\n|3 2|   |\n+---+---+\n|   |1 4|\n|   |3  |\n+---+---+'

        """
        return self.to_ascii()

    def _latex_(self):
        r"""nodetex
        Return a `\LaTeX` representation of a Sudoku puzzle as an array environment.

        EXAMPLES::

            sage: s = Sudoku('.4..32....14..3.')
            sage: s._latex_()
            '\\begin{array}{|*{2}{*{2}{r}|}}\\hline\n &4& & \\\\\n3&2& & \\\\\\hline\n & &1&4\\\\\n & &3& \\\\\\hline\n\\end{array}'
        """
        return self.to_latex()

    def _matrix_(self, R=None):
        r"""
        Return the puzzle as a matrix to support Sage's
        :func:`~sage.matrix.constructor.matrix` constructor.

        The base ring will be `\ZZ` if ``None`` is provided,
        and it is an error to specify any other base ring.

        EXAMPLES::

            sage: k = Sudoku('.4..32....14..3.')
            sage: matrix(k) # indirect doctest
            [0 4 0 0]
            [3 2 0 0]
            [0 0 1 4]
            [0 0 3 0]
            sage: matrix(ZZ,k)
            [0 4 0 0]
            [3 2 0 0]
            [0 0 1 4]
            [0 0 3 0]
            sage: matrix(QQ,k)
            Traceback (most recent call last):
            ...
            ValueError: Sudoku puzzles only convert to matrices over Integer Ring, not Rational Field
        """
        from sage.rings.integer_ring import ZZ, IntegerRing_class
        if R and not(isinstance(R, IntegerRing_class)):
            raise ValueError('Sudoku puzzles only convert to matrices over %s, not %s' % (ZZ, R))
        return self.to_matrix()

    def to_string(self):
        r"""
        Construct a string representing a Sudoku puzzle.

        Blank entries are represented as periods, single
        digits are not converted and two digit entries are
        converted to lower-case letters where ``10 = a``,
        ``11 = b``, etc.  This scheme limits puzzles to
        at most 36 symbols.

        EXAMPLES::

            sage: b = matrix(ZZ, 9, 9, [ [0,0,0,0,1,0,9,0,0], [8,0,0,4,0,0,0,0,0], [2,0,0,0,0,0,0,0,0], [0,7,0,0,3,0,0,0,0], [0,0,0,0,0,0,2,0,4], [0,0,0,0,0,0,0,5,8], [0,6,0,0,0,0,1,3,0], [7,0,0,2,0,0,0,0,0], [0,0,0,8,0,0,0,0,0] ])
            sage: Sudoku(b).to_string()
            '....1.9..8..4.....2.........7..3..........2.4.......58.6....13.7..2........8.....'

        TESTS:

        This tests the conversion of alphabetic characters as well as the
        input and output of Sudoku puzzles as strings. ::

            sage: j = Sudoku([0, 0, 0, 0, 10, 0, 0, 6, 9, 0, 3, 0, 0, 0, 0, 1, 13, 0, 2, 0, 0, 0, 8, 0, 0, 0, 0, 14, 0, 4, 0, 0, 0, 0, 11, 0, 0, 0, 0, 5, 0, 0, 12, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, 15, 0, 0, 0, 0, 1, 0, 14, 0, 0, 2, 0, 11, 0, 8, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 13, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15, 0, 0, 7, 0, 16, 0, 0, 9, 0, 10, 0, 0, 12, 0, 0, 0, 5, 0, 0, 0, 0, 0, 8, 0, 0, 15, 0, 0, 0, 0, 0, 1, 0, 0, 14, 0, 7, 9, 0, 12, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0, 2, 0, 0, 0, 6, 0, 0, 0, 0, 0, 16, 0, 7, 0, 0, 0, 0, 0, 0, 8, 4, 0, 0, 0, 0, 3, 0, 13, 0, 0, 10, 0, 5, 0, 0, 0, 0, 5, 0, 0, 0, 7, 0, 0, 14, 0, 0, 0, 12, 10, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0, 0, 15, 0, 0, 0, 0, 4, 0, 0, 0, 13, 0, 0, 14, 0, 0, 16, 0, 9, 2, 0, 6, 0, 0, 8, 0, 0, 0, 0])
            sage: st = j.to_string()
            sage: st
            '....a..69.3....1d.2...8....e.4....b....5..c.......7.......g...f....1.e..2.b.8..3.......4.d.....6.........f..7.g..9.a..c...5.....8..f.....1..e.79.c....b.....2...6.....g.7......84....3.d..a.5....5...7..e...ca.....3.1.......b......f....4...d..e..g.92.6..8....'
            sage: st == Sudoku(st).to_string()
            True

        A `49\times 49` puzzle with all entries equal to 40,
        which doesn't convert to a letter. ::

            sage: empty = [40]*2401
            sage: Sudoku(empty).to_string()
            Traceback (most recent call last):
            ...
            ValueError: Sudoku string representation is only valid for puzzles of size 36 or smaller
        """
        encoded = []
        for x in self.puzzle:
            if x == 0:
                encoded.append('.')
            elif 1 <= x <= 9:
                encoded.append(str(x))
            elif x <= 36:
                encoded.append(chr(x-10+ord('a')))
            else:
                raise ValueError('Sudoku string representation is only valid for puzzles of size 36 or smaller')
        return ''.join(encoded)

    def to_list(self):
        r"""
        Construct a list representing a Sudoku puzzle, in row-major order.

        EXAMPLES::

            sage: s = Sudoku('1.......2.9.4...5...6...7...5.9.3.......7.......85..4.7.....6...3...9.8...2.....1')
            sage: s.to_list()
            [1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 9, 0, 4, 0, 0, 0, 5, 0, 0, 0, 6, 0, 0, 0, 7, 0, 0, 0, 5, 0, 9, 0, 3, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 8, 5, 0, 0, 4, 0, 7, 0, 0, 0, 0, 0, 6, 0, 0, 0, 3, 0, 0, 0, 9, 0, 8, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1]

        TESTS:

        This tests the input and output of Sudoku puzzles as lists. ::

            sage: alist = [0, 4, 0, 0, 3, 2, 0, 0, 0, 0, 1, 4, 0, 0, 3, 0]
            sage: alist == Sudoku(alist).to_list()
            True
        """
        return list(self.puzzle)

    def to_matrix(self):
        r"""
        Construct a Sage matrix over `\ZZ` representing a Sudoku puzzle.

        EXAMPLES::

            sage: s = Sudoku('.4..32....14..3.')
            sage: s.to_matrix()
            [0 4 0 0]
            [3 2 0 0]
            [0 0 1 4]
            [0 0 3 0]

        TESTS:

        This tests the input and output of Sudoku puzzles as matrices over `\ZZ`. ::

            sage: g = matrix(ZZ, 9, 9, [ [1,0,0,0,0,7,0,9,0], [0,3,0,0,2,0,0,0,8], [0,0,9,6,0,0,5,0,0], [0,0,5,3,0,0,9,0,0], [0,1,0,0,8,0,0,0,2], [6,0,0,0,0,4,0,0,0], [3,0,0,0,0,0,0,1,0], [0,4,0,0,0,0,0,0,7], [0,0,7,0,0,0,3,0,0] ])
            sage: g == Sudoku(g).to_matrix()
            True
        """
        from sage.rings.integer_ring import ZZ
        from sage.matrix.constructor import matrix
        return matrix(ZZ, self.n*self.n, self.puzzle)


    def to_ascii(self):
        r"""
        Construct an ASCII-art version of a Sudoku puzzle.
        This is a modified version of the ASCII version of a subdivided matrix.

        EXAMPLES::

            sage: s = Sudoku('.4..32....14..3.')
            sage: print(s.to_ascii())
            +---+---+
            |  4|   |
            |3 2|   |
            +---+---+
            |   |1 4|
            |   |3  |
            +---+---+
            sage: s.to_ascii()
            '+---+---+\n|  4|   |\n|3 2|   |\n+---+---+\n|   |1 4|\n|   |3  |\n+---+---+'
        """
        from re import compile
        n = self.n
        nsquare = n*n
        m = self.to_matrix()
        m.subdivide(list(range(0,nsquare+1,n)), list(range(0,nsquare+1,n)))
        naked_zero = compile(r'([\|, ]+)0')
        blanked = naked_zero.sub(lambda x: x.group(1)+' ', m.str())
        brackets = compile(r'[\[,\]]')
        return brackets.sub('', blanked)

    def to_latex(self):
        r"""
        Create a string of `\LaTeX` code representing a Sudoku puzzle or solution.

        EXAMPLES::

            sage: s = Sudoku('.4..32....14..3.')
            sage: print(s.to_latex())
            \begin{array}{|*{2}{*{2}{r}|}}\hline
            &4& & \\
            3&2& & \\\hline
            & &1&4\\
            & &3& \\\hline
            \end{array}

        TESTS::

            sage: s = Sudoku('.4..32....14..3.')
            sage: s.to_latex()
            '\\begin{array}{|*{2}{*{2}{r}|}}\\hline\n &4& & \\\\\n3&2& & \\\\\\hline\n & &1&4\\\\\n & &3& \\\\\\hline\n\\end{array}'
        """
        n = self.n
        nsquare = n*n
        array = []
        array.append('\\begin{array}{|*{%s}{*{%s}{r}|}}\\hline\n' % (n, n))
        gen = (x for x in self.puzzle)
        for row in range(nsquare):
            for col in range(nsquare):
                entry = next(gen)
                array.append((str(entry) if entry else ' '))
                array.append(('' if col == nsquare - 1 else '&'))
            array.append(('\\\\\n' if (row+1) % n else '\\\\\\hline\n'))
        array.append('\\end{array}')
        return ''.join(array)

    def solve(self, algorithm='dlx'):
        r"""
        Return a generator object for the solutions of a Sudoku puzzle.

        INPUT:

        - algorithm -- default = ``'dlx'``, specify choice of solution algorithm. The
          two possible algorithms are ``'dlx'`` and ``'backtrack'``.

        OUTPUT:

        A generator that provides all solutions, as objects of
        the :class:`~sage.games.sudoku.Sudoku` class.

        Calling ``next()`` on the returned generator just once will find
        a solution, presuming it exists, otherwise it will return a
        ``StopIteration`` exception.  The generator may be used for
        iteration or wrapping the generator with ``list()`` will return
        all of the solutions as a list.  Solutions are returned as
        new objects of the :class:`~sage.games.sudoku.Sudoku` class,
        so may be printed or converted using other methods in this class.

        Generally, the DLX algorithm is very fast and very consistent.
        The backtrack algorithm is very variable in its performance,
        on some occasions markedly faster than DLX but usually slower
        by a similar factor, with the potential to be orders of magnitude
        slower.  See the docstrings for the
        :meth:`~sage.games.sudoku.Sudoku.dlx` and
        :meth:`~sage.games.sudoku.Sudoku.backtrack_all`
        methods for further discussions and examples of performance.
        Note that the backtrack algorithm is limited to puzzles of
        size `16\times 16` or smaller.

        EXAMPLES:

        This puzzle has 5 solutions, but the first one returned by each algorithm are identical. ::

            sage: h = Sudoku('8..6..9.5.............2.31...7318.6.24.....73...........279.1..5...8..36..3......')
            sage: h
            +-----+-----+-----+
            |8    |6    |9   5|
            |     |     |     |
            |     |  2  |3 1  |
            +-----+-----+-----+
            |    7|3 1 8|  6  |
            |2 4  |     |  7 3|
            |     |     |     |
            +-----+-----+-----+
            |    2|7 9  |1    |
            |5    |  8  |  3 6|
            |    3|     |     |
            +-----+-----+-----+
            sage: next(h.solve(algorithm='backtrack'))
            +-----+-----+-----+
            |8 1 4|6 3 7|9 2 5|
            |3 2 5|1 4 9|6 8 7|
            |7 9 6|8 2 5|3 1 4|
            +-----+-----+-----+
            |9 5 7|3 1 8|4 6 2|
            |2 4 1|9 5 6|8 7 3|
            |6 3 8|2 7 4|5 9 1|
            +-----+-----+-----+
            |4 6 2|7 9 3|1 5 8|
            |5 7 9|4 8 1|2 3 6|
            |1 8 3|5 6 2|7 4 9|
            +-----+-----+-----+
            sage: next(h.solve(algorithm='dlx'))
            +-----+-----+-----+
            |8 1 4|6 3 7|9 2 5|
            |3 2 5|1 4 9|6 8 7|
            |7 9 6|8 2 5|3 1 4|
            +-----+-----+-----+
            |9 5 7|3 1 8|4 6 2|
            |2 4 1|9 5 6|8 7 3|
            |6 3 8|2 7 4|5 9 1|
            +-----+-----+-----+
            |4 6 2|7 9 3|1 5 8|
            |5 7 9|4 8 1|2 3 6|
            |1 8 3|5 6 2|7 4 9|
            +-----+-----+-----+

        Gordon Royle maintains a list of 48072 Sudoku puzzles that each has
        a unique solution and exactly 17 "hints" (initially filled boxes).
        At this writing (May 2009) there is no known 16-hint puzzle with
        exactly one solution.  [sudoku:royle]_  This puzzle is number 3000
        in his database.  We solve it twice. ::

            sage: b = Sudoku('8..6..9.5.............2.31...7318.6.24.....73...........279.1..5...8..36..3......')
            sage: next(b.solve(algorithm='dlx')) == next(b.solve(algorithm='backtrack'))
            True


        These are the first 10 puzzles in a list of "Top 95" puzzles,
        [sudoku:top95]_ which we use to show that the two available algorithms obtain
        the same solution for each. ::

            sage: top =['4.....8.5.3..........7......2.....6.....8.4......1.......6.3.7.5..2.....1.4......',
            ....:       '52...6.........7.13...........4..8..6......5...........418.........3..2...87.....',
            ....:       '6.....8.3.4.7.................5.4.7.3..2.....1.6.......2.....5.....8.6......1....',
            ....:       '48.3............71.2.......7.5....6....2..8.............1.76...3.....4......5....',
            ....:       '....14....3....2...7..........9...3.6.1.............8.2.....1.4....5.6.....7.8...',
            ....:       '......52..8.4......3...9...5.1...6..2..7........3.....6...1..........7.4.......3.',
            ....:       '6.2.5.........3.4..........43...8....1....2........7..5..27...........81...6.....',
            ....:       '.524.........7.1..............8.2...3.....6...9.5.....1.6.3...........897........',
            ....:       '6.2.5.........4.3..........43...8....1....2........7..5..27...........81...6.....',
            ....:       '.923.........8.1...........1.7.4...........658.........6.5.2...4.....7.....9.....']
            sage: p = [Sudoku(top[i]) for i in range(10)]
            sage: verify = [next(p[i].solve(algorithm='dlx')) == next(p[i].solve(algorithm='backtrack')) for i in range(10)]
            sage: verify == [True]*10
            True

        TESTS:

        A `25\times 25` puzzle that the backtrack algorithm is not equipped to handle.  Since ``solve`` returns a generator this test will not go boom until we ask for a solution with ``next``. ::

            sage: too_big = Sudoku([0]*625)
            sage: next(too_big.solve(algorithm='backtrack'))
            Traceback (most recent call last):
            ...
            ValueError: The Sudoku backtrack algorithm is limited to puzzles of size 16 or smaller.

        An attempt to use a non-existent algorithm. ::

            sage: next(Sudoku([0]).solve(algorithm='bogus'))
            Traceback (most recent call last):
            ...
            NotImplementedError: bogus is not an algorithm for Sudoku puzzles
        """
        if algorithm == 'backtrack':
            if self.n > 4:
                raise ValueError('The Sudoku backtrack algorithm is limited to puzzles of size 16 or smaller.')
            else:
                gen = self.backtrack()
        elif algorithm == 'dlx':
            gen = self.dlx()
        else:
            raise NotImplementedError('%s is not an algorithm for Sudoku puzzles' % algorithm)
        for soln in gen:
            yield Sudoku(soln, verify_input = 'False')

    def backtrack(self):
        r"""
        Return a generator which iterates through all solutions of a Sudoku puzzle.

        This function is intended to be called from the
        :func:`~sage.games.sudoku.Sudoku.solve` method
        when the ``algorithm='backtrack'`` option is specified.
        However it may be called directly as a method of an
        instance of a Sudoku puzzle.

        At this point, this method calls
        :func:`~sage.games.sudoku_backtrack.backtrack_all` which
        constructs *all* of the solutions as a list.  Then the
        present method just returns the items of the list one at
        a time.  Once Cython supports closures and a yield statement
        is supported, then the contents of ``backtrack_all()``
        may be subsumed into this method and the
        :mod:`sage.games.sudoku_backtrack` module can be removed.

        This routine can have wildly variable performance, with a
        factor of 4000 observed between the fastest and slowest
        `9\times 9` examples tested. Examples designed to perform
        poorly for naive backtracking, will do poorly
        (such as ``d`` below).  However, examples meant to be
        difficult for humans often do very well, with a factor
        of 5 improvement over the `DLX` algorithm.

        Without dynamically allocating arrays in the Cython version,
        we have limited this function to `16\times 16` puzzles.
        Algorithmic details are in the
        :mod:`sage.games.sudoku_backtrack` module.

        EXAMPLES:

        This example was reported to be very difficult for human solvers.
        This algorithm works very fast on it, at about half the time
        of the DLX solver. [sudoku:escargot]_ ::

            sage: g = Sudoku('1....7.9..3..2...8..96..5....53..9...1..8...26....4...3......1..4......7..7...3..')
            sage: print(g)
            +-----+-----+-----+
            |1    |    7|  9  |
            |  3  |  2  |    8|
            |    9|6    |5    |
            +-----+-----+-----+
            |    5|3    |9    |
            |  1  |  8  |    2|
            |6    |    4|     |
            +-----+-----+-----+
            |3    |     |  1  |
            |  4  |     |    7|
            |    7|     |3    |
            +-----+-----+-----+
            sage: print(next(g.solve(algorithm='backtrack')))
            +-----+-----+-----+
            |1 6 2|8 5 7|4 9 3|
            |5 3 4|1 2 9|6 7 8|
            |7 8 9|6 4 3|5 2 1|
            +-----+-----+-----+
            |4 7 5|3 1 2|9 8 6|
            |9 1 3|5 8 6|7 4 2|
            |6 2 8|7 9 4|1 3 5|
            +-----+-----+-----+
            |3 5 6|4 7 8|2 1 9|
            |2 4 1|9 3 5|8 6 7|
            |8 9 7|2 6 1|3 5 4|
            +-----+-----+-----+

        This example has no entries in the top row and a half,
        and the top row of the solution is ``987654321`` and
        therefore a backtracking approach is slow, taking about
        750 times as long as the DLX solver. [sudoku:wikipedia]_ ::

            sage: c = Sudoku('..............3.85..1.2.......5.7.....4...1...9.......5......73..2.1........4...9')
            sage: print(c)
            +-----+-----+-----+
            |     |     |     |
            |     |    3|  8 5|
            |    1|  2  |     |
            +-----+-----+-----+
            |     |5   7|     |
            |    4|     |1    |
            |  9  |     |     |
            +-----+-----+-----+
            |5    |     |  7 3|
            |    2|  1  |     |
            |     |  4  |    9|
            +-----+-----+-----+
            sage: print(next(c.solve(algorithm='backtrack')))
            +-----+-----+-----+
            |9 8 7|6 5 4|3 2 1|
            |2 4 6|1 7 3|9 8 5|
            |3 5 1|9 2 8|7 4 6|
            +-----+-----+-----+
            |1 2 8|5 3 7|6 9 4|
            |6 3 4|8 9 2|1 5 7|
            |7 9 5|4 6 1|8 3 2|
            +-----+-----+-----+
            |5 1 9|2 8 6|4 7 3|
            |4 7 2|3 1 9|5 6 8|
            |8 6 3|7 4 5|2 1 9|
            +-----+-----+-----+
        """
        from .sudoku_backtrack import backtrack_all
        solutions = backtrack_all(self.n, self.puzzle)
        for soln in solutions:
            yield soln

    def dlx(self, count_only=False):
        r"""
        Return a generator that iterates through all solutions of a Sudoku puzzle.

        INPUT:

        - count_only -- boolean, default = False.
          If set to ``True`` the generator returned as output will
          simply generate ``None`` for each solution, so the
          calling routine can count these.

        OUTPUT:

        A generator that iterates over all the solutions.

        This function is intended to be called from the
        :func:`~sage.games.sudoku.Sudoku.solve` method
        with the ``algorithm='dlx'`` option. However it
        may be called directly as a method of an instance
        of a Sudoku puzzle if speed is important and you
        do not need automatic conversions on the output
        (or even just want to count solutions without looking
        at them).  In this case, inputting a puzzle as a list,
        with ``verify_input=False`` is the fastest way to
        create a puzzle.

        Or if only one solution is needed it can be obtained with
        one call to ``next()``, while the existence of a solution
        can be tested by catching the ``StopIteration`` exception
        with a ``try``.  Calling this particular method returns
        solutions as lists, in row-major order. It is up to you
        to work with this list for your own purposes. If you want
        fancier formatting tools, use the
        :func:`~sage.games.sudoku.Sudoku.solve` method, which
        returns a generator that creates
        :class:`sage.games.sudoku.Sudoku` objects.

        EXAMPLES:

        A `9\times 9` known to have one solution.  We get the one
        solution and then check to see if there are more or not. ::

            sage: e = Sudoku('4.....8.5.3..........7......2.....6.....8.4......1.......6.3.7.5..2.....1.4......')
            sage: print(next(e.dlx()))
            [4, 1, 7, 3, 6, 9, 8, 2, 5, 6, 3, 2, 1, 5, 8, 9, 4, 7, 9, 5, 8, 7, 2, 4, 3, 1, 6, 8, 2, 5, 4, 3, 7, 1, 6, 9, 7, 9, 1, 5, 8, 6, 4, 3, 2, 3, 4, 6, 9, 1, 2, 7, 5, 8, 2, 8, 9, 6, 4, 3, 5, 7, 1, 5, 7, 3, 2, 9, 1, 6, 8, 4, 1, 6, 4, 8, 7, 5, 2, 9, 3]
            sage: len(list(e.dlx()))
            1

        A `9\times 9` puzzle with multiple solutions.
        Once with actual solutions, once just to count. ::

            sage: h = Sudoku('8..6..9.5.............2.31...7318.6.24.....73...........279.1..5...8..36..3......')
            sage: len(list(h.dlx()))
            5
            sage: len(list(h.dlx(count_only=True)))
            5

        A larger puzzle, with multiple solutions, but we just get one. ::

            sage: j = Sudoku('....a..69.3....1d.2...8....e.4....b....5..c.......7.......g...f....1.e..2.b.8..3.......4.d.....6.........f..7.g..9.a..c...5.....8..f.....1..e.79.c....b.....2...6.....g.7......84....3.d..a.5....5...7..e...ca.....3.1.......b......f....4...d..e..g.92.6..8....')
            sage: print(next(j.dlx()))
            [5, 15, 16, 14, 10, 13, 7, 6, 9, 2, 3, 4, 11, 8, 12, 1, 13, 3, 2, 12, 11, 16, 8, 15, 1, 6, 7, 14, 10, 4, 9, 5, 1, 10, 11, 6, 9, 4, 3, 5, 15, 8, 12, 13, 16, 7, 14, 2, 9, 8, 7, 4, 12, 2, 1, 14, 10, 5, 16, 11, 6, 3, 15, 13, 12, 16, 4, 1, 13, 14, 9, 10, 2, 7, 11, 6, 8, 15, 5, 3, 3, 14, 5, 7, 16, 11, 15, 4, 12, 13, 8, 9, 1, 2, 10, 6, 2, 6, 13, 11, 1, 8, 5, 3, 4, 15, 14, 10, 7, 9, 16, 12, 15, 9, 8, 10, 2, 6, 12, 7, 3, 16, 5, 1, 4, 14, 13, 11, 8, 11, 3, 15, 5, 10, 4, 2, 13, 1, 6, 12, 14, 16, 7, 9, 16, 12, 14, 13, 7, 15, 11, 1, 8, 9, 4, 5, 2, 6, 3, 10, 6, 2, 10, 5, 14, 12, 16, 9, 7, 11, 15, 3, 13, 1, 4, 8, 4, 7, 1, 9, 8, 3, 6, 13, 16, 14, 10, 2, 5, 12, 11, 15, 11, 5, 9, 8, 6, 7, 13, 16, 14, 3, 1, 15, 12, 10, 2, 4, 7, 13, 15, 3, 4, 1, 10, 8, 5, 12, 2, 16, 9, 11, 6, 14, 10, 1, 6, 2, 15, 5, 14, 12, 11, 4, 9, 7, 3, 13, 8, 16, 14, 4, 12, 16, 3, 9, 2, 11, 6, 10, 13, 8, 15, 5, 1, 7]

        The puzzle ``h`` from above, but purposely made unsolvable with addition in second entry. ::

            sage: hbad = Sudoku('82.6..9.5.............2.31...7318.6.24.....73...........279.1..5...8..36..3......')
            sage: len(list(hbad.dlx()))
            0
            sage: next(hbad.dlx())
            Traceback (most recent call last):
            ...
            StopIteration

        A stupidly small puzzle to test the lower limits of arbitrary sized input. ::

            sage: s = Sudoku('.')
            sage: print(next(s.solve(algorithm='dlx')))
            +-+
            |1|
            +-+

        ALGORITHM:

        The ``DLXCPP`` solver finds solutions to the exact-cover
        problem with a "Dancing Links" backtracking algorithm.
        Given a `0-1` matrix, the solver finds a subset of the
        rows that sums to the all `1`'s vector.  The columns
        correspond to conditions, or constraints, that must be
        met by a solution, while the rows correspond to some
        collection of choices, or decisions.  A `1` in a row
        and column indicates that the choice corresponding to
        the row meets the condition corresponding to the column.

        So here, we code the notion of a Sudoku puzzle, and the
        hints already present, into such a `0-1` matrix.  Then the
        :class:`sage.combinat.matrices.dlxcpp.DLXCPP` solver makes
        the choices for the blank entries.
        """
        from sage.combinat.matrices.dlxcpp import DLXCPP

        n = self.n
        nsquare = n*n
        nfour = nsquare*nsquare

        # Boxes of the grid are numbered in row-major order
        # ``rcbox`` simply maps a row-column index pair to the box number it lives in
        rcbox   = [ [i//n + n*(j//n) for i in range(nsquare)] for j in range(nsquare)]

        # Every entry in a Sudoku puzzle satisfies four constraints
        # Every location has a single entry, and each row, column and box has each symbol once
        # These arrays can be thought of as assigning ID numbers to these constraints,
        # and correspond to column numbers of the `0-1` matrix describing the exact cover
        rows    = [[i+j for i in range(nsquare)] for j in range(0, nfour, nsquare)]
        cols    = [[i+j for i in range(nsquare)] for j in range(nfour, 2*nfour, nsquare)]
        boxes   = [[i+j for i in range(nsquare)] for j in range(2*nfour, 3*nfour, nsquare)]
        rowcol  = [[i+j for i in range(nsquare)] for j in range(3*nfour, 4*nfour, nsquare)]

        def make_row(row, col, entry):
            r"""
            Constructs a row of the `0-1` matrix describing
            the exact cover constraints for a Sudoku puzzle.

            If a (zero-based) ``entry`` is placed in location
            ``(row, col)`` of a Sudoku puzzle, then exactly four
            constraints are satisfied: the location has this entry,
            and the entry occurs in the row, column and box.
            This method looks up the constraint IDs for each of
            these four constraints, and returns a list of these four IDs.

            TESTS::

                sage: h = Sudoku('8..6..9.5.............2.31...7318.6.24.....73...........279.1..5...8..36..3......')
                sage: len(list(h.solve(algorithm='dlx')))  # indirect doctest
                5
            """
            box = rcbox[row][col]
            return [rows[row][entry], cols[col][entry], boxes[box][entry], rowcol[row][col]]

        # Construct the sparse `0-1` matrix for the exact cover formulation as the ``ones`` array
        # ``rowinfo`` remembers the location and entry that led to the row being added to the matrix
        rowinfo = []
        ones = []
        gen = (entry for entry in self.puzzle)
        for row in range(nsquare):
            for col in range(nsquare):
                puzz = next(gen)
                # All (zero-based) entries are possible, or only one is possible
                entries = ([puzz-1] if puzz else range(nsquare))
                for entry in entries:
                    ones.append(make_row(row, col, entry))
                    rowinfo.append((row, col, entry))

        # ``DLXCPP`` will solve the exact cover problem for the ``ones`` matrix
        # ``cover`` will contain a subset of the row indices so that the sum of the rows is the all `1`'s vector
        # These rows will represent the original hints, plus a single entry in every other location,
        # consistent with the requirements imposed on a solution to a Sudoku puzzle
        for cover in DLXCPP(ones):
            if not(count_only):
                solution = [0]*nfour
                for r in cover:
                    row, col, entry = rowinfo[r]
                    solution[row*nsquare+col] = entry+1
                yield solution
            else:
                yield None
