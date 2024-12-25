r"""
Combinatorial diagrams

A combinatorial diagram is a collection of cells `(i,j)` indexed by pairs of
natural numbers.

For arbitrary diagrams, see :class:`Diagram`. There are also two other specific
types of diagrams implemented here. They are northwest diagrams
(:class:`NorthwestDiagram`) and Rothe diagrams (:func:`RotheDiagram`, a special
kind of northwest diagram).

AUTHORS:

- Trevor K. Karn (2022-08-01): initial version
"""

# ****************************************************************************
#       Copyright (C) 2022 Trevor K. Karn <karnx018 (at) umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from itertools import product

from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.composition import Composition
from sage.combinat.partition import Partition
from sage.combinat.permutation import Permutations
from sage.combinat.skew_partition import SkewPartition
from sage.combinat.skew_tableau import SkewTableaux
from sage.combinat.tableau import Tableau
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_import import lazy_import
from sage.structure.element import Matrix
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

lazy_import('sage.combinat.tiling', 'Polyomino')


class Diagram(ClonableArray, metaclass=InheritComparisonClasscallMetaclass):
    r"""
    Combinatorial diagrams with positions indexed by rows and columns.

    The positions are indexed by rows and columns as in a matrix. For example,
    a Ferrers diagram is a diagram obtained from a partition
    `\lambda = (\lambda_0, \lambda_1, \ldots, \lambda_{\ell})`, where the
    cells are in rows `i` for `0 \leq i \leq \ell` and the cells in row `i`
    consist of `(i,j)` for `0 \leq j < \lambda_i`. In English notation, the
    indices are read from top left to bottom right as in a matrix.

    Indexing conventions are the same as
    :class:`~sage.combinat.partition.Partition`. Printing the diagram of a
    partition, however, will always be in English notation.

    EXAMPLES:

    To create an arbitrary diagram, pass a list of all cells::

        sage: from sage.combinat.diagram import Diagram
        sage: cells = [(0,0), (0,1), (1,0), (1,1), (4,4), (4,5), (4,6), (5,4), (7, 6)]
        sage: D = Diagram(cells); D
        [(0, 0), (0, 1), (1, 0), (1, 1), (4, 4), (4, 5), (4, 6), (5, 4), (7, 6)]

    We can visualize the diagram by printing ``O``'s and ``.``'s. ``O``'s are
    present in the cells which are present in the diagram and a ``.`` represents
    the absence of a cell in the diagram::

        sage: D.pp()
        O O . . . . .
        O O . . . . .
        . . . . . . .
        . . . . . . .
        . . . . O O O
        . . . . O . .
        . . . . . . .
        . . . . . . O

    We can also check if certain cells are contained in a given diagram::

        sage: (1, 0) in D
        True
        sage: (2, 2) in D
        False

    If you know that there are entire empty rows or columns at the end of the
    diagram, you can manually pass them with keyword arguments ``n_rows=`` or
    ``n_cols=``::

        sage: Diagram([(0,0), (0,3), (2,2), (2,4)]).pp()
        O . . O .
        . . . . .
        . . O . O
        sage: Diagram([(0,0), (0,3), (2,2), (2,4)], n_rows=6, n_cols=6).pp()
        O . . O . .
        . . . . . .
        . . O . O .
        . . . . . .
        . . . . . .
        . . . . . .
    """
    @staticmethod
    def __classcall_private__(self, cells, n_rows=None, n_cols=None, check=True):
        r"""
        Normalize the input so that it lives in the correct parent.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,0), (0,3), (2,2), (2,4)])
            sage: D.parent()
            Combinatorial diagrams
        """
        return Diagrams()(cells, n_rows, n_cols, check)

    def __init__(self, parent, cells, n_rows=None, n_cols=None, check=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D1 = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: D1.cells()
            [(0, 2), (0, 3), (1, 1), (3, 2)]
            sage: D1.nrows()
            4
            sage: D1.ncols()
            4
            sage: TestSuite(D1).run()

        We can specify the number of rows and columns explicitly,
        in case they are supposed to be empty::

            sage: D2 = Diagram([(0,2),(0,3),(1,1),(3,2)], n_cols=5)
            sage: D2.cells()
            [(0, 2), (0, 3), (1, 1), (3, 2)]
            sage: D2.ncols()
            5
            sage: D2.pp()
            . . O O .
            . O . . .
            . . . . .
            . . O . .
            sage: TestSuite(D2).run()
        """
        self._cells = frozenset(cells)

        if self._cells:
            # minimum possible number of rows/cols
            N_rows = max(c[0] for c in self._cells)
            N_cols = max(c[1] for c in self._cells)
        else:  # if there are no cells
            N_rows = -1
            N_cols = -1

        if n_rows is not None:
            if n_rows <= N_rows:
                raise ValueError('n_rows is too small')
            self._n_rows = n_rows
        else:
            self._n_rows = N_rows + 1
        if n_cols is not None:
            if n_cols <= N_cols:
                raise ValueError('n_cols is too small')
            self._n_cols = n_cols
        else:
            self._n_cols = N_cols + 1

        self._n_nonempty_rows = len(set(i for i, j in self._cells))
        self._n_nonempty_cols = len(set(j for i, j in self._cells))

        ClonableArray.__init__(self, parent, sorted(cells), check)

    def pp(self):
        r"""
        Return a visualization of the diagram.

        Cells which are present in the
        diagram are filled with a ``O``. Cells which are not present in the
        diagram are filled with a ``.``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: Diagram([(0,0), (0,3), (2,2), (2,4)]).pp()
            O . . O .
            . . . . .
            . . O . O
            sage: Diagram([(0,0), (0,3), (2,2), (2,4)], n_rows=6, n_cols=6).pp()
            O . . O . .
            . . . . . .
            . . O . O .
            . . . . . .
            . . . . . .
            . . . . . .
            sage: Diagram([]).pp()
            -
        """
        if self._n_rows == 0 or self._n_cols == 0:
            print('-')
            return
        print("\n".join(self._pretty_print()))

    def _ascii_art_(self):
        r"""
        Return a visualization of the diagram.

        Cells which are present in the
        diagram are filled with a ``O``. Cells which are not present in the
        diagram are filled with a ``.``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: ascii_art(Diagram([(0,0), (0,3), (2,2), (2,4)]))
            O . . O .
            . . . . .
            . . O . O
            sage: ascii_art(Diagram([(0,0), (0,3), (2,2), (2,4)], n_rows=6, n_cols=6))
            O . . O . .
            . . . . . .
            . . O . O .
            . . . . . .
            . . . . . .
            . . . . . .
            sage: ascii_art(Diagram([]))
            -
        """
        from sage.typeset.ascii_art import ascii_art
        if self._n_rows == 0 or self._n_cols == 0:
            return ascii_art("-")
        return ascii_art("\n".join(self._pretty_print()))

    def _unicode_art_(self):
        r"""
        Return a unicode visualization of the diagram.

        Cells which are present in the
        diagram are filled with a crossed box. Cells which are not present in the
        diagram are filled with an empty box.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: unicode_art(Diagram([(0,0), (0,3), (2,2), (2,4)]))
            ┌─┬─┬─┬─┬─┐
            │X│ │ │X│ │
            ├─┼─┼─┼─┼─┤
            │ │ │ │ │ │
            ├─┼─┼─┼─┼─┤
            │ │ │X│ │X│
            └─┴─┴─┴─┴─┘
            sage: unicode_art(Diagram([(0,0), (0,3), (2,2), (2,4)], n_rows=6, n_cols=6))
            ┌─┬─┬─┬─┬─┬─┐
            │X│ │ │X│ │ │
            ├─┼─┼─┼─┼─┼─┤
            │ │ │ │ │ │ │
            ├─┼─┼─┼─┼─┼─┤
            │ │ │X│ │X│ │
            ├─┼─┼─┼─┼─┼─┤
            │ │ │ │ │ │ │
            ├─┼─┼─┼─┼─┼─┤
            │ │ │ │ │ │ │
            ├─┼─┼─┼─┼─┼─┤
            │ │ │ │ │ │ │
            └─┴─┴─┴─┴─┴─┘
            sage: unicode_art(Diagram([]))
            ∅
        """
        from sage.typeset.unicode_art import unicode_art
        if self._n_rows == 0 or self._n_cols == 0:
            return unicode_art("∅")

        ndivs = self._n_cols - 1
        cell = "│X"
        empty = "│ "
        it = self._pretty_print(cell, empty)
        ret = "┌─" + "┬─"*ndivs + "┐"
        ret += "\n" + next(it) + "│"
        for row in it:
            ret += "\n├─" + "┼─"*ndivs + "┤"
            ret += "\n" + row + "│"
        ret += "\n└─" + "┴─"*ndivs + "┘"
        return unicode_art(ret)

    def _pretty_print(self, cell='O ', empty='. '):
        r"""
        Return a visualization of the diagram.

        Cells which are present in the
        diagram are filled with ``cell``. Cells which are not present in the
        diagram are filled with ``empty``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: "\n".join(Diagram([(0,0), (0,3), (2,2), (2,4)])._pretty_print('x ','. '))
            'x . . x . \n. . . . . \n. . x . x '
            sage: "\n".join(Diagram([(0,0), (0,3), (2,2), (2,4)], n_rows=6, n_cols=6)._pretty_print('x ','. '))
            'x . . x . . \n. . . . . . \n. . x . x . \n. . . . . . \n. . . . . . \n. . . . . . '
        """
        for i in range(self._n_rows):
            output_str = ''
            for j in range(self._n_cols):
                if (i, j) in self:
                    output_str += cell
                else:
                    output_str += empty
            yield output_str

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: latex(Diagram([]))
            {\emptyset}
            sage: latex(Diagram([(0,0), (0,3), (2,2), (2,4)]))
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{5}{p{0.6ex}}}\cline{1-1}\cline{4-4}
            \lr{\phantom{x}}&&&\lr{\phantom{x}}&\\\cline{1-1}\cline{4-4}
            &&&&\\\cline{3-3}\cline{5-5}
            &&\lr{\phantom{x}}&&\lr{\phantom{x}}\\\cline{3-3}\cline{5-5}
            \end{array}$}
            }
        """
        if self._n_rows == 0 or self._n_cols == 0:
            return "{\\emptyset}"

        lr = r'\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}'

        array = [[("\\phantom{x}" if (i, j) in self else None)
                  for j in range(self._n_cols)]
                 for i in range(self._n_rows)]

        def end_line(r):
            # give the line ending to row ``r``
            if r == 0:
                return "".join(r'\cline{%s-%s}' % (i+1, i+1)
                               for i, j in enumerate(array[0]) if j is not None)
            elif r == len(array):
                return r"\\" + "".join(r'\cline{%s-%s}' % (i+1, i+1)
                                       for i, j in enumerate(array[r-1]) if j is not None)
            else:
                out = r"\\" + "".join(r'\cline{%s-%s}' % (i+1, i+1)
                                      for i, j in enumerate(array[r-1]) if j is not None)
                out += "".join(r'\cline{%s-%s}' % (i+1, i+1)
                               for i, j in enumerate(array[r]) if j is not None)
                return out

        tex = r'\raisebox{-.6ex}{$\begin{array}[b]{*{%s}{p{0.6ex}}}' % (max(map(len, array)))
        tex += end_line(0)+'\n'
        for r in range(len(array)):
            tex += '&'.join('' if c is None else r'\lr{%s}' % (c,) for c in array[r])
            tex += end_line(r+1)+'\n'
        return '{%s\n%s\n}' % (lr, tex+r'\end{array}$}')

    def number_of_rows(self):
        r"""
        Return the total number of rows of ``self``.

        EXAMPLES:

        The following example has three rows which are filled, but they
        are contained in rows 0 to 3 (for a total of four)::

            sage: from sage.combinat.diagram import Diagram
            sage: D1 = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: D1.number_of_rows()
            4
            sage: D1.nrows()
            4

        The total number of rows includes including those which are empty.
        We can also include empty rows at the end::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,2),(0,3),(1,1),(3,2)], n_rows=6)
            sage: D.number_of_rows()
            6
            sage: D.pp()
            . . O O
            . O . .
            . . . .
            . . O .
            . . . .
            . . . .
        """
        return self._n_rows

    nrows = number_of_rows

    def number_of_cols(self):
        r"""
        Return the total number of rows of ``self``.

        EXAMPLES:

        The following example has three columns which are filled, but they
        are contained in rows 0 to 3 (for a total of four)::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: D.number_of_cols()
            4
            sage: D.ncols()
            4

        We can also include empty columns at the end::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,2),(0,3),(1,1),(3,2)], n_cols=6)
            sage: D.number_of_cols()
            6
            sage: D.pp()
            . . O O . .
            . O . . . .
            . . . . . .
            . . O . . .
        """
        return self._n_cols

    ncols = number_of_cols

    def cells(self):
        r"""
        Return a ``list`` of the cells contained in the diagram ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D1 = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: D1.cells()
            [(0, 2), (0, 3), (1, 1), (3, 2)]
        """
        return sorted(self._cells)

    def number_of_cells(self):
        r"""
        Return the total number of cells contained in the diagram ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D1 = Diagram([(0,2),(0,3),(1,1),(3,2)])
            sage: D1.number_of_cells()
            4
            sage: D1.n_cells()
            4
        """
        return len(self._cells)

    n_cells = number_of_cells

    size = number_of_cells

    def check(self):
        r"""
        Check that this is a valid diagram.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,0), (0,3), (2,2), (2,4)])
            sage: D.check()

        In the next two examples, a bad diagram is passed.
        The first example fails because one cell is indexed by negative
        integers::

            sage: D = Diagram([(0,0), (0,-3), (2,2), (2,4)])
            Traceback (most recent call last):
            ...
            ValueError: diagrams must be indexed by nonnegative integers

        The next example fails because one cell is indexed by rational
        numbers::

            sage: D = Diagram([(0,0), (0,3), (2/3,2), (2,4)])
            Traceback (most recent call last):
            ...
            ValueError: diagrams must be indexed by nonnegative integers
        """
        from sage.sets.non_negative_integers import NonNegativeIntegers
        NN = NonNegativeIntegers()
        if not all(i in NN for c in self._cells for i in c):
            raise ValueError("diagrams must be indexed by nonnegative integers")

    def specht_module(self, base_ring=None):
        r"""
        Return the Specht module corresponding to ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,0), (1,1), (2,2), (2,3)])
            sage: SM = D.specht_module(QQ)                                              # needs sage.modules
            sage: s = SymmetricFunctions(QQ).s()                                        # needs sage.modules
            sage: s(SM.frobenius_image())                                               # needs sage.modules
            s[2, 1, 1] + s[2, 2] + 2*s[3, 1] + s[4]
        """
        from sage.combinat.specht_module import SpechtModule
        from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ
        R = SymmetricGroupAlgebra(base_ring, len(self))
        return SpechtModule(R, self)

    def specht_module_dimension(self, base_ring=None):
        r"""
        Return the dimension of the Specht module corresponding to ``self``.

        INPUT:

        - ``base_ring`` -- (default: `\QQ`) the base ring

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,0), (1,1), (2,2), (2,3)])
            sage: D.specht_module_dimension()                                           # needs sage.modules
            12
            sage: D.specht_module(QQ).dimension()                                       # needs sage.modules
            12
        """
        from sage.combinat.specht_module import specht_module_rank
        return specht_module_rank(self, base_ring)


class Diagrams(UniqueRepresentation, Parent):
    r"""
    The class of combinatorial diagrams.

    A *combinatorial diagram* is a set of cells indexed by pairs of natural
    numbers. Calling an instance of :class:`Diagrams` is one way to construct
    diagrams.

    EXAMPLES::

        sage: from sage.combinat.diagram import Diagrams
        sage: Dgms = Diagrams()
        sage: D = Dgms([(0,0), (0,3), (2,2), (2,4)])
        sage: D.parent()
        Combinatorial diagrams
    """

    def __init__(self, category=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagrams
            sage: Dgms = Diagrams(); Dgms
            Combinatorial diagrams

        TESTS::

            sage: TestSuite(Dgms).run()
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets().or_subcategory(category))

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagrams
            sage: I = iter(Diagrams())
            sage: for i in range(10):
            ....:   print(next(I))
            []
            [(0, 0)]
            [(1, 0)]
            [(0, 0), (1, 0)]
            [(0, 1)]
            [(0, 0), (0, 1)]
            [(0, 1), (1, 0)]
            [(0, 0), (0, 1), (1, 0)]
            [(2, 0)]
            [(0, 0), (2, 0)]
            sage: next(I).parent()
            Combinatorial diagrams

            sage: from sage.combinat.diagram import NorthwestDiagrams
            sage: I = iter(NorthwestDiagrams())
            sage: for i in range(20):
            ....:   print(next(I))
            []
            [(0, 0)]
            [(1, 0)]
            [(0, 0), (1, 0)]
            [(0, 1)]
            [(0, 0), (0, 1)]
            [(0, 0), (0, 1), (1, 0)]
            [(2, 0)]
            [(0, 0), (2, 0)]
            [(1, 0), (2, 0)]
            [(0, 0), (1, 0), (2, 0)]
            [(0, 0), (0, 1), (2, 0)]
            [(0, 0), (0, 1), (1, 0), (2, 0)]
            [(1, 1)]
            [(0, 0), (1, 1)]
            [(1, 0), (1, 1)]
            [(0, 0), (1, 0), (1, 1)]
            [(0, 1), (1, 1)]
            [(0, 0), (0, 1), (1, 1)]
            [(0, 0), (0, 1), (1, 0), (1, 1)]
        """
        from sage.sets.non_negative_integers import NonNegativeIntegers
        from sage.categories.cartesian_product import cartesian_product
        from sage.combinat.subset import subsets
        # the product of positive integers automatically implements an
        # an enumeration which allows us to get out of the first column
        N = NonNegativeIntegers()
        NxN = cartesian_product([N, N])
        X = subsets(NxN)
        while True:
            cells = next(X)
            try:
                yield self.element_class(self, tuple((i, j) for i, j in cells))
            except ValueError:
                # if cells causes the .check method of a
                # subclass to fail, just go to the next one
                pass

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagrams
            sage: Dgms = Diagrams(); Dgms
            Combinatorial diagrams
        """
        return 'Combinatorial diagrams'

    def _element_constructor_(self, cells, n_rows=None, n_cols=None, check=True):
        r"""
        Cosntruct an element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagrams
            sage: Dgms = Diagrams()
            sage: Dgms([(0,1),(2,2)]).pp()
            . O .
            . . .
            . . O


            sage: from sage.combinat.tiling import Polyomino                            # needs sage.modules
            sage: p = Polyomino([(0,0),(1,0),(1,1),(1,2)])                              # needs sage.modules
            sage: Dgms(p).pp()                                                          # needs sage.modules
            O . .
            O O O

            sage: from sage.combinat.composition import Composition
            sage: a = Composition([4,2,0,2,4])
            sage: Dgms(a).pp()
            O O O O
            O O . .
            . . . .
            O O . .
            O O O O

            sage: M = Matrix([[1,1,1,1],[1,1,0,0],[0,0,0,0],[1,1,0,0],[1,1,1,1]])       # needs sage.modules
            sage: Dgms(M).pp()                                                          # needs sage.modules
            O O O O
            O O . .
            . . . .
            O O . .
            O O O O

        TESTS::

            sage: TestSuite(Dgms).run()
        """
        if isinstance(cells, Polyomino):
            return self.from_polyomino(cells)
        if isinstance(cells, Composition):
            return self.from_composition(cells)
        if isinstance(cells, Matrix):
            return self.from_zero_one_matrix(cells)

        return self.element_class(self, cells, n_rows, n_cols, check)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import Diagrams
            sage: Dgms = Diagrams()
            sage: D = Dgms.an_element(); D
            [(0, 2), (1, 1), (2, 3)]
            sage: D.pp()
            . . O .
            . O . .
            . . . O
        """
        return self([(0, 2), (1, 1), (2, 3)])

    def from_polyomino(self, p):
        r"""
        Create the diagram corresponding to a 2d
        :class:`~sage.combinat.tiling.Polyomino`

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino                            # needs sage.modules
            sage: p = Polyomino([(0,0),(1,0),(1,1),(1,2)])                              # needs sage.modules
            sage: from sage.combinat.diagram import Diagrams
            sage: Diagrams()(p).pp()                                                    # needs sage.modules
            O . .
            O O O

        We can also call this method directly::

            sage: Diagrams().from_polyomino(p).pp()                                     # needs sage.modules
            O . .
            O O O

        This only works for a 2d :class:`~sage.combinat.tiling.Polyomino`::

            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')     # needs sage.modules
            sage: Diagrams().from_polyomino(p)                                          # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: the polyomino must be 2 dimensional
        """
        if not p._dimension == 2:
            raise ValueError("the polyomino must be 2 dimensional")
        cells = list(map(tuple, p))
        return self.element_class(self, cells)

    def from_composition(self, alpha):
        r"""
        Create the diagram corresponding to a weak composition `\alpha \vDash n`.

        EXAMPLES::

            sage: alpha = Composition([3,0,2,1,4,4])
            sage: from sage.combinat.diagram import Diagrams
            sage: Diagrams()(alpha).pp()
            O O O .
            . . . .
            O O . .
            O . . .
            O O O O
            O O O O
            sage: Diagrams().from_composition(alpha).pp()
            O O O .
            . . . .
            O O . .
            O . . .
            O O O O
            O O O O
        """
        cells = []
        for i, n in enumerate(alpha):
            cells.extend((i, j) for j in range(n))
        return self.element_class(self, cells, check=False)

    def from_zero_one_matrix(self, M, check=True):
        r"""
        Get a diagram from a matrix with entries in `\{0, 1\}`, where
        positions of cells are indicated by the `1`'s.

        EXAMPLES::

            sage: M = matrix([[1,0,1,1],[0,1,1,0]])                                     # needs sage.modules
            sage: from sage.combinat.diagram import Diagrams
            sage: Diagrams()(M).pp()                                                    # needs sage.modules
            O . O O
            . O O .
            sage: Diagrams().from_zero_one_matrix(M).pp()                               # needs sage.modules
            O . O O
            . O O .

            sage: M = matrix([[1, 0, 0], [1, 0, 0], [0, 0, 0]])                         # needs sage.modules
            sage: Diagrams()(M).pp()                                                    # needs sage.modules
            O . .
            O . .
            . . .
        """
        # check matrix is zero-one
        n_rows, n_cols = M.dimensions()

        if check:
            zero = M.base_ring().zero()
            one = M.base_ring().one()
            for i in range(n_rows):
                for j in range(n_cols):
                    if not (M[i, j] == zero or M[i, j] == one):
                        raise ValueError("matrix entries must be 0 or 1")
        cells = [(i, j) for i in range(n_rows) for j in range(n_cols) if M[i, j]]

        return self.element_class(self, cells, n_rows, n_cols, check=False)

    Element = Diagram


####################
# Northwest diagrams
####################

class NorthwestDiagram(Diagram, metaclass=InheritComparisonClasscallMetaclass):
    r"""
    Diagrams with the northwest property.

    A diagram is a set of cells indexed by natural numbers. Such a diagram
    has the *northwest property* if the presence of cells `(i1, j1)` and
    `(i2, j2)` implies the presence of the cell
    `(\min(i1, i2), \min(j1, j2))`. Diagrams with the northwest property are
    called *northwest diagrams*.

    For general diagrams see :class:`Diagram`.

    EXAMPLES::

        sage: from sage.combinat.diagram import NorthwestDiagram
        sage: N = NorthwestDiagram([(0,0), (0, 2), (2,0)])

    To visualize them, use the ``.pp()`` method::

        sage: N.pp()
        O . O
        . . .
        O . .
    """
    @staticmethod
    def __classcall_private__(self, cells, n_rows=None, n_cols=None, check=True):
        """
        Normalize input to ensure a correct parent. This method also allows
        one to specify whether or not to check the northwest property for the
        provided cells.

        EXAMPLES::

            sage: from sage.combinat.diagram import NorthwestDiagram, NorthwestDiagrams
            sage: N1 = NorthwestDiagram([(0,1), (0,2)])
            sage: N2 = NorthwestDiagram([(0,1), (0,3)])
            sage: N1.parent() is N2.parent()
            True
            sage: N3 = NorthwestDiagrams()([(0,1), (0,2)])
            sage: N3.parent() is NorthwestDiagrams()
            True
            sage: N1.parent() is NorthwestDiagrams()
            True
        """
        return NorthwestDiagrams()(cells, n_rows, n_cols, check)

    def check(self):
        r"""
        A diagram has the northwest property if the presence of cells
        `(i1, j1)` and `(i2, j2)` implies the presence of the cell
        `(min(i1, i2), min(j1, j2))`. This method checks if the northwest
        property is satisfied for ``self``

        EXAMPLES::

            sage: from sage.combinat.diagram import NorthwestDiagram
            sage: N = NorthwestDiagram([(0,0), (0,3), (3,0)])
            sage: N.check()

        Here is a non-example::

            sage: notN = NorthwestDiagram([(0,1), (1,0)])  #.check() is implicit
            Traceback (most recent call last):
            ...
            ValueError: diagram is not northwest

        TESTS::

            sage: NorthwestDiagram([(0,1/2)])
            Traceback (most recent call last):
            ...
            ValueError: diagrams must be indexed by nonnegative integers
        """
        from itertools import combinations
        Diagram.check(self)
        if not all((min(i1, i2), min(j1, j2)) in self
                   for (i1, j1), (i2, j2) in combinations(self._cells, 2)):
            raise ValueError("diagram is not northwest")

    def peelable_tableaux(self):
        r"""
        Return the set of peelable tableaux whose diagram is ``self``.

        For a fixed northwest diagram `D`, we say that a Young tableau `T` is
        `D`-peelable if:

        1. the row indices of the cells in the first column of `D` are
           the entries in an initial segment in the first column of `T` and
        2. the tableau `Q` obtained by removing those cells from `T` and playing
           jeu de taquin is `(D-C)`-peelable, where `D-C` is the diagram formed
           by forgetting the first column of `D`.

        Reiner and Shimozono [RS1995]_ showed that the number
        `\operatorname{red}(w)` of reduced words of a permutation `w` may be
        computed using the peelable tableaux of the Rothe diagram `D(w)`.
        Explicitly,

        .. MATH::

            \operatorname{red}(w) = \sum_{T} f_{\operatorname{shape} T},

        where the sum runs over the `D(w)`-peelable tableaux `T` and `f_\lambda`
        is the number of standard Young tableaux of shape `\lambda` (which may
        be computed using the hook-length formula).

        EXAMPLES:

        We can compute the `D`-peelable diagrams for a northwest diagram `D`::

            sage: from sage.combinat.diagram import NorthwestDiagram
            sage: cells = [(0,0), (0,1), (0,2), (1,0), (2,0), (2,2), (2,4),
            ....:          (4,0), (4,2)]
            sage: D = NorthwestDiagram(cells); D.pp()
            O O O . .
            O . . . .
            O . O . O
            . . . . .
            O . O . .
            sage: D.peelable_tableaux()
            {[[1, 1, 1], [2, 3, 3], [3, 5], [5]],
            [[1, 1, 1, 3], [2, 3], [3, 5], [5]]}

        EXAMPLES:

        If the diagram is only one column, there is only one peelable tableau::

            sage: from sage.combinat.diagram import NorthwestDiagram
            sage: NWD = NorthwestDiagram([(0,0), (2,0)])
            sage: NWD.peelable_tableaux()
            {[[1], [3]]}

        From [RS1995]_, we know that there is only one peelable tableau for the
        Rothe diagram of the permutation (in one line notation) `251643`::

            sage: D = NorthwestDiagram([(1, 2), (1, 3), (3, 2), (3, 3), (4, 2)])
            sage: D.pp()
            . . . .
            . . O O
            . . . .
            . . O O
            . . O .

            sage: D.peelable_tableaux()
            {[[2, 2], [4, 4], [5]]}

        Here are all the intermediate steps to compute the peelables for the
        Rothe diagram of (in one-line notation) `64817235`. They are listed from
        deepest in the recursion to the final step. The recursion has depth five
        in this case so we will label the intermediate tableaux by `D_i` where
        `i` is the step in the recursion at which they appear.

        Start with the one that has a single column::

            sage: D5 = NorthwestDiagram([(2,0)]); D5.pp()
            .
            .
            O
            sage: D5.peelable_tableaux()
            {[[3]]}

        Now we know all of the `D_5` peelables, so we can compute the `D_4`
        peelables::

            sage: D4 = NorthwestDiagram([(0, 0), (2,0), (4, 0), (2, 2)])
            sage: D4.pp()
            O . .
            . . .
            O . O
            . . .
            O . .

            sage: D4.peelable_tableaux()
            {[[1, 3], [3], [5]]}

        There is only one `D_4` peelable, so we can compute the `D_3`
        peelables::

            sage: D3 = NorthwestDiagram([(0,0), (0,1), (2, 1), (2, 3), (4,1)])
            sage: D3.pp()
            O O . .
            . . . .
            . O . O
            . . . .
            . O . .

            sage: D3.peelable_tableaux()
            {[[1, 1], [3, 3], [5]], [[1, 1, 3], [3], [5]]}

        Now compute the `D_2` peelables::

            sage: cells = [(0,0), (0,1), (0,2), (1,0), (2,0), (2,2), (2,4),
            ....:          (4,0), (4,2)]
            sage: D2 = NorthwestDiagram(cells); D2.pp()
            O O O . .
            O . . . .
            O . O . O
            . . . . .
            O . O . .

            sage: D2.peelable_tableaux()
            {[[1, 1, 1], [2, 3, 3], [3, 5], [5]],
            [[1, 1, 1, 3], [2, 3], [3, 5], [5]]}

        And the `D_1` peelables::

            sage: cells = [(0,0), (0,1), (0,2), (0,3), (1,0), (1,1), (2,0),
            ....:          (2,1), (2,3), (2,5), (4,0), (4,1), (4,3)]
            sage: D1 = NorthwestDiagram(cells); D1.pp()
            O O O O . .
            O O . . . .
            O O . O . O
            . . . . . .
            O O . O . .

            sage: D1.peelable_tableaux()
            {[[1, 1, 1, 1], [2, 2, 3, 3], [3, 3, 5], [5, 5]],
            [[1, 1, 1, 1, 3], [2, 2, 3], [3, 3, 5], [5, 5]]}

        Which we can use to get the `D` peelables::

            sage: cells = [(0,0), (0,1), (0,2), (0,3), (0,4),
            ....:          (1,0), (1,1), (1,2),
            ....:          (2,0), (2,1), (2,2), (2,4), (2,6),
            ....:                 (4,1), (4,2), (4,4)]
            sage: D = NorthwestDiagram(cells); D.pp()
            O O O O O . .
            O O O . . . .
            O O O . O . O
            . . . . . . .
            . O O . O . .
            sage: D.peelable_tableaux()
            {[[1, 1, 1, 1, 1], [2, 2, 2, 3, 3], [3, 3, 3], [5, 5, 5]],
             [[1, 1, 1, 1, 1], [2, 2, 2, 3, 3], [3, 3, 3, 5], [5, 5]],
             [[1, 1, 1, 1, 1, 3], [2, 2, 2, 3], [3, 3, 3], [5, 5, 5]],
             [[1, 1, 1, 1, 1, 3], [2, 2, 2, 3], [3, 3, 3, 5], [5, 5]]}

        ALGORITHM:

        This implementation uses the algorithm suggested in Remark 25
        of [RS1995]_.

        TESTS:

        Corner case::

            sage: from sage.combinat.diagram import NorthwestDiagram
            sage: D = NorthwestDiagram([])
            sage: D.peelable_tableaux()
            {[]}
        """
        # TODO: There is a condition on the first column (if the rows in Dhat
        # are a subset of the rows in the first column) which simplifies the
        # description without performing JDT, so we should implement that

        # empty diagram case
        if not self:
            return set([Tableau([])])

        # if there is a single column in the diagram then there is only
        # one posslbe peelable tableau.
        if self._n_nonempty_cols == 1:
            return set([Tableau([[i+1] for i, j in self.cells()])])

        first_col = min(j for i, j in self._cells)

        dhat_cells = []
        new_vals_cells = []
        for i, j in self._cells:
            if j != first_col:
                dhat_cells.append((i, j))
            else:
                new_vals_cells.append(i + 1)

        new_vals = sorted(new_vals_cells)

        Dhat = NorthwestDiagram(dhat_cells)
        k = self.n_cells() - Dhat.n_cells()

        peelables = set()

        for Q in Dhat.peelable_tableaux():
            # get the vertical strips
            mu = Q.shape()
            vertical_strip_cells = mu.vertical_border_strip_cells(k)
            for s in vertical_strip_cells:
                sQ = SkewTableaux()(Q)  # sQ is skew - get it?
                # perform the jeu de taquin slides
                for c in s:
                    sQ = sQ.backward_slide(c)
                # create the new tableau by filling the columns
                sQ_new = sQ.to_list()
                for n in range(k):
                    sQ_new[n][0] = new_vals[n]

                T = Tableau(sQ_new)
                if T.is_column_strict():
                    peelables.add(T)

        return peelables


class NorthwestDiagrams(Diagrams):
    r"""
    Diagrams satisfying the northwest property.

    A diagram `D` is a *northwest diagram* if for every two cells `(i_1, j_1)`
    and `(i_2, j_2)` in `D` then there exists the cell
    `(\min(i_1, i_2), \min(j_1, j_2)) \in D`.

    EXAMPLES::

        sage: from sage.combinat.diagram import NorthwestDiagram
        sage: N = NorthwestDiagram([(0,0), (0, 10), (5,0)]); N.pp()
        O . . . . . . . . . O
        . . . . . . . . . . .
        . . . . . . . . . . .
        . . . . . . . . . . .
        . . . . . . . . . . .
        O . . . . . . . . . .

    Note that checking whether or not the northwest property is satisfied is
    automatically checked. The diagram found by adding the cell `(1,1)` to the
    diagram above is *not* a northwest diagram. The cell `(1,0)` should be
    present due to the presence of `(5,0)` and `(1,1)`::

        sage: from sage.combinat.diagram import Diagram
        sage: Diagram([(0, 0), (0, 10), (5, 0), (1, 1)]).pp()
        O . . . . . . . . . O
        . O . . . . . . . . .
        . . . . . . . . . . .
        . . . . . . . . . . .
        . . . . . . . . . . .
        O . . . . . . . . . .
        sage: NorthwestDiagram([(0, 0), (0, 10), (5, 0), (1, 1)])
        Traceback (most recent call last):
        ...
        ValueError: diagram is not northwest

    However, this behavior can be turned off if you are confident that
    you are providing a northwest diagram::

        sage: N = NorthwestDiagram([(0, 0), (0, 10), (5, 0),
        ....:                      (1, 1), (0, 1), (1, 0)],
        ....:                      check=False)
        sage: N.pp()
        O O . . . . . . . . O
        O O . . . . . . . . .
        . . . . . . . . . . .
        . . . . . . . . . . .
        . . . . . . . . . . .
        O . . . . . . . . . .

    Note that arbitrary diagrams which happen to be northwest diagrams
    only live in the parent of :class:`Diagrams`::

        sage: D = Diagram([(0, 0), (0, 10), (5, 0), (1, 1), (0, 1), (1, 0)])
        sage: D.pp()
        O O . . . . . . . . O
        O O . . . . . . . . .
        . . . . . . . . . . .
        . . . . . . . . . . .
        . . . . . . . . . . .
        O . . . . . . . . . .
        sage: from sage.combinat.diagram import NorthwestDiagrams
        sage: D in NorthwestDiagrams()
        False

    Here are some more examples::

        sage: from sage.combinat.diagram import NorthwestDiagram, NorthwestDiagrams
        sage: D = NorthwestDiagram([(0,1), (0,2), (1,1)]); D.pp()
        . O O
        . O .
        sage: NWDgms = NorthwestDiagrams()
        sage: D = NWDgms([(1,1), (1,2), (2,1)]); D.pp()
        . . .
        . O O
        . O .
        sage: D.parent()
        Combinatorial northwest diagrams

    Additionally, there are natural constructions of a northwest diagram
    given the data of a permutation (Rothe diagrams are the prototypical example
    of northwest diagrams), or the data of a partition of an integer, or a
    skew partition.

    The Rothe diagram `D(\omega)` of a permutation `\omega` is specified by
    the cells

    .. MATH::

        D(\omega) = \{(\omega_j, i) : i<j,\, \omega_i > \omega_j \}.

    We can construct one by calling :meth:`rothe_diagram` method on the set
    of all :class:`~sage.combinat.diagram.NorthwestDiagrams`::

        sage: w = Permutations(4)([4,3,2,1])
        sage: NorthwestDiagrams().rothe_diagram(w).pp()
        O O O .
        O O . .
        O . . .
        . . . .

    To turn a Ferrers diagram into a northwest diagram, we may call
    :meth:`from_partition`. This will return a Ferrer's diagram in the
    set of all northwest diagrams. For many use-cases it is probably better
    to get Ferrer's diagrams by the corresponding method on partitions, namely
    :meth:`sage.combinat.partitions.Partitions.ferrers_diagram`::

        sage: mu = Partition([7,3,1,1])
        sage: mu.pp()
        *******
        ***
        *
        *
        sage: NorthwestDiagrams().from_partition(mu).pp()
        O O O O O O O
        O O O . . . .
        O . . . . . .
        O . . . . . .

    It is also possible to turn a Ferrers diagram of a skew partition into a
    northwest diagram, altough it is more subtle than just using the skew
    diagram itself. One must first reflect the partition about a vertical axis
    so that the skew partition looks "backwards"::

        sage: mu, nu = Partition([5,4,3,2,1]), Partition([3,2,1])
        sage: s = mu/nu; s.pp()
           **
          **
         **
        **
        *
        sage: NorthwestDiagrams().from_skew_partition(s).pp()
        O O . . .
        . O O . .
        . . O O .
        . . . O O
        . . . . O
    """

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import NorthwestDiagrams
            sage: NWDgms = NorthwestDiagrams(); NWDgms
            Combinatorial northwest diagrams
        """
        return 'Combinatorial northwest diagrams'

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.diagram import NorthwestDiagrams
            sage: NWDgms = NorthwestDiagrams()
            sage: NWD = NWDgms.an_element(); NWD
            [(0, 1), (0, 2), (1, 1), (2, 3)]
            sage: NWD.pp()
            . O O .
            . O . .
            . . . O
            sage: NWD.parent() is NWDgms
            True
        """
        return self([(0, 1), (0, 2), (1, 1), (2, 3)])

    def rothe_diagram(self, w):
        r"""
        Return the Rothe diagram of ``w``.

        We construct a northwest diagram from a permutation by
        constructing its Rothe diagram. Formally, if `\omega` is
        a :class:`~sage.combinat.permutation.Permutation`
        then the Rothe diagram `D(\omega)` is the diagram whose cells are

        .. MATH::

            D(\omega) = \{(\omega_j, i) : i<j,\, \omega_i > \omega_j \}.

        Informally, one can construct the Rothe diagram by starting with all
        `n^2` possible cells, and then deleting the cells `(i, \omega(i))` as
        well as all cells to the right and below. (These are sometimes called
        "death rays".)

        .. SEEALSO::

            :func:`~sage.combinat.diagram.RotheDiagram`

        EXAMPLES::

            sage: from sage.combinat.diagram import NorthwestDiagrams
            sage: w = Permutations(3)([2,1,3])
            sage: NorthwestDiagrams().rothe_diagram(w).pp()
            O . .
            . . .
            . . .
            sage: NorthwestDiagrams().from_permutation(w).pp()
            O . .
            . . .
            . . .

            sage: w = Permutations(8)([2,5,4,1,3,6,7,8])
            sage: NorthwestDiagrams().rothe_diagram(w).pp()
            O . . . . . . .
            O . O O . . . .
            O . O . . . . .
            . . . . . . . .
            . . . . . . . .
            . . . . . . . .
            . . . . . . . .
            . . . . . . . .
        """
        return RotheDiagram(w)

    from_permutation = rothe_diagram

    def from_partition(self, mu):
        r"""
        Return the Ferrer's diagram of ``mu`` as a northwest diagram.

        EXAMPLES::

            sage: mu = Partition([5,2,1]); mu.pp()
            *****
            **
            *
            sage: mu.parent()
            Partitions
            sage: from sage.combinat.diagram import NorthwestDiagrams
            sage: D = NorthwestDiagrams().from_partition(mu)
            sage: D.pp()
            O O O O O
            O O . . .
            O . . . .
            sage: D.parent()
            Combinatorial northwest diagrams

        This will print in English notation even if the notation is set to
        French for the partition::

            sage: Partitions.options.convention="french"
            sage: mu.pp()
            *
            **
            *****
            sage: D.pp()
            O O O O O
            O O . . .
            O . . . .

        TESTS::

            sage: from sage.combinat.diagram import NorthwestDiagrams
            sage: mu = [5, 2, 1]
            sage: D = NorthwestDiagrams().from_partition(mu)
            Traceback (most recent call last):
            ...
            ValueError: mu must be a Partition
        """
        if not isinstance(mu, Partition):
            raise ValueError("mu must be a Partition")
        return self.element_class(self, mu.cells(), check=False)

    def from_skew_partition(self, s):
        r"""
        Get the northwest diagram found by reflecting a skew shape across
        a vertical plane.

        EXAMPLES::

            sage: mu, nu = Partition([3,2,1]), Partition([2,1])
            sage: s = mu/nu; s.pp()
              *
             *
            *
            sage: from sage.combinat.diagram import NorthwestDiagrams
            sage: D = NorthwestDiagrams().from_skew_partition(s)
            sage: D.pp()
            O . .
            . O .
            . . O

            sage: mu, nu = Partition([3,3,2]), Partition([2,2,2])
            sage: s = mu/nu; s.pp()
              *
              *
            sage: NorthwestDiagrams().from_skew_partition(s).pp()
            O . .
            O . .
            . . .

        TESTS::

            sage: mu = Partition([3,2,1])
            sage: NorthwestDiagrams().from_skew_partition(mu)
            Traceback (most recent call last):
            ...
            ValueError: mu must be a SkewPartition
        """
        if not isinstance(s, SkewPartition):
            raise ValueError("mu must be a SkewPartition")

        n_cols = s.outer()[0]
        n_rows = len(s.outer())

        cells = [(i, n_cols - 1 - j) for i, j in s.cells()]

        return self.element_class(self, cells, n_rows, n_cols, check=False)

    def from_parallelogram_polyomino(self, p):
        r"""
        Create the diagram corresponding to a
        :class:`~sage.combinat.parallelogram_polyomino.ParallelogramPolyomino`.

        EXAMPLES::

            sage: p = ParallelogramPolyomino([[0, 0, 1, 0, 0, 0, 1, 1],                 # needs sage.modules
            ....:                              [1, 1, 0, 1, 0, 0, 0, 0]])
            sage: from sage.combinat.diagram import NorthwestDiagrams
            sage: NorthwestDiagrams().from_parallelogram_polyomino(p).pp()              # needs sage.modules
            O O .
            O O O
            . O O
            . O O
            . O O
        """
        from sage.matrix.constructor import Matrix
        M = Matrix(p.get_array())
        return self.from_zero_one_matrix(M)

    Element = NorthwestDiagram


def RotheDiagram(w):
    r"""
    The Rothe diagram of a permutation ``w``.

    EXAMPLES::

        sage: w = Permutations(9)([1, 7, 4, 5, 9, 3, 2, 8, 6])
        sage: from sage.combinat.diagram import RotheDiagram
        sage: D = RotheDiagram(w); D.pp()
        . . . . . . . . .
        . O O O O O . . .
        . O O . . . . . .
        . O O . . . . . .
        . O O . . O . O .
        . O . . . . . . .
        . . . . . . . . .
        . . . . . O . . .
        . . . . . . . . .

    The Rothe diagram is a northwest diagram::

        sage: D.parent()
        Combinatorial northwest diagrams

    Some other examples::

        sage: RotheDiagram([2, 1, 4, 3]).pp()
        O . . .
        . . . .
        . . O .
        . . . .

        sage: RotheDiagram([4, 1, 3, 2]).pp()
        O O O .
        . . . .
        . O . .
        . . . .

    Currently, only elements of the set of
    :class:`sage.combinat.permutations.Permutations` are supported. In
    particular, elements of permutation groups are not supported::

        sage: w = SymmetricGroup(9).an_element()                                        # needs sage.groups
        sage: RotheDiagram(w)                                                           # needs sage.groups
        Traceback (most recent call last):
        ...
        ValueError: w must be a permutation

    TESTS::

        sage: w = Permutations(5)([1,2,3,4,5])
        sage: from sage.combinat.diagram import RotheDiagram
        sage: RotheDiagram(w).pp()
        . . . . .
        . . . . .
        . . . . .
        . . . . .
        . . . . .
    """
    P = Permutations()
    if w not in P:
        raise ValueError('w must be a permutation')
    w = P(w)

    N = w.size()
    winv = w.inverse()
    cells = [c for c in product(range(N), range(N))
             if c[0] + 1 < winv(c[1] + 1) and c[1] + 1 < w(c[0] + 1)]

    return NorthwestDiagram(cells, n_rows=N, n_cols=N, check=False)
