r"""
Tiling Solver

Tiling a n-dimensional box into non-intersecting n-dimensional polyominoes.

This uses dancing links code which is in Sage.  Dancing links were
originally introduced by Donald Knuth in 2000 [Knuth1]_. In
particular, Knuth used dancing links to solve tilings of a region by
2d pentaminoes.  Here we extend the method to any dimension.

In particular, the :mod:`sage.games.quantumino` module is based on
the Tiling Solver and allows to solve the 3d Quantumino puzzle.

This module defines two classes:

- :class:`sage.combinat.tiling.Polyomino` class, to represent polyominoes
  in arbitrary dimension. The goal of this class is to return all the
  rotated, reflected and/or translated copies of a polyomino that are
  contained in a certain box.

- :class:`sage.combinat.tiling.TilingSolver` class, to solve the general
  problem of tiling a rectangular `n`-dimensional box with a set of
  `n`-dimensional polyominoes. One can specify if rotations and reflections
  are allowed or not and if pieces can be reused or not. This class convert
  the tiling data into rows of a matrix that are passed to the DLX solver.
  It also allows to compute the number of solutions.

AUTHOR:

    - Sebastien Labbe, June 2011

EXAMPLES:

2d Easy Example
---------------

Here is a 2d example. Let us try to fill the `3 \times 2` rectangle with a
`1 \times 2` rectangle and a `2 \times 2` square. Obviously, there are two
solutions::

    sage: from sage.combinat.tiling import TilingSolver, Polyomino
    sage: p = Polyomino([(0,0), (0,1)])
    sage: q = Polyomino([(0,0), (0,1), (1,0), (1,1)])
    sage: T = TilingSolver([p,q], box=[3,2])
    sage: it = T.solve()
    sage: next(it)
    [Polyomino: [(0, 0), (0, 1), (1, 0), (1, 1)], Color: gray, Polyomino: [(2, 0), (2, 1)], Color: gray]
    sage: next(it)
    [Polyomino: [(1, 0), (1, 1), (2, 0), (2, 1)], Color: gray, Polyomino: [(0, 0), (0, 1)], Color: gray]
    sage: next(it)
    Traceback (most recent call last):
    ...
    StopIteration
    sage: T.number_of_solutions()
    2

1d Easy Example
---------------

Here is an easy one dimensional example where we try to tile a stick of
length 6 with three sticks of length 1, 2 and 3. There are six solutions::

    sage: p = Polyomino([[0]])
    sage: q = Polyomino([[0],[1]])
    sage: r = Polyomino([[0],[1],[2]])
    sage: T = TilingSolver([p,q,r], box=[6])
    sage: len(T.rows())
    15
    sage: it = T.solve()
    sage: next(it)
    [Polyomino: [(0,)], Color: gray, Polyomino: [(1,), (2,)], Color: gray, Polyomino: [(3,), (4,), (5,)], Color: gray]
    sage: next(it)
    [Polyomino: [(0,)], Color: gray, Polyomino: [(1,), (2,), (3,)], Color: gray, Polyomino: [(4,), (5,)], Color: gray]
    sage: T.number_of_solutions()
    6

2d Puzzle allowing reflections
------------------------------

The following is a puzzle owned by Florent Hivert::

    sage: from sage.combinat.tiling import Polyomino, TilingSolver
    sage: L = []
    sage: L.append(Polyomino([(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,3)], 'yellow'))
    sage: L.append(Polyomino([(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2)], "black"))
    sage: L.append(Polyomino([(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,3)], "gray"))
    sage: L.append(Polyomino([(0,0),(0,1),(0,2),(0,3),(1,0),(1,3)],"cyan"))
    sage: L.append(Polyomino([(0,0),(0,1),(0,2),(0,3),(1,0),(1,1)],"red"))
    sage: L.append(Polyomino([(0,0),(0,1),(0,2),(0,3),(1,1),(1,2)],"blue"))
    sage: L.append(Polyomino([(0,0),(0,1),(0,2),(0,3),(1,1),(1,3)],"green"))
    sage: L.append(Polyomino([(0,1),(0,2),(0,3),(1,0),(1,1),(1,3)],"magenta"))
    sage: L.append(Polyomino([(0,1),(0,2),(0,3),(1,0),(1,1),(1,2)],"orange"))
    sage: L.append(Polyomino([(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)],"pink"))

By default, rotations are allowed and reflections are not. In this case,
there are no solution::

    sage: T = TilingSolver(L, (8,8))
    sage: T.number_of_solutions()                       # long time (2.5 s)
    0

If reflections are allowed, there are solutions. Solve the puzzle and show
one solution::

    sage: T = TilingSolver(L, (8,8), reflection=True)
    sage: solution = next(T.solve())
    sage: G = sum([piece.show2d() for piece in solution], Graphics())
    sage: G.show(aspect_ratio=1, axes=False)

Compute the number of solutions::

    sage: T.number_of_solutions()                         # long time (2.6s)
    328

Create a animation of all the solutions::

    sage: a = T.animate()                            # not tested
    sage: a                                          # not tested
    Animation with 328 frames

3d Puzzle
---------

The same thing done in 3d *without* allowing reflections this time::

    sage: from sage.combinat.tiling import Polyomino, TilingSolver
    sage: L = []
    sage: L.append(Polyomino([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0),(1,2,0),(1,3,0)]))
    sage: L.append(Polyomino([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0),(1,2,0)]))
    sage: L.append(Polyomino([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0),(1,3,0)]))
    sage: L.append(Polyomino([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,3,0)]))
    sage: L.append(Polyomino([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0)]))
    sage: L.append(Polyomino([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,1,0),(1,2,0)]))
    sage: L.append(Polyomino([(0,0,0),(0,1,0),(0,2,0),(0,3,0),(1,1,0),(1,3,0)]))
    sage: L.append(Polyomino([(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0),(1,3,0)]))
    sage: L.append(Polyomino([(0,1,0),(0,2,0),(0,3,0),(1,0,0),(1,1,0),(1,2,0)]))
    sage: L.append(Polyomino([(0,0,0),(0,1,0),(0,2,0),(1,0,0),(1,1,0),(1,2,0)]))

Solve the puzzle and show one solution::

    sage: T = TilingSolver(L, (8,8,1))
    sage: solution = next(T.solve())
    sage: G = sum([piece.show3d(size=0.85) for piece in solution], Graphics())
    sage: G.show(aspect_ratio=1, viewer='tachyon')

Let us compute the number of solutions::

    sage: T.number_of_solutions()                              # long time (3s)
    328

Donald Knuth example : the Y pentamino
--------------------------------------

Donald Knuth [Knuth1]_ considered the problem of packing 45 Y pentaminoes into a
`15 \times 15` square::

    sage: from sage.combinat.tiling import Polyomino, TilingSolver
    sage: y = Polyomino([(0,0),(1,0),(2,0),(3,0),(2,1)])
    sage: T = TilingSolver([y], box=(5,10), reusable=True, reflection=True)
    sage: T.number_of_solutions()
    10
    sage: solution = next(T.solve())
    sage: G = sum([piece.show2d() for piece in solution], Graphics())
    sage: G.show(aspect_ratio=1)

::

    sage: T = TilingSolver([y], box=(15,15), reusable=True, reflection=True)
    sage: T.number_of_solutions()                      #not tested
    212

Animation of Donald Knuth's dancing links
-----------------------------------------

Animation of the solutions::

    sage: from sage.combinat.tiling import Polyomino, TilingSolver
    sage: Y = Polyomino([(0,0),(1,0),(2,0),(3,0),(2,1)], color='yellow')
    sage: T = TilingSolver([Y], box=(15,15), reusable=True, reflection=True)
    sage: a = T.animate(stop=40)            # long time  # optional -- ImageMagick
    sage: a                                 # long time  # optional -- ImageMagick
    Animation with 40 frames

Incremental animation of the solutions (one piece is removed/added at a time)::

    sage: a = T.animate('incremental', stop=40)   # long time  # optional -- ImageMagick
    sage: a                                       # long time  # optional -- ImageMagick
    Animation with 40 frames
    sage: a.show(delay=50, iterations=1)          # long time  # optional -- ImageMagick

5d Easy Example
---------------

Here is a 5d example. Let us try to fill the `2 \times 2 \times 2 \times 2
\times 2` rectangle with reusable `1 \times 1 \times 1 \times 1 \times 1`
rectangles. Obviously, there is one solution::

    sage: from sage.combinat.tiling import Polyomino, TilingSolver
    sage: p = Polyomino([(0,0,0,0,0)])
    sage: T = TilingSolver([p], box=(2,2,2,2,2), reusable=True)
    sage: rows = T.rows()                               # long time (3s)
    sage: rows                                          # long time (fast)
    [[0], [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12], [13], [14], [15], [16], [17], [18], [19], [20], [21], [22], [23], [24], [25], [26], [27], [28], [29], [30], [31]]
    sage: T.number_of_solutions()                       # long time (fast)
    1

REFERENCES:

.. [Knuth1] Knuth, Donald (2000). "Dancing links". :arxiv:`cs/0011047`.

"""
#*****************************************************************************
#       Copyright (C) 2011 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import itertools
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method, cached_function
from sage.misc.all import prod
from sage.combinat.all import WeylGroup
from sage.plot.all import Graphics
from sage.plot.polygon import polygon
from sage.plot.line import line
from sage.plot.circle import circle
from sage.modules.free_module_element import vector
from sage.plot.plot3d.platonic import cube
from sage.plot.animate import Animation
from sage.misc.mrange import xmrange

############################
# Orthogonal transformations
############################


@cached_function
def orthogonal_transformation(n, orientation_preserving=True):
    r"""
    Return the list of orthogonal transformation matrices in the
    `n`-dimensional vector space.

    INPUT:

    - ``n`` - positive integer, dimension of the space
    - ``orientation_preserving`` - bool (optional, default: ``True``),
      whether the orientation is preserved

    OUTPUT:

        list of matrices

    EXAMPLES::

        sage: from sage.combinat.tiling import orthogonal_transformation
        sage: orthogonal_transformation(2)
        [
        [1 0]  [ 0  1]  [-1  0]  [ 0 -1]
        [0 1], [-1  0], [ 0 -1], [ 1  0]
        ]
        sage: orthogonal_transformation(2, orientation_preserving=False)
        [
        [1 0]  [ 0 -1]  [ 1  0]  [ 0  1]  [0 1]  [-1  0]  [ 0 -1]  [-1  0]
        [0 1], [-1  0], [ 0 -1], [-1  0], [1 0], [ 0 -1], [ 1  0], [ 0  1]
        ]
        sage: orthogonal_transformation(3)
        [
        [1 0 0]  [ 1  0  0]  [-1  0  0]  [-1  0  0]  [0 0 1]  [ 0  0 -1]
        [0 1 0]  [ 0 -1  0]  [ 0  1  0]  [ 0 -1  0]  [1 0 0]  [ 1  0  0]
        [0 0 1], [ 0  0 -1], [ 0  0 -1], [ 0  0  1], [0 1 0], [ 0 -1  0],
        <BLANKLINE>
        [ 0  0 -1]  [ 0  0  1]  [0 1 0]  [ 0 -1  0]  [ 0  1  0]  [ 0 -1  0]
        [-1  0  0]  [-1  0  0]  [0 0 1]  [ 0  0 -1]  [ 0  0 -1]  [ 0  0  1]
        [ 0  1  0], [ 0 -1  0], [1 0 0], [ 1  0  0], [-1  0  0], [-1  0  0],
        <BLANKLINE>
        [ 0  1  0]  [ 0 -1  0]  [ 0  1  0]  [ 0 -1  0]  [ 1  0  0]  [ 1  0  0]
        [ 1  0  0]  [ 1  0  0]  [-1  0  0]  [-1  0  0]  [ 0  0 -1]  [ 0  0  1]
        [ 0  0 -1], [ 0  0  1], [ 0  0  1], [ 0  0 -1], [ 0  1  0], [ 0 -1  0],
        <BLANKLINE>
        [-1  0  0]  [-1  0  0]  [ 0  0 -1]  [ 0  0  1]  [ 0  0  1]  [ 0  0 -1]
        [ 0  0  1]  [ 0  0 -1]  [ 0  1  0]  [ 0 -1  0]  [ 0  1  0]  [ 0 -1  0]
        [ 0  1  0], [ 0 -1  0], [ 1  0  0], [ 1  0  0], [-1  0  0], [-1  0  0]
        ]

    TESTS::

        sage: orthogonal_transformation(1)
        [[1]]
        sage: orthogonal_transformation(0)
        Traceback (most recent call last):
        ...
        ValueError: ['B', 0] is not a valid Cartan type
    """
    if orientation_preserving:
        return [w.matrix() for w in WeylGroup(['B', n])
                if w.matrix().det() == 1]
    else:
        return [w.matrix() for w in WeylGroup(['B', n])]

##############################
# Class Polyomino
##############################


class Polyomino(SageObject):
    r"""
    Return the polyomino defined by a set of coordinates.

    The polyomino is the union of the unit square (or cube, or n-cube)
    centered at those coordinates. Such an object should be connected, but
    the code do not make this assumption.

    INPUT:

    - ``coords`` - iterable of tuple
    - ``color`` - string (optional, default: ``'gray'``), the color

    EXAMPLES::

        sage: from sage.combinat.tiling import Polyomino
        sage: Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
        Polyomino: [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)], Color: blue
    """
    def __init__(self, coords, color='gray'):
        r"""
        INPUT:

        - ``coords`` - iterable of tuple
        - ``color`` - string (optional, default: ``'gray'``), the color

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            Polyomino: [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)], Color: blue

        ::

            sage: from sage.combinat.tiling import Polyomino
            sage: Polyomino([(0,0), (1,0), (2,0)])
            Polyomino: [(0, 0), (1, 0), (2, 0)], Color: gray
        """
        assert isinstance(color, str)
        self._color = color
        self._blocs = frozenset(tuple(c) for c in coords)
        assert len(self._blocs) != 0, "Polyomino must be non empty"
        dimension_set = set(len(a) for a in self._blocs)
        assert len(dimension_set) <= 1, "coord must be all of the same dimension"
        self._dimension = dimension_set.pop()

    def __repr__(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='red')
            Polyomino: [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)], Color: red
        """
        s = "Polyomino: %s, " % sorted(self._blocs)
        s += "Color: %s" % self._color
        return s

    def __hash__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: hash(p)     # random
            2059134902
        """
        return hash(self._blocs)

    def __eq__(self, other):
        r"""
        Return whether self is equal to other.

        INPUT:

        - ``other`` - a polyomino

        OUTPUT:

            boolean

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: q = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='red')
            sage: p == q
            True
            sage: r = Polyomino([(0,0,0), (0,1,0), (1,1,0)], color='blue')
            sage: p == r
            False
        """
        return self._blocs == other._blocs

    def __ne__(self, other):
        r"""
        Return whether self is not equal to other.

        INPUT:

        - ``other`` - a polyomino

        OUTPUT:

            boolean

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: q = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='red')
            sage: p != q
            False
            sage: r = Polyomino([(0,0,0), (0,1,0), (1,1,0)], color='blue')
            sage: p != r
            True
        """
        return self._blocs != other._blocs

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: it = iter(p)
            sage: next(it)
            (1, 1, 0)
        """
        return iter(self._blocs)

    def color(self):
        r"""
        Return the color of the polyomino.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: p.color()
            'blue'
        """
        return self._color

    def __len__(self):
        r"""
        Return the size of the polyomino, i.e. the number of n-dimensional
        unit cubes.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: len(p)
            4
        """
        return len(self._blocs)

    def orthogonals(self, orientation_preserving=True):
        r"""
        Iterator over the images of self under orthogonal transformations.

        .. NOTE::

            No guarantee of uniqueness.

        INPUT:

        - ``orientation_preserving`` - bool (optional, default: ``True``),
          whether the orientation is preserved

        OUTPUT:

            iterator of Polyomino

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: L = list(p.orthogonals())
            sage: len(L)
            24
            sage: L = list(p.orthogonals(False))
            sage: len(L)
            48
        """
        return (m * self for m in
                orthogonal_transformation(self._dimension,
                                          orientation_preserving))

    def canonical_orthogonals(self, orientation_preserving=True):
        r"""
        Iterator over the image of self under orthogonal transformations
        where the coordinates are all positive and minimal.

        .. NOTE::

            No guarantee of uniqueness.

        INPUT:

        - ``orientation_preserving`` - bool (optional, default: ``True``),
          whether the orientation is preserved

        OUTPUT:

            iterator of Polyomino

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: L = list(p.canonical_orthogonals())
            sage: len(L)
            24

        They might not be all different::

            sage: s = set(p.canonical_orthogonals())
            sage: len(s)
            12

        With the non orientation-preserving::

            sage: s = set(p.canonical_orthogonals(False))
            sage: len(s)
            24
        """
        for q in self.orthogonals(orientation_preserving):
            yield q.canonical()

    def canonical(self):
        r"""
        Returns the translated copy of self having minimal and positive
        coordinates

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: p
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            sage: p.canonical()
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink

        TESTS::

            sage: p
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            sage: p + (3,4,5)
            Polyomino: [(3, 4, 5), (4, 4, 5), (4, 5, 5), (4, 5, 6), (4, 6, 5)], Color: deeppink
            sage: (p + (3,4,5)).canonical()
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
        """
        minxyz, maxxyz = self.bounding_box()
        return self - minxyz

    def __sub__(self, v):
        r"""
        Return a translated copy of self by the opposite of the
        vector v.

        INPUT:

        - ``v`` - tuple

        OUPUT:

            polyomino

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: p - (2,2,2)
            Polyomino: [(-2, -2, -2), (-1, -2, -2), (-1, -1, -2), (-1, -1, -1), (-1, 0, -2)], Color: deeppink
        """
        if not len(v) == self._dimension:
            raise ValueError("Dimension of input vector must match the "
                             "dimension of the polyomino")
        v = vector(v)
        return Polyomino([vector(p)-v for p in self], color=self._color)

    def __add__(self, v):
        r"""
        Return a translated copy of self by the vector v.

        INPUT:

        - ``v`` - tuple

        OUPUT:

            polyomino

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: p + (2,2,2)
            Polyomino: [(2, 2, 2), (3, 2, 2), (3, 3, 2), (3, 3, 3), (3, 4, 2)], Color: deeppink
        """
        if not len(v) == self._dimension:
            raise ValueError("Dimension of input vector must match "
                             "the dimension of the polyomino")
        v = vector(v)
        return Polyomino([vector(p)+v for p in self], color=self._color)

    def __rmul__(self, m):
        r"""
        Return the image of the polyomino under the application of the
        matrix m.

        INPUT:

        - ``m`` - square matrix, matching the dimension of self.

        OUPUT:

            Polyomino

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: m = matrix(3, [1,0,0,0,1,0,0,0,1])
            sage: m * p
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            sage: m = matrix(3, [1,0,0,0,0,-1,0,1,0])
            sage: m * p
            Polyomino: [(0, 0, 0), (1, -1, 1), (1, 0, 0), (1, 0, 1), (1, 0, 2)], Color: deeppink

        TESTS::

            sage: m = matrix(2, [1,0,0,1])
            sage: m * p
            Traceback (most recent call last):
            ...
            ValueError: Dimension of input matrix must match the dimension of the polyomino
        """
        if not m.nrows() == m.ncols() == self._dimension:
            raise ValueError("Dimension of input matrix must match the "
                             "dimension of the polyomino")
        return Polyomino([m * vector(p) for p in self], color=self._color)

    def bounding_box(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: p.bounding_box()
            [[0, 0, 0], [1, 2, 1]]
        """
        zipped_coords = zip(*self)
        return [map(min, zipped_coords), map(max, zipped_coords)]

    def translated(self, box):
        r"""
        Returns an iterator over the translated images of self inside a
        box.

        INPUT:

        - ``box`` - tuple, size of the box

        OUTPUT:

            iterator of 3d polyominoes

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: for t in p.translated(box=(5,8,2)): t
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            Polyomino: [(0, 1, 0), (1, 1, 0), (1, 2, 0), (1, 2, 1), (1, 3, 0)], Color: deeppink
            Polyomino: [(0, 2, 0), (1, 2, 0), (1, 3, 0), (1, 3, 1), (1, 4, 0)], Color: deeppink
            Polyomino: [(0, 3, 0), (1, 3, 0), (1, 4, 0), (1, 4, 1), (1, 5, 0)], Color: deeppink
            Polyomino: [(0, 4, 0), (1, 4, 0), (1, 5, 0), (1, 5, 1), (1, 6, 0)], Color: deeppink
            Polyomino: [(0, 5, 0), (1, 5, 0), (1, 6, 0), (1, 6, 1), (1, 7, 0)], Color: deeppink
            Polyomino: [(1, 0, 0), (2, 0, 0), (2, 1, 0), (2, 1, 1), (2, 2, 0)], Color: deeppink
            Polyomino: [(1, 1, 0), (2, 1, 0), (2, 2, 0), (2, 2, 1), (2, 3, 0)], Color: deeppink
            Polyomino: [(1, 2, 0), (2, 2, 0), (2, 3, 0), (2, 3, 1), (2, 4, 0)], Color: deeppink
            Polyomino: [(1, 3, 0), (2, 3, 0), (2, 4, 0), (2, 4, 1), (2, 5, 0)], Color: deeppink
            Polyomino: [(1, 4, 0), (2, 4, 0), (2, 5, 0), (2, 5, 1), (2, 6, 0)], Color: deeppink
            Polyomino: [(1, 5, 0), (2, 5, 0), (2, 6, 0), (2, 6, 1), (2, 7, 0)], Color: deeppink
            Polyomino: [(2, 0, 0), (3, 0, 0), (3, 1, 0), (3, 1, 1), (3, 2, 0)], Color: deeppink
            Polyomino: [(2, 1, 0), (3, 1, 0), (3, 2, 0), (3, 2, 1), (3, 3, 0)], Color: deeppink
            Polyomino: [(2, 2, 0), (3, 2, 0), (3, 3, 0), (3, 3, 1), (3, 4, 0)], Color: deeppink
            Polyomino: [(2, 3, 0), (3, 3, 0), (3, 4, 0), (3, 4, 1), (3, 5, 0)], Color: deeppink
            Polyomino: [(2, 4, 0), (3, 4, 0), (3, 5, 0), (3, 5, 1), (3, 6, 0)], Color: deeppink
            Polyomino: [(2, 5, 0), (3, 5, 0), (3, 6, 0), (3, 6, 1), (3, 7, 0)], Color: deeppink
            Polyomino: [(3, 0, 0), (4, 0, 0), (4, 1, 0), (4, 1, 1), (4, 2, 0)], Color: deeppink
            Polyomino: [(3, 1, 0), (4, 1, 0), (4, 2, 0), (4, 2, 1), (4, 3, 0)], Color: deeppink
            Polyomino: [(3, 2, 0), (4, 2, 0), (4, 3, 0), (4, 3, 1), (4, 4, 0)], Color: deeppink
            Polyomino: [(3, 3, 0), (4, 3, 0), (4, 4, 0), (4, 4, 1), (4, 5, 0)], Color: deeppink
            Polyomino: [(3, 4, 0), (4, 4, 0), (4, 5, 0), (4, 5, 1), (4, 6, 0)], Color: deeppink
            Polyomino: [(3, 5, 0), (4, 5, 0), (4, 6, 0), (4, 6, 1), (4, 7, 0)], Color: deeppink

        This method is independant of the translation of the polyomino::

            sage: q = Polyomino([(0,0,0), (1,0,0)])
            sage: list(q.translated((2,2,1)))
            [Polyomino: [(0, 0, 0), (1, 0, 0)], Color: gray, Polyomino: [(0, 1, 0), (1, 1, 0)], Color: gray]
            sage: q = Polyomino([(34,7,-9), (35,7,-9)])
            sage: list(q.translated((2,2,1)))
            [Polyomino: [(0, 0, 0), (1, 0, 0)], Color: gray, Polyomino: [(0, 1, 0), (1, 1, 0)], Color: gray]

        Inside smaller boxes::

            sage: list(p.translated(box=(2,2,3)))
            []
            sage: list(p.translated(box=(2,3,2)))
            [Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink]
            sage: list(p.translated(box=(3,2,2)))
            []
            sage: list(p.translated(box=(1,1,1)))
            []
            sage: list(p.translated(box=(1,1,-1)))
            []
        """
        if not len(box) == self._dimension:
            raise ValueError("Dimension of input box must match the "
                             "dimension of the polyomino")
        minxyz, maxxyz = map(vector, self.bounding_box())
        size = maxxyz - minxyz
        cano = self.canonical()
        for v in xmrange(vector(box) - vector(size), tuple):
            yield cano + v

    def translated_orthogonals(self, box, orientation_preserving=True):
        r"""
        Return the translated and rotated of self that lies in the box.

        INPUT:

        - ``box`` - tuple of size three, size of the box
        - ``orientation_preserving`` - bool (optional, default: ``True``),
          whether the orientation is preserved

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: L = list(p.translated_orthogonals(box=(5,8,2)))
            sage: len(L)
            360

        ::

            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,2,0),(1,2,1)], color='orange')
            sage: L = list(p.translated_orthogonals(box=(5,8,2)))
            sage: len(L)
            180

        ::

            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,2,0),(1,2,1)], color='orange')
            sage: L = list(p.translated_orthogonals((5,8,2), False))
            sage: len(L)
            360
        """
        if not len(box) == self._dimension:
            raise ValueError("Dimension of input box must match the "
                             "dimension of the polyomino")
        all_distinct_cano = set(self.canonical_orthogonals(orientation_preserving))
        for cano in all_distinct_cano:
            for t in cano.translated(box=box):
                yield t

    def neighbor_edges(self):
        r"""
        Return an iterator over the pairs of neighbor coordinates of the
        polyomino.

        Two points `P` and `Q` are neighbor if `P - Q` has one coordinate
        equal to `+1` or `-1` and zero everywhere else.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(0,0,1)])
            sage: list(sorted(edge) for edge in p.neighbor_edges())
            [[(0, 0, 0), (0, 0, 1)]]

        In 3d::

            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: L = sorted(sorted(edge) for edge in p.neighbor_edges())
            sage: for a in L: a
            [(0, 0, 0), (1, 0, 0)]
            [(1, 0, 0), (1, 1, 0)]
            [(1, 1, 0), (1, 1, 1)]
            [(1, 1, 0), (1, 2, 0)]

        In 2d::

            sage: p = Polyomino([(0,0),(1,0),(1,1),(1,2)])
            sage: L = sorted(sorted(edge) for edge in p.neighbor_edges())
            sage: for a in L: a
            [(0, 0), (1, 0)]
            [(1, 0), (1, 1)]
            [(1, 1), (1, 2)]
        """
        for P, Q in itertools.combinations(self, 2):
            P, Q = vector(P), vector(Q)
            s = sorted(map(abs, Q-P))
            firsts = s[:-1]
            last = s[-1]
            if last == 1 and all(f == 0 for f in firsts):
                yield P, Q

    def center(self):
        r"""
        Return the center of the polyomino.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(0,0,1)])
            sage: p.center()
            (0, 0, 1/2)

        In 3d::

            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: p.center()
            (4/5, 4/5, 1/5)

        In 2d::

            sage: p = Polyomino([(0,0),(1,0),(1,1),(1,2)])
            sage: p.center()
            (3/4, 3/4)
        """
        return sum(vector(t) for t in self) / len(self)

    def boundary(self):
        r"""
        Return the boundary of a 2d polyomino.

        INPUT:

        - ``self`` - a 2d polyomino

        OUTPUT:

        - list of edges (an edge is a pair of adjacent 2d coordinates)

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0), (1,0), (0,1), (1,1)])
            sage: p.boundary()
            [((0.5, 1.5), (1.5, 1.5)), ((-0.5, -0.5), (0.5, -0.5)), ((0.5, -0.5), (1.5, -0.5)), ((-0.5, 1.5), (0.5, 1.5)), ((-0.5, 0.5), (-0.5, 1.5)), ((-0.5, -0.5), (-0.5, 0.5)), ((1.5, 0.5), (1.5, 1.5)), ((1.5, -0.5), (1.5, 0.5))]
            sage: len(_)
            8
            sage: p = Polyomino([(5,5)])
            sage: p.boundary()
            [((4.5, 5.5), (5.5, 5.5)), ((4.5, 4.5), (5.5, 4.5)), ((4.5, 4.5), (4.5, 5.5)), ((5.5, 4.5), (5.5, 5.5))]
        """
        if self._dimension != 2:
            raise NotImplementedError("The method boundary is currently "
                                      "implemented "
                                      "only for dimension 2")
        from collections import defaultdict
        horizontal = defaultdict(int)
        vertical = defaultdict(int)
        for a in self:
            x, y = a = tuple(a)
            horizontal[a] += 1
            vertical[a] += 1
            horizontal[(x, y+1)] -= 1
            vertical[(x+1, y)] -= 1
        edges = []
        h = 0.5
        for (x, y), coeff in horizontal.iteritems():
            if coeff != 0:
                edges.append(((x-h, y-h), (x+h, y-h)))
        for (x, y), coeff in vertical.iteritems():
            if coeff != 0:
                edges.append(((x-h, y-h), (x-h, y+h)))
        return edges

    def show3d(self, size=1):
        r"""
        Returns a 3d Graphic object representing the polyomino.

        INPUT:

        - ``self`` - a polyomino of dimension 3
        - ``size`` - number (optional, default: ``1``), the size of each
          ``1 \times 1 \times 1`` cube. This does a homothety with respect
          to the center of the polyomino.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: p.show3d()
            Graphics3d Object
        """
        assert self._dimension == 3, "Dimension of the polyomino must be 3."
        G = Graphics()
        for p in self:
            G += cube(p, color=self._color)
        center = self.center()
        G = G.translate(-center)
        G = G.scale(size)
        G = G.translate(center)
        return G

    def show2d(self, size=0.7, color='black', thickness=1):
        r"""
        Returns a 2d Graphic object representing the polyomino.

        INPUT:

        - ``self`` - a polyomino of dimension 2
        - ``size`` - number (optional, default: ``0.7``), the size of each
          square.
        - ``color`` - color (optional, default: ``'black'``), color of
          the boundary line.
        - ``thickness`` - number (optional, default: ``1``), how thick the
          boundary line is.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0),(1,0),(1,1),(1,2)], color='deeppink')
            sage: p.show2d()              # long time (0.5s)
            Graphics object consisting of 17 graphics primitives
        """
        assert self._dimension == 2, "Dimension of the polyomino must be 2."
        h = size / 2.0
        G = Graphics()
        for a, b in self:
            G += circle((a, b), h, fill=True, color=self._color)
        k = h / 2.0
        for P, Q in self.neighbor_edges():
            a, b = (P + Q) / 2.0
            G += polygon([(a-k, b-k), (a+k, b-k), (a+k, b+k), (a-k, b+k),
                          (a-k, b-k)], color=self._color)
        for edge in self.boundary():
            G += line(edge, color=color, thickness=thickness)
        return G

#######################
# General tiling solver
#######################


class TilingSolver(SageObject):
    r"""
    Tiling solver

    Solve the problem of tiling a rectangular box with a certain number
    of pieces, called polyominoes, where each polyomino must be used
    exactly once.

    INPUT:

    - ``pieces`` - iterable of Polyominoes
    - ``box`` - tuple, size of the box
    - ``rotation`` - bool (optional, default: ``True``), whether to allow
      rotations
    - ``reflection`` - bool (optional, default: ``False``), whether to allow
      reflections
    - ``reusable`` - bool (optional, default: ``False``), whether to allow
      the pieces to be reused

    EXAMPLES:

    By default, rotations are allowed and reflections are not allowed::

        sage: from sage.combinat.tiling import TilingSolver, Polyomino
        sage: p = Polyomino([(0,0,0)])
        sage: q = Polyomino([(0,0,0), (0,0,1)])
        sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
        sage: T = TilingSolver([p,q,r], box=(1,1,6))
        sage: T
        Tiling solver of 3 pieces into the box (1, 1, 6)
        Rotation allowed: True
        Reflection allowed: False
        Reusing pieces allowed: False

    Solutions are given by an iterator::

        sage: it = T.solve()
        sage: for p in next(it): p
        Polyomino: [(0, 0, 0)], Color: gray
        Polyomino: [(0, 0, 1), (0, 0, 2)], Color: gray
        Polyomino: [(0, 0, 3), (0, 0, 4), (0, 0, 5)], Color: gray

    Another solution::

        sage: for p in next(it): p
        Polyomino: [(0, 0, 0)], Color: gray
        Polyomino: [(0, 0, 1), (0, 0, 2), (0, 0, 3)], Color: gray
        Polyomino: [(0, 0, 4), (0, 0, 5)], Color: gray

    TESTS::

        sage: T = TilingSolver([p,q,r], box=(1,1,6), rotation=False, reflection=True)
        Traceback (most recent call last):
        ...
        NotImplementedError: When reflection is allowed and rotation is not allowed
    """
    def __init__(self, pieces, box, rotation=True,
                 reflection=False, reusable=False):
        r"""
        Constructor.

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: T
            Tiling solver of 3 pieces into the box (1, 1, 6)
            Rotation allowed: True
            Reflection allowed: False
            Reusing pieces allowed: False
        """
        self._pieces = pieces
        self._box = box
        self._rotation = rotation
        self._reflection = reflection
        if not self._rotation and self._reflection:
            raise NotImplementedError("When reflection is allowed and "
                                      "rotation is not allowed")
        self._reusable = reusable
        self._starting_rows = None    # the starting row of each piece

    def __repr__(self):
        r"""
        String representation

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: TilingSolver([p,q,r], box=(1,1,6))
            Tiling solver of 3 pieces into the box (1, 1, 6)
            Rotation allowed: True
            Reflection allowed: False
            Reusing pieces allowed: False

        """
        N = len(self._pieces)
        s = "Tiling solver of %s pieces into the box %s\n" % (N, self._box)
        s += "Rotation allowed: %s\n" % self._rotation
        s += "Reflection allowed: %s\n" % self._reflection
        s += "Reusing pieces allowed: %s" % self._reusable
        return s

    def is_suitable(self):
        r"""
        Return whether the volume of the box is equal to sum of the volume
        of the polyominoes and the number of rows sent to the DLX solver is
        larger than zero.

        If these conditions are not verified, then the problem is not suitable
        in the sense that there are no solution.

        .. NOTE::

            The DLX solver throws a Segmentation Fault when the
            number of rows is zero::

                sage: from sage.combinat.matrices.dancing_links import dlx_solver
                sage: rows = []
                sage: x = dlx_solver(rows)
                sage: x.search()        # not tested
                BOOM !!!

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: T.is_suitable()
            True
            sage: T = TilingSolver([p,q,r], box=(1,1,7))
            sage: T.is_suitable()
            False
        """
        if self._reusable:
            return len(self.rows()) != 0
        else:
            return (sum(len(p) for p in self.pieces()) == prod(self._box)
                    and len(self.rows()) != 0)

    def pieces(self):
        r"""
        Return the list of pieces.

        OUTPUT:

            list of 3d polyominoes

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: for p in T._pieces: p
            Polyomino: [(0, 0, 0)], Color: gray
            Polyomino: [(0, 0, 0), (0, 0, 1)], Color: gray
            Polyomino: [(0, 0, 0), (0, 0, 1), (0, 0, 2)], Color: gray
        """
        return self._pieces

    def space(self):
        r"""
        Returns an iterator over all the non negative integer coordinates
        contained in the box.

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: list(T.space())
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 0, 4), (0, 0, 5)]
        """
        return xmrange(self._box, tuple)

    @cached_method
    def coord_to_int_dict(self):
        r"""
        Returns a dictionary mapping coordinates to integers.

        OUTPUT:

            dict

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: A = T.coord_to_int_dict()
            sage: sorted(A.iteritems())
            [((0, 0, 0), 3), ((0, 0, 1), 4), ((0, 0, 2), 5), ((0, 0, 3), 6), ((0, 0, 4), 7), ((0, 0, 5), 8)]

        Reusable pieces::

            sage: p = Polyomino([(0,0), (0,1)])
            sage: q = Polyomino([(0,0), (0,1), (1,0), (1,1)])
            sage: T = TilingSolver([p,q], box=[3,2], reusable=True)
            sage: B = T.coord_to_int_dict()
            sage: sorted(B.iteritems())
            [((0, 0), 0), ((0, 1), 1), ((1, 0), 2), ((1, 1), 3), ((2, 0), 4), ((2, 1), 5)]
        """
        if self._reusable:
            return dict((c, i) for i, c in enumerate(self.space()))
        else:
            number_of_pieces = len(self._pieces)
            return dict((c, i+number_of_pieces)
                        for i, c in enumerate(self.space()))

    @cached_method
    def int_to_coord_dict(self):
        r"""
        Returns a dictionary mapping integers to coordinates.

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: B = T.int_to_coord_dict()
            sage: sorted(B.iteritems())
            [(3, (0, 0, 0)), (4, (0, 0, 1)), (5, (0, 0, 2)), (6, (0, 0, 3)), (7, (0, 0, 4)), (8, (0, 0, 5))]

        Reusable pieces::

            sage: from sage.combinat.tiling import Polyomino, TilingSolver
            sage: p = Polyomino([(0,0), (0,1)])
            sage: q = Polyomino([(0,0), (0,1), (1,0), (1,1)])
            sage: T = TilingSolver([p,q], box=[3,2], reusable=True)
            sage: B = T.int_to_coord_dict()
            sage: sorted(B.iteritems())
            [(0, (0, 0)), (1, (0, 1)), (2, (1, 0)), (3, (1, 1)), (4, (2, 0)), (5, (2, 1))]

        TESTS:

        The methods ``int_to_coord_dict`` and ``coord_to_int_dict`` returns
        dictionary that are inverse of each other::

            sage: A = T.coord_to_int_dict()
            sage: B = T.int_to_coord_dict()
            sage: all(A[B[i]] == i for i in B)
            True
            sage: all(B[A[i]] == i for i in A)
            True

        """
        if self._reusable:
            return dict((i, c) for i, c in enumerate(self.space()))
        else:
            number_of_pieces = len(self._pieces)
            return dict((i+number_of_pieces, c)
                        for i, c in enumerate(self.space()))

    @cached_method
    def rows(self):
        r"""
        Creation of the rows

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: rows = T.rows()
            sage: for row in rows: row
            [0, 3]
            [0, 4]
            [0, 5]
            [0, 6]
            [0, 7]
            [0, 8]
            [1, 3, 4]
            [1, 4, 5]
            [1, 5, 6]
            [1, 6, 7]
            [1, 8, 7]
            [2, 3, 4, 5]
            [2, 4, 5, 6]
            [2, 5, 6, 7]
            [2, 8, 6, 7]
        """
        coord_to_int = self.coord_to_int_dict()
        rows = []
        self._starting_rows = []  # indices of the first row for each piece
        for i, p in enumerate(self._pieces):
            self._starting_rows.append(len(rows))
            if self._rotation and self._reflection:
                it = p.translated_orthogonals(self._box,
                                              orientation_preserving=False)
            elif self._rotation and not self._reflection:
                it = p.translated_orthogonals(self._box,
                                              orientation_preserving=True)
            elif not self._rotation and self._reflection:
                raise NotImplementedError("Reflection allowed, Rotation not "
                                          "allowed is not implemented")
            else:
                it = p.translated(self._box)
            if self._reusable:
                for q in it:
                    rows.append([coord_to_int[coord] for coord in q])
            else:
                for q in it:
                    rows.append([i] + [coord_to_int[coord] for coord in q])
        self._starting_rows.append(len(rows))
        return rows

    def nrows_per_piece(self):
        r"""
        Return the number of rows necessary by each piece.

        OUPUT:

            list

        EXAMPLES::

            sage: from sage.games.quantumino import QuantuminoSolver
            sage: q = QuantuminoSolver(0)
            sage: T = q.tiling_solver()
            sage: T.nrows_per_piece()             # long time (10s)
            [360, 360, 360, 360, 360, 180, 180, 672, 672, 360, 360, 180, 180, 360, 360, 180]
        """
        if self._starting_rows is None:
            self.rows()
        L = self._starting_rows
        return [L[i+1] - L[i] for i in xrange(len(L)-1)]

    def dlx_solver(self):
        r"""
        Return the sage DLX solver of that 3d tiling problem.

        OUTPUT:

            DLX Solver

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: x = T.dlx_solver()
            sage: x
            <sage.combinat.matrices.dancing_links.dancing_linksWrapper object at ...>
        """
        from sage.combinat.matrices.dancing_links import dlx_solver
        rows = self.rows()
        assert len(rows) != 0, "Number of rows given to the DLX solver must not be zero"
        x = dlx_solver(rows)
        return x

    def dlx_solutions(self):
        r"""
        Return an iterator over the row indices of the solutions.

        OUPUT:

            iterator

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: list(T.dlx_solutions())
            [[0, 7, 14], [0, 12, 10], [6, 13, 5], [6, 14, 2], [11, 9, 5], [11, 10, 3]]
        """
        if len(self.rows()) == 0:
            raise StopIteration
        x = self.dlx_solver()
        while x.search() == 1:
            yield x.get_solution()

    def dlx_common_prefix_solutions(self):
        r"""
        Return an iterator over the row indices of solutions and of partial
        solutions, i.e. the common prefix of two consecutive solutions.

        The purpose is to illustrate the backtracking and construct an
        animation of the evolution of solutions.

        OUPUT:

            iterator

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: list(T.dlx_solutions())
            [[0, 7, 14], [0, 12, 10], [6, 13, 5], [6, 14, 2], [11, 9, 5], [11, 10, 3]]
            sage: list(T.dlx_common_prefix_solutions())
            [[0, 7, 14], [0], [0, 12, 10], [], [6, 13, 5], [6], [6, 14, 2], [], [11, 9, 5], [11], [11, 10, 3]]

        ::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: y = Polyomino([(0,0),(1,0),(2,0),(3,0),(2,1)], color='yellow')
            sage: T = TilingSolver([y], box=(5,10), reusable=True, reflection=True)
            sage: for a in T.dlx_common_prefix_solutions(): a
            [64, 83, 149, 44, 179, 62, 35, 162, 132, 101]
            [64, 83, 149, 44, 179]
            [64, 83, 149, 44, 179, 154, 35, 162, 132, 175]
            [64, 83, 149]
            [64, 83, 149, 97, 39, 162, 35, 62, 48, 106]
            [64]
            [64, 157, 149, 136, 179, 62, 35, 162, 132, 101]
            [64, 157, 149, 136, 179]
            [64, 157, 149, 136, 179, 154, 35, 162, 132, 175]
            []
            [82, 119, 58, 97, 38, 87, 8, 63, 48, 107]
            [82, 119, 58, 97, 38]
            [82, 119, 58, 97, 38, 161, 8, 63, 140, 107]
            [82, 119]
            [82, 119, 150, 136, 180, 63, 8, 161, 131, 175]
            [82, 119, 150]
            [82, 119, 150, 171, 38, 87, 8, 63, 48, 107]
            [82, 119, 150, 171, 38]
            [82, 119, 150, 171, 38, 161, 8, 63, 140, 107]
        """
        it = self.dlx_solutions()
        B = next(it)
        while True:
            yield B
            A, B = B, next(it)
            common_prefix = []
            for a, b in itertools.izip(A, B):
                if a == b:
                    common_prefix.append(a)
                else:
                    break
            yield common_prefix

    def dlx_incremental_solutions(self):
        r"""
        Return an iterator over the row indices of the incremental
        solutions.

        Between two incremental solution, either one piece is added or one
        piece is removed.

        The purpose is to illustrate the backtracking and construct an
        animation of the evolution of solutions.

        OUPUT:

            iterator

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: list(T.dlx_solutions())
            [[0, 7, 14], [0, 12, 10], [6, 13, 5], [6, 14, 2], [11, 9, 5], [11, 10, 3]]
            sage: list(T.dlx_incremental_solutions())
            [[0, 7, 14], [0, 7], [0], [0, 12], [0, 12, 10], [0, 12], [0], [], [6], [6, 13], [6, 13, 5], [6, 13], [6], [6, 14], [6, 14, 2], [6, 14], [6], [], [11], [11, 9], [11, 9, 5], [11, 9], [11], [11, 10], [11, 10, 3]]

        ::

            sage: y = Polyomino([(0,0),(1,0),(2,0),(3,0),(2,1)], color='yellow')
            sage: T = TilingSolver([y], box=(5,10), reusable=True, reflection=True)
            sage: for a in T.dlx_solutions(): a
            [64, 83, 149, 44, 179, 62, 35, 162, 132, 101]
            [64, 83, 149, 44, 179, 154, 35, 162, 132, 175]
            [64, 83, 149, 97, 39, 162, 35, 62, 48, 106]
            [64, 157, 149, 136, 179, 62, 35, 162, 132, 101]
            [64, 157, 149, 136, 179, 154, 35, 162, 132, 175]
            [82, 119, 58, 97, 38, 87, 8, 63, 48, 107]
            [82, 119, 58, 97, 38, 161, 8, 63, 140, 107]
            [82, 119, 150, 136, 180, 63, 8, 161, 131, 175]
            [82, 119, 150, 171, 38, 87, 8, 63, 48, 107]
            [82, 119, 150, 171, 38, 161, 8, 63, 140, 107]
            sage: len(list(T.dlx_incremental_solutions()))
            123
        """
        it = self.dlx_solutions()
        B = next(it)
        while True:
            yield B
            A, B = B, next(it)
            common_prefix = 0
            for a, b in itertools.izip(A, B):
                if a == b:
                    common_prefix += 1
                else:
                    break
            for i in xrange(1, len(A)-common_prefix):
                yield A[:-i]
            for j in xrange(common_prefix, len(B)):
                yield B[:j]

    def solve(self, partial=None):
        r"""
        Returns an iterator of list of 3d polyominoes that are an exact
        cover of the box.

        INPUT:

        - ``partial`` - string (optional, default: ``None``), whether to
          include partial (incomplete) solutions. It can be one of the
          following:

          - ``None`` - include only complete solution
          - ``'common_prefix'`` - common prefix between two consecutive solutions
          - ``'incremental'`` - one piece change at a time

        OUTPUT:

            iterator of list of 3d polyominoes

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: it = T.solve()
            sage: for p in next(it): p
            Polyomino: [(0, 0, 0)], Color: gray
            Polyomino: [(0, 0, 1), (0, 0, 2)], Color: gray
            Polyomino: [(0, 0, 3), (0, 0, 4), (0, 0, 5)], Color: gray
            sage: for p in next(it): p
            Polyomino: [(0, 0, 0)], Color: gray
            Polyomino: [(0, 0, 1), (0, 0, 2), (0, 0, 3)], Color: gray
            Polyomino: [(0, 0, 4), (0, 0, 5)], Color: gray
            sage: for p in next(it): p
            Polyomino: [(0, 0, 0), (0, 0, 1)], Color: gray
            Polyomino: [(0, 0, 2), (0, 0, 3), (0, 0, 4)], Color: gray
            Polyomino: [(0, 0, 5)], Color: gray

        Including the partial solutions::

            sage: it = T.solve(partial='common_prefix')
            sage: for p in next(it): p
            Polyomino: [(0, 0, 0)], Color: gray
            Polyomino: [(0, 0, 1), (0, 0, 2)], Color: gray
            Polyomino: [(0, 0, 3), (0, 0, 4), (0, 0, 5)], Color: gray
            sage: for p in next(it): p
            Polyomino: [(0, 0, 0)], Color: gray
            sage: for p in next(it): p
            Polyomino: [(0, 0, 0)], Color: gray
            Polyomino: [(0, 0, 1), (0, 0, 2), (0, 0, 3)], Color: gray
            Polyomino: [(0, 0, 4), (0, 0, 5)], Color: gray
            sage: for p in next(it): p
            sage: for p in next(it): p
            Polyomino: [(0, 0, 0), (0, 0, 1)], Color: gray
            Polyomino: [(0, 0, 2), (0, 0, 3), (0, 0, 4)], Color: gray
            Polyomino: [(0, 0, 5)], Color: gray

        Colors are preserved when the polyomino can be reused::

            sage: p = Polyomino([(0,0,0)], color='yellow')
            sage: q = Polyomino([(0,0,0), (0,0,1)], color='yellow')
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)], color='yellow')
            sage: T = TilingSolver([p,q,r], box=(1,1,6), reusable=True)
            sage: it = T.solve()
            sage: for p in next(it): p
            Polyomino: [(0, 0, 0)], Color: yellow
            Polyomino: [(0, 0, 1)], Color: yellow
            Polyomino: [(0, 0, 2)], Color: yellow
            Polyomino: [(0, 0, 3)], Color: yellow
            Polyomino: [(0, 0, 4)], Color: yellow
            Polyomino: [(0, 0, 5)], Color: yellow

        TESTS::

            sage: T = TilingSolver([p,q,r], box=(1,1,7))
            sage: next(T.solve())
            Traceback (most recent call last):
            ...
            StopIteration

        """
        if not self.is_suitable():
            raise StopIteration
        int_to_coord = self.int_to_coord_dict()
        rows = self.rows()
        if partial is None:
            it = self.dlx_solutions()
        elif partial == 'common_prefix':
            it = self.dlx_common_prefix_solutions()
        elif partial == 'incremental':
            it = self.dlx_incremental_solutions()
        else:
            raise ValueError("Unknown value for partial (=%s)" % partial)
        for solution in it:
            pieces = []
            for row_number in solution:
                row = rows[row_number]
                if self._reusable:
                    no = -1
                    while self._starting_rows[no] < row_number:
                        no += 1
                    coords = [int_to_coord[i] for i in row]
                    p = Polyomino(coords, color=self._pieces[no].color())
                else:
                    no = row[0]
                    coords = [int_to_coord[i] for i in row[1:]]
                    p = Polyomino(coords, color=self._pieces[no].color())
                pieces.append(p)
            yield pieces

    def number_of_solutions(self):
        r"""
        Return the number of distinct solutions.

        OUPUT:

            integer

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0)])
            sage: q = Polyomino([(0,0), (0,1)])
            sage: r = Polyomino([(0,0), (0,1), (0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,6))
            sage: T.number_of_solutions()
            6

        ::

            sage: T = TilingSolver([p,q,r], box=(1,7))
            sage: T.number_of_solutions()
            0
        """
        if not self.is_suitable():
            return 0
        x = self.dlx_solver()
        N = 0
        while x.search() == 1:
            N += 1
        return N

    def animate(self, partial=None, stop=None, size=0.75, axes=False):
        r"""
        Return an animation of evolving solutions.

        INPUT:

        - ``partial`` - string (optional, default: ``None``), whether to
          include partial (incomplete) solutions. It can be one of the
          following:

          - ``None`` - include only complete solutions
          - ``'common_prefix'`` - common prefix between two consecutive solutions
          - ``'incremental'`` - one piece change at a time

        - ``stop`` - integer (optional, default:``None``), number of frames

        - ``size`` - number (optional, default: ``0.75``), the size of each
          ``1 \times 1`` square. This does a homothety with respect
          to the center of each polyomino.

        - ``axes`` - bool (optional, default:``False``), whether the x and
          y axes are shown.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino, TilingSolver
            sage: y = Polyomino([(0,0),(1,0),(2,0),(3,0),(2,1)], color='cyan')
            sage: T = TilingSolver([y], box=(5,10), reusable=True, reflection=True)
            sage: a = T.animate()
            sage: a             # optional -- ImageMagick
            Animation with 10 frames

        Include partial solutions (common prefix between two consecutive
        solutions)::

            sage: a = T.animate('common_prefix')
            sage: a             # optional -- ImageMagick
            Animation with 19 frames

        Incremental solutions (one piece removed or added at a time)::

            sage: a = T.animate('incremental')      # long time (2s)
            sage: a                                 # long time (2s)  # optional -- ImageMagick
            Animation with 123 frames

        ::

            sage: a.show()                 # optional -- ImageMagick

        The ``show`` function takes arguments to specify the delay between
        frames (measured in hundredths of a second, default value 20) and
        the number of iterations (default value 0, which means to iterate
        forever). To iterate 4 times with half a second between each frame::

            sage: a.show(delay=50, iterations=4) # optional

        Limit the number of frames::

            sage: a = T.animate('incremental', stop=13)     # not tested
            sage: a                                         # not tested
            Animation with 13 frames
        """
        dimension = len(self._box)
        if dimension == 2:
            it = self.solve(partial=partial)
            it = itertools.islice(it, stop)
            L = [sum([piece.show2d(size)
                      for piece in solution], Graphics()) for solution in it]
            xmax, ymax = self._box
            xmax = xmax-0.5
            ymax = ymax-0.5
            a = Animation(L, xmin=-0.5, ymin=-0.5,
                          xmax=xmax, ymax=ymax, aspect_ratio=1, axes=axes)
            return a
        elif dimension == 3:
            raise NotImplementedError("3d Animation must be implemented "
                                      "in Jmol first")
        else:
            raise NotImplementedError("Dimension must be 2 or 3 in order "
                                      "to make an animation")
