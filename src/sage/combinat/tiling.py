# sage.doctest: needs sage.combinat sage.modules
r"""
Tiling Solver

Tiling a n-dimensional polyomino with n-dimensional polyominoes.

This module defines two classes:

- :class:`sage.combinat.tiling.Polyomino` class, to represent polyominoes
  in arbitrary dimension. The goal of this class is to return all the
  rotated, reflected and/or translated copies of a polyomino that are
  contained in a certain box.

- :class:`sage.combinat.tiling.TilingSolver` class, to solve the problem of
  tiling a `n`-dimensional polyomino with a set of `n`-dimensional
  polyominoes. One can specify if rotations and reflections are allowed or
  not and if pieces can be reused or not. This class convert the tiling
  data into rows of a matrix that are passed to the DLX solver. It also
  allows to compute the number of solutions.

This uses dancing links code which is in Sage. Dancing links were
originally introduced by Donald Knuth in 2000 [Knuth1]_. Knuth used dancing
links to solve tilings of a region by 2d pentaminoes. Here we extend the
method to any dimension.

In particular, the :mod:`sage.games.quantumino` module is based on
the Tiling Solver and allows to solve the 3d Quantumino puzzle.

AUTHOR:

- Sébastien Labbé, June 2011, initial version
- Sébastien Labbé, July 2015, count solutions up to rotations
- Sébastien Labbé, April 2017, tiling a polyomino, not only a rectangular box

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
    [Polyomino: [(0, 0), (0, 1), (1, 0), (1, 1)], Color: gray,
     Polyomino: [(2, 0), (2, 1)], Color: gray]
    sage: next(it)
    [Polyomino: [(1, 0), (1, 1), (2, 0), (2, 1)], Color: gray,
     Polyomino: [(0, 0), (0, 1)], Color: gray]
    sage: next(it)
    Traceback (most recent call last):
    ...
    StopIteration
    sage: T.number_of_solutions()
    2

Scott's pentamino problem
-------------------------

As mentioned in the introduction of [Knuth1]_, Scott's pentamino problem
consists in tiling a chessboard leaving the center four squares vacant with
the 12 distinct pentaminoes.

The 12 pentaminoes::

    sage: from sage.combinat.tiling import Polyomino
    sage: I = Polyomino([(0,0),(1,0),(2,0),(3,0),(4,0)], color='brown')
    sage: N = Polyomino([(1,0),(1,1),(1,2),(0,2),(0,3)], color='yellow')
    sage: L = Polyomino([(0,0),(1,0),(0,1),(0,2),(0,3)], color='magenta')
    sage: U = Polyomino([(0,0),(1,0),(0,1),(0,2),(1,2)], color='violet')
    sage: X = Polyomino([(1,0),(0,1),(1,1),(1,2),(2,1)], color='pink')
    sage: W = Polyomino([(2,0),(2,1),(1,1),(1,2),(0,2)], color='green')
    sage: P = Polyomino([(1,0),(2,0),(0,1),(1,1),(2,1)], color='orange')
    sage: F = Polyomino([(1,0),(1,1),(0,1),(2,1),(2,2)], color='gray')
    sage: Z = Polyomino([(0,0),(1,0),(1,1),(1,2),(2,2)], color='yellow')
    sage: T = Polyomino([(0,0),(0,1),(1,1),(2,1),(0,2)], color='red')
    sage: Y = Polyomino([(0,0),(1,0),(2,0),(3,0),(2,1)], color='green')
    sage: V = Polyomino([(0,0),(0,1),(0,2),(1,0),(2,0)], color='blue')

A `8 \times 8` chessboard leaving the center four squares vacant::

    sage: import itertools
    sage: s = set(itertools.product(range(8), repeat=2))
    sage: s.difference_update([(3,3), (3,4), (4,3), (4,4)])
    sage: chessboard = Polyomino(s)
    sage: len(chessboard)
    60

This problem is represented by a matrix made of 1568 rows and 72 columns.
It has 65 different solutions up to isometries::

    sage: from sage.combinat.tiling import TilingSolver
    sage: T = TilingSolver([I,N,L,U,X,W,P,F,Z,T,Y,V], box=chessboard, reflection=True)
    sage: T
    Tiling solver of 12 pieces into a box of size 60
    Rotation allowed: True
    Reflection allowed: True
    Reusing pieces allowed: False
    sage: len(T.rows())                # long time
    1568
    sage: T.number_of_solutions()      # long time
    520
    sage: 520 / 8
    65

Showing one solution::

    sage: solution = next(T.solve())                                    # long time
    sage: G = sum([piece.show2d() for piece in solution], Graphics())   # long time, needs sage.plot
    sage: G.show(aspect_ratio=1, axes=False)                            # long time, needs sage.plot

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
    [Polyomino: [(0)], Color: gray,
     Polyomino: [(1), (2)], Color: gray,
     Polyomino: [(3), (4), (5)], Color: gray]
    sage: next(it)
    [Polyomino: [(0)], Color: gray,
     Polyomino: [(1), (2), (3)], Color: gray,
     Polyomino: [(4), (5)], Color: gray]
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
there are no solutions for tiling a `8 \times 8` rectangular box::

    sage: T = TilingSolver(L, box=(8,8))
    sage: T.number_of_solutions()                                       # long time (2.5s)
    0

If reflections are allowed, there are solutions. Solve the puzzle and show
one solution::

    sage: T = TilingSolver(L, box=(8,8), reflection=True)
    sage: solution = next(T.solve())                                    # long time (7s)
    sage: G = sum([piece.show2d() for piece in solution], Graphics())   # long time (<1s), needs sage.plot
    sage: G.show(aspect_ratio=1, axes=False)                            # long time (2s), needs sage.plot

Compute the number of solutions::

    sage: T.number_of_solutions()                                       # long time (2.6s)
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

    sage: T = TilingSolver(L, box=(8,8,1))
    sage: solution = next(T.solve())                                    # long time (8s)
    sage: G = sum([p.show3d(size=0.85) for p in solution], Graphics())  # long time (<1s), needs sage.plot
    sage: G.show(aspect_ratio=1, viewer='tachyon')                      # long time (2s), needs sage.plot

Let us compute the number of solutions::

    sage: T.number_of_solutions()                                       # long time (3s)
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
    sage: G = sum([p.show2d() for p in solution], Graphics())                           # needs sage.plot
    sage: G.show(aspect_ratio=1)                                        # long time (2s), needs sage.plot

::

    sage: T = TilingSolver([y], box=(15,15), reusable=True, reflection=True)
    sage: T.number_of_solutions()                      # not tested
    1696

Up to the symmetries of the square, there are 212 distinct solutions::

    sage: 1696 // 8
    212

Animation of Donald Knuth's dancing links
-----------------------------------------

Animation of the solutions::

    sage: from sage.combinat.tiling import Polyomino, TilingSolver
    sage: Y = Polyomino([(0,0),(1,0),(2,0),(3,0),(2,1)], color='yellow')
    sage: T = TilingSolver([Y], box=(15,15), reusable=True, reflection=True)
    sage: a = T.animate(stop=40); a                     # long time, optional - imagemagick, needs sage.plot
    Animation with 40 frames

Incremental animation of the solutions (one piece is removed/added at a time)::

    sage: a = T.animate('incremental', stop=40); a      # long time, optional - imagemagick, needs sage.plot
    Animation with 40 frames
    sage: a.show(delay=50, iterations=1)                # long time, optional - imagemagick, needs sage.plot

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
    [[0], [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12],
     [13], [14], [15], [16], [17], [18], [19], [20], [21], [22], [23], [24],
     [25], [26], [27], [28], [29], [30], [31]]
    sage: T.number_of_solutions()                       # long time (fast)
    1

REFERENCES:

.. [Knuth1] Knuth, Donald (2000). "Dancing links". :arxiv:`cs/0011047`.
"""
# ****************************************************************************
#       Copyright (C) 2011-2015 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools
from sage.structure.sage_object import SageObject
from sage.modules.free_module_element import vector
from sage.misc.cachefunc import cached_method, cached_function


#######################################
# n-cube isometry group transformations
#######################################


def ncube_isometry_group(n, orientation_preserving=True):
    r"""
    Return the isometry group of the `n`-cube as a list of matrices.

    INPUT:

    - ``n`` -- positive integer, dimension of the space
    - ``orientation_preserving`` -- boolean (default: ``True``);
      whether the orientation is preserved

    OUTPUT: list of matrices

    EXAMPLES::

        sage: from sage.combinat.tiling import ncube_isometry_group
        sage: ncube_isometry_group(2)
        [
        [1 0]  [ 0  1]  [-1  0]  [ 0 -1]
        [0 1], [-1  0], [ 0 -1], [ 1  0]
        ]
        sage: ncube_isometry_group(2, orientation_preserving=False)
        [
        [1 0]  [ 0 -1]  [ 1  0]  [ 0  1]  [0 1]  [-1  0]  [ 0 -1]  [-1  0]
        [0 1], [-1  0], [ 0 -1], [-1  0], [1 0], [ 0 -1], [ 1  0], [ 0  1]
        ]

    There are 24 orientation preserving isometries of the 3-cube::

        sage: ncube_isometry_group(3)
        [
        [1 0 0]  [ 1  0  0]  [ 1  0  0]  [ 0  1  0]  [0 1 0]  [ 0  0  1]
        [0 1 0]  [ 0  0  1]  [ 0  0 -1]  [-1  0  0]  [0 0 1]  [ 0 -1  0]
        [0 0 1], [ 0 -1  0], [ 0  1  0], [ 0  0  1], [1 0 0], [ 1  0  0],
        <BLANKLINE>
        [-1  0  0]  [ 0 -1  0]  [-1  0  0]  [-1  0  0]  [ 0 -1  0]  [ 0  0 -1]
        [ 0 -1  0]  [ 0  0 -1]  [ 0  0 -1]  [ 0  1  0]  [ 0  0  1]  [ 1  0  0]
        [ 0  0  1], [ 1  0  0], [ 0 -1  0], [ 0  0 -1], [-1  0  0], [ 0 -1  0],
        <BLANKLINE>
        [ 0  1  0]  [ 0  0  1]  [0 0 1]  [ 0 -1  0]  [ 0  0 -1]  [-1  0  0]
        [ 1  0  0]  [ 0  1  0]  [1 0 0]  [ 1  0  0]  [ 0  1  0]  [ 0  0  1]
        [ 0  0 -1], [-1  0  0], [0 1 0], [ 0  0  1], [ 1  0  0], [ 0  1  0],
        <BLANKLINE>
        [ 0 -1  0]  [ 0  0 -1]  [ 0  0  1]  [ 1  0  0]  [ 0  0 -1]  [ 0  1  0]
        [-1  0  0]  [-1  0  0]  [-1  0  0]  [ 0 -1  0]  [ 0 -1  0]  [ 0  0 -1]
        [ 0  0 -1], [ 0  1  0], [ 0 -1  0], [ 0  0 -1], [-1  0  0], [-1  0  0]
        ]

    TESTS::

        sage: ncube_isometry_group(1)
        [[1]]
        sage: ncube_isometry_group(0)
        Traceback (most recent call last):
        ...
        ValueError: ['B', 0] is not a valid Cartan type
    """
    from sage.combinat.root_system.weyl_group import WeylGroup
    L = [w.matrix() for w in WeylGroup(['B', n])]
    if orientation_preserving:
        return [m for m in L if m.det() == 1]
    else:
        return L


@cached_function
def ncube_isometry_group_cosets(n, orientation_preserving=True):
    r"""
    Return the quotient of the isometry group of the `n`-cube by the
    the isometry group of the rectangular parallelepiped.

    INPUT:

    - ``n`` -- positive integer, dimension of the space
    - ``orientation_preserving`` -- boolean (default: ``True``);
      whether the orientation is preserved

    OUTPUT: list of cosets, each coset being a sorted list of matrices

    EXAMPLES::

        sage: from sage.combinat.tiling import ncube_isometry_group_cosets
        sage: sorted(ncube_isometry_group_cosets(2))
        [[
        [-1  0]  [1 0]
        [ 0 -1], [0 1]
        ], [
        [ 0 -1]  [ 0  1]
        [ 1  0], [-1  0]
        ]]
        sage: sorted(ncube_isometry_group_cosets(2, False))
        [[
        [-1  0]  [-1  0]  [ 1  0]  [1 0]
        [ 0 -1], [ 0  1], [ 0 -1], [0 1]
        ], [
        [ 0 -1]  [ 0 -1]  [ 0  1]  [0 1]
        [-1  0], [ 1  0], [-1  0], [1 0]
        ]]

    ::

        sage: sorted(ncube_isometry_group_cosets(3))
        [[
        [-1  0  0]  [-1  0  0]  [ 1  0  0]  [1 0 0]
        [ 0 -1  0]  [ 0  1  0]  [ 0 -1  0]  [0 1 0]
        [ 0  0  1], [ 0  0 -1], [ 0  0 -1], [0 0 1]
        ], [
        [-1  0  0]  [-1  0  0]  [ 1  0  0]  [ 1  0  0]
        [ 0  0 -1]  [ 0  0  1]  [ 0  0 -1]  [ 0  0  1]
        [ 0 -1  0], [ 0  1  0], [ 0  1  0], [ 0 -1  0]
        ], [
        [ 0 -1  0]  [ 0 -1  0]  [ 0  1  0]  [ 0  1  0]
        [-1  0  0]  [ 1  0  0]  [-1  0  0]  [ 1  0  0]
        [ 0  0 -1], [ 0  0  1], [ 0  0  1], [ 0  0 -1]
        ], [
        [ 0 -1  0]  [ 0 -1  0]  [ 0  1  0]  [0 1 0]
        [ 0  0 -1]  [ 0  0  1]  [ 0  0 -1]  [0 0 1]
        [ 1  0  0], [-1  0  0], [-1  0  0], [1 0 0]
        ], [
        [ 0  0 -1]  [ 0  0 -1]  [ 0  0  1]  [0 0 1]
        [-1  0  0]  [ 1  0  0]  [-1  0  0]  [1 0 0]
        [ 0  1  0], [ 0 -1  0], [ 0 -1  0], [0 1 0]
        ], [
        [ 0  0 -1]  [ 0  0 -1]  [ 0  0  1]  [ 0  0  1]
        [ 0 -1  0]  [ 0  1  0]  [ 0 -1  0]  [ 0  1  0]
        [-1  0  0], [ 1  0  0], [ 1  0  0], [-1  0  0]
        ]]

    TESTS::

        sage: cosets = ncube_isometry_group_cosets(3, False)
        sage: len(cosets)
        6
        sage: [len(c) for c in cosets]
        [8, 8, 8, 8, 8, 8]
        sage: type(cosets[0][0])
        <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
    """
    from sage.misc.misc_c import prod
    from sage.matrix.constructor import diagonal_matrix
    G = ncube_isometry_group(n, orientation_preserving)

    # Construct the subgroup H of G of diagonal matrices
    it = itertools.product((1, -1), repeat=n)
    if orientation_preserving:
        H = [diagonal_matrix(L) for L in it if prod(L) == 1]
    else:
        H = [diagonal_matrix(L) for L in it]

    G_todo = set(G)
    # Make sure that H is a subset of G
    for h in H:
        h.set_immutable()
    assert all(h in G_todo for h in H), "H must be a subset of G"

    # Construct the cosets
    cosets = []
    for g in G:
        if g not in G_todo:
            continue
        left_coset = sorted(h*g for h in H)
        right_coset = sorted(g*h for h in H)
        assert left_coset == right_coset, "H must be a normal subgroup of G"
        for c in left_coset:
            c.set_immutable()
        G_todo.difference_update(left_coset)
        cosets.append(left_coset)
    return cosets

##############################
# Class Polyomino
##############################


class Polyomino(SageObject):
    r"""
    A polyomino in `\ZZ^d`.

    The polyomino is the union of the unit square (or cube, or n-cube)
    centered at those coordinates. Such an object should be connected, but
    the code does not make this assumption.

    INPUT:

    - ``coords`` -- iterable of integer coordinates in `\ZZ^d`
    - ``color`` -- string (default: ``'gray'``); color for display
    - ``dimension`` -- integer (default: ``None``); dimension of the space,
      if ``None``, it is guessed from the ``coords`` if ``coords`` is non
      empty

    EXAMPLES::

        sage: from sage.combinat.tiling import Polyomino
        sage: Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
        Polyomino: [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)], Color: blue
    """

    def __init__(self, coords, color='gray', dimension=None):
        r"""
        Constructor.

        See :mod:`Polyomino` for full documentation.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            Polyomino: [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)], Color: blue

        ::

            sage: Polyomino([(0,0), (1,0), (2,0)])
            Polyomino: [(0, 0), (1, 0), (2, 0)], Color: gray

        TESTS::

            sage: Polyomino([], dimension=2)
            Polyomino: [], Color: gray
        """
        from sage.modules.free_module import FreeModule
        from sage.rings.integer_ring import ZZ

        if not isinstance(color, str):
            raise TypeError("color = ({!r}) must be a string".format(color))
        self._color = color

        if not isinstance(coords, (tuple,list)):
            coords = list(coords)

        if dimension is None:
            if coords:
                self._dimension = ZZ(len(coords[0]))
            else:
                raise ValueError("dimension(={}) must be provided for"
                                 " the empty polyomino".format(dimension))
        else:
            self._dimension = dimension
        self._free_module = FreeModule(ZZ, self._dimension)

        self._blocs = coords
        self._blocs = [self._free_module(bloc) for bloc in self._blocs]
        for b in self._blocs:
            b.set_immutable()
        self._blocs = frozenset(self._blocs)

    def _repr_(self):
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

    def color(self, color=None):
        r"""
        Return or change the color of the polyomino.

        INPUT:

        - ``color`` -- string, RBG tuple or ``None`` (default: ``None``),
          if ``None``, it returns the current color

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: p.color()
            'blue'
        """
        if color is None:
            return self._color
        else:
            self._color = color

    def frozenset(self):
        r"""
        Return the elements of `\ZZ^d` in the polyomino as a frozenset.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='red')
            sage: p.frozenset()
            frozenset({(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)})
        """
        return self._blocs

    @cached_method
    def sorted_list(self):
        r"""
        Return the color of the polyomino.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: p.sorted_list()
            [(0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 1, 1)]
        """
        return sorted(self.frozenset())

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
        return len(self.frozenset())

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: it = iter(p)
            sage: next(it)
            (0, 0, 0)
        """
        return iter(self.sorted_list())

    @cached_method
    def bounding_box(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: p.bounding_box()
            [[0, 0, 0], [1, 2, 1]]
        """
        return [[min(w) for w in zip(*self)],
                [max(w) for w in zip(*self)]]

    def __hash__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: hash(p)     # random
            2059134902
        """
        return hash(self.frozenset())

    def __eq__(self, other):
        r"""
        Return whether ``self`` is equal to ``other``.

        INPUT:

        - ``other`` -- a polyomino

        OUTPUT: boolean

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
        return isinstance(other, Polyomino) and self.frozenset() == other.frozenset()

    def __ne__(self, other):
        r"""
        Return whether ``self`` is not equal to ``other``.

        INPUT:

        - ``other`` -- a polyomino

        OUTPUT: boolean

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
        return not self.frozenset() == other.frozenset()

    def __le__(self, other):
        r"""
        Return whether ``self`` is inside of ``other``.

        INPUT:

        - ``other`` -- a polyomino

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0)])
            sage: b = Polyomino([(0,0), (0,1), (1,1), (2,1)])
            sage: p <= b
            True
            sage: b <= p
            False
        """
        return isinstance(other, Polyomino) and self.frozenset() <= other.frozenset()

    def __ge__(self, other):
        r"""
        Return whether ``self`` contains ``other``.

        INPUT:

        - ``other`` -- a polyomino

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0)])
            sage: b = Polyomino([(0,0), (0,1), (1,1), (2,1)])
            sage: p >= b
            False
            sage: b >= p
            True
        """
        return isinstance(other, Polyomino) and self.frozenset() >= other.frozenset()

    def __lt__(self, other):
        r"""
        Return whether ``self`` is strictly inside of ``other``.

        INPUT:

        - ``other`` -- a polyomino

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0)])
            sage: b = Polyomino([(0,0), (0,1), (1,1), (2,1)])
            sage: p < b
            True
            sage: b < p
            False
        """
        return isinstance(other, Polyomino) and self.frozenset() < other.frozenset()

    def __gt__(self, other):
        r"""
        Return whether ``self`` strictly contains ``other``.

        INPUT:

        - ``other`` -- a polyomino

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0)])
            sage: b = Polyomino([(0,0), (0,1), (1,1), (2,1)])
            sage: p > b
            False
            sage: b > p
            True
        """
        return isinstance(other, Polyomino) and self.frozenset() > other.frozenset()

    def intersection(self, other):
        r"""
        Return the intersection of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a polyomino

        OUTPUT: polyomino

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: a = Polyomino([(0,0)])
            sage: b = Polyomino([(0,0), (0,1), (1,1), (2,1)])
            sage: a.intersection(b)
            Polyomino: [(0, 0)], Color: gray
            sage: a.intersection(b+(1,1))
            Polyomino: [], Color: gray
        """
        if not isinstance(other, Polyomino):
            raise TypeError("other(={}) must be a polyomino".format(other))
        return Polyomino(self.frozenset() & other.frozenset(),
                         color=self._color, dimension=self._dimension)

    def __sub__(self, v):
        r"""
        Return a translated copy of ``self`` by the opposite of the
        vector v.

        INPUT:

        - ``v`` -- tuple

        OUTPUT: polyomino

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: p - (2,2,2)
            Polyomino: [(-2, -2, -2), (-1, -2, -2), (-1, -1, -2), (-1, -1, -1), (-1, 0, -2)], Color: deeppink
        """
        v = self._free_module(v)
        return Polyomino([p-v for p in self], color=self._color)

    def __add__(self, v):
        r"""
        Return a translated copy of ``self`` by the vector v.

        INPUT:

        - ``v`` -- tuple

        OUTPUT: polyomino

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: p + (2,2,2)
            Polyomino: [(2, 2, 2), (3, 2, 2), (3, 3, 2), (3, 3, 3), (3, 4, 2)], Color: deeppink
        """
        v = self._free_module(v)
        return Polyomino([p+v for p in self], color=self._color)

    def __rmul__(self, m):
        r"""
        Return the image of the polyomino under the application of the
        matrix m.

        INPUT:

        - ``m`` -- square matrix, matching the dimension of ``self``

        OUTPUT: polyomino

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
        return Polyomino([m * p for p in self], color=self._color)

    def canonical(self):
        r"""
        Return the translated copy of ``self`` having minimal and nonnegative
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
        minxyz, _ = self.bounding_box()
        return self - minxyz

    def canonical_isometric_copies(self, orientation_preserving=True,
                                   mod_box_isometries=False):
        r"""
        Return the list of image of ``self`` under isometries of the `n`-cube
        where the coordinates are all nonnegative and minimal.

        INPUT:

        - ``orientation_preserving`` -- boolean (default: ``True``);
          if ``True``, the group of isometries of the `n`-cube is restricted
          to those that preserve the orientation, i.e. of determinant 1.

        - ``mod_box_isometries`` -- boolean (default: ``False``); whether to
          quotient the group of isometries of the `n`-cube by the
          subgroup of isometries of the `a_1\times a_2\cdots \times a_n`
          rectangular box where are the `a_i` are assumed to be distinct.

        OUTPUT: set of Polyomino

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: s = p.canonical_isometric_copies()
            sage: len(s)
            12

        With the non orientation-preserving::

            sage: s = p.canonical_isometric_copies(orientation_preserving=False)
            sage: len(s)
            24

        Modulo rotation by angle 180 degrees::

            sage: s = p.canonical_isometric_copies(mod_box_isometries=True)
            sage: len(s)
            3

        TESTS::

            sage: from sage.games.quantumino import pentaminos
            sage: [len(p.canonical_isometric_copies((5,8,2), mod_box_isometries=False)) for p in pentaminos]
            [24, 24, 24, 24, 24, 24, 12, 12, 24, 24, 24, 24, 12, 12, 24, 24, 12]
            sage: [len(p.canonical_isometric_copies((5,8,2), mod_box_isometries=True)) for p in pentaminos]
            [6, 6, 6, 6, 6, 6, 3, 3, 6, 6, 6, 6, 3, 3, 6, 6, 3]
        """
        if mod_box_isometries:
            L = ncube_isometry_group_cosets(self._dimension, orientation_preserving)
            P_cosets = {frozenset((m * self).canonical() for m in coset)
                        for coset in L}
            P_cosets_representents = [min(s, key=lambda a: a.sorted_list()) for s in P_cosets]
            return sorted(P_cosets_representents, key=lambda a: a.sorted_list())
        else:
            L = ncube_isometry_group(self._dimension, orientation_preserving)
            P_images = {(m * self).canonical() for m in L}
            return sorted(P_images, key=lambda a: a.sorted_list())

    def translated_copies(self, box):
        r"""
        Return an iterator over the translated images of ``self`` inside a
        polyomino.

        INPUT:

        - ``box`` -- polyomino or tuple of integers (size of a box)

        OUTPUT: iterator of 3d polyominoes

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: for t in p.translated_copies(box=(5,8,2)): t
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

        This method is independent of the translation of the polyomino::

            sage: q = Polyomino([(0,0,0), (1,0,0)])
            sage: list(q.translated_copies((2,2,1)))
            [Polyomino: [(0, 0, 0), (1, 0, 0)], Color: gray,
             Polyomino: [(0, 1, 0), (1, 1, 0)], Color: gray]
            sage: q = Polyomino([(34,7,-9), (35,7,-9)])
            sage: list(q.translated_copies((2,2,1)))
            [Polyomino: [(0, 0, 0), (1, 0, 0)], Color: gray,
             Polyomino: [(0, 1, 0), (1, 1, 0)], Color: gray]

        Inside smaller boxes::

            sage: list(p.translated_copies(box=(2,2,3)))
            []
            sage: list(p.translated_copies(box=(2,3,2)))
            [Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink]
            sage: list(p.translated_copies(box=(3,2,2)))
            []
            sage: list(p.translated_copies(box=(1,1,1)))
            []

        Using a Polyomino as input::

            sage: b = Polyomino([(0,0), (0,1), (0,2), (1,0), (1,1), (1,2)])
            sage: p = Polyomino([(0,0)])
            sage: list(p.translated_copies(b))
            [Polyomino: [(0, 0)], Color: gray,
             Polyomino: [(0, 1)], Color: gray,
             Polyomino: [(0, 2)], Color: gray,
             Polyomino: [(1, 0)], Color: gray,
             Polyomino: [(1, 1)], Color: gray,
             Polyomino: [(1, 2)], Color: gray]

        ::

            sage: p = Polyomino([(0,0), (1,0), (0,1)])
            sage: b = Polyomino([(0,0), (1,0), (2,0), (0,1), (1,1), (0,2)])
            sage: list(p.translated_copies(b))
            [Polyomino: [(0, 0), (0, 1), (1, 0)], Color: gray,
             Polyomino: [(0, 1), (0, 2), (1, 1)], Color: gray,
             Polyomino: [(1, 0), (1, 1), (2, 0)], Color: gray]
        """
        if not isinstance(box, Polyomino):
            ranges = [range(a) for a in box]
            box = Polyomino(itertools.product(*ranges))
        if not box._dimension == self._dimension:
            raise ValueError("Dimension of input box must match the "
                             "dimension of the polyomino")
        minxyz, maxxyz = self.bounding_box()
        minxyz, maxxyz = vector(minxyz), vector(maxxyz)
        size = maxxyz - minxyz
        boxminxyz, boxmaxxyz = box.bounding_box()
        ranges = [range(a, b-c+1) for (a,b,c) in zip(boxminxyz,
                                                   boxmaxxyz,
                                                   size)]
        cano = self.canonical()
        for v in itertools.product(*ranges):
            translated = cano + v
            if translated <= box:
                yield translated

    def translated_copies_intersection(self, box):
        r"""
        Return the set of non empty intersections of translated images of
        ``self`` with a polyomino.

        INPUT:

        - ``box`` -- polyomino or tuple of integers (size of a box)

        OUTPUT: set of 3d polyominoes

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0),(1,0)], color='deeppink')
            sage: sorted(sorted(a.frozenset())
            ....:        for a in p.translated_copies_intersection(box=(2,3)))
            [[(0, 0)],
             [(0, 0), (1, 0)],
             [(0, 1)],
             [(0, 1), (1, 1)],
             [(0, 2)],
             [(0, 2), (1, 2)],
             [(1, 0)],
             [(1, 1)],
             [(1, 2)]]

        Using a Polyomino as input::

            sage: b = Polyomino([(0,0), (0,1), (0,2), (1,0), (2,0)])
            sage: p = Polyomino([(0,0), (1,0)])
            sage: sorted(sorted(a.frozenset())
            ....:        for a in p.translated_copies_intersection(b))
            [[(0, 0)], [(0, 0), (1, 0)], [(0, 1)], [(0, 2)], [(1, 0), (2, 0)], [(2, 0)]]
        """
        if not isinstance(box, Polyomino):
            ranges = [range(a) for a in box]
            box = Polyomino(itertools.product(*ranges))
        if not box._dimension == self._dimension:
            raise ValueError("Dimension of input box must match the "
                             "dimension of the polyomino")
        minxyz, maxxyz = self.bounding_box()
        minxyz, maxxyz = vector(minxyz), vector(maxxyz)
        size = maxxyz - minxyz
        boxminxyz, boxmaxxyz = box.bounding_box()
        ranges = [range(a-c, b+1) for (a,b,c) in zip(boxminxyz,
                                                   boxmaxxyz,
                                                   size)]
        S = set()
        cano = self.canonical()
        for v in itertools.product(*ranges):
            translated = cano + v
            intersected = translated.intersection(box)
            if intersected:
                S.add(intersected)
        return S

    def isometric_copies(self, box, orientation_preserving=True,
            mod_box_isometries=False):
        r"""
        Return the translated and isometric images of ``self`` that lies in the box.

        INPUT:

        - ``box`` -- polyomino or tuple of integers (size of a box)

        - ``orientation_preserving`` -- boolean (default: ``True``);
          if ``True``, the group of isometries of the `n`-cube is restricted
          to those that preserve the orientation, i.e. of determinant 1.

        - ``mod_box_isometries`` -- boolean (default: ``False``); whether to
          quotient the group of isometries of the `n`-cube by the
          subgroup of isometries of the `a_1\times a_2\cdots \times a_n`
          rectangular box where are the `a_i` are assumed to be distinct.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,1,1),(1,2,0)], color='deeppink')
            sage: L = list(p.isometric_copies(box=(5,8,2)))
            sage: len(L)
            360

        ::

            sage: p = Polyomino([(0,0,0),(1,0,0),(1,1,0),(1,2,0),(1,2,1)], color='orange')
            sage: L = list(p.isometric_copies(box=(5,8,2)))
            sage: len(L)
            180
            sage: L = list(p.isometric_copies((5,8,2), False))
            sage: len(L)
            360
            sage: L = list(p.isometric_copies((5,8,2), mod_box_isometries=True))
            sage: len(L)
            45

        ::

            sage: p = Polyomino([(0,0), (1,0), (0,1)])
            sage: b = Polyomino([(0,0), (1,0), (2,0), (0,1), (1,1), (0,2)])
            sage: sorted(p.isometric_copies(b), key=lambda p: p.sorted_list())
            [Polyomino: [(0, 0), (0, 1), (1, 0)], Color: gray,
             Polyomino: [(0, 0), (0, 1), (1, 1)], Color: gray,
             Polyomino: [(0, 0), (1, 0), (1, 1)], Color: gray,
             Polyomino: [(0, 1), (0, 2), (1, 1)], Color: gray,
             Polyomino: [(0, 1), (1, 0), (1, 1)], Color: gray,
             Polyomino: [(1, 0), (1, 1), (2, 0)], Color: gray]
        """
        if not isinstance(box, Polyomino):
            ranges = [range(a) for a in box]
            box = Polyomino(itertools.product(*ranges))
        if not box._dimension == self._dimension:
            raise ValueError("Dimension of input box must match the "
                             "dimension of the polyomino")
        box_min_coords, box_max_coords = box.bounding_box()
        if mod_box_isometries and len({b - a for a, b in zip(box_min_coords,
                                                             box_max_coords)}) < box._dimension:
            raise NotImplementedError("The code below assumes that the"
                    " sizes of the box (={}) are all distinct when"
                    " argument `mod_box_isometries` is True.".format(box))
        all_distinct_cano = self.canonical_isometric_copies(orientation_preserving,
                                                            mod_box_isometries)
        for cano in all_distinct_cano:
            yield from cano.translated_copies(box=box)

    def isometric_copies_intersection(self, box, orientation_preserving=True):
        r"""
        Return the set of non empty intersections of isometric images of
        ``self`` with a polyomino.

        INPUT:

        - ``box`` -- polyomino or tuple of integers (size of a box)

        - ``orientation_preserving`` -- boolean (default: ``True``);
          if ``True``, the group of isometries of the `n`-cube is restricted
          to those that preserve the orientation, i.e. of determinant 1.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0),(1,0)], color='deeppink')
            sage: sorted(sorted(a.frozenset()) for a in p.isometric_copies_intersection(box=(2,3)))
            [[(0, 0)],
             [(0, 0), (0, 1)],
             [(0, 0), (1, 0)],
             [(0, 1)],
             [(0, 1), (0, 2)],
             [(0, 1), (1, 1)],
             [(0, 2)],
             [(0, 2), (1, 2)],
             [(1, 0)],
             [(1, 0), (1, 1)],
             [(1, 1)],
             [(1, 1), (1, 2)],
             [(1, 2)]]
        """
        all_distinct_cano = self.canonical_isometric_copies(orientation_preserving,
                                                            mod_box_isometries=False)
        return {t for cano in all_distinct_cano
                for t in cano.translated_copies_intersection(box=box)}

    def neighbor_edges(self):
        r"""
        Return an iterator over the pairs of neighbor coordinates inside of
        the polyomino.

        Two points `P` and `Q` in the polyomino are neighbor if `P - Q` has
        one coordinate equal to `+1` or `-1` and zero everywhere else.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0),(0,0,1)])
            sage: [sorted(edge) for edge in p.neighbor_edges()]
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
        return sum(self) / len(self)

    def boundary(self):
        r"""
        Return the boundary of a 2d polyomino.

        INPUT:

        - ``self`` -- a 2d polyomino

        OUTPUT:

        - list of edges (an edge is a pair of adjacent 2d coordinates)

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0), (1,0), (0,1), (1,1)])
            sage: sorted(p.boundary())
            [((-0.5, -0.5), (-0.5, 0.5)), ((-0.5, -0.5), (0.5, -0.5)),
             ((-0.5, 0.5), (-0.5, 1.5)), ((-0.5, 1.5), (0.5, 1.5)),
             ((0.5, -0.5), (1.5, -0.5)), ((0.5, 1.5), (1.5, 1.5)),
             ((1.5, -0.5), (1.5, 0.5)), ((1.5, 0.5), (1.5, 1.5))]
            sage: len(_)
            8
            sage: p = Polyomino([(5,5)])
            sage: sorted(p.boundary())
            [((4.5, 4.5), (4.5, 5.5)), ((4.5, 4.5), (5.5, 4.5)),
             ((4.5, 5.5), (5.5, 5.5)), ((5.5, 4.5), (5.5, 5.5))]
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
        for (x, y), coeff in horizontal.items():
            if coeff:
                edges.append(((x-h, y-h), (x+h, y-h)))
        for (x, y), coeff in vertical.items():
            if coeff:
                edges.append(((x-h, y-h), (x-h, y+h)))
        return edges

    def show3d(self, size=1):
        r"""
        Return a 3d Graphic object representing the polyomino.

        INPUT:

        - ``self`` -- a polyomino of dimension 3
        - ``size`` -- number (default: ``1``); the size of each
          ``1 \times 1 \times 1`` cube. This does a homothety with respect
          to the center of the polyomino.

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1)], color='blue')
            sage: p.show3d()                    # long time (2s)                        # needs sage.plot
            Graphics3d Object
        """
        assert self._dimension == 3, "Dimension of the polyomino must be 3."
        from sage.plot.graphics import Graphics
        from sage.plot.plot3d.platonic import cube
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
        Return a 2d Graphic object representing the polyomino.

        INPUT:

        - ``self`` -- a polyomino of dimension 2
        - ``size`` -- number (default: ``0.7``); the size of each square
        - ``color`` -- color (default: ``'black'``); color of the boundary line
        - ``thickness`` -- number (default: ``1``); how thick the boundary line is

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: p = Polyomino([(0,0),(1,0),(1,1),(1,2)], color='deeppink')
            sage: p.show2d()                    # long time (0.5s)                      # needs sage.plot
            Graphics object consisting of 17 graphics primitives
        """
        assert self._dimension == 2, "Dimension of the polyomino must be 2."
        from sage.plot.graphics import Graphics
        from sage.plot.circle import circle
        from sage.plot.line import line
        from sage.plot.polygon import polygon
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

    def self_surrounding(self, radius, remove_incomplete_copies=True,
                         ncpus=None):
        r"""
        Return a list of isometric copies of ``self`` surrounding it with an
        annulus of given radius.

        INPUT:

        - ``self`` -- a polyomino of dimension 2
        - ``radius`` -- integer
        - ``remove_incomplete_copies`` -- boolean (default: ``True``); whether
          to keep only complete copies of ``self`` in the output
        - ``ncpus`` -- integer (default: ``None``); maximal number of
          subprocesses to use at the same time. If ``None``, it detects the
          number of effective CPUs in the system using
          :func:`sage.parallel.ncpus.ncpus()`.
          If ``ncpus=1``, the first solution is searched serially.

        OUTPUT: list of polyominoes

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino
            sage: H = Polyomino([(-1, 1), (-1, 4), (-1, 7), (0, 0), (0, 1), (0, 2),
            ....: (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (1, 1), (1, 2),
            ....: (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (2, 0), (2, 2),
            ....: (2, 3), (2, 5), (2, 6), (2, 8)])
            sage: solution = H.self_surrounding(8)
            sage: G = sum([p.show2d() for p in solution], Graphics())                   # needs sage.plot

        ::

            sage: solution = H.self_surrounding(8, remove_incomplete_copies=False)
            sage: G = sum([p.show2d() for p in solution], Graphics())                   # needs sage.plot
        """
        # Define the box to tile
        minxyz, maxxyz = self.bounding_box()
        minxyz, maxxyz = vector(minxyz), vector(maxxyz)
        v = vector([radius for _ in range(self._dimension)])
        ranges = [range(a,b) for a,b in zip(minxyz-v, maxxyz+v)]
        box = Polyomino(itertools.product(*ranges))

        # Get the rows for this problem
        T = TilingSolver([self], box=box, reusable=True,
                reflection=True, rotation=True, outside=True)
        rows = T.rows()

        # Add one row to force the placement of the central tile
        coord_to_int = T.coord_to_int_dict()
        new_row = [coord_to_int[coord] for coord in self]
        new_row.append(len(coord_to_int))  # to force this row in the solution
        forced_row_number = len(rows)
        rows.append(new_row)

        # Construct the dancing links solver
        from sage.combinat.matrices.dancing_links import dlx_solver
        d = dlx_solver(rows)

        # Solve
        solution = d.one_solution(ncpus=ncpus)
        if solution is None:
            raise ValueError('No solution was found with radius={}, '
            'this tile can not be surrounded by itself'.format(radius))

        # Recover the polyominoes
        assert forced_row_number in solution
        solution.remove(forced_row_number)
        polyominoes = [T.row_to_polyomino(v) for v in solution]
        if remove_incomplete_copies:
            polyominoes = [p for p in polyominoes if len(p) == len(self)]

        # Recolor randomly the polyominoes
        from sage.plot.colors import Color
        from random import random
        for p in polyominoes:
            random_color = Color(tuple(random() for _ in range(3)))
            p.color(random_color)

        return polyominoes


#######################
# General tiling solver
#######################
class TilingSolver(SageObject):
    r"""
    Tiling solver.

    Solve the problem of tiling a polyomino with a certain number
    of polyominoes.

    INPUT:

    - ``pieces`` -- iterable of Polyominoes
    - ``box`` -- polyomino or tuple of integers (size of a box)
    - ``rotation`` -- boolean (default: ``True``); whether to allow
      rotations
    - ``reflection`` -- boolean (default: ``False``); whether to allow
      reflections
    - ``reusable`` -- boolean (default: ``False``); whether to allow
      the pieces to be reused
    - ``outside`` -- boolean (default: ``False``); whether to allow
      pieces to partially go outside of the box (all non-empty intersection
      of the pieces with the box are considered)

    EXAMPLES:

    By default, rotations are allowed and reflections are not allowed::

        sage: from sage.combinat.tiling import TilingSolver, Polyomino
        sage: p = Polyomino([(0,0,0)])
        sage: q = Polyomino([(0,0,0), (0,0,1)])
        sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
        sage: T = TilingSolver([p,q,r], box=(1,1,6))
        sage: T
        Tiling solver of 3 pieces into a box of size 6
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

    Tiling of a polyomino by polyominoes::

        sage: b = Polyomino([(0,0), (1,0), (1,1), (2,1), (1,2), (2,2), (0,3), (1,3)])
        sage: p = Polyomino([(0,0), (1,0)])
        sage: T = TilingSolver([p], box=b, reusable=True)
        sage: T.number_of_solutions()
        2

    TESTS::

        sage: T = TilingSolver([p,q,r], box=(1,1,6), rotation=False, reflection=True)
        Traceback (most recent call last):
        ...
        NotImplementedError: When reflection is allowed and rotation is not allowed
    """

    def __init__(self, pieces, box, rotation=True,
                 reflection=False, reusable=False, outside=False):
        r"""
        Constructor.

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: T
            Tiling solver of 3 pieces into a box of size 6
            Rotation allowed: True
            Reflection allowed: False
            Reusing pieces allowed: False
        """
        self._pieces = pieces
        self._free_module = self._pieces[0]._free_module
        if isinstance(box, Polyomino):
            self._box = box
        else:
            ranges = [range(a) for a in box]
            self._box = Polyomino(itertools.product(*ranges))
        self._rotation = rotation
        self._reflection = reflection
        if not self._rotation and self._reflection:
            raise NotImplementedError("When reflection is allowed and "
                                      "rotation is not allowed")
        self._reusable = reusable
        self._outside = outside

    def _repr_(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: TilingSolver([p,q,r], box=(1,1,6))
            Tiling solver of 3 pieces into a box of size 6
            Rotation allowed: True
            Reflection allowed: False
            Reusing pieces allowed: False
        """
        s = "Tiling solver of %s pieces " % len(self._pieces)
        s += "into a box of size %s\n" % len(self._box)
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
            return (sum(len(p) for p in self.pieces()) == len(self._box)
                    and len(self.rows()) != 0)

    def pieces(self):
        r"""
        Return the list of pieces.

        OUTPUT: list of 3d polyominoes

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
        Return an iterator over all the nonnegative integer coordinates
        contained in the space to tile.

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: list(T.space())
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 0, 4), (0, 0, 5)]
        """
        return iter(self._box)

    @cached_method
    def coord_to_int_dict(self):
        r"""
        Return a dictionary mapping coordinates to integers.

        OUTPUT: dictionary

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: A = T.coord_to_int_dict()
            sage: sorted(A.items())
            [((0, 0, 0), 3), ((0, 0, 1), 4), ((0, 0, 2), 5),
             ((0, 0, 3), 6), ((0, 0, 4), 7), ((0, 0, 5), 8)]

        Reusable pieces::

            sage: p = Polyomino([(0,0), (0,1)])
            sage: q = Polyomino([(0,0), (0,1), (1,0), (1,1)])
            sage: T = TilingSolver([p,q], box=[3,2], reusable=True)
            sage: B = T.coord_to_int_dict()
            sage: sorted(B.items())
            [((0, 0), 0), ((0, 1), 1), ((1, 0), 2), ((1, 1), 3),
             ((2, 0), 4), ((2, 1), 5)]
        """
        if self._reusable:
            return {c: i for i, c in enumerate(self.space())}

        number_of_pieces = len(self._pieces)
        return {c: i + number_of_pieces for i, c in enumerate(self.space())}

    @cached_method
    def int_to_coord_dict(self):
        r"""
        Return a dictionary mapping integers to coordinates.

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: B = T.int_to_coord_dict()
            sage: sorted(B.items())
            [(3, (0, 0, 0)), (4, (0, 0, 1)), (5, (0, 0, 2)),
             (6, (0, 0, 3)), (7, (0, 0, 4)), (8, (0, 0, 5))]

        Reusable pieces::

            sage: from sage.combinat.tiling import Polyomino, TilingSolver
            sage: p = Polyomino([(0,0), (0,1)])
            sage: q = Polyomino([(0,0), (0,1), (1,0), (1,1)])
            sage: T = TilingSolver([p,q], box=[3,2], reusable=True)
            sage: B = T.int_to_coord_dict()
            sage: sorted(B.items())
            [(0, (0, 0)), (1, (0, 1)), (2, (1, 0)),
             (3, (1, 1)), (4, (2, 0)), (5, (2, 1))]

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
            return dict(enumerate(self.space()))
        return dict(enumerate(self.space(), start=len(self._pieces)))

    @cached_method
    def rows_for_piece(self, i, mod_box_isometries=False):
        r"""
        Return the rows for the `i`-th piece.

        INPUT:

        - ``i`` -- integer; the `i`-th piece

        - ``mod_box_isometries`` -- boolean (default: ``False``); whether to
          consider only rows for positions up to the action of the
          quotient the group of isometries of the `n`-cube by the
          subgroup of isometries of the `a_1\times a_2\cdots \times a_n`
          rectangular box where are the `a_i` are assumed to be distinct.

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: T.rows_for_piece(0)
            [[0, 3], [0, 4], [0, 5], [0, 6], [0, 7], [0, 8]]
            sage: T.rows_for_piece(1)
            [[1, 3, 4], [1, 4, 5], [1, 5, 6], [1, 6, 7], [1, 7, 8]]
            sage: T.rows_for_piece(2)
            [[2, 3, 4, 5], [2, 4, 5, 6], [2, 5, 6, 7], [2, 6, 7, 8]]

        Less rows when using ``mod_box_isometries=True``::

            sage: a = Polyomino([(0,0,0), (0,0,1), (1,0,0)])
            sage: b = Polyomino([(0,0,0), (1,0,0), (0,1,0)])
            sage: T = TilingSolver([a,b], box=(2,1,3))
            sage: T.rows_for_piece(0)
            [[0, 2, 3, 5],
             [0, 3, 4, 6],
             [0, 2, 3, 6],
             [0, 3, 4, 7],
             [0, 2, 5, 6],
             [0, 3, 6, 7],
             [0, 3, 5, 6],
             [0, 4, 6, 7]]
            sage: T.rows_for_piece(0, mod_box_isometries=True)
            [[0, 2, 3, 5], [0, 3, 4, 6]]
            sage: T.rows_for_piece(1, mod_box_isometries=True)
            [[1, 2, 3, 5], [1, 3, 4, 6]]
        """
        p = self._pieces[i]
        if self._rotation:
            if self._reflection:
                orientation_preserving = False
            else:
                orientation_preserving = True
            if self._outside:
                it = p.isometric_copies_intersection(self._box,
                          orientation_preserving=orientation_preserving)
            else:
                it = p.isometric_copies(self._box,
                          orientation_preserving=orientation_preserving,
                          mod_box_isometries=mod_box_isometries)
        else:
            if self._reflection:
                raise NotImplementedError("Reflection allowed, Rotation not "
                                          "allowed is not implemented")
            else:
                if self._outside:
                    it = p.translated_copies_intersection(self._box)
                else:
                    it = p.translated_copies(self._box)
        coord_to_int = self.coord_to_int_dict()
        rows = []
        for q in it:
            L = [] if self._reusable else [i]
            L.extend(coord_to_int[coord] for coord in q)
            rows.append(L)
        return rows

    @cached_method
    def rows(self):
        r"""
        Creation of the rows.

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
            [1, 7, 8]
            [2, 3, 4, 5]
            [2, 4, 5, 6]
            [2, 5, 6, 7]
            [2, 6, 7, 8]
        """
        rows = []
        for i in range(len(self._pieces)):
            rows.extend(self.rows_for_piece(i))
        return rows

    def _rows_mod_box_isometries(self, i):
        r"""
        Return a list of rows representing the solutions up to isometries of
        the box.

        The positions of the `i`-th pieces are chosen up to isometries of
        the box. In dimension 3, there are four times less rows for that
        piece.

        It is currently implemented only when the pieces are not reusable.

        INPUT:

        - ``i`` -- integer; the `i`-th piece to consider, that piece must not
          be isometric to itself by a isometry that preserve the box

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0), (1,0,0), (1,1,0), (1,0,1), (2,0,1)],
            ....:               color='red')
            sage: T = TilingSolver([p], box=(3,4,2))
            sage: T._rows_mod_box_isometries(0)
            [[0, 1, 3, 4, 11, 13],
             [0, 3, 5, 6, 13, 15],
             [0, 9, 11, 12, 19, 21],
             [0, 11, 13, 14, 21, 23],
             [0, 1, 9, 10, 11, 18],
             [0, 3, 11, 12, 13, 20],
             [0, 5, 13, 14, 15, 22],
             [0, 2, 3, 4, 5, 11],
             [0, 4, 5, 6, 7, 13],
             [0, 10, 11, 12, 13, 19],
             [0, 12, 13, 14, 15, 21],
             [0, 2, 9, 10, 12, 20],
             [0, 4, 11, 12, 14, 22],
             [0, 6, 13, 14, 16, 24]]

        We test that there are four times less rows for that polyomino::

            sage: len(T.rows()) == 4 * len(T._rows_mod_box_isometries(0))
            True

        Now, a real use case. A solution of the game Quantumino is a tiling
        of a `5 \times 8 \times 2` box. Since a `5 \times 8 \times 2` box
        has four orientation preserving isometries, each solution up to
        rotation is counted four times by this dancing links solver::

            sage: from sage.games.quantumino import QuantuminoSolver
            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: q = QuantuminoSolver(0)
            sage: T = q.tiling_solver()
            sage: dlx_solver(T.rows())                  # long time (10s)
            Dancing links solver for 96 columns and 5484 rows

        It is possible to avoid to compute 4 times each solution up to
        rotations. This is done by choosing a piece (here the 0th) and
        considering 4 times less positions for that piece. To be precise,
        90 positions instead of 360, therefore the dancing links solver
        below has 270 less rows::

            sage: dlx_solver(T._rows_mod_box_isometries(0))  # long time (10s)
            Dancing links solver for 96 columns and 5214 rows
        """
        assert not self._reusable, ("this code assumes the pieces are not reusable")
        len_pieces = len(self._pieces)
        if not 0 <= i < len_pieces:
            raise ValueError("i(={}) must be 0 <= i < {}".format(i,len_pieces))
        rows = []
        for j in range(len_pieces):
            if j == i:
                rows.extend(self.rows_for_piece(j, mod_box_isometries=True))
            else:
                rows.extend(self.rows_for_piece(j))
        return rows

    def nrows_per_piece(self):
        r"""
        Return the number of rows necessary by each piece.

        OUTPUT: list

        EXAMPLES::

            sage: from sage.games.quantumino import QuantuminoSolver
            sage: q = QuantuminoSolver(0)
            sage: T = q.tiling_solver()
            sage: T.nrows_per_piece()             # long time (10s)
            [360, 360, 360, 360, 360, 180, 180, 672, 672, 360, 360, 180, 180, 360, 360, 180]
        """
        return [len(self.rows_for_piece(i)) for i in range(len(self._pieces))]

    def starting_rows(self):
        r"""
        Return the starting rows for each piece.

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: T.starting_rows()
            [0, 6, 11, 15]
        """
        s = 0
        S = [s]
        for a in self.nrows_per_piece():
            s += a
            S.append(s)
        return S

    def row_to_polyomino(self, row_number):
        r"""
        Return a polyomino associated to a row.

        INPUT:

        - ``row_number`` -- integer; the `i`-th row

        OUTPUT: polyomino

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: a = Polyomino([(0,0,0), (0,0,1), (1,0,0)], color='blue')
            sage: b = Polyomino([(0,0,0), (1,0,0), (0,1,0)], color='red')
            sage: T = TilingSolver([a,b], box=(2,1,3))
            sage: len(T.rows())
            16

        ::

            sage: T.row_to_polyomino(7)
            Polyomino: [(0, 0, 2), (1, 0, 1), (1, 0, 2)], Color: blue

        ::

            sage: T.row_to_polyomino(13)
            Polyomino: [(0, 0, 1), (1, 0, 1), (1, 0, 2)], Color: red

        TESTS:

        We check that issue :issue:`32252` is fixed and that colors of
        polyominoes are properly recovered::

            sage: v = Polyomino([(0, 0), (0, 1)], color='blue')
            sage: h = Polyomino([(0, 0), (1, 0)], color='red')
            sage: T = TilingSolver(pieces=[v, h], box=(2, 2),
            ....:                  rotation=False, reflection=False, reusable=True)
            sage: for i in range(4): print(i,T.row_to_polyomino(i))
            0 Polyomino: [(0, 0), (0, 1)], Color: blue
            1 Polyomino: [(1, 0), (1, 1)], Color: blue
            2 Polyomino: [(0, 0), (1, 0)], Color: red
            3 Polyomino: [(0, 1), (1, 1)], Color: red
        """
        rows = self.rows()
        row = rows[row_number]
        if self._reusable:
            if row_number < 0:
                row_number += len(rows)
            from bisect import bisect
            no = bisect(self.starting_rows(), row_number) - 1
            indices = row
        else:
            no = row[0]
            indices = row[1:]
        int_to_coord = self.int_to_coord_dict()
        coords = [int_to_coord[i] for i in indices]
        color = self._pieces[no].color()
        return Polyomino(coords, color=color)

    def dlx_solver(self):
        r"""
        Return the sage DLX solver of that tiling problem.

        OUTPUT: dLX Solver

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: T.dlx_solver()
            Dancing links solver for 9 columns and 15 rows
        """
        from sage.combinat.matrices.dancing_links import dlx_solver
        return dlx_solver(self.rows())

    def _dlx_solutions_iterator(self):
        r"""
        Return an iterator over the row indices of the solutions.

        OUTPUT: iterator

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: list(T._dlx_solutions_iterator())
            [[0, 7, 14], [0, 12, 10], [6, 13, 5], [6, 14, 2], [11, 9, 5], [11, 10, 3]]
        """
        if len(self.rows()) == 0:
            return

        x = self.dlx_solver()
        while x.search() == 1:
            yield x.get_solution()

    def _dlx_common_prefix_solutions_iterator(self):
        r"""
        Return an iterator over the row indices of solutions and of partial
        solutions, i.e. the common prefix of two consecutive solutions.

        The purpose is to illustrate the backtracking and construct an
        animation of the evolution of solutions.

        OUTPUT: iterator

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: list(T._dlx_solutions_iterator())
            [[0, 7, 14], [0, 12, 10], [6, 13, 5], [6, 14, 2], [11, 9, 5], [11, 10, 3]]
            sage: list(T._dlx_common_prefix_solutions_iterator())
            [[0, 7, 14], [0], [0, 12, 10], [], [6, 13, 5], [6], [6, 14, 2], [], [11, 9, 5], [11], [11, 10, 3]]

        ::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: y = Polyomino([(0,0),(1,0),(2,0),(3,0),(2,1)], color='yellow')
            sage: T = TilingSolver([y], box=(5,10), reusable=True, reflection=True)
            sage: for a in T._dlx_common_prefix_solutions_iterator(): a
            [0, 83, 114, 43, 158, 5, 128, 183, 168, 25]
            [0, 83, 114, 43, 158]
            [0, 83, 114, 43, 158, 33, 128, 183, 104, 25]
            [0, 83, 114]
            [0, 83, 114, 100, 52, 183, 128, 33, 95, 47]
            [0, 83]
            [0, 83, 178, 15, 158, 5, 128, 183, 168, 25]
            [0, 83, 178, 15, 158]
            [0, 83, 178, 15, 158, 33, 128, 183, 104, 25]
            []
            [56, 1, 113, 15, 159, 34, 155, 182, 168, 24]
            [56, 1, 113]
            [56, 1, 113, 164, 51, 118, 155, 34, 96, 47]
            [56, 1, 113, 164, 51]
            [56, 1, 113, 164, 51, 182, 155, 34, 96, 19]
            [56]
            [56, 29, 113, 100, 51, 118, 155, 34, 96, 47]
            [56, 29, 113, 100, 51]
            [56, 29, 113, 100, 51, 182, 155, 34, 96, 19]
        """
        it = self._dlx_solutions_iterator()
        B = next(it)
        while True:
            yield B
            try:
                A, B = B, next(it)
            except StopIteration:
                return
            common_prefix = []
            for a, b in zip(A, B):
                if a == b:
                    common_prefix.append(a)
                else:
                    break
            yield common_prefix

    def _dlx_incremental_solutions_iterator(self):
        r"""
        Return an iterator over the row indices of the incremental
        solutions.

        Between two incremental solution, either one piece is added or one
        piece is removed.

        The purpose is to illustrate the backtracking and construct an
        animation of the evolution of solutions.

        OUTPUT: iterator

        EXAMPLES::

            sage: from sage.combinat.tiling import TilingSolver, Polyomino
            sage: p = Polyomino([(0,0,0)])
            sage: q = Polyomino([(0,0,0), (0,0,1)])
            sage: r = Polyomino([(0,0,0), (0,0,1), (0,0,2)])
            sage: T = TilingSolver([p,q,r], box=(1,1,6))
            sage: list(T._dlx_solutions_iterator())
            [[0, 7, 14], [0, 12, 10], [6, 13, 5], [6, 14, 2], [11, 9, 5], [11, 10, 3]]
            sage: list(T._dlx_incremental_solutions_iterator())
            [[0, 7, 14], [0, 7], [0], [0, 12], [0, 12, 10], [0, 12], [0], [], [6], [6, 13], [6, 13, 5], [6, 13], [6], [6, 14], [6, 14, 2], [6, 14], [6], [], [11], [11, 9], [11, 9, 5], [11, 9], [11], [11, 10], [11, 10, 3]]

        ::

            sage: y = Polyomino([(0,0),(1,0),(2,0),(3,0),(2,1)], color='yellow')
            sage: T = TilingSolver([y], box=(5,10), reusable=True, reflection=True)
            sage: for a in T._dlx_solutions_iterator(): a
            [0, 83, 114, 43, 158, 5, 128, 183, 168, 25]
            [0, 83, 114, 43, 158, 33, 128, 183, 104, 25]
            [0, 83, 114, 100, 52, 183, 128, 33, 95, 47]
            [0, 83, 178, 15, 158, 5, 128, 183, 168, 25]
            [0, 83, 178, 15, 158, 33, 128, 183, 104, 25]
            [56, 1, 113, 15, 159, 34, 155, 182, 168, 24]
            [56, 1, 113, 164, 51, 118, 155, 34, 96, 47]
            [56, 1, 113, 164, 51, 182, 155, 34, 96, 19]
            [56, 29, 113, 100, 51, 118, 155, 34, 96, 47]
            [56, 29, 113, 100, 51, 182, 155, 34, 96, 19]
            sage: len(list(T._dlx_incremental_solutions_iterator()))
            123
        """
        it = self._dlx_solutions_iterator()
        B = next(it)
        while True:
            yield B
            try:
                A, B = B, next(it)
            except StopIteration:
                return
            common_prefix = 0
            for a, b in zip(A, B):
                if a == b:
                    common_prefix += 1
                else:
                    break
            for i in range(1, len(A)-common_prefix):
                yield A[:-i]
            for j in range(common_prefix, len(B)):
                yield B[:j]

    def solve(self, partial=None):
        r"""
        Return an iterator of list of polyominoes that are an exact cover
        of the box.

        INPUT:

        - ``partial`` -- string (default: ``None``); whether to
          include partial (incomplete) solutions. It can be one of the
          following:

          - ``None`` -- include only complete solution
          - ``'common_prefix'`` -- common prefix between two consecutive solutions
          - ``'incremental'`` -- one piece change at a time

        OUTPUT: iterator of list of polyominoes

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
            return
        if partial is None:
            it = self._dlx_solutions_iterator()
        elif partial == 'common_prefix':
            it = self._dlx_common_prefix_solutions_iterator()
        elif partial == 'incremental':
            it = self._dlx_incremental_solutions_iterator()
        else:
            raise ValueError("Unknown value for partial (=%s)" % partial)
        for solution in it:
            yield [self.row_to_polyomino(v) for v in solution]

    def number_of_solutions(self):
        r"""
        Return the number of distinct solutions.

        OUTPUT: integer

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

        - ``partial`` -- string (default: ``None``); whether to
          include partial (incomplete) solutions. It can be one of the
          following:

          - ``None`` -- include only complete solutions
          - ``'common_prefix'`` -- common prefix between two consecutive solutions
          - ``'incremental'`` -- one piece change at a time

        - ``stop`` -- integer (default: ``None``); number of frames

        - ``size`` -- number (default: ``0.75``); the size of each
          ``1 \times 1`` square. This does a homothety with respect
          to the center of each polyomino.

        - ``axes`` -- boolean (default: ``False``); whether the x and
          y axes are shown

        EXAMPLES::

            sage: from sage.combinat.tiling import Polyomino, TilingSolver
            sage: y = Polyomino([(0,0),(1,0),(2,0),(3,0),(2,1)], color='cyan')
            sage: T = TilingSolver([y], box=(5,10), reusable=True, reflection=True)
            sage: a = T.animate()                                                       # needs sage.plot
            sage: a                             # long time, optional - imagemagick, needs sage.plot
            Animation with 10 frames

        Include partial solutions (common prefix between two consecutive
        solutions)::

            sage: a = T.animate('common_prefix')                                        # needs sage.plot
            sage: a                             # long time, optional - imagemagick, needs sage.plot
            Animation with 19 frames

        Incremental solutions (one piece removed or added at a time)::

            sage: a = T.animate('incremental')  # long time (2s)                        # needs sage.plot
            sage: a                             # long time (2s), optional - imagemagick, needs sage.plot
            Animation with 123 frames

        ::

            sage: a.show()                      # long time, optional - imagemagick, needs sage.plot

        The ``show`` function takes arguments to specify the delay between
        frames (measured in hundredths of a second, default value 20) and
        the number of iterations (default: 0, which means to iterate
        forever). To iterate 4 times with half a second between each frame::

            sage: a.show(delay=50, iterations=4)        # long time, optional - imagemagick, needs sage.plot

        Limit the number of frames::

            sage: a = T.animate('incremental', stop=13); a      # not tested            # needs sage.plot
            Animation with 13 frames
        """
        dimension = self._box._dimension
        if dimension == 2:
            from sage.plot.graphics import Graphics
            from sage.plot.animate import Animation
            it = self.solve(partial=partial)
            it = itertools.islice(it, stop)
            L = [sum([piece.show2d(size)
                      for piece in solution], Graphics()) for solution in it]
            (xmin,ymin), (xmax,ymax) = self._box.bounding_box()
            xmax = xmax+0.5
            ymax = ymax+0.5
            a = Animation(L, xmin=xmin-0.5, ymin=ymin-0.5,
                          xmax=xmax, ymax=ymax, aspect_ratio=1, axes=axes)
            return a
        elif dimension == 3:
            raise NotImplementedError("3d Animation must be implemented "
                                      "in Jmol first")
        else:
            raise NotImplementedError("Dimension must be 2 or 3 in order "
                                      "to make an animation")
