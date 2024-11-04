r"""
Family Games America's Quantumino solver

This module allows to solve the `Quantumino puzzle
<http://familygamesamerica.com/mainsite/consumers/productview.php?pro_id=274&search=quantumino>`_
made by Family Games America (see also `this video
<http://www.youtube.com/watch?v=jX_VKzakZi8>`_ on Youtube). This puzzle was
left at the dinner room of the Laboratoire de Combinatoire Informatique
Mathématique in Montreal by Franco Saliola during winter 2011.

The solution uses the dancing links code which is in Sage and is based on
the more general code available in the module :mod:`sage.combinat.tiling`.
Dancing links were originally introduced by Donald Knuth in 2000
(:arxiv:`cs/0011047`). In particular,
Knuth used dancing links to solve tilings of a region by 2D pentaminos.
Here we extend the method for 3D pentaminos.

This module defines two classes :

- :class:`sage.games.quantumino.QuantuminoState` class, to represent a
  state of the Quantumino game, i.e. a solution or a partial solution.

- :class:`sage.games.quantumino.QuantuminoSolver` class, to find, enumerate
  and count the number of solutions of the Quantumino game where one of the
  piece is put aside.

AUTHOR:

- Sébastien Labbé, April 28th, 2011

DESCRIPTION (from [1]):

    "
    Pentamino games have been taken to a whole different level; a 3-D
    level, with this colorful creation! Using the original pentamino
    arrangements of 5 connected squares which date from 1907, players are
    encouraged to "think inside the box" as they try to fit 16 of the 17
    3-D pentamino pieces inside the playing perimeters. Remove a different
    piece each time you play for an entirely new challenge! Thousands of
    solutions to be found!
    Quantumino hands-on educational tool where players learn how shapes
    can be transformed or arranged into predefined shapes and spaces.
    Includes:
    1 wooden frame, 17 wooden blocks, instruction booklet.
    Age: 8+
    "

EXAMPLES:

Here are the 17 wooden blocks of the Quantumino puzzle numbered from 0 to 16 in
the following 3d picture. They will show up in 3D in your default (=Jmol)
viewer::

    sage: from sage.games.quantumino import show_pentaminos
    sage: show_pentaminos()
    Graphics3d Object

To solve the puzzle where the pentamino numbered 12 is put aside::

    sage: from sage.games.quantumino import QuantuminoSolver
    sage: s = next(QuantuminoSolver(12).solve())          # long time (10 s)
    sage: s                                               # long time (<1s)
    Quantumino state where the following pentamino is put aside :
    Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)], Color: blue
    sage: s.show3d()                                      # long time (<1s)
    Graphics3d Object

To remove the frame::

    sage: s.show3d().show(frame=False)                    # long time (<1s)

To solve the puzzle where the pentamino numbered 7 is put aside::

    sage: s = next(QuantuminoSolver(7).solve())           # long time (10 s)
    sage: s                                               # long time (<1s)
    Quantumino state where the following pentamino is put aside :
    Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange
    sage: s.show3d()                                      # long time (<1s)
    Graphics3d Object

The solution is iterable. This may be used to explicitly list the positions of each
pentamino::

    sage: for p in s: p                                   # long time (<1s)
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...
    Polyomino: [(...), (...), (...), (...), (...)], Color: ...

To get all the solutions, use the iterator returned by the ``solve``
method. Note that finding the first solution is the most time consuming
because it needs to create the complete data to describe the problem::

    sage: it = QuantuminoSolver(7).solve()
    sage: next(it)                                     # not tested (10s)
    Quantumino state where the following pentamino is put aside :
    Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange
    sage: next(it)                                     # not tested (0.001s)
    Quantumino state where the following pentamino is put aside :
    Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange
    sage: next(it)                                     # not tested (0.001s)
    Quantumino state where the following pentamino is put aside :
    Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange

To get the solution inside other boxes::

    sage: s = next(QuantuminoSolver(7, box=(4,4,5)).solve())        # not tested (2s)
    sage: s.show3d()                                                # not tested (<1s)

::

    sage: s = next(QuantuminoSolver(7, box=(2,2,20)).solve())       # not tested (1s)
    sage: s.show3d()                                                # not tested (<1s)

If there are no solution, a StopIteration error is raised::

    sage: next(QuantuminoSolver(7, box=(3,3,3)).solve())
    Traceback (most recent call last):
    ...
    StopIteration

The implementation allows a lot of introspection. From the
:class:`~sage.combinat.tiling.TilingSolver` object,
it is possible to retrieve the rows that are passed to the DLX
solver and count them. It is also possible to get an instance of the DLX
solver to play with it::

    sage: q = QuantuminoSolver(0)
    sage: T = q.tiling_solver()
    sage: T
    Tiling solver of 16 pieces into a box of size 80
    Rotation allowed: True
    Reflection allowed: False
    Reusing pieces allowed: False
    sage: rows = T.rows()                            # not tested (10 s)
    sage: len(rows)                                  # not tested (but fast)
    5484
    sage: x = T.dlx_solver()                         # long time (10 s)
    sage: x                                          # long time (fast)
    Dancing links solver for 96 columns and 5484 rows

TESTS:

We check that all pentaminos are equal to their canonical translate::

    sage: from sage.games.quantumino import pentaminos
    sage: all(p == p.canonical() for p in pentaminos)
    True

REFERENCES:

- [1] `Family Games America's Quantumino
  <http://familygamesamerica.com/mainsite/consumers/productview.php?pro_id=274&search=quantumino>`_
- [2] `Quantumino - How to Play <http://www.youtube.com/watch?v=jX_VKzakZi8>`_ on Youtube
- [3] Knuth, Donald (2000). *Dancing links*. :arxiv:`cs/0011047`.
"""
# ****************************************************************************
#       Copyright (C) 2011 Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.all", "Graphics")
lazy_import("sage.plot.plot3d.platonic", "cube")
lazy_import("sage.plot.plot3d.shapes2", "text3d")
from sage.modules.free_module_element import vector
from sage.combinat.tiling import Polyomino, TilingSolver

################################################
# Example:  The family games america: Quantumino
################################################
pentaminos = []
pentaminos.append(Polyomino([(0,0,0), (1,0,0), (1,1,0), (1,2,0), (1,1,1)], color='deeppink'))
pentaminos.append(Polyomino([(0,0,0), (1,0,0), (1,0,1), (2,0,0), (2,1,0)], color='deeppink'))
pentaminos.append(Polyomino([(0,0,0), (1,0,0), (1,1,0), (1,2,0), (0,0,1)], color='green'))
pentaminos.append(Polyomino([(0,0,0), (0,1,0), (0,2,0), (1,0,0), (1,0,1)], color='green'))
pentaminos.append(Polyomino([(0,1,0), (1,0,1), (1,1,0), (1,1,1), (1,2,0)], color='red'))
pentaminos.append(Polyomino([(0,0,0), (1,0,0), (1,1,0), (1,0,1), (2,0,1)], color='red'))
pentaminos.append(Polyomino([(0,0,0), (1,0,0), (1,1,0), (1,2,0), (1,2,1)], color='orange'))
pentaminos.append(Polyomino([(0,0,0), (1,0,0), (0,1,0), (0,2,0), (0,2,1)], color='orange'))
pentaminos.append(Polyomino([(0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1)], color='yellow'))
pentaminos.append(Polyomino([(0,0,0), (1,0,0), (1,1,0), (1,1,1), (0,0,1)], color='yellow'))
pentaminos.append(Polyomino([(0,0,0), (0,1,0), (1,1,0), (0,2,0), (1,1,1)], color='midnightblue'))
pentaminos.append(Polyomino([(0,0,0), (1,0,0), (1,1,0), (1,0,1), (1,2,0)], color='darkblue'))
pentaminos.append(Polyomino([(0,0,0), (1,0,0), (1,1,0), (1,1,1), (2,1,1)], color='blue'))
pentaminos.append(Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,1,1), (1,2,1)], color='blue'))
pentaminos.append(Polyomino([(0,0,0), (1,0,0), (1,1,0), (2,1,0), (2,1,1)], color='purple'))
pentaminos.append(Polyomino([(0,0,0), (0,1,0), (1,1,0), (1,2,0), (1,2,1)], color='purple'))
pentaminos.append(Polyomino([(0,1,0), (1,0,0), (1,1,0), (1,1,1), (1,2,0)], color='gray'))


def show_pentaminos(box=(5,8,2)):
    r"""
    Show the 17 3-D pentaminos included in the game and the `5 \times 8
    \times 2` box where 16 of them must fit.

    INPUT:

    - ``box`` -- tuple of size three (default: ``(5,8,2)``),
      size of the box

    OUTPUT: 3D Graphic object

    EXAMPLES::

        sage: from sage.games.quantumino import show_pentaminos
        sage: show_pentaminos()    # not tested (1s)

    To remove the frame do::

        sage: show_pentaminos().show(frame=False)  # not tested (1s)
    """
    G = Graphics()
    for i, p in enumerate(pentaminos):
        x = 4 * (i % 4)
        y = 4 * (i // 4)
        q = p + (x, y, 0)
        G += q.show3d()
        G += text3d(str(i), (x, y, 2))
    G += cube(color='gray',opacity=0.5).scale(box).translate((17, 6, 0))

    # hack to set the aspect ratio to 1
    a, b = G.bounding_box()
    a, b = map(vector, (a, b))
    G.frame_aspect_ratio(tuple(b - a))

    return G

##############################
# Class QuantuminoState
##############################


class QuantuminoState(SageObject):
    r"""
    A state of the Quantumino puzzle.

    Used to represent a solution or a partial solution of the Quantumino
    puzzle.

    INPUT:

    - ``pentos`` -- list of 16 3d pentamino representing the (partial)
      solution
    - ``aside`` -- 3d polyomino, the unused 3D pentamino
    - ``box`` -- tuple of size three (default: ``(5,8,2)``),
      size of the box

    EXAMPLES::

        sage: from sage.games.quantumino import pentaminos, QuantuminoState
        sage: p = pentaminos[0]
        sage: q = pentaminos[5]
        sage: r = pentaminos[11]
        sage: S = QuantuminoState([p,q], r)
        sage: S
        Quantumino state where the following pentamino is put aside :
        Polyomino: [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 2, 0)], Color: darkblue

    ::

        sage: from sage.games.quantumino import QuantuminoSolver
        sage: next(QuantuminoSolver(3).solve())      # not tested (1.5s)
        Quantumino state where the following pentamino is put aside :
        Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (1, 0, 0), (1, 0, 1)], Color: green
    """
    def __init__(self, pentos, aside, box=(5,8,2)):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import pentaminos, QuantuminoState
            sage: p = pentaminos[0]
            sage: q = pentaminos[5]
            sage: QuantuminoState([p], q)
            Quantumino state where the following pentamino is put aside :
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (2, 0, 1)], Color: red
        """
        assert all(isinstance(p, Polyomino) for p in pentos), "pentos must be an iterable of Polyomino"
        assert isinstance(aside, Polyomino), "aside must be a Polyomino"
        self._pentos = pentos
        self._aside = aside
        self._box = box

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import pentaminos, QuantuminoState
            sage: p = pentaminos[0]
            sage: q = pentaminos[5]
            sage: QuantuminoState([p], q)
            Quantumino state where the following pentamino is put aside :
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (2, 0, 1)], Color: red
        """
        return "Quantumino state where the following pentamino is put aside :\n%s" % self._aside

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from sage.games.quantumino import pentaminos, QuantuminoState
            sage: p = pentaminos[0]
            sage: q = pentaminos[5]
            sage: r = pentaminos[11]
            sage: S = QuantuminoState([p,q], r)
            sage: for a in S: a
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (2, 0, 1)], Color: red
        """
        return iter(self._pentos)

    def list(self):
        r"""
        Return the list of 3d polyomino making the solution.

        EXAMPLES::

            sage: from sage.games.quantumino import pentaminos, QuantuminoState
            sage: p = pentaminos[0]
            sage: q = pentaminos[5]
            sage: r = pentaminos[11]
            sage: S = QuantuminoState([p,q], r)
            sage: L = S.list()
            sage: L[0]
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            sage: L[1]
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (2, 0, 1)], Color: red
        """
        return list(self)

    def show3d(self, size=0.85):
        r"""
        Return the solution as a 3D Graphic object.

        OUTPUT: 3D Graphic Object

        EXAMPLES::

            sage: from sage.games.quantumino import QuantuminoSolver
            sage: s = next(QuantuminoSolver(0).solve())    # not tested (1.5s)
            sage: G = s.show3d()                            # not tested (<1s)
            sage: type(G)                                   # not tested
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>

        To remove the frame::

            sage: G.show(frame=False) # not tested

        To see the solution with Tachyon viewer::

            sage: G.show(viewer='tachyon', frame=False) # not tested
        """
        G = Graphics()
        for p in self:
            G += p.show3d(size=size)
        aside_pento = self._aside.canonical() + (2,-4,0)
        G += aside_pento.show3d(size=size)

        # the box to fill
        half_box = tuple(a/2 for a in self._box)
        b = cube(color='gray',opacity=0.2).scale(self._box).translate(half_box)
        b = b.translate((0, -.5, -.5))
        G += b

        # hack to set the aspect ratio to 1
        a,b = G.bounding_box()
        a,b = map(vector, (a,b))
        G.frame_aspect_ratio(tuple(b-a))

        return G

##############################
# Class QuantuminoSolver
##############################


class QuantuminoSolver(SageObject):
    r"""
    Return the Quantumino solver for the given box where one of the
    pentamino is put aside.

    INPUT:

    - ``aside`` -- integer, from 0 to 16, the aside pentamino
    - ``box`` -- tuple of size three (default: ``(5,8,2)``),
      size of the box

    EXAMPLES::

        sage: from sage.games.quantumino import QuantuminoSolver
        sage: QuantuminoSolver(9)
        Quantumino solver for the box (5, 8, 2)
        Aside pentamino number: 9
        sage: QuantuminoSolver(12, box=(5,4,4))
        Quantumino solver for the box (5, 4, 4)
        Aside pentamino number: 12
    """
    def __init__(self, aside, box=(5, 8, 2)):
        r"""
        Constructor.

        EXAMPLES::

            sage: from sage.games.quantumino import QuantuminoSolver
            sage: QuantuminoSolver(9)
            Quantumino solver for the box (5, 8, 2)
            Aside pentamino number: 9
        """
        if not (0 <= aside < 17):
            raise ValueError(f"aside (={aside}) must be between 0 and 16")
        self._aside = aside
        self._box = box

    def __repr__(self):
        r"""
        String representation.

        EXAMPLES::

            sage: from sage.games.quantumino import QuantuminoSolver
            sage: QuantuminoSolver(0)
            Quantumino solver for the box (5, 8, 2)
            Aside pentamino number: 0
        """
        s = "Quantumino solver for the box %s\n" % (self._box, )
        s += "Aside pentamino number: %s" % self._aside
        return s

    def tiling_solver(self):
        r"""
        Return the Tiling solver of the Quantumino Game where one of the
        pentamino is put aside.

        EXAMPLES::

            sage: from sage.games.quantumino import QuantuminoSolver
            sage: QuantuminoSolver(0).tiling_solver()
            Tiling solver of 16 pieces into a box of size 80
            Rotation allowed: True
            Reflection allowed: False
            Reusing pieces allowed: False
            sage: QuantuminoSolver(14).tiling_solver()
            Tiling solver of 16 pieces into a box of size 80
            Rotation allowed: True
            Reflection allowed: False
            Reusing pieces allowed: False
            sage: QuantuminoSolver(14, box=(5,4,4)).tiling_solver()
            Tiling solver of 16 pieces into a box of size 80
            Rotation allowed: True
            Reflection allowed: False
            Reusing pieces allowed: False
        """
        pieces = pentaminos[:self._aside] + pentaminos[self._aside+1:]
        return TilingSolver(pieces, box=self._box)

    def solve(self, partial=None):
        r"""
        Return an iterator over the solutions where one of the pentamino is
        put aside.

        INPUT:

        - ``partial`` -- string (default: ``None``), whether to
          include partial (incomplete) solutions. It can be one of the
          following:

          - ``None`` -- include only complete solution
          - ``'common'`` -- common part between two consecutive solutions
          - ``'incremental'`` -- one piece change at a time

        OUTPUT: iterator of :class:`QuantuminoState`

        EXAMPLES:

        Get one solution::

            sage: from sage.games.quantumino import QuantuminoSolver
            sage: s = next(QuantuminoSolver(8).solve())          # long time (9s)
            sage: s                                              # long time (fast)
            Quantumino state where the following pentamino is put aside :
            Polyomino: [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (1, 1, 0)], Color: yellow
            sage: s.show3d()                                     # long time (< 1s)
            Graphics3d Object

        The explicit solution::

            sage: for p in s: p                                  # long time (fast)
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...
            Polyomino: [(...), (...), (...), (...), (...)], Color: ...


        Enumerate the solutions::

            sage: it = QuantuminoSolver(0).solve()
            sage: next(it)                                          # not tested
            Quantumino state where the following pentamino is put aside :
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink
            sage: next(it)                                          # not tested
            Quantumino state where the following pentamino is put aside :
            Polyomino: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 2, 0)], Color: deeppink

        With the partial solutions included, one can see the evolution
        between consecutive solutions (an animation would be better)::

            sage: it = QuantuminoSolver(0).solve(partial='common')
            sage: next(it).show3d()               # not tested (2s)
            sage: next(it).show3d()               # not tested (< 1s)
            sage: next(it).show3d()               # not tested (< 1s)

        Generalizations of the game inside different boxes::

            sage: next(QuantuminoSolver(7, (4,4,5)).solve())       # long time (2s)
            Quantumino state where the following pentamino is put aside :
            Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange
            sage: next(QuantuminoSolver(7, (2,2,20)).solve())      # long time (1s)
            Quantumino state where the following pentamino is put aside :
            Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 2, 1), (1, 0, 0)], Color: orange
            sage: next(QuantuminoSolver(3, (2,2,20)).solve())      # long time (1s)
            Quantumino state where the following pentamino is put aside :
            Polyomino: [(0, 0, 0), (0, 1, 0), (0, 2, 0), (1, 0, 0), (1, 0, 1)], Color: green

        If the volume of the box is not 80, there is no solution::

            sage: next(QuantuminoSolver(7, box=(3,3,9)).solve())
            Traceback (most recent call last):
            ...
            StopIteration

        If the box is too small, there is no solution::

            sage: next(QuantuminoSolver(4, box=(40,2,1)).solve())
            Traceback (most recent call last):
            ...
            StopIteration
        """
        T = self.tiling_solver()
        aside = pentaminos[self._aside]
        for pentos in T.solve(partial=partial):
            yield QuantuminoState(pentos, aside, self._box)

    def number_of_solutions(self):
        r"""
        Return the number of solutions.

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.games.quantumino import QuantuminoSolver
            sage: QuantuminoSolver(4, box=(3,2,2)).number_of_solutions()
            0

        This computation takes several days::

            sage: QuantuminoSolver(0).number_of_solutions()                # not tested
            ??? hundreds of millions ???
        """
        return self.tiling_solver().number_of_solutions()
