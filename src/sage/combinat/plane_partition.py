# -*- coding: utf-8 -*-
r"""
Plane Partitions

AUTHORS:

- Jang Soo Kim (2016): Initial implementation
- Jessica Striker (2016): Added additional methods
"""
# ****************************************************************************
#       Copyright (C) 2016 Jang Soo Kim <jangsookim@skku.edu>,
#                     2016 Jessica Striker <jessicapalencia@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations
from typing import Iterator

from sage.structure.list_clone import ClonableArray
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.posets.posets import Poset
from sage.rings.integer import Integer
from sage.misc.misc_c import prod
from sage.combinat.tableau import Tableau
from sage.plot.plot3d.platonic import cube


class PlanePartition(ClonableArray,
        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A plane partition.

    A *plane partition* is a stack of cubes in the positive orthant.

    INPUT:

    - ``PP`` -- a list of lists which represents a tableau

    - ``box_size`` -- (optional) a list ``[A, B, C]`` of 3 positive integers,
      where ``A``, ``B``, ``C`` are the lengths of the box in the `x`-axis,
      `y`-axis, `z`-axis, respectively; if this is not given, it is
      determined by the smallest box bounding ``PP``

    OUTPUT:

    The plane partition whose tableau representation is ``PP``.

    EXAMPLES::

        sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
        sage: PP
        Plane partition [[4, 3, 3, 1], [2, 1, 1], [1, 1]]

    TESTS::

        sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
        sage: TestSuite(PP).run()
    """
    @staticmethod
    def __classcall_private__(cls, PP, box_size=None):
        """
        Construct a plane partition with the appropriate parent.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1], [2,1,1], [1,1]])
            sage: PP.parent() is PlanePartitions((3,4,4))
            True
        """
        if box_size is None:
            if PP:
                box_size = (len(PP), len(PP[0]), PP[0][0])
            else:
                box_size = (0, 0, 0)
        return PlanePartitions(box_size)(PP)

    def __init__(self, parent, PP, check=True):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: TestSuite(PP).run()
        """
        ClonableArray.__init__(self, parent, PP, check=check)
        self._max_x = parent._box[0]
        self._max_y = parent._box[1]
        self._max_z = parent._box[2]

    def check(self):
        """
        Check to see that ``self`` is a valid plane partition.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.check()
        """
        if len(self) == 0:
            return
        if len(self) > self.parent()._box[0]:
            raise ValueError("too big in z direction")
        if len(self[0]) > self.parent()._box[1]:
            raise ValueError("too big in y direction")
        if self[0][0] > self.parent()._box[2]:
            raise ValueError("too big in x direction")

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            Plane partition [[4, 3, 3, 1], [2, 1, 1], [1, 1]]
        """
        return "Plane partition {}".format(list(self))

    def to_tableau(self) -> Tableau:
        r"""
        Return the tableau class of ``self``.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.to_tableau()
            [[4, 3, 3, 1], [2, 1, 1], [1, 1]]
        """
        return Tableau(self)  # type: ignore

    def z_tableau(self):
        r"""
        Return the projection of ``self`` in the `z` direction.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.z_tableau()
            [[4, 3, 3, 1], [2, 1, 1, 0], [1, 1, 0, 0]]
        """
        Z = [[0 for i in range(self._max_y)] for j in range(self._max_x)]
        for C in self.cells():
            Z[C[0]][C[1]] += 1
        return Z

    def y_tableau(self):
        r"""
        Return the projection of ``self`` in the `y` direction.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.y_tableau()
            [[4, 3, 2], [3, 1, 0], [3, 0, 0], [1, 0, 0]]
        """
        Y = [[0 for i in range(self._max_x)] for j in range(self._max_z)]
        for C in self.cells():
            Y[C[2]][C[0]] += 1
        return Y

    def x_tableau(self):
        r"""
        Return the projection of ``self`` in the `x` direction.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.x_tableau()
            [[3, 2, 1, 1], [3, 1, 1, 0], [2, 1, 1, 0], [1, 0, 0, 0]]
        """
        X = [[0 for i in range(self._max_z)] for j in range(self._max_y)]
        for C in self.cells():
            X[C[1]][C[2]] += 1
        return X

    def cells(self) -> list:
        r"""
        Return the list of cells inside ``self``.

        EXAMPLES::

            sage: PP = PlanePartition([[3,1],[2]])
            sage: PP.cells()
            [[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, 0], [1, 0, 0], [1, 0, 1]]
        """
        L = []
        for r in range(len(self)):
            for c in range(len(self[r])):
                for h in range(self[r][c]):
                    L.append([r, c, h])
        return L

    def _repr_diagram(self, show_box=False, use_unicode=False) -> str:
        r"""
        Return a string of the 3D diagram of ``self``.

        INPUT:

        - ``show_box`` -- boolean (default: ``False``); if ``True``,
          also shows the visible tiles on the `xy`-, `yz`-, `zx`-planes
        - ``use_unicode`` -- boolean (default: ``False``); use unicode

        OUTPUT:

        A string of the 3D diagram of the plane partition.

        EXAMPLES::

            sage: print(PlanePartition([[4,3,3,1],[2,1,1],[1,1]])._repr_diagram())
                    __
                   /\_\
                __/\/_/
             __/\_\/\_\
            /\_\/_/\/\_\
            \/\_\_\/\/_/
             \/_/\_\/_/
                \/_/\_\
                   \/_/
            sage: print(PlanePartition([[4,3,3,1],[2,1,1],[1,1]])._repr_diagram(True))
                ______
               /_/_/\_\
              /_/_/\/_/\
             /_/\_\/\_\/\
            /\_\/_/\/\_\/\
            \/\_\_\/\/_/\/
             \/_/\_\/_/\/
              \_\/_/\_\/
               \_\_\/_/
        """
        x = self._max_x
        y = self._max_y
        z = self._max_z

        drawing = [[" " for i in range(2 * x + y + z)]
                   for j in range(y + z + 1)]

        hori = u"_" if use_unicode else "_"
        down = u"╲" if use_unicode else "\\"
        up = u"╱" if use_unicode else "/"

        def superpose(l, c, letter):
            # add the given letter at line l and column c
            exist = drawing[l][c]
            if exist == " " or exist == "_":
                drawing[l][c] = letter

        def add_topside(i, j, k):
            X = z + j - k
            Y = 2 * x - 2 * i + j + k
            superpose(X, Y - 2, hori)
            superpose(X, Y - 1, hori)
            superpose(X + 1, Y - 2, down)
            superpose(X + 1, Y - 1, hori)
            superpose(X + 1, Y, down)

        def add_rightside(i, j, k):
            X = z + j - k
            Y = 2 * x - 2 * i + j + k
            superpose(X - 1, Y - 1, hori)
            superpose(X - 1, Y, hori)
            superpose(X, Y - 2, up)
            superpose(X, Y - 1, hori)
            superpose(X, Y, up)

        def add_leftside(i, j, k):
            X = z + j - k
            Y = 2 * x - 2 * i + j + k
            superpose(X, Y, up)
            superpose(X, Y + 1, down)
            superpose(X + 1, Y + 1, up)
            superpose(X + 1, Y, down)

        tab = self.z_tableau()
        for r in range(len(tab)):
            for c in range(len(tab[r])):
                if tab[r][c] > 0 or show_box:
                    add_topside(r, c, tab[r][c])

        tab = self.y_tableau()
        for r in range(len(tab)):
            for c in range(len(tab[r])):
                if self.y_tableau()[r][c] > 0 or show_box:
                    add_rightside(c, tab[r][c], r)

        tab = self.x_tableau()
        for r in range(len(tab)):
            for c in range(len(tab[r])):
                if self.x_tableau()[r][c] > 0 or show_box:
                    add_leftside(tab[r][c], r, c)

        check = not show_box
        while check:
            if drawing and all(char == " " for char in drawing[-1]):
                drawing.pop()
            else:
                check = False

        if not drawing:
            return u"∅" if use_unicode else ""

        if use_unicode:
            return u'\n'.join(u"".join(s for s in row) for row in drawing)
        return '\n'.join("".join(s for s in row) for row in drawing)

    def _ascii_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: ascii_art(PP)
                    __
                   /\_\
                __/\/_/
             __/\_\/\_\
            /\_\/_/\/\_\
            \/\_\_\/\/_/
             \/_/\_\/_/
                \/_/\_\
                   \/_/
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self._repr_diagram().splitlines(), baseline=0)

    def _unicode_art_(self):
        r"""
        Return a unicode representation of ``self``.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: unicode_art(PP)
                    __
                   ╱╲_╲
                __╱╲╱_╱
             __╱╲_╲╱╲_╲
            ╱╲_╲╱_╱╲╱╲_╲
            ╲╱╲_╲_╲╱╲╱_╱
             ╲╱_╱╲_╲╱_╱
                ╲╱_╱╲_╲
                   ╲╱_╱
        """
        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt(self._repr_diagram(use_unicode=True).splitlines(), baseline=0)

    def pp(self, show_box=False):
        r"""
        Return a pretty print of the plane partition.

        INPUT:

        - ``show_box`` -- boolean (default: ``False``); if ``True``,
          also shows the visible tiles on the `xy`-, `yz`-, `zx`-planes

        OUTPUT:

        A pretty print of the plane partition.

        EXAMPLES::

            sage: PlanePartition([[4,3,3,1],[2,1,1],[1,1]]).pp()
                    __
                   /\_\
                __/\/_/
             __/\_\/\_\
            /\_\/_/\/\_\
            \/\_\_\/\/_/
             \/_/\_\/_/
                \/_/\_\
                   \/_/
            sage: PlanePartition([[4,3,3,1],[2,1,1],[1,1]]).pp(True)
                ______
               /_/_/\_\
              /_/_/\/_/\
             /_/\_\/\_\/\
            /\_\/_/\/\_\/\
            \/\_\_\/\/_/\/
             \/_/\_\/_/\/
              \_\/_/\_\/
               \_\_\/_/
        """
        print(self._repr_diagram(show_box))

    def _latex_(self, show_box=False,
                colors=["white", "lightgray", "darkgray"]) -> str:
        r"""
        Return latex code for ``self``, which uses TikZ package to draw
        the plane partition.

        INPUT:

        - ``show_box`` -- boolean (default: ``False``); if ``True``,
          also shows the visible tiles on the `xy`-, `yz`-, `zx`-planes

        - ``colors`` -- (default: ``["white", "lightgray", "darkgray"]``)
          list ``[A, B, C]`` of 3 strings representing colors

        OUTPUT:

        Latex code for drawing the plane partition.

        EXAMPLES::

            sage: PP = PlanePartition([[1]])
            sage: latex(PP)
            \begin{tikzpicture}
            \draw[fill=white,shift={(210:0)},shift={(-30:0)},shift={(90:1)}]
            (0,0)--(-30:1)--(0,-1)--(210:1)--(0,0);
            \draw[fill=darkgray,shift={(210:0)},shift={(-30:1)},shift={(90:0)}]
            (0,0)--(210:1)--(150:1)--(0,1)--(0,0);
            \draw[fill=lightgray,shift={(210:1)},shift={(-30:0)},shift={(90:0)}]
            (0,0)--(0,1)--(30:1)--(-30:1)--(0,0);
            \end{tikzpicture}
        """
        from sage.graphs.graph_latex import setup_latex_preamble
        setup_latex_preamble()

        ret = "\\begin{tikzpicture}\n"

        def add_topside(i, j, k):
            return "\\draw[fill={},shift={{(210:{})}},shift={{(-30:{})}},shift={{(90:{})}}]\n(0,0)--(-30:1)--(0,-1)--(210:1)--(0,0);\n".format(colors[0], i, j, k)

        def add_leftside(j, k, i):
            return "\\draw[fill={},shift={{(210:{})}},shift={{(-30:{})}},shift={{(90:{})}}]\n(0,0)--(0,1)--(30:1)--(-30:1)--(0,0);\n".format(colors[1], i, j, k)

        def add_rightside(k, i, j):
            return "\\draw[fill={},shift={{(210:{})}},shift={{(-30:{})}},shift={{(90:{})}}]\n(0,0)--(210:1)--(150:1)--(0,1)--(0,0);\n".format(colors[2], i, j, k)
        funcs = [add_topside, add_rightside, add_leftside]
        tableaux = [self.z_tableau(), self.y_tableau(), self.x_tableau()]
        for i in range(3):
            f = funcs[i]
            tab = tableaux[i]
            for r in range(len(tab)):
                for c in range(len(tab[r])):
                    if tab[r][c] > 0 or show_box:
                        ret += f(r, c, tab[r][c])
        return ret + "\\end{tikzpicture}"

    def plot(self, show_box=False, colors=None):
        r"""
        Return a plot of ``self``.

        INPUT:

        - ``show_box`` -- boolean (default: ``False``); if ``True``,
          also shows the visible tiles on the `xy`-, `yz`-, `zx`-planes

        - ``colors`` -- (default: ``["white", "lightgray", "darkgray"]``)
          list ``[A, B, C]`` of 3 strings representing colors

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.plot()
            Graphics object consisting of 27 graphics primitives
        """
        from sage.functions.trig import cos, sin
        from sage.plot.polygon import polygon
        from sage.symbolic.constants import pi
        from sage.plot.plot import plot
        if colors is None:
            colors = ["white", "lightgray", "darkgray"]
        Uside = [[0, 0], [cos(-pi / 6), sin(-pi / 6)],
                 [0, -1], [cos(7 * pi / 6), sin(7 * pi / 6)]]
        Lside = [[0, 0], [cos(-pi / 6), sin(-pi / 6)],
                 [cos(pi / 6), sin(pi / 6)], [0, 1]]
        Rside = [[0, 0], [0, 1], [cos(5 * pi / 6), sin(5 * pi / 6)],
                 [cos(7 * pi / 6), sin(7 * pi / 6)]]
        Xdir = [cos(7 * pi / 6), sin(7 * pi / 6)]
        Ydir = [cos(-pi / 6), sin(-pi / 6)]
        Zdir = [0, 1]

        def move(side, i, j, k):
            return [[P[0] + i * Xdir[0] + j * Ydir[0] + k * Zdir[0],
                     P[1] + i * Xdir[1] + j * Ydir[1] + k * Zdir[1]]
                    for P in side]

        def add_topside(i, j, k):
            return polygon(move(Uside, i, j, k), edgecolor="black",
                           color=colors[0])

        def add_leftside(i, j, k):
            return polygon(move(Lside, i, j, k), edgecolor="black",
                           color=colors[1])

        def add_rightside(i, j, k):
            return polygon(move(Rside, i, j, k), edgecolor="black",
                           color=colors[2])
        TP = plot([])
        for r in range(len(self.z_tableau())):
            for c in range(len(self.z_tableau()[r])):
                if self.z_tableau()[r][c] > 0 or show_box:
                    TP += add_topside(r, c, self.z_tableau()[r][c])
        for r in range(len(self.y_tableau())):
            for c in range(len(self.y_tableau()[r])):
                if self.y_tableau()[r][c] > 0 or show_box:
                    TP += add_rightside(c, self.y_tableau()[r][c], r)
        for r in range(len(self.x_tableau())):
            for c in range(len(self.x_tableau()[r])):
                if self.x_tableau()[r][c] > 0 or show_box:
                    TP += add_leftside(self.x_tableau()[r][c], r, c)
        TP.axes(show=False)
        return TP

    def plot3d(self, colors=None):
        r"""
        Return a 3D-plot of ``self``.

        INPUT:

        - ``colors`` -- (default: ``["white", "lightgray", "darkgray"]``)
          list ``[A, B, C]`` of 3 strings representing colors

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.plot3d()
            Graphics3d Object
        """
        if colors is None:
            colors = ["white", "lightgray", "darkgray"]
        return sum(cube(c, color=colors, frame_thickness=2,
                        frame_color='black', frame=False)
                   for c in self.cells())

    def complement(self, tableau_only=False):
        r"""
        Return the complement of ``self``.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.complement()
            Plane partition [[4, 4, 3, 3], [4, 3, 3, 2], [3, 1, 1, 0]]
            sage: PP.complement(True)
            [[4, 4, 3, 3], [4, 3, 3, 2], [3, 1, 1, 0]]
        """
        A = self._max_x
        B = self._max_y
        C = self._max_z
        T = [[C for i in range(B)] for j in range(A)]
        z_tab = self.z_tableau()
        for r in range(A):
            for c in range(B):
                T[A - 1 - r][B - 1 - c] = C - z_tab[r][c]
        if tableau_only:
            return T
        else:
            return type(self)(self.parent(), T, check=False)

    def transpose(self, tableau_only=False):
        r"""
        Return the transpose of ``self``.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.transpose()
            Plane partition [[4, 2, 1], [3, 1, 1], [3, 1, 0], [1, 0, 0]]
            sage: PP.transpose(True)
            [[4, 2, 1], [3, 1, 1], [3, 1, 0], [1, 0, 0]]
        """
        T = [[0 for i in range(self._max_x)] for j in range(self._max_y)]
        z_tab = self.z_tableau()
        for r in range(len(z_tab)):
            for c in range(len(z_tab[r])):
                T[c][r] = z_tab[r][c]
        if tableau_only:
            return T
        else:
            return type(self)(self.parent(), T, check=False)

    def is_SPP(self) -> bool:
        r"""
        Return whether ``self`` is a symmetric plane partition.

        A plane partition is symmetric if the corresponding tableau is
        symmetric about the diagonal.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_SPP()
            False
            sage: PP = PlanePartition([[3,3,2],[3,3,2],[2,2,2]])
            sage: PP.is_SPP()
            True
            sage: PP = PlanePartition([[3,2,1],[2,0,0]])
            sage: PP.is_SPP()
            False
            sage: PP = PlanePartition([[3,2,0],[2,0,0]])
            sage: PP.is_SPP()
            True
            sage: PP = PlanePartition([[3,2],[2,0],[1,0]])
            sage: PP.is_SPP()
            False
            sage: PP = PlanePartition([[3,2],[2,0],[0,0]])
            sage: PP.is_SPP()
            True
        """
        Z = self.z_tableau()
        c1 = len(Z)
        c2 = len(Z[0])
        size = max(c1, c2)
        T = [[0 for i in range(size)] for j in range(size)]
        for i in range(c1):
            for j in range(c2):
                T[i][j] = Z[i][j]
        return all(T[r][c] == T[c][r]
            for r in range(size)
            for c in range(r, size))

    def is_CSPP(self) -> bool:
        r"""
        Return whether ``self`` is a cyclically symmetric plane partition.

        A plane partition is cyclically symmetric if its `x`, `y`, and `z`
        tableaux are all equal.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_CSPP()
            False
            sage: PP = PlanePartition([[3,2,2],[3,1,0],[1,1,0]])
            sage: PP.is_CSPP()
            True
        """
        if self.z_tableau() == self.y_tableau():
            return True
        return False

    def is_TSPP(self) -> bool:
        r"""
        Return whether ``self`` is a totally symmetric plane partition.

        A plane partition is totally symmetric if it is both symmetric and
        cyclically symmetric.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_TSPP()
            False
            sage: PP = PlanePartition([[3,3,3],[3,3,2],[3,2,1]])
            sage: PP.is_TSPP()
            True
        """
        return self.is_CSPP() and self.is_SPP()

    def is_SCPP(self) -> bool:
        r"""
        Return whether ``self`` is a self-complementary plane partition.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_SCPP()
            False
            sage: PP = PlanePartition([[4,4,4,4],[4,4,2,0],[4,2,0,0],[0,0,0,0]])
            sage: PP.is_SCPP()
            True
        """
        return self.z_tableau() == self.complement(True)

    def is_TCPP(self) -> bool:
        r"""
        Return whether ``self`` is a transpose-complementary plane partition.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_TCPP()
            False
            sage: PP = PlanePartition([[4,4,3,2],[4,4,2,1],[4,2,0,0],[2,0,0,0]])
            sage: PP.is_TCPP()
            True
        """
        return self.transpose(True) == self.complement(True)

    def is_SSCPP(self) -> bool:
        r"""
        Return whether ``self`` is a symmetric, self-complementary
        plane partition.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_SSCPP()
            False
            sage: PP = PlanePartition([[4,3,3,2],[3,2,2,1],[3,2,2,1],[2,1,1,0]])
            sage: PP.is_SSCPP()
            True
            sage: PP = PlanePartition([[2,1],[1,0]])
            sage: PP.is_SSCPP()
            True
            sage: PP = PlanePartition([[4,3,2],[3,2,1],[2,1,0]])
            sage: PP.is_SSCPP()
            True
            sage: PP = PlanePartition([[4,2,2,2],[2,2,2,2],[2,2,2,2],[2,2,2,0]])
            sage: PP.is_SSCPP()
            True
        """
        return self.is_SPP() and self.is_SCPP()

    def is_CSTCPP(self) -> bool:
        r"""
        Return whether ``self`` is a cyclically symmetric and
        transpose-complementary plane partition.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_CSTCPP()
            False
            sage: PP = PlanePartition([[4,4,3,2],[4,3,2,1],[3,2,1,0],[2,1,0,0]])
            sage: PP.is_CSTCPP()
            True
        """
        return self.is_CSPP() and self.is_TCPP()

    def is_CSSCPP(self) -> bool:
        r"""
        Return whether ``self`` is a cyclically symmetric and
        self-complementary plane partition.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_CSSCPP()
            False
            sage: PP = PlanePartition([[4,4,4,1],[3,3,2,1],[3,2,1,1],[3,0,0,0]])
            sage: PP.is_CSSCPP()
            True
        """
        return self.is_CSPP() and self.is_SCPP()

    def is_TSSCPP(self) -> bool:
        r"""
        Return whether ``self`` is a totally symmetric self-complementary
        plane partition.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_TSSCPP()
            False
            sage: PP = PlanePartition([[4,4,3,2],[4,3,2,1],[3,2,1,0],[2,1,0,0]])
            sage: PP.is_TSSCPP()
            True
        """
        return self.is_TSPP() and self.is_SCPP()


PP = PlanePartition


class PlanePartitions(UniqueRepresentation, Parent):
    r"""
    All plane partitions inside a rectangular box of given side lengths.

    INPUT:

    - ``box_size`` -- a triple of positive integers indicating the size
      of the box containing the plane partition

    EXAMPLES:

    This will create an instance to manipulate the plane partitions
    in a `4 \times 3 \times 2` box::

        sage: P = PlanePartitions((4,3,2))
        sage: P
        Plane partitions inside a 4 x 3 x 2 box
        sage: P.cardinality()
        490

    .. SEEALSO::

        :class:`PlanePartition`
    """
    @staticmethod
    def __classcall_private__(cls, box_size):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: P1 = PlanePartitions((4,3,2))
            sage: P2 = PlanePartitions([4,3,2])
            sage: P1 is P2
            True
        """
        return super(PlanePartitions, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        r"""
        Initialize ``self``

        EXAMPLES::

            sage: PP = PlanePartitions((4,3,2))
            sage: TestSuite(PP).run()
        """
        if len(box_size) != 3:
            raise ValueError("invalid box size")
        self._box = box_size
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PlanePartitions((4,3,2))
            Plane partitions inside a 4 x 3 x 2 box
        """
        return "Plane partitions inside a {} x {} x {} box".format(
            self._box[0], self._box[1], self._box[2])

    def __iter__(self) -> Iterator[PP]:
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((1,2,1)))
            [Plane partition [[0, 0]],
             Plane partition [[1, 0]],
             Plane partition [[1, 1]]]
        """
        A = self._box[0]
        B = self._box[1]
        C = self._box[2]
        from sage.combinat.tableau import SemistandardTableaux as SST
        for T in SST([B for i in range(A)], max_entry=C + A):  # type:ignore
            PP = [[0 for _ in range(B)] for _ in range(A)]
            for r in range(A):
                for c in range(B):
                    PP[A - 1 - r][B - 1 - c] = T[r][c] - r - 1
            yield self.element_class(self, PP, check=False)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of plane partitions inside an `a \times b \times c`
        box is equal to

        .. MATH::

            \prod_{i=1}^{a} \prod_{j=1}^{b} \prod_{k=1}^{c}
            \frac{i+j+k-1}{i+j+k-2}.

        EXAMPLES::

            sage: P = PlanePartitions((4,3,5))
            sage: P.cardinality()
            116424
        """
        A = self._box[0]
        B = self._box[1]
        C = self._box[2]
        return Integer(prod(Integer(i + j + k - 1) / Integer(i + j + k - 2)
                            for i in range(1, A + 1)
                            for j in range(1, B + 1)
                            for k in range(1, C + 1)))

    def box(self) -> tuple:
        """
        Return the sizes of the box of the plane partitions of ``self``
        are contained in.

        EXAMPLES::

            sage: P = PlanePartitions((4,3,5))
            sage: P.box()
            (4, 3, 5)
        """
        return self._box

    def random_element(self) -> PP:
        r"""
        Return a uniformly random element of ``self``.

        ALGORITHM:

        This uses the
        :meth:`~sage.combinat.posets.posets.FinitePoset.random_order_ideal`
        method and the natural bijection with plane partitions.

        EXAMPLES::

            sage: P = PlanePartitions((4,3,5))
            sage: p = P.random_element()
            sage: p.parent() is P
            True
        """
        def leq(thing1, thing2):
            return all(thing1[i] <= thing2[i] for i in range(len(thing1)))
        elem = [(i, j, k) for i in range(self._box[0])
                for j in range(self._box[1])
                for k in range(self._box[2])]
        myposet = Poset((elem, leq))
        R = myposet.random_order_ideal()
        Z = [[0 for i in range(self._box[1])] for j in range(self._box[0])]
        for C in R:
            Z[C[0]][C[1]] += 1
        return self.element_class(self, Z, check=False)

    Element = PlanePartition
