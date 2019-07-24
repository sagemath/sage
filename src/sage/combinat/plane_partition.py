# -*- coding: utf-8 -*-
r"""
Plane Partitions

AUTHORS:

- Jang Soo Kim (2016): Initial implementation
- Jessica Striker (2016): Added additional methods
"""
#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function, absolute_import
from six.moves import range
from six import add_metaclass

from sage.structure.list_clone import ClonableArray
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.posets.posets import Poset
from sage.rings.integer import Integer
from sage.misc.all import prod
from sage.combinat.tableau import Tableau


#@add_metaclass(InheritComparisonClasscallMetaclass)
class PlanePartition(ClonableArray):
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
    def __classcall_private__(cls, PP):
        """
        Construct a plane partition with the appropriate parent.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1], [2,1,1], [1,1]])
            sage: PP.parent() is PlanePartitions((3,4,4))
            True
        """
	if isinstance(PP,PlanePartition):
            return PP
        if box_size is None:
            if PP:
                box_size = (len(PP), len(PP[0]), PP[0][0])
            else:
                box_size = (0, 0, 0)
        return PlanePartitions(box_size, symmetry=symmetry)(PP)

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

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            Plane partition [[4, 3, 3, 1], [2, 1, 1], [1, 1]]
        """
        return "Plane partition {}".format(list(self))

    def to_tableau(self):
        r"""
        Return the tableau class of ``self``.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.to_tableau()
            [[4, 3, 3, 1], [2, 1, 1], [1, 1]]
        """
        return Tableau(self)

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

    def cells(self):
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
                    L.append([r,c,h])
        return L

    def _repr_diagram(self, show_box=False, use_unicode=False):
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
            superpose(X, Y-2, hori)
            superpose(X, Y-1, hori)
            superpose(X + 1, Y-2, down)
            superpose(X + 1, Y-1, hori)
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

    def _latex_(self, show_box=False, colors=["white","lightgray","darkgray"]):
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

        def add_topside(i,j,k):
            return "\\draw[fill={},shift={{(210:{})}},shift={{(-30:{})}},shift={{(90:{})}}]\n(0,0)--(-30:1)--(0,-1)--(210:1)--(0,0);\n".format(colors[0],i,j,k)

        def add_leftside(j,k,i):
            return "\\draw[fill={},shift={{(210:{})}},shift={{(-30:{})}},shift={{(90:{})}}]\n(0,0)--(0,1)--(30:1)--(-30:1)--(0,0);\n".format(colors[1],i,j,k)

        def add_rightside(k,i,j):
            return "\\draw[fill={},shift={{(210:{})}},shift={{(-30:{})}},shift={{(90:{})}}]\n(0,0)--(210:1)--(150:1)--(0,1)--(0,0);\n".format(colors[2],i,j,k)
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

    def plot(self, show_box=False, colors=["white","lightgray","darkgray"]):
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
        Uside = [[0,0], [cos(-pi/6),sin(-pi/6)], [0,-1], [cos(7*pi/6),sin(7*pi/6)]]
        Lside = [[0,0], [cos(-pi/6),sin(-pi/6)], [cos(pi/6),sin(pi/6)], [0,1]]
        Rside = [[0,0], [0,1], [cos(5*pi/6),sin(5*pi/6)], [cos(7*pi/6),sin(7*pi/6)]]
        Xdir = [cos(7*pi/6), sin(7*pi/6)]
        Ydir = [cos(-pi/6), sin(-pi/6)]
        Zdir = [0, 1]

        def move(side, i, j, k):
            return [[P[0]+i*Xdir[0]+j*Ydir[0]+k*Zdir[0],
                     P[1]+i*Xdir[1]+j*Ydir[1]+k*Zdir[1]]
                    for P in side]

        def add_topside(i, j, k):
            return polygon(move(Uside,i,j,k), edgecolor="black", color=colors[0])
        def add_leftside(i, j, k):
            return polygon(move(Lside,i,j,k), edgecolor="black", color=colors[1])
        def add_rightside(i, j, k):
            return polygon(move(Rside,i,j,k), edgecolor="black", color=colors[2])
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
                T[A-1-r][B-1-c] = C - z_tab[r][c]
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

    def is_SPP(self):
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
        """
        z_tab = self.z_tableau()
        return all(z_tab[r][c] == z_tab[c][r]
                   for r in range(len(z_tab))
                   for c in range(r, len(z_tab[r])))

    def is_CSPP(self):
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

    def is_TSPP(self):
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

    def is_SCPP(self):
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

    def is_TCPP(self):
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

    def is_SSCPP(self):
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
        """
        return self.is_SPP() and self.is_SCPP()

    def is_CSTCPP(self):
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

    def is_CSSCPP(self):
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

    def is_TSSCPP(self):
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


class PlanePartitions(UniqueRepresentation, Parent):
    r"""
    A factory class for plane partitions.

    PlanePartitions([a,b,c]) returns the class of plane partitions that fit
    inside an a \times b \times c box.

    Optional keyword is 'symmetry'.

    Describe options.

    """
    @staticmethod
    def __classcall_private__(cls, box_size, symmetry=None):
        """
        Return correct parent based on input.

        EXAMPLES::

            sage: 1+1
            2
        """
        if len(box_size) != 3:
            raise ValueError("invalid box size")
        if symmetry == None:
            return PlanePartitions_all(box_size)
        elif symmetry == 'TSPP':
            return PlanePartitions_TSPP(box_size)
        elif symmetry == 'SPP':
            return PlanePartitions_SPP(box_size)
        elif symmetry == 'CSPP':
            return PlanePartitions_CSPP(box_size)
        elif symmetry == 'TSSCPP':
            return PlanePartitions_TSSCPP(box_size)
        else:
            raise ValueError("invalid symmetry class option; must be None, 'TSPP', 'SPP', 'CSPP', or 'TSCCPP' ")




    Element = PlanePartition


    def __init__(self, box_size, symmetry=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = PlanePartitions()
            sage: TestSuite(S).run()  # long time
        """
        self.symmetry=symmetry
        self._box = box_size
 #       Tableaux.__init__(self, **kwds)
 #       ClonableArray.__init__(self)


class PlanePartitions_all(PlanePartitions):
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
        return super(PlanePartitions_all, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        r"""
        Initialize ``self``

        EXAMPLES::

            sage: PP = PlanePartitions((4,3,2))
            sage: TestSuite(PP).run()
        """
        super(PlanePartitions_all,self ).__init__(box_size)
        self._box = box_size

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PlanePartitions((4,3,2))
            Plane partitions inside a 4 x 3 x 2 box
        """
        return "Plane partitions inside a %s x %s x %s box" % (self._box[0], self._box[1], self._box[2])
#        return "Plane partitions inside a box"

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((1,2,1)))
            [Plane partition [[0, 0]],
             Plane partition [[1, 0]],
             Plane partition [[1, 1]]]
        """
#        A = self._box[0]
#        B = self._box[1]
#        C = self._box[2]
#        from sage.combinat.tableau import SemistandardTableaux
#        for T in SemistandardTableaux([B for i in range(A)], max_entry=C+A):
#            PP = [[0 for i in range(B)] for j in range(A)]
#            for r in range(A):
#                for c in range(B):
#                    PP[A-1-r][B-1-c] = T[r][c] - r - 1
#            yield self.element_class(self, PP, check=False)
        def componentwise_comparer(thing1,thing2):
            if len(thing1) == len(thing2):
                if all(thing1[i] <= thing2[i] for i in range(len(thing1))):
                    return True
            return False
        def product_of_chains_poset(list_of_chain_lengths):
            elem = cartesian_product([range(chain_length) for chain_length in list_of_chain_lengths])
            return Poset((elem, componentwise_comparer))

        a = self._box[0]
        b = self._box[1]
        c = self._box[2]

        pocp = product_of_chains_poset([a,b,c])

        matrixList = [] #list of all PlaneParitions with parameters(a,b,c)

        #iterate through each antichain of product of chains poset with paramaters (a,b,c)
        for acl in pocp.antichains_iterator():
            ppMatrix = [[0] * (c) for i in range(b)] #creates a matrix for the plane parition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]

            #ac format ex: [x,y,z]
            #iterate through each antichain, assigning the y,z position in ppMatrix = the height of the stack (x + 1)
            for ac in acl:
                x = ac[0]
                y = ac[1]
                z = ac[2]
                ppMatrix[y][z] = (x+1)

            #for each value in current antichain, fill in the rest of the matrix by rule M[y,z] = Max(M[y+1,z], M[y,z+1]) antichiain is now in plane partitian format
            if acl != []:
                for i in range(b):
                    i = b-(i+1)
                    for j in range(c):
                        j = c-(j+1)
                        if (ppMatrix[i][j] == 0):
                            iValue = 0
                            jValue = 0
                            if i < b-1:
                                iValue = ppMatrix[i+1][j]
                            if j < c-1:
                                jValue = ppMatrix[i][j+1]
                            ppMatrix[i][j] = max(iValue,jValue)

            matrixList.append(ppMatrix) #add PlanePartition to list of plane partitions

        matrixList.sort()
        current = 0
        while current < len(matrixList):
            yield self.element_class(self, matrixList[current])
            current += 1


    def cardinality(self):
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
        return Integer(prod( Integer(i+j+k-1) / Integer(i+j+k-2)
                             for i in range(1, A+1) for j in range(1, B+1)
                             for k in range(1, C+1) ))

    def box(self):
        """
        Return the sizes of the box of the plane partitions of ``self``
        are contained in.

        EXAMPLES::

            sage: P = PlanePartitions((4,3,5))
            sage: P.box()
            (4, 3, 5)
        """
        return self._box

    def random_element(self):
        r"""
        Return a uniformly random element of ``self``.

        ALGORITHM:

        This uses the
        :meth:`~sage.combinat.posets.posets.FinitePoset.random_order_ideal`
        method and the natural bijection with plane partitions.

        EXAMPLES::

            sage: P = PlanePartitions((4,3,5))
            sage: P.random_element()
            Plane partition [[4, 3, 3], [4, 0, 0], [2, 0, 0], [0, 0, 0]]
        """
        def leq(thing1, thing2):
            return all(thing1[i] <= thing2[i] for i in range(len(thing1)))
        elem = [(i,j,k) for i in range(self._box[0]) for j in range(self._box[1])
                for k in range(self._box[2])]
        myposet = Poset((elem, leq))
        R = myposet.random_order_ideal()
        Z = [[0 for i in range(self._box[1])] for j in range(self._box[0])]
        for C in R:
            Z[C[0]][C[1]] += 1
        return self.element_class(self, Z, check=False)

class PlanePartitions_TSPP(PlanePartitions):

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
        return super(PlanePartitions_TSPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([3,3,3], symmetry=TSPP)
            sage: TestSuite(PP).run()
        """
        super(PlanePartitions_TSPP, self).__init__(box_size=box_size)
        self._box=box_size

    def _repr_(self):
        return " Totally Symmetric Plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])


    def __iter__(self):
        """
        EXAMPLES:
        


        """
        def componentwise_comparer(thing1,thing2):
            if len(thing1) == len(thing2):
                if all(thing1[i] <= thing2[i] for i in range(len(thing1))):
                    return True
            return False
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        pl = []
        for x in range(0,a):
            for y in range(x, b):
                    for z in range(y,c):
                        pl.append((x,y,z))

        pocp = Poset((pl,componentwise_comparer))
        matrixList = [] #list of all PlaneParitions with parameters(a,b,c)
        #iterate through each antichain of product of chains poset with paramaters (a,b,c)
        for acl in pocp.antichains_iterator():
            ppMatrix = [[0] * (c) for i in range(b)] #creates a matrix for the plane parition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]

            #ac format ex: [x,y,z]
            for ac in acl:
                x = ac[0]
                y = ac[1]
                z = ac[2]
                ppMatrix[y][z] = (x+1) #x,y,z
                ppMatrix[z][x] = (y+1) #y,z,x
                ppMatrix[x][y] = (z+1) #z,x,y
                ppMatrix[z][y] = (x+1) #x,z,y
                ppMatrix[x][z] = (y+1) #y,x,z
                ppMatrix[y][x] = (z+1) #z,y,x


            #for each value in current antichain, fill in the rest of the matrix by rule M[y,z] = Max(M[y+1,z], M[y,z+1]) antichiain is now in plane partitian format
            if acl != []:
                for i in range(b):
                    i = b-(i+1)
                    for j in range(c):
                        j = c-(j+1)
                        if (ppMatrix[i][j] == 0):
                            iValue = 0
                            jValue = 0
                            if i < b-1:
                                iValue = ppMatrix[i+1][j]
                            if j < c-1:
                                jValue = ppMatrix[i][j+1]
                            ppMatrix[i][j] = max(iValue,jValue)

            matrixList.append(ppMatrix) #add PlanePartition to list of plane partitions

        matrixList.sort()


        current = 0
        while current < len(matrixList):
            yield self.element_class(self, matrixList[current])
            current += 1



class PlanePartitions_SPP(PlanePartitions):

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
        return super(PlanePartitions_SPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([3,3,3], symmetry=TSPP)
            sage: TestSuite(PP).run()
        """
        if box_size[1] != box_size[2]:
            raise ValueError("invalid box size")
        super(PlanePartitions_SPP, self).__init__(box_size=box_size)
        self._box=box_size

    def _repr_(self):
        return " Symmetric Plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])


    def __iter__(self):
        def componentwise_comparer(thing1,thing2):
            if len(thing1) == len(thing2):
                if all(thing1[i] <= thing2[i] for i in range(len(thing1))):
                    return True
            return False
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        pl = []
        for x in range(0,a):
            for y in range(0, b):
                    for z in range(0,c):
                        if z <= y:
                            pl.append((x,y,z))


        pocp = Poset((pl,componentwise_comparer))

        matrixList = [] #list of all PlaneParitions with parameters(a,b,c)
            #iterate through each antichain of product of chains poset with paramaters (a,b,c)
        for acl in pocp.antichains_iterator():
            ppMatrix = [[0] * (c) for i in range(b)] #creates a matrix for the plane parition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]
                #ac format ex: [x,y,z]
                #iterate through each antichain, assigning the y,z position in ppMatrix = the height of the stack (x + 1)
            for ac in acl:
                x = ac[0]
                y = ac[1]
                z = ac[2]
                ppMatrix[y][z] = (x+1)

                #for each value in current antichain, fill in the rest of the matrix by rule M[y,z] = Max(M[y+1,z], M[y,z+1]) antichiain is now in plane partitian format
            if acl != []:
                for i in range(b):
                    i = b-(i+1)
                    for j in range(c):
                        j = c-(j+1)
                        if (ppMatrix[i][j] == 0) and i>=j:
                            iValue = 0
                            jValue = 0
                            if i < b-1:
                                iValue = ppMatrix[i+1][j]
                            if j < c-1:
                                jValue = ppMatrix[i][j+1]
                            ppMatrix[i][j] = max(iValue,jValue)
                        elif j>i:
                            ppMatrix[i][j] = ppMatrix[j][i]

            matrixList.append(ppMatrix) #add PlanePartition to list of plane partitions
        matrixList.sort()
        current = 0
        while current < len(matrixList):
            yield self.element_class(self, matrixList[current])
            current += 1


class PlanePartitions_CSPP(PlanePartitions):

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
        return super(PlanePartitions_CSPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([3,3,3], symmetry=TSPP)
            sage: TestSuite(PP).run()
        """
        if (box_size[0] != box_size[1] or box_size[1] != box_size[2]):
            raise ValueError("invalid box size; must be (r,r,r)")
        super(PlanePartitions_CSPP, self).__init__(box_size=box_size)
        self._box=box_size

    def _repr_(self):
        return " Cyclically Symmetric Plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])


    def __iter__(self):
        def componentwise_comparer2(thing1,thing2):
            x = thing2[0]
            y = thing2[1]
            z = thing2[2]

            if componentwise_comparer(thing1,(x,y,z)) or componentwise_comparer(thing1,(z,x,y)) or componentwise_comparer(thing1,(y,z,x)):
                return True
            return False

        pl = []
        for x in range(0,a):
            for y in range(0, b):
                    for z in range(x,c):
                        if y <= z  and (x != z or y == x):
                            pl.append((x,y,z))

        pocp = Poset((pl,componentwise_comparer2))
        matrixList = [] #list of all PlaneParitions with parameters(a,b,c)
        #iterate through each antichain of product of chains poset with paramaters (a,b,c)
        for acl in pocp.antichains_iterator():
            ppMatrix = [[0] * (c) for i in range(b)] #creates a matrix for the plane parition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]

            #ac format ex: [x,y,z]
            for ac in acl:
                x = ac[0]
                y = ac[1]
                z = ac[2]
                ppMatrix[y][z] = (x+1)
                ppMatrix[z][x] = (y+1)
                ppMatrix[x][y] = (z+1)


            #for each value in current antichain, fill in the rest of the matrix by rule M[y,z] = Max(M[y+1,z], M[y,z+1]) antichiain is now in plane partitian format
            if acl != []:
                for i in range(b):
                    i = b-(i+1)
                    for j in range(c):
                        j = c-(j+1)
                        if (ppMatrix[i][j] == 0):
                            iValue = 0
                            jValue = 0
                            if i < b-1:
                                iValue = ppMatrix[i+1][j]
                            if j < c-1:
                                jValue = ppMatrix[i][j+1]
                            ppMatrix[i][j] = max(iValue,jValue)

            matrixList.append(ppMatrix) #add PlanePartition to list of plane partitions

        matrixList.sort()

        current = 0
        while current < len(matrixList):
            yield self.element_class(self, matrixList[current])
            current += 1

class PlanePartitions_TSSCPP(PlanePartitions):

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
        return super(PlanePartitions_TSSCPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([3,3,3], symmetry=TSPP)
            sage: TestSuite(PP).run()
        """
        if (box_size[0] != box_size[1] or (box_size[1] != box_size[2]) or box_size[0] % 2 != 0):
            raise ValueError("invalid box size; must be (2r,2r,2r)")
        super(PlanePartitions_TSSCPP, self).__init__(box_size=box_size)
        self._box=box_size

    def _repr_(self):
        return " Totally Symmetric Self-Complementary Plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __iter__(self):
        def componentwise_comparer(thing1,thing2):
            if len(thing1) == len(thing2):
                if all(thing1[i] <= thing2[i] for i in range(len(thing1))):
                    return True
            return False
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        n = a
        b = n
        c = n

        pl = []
        for x in range(0,n/2 - 2 + 1):
            for y in range(x, n/2 - 2 + 1):
                    for z in range(0,n/2 - 2 + 1):
                        if z <= n/2 - 2 - y:
                            pl.append((x,y,z))

        pocp = Poset((pl,componentwise_comparer))

        matrixList = [] #list of all PlaneParitions with parameters(a,b,c)
        #iterate through each antichain of product of chains poset with paramaters (a,b,c)
        for acl in pocp.antichains_iterator():
            #ac format ex: [x,y,z]
            ppMatrix = [[0] * (c) for i in range(b)] #creates a matrix for the plane parition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]
            width = n/2 - 1
            height = n/2 - 1

            #generate inner triagle
            for i in range(width):
                for j in range(height):
                    if(i <= j):
                        for ac in acl:
                            if ac[0] == i and ac[1] == j:
                                zVal = ac[2]
                                matrixVal = ppMatrix[j +(n/2)] [i+ (n/2)]
                                if zVal + 1 > matrixVal:
                                    ppMatrix[j +(n/2)] [i+ (n/2)]= zVal + 1

            #fill back
            for i in range(width):
                i = width-(i+1)
                i = i + n/2
                for j in range(height):
                    j = height-(j+1)
                    j = j + n/2
                    if (ppMatrix[i][j] == 0):
                        if i >= j:
                            iValue = 0
                            jValue = 0
                            if i < n:
                                iValue = ppMatrix[i+1][j]
                            if j < n:
                               jValue = ppMatrix[i][j+1]
                            ppMatrix[i][j] = max(iValue,jValue)


            #fill half of triangle symmetrically
            for i in range(width):
                i = i + n/2
                for j in range(height):
                    j = j + n/2
                    if i >= j:
                        ppMatrix[j][i] = ppMatrix[i][j]

            #upper left box
            for i in range(n/2):
                for j in range(n/2):
                    ppMatrix[i][j] = n - ppMatrix[n-(i+1)][n-(j+1)]


            #fill in lower left cube with values n/2
            for i in range(n/2):
                for j in range(n/2):
                    x = i
                    y = j
                    if(ppMatrix[x][y+(n/2)]) == 0:
                        ppMatrix[x][y+(n/2)] = n/2
                    if(ppMatrix[x+(n/2)][y]) == 0:
                        ppMatrix[x+(n/2)][y] = n/2


            #add and subtract values from lower left cube to be rotation of lower right cube
            for i in range(n/2):
                for j in range(n/2):
                    x = i+(n/2)
                    y = j+(n/2)
                    if ppMatrix[x][y] > 0:
                        z = ppMatrix[x][y]
                        for cVal in range(z):
                            #build onto lower left cube
                            ppMatrix[x][0+cVal] += 1
                            #carve out of lower left cube
                            ppMatrix[n-(1+cVal)][(n/2)-(j+1)] -=1

            #fill in upper right cube symmetrically with lower left
            for i in range(n/2):
                for j in range(n/2):
                    ppMatrix[j][i+(n/2)] = ppMatrix[i+(n/2)][j]

            matrixList.append(ppMatrix) #add PlanePartition to list of plane partitions

        matrixList.sort()

        current = 0
        while current < len(matrixList):
            yield self.element_class(self, matrixList[current])
            current += 1
