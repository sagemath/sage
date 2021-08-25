# -*- coding: utf-8 -*-
r"""
Plane Partitions

AUTHORS:

- Jang Soo Kim (2016): Initial implementation
- Jessica Striker (2016): Added additional methods
- Kevin Dilks (2021): Added symmetry classes
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
from typing import NewType, Iterator, Tuple

from sage.structure.list_clone import ClonableList, ClonableArray
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.posets.posets import Poset
from sage.combinat.posets.poset_examples import posets
from sage.categories.cartesian_product import cartesian_product
from sage.rings.integer import Integer
from sage.misc.all import prod
from sage.combinat.tableau import Tableau
from sage.arith.misc import Sigma
from sage.functions.other import floor, ceil, binomial, factorial



PP = NewType('PP', 'PlanePartition')


#class PlanePartition(ClonableArray,
#        metaclass=InheritComparisonClasscallMetaclass):
#    r"""
#    A plane partition.

#    A *plane partition* is a stack of cubes in the positive orthant.

#    INPUT:

#    - ``PP`` -- a list of lists which represents a tableau

#    - ``box_size`` -- (optional) a list ``[A, B, C]`` of 3 positive integers,
#      where ``A``, ``B``, ``C`` are the lengths of the box in the `x`-axis,
#      `y`-axis, `z`-axis, respectively; if this is not given, it is
#      determined by the smallest box bounding ``PP``

#    OUTPUT:

#    The plane partition whose tableau representation is ``PP``.

#    EXAMPLES::

#        sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
#        sage: PP
#        Plane partition [[4, 3, 3, 1], [2, 1, 1], [1, 1]]

#    TESTS::

#        sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
#        sage: TestSuite(PP).run()
#    """

#@add_metaclass(InheritComparisonClasscallMetaclass)
class PlanePartition(ClonableList, metaclass=InheritComparisonClasscallMetaclass):
    @staticmethod
    def __classcall_private__(cls, PP):
        """
        Construct a plane partition with the appropriate parent.

        EXAMPLES::

            sage: p = PlanePartition([[2,1],[1]])
            sage: TestSuite(p).run()

            sage: p.parent()
            Plane partitions
            sage: p.category()
            Category of elements of Plane partitions
            sage: type(p)
            <class 'sage.combinat.plane_partition.PlanePartitions_all_with_category.element_class'>
        """
        if isinstance(PP,PlanePartition):
            return PP
        pp = PlanePartitions()
        return pp.element_class(pp, PP)  # The check() will raise the appropriate error

    def __init__(self, parent, pp, check=True):
        r"""
        Initialize a plane partition.

        TESTS::

            sage: a = PlanePartitions()([[2,1],[1]])
            sage: b = PlanePartitions([2,2,2])([[2,1],[1]])
            sage: c = PlanePartitions(4)([[2,1],[1]])

        Add more tests to show which parent a,b,c receive, check that a==b, and b==c, but a is not b, and b is not c.

        """
        if isinstance(pp, PlanePartition):
            ClonableList.__init__(self, parent, pp, check=False)
        else:
            pp = [list(_) for _ in pp]
            if pp:
                for i in reversed(range(len(pp))):
                    while pp[i] and not pp[i][-1]:
                        del pp[i][-1]
                    if not pp[i]:
                        pp.pop(i)        
            ClonableList.__init__(self, parent, pp, check=check)
        if self.parent()._box is None:
            if pp:
                self._max_x = len(pp)
                self._max_y = len(pp[0])
                self._max_z = pp[0][0]
            else:
                self._max_x = 0
                self._max_y = 0
                self._max_z = 0
        else:
            (self._max_x, self._max_y, self._max_z) = self.parent()._box

    def check(self):
        """
        Check to see that ``self`` is a valid plane partition.

        EXAMPLES::

            sage: a = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: a.check()
            sage: b = PlanePartition([[1,2],[1]])
            Traceback (most recent call last):
            ...          
            ValueError: Not weakly decreasing along rows
            sage: c = PlanePartition([[1,1],[2]])
            Traceback (most recent call last):
            ...
            ValueError: Not weakly decreasing along columns
        """
        for row in self:
            if not all(c >= 0 for c in row):
                raise ValueError("Entries not all nonnegative")
            if not all(row[i] >= row[i+1] for i in range(len(row)-1)):
                raise ValueError("Not weakly decreasing along rows")
        for row, next in zip(self, self[1:]):
            if not all(row[c] >= next[c] for c in range(len(next))):
                raise ValueError("Not weakly decreasing along columns")

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
        return Tableau(self)

    def z_tableau(self, tableau=True) -> Tableau:
        r"""
        Return the projection of ``self`` in the `z` direction.

        If ``tableau`` is set to ``False``, then only the list of lists
        consisting of the projection of boxes size onto the xy-plane
        is returned instead of a Tableau object. This output will
        not have empty trailing rows or trailing zeros removed.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.z_tableau()
            [[4, 3, 3, 1], [2, 1, 1, 0], [1, 1, 0, 0]]
        """
        Z = [[0 for i in range(self._max_y)] for j in range(self._max_x)]
        for C in self.cells():
            Z[C[0]][C[1]] += 1
        if tableau:
            return Tableau(Z)
        return Z

    def y_tableau(self, tableau=True) -> Tableau:
        r"""
        Return the projection of ``self`` in the `y` direction.

        If ``tableau`` is set to ``False``, then only the list of lists
        consisting of the projection of boxes size onto the xz-plane
        is returned instead of a Tableau object. This output will
        not have empty trailing rows or trailing zeros removed.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.y_tableau()
            [[4, 3, 2], [3, 1, 0], [3, 0, 0], [1, 0, 0]]
        """
        Y = [[0 for i in range(self._max_x)] for j in range(self._max_z)]
        for C in self.cells():
            Y[C[2]][C[0]] += 1
        if tableau:
            return Tableau(Y)
        return Y

    def x_tableau(self, tableau=True) -> Tableau:
        r"""
        Return the projection of ``self`` in the `x` direction.

        If ``tableau`` is set to ``False``, then only the list of lists
        consisting of the projection of boxes size onto the yz-plane
        is returned instead of a Tableau object. This output will
        not have empty trailing rows or trailing zeros removed.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.x_tableau()
            [[3, 2, 1, 1], [3, 1, 1, 0], [2, 1, 1, 0], [1, 0, 0, 0]]
        """
        X = [[0 for i in range(self._max_z)] for j in range(self._max_y)]
        for C in self.cells():
            X[C[1]][C[2]] += 1
        if tableau:
            return Tableau(X)
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

    def number_of_boxes(self) -> Integer:
        r"""
        Return the number of boxes in the plane partition.

        EXAMPLES::

            sage: PP = PlanePartition([[3,1],[2]])
            sage: PP.number_of_boxes()
            6
        """
        return sum(sum(k for k in row) for row in self)

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

    def plot(self, show_box=False, colors=["white", "lightgray", "darkgray"]):
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
        
    def contains(self, PP) -> bool: 
        """
        
        Return ``True`` if ``PP`` is a plane partition that fits
        inside ``self``.
        
        Specifically, ``self`` contains ``PP`` if, for all `i`, `j`,
        the height of ``PP`` at `ij` is less than or equal to the
        height of ``self`` at `ij`.
        
        EXAMPLES::
        
            sage: P1 = PlanePartition([[5,4,3],[3,2,2],[1]])
            sage: P2 = PlanePartition([[3,2],[1,1],[0,0],[0,0]])
            sage: P3 = PlanePartition([[5,5,5],[2,1,0]])
            sage: P1.contains(P2)
            True
            sage: P2.contains(P1)
            False
            sage: P1.contains(P3)
            False
            sage: P3.contains(P2)
            True
        """
        PP = PlanePartition(PP)
        if len(self) < len(PP):
            return False
        
        for i in range(len(PP)):
            if len(self[i]) < len(PP[i]):
                return False
        
        return all([self[i][j] >= PP[i][j] for j in range(len(PP[i])) for i in range(len(PP))])        
        

    def complement(self, tableau_only=False) -> PP:
        # Complement needs to be more intelligent about which parent to return
        r"""
        Return the complement of ``self``.

        If the parent of ``self`` consists only of partitions inside a given
        box, then the complement is taken in this box. Otherwise, the
        complement is taken in the smallest box containing the plane partition.
        The empty plane partition with no box specified is its own complement.

        If ``tableau_only`` is set to ``True``, then only the tableau
        consisting of the projection of boxes size onto the xy-plane
        is returned instead of a PlanePartition object. This output will
        not have empty trailing rows or trailing zeros removed.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.complement()
            Plane partition [[4, 4, 3, 3], [4, 3, 3, 2], [3, 1, 1]]
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
        elif not self.parent()._box:
            pp = PlanePartitions()
            return pp.element_class(pp, T)
        else:            
            return type(self)(self.parent(), T, check=False)

    def transpose(self, tableau_only=False) -> PP:
        r"""
        Return the transpose of ``self``.

        If ``tableau_only`` is set to ``True``, then only the tableau
        consisting of the projection of boxes size onto the xy-plane
        is returned instead of a PlanePartition object. This will
        not necessarily have trailing rows or trailing zeros removed.


        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.transpose()
            Plane partition [[4, 2, 1], [3, 1, 1], [3, 1], [1]]
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
        elif self.parent()._box == None or self.parent()._box[0] == self.parent()._box[1]:
            return type(self)(self.parent(), T, check=False)
        new_box = (self.parent()._box[1],self.parent()._box[0],self.parent()._box[2])
        return PlanePartitions(new_box,symmetry=self.parent._symmetry)
 #       elif self._max_x != self._max_y:
 #           raise ValueError("Tranpose only supports parents with symmetric dimensions")           
 #       else:
 #           return type(self)(self.parent(), T, check=False)

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
            False
            sage: PP = PlanePartitions([4,4,4])([[4,4,4,4],[4,4,2,0],[4,2,0,0],[0,0,0,0]])
            sage: PP.is_SCPP()
            True
        """
        return self.z_tableau() == self.complement(tableau_only=True)

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

    def to_order_ideal(self): 
        """
        Return the order ideal corresponding to ``self``.

        TODO: As many families of symmetric plane partitions are in bijection
              with order ideals in an associated poset, this function could
              feasibly have options to send symmetric plane partitions
              to the associated order ideal in that poset, instead.
        
        EXAMPLES::
        
            sage: PlanePartition([[3,2,1],[2,2],[2]]).to_order_ideal()
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 0), (0, 1, 1), (0, 2, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1), (2, 0, 0), (2, 0, 1)]
            sage: PlanePartition([[2,1],[1],[1]]).to_order_ideal()
            [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (2, 0, 0)]
        """
        (a, b, c) = (self._max_x, self._max_y, self._max_z)
        Q = posets.ProductOfChains([a,b,c])
        count = 0
        generate = []
        for i in range(len(self)):
            for j in range(len(self[i])):
                if (self[i][j] > 0):
                    generate.append((i,j,self[i][j]-1))
            count += 1
        oi = Q.order_ideal(generate)
        return oi
        


    def maximal_boxes(self) -> list:
        """
        Return the coordinates of the maximal boxes of ``self``.
        
        EXAMPLES::
        
            sage: PlanePartition([[3,2,1],[2,2],[2]]).maximal_boxes()
            [[0, 2, 0], [1, 1, 1], [0, 0, 2], [2, 0, 1]]
            sage: PlanePartition([[2,1],[1],[1]]).maximal_boxes()
            [[0, 1, 0], [0, 0, 1], [2, 0, 0]]
        """
        (a, b, c) = (self._max_x, self._max_y, self._max_z)
        Q = posets.ProductOfChains([a,b,c])
        count = 0
        generate = []
        for i in range(len(self)):
            for j in range(len(self[i])):
                if (self[i][j] > 0):
                    generate.append((i,j,self[i][j]-1))
            count+=1
        oi = Q.order_ideal_generators(generate)
        return [list(oi_elem) for oi_elem in oi]
  
    def cyclically_rotate(self, preserve_parent=False) -> PP:
        """
        Return the cyclic rotation of ``self``.

        By default, if the parent of ``self`` consists of plane
        partitions inside an `a \times b \times c` box, the result 
        will have a parent consisting of partitions inside
        a `c \times a \times b` box, unless the optional parameter
        ``preserve_parents`` is set to ``True``. Enabling this setting
        may give an element that is NOT an element of its parent.
        
        EXAMPLES::
        
            sage: PlanePartition([[3,2,1],[2,2],[2]]).cyclically_rotate()
            Plane partition [[3, 3, 1], [2, 2], [1]]
            sage: PP = PlanePartition([[4,1],[1],[1]])
            sage: PP.cyclically_rotate()
            Plane partition [[3, 1, 1, 1], [1]]
            sage: PP == PP.cyclically_rotate().cyclically_rotate().cyclically_rotate()
            True
            sage: PP = PlanePartitions([4,3,2]).random_element()
            sage: PP.cyclically_rotate().parent()
            Plane partitions inside a 2 x 4 x 3 box
            sage: PP = PlanePartitions([3,4,2])([[2,2,2,2],[2,2,2,2],[2,2,2,2]])
            sage: PP_rotated = PP.cyclically_rotate(preserve_parent=True)
            sage: PP_rotated in PP_rotated.parent()
            False
        """
        (a, b, c) = (self._max_x, self._max_y, self._max_z)
        new_antichain = []
        for elem in self.maximal_boxes():
            new = (elem[1], elem[2], elem[0])
            new_antichain.append(new)
        ppMatrix = [[0] * (c) for i in range(b)]
        for box in new_antichain:
            y = box[0]
            z = box[1]
            x = box[2]
            ppMatrix[y][z] = (x+1)
        if new_antichain != []:
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
        # Start code for determining correct parent
        if self.parent()._box == None or preserve_parent == True or (self.parent()._box[0] == self.parent()._box[1] == self.parent()._box[2]):
            return type(self)(self.parent(), ppMatrix, check=False)
        new_box = (self.parent()._box[2],self.parent()._box[0],self.parent()._box[1])
        return PlanePartitions(new_box,symmetry=self.parent()._symmetry)(ppMatrix)
        #return PlanePartition(ppMatrix)

class PlanePartitions(UniqueRepresentation, Parent):
    r"""
    A factory class for plane partitions.

    PlanePartitions([a,b,c]) returns the class of plane partitions that fit
    inside an a \times b \times c box.

    Optional keyword is 'symmetry'.

    Describe options.

    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwds):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`PlanePartitions`
        for more information.

        TESTS::

            sage: PlanePartitions()
            Plane partitions
            sage: PlanePartitions([3,3,3])
            Plane partitions inside a 3 x 3 x 3 box
            sage: PlanePartitions(3)
            Plane partitions of size 3
            sage: PlanePartitions([4,4,4], symmetry='TSSCPP')
            Totally symmetric self-complementary plane partitions inside a 4 x 4 x 4 box
        """
        symmetry = kwds.get('symmetry', None)
        if not args:
            return PlanePartitions_all()
        else:
            box_size = None
            if args:
                # The first arg could be either a size or a box size
                if isinstance(args[0], (int, Integer)):
                    return PlanePartitions_n(args[0])
                else:
                    box_size = args[0]
            if symmetry == None:
                return PlanePartitions_box(box_size)
            elif symmetry == 'SPP':
                return PlanePartitions_SPP(box_size)
            elif symmetry == 'CSPP':
                return PlanePartitions_CSPP(box_size)
            elif symmetry == 'TSPP':
                return PlanePartitions_TSPP(box_size)
            elif symmetry == 'SCPP':
                return PlanePartitions_SCPP(box_size)
            elif symmetry == 'TCPP':
                return PlanePartitions_TCPP(box_size)
            elif symmetry == 'SSCPP':
                return PlanePartitions_SSCPP(box_size)
            elif symmetry == 'CSTCPP':
                return PlanePartitions_CSTCPP(box_size)
            elif symmetry == 'CSSCPP':
                return PlanePartitions_CSSCPP(box_size)
            elif symmetry == 'TSSCPP':
                return PlanePartitions_TSSCPP(box_size)
            else:
                raise ValueError("invalid symmetry class option; must be None, 'SPP', 'CSPP', 'TSPP', 'SCPP', 'TCPP', 'SSCPP', 'CSTCPP', 'CSSCPP', or 'TSSCPP' ")


    Element = PlanePartition


    def __contains__(self, pp):
        """
        Check to see that ``self`` is a valid plane partition.

        .. TODO:
            
            Figure out how redundant this is, given that the check function
            exists for the factor class. Maybe only need __contains__
            on the fixed size and symmetry classes?
        """
        for row in pp:
            if not all(c >= 0 for c in row):
                return False
            if not all(row[i] >= row[i+1] for i in range(len(row)-1)):
                return False
        for row, next in zip(pp, pp[1:]):
            if not all(row[c] >= next[c] for c in range(len(next))):
                return False
        return True





class PlanePartitions_all(PlanePartitions):
    r"""
    All plane partitions.

    .. TODO:

    Consider giving this the structure of disjoint union of the classes
    PlanePartitions(n) for n an integer.
    """


    def __init__(self):
        r"""
        Initializes the class of all plane partitions.

        .. WARNING::

            Input is not checked; please use :class:`PlanePartitions` to
            ensure the options are properly parsed.

        TESTS::

            sage: from sage.combinat.plane_partition import PlanePartitions_all
            sage: P = PlanePartitions_all()
            sage: TestSuite(P).run()  # long time
        """
        self._box = None
        self._symmetry = None
        super(PlanePartitions_all, self).__init__(category=InfiniteEnumeratedSets())

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PlanePartitions()
            Plane partitions
        """
        return "Plane partitions"




class PlanePartitions_box(PlanePartitions):
    r"""
    All plane partitions that fit inside a box of a specified size.
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
        return super(PlanePartitions_box, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        r"""
        Initializes the class of plane partitions that fit in a box of a 
        specified size.

        EXAMPLES::

            sage: PP = PlanePartitions((4,3,2))
            sage: TestSuite(PP).run()
        """
        super(PlanePartitions_box,self).__init__(category=FiniteEnumeratedSets())
        self._box = box_size
        self._symmetry = None

    def _repr_(self) -> str:

        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PlanePartitions((4,3,2))
            Plane partitions inside a 4 x 3 x 2 box
        """

        return "Plane partitions inside a {} x {} x {} box".format(
            self._box[0], self._box[1], self._box[2])

    def __contains__(self, x):
        if len(x) == 0:
            return True
        return PlanePartitions.__contains__(self, x) and len(x) <= self._box[0] and len(x[0]) <= self._box[1] and x[0][0] <= self._box[2]

    def to_poset(self):
        r"""
        Returns the product of three chains poset, whose order ideals are
        naturally in bijection with plane partitions inside a box.
        """
        a=self._box[0]
        b=self._box[1]    
        c=self._box[2]
        return posets.ProductOfChains([a,b,c])

    def from_order_ideal(self, I) -> PP:
        r"""
        Return the plane partition corresponding to an order ideal in the
        poset given in :meth:`to_poset`.

        Note: input may not be checked ? Optional check parameter if too much overhead?
        """
        return self.from_antichain(self.to_poset().order_ideal_generators(I))

    def from_antichain(self, A) -> PP:
        r"""
        Return the plane partition corresponding to an antichain in the poset
        given in :meth:`to_poset`.

        Note: input may not be checked? Optional parameter if too much overhead?
        """
        a = self._box[0]
        b = self._box[1]
        c = self._box[2]
        ppMatrix = [[0] * (b) for i in range(a)] #creates a matrix for the plane partition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]

        #ac format ex: [x,y,z]
        #iterate through each antichain, assigning the y,z position in ppMatrix = the height of the stack (x + 1)
        for ac in A:
            x = ac[0]
            y = ac[1]
            z = ac[2]
            ppMatrix[x][y] = (z+1)

        #for each value in current antichain, fill in the rest of the matrix by rule M[y,z] = Max(M[y+1,z], M[y,z+1]) antichiain is now in plane partitian format
        if A != []:
            for i in range(a):
                i = a-(i+1)
                for j in range(b):
                    j = b-(j+1)
                    if (ppMatrix[i][j] == 0):
                        iValue = 0
                        jValue = 0
                        if i < a-1:
                            iValue = ppMatrix[i+1][j]
                        if j < b-1:
                            jValue = ppMatrix[i][j+1]
                        ppMatrix[i][j] = max(iValue,jValue)
        return self.element_class(self, ppMatrix)            
        

    def __iter__(self) -> Iterator:
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((1,2,1)))
            [Plane partition [], Plane partition [[1, 1]], Plane partition [[1]]]
        """

        A = self._box[0]
        B = self._box[1]
        C = self._box[2]
        from sage.combinat.tableau import SemistandardTableaux
        for T in SemistandardTableaux([B for i in range(A)], max_entry=C + A):
            PP = [[0 for i in range(B)] for j in range(A)]
            for r in range(A):
                for c in range(B):
                    PP[A - 1 - r][B - 1 - c] = T[r][c] - r - 1
            yield self.element_class(self, PP, check=False)

#        a = self._box[0]
#        b = self._box[1]
#        c = self._box[2]
#        pocp = posets.ProductOfChains([a,b,c])
#        matrixList = [] #list of all PlanePartitions with parameters(a,b,c)
        #iterate through each antichain of product of chains poset with parameters (a,b,c)
#        for acl in self.to_poset().antichains_iterator():
#            yield self.from_antichain(acl)



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

    def box(self) -> Tuple:
        """
        Return the size of the box of the plane partition of ``self``
        is contained in.

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
            sage: P.random_element()
            Plane partition [[4, 3, 3], [4], [2]]
        """
        Z = self.from_order_ideal(self.to_poset().random_order_ideal())
        return self.element_class(self, Z, check=False)



class PlanePartitions_n(PlanePartitions):
    """
    Plane partitions with a fixed number of boxes.
    """

    def __init__(self, n):
        r"""
        Initializes the class of plane partitions with ``n`` boxes.

        .. WARNING::

            Input is not checked; please use :class:`PlanePartitions` to
            ensure the options are properly parsed. (CHECK: does this make sense?)

        TESTS::

            sage: PP = PlanePartitions(4)
            sage: type(PP)
            <class 'sage.combinat.plane_partition.PlanePartitions_n_with_category'>
            sage: TestSuite(PP).run()
        """
        super(PlanePartitions_n, self).__init__(category=FiniteEnumeratedSets())
        self._n = n
#        self._box = (0,0,0)
        self._box = None
        self._symmetry = None

    def _repr_(self) -> str:
        """
        TESTS::

            sage: PlanePartitions(3)
            Plane partitions of size 3
        """
        return "Plane partitions of size {}".format(self._n)

    def __contains__(self, x) -> bool:
        return PlanePartitions.__contains__(self, x) and x.number_of_boxes() ==  self._n

    def __iter__(self) -> Iterator:
        r"""
        An iterator to generate all plane partitions of a fixed size.

        EXAMPLES::

            sage: list(PlanePartitions(2))
            [Plane partition [[2]], Plane partition [[1, 1]], Plane partition [[1], [1]]]
        """
        from sage.combinat.partition import Partitions
        def PP_first_row_iter(n, la):
            m = n-sum(la)
            if m < 0:
                yield
                return
            if m==0:
                yield [la]
                return
            for k in range(m,0,-1):
                for mu in P_in_shape_iter(k,la):
                    if mu is not None:
                        for PP in PP_first_row_iter(m, mu):
                            if PP is not None:
                                yield [la] + PP


        def P_in_shape_iter(n, la):
            if n<0 or sum(la)<n:
                yield
                return
            if n==0:
                yield []
                return
            if len(la)==1:
                if la[0]>=n:
                    yield [n]
                    return
                else:
                    yield
                    return
            if sum(la)==n:
                yield la
                return
            for mu_0 in range(min(n,la[0]),0,-1):
                new_la = [min(mu_0,la[i]) for i in range(1,len(la))]
                for mu in P_in_shape_iter(n-mu_0, new_la):
                    if mu is not None:
                        yield [mu_0]+mu
        n = self._n
        if n==0:
            yield PlanePartition([])
            return

        for m in range(n,0,-1):
            for la in Partitions(m):
                for a in PP_first_row_iter(n,la):
                    yield self.element_class(self, a, check=False)
            
    def cardinality(self) -> Integer:
        r"""
        Return the number of plane partitions with ``n`` boxes.

        Calculated using the recurrence relation

        .. MATH:

        PL(n) = \sum_{k=1}^n PL(n-k)\sigma_2(k)

        where ``\sigma_k(n)`` is the sum of the kth powers of
        divisors of n.

        """
        PPn = [1]
        for i in range(1,1+self._n):
            nextPPn = sum(PPn[i-k]*Sigma()(k,2) for k in range(1,i+1))/i
            PPn.append(nextPPn)
        return(Integer(PPn[-1]))

#Symmetry classes are enumerated and labelled in order as in Proofs and
#Confirmations/Stanley (with all plane partitions being the first class)




#Class 2

class PlanePartitions_SPP(PlanePartitions):

#How to handle when a!=b ?

    @staticmethod
    def __classcall_private__(cls, box_size):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: P1 = PlanePartitions((3,3,2), symmetry='SPP')
            sage: P2 = PlanePartitions([3,3,2], symmetry='SPP')
            sage: P1 is P2
            True
        """
        return super(PlanePartitions_SPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([3,3,2], symmetry='SPP')
            sage: TestSuite(PP).run()
        """
        if box_size[0] != box_size[1]:
            raise ValueError("x and y dimensions must match")
        super(PlanePartitions_SPP, self).__init__(category=FiniteEnumeratedSets())
        self._box = box_size
        self._symmetry = 'SPP'

    def _repr_(self) -> str:
        return "Symmetric plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __contains__(self, x) -> bool:
        return PlanePartitions.__contains__(self, x) and x.is_SPP()


    def to_poset(self):
        r"""
        Returns a poset whose order ideals are in bijection with
        symmetric plane partitions.
        """
        a=self._box[0]
        b=self._box[1]    
        c=self._box[2]
        pl = []
        cmp = lambda x,y : all(x[i]<=y[i] for i in range(len(x)))
        for x in range(0,a):
            for y in range(0,x+1):
                for z in range(0,c):
                    pl.append((x,y,z))
        return Poset((pl,cmp))

    def from_order_ideal(self, I) -> PP:
        r"""
        Return the plane partition corresponding to an order ideal in the
        poset given in :meth:`to_poset()`.

        Note: input may not be checked ? Optional check parameter if too much overhead?
        """
        return self.from_antichain(self.to_poset().order_ideal_generators(I))

    def from_antichain(self, A) -> PP:
        r"""
        Return the plane partition corresponding to an antichain in the poset
        given in :meth:`to_poset()`.

        Note: input may not be checked? Optional parameter if too much overhead?
        """
        #Initialize an empty plane partition
        a=self._box[0]
        b=self._box[1]    
        c=self._box[2]
        ppMatrix = [[0] * (b) for i in range(a)]
        #Antichain indicates where the 'corners' will be in the
        #plane partition
        for ac in A:
            x=ac[0]
            y=ac[1]
            z=ac[2]
            ppMatrix[x][y] = (z+1)
        #Fill out the rest of the plane partition using symmetry and the
        #rule pp[i][j]=max(pp[i][j+1],pp[i+1][j])
        if A!= []:
            for i in range(a):
                i = a-(i+1)
                for j in range(b):
                    j = b-(j+1)
                    if (ppMatrix[i][j] == 0) and i>=j:
                        iValue = 0
                        jValue = 0
                        if i < a-1:
                            iValue = ppMatrix[i+1][j]
                        if j < b-1:
                            jValue = ppMatrix[i][j+1]
                        ppMatrix[i][j] = max(iValue,jValue)
                    elif j>i:
                        ppMatrix[i][j] = ppMatrix[j][i]
        return self.element_class(self, ppMatrix)

    def __iter__(self) -> Iterator:
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((2,2,1), symmetry='SPP'))
            [Plane partition [],
            Plane partition [[1, 1], [1, 1]],
            Plane partition [[1, 1], [1]],
            Plane partition [[1]]]
        """
        for acl in self.to_poset().antichains_iterator():
            yield self.from_antichain(acl)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of symmetric plane partitions inside an `r \times r \times t`
        box is equal to

        .. MATH::

            \prod_{i=1}^{r} \frac{2*i + t - 1}{2*i - 1}  *
            \prod_{1 \leq i \leq j \leq r} \frac{i+j+t-1}{i+j-1}

        EXAMPLES::

            sage: P = PlanePartitions((3,3,2), symmetry='SPP')
            sage: P.cardinality()
            35
        """
        A = self._box[0]
        B = self._box[1]
        C = self._box[2]
        leftProduct = (prod( (2*i + C - 1) / (2*i - 1) for i in range(1,A+1)))
        rightProduct = (prod( (i + j + C - 1) / (i + j - 1) for j in range(1, A+1) for i in range(1, j) ))
        return Integer(leftProduct * rightProduct)

    def random_element(self) -> PP:
        r"""
        Return a uniformly random element of ``self``.

        ALGORITHM:

        This uses the
        :meth:`~sage.combinat.posets.posets.FinitePoset.random_order_ideal`
        method and the natural bijection between symmetric plane partitions
        and antichains in an associated poset. (FIX EXAMPLES)

        EXAMPLES::

            sage: P = PlanePartitions((4,3,5))
            sage: P.random_element()
            Plane partition [[4, 3, 3], [4], [2]]
        """
        Z = self.from_order_ideal(self.to_poset().random_order_ideal())
        return self.element_class(self, Z, check=False)

#Class 3

class PlanePartitions_CSPP(PlanePartitions):

#For some reason works if not(a == b == c)

    @staticmethod
    def __classcall_private__(cls, box_size):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: P1 = PlanePartitions((4,4,4), symmetry='CSPP')
            sage: P2 = PlanePartitions([4,4,4], symmetry='CSPP')
            sage: P1 is P2
            True
        """
        return super(PlanePartitions_CSPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([3,3,3], symmetry='CSPP')
            sage: TestSuite(PP).run()
        """
        if box_size[0] != box_size[1] or box_size[1] != box_size[2]:
            raise ValueError("x, y, and z dimensions must match")
        super(PlanePartitions_CSPP, self).__init__(category=FiniteEnumeratedSets())
        self._box = box_size
        self._symmetry = 'CSPP'

    def _repr_(self) -> str:
        return "Cyclically symmetric plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def to_poset(self):
        """
        MAKE SURE DOC EXPLICITLY DESCRIBES WHAT FUNDAMENTAL DOMAIN LOOKS LIKE
        """
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        cmp = lambda x,y : all(x[i]<= y[i] for i in range(len(x)))
        cmp2 = lambda x,y : cmp(x,y) or cmp(x,(y[2],y[0],y[1])) or cmp(x,(y[1],y[2],y[0]))
        pl = []
        for x in range(0,a):
            for y in range(0, b):
                    for z in range(x,c):
                        if y <= z  and (x != z or y == x):
                            pl.append((x,y,z))
        return Poset((pl, cmp2))

    def from_antichain(self, acl):
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        ppMatrix = [[0] * (c) for i in range(b)] 
        #creates a matrix for the plane parition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]
        #ac format ex: [x,y,z]
        for ac in acl:
            x = ac[0]
            y = ac[1]
            z = ac[2]
            ppMatrix[y][z] = (x+1)
            ppMatrix[z][x] = (y+1)
            ppMatrix[x][y] = (z+1)


        #for each value in current antichain, fill in the rest of the 
        #matrix by rule M[y,z] = Max(M[y+1,z], M[y,z+1]) antichiain is 
        #now in plane partition format
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
        return self.element_class(self, ppMatrix)

    def from_order_ideal(self, I) -> PP:
        r"""
        Return the plane partition corresponding to an order ideal in the
        poset given in :meth:`to_poset`.

        Note: input may not be checked ? Optional check parameter if too much overhead?
        """
        return self.from_antichain(self.to_poset().order_ideal_generators(I))    

    def __iter__(self) -> Iterator:
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((2,2,2), symmetry='CSPP'))
            [Plane partition [],
            Plane partition [[2, 2], [2, 2]],
            Plane partition [[2, 2], [2, 1]],
            Plane partition [[2, 1], [1]],
            Plane partition [[1]]]
        """
        for acl in self.to_poset().antichains_iterator():
            yield self.from_antichain(acl)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of cyclically symmetric plane partitions inside an 
        `r \times r \times r` box is equal to

        .. MATH::

            \left(\prod_{i=1}^{r} \frac{3i - 1}{3i - 2}\right)  
            \left(\prod_{1 \leq i \leq j \leq r} \frac{i+j+r-1}{2i+j-1}\right)

        EXAMPLES::

            sage: P = PlanePartitions((4,4,4), symmetry='CSPP')
            sage: P.cardinality()
            132
        """
        A = self._box[0]
        B = self._box[1]
        C = self._box[2]
        numerator = prod(3*i-1 for i in range(1, A+1)) * prod( (i+j+A-1) for j in range(1,A+1) for i in range(1,j+1))
        denominator = prod(3*i-2 for i in range(1, A+1)) * prod( (2*i+j-1) for j in range(1,A+1) for i in range(1,j+1))
        return Integer(numerator/denominator)


# Class 4


class PlanePartitions_TSPP(PlanePartitions):
# Totally symmetric plane partitions

# Make sure inputs checked, code doesn't have a,b,c treated properly

    @staticmethod
    def __classcall_private__(cls, box_size):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: P1 = PlanePartitions((4,4,4), symmetry='TSPP')
            sage: P2 = PlanePartitions([4,4,4], symmetry='TSPP')
            sage: P1 is P2
            True
        """
        return super(PlanePartitions_TSPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([3,3,3], symmetry='TSPP')
            sage: TestSuite(PP).run()
        """
        if box_size[0] != box_size[1] or box_size[1] != box_size[2]:
            raise ValueError("x, y, and z dimensions must match")
        super(PlanePartitions_TSPP, self).__init__(category=FiniteEnumeratedSets())
        self._box = box_size
        self._symmetry = 'TSPP'

    def _repr_(self) -> str:
        return "Totally symmetric plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def to_poset(self):
        """
        MAKE SURE DOC EXPLICITLY DESCRIBES WHAT FUNDAMENTAL DOMAIN LOOKS LIKE
        """
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        cmp = lambda x,y : all(x[i]<= y[i] for i in range(len(x)))
        pl = []
        for x in range(0,a):
            for y in range(x, b):
                    for z in range(y,c):
                        pl.append((x,y,z))
        return Poset((pl,cmp))

    def from_antichain(self, acl) -> PP:
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        ppMatrix = [[0] * (c) for i in range(b)] #creates a matrix for the plane partition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]
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
        return self.element_class(self, ppMatrix)

    def random_element(self) -> PP:
        return self.to_poset().from_antichain(self.to_poset().random_antichain())

    def __iter__(self) -> Iterator:
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((4,3,2), symmetry='TSPP'))
            [Plane partition [],
            Plane partition [[2, 2], [2, 2]],
            Plane partition [[2, 2], [2, 1]],
            Plane partition [[2, 1], [1]],
            Plane partition [[1]]]
        """
        for A in self.to_poset().antichains_iterator():
            yield self.from_antichain(A)


    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of totally symmetric plane partitions inside an 
        `r \times r \times r` box is equal to

        .. MATH::

            \prod_{1 \leq i \leq j \leq r} \frac{i+j+r-1}{i+2j-2}

        EXAMPLES::

            sage: P = PlanePartitions((4,4,4), symmetry='TSPP')
            sage: P.cardinality()
            66
        """
        A = self._box[0]
        B = self._box[1]
        C = self._box[2]
        return Integer(prod((i + j + A - 1) / (i + 2*j - 2) for j in range(1,A+1) for i in range(1,j+1) ))



# Class 5

class PlanePartitions_SCPP(PlanePartitions):
# Self-complementary plane partitions
    @staticmethod
    def __classcall_private__(cls, box_size):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: P1 = PlanePartitions((4,3,2), symmetry='SCPP')
            sage: P2 = PlanePartitions([4,3,2], symmetry='SCPP')
            sage: P1 is P2
            True
        """
        return super(PlanePartitions_SCPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([4,3,2], symmetry='SCPP')
            sage: TestSuite(PP).run()
        """
        if (box_size[0] % 2 == 1 and box_size[1] % 2 == 1 and box_size[2] % 2 == 1):
            raise ValueError("box sides cannot all be odd")
        super(PlanePartitions_SCPP, self).__init__(category=FiniteEnumeratedSets())
        self._box = box_size
        self._symmetry = 'SCPP'

    def _repr_(self) -> str:
        return "Self-complementary plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __iter__(self) -> Iterator:
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((2,2,2), symmetry='SCPP'))
            [Plane partition [],
            Plane partition [[2, 2], [2, 2]],
            Plane partition [[2, 2], [2, 1]],
            Plane partition [[2, 1], [1]],
            Plane partition [[1]]]
        """
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        def Partitions_inside_lambda(la):
            "Returns the list of partitions contained in la with the same number of parts including 0s."
            if len(la)==0:
                yield []
                return
            for mu_0 in range(la[0],0,-1):
                new_la = [min(mu_0,la[i]) for i in range(1,len(la))]
                for mu in Partitions_inside_lambda(new_la):
                    yield [mu_0]+mu
            yield [0 for i in la]
            return

        def Partitions_inside_lambda_with_smallest_at_least_k(la,k):
            "Returns the list of partitions contained in la with the smallest entry at least k"
            if len(la)==0:
                yield []
                return
            if la[-1] < k:
                yield
                return

            for mu in Partitions_inside_lambda([la[i]-k for i in range(len(la))]):
                yield ([mu[i]+k for i in range(len(la))])
            return

        def possible_middle_row_for_b_odd(a,c):
            "Returns the list of possible middle row for SCPP inside box(a,b,c) when b is odd"
            if a*c % 2 == 1:
                yield
                return
            for mu in Partitions_inside_lambda([floor(c/2) for i in range(floor(a/2))]):
                nu = [c-mu[len(mu)-1-i] for i in range(len(mu))]
                if a % 2 ==0:
                    la = nu + mu
                else:
                    la = nu + [c/2] + mu
                yield (la)
            return

        def possible_middle_row_for_b_even(a,c):
            "Returns the list of possible middle ((b/2)+1)st row for SCPP inside box(a,b,c) when b is even"
#            for mu in Partitions_inside_lambda([floor(c/2) for i in range((a+1)/2)]):
            for mu in Partitions_inside_lambda([floor(c/2) for i in range(floor((a+1)/2))]):
                nu = [c-mu[len(mu)-1-i] for i in range(floor(a/2))]
                for tau in Partitions_inside_lambda_with_smallest_at_least_k(nu,mu[0]):
                    la = tau + mu
                    yield (la)
            return



        def PPs_with_first_row_la_and_with_k_rows(la,k):
            "Returns PPs with first row la and with k rows in total"
            if k == 0:
                yield []
                return
            if k == 1:
                yield [la]
                return
            for mu in Partitions_inside_lambda(la):
                for PP in PPs_with_first_row_la_and_with_k_rows(mu,k-1):
                    yield ([la]+PP)

            return

        def complement(PP,c):
            "Returns the complement of PP with respect to height c"
            if len(PP) == 0:
                return []
            b = len(PP)
            a = len(PP[0])
            return [[c-PP[b-1-i][a-1-j] for j in range(a)] for i in range(b)]

        if a*b*c % 2 == 1:
            return

        if b % 2 == 1:
            for la in possible_middle_row_for_b_odd(a,c):   # la is the middle row of SCPP
                for PP in PPs_with_first_row_la_and_with_k_rows(la,(b+1)/2):
                    PP_below = PP[1:]
                    PP_above = complement(PP_below,c)
                    yield self.element_class(self, PP_above+[la]+PP_below)
        else:
            for la in possible_middle_row_for_b_even(a,c):   # la is the middle ((a/2)+1)st row of SCPP
                for PP in PPs_with_first_row_la_and_with_k_rows(la,b/2):
                    PP_below = PP
                    PP_above = complement(PP_below,c)
                    yield self.element_class(self, PP_above+PP_below)
        return

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of self complementary plane partitions inside a 
        `2a \times 2b \times 2c` box is equal to

        .. MATH::

            \left(\prod_{i=1}^{r}\prod_{j=1}^{b} \frac{i + j + c - 1}{i + j - 1}\right)^2



        The number of self complementary plane partitions inside an 
        `(2a+1) \times 2b \times 2c` box is equal to

        .. MATH::

            \left(\prod_{i=1}^{a}\prod_{j=1}^{b} \frac{i+j+c-1}{i+j-1} \right)
            \left(\prod_{i=1}^{a+1}\prod_{j=1}^{b} \frac{i+j+c-1}{i+j-1} \right)


        The number of self complementary plane partitions inside an 
        `(2a+1) \times (2b+1) \times 2c` box is equal to

        .. MATH::

            \left(\prod_{i=1}^{a+1}\prod_{j=1}^{b} \frac{i+j+c-1}{i+j-1} \right)
            \left(\prod_{i=1}^{a}\prod_{j=1}^{b+1} \frac{i+j+c-1}{i+j-1} \right)

        EXAMPLES::

            sage: P = PlanePartitions((4,4,4), symmetry='SCPP')
            sage: P.cardinality()
            400

            sage: P = PlanePartitions((5,4,4), symmetry='SCPP')
            sage: P.cardinality()
            1000

            sage: P = PlanePartitions((5,5,4), symmetry='SCPP')
            sage: P.cardinality()
            2500
        """
        r=self._box[0]
        s=self._box[1]
        t=self._box[2]
        if r % 2 == 0 and s % 2 == 0 and t % 2 == 0:
            return Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k -2 ) for i in range(1,1+r/2) for j in range(1,1+s/2) for k in range(1,1+t/2)))) * Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k -2 ) for i in range(1,1+r/2) for j in range(1,1+s/2) for k in range(1,1+t/2))))
        if r % 2 == 1 and s % 2 == 0 and t % 2 == 0:
            return Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+(r-1)/2) for j in range(1,1+s/2) for k in range(1,1+t/2)))) * Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+(r-1)/2+1) for j in range(1,1+s/2) for k in range(1,1+t/2))))
        if r % 2 == 0 and s % 2 == 1 and t % 2 == 0:
            return Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+r/2) for j in range(1,1+(s-1)/2) for k in range(1,1+t/2)))) * Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+r/2) for j in range(1,1+(s-1)/2+1) for k in range(1,1+t/2))))
        if r % 2 == 0 and s % 2 == 0 and t % 2 == 1:
            return Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+r/2) for j in range(1,1+s/2) for k in range(1,1+(t-1)/2)))) * Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+r/2) for j in range(1,1+s/2) for k in range(1,1+(t-1)/2+1))))
        if r % 2 == 1 and s % 2 == 1 and t % 2 == 0:
            return Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+(r-1)/2+1) for j in range(1,1+(s-1)/2) for k in range(1,1+t/2)))) * Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+(r-1)/2) for j in range(1,1+(s-1)/2+1) for k in range(1,1+t/2))))
        if r % 2 == 1 and s % 2 == 0 and t % 2 == 1:
            return Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+(r-1)/2+1) for j in range(1,1+s/2) for k in range(1,1+(t-1)/2)))) * Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+(r-1)/2) for j in range(1,1+s/2) for k in range(1,1+(t-1)/2+1))))
        if r % 2 == 0 and s % 2 == 1 and t % 2 == 1:
            return Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+r/2) for j in range(1,1+(s-1)/2+1) for k in range(1,1+(t-1)/2)))) * Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k - 2) for i in range(1,1+r/2) for j in range(1,1+(s-1)/2) for k in range(1,1+(t-1)/2+1))))
        if r % 2 == 1 and s % 2 == 1 and t % 2 == 1:
            return Integer(0)


#Class 6

class PlanePartitions_TCPP(PlanePartitions):
    @staticmethod
    def __classcall_private__(cls, box_size):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: P1 = PlanePartitions((4,3,2), symmetry='TCPP')
            sage: P2 = PlanePartitions([4,3,2], symmetry='TCPP')
            sage: P1 is P2
            True
        """
        return super(PlanePartitions_TCPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([3,3,3], symmetry='TCPP')
            sage: TestSuite(PP).run()
        """
        if (box_size[0] % 2 == 1 and box_size[1] % 2 == 1 and box_size[2] % 2 == 1):
            raise ValueError("x, y, and z dimensions cannot all be odd")
        super(PlanePartitions_TCPP, self).__init__(category=FiniteEnumeratedSets())
        self._box = box_size
        self._symmetry = 'TCPP'

    def _repr_(self) -> str:
        return "Transpose complement plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __iter__(self) -> Iterator:
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((3,3,2), symmetry='TCPP'))
            [Plane partition [[2, 2, 1], [2, 1], [1]],
            Plane partition [[2, 1, 1], [2, 1, 1], [1]],
            Plane partition [[2, 2, 1], [1, 1], [1, 1]],
            Plane partition [[2, 1, 1], [1, 1, 1], [1, 1]],
            Plane partition [[1, 1, 1], [1, 1, 1], [1, 1, 1]]]
        """
        for p in PlanePartitions(self._box):
            if p.is_TCPP():
                yield self.element_class(self,p)
        return



    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of transpose complement plane partitions inside an 
        `a \times a \times 2b` box is equal to

        .. MATH::

            \binom{b+1-1}{a-1}\prod_{1\leq i\leq j a-2} \frac{i + j + 2b - 1}{i + j - 1}
        """

        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        if a==b and c % 2 == 0:
            return Integer(binomial(c/2+a-1, a-1) * prod((c + i + j + 1)/(i + j + 1) for j in range(1,1+a-2) for i in range(1,1+j)))
        else:
            return Integer(0)

#Class 7

class PlanePartitions_SSCPP(PlanePartitions):
#Symmetric self-complementary plane partitions
    @staticmethod
    def __classcall_private__(cls, box_size):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: P1 = PlanePartitions((4,4,2), symmetry='SSCPP')
            sage: P2 = PlanePartitions([4,4,2], symmetry='SSCPP')
            sage: P1 is P2
            True
        """
        return super(PlanePartitions_SSCPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([4,4,2], symmetry='SSCPP')
            sage: TestSuite(PP).run()
        """


        if box_size[0]!=box_size[1]:
            raise ValueError("x and y dimensions must be equal")
        if (box_size[0] % 2 == 1 and box_size[1] % 2 == 1 and box_size[2] % 2 == 1):
            raise ValueError("x, y, and z dimensions cannot all be odd")
        super(PlanePartitions_SSCPP, self).__init__(category=FiniteEnumeratedSets())
        self._box = box_size
        self._symmetry = 'SSCPP'

    def _repr_(self) -> str:
        return "Symmetric self-complementary plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __iter__(self) -> Iterator:
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((3,3,2), symmetry='SSCPP'))
            [Plane partition [[2, 2, 1], [2, 1], [1]],
            Plane partition [[2, 1, 1], [1, 1, 1], [1, 1]],
            Plane partition [[1, 1, 1], [1, 1, 1], [1, 1, 1]]]
        """
        for p in PlanePartitions(self._box):
            if p.is_SSCPP():
                yield self.element_class(self,p)
        return

    def cardinality(self) -> Integer:
        r=self._box[0]
        s=self._box[1]
        t=self._box[2]
        if r % 2 == 0 and s % 2 == 0 and t % 2 == 0:
            return Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k -2 ) for i in range(1,1+r/2) for j in range(1,1+r/2) for k in range(1,1+t/2))))
        if r % 2 == 1 and s % 2 == 1 and t % 2 == 0:
            return Integer((prod( Integer(i + j + k - 1) / Integer(i + j + k -2 ) for i in range(1,1+(r-1)/2) for j in range(1,1+((r-1)/2)+1) for k in range(1,1+t/2))))
        return Integer(0)

#Class 8

class PlanePartitions_CSTCPP(PlanePartitions):
#Cyclically symmetric transpose complement partitions
    @staticmethod
    def __classcall_private__(cls, box_size):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: P1 = PlanePartitions((3,3,3), symmetry='CSTCPP')
            sage: P2 = PlanePartitions([3,3,3], symmetry='CSTCPP')
            sage: P1 is P2
            True
        """
        return super(PlanePartitions_CSTCPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([2,2,2], symmetry='CSTCPP')
            sage: TestSuite(PP).run()
        """
        if box_size[0] != box_size[1] or box_size[1] != box_size[2]:
            raise ValueError("x, y, and z dimensions must match")
        super(PlanePartitions_CSTCPP, self).__init__(category=FiniteEnumeratedSets())
        self._box = box_size
        self._symmetry = 'CSTCPP'

    def _repr_(self) -> str:
        return "Cyclically symmetric transpose complement partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __iter__(self) -> Iterator:
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((2,2,2), symmetry='CSTCPP'))
            [Plane partition [[2, 1], [1]]]
        """
        for p in PlanePartitions(self._box):
            if p.is_CSTCPP():
                yield self.element_class(self,p)
        return

    def cardinality(self) -> Integer:
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        if a == b and b == c and a % 2 == 0:
            return Integer(prod( ((3*i+1)*factorial(6*i)*factorial(2*i))/(factorial(4*i+1)*factorial(4*i)) for i in range(1+(a/2)-1)))
        return Integer(0)

# Class 9


class PlanePartitions_CSSCPP(PlanePartitions):
# Cyclically symmetric self-complementary plane partitions
    @staticmethod
    def __classcall_private__(cls, box_size):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: P1 = PlanePartitions((3,3,3), symmetry='CSSCPP')
            sage: P2 = PlanePartitions([3,3,3], symmetry='CSSCPP')
            sage: P1 is P2
            True
        """
        return super(PlanePartitions_CSSCPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([3,3,3], symmetry='CSSCPP')
            sage: TestSuite(PP).run()
        """
        if box_size[0] != box_size[1] or box_size[1] != box_size[2]:
            raise ValueError("x, y, and z dimensions must match")
        if (box_size[0] % 2 == 1 and box_size[1] % 2 == 1 and box_size[2] % 2 == 1):
            raise ValueError("x, y, and z dimensions cannot all be odd")  
        super(PlanePartitions_CSSCPP, self).__init__(category=FiniteEnumeratedSets())
        self._box = box_size
        self._symmetry = 'CSSCPP'

    def _repr_(self) -> str:
        return "Cyclically symmetric self-complementary plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __iter__(self) -> Iterator:
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((2,2,2), symmetry='CSSCPP'))
            [Plane partition [[2, 1], [1]]]
        """
        for p in PlanePartitions(self._box):
            if p.is_CSSCPP():
                yield self.element_class(self,p)
        return

    def cardinality(self) -> Integer:
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        if a == b and b == c and a % 2 == 0:
            return Integer(prod( ((factorial(3*i+1)**2/(factorial(a+i)**2) for i in range((a/2))))))
        return Integer(0)



#Class 10

class PlanePartitions_TSSCPP(PlanePartitions):
#Totally symmetric self-complementary plane partitions
    @staticmethod
    def __classcall_private__(cls, box_size):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: P1 = PlanePartitions((4,4,4), symmetry='TSSCPP')
            sage: P2 = PlanePartitions([4,4,4], symmetry='TSSCPP')
            sage: P1 is P2
            True
        """
        return super(PlanePartitions_TSSCPP, cls).__classcall__(cls, tuple(box_size))

    def __init__(self, box_size):
        """
        TESTS::
    
            sage: PP = PlanePartitions([4,4,4], symmetry='TSSCPP')
            sage: TestSuite(PP).run()
        """
        if (box_size[0] != box_size[1] or (box_size[1] != box_size[2]) or box_size[0] % 2 != 0):
            raise ValueError("x, y, and z dimensions must be (2r,2r,2r)")
        super(PlanePartitions_TSSCPP, self).__init__(category=FiniteEnumeratedSets())
        self._box = box_size
        self._symmetry = 'TSSCPP'

    def _repr_(self) -> str:
        return "Totally symmetric self-complementary plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def to_poset(self):
        cmp = lambda x,y : all(x[i] <= y[i] for i in range(len(x)))
        def componentwise_comparer(thing1,thing2):
            if len(thing1) == len(thing2):
                if all(thing1[i] <= thing2[i] for i in range(len(thing1))):
                    return True
            return False
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        if a != b or b != c or a != c:
            return

        pl = []
        for x in range(0,a/2 - 2 + 1):
            for y in range(x, a/2 - 2 + 1):
                    for z in range(0,a/2 - 2 + 1):
                        if z <= a/2 - 2 - y:
                            pl.append((x,y,z))

        return Poset((pl,cmp))
        
    def from_antichain(self, acl) -> PP:
        #ac format ex: [x,y,z]
        a=self._box[0]
        b=self._box[1]
        c=self._box[2]
        n=a
        ppMatrix = [[0] * (c) for i in range(b)] #creates a matrix for the plane parition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]
        width = n/2 - 1
        height = n/2 - 1

        #generate inner triangle
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
        return self.element_class(self, ppMatrix)
        




    def __iter__(self) -> Iterator:
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: list(PlanePartitions((4,4,4), symmetry='TSSCPP'))
            [Plane partition [[4, 4, 2, 2], [4, 4, 2, 2], [2, 2], [2, 2]],
            Plane partition [[4, 4, 3, 2], [4, 3, 2, 1], [3, 2, 1], [2, 1]]]
        """
#        def componentwise_comparer(thing1,thing2):
#            if len(thing1) == len(thing2):
#                if all(thing1[i] <= thing2[i] for i in range(len(thing1))):
#                    return True
#            return False
#        a=self._box[0]
#        b=self._box[1]
#        c=self._box[2]
#        n = a
#        b = n
#        c = n

#        pl = []
#        for x in range(0,n/2 - 2 + 1):
#            for y in range(x, n/2 - 2 + 1):
#                    for z in range(0,n/2 - 2 + 1):
#                        if z <= n/2 - 2 - y:
#                            pl.append((x,y,z))

#        pocp = Poset((pl,componentwise_comparer))
#        cmp = lambda x,y : all(x[i] <= y[i] for i in range(len(x)))
#        pocp = Poset((pl,cmp))

#        matrixList = [] #list of all PlaneParitions with parameters(a,b,c)
        #iterate through each antichain of product of chains poset with paramaters (a,b,c)
        for acl in self.to_poset().antichains_iterator():
            yield self.from_antichain(acl)
        return
#            #ac format ex: [x,y,z]
#            ppMatrix = [[0] * (c) for i in range(b)] #creates a matrix for the plane parition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]
#            width = n/2 - 1
#            height = n/2 - 1

#            #generate inner triagle
#            for i in range(width):
#                for j in range(height):
#                    if(i <= j):
#                        for ac in acl:
#                            if ac[0] == i and ac[1] == j:
#                                zVal = ac[2]
#                                matrixVal = ppMatrix[j +(n/2)] [i+ (n/2)]
#                                if zVal + 1 > matrixVal:
#                                    ppMatrix[j +(n/2)] [i+ (n/2)]= zVal + 1

#            #fill back
#            for i in range(width):
#                i = width-(i+1)
#                i = i + n/2
#                for j in range(height):
#                    j = height-(j+1)
#                    j = j + n/2
#                    if (ppMatrix[i][j] == 0):
#                        if i >= j:
#                            iValue = 0
#                            jValue = 0
#                            if i < n:
#                                iValue = ppMatrix[i+1][j]
#                            if j < n:
#                               jValue = ppMatrix[i][j+1]
#                            ppMatrix[i][j] = max(iValue,jValue)


#            #fill half of triangle symmetrically
#            for i in range(width):
#                i = i + n/2
#                for j in range(height):
#                    j = j + n/2
#                    if i >= j:
#                        ppMatrix[j][i] = ppMatrix[i][j]

#            #upper left box
#            for i in range(n/2):
#                for j in range(n/2):
#                    ppMatrix[i][j] = n - ppMatrix[n-(i+1)][n-(j+1)]


#            #fill in lower left cube with values n/2
#            for i in range(n/2):
#                for j in range(n/2):
#                    x = i
#                    y = j
#                    if(ppMatrix[x][y+(n/2)]) == 0:
#                        ppMatrix[x][y+(n/2)] = n/2
#                    if(ppMatrix[x+(n/2)][y]) == 0:
#                        ppMatrix[x+(n/2)][y] = n/2


#            #add and subtract values from lower left cube to be rotation of lower right cube
#            for i in range(n/2):
#                for j in range(n/2):
#                    x = i+(n/2)
#                    y = j+(n/2)
#                    if ppMatrix[x][y] > 0:
#                        z = ppMatrix[x][y]
#                        for cVal in range(z):
#                            #build onto lower left cube
#                            ppMatrix[x][0+cVal] += 1
#                            #carve out of lower left cube
#                            ppMatrix[n-(1+cVal)][(n/2)-(j+1)] -=1

#            #fill in upper right cube symmetrically with lower left
#            for i in range(n/2):
#                for j in range(n/2):
#                    ppMatrix[j][i+(n/2)] = ppMatrix[i+(n/2)][j]
#            yield self.element_class(self, ppMatrix)


