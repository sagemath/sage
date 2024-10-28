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
#                     2021 Kevin Dilks <kdilks@gmail.com>
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
from typing import NewType
from collections.abc import Iterator

from sage.structure.richcmp import richcmp, richcmp_method
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.tableau import Tableau
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.rings.integer import Integer
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.arith.misc import Sigma, binomial, factorial
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers

lazy_import('sage.modules.free_module_element', 'vector')


@richcmp_method
class PlanePartition(ClonableArray,
                     metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A plane partition.

    A *plane partition* is a stack of cubes in the positive orthant.

    INPUT:

    - ``PP`` -- list of lists which represents a tableau
    - ``box_size`` -- (optional) a list ``[A, B, C]`` of 3 positive integers,
      where ``A``, ``B``, ``C`` are the lengths of the box in the `x`-axis,
      `y`-axis, `z`-axis, respectively; if this is not given, it is
      determined by the smallest box bounding ``PP``

    OUTPUT: the plane partition whose tableau representation is ``PP``

    EXAMPLES::

        sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
        sage: PP
        Plane partition [[4, 3, 3, 1], [2, 1, 1], [1, 1]]

    TESTS::

        sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
        sage: TestSuite(PP).run()
        sage: hash(PP) # random
    """
    @staticmethod
    def __classcall_private__(cls, PP, box_size=None):
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
        if isinstance(PP, PlanePartition) and box_size is None:
            return PP
        pp = PlanePartitions(box_size=box_size)
        return pp.element_class(pp, PP)  # The check() will raise the appropriate error

    def __init__(self, parent, pp, check=True):
        r"""
        Initialize a plane partition.

        TESTS::

            sage: a = PlanePartitions()([[2,1],[1]])
            sage: b = PlanePartitions([2,2,2])([[2,1],[1]])
            sage: c = PlanePartitions(4)([[2,1],[1]])
            sage: a == b
            True
            sage: a is b
            False
            sage: a == c
            True
            sage: a is c
            False
        """
        if isinstance(pp, PlanePartition):
            ClonableArray.__init__(self, parent, pp, check=False)
        else:
            pp = [list(row) for row in pp]
            if pp:
                for i in reversed(range(len(pp))):
                    while pp[i] and not pp[i][-1]:
                        del pp[i][-1]
                    if not pp[i]:
                        pp.pop(i)
            pp = [tuple(row) for row in pp]
            ClonableArray.__init__(self, parent, pp, check=check)
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

    def __richcmp__(self, other, op):
        r"""
        Compare ``self`` to ``other``.

        .. TODO::

            This overwrites the comparison check of
            :class:`~sage.structure.list_clone.ClonableArray`
            in order to circumvent the coercion framework.
            Eventually this should be solved more elegantly,
            for example along the lines of what was done for
            `k`-tableaux.

            For now, this compares two elements by their underlying
            defining lists.

        INPUT:

        - ``other`` -- the element that ``self`` is compared to

        OUTPUT: boolean

        TESTS::

            sage: t = PlanePartition([[2,1],[1]])
            sage: t == 0
            False
            sage: t == PlanePartitions(4)([[2,1],[1]])
            True

            sage: s = PlanePartition([[3,1],[1]])
            sage: s != []
            True

            sage: t < s
            True
            sage: s < t
            False
            sage: s > t
            True
        """
        if isinstance(other, PlanePartition):
            return self._richcmp_(other, op)

        return richcmp(list(self), other, op)

    def check(self):
        """
        Check to see that ``self`` is a valid plane partition.

        EXAMPLES::

            sage: a = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: a.check()
            sage: b = PlanePartition([[1,2],[1]])
            Traceback (most recent call last):
            ...
            ValueError: not weakly decreasing along rows
            sage: c = PlanePartition([[1,1],[2]])
            Traceback (most recent call last):
            ...
            ValueError: not weakly decreasing along columns
            sage: d = PlanePartition([[2,-1],[-2]])
            Traceback (most recent call last):
            ...
            ValueError: entries not all nonnegative
            sage: e = PlanePartition([[3/2,1],[.5]])
            Traceback (most recent call last):
            ...
            ValueError: entries not all integers
        """
        if not all(a in ZZ for b in self for a in b):
            raise ValueError("entries not all integers")
        for row in self:
            if not all(c >= 0 for c in row):
                raise ValueError("entries not all nonnegative")
            if not all(row[i] >= row[i+1] for i in range(len(row)-1)):
                raise ValueError("not weakly decreasing along rows")
        for row, next in zip(self, self[1:]):
            if not all(row[c] >= next[c] for c in range(len(next))):
                raise ValueError("not weakly decreasing along columns")

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            Plane partition [[4, 3, 3, 1], [2, 1, 1], [1, 1]]
        """
        return "Plane partition {}".format([list(row) for row in self])

    def to_tableau(self) -> Tableau:
        r"""
        Return the tableau class of ``self``.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.to_tableau()
            [[4, 3, 3, 1], [2, 1, 1], [1, 1]]
        """
        return Tableau(self)  # type:ignore

    def z_tableau(self, tableau=True) -> Tableau:
        r"""
        Return the projection of ``self`` in the `z` direction.

        If ``tableau`` is set to ``False``, then only the list of lists
        consisting of the projection of boxes size onto the `xy`-plane
        is returned instead of a :class:`Tableau` object. This output will
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
        consisting of the projection of boxes size onto the `xz`-plane
        is returned instead of a :class:`Tableau` object. This output will
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
        consisting of the projection of boxes size onto the `yz`-plane
        is returned instead of a :class:`Tableau` object. This output will
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

    def cells(self) -> list[tuple[int, int, int]]:
        r"""
        Return the list of cells inside ``self``.

        Each cell is a tuple.

        EXAMPLES::

            sage: PP = PlanePartition([[3,1],[2]])
            sage: PP.cells()
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 0), (1, 0, 0), (1, 0, 1)]
        """
        return [(r, c, h)
                for r in range(len(self))
                for c in range(len(self[r]))
                for h in range(self[r][c])]

    def number_of_boxes(self) -> Integer:
        r"""
        Return the number of boxes in the plane partition.

        EXAMPLES::

            sage: PP = PlanePartition([[3,1],[2]])
            sage: PP.number_of_boxes()
            6
        """
        return sum(sum(row) for row in self)

    def _repr_diagram(self, show_box=False, use_unicode=False) -> str:
        r"""
        Return a string of the 3D diagram of ``self``.

        INPUT:

        - ``show_box`` -- boolean (default: ``False``); if ``True``,
          also shows the visible tiles on the `xy`-, `yz`-, `zx`-planes
        - ``use_unicode`` -- boolean (default: ``False``); use unicode

        OUTPUT: string of the 3D diagram of the plane partition

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

        hori = "_" if use_unicode else "_"
        down = "╲" if use_unicode else "\\"
        up = "╱" if use_unicode else "/"

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
            return "∅" if use_unicode else ""

        if use_unicode:
            return '\n'.join("".join(row) for row in drawing)
        return '\n'.join("".join(row) for row in drawing)

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

        OUTPUT: a pretty print of the plane partition

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

    def _repr_svg_(self) -> str:
        """
        Return the svg picture of a plane partition.

        This can be displayed by Jupyter.

        EXAMPLES::

            sage: PP = PlanePartition([[2, 1, 1], [1, 1]])
            sage: PP._repr_svg_()                                                       # needs sage.modules
            '<?xml...</g></svg>'
        """
        colors = ["snow", "tomato", "steelblue"]

        resu = '<?xml version=\"1.0\" standalone=\"no\"?>'
        resu += '<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" '
        resu += '\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">'
        resu += '<svg xmlns=\"http://www.w3.org/2000/svg\" '
        resu += 'xmlns:xlink=\"http://www.w3.org/1999/xlink\" width=\"300\" viewBox='

        resu1 = '<defs><polygon points=\"0, 0 -0.866, 0.5 0, 1 0.866, 0.5\" '
        resu1 += f'id=\"cz\" style=\"fill:{colors[0]}\"/>'
        resu1 += '<polygon points=\"0, 0 0.866, 0.5 0.866, -0.5 0, -1\" '
        resu1 += f'id=\"cx\" style=\"fill:{colors[1]}\"/>'
        resu1 += '<polygon points=\"0, 0 0, -1 -0.866, -0.5 -0.866, 0.5\" '
        resu1 += f'id=\"cy\" style=\"fill:{colors[2]}\"/></defs>'
        resu1 += '<g style=\"stroke-width:0.01;stroke-linejoin:bevel; '
        resu1 += 'stroke-linecap:butt; stroke:black; fill:red\">'

        vx = -vector([0.866, -0.5])
        vy = -vector([-0.866, -0.5])
        vz = -vector([0, 1])
        # Since we currently don't display the bounding box, just
        #   use the smallest one possible.
        Nx, Ny, Nz = self.bounding_box()

        resu += '\"%.3f %.3f %.3f %.3f \">' % (-0.866 * Nx, -Nz,
                                               0.866 * Nx + 0.866 * Ny,
                                               Nz + 0.5 * (Nx + Ny))
        resu += resu1

        mat = self.z_tableau()
        for i in range(Nx):
            for j in range(Ny):
                if mat[i][j]:
                    v = i * vx + j * vy + mat[i][j] * vz
                    resu += '<use transform=\"translate(%.3f, %.3f)' % (v[0], v[1])
                    resu += '\" xlink:href=\"#cz\" />'

        mat = self.y_tableau()
        for j in range(Nz):
            for k in range(Nx):
                if mat[j][k]:
                    v = j * vz + k * vx + mat[j][k] * vy
                    resu += '<use transform=\"translate(%.3f, %.3f)' % (v[0], v[1])
                    resu += '\" xlink:href=\"#cy\" />'

        mat = self.x_tableau()
        for k in range(Ny):
            for i in range(Nz):
                if mat[k][i]:
                    v = k * vy + i * vz + mat[k][i] * vx
                    resu += '<use transform=\"translate(%.3f, %.3f)' % (v[0], v[1])
                    resu += '\" xlink:href=\"#cx\" />'
        return resu + '</g></svg>'

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

        OUTPUT: latex code for drawing the plane partition

        EXAMPLES::

            sage: PP = PlanePartition([[1]])
            sage: latex(PP)                                                             # needs sage.graphs
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
            sage: PP.plot()                                                             # needs sage.plot
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
            return polygon(move(Uside, i, j, k), edgecolor='black',
                           color=colors[0])

        def add_leftside(i, j, k):
            return polygon(move(Lside, i, j, k), edgecolor='black',
                           color=colors[1])

        def add_rightside(i, j, k):
            return polygon(move(Rside, i, j, k), edgecolor='black',
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
        r"""
        Return ``True`` if ``PP`` is a plane partition that fits
        inside ``self``.

        Specifically, ``self`` contains ``PP`` if, for all `i`, `j`,
        the height of ``PP`` at `ij` is less than or equal to the
        height of ``self`` at `ij`.

        EXAMPLES::

            sage: P1 = PlanePartition([[5,4,3], [3,2,2], [1]])
            sage: P2 = PlanePartition([[3,2], [1,1], [0,0], [0,0]])
            sage: P3 = PlanePartition([[5,5,5], [2,1,0]])
            sage: P1.contains(P2)
            True
            sage: P2.contains(P1)
            False
            sage: P1.contains(P3)
            False
            sage: P3.contains(P2)
            True
        """
        if not isinstance(PP, PlanePartition):
            PP = PlanePartition(PP)
        if len(self) < len(PP):
            return False
        for rowself, rowPP in zip(self, PP):
            if len(rowself) < len(rowPP):
                return False
            if any(valself < valPP for valself, valPP in zip(rowself, rowPP)):
                return False
        return True

    def plot3d(self, colors=None):
        r"""
        Return a 3D-plot of ``self``.

        INPUT:

        - ``colors`` -- (default: ``["white", "lightgray", "darkgray"]``)
          list ``[A, B, C]`` of 3 strings representing colors

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.plot3d()                                                           # needs sage.plot
            Graphics3d Object
        """
        if colors is None:
            colors = ["white", "lightgray", "darkgray"]
        from sage.plot.plot3d.platonic import cube
        return sum(cube(c, color=colors, frame_thickness=2,
                        frame_color='black', frame=False)
                   for c in self.cells())

    def complement(self, tableau_only=False) -> PP:
        r"""
        Return the complement of ``self``.

        If the parent of ``self`` consists only of partitions inside a given
        box, then the complement is taken in this box. Otherwise, the
        complement is taken in the smallest box containing the plane partition.
        The empty plane partition with no box specified is its own complement.

        If ``tableau_only`` is set to ``True``, then only the tableau
        consisting of the projection of boxes size onto the `xy`-plane
        is returned instead of a :class:`PlanePartition`. This output will
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
                T[A-1-r][B-1-c] = C - z_tab[r][c]
        if tableau_only:
            return T
        P = self.parent()
        if not P._box:
            pp = PlanePartitions()
            return pp.element_class(pp, T)
        return P.element_class(P, T, check=False)

    def transpose(self, tableau_only=False) -> PP:
        r"""
        Return the transpose of ``self``.

        If ``tableau_only`` is set to ``True``, then only the tableau
        consisting of the projection of boxes size onto the `xy`-plane
        is returned instead of a :class:`PlanePartition`. This will
        not necessarily have trailing rows or trailing zeros removed.

        EXAMPLES::

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.transpose()
            Plane partition [[4, 2, 1], [3, 1, 1], [3, 1], [1]]
            sage: PP.transpose(True)
            [[4, 2, 1], [3, 1, 1], [3, 1, 0], [1, 0, 0]]

            sage: PPP = PlanePartitions([1, 2, 3])
            sage: PP = PPP([[1, 1]])
            sage: PT = PP.transpose(); PT
            Plane partition [[1], [1]]
            sage: PT.parent()
            Plane partitions inside a 2 x 1 x 3 box
        """
        T = [[0 for i in range(self._max_x)] for j in range(self._max_y)]
        z_tab = self.z_tableau()
        for r in range(len(z_tab)):
            for c in range(len(z_tab[r])):
                T[c][r] = z_tab[r][c]
        P = self.parent()
        if tableau_only:
            return T
        elif P._box is None or P._box[0] == P._box[1]:
            return P.element_class(P, T, check=False)
        new_box = (P._box[1], P._box[0], P._box[2])
        newP = PlanePartitions(new_box, symmetry=P._symmetry)
        return newP.element_class(newP, T)

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

        TESTS::

            sage: PlanePartition([]).is_SPP()
            True
        """
        if not self:
            return True
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

        TESTS::

            sage: PlanePartition([]).is_CSPP()
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

        TESTS::

            sage: PlanePartition([]).is_TSPP()
            True
        """
        return self.is_CSPP() and self.is_SPP()

    def is_SCPP(self) -> bool:
        r"""
        Return whether ``self`` is a self-complementary plane partition.

        Note that the complement of a plane partition (and thus the property of
        being self-complementary) is dependent on the choice of a box that it is
        contained in. If no parent/bounding box is specified,  the box is taken
        to be the smallest box that contains the plane partition.

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

        TESTS::

            sage: PlanePartition([]).is_SCPP()
            True
        """
        return self.z_tableau(tableau=False) == self.complement(tableau_only=True)

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

        TESTS::

            sage: PlanePartition([]).is_TCPP()
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

        TESTS::

            sage: PlanePartition([]).is_SSCPP()
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

        TESTS::

            sage: PlanePartition([]).is_CSTCPP()
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

        TESTS::

            sage: PlanePartition([]).is_CSSCPP()
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

        TESTS::

            sage: PlanePartition([]).is_TSSCPP()
            True
        """
        return self.is_TSPP() and self.is_SCPP()

    def to_order_ideal(self):
        r"""
        Return the order ideal corresponding to ``self``.

        .. TODO::

            As many families of symmetric plane partitions are in bijection
            with order ideals in an associated poset, this function could
            feasibly have options to send symmetric plane partitions
            to the associated order ideal in that poset, instead.

        EXAMPLES::

            sage: PlanePartition([[3,2,1],[2,2],[2]]).to_order_ideal()                  # needs sage.graphs sage.modules
            [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 0), (0, 1, 1), (0, 2, 0),
             (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1), (2, 0, 0), (2, 0, 1)]
            sage: PlanePartition([[2,1],[1],[1]]).to_order_ideal()                      # needs sage.graphs sage.modules
            [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (2, 0, 0)]
        """
        from sage.combinat.posets.poset_examples import posets
        (a, b, c) = (self._max_x, self._max_y, self._max_z)
        Q = posets.ProductOfChains([a, b, c])
        count = 0
        generate = []
        for i, row in enumerate(self):
            for j, val in enumerate(row):
                if val > 0:
                    generate.append((i, j, val-1))
            count += 1
        oi = Q.order_ideal(generate)
        return oi

    def maximal_boxes(self) -> list:
        r"""
        Return the coordinates of the maximal boxes of ``self``.

        The maximal boxes of a plane partitions are the boxes that can be
        removed from a plane partition and still yield a valid plane partition.

        EXAMPLES::

            sage: sorted(PlanePartition([[3,2,1],[2,2],[2]]).maximal_boxes())
            [[0, 0, 2], [0, 2, 0], [1, 1, 1], [2, 0, 1]]
            sage: sorted(PlanePartition([[2,1],[1],[1]]).maximal_boxes())
            [[0, 0, 1], [0, 1, 0], [2, 0, 0]]
        """
        generate = []
        for i, row in enumerate(self):
            for j, entry in enumerate(row):
                if (i == len(self)-1 or len(self[i+1])-1 < j or self[i+1][j] < entry) and (j == len(row)-1 or row[j+1] < entry):
                    generate.append([i, j, entry-1])
        return generate

    def cyclically_rotate(self, preserve_parent=False) -> PP:
        r"""
        Return the cyclic rotation of ``self``.

        By default, if the parent of ``self`` consists of plane
        partitions inside an `a \times b \times c` box, the result
        will have a parent consisting of partitions inside
        a `c \times a \times b` box, unless the optional parameter
        ``preserve_parent`` is set to ``True``. Enabling this setting
        may give an element that is **not** an element of its parent.

        EXAMPLES::

            sage: PlanePartition([[3,2,1],[2,2],[2]]).cyclically_rotate()
            Plane partition [[3, 3, 1], [2, 2], [1]]
            sage: PP = PlanePartition([[4,1],[1],[1]])
            sage: PP.cyclically_rotate()
            Plane partition [[3, 1, 1, 1], [1]]
            sage: PP == PP.cyclically_rotate().cyclically_rotate().cyclically_rotate()
            True

            sage: # needs sage.graphs sage.modules
            sage: PP = PlanePartitions([4,3,2]).random_element()
            sage: PP.cyclically_rotate().parent()
            Plane partitions inside a 2 x 4 x 3 box
            sage: PP = PlanePartitions([3,4,2])([[2,2,2,2],[2,2,2,2],[2,2,2,2]])
            sage: PP_rotated = PP.cyclically_rotate(preserve_parent=True)
            sage: PP_rotated in PP_rotated.parent()
            False
        """
        b = self._max_y
        c = self._max_z
        new_antichain = []
        for elem in self.maximal_boxes():
            new = (elem[1], elem[2], elem[0])
            new_antichain.append(new)
        pp_matrix = [[0] * (c) for i in range(b)]
        for box in new_antichain:
            y = box[0]
            z = box[1]
            x = box[2]
            pp_matrix[y][z] = x + 1
        if new_antichain:
            for i in range(b):
                i = b - (i+1)
                for j in range(c):
                    j = c - (j+1)
                    if pp_matrix[i][j] == 0:
                        iValue = 0
                        jValue = 0
                        if i < b-1:
                            iValue = pp_matrix[i+1][j]
                        if j < c-1:
                            jValue = pp_matrix[i][j+1]
                        pp_matrix[i][j] = max(iValue, jValue)
        # Start code for determining correct parent
        P = self.parent()
        if P._box is None or preserve_parent or (P._box[0] == P._box[1] == P._box[2]):
            return P.element_class(P, pp_matrix, check=preserve_parent)
        new_box = (P._box[2], P._box[0], P._box[1])
        newP = PlanePartitions(new_box, symmetry=P._symmetry)
        return newP.element_class(newP, pp_matrix)

    def bounding_box(self):
        r"""
        Return the smallest box `(a, b, c)` that ``self`` is contained in.

        EXAMPLES::

            sage: PP = PlanePartition([[5,2,1,1], [2,2], [2]])
            sage: PP.bounding_box()
            (3, 4, 5)
        """
        if not self:
            return (0, 0, 0)
        return (len(self), len(self[0]), self[0][0])


PP = NewType('PP', PlanePartition)


class PlanePartitions(UniqueRepresentation, Parent):
    r"""
    Plane partitions.

    ``PlanePartitions()`` returns the class of all plane partitions.

    ``PlanePartitions(n)`` return the class of all plane partitions with
    precisely `n` boxes.

    ``PlanePartitions([a, b, c])`` returns the class of plane partitions
    that fit inside an `a \times b \times c` box.

    ``PlanePartitions([a, b, c])`` has the optional keyword ``symmetry``, which
    restricts the plane partitions inside a box of the specified size satisfying
    certain symmetry conditions.

    - ``symmetry='SPP'`` gives the class of symmetric plane partitions. which
      is all plane partitions fixed under reflection across the diagonal.
      Requires that `a = b`.

    - ``symmetry='CSPP'`` gives the class of cyclic plane partitions, which
      is all plane partitions fixed under cyclic rotation of coordinates.
      Requires that `a = b = c`.

    - ``symmetry='TSPP'`` gives the class of totally symmetric plane partitions,
      which is all plane partitions fixed under any interchanging of coordinates.
      Requires that `a = b = c`.

    - ``symmetry='SCPP'`` gives the class of self-complementary plane partitions.
      which is all plane partitions that are equal to their own complement
      in the specified box. Requires at least one of `a,b,c` be even.

    - ``symmetry='TCPP'`` gives the class of transpose complement plane
      partitions, which is all plane partitions whose complement in the box
      of the specified size is equal to their transpose. Requires `a = b` and
      at least one of `a, b, c` be even.

    - ``symmetry='SSCPP'`` gives the class of symmetric self-complementary
      plane partitions, which is all plane partitions that are both
      symmetric and self-complementary. Requires `a = b` and at least one of
      `a, b, c` be even.

    - ``symmetry='CSTCPP'`` gives the class of cyclically symmetric transpose
      complement plane partitions, which is all plane partitions that are
      both symmetric and equal to the transpose of their complement. Requires
      `a = b = c`.

    - ``symmetry='CSSCPP'`` gives the class of cyclically symmetric
      self-complementary plane partitions, which is all plane partitions that
      are both cyclically symmetric and self-complementary. Requires `a = b = c`
      and all `a, b, c` be even.

    - ``symmetry='TSSCPP'`` gives the class of totally symmetric
      self-complementary plane partitions, which is all plane partitions that
      are totally symmetric and also self-complementary. Requires `a = b = c`
      and all `a, b, c` be even.

    EXAMPLES:

    If no arguments are passed, then the class of all plane partitions
    is returned::

        sage: PlanePartitions()
        Plane partitions
        sage: [[2,1],[1]] in PlanePartitions()
        True

    If an integer `n` is passed, then the class of plane partitions of `n`
    is returned::

        sage: PlanePartitions(3)
        Plane partitions of size 3
        sage: PlanePartitions(3).list()
        [Plane partition [[3]],
         Plane partition [[2, 1]],
         Plane partition [[1, 1, 1]],
         Plane partition [[2], [1]],
         Plane partition [[1, 1], [1]],
         Plane partition [[1], [1], [1]]]

    If a three-element tuple or list `[a,b,c]` is passed, then the class of all
    plane partitions that fit inside and `a \times b \times c` box is returned::

        sage: PlanePartitions([2,2,2])
        Plane partitions inside a 2 x 2 x 2 box
        sage: [[2,1],[1]] in PlanePartitions([2,2,2])
        True

    If an additional keyword ``symmetry`` is pass along with a three-element
    tuple or list `[a, b,c ]`, then the class of all plane partitions that fit
    inside an `a \times b \times c` box with the specified symmetry is returned::

        sage: PlanePartitions([2,2,2], symmetry='CSPP')
        Cyclically symmetric plane partitions inside a 2 x 2 x 2 box
        sage: [[2,1],[1]] in PlanePartitions([2,2,2], symmetry='CSPP')
        True

    .. SEEALSO::

        - :class:`PlanePartition`
        - :class:`PlanePartitions_all`
        - :class:`PlanePartitions_n`
        - :class:`PlanePartitions_box`
        - :class:`PlanePartitions_SPP`
        - :class:`PlanePartitions_CSPP`
        - :class:`PlanePartitions_TSPP`
        - :class:`PlanePartitions_SCPP`
        - :class:`PlanePartitions_TCPP`
        - :class:`PlanePartitions_SSCPP`
        - :class:`PlanePartitions_CSTCPP`
        - :class:`PlanePartitions_CSSCPP`
        - :class:`PlanePartitions_TSSCPP`
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwds):
        r"""
        Return the appropriate parent based on arguments.

        See the documentation for :class:`PlanePartitions` for more information.

        TESTS::

            sage: PlanePartitions()
            Plane partitions
            sage: PlanePartitions([3,3,3])
            Plane partitions inside a 3 x 3 x 3 box
            sage: PlanePartitions(3)
            Plane partitions of size 3
            sage: PlanePartitions([4,4,4], symmetry='TSSCPP')
            Totally symmetric self-complementary plane partitions inside a 4 x 4 x 4 box
            sage: PlanePartitions(4, symmetry='TSSCPP')
            Traceback (most recent call last):
            ...
            ValueError: the number of boxes may only be specified if no symmetry is required
        """
        symmetry = kwds.get('symmetry', None)
        box_size = kwds.get('box_size', None)

        if not args and symmetry is None and box_size is None:
            return PlanePartitions_all()

        if args and box_size is None:
            # The first arg could be either a size or a box size
            if isinstance(args[0], (int, Integer)):
                if symmetry is None:
                    return PlanePartitions_n(args[0])
                else:
                    raise ValueError("the number of boxes may only be specified if no symmetry is required")
            box_size = args[0]

        box_size = tuple(box_size)
        if symmetry is None:
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

        raise ValueError("invalid symmetry class option")

    def __init__(self, box_size=None, symmetry=None, category=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: PP = PlanePartitions(box_size=[2,2,1])
            sage: TestSuite(PP).run()                                                   # needs sage.modules
        """
        if box_size is not None and len(box_size) != 3:
            raise ValueError("invalid box size")
        self._box = box_size
        self._symmetry = symmetry
        Parent.__init__(self, category=category)

    Element = PlanePartition

    def __contains__(self, pp):
        """
        Check to see that ``pp`` is a valid plane partition.

        EXAMPLES::

            sage: [[3,2,1],[2,1]] in PlanePartitions()
            True
            sage: [[3,2,1],[1,2]] in PlanePartitions()
            False
            sage: [[3,2,1],[3,3]] in PlanePartitions()
            False
        """
        if isinstance(pp, PlanePartition):
            return True
        if isinstance(pp, (list, tuple)):
            if not pp:
                return True
            if not all(a in ZZ for b in pp for a in b):
                return False
            for row in pp:
                if not all(c >= 0 for c in row):
                    return False
                if not all(row[i] >= row[i+1] for i in range(len(row)-1)):
                    return False
            for row, nxt in zip(pp, pp[1:]):
                if not all(row[c] >= nxt[c] for c in range(len(nxt))):
                    return False
            return True
        return False

    def box(self) -> tuple:
        """
        Return the size of the box of the plane partition of ``self``
        is contained in.

        EXAMPLES::

            sage: P = PlanePartitions([4,3,5])
            sage: P.box()
            (4, 3, 5)

            sage: PP = PlanePartitions()
            sage: PP.box() is None
            True
        """
        return self._box

    def symmetry(self) -> str:
        """
        Return the symmetry class of ``self``.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,2], symmetry='SPP')
            sage: PP.symmetry()
            'SPP'
            sage: PP = PlanePartitions()
            sage: PP.symmetry() is None
            True
        """
        return self._symmetry


class PlanePartitions_all(PlanePartitions, DisjointUnionEnumeratedSets):
    r"""
    All plane partitions.
    """
    def __init__(self):
        r"""
        Initialize the class of all plane partitions.

        .. WARNING::

            Input is not checked; please use :class:`PlanePartitions` to
            ensure the options are properly parsed.

        TESTS::

            sage: from sage.combinat.plane_partition import PlanePartitions_all
            sage: P = PlanePartitions_all()
            sage: TestSuite(P).run()
        """
        # We manually set these here rather than invoking the super().__init__().
        # This is so DisjointUnionEnumeratedSets can make the Parent.__init__() call.
        self._box = None
        self._symmetry = None
        # super(PlanePartitions_all, self).__init__(category=InfiniteEnumeratedSets())

        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(NonNegativeIntegers(),
                                                    PlanePartitions_n),
                                             facade=True,
                                             keepkey=False)

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PlanePartitions()
            Plane partitions
        """
        return "Plane partitions"

    def an_element(self):
        r"""
        Return a particular element of the class.

        TESTS::

            sage: P = PlanePartitions()
            sage: P.an_element()
            Plane partition [[2, 1], [1]]
        """
        return self.element_class(self, [[2, 1], [1]])


class PlanePartitions_box(PlanePartitions):
    r"""
    All plane partitions that fit inside a box of a specified size.

    By convention, a plane partition in an `a \times b \times c` box
    will have at most `a` rows, of lengths at most `b`, with entries
    at most `c`.
    """
    def __init__(self, box_size):
        r"""
        Initialize the class of plane partitions that fit in a box of a
        specified size.

        EXAMPLES::

            sage: PP = PlanePartitions([4,3,2])
            sage: TestSuite(PP).run()           # long time                             # needs sage.modules
        """
        super().__init__(box_size, category=FiniteEnumeratedSets())

    def _repr_(self) -> str:
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PlanePartitions([4,3,2])
            Plane partitions inside a 4 x 3 x 2 box
        """
        return "Plane partitions inside a {} x {} x {} box".format(
            self._box[0], self._box[1], self._box[2])

    def __contains__(self, x):
        """
        TESTS::

            sage: [[2,1],[1]] in PlanePartitions([2,2,2])
            True
            sage: [[3,1],[1]] in PlanePartitions([2,2,2])
            False
            sage: [[2,1],[1],[1]] in PlanePartitions([2,2,2])
            False
            sage: [[2,1,1],[1]] in PlanePartitions([2,2,2])
            False
        """
        if len(x) == 0:
            return True
        return PlanePartitions.__contains__(self, x) and len(x) <= self._box[0] and len(x[0]) <= self._box[1] and x[0][0] <= self._box[2]

    def to_poset(self):
        r"""
        Return the product of three chains poset, whose order ideals are
        naturally in bijection with plane partitions inside a box.

        EXAMPLES::

            sage: PlanePartitions([2,2,2]).to_poset()                                   # needs sage.graphs sage.modules
            Finite lattice containing 8 elements
        """
        a = self._box[0]
        b = self._box[1]
        c = self._box[2]
        from sage.combinat.posets.poset_examples import posets
        return posets.ProductOfChains([a, b, c])

    def from_order_ideal(self, I) -> PP:
        r"""
        Return the plane partition corresponding to an order ideal in the
        poset given in :meth:`to_poset`.

        EXAMPLES::

            sage: I = [(1, 0, 0), (1, 0, 1), (1, 1, 0), (0, 1, 0),
            ....:      (0, 0, 0), (0, 0, 1), (0, 1, 1)]
            sage: PlanePartitions([2,2,2]).from_order_ideal(I)                          # needs sage.graphs sage.modules
            Plane partition [[2, 2], [2, 1]]
        """
        return self.from_antichain(self.to_poset().order_ideal_generators(I))

    def from_antichain(self, A) -> PP:
        r"""
        Return the plane partition corresponding to an antichain in the poset
        given in :meth:`to_poset`.

        EXAMPLES::

            sage: A = [(1,0,1), (0,1,1), (1,1,0)]
            sage: PlanePartitions([2,2,2]).from_antichain(A)
            Plane partition [[2, 2], [2, 1]]
        """
        a = self._box[0]
        b = self._box[1]
        # Creates a matrix for the plane partition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]
        pp_matrix = [[0] * (b) for i in range(a)]

        # ac format ex: [x,y,z]
        # iterate through each antichain, assigning the y,z position in pp_matrix = the height of the stack (x + 1)
        for ac in A:
            x = ac[0]
            y = ac[1]
            z = ac[2]
            pp_matrix[x][y] = z + 1

        # For each value in current antichain, fill in the rest of the matrix by
        # rule M[y,z] = Max(M[y+1,z], M[y,z+1]) antichain is now in plane partition format
        if A:
            for i in range(a):
                i = a - (i + 1)
                for j in range(b):
                    j = b - (j + 1)
                    if pp_matrix[i][j] == 0:
                        iValue = 0
                        jValue = 0
                        if i < a-1:
                            iValue = pp_matrix[i+1][j]
                        if j < b-1:
                            jValue = pp_matrix[i][j+1]
                        pp_matrix[i][j] = max(iValue, jValue)
        return self.element_class(self, pp_matrix)

    def __iter__(self) -> Iterator:
        r"""
        Iterate over all partitions that fit inside a box.

        EXAMPLES::

            sage: list(PlanePartitions([1,2,1]))                                        # needs sage.modules
            [Plane partition [], Plane partition [[1]], Plane partition [[1, 1]]]

        TESTS::

            sage: all(len(set(PP)) == PP.cardinality()
            ....:     for b in cartesian_product([range(4)]*3)
            ....:     if (PP := PlanePartitions(b)))
            True
        """
        A = self._box[0]
        B = self._box[1]
        C = self._box[2]
        if not A:
            yield self.element_class(self, [], check=False)
            return
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

            sage: P = PlanePartitions([4,3,5])
            sage: P.cardinality()
            116424
        """
        A = self._box[0]
        B = self._box[1]
        C = self._box[2]
        return Integer(prod(i + j + k - 1
                            for i in range(1, A + 1)
                            for j in range(1, B + 1)
                            for k in range(1, C + 1)) //
                       prod(i + j + k - 2
                            for i in range(1, A + 1)
                            for j in range(1, B + 1)
                            for k in range(1, C + 1)))

    def random_element(self) -> PP:
        r"""
        Return a uniformly random plane partition inside a box.

        ALGORITHM:

        This uses the
        :meth:`~sage.combinat.posets.posets.FinitePoset.random_order_ideal`
        method and the natural bijection with plane partitions.

        EXAMPLES::

            sage: P = PlanePartitions([4,3,5])
            sage: P.random_element()  # random                                          # needs sage.graphs sage.modules
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
        Initialize the class of plane partitions with ``n`` boxes.

        .. WARNING::

            Input is not checked; please use :class:`PlanePartitions` to
            ensure the options are properly parsed.

        TESTS::

            sage: PP = PlanePartitions(4)
            sage: type(PP)
            <class 'sage.combinat.plane_partition.PlanePartitions_n_with_category'>
            sage: TestSuite(PP).run()
        """
        super().__init__(category=FiniteEnumeratedSets())
        self._n = n

    def _repr_(self) -> str:
        """
        TESTS::

            sage: PlanePartitions(3)
            Plane partitions of size 3
        """
        return "Plane partitions of size {}".format(self._n)

    def __contains__(self, x) -> bool:
        """
        TESTS::

            sage: [[2,1],[1]] in PlanePartitions(4)
            True
            sage: [[2,1],[1]] in PlanePartitions(3)
            False
        """
        return PlanePartitions.__contains__(self, x) and PlanePartition(x).number_of_boxes() == self._n

    def __iter__(self) -> Iterator:
        r"""
        Iterate over all plane partitions of a fixed size.

        EXAMPLES::

            sage: list(PlanePartitions(2))
            [Plane partition [[2]], Plane partition [[1, 1]], Plane partition [[1], [1]]]

        TESTS::

            sage: all(len(set(PP)) == PP.cardinality() for n in range(9) if (PP := PlanePartitions(n)))
            True
        """
        from sage.combinat.partition import Partitions

        def P_in_shape_iter(n, la):
            if n < 0 or sum(la) < n:
                return
            if n == 0:
                yield []
                return
            if len(la) == 1:
                if la[0] >= n:
                    yield [n]
                return
            if sum(la) == n:
                yield la
                return
            for mu_0 in range(min(n, la[0]), 0, -1):
                new_la = [min(mu_0, la[i]) for i in range(1, len(la))]
                for mu in P_in_shape_iter(n-mu_0, new_la):
                    yield [mu_0] + mu

        def PP_first_row_iter(n, la):
            m = n - sum(la)
            if m < 0:
                return
            if m == 0:
                yield [la]
                return
            for k in range(m, 0, -1):
                for mu in P_in_shape_iter(k, la):
                    for PP in PP_first_row_iter(m, mu):
                        yield [la] + PP

        n = self._n
        if not n:
            yield PlanePartition([])
            return

        for m in range(n, 0, -1):
            for la in Partitions(m):
                for a in PP_first_row_iter(n, la):
                    yield self.element_class(self, a, check=False)

    def cardinality(self) -> Integer:
        r"""
        Return the number of plane partitions with ``n`` boxes.

        Calculated using the recurrence relation

        .. MATH::

            PL(n) = \sum_{k=1}^n PL(n-k) \sigma_2(k),

        where `\sigma_k(n)` is the sum of the `k`-th powers of
        divisors of `n`.

        EXAMPLES::

            sage: P = PlanePartitions(17)
            sage: P.cardinality()
            18334
        """
        PPn = [1]
        for i in range(1, 1+self._n):
            nextPPn = sum(PPn[i-k] * Sigma()(k, 2) for k in range(1, i+1)) / i
            PPn.append(nextPPn)
        return Integer(PPn[-1])


# Symmetry classes are enumerated and labelled in order as in Proofs and
# Confirmations/Stanley (with all plane partitions being the first class)

# Class 2
# Symmetric Plane Partitions
class PlanePartitions_SPP(PlanePartitions):
    r"""
    Plane partitions that fit inside a box of a specified size that are
    symmetric.
    """
    def __init__(self, box_size):
        """
        Initialize ``self``.

        TESTS::

            sage: PP = PlanePartitions([3,3,2], symmetry='SPP')
            sage: TestSuite(PP).run()                                                   # needs sage.graphs sage.modules
            sage: PlanePartitions([4,3,2], symmetry='SPP')
            Traceback (most recent call last):
            ...
            ValueError: x and y dimensions (4 and 3) must be equal
        """
        if box_size[0] != box_size[1]:
            raise ValueError("x and y dimensions ({} and {}) must be equal".format(box_size[0], box_size[1]))
        super().__init__(box_size, "SPP", category=FiniteEnumeratedSets())

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: PlanePartitions([3,3,2], symmetry='SPP')
            Symmetric plane partitions inside a 3 x 3 x 2 box
        """
        return "Symmetric plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __contains__(self, x) -> bool:
        """
        TESTS::

            sage: [[2,1],[1]] in PlanePartitions([2,2,2], symmetry='SPP')
            True
            sage: [[2,1],[1]] in PlanePartitions([1,1,1], symmetry='SPP')
            False
            sage: [[2,1],[2]] in PlanePartitions([2,2,2], symmetry='SPP')
            False
        """
        P = PlanePartition(x)
        max = (P._max_x, P._max_y, P._max_z)
        return (PlanePartitions.__contains__(self, x)
                and P.is_SPP()
                and all(a <= b for a, b in zip(max, self._box)))

    def to_poset(self):
        r"""
        Return a poset whose order ideals are in bijection with
        symmetric plane partitions.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,2], symmetry='SPP')
            sage: PP.to_poset()                                                         # needs sage.graphs
            Finite poset containing 12 elements
            sage: PP.to_poset().order_ideals_lattice().cardinality() == PP.cardinality()            # needs sage.graphs sage.modules sage.rings.finite_rings
            True
        """
        a = self._box[0]
        c = self._box[2]

        def comp(x, y):
            return all(a <= b for a, b in zip(x, y))

        pl = [(x, y, z) for x in range(a) for y in range(x + 1)
              for z in range(c)]
        from sage.combinat.posets.posets import Poset
        return Poset((pl, comp))

    def from_order_ideal(self, I) -> PP:
        r"""
        Return the symmetric plane partition corresponding to an order ideal
        in the poset given in :meth:`to_poset()`.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,2], symmetry='SPP')
            sage: I = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (2, 0, 0)]
            sage: PP.from_order_ideal(I)                                                # needs sage.graphs
            Plane partition [[1, 1, 1], [1, 1], [1]]
        """
        return self.from_antichain(self.to_poset().order_ideal_generators(I))

    def from_antichain(self, A) -> PP:
        r"""
        Return the symmetric plane partition corresponding to an antichain
        in the poset given in :meth:`to_poset()`.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,2], symmetry='SPP')
            sage: A = [(2, 2, 0), (1, 0, 1), (1, 1, 0)]
            sage: PP.from_antichain(A)
            Plane partition [[2, 2, 1], [2, 1, 1], [1, 1, 1]]
        """
        # Initialize an empty plane partition
        a = self._box[0]
        b = self._box[1]
        pp_matrix = [[0] * (b) for i in range(a)]
        # Antichain indicates where the 'corners' will be in the plane partition
        for ac in A:
            x = ac[0]
            y = ac[1]
            z = ac[2]
            pp_matrix[x][y] = z + 1
        # Fill out the rest of the plane partition using symmetry and the
        # rule pp[i][j]=max(pp[i][j+1],pp[i+1][j])
        if A:
            for i in range(a):
                i = a - (i + 1)
                for j in range(b):
                    j = b - (j + 1)
                    if pp_matrix[i][j] == 0 and i >= j:
                        iValue = 0
                        jValue = 0
                        if i < a - 1:
                            iValue = pp_matrix[i+1][j]
                        if j < b - 1:
                            jValue = pp_matrix[i][j+1]
                        pp_matrix[i][j] = max(iValue, jValue)
                    elif j > i:
                        pp_matrix[i][j] = pp_matrix[j][i]
        return self.element_class(self, pp_matrix)

    def __iter__(self) -> Iterator:
        """
        Iterate over all symmetric plane partitions.

        EXAMPLES::

            sage: list(PlanePartitions([2,2,1], symmetry='SPP'))                        # needs sage.graphs sage.modules sage.rings.finite_rings
            [Plane partition [],
             Plane partition [[1, 1], [1, 1]],
             Plane partition [[1, 1], [1]],
             Plane partition [[1]]]

        TESTS::

            sage: all(len(set(PP)) == PP.cardinality()                                  # needs sage.graphs sage.modules
            ....:     for a, b in cartesian_product([range(4)]*2)
            ....:     if (PP := PlanePartitions([a, a, b], symmetry='SPP')))
            True
        """
        for acl in self.to_poset().antichains_iterator():
            yield self.from_antichain(acl)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of symmetric plane partitions inside an `a \times a \times b`
        box is equal to

        .. MATH::

            \left(\prod_{i=1}^{a} \frac{2i + b - 1}{2i - 1}\right)
            \left(\prod_{1 \leq i < j \leq a} \frac{i+j+b-1}{i+j-1}\right).

        EXAMPLES::

            sage: P = PlanePartitions([3,3,2], symmetry='SPP')
            sage: P.cardinality()
            35
        """
        a = self._box[0]
        c = self._box[2]
        left_prod_num = prod(2*i + c - 1 for i in range(1, a+1))
        left_prod_den = prod(2*i - 1 for i in range(1, a+1))
        right_prod_num = prod(i + j + c - 1
                              for j in range(1, a+1)
                              for i in range(1, j))
        right_prod_den = prod(i + j - 1
                              for j in range(1, a+1)
                              for i in range(1, j))
        return Integer(left_prod_num * right_prod_num // left_prod_den // right_prod_den)

    def random_element(self) -> PP:
        r"""
        Return a uniformly random element of ``self``.

        ALGORITHM:

        This uses the
        :meth:`~sage.combinat.posets.posets.FinitePoset.random_order_ideal`
        method and the natural bijection between symmetric plane partitions
        and order ideals in an associated poset.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,2], symmetry='SPP')
            sage: PP.random_element()  # random                                         # needs sage.graphs
            Plane partition [[2, 2, 2], [2, 2], [2]]
        """
        Z = self.from_order_ideal(self.to_poset().random_order_ideal())
        return self.element_class(self, Z, check=False)


# Class 3
# Cyclically Symmetric Plane Partitions
class PlanePartitions_CSPP(PlanePartitions):
    r"""
    Plane partitions that fit inside a box of a specified size that are
    cyclically symmetric.
    """
    def __init__(self, box_size):
        """
        Initialize ``self``.

        TESTS::

            sage: PP = PlanePartitions([3,3,3], symmetry='CSPP')
            sage: TestSuite(PP).run()                                                   # needs sage.graphs sage.modules sage.rings.finite_rings
            sage: PlanePartitions([4,3,2], symmetry='CSPP')
            Traceback (most recent call last):
            ...
            ValueError: x, y, and z dimensions (4,3,2) must all be equal
        """
        if box_size[0] != box_size[1] or box_size[1] != box_size[2]:
            raise ValueError("x, y, and z dimensions ({},{},{}) must all be equal".format(*box_size))
        super().__init__(box_size, "CSPP", category=FiniteEnumeratedSets())

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: PlanePartitions([3,3,3], symmetry='CSPP')
            Cyclically symmetric plane partitions inside a 3 x 3 x 3 box
        """
        return "Cyclically symmetric plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __contains__(self, x) -> bool:
        """
        TESTS::

            sage: [[2,1],[1]] in PlanePartitions([2,2,2], symmetry='CSPP')
            True
            sage: [[2,1],[1]] in PlanePartitions([1,1,1], symmetry='CSPP')
            False
            sage: [[2,1],[2]] in PlanePartitions([2,2,2], symmetry='CSPP')
            False
        """
        P = PlanePartition(x)
        max = (P._max_x, P._max_y, P._max_z)
        return (PlanePartitions.__contains__(self, x)
                and P.is_CSPP()
                and all(a <= b for a, b in zip(max, self._box)))

    def to_poset(self):
        """
        Return a partially ordered set whose order ideals are in bijection with
        cyclically symmetric plane partitions.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,3], symmetry='CSPP')
            sage: PP.to_poset()                                                         # needs sage.graphs
            Finite poset containing 11 elements
            sage: PP.to_poset().order_ideals_lattice().cardinality() == PP.cardinality()            # needs sage.graphs
            True
        """
        a = self._box[0]
        b = self._box[1]
        c = self._box[2]

        def comp(x, y):
            return all(a <= b for a, b in zip(x, y))

        def comp2(x, y):
            return comp(x, y) or comp(x, (y[2], y[0], y[1])) or comp(x, (y[1], y[2], y[0]))

        pl = [(x, y, z) for x in range(a) for y in range(b) for z in range(x, c)
              if y <= z and (x != z or y == x)]
        from sage.combinat.posets.posets import Poset
        return Poset((pl, comp2))

    def from_antichain(self, acl) -> PP:
        r"""
        Return the cyclically symmetric plane partition corresponding to an
        antichain in the poset given in :meth:`to_poset()`.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,3], symmetry='CSPP')
            sage: A = [(0, 2, 2), (1, 1, 1)]
            sage: PP.from_antichain(A)
            Plane partition [[3, 3, 3], [3, 2, 1], [3, 1, 1]]
        """
        b = self._box[1]
        c = self._box[2]
        pp_matrix = [[0] * (c) for i in range(b)]
        # creates a matrix for the plane partition populated by 0s
        # EX: [[0,0,0], [0,0,0], [0,0,0]]
        # ac format ex: [x,y,z]
        for ac in acl:
            x = ac[0]
            y = ac[1]
            z = ac[2]
            pp_matrix[y][z] = (x+1)
            pp_matrix[z][x] = (y+1)
            pp_matrix[x][y] = (z+1)

        # For each value in current antichain, fill in the rest of the
        # matrix by rule M[y,z] = Max(M[y+1,z], M[y,z+1]) antichain is
        # now in plane partition format.
        if acl != []:
            for i in range(b):
                i = b - (i + 1)
                for j in range(c):
                    j = c - (j + 1)
                    if pp_matrix[i][j] == 0:
                        iValue = 0
                        jValue = 0
                        if i < b - 1:
                            iValue = pp_matrix[i+1][j]
                        if j < c - 1:
                            jValue = pp_matrix[i][j+1]
                        pp_matrix[i][j] = max(iValue, jValue)
        return self.element_class(self, pp_matrix)

    def from_order_ideal(self, I) -> PP:
        r"""
        Return the cylically symmetric plane partition corresponding
        to an order ideal in the poset given in :meth:`to_poset`.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,3], symmetry='CSPP')
            sage: I = [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1), (0, 1, 2),
            ....:      (1, 0, 2), (0, 2, 2), (1, 1, 1), (1, 1, 2), (1, 2, 2)]
            sage: PP.from_order_ideal(I)                                                # needs sage.graphs
            Plane partition [[3, 3, 3], [3, 3, 3], [3, 3, 2]]
        """
        return self.from_antichain(self.to_poset().order_ideal_generators(I))

    def random_element(self) -> PP:
        r"""
        Return a uniformly random element of ``self``.

        ALGORITHM:

        This uses the
        :meth:`~sage.combinat.posets.posets.FinitePoset.random_order_ideal`
        method and the natural bijection between cyclically symmetric plane
        partitions and order ideals in an associated poset.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,3], symmetry='CSPP')
            sage: PP.random_element()  # random                                         # needs sage.graphs
            Plane partition [[3, 2, 2], [3, 1], [1, 1]]
        """
        Z = self.from_order_ideal(self.to_poset().random_order_ideal())
        return self.element_class(self, Z, check=False)

    def __iter__(self) -> Iterator:
        """
        Iterate over all cyclically symmetric plane partitions.

        EXAMPLES::

            sage: list(PlanePartitions([2,2,2], symmetry='CSPP'))                       # needs sage.graphs sage.modules
            [Plane partition [],
             Plane partition [[2, 2], [2, 2]],
             Plane partition [[2, 2], [2, 1]],
             Plane partition [[2, 1], [1]],
             Plane partition [[1]]]

        TESTS::

            sage: all(len(set(PP)) == PP.cardinality() for n in range(5)                # needs sage.graphs sage.modules
            ....:     if (PP := PlanePartitions([n]*3, symmetry='CSPP')))
            True
        """
        for acl in self.to_poset().antichains_iterator():
            yield self.from_antichain(acl)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of cyclically symmetric plane partitions inside an
        `a \times a \times a` box is equal to

        .. MATH::

            \left(\prod_{i=1}^{a} \frac{3i - 1}{3i - 2}\right)
            \left(\prod_{1 \leq i < j \leq a} \frac{i+j+a-1}{2i+j-1}\right).

        EXAMPLES::

            sage: P = PlanePartitions([4,4,4], symmetry='CSPP')
            sage: P.cardinality()
            132
        """
        a = self._box[0]
        num = (prod(3*i - 1 for i in range(1, a + 1))
               * prod(i + j + a - 1 for j in range(1, a + 1)
                      for i in range(1, j + 1)))
        den = (prod(3*i - 2 for i in range(1, a + 1))
               * prod(2*i + j - 1 for j in range(1, a + 1)
                      for i in range(1, j + 1)))
        return Integer(num // den)


# Class 4
# Totally Symmetric Plane Partitions
class PlanePartitions_TSPP(PlanePartitions):
    r"""
    Plane partitions that fit inside a box of a specified size that are
    totally symmetric.
    """
    def __init__(self, box_size):
        """
        Initialize ``self``.

        TESTS::

            sage: PP = PlanePartitions([3,3,3], symmetry='TSPP')
            sage: TestSuite(PP).run()                                                   # needs sage.graphs sage.modules
            sage: PlanePartitions([4,3,2], symmetry='TSPP')
            Traceback (most recent call last):
            ...
            ValueError: x, y, and z dimensions (4,3,2) must all be equal
        """
        if box_size[0] != box_size[1] or box_size[1] != box_size[2]:
            raise ValueError("x, y, and z dimensions ({},{},{}) must all be equal".format(box_size[0], box_size[1], box_size[2]))
        super().__init__(box_size, "TSPP", category=FiniteEnumeratedSets())

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: PlanePartitions([3,3,3], symmetry='TSPP')
            Totally symmetric plane partitions inside a 3 x 3 x 3 box
        """
        return "Totally symmetric plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __contains__(self, x) -> bool:
        """
        TESTS::

            sage: [[2,1],[1]] in PlanePartitions([2,2,2], symmetry='TSPP')
            True
            sage: [[2,1],[1]] in PlanePartitions([1,1,1], symmetry='TSPP')
            False
            sage: [[2,1],[2]] in PlanePartitions([2,2,2], symmetry='TSPP')
            False
        """
        P = PlanePartition(x)
        maxval = (P._max_x, P._max_y, P._max_z)
        return (PlanePartitions.__contains__(self, x) and P.is_TSPP()
                and all(a <= b for a, b in zip(maxval, self._box)))

    def to_poset(self):
        r"""
        Return a poset whose order ideals are in bijection with totally
        symmetric plane partitions.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,3], symmetry='TSPP')
            sage: PP.to_poset()                                                         # needs sage.graphs
            Finite poset containing 10 elements
            sage: (PP.to_poset().order_ideals_lattice().cardinality()                   # needs sage.graphs sage.modules sage.rings.finite_rings
            ....:     == PP.cardinality())
            True
        """
        a = self._box[0]
        b = self._box[1]
        c = self._box[2]

        def comp(x, y):
            return all(a <= b for a, b in zip(x, y))

        pl = [(x, y, z) for x in range(a) for y in range(x, b) for z in range(y, c)]
        from sage.combinat.posets.posets import Poset
        return Poset((pl, comp))

    def from_antichain(self, acl) -> PP:
        r"""
        Return the totally symmetric plane partition corresponding to an
        antichain in the poset given in :meth:`to_poset()`.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,3], symmetry='TSPP')
            sage: A = [(0, 0, 2), (0, 1, 1)]
            sage: PP.from_antichain(A)
            Plane partition [[3, 2, 1], [2, 1], [1]]
        """
        b = self._box[1]
        c = self._box[2]
        pp_matrix = [[0] * (c) for i in range(b)]  # creates a matrix for the plane
        # partition populated by 0s EX: [[0,0,0], [0,0,0], [0,0,0]]
        for ac in acl:
            x = ac[0]
            y = ac[1]
            z = ac[2]

            pp_matrix[y][z] = x + 1  # x,y,z
            pp_matrix[z][x] = y + 1  # y,z,x
            pp_matrix[x][y] = z + 1  # z,x,y

            pp_matrix[z][y] = x + 1  # x,z,y
            pp_matrix[x][z] = y + 1  # y,x,z
            pp_matrix[y][x] = z + 1  # z,y,x

        # for each value in current antichain, fill in the rest of the matrix by
        # rule M[y,z] = Max(M[y+1,z], M[y,z+1]) antichain is now in plane partition format
        if acl != []:
            for i in range(b):
                i = b - (i + 1)
                for j in range(c):
                    j = c - (j + 1)
                    if pp_matrix[i][j] == 0:
                        iValue = 0
                        jValue = 0
                        if i < b - 1:
                            iValue = pp_matrix[i+1][j]
                        if j < c - 1:
                            jValue = pp_matrix[i][j+1]
                        pp_matrix[i][j] = max(iValue, jValue)
        return self.element_class(self, pp_matrix)

    def from_order_ideal(self, I) -> PP:
        r"""
        Return the totally symmetric plane partition corresponding
        to an order ideal in the poset given in :meth:`to_poset`.

        EXAMPLES::

            sage: PP = PlanePartitions([3,3,3], symmetry='TSPP')
            sage: I = [(0, 0, 0), (0, 0, 1), (0, 0, 2), (0, 1, 1)]
            sage: PP.from_order_ideal(I)                                                # needs sage.graphs
            Plane partition [[3, 2, 1], [2, 1], [1]]
        """
        return self.from_antichain(self.to_poset().order_ideal_generators(I))

    def __iter__(self) -> Iterator:
        """
        An iterator for totally symmetric plane partitions.

        EXAMPLES::

            sage: list(PlanePartitions([2,2,2], symmetry='TSPP'))                       # needs sage.graphs sage.modules
            [Plane partition [],
             Plane partition [[2, 2], [2, 2]],
             Plane partition [[2, 2], [2, 1]],
             Plane partition [[2, 1], [1]],
             Plane partition [[1]]]

        TESTS::

            sage: all(len(set(PP)) == PP.cardinality() for n in range(5) if (PP := PlanePartitions([n]*3, symmetry='TSPP')))
            True
        """
        for acl in self.to_poset().antichains_iterator():
            yield self.from_antichain(acl)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of totally symmetric plane partitions inside an
        `a \times a \times a` box is equal to

        .. MATH::

            \prod_{1 \leq i \leq j \leq a} \frac{i+j+a-1}{i+2j-2}.

        EXAMPLES::

            sage: P = PlanePartitions([4,4,4], symmetry='TSPP')
            sage: P.cardinality()
            66
        """
        a = self._box[0]
        num = prod(i + j + a - 1 for j in range(1, a + 1) for i in range(1, j + 1))
        den = prod(i + 2*j - 2 for j in range(1, a + 1) for i in range(1, j + 1))
        return Integer(num // den)


# Class 5
# Self-complementary Plane Partitions
class PlanePartitions_SCPP(PlanePartitions):
    r"""
    Plane partitions that fit inside a box of a specified size that are
    self-complementary.
    """
    def __init__(self, box_size):
        """
        Initialize ``self``.

        TESTS::

            sage: PP = PlanePartitions([4,3,2], symmetry='SCPP')
            sage: TestSuite(PP).run()
            sage: PlanePartitions([5,3,1], symmetry='SCPP')
            Traceback (most recent call last):
            ...
            ValueError: dimensions (5,3,1) cannot all be odd
        """
        if (box_size[0] % 2 == 1 and box_size[1] % 2 == 1 and box_size[2] % 2 == 1):
            raise ValueError("dimensions ({},{},{}) cannot all be odd".format(*box_size))
        super().__init__(box_size, "SCPP", category=FiniteEnumeratedSets())

    def __contains__(self, x) -> bool:
        """
        TESTS::

            sage: [[2,1],[1]] in PlanePartitions([2,2,2], symmetry='SCPP')
            True
            sage: [[2,1],[1]] in PlanePartitions([3,2,2], symmetry='SCPP')
            False
            sage: [[2,1],[1]] in PlanePartitions([2,1,1], symmetry='SCPP')
            False
            sage: [[2,1],[2]] in PlanePartitions([2,2,2], symmetry='SCPP')
            False
        """
        # P = PlanePartitions(self._box)(x)
        # max = (P._max_x, P._max_y, P._max_z)
        # return PlanePartitions.__contains__(self, x) and P.is_SCPP() and all( a<=b for a,b in zip(max,self._box))
        return x in PlanePartitions(self._box) and PlanePartitions(self._box)(x).is_SCPP()

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: PlanePartitions([4,3,2], symmetry='SCPP')
            Self-complementary plane partitions inside a 4 x 3 x 2 box
        """
        return "Self-complementary plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __iter__(self) -> Iterator:
        """
        An iterator for self-complementary plane partitions.

        EXAMPLES::

            sage: list(PlanePartitions([3,2,2], symmetry='SCPP'))
            [Plane partition [[1, 1], [1, 1], [1, 1]],
             Plane partition [[2, 1], [1, 1], [1]],
             Plane partition [[2, 2], [1, 1]],
             Plane partition [[2], [2], [2]],
             Plane partition [[2, 1], [2], [1]],
             Plane partition [[2, 2], [2]]]

        TESTS::

            sage: PP = PlanePartitions([3,4,5], symmetry='SCPP')
            sage: len(set(PP)) == PP.cardinality()
            True

            sage: all(len(set(PP)) == PP.cardinality()
            ....:     for b in cartesian_product([range(4)]*3)
            ....:     if is_even(prod(b)) and (PP := PlanePartitions(b, symmetry='SCPP')))
            True
        """
        b = self._box[0]
        a = self._box[1]
        c = self._box[2]

        def Partitions_inside_lambda(la):
            """
            Iterate over all partitions contained in la with the same number
            of parts including 0s.
            """
            from sage.combinat.partition import Partitions
            for k in range(sum(la), -1, -1):
                for mu in Partitions(k, outer=la):
                    yield mu + [0]*(len(la)-len(mu))

        def Partitions_inside_lambda_with_smallest_at_least_k(la, k):
            """
            Iterate over all partitions contained in la with the smallest
            entry at least k.
            """
            for mu in Partitions_inside_lambda([val - k for val in la]):
                yield [mu[i] + k for i in range(len(la))]

        def possible_middle_row_for_b_odd(a, c):
            """
            Iterate over all possible middle row for SCPP inside box(a,b,c)
            when b is odd.
            """
            if a * c % 2 == 1:
                yield
                return
            for mu in Partitions_inside_lambda([c // 2 for i in range(a // 2)]):
                nu = [c - mu[len(mu)-1-i] for i in range(len(mu))]
                if not a % 2:
                    la = nu + mu
                else:
                    la = nu + [c // 2] + mu
                yield la

        def possible_middle_row_for_b_even(a, c):
            """
            Iterate over all possible middle ((b/2)+1)st row for SCPP inside
            box(a,b,c) when b is even.
            """
            for mu in Partitions_inside_lambda([c // 2 for i in range((a+1) // 2)]):
                if not mu:
                    yield []
                    continue
                nu = [c - mu[len(mu)-1-i] for i in range(a // 2)]
                for tau in Partitions_inside_lambda_with_smallest_at_least_k(nu, mu[0]):
                    la = tau + mu
                    yield la

        def PPs_with_first_row_la_and_with_k_rows(la, k):
            "Iterate over PPs with first row la and with k rows in total."
            if k == 0:
                yield []
                return
            if k == 1:
                yield [la]
                return
            for mu in Partitions_inside_lambda(la):
                for PP in PPs_with_first_row_la_and_with_k_rows(mu, k-1):
                    yield [la] + PP

        def complement(PP, c):
            "Return the complement of PP with respect to height c"
            b = len(PP)
            if not b:
                return []
            a = len(PP[0])
            return [[c - PP[b-1-i][a-1-j] for j in range(a)] for i in range(b)]

        if b % 2 == 1:
            # la is the middle row of SCPP
            for la in possible_middle_row_for_b_odd(a, c):
                for PP in PPs_with_first_row_la_and_with_k_rows(la, (b+1) // 2):
                    PP_below = PP[1:]
                    PP_above = complement(PP_below, c)
                    yield self.element_class(self, PP_above + [la] + PP_below)
        else:
            # la is the middle ((a/2)+1)st row of SCPP
            for la in possible_middle_row_for_b_even(a, c):
                for PP in PPs_with_first_row_la_and_with_k_rows(la, b // 2):
                    PP_below = PP
                    PP_above = complement(PP_below, c)
                    yield self.element_class(self, PP_above + PP_below)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of self complementary plane partitions inside a
        `2a \times 2b \times 2c` box is equal to

        .. MATH::

            \left(\prod_{i=1}^{r}\prod_{j=1}^{b}
            \frac{i + j + c - 1}{i + j - 1}\right)^2.

        The number of self complementary plane partitions inside an
        `(2a+1) \times 2b \times 2c` box is equal to

        .. MATH::

            \left(\prod_{i=1}^{a} \prod_{j=1}^{b} \frac{i+j+c-1}{i+j-1} \right)
            \left(\prod_{i=1}^{a+1} \prod_{j=1}^{b} \frac{i+j+c-1}{i+j-1} \right).

        The number of self complementary plane partitions inside an
        `(2a+1) \times (2b+1) \times 2c` box is equal to

        .. MATH::

            \left(\prod_{i=1}^{a+1} \prod_{j=1}^{b} \frac{i+j+c-1}{i+j-1} \right)
            \left(\prod_{i=1}^{a} \prod_{j=1}^{b+1} \frac{i+j+c-1}{i+j-1} \right).

        EXAMPLES::

            sage: P = PlanePartitions([4,4,4], symmetry='SCPP')
            sage: P.cardinality()
            400

            sage: P = PlanePartitions([5,4,4], symmetry='SCPP')
            sage: P.cardinality()
            1000
            sage: P = PlanePartitions([4,5,4], symmetry='SCPP')
            sage: P.cardinality()
            1000
            sage: P = PlanePartitions([4,4,5], symmetry='SCPP')
            sage: P.cardinality()
            1000

            sage: P = PlanePartitions([5,5,4], symmetry='SCPP')
            sage: P.cardinality()
            2500
            sage: P = PlanePartitions([5,4,5], symmetry='SCPP')
            sage: P.cardinality()
            2500
            sage: P = PlanePartitions([4,5,5], symmetry='SCPP')
            sage: P.cardinality()
            2500
        """
        r = self._box[0]
        s = self._box[1]
        t = self._box[2]
        if r % 2 == 0:
            R = r // 2
            if s % 2 == 0:
                S = s // 2
                if t % 2 == 0:
                    T = t // 2
                    return Integer(prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                        for i in range(1, R+1) for j in range(1, S+1) for k in range(1, T+1))
                                   * prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                          for i in range(1, R+1) for j in range(1, S+1) for k in range(1, T+1)))
                else:
                    T = (t-1) // 2
                    return Integer(prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                        for i in range(1, R+1) for j in range(1, S+1) for k in range(1, T+1))
                                   * prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                          for i in range(1, R+1) for j in range(1, S+1) for k in range(1, T+2)))
            else:
                S = (s-1) // 2
                if t % 2 == 0:
                    T = t // 2
                    return Integer(prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                        for i in range(1, R+1) for j in range(1, S+1) for k in range(1, T+1))
                                   * prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                          for i in range(1, R+1) for j in range(1, S+2) for k in range(1, T+1)))
                else:
                    T = (t-1) // 2
                    return Integer(prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                        for i in range(1, R+1) for j in range(1, S+2) for k in range(1, T+1))
                                   * prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                          for i in range(1, R+1) for j in range(1, S+1) for k in range(1, T+2)))
        # r is odd
        R = (r-1) // 2
        if s % 2 == 0:
            S = s // 2
            if t % 2 == 0:
                T = t // 2
                return Integer(prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                    for i in range(1, R+1) for j in range(1, S+1) for k in range(1, T+1))
                               * prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                      for i in range(1, R+2) for j in range(1, S+1) for k in range(1, T+1)))
            else:
                T = (t-1) // 2
                return Integer(prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                    for i in range(1, R+2) for j in range(1, S+1) for k in range(1, T+1))
                               * prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                      for i in range(1, R+1) for j in range(1, S+1) for k in range(1, T+2)))
        # r and s are both odd
        S = (s-1) // 2
        if t % 2 == 0:
            T = t // 2
            return Integer(prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                for i in range(1, R+2) for j in range(1, S+1) for k in range(1, T+1))
                           * prod(Integer(i+j+k-1) / Integer(i+j+k-2)
                                  for i in range(1, R+1) for j in range(1, S+2) for k in range(1, T+1)))

        # Should never reach here as r, s, t are all odd, which the constructor should reject
        return Integer(0)


# Class 6
# Transpose-complement Plane Partitions
class PlanePartitions_TCPP(PlanePartitions):
    r"""
    Plane partitions that fit inside a box of a specified size that are
    transpose-complement.
    """
    def __init__(self, box_size):
        """
        Initialize ``self``.

        TESTS::

            sage: PP = PlanePartitions([3,3,2], symmetry='TCPP')
            sage: TestSuite(PP).run()                                                   # needs sage.graphs sage.modules

            sage: PlanePartitions([3,3,3], symmetry='TCPP')
            Traceback (most recent call last):
            ...
            ValueError: z dimension (3) must be even

            sage: PlanePartitions([4,3,2], symmetry='TCPP')
            Traceback (most recent call last):
            ...
            ValueError: x and y dimensions (4 and 3) must be equal
        """
        if box_size[2] % 2 == 1:
            raise ValueError("z dimension ({}) must be even".format(box_size[2]))
        if box_size[0] != box_size[1]:
            raise ValueError("x and y dimensions ({} and {}) must be equal".format(box_size[0], box_size[1]))
        super().__init__(box_size, "TCPP", category=FiniteEnumeratedSets())

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: PlanePartitions([3,3,2], symmetry='TCPP')
            Transpose complement plane partitions inside a 3 x 3 x 2 box
        """
        return "Transpose complement plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __iter__(self) -> Iterator:
        r"""
        Iterate over all transpose complement plane partitions.

        EXAMPLES::

            sage: list(PlanePartitions([3,3,2], symmetry='TCPP'))                       # needs sage.modules
            [Plane partition [[2, 2, 1], [2, 1], [1]],
            Plane partition [[2, 1, 1], [2, 1, 1], [1]],
            Plane partition [[2, 2, 1], [1, 1], [1, 1]],
            Plane partition [[2, 1, 1], [1, 1, 1], [1, 1]],
            Plane partition [[1, 1, 1], [1, 1, 1], [1, 1, 1]]]
        """
        for p in PlanePartitions(self._box):
            if p.is_TCPP():
                yield self.element_class(self, p)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of transpose complement plane partitions inside an
        `a \times a \times 2b` box is equal to

        .. MATH::

            \binom{b+1-1}{a-1} \prod_{1\leq i,j \leq a-2}
            \frac{i + j + 2b - 1}{i + j - 1}.

        EXAMPLES::

            sage: P = PlanePartitions([3,3,2], symmetry='TCPP')
            sage: P.cardinality()
            5
        """
        a = self._box[0]
        c = self._box[2]
        return Integer(binomial(c // 2 + a - 1, a - 1)
                       * prod(c + i + j + 1
                              for j in range(1, a - 1) for i in range(1, 1 + j))
                       // prod(i + j + 1
                               for j in range(1, a - 1) for i in range(1, 1 + j)))


# Class 7
# Symmetric Self-complementary Plane Partitions
class PlanePartitions_SSCPP(PlanePartitions):
    r"""
    Plane partitions that fit inside a box of a specified size that are
    symmetric self-complementary.
    """
    def __init__(self, box_size):
        """
        Initialize ``self``.

        TESTS::

            sage: PP = PlanePartitions([2, 2, 4], symmetry='SSCPP')
            sage: TestSuite(PP).run()                                                   # needs sage.modules

            sage: PP = PlanePartitions([4, 4, 2], symmetry='SSCPP')
            sage: TestSuite(PP).run()           # long time                             # needs sage.modules

            sage: PlanePartitions([4, 2, 2], symmetry='SSCPP')
            Traceback (most recent call last):
            ...
            ValueError: x and y dimensions (4 and 2) must be equal

            sage: PlanePartitions([4, 4, 3], symmetry='SSCPP')
            Traceback (most recent call last):
            ...
            ValueError: z dimension (3) must be even
        """
        if box_size[0] != box_size[1]:
            raise ValueError("x and y dimensions ({} and {}) must be equal".format(box_size[0], box_size[1]))
        if (box_size[2] % 2 == 1):
            raise ValueError("z dimension ({}) must be even".format(box_size[2]))
        super().__init__(box_size, "SSCPP", category=FiniteEnumeratedSets())

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: PlanePartitions([4, 4, 2], symmetry='SSCPP')
            Symmetric self-complementary plane partitions inside a 4 x 4 x 2 box
        """
        return "Symmetric self-complementary plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __iter__(self) -> Iterator:
        """
        Iterate over all symmetric self-complementary plane partitions.

        EXAMPLES::

            sage: list(PlanePartitions([4,4,2], symmetry='SSCPP'))                      # needs sage.modules
            [Plane partition [[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]],
             Plane partition [[2, 2, 2, 1], [2, 1, 1], [2, 1, 1], [1]],
             Plane partition [[2, 2, 1, 1], [2, 2, 1, 1], [1, 1], [1, 1]],
             Plane partition [[2, 2, 2, 1], [2, 2, 1], [2, 1], [1]],
             Plane partition [[2, 2, 1, 1], [2, 1, 1, 1], [1, 1, 1], [1, 1]],
             Plane partition [[2, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1]]]

        TESTS::

            sage: all(len(set(PP)) == PP.cardinality()
            ....:     for a, b in cartesian_product([range(5), range(0, 5, 2)])
            ....:     if (PP := PlanePartitions([a, a, b], symmetry='SSCPP')))
            True
        """
        # any SSCPP is a SPP
        for p in PlanePartitions(self._box, symmetry='SPP'):
            if p.is_SSCPP():
                yield self.element_class(self, p)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of symmetric self-complementary plane partitions inside a
        `2a \times 2a \times 2b` box is equal to

        .. MATH::

            \prod_{i=1}^a \prod_{j=1}^a \frac{i + j + b - 1}{i + j - 1}.

        The number of symmetric self-complementary plane partitions inside a
        `(2a+1) \times (2a+1) \times 2b` box is equal to

        .. MATH::

            \prod_{i=1}^a \prod_{j=1}^{a+1} \frac{i + j + b - 1}{i + j - 1}.

        EXAMPLES::

            sage: P = PlanePartitions([4,4,2], symmetry='SSCPP')
            sage: P.cardinality()
            6
            sage: Q = PlanePartitions([3,3,2], symmetry='SSCPP')
            sage: Q.cardinality()
            3
        """
        a = self._box[0]
        c = self._box[2]
        num = prod(i + j + k - 1
                   for i in range(1, 1 + a // 2)
                   for j in range(1, 1 + (a + 1) // 2)
                   for k in range(1, 1 + c // 2))
        den = prod(i + j + k - 2
                   for i in range(1, 1 + a // 2)
                   for j in range(1, 1 + (a + 1) // 2)
                   for k in range(1, 1 + c // 2))
        return Integer(num // den)


# Class 8
# Cyclically Symmetric Transpose-complement Partitions
class PlanePartitions_CSTCPP(PlanePartitions):
    r"""
    Plane partitions that fit inside a box of a specified size that are
    cyclically symmetric and transpose-complement.
    """
    def __init__(self, box_size):
        """
        TESTS::

            sage: PP = PlanePartitions([2,2,2], symmetry='CSTCPP')
            sage: TestSuite(PP).run()                                                   # needs sage.modules

            sage: PlanePartitions([4,3,2], symmetry='CSTCPP')
            Traceback (most recent call last):
            ...
            ValueError: x, y, and z dimensions (4,3,2) must all be equal

            sage: PlanePartitions([3,3,3], symmetry='CSTCPP')
            Traceback (most recent call last):
            ...
            ValueError: x, y, and z dimensions (3,3,3) must all be even
        """
        if box_size[0] != box_size[1] or box_size[1] != box_size[2]:
            raise ValueError("x, y, and z dimensions ({},{},{}) must all be equal".format(*box_size))
        if box_size[0] % 2 == 1:
            raise ValueError("x, y, and z dimensions ({},{},{}) must all be even".format(*box_size))
        super().__init__(box_size, "CSTPP", category=FiniteEnumeratedSets())

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: PlanePartitions([4,4,4], symmetry='CSTCPP')
            Cyclically symmetric transpose complement plane partitions inside a 4 x 4 x 4 box
        """
        return "Cyclically symmetric transpose complement plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __iter__(self) -> Iterator:
        """
        Iterate over all cyclically symmetry transpose complement plane partitions.

        EXAMPLES::

            sage: list(PlanePartitions([2,2,2], symmetry='CSTCPP'))                     # needs sage.modules
            [Plane partition [[2, 1], [1]]]

        TESTS::

            sage: all(len(set(PP)) == PP.cardinality()
            ....:     for n in range(0, 5, 2)
            ....:     if (PP := PlanePartitions([n]*3, symmetry='CSTCPP')))
            True
        """
        # any CSTCPP is a TSPP, a SSCPP and a CSSCPP
        for p in PlanePartitions(self._box, symmetry='TSPP'):
            if p.is_CSTCPP():
                yield self.element_class(self, p)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of cyclically symmetric transpose complement plane partitions
        inside a `2a \times 2a \times 2a` box is equal to

        .. MATH::

            \prod_{i=0}^{a-1} \frac{(3i+1)(6i)!(2i)!}{(4i+1)!(4i)!}.

        EXAMPLES::

            sage: P = PlanePartitions([6,6,6], symmetry='CSTCPP')
            sage: P.cardinality()
            11
        """
        a = self._box[0] // 2
        num = prod((3*i + 1) * factorial(6*i) * factorial(2*i) for i in range(a))
        den = prod((factorial(4*i + 1) * factorial(4*i)) for i in range(a))
        return Integer(num // den)


# Class 9
# Cyclically Symmetric Self-complementary Plane Partitions
class PlanePartitions_CSSCPP(PlanePartitions):
    r"""
    Plane partitions that fit inside a box of a specified size that are
    cyclically symmetric self-complementary.
    """
    def __init__(self, box_size):
        r"""
        Initialize ``self``.

        TESTS::

            sage: PP = PlanePartitions([2,2,2], symmetry='CSSCPP')
            sage: TestSuite(PP).run()                                                   # needs sage.modules
            sage: PlanePartitions([4,3,2], symmetry='CSSCPP')
            Traceback (most recent call last):
            ...
            ValueError: x, y, and z dimensions (4,3,2) must all be equal
            sage: PlanePartitions([3,3,3], symmetry='CSSCPP')
            Traceback (most recent call last):
            ...
            ValueError: x, y, and z dimensions (3,3,3) must all be even
        """
        if box_size[0] != box_size[1] or box_size[1] != box_size[2]:
            raise ValueError("x, y, and z dimensions ({},{},{}) must all be equal".format(*box_size))
        if box_size[0] % 2 == 1:
            raise ValueError("x, y, and z dimensions ({},{},{}) must all be even".format(*box_size))
        super().__init__(box_size, "CSSCPP", category=FiniteEnumeratedSets())

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: PlanePartitions([4,4,4], symmetry='CSSCPP')
            Cyclically symmetric self-complementary plane partitions inside a 4 x 4 x 4 box
        """
        return "Cyclically symmetric self-complementary plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def __iter__(self) -> Iterator:
        """
        Iterate over all cyclically symmetric self-complementary plane partitions.

        EXAMPLES::

            sage: list(PlanePartitions([2,2,2], symmetry='CSSCPP'))                     # needs sage.modules
            [Plane partition [[2, 1], [1]]]
        """
        # any CSSCPP is a SCPP and an CSPP, there are much fewer CSPP
        for p in PlanePartitions(self._box, symmetry='CSPP'):
            if p.is_CSSCPP():
                yield self.element_class(self, p)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of cyclically symmetric self-complementary plane partitions
        inside a `2a \times 2a \times 2a` box is equal to

        .. MATH::

            \left( \prod_{i=0}^{a-1} \frac{(3i+1)!}{(a+i)!} \right)^2.

        EXAMPLES::

            sage: P = PlanePartitions([6,6,6], symmetry='CSSCPP')
            sage: P.cardinality()
            49
        """
        a = self._box[0] // 2
        num = prod(factorial(3*i + 1)**2 for i in range(a))
        den = prod(factorial(a + i)**2 for i in range(a))
        return Integer(num // den)


# Class 10
# Totally Symmetric Self-complementary Plane Partitions
class PlanePartitions_TSSCPP(PlanePartitions):
    r"""
    Plane partitions that fit inside a box of a specified size that are
    totally symmetric self-complementary.
    """
    def __init__(self, box_size):
        """
        TESTS::

            sage: PP = PlanePartitions([4,4,4], symmetry='TSSCPP')
            sage: TestSuite(PP).run()                                                   # needs sage.modules
            sage: PlanePartitions([4,3,2], symmetry='TSSCPP')
            Traceback (most recent call last):
            ...
            ValueError: x, y, and z dimensions (4,3,2) must all be equal
            sage: PlanePartitions([3,3,3], symmetry='TSSCPP')
            Traceback (most recent call last):
            ...
            ValueError: x, y, and z dimensions (3,3,3) must all be even
        """
        if box_size[0] != box_size[1] or box_size[1] != box_size[2]:
            raise ValueError("x, y, and z dimensions ({},{},{}) must all be equal".format(*box_size))
        if box_size[0] % 2 == 1:
            raise ValueError("x, y, and z dimensions ({},{},{}) must all be even".format(*box_size))
        super().__init__(box_size, "TSSCPP", category=FiniteEnumeratedSets())

    def _repr_(self) -> str:
        r"""
        EXAMPLES::

            sage: PlanePartitions([4,4,4], symmetry='TSSCPP')
            Totally symmetric self-complementary plane partitions inside a 4 x 4 x 4 box
        """
        return "Totally symmetric self-complementary plane partitions inside a {} x {} x {} box".format(
                    self._box[0], self._box[1], self._box[2])

    def to_poset(self):
        r"""
        Return a poset whose order ideals are in bijection with
        totally symmetric self-complementary plane partitions.

        EXAMPLES::

            sage: PP = PlanePartitions([6,6,6], symmetry='TSSCPP')
            sage: PP.to_poset()                                                         # needs sage.graphs sage.modules
            Finite poset containing 4 elements
            sage: PP.to_poset().order_ideals_lattice().cardinality() == PP.cardinality()            # needs sage.graphs sage.modules
            True
        """
        from sage.combinat.posets.posets import Poset
        a = self._box[0]
        b = self._box[1]
        c = self._box[2]
        if a != b or b != c or a != c:
            return Poset()

        def comp(x, y):
            return all(xx <= yy for xx, yy in zip(x, y))

        A = a // 2
        pl = [(x, y, z) for x in range(A-1) for y in range(x, A-1)
              for z in range(A-1) if z <= A - 2 - y]
        return Poset((pl, comp))

    def from_antichain(self, acl) -> PP:
        r"""
        Return the totally symmetric self-complementary plane partition
        corresponding to an antichain in the poset given in :meth:`to_poset()`.

        EXAMPLES::

            sage: PP = PlanePartitions([6,6,6], symmetry='TSSCPP')
            sage: A = [(0, 0, 1), (1, 1, 0)]
            sage: PP.from_antichain(A)
            Plane partition [[6, 6, 6, 5, 5, 3], [6, 5, 5, 4, 3, 1], [6, 5, 4, 3, 2, 1],
                             [5, 4, 3, 2, 1], [5, 3, 2, 1, 1], [3, 1, 1]]
        """
        # ac format ex: [x,y,z]
        a = self._box[0]
        b = self._box[1]
        c = self._box[2]
        n = a
        N = n // 2
        pp_matrix = [[0] * (c) for i in range(b)]
        # creates a matrix for the plane partition populated by 0s
        # EX: [[0,0,0], [0,0,0], [0,0,0]]
        width = N - 1
        height = N - 1

        # generate inner triangle
        # FIXME: Make this iterator more efficient
        for i in range(width):
            for j in range(i, height):
                for ac in acl:
                    if ac[0] == i and ac[1] == j:
                        zVal = ac[2]
                        matrixVal = pp_matrix[j+N][i+N]
                        if zVal + 1 > matrixVal:
                            pp_matrix[j+N][i+N] = zVal + 1

        # fill back
        for i in range(width):
            i = width - (i + 1)
            i = i + N
            for j in range(height):
                j = height - (j + 1)
                j = j + N
                if pp_matrix[i][j] == 0:
                    if i >= j:
                        iValue = 0
                        jValue = 0
                        if i < n:
                            iValue = pp_matrix[i+1][j]
                        if j < n:
                            jValue = pp_matrix[i][j+1]
                        pp_matrix[i][j] = max(iValue, jValue)

        # fill half of triangle symmetrically
        for i in range(width):
            i += N
            for j in range(height):
                j += N
                if i >= j:
                    pp_matrix[j][i] = pp_matrix[i][j]

        # upper left box
        for i in range(N):
            for j in range(N):
                pp_matrix[i][j] = n - pp_matrix[n-(i+1)][n-(j+1)]

        # fill in lower left cube with values n/2
        for i in range(N):
            for j in range(N):
                x = i
                y = j
                if pp_matrix[x][y+N] == 0:
                    pp_matrix[x][y+N] = N
                if pp_matrix[x+N][y] == 0:
                    pp_matrix[x+N][y] = N

        # add and subtract values from lower left cube to be rotation of lower right cube
        for i in range(N):
            for j in range(N):
                x = i + N
                y = j + N
                if pp_matrix[x][y] > 0:
                    z = pp_matrix[x][y]
                    for cVal in range(z):
                        # build onto lower left cube
                        pp_matrix[x][0+cVal] += 1
                        # carve out of lower left cube
                        pp_matrix[n-(1+cVal)][N-(j+1)] -= 1

        # fill in upper right cube symmetrically with lower left
        for i in range(N):
            for j in range(N):
                pp_matrix[j][i+N] = pp_matrix[i+N][j]
        return self.element_class(self, pp_matrix)

    def from_order_ideal(self, I) -> PP:
        r"""
        Return the totally symmetric self-complementary plane partition
        corresponding to an order ideal in the poset given in :meth:`to_poset`.

        EXAMPLES::

            sage: PP = PlanePartitions([6,6,6], symmetry='TSSCPP')                      # needs sage.graphs
            sage: I = [(0, 0, 0), (0, 1, 0), (1, 1, 0)]
            sage: PP.from_order_ideal(I)                                                # needs sage.graphs
            Plane partition [[6, 6, 6, 5, 5, 3], [6, 5, 5, 3, 3, 1], [6, 5, 5, 3, 3, 1],
                             [5, 3, 3, 1, 1], [5, 3, 3, 1, 1], [3, 1, 1]]
        """
        return self.from_antichain(self.to_poset().order_ideal_generators(I))

    def __iter__(self) -> Iterator:
        """
        Iterate over all totally symmetric self-complementary plane partitions.

        EXAMPLES::

            sage: list(PlanePartitions([4,4,4], symmetry='TSSCPP'))                     # needs sage.graphs sage.modules
            [Plane partition [[4, 4, 2, 2], [4, 4, 2, 2], [2, 2], [2, 2]],
             Plane partition [[4, 4, 3, 2], [4, 3, 2, 1], [3, 2, 1], [2, 1]]]

        TESTS::

            sage: all(len(set(PP)) == PP.cardinality() for n in range(0,11,2)           # needs sage.graphs sage.modules
            ....:     if (PP := PlanePartitions([n]*3, symmetry='TSSCPP')))
            True
        """
        for acl in self.to_poset().antichains_iterator():
            yield self.from_antichain(acl)

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of ``self``.

        The number of totally symmetric self-complementary plane partitions
        inside a `2a \times 2a \times 2a` box is equal to

        .. MATH::

            \prod_{i=0}^{a-1} \frac{(3i+1)!}{(a+i)!}.

        EXAMPLES::

            sage: P = PlanePartitions([6,6,6], symmetry='TSSCPP')
            sage: P.cardinality()
            7
        """
        a = self._box[0] // 2
        num = prod(factorial(3*i + 1) for i in range(a))
        den = prod(factorial(a + i) for i in range(a))
        return Integer(num // den)
