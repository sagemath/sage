"""
Newton Polygons

This module implements finite Newton polygons and
infinite Newton polygons having a finite number of
slopes (and hence a last infinite slope).
"""

#############################################################################
#    Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
#############################################################################

import sage.geometry.abc

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import op_EQ, op_NE, op_LE, op_GE, op_LT
from sage.misc.cachefunc import cached_method

from sage.rings.infinity import Infinity
from sage.geometry.polyhedron.constructor import Polyhedron


class NewtonPolygon_element(Element):
    """
    Class for infinite Newton polygons with last slope.
    """
    def __init__(self, polyhedron, parent):
        """
        Initialize a Newton polygon.

        INPUT:

        - ``polyhedron`` -- a polyhedron defining the Newton polygon

        TESTS:

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NewtonPolygon([ (0,0), (1,1), (3,5) ])
            Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (3, 5)

            sage: NewtonPolygon([ (0,0), (1,1), (2,8), (3,5) ], last_slope=3)
            Infinite Newton polygon with 3 vertices: (0, 0), (1, 1), (3, 5) ending by an infinite line of slope 3

        ::

            sage: TestSuite(NewtonPolygon).run()
        """
        Element.__init__(self, parent)
        self._polyhedron = polyhedron
        self._vertices = None
        if polyhedron.is_mutable():
            polyhedron._add_dependent_object(self)

    def _repr_(self):
        """
        Return a string representation of this Newton polygon.

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP = NewtonPolygon([ (0,0), (1,1), (2,5) ]); NP
            Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (2, 5)

            sage: NP._repr_()
            'Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (2, 5)'
        """
        vertices = self.vertices()
        length = len(vertices)
        if self.last_slope() is Infinity:
            if length == 0:
                return "Empty Newton polygon"
            elif length == 1:
                return "Finite Newton polygon with 1 vertex: %s" % str(vertices[0])
            else:
                return "Finite Newton polygon with %s vertices: %s" % (length, str(vertices)[1:-1])
        else:
            if length == 1:
                return "Newton Polygon consisting of a unique infinite line of slope %s starting at %s" % (self.last_slope(), str(vertices[0]))
            else:
                return "Infinite Newton polygon with %s vertices: %s ending by an infinite line of slope %s" % (length, str(vertices)[1:-1], self.last_slope())

    def vertices(self, copy=True):
        """
        Return the list of vertices of this Newton polygon.

        INPUT:

        - ``copy`` -- boolean (default: ``True``)

        OUTPUT: the list of vertices of this Newton polygon (or a copy of it
        if ``copy`` is set to True)

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP = NewtonPolygon([ (0,0), (1,1), (2,5) ]); NP
            Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (2, 5)

            sage: v = NP.vertices(); v
            [(0, 0), (1, 1), (2, 5)]

        TESTS:

            sage: del v[0]
            sage: v
            [(1, 1), (2, 5)]
            sage: NP.vertices()
            [(0, 0), (1, 1), (2, 5)]
        """
        if self._vertices is None:
            self._vertices = [ tuple(v) for v in self._polyhedron.vertices() ]
            self._vertices.sort()
        if copy:
            return list(self._vertices)
        else:
            return self._vertices

    @cached_method
    def last_slope(self):
        """
        Return the last (infinite) slope of this Newton polygon
        if it is infinite and ``+Infinity`` otherwise.

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP1 = NewtonPolygon([ (0,0), (1,1), (2,8), (3,5) ], last_slope=3)
            sage: NP1.last_slope()
            3

            sage: NP2 = NewtonPolygon([ (0,0), (1,1), (2,5) ])
            sage: NP2.last_slope()
            +Infinity

        We check that the last slope of a sum (resp. a product) is the
        minimum of the last slopes of the summands (resp. the factors)::

            sage: (NP1 + NP2).last_slope()
            3
            sage: (NP1 * NP2).last_slope()
            3
        """
        rays = self._polyhedron.rays()
        for r in rays:
            if r[0] > 0:
                return r[1]/r[0]
        return Infinity

    def slopes(self, repetition=True):
        """
        Return the slopes of this Newton polygon.

        INPUT:

        - ``repetition`` -- boolean (default: ``True``)

        OUTPUT:

        The consecutive slopes (not including the last slope
        if the polygon is infinity) of this Newton polygon.

        If ``repetition`` is True, each slope is repeated a number of
        times equal to its length. Otherwise, it appears only one time.

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP = NewtonPolygon([ (0,0), (1,1), (3,6) ]); NP
            Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (3, 6)

            sage: NP.slopes()
            [1, 5/2, 5/2]

            sage: NP.slopes(repetition=False)
            [1, 5/2]
        """
        slopes = [ ]
        vertices = self.vertices(copy=False)
        for i in range(1,len(vertices)):
            dx = vertices[i][0] - vertices[i-1][0]
            dy = vertices[i][1] - vertices[i-1][1]
            slope = dy/dx
            if repetition:
                slopes.extend(dx * [slope])
            else:
                slopes.append(slope)
        return slopes

    def _add_(self, other):
        """
        Return the convex hull of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a Newton polygon

        OUTPUT: the Newton polygon, which is the convex hull of this Newton polygon and ``other``

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP1 = NewtonPolygon([ (0,0), (1,1), (2,6) ]); NP1
            Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (2, 6)
            sage: NP2 = NewtonPolygon([ (0,0), (1,3/2) ], last_slope=2); NP2
            Infinite Newton polygon with 2 vertices: (0, 0), (1, 3/2) ending by an infinite line of slope 2

            sage: NP1 + NP2
            Infinite Newton polygon with 2 vertices: (0, 0), (1, 1) ending by an infinite line of slope 2
        """
        polyhedron = self._polyhedron.convex_hull(other._polyhedron)
        return self.parent()(polyhedron)

    def _mul_(self, other):
        """
        Return the Minkowski sum of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a Newton polygon

        OUTPUT: the Newton polygon, which is the Minkowski sum of this Newton polygon and ``other``

        .. NOTE::

            If ``self`` and ``other`` are respective Newton polygons of some polynomials
            `f` and `g` the self*other is the Newton polygon of the product `fg`

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP1 = NewtonPolygon([ (0,0), (1,1), (2,6) ]); NP1
            Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (2, 6)
            sage: NP2 = NewtonPolygon([ (0,0), (1,3/2) ], last_slope=2); NP2
            Infinite Newton polygon with 2 vertices: (0, 0), (1, 3/2) ending by an infinite line of slope 2

            sage: NP = NP1 * NP2; NP
            Infinite Newton polygon with 3 vertices: (0, 0), (1, 1), (2, 5/2) ending by an infinite line of slope 2

        The slopes of ``NP`` is the union of those of ``NP1`` and those of ``NP2``
        which are less than the last slope::

            sage: NP1.slopes()
            [1, 5]
            sage: NP2.slopes()
            [3/2]
            sage: NP.slopes()
            [1, 3/2]
        """
        polyhedron = self._polyhedron.minkowski_sum(other._polyhedron)
        return self.parent()(polyhedron)

    def __pow__(self, exp, ignored=None):
        """
        Return ``self`` dilated by ``exp``.

        INPUT:

        - ``exp`` -- positive integer

        OUTPUT: this Newton polygon scaled by a factor ``exp``

        .. NOTE::

            If ``self`` is the Newton polygon of a polynomial `f`, then
            ``self^exp`` is the Newton polygon of `f^{exp}`.

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP = NewtonPolygon([ (0,0), (1,1), (2,6) ]); NP
            Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (2, 6)

            sage: NP^10
            Finite Newton polygon with 3 vertices: (0, 0), (10, 10), (20, 60)
        """
        polyhedron = self._polyhedron.dilation(exp)
        return self.parent()(polyhedron)

    def __lshift__(self, i):
        """
        Return ``self`` shifted by `(0,i)`.

        INPUT:

        - ``i`` -- a rational number

        OUTPUT: this Newton polygon shifted by the vector `(0,i)`

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP = NewtonPolygon([ (0,0), (1,1), (2,6) ]); NP
            Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (2, 6)

            sage: NP << 2
            Finite Newton polygon with 3 vertices: (0, 2), (1, 3), (2, 8)
        """
        polyhedron = self._polyhedron.translation((0,i))
        return self.parent()(polyhedron)

    def __rshift__(self, i):
        """
        Return ``self`` shifted by `(0,-i)`.

        INPUT:

        - ``i`` -- a rational number

        OUTPUT: this Newton polygon shifted by the vector `(0,-i)`

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP = NewtonPolygon([ (0,0), (1,1), (2,6) ]); NP
            Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (2, 6)

            sage: NP >> 2
            Finite Newton polygon with 3 vertices: (0, -2), (1, -1), (2, 4)
        """
        polyhedron = self._polyhedron.translation((0,-i))
        return self.parent()(polyhedron)

    def __call__(self, x):
        """
        Return `self(x)`.

        INPUT:

        - ``x`` -- a real number

        OUTPUT: the value of this Newton polygon at abscissa `x`

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP = NewtonPolygon([ (0,0), (1,1), (3,6) ]); NP
            Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (3, 6)

            sage: [ NP(i) for i in range(4) ]
            [0, 1, 7/2, 6]
        """
        # complexity: O(log(n))
        vertices = self.vertices()
        lastslope = self.last_slope()
        if len(vertices) == 0 or x < vertices[0][0]:
            return Infinity
        if x == vertices[0][0]:
            return vertices[0][1]
        if x == vertices[-1][0]:
            return vertices[-1][1]
        if x > vertices[-1][0]:
            return vertices[-1][1] + lastslope * (x - vertices[-1][0])
        a = 0
        b = len(vertices)
        while b - a > 1:
            c = (a + b) // 2
            if vertices[c][0] < x:
                a = c
            else:
                b = c
        xg, yg = vertices[a]
        xd, yd = vertices[b]
        return ((x-xg)*yd + (xd-x)*yg) / (xd-xg)

    def _richcmp_(self, other, op):
        r"""
        Comparisons of two Newton polygons.

        TESTS::

            sage: from sage.geometry.newton_polygon import NewtonPolygon

            sage: NP1 = NewtonPolygon([ (0,0), (1,1), (3,6) ])
            sage: NP2 = NewtonPolygon([ (0,0), (1,1), (2,6), (3,6) ])
            sage: NP1 == NP2
            True
            sage: NP1 != NP2
            False

            sage: NP1 >= NP1 and NP2 >= NP2
            True
            sage: NP1 > NP1 or NP2 > NP2
            False

            sage: NP1 = NewtonPolygon([ (0,0), (1,1), (2,6) ])
            sage: NP2 = NewtonPolygon([ (0,0), (1,3/2) ], last_slope=2)
            sage: NP3 = NP1 + NP2

            sage: NP1 <= NP2
            False
            sage: NP3 <= NP1
            True
            sage: NP3 <= NP2
            True

            sage: NP1 < NP1
            False
            sage: NP1 < NP2
            False

            sage: NP1 >= NP2
            False

            sage: NP1 >= NP3
            True

            sage: NP1 > NP1
            False
            sage: NP1 > NP2
            False

            sage: NP1 >= NP3 and NP2 >= NP3 and NP3 <= NP1 and NP3 <= NP2
            True
            sage: NP1 > NP3 and NP2 > NP3
            True
            sage: NP3 < NP2 and NP3 < NP1
            True
        """
        if self._polyhedron == other._polyhedron:
            return op == op_EQ or op == op_LE or op == op_GE
        elif op == op_NE:
            return True
        elif op == op_EQ:
            return False

        if op == op_LT or op == op_LE:
            if self.last_slope() > other.last_slope():
                return False
            return all(v in self._polyhedron for v in other.vertices())

        else:
            if self.last_slope() < other.last_slope():
                return False
            return all(v in other._polyhedron for v in self.vertices())

    def plot(self, **kwargs):
        """
        Plot this Newton polygon.

        .. NOTE::

            All usual rendering options (color, thickness, etc.) are available.

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP = NewtonPolygon([ (0,0), (1,1), (2,6) ])
            sage: polygon = NP.plot()                                                   # needs sage.plot
        """
        vertices = self.vertices()
        if len(vertices) == 0:
            from sage.plot.graphics import Graphics
            return Graphics()
        else:
            from sage.plot.line import line
            (xstart,ystart) = vertices[0]
            (xend,yend) = vertices[-1]
            if self.last_slope() is Infinity:
                return line([(xstart, ystart+1), (xstart,ystart+0.5)], linestyle='--', **kwargs) \
                     + line([(xstart, ystart+0.5)] + vertices + [(xend, yend+0.5)], **kwargs) \
                     + line([(xend, yend+0.5), (xend, yend+1)], linestyle='--', **kwargs)
            else:
                return line([(xstart, ystart+1), (xstart,ystart+0.5)], linestyle='--', **kwargs) \
                     + line([(xstart, ystart+0.5)] + vertices + [(xend+0.5, yend + 0.5*self.last_slope())], **kwargs) \
                     + line([(xend+0.5, yend + 0.5*self.last_slope()), (xend+1, yend+self.last_slope())], linestyle='--', **kwargs)

    def reverse(self, degree=None):
        r"""
        Return the symmetric of ``self``.

        INPUT:

        - ``degree`` -- integer (default: the top right abscissa of
          this Newton polygon)

        OUTPUT:

        The image this Newton polygon under the symmetry
        '(x,y) \mapsto (degree-x, y)`.

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NP = NewtonPolygon([ (0,0), (1,1), (2,5) ])
            sage: NP2 = NP.reverse(); NP2
            Finite Newton polygon with 3 vertices: (0, 5), (1, 1), (2, 0)

        We check that the slopes of the symmetric Newton polygon are
        the opposites of the slopes of the original Newton polygon::

            sage: NP.slopes()
            [1, 4]
            sage: NP2.slopes()
            [-4, -1]
        """
        if self.last_slope() is not Infinity:
            raise ValueError("Can only reverse *finite* Newton polygons")
        if degree is None:
            degree = self.vertices()[-1][0]
        vertices = [ (degree-x,y) for (x,y) in self.vertices() ]
        vertices.reverse()
        parent = self.parent()
        polyhedron = Polyhedron(base_ring=parent.base_ring(), vertices=vertices, rays=[(0,1)])
        return parent(polyhedron)


class ParentNewtonPolygon(Parent, UniqueRepresentation):
    r"""
    Construct a Newton polygon.

    INPUT:

    - ``arg`` -- list/tuple/iterable of vertices or of
      slopes. Currently, slopes must be rational numbers

    - ``sort_slopes`` -- boolean (default: ``True``);  whether slopes must be
      first sorted

    - ``last_slope`` -- rational or infinity (default:
      ``Infinity``); the last slope of the Newton polygon

    OUTPUT: the corresponding Newton polygon

    .. NOTE::

        By convention, a Newton polygon always contains the point
        at infinity `(0, \infty)`. These polygons are attached to
        polynomials or series over discrete valuation rings (e.g. padics).

    EXAMPLES:

    We specify here a Newton polygon by its vertices::

        sage: from sage.geometry.newton_polygon import NewtonPolygon
        sage: NewtonPolygon([ (0,0), (1,1), (3,5) ])
        Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (3, 5)

    We note that the convex hull of the vertices is automatically
    computed::

        sage: NewtonPolygon([ (0,0), (1,1), (2,8), (3,5) ])
        Finite Newton polygon with 3 vertices: (0, 0), (1, 1), (3, 5)

    Note that the value ``+Infinity`` is allowed as the second coordinate
    of a vertex::

        sage: NewtonPolygon([ (0,0), (1,Infinity), (2,8), (3,5) ])
        Finite Newton polygon with 2 vertices: (0, 0), (3, 5)

    If last_slope is set, the returned Newton polygon is infinite
    and ends with an infinite line having the specified slope::

        sage: NewtonPolygon([ (0,0), (1,1), (2,8), (3,5) ], last_slope=3)
        Infinite Newton polygon with 3 vertices: (0, 0), (1, 1), (3, 5) ending by an infinite line of slope 3

    Specifying a last slope may discard some vertices::

        sage: NewtonPolygon([ (0,0), (1,1), (2,8), (3,5) ], last_slope=3/2)
        Infinite Newton polygon with 2 vertices: (0, 0), (1, 1) ending by an infinite line of slope 3/2

    Next, we define a Newton polygon by its slopes::

        sage: NP = NewtonPolygon([0, 1/2, 1/2, 2/3, 2/3, 2/3, 1, 1])
        sage: NP
        Finite Newton polygon with 5 vertices: (0, 0), (1, 0), (3, 1), (6, 3), (8, 5)
        sage: NP.slopes()
        [0, 1/2, 1/2, 2/3, 2/3, 2/3, 1, 1]

    By default, slopes are automatically sorted::

        sage: NP2 = NewtonPolygon([0, 1, 1/2, 2/3, 1/2, 2/3, 1, 2/3])
        sage: NP2
        Finite Newton polygon with 5 vertices: (0, 0), (1, 0), (3, 1), (6, 3), (8, 5)
        sage: NP == NP2
        True

    except if the contrary is explicitly mentioned::

        sage: NewtonPolygon([0, 1, 1/2, 2/3, 1/2, 2/3, 1, 2/3], sort_slopes=False)
        Finite Newton polygon with 4 vertices: (0, 0), (1, 0), (6, 10/3), (8, 5)

    Slopes greater that or equal last_slope (if specified) are discarded::

        sage: NP = NewtonPolygon([0, 1/2, 1/2, 2/3, 2/3, 2/3, 1, 1], last_slope=2/3)
        sage: NP
        Infinite Newton polygon with 3 vertices: (0, 0), (1, 0), (3, 1) ending by an infinite line of slope 2/3
        sage: NP.slopes()
        [0, 1/2, 1/2]

    Be careful, do not confuse Newton polygons provided by this class
    with Newton polytopes. Compare::

        sage: NP = NewtonPolygon([ (0,0), (1,45), (3,6) ]); NP
        Finite Newton polygon with 2 vertices: (0, 0), (3, 6)

        sage: x, y = polygen(QQ,'x, y')
        sage: p = 1 + x*y**45 + x**3*y**6
        sage: p.newton_polytope()
        A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
        sage: p.newton_polytope().vertices()
        (A vertex at (0, 0), A vertex at (1, 45), A vertex at (3, 6))
    """

    Element = NewtonPolygon_element

    def __init__(self):
        """
        Parent class for all Newton polygons.

        EXAMPLES::

            sage: from sage.geometry.newton_polygon import ParentNewtonPolygon
            sage: ParentNewtonPolygon()
            Parent for Newton polygons

        TESTS:

        This class is a singleton::

            sage: ParentNewtonPolygon() is ParentNewtonPolygon()
            True

        ::

            sage: TestSuite(ParentNewtonPolygon()).run()
        """
        from sage.categories.semirings import Semirings
        from sage.rings.rational_field import QQ
        Parent.__init__(self, category=Semirings(), base=QQ)

    def _repr_(self):
        """
        Return the string representation of this parent,
        which is ``Parent for Newton polygons``.

        TESTS:

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NewtonPolygon
            Parent for Newton polygons

            sage: NewtonPolygon._repr_()
            'Parent for Newton polygons'
        """
        return "Parent for Newton polygons"

    def _an_element_(self):
        """
        Return a Newton polygon (which is the empty one).

        TESTS:

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NewtonPolygon._an_element_()
            Empty Newton polygon
        """
        return self(Polyhedron(base_ring=self.base_ring(), ambient_dim=2))

    def _element_constructor_(self, arg, sort_slopes=True, last_slope=Infinity):
        r"""
        INPUT:

        - ``arg`` -- an argument describing the Newton polygon

        - ``sort_slopes`` -- boolean (default: ``True``); whether
          slopes must be first sorted

        - ``last_slope`` -- rational or infinity (default:
          ``Infinity``). The last slope of the Newton polygon

        The first argument ``arg`` can be either:

        - a polyhedron in `\QQ^2`

        - the element ``0`` (corresponding to the empty Newton polygon)

        - the element ``1`` (corresponding to the Newton polygon of the
          constant polynomial equal to 1)

        - a list/tuple/iterable of vertices

        - a list/tuple/iterable of slopes

        OUTPUT: the corresponding Newton polygon

        For more informations, see :class:`ParentNewtonPolygon`.

        TESTS:

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: NewtonPolygon(0)
            Empty Newton polygon
            sage: NewtonPolygon(1)
            Finite Newton polygon with 1 vertex: (0, 0)
        """
        if isinstance(arg, sage.geometry.abc.Polyhedron):
            return self.element_class(arg, parent=self)
        if arg == 0:
            polyhedron = Polyhedron(base_ring=self.base_ring(), ambient_dim=2)
            return self.element_class(polyhedron, parent=self)
        if arg == 1:
            polyhedron = Polyhedron(base_ring=self.base_ring(),
                                    vertices=[(0, 0)], rays=[(0, 1)])
            return self.element_class(polyhedron, parent=self)
        if not isinstance(arg, list):
            try:
                arg = list(arg)
            except TypeError:
                raise TypeError("argument must be a list of coordinates or a list of (rational) slopes")
        if arg and arg[0] in self.base_ring():
            if sort_slopes:
                arg.sort()
            x = y = 0
            vertices = [(x, y)]
            for slope in arg:
                if slope not in self.base_ring():
                    raise TypeError("argument must be a list of coordinates or a list of (rational) slopes")
                x += 1
                y += slope
                vertices.append((x, y))
        else:
            vertices = [(x, y) for (x, y) in arg if y is not Infinity]

        if len(vertices) == 0:
            polyhedron = Polyhedron(base_ring=self.base_ring(), ambient_dim=2)
        else:
            rays = [(0, 1)]
            if last_slope is not Infinity:
                rays.append((1, last_slope))
            polyhedron = Polyhedron(base_ring=self.base_ring(), vertices=vertices, rays=rays)
        return self.element_class(polyhedron, parent=self)


NewtonPolygon = ParentNewtonPolygon()
