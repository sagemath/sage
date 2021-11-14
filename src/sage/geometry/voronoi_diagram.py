r"""
Voronoi diagram

This module provides the class :class:`VoronoiDiagram` for computing the
Voronoi diagram of a finite list of points in `\RR^d`.
"""

#*****************************************************************************
#       Copyright (C) 2012 Moritz Firsching <moritz@math.fu-berlin.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.rings.qqbar import AA
from sage.rings.rational_field import QQ
import sage.rings.abc
from sage.geometry.triangulation.point_configuration import PointConfiguration
from sage.modules.free_module_element import vector
from sage.plot.all import line, point, rainbow, plot


class VoronoiDiagram(SageObject):
    r"""
    Base class for the  Voronoi diagram.

    Compute the Voronoi diagram of a list of points.

    INPUT:

    - ``points`` -- a list of points. Any valid input for the
      :class:`PointConfiguration` will do.

    OUTPUT:

    An instance of the VoronoiDiagram class.

    EXAMPLES:

    Get the Voronoi diagram for some points in `\RR^3`::

        sage: V = VoronoiDiagram([[1, 3, .3], [2, -2, 1], [-1, 2, -.1]]); V
        The Voronoi diagram of 3 points of dimension 3 in the Real Double Field

        sage: VoronoiDiagram([])
        The empty Voronoi diagram.

    Get the Voronoi diagram of a regular pentagon in ``AA^2``.
    All cells meet at the origin::

        sage: DV = VoronoiDiagram([[AA(c) for c in v] for v in polytopes.regular_polygon(5).vertices_list()]); DV  # optional - sage.rings.number_field
        The Voronoi diagram of 5 points of dimension 2 in the Algebraic Real Field
        sage: all(P.contains([0, 0]) for P in DV.regions().values())                                               # optional - sage.rings.number_field
        True
        sage: any(P.interior_contains([0, 0]) for P in DV.regions().values())                                      # optional - sage.rings.number_field
        False

    If the vertices are not converted to ``AA`` before, the method throws an error::

        sage: polytopes.dodecahedron().vertices_list()[0][0].parent()                                              # optional - sage.rings.number_field
        Number Field in sqrt5 with defining polynomial x^2 - 5 with sqrt5 = 2.236067977499790?
        sage: VoronoiDiagram(polytopes.dodecahedron().vertices_list())                                             # optional - sage.rings.number_field
        Traceback (most recent call last):
        ...
        NotImplementedError: Base ring of the Voronoi diagram must be
        one of QQ, RDF, AA.

    ALGORITHM:

    We use hyperplanes tangent to the paraboloid one dimension higher to
    get a convex polyhedron and then project back to one dimension lower.

    .. TODO::

     - The dual construction: Delaunay triangulation
     - improve 2d-plotting
     - implement 3d-plotting
     - more general constructions, like Voroi diagrams with weights (power diagrams)

    REFERENCES:

     - [Mat2002]_ Ch.5.7, p.118.

    AUTHORS:

    - Moritz Firsching (2012-09-21)
    """
    def __init__(self, points):
        r"""
        See ``VoronoiDiagram`` for full documentation.

        EXAMPLES::

            sage: V = VoronoiDiagram([[1, 3, 3], [2, -2, 1], [-1 ,2, -1]]); V
            The Voronoi diagram of 3 points of dimension 3 in the Rational Field
        """
        self._P = {}
        self._points = PointConfiguration(points)
        self._n = self._points.n_points()
        if not self._n or self._points.base_ring().is_subring(QQ):
            self._base_ring = QQ
        elif isinstance(self._points.base_ring(), sage.rings.abc.RealDoubleField) or self._points.base_ring() == AA:
            self._base_ring = self._points.base_ring()
        elif isinstance(self._points.base_ring(), sage.rings.abc.RealField):
            from sage.rings.real_double import RDF
            self._base_ring = RDF
            self._points = PointConfiguration([[RDF(cor) for cor in poi]
                                               for poi in self._points])
        else:
            raise NotImplementedError('Base ring of the Voronoi diagram must '
                                      'be one of QQ, RDF, AA.')

        if self._n > 0:
            self._d = self._points.ambient_dim()
            e = [([sum(vector(i)[k] ** 2
                       for k in range(self._d))] +
                  [(-2) * vector(i)[l] for l in range(self._d)] + [1])
                 for i in self._points]
            # we attach hyperplane to the paraboloid

            e = [[self._base_ring(i) for i in k] for k in e]
            p = Polyhedron(ieqs=e, base_ring=self._base_ring)
            # To understand the reordering that takes place when
            # defining a rational polyhedron, we generate two sorted
            # lists, that are used a few lines below
            if self.base_ring() == QQ:
                enormalized = []
                for ineq in e:
                    if ineq[0] == 0:
                        enormalized.append(ineq)
                    else:
                        enormalized.append([i / ineq[0] for i in ineq[1:]])
                # print enormalized
                hlist = [list(ineq) for ineq in p.Hrepresentation()]
                hlistnormalized = []
                for ineq in hlist:
                    if ineq[0] == 0:
                        hlistnormalized.append(ineq)
                    else:
                        hlistnormalized.append([i / ineq[0] for i in ineq[1:]])
                # print hlistnormalized

        for i in range(self._n):
            # for base ring RDF and AA, Polyhedron keeps the order of the
            # points in the input, for QQ we resort
            if self.base_ring() == QQ:
                equ = p.Hrepresentation(hlistnormalized.index(enormalized[i]))
            else:
                equ = p.Hrepresentation(i)
            pvert = [[u[k] for k in range(self._d)] for u in equ.incident()
                     if u.is_vertex()]
            prays = [[u[k] for k in range(self._d)] for u in equ.incident()
                     if u.is_ray()]
            pline = [[u[k] for k in range(self._d)] for u in equ.incident()
                     if u.is_line()]
            (self._P)[self._points[i]] = Polyhedron(vertices=pvert,
                                                    lines=pline, rays=prays,
                                                    base_ring=self._base_ring)

    def points(self):
        r"""
        Return the input points (as a PointConfiguration).

        EXAMPLES::

            sage: V = VoronoiDiagram([[.5, 3], [2, 5], [4, 5], [4, -1]]); V.points()
            A point configuration in affine 2-space over Real Field
            with 53 bits of precision consisting of 4 points.
            The triangulations of this point configuration are
            assumed to be connected, not necessarily fine,
            not necessarily regular.
        """
        return self._points

    def ambient_dim(self):
        r"""
        Return the ambient dimension of the points.

        EXAMPLES::

            sage: V = VoronoiDiagram([[.5, 3], [2, 5], [4, 5], [4, -1]])
            sage: V.ambient_dim()
            2
            sage: V = VoronoiDiagram([[1, 2, 3, 4, 5, 6]]); V.ambient_dim()
            6
        """
        return self._d

    def regions(self):
        r"""
        Return the Voronoi regions of the Voronoi diagram as a
        dictionary of polyhedra.

        EXAMPLES::

            sage: V = VoronoiDiagram([[1, 3, .3], [2, -2, 1], [-1, 2, -.1]])
            sage: P = V.points()
            sage: V.regions() == {P[0]: Polyhedron(base_ring=RDF, lines=[(-RDF(0.375), RDF(0.13888888890000001), RDF(1.5277777779999999))],
            ....:                                                 rays=[(RDF(9), -RDF(1), -RDF(20)), (RDF(4.5), RDF(1), -RDF(25))],
            ....:                                                 vertices=[(-RDF(1.1074999999999999), RDF(1.149444444), RDF(9.0138888890000004))]),
            ....:                 P[1]: Polyhedron(base_ring=RDF, lines=[(-RDF(0.375), RDF(0.13888888890000001), RDF(1.5277777779999999))],
            ....:                                                 rays=[(RDF(9), -RDF(1), -RDF(20)), (-RDF(2.25), -RDF(1), RDF(2.5))],
            ....:                                                  vertices=[(-RDF(1.1074999999999999), RDF(1.149444444), RDF(9.0138888890000004))]),
            ....:                 P[2]: Polyhedron(base_ring=RDF, lines=[(-RDF(0.375), RDF(0.13888888890000001), RDF(1.5277777779999999))],
            ....:                                                 rays=[(RDF(4.5), RDF(1), -RDF(25)), (-RDF(2.25), -RDF(1), RDF(2.5))],
            ....:                                                 vertices=[(-RDF(1.1074999999999999), RDF(1.149444444), RDF(9.0138888890000004))])}
            True
        """
        return self._P

    def base_ring(self):
        r"""
        Return the base_ring of the regions of the Voronoi diagram.

        EXAMPLES::

            sage: V = VoronoiDiagram([[1, 3, 1], [2, -2, 1], [-1, 2, 1/2]]); V.base_ring()
            Rational Field
            sage: V = VoronoiDiagram([[1, 3.14], [2, -2/3], [-1, 22]]); V.base_ring()
            Real Double Field
            sage: V = VoronoiDiagram([[1, 3], [2, 4]]); V.base_ring()
            Rational Field
        """
        return self._base_ring

    def _repr_(self):
        r"""
        Return a description of the Voronoi diagram.

        EXAMPLES::

            sage: V = VoronoiDiagram(polytopes.regular_polygon(3).vertices()); V         # optional - sage.rings.number_field
            The Voronoi diagram of 3 points of dimension 2 in the Algebraic Real Field
            sage: VoronoiDiagram([])                                                     # optional - sage.rings.number_field
            The empty Voronoi diagram.
        """
        if self._n:
            desc = 'The Voronoi diagram of ' + str(self._n)
            desc += ' points of dimension ' + str(self.ambient_dim())
            desc += ' in the ' + str(self.base_ring())
            return desc

        return 'The empty Voronoi diagram.'

    def plot(self, cell_colors=None, **kwds):
        """
        Return a graphical representation for 2-dimensional Voronoi diagrams.

        INPUT:

        - ``cell_colors`` -- (default: ``None``) provide the colors for the cells, either as
          dictionary. Randomly colored cells are provided with ``None``.
        - ``**kwds`` -- optional keyword parameters, passed on as arguments for
          plot().

        OUTPUT:

        A graphics object.

        EXAMPLES::

            sage: P = [[0.671, 0.650], [0.258, 0.767], [0.562, 0.406], [0.254, 0.709], [0.493, 0.879]]

            sage: V = VoronoiDiagram(P); S=V.plot()                                            # optional - sage.plot
            sage: show(S, xmin=0, xmax=1, ymin=0, ymax=1, aspect_ratio=1, axes=false)          # optional - sage.plot

            sage: S=V.plot(cell_colors={0:'red', 1:'blue', 2:'green', 3:'white', 4:'yellow'})  # optional - sage.plot
            sage: show(S, xmin=0, xmax=1, ymin=0, ymax=1, aspect_ratio=1, axes=false)          # optional - sage.plot

            sage: S=V.plot(cell_colors=['red','blue','red','white', 'white'])                  # optional - sage.plot
            sage: show(S, xmin=0, xmax=1, ymin=0, ymax=1, aspect_ratio=1, axes=false)          # optional - sage.plot

            sage: S=V.plot(cell_colors='something else')                                       # optional - sage.plot
            Traceback (most recent call last):
            ...
            AssertionError: 'cell_colors' must be a list or a dictionary


        Trying to plot a Voronoi diagram of dimension other than 2 gives an
        error::

            sage: VoronoiDiagram([[1, 2, 3], [6, 5, 4]]).plot()                                # optional - sage.plot
            Traceback (most recent call last):
            ...
            NotImplementedError: Plotting of 3-dimensional Voronoi diagrams not
            implemented
        """

        if self.ambient_dim() == 2:
            S = line([])

            if cell_colors is None:
                from random import shuffle
                cell_colors = rainbow(self._n)
                shuffle(cell_colors)
            else:
                if not (isinstance(cell_colors, list) or (isinstance(cell_colors, dict))):
                    raise AssertionError("'cell_colors' must be a list or a dictionary")
            for i, p in enumerate(self._P):
                col = cell_colors[i]
                S += (self.regions()[p]).render_solid(color=col, zorder=1)
                S += point(p, color=col, pointsize=10, zorder=3)
                S += point(p, color='black', pointsize=20, zorder=2)
            return plot(S, **kwds)
        raise NotImplementedError('Plotting of ' + str(self.ambient_dim()) +
                                  '-dimensional Voronoi diagrams not' +
                                  ' implemented')

    def _are_points_in_regions(self):
        """
        Check if all points are contained in their regions.

        EXAMPLES::

            sage: py_trips = [[a, b] for a in range(1, 50) for b in range(1, 50) if ZZ(a^2 + b^2).is_square()]
            sage: v = VoronoiDiagram(py_trips)
            sage: v._are_points_in_regions()
            True
        """
        return all(self.regions()[p].contains(p) for p in self.points())
