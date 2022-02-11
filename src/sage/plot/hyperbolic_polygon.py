"""
Polygons and triangles in hyperbolic geometry

AUTHORS:

- Hartmut Monien (2011-08)
- Vincent Delecroix (2014-11)
"""
#*****************************************************************************
#       Copyright (C) 2011 Hartmut Monien <monien@th.physik.uni-bonn.de>,
#                     2014 Vincent Delecroix <20100.delecroix@gmail.com>,
#                     2015 Stefan Kraemer <skraemer@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.plot.bezier_path import BezierPath
from sage.misc.decorators import options, rename_keyword
from sage.rings.cc import CC


class HyperbolicPolygon(BezierPath):
    """
    Primitive class for hyperbolic polygon type.

    See ``hyperbolic_polygon?`` for information about plotting a hyperbolic
    polygon in the complex plane.

    INPUT:

    - ``pts`` -- coordinates of the polygon (as complex numbers)

    - ``options`` -- dict of valid plot options to pass to constructor

    EXAMPLES:

    Note that constructions should use :func:`hyperbolic_polygon` or
    :func:`hyperbolic_triangle`::

         sage: from sage.plot.hyperbolic_polygon import HyperbolicPolygon
         sage: print(HyperbolicPolygon([0, 1/2, I], {}))
         Hyperbolic polygon (0.000000000000000, 0.500000000000000, 1.00000000000000*I)
    """
    def __init__(self, pts, options):
        """
        Initialize HyperbolicPolygon.

        EXAMPLES::

            sage: from sage.plot.hyperbolic_polygon import HyperbolicPolygon
            sage: print(HyperbolicPolygon([0, 1/2, I], {}))
            Hyperbolic polygon (0.000000000000000, 0.500000000000000, 1.00000000000000*I)
        """
        pts = [CC(_) for _ in pts]
        self.path = []
        self._hyperbolic_arc(pts[0], pts[1], True)
        for i in range(1, len(pts) - 1):
            self._hyperbolic_arc(pts[i], pts[i + 1])
        self._hyperbolic_arc(pts[-1], pts[0])
        BezierPath.__init__(self, self.path, options)
        self._pts = pts

    def _repr_(self):
        """
        String representation of HyperbolicPolygon.

        TESTS::

            sage: from sage.plot.hyperbolic_polygon import HyperbolicPolygon
            sage: HyperbolicPolygon([0, 1/2, I], {})._repr_()
            'Hyperbolic polygon (0.000000000000000, 0.500000000000000, 1.00000000000000*I)'
        """
        return "Hyperbolic polygon ({})".format(", ".join(map(str, self._pts)))

    def _hyperbolic_arc(self, z0, z3, first=False):
        """
        Function to construct Bezier path as an approximation to
        the hyperbolic arc between the complex numbers z0 and z3 in the
        hyperbolic plane.
        """
        z0, z3 = (CC(z0), CC(z3))
        p = (abs(z0)*abs(z0)-abs(z3)*abs(z3))/(z0-z3).real()/2
        r = abs(z0-p)

        if abs(z3-z0)/r < 0.1:
            self.path.append([(z0.real(), z0.imag()), (z3.real(), z3.imag())])
            return

        if z0.imag() == 0 and z3.imag() == 0:
            p = (z0.real()+z3.real())/2
            zm = CC(p, r)
            self._hyperbolic_arc(z0, zm, first)
            self._hyperbolic_arc(zm, z3)
            return
        else:
            zm = ((z0+z3)/2-p)/abs((z0+z3)/2-p)*r+p
            t = (8*zm-4*(z0+z3)).imag()/3/(z3-z0).real()
            z1 = z0 + t*CC(z0.imag(), (p-z0.real()))
            z2 = z3 - t*CC(z3.imag(), (p-z3.real()))
        if first:
            self.path.append([(z0.real(), z0.imag()),
                              (z1.real(), z1.imag()),
                              (z2.real(), z2.imag()),
                              (z3.real(), z3.imag())])
            first = False
        else:
            self.path.append([(z1.real(), z1.imag()),
                              (z2.real(), z2.imag()),
                              (z3.real(), z3.imag())])


@rename_keyword(color='rgbcolor')
@options(alpha=1, fill=False, thickness=1, rgbcolor="blue", zorder=2,
         linestyle='solid')
def hyperbolic_polygon(pts, **options):
    r"""
    Return a hyperbolic polygon in the hyperbolic plane with vertices ``pts``.

    Type ``?hyperbolic_polygon`` to see all options.

    INPUT:

    - ``pts`` -- a list or tuple of complex numbers

    OPTIONS:

    - ``alpha`` -- default: 1

    - ``fill`` -- default: ``False``

    - ``thickness`` -- default: 1

    - ``rgbcolor`` -- default: ``'blue'``

    - ``linestyle`` -- (default: ``'solid'``) The style of the line, which is
      one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or ``'--'``,
      ``':'``, ``'-'``, ``'-.'``, respectively.

    EXAMPLES:

    Show a hyperbolic polygon with coordinates `-1`, `3i`, `2+2i`, `1+i`::

        sage: hyperbolic_polygon([-1,3*I,2+2*I,1+I])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_polygon([-1,3*I,2+2*I,1+I])
        sphinx_plot(P)

    With more options::

        sage: hyperbolic_polygon([-1,3*I,2+2*I,1+I], fill=True, color='red')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_polygon([-1,3*I,2+2*I,1+I], fill=True, color='red')
        sphinx_plot(P)

    """
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(g._extract_kwds_for_show(options))
    g.add_primitive(HyperbolicPolygon(pts, options))
    g.set_aspect_ratio(1)
    return g


def hyperbolic_triangle(a, b, c, **options):
    r"""
    Return a hyperbolic triangle in the hyperbolic plane with vertices ``(a,b,c)``.

    Type ``?hyperbolic_polygon`` to see all options.

    INPUT:

    - ``a, b, c`` -- complex numbers in the upper half complex plane

    OPTIONS:

    - ``alpha`` -- default: 1

    - ``fill`` -- default: ``False``

    - ``thickness`` -- default: 1

    - ``rgbcolor`` -- default: ``'blue'``

    - ``linestyle`` - (default: ``'solid'``) The style of the line, which is
      one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or ``'--'``,
      ``':'``, ``'-'``, ``'-.'``, respectively.

    EXAMPLES:

    Show a hyperbolic triangle with coordinates `0, 1/2+i\sqrt{3}/2` and
    `-1/2+i\sqrt{3}/2`::

         sage: hyperbolic_triangle(0, -1/2+I*sqrt(3)/2, 1/2+I*sqrt(3)/2)
         Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_triangle(0, 0.5*(-1+I*sqrt(3)), 0.5*(1+I*sqrt(3)))
        sphinx_plot(P)

    A hyperbolic triangle with coordinates `0, 1` and `2+i` and a dashed line::

         sage: hyperbolic_triangle(0, 1, 2+i, fill=true, rgbcolor='red', linestyle='--')
         Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_triangle(0, 1, 2+i, fill=true, rgbcolor='red', linestyle='--')
        sphinx_plot(P)

    """
    return hyperbolic_polygon((a, b, c), **options)
