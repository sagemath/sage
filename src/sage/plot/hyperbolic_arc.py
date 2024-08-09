# sage.doctest: needs sage.symbolic
r"""
Arcs in hyperbolic geometry

AUTHORS:

- Hartmut Monien (2011 - 08)

Three models of the hyperbolic plane are implemented:
Upper Half Plane, Poincaré Disk, and Klein Disk, each
with its different domain and metric tensor.

.. RUBRIC:: Upper half plane (UHP)

In this model, hyperbolic points are described by two coordinates, which
we will represent by a complex number in the domain

.. MATH::

    H = \{ z \in \CC \mid \Im(z)>0 \}

with the corresponding metric tensor

.. MATH::

   ds^2 = \frac{dzd\bar{z}}{\Im(z)^2}.

.. RUBRIC:: Poincaré disk (PD)

In this model, hyperbolic points are described by two coordinates, which we
will represent by a complex number within the unit circle, having therefore
the following domain

.. MATH::

    D = \{ z \in \CC \mid \lvert z \rvert < 1\}

with the corresponding metric tensor

.. MATH::

   ds^2 = 4 \frac{dzd\bar{z}}{(1-\lvert z \rvert^2)^2}.

.. RUBRIC:: Klein disk (KM)

In this model, the domain is again complex numbers within the unit circle as
in the Poincaré disk model, but the corresponding metric tensor is different:

.. MATH::

    ds^2 = \frac{dzd\bar{z}}{1-\lvert z \rvert^2}
    + \frac{(z \cdot dz)^2}{(1-\lvert z \rvert^2)^2}.

.. SEEALSO::

   :mod:`sage.geometry.hyperbolic_space.hyperbolic_geodesic`

REFERENCES:

For additional models of the hyperbolic plane and its relationship
see [CFKP1997]_. For a more detailed explanation on hyperbolic arcs
see [Sta1993]_.
"""
# *****************************************************************************
#       Copyright (C) 2011 Hartmut Monien <monien@th.physik.uni-bonn.de>,
#                     2015 Stefan Kraemer <skraemer@th.physik.uni-bonn.de>
#                     2016 Javier Honrubia <jhonrubia6@uned.es>
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
# *****************************************************************************

from sage.plot.bezier_path import BezierPath
from sage.plot.circle import circle
from sage.misc.decorators import options, rename_keyword
from sage.rings.cc import CC
from sage.geometry.hyperbolic_space.hyperbolic_constants import EPSILON


class HyperbolicArcCore(BezierPath):
    """
    Base class for Hyperbolic arcs and hyperbolic polygons in the
    hyperbolic plane.

    The Upper Half Model, Poincaré Disk Model, and Klein Disk model
    are supported.
    """
    def _bezier_path(self, z0, z1, model, first=False):
        """
        Construct a bezier path from a given arc object and store it
        in the ``path`` attribute.

        INPUT:

        - ``z0``, ``z1`` -- the points of the arc
        - ``model`` -- an model for the hyperbolic arc

        TESTS::

            sage: from sage.plot.hyperbolic_arc import HyperbolicArc
            sage: HyperbolicArc(0, 1/2+I*sqrt(3)/2, "UHP", {})  # indirect doctest
            Hyperbolic arc (0.000000000000000, 0.500000000000000 + 0.866025403784439*I)
        """
        import numpy as np
        from sage.rings.infinity import infinity
        EPSILON = 10 ** -5

        arc0 = model.get_geodesic(z0, z1).plot()[0]

        if isinstance(arc0, BezierPath):
            points = arc0.vertices
        else:
            points = arc0.bezier_path()[0].vertices
        if (
            ((z0.is_infinity() or z0 == infinity)
             and abs(CC(points[0][0], points[0][1]) - z1) < EPSILON)
            or ((z1.is_infinity() or z1 == infinity)
                and abs(CC(points[1][0], points[1][1]) - z0) < EPSILON)
            or (abs(CC(points[0][0], points[0][1]) - z0) >= EPSILON
                and not (z0.is_infinity() or z0 == infinity or z1.is_infinity()
                         or z1 == infinity))
            ):
            points = np.flipud(points)  # order is important

        if first:
            self.path.append(points[0:4])  # if it is a line it will append only two control points
            if isinstance(arc0, BezierPath):
                self.last_plotted = "line"
            else:
                N = 4
                # Add new triplets
                while N < len(points):
                    self.path.append(points[N: N + 3])
                    N += 3
                self.last_plotted = "arc"
        else:
            # the first point is equal to the last of the previous arc
            points = np.delete(points, 0, 0)
            N = 0
            if isinstance(arc0, BezierPath):
                self.path.append(points[0:1])
                self.last_plotted = "line"
            elif self.last_plotted == "line":  # actual segment is an arc
                # Add new triplets
                while N < len(points):
                    self.path.append(points[N: N + 3])
                    N += 3
                self.last_plotted = "arc"
            else:
                # Complete the last tuple of control points
                tail = self.path[-1]
                ltail = len(tail)
                while ltail < 3:
                    self.path[-1].append(points[N])
                    ltail += 1
                    N += 1
                # Add new triplets
                while N < len(points):
                    self.path.append(points[N: N + 3])
                    N += 3
                self.last_plotted = "arc"
        return


class HyperbolicArc(HyperbolicArcCore):
    r"""
    Primitive class for hyberbolic arc type.

    See ``hyperbolic_arc?`` for information about plotting a hyperbolic
    arc in the complex plane.

    INPUT:

    - ``A``, ``B`` -- end points of the hyperbolic arc
    - ``model`` -- the hyperbolic model used, which is one of the following:

      * ``'UHP'`` -- upper half plane
      * ``'PD'`` -- Poincaré disk
      * ``'KM'`` -- Klein disk

    TESTS::

         sage: from sage.plot.hyperbolic_arc import HyperbolicArc
         sage: HyperbolicArc(0, 1/2+I*sqrt(3)/2, "UHP", {})
         Hyperbolic arc (0.000000000000000, 0.500000000000000 + 0.866025403784439*I)
    """
    def __init__(self, A, B, model, options):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.plot.hyperbolic_arc import HyperbolicArc
            sage: arc = HyperbolicArc(0, 1/2+I*sqrt(3)/2, "UHP", {})
            sage: TestSuite(arc).run(skip='_test_pickling')  # no equality implemented
        """
        if model == "HM":
            raise ValueError("the hyperboloid model is not supported")
        from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicPlane
        HP = HyperbolicPlane()
        M = getattr(HP, model)()
        self.A = CC(A)
        self.B = CC(B)
        self.model = model
        M.point_test(self.A)
        M.point_test(self.B)
        self.path = []
        self._bezier_path(self.A, self.B, M, True)
        BezierPath.__init__(self, self.path, options)

    def _repr_(self):
        """
        String representation of HyperbolicArc.

        TESTS::

            sage: from sage.plot.hyperbolic_arc import HyperbolicArc
            sage: HyperbolicArc(0, 1/2+I*sqrt(3)/2, "UHP", {})
            Hyperbolic arc (0.000000000000000, 0.500000000000000 + 0.866025403784439*I)
        """
        return f"Hyperbolic arc ({self.A}, {self.B})"


@rename_keyword(color='rgbcolor')
@options(alpha=1, fill=False, thickness=1, rgbcolor='blue', zorder=2, linestyle='solid')
def hyperbolic_arc(a, b, model='UHP', **options):
    r"""
    Plot an arc from ``a`` to ``b`` in hyperbolic plane.

    INPUT:

    - ``a``, ``b`` -- complex numbers connected by a hyperbolic arc

    - ``model`` -- (default: ``'UHP'``) hyperbolic model used,
      which is one of the following:

      * ``'UHP'`` -- upper half plane
      * ``'PD'`` -- Poincaré disk
      * ``'KM'`` -- Klein disk
      * ``'HM'`` -- hyperboloid model

    OPTIONS:

    - ``alpha`` -- (default: 1)

    - ``thickness`` -- (default: 1)

    - ``rgbcolor`` -- (default: ``'blue'``)

    - ``linestyle`` -- (default: ``'solid'``) the style of the line, which
      is one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``,
      or ``'--'``, ``':'``, ``'-'``, ``'-.'``, respectively

    EXAMPLES:

    Show a hyperbolic arc from `0` to `1`::

        sage: hyperbolic_arc(0, 1)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(hyperbolic_arc(0,1))

    Show a hyperbolic arc from `1/2` to `i` with a red thick line::

        sage: hyperbolic_arc(0.5, I,color='red', thickness=2)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(hyperbolic_arc(0.5, I, color='red', thickness=2))

    Show a hyperbolic arc from `1+i` to `1+2i` with dashed line::

        sage: hyperbolic_arc(1+I, 1+2*I, linestyle='dashed', color='green')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(hyperbolic_arc(CC(1,1), CC(1,2), linestyle='dashed', color='green'))

    ::

         sage: hyperbolic_arc(-1+I, 1+2*I, linestyle='--', color='orange')
         Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(hyperbolic_arc(CC(-1,1), CC(1,2), linestyle='dashed'))

    Show a hyperbolic arc from a `1+i` to infinity::

        sage: hyperbolic_arc(1 + I, infinity, color='brown')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        from sage.rings.infinity import infinity
        sphinx_plot(hyperbolic_arc(CC(1,1), infinity, color='brown'))

    We can also plot hyperbolic arcs in other models.

    We show a hyperbolic arc from `i` to `-1` in red, another hyperbolic arc
    from `e^{i\pi/3}` to `0.6 \cdot e^{i 3\pi/4}` with dashed style in green,
    and finally a hyperbolic arc from `-0.5+0.5i` to `0.5-0.5i` together
    with the disk frontier in the Poincaré disk model::

        sage: z1 = CC(0,1)
        sage: z2 = CC(-1,0)
        sage: z3 = CC((cos(pi/3),sin(pi/3)))
        sage: z4 = CC((0.6*cos(3*pi/4),0.6*sin(3*pi/4)))
        sage: z5 = CC(-0.5,0.5)
        sage: z6 = CC(0.5,-0.5)
        sage: a1 = hyperbolic_arc(z1, z2, model='PD', color='red')
        sage: a2 = hyperbolic_arc(z3, z4, model='PD', color='green')
        sage: a3 = hyperbolic_arc(z5, z6, model='PD', linestyle='--')
        sage: a1 + a2 + a3
        Graphics object consisting of 6 graphics primitives

    .. PLOT::

        z1 = CC(0,1)
        z2 = CC(-1,0)
        z3 = CC((cos(pi/3),sin(pi/3)))
        z4 = CC((0.6*cos(3*pi/4),0.6*sin(3*pi/4)))
        z5 = CC(-0.5,0.5)
        z6 = CC(0.5,-0.5)
        a1 = hyperbolic_arc(z1, z2, model='PD', color='red')
        a2 = hyperbolic_arc(z3, z4, model='PD', color='green')
        a3 = hyperbolic_arc(z5, z6, model='PD', linestyle='--')
        P = a1 + a2 + a3
        sphinx_plot(P)

    We show the arcs defined by the same endpoints in the Klein disk
    model (note that these are *not* the image of those arcs when
    changing between the models)::

        sage: a1 = hyperbolic_arc(z1, z2, model='KM', color='red')
        sage: a2 = hyperbolic_arc(z3, z4, model='KM', color='green')
        sage: a3 = hyperbolic_arc(z5, z6, model='KM', linestyle='--')
        sage: a1 + a2 + a3
        Graphics object consisting of 6 graphics primitives

    .. PLOT::

        z1 = CC(0,1)
        z2 = CC(-1,0)
        z3 = CC((cos(pi/3),sin(pi/3)))
        z4 = CC((0.6*cos(3*pi/4),0.6*sin(3*pi/4)))
        z5 = CC(-0.5,0.5)
        z6 = CC(0.5,-0.5)
        a1 = hyperbolic_arc(z1, z2, model='KM', color='red')
        a2 = hyperbolic_arc(z3, z4, model='KM', color='green')
        a3 = hyperbolic_arc(z5, z6, model='KM', linestyle='--')
        P = a1 + a2 + a3
        sphinx_plot(P)

    Show a hyperbolic arc from `(1,2,\sqrt(6))` to `(-2,-3,\sqrt(14))`
    in the hyperboloid model::

        sage: a = (1,2,sqrt(6))
        sage: b = (-2,-3,sqrt(14))
        sage: hyperbolic_arc(a, b, model='HM')
        Graphics3d Object

    .. PLOT::

       a = (1,2,sqrt(6))
       b = (-2,-3,sqrt(14))
       sphinx_plot(hyperbolic_arc(a, b, model='HM'))
    """
    from sage.plot.graphics import Graphics

    g = Graphics()
    g._set_extra_kwds(g._extract_kwds_for_show(options))
    if model == "HM":
        # since KM is 3d we can not use HyperbolicArc class we plot it directly
        # and we also handle the hyperbolic_polygon in direct way
        from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicPlane

        # Check for valid points
        if a[2] < 0 or a[0]**2+a[1]**2-a[2]**2 + 1 > EPSILON:
            raise ValueError(f"{a} is not a valid point in the HM model")
        if b[2] < 0 or b[0]**2+b[1]**2-b[2]**2 + 1 > EPSILON:
            raise ValueError(f"{b} is not a valid point in the HM model")

        HM = HyperbolicPlane().HM()
        geodesic = HM.get_geodesic(a, b)
        g = g + geodesic.plot(show_hyperboloid=True, graphics_options=options)
    else:
        g.add_primitive(HyperbolicArc(a, b, model, options))
        if model == "PD" or model == "KM":
            g = g + circle((0, 0), 1, axes=False, color='black')
            g.set_aspect_ratio(1)
    return g
