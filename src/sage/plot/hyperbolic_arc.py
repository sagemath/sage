"""
Arcs in hyperbolic geometry

AUTHORS:

- Hartmut Monien (2011 - 08)
"""
#*****************************************************************************
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
#*****************************************************************************

from sage.plot.bezier_path import BezierPath
from sage.misc.decorators import options, rename_keyword
from sage.rings.cc import CC


class HyperbolicArc(BezierPath):
    """
    Primitive class for hyperbolic arc type. See ``hyperbolic_arc?`` for
    information about plotting a hyperbolic arc in the complex plane.

    INPUT:

    - ``a, b`` - coordinates of the hyperbolic arc in the complex plane
    - ``options`` - dict of valid plot options to pass to constructor

    EXAMPLES:

    Note that constructions should use ``hyperbolic_arc``::

         sage: from sage.plot.hyperbolic_arc import HyperbolicArc

         sage: print(HyperbolicArc(0, 1/2+I*sqrt(3)/2, {}))
         Hyperbolic arc (0.000000000000000, 0.500000000000000 + 0.866025403784439*I)
    """

    def __init__(self, A, B, options):
        A, B = (CC(A), CC(B))
        if A.imag()<0:
            raise ValueError("%s is not a valid point in the UHP model"%(A))
        if B.imag()<0:
            raise ValueError("%s is not a valid point in the UHP model"%(B))
        self.path = []
        self._hyperbolic_arc(A, B, True)
        BezierPath.__init__(self, self.path, options)
        self.A, self.B = (A, B)

    def _repr_(self):
        """
        String representation of HyperbolicArc.

        TESTS::

            sage: from sage.plot.hyperbolic_arc import HyperbolicArc
            sage: HyperbolicArc(0, 1/2+I*sqrt(3)/2, {})._repr_()
            'Hyperbolic arc (0.000000000000000, 0.500000000000000 + 0.866025403784439*I)'
        """

        return "Hyperbolic arc (%s, %s)" % (self.A, self.B)

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
            self.path.append([(z0.real(),z0.imag()), (z3.real(),z3.imag())])
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
@options(alpha=1, fill=False, thickness=1, rgbcolor="blue", zorder=2, linestyle='solid')
def hyperbolic_arc(a, b, **options):
    """
    Plot an arc from a to b in hyperbolic geometry in the complex upper
    half plane.

    INPUT:

    - ``a, b`` - complex numbers in the upper half complex plane
      connected bye the arc

    OPTIONS:

    - ``alpha`` - default: 1

    - ``thickness`` - default: 1

    - ``rgbcolor`` - default: 'blue'

    - ``linestyle`` - (default: ``'solid'``) The style of the line, which is one
      of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or ``'--'``,
      ``':'``, ``'-'``, ``'-.'``, respectively.

    Examples:

    Show a hyperbolic arc from 0 to 1::

         sage: hyperbolic_arc(0, 1)
         Graphics object consisting of 1 graphics primitive

    Show a hyperbolic arc from 1/2 to `i` with a red thick line::

         sage: hyperbolic_arc(1/2, I, color='red', thickness=2)
         Graphics object consisting of 1 graphics primitive

    Show a hyperbolic arc form `i` to `2 i` with dashed line::

         sage: hyperbolic_arc(I, 2*I, linestyle='dashed')
         Graphics object consisting of 1 graphics primitive
         sage: hyperbolic_arc(I, 2*I, linestyle='--')
         Graphics object consisting of 1 graphics primitive
    """
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(g._extract_kwds_for_show(options))
    g.add_primitive(HyperbolicArc(a, b, options))
    g.set_aspect_ratio(1)
    return g
