# sage_setup: distribution = sagemath-modules
# sage.doctest: needs numpy
"""
Complex Interpolation

AUTHORS:

- Ethan Van Andel (2009): initial version

Development supported by NSF award No. 0702939.
"""

# ****************************************************************************
#       Copyright (C) 2009 Ethan Van Andel <evlutte@gmail.com>,
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
# ****************************************************************************

import numpy as np
cimport numpy as np

if int(np.version.short_version[0]) > 1:
    np.set_printoptions(legacy="1.25")

from math import pi
cdef double TWOPI = 2*pi


def polygon_spline(pts):
    """
    Create a polygon from a set of complex or `(x,y)` points. The polygon
    will be a parametric curve from 0 to 2*pi. The returned values will be
    complex, not `(x,y)`.

    INPUT:

    - ``pts`` -- list or array of complex numbers of tuples of the form `(x,y)`

    EXAMPLES:

    A simple square::

        sage: pts = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        sage: ps = polygon_spline(pts)
        sage: fx = lambda x: ps.value(x).real
        sage: fy = lambda x: ps.value(x).imag
        sage: show(parametric_plot((fx, fy), (0, 2*pi)))                                # needs sage.plot sage.symbolic
        sage: m = Riemann_Map([lambda x: ps.value(real(x))],
        ....:                 [lambda x: ps.derivative(real(x))], 0)
        sage: show(m.plot_colored() + m.plot_spiderweb())                               # needs sage.plot

    Polygon approximation of a circle::

        sage: # needs sage.symbolic
        sage: pts = [e^(I*t / 25) for t in range(25)]
        sage: ps = polygon_spline(pts)
        sage: ps.derivative(2)
        (-0.0470303661...+0.1520363883...j)
    """
    return PSpline(pts)


cdef class PSpline:
    """
    A ``CCSpline`` object contains a polygon interpolation of a figure
    in the complex plane.

    EXAMPLES:

    A simple ``square``::

        sage: pts = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        sage: ps = polygon_spline(pts)
        sage: ps.value(0)
        (-1-1j)
        sage: ps.derivative(0)
        (1.27323954...+0j)
    """
    cdef int N
    cdef np.ndarray pts

    def __init__(self, pts):
        """
        TESTS::

            sage: pts = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
            sage: ps = polygon_spline(pts)
        """
        if isinstance(pts[0], tuple):
            self.pts = np.array(
                [complex(i[0], i[1]) for i in pts], dtype=np.complex128)
        else:
            self.pts = np.array(pts, dtype=np.complex128)
        self.N = len(pts)

    def value(self, double t):
        """
        Return the derivative (speed and direction of the curve) of a
        given point from the parameter ``t``.

        INPUT:

        - ``t`` -- double; the parameter value for the parameterized curve,
          between 0 and 2*pi

        OUTPUT:

        A complex number representing the point on the polygon corresponding
        to the input ``t``.

        EXAMPLES:

        ::

            sage: pts = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
            sage: ps = polygon_spline(pts)
            sage: ps.value(.5)
            (-0.363380227632...-1j)
            sage: ps.value(0) - ps.value(2*RDF.pi())
            0j
            sage: ps.value(10)
            (0.26760455264...+1j)
        """
        cdef double t1 = ((t / TWOPI*self.N) % self.N)
        pt1 = self.pts[int(t1)]
        pt2 = self.pts[(int(t1) + 1) % self.N]
        return pt1 + (pt2 - pt1) * (t1 - int(t1))

    def derivative(self, double t):
        """
        Return the derivative (speed and direction of the curve) of a
        given point from the parameter ``t``.

        INPUT:

        - ``t`` -- double; the parameter value for the parameterized curve,
          between 0 and 2*pi

        OUTPUT:

        A complex number representing the derivative at the point on the
        polygon corresponding to the input ``t``.

        EXAMPLES:

        ::

            sage: pts = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
            sage: ps = polygon_spline(pts)
            sage: ps.derivative(1 / 3)
            (1.27323954473...+0j)
            sage: from math import pi
            sage: ps.derivative(0) - ps.derivative(2*pi)
            0j
            sage: ps.derivative(10)
            (-1.27323954473...+0j)
        """
        cdef double t1 = ((t / TWOPI*self.N) % self.N)
        pt1 = self.pts[int(t1)]
        pt2 = self.pts[(int(t1) + 1) % self.N]
        return (pt2 - pt1) * self.N / TWOPI


def complex_cubic_spline(pts):
    """
    Create a cubic spline interpolated figure from a set of complex or
    `(x,y)` points. The figure will be a parametric curve from 0 to 2*pi.
    The returned values will be complex, not `(x,y)`.

    INPUT:

    - ``pts`` -- list or array of complex numbers, or tuples of the form `(x,y)`

    EXAMPLES:

    A simple ``square``::

        sage: pts = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        sage: cs = complex_cubic_spline(pts)
        sage: fx = lambda x: cs.value(x).real
        sage: fy = lambda x: cs.value(x).imag
        sage: from math import pi
        sage: show(parametric_plot((fx, fy), (0, 2*pi)))                                # needs sage.plot
        sage: m = Riemann_Map([lambda x: cs.value(real(x))],
        ....:                 [lambda x: cs.derivative(real(x))], 0)
        sage: show(m.plot_colored() + m.plot_spiderweb())                               # needs sage.plot

    Polygon approximation of a circle::

        sage: from cmath import exp
        sage: pts = [exp(1j * t / 25) for t in range(25)]
        sage: cs = complex_cubic_spline(pts)
        sage: cs.derivative(2)
        (-0.0497765406583...+0.151095006434...j)
    """
    return CCSpline(pts)


cdef class CCSpline:
    """
    A ``CCSpline`` object contains a cubic interpolation of a figure
    in the complex plane.

    EXAMPLES:

    A simple ``square``::

        sage: pts = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        sage: cs = complex_cubic_spline(pts)
        sage: cs.value(0)
        (-1-1j)
        sage: cs.derivative(0)
        (0.9549296...-0.9549296...j)
    """
    cdef int N
    cdef np.ndarray avec, bvec, cvec, dvec

    #standard cubic interpolation method
    def __init__(self, pts):
        """
        TESTS::

            sage: pts = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
            sage: cs = complex_cubic_spline(pts)
        """
        if isinstance(pts[0], tuple):
            pts = np.array(
                [complex(pt[0], pt[1]) for pt in pts], dtype=np.complex128)
        cdef int N, i
        N = len(pts)
        yvec = np.zeros(N, dtype=np.complex128)
        for i in range(N):
            yvec[i] = 3 * (pts[(i - 1) % N] - 2*pts[i] + pts[(i + 1) % N])
        bmat = np.zeros([N, N], dtype=np.complex128)
        for i in range(N):
            bmat[i, i] = 4
            bmat[(i - 1) % N, i] = 1
            bmat[(i + 1) % N, i] = 1
        bvec = (np.linalg.solve(bmat, yvec))
        cvec = np.zeros(N, dtype=np.complex128)
        for i in range(N):
            cvec[i] = (pts[(i + 1) % N] - pts[i] - 1.0/3.0 *
                       bvec[(i + 1) % N] - 2./3. * bvec[i])
        dvec = np.array(pts, dtype=np.complex128)
        avec = np.zeros(N, dtype=np.complex128)
        for i in range(N):
            avec[i] = 1.0/3.0 * (bvec[(i + 1) % N] - bvec[i])
        self.avec = avec
        self.bvec = bvec
        self.cvec = cvec
        self.dvec = dvec
        self.N = N

    def value(self, double t):
        """
        Return the location of a given point from the parameter ``t``.

        INPUT:

        - ``t`` -- double; the parameter value for the parameterized curve,
          between 0 and 2*pi

        OUTPUT:

        A complex number representing the point on the figure corresponding
        to the input ``t``.

        EXAMPLES:

        ::

            sage: pts = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
            sage: cs = complex_cubic_spline(pts)
            sage: cs.value(4 / 7)
            (-0.303961332787...-1.34716728183...j)
            sage: from math import pi
            sage: cs.value(0) - cs.value(2*pi)
            0j
            sage: cs.value(-2.73452)
            (0.934561222231...+0.881366116402...j)
        """
        cdef double t1 = (t / TWOPI * self.N) % self.N
        cdef int j = int(t1)
        return (self.avec[j] * (t1 - j)**3 + self.bvec[j] * (t1 - j)**2 +
                self.cvec[j] * (t1 - j) + self.dvec[j])

    def derivative(self, double t):
        """
        Return the derivative (speed and direction of the curve) of a
        given point from the parameter ``t``.

        INPUT:

        - ``t`` -- double; the parameter value for the parameterized curve,
          between 0 and 2*pi

        OUTPUT:

        A complex number representing the derivative at the point on the
        figure corresponding to the input ``t``.

        EXAMPLES::

            sage: pts = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
            sage: cs = complex_cubic_spline(pts)
            sage: cs.derivative(3 / 5)
            (1.40578892327...-0.225417136326...j)
            sage: from math import pi
            sage: cs.derivative(0) - cs.derivative(2 * pi)
            0j
            sage: cs.derivative(-6)
            (2.52047692949...-1.89392588310...j)
        """
        cdef double t1 = (t / TWOPI * self.N) % self.N
        cdef int j = int(t1)
        return (3 * self.avec[j] * (t1 - j)**2 + 2 * self.bvec[j] * (t1 - j) +
                self.cvec[j]) / TWOPI * self.N
