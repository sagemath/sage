# sage_setup: distribution = sagemath-modules
# sage.doctest: needs numpy sage.symbolic
"""
Riemann Mapping

AUTHORS:

- Ethan Van Andel (2009-2011): initial version and upgrades

- Robert Bradshaw (2009): his "complex_plot" was adapted for plot_colored

Development supported by NSF award No. 0702939.
"""

# ****************************************************************************
#       Copyright (C) 2011 Ethan Van Andel <evlutte@gmail.com>,
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

from cysignals.signals cimport sig_on, sig_off

from sage.misc.decorators import options

from sage.ext.fast_callable import fast_callable

from sage.rings.complex_double import CDF

from sage.arith.srange import srange

from sage.calculus.interpolation import spline

from sage.calculus.integration import numerical_integral


import numpy as np
cimport numpy as np

if int(np.version.short_version[0]) > 1:
    np.set_printoptions(legacy="1.25")

from math import pi
from math import sin
from math import cos

from math import log # used for complex plot lightness
from math import atan

from cmath import exp
from cmath import phase

FLOAT = np.float64
ctypedef np.float64_t FLOAT_T

COMPLEX = np.complex128
ctypedef np.complex128_t COMPLEX_T

cdef FLOAT_T PI = pi
cdef FLOAT_T TWOPI = 2*PI
cdef COMPLEX_T I = complex(0,1)

cdef class Riemann_Map:
    r"""
    The ``Riemann_Map`` class computes an interior or exterior Riemann map,
    or an Ahlfors map of a region given by the supplied boundary curve(s)
    and center point. The class also provides various methods to
    evaluate, visualize, or extract data from the map.

    A Riemann map conformally maps a simply connected region in
    the complex plane to the unit disc. The Ahlfors map does the same thing
    for multiply connected regions.

    Note that all the methods are numerical. As a result all answers have
    some imprecision. Moreover, maps computed with small number of
    collocation points, or for unusually shaped regions, may be very
    inaccurate. Error computations for the ellipse can be found in the
    documentation for :meth:`analytic_boundary` and :meth:`analytic_interior`.

    [BSV2010]_ provides an overview of the Riemann map and discusses the research
    that lead to the creation of this module.

    INPUT:

    - ``fs`` -- list of the boundaries of the region, given as
      complex-valued functions with domain `0` to `2*pi`. Note that the
      outer boundary must be parameterized counter clockwise
      (i.e. ``e^(I*t)``) while the inner boundaries must be clockwise
      (i.e. ``e^(-I*t)``).

    - ``fprimes`` -- list of the derivatives of the boundary functions
      (Must be in the same order as ``fs``)

    - ``a`` -- complex, the center of the Riemann map. Will be mapped to
      the origin of the unit disc. Note that ``a`` MUST be within
      the region in order for the results to be mathematically valid.

    The following inputs may be passed in as named parameters:

    - ``N`` -- integer (default: `500`); the number of collocation points
      used to compute the map. More points will give more accurate results,
      especially near the boundaries, but will take longer to compute.

    - ``exterior`` -- boolean (default: ``False``); if set to ``True``, the
      exterior map will be computed, mapping the exterior of the region to the
      exterior of the unit circle.

    The following inputs may be passed as named parameters in unusual
    circumstances:

    - ``ncorners`` -- integer (default: `4`); if mapping a figure with
      (equally t-spaced) corners -- corners that make a significant change in
      the direction of the boundary -- better results may be sometimes obtained by
      accurately giving this parameter. Used to add the proper constant to
      the theta correspondence function.

    - ``opp`` -- boolean (default: ``False``); set to ``True`` in very rare
      cases where the theta correspondence function is off by ``pi``, that
      is, if red is mapped left of the origin in the color plot.

    EXAMPLES:

    The unit circle identity map::

        sage: f(t) = e^(I*t)
        sage: fprime(t) = I*e^(I*t)
        sage: m = Riemann_Map([f], [fprime], 0)  # long time (4 sec)
        sage: m.plot_colored() + m.plot_spiderweb()  # long time
        Graphics object consisting of 22 graphics primitives

    The exterior map for the unit circle::

        sage: m = Riemann_Map([f], [fprime], 0, exterior=True)  # long time (4 sec)
        sage: #spiderwebs are not supported for exterior maps
        sage: m.plot_colored() # long time
        Graphics object consisting of 1 graphics primitive

    The unit circle with a small hole::

        sage: f(t) = e^(I*t)
        sage: fprime(t) = I*e^(I*t)
        sage: hf(t) = 0.5*e^(-I*t)
        sage: hfprime(t) = 0.5*-I*e^(-I*t)
        sage: m = Riemann_Map([f, hf], [fprime, hfprime], 0.5 + 0.5*I)
        sage: #spiderweb and color plots cannot be added for multiply
        sage: #connected regions. Instead we do this.
        sage: m.plot_spiderweb(withcolor = True)  # long time
        Graphics object consisting of 3 graphics primitives

    A square::

        sage: ps = polygon_spline([(-1, -1), (1, -1), (1, 1), (-1, 1)])
        sage: f = lambda t: ps.value(real(t))
        sage: fprime = lambda t: ps.derivative(real(t))
        sage: m = Riemann_Map([f], [fprime], 0.25, ncorners=4)
        sage: m.plot_colored() + m.plot_spiderweb()  # long time
        Graphics object consisting of 22 graphics primitives

    Compute rough error for this map::

        sage: x = 0.75  # long time
        sage: print("error = {}".format(m.inverse_riemann_map(m.riemann_map(x)) - x))  # long time
        error = (-0.000...+0.0016...j)

    A fun, complex region for demonstration purposes::

        sage: f(t) = e^(I*t)
        sage: fp(t) = I*e^(I*t)
        sage: ef1(t) = .2*e^(-I*t) +.4+.4*I
        sage: ef1p(t) = -I*.2*e^(-I*t)
        sage: ef2(t) = .2*e^(-I*t) -.4+.4*I
        sage: ef2p(t) = -I*.2*e^(-I*t)
        sage: pts = [(-.5,-.15-20/1000),(-.6,-.27-10/1000),(-.45,-.45),(0,-.65+10/1000),(.45,-.45),(.6,-.27-10/1000),(.5,-.15-10/1000),(0,-.43+10/1000)]
        sage: pts.reverse()
        sage: cs = complex_cubic_spline(pts)
        sage: mf = lambda x:cs.value(x)
        sage: mfprime = lambda x: cs.derivative(x)
        sage: m = Riemann_Map([f,ef1,ef2,mf],[fp,ef1p,ef2p,mfprime],0,N = 400) # long time
        sage: p = m.plot_colored(plot_points = 400) # long time

    ALGORITHM:

    This class computes the Riemann Map via the Szego kernel using an
    adaptation of the method described by [KT1986]_.
    """
    cdef int N, B, ncorners
    cdef f
    cdef opp
    cdef COMPLEX_T a
    cdef np.ndarray tk, tk2
    cdef np.ndarray cps, dps, szego, p_vector, pre_q_vector
    cdef np.ndarray p_vector_inverse, sinalpha, cosalpha, theta_array
    cdef x_range, y_range
    cdef exterior

    def __init__(self, fs, fprimes, COMPLEX_T a, int N=500, int ncorners=4,
        opp=False, exterior = False):
        """
        Initialize the ``Riemann_Map`` class. See the class :class:`Riemann_Map`
        for full documentation on the input of this initialization method.

        TESTS::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0, N = 10)
        """
        cdef double xmax, xmin, ymax, ymin, space
        cdef int i, k
        if N <= 2:
            raise ValueError(
                "The number of collocation points must be > 2.")
        if exterior and a!=0:
            raise ValueError("The exterior map requires a=0")
        try:
            fs = [fast_callable(f, domain=CDF) for f in fs]
            fprimes = [fast_callable(f, domain=CDF) for f in fprimes]
        except AttributeError:
            pass
        self.f = fs[0]
        self.a = a
        self.ncorners = ncorners
        self.N = N  # Number of collocation pts
        self.opp = opp
        self.exterior = exterior
        self.tk = np.array(np.arange(N) * TWOPI / N + 0.001 / N,
                           dtype=FLOAT)
        self.tk2 = np.zeros(N + 1, dtype=FLOAT)
        for i in range(N):
            self.tk2[i] = self.tk[i]
        self.tk2[N] = TWOPI
        self.B = len(fs) # number of boundaries of the figure
        if self.exterior and (self.B > 1):
            raise ValueError(
                "The exterior map is undefined for multiply connected domains")
        cdef np.ndarray[COMPLEX_T,ndim=2] cps = np.zeros([self.B, N],
            dtype=COMPLEX)
        cdef np.ndarray[COMPLEX_T,ndim=2] dps = np.zeros([self.B, N],
            dtype=COMPLEX)
        # Find the points on the boundaries and their derivatives.
        if self.exterior:
            for k in range(self.B):
                for i in range(N):
                    fk = fs[k](self.tk[N-i-1])
                    cps[k, i] = complex(1/fk)
                    dps[k, i] = complex(1/fk**2*fprimes[k](self.tk[N-i-1]))
        else:
            for k in range(self.B):
                for i in range(N):
                    cps[k, i] = complex(fs[k](self.tk[i]))
                    dps[k, i] = complex(fprimes[k](self.tk[i]))
        if self.exterior:
            xmax = (1/cps).real.max()
            xmin = (1/cps).real.min()
            ymax = (1/cps).imag.max()
            ymin = (1/cps).imag.min()
        else:
            xmax = cps.real.max()
            xmin = cps.real.min()
            ymax = cps.imag.max()
            ymin = cps.imag.min()
        space = 0.1 * max(xmax - xmin, ymax - ymin)
        #The default plotting window for this map.
        self.cps = cps
        self.dps = dps
        self.x_range = (xmin - space, xmax + space)
        self.y_range = (ymin - space, ymax + space)
        # Now we do some more computation
        self._generate_theta_array()
        self._generate_interior_mapper()
        self._generate_inverse_mapper()

    def _repr_(self):
        """
        Return a string representation of this :class:`Riemann_Map` object.

        TESTS::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: isinstance(Riemann_Map([f], [fprime], 0, N = 10)._repr_(), str)  # long time
            True
        """
        return "A Riemann or Ahlfors mapping of a figure to the unit circle."

    cdef _generate_theta_array(self):
        """
        Generates the essential data for the Riemann map, primarily the
        Szegő kernel and boundary correspondence.

        See [KT1986]_ for the algorithm.

        TESTS::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0, N = 10)
        """
        cdef np.ndarray[COMPLEX_T, ndim=1] cp = self.cps.flatten()
        cdef np.ndarray[COMPLEX_T, ndim=1] dp = self.dps.flatten()
        cdef int N = self.N
        cdef int NB = N * self.B
        cdef int B = self.B
        cdef int i, k
        cdef FLOAT_T saa, t0
        cdef np.ndarray[FLOAT_T, ndim=1] adp, sadp
        cdef np.ndarray[COMPLEX_T, ndim=1] h, hconj, g, normalized_dp, C, phi
        cdef np.ndarray[COMPLEX_T, ndim=2] K
        cdef np.ndarray[FLOAT_T, ndim=2] theta_array
        # Setting things up to use the Nystrom method
        adp = abs(dp)
        sadp = np.sqrt(adp)
        h = 1 / (TWOPI * I) * ((dp / adp) / (self.a - cp))
        hconj = h.conjugate()
        g = -sadp * hconj
        normalized_dp=dp/adp
        C = I / N * sadp # equivalent to -TWOPI / N * 1 / (TWOPI * I) * sadp
        errinvalid = np.geterr()['invalid'] # checks the current error handling for invalid
        errdivide = np.geterr()['divide'] # checks the current error handling for divide
        np.seterr(divide='ignore',invalid='ignore')
        K = np.array([C * sadp[t] * (normalized_dp/(cp-cp[t]) -
             (normalized_dp[t]/(cp-cp[t])).conjugate())
              for t in np.arange(NB)], dtype=np.complex128)
        np.seterr(divide=errdivide,invalid=errinvalid) # resets the error handling
        for i in range(NB):
            K[i, i] = 1
        # Nystrom Method for solving 2nd kind integrals
        phi = np.linalg.solve(K, g) / NB * TWOPI
        # the all-important Szego kernel
        szego = np.array(phi.flatten() / np.sqrt(dp), dtype=COMPLEX)
        self.szego = szego.reshape([B, N])
        start = 0
        # Finding the theta correspondence using phase. Misbehaves for some
        # regions.
        if B != 1:
            theta_array = np.zeros([1, NB])
            for i in range(NB):
                theta_array[0, i] = phase(-I * np.power(phi[i], 2) * dp[i])
            self.theta_array = np.concatenate(
                [theta_array.reshape([B, N]), np.zeros([B, 1])], axis=1)
            for k in range(B):
                self.theta_array[k, N] = self.theta_array[k, 0] + TWOPI
        # Finding the theta correspondence using abs. Well behaved, but
        # doesn't work on multiply connected domains.
        else:
            phi2 = phi.reshape([self.B, N])
            theta_array = np.zeros([B, N + 1], dtype=np.float64)
            for k in range(B):
                phik = phi2[k]
                saa = (np.dot(abs(phi), abs(phi))) * TWOPI / NB
                theta_array[k, 0] = 0
                for i in range(1, N):
                    theta_array[k, i] = (
                        theta_array[k, i - 1] +
                        ((TWOPI / NB * TWOPI *
                          abs(np.power(phi[1 * i], 2)) / saa +
                          TWOPI / NB * TWOPI *
                          abs(np.power(phi[1 * i - 1], 2)) / saa)) / 2)
                tmax = int(0.5 * N / self.ncorners)
                # Finding the initial value of the theta function.
                phimax = -I * phik[tmax]**2 * self.dps[k, tmax]
                if self.opp:
                    t0 = theta_array[k, tmax] + phase(phimax)
                else:
                    t0 = theta_array[k, tmax] - phase(phimax)
                for i in range(N):
                    theta_array[k, i] = theta_array[k, i] - t0
                theta_array[k, N] = TWOPI + theta_array[k, 0]
            self.theta_array = theta_array

    def get_szego(self, int boundary=-1, absolute_value=False):
        """
        Return a discretized version of the Szego kernel for each boundary
        function.

        INPUT:

        The following inputs may be passed in as named parameters:

        - ``boundary`` -- integer (default: `-1`); if < 0,
          :meth:`get_theta_points` will return the points for all boundaries.
          If >= 0, :meth:`get_theta_points` will return only the points for
          the boundary specified.

        - ``absolute_value`` -- boolean (default: ``False``); if ``True``, will
          return the absolute value of the (complex valued) Szego kernel
          instead of the kernel itself. Useful for plotting.

        OUTPUT:

        A list of points of the form
        ``[t value, value of the Szego kernel at that t]``.

        EXAMPLES:

        Generic use::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
            sage: sz = m.get_szego(boundary=0)
            sage: points = m.get_szego(absolute_value=True)
            sage: list_plot(points)                                                     # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        Extending the points by a spline::

            sage: s = spline(points)
            sage: s(3*pi / 4)
            0.0012158...
            sage: plot(s,0,2*pi)  # plot the kernel                                     # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        The unit circle with a small hole::

            sage: f(t) = e^(I*t)
            sage: fprime(t) = I*e^(I*t)
            sage: hf(t) = 0.5*e^(-I*t)
            sage: hfprime(t) = 0.5*-I*e^(-I*t)
            sage: m = Riemann_Map([f, hf], [fprime, hfprime], 0.5 + 0.5*I)

        Getting the szego for a specific boundary::

            sage: sz0 = m.get_szego(boundary=0)
            sage: sz1 = m.get_szego(boundary=1)
        """
        cdef int k, B
        if boundary < 0:
            temptk = self.tk
            for i in range(self.B - 1):
                temptk = np.concatenate([temptk, self.tk])
            if absolute_value:
                return np.column_stack(
                    [temptk, abs(self.szego.flatten())]).tolist()
            else:
                return np.column_stack([temptk, self.szego.flatten()]).tolist()
        else:
            if absolute_value:
                return np.column_stack(
                    [self.tk, abs(self.szego[boundary])]).tolist()
            else:
                return np.column_stack(
                    [self.tk, self.szego[boundary]]).tolist()

    def get_theta_points(self, int boundary=-1):
        """
        Return an array of points of the form
        ``[t value, theta in e^(I*theta)]``, that is, a discretized version
        of the theta/boundary correspondence function. In other words, a point
        in this array [t1, t2] represents that the boundary point given by f(t1)
        is mapped to a point on the boundary of the unit circle given by e^(I*t2).

        For multiply connected domains, ``get_theta_points`` will list the
        points for each boundary in the order that they were supplied.

        INPUT:

        The following input must all be passed in as named parameters:

        - ``boundary`` -- integer (default: `-1`); if < 0,
          ``get_theta_points()`` will return the points for all boundaries.
          If >= 0, ``get_theta_points()`` will return only the points for
          the boundary specified.

        OUTPUT:

        A list of points of the form ``[t value, theta in e^(I*theta)]``.

        EXAMPLES:

        Getting the list of points, extending it via a spline, getting the
        points for only the outside of a multiply connected domain::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
            sage: points = m.get_theta_points()
            sage: list_plot(points)                                                     # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        Extending the points by a spline::

            sage: s = spline(points)
            sage: s(3*pi / 4)
            1.627660...

        The unit circle with a small hole::

            sage: f(t) = e^(I*t)
            sage: fprime(t) = I*e^(I*t)
            sage: hf(t) = 0.5*e^(-I*t)
            sage: hfprime(t) = 0.5*-I*e^(-I*t)
            sage: m = Riemann_Map([f, hf], [hf, hfprime], 0.5 + 0.5*I)

        Getting the boundary correspondence for a specific boundary::

            sage: tp0 = m.get_theta_points(boundary=0)
            sage: tp1 = m.get_theta_points(boundary=1)
        """
        if boundary < 0:
            temptk = self.tk2
            for i in range(self.B - 1):
                temptk = np.concatenate([temptk, self.tk2])
            return np.column_stack(
                [temptk, self.theta_array.flatten()]).tolist()
        else:
            return np.column_stack(
                [self.tk2, self.theta_array[boundary]]).tolist()

    cdef _generate_interior_mapper(self):
        """
        Generates the data necessary to use the :meth:`riemann_map` function.
        As much setup as possible is done here to minimize the computation
        that must be done in ``riemann_map``.

        TESTS::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0, N = 10) # indirect doctest
        """
        cdef int N = self.N
        cdef double complex coeff = 3*I / (8*N)
        dps = self.dps
        theta_array = self.theta_array
        cdef np.ndarray[double complex, ndim=2] p_vector = np.zeros(
            [self.B, N + 1], dtype=np.complex128)
        cdef int k, i
        # Lots of setup for Simpson's method of integration.
        for k in range(self.B):
            for i in range(N // 3):
                p_vector[k, 3*i] = (2*coeff * dps[k, 3*i] *
                                    exp(I * theta_array[k, 3*i]))
                p_vector[k, 3*i + 1] = (3*coeff * dps[k, 3*i + 1] *
                                        exp(I * theta_array[k, 3*i + 1]))
                p_vector[k, 3*i + 2] = (3*coeff * dps[k, 3*i + 2] *
                                        exp(I * theta_array[k, 3*i + 2]))
            p_vector[k, 0] = 1*coeff * dps[k, 0] * exp(I * theta_array[k, 0])
            if N % 3 == 0:
                p_vector[k, N] = 1*coeff * dps[k, 0] * exp(I*theta_array[k, 0])
            elif (N - 2) % 3 == 0:
                p_vector[k, N - 2] = ((coeff + I/(3*N)) * dps[k, N - 2] *
                                      exp(I * theta_array[k, N - 2]))
                p_vector[k, N - 1] = (4*I / (3*N) * dps[k, N - 1] *
                                      exp(I * theta_array[k, N - 1]))
                p_vector[k, N] = (I / (3*N) * dps[k, 0] *
                                  exp(I * theta_array[k, 0]))
            else:
                p_vector[k, N - 4] = ((coeff + I / (3*N)) * dps[k, N - 4] *
                                      exp(I * theta_array[k, N - 4]))
                p_vector[k, N - 3] = (4*I / (3*N) * dps[k, N - 3] *
                                      exp(I * theta_array[k, N - 3]))
                p_vector[k, N - 2] = (2*I / (3*N) * dps[k, N - 2] *
                                      exp(I * theta_array[k, N - 2]))
                p_vector[k, N - 1] = (4*I / (3*N) * dps[k, N - 1] *
                                      exp(I * theta_array[k, N - 1]))
                p_vector[k, N] = (I / (3*N) * dps[k, 0] *
                                  exp(I * theta_array[k, 0]))
        self.p_vector = p_vector.flatten()
        cdef np.ndarray[double complex, ndim=1] pq = self.cps[:,list(range(N))+[0]].flatten()
        self.pre_q_vector = pq

    cpdef riemann_map(self, COMPLEX_T pt):
        """
        Return the Riemann mapping of a point.

        That is, given ``pt`` on the interior of the mapped region,
        ``riemann_map`` will return the point on the unit disk that
        ``pt`` maps to. Note that this method only works for interior
        points; accuracy breaks down very close to the boundary. To
        get boundary correspondence, use :meth:`get_theta_points`.

        INPUT:

        - ``pt`` -- a complex number representing the point to be
          inverse mapped

        OUTPUT:

        A complex number representing the point on the unit circle that
        the input point maps to.

        EXAMPLES:

        Can work for different types of complex numbers::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
            sage: m.riemann_map(0.25 + sqrt(-0.5))
            (0.137514...+0.876696...j)
            sage: I = CDF.gen()
            sage: m.riemann_map(1.3*I)
            (-1.56...e-05+0.989694...j)
            sage: m.riemann_map(0.4)
            (0.73324...+3.2...e-06j)
            sage: m.riemann_map(complex(-3, 0.0001))
            (1.405757...e-05+8.06...e-10j)
        """

        cdef COMPLEX_T pt1
        cdef np.ndarray[COMPLEX_T, ndim=1] q_vector
        if self.exterior:
            pt1 = 1/pt
            q_vector = 1 / (
                self.pre_q_vector - pt1)
            return -1/np.dot(self.p_vector, q_vector)
        else:
            pt1 = pt
            q_vector = 1 / (
                self.pre_q_vector - pt1)
            return -np.dot(self.p_vector, q_vector)

    cdef _generate_inverse_mapper(self):
        """
        Generates the data necessary to use the
        :meth:`inverse_riemann_map` function. As much setup as possible is
        done here to minimize the computation that must be done in
        ``inverse_riemann_map``.

        TESTS::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0, N = 10) # indirect doctest
        """
        cdef int N = self.N
        cdef int B = self.B
        cdef float di
        theta_array = self.theta_array
        self.p_vector_inverse = np.zeros([B, N], dtype=np.complex128)
        # Setup for trapezoid integration because integration points are
        # not equally spaced.
        for k in range(B):
            for i in range(N):
                di = theta_array[k, (i + 1) % N] - theta_array[k, (i - 1) % N]
                if di > PI:
                    di = di - TWOPI
                elif di < -PI:
                    di = di + TWOPI
                self.p_vector_inverse[k, i] = di / 2
        self.sinalpha = np.zeros([B, N], dtype=np.float64)
        for k in range(B):
            for i in range(N):
                self.sinalpha[k, i] = sin(-theta_array[k, i])
        self.cosalpha = np.zeros([B, N], dtype=np.float64)
        for k in range(B):
            for i in range(N):
                self.cosalpha[k, i] = cos(-theta_array[k, i])

    cpdef inverse_riemann_map(self, COMPLEX_T pt):
        """
        Return the inverse Riemann mapping of a point.

        That is, given ``pt`` on the interior of the unit disc,
        ``inverse_riemann_map()`` will return the point on the
        original region that would be Riemann mapped to ``pt``. Note
        that this method does not work for multiply connected domains.

        INPUT:

        - ``pt`` -- a complex number (usually with absolute value <= 1)
          representing the point to be inverse mapped

        OUTPUT: the point on the region that Riemann maps to the input point

        EXAMPLES:

        Can work for different types of complex numbers::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
            sage: m.inverse_riemann_map(0.5 + sqrt(-0.5))
            (0.406880...+0.3614702...j)
            sage: m.inverse_riemann_map(0.95)
            (0.486319...-4.90019052...j)
            sage: m.inverse_riemann_map(0.25 - 0.3*I)
            (0.1653244...-0.180936...j)
            sage: m.inverse_riemann_map(complex(-0.2, 0.5))
            (-0.156280...+0.321819...j)
        """
        if self.exterior:
            pt = 1/pt
        r = abs(pt)
        if r == 0:
            stheta = 0
            ctheta = 0
        else:
            stheta = pt.imag / r
            ctheta = pt.real / r
        k = 0
        mapped = (1 - r**2) / TWOPI * np.dot(
                self.p_vector_inverse[k] * self.cps[k],
                1 / (1 + r**2 - 2*r *
                     (ctheta * self.cosalpha[k] - stheta * self.sinalpha[k])))
        if self.exterior:
            return 1/mapped
        else:
            return mapped

    def plot_boundaries(self, plotjoined=True, rgbcolor=[0,0,0], thickness=1):
        """
        Plots the boundaries of the region for the Riemann map. Note that
        this method DOES work for multiply connected domains.

        INPUT:

        The following inputs may be passed in as named parameters:

        - ``plotjoined`` -- boolean (default: ``True``); if ``False``,
          discrete points will be drawn; otherwise they will be connected
          by lines. In this case, if ``plotjoined=False``, the points shown
          will be the original collocation points used to generate the
          Riemann map.

        - ``rgbcolor`` -- float array (default: ``[0,0,0]``) the
          red-green-blue color of the boundary

        - ``thickness`` -- positive float (default: ``1``) the thickness of
          the lines or points in the boundary

        EXAMPLES:

        General usage::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)

        Default plot::

            sage: m.plot_boundaries()                                                   # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        Big blue collocation points::

            sage: m.plot_boundaries(plotjoined=False, rgbcolor=[0,0,1], thickness=6)    # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.all import list_plot

        plots = list(range(self.B))
        for k in range(self.B):
            # This conditional should be eliminated when the thickness/pointsize
            # issue is resolved later. Same for the others in plot_spiderweb().
            if plotjoined:
                plots[k] = list_plot(
                    comp_pt(self.cps[k], 1), plotjoined=True,
                    rgbcolor=rgbcolor, thickness=thickness)
            else:
                plots[k] = list_plot(
                    comp_pt(self.cps[k], 1), rgbcolor=rgbcolor,
                    pointsize=thickness)
        return sum(plots)

    cpdef compute_on_grid(self, plot_range, int x_points):
        """
        Compute the Riemann map on a grid of points.

        Note that these points are complex of the form z = x + y*i.

        INPUT:

        - ``plot_range`` -- tuple of the form ``[xmin, xmax, ymin, ymax]``;
          if the value is ``[]``, the default plotting window of the map will
          be used

        - ``x_points`` -- integer; the size of the grid in the x direction;
          the number of points in the y direction is scaled accordingly

        OUTPUT:

        - a tuple containing ``[z_values, xmin, xmax, ymin, ymax]``
          where ``z_values`` is the evaluation of the map on the specified grid.

        EXAMPLES:

        General usage::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
            sage: data = m.compute_on_grid([],5)
            sage: data[0][8,1]
            (-0.0879...+0.9709...j)
        """
        cdef FLOAT_T xmin, xmax, xstep, ymin, ymax, ystep
        cdef int y_points
        if plot_range == []:
            xmin = self.x_range[0]
            xmax = self.x_range[1]
            ymin = self.y_range[0]
            ymax = self.y_range[1]
        else:
            xmin = plot_range[0]
            xmax = plot_range[1]
            ymin = plot_range[2]
            ymax = plot_range[3]
        xstep = (xmax - xmin) / x_points
        ystep = xstep
        y_points = int((ymax-ymin)/ ystep)
        cdef Py_ssize_t i, j
        cdef COMPLEX_T pt
        cdef np.ndarray[COMPLEX_T] pre_q_vector = self.pre_q_vector
        cdef np.ndarray[COMPLEX_T] p_vector = self.p_vector
        cdef np.ndarray[COMPLEX_T, ndim=2] z_values = np.empty(
            [y_points, x_points], dtype=np.complex128)
        if self.exterior:
            for i in range(x_points):
                for j in range(y_points):
                    pt = 1/(xmin + 0.5*xstep + i*xstep + I*(ymin + 0.5*ystep + j*ystep))
                    z_values[j, i] = 1/(-np.dot(p_vector,1/(pre_q_vector - pt)))
        else:
            for i in range(x_points):
                for j in range(y_points):
                    pt = xmin + 0.5*xstep + i*xstep + I*(ymin + 0.5*ystep + j*ystep)
                    z_values[j, i] = -np.dot(p_vector,1/(pre_q_vector - pt))
        return z_values, xmin, xmax, ymin, ymax

    @options(interpolation='catrom')
    def plot_spiderweb(self, spokes=16, circles=4, pts=32, linescale=0.99,
            rgbcolor=[0, 0, 0], thickness=1, plotjoined=True, withcolor=False,
            plot_points=200, min_mag=0.001, **options):
        """
        Generate a traditional "spiderweb plot" of the Riemann map.

        This shows what concentric circles and radial lines map to.
        The radial lines may exhibit erratic behavior near the
        boundary; if this occurs, decreasing ``linescale`` may
        mitigate the problem.

        For multiply connected domains the spiderweb is by necessity
        generated using the forward mapping. This method is more
        computationally intensive. In addition, these spiderwebs cannot
        be ``added`` to color plots. Instead the ``withcolor`` option
        must be used.

        In addition, spiderweb plots are not currently supported for
        exterior maps.

        INPUT:

        The following inputs may be passed in as named parameters:

        - ``spokes`` -- integer (default: 16); the number of equally
          spaced radial lines to plot

        - ``circles`` -- integer (default: 4); the number of equally
          spaced circles about the center to plot

        - ``pts`` -- integer (default: 32); the number of points to
          plot. Each radial line is made by ``1*pts`` points, each circle
          has ``2*pts`` points. Note that high values may cause erratic
          behavior of the radial lines near the boundaries.
          - only for simply connected domains

        - ``linescale`` -- float between 0 and 1; shrinks the radial lines
          away from the boundary to reduce erratic behavior
          - only for simply connected domains

        - ``rgbcolor`` -- float array (default: ``[0,0,0]``); the
          red-green-blue color of the spiderweb

        - ``thickness`` -- positive float (default: 1); the thickness of
          the lines or points in the spiderweb

        - ``plotjoined`` -- boolean (default: ``True``); if ``False``,
          discrete points will be drawn; otherwise they will be connected
          by lines
          - only for simply connected domains

        - ``withcolor`` -- boolean (default: ``False``); if ``True``,
          the spiderweb will be overlaid on the basic color plot

        - ``plot_points`` -- integer (default: 200); the size of the grid
          in the x direction. The number of points in the y_direction is scaled
          accordingly. Note that very large values can cause this function to
          run slowly.
          - only for multiply connected domains

        - ``min_mag`` -- float (default: 0.001); the magnitude cutoff
          below which spiderweb points are not drawn. This only applies
          to multiply connected domains and is designed to prevent
          "fuzz" at the edge of the domain. Some complicated multiply
          connected domains (particularly those with corners) may
          require a larger value to look clean outside.

        EXAMPLES:

        General usage::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)

        Default plot::

            sage: m.plot_spiderweb()                                                    # needs sage.plot
            Graphics object consisting of 21 graphics primitives

        Simplified plot with many discrete points::

            sage: m.plot_spiderweb(spokes=4, circles=1, pts=400,                        # needs sage.plot
            ....:                  linescale=0.95, plotjoined=False)
            Graphics object consisting of 6 graphics primitives

        Plot with thick, red lines::

            sage: m.plot_spiderweb(rgbcolor=[1,0,0], thickness=3)                       # needs sage.plot
            Graphics object consisting of 21 graphics primitives

        To generate the unit circle map, it's helpful to see what the
        original spiderweb looks like::

            sage: f(t) = e^(I*t)
            sage: fprime(t) = I*e^(I*t)
            sage: m = Riemann_Map([f], [fprime], 0, 1000)
            sage: m.plot_spiderweb()                                                    # needs sage.plot
            Graphics object consisting of 21 graphics primitives

        A multiply connected region with corners. We set ``min_mag`` higher
        to remove "fuzz" outside the domain::

            sage: ps = polygon_spline([(-4,-2),(4,-2),(4,2),(-4,2)])
            sage: z1 = lambda t: ps.value(t); z1p = lambda t: ps.derivative(t)
            sage: z2(t) = -2+exp(-I*t); z2p(t) = -I*exp(-I*t)
            sage: z3(t) = 2+exp(-I*t); z3p(t) = -I*exp(-I*t)
            sage: m = Riemann_Map([z1,z2,z3], [z1p,z2p,z3p], 0,             # long time
            ....:                 ncorners=4)
            sage: p = m.plot_spiderweb(withcolor=True, plot_points=500,         # long time, needs sage.plot
            ....:                      thickness=2.0, min_mag=0.1)
        """
        from sage.plot.complex_plot import ComplexPlot
        from sage.plot.all import list_plot, Graphics

        cdef int k, i
        if self.exterior:
            raise ValueError(
                "Spiderwebs for exterior maps are not currently    supported")
        if self.B == 1:  # The efficient simply connected
            edge = self.plot_boundaries(plotjoined=plotjoined,
                rgbcolor=rgbcolor, thickness=thickness)
            circle_list = list(range(circles))
            theta_array = self.theta_array[0]
            s = spline(np.column_stack([self.theta_array[0], self.tk2]).tolist())
            tmax = self.theta_array[0, self.N]
            tmin = self.theta_array[0, 0]
            for k in range(circles):
                temp = list(range(pts*2))
                for i in range(2*pts):
                    temp[i] = self.inverse_riemann_map(
                        (k + 1) / (circles + 1.0) * exp(I*i * TWOPI / (2*pts)))
                if plotjoined:
                    circle_list[k] = list_plot(comp_pt(temp, 1),
                        rgbcolor=rgbcolor, thickness=thickness, plotjoined=True)
                else:
                    circle_list[k] = list_plot(comp_pt(temp, 1),
                        rgbcolor=rgbcolor, pointsize=thickness)
            line_list = list(range(spokes))
            for k in range(spokes):
                temp = list(range(pts))
                angle = (k*1.0) / spokes * TWOPI
                if angle >= tmax:
                    angle -= TWOPI
                elif angle <= tmin:
                    angle += TWOPI
                for i in range(pts - 1):
                    temp[i] = self.inverse_riemann_map(
                        (i * 1.0) / (pts * 1.0) * exp(I * angle) * linescale)
                temp[pts - 1] = complex(
                    self.f(s(angle)) if angle <= tmax else self.f(s(angle-TWOPI)))
                if plotjoined:
                    line_list[k] = list_plot(
                        comp_pt(temp, 0), rgbcolor=rgbcolor, thickness=thickness,
                        plotjoined=True)
                else:
                    line_list[k] = list_plot(
                        comp_pt(temp, 0), rgbcolor=rgbcolor, pointsize=thickness)
            if withcolor:
                return edge + sum(circle_list) + sum(line_list) + \
                    self.plot_colored(plot_points=plot_points)
            else:
                return edge + sum(circle_list) + sum(line_list)
        else: # The more difficult multiply connected
            z_values, xmin, xmax, ymin, ymax = self.compute_on_grid([],
                plot_points)
            xstep = (xmax-xmin)/plot_points
            ystep = (ymax-ymin)/plot_points
            dr, dtheta= get_derivatives(z_values, xstep, ystep) # clean later

            g = Graphics()
            g.add_primitive(ComplexPlot(complex_to_spiderweb(z_values,dr,dtheta,
                spokes, circles, rgbcolor,thickness, withcolor, min_mag),
                (xmin, xmax), (ymin, ymax),options))
            return g + self.plot_boundaries(thickness = thickness)

    @options(interpolation='catrom')
    def plot_colored(self, plot_range=[], int plot_points=100, **options):
        """
        Generates a colored plot of the Riemann map. A red point on the
        colored plot corresponds to a red point on the unit disc.

        INPUT:

        The following inputs may be passed in as named parameters:

        - ``plot_range`` -- (default: ``[]``) list of 4 values
          ``(xmin, xmax, ymin, ymax)``. Declare if you do not want the plot
          to use the default range for the figure.

        - ``plot_points`` -- integer (default: 100); number of points to
          plot in the x direction. Points in the y direction are scaled
          accordingly. Note that very large values can cause this function to
          run slowly.

        EXAMPLES:

        Given a Riemann map m, general usage::

            sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
            sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
            sage: m = Riemann_Map([f], [fprime], 0)
            sage: m.plot_colored()                                                      # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        Plot zoomed in on a specific spot::

            sage: m.plot_colored(plot_range=[0,1,.25,.75])                              # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        High resolution plot::

            sage: m.plot_colored(plot_points=1000)      # long time (29s on sage.math, 2012), needs sage.plot
            Graphics object consisting of 1 graphics primitive

        To generate the unit circle map, it's helpful to see what the
        colors correspond to::

            sage: f(t) = e^(I*t)
            sage: fprime(t) = I*e^(I*t)
            sage: m = Riemann_Map([f], [fprime], 0, 1000)
            sage: m.plot_colored()                                                      # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.complex_plot import ComplexPlot
        from sage.plot.all import Graphics

        z_values, xmin, xmax, ymin, ymax = self.compute_on_grid(plot_range,
            plot_points)
        g = Graphics()
        g.add_primitive(ComplexPlot(complex_to_rgb(z_values), (xmin, xmax),
            (ymin, ymax),options))
        return g

cdef comp_pt(clist, loop=True):
    """
    Utility function to convert the list of complex numbers
    ``xderivs = get_derivatives(z_values, xstep, ystep)[0]`` to the plottable
    `(x,y)` form. If ``loop=True``, then the first point will be
    added as the last, i.e. to plot a closed circle.

    INPUT:

    - ``clist`` -- list of complex numbers

    - ``loop`` -- boolean (default: ``True``); controls whether or not the
      first point will be added as the last to plot a closed circle

    EXAMPLES:

    This tests it implicitly::

        sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
        sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
        sage: m = Riemann_Map([f], [fprime], 0)
        sage: m.plot_spiderweb()                                                        # needs sage.plot
        Graphics object consisting of 21 graphics primitives
    """
    list2 = [(c.real, c.imag) for c in clist]
    if loop:
        list2.append(list2[0])
    return list2

cpdef get_derivatives(np.ndarray[COMPLEX_T, ndim=2] z_values, FLOAT_T xstep,
    FLOAT_T ystep):
    """
    Compute the r*e^(I*theta) form of derivatives from the grid of points. The
    derivatives are computed using quick-and-dirty taylor expansion and
    assuming analyticity. As such ``get_derivatives`` is primarily intended
    to be used for comparisons in ``plot_spiderweb`` and not for
    applications that require great precision.

    INPUT:

    - ``z_values`` -- the values for a complex function evaluated on a grid
      in the complex plane, usually from ``compute_on_grid``

    - ``xstep`` -- float, the spacing of the grid points in the real direction

    OUTPUT:

    - A tuple of arrays, [``dr``, ``dtheta``], with each array 2 less in both
      dimensions than ``z_values``

      - ``dr`` -- the abs of the derivative of the function in the +r direction
      - ``dtheta`` -- the rate of accumulation of angle in the +theta direction

    EXAMPLES:

    Standard usage with compute_on_grid::

        sage: from sage.calculus.riemann import get_derivatives
        sage: f(t) = e^(I*t) - 0.5*e^(-I*t)
        sage: fprime(t) = I*e^(I*t) + 0.5*I*e^(-I*t)
        sage: m = Riemann_Map([f], [fprime], 0)
        sage: data = m.compute_on_grid([],19)
        sage: xstep = (data[2]-data[1])/19
        sage: ystep = (data[4]-data[3])/19
        sage: dr, dtheta = get_derivatives(data[0],xstep,ystep)
        sage: dr[8,8]
        0.241...
        sage: dtheta[5,5]
        5.907...
    """
    cdef np.ndarray[COMPLEX_T, ndim=2] xderiv
    cdef np.ndarray[FLOAT_T, ndim = 2] dr, dtheta, zabs
    imax = len(z_values)-2
    jmax = len(z_values[0])-2
    #(f(x+delta)-f(x-delta))/2delta
    xderiv = (z_values[1:-1,2:]-z_values[1:-1,:-2])/(2*xstep)
    #b/c the function is analytic, we know the magnitude of its
    #derivative is equal in all directions
    dr = np.abs(xderiv)
    # the abs(derivative) scaled by distance from origin
    zabs = np.abs(z_values[1:-1,1:-1])
    dtheta = np.divide(dr,zabs)
    return dr, dtheta

cpdef complex_to_spiderweb(np.ndarray[COMPLEX_T, ndim = 2] z_values,
    np.ndarray[FLOAT_T, ndim = 2] dr, np.ndarray[FLOAT_T, ndim = 2] dtheta,
    spokes, circles, rgbcolor, thickness, withcolor, min_mag):
    """
    Convert a grid of complex numbers into a matrix containing rgb data
    for the Riemann spiderweb plot.

    INPUT:

    - ``z_values`` -- a grid of complex numbers, as a list of lists

    - ``dr`` -- grid of floats, the r derivative of ``z_values``
      Used to determine precision

    - ``dtheta`` -- grid of floats, the theta derivative of ``z_values``
      Used to determine precision

    - ``spokes`` -- integer; the number of equally spaced radial lines to plot

    - ``circles`` -- integer; the number of equally spaced circles about the
      center to plot

    - ``rgbcolor`` -- float array; the red-green-blue color of the
      lines of the spiderweb

    - ``thickness`` -- positive float; the thickness of the lines or points
      in the spiderweb

    - ``withcolor`` -- boolean; if ``True`` the spiderweb will be overlaid
      on the basic color plot

    - ``min_mag`` -- float; the magnitude cutoff below which spiderweb
      points are not drawn. This only applies to multiply connected
      domains and is designed to prevent "fuzz" at the edge of the
      domain. Some complicated multiply connected domains (particularly
      those with corners) may require a larger value to look clean
      outside.

    OUTPUT:

    An `N x M x 3` floating point Numpy array ``X``, where
    ``X[i,j]`` is an (r,g,b) tuple.

    EXAMPLES::

        sage: from sage.calculus.riemann import complex_to_spiderweb
        sage: import numpy
        sage: zval = numpy.array([[0,1,1000], [.2+.3j,1,-.3j], [0,0,0]],
        ....:                    dtype=numpy.complex128)
        sage: deriv = numpy.array([[.1]],dtype = numpy.float64)
        sage: complex_to_spiderweb(zval, deriv, deriv, 4, 4, [0,0,0], 1, False, 0.001)
        array([[[1., 1., 1.],
                [1., 1., 1.],
                [1., 1., 1.]],
        <BLANKLINE>
               [[1., 1., 1.],
                [0., 0., 0.],
                [1., 1., 1.]],
        <BLANKLINE>
               [[1., 1., 1.],
                [1., 1., 1.],
                [1., 1., 1.]]])

        sage: complex_to_spiderweb(zval, deriv, deriv, 4, 4, [0,0,0], 1, True, 0.001)
        array([[[1.        , 1.        , 1.        ],
                [1.        , 0.05558355, 0.05558355],
                [0.17301243, 0.        , 0.        ]],
        <BLANKLINE>
               [[1.        , 0.96804683, 0.48044583],
                [0.        , 0.        , 0.        ],
                [0.77351965, 0.5470393 , 1.        ]],
        <BLANKLINE>
               [[1.        , 1.        , 1.        ],
                [1.        , 1.        , 1.        ],
                [1.        , 1.        , 1.        ]]])
     """
    cdef Py_ssize_t i, j, imax, jmax
    cdef FLOAT_T x, y, mag, arg, width, target, precision, dmag, darg
    cdef COMPLEX_T z
    cdef FLOAT_T DMAX = 70 # change to adjust rate_of_change cutoff below
    precision = thickness/150.0
    imax = len(z_values)
    jmax = len(z_values[0])
    cdef np.ndarray[FLOAT_T, ndim=3, mode="c"] rgb
    if withcolor:
        rgb = complex_to_rgb(z_values)
    else:
        rgb = np.zeros(dtype=FLOAT, shape=(imax, jmax, 3))
        rgb += 1
    if circles != 0:
        circ_radii = srange(0,1.0,1.0/circles)
    else:
        circ_radii = []
    if spokes != 0:
        # both -pi and pi are included
        spoke_angles = srange(-PI,PI+TWOPI/spokes,TWOPI/spokes)
    else:
        spoke_angles = []
    for i in range(imax-2): # the d arrays are 1 smaller on each side
        for j in range(jmax-2):
            z = z_values[i+1,j+1]
            mag = abs(z)
            arg = phase(z)
            dmag = dr[i,j]
            darg = dtheta[i,j]
            #points that change too rapidly are presumed to be borders
            #points that are too small are presumed to be outside
            if darg < DMAX and mag > min_mag:
                for target in circ_radii:
                    if abs(mag - target)/dmag < precision:
                        rgb[i+1,j+1] = rgbcolor
                        break
                for target in spoke_angles:
                    if abs(arg - target)/darg < precision:
                        rgb[i+1,j+1] = rgbcolor
                        break
    return rgb


cpdef complex_to_rgb(np.ndarray[COMPLEX_T, ndim = 2] z_values):
    r"""
    Convert from a (Numpy) array of complex numbers to its corresponding
    matrix of RGB values.  For internal use of :meth:`~Riemann_Map.plot_colored`
    only.

    INPUT:

    - ``z_values`` -- numpy array of complex numbers

    OUTPUT:

    An `N \times M \times 3` floating point Numpy array ``X``, where
    ``X[i,j]`` is an (r,g,b) tuple.

    EXAMPLES::

        sage: from sage.calculus.riemann import complex_to_rgb
        sage: import numpy
        sage: complex_to_rgb(numpy.array([[0, 1, 1000]], dtype=numpy.complex128))
        array([[[1.        , 1.        , 1.        ],
                [1.        , 0.05558355, 0.05558355],
                [0.17301243, 0.        , 0.        ]]])

        sage: complex_to_rgb(numpy.array([[0, 1j, 1000j]], dtype=numpy.complex128))
        array([[[1.        , 1.        , 1.        ],
                [0.52779177, 1.        , 0.05558355],
                [0.08650622, 0.17301243, 0.        ]]])


    TESTS::

        sage: complex_to_rgb([[0, 1, 10]])
        Traceback (most recent call last):
        ...
        TypeError: Argument 'z_values' has incorrect type (expected numpy.ndarray, got list)
    """
    cdef Py_ssize_t i, j, imax, jmax
    cdef FLOAT_T x, y, mag, arg
    cdef FLOAT_T lightness, hue, top, bot
    cdef FLOAT_T r, g, b
    cdef int ihue
    cdef COMPLEX_T z

    imax = len(z_values)
    jmax = len(z_values[0])
    cdef np.ndarray[FLOAT_T, ndim=3, mode="c"] rgb = np.empty(
        dtype=FLOAT, shape=(imax, jmax, 3))

    sig_on()
    for i in range(imax):
        row = z_values[i]
        for j in range(jmax):
            z = row[j]
            mag = abs(z)
            arg = phase(z)
            # tweak these levels to adjust how bright/dark the colors appear
            # output can range from -1 (black) to 1 (white)
            lightness = -(atan(log(mag*1.5+1)) * (4/PI) - 1)
            # in hsv, variable value, full saturation (s=1, v=1+lightness)
            if lightness < 0:
                bot = 0
                top = (1 + lightness)
            # in hsv, variable saturation, full value (v=1, s=1-lightness)
            else:
                bot = lightness
                top = 1
            # Note that does same thing as colorsys module hsv_to_rgb for
            # this setup, but in Cython.
            hue = 3*arg / PI
            # usual hsv hue is thus h=arg/(2*pi) for positive,
            # h=arg/(2*PI)+1 for negative
            if hue < 0:
                hue += 6
            ihue = <int>hue
            if ihue == 0:
                r = top
                g = bot + hue * (top - bot)
                b = bot
            elif ihue == 1:
                r = bot + (2 - hue) * (top - bot)
                g = top
                b = bot
            elif ihue == 2:
                r = bot
                g = top
                b = bot + (hue - 2) * (top - bot)
            elif ihue == 3:
                r = bot
                g = bot + (4 - hue) * (top - bot)
                b = top
            elif ihue == 4:
                r = bot + (hue - 4) * (top - bot)
                g = bot
                b = top
            else:
                r = top
                g = bot
                b = bot + (6 - hue) * (top - bot)
            rgb[i, j, 0] = r
            rgb[i, j, 1] = g
            rgb[i, j, 2] = b
    sig_off()
    return rgb

cpdef analytic_boundary(FLOAT_T t, int n, FLOAT_T epsilon):
    r"""
    Provides an exact (for `n = \infty`) Riemann boundary
    correspondence for the ellipse with axes `1 + \epsilon` and `1 - \epsilon`.
    The boundary is therefore given by `\exp(I t)+\epsilon\exp(-I t)`. It is
    primarily useful for testing the accuracy of the numerical
    :class:`Riemann_Map`.

    INPUT:

    - ``t`` -- the boundary parameter, from `0` to `2 \pi`

    - ``n`` -- integer; the number of terms to include
      (10 is fairly accurate, 20 is very accurate)

    - ``epsilon`` -- float; the skew of the ellipse (0 is circular)

    OUTPUT:

    A theta value from 0 to 2*pi, corresponding to the point on the
    circle e^(I*theta)

    TESTS:

    Checking the accuracy of this function for different n values::

        sage: from sage.calculus.riemann import analytic_boundary
        sage: t100 = analytic_boundary(pi/2, 100, .3)
        sage: abs(analytic_boundary(pi/2, 10, .3) - t100) < 10^-8
        True
        sage: abs(analytic_boundary(pi/2, 20, .3) - t100) < 10^-15
        True

    Using this to check the accuracy of the Riemann_Map boundary::

        sage: f(t) = e^(I*t)+.3*e^(-I*t)
        sage: fp(t) = I*e^(I*t)-I*.3*e^(-I*t)
        sage: m = Riemann_Map([f], [fp],0,200)
        sage: s = spline(m.get_theta_points())
        sage: test_pt = uniform(0,2*pi)
        sage: s(test_pt) - analytic_boundary(test_pt,20, .3) < 10^-3
        True
    """
    cdef FLOAT_T i
    cdef FLOAT_T result = t
    for i in range(1, n + 1):
        result += (2*(-1)**i/i)*(epsilon**i/(1+epsilon**(2*i)))*sin(2*i*t)
    return result


cpdef cauchy_kernel(t, args):
    """
    Intermediate function for the integration in :meth:`~Riemann_Map.analytic_interior`.

    INPUT:

    - ``t`` -- the boundary parameter, meant to be integrated over

    - ``args`` -- tuple containing:

      - ``epsilon`` -- float; the skew of the ellipse (0 is circular)

      - ``z`` -- complex; the point to be mapped

      - ``n`` -- integer; the number of terms to include
        (10 is fairly accurate, 20 is very accurate)

      - ``part`` -- will return the real ('r'), imaginary ('i') or
        complex ('c') value of the kernel

    TESTS:

    This is primarily tested implicitly by :meth:`~Riemann_Map.analytic_interior`.
    Here is a simple test::

        sage: from sage.calculus.riemann import cauchy_kernel
        sage: cauchy_kernel(.5,(.3, .1+.2*I, 10,'c'))
        (-0.584136405997...+0.5948650858950...j)
    """
    cdef COMPLEX_T result
    cdef FLOAT_T epsilon = args[0]
    cdef COMPLEX_T z = args[1]
    cdef int n = args[2]
    part = args[3]
    result = exp(I*analytic_boundary(t,n, epsilon))/(exp(I*t)+epsilon*exp(-I*t)-z) *  \
        (I*exp(I*t)-I*epsilon*exp(-I*t))
    if part == 'c':
        return result
    elif part == 'r':
        return result.real
    elif part == 'i':
        return result.imag
    else:
        return None


cpdef analytic_interior(COMPLEX_T z, int n, FLOAT_T epsilon):
    """
    Provides a nearly exact computation of the Riemann Map of an interior
    point of the ellipse with axes 1 + epsilon and 1 - epsilon. It is
    primarily useful for testing the accuracy of the numerical Riemann Map.

    INPUT:

    - ``z`` -- complex; the point to be mapped

    - ``n`` -- integer; the number of terms to include
      (10 is fairly accurate, 20 is very accurate)

    TESTS:

    Testing the accuracy of :class:`Riemann_Map`::

        sage: from sage.calculus.riemann import analytic_interior
        sage: f(t) = e^(I*t)+.3*e^(-I*t)
        sage: fp(t) = I*e^(I*t)-I*.3*e^(-I*t)
        sage: m = Riemann_Map([f],[fp],0,200)
        sage: abs(m.riemann_map(.5)-analytic_interior(.5, 20, .3)) < 10^-4
        True
        sage: m = Riemann_Map([f],[fp],0,2000)
        sage: abs(m.riemann_map(.5)-analytic_interior(.5, 20, .3)) < 10^-6
        True
    """
    # evaluates the Cauchy integral of the boundary, split into the real
    # and imaginary results because numerical_integral can't handle complex data.
    rp = 1/(TWOPI)*numerical_integral(cauchy_kernel,0,2*pi,
        params = [epsilon,z,n,'i'])[0]
    ip = 1/(TWOPI*I)*numerical_integral(cauchy_kernel,0,2*pi,
        params = [epsilon,z,n,'r'])[0]
    return rp + ip
