r"""
Class for computing sums over zeros of motivic `L`-functions

All computations are done to double precision.

AUTHORS:

- Simon Spicer (2014-10): first version
"""
# ****************************************************************************
#       Copyright (C) 2014 Simon Spicer <mlungu@uw.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import next_prime
from sage.functions.log import log, exp
from sage.functions.other import real, imag
from sage.libs.flint.ulong_extras cimport n_is_prime
from sage.misc.verbose import verbose
from sage.parallel.decorate import parallel
from sage.parallel.ncpus import ncpus as num_cpus
from sage.rings.complex_double import CDF
from sage.rings.infinity import PlusInfinity
from sage.rings.integer_ring import ZZ
from sage.rings.real_double import RDF
from sage.structure.sage_object cimport SageObject
from sage.symbolic.constants import euler_gamma, pi


cdef extern from "<math.h>":
    double c_exp "exp"(double)
    double c_log "log"(double)
    double c_cos "cos"(double)
    double c_acos "acos"(double)
    double c_sqrt "sqrt"(double)

# Global variable determining the number of CPUs to use for parallel computations
cdef NCPUS

cdef class LFunctionZeroSum_abstract(SageObject):
    r"""
    Abstract class for computing certain sums over zeros of a motivic `L`-function
    without having to determine the zeros themselves.
    """
    cdef _pi            # Pi to 64 bits
    cdef _euler_gamma   # Euler-Mascheroni constant = 0.5772...
    cdef _level         # The level of the form attached to self
    cdef _k             # The weight of the form attached to self
    cdef _C1            # = log(N)/2 - log(2*pi)
    cdef _C0            # = C1 - euler_gamma
    cdef _ncpus         # The number of CPUs to use for parallel computations

    def ncpus(self, n=None):
        r"""
        Set or return the number of CPUs to be used in parallel computations.

        If called with no input, the number of CPUs currently set is returned;
        else this value is set to `n`. If `n` is 0 then the number of CPUs is set
        to the max available.

        INPUT:

        - ``n`` -- (default: ``None``) if not ``None``, a nonnegative integer

        OUTPUT: if `n` is not ``None``, returns a positive integer

        EXAMPLES::

            sage: Z = LFunctionZeroSum(EllipticCurve("389a"))
            sage: Z.ncpus()
            1
            sage: Z.ncpus(2)
            sage: Z.ncpus()
            2

        The following output will depend on the system that Sage is running on.

        ::

            sage: Z.ncpus(0)
            sage: Z.ncpus()            # random
            4
        """
        if n is None:
            return self._ncpus
        elif n < 0:
            raise ValueError("input must be positive integer")
        elif n == 0:
            self._ncpus = num_cpus()
            NCPUS = self._ncpus
        else:
            self._ncpus = n
            NCPUS = self._ncpus

    def level(self):
        r"""
        Return the level of the form attached to ``self``.

        If ``self`` was constructed from an elliptic curve,
        then this is equal to the conductor of `E`.

        EXAMPLES::

            sage: E = EllipticCurve("389a")
            sage: Z = LFunctionZeroSum(E)
            sage: Z.level()
            389
        """
        return self._level

    def weight(self):
        r"""
        Return the weight of the form attached to ``self``.

        If ``self`` was constructed
        from an elliptic curve, then this is 2.

        EXAMPLES::

            sage: E = EllipticCurve("389a")
            sage: Z = LFunctionZeroSum(E)
            sage: Z.weight()
            2
        """
        return self._k

    def C0(self, include_euler_gamma=True):
        r"""
        Return the constant term of the logarithmic derivative of the
        completed `L`-function attached to ``self``.

        This is equal to
        `-\eta + \log(N)/2 - \log(2\pi)`, where `\eta` is the
        Euler-Mascheroni constant `= 0.5772...`
        and `N` is the level of the form attached to ``self``.

        INPUT:

        - ``include_euler_gamma`` -- boolean (default: ``True``); if set to
          ``False``, return the constant `\log(N)/2 - \log(2\pi)`, i.e., do
          not subtract the Euler-Mascheroni constant

        EXAMPLES::

            sage: E = EllipticCurve("389a")
            sage: Z = LFunctionZeroSum(E)
            sage: Z.C0() # tol 1.0e-13
            0.5666969404983447
            sage: Z.C0(include_euler_gamma=False) # tol 1.0e-13
            1.1439126053998776
        """
        # Computed at initialization
        return self._C1 if not include_euler_gamma else self._C0

    def cnlist(self, n, python_floats=False):
        r"""
        Return a list of Dirichlet coefficient of the logarithmic
        derivative of the `L`-function attached to ``self``, shifted so that
        the critical line lies on the imaginary axis, up to and
        including `n`.

        The `i`-th element of the returned list is ``a[i]``.

        INPUT:

        - ``n`` -- nonnegative integer

        - ``python_floats`` -- boolean (default: ``False``); if ``True`` return
          a list of Python floats instead of Sage Real Double Field elements

        OUTPUT: list of real numbers

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_EllipticCurve.cn`

        .. TODO::

            Speed this up; make more efficient

        EXAMPLES::

            sage: E = EllipticCurve("11a")
            sage: Z = LFunctionZeroSum(E)
            sage: cnlist = Z.cnlist(11)
            sage: for n in range(12): print((n, cnlist[n])) # tol 1.0e-13
            (0, 0.0)
            (1, 0.0)
            (2, 0.6931471805599453)
            (3, 0.3662040962227033)
            (4, 0.0)
            (5, -0.32188758248682003)
            (6, 0.0)
            (7, 0.555974328301518)
            (8, -0.34657359027997264)
            (9, 0.6103401603711721)
            (10, 0.0)
            (11, -0.21799047934530644)
        """
        if not python_floats:
            return [self.cn(i) for i in range(n + 1)]
        return [float(self.cn(i)) for i in range(n + 1)]

    def digamma(self, s, include_constant_term=True):
        r"""
        Return the digamma function `\digamma(s)` on the complex input s.

        This is given by
        `\digamma(s) = -\eta + \sum_{k=1}^{\infty} \frac{s-1}{k(k+s-1)}`,
        where `\eta` is the Euler-Mascheroni constant `=0.5772156649\ldots`.

        This function is needed in the computing the logarithmic derivative
        of the `L`-function attached to ``self``.

        INPUT:

        - ``s`` -- complex number

        - ``include_constant_term`` -- boolean (default: ``True``); if set
          to ``False``, only the value of the sum over `k` is returned without
          subtracting the Euler-Mascheroni constant, i.e., the returned value
          is equal to `\sum_{k=1}^{\infty} \frac{s-1}{k(k+s-1)}`

        OUTPUT:

        A real double precision number if the input is real and not a negative
        integer; Infinity if the input is a negative integer, and a complex
        number otherwise.

        EXAMPLES::

            sage: Z = LFunctionZeroSum(EllipticCurve("37a"))
            sage: Z.digamma(3.2) # tol 1.0e-13
            0.9988388912865993
            sage: Z.digamma(3.2,include_constant_term=False) # tol 1.0e-13
            1.576054556188132
            sage: Z.digamma(1+I) # tol 1.0e-13
            0.09465032062247625 + 1.076674047468581*I
            sage: Z.digamma(-2)
            +Infinity

        Evaluating the sum without the constant term at the positive integers n
        returns the (n-1)th harmonic number.

        ::

            sage: Z.digamma(3,include_constant_term=False)
            1.5
            sage: Z.digamma(6,include_constant_term=False)
            2.283333333333333
        """
        # imported here so as to avoid importing Numpy on Sage startup
        from scipy.special import psi

        if real(s) < 0 and imag(s) == 0:
            try:
                z = ZZ(s)
                return PlusInfinity()
            except Exception:
                pass

        if imag(s) == 0:
            F = RDF
        else:
            F = CDF
        # Cheating: SciPy already has this function implemented for complex inputs
        z = F(psi(F(s)))
        if include_constant_term:
            return z
        else:
            return z + self._euler_gamma

    def logarithmic_derivative(self, s, num_terms=10000, as_interval=False):
        r"""
        Compute the value of the logarithmic derivative
        `\frac{L^{\prime}}{L}` at the point s to *low* precision, where `L`
        is the `L`-function attached to ``self``.

        .. WARNING::

            The value is computed naively by evaluating the Dirichlet series
            for `\frac{L^{\prime}}{L}`; convergence is controlled by the
            distance of s from the critical strip `0.5<=\Re(s)<=1.5`.
            You may use this method to attempt to compute values inside the
            critical strip; however, results are then *not* guaranteed
            to be correct to any number of digits.

        INPUT:

        - ``s`` -- real or complex value

        - ``num_terms`` -- integer (default: 10000); the maximum number of
          terms summed in the Dirichlet series

        OUTPUT:

        A tuple (z,err), where z is the computed value, and err is an
        upper bound on the truncation error in this value introduced
        by truncating the Dirichlet sum.

        .. NOTE::

            For the default term cap of 10000, a value accurate to all 53
            bits of a double precision floating point number is only
            guaranteed when `|\Re(s-1)|>4.58`, although in practice inputs
            closer to the critical strip will still yield computed values
            close to the true value.

        EXAMPLES::

            sage: E = EllipticCurve([23,100])
            sage: Z = LFunctionZeroSum(E)
            sage: Z.logarithmic_derivative(10) # tol 1.0e-13
            (5.648066742632698e-05, 1.0974102859764345e-34)
            sage: Z.logarithmic_derivative(2.2) # tol 1.0e-13
            (0.5751257063594758, 0.024087912696974387)

        Increasing the number of terms should see the truncation error
        decrease.

        ::

            sage: Z.logarithmic_derivative(2.2,num_terms=50000) # long time # rel tol 1.0e-14
            (0.5751579645060139, 0.008988775519160675)

        Attempting to compute values inside the critical strip
        gives infinite error.

        ::

            sage: Z.logarithmic_derivative(1.3) # tol 1.0e-13
            (5.442994413920786, +Infinity)

        Complex inputs and inputs to the left of the critical strip
        are allowed.

        ::

            sage: Z.logarithmic_derivative(complex(3,-1)) # tol 1.0e-13
            (0.04764548578052381 + 0.16513832809989326*I, 6.584671359095225e-06)
            sage: Z.logarithmic_derivative(complex(-3,-1.1)) # tol 1.0e-13
            (-13.908452173241546 + 2.591443099074753*I, 2.7131584736258447e-14)

        The logarithmic derivative has poles at the negative integers.

        ::

            sage: Z.logarithmic_derivative(-3) # tol 1.0e-13
            (-Infinity, 2.7131584736258447e-14)
        """
        if imag(s) == 0:
            F = RDF
        else:
            F = CDF
        # Inputs left of the critical line are handled via the functional
        # equation of the logarithmic derivative
        if real(s - 1) < 0:
            a = -2 * self._C1 - self.digamma(s) - self.digamma(2 - s)
            b, err = self.logarithmic_derivative(2 - s, num_terms=num_terms)
            return (a + b, err)

        z = s - 1
        sigma = RDF(real(z))
        # Compute maximum possible Dirichlet series truncation error
        # When s is in the critical strip: no guaranteed precision
        if abs(sigma) <= 0.5:
            err = PlusInfinity()
        else:
            a = RDF(sigma) - RDF(0.5)
            b = log(RDF(num_terms)) * a
            err = (b + 1) * exp(-b) / a**2

        y = F.zero()
        n = ZZ(2)
        while n <= num_terms:
            if n.is_prime_power():
                cn = self.cn(n)
                y += cn / F(n)**z
            n += 1

        return (y, err)

    def completed_logarithmic_derivative(self, s, num_terms=10000):
        r"""
        Compute the value of the completed logarithmic derivative
        `\frac{\Lambda^{\prime}}{\Lambda}` at the point s to *low*
        precision, where `\Lambda = N^{s/2}(2\pi)^s \Gamma(s) L(s)`
        and `L` is the `L`-function attached to ``self``.

        .. WARNING::

            This is computed naively by evaluating the Dirichlet series
            for `\frac{L^{\prime}}{L}`; the convergence thereof is
            controlled by the distance of s from the critical strip
            `0.5<=\Re(s)<=1.5`.
            You may use this method to attempt to compute values inside the
            critical strip; however, results are then *not* guaranteed
            to be correct to any number of digits.

        INPUT:

        - ``s`` -- real or complex value

        - ``num_terms`` -- integer (default: 10000); the maximum number of
          terms summed in the Dirichlet series

        OUTPUT:

        A tuple (z,err), where z is the computed value, and err is an
        upper bound on the truncation error in this value introduced
        by truncating the Dirichlet sum.

        .. NOTE::

            For the default term cap of 10000, a value accurate to all 53
            bits of a double precision floating point number is only
            guaranteed when `|\Re(s-1)|>4.58`, although in practice inputs
            closer to the critical strip will still yield computed values
            close to the true value.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_EllipticCurve.logarithmic_derivative`

        EXAMPLES::

            sage: E = EllipticCurve([23,100])
            sage: Z = LFunctionZeroSum(E)
            sage: Z.completed_logarithmic_derivative(3) # tol 1.0e-13
            (6.64372066048195, 6.584671359095225e-06)

        Complex values are handled. The function is odd about s=1, so
        the value at 2-s should be minus the value at s.

        ::

            sage: Z.completed_logarithmic_derivative(complex(-2.2,1)) # tol 1.0e-13
            (-6.898080633125154 + 0.22557015394248361*I, 5.623853049808912e-11)
            sage: Z.completed_logarithmic_derivative(complex(4.2,-1)) # tol 1.0e-13
            (6.898080633125154 - 0.22557015394248361*I, 5.623853049808912e-11)
        """
        if real(s - 1) >= 0:
            Ls = self.logarithmic_derivative(s, num_terms)
            return (self._C1 + self.digamma(s) + Ls[0], Ls[1])
        else:
            Ls = self.logarithmic_derivative(2 - s, num_terms)
            return (-self._C1 - self.digamma(2 - s) - Ls[0], Ls[1])

    def zerosum(self, Delta=1, tau=0, function='sincsquared_fast', ncpus=None):
        r"""
        Bound from above the analytic rank of the form attached to ``self``.

        This bound is obtained by computing `\sum_{\gamma} f(\Delta(\gamma-\tau))`,
        where `\gamma` ranges over the imaginary parts of the zeros of `L_E(s)`
        along the critical strip, and `f(x)` is an appropriate even continuous
        `L_2` function such that `f(0)=1`.

        If `\tau=0`, then as `\Delta` increases this sum converges from above to
        the analytic rank of the `L`-function, as `f(0) = 1` is counted with
        multiplicity `r`, and the other terms all go to 0 uniformly.

        INPUT:

        - ``Delta`` -- positive real number (default: 1) parameter denoting the
          tightness of the zero sum

        - ``tau`` -- real parameter (default: 0) denoting the offset of the sum
          to be computed. When `\tau=0` the sum will converge to the analytic rank
          of the `L`-function as `\Delta` is increased. If `\tau` is the value
          of the imaginary part of a noncentral zero, the limit will be 1
          (assuming the zero is simple); otherwise, the limit will be 0.
          Currently only implemented for the sincsquared and cauchy functions;
          otherwise ignored.

        - ``function`` -- string (default: ``'sincsquared_fast'``); the function
          `f(x)` as described above. Currently implemented options for `f` are

          - ``sincsquared`` -- `f(x) = \left(\frac{\sin(\pi x)}{\pi x}\right)^2`

          - ``gaussian`` -- `f(x) = e^{-x^2}`

          - ``sincsquared_fast`` -- same as "sincsquared", but implementation
            optimized for elliptic curve `L`-functions, and tau must be 0. self
            must be attached to an elliptic curve over `\QQ` given by its global
            minimal model, otherwise the returned result will be incorrect.

          - ``sincsquared_parallel`` -- same as "sincsquared_fast", but optimized
            for parallel computation with large (>2.0) `\Delta` values. ``self`` must
            be attached to an elliptic curve over `\QQ` given by its global minimal
            model, otherwise the returned result will be incorrect.

          - ``cauchy`` -- `f(x) = \frac{1}{1+x^2}`; this is only computable to
            low precision, and only when `\Delta < 2`

        - ``ncpus`` -- (default: ``None``) if not ``None``, a positive integer
          defining the number of CPUs to be used for the computation. If left as
          ``None``, the maximum available number of CPUs will be used.
          Only implemented for algorithm="sincsquared_parallel"; ignored
          otherwise.

        .. WARNING::

            Computation time is exponential in `\Delta`, roughly doubling for
            every increase of 0.1 thereof. Using `\Delta=1` will yield a
            computation time of a few milliseconds; `\Delta=2` takes a few
            seconds, and `\Delta=3` takes upwards of an hour. Increase at your
            own risk beyond this!

        OUTPUT:

        A positive real number that bounds from above the number of zeros with
        imaginary part equal to `\tau`. When `\tau=0` this is an upper bound for
        the `L`-function's analytic rank.

        .. SEEALSO::

            :meth:`~sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field.analytic_rank_bound`
            for more documentation and examples on calling this method on elliptic curve
            `L`-functions.

        EXAMPLES::

            sage: E = EllipticCurve("389a"); E.rank()
            2
            sage: Z = LFunctionZeroSum(E)
            sage: E.lseries().zeros(3)
            [0.000000000, 0.000000000, 2.87609907]
            sage: Z.zerosum(Delta=1,function='sincsquared_fast') # tol 1.0e-13
            2.037500084595065
            sage: Z.zerosum(Delta=1,function='sincsquared_parallel') # tol 1.0e-11
            2.037500084595065
            sage: Z.zerosum(Delta=1,function='sincsquared') # tol 1.0e-13
            2.0375000845950644
            sage: Z.zerosum(Delta=1,tau=2.876,function='sincsquared') # tol 1.0e-13
            1.075551295651154
            sage: Z.zerosum(Delta=1,tau=1.2,function='sincsquared') # tol 1.0e-13
            0.10831555377490683
            sage: Z.zerosum(Delta=1,function='gaussian') # tol 1.0e-13
            2.056890425029435
        """
        # If Delta>6.95, then exp(2*pi*Delta)>sys.maxsize, so we get overflow
        # when summing over the logarithmic derivative coefficients
        if Delta > 6.95:
            raise ValueError("Delta value too large; will result in overflow")

        if function == "sincsquared_parallel":
            return self._zerosum_sincsquared_parallel(Delta=Delta, ncpus=ncpus)
        elif function == "sincsquared_fast":
            return self._zerosum_sincsquared_fast(Delta=Delta)
        elif function == "sincsquared":
            return self._zerosum_sincsquared(Delta=Delta, tau=tau)
        elif function == "gaussian":
            return self._zerosum_gaussian(Delta=Delta)
        elif function == "cauchy":
            return self._zerosum_cauchy(Delta=Delta, tau=tau)
        else:
            raise ValueError("input function not recognized")

    def _zerosum_sincsquared(self, Delta=1, tau=0):
        r"""
        Bound from above the analytic rank of the form attached to ``self``.

        This bound is obtained by computing `\sum_{\gamma} f(\Delta \cdot (\gamma-\tau))`,
        where `\gamma` ranges over the imaginary parts of the zeros of `L_E(s)`
        along the critical strip, and `f(x) = \sin(\pi x)/(\pi x)`

        If `\tau=0`, then as `\Delta` increases this sum limits from above to
        the analytic rank of the `L`-function, as `f(0) = 1` is counted with
        multiplicity `r`, and the other terms all go to 0 uniformly.

        INPUT:

        - ``Delta`` -- positive real number (default: 1) parameter denoting the
          tightness of the zero sum

        - ``tau`` -- real parameter (default: 0) denoting the offset of the sum
          to be computed. When tau=0 the sum will converge from above to the
          analytic rank of the `L`-function as `\Delta` is increased. If tau
          is the value of the imaginary part of a noncentral zero, the limit
          will be 1 (assuming GRH, the zero is simple); otherwise the limit
          will be 0.

        .. WARNING::

            Computation time is exponential in `\Delta`, roughly doubling for
            every increase of 0.1 thereof. Using `\Delta=1` will yield a
            computation time of a few milliseconds; `\Delta=2` takes a few
            seconds, and `\Delta=3` takes upwards of an hour. Increase at your
            own risk beyond this!

        .. WARNING::

            This method has *not* been optimized with Cython; as such
            computing with this method is slower than the central sum
            (i.e `\tau=0`) versions of self._zerosum_sincsquared_fast() and
            self._zerosum_sincsquared_parallel().

        OUTPUT:

        A positive real number that bounds from above the number of zeros with
        imaginary part equal to tau. When tau=0 this is an upper bound for the
        `L`-function's analytic rank.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve("37a"); E.rank()
            1
            sage: Z = LFunctionZeroSum(E)
            sage: E.lseries().zeros(2)
            [0.000000000, 5.00317001]

        E is a rank 1 curve; the lowest noncentral zero has imaginary part
        ~5.003. The zero sum with tau=0 indicates the probable existence of
        a zero at or very close to the central point.

        ::

            sage: Z._zerosum_sincsquared(Delta=1,tau=0) # tol 1.0e-13
            1.0103840698356257

        The zero sum also detects a zero at or near 5.003, as expected.

        ::

            sage: Z._zerosum_sincsquared(Delta=1,tau=5.003) # tol 1.0e-13
            1.0168124546878288

        However, there is definitely no zero with imaginary part near 2.5,
        as the sum would have to be at least 1.

        ::

            sage: Z._zerosum_sincsquared(Delta=1,tau=2.5) # tol 1.0e-13
            0.058058210806477814
        """
        npi = self._pi
        twopi = 2 * npi

        t = RDF(Delta * twopi)
        expt = RDF(exp(t))

        u = t * self.C0()

        # No offset: formulae are simpler
        if tau == 0:
            w = npi**2 / 6 - (RDF.one() / expt).dilog()

            y = RDF(0)
            n = int(1)
            while n < expt:
                cn = self.cn(n)
                if cn != 0:
                    logn = RDF(n).log()
                    y += cn * (t - logn)
                n += 1
        # When offset is nonzero, the digamma transform (w) must
        # be computed as an infinite sum
        else:
            tau = RDF(tau)
            cos_tau_t = (tau * t).cos()
            sin_tau_t = (tau * t).sin()
            w = RDF.zero()
            for k in range(1, 1001):
                a1 = tau**2 / (k * (k**2 + tau**2))
                a2 = (k**2 - tau**2) / (k**2 + tau**2)**2
                a3 = (2 * k * tau) / (k**2 + tau**2)**2

                w0 = a1 * t + a2
                w0 -= (a2 * cos_tau_t - a3 * sin_tau_t) * exp(-k * t)
                w += w0

            y = RDF.zero()
            n = int(1)
            while n < expt:
                cn = self.cn(n)
                if cn != 0:
                    logn = RDF(n).log()
                    y += cn * (t - logn) * (tau * logn).cos()
                n += 1

        return (u + w + y) * 2 / (t**2)

    def _zerosum_gaussian(self, Delta=1):
        r"""
        Return an upper bound on the analytic rank of the `L`-series attached
        to ``self`` by computing `\sum_{\gamma} f(\Delta*\gamma)`,
        where `\gamma` ranges over the imaginary parts of the zeros of `L_E(s)`
        along the critical strip, and `f(x) = \exp(-x^2)`.

        As `\Delta` increases this sum limits from above to the analytic rank
        of the form, as `f(0) = 1` is counted with multiplicity `r`, and the
        other terms all go to 0 uniformly.

        INPUT:

        - ``Delta`` -- positive real number (default: 1) parameter defining the
          tightness of the zero sum, and thus the closeness of the returned
          estimate to the actual analytic rank of the form attached to ``self``.

        .. WARNING::

            Computation time is exponential in `\Delta`, roughly doubling for
            every increase of 0.1 thereof. Using `\Delta=1` will yield a
            computation time of a few milliseconds; `\Delta=2` takes a few
            seconds, and `\Delta=3` takes upwards of an hour. Increase at your
            own risk beyond this!

        OUTPUT:

        A positive real number that bounds the analytic rank of the modular form
        attached to ``self`` from above.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve("37a"); E.rank()
            1
            sage: Z = LFunctionZeroSum(E)
            sage: Z._zerosum_gaussian(Delta=1) # tol 1.0e-13
            1.0563950773441664
        """
        # imported here so as to avoid importing Numpy on Sage startup
        from scipy.special import erfcx

        npi = self._pi
        Deltasqrtpi = Delta * npi.sqrt()

        t = RDF(Delta * npi * 2)
        expt = t.exp()

        u = self.C0()

        w = RDF.zero()
        for k in range(1, 1001):
            w += RDF.one() / k - erfcx(Delta * k) * Deltasqrtpi

        y = RDF.zero()
        n = int(1)

        # TO DO: Error analysis to make sure this bound is good enough to
        # avoid non-negligible truncation error
        while n < expt:
            cn = self.cn(n)
            if cn != 0:
                logn = log(RDF(n))
                y += cn * exp(-(logn / (2 * Delta))**2)
            n += 1
        # y is the truncation of an infinite sum, so we must add a value which
        # exceeds the max amount we could have left out.
        # WARNING: Truncation error analysis has *not* been done; the value of
        # 0.1 is empirical.
        return RDF(u + w + y + 0.1) / Deltasqrtpi

    def _zerosum_cauchy(self, Delta=1, tau=0, num_terms=None):
        r"""
        Bound from above the analytic rank of the form attached to ``self``
        by computing `\sum_{\gamma} f(\Delta*(\gamma-\tau))`,
        where `\gamma` ranges over the imaginary parts of the zeros of `L_E(s)`
        along the critical strip, and `f(x) = \frac{1}{1+x^2}`.

        If `\tau=0`, then as `\Delta` increases this sum converges from above
        to the analytic rank of the `L`-function, as `f(0) = 1` is counted
        with multiplicity `r`, and the other terms all go to 0 uniformly.

        INPUT:

        - ``Delta`` -- positive real number (default: 1); parameter denoting the
          tightness of the zero sum

        - ``tau`` -- real parameter (default: 0) denoting the offset of the sum
          to be computed. When tau=0 the sum will converge from above to the
          analytic rank of the `L`-function as `\Delta` is increased. If tau is
          the value of the imaginary part of a noncentral zero, the limit will
          be 1 (assuming GRH, the zero is simple); otherwise the limit will
          be 0.

        - ``num_terms`` -- positive integer (default: ``None``); the number of
          terms computed in the truncated Dirichlet series for the `L`-function
          attached to ``self``. If left at ``None``, this is set to
          `\ceil(e^{2 \pi \Delta})`, the same number of terms used in the other
          zero sum methods for this value of Delta.
          Increase num_terms to get more accuracy.

        .. WARNING::

            This value can only be provably computed when Delta < 2; an error
            will be thrown if a Delta value larger than 2 is supplied.
            Furthermore, beware that computation time is exponential in
            `\Delta`, roughly doubling for every increase of 0.1 thereof.
            Using `\Delta=1` will yield a computation time of a few
            milliseconds, while `\Delta=2` takes a few seconds.

        OUTPUT:

        A positive real number that bounds from above the number of zeros with
        imaginary part equal to tau. When tau=0 this is an upper bound for the
        `L`-function's analytic rank.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve("11a")
            sage: zeros = E.lseries().zeros(2)
            sage: zeros[0] # abs tol 1e-8
            6.36261389
            sage: zeros[1] # abs tol 1e-8
            8.60353962

        E is a rank zero curve; the lowest zero has imaginary part ~6.36. The
        zero sum with tau=0 indicates that there are no zeros at the central
        point (otherwise the returned value would be at least 1).

        ::

            sage: Z = LFunctionZeroSum(E)
            sage: Z._zerosum_cauchy(Delta=1,tau=0) # tol 1.0e-13
            0.9701073984459051

        The zero sum with tau=6.36 indicates there might be a zero in the
        vicinity.

        ::

            sage: Z._zerosum_cauchy(Delta=1,tau=6.36261389) # tol 1.0e-13
            2.180904626331156

        However, there are no zeros with imaginary part close to 1.5.

        ::

            sage: Z._zerosum_cauchy(Delta=1,tau=1.5) # tol 1.0e-13
            0.9827072037553375

        Because of the weak convergence of the Dirichlet series close to the
        critical line, the bound will in general get *worse* for larger Delta.
        This can be mitigated somewhat by increasing the number of terms.

        ::

            sage: Z._zerosum_cauchy(Delta=1.5) # tol 1.0e-13
            12.93835258975716
            sage: Z._zerosum_cauchy(Delta=1.5,num_terms=100000) # tol 1.0e-13
            10.395183960836599

        An error will be thrown if a Delta value >= 2 is passed.

        ::

            sage: Z._zerosum_cauchy(Delta=2)
            Traceback (most recent call last):
            ...
            ValueError: bound not provably computable for Delta >= 2
        """
        if Delta >= 2:
            raise ValueError("bound not provably computable for Delta >= 2")
        Del = RDF(Delta)
        if num_terms is None:
            num_terms = int(exp(2 * self._pi * Del))

        if tau == 0:
            one = RDF.one()
            s = one / Del + one
            u, err = self.completed_logarithmic_derivative(s, num_terms)
        else:
            one = CDF.one()
            s = CDF(one / Del + one, tau)
            u, err = self.completed_logarithmic_derivative(s, num_terms)
            u = u.real()

        return (u + err) / Del


cdef class LFunctionZeroSum_EllipticCurve(LFunctionZeroSum_abstract):
    r"""
    Subclass for computing certain sums over zeros of an elliptic curve `L`-function
    without having to determine the zeros themselves.
    """
    cdef _E     # The Elliptic curve attached to self
    cdef _e     # PARI ellcurve object used to compute a_p values

    def __init__(self, E, N=None, ncpus=1):
        r"""
        Initialize ``self``.

        INPUT:

        - ``E`` -- an elliptic curve defined over the rational numbers

        - ``N`` -- (default: ``None``) if not ``None``, a positive integer equal to
          the conductor of E. This is passable so that rank estimation
          can be done for curves whose (large) conductor has been precomputed.

        - ``ncpus`` -- (default: 1) the number of CPUs to use for computations;
          if set to ``None``, the max available amount will be used

        EXAMPLES::

            sage: from sage.lfunctions.zero_sums import LFunctionZeroSum_EllipticCurve
            sage: E = EllipticCurve([1,0,0,3,-4])
            sage: Z = LFunctionZeroSum_EllipticCurve(E); Z
            Zero sum estimator for L-function attached to Elliptic Curve defined by y^2 + x*y = x^3 + 3*x - 4 over Rational Field
            sage: E = EllipticCurve("5077a")
            sage: Z = LFunctionZeroSum_EllipticCurve(E); Z
            Zero sum estimator for L-function attached to Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
        """
        self._k = ZZ(2)
        self._E = E
        if N is not None:
            self._level = N
        else:
            self._level = E.conductor()
        # PARI minicurve for computing a_p coefficients
        self._e = E.pari_mincurve()

        self._pi = RDF(pi)
        self._euler_gamma = RDF(euler_gamma)

        # These constants feature in most (all?) sums over the L-function's zeros
        self._C1 = RDF(self._level).log() / 2 - log(self._pi * 2)
        self._C0 = self._C1 - self._euler_gamma

        # Number of CPUs to use for computations
        if ncpus is None:
            self._ncpus = num_cpus()
            NCPUS = self._ncpus
        else:
            self._ncpus = ncpus
            NCPUS = self._ncpus

    def __repr__(self):
        r"""
        Representation of ``self``.

        EXAMPLES::

            sage: Z = LFunctionZeroSum(EllipticCurve("37a")); Z
            Zero sum estimator for L-function attached to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        s = "Zero sum estimator for L-function attached to "
        return s + str(self._E)

    def elliptic_curve(self):
        r"""
        Return the elliptic curve associated with ``self``.

        EXAMPLES::

            sage: E = EllipticCurve([23,100])
            sage: Z = LFunctionZeroSum(E)
            sage: Z.elliptic_curve()
            Elliptic Curve defined by y^2 = x^3 + 23*x + 100 over Rational Field
        """
        return self._E

    def lseries(self):
        r"""
        Return the `L`-series associated with ``self``.

        EXAMPLES::

            sage: E = EllipticCurve([23,100])
            sage: Z = LFunctionZeroSum(E)
            sage: Z.lseries()
            Complex L-series of the Elliptic Curve defined by y^2 = x^3 + 23*x + 100 over Rational Field
        """
        return self._E.lseries()

    def cn(self, n):
        r"""
        Return the `n`-th Dirichlet coefficient of the logarithmic
        derivative of the `L`-function attached to ``self``, shifted so that
        the critical line lies on the imaginary axis.

        The returned value is
        zero if `n` is not a perfect prime power;
        when `n=p^e` for `p` a prime of bad reduction it is `-a_p^e log(p)/p^e`,
        where `a_p` is `+1, -1` or `0` according to the reduction type of `p`;
        and when `n=p^e` for a prime `p` of good reduction, the value
        is `-(\alpha_p^e + \beta_p^e) \log(p)/p^e`, where `\alpha_p`
        and `\beta_p` are the two complex roots of the characteristic equation
        of Frobenius at `p` on `E`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT:

        A real number which (by Hasse's Theorem) is at
        most `2\frac{log(n)}{\sqrt{n}}` in magnitude.

        EXAMPLES::

            sage: E = EllipticCurve("11a")
            sage: Z = LFunctionZeroSum(E)
            sage: for n in range(12): print((n, Z.cn(n))) # tol 1.0e-13
            (0, 0.0)
            (1, 0.0)
            (2, 0.6931471805599453)
            (3, 0.3662040962227033)
            (4, 0.0)
            (5, -0.32188758248682003)
            (6, 0.0)
            (7, 0.555974328301518)
            (8, -0.34657359027997264)
            (9, 0.6103401603711721)
            (10, 0.0)
            (11, -0.21799047934530644)
        """
        n = ZZ(n)
        if n == 0 or n == 1:
            return RDF.zero()
        if not n.is_prime_power():
            return RDF.zero()

        n_float = RDF(n)
        if n.is_prime():
            logn = log(n_float)
            ap = self._E.ap(n)
            return -ap * logn / n_float
        else:
            p, e = n.perfect_power()
            ap = self._E.ap(p)
            logp = RDF(p).log()
            if p.divides(self._level):
                return - ap**e * logp / n_float
            a, b = ap, 2
            # Coefficients for higher powers obey recursion relation
            for n in range(2, e + 1):
                a, b = ap * a - p * b, a
            return -a * logp / n_float

    cdef double _sincsquared_summand_1(self,
                                       unsigned long n,
                                       double t,
                                       int ap,
                                       double p,
                                       double logp,
                                       double thetap,
                                       double sqrtp,
                                       double logq,
                                       double thetaq,
                                       double sqrtq,
                                       double z) noexcept:
        r"""
        Private cdef method to compute the logarithmic derivative
        summand for the sinc^2 sum at prime values for when
        n <= sqrt(bound), bound = exp(t)
        Called in self._zerosum_sincsquared_fast() method
        """
        ap = self._e.ellap(n)
        p = n
        sqrtp = c_sqrt(p)
        thetap = c_acos(ap / (2 * sqrtp))
        logp = c_log(p)

        sqrtq = 1
        thetaq = 0
        logq = logp

        z = 0
        while logq < t:
            sqrtq *= sqrtp
            thetaq += thetap
            z += 2 * c_cos(thetaq) * (t - logq) / sqrtq
            logq += logp
        return -z * logp

    cdef double _sincsquared_summand_2(self,
                                       unsigned long n,
                                       double t,
                                       int ap,
                                       double p,
                                       double logp) noexcept:
        r"""
        Private cdef method to compute the logarithmic derivative
        summand for the sinc^2 sum at prime values for when
        sqrt(bound) < n < bound, bound = exp(t)
        Called in self._zerosum_sincsquared_fast() method
        """
        ap = self._e.ellap(n)
        p = n
        logp = c_log(p)
        return -(t - logp) * (logp / p) * ap

    cpdef _zerosum_sincsquared_fast(self, Delta=1, bad_primes=None):
        r"""
        A faster cythonized implementation of self._zerosum_sincsquared().

        .. NOTE::

            This will only produce correct output if self._E is given by its
            global minimal model, i.e., if self._E.is_minimal() is ``True``.

        INPUT:

        - ``Delta`` -- positive real parameter defining the
          tightness of the zero sum, and thus the closeness of the returned
          estimate to the actual analytic rank of the form attached to ``self``

        - ``bad_primes`` -- (default: ``None``) if not ``None``, a list of primes dividing
          the level of the form attached to ``self``. This is passable so that this
          method can be run on curves whose conductor is large enough to warrant
          precomputing bad primes.

        OUTPUT:

        A positive real number that bounds the analytic rank of the modular form
        attached to ``self`` from above.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum_sincsquared`
            for the more general but slower version of this method.

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve("37a")
            sage: Z = LFunctionZeroSum(E)
            sage: print((E.rank(),Z._zerosum_sincsquared_fast(Delta=1))) # tol 1.0e-13
            (1, 1.0103840698356263)
            sage: E = EllipticCurve("121a")
            sage: Z = LFunctionZeroSum(E)
            sage: print((E.rank(),Z._zerosum_sincsquared_fast(Delta=1.5))) # tol 1.0e-13
            (0, 0.0104712060086507)
        """
        # If Delta>6.619, then we will most likely get overflow: some ap values
        # will be too large to fit into a c int
        if Delta > 6.619:
            raise ValueError("Delta value too large; will result in overflow")

        cdef double npi = self._pi
        cdef double twopi = npi * 2
        cdef double eg = self._euler_gamma

        cdef double t, u, w, y, z, expt, bound1, logp, logq
        cdef double thetap, thetaq, sqrtp, sqrtq, p, q
        cdef int ap, aq

        cdef unsigned long n
        cdef double N_double = self._level

        t = twopi * Delta
        expt = c_exp(t)

        u = t * (-eg + c_log(N_double) / 2 - c_log(twopi))
        w = npi**2 / 6 - (RDF.one() / expt).dilog()

        y = 0
        # Do bad primes first. Add correct contributions and subtract
        # incorrect contribution, since we'll add them back later on.
        if bad_primes is None:
            bad_primes = self._level.prime_divisors()
        bad_primes = [prime for prime in bad_primes if prime < expt]
        for prime in bad_primes:
            n = prime
            ap = self._e.ellap(n)
            p = n
            sqrtp = c_sqrt(p)
            thetap = c_acos(ap / (2 * sqrtp))
            logp = c_log(p)

            q = 1
            sqrtq = 1
            aq = 1
            thetaq = 0
            logq = logp

            z = 0
            while logq < t:
                q *= p
                sqrtq *= sqrtp
                aq *= ap
                thetaq += thetap
                # Actual value of this term
                z += (aq / q) * (t - logq)
                # Incorrect value of this term to be removed below
                z -= 2 * c_cos(thetaq) * (t - logq) / sqrtq
                logq += logp
            y -= z * logp

        # Good prime case. Bad primes are treated as good primes, but their
        # contribution here is cancelled out above; this way we don't
        # have to check if each prime divides the level or not.

        # Must deal with n=2,3,5 separately
        for m in [2, 3, 5]:
            n = m
            if n < expt:
                y += self._sincsquared_summand_1(n, t, ap, p, logp, thetap,
                                                 sqrtp, logq, thetaq, sqrtq, z)
        # Now iterate only over those n that are 1 or 5 mod 6
        n = 11
        # First: those n that are <= sqrt(bound)
        bound1 = c_exp(t / 2)
        while n <= bound1:
            if n_is_prime(n - 4):
                y += self._sincsquared_summand_1(n - 4, t, ap, p, logp, thetap,
                                                 sqrtp, logq, thetaq, sqrtq, z)
            if n_is_prime(n):
                y += self._sincsquared_summand_1(n, t, ap, p, logp, thetap,
                                                 sqrtp, logq, thetaq, sqrtq, z)
            n += 6
        # Unlucky split case where n-4 <= sqrt(bound) but n isn't
        if n - 4 <= bound1 < n:
            if n_is_prime(n - 4):
                y += self._sincsquared_summand_1(n - 4, t, ap, p, logp, thetap,
                                                 sqrtp, logq, thetaq, sqrtq, z)
            if n <= expt and n_is_prime(n):
                y += self._sincsquared_summand_2(n, t, ap, p, logp)
            n += 6
        # Now sqrt(bound)< n < bound, so we don't need to consider higher
        # prime power logarithmic derivative coefficients
        while n <= expt:
            if n_is_prime(n - 4):
                y += self._sincsquared_summand_2(n - 4, t, ap, p, logp)
            if n_is_prime(n):
                y += self._sincsquared_summand_2(n, t, ap, p, logp)
            n += 6
        # Case where n-4 <= t but n isn't
        n = n - 4
        if n <= expt and n_is_prime(n):
            y += self._sincsquared_summand_2(n, t, ap, p, logp)

        return RDF(2 * (u + w + y) / (t**2))

    def _get_residue_data(self, n):
        r"""
        Method called by ``self._zerosum_sincsquared_parallel()`` to determine
        the optimal residue class breakdown when sieving for primes.

        Return a list of small primes, the product thereof, and a list of
        residues coprime to the product.

        INPUT:

        - ``n`` -- positive integer denoting the number of required chunks

        OUTPUT:

        A triple ``(small_primes, M, residue_chunks)`` such that

          - ``small_primes`` -- list of small primes

          - ``modulus`` -- the product of the small primes

          - ``residue_chunks`` -- list of lists consisting of all integers
             less than the modulus that are coprime to it, broken into `n`
             sublists of approximately equal size

        EXAMPLES::

            sage: E = EllipticCurve("37a"); Z = LFunctionZeroSum(E)
            sage: Z._get_residue_data(8)
            ([2, 3, 5, 7],
             210,
             [[1, 37, 71, 107, 143, 179],
              [11, 41, 73, 109, 149, 181],
              [13, 43, 79, 113, 151, 187],
              [17, 47, 83, 121, 157, 191],
              [19, 53, 89, 127, 163, 193],
              [23, 59, 97, 131, 167, 197],
              [29, 61, 101, 137, 169, 199],
              [31, 67, 103, 139, 173, 209]])
        """
        # If n <=48, primes are sieved for modulo 210
        if n <= 48:
            small_primes = [2, 3, 5, 7]
            modulus = 210
            residue_list = [1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                            53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,
                            107, 109, 113, 121, 127, 131, 137, 139, 143, 149,
                            151, 157, 163, 167, 169, 173, 179, 181, 187, 191,
                            193, 197, 199, 209]
        # General case for n > 480
        else:
            from sage.rings.finite_rings.integer_mod import mod

            modulus, p = 2, 2
            small_primes, residue_list = [2], [1]
            num_residues = 1
            # Enlarge residue_list by repeatedly applying Chinese Remainder
            # Theorem
            while num_residues < n:
                p = next_prime(p)
                small_primes.append(p)
                g, h = (mod(p, modulus)**(-1)).lift(), (mod(modulus, p)**(-1)).lift()
                residue_list = [(a * p * g + b * modulus * h) % (modulus * p)
                                for a in residue_list
                                for b in range(1, p)]
                num_residues = num_residues * (p - 1)
                modulus *= p
            residue_list.sort()

        # Break residue_list into n chunks
        residue_chunks = [[residue_list[i] for i in range(len(residue_list))
                           if i % n == k] for k in range(n)]

        return small_primes, modulus, residue_chunks

    @parallel(ncpus=NCPUS)
    def _sum_over_residues(self, residue_sum_data):
        r"""
        Return the p-power sum over residues in a residue chunk
        """
        modulus, residues = residue_sum_data[0], residue_sum_data[1]

        cdef double y = 0
        cdef unsigned long n = residues[0]
        cdef double t = residue_sum_data[2]
        cdef double expt = residue_sum_data[3]
        cdef double bound1 = residue_sum_data[4]

        cdef double z = 0
        cdef double p = 0
        cdef double sqrtp = 0
        cdef double sqrtq = 0
        cdef double logp = 0
        cdef double logq = 0
        cdef double thetap = 0
        cdef double thetaq = 0
        cdef int ap = 0

        # Generate a list of increments so that n iterates over integers with
        # residues in the residue list
        len_increment_list = len(residues)
        increments = [residues[i + 1] - residues[i]
                      for i in range(len_increment_list - 1)]
        increments.append(modulus + residues[0] - residues[-1])

        i = 0
        # up to bound1=sqrt(expt), higher powers of p must be summed over too
        while n < bound1:
            if n_is_prime(n):
                y += self._sincsquared_summand_1(n, t, ap, p, logp,
                                                 thetap, sqrtp, logq,
                                                 thetaq, sqrtq, z)
            n += increments[i]
            # cycle over increments
            i += 1
            if i >= len_increment_list:
                i = 0

        # when bound1 <= n < expt, we don't need to consider higher powers of p
        while n < expt:
            if n_is_prime(n):
                y += self._sincsquared_summand_2(n, t, ap, p, logp)
            n += increments[i]
            # cycle over increments
            i += 1
            if i >= len_increment_list:
                i = 0

        return y

    def _zerosum_sincsquared_parallel(self,
                                      Delta=1,
                                      bad_primes=None,
                                      ncpus=None):
        r"""
        Parallelized implementation of self._zerosum_sincsquared_fast().

        Faster than self._zerosum_sincsquared_fast() when Delta >= ~1.75.

        .. NOTE::

            This will only produce correct output if self._E is given by its
            global minimal model, i.e. if self._E.is_minimal() is ``True``.

        INPUT:

        - ``Delta`` -- positive real parameter defining the
          tightness of the zero sum, and thus the closeness of the returned
          estimate to the actual analytic rank of the form attached to ``self``.

        - ``bad_primes`` -- (default: ``None``) if not ``None``, a list of primes dividing
          the level of the form attached to ``self``. This is passable so that this
          method can be run on curves whose conductor is large enough to warrant
          precomputing bad primes.

        - ``ncpus`` -- (default: ``None``) if not ``None``, a positive integer
          defining the number of CPUs to be used for the computation. If left as
          ``None``, the maximum available number of CPUs will be used.

        OUTPUT:

        A positive real number that bounds the analytic rank of the modular form
        attached to ``self`` from above.

        .. SEEALSO::

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum_sincsquared`
            for the more general but slower version of this method.

            :meth:`~sage.lfunctions.zero_sums.LFunctionZeroSum_abstract.zerosum`
            for the public method that calls this private method.

        EXAMPLES::

            sage: E = EllipticCurve("37a"); print(E.rank())
            1
            sage: Z = LFunctionZeroSum(E)
            sage: print(Z._zerosum_sincsquared_parallel(Delta=1)) # tol 1.0e-11
            1.0103840698356263
            sage: E = EllipticCurve("121a"); print(E.rank())
            0
            sage: Z = LFunctionZeroSum(E)
            sage: print(Z._zerosum_sincsquared_parallel(Delta=1.5,ncpus=2)) # tol 1.0e-11
            0.01047120600865063
        """
        # If Delta>6.619, then we will most likely get overflow: some ap values
        # will be too large to fit into a c int
        if Delta > 6.619:
            raise ValueError("Delta value too large; will result in overflow")

        cdef double npi = self._pi
        cdef double twopi = npi * 2
        cdef double eg = self._euler_gamma
        cdef double N_double = self._level

        cdef double t, u, w, y, z, expt, bound1, logp, logq
        cdef double thetap, thetaq, sqrtp, sqrtq, p, q
        cdef int ap, aq
        cdef unsigned long n

        # Compute bounds and smooth part of sum
        t = twopi * Delta
        expt = c_exp(t)
        bound1 = c_exp(t / 2)
        u = t * (-eg + c_log(N_double) / 2 - c_log(twopi))
        w = npi**2 / 6 - (RDF.one() / expt).dilog()

        # Oscillating part of sum
        y = 0
        # Do bad primes first. Add correct contributions and subtract
        # incorrect contribution; the latter are added back later on.
        if bad_primes is None:
            bad_primes = self._level.prime_divisors()
        bad_primes = [prime for prime in bad_primes if prime < expt]
        for prime in bad_primes:
            n = prime
            ap = self._e.ellap(n)
            p = n
            sqrtp = c_sqrt(p)
            thetap = c_acos(ap / (2 * sqrtp))
            logp = c_log(p)

            q = 1
            sqrtq = 1
            aq = 1
            thetaq = 0
            logq = logp

            z = 0
            while logq < t:
                q *= p
                sqrtq *= sqrtp
                aq *= ap
                thetaq += thetap
                # Actual value of this term
                z += (aq / q) * (t - logq)
                # Incorrect value of this term to be removed below
                z -= 2 * c_cos(thetaq) * (t - logq) / sqrtq
                logq += logp
            y -= z * logp

        # Good prime case. Bad primes are treated as good primes, but their
        # contribution here is cancelled out above; this way we don't
        # have to check if each prime divides the level or not.
        if ncpus is not None:
            self._ncpus = ncpus
        else:
            ncpus = self._ncpus
        small_primes, modulus, residue_chunks = self._get_residue_data(ncpus)

        # We must deal with small primes separately
        for m in small_primes:
            n = m
            if n < expt:
                y += self._sincsquared_summand_1(n, t, ap, p, logp,
                                                 thetap, sqrtp, logq,
                                                 thetaq, sqrtq, z)

        # Now the rest of the sum via _sum_over_residues(), which is parallelized
        residue_sum_data = []
        for residues in residue_chunks:
            residue_sum_data.append([modulus, residues, t, expt, bound1])
        # residue_data = [[modulus]+residue_chunks[i] for i in range(len(residue_chunks))]
        # for summand in _sum_over_residues(residue_data):
        for summand in self._sum_over_residues(residue_sum_data):
            y += summand[1]

        return RDF(2 * (u + w + y) / (t**2))

    def analytic_rank_upper_bound(self,
                                  max_Delta=None,
                                  adaptive=True,
                                  root_number='compute',
                                  bad_primes=None,
                                  ncpus=None):
        r"""
        Return an upper bound for the analytic rank of the `L`-function
        `L_E(s)` attached to ``self``, conditional on the Generalized Riemann
        Hypothesis, via computing the zero sum `\sum_{\gamma} f(\Delta\gamma)`,
        where `\gamma` ranges over the imaginary parts of the zeros of `L(E,s)`
        along the critical strip, `f(x) = \left(\frac{\sin(\pi x)}{\pi x}\right)^2`,
        and `\Delta` is the tightness parameter whose maximum value is
        specified by max_Delta.

        This computation can be run on curves with
        very large conductor (so long as the conductor is known or quickly
        computable) when Delta is not too large (see below).

        Uses Bober's rank bounding method as described in [Bob2013]_.

        INPUT:

        - ``max_Delta`` -- (default: ``None``) if not ``None``, a positive real value
          specifying the maximum Delta value used in the zero sum; larger
          values of Delta yield better bounds - but runtime is exponential in
          Delta. If left as ``None``, Delta is set
          to `\min\left\{\frac{1}{\pi}\left(\log(N+1000)/2-\log(2\pi)-\eta\right), 2.5\right\}`,
          where `N` is the conductor of the curve attached to ``self``, and `\eta`
          is the Euler-Mascheroni constant `= 0.5772...`; the crossover
          point is at conductor ~8.3*10^8. For the former value, empirical
          results show that for about 99.7% of all curves the returned value
          is the actual analytic rank.

        - ``adaptive`` -- boolean (default: ``True``)

          - If ``True``, the computation is first run with small and then
            successively larger Delta values up to max_Delta. If at any
            point the computed bound is 0 (or 1 when root_number is -1
            or ``True``), the computation halts and that value is returned;
            otherwise the minimum of the computed bounds is returned.
          - If ``False``, the computation is run a single time with
            Delta=max_Delta, and the resulting bound returned.

        - ``root_number`` -- (default: ``'compute'``) string or integer

          - ``'compute'`` -- the root number of ``self`` is computed and used to
            (possibly) lower the analytic rank estimate by 1
          - ``'ignore'`` -- the above step is omitted
          - ``1`` -- this value is assumed to be the root number of
            ``self``. This is passable so that rank estimation can be done for
            curves whose root number has been precomputed.
          - ``-1`` -- this value is assumed to be the root number of
            ``self``. This is passable so that rank estimation can be done for
            curves whose root number has been precomputed.

        - ``bad_primes`` -- (default: ``None``) if not ``None``, a list of the primes
          of bad reduction for the curve attached to ``self``. This is passable
          so that rank estimation can be done for curves of large conductor
          whose bad primes have been precomputed.

        - ``ncpus`` -- (default: ``None``) if not ``None``, a positive integer
          defining the maximum number of CPUs to be used for the computation.
          If left as ``None``, the maximum available number of CPUs will be used.
          Note: Multiple processors will only be used for Delta values >= 1.75.

        .. NOTE::

            Output will be incorrect if the incorrect root number is specified.

        .. WARNING::

            Zero sum computation time is exponential in the tightness parameter
            `\Delta`, roughly doubling for every increase of 0.1 thereof.
            Using `\Delta=1` (and adaptive=False) will yield a runtime of a few
            milliseconds; `\Delta=2` takes a few seconds, and `\Delta=3` may
            take upwards of an hour. Increase beyond this at your own risk!

        OUTPUT:

        A nonnegative integer greater than or equal to the analytic rank of
        ``self``. If the returned value is 0 or 1 (the latter if parity is not
        ``False``), then this is the true analytic rank of ``self``.

        .. NOTE::

            If you use set_verbose(1), extra information about the computation
            will be printed.

        .. SEEALSO::

            :func:`LFunctionZeroSum`
            :meth:`EllipticCurve.root_number`
            :func:`~sage.misc.verbose.set_verbose`

        EXAMPLES:

        For most elliptic curves with small conductor the central zero(s)
        of `L_E(s)` are fairly isolated, so small values of `\Delta`
        will yield tight rank estimates.

        ::

            sage: E = EllipticCurve("11a")
            sage: E.rank()
            0
            sage: Z = LFunctionZeroSum(E)
            sage: Z.analytic_rank_upper_bound(max_Delta=1,ncpus=1)
            0

            sage: E = EllipticCurve([-39,123])
            sage: E.rank()
            1
            sage: Z = LFunctionZeroSum(E)
            sage: Z.analytic_rank_upper_bound(max_Delta=1)
            1

        This is especially true for elliptic curves with large rank.

        ::

            sage: for r in range(9):
            ....:     E = elliptic_curves.rank(r)[0]
            ....:     print((r, E.analytic_rank_upper_bound(max_Delta=1,
            ....:     adaptive=False,root_number='ignore')))
            (0, 0)
            (1, 1)
            (2, 2)
            (3, 3)
            (4, 4)
            (5, 5)
            (6, 6)
            (7, 7)
            (8, 8)

        However, some curves have `L`-functions with low-lying zeroes, and for these
        larger values of `\Delta` must be used to get tight estimates.

        ::

            sage: E = EllipticCurve("974b1")
            sage: r = E.rank(); r
            0
            sage: Z = LFunctionZeroSum(E)
            sage: Z.analytic_rank_upper_bound(max_Delta=1,root_number='ignore')
            1
            sage: Z.analytic_rank_upper_bound(max_Delta=1.3,root_number='ignore')
            0

        Knowing the root number of E allows us to use smaller Delta values
        to get tight bounds, thus speeding up runtime considerably.

        ::

            sage: Z.analytic_rank_upper_bound(max_Delta=0.6,root_number='compute')
            0

        The are a small number of curves which have pathologically low-lying
        zeroes. For these curves, this method will produce a bound that is
        strictly larger than the analytic rank, unless very large values of
        Delta are used. The following curve ("256944c1" in the Cremona tables)
        is a rank 0 curve with a zero at 0.0256...; the smallest Delta value
        for which the zero sum is strictly less than 2 is ~2.815.

        ::

            sage: E = EllipticCurve([0, -1, 0, -7460362000712, -7842981500851012704])
            sage: N,r = E.conductor(),E.analytic_rank(); N, r
            (256944, 0)
            sage: E.analytic_rank_upper_bound(max_Delta=1,adaptive=False)
            2
            sage: E.analytic_rank_upper_bound(max_Delta=2,adaptive=False)
            2

        This method is can be called on curves with large conductor.

        ::

            sage: E = EllipticCurve([-2934,19238])
            sage: Z = LFunctionZeroSum(E)
            sage: Z.analytic_rank_upper_bound()
            1

        And it can bound rank on curves with *very* large conductor, so long as
        you know beforehand/can easily compute the conductor and primes of bad
        reduction less than `e^{2\pi\Delta}`. The example below is of the rank
        28 curve discovered by Elkies that is the elliptic curve of (currently)
        largest known rank.

        ::

            sage: a4 = -20067762415575526585033208209338542750930230312178956502
            sage: a6 = 34481611795030556467032985690390720374855944359319180361266008296291939448732243429
            sage: E = EllipticCurve([1,-1,1,a4,a6])
            sage: bad_primes = [2,3,5,7,11,13,17,19,48463]
            sage: N = 3455601108357547341532253864901605231198511505793733138900595189472144724781456635380154149870961231592352897621963802238155192936274322687070
            sage: Z = LFunctionZeroSum(E,N)
            sage: Z.analytic_rank_upper_bound(max_Delta=2.37,adaptive=False, # long time
            ....: root_number=1,bad_primes=bad_primes,ncpus=2)
            32
        """
        # Helper function: compute zero sum and apply parity if not False
        def run_computation(Delta):
            verbose("Computing zero sum with Delta = %s" % Delta)
            # Empirically, the non-parallelized zero sum method runs faster
            # for Delta <= 1.75, regardless of the number of available CPUs.
            if Delta <= 1.75:
                bound = self._zerosum_sincsquared_fast(Delta=Delta,
                                                       bad_primes=bad_primes)
            else:
                bound = self._zerosum_sincsquared_parallel(Delta=Delta,
                                                           bad_primes=bad_primes,
                                                           ncpus=ncpus)
            verbose("Sum value is %s" % bound)
            bound = bound.floor()
            # parity is set to -1 when we're not taking root number into
            # account
            if parity == -1:
                verbose("Without invoking parity, rank bound is %s" % bound)
                return bound
            # parity is 0 if E has even analytic rank, and 1 if odd
            # analytic rank. The returned value must have the same parity
            # as the parity parameter.
            if bound % 2 != parity:
                bound -= 1
            verbose("Invoking parity, rank bound is %s" % bound)
            return bound

        # Get/compute parity
        if root_number == 1 or root_number == -1:
            parity = (1 - root_number) // 2
            verbose("Parity set to %s." % parity)
        elif root_number == "compute":
            verbose("Computing curve parity...")
            parity = (1 - self._e.ellrootno()) // 2
            verbose("Curve has parity %s." % parity)
        elif root_number == "ignore":
            verbose("Curve parity ignored.")
            parity = -1
        else:
            raise ValueError("root_number parameter not recognized")
        if parity == 1:
            halt_bound = 1
            verbose("Computation will halt if at any point bound is <= 1.")
        else:
            halt_bound = 0
            verbose("Computation will halt if at any point bound is 0.")

        # Compute max_Delta if necessary
        if max_Delta is None:
            verbose("Computing maximum Delta value")
            pi, eg = self._pi, self._euler_gamma
            # 1000 is arbitrary - increases Delta for small N
            max_Delta = (RDF(self._level + 1000).log() / 2 - log(2 * pi) - eg) / pi
            if max_Delta > 2.5:
                max_Delta = 2.5
                verbose("Computed max Delta value too big; setting to 2.5")
            else:
                verbose("Maximum Delta value to be used set at %s" % max_Delta)
        else:
            verbose("Maximum Delta value to be used set at %s" % max_Delta)

        # When max_Delta <= 1 it's not worth running the computation
        # multiple times, as it's so quick anyway
        if not adaptive or max_Delta <= 1:
            return run_computation(max_Delta)
        else:
            bound_list = []
            # Find starting value. This loop won't ever take long,
            # since max_Delta is never > 7.
            Delta = max_Delta
            while Delta > 1:
                Delta -= 0.2
            # Now go back up the sequence of Deltas
            while Delta <= max_Delta:
                bound = run_computation(Delta)
                if bound <= halt_bound:
                    verbose("computed bound <= halt_bound, so halting")
                    return bound
                else:
                    bound_list.append(bound)
                # Incrementing Delta by 0.2 each step means runtime
                # will increase by a factor of ~3.7
                Delta += 0.2

            # Since the zero sum is not strictly decreasing in Delta,
            # the last value is not necessarily the smallest
            smallest_bound = min(bound_list)
            verbose("Smallest bound computed is %s" % smallest_bound)
            return smallest_bound


def LFunctionZeroSum(X, *args, **kwds):
    r"""
    Constructor for the LFunctionZeroSum class.

    INPUT:

    - ``X`` -- a motivic object; currently only implemented for X = an elliptic curve
      over the rational numbers

    OUTPUT: an LFunctionZeroSum object

    EXAMPLES::

        sage: E = EllipticCurve("389a")
        sage: Z = LFunctionZeroSum(E); Z
        Zero sum estimator for L-function attached to Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field

    TESTS::

        sage: E = EllipticCurve([1.2,3.8])
        sage: LFunctionZeroSum(E)
        Traceback (most recent call last):
        ...
        NotImplementedError: currently only implemented for elliptic curves over QQ

        sage: f = Newforms(46)[0]
        sage: LFunctionZeroSum(f)
        Traceback (most recent call last):
        ...
        NotImplementedError: currently only implemented for elliptic curves over QQ
    """
    # Here to avoid import recursion
    from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field

    if isinstance(X, EllipticCurve_rational_field):
        return LFunctionZeroSum_EllipticCurve(X, *args, **kwds)

    raise NotImplementedError("currently only implemented for elliptic curves over QQ")
