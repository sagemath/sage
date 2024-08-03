"""
`L`-series for elliptic curves

AUTHORS:

- Simon Spicer (2014-08-15): Added LFunctionZeroSum class interface method

- Jeroen Demeyer (2013-10-17): Compute L series with arbitrary precision
  instead of floats.

- William Stein et al. (2005 and later)
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein
#       Copyright (C) 2013 Jeroen Demeyer
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from math import ceil, log, sqrt

from sage.misc.cachefunc import cached_method
from sage.misc.verbose import verbose
from sage.rings.rational_field import RationalField
from sage.rings.real_mpfr import RealField
from sage.structure.sage_object import SageObject


class Lseries_ell(SageObject):
    """
    An elliptic curve `L`-series.
    """
    def __init__(self, E):
        r"""
        Create an elliptic curve `L`-series.

        EXAMPLES::

            sage: EllipticCurve([1..5]).lseries()
            Complex L-series of the Elliptic Curve
             defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field
        """
        self.__E = E

    def elliptic_curve(self):
        r"""
        Return the elliptic curve that this `L`-series is attached to.

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: L = E.lseries()
            sage: L.elliptic_curve ()
            Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
        """
        return self.__E

    def taylor_series(self, a=1, prec=53, series_prec=6, var='z'):
        r"""
        Return the Taylor series of this `L`-series about `a` to
        the given precision (in bits) and the number of terms.

        The output is a series in ``var``, where you should view ``var`` as
        equal to `s-a`.  Thus this function returns the formal power
        series whose coefficients are `L^{(n)}(a)/n!`.

        INPUT:

        - ``a`` -- complex number
        - ``prec`` -- integer; precision in bits (default: 53)
        - ``series_prec`` -- integer (default: 6)
        - ``var`` -- variable (default: ``'z'``)

        EXAMPLES::

            sage: E = EllipticCurve('389a')
            sage: L = E.lseries()
            sage: L.taylor_series(series_prec=3)   # abs tol 1e-14
            -1.27685190980159e-23 + (7.23588070754027e-24)*z + 0.759316500288427*z^2 + O(z^3)  # 32-bit
            1.34667664606157e-19 + (-7.63157535163667e-20)*z + 0.759316500288427*z^2 + O(z^3)  # 64-bit
        """
        D = self.dokchitser(prec)
        return D.taylor_series(a, series_prec, var)

    def _repr_(self):
        r"""
        Return string representation of this `L`-series.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: L = E.lseries()
            sage: L._repr_()
            'Complex L-series of the Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field'
        """
        return "Complex L-series of the %s" % self.__E

    def dokchitser(self, prec=53,
                   max_imaginary_part=0,
                   max_asymp_coeffs=40,
                   algorithm=None):
        r"""
        Return an interface for computing with the `L`-series
        of this elliptic curve.

        This provides a way to compute Taylor expansions and higher
        derivatives of `L`-series.

        INPUT:

        - ``prec`` -- integer (default: 53); bits precision

        - ``max_imaginary_part`` -- real number (default: 0)

        - ``max_asymp_coeffs`` -- integer (default: 40)

        - ``algorithm`` -- string; ``'gp'`` (default), ``'pari'``, or ``'magma'``

        If algorithm is ``'gp'``, this returns an interface to Tim
        Dokchitser's program for computing with the `L`-functions.

        If algorithm is "pari", this returns instead an interface to Pari's
        own general implementation of `L`-functions.

        .. NOTE::

            If algorithm='magma', then the precision is in digits rather
            than bits and the object returned is a Magma `L`-series, which has
            different functionality from the Sage `L`-series.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: L = E.lseries().dokchitser()
            sage: L(2)
            0.381575408260711
            sage: L = E.lseries().dokchitser(algorithm='magma')         # optional - magma
            sage: L.Evaluate(2)                                         # optional - magma
            0.38157540826071121129371040958008663667709753398892116

        If the curve has too large a conductor, it is not possible to
        compute with the `L`-series using this command.  Instead a
        :exc:`RuntimeError` is raised::

            sage: e = EllipticCurve([1,1,0,-63900,-1964465932632])
            sage: L = e.lseries().dokchitser(15, algorithm='gp')
            Traceback (most recent call last):
            ...
            RuntimeError: unable to create L-series, due to precision or other limits in PARI

        Using the "pari" algorithm::

            sage: E = EllipticCurve('37a')
            sage: L = E.lseries().dokchitser(algorithm='pari')
            sage: L(2)
            0.381575408260711
        """
        if algorithm is None:
            algorithm = 'pari'

        if algorithm == 'magma':
            from sage.interfaces.magma import magma
            return magma(self.__E).LSeries(Precision=prec)

        if algorithm == 'pari':
            from sage.lfunctions.pari import LFunction, lfun_elliptic_curve
            L = LFunction(lfun_elliptic_curve(self.__E), prec=prec)
            L.rename('PARI L-function associated to %s' % self.__E)
            return L

        if algorithm == 'gp':
            from sage.lfunctions.all import Dokchitser
            key = (prec, max_imaginary_part, max_asymp_coeffs)
            try:
                return self.__dokchitser[key]
            except KeyError:
                pass
            except AttributeError:
                self.__dokchitser = {}
            L = Dokchitser(conductor=self.__E.conductor(),
                           gammaV=[0, 1],
                           weight=2,
                           eps=self.__E.root_number(),
                           poles=[],
                           prec=prec)
            s = 'e = ellinit(%s);' % list(self.__E.minimal_model().a_invariants())
            s += 'a(k) = ellak(e, k);'
            L.init_coeffs('a(k)', 1, pari_precode=s,
                          max_imaginary_part=max_imaginary_part,
                          max_asymp_coeffs=max_asymp_coeffs)
            L.rename('Dokchitser L-function associated to %s' % self.__E)
            self.__dokchitser[key] = L
            return L

        raise ValueError('algorithm must be "gp", "pari" or "magma"')

    def sympow(self, n, prec):
        r"""
        Return `L( Sym^{(n)}(E, \text{edge}))` to ``prec`` digits
        of precision.

        INPUT:

        - ``n`` -- integer

        - ``prec`` -- integer

        OUTPUT: string; real number to ``prec`` digits of precision as a string

        .. NOTE::

            Before using this function for the first time for
            a given ``n``, you may have to type ``sympow('-new_data <n>')``,
            where ``<n>`` is replaced by your value of ``n``.  This
            command takes a long time to run.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: a = E.lseries().sympow(2,16)   # not tested - requires precomputing "sympow('-new_data 2')"
            sage: a                              # not tested
            '2.492262044273650E+00'
            sage: RR(a)                          # not tested
            2.49226204427365
        """
        from sage.lfunctions.sympow import sympow
        return sympow.L(self.__E, n, prec)

    def sympow_derivs(self, n, prec, d):
        r"""
        Return `0`-th to `d`-th derivatives of `L( Sym^{(n)}(E,
        \text{edge}))` to ``prec`` digits of precision.

        INPUT:

        - ``n`` -- integer

        - ``prec`` -- integer

        - ``d`` -- integer

        OUTPUT: string; exactly as output by sympow

        .. NOTE::

            To use this function you may have to run a few commands
            like ``sympow('-new_data 1d2')``, each which takes a few
            minutes.  If this function fails it will indicate what
            commands have to be run.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: print(E.lseries().sympow_derivs(1,16,2))    # not tested -- requires precomputing "sympow('-new_data 2')"
            sympow 1.018 RELEASE  (c) Mark Watkins --- see README and COPYING for details
            Minimal model of curve  is [0,0,1,-1,0]
            At 37: Inertia Group is  C1 MULTIPLICATIVE REDUCTION
            Conductor is 37
            sp 1: Conductor at 37 is 1+0, root number is 1
            sp 1: Euler factor at 37 is 1+1*x
            1st sym power conductor is 37, global root number is -1
            NT 1d0: 35
            NT 1d1: 32
            NT 1d2: 28
            Maximal number of terms is 35
            Done with small primes 1049
            Computed:  1d0  1d1  1d2
            Checked out:  1d1
             1n0: 3.837774351482055E-01
             1w0: 3.777214305638848E-01
             1n1: 3.059997738340522E-01
             1w1: 3.059997738340524E-01
             1n2: 1.519054910249753E-01
             1w2: 1.545605024269432E-01
        """
        from sage.lfunctions.sympow import sympow
        return sympow.Lderivs(self.__E, n, prec, d)

    def zeros(self, n):
        r"""
        Return the imaginary parts of the first `n` nontrivial zeros
        on the critical line of the `L`-function in the upper half
        plane, as 32-bit reals.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.lseries().zeros(2)
            [0.000000000, 5.00317001]

            sage: a = E.lseries().zeros(20)             # long time
            sage: point([(1,x) for x in a])             # graph  (long time)
            Graphics object consisting of 1 graphics primitive

        AUTHORS: Uses Rubinstein's L-functions calculator.
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.zeros(n, L=self.__E)

    def zeros_in_interval(self, x, y, stepsize):
        r"""
        Return the imaginary parts of (most of) the nontrivial zeros
        on the critical line `\Re(s)=1` with positive imaginary part
        between ``x`` and ``y``, along with a technical quantity for each.

        INPUT:

        - ``x`` -- positive floating point number
        - ``y`` -- positive floating point number
        - ``stepsize`` -- positive floating point number

        OUTPUT:

        - list of pairs ``(zero, S(T))``.

        Rubinstein writes: The first column outputs the imaginary part
        of the zero, the second column a quantity related to ``S(T)`` (it
        increases roughly by 2 whenever a sign change, i.e. pair of
        zeros, is missed). Higher up the critical strip you should use
        a smaller stepsize so as not to miss zeros.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.lseries().zeros_in_interval(6, 10, 0.1)      # long time
            [(6.87039122, 0.248922780), (8.01433081, -0.140168533), (9.93309835, -0.129943029)]
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.zeros_in_interval(x, y, stepsize, L=self.__E)

    def values_along_line(self, s0, s1, number_samples):
        r"""
        Return values of `L(E, s)` at ``number_samples``
        equally-spaced sample points along the line from `s_0` to
        `s_1` in the complex plane.

        .. NOTE::

            The `L`-series is normalized so that the center of the
            critical strip is 1.

        INPUT:

        - ``s0``, ``s1`` -- complex numbers
        - ``number_samples`` -- integer

        OUTPUT:

        list -- list of pairs (`s`, `L(E,s)`), where the `s` are
                equally spaced sampled points on the line from
                ``s0`` to ``s1``.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.lseries().values_along_line(1, 0.5 + 20*I, 5)
            [(0.500000000, ...),
             (0.400000000 + 4.00000000*I, 3.31920245 - 2.60028054*I),
             (0.300000000 + 8.00000000*I, -0.886341185 - 0.422640337*I),
             (0.200000000 + 12.0000000*I, -3.50558936 - 0.108531690*I),
             (0.100000000 + 16.0000000*I, -3.87043288 - 1.88049411*I)]
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.values_along_line(s0-RationalField()('1/2'),
                                       s1-RationalField()('1/2'),
                                       number_samples, L=self.__E)

    def twist_values(self, s, dmin, dmax):
        r"""
        Return values of `L(E, s, \chi_d)` for each quadratic
        character `\chi_d` for `d_{\min} \leq d \leq d_{\max}`.

        .. NOTE::

            The `L`-series is normalized so that the center of the
            critical strip is 1.

        INPUT:

        - ``s`` -- complex numbers

        - ``dmin`` -- integer

        - ``dmax`` -- integer

        OUTPUT:

        - list of pairs `(d, L(E, s, \chi_d))`

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: vals = E.lseries().twist_values(1, -12, -4)
            sage: vals[0][0]
            -11
            sage: vals[0][1] # abs tol 1e-8
            1.47824342 + 0.0*I
            sage: vals[1][0]
            -8
            sage: vals[1][1] # abs tol 1e-8
            0.0 + 0.0*I
            sage: vals[2][0]
            -7
            sage: vals[2][1] # abs tol 1e-8
            1.85307619 + 0.0*I
            sage: vals[3][0]
            -4
            sage: vals[3][1] # abs tol 1e-8
            2.45138938 + 0.0*I
            sage: F = E.quadratic_twist(-8)
            sage: F.rank()
            1
            sage: F = E.quadratic_twist(-7)
            sage: F.rank()
            0
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.twist_values(s - RationalField()('1/2'), dmin, dmax, L=self.__E)

    def twist_zeros(self, n, dmin, dmax):
        r"""
        Return first `n` real parts of nontrivial zeros of
        `L(E,s,\chi_d)` for each quadratic character `\chi_d` with
        `d_{\min} \leq d \leq d_{\max}`.

        .. NOTE::

            The `L`-series is normalized so that the center of the
            critical strip is 1.

        INPUT:

        - ``n`` -- integer

        - ``dmin`` -- integer

        - ``dmax`` -- integer

        OUTPUT:

        - ``dict`` -- keys are the discriminants `d`, and
          values are list of corresponding zeros

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.lseries().twist_zeros(3, -4, -3)         # long time
            {-4: [1.60813783, 2.96144840, 3.89751747], -3: [2.06170900, 3.48216881, 4.45853219]}
        """
        from sage.lfunctions.lcalc import lcalc
        return lcalc.twist_zeros(n, dmin, dmax, L=self.__E)

    def at1(self, k=None, prec=None):
        r"""
        Compute `L(E,1)` using `k` terms of the series for `L(E,1)` as
        explained in Section 7.5.3 of Henri Cohen's book *A Course in
        Computational Algebraic Number Theory*.  If the argument `k`
        is not specified, then it defaults to `\sqrt{N}`, where `N` is
        the conductor.

        INPUT:

        - ``k`` -- number of terms of the series; if zero or ``None``,
          use `k = \sqrt{N}`, where `N` is the conductor

        - ``prec`` -- numerical precision in bits; if zero or ``None``,
          use a reasonable automatic default

        OUTPUT:

        A tuple of real numbers ``(L, err)`` where ``L`` is an
        approximation for `L(E,1)` and ``err`` is a bound on the error
        in the approximation.

        This function is disjoint from the PARI :pari:`elllseries`
        command, which is for a similar purpose.  To use that command
        (via the PARI C library), simply type
        ``E.pari_mincurve().elllseries(1)``.

        ALGORITHM:

        - Compute the root number eps.  If it is -1, return 0.

        - Compute the Fourier coefficients `a_n`, for `n` up to and
          including `k`.

        - Compute the sum

          .. MATH::

              2 \cdot \sum_{n=1}^{k} \frac{a_n}{n} \cdot \exp(-2*pi*n/\sqrt{N}),

          where `N` is the conductor of `E`.

        - Compute a bound on the tail end of the series, which is

          .. MATH::

                2 e^{-2 \pi (k+1) / \sqrt{N}} / (1 - e^{-2 \pi/\sqrt{N}}).

          For a proof see [Grigov-Jorza-Patrascu-Patrikis-Stein].

        EXAMPLES::

            sage: L, err = EllipticCurve('11a1').lseries().at1()
            sage: L, err
            (0.253804, 0.000181444)
            sage: parent(L)
            Real Field with 24 bits of precision
            sage: E = EllipticCurve('37b')
            sage: E.lseries().at1()
            (0.7257177, 0.000800697)
            sage: E.lseries().at1(100)
            (0.7256810619361527823362055410263965487367603361763, 1.52469e-45)
            sage: L,err = E.lseries().at1(100, prec=128)
            sage: L
            0.72568106193615278233620554102639654873
            sage: parent(L)
            Real Field with 128 bits of precision
            sage: err
            1.70693e-37
            sage: parent(err)
            Real Field with 24 bits of precision and rounding RNDU

        Rank 1 through 3 elliptic curves::

            sage: E = EllipticCurve('37a1')
            sage: E.lseries().at1()
            (0.0000000, 0.000000)
            sage: E = EllipticCurve('389a1')
            sage: E.lseries().at1()
            (-0.001769566, 0.00911776)
            sage: E = EllipticCurve('5077a1')
            sage: E.lseries().at1()
            (0.0000000, 0.000000)
        """
        sqrtN = sqrt(self.__E.conductor())
        if k:
            k = int(k)
        else:
            k = int(ceil(sqrtN))

        if prec:
            prec = int(prec)
        else:
            # Use the same precision as deriv_at1() below for
            # consistency
            prec = int(9.065*k/sqrtN + 1.443*log(k)) + 12
        R = RealField(prec)
        # Compute error term with bounded precision of 24 bits and
        # round towards +infinity
        Rerror = RealField(24, rnd='RNDU')

        if self.__E.root_number() == -1:
            return (R.zero(), Rerror.zero())

        an = self.__E.anlist(k)  # list of Sage Integers
        pi = R.pi()
        sqrtN = R(self.__E.conductor()).sqrt()

        z = (-2*pi/sqrtN).exp()
        zpow = z
        # Compute series sum and accumulate floating point errors
        L = R.zero()
        error = Rerror.zero()

        for n in range(1, k + 1):
            term = (zpow * an[n])/n
            zpow *= z
            L += term
            # We express relative error in units of epsilon, where
            # epsilon is a number divided by 2^precision.
            # Instead of multiplying the error by 2 after the loop
            # (to account for L *= 2), we already multiply it now.
            #
            # For multiplication and division, the relative error
            # in epsilons is bounded by (1+e)^n - 1, where n is the
            # number of operations (assuming exact inputs).
            # exp(x) additionally multiplies this error by abs(x) and
            # adds one epsilon. The inputs pi and sqrtN each contribute
            # another epsilon.
            # Assuming that 2*pi/sqrtN <= 2, the relative error for z is
            # 7 epsilon. This implies a relative error of (8n-1) epsilon
            # for zpow. We add 2 for the computation of term and 1/2 to
            # compensate for the approximation (1+e)^n = 1+ne.
            #
            # The error of the addition is at most half an ulp of the
            # result.
            #
            # Multiplying everything by two gives:
            error += term.epsilon(Rerror)*(16*n + 3) + L.ulp(Rerror)
        L *= 2

        # Add series error (we use (-2)/(z-1) instead of 2/(1-z)
        # because this causes 1/(1-z) to be rounded up)
        error += ((-2)*Rerror(zpow)) / Rerror(z - 1)
        return (L, error)

    def deriv_at1(self, k=None, prec=None):
        r"""
        Compute `L'(E,1)` using `k` terms of the series for `L'(E,1)`,
        under the assumption that `L(E,1) = 0`.

        The algorithm used is from Section 7.5.3 of Henri Cohen's book
        *A Course in Computational Algebraic Number Theory*.

        INPUT:

        - ``k`` -- number of terms of the series; if zero or ``None``,
          use `k = \sqrt{N}`, where `N` is the conductor

        - ``prec`` -- numerical precision in bits; if zero or ``None``,
          use a reasonable automatic default

        OUTPUT:

        A tuple of real numbers ``(L1, err)`` where ``L1`` is an
        approximation for `L'(E,1)` and ``err`` is a bound on the error
        in the approximation.

        .. WARNING::

            This function only makes sense if `L(E)` has positive order
            of vanishing at 1, or equivalently if `L(E,1) = 0`.

        ALGORITHM:

        - Compute the root number eps.  If it is 1, return 0.

        - Compute the Fourier coefficients `a_n`, for `n` up to and
          including `k`.

        - Compute the sum

          .. MATH::

                 2 \cdot \sum_{n=1}^{k} (a_n / n) \cdot E_1(2 \pi n/\sqrt{N}),

          where `N` is the conductor of `E`, and `E_1` is the
          exponential integral function.

        - Compute a bound on the tail end of the series, which is

          .. MATH::

                 2 e^{-2 \pi (k+1) / \sqrt{N}} / (1 - e^{-2 \pi/\sqrt{N}}).

          For a proof see [Grigorov-Jorza-Patrascu-Patrikis-Stein].  This
          is exactly the same as the bound for the approximation to
          `L(E,1)` produced by :meth:`at1`.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.lseries().deriv_at1()                                               # needs sage.symbolic
            (0.3059866, 0.000801045)
            sage: E.lseries().deriv_at1(100)                                            # needs sage.symbolic
            (0.3059997738340523018204836833216764744526377745903, 1.52493e-45)
            sage: E.lseries().deriv_at1(1000)                                           # needs sage.symbolic
            (0.305999773834052301820483683321676474452637774590771998..., 2.75031e-449)

        With less numerical precision, the error is bounded by numerical accuracy::

            sage: # needs sage.symbolic
            sage: L, err = E.lseries().deriv_at1(100, prec=64)
            sage: L, err
            (0.305999773834052302, 5.55318e-18)
            sage: parent(L)
            Real Field with 64 bits of precision
            sage: parent(err)
            Real Field with 24 bits of precision and rounding RNDU

        Rank 2 and rank 3 elliptic curves::

            sage: E = EllipticCurve('389a1')
            sage: E.lseries().deriv_at1()                                               # needs sage.symbolic
            (0.0000000, 0.000000)
            sage: E = EllipticCurve((1, 0, 1, -131, 558))  # curve 59450i1
            sage: E.lseries().deriv_at1()                                               # needs sage.symbolic
            (-0.00010911444, 0.142428)
            sage: E.lseries().deriv_at1(4000)                                           # needs sage.symbolic
            (6.990...e-50, 1.31318e-43)
        """
        sqrtN = sqrt(self.__E.conductor())
        if k:
            k = int(k)
        else:
            k = int(ceil(sqrtN))

        if prec:
            prec = int(prec)
        else:
            # Estimate number of bits for the computation, based on error
            # estimate below (the denominator of that error is close enough
            # to 1 that we can ignore it).
            # 9.065 = 2*Pi/log(2)
            # 1.443 = 1/log(2)
            # 12 is an arbitrary extra number of bits (it is chosen
            #    such that the precision is 24 bits when the conductor
            #    equals 11 and k is the default value 4)
            prec = int(9.065*k/sqrtN + 1.443*log(k)) + 12
        R = RealField(prec)
        # Compute error term with bounded precision of 24 bits and
        # round towards +infinity
        Rerror = RealField(24, rnd='RNDU')

        if self.__E.root_number() == 1:
            # Order of vanishing at 1 of L(E) is even and assumed to be
            # positive, so L'(E,1) = 0.
            return (R.zero(), Rerror.zero())

        from sage.functions.exp_integral import exponential_integral_1

        an = self.__E.anlist(k)  # list of Sage Integers
        pi = R.pi()
        sqrtN = R(self.__E.conductor()).sqrt()
        v = exponential_integral_1(2*pi/sqrtN, k)

        # Compute series sum and accumulate floating point errors
        L = R.zero()
        error = Rerror.zero()
        # Sum of |an[n]|/n
        sumann = Rerror.zero()

        for n in range(1, k + 1):
            term = (v[n-1] * an[n])/n
            L += term
            error += term.epsilon(Rerror)*5 + L.ulp(Rerror)
            sumann += Rerror(an[n].abs())/n
        L *= 2

        # Add error term for exponential_integral_1() errors.
        # Absolute error for 2*v[i] is 4*max(1, v[0])*2^-prec
        if v[0] > 1.0:
            sumann *= Rerror(v[0])
        error += (sumann >> (prec - 2))

        # Add series error (we use (-2)/(z-1) instead of 2/(1-z)
        # because this causes 1/(1-z) to be rounded up)
        z = (-2*pi/sqrtN).exp()
        zpow = ((-2*(k+1))*pi/sqrtN).exp()
        error += ((-2)*Rerror(zpow)) / Rerror(z - 1)
        return (L, error)

    def __call__(self, s):
        r"""
        Return the value of the `L`-series of the elliptic curve E at s, where s
        must be a real number.

        .. NOTE::

            If the conductor of the curve is large, say `>10^{12}`,
            then this function will take a very long time, since it
            uses an `O(\sqrt{N})` algorithm.

        EXAMPLES::

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: L = E.lseries()
            sage: L(1)
            0.000000000000000
            sage: L(1.1)
            0.285491007678148
            sage: L(1.1 + I)
            0.174851377216615 + 0.816965038124457*I
        """
        return self.dokchitser()(s)

    def L1_vanishes(self):
        r"""
        Return whether or not `L(E,1) = 0`. The result is provably
        correct if the Manin constant of the associated optimal
        quotient is <= 2.  This hypothesis on the Manin constant
        is true for all curves of conductor <= 40000 (by Cremona) and
        all semistable curves (i.e., squarefree conductor).

        ALGORITHM: see :meth:`L_ratio`.

        EXAMPLES::

            sage: E = EllipticCurve([0, -1, 1, -10, -20])   # 11A  = X_0(11)
            sage: E.lseries().L1_vanishes()
            False
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.lseries().L1_vanishes()
            False
            sage: E = EllipticCurve([0, 0, 1, -1, 0])       # 37A  (rank 1)
            sage: E.lseries().L1_vanishes()
            True
            sage: E = EllipticCurve([0, 1, 1, -2, 0])       # 389A (rank 2)
            sage: E.lseries().L1_vanishes()
            True
            sage: E = EllipticCurve([0, 0, 1, -38, 90])     # 361A (CM curve))
            sage: E.lseries().L1_vanishes()
            True
            sage: E = EllipticCurve([0,-1,1,-2,-1])         # 141C (13-isogeny)
            sage: E.lseries().L1_vanishes()
            False

        AUTHORS: William Stein, 2005-04-20.
        """
        return self.L_ratio() == 0

    @cached_method
    def L_ratio(self):
        r"""
        Return the ratio `L(E,1) / \Omega` as an exact rational number.

        The result is *provably* correct if the Manin
        constant of the associated optimal quotient is `\leq 2`.  This
        hypothesis on the Manin constant is true for all semistable
        curves (i.e., squarefree conductor), by a theorem of Mazur
        from his *Rational Isogenies of Prime Degree* paper.

        EXAMPLES::

            sage: E = EllipticCurve([0, -1, 1, -10, -20])   # 11A  = X_0(11)
            sage: E.lseries().L_ratio()
            1/5
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.lseries().L_ratio()
            1/25
            sage: E = EllipticCurve([0, 0, 1, -1, 0])       # 37A  (rank 1)
            sage: E.lseries().L_ratio()
            0
            sage: E = EllipticCurve([0, 1, 1, -2, 0])       # 389A (rank 2)
            sage: E.lseries().L_ratio()
            0
            sage: E = EllipticCurve([0, 0, 1, -38, 90])     # 361A (CM curve))
            sage: E.lseries().L_ratio()
            0
            sage: E = EllipticCurve([0,-1,1,-2,-1])         # 141C (13-isogeny)
            sage: E.lseries().L_ratio()
            1
            sage: E = EllipticCurve(RationalField(), [1, 0, 0, 1/24624, 1/886464])
            sage: E.lseries().L_ratio()
            2

        See :issue:`3651` and :issue:`15299`::

            sage: EllipticCurve([0,0,0,-193^2,0]).sha().an()
            4
            sage: EllipticCurve([1, 0, 1, -131, 558]).sha().an()  # long time
            1.00000000000000

        ALGORITHM: Compute the root number.  If it is -1 then `L(E,s)`
        vanishes to odd order at 1, hence vanishes.  If it is +1, use
        a result about modular symbols and Mazur's *Rational Isogenies*
        paper to determine a provably correct bound (assuming Manin
        constant is <= 2) so that we can determine whether `L(E,1) = 0`.

        AUTHORS: William Stein, 2005-04-20.
        """
        if not self.__E.is_minimal():
            return self.__E.minimal_model().lseries().L_ratio()

        QQ = RationalField()
        if self.__E.root_number() == -1:
            return QQ.zero()

        # Even root number.  Decide if L(E,1) = 0.  If E is a modular
        # *OPTIMAL* quotient of J_0(N) elliptic curve, we know that T *
        # L(E,1)/omega is an integer n, where T is the order of the
        # image of the rational torsion point (0)-(oo) in E(Q), and
        # omega is the least real Neron period.  (This is proved in my
        # Ph.D. thesis, but is probably well known.)  We can easily
        # compute omega to very high precision using AGM.  So to prove
        # that L(E,1) = 0 we compute T/omega * L(E,1) to sufficient
        # precision to determine it as an integer.  If eps is the
        # error in computation of L(E,1), then the error in computing
        # the product is (2T/Omega_E) * eps, and we need this to be
        # less than 0.5, i.e.,
        #          (2T/Omega_E) * eps < 0.5,
        # so
        #          eps < 0.5 * Omega_E / (2T) = Omega_E / (4*T).
        #
        # Since in general E need not be optimal, we have to choose
        # eps = Omega_E/(8*t*B), where t is the exponent of E(Q)_tor,
        # and is a multiple of the degree of an isogeny between E
        # and the optimal curve.
        #
        # NOTE: We *do* have to worry about the Manin constant, since
        # we are using the Neron model to compute omega, not the
        # newform.  My theorem replaces the omega above by omega/c,
        # where c is the Manin constant, and the bound must be
        # correspondingly smaller.  If the level is square free, then
        # the Manin constant is 1 or 2, so there's no problem (since
        # we took 8 instead of 4 in the denominator).  If the level
        # is divisible by a square, then the Manin constant could
        # be a divisible by an arbitrary power of that prime, except
        # that Edixhoven claims the primes that appear are <= 7.

        t = self.__E.torsion_subgroup().order()
        omega = self.__E.period_lattice().basis()[0]
        d = self.__E._multiple_of_degree_of_isogeny_to_optimal_curve()
        C = 8 * d * t
        eps = omega / C

        sqrtN = 2 * self.__E.conductor().isqrt()
        k = sqrtN + 10
        while True:
            L1, error_bound = self.at1(k)
            if error_bound < eps:
                n = (L1 * C / omega).round()
                quo = QQ((n, C))
                return quo / self.__E.real_components()
            k += sqrtN
            verbose("Increasing precision to %s terms." % k)

    def zero_sums(self, N=None):
        r"""
        Return an ``LFunctionZeroSum`` class object for efficient computation
        of sums over the zeros of ``self``.

        This can be used to bound analytic rank from above without
        having to compute with the `L`-series directly.

        INPUT:

        - ``N`` -- (default: ``None``) if not ``None``, the conductor of the
          elliptic curve attached to ``self``. This is passable so that zero
          sum computations can be done on curves for which the conductor
          has been precomputed.

        OUTPUT: a ``LFunctionZeroSum_EllipticCurve`` instance

        EXAMPLES::

            sage: E = EllipticCurve("5077a")
            sage: E.lseries().zero_sums()
            Zero sum estimator for L-function attached to
             Elliptic Curve defined by y^2 + y = x^3 - 7*x + 6 over Rational Field
        """
        from sage.lfunctions.zero_sums import LFunctionZeroSum
        return LFunctionZeroSum(self.__E, N=N)
