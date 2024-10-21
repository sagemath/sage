# sage_setup: distribution = sagemath-categories
"""
Number-theoretic functions
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
import math
import sys

from sage.misc.lazy_import import lazy_import
from sage.misc.misc import increase_recursion_limit
from sage.rings.integer_ring import ZZ
from sage.symbolic.function import GinacFunction, BuiltinFunction

lazy_import('sage.functions.gamma', 'psi')
lazy_import('sage.functions.other', 'factorial')

lazy_import('sage.combinat.combinat', 'bernoulli_polynomial')
lazy_import('sage.rings.cc', 'CC')
lazy_import('sage.rings.complex_mpfr', ['ComplexField', 'ComplexNumber'])
lazy_import('sage.rings.polynomial.polynomial_real_mpfr_dense', 'PolynomialRealDense')
lazy_import('sage.rings.real_double', 'RDF')
lazy_import('sage.rings.real_mpfr', ['RR', 'RealField', 'RealNumber'])

lazy_import('sage.libs.mpmath.sage_utils', 'call', as_='_mpmath_call')
lazy_import('sage.libs.mpmath.all', 'zeta', as_='_mpmath_zeta')


class Function_zeta(GinacFunction):
    def __init__(self):
        r"""
        Riemann zeta function at s with s a real or complex number.

        INPUT:

        - ``s`` -- real or complex number

        If s is a real number, the computation is done using the MPFR
        library. When the input is not real, the computation is done using
        the PARI C library.

        EXAMPLES::

            sage: RR = RealField(200)                                                   # needs sage.rings.real_mpfr
            sage: zeta(RR(2))                                                           # needs sage.rings.real_mpfr
            1.6449340668482264364724151666460251892189499012067984377356

            sage: # needs sage.symbolic
            sage: zeta(x)
            zeta(x)
            sage: zeta(2)
            1/6*pi^2
            sage: zeta(2.)
            1.64493406684823
            sage: zeta(I)
            zeta(I)
            sage: zeta(I).n()
            0.00330022368532410 - 0.418155449141322*I
            sage: zeta(sqrt(2))
            zeta(sqrt(2))
            sage: zeta(sqrt(2)).n()  # rel tol 1e-10
            3.02073767948603

        It is possible to use the ``hold`` argument to prevent
        automatic evaluation::

            sage: zeta(2, hold=True)                                                    # needs sage.symbolic
            zeta(2)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = zeta(2, hold=True); a.simplify()                                  # needs sage.symbolic
            1/6*pi^2

        The Laurent expansion of `\zeta(s)` at `s=1` is
        implemented by means of the
        :wikipedia:`Stieltjes constants <Stieltjes_constants>`::

            sage: s = SR('s')                                                           # needs sage.symbolic
            sage: zeta(s).series(s==1, 2)                                               # needs sage.symbolic
            1*(s - 1)^(-1) + euler_gamma + (-stieltjes(1))*(s - 1) + Order((s - 1)^2)

        Generally, the Stieltjes constants occur in the Laurent
        expansion of `\zeta`-type singularities::

            sage: zeta(2*s/(s+1)).series(s==1, 2)                                       # needs sage.symbolic
            2*(s - 1)^(-1) + (euler_gamma + 1) + (-1/2*stieltjes(1))*(s - 1) + Order((s - 1)^2)


        TESTS::

            sage: # needs sage.symbolic
            sage: latex(zeta(x))
            \zeta(x)
            sage: a = loads(dumps(zeta(x)))
            sage: a.operator() == zeta
            True
            sage: zeta(x)._sympy_()                                                     # needs sympy
            zeta(x)

            sage: zeta(1)                                                               # needs sage.symbolic
            Infinity
            sage: zeta(x).subs(x=1)                                                     # needs sage.symbolic
            Infinity

        Check that :issue:`19799` is resolved::

            sage: zeta(pi)                                                              # needs sage.symbolic
            zeta(pi)
            sage: zeta(pi).n()  # rel tol 1e-10                                         # needs sage.symbolic
            1.17624173838258

        Check that :issue:`20082` is fixed::

            sage: zeta(x).series(x==pi, 2)                                              # needs sage.symbolic
            (zeta(pi)) + (zetaderiv(1, pi))*(-pi + x) + Order((pi - x)^2)
            sage: (zeta(x) * 1/(1 - exp(-x))).residue(x==2*pi*I)                        # needs sage.symbolic
            zeta(2*I*pi)

        Check that :issue:`20102` is fixed::

            sage: (zeta(x)^2).series(x==1, 1)                                           # needs sage.symbolic
            1*(x - 1)^(-2) + (2*euler_gamma)*(x - 1)^(-1)
            + (euler_gamma^2 - 2*stieltjes(1)) + Order(x - 1)
            sage: (zeta(x)^4).residue(x==1)                                             # needs sage.symbolic
            4/3*euler_gamma*(3*euler_gamma^2 - 2*stieltjes(1))
            - 28/3*euler_gamma*stieltjes(1) + 2*stieltjes(2)

        Check that the right infinities are returned (:issue:`19439`)::

            sage: zeta(1.0)                                                             # needs sage.symbolic
            +infinity
            sage: zeta(SR(1.0))                                                         # needs sage.symbolic
            Infinity

        Fixed conversion::

            sage: zeta(3)._maple_init_()                                                # needs sage.symbolic
            'Zeta(3)'
            sage: zeta(3)._maple_().sage()      # optional - maple                      # needs sage.symbolic
            zeta(3)
        """
        GinacFunction.__init__(self, 'zeta',
                               conversions={'giac': 'Zeta',
                                            'maple': 'Zeta',
                                            'sympy': 'zeta',
                                            'mathematica': 'Zeta'})


zeta = Function_zeta()


class Function_stieltjes(GinacFunction):
    def __init__(self):
        r"""
        Stieltjes constant of index ``n``.

        ``stieltjes(0)`` is identical to the Euler-Mascheroni constant
        (:class:`sage.symbolic.constants.EulerGamma`). The Stieltjes
        constants are used in the series expansions of `\zeta(s)`.

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: _ = var('n')
            sage: stieltjes(n)
            stieltjes(n)
            sage: stieltjes(0)
            euler_gamma
            sage: stieltjes(2)
            stieltjes(2)
            sage: stieltjes(int(2))
            stieltjes(2)
            sage: stieltjes(2).n(100)
            -0.0096903631928723184845303860352
            sage: RR = RealField(200)                                                   # needs sage.rings.real_mpfr
            sage: stieltjes(RR(2))                                                      # needs sage.rings.real_mpfr
            -0.0096903631928723184845303860352125293590658061013407498807014

        It is possible to use the ``hold`` argument to prevent
        automatic evaluation::

            sage: stieltjes(0, hold=True)                                               # needs sage.symbolic
            stieltjes(0)

            sage: # needs sage.symbolic
            sage: latex(stieltjes(n))
            \gamma_{n}
            sage: a = loads(dumps(stieltjes(n)))
            sage: a.operator() == stieltjes
            True
            sage: stieltjes(x)._sympy_()                                                # needs sympy
            stieltjes(x)

            sage: stieltjes(x).subs(x==0)                                               # needs sage.symbolic
            euler_gamma
        """
        GinacFunction.__init__(self, "stieltjes", nargs=1,
                            conversions=dict(mathematica='StieltjesGamma',
                                sympy='stieltjes'),
                            latex_name=r'\gamma')


stieltjes = Function_stieltjes()


class Function_HurwitzZeta(BuiltinFunction):
    def __init__(self):
        r"""
        TESTS::

            sage: latex(hurwitz_zeta(x, 2))                                             # needs sage.symbolic
            \zeta\left(x, 2\right)
            sage: hurwitz_zeta(x, 2)._sympy_()                                          # needs sympy sage.symbolic
            zeta(x, 2)
        """
        BuiltinFunction.__init__(self, 'hurwitz_zeta', nargs=2,
                                 conversions=dict(mathematica='HurwitzZeta',
                                                  sympy='zeta'),
                                 latex_name=r'\zeta')

    def _eval_(self, s, x):
        r"""
        TESTS::

            sage: # needs sage.symbolic
            sage: hurwitz_zeta(x, 1)
            zeta(x)
            sage: hurwitz_zeta(4, 3)
            1/90*pi^4 - 17/16
            sage: hurwitz_zeta(-4, x)
            -1/5*x^5 + 1/2*x^4 - 1/3*x^3 + 1/30*x
            sage: hurwitz_zeta(0, x)
            -x + 1/2

            sage: hurwitz_zeta(3, 0.5)                                                  # needs mpmath
            8.41439832211716
        """
        if x == 1:
            return zeta(s)
        if s in ZZ and s > 1:
            return ((-1) ** s) * psi(s - 1, x) / factorial(s - 1)
        elif s in ZZ and s <= 0:
            return -bernoulli_polynomial(x, -s + 1) / (-s + 1)
        else:
            return

    def _evalf_(self, s, x, parent=None, algorithm=None):
        r"""
        TESTS::

            sage: hurwitz_zeta(11/10, 1/2).n()                                          # needs mpmath sage.symbolic
            12.1038134956837
            sage: hurwitz_zeta(11/10, 1/2).n(100)                                       # needs mpmath sage.symbolic
            12.103813495683755105709077413
            sage: hurwitz_zeta(11/10, 1 + 1j).n()                                       # needs mpmath sage.rings.real_mpfr
            9.85014164287853 - 1.06139499403981*I
        """
        return _mpmath_call(_mpmath_zeta, s, x, parent=parent)

    def _derivative_(self, s, x, diff_param):
        r"""
        TESTS::

            sage: y = var('y')                                                          # needs sage.symbolic
            sage: diff(hurwitz_zeta(x, y), y)                                           # needs sage.symbolic
            -x*hurwitz_zeta(x + 1, y)
        """
        if diff_param == 1:
            return -s * hurwitz_zeta(s + 1, x)
        else:
            raise NotImplementedError('derivative with respect to first '
                                      'argument')


hurwitz_zeta_func = Function_HurwitzZeta()


def hurwitz_zeta(s, x, **kwargs):
    r"""
    The Hurwitz zeta function `\zeta(s, x)`, where `s` and `x` are complex.

    The Hurwitz zeta function is one of the many zeta functions. It
    is defined as

    .. MATH::

             \zeta(s, x) = \sum_{k=0}^{\infty} (k + x)^{-s}.


    When `x = 1`, this coincides with Riemann's zeta function.
    The Dirichlet `L`-functions may be expressed as linear combinations
    of Hurwitz zeta functions.

    EXAMPLES:

    Symbolic evaluations::

        sage: # needs sage.symbolic
        sage: hurwitz_zeta(x, 1)
        zeta(x)
        sage: hurwitz_zeta(4, 3)
        1/90*pi^4 - 17/16
        sage: hurwitz_zeta(-4, x)
        -1/5*x^5 + 1/2*x^4 - 1/3*x^3 + 1/30*x
        sage: hurwitz_zeta(7, -1/2)
        127*zeta(7) - 128
        sage: hurwitz_zeta(-3, 1)
        1/120

    Numerical evaluations::

        sage: hurwitz_zeta(3, 1/2).n()                                                  # needs mpmath
        8.41439832211716
        sage: hurwitz_zeta(11/10, 1/2).n()                                              # needs sage.symbolic
        12.1038134956837
        sage: hurwitz_zeta(3, x).series(x, 60).subs(x=0.5).n()                          # needs sage.symbolic
        8.41439832211716
        sage: hurwitz_zeta(3, 0.5)                                                      # needs mpmath
        8.41439832211716

    REFERENCES:

    - :wikipedia:`Hurwitz_zeta_function`
    """
    return hurwitz_zeta_func(s, x, **kwargs)


class Function_zetaderiv(GinacFunction):
    def __init__(self):
        r"""
        Derivatives of the Riemann zeta function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: zetaderiv(1, x)
            zetaderiv(1, x)
            sage: zetaderiv(1, x).diff(x)
            zetaderiv(2, x)
            sage: var('n')
            n
            sage: zetaderiv(n, x)
            zetaderiv(n, x)
            sage: zetaderiv(1, 4).n()
            -0.0689112658961254
            sage: import mpmath; mpmath.diff(lambda x: mpmath.zeta(x), 4)               # needs mpmath
            mpf('-0.068911265896125382')

        TESTS::

            sage: latex(zetaderiv(2, x))                                                # needs sage.symbolic
            \zeta^\prime\left(2, x\right)
            sage: a = loads(dumps(zetaderiv(2, x)))                                     # needs sage.symbolic
            sage: a.operator() == zetaderiv                                             # needs sage.symbolic
            True

            sage: b = RBF(3/2, 1e-10)                                                   # needs sage.libs.flint
            sage: zetaderiv(1, b, hold=True)                                            # needs sage.libs.flint sage.symbolic
            zetaderiv(1, [1.500000000 +/- 1.01e-10])
            sage: zetaderiv(b, 1)                                                       # needs sage.libs.flint sage.symbolic
            zetaderiv([1.500000000 +/- 1.01e-10], 1)
        """
        GinacFunction.__init__(self, "zetaderiv", nargs=2,
                               conversions=dict(maple='Zeta'))

    def _evalf_(self, n, x, parent=None, algorithm=None):
        r"""
        TESTS::

            sage: zetaderiv(0, 3, hold=True).n() == zeta(3).n()                         # needs sage.symbolic
            True
            sage: zetaderiv(2, 3 + I).n()                                               # needs sage.symbolic
            0.0213814086193841 - 0.174938812330834*I
        """
        return _mpmath_call(_mpmath_zeta, x, 1, n, parent=parent)

    def _method_arguments(self, k, x, **args):
        r"""
        TESTS::

            sage: zetaderiv(1, RBF(3/2, 0.0001))                                        # needs sage.libs.flint
            [-3.93 +/- ...e-3]
        """
        return [x, k]


zetaderiv = Function_zetaderiv()


def zeta_symmetric(s):
    r"""
    Completed function `\xi(s)` that satisfies
    `\xi(s) = \xi(1-s)` and has zeros at the same points as the
    Riemann zeta function.

    INPUT:

    - ``s`` -- real or complex number

    If s is a real number the computation is done using the MPFR
    library. When the input is not real, the computation is done using
    the PARI C library.

    More precisely,

    .. MATH::

                xi(s) = \gamma(s/2 + 1) * (s-1) * \pi^{-s/2} * \zeta(s).

    EXAMPLES::

        sage: # needs sage.rings.real_mpfr
        sage: RR = RealField(200)
        sage: zeta_symmetric(RR(0.7))
        0.49758041465112690357779107525638385212657443284080589766062

        sage: # needs sage.libs.pari sage.rings.real_mpfr
        sage: zeta_symmetric(0.7)
        0.497580414651127
        sage: zeta_symmetric(1 - 0.7)
        0.497580414651127
        sage: C.<i> = ComplexField()
        sage: zeta_symmetric(0.5 + i*14.0)
        0.000201294444235258 + 1.49077798716757e-19*I
        sage: zeta_symmetric(0.5 + i*14.1)
        0.0000489893483255687 + 4.40457132572236e-20*I
        sage: zeta_symmetric(0.5 + i*14.2)
        -0.0000868931282620101 + 7.11507675693612e-20*I

    REFERENCE:

    - I copied the definition of xi from
      http://web.viu.ca/pughg/RiemannZeta/RiemannZetaLong.html
    """
    if not (isinstance(s, ComplexNumber) or isinstance(s, RealNumber)):
        s = ComplexField()(s)

    R = s.parent()
    if s == 1:  # deal with poles, hopefully
        return R(0.5)

    return (s/2 + 1).gamma() * (s-1) * (R.pi()**(-s/2)) * s.zeta()


class DickmanRho(BuiltinFunction):
    r"""
    Dickman's function is the continuous function satisfying the
    differential equation

    .. MATH::

         x \rho'(x) + \rho(x-1) = 0

    with initial conditions `\rho(x)=1` for
    `0 \le x \le 1`. It is useful in estimating the frequency
    of smooth numbers as asymptotically

    .. MATH::

         \Psi(a, a^{1/s}) \sim a \rho(s)

    where `\Psi(a,b)` is the number of `b`-smooth
    numbers less than `a`.

    ALGORITHM:

    Dickmans's function is analytic on the interval
    `[n,n+1]` for each integer `n`. To evaluate at
    `n+t, 0 \le t < 1`, a power series is recursively computed
    about `n+1/2` using the differential equation stated above.
    As high precision arithmetic may be needed for intermediate results
    the computed series are cached for later use.

    Simple explicit formulas are used for the intervals [0,1] and
    [1,2].

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: dickman_rho(2)
        0.306852819440055
        sage: dickman_rho(10)
        2.77017183772596e-11
        sage: dickman_rho(10.00000000000000000000000000000000000000)
        2.77017183772595898875812120063434232634e-11
        sage: plot(log(dickman_rho(x)), (x, 0, 15))                                     # needs sage.plot
        Graphics object consisting of 1 graphics primitive

    AUTHORS:

    - Robert Bradshaw (2008-09)

    REFERENCES:

    - G. Marsaglia, A. Zaman, J. Marsaglia. "Numerical
      Solutions to some Classical Differential-Difference Equations."
      Mathematics of Computation, Vol. 53, No. 187 (1989).
    """
    def __init__(self):
        """
        Construct an object to represent Dickman's rho function.

        TESTS::

            sage: dickman_rho(x)                                                        # needs sage.symbolic
            dickman_rho(x)
            sage: dickman_rho(3)                                                        # needs sage.symbolic
            0.0486083882911316
            sage: dickman_rho(pi)                                                       # needs sage.symbolic
            0.0359690758968463
        """
        self._cur_prec = 0
        BuiltinFunction.__init__(self, "dickman_rho", 1)

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: [dickman_rho(n) for n in [1..10]]                                     # needs sage.symbolic
            [1.00000000000000, 0.306852819440055, 0.0486083882911316,
             0.00491092564776083, 0.000354724700456040, 0.0000196496963539553,
             8.74566995329392e-7, 3.23206930422610e-8, 1.01624828273784e-9,
             2.77017183772596e-11]
            sage: dickman_rho(0)                                                        # needs sage.symbolic
            1.00000000000000
        """
        if not isinstance(x, RealNumber):
            try:
                x = RR(x)
            except (TypeError, ValueError):
                return None
        if x < 0:
            return x.parent()(0)
        elif x <= 1:
            return x.parent()(1)
        elif x <= 2:
            return 1 - x.log()
        n = x.floor()
        if self._cur_prec < x.parent().prec() or n not in self._f:
            self._cur_prec = rel_prec = x.parent().prec()
            # Go a bit beyond so we're not constantly re-computing.
            max = x.parent()(1.1)*x + 10
            abs_prec = (-self.approximate(max).log2() + rel_prec + 2*max.log2()).ceil()
            self._f = {}
            with increase_recursion_limit(int(max)):
                self._compute_power_series(max.floor(), abs_prec, cache_ring=x.parent())
        return self._f[n](2*(x-n-x.parent()(0.5)))

    def power_series(self, n, abs_prec):
        """
        This function returns the power series about `n+1/2` used
        to evaluate Dickman's function. It is scaled such that the interval
        `[n,n+1]` corresponds to `x` in `[-1,1]`.

        INPUT:

        - ``n`` -- the lower endpoint of the interval for which
          this power series holds

        - ``abs_prec`` -- the absolute precision of the
          resulting power series

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: f = dickman_rho.power_series(2, 20); f
            -9.9376e-8*x^11 + 3.7722e-7*x^10 - 1.4684e-6*x^9 + 5.8783e-6*x^8
             - 0.000024259*x^7 + 0.00010341*x^6 - 0.00045583*x^5 + 0.0020773*x^4
             - 0.0097336*x^3 + 0.045224*x^2 - 0.11891*x + 0.13032
            sage: f(-1), f(0), f(1)
            (0.30685, 0.13032, 0.048608)
            sage: dickman_rho(2), dickman_rho(2.5), dickman_rho(3)
            (0.306852819440055, 0.130319561832251, 0.0486083882911316)
        """
        return self._compute_power_series(n, abs_prec, cache_ring=None)

    def _compute_power_series(self, n, abs_prec, cache_ring=None):
        """
        Compute the power series giving Dickman's function on `[n, n+1]`, by
        recursion in `n`. For internal use; ``self.power_series()`` is a wrapper
        around this intended for the user.

        INPUT:

        - ``n`` -- the lower endpoint of the interval for which
          this power series holds

        - ``abs_prec`` -- the absolute precision of the
          resulting power series

        - ``cache_ring`` -- for internal use, caches the power
          series at this precision

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: f = dickman_rho.power_series(2, 20); f
            -9.9376e-8*x^11 + 3.7722e-7*x^10 - 1.4684e-6*x^9 + 5.8783e-6*x^8
             - 0.000024259*x^7 + 0.00010341*x^6 - 0.00045583*x^5 + 0.0020773*x^4
             - 0.0097336*x^3 + 0.045224*x^2 - 0.11891*x + 0.13032
        """
        if n <= 1:
            if n <= -1:
                return PolynomialRealDense(RealField(abs_prec)['x'])
            if n == 0:
                return PolynomialRealDense(RealField(abs_prec)['x'], [1])
            elif n == 1:
                nterms = (RDF(abs_prec) * RDF(2).log()/RDF(3).log()).ceil()
                R = RealField(abs_prec)
                neg_three = ZZ(-3)
                coeffs = [1 - R(1.5).log()] + [neg_three**-k/k for k in range(1, nterms)]
                f = PolynomialRealDense(R['x'], coeffs)
                if cache_ring is not None:
                    self._f[n] = f.truncate_abs(f[0] >> (cache_ring.prec()+1)).change_ring(cache_ring)
                return f
        else:
            f = self._compute_power_series(n-1, abs_prec, cache_ring)
            # integrand = f / (2n+1 + x)
            # We calculate this way because the most significant term is the constant term,
            # and so we want to push the error accumulation and remainder out to the least
            # significant terms.
            integrand = f.reverse().quo_rem(PolynomialRealDense(f.parent(), [1, 2*n+1]))[0].reverse()
            integrand = integrand.truncate_abs(RR(2)**-abs_prec)
            iintegrand = integrand.integral()
            ff = PolynomialRealDense(f.parent(), [f(1) + iintegrand(-1)]) - iintegrand
            i = 0
            while abs(f[i]) < abs(f[i+1]):
                i += 1
            rel_prec = int(abs_prec + abs(RR(f[i])).log2())
            if cache_ring is not None:
                self._f[n] = ff.truncate_abs(ff[0] >> (cache_ring.prec()+1)).change_ring(cache_ring)
            return ff.change_ring(RealField(rel_prec))

    def approximate(self, x, parent=None):
        r"""
        Approximate using de Bruijn's formula.

        .. MATH::

             \rho(x) \sim \frac{exp(-x \xi + Ei(\xi))}{\sqrt{2\pi x}\xi}

        which is asymptotically equal to Dickman's function, and is much
        faster to compute.

        REFERENCES:

        - N. De Bruijn, "The Asymptotic behavior of a function
          occurring in the theory of primes." J. Indian Math Soc. v 15.
          (1951)

        EXAMPLES::

            sage: dickman_rho.approximate(10)                                           # needs sage.rings.real_mpfr
            2.41739196365564e-11
            sage: dickman_rho(10)                                                       # needs sage.symbolic
            2.77017183772596e-11
            sage: dickman_rho.approximate(1000)                                         # needs sage.rings.real_mpfr
            4.32938809066403e-3464
        """
        log, exp, sqrt, pi = math.log, math.exp, math.sqrt, math.pi
        x = float(x)
        xi = log(x)
        y = (exp(xi)-1.0)/xi - x
        while abs(y) > 1e-12:
            dydxi = (exp(xi)*(xi-1.0) + 1.0)/(xi*xi)
            xi -= y/dydxi
            y = (exp(xi)-1.0)/xi - x
        return (-x*xi + RR(xi).eint()).exp() / (sqrt(2*pi*x)*xi)


dickman_rho = DickmanRho()
