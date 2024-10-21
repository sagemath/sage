# sage_setup: distribution = sagemath-categories
r"""
Orthogonal polynomials

Chebyshev polynomials
---------------------

The Chebyshev polynomial of the first kind arises as a solution
to the differential equation

.. MATH::

    (1-x^2)\,y'' - x\,y' + n^2\,y = 0

and those of the second kind as a solution to

.. MATH::

    (1-x^2)\,y'' - 3x\,y' + n(n+2)\,y = 0.

The Chebyshev polynomials of the first kind are defined by the
recurrence relation

.. MATH::

    T_0(x) = 1, \qquad T_1(x) = x, \qquad T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x).

The Chebyshev polynomials of the second kind are defined by the
recurrence relation

.. MATH::

    U_0(x) = 1, \qquad U_1(x) = 2x, \qquad U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x).

For integers `m,n`, they satisfy the orthogonality relations

.. MATH::

    \int_{-1}^1 T_n(x)T_m(x)\,\frac{dx}{\sqrt{1-x^2}} =
    \left\{ \begin{array}{cl} 0 & \text{if } n\neq m, \\ \pi & \text{if } n=m=0, \\ \pi/2 & \text{if } n = m \neq 0, \end{array} \right.

and

.. MATH::

    \int_{-1}^1 U_n(x)U_m(x)\sqrt{1-x^2}\,dx =\frac{\pi}{2}\delta_{m,n}.

They are named after Pafnuty Chebyshev (1821-1894, alternative
transliterations: Tchebyshef or Tschebyscheff).

Hermite polynomials
-------------------

The *Hermite polynomials* are defined either by

.. MATH::

    H_n(x) = (-1)^n e^{x^2/2} \frac{d^n}{dx^n} e^{-x^2/2}

(the "probabilists' Hermite polynomials"), or by

.. MATH::

    H_n(x) = (-1)^n e^{x^2}\frac{d^n}{dx^n}e^{-x^2}

(the "physicists' Hermite polynomials"). Sage (via Maxima) implements
the latter flavor. These satisfy the orthogonality relation

.. MATH::

    \int_{-\infty}^{\infty} H_n(x) H_m(x) \, e^{-x^2} \, dx
    = \sqrt{\pi} n! 2^n \delta_{nm}.

They are named in honor of Charles Hermite (1822-1901), but were first
introduced by Laplace in 1810 and also studied by Chebyshev in 1859.

Legendre polynomials
--------------------

Each *Legendre polynomial* `P_n(x)` is an `n`-th degree polynomial.
It may be expressed using Rodrigues' formula:

.. MATH::

    P_n(x) = (2^n n!)^{-1} {\frac{d^n}{dx^n} } \left[ (x^2 -1)^n \right].

These are solutions to Legendre's differential equation:

.. MATH::

    \frac{d}{dx} \left[ (1-x^2) {\frac{d}{dx}} P(x) \right] + n(n+1)P(x) = 0

and satisfy the orthogonality relation

.. MATH::

    \int_{-1}^{1} P_m(x) P_n(x)\,dx = {\frac{2}{2n + 1}} \delta_{mn}.

The *Legendre function of the second kind* `Q_n(x)` is another
(linearly independent) solution to the Legendre differential equation.
It is not an "orthogonal polynomial" however.

The associated Legendre functions of the first kind `P_\ell^m(x)` can
be given in terms of the "usual" Legendre polynomials by

.. MATH::

    \begin{aligned}
    P_{\ell}^m(x) &= (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}P_\ell(x) \\
    & = \frac{(-1)^m}{2^\ell \ell!} (1-x^2)^{m/2}\frac{d^{\ell+m}}{dx^{\ell+m}}(x^2-1)^{\ell}.
    \end{aligned}

Assuming `0 \le m \le \ell`, they satisfy the orthogonality relation:

.. MATH::

    \int_{-1}^{1} P_k^{(m)} P_{\ell}^{(m)} dx
    = \frac{2(\ell+m)!}{(2\ell+1)(\ell-m)!}\ \delta _{k,\ell},

where `\delta _{k,\ell}` is the Kronecker delta.

The associated Legendre functions of the second kind
`Q_\ell^m(x)` can be given in terms of the "usual"
Legendre polynomials by

.. MATH::

    Q_{\ell}^m(x) = (-1)^m (1-x^2)^{m/2} \frac{d^m}{dx^m} Q_{\ell}(x).

They are named after Adrien-Marie Legendre (1752-1833).

Laguerre polynomials
--------------------

*Laguerre polynomials* may be defined by the Rodrigues formula

.. MATH::

    L_n(x) = \frac{e^x}{n!} \frac{d^n}{dx^n} \left( e^{-x} x^n \right).

They are solutions of Laguerre's equation:

.. MATH::

    x\,y'' + (1 - x)\,y' + n\,y = 0

and satisfy the orthogonality relation

.. MATH::

    \int_0^{\infty} L_m(x) L_n(x) e^{-x} \, dx = \delta_{mn}.

The generalized Laguerre polynomials may be defined by the Rodrigues formula:

.. MATH::

   L_n^{(\alpha)}(x) = \frac{x^{-\alpha} e^x}{n!} \frac{d^n}{dx^n}
   \left(e^{-x} x^{n+\alpha}\right).

(These are also sometimes called the associated Laguerre
polynomials.) The simple Laguerre polynomials are recovered from
the generalized polynomials by setting `\alpha = 0`.

They are named after Edmond Laguerre (1834-1886).

Jacobi polynomials
------------------

*Jacobi polynomials* are a class of orthogonal polynomials. They
are obtained from hypergeometric series in cases where the series
is in fact finite:

.. MATH::

    P_n^{(\alpha,\beta)}(z) = \frac{(\alpha+1)_n}{n!}
    \,_2F_1\left(-n,1+\alpha+\beta+n; \alpha+1; \frac{1-z}{2}\right),

where `()_n` is Pochhammer's symbol (for the rising factorial),
(Abramowitz and Stegun p561.) and thus have the explicit expression

.. MATH::

    P_n^{(\alpha,\beta)} (z) = \frac{\Gamma(\alpha+n+1)}{n!\Gamma(\alpha+\beta+n+1)}
    \sum_{m=0}^n \binom{n}{m} \frac{\Gamma(\alpha+\beta+n+m+1)}{\Gamma(\alpha+m+1)}
    \left(\frac{z-1}{2}\right)^m.

They are named after Carl Gustav Jaboc Jacobi (1804-1851).

Gegenbauer polynomials
----------------------

*Ultraspherical* or *Gegenbauer polynomials* are given in terms of
the Jacobi polynomials `P_n^{(\alpha,\beta)}(x)` with
`\alpha = \beta = a - 1/2` by

.. MATH::

    C_n^{(a)}(x) = \frac{\Gamma(a+1/2)}{\Gamma(2a)}
    \frac{\Gamma(n+2a)}{\Gamma(n+a+1/2)} P_n^{(a-1/2,a-1/2)}(x).

They satisfy the orthogonality relation

.. MATH::

    \int_{-1}^1(1-x^2)^{a-1/2}C_m^{(a)}(x)C_n^{(a)}(x)\, dx
    = \delta_{mn}2^{1-2a}\pi \frac{\Gamma(n+2a)}{(n+a)\Gamma^2(a)\Gamma(n+1)},

for `a > -1/2`. They are obtained from hypergeometric series
in cases where the series is in fact finite:

.. MATH::

    C_n^{(a)}(z) = \frac{(2a)^{\underline{n}}}{n!}
    \,_2F_1\left(-n,2a+n; a+\frac{1}{2}; \frac{1-z}{2}\right)

where `\underline{n}` is the falling factorial. (See
Abramowitz and Stegun p561.)

They are named for Leopold Gegenbauer (1849-1903).

Krawtchouk polynomials
----------------------

The *Krawtchouk polynomials* are discrete orthogonal polynomials that
are given by the hypergeometric series

.. MATH::

    K_j(x; n, p) = (-1)^j \binom{n}{j} p^j
    \,_{2}F_1\left(-j,-x; -n; p^{-1}\right).

Since they are discrete orthogonal polynomials, they satisfy an orthogonality
relation defined on a discrete (in this case finite) set of points:

.. MATH::

    \sum_{m=0}^n K_i(m; n, p) K_j(m; n, p) \, \binom{n}{m} p^m q^{n-m}
    = \binom{n}{j} (pq)^j \delta_{ij},

where `q = 1 - p`. They can also be described by the recurrence relation

.. MATH::

    j K_j(x; n, p) = (x - (n-j+1) p - (j-1) q) K_{j-1}(x; n, p)
    - p q (n - j + 2) K_{j-2}(x; n, p),

where `K_0(x; n, p) = 1` and `K_1(x; n, p) = x - n p`.

They are named for Mykhailo Krawtchouk (1892-1942).


Meixner polynomials
-------------------

The *Meixner polynomials* are discrete orthogonal polynomials that
are given by the hypergeometric series

.. MATH::

    M_n(x; n, p) = (-1)^j \binom{n}{j} p^j
    \,_{2}F_1\left(-j,-x; -n; p^{-1}\right).

They satisfy an orthogonality relation:

.. MATH::

    \sum_{k=0}^{\infty} \tilde{M}_n(k; b, c) \tilde{M}_m(k; b, c) \, \frac{(b)_k}{k!} c^k
    = \frac{c^{-n} n!}{(b)_n (1-c)^b} \delta_{mn},

where `\tilde{M}_n(x; b, c) = M_n(x; b, c) / (b)_x`, for `b > 0 ` and
`0 < c < 1`. They can also be described by the recurrence relation

.. MATH::

   \begin{aligned}
    c (n-1+b) M_n(x; b, c) & = ((c-1) x + n-1 + c (n-1+b)) (b+n-1) M_{n-1}(x; b, c)
    \\ & \qquad - (b+n-1) (b+n-2) (n-1) M_{n-2}(x; b, c),
    \end{aligned}

where `M_0(x; b, c) = 0` and `M_1(x; b, c) = (1 - c^{-1}) x + b`.

They are named for Josef Meixner (1908-1994).

Hahn polynomials
----------------

The *Hahn polynomials* are discrete orthogonal polynomials that
are given by the hypergeometric series

.. MATH::

    Q_k(x; a, b, n) = \,_{3}F_2\left(-k,k+a+b+1,-x; a+1,-n; 1\right).

They satisfy an orthogonality relation:

.. MATH::

    \sum_{k=0}^{n-1} Q_i(k; a, b, n) Q_j(k; a, b, n) \, \rho(k)
    = \frac{\delta_{ij}}{\pi_i},

where

.. MATH::

    \begin{aligned}
    \rho(k) &= \binom{a+k}{k} \binom{b+n-k}{n-k},
    \\
    \pi_i &= \delta_{ij} \frac{(-1)^i i! (b+1)_i (i+a+b+1)_{n+1}}{n! (2i+a+b+1) (-n)_i (a+1)_i}.
    \end{aligned}

They can also be described by the recurrence relation

.. MATH::

    A Q_k(x; a,b,n) = (-x + A + C) Q_{k-1}(x; a,b,n) - C Q_{k-2}(x; a,b,n),

where `Q_0(x; a,b,n) = 1` and `Q_1(x; a,b,n) = 1 - \frac{a+b+2}{(a+1)n} x` and

.. MATH::

    A = \frac{(k+a+b) (k+a) (n-k+1)}{(2k+a+b-1) (2k+a+b)},
    \qquad
    C = \frac{(k-1) (k+b-1) (k+a+b+n)}{(2k+a+b-2) (2k+a+b-1)}.

They are named for Wolfgang Hahn (1911-1998), although they were first
introduced by Chebyshev in 1875.


Pochhammer symbol
-----------------

For completeness, the *Pochhammer symbol*, introduced by Leo August
Pochhammer, `(x)_n`, is used in the theory of special
functions to represent the "rising factorial" or "upper factorial"

.. MATH::

    (x)_n = x(x+1)(x+2) \cdots (x+n-1) = \frac{(x+n-1)!}{(x-1)!}.

On the other hand, the *falling factorial* or *lower factorial* is

.. MATH::

    x^{\underline{n}} = \frac{x!}{(x-n)!},

in the notation of Ronald L. Graham, Donald E. Knuth and Oren Patashnik
in their book Concrete Mathematics.

.. TODO::

    Implement Zernike polynomials.
    :wikipedia:`Zernike_polynomials`

REFERENCES:

- [AS1964]_
- :wikipedia:`Chebyshev_polynomials`
- :wikipedia:`Legendre_polynomials`
- :wikipedia:`Hermite_polynomials`
- http://mathworld.wolfram.com/GegenbauerPolynomial.html
- :wikipedia:`Jacobi_polynomials`
- :wikipedia:`Laguerre_polynomia`
- :wikipedia:`Associated_Legendre_polynomials`
- :wikipedia:`Kravchuk_polynomials`
- :wikipedia:`Meixner_polynomials`
- :wikipedia:`Hahn_polynomials`
- Roelof Koekeok and René F. Swarttouw, :arxiv:`math/9602214`
- [Koe1999]_

AUTHORS:

- David Joyner (2006-06)
- Stefan Reiterer (2010-)
- Ralf Stephan (2015-)

The original module wrapped some of the orthogonal/special functions
in the Maxima package "orthopoly" and was written by Barton
Willis of the University of Nebraska at Kearney.
"""
# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdj@usna.edu>
#                     2010 Stefan Reiterer <maldun.finsterschreck@gmail.com>
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

import warnings

import sage.rings.abc

from sage.arith.misc import rising_factorial
from sage.misc.lazy_import import lazy_import
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.symbolic.function import BuiltinFunction, GinacFunction
from sage.structure.element import Expression, parent

lazy_import('sage.functions.other', ['factorial', 'binomial'])

lazy_import('sage.misc.latex', 'latex')
lazy_import('sage.rings.cc', 'CC')
lazy_import('sage.rings.polynomial.polynomial_ring_constructor', 'PolynomialRing')
lazy_import('sage.rings.real_mpfr', 'RR')

lazy_import('sage.symbolic.ring', 'SR')
lazy_import('sage.calculus.calculus', 'maxima', as_='_maxima')

lazy_import('sage.libs.mpmath.sage_utils', 'call', as_='_mpmath_call')
lazy_import('sage.libs.mpmath.all', 'chebyt', as_='_mpmath_chebyt')
lazy_import('sage.libs.mpmath.all', 'chebyu', as_='_mpmath_chebyu')
lazy_import('sage.libs.mpmath.all', 'laguerre', as_='_mpmath_laguerre')
lazy_import('sage.libs.mpmath.all', 'legenp', as_='_mpmath_legenp')
lazy_import('sage.libs.mpmath.all', 'legenq', as_='_mpmath_legenq')

lazy_import('scipy.special', 'eval_chebyu', as_='_scipy_chebyu')


class OrthogonalFunction(BuiltinFunction):
    """
    Base class for orthogonal polynomials.

    This class is an abstract base class for all orthogonal polynomials since
    they share similar properties. The evaluation as a polynomial
    is either done via maxima, or with pynac.

    Convention: The first argument is always the order of the polynomial,
    the others are other values or parameters where the polynomial is
    evaluated.
    """
    def __init__(self, name, nargs=2, latex_name=None, conversions=None):
        """
        :class:`OrthogonalFunction` class needs the same input parameter as
        its parent class.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import OrthogonalFunction
            sage: new = OrthogonalFunction('testo_P')
            sage: new
            testo_P
        """
        self._maxima_name = None
        if conversions:
            try:
                self._maxima_name = conversions['maxima']
            except KeyError:
                pass
        super().__init__(name=name, nargs=nargs,
                                                 latex_name=latex_name,
                                                 conversions=conversions)

    def eval_formula(self, *args):
        """
        Evaluate this polynomial using an explicit formula.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import OrthogonalFunction
            sage: P = OrthogonalFunction('testo_P')
            sage: P.eval_formula(1, 2.0)
            Traceback (most recent call last):
            ...
            NotImplementedError: no explicit calculation of values implemented
        """
        raise NotImplementedError("no explicit calculation of values implemented")

    def _eval_special_values_(self, *args):
        """
        Evaluate the polynomial explicitly for special values.

        EXAMPLES::

            sage: var('n')                                                              # needs sage.symbolic
            n
            sage: chebyshev_T(n, -1)                                                    # needs sage.symbolic
            (-1)^n
        """
        raise ValueError("no special values known")

    def _eval_(self, n, *args):
        """
        The :meth:`_eval_()` method decides which evaluation suits best
        for the given input, and returns a proper value.

        EXAMPLES::

            sage: var('n,x')                                                            # needs sage.symbolic
            (n, x)
            sage: chebyshev_T(5, x)                                                     # needs sage.symbolic
            16*x^5 - 20*x^3 + 5*x
        """
        return None

    def __call__(self, *args, **kwds):
        """
        This overides the call method from SageObject to avoid problems with coercions,
        since the _eval_ method is able to handle more data types than symbolic functions
        would normally allow.
        Thus we have the distinction between algebraic objects (if n is an integer),
        and else as symbolic function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: chebyshev_T(5, x)
            16*x^5 - 20*x^3 + 5*x
            sage: chebyshev_T(5, x, algorithm='pari')                                   # needs sage.libs.pari
            16*x^5 - 20*x^3 + 5*x
            sage: chebyshev_T(5, x, algorithm='maxima')
            16*x^5 - 20*x^3 + 5*x
            sage: chebyshev_T(5, x, algorithm='recursive')
            16*x^5 - 20*x^3 + 5*x
        """
        algorithm = kwds.get('algorithm', None)
        if algorithm == 'pari':
            return self.eval_pari(*args, **kwds)
        elif algorithm == 'recursive':
            return self.eval_recursive(*args, **kwds)
        elif algorithm == 'maxima':
            kwds['hold'] = True
            return _maxima(self._eval_(*args, **kwds))._sage_()

        return super().__call__(*args, **kwds)


class ChebyshevFunction(OrthogonalFunction):
    """
    Abstract base class for Chebyshev polynomials of the first and second kind.

    EXAMPLES::

        sage: chebyshev_T(3, x)                                                         # needs sage.symbolic
        4*x^3 - 3*x
    """
    def __call__(self, n, *args, **kwds):
        """
        This overides the call method from :class:`SageObject` to avoid problems with coercions,
        since the ``_eval_`` method is able to handle more data types than symbolic functions
        would normally allow.
        Thus we have the distinction between algebraic objects (if n is an integer),
        and else as symbolic function.

        EXAMPLES::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = NumberField(x^3 - x - 1)                                      # needs sage.rings.number_field
            sage: chebyshev_T(5, a)                                                     # needs sage.rings.number_field
            16*a^2 + a - 4
            sage: chebyshev_T(5, MatrixSpace(ZZ, 2)([1, 2, -4, 7]))                     # needs sage.modules
            [-40799  44162]
            [-88324  91687]
            sage: R.<x> = QQ[]
            sage: parent(chebyshev_T(5, x))
            Univariate Polynomial Ring in x over Rational Field
            sage: chebyshev_T(5, 2, hold=True)                                          # needs sage.symbolic
            chebyshev_T(5, 2)
            sage: chebyshev_T(1, 2, 3)
            Traceback (most recent call last):
            ...
            TypeError: Symbolic function chebyshev_T takes exactly 2 arguments (3 given)
        """
        # If n is an integer: consider the polynomial as an algebraic (not symbolic) object
        if n in ZZ and not kwds.get('hold', False):
            try:
                return self._eval_(n, *args)
            except Exception:
                pass

        return super().__call__(n, *args, **kwds)

    def _eval_(self, n, x):
        """
        The :meth:`_eval_()` method decides which evaluation suits best
        for the given input, and returns a proper value.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: var('n,x')
            (n, x)
            sage: chebyshev_T(5, x)
            16*x^5 - 20*x^3 + 5*x
            sage: chebyshev_T(64, x)
            2*(2*(2*(2*(2*(2*x^2 - 1)^2 - 1)^2 - 1)^2 - 1)^2 - 1)^2 - 1
            sage: chebyshev_T(n, -1)
            (-1)^n
            sage: chebyshev_T(-7, x)
            64*x^7 - 112*x^5 + 56*x^3 - 7*x
            sage: chebyshev_T(3/2, x)
            chebyshev_T(3/2, x)

            sage: R.<t> = QQ[]
            sage: chebyshev_T(2, t)
            2*t^2 - 1
            sage: chebyshev_U(2, t)
            4*t^2 - 1
            sage: parent(chebyshev_T(4, RIF(5)))                                        # needs sage.rings.real_interval_field
            Real Interval Field with 53 bits of precision
            sage: RR2 = RealField(5)                                                    # needs sage.rings.real_mpfr
            sage: chebyshev_T(100000, RR2(2))                                           # needs sage.rings.real_mpfr
            8.9e57180
            sage: chebyshev_T(5, Qp(3)(2))                                              # needs sage.rings.padics
            2 + 3^2 + 3^3 + 3^4 + 3^5 + O(3^20)
            sage: chebyshev_T(100001/2, 2)                                              # needs sage.symbolic
            ...chebyshev_T(100001/2, 2)
            sage: chebyshev_U._eval_(1.5, Mod(8,9)) is None                             # needs mpmath
            True
        """
        # n is an integer => evaluate algebraically (as polynomial)
        if n in ZZ:
            n = ZZ(n)
            # Expanded symbolic expression only for small values of n
            if isinstance(x, Expression) and n.abs() < 32:
                return self.eval_formula(n, x)
            return self.eval_algebraic(n, x)

        if isinstance(x, Expression) or isinstance(n, Expression):
            # Check for known identities
            try:
                return self._eval_special_values_(n, x)
            except ValueError:
                # Don't evaluate => keep symbolic
                return None

        # n is not an integer and neither n nor x is symbolic.
        # We assume n and x are real/complex and evaluate numerically
        import sage.libs.mpmath.all as mpmath
        try:
            return self._evalf_(n, x)
        except mpmath.NoConvergence:
            warnings.warn("mpmath failed, keeping expression unevaluated",
                          RuntimeWarning)
            return None
        except Exception:
            # Numerical evaluation failed => keep symbolic
            return None


class Func_chebyshev_T(ChebyshevFunction):
    """
    Chebyshev polynomials of the first kind.

    REFERENCE:

    - [AS1964]_ 22.5.31 page 778 and 6.1.22 page 256.

    EXAMPLES::

       sage: chebyshev_T(5, x)                                                          # needs sage.symbolic
       16*x^5 - 20*x^3 + 5*x
       sage: var('k')                                                                   # needs sage.symbolic
       k
       sage: test = chebyshev_T(k, x); test                                             # needs sage.symbolic
       chebyshev_T(k, x)
    """
    def __init__(self):
        """
        Init method for the chebyshev polynomials of the first kind.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: var('n, x')
            (n, x)
            sage: from sage.functions.orthogonal_polys import Func_chebyshev_T
            sage: chebyshev_T2 = Func_chebyshev_T()
            sage: chebyshev_T2(1, x)
            x
            sage: chebyshev_T(x, x)._sympy_()                                           # needs sympy
            chebyshevt(x, x)
            sage: maxima(chebyshev_T(1, x, hold=True))
            _SAGE_VAR_x
            sage: maxima(chebyshev_T(n, chebyshev_T(n, x)))
            chebyshev_t(_SAGE_VAR_n,chebyshev_t(_SAGE_VAR_n,_SAGE_VAR_x))
        """
        ChebyshevFunction.__init__(self, 'chebyshev_T', nargs=2,
                                   conversions=dict(maxima='chebyshev_t',
                                                    mathematica='ChebyshevT',
                                                    sympy='chebyshevt',
                                                    giac='tchebyshev1'))

    def _latex_(self):
        r"""
        TESTS::

            sage: latex(chebyshev_T)
            T_n
        """
        return r"T_n"

    def _print_latex_(self, n, z):
        r"""
        TESTS::

            sage: latex(chebyshev_T(3, x, hold=True))                                   # needs sage.symbolic
            T_{3}\left(x\right)
        """
        return r"T_{{{}}}\left({}\right)".format(latex(n), latex(z))

    def _eval_special_values_(self, n, x):
        """
        Values known for special values of x.
        For details see [AS1964]_ 22.4 (p. 777)

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: var('n')
            n
            sage: chebyshev_T(n, 1)
            1
            sage: chebyshev_T(n, 0)
            1/2*(-1)^(1/2*n)*((-1)^n + 1)
            sage: chebyshev_T(n, -1)
            (-1)^n
            sage: chebyshev_T._eval_special_values_(3/2, x)
            Traceback (most recent call last):
            ...
            ValueError: no special value found
            sage: chebyshev_T._eval_special_values_(n, 0.1)
            Traceback (most recent call last):
            ...
            ValueError: no special value found
        """
        if x == 1:
            return x

        if x == -1:
            return x**n

        if x == 0:
            return (1+(-1)**n)*(-1)**(n/2)/2

        raise ValueError("no special value found")

    def _evalf_(self, n, x, **kwds):
        """
        Evaluates :class:`chebyshev_T` numerically with mpmath.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: chebyshev_T._evalf_(10, 3)
            2.26195370000000e7
            sage: chebyshev_T._evalf_(10, 3, parent=RealField(75))
            2.261953700000000000000e7
            sage: chebyshev_T._evalf_(10, I)                                            # needs sage.symbolic
            -3363.00000000000
            sage: chebyshev_T._evalf_(5, 0.3)
            0.998880000000000
            sage: chebyshev_T(1/2, 0)
            0.707106781186548
            sage: chebyshev_T(1/2, 3/2)
            1.11803398874989
            sage: chebyshev_T._evalf_(1.5, Mod(8,9))
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate chebyshev_T with parent Ring of integers modulo 9

        This simply evaluates using :class:`RealField` or :class:`ComplexField`::

            sage: chebyshev_T(1234.5, RDF(2.1))                                         # needs sage.rings.real_mpfr
            5.48174256255782e735
            sage: chebyshev_T(1234.5, I)                                                # needs sage.rings.real_mpfr sage.symbolic
            -1.21629397684152e472 - 1.21629397684152e472*I

        For large values of ``n``, mpmath fails (but the algebraic formula
        still works)::

            sage: chebyshev_T._evalf_(10^6, 0.1)                                        # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            NoConvergence: Hypergeometric series converges too slowly.
            Try increasing maxterms.
            sage: chebyshev_T(10^6, 0.1)                                                # needs sage.rings.real_mpfr
            0.636384327171504
        """
        try:
            real_parent = kwds['parent']
        except KeyError:
            real_parent = parent(x)

            if not isinstance(real_parent, (sage.rings.abc.RealField, sage.rings.abc.ComplexField)):
                # parent is not a real or complex field: figure out a good parent
                if x in RR:
                    x = RR(x)
                    real_parent = RR
                elif x in CC:
                    x = CC(x)
                    real_parent = CC

        if not isinstance(real_parent, (sage.rings.abc.RealField, sage.rings.abc.ComplexField)):
            raise TypeError("cannot evaluate chebyshev_T with parent {}".format(real_parent))

        return _mpmath_call(_mpmath_chebyt, n, x, parent=real_parent)

    def eval_formula(self, n, x):
        """
        Evaluate ``chebyshev_T`` using an explicit formula.
        See [AS1964]_ 227 (p. 782) for details for the recursions.
        See also [Koe1999]_ for fast evaluation techniques.

        INPUT:

        - ``n`` -- integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: chebyshev_T.eval_formula(-1, x)
            x
            sage: chebyshev_T.eval_formula(0, x)
            1
            sage: chebyshev_T.eval_formula(1, x)
            x
            sage: chebyshev_T.eval_formula(10, x)
            512*x^10 - 1280*x^8 + 1120*x^6 - 400*x^4 + 50*x^2 - 1
            sage: chebyshev_T.eval_algebraic(10, x).expand()
            512*x^10 - 1280*x^8 + 1120*x^6 - 400*x^4 + 50*x^2 - 1

            sage: chebyshev_T.eval_formula(2, 0.1) == chebyshev_T._evalf_(2, 0.1)       # needs sage.rings.complex_double
            True
        """
        if n < 0:
            return self.eval_formula(-n, x)
        elif n == 0:
            return parent(x).one()

        res = parent(x).zero()
        for j in range(n // 2 + 1):
            f = factorial(n-1-j) / factorial(j) / factorial(n-2*j)
            res += (-1)**j * (2*x)**(n-2*j) * f
        res *= n/2
        return res

    def eval_algebraic(self, n, x):
        """
        Evaluate :class:`chebyshev_T` as polynomial, using a recursive
        formula.

        INPUT:

        - ``n`` -- integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

            sage: chebyshev_T.eval_algebraic(5, x)                                      # needs sage.symbolic
            2*(2*(2*x^2 - 1)*x - x)*(2*x^2 - 1) - x
            sage: chebyshev_T(-7, x) - chebyshev_T(7, x)                                # needs sage.symbolic
            0
            sage: R.<t> = ZZ[]
            sage: chebyshev_T.eval_algebraic(-1, t)
            t
            sage: chebyshev_T.eval_algebraic(0, t)
            1
            sage: chebyshev_T.eval_algebraic(1, t)
            t
            sage: chebyshev_T(7^100, 1/2)
            1/2
            sage: chebyshev_T(7^100, Mod(2,3))
            2
            sage: n = 97; x = RIF(pi/2/n)                                               # needs sage.symbolic
            sage: chebyshev_T(n, cos(x)).contains_zero()                                # needs sage.symbolic
            True

            sage: # needs sage.rings.padics
            sage: R.<t> = Zp(2, 8, 'capped-abs')[]
            sage: chebyshev_T(10^6 + 1, t)
            (2^7 + O(2^8))*t^5 + O(2^8)*t^4 + (2^6 + O(2^8))*t^3 + O(2^8)*t^2
             + (1 + 2^6 + O(2^8))*t + O(2^8)
        """
        if n == 0:
            return parent(x).one()
        if n < 0:
            return self._eval_recursive_(-n, x)[0]
        return self._eval_recursive_(n, x)[0]

    def _eval_recursive_(self, n, x, both=False):
        """
        If ``both=True``, compute ``(T(n,x), T(n-1,x))`` using a
        recursive formula.
        If ``both=False``, return instead a tuple ``(T(n,x), False)``.

        EXAMPLES::

            sage: chebyshev_T._eval_recursive_(5, x)                                    # needs sage.symbolic
            (2*(2*(2*x^2 - 1)*x - x)*(2*x^2 - 1) - x, False)
            sage: chebyshev_T._eval_recursive_(5, x, True)                              # needs sage.symbolic
            (2*(2*(2*x^2 - 1)*x - x)*(2*x^2 - 1) - x, 2*(2*x^2 - 1)^2 - 1)
        """
        if n == 1:
            return x, parent(x).one()

        assert n >= 2
        a, b = self._eval_recursive_((n+1)//2, x, both or n % 2)
        if n % 2 == 0:
            return 2*a*a - 1, both and 2*a*b - x
        else:
            return 2*a*b - x, both and 2*b*b - 1

    def _eval_numpy_(self, n, x):
        """
        Evaluate ``self`` using numpy.

        EXAMPLES::

            sage: # needs numpy scipy
            sage: import numpy
            sage: z = numpy.array([1,2])
            sage: z2 = numpy.array([[1,2],[1,2]])
            sage: z3 = numpy.array([1,2,3.])
            sage: chebyshev_T(1,z)
            array([1., 2.])
            sage: chebyshev_T(1,z2)
            array([[1., 2.],
                   [1., 2.]])
            sage: chebyshev_T(1,z3)
            array([1., 2., 3.])
            sage: chebyshev_T(z,0.1)
            array([ 0.1 , -0.98])
        """
        from scipy.special import eval_chebyt
        return eval_chebyt(n, x)

    def _derivative_(self, n, x, diff_param):
        """
        Return the derivative of :class:`chebyshev_T` in form of the Chebyshev
        polynomial of the second kind :class:`chebyshev_U`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: var('k')
            k
            sage: derivative(chebyshev_T(k, x), x)
            k*chebyshev_U(k - 1, x)
            sage: derivative(chebyshev_T(3, x), x)
            12*x^2 - 3
            sage: derivative(chebyshev_T(k, x), k)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative w.r.t. to the index is not supported yet
        """
        if diff_param == 0:
            raise NotImplementedError("derivative w.r.t. to the index is not supported yet")
        elif diff_param == 1:
            return n*chebyshev_U(n-1, x)
        raise ValueError("illegal differentiation parameter {}".format(diff_param))


chebyshev_T = Func_chebyshev_T()


class Func_chebyshev_U(ChebyshevFunction):
    """
    Class for the Chebyshev polynomial of the second kind.

    REFERENCE:

    - [AS1964]_ 22.8.3 page 783 and 6.1.22 page 256.

    EXAMPLES::

        sage: R.<t> = QQ[]
        sage: chebyshev_U(2, t)
        4*t^2 - 1
        sage: chebyshev_U(3, t)
        8*t^3 - 4*t
    """
    def __init__(self):
        """
        Init method for the chebyshev polynomials of the second kind.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: var('n, x')
            (n, x)
            sage: from sage.functions.orthogonal_polys import Func_chebyshev_U
            sage: chebyshev_U2 = Func_chebyshev_U()
            sage: chebyshev_U2(1, x)
            2*x
            sage: chebyshev_U(x, x)._sympy_()                                           # needs sympy
            chebyshevu(x, x)
            sage: maxima(chebyshev_U(2,x, hold=True))
            3*(...-...(8*(1-_SAGE_VAR_x))/3)+(4*(1-_SAGE_VAR_x)^2)/3+1)
            sage: maxima(chebyshev_U(n,x, hold=True))
            chebyshev_u(_SAGE_VAR_n,_SAGE_VAR_x)
        """
        ChebyshevFunction.__init__(self, 'chebyshev_U', nargs=2,
                                   conversions=dict(maxima='chebyshev_u',
                                                    mathematica='ChebyshevU',
                                                    sympy='chebyshevu',
                                                    giac='tchebyshev2'))

    def _latex_(self):
        r"""
        TESTS::

            sage: latex(chebyshev_U)
            U_n
        """
        return r"U_n"

    def _print_latex_(self, n, z):
        r"""
        TESTS::

            sage: latex(chebyshev_U(3, x, hold=True))                                   # needs sage.symbolic
            U_{3}\left(x\right)
        """
        return r"U_{{{}}}\left({}\right)".format(latex(n), latex(z))

    def eval_formula(self, n, x):
        """
        Evaluate ``chebyshev_U`` using an explicit formula.

        See [AS1964]_ 227 (p. 782) for details on the recursions.
        See also [Koe1999]_ for the recursion formulas.

        INPUT:

        - ``n`` -- integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: chebyshev_U.eval_formula(10, x)
            1024*x^10 - 2304*x^8 + 1792*x^6 - 560*x^4 + 60*x^2 - 1
            sage: chebyshev_U.eval_formula(-2, x)
            -1
            sage: chebyshev_U.eval_formula(-1, x)
            0
            sage: chebyshev_U.eval_formula(0, x)
            1
            sage: chebyshev_U.eval_formula(1, x)
            2*x
            sage: chebyshev_U.eval_formula(2, 0.1) == chebyshev_U._evalf_(2, 0.1)
            True
        """
        if n < -1:
            return -self.eval_formula(-n-2, x)

        res = parent(x).zero()
        for j in range(n // 2 + 1):
            f = binomial(n-j, j)
            res += (-1)**j * (2*x)**(n-2*j) * f
        return res

    def eval_algebraic(self, n, x):
        """
        Evaluate :class:`chebyshev_U` as polynomial, using a recursive
        formula.

        INPUT:

        - ``n`` -- integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

            sage: chebyshev_U.eval_algebraic(5, x)                                      # needs sage.symbolic
            -2*((2*x + 1)*(2*x - 1)*x - 4*(2*x^2 - 1)*x)*(2*x + 1)*(2*x - 1)
            sage: parent(chebyshev_U(3, Mod(8,9)))
            Ring of integers modulo 9
            sage: parent(chebyshev_U(3, Mod(1,9)))
            Ring of integers modulo 9
            sage: chebyshev_U(-3, x) + chebyshev_U(1, x)                                # needs sage.symbolic
            0
            sage: chebyshev_U(-1, Mod(5,8))
            0
            sage: parent(chebyshev_U(-1, Mod(5,8)))
            Ring of integers modulo 8
            sage: R.<t> = ZZ[]
            sage: chebyshev_U.eval_algebraic(-2, t)
            -1
            sage: chebyshev_U.eval_algebraic(-1, t)
            0
            sage: chebyshev_U.eval_algebraic(0, t)
            1
            sage: chebyshev_U.eval_algebraic(1, t)
            2*t
            sage: n = 97; x = RIF(pi/n)                                                 # needs sage.symbolic
            sage: chebyshev_U(n - 1, cos(x)).contains_zero()                            # needs sage.symbolic
            True

            sage: # needs sage.rings.padics
            sage: R.<t> = Zp(2, 6, 'capped-abs')[]
            sage: chebyshev_U(10^6 + 1, t)
            (2 + O(2^6))*t + O(2^6)
        """
        if n == -1:
            return parent(x).zero()
        if n < 0:
            return -self._eval_recursive_(-n-2, x)[0]
        return self._eval_recursive_(n, x)[0]

    def _eval_recursive_(self, n, x, both=False):
        """
        If ``both=True``, compute ``(U(n,x), U(n-1,x))`` using a
        recursive formula.
        If ``both=False``, return instead a tuple ``(U(n,x), False)``.

        EXAMPLES::

            sage: chebyshev_U._eval_recursive_(3, x)                                    # needs sage.symbolic
            (4*((2*x + 1)*(2*x - 1) - 2*x^2)*x, False)
            sage: chebyshev_U._eval_recursive_(3, x, True)                              # needs sage.symbolic
            (4*((2*x + 1)*(2*x - 1) - 2*x^2)*x,
             ((2*x + 1)*(2*x - 1) + 2*x)*((2*x + 1)*(2*x - 1) - 2*x))
        """
        if n == 0:
            return parent(x).one(), 2*x

        assert n >= 1
        a, b = self._eval_recursive_((n-1)//2, x, True)
        if n % 2 == 0:
            return (b+a)*(b-a), both and 2*b*(x*b-a)
        else:
            return 2*a*(b-x*a), both and (b+a)*(b-a)

    def _evalf_(self, n, x, **kwds):
        """
        Evaluate :class:`chebyshev_U` numerically with mpmath.

        EXAMPLES::

            sage: chebyshev_U(5, -4 + 3.*I)                                             # needs sage.symbolic
            98280.0000000000 - 11310.0000000000*I
            sage: chebyshev_U(10, 3).n(75)                                              # needs sage.symbolic
            4.661117900000000000000e7
            sage: chebyshev_U._evalf_(1.5, Mod(8,9))                                    # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate chebyshev_U with parent Ring of integers modulo 9
        """
        try:
            real_parent = kwds['parent']
        except KeyError:
            real_parent = parent(x)

            if not isinstance(real_parent, (sage.rings.abc.RealField, sage.rings.abc.ComplexField)):
                # parent is not a real or complex field: figure out a good parent
                if x in RR:
                    x = RR(x)
                    real_parent = RR
                elif x in CC:
                    x = CC(x)
                    real_parent = CC

        if not isinstance(real_parent, (sage.rings.abc.RealField, sage.rings.abc.ComplexField)):
            raise TypeError("cannot evaluate chebyshev_U with parent {}".format(real_parent))

        return _mpmath_call(_mpmath_chebyu, n, x, parent=real_parent)

    def _eval_special_values_(self, n, x):
        """
        Values known for special values of x.
        See [AS1964]_ 22.4 (p.777).

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: var('n')
            n
            sage: chebyshev_U(n, 1)
            n + 1
            sage: chebyshev_U(n, 0)
            1/2*(-1)^(1/2*n)*((-1)^n + 1)
            sage: chebyshev_U(n, -1)
            (-1)^n*(n + 1)
            sage: chebyshev_U._eval_special_values_(n, 2)
            Traceback (most recent call last):
            ...
            ValueError: no special value found
        """
        if x == 1:
            return x*(n+1)

        if x == -1:
            return x**n*(n+1)

        if x == 0:
            return (1+(-1)**n)*(-1)**(n/2)/2

        raise ValueError("no special value found")

    def _eval_numpy_(self, n, x):
        """
        Evaluate ``self`` using numpy.

        EXAMPLES::

            sage: # needs numpy scipy
            sage: import numpy
            sage: z = numpy.array([1,2])
            sage: z2 = numpy.array([[1,2],[1,2]])
            sage: z3 = numpy.array([1,2,3.])
            sage: chebyshev_U(1,z)
            array([2., 4.])
            sage: chebyshev_U(1,z2)
            array([[2., 4.],
                   [2., 4.]])
            sage: chebyshev_U(1,z3)
            array([2., 4., 6.])
            sage: chebyshev_U(z,0.1)
            array([ 0.2 , -0.96])
        """
        return _scipy_chebyu(n, x)

    def _derivative_(self, n, x, diff_param):
        """
        Return the derivative of :class:`chebyshev_U` in form of the Chebyshev
        polynomials of the first and second kind.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: var('k')
            k
            sage: derivative(chebyshev_U(k,x), x)
            ((k + 1)*chebyshev_T(k + 1, x) - x*chebyshev_U(k, x))/(x^2 - 1)
            sage: derivative(chebyshev_U(3,x), x)
            24*x^2 - 4
            sage: derivative(chebyshev_U(k,x), k)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative w.r.t. to the index is not supported yet
        """
        if diff_param == 0:
            raise NotImplementedError("derivative w.r.t. to the index is not supported yet")
        elif diff_param == 1:
            return ((n+1)*chebyshev_T(n+1, x) - x*chebyshev_U(n, x)) / (x*x-1)
        raise ValueError("illegal differentiation parameter {}".format(diff_param))


chebyshev_U = Func_chebyshev_U()


class Func_legendre_P(GinacFunction):
    r"""
    EXAMPLES::

        sage: # needs sage.symbolic
        sage: legendre_P(4, 2.0)
        55.3750000000000
        sage: legendre_P(1, x)
        x
        sage: legendre_P(4, x + 1)
        35/8*(x + 1)^4 - 15/4*(x + 1)^2 + 3/8
        sage: legendre_P(1/2, I+1.)
        1.05338240025858 + 0.359890322109665*I
        sage: legendre_P(0, SR(1)).parent()
        Symbolic Ring

        sage: legendre_P(0, 0)                                                          # needs sage.symbolic
        1
        sage: legendre_P(1, x)                                                          # needs sage.symbolic
        x

        sage: # needs sage.symbolic
        sage: legendre_P(4, 2.)
        55.3750000000000
        sage: legendre_P(5.5, 1.00001)
        1.00017875754114
        sage: legendre_P(1/2, I + 1).n()
        1.05338240025858 + 0.359890322109665*I
        sage: legendre_P(1/2, I + 1).n(59)
        1.0533824002585801 + 0.35989032210966539*I
        sage: legendre_P(42, RR(12345678))
        2.66314881466753e309
        sage: legendre_P(42, Reals(20)(12345678))
        2.6632e309
        sage: legendre_P(201/2, 0).n()
        0.0561386178630179
        sage: legendre_P(201/2, 0).n(100)
        0.056138617863017877699963095883

        sage: # needs sage.symbolic
        sage: R.<x> = QQ[]
        sage: legendre_P(4, x)
        35/8*x^4 - 15/4*x^2 + 3/8
        sage: legendre_P(10000, x).coefficient(x, 1)
        0
        sage: var('t,x')
        (t, x)
        sage: legendre_P(-5, t)
        35/8*t^4 - 15/4*t^2 + 3/8
        sage: legendre_P(4, x + 1)
        35/8*(x + 1)^4 - 15/4*(x + 1)^2 + 3/8
        sage: legendre_P(4, sqrt(2))
        83/8
        sage: legendre_P(4, I*e)
        35/8*e^4 + 15/4*e^2 + 3/8

        sage: # needs sage.symbolic
        sage: n = var('n')
        sage: derivative(legendre_P(n,x), x)
        (n*x*legendre_P(n, x) - n*legendre_P(n - 1, x))/(x^2 - 1)
        sage: derivative(legendre_P(3,x), x)
        15/2*x^2 - 3/2
        sage: derivative(legendre_P(n,x), n)
        Traceback (most recent call last):
        ...
        RuntimeError: derivative w.r.t. to the index is not supported yet

    TESTS:

    Verify that :issue:`33962` is fixed::

        sage: [legendre_P(n, 0) for n in range(-10, 10)]                                # needs sage.symbolic
        [0, 35/128, 0, -5/16, 0, 3/8, 0, -1/2, 0, 1,
         1, 0, -1/2, 0, 3/8, 0, -5/16, 0, 35/128, 0]

    Verify that :issue:`33963` is fixed::

        sage: # needs sage.symbolic
        sage: n = var("n")
        sage: assume(n, "integer")
        sage: assume(n, "even")
        sage: legendre_P(n, 0)
        2^(-n + 2)*(-1)^(1/2*n)*gamma(n)/(n*gamma(1/2*n)^2)
        sage: forget()
    """
    def __init__(self):
        r"""
        Init method for the Legendre polynomials of the first kind.

        EXAMPLES::

            sage: loads(dumps(legendre_P))
            legendre_P
        """
        BuiltinFunction.__init__(self, 'legendre_P', nargs=2, latex_name=r"P",
                                 conversions={'maxima': 'legendre_p',
                                              'mathematica': 'LegendreP',
                                              'maple': 'LegendreP',
                                              'giac': 'legendre'})


legendre_P = Func_legendre_P()


class Func_legendre_Q(BuiltinFunction):
    def __init__(self):
        r"""
        EXAMPLES::

            sage: loads(dumps(legendre_Q))
            legendre_Q
            sage: maxima(legendre_Q(20, x, hold=True))._sage_().coefficient(x, 10)      # needs sage.symbolic
            -29113619535/131072*log(-(x + 1)/(x - 1))
        """
        BuiltinFunction.__init__(self, "legendre_Q", nargs=2, latex_name=r"Q",
                conversions={'maxima': 'legendre_q',
                             'mathematica': 'LegendreQ',
                             'maple': 'LegendreQ'})

    def _eval_(self, n, x, *args, **kwds):
        r"""
        Return an evaluation of this Legendre Q expression.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: legendre_Q(2,x)
            1/4*(3*x^2 - 1)*(log(x + 1) - log(-x + 1)) - 3/2*x
            sage: legendre_Q(5, 0)
            -8/15
            sage: legendre_Q(2, 2*x)
            1/4*(12*x^2 - 1)*(log(2*x + 1) - log(-2*x + 1)) - 3*x
            sage: legendre_Q(1/2, I+1.)
            -0.511424110789061 + 1.34356195297194*I
            sage: legendre_Q(-1, x)
            Infinity
        """
        ret = self._eval_special_values_(n, x)
        if ret is not None:
            return ret
        if n in ZZ:
            if n < 0:
                from sage.rings.infinity import unsigned_infinity
                return SR(unsigned_infinity)
            return self.eval_formula(n, x)

    def _eval_special_values_(self, n, x):
        """
        Special values known.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: var('n')
            n
            sage: legendre_Q(n, 0)
            -1/2*sqrt(pi)*gamma(1/2*n + 1/2)*sin(1/2*pi*n)/gamma(1/2*n + 1)
            sage: legendre_Q(-1., 0.)
            +infinity
            sage: legendre_Q(-1/2, 2)
            elliptic_kc(3/2)
        """
        if n == QQ(-1)/2:
            from sage.functions.special import elliptic_kc
            return elliptic_kc((x+1)/2)

        if x == 1:
            from sage.rings.infinity import unsigned_infinity
            return SR(unsigned_infinity)

        if x == -1:
            from sage.rings.infinity import unsigned_infinity
            return SR(unsigned_infinity)

        if x == 0:
            from .gamma import gamma
            from .other import sqrt
            from .trig import sin
            try:
                gam = gamma((n+1)/2)/gamma(n/2 + 1)
                if gam.is_infinity():
                    return gam
                return -(sqrt(SR.pi()))/2 * sin(SR.pi()/2*n) * gam
            except TypeError:
                pass

    def _evalf_(self, n, x, parent=None, **kwds):
        """
        Float evaluation of Legendre Q(n, x) function.

        EXAMPLES::

            sage: legendre_Q(4, 2.)                                                     # needs mpmath
            0.00116107583162041 - 86.9828465962674*I
            sage: legendre_Q(1/2, I+1.)                                                 # needs sage.symbolic
            -0.511424110789061 + 1.34356195297194*I
            sage: legendre_Q(1/2, I+1).n(59)                                            # needs sage.symbolic
            -0.51142411078906080 + 1.3435619529719394*I
        """
        ret = self._eval_special_values_(n, x)
        if ret is not None:
            return ret

        return _mpmath_call(_mpmath_legenq, n, 0, x, parent=parent)

    def eval_recursive(self, n, arg, **kwds):
        """
        Return expanded Legendre Q(n, arg) function expression.

        EXAMPLES::

            sage: legendre_Q.eval_recursive(2, x)                                       # needs sage.symbolic
            3/4*x^2*(log(x + 1) - log(-x + 1)) - 3/2*x - 1/4*log(x + 1) + 1/4*log(-x + 1)
            sage: legendre_Q.eval_recursive(20, x).expand().coefficient(x, 10)          # needs sage.symbolic
            -29113619535/131072*log(x + 1) + 29113619535/131072*log(-x + 1)
        """
        from sage.functions.log import ln
        if n == 0:
            return (ln(1+arg)-ln(1-arg))/2
        elif n == 1:
            return arg/2*(ln(1+arg)-ln(1-arg))-1

        x, l = PolynomialRing(QQ, 'x,l').gens()
        help1 = l / 2
        help2 = x / 2 * l - 1
        for j in range(1, n):
            help3 = (2 * j + 1) * x * help2 - j * help1
            help3 = help3 / (j + 1)
            help1 = help2
            help2 = help3

        sum1 = sum(help3.monomial_coefficient(mon)*arg**(mon.exponents()[0][0])
                   for mon in help3.monomials() if not l.divides(mon))
        sum2 = sum(help3.monomial_coefficient(mon)*arg**(mon.exponents()[0][0])*(ln(1+arg)-ln(1-arg))
                   for mon in help3.monomials() if l.divides(mon))
        return sum1 + sum2

    def eval_formula(self, n, arg, **kwds):
        """
        Return expanded Legendre ``Q(n, arg)`` function expression.

        REFERENCE:

        - T. M. Dunster, Legendre and Related Functions, https://dlmf.nist.gov/14.7#E2

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: legendre_Q.eval_formula(1, x)
            1/2*x*(log(x + 1) - log(-x + 1)) - 1
            sage: legendre_Q.eval_formula(2, x).expand().collect(log(1+x)).collect(log(1-x))
            1/4*(3*x^2 - 1)*log(x + 1) - 1/4*(3*x^2 - 1)*log(-x + 1) - 3/2*x
            sage: legendre_Q.eval_formula(20, x).coefficient(x, 10)
            -29113619535/131072*log(x + 1) + 29113619535/131072*log(-x + 1)
            sage: legendre_Q(0, 2)
            -1/2*I*pi + 1/2*log(3)

            sage: legendre_Q(0, 2.)                                                     # needs mpmath
            0.549306144334055 - 1.57079632679490*I
        """
        from sage.functions.log import ln
        if n == 0:
            return (ln(1+arg)-ln(1-arg))/2
        elif n == 1:
            return arg/2*(ln(1+arg)-ln(1-arg))-1

        arg = SR(arg)
        return legendre_P(n, arg)*(ln(1+arg)-ln(1-arg))/2 - self._Wfunc(n, arg)

    def _Wfunc(self, n, arg):
        """
        Helper function for ``eval_formula()``.

        EXAMPLES::

            sage: legendre_Q._Wfunc(2, x)                                               # needs sage.symbolic
            3/2*x
            sage: legendre_Q._Wfunc(7, x)                                               # needs sage.symbolic
            429/16*x^6 - 275/8*x^4 + 849/80*x^2 - 16/35
        """
        if n == 0:
            return 0
        if n == 1:
            return 1
        x = PolynomialRing(QQ, 'x').gen()
        help1 = 0
        help2 = 1
        for j in range(2, n + 1):
            help3 = (2 * j - 1) * x * help2 - (j - 1) * help1
            help3 = help3 / j
            help1 = help2
            help2 = help3

        return sum(b * arg**a for a, b in enumerate(help3))

    def _derivative_(self, n, x, *args, **kwds):
        """
        Return the derivative of legendre_Q.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: n = var('n')
            sage: derivative(legendre_Q(n,x), x)
            (n*x*legendre_Q(n, x) - n*legendre_Q(n - 1, x))/(x^2 - 1)
            sage: ex1 = legendre_Q(5, x, hold=True).diff(x).expand().simplify_full()
            sage: ex2 = legendre_Q(5, x).diff(x).expand().simplify_full()
            sage: ex1.subs(x=7).n() == ex2.subs(x=7).n()
            True
            sage: derivative(legendre_Q(n, x), n)
            Traceback (most recent call last):
            ...
            NotImplementedError: Derivative w.r.t. to the index is not supported.
        """
        diff_param = kwds['diff_param']
        if diff_param == 0:
            raise NotImplementedError("Derivative w.r.t. to the index is not supported.")
        else:
            return (n*x*legendre_Q(n, x) - n*legendre_Q(n-1, x))/(x**2 - 1)


legendre_Q = Func_legendre_Q()


class Func_assoc_legendre_P(BuiltinFunction):
    r"""
    Return the Ferrers function `\mathtt{P}_n^m(x)` of first kind for
    `x \in (-1,1)` with general order `m` and general degree `n`.

    Ferrers functions of first kind are one of two linearly independent
    solutions of the associated Legendre differential equation

    .. MATH::

        (1-x^2) \frac{\mathrm{d}^2 w}{\mathrm{d}x^2} -
            2x \frac{\mathrm{d} w}{\mathrm{d}x} +
            \left(n(n+1) - \frac{m^2}{1-x^2}\right) w = 0

    on the interval `x \in (-1, 1)` and are usually denoted by
    `\mathtt{P}_n^m(x)`.

    .. SEEALSO ::

        The other linearly independent solution is called *Ferrers function of
        second kind* and denoted by `\mathtt{Q}_n^m(x)`,
        see :class:`Func_assoc_legendre_Q`.

    .. WARNING::

        Ferrers functions must be carefully distinguished from associated
        Legendre functions which are defined on `\CC \setminus (- \infty, 1]`
        and have not yet been implemented.

    EXAMPLES:

    We give the first Ferrers functions for nonnegative integers
    `n` and `m` in the interval `-1<x<1`::

        sage: for n in range(4):                                                        # needs sage.symbolic
        ....:     for m in range(n+1):
        ....:         print(f"P_{n}^{m}({x}) = {gen_legendre_P(n, m, x)}")
        P_0^0(x) = 1
        P_1^0(x) = x
        P_1^1(x) = -sqrt(-x^2 + 1)
        P_2^0(x) = 3/2*x^2 - 1/2
        P_2^1(x) = -3*sqrt(-x^2 + 1)*x
        P_2^2(x) = -3*x^2 + 3
        P_3^0(x) = 5/2*x^3 - 3/2*x
        P_3^1(x) = -3/2*(5*x^2 - 1)*sqrt(-x^2 + 1)
        P_3^2(x) = -15*(x^2 - 1)*x
        P_3^3(x) = -15*(-x^2 + 1)^(3/2)

    These expressions for nonnegative integers are computed by the
    Rodrigues-type given in :meth:`eval_gen_poly`. Negative values for `n` are
    obtained by the following identity:

    .. MATH::

        P^{m}_{-n}(x) = P^{m}_{n-1}(x).

    For `n` being a nonnegative integer, negative values for `m` are
    obtained by

    .. MATH::

        P^{-|m|}_n(x) = (-1)^{|m|} \frac{(n-|m|)!}{(n+|m|)!} P_n^{|m|}(x),

    where `|m| \leq n`.

    Here are some specific values with negative integers::

        sage: # needs sage.symbolic
        sage: gen_legendre_P(-2, -1, x)
        1/2*sqrt(-x^2 + 1)
        sage: gen_legendre_P(2, -2, x)
        -1/8*x^2 + 1/8
        sage: gen_legendre_P(3, -2, x)
        -1/8*(x^2 - 1)*x
        sage: gen_legendre_P(1, -2, x)
        0

    Here are some other random values with floating numbers::

        sage: # needs sage.symbolic
        sage: m = var('m'); assume(m, 'integer')
        sage: gen_legendre_P(m, m, .2)
        0.960000000000000^(1/2*m)*(-1)^m*factorial(2*m)/(2^m*factorial(m))
        sage: gen_legendre_P(.2, m, 0)
        sqrt(pi)*2^m/(gamma(-1/2*m + 1.10000000000000)*gamma(-1/2*m + 0.400000000000000))
        sage: gen_legendre_P(.2, .2, .2)
        0.757714892929573

    TESTS:

    Some consistency checks::

        sage: gen_legendre_P(1, 1, x)                                                   # needs sage.symbolic
        -sqrt(-x^2 + 1)
        sage: gen_legendre_P.eval_gen_poly(1, 1, x)                                     # needs sage.symbolic
        -sqrt(-x^2 + 1)
        sage: gen_legendre_P(1, 1, 0.5)  # abs tol 1e-14                                # needs mpmath
        -0.866025403784439
        sage: gen_legendre_P.eval_gen_poly(1, 1, 0.5)  # abs tol 1e-14                  # needs sage.rings.real_mpfr
        -0.866025403784439
        sage: gen_legendre_P._evalf_(1, 1, 0.5)  # abs tol 1e-14                        # needs mpmath
        -0.866025403784439
        sage: gen_legendre_P(2/3, 1, 0.)  # abs tol 1e-14                               # needs mpmath
        -0.773063511309286
        sage: gen_legendre_P._eval_special_values_(2/3, 1, 0.).n()  # abs tol 1e-14     # needs sage.symbolic
        -0.773063511309286

    REFERENCES:

    - [DLMF-Legendre]_
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: loads(dumps(gen_legendre_P))
            gen_legendre_P
            sage: maxima(gen_legendre_P(20, 6, x, hold=True))._sage_().expand().coefficient(x,10)   # needs sage.symbolic
            2508866163428625/128

        TESTS::

            sage: fricas(gen_legendre_P(2, 1/2, x))                             # optional - fricas, needs sage.symbolic
                        1
            legendreP(2,-,x)
                        2

            sage: gen_legendre_P(3, 0, x)                                               # needs sage.symbolic
            5/2*x^3 - 3/2*x
            sage: fricas.legendreP(3, x)                                        # optional - fricas, needs sage.symbolic
            5  3   3
            - x  - - x
            2      2
        """
        BuiltinFunction.__init__(self, "gen_legendre_P", nargs=3,
                                 latex_name=r"\mathtt{P}",
                                 conversions={'maxima': 'assoc_legendre_p',
                                              'mathematica': 'LegendreP',
                                              'fricas': 'legendreP',
                                              'maple': 'LegendreP'})

    def _eval_(self, n, m, x, *args, **kwds):
        r"""
        Return an evaluation of this Legendre P(n, m, x) expression.

        EXAMPLES::

            sage: gen_legendre_P(13/2, 2, 0)                                            # needs sage.symbolic
            4*sqrt(pi)/(gamma(13/4)*gamma(-15/4))
            sage: gen_legendre_P(3, 2, x)                                               # needs sage.symbolic
            -15*(x^2 - 1)*x
        """
        ret = self._eval_special_values_(n, m, x)
        if ret is not None:
            return ret
        if n in ZZ and m in ZZ and (x in ZZ or not SR(x).is_numeric()):
            return self._eval_int_ord_deg_(n, m, x)

    def _eval_special_values_(self, n, m, x):
        """
        Special values known.

        EXAMPLES:

        Case `|m| > |n|` for integers::

            sage: gen_legendre_P(2, 3, 4)                                               # needs mpmath
            0

        Case `x = 0`::

            sage: # needs sage.symbolic
            sage: gen_legendre_P(13/2, 2, 0)
            4*sqrt(pi)/(gamma(13/4)*gamma(-15/4))
            sage: m, n = var('m,n')
            sage: gen_legendre_P(n, m, 0)
            sqrt(pi)*2^m/(gamma(-1/2*m + 1/2*n + 1)*gamma(-1/2*m - 1/2*n + 1/2))
            sage: gen_legendre_P(n, 3, 0)
            8*sqrt(pi)/(gamma(1/2*n - 1/2)*gamma(-1/2*n - 1))
            sage: gen_legendre_P(3, m, 0)
            sqrt(pi)*2^m/(gamma(-1/2*m + 5/2)*gamma(-1/2*m - 1))

        Case `m = n` for integers::

            sage: # needs sage.symbolic
            sage: m = var('m')
            sage: assume(m, 'integer')
            sage: gen_legendre_P(m, m, x)
            (-1)^m*(-x^2 + 1)^(1/2*m)*factorial(2*m)/(2^m*factorial(m))
            sage: gen_legendre_P(m, m, .2)
            0.960000000000000^(1/2*m)*(-1)^m*factorial(2*m)/(2^m*factorial(m))
            sage: gen_legendre_P(2, 2, x)
            -3*x^2 + 3

        Case `n = 0`::

            sage: gen_legendre_P(m, 0, x)                                               # needs sage.symbolic
            legendre_P(m, x)
            sage: gen_legendre_P(2, 0, 4) == legendre_P(2, 4)                           # needs sage.symbolic
            True
        """
        if m == 0:
            # https://dlmf.nist.gov/14.7#E1
            return legendre_P(n, x)
        if x == 0:
            from .gamma import gamma
            from .other import sqrt
            # https://dlmf.nist.gov/14.5#E1
            return 2**m*sqrt(SR.pi())/gamma(n/2-m/2+1)/gamma(QQ(1/2)-n/2-m/2)
        if m.is_integer() and n.is_integer():
            if abs(m) > abs(n):
                # https://dlmf.nist.gov/14.7#E10 and https://dlmf.nist.gov/14.9#E3
                # and https://dlmf.nist.gov/14.9#E5
                return ZZ.zero()
            if m == n:
                # http://dlmf.nist.gov/14.5.iv and https://dlmf.nist.gov/14.9#E3
                return (-1)**m*factorial(2*m)/(2**m*factorial(m)) * (1-x**2)**(m/2)

    def _eval_int_ord_deg_(self, n, m, x):
        r"""
        Evaluate the Ferrers function `P(n, m, x)` for `m` and `n` being
        concrete integers.

        TESTS::

            sage: gen_legendre_P._eval_int_ord_deg_(-2, 1, x)                           # needs sage.symbolic
            -sqrt(-x^2 + 1)
            sage: gen_legendre_P._eval_int_ord_deg_(2, -1, x)                           # needs sage.symbolic
            1/2*sqrt(-x^2 + 1)*x
            sage: gen_legendre_P._eval_int_ord_deg_(-2, -1, x)                          # needs sage.symbolic
            1/2*sqrt(-x^2 + 1)
        """
        # use connection formulas to fall back on nonnegative n and m:
        if n < 0:
            # https://dlmf.nist.gov/14.9#E5
            return self._eval_int_ord_deg_(-n-1, m, x)
        if m < 0:
            # https://dlmf.nist.gov/14.9#E3
            return (-1)**(-m)*factorial(n+m)/factorial(n-m) * self._eval_int_ord_deg_(n, -m, x)
        # apply Rodrigues formula:
        return self.eval_gen_poly(n, m, x)

    def _evalf_(self, n, m, x, parent=None, **kwds):
        """
        Float evaluation of Ferrers function P(n, m, x).

        EXAMPLES::

            sage: gen_legendre_P(10, 2, 3).n()  # abs tol 1e-14                         # needs sage.symbolic
            -7.19496360000000e8
            sage: gen_legendre_P(5/2, 2, 1. + I)                                        # needs sage.symbolic
            14.3165258449040 - 12.7850496155152*I
            sage: gen_legendre_P(5/2, 2, ComplexField(70)(1+I))                         # needs sage.rings.real_mpfr sage.symbolic
            14.316525844904028532 - 12.785049615515157033*I
            sage: gen_legendre_P(2/3, 1, 0.)                                            # needs mpmath
            -0.773063511309286
        """
        return _mpmath_call(_mpmath_legenp, n, m, x, parent=parent)

    def eval_gen_poly(self, n, m, arg, **kwds):
        r"""
        Return the Ferrers function of first kind `\mathtt{P}_n^m(x)` for
        integers `n > -1, m > -1` given by the following Rodrigues-type
        formula:

        .. MATH::

            \mathtt{P}_n^m(x) = (-1)^{m+n} \frac{(1-x^2)^{m/2}}{2^n n!}
                \frac{\mathrm{d}^{m+n}}{\mathrm{d}x^{m+n}} (1-x^2)^n.

        INPUT:

        - ``n`` -- integer degree
        - ``m`` -- integer order
        - ``x`` -- either an integer or a non-numerical symbolic expression

        EXAMPLES::

            sage: gen_legendre_P(7, 4, x)                                               # needs sage.symbolic
            3465/2*(13*x^3 - 3*x)*(x^2 - 1)^2
            sage: gen_legendre_P(3, 1, sqrt(x))                                         # needs sage.symbolic
            -3/2*(5*x - 1)*sqrt(-x + 1)

        REFERENCE:

        - [DLMF-Legendre]_, Section 14.7 eq. 10 (https://dlmf.nist.gov/14.7#E10)
        """
        if n < 0 or m < 0:
            return
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        p = (1-x**2)**ZZ(n)
        for _ in range(m + n):
            p = p.diff(x)
        ex1 = (1-arg**2)**(QQ(m)/2)/2**n/factorial(ZZ(n))
        ex2 = sum(b * arg**a for a, b in enumerate(p))
        return (-1)**(m+n)*ex1*ex2

    from sage.misc.superseded import deprecated_function_alias
    eval_poly = deprecated_function_alias(25034, eval_gen_poly)

    def _derivative_(self, n, m, x, *args, **kwds):
        """
        Return the derivative of ``gen_legendre_P(n,m,x)``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: m, n = var('m,n')
            sage: derivative(gen_legendre_P(n,m,x), x)
            -((n + 1)*x*gen_legendre_P(n, m, x)
             + (m - n - 1)*gen_legendre_P(n + 1, m, x))/(x^2 - 1)
            sage: gen_legendre_P(3, 2, x, hold=True).diff(x).expand().simplify_full()
            -45*x^2 + 15
            sage: derivative(gen_legendre_P(n,m,x), n)
            Traceback (most recent call last):
            ...
            NotImplementedError: Derivative w.r.t. to the index is not supported.
        """
        diff_param = kwds['diff_param']
        if diff_param == 0:
            raise NotImplementedError("Derivative w.r.t. to the index is not supported.")
        else:
            # https://dlmf.nist.gov/14.10#E4
            return ((m-n-1)*gen_legendre_P(n+1, m, x) + (n+1)*x*gen_legendre_P(n, m, x))/(1 - x**2)


gen_legendre_P = Func_assoc_legendre_P()


class Func_assoc_legendre_Q(BuiltinFunction):
    def __init__(self):
        r"""
        EXAMPLES::

            sage: loads(dumps(gen_legendre_Q))
            gen_legendre_Q
            sage: maxima(gen_legendre_Q(2, 1, 3, hold=True))._sage_().simplify_full()   # needs sage.symbolic
            1/4*sqrt(2)*(36*pi - 36*I*log(2) + 25*I)
        """
        BuiltinFunction.__init__(self, "gen_legendre_Q", nargs=3, latex_name=r"Q",
                                 conversions={'maxima': 'assoc_legendre_q',
                                              'mathematica': 'LegendreQ',
                                              'maple': 'LegendreQ'})

    def _eval_(self, n, m, x, *args, **kwds):
        r"""
        Return an evaluation of this Legendre Q(n, m, x) expression.

        EXAMPLES::

            sage: gen_legendre_Q(2, 1, 3)                                               # needs sage.symbolic
            -1/4*sqrt(-2)*(-36*I*pi + 36*log(2) - 25)
        """
        ret = self._eval_special_values_(n, m, x)
        if ret is not None:
            return ret
        if (n in ZZ and m in ZZ
                and n >= 0 and m >= 0
                and (x in ZZ or not SR(x).is_numeric())):
            return self.eval_recursive(n, m, x)

    def _eval_special_values_(self, n, m, x):
        """
        Special values known.

        EXAMPLES::

            sage: n, m = var('n m')                                                     # needs sage.symbolic
            sage: gen_legendre_Q(n, m, 0)                                               # needs sage.symbolic
            -sqrt(pi)*2^(m - 1)*gamma(1/2*m + 1/2*n + 1/2)*sin(1/2*pi*m + 1/2*pi*n)/gamma(-1/2*m + 1/2*n + 1)
        """
        if m == 0:
            return legendre_Q(n, x)
        if x.is_zero():
            from .gamma import gamma
            from .other import sqrt
            from .trig import sin
            if m in QQ and n in QQ:
                return -(sqrt(SR.pi()))*sin(SR.pi()/2*(m+n))*gamma(QQ(m+n+1)/2)/gamma(QQ(n-m)/2 + 1)*2**(m-1)
            elif isinstance(n, Expression) or isinstance(m, Expression):
                return -(sqrt(SR.pi()))*sin(SR.pi()/2*(m+n))*gamma((m+n+1)/2)/gamma((n-m)/2 + 1)*2**(m-1)

    def _evalf_(self, n, m, x, parent=None, **kwds):
        """
        Float evaluation of Legendre Q(n, m, x) function.

        EXAMPLES::

            sage: gen_legendre_Q(2, 1, 3.)                                              # needs mpmath
            -39.9859464434253 + 0.0165114736149193*I
            sage: gen_legendre_Q(2, 1, ComplexField(70)(3))                             # needs sage.rings.real_mpfr
            -39.985946443425296223 + 0.016511473614919329585*I
        """
        ret = self._eval_special_values_(n, m, x)
        if ret is not None:
            return ret

        return _mpmath_call(_mpmath_legenq, n, m, x, parent=parent)

    def eval_recursive(self, n, m, x, **kwds):
        """
        Return the associated Legendre Q(n, m, arg) function for integers `n > -1, m > -1`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: gen_legendre_Q(3, 4, x)
            48/(x^2 - 1)^2
            sage: gen_legendre_Q(4, 5, x)
            -384/((x^2 - 1)^2*sqrt(-x^2 + 1))
            sage: gen_legendre_Q(0, 1, x)
            -1/sqrt(-x^2 + 1)
            sage: gen_legendre_Q(0, 2, x)
            -1/2*((x + 1)^2 - (x - 1)^2)/(x^2 - 1)
            sage: gen_legendre_Q(2, 2, x).subs(x=2).expand()
            9/2*I*pi - 9/2*log(3) + 14/3
        """
        from sage.misc.functional import sqrt
        if m == n + 1 or n == 0:
            if m.mod(2).is_zero():
                denom = (1 - x**2)**(m/2)
            else:
                denom = sqrt(1 - x**2)*(1 - x**2)**((m-1)/2)
            if m == n + 1:
                return (-1)**m*(m-1).factorial()*2**n/denom
            else:
                return (-1)**m*(m-1).factorial()*((x+1)**m - (x-1)**m)/(2*denom)
        else:
            return ((n-m+1)*x*gen_legendre_Q(n, m-1, x)-(n+m-1)*gen_legendre_Q(n-1, m-1, x))/sqrt(1-x**2)

    def _derivative_(self, n, m, x, *args, **kwds):
        """
        Return the derivative of ``gen_legendre_Q(n,m,x)``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: m, n = var('m,n')
            sage: derivative(gen_legendre_Q(n,m,x), x)
            -((n + 1)*x*gen_legendre_Q(n, m, x)
             + (m - n - 1)*gen_legendre_Q(n + 1, m, x))/(x^2 - 1)
            sage: ex1 = gen_legendre_Q(3, 2, x, hold=True).diff(x).expand().simplify_full()
            sage: ex2 = gen_legendre_Q(3, 2, x).diff(x).expand().simplify_full()
            sage: ex1.subs(x=5).n() == ex2.subs(x=5).n()
            True
            sage: derivative(gen_legendre_Q(n,m,x), n)
            Traceback (most recent call last):
            ...
            NotImplementedError: Derivative w.r.t. to the index is not supported.
        """
        diff_param = kwds['diff_param']
        if diff_param == 0:
            raise NotImplementedError("Derivative w.r.t. to the index is not supported.")
        else:
            return ((n-m+1)*gen_legendre_Q(n+1, m, x) - (n+1)*x*gen_legendre_Q(n, m, x))/(x**2 - 1)


gen_legendre_Q = Func_assoc_legendre_Q()


class Func_hermite(GinacFunction):
    r"""
    Return the Hermite polynomial for integers `n > -1`.

    REFERENCE:

    - [AS1964]_ 22.5.40 and 22.5.41, page 779.

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: hermite(2, x)
        4*x^2 - 2
        sage: hermite(3, x)
        8*x^3 - 12*x
        sage: hermite(3, 2)
        40
        sage: S.<y> = PolynomialRing(RR)
        sage: hermite(3, y)
        8.00000000000000*y^3 - 12.0000000000000*y
        sage: R.<x,y> = QQ[]
        sage: hermite(3, y^2)
        8*y^6 - 12*y^2
        sage: w = var('w')
        sage: hermite(3, 2*w)
        64*w^3 - 24*w
        sage: hermite(5, 3.1416)
        5208.69733891963
        sage: hermite(5, RealField(100)(pi))
        5208.6167627118104649470287166

    Check that :issue:`17192` is fixed::

        sage: # needs sage.symbolic
        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: hermite(0, x)
        1
        sage: hermite(-1, x)
        Traceback (most recent call last):
        ...
        RuntimeError: hermite_eval: The index n must be a nonnegative integer
        sage: hermite(-7, x)
        Traceback (most recent call last):
        ...
        RuntimeError: hermite_eval: The index n must be a nonnegative integer
        sage: m, x = SR.var('m,x')
        sage: hermite(m, x).diff(m)
        Traceback (most recent call last):
        ...
        RuntimeError: derivative w.r.t. to the index is not supported yet
    """
    def __init__(self):
        r"""
        Init method for the Hermite polynomials.

        EXAMPLES::

            sage: loads(dumps(hermite))
            hermite
            sage: hermite(x, x)._sympy_()                                               # needs sympy sage.symbolic
            hermite(x, x)

        TESTS::

            sage: fricas(hermite(x, 5))         # optional - fricas                     # needs sage.symbolic
            hermiteH(x,5)

            sage: hermite(5, x)                                                         # needs sage.symbolic
            32*x^5 - 160*x^3 + 120*x
            sage: fricas.hermiteH(5, x)         # optional - fricas                     # needs sage.symbolic
                5        3
            32 x  - 160 x  + 120 x
        """
        GinacFunction.__init__(self, "hermite", nargs=2, latex_name=r"H",
                conversions={'maxima': 'hermite',
                             'mathematica': 'HermiteH',
                             'maple': 'HermiteH',
                             'fricas': 'hermiteH',
                             'sympy': 'hermite'}, preserved_arg=2)


hermite = Func_hermite()


class Func_jacobi_P(OrthogonalFunction):
    r"""
    Return the Jacobi polynomial `P_n^{(a,b)}(x)` for integers
    `n > -1` and a and b symbolic or `a > -1` and `b > -1`.

    The Jacobi polynomials are actually defined for all `a` and `b`.
    However, the Jacobi polynomial weight `(1-x)^a(1+x)^b` is not
    integrable for `a \leq -1` or `b \leq -1`.

    REFERENCE:

    - Table on page 789 in [AS1964]_.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: jacobi_P(2, 0, 0, x)                                                      # needs sage.libs.flint sage.symbolic
        3/2*x^2 - 1/2
        sage: jacobi_P(2, 1, 2, 1.2)                                                    # needs sage.libs.flint
        5.01000000000000
    """
    def __init__(self):
        r"""
        Init method for the Jacobi polynomials.

        EXAMPLES::

            sage: n, a, b, x = SR.var('n,a,b,x')                                        # needs sage.symbolic
            sage: loads(dumps(jacobi_P))
            jacobi_P
            sage: jacobi_P(n, a, b, x, hold=True)._sympy_()                             # needs sympy sage.symbolic
            jacobi(n, a, b, x)

        TESTS::

            sage: fricas(jacobi_P(1/2, 4, 1/3, x))                              # optional - fricas, needs sage.symbolic
                    1   1
            jacobiP(-,4,-,x)
                    2   3

            sage: jacobi_P(1, 2, 3, x)                                                  # needs sage.symbolic
            7/2*x - 1/2
            sage: fricas.jacobiP(1, 2, 3, x)                                    # optional - fricas, needs sage.symbolic
            7 x - 1
            -------
               2
        """
        OrthogonalFunction.__init__(self, "jacobi_P", nargs=4, latex_name=r"P",
                conversions={'maxima': 'jacobi_p',
                             'mathematica': 'JacobiP',
                             'maple': 'JacobiP',
                             'fricas': 'jacobiP',
                             'sympy': 'jacobi'})

    def _eval_(self, n, a, b, x):
        """
        EXAMPLES::

            sage: # needs sage.symbolic
            sage: n, a, b, x = SR.var('n,a,b,x')
            sage: jacobi_P(1, n, n, n)
            (n + 1)*n
            sage: jacobi_P(2, n, n, n)
            1/4*(2*n - 1)*(n + 2)*(n + 1)^2
            sage: jacobi_P(1, n, n, x)
            (n + 1)*x
            sage: jacobi_P(3, 2, 1, x)
            21/2*x^3 + 7/2*x^2 - 7/2*x - 1/2
            sage: jacobi_P(1, a, b, x)
            1/2*a*x + 1/2*b*x + 1/2*a - 1/2*b + x

        TESTS:

        Check that :issue:`17192` is fixed::

            sage: x = PolynomialRing(QQ, 'x').gen()
            sage: jacobi_P(0, 0, 0, x)                                                  # needs sage.libs.flint sage.symbolic
            1
            sage: jacobi_P(-1, 0, 0, x)                                                 # needs sage.libs.flint sage.symbolic
            1
            sage: jacobi_P(-1, 1, 1, x)                                                 # needs sage.libs.flint sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: n must be greater than -1, got n = -1

            sage: jacobi_P(-7, 0, 0, x)                                                 # needs sage.libs.flint sage.symbolic
            231/16*x^6 - 315/16*x^4 + 105/16*x^2 - 5/16
            sage: jacobi_P(-7, 0, 2, x)                                                 # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: n must be greater than -1, got n = -7
        """
        if SR(a).is_trivial_zero() and SR(b).is_trivial_zero():
            return legendre_P(n, x)
        if SR(n).is_numeric() and not (n > -1):
            raise ValueError("n must be greater than -1, got n = {0}".format(n))
        if n not in ZZ:
            return
        from .gamma import gamma
        s = sum(binomial(n, m) * gamma(a+b+n+m+1) / gamma(a+m+1) * ((x-1)/2)**m for m in range(n+1))
        r = gamma(a+n+1) / factorial(n) / gamma(n+a+b+1) * s
        return r.to_gamma().gamma_normalize().normalize()

    def _evalf_(self, n, a, b, x, **kwds):
        """
        EXAMPLES::

            sage: jacobi_P(2, 1, 2, 1.2)                                                # needs sage.symbolic
            5.01000000000000
            sage: jacobi_P(2, 1, 2, 1.2, hold=True).n(20)                               # needs sage.symbolic
            5.0100
            sage: jacobi_P(2, 1, 2, pi + I, hold=True).n(100)                           # needs sage.symbolic
            41.103034125334442891187112674 + 31.486722862692829003857755524*I
        """
        from sage.rings.complex_arb import ComplexBallField as CBF
        the_parent = kwds.get('parent', None)
        if the_parent is None:
            the_parent = parent(x)
        prec = the_parent.precision()
        BF = CBF(prec+5)
        ret = BF(x).jacobi_P(BF(n), BF(a), BF(b))
        return SR(ret)._eval_self(the_parent)


jacobi_P = Func_jacobi_P()


class Func_ultraspherical(GinacFunction):
    r"""
    Return the ultraspherical (or Gegenbauer) polynomial ``gegenbauer(n,a,x)``,

    .. MATH::

        C_n^{a}(x) = \sum_{k=0}^{\lfloor n/2\rfloor} (-1)^k
        \frac{\Gamma(n-k+a)}{\Gamma(a)k!(n-2k)!} (2x)^{n-2k}.

    When `n` is a nonnegative integer, this formula gives a
    polynomial in `z` of degree `n`, but all parameters are
    permitted to be complex numbers. When `a = 1/2`, the
    Gegenbauer polynomial reduces to a Legendre polynomial.

    Computed using Pynac.

    For numerical evaluation, consider using the `mpmath library
    <http://mpmath.org/doc/current/functions/orthogonal.html#gegenbauer-polynomials>`_,
    as it also allows complex numbers (and negative `n` as well);
    see the examples below.

    REFERENCE:

    - [AS1964]_ 22.5.27

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: ultraspherical(8, 101/11, x)
        795972057547264/214358881*x^8 - 62604543852032/19487171*x^6...
        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: ultraspherical(2, 3/2, x)
        15/2*x^2 - 3/2
        sage: ultraspherical(1, 1, x)
        2*x
        sage: t = PolynomialRing(RationalField(), "t").gen()
        sage: gegenbauer(3, 2, t)
        32*t^3 - 12*t
        sage: x = SR.var('x')
        sage: n = ZZ.random_element(5, 5001)
        sage: a = QQ.random_element().abs() + 5
        sage: s = (  (n + 1)*ultraspherical(n + 1, a, x)
        ....:      - 2*x*(n + a)*ultraspherical(n, a, x)
        ....:      + (n + 2*a - 1)*ultraspherical(n - 1, a, x) )
        sage: s.expand().is_zero()
        True
        sage: ultraspherical(5, 9/10, 3.1416)
        6949.55439044240
        sage: ultraspherical(5, 9/10, RealField(100)(pi))                               # needs sage.rings.real_mpfr
        6949.4695419382702451843080687

        sage: # needs sage.symbolic
        sage: a, n = SR.var('a,n')
        sage: gegenbauer(2, a, x)
        2*(a + 1)*a*x^2 - a
        sage: gegenbauer(3, a, x)
        4/3*(a + 2)*(a + 1)*a*x^3 - 2*(a + 1)*a*x
        sage: gegenbauer(3, a, x).expand()
        4/3*a^3*x^3 + 4*a^2*x^3 + 8/3*a*x^3 - 2*a^2*x - 2*a*x
        sage: gegenbauer(10, a, x).expand().coefficient(x, 2)
        1/12*a^6 + 5/4*a^5 + 85/12*a^4 + 75/4*a^3 + 137/6*a^2 + 10*a
        sage: ex = gegenbauer(100, a, x)
        sage: (ex.subs(a==55/98) - gegenbauer(100, 55/98, x)).is_trivial_zero()
        True

        sage: # needs sage.symbolic
        sage: gegenbauer(2, -3, x)
        12*x^2 + 3
        sage: gegenbauer(120, -99/2, 3)
        1654502372608570682112687530178328494861923493372493824
        sage: gegenbauer(5, 9/2, x)
        21879/8*x^5 - 6435/4*x^3 + 1287/8*x
        sage: gegenbauer(15, 3/2, 5)
        3903412392243800

        sage: derivative(gegenbauer(n, a, x), x)                                        # needs sage.symbolic
        2*a*gegenbauer(n - 1, a + 1, x)
        sage: derivative(gegenbauer(3, a, x), x)                                        # needs sage.symbolic
        4*(a + 2)*(a + 1)*a*x^2 - 2*(a + 1)*a
        sage: derivative(gegenbauer(n, a, x), a)                                        # needs sage.symbolic
        Traceback (most recent call last):
        ...
        RuntimeError: derivative w.r.t. to the second index is not supported yet

    Numerical evaluation with the mpmath library::

        sage: # needs mpmath
        sage: from mpmath import gegenbauer as gegenbauer_mp
        sage: from mpmath import mp
        sage: print(gegenbauer_mp(-7,0.5,0.3))
        0.1291811875
        sage: with mp.workdps(25):
        ....:     print(gegenbauer_mp(2+3j, -0.75, -1000j))
        (-5038991.358609026523401901 + 9414549.285447104177860806j)

    TESTS:

    Check that :issue:`17192` is fixed::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: ultraspherical(0, 1, x)                                                   # needs sage.symbolic
        1

        sage: ultraspherical(-1, 1, x)                                                  # needs sage.symbolic
        Traceback (most recent call last):
        ...
        RuntimeError: gegenb_eval: The index n must be a nonnegative integer

        sage: ultraspherical(-7, 1, x)                                                  # needs sage.symbolic
        Traceback (most recent call last):
        ...
        RuntimeError: gegenb_eval: The index n must be a nonnegative integer
    """
    def __init__(self):
        r"""
        Init method for the ultraspherical polynomials.

        EXAMPLES::

            sage: loads(dumps(ultraspherical))
            gegenbauer
            sage: ultraspherical(x, x, x)._sympy_()                                     # needs sympy sage.symbolic
            gegenbauer(x, x, x)
        """
        GinacFunction.__init__(self, "gegenbauer", nargs=3, latex_name=r"C",
                               conversions={'maxima': 'ultraspherical',
                                            'mathematica': 'GegenbauerC',
                                            'maple': 'GegenbauerC',
                                            'sympy': 'gegenbauer'})


ultraspherical = Func_ultraspherical()
gegenbauer = Func_ultraspherical()


class Func_laguerre(OrthogonalFunction):
    """
    REFERENCE:

    - [AS1964]_ 22.5.16, page 778 and page 789.
    """
    def __init__(self):
        r"""
        Init method for the Laguerre polynomials.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: n, x = var('n,x')
            sage: laguerre(x, x)._sympy_()                                              # needs sympy
            laguerre(x, x)
            sage: maxima(laguerre(1, x, hold=True))
            1-_SAGE_VAR_x
            sage: maxima(laguerre(n, laguerre(n, x)))
            laguerre(_SAGE_VAR_n,laguerre(_SAGE_VAR_n,_SAGE_VAR_x))

        TESTS::

            sage: loads(dumps(laguerre))
            laguerre
        """
        OrthogonalFunction.__init__(self, "laguerre", nargs=2, latex_name=r"L",
                conversions={'maxima': 'laguerre',
                             'mathematica': 'LaguerreL',
                             # 'fricas': 'laguerreL',  3 arguments ?
                             'maple': 'LaguerreL',
                             'sympy': 'laguerre'})

    def _eval_(self, n, x, *args, **kwds):
        r"""
        Return an evaluation of this Laguerre polynomial expression.

        EXAMPLES::

            sage: x = PolynomialRing(QQ, 'x').gen()
            sage: laguerre(2, x)                                                        # needs mpmath
            1/2*x^2 - 2*x + 1
            sage: laguerre(3, x)                                                        # needs mpmath
            -1/6*x^3 + 3/2*x^2 - 3*x + 1
            sage: laguerre(2, 2)                                                        # needs mpmath
            -1
            sage: laguerre(-1, x)                                                       # needs sage.symbolic
            e^x
            sage: laguerre(-6, x)                                                       # needs sage.symbolic
            1/120*(x^5 + 25*x^4 + 200*x^3 + 600*x^2 + 600*x + 120)*e^x
            sage: laguerre(-9,2)                                                        # needs sage.symbolic
            66769/315*e^2
        """
        from sage.rings.integer import Integer
        from sage.functions.log import exp
        ret = self._eval_special_values_(n, x)
        if ret is not None:
            return ret
        if isinstance(n, (Integer, int)):
            if n >= 0 and not hasattr(x, 'prec'):
                return self._pol_laguerre(n, x)
            elif n < 0:
                return exp(x)*laguerre(-n-1, -x)

    def _eval_special_values_(self, n, x):
        """
        Special values known.

        EXAMPLES::

            sage: laguerre(0, 0)                                                        # needs mpmath
            1
            sage: laguerre(1, x)                                                        # needs sage.symbolic
            -x + 1
        """
        if n == 0 or x == 0:
            return ZZ(1)
        if n == 1:
            return ZZ(1) - x

    def _pol_laguerre(self, n, x):
        """
        Fast creation of Laguerre polynomial.

        EXAMPLES::

            sage: laguerre(3, sin(x))                                                   # needs sage.symbolic
            -1/6*sin(x)^3 + 3/2*sin(x)^2 - 3*sin(x) + 1
            sage: R.<x> = PolynomialRing(QQ, 'x')
            sage: laguerre(4, x)                                                        # needs mpmath
            1/24*x^4 - 2/3*x^3 + 3*x^2 - 4*x + 1
            sage: laguerre(4, x + 1)                                                    # needs mpmath
            1/24*(x + 1)^4 - 2/3*(x + 1)^3 + 3*(x + 1)^2 - 4*x - 3
            sage: laguerre(10, 1 + I)                                                   # needs sage.symbolic
            142511/113400*I + 95867/22680
        """
        if hasattr(x, 'pyobject'):
            try:
                x = x.pyobject()
            except TypeError:
                pass
        return SR(sum(binomial(n, k) * (-1)**k / factorial(k) * x**k
                      for k in range(n + 1)))

    def _evalf_(self, n, x, **kwds):
        """
        Return the evaluation of `laguerre(n,x)` with floating point `x`.

        EXAMPLES::

            sage: laguerre(100, RealField(300)(pi))                                     # needs sage.symbolic
            -0.638322077840648311606324...
            sage: laguerre(10, 1. + I)                                                  # needs sage.symbolic
            4.22694003527337 + 1.25671075837743*I
            sage: laguerre(-9, 2.)                                                      # needs sage.symbolic
            1566.22186244286
        """
        the_parent = kwds.get('parent', None)
        if the_parent is None:
            the_parent = parent(x)
        if n < 0:
            # work around mpmath issue 307
            from sage.functions.log import exp
            return exp(x) * _mpmath_call(_mpmath_laguerre, -n-1, 0, -x, parent=the_parent)
        else:
            return _mpmath_call(_mpmath_laguerre, n, 0, x, parent=the_parent)

    def _derivative_(self, n, x, *args, **kwds):
        """
        Return the derivative of `laguerre(n,x)`.

        EXAMPLES::

            sage: n = var('n')                                                          # needs sage.symbolic
            sage: diff(laguerre(n, x), x)                                               # needs sage.symbolic
            -gen_laguerre(n - 1, 1, x)

        TESTS::

            sage: diff(laguerre(x, x))                                                  # needs sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: Derivative w.r.t. to the index is not supported.
        """
        diff_param = kwds['diff_param']
        if diff_param == 0:
            raise NotImplementedError("Derivative w.r.t. to the index is not supported.")
        if diff_param == 1:
            return -gen_laguerre(n-1, 1, x)
        raise ValueError(f"illegal differentiation parameter {diff_param}")


laguerre = Func_laguerre()


class Func_gen_laguerre(OrthogonalFunction):
    """
    REFERENCE:

    - [AS1964]_ 22.5.16, page 778 and page 789.
    """
    def __init__(self):
        r"""
        Init method for the Laguerre polynomials.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: a, n, x = var('a, n, x')
            sage: gen_laguerre(x, x, x)._sympy_()                                       # needs sympy
            assoc_laguerre(x, x, x)
            sage: maxima(gen_laguerre(1, 2, x, hold=True))
            3*(1-_SAGE_VAR_x/3)
            sage: maxima(gen_laguerre(n, a, gen_laguerre(n, a, x)))
            gen_laguerre(_SAGE_VAR_n,_SAGE_VAR_a, gen_laguerre(_SAGE_VAR_n,_SAGE_VAR_a,_SAGE_VAR_x))

        TESTS::

            sage: loads(dumps(gen_laguerre))
            gen_laguerre
        """
        OrthogonalFunction.__init__(self, "gen_laguerre", nargs=3, latex_name=r"L",
                                    conversions={'maxima': 'gen_laguerre',
                                                 'mathematica': 'LaguerreL',
                                                 'maple': 'LaguerreL',
                                                 'sympy': 'assoc_laguerre'})

    def _eval_(self, n, a, x, *args, **kwds):
        r"""
        Return an evaluation of this Laguerre polynomial expression.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: gen_laguerre(2, 1, x)
            1/2*x^2 - 3*x + 3
            sage: gen_laguerre(2, 1/2, x)
            1/2*x^2 - 5/2*x + 15/8
            sage: gen_laguerre(2, -1/2, x)
            1/2*x^2 - 3/2*x + 3/8
            sage: gen_laguerre(2, 0, x)
            1/2*x^2 - 2*x + 1
            sage: gen_laguerre(3, 0, x)
            -1/6*x^3 + 3/2*x^2 - 3*x + 1
        """
        from sage.rings.integer import Integer
        ret = self._eval_special_values_(n, a, x)
        if ret is not None:
            return ret
        if isinstance(n, Integer):
            if n >= 0 and not hasattr(x, 'prec'):
                return self._pol_gen_laguerre(n, a, x)

    def _eval_special_values_(self, n, a, x):
        """
        Special values known.

        EXAMPLES::

            sage: gen_laguerre(0, 1, pi)                                                # needs sage.symbolic
            1
            sage: gen_laguerre(1, 2, x)                                                 # needs sage.symbolic
            -x + 3
            sage: gen_laguerre(3, 4, 0)                                                 # needs mpmath
            35
        """
        if n == 0:
            return ZZ(1)
        if n == 1:
            return ZZ(1) + a - x
        if a == 0:
            return laguerre(n, x)
        if x == 0:
            from sage.arith.misc import binomial
            return binomial(n+a, n)

    def _pol_gen_laguerre(self, n, a, x):
        """
        EXAMPLES::

            sage: gen_laguerre(3, 1/2, sin(x))                                          # needs sage.symbolic
            -1/6*sin(x)^3 + 7/4*sin(x)^2 - 35/8*sin(x) + 35/16
            sage: R.<x> = PolynomialRing(QQ, 'x')
            sage: gen_laguerre(4, -1/2, x)                                              # needs mpmath
            1/24*x^4 - 7/12*x^3 + 35/16*x^2 - 35/16*x + 35/128
            sage: gen_laguerre(4, -1/2, x + 1)                                          # needs mpmath
            1/24*(x + 1)^4 - 7/12*(x + 1)^3 + 35/16*(x + 1)^2 - 35/16*x - 245/128
            sage: gen_laguerre(10, 1, 1 + I)                                            # needs sage.symbolic
            25189/2100*I + 11792/2835
        """
        return sum(binomial(n + a, n - k) * (-1)**k / factorial(k) * x**k
                   for k in range(n + 1))

    def _evalf_(self, n, a, x, **kwds):
        """
        EXAMPLES::

            sage: gen_laguerre(100, 1, RealField(300)(pi))                              # needs sage.symbolic
            -0.89430788373354541911...
            sage: gen_laguerre(10, 1/2, 1. + I)                                         # needs sage.symbolic
            5.34469635574906 + 5.23754057922902*I
        """
        the_parent = kwds.get('parent', None)
        if the_parent is None:
            the_parent = parent(x)
        return _mpmath_call(_mpmath_laguerre, n, a, x, parent=the_parent)

    def _derivative_(self, n, a, x, diff_param):
        """
        Return the derivative of `gen_laguerre(n,a,x)`.

        EXAMPLES::

            sage: a, n = var('a,n')                                                     # needs sage.symbolic
            sage: diff(gen_laguerre(n,a,x), x)                                          # needs sage.symbolic
            -gen_laguerre(n - 1, a + 1, x)
            sage: gen_laguerre(n,a,x).diff(a)                                           # needs sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: Derivative w.r.t. to the second index is not supported.

        TESTS::

            sage: diff(gen_laguerre(n,a,x), n)                                          # needs sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: Derivative w.r.t. to the index is not supported.
        """
        if diff_param == 0:
            raise NotImplementedError("Derivative w.r.t. to the index is not supported.")
        elif diff_param == 1:
            raise NotImplementedError("Derivative w.r.t. to the second index is not supported.")
        elif diff_param == 2:
            return -gen_laguerre(n - 1, a + 1, x)
        else:
            raise ValueError("illegal differentiation parameter {}".format(diff_param))


gen_laguerre = Func_gen_laguerre()


class Func_krawtchouk(OrthogonalFunction):
    r"""
    Krawtchouk polynomials `K_j(x; n, p)`.

    INPUT:

    - ``j`` -- the degree
    - ``x`` -- the independent variable `x`
    - ``n`` -- the number of discrete points
    - ``p`` -- the parameter `p`

    .. SEEALSO::

        :func:`sage.coding.delsarte_bounds.krawtchouk`
        `\bar{K}^{n,q}_l(x)`, which are related by

        .. MATH::

            (-q)^j \bar{K}^{n,q^{-1}}_j(x) = K_j(x; n, 1-q).

    EXAMPLES:

    We verify the orthogonality for `n = 4`::

        sage: n = 4
        sage: p = SR.var('p')                                                           # needs sage.symbolic
        sage: matrix([[sum(binomial(n,m) * p**m * (1-p)**(n-m)                          # needs sage.symbolic
        ....:              * krawtchouk(i,m,n,p) * krawtchouk(j,m,n,p)
        ....:              for m in range(n+1)).expand().factor()
        ....:          for i in range(n+1)] for j in range(n+1)])
        [               1                0                0                0                0]
        [               0     -4*(p - 1)*p                0                0                0]
        [               0                0  6*(p - 1)^2*p^2                0                0]
        [               0                0                0 -4*(p - 1)^3*p^3                0]
        [               0                0                0                0    (p - 1)^4*p^4]

    We verify the relationship between the Krawtchouk implementations::

        sage: q = SR.var('q')                                                           # needs sage.symbolic
        sage: all(codes.bounds.krawtchouk(n, 1/q, j, x)*(-q)^j                          # needs sage.symbolic
        ....:     == krawtchouk(j, x, n, 1-q) for j in range(n+1))
        True
    """
    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: k, x, n, p = var('k,x,n,p')                                           # needs sage.symbolic
            sage: TestSuite(krawtchouk).run()
            sage: TestSuite(krawtchouk(k, x, n, p)).run()                               # needs sage.symbolic
            sage: TestSuite(krawtchouk(3, x, n, p)).run()                               # needs sage.symbolic
        """
        super().__init__(name='krawtchouk', nargs=4, latex_name='K')

    def eval_formula(self, k, x, n, p):
        r"""
        Evaluate ``self`` using an explicit formula.

        EXAMPLES::

            sage: x, n, p = var('x,n,p')                                                # needs sage.symbolic
            sage: krawtchouk.eval_formula(3, x, n, p).expand().collect(x)               # needs sage.symbolic
            -1/6*n^3*p^3 + 1/2*n^2*p^3 - 1/3*n*p^3 - 1/2*(n*p - 2*p + 1)*x^2
             + 1/6*x^3 + 1/6*(3*n^2*p^2 - 9*n*p^2 + 3*n*p + 6*p^2 - 6*p + 2)*x
        """
        q = 1 - p
        return sum((-1)**(k-i) * binomial(n-x, k-i) * binomial(x, i) * p**(k-i) * q**i
                   for i in range(k+1))

    def _eval_(self, j, x, n, p, *args, **kwds):
        r"""
        Return an evaluation of the Krawtchouk polynomial `K_j(x; n, p)`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: k, x, n, p = var('k,x,n,p')
            sage: krawtchouk(3, x, 5, p).expand()
            -10*p^3 + 6*p^2*x - 3/2*p*x^2 + 1/6*x^3 + 3/2*p*x - 1/2*x^2 + 1/3*x
            sage: krawtchouk(k, x, 5, p)
            (-1)^k*p^k*binomial(5, k)*hypergeometric((-k, -x), (-5,), 1/p)
            sage: krawtchouk(2, x, n, p).collect(x)
            1/2*n^2*p^2 - 1/2*n*p^2 - 1/2*(2*n*p - 2*p + 1)*x + 1/2*x^2
            sage: krawtchouk(k, x, n, p)
            (-1)^k*p^k*binomial(n, k)*hypergeometric((-k, -x), (-n,), 1/p)

            sage: k3_hypergeo = krawtchouk(k,x,n,p)(k=3).simplify_hypergeometric()      # needs sage.symbolic
            sage: bool(k3_hypergeo == krawtchouk(3,x,n,p))                              # needs sage.symbolic
            True

            sage: krawtchouk(2, x, n, p, hold=True)                                     # needs sage.symbolic
            krawtchouk(2, x, n, p)
        """
        if kwds.get('hold', False):
            return None
        if j not in ZZ or j < 0:
            from sage.functions.hypergeometric import hypergeometric
            return (-1)**j * binomial(n, j) * p**j * hypergeometric([-j, -x], [-n], 1/p)
        try:
            return self.eval_formula(j, x, n, p)
        except (TypeError, ValueError):
            return self.eval_recursive(j, x, n, p)

    def eval_recursive(self, j, x, n, p, *args, **kwds):
        r"""
        Return the Krawtchouk polynomial `K_j(x; n, p)` using the
        recursive formula.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: x, n, p = var('x,n,p')
            sage: krawtchouk.eval_recursive(0, x, n, p)
            1
            sage: krawtchouk.eval_recursive(1, x, n, p)
            -n*p + x
            sage: krawtchouk.eval_recursive(2, x, n, p).collect(x)
            1/2*n^2*p^2 + 1/2*n*(p - 1)*p - n*p^2 + 1/2*n*p
             - 1/2*(2*n*p - 2*p + 1)*x + 1/2*x^2
            sage: bool(krawtchouk.eval_recursive(2, x, n, p) == krawtchouk(2, x, n, p))
            True
            sage: bool(krawtchouk.eval_recursive(3, x, n, p) == krawtchouk(3, x, n, p))
            True
            sage: bool(krawtchouk.eval_recursive(4, x, n, p) == krawtchouk(4, x, n, p))
            True

            sage: M = matrix([[-1/2, -1], [1, 0]])                                      # needs sage.modules
            sage: krawtchouk.eval_recursive(2, M, 3, 1/2)                               # needs sage.modules
            [ 9/8  7/4]
            [-7/4  1/4]
        """
        if j == 0:
            return parent(x).one()
        elif j == 1:
            return x - n * p
        q = 1 - p
        tm2 = p * q * (n - (j-1) + 1) * krawtchouk.eval_recursive(j-2, x, n, p)
        tm1 = (x - p*(n-(j-1)) - (j-1)*q) * krawtchouk.eval_recursive(j-1, x, n, p)
        return (tm1 - tm2) / j


krawtchouk = Func_krawtchouk()


class Func_meixner(OrthogonalFunction):
    r"""
    Meixner polynomials `M_n(x; b, c)`.

    INPUT:

    - ``n`` -- the degree
    - ``x`` -- the independent variable `x`
    - ``b``, ``c`` -- the parameters `b`, `c`
    """
    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: n, x, b, c = var('n,x,b,c')                                           # needs sage.symbolic
            sage: TestSuite(meixner).run()
            sage: TestSuite(meixner(3, x, b, c)).run()                                  # needs sage.symbolic
            sage: TestSuite(meixner(n, x, b, c)).run()                                  # needs sage.symbolic
        """
        super().__init__(name='meixner', nargs=4, latex_name='M')

    def eval_formula(self, n, x, b, c):
        r"""
        Evaluate ``self`` using an explicit formula.

        EXAMPLES::

            sage: x, b, c = var('x,b,c')                                                # needs sage.symbolic
            sage: meixner.eval_formula(3, x, b, c).expand().collect(x)                  # needs sage.symbolic
            -x^3*(3/c - 3/c^2 + 1/c^3 - 1) + b^3
             + 3*(b - 2*b/c + b/c^2 - 1/c - 1/c^2 + 1/c^3 + 1)*x^2 + 3*b^2
             + (3*b^2 + 6*b - 3*b^2/c - 3*b/c - 3*b/c^2 - 2/c^3 + 2)*x + 2*b
        """
        from sage.misc.misc_c import prod

        def P(val, k):
            return prod(val + j for j in range(k))
        return sum((-1)**k * binomial(n, k) * binomial(x, k) * factorial(k)
                   * P(x + b, n - k) * c**-k
                   for k in range(n+1))

    def _eval_(self, n, x, b, c, *args, **kwds):
        r"""
        Return an evaluation of the Meixner polynomial `M_n(x; b, c)`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: n, x, b, c = var('n,x,b,c')
            sage: meixner(2, x, b, c).collect(x)
            -x^2*(2/c - 1/c^2 - 1) + b^2 + (2*b - 2*b/c - 1/c^2 + 1)*x + b
            sage: meixner(3, x, b, c).factor().collect(x)
            -x^3*(3/c - 3/c^2 + 1/c^3 - 1) + b^3
             + 3*(b - 2*b/c + b/c^2 - 1/c - 1/c^2 + 1/c^3 + 1)*x^2 + 3*b^2
             + (3*b^2 + 6*b - 3*b^2/c - 3*b/c - 3*b/c^2 - 2/c^3 + 2)*x + 2*b
            sage: meixner(n, x, b, c)
            gamma(b + n)*hypergeometric((-n, -x), (b,), -1/c + 1)/gamma(b)

            sage: # needs sage.symbolic
            sage: n3_hypergeo = meixner(n, x, b, c)(n=3).simplify_hypergeometric()
            sage: n3_hypergeo = n3_hypergeo.simplify_full()
            sage: bool(n3_hypergeo == meixner(3, x, b, c))
            True
            sage: n4_hypergeo = meixner(n, x, b, c)(n=4).simplify_hypergeometric()
            sage: n4_hypergeo = n4_hypergeo.simplify_full()
            sage: bool(n4_hypergeo == meixner(4, x, b, c))
            True

            sage: meixner(2, x, b, c, hold=True)                                        # needs sage.symbolic
            meixner(2, x, b, c)
        """
        if kwds.get('hold', False):
            return None
        if n not in ZZ or n < 0:
            from sage.functions.hypergeometric import hypergeometric
            from sage.functions.gamma import gamma
            return gamma(b + n) / gamma(b) * hypergeometric([-n, -x], [b], 1 - 1/c)
        try:
            return self.eval_formula(n, x, b, c)
        except (TypeError, ValueError):
            return self.eval_recursive(n, x, b, c)

    def eval_recursive(self, n, x, b, c, *args, **kwds):
        r"""
        Return the Meixner polynomial `M_n(x; b, c)` using the
        recursive formula.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: x, b, c = var('x,b,c')
            sage: meixner.eval_recursive(0, x, b, c)
            1
            sage: meixner.eval_recursive(1, x, b, c)
            -x*(1/c - 1) + b
            sage: meixner.eval_recursive(2, x, b, c).simplify_full().collect(x)
            -x^2*(2/c - 1/c^2 - 1) + b^2 + (2*b - 2*b/c - 1/c^2 + 1)*x + b
            sage: bool(meixner(2, x, b, c) == meixner.eval_recursive(2, x, b, c))
            True
            sage: bool(meixner(3, x, b, c) == meixner.eval_recursive(3, x, b, c))
            True
            sage: bool(meixner(4, x, b, c) == meixner.eval_recursive(4, x, b, c))
            True
            sage: M = matrix([[-1/2, -1], [1, 0]])
            sage: ret = meixner.eval_recursive(2, M, b, c).simplify_full().factor()
            sage: for i in range(2):  # make the output polynomials in 1/c
            ....:     for j in range(2):
            ....:         ret[i, j] = ret[i, j].collect(c)
            sage: ret
            [b^2 + 1/2*(2*b + 3)/c - 1/4/c^2 - 5/4    -2*b + (2*b - 1)/c + 3/2/c^2 - 1/2]
            [    2*b - (2*b - 1)/c - 3/2/c^2 + 1/2             b^2 + b + 2/c - 1/c^2 - 1]
        """
        if n == 0:
            return parent(x).one()
        elif n == 1:
            return (1 - 1/c) * x + b
        tm2 = (b+n-1) * (b+n-2) * (n - 1) * meixner.eval_recursive(n-2, x, b, c)
        tm1 = (b+n-1) * ((c-1) * x + n-1 + (n-1+b) * c) * meixner.eval_recursive(n-1, x, b, c)
        return (tm1 - tm2) / (c * (n - 1 + b))


meixner = Func_meixner()


class Func_hahn(OrthogonalFunction):
    r"""
    Hahn polynomials `Q_k(x; a, b, n)`.

    INPUT:

    - ``k`` -- the degree
    - ``x`` -- the independent variable `x`
    - ``a``, ``b`` -- the parameters `a`, `b`
    - ``n`` -- the number of discrete points

    EXAMPLES:

    We verify the orthogonality for `n = 3`::

        sage: # needs sage.symbolic
        sage: n = 2
        sage: a, b = SR.var('a,b')
        sage: def rho(k, a, b, n):
        ....:     return binomial(a + k, k) * binomial(b + n - k, n - k)
        sage: M = matrix([[sum(rho(k, a, b, n)
        ....:                  * hahn(i, k, a, b, n) * hahn(j, k, a, b, n)
        ....:                  for k in range(n + 1)).expand().factor()
        ....:              for i in range(n+1)] for j in range(n+1)])
        sage: M = M.factor()
        sage: P = rising_factorial
        sage: def diag(i, a, b, n):
        ....:    return ((-1)^i * factorial(i) * P(b + 1, i) * P(i + a + b + 1, n + 1)
        ....:            / (factorial(n) * (2*i + a + b + 1) * P(-n, i) * P(a + 1, i)))
        sage: all(M[i,i] == diag(i, a, b, n) for i in range(3))
        True
        sage: all(M[i,j] == 0 for i in range(3) for j in range(3) if i != j)
        True
    """
    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: k, x, a, b, n = var('k,x,a,b,n')                                      # needs sage.symbolic
            sage: TestSuite(hahn).run()
            sage: TestSuite(hahn(3, x, a, b, n)).run()                                  # needs sage.symbolic
            sage: TestSuite(hahn(k, x, a, b, n)).run(skip='_test_category')             # needs sage.symbolic
        """
        super().__init__(name='hahn', nargs=5, latex_name='Q')

    def eval_formula(self, k, x, a, b, n):
        r"""
        Evaluate ``self`` using an explicit formula.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: k, x, a, b, n = var('k,x,a,b,n')
            sage: Q2 = hahn.eval_formula(2, x, a, b, n).simplify_full()
            sage: Q2.coefficient(x^2).factor()
            (a + b + 4)*(a + b + 3)/((a + 2)*(a + 1)*(n - 1)*n)
            sage: Q2.coefficient(x).factor()
            -(2*a*n - a + b + 4*n)*(a + b + 3)/((a + 2)*(a + 1)*(n - 1)*n)
            sage: Q2(x=0)
            1
        """
        P = rising_factorial
        return sum(P(-k, i) * P(k+a+b+1, i) * P(-x, i) / (P(a+1, i) * P(-n, i) * factorial(i))
                   for i in range(k+1))

    def _eval_(self, k, x, a, b, n, *args, **kwds):
        r"""
        Return an evaluation of the Hahn polynomial `Q_k(x; a, b, n)`.

        EXAMPLES::

            sage: k, x, a, b, n = var('k,x,a,b,n')                                      # needs sage.symbolic
            sage: hahn(1, x, a, b, n).collect(x)                                        # needs sage.symbolic
            -(a + b + 2)*x/((a + 1)*n) + 1
            sage: hahn(k, x, a, b, n)                                                   # needs sage.symbolic
            hypergeometric((-k, a + b + k + 1, -x), (a + 1, -n), 1)

            sage: # needs sage.symbolic
            sage: k2_hypergeo = hahn(k, x, a, b, n)(k=2).simplify_hypergeometric()
            sage: bool(k2_hypergeo == hahn(2, x, a, b, n))
            True
            sage: k3_hypergeo = hahn(k, x, a, b, n)(k=3).simplify_hypergeometric()
            sage: bool(k3_hypergeo == hahn(3, x, a, b, n))
            True

            sage: hahn(2, x, a, b, n, hold=True)                                        # needs sage.symbolic
            hahn(2, x, a, b, n)
        """
        if kwds.get('hold', False):
            return None
        if k not in ZZ or k < 0:
            from sage.functions.hypergeometric import hypergeometric
            return hypergeometric([-k, k+a+b+1, -x], [a+1, -n], 1)
        try:
            return self.eval_formula(k, x, a, b, n)
        except (TypeError, ValueError):
            return self.eval_recursive(k, x, a, b, n)

    def eval_recursive(self, k, x, a, b, n, *args, **kwds):
        r"""
        Return the Hahn polynomial `Q_k(x; a, b, n)` using the
        recursive formula.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: x, a, b, n = var('x,a,b,n')
            sage: hahn.eval_recursive(0, x, a, b, n)
            1
            sage: hahn.eval_recursive(1, x, a, b, n)
            -(a + b + 2)*x/((a + 1)*n) + 1
            sage: bool(hahn(2, x, a, b, n) == hahn.eval_recursive(2, x, a, b, n))
            True
            sage: bool(hahn(3, x, a, b, n) == hahn.eval_recursive(3, x, a, b, n))
            True
            sage: bool(hahn(4, x, a, b, n) == hahn.eval_recursive(4, x, a, b, n))
            True
            sage: M = matrix([[-1/2, -1], [1, 0]])                                      # needs sage.modules
            sage: ret = hahn.eval_recursive(2, M, 1, 2, n).simplify_full().factor()     # needs sage.modules
            sage: ret                                                                   # needs sage.modules
            [1/4*(4*n^2 + 8*n - 19)/((n - 1)*n)          3/2*(4*n + 3)/((n - 1)*n)]
            [        -3/2*(4*n + 3)/((n - 1)*n)          (n^2 - n - 7)/((n - 1)*n)]
        """
        if k == 0:
            return parent(x).one()
        elif k == 1:
            return -(a+b+2) / ((a+1)*n) * x + 1
        A = (k+a+b) * (k+a) * (n-k+1) / ((2*k+a+b-1) * (2*k+a+b))
        C = (k-1) * (k+b-1) * (k+a+b+n) / ((2*k+a+b-2) * (2*k+a+b-1))
        Hm1 = (-x + A + C) * hahn.eval_recursive(k-1, x, a, b, n)
        Hm2 = C * hahn.eval_recursive(k-2, x, a, b, n)
        return (Hm1 - Hm2) / A


hahn = Func_hahn()
