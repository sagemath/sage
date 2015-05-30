r"""
Orthogonal Polynomials

This module wraps some of the orthogonal/special functions in the
Maxima package "orthopoly". This package was written by Barton
Willis of the University of Nebraska at Kearney. It is released
under the terms of the General Public License (GPL). Send
Maxima-related bug reports and comments on this module to
willisb@unk.edu. In your report, please include Maxima and specfun
version information.


-  The Chebyshev polynomial of the first kind arises as a solution
   to the differential equation

   .. math::

         (1-x^2)\,y'' - x\,y' + n^2\,y = 0


   and those of the second kind as a solution to

   .. math::

         (1-x^2)\,y'' - 3x\,y' + n(n+2)\,y = 0.


   The Chebyshev polynomials of the first kind are defined by the
   recurrence relation

   .. math::

     T_0(x) = 1 \, T_1(x) = x \, T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x). \,


   The Chebyshev polynomials of the second kind are defined by the
   recurrence relation

   .. math::

     U_0(x) = 1 \, U_1(x) = 2x \, U_{n+1}(x) = 2xU_n(x) - U_{n-1}(x). \,



   For integers `m,n`, they satisfy the orthogonality
   relations

   .. math::

     \int_{-1}^1 T_n(x)T_m(x)\,\frac{dx}{\sqrt{1-x^2}} =\left\{ \begin{matrix} 0 &: n\ne m~~~~~\\ \pi &: n=m=0\\ \pi/2 &: n=m\ne 0 \end{matrix} \right.


   and


   .. math::

     \int_{-1}^1 U_n(x)U_m(x)\sqrt{1-x^2}\,dx =\frac{\pi}{2}\delta_{m,n}.



   They are named after Pafnuty Chebyshev (alternative
   transliterations: Tchebyshef or Tschebyscheff).

-  The Hermite polynomials are defined either by

   .. math::

     H_n(x)=(-1)^n e^{x^2/2}\frac{d^n}{dx^n}e^{-x^2/2}


   (the "probabilists' Hermite polynomials"), or by


   .. math::

     H_n(x)=(-1)^n e^{x^2}\frac{d^n}{dx^n}e^{-x^2}


   (the "physicists' Hermite polynomials"). Sage (via Maxima)
   implements the latter flavor. These satisfy the orthogonality
   relation

   .. math::

     \int_{-\infty}^\infty H_n(x)H_m(x)\,e^{-x^2}\,dx ={n!2^n}{\sqrt{\pi}}\delta_{nm}



   They are named in honor of Charles Hermite.

-  Each *Legendre polynomial* `P_n(x)` is an `n`-th degree polynomial.
   It may be expressed using Rodrigues' formula:

   .. math::

      P_n(x) = (2^n n!)^{-1} {\frac{d^n}{dx^n} } \left[ (x^2 -1)^n \right].

   These are solutions to Legendre's differential equation:

   .. math::

      {\frac{d}{dx}} \left[ (1-x^2) {\frac{d}{dx}} P(x) \right] + n(n+1)P(x) = 0.

   and satisfy the orthogonality relation

   .. math::

      \int_{-1}^{1} P_m(x) P_n(x)\,dx = {\frac{2}{2n + 1}} \delta_{mn}

   The *Legendre function of the second kind* `Q_n(x)` is another
   (linearly independent) solution to the Legendre differential equation.
   It is not an "orthogonal polynomial" however.

   The associated Legendre functions of the first kind
   `P_\ell^m(x)` can be given in terms of the "usual"
   Legendre polynomials by

   .. math::

     \begin{array}{ll} P_\ell^m(x)    &=  (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}P_\ell(x) \\ &=  \frac{(-1)^m}{2^\ell \ell!} (1-x^2)^{m/2}\frac{d^{\ell+m}}{dx^{\ell+m}}(x^2-1)^\ell. \end{array}


   Assuming `0 \le m \le \ell`, they satisfy the orthogonality
   relation:

   .. math::

      \int_{-1}^{1} P_k ^{(m)} P_\ell ^{(m)} dx  = \frac{2 (\ell+m)!}{(2\ell+1)(\ell-m)!}\ \delta _{k,\ell},


   where `\delta _{k,\ell}` is the Kronecker delta.

   The associated Legendre functions of the second kind
   `Q_\ell^m(x)` can be given in terms of the "usual"
   Legendre polynomials by


   .. math::

     Q_\ell^m(x)   =  (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}Q_\ell(x).



   They are named after Adrien-Marie Legendre.

-  Laguerre polynomials may be defined by the Rodrigues formula

   .. math::

      L_n(x)=\frac{e^x}{n!}\frac{d^n}{dx^n}\left(e^{-x} x^n\right).


   They are solutions of Laguerre's equation:


   .. math::

      x\,y'' + (1 - x)\,y' + n\,y = 0\,

   and satisfy the orthogonality relation


   .. math::

      \int_0^\infty L_m(x) L_n(x) e^{-x}\,dx = \delta_{mn}.



   The generalized Laguerre polynomials may be defined by the
   Rodrigues formula:


   .. math::

       L_n^{(\alpha)}(x)   = {\frac{x^{-\alpha} e^x}{n!}}{\frac{d^n}{dx^n}} \left(e^{-x} x^{n+\alpha}\right) .


   (These are also sometimes called the associated Laguerre
   polynomials.) The simple Laguerre polynomials are recovered from
   the generalized polynomials by setting `\alpha =0`.

   They are named after Edmond Laguerre.

-  Jacobi polynomials are a class of orthogonal polynomials. They
   are obtained from hypergeometric series in cases where the series
   is in fact finite:

   .. math::

     P_n^{(\alpha,\beta)}(z) =\frac{(\alpha+1)_n}{n!} \,_2F_1\left(-n,1+\alpha+\beta+n;\alpha+1;\frac{1-z}{2}\right) ,


   where `()_n` is Pochhammer's symbol (for the rising
   factorial), (Abramowitz and Stegun p561.) and thus have the
   explicit expression


   .. math::

     P_n^{(\alpha,\beta)} (z) = \frac{\Gamma (\alpha+n+1)}{n!\Gamma (\alpha+\beta+n+1)} \sum_{m=0}^n {n\choose m} \frac{\Gamma (\alpha + \beta + n + m + 1)}{\Gamma (\alpha + m + 1)} \left(\frac{z-1}{2}\right)^m .



   They are named after Carl Jacobi.

-  Ultraspherical or Gegenbauer polynomials are given in terms of
   the Jacobi polynomials `P_n^{(\alpha,\beta)}(x)` with
   `\alpha=\beta=a-1/2` by


   .. math::

     C_n^{(a)}(x)= \frac{\Gamma(a+1/2)}{\Gamma(2a)}\frac{\Gamma(n+2a)}{\Gamma(n+a+1/2)} P_n^{(a-1/2,a-1/2)}(x).


   They satisfy the orthogonality relation

   .. math::

     \int_{-1}^1(1-x^2)^{a-1/2}C_m^{(a)}(x)C_n^{(a)}(x)\, dx =\delta_{mn}2^{1-2a}\pi \frac{\Gamma(n+2a)}{(n+a)\Gamma^2(a)\Gamma(n+1)} ,


   for `a>-1/2`. They are obtained from hypergeometric series
   in cases where the series is in fact finite:


   .. math::

     C_n^{(a)}(z) =\frac{(2a)^{\underline{n}}}{n!} \,_2F_1\left(-n,2a+n;a+\frac{1}{2};\frac{1-z}{2}\right)


   where `\underline{n}` is the falling factorial. (See
   Abramowitz and Stegun p561)

   They are named for Leopold Gegenbauer (1849-1903).


For completeness, the Pochhammer symbol, introduced by Leo August
Pochhammer, `(x)_n`, is used in the theory of special
functions to represent the "rising factorial" or "upper factorial"

.. math::

         (x)_n=x(x+1)(x+2)\cdots(x+n-1)=\frac{(x+n-1)!}{(x-1)!}.


On the other hand, the "falling factorial" or "lower factorial" is

.. math::

     x^{\underline{n}}=\frac{x!}{(x-n)!} ,


in the notation of Ronald L. Graham, Donald E. Knuth and Oren
Patashnik in their book Concrete Mathematics.

.. note::

   The first call of any of these will usually cost a bit extra
   (it loads "specfun", but I'm not sure if that is the real reason).
   The next call is usually faster but not always.

.. TODO::

    Implement associated Legendre polynomials and Zernike
    polynomials. (Neither is in Maxima.)
    :wikipedia:`Associated_Legendre_polynomials`
    :wikipedia:`Zernike_polynomials`

REFERENCES:

.. [ASHandbook] Abramowitz and Stegun: Handbook of Mathematical Functions,
    http://www.math.sfu.ca/ cbm/aands/

.. :wikipedia:`Chebyshev_polynomials`

.. :wikipedia:`Legendre_polynomials`

.. :wikipedia:`Hermite_polynomials`

.. http://mathworld.wolfram.com/GegenbauerPolynomial.html

.. :wikipedia:`Jacobi_polynomials`

.. :wikipedia:`Laguerre_polynomia`

.. :wikipedia:`Associated_Legendre_polynomials`

.. [EffCheby] Wolfram Koepf: Effcient Computation of Chebyshev Polynomials
    in Computer Algebra
    Computer Algebra Systems: A Practical Guide.
    John Wiley, Chichester (1999): 79-99.

AUTHORS:

- David Joyner (2006-06)
- Stefan Reiterer (2010-)
"""

#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import warnings

from sage.misc.sage_eval import sage_eval
from sage.rings.all import ZZ, RR, CC
from sage.rings.real_mpfr import is_RealField
from sage.rings.complex_field import is_ComplexField
from sage.calculus.calculus import maxima


from sage.symbolic.function import BuiltinFunction
from sage.symbolic.expression import is_Expression
from sage.functions.other import factorial, binomial
from sage.structure.all import parent

_done = False
def _init():
    """
    Internal function which checks if Maxima has loaded the
    "orthopoly" package.  All functions using this in this
    file should call this function first.

    TEST:

    The global starts ``False``::

        sage: sage.functions.orthogonal_polys._done
        False

    Then after using one of these functions, it changes::

        sage: from sage.functions.orthogonal_polys import legendre_P
        sage: legendre_P(2,x)
        3/2*(x - 1)^2 + 3*x - 2
        sage: sage.functions.orthogonal_polys._done
        True


    Note that because here we use a Pynac variable ``x``,
    the representation of the function is different from
    its actual doctest, where a polynomial indeterminate
    ``x`` is used.
    """
    global _done
    if _done:
        return
    maxima.eval('load("orthopoly");')
    # TODO -- make it possible to use the intervals returned
    # instead of just discarding this info!
    maxima.eval('orthopoly_returns_intervals:false;')
    _done = True


class OrthogonalPolynomial(BuiltinFunction):
    """
    Base class for orthogonal polynomials.

    This class is an abstract base class for all orthogonal polynomials since
    they share similar properties. The evaluation as a polynomial
    is either done via maxima, or with pynac.

    Convention: The first argument is always the order of the polynomial,
    the others are other values or parameters where the polynomial is
    evaluated.
    """
    def __init__(self, name, nargs=2, latex_name=None, conversions={}):
        """
        :class:`OrthogonalPolynomial` class needs the same input parameter as
        it's parent class.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import OrthogonalPolynomial
            sage: new = OrthogonalPolynomial('testo_P')
            sage: new
            testo_P
        """
        try:
            self._maxima_name = conversions['maxima']
        except KeyError:
            self._maxima_name = None

        super(OrthogonalPolynomial,self).__init__(name=name, nargs=nargs,
                                 latex_name=latex_name, conversions=conversions)

    def _maxima_init_evaled_(self, *args):
        r"""
        Return a string which represents this function evaluated at
        ``n, x`` in Maxima.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import OrthogonalPolynomial
            sage: P = OrthogonalPolynomial('testo_P')
            sage: P._maxima_init_evaled_(2, 5) is None
            True
        """
        return None

    def eval_formula(self, *args):
        """
        Evaluate this polynomial using an explicit formula.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import OrthogonalPolynomial
            sage: P = OrthogonalPolynomial('testo_P')
            sage: P.eval_formula(1,2.0)
            Traceback (most recent call last):
            ...
            NotImplementedError: no explicit calculation of values implemented
        """
        raise NotImplementedError("no explicit calculation of values implemented")

    def _eval_special_values_(self, *args):
        """
        Evaluate the polynomial explicitly for special values.

        EXAMPLES::

            sage: var('n')
            n
            sage: chebyshev_T(n,-1)
            (-1)^n
        """
        raise ValueError("no special values known")

    def _eval_(self, n, *args):
        """
        The :meth:`_eval_()` method decides which evaluation suits best
        for the given input, and returns a proper value.

        EXAMPLES::

            sage: var('n,x')
            (n, x)
            sage: chebyshev_T(5,x)
            16*x^5 - 20*x^3 + 5*x
        """
        return None

    def __call__(self, n, *args, **kwds):
        """
        This overides the call method from SageObject to avoid problems with coercions,
        since the _eval_ method is able to handle more data types than symbolic functions
        would normally allow.
        Thus we have the distinction between algebraic objects (if n is an integer),
        and else as symbolic function.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-x-1)
            sage: chebyshev_T(5, a)
            16*a^2 + a - 4
        """
        return super(OrthogonalPolynomial,self).__call__(n, *args, **kwds)

class ChebyshevPolynomial(OrthogonalPolynomial):
    """
    Abstract base class for Chebyshev polynomials of the first and second kind.

    EXAMPLES::

        sage: chebyshev_T(3,x)
        4*x^3 - 3*x
    """
    def __call__(self, n, *args, **kwds):
        """
        This overides the call method from SageObject to avoid problems with coercions,
        since the _eval_ method is able to handle more data types than symbolic functions
        would normally allow.
        Thus we have the distinction between algebraic objects (if n is an integer),
        and else as symbolic function.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-x-1)
            sage: chebyshev_T(5, a)
            16*a^2 + a - 4
            sage: chebyshev_T(5,MatrixSpace(ZZ, 2)([1, 2, -4, 7]))
            [-40799  44162]
            [-88324  91687]
            sage: R.<x> = QQ[]
            sage: parent(chebyshev_T(5, x))
            Univariate Polynomial Ring in x over Rational Field
            sage: chebyshev_T(5, 2, hold=True)
            chebyshev_T(5, 2)
            sage: chebyshev_T(1,2,3)
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

        return super(ChebyshevPolynomial,self).__call__(n, *args, **kwds)

    def _eval_(self, n, x):
        """
        The :meth:`_eval_()` method decides which evaluation suits best
        for the given input, and returns a proper value.

        EXAMPLES::

            sage: var('n,x')
            (n, x)
            sage: chebyshev_T(5,x)
            16*x^5 - 20*x^3 + 5*x
            sage: chebyshev_T(64, x)
            2*(2*(2*(2*(2*(2*x^2 - 1)^2 - 1)^2 - 1)^2 - 1)^2 - 1)^2 - 1
            sage: chebyshev_T(n,-1)
            (-1)^n
            sage: chebyshev_T(-7,x)
            64*x^7 - 112*x^5 + 56*x^3 - 7*x
            sage: chebyshev_T(3/2,x)
            chebyshev_T(3/2, x)
            sage: R.<t> = QQ[]
            sage: chebyshev_T(2,t)
            2*t^2 - 1
            sage: chebyshev_U(2,t)
            4*t^2 - 1
            sage: parent(chebyshev_T(4, RIF(5)))
            Real Interval Field with 53 bits of precision
            sage: RR2 = RealField(5)
            sage: chebyshev_T(100000,RR2(2))
            8.9e57180
            sage: chebyshev_T(5,Qp(3)(2))
            2 + 3^2 + 3^3 + 3^4 + 3^5 + O(3^20)
            sage: chebyshev_T(100001/2, 2)
            doctest:...: RuntimeWarning: mpmath failed, keeping expression unevaluated
            chebyshev_T(100001/2, 2)
            sage: chebyshev_U._eval_(1.5, Mod(8,9)) is None
            True
        """
        # n is an integer => evaluate algebraically (as polynomial)
        if n in ZZ:
            n = ZZ(n)
            # Expanded symbolic expression only for small values of n
            if is_Expression(x) and n.abs() < 32:
                return self.eval_formula(n, x)
            return self.eval_algebraic(n, x)

        if is_Expression(x) or is_Expression(n):
            # Check for known identities
            try:
                return self._eval_special_values_(n, x)
            except ValueError:
                # Don't evaluate => keep symbolic
                return None

        # n is not an integer and neither n nor x is symbolic.
        # We assume n and x are real/complex and evaluate numerically
        try:
            import sage.libs.mpmath.all as mpmath
            return self._evalf_(n, x)
        except mpmath.NoConvergence:
            warnings.warn("mpmath failed, keeping expression unevaluated",
                          RuntimeWarning)
            return None
        except Exception:
            # Numerical evaluation failed => keep symbolic
            return None

    
class Func_chebyshev_T(ChebyshevPolynomial):
    """
    Chebyshev polynomials of the first kind.

    REFERENCE:

    - [ASHandbook]_ 22.5.31 page 778 and 6.1.22 page 256.

    EXAMPLES::

       sage: chebyshev_T(5,x)
       16*x^5 - 20*x^3 + 5*x
       sage: var('k')
       k
       sage: test = chebyshev_T(k,x)
       sage: test
       chebyshev_T(k, x)
    """
    def __init__(self):
        """
        Init method for the chebyshev polynomials of the first kind.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import Func_chebyshev_T
            sage: chebyshev_T2 = Func_chebyshev_T()
            sage: chebyshev_T2(1,x)
            x
        """
        ChebyshevPolynomial.__init__(self, "chebyshev_T", nargs=2,
                                     conversions=dict(maxima='chebyshev_t',
                                                      mathematica='ChebyshevT'))

    def _eval_special_values_(self, n, x):
        """
        Values known for special values of x.
        For details see [ASHandbook]_ 22.4 (p. 777)

        EXAMPLES:

            sage: var('n')
            n
            sage: chebyshev_T(n,1)
            1
            sage: chebyshev_T(n,0)
            1/2*(-1)^(1/2*n)*((-1)^n + 1)
            sage: chebyshev_T(n,-1)
            (-1)^n
            sage: chebyshev_T._eval_special_values_(3/2,x)
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

            sage: chebyshev_T._evalf_(10,3)
            2.26195370000000e7
            sage: chebyshev_T._evalf_(10,3,parent=RealField(75))
            2.261953700000000000000e7
            sage: chebyshev_T._evalf_(10,I)
            -3363.00000000000
            sage: chebyshev_T._evalf_(5,0.3)
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

            sage: chebyshev_T(1234.5, RDF(2.1))
            5.48174256255782e735
            sage: chebyshev_T(1234.5, I)
            -1.21629397684152e472 - 1.21629397684152e472*I

        For large values of ``n``, mpmath fails (but the algebraic formula
        still works)::

            sage: chebyshev_T._evalf_(10^6, 0.1)
            Traceback (most recent call last):
            ...
            NoConvergence: Hypergeometric series converges too slowly. Try increasing maxterms.
            sage: chebyshev_T(10^6, 0.1)
            0.636384327171504
        """
        try:
            real_parent = kwds['parent']
        except KeyError:
            real_parent = parent(x)

            if not is_RealField(real_parent) and not is_ComplexField(real_parent):
                # parent is not a real or complex field: figure out a good parent
                if x in RR:
                    x = RR(x)
                    real_parent = RR
                elif x in CC:
                    x = CC(x)
                    real_parent = CC

        if not is_RealField(real_parent) and not is_ComplexField(real_parent):
            raise TypeError("cannot evaluate chebyshev_T with parent {}".format(real_parent))

        from sage.libs.mpmath.all import call as mpcall
        from sage.libs.mpmath.all import chebyt as mpchebyt

        return mpcall(mpchebyt, n, x, parent=real_parent)

    def _maxima_init_evaled_(self, n, x):
        """
        Evaluate the Chebyshev polynomial ``self`` with maxima.

        EXAMPLES::

            sage: var('n, x')
            (n, x)
            sage: chebyshev_T._maxima_init_evaled_(1,x)
            '_SAGE_VAR_x'
            sage: maxima(chebyshev_T(n, chebyshev_T(n, x)))
            chebyshev_t(_SAGE_VAR_n,chebyshev_t(_SAGE_VAR_n,_SAGE_VAR_x))
        """
        return maxima.eval('chebyshev_t({0},{1})'.format(n._maxima_init_(), x._maxima_init_()))

    def eval_formula(self, n, x):
        """
        Evaluate ``chebyshev_T`` using an explicit formula.
        See [ASHandbook]_ 227 (p. 782) for details for the recurions.
        See also [EffCheby]_ for fast evaluation techniques.

        INPUT:

        - ``n`` -- an integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

            sage: chebyshev_T.eval_formula(-1,x)
            x
            sage: chebyshev_T.eval_formula(0,x)
            1
            sage: chebyshev_T.eval_formula(1,x)
            x
            sage: chebyshev_T.eval_formula(2,0.1) == chebyshev_T._evalf_(2,0.1)
            True
            sage: chebyshev_T.eval_formula(10,x)
            512*x^10 - 1280*x^8 + 1120*x^6 - 400*x^4 + 50*x^2 - 1
            sage: chebyshev_T.eval_algebraic(10,x).expand()
            512*x^10 - 1280*x^8 + 1120*x^6 - 400*x^4 + 50*x^2 - 1
        """
        if n < 0:
            return self.eval_formula(-n, x)
        elif n == 0:
            return parent(x).one()

        res = parent(x).zero()
        for j in xrange(0, n//2+1):
            f = factorial(n-1-j) / factorial(j) / factorial(n-2*j)
            res += (-1)**j * (2*x)**(n-2*j) * f
        res *= n/2
        return res

    def eval_algebraic(self, n, x):
        """
        Evaluate :class:`chebyshev_T` as polynomial, using a recursive
        formula.

        INPUT:

        - ``n`` -- an integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

            sage: chebyshev_T.eval_algebraic(5, x)
            2*(2*(2*x^2 - 1)*x - x)*(2*x^2 - 1) - x
            sage: chebyshev_T(-7, x) - chebyshev_T(7,x)
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
            sage: n = 97; x = RIF(pi/2/n)
            sage: chebyshev_T(n, cos(x)).contains_zero()
            True
            sage: R.<t> = Zp(2, 8, 'capped-abs')[]
            sage: chebyshev_T(10^6+1, t)
            (2^7 + O(2^8))*t^5 + (O(2^8))*t^4 + (2^6 + O(2^8))*t^3 + (O(2^8))*t^2 + (1 + 2^6 + O(2^8))*t + (O(2^8))
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

            sage: chebyshev_T._eval_recursive_(5, x)
            (2*(2*(2*x^2 - 1)*x - x)*(2*x^2 - 1) - x, False)
            sage: chebyshev_T._eval_recursive_(5, x, True)
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

            sage: import numpy
            sage: z = numpy.array([1,2])
            sage: z2 = numpy.array([[1,2],[1,2]])
            sage: z3 = numpy.array([1,2,3.])
            sage: chebyshev_T(1,z)
            array([ 1.,  2.])
            sage: chebyshev_T(1,z2)
            array([[ 1.,  2.],
                   [ 1.,  2.]])
            sage: chebyshev_T(1,z3)
            array([ 1.,  2.,  3.])
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

            sage: var('k')
            k
            sage: derivative(chebyshev_T(k,x),x)
            k*chebyshev_U(k - 1, x)
            sage: derivative(chebyshev_T(3,x),x)
            12*x^2 - 3
            sage: derivative(chebyshev_T(k,x),k)
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

class Func_chebyshev_U(ChebyshevPolynomial):
    """
    Class for the Chebyshev polynomial of the second kind.

    REFERENCE:

    - [ASHandbook]_ 22.8.3 page 783 and 6.1.22 page 256.

    EXAMPLES::

        sage: R.<t> = QQ[]
        sage: chebyshev_U(2,t)
        4*t^2 - 1
        sage: chebyshev_U(3,t)
        8*t^3 - 4*t
    """
    def __init__(self):
        """
        Init method for the chebyshev polynomials of the second kind.

        EXAMPLES::

            sage: from sage.functions.orthogonal_polys import Func_chebyshev_U
            sage: chebyshev_U2 = Func_chebyshev_U()
            sage: chebyshev_U2(1,x)
            2*x
        """
        ChebyshevPolynomial.__init__(self, "chebyshev_U", nargs=2,
                                     conversions=dict(maxima='chebyshev_u',
                                                      mathematica='ChebyshevU'))

    def eval_formula(self, n, x):
        """
        Evaluate ``chebyshev_U`` using an explicit formula.
        See [ASHandbook]_ 227 (p. 782) for details on the recurions.
        See also [EffCheby]_ for the recursion formulas.

        INPUT:

        - ``n`` -- an integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

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
            sage: chebyshev_U.eval_formula(2,0.1) == chebyshev_U._evalf_(2,0.1)
            True
        """
        if n < -1:
            return -self.eval_formula(-n-2, x)

        res = parent(x).zero()
        for j in xrange(0, n//2+1):
            f = binomial(n-j, j)
            res += (-1)**j * (2*x)**(n-2*j) * f
        return res

    def eval_algebraic(self, n, x):
        """
        Evaluate :class:`chebyshev_U` as polynomial, using a recursive
        formula.

        INPUT:

        - ``n`` -- an integer

        - ``x`` -- a value to evaluate the polynomial at (this can be
          any ring element)

        EXAMPLES::

            sage: chebyshev_U.eval_algebraic(5,x)
            -2*((2*x + 1)*(2*x - 1)*x - 4*(2*x^2 - 1)*x)*(2*x + 1)*(2*x - 1)
            sage: parent(chebyshev_U(3, Mod(8,9)))
            Ring of integers modulo 9
            sage: parent(chebyshev_U(3, Mod(1,9)))
            Ring of integers modulo 9
            sage: chebyshev_U(-3,x) + chebyshev_U(1,x)
            0
            sage: chebyshev_U(-1,Mod(5,8))
            0
            sage: parent(chebyshev_U(-1,Mod(5,8)))
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
            sage: n = 97; x = RIF(pi/n)
            sage: chebyshev_U(n-1, cos(x)).contains_zero()
            True
            sage: R.<t> = Zp(2, 6, 'capped-abs')[]
            sage: chebyshev_U(10^6+1, t)
            (2 + O(2^6))*t + (O(2^6))
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

            sage: chebyshev_U._eval_recursive_(3, x)
            (4*((2*x + 1)*(2*x - 1) - 2*x^2)*x, False)
            sage: chebyshev_U._eval_recursive_(3, x, True)
            (4*((2*x + 1)*(2*x - 1) - 2*x^2)*x, ((2*x + 1)*(2*x - 1) + 2*x)*((2*x + 1)*(2*x - 1) - 2*x))
        """
        if n == 0:
            return parent(x).one(), 2*x

        assert n >= 1
        a, b = self._eval_recursive_((n-1)//2, x, True)
        if n % 2 == 0:
            return (b+a)*(b-a), both and 2*b*(x*b-a)
        else:
            return 2*a*(b-x*a), both and (b+a)*(b-a)

    def _maxima_init_evaled_(self, n, x):
        """
        Uses maxima to evaluate ``self``.

        EXAMPLES::

            sage: var('n, x')
            (n, x)
            sage: maxima(chebyshev_U(5,x))
            32*_SAGE_VAR_x^5-32*_SAGE_VAR_x^3+6*_SAGE_VAR_x
            sage: maxima(chebyshev_U(n,x))
            chebyshev_u(_SAGE_VAR_n,_SAGE_VAR_x)
            sage: maxima(chebyshev_U(2,x))
            4*_SAGE_VAR_x^2-1
        """
        return maxima.eval('chebyshev_u({0},{1})'.format(n._maxima_init_(), x._maxima_init_()))

    def _evalf_(self, n, x, **kwds):
        """
        Evaluate :class:`chebyshev_U` numerically with mpmath.

        EXAMPLES::

            sage: chebyshev_U(5,-4+3.*I)
            98280.0000000000 - 11310.0000000000*I
            sage: chebyshev_U(10,3).n(75)
            4.661117900000000000000e7
            sage: chebyshev_U._evalf_(1.5, Mod(8,9))
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate chebyshev_U with parent Ring of integers modulo 9
        """
        try:
            real_parent = kwds['parent']
        except KeyError:
            real_parent = parent(x)

            if not is_RealField(real_parent) and not is_ComplexField(real_parent):
                # parent is not a real or complex field: figure out a good parent
                if x in RR:
                    x = RR(x)
                    real_parent = RR
                elif x in CC:
                    x = CC(x)
                    real_parent = CC

        if not is_RealField(real_parent) and not is_ComplexField(real_parent):
            raise TypeError("cannot evaluate chebyshev_U with parent {}".format(real_parent))

        from sage.libs.mpmath.all import call as mpcall
        from sage.libs.mpmath.all import chebyu as mpchebyu

        return mpcall(mpchebyu, n, x, parent=real_parent)

    def _eval_special_values_(self, n, x):
        """
        Values known for special values of x.
        See [ASHandbook]_ 22.4 (p.777).

        EXAMPLES::

            sage: var('n')
            n
            sage: chebyshev_U(n,1)
            n + 1
            sage: chebyshev_U(n,0)
            1/2*(-1)^(1/2*n)*((-1)^n + 1)
            sage: chebyshev_U(n,-1)
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

            sage: import numpy
            sage: z = numpy.array([1,2])
            sage: z2 = numpy.array([[1,2],[1,2]])
            sage: z3 = numpy.array([1,2,3.])
            sage: chebyshev_U(1,z)
            array([ 2.,  4.])
            sage: chebyshev_U(1,z2)
            array([[ 2.,  4.],
                   [ 2.,  4.]])
            sage: chebyshev_U(1,z3)
            array([ 2.,  4.,  6.])
            sage: chebyshev_U(z,0.1)
            array([ 0.2 , -0.96])
        """
        from scipy.special import eval_chebyu
        return eval_chebyu(n, x)

    def _derivative_(self, n, x, diff_param):
        """
        Return the derivative of :class:`chebyshev_U` in form of the Chebyshev
        polynomials of the first and second kind.

        EXAMPLES::

            sage: var('k')
            k
            sage: derivative(chebyshev_U(k,x),x)
            ((k + 1)*chebyshev_T(k + 1, x) - x*chebyshev_U(k, x))/(x^2 - 1)
            sage: derivative(chebyshev_U(3,x),x)
            24*x^2 - 4
            sage: derivative(chebyshev_U(k,x),k)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative w.r.t. to the index is not supported yet
        """
        if diff_param == 0:
            raise NotImplementedError("derivative w.r.t. to the index is not supported yet")
        elif diff_param == 1:
            return ((n+1)*chebyshev_T(n+1, x) - x*chebyshev_U(n,x)) / (x*x-1)
        raise ValueError("illegal differentiation parameter {}".format(diff_param))

chebyshev_U = Func_chebyshev_U()


def gen_laguerre(n,a,x):
    """
    Returns the generalized Laguerre polynomial for integers `n > -1`.
    Typically, `a = 1/2` or `a = -1/2`.

    REFERENCES:

    - Table on page 789 in [ASHandbook]_.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: gen_laguerre(2,1,x)
        1/2*x^2 - 3*x + 3
        sage: gen_laguerre(2,1/2,x)
        1/2*x^2 - 5/2*x + 15/8
        sage: gen_laguerre(2,-1/2,x)
        1/2*x^2 - 3/2*x + 3/8
        sage: gen_laguerre(2,0,x)
        1/2*x^2 - 2*x + 1
        sage: gen_laguerre(3,0,x)
        -1/6*x^3 + 3/2*x^2 - 3*x + 1

    Check that :trac:`17192` is fixed::
        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: gen_laguerre(0,1,x)
        1

        sage: gen_laguerre(-1,1,x)
        Traceback (most recent call last):
        ...
        ValueError: n must be greater than -1, got n = -1

        sage: gen_laguerre(-7,1,x)
        Traceback (most recent call last):
        ...
        ValueError: n must be greater than -1, got n = -7
    """
    if not (n > -1):
        raise ValueError("n must be greater than -1, got n = {0}".format(n))

    _init()
    return sage_eval(maxima.eval('gen_laguerre(%s,%s,x)'%(ZZ(n),a)), locals={'x':x})

def gen_legendre_P(n,m,x):
    r"""
    Returns the generalized (or associated) Legendre function of the
    first kind.

    The awkward code for when m is odd and 1 results from the fact that
    Maxima is happy with, for example, `(1 - t^2)^3/2`, but
    Sage is not. For these cases the function is computed from the
    (m-1)-case using one of the recursions satisfied by the Legendre
    functions.

    REFERENCE:

    - Gradshteyn and Ryzhik 8.706 page 1000.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: gen_legendre_P(2, 0, t)
        3/2*t^2 - 1/2
        sage: gen_legendre_P(2, 0, t) == legendre_P(2, t)
        True
        sage: gen_legendre_P(3, 1, t)
        -3/2*(5*t^2 - 1)*sqrt(-t^2 + 1)
        sage: gen_legendre_P(4, 3, t)
        105*(t^2 - 1)*sqrt(-t^2 + 1)*t
        sage: gen_legendre_P(7, 3, I).expand()
        -16695*sqrt(2)
        sage: gen_legendre_P(4, 1, 2.5)
        -583.562373654533*I
    """
    from sage.functions.all import sqrt
    _init()
    if m.mod(2).is_zero() or m.is_one():
        return sage_eval(maxima.eval('assoc_legendre_p(%s,%s,x)'%(ZZ(n),ZZ(m))), locals={'x':x})
    else:
        return sqrt(1-x**2)*(((n-m+1)*x*gen_legendre_P(n,m-1,x)-(n+m-1)*gen_legendre_P(n-1,m-1,x))/(1-x**2))

def gen_legendre_Q(n,m,x):
    """
    Returns the generalized (or associated) Legendre function of the
    second kind.

    Maxima restricts m = n. Hence the cases m n are computed using the
    same recursion used for gen_legendre_P(n,m,x) when m is odd and
    1.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: gen_legendre_Q(2,0,t)
        3/4*t^2*log(-(t + 1)/(t - 1)) - 3/2*t - 1/4*log(-(t + 1)/(t - 1))
        sage: gen_legendre_Q(2,0,t) - legendre_Q(2, t)
        0
        sage: gen_legendre_Q(3,1,0.5)
        2.49185259170895
        sage: gen_legendre_Q(0, 1, x)
        -1/sqrt(-x^2 + 1)
        sage: gen_legendre_Q(2, 4, x).factor()
        48*x/((x + 1)^2*(x - 1)^2)
    """
    from sage.functions.all import sqrt
    if m <= n:
        _init()
        return sage_eval(maxima.eval('assoc_legendre_q(%s,%s,x)'%(ZZ(n),ZZ(m))), locals={'x':x})
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
        return ((n-m+1)*x*gen_legendre_Q(n,m-1,x)-(n+m-1)*gen_legendre_Q(n-1,m-1,x))/sqrt(1-x**2)

def hermite(n,x):
    """
    Returns the Hermite polynomial for integers `n > -1`.

    REFERENCE:

    - [ASHandbook]_ 22.5.40 and 22.5.41, page 779.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: hermite(2,x)
        4*x^2 - 2
        sage: hermite(3,x)
        8*x^3 - 12*x
        sage: hermite(3,2)
        40
        sage: S.<y> = PolynomialRing(RR)
        sage: hermite(3,y)
        8.00000000000000*y^3 - 12.0000000000000*y
        sage: R.<x,y> = QQ[]
        sage: hermite(3,y^2)
        8*y^6 - 12*y^2
        sage: w = var('w')
        sage: hermite(3,2*w)
        8*(8*w^2 - 3)*w

    Check that :trac:`17192` is fixed::
        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: hermite(0,x)
        1

        sage: hermite(-1,x)
        Traceback (most recent call last):
        ...
        ValueError: n must be greater than -1, got n = -1

        sage: hermite(-7,x)
        Traceback (most recent call last):
        ...
        ValueError: n must be greater than -1, got n = -7
    """
    if not (n > -1):
        raise ValueError("n must be greater than -1, got n = {0}".format(n))

    _init()
    return sage_eval(maxima.eval('hermite(%s,x)'%ZZ(n)), locals={'x':x})

def jacobi_P(n,a,b,x):
    r"""
    Returns the Jacobi polynomial `P_n^{(a,b)}(x)` for
    integers `n > -1` and a and b symbolic or `a > -1`
    and `b > -1`. The Jacobi polynomials are actually defined
    for all a and b. However, the Jacobi polynomial weight
    `(1-x)^a(1+x)^b` isn't integrable for `a \leq -1`
    or `b \leq -1`.

    REFERENCE:

    - Table on page 789 in [ASHandbook]_.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: jacobi_P(2,0,0,x)
        3/2*x^2 - 1/2
        sage: jacobi_P(2,1,2,1.2)        # random output of low order bits
        5.009999999999998

    Check that :trac:`17192` is fixed::
        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: jacobi_P(0,0,0,x)
        1

        sage: jacobi_P(-1,0,0,x)
        Traceback (most recent call last):
        ...
        ValueError: n must be greater than -1, got n = -1

        sage: jacobi_P(-7,0,0,x)
        Traceback (most recent call last):
        ...
        ValueError: n must be greater than -1, got n = -7
    """
    if not (n > -1):
        raise  ValueError("n must be greater than -1, got n = {0}".format(n))

    _init()
    return sage_eval(maxima.eval('jacobi_p(%s,%s,%s,x)'%(ZZ(n),a,b)), locals={'x':x})

def laguerre(n,x):
    """
    Return the Laguerre polynomial for integers `n > -1`.

    REFERENCE:

    - [ASHandbook]_ 22.5.16, page 778 and page 789.

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: laguerre(2,x)
        1/2*x^2 - 2*x + 1
        sage: laguerre(3,x)
        -1/6*x^3 + 3/2*x^2 - 3*x + 1
        sage: laguerre(2,2)
        -1

    Check that :trac:`17192` is fixed::
        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: laguerre(0,x)
        1

        sage: laguerre(-1,x)
        Traceback (most recent call last):
        ...
        ValueError: n must be greater than -1, got n = -1

        sage: laguerre(-7,x)
        Traceback (most recent call last):
        ...
        ValueError: n must be greater than -1, got n = -7
    """
    if not (n > -1):
        raise ValueError("n must be greater than -1, got n = {0}".format(n))

    _init()
    return sage_eval(maxima.eval('laguerre(%s,x)'%ZZ(n)), locals={'x':x})

def legendre_P(n,x):
    """
    Returns the Legendre polynomial of the first kind.

    REFERENCE:

    - [ASHandbook]_ 22.5.35 page 779.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: legendre_P(2,t)
        3/2*t^2 - 1/2
        sage: legendre_P(3, 1.1)
        1.67750000000000
        sage: legendre_P(3, 1 + I)
        7/2*I - 13/2
        sage: legendre_P(3, MatrixSpace(ZZ, 2)([1, 2, -4, 7]))
        [-179  242]
        [-484  547]
        sage: legendre_P(3, GF(11)(5))
        8
    """
    _init()
    return sage_eval(maxima.eval('legendre_p(%s,x)'%ZZ(n)), locals={'x':x})

def legendre_Q(n,x):
    """
    Returns the Legendre function of the second kind.

    Computed using Maxima.

    EXAMPLES::

        sage: P.<t> = QQ[]
        sage: legendre_Q(2, t)
        3/4*t^2*log(-(t + 1)/(t - 1)) - 3/2*t - 1/4*log(-(t + 1)/(t - 1))
        sage: legendre_Q(3, 0.5)
        -0.198654771479482
        sage: legendre_Q(4, 2)
        443/16*I*pi + 443/16*log(3) - 365/12
        sage: legendre_Q(4, 2.0)
        0.00116107583162324 + 86.9828465962674*I
    """
    _init()
    return sage_eval(maxima.eval('legendre_q(%s,x)'%ZZ(n)), locals={'x':x})

def ultraspherical(n,a,x):
    """
    Returns the ultraspherical (or Gegenbauer) polynomial for integers
    `n > -1`.

    Computed using Maxima.

    REFERENCE:

    - [ASHandbook]_ 22.5.27

    EXAMPLES::

        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: ultraspherical(2,3/2,x)
        15/2*x^2 - 3/2
        sage: ultraspherical(2,1/2,x)
        3/2*x^2 - 1/2
        sage: ultraspherical(1,1,x)
        2*x
        sage: t = PolynomialRing(RationalField(),"t").gen()
        sage: gegenbauer(3,2,t)
        32*t^3 - 12*t

    Check that :trac:`17192` is fixed::
        sage: x = PolynomialRing(QQ, 'x').gen()
        sage: ultraspherical(0,1,x)
        1

        sage: ultraspherical(-1,1,x)
        Traceback (most recent call last):
        ...
        ValueError: n must be greater than -1, got n = -1

        sage: ultraspherical(-7,1,x)
        Traceback (most recent call last):
        ...
        ValueError: n must be greater than -1, got n = -7
    """
    if not (n > -1):
        raise ValueError("n must be greater than -1, got n = {0}".format(n))

    _init()
    return sage_eval(maxima.eval('ultraspherical(%s,%s,x)'%(ZZ(n),a)), locals={'x':x})

gegenbauer = ultraspherical
