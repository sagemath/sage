"""
Dense univariate polynomials over `\ZZ`, implemented using FLINT.

AUTHORS:

- David Harvey: rewrote to talk to NTL directly, instead of via ntl.pyx
  (2007-09); a lot of this was based on Joel Mohler's recent rewrite of the NTL
  wrapper
- David Harvey: split off from polynomial_element_generic.py (2007-09)
- Burcin Erocal: rewrote to use FLINT (2008-06-16)
"""

################################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2008-2010 Burcin Erocal <burcin@erocal.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
################################################################################

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/gmp.pxi"
include "sage/libs/ntl/decl.pxi"

from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.structure.element cimport ModuleElement, RingElement
from sage.structure.element import coerce_binop

from sage.rings.polynomial.polynomial_element cimport is_Polynomial

from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.libs.all import pari, pari_gen
from sage.structure.factorization import Factorization

from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.arith import lcm

from sage.libs.flint.fmpz_poly cimport fmpz_poly_reverse
from sage.libs.flint.ntl_interface cimport fmpz_poly_set_ZZX, fmpz_poly_get_ZZX
from sage.libs.ntl.ntl_ZZX_decl cimport *, vec_pair_ZZX_long_c

cdef extern from "limits.h":
    long LONG_MAX

cdef extern from "flint/flint.h":
    int FLINT_BITS

cdef class Polynomial_integer_dense_flint(Polynomial):
    r"""
    A dense polynomial over the integers, implemented via FLINT.
    """

    def __cinit__(self):
        r"""
        This calls the underlying FLINT fmpz_poly constructor
        """
        fmpz_poly_init(self.__poly)


    def __dealloc__(self):
        r"""
        calls the underlying FLINT fmpz_poly destructor
        """
        fmpz_poly_clear(self.__poly)

    cdef Polynomial_integer_dense_flint _new(self):
        r"""
        Quickly creates a new initialized Polynomial_integer_dense_flint
        with the correct parent and _is_gen == 0.
        """
        cdef Polynomial_integer_dense_flint x = PY_NEW(Polynomial_integer_dense_flint)
        x._parent = self._parent
        x._is_gen = 0
        return x

    cpdef Polynomial _new_constant_poly(self, a, Parent P):
        r"""
        Quickly creates a new constant polynomial with value a in parent P

        ASSUMPTION:

        The given value has to be in the base ring of P. This assumption is not
        verified.

        EXAMPLE::

            sage: R.<x> = ZZ[]
            sage: x._new_constant_poly(2,R)
            2

        """
        cdef Polynomial_integer_dense_flint x = PY_NEW(Polynomial_integer_dense_flint)
        x._parent = P
        x._is_gen = 0
        if not PY_TYPE_CHECK(a, Integer):
            a = ZZ(a)
        fmpz_poly_set_coeff_mpz(x.__poly, 0, (<Integer>a).value)
        return x

    def __init__(self, parent, x=None, check=True, is_gen=False,
            construct=False):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: x
            x

        Construct from list::

            sage: R([])
            0
            sage: R([1, -2, 3])
            3*x^2 - 2*x + 1

        Coercions from other rings are attempted automatically::

            sage: R([1, -6/3, 3])
            3*x^2 - 2*x + 1
            sage: R([1, 5/2, 2])
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer

        Construct from constant::

            sage: R(3)
            3

        Coercion from PARI polynomial::

            sage: f = R([-1, 2, 5]); f
            5*x^2 + 2*x - 1
            sage: type(f)
            <type 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
            sage: type(pari(f))
            <type 'sage.libs.pari.gen.gen'>
            sage: type(R(pari(f)))
            <type 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
            sage: R(pari(f))
            5*x^2 + 2*x - 1

        Coercion from NTL polynomial::

            sage: f = ntl.ZZX([1, 2, 3])
            sage: print R(f)
            3*x^2 + 2*x + 1

        Coercion from dictionary::

            sage: f = R({2: -4, 3: 47}); f
            47*x^3 - 4*x^2

        Coercion from fraction field element with trivial denominator::

            sage: f = (x^3 - 1) / (x - 1)
            sage: type(f)
            <type 'sage.rings.fraction_field_element.FractionFieldElement'>
            sage: g = R(f); g
            x^2 + x + 1

            sage: ZZ['x']({2^3: 1})
            x^8

        """
        Polynomial.__init__(self, parent, is_gen=is_gen)

        cdef Py_ssize_t degree
        cdef Py_ssize_t i

        if x is None:
            return         # leave initialized to 0 polynomial.

        if isinstance(x, Polynomial):
            if x.parent() is self.parent():
                sig_on()
                fmpz_poly_set(self.__poly, \
                        (<Polynomial_integer_dense_flint>x).__poly)
                sig_off()
                return
            else:
                # coerce coefficients into Sage integers
                x = [Integer(a) for a in x.list()]
                check = False

        elif isinstance(x, dict):
            x = x.items()
            degree = 0
            # find max degree to allocate only once
            for ii, a in x:
                # mpoly dict style has tuple keys
                i = ii[0] if PY_TYPE_CHECK_EXACT(ii, tuple) else ii
                if i < 0:
                    raise ValueError, "Negative monomial degrees not allowed: %s" % i
                elif i > degree:
                    degree = i
            try:
                sig_on()
                fmpz_poly_realloc(self.__poly, degree + 1)
                sig_off()
            except RuntimeError:
                raise OverflowError, "Cannot allocate memory!"
            # now fill them in
            for ii, a in x:
                i = ii[0] if PY_TYPE_CHECK_EXACT(ii, tuple) else ii
                if PY_TYPE_CHECK_EXACT(a, int):
                    sig_on()
                    fmpz_poly_set_coeff_si(self.__poly, i, a)
                    sig_off()
                else:
                    if not PY_TYPE_CHECK(a, Integer):
                        a = ZZ(a)
                    sig_on()
                    fmpz_poly_set_coeff_mpz(self.__poly, i, (<Integer>a).value)
                    sig_off()
            return

        elif isinstance(x, pari_gen):
            x = [Integer(w) for w in x.list()]
            check = False

        elif isinstance(x, ntl_ZZX):    # coercion from ntl.pyx object
            fmpz_poly_set_ZZX(self.__poly, (<ntl_ZZX>x).x)
            return

        elif isinstance(x, FractionFieldElement) and \
                 isinstance(x.numerator(), Polynomial_integer_dense_flint):
            if x.denominator() == 1:
                # fraction of the form f(x)/1
                sig_on()
                fmpz_poly_set(self.__poly,
                        (<Polynomial_integer_dense_flint>x.numerator()).__poly)
                sig_off()
                return

        elif not isinstance(x, list):
            x = [x]   # constant polynomials

        sig_on()
        fmpz_poly_realloc(self.__poly, len(x))
        sig_off()
        for i from 0 <= i < len(x):
            a = x[i]
            if PY_TYPE_CHECK_EXACT(a, int):
                sig_on()
                fmpz_poly_set_coeff_si(self.__poly, i, a)
                sig_off()
            else:
                if not PY_TYPE_CHECK(a, Integer):
                    a = ZZ(a)
                sig_on()
                fmpz_poly_set_coeff_mpz(self.__poly, i, (<Integer>a).value)
                sig_off()

    def __call__(self, *x, **kwds):
        """
        Calls this polynomial with the given parameters, which can be
        interpreted as polynomial composition or evaluation by this
        method.

        If the argument is not simply an integer (``int``, ``long`` or
        ``Integer``) or a polynomial (of the same type as ``self``),
        the call is passed on to the generic implementation in the
        ``Polynomial`` class.

        EXAMPLES:

        The first example illustrates polynomial composition::

            sage: R.<t> = ZZ[]
            sage: f = t^2 - 1
            sage: g = t + 1
            sage: f(g)          # indirect doctest
            t^2 + 2*t

        Now we illustrate how a polynomial can be evaluated at an
        integer::

            sage: f(2)          # indirect doctest
            3
        """
        cdef Polynomial_integer_dense_flint f
        cdef Integer a, z

        cdef unsigned long limbs
        cdef fmpz_t a_fmpz
        cdef fmpz_t z_fmpz

        if len(x) == 1:
            x0 = x[0]
            if isinstance(x, Polynomial_integer_dense_flint):
                f = self._new()
                sig_on()
                fmpz_poly_compose(f.__poly, self.__poly, \
                    (<Polynomial_integer_dense_flint> x0).__poly)
                sig_off()
                return f
            if isinstance(x0, (int, long)):
                x0 = Integer(x0)
            if isinstance(x0, Integer):
                a = <Integer> x0

                if fmpz_poly_length(self.__poly) == 0:
                    return ZZ.zero()
                if mpz_sgn(a.value) == 0:
                    return self[0]

                z = PY_NEW(Integer)

                sig_on()
                fmpz_init(a_fmpz)
                fmpz_init(z_fmpz)
                fmpz_set_mpz(a_fmpz, a.value)
                fmpz_poly_evaluate_fmpz(z_fmpz, self.__poly, a_fmpz)
                fmpz_get_mpz(z.value, z_fmpz)
                fmpz_clear(a_fmpz)
                fmpz_clear(z_fmpz)
                sig_off()

                return z

        return Polynomial.__call__(self, *x, **kwds)

    cpdef Integer content(self):
        r"""
        Return the greatest common divisor of the coefficients of this
        polynomial.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: (2*x^2 - 4*x^4 + 14*x^7).content()
            2
            sage: x.content()
            1
            sage: R(1).content()
            1
            sage: R(0).content()
            0

        TESTS::

            sage: t = x^2+x+1
            sage: t.content()
            1
            sage: (123456789123456789123456789123456789123456789*t).content()
            123456789123456789123456789123456789123456789
        """
        cdef fmpz_t c
        fmpz_init(c)
        fmpz_poly_content(c, self.__poly)
        cdef Integer z = PY_NEW(Integer)
        fmpz_get_mpz(z.value, c)
        fmpz_clear(c)
        return z

    def __reduce__(self):
        r"""
        Used for pickling.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: loads(dumps(x)) == x
            True
            sage: f = 2*x + 3
            sage: loads(dumps(f)) == f
            True
        """
        return Polynomial_integer_dense_flint, \
               (self.parent(), self.list(), False, self.is_gen())

    def __getitem__(self, n):
        r"""
        Returns coefficient of x^n, or zero if n is negative.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 2*x^2 - 3
            sage: f[0]
            -3
            sage: f[1]
            0
            sage: f[2]
            2
            sage: f[3]
            0
            sage: f[-1]
            0
            sage: f = 1 + x + 2*x^2 + 3*x^3 + 4*x^4 + 5*x^5
            sage: f[2:4]
            3*x^3 + 2*x^2
            sage: f[-2:4]
            3*x^3 + 2*x^2 + x + 1
            sage: f[4:100]
            5*x^5 + 4*x^4
        """
        cdef long k
        cdef Integer z = PY_NEW(Integer)
        if isinstance(n, slice):
            start = max(0, n.start)
            stop = n.stop
            if stop is None or stop > self.degree()+1:
                stop = self.degree() + 1
            v = [self[k] for k from start <= k < stop]
            P = self.parent()
            return P([0] * int(start) + v)
        else:
            if n < 0 or n > fmpz_poly_degree(self.__poly):
                return z
            else:
                fmpz_poly_get_coeff_mpz(z.value, self.__poly, n)
                return z

    def _repr(self, name=None, bint latex=False):
        """
        Return string representation of this polynomial.

        EXAMPLES::

            sage: R.<x> = ZZ['x']
            sage: (-x+1)^5
            -x^5 + 5*x^4 - 10*x^3 + 10*x^2 - 5*x + 1
            sage: ((-x+1)^5)._repr()
            '-x^5 + 5*x^4 - 10*x^3 + 10*x^2 - 5*x + 1'
            sage: ((-x+1)^5)._repr(name='y')
            '-y^5 + 5*y^4 - 10*y^3 + 10*y^2 - 5*y + 1'
        """
        if name is None:
            name = self.parent().variable_name()
        cdef long i
        cdef mpz_t coef
        mpz_init(coef)
        all = []
        for i from fmpz_poly_degree(self.__poly) >= i >= 0:
            fmpz_poly_get_coeff_mpz(coef, self.__poly, i)
            sign = mpz_sgn(coef)
            if sign:
                if sign > 0:
                    sign_str = '+'
                    coeff_str = mpz_to_str(coef)
                else:
                    sign_str = '-'
                    coeff_str = mpz_to_str(coef)[1:]
                if i > 0:
                    if coeff_str == '1':
                        coeff_str = ''
                    elif not latex:
                        coeff_str = coeff_str + '*'
                if i > 1:
                    if latex:
                        PyList_Append(all, " %s %s%s^{%s}" % (sign_str,
                            coeff_str, name, i))
                    else:
                        PyList_Append(all, " %s %s%s^%s" % (sign_str,
                            coeff_str, name, i))
                elif i == 1:
                    PyList_Append(all, " %s %s%s" % (sign_str, coeff_str, name))
                else:
                    PyList_Append(all, " %s %s" % (sign_str, coeff_str))
        mpz_clear(coef)
        if len(all) == 0:
            return '0'
        leading = all[0]
        if leading[1] == '+':
            all[0] = leading[3:]
        else:
            all[0] = '-' + leading[3:]
        return ''.join(all)

    def _latex_(self, name=None):
        """
        Return the latex representation of this polynomial.

        EXAMPLES::

            sage: R.<t> = ZZ['t']
            sage: latex(t^10-t^2-5*t+1)
            t^{10} - t^{2} - 5t + 1
            sage: cyclotomic_polynomial(10^5)._latex_()
            'x^{40000} - x^{30000} + x^{20000} - x^{10000} + 1'
            sage: cyclotomic_polynomial(10^5)._latex_(name='y')
            'y^{40000} - y^{30000} + y^{20000} - y^{10000} + 1'
        """
        if name is None:
            name = self.parent().latex_variable_names()[0]
        return self._repr(name=name, latex=True)

    cpdef ModuleElement _add_(self, ModuleElement right):
        r"""
        Returns self plus right.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 2*x + 1
            sage: g = -3*x^2 + 6
            sage: f + g
            -3*x^2 + 2*x + 7
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_add(x.__poly, self.__poly,
                (<Polynomial_integer_dense_flint>right).__poly)
        sig_off()
        return x


    cpdef ModuleElement _sub_(self, ModuleElement right):
        r"""
        Return self minus right.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 2*x + 1
            sage: g = -3*x^2 + 6
            sage: f - g
            3*x^2 + 2*x - 5
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_sub(x.__poly, self.__poly,
                (<Polynomial_integer_dense_flint>right).__poly)
        sig_off()
        return x


    cpdef ModuleElement _neg_(self):
        r"""
        Returns negative of self.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 2*x - 1
            sage: -f
            -2*x + 1
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_neg(x.__poly, self.__poly)
        sig_off()
        return x

    @coerce_binop
    def quo_rem(self, Polynomial_integer_dense_flint right):
        r"""
        Attempts to divide self by right, and return a quotient and remainder.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = R(range(10)); g = R([-1, 0, 1])
            sage: q, r = f.quo_rem(g)
            sage: q, r
            (9*x^7 + 8*x^6 + 16*x^5 + 14*x^4 + 21*x^3 + 18*x^2 + 24*x + 20, 25*x + 20)
            sage: q*g + r == f
            True

            sage: f = x^2
            sage: f.quo_rem(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero polynomial

            sage: f = (x^2 + 3) * (2*x - 1)
            sage: f.quo_rem(2*x - 1)
            (x^2 + 3, 0)

            sage: f = x^2
            sage: f.quo_rem(2*x - 1)
            (0, x^2)

        TESTS::

            sage: z = R(0)
            sage: z.quo_rem(1)
            (0, 0)
            sage: z.quo_rem(x)
            (0, 0)
            sage: z.quo_rem(2*x)
            (0, 0)

        Ticket #383, make sure things get coerced correctly::

            sage: f = x+1; parent(f)
            Univariate Polynomial Ring in x over Integer Ring
            sage: g = x/2; parent(g)
            Univariate Polynomial Ring in x over Rational Field
            sage: f.quo_rem(g)
            (2, 1)
            sage: g.quo_rem(f)
            (1/2, -1/2)
            sage: parent(f.quo_rem(g)[0])
            Univariate Polynomial Ring in x over Rational Field
            sage: f.quo_rem(3)
            (0, x + 1)
            sage: (5*x+7).quo_rem(3)
            (x + 2, 2*x + 1)
        """
        if right.is_zero():
            raise ZeroDivisionError, "division by zero polynomial"

        if self.is_zero():
            return self, self

        cdef Polynomial_integer_dense_flint qq = self._new()
        cdef Polynomial_integer_dense_flint rr = self._new()

        sig_on()
        fmpz_poly_divrem(qq.__poly, rr.__poly, self.__poly, right.__poly)
        sig_off()
        return qq, rr

    cpdef bint is_zero(self):
        """
        Returns True if self is equal to zero.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: R(0).is_zero()
            True
            sage: R(1).is_zero()
            False
            sage: x.is_zero()
            False
        """
        return (fmpz_poly_degree(self.__poly) == -1)

    def __nonzero__(self):
        """
        Check if self is not zero.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: R(0).__nonzero__()
            False
            sage: R(1).__nonzero__()
            True
            sage: x.__nonzero__()
            True
        """
        return not (fmpz_poly_degree(self.__poly) == -1)

    @coerce_binop
    def gcd(self, right):
        r"""
        Return the GCD of self and right.  The leading
        coefficient need not be 1.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = (6*x + 47)*(7*x^2 - 2*x + 38)
            sage: g = (6*x + 47)*(3*x^3 + 2*x + 1)
            sage: f.gcd(g)
            6*x + 47
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_gcd(x.__poly, self.__poly,
                (<Polynomial_integer_dense_flint>right).__poly)
        sig_off()
        return x


    @coerce_binop
    def lcm(self, right):
        """
        Return the LCM of self and right.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = (6*x + 47)*(7*x^2 - 2*x + 38)
            sage: g = (6*x + 47)*(3*x^3 + 2*x + 1)
            sage: h = f.lcm(g); h
            126*x^6 + 951*x^5 + 486*x^4 + 6034*x^3 + 585*x^2 + 3706*x + 1786
            sage: h == (6*x + 47)*(7*x^2 - 2*x + 38)*(3*x^3 + 2*x + 1)
            True
        """
        g = self.gcd(right)
        return (self//g)*right


    @coerce_binop
    def xgcd(self, right):
        """
        This function can't in general return ``(g,s,t)`` as above,
        since they need not exist.  Instead, over the integers, we
        first multiply `g` by a divisor of the resultant of `a/g` and
        `b/g`, up to sign, and return ``g, u, v`` such that
        ``g = s*self + s*right``.  But note that this `g` may be a
        multiple of the gcd.

        If ``self`` and ``right`` are coprime as polynomials over the
        rationals, then ``g`` is guaranteed to be the resultant of
        self and right, as a constant polynomial.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(ZZ)
            sage: F = (x^2 + 2)*x^3; G = (x^2+2)*(x-3)
            sage: g, u, v = F.xgcd(G)
            sage: g, u, v
            (27*x^2 + 54, 1, -x^2 - 3*x - 9)
            sage: u*F + v*G
            27*x^2 + 54
            sage: x.xgcd(P(0))
            (x, 1, 0)
            sage: f = P(0)
            sage: f.xgcd(x)
            (x, 0, 1)
            sage: F = (x-3)^3; G = (x-15)^2
            sage: g, u, v = F.xgcd(G)
            sage: g, u, v
            (2985984, -432*x + 8208, 432*x^2 + 864*x + 14256)
            sage: u*F + v*G
            2985984
        """
        cdef Polynomial_integer_dense_flint ss = self._new()
        cdef Polynomial_integer_dense_flint tt = self._new()
        cdef fmpz_t r
        fmpz_init(r)

        sig_on()
        fmpz_poly_xgcd(r, ss.__poly, tt.__poly, self.__poly,
                (<Polynomial_integer_dense_flint>right).__poly)
        sig_off()
        cdef Integer rr = PY_NEW(Integer)
        fmpz_get_mpz(rr.value, r)
        fmpz_clear(r)

        if rr.is_zero():
            f = self.base_extend(QQ)
            g, u, v = f.xgcd(right.base_extend(QQ))
            d = lcm([g.denominator(), u.denominator(), v.denominator()])
            R = self.parent()
            return R(d*g), R(d*u), R(d*v)
        else:
            return self._parent(rr), ss, tt


    cpdef RingElement _mul_(self, RingElement right):
        r"""
        Returns self multiplied by right.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: (x - 2)*(x^2 - 8*x + 16)
            x^3 - 10*x^2 + 32*x - 32
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_mul(x.__poly, self.__poly,
                (<Polynomial_integer_dense_flint>right).__poly)
        sig_off()
        return x


    cpdef ModuleElement _lmul_(self, RingElement right):
        r"""
        Returns self multiplied by right, where right is a scalar (integer).

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: x*3
            3*x
            sage: (2*x^2 + 4)*3
            6*x^2 + 12
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_scalar_mul_mpz(x.__poly, self.__poly, (<Integer>right).value)
        sig_off()
        return x


    cpdef ModuleElement _rmul_(self, RingElement right):
        r"""
        Returns self multiplied by right, where right is a scalar (integer).

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: 3*x
            3*x
            sage: 3*(2*x^2 + 4)
            6*x^2 + 12
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_scalar_mul_mpz(x.__poly, self.__poly, (<Integer>right).value)
        sig_off()
        return x

    def __pow__(Polynomial_integer_dense_flint self, int exp, ignored):
        """

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: r = 2*x + 2
            sage: r^0
            1
            sage: r^2
            4*x^2 + 8*x + 4
            sage: r^-2
            1/(4*x^2 + 8*x + 4)

            sage: x^(2^20)
            x^1048576

        TESTS::

            sage: z = R(0)
            sage: z^0
            1
            sage: z^1
            0
            sage: z^-1
            Traceback (most recent call last):
            ...
            ZeroDivisionError: negative exponent in power of zero
        """
        cdef long nn = exp
        cdef Polynomial_integer_dense_flint res = self._new()
        if self.is_zero():
            if exp == 0:
                fmpz_poly_set_coeff_si(res.__poly, 0, 1)
                return res
            elif exp < 0:
                raise ZeroDivisionError, "negative exponent in power of zero"
            else:
                return res
        if exp < 0:
            sig_on()
            fmpz_poly_pow(res.__poly, self.__poly, -nn)
            sig_off()
            return ~res
        else:
            if self is self._parent.gen():
                sig_on()
                fmpz_poly_set_coeff_ui(res.__poly, exp, 1)
                sig_off()
            else:
                sig_on()
                fmpz_poly_pow(res.__poly, self.__poly, nn)
                sig_off()
            return res

    def __floordiv__(Polynomial_integer_dense_flint self, right):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: (x^2+1)//x
            x
            sage: (5*x^2+1)//(2*x)
            2*x

        Divide by a scalar.

        ::

            sage: (5*x^3 + 5*x + 10)//5
            x^3 + x + 2

        TESTS::

            sage: x//0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero

            sage: (x^2 + 13*x + 169) // 13
            x + 13
        """
        cdef Polynomial_integer_dense_flint res = self._new()
        cdef Polynomial
        cdef long t
        if right == 0:
            raise ZeroDivisionError, "division by zero"
        if not PY_TYPE_CHECK(right, Polynomial_integer_dense_flint):
            if right in ZZ:
                sig_on()
                fmpz_poly_scalar_fdiv_mpz(res.__poly, self.__poly,
                        (<Integer>ZZ(right)).value)
                sig_off()
                return res
        if self._parent is not right.parent():
            right = self._parent(right)
        sig_on()
        fmpz_poly_div(res.__poly, self.__poly,
                (<Polynomial_integer_dense_flint>right).__poly)
        sig_off()
        return res

    cpdef _unsafe_mutate(self, long n, value):
        r"""
        Sets coefficient of `x^n` to value.

        This is very unsafe, because Sage polynomials are supposed to be
        immutable. (Shhhh don't tell anyone!)

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 2*x^2 + 3; f
            2*x^2 + 3
            sage: f._unsafe_mutate(1, 42); f
            2*x^2 + 42*x + 3

            sage: f._unsafe_mutate(1, int(5)); f
            2*x^2 + 5*x + 3
            sage: f._unsafe_mutate(1, Zmod(15)(7)); f
            2*x^2 + 7*x + 3
        """
        if n < 0:
            raise IndexError, "n must be >= 0"
        if PY_TYPE_CHECK(value, int):
            sig_on()
            fmpz_poly_set_coeff_si(self.__poly, n, value)
            sig_off()
        elif PY_TYPE_CHECK(value, Integer):
            sig_on()
            fmpz_poly_set_coeff_mpz(self.__poly, n, (<Integer>value).value)
            sig_off()
        else:
            value = Integer(value)
            sig_on()
            fmpz_poly_set_coeff_mpz(self.__poly, n, (<Integer>value).value)
            sig_off()

    def real_root_intervals(self):
        """
        Returns isolating intervals for the real roots of this
        polynomial.

        EXAMPLE:
        We compute the roots of the characteristic polynomial of some
        Salem numbers::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 1 - x^2 - x^3 - x^4 + x^6
            sage: f.real_root_intervals()
            [((1/2, 3/4), 1), ((1, 3/2), 1)]
        """

        from sage.rings.polynomial.real_roots import real_roots

        return real_roots(self)

##     def __copy__(self):
##         f = Polynomial_integer_dense(self.parent())
##         f.__poly = self.__poly.copy()
##         return f


    def degree(self, gen=None):
        """
        Return the degree of this polynomial.  The zero polynomial
        has degree -1.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: x.degree()
            1
            sage: (x^2).degree()
            2
            sage: R(1).degree()
            0
            sage: R(0).degree()
            -1
        """
        return fmpz_poly_degree(self.__poly)

    def pseudo_divrem(self, B):
        """
        Write ``A = self``.  This function computes polynomials `Q` and `R`
        and an integer `d` such that

        .. math::

             \mathop{\mathrm{lead}}(B)^d A = B Q + R

        where R has degree less than that of B.

        INPUT:

        - ``B`` -- a polynomial over `\ZZ`

        OUTPUT:

        - ``Q``, ``R`` -- polynomials
        - ``d`` -- nonnegative integer

        EXAMPLES::

            sage: R.<x> = ZZ['x']
            sage: A = R(range(10)); B = 3*R([-1, 0, 1])
            sage: Q, R, d = A.pseudo_divrem(B)
            sage: Q, R, d
            (9*x^7 + 8*x^6 + 16*x^5 + 14*x^4 + 21*x^3 + 18*x^2 + 24*x + 20, 75*x + 60, 1)
            sage: B.leading_coefficient()^d * A == B*Q + R
            True
        """
        cdef Polynomial_integer_dense_flint Q = self._new(), R = self._new(), _B = B
        cdef unsigned long d
        fmpz_poly_pseudo_divrem(Q.__poly, R.__poly, &d, self.__poly, _B.__poly)
        return Q, R, Integer(d)

    def discriminant(self, proof=True):
        r"""
        Return the discriminant of self, which is by definition

        .. math::

            (-1)^{m(m-1)/2} \mathop{\mathrm{resultant}}(a, a')/\mathop{\mathrm{lc}}(a),

        where `m = \mathop{\mathrm{deg}}(a)`, and `\mathop{\mathrm{lc}}(a)` is
        the leading coefficient of a. If ``proof`` is False (the default is
        True), then this function may use a randomized strategy that errors
        with probability no more than `2^{-80}`.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = 3*x^3 + 2*x + 1
            sage: f.discriminant()
            -339
            sage: f.discriminant(proof=False)
            -339
        """
        cdef ZZX_c ntl_poly
        cdef ZZ_c* temp
        cdef Integer x
        fmpz_poly_get_ZZX(ntl_poly, self.__poly)

        temp = ZZX_discriminant(&ntl_poly, proof)
        x = PY_NEW(Integer)
        ZZ_to_mpz(&x.value, temp)
        ZZ_delete(temp)

        return x


    def _pari_(self, variable=None):
        """
        EXAMPLES::

            sage: t = PolynomialRing(ZZ,"t").gen()
            sage: f = t^3 + 3*t - 17
            sage: pari(f)
            t^3 + 3*t - 17
            sage: f._pari_(variable='y')
            y^3 + 3*y - 17
        """
        if variable is None:
            variable = self.parent().variable_name()
        return pari(self.list()).Polrev(variable)


    def squarefree_decomposition(Polynomial_integer_dense_flint self):
        """
        Return the square-free decomposition of self.  This is
        a partial factorization of self into square-free, relatively
        prime polynomials.

        This is a wrapper for the NTL function SquareFreeDecomp.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: p = (x-1)^2 * (x-2)^2 * (x-3)^3 * (x-4)
            sage: p.squarefree_decomposition()
            (x - 4) * (x^2 - 3*x + 2)^2 * (x - 3)^3
            sage: p = 37 * (x-1)^2 * (x-2)^2 * (x-3)^3 * (x-4)
            sage: p.squarefree_decomposition()
            (37) * (x - 4) * (x^2 - 3*x + 2)^2 * (x - 3)^3
        """
        cdef ZZX_c** v
        cdef long* e
        cdef long i, n
        cdef fmpz_poly_t ppart
        cdef ZZX_c ntl_poly
        cdef Integer z
        cdef Polynomial_integer_dense_flint fac

        z = self.content()
        if not z.is_one():
            fmpz_poly_init(ppart)

            fmpz_poly_primitive_part(ppart, self.__poly)

            fmpz_poly_get_ZZX(ntl_poly, ppart)
            fmpz_poly_clear(ppart)
        else:
            fmpz_poly_get_ZZX(ntl_poly, self.__poly)

        ZZX_squarefree_decomposition(&v, &e, &n, &ntl_poly)

        F = []
        for i from 0 <= i < n:
            fac = self._new()
            fmpz_poly_set_ZZX(fac.__poly, v[i][0])
            F.append( (fac,e[i]) )
            ZZX_delete(v[i])
        sage_free(v)
        sage_free(e)

        return Factorization(F, unit=z, sort=False)

    def _factor_pari(self):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = (x^2-2)*(x^5-3)^2
            sage: f._factor_pari()
            (x^2 - 2) * (x^5 - 3)^2
            sage: (1234567898765432123456789876543212345678987*f)._factor_pari()
            1234567898765432123456789876543212345678987 * (x^2 - 2) * (x^5 - 3)^2
        """
        return Polynomial.factor(self) # uses pari for integers over ZZ

    def _factor_ntl(self):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = (x^2-2)*(x^5-3)^2
            sage: f._factor_ntl()
            (x^2 - 2) * (x^5 - 3)^2
            sage: (12345678987654321234567898765432123456789876*f)._factor_ntl()
            12345678987654321234567898765432123456789876 * (x^2 - 2) * (x^5 - 3)^2
        """
        cdef Polynomial_integer_dense_flint fac_py
        cdef fmpz_t tcontent
        cdef ZZX_c ntl_poly
        cdef ZZ_c content
        cdef vec_pair_ZZX_long_c factors
        cdef long i
        cdef int sig_me = fmpz_poly_degree(self.__poly)

        fmpz_poly_get_ZZX(ntl_poly, self.__poly)

        if sig_me > 10:
            sig_on()
        ZZX_factor(content, factors, ntl_poly, 0, 0)
        if sig_me > 10:
            sig_off()

        results = []
        unit = None

        if ZZ_sign(content) < 0:
            unit = Integer(-1)
            ZZ_abs(content, content)

        if not ZZ_IsOne(content):
            fac_py = self._new()
            fmpz_init(tcontent)
            fmpz_set_ZZ(tcontent, content)
            fmpz_poly_set_coeff_fmpz(fac_py.__poly, 0, tcontent)
            results.append( (fac_py,1) )
            fmpz_clear(tcontent)

        for i from 0 <= i < factors.length():
            fac_py = self._new()
            fmpz_poly_set_ZZX(fac_py.__poly, factors.RawGet(i).a)
            results.append( (fac_py,factors.RawGet(i).b) )
        return Factorization(results, unit = unit)

    def factor(self):
        """
        This function overrides the generic polynomial factorization to
        make a somewhat intelligent decision to use Pari or NTL based on
        some benchmarking.

        Note: This function factors the content of the polynomial,
        which can take very long if it's a really big integer.  If you
        do not need the content factored, divide it out of your
        polynomial before calling this function.

        EXAMPLES::

            sage: R.<x>=ZZ[]
            sage: f=x^4-1
            sage: f.factor()
            (x - 1) * (x + 1) * (x^2 + 1)
            sage: f=1-x
            sage: f.factor()
            (-1) * (x - 1)
            sage: f.factor().unit()
            -1
            sage: f = -30*x; f.factor()
            (-1) * 2 * 3 * 5 * x
        """
        cdef int i
        cdef long deg = fmpz_poly_degree(self.__poly)
        # it appears that pari has a window from about degrees 30 and 300
        # in which it beats NTL.
        c = self.content()
        g = self//c
        if deg < 30 or deg > 300:
            return c.factor()*g._factor_ntl()
        else:
            return c.factor()*g._factor_pari()

    def factor_mod(self, p):
        """
        Return the factorization of self modulo the prime `p`.

        INPUT:

        - ``p`` -- prime

        OUTPUT:

        factorization of self reduced modulo p.

        EXAMPLES::

            sage: R.<x> = ZZ['x']
            sage: f = -3*x*(x-2)*(x-9) + x
            sage: f.factor_mod(3)
            x
            sage: f = -3*x*(x-2)*(x-9)
            sage: f.factor_mod(3)
            Traceback (most recent call last):
            ...
            ValueError: factorization of 0 not defined

            sage: f = 2*x*(x-2)*(x-9)
            sage: f.factor_mod(7)
            (2) * x * (x + 5)^2
        """
        from sage.rings.finite_rings.constructor import FiniteField
        p = Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"
        if all([c%p==0 for c in self.coefficients()]):
            raise ValueError, "factorization of 0 not defined"
        f = self._pari_()
        G = f.factormod(p)
        k = FiniteField(p)
        R = k[self.parent().variable_name()]
        return R(1)._factor_pari_helper(G, unit=R(self).leading_coefficient())


    def factor_padic(self, p, prec=10):
        """
        Return `p`-adic factorization of self to given precision.

        INPUT:

        - ``p`` -- prime

        - ``prec`` -- integer; the precision

        OUTPUT:

        - factorization of ``self`` over the completion at `p`.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = x^2 + 1
            sage: f.factor_padic(5, 4)
            ((1 + O(5^4))*x + (2 + 5 + 2*5^2 + 5^3 + O(5^4))) * ((1 + O(5^4))*x + (3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4)))

        A more difficult example::

            sage: f = 100 * (5*x + 1)^2 * (x + 5)^2
            sage: f.factor_padic(5, 10)
            (4 + O(5^10)) * ((5 + O(5^11)))^2 * ((1 + O(5^10))*x + (5 + O(5^10)))^2 * ((5 + O(5^10))*x + (1 + O(5^10)))^2

        """
        from sage.rings.padics.factory import Zp

        p = Integer(p)
        prec = Integer(prec)

        # Parent field for coefficients and polynomial
        K = Zp(p, prec, type='capped-rel')
        R = K[self.parent().variable_name()]

        # Factor the *exact* polynomial using factorpadic()
        G = self._pari_with_name().factorpadic(p, prec)

        from sage.rings.polynomial.padics.polynomial_padic import _pari_padic_factorization_to_sage
        return _pari_padic_factorization_to_sage(G, R, self.leading_coefficient())

    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of self.

        EXAMPLES::

            sage: x = PolynomialRing(ZZ,'x').0
            sage: f = x^3 + 3*x - 17
            sage: f.list()
            [-17, 3, 0, 1]
            sage: f = PolynomialRing(ZZ,'x')(0)
            sage: f.list()
            []
        """
        return [self[i] for i in range(self.degree()+1)]


    @coerce_binop
    def resultant(self, other, proof=True):
        """
        Returns the resultant of self and other, which must lie in the same
        polynomial ring.

        If ``proof = False`` (the default is ``proof=True``), then this function may
        use a randomized strategy that errors with probability no more than
        `2^{-80}`.

        INPUT:

        - other -- a polynomial

        OUTPUT:

        an element of the base ring of the polynomial ring

        EXAMPLES::

            sage: x = PolynomialRing(ZZ,'x').0
            sage: f = x^3 + x + 1;  g = x^3 - x - 1
            sage: r = f.resultant(g); r
            -8
            sage: r.parent() is ZZ
            True
        """
        if not isinstance(other, Polynomial_integer_dense_flint):
            other = self.parent()(other)
        elif self.parent() is not other.parent():
            raise TypeError

        cdef fmpz_t res
        fmpz_init(res)
        cdef Integer x = PY_NEW(Integer)

        sig_on()
        fmpz_poly_resultant(res, self.__poly,
                (<Polynomial_integer_dense_flint>other).__poly)
        sig_off()
        fmpz_get_mpz(x.value, res)
        fmpz_clear(res)
        return x

    def reverse(self, degree=None):
        """
        Return a polynomial with the coefficients of this polynomial reversed.

        If an optional degree argument is given the coefficient list will be
        truncated or zero padded as necessary and the reverse polynomial will
        have the specified degree.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: p = R([1,2,3,4]); p
            4*x^3 + 3*x^2 + 2*x + 1
            sage: p.reverse()
            x^3 + 2*x^2 + 3*x + 4
            sage: p.reverse(degree=6)
            x^6 + 2*x^5 + 3*x^4 + 4*x^3
            sage: p.reverse(degree=2)
            x^2 + 2*x + 3

        TESTS::

            sage: p.reverse(degree=1.5r)
            Traceback (most recent call last):
            ...
            ValueError: degree argument must be a non-negative integer, got 1.5
        """
        cdef Polynomial_integer_dense_flint res = self._new()
        cdef unsigned long d
        if degree:
            d = degree
            if d != degree:
                raise ValueError, "degree argument must be a non-negative integer, got %s"%(degree)
            # FLINT expects length
            fmpz_poly_reverse(res.__poly, self.__poly, d+1)
        else:
            fmpz_poly_reverse(res.__poly, self.__poly,
                    fmpz_poly_length(self.__poly))
        return res
