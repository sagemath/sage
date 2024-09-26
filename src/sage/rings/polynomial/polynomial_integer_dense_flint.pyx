# distutils: libraries = NTL_LIBRARIES gmp
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
r"""
Dense univariate polynomials over `\ZZ`, implemented using FLINT

AUTHORS:

- David Harvey: rewrote to talk to NTL directly, instead of via ntl.pyx
  (2007-09); a lot of this was based on Joel Mohler's recent rewrite of the NTL
  wrapper
- David Harvey: split off from polynomial_element_generic.py (2007-09)
- Burcin Erocal: rewrote to use FLINT (2008-06-16)

TESTS:

We check that the buggy gcd is fixed (see :issue:`17816`)::

    sage: R.<q> = ZZ[]
    sage: X = 3*q^12 - 8*q^11 - 24*q^10 - 48*q^9 - 84*q^8 - 92*q^7 - 92*q^6 - 70*q^5 - 50*q^4 - 27*q^3 - 13*q^2 - 4*q - 1
    sage: Y = q^13 - 2*q^12 + 2*q^10 - q^9
    sage: gcd(X,Y)
    1
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2008-2010 Burcin Erocal <burcin@erocal.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.memory cimport sig_free
from cysignals.signals cimport sig_on, sig_off

include "sage/libs/ntl/decl.pxi"

from cpython.long cimport PyLong_AsLong
from sage.libs.gmp.mpz cimport *
from sage.arith.long cimport pyobject_to_long, is_small_python_int

from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.structure.element cimport Element
from sage.structure.element import coerce_binop

from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from sage.libs.pari.all import pari, pari_gen
from sage.structure.factorization import Factorization

from sage.rings.fraction_field_element import FractionFieldElement
from sage.arith.functions import lcm

from sage.libs.flint.arb_fmpz_poly cimport arb_fmpz_poly_evaluate_arb, arb_fmpz_poly_evaluate_acb
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.fmpz_poly_sage cimport *
from sage.libs.flint.types cimport ulong, fmpz_poly_t
from sage.libs.flint.ntl_interface cimport fmpz_set_ZZ, fmpz_poly_set_ZZX, fmpz_poly_get_ZZX
from sage.libs.ntl.ZZX cimport *
from sage.rings.complex_arb cimport ComplexBall
from sage.rings.integer cimport Integer, smallInteger
from sage.rings.real_arb cimport RealBall
from sage.rings.real_mpfr cimport RealNumber
from sage.rings.real_mpfi cimport RealIntervalFieldElement

from sage.rings.polynomial.evaluation_flint cimport fmpz_poly_evaluation_mpfr, fmpz_poly_evaluation_mpfi


cdef class Polynomial_integer_dense_flint(Polynomial):
    r"""
    A dense polynomial over the integers, implemented via FLINT.

    TESTS::

        sage: f = ZZ['x'].random_element()
        sage: from sage.rings.polynomial.polynomial_integer_dense_flint import Polynomial_integer_dense_flint
        sage: isinstance(f, Polynomial_integer_dense_flint)
        True

    .. automethod:: _add_
    .. automethod:: _sub_
    .. automethod:: _lmul_
    .. automethod:: _rmul_
    .. automethod:: _mul_
    .. automethod:: _mul_trunc_
    """

    def __cinit__(self):
        r"""
        This calls the underlying FLINT fmpz_poly constructor
        """
        fmpz_poly_init(self._poly)

    def __dealloc__(self):
        r"""
        Call the underlying FLINT fmpz_poly destructor
        """
        fmpz_poly_clear(self._poly)

    cdef Polynomial_integer_dense_flint _new(self):
        r"""
        Quickly create a new initialized ``Polynomial_integer_dense_flint``
        with the correct parent and ``_is_gen == 0``.
        """
        cdef Polynomial_integer_dense_flint x = Polynomial_integer_dense_flint.__new__(Polynomial_integer_dense_flint)
        x._parent = self._parent
        x._is_gen = 0
        return x

    cpdef Polynomial _new_constant_poly(self, a, Parent P):
        r"""
        Quickly create a new constant polynomial with value ``a`` in parent ``P``.

        ASSUMPTION:

        The given value has to be in the base ring of ``P``. This assumption is
        not verified.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: x._new_constant_poly(2,R)
            2
        """
        cdef Polynomial_integer_dense_flint x = Polynomial_integer_dense_flint.__new__(Polynomial_integer_dense_flint)
        x._parent = P
        x._is_gen = 0
        if not isinstance(a, Integer):
            a = ZZ(a)
        fmpz_poly_set_coeff_mpz(x._poly, 0, (<Integer>a).value)
        return x

    def __init__(self, parent, x=None, check=True, is_gen=False,
            construct=False):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: x
            x

        Construct from list and tuple::

            sage: R([])
            0
            sage: R([1, -2, 3])
            3*x^2 - 2*x + 1

            sage: R(())
            0
            sage: R((1, -2, 3))
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
            <class 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
            sage: type(pari(f))
            <class 'cypari2.gen.Gen'>
            sage: type(R(pari(f)))
            <class 'sage.rings.polynomial.polynomial_integer_dense_flint.Polynomial_integer_dense_flint'>
            sage: R(pari(f))
            5*x^2 + 2*x - 1

        Coercion from NTL polynomial::

            sage: f = ntl.ZZX([1, 2, 3])
            sage: print(R(f))
            3*x^2 + 2*x + 1

        Coercion from dictionary::

            sage: f = R({2: -4, 3: 47}); f
            47*x^3 - 4*x^2

        Coercion from fraction field element with trivial denominator::

            sage: f = (x^3 - 1) / (x - 1)
            sage: type(f)
            <class 'sage.rings.fraction_field_element.FractionFieldElement'>
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
                fmpz_poly_set(self._poly,
                              (<Polynomial_integer_dense_flint>x)._poly)
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
                i = ii[0] if type(ii) is tuple else ii
                if i < 0:
                    raise ValueError("Negative monomial degrees not allowed: %s" % i)
                elif i > degree:
                    degree = i
            try:
                sig_on()
                fmpz_poly_realloc(self._poly, degree + 1)
                sig_off()
            except RuntimeError:
                raise OverflowError("Cannot allocate memory!")
            # now fill them in
            for ii, a in x:
                i = ii[0] if type(ii) is tuple else ii
                if is_small_python_int(a):
                    sig_on()
                    fmpz_poly_set_coeff_si(self._poly, i, a)
                    sig_off()
                else:
                    if not isinstance(a, Integer):
                        a = ZZ(a)
                    sig_on()
                    fmpz_poly_set_coeff_mpz(self._poly, i, (<Integer>a).value)
                    sig_off()
            return

        elif isinstance(x, pari_gen):
            x = [Integer(w) for w in x.list()]
            check = False

        elif isinstance(x, ntl_ZZX):    # coercion from ntl.pyx object
            fmpz_poly_set_ZZX(self._poly, (<ntl_ZZX>x).x)
            return

        elif isinstance(x, FractionFieldElement) and \
                 isinstance(x.numerator(), Polynomial_integer_dense_flint):
            if x.denominator() == 1:
                # fraction of the form f(x)/1
                sig_on()
                fmpz_poly_set(self._poly,
                        (<Polynomial_integer_dense_flint>x.numerator())._poly)
                sig_off()
                return

        elif not isinstance(x, (list, tuple)):
            x = [x]   # constant polynomials

        sig_on()
        fmpz_poly_realloc(self._poly, len(x))
        sig_off()
        for i from 0 <= i < len(x):
            a = x[i]
            if is_small_python_int(a):
                sig_on()
                fmpz_poly_set_coeff_si(self._poly, i, a)
                sig_off()
            else:
                if not isinstance(a, Integer):
                    a = ZZ(a)
                sig_on()
                fmpz_poly_set_coeff_mpz(self._poly, i, (<Integer>a).value)
                sig_off()

    def _eval_mpfr_(self, RealNumber a):
        r"""
        Evaluate this polynomial on the real number element ``a``.

        This method uses Horner's rule and might not be appropriate for
        polynomials of large degree.

        TESTS::

            sage: R.<x> = PolynomialRing(ZZ, implementation='FLINT')
            sage: (x+1)._eval_mpfr_(RR(1.2))
            2.20000000000000
            sage: (x^2)._eval_mpfr_(RR(2.2))
            4.84000000000000
            sage: R.zero()._eval_mpfr_(RR(2.1))
            0.000000000000000
            sage: R.one()._eval_mpfr_(RR(2.1))
            1.00000000000000

            sage: p = x^3 - 2*x^2 + x -1
            sage: p._eval_mpfr_(RR(1.3))
            -0.883000000000000
        """
        cdef RealNumber res = a._new()
        sig_on()
        fmpz_poly_evaluation_mpfr(res.value, self._poly, a.value)
        sig_off()
        return res

    def _eval_mpfi_(self, RealIntervalFieldElement a):
        r"""
        Evaluate this polynomial on the real interval ``a``.

        This method uses Horner's rule and might not be appropriate for
        polynomials of large degree.

        TESTS::

            sage: R.<x> = PolynomialRing(ZZ, implementation='FLINT')
            sage: (x+1)._eval_mpfi_(RIF(1.5))
            2.5000000000000000?
            sage: (x^2)._eval_mpfi_(RIF(1.333,1.334))
            1.78?
            sage: R.zero()._eval_mpfi_(RIF(2.1))
            0
            sage: R.one()._eval_mpfi_(RIF(2.1))
            1

            sage: p = x^3 - x^2 - x - 1
            sage: r = p.roots(RIF, multiplicities=False)[0]                             # needs sage.libs.linbox
            sage: p._eval_mpfi_(r)                                                      # needs sage.libs.linbox
            0.?e-27
        """
        cdef RealIntervalFieldElement res = a._new()
        sig_on()
        fmpz_poly_evaluation_mpfi(res.value, self._poly, a.value)
        sig_off()
        return res

    def __call__(self, *x, **kwds):
        """
        Call this polynomial with the given parameters, which can be
        interpreted as polynomial composition or evaluation by this
        method.

        If the argument is not simply an integer (``int`` or
        ``Integer``) a real number (``RealNumber``) a real interval
        (``RealIntervalFieldElement``) or a polynomial (of the same type as
        ``self``), the call is passed on to the generic implementation in the
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

        TESTS:

            sage: t(-sys.maxsize-1r) == t(-sys.maxsize-1)
            True
            sage: (t^2+3)(RealBallField(100)(1/3))
            [3.1111111111111111111111111111...]
            sage: (t^2+3)(ComplexBallField(10)(i))                                      # needs sage.symbolic
            2.00
        """
        cdef Polynomial_integer_dense_flint f
        cdef Integer a, z
        cdef RealBall arb_a, arb_z
        cdef ComplexBall acb_a, acb_z

        cdef unsigned long limbs
        cdef fmpz_t a_fmpz
        cdef fmpz_t z_fmpz

        if len(x) == 1:
            x0 = x[0]
            if isinstance(x, Polynomial_integer_dense_flint):
                f = self._new()
                sig_on()
                fmpz_poly_compose(f._poly, self._poly,
                                  (<Polynomial_integer_dense_flint> x0)._poly)
                sig_off()
                return f
            if is_small_python_int(x0):
                z = Integer.__new__(Integer)
                sig_on()
                fmpz_init(a_fmpz)
                fmpz_init(z_fmpz)
                fmpz_set_si(a_fmpz, PyLong_AsLong(x0))
                fmpz_poly_evaluate_fmpz(z_fmpz, self._poly, a_fmpz)
                fmpz_get_mpz(z.value, z_fmpz)
                fmpz_clear(a_fmpz)
                fmpz_clear(z_fmpz)
                sig_off()
                return z
            if isinstance(x0, int):
                x0 = Integer(x0)
            if isinstance(x0, Integer):
                a = <Integer> x0

                if fmpz_poly_length(self._poly) == 0:
                    return ZZ.zero()
                if mpz_sgn(a.value) == 0:
                    return self[0]

                z = Integer.__new__(Integer)

                sig_on()
                fmpz_init(a_fmpz)
                fmpz_init(z_fmpz)
                fmpz_set_mpz(a_fmpz, a.value)
                fmpz_poly_evaluate_fmpz(z_fmpz, self._poly, a_fmpz)
                fmpz_get_mpz(z.value, z_fmpz)
                fmpz_clear(a_fmpz)
                fmpz_clear(z_fmpz)
                sig_off()

                return z

            if isinstance(x0, RealNumber):
                return self._eval_mpfr_(<RealNumber> x0)
            if isinstance(x0, RealIntervalFieldElement):
                return self._eval_mpfi_(<RealIntervalFieldElement> x0)
            if isinstance(x0, RealBall):
                arb_a = <RealBall> x0
                arb_z = arb_a._new()
                sig_on()
                arb_fmpz_poly_evaluate_arb(arb_z.value, self._poly, arb_a.value, arb_a._parent._prec)
                sig_off()
                return arb_z
            if isinstance(x0, ComplexBall):
                acb_a = <ComplexBall> x0
                acb_z = acb_a._new()
                sig_on()
                arb_fmpz_poly_evaluate_acb(acb_z.value, self._poly, acb_a.value, acb_a._parent._prec)
                sig_off()
                return acb_z

        return Polynomial.__call__(self, *x, **kwds)

    cpdef Integer content(self):
        r"""
        Return the greatest common divisor of the coefficients of this
        polynomial. The sign is the sign of the leading coefficient.  The
        content of the zero polynomial is zero.

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

        Verify that :issue:`13053` has been resolved::

            sage: R(-1).content()
            -1
        """
        if self.is_zero():
            return ZZ.zero()

        cdef fmpz_t c
        fmpz_init(c)
        fmpz_poly_get_coeff_fmpz(c, self._poly, fmpz_poly_degree(self._poly))
        cdef int sign = fmpz_sgn(c)

        fmpz_poly_content(c, self._poly)

        cdef Integer z = Integer.__new__(Integer)
        fmpz_get_mpz(z.value, c)
        fmpz_clear(c)
        return z if sign == 1 else -z

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

    cdef get_unsafe(self, Py_ssize_t n):
        """
        Return the `n`-th coefficient of ``self``.

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
            sage: f[:4]
            3*x^3 + 2*x^2 + x + 1
            sage: f[:100]
            5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x + 1
        """
        cdef Integer z = Integer.__new__(Integer)
        fmpz_poly_get_coeff_mpz(z.value, self._poly, n)
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
        cdef Integer coef = Integer.__new__(Integer)
        cdef list all = []
        for i from fmpz_poly_degree(self._poly) >= i >= 0:
            fmpz_poly_get_coeff_mpz(coef.value, self._poly, i)
            if coef:
                if coef > 0:
                    sign_str = '+'
                    coeff_str = str(coef)
                else:
                    sign_str = '-'
                    coeff_str = str(coef)[1:]
                if i > 0:
                    if coeff_str == '1':
                        coeff_str = ''
                    elif not latex:
                        coeff_str = coeff_str + '*'
                if i > 1:
                    if latex:
                        all.append(" %s %s%s^{%s}" % (sign_str,
                            coeff_str, name, i))
                    else:
                        all.append(" %s %s%s^%s" % (sign_str,
                            coeff_str, name, i))
                elif i == 1:
                    all.append(" %s %s%s" % (sign_str, coeff_str, name))
                else:
                    all.append(" %s %s" % (sign_str, coeff_str))
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

    cpdef _add_(self, right):
        r"""
        Return ``self`` plus ``right``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 2*x + 1
            sage: g = -3*x^2 + 6
            sage: f + g
            -3*x^2 + 2*x + 7
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_add(x._poly, self._poly,
                (<Polynomial_integer_dense_flint>right)._poly)
        sig_off()
        return x

    cpdef _sub_(self, right):
        r"""
        Return ``self`` minus ``right``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 2*x + 1
            sage: g = -3*x^2 + 6
            sage: f - g
            3*x^2 + 2*x - 5
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_sub(x._poly, self._poly,
                (<Polynomial_integer_dense_flint>right)._poly)
        sig_off()
        return x

    cpdef _neg_(self):
        r"""
        Return negative of ``self``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 2*x - 1
            sage: -f
            -2*x + 1
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_neg(x._poly, self._poly)
        sig_off()
        return x

    @coerce_binop
    def quo_rem(self, Polynomial_integer_dense_flint right):
        r"""
        Attempts to divide ``self`` by ``right``, and return a quotient and remainder.

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

        :issue:`383`, make sure things get coerced correctly::

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
            raise ZeroDivisionError("division by zero polynomial")

        if self.is_zero():
            return self, self

        cdef Polynomial_integer_dense_flint qq = self._new()
        cdef Polynomial_integer_dense_flint rr = self._new()

        sig_on()
        fmpz_poly_divrem(qq._poly, rr._poly, self._poly, right._poly)
        sig_off()
        return qq, rr

    cpdef bint is_zero(self) except -1:
        """
        Return ``True`` if ``self`` is equal to zero.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: R(0).is_zero()
            True
            sage: R(1).is_zero()
            False
            sage: x.is_zero()
            False
        """
        return (fmpz_poly_degree(self._poly) == -1)

    cpdef bint is_one(self) except -1:
        """
        Return ``True`` if ``self`` is equal to one.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: R(0).is_one()
            False
            sage: R(1).is_one()
            True
            sage: x.is_one()
            False
        """
        return fmpz_poly_is_one(self._poly)

    def __bool__(self):
        """
        Check if ``self`` is not zero.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: bool(R(0))
            False
            sage: bool(R(1))
            True
            sage: bool(x)
            True
        """
        return not (fmpz_poly_degree(self._poly) == -1)

    @coerce_binop
    def gcd(self, right):
        r"""
        Return the GCD of ``self`` and ``right``.  The leading
        coefficient need not be 1.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = (6*x + 47) * (7*x^2 - 2*x + 38)
            sage: g = (6*x + 47) * (3*x^3 + 2*x + 1)
            sage: f.gcd(g)
            6*x + 47
        """
        cdef Polynomial_integer_dense_flint _right = <Polynomial_integer_dense_flint> right
        if self.is_zero() or _right.is_one():
            return right
        elif self.is_one() or _right.is_zero():
            return self
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_gcd(x._poly, self._poly,
                (<Polynomial_integer_dense_flint>right)._poly)
        sig_off()
        return x

    @coerce_binop
    def lcm(self, right):
        """
        Return the LCM of ``self`` and ``right``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = (6*x + 47) * (7*x^2 - 2*x + 38)
            sage: g = (6*x + 47) * (3*x^3 + 2*x + 1)
            sage: h = f.lcm(g); h
            126*x^6 + 951*x^5 + 486*x^4 + 6034*x^3 + 585*x^2 + 3706*x + 1786
            sage: h == (6*x + 47) * (7*x^2 - 2*x + 38) * (3*x^3 + 2*x + 1)
            True

        TESTS:

        Check that :issue:`32033` has been fixed::

            sage: R.<t> = ZZ[]
            sage: lcm(R(0), R(0))
            0
        """
        if self.is_zero() or right.is_zero():
            P = self.parent()
            return P.zero()
        g = self.gcd(right)
        return (self//g)*right

    @coerce_binop
    def xgcd(self, right):
        r"""
        Return a triple `(g,s,t)` such that `g = s\cdot{}` ``self`` + `t\cdot{}` ``right`` and such
        that `g` is the gcd of ``self`` and ``right`` up to a divisor of the
        resultant of ``self`` and ``other``.

        As integer polynomials do not form a principal ideal domain, it is not
        always possible given `a` and `b` to find a pair `s,t` such that
        `gcd(a,b) = sa + tb`. Take `a=x+2` and `b=x+4` as an example for which the
        gcd is `1` but the best you can achieve in the Bezout identity is `2`.

        If ``self`` and ``right`` are coprime as polynomials over the
        rationals, then ``g`` is guaranteed to be the resultant of
        ``self`` and ``right``, as a constant polynomial.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(ZZ)

            sage: (x + 2).xgcd(x + 4)
            (2, -1, 1)
            sage: (x + 2).resultant(x + 4)
            2
            sage: (x + 2).gcd(x + 4)
            1

            sage: F = (x^2 + 2)*x^3; G = (x^2 + 2) * (x - 3)
            sage: g, u, v = F.xgcd(G)
            sage: g, u, v
            (27*x^2 + 54, 1, -x^2 - 3*x - 9)
            sage: u*F + v*G
            27*x^2 + 54

            sage: zero = P(0)
            sage: x.xgcd(zero)
            (x, 1, 0)
            sage: zero.xgcd(x)
            (x, 0, 1)

            sage: F = (x - 3)^3; G = (x - 15)^2
            sage: g, u, v = F.xgcd(G)
            sage: g, u, v
            (2985984, -432*x + 8208, 432*x^2 + 864*x + 14256)
            sage: u*F + v*G
            2985984

        TESTS:

        Check that :issue:`17675` is fixed::

            sage: R.<x> = ZZ['x']
            sage: R(2).xgcd(R(2))
            (2, 0, 1)
            sage: R.zero().xgcd(R(2))
            (2, 0, 1)
            sage: R(2).xgcd(R.zero())
            (2, 1, 0)
        """
        # trivial cases first
        if self.is_zero():
            return (right, Integer(0), Integer(1))
        elif right.is_zero():
            return (self, Integer(1), Integer(0))
        elif self.is_constant() and right.is_constant():
            # this case is needed as the behavior of the function fmpz_poly_xgcd
            # in FLINT is not defined with constant input
            return self.constant_coefficient().xgcd(right.constant_coefficient())

        cdef Polynomial_integer_dense_flint ss = self._new()
        cdef Polynomial_integer_dense_flint tt = self._new()
        cdef fmpz_t r
        fmpz_init(r)

        sig_on()
        fmpz_poly_xgcd(r, ss._poly, tt._poly, self._poly,
                (<Polynomial_integer_dense_flint>right)._poly)
        sig_off()
        cdef Integer rr = Integer.__new__(Integer)
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

    cpdef _mul_(self, right):
        r"""
        Return ``self`` multiplied by ``right``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: (x - 2)*(x^2 - 8*x + 16)
            x^3 - 10*x^2 + 32*x - 32
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_mul(x._poly, self._poly,
                (<Polynomial_integer_dense_flint>right)._poly)
        sig_off()
        return x

    cpdef Polynomial _mul_trunc_(self, Polynomial right, long n):
        r"""
        Truncated multiplication.

        .. SEEALSO::

            :meth:`_mul_` for standard multiplication

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: p1 = 1 + x + x^2 + x^4
            sage: p2 = -2 + 3*x^2 + 5*x^4
            sage: p1._mul_trunc_(p2, 4)
            3*x^3 + x^2 - 2*x - 2
            sage: (p1*p2).truncate(4)
            3*x^3 + x^2 - 2*x - 2
            sage: p1._mul_trunc_(p2, 6)
            5*x^5 + 6*x^4 + 3*x^3 + x^2 - 2*x - 2
        """
        if n <= 0:
            raise ValueError("length must be > 0")

        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_mullow(x._poly, self._poly,
                    (<Polynomial_integer_dense_flint>right)._poly,
                    n)
        sig_off()
        return x

    cpdef _lmul_(self, Element right):
        r"""
        Return ``self`` multiplied by ``right``, where ``right`` is a scalar (integer).

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: x*3
            3*x
            sage: (2*x^2 + 4)*3
            6*x^2 + 12
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_scalar_mul_mpz(x._poly, self._poly, (<Integer>right).value)
        sig_off()
        return x

    cpdef _rmul_(self, Element right):
        r"""
        Return ``self`` multiplied by ``right``, where right is a scalar (integer).

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: 3*x
            3*x
            sage: 3*(2*x^2 + 4)
            6*x^2 + 12
        """
        cdef Polynomial_integer_dense_flint x = self._new()
        sig_on()
        fmpz_poly_scalar_mul_mpz(x._poly, self._poly, (<Integer>right).value)
        sig_off()
        return x

    def __pow__(Polynomial_integer_dense_flint self, exp, ignored):
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

        Check that :issue:`18278` is fixed::

            sage: R.<x> = ZZ[]
            sage: x^(1/2)
            Traceback (most recent call last):
            ...
            ValueError: not a 2nd power

            sage: x^(2^100)
            Traceback (most recent call last):
            ...
            OverflowError: Sage Integer too large to convert to C long

        Test fractional powers (:issue:`20086`)::

            sage: P.<R> = ZZ[]
            sage: (R^3 + 6*R^2 + 12*R + 8)^(1/3)
            R + 2
            sage: _.parent()
            Univariate Polynomial Ring in R over Integer Ring
            sage: P(4)^(1/2)
            2
            sage: _.parent()
            Univariate Polynomial Ring in R over Integer Ring

            sage: (R^2 + 3)^(1/2)
            Traceback (most recent call last):
            ...
            ValueError: 3 is not a 2nd power

            Ring in R over Integer Ring
            sage: P(2)^P(2)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Univariate Polynomial
            Ring in R over Integer Ring to Rational Field
            sage: (R + 1)^P(2)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Univariate Polynomial
            Ring in R over Integer Ring to Rational Field
            sage: (R + 1)^R
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Univariate Polynomial
            Ring in R over Integer Ring to Rational Field
            sage: 2^R
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Univariate Polynomial
            Ring in R over Integer Ring to Rational Field
        """
        cdef long nn
        cdef Polynomial_integer_dense_flint res

        try:
            nn = pyobject_to_long(exp)
        except TypeError:
            n = QQ.coerce(exp)
            num = n.numerator()
            den = n.denominator()

            if fmpz_poly_degree(self._poly) == 0:
                return self.parent()(self[0].nth_root(den) ** num)

            return self.nth_root(den) ** num

        else:
            res = self._new()

            if self.is_zero():
                if nn == 0:
                    fmpz_poly_set_coeff_si(res._poly, 0, 1)
                    return res
                elif nn < 0:
                    raise ZeroDivisionError("negative exponent in power of zero")
                else:
                    return res
            if nn < 0:
                sig_on()
                fmpz_poly_pow(res._poly, self._poly, -nn)
                sig_off()
                return ~res
            else:
                if self is self._parent.gen():
                    sig_on()
                    fmpz_poly_set_coeff_ui(res._poly, nn, 1)
                    sig_off()
                else:
                    sig_on()
                    fmpz_poly_pow(res._poly, self._poly, nn)
                    sig_off()
                return res

    cpdef Polynomial _power_trunc(self, unsigned long n, long prec):
        r"""
        Truncated power.

        TESTS::

            sage: R.<x> = ZZ[]
            sage: (x**2 - x + 1).power_trunc(100, 5)
            4411275*x^4 - 171600*x^3 + 5050*x^2 - 100*x + 1
            sage: R.zero().power_trunc(0, 1)
            1

            sage: x._power_trunc(2, -1)
            0
            sage: parent(_) is R
            True
        """
        if prec <= 0:
            # NOTE: flint crashes for prec < 0
            return self._parent.zero()

        cdef Polynomial_integer_dense_flint res
        res = self._new()
        fmpz_poly_pow_trunc(res._poly, self._poly, n, prec)
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
        cdef Polynomial_integer_dense_flint _right
        cdef Polynomial_integer_dense_flint res
        if not isinstance(right, Polynomial_integer_dense_flint):
            if isinstance(right, Integer):
                if (<Integer> right).is_zero():
                    raise ZeroDivisionError("division by zero")
                else:
                    res = self._new()
                    sig_on()
                    fmpz_poly_scalar_fdiv_mpz(res._poly, self._poly,
                            (<Integer> right).value)
                    sig_off()
                return res
            elif right in ZZ:
                res = self._new()
                sig_on()
                fmpz_poly_scalar_fdiv_mpz(res._poly, self._poly,
                        (<Integer>ZZ(right)).value)
                sig_off()
                return res
            else:
                right = self._parent(right)
        _right = <Polynomial_integer_dense_flint> right
        if _right.is_zero():
            raise ZeroDivisionError("division by zero")
        elif _right.is_one():
            # quite typical when removing gcds...
            return self
        else:
            res = self._new()
            sig_on()
            fmpz_poly_div(res._poly, self._poly, _right._poly)
            sig_off()
            return res

    cpdef Polynomial inverse_series_trunc(self, long prec):
        r"""
        Return a polynomial approximation of precision ``prec`` of the inverse
        series of this polynomial.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: p = 1 + x + 2*x^2
            sage: q5 = p.inverse_series_trunc(5)
            sage: q5
            -x^4 + 3*x^3 - x^2 - x + 1
            sage: p*q5
            -2*x^6 + 5*x^5 + 1

            sage: (x-1).inverse_series_trunc(5)
            -x^4 - x^3 - x^2 - x - 1

            sage: q100 = p.inverse_series_trunc(100)
            sage: (q100 * p).truncate(100)
            1

        TESTS::

            sage: ZZ['x'].zero().inverse_series_trunc(4)
            Traceback (most recent call last):
            ...
            ValueError: constant term is zero
            sage: ZZ['x'](2).inverse_series_trunc(4)
            Traceback (most recent call last):
            ...
            ValueError: constant term 2 is not a unit
            sage: x = polygen(ZZ)
            sage: (x+1).inverse_series_trunc(0)
            Traceback (most recent call last):
            ...
            ValueError: the precision must be positive, got 0
        """
        if prec <= 0:
            raise ValueError("the precision must be positive, got {}".format(prec))

        if fmpz_poly_degree(self._poly) == -1:
            raise ValueError("constant term is zero")
        cdef fmpz_t c = fmpz_poly_get_coeff_ptr(self._poly, 0)
        if fmpz_cmp_si(c, 1) and fmpz_cmp_si(c, -1):
            raise ValueError("constant term {} is not a unit".format(self[0]))

        cdef Polynomial_integer_dense_flint res = self._new()
        if prec <= 0:
            return res
        sig_on()
        fmpz_poly_inv_series(res._poly, self._poly, prec)
        sig_off()
        return res

    cpdef _unsafe_mutate(self, long n, value):
        r"""
        Set coefficient of `x^n` to value.

        This is very unsafe, because Sage polynomials are supposed to be
        immutable.

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
            raise IndexError("n must be >= 0")
        if isinstance(value, int):
            sig_on()
            fmpz_poly_set_coeff_si(self._poly, n, value)
            sig_off()
        elif isinstance(value, Integer):
            sig_on()
            fmpz_poly_set_coeff_mpz(self._poly, n, (<Integer>value).value)
            sig_off()
        else:
            value = Integer(value)
            sig_on()
            fmpz_poly_set_coeff_mpz(self._poly, n, (<Integer>value).value)
            sig_off()

    def real_root_intervals(self):
        """
        Return isolating intervals for the real roots of this
        polynomial.

        EXAMPLES:
        We compute the roots of the characteristic polynomial of some
        Salem numbers::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = 1 - x^2 - x^3 - x^4 + x^6
            sage: f.real_root_intervals()                                               # needs sage.libs.linbox
            [((1/2, 3/4), 1), ((1, 3/2), 1)]
        """

        from sage.rings.polynomial.real_roots import real_roots

        return real_roots(self)

    #     def __copy__(self):
    #         f = Polynomial_integer_dense(self.parent())
    #         f._poly = self._poly.copy()
    #         return f

    def degree(self, gen=None):
        """
        Return the degree of this polynomial.

        The zero polynomial has degree `-1`.

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

        TESTS::

            sage: type(x.degree())
            <class 'sage.rings.integer.Integer'>
        """
        return smallInteger(fmpz_poly_degree(self._poly))

    def pseudo_divrem(self, B):
        r"""
        Write ``A = self``.  This function computes polynomials `Q` and `R`
        and an integer `d` such that

        .. MATH::

             \mathop{\mathrm{lead}}(B)^d A = B Q + R

        where R has degree less than that of B.

        INPUT:

        - ``B`` -- a polynomial over `\ZZ`

        OUTPUT:

        - ``Q``, ``R`` -- polynomials
        - ``d`` -- nonnegative integer

        EXAMPLES::

            sage: R.<x> = ZZ['x']
            sage: A = R(range(10))
            sage: B = 3*R([-1, 0, 1])
            sage: Q, R, d = A.pseudo_divrem(B)
            sage: Q, R, d
            (9*x^7 + 8*x^6 + 16*x^5 + 14*x^4 + 21*x^3 + 18*x^2 + 24*x + 20, 75*x + 60, 1)
            sage: B.leading_coefficient()^d * A == B*Q + R
            True
        """
        cdef Polynomial_integer_dense_flint Q = self._new(), R = self._new(), _B = B
        cdef ulong d
        fmpz_poly_pseudo_divrem(Q._poly, R._poly, &d, self._poly, _B._poly)
        return Q, R, Integer(d)

    def discriminant(self, proof=True):
        r"""
        Return the discriminant of ``self``, which is by definition

        .. MATH::

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

        TESTS:

        Confirm that :issue:`17603` has been applied::

            sage: f.disc()
            -339
        """
        cdef ZZX_c ntl_poly
        cdef ZZ_c* temp
        fmpz_poly_get_ZZX(ntl_poly, self._poly)

        temp = ZZX_discriminant(&ntl_poly, proof)
        cdef Integer x = Integer.__new__(Integer)
        ZZ_to_mpz(x.value, temp)
        del temp

        return x

    # Alias for discriminant
    disc = discriminant

    def __pari__(self, variable=None):
        """
        EXAMPLES::

            sage: t = PolynomialRing(ZZ,"t").gen()
            sage: f = t^3 + 3*t - 17
            sage: pari(f)
            t^3 + 3*t - 17
            sage: f.__pari__(variable='y')
            y^3 + 3*y - 17
        """
        if variable is None:
            variable = self.parent().variable_name()
        return pari(self.list()).Polrev(variable)

    def squarefree_decomposition(Polynomial_integer_dense_flint self):
        """
        Return the square-free decomposition of ``self``.  This is
        a partial factorization of ``self`` into square-free, relatively
        prime polynomials.

        This is a wrapper for the NTL function ``SquareFreeDecomp``.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: p = (x-1)^2 * (x-2)^2 * (x-3)^3 * (x-4)
            sage: p.squarefree_decomposition()
            (x - 4) * (x^2 - 3*x + 2)^2 * (x - 3)^3
            sage: p = 37 * (x-1)^2 * (x-2)^2 * (x-3)^3 * (x-4)
            sage: p.squarefree_decomposition()
            (37) * (x - 4) * (x^2 - 3*x + 2)^2 * (x - 3)^3

        TESTS:

        Verify that :issue:`13053` has been resolved::

            sage: R.<x> = PolynomialRing(ZZ, implementation='FLINT')
            sage: f=-x^2
            sage: f.squarefree_decomposition()
            (-1) * x^2
        """
        cdef ZZX_c** v
        cdef long* e
        cdef long i, n
        cdef fmpz_poly_t ppart
        cdef ZZX_c ntl_poly
        cdef Integer z
        cdef Polynomial_integer_dense_flint fac

        # the sign of the content is the sign of the leading coefficient
        z = self.content()
        if not z.is_one():
            fmpz_poly_init(ppart)

            # the primitive part returned by FLINT has positive leading
            # coefficient
            fmpz_poly_primitive_part(ppart, self._poly)

            fmpz_poly_get_ZZX(ntl_poly, ppart)
            fmpz_poly_clear(ppart)
        else:
            fmpz_poly_get_ZZX(ntl_poly, self._poly)

        # input is primitive, with positive leading coefficient
        ZZX_squarefree_decomposition(&v, &e, &n, &ntl_poly)

        F = []
        for i from 0 <= i < n:
            fac = self._new()
            fmpz_poly_set_ZZX(fac._poly, v[i][0])
            F.append( (fac,e[i]) )
            del v[i]
        sig_free(v)
        sig_free(e)

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
        cdef int sig_me = fmpz_poly_degree(self._poly)

        fmpz_poly_get_ZZX(ntl_poly, self._poly)

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
            fmpz_poly_set_coeff_fmpz(fac_py._poly, 0, tcontent)
            results.append( (fac_py,1) )
            fmpz_clear(tcontent)

        for i from 0 <= i < factors.length():
            fac_py = self._new()
            fmpz_poly_set_ZZX(fac_py._poly, factors.RawGet(i).a)
            results.append( (fac_py,factors.RawGet(i).b) )
        return Factorization(results, unit = unit)

    def factor(self):
        """
        This function overrides the generic polynomial factorization to
        make a somewhat intelligent decision to use PARI or NTL based on
        some benchmarking.

        Note: This function factors the content of the polynomial,
        which can take very long if it's a really big integer.  If you
        do not need the content factored, divide it out of your
        polynomial before calling this function.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: f = x^4 - 1
            sage: f.factor()
            (x - 1) * (x + 1) * (x^2 + 1)
            sage: f = 1 - x
            sage: f.factor()
            (-1) * (x - 1)
            sage: f.factor().unit()
            -1
            sage: f = -30*x; f.factor()
            (-1) * 2 * 3 * 5 * x
        """
        cdef int i
        cdef long deg = fmpz_poly_degree(self._poly)
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
        Return the factorization of ``self`` modulo the prime `p`.

        INPUT:

        - ``p`` -- prime

        OUTPUT: factorization of ``self`` reduced modulo `p`

        EXAMPLES::

            sage: R.<x> = ZZ['x']
            sage: f = -3*x*(x-2)*(x-9) + x
            sage: f.factor_mod(3)
            x
            sage: f = -3 * x * (x - 2) * (x - 9)
            sage: f.factor_mod(3)
            Traceback (most recent call last):
            ...
            ArithmeticError: factorization of 0 is not defined

            sage: f = 2 * x * (x - 2) * (x - 9)
            sage: f.factor_mod(7)
            (2) * x * (x + 5)^2
        """
        from sage.rings.finite_rings.finite_field_constructor import FiniteField
        p = Integer(p)
        if not p.is_prime():
            raise ValueError("p must be prime")
        if not any(c % p for c in self.coefficients()):
            raise ArithmeticError("factorization of 0 is not defined")
        f = self.__pari__()
        G = f.factormod(p)
        k = FiniteField(p)
        R = k[self.parent().variable_name()]
        return R(1)._factor_pari_helper(G, unit=R(self).leading_coefficient())

    def factor_padic(self, p, prec=10):
        """
        Return `p`-adic factorization of ``self`` to given precision.

        INPUT:

        - ``p`` -- prime

        - ``prec`` -- integer; the precision

        OUTPUT: factorization of ``self`` over the completion at `p`

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: f = x^2 + 1
            sage: f.factor_padic(5, 4)
            ((1 + O(5^4))*x + 2 + 5 + 2*5^2 + 5^3 + O(5^4))
            * ((1 + O(5^4))*x + 3 + 3*5 + 2*5^2 + 3*5^3 + O(5^4))

        A more difficult example::

            sage: f = 100 * (5*x + 1)^2 * (x + 5)^2
            sage: f.factor_padic(5, 10)
            (4 + O(5^10)) * (5 + O(5^11))^2 * ((1 + O(5^10))*x + 5 + O(5^10))^2
            * ((5 + O(5^10))*x + 1 + O(5^10))^2
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

    cpdef list list(self, bint copy=True):
        """
        Return a new copy of the list of the underlying
        elements of ``self``.

        EXAMPLES::

            sage: x = PolynomialRing(ZZ,'x').0
            sage: f = x^3 + 3*x - 17
            sage: f.list()
            [-17, 3, 0, 1]
            sage: f = PolynomialRing(ZZ,'x')(0)
            sage: f.list()
            []
        """
        return [self.get_unsafe(i) for i in range(self.degree()+1)]

    @coerce_binop
    def resultant(self, other, proof=True):
        """
        Return the resultant of ``self`` and ``other``, which must lie in the same
        polynomial ring.

        If ``proof=False`` (the default is ``proof=True``), then this function may
        use a randomized strategy that errors with probability no more than
        `2^{-80}`.

        INPUT:

        - ``other`` -- a polynomial

        OUTPUT: an element of the base ring of the polynomial ring

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
        cdef Integer x = Integer.__new__(Integer)

        sig_on()
        fmpz_poly_resultant(res, self._poly,
                (<Polynomial_integer_dense_flint>other)._poly)
        sig_off()
        fmpz_get_mpz(x.value, res)
        fmpz_clear(res)
        return x

    def reverse(self, degree=None):
        """
        Return a polynomial with the coefficients of this polynomial reversed.

        If an optional degree argument is given the coefficient list will be
        truncated or zero padded as necessary before computing the reverse.

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
            ValueError: degree argument must be a nonnegative integer, got 1.5

        Check that this implementation is compatible with the generic one::

            sage: p = R([0,1,0,2])
            sage: all(p.reverse(d) == Polynomial.reverse(p, d)
            ....:     for d in [None, 0, 1, 2, 3, 4])
            True
        """
        cdef Polynomial_integer_dense_flint res = self._new()
        cdef unsigned long d
        if degree is not None:
            if degree < 0:
                raise ValueError("degree argument must be a nonnegative integer, got %s" % (degree))
            d = degree
            if d != degree:
                raise ValueError("degree argument must be a nonnegative integer, got %s" % (degree))
            # FLINT expects length
            fmpz_poly_reverse(res._poly, self._poly, d+1)
        else:
            fmpz_poly_reverse(res._poly, self._poly,
                    fmpz_poly_length(self._poly))
        return res

    def revert_series(self, n):
        r"""
        Return a polynomial `f` such that `f(` ``self`` `(x)) =` ``self`` `(f(x)) = x` (mod `x^n`).

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: f = t - t^3 + t^5
            sage: f.revert_series(6)
            2*t^5 + t^3 + t

            sage: f.revert_series(-1)
            Traceback (most recent call last):
            ...
            ValueError: argument n must be a nonnegative integer, got -1

            sage: g = - t^3 + t^5
            sage: g.revert_series(6)
            Traceback (most recent call last):
            ...
            ValueError: self must have constant coefficient 0 and a unit for coefficient t^1
        """
        cdef Polynomial_integer_dense_flint res = self._new()
        cdef unsigned long m
        if n < 0:
            raise ValueError("argument n must be a nonnegative integer, got {}".format(n))
        m = n
        if not self[0].is_zero() or not self[1].is_unit():
            raise ValueError("self must have constant coefficient 0 and a unit for coefficient {}^1".format(self.parent().gen()))

        sig_on()
        fmpz_poly_revert_series(res._poly, self._poly, m)
        sig_off()

        return res
