# -*- coding: utf-8
r"""
Univariate polynomials over `\CC` with interval coefficients using Arb.

This is a binding to the `Arb library <http://fredrikj.net/arb/>`_; it
may be useful to refer to its documentation for more details.

Parts of the documentation for this module are copied or adapted from
Arb's own documentation, licenced under the GNU General Public License
version 2, or later.

.. SEEALSO::

    - :mod:`Complex balls using Arb <sage.rings.complex_arb>`

TESTS:

    sage: type(polygen(ComplexBallField(140)))
    <type 'sage.rings.polynomial.polynomial_complex_arb.Polynomial_complex_arb'>
    sage: Pol.<x> = CBF[]
    sage: (x+1/2)^3
    x^3 + 1.500000000000000*x^2 + 0.7500000000000000*x + 0.1250000000000000

"""

include "cysignals/signals.pxi"

from sage.libs.arb.acb cimport *
from sage.rings.integer cimport Integer, smallInteger
from sage.rings.complex_arb cimport ComplexBall

from sage.rings.complex_arb import ComplexBallField
from sage.structure.element import coerce_binop, have_same_parent

cdef inline long prec(Polynomial_complex_arb pol):
    return pol._parent._base._prec

cdef class Polynomial_complex_arb(Polynomial):
    r"""
    Wrapper for `Arb <http://fredrikj.net/arb/>`_ polynomials of type
    ``acb_poly_t``

    EXAMPLES::

        sage: Pol.<x> = CBF[]
        sage: type(x)
        <type 'sage.rings.polynomial.polynomial_complex_arb.Polynomial_complex_arb'>

        sage: Pol(), Pol(1), Pol([0,1,2]), Pol({1: pi, 3: i})
        (0,
         1.000000000000000,
         2.000000000000000*x^2 + x,
         I*x^3 + ([3.141592653589793 +/- 5.61e-16])*x)

        sage: Pol("x - 2/3")
        x + [-0.666666666666667 +/- 4.82e-16]
        sage: Pol(polygen(QQ))
        x

        sage: [Pol.has_coerce_map_from(P) for P in
        ....: QQ['x'], QuadraticField(-1), RealBallField(100)]
        [True, True, True]
        sage: [Pol.has_coerce_map_from(P) for P in
        ....: QQ['y'], RR, CC, RDF, CDF, RIF, CIF, RealBallField(20)]
        [False, False, False, False, False, False, False, False]
    """

    # Memory management and initialization

    def __cinit__(self):
        r"""
        TESTS::

            sage: ComplexBallField(2)['y']()
            0
        """
        acb_poly_init(self.__poly)

    def __dealloc__(self):
        r"""
        TESTS::

            sage: pol = CBF['x']()
            sage: del pol
        """
        acb_poly_clear(self.__poly)

    cdef Polynomial_complex_arb _new(self):
        r"""
        Return a new polynomial with the same parent as this one.
        """
        cdef Polynomial_complex_arb res = Polynomial_complex_arb.__new__(Polynomial_complex_arb)
        res._parent = self._parent
        res._is_gen = 0
        return res

    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        r"""
        Initialize this polynomial to the specified value.

        TESTS::

            sage: from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
            sage: Pol = CBF['x']
            sage: Polynomial_complex_arb(Pol)
            0
            sage: Polynomial_complex_arb(Pol, is_gen=True)
            x
            sage: Polynomial_complex_arb(Pol, 42, is_gen=True)
            x
            sage: Polynomial_complex_arb(Pol, CBF(1))
            1.000000000000000
            sage: Polynomial_complex_arb(Pol, [])
            0
            sage: Polynomial_complex_arb(Pol, [0])
            0
            sage: Polynomial_complex_arb(Pol, [0, 2, 0])
            2.000000000000000*x
            sage: Polynomial_complex_arb(Pol, (1,))
            1.000000000000000
            sage: Polynomial_complex_arb(Pol, (CBF(i), 1))
            x + I
            sage: Polynomial_complex_arb(Pol, polygen(QQ,'y')+2)
            x + 2.000000000000000
            sage: Polynomial_complex_arb(Pol, QQ['x'](0))
            0
            sage: Polynomial_complex_arb(Pol, {10: pi})
            ([3.141592653589793 +/- 5.61e-16])*x^10
            sage: Polynomial_complex_arb(Pol, pi)
            [3.141592653589793 +/- 5.61e-16]
        """
        cdef ComplexBall ball
        cdef Polynomial pol
        cdef list lst
        cdef tuple tpl
        cdef dict dct
        cdef long length, i

        Polynomial.__init__(self, parent, is_gen=is_gen)

        if is_gen:
            acb_poly_set_coeff_si(self.__poly, 1, 1)
        elif x is None:
            acb_poly_zero(self.__poly)
        elif isinstance(x, Polynomial_complex_arb):
            acb_poly_set(self.__poly, (<Polynomial_complex_arb> x).__poly)
        elif isinstance(x, ComplexBall):
            acb_poly_set_coeff_acb(self.__poly, 0, (<ComplexBall> x).value)
        else:
            Coeff = parent.base_ring()
            if isinstance(x, list):
                lst = <list> x
                length = len(lst)
                sig_on(); acb_poly_fit_length(self.__poly, length); sig_off()
                for i in range(length):
                    ball = Coeff(lst[i])
                    acb_poly_set_coeff_acb(self.__poly, i, ball.value)
            elif isinstance(x, tuple):
                tpl = <tuple> x
                length = len(tpl)
                sig_on(); acb_poly_fit_length(self.__poly, length); sig_off()
                for i in range(length):
                    ball = Coeff(tpl[i])
                    acb_poly_set_coeff_acb(self.__poly, i, ball.value)
            elif isinstance(x, Polynomial):
                pol = <Polynomial> x
                length = pol.degree() + 1
                sig_on(); acb_poly_fit_length(self.__poly, length); sig_off()
                for i in range(length):
                    ball = Coeff(pol.get_unsafe(i))
                    acb_poly_set_coeff_acb(self.__poly, i, ball.value)
            elif isinstance(x, dict):
                dct = <dict> x
                if len(dct) == 0:
                    acb_poly_zero(self.__poly)
                else:
                    length = max(int(i) for i in dct) + 1
                    sig_on(); acb_poly_fit_length(self.__poly, length); sig_off()
                    for i, c in dct.iteritems():
                        ball = Coeff(c)
                        acb_poly_set_coeff_acb(self.__poly, i, ball.value)
            else:
                ball = Coeff(x)
                acb_poly_set_coeff_acb(self.__poly, 0, ball.value)

    # Access

    def degree(self):
        r"""
        Return the (apparent) degree of this polynomial.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x^2 + 1).degree()
            2
            sage: pol = (x/3 + 1) - x/3; pol
            ([+/- 1.12e-16])*x + 1.000000000000000
            sage: pol.degree()
            1
            sage: Pol([1, 0, 0, 0]).degree()
            0
        """
        return smallInteger(acb_poly_degree(self.__poly))

    cdef get_unsafe(self, Py_ssize_t n):
        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        res._parent = self._parent._base
        acb_poly_get_coeff_acb(res.value, self.__poly, n)
        return res

    def list(self):
        r"""
        Return the coefficient list of this polynomial.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x^2/3).list()
            [0, 0, [0.3333333333333333 +/- 7.04e-17]]
            sage: Pol(0).list()
            []
            sage: Pol([0, 1, RBF(0, rad=.1), 0]).list()
            [0, 1.000000000000000, [+/- 0.101]]
        """
        cdef unsigned long length = acb_poly_length(self.__poly)
        return [self.get_unsafe(n) for n in range(length)]

    def __nonzero__(self):
        r"""
        Return ``False`` if this polynomial is exactly zero, ``True`` otherwise.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: bool(Pol(0))
            False
            sage: z = Pol(1/3) - 1/3
            sage: bool(z)
            True
        """
        return acb_poly_length(self.__poly)

    # Ring and Euclidean arithmetic

    cpdef _add_(self, other):
        r"""
        Return the sum of two polynomials.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x + 1) + (x/3 - 2)
            ([1.333333333333333 +/- 5.37e-16])*x - 1.000000000000000
        """
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_add(
                res.__poly,
                self.__poly,
                (<Polynomial_complex_arb> other).__poly,
                prec(self))
        sig_off()
        return res

    cpdef _neg_(self):
        r"""
        Return the opposite of this polynomial.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: -(x/3 - 2)
            ([-0.3333333333333333 +/- 7.04e-17])*x + 2.000000000000000
        """
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_neg(res.__poly, self.__poly)
        sig_off()
        return res

    cpdef _sub_(self, other):
        r"""
        Return the difference of two polynomials.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x + 1) - (x/3 - 2)
            ([0.666666666666667 +/- 5.37e-16])*x + 3.000000000000000
        """
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_sub(
                res.__poly,
                self.__poly,
                (<Polynomial_complex_arb> other).__poly,
                prec(self))
        sig_off()
        return res

    cpdef _mul_(self, other):
        r"""
        Return the product of two polynomials.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x + 1)*(x/3 - 2)
            ([0.3333333333333333 +/- 7.04e-17])*x^2
            + ([-1.666666666666667 +/- 7.59e-16])*x - 2.000000000000000
        """
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_mul(
                res.__poly,
                self.__poly,
                (<Polynomial_complex_arb> other).__poly,
                prec(self))
        sig_off()
        return res

    @coerce_binop
    def quo_rem(self, divisor):
        r"""
        Compute the Euclidean division of this ball polynomial by ``divisor``.

        Raises a ``ZeroDivisionError`` when the divisor is zero or its leading
        coefficient contains zero. Returns a pair (quotient, remainder)
        otherwise.

        EXAMPLES::

            sage: Pol.<x> = CBF[]

            sage: (x^3/7 - CBF(i)).quo_rem(x + CBF(pi))
            (([0.1428571428571428 +/- 7.70e-17])*x^2 +
            ([-0.448798950512828 +/- 6.74e-16])*x
            + [1.40994348586991 +/- 3.34e-15],
            [-4.42946809718569 +/- 9.00e-15] - I)

            sage: Pol(0).quo_rem(x + 1)
            (0, 0)

            sage: (x + 1).quo_rem(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ('cannot divide by this polynomial', 0)

            sage: div = (x^2/3 + x + 1) - x^2/3; div
            ([+/- 1.12e-16])*x^2 + x + 1.000000000000000
            sage: (x + 1).quo_rem(div)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ('cannot divide by this polynomial',
            ([+/- 1.12e-16])*x^2 + x + 1.000000000000000)
        """
        cdef Polynomial_complex_arb div = <Polynomial_complex_arb> divisor
        cdef Polynomial_complex_arb quo = self._new()
        cdef Polynomial_complex_arb rem = self._new()
        sig_on()
        cdef bint success = acb_poly_divrem(quo.__poly, rem.__poly, self.__poly,
                div.__poly, prec(self))
        sig_off()
        if success:
            return quo, rem
        else:
            raise ZeroDivisionError("cannot divide by this polynomial", divisor)

    # Syntactic transformations

    cpdef Polynomial truncate(self, long n):
        r"""
        Return the truncation to degree `n - 1` of this polynomial.

        EXAMPLES::

            sage: pol = CBF['x'](range(1,5)); pol
            4.000000000000000*x^3 + 3.000000000000000*x^2 + 2.000000000000000*x + 1.000000000000000
            sage: pol.truncate(2)
            2.000000000000000*x + 1.000000000000000
            sage: pol.truncate(0)
            0
            sage: pol.truncate(-1)
            0

        TESTS::

            sage: pol.truncate(6)
            4.000000000000000*x^3 + 3.000000000000000*x^2 + 2.000000000000000*x + 1.000000000000000
            sage: pol.truncate(4)
            4.000000000000000*x^3 + 3.000000000000000*x^2 + 2.000000000000000*x + 1.000000000000000
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_set(res.__poly, self.__poly)
        acb_poly_truncate(res.__poly, n)
        sig_off()
        return res

    cdef _inplace_truncate(self, long n):
        if n < 0:
            n = 0
        acb_poly_truncate(self.__poly, n)
        return self

    def __lshift__(val, n):
        r"""
        Shift ``val`` to the left, i.e. multiply it by `x^n`, throwing away
        coefficients if `n < 0`.

        EXAMPLES::

            sage: pol = CBF['x'](range(1,5)); pol
            4.000000000000000*x^3 + 3.000000000000000*x^2 + 2.000000000000000*x + 1.000000000000000
            sage: pol << 2
            4.000000000000000*x^5 + 3.000000000000000*x^4 + 2.000000000000000*x^3 + x^2
            sage: pol << (-2)
            4.000000000000000*x + 3.000000000000000

        TESTS::

            sage: 1 << pol
            Traceback (most recent call last):
            ...
            TypeError: unsupported operands for <<: 1, 4.000000000000000*x^3 + 3.000000000000000*x^2 + 2.000000000000000*x + 1.000000000000000
        """
        if not isinstance(val, Polynomial_complex_arb):
            raise TypeError("unsupported operand type(s) for <<: '{}' and '{}'"
                            .format(type(val).__name__, type(n).__name__))
        if n < 0:
            return val.__rshift__(-n)
        cdef Polynomial_complex_arb self = (<Polynomial_complex_arb> val)
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_shift_left(res.__poly, self.__poly, n)
        sig_off()
        return res

    def __rshift__(val, n):
        r"""
        Shift ``val`` to the left, i.e. divide it by `x^n`, throwing away
        coefficients if `n > 0`.

        EXAMPLES::

            sage: pol = CBF['x'](range(1,5)); pol
            4.000000000000000*x^3 + 3.000000000000000*x^2 + 2.000000000000000*x + 1.000000000000000
            sage: pol >> 2
            4.000000000000000*x + 3.000000000000000
            sage: pol >> -2
            4.000000000000000*x^5 + 3.000000000000000*x^4 + 2.000000000000000*x^3 + x^2

        TESTS::

            sage: 1 >> pol
            Traceback (most recent call last):
            ...
            TypeError: unsupported operands for >>: 1, 4.000000000000000*x^3 + 3.000000000000000*x^2 + 2.000000000000000*x + 1.000000000000000
        """
        if not isinstance(val, Polynomial_complex_arb):
            raise TypeError("unsupported operand type(s) for <<: '{}' and '{}'"
                            .format(type(val).__name__, type(n).__name__))
        if n < 0:
            return val.__lshift__(-n)
        cdef Polynomial_complex_arb self = (<Polynomial_complex_arb> val)
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_shift_right(res.__poly, self.__poly, n)
        sig_off()
        return res

    # Truncated and power series arithmetic

    cpdef Polynomial _mul_trunc_(self, Polynomial other, long n):
        r"""
        Return the product of ``self`` and ``other``, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x + 1)._mul_trunc_(x + 2, 2)
            3.000000000000000*x + 2.000000000000000
            sage: (x + 1)._mul_trunc_(x + 2, 0)
            0
            sage: (x + 1)._mul_trunc_(x + 2, -1)
            0

        TESTS::

            sage: (x + 1)._mul_trunc_(x + 2, 4)
            x^2 + 3.000000000000000*x + 2.000000000000000
        """
        cdef Polynomial_complex_arb my_other = <Polynomial_complex_arb> other
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_mullow(res.__poly, self.__poly, my_other.__poly, n, prec(self))
        sig_off()
        return res

    cpdef Polynomial inverse_series_trunc(self, long n):
        r"""
        Return the power series expansion at 0 of the inverse of this
        polynomial, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (1 - x/3).inverse_series_trunc(3)
            ([0.1111111111111111 +/- 5.99e-17])*x^2 + ([0.3333333333333333 +/- 7.04e-17])*x + 1.000000000000000
            sage: x.inverse_series_trunc(1)
            [+/- inf]
            sage: Pol(0).inverse_series_trunc(2)
            (nan + nan*I)*x + nan + nan*I

        TESTS::

            sage: Pol(0).inverse_series_trunc(-1)
            0
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_inv_series(res.__poly, self.__poly, n, prec(self))
        sig_off()
        return res

    cpdef Polynomial _power_trunc(self, unsigned long expo, long n):
        r"""
        Return a power of this polynomial, truncated before degree `n`.

        INPUT:

        - ``expo`` - non-negative integer exponent
        - ``n`` - truncation order

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x^2 + 1)._power_trunc(10^9, 3)
            1000000000.000000*x^2 + 1.000000000000000
            sage: (x^2 + 1)._power_trunc(10^20, 0)
            Traceback (most recent call last):
                ...
            OverflowError: long int too large to convert

        TESTS::

            sage: (x^2 + 1)._power_trunc(10, -3)
            0
            sage: (x^2 + 1)._power_trunc(-1, 0)
            Traceback (most recent call last):
            ...
            OverflowError: can't convert negative value to unsigned long
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_pow_ui_trunc_binexp(res.__poly, self.__poly, expo, n, prec(self))
        sig_off()
        return res

    def _log_series(self, long n):
        r"""
        Return the power series expansion at 0 of the logarithm of this
        polynomial, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (1 + x/3)._log_series(3)
            ([-0.0555555555555555 +/- 7.10e-17])*x^2 + ([0.3333333333333333 +/- 7.04e-17])*x
            sage: (-1 + x)._log_series(3)
            -0.5000000000000000*x^2 - x + [3.141592653589793 +/- 5.61e-16]*I

        An example where the constant term crosses the branch cut of the
        logarithm::

            sage: pol = CBF(-1, RBF(0, rad=.01)) + x; pol
            x - 1.000000000000000 + [+/- 0.0101]*I
            sage: pol._log_series(2)
            ([-1.000 +/- 1.01e-4] + [+/- 0.0101]*I)*x + [+/- 5.01e-5] + [+/- 3.15]*I

        Some cases where the result is not defined::

            sage: x._log_series(1)
            nan + nan*I
            sage: Pol(0)._log_series(1)
            nan + nan*I
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_log_series(res.__poly, self.__poly, n, prec(self))
        sig_off()
        return res

    def _exp_series(self, long n):
        r"""
        Return the power series expansion at 0 of the exponential of this
        polynomial, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: x._exp_series(3)
            0.5000000000000000*x^2 + x + 1.000000000000000
            sage: (1 + x/3)._log_series(3)._exp_series(3)
            ([+/- 5.09e-17])*x^2 + ([0.3333333333333333 +/- 7.04e-17])*x + 1.000000000000000
            sage: (CBF(0, pi) + x)._exp_series(4)
            ([-0.166...] + [+/- ...]*I)*x^3 + ([-0.500...] + [+/- ...]*I)*x^2
            + ([-1.000...] + [+/- ...]*I)*x + [-1.000...] + [+/- ...]*I
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_exp_series(res.__poly, self.__poly, n, prec(self))
        sig_off()
        return res

    def _sqrt_series(self, long n):
        r"""
        Return the power series expansion at 0 of the square root of this
        polynomial, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (1 + x)._sqrt_series(3)
            -0.1250000000000000*x^2 + 0.5000000000000000*x + 1.000000000000000
            sage: pol = CBF(-1, RBF(0, rad=.01)) + x; pol
            x - 1.000000000000000 + [+/- 0.0101]*I
            sage: pol._sqrt_series(2)
            ([+/- 7.51e-3] + [+/- 0.501]*I)*x + [+/- 5.01e-3] + [+/- 1.01]*I
            sage: x._sqrt_series(2)
            ([+/- inf] + [+/- inf]*I)*x
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_sqrt_series(res.__poly, self.__poly, n, prec(self))
        sig_off()
        return res

    def compose_trunc(self, Polynomial other, long n):
        r"""
        Return the composition of ``self`` and ``other``, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: Pol.<x> = CBF[]
            sage: pol = x*(x-1)^2
            sage: pol.compose_trunc(x + x^2, 4)
            -3.000000000000000*x^3 - x^2 + x
            sage: pol.compose_trunc(1 + x, 4)
            x^3 + x^2
            sage: pol.compose_trunc(2 + x/3, 2)
            ([1.666666666666667 +/- 9.81e-16])*x + 2.000000000000000
            sage: pol.compose_trunc(2 + x/3, 0)
            0
            sage: pol.compose_trunc(2 + x/3, -1)
            0
        """
        if n < 0:
            n = 0
        if not isinstance(other, Polynomial_complex_arb):
            return self(other).truncate(n)
        cdef Polynomial_complex_arb other1 = <Polynomial_complex_arb> other
        cdef Polynomial_complex_arb res = self._new()
        cdef acb_poly_t self_ts, other_ts, lin
        cdef acb_ptr cc
        if acb_poly_length(other1.__poly) > 0:
            cc = acb_poly_get_coeff_ptr(other1.__poly, 0)
            if not acb_is_zero(cc):
                sig_on()
                try:
                    acb_poly_init(self_ts)
                    acb_poly_init(other_ts)
                    ### Not yet supported in sage's version of arb
                    #acb_poly_taylor_shift(self_ts, self.__poly, cc, prec(self))
                    acb_poly_init(lin)
                    acb_poly_set_coeff_acb(lin, 0, cc)
                    acb_poly_set_coeff_si(lin, 1, 1)
                    acb_poly_compose(self_ts, self.__poly, lin, prec(self))
                    ###
                    acb_poly_set(other_ts, other1.__poly)
                    acb_zero(acb_poly_get_coeff_ptr(other_ts, 0))
                    acb_poly_compose_series(res.__poly, self_ts, other_ts, n, prec(self))
                finally:
                    ###
                    acb_poly_clear(lin)
                    ###
                    acb_poly_clear(other_ts)
                    acb_poly_clear(self_ts)
                    sig_off()
                return res
        sig_on()
        acb_poly_compose_series(res.__poly, self.__poly, other1.__poly, n, prec(self))
        sig_off()
        return res

    # Evaluation

    def __call__(self, *x, **kwds):
        r"""
        Evaluate this polynomial.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: pol = x^2 - 1
            sage: pol(CBF(pi))
            [8.86960440108936 +/- 8.36e-15]
            sage: pol(x^3 + 1)
            x^6 + 2.000000000000000*x^3
            sage: pol(matrix([[1,2],[3,4]]))
            [6.000000000000000 10.00000000000000]
            [15.00000000000000 21.00000000000000]
        """
        cdef ComplexBall ball
        cdef Polynomial_complex_arb poly
        if len(x) == 1 and not kwds:
            point = x[0]
            if isinstance(point, ComplexBall):
                # parent of result = base ring of self (not parent of point)
                ball = ComplexBall.__new__(ComplexBall)
                ball._parent = self._parent._base
                sig_on()
                acb_poly_evaluate(ball.value, self.__poly,
                        (<ComplexBall> point).value, prec(self))
                sig_off()
                return ball
            elif isinstance(point, Polynomial_complex_arb):
                poly = self._new()
                sig_on()
                acb_poly_compose(poly.__poly, self.__poly,
                        (<Polynomial_complex_arb> point).__poly, prec(self))
                sig_off()
                return poly
            # TODO: perhaps add more special cases, e.g. for real ball,
            # integers and rationals
        return Polynomial.__call__(self, *x, **kwds)
