# sage_setup: distribution = sagemath-flint
r"""
Univariate polynomials over `\CC` with Arb ball coefficients.

This is a binding to the `acb_poly module of FLINT <https://flintlib.org/doc/acb_poly.html>`_;
it may be useful to refer to its documentation for more details.

Parts of the documentation for this module are copied or adapted from Arb's
(now FLINT's) own documentation, licenced (at the time) under the GNU General
Public License version 2, or later.

.. SEEALSO::

    - :mod:`Complex balls using Arb <sage.rings.complex_arb>`

TESTS:

    sage: type(polygen(ComplexBallField(140)))
    <class 'sage.rings.polynomial.polynomial_complex_arb.Polynomial_complex_arb'>
    sage: Pol.<x> = CBF[]
    sage: (x+1/2)^3
    x^3 + 1.500000000000000*x^2 + 0.7500000000000000*x + 0.1250000000000000
"""

from cysignals.signals cimport sig_on, sig_off

from sage.libs.flint.acb cimport *
from sage.libs.flint.fmpz cimport *
from sage.rings.integer cimport Integer, smallInteger
from sage.rings.complex_arb cimport ComplexBall
from sage.structure.element cimport Element

from sage.structure.element import coerce_binop

cdef inline long prec(Polynomial_complex_arb pol) noexcept:
    return pol._parent._base._prec


cdef class Polynomial_complex_arb(Polynomial):
    r"""
    Wrapper for `FLINT <https://flintlib.org>`_ polynomials of type
    ``acb_poly_t``

    EXAMPLES::

        sage: Pol.<x> = CBF[]
        sage: type(x)
        <class 'sage.rings.polynomial.polynomial_complex_arb.Polynomial_complex_arb'>

        sage: Pol(), Pol(1), Pol([0,1,2]), Pol({1: pi, 3: i})                           # needs sage.symbolic
        (0,
         1.000000000000000,
         2.000000000000000*x^2 + x,
         I*x^3 + ([3.141592653589793 +/- ...e-16])*x)

        sage: Pol("x - 2/3")
        x + [-0.666666666666667 +/- ...e-16]
        sage: Pol(polygen(QQ))
        x

        sage: all(Pol.has_coerce_map_from(P) for P in
        ....:     (QQ['x'], QuadraticField(-1), RealBallField(100)))
        True
        sage: any(Pol.has_coerce_map_from(P) for P in
        ....:     (QQ['y'], RR, CC, RDF, CDF, RIF, CIF, RealBallField(20)))
        False
    """

    # Memory management and initialization

    def __cinit__(self):
        r"""
        TESTS::

            sage: ComplexBallField(2)['y']()
            0
        """
        acb_poly_init(self._poly)

    def __dealloc__(self):
        r"""
        TESTS::

            sage: pol = CBF['x']()
            sage: del pol
        """
        acb_poly_clear(self._poly)

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
            sage: Polynomial_complex_arb(Pol, (CBF(i), 1))                              # needs sage.symbolic
            x + I
            sage: Polynomial_complex_arb(Pol, polygen(QQ,'y')+2)
            x + 2.000000000000000
            sage: Polynomial_complex_arb(Pol, QQ['x'](0))
            0
            sage: Polynomial_complex_arb(Pol, {10: pi})                                 # needs sage.symbolic
            ([3.141592653589793 +/- ...e-16])*x^10
            sage: Polynomial_complex_arb(Pol, pi)                                       # needs sage.symbolic
            [3.141592653589793 +/- ...e-16]
        """
        cdef ComplexBall ball
        cdef Polynomial pol
        cdef list lst
        cdef tuple tpl
        cdef dict dct
        cdef long length, i

        Polynomial.__init__(self, parent, is_gen=is_gen)

        if is_gen:
            acb_poly_set_coeff_si(self._poly, 1, 1)
        elif x is None:
            acb_poly_zero(self._poly)
        elif isinstance(x, Polynomial_complex_arb):
            acb_poly_set(self._poly, (<Polynomial_complex_arb> x)._poly)
        elif isinstance(x, ComplexBall):
            acb_poly_set_coeff_acb(self._poly, 0, (<ComplexBall> x).value)
        else:
            Coeff = parent.base_ring()
            if isinstance(x, list):
                lst = <list> x
                length = len(lst)
                sig_on()
                acb_poly_fit_length(self._poly, length)
                sig_off()
                for i in range(length):
                    ball = Coeff(lst[i])
                    acb_poly_set_coeff_acb(self._poly, i, ball.value)
            elif isinstance(x, tuple):
                tpl = <tuple> x
                length = len(tpl)
                sig_on()
                acb_poly_fit_length(self._poly, length)
                sig_off()
                for i in range(length):
                    ball = Coeff(tpl[i])
                    acb_poly_set_coeff_acb(self._poly, i, ball.value)
            elif isinstance(x, Polynomial):
                pol = <Polynomial> x
                length = pol.degree() + 1
                sig_on()
                acb_poly_fit_length(self._poly, length)
                sig_off()
                for i in range(length):
                    ball = Coeff(pol.get_unsafe(i))
                    acb_poly_set_coeff_acb(self._poly, i, ball.value)
            elif isinstance(x, dict):
                dct = <dict> x
                if len(dct) == 0:
                    acb_poly_zero(self._poly)
                else:
                    length = max(int(i) for i in dct) + 1
                    sig_on()
                    acb_poly_fit_length(self._poly, length)
                    sig_off()
                    for i, c in dct.iteritems():
                        ball = Coeff(c)
                        acb_poly_set_coeff_acb(self._poly, i, ball.value)
            else:
                ball = Coeff(x)
                acb_poly_set_coeff_acb(self._poly, 0, ball.value)

    def __reduce__(self):
        r"""
        Serialize a polynomial for pickling.

        TESTS::

            sage: # needs sage.symbolic
            sage: Pol.<x> = ComplexBallField(42)[]
            sage: pol = (x + i)/3
            sage: pol2 = loads(dumps(pol))
            sage: pol.degree() == pol2.degree()
            True
            sage: all(a.identical(b) for (a, b) in zip(pol, pol2))
            True
        """
        return (self.__class__,
               (self.parent(), self.list(), False, self.is_gen()))

    # Access

    def degree(self):
        r"""
        Return the (apparent) degree of this polynomial.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x^2 + 1).degree()
            2
            sage: pol = (x/3 + 1) - x/3; pol
            ([+/- ...e-16])*x + 1.000000000000000
            sage: pol.degree()
            1
            sage: Pol([1, 0, 0, 0]).degree()
            0
        """
        return smallInteger(acb_poly_degree(self._poly))

    cdef get_unsafe(self, Py_ssize_t n):
        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        res._parent = self._parent._base
        acb_poly_get_coeff_acb(res.value, self._poly, n)
        return res

    cpdef list list(self, bint copy=True):
        r"""
        Return the coefficient list of this polynomial.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x^2/3).list()
            [0, 0, [0.3333333333333333 +/- ...e-17]]
            sage: Pol(0).list()
            []
            sage: Pol([0, 1, RBF(0, rad=.1), 0]).list()
            [0, 1.000000000000000, [+/- 0.101]]
        """
        cdef unsigned long length = acb_poly_length(self._poly)
        return [self.get_unsafe(n) for n in range(length)]

    def __bool__(self):
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
        return acb_poly_length(self._poly)

    # Ring and Euclidean arithmetic

    cpdef _add_(self, other):
        r"""
        Return the sum of two polynomials.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x + 1) + (x/3 - 2)
            ([1.333333333333333 +/- ...e-16])*x - 1.000000000000000
        """
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_add(
                res._poly,
                self._poly,
                (<Polynomial_complex_arb> other)._poly,
                prec(self))
        sig_off()
        return res

    cpdef _neg_(self):
        r"""
        Return the opposite of this polynomial.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: -(x/3 - 2)
            ([-0.3333333333333333 +/- ...e-17])*x + 2.000000000000000
        """
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_neg(res._poly, self._poly)
        sig_off()
        return res

    cpdef _sub_(self, other):
        r"""
        Return the difference of two polynomials.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x + 1) - (x/3 - 2)
            ([0.666666666666667 +/- ...e-16])*x + 3.000000000000000
        """
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_sub(
                res._poly,
                self._poly,
                (<Polynomial_complex_arb> other)._poly,
                prec(self))
        sig_off()
        return res

    cpdef _mul_(self, other):
        r"""
        Return the product of two polynomials.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x + 1)*(x/3 - 2)
            ([0.3333333333333333 +/- ...e-17])*x^2
            + ([-1.666666666666667 +/- ...e-16])*x - 2.000000000000000
        """
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_mul(
                res._poly,
                self._poly,
                (<Polynomial_complex_arb> other)._poly,
                prec(self))
        sig_off()
        return res

    cpdef _lmul_(self, Element a):
        r"""
        TESTS::

            sage: Pol.<x> = CBF[]
            sage: (x + 1)._lmul_(CBF(3))
            3.000000000000000*x + 3.000000000000000
            sage: (1 + x)*(1/3)
            ([0.3333333333333333 +/- ...e-17])*x + [0.3333333333333333 +/- ...e-17]
            sage: (1 + x)*GF(2)(1)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s)...
        """
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_scalar_mul(res._poly, self._poly, (<ComplexBall> a).value, prec(self))
        sig_off()
        return res

    cpdef _rmul_(self, Element a):
        r"""
        TESTS::

            sage: Pol.<x> = CBF[]
            sage: (x + 1)._rmul_(CBF(3))
            3.000000000000000*x + 3.000000000000000
            sage: (1/3)*(1 + x)
            ([0.3333333333333333 +/- ...e-17])*x + [0.3333333333333333 +/- ...e-17]
        """
        return self._lmul_(a)

    @coerce_binop
    def quo_rem(self, divisor):
        r"""
        Compute the Euclidean division of this ball polynomial by ``divisor``.

        Raises a :exc:`ZeroDivisionError` when the divisor is zero or its leading
        coefficient contains zero. Returns a pair (quotient, remainder)
        otherwise.

        EXAMPLES::

            sage: Pol.<x> = CBF[]

            sage: (x^3/7 - CBF(i)).quo_rem(x + CBF(pi))                                 # needs sage.symbolic
            (([0.1428571428571428 +/- ...e-17])*x^2
               + ([-0.448798950512828 +/- ...e-16])*x
               + [1.409943485869908 +/- ...e-16],
             [-4.42946809718569 +/- ...e-15] - I)

            sage: Pol(0).quo_rem(x + 1)
            (0, 0)

            sage: (x + 1).quo_rem(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ('cannot divide by this polynomial', 0)

            sage: div = (x^2/3 + x + 1) - x^2/3; div
            ([+/- ...e-16])*x^2 + x + 1.000000000000000
            sage: (x + 1).quo_rem(div)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ('cannot divide by this polynomial',
            ([+/- ...e-16])*x^2 + x + 1.000000000000000)
        """
        cdef Polynomial_complex_arb div = <Polynomial_complex_arb> divisor
        cdef Polynomial_complex_arb quo = self._new()
        cdef Polynomial_complex_arb rem = self._new()
        sig_on()
        cdef bint success = acb_poly_divrem(quo._poly, rem._poly, self._poly,
                div._poly, prec(self))
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
        acb_poly_set(res._poly, self._poly)
        acb_poly_truncate(res._poly, n)
        sig_off()
        return res

    cdef _inplace_truncate(self, long n):
        if n < 0:
            n = 0
        acb_poly_truncate(self._poly, n)
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
        acb_poly_shift_left(res._poly, self._poly, n)
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
        acb_poly_shift_right(res._poly, self._poly, n)
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
        acb_poly_mullow(res._poly, self._poly, my_other._poly, n, prec(self))
        sig_off()
        return res

    cpdef Polynomial inverse_series_trunc(self, long n):
        r"""
        Return the power series expansion at 0 of the inverse of this
        polynomial, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (1 - x/3).inverse_series_trunc(3)
            ([0.1111111111111111 +/- ...e-17])*x^2 + ([0.3333333333333333 +/- ...e-17])*x + 1.000000000000000
            sage: x.inverse_series_trunc(1)
            nan
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
        acb_poly_inv_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    cpdef Polynomial _power_trunc(self, unsigned long expo, long n):
        r"""
        Return a power of this polynomial, truncated before degree `n`.

        INPUT:

        - ``expo`` -- nonnegative integer exponent
        - ``n`` -- truncation order

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (x^2 + 1)._power_trunc(10^9, 3)
            1000000000.000000*x^2 + 1.000000000000000
            sage: (x^2 + 1)._power_trunc(10^20, 0)
            Traceback (most recent call last):
                ...
            OverflowError: ... int too large to convert...

        TESTS::

            sage: (x^2 + 1)._power_trunc(10, -3)
            0
            sage: (x^2 + 1)._power_trunc(-1, 0)
            Traceback (most recent call last):
            ...
            OverflowError: can...t convert negative value to unsigned long
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_pow_ui_trunc_binexp(res._poly, self._poly, expo, n, prec(self))
        sig_off()
        return res

    def _log_series(self, long n):
        r"""
        Return the power series expansion at 0 of the logarithm of this
        polynomial, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (1 + x/3)._log_series(3)
            ([-0.0555555555555555 +/- ...e-17])*x^2 + ([0.3333333333333333 +/- ...e-17])*x
            sage: (-1 + x)._log_series(3)
            -0.5000000000000000*x^2 - x + [3.141592653589793 +/- ...e-16]*I

        An example where the constant term crosses the branch cut of the
        logarithm::

            sage: pol = CBF(-1, RBF(0, rad=.01)) + x; pol
            x - 1.000000000000000 + [+/- 0.0101]*I
            sage: pol._log_series(2)
            ([-1.000 +/- ...e-4] + [+/- 0.0101]*I)*x + [+/- ...e-5] + [+/- 3.15]*I

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
        acb_poly_log_series(res._poly, self._poly, n, prec(self))
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
            ([+/- ...e-17])*x^2 + ([0.3333333333333333 +/- ...e-17])*x + 1.000000000000000
            sage: (CBF(0, pi) + x)._exp_series(4)                                       # needs sage.symbolic
            ([-0.166...] + [+/- ...]*I)*x^3 + ([-0.500...] + [+/- ...]*I)*x^2
            + ([-1.000...] + [+/- ...]*I)*x + [-1.000...] + [+/- ...]*I
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_exp_series(res._poly, self._poly, n, prec(self))
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
            ([+/- ...e-3] + [+/- 0.501]*I)*x + [+/- ...e-3] + [+/- 1.01]*I
            sage: x._sqrt_series(2)
            (nan + nan*I)*x
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_sqrt_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    def _gamma_series(self, long n):
        r"""
        Return the series expansion of the gamma function composed
        with this polynomial, truncated before degree ``n``.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (1 + x)._gamma_series(3)
            ([0.98905599532797...])*x^2 + ([-0.57721566490153...])*x + 1.000000000000000
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_gamma_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    def _lgamma_series(self, long n):
        r"""
        Return the series expansion of the log-gamma function composed
        with this polynomial, truncated before degree ``n``.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (1000 + x)._lgamma_series(3)
            ([0.00050025008333331...])*x^2 + ([6.9072551956488...])*x + [5905.2204232091...]
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_lgamma_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    def _rgamma_series(self, long n):
        r"""
        Return the series expansion of the reciprocal gamma function composed
        with this polynomial, truncated before degree ``n``.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (-1 + x)._rgamma_series(4)
            ([1.23309373642178...])*x^3 + ([0.422784335098467...])*x^2 - x
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_rgamma_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    def _lambert_w_series(self, long n, branch=0):
        r"""
        Return the series expansion of the specified branch of the Lambert W
        function composed with this polynomial, truncated before degree ``n``.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (1 + x)._lambert_w_series(3)
            ([-0.10727032...])*x^2 + ([0.36189625...])*x + [0.56714329...]
            sage: (CBF(1, 1) + x)._lambert_w_series(2)
            ([0.26651990...] + [-0.15238505...]*I)*x + [0.65696606...] + [0.32545033...]*I
            sage: (1 + x)._lambert_w_series(2, branch=3)
            ([1.00625557...] + [0.05775573...]*I)*x + [-2.85358175...] + [17.1135355...]*I
            sage: (1 + x)._lambert_w_series(2, branch=-3)
            ([1.00625557...] + [-0.05775573...]*I)*x + [-2.85358175...] + [-17.1135355...]*I
            sage: (1 + x)._lambert_w_series(2, branch=2^100)
            ([1.00000000...] + [1.25551112...]*I)*x + [-71.1525951...] + [7.96488362...]*I
        """
        cdef fmpz_t _branch
        fmpz_init(_branch)
        fmpz_set_mpz(_branch, (<Integer> Integer(branch)).value)
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        acb_poly_lambertw_series(res._poly, self._poly, _branch, 0, n, prec(self))
        sig_off()
        fmpz_clear(_branch)
        return res

    def _zeta_series(self, long n, a=1, deflate=False):
        r"""
        Return the series expansion of the Hurwitz zeta function composed
        with this polynomial, truncated before degree ``n``.

        For ``a = 1``, this computes the usual Riemann zeta function.

        If ``deflate`` is True, evaluate ζ(s,a) + 1/(1-s), see the FLINT
        documentation for details.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: (CBF(1/2, 1) + x)._zeta_series(2)
            ([0.55898247...] + [-0.64880821...]*I)*x + [0.14393642...] + [-0.72209974...]*I
            sage: (1/2 + x^2)._zeta_series(3, a=1/3)
            ([-2.13199508...])*x^2 + [-0.11808332...]
            sage: (1 + x)._zeta_series(2, deflate=True)
            ([0.07281584...])*x + [0.57721566...]
        """
        if n < 0:
            n = 0
        cdef ComplexBall _a = <ComplexBall> (self._parent._base.coerce(a))
        cdef Polynomial_complex_arb res = self._new()
        sig_on()
        acb_poly_zeta_series(res._poly, self._poly, _a.value, deflate, n, prec(self))
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
            ([1.666666666666667 +/- ...e-16])*x + 2.000000000000000
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
        cdef acb_poly_t self_ts, other_ts
        cdef acb_ptr cc
        if acb_poly_length(other1._poly) > 0:
            cc = acb_poly_get_coeff_ptr(other1._poly, 0)
            if not acb_is_zero(cc):
                sig_on()
                try:
                    acb_poly_init(self_ts)
                    acb_poly_init(other_ts)
                    acb_poly_taylor_shift(self_ts, self._poly, cc, prec(self))
                    acb_poly_set(other_ts, other1._poly)
                    acb_zero(acb_poly_get_coeff_ptr(other_ts, 0))
                    acb_poly_compose_series(res._poly, self_ts, other_ts, n, prec(self))
                finally:
                    acb_poly_clear(other_ts)
                    acb_poly_clear(self_ts)
                    sig_off()
                return res
        sig_on()
        acb_poly_compose_series(res._poly, self._poly, other1._poly, n, prec(self))
        sig_off()
        return res

    def revert_series(self, long n):
        r"""
        Return a polynomial ``f`` such that
        ``f(self(x)) = self(f(x)) = x mod x^n``.

        EXAMPLES::

            sage: Pol.<x> = CBF[]

            sage: (2*x).revert_series(5)
            0.5000000000000000*x

            sage: (x + x^3/6 + x^5/120).revert_series(6)
            ([0.075000000000000 +/- ...e-17])*x^5 + ([-0.166666666666667 +/- ...e-16])*x^3 + x

            sage: (1 + x).revert_series(6)
            Traceback (most recent call last):
            ...
            ValueError: the constant coefficient must be zero

            sage: (x^2).revert_series(6)
            Traceback (most recent call last):
            ...
            ValueError: the linear term must be nonzero
        """
        cdef Polynomial_complex_arb res = self._new()
        if n < 0:
            n = 0
        if not acb_is_zero(acb_poly_get_coeff_ptr(self._poly, 0)):
            raise ValueError("the constant coefficient must be zero")
        if acb_contains_zero(acb_poly_get_coeff_ptr(self._poly, 1)):
            raise ValueError("the linear term must be nonzero")
        sig_on()
        acb_poly_revert_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    # Evaluation

    def __call__(self, *x, **kwds):
        r"""
        Evaluate this polynomial.

        EXAMPLES::

            sage: Pol.<x> = CBF[]
            sage: pol = x^2 - 1
            sage: pol(CBF(pi))                                                          # needs sage.symbolic
            [8.86960440108936 +/- ...e-15]
            sage: pol(x^3 + 1)
            x^6 + 2.000000000000000*x^3
            sage: pol(matrix([[1,2],[3,4]]))
            [6.000000000000000 10.00000000000000]
            [15.00000000000000 21.00000000000000]

        TESTS::

            sage: P.<x> = CBF[]
            sage: Q.<y> = CBF[]
            sage: x(y)
            y
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
                acb_poly_evaluate(ball.value, self._poly,
                        (<ComplexBall> point).value, prec(self))
                sig_off()
                return ball
            elif isinstance(point, Polynomial_complex_arb):
                poly = (<Polynomial_complex_arb> point)._new()
                sig_on()
                acb_poly_compose(poly._poly, self._poly,
                        (<Polynomial_complex_arb> point)._poly, prec(self))
                sig_off()
                return poly
            # TODO: perhaps add more special cases, e.g. for real ball,
            # integers and rationals
        return Polynomial.__call__(self, *x, **kwds)
