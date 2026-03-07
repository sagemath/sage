r"""
Univariate polynomials over `\RR` with Arb ball coefficients.

This is a binding to the `arb_poly module of FLINT <https://flintlib.org/doc/arb_poly.html>`_;
it may be useful to refer to its documentation for more details.

Parts of the documentation for this module are copied or adapted from Arb's
(now FLINT's) own documentation, licenced (at the time) under the GNU General
Public License version 2, or later.

.. SEEALSO::

    - :mod:`Real balls using Arb <sage.rings.real_arb>`

TESTS:

    sage: type(polygen(RealBallField(140)))
    <class 'sage.rings.polynomial.polynomial_real_arb.Polynomial_real_arb'>
    sage: Pol.<x> = RBF[]
    sage: (x+1/2)^3
    x^3 + 1.500000000000000*x^2 + 0.7500000000000000*x + 0.1250000000000000
"""

from cysignals.signals cimport sig_on, sig_off

from sage.libs.flint.arb cimport *
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.arb_poly_macros cimport arb_poly_get_coeff_ptr
from sage.rings.integer cimport Integer, smallInteger
from sage.rings.real_arb cimport RealBall
from sage.structure.element cimport Element

from sage.structure.element import coerce_binop

cdef inline long prec(Polynomial_real_arb pol) noexcept:
    return pol._parent._base._prec


cdef class Polynomial_real_arb(Polynomial):
    r"""
    Wrapper for `FLINT <https://flintlib.org>`_ polynomials of type
    ``arb_poly_t``

    EXAMPLES::

        sage: Pol.<x> = RBF[]
        sage: type(x)
        <class 'sage.rings.polynomial.polynomial_real_arb.Polynomial_real_arb'>

        sage: Pol(), Pol(1), Pol([0,1,2]), Pol({1: pi, 3: e})                          # needs sage.symbolic  # abs tol 1e-15
        (0,
         1.000000000000000,
         2.000000000000000*x^2 + x,
         [2.718281828459045 +/- 5.41e-16]*x^3 + [3.141592653589793 +/- 3.39e-16]*x)

        sage: Pol("x - 2/3")
        x + [-0.666666666666667 +/- ...e-16]
        sage: Pol(polygen(QQ))
        x

        sage: all(Pol.has_coerce_map_from(P) for P in
        ....:     (QQ['x'], QuadraticField(2), RealBallField(100)))
        True
        sage: any(Pol.has_coerce_map_from(P) for P in
        ....:     (QQ['y'], RR, CC, RDF, CDF, RIF, CIF, RealBallField(20)))
        False
    """

    # Memory management and initialization

    def __cinit__(self):
        r"""
        TESTS::

            sage: RealBallField(2)['y']()
            0
        """
        arb_poly_init(self._poly)

    def __dealloc__(self):
        r"""
        TESTS::

            sage: pol = RBF['x']()
            sage: del pol
        """
        arb_poly_clear(self._poly)

    cdef Polynomial_real_arb _new(self):
        r"""
        Return a new polynomial with the same parent as this one.
        """
        cdef Polynomial_real_arb res = Polynomial_real_arb.__new__(Polynomial_real_arb)
        res._parent = self._parent
        res._is_gen = 0
        return res

    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        r"""
        Initialize this polynomial to the specified value.

        TESTS::

            sage: from sage.rings.polynomial.polynomial_real_arb import Polynomial_real_arb
            sage: Pol = RBF['x']
            sage: Polynomial_real_arb(Pol)
            0
            sage: Polynomial_real_arb(Pol, is_gen=True)
            x
            sage: Polynomial_real_arb(Pol, 42, is_gen=True)
            x
            sage: Polynomial_real_arb(Pol, RBF(1))
            1.000000000000000
            sage: Polynomial_real_arb(Pol, [])
            0
            sage: Polynomial_real_arb(Pol, [0])
            0
            sage: Polynomial_real_arb(Pol, [0, 2, 0])
            2.000000000000000*x
            sage: Polynomial_real_arb(Pol, (1,))
            1.000000000000000
            sage: Polynomial_real_arb(Pol, (RBF(i), 1))                              # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: nonzero imaginary part
            sage: Polynomial_real_arb(Pol, polygen(QQ,'y')+2)
            x + 2.000000000000000
            sage: Polynomial_real_arb(Pol, QQ['x'](0))
            0
            sage: Polynomial_real_arb(Pol, {10: pi})                                 # needs sage.symbolic  # abs tol 1e-15
            [3.141592653589793 +/- 3.39e-16]*x^10
            sage: Polynomial_real_arb(Pol, pi)                                       # needs sage.symbolic
            [3.141592653589793 +/- ...e-16]
        """
        cdef RealBall ball
        cdef Polynomial pol
        cdef list lst
        cdef tuple tpl
        cdef dict dct
        cdef long length, i

        Polynomial.__init__(self, parent, is_gen=is_gen)

        if is_gen:
            arb_poly_set_coeff_si(self._poly, 1, 1)
        elif x is None:
            arb_poly_zero(self._poly)
        elif isinstance(x, Polynomial_real_arb):
            arb_poly_set(self._poly, (<Polynomial_real_arb> x)._poly)
        elif isinstance(x, RealBall):
            arb_poly_set_coeff_arb(self._poly, 0, (<RealBall> x).value)
        else:
            Coeff = parent.base_ring()
            if isinstance(x, list):
                lst = <list> x
                length = len(lst)
                sig_on()
                arb_poly_fit_length(self._poly, length)
                sig_off()
                for i in range(length):
                    ball = Coeff(lst[i])
                    arb_poly_set_coeff_arb(self._poly, i, ball.value)
            elif isinstance(x, tuple):
                tpl = <tuple> x
                length = len(tpl)
                sig_on()
                arb_poly_fit_length(self._poly, length)
                sig_off()
                for i in range(length):
                    ball = Coeff(tpl[i])
                    arb_poly_set_coeff_arb(self._poly, i, ball.value)
            elif isinstance(x, Polynomial):
                pol = <Polynomial> x
                length = pol.degree() + 1
                sig_on()
                arb_poly_fit_length(self._poly, length)
                sig_off()
                for i in range(length):
                    ball = Coeff(pol.get_unsafe(i))
                    arb_poly_set_coeff_arb(self._poly, i, ball.value)
            elif isinstance(x, dict):
                dct = <dict> x
                if not dct:
                    arb_poly_zero(self._poly)
                else:
                    length = max(int(i) for i in dct) + 1
                    sig_on()
                    arb_poly_fit_length(self._poly, length)
                    sig_off()
                    for i, c in dct.items():
                        ball = Coeff(c)
                        arb_poly_set_coeff_arb(self._poly, i, ball.value)
            else:
                ball = Coeff(x)
                arb_poly_set_coeff_arb(self._poly, 0, ball.value)

    def __reduce__(self):
        r"""
        Serialize a polynomial for pickling.

        TESTS::

            sage: # needs sage.symbolic
            sage: Pol.<x> = RealBallField(42)[]
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

            sage: Pol.<x> = RBF[]
            sage: (x^2 + 1).degree()
            2
            sage: pol = (x/3 + 1) - x/3; pol                                          # abs tol 1e-15
            [+/- 1.12e-16]*x + 1.000000000000000
            sage: pol.degree()
            1
            sage: Pol([1, 0, 0, 0]).degree()
            0
        """
        return smallInteger(arb_poly_degree(self._poly))

    cdef get_unsafe(self, Py_ssize_t n):
        cdef RealBall res = RealBall.__new__(RealBall)
        res._parent = self._parent._base
        arb_poly_get_coeff_arb(res.value, self._poly, n)
        return res

    cpdef list list(self, bint copy=True):
        r"""
        Return the coefficient list of this polynomial.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (x^2/3).list()
            [0, 0, [0.3333333333333333 +/- ...e-17]]
            sage: Pol(0).list()
            []
            sage: Pol([0, 1, RBF(0, rad=.1), 0]).list()
            [0, 1.000000000000000, [+/- 0.101]]
        """
        cdef unsigned long length = arb_poly_length(self._poly)
        return [self.get_unsafe(n) for n in range(length)]

    def __bool__(self):
        r"""
        Return ``False`` if this polynomial is exactly zero, ``True`` otherwise.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: bool(Pol(0))
            False
            sage: z = Pol(1/3) - 1/3
            sage: bool(z)
            True
        """
        return arb_poly_length(self._poly)

    # Ring and Euclidean arithmetic

    cpdef _add_(self, other):
        r"""
        Return the sum of two polynomials.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (x + 1) + (x/3 - 2)                                                 # abs tol 1e-15
            [1.333333333333333 +/- 5.37e-16]*x - 1.000000000000000
        """
        cdef Polynomial_real_arb res = self._new()
        sig_on()
        arb_poly_add(
                res._poly,
                self._poly,
                (<Polynomial_real_arb> other)._poly,
                prec(self))
        sig_off()
        return res

    cpdef _neg_(self):
        r"""
        Return the opposite of this polynomial.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: -(x/3 - 2)                                                         # abs tol 1e-15
            [-0.3333333333333333 +/- 7.04e-17]*x + 2.000000000000000
        """
        cdef Polynomial_real_arb res = self._new()
        sig_on()
        arb_poly_neg(res._poly, self._poly)
        sig_off()
        return res

    cpdef _sub_(self, other):
        r"""
        Return the difference of two polynomials.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (x + 1) - (x/3 - 2)                                                 # abs tol 1e-15
            [0.666666666666667 +/- 5.37e-16]*x + 3.000000000000000
        """
        cdef Polynomial_real_arb res = self._new()
        sig_on()
        arb_poly_sub(
                res._poly,
                self._poly,
                (<Polynomial_real_arb> other)._poly,
                prec(self))
        sig_off()
        return res

    cpdef _mul_(self, other):
        r"""
        Return the product of two polynomials.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (x + 1)*(x/3 - 2)                                                   # abs tol 1e-15
            [0.3333333333333333 +/- 7.04e-17]*x^2 + [-1.666666666666667 +/- 7.59e-16]*x - 2.000000000000000
        """
        cdef Polynomial_real_arb res = self._new()
        sig_on()
        arb_poly_mul(
                res._poly,
                self._poly,
                (<Polynomial_real_arb> other)._poly,
                prec(self))
        sig_off()
        return res

    cpdef _lmul_(self, Element a):
        r"""
        TESTS::

            sage: Pol.<x> = RBF[]
            sage: (x + 1)._lmul_(RBF(3))
            3.000000000000000*x + 3.000000000000000
            sage: (1 + x)*(1/3)                                                       # abs tol 1e-15
            [0.3333333333333333 +/- 7.04e-17]*x + [0.3333333333333333 +/- 7.04e-17]
            sage: (1 + x)*GF(2)(1)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s)...
        """
        cdef Polynomial_real_arb res = self._new()
        sig_on()
        arb_poly_scalar_mul(res._poly, self._poly, (<RealBall> a).value, prec(self))
        sig_off()
        return res

    cpdef _rmul_(self, Element a):
        r"""
        TESTS::

            sage: Pol.<x> = RBF[]
            sage: (x + 1)._rmul_(RBF(3))
            3.000000000000000*x + 3.000000000000000
            sage: (1/3)*(1 + x)                                                       # abs tol 1e-15
            [0.3333333333333333 +/- 7.04e-17]*x + [0.3333333333333333 +/- 7.04e-17]
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

            sage: Pol.<x> = RBF[]

            sage: (x^3/7 - 2).quo_rem(x + 3)                                           # abs tol 1e-14
            ([0.1428571428571428 +/- 7.70e-17]*x^2 + [-0.428571428571429 +/- 5.36e-16]*x + [1.285714285714286 +/- 8.85e-16],
             [-5.85714285714286 +/- 4.66e-15])

            sage: Pol(0).quo_rem(x + 1)
            (0, 0)

            sage: (x + 1).quo_rem(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ('cannot divide by this polynomial', 0)

            sage: div = (x^2/3 + x + 1) - x^2/3; div                                  # abs tol 1e-15
            [+/- 1.12e-16]*x^2 + x + 1.000000000000000
            sage: (x + 1).quo_rem(div)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: ('cannot divide by this polynomial', [+/- 1.12e-16]*x^2 + x + 1.000000000000000)
        """
        cdef Polynomial_real_arb div = <Polynomial_real_arb> divisor
        cdef Polynomial_real_arb quo = self._new()
        cdef Polynomial_real_arb rem = self._new()
        sig_on()
        cdef bint success = arb_poly_divrem(quo._poly, rem._poly, self._poly,
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

            sage: pol = RBF['x'](range(1,5)); pol
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
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        arb_poly_set(res._poly, self._poly)
        arb_poly_truncate(res._poly, n)
        sig_off()
        return res

    cdef _inplace_truncate(self, long n):
        if n < 0:
            n = 0
        arb_poly_truncate(self._poly, n)
        return self

    def __lshift__(val, n):
        r"""
        Shift ``val`` to the left, i.e. multiply it by `x^n`, throwing away
        coefficients if `n < 0`.

        EXAMPLES::

            sage: pol = RBF['x'](range(1,5)); pol
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
        if not isinstance(val, Polynomial_real_arb):
            raise TypeError("unsupported operand type(s) for <<: '{}' and '{}'"
                            .format(type(val).__name__, type(n).__name__))
        if n < 0:
            return val.__rshift__(-n)
        cdef Polynomial_real_arb self = (<Polynomial_real_arb> val)
        cdef Polynomial_real_arb res = self._new()
        sig_on()
        arb_poly_shift_left(res._poly, self._poly, n)
        sig_off()
        return res

    def __rshift__(val, n):
        r"""
        Shift ``val`` to the left, i.e. divide it by `x^n`, throwing away
        coefficients if `n > 0`.

        EXAMPLES::

            sage: pol = RBF['x'](range(1,5)); pol
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
        if not isinstance(val, Polynomial_real_arb):
            raise TypeError("unsupported operand type(s) for <<: '{}' and '{}'"
                            .format(type(val).__name__, type(n).__name__))
        if n < 0:
            return val.__lshift__(-n)
        cdef Polynomial_real_arb self = (<Polynomial_real_arb> val)
        cdef Polynomial_real_arb res = self._new()
        sig_on()
        arb_poly_shift_right(res._poly, self._poly, n)
        sig_off()
        return res

    # Truncated and power series arithmetic

    cpdef Polynomial _mul_trunc_(self, Polynomial other, long n):
        r"""
        Return the product of ``self`` and ``other``, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
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
        cdef Polynomial_real_arb my_other = <Polynomial_real_arb> other
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        arb_poly_mullow(res._poly, self._poly, my_other._poly, n, prec(self))
        sig_off()
        return res

    cpdef Polynomial inverse_series_trunc(self, long n):
        r"""
        Return the power series expansion at 0 of the inverse of this
        polynomial, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (1 - x/3).inverse_series_trunc(3)                                   # abs tol 1e-15
            [0.1111111111111111 +/- 5.99e-17]*x^2 + [0.3333333333333333 +/- 7.04e-17]*x + 1.000000000000000
            sage: x.inverse_series_trunc(1)
            nan
            sage: Pol(0).inverse_series_trunc(2)
            nan*x + nan

        TESTS::

            sage: Pol(0).inverse_series_trunc(-1)
            0
        """
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        arb_poly_inv_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    cpdef Polynomial _power_trunc(self, unsigned long expo, long n):
        r"""
        Return a power of this polynomial, truncated before degree `n`.

        INPUT:

        - ``expo`` -- nonnegative integer exponent
        - ``n`` -- truncation order

        EXAMPLES::

            sage: Pol.<x> = RBF[]
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
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        arb_poly_pow_ui_trunc_binexp(res._poly, self._poly, expo, n, prec(self))
        sig_off()
        return res

    def _log_series(self, long n):
        r"""
        Return the power series expansion at 0 of the logarithm of this
        polynomial, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (1 + x/3)._log_series(3)                                            # abs tol 1e-15
            [-0.0555555555555555 +/- 7.10e-17]*x^2 + [0.3333333333333333 +/- 7.04e-17]*x

        Some cases where the result is not defined::
            sage: (-1 + x)._log_series(3)  # the result's constant term is pi*I, not real
            nan*x^2 + nan*x + nan
            sage: x._log_series(1)
            nan
            sage: Pol(0)._log_series(1)
            nan
        """
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        arb_poly_log_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    def _exp_series(self, long n):
        r"""
        Return the power series expansion at 0 of the exponential of this
        polynomial, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: x._exp_series(3)
            0.5000000000000000*x^2 + x + 1.000000000000000
            sage: (1 + x/3)._log_series(3)._exp_series(3)                              # abs tol 1e-15
            [+/- 4.79e-17]*x^2 + [0.3333333333333333 +/- 7.04e-17]*x + 1.000000000000000
        """
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        arb_poly_exp_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    def _sqrt_series(self, long n):
        r"""
        Return the power series expansion at 0 of the square root of this
        polynomial, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (1 + x)._sqrt_series(3)
            -0.1250000000000000*x^2 + 0.5000000000000000*x + 1.000000000000000
            sage: x._sqrt_series(2)
            nan*x
        """
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        arb_poly_sqrt_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    def _gamma_series(self, long n):
        r"""
        Return the series expansion of the gamma function composed
        with this polynomial, truncated before degree ``n``.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (1 + x)._gamma_series(3)                                            # abs tol 1e-12
            [0.989055995327973 +/- 6.12e-16]*x^2 + [-0.577215664901533 +/- 3.58e-16]*x + 1.000000000000000
        """
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        arb_poly_gamma_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    def _lgamma_series(self, long n):
        r"""
        Return the series expansion of the log-gamma function composed
        with this polynomial, truncated before degree ``n``.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (1000 + x)._lgamma_series(3)                                        # abs tol 1e-15
            [0.000500250083333317 +/- 4.84e-19]*x^2 + [6.90725519564881 +/- 2.90e-15]*x + [5905.220423209181 +/- 2.49e-13]
        """
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        arb_poly_lgamma_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    def _rgamma_series(self, long n):
        r"""
        Return the series expansion of the reciprocal gamma function composed
        with this polynomial, truncated before degree ``n``.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (-1 + x)._rgamma_series(4)                                          # abs tol 1e-15
            [1.233093736421787 +/- 5.71e-16]*x^3 + [0.4227843350984671 +/- 9.09e-17]*x^2 - x
        """
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        arb_poly_rgamma_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    def _lambert_w_series(self, long n, branch=0):
        r"""
        Return the series expansion of the specified branch of the Lambert W
        function composed with this polynomial, truncated before degree ``n``.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (1 + x)._lambert_w_series(3)  # abs tol 1e-14
            [-0.107270323141072 +/- 7.79e-16]*x^2 + [0.361896256634889 +/- 2.99e-16]*x + [0.567143290409784 +/- 2.72e-16]
            sage: (-1/4 + x)._lambert_w_series(2, branch=-1)  # abs tol 1e-14
            [-7.46833129610253 +/- 3.12e-15]*x + [-2.153292364110349 +/- 8.59e-16]
        """
        if branch == 0:
            flags = 0
        elif branch == -1:
            flags = 1
        else:
            raise ValueError("for other branches, use polynomials over CBF")
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        sig_on()
        arb_poly_lambertw_series(res._poly, self._poly, flags, n, prec(self))
        sig_off()
        return res

    def _zeta_series(self, long n, a=1, deflate=False):
        r"""
        Return the series expansion of the Hurwitz zeta function composed
        with this polynomial, truncated before degree ``n``.

        For ``a = 1``, this computes the usual Riemann zeta function.

        If ``deflate`` is True, evaluate ζ(s,a) + 1/(1-s), see the FLINT
        documentation for details.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: (1/2 + x^2)._zeta_series(3, a=1/3)                                  # abs tol 1e-15
            [-2.13199508553675 +/- 2.90e-15]*x^2 + [-0.118083327934222 +/- 5.54e-16]
            sage: (1 + x)._zeta_series(2, deflate=True)                               # abs tol 1e-15
            [0.0728158454836767 +/- 2.97e-17]*x + [0.5772156649015329 +/- 4.09e-17]
        """
        if n < 0:
            n = 0
        cdef RealBall _a = <RealBall> (self._parent._base.coerce(a))
        cdef Polynomial_real_arb res = self._new()
        sig_on()
        arb_poly_zeta_series(res._poly, self._poly, _a.value, deflate, n, prec(self))
        sig_off()
        return res

    def compose_trunc(self, Polynomial other, long n):
        r"""
        Return the composition of ``self`` and ``other``, truncated before degree `n`.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: Pol.<x> = RBF[]
            sage: pol = x*(x-1)^2
            sage: pol.compose_trunc(x + x^2, 4)
            -3.000000000000000*x^3 - x^2 + x
            sage: pol.compose_trunc(1 + x, 4)
            x^3 + x^2
            sage: pol.compose_trunc(2 + x/3, 2)                                       # abs tol 1e-15
            [1.666666666666667 +/- 9.81e-16]*x + 2.000000000000000
            sage: pol.compose_trunc(2 + x/3, 0)
            0
            sage: pol.compose_trunc(2 + x/3, -1)
            0
        """
        if n < 0:
            n = 0
        if not isinstance(other, Polynomial_real_arb):
            return self(other).truncate(n)
        cdef Polynomial_real_arb other1 = <Polynomial_real_arb> other
        cdef Polynomial_real_arb res = self._new()
        cdef arb_poly_t self_ts, other_ts
        cdef arb_ptr cc
        if arb_poly_length(other1._poly) > 0:
            cc = arb_poly_get_coeff_ptr(other1._poly, 0)
            if not arb_is_zero(cc):
                sig_on()
                try:
                    arb_poly_init(self_ts)
                    arb_poly_init(other_ts)
                    arb_poly_taylor_shift(self_ts, self._poly, cc, prec(self))
                    arb_poly_set(other_ts, other1._poly)
                    arb_zero(arb_poly_get_coeff_ptr(other_ts, 0))
                    arb_poly_compose_series(res._poly, self_ts, other_ts, n, prec(self))
                finally:
                    arb_poly_clear(other_ts)
                    arb_poly_clear(self_ts)
                    sig_off()
                return res
        sig_on()
        arb_poly_compose_series(res._poly, self._poly, other1._poly, n, prec(self))
        sig_off()
        return res

    def revert_series(self, long n):
        r"""
        Return a polynomial ``f`` such that
        ``f(self(x)) = self(f(x)) = x mod x^n``.

        EXAMPLES::

            sage: Pol.<x> = RBF[]

            sage: (2*x).revert_series(5)
            0.5000000000000000*x

            sage: (x + x^3/6 + x^5/120).revert_series(6)                              # abs tol 1e-15
            [0.075000000000000 +/- 9.96e-17]*x^5 + [-0.166666666666667 +/- 4.45e-16]*x^3 + x

            sage: (1 + x).revert_series(6)
            Traceback (most recent call last):
            ...
            ValueError: the constant coefficient must be zero

            sage: (x^2).revert_series(6)
            Traceback (most recent call last):
            ...
            ValueError: the linear term must be nonzero
        """
        cdef Polynomial_real_arb res = self._new()
        if n < 0:
            n = 0
        if not arb_is_zero(arb_poly_get_coeff_ptr(self._poly, 0)):
            raise ValueError("the constant coefficient must be zero")
        if arb_contains_zero(arb_poly_get_coeff_ptr(self._poly, 1)):
            raise ValueError("the linear term must be nonzero")
        sig_on()
        arb_poly_revert_series(res._poly, self._poly, n, prec(self))
        sig_off()
        return res

    # Evaluation

    def __call__(self, *x, **kwds):
        r"""
        Evaluate this polynomial.

        EXAMPLES::

            sage: Pol.<x> = RBF[]
            sage: pol = x^2 - 1
            sage: pol(RBF(pi))                                                          # needs sage.symbolic
            [8.86960440108936 +/- ...e-15]
            sage: pol(x^3 + 1)
            x^6 + 2.000000000000000*x^3
            sage: pol(matrix([[1,2],[3,4]]))
            [6.000000000000000 10.00000000000000]
            [15.00000000000000 21.00000000000000]

        TESTS::

            sage: P.<x> = RBF[]
            sage: Q.<y> = RBF[]
            sage: x(y)
            y
        """
        cdef RealBall ball
        cdef Polynomial_real_arb poly
        if len(x) == 1 and not kwds:
            point = x[0]
            if isinstance(point, RealBall):
                # parent of result = base ring of self (not parent of point)
                ball = RealBall.__new__(RealBall)
                ball._parent = self._parent._base
                sig_on()
                arb_poly_evaluate(ball.value, self._poly,
                        (<RealBall> point).value, prec(self))
                sig_off()
                return ball
            elif isinstance(point, Polynomial_real_arb):
                poly = (<Polynomial_real_arb> point)._new()
                sig_on()
                arb_poly_compose(poly._poly, self._poly,
                        (<Polynomial_real_arb> point)._poly, prec(self))
                sig_off()
                return poly
            # TODO: perhaps add more special cases, e.g. for real ball,
            # integers and rationals
        return Polynomial.__call__(self, *x, **kwds)
