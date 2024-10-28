r"""
Dense univariate polynomials over `\RR`, implemented using MPFR

TESTS:

Check that operations with numpy elements work well (see :issue:`18076` and
:issue:`8426`)::

    sage: # needs numpy
    sage: import numpy
    sage: x = polygen(RR)
    sage: x * numpy.int32('1')
    x
    sage: numpy.int32('1') * x
    x
    sage: x * numpy.int64('1')
    x
    sage: numpy.int64('1') * x
    x
    sage: x * numpy.float32('1.5')
    1.50000000000000*x
    sage: numpy.float32('1.5') * x
    1.50000000000000*x
"""

from cysignals.memory cimport check_allocarray, check_reallocarray, sig_free
from cysignals.signals cimport sig_on, sig_off

from cpython.long cimport PyLong_AsLong
from cpython.float cimport PyFloat_AS_DOUBLE

from sage.structure.parent cimport Parent
from sage.rings.polynomial.polynomial_element cimport Polynomial, _dict_to_list
from sage.rings.real_mpfr cimport RealField_class, RealNumber
from sage.rings.integer cimport Integer, smallInteger
from sage.rings.rational cimport Rational

from sage.structure.element cimport Element
from sage.structure.element cimport parent
from sage.structure.element import coerce_binop
from sage.libs.mpfr cimport *

try:
    from sage.libs.pari.all import pari_gen
except ImportError:
    pari_gen = ()


cdef class PolynomialRealDense(Polynomial):
    r"""

    TESTS::

        sage: f = RR['x'].random_element()
        sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
        sage: isinstance(f, PolynomialRealDense)
        True
    """

    cdef Py_ssize_t _degree
    cdef mpfr_t* _coeffs
    cdef RealField_class _base_ring

    def __cinit__(self):
        """
        TESTS::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: PolynomialRealDense(RR['x'])
            0
        """
        self._coeffs = NULL

    def __init__(self, Parent parent, x=0, check=None, bint is_gen=False, construct=None):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: PolynomialRealDense(RR['x'], [1, int(2), RR(3), 4/1, pi])             # needs sage.symbolic
            3.14159265358979*x^4 + 4.00000000000000*x^3 + 3.00000000000000*x^2 + 2.00000000000000*x + 1.00000000000000
            sage: PolynomialRealDense(RR['x'], None)
            0

        TESTS:

        Check that errors and interrupts are handled properly (see :issue:`10100`)::

            sage: a = var('a')                                                          # needs sage.symbolic
            sage: PolynomialRealDense(RR['x'], [1,a])                                   # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate symbolic expression to a numeric value
            sage: R.<x> = SR[]                                                          # needs sage.symbolic
            sage: (x-a).change_ring(RR)                                                 # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: cannot evaluate symbolic expression to a numeric value
            sage: sig_on_count()
            0

        Test that we don't clean up uninitialized coefficients (:issue:`9826`)::

            sage: k.<a> = GF(7^3)                                                       # needs sage.rings.finite_rings
            sage: P.<x> = PolynomialRing(k)                                             # needs sage.rings.finite_rings
            sage: (a*x).complex_roots()                                                 # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'a' to a real number

        Check that :issue:`17190` is fixed::

            sage: RR['x']({})
            0
        """
        Polynomial.__init__(self, parent, is_gen=is_gen)
        self._base_ring = parent._base
        cdef Py_ssize_t i, degree
        cdef int prec = self._base_ring._prec
        cdef mpfr_rnd_t rnd = self._base_ring.rnd
        if x is None:
            self._coeffs = <mpfr_t*>check_allocarray(1, sizeof(mpfr_t)) # degree zero
            mpfr_init2(self._coeffs[0], prec)
            mpfr_set_si(self._coeffs[0], 0, rnd)
            self._normalize()
            return
        if is_gen:
            x = [0, 1]
        elif isinstance(x, (int, float, Integer, Rational, RealNumber)):
            x = [x]
        elif isinstance(x, dict):
            x = _dict_to_list(x, self._base_ring.zero())
        elif isinstance(x, pari_gen):
            x = [self._base_ring(w) for w in x.list()]
        elif not isinstance(x, list):
            try:
                x = list(x)
            except TypeError:  # x is not iterable
                x = [self._base_ring(x)]

        sig_on()
        degree = len(x) - 1
        self._degree = -1
        cdef mpfr_t* coeffs
        coeffs = <mpfr_t*>check_allocarray(degree+1, sizeof(mpfr_t))
        try:  # We might get Python exceptions here
            for i from 0 <= i <= degree:
                mpfr_init2(coeffs[i], prec)
                self._degree += 1
                a = x[i]
                if type(a) is RealNumber:
                    mpfr_set(coeffs[i], (<RealNumber>a).value, rnd)
                elif type(a) is int:
                    mpfr_set_si(coeffs[i], PyLong_AsLong(a), rnd)
                elif type(a) is float:
                    mpfr_set_d(coeffs[i], PyFloat_AS_DOUBLE(a), rnd)
                elif type(a) is Integer:
                    mpfr_set_z(coeffs[i], (<Integer>a).value, rnd)
                elif type(a) is Rational:
                    mpfr_set_q(coeffs[i], (<Rational>a).value, rnd)
                else:
                    a = self._base_ring(a)
                    mpfr_set(coeffs[i], (<RealNumber>a).value, rnd)
        finally:
            sig_off()
        self._coeffs = coeffs
        self._normalize()

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._coeffs != NULL:
            for i from 0 <= i <= self._degree:
                mpfr_clear(self._coeffs[i])
            sig_free(self._coeffs)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [-2, 0, 1])
            sage: loads(dumps(f)) == f
            True
        """
        return make_PolynomialRealDense, (self._parent, self.list())

    cdef _normalize(self):
        """
        Remove all leading 0's.
        """
        cdef Py_ssize_t i
        if self._degree >= 0 and mpfr_zero_p(self._coeffs[self._degree]):
            i = self._degree
            while i >= 0 and mpfr_zero_p(self._coeffs[i]):
                mpfr_clear(self._coeffs[i])
                i -= 1
            self._coeffs = <mpfr_t*>check_reallocarray(self._coeffs, i+1, sizeof(mpfr_t))
            self._degree = i

    cdef get_unsafe(self, Py_ssize_t i):
        """
        Return the `i`-th coefficient of ``self``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], range(5)); f
            4.00000000000000*x^4 + 3.00000000000000*x^3 + 2.00000000000000*x^2 + x
            sage: f[0]
            0.000000000000000
            sage: f[3]
            3.00000000000000
            sage: f[5]
            0.000000000000000

        Test slices::

            sage: R.<x> = RealField(10)[]
            sage: f = (x+1)^5; f
            x^5 + 5.0*x^4 + 10.*x^3 + 10.*x^2 + 5.0*x + 1.0
            sage: f[:3]
            10.*x^2 + 5.0*x + 1.0
        """
        cdef RealNumber r = <RealNumber>RealNumber.__new__(RealNumber, self._base_ring)
        mpfr_set(r.value, self._coeffs[i], self._base_ring.rnd)
        return r

    cdef PolynomialRealDense _new(self, Py_ssize_t degree):
        cdef Py_ssize_t i
        cdef int prec = self._base_ring._prec
        cdef PolynomialRealDense f = <PolynomialRealDense>PolynomialRealDense.__new__(PolynomialRealDense)
        f._parent = self._parent
        f._base_ring = self._base_ring
        f._degree = degree
        if degree >= 0:
            f._coeffs = <mpfr_t*>check_allocarray(degree+1, sizeof(mpfr_t))
            for i from 0 <= i <= degree:
                mpfr_init2(f._coeffs[i], prec)
        return f

    def degree(self):
        """
        Return the degree of the polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [1, 2, 3]); f
            3.00000000000000*x^2 + 2.00000000000000*x + 1.00000000000000
            sage: f.degree()
            2

        TESTS::

            sage: type(f.degree())
            <class 'sage.rings.integer.Integer'>
        """
        return smallInteger(self._degree)

    cpdef Polynomial truncate(self, long n):
        r"""
        Return the polynomial of degree `< n` which is equivalent to ``self``
        modulo `x^n`.

        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RealField(10)['x'], [1, 2, 4, 8])
            sage: f.truncate(3)
            4.0*x^2 + 2.0*x + 1.0
            sage: f.truncate(100)
            8.0*x^3 + 4.0*x^2 + 2.0*x + 1.0
            sage: f.truncate(1)
            1.0
            sage: f.truncate(0)
            0
        """
        if n <= 0:
            return self._new(-1)
        if n > self._degree:
            return self
        cdef PolynomialRealDense f = self._new(n-1)
        cdef Py_ssize_t i
        for i from 0 <= i < n:
            mpfr_set(f._coeffs[i], self._coeffs[i], self._base_ring.rnd)
        f._normalize()
        return f

    def truncate_abs(self, RealNumber bound):
        """
        Truncate all high order coefficients below ``bound``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RealField(10)['x'], [10^-k for k in range(10)])
            sage: f
            1.0e-9*x^9 + 1.0e-8*x^8 + 1.0e-7*x^7 + 1.0e-6*x^6 + 0.000010*x^5
             + 0.00010*x^4 + 0.0010*x^3 + 0.010*x^2 + 0.10*x + 1.0
            sage: f.truncate_abs(0.5e-6)
            1.0e-6*x^6 + 0.000010*x^5 + 0.00010*x^4 + 0.0010*x^3 + 0.010*x^2 + 0.10*x + 1.0
            sage: f.truncate_abs(10.0)
            0
            sage: f.truncate_abs(1e-100) == f
            True
        """
        cdef Py_ssize_t i
        for i from self._degree >= i >= 0:
            if mpfr_cmpabs(self._coeffs[i], bound.value) >= 0:
                return self.truncate(i+1)
        return self._new(-1)

    cpdef shift(self, Py_ssize_t n):
        r"""
        Return this polynomial multiplied by the power `x^n`. If `n`
        is negative, terms below `x^n` will be discarded. Does not
        change this polynomial.

        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [1, 2, 3]); f
            3.00000000000000*x^2 + 2.00000000000000*x + 1.00000000000000
            sage: f.shift(10)
            3.00000000000000*x^12 + 2.00000000000000*x^11 + x^10
            sage: f.shift(-1)
            3.00000000000000*x + 2.00000000000000
            sage: f.shift(-10)
            0

        TESTS::

            sage: f = RR['x'](0)
            sage: f.shift(3).is_zero()
            True
            sage: f.shift(-3).is_zero()
            True
        """
        cdef Py_ssize_t i
        cdef Py_ssize_t nn = 0 if n < 0 else n
        cdef PolynomialRealDense f
        if n == 0 or self._degree < 0:
            return self
        elif self._degree < -n:
            return self._new(-1)
        else:
            f = self._new(self._degree + n)
            for i from 0 <= i < n:
                mpfr_set_ui(f._coeffs[i], 0, self._base_ring.rnd)
            for i from nn <= i <= self._degree + n:
                mpfr_set(f._coeffs[i], self._coeffs[i-n], self._base_ring.rnd)
        return f

    cpdef list list(self, bint copy=True):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [1, 0, -2]); f
            -2.00000000000000*x^2 + 1.00000000000000
            sage: f.list()
            [1.00000000000000, 0.000000000000000, -2.00000000000000]
        """
        cdef RealNumber r
        cdef Py_ssize_t i
        cdef list L = []
        for i from 0 <= i <= self._degree:
            r = <RealNumber>RealNumber(self._base_ring)
            mpfr_set(r.value, self._coeffs[i], self._base_ring.rnd)
            L.append(r)
        return L

    def __neg__(self):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [-2,0,1])
            sage: -f
            -x^2 + 2.00000000000000
        """
        cdef Py_ssize_t i
        cdef mpfr_rnd_t rnd = self._base_ring.rnd
        cdef PolynomialRealDense f = self._new(self._degree)
        for i from 0 <= i <= f._degree:
            mpfr_neg(f._coeffs[i], self._coeffs[i], rnd)
        return f

    cpdef _add_(left, _right):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [-2,0,1]); f
            x^2 - 2.00000000000000
            sage: g = PolynomialRealDense(RR['x'], range(5)); g
            4.00000000000000*x^4 + 3.00000000000000*x^3 + 2.00000000000000*x^2 + x
            sage: f+g
            4.00000000000000*x^4 + 3.00000000000000*x^3 + 3.00000000000000*x^2 + x - 2.00000000000000
            sage: g + f == f + g
            True
            sage: f + (-f)
            0
        """
        cdef Py_ssize_t i
        cdef mpfr_rnd_t rnd = left._base_ring.rnd
        cdef PolynomialRealDense right = _right
        cdef Py_ssize_t min = left._degree if left._degree < right._degree else right._degree
        cdef Py_ssize_t max = left._degree if left._degree > right._degree else right._degree
        cdef PolynomialRealDense f = left._new(max)
        for i from 0 <= i <= min:
            mpfr_add(f._coeffs[i], left._coeffs[i], right._coeffs[i], rnd)
        if left._degree < right._degree:
            for i from min < i <= max:
                mpfr_set(f._coeffs[i], right._coeffs[i], rnd)
        else:
            for i from min < i <= max:
                mpfr_set(f._coeffs[i], left._coeffs[i], rnd)
        f._normalize()
        return f

    cpdef _sub_(left, _right):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [-3,0,1]); f
            x^2 - 3.00000000000000
            sage: g = PolynomialRealDense(RR['x'], range(4)); g
            3.00000000000000*x^3 + 2.00000000000000*x^2 + x
            sage: f-g
            -3.00000000000000*x^3 - x^2 - x - 3.00000000000000
            sage: (f-g) == -(g-f)
            True
        """
        cdef Py_ssize_t i
        cdef mpfr_rnd_t rnd = left._base_ring.rnd
        cdef PolynomialRealDense right = _right
        cdef Py_ssize_t min = left._degree if left._degree < right._degree else right._degree
        cdef Py_ssize_t max = left._degree if left._degree > right._degree else right._degree
        cdef PolynomialRealDense f = left._new(max)
        for i from 0 <= i <= min:
            mpfr_sub(f._coeffs[i], left._coeffs[i], right._coeffs[i], rnd)
        if left._degree < right._degree:
            for i from min < i <= max:
                mpfr_neg(f._coeffs[i], right._coeffs[i], rnd)
        else:
            for i from min < i <= max:
                mpfr_set(f._coeffs[i], left._coeffs[i], rnd)
        f._normalize()
        return f

    cpdef _lmul_(self, Element c):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [-5,0,0,1]); f
            x^3 - 5.00000000000000
            sage: 4.0 * f
            4.00000000000000*x^3 - 20.0000000000000
            sage: f * -0.2
            -0.200000000000000*x^3 + 1.00000000000000
        """
        cdef Py_ssize_t i
        cdef RealNumber a = c
        if mpfr_zero_p(a.value):
            return self._new(-1)
        cdef mpfr_rnd_t rnd = self._base_ring.rnd
        cdef PolynomialRealDense f = self._new(self._degree)
        for i from 0 <= i <= self._degree:
            mpfr_mul(f._coeffs[i], self._coeffs[i], a.value, rnd)
        return f

    cpdef _mul_(left, _right):
        """
        Here we use the naive `O(n^2)` algorithm, as asymptotically faster algorithms such
        as Karatsuba can have very inaccurate results due to intermediate rounding errors.

        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [1e20, 1])
            sage: g = PolynomialRealDense(RR['x'], [1e30, 1])
            sage: f*g
            x^2 + 1.00000000010000e30*x + 1.00000000000000e50
            sage: f._mul_karatsuba(g,0)
            x^2 + 1.00000000000000e50
            sage: f = PolynomialRealDense(RR['x'], range(5))
            sage: g = PolynomialRealDense(RR['x'], range(3))
            sage: f*g
            8.00000000000000*x^6 + 10.0000000000000*x^5 + 7.00000000000000*x^4 + 4.00000000000000*x^3 + x^2
        """
        cdef Py_ssize_t i, j
        cdef mpfr_rnd_t rnd = left._base_ring.rnd
        cdef PolynomialRealDense right = _right
        cdef PolynomialRealDense f
        cdef mpfr_t tmp
        if left._degree < 0 or right._degree < 0:
            f = left._new(-1)
        else:
            f = left._new(left._degree + right._degree)
        sig_on()
        mpfr_init2(tmp, left._base_ring._prec)
        for i from 0 <= i <= f._degree:
            # Yes, we could make this more efficient by initializing with
            # a multiple of left rather than all zeros...
            mpfr_set_ui(f._coeffs[i], 0, rnd)
        for i from 0 <= i <= left._degree:
            for j from 0 <= j <= right._degree:
                mpfr_mul(tmp, left._coeffs[i], right._coeffs[j], rnd)
                mpfr_add(f._coeffs[i+j], f._coeffs[i+j], tmp, rnd)
        mpfr_clear(tmp)
        sig_off()
        return f

    def _derivative(self, var=None):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [pi, 0, 2, 1])                       # needs sage.symbolic
            sage: f.derivative()                                                        # needs sage.symbolic
            3.00000000000000*x^2 + 4.00000000000000*x

        TESTS::

            sage: x, y = var('x,y')                                                     # needs sage.symbolic
            sage: f.derivative(x)                                                       # needs sage.symbolic
            3.00000000000000*x^2 + 4.00000000000000*x
            sage: f.derivative(y)                                                       # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate with respect to y
        """
        if var is not None and var != self._parent.gen():
            raise ValueError("cannot differentiate with respect to {}".format(var))

        cdef mpfr_rnd_t rnd = self._base_ring.rnd
        cdef PolynomialRealDense f = self._new(self._degree-1)
        for i from 0 <= i < self._degree:
            mpfr_mul_ui(f._coeffs[i], self._coeffs[i+1], i+1, rnd)
        return f

    def integral(self):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [3, pi, 1])                          # needs sage.symbolic
            sage: f.integral()                                                          # needs sage.symbolic
            0.333333333333333*x^3 + 1.57079632679490*x^2 + 3.00000000000000*x
        """
        cdef mpfr_rnd_t rnd = self._base_ring.rnd
        cdef PolynomialRealDense f = self._new(self._degree+1)
        mpfr_set_ui(f._coeffs[0], 0, rnd)
        for i from 0 <= i <= self._degree:
            mpfr_div_ui(f._coeffs[i+1], self._coeffs[i], i+1, rnd)
        return f

    def reverse(self, degree=None):
        """
        Return reverse of the input polynomial thought as a polynomial of
        degree ``degree``.

        If `f` is a degree-`d` polynomial, its reverse is `x^d f(1/x)`.

        INPUT:

        - ``degree`` -- ``None`` or an integer; if specified, truncate or zero
          pad the list of coefficients to this degree before reversing it

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: f = RR['x']([-3, pi, 0, 1])
            sage: f.reverse()
            -3.00000000000000*x^3 + 3.14159265358979*x^2 + 1.00000000000000
            sage: f.reverse(2)
            -3.00000000000000*x^2 + 3.14159265358979*x
            sage: f.reverse(5)
            -3.00000000000000*x^5 + 3.14159265358979*x^4 + x^2

        TESTS:

        We check that this implementation is compatible with the generic one::

            sage: all(f.reverse(d) == Polynomial.reverse(f, d)                          # needs sage.symbolic
            ....:     for d in [None, 0, 1, 2, 3, 4, 5])
            True
        """
        if degree is None: degree = self._degree

        cdef mpfr_rnd_t rnd = self._base_ring.rnd
        cdef PolynomialRealDense f = self._new(degree)

        cdef int i
        for i in range(1+min(degree, self._degree)):
            mpfr_set(f._coeffs[degree-i], self._coeffs[i], rnd)
        for i in range(1 + self._degree, 1 + degree):
            mpfr_set_si(f._coeffs[degree-i], 0, rnd)
        f._normalize()
        return f

    @coerce_binop
    def quo_rem(self, PolynomialRealDense other):
        """
        Return the quotient with remainder of ``self`` by ``other``.

        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [-2, 0, 1])
            sage: g = PolynomialRealDense(RR['x'], [5, 1])
            sage: q, r = f.quo_rem(g)
            sage: q
            x - 5.00000000000000
            sage: r
            23.0000000000000
            sage: q*g + r == f
            True
            sage: fg = f*g
            sage: fg.quo_rem(f)
            (x + 5.00000000000000, 0)
            sage: fg.quo_rem(g)
            (x^2 - 2.00000000000000, 0)

            sage: # needs sage.symbolic
            sage: f = PolynomialRealDense(RR['x'], range(5))
            sage: g = PolynomialRealDense(RR['x'], [pi,3000,4])
            sage: q, r = f.quo_rem(g)
            sage: g*q + r == f
            True

        TESTS:

        Check that :issue:`18467` is fixed::

            sage: S.<x> = RR[]
            sage: z = S.zero()
            sage: z.degree()
            -1
            sage: q, r = z.quo_rem(x)
            sage: q.degree()
            -1
        """
        if other._degree < 0:
            raise ZeroDivisionError("other must be nonzero")
        elif other._degree == 0:
            return self * ~other[0], self._parent.zero()
        elif other._degree > self._degree:
            return self._parent.zero(), self
        cdef mpfr_rnd_t rnd = self._base_ring.rnd
        cdef PolynomialRealDense q, r
        cdef Py_ssize_t i, j
        cdef mpfr_t tmp
        # Make divisor monic for simplicity
        leading = other[other._degree]
        other = other * ~leading
        r = self * ~leading
        q = self._new(self._degree - other._degree)
        # This is the standard division algorithm
        sig_on()
        mpfr_init2(tmp, self._base_ring._prec)
        for i from self._degree >= i >= other._degree:
            mpfr_set(q._coeffs[i-other._degree], r._coeffs[i], rnd)
            for j from 0 <= j < other._degree:
                mpfr_mul(tmp, r._coeffs[i], other._coeffs[j], rnd)
                mpfr_sub(r._coeffs[i-other._degree+j], r._coeffs[i-other._degree+j], tmp, rnd)
            r._degree -= 1
            mpfr_clear(r._coeffs[i])
        mpfr_clear(tmp)
        sig_off()
        r._normalize()
        return q, r * leading

    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [-2, 0, 1])
            sage: f(10)
            98.0000000000000
            sage: f(CC.0)
            -3.00000000000000
            sage: f(2.0000000000000000000000000000000000000000000)
            2.00000000000000
            sage: f(RealField(10)(2))
            2.0
            sage: f(pi)                                                                 # needs sage.symbolic
            1.00000000000000*pi^2 - 2.00000000000000


            sage: f = PolynomialRealDense(RR['x'], range(5))
            sage: f(1)
            10.0000000000000
            sage: f(-1)
            2.00000000000000
            sage: f(0)
            0.000000000000000
            sage: f = PolynomialRealDense(RR['x'])
            sage: f(12)
            0.000000000000000

        TESTS::

            sage: R.<x> = RR[]  # Issue #17311
            sage: (x^2+1)(x=5)
            26.0000000000000
        """
        if len(args) == 1:
            xx = args[0]
        else:
            return Polynomial.__call__(self, *args, **kwds)

        if not isinstance(xx, RealNumber):
            if self._base_ring.has_coerce_map_from(parent(xx)):
                xx = self._base_ring(xx)
            else:
                return Polynomial.__call__(self, xx)

        cdef Py_ssize_t i
        cdef mpfr_rnd_t rnd = self._base_ring.rnd
        cdef RealNumber x = <RealNumber>xx
        cdef RealNumber res

        if (<RealField_class>x._parent)._prec < self._base_ring._prec:
            res = RealNumber(x._parent)
        else:
            res = RealNumber(self._base_ring)
        # Optimize some very useful and common cases:
        if self._degree < 0:
            mpfr_set_ui(res.value, 0, rnd)
        elif mpfr_zero_p(x.value):
            mpfr_set(res.value, self._coeffs[0], rnd)
        elif mpfr_cmp_ui(x.value, 1) == 0:
            mpfr_set(res.value, self._coeffs[0], rnd)
            for i from 0 < i <= self._degree:
                mpfr_add(res.value, res.value, self._coeffs[i], rnd)
        elif mpfr_cmp_si(x.value, -1) == 0:
            mpfr_set(res.value, self._coeffs[0], rnd)
            for i from 2 <= i <= self._degree by 2:
                mpfr_add(res.value, res.value, self._coeffs[i], rnd)
            for i from 1 <= i <= self._degree by 2:
                mpfr_sub(res.value, res.value, self._coeffs[i], rnd)
        else:
            mpfr_set(res.value, self._coeffs[self._degree], rnd)
            for i from self._degree > i >= 0:
                mpfr_mul(res.value, res.value, x.value, rnd)
                mpfr_add(res.value, res.value, self._coeffs[i], rnd)
        return res

    def change_ring(self, R):
        """
        EXAMPLES::

            sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import PolynomialRealDense
            sage: f = PolynomialRealDense(RR['x'], [-2, 0, 1.5])
            sage: f.change_ring(QQ)
            3/2*x^2 - 2
            sage: f.change_ring(RealField(10))
            1.5*x^2 - 2.0
            sage: f.change_ring(RealField(100))
            1.5000000000000000000000000000*x^2 - 2.0000000000000000000000000000
        """
        cdef Py_ssize_t i
        cdef mpfr_rnd_t rnd = self._base_ring.rnd
        cdef PolynomialRealDense f
        if isinstance(R, RealField_class):
            f = PolynomialRealDense(R[self.variable_name()])
            f = f._new(self._degree)
            for i from 0 <= i <= self._degree:
                mpfr_set(f._coeffs[i], self._coeffs[i], rnd)
            return f
        else:
            return Polynomial.change_ring(self, R)


def make_PolynomialRealDense(parent, data):
    """
    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_real_mpfr_dense import make_PolynomialRealDense
        sage: make_PolynomialRealDense(RR['x'], [1,2,3])
        3.00000000000000*x^2 + 2.00000000000000*x + 1.00000000000000
    """
    return PolynomialRealDense(parent, data)
