r"""
Univariate polynomials over `\QQ` implemented via FLINT

AUTHOR:

- Sebastian Pancratz
"""

###############################################################################
#          Copyright (C) 2010 Sebastian Pancratz <sfp@pancratz.org>           #
#                                                                             #
#     Distributed under the terms of the GNU General Public License (GPL)     #
#                                                                             #
#                        http://www.gnu.org/licenses/                         #
###############################################################################

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/gmp.pxi"
include "sage/libs/ntl/decl.pxi"

include "sage/ext/cdefs.pxi"
include "sage/libs/flint/fmpz.pxi"
include "sage/libs/flint/fmpz_poly.pxi"
include "sage/libs/flint/fmpq_poly.pxd"

from sage.interfaces.all import singular as singular_default

from sage.libs.all import pari, pari_gen
from sage.libs.flint.ntl_interface cimport *
from sage.libs.flint.fmpz_poly cimport fmpz_poly_set

from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.rational cimport Rational
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint

from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.structure.element import coerce_binop
from sage.structure.factorization import Factorization


cdef inline bint _do_sig(fmpq_poly_t op):
    """
    Returns 1 when signal handling should be carried out for an operation
    on this polynomial and 0 otherwise.

    Strictly speaking, whether or not signal handling should be carried
    ought to depend on the operation as well as the operands in question.
    For simplicity we carry out signal handling for all but the simplest
    of operands regardless of the operation.

    TESTS::

        sage: R.<t> = QQ[]
        sage: f = 1 + t/2
        sage: g = 2/3 + t^2
        sage: _ = f * g      # indirect doctest
    """
    # Trac #12173: check that the degree is greater than 1000 before computing
    # the max limb size
    return fmpq_poly_length(op) > 0 and \
       (fmpq_poly_degree(op) > 1000 or
        _fmpz_vec_max_limbs(fmpq_poly_numref(op), fmpq_poly_length(op)) > 1)

cdef class Polynomial_rational_flint(Polynomial):
    """
    Univariate polynomials over the rationals, implemented via FLINT.

    Internally, we represent rational polynomial as the quotient of an integer
    polynomial and a positive denominator which is coprime to the content of
    the numerator.
    """

    ###########################################################################
    # Allocation & initialisation                                             #
    ###########################################################################

    cdef Polynomial_rational_flint _new(self):
        """
        Quickly creates a new polynomial object in this class.

        OUTPUT:

        - Polynomial of type Polynomial_rational_flint

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 2/3*t^2
            sage: g = -1/2*t + 2
            sage: f + g           # indirect doctest
            2/3*t^2 - 1/2*t + 2
        """
        cdef Polynomial_rational_flint res = PY_NEW(Polynomial_rational_flint)
        res._parent = self._parent
        res._is_gen = 0
        return res

    cpdef Polynomial _new_constant_poly(self, x, Parent P):
        r"""
        Quickly creates a new constant polynomial with value x in parent P

        ASSUMPTION:

        x must be a rational or convertible to an int.

        EXAMPLE::

        sage: R.<x> = QQ[]
            sage: x._new_constant_poly(2/1,R)
            2
            sage: x._new_constant_poly(2,R)
            2
            sage: x._new_constant_poly("2",R)
            2
            sage: x._new_constant_poly("2.1",R)
            Traceback (most recent call last):
            ...
            ValueError: invalid literal for int() with base 10: '2.1'

        """
        cdef Polynomial_rational_flint res = PY_NEW(Polynomial_rational_flint)
        res._parent = P
        res._is_gen = <char>0
        if PY_TYPE_CHECK(x, int):
            fmpq_poly_set_si(res.__poly, <int> x)

        elif PY_TYPE_CHECK(x, Integer):
            fmpq_poly_set_mpz(res.__poly, (<Integer> x).value)

        elif PY_TYPE_CHECK(x, Rational):
            fmpq_poly_set_mpq(res.__poly, (<Rational> x).value)

        else:
            fmpq_poly_set_si(res.__poly, int(x))
        return res


    def __cinit__(self):
        """
        Initialises the underlying data structure.

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 2/3 * t - 7  #indirect doctest
        """
        fmpq_poly_init(self.__poly)

    def __dealloc__(self):
        """
        Deallocates the underlying data structure.

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 1/3 * t
            sage: del f        # untested
        """
        fmpq_poly_clear(self.__poly)

    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        Initialises the associated data for the polynomial self.

        INPUT:

        - ``parent`` - Polynomial ring, the parent of self
        - ``x`` - Data for the new polynomial self, e.g. a polynomial, an
          integer, a rational, a list of rationals, a dictionary with keys
          the degrees and the rational coefficients, etc (default: ``None``)
        - `check`` - Whether the integrity of the data needs to be verified,
          largely ignored by this method (default: ``True``)
        - ``is_gen`` - Whether self shall be initialised as the generator of
          the parent polynomial ring
        - ``construct`` - Whether the element shall always be constructed
          as an independent copy of any input data (default: ``False``)

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = -4 * t^2 + 1/3 * t - 1/7  # indirect doctest

            sage: f = ZZ['x']([1..10^6])
            sage: g = f.change_ring(QQ)
            sage: g[:10]
            10*x^9 + 9*x^8 + 8*x^7 + 7*x^6 + 6*x^5 + 5*x^4 + 4*x^3 + 3*x^2 + 2*x + 1
        """
        cdef long deg
        cdef unsigned long n
        cdef Rational c
        cdef list L1
        cdef mpq_t * L2

        Polynomial.__init__(self, parent, is_gen=is_gen)

        if is_gen:
            fmpq_poly_set_coeff_si(self.__poly, 1, 1)

        elif PY_TYPE_CHECK(x, Polynomial_rational_flint):
            fmpq_poly_set(self.__poly, (<Polynomial_rational_flint> x).__poly)

        elif PY_TYPE_CHECK(x, int):
            fmpq_poly_set_si(self.__poly, <int> x)

        elif PY_TYPE_CHECK(x, Integer):
            fmpq_poly_set_mpz(self.__poly, (<Integer> x).value)

        elif PY_TYPE_CHECK(x, Rational):
            fmpq_poly_set_mpq(self.__poly, (<Rational> x).value)

        elif PY_TYPE_CHECK(x, list) or PY_TYPE_CHECK(x, tuple):

            if len(x) == 0:
                return
            elif len(x) == 1:
                Polynomial_rational_flint.__init__(self, parent, x[0], \
                                check=check, is_gen=False, construct=construct)
                return

            L1 = [e if isinstance(e, Rational) else Rational(e) for e in x]
            n  = <unsigned long> len(x)
            sig_on()
            L2 = <mpq_t *> sage_malloc(n * sizeof(mpq_t))
            for deg from 0 <= deg < n:
                mpq_init(L2[deg])
                mpq_set(L2[deg], (<Rational> L1[deg]).value)
            fmpq_poly_set_array_mpq(self.__poly, L2, n)
            for deg from 0 <= deg < n:
                mpq_clear(L2[deg])
            sage_free(L2)
            sig_off()

#           deg = 0
#           for e in x:
#               c = Rational(e)
#               fmpq_poly_set_coeff_mpq(self.__poly, deg, c.value)
#               deg += 1

        elif PY_TYPE_CHECK(x, dict):
            for deg, e in x.iteritems():
                c = Rational(e)
                fmpq_poly_set_coeff_mpq(self.__poly, deg, c.value)

        elif PY_TYPE_CHECK(x, pari_gen):
            k = self._parent.base_ring()
            x = [k(w) for w in x.list()]
            Polynomial_rational_flint.__init__(self, parent, x, check=True, \
                                             is_gen=False, construct=construct)

        elif PY_TYPE_CHECK(x, Polynomial_integer_dense_flint):
            fmpq_poly_set_fmpz_poly(self.__poly, (<Polynomial_integer_dense_flint>x).__poly)

        elif PY_TYPE_CHECK(x, Polynomial):
            k = self._parent.base_ring()
            x = [k(w) for w in list(x)]
            Polynomial_rational_flint.__init__(self, parent, x, check=True, \
                                             is_gen=False, construct=construct)

        elif PY_TYPE_CHECK(x, FractionFieldElement) and (x.parent().base() is parent or x.parent().base() == parent) and x.denominator() == 1:
            x = x.numerator()
            Polynomial_rational_flint.__init__(self, parent, x, check=check, \
                                            is_gen=is_gen, construct=construct)

        else:
            x = parent.base_ring()(x)
            Polynomial_rational_flint.__init__(self, parent, x, check=check, \
                                            is_gen=is_gen, construct=construct)

    def __reduce__(self):
        """
        This is used when pickling polynomials.

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 2/3 * t^2 + 1
            sage: r = f.__reduce__(); r
            (<type 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint'>, (Univariate Polynomial Ring in t over Rational Field, [1, 0, 2/3], False, False))
            sage: r[0](*r[1])
            2/3*t^2 + 1
            sage: loads(dumps(f)) == f
            True
        """
        return Polynomial_rational_flint, \
               (self.parent(), self.list(), False, self.is_gen())

    def __copy__(self):
        """
        Returns a copy of self.

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 4/5 * t^3 - 1/17
            sage: copy(f) == f
            True
        """
        cdef Polynomial_rational_flint res = self._new()
        fmpq_poly_set(res.__poly, self.__poly)
        return res

    def _singular_(self, singular=singular_default, have_ring=False):
        """
        Returns a Singular representation of self.

        INPUT:

        - ``singular`` - Singular interpreter (default: default interpreter)
        - ``have_ring`` - set to True if the ring was already set in Singular

        EXAMPLES::

            sage: P.<x> = PolynomialRing(QQ)
            sage: f = 3*x^2 + 2*x + 5
            sage: singular(f)
            3*x^2+2*x+5
        """
        if not have_ring:
            self._parent._singular_(singular).set_ring()  # Expensive!
        return singular(self._singular_init_())

    def list(self):
        """
        Returns a list with the coefficients of self.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2/2 + t^3/3 + t^4/4
            sage: f.list()
            [1, 1, 1/2, 1/3, 1/4]
            sage: g = R(0)
            sage: g.list()
            []
        """
        cdef unsigned long length = fmpq_poly_length(self.__poly)
        return [self[n] for n in range(length)]

    ###########################################################################
    # Basis access                                                            #
    ###########################################################################

    def degree(self):
        """
        Returns the degree of self.

        By convention, the degree of the zero polynomial is -1.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2/2 + t^3/3 + t^4/4
            sage: f.degree()
            4
            sage: g = R(0)
            sage: g.degree()
            -1
        """
        cdef Integer deg = PY_NEW(Integer)
        mpz_set_si(deg.value, fmpq_poly_degree(self.__poly))
        return deg

    def __getitem__(self, n):
        """
        Returns coefficient of the monomial of degree `n` if `n` is an integer,
        returns the monomials of self of degree in slice `n` if `n` is a slice.

        INPUT:

        - ``n`` - Degree of the monomial whose coefficient is to be returned
                  or a slice.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2/2 + t^3/3 + t^4/4
            sage: f[-1], f[0], f[3], f[5]            # indirect doctest
            (0, 1, 1/3, 0)
            sage: f[1:3]                             # indirect doctest
            1/2*t^2 + t
        """
        cdef Rational z = PY_NEW(Rational)
        cdef Polynomial_rational_flint res = self._new()
        cdef bint do_sig = _do_sig(self.__poly)
        if isinstance(n, slice):
            start, stop, step = n.indices(self.degree() + 1)
            if do_sig: sig_on()
            fmpq_poly_get_slice(res.__poly, self.__poly, start, stop)
            if do_sig: sig_off()
            return res
        else:
            if 0 <= n and n < fmpq_poly_length(self.__poly):
                fmpq_poly_get_coeff_mpq(z.value, self.__poly, n)
            return z

    cpdef _unsafe_mutate(self, unsigned long n, value):
        """
        Sets the `n`th coefficient of self to value.

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2/2 + t^3/3 + t^4/4
            sage: f._unsafe_mutate(4, 1/5)
            sage: f
            1/5*t^4 + 1/3*t^3 + 1/2*t^2 + t + 1

        WARNING:

        Polynomials in Sage are meant to be immutable, and some methods may
        rely on this convention.  This method should be used only with the
        utmost care.
        """
        cdef bint do_sig = _do_sig(self.__poly)

        if PY_TYPE_CHECK(value, int):
            if do_sig: sig_on()
            fmpq_poly_set_coeff_si(self.__poly, n, value)
            if do_sig: sig_off()
        elif PY_TYPE_CHECK(value, Integer):
            if do_sig: sig_on()
            fmpq_poly_set_coeff_mpz(self.__poly, n, (<Integer> value).value)
            if do_sig: sig_off()
        elif PY_TYPE_CHECK(value, Rational):
            if do_sig: sig_on()
            fmpq_poly_set_coeff_mpq(self.__poly, n, (<Rational> value).value)
            if do_sig: sig_off()
        else:
            value = Rational(value)
            if do_sig: sig_on()
            fmpq_poly_set_coeff_mpq(self.__poly, n, (<Rational> value).value)
            if do_sig: sig_off()

    def __call__(self, *x, **kwds):
        """
        Calls this polynomial with the given parameters, which can be
        interpreted as polynomial composition or evaluation by this
        method.

        If the argument is not simply an integer, a rational, or a
        polynomial, the call is passed on to the generic implementation
        in the Polynomial class.

        EXAMPLES:

        The first example illustrates polynomial composition::

            sage: R.<t> = QQ[]
            sage: f = t^2 - 1
            sage: g = t + 1
            sage: f(g)          # indirect doctest
            t^2 + 2*t

        Now we illustrate how a polynomial can be evaluated at a rational
        number::

            sage: f(-2/3)       # indirect doctest
            -5/9
        """
        cdef Polynomial_rational_flint f
        cdef Rational r

        if len(x) == 1:
            a = x[0]
            if isinstance(a, Polynomial_rational_flint):
                f = (<Polynomial_rational_flint> a)._new()
                sig_on()
                fmpq_poly_compose(f.__poly, self.__poly, \
                    (<Polynomial_rational_flint> a).__poly)
                sig_off()
                return f
            if isinstance(a, Rational):
                r = PY_NEW(Rational)
                sig_on()
                fmpq_poly_evaluate_mpq(r.value, self.__poly, (<Rational> a).value)
                sig_off()
                return r
            if isinstance(a, Integer):
                r = PY_NEW(Rational)
                sig_on()
                fmpq_poly_evaluate_mpz(r.value, self.__poly, (<Integer> a).value)
                sig_off()
                return r

        return Polynomial.__call__(self, *x, **kwds)

    cpdef Polynomial truncate(self, long n):
        """
        Returns self truncated modulo `t^n`.

        INPUT:

        - ``n`` - The power of `t` modulo which self is truncated

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 - t + 1/2*t^2 - 1/3*t^3
            sage: f.truncate(0)
            0
            sage: f.truncate(2)
            -t + 1
        """
        cdef Polynomial_rational_flint res
        cdef bint do_sig

        if (n >= fmpq_poly_length(self.__poly)):
            return self
        else:
            res = self._new()
            if n > 0:
                do_sig = _do_sig(self.__poly)
                if do_sig: sig_on()
                fmpq_poly_get_slice(res.__poly, self.__poly, 0, n)
                if do_sig: sig_off()
            return res

    def reverse(self, n = None):
        """
        Reverses the coefficients of self - thought of as a polynomial of
        length `n`.

        Assumes that whenever `n` is not ``None`` it is an integral value
        that fits into an unsigned long.  Otherwise, a ValueError is raised.

        INPUT:

        - ``n`` - Integral value that fits in an unsigned long:  the power
          of `t` modulo which we consider self before reversing it
          (default:  ``None``, interpreted as the length of self)

        OUTPUT:

        - The reversed polynomial as a Polynomial_rational_flint

        EXAMPLES:

        We first consider the simplest case, where we reverse all coefficients
        of a polynomial and obtain a polynomial of the same degree::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2 / 2 + t^3 / 3 + t^4 / 4
            sage: f.reverse()
            t^4 + t^3 + 1/2*t^2 + 1/3*t + 1/4

        Next, an example we the returned polynomial has lower degree because
        the original polynomial has low coefficients equal to zero::

            sage: R.<t> = QQ[]
            sage: f = 3/4*t^2 + 6*t^7
            sage: f.reverse()
            3/4*t^5 + 6

        The next example illustrates the passing of a value for `n` less than
        the length of self, notationally resulting in truncation prior to
        reversing::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2 / 2 + t^3 / 3 + t^4 / 4
            sage: f.reverse(3)
            t^2 + t + 1/2

        Now we illustrate the passing of a value for `n` greater than the
        length of self, notationally resulting in zero padding at the top
        end prior to reversing::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2 / 2 + t^3 / 3
            sage: f.reverse(5)
            t^4 + t^3 + 1/2*t^2 + 1/3*t

        TESTS::

        We illustrate two ways in which the interpretation of `n` as an
        unsigned long int may fail.  Firstly, an integral value which is
        too large, yielding an OverflowError::

            sage: R.<t> = QQ[]
            sage: f = 1 + t/2
            sage: f.reverse(2**64)
            Traceback (most recent call last):
            ...
            OverflowError: long int too large to convert

        Secondly, a value which cannot be converted to an integral value,
        resulting in a ValueError::

            sage: R.<t> = QQ[]
            sage: f = 1 + t/2
            sage: f.reverse(I)
            Traceback (most recent call last):
            ...
            ValueError: cannot convert I to int
        """
        cdef unsigned long len
        cdef Polynomial_rational_flint res
        cdef bint do_sig

        if n is None:
            len = fmpq_poly_length(self.__poly)
        else:
            len = <unsigned long> n

        res = self._new()
        do_sig = _do_sig(self.__poly)
        if do_sig: sig_on()
        fmpq_poly_reverse(res.__poly, self.__poly, len)
        if do_sig: sig_off()
        return res

    ###########################################################################
    # Comparisons                                                             #
    ###########################################################################

    def is_zero(self):
        """
        Returns whether or not self is the zero polynomial.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 - t + 1/2*t^2 - 1/3*t^3
            sage: f.is_zero()
            False
            sage: R(0).is_zero()
            True
        """
        return bool(fmpq_poly_is_zero(self.__poly))

    def __nonzero__(self):
        """
        Returns whether or not self is non-zero.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 - t + 1/2*t^2 - 1/3*t^3
            sage: bool(f)
            True
            sage: bool(R(0))
            False
        """
        return not fmpq_poly_is_zero(self.__poly)

    ###########################################################################
    # Shifting                                                                #
    ###########################################################################

    def __lshift__(self, n):
        """
        Notationally multiplies self by `t^n`.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: t << 10                     # indirect doctest
            t^11

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = R.random_element(1000)
            sage: (f << 23) >> 23 == f        # indirect doctest
            True
        """
        cdef unsigned long k = <unsigned long> n
        cdef Polynomial_rational_flint f = <Polynomial_rational_flint> self
        cdef Polynomial_rational_flint res
        cdef bint do_sig

        if k == 0 or fmpq_poly_is_zero(f.__poly):
            return self
        else:
            res = f._new()
            do_sig = fmpq_poly_length(f.__poly) > 5000 or n > 5000

            if do_sig: sig_on()
            fmpq_poly_shift_left(res.__poly, f.__poly, k)
            if do_sig: sig_off()
            return res

    def __rshift__(self, n):
        """
        Notationally returns the quotient of Euclidean division of self
        by `t^n`.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2/2 + t^3/3 + t^4/4
            sage: f >> 2
            1/4*t^2 + 1/3*t + 1/2
        """
        cdef unsigned long k = <unsigned long> n
        cdef Polynomial_rational_flint f = <Polynomial_rational_flint> self
        cdef Polynomial_rational_flint res
        cdef bint do_sig

        if k == 0 or fmpq_poly_is_zero(f.__poly):
            return self
        else:
            res = f._new()
            do_sig = _do_sig(f.__poly)

            if do_sig: sig_on()
            fmpq_poly_shift_right(res.__poly, f.__poly, k)
            if do_sig: sig_off()
            return res

    ###########################################################################
    # Arithmetic                                                              #
    ###########################################################################

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Returns the sum of two rational polynomials.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 2/3 + t + 2*t^3
            sage: g = -1 + t/3 - 10/11*t^4
            sage: f + g
            -10/11*t^4 + 2*t^3 + 4/3*t - 1/3

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = R.random_element(2000)
            sage: f + f == 2 * f              # indirect doctest
            True
        """
        cdef Polynomial_rational_flint op2 = <Polynomial_rational_flint> right
        cdef Polynomial_rational_flint res = self._new()
        cdef bint do_sig = _do_sig(self.__poly) or _do_sig(op2.__poly)

        if do_sig: sig_on()
        fmpq_poly_add(res.__poly, self.__poly, op2.__poly)
        if do_sig: sig_off()
        return res

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Returns the difference of two rational polynomials.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = -10/11*t^4 + 2*t^3 + 4/3*t - 1/3
            sage: g = 2*t^3
            sage: f - g                                 # indirect doctest
            -10/11*t^4 + 4/3*t - 1/3

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = R.random_element(2000)
            sage: f - f/2 == 1/2 * f          # indirect doctest
            True
            sage: f[:1000] == f - f[1000:]    # indirect doctest
            True
        """
        cdef Polynomial_rational_flint op2 = <Polynomial_rational_flint> right
        cdef Polynomial_rational_flint res = self._new()
        cdef bint do_sig = _do_sig(self.__poly) or _do_sig(op2.__poly)

        if do_sig: sig_on()
        fmpq_poly_sub(res.__poly, self.__poly, op2.__poly)
        if do_sig: sig_off()
        return res

    cpdef ModuleElement _neg_(self):
        """
        Returns the difference of two rational polynomials.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 3*t/2
            sage: -f            # indirect doctest
            -3/2*t

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = R.random_element(2000)
            sage: f + (-f) == 0               # indirect doctest
            True
        """
        cdef Polynomial_rational_flint res = self._new()
        cdef bint do_sig = _do_sig(self.__poly)

        if do_sig: sig_on()
        fmpq_poly_neg(res.__poly, self.__poly)
        if do_sig: sig_off()
        return res

    @coerce_binop
    def quo_rem(self, right):
        """
        Returns the quotient and remainder of the Euclidean division of
        self and right.

        Raises a ZerodivisionError if right is zero.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = R.random_element(2000)
            sage: g = R.random_element(1000)
            sage: q, r = f.quo_rem(g)
            sage: f == q*g + r
            True
        """
        if right.is_zero():
            raise ZeroDivisionError, "division by zero polynomial"
        if self.is_zero():
            return self, self

        cdef Polynomial_rational_flint qq = self._new()
        cdef Polynomial_rational_flint rr = self._new()

        sig_on()
        fmpq_poly_divrem(qq.__poly, rr.__poly, self.__poly, \
                         (<Polynomial_rational_flint> right).__poly)
        sig_off()
        return qq, rr

    @coerce_binop
    def gcd(self, right):
        """
        Returns the (monic) greatest common divisor of self and right.

        Corner cases:  if self and right are both zero, returns zero.  If
        only one of them is zero, returns the other polynomial, up to
        normalisation.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = -2 + 3*t/2 + 4*t^2/7 - t^3
            sage: g = 1/2 + 4*t + 2*t^4/3
            sage: f.gcd(g)
            1
            sage: f = (-3*t + 1/2) * f
            sage: g = (-3*t + 1/2) * (4*t^2/3 - 1) * g
            sage: f.gcd(g)
            t - 1/6
        """
        cdef Polynomial_rational_flint res = self._new()

        sig_on()
        fmpq_poly_gcd(res.__poly, self.__poly, \
                (<Polynomial_rational_flint> right).__poly)
        sig_off()
        return res

    @coerce_binop
    def lcm(self, right):
        """
        Returns the monic (or zero) least common multiple of self and right.

        Corner cases:  if either of self and right are zero, returns zero.
        This behaviour is ensures that the relation lcm(a,b) gcd(a,b) == a b
        holds up to multiplication by rationals.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = -2 + 3*t/2 + 4*t^2/7 - t^3
            sage: g = 1/2 + 4*t + 2*t^4/3
            sage: f.lcm(g)
            t^7 - 4/7*t^6 - 3/2*t^5 + 8*t^4 - 75/28*t^3 - 66/7*t^2 + 87/8*t + 3/2
            sage: f.lcm(g) * f.gcd(g) // (f * g)
            -3/2
        """
        cdef Polynomial_rational_flint res = self._new()

        sig_on()
        fmpq_poly_lcm(res.__poly, self.__poly, \
                      (<Polynomial_rational_flint> right).__poly)
        sig_off()
        return res

    @coerce_binop
    def xgcd(self, right):
        """
        Returns polynomials d, s, and t such that d == s * self + t * right,
        where d is the (monic) greatest common divisor of self and right.
        The choice of s and t is not specified any further.

        Corner cases:  if self and right are zero, returns zero polynomials.
        Otherwise, if only self is zero, returns (d, s, t) = (right, 0, 1) up
        to normalisation, and similarly if only right is zero.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 2/3 + 3/4 * t - t^2
            sage: g = -3 + 1/7 * t
            sage: f.xgcd(g)
            (1, -12/5095, -84/5095*t - 1701/5095)

        TESTS:

        The following example used to crash (cf. #11771)::

            sage: R.<t> = QQ[]
            sage: f = 10**383 * (t+1)
            sage: g = 10**445 * t^2 + 1
            sage: r = f.xgcd(g)
            sage: r[0] == f.gcd(g)
            True
            sage: r[1]*f + r[2]*g == r[0]
            True
        """
        cdef Polynomial_rational_flint d = self._new()
        cdef Polynomial_rational_flint s = self._new()
        cdef Polynomial_rational_flint t = self._new()

        sig_on()
        fmpq_poly_xgcd(d.__poly, s.__poly, t.__poly, self.__poly, (<Polynomial_rational_flint>right).__poly)
        sig_off()
        return d, s, t

    cpdef RingElement _mul_(self, RingElement right):
        """
        Returns the product of self and right.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = -1 + 3*t/2 - t^3
            sage: g = 2/3 + 7/3*t + 3*t^2
            sage: f * g                           # indirect doctest
            -3*t^5 - 7/3*t^4 + 23/6*t^3 + 1/2*t^2 - 4/3*t - 2/3

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = R.random_element(2000)
            sage: g = R.random_element(2000)
            sage: (f + g) * (f - g) == f^2 - g^2  # indirect doctest
            True
        """
        cdef Polynomial_rational_flint op2 = <Polynomial_rational_flint> right
        cdef Polynomial_rational_flint res = self._new()
        cdef bint do_sig = _do_sig(self.__poly) or _do_sig(op2.__poly)

        if do_sig: sig_on()
        fmpq_poly_mul(res.__poly, self.__poly, op2.__poly)
        if do_sig: sig_off()
        return res

    cpdef ModuleElement _rmul_(self, RingElement left):
        r"""
        Returns left * self, where left is a rational number.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 3/2*t^3 - t + 1/3
            sage: 6 * f                  # indirect doctest
            9*t^3 - 6*t + 2
        """
        cdef Polynomial_rational_flint res = self._new()
        cdef bint do_sig = _do_sig(self.__poly)

        if do_sig: sig_on()
        fmpq_poly_scalar_mul_mpq(res.__poly, self.__poly, \
                                 (<Rational> left).value)
        if do_sig: sig_off()
        return res

    cpdef ModuleElement _lmul_(self, RingElement right):
        r"""
        Returns self * right, where right is a rational number.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 3/2*t^3 - t + 1/3
            sage: f * 6                   # indirect doctest
            9*t^3 - 6*t + 2
        """
        cdef Polynomial_rational_flint res = self._new()
        cdef bint do_sig = _do_sig(self.__poly)

        if do_sig: sig_on()
        fmpq_poly_scalar_mul_mpq(res.__poly, self.__poly, \
                                 (<Rational> right).value)
        if do_sig: sig_off()
        return res

    def __pow__(Polynomial_rational_flint self, exp, ignored):
        """
        Returns self raised to the power of exp.

        The corner case of exp == 0 is handled by returning the constant
        polynomial 1.  Note that this includes the case 0^0 == 1.

        This method only supports integral values for exp that fit into
        a signed long int.

        INPUT:

        - exp - Exponent

        OUTPUT:

        - Polynomial; self raised to the power of exp

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1/2 + 2*t - t^2/3
            sage: f^0
            1
            sage: f^3
            -1/27*t^6 + 2/3*t^5 - 23/6*t^4 + 6*t^3 + 23/4*t^2 + 3/2*t + 1/8

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 1/2 + t
            sage: f^0
            1
            sage: R(0)^0
            1

        We verify the size condition on the exponent::

            sage: R.<t> = QQ[]
            sage: (1 + t)^(2^63)
            Traceback (most recent call last):
            ...
            OverflowError: Python int too large to convert to C long
        """
        cdef long n
        cdef Polynomial_rational_flint res

        if not PY_TYPE_CHECK(exp, Integer):
            try:
                exp = Integer(exp)
            except TypeError:
                raise TypeError, "non-integral exponents not supported"

        n = <long> exp

        if n < 0:
            if fmpq_poly_is_zero(self.__poly):
                raise ZeroDivisionError, "negative exponent in power of zero"
            res = self._new()
            sig_on()
            fmpq_poly_pow(res.__poly, self.__poly, -n)
            sig_off()
            return ~res
        else:
            res = self._new()
            if self._is_gen:
                fmpq_poly_set_coeff_si(res.__poly, n, 1)
            else:
                sig_on()
                fmpq_poly_pow(res.__poly, self.__poly, n)
                sig_off()
            return res

    def __floordiv__(Polynomial_rational_flint self, right):
        """
        Returns the quotient of self and right obtain by Euclidean division.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = t^3 - t/2 + 1/5
            sage: g = 2/3*t - 1
            sage: f // g                       # indirect doctest
            3/2*t^2 + 9/4*t + 21/8

        TESTS::

            sage: R.<t> = QQ[]
            sage: f  = R.random_element(1000)
            sage: g  = R.random_element(500)
            sage: if g == 0: g = R(1)
            sage: qr = f.quo_rem(g)
            sage: q  = f // g                  # indirect doctest
            sage: qr[0] == q
            True
        """
        cdef Polynomial_rational_flint res
        cdef bint do_sig

        if right == 0:
            raise ZeroDivisionError, "division by zero polynomial"

        if not PY_TYPE_CHECK(right, Polynomial_rational_flint):
            if right in QQ:
                res = self._new()
                do_sig = _do_sig(self.__poly)

                if do_sig: sig_on()
                fmpq_poly_scalar_div_mpq(res.__poly, self.__poly,
                                                  (<Rational> QQ(right)).value)
                if do_sig: sig_off()
                return res

            right = self._parent(right)

        res = self._new()
        sig_on()
        fmpq_poly_div(res.__poly, self.__poly,
                                     (<Polynomial_rational_flint>right).__poly)
        sig_off()
        return res

    def __mod__(Polynomial_rational_flint self, right):
        """
        Returns the remainder of self and right obtain by Euclidean division.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = t^3 - t/2 + 1/5
            sage: g = 2/3*t - 1
            sage: f % g                        # indirect doctest
            113/40

        TESTS::

            sage: R.<t> = QQ[]
            sage: f  = R.random_element(1000)
            sage: g  = R.random_element(500)
            sage: if g == 0: g = R(1)
            sage: qr = f.quo_rem(g)
            sage: r  = f % g                   # indirect doctest
            sage: qr[1] == r
            True
        """
        cdef Polynomial_rational_flint res

        if right == 0:
            raise ZeroDivisionError, "division by zero polynomial"

        if not PY_TYPE_CHECK(right, Polynomial_rational_flint):
            right = self._parent(right)

        res = self._new()
        sig_on()
        fmpq_poly_rem(res.__poly, self.__poly,
                                     (<Polynomial_rational_flint>right).__poly)
        sig_off()
        return res

    ###########################################################################
    # Further methods                                                         #
    ###########################################################################

    def numerator(self):
        """
        Returns the numerator of self.

        Representing self as the quotient of an integer polynomial and
        a positive integer denominator (coprime to the content of the
        polynomial), returns the integer polynomial.

        EXAMPLE::

            sage: R.<t> = QQ[]
            sage: f = (3 * t^3 + 1) / -3
            sage: f.numerator()
            -3*t^3 - 1
        """
        cdef Polynomial_integer_dense_flint num = \
                                         PY_NEW(Polynomial_integer_dense_flint)
        parent = ZZ[self.variable_name()]
        sig_on()
        Polynomial_integer_dense_flint.__init__(num, parent, x=None, \
                                    check=False, is_gen=False, construct=False)
        fmpq_poly_get_numerator(num.__poly, self.__poly)
        sig_off()
        return num

    def denominator(self):
        """
        Returns the denominator of self.

        EXAMPLE::

            sage: R.<t> = QQ[]
            sage: f = (3 * t^3 + 1) / -3
            sage: f.denominator()
            3
        """
        cdef Integer den = PY_NEW(Integer)
        if fmpq_poly_denref(self.__poly) is NULL:
            mpz_set_ui(den.value, 1)
        else:
            fmpz_get_mpz(den.value, <fmpz *> fmpq_poly_denref(self.__poly))
        return den

    def _derivative(self, var = None):
        """
        Returns the derivative of self with respect to ``var``.

        INPUT:

        -  ``var`` - Must be either (equal to) the generator of the polynomial
           ring to which this polynomial belongs, or ``None``; either way the
           behaviour is the same.

        OUTPUT:

        -  Derivative as a ``Polynomial_rational_flint``

        .. seealso:: :meth:`.derivative`

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^4 - x - 1
            sage: f._derivative()
            4*x^3 - 1
            sage: f._derivative(None)
            4*x^3 - 1
            sage: f._derivative(2*x)
            Traceback (most recent call last):
            ...
            ValueError: Cannot differentiate with respect to 2*x
            sage: y = var("y")
            sage: f._derivative(y)
            Traceback (most recent call last):
            ...
            ValueError: Cannot differentiate with respect to y
        """
        cdef Polynomial_rational_flint der
        cdef bint do_sig

        if var is not None and var != self._parent.gen():
            raise ValueError, "Cannot differentiate with respect to %s" %var

        der = self._new()
        do_sig = _do_sig(self.__poly)

        if do_sig: sig_on()
        fmpq_poly_derivative(der.__poly, self.__poly)
        if do_sig: sig_off()
        return der

    def real_root_intervals(self):
        """
        Returns isolating intervals for the real roots of self.

        EXAMPLES:

        We compute the roots of the characteristic polynomial of some
        Salem numbers::

            sage: R.<t> = QQ[]
            sage: f = 1 - t^2 - t^3 - t^4 + t^6
            sage: f.real_root_intervals()
            [((1/2, 3/4), 1), ((1, 3/2), 1)]
        """
        from sage.rings.polynomial.real_roots import real_roots
        return real_roots(self)

    @coerce_binop
    def resultant(Polynomial_rational_flint self, right):
        r"""
        Returns the resultant of self and right.

        Enumerating the roots over `\QQ` as `r_1, \cdots, r_m` and
        `s_1, \cdots, s_n` and letting `x` and `y` denote the leading
        coefficients of `f` and `g`, the resultant of the two polynomials
        is defined by

        .. math::

            x^{\deg g} y^{\deg f} \prod_{i,j} (r_i - s_j).

        Corner cases:  if one of the polynomials is zero, the resultant
        is zero.  Note that otherwise if one of the polynomials is constant,
        the last term in the above is the empty product.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = (t - 2/3) * (t + 4/5) * (t - 1)
            sage: g = (t - 1/3) * (t + 1/2) * (t + 1)
            sage: f.resultant(g)
            119/1350
            sage: h = (t - 1/3) * (t + 1/2) * (t - 1)
            sage: f.resultant(h)
            0
        """
        cdef Rational res = PY_NEW(Rational)
        cdef fmpq_t t
        fmpq_init(t)
        sig_on()
        fmpq_poly_resultant(t, self.__poly, \
                            (<Polynomial_rational_flint>right).__poly)
        fmpq_get_mpq(res.value, t)
        sig_off()
        fmpq_clear(t)
        return res

    def is_irreducible(self):
        r"""
        Returns whether self is irreducible.

        This method computes the primitive part as an element of `\ZZ[t]` and
        calls the method ``is_irreducible`` for elements of that polynomial
        ring.

        By definition, over any integral domain, an element `r` is irreducible
        if and only if it is non-zero, not a unit and whenever `r = ab` then
        `a` or `b` is a unit.

        OUTPUT:

        -  ``bool`` - Whether this polynomial is irreducible

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: (t^2 + 2).is_irreducible()
            True
            sage: (t^2 - 1).is_irreducible()
            False

        TESTS::

            sage: R.<t> = QQ[]
            sage: R(0).is_irreducible()
            False
            sage: R(-1/2).is_irreducible()
            False
            sage: (t+1).is_irreducible()
            True
        """
        cdef Polynomial_integer_dense_flint primitive
        cdef unsigned long length = fmpq_poly_length(self.__poly)

        if length < 2:
            return False
        elif length == 2:
            return True
        else:
            primitive = PY_NEW(Polynomial_integer_dense_flint)
            parent = ZZ[self.variable_name()]
            sig_on()
            Polynomial_integer_dense_flint.__init__(primitive, parent, \
                             x=None, check=True, is_gen=False, construct=False)

            fmpq_poly_get_numerator(primitive.__poly, self.__poly)
            fmpz_poly_primitive_part(primitive.__poly, primitive.__poly)

            sig_off()
            return primitive.is_irreducible()

    ###########################################################################
    # Methods using PARI                                                      #
    ###########################################################################

    def galois_group(self, pari_group = False, algorithm = 'pari'):
        """
        Returns the Galois group of self as a permutation group.

        INPUT:

        -  ``self`` - Irreducible polynomial

        -  ``pari_group`` - bool (default: ``False``); if ``True`` instead
           return the Galois group as a PARI group.  This has a useful label
           in it, and may be slightly faster since it doesn't require looking
           up a group in Gap.  To get a permutation group from a PARI
           group ``P``, type ``PermutationGroup(P)``.

        -  ``algorithm`` - ``'pari'``, ``'kash'``, ``'magma'`` (default:
           ``'pari'``, except when the degree is at least 12 in which case
           ``'kash'`` is tried).

        OUTPUT:

        -  Galois group

        ALGORITHM:

        The Galois group is computed using PARI in C library mode, or possibly
        KASH or MAGMA.

        .. note::

            The PARI documentation contains the following warning: The method
            used is that of resolvent polynomials and is sensitive to the
            current precision. The precision is updated internally but, in very
            rare cases, a wrong result may be returned if the initial precision
            was not sufficient.

            MAGMA does not return a provably correct result.  Please see the
            MAGMA documentation for how to obtain a provably correct result.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: G = f.galois_group(); G            # optional - database_gap
            Transitive group number 5 of degree 4
            sage: G.gens()                           # optional - database_gap
            [(1,2,3,4), (1,2)]
            sage: G.order()                          # optional - database_gap
            24

        It is potentially useful to instead obtain the corresponding PARI
        group, which is little more than a 4-tuple.  See the PARI manual for
        the exact details.  (Note that the third entry in the tuple is in the
        new standard ordering.)

        ::

            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: G = f.galois_group(pari_group=True); G
            PARI group [24, -1, 5, "S4"] of degree 4
            sage: PermutationGroup(G)                # optional - database_gap
            Transitive group number 5 of degree 4

        You can use KASH to compute Galois groups as well.  The advantage is
        that KASH can compute Galois groups of fields up to degree 23, whereas
        PARI only goes to degree 11.  (In my not-so-thorough experiments PARI
        is faster than KASH.)

        ::

            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: f.galois_group(algorithm='kash')   # optional - kash
            Transitive group number 5 of degree 4

            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: f.galois_group(algorithm='magma')  # optional - magma, database_gap
            Transitive group number 5 of degree 4

        TESTS:

        We illustrate the behaviour in the case of reducible polynomials::

            sage: R.<t> = QQ[]
            sage: f = (1 + t)^2
            sage: f.galois_group()
            Traceback (most recent call last):
            ...
            ValueError: The polynomial must be irreducible
        """
        from sage.groups.all import PariGroup, PermutationGroup, TransitiveGroup

        if not self.is_irreducible():
            raise ValueError, "The polynomial must be irreducible"

        if self.degree() > 11 and algorithm == 'pari':
            algorithm = 'kash'

        if self.degree() > 23 and algorithm == 'kash':
            raise NotImplementedError, "Galois group computation is " + \
                "supported for degrees up to 11 using Pari, or up to 23 " + \
                "if the optional package KASH is installed.  Try " + \
                "algorithm='magma' if you have magma."

        if algorithm == 'pari':
            G = self._pari_().Polrev().polgalois()
            H = PariGroup(G, self.degree())
            if pari_group:
                return H
            else:
                return PermutationGroup(H)

        elif algorithm == 'kash':
            try:
                from sage.interfaces.all import kash
                kash.eval('X := PolynomialRing(RationalField()).1')
                s = self._repr(name='X')
                G = kash('Galois(%s)'%s)
                d = int(kash.eval('%s.ext1'%G.name()))
                n = int(kash.eval('%s.ext2'%G.name()))
                return TransitiveGroup(d, n)
            except RuntimeError, msg:
                raise NotImplementedError, (str(msg) + "\nSorry, " +
                    "computation of Galois groups of fields of degree " +
                    "bigger than 11 is not yet implemented.  Try installing " +
                    "the optional free (closed source) KASH package, which " +
                    "supports degrees up to 23, or use algorithm='magma' if " +
                    "you have magma.")

        elif algorithm == 'magma':
            from sage.interfaces.all import magma
            X = magma(self).GaloisGroup()
            try:
                n, d = X.TransitiveGroupIdentification(nvals=2)
                d = int(d)
                n = int(n)
            except RuntimeError, msg:
                raise RuntimeError, (str(msg) + "\nUnable to lookup " +
                    "description of Galois group as a transitive " +
                    "group.\n%s" %X)
            return TransitiveGroup(d, n)

        else:
            raise ValueError, "Algorithm %s not supported." %algorithm

    def factor_mod(self, p):
        """
        Returns the factorization of self modulo the prime ``p``.

        Assumes that the degree of this polynomial is at least one, and raises
        a ``ValueError`` otherwise.

        INPUT:

        -  ``p`` - Prime number

        OUTPUT:

        -  Factorization of this polynomial  modulo ``p``

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: (x^5 + 17*x^3 + x+ 3).factor_mod(3)
            x * (x^2 + 1)^2
            sage: (x^5 + 2).factor_mod(5)
            (x + 2)^5
        """
        from sage.rings.finite_rings.constructor import FiniteField

        p = Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"

        if self.degree() < 1:
            raise ValueError, "The polynomial must have degree at least 1"

        G = self._pari_().factormod(p)
        K = FiniteField(p)
        R = K[self.parent().variable_name()]
        return R(1)._factor_pari_helper(G, unit=R(self).leading_coefficient())

    def factor_padic(self, p, prec=10):
        """
        Return the `p`-adic factorization of this polynomial to the given
        precision.

        INPUT:

        -  ``p`` - Prime number

        -  ``prec`` - Integer; the precision

        OUTPUT:

        - factorization of ``self`` viewed as a `p`-adic polynomial

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^3 - 2
            sage: f.factor_padic(2)
            (1 + O(2^10))*x^3 + (2 + 2^2 + 2^3 + 2^4 + 2^5 + 2^6 + 2^7 + 2^8 + 2^9 + O(2^10))
            sage: f.factor_padic(3)
            (1 + O(3^10))*x^3 + (1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10))
            sage: f.factor_padic(5)
            ((1 + O(5^10))*x + (2 + 4*5 + 2*5^2 + 2*5^3 + 5^4 + 3*5^5 + 4*5^7 + 2*5^8 + 5^9 + O(5^10))) * ((1 + O(5^10))*x^2 + (3 + 2*5^2 + 2*5^3 + 3*5^4 + 5^5 + 4*5^6 + 2*5^8 + 3*5^9 + O(5^10))*x + (4 + 5 + 2*5^2 + 4*5^3 + 4*5^4 + 3*5^5 + 3*5^6 + 4*5^7 + 4*5^9 + O(5^10)))

        The input polynomial is considered to have "infinite" precision,
        therefore the `p`-adic factorization of the polynomial is not
        the same as first coercing to `Q_p` and then factoring
        (see also :trac:`15422`)::

            sage: f = x^2 - 3^6
            sage: f.factor_padic(3,5)
            ((1 + O(3^5))*x + (3^3 + O(3^5))) * ((1 + O(3^5))*x + (2*3^3 + 2*3^4 + O(3^5)))
            sage: f.change_ring(Qp(3,5)).factor()
            Traceback (most recent call last):
            ...
            PrecisionError: p-adic factorization not well-defined since the discriminant is zero up to the requestion p-adic precision

        A more difficult example::

            sage: f = 100 * (5*x + 1)^2 * (x + 5)^2
            sage: f.factor_padic(5, 10)
            (4*5^4 + O(5^14)) * ((1 + O(5^9))*x + (5^-1 + O(5^9)))^2 * ((1 + O(5^10))*x + (5 + O(5^10)))^2

        Try some bogus inputs::

            sage: f.factor_padic(3,-1)
            Traceback (most recent call last):
            ...
            ValueError: prec_cap must be non-negative.
            sage: f.factor_padic(6,10)
            Traceback (most recent call last):
            ...
            ValueError: p must be prime
            sage: f.factor_padic('hello', 'world')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x (=hello) to an integer
        """
        from sage.rings.padics.factory import Qp

        p = Integer(p)
        prec = Integer(prec)

        # Parent field for coefficients and polynomial
        K = Qp(p, prec, type='capped-rel')
        R = K[self.parent().variable_name()]

        # Factor the *exact* polynomial using factorpadic()
        G = self._pari_with_name().factorpadic(p, prec)

        from sage.rings.polynomial.padics.polynomial_padic import _pari_padic_factorization_to_sage
        return _pari_padic_factorization_to_sage(G, R, self.leading_coefficient())

    def hensel_lift(self, p, e):
        r"""
        Assuming that this polynomial factors modulo `p` into distinct
        monic factors, computes the Hensel lifts of these factors modulo
        `p^e`. We assume that ``self`` has integer coefficients.

        Returns an empty list if this polynomial has degree less than one.

        INPUT:

        -  ``p`` - Prime number; coerceable to ``Integer``
        -  ``e`` - Exponent; coerceable to ``Integer``

        OUTPUT:

        -  Hensel lifts; list of polynomials over `\ZZ / p^e \ZZ`

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R((x-1)*(x+1)).hensel_lift(7, 2)
            [x + 1, x + 48]

        If the input polynomial `f` is not monic, we get a factorization of
        `f / lc(f)`::

            sage: R(2*x^2 - 2).hensel_lift(7, 2)
            [x + 1, x + 48]

        TESTS::

            sage: R.<x> = QQ[]
            sage: R(0).hensel_lift(7, 2)
            []
            sage: R(x).hensel_lift(7, 2)
            [x]
            sage: R(x-1).hensel_lift(7, 2)
            [x + 48]
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

        p = Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"
        e = Integer(e)
        if e < 1:
            raise ValueError, "e must be at least 1"

        # The relevant PARI method doesn't seem to play well with constant and
        # linear polynomials, so we handle these separately.
        #
        if self.degree() < 1:
            return [ ]
        elif self.degree() == 1:
            R = IntegerModRing(p**e)
            S = R[self.parent().variable_name()]
            return [S(self)]

        F = self.factor_mod(p)
        y = []
        for g, n in F:
            if n > 1:
                raise ArithmeticError, ("The polynomial must be square free " +
                    "modulo p")
            y.append(self.parent()(g))
        H = self._pari_().polhensellift(y, p, e)
        R = IntegerModRing(p**e)
        S = R[self.parent().variable_name()]
        return [S(m) for m in H]

    def discriminant(self):
        r"""
        Returns the discriminant of this polynomial.

        The discriminant `R_n` is defined as

        .. math::

            R_n = a_n^{2 n-2} \prod_{1 \le i < j \le n} (r_i - r_j)^2,

        where `n` is the degree of this polynomial, `a_n` is the leading
        coefficient and the roots over `\QQbar` are `r_1, \ldots, r_n`.

        The discriminant of constant polynomials is defined to be 0.

        OUTPUT:

        -  Discriminant, an element of the base ring of the polynomial ring

        .. note::

            Note the identity `R_n(f) := (-1)^(n (n-1)/2) R(f,f') a_n^(n-k-2)`,
            where `n` is the degree of this polynomial, `a_n` is the leading
            coefficient, `f'` is the derivative of `f`, and `k` is the degree
            of `f'`.  Calls :meth:`.resultant`.

        ALGORITHM:

        Use PARI.

        EXAMPLES:

        In the case of elliptic curves in special form, the discriminant is
        easy to calculate::

            sage: R.<t> = QQ[]
            sage: f = t^3 + t + 1
            sage: d = f.discriminant(); d
            -31
            sage: d.parent() is QQ
            True
            sage: EllipticCurve([1, 1]).discriminant() / 16
            -31

        ::

            sage: R.<t> = QQ[]
            sage: f = 2*t^3 + t + 1
            sage: d = f.discriminant(); d
            -116

        ::

            sage: R.<t> = QQ[]
            sage: f = t^3 + 3*t - 17
            sage: f.discriminant()
            -7911

        TESTS::

            sage: R.<t> = QQ[]
            sage: R(0).discriminant()
            0
            sage: R(2/3).discriminant()
            0
            sage: (t + 1/2).discriminant()
            1
        """
        return QQ(self._pari_().poldisc())

    # Alias for discriminant
    disc = discriminant


