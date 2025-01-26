# distutils: libraries = NTL_LIBRARIES gmp
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
r"""
Univariate polynomials over `\QQ` implemented via FLINT

AUTHOR:

- Sebastian Pancratz
"""

# ****************************************************************************
#          Copyright (C) 2010 Sebastian Pancratz <sfp@pancratz.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_str, sig_off

from cpython.long cimport PyLong_AsLong
from sage.arith.long cimport pyobject_to_long

from sage.libs.flint.acb cimport acb_div_fmpz
from sage.libs.flint.arb cimport arb_div_fmpz
from sage.libs.flint.arb_fmpz_poly cimport _arb_fmpz_poly_evaluate_arb, _arb_fmpz_poly_evaluate_acb
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpq cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.fmpq_poly cimport *
from sage.libs.flint.fmpq_poly_sage cimport *
from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *

from cypari2.gen import Gen as pari_gen

from sage.rings.complex_arb cimport ComplexBall
from sage.rings.integer cimport Integer, smallInteger
from sage.rings.integer_ring import ZZ
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.rational cimport Rational
from sage.rings.rational_field import QQ
from sage.rings.real_arb cimport RealBall
from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint

from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.structure.element import coerce_binop

from sage.misc.cachefunc import cached_method

cdef inline bint _do_sig(fmpq_poly_t op) noexcept:
    """
    Return 1 when signal handling should be carried out for an operation
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
    # Issue #12173: check that the degree is greater than 1000 before computing
    # the max limb size
    return (fmpq_poly_length(op) > 0 and
            (fmpq_poly_degree(op) > 1000 or
                sage_fmpq_poly_max_limbs(op) > 1))

cdef class Polynomial_rational_flint(Polynomial):
    """
    Univariate polynomials over the rationals, implemented via FLINT.

    Internally, we represent rational polynomial as the quotient of an integer
    polynomial and a positive denominator which is coprime to the content of
    the numerator.

    TESTS::

        sage: f = QQ['x'].random_element()
        sage: from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint
        sage: isinstance(f, Polynomial_rational_flint)
        True

    .. automethod:: _add_
    .. automethod:: _sub_
    .. automethod:: _lmul_
    .. automethod:: _rmul_
    .. automethod:: _mul_
    .. automethod:: _mul_trunc_
    """

    ###########################################################################
    # Allocation & initialisation                                             #
    ###########################################################################

    cdef Polynomial_rational_flint _new(self):
        """
        Quickly create a new polynomial object in this class.

        OUTPUT: polynomial of type ``Polynomial_rational_flint``

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 2/3*t^2
            sage: g = -1/2*t + 2
            sage: f + g           # indirect doctest
            2/3*t^2 - 1/2*t + 2
        """
        cdef Polynomial_rational_flint res = Polynomial_rational_flint.__new__(Polynomial_rational_flint)
        res._parent = self._parent
        res._is_gen = 0
        return res

    cpdef Polynomial _new_constant_poly(self, x, Parent P):
        r"""
        Quickly create a new constant polynomial with value x in parent P.

        ASSUMPTION:

        x must be a rational or convertible to an int.

        EXAMPLES::

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
        cdef Polynomial_rational_flint res = Polynomial_rational_flint.__new__(Polynomial_rational_flint)
        res._parent = P
        res._is_gen = <char>0
        if isinstance(x, int):
            fmpq_poly_set_si(res._poly, <int> x)

        elif isinstance(x, Integer):
            fmpq_poly_set_mpz(res._poly, (<Integer> x).value)

        elif isinstance(x, Rational):
            fmpq_poly_set_mpq(res._poly, (<Rational> x).value)

        else:
            fmpq_poly_set_si(res._poly, int(x))
        return res

    def __cinit__(self):
        """
        Initialises the underlying data structure.

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 2/3 * t - 7  #indirect doctest
        """
        fmpq_poly_init(self._poly)

    def __dealloc__(self):
        """
        Deallocate the underlying data structure.

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 1/3 * t
            sage: del f
        """
        fmpq_poly_clear(self._poly)

    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        Initialises the associated data for the polynomial ``self``.

        INPUT:

        - ``parent`` -- polynomial ring, the parent of ``self``
        - ``x`` -- data for the new polynomial ``self``, e.g. a polynomial, an
          integer, a rational, a list of rationals, a dictionary with keys
          the degrees and the rational coefficients, etc (default: ``None``)
        - ``check`` -- whether the integrity of the data needs to be verified,
          largely ignored by this method (default: ``True``)
        - ``is_gen`` -- whether ``self`` shall be initialised as the generator of
          the parent polynomial ring
        - ``construct`` -- whether the element shall always be constructed
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
        cdef fmpq_t q

        Polynomial.__init__(self, parent, is_gen=is_gen)

        if is_gen:
            fmpq_poly_set_coeff_si(self._poly, 1, 1)

        elif isinstance(x, Polynomial_rational_flint):
            fmpq_poly_set(self._poly, (<Polynomial_rational_flint> x)._poly)

        elif isinstance(x, int):
            fmpq_poly_set_si(self._poly, <int> x)

        elif isinstance(x, Integer):
            fmpq_poly_set_mpz(self._poly, (<Integer> x).value)

        elif isinstance(x, Rational):
            fmpq_poly_set_mpq(self._poly, (<Rational> x).value)

        elif isinstance(x, (list, tuple)):

            if len(x) == 0:
                return
            elif len(x) == 1:
                Polynomial_rational_flint.__init__(self, parent, x[0],
                                check=check, is_gen=False, construct=construct)
                return

            L1 = [e if isinstance(e, Rational) else Rational(e) for e in x]
            n  = <unsigned long> len(x)
            sig_on()
            fmpq_poly_fit_length(self._poly, n)
            for deg from 0 <= deg < n:
                fmpq_init_set_readonly(q, (<Rational> L1[deg]).value)
                fmpq_poly_set_coeff_fmpq(self._poly, deg, q)
                fmpq_clear_readonly(q)
            sig_off()

#           deg = 0
#           for e in x:
#               c = Rational(e)
#               fmpq_poly_set_coeff_mpq(self._poly, deg, c.value)
#               deg += 1

        elif isinstance(x, dict):
            for deg, e in x.iteritems():
                c = Rational(e)
                fmpq_poly_set_coeff_mpq(self._poly, deg, c.value)

        elif isinstance(x, pari_gen):
            k = self._parent.base_ring()
            x = [k(w) for w in x.list()]
            Polynomial_rational_flint.__init__(self, parent, x, check=True,
                                             is_gen=False, construct=construct)

        elif isinstance(x, Polynomial_integer_dense_flint):
            fmpq_poly_set_fmpz_poly(self._poly, (<Polynomial_integer_dense_flint>x)._poly)

        elif isinstance(x, Polynomial):
            k = self._parent.base_ring()
            x = [k(w) for w in list(x)]
            Polynomial_rational_flint.__init__(self, parent, x, check=True,
                                             is_gen=False, construct=construct)

        elif isinstance(x, FractionFieldElement) and (x.parent().base() is parent or x.parent().base() == parent) and x.denominator() == 1:
            x = x.numerator()
            Polynomial_rational_flint.__init__(self, parent, x, check=check,
                                            is_gen=is_gen, construct=construct)

        else:
            x = parent.base_ring()(x)
            Polynomial_rational_flint.__init__(self, parent, x, check=check,
                                            is_gen=is_gen, construct=construct)

    def __reduce__(self):
        """
        This is used when pickling polynomials.

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 2/3 * t^2 + 1
            sage: r = f.__reduce__(); r
            (<class 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint'>, (Univariate Polynomial Ring in t over Rational Field, [1, 0, 2/3], False, False))
            sage: r[0](*r[1])
            2/3*t^2 + 1
            sage: loads(dumps(f)) == f
            True
        """
        return (Polynomial_rational_flint,
               (self.parent(), self.list(), False, self.is_gen()))

    def __copy__(self):
        """
        Return a copy of ``self``.

        TESTS::

            sage: R.<t> = QQ[]
            sage: f = 4/5 * t^3 - 1/17
            sage: copy(f) == f
            True
        """
        cdef Polynomial_rational_flint res = self._new()
        fmpq_poly_set(res._poly, self._poly)
        return res

    def _singular_(self, singular=None):
        """
        Return a Singular representation of ``self``.

        INPUT:

        - ``singular`` -- Singular interpreter (default: default interpreter)

        EXAMPLES::

            sage: P.<x> = PolynomialRing(QQ)
            sage: f = 3*x^2 + 2*x + 5
            sage: singular(f)                                                           # needs sage.libs.singular
            3*x^2+2*x+5
        """
        if singular is None:
            from sage.interfaces.singular import singular
        self._parent._singular_(singular).set_ring()  # Expensive!
        return singular(self._singular_init_())

    cpdef list list(self, bint copy=True):
        """
        Return a list with the coefficients of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2/2 + t^3/3 + t^4/4
            sage: f.list()
            [1, 1, 1/2, 1/3, 1/4]
            sage: g = R(0)
            sage: g.list()
            []
        """
        cdef unsigned long length = fmpq_poly_length(self._poly)
        return [self.get_unsafe(n) for n in range(length)]

    ###########################################################################
    # Basis access                                                            #
    ###########################################################################

    def degree(self):
        """
        Return the degree of ``self``.

        By convention, the degree of the zero polynomial is `-1`.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2/2 + t^3/3 + t^4/4
            sage: f.degree()
            4
            sage: g = R(0)
            sage: g.degree()
            -1

        TESTS::

            sage: type(f.degree())
            <class 'sage.rings.integer.Integer'>
        """
        return smallInteger(fmpq_poly_degree(self._poly))

    cdef get_unsafe(self, Py_ssize_t n):
        """
        Return the `n`-th coefficient of ``self``.

        INPUT:

        - ``n`` -- degree of the monomial whose coefficient is to be
          returned

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2/2 + t^3/3 + t^4/4
            sage: f[-1], f[0], f[3], f[5]           # indirect doctest
            (0, 1, 1/3, 0)
            sage: f[:3]                             # indirect doctest
            1/2*t^2 + t + 1
        """
        cdef Rational z = Rational.__new__(Rational)
        fmpq_poly_get_coeff_mpq(z.value, self._poly, n)
        return z

    cpdef _unsafe_mutate(self, unsigned long n, value):
        """
        Set the `n`-th coefficient of ``self`` to ``value``.

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
        cdef bint do_sig = _do_sig(self._poly)
        cdef fmpz_t tmpfz

        if isinstance(value, int):
            if do_sig: sig_str("FLINT exception")
            fmpq_poly_set_coeff_si(self._poly, n, value)
            if do_sig: sig_off()
        elif isinstance(value, Integer):
            if do_sig: sig_str("FLINT exception")
            fmpz_init_set_readonly(tmpfz, (<Integer> value).value)
            fmpq_poly_set_coeff_fmpz(self._poly, n, tmpfz)
            fmpz_clear_readonly(tmpfz)
            if do_sig: sig_off()
        elif isinstance(value, Rational):
            if do_sig: sig_str("FLINT exception")
            fmpq_poly_set_coeff_mpq(self._poly, n, (<Rational> value).value)
            if do_sig: sig_off()
        else:
            value = Rational(value)
            if do_sig: sig_str("FLINT exception")
            fmpq_poly_set_coeff_mpq(self._poly, n, (<Rational> value).value)
            if do_sig: sig_off()

    def __call__(self, *x, **kwds):
        """
        Call this polynomial with the given parameters, which can be
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

        TESTS:

            sage: t(-sys.maxsize-1r) == t(-sys.maxsize-1)
            True
            sage: (t/3)(RealBallField(100)(1))
            [0.33333333333333333333333333333...]
            sage: (t/3)(ComplexBallField(10)(1+i))                                      # needs sage.symbolic
            [0.33...] + [0.33...]*I
        """
        cdef Polynomial_rational_flint f
        cdef Rational r
        cdef fmpz_t tmpfz
        cdef fmpq_t tmpfq, tmpfq1
        cdef RealBall arb_a, arb_z
        cdef ComplexBall acb_a, acb_z

        if len(x) == 1:
            a = x[0]
            if isinstance(a, Polynomial_rational_flint):
                f = (<Polynomial_rational_flint> a)._new()
                sig_str("FLINT exception")
                fmpq_poly_compose(f._poly, self._poly,
                    (<Polynomial_rational_flint> a)._poly)
                sig_off()
                return f
            elif isinstance(a, Rational):
                r = Rational.__new__(Rational)
                sig_str("FLINT exception")
                fmpq_init_set_readonly(tmpfq, (<Rational> a).value)
                fmpq_init(tmpfq1)
                fmpq_poly_evaluate_fmpq(tmpfq1, self._poly, tmpfq)
                fmpq_get_mpq(r.value, tmpfq1)
                fmpq_clear(tmpfq1)
                fmpq_clear_readonly(tmpfq)
                sig_off()
                return r
            elif isinstance(a, Integer):
                r = Rational.__new__(Rational)
                sig_str("FLINT exception")
                fmpz_init_set_readonly(tmpfz, (<Integer> a).value)
                fmpq_init(tmpfq)
                fmpq_poly_evaluate_fmpz(tmpfq, self._poly, tmpfz)
                fmpq_get_mpq(r.value, tmpfq)
                fmpq_clear(tmpfq)
                fmpz_clear_readonly(tmpfz)
                sig_off()
                return r
            elif isinstance(a, int):
                r = Rational.__new__(Rational)
                sig_str("FLINT exception")
                fmpz_init(tmpfz)
                fmpq_init(tmpfq)
                fmpz_set_si(tmpfz, PyLong_AsLong(a))
                fmpq_poly_evaluate_fmpz(tmpfq, self._poly, tmpfz)
                fmpq_get_mpq(r.value, tmpfq)
                fmpq_clear(tmpfq)
                fmpz_clear(tmpfz)
                sig_off()
                return r
            if isinstance(a, RealBall):
                arb_a = <RealBall> a
                arb_z = arb_a._new()
                sig_on()
                _arb_fmpz_poly_evaluate_arb(arb_z.value, fmpq_poly_numref(self._poly),
                        fmpq_poly_length(self._poly), arb_a.value, arb_a._parent._prec)
                arb_div_fmpz(arb_z.value, arb_z.value, fmpq_poly_denref(self._poly), arb_a._parent._prec)
                sig_off()
                return arb_z
            if isinstance(a, ComplexBall):
                acb_a = <ComplexBall> a
                acb_z = acb_a._new()
                sig_on()
                _arb_fmpz_poly_evaluate_acb(acb_z.value, fmpq_poly_numref(self._poly),
                        fmpq_poly_length(self._poly), acb_a.value, acb_a._parent._prec)
                acb_div_fmpz(acb_z.value, acb_z.value, fmpq_poly_denref(self._poly), acb_a._parent._prec)
                sig_off()
                return acb_z

        return Polynomial.__call__(self, *x, **kwds)

    cpdef Polynomial truncate(self, long n):
        """
        Return ``self`` truncated modulo `t^n`.

        INPUT:

        - ``n`` -- the power of `t` modulo which ``self`` is truncated

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

        if (n >= fmpq_poly_length(self._poly)):
            return self
        else:
            res = self._new()
            if n > 0:
                do_sig = _do_sig(self._poly)
                if do_sig: sig_str("FLINT exception")
                fmpq_poly_get_slice(res._poly, self._poly, 0, n)
                if do_sig: sig_off()
            return res

    def reverse(self, degree=None):
        """
        Reverse the coefficients of this polynomial (thought of as a polynomial
        of degree ``degree``).

        INPUT:

        - ``degree`` -- ``None`` or integral value that fits in an ``unsigned
          long`` (default: degree of ``self``); if specified, truncate or zero
          pad the list of coefficients to this degree before reversing it

        EXAMPLES:

        We first consider the simplest case, where we reverse all coefficients
        of a polynomial and obtain a polynomial of the same degree::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2 / 2 + t^3 / 3 + t^4 / 4
            sage: f.reverse()
            t^4 + t^3 + 1/2*t^2 + 1/3*t + 1/4

        Next, an example where the returned polynomial has lower degree because
        the original polynomial has low coefficients equal to zero::

            sage: R.<t> = QQ[]
            sage: f = 3/4*t^2 + 6*t^7
            sage: f.reverse()
            3/4*t^5 + 6

        The next example illustrates the passing of a value for ``degree`` less
        than the length of ``self``, notationally resulting in truncation prior to
        reversing::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2 / 2 + t^3 / 3 + t^4 / 4
            sage: f.reverse(2)
            t^2 + t + 1/2

        Now we illustrate the passing of a value for ``degree`` greater than
        the length of ``self``, notationally resulting in zero padding at the top
        end prior to reversing::

            sage: R.<t> = QQ[]
            sage: f = 1 + t + t^2 / 2 + t^3 / 3
            sage: f.reverse(4)
            t^4 + t^3 + 1/2*t^2 + 1/3*t

        TESTS:

        We illustrate two ways in which the interpretation of ``degree`` as an
        unsigned long int may fail.  Firstly, an integral value which is
        too large, yielding an :exc:`OverflowError`::

            sage: R.<t> = QQ[]
            sage: f = 1 + t/2
            sage: f.reverse(2**64 - 1)
            Traceback (most recent call last):
            ...
            OverflowError:... int too large to convert...

        Secondly, a value which cannot be converted to an integral value,
        resulting in a ValueError::

            sage: R.<t> = QQ[]
            sage: f = 1 + t/2
            sage: f.reverse(I)
            Traceback (most recent call last):
            ...
            ValueError: degree must be convertible to long

        We check that this specialized implementation is compatible with the
        generic one::

            sage: all((t + 2*t^2).reverse(degree=d)
            ....:     == Polynomial.reverse(t + 2*t^2, degree=d)
            ....:     for d in [None, 0, 1, 2, 3, 4, 5])
            True
        """
        cdef unsigned long len
        cdef Polynomial_rational_flint res
        cdef bint do_sig

        if degree is None:
            len = fmpq_poly_length(self._poly)
        else:
            try:
                len = <unsigned long> (degree + 1)
            except (TypeError, ValueError):
                raise ValueError('degree must be convertible to long')

        res = self._new()
        do_sig = _do_sig(self._poly)
        if do_sig: sig_str("FLINT exception")
        fmpq_poly_reverse(res._poly, self._poly, len)
        if do_sig: sig_off()
        return res

    def revert_series(self, n):
        r"""
        Return a polynomial `f` such that ``f(self(x)) = self(f(x)) = x mod x^n``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = t - t^3/6 + t^5/120
            sage: f.revert_series(6)
            3/40*t^5 + 1/6*t^3 + t

            sage: f.revert_series(-1)
            Traceback (most recent call last):
            ValueError: argument n must be a nonnegative integer, got -1

            sage: g = - t^3/3 + t^5/5
            sage: g.revert_series(6)
            Traceback (most recent call last):
            ...
            ValueError: self must have constant coefficient 0 and a unit for coefficient t^1
        """

        cdef Polynomial_rational_flint res = self._new()
        cdef unsigned long m
        if n < 0:
            raise ValueError("argument n must be a nonnegative integer, got {}".format(n))
        m = n
        if not self[0].is_zero() or not self[1].is_unit():
            raise ValueError("self must have constant coefficient 0 and a unit for coefficient {}^1".format(self.parent().gen()))

        sig_str("FLINT exception")
        fmpq_poly_revert_series(res._poly, self._poly, m)
        sig_off()

        return res

    ###########################################################################
    # Comparisons                                                             #
    ###########################################################################

    cpdef bint is_zero(self) except -1:
        """
        Return whether or not ``self`` is the zero polynomial.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 - t + 1/2*t^2 - 1/3*t^3
            sage: f.is_zero()
            False
            sage: R(0).is_zero()
            True
        """
        return fmpq_poly_is_zero(self._poly)

    cpdef bint is_one(self) except -1:
        r"""
        Return whether or not this polynomial is one.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R([0,1]).is_one()
            False
            sage: R([1]).is_one()
            True
            sage: R([0]).is_one()
            False
            sage: R([-1]).is_one()
            False
            sage: R([1,1]).is_one()
            False
        """
        return fmpq_poly_is_one(self._poly)

    def __bool__(self):
        """
        Return whether or not ``self`` is nonzero.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1 - t + 1/2*t^2 - 1/3*t^3
            sage: bool(f)
            True
            sage: bool(R(0))
            False
        """
        return not fmpq_poly_is_zero(self._poly)

    ###########################################################################
    # Shifting                                                                #
    ###########################################################################

    def __lshift__(self, n):
        """
        Notationally multiply ``self`` by `t^n`.

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

        if k == 0 or fmpq_poly_is_zero(f._poly):
            return self
        else:
            res = f._new()
            do_sig = fmpq_poly_length(f._poly) > 5000 or n > 5000

            if do_sig: sig_str("FLINT exception")
            fmpq_poly_shift_left(res._poly, f._poly, k)
            if do_sig: sig_off()
            return res

    def __rshift__(self, n):
        """
        Notationally return the quotient of Euclidean division of ``self``
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

        if k == 0 or fmpq_poly_is_zero(f._poly):
            return self
        else:
            res = f._new()
            do_sig = _do_sig(f._poly)

            if do_sig: sig_str("FLINT exception")
            fmpq_poly_shift_right(res._poly, f._poly, k)
            if do_sig: sig_off()
            return res

    ###########################################################################
    # Arithmetic                                                              #
    ###########################################################################

    cpdef _add_(self, right):
        """
        Return the sum of two rational polynomials.

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
        cdef bint do_sig = _do_sig(self._poly) or _do_sig(op2._poly)

        if do_sig: sig_str("FLINT exception")
        fmpq_poly_add(res._poly, self._poly, op2._poly)
        if do_sig: sig_off()
        return res

    cpdef _sub_(self, right):
        """
        Return the difference of two rational polynomials.

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
        """
        cdef Polynomial_rational_flint op2 = <Polynomial_rational_flint> right
        cdef Polynomial_rational_flint res = self._new()
        cdef bint do_sig = _do_sig(self._poly) or _do_sig(op2._poly)

        if do_sig: sig_str("FLINT exception")
        fmpq_poly_sub(res._poly, self._poly, op2._poly)
        if do_sig: sig_off()
        return res

    cpdef _neg_(self):
        """
        Return the difference of two rational polynomials.

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
        cdef bint do_sig = _do_sig(self._poly)

        if do_sig: sig_str("FLINT exception")
        fmpq_poly_neg(res._poly, self._poly)
        if do_sig: sig_off()
        return res

    @coerce_binop
    def quo_rem(self, right):
        r"""
        Return the quotient and remainder of the Euclidean division of
        ``self`` and ``right``.

        Raises a :exc:`ZeroDivisionError` if ``right`` is zero.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = R.random_element(2000)
            sage: g = R.random_element(1000)
            sage: q, r = f.quo_rem(g)
            sage: f == q*g + r
            True
        """
        if right.is_zero():
            raise ZeroDivisionError("division by zero polynomial")
        if self.is_zero():
            return self, self

        cdef Polynomial_rational_flint qq = self._new()
        cdef Polynomial_rational_flint rr = self._new()

        sig_str("FLINT exception")
        fmpq_poly_divrem(qq._poly, rr._poly, self._poly,
                         (<Polynomial_rational_flint> right)._poly)
        sig_off()
        return qq, rr

    @coerce_binop
    def gcd(self, right):
        r"""
        Return the (monic) greatest common divisor of ``self`` and ``right``.

        Corner cases:  if ``self`` and ``right`` are both zero, returns zero.  If
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

        sig_str("FLINT exception")
        fmpq_poly_gcd(res._poly, self._poly,
                (<Polynomial_rational_flint> right)._poly)
        sig_off()
        return res

    @coerce_binop
    def lcm(self, right):
        r"""
        Return the monic (or zero) least common multiple of ``self`` and ``right``.

        Corner cases:  if either of ``self`` and ``right`` are zero, returns zero.
        This behaviour is ensures that the relation `\lcm(a,b)\cdot \gcd(a,b) = a\cdot b`
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

        sig_str("FLINT exception")
        fmpq_poly_lcm(res._poly, self._poly,
                      (<Polynomial_rational_flint> right)._poly)
        sig_off()
        return res

    @coerce_binop
    def xgcd(self, right):
        r"""
        Return polynomials `d`, `s`, and `t` such that ``d == s * self + t * right``,
        where `d` is the (monic) greatest common divisor of ``self`` and ``right``.
        The choice of `s` and `t` is not specified any further.

        Corner cases:  if ``self`` and ``right`` are zero, returns zero polynomials.
        Otherwise, if only ``self`` is zero, returns ``(d, s, t) = (right, 0, 1)`` up
        to normalisation, and similarly if only ``right`` is zero.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 2/3 + 3/4 * t - t^2
            sage: g = -3 + 1/7 * t
            sage: f.xgcd(g)
            (1, -12/5095, -84/5095*t - 1701/5095)

        TESTS:

        The following example used to crash (cf. :issue:`11771`)::

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

        sig_str("FLINT exception")
        fmpq_poly_xgcd(d._poly, s._poly, t._poly, self._poly, (<Polynomial_rational_flint>right)._poly)
        sig_off()
        return d, s, t

    cpdef _mul_(self, right):
        """
        Return the product of ``self`` and ``right``.

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
        cdef bint do_sig = _do_sig(self._poly) or _do_sig(op2._poly)

        if do_sig: sig_str("FLINT exception")
        fmpq_poly_mul(res._poly, self._poly, op2._poly)
        if do_sig: sig_off()
        return res

    cpdef Polynomial _mul_trunc_(self, Polynomial right, long n):
        r"""
        Truncated multiplication.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: p1 = 1/2 - 3*x + 2/7*x**3
            sage: p2 = x + 2/5*x**5 + x**7
            sage: p1._mul_trunc_(p2, 5)
            2/7*x^4 - 3*x^2 + 1/2*x
            sage: (p1*p2).truncate(5)
            2/7*x^4 - 3*x^2 + 1/2*x

            sage: p1._mul_trunc_(p2, 1)
            0
            sage: p1._mul_trunc_(p2, 0)
            Traceback (most recent call last):
            ...
            ValueError: n must be > 0

        ALGORITHM:

        Call the FLINT method ``fmpq_poly_mullow``.
        """
        cdef Polynomial_rational_flint op2 = <Polynomial_rational_flint> right
        cdef Polynomial_rational_flint res = self._new()
        cdef bint do_sig = _do_sig(self._poly) or _do_sig(op2._poly)

        if n <= 0:
            raise ValueError("n must be > 0")

        if do_sig: sig_str("FLINT exception")
        fmpq_poly_mullow(res._poly, self._poly, op2._poly, n)
        if do_sig: sig_off()
        return res

    cpdef _rmul_(self, Element left):
        r"""
        Return ``left * self``, where ``left`` is a rational number.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 3/2*t^3 - t + 1/3
            sage: 6 * f                  # indirect doctest
            9*t^3 - 6*t + 2
        """
        cdef Polynomial_rational_flint res = self._new()
        cdef bint do_sig = _do_sig(self._poly)

        if do_sig: sig_str("FLINT exception")
        fmpq_poly_scalar_mul_mpq(res._poly, self._poly,
                                 (<Rational> left).value)
        if do_sig: sig_off()
        return res

    cpdef _lmul_(self, Element right):
        r"""
        Return ``self * right``, where ``right`` is a rational number.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 3/2*t^3 - t + 1/3
            sage: f * 6                   # indirect doctest
            9*t^3 - 6*t + 2
        """
        cdef Polynomial_rational_flint res = self._new()
        cdef bint do_sig = _do_sig(self._poly)

        if do_sig: sig_str("FLINT exception")
        fmpq_poly_scalar_mul_mpq(res._poly, self._poly,
                                 (<Rational> right).value)
        if do_sig: sig_off()
        return res

    def __pow__(Polynomial_rational_flint self, exp, mod):
        """
        Return ``self`` raised to the power of ``exp``.

        The corner case of ``exp == 0`` is handled by returning the constant
        polynomial 1.  Note that this includes the case ``0^0 == 1``.

        This method only supports integral values for ``exp`` that fit into
        a signed long int (except when this is a constant polynomial).

        INPUT:

        - ``exp`` -- exponent

        OUTPUT: polynomial; this polynomial raised to the power of ``exp``

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = 1/2 + 2*t - t^2/3
            sage: f^0
            1
            sage: f^3
            -1/27*t^6 + 2/3*t^5 - 23/6*t^4 + 6*t^3 + 23/4*t^2 + 3/2*t + 1/8
            sage: f^(-3)
            -27/(t^6 - 18*t^5 + 207/2*t^4 - 162*t^3 - 621/4*t^2 - 81/2*t - 27/8)

        TESTS::

            sage: R.<t> = QQ[]
            sage: t^0
            1
            sage: R(0)^0
            1
            sage: R(0)^(-1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: negative exponent in power of zero

        We verify the checking of the exponent::

            sage: R.<t> = QQ[]
            sage: (1 + t)^(2/3)
            Traceback (most recent call last):
            ...
            ValueError: not a 3rd power
            sage: (1 + t)^(2^63)
            Traceback (most recent call last):
            ...
            OverflowError: Sage Integer too large to convert to C long

        FLINT memory errors do not crash Sage (:issue:`17629`)::

            sage: t^(sys.maxsize//2)
            Traceback (most recent call last):
            ...
            RuntimeError: FLINT exception

        Flush the output buffer to get rid of stray output -- see
        :issue:`28649`::

            sage: from sage.misc.misc_c import cyflush
            sage: cyflush()
            ...

        Test fractional powers (:issue:`20086`)::

            sage: P.<R> = QQ[]
            sage: (1/27*R^3 + 2/3*R^2 + 4*R + 8)^(2/3)
            1/9*R^2 + 4/3*R + 4
            sage: _.parent()
            Univariate Polynomial Ring in R over Rational Field
            sage: P(1/4)^(1/2)
            1/2
            sage: _.parent()
            Univariate Polynomial Ring in R over Rational Field

            sage: (R+2)^(2/5)
            Traceback (most recent call last):
            ...
            ValueError: not a 5th power

            sage: P(1/3)^(1/2)
            Traceback (most recent call last):
            ...
            ValueError: not a perfect 2nd power
            sage: P(4)^P(1/2)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Univariate Polynomial
            Ring in R over Rational Field to Rational Field
            sage: (R + 1)^P(2)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Univariate Polynomial
            Ring in R over Rational Field to Rational Field
            sage: (R + 1)^R
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Univariate Polynomial
            Ring in R over Rational Field to Rational Field
            sage: 2^R
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Univariate Polynomial
            Ring in R over Rational Field to Rational Field

        Check that using third argument raises an error::

            sage: R.<x> = PolynomialRing(QQ)
            sage: pow(x, 2, x)
            Traceback (most recent call last):
            ...
            NotImplementedError: pow() with a modulus is not implemented for this ring
        """
        cdef Polynomial_rational_flint res
        cdef long n

        if mod is not None:
            raise NotImplementedError(
                "pow() with a modulus is not implemented for this ring"
            )

        try:
            n = pyobject_to_long(exp)
        except TypeError:
            r = QQ.coerce(exp)
            num = r.numerator()
            den = r.denominator()

            if fmpq_poly_degree(self._poly) == 0:
                return self.parent()(self[0].nth_root(den) ** num)

            return self.nth_root(den) ** num

        else:
            if n < 0:
                if fmpq_poly_is_zero(self._poly):
                    raise ZeroDivisionError("negative exponent in power of zero")
                res = self._new()
                sig_str("FLINT exception")
                fmpq_poly_pow(res._poly, self._poly, -n)
                sig_off()
                return ~res
            else:
                res = self._new()
                sig_str("FLINT exception")
                if self._is_gen:
                    fmpq_poly_set_coeff_si(res._poly, n, 1)
                else:
                    fmpq_poly_pow(res._poly, self._poly, n)
                sig_off()
                return res

    def __floordiv__(Polynomial_rational_flint self, right):
        """
        Return the quotient of ``self`` and ``right`` obtained by Euclidean division.

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
        cdef fmpq_t tmpfq

        if right == 0:
            raise ZeroDivisionError("division by zero polynomial")

        if not isinstance(right, Polynomial_rational_flint):
            if right in QQ:
                res = self._new()
                do_sig = _do_sig(self._poly)

                if do_sig: sig_str("FLINT exception")
                fmpq_init_set_readonly(tmpfq, (<Rational> QQ(right)).value)
                fmpq_poly_scalar_div_fmpq(res._poly, self._poly, tmpfq)
                fmpq_clear_readonly(tmpfq)
                if do_sig: sig_off()
                return res

            right = self._parent(right)

        res = self._new()
        sig_str("FLINT exception")
        fmpq_poly_div(res._poly, self._poly,
                                     (<Polynomial_rational_flint>right)._poly)
        sig_off()
        return res

    cpdef Polynomial inverse_series_trunc(self, long prec):
        r"""
        Return a polynomial approximation of precision ``prec`` of the inverse
        series of this polynomial.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: p = 2 + x - 3/5*x**2
            sage: q5 = p.inverse_series_trunc(5)
            sage: q5
            151/800*x^4 - 17/80*x^3 + 11/40*x^2 - 1/4*x + 1/2
            sage: q5 * p
            -453/4000*x^6 + 253/800*x^5 + 1

            sage: q100 = p.inverse_series_trunc(100)
            sage: (q100 * p).truncate(100)
            1

        TESTS::

            sage: (0*x).inverse_series_trunc(4)
            Traceback (most recent call last):
            ...
            ValueError: constant term is zero
            sage: x.inverse_series_trunc(4)
            Traceback (most recent call last):
            ...
            ValueError: constant term is zero
            sage: (x+1).inverse_series_trunc(0)
            Traceback (most recent call last):
            ...
            ValueError: the precision must be positive, got 0
        """
        if prec <= 0:
            raise ValueError("the precision must be positive, got {}".format(prec))
        if fmpq_poly_degree(self._poly) == -1 or \
           fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term is zero")

        cdef Polynomial_rational_flint res = self._new()
        if prec <= 0:
            return res
        sig_on()
        fmpq_poly_inv_series(res._poly, self._poly, prec)
        sig_off()
        return res

    cpdef _mod_(self, right):
        """
        Return the remainder of ``self`` and ``right`` obtain by Euclidean division.

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
            raise ZeroDivisionError("division by zero polynomial")

        res = self._new()
        sig_str("FLINT exception")
        fmpq_poly_rem(res._poly, self._poly,
                                     (<Polynomial_rational_flint>right)._poly)
        sig_off()
        return res

    ###########################################################################
    # Further methods                                                         #
    ###########################################################################

    def numerator(self):
        """
        Return the numerator of ``self``.

        Representing ``self`` as the quotient of an integer polynomial and
        a positive integer denominator (coprime to the content of the
        polynomial), returns the integer polynomial.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = (3 * t^3 + 1) / -3
            sage: f.numerator()
            -3*t^3 - 1
        """
        cdef Polynomial_integer_dense_flint num
        num = Polynomial_integer_dense_flint.__new__(Polynomial_integer_dense_flint)
        parent = ZZ[self.variable_name()]
        Polynomial_integer_dense_flint.__init__(num, parent, x=None,
                                    check=False, is_gen=False, construct=False)
        sig_str("FLINT exception")
        fmpq_poly_get_numerator(num._poly, self._poly)
        sig_off()
        return num

    def denominator(self):
        """
        Return the denominator of ``self``.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f = (3 * t^3 + 1) / -3
            sage: f.denominator()
            3
        """
        cdef Integer den = Integer.__new__(Integer)
        if fmpq_poly_denref(self._poly) is NULL:
            mpz_set_ui(den.value, 1)
        else:
            fmpz_get_mpz(den.value, <fmpz *> fmpq_poly_denref(self._poly))
        return den

    def _derivative(self, var=None):
        """
        Return the derivative of this polynomial with respect to ``var``.

        INPUT:

        - ``var`` -- must be either (equal to) the generator of the polynomial
          ring to which this polynomial belongs, or ``None``; either way the
          behaviour is the same

        OUTPUT: derivative as a ``Polynomial_rational_flint``

        .. SEEALSO:: :meth:`~Polynomial.derivative`

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
            ValueError: cannot differentiate with respect to 2*x

        Check that :issue:`28187` is fixed::

            sage: x = var("x")                                                          # needs sage.symbolic
            sage: f._derivative(x)                                                      # needs sage.symbolic
            4*x^3 - 1
        """
        cdef Polynomial_rational_flint der
        cdef bint do_sig

        if var is not None and var != self._parent.gen():
            raise ValueError("cannot differentiate with respect to {}".format(var))

        der = self._new()
        do_sig = _do_sig(self._poly)

        if do_sig: sig_str("FLINT exception")
        fmpq_poly_derivative(der._poly, self._poly)
        if do_sig: sig_off()
        return der

    def real_root_intervals(self):
        """
        Return isolating intervals for the real roots of ``self``.

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
        Return the resultant of ``self`` and ``right``.

        Enumerating the roots over `\QQ` as `r_1, \dots, r_m` and
        `s_1, \dots, s_n` and letting `x` and `y` denote the leading
        coefficients of `f` and `g`, the resultant of the two polynomials
        is defined by

        .. MATH::

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
        cdef Rational res = Rational.__new__(Rational)
        cdef fmpq_t t
        fmpq_init(t)
        sig_str("FLINT exception")
        fmpq_poly_resultant(t, self._poly,
                            (<Polynomial_rational_flint>right)._poly)
        fmpq_get_mpq(res.value, t)
        sig_off()
        fmpq_clear(t)
        return res

    @cached_method
    def is_irreducible(self):
        r"""
        Return whether this polynomial is irreducible.

        This method computes the primitive part as an element of `\ZZ[t]` and
        calls the method ``is_irreducible`` for elements of that polynomial
        ring.

        By definition, over any integral domain, an element `r` is irreducible
        if and only if it is nonzero, not a unit and whenever `r = ab` then
        `a` or `b` is a unit.

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
            sage: (t + 1).is_irreducible()
            True

        Test that caching works::

           sage: R.<t> = QQ[]
           sage: f = t + 1
           sage: f.is_irreducible()
           True
           sage: f.is_irreducible.cache
           True
        """
        cdef Polynomial_integer_dense_flint primitive
        cdef unsigned long length = fmpq_poly_length(self._poly)

        if length < 2:
            return False
        elif length == 2:
            return True
        else:
            primitive = Polynomial_integer_dense_flint.__new__(Polynomial_integer_dense_flint)
            parent = ZZ[self.variable_name()]
            Polynomial_integer_dense_flint.__init__(primitive, parent,
                             x=None, check=True, is_gen=False, construct=False)

            sig_str("FLINT exception")
            fmpq_poly_get_numerator(primitive._poly, self._poly)
            fmpz_poly_primitive_part(primitive._poly, primitive._poly)
            sig_off()
            return primitive.is_irreducible()

    #######################################################
    # Transcendental functions (return truncated series)  #
    #######################################################

    def _log_series(self, long prec):
        r"""
        Return the logarithm of this polynomial up to precision ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: (1+x)._log_series(5)
            -1/4*x^4 + 1/3*x^3 - 1/2*x^2 + x

            sage: (1/3*x^3 - 2*x^2 + x + 1)._log_series(10)._exp_series(10)
            1/3*x^3 - 2*x^2 + x + 1

        TESTS::

            sage: x._log_series(5)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 1 in order to take logarithm
            sage: (0*x)._log_series(5)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 1 in order to take logarithm
        """
        if fmpq_poly_degree(self._poly) == -1 or \
           fmpz_cmp(fmpq_poly_numref(self._poly),
                    fmpq_poly_denref(self._poly)):
            raise ValueError("constant term should be 1 in order to take logarithm")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_log_series(res._poly, self._poly, prec)
        sig_off()
        return res

    def _exp_series(self, long prec):
        r"""
        Return the exponential of this polynomial up to precision ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: x._exp_series(5)
            1/24*x^4 + 1/6*x^3 + 1/2*x^2 + x + 1
            sage: (1/3*x^4 - 3*x^2 - 1/2*x)._exp_series(5)._log_series(5)
            1/3*x^4 - 3*x^2 - 1/2*x

        TESTS::

            sage: (x+1)._exp_series(5)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 0 in order to take exponential
            sage: (0*x)._exp_series(5)
            1
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if fmpq_poly_degree(self._poly) == -1:
            return self._parent.one()
        elif not fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term should be 0 in order to take exponential")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_exp_series(res._poly, self._poly, prec)
        sig_off()
        return res

    def _atan_series(self, long prec):
        r"""
        Return the arctangent of this polynomial up to precision ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: x._atan_series(7)
            1/5*x^5 - 1/3*x^3 + x
            sage: (1/5*x^3 + 2*x^2 - x)._atan_series(10)._tan_series(10)
            1/5*x^3 + 2*x^2 - x

        TESTS::

            sage: (1+x)._atan_series(3)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 0 in order to take arctangent
            sage: (0*x)._atan_series(10)
            0
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if fmpq_poly_degree(self._poly) == -1:
            return self._parent.zero()
        elif not fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term should be 0 in order to take arctangent")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_atan_series(res._poly, self._poly, prec)
        sig_off()
        return res

    def _atanh_series(self, long prec):
        r"""
        Return the hyperbolic arctangent of this polynomial up to precision
        ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: x._atanh_series(7)
            1/5*x^5 + 1/3*x^3 + x
            sage: (1/5*x^3 + 2*x^2 - x)._atanh_series(10)._tanh_series(10)
            1/5*x^3 + 2*x^2 - x

        TESTS::

            sage: (0*x)._atanh_series(10)
            0
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if fmpq_poly_degree(self._poly) == -1:
            return self._parent.zero()
        elif not fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term should be 0 in order to take hyperbolic arctangent")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_atanh_series(res._poly, self._poly, prec)
        sig_off()
        return res

    def _asin_series(self, long prec):
        r"""
        Return the arcsine of this polynomial up to precision ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: x._asin_series(7)
            3/40*x^5 + 1/6*x^3 + x
            sage: (1/5*x^3 + 2*x^2 - x)._asin_series(10)._sin_series(10)
            1/5*x^3 + 2*x^2 - x

        TESTS::

            sage: (x+1)._asin_series(5)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 0 in order to take arcsine
            sage: (0*x)._asin_series(5)
            0
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if fmpq_poly_degree(self._poly) == -1:
            return self._parent.zero()
        elif not fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term should be 0 in order to take arcsine")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_asin_series(res._poly, self._poly, prec)
        sig_off()
        return res

    def _asinh_series(self, long prec):
        r"""
        Return the hyperbolic arcsine of this polynomial up to precision
        ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: x._asinh_series(7)
            3/40*x^5 - 1/6*x^3 + x
            sage: (1/5*x^3 + 2*x^2 - x)._asinh_series(10)._sinh_series(10)
            1/5*x^3 + 2*x^2 - x

        TESTS::

            sage: (x+1)._asinh_series(5)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 0 in order to take hyperbolic arcsine
            sage: (0*x)._asinh_series(5)
            0
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if fmpq_poly_degree(self._poly) == -1:
            return self._parent.zero()
        elif not fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term should be 0 in order to take hyperbolic arcsine")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_asinh_series(res._poly, self._poly, prec)
        sig_off()
        return res

    def _tan_series(self, long prec):
        r"""
        Return the tangent of this polynomial up to precision ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: x._tan_series(8)
            17/315*x^7 + 2/15*x^5 + 1/3*x^3 + x
            sage: (1/5*x^3 + 2*x^2 - x)._tan_series(10)._atan_series(10)
            1/5*x^3 + 2*x^2 - x

        TESTS::

            sage: (x+1)._tan_series(10)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 0 in order to take tangent
            sage: (0*x)._tan_series(5)
            0
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if fmpq_poly_degree(self._poly) == -1:
            return self._parent.zero()
        elif not fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term should be 0 in order to take tangent")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_tan_series(res._poly, self._poly, prec)
        sig_off()
        return res

    def _sin_series(self, long prec):
        r"""
        Return the sine of this polynomial up to precision ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: x._sin_series(8)
            -1/5040*x^7 + 1/120*x^5 - 1/6*x^3 + x
            sage: (1/5*x^3 - 2*x^2 + 1/2*x)._sin_series(10)._asin_series(10)
            1/5*x^3 - 2*x^2 + 1/2*x

        TESTS::

            sage: (x+1)._sin_series(10)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 0 in order to take sine
            sage: (0*x)._sin_series(5)
            0
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if fmpq_poly_degree(self._poly) == -1:
            return self._parent.zero()
        elif not fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term should be 0 in order to take sine")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_sin_series(res._poly, self._poly, prec)
        sig_off()
        return res

    def _cos_series(self, long prec):
        r"""
        Return the cosine of this polynomial up to precision ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: x._cos_series(10)
            1/40320*x^8 - 1/720*x^6 + 1/24*x^4 - 1/2*x^2 + 1

        TESTS::

            sage: (x+1)._cos_series(5)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 0 in order to take cosine
            sage: (0*x)._cos_series(5)
            0
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if fmpq_poly_degree(self._poly) == -1:
            return self._parent.zero()
        elif not fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term should be 0 in order to take cosine")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_cos_series(res._poly, self._poly, prec)
        sig_off()
        return res

    def _sinh_series(self, long prec):
        r"""
        Return the hyperbolic sine of this polynomial up to precision ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: x._sinh_series(8)
            1/5040*x^7 + 1/120*x^5 + 1/6*x^3 + x

        TESTS::

            sage: (x+1)._sinh_series(5)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 0 in order to take hyperbolic sine
            sage: (0*x)._sinh_series(5)
            0
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if fmpq_poly_degree(self._poly) == -1:
            return self._parent.zero()
        elif not fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term should be 0 in order to take hyperbolic sine")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_sinh_series(res._poly, self._poly, prec)
        sig_off()
        return res

    def _cosh_series(self, long prec):
        r"""
        Return the hyperbolic cosine of this polynomial up to precision
        ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: x._cosh_series(8)
            1/720*x^6 + 1/24*x^4 + 1/2*x^2 + 1

        A trigonometric identity::

            sage: x._cosh_series(8) + x._sinh_series(8) == x._exp_series(8)
            True

        TESTS::

            sage: (x+1)._cosh_series(5)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 0 in order to take hyperbolic cosine
            sage: (0*x)._cosh_series(5)
            1
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if fmpq_poly_degree(self._poly) == -1:
            return self._parent.one()
        elif not fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term should be 0 in order to take hyperbolic cosine")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_cosh_series(res._poly, self._poly, prec)
        sig_off()
        return res

    def _tanh_series(self, long prec):
        r"""
        Return the hyperbolic tangent of this polynomial up to precision
        ``prec``.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: x._tanh_series(8)
            -17/315*x^7 + 2/15*x^5 - 1/3*x^3 + x

        TESTS::

            sage: (x+1)._tanh_series(5)
            Traceback (most recent call last):
            ...
            ValueError: constant term should be 0 in order to take hyperbolic tangent
            sage: (0*x)._tanh_series(5)
            0
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        if fmpq_poly_degree(self._poly) == -1:
            return self._parent.zero()
        elif not fmpz_is_zero(fmpq_poly_numref(self._poly)):
            raise ValueError("constant term should be 0 in order to take hyperbolic tangent")

        cdef Polynomial_rational_flint res = self._new()
        sig_on()
        fmpq_poly_tanh_series(res._poly, self._poly, prec)
        sig_off()
        return res

    ###########################################################################
    # Methods using PARI                                                      #
    ###########################################################################

    def galois_group(self, pari_group=False, algorithm='pari'):
        """
        Return the Galois group of this polynomial as a permutation group.

        INPUT:

        - ``self`` -- irreducible polynomial

        - ``pari_group`` -- boolean (default: ``False``); if ``True`` instead
          return the Galois group as a PARI group.  This has a useful label
          in it, and may be slightly faster since it doesn't require looking
          up a group in GAP.  To get a permutation group from a PARI
          group ``P``, type ``PermutationGroup(P)``.

        - ``algorithm`` -- ``'pari'``, ``'gap'``, ``'kash'``, ``'magma'`` (default:
          ``'pari'``, for degrees is at most 11;
          ``'gap'``, for degrees from 12 to 15;
          ``'kash'``, for degrees from 16 or more).

        OUTPUT: Galois group

        ALGORITHM:

        The Galois group is computed using PARI in C library mode, or possibly
        GAP, KASH, or MAGMA.

        .. NOTE::

            The PARI documentation contains the following warning: The method
            used is that of resolvent polynomials and is sensitive to the
            current precision. The precision is updated internally but, in very
            rare cases, a wrong result may be returned if the initial precision
            was not sufficient.

            GAP uses the "Transitive Groups Libraries" from the "TransGrp"
            GAP package which comes installed with the "gap" Sage package.

            MAGMA does not return a provably correct result.  Please see the
            MAGMA documentation for how to obtain a provably correct result.

        EXAMPLES::

            sage: # needs sage.groups sage.libs.pari
            sage: R.<x> = QQ[]
            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: G = f.galois_group(); G
            Transitive group number 5 of degree 4
            sage: G.gens()
            ((1,2,3,4), (1,2))
            sage: G.order()
            24

        It is potentially useful to instead obtain the corresponding PARI
        group, which is little more than a 4-tuple.  See the PARI manual for
        the exact details.  (Note that the third entry in the tuple is in the
        new standard ordering.)

        ::

            sage: # needs sage.groups sage.libs.pari
            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: G = f.galois_group(pari_group=True); G
            PARI group [24, -1, 5, "S4"] of degree 4
            sage: PermutationGroup(G)
            Transitive group number 5 of degree 4

        You can use KASH or GAP to compute Galois groups as well.  The advantage is
        that KASH (resp. GAP) can compute Galois groups of fields up to
        degree 23 (resp. 15), whereas PARI only goes to degree 11.
        (In my not-so-thorough experiments PARI is faster than KASH.)

        ::

            sage: R.<x> = QQ[]
            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: f.galois_group(algorithm='kash')   # optional - kash
            Transitive group number 5 of degree 4

            sage: # needs sage.libs.gap
            sage: f = x^4 - 17*x^3 - 2*x + 1
            sage: f.galois_group(algorithm='gap')
            Transitive group number 5 of degree 4
            sage: f = x^13 - 17*x^3 - 2*x + 1
            sage: f.galois_group(algorithm='gap')
            Transitive group number 9 of degree 13
            sage: f = x^12 - 2*x^8 - x^7 + 2*x^6 + 4*x^4 - 2*x^3 - x^2 - x + 1
            sage: f.galois_group(algorithm='gap')
            Transitive group number 183 of degree 12

            sage: f.galois_group(algorithm='magma')  # optional - magma
            Transitive group number 5 of degree 4

        TESTS:

        We illustrate the behaviour in the case of reducible polynomials::

            sage: R.<t> = QQ[]
            sage: f = (1 + t)^2
            sage: f.galois_group()
            Traceback (most recent call last):
            ...
            ValueError: The polynomial must be irreducible

        Variable names that are reserved in PARI, such as ``zeta``,
        are supported (see :issue:`20631`)::

            sage: R.<zeta> = QQ[]
            sage: (zeta^2 + zeta + 1).galois_group(pari_group=True)                     # needs sage.libs.pari
            PARI group [2, -1, 1, "S2"] of degree 2
        """
        from sage.groups.pari_group import PariGroup
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.groups.perm_gps.permgroup_named import TransitiveGroup

        if not self.is_irreducible():
            raise ValueError("The polynomial must be irreducible")

        if self.degree() > 11 and algorithm == 'pari':
            if self.degree() < 16:
                algorithm = 'gap'
            else:
                algorithm = 'kash'

        if self.degree() > 21 and algorithm == 'kash':
            raise NotImplementedError("Galois group computation is "
                "supported for degrees up to 11 using PARI, or up to 21 "
                "if KASH is installed.  Try "
                "algorithm='magma' if you have magma.")

        if algorithm == 'pari':
            G = self._pari_with_name().Polrev().polgalois()
            H = PariGroup(G, self.degree())
            if pari_group:
                return H
            else:
                return PermutationGroup(H)

        elif algorithm == 'kash':
            try:
                from sage.interfaces.kash import kash
                kash.eval('X := PolynomialRing(RationalField()).1')
                s = self._repr(name='X')
                G = kash('Galois(%s)'%s)
                d = int(kash.eval('%s.ext1'%G.name()))
                n = int(kash.eval('%s.ext2'%G.name()))
                return TransitiveGroup(d, n)
            except RuntimeError as msg:
                raise NotImplementedError(str(msg) + "\nSorry, " +
                    "computation of Galois groups of fields of degree " +
                    "bigger than 11 is not yet implemented.  Try installing " +
                    "the optional free (closed source) KASH software, which " +
                    "supports degrees up to 21, or use algorithm='magma' if " +
                    "you have magma.")

        elif algorithm == 'gap':
            if self.degree() > 15:
                raise NotImplementedError("Galois group computation is " +
                    "supported for degrees up to 15 using GAP. Try " +
                    "algorithm='kash'.")
            from sage.libs.gap.libgap import libgap
            fgap = libgap(self)
            return TransitiveGroup(self.degree(), fgap.GaloisType())

        elif algorithm == 'magma':
            from sage.interfaces.magma import magma
            X = magma(self).GaloisGroup()
            try:
                n, d = X.TransitiveGroupIdentification(nvals=2)
                d = int(d)
                n = int(n)
            except RuntimeError as msg:
                raise RuntimeError(str(msg) + "\nUnable to lookup " +
                    "description of Galois group as a transitive " +
                    "group.\n%s" % X)
            return TransitiveGroup(d, n)

        else:
            raise ValueError("Algorithm %s not supported." % algorithm)

    def factor_mod(self, p):
        """
        Return the factorization of ``self`` modulo the prime `p`.

        Assumes that the degree of this polynomial is at least one, and raises
        a :exc:`ValueError` otherwise.

        INPUT:

        - ``p`` -- prime number

        OUTPUT: factorization of this polynomial  modulo `p`

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: (x^5 + 17*x^3 + x + 3).factor_mod(3)
            x * (x^2 + 1)^2
            sage: (x^5 + 2).factor_mod(5)
            (x + 2)^5

        Variable names that are reserved in PARI, such as ``zeta``,
        are supported (see :issue:`20631`)::

            sage: R.<zeta> = QQ[]
            sage: (zeta^2 + zeta + 1).factor_mod(7)
            (zeta + 3) * (zeta + 5)
        """
        from sage.rings.finite_rings.finite_field_constructor import FiniteField

        p = Integer(p)
        if not p.is_prime():
            raise ValueError("p must be prime")

        if self.degree() < 1:
            raise ValueError("The polynomial must have degree at least 1")

        G = self._pari_with_name().factormod(p)
        K = FiniteField(p)
        R = K[self.parent().variable_name()]
        return R(1)._factor_pari_helper(G, unit=R(self).leading_coefficient())

    def factor_padic(self, p, prec=10):
        r"""
        Return the `p`-adic factorization of this polynomial to the given
        precision.

        INPUT:

        - ``p`` -- prime number

        - ``prec`` -- integer; the precision

        OUTPUT: factorization of ``self`` viewed as a `p`-adic polynomial

        EXAMPLES::

            sage: # needs sage.rings.padic
            sage: R.<x> = QQ[]
            sage: f = x^3 - 2
            sage: f.factor_padic(2)
            (1 + O(2^10))*x^3 + O(2^10)*x^2 + O(2^10)*x
            + 2 + 2^2 + 2^3 + 2^4 + 2^5 + 2^6 + 2^7 + 2^8 + 2^9 + O(2^10)
            sage: f.factor_padic(3)
            (1 + O(3^10))*x^3 + O(3^10)*x^2 + O(3^10)*x
            + 1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10)
            sage: f.factor_padic(5)
            ((1 + O(5^10))*x
             + 2 + 4*5 + 2*5^2 + 2*5^3 + 5^4 + 3*5^5 + 4*5^7 + 2*5^8 + 5^9 + O(5^10))
            * ((1 + O(5^10))*x^2
               + (3 + 2*5^2 + 2*5^3 + 3*5^4 + 5^5 + 4*5^6 + 2*5^8 + 3*5^9 + O(5^10))*x
               + 4 + 5 + 2*5^2 + 4*5^3 + 4*5^4 + 3*5^5 + 3*5^6 + 4*5^7 + 4*5^9 + O(5^10))

        The input polynomial is considered to have "infinite" precision,
        therefore the `p`-adic factorization of the polynomial is not
        the same as first coercing to `\QQ_p` and then factoring
        (see also :issue:`15422`)::

            sage: # needs sage.rings.padic
            sage: f = x^2 - 3^6
            sage: f.factor_padic(3, 5)
            ((1 + O(3^5))*x + 3^3 + O(3^5)) * ((1 + O(3^5))*x + 2*3^3 + 2*3^4 + O(3^5))
            sage: f.change_ring(Qp(3,5)).factor()
            Traceback (most recent call last):
            ...
            PrecisionError: p-adic factorization not well-defined since
            the discriminant is zero up to the requestion p-adic precision

        A more difficult example::

            sage: R.<x> = QQ[]
            sage: f = 100 * (5*x + 1)^2 * (x + 5)^2
            sage: f.factor_padic(5, 10)                                                 # needs sage.rings.padic
            (4*5^4 + O(5^14)) * ((1 + O(5^9))*x + 5^-1 + O(5^9))^2
            * ((1 + O(5^10))*x + 5 + O(5^10))^2

        Try some bogus inputs::

            sage: # needs sage.rings.padic
            sage: f.factor_padic(3, -1)
            Traceback (most recent call last):
            ...
            ValueError: prec_cap must be nonnegative
            sage: f.factor_padic(6, 10)
            Traceback (most recent call last):
            ...
            ValueError: p must be prime
            sage: f.factor_padic('hello', 'world')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'hello' to an integer
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

        Return an empty list if this polynomial has degree less than one.

        INPUT:

        - ``p`` -- prime number; coerceable to :class:`Integer`
        - ``e`` -- exponent; coerceable to :class:`Integer`

        OUTPUT: Hensel lifts; list of polynomials over `\ZZ / p^e \ZZ`

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
            sage: R(x - 1).hensel_lift(7, 2)
            [x + 48]

        Variable names that are reserved in PARI, such as ``I``, are
        supported (see :issue:`20631`)::

            sage: R.<I> = QQ[]
            sage: (I^2 + 1).hensel_lift(5, 3)
            [I + 57, I + 68]
            sage: (I^2 + 1).hensel_lift(2, 3)
            Traceback (most recent call last):
            ...
            ValueError: I^2 + 1 is not square-free modulo 2
        """
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

        p = Integer(p)
        if not p.is_prime():
            raise ValueError("p must be prime")
        e = Integer(e)
        if e < 1:
            raise ValueError("e must be at least 1")

        # The relevant PARI method doesn't seem to play well with constant and
        # linear polynomials, so we handle these separately.
        #
        if self.degree() < 1:
            return [ ]
        elif self.degree() == 1:
            R = IntegerModRing(p**e)
            S = R[self.parent().variable_name()]
            return [S(self)]

        f = self._pari_with_name()
        F = f.factormod(p)
        if any(n > 1 for n in F[1]):
            raise ValueError("{} is not square-free modulo {}".format(self, p))
        H = f.polhensellift(F[0].liftint(), p, e)
        R = IntegerModRing(p**e)
        S = R[self.parent().variable_name()]
        return [S(m) for m in H]

    def discriminant(self):
        r"""
        Return the discriminant of this polynomial.

        The discriminant `R_n` is defined as

        .. MATH::

            R_n = a_n^{2 n-2} \prod_{1 \le i < j \le n} (r_i - r_j)^2,

        where `n` is the degree of this polynomial, `a_n` is the leading
        coefficient and the roots over `\QQbar` are `r_1, \ldots, r_n`.

        The discriminant of constant polynomials is defined to be 0.

        OUTPUT: discriminant, an element of the base ring of the polynomial ring

        .. NOTE::

            Note the identity `R_n(f) := (-1)^{(n (n-1)/2)} R(f,f') a_n^{(n-k-2)}`,
            where `n` is the degree of this polynomial, `a_n` is the leading
            coefficient, `f'` is the derivative of `f`, and `k` is the degree
            of `f'`.  Calls :meth:`resultant`.

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
            sage: EllipticCurve([1, 1]).discriminant() / 16                             # needs sage.schemes
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

        Variable names that are reserved in PARI, such as ``I``, are
        supported (see :issue:`20631`)::

            sage: R.<I> = QQ[]
            sage: (I^2 + 1).discriminant()
            -4
        """
        return QQ(self._pari_with_name().poldisc())

    # Alias for discriminant
    disc = discriminant

    def galois_group_davenport_smith_test(self, num_trials=50, assume_irreducible=False):
        """
        Use the Davenport-Smith test to attempt to certify that `f` has Galois group `A_n` or `S_n`.

        Return 1 if the Galois group is certified as `S_n`, 2 if `A_n`, or 0 if no conclusion is reached.

        By default, we first check that `f` is irreducible. For extra efficiency, one can override this
        by specifying ``assume_irreducible=True``; this yields undefined results if `f` is not irreducible.

        A corresponding function in Magma is ``IsEasySnAn``.

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: u = x^7 + x + 1
            sage: u.galois_group_davenport_smith_test()
            1
            sage: u = x^7 - x^4 - x^3 + 3*x^2 - 1
            sage: u.galois_group_davenport_smith_test()
            2
            sage: u = x^7 - 2
            sage: u.galois_group_davenport_smith_test()
            0
        """
        from sage.arith.misc import primes_first_n
        from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

        if not assume_irreducible and not self.is_irreducible():
            return 0
        d = self.degree()
        for p in primes_first_n(num_trials):
            fp = self.change_ring(IntegerModRing(p))
            g = fp.factor()[-1][0]
            d1 = g.degree()
            # Here we use the fact that a transitive permutation representation with a long prime cycle
            # must have image at least as big as A_n.
            if (d1 <= 7 and (d,d1) in ((1,1),(2,2),(3,2),(3,3),(4,3),(5,3),(5,4),(6,5),(7,5))) or\
                (d1 > d/2 and d1 < d-2 and d1.is_prime()):
                return (2 if self.disc().is_square() else 1)
        return 0
