# distutils: extra_compile_args = -D_XPG6
"""
FLINT fmpz_poly class wrapper

AUTHORS:

- Robert Bradshaw (2007-09-15) Initial version.
- William Stein (2007-10-02) update for new flint; add arithmetic and creation
  of coefficients of arbitrary size.
"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.sequence cimport *

from cysignals.memory cimport sig_free

from sage.arith.long cimport pyobject_to_long
from sage.cpython.string cimport char_to_str, str_to_bytes
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.structure.sage_object cimport SageObject
from sage.rings.integer cimport Integer

cdef class Fmpz_poly(SageObject):

    def __cinit__(self):
        fmpz_poly_init(self.poly)

    def __init__(self, v):
        """
        Construct a new fmpz_poly from a sequence, constant coefficient,
        or string (in the same format as it prints).

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: Fmpz_poly([1,2,3])
            3  1 2 3
            sage: Fmpz_poly(5)
            1  5
            sage: Fmpz_poly(str(Fmpz_poly([3,5,7])))
            3  3 5 7
        """
        cdef Py_ssize_t i
        cdef long c
        cdef Integer w
        if isinstance(v, str):
            if not fmpz_poly_set_str(self.poly, str_to_bytes(v)):
                return
            else:
                raise ValueError("Unable to create Fmpz_poly from that string.")
        if not PySequence_Check(v):
            v = [v]
        try:
            fmpz_poly_set_coeff_si(self.poly, 0, 1)
            fmpz_poly_set_coeff_si(self.poly, 0, 0)
            for i from 0 <= i < len(v):
                #fmpz_poly_set_coeff_si(self.poly, i, v[i])
                w = Integer(v[i])
                fmpz_poly_set_coeff_mpz(self.poly, i, w.value)
        except OverflowError:
            raise ValueError("No fmpz_poly_set_coeff_mpz() method.")

    def __dealloc__(self):
        fmpz_poly_clear(self.poly)

    def __setitem__(self, i, value):
        """
        Set the `i`-th item of self, which is the coefficient of the `x^i` term.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly(range(10))
            sage: f[7] = 100; f
            10  0 1 2 3 4 5 6 100 8 9
            sage: f[2] = 10**100000
            sage: f[2] == 10**100000
            True
        """
        if isinstance(value, Integer):
            fmpz_poly_set_coeff_mpz(self.poly, i, (<Integer>value).value)
        else:
            fmpz_poly_set_coeff_si(self.poly, i, value)

    def __getitem__(self, i):
        """
        Return the `i`-th item of self, which is the coefficient of the `x^i` term.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly(range(100))
            sage: f[13]
            13
            sage: f[200]
            0
        """
        cdef Integer res = Integer.__new__(Integer)
        fmpz_poly_get_coeff_mpz(res.value, self.poly, i)
        return res

    def __repr__(self):
        """
        Print ``self`` according to the native FLINT format.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([0,1]); f^7
            8  0 0 0 0 0 0 0 1
        """
        cdef char* ss = fmpz_poly_get_str(self.poly)
        cdef object s = char_to_str(ss)
        sig_free(ss)
        return s

    def degree(self):
        """
        The degree of ``self``.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([1,2,3]); f
            3  1 2 3
            sage: f.degree()
            2
            sage: Fmpz_poly(range(1000)).degree()
            999
            sage: Fmpz_poly([2,0]).degree()
            0
        """
        return fmpz_poly_degree(self.poly)

    def list(self):
        """
        Return ``self`` as a list of coefficients, lowest terms first.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([2,1,0,-1])
            sage: f.list()
            [2, 1, 0, -1]
        """
        return [self[i] for i in range(self.degree() + 1)]

    def __add__(left, right):
        """
        Add together two Flint polynomials.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: Fmpz_poly([1,2,3]) + Fmpz_poly(range(6))
            6  1 3 5 3 4 5
        """
        if not isinstance(left, Fmpz_poly) or not isinstance(right, Fmpz_poly):
            raise TypeError
        cdef Fmpz_poly res = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        fmpz_poly_add(res.poly, (<Fmpz_poly>left).poly, (<Fmpz_poly>right).poly)
        return res

    def __sub__(left, right):
        """
        Subtract two Flint polynomials.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: Fmpz_poly([10,2,3]) - Fmpz_poly([4,-2,1])
            3  6 4 2
        """
        if not isinstance(left, Fmpz_poly) or not isinstance(right, Fmpz_poly):
            raise TypeError
        cdef Fmpz_poly res = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        fmpz_poly_sub(res.poly, (<Fmpz_poly>left).poly, (<Fmpz_poly>right).poly)
        return res

    def __neg__(self):
        """
        Return the negative of ``self``.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: -Fmpz_poly([2,10,2,3,18,-5])
            6  -2 -10 -2 -3 -18 5
        """
        cdef Fmpz_poly res = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        fmpz_poly_neg(res.poly, self.poly)
        return res

    def __mul__(left, right):
        """
        Return the product of left and right.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([0,1]); g = Fmpz_poly([2,3,4])
            sage: f*g
            4  0 2 3 4
            sage: f = Fmpz_poly([1,0,-1]); g = Fmpz_poly([2,3,4])
            sage: f*g
            5  2 3 2 -3 -4

            Scalar multiplication
            sage: f * 3
            3  3 0 -3
            sage: f * 5r
            3  5 0 -5
        """
        cdef Fmpz_poly res = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        if not isinstance(left, Fmpz_poly) or not isinstance(right, Fmpz_poly):
            if isinstance(left, int):
                fmpz_poly_scalar_mul_si(res.poly, (<Fmpz_poly>right).poly, left)
            elif isinstance(left, Integer):
                fmpz_poly_scalar_mul_mpz(res.poly, (<Fmpz_poly>right).poly, (<Integer>left).value)
            elif isinstance(right, int):
                fmpz_poly_scalar_mul_si(res.poly, (<Fmpz_poly>left).poly, right)
            elif isinstance(right, Integer):
                fmpz_poly_scalar_mul_mpz(res.poly, (<Fmpz_poly>left).poly, (<Integer>right).value)
            else:
                raise TypeError
        else:
            fmpz_poly_mul(res.poly, (<Fmpz_poly>left).poly, (<Fmpz_poly>right).poly)
        return res

    def __pow__(self, n, dummy):
        """
        Return ``self`` raised to the power of `n`.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([1,1])
            sage: f**6
            7  1 6 15 20 15 6 1
            sage: f = Fmpz_poly([2])
            sage: f^150
            1  1427247692705959881058285969449495136382746624
            sage: 2^150
            1427247692705959881058285969449495136382746624

            sage: f**(3/2)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert rational 3/2 to an integer
        """
        cdef long nn = pyobject_to_long(n)
        if not isinstance(self, Fmpz_poly):
            raise TypeError
        cdef Fmpz_poly res = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        fmpz_poly_pow(res.poly, (<Fmpz_poly>self).poly, nn)
        return res

    def pow_truncate(self, exp, n):
        """
        Return ``self`` raised to the power of ``exp`` mod `x^n`.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([1,2])
            sage: f.pow_truncate(10,3)
            3  1 20 180
            sage: f.pow_truncate(1000,3)
            3  1 2000 1998000
        """
        if exp < 0:
            raise ValueError("Exponent must be at least 0")
        if n < 0:
            raise ValueError("Exponent must be at least 0")
        cdef long exp_c = exp, nn = n
        cdef Fmpz_poly res = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        fmpz_poly_pow_trunc(res.poly, (<Fmpz_poly>self).poly, exp_c, nn)
        return res

    def __floordiv__(left, right):
        """
        Return left // right, truncated.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([3,4,5])
            sage: g = f^5; g
            11  243 1620 6345 16560 32190 47224 53650 46000 29375 12500 3125
            sage: g // f
            9  81 432 1404 2928 4486 4880 3900 2000 625
            sage: f^4
            9  81 432 1404 2928 4486 4880 3900 2000 625
        """
        if not isinstance(left, Fmpz_poly) or not isinstance(right, Fmpz_poly):
            raise TypeError
        cdef Fmpz_poly res = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        fmpz_poly_div(res.poly, (<Fmpz_poly>left).poly, (<Fmpz_poly>right).poly)
        return res

    def div_rem(self, Fmpz_poly other):
        """
        Return ``self / other, self % other``.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([1,3,4,5])
            sage: g = f^23
            sage: g.div_rem(f)[1]
            0
            sage: g.div_rem(f)[0] - f^22
            0
            sage: f = Fmpz_poly([1..10])
            sage: g = Fmpz_poly([1,3,5])
            sage: q, r = f.div_rem(g)
            sage: q*f+r
            17  1 2 3 4 4 4 10 11 17 18 22 26 30 23 26 18 20
            sage: g
            3  1 3 5
            sage: q*g+r
            10  1 2 3 4 5 6 7 8 9 10
        """
        cdef Fmpz_poly Q = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        cdef Fmpz_poly R = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        fmpz_poly_divrem(Q.poly, R.poly, self.poly, other.poly)
        return Q, R

    def left_shift(self, unsigned long n):
        """
        Left shift ``self`` by `n`.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([1,2])
            sage: f.left_shift(1).list() == [0,1,2]
            True
        """
        cdef Fmpz_poly res = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)

        fmpz_poly_shift_left(res.poly, self.poly, n)

        return res

    def right_shift(self, unsigned long n):
        """
        Right shift ``self`` by `n`.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([1,2])
            sage: f.right_shift(1).list() == [2]
            True
        """
        cdef Fmpz_poly res = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)

        fmpz_poly_shift_right(res.poly, self.poly, n)

        return res

    def pseudo_div(self, Fmpz_poly other):
        cdef ulong d
        cdef Fmpz_poly Q = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        fmpz_poly_pseudo_div(Q.poly, &d, self.poly, other.poly)
        return Q, d

    def pseudo_div_rem(self, Fmpz_poly other):
        cdef ulong d
        cdef Fmpz_poly Q = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        cdef Fmpz_poly R = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        fmpz_poly_pseudo_divrem(Q.poly, R.poly, &d, self.poly, other.poly)
        return Q, R, d

    def derivative(self):
        """
        Return the derivative of ``self``.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([1,2,6])
            sage: f.derivative().list() == [2, 12]
            True
        """
        cdef Fmpz_poly res = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)

        fmpz_poly_derivative(res.poly, self.poly)

        return res

    def __copy__(self):
        cdef Fmpz_poly res = <Fmpz_poly>Fmpz_poly.__new__(Fmpz_poly)
        fmpz_poly_set(res.poly, self.poly)
        return res

    def truncate(self, n):
        """
        Return the truncation of ``self`` at degree `n`.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([1,1])
            sage: g = f**10; g
            11  1 10 45 120 210 252 210 120 45 10 1
            sage: g.truncate(5)
            5  1 10 45 120 210
        """
        cdef Fmpz_poly g = self.__copy__()
        fmpz_poly_truncate(g.poly, n)
        return g

    def _unsafe_mutate_truncate(self, n):
        """
        Return the truncation of ``self`` at degree `n`.

        Don't do this unless you know there are no other references to
        this polynomial!!!!!

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([1,1])
            sage: g = f**10; g
            11  1 10 45 120 210 252 210 120 45 10 1
            sage: g._unsafe_mutate_truncate(5); g
            5  1 10 45 120 210
        """
        cdef long nn = n
        fmpz_poly_truncate(self.poly, nn) # mutating!

    def _sage_(self, var='x'):
        """
        Return ``self`` as an element of the sage ``ZZ[var]``.

        EXAMPLES::

            sage: from sage.libs.flint.fmpz_poly_sage import Fmpz_poly
            sage: f = Fmpz_poly([1,1])
            sage: f._sage_('t')
            t + 1
            sage: Fmpz_poly([-1,0,0,1])._sage_()
            x^3 - 1
        """
        from sage.rings.integer_ring import ZZ
        return ZZ[var](self.list())


# Functions removed from flint but still needed in Sage. Code adapted from
# earlier versions of flint.

cdef void fmpz_poly_scalar_mul_mpz(fmpz_poly_t rop, const fmpz_poly_t op, const mpz_t c) noexcept:
    cdef fmpz_t f
    fmpz_init_set_readonly(f, c)
    fmpz_poly_scalar_mul_fmpz(rop, op, f)
    fmpz_clear_readonly(f)

cdef void fmpz_poly_scalar_divexact_mpz(fmpz_poly_t rop, const fmpz_poly_t op, const mpz_t c) noexcept:
    cdef fmpz_t f
    fmpz_init_set_readonly(f, c)
    fmpz_poly_scalar_divexact_fmpz(rop, op, f)
    fmpz_clear_readonly(f)

cdef void fmpz_poly_scalar_fdiv_mpz(fmpz_poly_t rop, const fmpz_poly_t op, const mpz_t c) noexcept:
    cdef fmpz_t f
    fmpz_init_set_readonly(f, c)
    fmpz_poly_scalar_fdiv_fmpz(rop, op, f)
    fmpz_clear_readonly(f)

cdef void fmpz_poly_set_coeff_mpz(fmpz_poly_t poly, slong n, const mpz_t x) noexcept:
    cdef fmpz_t t
    fmpz_init_set_readonly(t, x)
    fmpz_poly_set_coeff_fmpz(poly, n, t)
    fmpz_clear_readonly(t)

cdef void fmpz_poly_get_coeff_mpz(mpz_t x, const fmpz_poly_t poly, slong n) noexcept:
    cdef fmpz_t t
    fmpz_init(t)
    fmpz_poly_get_coeff_fmpz(t, poly, n)
    fmpz_get_mpz(x, t)
    fmpz_clear(t)

cdef void fmpz_poly_set_mpz(fmpz_poly_t poly, const mpz_t x) noexcept:
    fmpz_poly_fit_length(poly, 1)
    fmpz_set_mpz(poly.coeffs, x)
    _fmpz_poly_set_length(poly, 1)
    _fmpz_poly_normalise(poly)
