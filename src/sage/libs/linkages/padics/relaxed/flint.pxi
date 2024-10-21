# sage_setup: distribution = sagemath-flint
r"""
This linkage file implements the relaxed padics API using flint.

AUTHOR:

- Xavier Caruso (2021-01): initial version
"""

# ****************************************************************************
#       Copyright (C) 2021 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.libs.flint.types cimport flint_rand_t
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.flint.fmpz_poly_sage cimport *

cdef extern from "flint_helper.c":
    cdef void flint_randseed(flint_rand_t state, ulong seed1, ulong seed2)
    cdef fmpz* get_coeff(fmpz_poly_t poly, slong i)
    cdef void get_slice(fmpz_poly_t slice, fmpz_poly_t poly, slong start, slong length)
    cdef void iadd_coeff(fmpz_poly_t poly, fmpz_t summand, slong i)
    cdef void isub_coeff(fmpz_poly_t poly, fmpz_t summand, slong i)
    cdef void iadd_shifted(fmpz_poly_t poly, fmpz_poly_t summand, slong shift)
    cdef void reduce_coeff(fmpz_poly_t poly, slong i, fmpz_t modulus)
    cdef void reducesmall_coeff(fmpz_poly_t poly, slong i, fmpz_t modulus)
    cdef void reduceneg_coeff(fmpz_poly_t poly, slong i, fmpz_t modulus)

from sage.rings.padics.pow_computer_flint cimport PowComputer_flint

from sage.ext.stdsage cimport PY_NEW

from sage.rings.finite_rings.finite_field_constructor import GF


# Operations on digits (intended to be small elements in the exact subring)
###########################################################################

cdef fmpz_t digit_zero
digit_init(digit_zero)

cdef inline void digit_init(fmpz_t a) noexcept:
    r"""
    Initialize a digit and set to it the value `0`.

    INPUT:

    - ``a`` -- the ``cdigit`` to initialize
    """
    fmpz_init(a)

cdef inline void digit_clear(fmpz_t a) noexcept:
    r"""
    Deallocate memory assigned to a digit.

    INPUT:

    - ``a`` -- the ``cdigit`` to deallocate
    """
    fmpz_clear(a)

# get and set

cdef inline Integer digit_get_sage(fmpz_t a):
    r"""
    Convert a digit to a Sage element.

    INPUT:

    - ``a`` -- a ``cdigit``

    OUTPUT:

    A Sage element representing the same digit.
    """
    cdef Integer elt = PY_NEW(Integer)
    fmpz_get_mpz(elt.value, a)
    return elt

cdef inline void digit_set(fmpz_t a, fmpz_t b) noexcept:
    r"""
    Set up a digit.

    INPUT:

    - ``a`` -- the ``cdigit`` to set up
    - ``b`` -- the ``cdigit`` containing the value to assign
    """
    fmpz_set(a, b)

cdef inline void digit_set_ui(fmpz_t a, slong b) noexcept:
    r"""
    Set an integral value of a digit.

    INPUT:

    - ``a`` -- the ``cdigit`` to set up
    - ``b`` -- integer; the value of assign
    """
    fmpz_set_ui(a, b)

cdef inline void digit_set_sage(fmpz_t a, Integer elt) noexcept:
    r"""
    Set the value of a digit.

    INPUT:

    - ``a`` -- the ``cdigit`` to be assigned
    - ``b`` -- a Sage element containing the value to assign
    """
    fmpz_set_mpz(a, elt.value)

# comparisons

cdef inline bint digit_equal(fmpz_t a, fmpz_t b) noexcept:
    r"""
    Comparison of two digits.

    INPUT:

    - ``a`` -- a ``cdigit``
    - ``b`` -- a ``cdigit``

    OUTPUT:

    `1` of `a` is equal to `b`, `0` otherwise
    """
    return fmpz_equal(a, b)

cdef inline bint digit_equal_ui(fmpz_t a, slong b) noexcept:
    r"""
    Comparison of a digit and an integer

    INPUT:

    - ``a`` -- a ``cdigit``
    - ``b`` -- integer

    OUTPUT:

    `1` of `a` is equal to `b`, `0` otherwise
    """
    return fmpz_equal_ui(a, b)

cdef inline bint digit_is_zero(fmpz_t a) noexcept:
    r"""
    Comparison to zero

    INPUT:

    - ``a`` -- a ``cdigit``

    OUTPUT:

    `1` of `a` vanishes, `0` otherwise
    """
    return fmpz_is_zero(a)

# random

cdef inline void digit_random_init(flint_rand_t generator, slong seed) noexcept:
    r"""
    Initialize the random generator with a new seed

    INPUT:

    - ``generator`` -- the generator to initialize
    - ``seed`` -- the seed
    """
    flint_randseed(generator, seed, seed*seed + 1)

cdef inline void digit_random(fmpz_t res, PowComputer_flint prime_pow, flint_rand_t generator) noexcept:
    r"""
    Set a digit to a random value in the distinguished set of representatives.

    INPUT:

    - ``res`` -- the ``cdigit`` to be assigned
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_randm(res, generator, prime_pow.fprime)

# operations

cdef inline void digit_add(fmpz_t res, fmpz_t a, fmpz_t b) noexcept:
    r"""
    Add two digits.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the first summand
    - ``b`` -- a ``cdigit``, the second summand
    """
    fmpz_add(res, a, b)

cdef inline void digit_sub(fmpz_t res, fmpz_t a, fmpz_t b) noexcept:
    r"""
    Subtract two digits.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the minuend
    - ``b`` -- a ``cdigit``, the subtrahend
    """
    fmpz_sub(res, a, b)

cdef inline void digit_mul(fmpz_t res, fmpz_t a, fmpz_t b) noexcept:
    r"""
    Multiply two digits.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the first factor
    - ``b`` -- a ``cdigit``, the second factor
    """
    fmpz_mul(res, a, b)

cdef inline void digit_mod(fmpz_t res, fmpz_t a, PowComputer_flint prime_pow) noexcept:
    r"""
    Reduce a digit modulo the uniformizer.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the digit to reduce
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_mod(res, a, prime_pow.fprime)

cdef inline void digit_quorem(fmpz_t quo, fmpz_t rem, fmpz_t a, PowComputer_flint prime_pow) noexcept:
    r"""
    Reduce a digit modulo the uniformizer and keep the carry.

    INPUT:

    - ``quo`` -- a ``cdigit`` to store the carry
    - ``rem`` -- a ``cdigit`` to store the reduction
    - ``a`` -- a ``cdigit``, the digit to reduce
    - ``prime_pow`` -- the PowComputer for the ring
    """
    fmpz_tdiv_qr(quo, rem, a, prime_pow.fprime)

cdef inline void digit_smallest(cdigit res, cdigit carry, cdigit a, PowComputer_flint prime_pow) noexcept:
    r"""
    Compute the smallest representative of a digit.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the smallest representative
    - ``carry`` -- a ``cdigit`` to store the carry
    - ``a`` -- a ``cdigit``, the digit to reduce
    - ``prime_pow`` -- the PowComputer for the ring

    .. NOTE::

        This function assumes that ``a`` is always reduced in the
        usual sense, that is belongs to the range `[0, p-1]`.
    """
    cdef fmpz_t b
    fmpz_init(b)
    fmpz_mul_ui(b, a, 2)
    if fmpz_cmp(b, prime_pow.fprime) > 0:
        fmpz_sub(res, a, prime_pow.fprime)
        fmpz_set_ui(carry, 1)
    else:
        fmpz_set(res, a)
        fmpz_set_ui(carry, 0)
    fmpz_clear(b)

cdef inline void digit_inv(fmpz_t res, fmpz_t a, PowComputer_flint prime_pow) noexcept:
    r"""
    Compute the multiplicative inverse of a digit modulo the uniformizer.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the digit to invert
    - ``prime_pow`` -- the PowComputer for the ring
    """
    cdef fmpz_t gcd
    fmpz_init(gcd)
    fmpz_gcdinv(gcd, res, a, prime_pow.fprime)
    fmpz_clear(gcd)

cdef bint digit_sqrt(fmpz_t ans, fmpz_t x, PowComputer_flint prime_pow) noexcept:
    r"""
    Compute the square root of a digit modulo the uniformizer.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the digit to square root
    - ``prime_pow`` -- the PowComputer for the ring
    """
    return not fmpz_sqrtmod(ans, x, prime_pow.fprime)


# Operations on elements (represented as series of digits)
##########################################################

cdef inline void element_init(fmpz_poly_t x) noexcept:
    r"""
    Initialize an element.

    INPUT:

    - ``x`` -- the ``celement`` to initialize
    """
    fmpz_poly_init(x)

cdef inline void element_clear(fmpz_poly_t x) noexcept:
    r"""
    Deallocate memory assigned to an element.

    INPUT:

    - ``x`` -- the ``celement`` to deallocate
    """
    fmpz_poly_clear(x)

# get and set

cdef inline Integer element_get_sage(fmpz_poly_t x, PowComputer_flint prime_pow):
    r"""
    Convert a digit to a Sage element.

    INPUT:

    - ``x`` -- a ``celement``
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT:

    A Sage integer representing the same element modulo the known precision
    """
    cdef fmpz_t value
    fmpz_init(value)
    fmpz_poly_evaluate_fmpz(value, x, prime_pow.fprime)
    cdef Integer ans = digit_get_sage(value)
    fmpz_clear(value)
    return ans

cdef inline void element_set(fmpz_poly_t x, fmpz_poly_t y) noexcept:
    r"""
    Set an element

    INPUT:

    - ``x`` -- the ``celement`` to be assigned
    - ``y`` -- the ``celement`` containing the value to assign
    """
    fmpz_poly_set(x, y)

# get and set digits

cdef inline fmpz* element_get_digit(fmpz_poly_t x, slong i) noexcept:
    r"""
    Return the `i`-th coefficient of `x`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``i`` -- integer
    """
    return get_coeff(x, i)

cdef inline void element_get_slice(fmpz_poly_t res, fmpz_poly_t x, slong start, slong length) noexcept:
    r"""
    Select a slice of an element.

    INPUT:

    - ``res`` -- a ``celement`` to store the slice
    - ``x`` -- a ``celement``, the element from which the slice is extracted
    - ``start`` -- integer; the start position of the slice
    - ``length`` -- integer; the length of the slice

    .. NOTE::

        The function only sets up a pointer to the requested slice
        (the slice is not copied). Hence any future modification
        of the slice ``res`` will affect the container ``x``.
    """
    get_slice(res, x, start, length)

cdef inline void element_set_digit(fmpz_poly_t x, fmpz_t a, slong i) noexcept:
    r"""
    Set `i`-th coefficient of `x` to the value `a`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``a`` -- a ``cdigit``
    - ``i`` -- integer
    """
    fmpz_poly_set_coeff_fmpz(x, i, a)

cdef inline void element_set_digit_ui(fmpz_poly_t x, slong a, slong i) noexcept:
    r"""
    Set `i`-th coefficient of `x` to the value `a`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``a`` -- integer
    - ``i`` -- integer
    """
    fmpz_poly_set_coeff_ui(x, i, a)

cdef inline void element_set_digit_sage(fmpz_poly_t x, Integer a, slong i) noexcept:
    r"""
    Set `i`-th coefficient of `x` to the value `a`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``a`` -- a Sage element
    - ``i`` -- integer
    """
    fmpz_poly_set_coeff_mpz(x, i, a.value)

# operations

cdef inline void element_iadd_digit(fmpz_poly_t x, fmpz_t a, slong i) noexcept:
    r"""
    Inplace addition:
    add `a` to the `i`-th coefficient of `x`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``a`` -- a ``cdigit``
    - ``i`` -- integer
    """
    iadd_coeff(x, a, i)

cdef inline void element_isub_digit(fmpz_poly_t x, fmpz_t a, slong i) noexcept:
    r"""
    Inplace subtraction:
    subtract `a` to the `i`-th coefficient of `x`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``a`` -- a ``cdigit``
    - ``i`` -- integer
    """
    isub_coeff(x, a, i)

cdef inline void element_iadd_slice(fmpz_poly_t x, fmpz_poly_t slice, slong start) noexcept:
    r"""
    Inplace addition:
    add a slice to an element

    INPUT:

    - ``x`` -- a ``celement``, the element to update
    - ``slice`` -- a ``celement``, the slice to be added
    - ``start`` -- the position from which the slice will be added
    """
    iadd_shifted(x, slice, start)

cdef inline void element_scalarmul(fmpz_poly_t res, fmpz_poly_t x, fmpz_t a) noexcept:
    r"""
    Scalar multiplication.

    INPUT:

    - ``res`` -- a ``celement``, the element to store the result
    - ``x`` -- a ``celement``, the multiplicand
    - ``a`` -- a ``cdigit``, the multiplier
    """
    fmpz_poly_scalar_mul_fmpz(res, x, a)

cdef inline void element_mul(fmpz_poly_t res, fmpz_poly_t x, fmpz_poly_t y) noexcept:
    r"""
    Multiplication.

    INPUT:

    - ``res`` -- a ``celement``, the element to store the result
    - ``x`` -- a ``celement``, the first factor
    - ``y`` -- a ``celement``, the second factor
    """
    fmpz_poly_mul(res, x, y)

cdef inline void element_reduce_digit(fmpz_poly_t x, slong i, PowComputer_flint prime_pow) noexcept:
    r"""
    Reduce the `i`-th digit of `x` and propagate carry.

    INPUT:

    - ``x`` -- a ``celement``, the element to update
    - ``i`` -- integer
    - ``prime_pow`` -- the PowComputer for the ring
    """
    reduce_coeff(x, i, prime_pow.fprime)

cdef inline void element_reducesmall_digit(fmpz_poly_t x, slong i, PowComputer_flint prime_pow) noexcept:
    r"""
    Reduce the `i`-th digit of `x` and propagate carry,
    assuming that `x` is between `0` and `2*p - 1`.

    INPUT:

    - ``x`` -- a ``celement``, the element to update
    - ``i`` -- integer
    - ``prime_pow`` -- the PowComputer for the ring
    """
    reducesmall_coeff(x, i, prime_pow.fprime)

cdef inline void element_reduceneg_digit(fmpz_poly_t x, slong i, PowComputer_flint prime_pow) noexcept:
    r"""
    Reduce the `i`-th digit of `x` and propagate carry,
    assuming that `x` is between `-p` and `p-1`.

    INPUT:

    - ``x`` -- a ``celement``, the element to update
    - ``i`` -- integer
    - ``prime_pow`` -- the PowComputer for the ring
    """
    reduceneg_coeff(x, i, prime_pow.fprime)

cdef inline void element_shift_right(fmpz_poly_t x) noexcept:
    r"""
    Remove the first digit of ``x``.

    INPUT:

    - ``x`` -- a ``celement``
    """
    fmpz_poly_shift_right(x, x, 1)
