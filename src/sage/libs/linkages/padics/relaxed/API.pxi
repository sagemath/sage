# sage_setup: distribution = sagemath-flint
r"""
This file defines the common API for relaxed `p`-adic numbers.

Let `K` be a given `p`-adic ring/field with uniformizer `\pi`.
In the relaxed model, we represent elements of `K` as (Laurent) polynomials over an
explicit exact subring of `K` (whose evaluation at `\pi` are good approximations
of the element we want to represent).

The coefficients of the aforementioned polynomials are called the *digits*.
They are intended to be chosen in a fixed set of small representatives of elements
of the residue field of `K`, but we do not impose this condition because we want
to keep stability under addition and multiplication.
On the contrary, we let the digits be arbitrary elements in the exact subring
but provide functions to reduce digits (i.e. replace them by the representative
of its reduction modulo `\pi`) and propagate the carry accordingly.

Our API is based on two template types (which have to be instantiated in
implementations):

- ``cdigit`` -- the type of a digit

- ``celement`` -- the type of an element (which is then a polynomial whose
  coefficients are elements of type ``cdigit``)

The remainder of the file gives the function signatures.

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


# Operations on digits
######################

cdef cdigit digit_zero
digit_init(digit_zero)

cdef inline void digit_init(cdigit a):
    r"""
    Initialize a digit and set to it the value `0`.

    INPUT:

    - ``a`` -- the ``cdigit`` to initialize
    """
    pass

cdef inline void digit_clear(cdigit a):
    r"""
    Deallocate memory assigned to a digit.

    INPUT:

    - ``a`` -- the ``cdigit`` to deallocate
    """
    pass

# get and set

cdef inline Element digit_get_sage(cdigit a):
    r"""
    Convert a digit to a Sage element.

    INPUT:

    - ``a`` -- a ``cdigit``

    OUTPUT:

    A Sage element representing the same digit.
    """
    pass

cdef inline void digit_set(cdigit a, cdigit b):
    r"""
    Set up a digit.

    INPUT:

    - ``a`` -- the ``cdigit`` to set up
    - ``b`` -- the ``cdigit`` containing the value to assign
    """
    pass

cdef inline void digit_set_ui(cdigit a, long b):
    r"""
    Set an integral value of a digit.

    INPUT:

    - ``a`` -- the ``cdigit`` to set up
    - ``b`` -- integer; the value of assign
    """
    pass

cdef inline void digit_set_sage(cdigit a, Element b):
    r"""
    Set the value of a digit.

    INPUT:

    - ``a`` -- the ``cdigit`` to be assigned
    - ``b`` -- a Sage element containing the value to assign
    """
    pass

# comparisons

cdef inline bint digit_equal(cdigit a, cdigit b):
    r"""
    Comparison of two digits.

    INPUT:

    - ``a`` -- a ``cdigit``
    - ``b`` -- a ``cdigit``

    OUTPUT:

    `1` of `a` is equal to `b`, `0` otherwise
    """
    pass

cdef inline bint digit_equal_ui(cdigit a, long b):
    r"""
    Comparison of a digit and an integer.

    INPUT:

    - ``a`` -- a ``cdigit``
    - ``b`` -- integer

    OUTPUT:

    `1` of `a` is equal to `b`, `0` otherwise
    """
    pass

cdef inline bint digit_is_zero(cdigit a):
    r"""
    Comparison to zero.

    INPUT:

    - ``a`` -- a ``cdigit``

    OUTPUT:

    `1` of `a` vanishes, `0` otherwise
    """
    pass

cdef inline void digit_random_init(randgen generator, long seed):
    r"""
    Initialize the random generator with a new seed.

    INPUT:

    - ``generator`` -- the generator to initialize
    - ``seed`` -- the seed
    """
    pass

cdef inline void digit_random(cdigit res, PowComputer_class prime_pow):
    r"""
    Set a digit to a random value in the distinguished set of representatives.

    INPUT:

    - ``res`` -- the ``cdigit`` to be assigned
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

# operations

cdef inline void digit_add(cdigit res, cdigit a, cdigit b):
    r"""
    Add two digits.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the first summand
    - ``b`` -- a ``cdigit``, the second summand
    """
    pass

cdef inline void digit_sub(cdigit res, cdigit a, cdigit b):
    r"""
    Subtract two digits.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the minuend
    - ``b`` -- a ``cdigit``, the subtrahend
    """
    pass

cdef inline void digit_mul(cdigit res, cdigit a, cdigit b):
    r"""
    Multiply two digits.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the first factor
    - ``b`` -- a ``cdigit``, the second factor
    """
    pass

cdef inline void digit_mod(cdigit res, cdigit a, PowComputer_class prime_pow):
    r"""
    Reduce a digit modulo the uniformizer.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the digit to reduce
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline void digit_quorem(cdigit quo, cdigit rem, cdigit a, PowComputer_class prime_pow):
    r"""
    Reduce a digit modulo the uniformizer and keep the carry.

    INPUT:

    - ``quo`` -- a ``cdigit`` to store the carry
    - ``rem`` -- a ``cdigit`` to store the reduction
    - ``a`` -- a ``cdigit``, the digit to reduce
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline void digit_smallest(cdigit res, cdigit carry, cdigit a, PowComputer_class prime_pow):
    r"""
    Compute the smallest representative of a digit.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the smallest representative
    - ``carry`` -- a ``cdigit`` to store the carry
    - ``a`` -- a ``cdigit``, the digit to reduce
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline void digit_inv(cdigit res, cdigit a, PowComputer_class prime_pow):
    r"""
    Compute the multiplicative inverse of a digit modulo the uniformizer.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the digit to invert
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef bint digit_sqrt(cdigit ans, cdigit x, PowComputer_class prime_pow):
    r"""
    Compute the square root of a digit modulo the uniformizer.

    INPUT:

    - ``res`` -- a ``cdigit`` to store the result
    - ``a`` -- a ``cdigit``, the digit to square root
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass


# Operations on elements
########################

cdef inline void element_init(celement x):
    r"""
    Initialize an element.

    INPUT:

    - ``x`` -- the ``celement`` to initialize
    """
    pass

cdef inline void element_clear(celement x):
    r"""
    Deallocate memory assigned to an element.

    INPUT:

    - ``x`` -- the ``celement`` to deallocate
    """
    pass

# get and set

cdef inline Element element_get_sage(celement x, PowComputer_class prime_pow):
    r"""
    Convert a digit to a Sage element.

    INPUT:

    - ``x`` -- a ``celement``
    - ``prime_pow`` -- the PowComputer for the ring

    OUTPUT:

    A Sage `p`-adic representing the same element.
    """
    pass

cdef inline void element_set(celement x, celement y):
    r"""
    Set an element.

    INPUT:

    - ``x`` -- the ``celement`` to be assigned
    - ``y`` -- the ``celement`` containing the value to assign
    """
    pass

# get and set digits

cdef inline cdigit element_get_digit(celement x, long i):
    r"""
    Return the `i`-th coefficient of `x`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``i`` -- integer
    """
    pass

cdef inline void element_get_slice(celement res, celement x, long start, long length):
    r"""
    Select a slice of an element.

    INPUT:

    - ``res`` -- a ``celement`` to store the slice
    - ``x`` -- a ``celement``, the element from which the slice is extracted
    - ``start`` -- integer; the start position of the slice
    - ``length`` -- integer; the length of the slice

    .. NOTE::

        This function only sets up a pointer to the requested slice
        (the slice is not copied). Hence any future modification
        of the slice ``res`` will affect the container ``x``.
    """
    pass

cdef inline void element_set_digit(celement x, cdigit a, long i):
    r"""
    Set `i`-th coefficient of `x` to the value `a`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``a`` -- a ``cdigit``
    - ``i`` -- integer
    """
    pass

cdef inline void element_set_digit_ui(celement x, long a, long i):
    r"""
    Set `i`-th coefficient of `x` to the value `a`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``a`` -- integer
    - ``i`` -- integer
    """
    pass

cdef inline void element_set_digit_sage(celement x, Element a, long i):
    r"""
    Set `i`-th coefficient of `x` to the value `a`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``a`` -- a Sage element
    - ``i`` -- integer
    """
    pass

# operations

cdef inline void element_iadd_digit(celement x, cdigit a, long i):
    r"""
    Inplace addition:
    add `a` to the `i`-th coefficient of `x`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``a`` -- a ``cdigit``
    - ``i`` -- integer
    """
    pass

cdef inline void element_isub_digit(celement x, cdigit a, long i):
    r"""
    Inplace subtraction:
    subtract `a` to the `i`-th coefficient of `x`.

    INPUT:

    - ``x`` -- a ``celement``
    - ``a`` -- a ``cdigit``
    - ``i`` -- integer
    """
    pass

cdef inline void element_iadd_slice(celement x, celement slice, long start):
    r"""
    Inplace addition:
    add a slice to an element

    INPUT:

    - ``x`` -- a ``celement``, the element to update
    - ``slice`` -- a ``celement``, the slice to be added
    - ``start`` -- the position from which the slice will be added
    """
    pass

cdef inline void element_scalarmul(celement res, celement x, cdigit a):
    r"""
    Scalar multiplication.

    INPUT:

    - ``res`` -- a ``celement``, the element to store the result
    - ``x`` -- a ``celement``, the multiplicand
    - ``a`` -- a ``cdigit``, the multiplier
    """
    pass

cdef inline void element_mul(celement res, celement x, celement y):
    r"""
    Multiplication.

    INPUT:

    - ``res`` -- a ``celement``, the element to store the result
    - ``x`` -- a ``celement``, the first factor
    - ``y`` -- a ``celement``, the second factor
    """
    pass

cdef inline void element_reduce_digit(celement x, long i, PowComputer_class prime_pow):
    r"""
    Reduce the `i`-th digit of `x` and propagate carry.

    INPUT:

    - ``x`` -- a ``celement``, the element to update
    - ``i`` -- integer
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline void element_reducesmall_digit(fmpz_poly_t x, long i, PowComputer_flint prime_pow):
    r"""
    Reduce the `i`-th digit of `x` and propagate carry,
    assuming that `x` is between `0` and `2*p - 1`.

    INPUT:

    - ``x`` -- a ``celement``, the element to update
    - ``i`` -- integer
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline void element_reduceneg_digit(fmpz_poly_t x, long i, PowComputer_flint prime_pow):
    r"""
    Reduce the `i`-th digit of `x` and propagate carry,
    assuming that `x` is between `-p` and `p-1`.

    INPUT:

    - ``x`` -- a ``celement``, the element to update
    - ``i`` -- integer
    - ``prime_pow`` -- the PowComputer for the ring
    """
    pass

cdef inline void element_shift_right(celement x):
    r"""
    Remove the first digit of ``x``.

    INPUT:

    - ``x`` -- a ``celement``
    """
    pass
