# sage_setup: distribution = sagemath-objects
"""
Various functions to deal with conversion mpz <-> Python int/long

For doctests, see :class:`Integer`.

AUTHORS:

- Gonzalo Tornaria (2006): initial version

- David Harvey (2007-08-18): added ``mpz_get_pyintlong`` function
  (:issue:`440`)

- Jeroen Demeyer (2015-02-24): moved from c_lib, rewritten using
  ``mpz_export`` and ``mpz_import`` (:issue:`17853`)
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from cpython.long cimport PyLong_FromLong
from cpython.longintrepr cimport _PyLong_New, py_long, digit, PyLong_SHIFT
from sage.cpython.pycore_long cimport (ob_digit, _PyLong_IsNegative,
    _PyLong_DigitCount, _PyLong_SetSignAndDigitCount)
from sage.libs.gmp.mpz cimport *

cdef extern from *:
    """
    /* Compatibility for python 3.8, can be removed later */
    #if PY_VERSION_HEX < 0x030900A4 && !defined(Py_SET_SIZE)
    static inline void _Py_SET_SIZE(PyVarObject *ob, Py_ssize_t size)
    { ob->ob_size = size; }
    #define Py_SET_SIZE(ob, size) _Py_SET_SIZE((PyVarObject*)(ob), size)
    #endif
    """
    void Py_SET_SIZE(object, Py_ssize_t)
    int hash_bits """
        #ifdef _PyHASH_BITS
        _PyHASH_BITS         /* Python 3 */
        #else
        (8 * sizeof(void*))  /* Python 2 */
        #endif
        """
    int limb_bits "(8 * sizeof(mp_limb_t))"


# Unused bits in every PyLong digit
cdef size_t PyLong_nails = 8*sizeof(digit) - PyLong_SHIFT


cdef mpz_get_pylong_large(mpz_srcptr z):
    """
    Convert a nonzero ``mpz`` to a Python ``long``.
    """
    cdef size_t nbits = mpz_sizeinbase(z, 2)
    cdef size_t pylong_size = (nbits + PyLong_SHIFT - 1) // PyLong_SHIFT
    cdef py_long L = _PyLong_New(pylong_size)
    mpz_export(ob_digit(L), NULL, -1, sizeof(digit), 0, PyLong_nails, z)
    _PyLong_SetSignAndDigitCount(L, mpz_sgn(z), pylong_size)
    return L


cdef mpz_get_pylong(mpz_srcptr z):
    """
    Convert an ``mpz`` to a Python ``long``.
    """
    if mpz_fits_slong_p(z):
        return PyLong_FromLong(mpz_get_si(z))
    return mpz_get_pylong_large(z)


cdef mpz_get_pyintlong(mpz_srcptr z):
    """
    Convert an ``mpz`` to a Python ``int`` if possible, or a ``long``
    if the value is too large.
    """
    if mpz_fits_slong_p(z):
        return PyLong_FromLong(mpz_get_si(z))
    return mpz_get_pylong_large(z)


cdef int mpz_set_pylong(mpz_ptr z, py_long L) except -1:
    """
    Convert a Python ``long`` `L` to an ``mpz``.
    """
    cdef Py_ssize_t pylong_size = _PyLong_DigitCount(L)
    mpz_import(z, pylong_size, -1, sizeof(digit), 0, PyLong_nails, ob_digit(L))
    if _PyLong_IsNegative(L):
        mpz_neg(z, z)


cdef Py_hash_t mpz_pythonhash(mpz_srcptr z) noexcept:
    """
    Hash an ``mpz``, where the hash value is the same as the hash value
    of the corresponding Python ``int`` or ``long``, except that we do
    not replace -1 by -2 (the Cython wrapper for ``__hash__`` does that).
    """
    if mpz_sgn(z) == 0:
        return 0

    # The hash value equals z % m where m = 2 ^ hash_bits - 1.
    #
    # Safely compute 2 ^ hash_bits - 1 without overflow
    cdef mp_limb_t modulus = (((<mp_limb_t>(1) << (hash_bits - 1)) - 1) * 2) + 1

    cdef mp_limb_t h = 0
    cdef mp_limb_t x, y
    cdef size_t i, n
    cdef unsigned int r
    n = mpz_size(z)
    for i in range(n):
        x = mpz_getlimbn(z, i)

        # Computing modulo 2 ^ hash_bits - 1 means that the bit at
        # position j is really moved to position (j % hash_bits).
        # We need to shift every bit of x left by (limb_bits * i)
        # and then put it in the right position to account for
        # the modulo operation. Store the result in y.
        if limb_bits == hash_bits:
            y = x
        else:
            r = (limb_bits * i) % hash_bits
            y = (x << r) & modulus
            y += (x >> (hash_bits - r)) & modulus
            # Only do this shift if we don't shift more than the size of the
            # type
            if r > 2 * hash_bits - limb_bits:
                y += (x >> (2 * hash_bits - r))
            # At this point, y <= 2 * modulus, so y did not overflow, but we
            # need y <= modulus. We use > instead of >= on the line below
            # because it generates more efficient code.
            if y > modulus:
                y -= modulus

        # Safely compute h = (h + y) % modulus knowing that h < modulus
        # and y <= modulus
        if h < modulus - y:
            h = h + y
        else:
            h = h - (modulus - y)

    # Special case for Python 2
    if limb_bits == hash_bits and h == 0:
        h = -1

    if mpz_sgn(z) < 0:
        return -h
    return h


assert hash_bits <= limb_bits <= 2 * hash_bits
