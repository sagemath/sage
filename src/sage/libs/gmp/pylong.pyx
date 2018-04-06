"""
Various functions to deal with conversion mpz <-> Python int/long

For doctests, see :class:`Integer`.

AUTHORS:

- Gonzalo Tornaria (2006): initial version

- David Harvey (2007-08-18): added ``mpz_get_pyintlong`` function
  (:trac:`440`)

- Jeroen Demeyer (2015-02-24): moved from c_lib, rewritten using
  ``mpz_export`` and ``mpz_import`` (:trac:`17853`)
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from cpython.object cimport Py_SIZE
from cpython.int cimport PyInt_FromLong
from cpython.long cimport PyLong_FromLong
from cpython.longintrepr cimport _PyLong_New, PyLongObject, digit, PyLong_SHIFT
from .mpz cimport *

cdef extern from *:
    Py_ssize_t* Py_SIZE_PTR "&Py_SIZE"(object)


# Unused bits in every PyLong digit
cdef size_t PyLong_nails = 8*sizeof(digit) - PyLong_SHIFT


cdef mpz_get_pylong_large(mpz_srcptr z):
    """
    Convert a non-zero ``mpz`` to a Python ``long``.
    """
    cdef size_t nbits = mpz_sizeinbase(z, 2)
    cdef size_t pylong_size = (nbits + PyLong_SHIFT - 1) // PyLong_SHIFT
    L = _PyLong_New(pylong_size)
    mpz_export((<PyLongObject*>L).ob_digit, NULL,
            -1, sizeof(digit), 0, PyLong_nails, z)
    if mpz_sgn(z) < 0:
        # Set correct size (use a pointer to hack around Cython's
        # non-support for lvalues).
        sizeptr = Py_SIZE_PTR(L)
        sizeptr[0] = -pylong_size
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
        return PyInt_FromLong(mpz_get_si(z))
    return mpz_get_pylong_large(z)


cdef int mpz_set_pylong(mpz_ptr z, L) except -1:
    """
    Convert a Python ``long`` `L` to an ``mpz``.
    """
    cdef Py_ssize_t pylong_size = Py_SIZE(L)
    if pylong_size < 0:
        pylong_size = -pylong_size
    mpz_import(z, pylong_size, -1, sizeof(digit), 0, PyLong_nails,
            (<PyLongObject*>L).ob_digit)
    if Py_SIZE(L) < 0:
        mpz_neg(z, z)


cdef Py_hash_t mpz_pythonhash(mpz_srcptr z):
    """
    Hash an ``mpz``, where the hash value is the same as the hash value
    of the corresponding Python ``long``, except that we do not replace
    -1 by -2 (the Cython wrapper for ``__hash__`` does that).
    """
    # Add all limbs, adding 1 for every carry
    cdef mp_limb_t h1 = 0
    cdef mp_limb_t h0
    cdef size_t i, n
    n = mpz_size(z)
    for i in range(n):
        h0 = h1
        h1 += mpz_getlimbn(z, i)
        # Add 1 on overflow
        if h1 < h0: h1 += 1

    cdef Py_hash_t h = h1
    if mpz_sgn(z) < 0:
        return -h
    return h
