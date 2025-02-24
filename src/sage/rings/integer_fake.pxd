# sage_setup: distribution = sagemath-objects
"""
Fake Integer interface

This exists to solve the problem of cyclic imports involving the
``Integer`` class. The problem is that ``Integer`` depends on the
coercion model and the coercion model depends on ``Integer``.

Therefore, this should only be used to implement things at a lower
level than ``Integer``, such as the coercion model.

This provides two functions:

- ``Integer_AS_MPZ(x)``: access the value of the Integer ``x`` as GMP
  ``mpz_t``.

- ``is_Integer(x)``: is ``x`` an Integer?

TESTS::

    sage: cython(                                                                       # needs sage.misc.cython
    ....: '''
    ....: from sage.rings.integer_fake cimport Integer_AS_MPZ, is_Integer
    ....: from sage.rings.integer cimport Integer
    ....: cdef Integer x = Integer(123456789)
    ....: assert is_Integer(x)
    ....: assert Integer_AS_MPZ(x) is x.value
    ....: ''')
"""

#*****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.ref cimport PyTypeObject, Py_TYPE
from sage.libs.gmp.types cimport mpz_ptr

cdef extern from "integer_fake.h":
    PyTypeObject* Integer       # Imported as needed
    mpz_ptr Integer_AS_MPZ(x)
    bint unlikely(bint c)       # Defined by Cython


cdef inline bint is_Integer(x) noexcept:
    global Integer
    if unlikely(Integer is NULL):
        import sage.rings.integer
        Integer = <PyTypeObject*>sage.rings.integer.Integer
    return Py_TYPE(x) is Integer
