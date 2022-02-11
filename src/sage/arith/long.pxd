r"""
Fast conversion of Python objects to C long
"""

# ****************************************************************************
#       Copyright (C) 2015 Vincent Delecroix <20100.delecroix@gmail.com>
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.limits cimport LONG_MIN, LONG_MAX

from cpython.object cimport Py_SIZE
from cpython.int cimport PyInt_AS_LONG
from cpython.long cimport PyLong_AsLong
from cpython.number cimport PyNumber_Index, PyIndex_Check
from cpython.longintrepr cimport py_long, PyLong_SHIFT, digit

from sage.libs.gmp.mpz cimport mpz_fits_slong_p, mpz_get_si
from sage.rings.integer_fake cimport is_Integer, Integer_AS_MPZ


cdef inline long pyobject_to_long(x) except? LONG_MIN:
    r"""
    Given a Python object ``x`` cast it quickly to a C long.

    A ``TypeError`` is raised if the input cannot be converted to an integer or
    an ``OverflowError`` is raised if it does not fit into a C long.

    TESTS:

    We test indirectly that ``Integer.__pow__`` works::

        sage: a = 10
        sage: a^10
        10000000000
        sage: a^(10r)
        10000000000
        sage: a^(10l)
        10000000000
        sage: a^(10/1)
        10000000000
        sage: a^(2**258)
        Traceback (most recent call last):
        ...
        OverflowError: exponent must be at most 2147483647           # 32-bit
        OverflowError: exponent must be at most 9223372036854775807  # 64-bit

    See :trac:`22319`::

        sage: a^pari(10)
        10000000000
    """
    cdef long value
    cdef int err
    if integer_check_long_py(x, &value, &err):
        if err:
            raise OverflowError("Python int too large to convert to C long")
        return value
    if is_Integer(x):
        z = Integer_AS_MPZ(x)
        if mpz_fits_slong_p(z):
            return mpz_get_si(z)
        else:
            raise OverflowError("Sage Integer too large to convert to C long")

    return PyNumber_Index(x)


# Error values for integer_check_long()
cdef enum:
    ERR_TYPE = 1
    ERR_INDEX = 2
    ERR_OVERFLOW = 3


cdef inline bint integer_check_long(x, long* value, int* err) except -1:
    """
    Return whether ``x`` is some integer type. This is true for the
    Python types ``int`` and ``long``, for Sage Integers and for types
    implementing ``__index__``.

    If possible, compute the value of this integer as C long and store
    it in ``*value``.

    Errors are returned as an error indicator ``*err`` (without raising
    any Python exception).

    Possible errors when returning ``True``:

    - ``0``: ``x`` was successfully converted to a C long and its value
      is stored in ``*value``.

    - ``ERR_OVERFLOW``: ``x`` is an integer type but too large to store
      in a C long.

    Possible errors when returning ``False``:

    - ``ERR_TYPE``: ``x`` is not an integer type of any kind.

    - ``ERR_INDEX``: ``x`` implements ``__index__`` but a ``TypeError``
      was raised calling ``__index__()``.

    - Other exceptions in ``__index__`` are simply propagated. This is
      the only way this function can raise an exception.

    EXAMPLES:

    We create a pure Python wrapper of this function::

        sage: cython('''
        ....: from sage.arith.long cimport *
        ....: from sage.rings.integer cimport smallInteger
        ....: def check_long(x):
        ....:     cdef long value
        ....:     cdef int err
        ....:     cdef bint c = integer_check_long(x, &value, &err)
        ....:     if c:
        ....:         if err == 0:
        ....:             return value
        ....:         elif err == ERR_OVERFLOW:
        ....:             raise OverflowError("integer_check_long: overflow")
        ....:     elif err == ERR_TYPE:
        ....:         raise TypeError("integer_check_long: wrong type")
        ....:     elif err == ERR_INDEX:
        ....:         raise TypeError("integer_check_long: bad __index__")
        ....:     assert False
        ....: from libc.limits cimport LONG_MIN, LONG_MAX
        ....: def long_min():
        ....:     return smallInteger(LONG_MIN)
        ....: def long_max():
        ....:     return smallInteger(LONG_MAX)
        ....: ''')
        sage: types = (ZZ, QQ, int)
        sage: L = [1, 12345, 10^9, 2^30, long_max()//9, long_max()//3, long_max()]
        sage: L += [-x for x in L] + [0, long_min()]
        sage: for v in L:
        ....:     for t in (Integer, int):
        ....:         assert check_long(t(v)) == v
        sage: check_long(2^100)
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow
        sage: check_long(long_max() + 1)
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow
        sage: check_long(long_min() - 1)
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow
        sage: check_long("hello")
        Traceback (most recent call last):
        ...
        TypeError: integer_check_long: wrong type
        sage: check_long(2/3)
        Traceback (most recent call last):
        ...
        TypeError: integer_check_long: bad __index__
    """
    cdef int c = integer_check_long_py(x, value, err)
    if c:
        return c
    if is_Integer(x):
        z = Integer_AS_MPZ(x)
        if mpz_fits_slong_p(z):
            value[0] = mpz_get_si(z)
            err[0] = 0
        else:
            err[0] = ERR_OVERFLOW
        return 1
    elif PyIndex_Check(x):
        err[0] = ERR_INDEX
        try:
            x = PyNumber_Index(x)
        except TypeError:
            return 0
        return integer_check_long_py(x, value, err)
    else:
        err[0] = ERR_TYPE
        return 0


cdef inline long dig(const digit* D, int n):
    # Convenient helper function for integer_check_long_py()
    return (<long>D[n]) << (n * PyLong_SHIFT)


cdef inline bint integer_check_long_py(x, long* value, int* err):
    """
    Part of ``integer_check_long`` in ``long.pxd``, checking only for
    Python objects of type ``int`` and ``long``. See that function for
    documentation and tests.
    """
    if not isinstance(x, long):
        if isinstance(x, int):
            # This can happen only on Python 2
            value[0] = PyInt_AS_LONG(x)
            err[0] = 0
            return 1
        err[0] = ERR_TYPE
        return 0

    # x is a Python "long" (called "int" on Python 3)
    cdef const digit* D = (<py_long>x).ob_digit
    cdef Py_ssize_t size = Py_SIZE(x)

    # We assume that PyLong_SHIFT is 15 on a 32-bit system and 30 on a
    # 64-bit system. This is not guaranteed by Python, but it is the
    # default configuration.
    #
    # This way, we know that 1 and 2 digits certainly fit in a C long
    # and 4 or more digits never fit. For 3 digits, we need an explicit
    # overflow check.
    cdef int BITS_IN_LONG = 8 * sizeof(long) - 1
    if not (2 * PyLong_SHIFT <= BITS_IN_LONG < 4 * PyLong_SHIFT):
        raise AssertionError

    cdef long lead
    cdef long lead_3_overflow = (<long>1) << (BITS_IN_LONG - 2 * PyLong_SHIFT)
    if size == 0:
        value[0] = 0
        err[0] = 0
    elif size == 1:
        value[0] = dig(D, 0)
        err[0] = 0
    elif size == -1:
        value[0] = -dig(D, 0)
        err[0] = 0
    elif size == 2:
        value[0] = dig(D, 0) + dig(D, 1)
        err[0] = 0
    elif size == -2:
        value[0] = -(dig(D, 0) + dig(D, 1))
        err[0] = 0
    elif size == 3:
        lead = D[2]
        if lead < lead_3_overflow:
            value[0] = dig(D, 0) + dig(D, 1) + dig(D, 2)
            err[0] = 0
        else:
            err[0] = ERR_OVERFLOW
    elif size == -3:
        lead = D[2]
        if lead < lead_3_overflow:
            value[0] = -(dig(D, 0) + dig(D, 1) + dig(D, 2))
            err[0] = 0
        elif D[0] == 0 and D[1] == 0 and lead == lead_3_overflow:
            # Special case for LONG_MIN
            value[0] = (<long>-1) << BITS_IN_LONG
            err[0] = 0
        else:
            err[0] = ERR_OVERFLOW
    else:
        # 4 digits or more: guaranteed overflow
        err[0] = ERR_OVERFLOW
    return 1


cdef inline bint is_small_python_int(obj):
    """
    Test whether Python object is a small Python integer.

    Meaning that it can be converted to a C long. In Python 2,
    this is equivalent to it being the ``int`` Python type. In Python
    3, the ``int`` Python type has unlimited precision so we need to
    check its range.

    EXAMPLES::

        sage: cython('''
        ....: from sage.arith.long cimport is_small_python_int
        ....: def is_small_wrapper(x):
        ....:     return is_small_python_int(x)
        ....: ''')
        sage: is_small_wrapper(int(3))
        True
        sage: is_small_wrapper(ZZ(3))   # not a Python int
        False
        sage: import sys
        sage: is_small_wrapper(int(sys.maxsize))   # does fit into C long
        True
        sage: is_small_wrapper(int(sys.maxsize + 1))   # does not fit into C long
        False
    """
    return (type(obj) is int) and (LONG_MIN <= obj <= LONG_MAX)
