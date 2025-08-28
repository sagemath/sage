# sage_setup: distribution = sagemath-objects
# sage.doctest: needs sage.misc.cython
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
from cpython.number cimport PyNumber_Index, PyIndex_Check
from cpython.longintrepr cimport py_long, PyLong_SHIFT, digit
from sage.cpython.pycore_long cimport (
    ob_digit, _PyLong_IsNegative, _PyLong_DigitCount)

from sage.libs.gmp.mpz cimport mpz_fits_slong_p, mpz_get_si
from sage.rings.integer_fake cimport is_Integer, Integer_AS_MPZ


cdef inline long pyobject_to_long(x) except? LONG_MIN:
    r"""
    Given a Python object ``x`` cast it quickly to a C long.

    A :exc:`TypeError` is raised if the input cannot be converted to
    an integer or
    an :exc:`OverflowError` is raised if it does not fit into a C long.

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

    See :issue:`22319`::

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
    Python type ``int``, for Sage Integers and for types
    implementing ``__index__``.

    If possible, compute the value of this integer as C long and store
    it in ``*value``.

    Errors are returned as an error indicator ``*err`` (without raising
    any Python exception).

    Possible errors when returning ``True``:

    - ``0`` -- ``x`` was successfully converted to a C long and its value
      is stored in ``*value``

    - ``ERR_OVERFLOW`` -- ``x`` is an integer type but too large to store
      in a C long

    Possible errors when returning ``False``:

    - ``ERR_TYPE`` -- ``x`` is not an integer type of any kind

    - ``ERR_INDEX`` -- ``x`` implements ``__index__`` but a :exc:`TypeError`
      was raised calling ``__index__()``

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
        ....:             raise OverflowError(f"integer_check_long: overflow ({x})")
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
        sage: L = [1, 12345, 10^9, 2^30, long_max()//9, long_max()//3, long_max()]
        sage: L += [-x for x in L] + [0, long_min()]
        sage: for v in L:
        ....:     for t in (Integer, int, QQ):
        ....:         assert check_long(t(v)) == v
        sage: check_long(2^100)
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow (...)
        sage: check_long(long_max() + 1)
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow (...)
        sage: check_long(long_min() - 1)
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow (...)
        sage: check_long("hello")
        Traceback (most recent call last):
        ...
        TypeError: integer_check_long: wrong type
        sage: check_long(2/3)
        Traceback (most recent call last):
        ...
        TypeError: integer_check_long: bad __index__

    Repeat the overflow tests with python integers:

        sage: check_long(int(2^100))
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow (...)
        sage: check_long(int(long_max() + 1))
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow (...)
        sage: check_long(int(long_min() - 1))
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow (...)

    And again with rationals:

        sage: check_long(QQ(2^100))
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow (...)
        sage: check_long(QQ(long_max() + 1))
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow (...)
        sage: check_long(QQ(long_min() - 1))
        Traceback (most recent call last):
        ...
        OverflowError: integer_check_long: overflow (...)
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


cdef inline long dig(const digit* D, int n) noexcept:
    # Convenient helper function for integer_check_long_py()
    return (<long>D[n]) << (n * PyLong_SHIFT)


cdef inline bint integer_check_long_py(x, long* value, int* err) noexcept:
    """
    Return whether ``x`` is a Python object of type ``int``.

    If possible, compute the value of this integer as C long and store
    it in ``*value``.

    Errors are returned as an error indicator ``*err`` (without raising
    any Python exception).

    Possible errors when returning ``True``:

    - ``0`` -- ``x`` was successfully converted to a C long and its value
      is stored in ``*value``

    - ``ERR_OVERFLOW`` -- ``x`` is a Python object of type ``int`` but
      too large to store in a C long

    Possible errors when returning ``False``:

    - ``ERR_TYPE`` -- ``x`` is not a Python object of type ``int``

    EXAMPLES:

    We create a pure Python wrapper of this function::

        sage: cython('''
        ....: from sage.arith.long cimport *
        ....: def check_long_py(x):
        ....:     cdef long value
        ....:     cdef int err
        ....:     cdef bint c = integer_check_long_py(x, &value, &err)
        ....:     if c:
        ....:         if err == 0:
        ....:             return value
        ....:         elif err == ERR_OVERFLOW:
        ....:             return f"Overflow ({x})"
        ....:     elif err == ERR_TYPE:
        ....:         return f"Bad type ({x})"
        ....:     return f"This should never happen ({x})"
        ....: from libc.limits cimport LONG_MIN, LONG_MAX
        ....: def long_min():
        ....:     return LONG_MIN
        ....: def long_max():
        ....:     return LONG_MAX
        ....: ''')
        sage: L = [1, 12345, 10^9, 2^30, long_max()//9, long_max()//3, long_max()]
        sage: L += [-x for x in L] + [0, long_min()]
        sage: for v in L:
        ....:     assert check_long_py(int(v)) == v
        sage: check_long_py(int(2^60))
        1152921504606846976                 # 64-bit
        'Overflow (...)'                    # 32-bit
        sage: check_long_py(int(2^61))
        2305843009213693952                 # 64-bit
        'Overflow (...)'                    # 32-bit
        sage: check_long_py(int(2^62))
        4611686018427387904                 # 64-bit
        'Overflow (...)'                    # 32-bit
        sage: check_long_py(int(2^63))
        'Overflow (...)'
        sage: check_long_py(int(2^100))
        'Overflow (...)'
        sage: check_long_py(int(long_max() + 1))
        'Overflow (...)'
        sage: check_long_py(int(long_min() - 1))
        'Overflow (...)'
        sage: check_long_py(389)
        'Bad type (...)'
        sage: check_long_py("hello")
        'Bad type (...)'
        sage: check_long_py(2/3)
        'Bad type (...)'
    """
    if not isinstance(x, int):
        err[0] = ERR_TYPE
        return 0

    # x is a Python "int" (aka PyLongObject or py_long in cython)
    cdef const digit* D = ob_digit(x)
    cdef Py_ssize_t size = _PyLong_DigitCount(x)

    if _PyLong_IsNegative(x):
        size = -size

    # We assume PyLong_SHIFT <= BITS_IN_LONG <= 3 * PyLong_SHIFT.
    # This is true in all the default configurations:
    # - BITS_IN_LONG = 63, PyLong_SHIFT = 30
    # - BITS_IN_LONG = 31, PyLong_SHIFT = 15 (python <= 3.10)
    # - BITS_IN_LONG = 31, PyLong_SHIFT = 30 (new in python 3.11)
    # cf. https://github.com/sagemath/sage/issues/33842#comment:130
    #
    # This way, we know that 1 digit certainly fits in a C long
    # and 4 or more digits never fit.
    # For 2 or 3 digits, we need an explicit overflow check.
    cdef int BITS_IN_LONG = 8 * sizeof(long) - 1
    if not (PyLong_SHIFT <= BITS_IN_LONG <= 3 * PyLong_SHIFT):
        raise AssertionError(
                f"PyLong_SHIFT = {PyLong_SHIFT}, "
                f"BITS_IN_LONG = {BITS_IN_LONG}")

    cdef long lead
    cdef long lead_2_overflow = (<long>1) << (BITS_IN_LONG - PyLong_SHIFT)
    cdef long lead_3_overflow
    if BITS_IN_LONG < 2 * PyLong_SHIFT:
        # in this case 3 digit is always overflow
        lead_3_overflow = 0
    else:
        lead_3_overflow = (<long>1) << (BITS_IN_LONG - 2 * PyLong_SHIFT)
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
        if BITS_IN_LONG < 2 * PyLong_SHIFT and D[1] >= lead_2_overflow:
            err[0] = ERR_OVERFLOW
            return 1
        value[0] = dig(D, 0) + dig(D, 1)
        err[0] = 0
    elif size == -2:
        if BITS_IN_LONG < 2 * PyLong_SHIFT and D[1] >= lead_2_overflow:
            if D[0] == 0 and D[1] == lead_2_overflow:
                # Special case for LONG_MIN
                value[0] = (<long>-1) << BITS_IN_LONG
                err[0] = 0
            else:
                err[0] = ERR_OVERFLOW
            return 1
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


cdef inline bint is_small_python_int(obj) noexcept:
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
