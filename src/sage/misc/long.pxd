r"""
Fast conversion macro to long
"""

from libc.limits cimport LONG_MIN

from cpython.int cimport PyInt_AS_LONG
from cpython.long cimport PyLong_AsLong
from cpython.number cimport PyNumber_Index

from sage.libs.gmp.mpz cimport mpz_fits_slong_p, mpz_get_si

from sage.rings.integer cimport Integer

cdef inline long pyobject_to_long(x) except? LONG_MIN:
    r"""
    Given a Python object ``x`` cast it quickly to a C long.

    A ``TypeError`` is raised if the input can not be converted to an integer or
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
        RuntimeError: exponent must be at most 2147483647            # 32-bit
        RuntimeError: exponent must be at most 9223372036854775807   # 64-bit
    """
    if isinstance(x, int):
        return PyInt_AS_LONG(x)
    elif type(x) is Integer:
        if mpz_fits_slong_p((<Integer>x).value):
            return mpz_get_si((<Integer>x).value)
        else:
            raise OverflowError("Sage Integer too large to convert to C long")
    elif isinstance(x, long):
        return PyLong_AsLong(x)

    return PyNumber_Index(x)
