from libc.stdint cimport uint32_t
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE

# Export this for use by Python modules
cdef extern from "Python.h":
    cpdef richcmp "PyObject_RichCompare"(object, object, int)


cdef class SageObject:
    pass


cpdef inline bint rich_to_bool(int op, int c):
    """
    Return the corresponding ``True`` or ``False`` value for a rich
    comparison, given the result of an ordinary comparison.

    INPUT:

    - ``op`` -- a rich comparison operation (e.g. ``Py_EQ``)

    - ``c`` -- the result of an ordinary comparison: -1, 0 or 1.

    OUTPUT: 1 or 0 (corresponding to ``True`` and ``False``)

    .. SEEALSO::

        ``rich_to_bool_sgn`` if ``c`` could be outside the [-1, 1]
        range.

    EXAMPLES::

        sage: from sage.structure.sage_object import (rich_to_bool,
        ....:    op_EQ, op_NE, op_LT, op_LE, op_GT, op_GE)
        sage: for op in (op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE):
        ....:     for c in (-1,0,1):
        ....:         print(rich_to_bool(op, c))
        True False False
        True True False
        False True False
        True False True
        False False True
        False True True

    Indirect tests using integers::

        sage: 0 < 5, 5 < 5, 5 < -8
        (True, False, False)
        sage: 0 <= 5, 5 <= 5, 5 <= -8
        (True, True, False)
        sage: 0 >= 5, 5 >= 5, 5 >= -8
        (False, True, True)
        sage: 0 > 5, 5 > 5, 5 > -8
        (False, False, True)
        sage: 0 == 5, 5 == 5, 5 == -8
        (False, True, False)
        sage: 0 != 5, 5 != 5, 5 != -8
        (True, False, True)

    TESTS::

        sage: from sage.structure.sage_object import py_rich_to_bool
        sage: py_rich_to_bool(op_EQ, 0)
        doctest:...: DeprecationWarning: py_rich_to_bool is deprecated. Please use sage.structure.sage_object.rich_to_bool instead.
        See http://trac.sagemath.org/21128 for details.
        True
    """
    # op is a value in [0,5], c a value in [-1,1]. We implement this
    # function very efficienly using a bitfield. Note that the masking
    # below implies we consider c mod 4, so c = -1 implicitly becomes
    # c = 3.

    # The 4 lines below involve just constants, so the compiler should
    # optimize them to just one constant value for "bits".
    cdef uint32_t less_bits = (1 << Py_LT) + (1 << Py_LE) + (1 << Py_NE)
    cdef uint32_t equal_bits = (1 << Py_LE) + (1 << Py_GE) + (1 << Py_EQ)
    cdef uint32_t greater_bits = (1 << Py_GT) + (1 << Py_GE) + (1 << Py_NE)
    cdef uint32_t bits = (less_bits << 24) + (equal_bits) + (greater_bits << 8)

    cdef int shift = 8*c + op

    # The shift masking (shift & 31) will likely be optimized away by
    # the compiler since shift and bit test instructions implicitly
    # mask their offset.
    return (bits >> (shift & 31)) & 1


cpdef inline bint rich_to_bool_sgn(int op, int c):
    """
    Same as ``rich_to_bool``, but allow any `c < 0` and `c > 0`
    instead of only `-1` and `1`.

    .. NOTE::

        This is in particular needed for ``mpz_cmp()``.
    """
    return rich_to_bool(op, (c > 0) - (c < 0))
