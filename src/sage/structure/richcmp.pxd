# sage_setup: distribution = sagemath-objects
from libc.stdint cimport uint32_t
from cpython.object cimport (Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE,
                             PyObject_RichCompare)


cpdef inline richcmp(x, y, int op):
    """
    Return the result of the rich comparison of ``x`` and ``y`` with
    operator ``op``.

    INPUT:

    - ``x``, ``y`` -- arbitrary Python objects

    - ``op`` -- comparison operator (one of ``op_LT``, ``op_LE``,
      ``op_EQ``, ``op_NE``, ``op_GT``, ``op_GE``)

    EXAMPLES::

        sage: from sage.structure.richcmp import *
        sage: richcmp(3, 4, op_LT)
        True
        sage: richcmp(x, x^2, op_EQ)                                                    # needs sage.symbolic
        x == x^2

    The two examples above are completely equivalent to ``3 < 4``
    and ``x == x^2``. For this reason, it only makes sense in practice
    to call ``richcmp`` with a non-constant value for ``op``.

    We can write a custom ``Element`` class which shows a more
    realistic example of how to use this::

        sage: from sage.structure.element import Element
        sage: class MyElement(Element):
        ....:     def __init__(self, parent, value):
        ....:         Element.__init__(self, parent)
        ....:         self.v = value
        ....:     def _richcmp_(self, other, op):
        ....:         return richcmp(self.v, other.v, op)
        sage: P = Parent()
        sage: x = MyElement(P, 3)
        sage: y = MyElement(P, 3)
        sage: x < y
        False
        sage: x == y
        True
        sage: x > y
        False
    """
    return PyObject_RichCompare(x, y, op)


cpdef richcmp_item(x, y, int op)


cpdef inline richcmp_not_equal(x, y, int op):
    """
    Like ``richcmp(x, y, op)`` but assuming that `x` is not equal to `y`.

    INPUT:

    - ``op`` -- a rich comparison operation (e.g. ``Py_EQ``)

    OUTPUT:

    If ``op`` is not ``op_EQ`` or ``op_NE``, the result of
    ``richcmp(x, y, op)``. If ``op`` is ``op_EQ``, return
    ``False``. If ``op`` is ``op_NE``, return ``True``.

    This is useful to compare lazily two objects A and B according to 2
    (or more) different parameters, say width and height for example.
    One could use::

        return richcmp((A.width(), A.height()), (B.width(), B.height()), op)

    but this will compute both width and height in all cases, even if
    A.width() and B.width() are enough to decide the comparison.

    Instead one can do::

        wA = A.width()
        wB = B.width()
        if wA != wB:
            return richcmp_not_equal(wA, wB, op)
        return richcmp(A.height(), B.height(), op)

    The difference with ``richcmp`` is that ``richcmp_not_equal``
    assumes that its arguments are not equal, which is excluding the case
    where the comparison cannot be decided so far, without
    knowing the rest of the parameters.

    EXAMPLES::

        sage: from sage.structure.richcmp import (richcmp_not_equal,
        ....:    op_EQ, op_NE, op_LT, op_LE, op_GT, op_GE)
        sage: for op in (op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE):
        ....:     print(richcmp_not_equal(3, 4, op))
        True
        True
        False
        True
        False
        False
        sage: for op in (op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE):
        ....:     print(richcmp_not_equal(5, 4, op))
        False
        False
        False
        True
        True
        True
    """
    if op == Py_EQ:
        return False
    elif op == Py_NE:
        return True
    return richcmp(x, y, op)


cpdef inline bint rich_to_bool(int op, int c) noexcept:
    """
    Return the corresponding ``True`` or ``False`` value for a rich
    comparison, given the result of an old-style comparison.

    INPUT:

    - ``op`` -- a rich comparison operation (e.g. ``Py_EQ``)

    - ``c`` -- the result of an old-style comparison: -1, 0 or 1

    OUTPUT: 1 or 0 (corresponding to ``True`` and ``False``)

    .. SEEALSO::

        :func:`rich_to_bool_sgn` if ``c`` could be outside the
        [-1, 0, 1] range.

    EXAMPLES::

        sage: from sage.structure.richcmp import (rich_to_bool,
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
    """
    # op is a value in [0,5], c a value in [-1,1]. We implement this
    # function very efficiently using a bitfield. Note that the masking
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


cpdef inline bint rich_to_bool_sgn(int op, Py_ssize_t c) noexcept:
    """
    Same as ``rich_to_bool``, but allow any `c < 0` and `c > 0`
    instead of only `-1` and `1`.

    .. NOTE::

        This is in particular needed for ``mpz_cmp()``.
    """
    return rich_to_bool(op, (c > 0) - (c < 0))


cpdef inline int revop(int op) noexcept:
    """
    Return the reverse operation of ``op``.

    For example, <= becomes >=, etc.

    EXAMPLES::

        sage: from sage.structure.richcmp import revop
        sage: [revop(i) for i in range(6)]
        [4, 5, 2, 3, 0, 1]
    """
    return (5 - op) ^ 1
