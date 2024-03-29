from cpython.number cimport PyNumber_TrueDivide
from sage.structure.element cimport Element


ctypedef fused ulong_or_object:
    unsigned long
    object


cpdef generic_power(a, n) noexcept
cdef generic_power_long(a, long n) noexcept
cdef generic_power_pos(a, ulong_or_object n) noexcept  # n > 0


cdef inline invert(a) noexcept:
    """
    Return ``a^(-1)``.
    """
    if isinstance(a, Element):
        return ~a
    return PyNumber_TrueDivide(type(a)(1), a)


cdef inline one(a) noexcept:
    """
    Return ``a^0``.
    """
    if isinstance(a, Element):
        return (<Element>a)._parent.one()
    return type(a)(1)
