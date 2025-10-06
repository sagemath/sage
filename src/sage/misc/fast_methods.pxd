# sage_setup: distribution = sagemath-objects
cdef extern from "Python.h":
    cdef size_t SIZEOF_VOID_P

cdef class FastHashable_class:
    cdef Py_ssize_t _hash

cdef inline long hash_by_id(void * p) noexcept:
    r"""
    This function is a copy paste from the default Python hash function.
    """
    cdef long x
    cdef size_t y = <size_t>p
    # bottom 3 or 4 bits are likely to be 0; rotate y by 4 to avoid
    # excessive hash collisions for dicts and sets
    y = (y >> 4) | (y << (8 * SIZEOF_VOID_P - 4))
    x = <long>y
    if x == -1:
        x = -2
    return x
