"""
Low-level memory allocation functions
"""

#*****************************************************************************
#       Copyright (C) 2011-2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include 'sage/ext/interrupt.pxi'
cimport cython
from libc.stdlib cimport malloc, calloc, realloc, free

cdef extern from *:
    int unlikely(int) nogil  # Defined by Cython


cdef inline void* sage_malloc "sage_malloc"(size_t n) nogil:
    sig_block()
    cdef void* ret = malloc(n)
    sig_unblock()
    return ret


cdef inline void* sage_realloc "sage_realloc"(void* ptr, size_t size) nogil:
    sig_block()
    cdef void* ret = realloc(ptr, size)
    sig_unblock()
    return ret


cdef inline void* sage_calloc "sage_calloc"(size_t nmemb, size_t size) nogil:
    sig_block()
    cdef void* ret = calloc(nmemb, size)
    sig_unblock()
    return ret


cdef inline void sage_free "sage_free"(void* ptr) nogil:
    sig_block()
    free(ptr)
    sig_unblock()


@cython.cdivision(True)
cdef inline size_t mul_overflowcheck(size_t a, size_t b) nogil:
    """
    Return a*b, checking for overflow. Assume that a > 0.
    If overflow occurs, return <size_t>(-1).
    We assume that malloc(<size_t>-1) always fails.
    """
    # If a and b both less than MUL_NO_OVERFLOW, no overflow can occur
    cdef size_t MUL_NO_OVERFLOW = ((<size_t>1) << (4*sizeof(size_t)))
    if a >= MUL_NO_OVERFLOW or b >= MUL_NO_OVERFLOW:
        if unlikely(b > (<size_t>-1) // a):
            return <size_t>(-1)
    return a*b


cdef inline void* check_allocarray(size_t nmemb, size_t size) except? NULL:
    """
    Allocate memory for ``nmemb`` elements of size ``size``.
    """
    if nmemb == 0:
        return NULL
    cdef size_t n = mul_overflowcheck(nmemb, size)
    cdef void* ret = sage_malloc(n)
    if unlikely(ret == NULL):
        raise MemoryError("failed to allocate %s * %s bytes" % (nmemb, size))
    return ret


cdef inline void* check_reallocarray(void* ptr, size_t nmemb, size_t size) except? NULL:
    """
    Re-allocate memory at ``ptr`` to hold ``nmemb`` elements of size
    ``size``. If ``ptr`` equals ``NULL``, this behaves as
    ``check_allocarray``.

    When ``nmemb`` equals 0, then free the memory at ``ptr``.
    """
    if nmemb == 0:
        sage_free(ptr)
        return NULL
    cdef size_t n = mul_overflowcheck(nmemb, size)
    cdef void* ret = sage_realloc(ptr, n)
    if unlikely(ret == NULL):
        raise MemoryError("failed to allocate %s * %s bytes" % (nmemb, size))
    return ret


cdef inline void* check_malloc(size_t n) except? NULL:
    """
    Allocate ``n`` bytes of memory.
    """
    if n == 0:
        return NULL
    cdef void* ret = sage_malloc(n)
    if unlikely(ret == NULL):
        raise MemoryError("failed to allocate %s bytes" % n)
    return ret


cdef inline void* check_realloc(void* ptr, size_t n) except? NULL:
    """
    Re-allocate memory at ``ptr`` to hold ``n`` bytes.
    If ``ptr`` equals ``NULL``, this behaves as ``check_malloc``.
    """
    if n == 0:
        sage_free(ptr)
        return NULL
    cdef void* ret = sage_realloc(ptr, n)
    if unlikely(ret == NULL):
        raise MemoryError("failed to allocate %s bytes" % n)
    return ret


cdef inline void* check_calloc(size_t nmemb, size_t size) except? NULL:
    """
    Allocate memory for ``nmemb`` elements of size ``size``. The
    resulting memory is zeroed.
    """
    if nmemb == 0:
        return NULL
    cdef void* ret = sage_calloc(nmemb, size)
    if unlikely(ret == NULL):
        raise MemoryError("failed to allocate %s * %s bytes" % (nmemb, size))
    return ret
