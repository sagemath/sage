# distutils: libraries = gmp
"""
Low-level memory allocation functions

TESTS:

Check that a ``MemoryError`` is raised if we try to allocate a
ridiculously large integer, see :trac:`15363`::

    sage: 2^(2^63-2)
    Traceback (most recent call last):
    ...
    RuntimeError: exponent must be at most 2147483647          # 32-bit
    MemoryError: failed to allocate 1152921504606847008 bytes  # 64-bit

AUTHORS:

- Jeroen Demeyer (2011-01-13): initial version (:trac:`10258`)

- Jeroen Demeyer (2014-12-14): add more functions (:trac:`10257`)

- Jeroen Demeyer (2015-03-02): move from ``c_lib`` to Cython (:trac:`17881`)
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


cdef extern from "Python.h":
    # Declare as returning void without except value
    void PyErr_Format(object exception, char *format, ...)

cdef extern from "gmp.h":
    void mp_set_memory_functions(
        void *(*) (size_t),
        void *(*) (void *, size_t, size_t),
        void (*) (void *, size_t))


cdef void alloc_error(size_t size) nogil:
    """
    Jump back to ``sig_on()``, raising a ``MemoryError``.
    """
    with gil:
        PyErr_Format(MemoryError, "failed to allocate %zu bytes", size)
    sig_error()


cdef void* sage_sig_malloc(size_t size) nogil:
    """
    ``malloc()`` function for the MPIR/GMP library.

    Out-of-memory errors are handled using the ``sig_error`` mechanism.
    """
    cdef void* p = sage_malloc(size)
    if unlikely(p == NULL):
        alloc_error(size)
    return p


cdef void* sage_sig_realloc(void *ptr, size_t old_size, size_t new_size) nogil:
    """
    ``realloc()`` function for the MPIR/GMP library.

    Out-of-memory errors are handled using the ``sig_error`` mechanism.
    """
    cdef void* p = sage_realloc(ptr, new_size)
    if unlikely(p == NULL):
        alloc_error(new_size)
    return p


cdef void sage_sig_free(void *ptr, size_t size) nogil:
    """
    ``free()`` function for the MPIR/GMP library.
    """
    sage_free(ptr)


def init_memory_functions():
    """
    Set the MPIR/GMP memory functions to the above functions.

    EXAMPLES::

        sage: from sage.ext.memory import init_memory_functions
        sage: init_memory_functions()
    """
    mp_set_memory_functions(sage_sig_malloc, sage_sig_realloc, sage_sig_free)

init_memory_functions()
