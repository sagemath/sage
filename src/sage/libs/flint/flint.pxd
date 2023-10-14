# distutils: libraries = flint
# distutils: depends = flint/flint.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    mp_limb_t FLINT_BIT_COUNT(mp_limb_t x)
    # Returns the number of binary bits required to represent an ``ulong x``.  If
    # `x` is zero, returns `0`.

    void * flint_malloc(size_t size)

    void * flint_realloc(void * ptr, size_t size)

    void * flint_calloc(size_t num, size_t size)

    void flint_free(void * ptr)

    flint_rand_s * flint_rand_alloc()
    # Allocates a ``flint_rand_t`` object to be used like a heap-allocated
    # ``flint_rand_t`` in external libraries.
    # The random state is not initialised.

    void flint_rand_free(flint_rand_s * state)
    # Frees a random state object as allocated using :func:`flint_rand_alloc`.

    void flint_randinit(flint_rand_t state)
    # Initialize a :type:`flint_rand_t`.

    void flint_randclear(flint_rand_t state)
    # Free all memory allocated by :func:`flint_rand_init`.

    void flint_set_num_threads(int num_threads)
    # Set up a thread pool of ``num_threads - 1`` worker threads (in addition
    # to the master thread) and set the maximum number of worker threads the
    # master thread can start to ``num_threads - 1``.
    # This function may only be called globally from the master thread. It can
    # also be called at a global level to change the size of the thread pool, but
    # an exception is raised if the thread pool is in use (threads have been
    # woken but not given back). The function cannot be called from inside
    # worker threads.

    int flint_get_num_threads()
    # When called at the global level, this function returns one more than the
    # number of worker threads in the Flint thread pool, i.e. it returns the
    # number of workers in the thread pool plus one for the master thread.
    # In general, this function returns one more than the number of additional
    # worker threads that can be started by the current thread.
    # Use :func:`thread_pool_wake` to set this number for a given worker thread.
    # See also: :func:`flint_get_num_available_threads`.

    int flint_set_num_workers(int num_workers)
    # Restricts the number of worker threads that can be started by the current
    # thread to ``num_workers``. This function can be called from any thread.
    # Assumes that the Flint thread pool is already set up.
    # The function returns the old number of worker threads that can be started.
    # The function can only be used to reduce the number of workers that can be
    # started from a thread. It cannot be used to increase the number. If a
    # higher number is passed, the function has no effect.
    # The number of workers must be restored to the original value by a call to
    # :func:`flint_reset_num_workers` before the thread is returned to the thread
    # pool.
    # The main use of this function and :func:`flint_reset_num_workers` is to cheaply
    # and temporarily restrict the number of workers that can be started, e.g. by
    # a function that one wishes to call from a thread, and cheaply restore the
    # number of workers to its original value before exiting the current thread.

    void flint_reset_num_workers(int num_workers)
    # After a call to :func:`flint_set_num_workers` this function must be called to
    # set the number of workers that may be started by the current thread back to
    # its original value.

    int flint_printf(const char * str, ...)
    int flint_fprintf(FILE * f, const char * str, ...)
    int flint_sprintf(char * s, const char * str, ...)
    # These are equivalent to the standard library functions ``printf``,
    # ``vprintf``, ``fprintf``, and ``sprintf`` with an additional length modifier
    # "w" for use with an :type:`mp_limb_t` type. This modifier can be used with
    # format specifiers "d", "x", or "u", thereby outputting the limb as a signed
    # decimal, hexadecimal, or unsigned decimal integer.

    int flint_scanf(const char * str, ...)
    int flint_fscanf(FILE * f, const char * str, ...)
    int flint_sscanf(const char * s, const char * str, ...)
    # These are equivalent to the standard library functions ``scanf``,
    # ``fscanf``, and ``sscanf`` with an additional length modifier "w" for
    # reading an :type:`mp_limb_t` type.
