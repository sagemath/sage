# distutils: libraries = flint
# distutils: depends = flint/thread_pool.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void thread_pool_init(thread_pool_t T, slong size)
    # Initialise ``T`` and create ``size`` sleeping
    # threads that are available to work.
    # If `size \le 0` no threads are created and future calls to
    # :func:`thread_pool_request` will return `0` (unless
    # :func:`thread_pool_set_size` has been called).

    slong thread_pool_get_size(thread_pool_t T)
    # Return the number of threads in ``T``.

    int thread_pool_set_size(thread_pool_t T, slong new_size)
    # If all threads in ``T`` are in the available state, resize ``T`` and return 1.
    # Otherwise, return ``0``.

    slong thread_pool_request(thread_pool_t T, thread_pool_handle * out, slong requested)
    # Put at most ``requested`` threads in the unavailable state and return
    # their handles. The handles are written to ``out`` and the number of
    # handles written is returned. These threads must be released by a call to
    # ``thread_pool_give_back``.

    void thread_pool_wake(thread_pool_t T, thread_pool_handle i, int max_workers, void (*f)(void*), void * a)
    # Wake up a sleeping thread ``i`` and have it work on ``f(a)``. The thread
    # being woken will be allowed to start ``max_workers`` additional worker
    # threads. Usually this value should be set to ``0``.

    void thread_pool_wait(thread_pool_t T, thread_pool_handle i)
    # Wait for thread ``i`` to finish working and go back to sleep.

    void thread_pool_give_back(thread_pool_t T, thread_pool_handle i)
    # Put thread ``i`` back in the available state. This thread should be sleeping
    # when this function is called.

    void thread_pool_clear(thread_pool_t T)
    # Release any resources used by ``T``. All threads should be given back before
    # this function is called.
