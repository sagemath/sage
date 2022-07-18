# distutils: libraries = flint
# distutils: depends = flint/thread_pool.h

#*****************************************************************************
#       Copyright (C) 2021 Vincent Delecroix
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.flint.types cimport slong

# flint/thread_pool.h
cdef extern from "flint/thread_pool.h":
    ctypedef struct thread_pool_entry_struct:
        pass

    ctypedef thread_pool_entry_struct thread_pool_entry_t[1]

    ctypedef struct thread_pool_struct:
        pass

    ctypedef thread_pool_struct thread_pool_t[1]
    ctypedef int thread_pool_handle


    extern thread_pool_t global_thread_pool
    extern int global_thread_pool_initialized

    void * thread_pool_idle_loop(void * varg)

    void thread_pool_init(thread_pool_t T, slong l)

    int thread_pool_set_affinity(thread_pool_t T,
                                                     int * cpus, slong length)

    int thread_pool_restore_affinity(thread_pool_t T)

    slong thread_pool_get_size(thread_pool_t T)

    int thread_pool_set_size(thread_pool_t T, slong new_size)

    slong thread_pool_request(thread_pool_t T,
                                    thread_pool_handle * out, slong requested)

    void thread_pool_wake(thread_pool_t T, thread_pool_handle i,
                                  int max_workers, void (*f)(void*), void * a)

    void thread_pool_wait(thread_pool_t T, thread_pool_handle i)

    void thread_pool_give_back(thread_pool_t T, thread_pool_handle i)

    void thread_pool_clear(thread_pool_t T)

