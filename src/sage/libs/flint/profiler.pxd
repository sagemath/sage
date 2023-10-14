# distutils: libraries = flint
# distutils: depends = flint/profiler.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void timeit_start(timeit_t t)
    void timeit_stop(timeit_t t)
    # Gives wall and user time - useful for parallel programming.
    # Example usage::
    # timeit_t t0;
    # // ...
    # timeit_start(t0);
    # // do stuff, take some time
    # timeit_stop(t0);
    # flint_printf("cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

    void start_clock(int n)

    void stop_clock(int n)

    double get_clock(int n)
    # Gives time based on cycle counter.
    # First one must ensure the processor speed in cycles per second
    # is set correctly in ``profiler.h``, in the macro definition
    # ``#define FLINT_CLOCKSPEED``.
    # One can access the cycle counter directly by :func:`get_cycle_counter`
    # which returns the current cycle counter as a ``double``.
    # A sample usage of clocks is::
    # init_all_clocks();
    # start_clock(n);
    # // do something
    # stop_clock(n);
    # flint_printf("Time in seconds is %f.3\n", get_clock(n));
    # where ``n`` is a clock number (from 0-19 by default). The number of
    # clocks can be changed by altering ``FLINT_NUM_CLOCKS``. One can also
    # initialise an individual clock with ``init_clock(n)``.

    void prof_repeat(double *min, double *max, profile_target_t target, void *arg)
    # Allows one to automatically time a given function. Here is a sample usage:
    # Suppose one has a function one wishes to profile::
    # void myfunc(ulong a, ulong b);
    # One creates a struct for passing arguments to our function::
    # typedef struct
    # {
    # ulong a, b;
    # } myfunc_t;
    # a sample function::
    # void sample_myfunc(void * arg, ulong count)
    # {
    # myfunc_t * params = (myfunc_t *) arg;
    # ulong a = params->a;
    # ulong b = params->b;
    # for (ulong i = 0; i < count; i++)
    # {
    # prof_start();
    # myfunc(a, b);
    # prof_stop();
    # }
    # }
    # Then we do the profile::
    # double min, max;
    # myfunc_t params;
    # params.a = 3;
    # params.b = 4;
    # prof_repeat(&min, &max, sample_myfunc, &params);
    # flint_printf("Min time is %lf.3s, max time is %lf.3s\n", min, max);
    # If either of the first two parameters to ``prof_repeat`` is
    # ``NULL``, that value is not stored.
    # One may set the minimum time in microseconds for a timing run by
    # adjusting ``DURATION_THRESHOLD`` and one may set a target duration
    # in microseconds by adjusting ``DURATION_TARGET`` in ``profiler.h``.

    void get_memory_usage(meminfo_t meminfo)
    # Obtains information about the memory usage of the current process.
    # The meminfo object contains the slots ``size`` (virtual memory size),
    # ``peak`` (peak virtual memory size), ``rss`` (resident set size),
    # ``hwm`` (peak resident set size). The values are stored in kilobytes
    # (1024 bytes). This function currently only works on Linux.
