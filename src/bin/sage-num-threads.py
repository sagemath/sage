#!/usr/bin/env python3
#
# Determine the number of threads to be used by Sage.
# This is a simplified version of SAGE_ROOT/build/bin/sage-build-num-threads.py
#
# Outputs three space-separated numbers:
# 1) The number of threads to use for Sage, based on MAKE, MAKEFLAGS
#    and SAGE_NUM_THREADS
# 2) The number of threads to use when parallel execution is explicitly
#    asked for (e.g. sage -tp)
# 3) The number of CPU cores in the system, as determined by
#    multiprocessing.cpu_count()
#
# AUTHOR: Jeroen Demeyer (2011-12-08): Github issue #12016
#
import os
import multiprocessing


def number_of_cores():
    """
    Try to determine the number of CPU cores in this system.
    If successful return that number. Otherwise return 1.
    """
    # If the environment variable SAGE_NUM_CORES exists, use that value.
    # This is useful for testing purposes.
    try:
        n = int(os.environ["SAGE_NUM_CORES"])
        if n > 0:
            return n
    except (ValueError, KeyError):
        pass

    try:
        n = multiprocessing.cpu_count()
        if n > 0:
            return n
    except NotImplementedError:
        pass

    try:  # Solaris fix
        from subprocess import Popen, PIPE
        p = Popen(['sysctl', '-n', 'hw.ncpu'],
                  stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
        n = int(p.stdout.read().strip())
        if n > 0:
            return n
    except (ValueError, OSError):
        pass

    return 1


def num_threads():
    """
    Determine the number of threads from the environment variable
    :envvar:`SAGE_NUM_THREADS`. If it is 0 or not provided, use a default
    of ``min(8, number_of_cores)``.

    OUTPUT:

    a tuple (num_threads, num_threads_parallel, num_cores)
    """
    num_cores = number_of_cores()

    num_threads = None

    # Number of threads to use when parallel execution is explicitly
    # asked for
    num_threads_parallel = num_threads
    if num_threads_parallel is None:
        num_threads_parallel = max(min(8, num_cores), 2)

    try:
        sage_num_threads = int(os.environ["SAGE_NUM_THREADS"])
        if sage_num_threads == 0:
            # If SAGE_NUM_THREADS is 0, use the default only
            # if none of the above variables specified anything.
            if num_threads is None:
                num_threads = min(8, num_cores)
        elif sage_num_threads > 0:
            # SAGE_NUM_THREADS overrides everything
            num_threads = sage_num_threads
            num_threads_parallel = sage_num_threads
    except (ValueError, KeyError):
        pass

    # If we still don't know, use 1 thread
    if num_threads is None:
        num_threads = 1

    # Finally, use SAGE_NUM_THREADS_PARALLEL if set.  A user isn't
    # likely to set this, but it ensures that sage-env is idempotent
    # if called more than once.
    try:
        sage_num_threads = int(os.environ["SAGE_NUM_THREADS_PARALLEL"])
        if sage_num_threads > 0:
            num_threads_parallel = sage_num_threads
    except (ValueError, KeyError):
        pass

    return (num_threads, num_threads_parallel, num_cores)


print(*num_threads())
