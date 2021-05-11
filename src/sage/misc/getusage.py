"""
Get resource usage of process

AUTHORS:

- William Stein (2006-03-04): initial version

- Jeroen Demeyer (2016-11-14): implement as thin wrapper over
  ``psutil`` package
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


import sys


def get_memory_usage(t=None):
    """
    Return the memory usage of the current process in megabytes.

    INPUT:

    - ``t`` -- a float (default: None); output of an earlier call.
      If this is given, return the current usage minus `t`.

    OUTPUT: a float representing the number of megabytes used.

    EXAMPLES::

        sage: t = get_memory_usage(); t  # random
        873.98046875
        sage: type(t)
        <... 'float'>
    """
    import psutil
    m = psutil.Process().memory_info().vms / 1048576
    if t is None:
        return m
    else:
        return m - t


def virtual_memory_limit():
    """
    Return the upper limit for virtual memory usage.

    This is the value set by ``ulimit -v`` at the command line or a
    practical limit if no limit is set. In any case, the value is
    bounded by ``sys.maxsize``.

    OUTPUT:

    Integer. The virtual memory limit in bytes.

    EXAMPLES::

        sage: from sage.misc.getusage import virtual_memory_limit
        sage: virtual_memory_limit() > 0
        True
        sage: virtual_memory_limit() <= sys.maxsize
        True
    """
    import resource
    try:
        vmax = resource.getrlimit(resource.RLIMIT_AS)[0]
    except resource.error:
        vmax = resource.RLIM_INFINITY
    if vmax == resource.RLIM_INFINITY:
        import psutil
        vmax = psutil.virtual_memory().total + psutil.swap_memory().total
    return min(vmax, sys.maxsize)
