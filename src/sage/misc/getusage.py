"""
Get resource usage of process

AUTHORS:

- William Stein (2006-03-04): initial version

- Jeroen Demeyer (2016-11-14): implement as thin wrapper over
  ``psutil`` package
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, division, print_function

import os
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


########################################################################
# Old deprecated stuff below
########################################################################

def top():
    """
    Return the 'top' or 'prstat' line that contains this running Sage
    process.
    For FreeBSD, return the line containing this running Sage process from
    'ps -axwww -o pid,user,vsz,rss,state,pri,nice,time,cpu,comm'.

    OUTPUT:

    - a string

    EXAMPLES::

        sage: top()              # random output
        '72373 python       0.0%  0:01.36   1    14+  1197   39M+   34M+   55M+  130M+'

    NOTES:

    The external command 'top' (http://www.unixtop.org/) is called on
    Linux, and most other operating systems. The output format of
    'top' is not consistent across all platforms and all versions of
    'top'. If the :func:`top` function does not work in Sage, you may
    need to install 'top'.

    The external command 'prstat' is called on the Solaris and
    OpenSolaris systems. That is part of Solaris, and will not need to
    be installed. The columns used in the 'prstat' output are::

        PID USERNAME  SIZE   RSS STATE  PRI NICE      TIME  CPU PROCESS/NLWP
    """
    from sage.misc.superseded import deprecation
    deprecation(21805, "the function top() is deprecated.")

    U = os.uname()[0].lower()
    pid = os.getpid()

    if U == 'linux':
        cmd = 'top -b -n 1 -p %s' % pid
    elif U == 'darwin':
        cmd = 'top -l 1 |grep "^ *%s "' % pid
    elif U == 'sunos':
        cmd = '/usr/bin/prstat -n 100000 1 1  | grep "^ *%s "' % pid
    elif U[:6] == 'cygwin':
        cmd = 'top -b -n 1 -p %s' % pid
    elif U == 'freebsd':
        cmd = 'ps -axwww -o pid,user,vsz,rss,state,pri,nice,time,cpu,comm | grep "^ *%s "' % pid
    else:
        raise NotImplementedError("top not implemented on platform %s" % U)

    r = os.popen(cmd).read()
    r = r.strip()
    i = r.rfind('\n')
    if i == -1:
        return r
    return r[i+1:]


########################################################################
# The following is adapted from
#   http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/286222
# Python Cookbook, by Jean Brouwers
########################################################################

_proc_status = '/proc/%d/status' % os.getpid()

def VmB(VmKey):
    """
    Function used internally by this module.
    """
    from sage.misc.superseded import deprecation
    deprecation(21805, "the function VmB() is deprecated.")

    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except Exception:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
    return float(v[1])/1024.0

def linux_memory_usage():
    """
    Return memory usage in megabytes.
    """
    return VmB('VmSize:')
