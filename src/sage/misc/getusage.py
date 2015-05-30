"""
Get resource usage of process

AUTHORS:

- William Stein (2006-03-04): initial version
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

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

def get_memory_usage(t=None):
    """
    Return memory usage.

    INPUT:

    - ``t`` - a float (default: None); output of an earlier call

    OUTPUT:

    - ``Linux`` - Returns float number (in megabytes)

    - ``OS X`` - Returns float number (in megabytes) that matches
      VSIZE column of ``top``

    - ``Solaris or OpenSolaris`` - Returns float number (in megabytes)
      that matches RSS column of ``prstat``. Depending on the memory
      usage, ``prstat`` will output the data in KB, MB or GB. In each
      case, the value returned by this function will always be in MB.

    - ``FreeBSD`` - Returns float number (in megabytes) that matches
      RSS column of ``ps -auxwww``

    - ``other`` - not implemented for any other operating systems

    EXAMPLES::

        sage: t = get_memory_usage(); t  # random
        873.98046875

    .. NOTE::

        * Currently, :func:`get_memory_usage` calls ``prstat`` on Solaris
          and OpenSolaris to get the data it requires. In the long term, a
          better solution would be to use Solaris system calls.

        * In some instances, ``top`` may be used on OS X. This may break
          if the memory usage is greater than 9999 MB. However, normally
          ``top`` is not used on OS X.
    """
    U = os.uname()[0].lower()
    if U == 'linux' or U[:6] == 'cygwin':
        m = linux_memory_usage()
    elif U == 'darwin':
        try:
            from sage.misc.darwin_utilities import darwin_memory_usage
            m = float(darwin_memory_usage()) / (1024 * 1024)
        except ImportError:
            # darwin_utilities is not supported on some versions of OS X.
            m = float(top().split()[-1].strip('M+'))
    elif U == 'sunos':
        # Sun's 'prstat' command appends K, M or G depending on whether
        # the memory usage is in KB. MB or GB. So we need to strip off
        # the letter, and convert to a consistent unit of MB.
        memory_in_KB_MB_or_GB = top().split()[3]
        if memory_in_KB_MB_or_GB.endswith("K"):
            m = float(memory_in_KB_MB_or_GB.strip("K")) / 1024
        elif memory_in_KB_MB_or_GB.endswith("M"):
            m = float(memory_in_KB_MB_or_GB.strip("M"))
        elif memory_in_KB_MB_or_GB.endswith("G"):
            m = float(memory_in_KB_MB_or_GB.strip("G")) * 1024
    elif U == 'freebsd':
        memory_in_KB = top().split()[3]
        m = float(memory_in_KB) / 1024
    else:
        raise NotImplementedError("memory usage not implemented on platform %s" % U)

    if t is None:
        return m
    else:
        return m - t



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
