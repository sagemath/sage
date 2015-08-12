"""
Context manager for safer forking

AUTHOR:

- Jeroen Demeyer (2015-06-20): initial version, :trac:`18741`

"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sys
from posix.unistd cimport getpid, _exit


cdef class ContainChildren(object):
    """
    Context manager which will ensure that all forked child processes
    will be forced to exit if they try to exit the context.

    This can be used as a guard against race conditions, where a child
    process wants to ``fork`` and ``exec`` but it gets interrupted after
    the ``fork`` but before the ``exec``. In such situations, the child
    would raise ``KeyboardInterrupt`` and execute code which is really
    meant for the parent process.

    INPUT:

    - ``exitcode`` -- (integer, default 0) exit code to use when a
      child process tries to exit the with block normally (not due to
      an exception)

    - ``exceptcode`` -- (integer, default 1) exit code to use when a
      child process tries to exit the with block due to an exception

    - ``silent`` -- (boolean, default ``False``) if ``False``, print
      exceptions raised by the child process.

    EXAMPLES::

        sage: from sage.parallel.safefork import ContainChildren
        sage: pid = os.getpid()
        sage: with ContainChildren():
        ....:     child = os.fork()
        sage: assert pid == os.getpid()  # We are the parent process

    By default, exceptions raised by the child process are printed::

        sage: with ContainChildren():
        ....:     child = os.fork()
        ....:     if child == 0:
        ....:         raise RuntimeError("Exception raised by child")
        ....:     _ = os.waitpid(child, 0r)
        Exception raised by child process with pid=...:
        Traceback (most recent call last):
        ...
        RuntimeError: Exception raised by child

    The same example with ``silent=True`` does not show the exception::

        sage: with ContainChildren(silent=True):
        ....:     child = os.fork()
        ....:     if child == 0:
        ....:         raise RuntimeError("Exception raised by child")
        ....:     _ = os.waitpid(child, 0r)
    """
    def __init__(self, exitcode=0, exceptcode=1, silent=False):
        """
        TESTS:

        Check that the exit codes work properly::

            sage: from sage.parallel.safefork import ContainChildren
            sage: with ContainChildren(exitcode=11):
            ....:     child = os.fork()
            sage: pid, st = os.waitpid(child, 0r)
            sage: os.WEXITSTATUS(st)
            11
            sage: with ContainChildren(exceptcode=12, silent=True):
            ....:     child = os.fork()
            ....:     if child == 0:
            ....:         raise RuntimeError("Exception raised by child")
            sage: pid, st = os.waitpid(child, 0r)
            sage: os.WEXITSTATUS(st)
            12
        """
        self.exitcode = exitcode
        self.exceptcode = exceptcode
        self.silent = silent

    def __enter__(self):
        """
        Store the current PID and flush the standard output and error
        streams.

        TESTS:

        The flushing solves the following double-output problem::

            sage: try:
            ....:     sys.stdout.write("X ")
            ....:     if os.fork() == 0:
            ....:         sys.stdout.write("Y ")
            ....:         sys.stdout.flush()
            ....:         os._exit(0)
            ....:     sleep(0.5)  # Give the child process time
            ....:     print("Z")
            ....: finally:
            ....:     pass
            X Y X Z

        With ``ContainChildren()``, no additional flushes are needed::

            sage: from sage.parallel.safefork import ContainChildren
            sage: try:
            ....:     sys.stdout.write("X ")
            ....:     with ContainChildren():
            ....:         sys.stdout.write("Y ")
            ....:     sleep(0.5)  # Give the child process time
            ....:     print("Z")
            ....: finally:
            ....:     pass
            X Y Z
        """
        self.parentpid = getpid()
        # Since we're probably forking, it's a good idea to flush
        # output streams now.
        sys.stdout.flush()
        sys.stderr.flush()

    def __exit__(self, *exc):
        """
        TESTS:

        Check that exceptions raised by the parent work normally::

            sage: from sage.parallel.safefork import ContainChildren
            sage: with ContainChildren():
            ....:     assert os.fork() == 0
            Traceback (most recent call last):
            ...
            AssertionError
        """
        cdef int pid = getpid()
        if pid == self.parentpid:
            # We are the parent process
            return

        # We are a child process!
        cdef int exitcode
        try:
            if exc[0] is None:  # No exception
                exitcode = self.exitcode
            else:               # Exception was raised
                exitcode = self.exceptcode
                if not self.silent:
                    sys.stderr.write("Exception raised by child process with pid=%s:\n"%pid)
                    import traceback
                    traceback.print_exception(*exc)
            sys.stdout.flush()
            sys.stderr.flush()
        finally:
            # This ensures that we can never exit this function alive.
            _exit(exitcode)
