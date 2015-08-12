"""
Sage wrapper around pexpect's ``spawn`` class

AUTHOR:

- Jeroen Demeyer (2015-02-01): initial version, :trac:`17686`

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

from pexpect import *

from cpython.ref cimport Py_INCREF
from libc.signal cimport *
from posix.signal cimport killpg
from posix.unistd cimport getpid, getpgid, close, fork

from time import sleep

from sage.parallel.safefork cimport ContainChildren


class SageSpawn(spawn):
    def __init__(self, *args, **kwds):
        """
        Spawn a subprocess in a pseudo-tty.

        - ``*args``, ``**kwds``: see :class:`pexpect.spawn`.

        - ``name`` -- human-readable name for this process, used for
          display purposes only.

        - ``quit_string`` -- (default: ``None``) if not ``None``, send
          this string to the child process before killing it.

        EXAMPLES::

            sage: from sage.interfaces.sagespawn import SageSpawn
            sage: SageSpawn("sleep 1", name="Sleeping Beauty")
            Sleeping Beauty with PID ... running ...
        """
        self.__name = kwds.pop("name", self.__class__.__name__)
        self.quit_string = kwds.pop("quit_string", None)

        with ContainChildren(silent=True):
            spawn.__init__(self, *args, **kwds)


    def __repr__(self):
        """
        A string representing this subprocess.

        EXAMPLES::

            sage: from sage.interfaces.sagespawn import SageSpawn
            sage: s = SageSpawn("true", name="stupid process")
            sage: s  # indirect doctest
            stupid process with PID ... running .../true
            sage: while s.isalive():  # Wait until the process finishes
            ....:     sleep(0.1)
            sage: s  # indirect doctest
            stupid process finished running .../true
        """
        try:
            cmd = " ".join(self.args)
            if not self.terminated:
                return "%s with PID %s running %s" % (self.__name, self.pid, cmd)
            else:
                return "%s finished running %s" % (self.__name, cmd)
        except Exception:
            return object.__repr__(self)

    def _keep_alive(self):
        """
        Ensure that the child process does not get terminated, even if
        all references to ``self`` are deleted.

        This is needed in forked Sage processes: a forked process
        should call this method and then delete all references to
        ``self``.

        EXAMPLES::

            sage: from sage.interfaces.sagespawn import SageSpawn
            sage: s = SageSpawn("sh", ["-c", "while true; do sleep 1; done"])
            sage: s._keep_alive()
            sage: pid = s.pid
            sage: del s
            sage: import gc
            sage: _ = gc.collect()

        All references to ``s`` have been deleted, but the process is
        still running::

            sage: from signal import SIGTERM
            sage: os.kill(pid, SIGTERM)
        """
        Py_INCREF(self)

    def expect_peek(self, *args, **kwds):
        r"""
        Like :meth:`expect` but restore the read buffer such that it
        looks like nothing was actually read. The next reading will
        continue at the current position.

        EXAMPLES::

            sage: from sage.interfaces.sagespawn import SageSpawn
            sage: E = SageSpawn("sh", ["-c", "echo hello world"])
            sage: _ = E.expect_peek("w")
            sage: E.read()
            'hello world\r\n'
        """
        ret = self.expect(*args, **kwds)
        self.buffer = self.before + self.after + self.buffer
        return ret

    def expect_upto(self, *args, **kwds):
        r"""
        Like :meth:`expect` but restore the read buffer starting from
        the matched string. The next reading will continue starting
        with the matched string.

        EXAMPLES::

            sage: from sage.interfaces.sagespawn import SageSpawn
            sage: E = SageSpawn("sh", ["-c", "echo hello world"])
            sage: _ = E.expect_upto("w")
            sage: E.read()
            'world\r\n'
        """
        ret = self.expect(*args, **kwds)
        self.buffer = self.after + self.buffer
        return ret

    def close(self):
        """
        Quit the child process: send the quit string, close the
        pseudo-tty and kill the process.

        This function returns immediately, it doesn't wait for the
        child process to die.

        EXAMPLES::

            sage: from sage.interfaces.sagespawn import SageSpawn
            sage: s = SageSpawn("sleep 1000")
            sage: s.close()
            sage: while s.isalive():  # long time (5 seconds)
            ....:     sleep(0.1)
        """
        cdef int fd = self.child_fd
        if fd != -1:
            if self.quit_string is not None:
                try:
                    # This can fail if the process already exited
                    self.sendline(self.quit_string)
                except OSError:
                    pass
            close(fd)
            self.child_fd = -1
            self.closed = 1
            self.terminate_async()

    def terminate_async(self, interval=5.0):
        """
        Terminate the child process group asynchronously.

        This function returns immediately, while the child is slowly
        being killed in the background.

        INPUT:

        - ``interval`` -- (default: 5) how much seconds to wait between
          sending two signals.

        EXAMPLES:

        Run an infinite loop in the shell::

            sage: from sage.interfaces.sagespawn import SageSpawn
            sage: s = SageSpawn("sh", ["-c", "while true; do sleep 1; done"])

        Check that the process eventually dies after calling
        ``terminate_async``::

            sage: s.terminate_async(interval=0.2)
            sage: while True:
            ....:     try:
            ....:         os.kill(s.pid, 0)
            ....:     except OSError:
            ....:         sleep(0.1)
            ....:     else:
            ....:         break  # process got killed
        """
        cdef int pg, counter = 0
        cdef int thispg = getpgid(getpid())
        assert thispg != -1

        with ContainChildren():
            if fork():
                # Parent process
                return

            # We need to avoid a race condition where the spawned
            # process has not started up completely yet: we need to
            # wait until the spawned process has changed its process
            # group. See Trac #18741.
            pg = getpgid(self.pid)
            while pg == thispg:
                counter += 1
                if counter >= 20:
                    # Something is seriously wrong, give up...
                    raise RuntimeError("%r is not starting up" % self)
                sleep(interval * 0.125)
                pg = getpgid(self.pid)

            # If we failed to determine the process group, probably
            # this child is already dead...
            if pg == -1:
                return

            # If any of these killpg() calls fail, it's most likely
            # because the process is actually killed.
            if killpg(pg, SIGCONT): return
            sleep(interval)
            if killpg(pg, SIGINT): return
            sleep(interval)
            if killpg(pg, SIGHUP): return
            sleep(interval)
            if killpg(pg, SIGTERM): return
            sleep(interval)
            if killpg(pg, SIGKILL): return
