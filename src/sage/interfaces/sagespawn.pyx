# sage.doctest: optional - ptyprocess
"""
Sage wrapper around pexpect's ``spawn`` class and
the ptyprocess's ``PtyProcess`` class.

AUTHOR:

- Jeroen Demeyer (2015-02-01): initial version, see :issue:`17686`.

- Jeroen Demeyer (2015-12-04): add support for pexpect 4 + ptyprocess,
  see :issue:`10295`.
"""

# ***************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from pexpect import *
from ptyprocess import PtyProcess

from cpython.ref cimport Py_INCREF
from libc.signal cimport *
from posix.signal cimport killpg
from posix.unistd cimport getpid, getpgid, fork

from time import sleep

from sage.cpython.string cimport str_to_bytes
from sage.interfaces.process cimport ContainChildren


class SageSpawn(spawn):
    def __init__(self, *args, **kwds):
        """
        Spawn a subprocess in a pseudo-tty.

        - ``*args``, ``**kwds`` -- see :class:`pexpect.spawn`

        - ``name`` -- human-readable name for this process, used for
          display purposes only

        - ``quit_string`` -- (default: ``None``) if not ``None``, send
          this string to the child process before killing it

        EXAMPLES::

            sage: from sage.interfaces.sagespawn import SageSpawn
            sage: SageSpawn("sleep 1", name="Sleeping Beauty")
            Sleeping Beauty with PID ... running ...
        """
        self.__name = kwds.pop("name", self.__class__.__name__)
        self.quit_string = kwds.pop("quit_string", None)

        kwds.setdefault("ignore_sighup", True)

        # Use a *serious* read buffer of 4MiB, not the ridiculous 2000 bytes
        # that pexpect uses by default.
        kwds.setdefault("maxread", 4194304)

        with ContainChildren(silent=True):
            spawn.__init__(self, *args, **kwds)

        self.delaybeforesend = None
        self.delayafterread = None

    def _spawnpty(self, args, **kwds):
        """
        Create an instance of :class:`SagePtyProcess`.

        EXAMPLES::

            sage: from sage.interfaces.sagespawn import SageSpawn
            sage: s = SageSpawn("sleep 1")
            sage: s.ptyproc
            SagePtyProcess.spawn(...)
        """
        ptyproc = SagePtyProcess.spawn(args, **kwds)
        ptyproc.quit_string = self.quit_string
        return ptyproc

    def __repr__(self):
        """
        A string representing this subprocess.

        EXAMPLES::

            sage: from sage.interfaces.sagespawn import SageSpawn
            sage: s = SageSpawn("true", name="stupid process")
            sage: s  # indirect doctest
            stupid process with PID ... running .../true
            sage: while s.isalive():  # Wait until the process finishes
            ....:     sleep(float(0.1))
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
            sage: E.read().decode('ascii')
            'hello world\r\n'
        """
        ret = self.expect(*args, **kwds)
        self._before = self.buffer_type()
        self._before.write(self.before + self.after + self.buffer)
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
            sage: E.read().decode('ascii')
            'world\r\n'
        """
        ret = self.expect(*args, **kwds)
        self._before = self.buffer_type()
        self._before.write(self.after + self.buffer)
        return ret


class SagePtyProcess(PtyProcess):
    def close(self, force=None):
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
            ....:     sleep(float(0.1))
        """
        if not self.closed:
            if self.quit_string is not None:
                try:
                    # This can fail if the process already exited
                    # PtyProcess.write takes bytes; ideally we would use
                    # an encoding picked specifically for the target process
                    # but the default (UTF-8) will do now, since I don't
                    # think we have any non-ASCII quit_strings anyways.
                    self.write(str_to_bytes(self.quit_string))
                except (OSError, IOError):
                    pass
            self.fileobj.close()
            self.fd = -1
            self.closed = 1
            self.terminate_async()

    def terminate_async(self, interval=5.0):
        """
        Terminate the child process group asynchronously.

        This function returns immediately, while the child is slowly
        being killed in the background.

        INPUT:

        - ``interval`` -- (default: 5) how much seconds to wait between
          sending two signals

        EXAMPLES:

        Run an infinite loop in the shell::

            sage: from sage.interfaces.sagespawn import SageSpawn
            sage: s = SageSpawn("sh", ["-c", "while true; do sleep 1; done"])

        Check that the process eventually dies after calling
        ``terminate_async``::

            sage: s.ptyproc.terminate_async(interval=float(0.2))
            sage: while True:
            ....:     try:
            ....:         os.kill(s.pid, 0)
            ....:     except OSError:
            ....:         sleep(float(0.1))
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
            # group. See Issue #18741.
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
            if killpg(pg, SIGCONT):
                return
            sleep(interval)
            if killpg(pg, SIGINT):
                return
            sleep(interval)
            if killpg(pg, SIGHUP):
                return
            sleep(interval)
            if killpg(pg, SIGTERM):
                return
            sleep(interval)
            if killpg(pg, SIGKILL):
                return
