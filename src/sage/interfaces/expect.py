# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.libs.gap sage.libs.pari sage.libs.singular sage.symbolic
"""
Common Interface Functionality through Pexpect

See the examples in the other sections for how to use specific
interfaces. The interface classes all derive from the generic
interface that is described in this section.

AUTHORS:

- William Stein (2005): initial version

- William Stein (2006-03-01): got rid of infinite loop on startup if
  client system missing

- Felix Lawrence (2009-08-21): edited ._sage_() to support lists and float exponents in foreign notation.

- Simon King (2010-09-25): Expect._local_tmpfile() depends on
  Expect.pid() and is cached; Expect.quit() clears that cache,
  which is important for forking.

- Jean-Pierre Flori (2010,2011): Split non Pexpect stuff into a parent class.

- Simon King (2010-11-23): Ensure that the interface is started again
  after a crash, when a command is executed in _eval_line. Allow
  synchronisation of the GAP interface.

- François Bissey, Bill Page, Jeroen Demeyer (2015-12-09): Upgrade to
  pexpect 4.0.1 + patches, see :issue:`10295`.
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import io
import os
import re
import shlex
import signal
import sys
import weakref
import time
import gc
from . import quit
from random import randrange

import pexpect
from pexpect import ExceptionPexpect
import sage.interfaces.abc
from sage.interfaces.sagespawn import SageSpawn
from sage.interfaces.interface import (Interface, InterfaceElement,
            InterfaceFunction, InterfaceFunctionElement)

from sage.structure.element import RingElement

from sage.env import SAGE_EXTCODE, LOCAL_IDENTIFIER
from sage.misc.object_multiplexer import Multiplex
from sage.misc.instancedoc import instancedoc

from sage.cpython.string import str_to_bytes, bytes_to_str

BAD_SESSION = -2

# The subprocess is a shared resource.  In a multi-threaded
# environment, there would have to be a lock to control access to the
# subprocess.  Fortunately, Sage does not use Python threads.
# Unfortunately, the combination of the garbage collector and __del__
# methods gives rise to the same issues.  So, in places where we need
# to do a sequence of operations on the subprocess and make sure
# nothing else intervenes (for example, when we write a command and
# then read back the result) we need to disable the garbage collector.
# See TRAC #955 for a more detailed description of this problem.

# To turn off the garbage collector for a particular region of code,
# do:
#   with gc_disabled():
#       ... your code goes here ...
# The garbage collector will be returned to its original state
# whenever the code exits by any means (falling off the end, executing
# "return", "break", or "continue", raising an exception, ...)


class gc_disabled():
    """
    This is a "with" statement context manager. Garbage collection is
    disabled within its scope. Nested usage is properly handled.

    EXAMPLES::

        sage: import gc
        sage: from sage.interfaces.expect import gc_disabled
        sage: gc.isenabled()
        True
        sage: with gc_disabled():
        ....:     print(gc.isenabled())
        ....:     with gc_disabled():
        ....:         print(gc.isenabled())
        ....:     print(gc.isenabled())
        False
        False
        False
        sage: gc.isenabled()
        True
    """
    def __enter__(self):
        self._enabled = gc.isenabled()
        gc.disable()

    def __exit__(self, ty, val, tb):
        if self._enabled:
            gc.enable()
        return False


class Expect(Interface):
    """
    Expect interface object.
    """
    def __init__(self, name, prompt, command=None, env={}, server=None,
                 server_tmpdir=None,
                 ulimit=None, maxread=None,
                 script_subdirectory=None, restart_on_ctrlc=False,
                 verbose_start=False, init_code=[], max_startup_time=None,
                 logfile=None, eval_using_file_cutoff=0,
                 do_cleaner=True, remote_cleaner=False, path=None,
                 terminal_echo=True):

        Interface.__init__(self, name)

        # Read environment variables
        env_name = 'SAGE_%s_{}' % self.name().upper()
        if server is None:
            server = os.getenv(env_name.format('SERVER'))
        if server_tmpdir is None:
            server_tmpdir = os.getenv(env_name.format('TMPDIR'))
        if script_subdirectory is None:
            script_subdirectory = os.getenv(env_name.format('SCRIPT_SUBDIRECTORY'))
        self.__is_remote = False
        self.__remote_ulimit = None
        self.__remote_cleaner = remote_cleaner
        self._expect = None
        self._eval_using_file_cutoff = eval_using_file_cutoff
        self.__verbose_start = verbose_start
        self.set_server_and_command(server, command, server_tmpdir, ulimit)
        self._env = env
        self.__do_cleaner = do_cleaner
        self._change_prompt(prompt)
        self._restart_on_ctrlc = restart_on_ctrlc
        if path is not None:
            self.__path = os.path.abspath(path)
        elif script_subdirectory is None:
            self.__path = os.getcwd()
        else:
            self.__path = os.path.join(SAGE_EXTCODE, name, script_subdirectory)
        if not os.path.isdir(self.__path):
            raise OSError("path %r does not exist" % self.__path)
        self.__initialized = False
        self.__seq = -1
        self._session_number = 0
        self.__init_code = init_code

        # Handle the log file
        if isinstance(logfile, str):
            self.__logfile = None
            self.__logfilename = logfile
        else:
            self.__logfile = logfile
            self.__logfilename = None

        quit.expect_objects.append(weakref.ref(self))
        self._available_vars = []
        self._terminal_echo = terminal_echo

    def set_server_and_command(self, server=None, command=None, server_tmpdir=None, ulimit=None):
        """
        Changes the server and the command to use for this interface.

        This raises a :exc:`RuntimeError` if the interface is already started.

        INPUT:

        - ``server`` -- string or ``None`` (default); name of a remote host to connect to using ``ssh``

        - ``command`` -- one of:

          - a string; command line passed to the shell

          - a sequence of an :class:`~sage.features.Executable` and strings, arguments to
            pass to the executable.

        EXAMPLES::

            sage: magma.set_server_and_command(server='remote', command='mymagma')  # indirect doctest
            No remote temporary directory (option server_tmpdir) specified, using /tmp/ on remote
            sage: magma.server()
            'remote'
            sage: magma.command()
            'ssh -t remote mymagma'
        """
        if self._expect:
            raise RuntimeError("interface has already started")
        self._server = server
        self.__remote_ulimit = ulimit
        if server is not None:
            self.__is_remote = True
            self._eval_using_file_cutoff = 0  # don't allow this!
            if self.__verbose_start:
                print("Using remote server")
            if server_tmpdir is None:
                # TO DO: Why default to /tmp/? Might be better to use the expect process itself to get a tmp folder
                print("No remote temporary directory (option server_tmpdir) specified, using /tmp/ on " + server)
                self.__remote_tmpdir = "/tmp/"
            else:
                self.__remote_tmpdir = server_tmpdir
        else:
            self.__is_remote = False
        self.__command = command

    def server(self):
        """
        Return the server used in this interface.

        EXAMPLES::

            sage: magma.set_server_and_command(server='remote')
            No remote temporary directory (option server_tmpdir) specified, using /tmp/ on remote
            sage: magma.server() # indirect doctest
            'remote'
        """
        return self._server

    def command(self):
        """
        Return the command used in this interface as a string.

        EXAMPLES::

            sage: magma.set_server_and_command(command='magma-2.19')
            sage: magma.command()  # indirect doctest
            'magma-2.19'
        """
        command = self.__command
        server = self.server()
        if command is None:
            env_name = 'SAGE_%s_{}' % self.name().upper()  # same as in __init__
            command = os.getenv(env_name.format('COMMAND'), self.name())
        elif not isinstance(command, str):
            executable = command[0]
            if server:
                executable = executable.name
            else:
                executable = executable.absolute_filename()
            command = ' '.join([shlex.quote(executable)]
                               + [shlex.quote(arg) for arg in command[1:]])
        if server:
            if self.__remote_ulimit:
                command = f"ulimit {self.__remote_ulimit}; {command}"
            command = f"ssh -t {shlex.quote(server)} {shlex.quote(command)}"

        return command

    def _get(self, wait=0.1, alternate_prompt=None):
        if self._expect is None:
            self._start()
        E = self._expect
        wait = float(wait)
        try:
            if alternate_prompt is None:
                E.expect(self._prompt, timeout=wait)
            else:
                E.expect(str_to_bytes(alternate_prompt), timeout=wait)
        except pexpect.TIMEOUT:
            # TODO: In case an unexpected error occurred it's possible that the
            # contents of E.before, if consisting of multi-byte encoded text,
            # may be incomplete and contain errors, so merely calling
            # bytes_to_str here is probably not sufficient.
            return False, self._before()
        except pexpect.EOF:
            return True, self._before()
        except Exception:   # weird major problem!
            return True, self._before()
        return True, self._before()

    def _send(self, cmd):
        if self._expect is None:
            self._start()
        E = self._expect
        self.__so_far = ''
        E.sendline(cmd)

    def is_running(self):
        """
        Return ``True`` if ``self`` is currently running.
        """
        if self._expect is None:
            return False
        try:
            os.kill(self._expect.pid, 0)
        except OSError:
            # This means the process is not running
            return False
        return True

    def _so_far(self, wait=0.1, alternate_prompt=None):
        """
        Return whether done and output so far and new output since last
        time called.
        """
        done, new = self._get(wait=wait, alternate_prompt=alternate_prompt)
        try:
            if done:
                # if new is not None:
                X = self.__so_far + new
                del self.__so_far
                return True, X, new
            # new = self._expect.before
            try:
                self.__so_far += new
            except (AttributeError, TypeError):
                self.__so_far = new
            return False, self.__so_far, new
        except AttributeError as msg:   # no __so_far
            raise RuntimeError(msg)

    def is_remote(self):
        return self.__is_remote

    def is_local(self):
        return not self.__is_remote

    def user_dir(self):
        return self.__path

    def _change_prompt(self, prompt):
        if isinstance(prompt, str):
            prompt = str_to_bytes(prompt)
        elif (isinstance(prompt, type(re.compile(''))) and
                isinstance(prompt.pattern, str)):
            prompt = re.compile(str_to_bytes(prompt.pattern),
                                prompt.flags & ~re.U)
        self._prompt = prompt

    def path(self):
        return self.__path

    def expect(self):
        if self._expect is None:
            self._start()
        return self._expect

    def pid(self):
        """
        Return the PID of the underlying sub-process.

        REMARK:

        If the interface terminates unexpectedly, the original
        PID will still be used. But if it was terminated using
        :meth:`quit`, a new sub-process with a new PID is
        automatically started.

        EXAMPLES::

            sage: pid = gap.pid()
            sage: gap.eval('quit;')
            ''
            sage: pid == gap.pid()
            True
            sage: gap.quit()
            sage: pid == gap.pid()
            False
        """
        if self._expect is None:
            self._start()
        return self._expect.pid

    def _install_hints(self):
        r"""
        Hints for installing needed slave program on your computer.

        There are no hints by default.
        """
        return ''

    def _install_hints_ssh(self):
        r"""
        Hints for installing passwordless authentication on your
        computer...
        """
        # Written by Paul-Olivier Dehaye 2007/08/23
        return """
In order for Sage (on "local") to launch a "slave" process on "remote", the following command needs to work from local's console, without the need to enter any password:

       "ssh -t remote slave",

where "slave" could be "math" (for text-mode Mathematica), "gap", "magma", "sage", "maple", etc.

This thus requires passwordless authentication to be setup, which can be done with commands like these:
        cd; ssh-keygen -t rsa; scp .ssh/id_rsa.pub remote:.ssh/authorized_keys2\n
(WARNING: this would overwrite your current list of authorized keys on "remote")

In many cases, the server that can actually run "slave" is not accessible from the internet directly, but you have to hop through an intermediate trusted server, say "gate".
If that is your case, get help with _install_hints_ssh_through_gate().
"""

    def _install_hints_ssh_through_gate(self):
        r"""
        Hints for installing passwordless authentication through a gate
        """
        # Written by Paul-Olivier Dehaye 2007/08/23
        return """

 We assume you would like to run a "slave" process  on a machine called "remote" from a machine running Sage called "local". We also assume "remote" can only be accessed from "local" by ssh'ing first to "gate" (this is a fairly common setup). Sometimes, "gate" and "remote" have a shared filesystem, and this helps a bit.

  Note: You cannot just create shell scripts on "local" and "gate" that would use two successive SSH connections to "remote" in order to simulate running "slave" locally. This is because Sage will sometimes use files (and scp)  to communicate with "remote", which shell scripts would not take care of.

You need to setup:
 * passwordless authentication to "gate" from "local"
 * add passwordless authentication to "remote" from "local",
   for instance by appending the file local:~/.ssh/id_rsa.pub to remote:~/.ssh/authorized_keys2 and logging in once
      (this is only needed if "remote" and "gate" don\'t share filesystem)
 * add a few lines to your local:~/.ssh/ssh_config. Mine look like

       Host remote_for_sage
            ProxyCommand ssh gate nc -w 1 remote 22

That's it, normally.

The last step tells ssh that whenever an ssh connection is required to
the host "remote_for_sage", it should tunnel it through "gate". Any
attempt to scp-connect to "remote_for_sage" will use ssh and thus
this configuration file, and properly channel those file transfers
through the tunnel.

A good test is to attempt an scp connection from the command-line
of "local" to "remote_for_sage" as if no tunnel through "gate" was
required. No password should be asked for the second time around.

Finally, we created the new name "remote_for_sage" for "remote",
but this name only exists locally. this is to avoid interfering
with any other program that might already ssh to "remote" in
their own way.

If this all works, you can then make calls like:
         math = Mathematica(server='remote_for_sage')
"""

    def _do_cleaner(self):
        try:
            return self.__do_cleaner
        except AttributeError:
            return False

    def _start(self, alt_message=None, block_during_init=True):
        if self.is_running():
            # In case one is already running. We check first because
            # quit() can reset the local temporary filename at an
            # unexpected time as the process is started "on demand."
            self.quit()

        self._session_number += 1

        if self.__logfile is None:
            # If the 'SAGE_PEXPECT_LOG' environment variable is set and
            # there is no logfile already defined, then create a
            # logfile in .sage/pexpect_logs/
            if self.__logfilename is None and 'SAGE_PEXPECT_LOG' in os.environ:
                from sage.env import DOT_SAGE
                logs = os.path.join(DOT_SAGE, 'pexpect_logs')
                os.makedirs(logs, exist_ok=True)

                filename = '{name}-{pid}-{id}-{session}'.format(
                    name=self.name(), pid=os.getpid(), id=id(self),
                    session=self._session_number)
                self.__logfilename = os.path.join(logs, filename)
            if self.__logfilename is not None:
                self.__logfile = open(self.__logfilename, 'wb')

        cmd = self.command()

        if self.__verbose_start:
            print(cmd)
            print("Starting %s" % cmd.split()[0])

        if self.__remote_cleaner and self._server:
            c = 'ssh %s "nohup sage -cleaner"  &' % self._server
            os.system(c)

        # Unset some environment variables for the children to
        # reduce the chances they do something complicated breaking
        # the terminal interface.
        # See Issue #12221 and #13859.
        pexpect_env = dict(os.environ)
        pexpect_env.update(self._env)
        pexpect_del_vars = ['TERM', 'COLUMNS']
        for i in pexpect_del_vars:
            try:
                del pexpect_env[i]
            except KeyError:
                pass

        # Run child from self.__path
        currentdir = os.getcwd()
        os.chdir(self.__path)
        try:
            try:
                self._expect = SageSpawn(cmd,
                        logfile=self.__logfile,
                        timeout=None,  # no timeout
                        env=pexpect_env,
                        name=self._repr_(),
                        echo=self._terminal_echo,
                        # Work around https://bugs.python.org/issue1652
                        preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL),
                        quit_string=self._quit_string())

                # Attempt to shutdown the running process gracefully
                # when sage terminates.
                import atexit
                atexit.register(self.quit)

            except (ExceptionPexpect, pexpect.EOF) as e:
                # Change pexpect errors to RuntimeError
                raise RuntimeError("unable to start %s because the command %r failed: %s\n%s" %
                        (self.name(), cmd, e, self._install_hints()))
        except BaseException:
            self._expect = None
            self._session_number = BAD_SESSION
            raise
        finally:
            os.chdir(currentdir)

        if self._do_cleaner():
            quit.register_spawned_process(self._expect.pid, cmd)

        try:
            self._expect.expect(self._prompt)
        except (pexpect.TIMEOUT, pexpect.EOF) as msg:
            self._expect = None
            self._session_number = BAD_SESSION
            raise RuntimeError("unable to start %s: %s" % (self.name(), msg))
        self._expect.timeout = None

        with gc_disabled():
            if block_during_init:
                for X in self.__init_code:
                    self.eval(X)
            else:
                for X in self.__init_code:
                    self._send(X)

    def _isalive(self):
        """
        Wrapper for pexpect's ``spawn.isalive()``.

        Handles an issue where if the underlying process disappear (died / was
        killed and wait()-ed by another process) before pexpect itself could
        wait() on it, then pexpect (really ptyprocess) raises an exception but
        does *not* mark the process as terminated.  The same exception results,
        then, from any attempt to close the pexpect process.

        See https://github.com/sagemath/sage/issues/28354
        """
        try:
            return self._expect is not None and self._expect.isalive()
        except ExceptionPexpect:
            self._expect.ptyproc.terminated = True
            return False

    def _close(self, force=True):
        """
        Wrapper for pexpect's ``spawn.close()``.

        Since the underlying method calls ``isalive()`` it is affected by the
        same issue described in ``_isalive()`` above.
        """
        try:
            if self._expect is not None:
                self._expect.close(force=force)
        except (ExceptionPexpect, OSError):
            self._expect.ptyproc.fd = -1
            self._expect.ptyproc.closed = True
            self._expect.child_fd = -1
            self._expect.closed = True

    def clear_prompts(self):
        while True:
            try:
                self._expect.expect(self._prompt, timeout=0.1)
            except pexpect.TIMEOUT:
                return

    def _reset_expect(self):
        """
        Delete ``self._expect`` and reset any state.

        This is called by :meth:`quit` and :meth:`detach`.

        EXAMPLES::

            sage: gp("eulerphi(49)")
            42
            sage: print(gp._expect)
            PARI/GP interpreter with PID ...
            sage: gp._reset_expect()
            sage: print(gp._expect)
            None
            sage: gp("eulerphi(49)")
            42
        """
        self._session_number += 1
        try:
            # Spaghetti alert: when running several computations in
            # parallel, the pexpect interface is reset in each one,
            # and this next line is needed to trigger the generation
            # of a new temporary file when otherwise the existing
            # member variable would be shared. That also means that
            # you can't use this method to clean up an existing
            # tmpfile, because you can delete the one that the parent
            # is using.
            del self.__local_tmpfile
        except AttributeError:
            pass
        self._expect = None

    def quit(self, verbose=False):
        """
        Quit the running subprocess.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``); whether to print a
          message when quitting the process

        EXAMPLES::

            sage: a = maxima('y')
            sage: maxima.quit(verbose=True)
            Exiting Maxima with PID ... running ...maxima...
            sage: a._check_valid()
            Traceback (most recent call last):
            ...
            ValueError: The maxima session in which this object was defined is no longer running.

        Calling ``quit()`` a second time does nothing::

            sage: maxima.quit(verbose=True)
        """
        if self._expect is not None:
            if verbose:
                if self.is_remote():
                    print("Exiting %r (running on %s)" % (self._expect, self._server))
                else:
                    print("Exiting %r" % (self._expect,))
            self._close()
        self._reset_expect()

    def detach(self):
        """
        Forget the running subprocess: keep it running but pretend that
        it's no longer running.

        EXAMPLES::

            sage: a = maxima('y')
            sage: saved_expect = maxima._expect  # Save this to close later
            sage: maxima.detach()
            sage: a._check_valid()
            Traceback (most recent call last):
            ...
            ValueError: The maxima session in which this object was defined is no longer running.
            sage: saved_expect.close()  # Close child process

        Calling ``detach()`` a second time does nothing::

            sage: maxima.detach()
        """
        try:
            self._expect._keep_alive()
        except AttributeError:
            pass
        self._reset_expect()

    def _quit_string(self):
        """
        Return the string which will be used to quit the application.

        EXAMPLES::

            sage: gp._quit_string()
            '\\q'
            sage: maxima._quit_string()
            'quit();'
        """
        return 'quit'

    def _send_interrupt(self):
        """
        Send an interrupt to the application.  This is used internally
        by :meth:`interrupt`.

        First CTRL-C to stop the current command, then quit.
        """
        self._expect.sendline(chr(3))
        self._expect.sendline(self._quit_string())

    def _local_tmpfile(self):
        """
        Return a filename that is used to buffer long command lines for this interface.

        INPUT:

        - ``e`` -- an expect interface instance

        OUTPUT:

        A string that provides a temporary filename and is unique for the
        given interface.

        TESTS:

        The filename is cached::

            sage: gap._local_tmpfile() is gap._local_tmpfile()
            True

        The following two problems were fixed in :issue:`10004`.

        1. Different interfaces have different temp-files::

            sage: gap._local_tmpfile() != singular._local_tmpfile()
            True

        2. Interface instances in different branches of a parallelised
           function have different temp-files::

            sage: @parallel
            ....: def f(n):
            ....:     return gap._local_tmpfile()
            sage: L = [t[1] for t in f(list(range(5)))]
            sage: len(set(L))
            5

        The following used to fail::

            sage: s = gap._local_tmpfile()
            sage: L = [t[1] for t in f(list(range(5)))]
            sage: len(set(L))
            5

        AUTHOR:

        - Simon King (2010-09): Making the tmp-file unique for the interface instance
        """
        try:
            return self.__local_tmpfile
        except AttributeError:
            pass

        import atexit
        from tempfile import NamedTemporaryFile
        # FriCAS uses the ".input" suffix, and the other
        # interfaces are suffix-agnostic, so using ".input" here
        # lets us avoid a subclass override for FriCAS.
        with NamedTemporaryFile(suffix='.input', delete=False) as f:
            self.__local_tmpfile = f.name
            atexit.register(lambda: os.remove(f.name))
        return self.__local_tmpfile

    def _remote_tmpdir(self):
        return self.__remote_tmpdir

    def _remote_tmpfile(self):
        try:
            return self.__remote_tmpfile
        except AttributeError:
            self.__remote_tmpfile = self._remote_tmpdir() + "/interface_%s:%s" % (LOCAL_IDENTIFIER, self.pid())
            return self.__remote_tmpfile

    def _send_tmpfile_to_server(self, local_file=None, remote_file=None):
        if local_file is None:
            local_file = self._local_tmpfile()
        if remote_file is None:
            remote_file = self._remote_tmpfile()
        cmd = 'scp "%s" %s:"%s" 1>&2 2>/dev/null' % (local_file, self._server, remote_file)
        os.system(cmd)

    def _get_tmpfile_from_server(self, local_file=None, remote_file=None):
        if local_file is None:
            local_file = self._local_tmpfile()
        if remote_file is None:
            remote_file = self._remote_tmpfile()
        cmd = 'scp %s:"%s" "%s" 1>&2 2>/dev/null' % (self._server, remote_file, local_file)
        os.system(cmd)

    def _remove_tmpfile_from_server(self):
        if self.__remote_tmpfile is not None:
            raise NotImplementedError

    def _eval_line_using_file(self, line, restart_if_needed=True):
        """
        Evaluate a line of commands, using a temporary file.

        REMARK:

        By default, this is called when a long command is
        evaluated in :meth:`eval`.

        If the command can not be evaluated since the interface
        has crashed, it is automatically restarted and tried
        again *once*.

        INPUT:

        - ``line`` -- string; a command
        - ``restart_if_needed`` -- boolean (default: ``True``);
          if it is ``True``, the command evaluation is evaluated
          a second time after restarting the interface, if an
          :exc:`EOFError` occurred.

        TESTS::

            sage: singular._eval_line_using_file('def a=3;')
            ''
            sage: singular('a')
            3
            sage: singular.eval('quit;')
            ''
            sage: singular._eval_line_using_file('def a=3;')
            Singular crashed -- automatically restarting.
            ''
            sage: singular('a')
            3
            sage: singular.eval('quit;')
            ''
            sage: singular._eval_line_using_file('def a=3;', restart_if_needed=False)
            Traceback (most recent call last):
            ...
            RuntimeError: Singular terminated unexpectedly while reading in a large line...

        We end by triggering a re-start of Singular, since otherwise
        the doc test of another method would fail by a side effect.
        ::

            sage: singular(3)
            Singular crashed -- automatically restarting.
            3
        """
        with open(self._local_tmpfile(), 'w') as F:
            F.write(line + '\n')

        tmp_to_use = self._local_tmpfile()
        if self.is_remote():
            self._send_tmpfile_to_server()
            tmp_to_use = self._remote_tmpfile()
        try:
            s = self._eval_line(self._read_in_file_command(tmp_to_use), allow_use_file=False, restart_if_needed=False)
        except pexpect.EOF:
            if self._quit_string() in line:
                # we expect to get an EOF if we're quitting.
                return ''
            elif restart_if_needed:  # the subprocess might have crashed
                try:
                    self._synchronize()
                    return self._post_process_from_file(self._eval_line_using_file(line, restart_if_needed=False))
                except RuntimeError as msg:
                    raise RuntimeError('%s terminated unexpectedly while reading in a large line:\n%s' % (self, msg.args[0]))
                except TypeError:
                    pass
            raise RuntimeError('%s terminated unexpectedly while reading in a large line' % self)
        except RuntimeError as msg:
            if self._quit_string() in line:
                if not self._isalive():
                    return ''
                raise
            if restart_if_needed and not self._isalive():
                try:
                    self._synchronize()
                    return self._post_process_from_file(self._eval_line_using_file(line, restart_if_needed=False))
                except TypeError:
                    pass
                except RuntimeError:
                    raise RuntimeError('%s terminated unexpectedly while reading in a large line' % self)
            if "Input/output error" in msg.args[0]:
                # This occurs on non-linux machines
                raise RuntimeError('%s terminated unexpectedly while reading in a large line' % self)
            raise RuntimeError('%s terminated unexpectedly while reading in a large line:\n%s' % (self, msg.args[0]))
        return self._post_process_from_file(s)

    def _post_process_from_file(self, s):
        return s

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True, restart_if_needed=True):
        """
        Evaluate a line of commands.

        REMARK:

        By default, a long command (length exceeding ``self._eval_using_file_cutoff``)
        is evaluated using :meth:`_eval_line_using_file`.

        If the command can not be evaluated since the interface
        has crashed, it is automatically restarted and tried
        again *once*.

        If the optional ``wait_for_prompt`` is ``False`` then even a very
        long line will not be evaluated by :meth:`_eval_line_using_file`,
        since this does not support the ``wait_for_prompt`` option.

        INPUT:

        - ``line`` -- string; a command
        - ``allow_use_file`` -- boolean (default: ``True``);
          allow to evaluate long commands using :meth:`_eval_line_using_file`.
        - ``wait_for_prompt`` -- boolean (default: ``True``);
          wait until the prompt appears in the sub-process' output.
        - ``restart_if_needed`` -- boolean (default: ``True``);
          if it is ``True``, the command evaluation is evaluated
          a second time after restarting the interface, if an
          :exc:`EOFError` occurred.

        TESTS::

            sage: from sage.interfaces.singular import singular
            sage: singular._eval_line('def a=3;')
            ''
            sage: singular('a')
            3
            sage: singular.eval('quit;')
            ''
            sage: singular._eval_line('def a=3;')
            Singular crashed -- automatically restarting.
            ''
            sage: singular('a')
            3
            sage: singular.eval('kill a')
            ''

        We are now sending a command that would run forever. But since
        we declare that we are not waiting for a prompt, we can interrupt
        it without a KeyboardInterrupt. At the same time, we test that
        the line is not forwarded to :meth:`_eval_line_using_file`, since
        that method would not support the ``wait_for_prompt`` option.
        For reasons which are currently not understood, the ``interrupt``
        test usually returns immediately, but sometimes it takes a very
        long time on the same system. ::

            sage: cutoff = singular._eval_using_file_cutoff
            sage: singular._eval_using_file_cutoff = 4
            sage: singular._eval_line('for(int i=1;i<=3;i++){i=1;};', wait_for_prompt=False)
            ''
            sage: singular.interrupt()
            True
            sage: singular._eval_using_file_cutoff = cutoff

        The interface still works after this interrupt::

            sage: singular('2+3')
            5

        Last, we demonstrate that by default the execution of a command
        is tried twice if it fails the first time due to a crashed
        interface::

            sage: singular.eval('quit;')
            ''
            sage: singular._eval_line_using_file('def a=3;', restart_if_needed=False)
            Traceback (most recent call last):
            ...
            RuntimeError: Singular terminated unexpectedly while reading in a large line...

        Since the test of the next method would fail, we re-start
        Singular now. ::

            sage: singular('2+3')
            Singular crashed -- automatically restarting.
            5
        """
        if allow_use_file and wait_for_prompt and self._eval_using_file_cutoff and len(line) > self._eval_using_file_cutoff:
            return self._eval_line_using_file(line)
        try:
            if self._expect is None:
                self._start()
            E = self._expect
            try:
                if len(line) >= 4096:
                    raise RuntimeError("Sending more than 4096 characters with %s on a line may cause a hang and you're sending %s characters" % (self, len(line)))
                E.sendline(line)
                if not wait_for_prompt:
                    return ''

            except OSError as msg:
                if restart_if_needed:
                    # The subprocess most likely crashed.
                    # If it's really still alive, we fall through
                    # and raise RuntimeError.
                    if sys.platform.startswith('sunos'):
                        # On (Open)Solaris, we might need to wait a
                        # while because the process might not die
                        # immediately. See Issue #14371.
                        for t in [0.5, 1.0, 2.0]:
                            if self._isalive():
                                time.sleep(t)
                            else:
                                break
                    if not self._isalive():
                        try:
                            self._synchronize()
                        except (TypeError, RuntimeError):
                            pass
                        return self._eval_line(line, allow_use_file=allow_use_file, wait_for_prompt=wait_for_prompt, restart_if_needed=False)
                raise RuntimeError("%s\nError evaluating %s in %s" % (msg, line, self))

            if line:
                try:
                    if isinstance(wait_for_prompt, str):
                        E.expect(str_to_bytes(wait_for_prompt))
                    else:
                        E.expect(self._prompt)
                except pexpect.EOF as msg:
                    try:
                        if self.is_local():
                            tmp_to_use = self._local_tmpfile()
                        else:
                            tmp_to_use = self._remote_tmpfile()
                        if self._read_in_file_command(tmp_to_use) in line:
                            raise pexpect.EOF(msg)
                    except NotImplementedError:
                        pass
                    if self._quit_string() in line:
                        # we expect to get an EOF if we're quitting.
                        return ''
                    elif restart_if_needed:  # the subprocess might have crashed
                        try:
                            self._synchronize()
                            return self._eval_line(line, allow_use_file=allow_use_file, wait_for_prompt=wait_for_prompt, restart_if_needed=False)
                        except (TypeError, RuntimeError):
                            pass
                    raise RuntimeError("%s\n%s crashed executing %s" % (msg, self, line))
                if self._terminal_echo:
                    out = self._before()
                else:
                    out = self._before().rstrip('\n\r')
            else:
                if self._terminal_echo:
                    out = '\n\r'
                else:
                    out = ''
        except KeyboardInterrupt:
            self._keyboard_interrupt()
            raise KeyboardInterrupt("Ctrl-c pressed while running %s" % self)
        if self._terminal_echo:
            i = out.find("\n")
            j = out.rfind("\r")
            return out[i + 1:j].replace('\r\n', '\n')
        else:
            return out.replace('\r\n', '\n')

    def _keyboard_interrupt(self):
        print("Interrupting %s..." % self)
        if self._restart_on_ctrlc:
            try:
                self._close()
            except pexpect.ExceptionPexpect as msg:
                raise pexpect.ExceptionPexpect("THIS IS A BUG -- PLEASE REPORT. This should never happen.\n" + msg)
            self._start()
            raise KeyboardInterrupt("Restarting %s (WARNING: all variables defined in previous session are now invalid)" % self)
        else:
            self._expect.sendline(chr(3))  # send ctrl-c
            self._expect.expect(self._prompt)
            self._expect.expect(self._prompt)
            raise KeyboardInterrupt("Ctrl-c pressed while running %s" % self)

    def interrupt(self, tries=5, timeout=2.0, quit_on_fail=True):
        E = self._expect
        if E is None:
            return True
        try:
            for i in range(tries):
                self._send_interrupt()
                try:
                    E.expect(self._prompt, timeout=timeout)
                except (pexpect.TIMEOUT, pexpect.EOF):
                    pass
                else:
                    return True  # Success
        except Exception:
            pass
        # Failed to interrupt...
        if quit_on_fail:
            self.quit()
        return False

    ###########################################################################
    # BEGIN Synchronization code.
    ###########################################################################

    def _before(self, encoding=None, errors=None) -> str:
        r"""
        Return the previous string that was sent through the interface.

        This returns a ``str`` object.

        The ``encoding`` and ``errors`` arguments are passed to
        :func:`sage.misc.cpython.bytes_to_str`.

        EXAMPLES::

            sage: singular('2+3')
            5
            sage: singular._before()
            '5\r\n'
        """
        return bytes_to_str(self._expect.before, encoding, errors)

    def _after(self, encoding=None, errors=None):
        r"""
        Return trailing data in the buffer after the text matched by the expect
        interface.

        When the ``spawn.after`` attribute contains bytes, this returns ``str``
        objects on both Python 2 and Python 3.  There are also cases (such as
        exceptions) where the ``.after`` attribute contains either an exception
        type or ``None``, in which case those values are returned.

        The ``encoding`` and ``errors`` arguments are passed to
        :func:`sage.misc.cpython.bytes_to_str`.

        EXAMPLES::

            sage: singular('2+3')
            5
            sage: singular._after()
            '> '
        """
        after = self._expect.after
        if isinstance(after, bytes):
            return bytes_to_str(after, encoding, errors)

        return after

    def _readline(self, size=-1, encoding=None, errors=None) -> str:
        r"""
        Wraps ``spawn.readline`` to pass the return values through
        ``bytes_to_str``, like `Expect._before` and `Expect._after`.

        EXAMPLES::

            sage: a = singular(1)
            sage: singular._expect.sendline('1+1;')
            5
            sage: singular._readline()
            '2\r\n'
        """
        return bytes_to_str(self._expect.readline(size=size), encoding, errors)

    def _interrupt(self):
        for i in range(15):
            try:
                self._sendstr('quit;\n' + chr(3))
                self._expect_expr(timeout=2)
            except pexpect.TIMEOUT:
                pass
            except pexpect.EOF:
                self._crash_msg()
                self.quit()
            else:
                return

    def _expect_expr(self, expr=None, timeout=None):
        r"""
        Wait for a given expression expr (which could be a regular
        expression or list of regular expressions) to appear in the output
        for at most timeout seconds.

        Use ``r._expect.before`` to see what was put in the
        output stream before the expression.

        INPUT:

        - ``expr`` -- ``None`` or a string or list of strings
          (default: ``None``)

        - ``timeout`` -- ``None`` or a number (default: ``None``)

        EXAMPLES:

        We test all of this using the Singular interface. First we put
        10 + 15 in the input stream::

            sage: singular._sendstr('def abc = 10 + 15;\n')

        Then we tell singular to print 10, which is an arbitrary number
        different from the expected result 35::

            sage: singular._sendstr('10;\n')

        Here an exception is raised because 25 hasn't appeared yet in the
        output stream. The key thing is that this doesn't lock, but instead
        quickly raises an exception::

            sage: t = walltime()
            sage: try:
            ....:    singular._expect_expr('25', timeout=float(0.4))
            ....: except Exception:
            ....:    print('Did not get expression')
            Did not get expression

        A quick consistency check on the time that the above took::

            sage: w = walltime(t); 0.3 < w < 10
            True

        We tell Singular to print abc, which equals 25::

            sage: singular._sendstr('abc;\n')

        Now 25 is in the output stream, so we can wait for it::

            sage: singular._expect_expr('25')

        This gives us everything before the 25, including the 10 we printed earlier::

            sage: singular._expect.before.decode('ascii')
            '...10\r\n> '

        We test interrupting ``_expect_expr`` using the GP interface,
        see :issue:`6661`.  Unfortunately, this test doesn't work reliably using
        Singular, see :issue:`9163` and the follow-up :issue:`10476`.
        The ``gp.eval('0')`` in this test makes sure that ``gp`` is
        running, so a timeout of 1 second should be sufficient. ::

            sage: print(sage0.eval("dummy=gp.eval('0'); alarm(1); gp._expect_expr('1')"))  # long time
            ...Interrupting PARI/GP interpreter. Please wait a few seconds...
            ...
            AlarmInterrupt:
        """
        if expr is None:
            # the following works around gap._prompt_wait not being defined
            expr = getattr(self, '_prompt_wait', None) or self._prompt
        if self._expect is None:
            self._start()
        try:
            if timeout:
                i = self._expect.expect(expr, timeout=timeout)
            else:
                i = self._expect.expect(expr)
            if i > 0:
                v = self._before()
                self.quit()
                raise ValueError("%s\nComputation failed due to a bug in %s -- NOTE: Had to restart." % (v, self))
        except KeyboardInterrupt:
            print("Control-C pressed. Interrupting %s. Please wait a few seconds..." % self)
            self.interrupt()
            raise

    def _sendstr(self, string):
        r"""
        Send a string to the pexpect interface, autorestarting the expect
        interface if anything goes wrong.

        INPUT:

        - ``string`` -- string

        EXAMPLES: We illustrate this function using the Singular interface::

            sage: singular._synchronize()
            sage: singular._sendstr('int i = 5;')
            sage: singular('i')
            5
        """
        if self._expect is None:
            self._start()
        try:
            os.write(self._expect.child_fd, str_to_bytes(string))
        except OSError:
            self._crash_msg()
            self.quit()
            self._sendstr(string)

    def _crash_msg(self):
        r"""
        Show a message if the interface crashed.

        EXAMPLES::

            sage: singular._crash_msg()
            Singular crashed -- automatically restarting.

        ::

            sage: singular('2+3')
            5
            sage: singular._sendstr('quit;\n')   # make it so that singular appears to die.
            sage: singular('2+3')
            Singular crashed -- automatically restarting.
            5
        """
        print("%s crashed -- automatically restarting." % self)

    def _synchronize(self, cmd='1+%s;\n'):
        """
        Synchronize pexpect interface.

        This put a random integer (plus one!) into the output stream, then
        waits for it, thus resynchronizing the stream. If the random
        integer doesn't appear within 1 second, the interface is sent
        interrupt signals.

        This way, even if you somehow left the interface in a busy state
        computing, calling _synchronize gets everything fixed.

        EXAMPLES: We observe nothing, just as it should be::

            sage: singular._synchronize()

        TESTS:

        This illustrates a synchronization bug being fixed (thanks
        to Simon King and David Joyner for tracking this down)::

            sage: R.<x> = QQ[]; f = x^3 + x + 1;  g = x^3 - x - 1; r = f.resultant(g); gap(ZZ); singular(R)
            Integers
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C
        """
        if self._expect is None:
            return
        rnd = randrange(2147483647)
        s = str(rnd + 1)
        cmd = cmd % rnd
        self._sendstr(cmd)
        try:
            self._expect_expr(timeout=0.5)
            if s not in self._before():
                self._expect_expr(s, timeout=0.5)
                self._expect_expr(timeout=0.5)
        except pexpect.TIMEOUT:
            self._interrupt()
        except pexpect.EOF:
            self._crash_msg()
            self.quit()

    ###########################################################################
    # END Synchronization code.
    ###########################################################################

    def eval(self, code, strip=True, synchronize=False, locals=None, allow_use_file=True,
             split_lines='nofile', **kwds):
        """
        INPUT:

        - ``code`` -- text to evaluate

        - ``strip`` -- boolean; whether to strip output prompts,
          etc. (ignored in the base class)

        - ``locals`` -- None (ignored); this is used for compatibility
          with the Sage notebook's generic system interface

        - ``allow_use_file`` -- boolean (default: ``True``); if ``True`` and
          ``code`` exceeds an interface-specific threshold then ``code`` will
          be communicated via a temporary file rather that the character-based
          interface. If ``False`` then the code will be communicated via the
          character interface.

        - ``split_lines`` -- Tri-state (default: ``'nofile'``); if "nofile"
          then ``code`` is sent line by line unless it gets communicated via a
          temporary file. If ``True`` then ``code`` is sent line by line, but
          some lines individually might be sent via temporary file. Depending
          on the interface, this may transform grammatical ``code`` into
          ungrammatical input. If ``False``, then the whole block of code is
          evaluated all at once.

        - ``**kwds`` -- all other arguments are passed onto the ``_eval_line``
          method. An often useful example is ``reformat=False``.
        """
        if synchronize:
            try:
                self._synchronize()
            except AttributeError:
                pass

        if strip:
            try:
                code = self._strip_prompt(code)
            except AttributeError:
                pass

        if not isinstance(code, str):
            raise TypeError('input code must be a string.')

        # Remove extra whitespace
        code = code.strip()

        try:
            with gc_disabled():
                if (split_lines == "nofile" and allow_use_file and
                        self._eval_using_file_cutoff and len(code) > self._eval_using_file_cutoff):
                    return self._eval_line_using_file(code)
                elif split_lines:
                    return '\n'.join(self._eval_line(L, allow_use_file=allow_use_file, **kwds)
                                     for L in code.split('\n') if L)
                else:
                    return self._eval_line(code, allow_use_file=allow_use_file, **kwds)
        # DO NOT CATCH KeyboardInterrupt, as it is being caught
        # by _eval_line
        # In particular, do NOT call self._keyboard_interrupt()
        except TypeError as s:
            raise TypeError('error evaluating "%s":\n%s' % (code, s))

    ############################################################
    #         Functions for working with variables.
    #  The first three must be overloaded by derived classes,
    #  and the definition depends a lot on the class.  But
    #  the functionality one gets from this is very nice.
    ############################################################

    def _object_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.expect import Expect
            sage: Expect._object_class(maxima)
            <class 'sage.interfaces.expect.ExpectElement'>
        """
        return ExpectElement

    def _function_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.expect import Expect
            sage: Expect._function_class(maxima)
            <class 'sage.interfaces.expect.ExpectFunction'>
        """
        return ExpectFunction

    def _function_element_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.expect import Expect
            sage: Expect._function_element_class(maxima)
            <class 'sage.interfaces.expect.FunctionElement'>
        """
        return FunctionElement


@instancedoc
class ExpectFunction(InterfaceFunction):
    """
    Expect function.
    """
    pass


@instancedoc
class FunctionElement(InterfaceFunctionElement):
    """
    Expect function element.
    """
    pass


def is_ExpectElement(x):
    """
    Return ``True`` if ``x`` is of type :class:`ExpectElement`.

    This function is deprecated; use :func:`isinstance`
    (of :class:`sage.interfaces.abc.ExpectElement`) instead.

    EXAMPLES::

        sage: from sage.interfaces.expect import is_ExpectElement
        sage: is_ExpectElement(2)
        doctest:...: DeprecationWarning: the function is_ExpectElement is deprecated; use isinstance(x, sage.interfaces.abc.ExpectElement) instead
        See https://github.com/sagemath/sage/issues/34804 for details.
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(34804, "the function is_ExpectElement is deprecated; use isinstance(x, sage.interfaces.abc.ExpectElement) instead")

    return isinstance(x, ExpectElement)


@instancedoc
class ExpectElement(InterfaceElement, sage.interfaces.abc.ExpectElement):
    """
    Expect element.
    """
    def __init__(self, parent, value, is_name=False, name=None):
        RingElement.__init__(self, parent)
        self._create = value
        if parent is None:
            return     # means "invalid element"
        # idea: Joe Wetherell -- try to find out if the output
        # is too long and if so get it using file, otherwise
        # don't.
        if isinstance(value, str) and parent._eval_using_file_cutoff and \
           parent._eval_using_file_cutoff < len(value):
            self._get_using_file = True

        if is_name:
            self._name = value
        else:
            try:
                self._name = parent._create(value, name=name)
            # Convert ValueError and RuntimeError to TypeError for
            # coercion to work properly.
            except (RuntimeError, ValueError) as x:
                self._session_number = -1
                raise TypeError(*x.args)
            except BaseException:
                self._session_number = -1
                raise
        self._session_number = parent._session_number

    def __hash__(self):
        """
        Return the hash of ``self``.

        This is a default implementation of hash
        which just takes the hash of the string of ``self``.
        """
        return hash('%s%s' % (self, self._session_number))

    def _check_valid(self):
        """
        Check that this object is valid, i.e., the session in which this
        object is defined is still running. This is relevant for
        interpreters that can't be interrupted via ctrl-C, hence get
        restarted.
        """
        try:
            P = self.parent()
            if P is None or P._session_number == BAD_SESSION or self._session_number == -1 or \
               P._session_number != self._session_number:
                raise ValueError("The %s session in which this object was defined is no longer running." % P.name())
        except AttributeError:
            raise ValueError("The session in which this object was defined is no longer running.")
        return P

    def __del__(self):
        try:
            self._check_valid()
        except ValueError:
            return
        try:
            if hasattr(self, '_name'):
                P = self.parent()
                if P is not None:
                    P.clear(self._name)

        except (RuntimeError, ExceptionPexpect):  # needed to avoid infinite loops in some rare cases
            pass

#    def _sage_repr(self):
# TO DO: this could use file transfers when self.is_remote()


class StdOutContext:
    """
    A context in which all communication between Sage and a subprocess
    interfaced via pexpect is printed to stdout.
    """
    def __init__(self, interface, silent=False, stdout=None):
        """
        Construct a new context in which all communication between Sage
        and a subprocess interfaced via pexpect is printed to stdout.

        INPUT:

        - ``interface`` -- the interface whose communication shall be dumped

        - ``silent`` -- if ``True`` this context does nothing

        - ``stdout`` -- (optional) parameter for alternative stdout device (default: ``None``)

        EXAMPLES::

            sage: from sage.interfaces.expect import StdOutContext
            sage: with StdOutContext(Gp()) as g:
            ....:     g('1+1')
            sage=...
        """
        self.interface = interface
        self.silent = silent
        self.stdout = stdout if stdout else sys.stdout

    def __enter__(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.expect import StdOutContext
            sage: with StdOutContext(singular):
            ....:     singular.eval('1+1')
            1+1;
            ...
        """
        if self.silent:
            return self.interface
        if self.interface._expect is None:
            self.interface._start()
        self._logfile_backup = self.interface._expect.logfile

        if isinstance(self.stdout, io.TextIOWrapper):
            stdout = self.stdout.buffer
        else:
            stdout = self.stdout

        if self.interface._expect.logfile:
            self.interface._expect.logfile = Multiplex(self.interface._expect.logfile, stdout)
        else:
            self.interface._expect.logfile = Multiplex(stdout)
        return self.interface

    def __exit__(self, typ, value, tb):
        r"""
        EXAMPLES::

            sage: from sage.interfaces.expect import StdOutContext
            sage: with StdOutContext(gap):
            ....:     gap('1+1')
            \$sage...
        """
        if self.silent:
            return
        self.interface._expect.logfile.flush()
        self.stdout.write("\n")
        self.interface._expect.logfile = self._logfile_backup
