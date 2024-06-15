r"""
Parallel Interface to the Sage interpreter

This is an expect interface to \emph{multiple} copy of the \sage
interpreter, which can all run simultaneous calculations.  A PSage
object does not work as well as the usual Sage object, but does have
the great property that when you construct an object in a PSage you
get back a prompt immediately.  All objects constructed for that
PSage print <<currently executing code>> until code execution
completes, when they print as normal.

\note{BUG -- currently non-idle PSage subprocesses do not stop when
\sage exits.  I would very much like to fix this but don't know how.}

EXAMPLES:

We illustrate how to factor 3 integers in parallel.
First start up 3 parallel Sage interfaces::

    sage: v = [PSage() for _ in range(3)]

Next, request factorization of one random integer in each copy. ::

    sage: w = [x('factor(2^%s-1)'% randint(250,310)) for x in v]  # long time (5s on sage.math, 2011)

Print the status::

    sage: w       # long time, random output (depends on timing)
    [3 * 11 * 31^2 * 311 * 11161 * 11471 * 73471 * 715827883 * 2147483647 * 4649919401 * 18158209813151 * 5947603221397891 * 29126056043168521,
     <<currently executing code>>,
     9623 * 68492481833 * 23579543011798993222850893929565870383844167873851502677311057483194673]

Note that at the point when we printed two of the factorizations had
finished but a third one hadn't.   A few seconds later all three have
finished::

    sage: w       # long time, random output
    [3 * 11 * 31^2 * 311 * 11161 * 11471 * 73471 * 715827883 * 2147483647 * 4649919401 * 18158209813151 * 5947603221397891 * 29126056043168521,
     23^2 * 47 * 89 * 178481 * 4103188409 * 199957736328435366769577 * 44667711762797798403039426178361,
     9623 * 68492481833 * 23579543011798993222850893929565870383844167873851502677311057483194673]
"""

import os
import time

from .sage0 import Sage, SageElement
from pexpect import ExceptionPexpect

number = 0


class PSage(Sage):
    def __init__(self, **kwds):
        if 'server' in kwds:
            raise NotImplementedError("PSage doesn't work on remote server yet.")
        Sage.__init__(self, **kwds)
        import sage.misc.misc
        T = sage.misc.temporary_file.tmp_dir('sage_smp')
        self.__tmp_dir = T
        self.__tmp = '%s/lock' % T
        self._unlock()
        self._unlock_code = "with open('%s', 'w') as f: f.write('__unlocked__')" % self.__tmp

        global number
        self._number = number
        number += 1

    def _repr_(self):
        """
        TESTS::

            sage: from sage.interfaces.psage import PSage
            sage: PSage()                                   # indirect doctest
            A running non-blocking (parallel) instance of Sage (number ...)
        """
        return 'A running non-blocking (parallel) instance of Sage (number %s)' % (self._number)

    def _unlock(self):
        self._locked = False
        with open(self.__tmp, 'w') as fobj:
            fobj.write('__unlocked__')

    def _lock(self):
        self._locked = True
        with open(self.__tmp, 'w') as fobj:
            fobj.write('__locked__')

    def _start(self):
        Sage._start(self)
        self.expect().timeout = 0.25
        self.expect().delaybeforesend = 0.01

    def is_locked(self) -> bool:
        try:
            with open(self.__tmp) as fobj:
                if fobj.read() != '__locked__':
                    return False
        except FileNotFoundError:
            # Directory may have already been deleted :issue:`30730`
            return False
        # looks like we are locked, but check health first
        try:
            self.expect().expect(self._prompt)
            self.expect().expect(self._prompt)
        except ExceptionPexpect:
            return False
        return True

    def __del__(self):
        """
        TESTS:

        Check that :issue:`29989` is fixed::

            sage: PSage().__del__()
        """
        try:
            files = os.listdir(self.__tmp_dir)
            for x in files:
                os.remove(os.path.join(self.__tmp_dir, x))
            os.removedirs(self.__tmp_dir)
        except OSError:
            pass

        if not (self._expect is None):
            cmd = 'kill -9 %s' % self._expect.pid
            os.system(cmd)

    def eval(self, x, strip=True, **kwds):
        """
        INPUT:

        - ``x`` -- code
        - ``strip`` --ignored
        """
        if self.is_locked():
            return "<<currently executing code>>"
        if self._locked:
            self._locked = False
            # self._expect.expect('__unlocked__')
            self.expect().send('\n')
            self.expect().expect(self._prompt)
            self.expect().expect(self._prompt)
        try:
            return Sage.eval(self, x, **kwds)
        except ExceptionPexpect:
            return "<<currently executing code>>"

    def get(self, var):
        """
        Get the value of the variable var.
        """
        try:
            return self.eval('print(%s)' % var)
        except ExceptionPexpect:
            return "<<currently executing code>>"

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s=%s' % (var, value)
        self._send_nowait(cmd)
        time.sleep(0.02)

    def _send_nowait(self, x):
        if x.find('\n') != -1:
            raise ValueError("x must not have any newlines")
        # Now we want the client Python process to execute two things.
        # The first is x and the second is c.  The effect of c
        # will be to unlock the lock.
        if self.is_locked():
            return "<<currently executing code>>"
        E = self.expect()
        self._lock()
        E.write(self.preparse(x) + '\n')
        try:
            E.expect(self._prompt)
        except ExceptionPexpect:
            pass
        E.write(self._unlock_code + '\n\n')

    def _object_class(self):
        return PSageElement


class PSageElement(SageElement):
    def is_locked(self):
        return self.parent().is_locked()
