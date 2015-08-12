"""
Interface to the ``pselect()`` system call

This module defines a class :class:`PSelecter` which can be used to
call the system call ``pselect()`` and which can also be used in a
``with`` statement to block given signals until
:meth:`PSelecter.pselect` is called.

AUTHORS:

- Jeroen Demeyer (2013-02-07): initial version (:trac:`14079`)

Waiting for subprocesses
------------------------

One possible use is to wait with a **timeout** until **any child process**
exits, as opposed to ``os.wait()`` which doesn't have a timeout or
``multiprocessing.Process.join()`` which waits for one specific process.

Since ``SIGCHLD`` is ignored by default, we first need to install a
signal handler for ``SIGCHLD``. It doesn't matter what it does, as long
as the signal isn't ignored::

    sage: import signal
    sage: def dummy_handler(sig, frame):
    ....:     pass
    ....:
    sage: _ = signal.signal(signal.SIGCHLD, dummy_handler)

We wait for a child created using the ``subprocess`` module::

    sage: from sage.ext.pselect import PSelecter
    sage: from subprocess import *
    sage: with PSelecter([signal.SIGCHLD]) as sel:
    ....:     p = Popen(["sleep", "1"])
    ....:     _ = sel.sleep()
    ....:
    sage: p.poll()  # p should be finished
    0

Now using the ``multiprocessing`` module::

    sage: from sage.ext.pselect import PSelecter
    sage: from multiprocessing import *
    sage: import time
    sage: with PSelecter([signal.SIGCHLD]) as sel:
    ....:     p = Process(target=time.sleep, args=(1,))
    ....:     p.start()
    ....:     _ = sel.sleep()
    ....:     p.is_alive()  # p should be finished
    ....:
    False
"""
#*****************************************************************************
#       Copyright (C) 2013 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "signals.pxi"
cimport libc.errno


cpdef int get_fileno(f) except -1:
    """
    Return the file descriptor of ``f``.

    INPUT:

    - ``f`` -- an object with a ``.fileno`` method or an integer,
      which is a file descriptor.

    OUTPUT: A C ``long`` representing the file descriptor.

    EXAMPLES::

        sage: from sage.ext.pselect import get_fileno
        sage: get_fileno(open(os.devnull))  # random
        5
        sage: get_fileno(42)
        42
        sage: get_fileno(None)
        Traceback (most recent call last):
        ...
        TypeError: an integer is required
        sage: get_fileno(-1)
        Traceback (most recent call last):
        ...
        ValueError: Invalid file descriptor
        sage: get_fileno(2^30)
        Traceback (most recent call last):
        ...
        ValueError: Invalid file descriptor
    """
    cdef int n
    try:
        n = f.fileno()
    except AttributeError:
        n = f
    if n < 0 or n >= FD_SETSIZE:
        raise ValueError("Invalid file descriptor")
    return n


cdef class PSelecter:
    """
    This class gives an interface to the ``pselect`` system call.

    It can be used in a ``with`` statement to block given signals
    such that they can only occur during the :meth:`pselect()` or
    :meth:`sleep()` calls.

    As an example, we block the ``SIGHUP`` and ``SIGALRM`` signals and
    then raise a ``SIGALRM`` signal. The interrupt will only be seen
    during the :meth:`sleep` call::

        sage: from sage.ext.pselect import PSelecter
        sage: import signal, time
        sage: with PSelecter([signal.SIGHUP, signal.SIGALRM]) as sel:
        ....:     os.kill(os.getpid(), signal.SIGALRM)
        ....:     time.sleep(0.5)  # Simply sleep, no interrupt detected
        ....:     try:
        ....:         _ = sel.sleep(1)  # Interrupt seen here
        ....:     except AlarmInterrupt:
        ....:         print("Interrupt OK")
        ....:
        Interrupt OK

    .. WARNING::

        If ``SIGCHLD`` is blocked inside the ``with`` block, then you
        should not use ``Popen().wait()`` or ``Process().join()``
        because those might block, even if the process has actually
        exited. Use non-blocking alternatives such as ``Popen.poll()``
        or ``multiprocessing.active_children()`` instead.
    """
    cdef sigset_t oldset
    cdef sigset_t blockset

    def __cinit__(self):
        """
        Store old signal mask, needed if this class is used *without*
        a ``with`` statement.

        EXAMPLES::

            sage: from sage.ext.pselect import PSelecter
            sage: PSelecter()
            <sage.ext.pselect.PSelecter ...>
        """
        cdef sigset_t emptyset
        sigemptyset(&emptyset)
        sigprocmask(SIG_BLOCK, &emptyset, &self.oldset)

    def __init__(self, block=[]):
        """
        Store list of signals to block during ``pselect()``.

        EXAMPLES::

            sage: from sage.ext.pselect import PSelecter
            sage: from signal import *
            sage: PSelecter([SIGINT, SIGSEGV])
            <sage.ext.pselect.PSelecter ...>
        """
        sigemptyset(&self.blockset)
        for sig in block:
            sigaddset(&self.blockset, sig)

    def __enter__(self):
        """
        Block signals chosen during :meth:`__init__` in this ``with`` block.

        OUTPUT: ``self``

        TESTS:

        Test nesting, where the inner ``with`` statements should have no
        influence, in particular they should not unblock signals which
        were already blocked upon entering::

            sage: from sage.ext.pselect import PSelecter
            sage: import signal
            sage: with PSelecter([signal.SIGALRM]) as sel:
            ....:     os.kill(os.getpid(), signal.SIGALRM)
            ....:     with PSelecter([signal.SIGFPE]) as sel2:
            ....:         _ = sel2.sleep(0.1)
            ....:     with PSelecter([signal.SIGALRM]) as sel3:
            ....:         _ = sel3.sleep(0.1)
            ....:     try:
            ....:         _ = sel.sleep(0.1)
            ....:     except AlarmInterrupt:
            ....:         print("Interrupt OK")
            ....:
            Interrupt OK
        """
        sigprocmask(SIG_BLOCK, &self.blockset, &self.oldset)
        return self

    def __exit__(self, type, value, traceback):
        """
        Reset signal mask to what it was before :meth:`__enter__`.

        EXAMPLES:

        Install a ``SIGCHLD`` handler::

            sage: import signal
            sage: def child_handler(sig, frame):
            ....:     global got_child
            ....:     got_child = 1
            ....:
            sage: _ = signal.signal(signal.SIGCHLD, child_handler)
            sage: got_child = 0

        Start a process which will cause a ``SIGCHLD`` signal::

            sage: import time
            sage: from multiprocessing import *
            sage: from sage.ext.pselect import PSelecter
            sage: w = PSelecter([signal.SIGCHLD])
            sage: with w:
            ....:     p = Process(target=time.sleep, args=(0.25,))
            ....:     t0 = walltime()
            ....:     p.start()
            ....:

        This ``sleep`` should be interruptible now::

            sage: time.sleep(1)
            sage: t = walltime(t0)
            sage: (0.2 <= t <= 0.9) or t
            True
            sage: got_child
            1
            sage: p.join()
        """
        sigprocmask(SIG_SETMASK, &self.oldset, NULL)

    def pselect(self, rlist=[], wlist=[], xlist=[], timeout=None):
        """
        Wait until one of the given files is ready, or a signal has
        been received, or until ``timeout`` seconds have past.

        INPUT:

        - ``rlist`` -- (default: ``[]``) a list of files to wait for
          reading.

        - ``wlist`` -- (default: ``[]``) a list of files to wait for
          writing.

        - ``xlist`` -- (default: ``[]``) a list of files to wait for
          exceptions.

        - ``timeout`` -- (default: ``None``) a timeout in seconds,
          where ``None`` stands for no timeout.

        OUTPUT: A 4-tuple ``(rready, wready, xready, tmout)`` where the
        first three are lists of file descriptors which are ready,
        that is a subset of ``(rlist, wlist, xlist)``. The fourth is a
        boolean which is ``True`` if and only if the command timed out.
        If ``pselect`` was interrupted by a signal, the output is
        ``([], [], [], False)``.

        .. SEEALSO::

            Use the :meth:`sleep` method instead if you don't care about
            file descriptors.

        EXAMPLES:

        The file ``/dev/null`` should always be available for reading
        and writing::

            sage: from sage.ext.pselect import PSelecter
            sage: f = open(os.devnull, "r+")
            sage: sel = PSelecter()
            sage: sel.pselect(rlist=[f])
            ([<open file '/dev/null', mode 'r+' at ...>], [], [], False)
            sage: sel.pselect(wlist=[f])
            ([], [<open file '/dev/null', mode 'r+' at ...>], [], False)

        A list of various files, all of them should be ready for
        reading. Also create a pipe, which should be ready for
        writing, but not reading (since nothing has been written)::

            sage: f = open(os.devnull, "r")
            sage: g = open(os.path.join(SAGE_LOCAL, 'bin', 'python'), "r")
            sage: (pr, pw) = os.pipe()
            sage: r, w, x, t = PSelecter().pselect([f,g,pr,pw], [pw], [pr,pw])
            sage: len(r), len(w), len(x), t
            (2, 1, 0, False)

        Checking for exceptions on the pipe should simply time out::

            sage: sel.pselect(xlist=[pr,pw], timeout=0.2)
            ([], [], [], True)

        TESTS:

        It is legal (but silly) to list the same file multiple times::

            sage: r, w, x, t = PSelecter().pselect([f,g,f,f,g])
            sage: len(r)
            5

        Invalid input::

            sage: PSelecter().pselect([None])
            Traceback (most recent call last):
            ...
            TypeError: an integer is required

        Open a file and close it, but save the (invalid) file
        descriptor::

            sage: f = open(os.devnull, "r")
            sage: n = f.fileno()
            sage: f.close()
            sage: PSelecter().pselect([n])
            Traceback (most recent call last):
            ...
            IOError: ...
        """
        # Convert given lists to fd_set
        cdef fd_set rfds, wfds, xfds
        FD_ZERO(&rfds)
        FD_ZERO(&wfds)
        FD_ZERO(&xfds)

        cdef int nfds = 0
        cdef int n
        for f in rlist:
            n = get_fileno(f)
            if (n >= nfds): nfds = n + 1
            FD_SET(n, &rfds)
        for f in wlist:
            n = get_fileno(f)
            if (n >= nfds): nfds = n + 1
            FD_SET(n, &wfds)
        for f in xlist:
            n = get_fileno(f)
            if (n >= nfds): nfds = n + 1
            FD_SET(n, &xfds)

        cdef double tm
        cdef timespec tv
        cdef timespec *ptv = NULL
        cdef int ret
        if timeout is not None:
            tm = timeout
            if tm < 0:
                tm = 0
            tv.tv_sec = <long>tm
            tv.tv_nsec = <long>(1e9 * (tm - <double>tv.tv_sec))
            ptv = &tv

        ret = pselect(nfds, &rfds, &wfds, &xfds, ptv, &self.oldset)
        # No file descriptors ready => timeout
        if ret == 0:
            return ([], [], [], True)

        # Error?
        cdef int err
        if ret < 0:
            # Save this value in case Python statements change it.
            err = libc.errno.errno
            if err == libc.errno.EINTR:
                return ([], [], [], False)
            import os
            if err == libc.errno.EBADF:
                raise IOError(err, os.strerror(err))
            else:
                raise OSError(err, os.strerror(err))

        # Figure out which file descriptors to return
        rready = []
        wready = []
        xready = []
        for f in rlist:
            n = get_fileno(f)
            if FD_ISSET(n, &rfds):
                rready.append(f)
        for f in wlist:
            n = get_fileno(f)
            if FD_ISSET(n, &wfds):
                wready.append(f)
        for f in xlist:
            n = get_fileno(f)
            if FD_ISSET(n, &xfds):
                xready.append(f)

        return (rready, wready, xready, False)

    def sleep(self, timeout=None):
        """
        Wait until a signal has been received, or until ``timeout``
        seconds have past.

        This is implemented as a special case of :meth:`pselect` with
        empty lists of file descriptors.

        INPUT:

        - ``timeout`` -- (default: ``None``) a timeout in seconds,
          where ``None`` stands for no timeout.

        OUTPUT: A boolean which is ``True`` if the call timed out,
        False if it was interrupted.

        EXAMPLES:

        A simple wait with timeout::

            sage: from sage.ext.pselect import PSelecter
            sage: sel = PSelecter()
            sage: sel.sleep(timeout=0.1)
            True

        0 or negative time-outs are allowed, ``sleep`` should then
        return immediately::

            sage: sel.sleep(timeout=0)
            True
            sage: sel.sleep(timeout=-123.45)
            True
        """
        return self.pselect(timeout=timeout)[3]
