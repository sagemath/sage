# distutils: libraries = CYGWIN_SQLITE3_LIBS
"""
Miscellaneous operating system functions
"""

from posix.unistd cimport dup, dup2, close
from cpython.exc cimport PyErr_SetFromErrno

import os
import contextlib


def have_program(program, path=None):
    """
    Return ``True`` if a ``program`` executable is found in the path
    given by ``path``.

    INPUT:

    - ``program`` - a string, the name of the program to check.

    - ``path`` - string or None. Paths to search for ``program``,
      separated by ``os.pathsep``. If ``None``, use the :envvar:`PATH`
      environment variable.

    OUTPUT: bool

    EXAMPLES::

        sage: from sage.misc.sage_ostools import have_program
        sage: have_program('ls')
        True
        sage: have_program('there_is_not_a_program_with_this_name')
        False
        sage: from sage.env import SAGE_VENV
        sage: have_program('sage', os.path.join(SAGE_VENV, 'bin'))
        True
        sage: have_program('sage', '/there_is_not_a_path_with_this_name')
        False
        sage: have_program('there_is_not_a_program_with_this_name', os.path.join(SAGE_VENV, 'bin'))
        False
    """
    if path is None:
        path = os.environ.get('PATH', "")
    for p in path.split(os.pathsep):
        try:
            if os.access(os.path.join(p, program), os.X_OK):
                return True
        except OSError:
            pass
    return False


@contextlib.contextmanager
def restore_cwd(chdir=None):
    """
    Context manager that restores the original working directory upon exiting.

    INPUT:

    - ``chdir`` -- optionally change directories to the given directory
      upon entering the context manager

    EXAMPLES::

        sage: import os
        sage: from sage.misc.sage_ostools import restore_cwd
        sage: from sage.misc.misc import SAGE_TMP
        sage: cwd = os.getcwd()
        sage: with restore_cwd(str(SAGE_TMP)):
        ....:     print(os.getcwd() == os.path.realpath(SAGE_TMP))
        True
        sage: cwd == os.getcwd()
        True
    """
    orig_cwd = os.getcwd()
    if chdir is not None:
        os.chdir(chdir)
    try:
        yield
    finally:
        os.chdir(orig_cwd)


cdef file_and_fd(x, int* fd):
    """
    If ``x`` is a file, return ``x`` and set ``*fd`` to its file
    descriptor. If ``x`` is an integer, return ``None`` and set
    ``*fd`` to ``x``. Otherwise, set ``*fd = -1`` and raise a
    ``TypeError``.
    """
    fd[0] = -1
    try:
        m = x.fileno
    except AttributeError:
        try:
            fd[0] = x
        except Exception:  # Some objects raise bizarre exceptions when calling int()
            raise TypeError(f"{x} must be a Python file or an integer")
        return None
    fd[0] = m()
    return x


cdef class redirection:
    r"""
    Context to implement redirection of files, analogous to the
    ``>file`` or ``1>&2`` syntax in POSIX shells.

    Unlike the ``redirect_stdout`` and ``redirect_stderr`` contexts in
    the Python 3.4 standard library, this acts on the OS level, not on
    the Python level. This implies that it only works for true files,
    not duck-type file objects such as ``StringIO``.

    INPUT:

    - ``source`` -- the file to be redirected

    - ``dest`` -- where the source file should be redirected to

    - ``close`` -- (boolean, default: ``True``) whether to close the
      destination file upon exiting the context. This is only supported
      if ``dest`` is a Python file.

    The ``source`` and ``dest`` arguments can be either Python files
    or file descriptors.

    EXAMPLES::

        sage: from sage.misc.sage_ostools import redirection
        sage: fn = tmp_filename()

    ::

        sage: with redirection(sys.stdout, open(fn, 'w')):
        ....:     print("hello world!")
        sage: with open(fn) as f:
        ....:     _ = sys.stdout.write(f.read())
        hello world!

    We can do the same using a file descriptor as source::

        sage: fd = sys.stdout.fileno()
        sage: with redirection(fd, open(fn, 'wb')):
        ....:     _ = os.write(fd, b"hello world!\n")
        sage: with open(fn) as f:
        ....:     _ = sys.stdout.write(f.read())
        hello world!

    The converse also works::

        sage: with open(fn, 'w') as f:
        ....:     _ = f.write("This goes to the file\n")
        ....:     with redirection(f, sys.stdout, close=False):
        ....:         _ = f.write("This goes to stdout\n")
        ....:     _ = f.write("This goes to the file again\n")
        This goes to stdout
        sage: with open(fn) as f:
        ....:     _ = sys.stdout.write(f.read())
        This goes to the file
        This goes to the file again

    The same :class:`redirection` instance can be reused multiple times,
    provided that ``close=False``::

        sage: f = open(fn, 'w+')
        sage: r = redirection(sys.stdout, f, close=False)
        sage: with r:
        ....:     print("Line 1")
        sage: with r:
        ....:     print("Line 2")
        sage: with f:
        ....:     _ = f.seek(0)
        ....:     _ = sys.stdout.write(f.read())
        Line 1
        Line 2

    The redirection also works for subprocesses::

        sage: import subprocess
        sage: with redirection(sys.stdout, open(fn, 'w')):
        ....:     _ = subprocess.call(["echo", "hello world"])
        sage: with open(fn) as f:
        ....:     _ = sys.stdout.write(f.read())
        hello world

    TESTS::

        sage: import io
        sage: redirection(sys.stdout, io.StringIO())
        Traceback (most recent call last):
        ...
        io.UnsupportedOperation: fileno

    The redirection is removed and the destination file is closed even
    in the case of errors::

        sage: f = open(os.devnull, 'w')
        sage: with redirection(sys.stdout, f):
        ....:     raise KeyboardInterrupt
        Traceback (most recent call last):
        ...
        KeyboardInterrupt
        sage: f.closed
        True

    Reusing a :class:`redirection` instance with ``close=True``::

        sage: f = open(fn, 'w+')
        sage: r = redirection(sys.stdout, f, close=True)
        sage: with r:
        ....:     print("Line 1")
        sage: with r:
        ....:     print("Line 2")
        Traceback (most recent call last):
        ...
        ValueError: invalid destination file

    Using ``close=True`` on a non-file::

        sage: redirection(sys.stdout, 0, close=True)
        Traceback (most recent call last):
        ...
        ValueError: close=True requires a Python destination file

    Passing invalid file descriptors::

        sage: with redirection(-123, open(os.devnull, 'w')):
        ....:     pass
        Traceback (most recent call last):
        ...
        OSError: [Errno 9] Bad file descriptor
        sage: with redirection(sys.stdout, -123):
        ....:     pass
        Traceback (most recent call last):
        ...
        OSError: [Errno 9] Bad file descriptor

    Nesting the same :class:`redirection` object is not allowed::

        sage: with redirection(sys.stdout, open(os.devnull, 'w')) as r:
        ....:     with r:
        ....:         pass
        Traceback (most recent call last):
        ...
        ValueError: source already redirected
    """
    cdef readonly source_file, dest_file
    cdef readonly int source_fd, dest_fd, dup_source_fd
    cdef bint close_dest

    def __cinit__(self):
        self.source_fd = -1
        self.dest_fd = -1
        self.dup_source_fd = -1

    def __init__(self, source, dest, close=None):
        self.source_file = file_and_fd(source, &self.source_fd)
        self.dest_file = file_and_fd(dest, &self.dest_fd)

        if self.dest_file is None:
            if close:
                raise ValueError("close=True requires a Python destination file")
        elif close is None:
            close = True
        self.close_dest = close

    cdef int flush(self) except -1:
        for f in (self.source_file, self.dest_file):
            if f is not None:
                f.flush()

    def __enter__(self):
        # Basic sanity checks
        if self.source_fd == -1:
            raise ValueError("invalid source file")
        if self.dest_fd == -1:
            raise ValueError("invalid destination file")
        if self.dup_source_fd != -1:
            raise ValueError("source already redirected")

        self.flush()

        dupsrc = dup(self.source_fd)
        if dupsrc == -1:
            PyErr_SetFromErrno(OSError)
        fd = dup2(self.dest_fd, self.source_fd)
        if fd == -1:
            try:
                PyErr_SetFromErrno(OSError)
            finally:
                close(dupsrc)
        self.dup_source_fd = dupsrc
        return self

    def __exit__(self, *args):
        try:
            self.flush()
        finally:
            fd = dup2(self.dup_source_fd, self.source_fd)
            if fd == -1:
                PyErr_SetFromErrno(OSError)
            close(self.dup_source_fd)
            self.dup_source_fd = -1
            if self.close_dest:
                self.dest_file.close()
                self.dest_fd = -1


IF PY_PLATFORM == 'cygwin':
    from libc.stddef cimport wchar_t

    cdef extern from "Windows.h":
        int SetDllDirectoryW(wchar_t* lpPathName)

    cdef extern from "sqlite3.h":
        int sqlite3_initialize()

    def fix_for_ticket_30157():
        """
        Cygwin-only workaround for an issue caused by the sqlite3 library.

        See :trac:`30157`.

        The issue here is that when the sqlite3 library is first initialized
        it modifies Windows' default DLL search path order, which can possibly
        break the correct search path for subsequent DLL loads.

        This workaround ensures that the sqlite3 library is initialized very
        early on (this does not add any significant overhead) and then
        immediately undoes its deleterious effects.  In particular, calling
        SetDllDirectoryW(NULL) restores the default DLL search path.

        To be clear, there's no reason sqlite3 needs this to function
        correctly; it's just a poorly-considered hack that attempted to work
        around a problem that doesn't affect us.

        Returns 0 if it succeeeds or a non-zero value if not.
        """

        ret = sqlite3_initialize()

        if ret != 0:
            # Library initialization failed for some reason
            return ret

        # SetDllDirectory returns 1 if it succeeds.
        return not SetDllDirectoryW(NULL)
