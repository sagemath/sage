"""
Miscellaneous operating system functions
"""

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
        sage: have_program('sage', path=SAGE_ROOT)
        True
        sage: have_program('ls', path=SAGE_ROOT)
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

    EXAMPLES:

        sage: import os
        sage: from sage.misc.sage_ostools import restore_cwd
        sage: from sage.misc.misc import SAGE_TMP
        sage: cwd = os.getcwd()
        sage: with restore_cwd(str(SAGE_TMP)):
        ....:     print(os.getcwd() == SAGE_TMP)
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
