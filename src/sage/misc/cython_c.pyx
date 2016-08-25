"On-the-fly generation of compiled extensions"

from sage.misc.temporary_file import tmp_filename
from sage.misc.cython import cython_import_all
from sage.repl.user_globals import get_globals


def cython_compile(code,
          verbose=False, compile_message=False,
          make_c_file_nice=False, use_cache=False):
    """
    Given a block of Cython code (as a text string), this function
    compiles it using a C compiler, and includes it into the global
    scope of the module that called this function.

    The following pragmas are available:

    - ``clang`` - may be either c or c++ (or C or C++) indicating
      whether a C or C++ compiler should be used.

    - ``clib`` - additional libraries to be linked in, the space
      separated list is split and passed to distutils.

    - ``cinclude`` - additional directories to search for header
      files. The space separated list is split and passed to
      distutils.

    - ``cfile`` - additional C or C++ files to be compiled. Also,
      ``$SAGE_ROOT`` is expanded, but other environment variables
      are not.

    - ``cargs`` - additional parameters passed to the compiler

    For example::

        #clang C++
        #clib givaro
        #cinclude /usr/local/include/
        #cargs -ggdb
        #cfile foo.c

    AUTHOR: William Stein, 2006-10-31

    .. warn:

        Only use this from Python code, not from extension code, since
        from extension code you would change the global scope (i.e.,
        of the Sage interpreter). And it would be stupid, since you're
        already writing Cython!

        Also, never use this in the standard Sage library.  Any code
        that uses this can only run on a system that has a C compiler
        installed, and we want to avoid making that assumption for
        casual Sage usage.  Also, any code that uses this in the
        library would greatly slow down startup time, since currently
        there is no caching.

    .. todo:

        Need to create a clever caching system so code only gets
        compiled once.
    """
    tmpfile = tmp_filename(ext=".spyx")
    open(tmpfile,'w').write(code)
    cython_import_all(tmpfile, get_globals(),
                      verbose=verbose, compile_message=compile_message,
                      use_cache=use_cache,
                      create_local_c_file=False)
