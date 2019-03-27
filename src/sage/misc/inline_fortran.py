"""
Fortran compiler
"""
from __future__ import absolute_import

import os
import shutil
import sys

import six

from six import iteritems

from sage.misc.temporary_file import tmp_dir


def _import_module_from_path(name, path=None):
    """
    Import the module named ``name`` by searching the given path entries (or
    `sys.path` by default).

    Returns a fully executed module object without inserting that module into
    `sys.modules`.

    EXAMPLES:

        sage: from sage.misc.inline_fortran import _import_module_from_path
        sage: modname = '___test__import_module_from_path'
        sage: tmpdir = tmp_dir()
        sage: filename = os.path.join(tmpdir, modname + '.py')
        sage: with open(filename, 'w') as fobj:
        ....:     _ = fobj.write('foo = "bar"')
        sage: mod = _import_module_from_path(modname)
        Traceback (most recent call last):
        ...
        ImportError: No module named ___test__import_module_from_path
        sage: mod = _import_module_from_path('DoEsNoTeXiSt', path=[tmpdir])
        Traceback (most recent call last):
        ...
        ImportError: No module named DoEsNoTeXiSt
        sage: mod = _import_module_from_path(modname, path=[tmpdir])
        sage: mod.foo
        'bar'
    """

    if path is None:
        path = sys.path

    return _import_module_from_path_impl(name, path)


if six.PY2:
    import imp

    def _import_module_from_path_impl(name, path):
        """Implement ``_import_module_from_path for Python 2."""

        # Note: Raises an ImportError if not found
        fileobj, pathname, description = imp.find_module(name, path)
        try:
            # Executes the module in fileobj using the appropriate loader and
            # returns the module
            return imp.load_module(name, fileobj, pathname, description)
        finally:
            fileobj.close()
else:
    import importlib

    def _import_module_from_path_impl(name, path):
        """Implement ``_import_module_from_path for Python 3.4+."""

        # This is remarkably tricky to do right, considering that the new
        # importlib is supposed to make direct interaction with the import
        # system easier.  I blame the ModuleSpec stuff...
        finder = importlib.machinery.PathFinder()
        spec = finder.find_spec(name, path=path)
        if spec is None:
            raise ImportError('No module named {}'.format(name))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod


class InlineFortran:
    def __init__(self, globals=None):
        # globals=None means: use user globals from REPL
        self.globs = globals
        self.library_paths=[]
        self.libraries=[]
        self.verbose = False

    def __repr__(self):
        return "Interface to Fortran compiler"

    def __call__(self, *args, **kwds):
        return self.eval(*args, **kwds)

    def eval(self, x, globals=None, locals=None):
        """
        Compile fortran code ``x`` and adds the functions in it to
        ``globals``.

        INPUT:

        - ``x`` -- Fortran code

        - ``globals`` -- a dict to which to add the functions from the
          fortran module

        - ``locals`` -- ignored

        EXAMPLES::

            sage: code = '''
            ....: C FILE: FIB1.F
            ....:       SUBROUTINE FIB(A,N)
            ....: C
            ....: C     CALCULATE FIRST N FIBONACCI NUMBERS
            ....: C
            ....:       INTEGER N
            ....:       REAL*8 A(N)
            ....:       DO I=1,N
            ....:          IF (I.EQ.1) THEN
            ....:             A(I) = 0.0D0
            ....:          ELSEIF (I.EQ.2) THEN
            ....:             A(I) = 1.0D0
            ....:          ELSE
            ....:             A(I) = A(I-1) + A(I-2)
            ....:          ENDIF
            ....:       ENDDO
            ....:       END
            ....: C END FILE FIB1.F
            ....: '''
            sage: fortran(code, globals())
            sage: import numpy
            sage: a = numpy.array(range(10), dtype=float)
            sage: fib(a, 10)
            sage: a
            array([  0.,   1.,   1.,   2.,   3.,   5.,   8.,  13.,  21.,  34.])

        TESTS::

            sage: os.chdir(SAGE_ROOT)
            sage: fortran.eval("SYNTAX ERROR !@#$")
            Traceback (most recent call last):
            ...
            RuntimeError: failed to compile Fortran code:...
            sage: os.getcwd() == SAGE_ROOT
            True
        """
        if globals is None:
            globals = self.globs
            if globals is None:
                from sage.repl.user_globals import get_globals
                globals = get_globals()

        from numpy import f2py

        # Create everything in a temporary directory
        mytmpdir = tmp_dir()

        try:
            old_cwd = os.getcwd()
            os.chdir(mytmpdir)

            name = "fortran_module"  # Python module name
            # if the first line has !f90 as a comment, gfortran will
            # treat it as Fortran 90 code
            if x.startswith('!f90'):
                fortran_file = name + '.f90'
            else:
                fortran_file = name + '.f'

            s_lib_path = ""
            s_lib = ""
            for s in self.library_paths:
                s_lib_path = s_lib_path + "-L%s "

            for s in self.libraries:
                s_lib = s_lib + "-l%s "%s

            log = name + ".log"
            extra_args = ('--quiet --f77exec=sage-inline-fortran '
                          '--f90exec=sage-inline-fortran {lib_path} {lib} '
                          '> {log} 2>&1'.format(lib_path=s_lib_path,
                                                lib=s_lib, log=log))

            f2py.compile(x, name, extra_args=extra_args,
                         source_fn=fortran_file)

            with open(log) as fobj:
                log_string = fobj.read()

            # Note that f2py() doesn't raise an exception if it fails.
            # In that case, the import below will fail.
            try:
                mod = _import_module_from_path(name, [mytmpdir])
            except ImportError:
                raise RuntimeError("failed to compile Fortran code:\n" +
                                   log_string)

            if self.verbose:
                print(log_string)
        finally:
            os.chdir(old_cwd)

            if sys.platform != 'cygwin':
                # Do not delete temporary DLLs on Cygwin; this will cause
                # future forks of this process to fail.  Instead temporary DLLs
                # will be cleaned up upon process exit
                try:
                    shutil.rmtree(mytmpdir)
                except OSError:
                    # This can fail for example over NFS
                    pass

        for k, x in iteritems(mod.__dict__):
            if k[0] != '_':
                globals[k] = x

    def add_library(self,s):
       self.libraries.append(s)

    def add_library_path(self,s):
       self.library_paths.append(s)

# An instance
fortran = InlineFortran()
