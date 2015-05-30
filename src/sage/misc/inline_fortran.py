"""
Fortran compiler
"""
import os, imp, shutil

from sage.misc.temporary_file import tmp_dir


class InlineFortran:
    def __init__(self, globals=None):
        if globals is None:
            self.globs = {}
        else:
            self.globs = globals  # Deprecated
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
        if len(x.splitlines()) == 1 and os.path.exists(x):
            from sage.misc.superseded import deprecation
            deprecation(2891, "Calling fortran() with a filename is deprecated, use fortran(open(f).read) instead")
            filename = x
            x = open(x).read()
            if filename.lower().endswith('.f90'):
                x = '!f90\n' + x

        if globals is None:
            globals = self.globs

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
            extra_args = '--quiet --f77exec=sage-inline-fortran --f90exec=sage-inline-fortran %s %s >"%s" 2>&1'%(
                s_lib_path, s_lib, log)

            f2py.compile(x, name, extra_args = extra_args, source_fn=fortran_file)
            log_string = open(log).read()

            # Note that f2py() doesn't raise an exception if it fails.
            # In that case, the import below will fail.
            try:
                file, pathname, description = imp.find_module(name, [mytmpdir])
            except ImportError:
                raise RuntimeError("failed to compile Fortran code:\n" + log_string)
            try:
                m = imp.load_module(name, file, pathname, description)
            finally:
                file.close()

            if self.verbose:
                print(log_string)
        finally:
            os.chdir(old_cwd)
            try:
                shutil.rmtree(mytmpdir)
            except OSError:
                # This can fail for example over NFS
                pass

        for k, x in m.__dict__.iteritems():
            if k[0] != '_':
                globals[k] = x

    def add_library(self,s):
       self.libraries.append(s)

    def add_library_path(self,s):
       self.library_paths.append(s)

# An instance
fortran = InlineFortran()
