"""
Cython support functions

AUTHORS:

- William Stein (2006-01-18): initial version
- William Stein (2007-07-28): update from sagex to cython
- Martin Albrecht & William Stein (2011-08): cfile & cargs
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import builtins
import os
import sys
import shutil

from sage.env import (SAGE_LOCAL, cython_aliases,
                      sage_include_directories)
from sage.misc.misc import SPYX_TMP, sage_makedirs
from .temporary_file import tmp_filename
from sage.repl.user_globals import get_globals
from sage.misc.sage_ostools import restore_cwd, redirection
from sage.cpython.string import str_to_bytes
from sage.misc.cachefunc import cached_function


@cached_function
def _standard_libs_libdirs_incdirs_aliases():
    r"""
    Return the list of libraries and library directories.

    EXAMPLES::

        sage: from sage.misc.cython import _standard_libs_libdirs_incdirs_aliases
        sage: _standard_libs_libdirs_incdirs_aliases()
        (['mpfr', 'gmp', 'gmpxx', 'pari', ...],
         [...],
         [...],
         {...})
    """
    aliases = cython_aliases()
    standard_libs = [
        'mpfr', 'gmp', 'gmpxx', 'pari', 'm',
        'ec', 'gsl',
    ] + aliases["CBLAS_LIBRARIES"] + [
        'ntl']
    standard_libdirs = []
    if SAGE_LOCAL:
        standard_libdirs.append(os.path.join(SAGE_LOCAL, "lib"))
    standard_libdirs.extend(aliases["CBLAS_LIBDIR"] + aliases["NTL_LIBDIR"])
    standard_incdirs = sage_include_directories() + aliases["CBLAS_INCDIR"] + aliases["NTL_INCDIR"]
    return standard_libs, standard_libdirs, standard_incdirs, aliases

################################################################
# If the user attaches a .spyx file and changes it, we have
# to reload an .so.
#
# PROBLEM: Python does not allow one to reload an .so extension module.
# Solution, we create a different .so file and load that one,
# overwriting the definitions of everything in the original .so file.
#
# HOW: By using a sequence_number for each .spyx file; we keep
# these sequence numbers in a dict.
#
################################################################


sequence_number = {}


# Only used on Cygwin; see note below about --enable-auto-image-base
CYTHON_IMAGE_BASE = 0x580000000


def cython(filename, verbose=0, compile_message=False,
           use_cache=False, create_local_c_file=False, annotate=True, sage_namespace=True,
           create_local_so_file=False):
    r"""
    Compile a Cython file. This converts a Cython file to a C (or C++ file),
    and then compiles that. The .c file and the .so file are
    created in a temporary directory.

    INPUT:

    - ``filename`` -- the name of the file to be compiled. Should end with
      'pyx'.

    - ``verbose`` (integer, default 0) -- level of verbosity. A negative
      value ensures complete silence.

    - ``compile_message`` (bool, default False) -- if True, print
      ``'Compiling <filename>...'`` to the standard error.

    - ``use_cache`` (bool, default False) -- if True, check the
      temporary build directory to see if there is already a
      corresponding .so file. If so, and if the .so file is newer than the
      Cython file, don't recompile, just reuse the .so file.

    - ``create_local_c_file`` (bool, default False) -- if True, save a
      copy of the ``.c`` or ``.cpp`` file in the current directory.

    - ``annotate`` (bool, default True) -- if True, create an html file which
      annotates the conversion from .pyx to .c. By default this is only created
      in the temporary directory, but if ``create_local_c_file`` is also True,
      then save a copy of the .html file in the current directory.

    - ``sage_namespace`` (bool, default True) -- if True, import
      ``sage.all``.

    - ``create_local_so_file`` (bool, default False) -- if True, save a
      copy of the compiled .so file in the current directory.

    OUTPUT: a tuple ``(name, dir)`` where ``name`` is the name
    of the compiled module and ``dir`` is the directory containing
    the generated files.

    TESTS:

    Before :trac:`12975`, it would have been needed to write ``#clang c++``,
    but upper case ``C++`` has resulted in an error.
    Using pkgconfig to find the libraries, headers and macros. This is a
    work around while waiting for :trac:`22461` which will offer a better
    solution::

        sage: code = [
        ....: "#clang C++",
        ....: "from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular",
        ....: "from sage.libs.singular.polynomial cimport singular_polynomial_pow",
        ....: "def test(MPolynomial_libsingular p):",
        ....: "    singular_polynomial_pow(&p._poly, p._poly, 2, p._parent_ring)"]
        sage: cython(os.linesep.join(code))

    The function ``test`` now manipulates internal C data of polynomials,
    squaring them::

        sage: P.<x,y>=QQ[]
        sage: test(x)
        sage: x
        x^2

    Check that compiling C++ code works::

        sage: cython("# distutils: language = c++\n"+
        ....:        "from libcpp.vector cimport vector\n"
        ....:        "cdef vector[int] * v = new vector[int](4)\n")

    Check that compiling C++ code works when creating a local C file,
    first moving to a tempdir to avoid clutter.  Before :trac:`22113`,
    the create_local_c_file argument was not tested for C++ code::

        sage: d = sage.misc.temporary_file.tmp_dir()
        sage: os.chdir(d)
        sage: with open("test.pyx", 'w') as f:
        ....:     _ = f.write("# distutils: language = c++\n"
        ....:       "from libcpp.vector cimport vector\n"
        ....:       "cdef vector[int] * v = new vector[int](4)\n")
        sage: output = sage.misc.cython.cython("test.pyx", create_local_c_file=True)

    Accessing a ``.pxd`` file from the current directory works::

        sage: d = sage.misc.temporary_file.tmp_dir()
        sage: os.chdir(d)
        sage: with open("helper.pxd", 'w') as f:
        ....:     _ = f.write("cdef inline int the_answer(): return 42")
        sage: cython('''
        ....: from helper cimport the_answer
        ....: print(the_answer())
        ....: ''')
        42

    Warning and error messages generated by Cython are properly
    handled. Warnings are only shown if verbose >= 0::

        sage: code = '''
        ....: def test_unreachable():
        ....:     raise Exception
        ....:     return 42
        ....: '''
        sage: cython(code, verbose=-1)
        sage: cython(code, verbose=0)
        warning: ...:4:4: Unreachable code

        sage: cython("foo = bar\n")
        Traceback (most recent call last):
        ...
        RuntimeError: Error compiling Cython file:
        ------------------------------------------------------------
        ...
        foo = bar
             ^
        ------------------------------------------------------------
        <BLANKLINE>
        ...:1:6: undeclared name not builtin: bar

        sage: cython("cdef extern from 'no_such_header_file': pass")
        Traceback (most recent call last):
        ...
        RuntimeError: ...

    As of :trac:`29139` the default is ``cdivision=True``::

        sage: cython('''
        ....: cdef size_t foo = 3/2
        ....: ''')
    """
    if not filename.endswith('pyx'):
        print("Warning: file (={}) should have extension .pyx".format(filename), file=sys.stderr)

    # base is the name of the .so module that we create. If we are
    # creating a local shared object file, we use a more natural
    # naming convention. If we are not creating a local shared object
    # file, the main constraint is that it is unique and determined by
    # the file that we're running Cython on, so that in some cases we
    # can cache the result (e.g., recompiling the same pyx file during
    # the same session).
    if create_local_so_file:
        base, ext = os.path.splitext(os.path.basename(filename))
    else:
        base = os.path.abspath(filename)
    base = sanitize(base)

    # This is the *temporary* directory where we store the pyx file.
    # This is deleted when Sage exits, which means pyx files must be
    # rebuilt every time Sage is restarted at present.
    target_dir = os.path.join(SPYX_TMP, base)

    # Build directory for Cython/distutils
    build_dir = os.path.join(target_dir, "build")

    if os.path.exists(target_dir):
        # There is already a module here. Maybe we do not have to rebuild?
        # Find the name.
        if use_cache:
            from sage.misc.sageinspect import loadable_module_extension
            prev_so = [F for F in os.listdir(target_dir) if F.endswith(loadable_module_extension())]
            if len(prev_so) > 0:
                prev_so = prev_so[0]     # should have length 1 because of deletes below
                if os.path.getmtime(filename) <= os.path.getmtime('%s/%s' % (target_dir, prev_so)):
                    # We do not have to rebuild.
                    return prev_so[:-len(loadable_module_extension())], target_dir

        # Delete all ordinary files in target_dir
        for F in os.listdir(target_dir):
            G = os.path.join(target_dir, F)
            if os.path.isdir(G):
                continue
            try:
                os.unlink(G)
            except OSError:
                pass
    else:
        sage_makedirs(target_dir)

    if create_local_so_file:
        name = base
    else:
        global sequence_number
        if base not in sequence_number:
            sequence_number[base] = 0
        name = '%s_%s' % (base, sequence_number[base])

        # increment the sequence number so will use a different one next time.
        sequence_number[base] += 1

    if compile_message:
        sys.stderr.write("Compiling {}...\n".format(filename))
        sys.stderr.flush()

    # Copy original file to the target directory.
    pyxfile = os.path.join(target_dir, name + ".pyx")
    shutil.copy(filename, pyxfile)

    # Add current working directory to includes. This is needed because
    # we cythonize from a different directory. See Trac #24764.
    standard_libs, standard_libdirs, standard_includes, aliases = _standard_libs_libdirs_incdirs_aliases()
    includes = [os.getcwd()] + standard_includes

    # Now do the actual build, directly calling Cython and distutils
    from Cython.Build import cythonize
    from Cython.Compiler.Errors import CompileError
    import Cython.Compiler.Options

    try:
        # Import setuptools before importing distutils, so that setuptools
        # can replace distutils by its own vendored copy.
        import setuptools
        from setuptools.dist import Distribution
        from setuptools.extension import Extension
    except ImportError:
        # Fall back to distutils (stdlib); note that it is deprecated
        # in Python 3.10, 3.11; https://www.python.org/dev/peps/pep-0632/
        from distutils.dist import Distribution
        from distutils.core import Extension

    from distutils.log import set_verbosity
    set_verbosity(verbose)

    Cython.Compiler.Options.annotate = annotate
    Cython.Compiler.Options.embed_pos_in_docstring = True
    Cython.Compiler.Options.pre_import = "sage.all" if sage_namespace else None

    extra_compile_args = ['-w']  # no warnings
    extra_link_args = []
    if sys.platform == 'cygwin':
        # Link using --enable-auto-image-base, reducing the likelihood of the
        # new DLL colliding with existing DLLs in memory.
        # Note: Cygwin locates --enable-auto-image-base DLLs in the range
        # 0x4:00000000 up to 0x6:00000000; this is documented in heap.cc in the
        # Cygwin sources, while 0x6:00000000 and up is reserved for the Cygwin
        # heap.
        # In practice, Sage has enough DLLs that when rebasing everything we
        # use up through, approximately, 0x4:80000000 though there is nothing
        # precise here.  When running 'rebase' it just start from 0x2:00000000
        # (below that is reserved for Cygwin's DLL and some internal
        # structures).
        # Therefore, to minimize the likelihood of collision with one of Sage's
        # standard DLLs, while giving ~2GB (should be more than enough) for
        # Sage to grow, we base these DLLs from 0x5:8000000, leaving again ~2GB
        # for temp DLLs which in normal use should be more than enough.
        # See https://trac.sagemath.org/ticket/28258
        # It should be noted, this is not a bulletproof solution; there is
        # still a potential for DLL overlaps with this.  But this reduces the
        # probability thereof, especially in normal practice.
        dll_filename = os.path.splitext(pyxfile)[0] + '.dll'
        image_base = _compute_dll_image_base(dll_filename)
        extra_link_args.extend(['-Wl,--disable-auto-image-base',
                                '-Wl,--image-base=0x{:x}'.format(image_base)])

    ext = Extension(name,
                    sources=[pyxfile],
                    extra_compile_args=extra_compile_args,
                    extra_link_args=extra_link_args,
                    libraries=standard_libs,
                    library_dirs=standard_libdirs)

    directives = dict(language_level=sys.version_info[0], cdivision=True)

    try:
        # Change directories to target_dir so that Cython produces the correct
        # relative path; https://trac.sagemath.org/ticket/24097
        with restore_cwd(target_dir):
            try:
                ext, = cythonize([ext],
                                 aliases=aliases,
                                 include_path=includes,
                                 compiler_directives=directives,
                                 quiet=(verbose <= 0),
                                 errors_to_stderr=False,
                                 use_listing_file=True)
            finally:
                # Read the "listing file" which is the file containing
                # warning and error messages generated by Cython.
                try:
                    with open(name + ".lis") as f:
                        cython_messages = f.read()
                except IOError:
                    cython_messages = "Error compiling Cython file"
    except CompileError:
        raise RuntimeError(cython_messages.strip())

    if verbose >= 0:
        sys.stderr.write(cython_messages)
        sys.stderr.flush()

    if create_local_c_file:
        shutil.copy(os.path.join(target_dir, ext.sources[0]),
                    os.curdir)
        if annotate:
            shutil.copy(os.path.join(target_dir, name + ".html"),
                        os.curdir)

    # This emulates running "setup.py build" with the correct options
    dist = Distribution()
    dist.ext_modules = [ext]
    dist.include_dirs = includes
    buildcmd = dist.get_command_obj("build")
    buildcmd.build_base = build_dir
    buildcmd.build_lib = target_dir

    try:
        # Capture errors from distutils and its child processes
        with open(os.path.join(target_dir, name + ".err"), 'w+') as errfile:
            try:
                # Redirect stderr to errfile.  We use the file descriptor
                # number "2" instead of "sys.stderr" because we really
                # want to redirect the messages from GCC. These are sent
                # to the actual stderr, regardless of what sys.stderr is.
                sys.stderr.flush()
                with redirection(2, errfile, close=False):
                    dist.run_command("build")
            finally:
                errfile.seek(0)
                distutils_messages = errfile.read()
    except Exception as msg:
        msg = str(msg) + "\n" + distutils_messages
        raise RuntimeError(msg.strip())

    if verbose >= 0:
        sys.stderr.write(distutils_messages)
        sys.stderr.flush()

    if create_local_so_file:
        # Copy module to current directory
        from sage.misc.sageinspect import loadable_module_extension
        shutil.copy(os.path.join(target_dir, name + loadable_module_extension()),
                    os.curdir)

    return name, target_dir


################################################################
# COMPILE
################################################################
def cython_lambda(vars, expr, verbose=0, **kwds):
    """
    Create a compiled function which evaluates ``expr`` assuming machine values
    for ``vars``.

    INPUT:

    - ``vars`` - list of pairs (variable name, c-data type), where the variable
      names and data types are strings, OR a string such as ``'double x, int y,
      int z'``

    - ``expr`` - an expression involving the vars and constants; you can access
      objects defined in the current module scope ``globals()`` using
      ``sage.object_name``.

    .. warning::

        Accessing ``globals()`` doesn't actually work, see :trac:`12446`.

    EXAMPLES:

    We create a Lambda function in pure Python (using the r to make sure the 3.2
    is viewed as a Python float)::

        sage: f = lambda x,y: x*x + y*y + x + y + 17r*x + 3.2r

    We make the same Lambda function, but in a compiled form. ::

        sage: g = cython_lambda('double x, double y', 'x*x + y*y + x + y + 17*x + 3.2')
        sage: g(2,3)
        55.2
        sage: g(0,0)
        3.2

    In order to access Sage globals, prefix them with ``sage.``::

        sage: f = cython_lambda('double x', 'sage.sin(x) + sage.a')
        sage: f(0)
        Traceback (most recent call last):
        ...
        NameError: global 'a' is not defined
        sage: a = 25
        sage: f(10)
        24.45597888911063
        sage: a = 50
        sage: f(10)
        49.45597888911063
    """
    if isinstance(vars, str):
        v = vars
    else:
        v = ', '.join(['%s %s' % (typ, var) for typ, var in vars])

    s = """
cdef class _s:
    cdef globals

    def __init__(self):
        from sage.repl.user_globals import get_globals
        self.globals = get_globals()

    def __getattr__(self, name):
        try:
            return self.globals[name]
        except KeyError:
            raise NameError("global {!r} is not defined".format(name))

sage = _s()

def f(%s):
    return %s
    """ % (v, expr)
    if verbose > 0:
        print(s)
    tmpfile = tmp_filename(ext=".pyx")
    with open(tmpfile, 'w') as f:
        f.write(s)

    d = {}
    cython_import_all(tmpfile, d, verbose=verbose, **kwds)
    return d['f']


################################################################
# IMPORT
################################################################
def cython_import(filename, **kwds):
    """
    Compile a file containing Cython code, then import and return the
    module.  Raises an ``ImportError`` if anything goes wrong.

    INPUT:

    - ``filename`` - a string; name of a file that contains Cython
      code

    See the function :func:`sage.misc.cython.cython` for documentation
    for the other inputs.

    OUTPUT:

    - the module that contains the compiled Cython code.
    """
    name, build_dir = cython(filename, **kwds)

    oldpath = sys.path
    try:
        sys.path.append(build_dir)
        return builtins.__import__(name)
    finally:
        sys.path = oldpath


def cython_import_all(filename, globals, **kwds):
    """
    Imports all non-private (i.e., not beginning with an underscore)
    attributes of the specified Cython module into the given context.
    This is similar to::

        from module import *

    Raises an ``ImportError`` exception if anything goes wrong.

    INPUT:

    - ``filename`` - a string; name of a file that contains Cython
      code
    """
    m = cython_import(filename, **kwds)
    for k, x in m.__dict__.items():
        if k[0] != '_':
            globals[k] = x


def sanitize(f):
    """
    Given a filename ``f``, replace it by a filename that is a valid Python
    module name.

    This means that the characters are all alphanumeric or ``_``'s and doesn't
    begin with a numeral.

    EXAMPLES::

        sage: from sage.misc.cython import sanitize
        sage: sanitize('abc')
        'abc'
        sage: sanitize('abc/def')
        'abc_def'
        sage: sanitize('123/def-hij/file.py')
        '_123_def_hij_file_py'
    """
    s = ''
    if f[0].isdigit():
        s += '_'
    for a in f:
        if a.isalnum():
            s += a
        else:
            s += '_'
    return s


def compile_and_load(code, **kwds):
    r"""
    INPUT:

    - ``code`` -- string containing code that could be in a .pyx file
      that is attached or put in a %cython block in the notebook.

    OUTPUT: a module, which results from compiling the given code and
    importing it

    EXAMPLES::

        sage: from sage.misc.cython import compile_and_load
        sage: module = compile_and_load("def f(int n):\n    return n*n")
        sage: module.f(10)
        100

    TESTS::

        sage: code = '''
        ....: from sage.rings.rational cimport Rational
        ....: from sage.rings.polynomial.polynomial_rational_flint cimport Polynomial_rational_flint
        ....: from sage.libs.flint.fmpq_poly cimport fmpq_poly_length, fmpq_poly_get_coeff_mpq, fmpq_poly_set_coeff_mpq
        ....:
        ....: def evaluate_at_power_of_gen(Polynomial_rational_flint f, unsigned long n):
        ....:     assert n >= 1
        ....:     cdef Polynomial_rational_flint res = f._new()
        ....:     cdef unsigned long k
        ....:     cdef Rational z = Rational(0)
        ....:     for k in range(fmpq_poly_length(f.__poly)):
        ....:         fmpq_poly_get_coeff_mpq(z.value, f.__poly, k)
        ....:         fmpq_poly_set_coeff_mpq(res.__poly, n*k, z.value)
        ....:     return res
        ....: '''
        sage: module = compile_and_load(code)  # long time
        sage: R.<x> = QQ[]
        sage: module.evaluate_at_power_of_gen(x^3 + x - 7, 5)  # long time
        x^15 + x^5 - 7
    """
    tmpfile = tmp_filename(ext=".pyx")
    with open(tmpfile, 'w') as f:
        f.write(code)
    return cython_import(tmpfile, **kwds)


def cython_compile(code, **kwds):
    """
    Given a block of Cython code (as a text string), this function
    compiles it using a C compiler, and includes it into the global
    namespace.

    AUTHOR: William Stein, 2006-10-31

    .. WARNING::

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

    .. TODO::

        Need to create a clever caching system so code only gets
        compiled once.
    """
    tmpfile = tmp_filename(ext=".pyx")
    with open(tmpfile, 'w') as f:
        f.write(code)
    return cython_import_all(tmpfile, get_globals(), **kwds)


# THe following utility functions are used on Cygwin only to work around a
# shortcoming in ld/binutils; see https://trac.sagemath.org/ticket/28258
def _strhash(s):
    """
    Implementation of the strhash function from binutils

    See
    https://sourceware.org/git/gitweb.cgi?p=binutils-gdb.git;a=blob;f=binutils/dllwrap.c;hb=d0e6d77b3fb64561ca66535f8ed6ca523eac923e#l449

    This is used to re-implement support for --enable-auto-image-base with a
    custom base address, which is currently broken in ld/binutils for PE32+
    binaries.

    TESTS::

        sage: from sage.misc.cython import _strhash
        sage: hex(_strhash('test.pyx'))
        '0x3d99a20'
    """

    s = str_to_bytes(s)
    h = 0
    l = len(s)

    for c in s:
        h += c + (c << 17)
        h ^= h >> 2

    h += l + (l << 17)
    h ^= h >> 2
    return h


def _compute_dll_image_base(filename, auto_image_base=None):
    """
    Compute a DLL image base address from a hash of its filename, using the
    same implementation as --enable-auto-image-base from binutils, but with
    support for a custom minimum base address.

    See
    https://sourceware.org/git/gitweb.cgi?p=binutils-gdb.git;a=blob;f=ld/emultempl/pe.em;h=c18cb266beb23d14108044571382e0c7f4dbedb3;hb=d0e6d77b3fb64561ca66535f8ed6ca523eac923e#l919
    which implements this for PE32 images (for PE32+ images the support for
    custom base addresses is broken).

    See also :trac:`28258`.

    TESTS::

        sage: from sage.misc.cython import _compute_dll_image_base
        sage: hex(_compute_dll_image_base('test.pyx'))
        '0x59a200000'

    Test that given a random filename it is always below 0x6:00000000; in
    principle this could fail randomly, but if it never fails then no problem.
    We do not test for collisions since there is always a random chance of
    collisions with two arbitrary filenames::

        sage: from sage.misc.temporary_file import tmp_filename
        sage: image_base = _compute_dll_image_base(tmp_filename(ext='.pyx'))
        sage: hex(image_base)  # random
        '0x59a200000'
        sage: image_base <= (0x600000000 - 0x10000)
        True
    """

    if auto_image_base is None:
        # The default, which can in principle be modified at runtime if needed
        auto_image_base = CYTHON_IMAGE_BASE

    # See https://sourceware.org/git/gitweb.cgi?p=binutils-gdb.git;a=blob;f=ld/emultempl/pep.em;h=dca36cc341aabfa2ec9db139b8ec94690165201a;hb=d0e6d77b3fb64561ca66535f8ed6ca523eac923e#l113
    # where move_default_addr_high=1 on Cygwin targets.
    # This is actually derived from taking the upper limit 0x6:00000000 minus
    # the lower limit 0x4:00000000 minus a gap of at least 0x10000.  I don't
    # know why but this is supposed to work: https://narkive.com/HCyWCS6h.23
    # Currently in binutils auto_image_base of 0x4:00000000 is hard-coded, so
    # this comes out to 0x1:ffff0000, but here we must compute it.
    NT_DLL_AUTO_IMAGE_MASK = 0x600000000 - auto_image_base - 0x10000

    h = _strhash(filename)
    return auto_image_base + ((h << 16) & NT_DLL_AUTO_IMAGE_MASK)
