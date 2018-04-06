"""
Cython support functions

AUTHORS:

- William Stein (2006-01-18): initial version
- William Stein (2007-07-28): update from sagex to cython
- Martin Albrecht & William Stein (2011-08): cfile & cargs
"""
#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function, absolute_import
from six.moves import builtins
from six import iteritems

import os
import sys
import re
import shutil
import pkgconfig

from textwrap import dedent

from sage.env import (SAGE_LOCAL, SAGE_SRC, cython_aliases,
        sage_include_directories)
from sage.misc.misc import SPYX_TMP, sage_makedirs
from .temporary_file import tmp_filename
from sage.misc.superseded import deprecated_function_alias
from sage.repl.user_globals import get_globals


# CBLAS can be one of multiple implementations
cblas_pc = pkgconfig.parse('cblas')
cblas_libs = list(cblas_pc['libraries'])
cblas_library_dirs = list(cblas_pc['library_dirs'])
cblas_include_dirs = list(cblas_pc['include_dirs'])

standard_libs = [
    'mpfr', 'gmp', 'gmpxx', 'stdc++', 'pari', 'm',
    'ec', 'gsl',
] + cblas_libs + [
    'ntl']


# Functions which used to be automatically declared.
# We list these here in order to give useful warnings.
old_pxi_names = {
    "cysignals.signals": [
        "sig_on", "sig_str", "sig_check", "sig_off",
        "sig_retry", "sig_error", "sig_block", "sig_unblock",
        "sig_on_no_except", "sig_str_no_except", "sig_check_no_except",
        "cython_check_exception",
    ],
    "sage.ext.stdsage": [
        "PY_NEW", "HAS_DICTIONARY",
    ],
    "cysignals.memory": [
        "sig_malloc", "sig_realloc", "sig_calloc", "sig_free",
        "check_allocarray", "check_reallocarray",
        "check_malloc", "check_realloc", "check_calloc",
    ],
    "libc.string": [
        "strlen", "strcpy", "memset", "memcpy", "memcmp",
    ],
    "libc.math": [
        "sqrt", "frexp", "ldexp",
    ],
    "libc.stdio": [
        "stdin", "stdout", "stderr",
        "FOPEN_MAX", "FILENAME_MAX",
        "fopen", "freopen", "fdopen", "fclose",
        "remove", "rename", "tmpfile",
        "setvbuf", "BUFSIZ", "setbuf",
        "fread", "fwrite", "fflush",
        "EOF", "clearerr", "feof", "ferror",
        "SEEK_SET", "SEEK_CUR", "SEEK_END",
        "fseek", "rewind", "ftell", "fgetpos", "fsetpos",
        "scanf", "sscanf", "fscanf",
        "printf", "sprintf", "snprintf", "fprintf",
        "perror", "gets", "fgets", "getchar", "fgetc", "getc", "ungetc",
        "puts", "fputs", "putchar", "fputc", "putc", "getline",
    ]
    }


def _parse_keywords(kwd, s):
    r"""
    Given a keyword ``kwd`` and a string ``s``, return a list of all arguments
    on the same line as that keyword in ``s``, as well as a new copy of ``s``
    in which each occurrence of ``kwd`` is in a comment. If a comment already
    occurs on the line containing ``kwd``, no words after the ``#`` are added
    to the list.

    EXAMPLES::

        sage: from sage.misc.cython import _parse_keywords
        sage: _parse_keywords('clib', " clib foo bar baz\n #cinclude bar\n")
        doctest:...: DeprecationWarning: the Sage-specific Cython pragma '#clib' is deprecated;
        use '# distutils: libraries' instead
        See http://trac.sagemath.org/24105 for details.
        (['foo', 'bar', 'baz'], ' #clib foo bar baz\n #cinclude bar\n')
        sage: _parse_keywords('clib', "# qux clib foo bar baz\n #cinclude bar\n")
        (['foo', 'bar', 'baz'], '# qux clib foo bar baz\n #cinclude bar\n')
        sage: _parse_keywords('clib', "# clib foo bar # baz\n #cinclude bar\n")
        (['foo', 'bar'], '# clib foo bar # baz\n #cinclude bar\n')

    TESTS::

        sage: from sage.misc.cython import parse_keywords
        sage: parse_keywords('kwd', "#kwd foo")
        doctest:...: DeprecationWarning: parse_keywords is deprecated. Please use sage.misc.cython._parse_keywords instead.
        See http://trac.sagemath.org/24105 for details.
        (['foo'], '#kwd foo')
    """
    j = 0
    v = []
    while True:
        # see if kwd occurs
        i = s[j:].find(kwd)
        if i == -1: break
        j = i + j

        # add a hash, if necessary
        last_hash = s[:j].rfind('#')
        last_newline = s[:j].rfind('\n')
        if last_hash > last_newline:
            j += len(kwd)
        else:
            s = s[:j] + '#' + s[j:]
            j += len(kwd) + 1

        # find all other words on this line
        k = s[j:].find('\n')
        if k == -1:
            k = len(s)

        # add them to our list, until we find a comment
        for X in s[j:j+k].split():
            if X[0] == '#':   # skip rest of line
                break
            v.append(X)

    if v:
        replacements = {
            "clang": "# distutils: language",
            "clib": "# distutils: libraries",
            "cfile": "# distutils: sources",
            "cinclude": "# distutils: include_dirs",
            "cargs": "# distutils: extra_compile_args",
        }
        if kwd in replacements:
            from sage.misc.superseded import deprecation
            deprecation(24105, "the Sage-specific Cython pragma {!r} is deprecated;\n"
                "use {!r} instead".format("#" + kwd, replacements[kwd]))
    return v, s


def _environ_parse(s):
    r"""
    Given a string s, find each substring of the form ``'\$ABC'``. If the
    environment variable :envvar:`$ABC` is set, replace ``'\$ABC'`` with its
    value and move on to the next such substring. If it is not set, stop
    parsing there.

    EXAMPLES::

        sage: from sage.misc.cython import _environ_parse
        sage: _environ_parse('$SAGE_LOCAL') == SAGE_LOCAL
        True
        sage: _environ_parse('$THIS_IS_NOT_DEFINED_ANYWHERE')
        '$THIS_IS_NOT_DEFINED_ANYWHERE'
        sage: os.environ['DEFINE_THIS'] = 'hello'
        sage: _environ_parse('$DEFINE_THIS/$THIS_IS_NOT_DEFINED_ANYWHERE/$DEFINE_THIS')
        'hello/$THIS_IS_NOT_DEFINED_ANYWHERE/$DEFINE_THIS'

    TESTS::

        sage: from sage.misc.cython import environ_parse
        sage: environ_parse('$SAGE_LOCAL') == SAGE_LOCAL
        doctest:...: DeprecationWarning: environ_parse is deprecated. Please use sage.misc.cython._environ_parse instead.
        See http://trac.sagemath.org/24105 for details.
        True
    """
    i = s.find('$')
    if i == -1:
        return s
    j = s[i:].find('/')
    if j == -1:
        j = len(s)
    else:
        j = i + j
    name = s[i+1:j]
    if name in os.environ:
        s = s[:i] + os.environ[name] + s[j:]
    else:
        return s
    return _environ_parse(s)


def _pyx_preparse(s):
    r"""
    Preparse a pyx file:

    * parse ``clang`` pragma (c or c++)
    * parse ``clib`` pragma (additional libraries to link in)
    * parse ``cinclude`` (additional include directories)
    * parse ``cfile`` (additional files to be included)
    * parse ``cargs`` (additional parameters passed to the compiler)

    The pragmas:

    - ``clang`` - may be either ``'c'`` or ``'c++'`` indicating whether a C or
      C++ compiler should be used

    - ``clib`` - additional libraries to be linked in, the space separated list
      is split and passed to distutils.

    - ``cinclude`` - additional directories to search for header files. The
      space separated list is split and passed to distutils.

    - ``cfile`` - additional C or C++ files to be compiled. Also,
      :envvar:`$SAGE_SRC` and :envvar:`$SAGE_LOCAL` are expanded, but other
      environment variables are not.

    - ``cargs`` - additional parameters passed to the compiler

    OUTPUT: preamble, libs, includes, language, files, args

    EXAMPLES::

        sage: from sage.misc.cython import _pyx_preparse
        sage: _pyx_preparse("")
        ('',
        ['mpfr',
        'gmp',
        'gmpxx',
        'stdc++',
        'pari',
        'm',
        'ec',
        'gsl',
        ...,
        'ntl'],
        ['.../include',
        '.../include/python...',
        '.../python.../numpy/core/include',
        '...',
        '.../sage/ext',
        '.../cysignals'],
        'c',
        [], ['-w', '-O2'],...)
        sage: s, libs, inc, lang, f, args, libdirs = _pyx_preparse("# clang c++\n #clib foo\n # cinclude bar\n")
        doctest:...: DeprecationWarning: the Sage-specific Cython pragma '#clang' is deprecated;
        use '# distutils: language' instead
        See http://trac.sagemath.org/24105 for details.
        doctest:...: DeprecationWarning: the Sage-specific Cython pragma '#clib' is deprecated;
        use '# distutils: libraries' instead
        See http://trac.sagemath.org/24105 for details.
        doctest:...: DeprecationWarning: the Sage-specific Cython pragma '#cinclude' is deprecated;
        use '# distutils: include_dirs' instead
        See http://trac.sagemath.org/24105 for details.
        sage: lang
        'c++'

        sage: libs
        ['foo', 'mpfr',
        'gmp', 'gmpxx',
        'stdc++',
        'pari',
        'm',
        'ec',
        'gsl',
        ...,
        'ntl']
        sage: libs[1:] == sage.misc.cython.standard_libs
        True

        sage: inc
        ['bar',
        '.../include',
        '.../include/python...',
        '.../python.../numpy/core/include',
        '...',
        '.../sage/ext',
        '.../cysignals']

        sage: s, libs, inc, lang, f, args, libdirs = _pyx_preparse("# cargs -O3 -ggdb\n")
        doctest:...: DeprecationWarning: the Sage-specific Cython pragma '#cargs' is deprecated;
        use '# distutils: extra_compile_args' instead
        See http://trac.sagemath.org/24105 for details.
        sage: args
        ['-w', '-O2', '-O3', '-ggdb']

    TESTS::

        sage: module = sage.misc.cython.import_test("trac11680")  # long time (7s on sage.math, 2012)
        sage: R.<x> = QQ[]
        sage: module.evaluate_at_power_of_gen(x^3 + x - 7, 5)  # long time
        x^15 + x^5 - 7

    ::

        sage: from sage.misc.cython import pyx_preparse
        sage: _ = pyx_preparse("")
        doctest:...: DeprecationWarning: pyx_preparse is deprecated. Please use sage.misc.cython._pyx_preparse instead.
        See http://trac.sagemath.org/24105 for details.
    """
    lang, s = _parse_keywords('clang', s)
    if lang:
        lang = lang[0].lower() # this allows both C++ and c++
    else:
        lang = "c"

    v, s = _parse_keywords('clib', s)
    libs = v + standard_libs

    additional_source_files, s = _parse_keywords('cfile', s)

    v, s = _parse_keywords('cinclude', s)
    inc = [_environ_parse(x.replace('"','').replace("'","")) for x in v] + sage_include_directories()
    args, s = _parse_keywords('cargs', s)
    args = ['-w','-O2'] + args
    libdirs = cblas_library_dirs

    # Add cysignals directory to includes
    for path in sys.path:
        cysignals_path = os.path.join(path, "cysignals")
        if os.path.isdir(cysignals_path):
            inc.append(cysignals_path)

    return s, libs, inc, lang, additional_source_files, args, libdirs


parse_keywords = deprecated_function_alias(24105, _parse_keywords)
environ_parse = deprecated_function_alias(24105, _environ_parse)
pyx_preparse = deprecated_function_alias(24105, _pyx_preparse)


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

    - ``verbose`` (integer, default 0) -- level of verbosity.

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

        sage: import sage.misc.cython
        sage: d = sage.misc.temporary_file.tmp_dir()
        sage: os.chdir(d)
        sage: with open("test.pyx", 'w') as f:
        ....:     _ = f.write("# distutils: language = c++\n"
        ....:       "from libcpp.vector cimport vector\n"
        ....:       "cdef vector[int] * v = new vector[int](4)\n")
        sage: output = sage.misc.cython.cython("test.pyx", create_local_c_file=True)

    Sage used to automatically include various ``.pxi`` files. Since
    :trac:`22805`, we no longer do this. But we make sure to give a
    useful message in case the ``.pxi`` files were needed::

        sage: cython("sig_malloc(0)")
        Traceback (most recent call last):
        ...
        RuntimeError: Error converting ... to C
        NOTE: Sage no longer automatically includes the deprecated files
        "cdefs.pxi", "signals.pxi" and "stdsage.pxi" in Cython files.
        You can fix your code by adding "from cysignals.memory cimport sig_malloc".
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
        base = sanitize(base)
    else:
        base = sanitize(os.path.abspath(filename))

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
                if os.path.getmtime(filename) <= os.path.getmtime('%s/%s'%(target_dir, prev_so)):
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
        name = '%s_%s'%(base, sequence_number[base])

        # increment the sequence number so will use a different one next time.
        sequence_number[base] += 1

    if compile_message:
        print("Compiling {}...".format(filename), file=sys.stderr)

    with open(filename) as f:
        (preparsed, libs, includes, language, additional_source_files,
         extra_args, libdirs) = _pyx_preparse(f.read())

    # New filename with preparsed code.
    # NOTE: if we ever stop preparsing, we should still copy the
    # original file to the target directory.
    pyxfile = os.path.join(target_dir, name + ".pyx")
    with open(pyxfile, 'w') as f:
        f.write(preparsed)

    extra_sources = []
    for fname in additional_source_files:
        fname = fname.replace("$SAGE_SRC", SAGE_SRC)
        fname = fname.replace("$SAGE_LOCAL", SAGE_LOCAL)
        extra_sources.append(fname)

    # Now do the actual build, directly calling Cython and distutils
    from Cython.Build import cythonize
    from Cython.Compiler.Errors import CompileError
    import Cython.Compiler.Options
    from distutils.dist import Distribution
    from distutils.core import Extension
    from distutils.log import set_verbosity
    set_verbosity(verbose)

    Cython.Compiler.Options.annotate = annotate
    Cython.Compiler.Options.embed_pos_in_docstring = True
    Cython.Compiler.Options.pre_import = "sage.all" if sage_namespace else None

    ext = Extension(name,
                    sources=[pyxfile] + extra_sources,
                    libraries=libs,
                    library_dirs=[os.path.join(SAGE_LOCAL, "lib")] + libdirs,
                    extra_compile_args=extra_args,
                    language=language)

    orig_cwd = os.getcwd()
    try:
        # Change directories to target_dir so that Cython produces the correct
        # relative path; https://trac.sagemath.org/ticket/24097
        os.chdir(target_dir)
        ext, = cythonize([ext],
                         aliases=cython_aliases(),
                         include_path=includes,
                         quiet=not verbose)
    except CompileError:
        # Check for names in old_pxi_names
        note = ''
        for pxd, names in old_pxi_names.items():
            for name in names:
                if re.search(r"\b{}\b".format(name), preparsed):
                    note += dedent(
                        """
                        NOTE: Sage no longer automatically includes the deprecated files
                        "cdefs.pxi", "signals.pxi" and "stdsage.pxi" in Cython files.
                        You can fix your code by adding "from {} cimport {}".
                        """.format(pxd, name))
        raise RuntimeError("Error converting {} to C".format(filename) +
                           note)
    finally:
        os.chdir(orig_cwd)

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
    dist.run_command("build")

    if create_local_so_file:
        # Copy module to current directory
        from sage.misc.sageinspect import loadable_module_extension
        shutil.copy(os.path.join(target_dir, name + loadable_module_extension()),
                    os.curdir)

    return name, target_dir


def subtract_from_line_numbers(s, n):
    r"""
    Given a string ``s`` and an integer ``n``, for any line of ``s`` which has
    the form ``'text:NUM:text'`` subtract ``n`` from NUM and return
    ``'text:(NUM-n):text'``. Return other lines of ``s`` without change.

    EXAMPLES::

        sage: from sage.misc.cython import subtract_from_line_numbers
        sage: subtract_from_line_numbers('hello:1234:hello', 3)
        doctest:...: DeprecationWarning: subtract_from_line_numbers is deprecated
        See http://trac.sagemath.org/22805 for details.
        'hello:1231:hello\n'
        sage: subtract_from_line_numbers('text:123\nhello:1234:', 3)
        'text:123\nhello:1231:\n'
    """
    from sage.misc.superseded import deprecation
    deprecation(22805, 'subtract_from_line_numbers is deprecated')

    ans = []
    for X in s.split('\n'):
        i = X.find(':')
        j = i+1 + X[i+1:].find(':')
        try:
            ans.append('%s:%s:%s\n'%(X[:i], int(X[i+1:j]) - n, X[j+1:]))
        except ValueError:
            ans.append(X)
    return '\n'.join(ans)


################################################################
# COMPILE
################################################################
def cython_lambda(vars, expr,
                 verbose=False,
                 compile_message=False,
                 use_cache=False):
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
        v = ', '.join(['%s %s'%(typ,var) for typ, var in vars])

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
    """%(v, expr)
    if verbose:
        print(s)
    tmpfile = tmp_filename(ext=".spyx")
    open(tmpfile,'w').write(s)

    d = {}
    cython_import_all(tmpfile, d,
                      verbose=verbose, compile_message=compile_message,
                      use_cache=use_cache,
                      create_local_c_file=False)
    return d['f']


def cython_create_local_so(filename):
    r"""
    Compile filename and make it available as a loadable shared object file.

    INPUT:

    - ``filename`` - string: a Cython (.spyx) file

    OUTPUT: None

    EFFECT: A compiled, python "importable" loadable shared object file is created.

    .. note::

        Shared object files are *not* reloadable. The intent is for
        imports in other scripts. A possible development cycle might
        go thus:

        - Attach a .spyx file
        - Interactively test and edit it to your satisfaction
        - Use ``cython_create_local_so`` to create the shared object file
        - Import the .so file in other scripts

    EXAMPLES::

        sage: curdir = os.path.abspath(os.curdir)
        sage: dir = tmp_dir(); os.chdir(dir)
        sage: f = open('hello.spyx', 'w')
        sage: s = "def hello():\n    print('hello')\n"
        sage: _ = f.write(s)
        sage: f.close()
        sage: cython_create_local_so('hello.spyx')
        Compiling hello.spyx...
        sage: sys.path.append('.')
        sage: import hello
        sage: hello.hello()
        hello
        sage: os.chdir(curdir)

    AUTHORS:

    - David Fu (2008-04-09): initial version
    """
    cython(filename, compile_message=True, use_cache=False, create_local_so_file=True)


################################################################
# IMPORT
################################################################
def cython_import(filename, verbose=False, compile_message=False,
                 use_cache=False, create_local_c_file=True, **kwds):
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
    name, build_dir = cython(filename, verbose=verbose,
                             compile_message=compile_message,
                             use_cache=use_cache,
                             create_local_c_file=create_local_c_file,
                             **kwds)
    sys.path.append(build_dir)
    return builtins.__import__(name)


def cython_import_all(filename, globals, verbose=False, compile_message=False,
                     use_cache=False, create_local_c_file=True):
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
    m = cython_import(filename, verbose=verbose, compile_message=compile_message,
                     use_cache=use_cache,
                     create_local_c_file=create_local_c_file)
    for k, x in iteritems(m.__dict__):
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


def compile_and_load(code):
    r"""
    INPUT:

    - ``code`` -- string containing code that could be in a .pyx file
      that is attached or put in a %cython block in the notebook.

    OUTPUT: a module, which results from compiling the given code and
    importing it

    EXAMPLES::

        sage: module = sage.misc.cython.compile_and_load("def f(int n):\n    return n*n")
        sage: module.f(10)
        100
    """
    file = tmp_filename(ext=".pyx")
    open(file,'w').write(code)
    return cython_import(file, create_local_c_file=False)


def cython_compile(code,
          verbose=False, compile_message=False,
          make_c_file_nice=False, use_cache=False):
    """
    Given a block of Cython code (as a text string), this function
    compiles it using a C compiler, and includes it into the global
    namespace.

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
    tmpfile = tmp_filename(ext=".pyx")
    open(tmpfile,'w').write(code)
    cython_import_all(tmpfile, get_globals(),
                      verbose=verbose, compile_message=compile_message,
                      use_cache=use_cache,
                      create_local_c_file=False)


TESTS = {
'trac11680':"""
from sage.rings.rational cimport Rational
from sage.rings.polynomial.polynomial_rational_flint cimport Polynomial_rational_flint
from sage.libs.flint.fmpq_poly cimport fmpq_poly_length, fmpq_poly_get_coeff_mpq, fmpq_poly_set_coeff_mpq

def evaluate_at_power_of_gen(Polynomial_rational_flint f, unsigned long n):
    assert n >= 1
    cdef Polynomial_rational_flint res = f._new()
    cdef unsigned long k
    cdef Rational z = Rational(0)
    for k in range(fmpq_poly_length(f.__poly)):
        fmpq_poly_get_coeff_mpq(z.value, f.__poly, k)
        fmpq_poly_set_coeff_mpq(res.__poly, n*k, z.value)
    return res
""",

'trac11680b':"""
def f(int a, int b, int c):
    return a+b+c
"""
}

def import_test(name):
    """
    This is used by the testing infrastructure to test building
    Cython programs.

    INPUT:

    - ``name`` -- string; name of a key to the TESTS dictionary above

    OUTPUT: a module, which results from compiling the given code and importing it

    EXAMPLES::

        sage: module = sage.misc.cython.import_test("trac11680b")
        sage: module.f(2,3,4)
        9
    """
    return compile_and_load(TESTS[name])

