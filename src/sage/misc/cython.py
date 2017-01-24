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

from __future__ import print_function
from __future__ import absolute_import

import os
import sys
import platform
from six.moves import builtins

from sage.env import SAGE_LOCAL, SAGE_SRC
from .misc import SPYX_TMP
from .temporary_file import tmp_filename
import pkgconfig


# CBLAS can be one of multiple implementations
cblas_pc = pkgconfig.parse('cblas')
cblas_libs = list(cblas_pc['libraries'])
cblas_library_dirs = list(cblas_pc['library_dirs'])
cblas_include_dirs = list(cblas_pc['include_dirs'])

# TODO: Remove Cygwin hack by installing a suitable cblas.pc
if os.path.exists('/usr/lib/libblas.dll.a'):
    cblas_libs = 'gslcblas'

standard_libs = [
    'mpfr', 'gmp', 'gmpxx', 'stdc++', 'pari', 'm', 
    'ec', 'gsl',
] + cblas_libs + [
    'ntl']


offset = 0

def parse_keywords(kwd, s):
    r"""
    Given a keyword ``kwd`` and a string ``s``, return a list of all arguments
    on the same line as that keyword in ``s``, as well as a new copy of ``s``
    in which each occurrence of ``kwd`` is in a comment. If a comment already
    occurs on the line containing ``kwd``, no words after the ``#`` are added
    to the list.

    EXAMPLES::

        sage: import sage.misc.cython
        sage: sage.misc.cython.parse_keywords('clib', " clib foo bar baz\n #cinclude bar\n")
        (['foo', 'bar', 'baz'], ' #clib foo bar baz\n #cinclude bar\n')

        sage: sage.misc.cython.parse_keywords('clib', "# qux clib foo bar baz\n #cinclude bar\n")
        (['foo', 'bar', 'baz'], '# qux clib foo bar baz\n #cinclude bar\n')
        sage: sage.misc.cython.parse_keywords('clib', "# clib foo bar # baz\n #cinclude bar\n")
        (['foo', 'bar'], '# clib foo bar # baz\n #cinclude bar\n')
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

    return v, s

def environ_parse(s):
    r"""
    Given a string s, find each substring of the form ``'\$ABC'``. If the
    environment variable :envvar:`$ABC` is set, replace ``'\$ABC'`` with its
    value and move on to the next such substring. If it is not set, stop
    parsing there.

    EXAMPLES::

        sage: from sage.misc.cython import environ_parse
        sage: environ_parse('$SAGE_LOCAL') == SAGE_LOCAL
        True
        sage: environ_parse('$THIS_IS_NOT_DEFINED_ANYWHERE')
        '$THIS_IS_NOT_DEFINED_ANYWHERE'
        sage: os.environ['DEFINE_THIS'] = 'hello'
        sage: environ_parse('$DEFINE_THIS/$THIS_IS_NOT_DEFINED_ANYWHERE/$DEFINE_THIS')
        'hello/$THIS_IS_NOT_DEFINED_ANYWHERE/$DEFINE_THIS'
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
    return environ_parse(s)

def pyx_preparse(s):
    r"""
    Preparse a pyx file:

    * include ``cdefs.pxi``, ``signals.pxi`` from ``cysignals``, ``stdsage.pxi``
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

        sage: from sage.misc.cython import pyx_preparse
        sage: pyx_preparse("")
        ('\ninclude "cysignals/signals.pxi"  # ctrl-c interrupt block support\ninclude "stdsage.pxi"\n\ninclude "cdefs.pxi"\n',
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
        sage: s, libs, inc, lang, f, args, libdirs = pyx_preparse("# clang c++\n #clib foo\n # cinclude bar\n")
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

        sage: s, libs, inc, lang, f, args, libdirs = pyx_preparse("# cargs -O3 -ggdb\n")
        sage: args
        ['-w', '-O2', '-O3', '-ggdb']

    TESTS::

        sage: module = sage.misc.cython.import_test("trac11680")  # long time (7s on sage.math, 2012)
        sage: R.<x> = QQ[]
        sage: module.evaluate_at_power_of_gen(x^3 + x - 7, 5)  # long time
        x^15 + x^5 - 7
    """
    from sage.env import sage_include_directories

    lang, s = parse_keywords('clang', s)
    if lang:
        lang = lang[0].lower() # this allows both C++ and c++
    else:
        lang = "c"

    v, s = parse_keywords('clib', s)
    libs = v + standard_libs

    additional_source_files, s = parse_keywords('cfile', s)

    v, s = parse_keywords('cinclude', s)
    inc = [environ_parse(x.replace('"','').replace("'","")) for x in v] + sage_include_directories()
    s = """\ninclude "cdefs.pxi"\n""" + s
    s = """\ninclude "cysignals/signals.pxi"  # ctrl-c interrupt block support\ninclude "stdsage.pxi"\n""" + s
    args, s = parse_keywords('cargs', s)
    args = ['-w','-O2'] + args
    libdirs = cblas_library_dirs

    # Add cysignals directory to includes
    for path in sys.path:
        cysignals_path = os.path.join(path, "cysignals")
        if os.path.isdir(cysignals_path):
            inc.append(cysignals_path)

    return s, libs, inc, lang, additional_source_files, args, libdirs

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

def cython(filename, verbose=False, compile_message=False,
           use_cache=False, create_local_c_file=False, annotate=True, sage_namespace=True,
           create_local_so_file=False):
    r"""
    Compile a Cython file. This converts a Cython file to a C (or C++ file),
    and then compiles that. The .c file and the .so file are
    created in a temporary directory.

    INPUT:

    - ``filename`` - the name of the file to be compiled. Should end with
      'pyx'.

    - ``verbose`` (bool, default False) - if True, print debugging
      information.

    - ``compile_message`` (bool, default False) - if True, print
      ``'Compiling <filename>...'`` to the standard error.

    - ``use_cache`` (bool, default False) - if True, check the
      temporary build directory to see if there is already a
      corresponding .so file. If so, and if the .so file is newer than the
      Cython file, don't recompile, just reuse the .so file.

    - ``create_local_c_file`` (bool, default False) - if True, save a
      copy of the .c file in the current directory.

    - ``annotate`` (bool, default True) - if True, create an html file which
      annotates the conversion from .pyx to .c. By default this is only created
      in the temporary directory, but if ``create_local_c_file`` is also True,
      then save a copy of the .html file in the current directory.

    - ``sage_namespace`` (bool, default True) - if True, import
      ``sage.all``.

    - ``create_local_so_file`` (bool, default False) - if True, save a
      copy of the compiled .so file in the current directory.

    TESTS:

    Before :trac:`12975`, it would have been needed to write ``#clang c++``,
    but upper case ``C++`` has resulted in an error::

        sage: code = [
        ....: "#clang C++",
        ....: "#clib m readline Singular givaro ntl gmpxx gmp",
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

    Check that compiling c++ code works::

        sage: cython("#clang C++\n"+
        ....:        "from libcpp.vector cimport vector\n"
        ....:        "cdef vector[int] * v = new vector[int](4)\n")
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

    # This is the *temporary* directory where we build the pyx file.
    # This is deleted when sage exits, which means pyx files must be
    # rebuilt every time Sage is restarted at present.
    build_dir = os.path.join(SPYX_TMP, base)

    if os.path.exists(build_dir):
        # There is already a module here. Maybe we do not have to rebuild?
        # Find the name.
        if use_cache:
            from sage.misc.sageinspect import loadable_module_extension
            prev_so = [F for F in os.listdir(build_dir) if F.endswith(loadable_module_extension())]
            if len(prev_so) > 0:
                prev_so = prev_so[0]     # should have length 1 because of deletes below
                if os.path.getmtime(filename) <= os.path.getmtime('%s/%s'%(build_dir, prev_so)):
                    # We do not have to rebuild.
                    return prev_so[:-len(loadable_module_extension())], build_dir
    else:
        os.makedirs(build_dir)
    for F in os.listdir(build_dir):
        G = '%s/%s'%(build_dir,F)
        try:
            if not os.path.isdir(G):
                os.unlink(G)
        except OSError:
            pass

    # Get the absolute path to the directory that contains the pyx file.
    # We will use this only to make some convenient symbolic links.
    abs_base = os.path.split(os.path.abspath(filename))[0]

    # bad things happen if the current directory is SAGE_SRC
    if not os.path.exists("%s/sage" % abs_base):
        cmd = 'cd "%s"; ln -sf "%s"/* .'%(build_dir, abs_base)
        os.system(cmd)
        if os.path.exists("%s/setup.py" % build_dir):
            os.unlink("%s/setup.py" % build_dir)

    if compile_message:
        print("Compiling {}...".format(filename), file=sys.stderr)

    F = open(filename).read()

    F, libs, includes, language, additional_source_files, extra_args, libdirs = pyx_preparse(F)

    # add the working directory to the includes so custom headers etc. work
    includes.append(os.path.split(os.path.splitext(filename)[0])[0])

    if language == 'c++':
        extension = "cpp"
    else:
        extension = "c"

    if create_local_so_file:
        name = base
    else:
        global sequence_number
        if base not in sequence_number:
            sequence_number[base] = 0
        name = '%s_%s'%(base, sequence_number[base])

        # increment the sequence number so will use a different one next time.
        sequence_number[base] += 1

    file_list = []
    for fname in additional_source_files:
        fname = fname.replace("$SAGE_SRC", SAGE_SRC)
        fname = fname.replace("$SAGE_LOCAL", SAGE_LOCAL)
        if fname.startswith(os.path.sep):
            file_list.append("'"+fname+"'")
        else:
            file_list.append("'"+os.path.abspath(os.curdir)+"/"+fname+"'")
    additional_source_files = ",".join(file_list)

    pyx = '%s/%s.pyx'%(build_dir, name)
    open(pyx,'w').write(F)
    setup="""
# Build using 'python setup.py'
import distutils.sysconfig, os, sys
from distutils.core import setup, Extension

from sage.env import SAGE_LOCAL

extra_compile_args = %s

ext_modules = [Extension('%s', sources=['%s.%s', %s],
                     libraries=%s,
                     library_dirs=[SAGE_LOCAL + '/lib/'] + %s,
                     extra_compile_args = extra_compile_args,
                     language = '%s' )]

setup(ext_modules = ext_modules,
      include_dirs = %s)
    """%(extra_args, name, name, extension, additional_source_files, libs, libdirs, language, includes)
    open('%s/setup.py'%build_dir,'w').write(setup)
    cython_include = ' '.join(["-I '%s'"%x for x in includes if len(x.strip()) > 0 ])

    options = ['-p']
    if annotate:
        options.append('-a')
    if sage_namespace:
        options.append('--pre-import sage.all')

    cmd = "cd '{DIR}' && cython {OPT} {INC} {LANG} '{NAME}.pyx' 1>log 2>err ".format(
        DIR=build_dir,
        OPT=' '.join(options),
        INC=cython_include,
        LANG='--cplus' if language=='c++' else '',
        NAME=name)

    if create_local_c_file:
        target_c = '%s/_%s.c'%(os.path.abspath(os.curdir), base)
        if language == 'c++':
            target_c = target_c + "pp"
        cmd += " && cp '%s.c' '%s'"%(name, target_c)
        if annotate:
            target_html = '%s/_%s.html'%(os.path.abspath(os.curdir), base)
            cmd += " && cp '%s.html' '%s'"%(name, target_html)

    if verbose:
        print(cmd)
    if os.system(cmd):
        log = open('%s/log'%build_dir).read()
        err = subtract_from_line_numbers(open('%s/err'%build_dir).read(), offset)
        raise RuntimeError("Error converting {} to C:\n{}\n{}".format(filename, log, err))

    cmd = 'cd %s && python setup.py build 1>log 2>err'%build_dir
    if verbose:
        print(cmd)
    if os.system(cmd):
        log = open('%s/log'%build_dir).read()
        err = open('%s/err'%build_dir).read()
        raise RuntimeError("Error compiling {}:\n{}\n{}".format(filename, log, err))

    # Move from lib directory.
    cmd = 'mv %s/build/lib.*/* %s'%(build_dir, build_dir)
    if verbose:
        print(cmd)
    if os.system(cmd):
        raise RuntimeError("Error copying extension module for {}".format(filename))

    if create_local_so_file:
        # Copy from lib directory into local directory
        from sage.misc.sageinspect import loadable_module_extension
        cmd = 'cp %s/%s%s %s'%(build_dir, name, loadable_module_extension(), os.path.abspath(os.curdir))
        if os.system(cmd):
            raise RuntimeError("Error making local copy of shared object library for {}".format(filename))

    return name, build_dir



def subtract_from_line_numbers(s, n):
    r"""
    Given a string ``s`` and an integer ``n``, for any line of ``s`` which has
    the form ``'text:NUM:text'`` subtract ``n`` from NUM and return
    ``'text:(NUM-n):text'``. Return other lines of ``s`` without change.

    EXAMPLES::

        sage: from sage.misc.cython import subtract_from_line_numbers
        sage: subtract_from_line_numbers('hello:1234:hello', 3)
        'hello:1231:hello\n'
        sage: subtract_from_line_numbers('text:123\nhello:1234:', 3)
        'text:123\nhello:1231:\n'
    """
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
        sage: f.write(s)
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
    for k, x in m.__dict__.iteritems():
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


TESTS = {
'trac11680':"""
#cargs -O3 -ggdb
#clib flint

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

