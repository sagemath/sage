r"""
Sage Runtime Environment

AUTHORS:

- \R. Andrew Ohana (2012): Initial version.

Verify that importing ``sage.all`` works in Sage's Python without any ``SAGE_``
environment variables, and has the same ``SAGE_ROOT`` and ``SAGE_LOCAL``
(see also :trac:`29446`)::

    sage: env = {k:v for (k,v) in os.environ.items() if not k.startswith("SAGE_")}
    sage: from subprocess import check_output
    sage: cmd = "from sage.all import SAGE_ROOT, SAGE_LOCAL; print((SAGE_ROOT, SAGE_LOCAL))"
    sage: out = check_output([sys.executable, "-c", cmd], env=env).decode().strip()   # long time
    sage: out == repr((SAGE_ROOT, SAGE_LOCAL))                                        # long time
    True
"""

# ****************************************************************************
#       Copyright (C) 2013 R. Andrew Ohana <andrew.ohana@gmail.com>
#       Copyright (C) 2019 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from typing import List, Optional
import sage
import os
import socket
import sys
import sysconfig
from . import version
from pathlib import Path
import subprocess


# All variables set by var() appear in this SAGE_ENV dict
SAGE_ENV = dict()


def join(*args):
    """
    Join paths like ``os.path.join`` except that the result is ``None``
    if any of the components is ``None``.

    EXAMPLES::

        sage: from sage.env import join
        sage: print(join("hello", "world"))
        hello/world
        sage: print(join("hello", None))
        None
    """
    if any(a is None for a in args):
        return None
    return os.path.join(*args)


def var(key: str, *fallbacks: Optional[str], force: bool = False) -> Optional[str]:
    """
    Set ``SAGE_ENV[key]`` and return the value.

    If ``key`` is an environment variable, this is the value.
    Otherwise, the ``fallbacks`` are tried until one is found which
    is not ``None``. If the environment variable is not set and all
    fallbacks are ``None``, then the final value is ``None``.

    INPUT:

    - ``key`` -- string.

    - ``fallbacks`` -- tuple containing ``str`` or ``None`` values.

    - ``force`` -- boolean (optional, default is ``False``). If
      ``True``, skip the environment variable and only use the
      fallbacks.

    OUTPUT:

    The value of the environment variable or its fallbacks.

    EXAMPLES::

        sage: import os, sage.env
        sage: sage.env.SAGE_ENV = dict()
        sage: os.environ['SAGE_FOO'] = 'foo'
        sage: sage.env.var('SAGE_FOO', 'unused')
        'foo'
        sage: sage.env.SAGE_FOO
        'foo'
        sage: sage.env.SAGE_ENV['SAGE_FOO']
        'foo'

    If the environment variable does not exist, the fallbacks (if any)
    are used. In most typical uses, there is exactly one fallback::

        sage: _ = os.environ.pop('SAGE_BAR', None)  # ensure that SAGE_BAR does not exist
        sage: sage.env.var('SAGE_BAR', 'bar')
        'bar'
        sage: sage.env.SAGE_BAR
        'bar'
        sage: sage.env.SAGE_ENV['SAGE_BAR']
        'bar'

    Test multiple fallbacks::

        sage: sage.env.var('SAGE_BAR', None, 'yes', 'no')
        'yes'
        sage: sage.env.SAGE_BAR
        'yes'

    If all fallbacks are ``None``, the result is ``None``::

        sage: sage.env.var('SAGE_BAR')
        sage: print(sage.env.SAGE_BAR)
        None
        sage: sage.env.var('SAGE_BAR', None)
        sage: print(sage.env.SAGE_BAR)
        None

    Test the ``force`` keyword::

        sage: os.environ['SAGE_FOO'] = 'foo'
        sage: sage.env.var('SAGE_FOO', 'forced', force=True)
        'forced'
        sage: sage.env.SAGE_FOO
        'forced'
        sage: sage.env.var('SAGE_FOO', 'forced', force=False)
        'foo'
        sage: sage.env.SAGE_FOO
        'foo'
    """
    if force:
        value = None
    else:
        value = os.environ.get(key)
    if value is None:
        try:
            import sage_conf
            value = getattr(sage_conf, key, None)
        except ImportError:
            pass
    # Try all fallbacks in order as long as we don't have a value
    for f in fallbacks:
        if value is not None:
            break
        value = f
    SAGE_ENV[key] = value
    globals()[key] = value
    return value


# system info
UNAME = var("UNAME", os.uname()[0])
HOSTNAME = var("HOSTNAME", socket.gethostname())
LOCAL_IDENTIFIER = var("LOCAL_IDENTIFIER", "{}.{}".format(HOSTNAME, os.getpid()))

# version info
SAGE_VERSION = var("SAGE_VERSION", version.version)
SAGE_DATE = var("SAGE_DATE", version.date)
SAGE_VERSION_BANNER = var("SAGE_VERSION_BANNER", version.banner)

# virtual environment where sagelib is installed
SAGE_VENV = var("SAGE_VENV", os.path.abspath(sys.prefix))
SAGE_LIB = var("SAGE_LIB", os.path.dirname(os.path.dirname(sage.__file__)))
SAGE_EXTCODE = var("SAGE_EXTCODE", join(SAGE_LIB, "sage", "ext_data"))
SAGE_VENV_SPKG_INST = var("SAGE_VENV_SPKG_INST", join(SAGE_VENV, "var", "lib", "sage", "installed"))

# prefix hierarchy where non-Python packages are installed
SAGE_LOCAL = var("SAGE_LOCAL", SAGE_VENV)
SAGE_SHARE = var("SAGE_SHARE", join(SAGE_LOCAL, "share"))
SAGE_DOC = var("SAGE_DOC", join(SAGE_SHARE, "doc", "sage"))
SAGE_SPKG_INST = var("SAGE_SPKG_INST", join(SAGE_LOCAL, "var", "lib", "sage", "installed"))

# source tree of the Sage distribution
SAGE_ROOT = var("SAGE_ROOT")  # no fallback for SAGE_ROOT
SAGE_SRC = var("SAGE_SRC", join(SAGE_ROOT, "src"), SAGE_LIB)
SAGE_DOC_SRC = var("SAGE_DOC_SRC", join(SAGE_ROOT, "src", "doc"), SAGE_DOC)
SAGE_PKGS = var("SAGE_PKGS", join(SAGE_ROOT, "build", "pkgs"))
SAGE_ROOT_GIT = var("SAGE_ROOT_GIT", join(SAGE_ROOT, ".git"))

# ~/.sage
DOT_SAGE = var("DOT_SAGE", join(os.environ.get("HOME"), ".sage"))
SAGE_STARTUP_FILE = var("SAGE_STARTUP_FILE", join(DOT_SAGE, "init.sage"))

# for sage_setup.setenv
SAGE_ARCHFLAGS = var("SAGE_ARCHFLAGS", "unset")
SAGE_PKG_CONFIG_PATH = var("SAGE_PKG_CONFIG_PATH")

# installation directories for various packages
CONWAY_POLYNOMIALS_DATA_DIR = var("CONWAY_POLYNOMIALS_DATA_DIR", join(SAGE_SHARE, "conway_polynomials"))
GRAPHS_DATA_DIR = var("GRAPHS_DATA_DIR", join(SAGE_SHARE, "graphs"))
ELLCURVE_DATA_DIR = var("ELLCURVE_DATA_DIR", join(SAGE_SHARE, "ellcurves"))
POLYTOPE_DATA_DIR = var("POLYTOPE_DATA_DIR", join(SAGE_SHARE, "reflexive_polytopes"))
GAP_ROOT_DIR = var("GAP_ROOT_DIR", join(SAGE_SHARE, "gap"))
THEBE_DIR = var("THEBE_DIR", join(SAGE_SHARE, "thebe"))
COMBINATORIAL_DESIGN_DATA_DIR = var("COMBINATORIAL_DESIGN_DATA_DIR", join(SAGE_SHARE, "combinatorial_designs"))
CREMONA_MINI_DATA_DIR = var("CREMONA_MINI_DATA_DIR", join(SAGE_SHARE, "cremona"))
CREMONA_LARGE_DATA_DIR = var("CREMONA_LARGE_DATA_DIR", join(SAGE_SHARE, "cremona"))
JMOL_DIR = var("JMOL_DIR", join(SAGE_SHARE, "jmol"))
MATHJAX_DIR = var("MATHJAX_DIR", join(SAGE_SHARE, "mathjax"))
MTXLIB = var("MTXLIB", join(SAGE_SHARE, "meataxe"))
THREEJS_DIR = var("THREEJS_DIR", join(SAGE_SHARE, "threejs-sage"))
PPLPY_DOCS = var("PPLPY_DOCS", join(SAGE_SHARE, "doc", "pplpy"))
MAXIMA = var("MAXIMA", "maxima")
MAXIMA_FAS = var("MAXIMA_FAS")
KENZO_FAS = var("KENZO_FAS")
SAGE_NAUTY_BINS_PREFIX = var("SAGE_NAUTY_BINS_PREFIX", "")
FOURTITWO_HILBERT = var("FOURTITWO_HILBERT")
FOURTITWO_MARKOV = var("FOURTITWO_MARKOV")
FOURTITWO_GRAVER = var("FOURTITWO_GRAVER")
FOURTITWO_ZSOLVE = var("FOURTITWO_ZSOLVE")
FOURTITWO_QSOLVE = var("FOURTITWO_QSOLVE")
FOURTITWO_RAYS = var("FOURTITWO_RAYS")
FOURTITWO_PPI = var("FOURTITWO_PPI")
FOURTITWO_CIRCUITS = var("FOURTITWO_CIRCUITS")
FOURTITWO_GROEBNER = var("FOURTITWO_GROEBNER")
ARB_LIBRARY = var("ARB_LIBRARY", "arb")
CBLAS_PC_MODULES = var("CBLAS_PC_MODULES", "cblas:openblas:blas")
ECL_CONFIG = var("ECL_CONFIG", "ecl-config")
NTL_INCDIR = var("NTL_INCDIR")
NTL_LIBDIR = var("NTL_LIBDIR")
LIE_INFO_DIR = var("LIE_INFO_DIR", join(SAGE_LOCAL, "lib", "LiE"))

# The path to libSingular, to be passed to dlopen(). This will
# typically be set to an absolute path in sage_conf, but the relative
# fallback path here works on systems where dlopen() searches the
# system's library locations.
LIBSINGULAR_PATH = var("LIBSINGULAR_PATH", "libSingular.so")

# OpenMP
OPENMP_CFLAGS = var("OPENMP_CFLAGS", "")
OPENMP_CXXFLAGS = var("OPENMP_CXXFLAGS", "")

# Make sure mpmath uses Sage types
os.environ['MPMATH_SAGE'] = '1'

# misc
SAGE_BANNER = var("SAGE_BANNER", "")
SAGE_IMPORTALL = var("SAGE_IMPORTALL", "yes")


def _get_shared_lib_path(*libnames: str) -> Optional[str]:
    """
    Return the full path to a shared library file installed in
    ``$SAGE_LOCAL/lib`` or the directories associated with the
    Python sysconfig.

    This can also be passed more than one library name (e.g. for cases where
    some library may have multiple names depending on the platform) in which
    case the first one found is returned.

    This supports most *NIX variants (in which ``lib<libname>.so`` is found
    under ``$SAGE_LOCAL/lib``), macOS (same, but with the ``.dylib``
    extension), and Cygwin (under ``$SAGE_LOCAL/bin/cyg<libname>.dll``,
    or ``$SAGE_LOCAL/bin/cyg<libname>-*.dll`` for versioned DLLs).

    For distributions like Debian that use a multiarch layout, we also try the
    multiarch lib paths (i.e. ``/usr/lib/<arch>/``).

    This returns ``None`` if no matching library file could be found.

    EXAMPLES::

        sage: from sage.env import _get_shared_lib_path
        sage: "gap" in _get_shared_lib_path("gap")
        True
        sage: _get_shared_lib_path("an_absurd_lib") is None
        True

    """

    for libname in libnames:
        search_directories: List[Path] = []
        patterns: List[str] = []
        if sys.platform == 'cygwin':
            # Later down we take the first matching DLL found, so search
            # SAGE_LOCAL first so that it takes precedence
            if SAGE_LOCAL:
                search_directories.append(Path(SAGE_LOCAL) / 'bin')
            search_directories.append(Path(sysconfig.get_config_var('BINDIR')))
            # Note: The following is not very robust, since if there are multible
            # versions for the same library this just selects one more or less
            # at arbitrary. However, practically speaking, on Cygwin, there
            # will only ever be one version
            patterns = [f'cyg{libname}.dll', f'cyg{libname}-*.dll']
        else:
            if sys.platform == 'darwin':
                ext = 'dylib'
            else:
                ext = 'so'

            if SAGE_LOCAL:
                search_directories.append(Path(SAGE_LOCAL) / 'lib')
            libdir = sysconfig.get_config_var('LIBDIR')
            if libdir is not None:
                libdir = Path(libdir)
                search_directories.append(libdir)

                multiarchlib = sysconfig.get_config_var('MULTIARCH')
                if multiarchlib is not None:
                    search_directories.append(libdir / multiarchlib),

            patterns = [f'lib{libname}.{ext}']

        for directory in search_directories:
            for pattern in patterns:
                path = next(directory.glob(pattern), None)
                if path is not None:
                    return str(path.resolve())

    # Just return None if no files were found
    return None

# locate libgap shared object
GAP_SO = var("GAP_SO", _get_shared_lib_path("gap", ""))

# post process
if DOT_SAGE is not None and ' ' in DOT_SAGE:
    if UNAME[:6] == 'CYGWIN':
        # on windows/cygwin it is typical for the home directory
        # to have a space in it.  Fortunately, users also have
        # write privileges to c:\cygwin\home, so we just put
        # .sage there.
        DOT_SAGE = var("DOT_SAGE", "/home/.sage", force=True)
    else:
        print("Your home directory has a space in it.  This")
        print("will probably break some functionality of Sage.  E.g.,")
        print("the GAP interface will not work. A workaround")
        print("is to set the environment variable HOME to a")
        print("directory with no spaces that you have write")
        print("permissions to before you start sage.")


def sage_include_directories(use_sources=False):
    """
    Return the list of include directories for compiling Sage extension modules.

    INPUT:

    -  ``use_sources`` -- (default: False) a boolean

    OUTPUT:

    a list of include directories to be used to compile sage code
    1. while building sage (use_sources='True')
    2. while using sage (use_sources='False')

    EXAMPLES:

    Expected output while using Sage::

        sage: import sage.env
        sage: sage.env.sage_include_directories()
        ['.../include/python...',
        '.../python.../numpy/core/include']

    To check that C/C++ files are correctly found, we verify that we can
    always find the include file ``sage/cpython/cython_metaclass.h``,
    with both values for ``use_sources``::

        sage: file = os.path.join("sage", "cpython", "cython_metaclass.h")
        sage: dirs = sage.env.sage_include_directories(use_sources=True)
        sage: any(os.path.isfile(os.path.join(d, file)) for d in dirs)
        True
        sage: dirs = sage.env.sage_include_directories(use_sources=False)
        sage: any(os.path.isfile(os.path.join(d, file)) for d in dirs)
        True
    """
    import distutils.sysconfig

    TOP = SAGE_SRC if use_sources else SAGE_LIB

    dirs = [TOP,
            distutils.sysconfig.get_python_inc()]
    try:
        import numpy
        dirs.append(numpy.get_include())
    except ModuleNotFoundError:
        pass
    return dirs

def get_cblas_pc_module_name() -> str:
    """
    Return the name of the BLAS libraries to be used.
    """
    import pkgconfig
    cblas_pc_modules = CBLAS_PC_MODULES.split(':')
    return next((blas_lib for blas_lib in cblas_pc_modules if pkgconfig.exists(blas_lib)))


default_required_modules = ('fflas-ffpack', 'givaro', 'gsl', 'linbox', 'Singular',
                            'libpng', 'gdlib', 'm4ri', 'zlib', 'cblas')


default_optional_modules = ('lapack',)


def cython_aliases(required_modules=None,
                   optional_modules=None):
    """
    Return the aliases for compiling Cython code. These aliases are
    macros which can occur in ``# distutils`` headers.

    INPUT:

    - ``required_modules`` -- (default: taken from ``default_required_modules``)
      iterable of ``str`` values.

    - ``optional_modules`` -- (default: taken from ``default_optional_modules``)
      iterable of ``str`` values.

    EXAMPLES::

        sage: from sage.env import cython_aliases
        sage: cython_aliases()
        {...}
        sage: sorted(cython_aliases().keys())
        ['ARB_LIBRARY',
         'CBLAS_CFLAGS',
         ...,
         'ZLIB_LIBRARIES']
        sage: cython_aliases(required_modules=('module-that-is-assumed-to-not-exist'))
        Traceback (most recent call last):
        ...
        PackageNotFoundError: ...
        sage: cython_aliases(required_modules=(), optional_modules=('module-that-is-assumed-to-not-exist'))
        {...}

    TESTS:

    We can use ``cython.parallel`` regardless of whether OpenMP is supported.
    This will run in parallel, if OpenMP is supported::

        sage: cython('''
        ....: #distutils: extra_compile_args = OPENMP_CFLAGS
        ....: #distutils: extra_link_args = OPENMP_CFLAGS
        ....: from cython.parallel import prange
        ....:
        ....: cdef int i
        ....: cdef int n = 30
        ....: cdef int sum = 0
        ....:
        ....: for i in prange(n, num_threads=4, nogil=True):
        ....:     sum += i
        ....:
        ....: print(sum)
        ....: ''')
        435
    """
    import pkgconfig
    import itertools

    if required_modules is None:
        required_modules = default_required_modules

    if optional_modules is None:
        optional_modules = default_optional_modules

    aliases = {}

    for lib, required in itertools.chain(((lib, True) for lib in required_modules),
                                         ((lib, False) for lib in optional_modules)):
        var = lib.upper().replace("-", "") + "_"
        if lib == 'cblas':
            lib = get_cblas_pc_module_name()
        if lib == 'zlib':
            aliases[var + "CFLAGS"] = ""
            try:
                pc = pkgconfig.parse('zlib')
                libs = pkgconfig.libs(lib)
            except pkgconfig.PackageNotFoundError:
                from collections import defaultdict
                pc = defaultdict(list, {'libraries': ['z']})
                libs = "-lz"
        else:
            try:
                aliases[var + "CFLAGS"] = pkgconfig.cflags(lib).split()
                pc = pkgconfig.parse(lib)
                libs = pkgconfig.libs(lib)
            except pkgconfig.PackageNotFoundError:
                if required:
                    raise
                else:
                    continue

        # It may seem that INCDIR is redundant because the -I options are also
        # passed in CFLAGS.  However, "extra_compile_args" are put at the end
        # of the compiler command line.  "include_dirs" go to the front; the
        # include search order matters.
        aliases[var + "INCDIR"] = pc['include_dirs']
        aliases[var + "LIBDIR"] = pc['library_dirs']
        aliases[var + "LIBEXTRA"] = list(filter(lambda s: not s.startswith(('-l','-L')), libs.split()))
        aliases[var + "LIBRARIES"] = pc['libraries']

    # uname-specific flags
    UNAME = os.uname()

    def uname_specific(name, value, alternative):
        if name in UNAME[0]:
            return value
        else:
            return alternative

    aliases["LINUX_NOEXECSTACK"] = uname_specific("Linux", ["-Wl,-z,noexecstack"],
                                                  [])
    aliases["CYGWIN_SQLITE3_LIBS"] = uname_specific("CYGWIN", ["sqlite3"],
                                                    [])

    # LinBox needs special care because it actually requires C++11 with
    # GNU extensions: -std=c++11 does not work, you need -std=gnu++11
    # (this is true at least with GCC 7.2.0).
    #
    # Further, note that LinBox does not add any C++11 flag in its .pc
    # file (possibly because of confusion between CFLAGS and CXXFLAGS?).
    # This is not a problem in practice since LinBox depends on
    # fflas-ffpack and fflas-ffpack does add such a C++11 flag.
    if "LINBOX_CFLAGS" in aliases:
        aliases["LINBOX_CFLAGS"].append("-std=gnu++11")

    aliases["ARB_LIBRARY"] = ARB_LIBRARY

    # TODO: Remove Cygwin hack by installing a suitable cblas.pc
    if os.path.exists('/usr/lib/libblas.dll.a'):
        aliases["CBLAS_LIBS"] = ['gslcblas']

    try:
        aliases["M4RI_CFLAGS"].remove("-pedantic")
    except (ValueError, KeyError):
        pass

    # Determine ecl-specific compiler arguments using the ecl-config script
    ecl_cflags = subprocess.run([ECL_CONFIG, "--cflags"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.split()
    aliases["ECL_CFLAGS"] = list(filter(lambda s: not s.startswith('-I'), ecl_cflags))
    aliases["ECL_INCDIR"] = list(map(lambda s: s[2:], filter(lambda s: s.startswith('-I'), ecl_cflags)))
    ecl_libs = subprocess.run([ECL_CONFIG, "--libs"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.split()
    aliases["ECL_LIBDIR"] = list(map(lambda s: s[2:], filter(lambda s: s.startswith('-L'), ecl_libs)))
    aliases["ECL_LIBRARIES"] = list(map(lambda s: s[2:], filter(lambda s: s.startswith('-l'), ecl_libs)))
    aliases["ECL_LIBEXTRA"] = list(filter(lambda s: not s.startswith(('-l','-L')), ecl_libs))

    # NTL
    aliases["NTL_CFLAGS"] = ['-std=c++11']
    aliases["NTL_INCDIR"] = [NTL_INCDIR] if NTL_INCDIR else []
    aliases["NTL_LIBDIR"] = [NTL_LIBDIR] if NTL_LIBDIR else []
    aliases["NTL_LIBRARIES"] = ['ntl']
    aliases["NTL_LIBEXTRA"] = []

    # OpenMP
    aliases["OPENMP_CFLAGS"] = OPENMP_CFLAGS.split()
    aliases["OPENMP_CXXFLAGS"] = OPENMP_CXXFLAGS.split()

    return aliases
