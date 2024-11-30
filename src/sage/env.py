# sage_setup: distribution = sagemath-environment
r"""
Sage Runtime Environment

Verify that importing ``sage.all`` works in Sage's Python without any ``SAGE_``
environment variables, and has the same ``SAGE_ROOT`` and ``SAGE_LOCAL``
(see also :issue:`29446`)::

    sage: env = {k:v for (k,v) in os.environ.items() if not k.startswith("SAGE_")}
    sage: from subprocess import check_output
    sage: module_name = "sage.all"   # hide .all import from the linter
    sage: cmd  = f"from {module_name} import SAGE_ROOT, SAGE_LOCAL;"
    sage: cmd +=  "from os.path import samefile;"
    sage: cmd += f"s1 = samefile(SAGE_ROOT, '{SAGE_ROOT}') if SAGE_ROOT else True;"
    sage: cmd += f"s2 = samefile(SAGE_LOCAL, '{SAGE_LOCAL}');"
    sage: cmd += "print(s1 and s2);"
    sage: out = check_output([sys.executable, "-c", cmd], env=env).decode().strip()   # long time
    sage: out == "True"                                                               # long time
    True

AUTHORS:

- \R. Andrew Ohana (2012): initial version
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

import os
import socket
import sys
import sysconfig
from typing import Optional

from . import version

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

    - ``key`` -- string

    - ``fallbacks`` -- tuple containing ``str`` or ``None`` values

    - ``force`` -- boolean (default: ``False``); if
      ``True``, skip the environment variable and only use the
      fallbacks

    OUTPUT: the value of the environment variable or its fallbacks

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
            try:
                import sage.config
                value = getattr(sage.config, key, None)
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
HOSTNAME = var("HOSTNAME", socket.gethostname())
LOCAL_IDENTIFIER = var("LOCAL_IDENTIFIER", "{}.{}".format(HOSTNAME, os.getpid()))

# version info
SAGE_VERSION = var("SAGE_VERSION", version.version)
SAGE_DATE = var("SAGE_DATE", version.date)
SAGE_VERSION_BANNER = var("SAGE_VERSION_BANNER", version.banner)

# virtual environment where sagelib is installed
SAGE_VENV = var("SAGE_VENV", os.path.abspath(sys.prefix))
SAGE_LIB = var("SAGE_LIB", os.path.dirname(os.path.dirname(__file__)))
SAGE_EXTCODE = var("SAGE_EXTCODE", join(SAGE_LIB, "sage", "ext_data"))
SAGE_VENV_SPKG_INST = var("SAGE_VENV_SPKG_INST", join(SAGE_VENV, "var", "lib", "sage", "installed"))

# prefix hierarchy where non-Python packages are installed
SAGE_LOCAL = var("SAGE_LOCAL", SAGE_VENV)
SAGE_SHARE = var("SAGE_SHARE", join(SAGE_LOCAL, "share"))
SAGE_DOC = var("SAGE_DOC", join(SAGE_SHARE, "doc", "sage"))
SAGE_LOCAL_SPKG_INST = var("SAGE_LOCAL_SPKG_INST", join(SAGE_LOCAL, "var", "lib", "sage", "installed"))
SAGE_SPKG_INST = var("SAGE_SPKG_INST", join(SAGE_LOCAL, "var", "lib", "sage", "installed"))  # deprecated

# source tree of the Sage distribution
SAGE_ROOT = var("SAGE_ROOT")  # no fallback for SAGE_ROOT
SAGE_SRC = var("SAGE_SRC", join(SAGE_ROOT, "src"), SAGE_LIB)
SAGE_DOC_SRC = var("SAGE_DOC_SRC", join(SAGE_ROOT, "src", "doc"), SAGE_DOC)
SAGE_PKGS = var("SAGE_PKGS", join(SAGE_ROOT, "build", "pkgs"))
SAGE_ROOT_GIT = var("SAGE_ROOT_GIT", join(SAGE_ROOT, ".git"))

# Sage doc server (local server with PORT if URL is not given)
SAGE_DOC_SERVER_URL = var("SAGE_DOC_SERVER_URL")
# The default port is 0 so that the system will assign a random unused port > 1024
SAGE_DOC_LOCAL_PORT = var("SAGE_DOC_LOCAL_PORT", "0")

# ~/.sage
DOT_SAGE = var("DOT_SAGE", join(os.environ.get("HOME"), ".sage"))
SAGE_STARTUP_FILE = var("SAGE_STARTUP_FILE", join(DOT_SAGE, "init.sage"))

# for sage_setup.setenv
SAGE_ARCHFLAGS = var("SAGE_ARCHFLAGS", "unset")
SAGE_PKG_CONFIG_PATH = var("SAGE_PKG_CONFIG_PATH")

# colon-separated search path for databases.
SAGE_DATA_PATH = var("SAGE_DATA_PATH",
                     os.pathsep.join(filter(None, [
                         join(DOT_SAGE, "db"),
                         join(SAGE_SHARE, "sagemath"),
                         SAGE_SHARE,
                         ])))

# database directories, the default is to search in SAGE_DATA_PATH
CREMONA_LARGE_DATA_DIR = var("CREMONA_LARGE_DATA_DIR")
CREMONA_MINI_DATA_DIR = var("CREMONA_MINI_DATA_DIR")
ELLCURVE_DATA_DIR = var("ELLCURVE_DATA_DIR")
GRAPHS_DATA_DIR = var("GRAPHS_DATA_DIR")
POLYTOPE_DATA_DIR = var("POLYTOPE_DATA_DIR")

# installation directories for various packages
JMOL_DIR = var("JMOL_DIR")
MATHJAX_DIR = var("MATHJAX_DIR", join(SAGE_SHARE, "mathjax"))
MTXLIB = var("MTXLIB", join(SAGE_SHARE, "meataxe"))
THREEJS_DIR = var("THREEJS_DIR")
PPLPY_DOCS = var("PPLPY_DOCS", join(SAGE_SHARE, "doc", "pplpy"))
MAXIMA = var("MAXIMA", "maxima")
MAXIMA_FAS = var("MAXIMA_FAS")
KENZO_FAS = var("KENZO_FAS")
SAGE_NAUTY_BINS_PREFIX = var("SAGE_NAUTY_BINS_PREFIX", "")
SAGE_ECMBIN = var("SAGE_ECMBIN", "ecm")
RUBIKS_BINS_PREFIX = var("RUBIKS_BINS_PREFIX", "")
FOURTITWO_HILBERT = var("FOURTITWO_HILBERT")
FOURTITWO_MARKOV = var("FOURTITWO_MARKOV")
FOURTITWO_GRAVER = var("FOURTITWO_GRAVER")
FOURTITWO_ZSOLVE = var("FOURTITWO_ZSOLVE")
FOURTITWO_QSOLVE = var("FOURTITWO_QSOLVE")
FOURTITWO_RAYS = var("FOURTITWO_RAYS")
FOURTITWO_PPI = var("FOURTITWO_PPI")
FOURTITWO_CIRCUITS = var("FOURTITWO_CIRCUITS")
FOURTITWO_GROEBNER = var("FOURTITWO_GROEBNER")
CBLAS_PC_MODULES = var("CBLAS_PC_MODULES", "cblas:openblas:blas")
ECL_CONFIG = var("ECL_CONFIG", "ecl-config")
NTL_INCDIR = var("NTL_INCDIR")
NTL_LIBDIR = var("NTL_LIBDIR")
LIE_INFO_DIR = var("LIE_INFO_DIR", join(SAGE_LOCAL, "lib", "LiE"))
SINGULAR_BIN = var("SINGULAR_BIN") or "Singular"

# OpenMP
OPENMP_CFLAGS = var("OPENMP_CFLAGS", "")
OPENMP_CXXFLAGS = var("OPENMP_CXXFLAGS", "")

# Make sure mpmath uses Sage types
os.environ['MPMATH_SAGE'] = '1'

# misc
SAGE_BANNER = var("SAGE_BANNER", "")
SAGE_IMPORTALL = var("SAGE_IMPORTALL", "yes")

# GAP memory and args

SAGE_GAP_MEMORY = var('SAGE_GAP_MEMORY', None)
SAGE_GAP_COMMAND = var('SAGE_GAP_COMMAND', None)

# The semicolon-separated search path for GAP packages. It is passed
# directly to GAP via the -l flag.
GAP_ROOT_PATHS = var("GAP_ROOT_PATHS",
                     ";".join([join(SAGE_LOCAL, "lib", "gap"),
                               join(SAGE_LOCAL, "share", "gap")]))

# post process
if DOT_SAGE is not None and ' ' in DOT_SAGE:
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

    - ``use_sources`` -- boolean (default: ``False``)

    OUTPUT:

    a list of include directories to be used to compile sage code
    1. while building sage (use_sources='True')
    2. while using sage (use_sources='False')

    EXAMPLES:

    Expected output while using Sage::

        sage: import sage.env
        sage: sage.env.sage_include_directories()
        ['...',
         '.../numpy/...core/include',
         '.../include/python...']

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
    if use_sources:
        dirs = [SAGE_SRC]
    else:
        import sage
        dirs = [os.path.dirname(directory)
                for directory in sage.__path__]
    try:
        import numpy
        dirs.append(numpy.get_include())
    except ModuleNotFoundError:
        pass

    dirs.append(sysconfig.get_config_var('INCLUDEPY'))

    return dirs
