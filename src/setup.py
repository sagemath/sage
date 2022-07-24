#!/usr/bin/env python

## This version of setup.py is used by the Sage distribution
## only when configure --enable-editable has been used.
##
## Distribution packaging should use build/pkgs/sagelib/src/setup.py
## instead.

import os
import platform
import sys
import time
from setuptools import setup, find_namespace_packages
from distutils import log
import multiprocessing.pool

# PEP 517 builds do not have . in sys.path
sys.path.insert(0, os.path.dirname(__file__))

import sage.misc.lazy_import_cache

from sage.misc.package import is_package_installed_and_updated
from sage_setup.command.sage_build_ext_minimal import sage_build_ext_minimal
from sage_setup.command.sage_install import sage_develop, sage_install
from sage_setup.find import filter_cython_sources
from sage_setup.cython_options import compiler_directives, compile_time_env_variables
from sage_setup.extensions import create_extension
from sage_setup.excepthook import excepthook

# Work around a Cython problem in Python 3.8.x on macOS
# https://github.com/cython/cython/issues/3262
if platform.system() == 'Darwin':
    import multiprocessing
    multiprocessing.set_start_method('fork', force=True)

# ########################################################
# ## Set source directory
# ########################################################

import sage.env
sage.env.SAGE_SRC = os.getcwd()
from sage.env import *

sys.excepthook = excepthook

from sage_setup.setenv import setenv
setenv()

# ########################################################
# ## Configuration
# ########################################################

if len(sys.argv) > 1 and (sys.argv[1] == "sdist" or sys.argv[1] == "egg_info"):
    sdist = True
else:
    sdist = False

# ########################################################
# ## Testing related stuff
# ########################################################

# Remove (potentially invalid) star import caches
if os.path.exists(sage.misc.lazy_import_cache.get_cache_file()):
    os.unlink(sage.misc.lazy_import_cache.get_cache_file())


# ########################################################
# ## Discovering Sources
# ########################################################
if sdist:
    extensions = []
    python_packages = []
else:
    log.info("Generating auto-generated sources")
    from sage_setup.autogen import autogen_all
    autogen_all()

    log.info("Discovering Python/Cython source code....")
    t = time.time()

    # Exclude a few files if the corresponding distribution is not loaded
    optional_packages = ['mcqd', 'bliss', 'tdlib',
                         'coxeter3', 'fes', 'sirocco', 'meataxe']
    not_installed_packages = [package for package in optional_packages
                              if not is_package_installed_and_updated(package)]

    distributions_to_exclude = [f"sagemath-{pkg}"
                                for pkg in not_installed_packages]
    files_to_exclude = filter_cython_sources(SAGE_SRC, distributions_to_exclude)

    log.debug(f"files_to_exclude = {files_to_exclude}")

    python_packages = find_namespace_packages(where=SAGE_SRC, include=['sage', 'sage.*'])
    log.debug(f"python_packages = {python_packages}")

    log.info(f"Discovered Python/Cython sources, time: {(time.time() - t):.2f} seconds.")

    # from sage_build_cython:
    import Cython.Compiler.Options
    Cython.Compiler.Options.embed_pos_in_docstring = True
    gdb_debug = os.environ.get('SAGE_DEBUG', None) != 'no'

    try:
        from Cython.Build import cythonize
        from sage.env import cython_aliases, sage_include_directories
        from sage.misc.package_dir import cython_namespace_package_support
        with cython_namespace_package_support():
            extensions = cythonize(
                ["sage/**/*.pyx"],
                exclude=files_to_exclude,
                include_path=sage_include_directories(use_sources=True) + ['.'],
                compile_time_env=compile_time_env_variables(),
                compiler_directives=compiler_directives(False),
                aliases=cython_aliases(),
                create_extension=create_extension,
                gdb_debug=gdb_debug,
                nthreads=4)
    except Exception as exception:
        log.warn(f"Exception while cythonizing source files: {repr(exception)}")
        raise

# ########################################################
# ## Distutils
# ########################################################
code = setup(
    packages=python_packages,
    cmdclass={
        "build_ext": sage_build_ext_minimal,
        "develop":   sage_develop,
        "install":   sage_install,
    },
    ext_modules=extensions
)
