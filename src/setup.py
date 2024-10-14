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
from setuptools.dist import Distribution
from distutils import log
import multiprocessing.pool

# PEP 517 builds do not have . in sys.path
sys.path.insert(0, os.path.dirname(__file__))

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

# setuptools plugins considered harmful:
# If build isolation is not in use and setuptools_scm is installed,
# then its file_finders entry point is invoked, which we don't need.
# And with setuptools_scm 8, we get more trouble:
# LookupError: pyproject.toml does not contain a tool.setuptools_scm section
# LookupError: setuptools-scm was unable to detect version ...
# We just remove all handling of "setuptools.finalize_distribution_options" entry points.
Distribution._removed = staticmethod(lambda ep: True)

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

if len(sys.argv) > 1 and (sys.argv[1] in ["sdist", "egg_info", "dist_info"]):
    sdist = True
else:
    sdist = False

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

    log.info("Discovering Python/Cython source code...")

    optional_packages = ['mcqd', 'bliss', 'tdlib',
                         'coxeter3', 'sirocco', 'meataxe', 'cmr']
    distributions_to_exclude = [f"sagemath-{pkg}"
                                for pkg in optional_packages]
    files_to_exclude = filter_cython_sources(SAGE_SRC, distributions_to_exclude)

    log.debug(f"files_to_exclude = {files_to_exclude}")

    python_packages = find_namespace_packages(where=SAGE_SRC, include=['sage', 'sage.*'])
    log.debug(f"python_packages = {python_packages}")

    log.info(f"Discovering Python/Cython source code... done")

    # from sage_build_cython:
    import Cython.Compiler.Options
    Cython.Compiler.Options.embed_pos_in_docstring = True
    gdb_debug = os.environ.get('SAGE_DEBUG', None) != 'no'

    aliases = cython_aliases()
    log.debug(f"aliases = {aliases}")
    include_path = sage_include_directories(use_sources=True) + ['.']
    log.debug(f"include_path = {include_path}")
    nthreads = sage_build_ext_minimal.get_default_number_build_jobs()
    log.info(f"Cythonizing with {nthreads} threads...")
    try:
        from Cython.Build import cythonize
        from sage.env import cython_aliases, sage_include_directories
        from sage.misc.package_dir import cython_namespace_package_support
        with cython_namespace_package_support():
            extensions = cythonize(
                ["sage/**/*.pyx"],
                exclude=files_to_exclude,
                include_path=include_path,
                compile_time_env=compile_time_env_variables(),
                compiler_directives=compiler_directives(False),
                aliases=aliases,
                create_extension=create_extension,
                gdb_debug=gdb_debug,
                nthreads=nthreads)
    except Exception as exception:
        log.warn(f"Exception while cythonizing source files: {repr(exception)}")
        raise
    log.info(f"Cythonizing with {nthreads} threads... done")

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
