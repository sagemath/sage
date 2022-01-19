#!/usr/bin/env python

import os
import sys
import time
# Import setuptools before importing distutils, so that setuptools
# can replace distutils by its own vendored copy.
import setuptools
from distutils import log
from setuptools import setup

# Work around a Cython problem in Python 3.8.x on macOS
# https://github.com/cython/cython/issues/3262
if os.uname().sysname == 'Darwin':
    import multiprocessing
    multiprocessing.set_start_method('fork', force=True)

#########################################################
### Set source directory
#########################################################

# PEP 517 builds do not have . in sys.path
sys.path.insert(0, os.path.dirname(__file__))

import sage.env
sage.env.SAGE_SRC = os.getcwd()
from sage.env import *

#########################################################
### Configuration
#########################################################

if len(sys.argv) > 1 and (sys.argv[1] == "sdist" or sys.argv[1] == "egg_info"):
    sdist = True
else:
    sdist = False

if sdist:
    cmdclass = {}
else:
    from sage_setup.excepthook import excepthook
    sys.excepthook = excepthook

    from sage_setup.setenv import setenv
    setenv()

    from sage_setup.command.sage_build import sage_build
    from sage_setup.command.sage_build_cython import sage_build_cython
    from sage_setup.command.sage_build_ext import sage_build_ext
    from sage_setup.command.sage_install import sage_install_and_clean

    cmdclass = dict(build=sage_build,
                    build_cython=sage_build_cython,
                    build_ext=sage_build_ext,
                    install=sage_install_and_clean)

#########################################################
### Testing related stuff
#########################################################

# Remove (potentially invalid) star import caches
import sage.misc.lazy_import_cache
if os.path.exists(sage.misc.lazy_import_cache.get_cache_file()):
    os.unlink(sage.misc.lazy_import_cache.get_cache_file())

#########################################################
### Discovering Sources
#########################################################

if sdist:
    # No need to compute distributions.  This avoids a dependency on Cython
    # just to make an sdist.
    distributions = None
    python_packages = []
    python_modules = []
    cython_modules = []
else:
    # TODO: This should be quiet by default
    print("Discovering Python/Cython source code....")
    t = time.time()
    from sage_setup.optional_extension import is_package_installed_and_updated
    distributions = ['']
    optional_packages_with_extensions = ['mcqd', 'bliss', 'tdlib',
                                         'coxeter3', 'fes', 'sirocco', 'meataxe']
    distributions += ['sagemath-{}'.format(pkg)
                      for pkg in optional_packages_with_extensions
                      if is_package_installed_and_updated(pkg)]
    log.warn('distributions = {0}'.format(distributions))
    from sage_setup.find import find_python_sources
    python_packages, python_modules, cython_modules = find_python_sources(
        SAGE_SRC, ['sage'], distributions=distributions)

    log.debug('python_packages = {0}'.format(python_packages))
    print("Discovered Python/Cython sources, time: %.2f seconds." % (time.time() - t))


#########################################################
### Distutils
#########################################################

code = setup(
      packages = python_packages,
      cmdclass = cmdclass,
      ext_modules = cython_modules)
