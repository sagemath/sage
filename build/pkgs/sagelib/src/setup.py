#!/usr/bin/env python

from __future__ import print_function

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

import sage.env
sage.env.SAGE_SRC = os.getcwd()
from sage.env import *

from sage_setup.excepthook import excepthook
sys.excepthook = excepthook

#########################################################
### Configuration
#########################################################

if len(sys.argv) > 1 and sys.argv[1] == "sdist":
    sdist = True
else:
    sdist = False

#########################################################
### Testing related stuff
#########################################################

# Remove (potentially invalid) star import caches
import sage.misc.lazy_import_cache
if os.path.exists(sage.misc.lazy_import_cache.get_cache_file()):
    os.unlink(sage.misc.lazy_import_cache.get_cache_file())


from sage_setup.command.sage_build import sage_build
from sage_setup.command.sage_build_cython import sage_build_cython
from sage_setup.command.sage_build_ext import sage_build_ext


#########################################################
### Discovering Sources
#########################################################

# TODO: This should be quiet by default
print("Discovering Python/Cython source code....")
t = time.time()

from sage_setup.optional_extension import is_package_installed_and_updated

if sdist:
    # No need to compute distributions.  This avoids a dependency on Cython
    # just to make an sdist.
    distributions = None
else:
    distributions = ['']
    optional_packages_with_extensions = ['mcqd', 'bliss', 'tdlib', 'primecount',
                                         'coxeter3', 'fes', 'sirocco', 'meataxe']
    distributions += ['sage-{}'.format(pkg)
                      for pkg in optional_packages_with_extensions
                      if is_package_installed_and_updated(pkg)]
    log.warn('distributions = {0}'.format(distributions))

from sage_setup.find import find_python_sources
python_packages, python_modules, cython_modules = find_python_sources(
    SAGE_SRC, ['sage', 'sage_setup'], distributions=distributions)

log.debug('python_packages = {0}'.format(python_packages))

print("Discovered Python/Cython sources, time: %.2f seconds." % (time.time() - t))


from sage_setup.command.sage_install import sage_install_and_clean

#########################################################
### Distutils
#########################################################

code = setup(
      packages = python_packages,
      package_data = {
          'sage.libs.gap': ['sage.gaprc'],
          'sage.interfaces': ['sage-maxima.lisp'],
          'sage.doctest':  ['tests/*'],
          'sage': ['ext_data/*',
                   'ext_data/kenzo/*',
                   'ext_data/singular/*',
                   'ext_data/singular/function_field/*',
                   'ext_data/images/*',
                   'ext_data/doctest/*',
                   'ext_data/doctest/invalid/*',
                   'ext_data/doctest/rich_output/*',
                   'ext_data/doctest/rich_output/example_wavefront/*',
                   'ext_data/gap/*',
                   'ext_data/gap/joyner/*',
                   'ext_data/mwrank/*',
                   'ext_data/notebook-ipython/*',
                   'ext_data/nbconvert/*',
                   'ext_data/graphs/*',
                   'ext_data/pari/*',
                   'ext_data/pari/dokchitser/*',
                   'ext_data/pari/buzzard/*',
                   'ext_data/pari/simon/*',
                   'ext_data/magma/*',
                   'ext_data/magma/latex/*',
                   'ext_data/magma/sage/*',
                   'ext_data/valgrind/*',
                   'ext_data/threejs/*']
      },
      scripts = [## The sage script
                 'bin/sage',
                 ## Other scripts that should be in the path also for OS packaging of sage:
                 'bin/sage-eval',
                 'bin/sage-runtests',          # because it is useful for doctesting user scripts too
                 'bin/sage-fixdoctests',       # likewise
                 'bin/sage-coverage',          # because it is useful for coverage-testing user scripts too
                 'bin/sage-cython',            # deprecated, might be used in user package install scripts
                 ## Helper scripts invoked by sage script
                 ## (they would actually belong to something like libexec)
                 'bin/sage-cachegrind',
                 'bin/sage-callgrind',
                 'bin/sage-massif',
                 'bin/sage-omega',
                 'bin/sage-valgrind',
                 'bin/sage-venv-config',
                 'bin/sage-version.sh',
                 'bin/sage-cleaner',
                 ## Only makes sense in sage-the-distribution. TODO: Move to another installation script.
                 'bin/sage-list-packages',
                 'bin/sage-location',
                 ## Uncategorized scripts in alphabetical order
                 'bin/math-readline',
                 'bin/sage-env',
                 # sage-env-config -- installed by sage_conf
                 # sage-env-config.in -- not to be installed',
                 'bin/sage-gdb-commands',
                 'bin/sage-grep',
                 'bin/sage-grepdoc',
                 'bin/sage-inline-fortran',
                 'bin/sage-ipynb2rst',
                 'bin/sage-ipython',
                 'bin/sage-native-execute',
                 'bin/sage-notebook',
                 'bin/sage-num-threads.py',
                 'bin/sage-open',
                 'bin/sage-preparse',
                 'bin/sage-python',
                 'bin/sage-rebase.bat',
                 'bin/sage-rebase.sh',
                 'bin/sage-rebaseall.bat',
                 'bin/sage-rebaseall.sh',
                 'bin/sage-rst2txt',
                 'bin/sage-run',
                 'bin/sage-run-cython',
                 'bin/sage-startuptime.py',
                 'bin/sage-update-src',
                 'bin/sage-update-version',
                 ],
      cmdclass = dict(build=sage_build,
                      build_cython=sage_build_cython,
                      build_ext=sage_build_ext,
                      install=sage_install_and_clean),
      ext_modules = cython_modules)
