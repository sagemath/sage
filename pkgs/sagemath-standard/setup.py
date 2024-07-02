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

# If build isolation is not in use and setuptools_scm is installed,
# then its file_finders entry point is invoked, which we don't need.
# Workaround from ​https://github.com/pypa/setuptools_scm/issues/190#issuecomment-351181286
try:
    import setuptools_scm.integration
    setuptools_scm.integration.find_files = lambda _: []
except ImportError:
    pass

# Different workaround: disable `walk_revctrl` in setuptools
# This is needed for setuptools_scm >= 8, should work for any version
import setuptools.command.egg_info
setuptools.command.egg_info.walk_revctrl = lambda: ()

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

from sage_setup.excepthook import excepthook
sys.excepthook = excepthook

from sage_setup.setenv import setenv
setenv()

from sage_setup.command.sage_build_cython import sage_build_cython
from sage_setup.command.sage_build_ext import sage_build_ext
from sage_setup.command.sage_install import sage_develop, sage_install

cmdclass = dict(build_cython=sage_build_cython,
                build_ext=sage_build_ext,
                develop=sage_develop,
                install=sage_install)

#########################################################
### Discovering Sources
#########################################################

if any(x in sys.argv
       for x in ['build', 'build_ext', 'bdist_wheel', 'install']):
    log.info("Generating auto-generated sources")
    from sage_setup.autogen import autogen_all
    autogen_all()

# TODO: This should be quiet by default
print("Discovering Python/Cython source code....")
t = time.time()
distributions = ['sagemath-categories',
                 'sagemath-environment',
                 'sagemath-objects',
                 'sagemath-repl',
                 '']
log.warn('distributions = {0}'.format(distributions))
from sage_setup.find import find_python_sources
python_packages, python_modules, cython_modules = find_python_sources(
    SAGE_SRC, ['sage'], distributions=distributions)

log.debug('python_packages = {0}'.format(python_packages))
log.debug('python_modules = {0}'.format(python_modules))
log.debug('cython_modules = {0}'.format(cython_modules))

print("Discovered Python/Cython sources, time: %.2f seconds." % (time.time() - t))

#########################################################
### Distutils
#########################################################

code = setup(
      packages=python_packages,
      cmdclass=cmdclass,
      ext_modules=cython_modules,
)
