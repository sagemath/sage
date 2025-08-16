#!/usr/bin/env python

import logging
import os
import sys
from setuptools import setup

logging.basicConfig(
    level=logging.INFO,
    format='%(message)s'  # keep as distutils.log the same
)
log = logging.getLogger(__name__)

# Work around a Cython problem in Python 3.8.x on macOS
# https://github.com/cython/cython/issues/3262
if os.uname().sysname == 'Darwin':
    import multiprocessing
    multiprocessing.set_start_method('fork', force=True)

# PEP 517 builds do not have . in sys.path
sys.path.insert(0, os.path.dirname(__file__))

from sage_setup.excepthook import excepthook
sys.excepthook = excepthook

from sage_setup.setenv import setenv
setenv()

import sage.env
sage.env.default_required_modules = sage.env.default_optional_modules = ()

from sage_setup.command.sage_build_cython import sage_build_cython
from sage_setup.command.sage_build_ext import sage_build_ext

cmdclass = dict(build_cython=sage_build_cython,
                build_ext=sage_build_ext)

from sage_setup.find import find_python_sources
python_packages, python_modules, cython_modules = find_python_sources(
    '.', ['sage'])   # for now, we do the filtering using MANIFEST

log.warning('python_packages = {0}'.format(python_packages))
log.warning('python_modules = {0}'.format(python_modules))
log.warning('cython_modules = {0}'.format(cython_modules))

setup(
    cmdclass = cmdclass,
    packages = python_packages,
    py_modules  = python_modules,
    ext_modules = cython_modules,
)
