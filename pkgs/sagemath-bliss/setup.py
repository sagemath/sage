#!/usr/bin/env python

from distutils import log
from setuptools import setup

# Work around a Cython problem in Python 3.8.x on macOS
# https://github.com/cython/cython/issues/3262
import os
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

# PEP 517 builds do not have . in sys.path
import sys
sys.path.insert(0, os.path.dirname(__file__))

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

    import sage.env
    sage.env.default_required_modules = sage.env.default_optional_modules = ()

    from sage_setup.command.sage_build_cython import sage_build_cython
    from sage_setup.command.sage_build_ext import sage_build_ext
    sage_build_cython.built_distributions = ['sagemath-bliss']

    cmdclass = dict(build_cython=sage_build_cython,
                    build_ext=sage_build_ext)

from sage_setup.find import find_python_sources
python_packages, python_modules, cython_modules = find_python_sources(
    '.', ['sage'], distributions=['sagemath-bliss'])

log.warn('python_packages = {0}'.format(python_packages))
log.warn('python_modules = {0}'.format(python_modules))
log.warn('cython_modules = {0}'.format(cython_modules))

setup(
    cmdclass = cmdclass,
    packages = python_packages,
    py_modules  = python_modules,
    ext_modules = cython_modules,
)
