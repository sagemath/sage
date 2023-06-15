#!/usr/bin/env python

from distutils import log
from setuptools import setup

# Work around a Cython problem in Python 3.8.x on macOS
# https://github.com/cython/cython/issues/3262
import os
if os.uname().sysname == 'Darwin':
    import multiprocessing
    multiprocessing.set_start_method('fork', force=True)

# PEP 517 builds do not have . in sys.path
import sys
sys.path.insert(0, os.path.dirname(__file__))

from sage_setup.excepthook import excepthook
sys.excepthook = excepthook

from sage_setup.setenv import setenv
setenv()

import sage.env
sage.env.default_required_modules = sage.env.default_optional_modules = ('gsl', )
print(f'##################### {sage.env.SAGE_SRC=}')

if any(x in sys.argv
       for x in ['build', 'bdist_wheel', 'install']):
    log.info("Generating auto-generated sources")
    # from sage_setup.autogen import autogen_all
    # autogen_all()
    from sage_setup.autogen.interpreters import rebuild
    rebuild(os.path.join("sage", "ext", "interpreters"),
            interpreters=['CDF', 'RDF', 'RR', 'CC'],
            distribution='sagemath-modules')

from sage_setup.command.sage_build_cython import sage_build_cython
from sage_setup.command.sage_build_ext import sage_build_ext

cmdclass = dict(build_cython=sage_build_cython,
                build_ext=sage_build_ext)

from sage_setup.find import find_python_sources
python_packages, python_modules, cython_modules = find_python_sources(
    '.', ['sage'])   # for now, we do the filtering using MANIFEST

log.debug('python_packages = {0}'.format(python_packages))
log.debug('python_modules = {0}'.format(python_modules))
log.debug('cython_modules = {0}'.format(cython_modules))

setup(
    cmdclass = cmdclass,
    packages = python_packages,
    py_modules  = python_modules,
    ext_modules = cython_modules,
)
