#!/usr/bin/env python

import time

from distutils import log
from setuptools import setup
from sage_setup.command.sage_build_ext_minimal import sage_build_ext_minimal
from sage_setup.cython_options import compiler_directives, compile_time_env_variables
from sage_setup.extensions import create_extension
from sage_setup.find import find_python_sources, find_extra_files

# Work around a Cython problem in Python 3.8.x on macOS
# https://github.com/cython/cython/issues/3262
import os
if os.uname().sysname == 'Darwin':
    import multiprocessing
    multiprocessing.set_start_method('fork', force=True)

import sys

from sage_setup.excepthook import excepthook
sys.excepthook = excepthook

from sage_setup.setenv import setenv
setenv()

import sage.env
sage.env.default_required_modules = sage.env.default_optional_modules = ()

cmdclass = dict(build_ext=sage_build_ext_minimal)

sdist = len(sys.argv) > 1 and (sys.argv[1] in ["sdist", "egg_info", "dist_info"])

# ########################################################
# ## Discovering Sources
# ########################################################
if sdist:
    extensions = []
    python_modules = []
    python_packages = []
    package_data = {}
else:
    log.info("Discovering Python/Cython source code....")
    t = time.time()

    python_packages, python_modules, cython_modules = find_python_sources(
        '.', ['sage'], distributions=['sagemath-schemes'])
    extra_files = find_extra_files(
        '.', ['sage'], '/doesnotexist', distributions=['sagemath-schemes'])
    package_data = {"": [f
                         for pkg, files in extra_files.items()
                         for f in files ]}
    package_data.update({})
    python_packages += list(package_data)

    log.debug('python_packages = {0}'.format(sorted(python_packages)))
    log.debug('python_modules = {0}'.format(sorted(m if isinstance(m, str) else m.name for m in python_modules)))
    log.debug('cython_modules = {0}'.format(sorted(m if isinstance(m, str) else m.name for m in cython_modules)))
    log.debug('package_data = {0}'.format(package_data))

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
                cython_modules,
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

setup(cmdclass=cmdclass,
      packages=python_packages,
      py_modules=python_modules,
      ext_modules=extensions,
      package_data=package_data,
)
