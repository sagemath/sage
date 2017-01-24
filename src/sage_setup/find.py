"""
Recursive Directory Contents
"""
#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import os


def find_python_sources(src_dir, modules=('sage',)):
    """
    Find all currently installed files

    INPUT:

    - ``src_dir`` -- string. The root directory for the sources.

    - ``module`` -- list/tuple/iterable of strings (default:
      ``('sage',)``). The top-level directory name(s) in ``src_dir``.

    OUTPUT:

    Pair consisting of the list of package names (directories with
    ``__init__.py``) and module names (remaining ``*.py`` files). Both
    use dot as separator.

    EXAMPLES::

        sage: from sage.env import SAGE_SRC
        sage: py_packages, py_modules = find_python_sources(SAGE_SRC)
        sage: examples = ['sage.structure', 'sage.structure.formal_sum',
        ....:             'sage.structure.sage_object', 'sage.doctest.tests']
        sage: [m in py_packages for m in examples]
        [True, False, False, False]
        sage: [m in py_modules for m in examples]
        [False, True, False, False]

        sage: timeit('find_python_sources(SAGE_SRC)',         # random output
        ....:        number=1, repeat=1)
        1 loops, best of 1: 18.8 ms per loop

        sage: find_python_sources(SAGE_SRC, modules=['sage_setup'])
        (['sage_setup', ...], [...'sage_setup.find'...])
    """
    PYMOD_EXT = os.path.extsep + 'py'
    INIT_FILE = '__init__' + PYMOD_EXT

    python_packages = []
    python_modules = []

    cwd = os.getcwd()
    try:
        os.chdir(src_dir)
        for module in modules:
            for dirpath, dirnames, filenames in os.walk(module):
                if INIT_FILE not in filenames:
                    continue
                dirpath = dirpath.replace(os.path.sep, '.')
                python_packages.append(dirpath)
                for filename in filenames:
                    base, ext = os.path.splitext(filename)
                    if ext == PYMOD_EXT and base != '__init__':
                        python_modules.append(dirpath + '.' + base)
    finally:
        os.chdir(cwd)
    return python_packages, python_modules


def find_extra_files(packages, src_dir, cythonized_dir, special_filenames=[]):
    """
    Find all extra files which should be installed.

    These are:

    1. From ``src_dir``: all .pxd and .pxi files and files listed in
       ``special_filenames``.
    2. From ``cythonized_dir``: all .h files (these are both the .h files
       from the Sage sources, as well as all Cython-generated .h files).

    INPUT:

    - ``packages`` -- a list of Python packages to be considered

    - ``src_dir`` -- the directory where to look for source files

    - ``cythonized_dir`` -- the directory where the Cython-generated
      files are

    - ``special_filenames`` -- a list of filenames to be installed from
      ``src_dir``

    EXAMPLES::

        sage: from sage_setup.find import find_extra_files
        sage: from sage.env import SAGE_SRC, SAGE_CYTHONIZED
        sage: find_extra_files(["sage.modular.arithgroup"], SAGE_SRC, SAGE_CYTHONIZED)
        [('sage/modular/arithgroup',
          ['.../src/sage/modular/arithgroup/farey.pxd', ...farey_symbol.h...])]
    """
    data_files = []

    for package in packages:
        dir = package.replace('.', os.path.sep)
        sdir = os.path.join(src_dir, dir)
        cydir = os.path.join(cythonized_dir, dir)

        files = [os.path.join(sdir, f) for f in os.listdir(sdir)
                if f.endswith((".pxd", ".pxi")) or f in special_filenames]
        if os.path.isdir(cydir):  # Not every directory contains Cython files
            files += [os.path.join(cydir, f) for f in os.listdir(cydir)
                    if f.endswith(".h")]

        if files:
            data_files.append((dir, files))

    return data_files


def installed_files_by_module(site_packages, modules=('sage',)):
    """
    Find all currently installed files

    INPUT:

    - ``site_packages`` -- string. The root Python path where the Sage
      library is being installed. If the path doesn't exist, returns
      an empty dictionary.

    - ``module`` -- list/tuple/iterable of strings (default:
      ``('sage',)``). The top-level directory name(s) in
      ``site_packages``.

    OUTPUT:

    A dictionary whose keys are module names (``'sage.module.foo'``)
    and values are list of corresponding file names
    ``['sage/module/foo.py', 'sage/module/foo.pyc']`` relative to
    ``site_packages``.

    EXAMPLES::

        sage: from site import getsitepackages
        sage: site_packages = getsitepackages()[0]
        sage: files_by_module = installed_files_by_module(site_packages)
        sage: from sage.misc.sageinspect import loadable_module_extension
        sage: 'sage/structure/sage_object' + loadable_module_extension() in \
        ....:     files_by_module['sage.structure.sage_object']
        True
        sage: sorted(files_by_module['sage.structure'])
        ['sage/structure/__init__.py', 'sage/structure/__init__.pyc']

    This takes about 30ms with warm cache:

        sage: timeit('installed_files_by_module(site_packages)',       # random output
        ....:        number=1, repeat=1)
        1 loops, best of 1: 29.6 ms per loop
    """
    module_files = dict()
    def add(module, filename):
        module_files.setdefault(module, set([filename])).add(filename)

    cwd = os.getcwd()
    try:
        os.chdir(site_packages)
    except OSError:
        return module_files
    try:
        for module in modules:
            for dirpath, dirnames, filenames in os.walk(module):
                module_dir = '.'.join(dirpath.split(os.path.sep))
                for filename in filenames:
                    base, ext = os.path.splitext(filename)
                    filename = os.path.join(dirpath, filename)
                    if base == '__init__':
                        add(module_dir, filename)
                    else:
                        add(module_dir + '.' + base, filename)
    finally:
        os.chdir(cwd)
    return module_files
