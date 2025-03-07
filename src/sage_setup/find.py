# sage.doctest: needs SAGE_SRC
"""
Recursive Directory Contents
"""
# ****************************************************************************
#       Copyright (C) 2014      Volker Braun <vbraun.name@gmail.com>
#                     2014      R. Andrew Ohana
#                     2015-2018 Jeroen Demeyer
#                     2017      Erik M. Bray
#                     2021      Tobias Diez
#                     2020-2022 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

import importlib.machinery
import importlib.util

import os

from collections import defaultdict

from sage.misc.package_dir import is_package_or_sage_namespace_package_dir as is_package_or_namespace_package_dir
from sage.misc.package_dir import read_distribution, SourceDistributionFilter

assert read_distribution  # unused in this file, re-export for compatibility


def find_python_sources(src_dir, modules=['sage'], distributions=None,
                        exclude_distributions=None):
    """
    Find all Python packages and Python/Cython modules in the sources.

    INPUT:

    - ``src_dir`` -- root directory for the sources

    - ``modules`` -- (default: ``['sage']``) sequence of strings:
      the top-level directories in ``src_dir`` to be considered

    - ``distributions`` -- (default: ``None``) if not ``None``,
      should be a sequence or set of strings: only find modules whose
      ``distribution`` (from a ``# sage_setup: distribution = PACKAGE``
      directive in the module source file) is an element of
      ``distributions``.

    - ``exclude_distributions`` -- (default: ``None``) if not ``None``,
      should be a sequence or set of strings: exclude modules whose
      ``distribution`` (from a ``# sage_setup: distribution = PACKAGE``
      directive in the module source file) is in ``exclude_distributions``.

    OUTPUT:

    Triple consisting of

    - the list of package names (corresponding to ordinary packages
      or namespace packages, according to
      :func:`sage.misc.namespace_package.is_package_or_sage_namespace_package_dir`)

    - Python module names (corresponding to ``*.py`` files in ordinary
      or namespace packages)

    - Cython extensions (corresponding to ``*.pyx`` files).

    Both use dot as separator.

    EXAMPLES::

        sage: from sage.env import SAGE_SRC
        sage: from sage_setup.find import find_python_sources
        sage: py_packages, py_modules, cy_modules = find_python_sources(SAGE_SRC)

    Ordinary package (with ``__init__.py``)::

        sage: ['sage.structure' in L for L in (py_packages, py_modules)]
        [True, False]

    Python module in an ordinary package::

        sage: ['sage.structure.formal_sum' in L for L in (py_packages, py_modules)]
        [False, True]

    Cython module in an ordinary package::

        sage: ['sage.structure.sage_object' in L for L in (py_packages, py_modules)]
        [False, False]

    Subdirectory without any Python files::

        sage: ['sage.doctest.tests' in L for L in (py_packages, py_modules)]
        [False, False]

    Another subdirectory that is neither an ordinary nor a namespace package::

        sage: ['sage.extdata' in L for L in (py_packages, py_modules)]
        [False, False]

    Package designated to become an implicit namespace package (no ``__init__.py``, PEP 420,
    but with an ``all.py`` file per Sage library conventions)::

        sage: ['sage.graphs.graph_decompositions' in L for L in (py_packages, py_modules)]
        [True, False]

    Python module in a package designated to become an implicit namespace package::

        sage: ['sage.graphs.graph_decompositions.modular_decomposition' in L for L in (py_packages, py_modules)]
        [False, True]

    Python file (not module) in a directory that is neither an ordinary nor a namespace package::

        sage: ['sage.ext_data.nbconvert.postprocess' in L for L in (py_packages, py_modules)]
        [False, False]

    Filtering by distribution (distribution package)::

        sage: find_python_sources(SAGE_SRC, distributions=['sagemath-tdlib'])
        ([], [], [<setuptools.extension.Extension('sage.graphs.graph_decompositions.tdlib')...>])

    Benchmarking::

        sage: timeit('find_python_sources(SAGE_SRC)',         # random output
        ....:        number=1, repeat=1)
        1 loops, best of 1: 30 ms per loop

        sage: timeit('find_python_sources(SAGE_SRC, distributions=[""])', # random output
        ....:        number=1, repeat=1)
        1 loops, best of 1: 850 ms per loop

        sage: find_python_sources(SAGE_SRC, modules=['sage_setup'])
        (['sage_setup', ...], [...'sage_setup.find'...], [])
    """
    from setuptools import Extension

    PYMOD_EXT = get_extensions('source')[0]

    python_packages = []
    python_modules = []
    cython_modules = []

    distribution_filter = SourceDistributionFilter(distributions, exclude_distributions)

    cwd = os.getcwd()
    try:
        os.chdir(src_dir)
        for module in modules:
            for dirpath, dirnames, filenames in os.walk(module):
                package = dirpath.replace(os.path.sep, '.')
                if not is_package_or_namespace_package_dir(dirpath):
                    # Skip any subdirectories
                    dirnames[:] = []
                    continue
                # Ordinary package or namespace package.
                if distributions is None or '' in distributions:
                    python_packages.append(package)

                for filename in filenames:
                    base, ext = os.path.splitext(filename)
                    filepath = os.path.join(dirpath, filename)
                    if ext == PYMOD_EXT and base != '__init__':
                        if filepath in distribution_filter:
                            python_modules.append(package + '.' + base)
                    if ext == '.pyx':
                        if filepath in distribution_filter:
                            cython_modules.append(Extension(package + '.' + base,
                                                            sources=[os.path.join(dirpath, filename)]))

    finally:
        os.chdir(cwd)
    return python_packages, python_modules, cython_modules


def filter_cython_sources(src_dir, distributions, exclude_distributions=None):
    """
    Find all Cython modules in the given source directory that belong to the
    given distributions.

    INPUT:

    - ``src_dir`` -- root directory for the sources

    - ``distributions`` -- a sequence or set of strings: only find modules whose
      ``distribution`` (from a ``# sage_setup: distribution = PACKAGE``
      directive in the module source file) is an element of
      ``distributions``.

    OUTPUT: list of absolute paths to Cython files (``*.pyx``)

    EXAMPLES::

        sage: from sage.env import SAGE_SRC
        sage: from sage_setup.find import filter_cython_sources
        sage: cython_modules = filter_cython_sources(SAGE_SRC, ["sagemath-tdlib"])

    Cython module relying on tdlib::

        sage: any(f.endswith('sage/graphs/graph_decompositions/tdlib.pyx') for f in cython_modules)
        True

    Cython module not relying on tdlib::

        sage: any(f.endswith('sage/structure/sage_object.pyx') for f in cython_modules)
        False

    Benchmarking::

        sage: timeit('filter_cython_sources(SAGE_SRC, ["sagemath-tdlib"])', # random output
        ....:        number=1, repeat=1)
        1 loops, best of 1: 850 ms per loop
    """
    files: list[str] = []
    distribution_filter = SourceDistributionFilter(distributions, exclude_distributions)
    for dirpath, dirnames, filenames in os.walk(src_dir):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            base, ext = os.path.splitext(filename)
            if ext == '.pyx' and filepath in distribution_filter:
                files.append(filepath)

    return files


def _cythonized_dir(src_dir=None, editable_install=None):
    """
    Return the path where Cython-generated files are placed by the build system.

    INPUT:

    - ``src_dir`` -- string or path (default: the value of ``SAGE_SRC``).  The
      root directory for the sources.

    - ``editable_install`` -- boolean (default: determined from the existing
      installation). Whether this is an editable install of the Sage library.

    EXAMPLES::

        sage: from sage_setup.find import _cythonized_dir
        sage: from sage.env import SAGE_SRC
        sage: _cythonized_dir(SAGE_SRC)
        PosixPath('...')
        sage: _cythonized_dir(SAGE_SRC, editable_install=False) # optional - sage_spkg
        PosixPath('.../build/cythonized')

    """
    from importlib import import_module
    from pathlib import Path
    from sage.env import SAGE_ROOT, SAGE_SRC
    if editable_install is None:
        if src_dir is None:
            src_dir = SAGE_SRC
        src_dir = Path(src_dir)
        # sage.cpython is an ordinary package, so it has __file__
        sage_cpython = import_module('sage.cpython')
        d = Path(sage_cpython.__file__).resolve().parent.parent.parent
        editable_install = d == src_dir.resolve()
    if editable_install:
        # Editable install: Cython generates files in the source tree
        return src_dir
    else:
        return Path(SAGE_ROOT) / "build" / "pkgs" / "sagelib" / "src" / "build" / "cythonized"


def find_extra_files(src_dir, modules, cythonized_dir, special_filenames=[], *,
                     distributions=None):
    """
    Find all extra files which should be installed.

    These are (for each ``module`` in ``modules``):

    1. From ``src_dir/module``: all .pyx, .pxd and .pxi files and files
       listed in ``special_filenames``.
    2. From ``cythonized_dir/module``: all .h files (both the .h files
       from the sources, as well as all Cython-generated .h files).

    The directories are searched recursively, but only package
    directories (containing ``__init__.py`` or a Cython equivalent)
    are considered.

    INPUT:

    - ``src_dir`` -- root directory for the sources

    - ``modules`` -- sequence of strings:
      the top-level directories in ``src_dir`` to be considered

    - ``cythonized_dir`` -- the directory where the Cython-generated
      files are

    - ``special_filenames`` -- a list of filenames to be installed from
      ``src_dir``

    - ``distributions`` -- (default: ``None``) if not ``None``,
      should be a sequence or set of strings: only find files whose
      ``distribution`` (from a ``# sage_setup: distribution = PACKAGE``
      directive in the file) is an element of ``distributions``.

    OUTPUT: dict with items ``{dir: files}`` where ``dir`` is a
    directory relative to ``src_dir`` and ``files`` is a list of
    filenames inside that directory.

    EXAMPLES::

        sage: from sage_setup.find import find_extra_files, _cythonized_dir
        sage: from sage.env import SAGE_SRC, SAGE_ROOT
        sage: cythonized_dir = _cythonized_dir(SAGE_SRC)
        sage: extras = find_extra_files(SAGE_SRC, ["sage"], cythonized_dir)
        sage: extras["sage/libs/mpfr"]
        [...sage/libs/mpfr/types.pxd...]
        sage: sorted(extras["sage/ext/interpreters"])
        ['.../sage/ext/interpreters/wrapper_cdf.h', ...wrapper_cdf.pxd...]
    """
    data_files = {}
    cy_exts = ('.pxd', '.pxi', '.pyx')

    cwd = os.getcwd()
    try:
        os.chdir(src_dir)
        for module in modules:
            for dir, dirnames, filenames in os.walk(module):
                if not is_package_or_namespace_package_dir(dir):
                    continue
                sdir = os.path.join(src_dir, dir)
                cydir = os.path.join(cythonized_dir, dir)

                files = [os.path.join(sdir, f) for f in filenames
                         if f.endswith(cy_exts) or f in special_filenames]
                if os.path.isdir(cydir):  # Not every directory contains Cython files
                    files += [os.path.join(cydir, f) for f in os.listdir(cydir)
                              if f.endswith(".h")]
                else:
                    files += [os.path.join(sdir, f) for f in filenames
                              if f.endswith(".h")]

                if distributions is not None:
                    files = [f for f in files
                             if read_distribution(f) in distributions]

                if files:
                    data_files[dir] = files
    finally:
        os.chdir(cwd)

    return data_files


def installed_files_by_module(site_packages, modules=('sage',)):
    """
    Find all currently installed files

    INPUT:

    - ``site_packages`` -- string. The root Python path where the Sage
      library is being installed. If the path doesn't exist, returns
      an empty dictionary.

    - ``modules`` -- list/tuple/iterable of strings (default:
      ``('sage',)``). The top-level directory name(s) in
      ``site_packages``.

    OUTPUT:

    A dictionary whose keys are module names (``'sage.module.foo'``)
    and values are list of corresponding file names
    ``['sage/module/foo.py', 'sage/module/foo.pyc']`` relative to
    ``site_packages``.

    EXAMPLES::

        sage: site_packages = os.path.dirname(os.path.dirname(os.path.dirname(sage.cpython.__file__)))
        sage: from sage_setup.find import installed_files_by_module
        sage: files_by_module = installed_files_by_module(site_packages)
        sage: (f,) = files_by_module['sage.structure.sage_object']; f
        'sage/structure/sage_object...'
        sage: (f1, f2) = sorted(files_by_module['sage.structure'])
        sage: f1
        'sage/structure/__init__.py'
        sage: f2
        'sage/structure/....pyc'

    This takes about 30ms with warm cache::

        sage: timeit('installed_files_by_module(site_packages)',       # random output
        ....:        number=1, repeat=1)
        1 loops, best of 1: 29.6 ms per loop
    """

    module_files = defaultdict(set)
    module_exts = get_extensions()

    def add(module, filename, dirpath):
        # Find the longest extension that matches the filename
        best_ext = ''

        for ext in module_exts:
            if filename.endswith(ext) and len(ext) > len(best_ext):
                best_ext = ext

        if not best_ext:
            return

        base = filename[:-len(best_ext)]
        filename = os.path.join(dirpath, filename)

        if base != '__init__':
            module += '.' + base

        module_files[module].add(filename)

        cache_filename = importlib.util.cache_from_source(filename)
        if os.path.exists(cache_filename):
            module_files[module].add(cache_filename)

    cwd = os.getcwd()
    try:
        os.chdir(site_packages)
    except OSError:
        return module_files
    try:
        for module in modules:
            for dirpath, dirnames, filenames in os.walk(module):
                module_dir = '.'.join(dirpath.split(os.path.sep))

                if os.path.basename(dirpath) == '__pycache__':
                    continue

                for filename in filenames:
                    add(module_dir, filename, dirpath)
    finally:
        os.chdir(cwd)
    return module_files


def get_extensions(type=None):
    """
    Returns the filename extensions for different types of Python module files.

    By default returns all extensions, but can be filtered by type.  The
    possible types are 'source' (for pure Python sources), 'bytecode' (for
    compiled bytecode files (i.e. pyc files), or 'extension' for C extension
    modules.

    INPUT:

    - ``type`` -- the module type ('source', 'bytecode', or 'extension') or
      None

    EXAMPLES::

        sage: from sage_setup.find import get_extensions
        sage: get_extensions()  # random - depends on Python version
        ['.so', 'module.so', '.py', '.pyc']
        sage: get_extensions('source')
        ['.py']
        sage: get_extensions('bytecode')
        ['.pyc']
        sage: get_extensions('extension')  # random - depends on Python version
        ['.so', 'module.so']
    """

    if type:
        type = type.lower()
        if type not in ('source', 'bytecode', 'extension'):
            raise ValueError(
                "type must by one of 'source' (for Python sources), "
                "'bytecode' (for compiled Python bytecoe), or 'extension' "
                "(for C extension modules).")

    # Note: There is at least one case, for extension modules, where the
    # 'extension' does not begin with '.', but rather with 'module', for cases
    # in Python's stdlib, for example, where an extension module can be named
    # like "<modname>module.so".  This breaks things for us if we have a Cython
    # module literally named "module".
    return [ext for ext in _get_extensions(type) if ext[0] == '.']


def _get_extensions(type):
    """
    Python 3.3+ implementation of ``get_extensions()`` using the
    `importlib.extensions` module.
    """

    if type:
        return {'source': importlib.machinery.SOURCE_SUFFIXES,
                'bytecode': importlib.machinery.BYTECODE_SUFFIXES,
                'extension': importlib.machinery.EXTENSION_SUFFIXES}[type]

    return importlib.machinery.all_suffixes()
