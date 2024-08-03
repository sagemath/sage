# sage_setup: distribution = sagemath-environment
"""
Recognizing package directories
"""
# ****************************************************************************
#       Copyright (C) 2020-2022 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import glob
import re
import sys

from collections import defaultdict
from contextlib import contextmanager


class SourceDistributionFilter:
    r"""
    A :class:`collections.abc.Container` for source files in distributions.

    INPUT:

    - ``include_distributions`` -- (default: ``None``) if not ``None``,
      should be a sequence or set of strings: include files whose
      ``distribution`` (from a ``# sage_setup:`` ``distribution = PACKAGE``
      directive in the source file) is an element of ``distributions``.

    - ``exclude_distributions`` -- (default: ``None``) if not ``None``,
      should be a sequence or set of strings: exclude files whose
      ``distribution`` (from a ``# sage_setup:`` ``distribution = PACKAGE``
      directive in the module source file) is in ``exclude_distributions``.

    EXAMPLES::

        sage: from sage.misc.package_dir import SourceDistributionFilter
        sage: F = SourceDistributionFilter()
        sage: sage.misc.package_dir.__file__ in F
        True
        sage: F = SourceDistributionFilter(include_distributions=['sagemath-environment'])
        sage: sage.misc.package_dir.__file__ in F
        True
        sage: F = SourceDistributionFilter(exclude_distributions=['sagemath-environment'])
        sage: sage.misc.package_dir.__file__ in F
        False
    """
    def __init__(self, include_distributions=None, exclude_distributions=None):
        r"""
        TESTS:

        ``exclude_distributions=None`` is normalized to the empty tuple::

            sage: from sage.misc.package_dir import SourceDistributionFilter
            sage: F = SourceDistributionFilter()
            sage: F._exclude_distributions
            ()
        """
        self._include_distributions = include_distributions
        if exclude_distributions is None:
            exclude_distributions = ()
        self._exclude_distributions = exclude_distributions

    def __contains__(self, filename):
        r"""
        TESTS:

        No file access is used when neither ``include_distributions`` nor
        ``exclude_distributions`` is given::

            sage: from sage.misc.package_dir import SourceDistributionFilter
            sage: F = SourceDistributionFilter()
            sage: '/doesnotexist' in F
            True

        ``exclude_distributions`` can also be an empty container::

            sage: F = SourceDistributionFilter(exclude_distributions=())
            sage: '/doesnotexist' in F
            True
        """
        if self._include_distributions is None and not self._exclude_distributions:
            return True
        distribution = read_distribution(filename)
        if self._include_distributions is not None:
            if distribution not in self._include_distributions:
                return False
            return distribution not in self._exclude_distributions


distribution_directive = re.compile(r"(\s*#?\s*)(sage_setup:\s*distribution\s*=\s*([-_A-Za-z0-9]*))")


def read_distribution(src_file):
    r"""
    Parse ``src_file`` for a ``# sage_setup:`` ``distribution = PKG`` directive.

    INPUT:

    - ``src_file`` -- file name of a Python or Cython source file

    OUTPUT:

    A string, the name of the distribution package (``PKG``), or the empty
    string if no directive was found.

    EXAMPLES::

        sage: # needs SAGE_SRC
        sage: from sage.env import SAGE_SRC
        sage: from sage.misc.package_dir import read_distribution
        sage: read_distribution(os.path.join(SAGE_SRC, 'sage', 'graphs', 'graph_decompositions', 'tdlib.pyx'))
        'sagemath-tdlib'
        sage: read_distribution(os.path.join(SAGE_SRC, 'sage', 'graphs', 'graph_decompositions', 'modular_decomposition.py'))
        ''
    """
    with open(src_file, encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            # Adapted from Cython's Build/Dependencies.py
            line = line.lstrip()
            if not line:
                continue
            if line.startswith('#') or line.startswith(';'):
                line = line[1:].lstrip()
            elif line.startswith('/*') or line.startswith('//') or line.startswith(';;'):
                line = line[2:].lstrip()
            else:
                break
            kind = "sage_setup:"
            if line.startswith(kind):
                key, _, value = (s.strip() for s in line[len(kind):].partition('='))
                if key == "distribution":
                    return value
    return ''


def update_distribution(src_file, distribution, *, verbose=False):
    r"""
    Add or update a ``# sage_setup:`` ``distribution = PKG`` directive in ``src_file``.

    For a Python or Cython file, if a ``distribution`` directive
    is not already present, it is added.

    For any other file, if a ``distribution`` directive is not already
    present, no action is taken.

    INPUT:

    - ``src_file`` -- file name of a source file

    EXAMPLES::

        sage: from sage.misc.package_dir import read_distribution, update_distribution
        sage: import tempfile
        sage: def test(filename, file_contents):
        ....:     with tempfile.TemporaryDirectory() as d:
        ....:         fname = os.path.join(d, filename)
        ....:         with open(fname, 'w') as f:
        ....:             f.write(file_contents)
        ....:         with open(fname, 'r') as f:
        ....:             print(f.read() + "====")
        ....:         update_distribution(fname, 'sagemath-categories')
        ....:         with open(fname, 'r') as f:
        ....:             print(f.read() + "====")
        ....:         update_distribution(fname, '')
        ....:         with open(fname, 'r') as f:
        ....:             print(f.read(), end="")
        sage: test('module.py', '# Python file\n')
        # Python file
        ====
        # sage_setup: distribution...= sagemath-categories
        # Python file
        ====
        # sage_setup: distribution...=
        # Python file
        sage: test('file.cpp', '// sage_setup: ' 'distribution=sagemath-modules\n'
        ....:                  '// C++ file with existing directive\n')
        // sage_setup: distribution...=sagemath-modules
        // C++ file with existing directive
        ====
        // sage_setup: distribution...= sagemath-categories
        // C++ file with existing directive
        ====
        // sage_setup: distribution...=
        // C++ file with existing directive
        sage: test('file.cpp', '// C++ file without existing directive\n')
        // C++ file without existing directive
        ====
        // C++ file without existing directive
        ====
        // C++ file without existing directive
    """
    if not distribution:
        distribution = ''
    directive = 'sage_setup: ' f'distribution = {distribution}'.rstrip()
    try:
        with open(src_file) as f:
            src_lines = f.read().splitlines()
    except UnicodeDecodeError:
        # Silently skip binary files
        return
    any_found = False
    any_change = False
    for i, line in enumerate(src_lines):
        if m := distribution_directive.search(line):
            old_distribution = m.group(3)
            if any_found:
                # Found a second distribution directive; remove it.
                if not (line := distribution_directive.sub(r'', line)):
                    line = None
            else:
                line = distribution_directive.sub(fr'\1{directive}', line)
            if line != src_lines[i]:
                src_lines[i] = line
                any_change = True
                if verbose:
                    print(f"{src_file}: changed 'sage_setup: " f"distribution' "
                          f"from {old_distribution!r} to {distribution!r}")
            any_found = True
    if not any_found:
        if any(src_file.endswith(ext)
               for ext in [".pxd", ".pxi", ".py", ".pyx", ".sage"]):
            src_lines.insert(0, f'# {directive}')
            any_change = True
            if verbose:
                print(f"{src_file}: added 'sage_setup: "
                      f"distribution = {distribution}' directive")
    if not any_change:
        return
    with open(src_file, 'w') as f:
        for line in src_lines:
            if line is not None:
                f.write(line + '\n')


def is_package_or_sage_namespace_package_dir(path, *, distribution_filter=None):
    r"""
    Return whether ``path`` is a directory that contains a Python package.

    Ordinary Python packages are recognized by the presence of ``__init__.py``.

    Implicit namespace packages (PEP 420) are only recognized if they
    follow the conventions of the Sage library, i.e., the directory contains
    a file ``all.py`` or a file matching the pattern ``all__*.py``
    such as ``all__sagemath_categories.py``.

    INPUT:

    - ``path`` -- a directory name

    - ``distribution_filter`` -- (default: ``None``)
      only consider ``all*.py`` files whose distribution (from a
      ``# sage_setup:`` ``distribution = PACKAGE`` directive in the source file)
      is an element of ``distribution_filter``.

    EXAMPLES:

    :mod:`sage.cpython` is an ordinary package::

        sage: from sage.misc.package_dir import is_package_or_sage_namespace_package_dir
        sage: directory = sage.cpython.__path__[0]; directory
        '.../sage/cpython'
        sage: is_package_or_sage_namespace_package_dir(directory)
        True

    :mod:`sage.libs.mpfr` only has an ``__init__.pxd`` file, but we consider
    it a package directory for consistency with Cython::

        sage: directory = os.path.join(sage.libs.__path__[0], 'mpfr'); directory
        '.../sage/libs/mpfr'
        sage: is_package_or_sage_namespace_package_dir(directory)       # known bug (seen in build.yml)
        True

    :mod:`sage` is designated to become an implicit namespace package::

        sage: directory = sage.__path__[0]; directory
        '.../sage'
        sage: is_package_or_sage_namespace_package_dir(directory)       # known bug (seen in build.yml)
        True

    Not a package::

        sage: directory = os.path.join(sage.symbolic.__path__[0], 'ginac'); directory   # needs sage.symbolic
        '.../sage/symbolic/ginac'
        sage: is_package_or_sage_namespace_package_dir(directory)                       # needs sage.symbolic
        False
    """
    if os.path.exists(os.path.join(path, '__init__.py')):                # ordinary package
        return True
    if os.path.exists(os.path.join(path, '__init__.pxd')):               # for consistency with Cython
        return True
    fname = os.path.join(path, 'all.py')
    if os.path.exists(fname):
        if distribution_filter is None or fname in distribution_filter:  # complete namespace package
            return True
    for fname in glob.iglob(os.path.join(path, 'all__*.py')):
        if distribution_filter is None or fname in distribution_filter:  # partial namespace package
            return True
    return False


@contextmanager
def cython_namespace_package_support():
    r"""
    Activate namespace package support in Cython 0.x.

    See https://github.com/cython/cython/issues/2918#issuecomment-991799049
    """
    import Cython.Build.Dependencies
    import Cython.Build.Cythonize
    import Cython.Utils
    orig_is_package_dir = Cython.Utils.is_package_dir
    Cython.Utils.is_package_dir = Cython.Build.Cythonize.is_package_dir = Cython.Build.Dependencies.is_package_dir = Cython.Utils.cached_function(is_package_or_sage_namespace_package_dir)
    try:
        yield
    finally:
        Cython.Utils.is_package_dir = Cython.Build.Cythonize.is_package_dir = Cython.Build.Dependencies.is_package_dir = orig_is_package_dir


def walk_packages(path=None, prefix='', onerror=None):
    r"""
    Yield :class:`pkgutil.ModuleInfo` for all modules recursively on ``path``.

    This version of the standard library function :func:`pkgutil.walk_packages`
    addresses https://github.com/python/cpython/issues/73444 by handling
    the implicit namespace packages in the package layout used by Sage;
    see :func:`is_package_or_sage_namespace_package_dir`.

    INPUT:

    - ``path`` -- list of paths to look for modules in or
      ``None`` (all accessible modules)

    - ``prefix`` -- string to output on the front of every module name
      on output

    - ``onerror`` -- a function which gets called with one argument (the
      name of the package which was being imported) if any exception
      occurs while trying to import a package.  If ``None``, ignore
      :exc:`ImportError` but propagate all other exceptions.

    EXAMPLES::

        sage: sorted(sage.misc.package_dir.walk_packages(sage.misc.__path__))  # a namespace package
        [..., ModuleInfo(module_finder=FileFinder('.../sage/misc'), name='package_dir', ispkg=False), ...]
    """
    # Adapted from https://github.com/python/cpython/blob/3.11/Lib/pkgutil.py

    def iter_modules(path=None, prefix=''):
        """
        Yield :class:`ModuleInfo` for all submodules on ``path``.
        """
        from pkgutil import get_importer, iter_importers, ModuleInfo

        if path is None:
            importers = iter_importers()
        elif isinstance(path, str):
            raise ValueError("path must be None or list of paths to look for modules in")
        else:
            importers = map(get_importer, path)

        yielded = {}
        for i in importers:
            for name, ispkg in iter_importer_modules(i, prefix):
                if name not in yielded:
                    yielded[name] = 1
                    yield ModuleInfo(i, name, ispkg)

    def iter_importer_modules(importer, prefix=''):
        r"""
        Yield :class:`ModuleInfo` for all modules of ``importer``.
        """
        from importlib.machinery import FileFinder

        if isinstance(importer, FileFinder):
            if importer.path is None or not os.path.isdir(importer.path):
                return

            yielded = {}
            import inspect
            try:
                filenames = os.listdir(importer.path)
            except OSError:
                # ignore unreadable directories like import does
                filenames = []
            filenames.sort()  # handle packages before same-named modules

            for fn in filenames:
                modname = inspect.getmodulename(fn)
                if modname and (modname in ['__init__', 'all']
                                or modname.startswith('all__')
                                or modname in yielded):
                    continue

                path = os.path.join(importer.path, fn)
                ispkg = False

                if not modname and os.path.isdir(path) and '.' not in fn:
                    modname = fn
                    if not (ispkg := is_package_or_sage_namespace_package_dir(path)):
                        continue

                if modname and '.' not in modname:
                    yielded[modname] = 1
                    yield prefix + modname, ispkg

        elif not hasattr(importer, 'iter_modules'):
            yield from []

        else:
            yield from importer.iter_modules(prefix)

    def seen(p, m={}):
        if p in m:
            return True
        m[p] = True

    for info in iter_modules(path, prefix):
        yield info

        if info.ispkg:
            try:
                __import__(info.name)
            except ImportError:
                if onerror is not None:
                    onerror(info.name)
            except Exception:
                if onerror is not None:
                    onerror(info.name)
                else:
                    raise
            else:
                path = getattr(sys.modules[info.name], '__path__', None) or []

                # don't traverse path items we've seen before
                path = [p for p in path if not seen(p)]

                yield from walk_packages(path, info.name + '.', onerror)


def _all_filename(distribution):
    if not distribution:
        return 'all.py'
    return f"all__{distribution.replace('-', '_')}.py"


def _distribution_from_all_filename(filename):
    if m := re.match('all(__(.*?))?[.]py', filename):
        if distribution_per_all_filename := m.group(2):
            return distribution_per_all_filename.replace('_', '-')
        return ''
    return False


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(prog="sage --fixdistributions",
                            description="Maintenance tool for distribution packages of the Sage library",
                            epilog="By default, '%(prog)s' shows the distribution of each file.")
    parser.add_argument('--add', metavar='DISTRIBUTION', type=str, default=None,
                        help=("add a 'sage_setup: DISTRIBUTION' directive to FILES; "
                              "do not change files that already have a nonempty directive"))
    parser.add_argument('--set', metavar='DISTRIBUTION', type=str, default=None,
                        help="add or update the 'sage_setup: DISTRIBUTION' directive in FILES")
    parser.add_argument('--from-egg-info', action='store_true', default=False,
                        help="take FILES from pkgs/DISTRIBUTION/DISTRIBUTION.egg-info/SOURCES.txt")
    parser.add_argument("filename", metavar='FILES', nargs='*', type=str,
                        help=("source files or directories (default: all files from SAGE_SRC, "
                              "unless --from-egg-info, --add, or --set are used)"))

    args = parser.parse_args()

    distribution = args.set or args.add or ''

    if distribution == 'all':
        distributions = ["sagemath-bliss",
                         "sagemath-brial",
                         "sagemath-coxeter3",
                         "sagemath-mcqd",
                         "sagemath-meataxe",
                         "sagemath-sirocco",
                         "sagemath-tdlib",
                         "sagemath-environment",
                         "sagemath-categories",
                         "sagemath-repl",
                         "sagemath-objects"]
    else:
        distributions = [distribution.replace('_', '-')]

    if args.from_egg_info:
        if not distribution:
            print("Switch '--from-egg-info' must be used with either "
                  "'--add DISTRIBUTION' or '--set DISTRIBUTION'")
            sys.exit(1)
    elif not args.filename:
        if distribution:
            print("Switches '--add' and '--set' require the switch '--from-egg-info' "
                  "or one or more file or directory names")
            sys.exit(1)
        from sage.env import SAGE_SRC
        if (not SAGE_SRC
                or not os.path.exists(os.path.join(SAGE_SRC, 'sage'))
                or not os.path.exists(os.path.join(SAGE_SRC, 'conftest_test.py'))):
            print(f'{SAGE_SRC=} does not seem to contain a copy of the Sage source tree')
            sys.exit(1)
        args.filename = [os.path.join(SAGE_SRC, 'sage')]

    ordinary_packages = set()
    package_distributions_per_directives = defaultdict(set)     # path -> set of strings (distributions)
    package_distributions_per_all_files = defaultdict(set)      # path -> set of strings (distributions)

    def handle_file(root, file):
        path = os.path.join(root, file)
        if args.set is not None:
            update_distribution(path, distribution, verbose=True)
            file_distribution = distribution
        elif args.add is not None:
            if not (file_distribution := read_distribution(path)):
                update_distribution(path, distribution, verbose=True)
                file_distribution = distribution
        else:
            file_distribution = read_distribution(path)
            print(f'{path}: file in distribution {file_distribution!r}')
        package_distributions_per_directives[root].add(file_distribution)
        if file.startswith('__init__.'):
            ordinary_packages.add(root)
        elif (distribution_per_all_filename := _distribution_from_all_filename(file)) is False:
            # Not an all*.py file.
            pass
        elif not distribution_per_all_filename:
            # An all.py file.
            if file_distribution:
                # The all.py is declared to belong to a named distribution, that's OK
                package_distributions_per_all_files[root].add(file_distribution)
            else:
                pass
        else:
            # An all__*.py file
            if distribution_per_all_filename != file_distribution:
                print(f'{path}: file should go in distribution {distribution_per_all_filename!r}, not {file_distribution!r}')
            package_distributions_per_all_files[root].add(distribution_per_all_filename)

    for distribution in distributions:

        paths = list(args.filename)

        if args.from_egg_info:
            from sage.env import SAGE_ROOT
            if not distribution:
                print("Switch '--from-egg-info' must be used with either "
                      "'--add DISTRIBUTION' or '--set DISTRIBUTION'")
                sys.exit(1)
            if not SAGE_ROOT:
                print(f'{SAGE_ROOT=} does not seem to contain a copy of the Sage source root')
                sys.exit(1)
            distribution_dir = os.path.join(SAGE_ROOT, 'pkgs', distribution)
            if not os.path.exists(distribution_dir):
                print(f'{distribution_dir} does not exist')
                sys.exit(1)
            distribution_underscore = distribution.replace('-', '_')
            try:
                with open(os.path.join(distribution_dir,
                                       f'{distribution_underscore}.egg-info', 'SOURCES.txt')) as f:
                    paths.extend(os.path.join(SAGE_ROOT, 'src', line.strip())
                                 for line in f
                                 if line.startswith('sage/'))
                print(f"sage --fixdistributions: found egg-info of distribution {distribution!r}")
            except FileNotFoundError:
                if len(distributions) > 1:
                    print(f"sage --fixdistributions: distribution {distribution!r} does not have egg-info, skipping it; "
                          f"run 'make {distribution_underscore}-sdist' or 'make {distribution_underscore}' to create it")
                    continue
                else:
                    print(f"sage --fixdistributions: distribution {distribution!r} does not have egg-info; "
                          f"run 'make {distribution_underscore}-sdist' or 'make {distribution_underscore}' to create it")
                    sys.exit(1)

        for path in paths:
            path = os.path.relpath(path)
            if os.path.isdir(path):
                if not is_package_or_sage_namespace_package_dir(path):
                    print(f'{path}: non-package directory')
                else:
                    for root, dirs, files in os.walk(path):
                        for dir in sorted(dirs):
                            path = os.path.join(root, dir)
                            if any(dir.startswith(prefix) for prefix in ['.', 'build', 'dist', '__pycache__', '_vendor', '.tox']):
                                # Silently skip
                                dirs.remove(dir)
                            elif not is_package_or_sage_namespace_package_dir(path):
                                print(f'{path}: non-package directory')
                                dirs.remove(dir)
                        for file in sorted(files):
                            if any(file.endswith(ext) for ext in [".pyc", ".pyo", ".bak", ".so", "~"]):
                                continue
                            handle_file(root, file)
            else:
                handle_file(*os.path.split(path))

    print(f"sage --fixdistributions: checking consistency")

    for package in ordinary_packages:
        if len(package_distributions_per_directives[package]) > 1:
            print(f'{package}: ordinary packages (with __init__.py) cannot be split in several distributions ('
                  + ', '.join(f'{dist!r}'
                              for dist in sorted(package_distributions_per_directives[package])) + ')')

    for package, distributions_per_directives in package_distributions_per_directives.items():
        if package in ordinary_packages:
            pass
        elif ((missing_all_files := distributions_per_directives - package_distributions_per_all_files[package])
                and not (missing_all_files == {''} and len(distributions_per_directives) < 2)):
            s = '' if len(missing_all_files) == 1 else 's'
            print(f'{package}: missing file{s} ' + ', '.join(_all_filename(distribution)
                                                             for distribution in missing_all_files))
