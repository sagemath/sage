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
import sys
from contextlib import contextmanager


def is_package_or_sage_namespace_package_dir(path):
    r"""
    Return whether ``path`` is a directory that contains a Python package.

    Ordinary Python packages are recognized by the presence of ``__init__.py``.

    Implicit namespace packages (PEP 420) are only recognized if they
    follow the conventions of the Sage library, i.e., the directory contains
    a file ``all.py``.

    INPUT:

    - ``path`` -- a directory name

    EXAMPLES:

    :mod:`sage.cpython` is an ordinary package::

        sage: # optional - !meson_editable
        sage: from sage.misc.package_dir import is_package_or_sage_namespace_package_dir
        sage: directory = sage.cpython.__path__[0]; directory
        '.../sage/cpython'
        sage: is_package_or_sage_namespace_package_dir(directory)
        True

    :mod:`sage.libs.mpfr` only has an ``__init__.pxd`` file, but we consider
    it a package directory for consistency with Cython::

        sage: # optional - !meson_editable
        sage: directory = os.path.join(sage.libs.__path__[0], 'mpfr'); directory
        '.../sage/libs/mpfr'
        sage: is_package_or_sage_namespace_package_dir(directory)       # known bug (seen in build.yml)
        True

    :mod:`sage` is designated to become an implicit namespace package::

        sage: # optional - !meson_editable
        sage: directory = sage.__path__[0]; directory
        '.../sage'
        sage: is_package_or_sage_namespace_package_dir(directory)
        True

    Not a package::

        sage: # optional - !meson_editable
        sage: directory = os.path.join(sage.symbolic.__path__[0], 'ginac'); directory   # needs sage.symbolic
        '.../sage/symbolic/ginac'
        sage: is_package_or_sage_namespace_package_dir(directory)                       # needs sage.symbolic
        False

    TESTS::

        sage: # optional - meson_editable
        sage: from sage.misc.package_dir import is_package_or_sage_namespace_package_dir
        sage: directory = os.path.dirname(sage.cpython.__file__); directory
        '.../sage/cpython'
        sage: is_package_or_sage_namespace_package_dir(directory)
        True

        sage: # optional - meson_editable
        sage: directory = os.path.join(os.path.dirname(sage.libs.__file__), 'mpfr'); directory
        '.../sage/libs/mpfr'
        sage: is_package_or_sage_namespace_package_dir(directory)
        True

        sage: # optional - meson_editable, sage.symbolic
        sage: directory = os.path.join(os.path.dirname(sage.symbolic.__file__), 'ginac'); directory
        '.../sage/symbolic/ginac'
        sage: is_package_or_sage_namespace_package_dir(directory)
        False
    """
    if os.path.exists(os.path.join(path, "__init__.py")):  # ordinary package
        return True
    if os.path.exists(
        os.path.join(path, "__init__.pxd")
    ):  # for consistency with Cython
        return True
    fname = os.path.join(path, "all.py")
    if os.path.exists(fname):
        return True  # namespace package
    return False


@contextmanager
def cython_namespace_package_support():
    r"""
    Activate namespace package support in Cython 0.x.

    See https://github.com/cython/cython/issues/2918#issuecomment-991799049
    """
    import Cython.Build.Cythonize
    import Cython.Build.Dependencies
    import Cython.Utils

    orig_is_package_dir = Cython.Utils.is_package_dir
    Cython.Utils.is_package_dir = Cython.Build.Cythonize.is_package_dir = (
        Cython.Build.Dependencies.is_package_dir
    ) = Cython.Utils.cached_function(is_package_or_sage_namespace_package_dir)
    try:
        yield
    finally:
        Cython.Utils.is_package_dir = Cython.Build.Cythonize.is_package_dir = (
            Cython.Build.Dependencies.is_package_dir
        ) = orig_is_package_dir


def walk_packages(path=None, prefix="", onerror=None):
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

        sage: # optional - !meson_editable
        sage: sorted(sage.misc.package_dir.walk_packages(sage.misc.__path__))  # a namespace package
        [..., ModuleInfo(module_finder=FileFinder('.../sage/misc'), name='package_dir', ispkg=False), ...]

    TESTS::

        sage: # optional - meson_editable
        sage: sorted(sage.misc.package_dir.walk_packages(sage.misc.__path__))
        [..., ModuleInfo(module_finder=<...MesonpyPathFinder object...>, name='package_dir', ispkg=False), ...]
    """
    # Adapted from https://github.com/python/cpython/blob/3.11/Lib/pkgutil.py

    def iter_modules(path=None, prefix=""):
        """
        Yield :class:`ModuleInfo` for all submodules on ``path``.
        """
        from pkgutil import ModuleInfo, get_importer, iter_importers

        if path is None:
            importers = iter_importers()
        elif isinstance(path, str):
            raise ValueError(
                "path must be None or list of paths to look for modules in"
            )
        else:
            importers = map(get_importer, path)

        yielded = {}
        for i in importers:
            for name, ispkg in iter_importer_modules(i, prefix):
                if name not in yielded:
                    yielded[name] = 1
                    yield ModuleInfo(i, name, ispkg)

    def _iter_importer_modules_helper(importer, prefix=""):
        r"""
        Helper function for :func:`iter_importer_modules`.
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
                path = os.path.join(importer.path, fn)
                ispkg = False

                if not modname and os.path.isdir(path) and "." not in fn:
                    modname = fn
                    if not (ispkg := is_package_or_sage_namespace_package_dir(path)):
                        continue

                if modname and "." not in modname:
                    yielded[modname] = 1
                    yield prefix + modname, ispkg

        elif not hasattr(importer, "iter_modules"):
            yield from []

        else:
            yield from importer.iter_modules(prefix)

    def iter_importer_modules(importer, prefix=""):
        r"""
        Yield :class:`ModuleInfo` for all modules of ``importer``.
        """
        for name, ispkg in sorted(_iter_importer_modules_helper(importer, prefix)):
            # we sort again for consistency of output ordering if importer is not
            # a FileFinder (needed in doctest of :func:`sage.misc.dev_tools/load_submodules`)
            modname = name.rsplit(".", 1)[-1]
            if modname in ["__init__", "all"] or modname.startswith("all__"):
                continue
            yield name, ispkg

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
                path = getattr(sys.modules[info.name], "__path__", None) or []

                # don't traverse path items we've seen before
                path = [p for p in path if not seen(p)]

                yield from walk_packages(path, info.name + ".", onerror)
