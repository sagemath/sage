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
from contextlib import contextmanager


class SourceDistributionFilter:
    r"""
    A :class:`collections.abc.Container` for source files in distributions.

    INPUT:

    - ``include_distributions`` -- (default: ``None``) if not ``None``,
      should be a sequence or set of strings: include files whose
      ``distribution`` (from a ``# sage_setup: distribution = PACKAGE``
      directive in the source file) is an element of ``distributions``.

    - ``exclude_distributions`` -- (default: ``None``) if not ``None``,
      should be a sequence or set of strings: exclude files whose
      ``distribution`` (from a ``# sage_setup: distribution = PACKAGE``
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


def read_distribution(src_file):
    """
    Parse ``src_file`` for a ``# sage_setup: distribution = PKG`` directive.

    INPUT:

    - ``src_file`` -- file name of a Python or Cython source file

    OUTPUT:

    - a string, the name of the distribution package (``PKG``); or the empty
      string if no directive was found.

    EXAMPLES::

        sage: from sage.env import SAGE_SRC
        sage: from sage.misc.package_dir import read_distribution
        sage: read_distribution(os.path.join(SAGE_SRC, 'sage', 'graphs', 'graph_decompositions', 'tdlib.pyx'))
        'sagemath-tdlib'
        sage: read_distribution(os.path.join(SAGE_SRC, 'sage', 'graphs', 'graph_decompositions', 'modular_decomposition.py'))
        ''
    """
    from Cython.Utils import open_source_file
    with open_source_file(src_file, error_handling='ignore') as fh:
        for line in fh:
            # Adapted from Cython's Build/Dependencies.py
            line = line.lstrip()
            if not line:
                continue
            if line[0] != '#':
                break
            line = line[1:].lstrip()
            kind = "sage_setup:"
            if line.startswith(kind):
                key, _, value = [s.strip() for s in line[len(kind):].partition('=')]
                if key == "distribution":
                    return value
    return ''


def is_package_or_sage_namespace_package_dir(path, *, distribution_filter=None):
    r"""
    Return whether ``path`` is a directory that contains a Python package.

    Ordinary Python packages are recognized by the presence of ``__init__.py``.

    Implicit namespace packages (PEP 420) are only recognized if they
    follow the conventions of the Sage library, i.e., the directory contains
    a file ``all.py`` or a file matching the pattern ``all__*.py``
    such as ``all__sagemath_categories.py``.

    INPUT:

    - ``path`` -- a directory name.

    - ``distribution_filter`` -- (optional, default: ``None``)
      only consider ``all*.py`` files whose distribution (from a
      ``# sage_setup: distribution = PACKAGE`` directive in the source file)
      is an element of ``distribution_filter``.

    EXAMPLES:

    :mod:`sage.cpython` is an ordinary package::

        sage: from sage.misc.package_dir import is_package_or_sage_namespace_package_dir
        sage: directory = os.path.dirname(sage.cpython.__file__); directory
        '.../sage/cpython'
        sage: is_package_or_sage_namespace_package_dir(directory)
        True

    :mod:`sage.libs.mpfr` only has an ``__init__.pxd`` file, but we consider
    it a package directory for consistency with Cython::

        sage: directory = os.path.join(os.path.dirname(sage.libs.all.__file__), 'mpfr'); directory
        '.../sage/libs/mpfr'
        sage: is_package_or_sage_namespace_package_dir(directory)
        True

    :mod:`sage` is designated to become an implicit namespace package::

        sage: directory = os.path.dirname(sage.env.__file__); directory
        '.../sage'
        sage: is_package_or_sage_namespace_package_dir(directory)
        True

    Not a package::

        sage: directory = os.path.join(os.path.dirname(sage.symbolic.__file__), 'ginac'); directory
        '.../sage/symbolic/ginac'
        sage: is_package_or_sage_namespace_package_dir(directory)
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
    Activate namespace package support in Cython 0.x

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
