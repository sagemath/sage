"""
Recognizing package directories
"""
import os
import glob

def is_package_or_sage_namespace_package_dir(path):
    r"""
    Return whether ``path`` is a directory that contains a Python package.

    Ordinary Python packages are recognized by the presence of `__init__.py`.

    Implicit namespace packages (PEP 420) are only recognized if they
    follow the conventions of the Sage library, i.e., the directory contains
    a file ``all.py`` or a file matching the pattern ``all__*.py``
    such as ``all__sagemath_categories.py``.

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
    if os.path.exists(os.path.join(path, '__init__.py')):   # ordinary package
        return True
    if os.path.exists(os.path.join(path, '__init__.pxd')):  # for consistency with Cython
        return True
    if os.path.exists(os.path.join(path, 'all.py')):        # complete namespace package
        return True
    for _ in glob.iglob(os.path.join(path, 'all__*.py')):
        return True                                         # partial namespace package
    return False
