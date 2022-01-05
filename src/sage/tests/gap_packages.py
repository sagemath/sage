"""
Test the optional GAP packages

TESTS::

    sage: from sage.tests.gap_packages import all_installed_packages, test_packages
    sage: pkgs = all_installed_packages(ignore_dot_gap=True)
    sage: test_packages(pkgs, only_failures=True)    # optional - gap_packages
      Status   Package   GAP Output
    +--------+---------+------------+

    sage: test_packages(['atlasrep', 'tomlib'])
      Status   Package    GAP Output
    +--------+----------+------------+
               atlasrep   true
               tomlib     true
"""

import os

from sage.libs.gap.libgap import libgap


def test_packages(packages, only_failures=False):
    """
    Return list of all installed packages.

    INPUT:

    - ``packages`` -- a list/tuple/iterable of strings. The names of
      GAP packages to try to import.

    - ``only_failures`` -- boolean, default ``False``. Whether to only
      include failures in the table.

    OUTPUT:

    A table of the installed packages and whether they load
    successfully.

    EXAMPLES::

        sage: from sage.tests.gap_packages import all_installed_packages, test_packages
        sage: test_packages(['GAPDoc'])
          Status   Package   GAP Output
        +--------+---------+------------+
                   GAPDoc    true

    All packages, including user-installed ones::

        sage: pkgs = all_installed_packages()
        sage: test_packages(pkgs)    # random output
          Status    Package      GAP Output
        +---------+------------+------------+
                   Alnuth       true
                   GAPDoc       true
                   sonata       true
                   tomlib       true
                   toric        true
    """
    rows = [['Status', 'Package', 'GAP Output']]
    for pkgdir in packages:
        # to allow weird suffixes e.g. 'qpa-version'
        pkg = pkgdir.split('-')[0]
        orig_warning_level = libgap.InfoLevel(libgap.InfoWarning)
        # Silence warnings about missing optional packages that might occur
        # when loading packages; they're not important for the purposes of this
        # test code
        libgap.SetInfoLevel(libgap.InfoWarning, 0)
        try:
            output = libgap.LoadPackage(pkg)
        finally:
            # Restore the original warning level
            libgap.SetInfoLevel(libgap.InfoWarning, orig_warning_level)

        ok = bool(output)
        status = '' if ok else 'Failure'
        if ok and only_failures:
            continue
        rows.append([status, pkg, str(output)])
    from sage.misc.table import table
    return table(rows, header_row=True)


def all_installed_packages(ignore_dot_gap=False, gap=None):
    """
    Return list of all installed packages.

    INPUT:

    - ``ignore_dot_gap`` -- Boolean (default: ``False``). Whether to
      ignore the `.gap/` directory (usually in the user home
      directory) when searching for packages.

    - ``gap`` -- The GAP interface to use (default: ``libgap``); can
      be either ``libgap`` or a pexpect ``Gap`` instance.

    OUTPUT:

    Tuple of strings in alphabetic order.

    EXAMPLES::

        sage: from sage.tests.gap_packages import all_installed_packages
        sage: all_installed_packages()
        (...'GAPDoc'...)
        sage: all_installed_packages(ignore_dot_gap=True) == all_installed_packages(gap=gap, ignore_dot_gap=True)
        True
    """
    if gap is None:
        gap = libgap

    if gap == libgap:
        paths = [str(p) for p in gap.eval('GAPInfo.RootPaths')]
    else:
        paths = [str(p) for p in gap('GAPInfo.RootPaths')]

    packages = []
    for path in paths:
        if ignore_dot_gap and path.endswith('/.gap/'):
            continue
        pkg_dir = os.path.join(path, 'pkg')
        if not os.path.exists(pkg_dir):
            continue
        for subdir in os.listdir(pkg_dir):
            if not os.path.isdir(os.path.join(pkg_dir, subdir)):
                continue
            packages.append(subdir.rstrip('-.0123456789'))
    packages.sort()
    return tuple(packages)
