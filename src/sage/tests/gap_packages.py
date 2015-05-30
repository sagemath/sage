"""
Test the optional GAP packages

TESTS::

    sage: from sage.tests.gap_packages import all_installed_packages, test_packages
    sage: pkgs = all_installed_packages(ignore_dot_gap=True)
    sage: test_packages(pkgs, only_failures=True)    # optional - gap_packages
      Status    Package    GAP Output
    +---------+----------+------------+
      Failure   HAPcryst   fail

These are packages in the ``database_gap`` package::

    sage: test_packages(['atlasrep', 'tomlib'])    # optional - database_gap
      Status   Package    GAP Output
    +--------+----------+------------+
               atlasrep   true
               tomlib     true
"""

import os
import os.path

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
          Failure   HAPcryst     fail
                    Hap          true
                    autpgrp      true
                    crime        true
                    ctbllib      true
                    design       true
                    factint      true
                    grape        true
                    guava        true
                    laguna       true
                    polycyclic   true
                    polymaking   true
                    sonata       true
                    toric        true
    """
    rows = [['Status', 'Package', 'GAP Output']]
    for pkg in packages:
        output = libgap.eval('LoadPackage("{0}")'.format(pkg))
        ok = bool(output)
        status = '' if ok else 'Failure'
        if ok and only_failures:
            continue
        rows.append([status, pkg, str(output)])
    from sage.misc.table import table
    return table(rows, header_row=True)



def all_installed_packages(ignore_dot_gap=False):
    """
    Return list of all installed packages.

    INPUT:

    - ``ignore_dot_gap`` -- Boolean (default: ``False``). Whether to
      ignore the `.gap/` directory (usually in the user home
      directory) when searching for packages.

    OUTPUT:

    Tuple of strings in alphabetic order.

    EXAMPLES::

        sage: from sage.tests.gap_packages import all_installed_packages
        sage: all_installed_packages()
        (...'GAPDoc',...)
    """
    packages = []
    for path in libgap.eval('GAP_ROOT_PATHS').sage():
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
