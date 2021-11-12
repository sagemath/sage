"""
Optional extensions

An "optional extension" is a Cython extension which is always
cythonized (i.e. converted to a .c or .cpp file), but which is only
compiled depending on some condition. Typically, this condition is a
package which must be installed.
"""

# ****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from setuptools.extension import Extension
from sage.misc.package import list_packages

all_packages = list_packages(local=True)


class CythonizeExtension(Extension):
    """
    A class for extensions which are only cythonized, but not built.

    The file ``src/setup.py`` contains some logic to check the
    ``skip_build`` attribute of extensions.

    EXAMPLES::

        sage: from sage_setup.optional_extension import CythonizeExtension
        sage: ext = CythonizeExtension("foo", ["foo.c"])
        sage: ext.skip_build
        True
    """
    skip_build = True


def is_package_installed_and_updated(pkg):
    from sage.misc.package import is_package_installed
    try:
        pkginfo = all_packages[pkg]
    except KeyError:
        # Might be an installed old-style package
        condition = is_package_installed(pkg)
    else:
        condition = (pkginfo.installed_version == pkginfo.remote_version)
    return condition


def OptionalExtension(*args, **kwds):
    """
    If some condition (see INPUT) is satisfied, return an ``Extension``.
    Otherwise, return a ``CythonizeExtension``.

    Typically, the condition is some optional package or something
    depending on the operating system.

    INPUT:

    - ``condition`` -- (boolean) the actual condition

    - ``package`` -- (string) the condition is that this package is
      installed and up-to-date (only used if ``condition`` is not given)

    EXAMPLES::

        sage: from sage_setup.optional_extension import OptionalExtension
        sage: ext = OptionalExtension("foo", ["foo.c"], condition=False)
        sage: print(ext.__class__.__name__)
        CythonizeExtension
        sage: ext = OptionalExtension("foo", ["foo.c"], condition=True)
        sage: print(ext.__class__.__name__)
        Extension
        sage: ext = OptionalExtension("foo", ["foo.c"], package="no_such_package")
        sage: print(ext.__class__.__name__)
        CythonizeExtension
        sage: ext = OptionalExtension("foo", ["foo.c"], package="gap")
        sage: print(ext.__class__.__name__)
        Extension
    """
    try:
        condition = kwds.pop("condition")
    except KeyError:
        pkg = kwds.pop("package")
        condition = is_package_installed_and_updated(pkg)
    if condition:
        return Extension(*args, **kwds)
    else:
        return CythonizeExtension(*args, **kwds)
