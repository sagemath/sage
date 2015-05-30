r"""
Sage package management commands

A Sage package has the extension .spkg. It is a tarball that is
(usually) bzip2 compressed that contains arbitrary data and an
spkg-install file. An Sage package typically has the following
components:


-  spkg-install - shell script that is run to install the package

-  Sage.txt - file that describes how the package was made, who
   maintains it, etc.

-  sage - directory with extra patched version of files that needed
   during the install


Use the ``install_package`` command to install a new
package, and use ``optional_packages`` to list all
optional packages available on the central Sage server. The
``upgrade`` command upgrades all *standard* packages -
there is no auto-upgrade command for optional packages.

All package management can also be done via the Sage command line.
"""

import os

__installed_packages = None


def install_all_optional_packages(force=True, dry_run=False):
    r"""
    Install all available optional spkg's in the official Sage spkg
    repository.  Returns a list of all spkg's that *fail* to install.

    INPUT:

    - ``force`` -- bool (default: ``True``); whether to force
      reinstall of spkg's that are already installed.

    - ``dry_run`` -- bool (default: ``False``); if ``True``, just list
      the packages that would be installed in order, but don't
      actually install them.

    OUTPUT:

        list of strings

    .. NOTE::

        This is designed mainly for testing purposes.  This also
        doesn't do anything with respect to dependencies -- the
        packages are installed in alphabetical order.  Dependency
        issues will be dealt with in a future version.


    AUTHOR:

        -- William Stein (2008-12)

    EXAMPLES::

        sage: sage.misc.package.install_all_optional_packages(dry_run=True)  # optional - internet
        Installing ...
        []
    """
    # Get list of all packages from the package server
    installed, not_installed = optional_packages()
    if force:
        v = installed + not_installed
    else:
        v = not_installed
    failed = []

    # sort the list of packages in alphabetical order
    v.sort()

    # install each one
    for pkg in v:
        try:
            print("Installing {}...".format(pkg))
            if not dry_run:
                # only actually do the install of the package if dry_run is False
                install_package(pkg, force=force)
        except ValueError as msg:
            # An error occurred -- catch exception and record this in
            # the list of failed installed.
            print("*"*70)
            print("FAILED to install '{}'".format(pkg))
            print("*"*70)
            failed.append(pkg)
    return failed

def install_package(package=None, force=False):
    """
    Install a package or return a list of all packages that have been
    installed into this Sage install.

    You must have an internet connection. Also, you will have to
    restart Sage for the changes to take affect.

    It is not needed to provide the version number.

    INPUT:

    -  ``package`` -- optional; if specified, install the
       given package. If not, list all installed packages.

    - ``force`` -- boolean (default: ``False``); if ``True``, reinstall
      the package given if it is already installed. (Otherwise, an already
      installed package doesn't get reinstalled, as with 'sage -i ...'.)

    EXAMPLES:

    With no arguments, list the installed packages::

        sage: install_package()
        [...'atlas...'python...]

    With an argument, install the named package::

        sage: install_package('chomp')  # not tested
        Attempting to download package chomp-20100213.p2
        ...

    IMPLEMENTATION:

    Calls 'sage -f ...' to (re)install the package if a package name is
    given. If no package name is given, simply list the contents of
    ``spkg/installed``.

    .. seealso:: :func:`optional_packages`, :func:`upgrade`
    """
    global __installed_packages
    if os.uname()[0][:6] == 'CYGWIN' and package is not None:
        print("install_package may not work correctly under Microsoft Windows")
        print("since you can't change an opened file.  Quit all")
        print("instances of Sage and use 'sage -i {}' instead or".format(package))
        print("use the force option to install_package().")

    if package is None:
        if __installed_packages is None:
            import sage.env
            __installed_packages = sorted(os.listdir(sage.env.SAGE_SPKG_INST))
        return __installed_packages



    import stopgap
    stopgap.stopgap("The Sage function 'install_packages' is currently broken: it does not correctly install new packages. Please use 'sage -i {}' from a shell prompt instead.".format(package), 16759)
    return



    # Get full package name / list of candidates:
    if force:
        # Also search packages already installed.
        S = [P for P in standard_packages()[0] if P.startswith(package)]
        O = [P for P in optional_packages()[0] if P.startswith(package)]
        E = [P for P in experimental_packages()[0] if P.startswith(package)]
    else:
        S,O,E = [], [], []
    S.extend([P for P in standard_packages()[1] if P.startswith(package)])
    O.extend([P for P in optional_packages()[1] if P.startswith(package)])
    E.extend([P for P in experimental_packages()[1] if P.startswith(package)])
    L = S+O+E
    if len(L) > 1:
        if force:
            print("Possible package names starting with '{}' are:".format(package))
        else:
            print("Possible names of non-installed packages starting with '{}':".format(package))
        for P in L:
            print(" ", P)
        raise ValueError("There is more than one package name starting with '{}'. Please specify!".format(package))
    if len(L) == 0:
        if not force:
            if is_package_installed(package):
                raise ValueError("Package is already installed. Try install_package('{}',force=True)".format(package))
        raise ValueError("There is no package name starting with '{}'.".format(package))
    # len(L) == 1, i.e. exactly one package matches the given one.
    os.system('sage -f "{}"'.format(L[0]))
    __installed_packages = None


def is_package_installed(package):
    """
    Return true if a package starting with the given string is
    installed.

    EXAMPLES::

        sage: is_package_installed('sage')
        True
    """
    return any(p.startswith(package) for p in install_package())

def _package_lists_from_sage_output(package_type):
    r"""
    Helper function for :func:`standard_packages`, :func:`optional_packages` and
    :func:`experimental_packages`.

    INPUT:

    - ``package_type`` (string) -- one of `"standard"`, `"optional"` or
      `"experimental"`

    OUTPUT:

    A pair of list ``(installed,not_installed)`` with the corresponding
    packages' name.

    EXAMPLE::

        sage: from sage.misc.package import standard_packages # indirect doctest
        sage: installed, not_installed = standard_packages()  # optional internet
    """
    if package_type not in ['standard','optional','experimental']:
        raise ValueError("'package_type' must be one of 'standard','optional','experimental'.")

    R = os.popen('sage -'+package_type).read()
    X = R.split('\n')
    try:
        X = X[2+X.index((package_type+' packages:').capitalize()):]
        X = X[:X.index('')]
    except ValueError as msg:
        print(R)
        print("The list of %s packages (shown above) appears to be currently" % package_type)
        print("not available or corrupted (network error?).")
        return [], []

    def pname_from_line(line):
        pname = line.split(' ')[0]
        if '-' in pname:
            return pname[:pname.rfind('-')]
        else:
            return pname

    installed     = [pname_from_line(l) for l in X
                     if ('already installed' in l or
                         'installed version' in l)]

    not_installed = [pname_from_line(l) for l in X
                     if 'not installed' in l]

    if len(X) != len(installed) + len(not_installed):
        raise RuntimeError("Some package is missing from the list")

    return installed, not_installed


def standard_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed standard packages that are available
    from the Sage repository. You must have an internet connection.

    OUTPUT:

    -  installed standard packages (as a list)

    -  NOT installed standard packages (as a list)

    Use ``install_package(package_name)`` to install or
    re-install a given package.

    .. seealso:: :func:`install_package`, :func:`upgrade`

    EXAMPLE::

        sage: from sage.misc.package import standard_packages
        sage: installed, not_installed = standard_packages() # optional internet
        sage: installed[0], installed[-1]                    # optional internet
        ('atlas', 'zn_poly')
        sage: 'mercurial' in not_installed                   # optional internet
        True

    """
    return _package_lists_from_sage_output('standard')

def optional_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed optional packages that are available
    from the Sage repository. You must have an internet connection.

    OUTPUT:

    -  installed optional packages (as a list)

    -  NOT installed optional packages (as a list)


    Use ``install_package(package_name)`` to install or
    re-install a given package.

    .. seealso:: :func:`install_package`, :func:`upgrade`

    EXAMPLE::

        sage: from sage.misc.package import optional_packages
        sage: installed, not_installed = optional_packages() # optional internet
        sage: min(installed+not_installed)                   # optional internet
        '4ti2'
        sage: max(installed+not_installed)                   # optional internet
        'zeromq'
    """
    return _package_lists_from_sage_output('optional')

def experimental_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed experimental packages that are available
    from the Sage repository. You must have an internet connection.

    OUTPUT:

    -  installed experimental packages (as a list)

    -  NOT installed experimental packages (as a list)


    Use ``install_package(package_name)`` to install or
    re-install a given package.

    .. seealso:: :func:`install_package`, :func:`upgrade`

    EXAMPLE::

        sage: from sage.misc.package import experimental_packages
        sage: installed, not_installed = experimental_packages() # optional internet
        sage: min(installed+not_installed)                   # optional internet
        'PyQt4'
        sage: max(installed+not_installed)                   # optional internet
        'yassl'
    """
    return _package_lists_from_sage_output('experimental')

#################################################################
# Upgrade to latest version of Sage
#################################################################


def upgrade():
    """
    Download and build the latest version of Sage.

    You must have an internet connection. Also, you will have to
    restart Sage for the changes to take affect.

    This upgrades to the latest version of core packages (optional
    packages are not automatically upgraded).

    This will not work on systems that don't have a C compiler.

    .. seealso:: :func:`install_package`, :func:`optional_packages`
    """
    global __installed_packages
    if os.uname()[0][:6] == 'CYGWIN':
        print("Upgrade may not work correctly under Microsoft Windows")
        print("since you can't change an opened file.  Quit all")
        print("instances of Sage and use 'sage -upgrade' instead.")
        return

    os.system('sage -upgrade')
    __installed_packages = None
    print("You should quit and restart Sage now.")


def package_mesg(package_name):
    mesg  = 'To install the package %s type install_package("%s")\n'%(package_name, package_name)
    mesg += 'at the sage prompt.  Note, the version number might\n'
    mesg += 'change; if so, type optional_packages() to see a list \n'
    mesg += 'of possibilities.   All this requires an internet connection.'
    mesg += 'For more help, type optional_packages?'
    return mesg

