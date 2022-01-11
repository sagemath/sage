r"""
Listing Sage packages

This module can be used to see which Sage packages are installed
and which packages are available for installation.

For more information about creating Sage packages, see the "Packaging
Third-Party Code" section of the Sage Developer's Guide.

Actually installing the packages should be done via the command
line, using the following commands:

- ``sage -i PACKAGE_NAME`` -- install the given package

- ``sage -f PACKAGE_NAME`` -- re-install the given package, even if it
  was already installed

To list the packages available, either use in a terminal one of ``sage
-standard``, ``sage -optional`` or ``sage -experimental``. Or the following
command inside Sage::

    sage: from sage.misc.package import list_packages
    sage: pkgs = list_packages(local=True)  # optional - build
    sage: sorted(pkgs.keys())  # optional - build, random
    ['4ti2',
     'alabaster',
     'arb',
     ...
     'zlib',
     'zn_poly']

Functions
---------
"""

# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from typing import Dict, List, NamedTuple, Optional, Union

import sage.env

import json
import os
import subprocess
import sys
from pathlib import Path
from urllib.request import urlopen
from urllib.error import URLError
from ssl import SSLContext

DEFAULT_PYPI = 'https://pypi.org/pypi'

def pkgname_split(name):
    r"""
    Split a pkgname into a list of strings, 'name, version'.

    For some packages, the version string might be empty.

    EXAMPLES::

        sage: from sage.misc.package import pkgname_split
        sage: pkgname_split('hello_world-1.2')
        ['hello_world', '1.2']
    """
    return (name.split('-',1) + [''])[:2]

def pip_remote_version(pkg, pypi_url=DEFAULT_PYPI, ignore_URLError=False):
    r"""
    Return the version of this pip package available on PyPI.

    INPUT:

    - ``pkg`` -- the package

    - ``pypi_url`` -- (string, default: standard PyPI url) an optional Python
      package repository to use

    - ``ignore_URLError`` -- (default: ``False``) if set to ``True`` then no
      error is raised if the connection fails and the function returns ``None``

    EXAMPLES:

    The following test does fail if there is no TLS support (see e.g.
    :trac:`19213`)::

        sage: from sage.misc.package import pip_remote_version
        sage: pip_remote_version('beautifulsoup4') # optional - internet # not tested
        '...'

    These tests are reliable since the tested package does not exist::

        sage: nap = 'hey_this_is_NOT_a_python_package'
        sage: pypi = 'http://this.is.not.pypi.com/'
        sage: pip_remote_version(nap, pypi_url=pypi, ignore_URLError=True) # optional - internet
        doctest:...: UserWarning: failed to fetch the version of
        pkg='hey_this_is_NOT_a_python_package' at
        http://this.is.not.pypi.com/.../json
        sage: pip_remote_version(nap, pypi_url=pypi, ignore_URLError=False) # optional - internet
        Traceback (most recent call last):
        ...
        HTTPError: HTTP Error 404: Not Found
    """
    url = '{pypi_url}/{pkg}/json'.format(pypi_url=pypi_url, pkg=pkg)

    try:
        f = urlopen(url, context=SSLContext())
        text = f.read()
        f.close()
    except URLError:
        if ignore_URLError:
            import warnings
            warnings.warn("failed to fetch the version of pkg={!r} at {}".format(pkg, url))
            return
        else:
            raise

    info = json.loads(text)
    stable_releases = [v for v in info['releases'] if 'a' not in v and 'b' not in v]
    return max(stable_releases)

def pip_installed_packages(normalization=None):
    r"""
    Return a dictionary `name->version` of installed pip packages.

    This command returns *all* pip-installed packages. Not only Sage packages.

    INPUT:

    - ``normalization`` -- (optional, default: ``None``) according to which rule to
      normalize the package name, either ``None`` (as is) or ``'spkg'`` (format
      as in the Sage distribution in ``build/pkgs/``), i.e., lowercased and
      dots and dashes replaced by underscores.

    EXAMPLES::

        sage: from sage.misc.package import pip_installed_packages
        sage: d = pip_installed_packages()  # optional - build
        sage: 'scipy' in d  # optional - build
        True
        sage: d['scipy']  # optional - build
        '...'
        sage: d['beautifulsoup4']   # optional - build beautifulsoup4
        '...'
        sage: d['prompt-toolkit']   # optional - build
        '...'
        sage: d = pip_installed_packages(normalization='spkg')  # optional - build
        sage: d['prompt_toolkit']   # optional - build
        '...'

    """
    with open(os.devnull, 'w') as devnull:
        proc = subprocess.Popen(
            [sys.executable, "-m", "pip", "list", "--no-index", "--format", "json"],
            stdout=subprocess.PIPE,
            stderr=devnull,
        )
        stdout = proc.communicate()[0].decode()

        def normalize(name: str) -> str:
            if normalization is None:
                return name
            elif normalization == 'spkg':
                return name.lower().replace('-', '_').replace('.', '_')
            else:
                raise NotImplementedError(f'normalization {normalization} is not implemented')
        try:
            return {normalize(package['name']): package['version']
                    for package in json.loads(stdout)}
        except json.decoder.JSONDecodeError:
            # Something went wrong while parsing the output from pip.
            # This may happen if pip is not correctly installed.
            return {}


class PackageInfo(NamedTuple):
    """Represents information about a package."""
    name: str
    type: Optional[str] = None
    source: Optional[str] = None
    installed_version: Optional[str] = None
    remote_version: Optional[str] = None

    def is_installed(self) -> bool:
        r"""
        Whether the package is installed in the system.
        """
        return self.installed_version is not None

    def __getitem__(self, key: Union[int, str]):
        r"""
        Only for backwards compatibility to allow dict-like access.

        TESTS::

            sage: from sage.misc.package import PackageInfo
            sage: package = PackageInfo("test_package")
            sage: package["name"]
            doctest:warning...
            dict-like access is deprecated, use pkg.name instead of pkg['name'], for example
            See https://trac.sagemath.org/31013 for details.
            'test_package'
            sage: package[0]
            'test_package'
        """
        if isinstance(key, str):
            from sage.misc.superseded import deprecation

            if key == "installed":
                deprecation(31013, "dict-like access via 'installed' is deprecated, use method is_installed instead")
                return self.is_installed()
            else:
                deprecation(31013, "dict-like access is deprecated, use pkg.name instead of pkg['name'], for example")
                return self.__getattribute__(key)
        else:
            return tuple.__getitem__(self, key)


def list_packages(*pkg_types: str, pkg_sources: List[str] = ['normal', 'pip', 'script'],
                  local: bool = False, ignore_URLError: bool = False, exclude_pip: bool = False) -> Dict[str, PackageInfo]:
    r"""
    Return a dictionary of information about each package.

    The keys are package names and values are named tuples with the following keys:

    - ``'type'``: either ``'base``, ``'standard'``, ``'optional'``, or ``'experimental'``
    - ``'source'``: either ``'normal', ``'pip'``, or ``'script'``
    - ``'installed'``: boolean
    - ``'installed_version'``: ``None`` or a string
    - ``'remote_version'``: string

    INPUT:

    - ``pkg_types`` -- (optional) a sublist of ``'base``, ``'standard'``, ``'optional'``,
      or ``'experimental'``.  If provided, list only the packages with the
      given type(s), otherwise list all packages.

    - ``pkg_sources`` -- (optional) a sublist of ``'normal', ``'pip'``, or ``'script'``.
      If provided, list only the packages with the given source(s), otherwise list all
      packages.

    - ``local`` -- (optional, default: ``False``) if set to ``True``, then do not
      consult remote (PyPI) repositories for package versions (only applicable for
      ``'pip'`` type)

    - ``exclude_pip`` -- (optional, default: ``False``) if set to ``True``, then
      pip packages are not considered.  This is the same as removing ``'pip'``
      from ``pkg_sources``.

    - ``ignore_URLError`` -- (default: ``False``) if set to ``True``, then
      connection errors will be ignored

    EXAMPLES::

        sage: from sage.misc.package import list_packages
        sage: L = list_packages('standard')    # optional - build
        sage: sorted(L.keys())                 # optional - build, random
        ['alabaster',
         'arb',
         'babel',
         ...
         'zn_poly']
        sage: sage_conf_info = L['sage_conf']  # optional - build
        sage: sage_conf_info.type              # optional - build
        'standard'
        sage: sage_conf_info.is_installed()    # optional - build
        True
        sage: sage_conf_info.source            # optional - build
        'script'

        sage: L = list_packages(pkg_sources=['pip'], local=True)  # optional - build internet
        sage: bs4_info = L['beautifulsoup4'] # optional - build internet
        sage: bs4_info.type                    # optional - build internet
        'optional'
        sage: bs4_info.source                  # optional - build internet
        'pip'

    Check the option ``exclude_pip``::

        sage: [p for p, d in list_packages('optional', exclude_pip=True).items()  # optional - build
        ....:  if d.source == 'pip']
        []
    """
    if not pkg_types:
        pkg_types = ('base', 'standard', 'optional', 'experimental')
    elif any(pkg_type not in ('base', 'standard', 'optional', 'experimental') for pkg_type in pkg_types):
        raise ValueError("Each pkg_type must be one of 'base', 'standard', 'optional', 'experimental'")

    if exclude_pip:
        pkg_sources = [s for s in pkg_sources if s != 'pip']

    pkgs = {p: PackageInfo(name=p, installed_version=v)
            for p, v in installed_packages('pip' not in pkg_sources).items()}

    # Add additional information based on Sage's package repository
    lp = []
    SAGE_PKGS = sage.env.SAGE_PKGS
    if not SAGE_PKGS:
        return pkgs

    try:
        lp = os.listdir(SAGE_PKGS)
    except FileNotFoundError:
        return pkgs

    for p in lp:

        try:
            f = open(os.path.join(SAGE_PKGS, p, "type"))
        except IOError:
            # Probably an empty directory => ignore
            continue

        with f:
            typ = f.read().strip()

        if os.path.isfile(os.path.join(SAGE_PKGS, p, "requirements.txt")):
            src = 'pip'
        elif os.path.isfile(os.path.join(SAGE_PKGS, p, "checksums.ini")):
            src = 'normal'
        else:
            src = 'script'

        if typ not in pkg_types or src not in pkg_sources:
            try:
                del pkgs[p]
            except KeyError:
                pass
            continue

        if src == 'pip':
            if not local:
                remote_version = pip_remote_version(p, ignore_URLError=ignore_URLError)
            else:
                remote_version = None
        elif src == 'normal':
            # If package-version.txt does not exist, that is an error
            # in the build system => we just propagate the exception
            package_filename = os.path.join(SAGE_PKGS, p, "package-version.txt")
            with open(package_filename) as f:
                remote_version = f.read().strip()
        else:
            remote_version = None

        pkg = pkgs.get(p, PackageInfo(name=p))
        pkgs[p] = PackageInfo(p, typ, src, pkg.installed_version, remote_version)

    return pkgs

def _spkg_inst_dirs():
    """
    Generator for the installation manifest directories as resolved paths.

    It yields first ``SAGE_SPKG_INST``, then ``SAGE_VENV_SPKG_INST``,
    if defined; but it both resolve to the same directory, it only yields
    one element.

    EXAMPLES::

        sage: from sage.misc.package import _spkg_inst_dirs
        sage: list(_spkg_inst_dirs())
        [...]

    """
    last_inst_dir = None
    for inst_dir in (sage.env.SAGE_SPKG_INST, sage.env.SAGE_VENV_SPKG_INST):
        if inst_dir:
            inst_dir = Path(inst_dir).resolve()
            if inst_dir.is_dir() and inst_dir != last_inst_dir:
                yield inst_dir
                last_inst_dir = inst_dir

def installed_packages(exclude_pip=True):
    """
    Return a dictionary of all installed packages, with version numbers.

    INPUT:

    - ``exclude_pip`` -- (optional, default: ``True``) whether "pip" packages
      are excluded from the list

    EXAMPLES::

        sage: sorted(installed_packages().keys())  # optional - build
        [...'gmpy2', ...'sage_conf', ...]
        sage: installed_packages()['gmpy2']  # optional - build, random
        '2.1.0b5'

    .. SEEALSO::

        :func:`sage.misc.package.list_packages`
    """
    installed = {}
    if not exclude_pip:
        installed.update(pip_installed_packages(normalization='spkg'))
    # Sage packages should override pip packages (Trac #23997)

    for inst_dir in _spkg_inst_dirs():
        try:
            lp = os.listdir(inst_dir)
            installed.update(pkgname_split(pkgname) for pkgname in lp
                             if not pkgname.startswith('.'))
        except FileNotFoundError:
            pass
    return installed


def is_package_installed(package, exclude_pip=True):
    """
    Return whether (any version of) ``package`` is installed.

    INPUT:

    - ``package`` -- the name of the package

    - ``exclude_pip`` -- (optional, default: ``True``) whether to consider pip
      type packages


    EXAMPLES::

        sage: is_package_installed('gap')  # optional - build
        True

    Giving just the beginning of the package name is not good enough::

        sage: is_package_installed('matplotli')  # optional - build
        False

    Otherwise, installing "pillow" would cause this function to think
    that "pil" is installed, for example.

    .. NOTE::

        Do not use this function to check whether you can use a feature from an
        external library. This only checks whether something was installed with
        ``sage -i`` but it may have been installed by other means (for example
        if this copy of Sage has been installed as part of a distribution.)
        Use the framework provided by :mod:`sage.features` to check
        whether a library is installed and functional.
    """
    return any(p == package for p in installed_packages(exclude_pip))


def package_versions(package_type, local=False):
    r"""
    Return version information for each Sage package.

    INPUT:

    - ``package_type`` -- (string) one of ``"standard"``, ``"optional"`` or
      ``"experimental"``

    - ``local`` -- (boolean, default: ``False``) only query local data (no internet needed)

    For packages of the given type, return a dictionary whose entries
    are of the form ``'package': (installed, latest)``, where
    ``installed`` is the installed version (or ``None`` if not
    installed) and ``latest`` is the latest available version. If the
    package has a directory in ``SAGE_ROOT/build/pkgs/``, then
    ``latest`` is determined by the file ``package-version.txt`` in
    that directory.  If ``local`` is ``False``, then Sage's servers are
    queried for package information.

    .. SEEALSO:: :func:`sage.misc.package.list_packages`

    EXAMPLES::

        sage: std = package_versions('standard', local=True)  # optional - build
        sage: 'gap' in std  # optional - build
        True
        sage: std['zn_poly']  # optional - build, random
        ('0.9.p12', '0.9.p12')
    """
    return {pkg.name: (pkg.installed_version, pkg.remote_version) for pkg in list_packages(package_type, local=local).values()}


def standard_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed standard packages that are available
    from the Sage repository.

    OUTPUT:

    - installed standard packages (as a list)

    - NOT installed standard packages (as a list)

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

    .. SEEALSO:: :func:`sage.misc.package.list_packages`

    EXAMPLES::

        sage: from sage.misc.package import standard_packages
        sage: installed, not_installed = standard_packages()  # optional - build
        doctest:...: DeprecationWarning: ...
        sage: 'numpy' in installed                            # optional - build
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(30747,
                'the functions standard_packages, optional_packages, experimental_packages '
                'are deprecated, use sage.features instead')
    pkgs = list_packages('standard', local=True).values()
    return (sorted(pkg.name for pkg in pkgs if pkg.is_installed()),
            sorted(pkg.name for pkg in pkgs if not pkg.is_installed()))


def optional_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed optional packages that are available
    from the Sage repository.

    OUTPUT:

    - installed optional packages (as a list)

    - NOT installed optional packages (as a list)

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

    .. SEEALSO:: :func:`sage.misc.package.list_packages`

    EXAMPLES::

        sage: from sage.misc.package import optional_packages
        sage: installed, not_installed = optional_packages()  # optional - build
        doctest:...: DeprecationWarning: ...
        sage: 'beautifulsoup4' in installed+not_installed  # optional - build
        True

        sage: 'beautifulsoup4' in installed   # optional - build beautifulsoup4
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(30747,
                'the functions standard_packages, optional_packages, experimental_packages '
                'are deprecated, use sage.features instead')
    pkgs = list_packages('optional', local=True)
    pkgs = pkgs.values()
    return (sorted(pkg.name for pkg in pkgs if pkg.is_installed()),
            sorted(pkg.name for pkg in pkgs if not pkg.is_installed()))


def experimental_packages():
    """
    Return two lists. The first contains the installed and the second
    contains the not-installed experimental packages that are available
    from the Sage repository.

    OUTPUT:

    - installed experimental packages (as a list)

    - NOT installed experimental packages (as a list)

    Run ``sage -i package_name`` from a shell to install a given
    package or ``sage -f package_name`` to re-install it.

    .. SEEALSO:: :func:`sage.misc.package.list_packages`

    EXAMPLES::

        sage: from sage.misc.package import experimental_packages
        sage: installed, not_installed = experimental_packages()  # optional - build
        doctest:...: DeprecationWarning: ...
    """
    from sage.misc.superseded import deprecation
    deprecation(30747,
                'the functions standard_packages, optional_packages, experimental_packages '
                'are deprecated, use sage.features instead')
    pkgs = list_packages('experimental', local=True).values()
    return (sorted(pkg.name for pkg in pkgs if pkg.is_installed()),
            sorted(pkg.name for pkg in pkgs if not pkg.is_installed()))

def package_manifest(package):
    """
    Return the manifest for ``package``.

    INPUT:

    - ``package`` -- package name

    The manifest is written in the file
    ``SAGE_SPKG_INST/package-VERSION``. It is a JSON file containing a
    dictionary with the package name, version, installation date, list
    of installed files, etc.

    EXAMPLES::

        sage: from sage.misc.package import package_manifest
        sage: sagetex_manifest = package_manifest('sagetex')  # optional - build
        sage: sagetex_manifest['package_name'] == 'sagetex'  # optional - build
        True
        sage: 'files' in sagetex_manifest  # optional - build
        True

    Test a nonexistent package::

        sage: package_manifest('dummy-package')  # optional - build
        Traceback (most recent call last):
        ...
        KeyError: 'dummy-package'
    """
    version = installed_packages()[package]
    for inst_dir in _spkg_inst_dirs():
        stamp_file = os.path.join(inst_dir,
                                  '{}-{}'.format(package, version))
        try:
            with open(stamp_file) as f:
                return json.load(f)
        except FileNotFoundError:
            pass
    raise RuntimeError('package manifest directory changed at runtime')

class PackageNotFoundError(RuntimeError):
    """
    This class defines the exception that should be raised when a
    function, method, or class cannot detect a Sage package that it
    depends on.

    This exception should be raised with a single argument, namely
    the name of the package.

    When a ``PackageNotFoundError`` is raised, this means one of the
    following:

    - The required optional package is not installed.

    - The required optional package is installed, but the relevant
      interface to that package is unable to detect the package.

    Raising a ``PackageNotFoundError`` is deprecated.  Use
    :class:`sage.features.FeatureNotPresentError` instead.

    User code can continue to catch ``PackageNotFoundError`` exceptions
    for compatibility with older versions of the Sage library.
    This does not cause deprecation warnings.

    EXAMPLES::

        sage: from sage.misc.package import PackageNotFoundError
        sage: try:
        ....:     pass
        ....: except PackageNotFoundError:
        ....:     pass

    """

    def __init__(self, *args):
        """
        TESTS::

            sage: from sage.misc.package import PackageNotFoundError
            sage: raise PackageNotFoundError("my_package")
            Traceback (most recent call last):
            ...
            PackageNotFoundError: the package 'my_package' was not found. You can install it by running 'sage -i my_package' in a shell
        """
        super().__init__(*args)
        # We do not deprecate the whole class because we want
        # to allow user code to handle this exception without causing
        # a deprecation warning.
        from sage.misc.superseded import deprecation
        deprecation(30607, "Instead of raising PackageNotFoundError, raise sage.features.FeatureNotPresentError")

    def __str__(self):
        """
        Return the actual error message.

        EXAMPLES::

            sage: from sage.misc.package import PackageNotFoundError
            sage: str(PackageNotFoundError("my_package"))
            doctest:warning...
            "the package 'my_package' was not found. You can install it by running 'sage -i my_package' in a shell"
        """
        return ("the package {0!r} was not found. "
            "You can install it by running 'sage -i {0}' in a shell"
            .format(self.args[0]))
