# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of package systems ``sage_spkg``, ``conda``, ``pip``, ``debian``, ``fedora`` etc.
"""

# *****************************************************************************
#       Copyright (C) 2021-2022 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import Feature


class PackageSystem(Feature):
    r"""
    A :class:`Feature` describing a system package manager.

    EXAMPLES::

        sage: from sage.features.pkg_systems import PackageSystem
        sage: PackageSystem('conda')
        Feature('conda')
    """
    def _is_present(self):
        r"""
        Test whether ``self`` appears in the list of available package systems.

        EXAMPLES::

            sage: from sage.features.pkg_systems import PackageSystem
            sage: debian = PackageSystem('debian')
            sage: debian.is_present()  # indirect doctest, random
            True
        """
        from . import package_systems
        return self in package_systems()

    def spkg_installation_hint(self, spkgs, *, prompt="  !", feature=None):
        r"""
        Return a string that explains how to install ``feature``.

        EXAMPLES::

            sage: from sage.features.pkg_systems import PackageSystem
            sage: homebrew = PackageSystem('homebrew')
            sage: homebrew.spkg_installation_hint('openblas')  # optional - SAGE_ROOT
            'To install openblas using the homebrew package manager, you can try to run:\n!brew install openblas'
        """
        if isinstance(spkgs, (tuple, list)):
            spkgs = ' '.join(spkgs)
        if feature is None:
            feature = spkgs
        return self._spkg_installation_hint(spkgs, prompt, feature)

    def _spkg_installation_hint(self, spkgs, prompt, feature):
        r"""
        Return a string that explains how to install ``feature``.

        Override this method in derived classes.

        EXAMPLES::

            sage: from sage.features.pkg_systems import PackageSystem
            sage: fedora = PackageSystem('fedora')
            sage: fedora.spkg_installation_hint('openblas')  # optional - SAGE_ROOT
            'To install openblas using the fedora package manager, you can try to run:\n!sudo yum install openblas-devel'
        """
        from subprocess import run, CalledProcessError, PIPE
        lines = []
        system = self.name
        try:
            proc = run(f'sage-get-system-packages {system} {spkgs}',
                       shell=True, capture_output=True, text=True, check=True)
            system_packages = proc.stdout.strip()
            print_sys = f'sage-print-system-package-command {system} --verbose --sudo --prompt="{prompt}"'
            command = f'{print_sys} update && {print_sys} install {system_packages}'
            proc = run(command, shell=True, capture_output=True, text=True, check=True)
            command = proc.stdout.strip()
            if command:
                lines.append(f'To install {feature} using the {system} package manager, you can try to run:')
                lines.append(command)
                return '\n'.join(lines)
        except CalledProcessError:
            pass
        return f'No equivalent system packages for {system} are known to Sage.'


class SagePackageSystem(PackageSystem):
    r"""
    A :class:`Feature` describing the package manager of the SageMath distribution.

    EXAMPLES::

        sage: from sage.features.pkg_systems import SagePackageSystem
        sage: SagePackageSystem()
        Feature('sage_spkg')
    """
    @staticmethod
    def __classcall__(cls):
        r"""
        Normalize initargs.

        TESTS::

            sage: from sage.features.pkg_systems import SagePackageSystem
            sage: SagePackageSystem() is SagePackageSystem()  # indirect doctest
            True
        """
        return PackageSystem.__classcall__(cls, "sage_spkg")

    def _is_present(self):
        r"""
        Test whether ``sage-spkg`` is available.

        EXAMPLES::

            sage: from sage.features.pkg_systems import SagePackageSystem
            sage: bool(SagePackageSystem().is_present())  # indirect doctest, optional - sage_spkg
            True
        """
        from subprocess import run, DEVNULL, CalledProcessError
        try:
            # "sage -p" is a fast way of checking whether sage-spkg is available.
            run('sage -p', shell=True, stdout=DEVNULL, stderr=DEVNULL, check=True)
        except CalledProcessError:
            return False
        # Check if there are any installation records.
        try:
            from sage.misc.package import installed_packages
        except ImportError:
            return False
        for pkg in installed_packages(exclude_pip=True):
            return True
        return False

    def _spkg_installation_hint(self, spkgs, prompt, feature):
        r"""
        Return a string that explains how to install ``feature``.

        EXAMPLES::

            sage: from sage.features.pkg_systems import SagePackageSystem
            sage: print(SagePackageSystem().spkg_installation_hint(['foo', 'bar'], prompt="### ", feature='foobarability'))  # indirect doctest
            To install foobarability using the Sage package manager, you can try to run:
            ### sage -i foo bar
        """
        lines = []
        lines.append(f'To install {feature} using the Sage package manager, you can try to run:')
        lines.append(f'{prompt}sage -i {spkgs}')
        return '\n'.join(lines)


class PipPackageSystem(PackageSystem):
    r"""
    A :class:`Feature` describing the Pip package manager.

    EXAMPLES::

        sage: from sage.features.pkg_systems import PipPackageSystem
        sage: PipPackageSystem()
        Feature('pip')
    """
    @staticmethod
    def __classcall__(cls):
        r"""
        Normalize initargs.

        TESTS::

            sage: from sage.features.pkg_systems import PipPackageSystem
            sage: PipPackageSystem() is PipPackageSystem()  # indirect doctest
            True
        """
        return PackageSystem.__classcall__(cls, "pip")

    def _is_present(self):
        r"""
        Test whether ``pip`` is available.

        EXAMPLES::

            sage: from sage.features.pkg_systems import PipPackageSystem
            sage: bool(PipPackageSystem().is_present())    # indirect doctest
            True
        """
        from subprocess import run, DEVNULL, CalledProcessError
        try:
            run('sage -pip --version', shell=True, stdout=DEVNULL, stderr=DEVNULL, check=True)
            return True
        except CalledProcessError:
            return False
