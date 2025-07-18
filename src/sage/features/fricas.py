# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of ``fricas``
"""

# *****************************************************************************
#       Copyright (C) 2023 Dima Pasechnik
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import subprocess
from . import Executable, FeatureTestResult
from packaging.version import Version


class FriCAS(Executable):
    r"""
    A :class:`~sage.features.Feature` which checks for the :ref:`fricas <fricas>` binary.

    EXAMPLES::

        sage: from sage.features.fricas import FriCAS
        sage: FriCAS().is_present()  # optional - fricas
        FeatureTestResult('fricas', True)
    """
    MINIMUM_VERSION = "1.3.8"

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.fricas import FriCAS
            sage: isinstance(FriCAS(), FriCAS)
            True
        """
        Executable.__init__(self, name='fricas', spkg='fricas',
                            executable='fricas',
                            url='https://fricas.github.io')

    def get_version(self):
        r"""
        Retrieve the installed FriCAS version

        EXAMPLES::
            sage: from sage.features.fricas import FriCAS
            sage: FriCAS().get_version() # optional - fricas
            '1.3...'
        """
        try:
            output = subprocess.check_output(['fricas', '--version'], stderr=subprocess.STDOUT)
            version_line = output.decode('utf-8').strip()
            version = version_line.split()[1]
            return version
        except subprocess.CalledProcessError:
            return None

    def is_functional(self):
        r"""
        Check whether ``fricas`` works on trivial input.

        EXAMPLES::

            sage: from sage.features.fricas import FriCAS
            sage: FriCAS().is_functional()  # optional - fricas
            FeatureTestResult('fricas', True)
        """
        command = ['fricas -nosman -eval ")quit"']
        try:
            lines = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                                     reason="Call `{command}` failed with exit code {e.returncode}".format(command=" ".join(command), e=e))

        expected = b"FriCAS"
        if lines.find(expected) == -1:
            return FeatureTestResult(self, False,
                                     reason="Call `{command}` did not produce output which contains `{expected}`".format(command=" ".join(command),
                                                                                                                         expected=expected))
        version = self.get_version()
        if version is None:
            return FeatureTestResult(self, False,
                                     reason="Could not determine FriCAS version")

        try:
            if Version(version) < Version(self.MINIMUM_VERSION):
                return FeatureTestResult(self, False, reason=f"FriCAS version {version} is too old; minimum required is {self.MINIMUM_VERSION}")
            return FeatureTestResult(self, True)
        except ValueError:
            return FeatureTestResult(self, False, reason="Invalid Version Format")


def all_features():
    return [FriCAS()]
