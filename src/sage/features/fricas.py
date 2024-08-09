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

import os
import subprocess
from . import Executable, FeatureTestResult

class FriCAS(Executable):
    r"""
    A :class:`~sage.features.Feature` which checks for the :ref:`fricas <fricas>` binary.

    EXAMPLES::

        sage: from sage.features.fricas import FriCAS
        sage: FriCAS().is_present()  # optional - fricas
        FeatureTestResult('fricas', True)
    """
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
                    reason="Call `{command}` did not produce output which contains `{expected}`".format(command=" ".join(command), expected=expected))

        return FeatureTestResult(self, True)

def all_features():
    return [FriCAS()]
