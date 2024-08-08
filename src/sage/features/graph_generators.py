# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of graph generator programs ``benzene``, ``buckygen``, ``plantri``
"""

# *****************************************************************************
#       Copyright (C) 2016 Julian Rüth
#                     2018 Jeroen Demeyer
#                     2019 Frédéric Chapoton
#                     2021 Matthias Koeppe
#                     2021 Kwankyu Lee
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import os
import subprocess

from . import Executable, FeatureTestResult


class Plantri(Executable):
    r"""
    A :class:`~sage.features.Feature` which checks for the :ref:`plantri <spkg_plantri>` binary.

    EXAMPLES::

        sage: from sage.features.graph_generators import Plantri
        sage: Plantri().is_present()  # optional - plantri
        FeatureTestResult('plantri', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graph_generators import Plantri
            sage: isinstance(Plantri(), Plantri)
            True
        """
        Executable.__init__(self, name='plantri', spkg='plantri',
                            executable='plantri',
                            url='http://users.cecs.anu.edu.au/~bdm/plantri/')

    def is_functional(self):
        r"""
        Check whether ``plantri`` works on trivial input.

        EXAMPLES::

            sage: from sage.features.graph_generators import Plantri
            sage: Plantri().is_functional()  # optional - plantri
            FeatureTestResult('plantri', True)
        """
        command = ["plantri", "4"]
        try:
            lines = subprocess.check_output(command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                    reason="Call `{command}` failed with exit code {e.returncode}".format(command=" ".join(command), e=e))

        expected = b"1 triangulations written"
        if lines.find(expected) == -1:
            return FeatureTestResult(self, False,
                    reason="Call `{command}` did not produce output which contains `{expected}`".format(command=" ".join(command), expected=expected))

        return FeatureTestResult(self, True)


class Buckygen(Executable):
    r"""
    A :class:`~sage.features.Feature` which checks for the :ref:`buckygen <spkg_buckygen>` binary.

    EXAMPLES::

        sage: from sage.features.graph_generators import Buckygen
        sage: Buckygen().is_present()  # optional - buckygen
        FeatureTestResult('buckygen', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graph_generators import Buckygen
            sage: isinstance(Buckygen(), Buckygen)
            True
        """
        Executable.__init__(self, name='buckygen', spkg='buckygen',
                            executable='buckygen',
                            url='http://caagt.ugent.be/buckygen/')

    def is_functional(self):
        r"""
        Check whether ``buckygen`` works on trivial input.

        EXAMPLES::

            sage: from sage.features.graph_generators import Buckygen
            sage: Buckygen().is_functional()  # optional - buckygen
            FeatureTestResult('buckygen', True)
        """
        command = ["buckygen", "-d", "22d"]
        try:
            lines = subprocess.check_output(command, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                    reason="Call `{command}` failed with exit code {e.returncode}".format(command=" ".join(command), e=e))

        expected = b"Number of fullerenes generated with 13 vertices: 0"
        if lines.find(expected) == -1:
            return FeatureTestResult(self, False,
                    reason="Call `{command}` did not produce output which contains `{expected}`".format(command=" ".join(command), expected=expected))

        return FeatureTestResult(self, True)


class Benzene(Executable):
    r"""
    A :class:`~sage.features.Feature` which checks for the :ref:`benzene <spkg_benzene>`
    binary.

    EXAMPLES::

        sage: from sage.features.graph_generators import Benzene
        sage: Benzene().is_present()  # optional - benzene
        FeatureTestResult('benzene', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graph_generators import Benzene
            sage: isinstance(Benzene(), Benzene)
            True
        """
        Executable.__init__(self, name='benzene', spkg='benzene',
                            executable='benzene',
                            url='http://www.grinvin.org/')

    def is_functional(self):
        r"""
        Check whether ``benzene`` works on trivial input.

        EXAMPLES::

            sage: from sage.features.graph_generators import Benzene
            sage: Benzene().is_functional()  # optional - benzene
            FeatureTestResult('benzene', True)
        """
        devnull = open(os.devnull, 'wb')
        command = ["benzene", "2", "p"]
        try:
            lines = subprocess.check_output(command, stderr=devnull)
        except subprocess.CalledProcessError as e:
            return FeatureTestResult(self, False,
                    reason="Call `{command}` failed with exit code {e.returncode}".format(command=" ".join(command), e=e))

        expected = b">>planar_code<<"
        if not lines.startswith(expected):
            return FeatureTestResult(self, False,
                    reason="Call `{command}` did not produce output that started with `{expected}`.".format(command=" ".join(command), expected=expected))

        return FeatureTestResult(self, True)


def all_features():
    return [Plantri(),
            Buckygen(),
            Benzene()]
