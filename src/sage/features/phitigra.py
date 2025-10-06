# sage_setup: distribution = sagemath-environment
r"""
Check for ``phitigra``
"""

# *****************************************************************************
#       Copyright (C) 2022 Jean-Florent Raymond
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import PythonModule


class Phitigra(PythonModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of :ref:`phitigra <spkg_phitigra>`.

    Phitigra is provided by an optional package in the Sage distribution.

    EXAMPLES::

        sage: from sage.features.phitigra import Phitigra
        sage: Phitigra().is_present()                     # optional - phitigra
        FeatureTestResult('phitigra', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.phitigra import Phitigra
            sage: isinstance(Phitigra(), Phitigra)
            True
        """
        PythonModule.__init__(self, 'phitigra', spkg='phitigra')


def all_features():
    return [Phitigra()]
