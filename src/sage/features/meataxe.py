# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of ``meataxe``
"""

# *****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#                     2021 Kwankyu Lee
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************


from . import PythonModule
from .join_feature import JoinFeature


class Meataxe(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the Sage modules
    that depend on the :ref:`meataxe <spkg_meataxe>` library.

    EXAMPLES::

        sage: from sage.features.meataxe import Meataxe
        sage: Meataxe().is_present()  # optional - meataxe
        FeatureTestResult('meataxe', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.meataxe import Meataxe
            sage: isinstance(Meataxe(), Meataxe)
            True
        """
        JoinFeature.__init__(self, 'meataxe',
                             [PythonModule('sage.matrix.matrix_gfpn_dense',
                                           spkg='sagemath_meataxe')])


def all_features():
    return [Meataxe()]
