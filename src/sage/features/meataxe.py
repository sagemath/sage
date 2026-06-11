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

from sage.config import meataxe_enabled
from sage.features import PythonModule
from sage.features.build_feature import BuildFeature

class Meataxe(BuildFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    the Sage modules that depend on the :ref:`meataxe <spkg_meataxe>`
    library.

    EXAMPLES::

        sage: from sage.features.meataxe import Meataxe
        sage: Meataxe().is_present()  # needs meataxe
        FeatureTestResult('meataxe', True)

    """
    _enabled_in_build = meataxe_enabled

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.meataxe import Meataxe
            sage: isinstance(Meataxe(), Meataxe)
            True

        """
        super().__init__("meataxe", spkg="meataxe")

    def is_present_at_runtime(self):
        r"""
        TESTS::

            sage: from sage.features import FeatureTestResult
            sage: from sage.features.meataxe import Meataxe
            sage: result = Meataxe().is_present_at_runtime()
            sage: isinstance(result, FeatureTestResult)
            True
            sage: result  # needs meataxe
            FeatureTestResult('meataxe', True)

        """
        modname = "sage.matrix.matrix_gfpn_dense"
        result = PythonModule(modname)._is_present()
        result.feature = self
        return result

def all_features():
    return [Meataxe()]
