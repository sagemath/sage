r"""
Features for testing the presence of ``tdlib``
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

from sage.config import tdlib_enabled
from sage.features import PythonModule
from sage.features.build_feature import BuildFeature

class Tdlib(BuildFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    the SageMath interface to the :ref:`tdlib <spkg_tdlib>` library.
    """
    _enabled_in_build = tdlib_enabled

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.tdlib import Tdlib
            sage: isinstance(Tdlib(), Tdlib)
            True

        """
        super().__init__("tdlib", spkg="tdlib")

    def is_present_at_runtime(self):
        r"""
        TESTS::

            sage: from sage.features import FeatureTestResult
            sage: from sage.features.tdlib import Tdlib
            sage: result = Tdlib().is_present_at_runtime()
            sage: isinstance(result, FeatureTestResult)
            True
            sage: result  # needs tdlib
            FeatureTestResult('tdlib', True)

        """
        modname = "sage.graphs.graph_decompositions.tdlib"
        result = PythonModule(modname)._is_present()
        result.feature = self
        return result

def all_features():
    return [Tdlib()]
