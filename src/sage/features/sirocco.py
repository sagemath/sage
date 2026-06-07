r"""
Features for testing the presence of ``sirocco``
"""

# *****************************************************************************
#       Copyright (C) 2016      Julian Rüth
#                     2018      Jeroen Demeyer
#                     2021-2024 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.config import sirocco_enabled
from sage.features import PythonModule
from sage.features.build_feature import BuildFeature

class Sirocco(BuildFeature):
    r"""
    A :class:`~sage.features.Feature` which describes whether the
    :mod:`sage.libs.sirocco` module is available in this installation
    of Sage.

    EXAMPLES::

        sage: from sage.features.sirocco import Sirocco
        sage: Sirocco().require()  # needs sirocco

    """
    _enabled_in_build = sirocco_enabled

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sirocco import Sirocco
            sage: Sirocco()
            Feature('sirocco')

        """
        super().__init__("sirocco", spkg="sirocco")

    def is_present_at_runtime(self):
        r"""
        TESTS::

            sage: from sage.features import FeatureTestResult
            sage: from sage.features.sirocco import Sirocco
            sage: result = Sirocco().is_present_at_runtime()
            sage: isinstance(result, FeatureTestResult)
            True
            sage: result  # needs sirocco
            FeatureTestResult('sirocco', True)

        """
        result = PythonModule("sage.libs.sirocco")._is_present()
        result.feature = self
        return result

def all_features():
    return [Sirocco()]
