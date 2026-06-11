r"""
Features for testing the presence of ``mcqd``
"""

# *****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.config import mcqd_enabled
from sage.features import PythonModule
from sage.features.build_feature import BuildFeature

class Mcqd(BuildFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    the :mod:`~sage.graphs.mcqd` module, which is the SageMath
    interface to the :ref:`mcqd <spkg_mcqd>` library

    EXAMPLES::

        sage: from sage.features.mcqd import Mcqd
        sage: Mcqd().is_present()  # needs mcqd
        FeatureTestResult('mcqd', True)

    """
    _enabled_in_build = mcqd_enabled

    def __init__(self):
        """
        TESTS::

            sage: from sage.features.mcqd import Mcqd
            sage: isinstance(Mcqd(), Mcqd)
            True

        """
        super().__init__("mcqd", spkg="mcqd")

    def is_present_at_runtime(self):
        r"""
        TESTS::

            sage: from sage.features import FeatureTestResult
            sage: from sage.features.mcqd import Mcqd
            sage: result = Mcqd().is_present_at_runtime()
            sage: isinstance(result, FeatureTestResult)
            True
            sage: result  # needs mcqd
            FeatureTestResult('mcqd', True)

        """
        result = PythonModule("sage.graphs.mcqd")._is_present()
        result.feature = self
        return result

def all_features():
    return [Mcqd()]
