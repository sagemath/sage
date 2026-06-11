r"""
Features for testing the presence of ``coxeter3``
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

from sage.config import coxeter3_enabled
from sage.features import PythonModule
from sage.features.build_feature import BuildFeature

class Coxeter3(BuildFeature):
    r"""
    A :class:`~sage.features.Feature` which describes whether the
    :mod:`sage.libs.coxeter3` module is available in this installation
    of Sage.

    EXAMPLES::

        sage: from sage.features.coxeter3 import Coxeter3
        sage: Coxeter3().require()  # needs coxeter3

    """
    _enabled_in_build = coxeter3_enabled

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.coxeter3 import Coxeter3
            sage: Coxeter3()
            Feature('coxeter3')

        """
        super().__init__("coxeter3", spkg="coxeter3")

    def is_present_at_runtime(self):
        r"""
        TESTS::

            sage: from sage.features import FeatureTestResult
            sage: from sage.features.coxeter3 import Coxeter3
            sage: result = Coxeter3().is_present_at_runtime()
            sage: isinstance(result, FeatureTestResult)
            True
            sage: result  # needs coxeter3
            FeatureTestResult('coxeter3', True)

        """
        # The build system installs the sage/libs/coxeter3 source code
        # even when coxeter3 support is disabled, so a naive import of
        # that module will actually succeed. We check for a
        # conditionally-compiled cython module instead.
        cython_modname = "sage.libs.coxeter3.coxeter"
        result = PythonModule(cython_modname)._is_present()
        result.feature = self
        return result

def all_features():
    return [Coxeter3()]
