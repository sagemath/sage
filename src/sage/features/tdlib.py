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
from sage.features.build_feature import BuildModule

class Tdlib(BuildModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    the SageMath interface to the treedec (formerly tdlib) library.

    EXAMPLES::

        sage: from sage.features.tdlib import Tdlib
        sage: Tdlib().is_present()  # needs tdlib
        FeatureTestResult('tdlib', True)
        sage: Tdlib().is_present()  # needs !tdlib
        FeatureTestResult('tdlib', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !tdlib`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.tdlib import Tdlib
        sage: Tdlib().is_present_at_runtime()  # needs tdlib
        FeatureTestResult('tdlib', True)

    """
    _enabled_in_build = tdlib_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.tdlib import Tdlib
            sage: Tdlib()
            Feature('tdlib')

        """
        module_name = "sage.graphs.graph_decompositions.tdlib"
        super().__init__("tdlib", module_name, spkg="tdlib")


def all_features():
    return [Tdlib()]
