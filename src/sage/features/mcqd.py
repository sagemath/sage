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
from sage.features.build_feature import BuildModule

class Mcqd(BuildModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    the ``sage.graphs.mcqd`` module, which is the SageMath
    interface to the :ref:`mcqd <spkg_mcqd>` library

    EXAMPLES::

        sage: from sage.features.mcqd import Mcqd
        sage: Mcqd().is_present()  # needs mcqd
        FeatureTestResult('mcqd', True)
        sage: Mcqd().is_present()  # needs !mcqd
        FeatureTestResult('mcqd', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !mcqd`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.mcqd import Mcqd
        sage: Mcqd().is_present_at_runtime()  # needs mcqd
        FeatureTestResult('mcqd', True)

    """
    _enabled_in_build = mcqd_enabled

    def __init__(self):
        """
        EXAMPLES::

            sage: from sage.features.mcqd import Mcqd
            sage: Mcqd()
            Feature('mcqd')

        """
        module_name = "sage.graphs.mcqd"
        super().__init__("mcqd", module_name)


def all_features():
    return [Mcqd()]
