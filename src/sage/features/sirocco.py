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
from sage.features.build_feature import BuildModule

class Sirocco(BuildModule):
    r"""
    A :class:`~sage.features.Feature` which describes whether the
    :mod:`sage.libs.sirocco` module is available in this installation
    of Sage.

    EXAMPLES::

        sage: from sage.features.sirocco import Sirocco
        sage: Sirocco().is_present()  # needs sirocco
        FeatureTestResult('sirocco', True)
        sage: Sirocco().is_present()  # needs !sirocco
        FeatureTestResult('sirocco', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !sirocco`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.sirocco import Sirocco
        sage: Sirocco().is_present_at_runtime()  # needs sirocco
        FeatureTestResult('sirocco', True)

    """
    _enabled_in_build = sirocco_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.sirocco import Sirocco
            sage: Sirocco()
            Feature('sirocco')

        """
        module_name = "sage.libs.sirocco"
        super().__init__("sirocco", module_name, spkg="sirocco")


def all_features():
    return [Sirocco()]
