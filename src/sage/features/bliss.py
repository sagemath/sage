r"""
Features for testing the presence of ``bliss``
"""

# *****************************************************************************
#       Copyright (C) 2016 Julian Rüth
#                     2018 Jeroen Demeyer
#                     2021 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.config import bliss_enabled
from sage.features.build_feature import BuildModule

class Bliss(BuildModule):
    r"""
    A :class:`~sage.features.Feature` which describes whether the
    :mod:`sage.graphs.bliss` module is available in this installation
    of Sage.

    EXAMPLES::

        sage: from sage.features.bliss import Bliss
        sage: Bliss().is_present()  # needs bliss
        FeatureTestResult('bliss', True)
        sage: Bliss().is_present()  # needs !bliss
        FeatureTestResult('bliss', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !bliss`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.bliss import Bliss
        sage: Bliss().is_present_at_runtime()  # needs bliss
        FeatureTestResult('bliss', True)

    """
    _enabled_in_build = bliss_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.bliss import Bliss
            sage: Bliss()
            Feature('bliss')

        """
        module_name = "sage.graphs.bliss"
        super().__init__("bliss",
                         module_name,
                         url='http://www.tcs.hut.fi/Software/bliss/')

def all_features():
    return [Bliss()]
