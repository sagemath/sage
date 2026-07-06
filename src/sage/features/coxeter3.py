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
from sage.features.build_feature import BuildModule

class Coxeter3(BuildModule):
    r"""
    A :class:`~sage.features.Feature` which describes whether the
    :mod:`sage.libs.coxeter3` module is available in this installation
    of Sage.

    EXAMPLES::

        sage: from sage.features.coxeter3 import Coxeter3
        sage: Coxeter3().is_present()  # needs coxeter3
        FeatureTestResult('coxeter3', True)
        sage: Coxeter3().is_present()  # needs !coxeter3
        FeatureTestResult('coxeter3', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !coxeter3`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.coxeter3 import Coxeter3
        sage: Coxeter3().is_present_at_runtime()  # needs coxeter3
        FeatureTestResult('coxeter3', True)

    """
    _enabled_in_build = coxeter3_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.coxeter3 import Coxeter3
            sage: Coxeter3()
            Feature('coxeter3')

        """
        module_name = "sage.libs.coxeter3.coxeter"
        super().__init__("coxeter3", module_name)


def all_features():
    return [Coxeter3()]
