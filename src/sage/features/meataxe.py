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
from sage.features.build_feature import BuildModule

class Meataxe(BuildModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    the Sage modules that depend on the :ref:`meataxe <spkg_meataxe>`
    library.

    EXAMPLES::

        sage: from sage.features.meataxe import Meataxe
        sage: Meataxe().is_present()  # needs meataxe
        FeatureTestResult('meataxe', True)
        sage: Meataxe().is_present()  # needs !meataxe
        FeatureTestResult('meataxe', False)

    A runtime check. We only check the "present" case because, if
    feature checks are _not_ deferred, the ``needs !meataxe`` can be
    satisfied (disabled at build time) at the same time we are able to
    import a module that was installed for a previous build of sage::

        sage: from sage.features.meataxe import Meataxe
        sage: Meataxe().is_present_at_runtime()  # needs meataxe
        FeatureTestResult('meataxe', True)

    """
    _enabled_in_build = meataxe_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.meataxe import Meataxe
            sage: Meataxe()
            Feature('meataxe')

        """
        module_name = "sage.matrix.matrix_gfpn_dense"
        super().__init__("meataxe", module_name)


def all_features():
    return [Meataxe()]
