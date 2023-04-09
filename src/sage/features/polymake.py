r"""
Feature for testing the presence of the Python interface to polymake
"""

# *****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import PythonModule
from .join_feature import JoinFeature


class JuPyMake(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the :mod:`JuPyMake`
    module, a Python interface to the polymake library.

    EXAMPLES::

        sage: from sage.features.polymake import JuPyMake
        sage: JuPyMake().is_present()  # optional - jupymake
        FeatureTestResult('jupymake', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.polymake import JuPyMake
            sage: isinstance(JuPyMake(), JuPyMake)
            True
        """
        JoinFeature.__init__(self, "jupymake",
                             [PythonModule("JuPyMake", spkg="jupymake")])


def all_features():
    return [JuPyMake()]
