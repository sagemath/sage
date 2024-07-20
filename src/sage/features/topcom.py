# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of topcom executables
"""

# *****************************************************************************
#       Copyright (C) 2022-2024 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import Executable
from .join_feature import JoinFeature


class TOPCOMExecutable(Executable):
    r"""
    A :class:`~sage.features.Feature` which checks for executables from the :ref:`TOPCOM <spkg_topcom>` package.

    EXAMPLES::

        sage: from sage.features.topcom import TOPCOMExecutable
        sage: TOPCOMExecutable('points2allfinetriangs').is_present()    # optional - topcom
        FeatureTestResult('topcom_points2allfinetriangs', True)
    """
    def __init__(self, name):
        r"""
        TESTS::

            sage: from sage.features.topcom import TOPCOMExecutable
            sage: isinstance(TOPCOMExecutable('points2finetriangs'), TOPCOMExecutable)
            True
        """
        Executable.__init__(self, name=f"topcom_{name}",
                            executable=name,
                            spkg="topcom")


class TOPCOM(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the executables
    which comes as a part of :ref:`TOPCOM <spkg_topcom>`.

    EXAMPLES::

        sage: from sage.features.topcom import TOPCOM
        sage: TOPCOM().is_present()                             # optional - topcom
        FeatureTestResult('topcom', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.topcom import TOPCOM
            sage: isinstance(TOPCOM(), TOPCOM)
            True
        """
        JoinFeature.__init__(self, "topcom",
                             [TOPCOMExecutable(name)
                              for name in ('points2allfinetriangs',)])


def all_features():
    return [TOPCOM()]
