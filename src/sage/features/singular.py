r"""
Features for testing the presence of Singular
"""

# *****************************************************************************
#       Copyright (C) 2022-2023 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import Executable, PythonModule
from .join_feature import JoinFeature
from sage.env import SINGULAR_BIN


class Singular(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the Singular executable.

    EXAMPLES::

        sage: from sage.features.singular import Singular
        sage: Singular().is_present()
        FeatureTestResult('singular', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.singular import Singular
            sage: isinstance(Singular(), Singular)
            True
        """
        Executable.__init__(self, "singular", SINGULAR_BIN,
                            spkg='singular')


class sage__libs__singular(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of :mod:`sage.libs.singular`
    (the library interface to Singular) and :mod:`sage.interfaces.singular` (the pexpect
    interface to Singular). By design, we do not distinguish between these two, in order
    to facilitate the conversion of code from the pexpect interface to the library
    interface.

    EXAMPLES::

        sage: from sage.features.singular import sage__libs__singular
        sage: sage__libs__singular().is_present()                       # optional - sage.libs.singular
        FeatureTestResult('sage.libs.singular', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.singular import sage__libs__singular
            sage: isinstance(sage__libs__singular(), sage__libs__singular)
            True
        """
        JoinFeature.__init__(self, 'sage.libs.singular',
                             [PythonModule('sage.libs.singular.singular'),
                              PythonModule('sage.interfaces.singular')])


def all_features():
    return [Singular(),
            sage__libs__singular()]
