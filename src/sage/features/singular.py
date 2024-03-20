# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of ``singular`` and the SageMath interfaces to it
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
from .sagemath import sage__libs__singular
from sage.env import SINGULAR_BIN


class Singular(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the :ref:`singular <spkg_singular>` executable.

    .. SEEALSO::

        :class:`Feature sage.libs.singular <~sage.features.sagemath.sage__libs__singular>`

    EXAMPLES::

        sage: from sage.features.singular import Singular
        sage: Singular().is_present()                                                   # needs singular
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
                            spkg='singular', type='standard')


def all_features():
    return [Singular()]
