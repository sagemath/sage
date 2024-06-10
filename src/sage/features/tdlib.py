# sage_setup: distribution = sagemath-environment
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

from . import PythonModule
from .join_feature import JoinFeature


class Tdlib(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the SageMath interface to the :ref:`tdlib <spkg_tdlib>` library.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.tdlib import Tdlib
            sage: isinstance(Tdlib(), Tdlib)
            True
        """
        JoinFeature.__init__(self, 'tdlib',
                             [PythonModule('sage.graphs.graph_decompositions.tdlib',
                                           spkg='sagemath_tdlib')])


def all_features():
    return [Tdlib()]
