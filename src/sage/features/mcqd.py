# sage_setup: distribution = sagemath-environment
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

from . import PythonModule
from .join_feature import JoinFeature


class Mcqd(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the :mod:`~sage.graphs.mcqd` module,
    which is the SageMath interface to the :ref:`mcqd <spkg_mcqd>` library

    EXAMPLES::

        sage: from sage.features.mcqd import Mcqd
        sage: Mcqd().is_present()  # optional - mcqd
        FeatureTestResult('mcqd', True)
    """

    def __init__(self):
        """
        TESTS::

            sage: from sage.features.mcqd import Mcqd
            sage: isinstance(Mcqd(), Mcqd)
            True
        """
        JoinFeature.__init__(self, 'mcqd',
                             [PythonModule('sage.graphs.mcqd',
                                           spkg='sagemath_mcqd')])


def all_features():
    return [Mcqd()]
