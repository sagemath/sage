# sage_setup: distribution = sagemath-environment
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

from . import PythonModule
from .join_feature import JoinFeature


class Coxeter3(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` which describes whether the :mod:`sage.libs.coxeter3`
    module is available in this installation of Sage.

    EXAMPLES::

        sage: from sage.features.coxeter3 import Coxeter3
        sage: Coxeter3().require()  # optional - coxeter3
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.coxeter3 import Coxeter3
            sage: Coxeter3()
            Feature('coxeter3')
        """
        JoinFeature.__init__(self, "coxeter3",
                             [PythonModule("sage.libs.coxeter3.coxeter",
                                           spkg='sagemath_coxeter3')])


def all_features():
    return [Coxeter3()]
