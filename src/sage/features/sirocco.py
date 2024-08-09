# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of ``sirocco``
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


class Sirocco(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` which describes whether the :mod:`sage.libs.sirocco`
    module is available in this installation of Sage.

    EXAMPLES::

        sage: from sage.features.sirocco import Sirocco
        sage: Sirocco().require()  # optional - sirocco
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sirocco import Sirocco
            sage: Sirocco()
            Feature('sirocco')
        """
        JoinFeature.__init__(self, "sirocco",
                             [PythonModule("sage.libs.sirocco",
                                           spkg='sagemath_sirocco')])


def all_features():
    return [Sirocco()]
