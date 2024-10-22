# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of ``bliss``
"""

# *****************************************************************************
#       Copyright (C) 2016 Julian Rüth
#                     2018 Jeroen Demeyer
#                     2021 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from . import PythonModule


class Bliss(PythonModule):
    r"""
    A :class:`~sage.features.Feature` which describes whether the :mod:`sage.graphs.bliss`
    module is available in this installation of Sage.

    EXAMPLES::

        sage: from sage.features.bliss import Bliss
        sage: Bliss().require()  # optional - bliss
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.bliss import Bliss
            sage: Bliss()
            Feature('sage.graphs.bliss')
        """
        PythonModule.__init__(
            self,
            "sage.graphs.bliss",
            spkg="sagemath_bliss",
            url="http://www.tcs.hut.fi/Software/bliss/",
        )


def all_features():
    return [Bliss()]
