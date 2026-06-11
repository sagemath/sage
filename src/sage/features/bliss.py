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

from sage.config import bliss_enabled
from sage.features import CythonFeature, PythonModule
from sage.features.build_feature import BuildFeature

TEST_CODE = """
# distutils: language=c++
# distutils: libraries=bliss

cdef extern from "bliss/graph.hh" namespace "bliss":
    cdef cppclass Graph:
        Graph(const unsigned int)

from cysignals.signals cimport sig_on, sig_off

sig_on()
Graph(1)
sig_off()
"""


class BlissLibrary(CythonFeature):
    r"""
    A :class:`~sage.features.Feature` which describes whether the :ref:`Bliss library <spkg_bliss>` is
    present and functional.

    EXAMPLES::

        sage: from sage.features.bliss import BlissLibrary
        sage: BlissLibrary().require()  # optional - libbliss
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.bliss import BlissLibrary
            sage: BlissLibrary()
            Feature('libbliss')
        """
        CythonFeature.__init__(self, "libbliss", test_code=TEST_CODE,
                               spkg='bliss',
                               url='http://www.tcs.hut.fi/Software/bliss/')


class Bliss(BuildFeature):
    r"""
    A :class:`~sage.features.Feature` which describes whether the
    :mod:`sage.graphs.bliss` module is available in this installation
    of Sage.

    EXAMPLES::

        sage: from sage.features.bliss import Bliss
        sage: Bliss().require()  # needs bliss
    """
    _enabled_in_build = bliss_enabled

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.bliss import Bliss
            sage: Bliss()
            Feature('bliss')

        """
        super().__init__("bliss",
                         url='http://www.tcs.hut.fi/Software/bliss/')

    def is_present_at_runtime(self):
        r"""
        TESTS::

            sage: from sage.features import FeatureTestResult
            sage: from sage.features.bliss import Bliss
            sage: result = Bliss().is_present_at_runtime()
            sage: isinstance(result, FeatureTestResult)
            True
            sage: result  # needs bliss
            FeatureTestResult('bliss', True)

        """
        result = PythonModule("sage.graphs.bliss")._is_present()
        result.feature = self
        return result

def all_features():
    return [Bliss()]
