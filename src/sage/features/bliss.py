# -*- coding: utf-8 -*-
r"""
Checks for bliss
"""
from . import CythonFeature, PythonModule
from .join_feature import JoinFeature


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
    A :class:`Feature` which describes whether the Bliss library is
    present and functional.

    EXAMPLES::

        sage: from sage.features.bliss import BlissLibrary
        sage: BlissLibrary().require()  # optional: bliss
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.bliss import BlissLibrary
            sage: BlissLibrary()
            Feature('Bliss')
        """
        CythonFeature.__init__(self, "Bliss", test_code=TEST_CODE,
                               spkg="bliss",
                               url="http://www.tcs.hut.fi/Software/bliss/")


class Bliss(JoinFeature):
    r"""
    A :class:`Feature` which describes whether the :mod:`sage.graphs.bliss`
    module is available in this installation of Sage.

    EXAMPLES::

        sage: from sage.features.bliss import Bliss
        sage: Bliss().require()  # optional: bliss
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.bliss import Bliss
            sage: Bliss()
            Feature('bliss')
        """
        # Currently part of sagemath_standard, conditionally built.
        # Will be changed to spkg='sagemath_bliss' later
        JoinFeature.__init__(self, "bliss",
                             [PythonModule("sage.graphs.bliss", spkg="bliss",
                                           url="http://www.tcs.hut.fi/Software/bliss/")])
