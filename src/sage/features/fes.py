# -*- coding: utf-8 -*-
r"""
Features for testing the presence of ``fes``
"""

from . import CythonFeature, PythonModule
from .join_feature import JoinFeature


TEST_CODE = """
# distutils: libraries=fes

from libc.stdint cimport uint64_t
cdef extern from "<fes_interface.h>":
    ctypedef int (*solution_callback_t)(void *, uint64_t)
    void exhaustive_search_wrapper(int n, int n_eqs, int degree, int ***coeffs, solution_callback_t callback, void* callback_state, int verbose)

solutions = 0

class InternalState:
    verbose = False
    sols = []
    max_sols = 0

cdef int report_solution(void *_state, uint64_t i):
    global solutions
    solutions += 1
    return 0

sig_on()
cdef int ***coeffs = <int ***> sig_calloc(1, sizeof(int **))
coeffs[0] = <int **> sig_calloc(3, sizeof(int *))
coeffs[0][0] = <int *> sig_calloc(1, sizeof(int))
coeffs[0][1] = <int *> sig_calloc(2, sizeof(int))
coeffs[0][2] = <int *> sig_calloc(1, sizeof(int))
coeffs[0][2][0] = 1  # x*y = 0
internal_state = InternalState()

exhaustive_search_wrapper(2, 1, 2, coeffs, report_solution, <void *>internal_state, 0)

sig_free(coeffs[0][2])
sig_free(coeffs[0][1])
sig_free(coeffs[0][0])
sig_free(coeffs[0])
sig_free(coeffs)
sig_off()

if solutions != 3:
    raise AssertionError("libFES did not find three solutions for x*y = 0")
"""


class LibFESLibrary(CythonFeature):
    r"""
    A :class:`~sage.features.Feature` which describes whether the FES library
    is present and functional.

    EXAMPLES::

        sage: from sage.features.fes import LibFESLibrary
        sage: LibFESLibrary().require()  # optional - fes
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.fes import LibFESLibrary
            sage: isinstance(LibFESLibrary(), LibFESLibrary)
            True
        """
        CythonFeature.__init__(self, "LibFES", test_code=TEST_CODE, spkg="fes",
                               url="http://www.lifl.fr/~bouillag/fes/")


class LibFES(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` which describes whether the :mod:`sage.libs.fes`
    module has been enabled for this build of Sage and is functional.

    EXAMPLES::

        sage: from sage.features.fes import LibFES
        sage: LibFES().require()  # optional - fes
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.fes import LibFES
            sage: isinstance(LibFES(), LibFES)
            True
        """
        JoinFeature.__init__(self, 'fes',
                             [PythonModule("sage.libs.fes")],
                             spkg="fes",
                             url="http://www.lifl.fr/~bouillag/fes/")
