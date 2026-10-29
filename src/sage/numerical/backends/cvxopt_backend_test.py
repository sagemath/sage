from __future__ import annotations
import pytest

from sage.numerical.backends.generic_backend_test import GenericBackendTests
from sage.numerical.mip import MixedIntegerLinearProgram
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sage.numerical.backends.generic_backend import GenericBackend
    from sage.structure.sage_object import SageObject


class TestCVXOPTBackend(GenericBackendTests):

    @pytest.fixture
    def backend(self) -> GenericBackend:
        return MixedIntegerLinearProgram(solver="CVXOPT").get_backend()

    def test_sage_unittest_testsuite(self, sage_object: SageObject):
        # TODO: Remove this test as soon as all old test methods are migrated
        from sage.misc.sage_unittest import TestSuite
        TestSuite(sage_object).run(verbose=True, raise_on_failure=True, skip=("_test_pickling","_test_solve","_test_solve_trac_18572"))
