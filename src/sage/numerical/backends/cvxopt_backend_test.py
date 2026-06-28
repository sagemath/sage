import pytest

from sage.structure.sage_object import SageObject
from sage.numerical.backends.generic_backend_test import GenericBackendTests
from sage.numerical.backends.generic_backend import GenericBackend
from sage.numerical.mip import MixedIntegerLinearProgram


class TestCVXOPTBackend(GenericBackendTests):

    @pytest.fixture
    def backend(self) -> GenericBackend:
        return MixedIntegerLinearProgram(solver="CVXOPT").get_backend()

    def test_sage_unittest_testsuite(self, sage_object: SageObject, run_test_suite):
        # TODO: Remove this test as soon as all old test methods are migrated
        run_test_suite(
            sage_object,
            verbose=True,
            skip=("_test_pickling", "_test_solve", "_test_solve_trac_18572"),
        )
