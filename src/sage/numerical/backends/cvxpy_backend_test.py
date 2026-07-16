import pytest
from sage.numerical.backends.generic_backend_test import GenericBackendTests
from sage.numerical.backends.generic_backend import GenericBackend
from sage.numerical.mip import MixedIntegerLinearProgram

# Skip the whole module if cvxpy is not installed. This must be a
# module-level statement: used as a class decorator, importorskip returns
# the imported module and the decoration fails with "module object is not
# callable" whenever cvxpy *is* present.
pytest.importorskip("cvxpy")


class TestCVXPYBackend(GenericBackendTests):

    @pytest.fixture
    def backend(self) -> GenericBackend:
        return MixedIntegerLinearProgram(solver="CVXPY").get_backend()
