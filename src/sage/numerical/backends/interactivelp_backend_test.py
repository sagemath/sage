import pytest
from sage.numerical.backends.generic_backend_test import GenericBackendTests
from sage.numerical.backends.generic_backend import GenericBackend
from sage.numerical.mip import MixedIntegerLinearProgram

class TestInteractiveLPBackend(GenericBackendTests):

    @pytest.fixture
    def backend(self) -> GenericBackend:
        return MixedIntegerLinearProgram(solver='InteractiveLP').get_backend()
