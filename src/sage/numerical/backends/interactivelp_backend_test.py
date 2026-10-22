from __future__ import annotations
import pytest
from sage.numerical.backends.generic_backend_test import GenericBackendTests
from sage.numerical.mip import MixedIntegerLinearProgram
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sage.numerical.backends.generic_backend import GenericBackend


class TestInteractiveLPBackend(GenericBackendTests):

    @pytest.fixture
    def backend(self) -> GenericBackend:
        return MixedIntegerLinearProgram(solver='InteractiveLP').get_backend()
