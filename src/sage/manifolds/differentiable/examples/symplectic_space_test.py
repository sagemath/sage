import sage.all
from sage.manifolds.differentiable.symplectic_form import SymplecticForm
from sage.manifolds.differentiable.examples.symplectic_space import (
    StandardSymplecticSpace,
)
import pytest


class TestR2VectorSpace:
    @pytest.fixture
    def M(self):
        return StandardSymplecticSpace(2, 'R2', symplectic_name='omega')

    @pytest.fixture
    def omega(self, M: StandardSymplecticSpace):
        return M.symplectic_form()

    def test_repr(self, M: StandardSymplecticSpace):
        assert str(M) == "Standard symplectic space R2"

    def test_display(self, omega: SymplecticForm):
        assert str(omega.display()) == r"omega = -dq∧dp"


class TestR4VectorSpace:
    @pytest.fixture
    def M(self):
        return StandardSymplecticSpace(4, 'R4', symplectic_name='omega')

    @pytest.fixture
    def omega(self, M: StandardSymplecticSpace):
        return M.symplectic_form()

    def test_repr(self, M: StandardSymplecticSpace):
        assert str(M) == "Standard symplectic space R4"

    def test_display(self, omega: SymplecticForm):
        assert str(omega.display()) == r"omega = -dq1∧dp1 - dq2∧dp2"
