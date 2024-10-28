# pylint: disable=missing-function-docstring
from _pytest.fixtures import FixtureRequest
import pytest

from sage.manifolds.manifold import Manifold
from sage.manifolds.differentiable.manifold import DifferentiableManifold
from sage.manifolds.differentiable.examples.sphere import Sphere
from sage.manifolds.differentiable.symplectic_form import SymplecticForm
from sage.manifolds.differentiable.examples.symplectic_space import (
    StandardSymplecticSpace,
)
from sage.symbolic.function_factory import function
from sage.manifolds.differentiable.scalarfield import DiffScalarField


class TestGenericSymplecticForm:
    @pytest.fixture
    def omega(self):
        # Generic 6-dimensional manifold
        M = Manifold(6, "M")
        return SymplecticForm(M, "omega")

    def test_repr(self, omega: SymplecticForm):
        assert (
            str(omega)
            == "Symplectic form omega on the 6-dimensional differentiable manifold M"
        )

    def test_new_instance_repr(self, omega: SymplecticForm):
        omega1 = omega._new_instance()  # type: ignore reportPrivateUsage
        assert (
            str(omega1)
            == "Symplectic form unnamed symplectic form on the 6-dimensional differentiable manifold M"
        )

    def test_new_instance_same_type(self, omega: SymplecticForm):
        omega1 = omega._new_instance()  # type: ignore reportPrivateUsage
        assert type(omega1) is type(omega)

    def test_new_instance_same_parent(self, omega: SymplecticForm):
        omega1 = omega._new_instance()  # type: ignore reportPrivateUsage
        assert omega1.parent() == omega.parent()

    def test_init_derived_sets_poisson_tensor(self, omega: SymplecticForm):
        omega._init_derived()  # type: ignore reportPrivateUsage
        assert omega._poisson is None  # type: ignore reportPrivateUsage

    def test_del_derived_resets_poisson_tensor(self, omega: SymplecticForm):
        omega._del_derived()  # type: ignore reportPrivateUsage
        assert omega._poisson is None  # type: ignore reportPrivateUsage


class TestCoherenceOfFormulas:
    r"""
    Test correctness of the implementation, by checking that equivalent formulas give the correct result.
    We check it for the examples of `\R^2` and `S^2`, which should be enough.
    """

    @pytest.fixture(params=['R2', 'S2'])
    def M(self, request: FixtureRequest):
        if request.param == 'R2':
            return StandardSymplecticSpace(2, 'R2', symplectic_name='omega')
        elif request.param == 'S2':
            # Init stereographic coordinates to get a complete atlas
            return Sphere(2, coordinates='stereographic')

    @pytest.fixture()
    def omega(self, M):
        if isinstance(M, StandardSymplecticSpace):
            return M.symplectic_form()
        else:
            return SymplecticForm.wrap(M.metric().volume_form(), "omega")

    def test_flat_of_hamiltonian_vector_field(
        self, M: DifferentiableManifold, omega: SymplecticForm
    ):
        H = generic_scalar_field(M, "H")
        assert omega.flat(omega.hamiltonian_vector_field(H)) == -H.differential()

    def test_hamiltonian_vector_field_contr_symplectic_form(
        self, M: DifferentiableManifold, omega: SymplecticForm
    ):
        H = generic_scalar_field(M, "H")
        # X_H \lrcorner \omega + \dif H = 0
        assert omega.contract(0, omega.hamiltonian_vector_field(H)) == -H.differential()

    def test_poisson_bracket_as_contraction_symplectic_form(
        self, M: DifferentiableManifold, omega: SymplecticForm
    ):
        f = generic_scalar_field(M, "f")
        g = generic_scalar_field(M, "g")
        # {f, g} = \omega(X_f, X_g)
        assert omega.poisson_bracket(f, g) == omega.contract(
            0, omega.hamiltonian_vector_field(f)
        ).contract(0, omega.hamiltonian_vector_field(g))

    def test_poisson_bracket_as_contraction_poisson_tensor(
        self, M: DifferentiableManifold, omega: SymplecticForm
    ):
        f = generic_scalar_field(M, "f")
        g = generic_scalar_field(M, "g")
        # {f, g} = \pi(\dif f, \dif g)
        assert omega.poisson_bracket(f, g) == omega.poisson().contract(
            0, f.exterior_derivative()
        ).contract(0, g.exterior_derivative())

    def test_poisson_bracket_as_contraction_hamiltonian_vector_field(
        self, M: DifferentiableManifold, omega: SymplecticForm
    ):
        f = generic_scalar_field(M, "f")
        g = generic_scalar_field(M, "g")
        # {f, g} =  X_f (g)
        assert omega.poisson_bracket(f, g) == omega.hamiltonian_vector_field(f)(g)
        # {f, g} =  -X_g(f)
        assert omega.poisson_bracket(f, g) == -omega.hamiltonian_vector_field(g)(f)

    def test_poisson_bracket_as_commutator_hamiltonian_vector_fields(
        self, M: DifferentiableManifold, omega: SymplecticForm
    ):
        f = generic_scalar_field(M, "f")
        g = generic_scalar_field(M, "g")
        # [X_f, X_g] = X_{{f,g}}
        assert omega.hamiltonian_vector_field(f).bracket(
            omega.hamiltonian_vector_field(g)
        ) == omega.hamiltonian_vector_field(omega.poisson_bracket(f, g))

    def test_hodge_star_of_one_is_volume(
        self, M: DifferentiableManifold, omega: SymplecticForm
    ):
        assert M.one_scalar_field().hodge_dual(omega) == omega.volume_form()

    def test_hodge_star_of_volume_is_one(
        self, M: DifferentiableManifold, omega: SymplecticForm
    ):
        assert omega.volume_form().hodge_dual(omega) == M.one_scalar_field()

    def test_trace_of_two_form_is_given_using_contraction_with_omega(
        self, M: DifferentiableManifold, omega: SymplecticForm
    ):
        a = M.diff_form(2)
        a[1,2] = 3
        assert a.trace(using=omega) == a.up(omega, 1).trace()

    def test_omega_on_forms_is_determinant_for_decomposables(
        self, M: DifferentiableManifold, omega: SymplecticForm
    ):
        a = M.one_form(1,2)
        b = M.one_form(3,4)
        c = M.one_form(5,6)
        d = M.one_form(7,8)

        assert omega.on_forms(a.wedge(b), c.wedge(d)) == omega.on_forms(a,c) * omega.on_forms(b, d) - omega.on_forms(a,d) * omega.on_forms(b,c)

    def test_omega_on_one_forms_is_omega_on_dual_vectors(
        self, M: DifferentiableManifold, omega: SymplecticForm
    ):
        a = M.one_form(1,2)
        b = M.one_form(3,4)
        assert omega.on_forms(a, b) == omega(a.up(omega), b.up(omega))


def generic_scalar_field(M: DifferentiableManifold, name: str) -> DiffScalarField:
    chart_functions = {chart: function(name)(*chart[:]) for chart in M.atlas()}
    return M.scalar_field(chart_functions, name=name)


class TestR2VectorSpace:
    @pytest.fixture
    def M(self):
        return StandardSymplecticSpace(2, 'R2', symplectic_name='omega')

    @pytest.fixture
    def omega(self, M):
        return M.symplectic_form()

    def test_display(self, omega: SymplecticForm):
        assert str(omega.display()) == r"omega = -dq∧dp"

    def test_poisson(self, omega: SymplecticForm):
        # TODO: Shouldn't this be written with wedge product?
        assert str(omega.poisson().display()) == r"poisson_omega = -e_q∧e_p"

    def test_hamiltonian_vector_field(
        self, M: StandardSymplecticSpace, omega: SymplecticForm
    ):
        H = generic_scalar_field(M, "H")
        XH = omega.hamiltonian_vector_field(H)
        assert str(XH.display()) == r"XH = d(H)/dp e_q - d(H)/dq e_p"

    def test_flat(self, M: StandardSymplecticSpace, omega: SymplecticForm):
        X = M.vector_field(1, 2, name='X')
        assert str(X.display()) == r"X = e_q + 2 e_p"
        assert str(omega.flat(X).display()) == r"X_flat = 2 dq - dp"

    def test_hodge_star(self, M: StandardSymplecticSpace, omega: SymplecticForm):
        # Standard basis
        e = M.one_form(0,1, name='e')
        f = M.one_form(1,0, name='f')
        assert e.wedge(f) == omega

        assert M.one_scalar_field().hodge_dual(omega) == omega
        assert e.hodge_dual(omega) == e
        assert f.hodge_dual(omega) == f
        assert omega.hodge_dual(omega) == M.one_scalar_field()

    def test_omega_on_one_forms(self, M: StandardSymplecticSpace, omega: SymplecticForm):
        # Standard basis
        e = M.one_form(0,1, name='e')
        f = M.one_form(1,0, name='f')
        assert e.wedge(f) == omega

        assert omega.on_forms(e, f) == 1

    def test_hodge_star_is_given_using_omega_on_forms(
        self, M: StandardSymplecticSpace, omega: SymplecticForm
    ):
        a = M.one_form(1,2)
        b = M.one_form(3,4)
        assert a.wedge(b.hodge_dual(omega)) == omega.on_forms(a, b) * omega.volume_form()
