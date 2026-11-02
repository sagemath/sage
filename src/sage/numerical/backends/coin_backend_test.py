from __future__ import annotations

import pytest

from sage.numerical.backends.generic_backend_test import GenericBackendTests

CoinBackend = pytest.importorskip('sage.numerical.backends.coin_backend').CoinBackend


class TestCoinBackend(GenericBackendTests):

    @pytest.fixture
    def backend(self):
        return CoinBackend()

    def test_sage_unittest_testsuite(self, sage_object):
        from sage.misc.sage_unittest import TestSuite

        # Keep the exclusions used by CoinBackend's existing doctest.
        TestSuite(sage_object).run(raise_on_failure=True,
                                  skip=['_test_pickling', '_test_solve'])

    def test_solver_selection(self):
        from sage.numerical.mip import MixedIntegerLinearProgram

        program = MixedIntegerLinearProgram(solver='Coin')
        assert isinstance(program.get_backend(), CoinBackend)

    @pytest.mark.parametrize('integer', [False, True])
    def test_constraint_generation_with_default_solver(self, integer):
        from sage.numerical.backends.generic_backend import default_mip_solver
        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

        former_solver = default_mip_solver()
        try:
            default_mip_solver('Coin')
            p = MixedIntegerLinearProgram(constraint_generation=True)
            assert isinstance(p.get_backend(), CoinBackend)
            x = p.new_variable(integer=integer, nonnegative=True)
            p.set_objective(x[0] + 2*x[1])
            p.add_constraint(x[0], max=2)
            p.add_constraint(x[1], max=2)
            assert p.solve() == pytest.approx(6)

            # Each solve must include constraints added since the previous solve.
            p.add_constraint(x[0] + x[1], max=2.5)
            assert p.solve() == pytest.approx(4 if integer else 4.5)
            p.add_constraint(x[1], max=0)
            assert p.solve() == pytest.approx(2)

            # Infeasibility and later recovery must also use the current model.
            p.add_constraint(x[0], min=3)
            with pytest.raises(MIPSolverException, match='proven infeasible'):
                p.solve()
            p.remove_constraint(p.number_of_constraints()-1)
            assert p.solve() == pytest.approx(2)
            assert default_mip_solver() == 'Coin'
        finally:
            default_mip_solver(former_solver)

    def test_feature(self):
        from sage.features.mip_backends import COIN

        assert COIN().name == 'cbc'
        assert COIN().spkg == 'cbc'
        assert COIN().is_present()

    def test_integer_infeasible_with_feasible_relaxation(self, backend):
        from sage.numerical.mip import MIPSolverException

        backend.add_variables(4, binary=True)
        # The LP relaxation is feasible, but a sum of binary variables
        # cannot equal 11/7. Checking the LP solver's status used to accept
        # its fractional solution and round it to an invalid assignment.
        backend.add_linear_constraint([(i, 7) for i in range(4)], 11, 11)
        with pytest.raises(MIPSolverException, match='proven infeasible'):
            backend.solve()

        # Solving again must use the new model and its integer solution.
        backend.remove_constraint(0)
        backend.add_linear_constraint([(i, 7) for i in range(4)], 14, 14)
        backend.solve()
        values = [backend.get_variable_value(i) for i in range(4)]
        assert all(v in (0, 1) for v in values)
        assert sum(values) == 2

    @pytest.mark.parametrize('maximization', [False, True])
    def test_integer_solution_and_objective(self, maximization):
        backend = CoinBackend(maximization=maximization)
        backend.add_variables(3, binary=True)
        backend.add_linear_constraint([(0, 3), (1, 4), (2, 5)], None, 7)
        sign = 1 if maximization else -1
        backend.set_objective([sign * c for c in (5, 6, 7)], 13)
        backend.solve()
        assert [backend.get_variable_value(i) for i in range(3)] == [1, 1, 0]
        assert backend.get_objective_value() == pytest.approx(13 + sign * 11)

    @pytest.mark.parametrize('integer', [False, True])
    def test_unbounded(self, backend, integer):
        from sage.numerical.mip import MIPSolverException

        backend.add_variable(integer=integer, obj=1)
        with pytest.raises(MIPSolverException, match='unbounded'):
            backend.solve()

    def test_basis_solution_after_solve(self, backend):
        backend.add_variables(2)
        backend.add_linear_constraint([(0, 2), (1, -3)], None, 6)
        backend.add_linear_constraint([(0, 3), (1, 2)], None, 6)
        backend.set_objective([1, 1], 7)
        backend.solve()
        assert backend.get_objective_value() == pytest.approx(10)

        # An explicitly selected basis need not be optimal. Its values
        # must take precedence over the previous branch-and-bound result.
        assert backend.set_basis_status([3, 3], [1, 1]) == 0
        assert [backend.get_variable_value(i) for i in range(2)] == [0, 0]
        assert backend.get_objective_value() == pytest.approx(7)
        assert backend.set_basis_status([1, 3], [1, 3]) == 0
        assert [backend.get_variable_value(i) for i in range(2)] == [2, 0]
        assert backend.get_objective_value() == pytest.approx(9)
        backend.solve()
        assert backend.get_objective_value() == pytest.approx(10)
