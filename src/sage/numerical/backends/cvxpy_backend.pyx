# sage_setup: distribution = sagemath-polyhedra
# sage.doctest: optional - cvxpy
r"""
CVXPY Backend

AUTHORS:

- Nathann Cohen (2010-10)      : generic_backend template
- Matthias Koeppe (2022-03)    : this backend
"""
# ****************************************************************************
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#       Copyright (C) 2022 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.numerical.mip import MIPSolverException
from sage.rings.real_double import RDF
from copy import copy
import cvxpy
from cvxpy.atoms.affine.add_expr import AddExpression
from cvxpy.expressions.constants import Constant
from cvxpy.constraints.zero import Equality

cdef class CVXPYBackend:
    """
    MIP Backend that delegates to CVXPY.

    CVXPY interfaces to various solvers, see
    https://www.cvxpy.org/install/index.html#install and
    https://www.cvxpy.org/tutorial/advanced/index.html#choosing-a-solver

    EXAMPLES::

        sage: import cvxpy
        sage: cvxpy.installed_solvers()                                                # random

    Using the default solver determined by CVXPY::

        sage: p = MixedIntegerLinearProgram(solver='CVXPY'); p.solve()
        0.0

    Using a specific solver::

        sage: p = MixedIntegerLinearProgram(solver='CVXPY/OSQP'); p.solve()
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/ECOS'); p.solve()
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/SCS'); p.solve()
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/SciPy/HiGHS'); p.solve()
        0.0

    Open-source solvers provided by optional packages::

        sage: p = MixedIntegerLinearProgram(solver='CVXPY/GLPK'); p.solve()             # needs cvxopt
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/GLPK_MI'); p.solve()          # needs cvxopt
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/CVXOPT'); p.solve()           # needs cvxopt
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/GLOP'); p.solve()            # optional - ortools
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/PDLP'); p.solve()            # optional - ortools
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/CBC'); p.solve()             # optional - cylp
        0.0

    Non-free solvers::

        sage: p = MixedIntegerLinearProgram(solver='CVXPY/Gurobi'); p.solve()          # optional - gurobi
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/CPLEX'); p.solve()           # optional - cplex
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/MOSEK'); p.solve()           # optional - mosek
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/SCIP'); p.solve()            # optional - pyscipopt
        0.0
        sage: p = MixedIntegerLinearProgram(solver='CVXPY/XPRESS'); p.solve()          # optional - xpress
        0.0
        sage: p = MixedIntegerLinearProgram(solver="CVXPY/NAG"); p.solve()             # optional - naginterfaces
        0.0
    """

    def __cinit__(self, maximization=True, base_ring=None, cvxpy_solver=None, cvxpy_solver_args=None):
        """
        Cython constructor.

        INPUT:

        - ``maximization``-- boolean (default: ``True``); whether this is a
          maximization or minimization problem

        - ``base_ring`` -- (optional) must be ``RDF`` if provided

        - ``cvxpy_solver`` -- (optional) passed to :meth:`cvxpy.Problem.solve`
          as the parameter ``solver``

        - ``cvxpy_solver_args`` -- dictionary (optional); passed to
          :meth:`cvxpy.Problem.solve` as additional keyword arguments

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
        """
        if base_ring is None:
            base_ring = RDF
        elif base_ring != RDF:
            raise ValueError('base_ring must be RDF')

        if cvxpy_solver_args is None:
            cvxpy_solver_args = {}

        if isinstance(cvxpy_solver, str):
            cvxpy_solver = cvxpy_solver.upper()
            if cvxpy_solver.startswith("SCIPY/"):
                cvxpy_solver_args['scipy_options'] = {"method": cvxpy_solver[len("SCIPY/"):]}
                cvxpy_solver = "SCIPY"
            cvxpy_solver = getattr(cvxpy, cvxpy_solver)
        self._cvxpy_solver = cvxpy_solver
        self._cvxpy_solver_args = cvxpy_solver_args

        self.set_verbosity(0)
        self.variables = []
        self.constraint_names = []
        self.Matrix = []
        self.row_lower_bound = []
        self.row_upper_bound = []
        self.col_lower_bound = []
        self.col_upper_bound = []
        self.objective_coefficients = []
        self.obj_constant_term = 0.0
        if maximization:
            objective = cvxpy.Maximize(0)
        else:
            objective = cvxpy.Minimize(0)
        self.problem = cvxpy.Problem(objective, ())

    cpdef __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver='CVXPY')
            sage: b = p.new_variable()
            sage: p.add_constraint(b[1] + b[2] <= 6)
            sage: p.set_objective(b[1] + b[2])
            sage: cp = copy(p.get_backend())
            sage: cp.solve()
            doctest:warning...:
            FutureWarning:...
            0
            sage: cp.get_objective_value()  # abs tol 1e-7
            6.0
        """
        cdef CVXPYBackend cp = type(self)(base_ring=self.base_ring())
        cp.problem = self.problem                   # it's considered immutable; so no need to copy.
        cp.variables = copy(self.variables)
        cp.constraint_names = copy(self.constraint_names)
        cp.Matrix = [row[:] for row in self.Matrix]
        cp.row_lower_bound = self.row_lower_bound[:]
        cp.row_upper_bound = self.row_upper_bound[:]
        cp.col_lower_bound = self.col_lower_bound[:]
        cp.col_upper_bound = self.col_upper_bound[:]
        cp.objective_coefficients = self.objective_coefficients[:]
        cp.obj_constant_term = self.obj_constant_term
        return cp

    cpdef cvxpy_problem(self):
        return self.problem

    def cvxpy_variables(self):
        return self.variables

    cpdef int add_variable(self, lower_bound=0, upper_bound=None,
                           binary=False, continuous=True, integer=False,
                           obj=None, name=None, coefficients=None) except -1:
        ## coefficients is an extension in this backend,
        ## and a proposed addition to the interface, to unify this with add_col.
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both nonnegative and real.

        INPUT:

        - ``lower_bound`` -- the lower bound of the variable (default: 0)

        - ``upper_bound`` -- the upper bound of the variable (default: ``None``)

        - ``binary`` -- ``True`` if the variable is binary (default: ``False``)

        - ``continuous`` -- ``True`` if the variable is continuous (default: ``True``)

        - ``integer`` -- ``True`` if the variable is integral (default: ``False``)

        - ``obj`` -- (optional) coefficient of this variable in the objective function (default: 0)

        - ``name`` -- an optional name for the newly added variable (default: ``None``)

        - ``coefficients`` -- (optional) an iterable of pairs ``(i, v)``; in each
          pair, ``i`` is a row index (integer) and ``v`` is a value (element of :meth:`base_ring`)

        OUTPUT: the index of the newly created variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable(continuous=True, integer=True)
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.add_variable(name='x',obj=1)
            1
            sage: p.col_name(1)
            'x'
            sage: p.objective_coefficient(1)
            1.0
        """
        cdef int vtype = int(binary) + int(continuous) + int(integer)
        if vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        for i in range(len(self.Matrix)):
            self.Matrix[i].append(0.0)
        if lower_bound is not None:
            lower_bound = float(lower_bound)
        self.col_lower_bound.append(lower_bound)
        if upper_bound is not None:
            upper_bound = float(upper_bound)
        self.col_upper_bound.append(upper_bound)
        if obj is None:
            obj = 0.0
        else:
            obj = float(obj)
        self.objective_coefficients.append(obj)

        if name is None:
            name = f'x_{self.ncols()}'

        if binary:
            variable = cvxpy.Variable(name=name, boolean=True)
        elif integer:
            variable = cvxpy.Variable(name=name, integer=True)
        else:
            variable = cvxpy.Variable(name=name)

        self.variables.append(variable)
        index = self.ncols() - 1

        if coefficients is not None:
            constraints = list(self.problem.constraints)
            for i, v in coefficients:
                if not isinstance(constraints[i], Equality):
                    raise NotImplementedError('adding coefficients to inequalities is ambiguous '
                                              'because cvxpy rewrites all inequalities as <=')
                constraints[i] = type(constraints[i])(constraints[i].args[0] + float(v) * variable,
                                                      constraints[i].args[1])
            self.problem = cvxpy.Problem(self.problem.objective, constraints)

        self.add_linear_constraint([(index, 1)], lower_bound, upper_bound)

        if obj:
            objective = type(self.problem.objective)(self.problem.objective.args[0]
                                                     + obj * variable)
            self.problem = cvxpy.Problem(objective, self.problem.constraints)

        return index

    cpdef set_verbosity(self, int level):
        """
        Set the log (verbosity) level.

        This is currently ignored.

        INPUT:

        - ``level`` -- integer; from 0 (no verbosity) to 3

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.set_verbosity(2)
        """
        pass

    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` -- an iterable of pairs ``(i, v)``. In each
          pair, ``i`` is a variable index (integer) and ``v`` is a
          value (element of :meth:`base_ring`).

        - ``lower_bound`` -- element of :meth:`base_ring` or
          ``None``; the lower bound

        - ``upper_bound`` -- element of :meth:`base_ring` or
          ``None``; the upper bound

        - ``name`` -- string or ``None``; optional name for this row

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_variables(5)
            4
            sage: index = p.nrows()
            sage: p.add_linear_constraint( zip(range(5), range(5)), 2, 2)
            sage: p.row(index)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(index)
            (2, 2)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1, 1, name='foo')
            sage: p.row_name(1)
            'constraint_1'
        """
        last = len(self.Matrix)
        self.Matrix.append([])
        for i in range(len(self.objective_coefficients)):
            self.Matrix[last].append(0.0)
        for a in coefficients:
            self.Matrix[last][a[0]] = a[1]

        self.row_lower_bound.append(lower_bound)
        self.row_upper_bound.append(upper_bound)

        terms = [v * self.variables[i] for i, v in coefficients]
        if terms:
            expr = AddExpression(terms)
        else:
            expr = Constant(0)
        constraints = list(self.problem.constraints)
        if lower_bound is not None and lower_bound == upper_bound:
            constraints.append(expr == upper_bound)
        elif lower_bound is not None:
            constraints.append(lower_bound <= expr)
        elif upper_bound is not None:
            constraints.append(expr <= upper_bound)
        self.constraint_names.append(name)
        self.problem = cvxpy.Problem(self.problem.objective, constraints)

    cpdef add_col(self, indices, coeffs):
        """
        Add a column.

        INPUT:

        - ``indices`` -- list of integers; this list contains the
          indices of the constraints in which the variable's
          coefficient is nonzero

        - ``coeffs`` -- list of real values; associates a coefficient
          to the variable in each of the constraints in which it
          appears. Namely, the i-th entry of ``coeffs`` corresponds to
          the coefficient of the variable in the constraint
          represented by the i-th entry in ``indices``.

        .. NOTE::

            ``indices`` and ``coeffs`` are expected to be of the same
            length.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(list(range(5)), list(range(5)))
            Traceback (most recent call last):
            ...
            NotImplementedError: adding coefficients to inequalities is ambiguous because cvxpy rewrites all inequalities as <=
            sage: p.nrows()
            5
        """
        for i in range(len(self.Matrix)):
            self.Matrix[i].append(0)
        for i in range(len(indices)):
            self.Matrix[indices[i]][len(self.Matrix[indices[i]]) - 1] = coeffs[i]

        self.col_lower_bound.append(None)
        self.col_upper_bound.append(None)
        #self.objective_coefficients.append(0) goes to "self.add_variable"
        self.add_variable(coefficients=zip(indices, coeffs))

    cpdef set_objective(self, list coeff, d=0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- list of real values, whose i-th element is the
          coefficient of the i-th variable in the objective function

        - ``d`` -- double; the constant term in the linear function (set to `0`
          by default)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1.0, 1.0, 2.0, 1.0, 3.0]
        """
        if self.variables:
            expr = AddExpression([c * x for c, x in zip(coeff, self.variables)])
            for i in range(len(coeff)):
                self.objective_coefficients[i] = float(coeff[i])
        else:
            expr = Constant(0)
        objective = type(self.problem.objective)(expr)
        constraints = list(self.problem.constraints)
        self.problem = cvxpy.Problem(objective, constraints)
        self.obj_constant_term = d

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` -- integer:

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        expr = self.problem.objective.args[0]
        if sense == 1:
            objective = cvxpy.Maximize(expr)
        else:
            objective = cvxpy.Minimize(expr)
        self.problem = cvxpy.Problem(objective, self.problem.constraints)

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective function.

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``coeff`` -- double; its coefficient or ``None`` for
          reading (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0.0
            sage: p.objective_coefficient(0, 2)
            sage: p.objective_coefficient(0)
            2
        """
        if coeff is not None:
            self.objective_coefficients[variable] = coeff
            expr = self.problem.objective.args[0] + coeff * self.variables[variable]
            objective = type(self.problem.objective)(expr)
            constraints = list(self.problem.constraints)
            self.problem = cvxpy.Problem(objective, constraints)
        else:
            if variable < len(self.objective_coefficients):
                return self.objective_coefficients[variable]
            else:
                return 0

    cpdef int solve(self) except -1:
        """
        Solve the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution cannot be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(list(range(5)), list(range(5)))
            Traceback (most recent call last):
            ...
            NotImplementedError: adding coefficients to inequalities is ambiguous because cvxpy rewrites all inequalities as <=
            sage: p.solve()
            0
            sage: p.objective_coefficient(0,1)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: ...
        """
        try:
            self.problem.solve(solver=self._cvxpy_solver, **self._cvxpy_solver_args)
        except Exception as e:
            raise MIPSolverException(f"cvxpy.Problem.solve raised exception: {e}")
        status = self.problem.status
        if 'optimal' in status:
            return 0
        if 'infeasible' in status:
            raise MIPSolverException(f"cvxpy.Problem.solve: Problem has no feasible solution")
        if 'unbounded' in status:
            raise MIPSolverException(f"cvxpy.Problem.solve: Problem is unbounded")
        raise MIPSolverException(f"cvxpy.Problem.solve reported an unknown problem status: {status}")

    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behavior is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0,1), (1,2)], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()  # abs tol 1e-7
            7.5
            sage: p.get_variable_value(0)  # abs tol 1e-7
            0.0
            sage: p.get_variable_value(1)  # abs tol 1e-7
            1.5
        """
        return self.problem.value + self.obj_constant_term

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behavior is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0,1), (1, 2)], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()  # abs tol 1e-7
            7.5
            sage: p.get_variable_value(0)  # abs tol 1e-7
            0.0
            sage: p.get_variable_value(1)  # abs tol 1e-7
            1.5
        """
        return float(self.variables[variable].value)

    cpdef int ncols(self) noexcept:
        """
        Return the number of columns/variables.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.ncols()
            0
            sage: p.add_variables(2)
            1
            sage: p.ncols()
            2
        """
        return len(self.variables)

    cpdef int nrows(self) noexcept:
        """
        Return the number of rows/constraints.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(2, 0, None)
            sage: p.nrows()
            2
        """
        return len(self.problem.constraints)

    cpdef bint is_maximization(self) noexcept:
        """
        Test whether the problem is a maximization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        return isinstance(self.problem.objective, cvxpy.Maximize)

    cpdef problem_name(self, name=None):
        """
        Return or define the problem's name.

        INPUT:

        - ``name`` -- string; the problem's name. When set to
          ``None`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.problem_name("There_once_was_a_french_fry")
            sage: print(p.problem_name())
            There_once_was_a_french_fry
        """
        if name is None:
            if self.prob_name is not None:
                return self.prob_name
            else:
                return ""
        else:
            self.prob_name = str(name)

    cpdef row(self, int i):
        """
        Return a row.

        INPUT:

        - ``index`` -- integer; the constraint's id

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_linear_constraint`` method.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_variables(5)
            4
            sage: index = p.nrows()
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row(index)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(index)
            (2, 2)
        """
        idx = []
        coef = []
        for j in range(len(self.Matrix[i])):
            if self.Matrix[i][j] != 0:
                idx.append(j)
                coef.append(self.Matrix[i][j])
        return (idx, coef)

    cpdef row_bounds(self, int index):
        """
        Return the bounds of a specific constraint.

        INPUT:

        - ``index`` -- integer; the constraint's id

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the constraint is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_variables(5)
            4
            sage: index = p.nrows()
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row(index)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(index)
            (2, 2)
        """
        return (self.row_lower_bound[index], self.row_upper_bound[index])

    cpdef col_bounds(self, int index):
        """
        Return the bounds of a specific variable.

        INPUT:

        - ``index`` -- integer; the variable's id

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the variable is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0.0, 5)
        """
        return (self.col_lower_bound[index], self.col_upper_bound[index])

    cpdef bint is_variable_binary(self, int index) noexcept:
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_binary(0)
            False
        """
        return bool(self.variables[index].boolean_idx)

    cpdef bint is_variable_integer(self, int index) noexcept:
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_integer(0)
            False
        """
        return bool(self.variables[index].integer_idx)

    cpdef bint is_variable_continuous(self, int index) noexcept:
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_continuous(0)
            True
        """
        return not self.is_variable_binary(index) and not self.is_variable_integer(index)

    cpdef row_name(self, int index):
        """
        Return the ``index``-th row name.

        INPUT:

        - ``index`` -- integer; the row's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_linear_constraint([], 2, 2)
            sage: p.row_name(0)
            'constraint_0'
        """
        return self.constraint_names[index] or ("constraint_" + repr(index))

    cpdef col_name(self, int index):
        """
        Return the ``index``-th col name.

        INPUT:

        - ``index`` -- integer; the col's id

        - ``name`` -- (``char *``) its name; when set to ``NULL``
          (default), the method returns the current name

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_variable()
            0
            sage: p.col_name(0)
            'x_0'
        """
        return self.variables[index].name()

    cpdef variable_upper_bound(self, int index, value=False):
        """
        Return or define the upper bound on a variable.

        INPUT:

        - ``index`` -- integer; the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``None``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0.0, 5)
            sage: p.variable_upper_bound(0, None)
            sage: p.col_bounds(0)
            (0.0, None)
        """
        if value is not False:
            self.col_upper_bound[index] = value
        else:
            return self.col_upper_bound[index]

    cpdef variable_lower_bound(self, int index, value=False):
        """
        Return or define the lower bound on a variable.

        INPUT:

        - ``index`` -- integer; the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``None``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='CVXPY')
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)
            sage: p.col_bounds(0)
            (5, None)
            sage: p.variable_lower_bound(0, None)
            sage: p.col_bounds(0)
            (None, None)
        """
        if value is not False:
            self.col_lower_bound[index] = value
        else:
            return self.col_lower_bound[index]
