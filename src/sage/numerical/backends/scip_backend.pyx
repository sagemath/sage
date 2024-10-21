# sage_setup: distribution = sagemath-polyhedra
# distutils: language = c++
# sage.doctest: optional - pyscipopt
"""
SCIP Backend

AUTHORS:

- Nathann Cohen (2010-10): generic backend
- Matthias Koeppe (2017): stubs
- Moritz Firsching (2018-04): rest
"""
# *****************************************************************************
#       Copyright (C) 2010      Nathann Cohen <nathann.cohen@gmail.com>
#                     2017-2022 Matthias Koeppe
#                     2018      Moritz Firsching
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from os.path import splitext
from sage.numerical.mip import MIPSolverException
from pyscipopt import Model


cdef class SCIPBackend(GenericBackend):
    """
    MIP Backend that uses the SCIP solver.

    TESTS:

    General backend testsuite::

        sage: p = MixedIntegerLinearProgram(solver='SCIP')
        sage: TestSuite(p.get_backend()).run(skip='_test_pickling')
    """

    def __cinit__(self, maximization=True):
        """
        Constructor.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='SCIP')
        """
        self.model = Model('')
        if maximization:
            self.set_sense(1)
        else:
            self.set_sense(-1)
        self.obj_constant_term = 0.0
        self.variables = []
        self.model.hideOutput()
        # always set this to None if the list of constraints may change
        self.constraints = None

    def get_constraints(self):
        """
        Get all constraints of the problem.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver='SCIP')
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(zip([0, 1, 2], [8, 6, 1]), None, 48)
            sage: lp.add_linear_constraint(zip([0, 1, 2], [2, 1.5, 0.5]), None, 8)

            sage: lp.get_constraints()
            [c1, c2]
            sage: lp.row(1)                                                     # indirect doctest
            ([0, 1, 2], [2.0, 1.5, 0.5])
        """
        if self.constraints is None:
            self.constraints = self.model.getConss()
        return self.constraints

    cpdef _get_model(self):
        """
        Get the model as a pyscipopt Model.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p._get_model()
            <pyscipopt.scip.Model object at ...
        """
        return self.model

    cpdef int add_variable(self, lower_bound=0.0, upper_bound=None,
                           binary=False, continuous=False, integer=False,
                           obj=0.0, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive, real and the coefficient in the
        objective function is 0.0.

        INPUT:

        - ``lower_bound`` -- the lower bound of the variable (default: 0)

        - ``upper_bound`` -- the upper bound of the variable (default: ``None``)

        - ``binary`` -- ``True`` if the variable is binary (default: ``False``)

        - ``continuous`` -- ``True`` if the variable is binary (default: ``True``)

        - ``integer`` -- ``True`` if the variable is binary (default: ``False``)

        - ``obj`` -- (optional) coefficient of this variable in the objective function (default: 0.0)

        - ``name`` -- an optional name for the newly added variable (default: ``None``)

        OUTPUT: the index of the newly created variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable(binary=True)
            1
            sage: p.add_variable(lower_bound=-2.0, integer=True)
            2
            sage: p.add_variable(continuous=True, integer=True)
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.add_variable(name='x', obj=1.0)
            3
            sage: p.col_name(3)
            'x'
            sage: p.objective_coefficient(3)
            1.0
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer'"
                             "and 'continuous' must be 'True'.")

        if name is None:
            vname = ''
        else:
            assert type(name) in [str, unicode]
            vname = name

        if continuous:
            vtypestr = 'C'
        if binary:
            vtypestr = 'B'
        if integer:
            vtypestr = 'I'

        v = self.model.addVar(name=vname, vtype=vtypestr, ub=upper_bound,
                              lb=lower_bound, obj=obj, pricedVar=False)
        index = v.getIndex()
        assert index == self.ncols()
        self.variables.append(v)

        return index

    cpdef set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable.

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``vtype`` -- integer:

            *  1  Integer
            *  0  Binary
            * -1  Real

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)
            True
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        vtypenames = {1: 'I', 0: 'B', -1: 'C'}
        self.model.chgVarType(var=self.variables[variable], vtype=vtypenames[vtype])

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` -- integer:

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        if sense == 1:
            self.model.setMaximize()
        elif sense == -1:
            self.model.setMinimize()
        else:
            raise AssertionError("sense must be either 1 or -1")

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective function.

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``coeff`` -- double; its coefficient or ``None`` for
          reading (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0.0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2.0
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        if coeff is None:
            return self.variables[variable].getObj()
        else:
            objexpr = self.model.getObjective()
            var = self.variables[variable]
            linfun = sum([e * c for e, c in objexpr.terms.iteritems() if e != var]) + var * coeff
            self.model.setObjective(linfun, sense=self.model.getObjectiveSense())

    cpdef problem_name(self, name=None):
        """
        Return or define the problem's name.

        INPUT:

        - ``name`` -- string; the problem's name. When set to
          ``None`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.problem_name("Nomen est omen")
            sage: p.problem_name()
            'Nomen est omen'
        """
        if name is None:
            return self.model.getProbName()
        else:
            self.model.setProbName(name)

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
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1.0, 1.0, 2.0, 1.0, 3.0]
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        linfun = sum([c * x for c, x in zip(coeff, self.variables)]) + d
        self.model.setObjective(linfun, sense=self.model.getObjectiveSense())

    cpdef set_verbosity(self, int level):
        """
        Set the verbosity level.

        INPUT:

        - ``level`` -- integer; from 0 (no verbosity) to 1

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.set_verbosity(1)

        TODOs::

            - Currently, the output is written to stdout, even when running
            Jupyter: https://github.com/SCIP-Interfaces/PySCIPOpt/issues/116 .
            This should be fixed upstream
            - get access to more verbosity levels (e.g. via parameter settings)
        """
        if level == 0:
            self.model.hideOutput()
        elif level == 1:
            self.model.hideOutput(False)
        else:
            raise AssertionError('level must be "0" or "1"')

    cpdef remove_constraint(self, int i):
        r"""
        Remove a constraint from ``self``.

        INPUT:

        - ``i`` -- index of the constraint to remove

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='SCIP')
            sage: x, y = p['x'], p['y']
            sage: p.add_constraint(2*x + 3*y <= 6)
            sage: p.add_constraint(3*x + 2*y <= 6)
            sage: p.add_constraint(x >= 0)
            sage: p.set_objective(x + y + 7)
            sage: p.set_integer(x); p.set_integer(y)
            sage: p.solve()
            9.0
            sage: p.remove_constraint(0)
            sage: p.solve()
            10.0

        Removing fancy constraints does not make Sage crash::

            sage: MixedIntegerLinearProgram(solver='SCIP').remove_constraint(-2)
            Traceback (most recent call last):
            ...
            ValueError: The constraint's index i must satisfy 0 <= i < number_of_constraints
        """
        if i < 0 or i >= self.nrows():
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
            self.constraints = None
        self.model.delCons(self.get_constraints()[i])
        self.constraints = None

    cpdef remove_constraints(self, constraints):
        r"""
        Remove several constraints.

        INPUT:

        - ``constraints`` -- an iterable containing the indices of the rows to remove

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='SCIP')
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0, 2), (1, 3)], None, 6)
            sage: p.add_linear_constraint([(0, 3), (1, 2)], None, 6)
            sage: p.row(0)
            ([0, 1], [2.0, 3.0])
            sage: p.remove_constraints([0, 1])
            sage: p.nrows()
            0
        """
        if isinstance(constraints, int):
            self.remove_constraint(constraints)
            return

        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
            self.constraints = None

        all_constraints = self.get_constraints()
        to_remove = [all_constraints[i] for i in constraints]
        for constraint in to_remove:
            self.model.delCons(constraint)
        self.constraints = None

    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` an iterable with ``(c, v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (real
          value).

        - ``lower_bound`` -- a lower bound, either a real value or ``None``

        - ``upper_bound`` -- an upper bound, either a real value or ``None``

        - ``name`` -- an optional name for this row (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint( zip(range(5), range(5)), 2.0, 2.0)
            sage: p.row_bounds(0)
            (2.0, 2.0)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1.0, 1.0, name='foo')
            sage: p.row_name(1)
            'foo'
        """
        if self.model.getStatus() != 'unknown':
            self.model.freeTransform()
        mvars = self.variables
        from pyscipopt.scip import quicksum
        linfun = quicksum([v * mvars[c] for c, v in coefficients])
        # we introduced patch 0001 for pyscipopt, in order to handle the case
        # when linfun is an empty expression.
        if name is None:
            name = ''

        if lower_bound is None:
            lower_bound = -self.model.infinity()
        if upper_bound is None:
            upper_bound = self.model.infinity()

        cons = lower_bound <= (linfun <= upper_bound)
        self.model.addCons(cons, name=name)
        self.constraints = None

    cpdef row(self, int index):
        r"""
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
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)
        """
        namedvars = [var.name for var in self.variables]
        valslinear = self.model.getValsLinear(self.get_constraints()[index])
        cdef list indices = []
        cdef list values = []
        for var, coeff in valslinear.iteritems():
            indices.append(namedvars.index(var))
            values.append(coeff)
        return (indices, values)

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
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row_bounds(0)
            (2.0, 2.0)
        """
        cons = self.get_constraints()[index]
        lhs = self.model.getLhs(cons)
        rhs = self.model.getRhs(cons)
        if lhs == -self.model.infinity():
            lhs = None
        if rhs == self.model.infinity():
            rhs = None
        return (lhs, rhs)

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
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0.0, 5.0)
        """
        var = self.variables[index]
        lb = var.getLbOriginal()
        if lb == -self.model.infinity():
            lb = None
        ub = var.getUbOriginal()
        if ub == self.model.infinity():
            ub = None
        return (lb, ub)

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
            sage: p = get_solver(solver = "SCIP")
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5
        """
        mcons = self.get_constraints()
        self.constraints = None  # is this necessary?
        # after update of
        index = self.add_variable(lower_bound=-self.model.infinity())
        var = self.variables[index]

        for i, coeff in zip(indices, coeffs):
            self.model.addConsCoeff(var=var, cons=mcons[i], coeff=coeff)

    cpdef int solve(self) except -1:
        """
        Solve the problem.

        Sage uses SCIP's implementation of the branch-and-cut algorithm
        to solve the mixed-integer linear program.
        (If all variables are continuous, the algorithm reduces to solving the
        linear program by the simplex method.)

        EXAMPLES::

            sage: lp = MixedIntegerLinearProgram(solver = 'SCIP', maximization = False)
            sage: x, y = lp[0], lp[1]
            sage: lp.add_constraint(-2*x + y <= 1)
            sage: lp.add_constraint(x - y <= 1)
            sage: lp.add_constraint(x + y >= 2)
            sage: lp.set_objective(x + y)
            sage: lp.set_integer(x)
            sage: lp.set_integer(y)
            sage: lp.solve()
            2.0
            sage: lp.get_values([x, y])
            [1.0, 1.0]

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.solve()
            0
            sage: p.objective_coefficient(0,1)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: ...

            sage: lp = MixedIntegerLinearProgram(solver = "SCIP")
            sage: v = lp.new_variable(nonnegative=True)
            sage: lp.add_constraint(v[1] +v[2] -2.0 *v[3], max=-1.0)
            sage: lp.add_constraint(v[0] -4.0/3 *v[1] +1.0/3 *v[2], max=-1.0/3)
            sage: lp.add_constraint(v[0] +0.5 *v[1] -0.5 *v[2] +0.25 *v[3], max=-0.25)
            sage: lp.solve()
            0.0

        Solving a LP within the acceptable gap. No exception is raised, even if
        the result is not optimal. To do this, we try to compute the maximum
        number of disjoint balls (of diameter 1) in a hypercube::

            sage: g = graphs.CubeGraph(9)
            sage: p = MixedIntegerLinearProgram(solver = "SCIP")

            sage: b = p.new_variable(binary=True)
            sage: p.set_objective(p.sum(b[v] for v in g))
            sage: for v in g:
            ....:     p.add_constraint(b[v] + p.sum(b[u] for u in g.neighbors(v)) <= 1)
            sage: p.add_constraint(b[v] == 1) # Force an easy non-0 solution
            sage: p.solver_parameter("limits/absgap", 100)
            sage: p.solve() # rel tol 100
            1
        """
        if self.model.getStatus() != 'unknown' or self.model.getStage() > 1:
            # This should actually be self.model.freeReoptSolve, but it seems to fail.
            self.model.freeTransform()
        self.model.optimize()

        status = self.model.getStatus()

        if status == 'unbounded':
            raise MIPSolverException("SCIP: Solution is unbounded")
        elif status == 'infeasible':
            raise MIPSolverException("SCIP: There is no feasible solution")
        # elif status == 'timelimit':
        #     raise MIPSolverException("SCIP: Time limit reached")
        return 0

    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            7.5
            sage: p.get_variable_value(0) # abs tol 1e-15
            0.0
            sage: p.get_variable_value(1)
            1.5
        """
        return self.model.getObjVal()

    cpdef best_known_objective_bound(self):
        r"""
        Return the value of the currently best known bound.

        This method returns the current best upper (resp. lower) bound on the
        optimal value of the objective function in a maximization
        (resp. minimization) problem. It is equal to the output of
        :meth:`get_objective_value` if the MILP found an optimal solution, but
        it can differ if it was interrupted manually or after a time limit (cf
        :meth:`solver_parameter`).

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: g = graphs.CubeGraph(9)
            sage: p = MixedIntegerLinearProgram(solver='SCIP')
            sage: p.solver_parameter("limits/gap", 100)
            sage: b = p.new_variable(binary=True)
            sage: p.set_objective(p.sum(b[v] for v in g))
            sage: for v in g:
            ....:     p.add_constraint(b[v]+p.sum(b[u] for u in g.neighbors(v)) <= 1)
            sage: p.add_constraint(b[v] == 1) # Force an easy non-0 solution
            sage: p.solve() # rel tol 100
            1.0
            sage: backend = p.get_backend()
            sage: backend.best_known_objective_bound() # random
            31.0
        """
        return self.model.getPrimalbound()

    cpdef get_relative_objective_gap(self):
        r"""
        Return the relative objective gap of the best known solution.

        For a minimization problem, this value is computed by
        `(\texttt{bestinteger} - \texttt{bestobjective}) / (1e-10 +
        |\texttt{bestobjective}|)`, where ``bestinteger`` is the value returned
        by :meth:`get_objective_value` and ``bestobjective`` is the value
        returned by :meth:`best_known_objective_bound`. For a maximization
        problem, the value is computed by `(\texttt{bestobjective} -
        \texttt{bestinteger}) / (1e-10 + |\texttt{bestobjective}|)`.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: g = graphs.CubeGraph(9)
            sage: p = MixedIntegerLinearProgram(solver='SCIP')
            sage: p.solver_parameter("limits/gap", 100)
            sage: b = p.new_variable(binary=True)
            sage: p.set_objective(p.sum(b[v] for v in g))
            sage: for v in g:
            ....:     p.add_constraint(b[v]+p.sum(b[u] for u in g.neighbors(v)) <= 1)
            sage: p.add_constraint(b[v] == 1) # Force an easy non-0 solution
            sage: p.solve() # rel tol 100
            1.0
            sage: backend = p.get_backend()
            sage: backend.get_relative_objective_gap() # random
            46.99999999999999

        TESTS:

        Just make sure that the variable *has* been defined, and is not just
        undefined::

            sage: backend.get_relative_objective_gap() > 1                              # needs sage.graphs
            True
        """
        return self.model.getGap()

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            7.5
            sage: p.get_variable_value(0) # abs tol 1e-15
            0.0
            sage: p.get_variable_value(1)
            1.5
        """
        return self.model.getVal(self.variables[variable])

    cpdef get_row_prim(self, int i):
        r"""
        Return the value of the auxiliary variable associated with i-th row.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "SCIP")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(zip([0, 1, 2], [8, 6, 1]), None, 48)
            sage: lp.add_linear_constraint(zip([0, 1, 2], [4, 2, 1.5]), None, 20)
            sage: lp.add_linear_constraint(zip([0, 1, 2], [2, 1.5, 0.5]), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: lp.solver_parameter('presolving/maxrounds',0)
            sage: lp.solve()
            0
            sage: lp.get_objective_value()
            280.0
            sage: lp.get_row_prim(0)
            24.0
            sage: lp.get_row_prim(1)
            20.0
            sage: lp.get_row_prim(2)
            8.0
        """
        return self.model.getActivity(self.get_constraints()[i])

    cpdef int ncols(self) noexcept:
        """
        Return the number of columns/variables.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
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
            sage: p = get_solver(solver = "SCIP")
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(2, 2, None)
            sage: p.nrows()
            2

        TESTS::

        After calling :meth:`remove_constraints` we know that
        `self.constraints is None`.  `SCIP` keeps track of the number
        of constraints, so we can do the optimization::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='SCIP')
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0, 2), (1, 3)], None, 6)
            sage: p.row(0)
            ([0, 1], [2.0, 3.0])
            sage: p.remove_constraints([0])
            sage: p.nrows()
            0
        """
        if self.constraints is None:
            return self.model.getNConss()
        return len(self.get_constraints())

    cpdef col_name(self, int index):
        """
        Return the ``index``-th col name.

        INPUT:

        - ``index`` -- integer; the col's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variable(name='I am a variable')
            0
            sage: p.col_name(0)
            'I am a variable'
        """
        return self.variables[index].name

    cpdef row_name(self, int index):
        """
        Return the ``index``-th row name.

        INPUT:

        - ``index`` -- integer; the row's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_linear_constraints(1, 2, None, names=['Empty constraint 1'])
            sage: p.row_name(0)
            'Empty constraint 1'
        """
        return self.get_constraints()[index].name

    cpdef bint is_variable_binary(self, int index) noexcept:
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,0)
            sage: p.is_variable_binary(0)
            True
        """
        return self.variables[index].vtype() == 'BINARY'

    cpdef bint is_variable_integer(self, int index) noexcept:
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)
            True
        """
        return self.variables[index].vtype() == 'INTEGER'

    cpdef bint is_variable_continuous(self, int index) noexcept:
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_continuous(0)
            True
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_continuous(0)
            False
        """
        return self.variables[index].vtype() == 'CONTINUOUS'

    cpdef bint is_maximization(self) noexcept:
        """
        Test whether the problem is a maximization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        return self.model.getObjectiveSense() != 'minimize'

    cpdef variable_upper_bound(self, int index, value=False):
        """
        Return or define the upper bound on a variable.

        INPUT:

        - ``index`` -- integer; the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0.0, 5.0)

        TESTS:

        :issue:`14581`::

            sage: P = MixedIntegerLinearProgram(solver='SCIP')
            sage: x = P["x"]
            sage: P.set_max(x, 0)
            sage: P.get_max(x)
            0.0

        Check that :issue:`10232` is fixed::

            sage: p = get_solver(solver='SCIP')
            sage: p.variable_upper_bound(2)
            Traceback (most recent call last):
            ...
            ValueError: The variable's id must satisfy 0 <= id < number_of_variables
            sage: p.variable_upper_bound(3, 5)
            Traceback (most recent call last):
            ...
            ValueError: The variable's id must satisfy 0 <= id < number_of_variables

            sage: p.add_variable()
            0
            sage: p.variable_upper_bound(0, 'hey!')
            Traceback (most recent call last):
            ...
            TypeError: must be real number, not str
        """
        if index < 0 or index >= self.ncols():
            raise ValueError("The variable's id must satisfy 0 <= id < number_of_variables")
        var = self.variables[index]
        if value is False:
            return var.getUbOriginal()
        else:
            self.model.chgVarUb(var=var, ub=value)

    cpdef variable_lower_bound(self, int index, value=False):
        """
        Return or define the lower bound on a variable.

        INPUT:

        - ``index`` -- integer; the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)
            sage: p.col_bounds(0)
            (5.0, None)

        TESTS:

        :issue:`14581`::

            sage: P = MixedIntegerLinearProgram(solver='SCIP')
            sage: x = P["x"]
            sage: P.set_min(x, 5)
            sage: P.set_min(x, 0)
            sage: P.get_min(x)
            0.0

        Check that :issue:`10232` is fixed::

            sage: p = get_solver(solver='SCIP')
            sage: p.variable_lower_bound(2)
            Traceback (most recent call last):
            ...
            ValueError: The variable's id must satisfy 0 <= id < number_of_variables
            sage: p.variable_lower_bound(3, 5)
            Traceback (most recent call last):
            ...
            ValueError: The variable's id must satisfy 0 <= id < number_of_variables

            sage: p.add_variable()
            0
            sage: p.variable_lower_bound(0, 'hey!')
            Traceback (most recent call last):
            ...
            TypeError: must be real number, not str
        """
        if index < 0 or index >= self.ncols():
            raise ValueError("The variable's id must satisfy 0 <= id < number_of_variables")
        var = self.variables[index]
        if value is False:
            return var.getLbOriginal()
        else:
            self.model.chgVarLb(var=var, lb=value)

    cpdef write_cip(self, filename):
        """
        Write the problem to a ``.cip`` file.

        INPUT:

        - ``filename`` -- string

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
            sage: p.set_objective([2, 5])
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(suffix='.cip') as f:
            ....:     p.write_cip(f.name)
            wrote problem to file ...
        """
        self.model.writeProblem(filename)

    cpdef write_lp(self, filename):
        """
        Write the problem to a ``.lp`` file.

        INPUT:

        - ``filename`` -- string

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
            sage: p.set_objective([2, 5])
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(suffix='.lp') as f:
            ....:     p.write_lp(f.name)
            wrote problem to file ...
        """
        filenamestr = filename
        fname, fext = splitext(filenamestr)

        if fext.lower() != 'lp':
            filenamestr = filename + '.lp'

        self.model.writeProblem(filenamestr)

    cpdef write_mps(self, filename, int modern):
        """
        Write the problem to a ``.mps`` file.

        INPUT:

        - ``filename`` -- string

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "SCIP")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
            sage: p.set_objective([2, 5])
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(suffix='.mps') as f:
            ....:     p.write_mps(f.name, 2)
            wrote problem to file ...
        """
        filenamestr = filename
        fname, fext = splitext(filenamestr)

        if fext.lower() != 'mps':
            filenamestr = filename + '.mps'

        self.model.writeProblem(filenamestr)

    cpdef __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "SCIP")
            sage: b = p.new_variable()
            sage: p.add_constraint(b[1] + b[2] <= 6)
            sage: p.set_objective(b[1] + b[2])
            sage: copy(p).solve()
            6.0
        """
        cdef SCIPBackend cp = type(self)(maximization=self.is_maximization())
        cp.model = Model(sourceModel=self.model, origcopy=True)
        cp.problem_name(self.problem_name())
        cp.obj_constant_term = self.obj_constant_term
        cp.variables = cp.model.getVars()
        return cp

    cpdef solver_parameter(self, name, value=None):
        """
        Return or define a solver parameter.

        INPUT:

        - ``name`` -- string; the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value

        EXAMPLES:

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='SCIP')
            sage: p.solver_parameter("limits/time", 1)
            sage: p.solver_parameter("limits/time")
            1.0
        """
        if value is not None:
            if name.lower() == 'timelimit':
                self.model.setRealParam("limits/time", float(value))
            else:
                self.model.setParam(name, value)
        else:
            return self.model.getParam(name)
