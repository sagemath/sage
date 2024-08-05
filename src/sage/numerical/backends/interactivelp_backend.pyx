r"""
InteractiveLP Backend

AUTHORS:

- Nathann Cohen (2010-10)      : generic_backend template
- Matthias Koeppe (2016-03)    : this backend
"""

# ****************************************************************************
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#       Copyright (C) 2016 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.numerical.mip import MIPSolverException
from sage.numerical.interactive_simplex_method import InteractiveLPProblem, default_variable_name
from sage.modules.free_module_element import vector
from copy import copy


cdef class InteractiveLPBackend:
    """
    MIP Backend that works with :class:`InteractiveLPProblem`.

    This backend should be used only for linear programs over general fields,
    or for educational purposes.  For fast computations with floating point
    arithmetic, use one of the numerical backends. For exact computations
    with rational numbers, use backend 'PPL'.

    There is no support for integer variables.

    EXAMPLES::

        sage: from sage.numerical.backends.generic_backend import get_solver
        sage: p = get_solver(solver = "InteractiveLP")
    """

    def __cinit__(self, maximization = True, base_ring = None):
        """
        Cython constructor.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")

        This backend can work with irrational algebraic numbers::

            sage: # needs sage.rings.number_field
            sage: poly = polytopes.dodecahedron(base_ring=AA)
            sage: lp, x = poly.to_linear_program(solver='InteractiveLP', return_variable=True)
            sage: lp.set_objective(x[0] + x[1] + x[2])
            sage: lp.solve()
            2.291796067500631?
            sage: lp.get_values(x[0], x[1], x[2])
            [0.763932022500211?, 0.763932022500211?, 0.763932022500211?]
            sage: lp.set_objective(x[0] - x[1] - x[2])
            sage: lp.solve()
            2.291796067500631?
            sage: lp.get_values(x[0], x[1], x[2])
            [0.763932022500211?, -0.763932022500211?, -0.763932022500211?]
        """

        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ

        self.lp = InteractiveLPProblem([], [], [], base_ring=base_ring)
        self.set_verbosity(0)

        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

        self.row_names = []

    cpdef __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "InteractiveLP")
            sage: b = p.new_variable()
            sage: p.add_constraint(b[1] + b[2] <= 6)
            sage: p.set_objective(b[1] + b[2])
            sage: cp = copy(p.get_backend())
            sage: cp.solve()
            0
            sage: cp.get_objective_value()
            6
        """
        cdef InteractiveLPBackend cp = type(self)(base_ring=self.base_ring())
        cp.lp = self.lp                   # it's considered immutable; so no need to copy.
        cp.row_names = copy(self.row_names)
        cp.prob_name = self.prob_name
        return cp

    cpdef base_ring(self):
        """
        Return the base ring.

        OUTPUT: a ring; the coefficients that the chosen solver supports

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.base_ring()
            Rational Field
        """
        return self.lp.base_ring()

    def _variable_type_from_bounds(self, lower_bound, upper_bound):
        """
        Return a string designating a variable type in `InteractiveLPProblem`.

        INPUT:

        - ``lower_bound`` -- the lower bound of the variable

        - ``upper_bound`` -- the upper bound of the variable

        OUTPUT: string, one of ``''``, ``'<='``, ``'>='``

        The function raises an error if this pair of bounds cannot be
        represented by an `InteractiveLPProblem` variable type.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p._variable_type_from_bounds(0, None)
            '>='
            sage: p._variable_type_from_bounds(None, 0)
            '<='
            sage: p._variable_type_from_bounds(None, None)
            ''
            sage: p._variable_type_from_bounds(None, 5)
            Traceback (most recent call last):
            ...
            NotImplementedError: General variable bounds not supported
        """
        if lower_bound is None:
            if upper_bound is None:
                return ""
            elif upper_bound == 0:
                return "<="
            else:
                raise NotImplementedError("General variable bounds not supported")
        elif lower_bound == 0:
            if upper_bound is None:
                return ">="
            else:
                raise NotImplementedError("General variable bounds not supported")
        else:
            raise NotImplementedError("General variable bounds not supported")

    cpdef int add_variable(self, lower_bound=0, upper_bound=None,
                           binary=False, continuous=True, integer=False,
                           obj=None, name=None, coefficients=None) except -1:
        ## coefficients is an extension in this backend,
        ## and a proposed addition to the interface, to unify this with add_col.
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both nonnegative and real.

        In this backend, variables are always continuous (real).
        If integer variables are requested via the parameters
        ``binary`` and ``integer``, an error will be raised.

        INPUT:

        - ``lower_bound`` -- the lower bound of the variable (default: 0)

        - ``upper_bound`` -- the upper bound of the variable (default: ``None``)

        - ``binary`` -- ``True`` if the variable is binary (default: ``False``)

        - ``continuous`` -- ``True`` if the variable is continuous (default: ``True``)

        - ``integer`` -- ``True`` if the variable is integral (default: ``False``)

        - ``obj`` -- (optional) coefficient of this variable in the objective function (default: 0)

        - ``name`` -- an optional name for the newly added variable (default: ``None``)

        - ``coefficients`` -- (optional) an iterable of pairs ``(i, v)``. In each
          pair, ``i`` is a variable index (integer) and ``v`` is a
          value (element of :meth:`base_ring`).

        OUTPUT: the index of the newly created variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
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
            1
        """
        A, b, c, x, constraint_types, variable_types, problem_type, ring, d = self._AbcxCVPRd()
        cdef int vtype = int(binary) + int(continuous) + int(integer)
        if vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")
        if not continuous:
            raise NotImplementedError("Integer variables are not supported")
        variable_types = variable_types + (self._variable_type_from_bounds(lower_bound, upper_bound),)
        col = vector(ring, self.nrows())
        if coefficients is not None:
            for (i, v) in coefficients:
                col[i] = v
        A = A.augment(col)
        if obj is None:
            obj = 0
        c = tuple(c) + (obj,)
        if name is None:
            var_names = default_variable_name("primal decision")
            name = "{}{:d}".format(var_names, self.ncols() + 1)
        x = tuple(x) + (name,)
        self.lp = InteractiveLPProblem(A, b, c, x,
                                       constraint_types, variable_types,
                                       problem_type, ring, objective_constant_term=d)
        return self.ncols() - 1

    cpdef  set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable.

        In this backend, variables are always continuous (real).
        If integer or binary variables are requested via the parameter
        ``vtype``, an error will be raised.

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``vtype`` -- integer:

            *  1  Integer
            *  0  Binary
            *  -1  Continuous

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,-1)
            sage: p.is_variable_continuous(0)
            True
        """
        if vtype == -1:
            pass
        else:
            raise NotImplementedError()

    def _AbcxCVPRd(self):
        """
        Retrieve all problem data from the LP.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p._AbcxCVPRd()
            ([], (), (), (), (), (), 'max', Rational Field, 0)
        """
        A, b, c, x = self.lp.Abcx()
        constraint_types = self.lp.constraint_types()
        variable_types = self.lp.variable_types()
        problem_type = self.lp.problem_type()
        base_ring = self.lp.base_ring()
        d = self.lp.objective_constant_term()
        return A, b, c, x, constraint_types, variable_types, problem_type, base_ring, d

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` -- integer:

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        A, b, c, x, constraint_types, variable_types, problem_type, ring, d = self._AbcxCVPRd()
        if sense == +1:
            problem_type = "max"
        else:
            problem_type = "min"
        self.lp = InteractiveLPProblem(A, b, c, x,
                                       constraint_types, variable_types,
                                       problem_type, ring, objective_constant_term=d)

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective
        function

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``coeff`` -- double; its coefficient

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2
        """
        if coeff is None:
            return self.lp.objective_coefficients()[variable]
        else:
            A, b, c, x, constraint_types, variable_types, problem_type, ring, d = self._AbcxCVPRd()
            c = list(c)
            c[variable] = coeff
            self.lp = InteractiveLPProblem(A, b, c, x,
                                           constraint_types, variable_types,
                                           problem_type, ring, objective_constant_term=d)

    cpdef objective_constant_term(self, d=None):
        """
        Set or get the constant term in the objective function.

        INPUT:

        - ``d`` -- double; its coefficient. If ``None`` (default), return the
          current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.objective_constant_term()
            0
            sage: p.objective_constant_term(42)
            sage: p.objective_constant_term()
            42
        """
        if d is None:
            return self.lp.objective_constant_term()
        else:
            A, b, c, x, constraint_types, variable_types, problem_type, ring, _ = self._AbcxCVPRd()
            self.lp = InteractiveLPProblem(A, b, c, x,
                                           constraint_types, variable_types,
                                           problem_type, ring, objective_constant_term=d)

    cpdef set_objective(self, list coeff, d=0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- list of real values, whose i-th element is the
          coefficient of the i-th variable in the objective function

        - ``d`` -- real; the constant term in the linear function (set to `0`
          by default)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1, 1, 2, 1, 3]

        Constants in the objective function are respected::

            sage: p = MixedIntegerLinearProgram(solver='InteractiveLP')
            sage: x,y = p[0], p[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: p.solve()
            47/5

        TESTS:

        Constants also work the right way for min::

            sage: p = MixedIntegerLinearProgram(solver='InteractiveLP', maximization=False)
            sage: x,y = p[0], p[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(-x - y - 7)
            sage: p.solve()
            -47/5
        """
        A, b, _, x, constraint_types, variable_types, problem_type, ring, _ = self._AbcxCVPRd()
        c = coeff
        self.lp = InteractiveLPProblem(A, b, c, x,
                                       constraint_types, variable_types,
                                       problem_type, ring, objective_constant_term=d)

    cpdef set_verbosity(self, int level):
        """
        Set the log (verbosity) level.

        INPUT:

        - ``level`` -- integer; from 0 (no verbosity) to 3

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.set_verbosity(2)
        """
        self.verbosity = level

    cpdef remove_constraint(self, int i):
        r"""
        Remove a constraint.

        INPUT:

        - ``i`` -- index of the constraint to remove

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='InteractiveLP')
            sage: v = p.new_variable(nonnegative=True)
            sage: x,y = v[0], v[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: p.solve()
            47/5
            sage: p.remove_constraint(0)
            sage: p.solve()
            10
            sage: p.get_values([x,y])
            [0, 3]
        """
        A, b, c, x, constraint_types, variable_types, problem_type, ring, d = self._AbcxCVPRd()
        A = A.delete_rows((i,))
        b = list(b)
        del b[i]
        constraint_types = list(constraint_types)
        del constraint_types[i]
        self.lp = InteractiveLPProblem(A, b, c, x,
                                       constraint_types, variable_types,
                                       problem_type, ring, objective_constant_term=d)

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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint( zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)
            (2, 2)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1, 1, name='foo')
            sage: p.row_name(1)
            'foo'
        """
        A, b, c, x, constraint_types, variable_types, problem_type, ring, d = self._AbcxCVPRd()
        if lower_bound is None:
            if upper_bound is None:
                raise ValueError("At least one of lower_bound and upper_bound must be provided")
            else:
                constraint_types = constraint_types + ("<=",)
                b = tuple(b) + (upper_bound,)
        else:
            if upper_bound is None:
                constraint_types = constraint_types + (">=",)
                b = tuple(b) + (lower_bound,)
            elif lower_bound == upper_bound:
                constraint_types = constraint_types + ("==",)
                b = tuple(b) + (lower_bound,)
            else:
                raise NotImplementedError("Ranged constraints are not supported")

        row = vector(ring, self.ncols())
        for (i, v) in coefficients:
            row[i] = v
        A = A.stack(row)

        self.row_names.append(name)

        self.lp = InteractiveLPProblem(A, b, c, x,
                                       constraint_types, variable_types,
                                       problem_type, ring, objective_constant_term=d)

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

            ``indices`` and ``coeffs`` are expected to be of the same length.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(list(range(5)), list(range(5)))
            sage: p.nrows()
            5
        """
        self.add_variable(coefficients=zip(indices, coeffs))

    cpdef int solve(self) except -1:
        """
        Solve the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution cannot be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(list(range(5)), list(range(5)))
            sage: p.solve()
            0
            sage: p.objective_coefficient(0,1)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: ...
        """
        ## FIXME: standard_form should allow to pass slack names (which we would take from row_names).
        ## FIXME: Perhaps also pass the problem name as objective name
        lp_std_form, transformation = self.lp.standard_form(transformation=True)
        self.lp_std_form, self.std_form_transformation = lp_std_form, transformation
        output = lp_std_form.run_revised_simplex_method()
        ## FIXME: Display output as a side effect if verbosity is high enough
        d = self.final_dictionary = lp_std_form.final_revised_dictionary()
        if d.is_optimal():
            if lp_std_form.auxiliary_variable() in d.basic_variables():
                raise MIPSolverException("InteractiveLP: Problem has no feasible solution")
            else:
                return 0
        else:
            raise MIPSolverException("InteractiveLP: Problem is unbounded")

    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behavior is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0,1), (1,2)], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            15/2
            sage: p.get_variable_value(0)
            0
            sage: p.get_variable_value(1)
            3/2
        """
        d = self.final_dictionary
        v = d.objective_value()
        if self.lp_std_form.is_negative():
            v = - v
        return v

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behavior is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0,1), (1, 2)], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            15/2
            sage: p.get_variable_value(0)
            0
            sage: p.get_variable_value(1)
            3/2
        """
        solution = self.std_form_transformation(self.final_dictionary.basic_solution())
        return solution[variable]

    cpdef int ncols(self) noexcept:
        """
        Return the number of columns/variables.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variables(2)
            1
            sage: p.ncols()
            2
        """
        return self.lp.n_variables()

    cpdef int nrows(self) noexcept:
        """
        Return the number of rows/constraints.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(2, 0, None)
            sage: p.nrows()
            2
        """
        return self.lp.n_constraints()

    cpdef bint is_maximization(self) noexcept:
        """
        Test whether the problem is a maximization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        return self.lp.problem_type() == "max"

    cpdef problem_name(self, name=None):
        """
        Return or define the problem's name.

        INPUT:

        - ``name`` -- string; the problem's name. When set to
          ``None`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 0, None)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
        """
        A, b, c, x = self.lp.Abcx()
        indices = []
        coeffs = []
        for j in range(self.ncols()):
            if A[i][j] != 0:
                indices.append(j)
                coeffs.append(A[i][j])
        return (indices, coeffs)

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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row_bounds(0)
            (2, 2)
        """
        A, b, c, x = self.lp.Abcx()
        constraint_types = self.lp.constraint_types()
        if constraint_types[index] == '>=':
            return (b[index], None)
        elif constraint_types[index] == '<=':
            return (None, b[index])
        elif constraint_types[index] == '==':
            return (b[index], b[index])
        else:
            raise ValueError("Bad constraint_type")

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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variable(lower_bound=None)
            0
            sage: p.col_bounds(0)
            (None, None)
            sage: p.variable_lower_bound(0, 0)
            sage: p.col_bounds(0)
            (0, None)
        """
        t = self.lp.variable_types()[index]
        if t == ">=":
            return (0, None)
        elif t == "<=":
            return (None, 0)
        elif t == "":
            return (None, None)
        else:
            raise ValueError("Bad _variable_types")

    cpdef bint is_variable_binary(self, int index) noexcept:
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_binary(0)
            False
        """
        return False

    cpdef bint is_variable_integer(self, int index) noexcept:
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_integer(0)
            False
        """
        return False

    cpdef bint is_variable_continuous(self, int index) noexcept:
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_continuous(0)
            True
        """
        return True

    cpdef row_name(self, int index):
        """
        Return the ``index``-th row name.

        INPUT:

        - ``index`` -- integer; the row's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_linear_constraints(1, 2, None, names=['Empty constraint 1'])
            sage: p.row_name(0)
            'Empty constraint 1'
        """
        return self.row_names[index]

    cpdef col_name(self, int index):
        """
        Return the ``index``-th column name.

        INPUT:

        - ``index`` -- integer; the column id

        - ``name`` -- (``char *``) its name; when set to ``NULL``
          (default), the method returns the current name

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variable(name='I_am_a_variable')
            0
            sage: p.col_name(0)
            'I_am_a_variable'
        """
        return str(self.lp.decision_variables()[index])

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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variable(lower_bound=None)
            0
            sage: p.col_bounds(0)
            (None, None)
            sage: p.variable_upper_bound(0) is None
            True
            sage: p.variable_upper_bound(0, 0)
            sage: p.col_bounds(0)
            (None, 0)
            sage: p.variable_upper_bound(0)
            0
            sage: p.variable_upper_bound(0, None)
            sage: p.variable_upper_bound(0) is None
            True
        """
        bounds = self.col_bounds(index)
        if value is False:
            return bounds[1]
        else:
            if value != bounds[1]:
                bounds = (bounds[0], value)
                A, b, c, x, constraint_types, variable_types, problem_type, ring, d = self._AbcxCVPRd()
                variable_types = list(variable_types)
                variable_types[index] = self._variable_type_from_bounds(*bounds)
                self.lp = InteractiveLPProblem(A, b, c, x,
                                               constraint_types, variable_types,
                                               problem_type, ring, objective_constant_term=d)

    cpdef variable_lower_bound(self, int index, value=False):
        """
        Return or define the lower bound on a variable.

        INPUT:

        - ``index`` -- integer; the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has no lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variable(lower_bound=None)
            0
            sage: p.col_bounds(0)
            (None, None)
            sage: p.variable_lower_bound(0) is None
            True
            sage: p.variable_lower_bound(0, 0)
            sage: p.col_bounds(0)
            (0, None)
            sage: p.variable_lower_bound(0)
            0
            sage: p.variable_lower_bound(0, None)
            sage: p.variable_lower_bound(0) is None
            True
        """
        bounds = self.col_bounds(index)
        if value is False:
            return bounds[0]
        else:
            if value != bounds[0]:
                bounds = (value, bounds[1])
                A, b, c, x, constraint_types, variable_types, problem_type, ring, d = self._AbcxCVPRd()
                variable_types = list(variable_types)
                variable_types[index] = self._variable_type_from_bounds(*bounds)
                self.lp = InteractiveLPProblem(A, b, c, x,
                                               constraint_types, variable_types,
                                               problem_type, ring, objective_constant_term=d)

    cpdef bint is_variable_basic(self, int index) noexcept:
        """
        Test whether the given variable is basic.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver='InteractiveLP')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_variable_basic(0)
            True
            sage: b.is_variable_basic(1)
            False
        """
        return self.lp_std_form.decision_variables()[index] in self.final_dictionary.basic_variables()

    cpdef bint is_variable_nonbasic_at_lower_bound(self, int index) noexcept:
        """
        Test whether the given variable is nonbasic at lower bound.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver='InteractiveLP')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_variable_nonbasic_at_lower_bound(0)
            False
            sage: b.is_variable_nonbasic_at_lower_bound(1)
            True
        """
        return self.lp_std_form.decision_variables()[index] in self.final_dictionary.nonbasic_variables()

    cpdef bint is_slack_variable_basic(self, int index) noexcept:
        """
        Test whether the slack variable of the given row is basic.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver='InteractiveLP')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_slack_variable_basic(0)
            True
            sage: b.is_slack_variable_basic(1)
            False
        """
        return self.lp_std_form.slack_variables()[index] in self.final_dictionary.basic_variables()

    cpdef bint is_slack_variable_nonbasic_at_lower_bound(self, int index) noexcept:
        """
        Test whether the given variable is nonbasic at lower bound.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver='InteractiveLP')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_slack_variable_nonbasic_at_lower_bound(0)
            False
            sage: b.is_slack_variable_nonbasic_at_lower_bound(1)
            True
        """
        return self.lp_std_form.slack_variables()[index] in self.final_dictionary.nonbasic_variables()

    cpdef dictionary(self):
        # Proposed addition to the general interface,
        # which would for other solvers return backend dictionaries (#18804)
        """
        Return a dictionary representing the current basis.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver='InteractiveLP')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: d = b.dictionary(); d
            LP problem dictionary ...
            sage: set(d.basic_variables())
            {x1, x3}
            sage: d.basic_solution()
            (17/8, 0)

        TESTS:

        Compare with a dictionary obtained in another way::

            sage: lp, basis = p.interactive_lp_problem()
            sage: lp.dictionary(*basis).basic_solution()
            (17/8, 0)
        """
        return self.final_dictionary

    cpdef interactive_lp_problem(self):
        """
        Return the :class:`InteractiveLPProblem` object associated with this backend.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver='InteractiveLP')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: b.interactive_lp_problem()
            LP problem ...
        """
        return self.lp
