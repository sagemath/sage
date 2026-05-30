# sage.doctest: optional - pplpy
"""
PPL Backend

AUTHORS:

- Risan (2012-02): initial implementation

- Jeroen Demeyer (2014-08-04) allow rational coefficients for
  constraints and objective function (:issue:`16755`)
"""

# ****************************************************************************
#       Copyright (C) 2010 Risan <ptrrsn.1@gmail.com>
#       Copyright (C) 2014 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.numerical.mip import MIPSolverException
from ppl import MIP_Problem, Variable, Variables_Set, Linear_Expression
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.numerical.backends.generic_backend cimport GenericBackend
from copy import copy


cdef class PPLBackend(GenericBackend):
    """
    MIP Backend that uses the exact MIP solver from the Parma Polyhedra Library.
    """

    cdef object mip
    cdef list Matrix
    cdef list row_lower_bound
    cdef list row_upper_bound
    cdef list col_lower_bound
    cdef list col_upper_bound
    cdef list objective_function
    cdef list row_name_var
    cdef list col_name_var
    cdef int is_maximize
    cdef str name
    cdef object integer_variables

    # Common denominator for objective function in self.mip (not for the constant term)
    cdef Integer obj_denominator

    def __cinit__(self, maximization=True, base_ring=None):
        """
        Constructor.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver = "PPL")

        TESTS:

        Raise an error if a ``base_ring`` is requested that is not supported::

            sage: p = MixedIntegerLinearProgram(solver='PPL', base_ring=AA)             # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            TypeError: The PPL backend only supports rational data.
        """

        if base_ring is not None:
            from sage.rings.rational_field import QQ
            if base_ring is not QQ:
                raise TypeError('The PPL backend only supports rational data.')

        self.Matrix = []
        self.row_lower_bound = []
        self.row_upper_bound = []
        self.col_lower_bound = []
        self.col_upper_bound = []
        self.objective_function = []
        self.row_name_var = []
        self.col_name_var = []
        self.name = ''
        self.obj_constant_term = Rational(0)
        self.obj_denominator = Integer(1)
        self.integer_variables = set()

        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

    cpdef base_ring(self):
        from sage.rings.rational_field import QQ
        return QQ

    cpdef zero(self):
        return self.base_ring()(0)

    cpdef __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "PPL")
            sage: b = p.new_variable()
            sage: p.add_constraint(b[1] + b[2] <= 6)
            sage: p.set_objective(b[1] + b[2])
            sage: cp = copy(p.get_backend())
            sage: cp.solve()
            0
            sage: cp.get_objective_value()
            6
        """
        cdef PPLBackend cp = type(self)()
        cp.Matrix = [row[:] for row in self.Matrix]
        cp.row_lower_bound = self.row_lower_bound[:]
        cp.row_upper_bound = self.row_upper_bound[:]
        cp.col_lower_bound = self.col_lower_bound[:]
        cp.col_upper_bound = self.col_upper_bound[:]
        cp.objective_function = self.objective_function[:]
        cp.row_name_var = self.row_name_var[:]
        cp.col_name_var = self.col_name_var[:]
        cp.name = self.name
        cp.obj_constant_term = self.obj_constant_term
        cp.obj_denominator = self.obj_denominator
        cp.integer_variables = copy(self.integer_variables)
        cp.is_maximize = self.is_maximize
        return cp

    def init_mip(self):
        """
        Converting the matrix form of the MIP Problem to PPL MIP_Problem.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='PPL')
            sage: p.base_ring()
            Rational Field
            sage: type(p.zero())
            <class 'sage.rings.rational.Rational'>
            sage: p.init_mip()
        """

        self.mip = MIP_Problem()

        # Common denominator (for objective function and every constraint)
        cdef Integer denom, newdenom

        self.mip.add_space_dimensions_and_embed(len(self.objective_function))

        # Integrality

        ivar = Variables_Set()
        for i in self.integer_variables:
            ivar.insert(Variable(i))
        self.mip.add_to_integer_space_dimensions(ivar)

        # Objective function
        mip_obj = Linear_Expression(0)
        denom = Integer(1)
        for i in range(len(self.objective_function)):
            coeff = self.objective_function[i] * denom
            newdenom = coeff.denominator()
            if newdenom != 1:
                assert newdenom >= 2
                denom *= newdenom
                mip_obj *= newdenom
                coeff *= newdenom
            mip_obj = mip_obj + Linear_Expression(coeff * Variable(i))
        self.mip.set_objective_function(mip_obj)
        self.obj_denominator = denom

        # Constraints
        for i in range(len(self.Matrix)):
            l = Linear_Expression(0)
            denom = Integer(1)
            for j in range(len(self.Matrix[i])):
                coeff = self.Matrix[i][j] * denom
                newdenom = coeff.denominator()
                if newdenom != 1:
                    assert newdenom >= 2
                    denom *= newdenom
                    l *= newdenom
                    coeff *= newdenom
                l = l + Linear_Expression(coeff * Variable(j))
            self._add_rational_constraint(l, denom, self.row_lower_bound[i], self.row_upper_bound[i])

        assert len(self.col_lower_bound) == len(self.col_upper_bound)
        for i in range(len(self.col_lower_bound)):
            self._add_rational_constraint(Variable(i), 1, self.col_lower_bound[i], self.col_upper_bound[i])

        if self.is_maximize == 1:
            self.mip.set_optimization_mode('maximization')
        else:
            self.mip.set_optimization_mode('minimization')

    def _add_rational_constraint(self, e, denom, lower, upper):
        """
        Helper function for adding constraints: add the constraint
        ``lower <= e/denom <= upper``.

        INPUT:

        - ``e`` -- a linear expression (type ``Linear_Expression``)

        - ``denom`` -- positive integer

        - ``lower``, ``upper`` -- a rational number or ``None``, where
          ``None`` means that there is no constraint

        TESTS:

        Create a linear system with only equalities as constraints::

            sage: p = MixedIntegerLinearProgram(solver='PPL')
            sage: x = p.new_variable(nonnegative=False)
            sage: n = 40
            sage: v = random_vector(QQ, n)
            sage: M = random_matrix(QQ, 2*n, n)
            sage: for j in range(2*n):  # indirect doctest
            ....:     lhs = p.sum(M[j,i]*x[i] for i in range(n))
            ....:     rhs = M.row(j).inner_product(v)
            ....:     p.add_constraint(lhs == rhs)
            sage: p.solve()  # long time
            0
        """
        if lower == upper:
            if lower is not None:
                rhs = Rational(lower * denom)
                self.mip.add_constraint(e * rhs.denominator() == rhs.numerator())
        else:
            if lower is not None:
                rhs = Rational(lower * denom)
                self.mip.add_constraint(e * rhs.denominator() >= rhs.numerator())
            if upper is not None:
                rhs = Rational(upper * denom)
                self.mip.add_constraint(e * rhs.denominator() <= rhs.numerator())

    cpdef int add_variable(self, lower_bound=0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

        It has not been implemented for selecting the variable type yet.

        INPUT:

        - ``lower_bound`` -- the lower bound of the variable (default: 0)

        - ``upper_bound`` -- the upper bound of the variable (default: ``None``)

        - ``binary`` -- ``True`` if the variable is binary (default: ``False``)

        - ``continuous`` -- ``True`` if the variable is continuous (default: ``True``)

        - ``integer`` -- ``True`` if the variable is integral (default: ``False``)

        - ``obj`` -- (optional) coefficient of this variable in the objective function (default: 0)

        - ``name`` -- an optional name for the newly added variable (default: ``None``)

        OUTPUT: the index of the newly created variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable(lower_bound=-2)
            1
            sage: p.add_variable(name='x',obj=2/3)
            2
            sage: p.col_name(2)
            'x'
            sage: p.objective_coefficient(2)
            2/3
            sage: p.add_variable(integer=True)
            3
        """
        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        for i in range(len(self.Matrix)):
            self.Matrix[i].append(0)
        self.col_lower_bound.append(lower_bound)
        self.col_upper_bound.append(upper_bound)
        self.objective_function.append(obj)
        self.col_name_var.append(name)

        n = len(self.objective_function) - 1
        if binary:
            self.set_variable_type(n,0)
        elif integer:
            self.set_variable_type(n,1)

        return n

    cpdef int add_variables(self, int n, lower_bound=0, upper_bound=None, binary=False, continuous=True, integer=False, obj=0, names=None) except -1:
        """
        Add ``n`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        It has not been implemented for selecting the variable type yet.

        INPUT:

        - ``n`` -- the number of new variables (must be > 0)

        - ``lower_bound`` -- the lower bound of the variable (default: 0)

        - ``upper_bound`` -- the upper bound of the variable (default: ``None``)

        - ``binary`` -- ``True`` if the variable is binary (default: ``False``)

        - ``continuous`` -- ``True`` if the variable is continuous (default: ``True``)

        - ``integer`` -- ``True`` if the variable is integral (default: ``False``)

        - ``obj`` -- coefficient of all variables in the objective function (default: 0)

        - ``names`` -- list of names (default: ``None``)

        OUTPUT: the index of the variable created last

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.ncols()
            0
            sage: p.add_variables(5)
            4
            sage: p.ncols()
            5
            sage: p.add_variables(2, lower_bound=-2.0, obj=42.0, names=['a','b'])
            6

        TESTS:

        Check that arguments are used::

            sage: p.col_bounds(5) # tol 1e-8
            (-2.0, None)
            sage: p.col_name(5)
            'a'
            sage: p.objective_coefficient(5) # tol 1e-8
            42.0
        """
        if binary or integer:
            raise NotImplementedError("The PPL backend in Sage only supports continuous variables")
        for k in range(n):
            for i in range(len(self.Matrix)):
                self.Matrix[i].append(0)
            self.col_lower_bound.append(lower_bound)
            self.col_upper_bound.append(upper_bound)
            self.objective_function.append(obj)
            if names is not None:
                self.col_name_var.append(names[k])
            else:
                self.col_name_var.append(None)
        return len(self.objective_function) - 1

    cpdef  set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable.

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``vtype`` -- integer:

            *  1  Integer
            *  0  Binary
            *  -1  Continuous

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.add_variables(5)
            4
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)
            True
            sage: p.set_variable_type(3,0)
            sage: p.is_variable_integer(3) or p.is_variable_binary(3)
            True
            sage: p.col_bounds(3) # tol 1e-6
            (0, 1)
            sage: p.set_variable_type(3, -1)
            sage: p.is_variable_continuous(3)
            True

        TESTS:

        Test that an exception is raised when an invalid type is passed::

            sage: p.set_variable_type(3, -2)
            Traceback (most recent call last):
            ...
            ValueError: ...
        """
        if vtype == -1:
            if variable in self.integer_variables:
                self.integer_variables.remove(variable)
        elif vtype == 0:
            self.integer_variables.add(variable)
            self.variable_lower_bound(variable, 0)
            self.variable_upper_bound(variable, 1)
        elif vtype == 1:
            self.integer_variables.add(variable)
        else:
            raise ValueError("Invalid variable type: {}".format(vtype))

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` -- integer:

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        if sense == 1:
            self.is_maximize = 1
        else:
            self.is_maximize = 0

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective
        function

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``coeff`` -- integer; its coefficient

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2
        """
        if coeff is not None:
            self.objective_function[variable] = coeff
        else:
            return self.objective_function[variable]

    cpdef set_objective(self, list coeff, d=0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- list of real values, whose i-th element is the
          coefficient of the i-th variable in the objective function

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='PPL')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(x[0]*5 + x[1]/11 <= 6)
            sage: p.set_objective(x[0])
            sage: p.solve()
            6/5
            sage: p.set_objective(x[0]/2 + 1)
            sage: p.show()
            Maximization:
              1/2 x_0 + 1
            <BLANKLINE>
            Constraints:
              constraint_0: 5 x_0 + 1/11 x_1 <= 6
            Variables:
              x_0 is a continuous variable (min=0, max=+oo)
              x_1 is a continuous variable (min=0, max=+oo)
            sage: p.solve()
            8/5

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1, 1, 2, 1, 3]
        """
        for i in range(len(coeff)):
            self.objective_function[i] = coeff[i]
        self.obj_constant_term = Rational(d)

    cpdef set_verbosity(self, int level):
        """
        Set the log (verbosity) level. Not Implemented.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.set_verbosity(0)
        """

    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` -- an iterable with ``(c,v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (real
          value).

        - ``lower_bound`` -- a lower bound, either a real value or ``None``

        - ``upper_bound`` -- an upper bound, either a real value or ``None``

        - ``name`` -- an optional name for this row (default: ``None``)

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='PPL')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(x[0]/2 + x[1]/3 <= 2/5)
            sage: p.set_objective(x[1])
            sage: p.solve()
            6/5
            sage: p.add_constraint(x[0] - x[1] >= 1/10)
            sage: p.solve()
            21/50
            sage: p.set_max(x[0], 1/2)
            sage: p.set_min(x[1], 3/8)
            sage: p.solve()
            2/5

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2.0, 2.0)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)
            (2.00000000000000, 2.00000000000000)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1.0, 1.0, name='foo')
            sage: p.row_name(-1)
            'foo'
        """
        last = len(self.Matrix)
        self.Matrix.append([])
        for i in range(len(self.objective_function)):
            self.Matrix[last].append(0)
        for a in coefficients:
            self.Matrix[last][a[0]] = a[1]

        self.row_lower_bound.append(lower_bound)
        self.row_upper_bound.append(upper_bound)
        self.row_name_var.append(name)

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
            sage: p = get_solver(solver = "PPL")
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(list(range(5)), list(range(5)))
            sage: p.nrows()
            5
        """
        for i in range(len(self.Matrix)):
            self.Matrix[i].append(0)
        for i in range(len(indices)):
            self.Matrix[indices[i]][len(self.Matrix[indices[i]]) - 1] = coeffs[i]

        self.col_lower_bound.append(None)
        self.col_upper_bound.append(None)
        self.objective_function.append(0)
        self.col_name_var.append(None)

    cpdef add_linear_constraints(self, int number, lower_bound, upper_bound, names=None):
        """
        Add constraints.

        INPUT:

        - ``number`` -- integer; the number of constraints to add

        - ``lower_bound`` -- a lower bound, either a real value or ``None``

        - ``upper_bound`` -- an upper bound, either a real value or ``None``

        - ``names`` -- an optional list of names (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraints(5, None, 2)
            sage: p.row(4)
            ([], [])
            sage: p.row_bounds(4)
            (None, 2)
        """
        for i in range(number):
            self.Matrix.append([])
            for j in range(len(self.objective_function)):
                self.Matrix[i].append(0)
            self.row_lower_bound.append(lower_bound)
            self.row_upper_bound.append(upper_bound)
            if names is not None:
                self.row_name_var.append(names[i])
            else:
                self.row_name_var.append(None)

    cpdef int solve(self) except -1:
        # integer example copied from cplex_backend.pyx
        """
        Solve the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution cannot be computed for any reason (none
            exists, or the solver was not able to find it, etc...)

        EXAMPLES:

        A linear optimization problem::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(list(range(5)), list(range(5)))
            sage: p.solve()
            0

        An unbounded problem::

            sage: p.objective_coefficient(0,1)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: ...

        An integer optimization problem::

            sage: p = MixedIntegerLinearProgram(solver='PPL')
            sage: x = p.new_variable(integer=True, nonnegative=True)
            sage: p.add_constraint(2*x[0] + 3*x[1], max = 6)
            sage: p.add_constraint(3*x[0] + 2*x[1], max = 6)
            sage: p.set_objective(x[0] + x[1] + 7)
            sage: p.solve()
            9
        """
        self.init_mip()

        ans = self.mip.solve()

        if ans['status'] == 'optimized':
            pass
        elif ans['status'] == 'unbounded':
            raise MIPSolverException("PPL : Solution is unbounded")
        elif ans['status'] == 'unfeasible':
            raise MIPSolverException("PPL : There is no feasible solution")

        return 0

    cpdef get_objective_value(self):
        """
        Return the exact value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='PPL')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(5/13*x[0] + x[1]/2 == 8/7)
            sage: p.set_objective(5/13*x[0] + x[1]/2)
            sage: p.solve()
            8/7

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
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
        ans = Rational(self.mip.optimal_value())
        return ans / self.obj_denominator + self.obj_constant_term

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
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
        g = self.mip.optimizing_point()
        return Integer(g.coefficient(Variable(variable))) / Integer(g.divisor())

    cpdef int ncols(self) noexcept:
        """
        Return the number of columns/variables.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.ncols()
            0
            sage: p.add_variables(2)
            1
            sage: p.ncols()
            2
        """
        return len(self.objective_function)

    cpdef int nrows(self) noexcept:
        """
        Return the number of rows/constraints.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(2, 2.0, None)
            sage: p.nrows()
            2
        """
        return len(self.Matrix)

    cpdef bint is_maximization(self) noexcept:
        """
        Test whether the problem is a maximization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        if self.is_maximize == 1:
            return 1
        else:
            return 0

    cpdef problem_name(self, name=None):
        """
        Return or define the problem's name.

        INPUT:

        - ``name`` -- string; the problem's name. When set to
          ``None`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.problem_name("There once was a french fry")
            sage: print(p.problem_name())
            There once was a french fry
        """
        if name is None:
            return self.name
        self.name = name

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
            sage: p = get_solver(solver = "PPL")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)
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
            sage: p = get_solver(solver = "PPL")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)
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
            sage: p = get_solver(solver = "PPL")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0, 5)
        """
        return (self.col_lower_bound[index], self.col_upper_bound[index])

    cpdef bint is_variable_binary(self, int index) noexcept:
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_binary(0)
            False
        """
        return index in self.integer_variables and self.col_bounds(index) == (0, 1)

    cpdef bint is_variable_integer(self, int index) noexcept:
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_integer(0)
            False
        """
        return index in self.integer_variables and self.col_bounds(index) != (0, 1)

    cpdef bint is_variable_continuous(self, int index) noexcept:
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_continuous(0)
            True
        """
        return index not in self.integer_variables

    cpdef row_name(self, int index):
        """
        Return the ``index``-th row name.

        INPUT:

        - ``index`` -- integer; the row's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.add_linear_constraints(1, 2, None, names=["Empty constraint 1"])
            sage: p.row_name(0)
            'Empty constraint 1'
        """
        if self.row_name_var[index] is not None:
            return self.row_name_var[index]
        return "constraint_" + repr(index)

    cpdef col_name(self, int index):
        """
        Return the ``index``-th col name.

        INPUT:

        - ``index`` -- integer; the col's id

        - ``name`` -- (``char *``) its name; when set to ``NULL``
          (default), the method returns the current name

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.add_variable(name="I am a variable")
            0
            sage: p.col_name(0)
            'I am a variable'
        """
        if self.col_name_var[index] is not None:
            return self.col_name_var[index]
        return "x_" + repr(index)

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
            sage: p = get_solver(solver = "PPL")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0, 5)
            sage: p.variable_upper_bound(0, None)
            sage: p.col_bounds(0)
            (0, None)
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
          variable has not lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "PPL")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0, None)
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
