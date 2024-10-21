# sage_setup: distribution = sagemath-glpk
"""
GLPK Backend

AUTHORS:

- Nathann Cohen (2010-10): initial implementation
- John Perry (2012-01): glp_simplex preprocessing
- John Perry and Raniere Gaia Silva (2012-03): solver parameters
- Christian Kuper (2012-10): Additions for sensitivity analysis
"""

# ****************************************************************************
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.float cimport DBL_MAX
from libc.limits cimport INT_MAX
from cysignals.memory cimport sig_malloc, sig_free
from cysignals.signals cimport sig_on, sig_off
from memory_allocator cimport MemoryAllocator

from sage.cpython.string cimport char_to_str, str_to_bytes
from sage.cpython.string import FS_ENCODING
from sage.numerical.mip import MIPSolverException
from sage.libs.glpk.constants cimport *
from sage.libs.glpk.lp cimport *


cdef class GLPKBackend(GenericBackend):
    """
    MIP Backend that uses the GLPK solver.
    """

    def __cinit__(self, maximization = True):
        """
        Constructor.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='GLPK')
        """
        self.lp = glp_create_prob()
        self.simplex_or_intopt = glp_simplex_then_intopt
        self.smcp = <glp_smcp* > sig_malloc(sizeof(glp_smcp))
        glp_init_smcp(self.smcp)
        self.iocp = <glp_iocp* > sig_malloc(sizeof(glp_iocp))
        glp_init_iocp(self.iocp)

        self.iocp.cb_func = glp_callback                      # callback function
        self.iocp.cb_info = <void *> &(self.search_tree_data) # callback data
        self.iocp.presolve = GLP_ON
        self.set_verbosity(0)
        self.obj_constant_term = 0.0

        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

    cpdef int add_variable(self, lower_bound=0.0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0.0, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive, real and the coefficient in the
        objective function is 0.0.

        INPUT:

        - ``lower_bound`` -- the lower bound of the variable (default: 0)

        - ``upper_bound`` -- the upper bound of the variable (default: ``None``)

        - ``binary`` -- ``True`` if the variable is binary (default: ``False``)

        - ``continuous`` -- ``True`` if the variable is continuous (default: ``True``)

        - ``integer`` -- ``True`` if the variable is integral (default: ``False``)

        - ``obj`` -- (optional) coefficient of this variable in the objective function (default: 0.0)

        - ``name`` -- an optional name for the newly added variable (default: ``None``)

        OUTPUT: the index of the newly created variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
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
        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        glp_add_cols(self.lp, 1)
        cdef int n_var = glp_get_num_cols(self.lp)

        self.variable_lower_bound(n_var - 1, lower_bound)
        self.variable_upper_bound(n_var - 1, upper_bound)

        if continuous:
            glp_set_col_kind(self.lp, n_var, GLP_CV)
        elif binary:
            glp_set_col_kind(self.lp, n_var, GLP_BV)
        elif integer:
            glp_set_col_kind(self.lp, n_var, GLP_IV)

        if name is not None:
            glp_set_col_name(self.lp, n_var, str_to_bytes(name))

        if obj:
            self.objective_coefficient(n_var - 1, obj)

        return n_var - 1

    cpdef int add_variables(self, int number, lower_bound=0.0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0.0, names=None) except -1:
        """
        Add ``number`` new variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive, real and their coefficient in
        the objective function is 0.0.

        INPUT:

        - ``n`` -- the number of new variables (must be > 0)

        - ``lower_bound`` -- the lower bound of the variable (default: 0)

        - ``upper_bound`` -- the upper bound of the variable (default: ``None``)

        - ``binary`` -- ``True`` if the variable is binary (default: ``False``)

        - ``continuous`` -- ``True`` if the variable is binary (default: ``True``)

        - ``integer`` -- ``True`` if the variable is binary (default: ``False``)

        - ``obj`` -- coefficient of all variables in the objective function (default: 0.0)

        - ``names`` -- list of names (default: ``None``)

        OUTPUT: the index of the variable created last

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variables(5)
            4
            sage: p.ncols()
            5
            sage: p.add_variables(2, lower_bound=-2.0, integer=True, obj=42.0, names=['a','b'])
            6

        TESTS:

        Check that arguments are used::

            sage: p.col_bounds(5) # tol 1e-8
            (-2.0, None)
            sage: p.is_variable_integer(5)
            True
            sage: p.col_name(5)
            'a'
            sage: p.objective_coefficient(5)
            42.0
        """
        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        glp_add_cols(self.lp, number)

        cdef int n_var
        n_var = glp_get_num_cols(self.lp)

        cdef int i

        for 0<= i < number:
            self.variable_lower_bound(n_var - i - 1, lower_bound)
            self.variable_upper_bound(n_var - i - 1, upper_bound)
            if continuous:
                glp_set_col_kind(self.lp, n_var - i, GLP_CV)
            elif binary:
                glp_set_col_kind(self.lp, n_var - i, GLP_BV)
            elif integer:
                glp_set_col_kind(self.lp, n_var - i, GLP_IV)

            if obj:
                self.objective_coefficient(n_var - i - 1, obj)

            if names is not None:
                glp_set_col_name(self.lp, n_var - i,
                                 str_to_bytes(names[number - i - 1]))

        return n_var - 1

    cpdef set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable.

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``vtype`` -- integer:

            *  1  Integer
            *  0  Binary
            * -1 Real

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)
            True

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.set_variable_type(2,0)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 2
        """
        if variable < 0 or variable > (self.ncols() - 1):
            raise ValueError("invalid variable index %d" % variable)

        if vtype==1:
            glp_set_col_kind(self.lp, variable+1, GLP_IV)

        elif vtype==0:
            glp_set_col_kind(self.lp, variable+1, GLP_BV)

        else:
            glp_set_col_kind(self.lp, variable+1, GLP_CV)

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` -- integer:

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        if sense == 1:
            glp_set_obj_dir(self.lp, GLP_MAX)
        else:
            glp_set_obj_dir(self.lp, GLP_MIN)

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective function.

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``coeff`` -- double; its coefficient or ``None`` for
          reading (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0.0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2.0

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.objective_coefficient(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 2
        """
        if variable < 0 or variable > (self.ncols() - 1):
            raise ValueError("invalid variable index %d" % variable)

        if coeff is None:
            return glp_get_obj_coef(self.lp, variable + 1)
        else:
            glp_set_obj_coef(self.lp, variable + 1, coeff)

    cpdef problem_name(self, name=None):
        """
        Return or define the problem's name.

        INPUT:

        - ``name`` -- string; the problem's name. When set to
          ``None`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.problem_name("There once was a french fry")
            sage: print(p.problem_name())
            There once was a french fry
        """
        cdef char * n

        if name is None:
            n =  <char *> glp_get_prob_name(self.lp)
            if n == NULL:
                return ""
            else:
                return char_to_str(n)

        else:
            name = str_to_bytes(name)
            if len(name) > 255:
                raise ValueError("Problem name for GLPK must not be longer than 255 characters.")
            glp_set_prob_name(self.lp, name)

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
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1.0, 1.0, 2.0, 1.0, 3.0]
        """
        cdef int i

        for i,v in enumerate(coeff):
            glp_set_obj_coef(self.lp, i+1, v)

        glp_set_obj_coef(self.lp, 0, d)

        self.obj_constant_term = d

    cpdef set_verbosity(self, int level):
        """
        Set the verbosity level.

        INPUT:

        - ``level`` -- integer; from 0 (no verbosity) to 3

        EXAMPLES::

            sage: p.<x> = MixedIntegerLinearProgram(solver='GLPK')
            sage: p.add_constraint(10 * x[0] <= 1)
            sage: p.add_constraint(5 * x[1] <= 1)
            sage: p.set_objective(x[0] + x[1])
            sage: p.solve()
            0.30000000000000004
            sage: p.get_backend().set_verbosity(3)
            sage: p.solver_parameter("simplex_or_intopt", "intopt_only")
            sage: p.solve()
            GLPK Integer Optimizer...
            2 rows, 2 columns, 2 non-zeros
            0 integer variables, none of which are binary
            Preprocessing...
            Objective value =   3.000000000e-01
            INTEGER OPTIMAL SOLUTION FOUND BY MIP PREPROCESSOR
            0.30000000000000004

        ::

            sage: p.<x> = MixedIntegerLinearProgram(solver='GLPK/exact')
            sage: p.add_constraint(10 * x[0] <= 1)
            sage: p.add_constraint(5 * x[1] <= 1)
            sage: p.set_objective(x[0] + x[1])
            sage: p.solve() # tol 1e-14
            0.3
            sage: p.get_backend().set_verbosity(2)
            sage: p.solve() # tol 1e-14
            *     2:   objval =                    0.3   (0)
            *     2:   objval =                    0.3   (0)
            0.3
            sage: p.get_backend().set_verbosity(3)
            sage: p.solve() # tol 1e-14
            glp_exact: 2 rows, 2 columns, 2 non-zeros
            ...
            *     2:   objval =                    0.3   (0)
            *     2:   objval =                    0.3   (0)
            OPTIMAL SOLUTION FOUND
            0.3
        """
        if level == 0:
            self.iocp.msg_lev = GLP_MSG_OFF
            self.smcp.msg_lev = GLP_MSG_OFF
        elif level == 1:
            self.iocp.msg_lev = GLP_MSG_ERR
            self.smcp.msg_lev = GLP_MSG_ERR
        elif level == 2:
            self.iocp.msg_lev = GLP_MSG_ON
            self.smcp.msg_lev = GLP_MSG_ON
        else:
            self.iocp.msg_lev = GLP_MSG_ALL
            self.smcp.msg_lev = GLP_MSG_ALL

    cpdef remove_constraint(self, int i):
        r"""
        Remove a constraint from ``self``.

        INPUT:

        - ``i`` -- index of the constraint to remove

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='GLPK')
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

            sage: MixedIntegerLinearProgram(solver = "GLPK").remove_constraint(-2)
            Traceback (most recent call last):
            ...
            ValueError: The constraint's index i must satisfy 0 <= i < number_of_constraints
        """
        cdef int rows[2]

        if i < 0 or i >= glp_get_num_rows(self.lp):
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")

        rows[1] = i + 1
        glp_del_rows(self.lp, 1, rows)
        glp_std_basis(self.lp)

    cpdef remove_constraints(self, constraints):
        r"""
        Remove several constraints.

        INPUT:

        - ``constraints`` -- an iterable containing the indices of the rows to remove

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='GLPK')
            sage: x, y = p['x'], p['y']
            sage: p.add_constraint(2*x + 3*y <= 6)
            sage: p.add_constraint(3*x + 2*y <= 6)
            sage: p.add_constraint(x >= 0)
            sage: p.set_objective(x + y + 7)
            sage: p.set_integer(x); p.set_integer(y)
            sage: p.solve()
            9.0
            sage: p.remove_constraints([0])
            sage: p.solve()
            10.0
            sage: p.get_values([x,y])
            [0.0, 3.0]

        TESTS:

        Removing fancy constraints does not make Sage crash::

            sage: MixedIntegerLinearProgram(solver= "GLPK").remove_constraints([0, -2])
            Traceback (most recent call last):
            ...
            ValueError: The constraint's index i must satisfy 0 <= i < number_of_constraints
        """
        cdef int i, c
        cdef int m = len(constraints)
        cdef int * rows = <int *>sig_malloc((m + 1) * sizeof(int *))
        cdef int nrows = glp_get_num_rows(self.lp)

        for i in range(m):

            c = constraints[i]
            if c < 0 or c >= nrows:
                sig_free(rows)
                raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")

            rows[i+1] = c + 1

        glp_del_rows(self.lp, m, rows)
        sig_free(rows)
        glp_std_basis(self.lp)

    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` an iterable with ``(c,v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (real
          value).

        - ``lower_bound`` -- a lower bound, either a real value or ``None``

        - ``upper_bound`` -- an upper bound, either a real value or ``None``

        - ``name`` -- an optional name for this row (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2.0, 2.0)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)
            sage: p.add_linear_constraint(zip(range(5), range(5)), 1.0, 1.0, name='foo')
            sage: p.row_name(1)
            'foo'

        TESTS:

        This used to crash Sage, but was fixed in :issue:`19525`::

            sage: p = MixedIntegerLinearProgram(solver='glpk')
            sage: q = MixedIntegerLinearProgram(solver='glpk')
            sage: q.add_constraint(p.new_variable()[0] <= 1)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 0
        """
        if lower_bound is None and upper_bound is None:
            raise ValueError("At least one of 'upper_bound' or 'lower_bound' must be set.")

        # We're going to iterate through this more than once.
        coefficients = list(coefficients)

        for (index, _) in coefficients:
            if index < 0 or index > (self.ncols() - 1):
                raise ValueError("invalid variable index %d" % index)

        glp_add_rows(self.lp, 1)
        cdef int n = glp_get_num_rows(self.lp)

        cdef MemoryAllocator mem = MemoryAllocator()
        cdef int * row_i
        cdef double * row_values

        cdef int n_coeff = len(coefficients)
        row_i = <int*>mem.allocarray(n_coeff + 1, sizeof(int))
        row_values = <double*>mem.allocarray(n_coeff + 1, sizeof(double))

        cdef Py_ssize_t i = 1
        for c,v in coefficients:
            row_i[i] = c+1
            row_values[i] = v
            i += 1

        sig_on()
        glp_set_mat_row(self.lp, n, n_coeff, row_i, row_values)
        sig_off()

        if upper_bound is not None and lower_bound is None:
            glp_set_row_bnds(self.lp, n, GLP_UP, upper_bound, upper_bound)
        elif lower_bound is not None and upper_bound is None:
            glp_set_row_bnds(self.lp, n, GLP_LO, lower_bound, lower_bound)
        elif upper_bound is not None and lower_bound is not None:
            if lower_bound == upper_bound:
                glp_set_row_bnds(self.lp, n, GLP_FX, lower_bound, upper_bound)
            else:
                glp_set_row_bnds(self.lp, n, GLP_DB, lower_bound, upper_bound)

        if name is not None:
            glp_set_row_name(self.lp, n, str_to_bytes(name))

    cpdef add_linear_constraints(self, int number, lower_bound, upper_bound, names=None):
        """
        Add ``'number`` linear constraints.

        INPUT:

        - ``number`` -- integer; the number of constraints to add

        - ``lower_bound`` -- a lower bound, either a real value or ``None``

        - ``upper_bound`` -- an upper bound, either a real value or ``None``

        - ``names`` -- an optional list of names (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraints(5, None, 2)
            sage: p.row(4)
            ([], [])
            sage: p.row_bounds(4)
            (None, 2.0)
            sage: p.add_linear_constraints(2, None, 2, names=['foo','bar'])
        """
        if lower_bound is None and upper_bound is None:
            raise ValueError("At least one of 'upper_bound' or 'lower_bound' must be set.")

        glp_add_rows(self.lp, number)
        cdef int n = glp_get_num_rows(self.lp)

        cdef int i
        for 0<= i < number:
            if upper_bound is not None and lower_bound is None:
                glp_set_row_bnds(self.lp, n-i, GLP_UP, upper_bound, upper_bound)
            elif lower_bound is not None and upper_bound is None:
                glp_set_row_bnds(self.lp, n-i, GLP_LO, lower_bound, lower_bound)
            elif upper_bound is not None and lower_bound is not None:
                if lower_bound == upper_bound:
                    glp_set_row_bnds(self.lp, n-i, GLP_FX, lower_bound, upper_bound)
                else:
                    glp_set_row_bnds(self.lp, n-i, GLP_DB, lower_bound, upper_bound)
            if names is not None:
                glp_set_row_name(self.lp, n-i,
                                 str_to_bytes(names[number-i-1]))

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
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.row(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid row index 2
        """
        if index < 0 or index > (self.nrows() - 1):
            raise ValueError("invalid row index %d" % index)

        cdef int n = glp_get_num_cols(self.lp)
        cdef MemoryAllocator mem = MemoryAllocator()
        cdef int * c_indices = <int*>mem.allocarray(n+1, sizeof(int))
        cdef double * c_values = <double*>mem.allocarray(n+1, sizeof(double))
        cdef list indices = []
        cdef list values = []
        cdef int i,j

        i = glp_get_mat_row(self.lp, index + 1, c_indices, c_values)
        for 0 < j <= i:
            indices.append(c_indices[j]-1)
            values.append(c_values[j])

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
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.row_bounds(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid row index 2
        """
        cdef double ub
        cdef double lb

        if index < 0 or index > (self.nrows() - 1):
            raise ValueError("invalid row index %d" % index)

        ub = glp_get_row_ub(self.lp, index + 1)
        lb = glp_get_row_lb(self.lp, index +1)

        return (
            (lb if lb != -DBL_MAX else None),
            (ub if ub != +DBL_MAX else None)
            )

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
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0.0, 5.0)

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.col_bounds(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid column index 2
        """

        cdef double ub
        cdef double lb

        if index < 0 or index > (self.ncols() - 1):
            raise ValueError("invalid column index %d" % index)

        ub = glp_get_col_ub(self.lp, index +1)
        lb = glp_get_col_lb(self.lp, index +1)

        return (
            (lb if lb != -DBL_MAX else None),
            (ub if ub != +DBL_MAX else None)
            )

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
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5
        """

        glp_add_cols(self.lp, 1)
        cdef int n = glp_get_num_cols(self.lp)

        cdef int * col_i
        cdef double * col_values

        col_i = <int *> sig_malloc((len(indices)+1) * sizeof(int))
        col_values = <double *> sig_malloc((len(indices)+1) * sizeof(double))

        for i,v in enumerate(indices):
            col_i[i+1] = v+1
        for i,v in enumerate(coeffs):
            col_values[i+1] = v

        glp_set_mat_col(self.lp, n, len(indices), col_i, col_values)
        glp_set_col_bnds(self.lp, n, GLP_LO, 0,0)
        sig_free(col_i)
        sig_free(col_values)

    cpdef int solve(self) except -1:
        """
        Solve the problem.

        Sage uses GLPK's implementation of the branch-and-cut
        algorithm (``glp_intopt``) to solve the mixed-integer linear
        program.  This algorithm can be requested explicitly by
        setting the solver parameter "simplex_or_intopt" to
        "intopt_only". By default, the simplex method will be used
        first to detect pathological problems that the integer solver
        cannot handle. If all variables are continuous, the integer
        algorithm reduces to solving the linear program by the simplex
        method.

        EXAMPLES::

            sage: lp = MixedIntegerLinearProgram(solver = 'GLPK', maximization = False)
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
            the solution cannot be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.solve()
            0
            sage: p.objective_coefficient(0,1)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: ...

        .. WARNING::

            GLPK's ``glp_intopt`` sometimes fails catastrophically
            when given a system it cannot solve (:issue:`12309`). It
            can loop indefinitely, or just plain segfault. Upstream
            considers this behavior "essentially innate" to the
            current design, and suggests preprocessing with
            ``glp_simplex``, which is what SageMath does by default.
            Set the ``simplex_or_intopt`` solver parameter to
            ``glp_intopt_only`` at your own risk.

        EXAMPLES::

            sage: lp = MixedIntegerLinearProgram(solver = "GLPK")
            sage: v = lp.new_variable(nonnegative=True)
            sage: lp.add_constraint(v[1] +v[2] -2.0 *v[3], max=-1.0)
            sage: lp.add_constraint(v[0] -4.0/3 *v[1] +1.0/3 *v[2], max=-1.0/3)
            sage: lp.add_constraint(v[0] +0.5 *v[1] -0.5 *v[2] +0.25 *v[3], max=-0.25)
            sage: lp.solve()
            0.0
            sage: lp.add_constraint(v[0] +4.0 *v[1] -v[2] +v[3], max=-1.0)
            sage: lp.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: GLPK: Problem has no feasible solution

        If we switch to "simplex_only", the integrality constraints are ignored,
        and we get an optimal solution to the continuous relaxation.

        EXAMPLES::

            sage: lp = MixedIntegerLinearProgram(solver = 'GLPK', maximization = False)
            sage: x, y = lp[0], lp[1]
            sage: lp.add_constraint(-2*x + y <= 1)
            sage: lp.add_constraint(x - y <= 1)
            sage: lp.add_constraint(x + y >= 2)
            sage: lp.set_objective(x + y)
            sage: lp.set_integer(x)
            sage: lp.set_integer(y)
            sage: lp.solver_parameter("simplex_or_intopt", "simplex_only") # use simplex only
            sage: lp.solve()
            2.0
            sage: lp.get_values([x, y])
            [1.5, 0.5]

        If one solves a linear program and wishes to access dual information
        (`get_col_dual` etc.) or tableau data (`get_row_stat` etc.),
        one needs to switch to "simplex_only" before solving.

        GLPK also has an exact rational simplex solver.  The only access
        to data is via double-precision floats, which means that
        rationals in the input data may be rounded before the exact
        solver sees them. Thus, it is unreasonable to expect that
        arbitrary LPs with rational coefficients are solved exactly.
        Once the LP has been read into the backend, it reconstructs
        rationals from doubles and does solve exactly over the rationals,
        but results are returned as as doubles.

        EXAMPLES::

            sage: lp.solver_parameter("simplex_or_intopt", "exact_simplex_only") # use exact simplex only
            sage: lp.solve()
            2.0
            sage: lp.get_values([x, y])
            [1.5, 0.5]

        If you need the rational solution, you need to retrieve the
        basis information via ``get_col_stat`` and ``get_row_stat``
        and calculate the corresponding basic solution.  Below we only
        test that the basis information is indeed available.
        Calculating the corresponding basic solution is left as an
        exercise.

        EXAMPLES::

            sage: lp.get_backend().get_row_stat(0)
            1
            sage: lp.get_backend().get_col_stat(0)
            1

        Below we test that integers that can be exactly represented by
        IEEE 754 double-precision floating point numbers survive the
        rational reconstruction done by ``glp_exact`` and the subsequent
        conversion to double-precision floating point numbers.

        EXAMPLES::

            sage: lp = MixedIntegerLinearProgram(solver = 'GLPK', maximization = True)
            sage: test = 2^53 - 43
            sage: lp.solver_parameter("simplex_or_intopt", "exact_simplex_only") # use exact simplex only
            sage: x = lp[0]
            sage: lp.add_constraint(x <= test)
            sage: lp.set_objective(x)
            sage: lp.solve() == test # yes, we want an exact comparison here
            True
            sage: lp.get_values(x) == test # yes, we want an exact comparison here
            True

        Below we test that GLPK backend can detect unboundedness in
        "simplex_only" mode (:issue:`18838`).

        EXAMPLES::

            sage: lp = MixedIntegerLinearProgram(maximization=True, solver = "GLPK")
            sage: lp.set_objective(lp[0])
            sage: lp.solver_parameter("simplex_or_intopt", "simplex_only")
            sage: lp.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: GLPK: Problem has unbounded solution
            sage: lp.set_objective(lp[1])
            sage: lp.solver_parameter("primal_v_dual", "GLP_DUAL")
            sage: lp.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: GLPK: Problem has unbounded solution
            sage: lp.solver_parameter("simplex_or_intopt", "simplex_then_intopt")
            sage: lp.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: GLPK: The LP (relaxation) problem has no dual feasible solution
            sage: lp.solver_parameter("simplex_or_intopt", "intopt_only")
            sage: lp.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: GLPK: The LP (relaxation) problem has no dual feasible solution
            sage: lp.set_max(lp[1],5)
            sage: lp.solve()
            5.0

        Solving a LP within the acceptable gap. No exception is raised, even if
        the result is not optimal. To do this, we try to compute the maximum
        number of disjoint balls (of diameter 1) in a hypercube::

            sage: # needs sage.graphs
            sage: g = graphs.CubeGraph(9)
            sage: p = MixedIntegerLinearProgram(solver='GLPK')
            sage: p.solver_parameter("mip_gap_tolerance",100)
            sage: b = p.new_variable(binary=True)
            sage: p.set_objective(p.sum(b[v] for v in g))
            sage: for v in g:
            ....:     p.add_constraint(b[v]+p.sum(b[u] for u in g.neighbors(v)) <= 1)
            sage: p.add_constraint(b[v] == 1) # Force an easy non-0 solution
            sage: p.solve() # rel tol 100
            1

        Same, now with a time limit::

            sage: # needs sage.graphs
            sage: p.solver_parameter("mip_gap_tolerance",1)
            sage: p.solver_parameter("timelimit",3.0)
            sage: p.solve() # rel tol 100
            1
        """

        cdef int solve_status
        cdef int solution_status
        global solve_status_msg
        global solution_status_msg

        if (self.simplex_or_intopt == glp_simplex_only
            or self.simplex_or_intopt == glp_simplex_then_intopt
            or self.simplex_or_intopt == glp_exact_simplex_only):
            if self.simplex_or_intopt == glp_exact_simplex_only:
                solve_status = glp_exact(self.lp, self.smcp)
            else:
                solve_status = glp_simplex(self.lp, self.smcp)
            solution_status = glp_get_status(self.lp)

        if ((self.simplex_or_intopt == glp_intopt_only)
            or (self.simplex_or_intopt == glp_simplex_then_intopt) and (solution_status != GLP_UNDEF) and (solution_status != GLP_NOFEAS)):
            sig_on()
            solve_status = glp_intopt(self.lp, self.iocp)
            solution_status = glp_mip_status(self.lp)
            sig_off()

        if solution_status == GLP_OPT:
            pass
        elif (solution_status == GLP_FEAS) and (solve_status == GLP_ETMLIM or solve_status == GLP_EITLIM
                                                or solve_status == GLP_EMIPGAP or solve_status == GLP_EOBJLL or solve_status == GLP_EOBJUL):
            # no exception when time limit or iteration limit or  mip gap tolerances or objective limits reached.
            pass
        elif solution_status == GLP_UNDEF:
            raise MIPSolverException("GLPK: "+solve_status_msg.get(solve_status, "unknown error during call to GLPK : "+str(solve_status)))
        else:
            raise MIPSolverException("GLPK: "+solution_status_msg.get(solution_status, "unknown error during call to GLPK : "+str(solution_status)))
        return 0

    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
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
        if (self.simplex_or_intopt != glp_simplex_only
            and self.simplex_or_intopt != glp_exact_simplex_only):
            return glp_mip_obj_val(self.lp)
        else:
            return glp_get_obj_val(self.lp)

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
            sage: p = MixedIntegerLinearProgram(solver='GLPK')
            sage: p.solver_parameter("mip_gap_tolerance",100)
            sage: b = p.new_variable(binary=True)
            sage: p.set_objective(p.sum(b[v] for v in g))
            sage: for v in g:
            ....:     p.add_constraint(b[v]+p.sum(b[u] for u in g.neighbors(v)) <= 1)
            sage: p.add_constraint(b[v] == 1) # Force an easy non-0 solution
            sage: p.solve() # rel tol 100
            1.0
            sage: backend = p.get_backend()
            sage: backend.best_known_objective_bound() # random
            48.0
        """
        return self.search_tree_data.best_bound

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
            sage: p = MixedIntegerLinearProgram(solver='GLPK')
            sage: p.solver_parameter("mip_gap_tolerance",100)
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
        return self.search_tree_data.mip_gap

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
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

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.get_variable_value(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 2
        """
        if variable < 0 or variable > (self.ncols() - 1):
            raise ValueError("invalid variable index %d" % variable)

        if (self.simplex_or_intopt != glp_simplex_only
            and self.simplex_or_intopt != glp_exact_simplex_only):
            return glp_mip_col_val(self.lp, variable + 1)
        else:
            return glp_get_col_prim(self.lp, variable + 1)

    cpdef get_row_prim(self, int i):
        r"""
        Return the value of the auxiliary variable associated with i-th row.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "GLPK")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
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

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.get_row_prim(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid row index 2
        """
        if i < 0 or i > (self.nrows() - 1):
            raise ValueError("invalid row index %d" % i)

        return glp_get_row_prim(self.lp, i+1)

    cpdef int ncols(self) noexcept:
        """
        Return the number of columns/variables.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variables(2)
            1
            sage: p.ncols()
            2
        """
        return glp_get_num_cols(self.lp)

    cpdef int nrows(self) noexcept:
        """
        Return the number of rows/constraints.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(2, 2, None)
            sage: p.nrows()
            2
        """

        return glp_get_num_rows(self.lp)

    cpdef col_name(self, int index):
        """
        Return the ``index``-th col name.

        INPUT:

        - ``index`` -- integer; the col's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable(name='I am a variable')
            0
            sage: p.col_name(0)
            'I am a variable'

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.col_name(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid column index 2
        """
        cdef char * s

        if index < 0 or index > (self.ncols() - 1):
            raise ValueError("invalid column index %d" % index)

        glp_create_index(self.lp)
        s = <char*> glp_get_col_name(self.lp, index + 1)

        if s != NULL:
            return char_to_str(s)
        else:
            return ""

    cpdef row_name(self, int index):
        """
        Return the ``index``-th row name.

        INPUT:

        - ``index`` -- integer; the row's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_linear_constraints(1, 2, None, names=['Empty constraint 1'])
            sage: p.row_name(0)
            'Empty constraint 1'

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.row_name(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid row index 2
        """
        cdef char *  s

        if index < 0 or index > (self.nrows() - 1):
            raise ValueError("invalid row index %d" % index)

        glp_create_index(self.lp)
        s = <char*> glp_get_row_name(self.lp, index + 1)

        if s != NULL:
            return char_to_str(s)
        else:
            return ""

    cpdef bint is_variable_binary(self, int index) noexcept:
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,0)
            sage: p.is_variable_binary(0)
            True

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.is_variable_binary(2)
            False
        """
        if index < 0 or index > (self.ncols() - 1):
            # This is how the other backends behave, and this method is
            # unable to raise a python exception as currently defined.
            return False

        return glp_get_col_kind(self.lp, index + 1) == GLP_BV

    cpdef bint is_variable_integer(self, int index) noexcept:
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)
            True

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.is_variable_integer(2)
            False
        """
        if index < 0 or index > (self.ncols() - 1):
            # This is how the other backends behave, and this method is
            # unable to raise a python exception as currently defined.
            return False

        return glp_get_col_kind(self.lp, index + 1) == GLP_IV

    cpdef bint is_variable_continuous(self, int index) noexcept:
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_continuous(0)
            True
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_continuous(0)
            False

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.is_variable_continuous(2)
            False
        """
        if index < 0 or index > (self.ncols() - 1):
            # This is how the other backends behave, and this method is
            # unable to raise a python exception as currently defined.
            return False

        return glp_get_col_kind(self.lp, index + 1) == GLP_CV

    cpdef bint is_maximization(self) noexcept:
        """
        Test whether the problem is a maximization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """

        return glp_get_obj_dir(self.lp) == GLP_MAX

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
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0.0, 5.0)

        TESTS:

        :issue:`14581`::

            sage: P = MixedIntegerLinearProgram(solver='GLPK')
            sage: x = P["x"]
            sage: P.set_max(x, 0)
            sage: P.get_max(x)
            0.0

        Check that :issue:`10232` is fixed::

            sage: p = get_solver(solver='GLPK')
            sage: p.variable_upper_bound(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 2

            sage: p.variable_upper_bound(-1)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index -1

            sage: p.variable_upper_bound(3, 5)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 3

            sage: p.add_variable()
            0
            sage: p.variable_upper_bound(0, 'hey!')
            Traceback (most recent call last):
            ...
            TypeError: must be real number, not str
        """
        cdef double x
        cdef double min
        cdef double dvalue

        if index < 0 or index > (self.ncols() - 1):
            raise ValueError("invalid variable index %d" % index)

        if value is False:
            sig_on()
            x = glp_get_col_ub(self.lp, index +1)
            sig_off()
            if x == DBL_MAX:
                return None
            else:
                return x
        else:
            sig_on()
            min = glp_get_col_lb(self.lp, index + 1)
            sig_off()

            if value is None:
                sig_on()
                if min == -DBL_MAX:
                    glp_set_col_bnds(self.lp, index + 1, GLP_FR, 0, 0)
                else:
                    glp_set_col_bnds(self.lp, index + 1, GLP_LO, min, 0)
                sig_off()
            else:
                dvalue = <double?> value

                sig_on()
                if min == -DBL_MAX:
                    glp_set_col_bnds(self.lp, index + 1, GLP_UP, 0, dvalue)

                elif min == dvalue:
                    glp_set_col_bnds(self.lp, index + 1, GLP_FX,  dvalue, dvalue)
                else:
                    glp_set_col_bnds(self.lp, index + 1, GLP_DB, min, dvalue)
                sig_off()

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
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)
            sage: p.col_bounds(0)
            (5.0, None)

        TESTS:

        :issue:`14581`::

            sage: P = MixedIntegerLinearProgram(solver='GLPK')
            sage: x = P["x"]
            sage: P.set_min(x, 5)
            sage: P.set_min(x, 0)
            sage: P.get_min(x)
            0.0

        Check that :issue:`10232` is fixed::

            sage: p = get_solver(solver='GLPK')
            sage: p.variable_lower_bound(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 2

            sage: p.variable_lower_bound(-1)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index -1

            sage: p.variable_lower_bound(3, 5)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 3

            sage: p.add_variable()
            0
            sage: p.variable_lower_bound(0, 'hey!')
            Traceback (most recent call last):
            ...
            TypeError: must be real number, not str
        """
        cdef double x
        cdef double max
        cdef double dvalue

        if index < 0 or index > (self.ncols() - 1):
            raise ValueError("invalid variable index %d" % index)

        if value is False:
            sig_on()
            x = glp_get_col_lb(self.lp, index +1)
            sig_off()
            if x == -DBL_MAX:
                return None
            else:
                return x
        else:
            sig_on()
            max = glp_get_col_ub(self.lp, index + 1)
            sig_off()

            if value is None:
                sig_on()
                if max == DBL_MAX:
                    glp_set_col_bnds(self.lp, index + 1, GLP_FR, 0.0, 0.0)
                else:
                    glp_set_col_bnds(self.lp, index + 1, GLP_UP, 0.0, max)
                sig_off()

            else:
                dvalue = <double?> value

                sig_on()
                if max == DBL_MAX:
                    glp_set_col_bnds(self.lp, index + 1, GLP_LO, value, 0.0)
                elif max == value:
                    glp_set_col_bnds(self.lp, index + 1, GLP_FX,  value, value)
                else:
                    glp_set_col_bnds(self.lp, index + 1, GLP_DB, value, max)
                sig_off()

    cpdef write_lp(self, filename):
        """
        Write the problem to a ``.lp`` file.

        INPUT:

        - ``filename`` -- string

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
            sage: p.set_objective([2, 5])
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(suffix='.lp') as f:
            ....:     _ = p.write_lp(f.name)
            ....:     len(f.readlines())
            ...
            9 lines were written
            9
        """
        filename = str_to_bytes(filename, FS_ENCODING, 'surrogateescape')
        glp_write_lp(self.lp, NULL, filename)

    cpdef write_mps(self, filename, int modern):
        """
        Write the problem to a ``.mps`` file.

        INPUT:

        - ``filename`` -- string

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
            sage: p.set_objective([2, 5])
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(suffix='mps') as f:
            ....:     _ = p.write_mps(f.name, 2)
            ....:     len(f.readlines())
            ...
            17 records were written
            17
        """
        filename = str_to_bytes(filename, FS_ENCODING, 'surrogateescape')
        glp_write_mps(self.lp, modern, NULL, filename)

    cpdef __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "GLPK")
            sage: b = p.new_variable()
            sage: p.add_constraint(b[1] + b[2] <= 6)
            sage: p.set_objective(b[1] + b[2])
            sage: copy(p).solve()
            6.0
        """
        cdef GLPKBackend p = type(self)(maximization = (1 if self.is_maximization() else -1))
        p.simplex_or_intopt = self.simplex_or_intopt
        p.iocp.tm_lim = self.iocp.tm_lim
        glp_copy_prob(p.lp, self.lp, 1)
        return p

    cpdef solver_parameter(self, name, value=None):
        """
        Return or define a solver parameter.

        INPUT:

        - ``name`` -- string; the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value

        You can supply the name of a parameter and its value using either a
        string or a ``glp_`` constant (which are defined as Cython variables of
        this module).

        In most cases, you can use the same name for a parameter as that given
        in the GLPK documentation, which is available by downloading GLPK from
        <http://www.gnu.org/software/glpk/>. The exceptions relate to parameters
        common to both methods; these require you to append ``_simplex`` or
        ``_intopt`` to the name to resolve ambiguity, since the interface allows
        access to both.

        We have also provided more meaningful names, to assist readability.

        Parameter **names** are specified in lower case.
        To use a constant instead of a string, prepend ``glp_`` to the name.
        For example, both ``glp_gmi_cuts`` or ``'gmi_cuts'`` control whether
        to solve using Gomory cuts.

        Parameter **values** are specified as strings in upper case,
        or as constants in lower case. For example, both ``glp_on`` and ``"GLP_ON"``
        specify the same thing.

        Naturally, you can use ``True`` and ``False`` in cases where ``glp_on`` and ``glp_off``
        would be used.

        A list of parameter names, with their possible values:

        **General-purpose parameters:**

        .. list-table::
         :widths: 30 70

         * - ``timelimit``

           - specify the time limit IN SECONDS.  This affects both simplex
             and intopt.

         * - ``timelimit_simplex`` and ``timelimit_intopt``

           - specify the time limit IN MILLISECONDS. (This is glpk's
             default.)

         * - ``simplex_or_intopt``

           - specify which solution routines in GLPK to use. Set this
             to either ``simplex_only``, ``exact_simplex_only``,
             ``intopt_only``, or ``simplex_then_intopt`` (the
             default). The ``simplex_then_intopt`` option does some
             extra work, but avoids hangs/crashes in GLPK on problems
             with no solution; SageMath will try simplex first, then
             perform integer optimization only if a solution of the LP
             relaxation exists. If you know that your system is not
             pathological, one of the other options will be faster.

         * - ``verbosity_intopt`` and ``verbosity_simplex``

           - one of ``GLP_MSG_OFF``, ``GLP_MSG_ERR``, ``GLP_MSG_ON``, or
             ``GLP_MSG_ALL``. The default is ``GLP_MSG_OFF``.

         * - ``output_frequency_intopt`` and ``output_frequency_simplex``

           - the output frequency, in milliseconds. Default is 5000.

         * - ``output_delay_intopt`` and ``output_delay_simplex``

           - the output delay, in milliseconds, regarding the use of the
             simplex method on the LP relaxation. Default is 10000.


        **intopt-specific parameters:**

        .. list-table::
         :widths: 30 70

         * - ``branching``
           - - ``GLP_BR_FFV``  first fractional variable
             - ``GLP_BR_LFV``  last fractional variable
             - ``GLP_BR_MFV``  most fractional variable
             - ``GLP_BR_DTH``  Driebeck-Tomlin heuristic (default)
             - ``GLP_BR_PCH``  hybrid pseudocost heuristic

         * - ``backtracking``
           - - ``GLP_BT_DFS``  depth first search
             - ``GLP_BT_BFS``  breadth first search
             - ``GLP_BT_BLB``  best local bound (default)
             - ``GLP_BT_BPH``  best projection heuristic

         * - ``preprocessing``
           - - ``GLP_PP_NONE``
             - ``GLP_PP_ROOT``  preprocessing only at root level
             - ``GLP_PP_ALL``   (default)


         * - ``feasibility_pump``

           - ``GLP_ON`` or ``GLP_OFF`` (default)

         * - ``gomory_cuts``

           - ``GLP_ON`` or ``GLP_OFF`` (default)

         * - ``mixed_int_rounding_cuts``

           - ``GLP_ON`` or ``GLP_OFF`` (default)

         * - ``mixed_cover_cuts``

           - ``GLP_ON`` or ``GLP_OFF`` (default)

         * - ``clique_cuts``

           - ``GLP_ON`` or ``GLP_OFF`` (default)

         * - ``absolute_tolerance``

           - (double) used to check if optimal solution to LP relaxation is
             integer feasible. GLPK manual advises, "do not change... without
             detailed understanding of its purpose."

         * - ``relative_tolerance``

           - (double) used to check if objective value in LP relaxation is not
             better than best known integer solution. GLPK manual advises, "do
             not change... without detailed understanding of its purpose."

         * - ``mip_gap_tolerance``

           - (double) relative mip gap tolerance. Default is 0.0.

         * - ``presolve_intopt``

           - ``GLP_ON`` (default) or ``GLP_OFF``

         * - ``binarize``

           - ``GLP_ON`` or ``GLP_OFF`` (default)


        **simplex-specific parameters:**

        .. list-table::
         :widths: 30 70

         * - ``primal_v_dual``

           - - ``GLP_PRIMAL``  (default)
             - ``GLP_DUAL``
             - ``GLP_DUALP``

         * - ``pricing``

           - - ``GLP_PT_STD``    standard (textbook)
             - ``GLP_PT_PSE``    projected steepest edge (default)

         * - ``ratio_test``

           - - ``GLP_RT_STD``  standard (textbook)
             - ``GLP_RT_HAR``  Harris' two-pass ratio test (default)

         * - ``tolerance_primal``

           - (double) tolerance used to check if basic solution is primal
             feasible. GLPK manual advises, "do not change... without
             detailed understanding of its purpose."

         * - ``tolerance_dual``

           - (double) tolerance used to check if basic solution is dual
             feasible. GLPK manual advises, "do not change... without
             detailed understanding of its purpose."

         * - ``tolerance_pivot``

           - (double) tolerance used to choose pivot. GLPK manual advises,
             "do not change... without detailed understanding of its
             purpose."

         * - ``obj_lower_limit``

           - (double) lower limit of the objective function.  The default is
             ``-DBL_MAX``.

         * - ``obj_upper_limit``

           - (double) upper limit of the objective function.  The default is
             ``DBL_MAX``.

         * - ``iteration_limit``

           - (int) iteration limit of the simplex algorithm.  The default is
             ``INT_MAX``.

         * - ``presolve_simplex``

           - ``GLP_ON`` or ``GLP_OFF`` (default).

        .. NOTE::

            The coverage for GLPK's control parameters for simplex and integer optimization
            is nearly complete. The only thing lacking is a wrapper for callback routines.

            To date, no attempt has been made to expose the interior point methods.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.solver_parameter("timelimit", 60)
            sage: p.solver_parameter("timelimit")
            60.0

        - Don't forget the difference between ``timelimit`` and ``timelimit_intopt``

        ::

            sage: p.solver_parameter("timelimit_intopt")
            60000

        If you don't care for an integer answer, you can ask for an LP
        relaxation instead.  The default solver performs integer optimization,
        but you can switch to the standard simplex algorithm through the
        ``glp_simplex_or_intopt`` parameter.

        EXAMPLES::

            sage: lp = MixedIntegerLinearProgram(solver = 'GLPK', maximization = False)
            sage: x, y = lp[0], lp[1]
            sage: lp.add_constraint(-2*x + y <= 1)
            sage: lp.add_constraint(x - y <= 1)
            sage: lp.add_constraint(x + y >= 2)
            sage: lp.set_integer(x); lp.set_integer(y)
            sage: lp.set_objective(x + y)
            sage: lp.solve()
            2.0
            sage: lp.get_values([x,y])
            [1.0, 1.0]
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: lp.solve()
            2.0
            sage: lp.get_values([x,y])
            [1.5, 0.5]

        You can get GLPK to spout all sorts of information at you.
        The default is to turn this off, but sometimes (debugging) it's very useful::

            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_then_intopt)
            sage: lp.solver_parameter(backend.glp_mir_cuts, backend.glp_on)
            sage: lp.solver_parameter(backend.glp_msg_lev_intopt, backend.glp_msg_all)
            sage: lp.solver_parameter(backend.glp_mir_cuts)
            1

        If you actually try to solve ``lp``, you will get a lot of detailed information.
        """

        if isinstance(name, str):
            if name == "out_frq" or name == "out_dly" or name == "tm_lim" or name == "msg_lev":
                raise ValueError("To set parameter " + name + " you must specify the solver. Append either _simplex or _intopt.")
            name = solver_parameter_names_dict[name]

        if isinstance(value, str):
            value = solver_parameter_values[value]

        if name == timelimit_intopt:
            if value is None: return self.iocp.tm_lim
            else: self.iocp.tm_lim = value

        if name == timelimit_seconds:
            if value is None: return self.iocp.tm_lim / 1000.0
            else:
                self.iocp.tm_lim = value * 1000.0
                self.smcp.tm_lim = value * 1000.0

        elif name == timelimit_simplex:
            if value is None: return self.smcp.tm_lim
            else: self.smcp.tm_lim = value

        elif name == simplex_or_intopt:
            if value is None:
                return self.simplex_or_intopt
            if value not in (simplex_only, intopt_only, simplex_then_intopt, exact_simplex_only):
                raise MIPSolverException("GLPK: invalid value for simplex_or_intopt; see documentation")
            self.simplex_or_intopt = value

        elif name == msg_lev_simplex:
            if value is None: return self.smcp.msg_lev
            else: self.smcp.msg_lev = value

        elif name == msg_lev_intopt:
            if value is None: return self.iocp.msg_lev
            else: self.iocp.msg_lev = value

        elif name == br_tech:
            if value is None: return self.iocp.br_tech
            else: self.iocp.br_tech = value

        elif name == bt_tech:
            if value is None: return self.iocp.bt_tech
            else: self.iocp.bt_tech = value

        elif name == pp_tech:
            if value is None: return self.iocp.pp_tech
            else: self.iocp.pp_tech = value

        elif name == fp_heur:
            if value is None: return self.iocp.fp_heur
            else:
                if value: value = GLP_ON
                else: value = GLP_OFF
                self.iocp.fp_heur = value

        elif name == gmi_cuts:
            if value is None: return self.iocp.gmi_cuts
            else:
                if value: value = GLP_ON
                else: value = GLP_OFF
                self.iocp.gmi_cuts = value

        elif name == mir_cuts:
            if value is None: return self.iocp.mir_cuts
            else:
                if value: value = GLP_ON
                else: value = GLP_OFF
                self.iocp.mir_cuts = value

        elif name == cov_cuts:
            if value is None: return self.iocp.cov_cuts
            else:
                if value: value = GLP_ON
                else: value = GLP_OFF
                self.iocp.cov_cuts = value

        elif name == clq_cuts:
            if value is None: return self.iocp.clq_cuts
            else:
                if value: value = GLP_ON
                else: value = GLP_OFF
                self.iocp.clq_cuts = value

        elif name == tol_int:
            if value is None: return self.iocp.tol_int
            else: self.iocp.tol_int = value

        elif name == tol_obj:
            if value is None: return self.iocp.tol_obj
            else: self.iocp.tol_obj = value

        elif name == mip_gap:
            if value is None: return self.iocp.mip_gap
            else: self.iocp.mip_gap = value

        elif name == tm_lim_intopt:
            if value is None: return self.iocp.tm_lim
            else: self.iocp.tm_lim = value

        elif name == out_frq_intopt:
            if value is None: return self.iocp.out_frq
            else: self.iocp.out_frq = value

        elif name == out_dly_intopt:
            if value is None: return self.iocp.out_dly
            else: self.iocp.out_dly = value

        elif name == presolve_intopt:
            if value is None: return self.iocp.presolve
            else:
                if value: value = GLP_ON
                else: value = GLP_OFF
                self.iocp.presolve = value

        elif name == binarize:
            if value is None: return self.iocp.binarize
            else:
                if value: value = GLP_ON
                else: value = GLP_OFF
                self.iocp.binarize = value

        elif name == msg_lev_simplex:
            if value is None: return self.smcp.msg_lev
            else: self.smcp.msg_lev = value

        elif name == meth:
            if value is None: return self.smcp.meth
            else: self.smcp.meth = value

        elif name == pricing:
            if value is None: return self.smcp.pricing
            else: self.smcp.pricing = value

        elif name == r_test:
            if value is None: return self.smcp.r_test
            else: self.smcp.r_test = value

        elif name == tol_bnd:
            if value is None: return self.smcp.tol_bnd
            else: self.smcp.tol_bnd = value

        elif name == tol_dj:
            if value is None: return self.smcp.tol_dj
            else: self.smcp.tol_dj = value

        elif name == tol_piv:
            if value is None: return self.smcp.tol_piv
            else: self.smcp.tol_piv = value

        elif name == obj_ll:
            if value is None: return self.smcp.obj_ll
            else: self.smcp.obj_ll = value

        elif name == obj_ul:
            if value is None: return self.smcp.obj_ul
            else: self.smcp.obj_ul = value

        elif name == it_lim:
            if value is None: return self.smcp.it_lim
            else: self.smcp.it_lim = value

        elif name == tm_lim_simplex:
            if value is None: return self.smcp.tm_lim
            else: self.smcp.tm_lim = value

        elif name == out_frq_simplex:
            if value is None: return self.smcp.out_frq
            else: self.smcp.out_frq = value

        elif name == out_dly_simplex:
            if value is None: return self.smcp.out_dly
            else: self.smcp.out_dly = value

        elif name == presolve_simplex:
            if value is None: return self.smcp.presolve
            else:
                if value: value = GLP_ON
                else: value = GLP_OFF
                self.smcp.presolve = value
        else:
            raise ValueError("This parameter is not available.")

    cpdef bint is_variable_basic(self, int index) noexcept:
        """
        Test whether the given variable is basic.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver='GLPK')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0
            sage: b.is_variable_basic(0)
            True
            sage: b.is_variable_basic(1)
            False
        """
        return self.get_col_stat(index) == GLP_BS

    cpdef bint is_variable_nonbasic_at_lower_bound(self, int index) noexcept:
        """
        Test whether the given variable is nonbasic at lower bound.
        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver='GLPK')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0
            sage: b.is_variable_nonbasic_at_lower_bound(0)
            False
            sage: b.is_variable_nonbasic_at_lower_bound(1)
            True
        """
        return self.get_col_stat(index) == GLP_NL

    cpdef bint is_slack_variable_basic(self, int index) noexcept:
        """
        Test whether the slack variable of the given row is basic.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver='GLPK')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0
            sage: b.is_slack_variable_basic(0)
            True
            sage: b.is_slack_variable_basic(1)
            False
        """
        return self.get_row_stat(index) == GLP_BS

    cpdef bint is_slack_variable_nonbasic_at_lower_bound(self, int index) noexcept:
        """
        Test whether the slack variable of the given row is nonbasic at lower bound.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver='GLPK')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(5.5 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0
            sage: b.is_slack_variable_nonbasic_at_lower_bound(0)
            False
            sage: b.is_slack_variable_nonbasic_at_lower_bound(1)
            True
        """
        return self.get_row_stat(index) == GLP_NU

    cpdef int print_ranges(self, filename=None) except -1:
        r"""
        Print results of a sensitivity analysis

        If no filename is given as an input the results of the
        sensitivity analysis are displayed on the screen. If a
        filename is given they are written to a file.

        INPUT:

        - ``filename`` -- (optional) name of the file

        OUTPUT: zero if the operations was successful otherwise nonzero

        .. NOTE::

            This method is only effective if an optimal solution has been found
            for the lp using the simplex algorithm. In all other cases an error
            message is printed.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint(list(zip([0, 1], [1, 2])), None, 3)
            sage: p.set_objective([2, 5])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: p.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: p.print_ranges()
            glp_print_ranges: optimal basic solution required
            1
            sage: p.solve()
            0
            sage: from tempfile import NamedTemporaryFile
            sage: with NamedTemporaryFile(mode='r+t', suffix='.tmp') as f:
            ....:     p.print_ranges(f.name)
            ....:     for ll in f.readlines():
            ....:         if ll: print(ll)
            ...
            GLPK ... - SENSITIVITY ANALYSIS REPORT                                                                         Page   1
            Problem:
            Objective:  7.5 (MAXimum)
               No. Row name     St      Activity         Slack   Lower bound       Activity      Obj coef  Obj value at Limiting
                                                      Marginal   Upper bound          range         range   break point variable
            ------ ------------ -- ------------- ------------- -------------  ------------- ------------- ------------- ------------
                 1              NU       3.00000        .               -Inf         .           -2.50000        .
                                                       2.50000       3.00000           +Inf          +Inf          +Inf
            GLPK ... - SENSITIVITY ANALYSIS REPORT                                                                         Page   2
            Problem:
            Objective:  7.5 (MAXimum)
               No. Column name  St      Activity      Obj coef   Lower bound       Activity      Obj coef  Obj value at Limiting
                                                      Marginal   Upper bound          range         range   break point variable
            ------ ------------ -- ------------- ------------- -------------  ------------- ------------- ------------- ------------
                 1              NL        .            2.00000        .                -Inf          -Inf          +Inf
                                                       -.50000          +Inf        3.00000       2.50000       6.00000
                 2              BS       1.50000       5.00000        .                -Inf       4.00000       6.00000
                                                        .               +Inf        1.50000          +Inf          +Inf
            End of report
        """
        if filename is None:
            import tempfile
            with tempfile.NamedTemporaryFile() as f:
                res = glp_print_ranges(self.lp, 0, 0, 0,
                                       str_to_bytes(f.name,
                                                    FS_ENCODING,
                                                    'surrogateescape'))
                if res == 0:
                    # Success; now write it to stdout.
                    with open(f.name) as fh:
                        for line in fh:
                            print(line, end=" ")
                    print("\n")
        else:
            res = glp_print_ranges(self.lp, 0, 0, 0,
                                   str_to_bytes(filename,
                                                FS_ENCODING,
                                                'surrogateescape'))
        return res

    cpdef double get_row_dual(self, int variable) noexcept:
        r"""
        Return the dual value of a constraint.

        The dual value of the i-th row is also the value of the i-th variable
        of the dual problem.

        The dual value of a constraint is the shadow price of the constraint.
        The shadow price is the amount by which the objective value will change
        if the constraints bounds change by one unit under the precondition
        that the basis remains the same.

        INPUT:

        - ``variable`` -- the number of the constraint

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.
           If the simplex algorithm has not been used for solving 0.0 will
           be returned.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "GLPK")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: lp.solve()
            0
            sage: lp.get_row_dual(0)   # tolerance 0.00001
            0.0
            sage: lp.get_row_dual(1)   # tolerance 0.00001
            10.0
        """

        if self.simplex_or_intopt == simplex_only:
            return glp_get_row_dual(self.lp, variable+1)
        else:
            return 0.0

    cpdef double get_col_dual(self, int variable) except? -1:
        """
        Return the dual value (reduced cost) of a variable

        The dual value is the reduced cost of a variable.
        The reduced cost is the amount by which the objective coefficient
        of a non basic variable has to change to become a basic variable.

        INPUT:

        - ``variable`` -- the number of the variable

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.
           If the simplex algorithm has not been used for solving just a
           0.0 will be returned.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(3)
            2
            sage: p.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: p.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: p.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: p.set_objective([60, 30, 20])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: p.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: p.solve()
            0
            sage: p.get_col_dual(1)
            -5.0

        TESTS:

        We sanity check the input that will be passed to GLPK::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: p.get_col_dual(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid column index 2
        """
        if variable < 0 or variable > (self.ncols() - 1):
            raise ValueError("invalid column index %d" % variable)

        if self.simplex_or_intopt == simplex_only:
            return glp_get_col_dual(self.lp, variable+1)
        else:
            return 0.0

    cpdef int get_row_stat(self, int i) except? -1:
        """
        Retrieve the status of a constraint.

        INPUT:

        - ``i`` -- the index of the constraint

        OUTPUT:

        Current status assigned to the auxiliary variable associated with i-th
        row:

            * GLP_BS = 1     basic variable
            * GLP_NL = 2     non-basic variable on lower bound
            * GLP_NU = 3     non-basic variable on upper bound
            * GLP_NF = 4     non-basic free (unbounded) variable
            * GLP_NS = 5     non-basic fixed variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "GLPK")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: lp.solve()
            0
            sage: lp.get_row_stat(0)
            1
            sage: lp.get_row_stat(1)
            3
            sage: lp.get_row_stat(-1)
            Traceback (most recent call last):
            ...
            ValueError: The constraint's index i must satisfy 0 <= i < number_of_constraints
        """
        if i < 0 or i >= glp_get_num_rows(self.lp):
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")
        return glp_get_row_stat(self.lp, i+1)

    cpdef int get_col_stat(self, int j) except? -1:
        """
        Retrieve the status of a variable.

        INPUT:

        - ``j`` -- the index of the variable

        OUTPUT:

        Current status assigned to the structural variable associated
        with j-th column:

            * GLP_BS = 1     basic variable
            * GLP_NL = 2     non-basic variable on lower bound
            * GLP_NU = 3     non-basic variable on upper bound
            * GLP_NF = 4     non-basic free (unbounded) variable
            * GLP_NS = 5     non-basic fixed variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "GLPK")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: lp.solve()
            0
            sage: lp.get_col_stat(0)
            1
            sage: lp.get_col_stat(1)
            2
            sage: lp.get_col_stat(100)
            Traceback (most recent call last):
            ...
            ValueError: The variable's index j must satisfy 0 <= j < number_of_variables
        """
        if j < 0 or j >= glp_get_num_cols(self.lp):
            raise ValueError("The variable's index j must satisfy 0 <= j < number_of_variables")

        return glp_get_col_stat(self.lp, j+1)

    cpdef set_row_stat(self, int i, int stat):
        r"""
        Set the status of a constraint.

        INPUT:

        - ``i`` -- the index of the constraint

        - ``stat`` -- the status to set to

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "GLPK")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: lp.solve()
            0
            sage: lp.get_row_stat(0)
            1
            sage: lp.set_row_stat(0, 3)
            sage: lp.get_row_stat(0)
            3
        """
        if i < 0 or i >= glp_get_num_rows(self.lp):
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")

        glp_set_row_stat(self.lp, i+1, stat)

    cpdef set_col_stat(self, int j, int stat):
        r"""
        Set the status of a variable.

        INPUT:

        - ``j`` -- the index of the constraint

        - ``stat`` -- the status to set to

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "GLPK")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: lp.solve()
            0
            sage: lp.get_col_stat(0)
            1
            sage: lp.set_col_stat(0, 2)
            sage: lp.get_col_stat(0)
            2
        """
        if j < 0 or j >= glp_get_num_cols(self.lp):
            raise ValueError("The variable's index j must satisfy 0 <= j < number_of_variables")

        glp_set_col_stat(self.lp, j+1, stat)

    cpdef int warm_up(self) noexcept:
        r"""
        Warm up the basis using current statuses assigned to rows and cols.

        OUTPUT: the warming up status

            * 0             The operation has been successfully performed.
            * GLP_EBADB     The basis matrix is invalid.
            * GLP_ESING     The basis matrix is singular within the working precision.
            * GLP_ECOND     The basis matrix is ill-conditioned.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "GLPK")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: lp.solve()
            0
            sage: lp.get_objective_value()
            280.0
            sage: lp.set_row_stat(0,3)
            sage: lp.set_col_stat(1,1)
            sage: lp.warm_up()
            0
        """
        return glp_warm_up(self.lp)

    cpdef eval_tab_row(self, int k):
        r"""
        Compute a row of the current simplex tableau.

        A row corresponds to some basic variable specified by the parameter
        ``k`` as follows:

        - if `0 \leq k \leq m-1`, the basic variable is `k`-th auxiliary
          variable,

        - if `m \leq k \leq m+n-1`, the basic variable is `(k-m)`-th structural
          variable,

        where `m` is the number of rows and `n` is the number of columns in the
        specified problem object.

        .. NOTE::

            The basis factorization must exist and the variable with
            index ``k`` must be basic. Otherwise, a :exc:`ValueError` is
            be raised.

        INPUT:

        - ``k`` -- integer; the id of the basic variable

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient in the computed row
        of the current simplex tableau.

        .. NOTE::

            Elements in ``indices`` have the same sense as index ``k``.
            All these variables are non-basic by definition.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "GLPK")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: lp.eval_tab_row(0)
            Traceback (most recent call last):
            ...
            ValueError: basis factorization does not exist
            sage: lp.solve()
            0
            sage: lp.eval_tab_row(0)
            ([1, 2, 4], [-2.0, 8.0, -2.0])
            sage: lp.eval_tab_row(3)
            ([1, 2, 4], [-0.5, 1.5, -1.25])
            sage: lp.eval_tab_row(5)
            ([1, 2, 4], [2.0, -4.0, 2.0])
            sage: lp.eval_tab_row(1)
            Traceback (most recent call last):
            ...
            ValueError: slack variable 1 is not basic
            sage: lp.eval_tab_row(-1)
            Traceback (most recent call last):
            ...
            ValueError: ...
        """
        cdef int m = self.nrows()
        cdef int n = self.ncols()
        cdef int i,j

        if k < 0 or k >= n + m:
            raise ValueError("k = %s; Variable number out of range" % k)

        if glp_bf_exists(self.lp) == 0:
            raise ValueError("basis factorization does not exist")

        if k < m:
            if not self.is_slack_variable_basic(k):
                raise ValueError("slack variable %d is not basic" % k)
        else:
            if not self.is_variable_basic(k-m):
                raise ValueError("variable %d is not basic" % (k-m) )

        cdef MemoryAllocator mem = MemoryAllocator()
        cdef int    * c_indices = <int*>mem.allocarray(n+1, sizeof(int))
        cdef double * c_values  = <double*>mem.allocarray(n+1, sizeof(double))

        i = glp_eval_tab_row(self.lp, k + 1, c_indices, c_values)

        indices = [c_indices[j+1] - 1 for j in range(i)]
        values  = [c_values[j+1]      for j in range(i)]
        return (indices, values)

    cpdef eval_tab_col(self, int k):
        r"""
        Compute a column of the current simplex tableau.

        A (column) corresponds to some non-basic variable specified by the
        parameter ``k`` as follows:

        - if `0 \leq k \leq m-1`, the non-basic variable is `k`-th auxiliary
          variable,

        - if `m \leq k \leq m+n-1`, the non-basic variable is `(k-m)`-th
          structural variable,

        where `m` is the number of rows and `n` is the number of columns
        in the specified problem object.

        .. NOTE::

            The basis factorization must exist and the variable with
            index ``k`` must not be basic. Otherwise, a :exc:`ValueError` is
            be raised.

        INPUT:

        - ``k`` -- integer; the id of the non-basic variable

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient in the computed column
        of the current simplex tableau.

        .. NOTE::

            Elements in ``indices`` have the same sense as index `k`.
            All these variables are basic by definition.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "GLPK")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: lp.eval_tab_col(1)
            Traceback (most recent call last):
            ...
            ValueError: basis factorization does not exist
            sage: lp.solve()
            0
            sage: lp.eval_tab_col(1)
            ([0, 5, 3], [-2.0, 2.0, -0.5])
            sage: lp.eval_tab_col(2)
            ([0, 5, 3], [8.0, -4.0, 1.5])
            sage: lp.eval_tab_col(4)
            ([0, 5, 3], [-2.0, 2.0, -1.25])
            sage: lp.eval_tab_col(0)
            Traceback (most recent call last):
            ...
            ValueError: slack variable 0 is basic
            sage: lp.eval_tab_col(-1)
            Traceback (most recent call last):
            ...
            ValueError: ...
        """
        cdef int m = self.nrows()
        cdef int n = self.ncols()
        cdef int i,j

        if k < 0 or k >= m + n:
            raise ValueError("k = %s; Variable number out of range" % k)

        if glp_bf_exists(self.lp) == 0:
            raise ValueError("basis factorization does not exist")

        if k < m:
            if self.is_slack_variable_basic(k):
                raise ValueError("slack variable %d is basic" % k)
        else:
            if self.is_variable_basic(k-m):
                raise ValueError("variable %d is basic" % (k-m) )

        cdef MemoryAllocator mem = MemoryAllocator()
        cdef int    * c_indices = <int*>mem.allocarray(m+1, sizeof(int))
        cdef double * c_values  = <double*>mem.allocarray(m+1, sizeof(double))

        i = glp_eval_tab_col(self.lp, k + 1, c_indices, c_values)

        indices = [c_indices[j+1] - 1 for j in range(i)]
        values  = [c_values[j+1]      for j in range(i)]
        return (indices, values)

    def __dealloc__(self):
        """
        Destructor
        """
        glp_delete_prob(self.lp)
        sig_free(self.iocp)
        sig_free(self.smcp)

cdef void glp_callback(glp_tree* tree, void* info) noexcept:
    r"""
    A callback routine called by glp_intopt

    This function fills the ``search_tree_data`` structure of a ``GLPKBackend``
    object. It is fed to ``glp_intopt`` (which calls it while the search tree is
    being built) through ``iocp.cb_func``.

    INPUT:

    - ``tree`` -- a pointer toward ``glp_tree``, which is GLPK's search tree

    - ``info`` -- a ``void *`` to let the function know *where* it should store
      the data we need. The value of ``info`` is equal to the one stored in
      iocp.cb_info.
    """
    cdef search_tree_data_t * data = <search_tree_data_t *> info
    data.mip_gap = glp_ios_mip_gap(tree)

    cdef int node_id = glp_ios_best_node(tree)
    data.best_bound = glp_ios_node_bound(tree, node_id)

# parameter names

cdef enum solver_parameter_names:
    timelimit_seconds, timelimit_simplex, timelimit_intopt, simplex_or_intopt,
    msg_lev_simplex,
    msg_lev_intopt, br_tech, bt_tech, pp_tech, fp_heur, gmi_cuts,
    mir_cuts, cov_cuts, clq_cuts, tol_int, tol_obj, mip_gap,
    tm_lim_intopt, out_frq_intopt, out_dly_intopt, presolve_intopt,
    binarize, meth, pricing, r_test, tol_bnd, tol_dj, tol_piv, obj_ll,
    obj_ul, it_lim, tm_lim_simplex, out_frq_simplex, out_dly_simplex,
    presolve_simplex

glp_tm_lim_simplex = timelimit_simplex
glp_tm_lim_intopt = timelimit_intopt
glp_simplex_or_intopt = simplex_or_intopt
glp_msg_lev_intopt = glp_verbosity_intopt = msg_lev_intopt
glp_msg_lev_simplex = glp_verbosity_simplex = msg_lev_simplex
glp_br_tech = glp_branching = br_tech
glp_bt_tech = glp_backtracking = bt_tech
glp_pp_tech = glp_preprocessing = pp_tech
glp_fp_heur = glp_feasibility_pump = fp_heur
glp_gmi_cuts = glp_gomory_cuts = gmi_cuts
glp_mir_cuts = glp_mixed_int_rounding_cuts = mir_cuts
glp_cov_cuts = glp_mixed_cover_cuts = cov_cuts
glp_clq_cuts = glp_clique_cuts = clq_cuts
glp_tol_int = glp_absolute_tolerance = tol_int
glp_tol_obj = glp_relative_tolerance = tol_obj
glp_mip_gap = glp_mip_gap_tolerance = mip_gap
glp_out_frq_intopt = glp_output_frequency_intopt = out_frq_intopt
glp_out_dly_intopt = glp_output_delay_intopt = out_dly_intopt
glp_presolve_intopt = presolve_intopt
glp_binarize = binarize
glp_meth = glp_primal_v_dual = meth
glp_pricing = pricing
glp_r_test = glp_ratio_test = r_test
glp_tol_bnd = glp_tolerance_primal = tol_bnd
glp_tol_dj = glp_tolerance_dual = tol_dj
glp_tol_piv = glp_tolerance_pivot = tol_piv
glp_obj_ll = glp_obj_lower_limit = obj_ll
glp_obj_ul = glp_obj_upper_limit = obj_ul
glp_it_lim = glp_iteration_limit = it_lim
glp_out_frq_simplex = glp_output_frequency_intopt = out_frq_simplex
glp_out_dly_simplex = glp_output_delay_simplex = out_dly_simplex
glp_presolve_simplex = presolve_simplex

solver_parameter_names_dict = {
  'timelimit': timelimit_seconds,
  'timelimit_intopt': timelimit_intopt,
  'tm_lim_intopt': timelimit_intopt,
  'timelimit_simplex': timelimit_simplex,
  'tm_lim_simplex': timelimit_simplex,
  'simplex_or_intopt': simplex_or_intopt,
  'msg_lev_simplex': msg_lev_simplex, 'verbosity_simplex': msg_lev_simplex,
  'msg_lev_intopt': msg_lev_intopt, 'verbosity_intopt': msg_lev_intopt,
  'br_tech': br_tech, 'branching': br_tech,
  'bt_tech': bt_tech, 'backtracking': bt_tech,
  'pp_tech': pp_tech, 'preprocessing': pp_tech,
  'fp_heur': fp_heur, 'feasibility_pump': fp_heur,
  'gmi_cuts': gmi_cuts, 'gomory_cuts': gmi_cuts,
  'mir_cuts': mir_cuts, 'mixed_int_rounding_cuts': mir_cuts,
  'cov_cuts': cov_cuts, 'mixed_cover_cuts': cov_cuts,
  'clq_cuts': clq_cuts, 'clique_cuts': clq_cuts,
  'tol_int': tol_int, 'absolute_tolerance': tol_int,
  'tol_obj': tol_obj, 'relative_tolerance': tol_obj,
  'mip_gap': mip_gap, 'mip_gap_tolerance': mip_gap,
  'out_frq_intopt': out_frq_intopt, 'output_frequency_intopt': out_frq_intopt,
  'out_dly_intopt': out_dly_intopt, 'output_delay_intopt': out_dly_intopt,
  'presolve_intopt': presolve_intopt, 'binarize': binarize,
  'meth': meth, 'primal_v_dual': meth,
  'pricing': pricing,
  'r_test': r_test, 'ratio_test': r_test,
  'tol_bnd': tol_bnd, 'tolerance_primal': tol_bnd,
  'tol_dj': tol_dj, 'tolerance_dual': tol_dj,
  'tol_piv': tol_piv, 'tolerance_pivot': tol_piv,
  'obj_ll': obj_ll, 'obj_lower_limit': obj_ll,
  'obj_ul': obj_ul, 'obj_upper_limit': obj_ul,
  'it_lim': it_lim, 'iteration_limit': it_lim,
  'out_frq_simplex': out_frq_simplex, 'output_frequency_simplex': out_frq_simplex,
  'out_dly_simplex': out_dly_simplex, 'output_delay_simplex': out_dly_simplex,
  'presolve_simplex': presolve_simplex
}

# parameter values

glp_msg_off = GLP_MSG_OFF
glp_msg_on = GLP_MSG_ON
glp_msg_err = GLP_MSG_ERR
glp_msg_all = GLP_MSG_ALL
glp_msg_dbg = GLP_MSG_DBG

glp_primal = GLP_PRIMAL
glp_dual = GLP_DUAL
glp_dualp = GLP_DUALP

glp_pt_std = GLP_PT_STD
glp_pt_pse = GLP_PT_PSE

glp_rt_std = GLP_RT_STD
glp_rt_har = GLP_RT_HAR

dbl_max = DBL_MAX
int_max = INT_MAX

glp_on = GLP_ON
glp_off = GLP_OFF

glp_br_ffv = GLP_BR_FFV
glp_br_lfv = GLP_BR_LFV
glp_br_mfv = GLP_BR_MFV
glp_br_dth = GLP_BR_DTH
glp_br_pch = GLP_BR_PCH

glp_bt_dfs = GLP_BT_DFS
glp_bt_bfs = GLP_BT_BFS
glp_bt_blb = GLP_BT_BLB
glp_bt_bph = GLP_BT_BPH

glp_pp_none = GLP_PP_NONE
glp_pp_root = GLP_PP_ROOT
glp_pp_all = GLP_PP_ALL

glp_max = GLP_MAX
glp_min = GLP_MIN
glp_up = GLP_UP
glp_fr = GLP_FR
glp_db = GLP_DB
glp_fx = GLP_FX
glp_lo = GLP_LO
glp_cv = GLP_CV
glp_iv = GLP_IV
glp_bv = GLP_BV
glp_mps_deck = GLP_MPS_DECK
glp_mps_file = GLP_MPS_FILE

glp_undef = GLP_UNDEF
glp_opt = GLP_OPT
glp_feas = GLP_FEAS
glp_nofeas = GLP_NOFEAS

glp_bs = GLP_BS
glp_nl = GLP_NL
glp_nu = GLP_NU
glp_nf = GLP_NF

cdef enum more_parameter_values:
    simplex_only, simplex_then_intopt, intopt_only, exact_simplex_only

glp_simplex_only = simplex_only
glp_simplex_then_intopt = simplex_then_intopt
glp_intopt_only = intopt_only
glp_exact_simplex_only = exact_simplex_only

# dictionaries for those who prefer to use strings

solver_parameter_values = {

  'simplex_only': simplex_only,
  'simplex_then_intopt': simplex_then_intopt,
  'intopt_only': intopt_only,
  'exact_simplex_only': exact_simplex_only,

  'GLP_MSG_OFF' : GLP_MSG_OFF,
  'GLP_MSG_ON' : GLP_MSG_ON,
  'GLP_MSG_ERR' : GLP_MSG_ERR,
  'GLP_MSG_ALL' : GLP_MSG_ALL,
  'GLP_MSG_DBG' : GLP_MSG_DBG,

  'GLP_PRIMAL' : GLP_PRIMAL,
  'GLP_DUAL' : GLP_DUAL,
  'GLP_DUALP' : GLP_DUALP,

  'GLP_PT_STD' : GLP_PT_STD,
  'GLP_PT_PSE' : GLP_PT_PSE,

  'GLP_RT_STD' : GLP_RT_STD,
  'GLP_RT_HAR' : GLP_RT_HAR,

  'DBL_MAX' : DBL_MAX,
  'INT_MAX' : INT_MAX,

  'GLP_ON' : GLP_ON,
  'GLP_OFF' : GLP_OFF,

  'GLP_BR_FFV' : GLP_BR_FFV,
  'GLP_BR_LFV' : GLP_BR_LFV,
  'GLP_BR_MFV' : GLP_BR_MFV,
  'GLP_BR_DTH' : GLP_BR_DTH,
  'GLP_BR_PCH' : GLP_BR_PCH,

  'GLP_BT_DFS' : GLP_BT_DFS,
  'GLP_BT_BFS' : GLP_BT_BFS,
  'GLP_BT_BLB' : GLP_BT_BLB,
  'GLP_BT_BPH' : GLP_BT_BPH,

  'GLP_PP_NONE' : GLP_PP_NONE,
  'GLP_PP_ROOT' : GLP_PP_ROOT,
  'GLP_PP_ALL' : GLP_PP_ALL,

  'GLP_MAX' : GLP_MAX,
  'GLP_MIN' : GLP_MIN,
  'GLP_UP' : GLP_UP,
  'GLP_FR' : GLP_FR,
  'GLP_DB' : GLP_DB,
  'GLP_FX' : GLP_FX,
  'GLP_LO' : GLP_LO,
  'GLP_CV' : GLP_CV,
  'GLP_IV' : GLP_IV,
  'GLP_BV' : GLP_BV,
  'GLP_MPS_DECK' : GLP_MPS_DECK,
  'GLP_MPS_FILE' : GLP_MPS_FILE,

  'GLP_UNDEF' : GLP_UNDEF,
  'GLP_OPT' : GLP_OPT,
  'GLP_FEAS' : GLP_FEAS,
  'GLP_NOFEAS' : GLP_NOFEAS

}

cdef dict solve_status_msg = {
    GLP_EBADB: "The initial basis specified in the problem object is invalid",
    GLP_ESING: "The basis matrix corresponding to the initial basis is singular within the working precision",
    GLP_ECOND: "The basis matrix corresponding to the initial basis is ill-conditioned, i.e. its condition number is too large",
    GLP_EBOUND: "Some variables (auxiliary or structural) have incorrect bounds",
    GLP_EFAIL: "Solver failure",
    GLP_EOBJLL: "The objective lower limit has been reached",
    GLP_EOBJUL: "The objective upper limit has been reached",
    GLP_EITLIM: "The iteration limit has been exceeded",
    GLP_ETMLIM: "The time limit has been exceeded",
    GLP_ENOPFS: "The LP (relaxation) problem has no primal feasible solution",
    GLP_ENODFS: "The LP (relaxation) problem has no dual feasible solution",
    GLP_EROOT:  "Optimal basis for initial LP relaxation is not provided",
    GLP_ESTOP:  "The search was prematurely terminated by application",
    GLP_EMIPGAP: "The relative mip gap tolerance has been reached",
}

cdef dict solution_status_msg = {
    GLP_UNDEF: "Solution is undefined",
    GLP_FEAS: "Feasible solution found, while optimality has not been proven",
    GLP_INFEAS: "Solution is infeasible",
    GLP_NOFEAS: "Problem has no feasible solution",
    GLP_OPT: "Solution is optimal",
    GLP_UNBND: "Problem has unbounded solution",
}
