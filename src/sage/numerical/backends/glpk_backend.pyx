"""
GLPK Backend

AUTHORS:

- Nathann Cohen (2010-10): initial implementation
- John Perry (2012-01): glp_simplex preprocessing
- John Perry and Raniere Gaia Silva (2012-03): solver parameters
- Christian Kuper (2012-10): Additions for sensitivity analysis
"""

##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.numerical.mip import MIPSolverException

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"

cdef class GLPKBackend(GenericBackend):

    def __cinit__(self, maximization = True):
        """
        Constructor

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver="GLPK")
        """
        self.lp = glp_create_prob()
        self.simplex_or_intopt = glp_intopt_only
        self.smcp = <c_glp_smcp* > sage_malloc(sizeof(c_glp_smcp))
        glp_init_smcp(self.smcp)
        self.iocp = <c_glp_iocp* > sage_malloc(sizeof(c_glp_iocp))
        glp_init_iocp(self.iocp)
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

        - ``lower_bound`` - the lower bound of the variable (default: 0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is binary (default: ``True``).

        - ``integer`` - ``True`` if the variable is binary (default: ``False``).

        - ``obj`` - (optional) coefficient of this variable in the objective function (default: 0.0)

        - ``name`` - an optional name for the newly added variable (default: ``None``).

        OUTPUT: The index of the newly created variable

        EXAMPLE::

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
            sage: p.add_variable(name='x',obj=1.0)
            3
            sage: p.col_name(3)
            'x'
            sage: p.objective_coefficient(3)
            1.0
        """
        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if  vtype == 0:
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
            glp_set_col_name(self.lp, n_var, name)

        if obj:
            self.objective_coefficient(n_var - 1, obj)

        return n_var - 1

    cpdef int add_variables(self, int number, lower_bound=0.0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0.0, names=None) except -1:
        """
        Add ``number`` new variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive, real and theor coefficient in
        the objective function is 0.0.

        INPUT:

        - ``n`` - the number of new variables (must be > 0)

        - ``lower_bound`` - the lower bound of the variable (default: 0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is binary (default: ``True``).

        - ``integer`` - ``True`` if the variable is binary (default: ``False``).

        - ``obj`` - (optional) coefficient of all variables in the objective function (default: 0.0)

        - ``names`` - optional list of names (default: ``None``)

        OUTPUT: The index of the variable created last.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variables(5)
            4
            sage: p.ncols()
            5
            sage: p.add_variables(2, lower_bound=-2.0, integer=True, names=['a','b'])
            6
        """
        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if  vtype == 0:
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
                glp_set_col_name(self.lp, n_var, names[number - i - 1])

        return n_var - 1

    cpdef set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``vtype`` (integer) :

            *  1  Integer
            *  0  Binary
            * -1 Real

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)
            True
        """

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

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

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
        Set or get the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient or ``None`` for
          reading (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0.0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2.0
        """
        if coeff is None:
            return glp_get_obj_coef(self.lp, variable + 1)
        else:
            glp_set_obj_coef(self.lp, variable + 1, coeff)

    cpdef problem_name(self, char * name = NULL):
        """
        Return or define the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.problem_name("There once was a french fry")
            sage: print p.problem_name()
            There once was a french fry
        """
        cdef char * n

        if name == NULL:
            n =  <char *> glp_get_prob_name(self.lp)
            if n == NULL:
                return ""
            else:
                return n

        else:
            if len(name) > 255:
                raise ValueError("Problem name for GLPK must not be longer than 255 characters.")
            glp_set_prob_name(self.lp, name)

    cpdef set_objective(self, list coeff, d = 0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` - a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        - ``d`` (double) -- the constant term in the linear function (set to `0` by default)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: map(lambda x :p.objective_coefficient(x), range(5))
            [1.0, 1.0, 2.0, 1.0, 3.0]
        """
        cdef int i

        for i,v in enumerate(coeff):
            glp_set_obj_coef(self.lp, i+1, v)

        glp_set_obj_coef(self.lp, 0, d)

        self.obj_constant_term = d

    cpdef set_verbosity(self, int level):
        """
        Set the verbosity level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.set_verbosity(2)

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
        Remove a constraint from self.

        INPUT:

        - ``i`` -- index of the constraint to remove

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver='GLPK')
            sage: x,y = p[0], p[1]
            doctest:...: DeprecationWarning: The default value of 'nonnegative' will change, to False instead of True. You should add the explicit 'nonnegative=True'.
            See http://trac.sagemath.org/15521 for details.
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
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

        - ``constraints`` -- an iterable containing the indices of the rows to remove.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver='GLPK')
            sage: x,y = p[0], p[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: p.set_integer(x); p.set_integer(y)
            sage: p.solve()
            9.0
            sage: p.remove_constraints([1])
            sage: p.solve()
            10.0
            sage: p.get_values([x,y])
            [3.0, 0.0]

        TESTS:

        Removing fancy constraints does not make Sage crash::

            sage: MixedIntegerLinearProgram(solver= "GLPK").remove_constraints([0, -2])
            Traceback (most recent call last):
            ...
            ValueError: The constraint's index i must satisfy 0 <= i < number_of_constraints
        """
        cdef int i, c
        cdef int m = len(constraints)
        cdef int * rows = <int *>sage_malloc((m + 1) * sizeof(int *))
        cdef int nrows = glp_get_num_rows(self.lp)

        for i in xrange(m):

            c = constraints[i]
            if c < 0 or c >= nrows:
                sage_free(rows)
                raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")

            rows[i+1] = c + 1

        glp_del_rows(self.lp, m, rows)
        sage_free(rows)
        glp_std_basis(self.lp)

    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` an iterable with ``(c,v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (real
          value).

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``name`` - an optional name for this row (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint( zip(range(5), range(5)), 2.0, 2.0)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1.0, 1.0, name='foo')
            sage: p.row_name(1)
            'foo'
        """
        if lower_bound is None and upper_bound is None:
            raise ValueError("At least one of 'upper_bound' or 'lower_bound' must be set.")

        glp_add_rows(self.lp, 1)
        cdef int n = glp_get_num_rows(self.lp)

        cdef int * row_i
        cdef double * row_values

        row_i = <int *> sage_malloc((len(coefficients)+1) * sizeof(int))
        row_values = <double *> sage_malloc((len(coefficients)+1) * sizeof(double))

        for i,(c,v) in enumerate(coefficients):
            row_i[i+1] = c+1
            row_values[i+1] = v

        glp_set_mat_row(self.lp, n, len(coefficients), row_i, row_values)

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
            glp_set_row_name(self.lp, n, name)

        sage_free(row_i)
        sage_free(row_values)

    cpdef add_linear_constraints(self, int number, lower_bound, upper_bound, names=None):
        """
        Add ``'number`` linear constraints.

        INPUT:

        - ``number`` (integer) -- the number of constraints to add.

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``names`` - an optional list of names (default: ``None``)

        EXAMPLE::

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
                glp_set_row_name(self.lp, n-i, names[number-i-1])

    cpdef row(self, int index):
        r"""
        Return a row

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_linear_constraint`` method.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)
        """
        cdef int n = glp_get_num_cols(self.lp)
        cdef int * c_indices = <int*> sage_malloc((n+1)*sizeof(int))
        cdef double * c_values = <double*> sage_malloc((n+1)*sizeof(double))
        cdef list indices = []
        cdef list values = []
        cdef int i,j

        i = glp_get_mat_row(self.lp, index + 1, c_indices, c_values)
        for 0 < j <= i:
            indices.append(c_indices[j]-1)
            values.append(c_values[j])

        sage_free(c_indices)
        sage_free(c_values)

        return (indices, values)

    cpdef row_bounds(self, int index):
        """
        Return the bounds of a specific constraint.

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the constraint is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)
        """
        cdef double ub
        cdef double lb

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

        - ``index`` (integer) -- the variable's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the variable is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0.0, 5.0)
        """

        cdef double ub
        cdef double lb

        ub = glp_get_col_ub(self.lp, index +1)
        lb = glp_get_col_lb(self.lp, index +1)

        return (
            (lb if lb != -DBL_MAX else None),
            (ub if ub != +DBL_MAX else None)
            )

    cpdef add_col(self, list indices, list coeffs):
        """
        Add a column.

        INPUT:

        - ``indices`` (list of integers) -- this list constains the
          indices of the constraints in which the variable's
          coefficient is nonzero

        - ``coeffs`` (list of real values) -- associates a coefficient
          to the variable in each of the constraints in which it
          appears. Namely, the ith entry of ``coeffs`` corresponds to
          the coefficient of the variable in the constraint
          represented by the ith entry in ``indices``.

        .. NOTE::

            ``indices`` and ``coeffs`` are expected to be of the same
            length.

        EXAMPLE::

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

        col_i = <int *> sage_malloc((len(indices)+1) * sizeof(int))
        col_values = <double *> sage_malloc((len(indices)+1) * sizeof(double))

        for i,v in enumerate(indices):
            col_i[i+1] = v+1
        for i,v in enumerate(coeffs):
            col_values[i+1] = v

        glp_set_mat_col(self.lp, n, len(indices), col_i, col_values)
        glp_set_col_bnds(self.lp, n, GLP_LO, 0,0)


    cpdef int solve(self) except -1:
        """
        Solve the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

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

            Sage uses GLPK's ``glp_intopt`` to find solutions.
            This routine sometimes FAILS CATASTROPHICALLY
            when given a system it cannot solve. (Ticket #12309.)
            Here, "catastrophic" can mean either "infinite loop" or
            segmentation fault. Upstream considers this behavior
            "essentially innate" to their design, and suggests
            preprocessing it with ``glpk_simplex`` first.
            Thus, if you suspect that your system is infeasible,
            set the ``preprocessing`` option first.

        EXAMPLE::

            sage: lp = MixedIntegerLinearProgram(solver = "GLPK")
            sage: lp.add_constraint( lp[1] +lp[2] -2.0 *lp[3] , max =-1.0)
            sage: lp.add_constraint( lp[0] -4.0/3 *lp[1] +1.0/3 *lp[2] , max =-1.0/3)
            sage: lp.add_constraint( lp[0] +0.5 *lp[1] -0.5 *lp[2] +0.25 *lp[3] , max =-0.25)
            sage: lp.solve()
            0.0
            sage: lp.add_constraint( lp[0] +4.0 *lp[1] -lp[2] +lp[3] , max =-1.0)
            sage: lp.solve()
            Traceback (most recent call last):
            ...
            RuntimeError: GLPK : Signal sent, try preprocessing option
            sage: lp.solver_parameter("simplex_or_intopt", "simplex_then_intopt")
            sage: lp.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: 'GLPK : Simplex cannot find a feasible solution'

        The user can ask sage to solve via ``simplex`` or ``intopt``.
        The default solver is ``intopt``, so we get integer solutions.

        EXAMPLE::

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

        If we switch to ``simplex``, we get continuous solutions.

        EXAMPLE::

            sage: lp.solver_parameter("simplex_or_intopt", "simplex_only") # use simplex only
            sage: lp.solve()
            2.0
            sage: lp.get_values([x, y])
            [1.5, 0.5]
        """

        cdef int status
        if self.simplex_or_intopt == glp_simplex_only or self.simplex_or_intopt == glp_simplex_then_intopt:
            status = glp_simplex(self.lp, self.smcp)
            status = glp_get_prim_stat(self.lp)
            if status == GLP_OPT or status == GLP_FEAS:
                pass
            elif status == GLP_UNDEF or status == GLP_NOFEAS:
                raise MIPSolverException("GLPK : Simplex cannot find a feasible solution")

        if self.simplex_or_intopt != glp_simplex_only:
          sig_str('GLPK : Signal sent, try preprocessing option')
          status = glp_intopt(self.lp, self.iocp)
          sig_off()
          # this is necessary to catch errors when certain options are enabled, e.g. tm_lim
          if status == GLP_ETMLIM: raise MIPSolverException("GLPK : The time limit was reached")
          elif status == GLP_EITLIM: raise MIPSolverException("GLPK : The iteration limit was reached")

          status = glp_mip_status(self.lp)
          if status == GLP_OPT:
              pass
          elif status == GLP_UNDEF:
              raise MIPSolverException("GLPK : Solution is undefined")
          elif status == GLP_FEAS:
              raise MIPSolverException("GLPK : Feasible solution found, while optimality has not been proven")
          elif status == GLP_INFEAS:
              raise MIPSolverException("GLPK : Solution is infeasible")
          elif status == GLP_UNBND:
              raise MIPSolverException("GLPK : Problem has unbounded solution")
          elif status == GLP_NOFEAS:
              raise MIPSolverException("GLPK : There is no feasible integer solution to this Linear Program")

        return 0

    cpdef get_objective_value(self):
        """
        Returns the value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

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
        if self.simplex_or_intopt != glp_simplex_only:
          return glp_mip_obj_val(self.lp)
        else:
          return glp_get_obj_val(self.lp)

    cpdef get_variable_value(self, int variable):
        """
        Returns the value of a variable given by the solver.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

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
        if self.simplex_or_intopt != glp_simplex_only:
          return glp_mip_col_val(self.lp, variable+1)
        else:
          return glp_get_col_prim(self.lp, variable+1)

    cpdef int ncols(self):
        """
        Return the number of columns/variables.

        EXAMPLE::

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

    cpdef int nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLE::

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
        Return the ``index`` th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable(name='I am a variable')
            0
            sage: p.col_name(0)
            'I am a variable'
        """
        cdef char * s

        glp_create_index(self.lp)
        s = <char*> glp_get_col_name(self.lp, index + 1)

        if s != NULL:
            return s
        else:
            return ""

    cpdef row_name(self, int index):
        """
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_linear_constraints(1, 2, None, names=['Empty constraint 1'])
            sage: p.row_name(0)
            'Empty constraint 1'
        """
        cdef char *  s

        glp_create_index(self.lp)
        s = <char*> glp_get_row_name(self.lp, index + 1)

        if s != NULL:
            return s
        else:
            return ""

    cpdef bint is_variable_binary(self, int index):
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,0)
            sage: p.is_variable_binary(0)
            True

        """
        return glp_get_col_kind(self.lp, index + 1) == GLP_BV

    cpdef bint is_variable_integer(self, int index):
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)
            True
        """
        return glp_get_col_kind(self.lp, index + 1) == GLP_IV

    cpdef bint is_variable_continuous(self, int index):
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

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

        """
        return glp_get_col_kind(self.lp, index + 1) == GLP_CV

    cpdef bint is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """

        return glp_get_obj_dir(self.lp) == GLP_MAX

    cpdef variable_upper_bound(self, int index, value = False):
        """
        Return or define the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

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

        :trac:`14581`::

            sage: P = MixedIntegerLinearProgram(solver="GLPK")
            sage: x = P["x"]
            sage: P.set_max(x, 0)
            sage: P.get_max(x)
            0.0
        """
        cdef double x
        cdef double min

        if value is False:
            x = glp_get_col_ub(self.lp, index +1)
            if x == DBL_MAX:
                return None
            else:
                return x
        else:
            min = glp_get_col_lb(self.lp, index + 1)

            if value is None and min == -DBL_MAX:
                glp_set_col_bnds(self.lp, index + 1, GLP_FR, 0, 0)

            elif value is None:
                glp_set_col_bnds(self.lp, index + 1, GLP_LO, min, 0)

            elif min == -DBL_MAX:
                glp_set_col_bnds(self.lp, index + 1, GLP_UP, 0, value)

            elif min == value:
                glp_set_col_bnds(self.lp, index + 1, GLP_FX,  value, value)

            else:
                glp_set_col_bnds(self.lp, index + 1, GLP_DB, min, value)

    cpdef variable_lower_bound(self, int index, value = False):
        """
        Return or define the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

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

        :trac:`14581`::

            sage: P = MixedIntegerLinearProgram(solver="GLPK")
            sage: x = P["x"]
            sage: P.set_min(x, 5)
            sage: P.set_min(x, 0)
            sage: P.get_min(x)
            0.0
        """
        cdef double x
        cdef double max

        if value is False:
            x = glp_get_col_lb(self.lp, index +1)
            if x == -DBL_MAX:
                return None
            else:
                return x
        else:
            max = glp_get_col_ub(self.lp, index + 1)

            if value is None and max == DBL_MAX:
                glp_set_col_bnds(self.lp, index + 1, GLP_FR, 0.0, 0.0)

            elif value is None:
                glp_set_col_bnds(self.lp, index + 1, GLP_UP, 0.0, max)

            elif max == DBL_MAX:
                glp_set_col_bnds(self.lp, index + 1, GLP_LO, value, 0.0)

            elif max == value:
                glp_set_col_bnds(self.lp, index + 1, GLP_FX,  value, value)

            else:
                glp_set_col_bnds(self.lp, index + 1, GLP_DB, value, max)

    cpdef write_lp(self, char * filename):
        """
        Write the problem to a .lp file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))
            Writing problem data to ...
            9 lines were written
        """
        glp_write_lp(self.lp, NULL, filename)

    cpdef write_mps(self, char * filename, int modern):
        """
        Write the problem to a .mps file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.write_mps(os.path.join(SAGE_TMP, "lp_problem.mps"), 2)
            Writing problem data to...
            17 records were written
        """
        glp_write_mps(self.lp, modern, NULL,  filename)

    cpdef GLPKBackend copy(self):
        """
        Returns a copy of self.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "GLPK")
            sage: b = p.new_variable()
            sage: p.add_constraint(b[1] + b[2] <= 6)
            sage: p.set_objective(b[1] + b[2])
            sage: copy(p).solve()
            6.0
        """
        cdef GLPKBackend p = GLPKBackend(maximization = (1 if self.is_maximization() else -1))
        p.simplex_or_intopt = self.simplex_or_intopt
        p.iocp.tm_lim = self.iocp.tm_lim
        glp_copy_prob(p.lp, self.lp, 1)
        return p


    cpdef solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter

        INPUT:

        - ``name`` (string) -- the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value.

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
        For example, both ``glp_gmi_cuts`` or ``"gmi_cuts"`` control whether
        to solve using Gomory cuts.

        Parameter **values** are specificed as strings in upper case,
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

           - whether to use the ``simplex`` or ``intopt`` routines in
             GLPK. This is controlled by using ``glp_simplex_only``,
             ``glp_intopt_only``, and ``glp_simplex_then_intopt``. The latter
             is useful to deal with a problem in GLPK where problems with no
             solution hang when using integer optimization; if you specify
             ``glp_simplex_then_intopt``, sage will try simplex first, then
             perform integer optimization only if a solution of the LP
             relaxation exists.

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

           - ``GLP_ON`` (default) or ``GLP_OFF``.

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

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.solver_parameter("timelimit", 60)
            sage: p.solver_parameter("timelimit")
            60.0

        - Don't forget the difference between ``timelimit`` and ``timelimit_intopt``

        ::

            sage: p.solver_parameter("timelimit_intopt")
            60000

        If you don't care for an integer answer, you can ask for an lp
        relaxation instead.  The default solver performs integer optimization,
        but you can switch to the standard simplex algorithm through the
        ``glp_simplex_or_intopt`` parameter.

        EXAMPLE::

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

        You can get glpk to spout all sorts of information at you.
        The default is to turn this off, but sometimes (debugging) it's very useful::

            sage: lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_then_intopt)
            sage: lp.solver_parameter(backend.glp_mir_cuts, backend.glp_on)
            sage: lp.solver_parameter(backend.glp_msg_lev_intopt, backend.glp_msg_all)
            sage: lp.solver_parameter(backend.glp_mir_cuts)
            1

        If you actually try to solve ``lp``, you will get a lot of detailed information.
        """

        if type(name) == str:
          if name == "out_frq" or name == "out_dly" or name == "tm_lim" or name == "msg_lev":
            raise ValueError("To set parameter " + name + " you must specify the solver. Append either _simplex or _intopt.")
          name = solver_parameter_names_dict[name]

        if type(value) == str: value = solver_parameter_values[value]

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
            if value is None: return self.simplex_or_intopt
            if not value in (simplex_only,intopt_only,simplex_then_intopt):
                raise MIPSolverException, "GLPK: invalid value for simplex_or_intopt; see documentation"
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


    cpdef int print_ranges(self, char * filename = NULL) except -1:
        r"""
        Print results of a sensitivity analysis

        If no filename is given as an input the results of the
        sensitivity analysis are displayed on the screen. If a
        filename is given they are written to a file.

        INPUT:

        - ``filename`` -- (optional) name of the file

        OUTPUT:

        Zero if the operations was successful otherwise nonzero.

        .. NOTE::

            This method is only effective if an optimal solution has been found
            for the lp using the simplex algorithm. In all other cases an error
            message is printed.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint(zip([0, 1], [1, 2]), None, 3)
            sage: p.set_objective([2, 5])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: p.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: p.print_ranges()
            glp_print_ranges: optimal basic solution required
            1
            sage: p.solve()
            0
            sage: p.print_ranges()
            Write sensitivity analysis report to...
            GLPK ... - SENSITIVITY ANALYSIS REPORT                                                                         Page   1
            <BLANKLINE>
            Problem:
            Objective:  7.5 (MAXimum)
            <BLANKLINE>
               No. Row name     St      Activity         Slack   Lower bound       Activity      Obj coef  Obj value at Limiting
                                                      Marginal   Upper bound          range         range   break point variable
            ------ ------------ -- ------------- ------------- -------------  ------------- ------------- ------------- ------------
                 1              NU       3.00000        .               -Inf         .           -2.50000        .
                                                       2.50000       3.00000           +Inf          +Inf          +Inf
            <BLANKLINE>
            GLPK ... - SENSITIVITY ANALYSIS REPORT                                                                         Page   2
            <BLANKLINE>
            Problem:
            Objective:  7.5 (MAXimum)
            <BLANKLINE>
               No. Column name  St      Activity      Obj coef   Lower bound       Activity      Obj coef  Obj value at Limiting
                                                      Marginal   Upper bound          range         range   break point variable
            ------ ------------ -- ------------- ------------- -------------  ------------- ------------- ------------- ------------
                 1              NL        .            2.00000        .                -Inf          -Inf          +Inf
                                                       -.50000          +Inf        3.00000       2.50000       6.00000
            <BLANKLINE>
                 2              BS       1.50000       5.00000        .                -Inf       4.00000       6.00000
                                                        .               +Inf        1.50000          +Inf          +Inf
            <BLANKLINE>
            End of report
            <BLANKLINE>
            0
        """

        from sage.misc.all import SAGE_TMP

        if filename == NULL:
            fname = SAGE_TMP+"/ranges.tmp"
        else:
            fname = filename

        res = glp_print_ranges(self.lp, 0, 0, 0, fname)

        if filename == NULL:
            if res == 0:
                with open(fname) as f:
                    for line in f:
                        print line,
                print

        return res

    cpdef double get_row_dual(self, int variable):
        r"""
        Returns the dual value of a constraint.

        The dual value of the ith row is also the value of the ith variable
        of the dual problem.

        The dual value of a constraint is the shadow price of the constraint.
        The shadow price is the amount by which the objective value will change
        if the constraints bounds change by one unit under the precondition
        that the basis remains the same.

        INPUT:

        - ``variable`` -- The number of the constraint

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.
           If the simplex algorithm has not been used for solving 0.0 will
           be returned.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "GLPK")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(zip([0, 1, 2], [8, 6, 1]), None, 48)
            sage: lp.add_linear_constraint(zip([0, 1, 2], [4, 2, 1.5]), None, 20)
            sage: lp.add_linear_constraint(zip([0, 1, 2], [2, 1.5, 0.5]), None, 8)
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

    cpdef double get_col_dual(self, int variable):
        """
        Returns the dual value (reduced cost) of a variable

        The dual value is the reduced cost of a variable.
        The reduced cost is the amount by which the objective coefficient
        of a non basic variable has to change to become a basic variable.

        INPUT:

        - ``variable`` -- The number of the variable

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.
           If the simplex algorithm has not been used for solving just a
           0.0 will be returned.


        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(3)
            2
            sage: p.add_linear_constraint(zip([0, 1, 2], [8, 6, 1]), None, 48)
            sage: p.add_linear_constraint(zip([0, 1, 2], [4, 2, 1.5]), None, 20)
            sage: p.add_linear_constraint(zip([0, 1, 2], [2, 1.5, 0.5]), None, 8)
            sage: p.set_objective([60, 30, 20])
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: p.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: p.solve()
            0
            sage: p.get_col_dual(1)
            -5.0

        """
        if self.simplex_or_intopt == simplex_only:
            return glp_get_col_dual(self.lp, variable+1)
        else:
            return 0.0


    def __dealloc__(self):
        """
        Destructor
        """
        glp_delete_prob(self.lp)
        sage_free(self.iocp)
        sage_free(self.smcp)

# parameter names

cdef enum solver_parameter_names:
  timelimit_seconds, timelimit_simplex, timelimit_intopt, simplex_or_intopt, msg_lev_simplex,
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
  'out_frq_simplex': out_frq_simplex, 'output_frequency_intopt': out_frq_simplex,
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

cdef enum more_parameter_values:
  simplex_only, simplex_then_intopt, intopt_only

glp_simplex_only = simplex_only
glp_simplex_then_intopt = simplex_then_intopt
glp_intopt_only = intopt_only

# dictionaries for those who prefer to use strings

solver_parameter_values = {

  'simplex_only': simplex_only,
  'simplex_then_intopt': simplex_then_intopt,
  'intopt_only': intopt_only,

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
