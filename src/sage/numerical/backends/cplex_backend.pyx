"""
CPLEX Backend

AUTHORS:

- Nathann Cohen (2010-10): initial implementation

"""

##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

include "cysignals/memory.pxi"

from sage.numerical.mip import MIPSolverException

cdef class CPLEXBackend(GenericBackend):

    """
    MIP Backend that uses the CPLEX solver.

    TESTS:

    General backend testsuite::

        sage: p = MixedIntegerLinearProgram(solver="CPLEX")                 # optional - CPLEX
        sage: TestSuite(p.get_backend()).run(skip="_test_pickling")         # optional - CPLEX

    """

    def __cinit__(self, maximization = True):
        """
        Constructor

        EXAMPLE::

        sage: from sage.numerical.mip import MixedIntegerLinearProgram      # optional - CPLEX
        sage: p = MixedIntegerLinearProgram(solver="CPLEX")                 # optional - CPLEX
        """

        cdef int status
        self.env = CPXopenCPLEX (&status)
        check(status)

        cdef char * tmp = ""
        self.lp = CPXcreateprob (self.env, &status, tmp);
        check(status)

        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

        self.obj_constant_term = 0.0
        self._logfilename = ''

    cpdef int add_variable(self, lower_bound=0.0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0.0, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

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

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")                  # optional - CPLEX
            sage: p.ncols()                                         # optional - CPLEX
            0
            sage: p.add_variable()                                  # optional - CPLEX
            0
            sage: p.ncols()                                         # optional - CPLEX
            1
            sage: p.add_variable(binary=True)                       # optional - CPLEX
            1
            sage: p.add_variable(lower_bound=-2.0, integer=True)    # optional - CPLEX
            2
            sage: p.add_variable(continuous=True, integer=True)     # optional - CPLEX
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.add_variable(name='x',obj=1.0)                  # optional - CPLEX
            3
            sage: p.col_name(3)                                     # optional - CPLEX
            'x'
            sage: p.objective_coefficient(3)                        # optional - CPLEX
            1.0

        """
        cdef char * c_name
        cdef double c_coeff = obj
        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if  vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        cdef int status
        status = CPXnewcols(self.env, self.lp, 1, NULL, NULL, NULL, NULL, NULL)
        check(status)

        cdef int n
        n = CPXgetnumcols(self.env, self.lp) - 1

        if lower_bound != 0.0:
            self.variable_lower_bound(n, lower_bound)
        if upper_bound is not None:
            self.variable_upper_bound(n, upper_bound)

        if binary:
            self.set_variable_type(n,0)
        elif integer:
            self.set_variable_type(n,1)

        if name is not None:
            c_name = name
            status = CPXchgcolname(self.env, self.lp, 1, &n, &c_name)
            check(status)

        if c_coeff:
            status = CPXchgobj(self.env, self.lp, 1, &n, &c_coeff)
            check(status)

        return n

    cpdef int add_variables(self, int number, lower_bound=0.0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0.0, names=None) except -1:
        """
        Add ``number`` new variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

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

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")                         # optional - CPLEX
            sage: p.ncols()                                                # optional - CPLEX
            0
            sage: p.add_variables(5)                                       # optional - CPLEX
            4
            sage: p.ncols()                                                # optional - CPLEX
            5
            sage: p.add_variables(2, lower_bound=-2.0, integer=True, obj=42.0, names=['a','b']) # optional - CPLEX
            6

        TESTS:

        Check that arguments are used::

            sage: p.col_bounds(5) # tol 1e-8, optional - CPLEX
            (-2.0, None)
            sage: p.is_variable_integer(5)   # optional - CPLEX
            True
            sage: p.col_name(5)              # optional - CPLEX
            'a'
            sage: p.objective_coefficient(5) # tol 1e-8, optional - CPLEX
            42.0
        """
        cdef char * c_name
        cdef double c_coeff = obj
        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if  vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        cdef int numcols_before
        numcols_before = CPXgetnumcols(self.env, self.lp)

        cdef int status
        status = CPXnewcols(self.env, self.lp, number, NULL, NULL, NULL, NULL, NULL)
        check(status)

        cdef int n
        n = CPXgetnumcols(self.env, self.lp) - 1

        cdef int i, j

        for 0<= i < number:
            j = numcols_before + i

            if lower_bound != 0.0:
                self.variable_lower_bound(j, lower_bound)
            if upper_bound is not None:
                self.variable_upper_bound(j, upper_bound)

            if binary:
                self.set_variable_type(j, 0)
            elif integer:
                self.set_variable_type(j, 1)

            if names:
                c_name = names[i]
                status = CPXchgcolname(self.env, self.lp, 1, &j, &c_name)
                check(status)

            if c_coeff:
                status = CPXchgobj(self.env, self.lp, 1, &j, &c_coeff)
                check(status)

        return n

    cpdef set_variable_type(self, int variable, int vtype):
        r"""
        Sets the type of a variable

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``vtype`` (integer) :

            *  1  Integer
            *  0  Binary
            * -1 Real

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")   # optional - CPLEX
            sage: p.ncols()                                        # optional - CPLEX
            0
            sage: p.add_variable()                                  # optional - CPLEX
            0
            sage: p.set_variable_type(0,1)                          # optional - CPLEX
            sage: p.is_variable_integer(0)                          # optional - CPLEX
            True
        """

        cdef int status

        cdef char type
        if vtype == 1:
            type = 'I'
        elif vtype == 0:
            type = 'B'
        else:
            type = 'C'

        status = CPXchgctype(self.env, self.lp, 1, &variable, &type)
        check(status)

    cpdef set_sense(self, int sense):
        r"""
        Sets the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.is_maximization()                              # optional - CPLEX
            True
            sage: p.set_sense(-1)                              # optional - CPLEX
            sage: p.is_maximization()                              # optional - CPLEX
            False
        """

        CPXchgobjsen(self.env, self.lp, -sense)

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient or ``None`` for
          reading (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")       # optional -- CPLEX
            sage: p.add_variable()                      # optional -- CPLEX
            0
            sage: p.objective_coefficient(0)            # optional -- CPLEX
            0.0
            sage: p.objective_coefficient(0,2)          # optional -- CPLEX
            sage: p.objective_coefficient(0)            # optional -- CPLEX
            2.0
        """

        cdef int status
        cdef double value

        if coeff is None:
            status = CPXgetobj(self.env, self.lp, &value, variable, variable)
            check(status)
            return value

        else:
            value = coeff
            status = CPXchgobj(self.env, self.lp, 1, &variable, &value)
            check(status)

    cpdef problem_name(self, char * name = NULL):
        r"""
        Returns or defines the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")   # optional - CPLEX
            sage: p.problem_name("There once was a french fry") # optional - CPLEX
            sage: print(p.problem_name())                       # optional - CPLEX
            There once was a french fry
        """

        cdef int status
        cdef int zero
        cdef char * n
        if name == NULL:

            n = <char*> sig_malloc(500*sizeof(char))
            status = CPXgetprobname(self.env, self.lp, n, 500, &zero)
            check(status)
            s = str(n)
            sig_free(n)
            return s


        else:
            status = CPXchgprobname(self.env, self.lp, name)
            check(status)


    cpdef set_objective(self, list coeff, d = 0.0):
        r"""
        Sets the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        - ``d`` (double) -- the constant term in the linear function (set to `0` by default)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")    # optional - CPLEX
            sage: p.add_variables(5)                                 # optional - CPLEX
            4
            sage: p.set_objective([1, 1, 2, 1, 3])                   # optional - CPLEX
            sage: map(lambda x :p.objective_coefficient(x), range(5))  # optional - CPLEX
            [1.0, 1.0, 2.0, 1.0, 3.0]

        Constants in the objective function are respected::

            sage: p = MixedIntegerLinearProgram(solver='CPLEX') # optional - CPLEX
            sage: var = p.new_variable(nonnegative=True)        # optional - CPLEX
            sage: x,y = var[0], var[1]                          # optional - CPLEX
            sage: p.add_constraint(2*x + 3*y, max = 6)          # optional - CPLEX
            sage: p.add_constraint(3*x + 2*y, max = 6)          # optional - CPLEX
            sage: p.set_objective(x + y + 7)                    # optional - CPLEX
            sage: p.set_integer(x); p.set_integer(y)            # optional - CPLEX
            sage: p.solve()                                     # optional - CPLEX
            9.0
        """

        cdef int status
        cdef int n = self.ncols()
        cdef double * c_coeff = <double *> sig_malloc(n * sizeof(double))
        cdef int * c_indices = <int *> sig_malloc(n * sizeof(int))

        for i,v in enumerate(coeff):
            c_coeff[i] = v
            c_indices[i] = i

        status = CPXchgobj(self.env, self.lp, n, c_indices, c_coeff)
        check(status)

        sig_free(c_coeff)
        sig_free(c_indices)

        self.obj_constant_term = d


    cpdef set_verbosity(self, int level):
        r"""
        Sets the log (verbosity) level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")   # optional - CPLEX
            sage: p.set_verbosity(2)                                # optional - CPLEX

        """

        cdef int status
        if level == 0:
            status = CPXsetintparam (self.env, CPX_PARAM_SCRIND, 0)
            check(status)
        else:
            status = CPXsetintparam (self.env, CPX_PARAM_SCRIND, CPX_ON)
            check(status)

    cpdef remove_constraint(self, int i):
        r"""
        Remove a constraint from self.

        INPUT:

        - ``i`` -- index of the constraint to remove

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver='CPLEX')# optional - CPLEX
            sage: var = p.new_variable(nonnegative=True)       # optional - CPLEX
            sage: x,y = var[0], var[1]                         # optional - CPLEX
            sage: p.add_constraint(2*x + 3*y, max = 6)         # optional - CPLEX
            sage: p.add_constraint(3*x + 2*y, max = 6)         # optional - CPLEX
            sage: p.set_objective(x + y + 7)                   # optional - CPLEX
            sage: p.set_integer(x); p.set_integer(y)           # optional - CPLEX
            sage: p.solve()                                    # optional - CPLEX
            9.0
            sage: p.remove_constraint(0)                       # optional - CPLEX
            sage: p.solve()                                    # optional - CPLEX
            10.0
            sage: p.get_values([x,y])                          # optional - CPLEX
            [0.0, 3.0]
        """
        cdef int status
        status = CPXdelrows(self.env, self.lp, i, i)
        check(status)

    cpdef add_linear_constraints(self, int number, lower_bound, upper_bound, names = None):
        """
        Add ``'number`` linear constraints.

        INPUT:

        - ``number`` (integer) -- the number of constraints to add.

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``names`` - an optional list of names (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")       # optional - CPLEX
            sage: p.add_variables(5)                     # optional - CPLEX
            4
            sage: p.add_linear_constraints(5, None, 2)   # optional - CPLEX
            sage: p.row(4)                               # optional - CPLEX
            ([], [])
            sage: p.row_bounds(4)                        # optional - CPLEX
            (None, 2.0)
            sage: p.add_linear_constraints(2, None, 2, names=['foo','bar']) # optional - CPLEX
        """
        if lower_bound is None and upper_bound is None:
            raise ValueError("At least one of 'upper_bound' or 'lower_bound' must be set.")

        cdef int status
        cdef char * sense = <char *> sig_malloc(number * sizeof(char))
        cdef double * bound = <double *> sig_malloc(number * sizeof(double))
        cdef double * rng = NULL
        cdef int i
        cdef char ** c_names = <char **> sig_malloc(number * sizeof(char *))

        if upper_bound == lower_bound:
            sense[0] = 'E'
            bound[0] = lower_bound

        elif upper_bound is not None and lower_bound is not None:
            if  upper_bound < lower_bound:
                raise ValueError("The upper bound must be at least equal to the lower bound !")

            rng = <double *> sig_malloc(number * sizeof(double))

            sense[0] = 'R'
            bound[0] = lower_bound
            rng[0] = upper_bound - lower_bound

        elif upper_bound is not None:
            sense[0] = 'L'
            bound[0] = upper_bound

        elif lower_bound is not None:
            sense[0] = 'G'
            bound[0] = lower_bound

        if names:
            c_names[0] = names[0]

        for 1<= i <number:
            sense[i] = sense[0]
            bound[i] = bound[0]
            if rng != NULL:
                rng[i] = rng[0]
            if names:
                c_names[i] = names[i]

        status = CPXnewrows(self.env, self.lp, number, bound, sense, rng, c_names if names else NULL)

        sig_free(sense)
        sig_free(bound)
        sig_free(c_names)
        check(status)

    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name = None):
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

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")                             # optional - CPLEX
            sage: p.add_variables(5)                                           # optional - CPLEX
            4
            sage: p.add_linear_constraint( zip(range(5), range(5)), 2.0, 2.0)  # optional - CPLEX
            sage: p.row(0)                                                     # optional - CPLEX
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                                              # optional - CPLEX
            (2.0, 2.0)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1.0, 1.0, name='foo') # optional - CPLEX
            sage: p.row_name(1)                                                           # optional - CPLEX
            'foo'

        """
        if lower_bound is None and upper_bound is None:
            raise ValueError("At least one of 'upper_bound' or 'lower_bound' must be set.")

        cdef int status
        cdef int i, j
        cdef int n = len(coefficients)
        cdef int nrows = self.nrows()
        cdef char sense

        cdef char * c_name

        cdef double * c_coeff
        cdef int * c_indices
        cdef int * c_row
        cdef double bound
        cdef double rng
        cdef double c

        c_coeff = <double *> sig_malloc(n * sizeof(double))
        c_indices = <int *> sig_malloc(n * sizeof(int))
        c_row = <int *> sig_malloc(n * sizeof(int))

        for i, (j, c) in enumerate(coefficients):
            c_coeff[i] = c
            c_indices[i] = j
            c_row[i] = nrows

        if upper_bound is None and lower_bound is None:
            pass

        elif upper_bound == lower_bound:
            sense = 'E'
            bound = lower_bound

        elif upper_bound is not None and lower_bound is not None:
            if  upper_bound < lower_bound:
                raise ValueError("The upper bound must be at least equal to the lower bound !")

            sense = 'R'
            bound = lower_bound
            rng = upper_bound - lower_bound

        elif upper_bound is not None:
            sense = 'L'
            bound = upper_bound

        elif lower_bound is not None:
            sense = 'G'
            bound = lower_bound

        if name:
            c_name = name

        status = CPXnewrows(self.env, self.lp, 1, &bound, &sense, &rng, NULL if (name is None) else &c_name)

        check(status)
        status = CPXchgcoeflist(self.env, self.lp, n, c_row, c_indices, c_coeff)
        check(status)

        # Free memory
        sig_free(c_coeff)
        sig_free(c_indices)
        sig_free(c_row)

    cpdef row(self, int index):
        r"""
        Returns a row

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_linear_constraint`` method.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variables(5)                               # optional - CPLEX
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)       # optional - CPLEX
            sage: p.row(0)                                     # optional - CPLEX
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                              # optional - CPLEX
            (2.0, 2.0)
        """

        cdef int status
        cdef int n,i
        cdef int zero
        cdef list indices = []
        cdef list values = []

        cdef double * c_coeff = <double *> sig_malloc((self.ncols()+10) * sizeof(double))
        cdef int * c_indices = <int *> sig_malloc((self.ncols()+10) * sizeof(int))

        status = CPXgetrows(self.env, self.lp, &n, &zero, c_indices, c_coeff, self.ncols()+3, &zero, index, index)

        check(status)

        for 0<= i<n:
            indices.append(c_indices[i])
            values.append(c_coeff[i])

        sig_free(c_coeff)
        sig_free(c_indices)

        return (indices, values)

    cpdef row_bounds(self, int index):
        r"""
        Returns the bounds of a specific constraint.

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the constraint is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variables(5)                               # optional - CPLEX
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)       # optional - CPLEX
            sage: p.row(0)                                     # optional - CPLEX
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                              # optional - CPLEX
            (2.0, 2.0)
        """

        cdef int status
        cdef double rng
        cdef double value
        status = CPXgetrhs(self.env, self.lp, &value, index, index)
        check(status)

        cdef char direction
        status = CPXgetsense(self.env, self.lp, &direction, index, index)
        check(status)

        if direction == 'L':
            return (None, value)
        elif direction == 'G':
            return (value, None)
        elif direction == 'E':
            return (value, value)
        elif direction == 'R':
            status = CPXgetrngval(self.env, self.lp, &rng, index, index)
            check(status)
            return (value, value + rng)

    cpdef col_bounds(self, int index):
        r"""
        Returns the bounds of a specific variable.

        INPUT:

        - ``index`` (integer) -- the variable's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the variable is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variable()                                 # optional - CPLEX
            0
            sage: p.col_bounds(0)                              # optional - CPLEX
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                         # optional - CPLEX
            sage: p.col_bounds(0)                              # optional - CPLEX
            (0.0, 5.0)
        """

        cdef int status
        cdef double ub
        cdef double lb

        status = CPXgetub(self.env, self.lp, &ub, index, index)
        check(status)

        status = CPXgetlb(self.env, self.lp, &lb, index, index)
        check(status)

        return (None if lb <= -int(CPX_INFBOUND) else lb,
                None if ub >= +int(CPX_INFBOUND) else ub)

    cpdef add_col(self, list indices, list coeffs):
        r"""
        Adds a column.

        INPUT:

        - ``indices`` (list of integers) -- this list contains the
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

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")                  # optional - CPLEX
            sage: p.ncols()                                       # optional - CPLEX
            0
            sage: p.nrows()                                       # optional - CPLEX
            0
            sage: p.add_linear_constraints(5, 0, None)                      # optional - CPLEX
            sage: p.add_col(range(5), range(5))                    # optional - CPLEX
            sage: p.nrows()                                       # optional - CPLEX
            5
        """

        cdef int status
        cdef int i
        cdef int n = len(indices)
        cdef int ncols = self.ncols()

        status = CPXnewcols(self.env, self.lp, 1, NULL, NULL, NULL, NULL, NULL)


        check(status)

        cdef double * c_coeff = <double *> sig_malloc(n * sizeof(double))
        cdef int * c_indices = <int *> sig_malloc(n * sizeof(int))
        cdef int * c_col = <int *> sig_malloc(n * sizeof(int))

        for 0<= i < n:
            c_coeff[i] = coeffs[i]
            c_indices[i] = indices[i]
            c_col[i] = ncols


        status = CPXchgcoeflist(self.env, self.lp, n, c_indices, c_col, c_coeff)
        check(status)

        sig_free(c_coeff)
        sig_free(c_indices)
        sig_free(c_col)

    cpdef int solve(self) except -1:
        r"""
        Solves the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE:

        A simple maximization problem::

            sage: p = MixedIntegerLinearProgram(solver='CPLEX')      # optional - CPLEX
            sage: x = p.new_variable(integer=True, nonnegative=True) # optional - CPLEX
            sage: p.add_constraint(2*x[0] + 3*x[1], max = 6)         # optional - CPLEX
            sage: p.add_constraint(3*x[0] + 2*x[1], max = 6)         # optional - CPLEX
            sage: p.set_objective(x[0] + x[1] + 7)                   # optional - CPLEX
            sage: p.solve()                                          # optional - CPLEX
            9.0

        A problem without feasible solution::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")                               # optional - CPLEX
            sage: p.add_linear_constraints(5, 0, None)                           # optional - CPLEX
            sage: p.add_col(range(5), range(5))                                  # optional - CPLEX
            sage: p.objective_coefficient(0,1)                                   # optional - CPLEX
            sage: p.solve()                                                      # optional - CPLEX
            Traceback (most recent call last):
            ...
            MIPSolverException: CPLEX: The problem is infeasible or unbounded
        """
        cdef int status
        cdef int ptype
        cdef int solnmethod_p, solntype_p, pfeasind_p, dfeasind_p

        ptype = CPXgetprobtype(self.env, self.lp)

        if ptype == 1:
            status = CPXmipopt(self.env, self.lp)
        elif ptype == 0:
            status = CPXlpopt(self.env, self.lp)
        else:
            raise MIPSolverException("CPLEX: Unknown problem type")

        check(status)

        stat = CPXgetstat(self.env, self.lp)
        if stat == CPX_STAT_OPTIMAL or stat == CPXMIP_OPTIMAL:
            return 0
        elif stat == CPX_STAT_INFEASIBLE or stat == CPXMIP_INFEASIBLE:
            raise MIPSolverException("CPLEX: The problem has no feasible solution")
        elif stat == CPX_STAT_UNBOUNDED or stat == CPXMIP_UNBOUNDED:
            raise MIPSolverException("CPLEX: The problem is unbounded")
        elif stat == CPX_STAT_INForUNBD or stat == CPXMIP_INForUNBD:
            raise MIPSolverException("CPLEX: The problem is infeasible or unbounded")
        else:
            # TODO: Many more stats to be handled.
            pass

        # No exception should be raised when CPX_STAT_ABORT_... or CPXMIP_ABORT_...
        # This is so that when a time limit etc. is reached, we obtain meaningful information.

        status = CPXsolninfo(self.env, self.lp, &solnmethod_p, &solntype_p, &pfeasind_p, &dfeasind_p)
        check(status)

        if solntype_p == CPX_NO_SOLN or not pfeasind_p:
            raise MIPSolverException("CPLEX: No solution known to be primal feasible is available")

        return 0

    cpdef get_objective_value(self):
        r"""
        Returns the value of the objective function.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variables(2)                               # optional - CPLEX
            1
            sage: p.add_linear_constraint([(0,1), (1,2)], None, 3)          # optional - CPLEX
            sage: p.set_objective([2, 5])                          # optional - CPLEX
            sage: p.solve()                                        # optional - CPLEX
            0
            sage: p.get_objective_value()                          # optional - CPLEX
            7.5
            sage: p.get_variable_value(0)                          # optional - CPLEX
            0.0
            sage: p.get_variable_value(1)                          # optional - CPLEX
            1.5
        """

        cdef int status
        cdef double value
        status = CPXgetobjval (self.env, self.lp, &value)
        check(status)

        return value + self.obj_constant_term


    cpdef best_known_objective_bound(self):
        r"""
        Return the value of the currently best known bound.

        This method returns the current best upper (resp. lower) bound on the
        optimal value of the objective function in a maximization
        (resp. minimization) problem. It is equal to the output of
        :meth:`~MixedIntegerLinearProgram.get_objective_value` if the MILP found
        an optimal solution, but it can differ if it was interrupted manually or
        after a time limit (cf :meth:`MixedIntegerLinearProgram.solver_parameter`).

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

           When using the ``CPX_PARAM_POPULATELIM`` parameter to get a pool of
           solutions, the value can exceed the optimal solution value if
           ``CPLEX`` has already solved the model to optimality but continues to
           search for additional solutions.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver="CPLEX")        # optional - CPLEX
            sage: b = p.new_variable(binary=True)                      # optional - CPLEX
            sage: for u,v in graphs.CycleGraph(5).edges(labels=False): # optional - CPLEX
            ....:     p.add_constraint(b[u]+b[v]<=1)                   # optional - CPLEX
            sage: p.set_objective(p.sum(b[x] for x in range(5)))       # optional - CPLEX
            sage: p.solve()                                            # optional - CPLEX
            2.0
            sage: pb = p.get_backend()                                 # optional - CPLEX
            sage: pb.get_objective_value()                             # optional - CPLEX
            2.0
            sage: pb.best_known_objective_bound()                      # optional - CPLEX
            2.0
        """
        cdef int status
        cdef double value
        status = CPXgetbestobjval (self.env, self.lp, &value)
        check(status)

        return value + self.obj_constant_term

    cpdef get_relative_objective_gap(self):
        r"""
        Return the relative objective gap of the best known solution.

        For a minimization problem, this value is computed by
        `(\texttt{bestinteger} - \texttt{bestobjective}) / (1e-10 +
        |\texttt{bestobjective}|)`, where ``bestinteger`` is the value returned
        by :meth:`~MixedIntegerLinearProgram.get_objective_value` and
        ``bestobjective`` is the value returned by
        :meth:`MixedIntegerLinearProgram.best_known_objective_bound`. For a
        maximization problem, the value is computed by `(\texttt{bestobjective}
        - \texttt{bestinteger}) / (1e-10 + |\texttt{bestobjective}|)`.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver="CPLEX")        # optional - CPLEX
            sage: b = p.new_variable(binary=True)                      # optional - CPLEX
            sage: for u,v in graphs.CycleGraph(5).edges(labels=False): # optional - CPLEX
            ....:     p.add_constraint(b[u]+b[v]<=1)                   # optional - CPLEX
            sage: p.set_objective(p.sum(b[x] for x in range(5)))       # optional - CPLEX
            sage: p.solve()                                            # optional - CPLEX
            2.0
            sage: pb = p.get_backend()                                 # optional - CPLEX
            sage: pb.get_objective_value()                             # optional - CPLEX
            2.0
            sage: pb.get_relative_objective_gap()                      # optional - CPLEX
            0.0
        """
        cdef int status
        cdef double value
        status = CPXgetmiprelgap (self.env, self.lp, &value)
        check(status)

        return value

    cpdef get_variable_value(self, int variable):
        r"""
        Returns the value of a variable given by the solver.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX") # optional - CPLEX
            sage: p.add_variables(2)                              # optional - CPLEX
            1
            sage: p.add_linear_constraint([(0,1), (1,2)], None, 3)         # optional - CPLEX
            sage: p.set_objective([2, 5])                         # optional - CPLEX
            sage: p.solve()                                       # optional - CPLEX
            0
            sage: p.get_objective_value()                         # optional - CPLEX
            7.5
            sage: p.get_variable_value(0)                         # optional - CPLEX
            0.0
            sage: p.get_variable_value(1)                         # optional - CPLEX
            1.5
        """

        cdef int status
        cdef int zero
        cdef char ctype
        cdef double value
        status = CPXgetx(self.env, self.lp, &value, variable, variable)
        check(status)

        status = CPXgetctype(self.env, self.lp, &ctype, variable, variable)

        return value if (status == 3003 or ctype=='C') else round(value)

    cpdef int ncols(self):
        r"""
        Returns the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.ncols()                                       # optional - CPLEX
            0
            sage: p.add_variables(2)                               # optional - CPLEX
            1
            sage: p.ncols()                                       # optional - CPLEX
            2
        """

        return CPXgetnumcols(self.env, self.lp)

    cpdef int nrows(self):
        r"""
        Returns the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX") # optional - CPLEX
            sage: p.nrows()                                      # optional - CPLEX
            0
            sage: p.add_linear_constraints(2, 2, None)                     # optional - CPLEX
            sage: p.nrows()                                      # optional - CPLEX
            2
        """

        return CPXgetnumrows(self.env, self.lp)

    cpdef row_name(self, int index):
        r"""
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver       # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")                                     # optional - CPLEX
            sage: p.add_linear_constraints(1, 2, None, names=['Empty constraint 1'])   # optional - CPLEX
            sage: p.row_name(0)                                                        # optional - CPLEX
            'Empty constraint 1'
        """

        cdef int status
        cdef int zero
        cdef char * n

        n = <char *>sig_malloc(500*sizeof(char))
        status = CPXgetrowname(self.env, self.lp, &n, n, 500, &zero, index, index)
        if status == 1219:
            sig_free(n)
            return ""
        check(status)

        s = str(n)
        sig_free(n)

        return s

    cpdef col_name(self, int index):
        r"""
        Returns the ``index`` th col name.

        INPUT:

        - ``index`` (integer) -- the col's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")         # optional - CPLEX
            sage: p.add_variable(name='I am a variable')   # optional - CPLEX
            0
            sage: p.col_name(0)                            # optional - CPLEX
            'I am a variable'
        """

        cdef int status
        cdef char * n
        cdef int zero

        n = <char *>sig_malloc(500*sizeof(char))
        status = CPXgetcolname(self.env, self.lp, &n, n, 500, &zero, index, index)
        if status == 1219:
            sig_free(n)
            return ""
        check(status)

        s = str(n)
        sig_free(n)
        return s

    cpdef bint is_variable_binary(self, int index):
        r"""
        Tests whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.ncols()                                       # optional - CPLEX
            0
            sage: p.add_variable()                                 # optional - CPLEX
            0
            sage: p.set_variable_type(0,0)                         # optional - CPLEX
            sage: p.is_variable_binary(0)                          # optional - CPLEX
            True

        """

        cdef int status
        cdef char ctype

        status = CPXgetctype(self.env, self.lp, &ctype, index, index)

        # status = 3003 when the problem is a LP and not a MILP
        if status == 3003:
            return False

        check(status)

        return ctype == 'B'


    cpdef bint is_variable_integer(self, int index):
        r"""
        Tests whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.ncols()                                       # optional - CPLEX
            0
            sage: p.add_variable()                                 # optional - CPLEX
            0
            sage: p.set_variable_type(0,1)                         # optional - CPLEX
            sage: p.is_variable_integer(0)                         # optional - CPLEX
            True
        """

        cdef int status
        cdef char ctype

        status = CPXgetctype(self.env, self.lp, &ctype, index, index)

        # status = 3003 when the problem is a LP and not a MILP
        if status == 3003:
            return False

        check(status)

        return ctype == 'I'


    cpdef bint is_variable_continuous(self, int index):
        r"""
        Tests whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.ncols()                                       # optional - CPLEX
            0
            sage: p.add_variable()                                 # optional - CPLEX
            0
            sage: p.is_variable_continuous(0)                      # optional - CPLEX
            True
            sage: p.set_variable_type(0,1)                         # optional - CPLEX
            sage: p.is_variable_continuous(0)                      # optional - CPLEX
            False

        """

        cdef int status
        cdef char ctype

        status = CPXgetctype(self.env, self.lp, &ctype, index, index)

        # status = 3003 when the problem is a LP and not a MILP
        if status == 3003:
            return True

        check(status)

        return ctype == 'C'


    cpdef bint is_maximization(self):
        r"""
        Tests whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX") # optional - CPLEX
            sage: p.is_maximization()                             # optional - CPLEX
            True
            sage: p.set_sense(-1)                             # optional - CPLEX
            sage: p.is_maximization()                             # optional - CPLEX
            False
        """

        return -1 == CPXgetobjsen(self.env, self.lp)

    cpdef variable_upper_bound(self, int index, value = False):
        r"""
        Returns or defines the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variable()                                 # optional - CPLEX
            0
            sage: p.col_bounds(0)                              # optional - CPLEX
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                         # optional - CPLEX
            sage: p.col_bounds(0)                              # optional - CPLEX
            (0.0, 5.0)

        TESTS:

        :trac:`14581`::

            sage: P = MixedIntegerLinearProgram(solver="CPLEX") # optional - CPLEX
            sage: var = P.new_variable(nonnegative=False)       # optional - CPLEX
            sage: x = var["x"]                                  # optional - CPLEX
            sage: P.set_max(x, 0)                               # optional - CPLEX
            sage: P.get_max(x)                                  # optional - CPLEX
            0.0
        """
        cdef int status
        cdef double ub
        cdef char x
        cdef double c_value

        if value is False:
            status = CPXgetub(self.env, self.lp, &ub, index, index)
            check(status)

            return ub if ub < int(CPX_INFBOUND) else None

        else:
            x = 'U'
            c_value = value if value is not None else +CPX_INFBOUND
            status = CPXchgbds(self.env, self.lp, 1, &index, &x, &c_value)
            check(status)

    cpdef variable_lower_bound(self, int index, value = False):
        r"""
        Returns or defines the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variable()                                 # optional - CPLEX
            0
            sage: p.col_bounds(0)                              # optional - CPLEX
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)                         # optional - CPLEX
            sage: p.col_bounds(0)                              # optional - CPLEX
            (5.0, None)

        TESTS:

        :trac:`14581`::

            sage: P = MixedIntegerLinearProgram(solver="CPLEX") # optional - CPLEX
            sage: var = P.new_variable(nonnegative=False)       # optional - CPLEX
            sage: x = var["x"]                                  # optional - CPLEX
            sage: P.set_min(x, 5)                               # optional - CPLEX
            sage: P.set_min(x, 0)                               # optional - CPLEX
            sage: P.get_min(x)                                  # optional - CPLEX
            0.0
        """
        cdef int status
        cdef double lb
        cdef char x
        cdef double c_value

        if value is False:
            status = CPXgetlb(self.env, self.lp, &lb, index, index)
            check(status)
            return None if lb <= int(-CPX_INFBOUND) else lb

        else:
            x = 'L'
            c_value = value if value is not None else -CPX_INFBOUND
            status = CPXchgbds(self.env, self.lp, 1, &index, &x, &c_value)
            check(status)

    cpdef write_lp(self, char * filename):
        r"""
        Writes the problem to a .lp file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variables(2)                               # optional - CPLEX
            1
            sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3)          # optional - CPLEX
            sage: p.set_objective([2, 5])                          # optional - CPLEX
            sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))            # optional - CPLEX
        """

        cdef int status
        cdef char * ext = "LP"
        status = CPXwriteprob(self.env, self.lp, filename, ext)
        check(status)

    cpdef write_mps(self, char * filename, int modern):
        r"""
        Writes the problem to a .mps file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variables(2)                               # optional - CPLEX
            1
            sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3)          # optional - CPLEX
            sage: p.set_objective([2, 5])                          # optional - CPLEX
            sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))            # optional - CPLEX
        """

        cdef int status
        cdef char * ext = "MPS"
        status = CPXwriteprob(self.env, self.lp, filename, ext)
        check(status)

    cpdef __copy__(self):
        r"""
        Returns a copy of self.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = MixedIntegerLinearProgram(solver = "CPLEX")        # optional - CPLEX
            sage: b = p.new_variable(nonnegative=True)         # optional - CPLEX
            sage: p.add_constraint(b[1] + b[2] <= 6)           # optional - CPLEX
            sage: p.set_objective(b[1] + b[2])                 # optional - CPLEX
            sage: copy(p).solve()                              # optional - CPLEX
            6.0
        """
        cdef CPLEXBackend p = type(self)()

        p.lp = CPXcloneprob(p.env, self.lp, &status)
        check(status)

        p._mixed = self._mixed
        p.current_sol = self.current_sol

        status = CPXchgprobtype(p.env, p.lp, CPXgetprobtype(self.env, self.lp))
        check(status)

        return p

    cpdef solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter

        INPUT:

        - ``name`` (string) -- the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value.

        .. NOTE::

           The list of available parameters is available at
           :meth:`sage.numerical.mip.MixedIntegerLinearProgram.solver_parameter`

        EXAMPLE:

        Set a computation time limit::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")         # optional - CPLEX
            sage: p.solver_parameter("timelimit", 60)      # optional - CPLEX
            sage: p.solver_parameter("timelimit")          # optional - CPLEX
            60.0

        Set the logfile (no log file by default)::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")           # optional - CPLEX
            sage: p.solver_parameter("logfile")              # optional - CPLEX
            ''
            sage: p.solver_parameter("logfile", '/dev/null') # optional - CPLEX
            sage: p.solver_parameter("logfile")              # optional - CPLEX
            '/dev/null'

        Disable the logfile::

            sage: from sage.numerical.backends.generic_backend import get_solver # optional - CPLEX
            sage: p = get_solver(solver = "CPLEX")         # optional - CPLEX
            sage: p.solver_parameter("logfile", '')        # optional - CPLEX
            sage: p.solver_parameter("logfile")            # optional - CPLEX
            ''

        TEST:

        Print the logfile's content (through :class:`MixedIntegerLinearProgram`)::

            sage: filename = tmp_filename(ext='.txt')                  # optional - CPLEX
            sage: p = MixedIntegerLinearProgram(solver="CPLEX")        # optional - CPLEX
            sage: b = p.new_variable(binary=True)                      # optional - CPLEX
            sage: for u,v in graphs.CycleGraph(5).edges(labels=False): # optional - CPLEX
            ....:     p.add_constraint(b[u]+b[v]<=1)                   # optional - CPLEX
            sage: p.set_objective(p.sum(b[x] for x in range(5)))       # optional - CPLEX
            sage: p.solver_parameter("logfile", filename)              # optional - CPLEX
            sage: p.solve()                                            # optional - CPLEX
            2.0
            sage: with open(filename,'r') as f:                        # optional - CPLEX
            ....:     print(f.read())                                  # optional - CPLEX
            Found incumbent of value ...
            Reduced MIP has 5 rows, 5 columns, and 10 nonzeros.
            ...
            Elapsed time = ... sec. (... ticks, tree = ... MB, solutions = ...)

        Bad name::

            sage: p.solver_parameter("Heyyyyyyyyyyy Joeeeeeeeee !!",-4) # optional - CPLEX
            Traceback (most recent call last):
            ...
            ValueError: This parameter is not available.
        """
        cdef int intv
        cdef double doublev
        cdef char * strv

        # Specific action for log file
        cdef FILE *ff
        if name.lower() == "logfile":
            if value is None: # Return logfile name
                return self._logfilename
            elif not value:   # Close current logfile and disable logs
                check( CPXsetlogfile(self.env, NULL) )
                self._logfilename = ''
            else:             # Set log file to logfilename
                ff = fopen(value, "a")
                if not ff:
                    raise ValueError("Unable to append file {}.".format(value))
                check( CPXsetlogfile(self.env, ff) )
                self._logfilename = value
            return


        # If the name has to be translated to a CPLEX parameter ID
        if name == "timelimit":
            name = "CPX_PARAM_TILIM"

        cdef int paramid = parameters.get(name,-1)

        if paramid == -1:
            raise ValueError("This parameter is not available.")

        # Type of the parameter. Can be INT (1), Double(2) or String(3)
        cdef int paramtype
        check(CPXgetparamtype(self.env, paramid, &paramtype))

        if value is None:
            if paramtype == 1:
                check(CPXgetintparam(self.env, paramid, &intv))
                return intv
            elif paramtype == 2:
                check(CPXgetdblparam(self.env, paramid, &doublev))
                return doublev
            else:
                strv = <char *>sig_malloc(500*sizeof(char))
                status = CPXgetstrparam(self.env, paramid, strv)
                s = str(strv)
                sig_free(strv)
                check(status)
                return s
        else:
            if paramtype == 1:
                check(CPXsetintparam(self.env, paramid, value))
            elif paramtype == 2:
                check(CPXsetdblparam(self.env, paramid, value))
            else:
                check(CPXsetstrparam(self.env, paramid, value))

    def __dealloc__(self):
        r"""
        Destructor for the class
        """
        cdef int status

        if self.lp != NULL:
            status = CPXfreeprob(self.env, &self.lp)

        if self.env != NULL:
            status = CPXcloseCPLEX(&self.env)

cdef check(number):
    r"""
    Given a number, raises the corresponding exception or does nothing
    if ``number == 0``.

    - ``number`` (integer) -- number corresponding to the error. If
      this number is unknown, the message contained in the raised
      exception will mention it.
    """

    # Below 1000 are 0 (no error), and some quality reports (but the
    # ERR* codes are above 1000)
    if number > 1000:
        from sage.numerical.mip import MIPSolverException
        default = "Error reported by the solver (unknown error number : "+str(number)+")"
        raise MIPSolverException("CPLEX: "+errors.get(number,default))

# Error codes
#
# Todo : when common error codes are returned, rewrite the message to
# be more meaningful

cdef dict errors

errors = {
    1001 : "CPXERR_NO_MEMORY",
    1002 : "CPXERR_NO_ENVIRONMENT",
    1003 : "CPXERR_BAD_ARGUMENT",
    1004 : "CPXERR_NULL_POINTER",
    1006 : "CPXERR_CALLBACK",
    1009 : "CPXERR_NO_PROBLEM",
    1012 : "CPXERR_LIMITS_TOO_BIG",
    1013 : "CPXERR_BAD_PARAM_NUM",
    1014 : "CPXERR_PARAM_TOO_SMALL",
    1015 : "CPXERR_PARAM_TOO_BIG",
    1016 : "CPXERR_RESTRICTED_VERSION",
    1017 : "CPXERR_NOT_FOR_MIP",
    1018 : "CPXERR_NOT_FOR_QP",
    1019 : "CPXERR_CHILD_OF_CHILD",
    1020 : "CPXERR_TOO_MANY_THREADS",
    1021 : "CPXERR_CANT_CLOSE_CHILD",
    1022 : "CPXERR_BAD_PROB_TYPE",
    1023 : "CPXERR_NOT_ONE_PROBLEM",
    1024 : "CPXERR_NOT_MILPCLASS",
    1026 : "CPXERR_STR_PARAM_TOO_LONG",
    1027 : "CPXERR_DECOMPRESSION",
    1028 : "CPXERR_BAD_PARAM_NAME",
    1029 : "CPXERR_NOT_MIQPCLASS",
    1031 : "CPXERR_NOT_FOR_QCP",
    1051 : "CPXERR_MSG_NO_CHANNEL",
    1052 : "CPXERR_MSG_NO_FILEPTR",
    1053 : "CPXERR_MSG_NO_FUNCTION",
    1101 : "CPXERR_PRESLV_INForUNBD",
    1103 : "CPXERR_PRESLV_NO_PROB",
    1106 : "CPXERR_PRESLV_ABORT",
    1107 : "CPXERR_PRESLV_BASIS_MEM",
    1108 : "CPXERR_PRESLV_COPYSOS",
    1109 : "CPXERR_PRESLV_COPYORDER",
    1110 : "CPXERR_PRESLV_SOLN_MIP",
    1111 : "CPXERR_PRESLV_SOLN_QP",
    1112 : "CPXERR_PRESLV_START_LP",
    1114 : "CPXERR_PRESLV_FAIL_BASIS",
    1115 : "CPXERR_PRESLV_NO_BASIS",
    1117 : "CPXERR_PRESLV_INF",
    1118 : "CPXERR_PRESLV_UNBD",
    1119 : "CPXERR_PRESLV_DUAL",
    1120 : "CPXERR_PRESLV_UNCRUSHFORM",
    1121 : "CPXERR_PRESLV_CRUSHFORM",
    1122 : "CPXERR_PRESLV_BAD_PARAM",
    1123 : "CPXERR_PRESLV_TIME_LIM",
    1200 : "CPXERR_INDEX_RANGE",
    1201 : "CPXERR_COL_INDEX_RANGE",
    1203 : "CPXERR_ROW_INDEX_RANGE",
    1205 : "CPXERR_INDEX_RANGE_LOW",
    1206 : "CPXERR_INDEX_RANGE_HIGH",
    1207 : "CPXERR_NEGATIVE_SURPLUS",
    1208 : "CPXERR_ARRAY_TOO_LONG",
    1209 : "CPXERR_NAME_CREATION",
    1210 : "CPXERR_NAME_NOT_FOUND",
    1211 : "CPXERR_NO_RHS_IN_OBJ",
    1215 : "CPXERR_BAD_SENSE",
    1216 : "CPXERR_NO_RNGVAL",
    1217 : "CPXERR_NO_SOLN",
    1219 : "CPXERR_NO_NAMES",
    1221 : "CPXERR_NOT_FIXED",
    1222 : "CPXERR_DUP_ENTRY",
    1223 : "CPXERR_NO_BARRIER_SOLN",
    1224 : "CPXERR_NULL_NAME",
    1225 : "CPXERR_NAN",
    1226 : "CPXERR_ARRAY_NOT_ASCENDING",
    1227 : "CPXERR_COUNT_RANGE",
    1228 : "CPXERR_COUNT_OVERLAP",
    1229 : "CPXERR_BAD_LUB",
    1230 : "CPXERR_NODE_INDEX_RANGE",
    1231 : "CPXERR_ARC_INDEX_RANGE",
    1232 : "CPXERR_NO_DUAL_SOLN",
    1233 : "CPXERR_DBL_MAX",
    1234 : "CPXERR_THREAD_FAILED",
    1251 : "CPXERR_INDEX_NOT_BASIC",
    1252 : "CPXERR_NEED_OPT_SOLN",
    1253 : "CPXERR_BAD_STATUS",
    1254 : "CPXERR_NOT_UNBOUNDED",
    1255 : "CPXERR_SBASE_INCOMPAT",
    1256 : "CPXERR_SINGULAR",
    1257 : "CPXERR_PRIIND",
    1258 : "CPXERR_NO_LU_FACTOR",
    1260 : "CPXERR_NO_SENSIT",
    1261 : "CPXERR_NO_BASIC_SOLN",
    1262 : "CPXERR_NO_BASIS",
    1263 : "CPXERR_ABORT_STRONGBRANCH",
    1264 : "CPXERR_NO_NORMS",
    1265 : "CPXERR_NOT_DUAL_UNBOUNDED",
    1266 : "CPXERR_TILIM_STRONGBRANCH",
    1267 : "CPXERR_BAD_PIVOT",
    1268 : "CPXERR_TILIM_CONDITION_NO",
    1292 : "CPXERR_BAD_METHOD",
    1421 : "CPXERR_NO_FILENAME",
    1422 : "CPXERR_FAIL_OPEN_WRITE",
    1423 : "CPXERR_FAIL_OPEN_READ",
    1424 : "CPXERR_BAD_FILETYPE",
    1425 : "CPXERR_XMLPARSE",
    1431 : "CPXERR_TOO_MANY_ROWS",
    1432 : "CPXERR_TOO_MANY_COLS",
    1433 : "CPXERR_TOO_MANY_COEFFS",
    1434 : "CPXERR_BAD_NUMBER",
    1435 : "CPXERR_BAD_EXPO_RANGE",
    1436 : "CPXERR_NO_OBJ_SENSE",
    1437 : "CPXERR_QCP_SENSE_FILE",
    1438 : "CPXERR_BAD_LAZY_UCUT",
    1439 : "CPXERR_BAD_INDCONSTR",
    1441 : "CPXERR_NO_NAME_SECTION",
    1442 : "CPXERR_BAD_SOS_TYPE",
    1443 : "CPXERR_COL_ROW_REPEATS",
    1444 : "CPXERR_RIM_ROW_REPEATS",
    1445 : "CPXERR_ROW_REPEATS",
    1446 : "CPXERR_COL_REPEATS",
    1447 : "CPXERR_RIM_REPEATS",
    1448 : "CPXERR_ROW_UNKNOWN",
    1449 : "CPXERR_COL_UNKNOWN",
    1453 : "CPXERR_NO_ROW_SENSE",
    1454 : "CPXERR_EXTRA_FX_BOUND",
    1455 : "CPXERR_EXTRA_FR_BOUND",
    1456 : "CPXERR_EXTRA_BV_BOUND",
    1457 : "CPXERR_BAD_BOUND_TYPE",
    1458 : "CPXERR_UP_BOUND_REPEATS",
    1459 : "CPXERR_LO_BOUND_REPEATS",
    1460 : "CPXERR_NO_BOUND_TYPE",
    1461 : "CPXERR_NO_QMATRIX_SECTION",
    1462 : "CPXERR_BAD_SECTION_ENDATA",
    1463 : "CPXERR_INT_TOO_BIG_INPUT",
    1464 : "CPXERR_NAME_TOO_LONG",
    1465 : "CPXERR_LINE_TOO_LONG",
    1471 : "CPXERR_NO_ROWS_SECTION",
    1472 : "CPXERR_NO_COLUMNS_SECTION",
    1473 : "CPXERR_BAD_SECTION_BOUNDS",
    1474 : "CPXERR_RANGE_SECTION_ORDER",
    1475 : "CPXERR_BAD_SECTION_QMATRIX",
    1476 : "CPXERR_NO_OBJECTIVE",
    1477 : "CPXERR_ROW_REPEAT_PRINT",
    1478 : "CPXERR_COL_REPEAT_PRINT",
    1479 : "CPXERR_RIMNZ_REPEATS",
    1480 : "CPXERR_EXTRA_INTORG",
    1481 : "CPXERR_EXTRA_INTEND",
    1482 : "CPXERR_EXTRA_SOSORG",
    1483 : "CPXERR_EXTRA_SOSEND",
    1484 : "CPXERR_TOO_MANY_RIMS",
    1485 : "CPXERR_TOO_MANY_RIMNZ",
    1486 : "CPXERR_NO_ROW_NAME",
    1487 : "CPXERR_BAD_OBJ_SENSE",
    1550 : "CPXERR_BAS_FILE_SHORT",
    1551 : "CPXERR_BAD_INDICATOR",
    1552 : "CPXERR_NO_ENDATA",
    1553 : "CPXERR_FILE_ENTRIES",
    1554 : "CPXERR_SBASE_ILLEGAL",
    1555 : "CPXERR_BAS_FILE_SIZE",
    1556 : "CPXERR_NO_VECTOR_SOLN",
    1560 : "CPXERR_NOT_SAV_FILE",
    1561 : "CPXERR_SAV_FILE_DATA",
    1562 : "CPXERR_SAV_FILE_WRITE",
    1563 : "CPXERR_FILE_FORMAT",
    1602 : "CPXERR_ADJ_SIGNS",
    1603 : "CPXERR_RHS_IN_OBJ",
    1604 : "CPXERR_ADJ_SIGN_SENSE",
    1605 : "CPXERR_QUAD_IN_ROW",
    1606 : "CPXERR_ADJ_SIGN_QUAD",
    1607 : "CPXERR_NO_OPERATOR",
    1608 : "CPXERR_NO_OP_OR_SENSE",
    1609 : "CPXERR_NO_ID_FIRST",
    1610 : "CPXERR_NO_RHS_COEFF",
    1611 : "CPXERR_NO_NUMBER_FIRST",
    1612 : "CPXERR_NO_QUAD_EXP",
    1613 : "CPXERR_QUAD_EXP_NOT_2",
    1614 : "CPXERR_NO_QP_OPERATOR",
    1615 : "CPXERR_NO_NUMBER",
    1616 : "CPXERR_NO_ID",
    1617 : "CPXERR_BAD_ID",
    1618 : "CPXERR_BAD_EXPONENT",
    1619 : "CPXERR_Q_DIVISOR",
    1621 : "CPXERR_NO_BOUND_SENSE",
    1622 : "CPXERR_BAD_BOUND_SENSE",
    1623 : "CPXERR_NO_NUMBER_BOUND",
    1627 : "CPXERR_NO_SOS_SEPARATOR",
    1650 : "CPXERR_INVALID_NUMBER",
    1660 : "CPXERR_PRM_DATA",
    1661 : "CPXERR_PRM_HEADER",
    1719 : "CPXERR_NO_CONFLICT",
    1720 : "CPXERR_CONFLICT_UNSTABLE",
    1801 : "CPXERR_WORK_FILE_OPEN",
    1802 : "CPXERR_WORK_FILE_READ",
    1803 : "CPXERR_WORK_FILE_WRITE",
    1804 : "CPXERR_IN_INFOCALLBACK",
    1805 : "CPXERR_MIPSEARCH_WITH_CALLBACKS",
    1806 : "CPXERR_LP_NOT_IN_ENVIRONMENT",
    1807 : "CPXERR_PARAM_INCOMPATIBLE",
    32000 : "CPXERR_LICENSE_MIN",
    32201 : "CPXERR_ILOG_LICENSE",
    32301 : "CPXERR_NO_MIP_LIC",
    32302 : "CPXERR_NO_BARRIER_LIC",
    32305 : "CPXERR_NO_MIQP_LIC",
    32018 : "CPXERR_BADLDWID",
    32023 : "CPXERR_BADPRODUCT",
    32024 : "CPXERR_ALGNOTLICENSED",
    32999 : "CPXERR_LICENSE_MAX",
    1701 : "CPXERR_IIS_NO_INFO",
    1702 : "CPXERR_IIS_NO_SOLN",
    1703 : "CPXERR_IIS_FEAS",
    1704 : "CPXERR_IIS_NOT_INFEAS",
    1705 : "CPXERR_IIS_OPT_INFEAS",
    1706 : "CPXERR_IIS_DEFAULT",
    1707 : "CPXERR_IIS_NO_BASIC",
    1709 : "CPXERR_IIS_NO_LOAD",
    1710 : "CPXERR_IIS_SUB_OBJ_LIM",
    1711 : "CPXERR_IIS_SUB_IT_LIM",
    1712 : "CPXERR_IIS_SUB_TIME_LIM",
    1713 : "CPXERR_IIS_NUM_BEST",
    1714 : "CPXERR_IIS_SUB_ABORT",
    3003 : "CPXERR_NOT_MIP",
    3006 : "CPXERR_BAD_PRIORITY",
    3007 : "CPXERR_ORDER_BAD_DIRECTION",
    3009 : "CPXERR_ARRAY_BAD_SOS_TYPE",
    3010 : "CPXERR_UNIQUE_WEIGHTS",
    3012 : "CPXERR_BAD_DIRECTION",
    3015 : "CPXERR_NO_SOS",
    3016 : "CPXERR_NO_ORDER",
    3018 : "CPXERR_INT_TOO_BIG",
    3019 : "CPXERR_SUBPROB_SOLVE",
    3020 : "CPXERR_NO_MIPSTART",
    3021 : "CPXERR_BAD_CTYPE",
    3023 : "CPXERR_NO_INT_X",
    3024 : "CPXERR_NO_SOLNPOOL",
    3301 : "CPXERR_MISS_SOS_TYPE",
    3412 : "CPXERR_NO_TREE",
    3413 : "CPXERR_TREE_MEMORY_LIMIT",
    3414 : "CPXERR_FILTER_VARIABLE_TYPE",
    3504 : "CPXERR_NODE_ON_DISK",
    3601 : "CPXERR_PTHREAD_MUTEX_INIT",
    3603 : "CPXERR_PTHREAD_CREATE",
    1212 : "CPXERR_UNSUPPORTED_CONSTRAINT_TYPE",
    1213 : "CPXERR_ILL_DEFINED_PWL",
    1530 : "CPXERR_NET_DATA",
    1531 : "CPXERR_NOT_MIN_COST_FLOW",
    1532 : "CPXERR_BAD_ROW_ID",
    1537 : "CPXERR_BAD_CHAR",
    1538 : "CPXERR_NET_FILE_SHORT",
    5002 : "CPXERR_Q_NOT_POS_DEF",
    5004 : "CPXERR_NOT_QP",
    5011 : "CPXERR_Q_DUP_ENTRY",
    5012 : "CPXERR_Q_NOT_SYMMETRIC",
    5014 : "CPXERR_Q_NOT_INDEF",
    6002 : "CPXERR_QCP_SENSE"
    }

cdef dict parameters
parameters = {
    "CPX_PARAMTYPE_NONE" : 0,
    "CPX_PARAMTYPE_INT" : 1,
    "CPX_PARAMTYPE_DOUBLE" : 2,
    "CPX_PARAMTYPE_STRING" : 3,
    "CPX_PARAM_ADVIND" : 1001,
    "CPX_PARAM_AGGFILL" : 1002,
    "CPX_PARAM_AGGIND" : 1003,
    "CPX_PARAM_BASINTERVAL" : 1004,
    "CPX_PARAM_CFILEMUL" : 1005,
    "CPX_PARAM_CLOCKTYPE" : 1006,
    "CPX_PARAM_CRAIND" : 1007,
    "CPX_PARAM_DEPIND" : 1008,
    "CPX_PARAM_DPRIIND" : 1009,
    "CPX_PARAM_PRICELIM" : 1010,
    "CPX_PARAM_EPMRK" : 1013,
    "CPX_PARAM_EPOPT" : 1014,
    "CPX_PARAM_EPPER" : 1015,
    "CPX_PARAM_EPRHS" : 1016,
    "CPX_PARAM_FASTMIP" : 1017,
    "CPX_PARAM_SIMDISPLAY" : 1019,
    "CPX_PARAM_ITLIM" : 1020,
    "CPX_PARAM_ROWREADLIM" : 1021,
    "CPX_PARAM_NETFIND" : 1022,
    "CPX_PARAM_COLREADLIM" : 1023,
    "CPX_PARAM_NZREADLIM" : 1024,
    "CPX_PARAM_OBJLLIM" : 1025,
    "CPX_PARAM_OBJULIM" : 1026,
    "CPX_PARAM_PERIND" : 1027,
    "CPX_PARAM_PERLIM" : 1028,
    "CPX_PARAM_PPRIIND" : 1029,
    "CPX_PARAM_PREIND" : 1030,
    "CPX_PARAM_REINV" : 1031,
    "CPX_PARAM_REVERSEIND" : 1032,
    "CPX_PARAM_RFILEMUL" : 1033,
    "CPX_PARAM_SCAIND" : 1034,
    "CPX_PARAM_SCRIND" : 1035,
    "CPX_PARAM_SINGLIM" : 1037,
    "CPX_PARAM_SINGTOL" : 1038,
    "CPX_PARAM_TILIM" : 1039,
    "CPX_PARAM_XXXIND" : 1041,
    "CPX_PARAM_PREDUAL" : 1044,
    "CPX_PARAM_EPOPT_H" : 1049,
    "CPX_PARAM_EPRHS_H" : 1050,
    "CPX_PARAM_PREPASS" : 1052,
    "CPX_PARAM_DATACHECK" : 1056,
    "CPX_PARAM_REDUCE" : 1057,
    "CPX_PARAM_PRELINEAR" : 1058,
    "CPX_PARAM_LPMETHOD" : 1062,
    "CPX_PARAM_QPMETHOD" : 1063,
    "CPX_PARAM_WORKDIR" : 1064,
    "CPX_PARAM_WORKMEM" : 1065,
    "CPX_PARAM_THREADS" : 1067,
    "CPX_PARAM_CONFLICTDISPLAY" : 1074,
    "CPX_PARAM_SIFTDISPLAY" : 1076,
    "CPX_PARAM_SIFTALG" : 1077,
    "CPX_PARAM_SIFTITLIM" : 1078,
    "CPX_PARAM_MPSLONGNUM" : 1081,
    "CPX_PARAM_MEMORYEMPHASIS" : 1082,
    "CPX_PARAM_NUMERICALEMPHASIS" : 1083,
    "CPX_PARAM_FEASOPTMODE" : 1084,
    "CPX_PARAM_PARALLELMODE" : 1109,
    "CPX_PARAM_TUNINGMEASURE" : 1110,
    "CPX_PARAM_TUNINGREPEAT" : 1111,
    "CPX_PARAM_TUNINGTILIM" : 1112,
    "CPX_PARAM_TUNINGDISPLAY" : 1113,
    "CPX_PARAM_WRITELEVEL" : 1114,
    "CPX_PARAM_ALL_MIN" : 1000,
    "CPX_PARAM_ALL_MAX" : 6000,
    "CPX_PARAM_BARDSTART" : 3001,
    "CPX_PARAM_BAREPCOMP" : 3002,
    "CPX_PARAM_BARGROWTH" : 3003,
    "CPX_PARAM_BAROBJRNG" : 3004,
    "CPX_PARAM_BARPSTART" : 3005,
    "CPX_PARAM_BARALG" : 3007,
    "CPX_PARAM_BARCOLNZ" : 3009,
    "CPX_PARAM_BARDISPLAY" : 3010,
    "CPX_PARAM_BARITLIM" : 3012,
    "CPX_PARAM_BARMAXCOR" : 3013,
    "CPX_PARAM_BARORDER" : 3014,
    "CPX_PARAM_BARSTARTALG" : 3017,
    "CPX_PARAM_BARCROSSALG" : 3018,
    "CPX_PARAM_BARQCPEPCOMP" : 3020,
    "CPX_PARAM_BRDIR" : 2001,
    "CPX_PARAM_BTTOL" : 2002,
    "CPX_PARAM_CLIQUES" : 2003,
    "CPX_PARAM_COEREDIND" : 2004,
    "CPX_PARAM_COVERS" : 2005,
    "CPX_PARAM_CUTLO" : 2006,
    "CPX_PARAM_CUTUP" : 2007,
    "CPX_PARAM_EPAGAP" : 2008,
    "CPX_PARAM_EPGAP" : 2009,
    "CPX_PARAM_EPINT" : 2010,
    "CPX_PARAM_MIPDISPLAY" : 2012,
    "CPX_PARAM_MIPINTERVAL" : 2013,
    "CPX_PARAM_INTSOLLIM" : 2015,
    "CPX_PARAM_NODEFILEIND" : 2016,
    "CPX_PARAM_NODELIM" : 2017,
    "CPX_PARAM_NODESEL" : 2018,
    "CPX_PARAM_OBJDIF" : 2019,
    "CPX_PARAM_MIPORDIND" : 2020,
    "CPX_PARAM_RELOBJDIF" : 2022,
    "CPX_PARAM_STARTALG" : 2025,
    "CPX_PARAM_SUBALG" : 2026,
    "CPX_PARAM_TRELIM" : 2027,
    "CPX_PARAM_VARSEL" : 2028,
    "CPX_PARAM_BNDSTRENIND" : 2029,
    "CPX_PARAM_HEURFREQ" : 2031,
    "CPX_PARAM_MIPORDTYPE" : 2032,
    "CPX_PARAM_CUTSFACTOR" : 2033,
    "CPX_PARAM_RELAXPREIND" : 2034,
    "CPX_PARAM_PRESLVND" : 2037,
    "CPX_PARAM_BBINTERVAL" : 2039,
    "CPX_PARAM_FLOWCOVERS" : 2040,
    "CPX_PARAM_IMPLBD" : 2041,
    "CPX_PARAM_PROBE" : 2042,
    "CPX_PARAM_GUBCOVERS" : 2044,
    "CPX_PARAM_STRONGCANDLIM" : 2045,
    "CPX_PARAM_STRONGITLIM" : 2046,
    "CPX_PARAM_FRACCAND" : 2048,
    "CPX_PARAM_FRACCUTS" : 2049,
    "CPX_PARAM_FRACPASS" : 2050,
    "CPX_PARAM_FLOWPATHS" : 2051,
    "CPX_PARAM_MIRCUTS" : 2052,
    "CPX_PARAM_DISJCUTS" : 2053,
    "CPX_PARAM_AGGCUTLIM" : 2054,
    "CPX_PARAM_MIPCBREDLP" : 2055    ,
    "CPX_PARAM_CUTPASS" : 2056,
    "CPX_PARAM_MIPEMPHASIS" : 2058,
    "CPX_PARAM_SYMMETRY" : 2059,
    "CPX_PARAM_DIVETYPE" : 2060,
    "CPX_PARAM_RINSHEUR" : 2061,
    "CPX_PARAM_SUBMIPNODELIM" : 2062,
    "CPX_PARAM_LBHEUR" : 2063,
    "CPX_PARAM_REPEATPRESOLVE" : 2064,
    "CPX_PARAM_PROBETIME" : 2065,
    "CPX_PARAM_POLISHTIME" : 2066,
    "CPX_PARAM_REPAIRTRIES" : 2067,
    "CPX_PARAM_EPLIN" : 2068,
    "CPX_PARAM_EPRELAX" : 2073,
    "CPX_PARAM_FPHEUR" : 2098,
    "CPX_PARAM_EACHCUTLIM" : 2102,
    "CPX_PARAM_SOLNPOOLCAPACITY" : 2103 ,
    "CPX_PARAM_SOLNPOOLREPLACE" : 2104 ,
    "CPX_PARAM_SOLNPOOLGAP" : 2105 ,
    "CPX_PARAM_SOLNPOOLAGAP" : 2106 ,
    "CPX_PARAM_SOLNPOOLINTENSITY" : 2107,
    "CPX_PARAM_POPULATELIM" : 2108,
    "CPX_PARAM_MIPSEARCH" : 2109,
    "CPX_PARAM_MIQCPSTRAT" : 2110,
    "CPX_PARAM_ZEROHALFCUTS" : 2111,
    "CPX_PARAM_POLISHAFTEREPAGAP" : 2126,
    "CPX_PARAM_POLISHAFTEREPGAP" : 2127,
    "CPX_PARAM_POLISHAFTERNODE" : 2128,
    "CPX_PARAM_POLISHAFTERINTSOL" : 2129,
    "CPX_PARAM_POLISHAFTERTIME" : 2130,
    "CPX_PARAM_MCFCUTS" : 2134,
    "CPX_PARAM_NETITLIM" : 5001,
    "CPX_PARAM_NETEPOPT" : 5002,
    "CPX_PARAM_NETEPRHS" : 5003,
    "CPX_PARAM_NETPPRIIND" : 5004,
    "CPX_PARAM_NETDISPLAY" : 5005,
    "CPX_PARAM_QPNZREADLIM" : 4001,
    "CPX_PARAM_QPMAKEPSDIND" : 4010,
    }
