
"""
COIN Backend

AUTHORS:

- Nathann Cohen (2010-10): initial implementation

- John Perry (2012-03): major modifications to the interface in order to update
  the CBC package to version 2.7.5
"""

##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from __future__ import print_function

include "cysignals/memory.pxi"
include "cysignals/signals.pxi"

from sage.numerical.mip import MIPSolverException
from copy import copy

cdef class CoinBackend(GenericBackend):

    """
    MIP Backend that uses the COIN solver (CBC).

    TESTS:

    General backend testsuite::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")                       # optional - cbc
            sage: TestSuite(p).run(skip="_test_pickling")               # optional - cbc
    """

    def __cinit__(self, maximization = True):
        """
        Cython constructor

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")                  # optional - cbc

        """

        # Coin devs seem to favor OsiClpSolverInterface
        self.si = new OsiClpSolverInterface()
        self.model =  new CbcModel(self.si[0])
        self.prob_name = None
        self.row_names = []
        self.col_names = []
        self.set_verbosity(0)

        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

        self.obj_constant_term = 0.0

    def __dealloc__(self):
        r"""
        Destructor function
        """
        del self.si
        del self.model

    cpdef int add_variable(self, lower_bound=0.0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0.0, name=None) except -1:
        r"""
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

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")                  # optional - cbc
            sage: p.ncols()                                         # optional - cbc
            0
            sage: p.add_variable()                                  # optional - cbc
            0
            sage: p.ncols()                                         # optional - cbc
            1
            sage: p.add_variable(binary=True)                       # optional - cbc
            1
            sage: p.add_variable(lower_bound=-2.0, integer=True)    # optional - cbc
            2
            sage: p.add_variable(continuous=True, integer=True)     # optional - cbc
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.add_variable(name='x',obj=1.0)                  # optional - cbc
            3
            sage: p.col_name(3)                                     # optional - cbc
            'x'
            sage: p.objective_coefficient(3)                        # optional - cbc
            1.0
        """

        # for some reason, Cython is not accepting the line below, which appeare
        #cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        cdef int vtype = int(binary) + int(continuous) + int(integer)
        if  vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        self.si.addCol(0, NULL, NULL, 0, self.si.getInfinity(), 0)

        cdef int n
        n = self.si.getNumCols() - 1

        if lower_bound != 0.0:
            self.variable_lower_bound(n, lower_bound)
        if upper_bound is not None:
            self.variable_upper_bound(n, upper_bound)

        if binary:
            self.set_variable_type(n,0)
        elif integer:
            self.set_variable_type(n,1)

        if name:
            self.col_names.append(name)
        else:
            self.col_names.append("")

        if obj:
            self.si.setObjCoeff(n, obj)

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

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")                         # optional - cbc
            sage: p.ncols()                                                # optional - cbc
            0
            sage: p.add_variables(5)                                       # optional - cbc
            4
            sage: p.ncols()                                                # optional - cbc
            5
            sage: p.add_variables(2, lower_bound=-2.0, integer=True, obj=42.0, names=['a','b']) # optional - cbc
            6

        TESTS:

        Check that arguments are used::

            sage: p.col_bounds(5) # tol 1e-8, optional - cbc
            (-2.0, None)
            sage: p.is_variable_integer(5)   # optional - cbc
            True
            sage: p.col_name(5)              # optional - cbc
            'a'
            sage: p.objective_coefficient(5) # tol 1e-8, optional - cbc
            42.0
        """
        #cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        cdef int vtype = int(binary) + int(continuous) + int(integer)
        if  vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        cdef int n
        n = self.si.getNumCols()

        cdef int i

        for 0<= i < number:

            self.si.addCol(0, NULL, NULL, 0, self.si.getInfinity(), 0)

            if lower_bound != 0.0:
                self.variable_lower_bound(n + i, lower_bound)
            if upper_bound is not None:
                self.variable_upper_bound(n + i, upper_bound)

            if binary:
                self.set_variable_type(n + i,0)
            elif integer:
                self.set_variable_type(n + i,1)

            if obj:
                self.si.setObjCoeff(n + i, obj)

        if names is not None:
            for name in names:
                self.col_names.append(name)
        else:
            self.col_names.extend(['' for i in range(number)])

        return n + number -1

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

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")   # optional - cbc
            sage: p.ncols()                                        # optional - cbc
            0
            sage: p.add_variable()                                  # optional - cbc
            0
            sage: p.set_variable_type(0,1)                          # optional - cbc
            sage: p.is_variable_integer(0)                          # optional - cbc
            True
        """

        if vtype == 1:
            self.si.setInteger(variable)
        elif vtype == 0:
            self.si.setColLower(variable, 0)
            self.si.setInteger(variable)
            self.si.setColUpper(variable, 1)
        else:
            self.si.setContinuous(variable)

    cpdef set_sense(self, int sense):
        r"""
        Sets the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.is_maximization()                              # optional - cbc
            True
            sage: p.set_sense(-1)                              # optional - cbc
            sage: p.is_maximization()                              # optional - cbc
            False
        """
        self.si.setObjSense(-sense)

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient or ``None`` for
          reading (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")       # optional -- cbc
            sage: p.add_variable()                      # optional -- cbc
            0
            sage: p.objective_coefficient(0)            # optional -- cbc
            0.0
            sage: p.objective_coefficient(0,2)          # optional -- cbc
            sage: p.objective_coefficient(0)            # optional -- cbc
            2.0
        """
        if coeff is not None:
            self.si.setObjCoeff(variable, coeff)
        else:
            return self.si.getObjCoefficients()[variable]

    cpdef set_objective(self, list coeff, d = 0.0):
        r"""
        Sets the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        - ``d`` (double) -- the constant term in the linear function (set to `0` by default)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")    # optional - cbc
            sage: p.add_variables(5)                                 # optional - cbc
            4
            sage: p.set_objective([1, 1, 2, 1, 3])                   # optional - cbc
            sage: map(lambda x :p.objective_coefficient(x), range(5))  # optional - cbc
            [1.0, 1.0, 2.0, 1.0, 3.0]

        Constants in the objective function are respected::

            sage: p = MixedIntegerLinearProgram(solver='Coin')  # optional - cbc
            sage: v = p.new_variable(nonnegative=True)          # optional - cbc
            sage: x,y = v[0], v[1]                              # optional - cbc
            sage: p.add_constraint(2*x + 3*y, max = 6)          # optional - cbc
            sage: p.add_constraint(3*x + 2*y, max = 6)          # optional - cbc
            sage: p.set_objective(x + y + 7)                    # optional - cbc
            sage: p.set_integer(x); p.set_integer(y)            # optional - cbc
            sage: p.solve()                                     # optional - cbc
            9.0
        """

        cdef int i

        for i,v in enumerate(coeff):
            self.si.setObjCoeff(i, v)

        self.obj_constant_term = d

    cpdef set_verbosity(self, int level):
        r"""
        Sets the log (verbosity) level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")   # optional - cbc
            sage: p.set_verbosity(2)                                # optional - cbc

        """
        self.model.setLogLevel(level)

    cpdef remove_constraint(self, int i):
        r"""
        Remove a constraint from self.

        INPUT:

        - ``i`` -- index of the constraint to remove

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver='Coin') # optional - cbc
            sage: v = p.new_variable(nonnegative=True)         # optional - cbc
            sage: x,y = v[0], v[1]                             # optional - cbc
            sage: p.add_constraint(2*x + 3*y, max = 6)         # optional - cbc
            sage: p.add_constraint(3*x + 2*y, max = 6)         # optional - cbc
            sage: p.set_objective(x + y + 7)                   # optional - cbc
            sage: p.set_integer(x); p.set_integer(y)           # optional - cbc
            sage: p.solve()                                    # optional - cbc
            9.0
            sage: p.remove_constraint(0)                       # optional - cbc
            sage: p.solve()                                    # optional - cbc
            10.0
            sage: p.get_values([x,y])                          # optional - cbc
            [0.0, 3.0]

        TESTS:

        Removing fancy constraints does not make Sage crash::

            sage: MixedIntegerLinearProgram(solver='Coin').remove_constraint(-2)  # optional - cbc
            Traceback (most recent call last):
            ...
            ValueError: The constraint's index i must satisfy 0 <= i < number_of_constraints
        """
        cdef int rows [1]

        if i < 0 or i >= self.si.getNumRows():
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")
        rows[0] = i
        self.si.deleteRows(1,rows)

    cpdef remove_constraints(self, constraints):
        r"""
        Remove several constraints.

        INPUT:

        - ``constraints`` -- an interable containing the indices of the rows to remove

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver='Coin') # optional - cbc
            sage: v = p.new_variable(nonnegative=True)         # optional - cbc
            sage: x,y = v[0], v[1]                             # optional - cbc
            sage: p.add_constraint(2*x + 3*y, max = 6)         # optional - cbc
            sage: p.add_constraint(3*x + 2*y, max = 6)         # optional - cbc
            sage: p.set_objective(x + y + 7)                   # optional - cbc
            sage: p.set_integer(x); p.set_integer(y)           # optional - cbc
            sage: p.solve()                                    # optional - cbc
            9.0
            sage: p.get_values(x)                              # optional - cbc
            2.0...
            sage: p.get_values(y)                              # optional - cbc
            0.0...
            sage: p.remove_constraints([0])                    # optional - cbc
            sage: p.solve()                                    # optional - cbc
            10.0
            sage: p.get_values([x,y])                          # optional - cbc
            [0.0, 3.0]

        TESTS:

        Removing fancy constraints do not make Sage crash::

            sage: MixedIntegerLinearProgram(solver='Coin').remove_constraints([0, -2])  # optional - cbc
            Traceback (most recent call last):
            ...
            ValueError: The constraint's index i must satisfy 0 <= i < number_of_constraints
        """
        cdef int i, c
        cdef int m = len(constraints)
        cdef int * rows = <int *>check_malloc(m * sizeof(int *))
        cdef int nrows = self.si.getNumRows()

        for i in xrange(m):

            c = constraints[i]
            if c < 0 or c >= nrows:
                sig_free(rows)
                raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")

            rows[i] = c

        self.si.deleteRows(m,rows)
        sig_free(rows)

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

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")                              # optional - cbc
            sage: p.add_variables(5)                                           # optional - cbc
            4
            sage: p.add_linear_constraint( zip(range(5), range(5)), 2.0, 2.0)  # optional - cbc
            sage: p.row(0)                                                     # optional - cbc
            ([0, 1, 2, 3, 4], [0.0, 1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                                              # optional - cbc
            (2.0, 2.0)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1.0, 1.0, name='foo') # optional - cbc
            sage: p.row_name(1)                                                           # optional - cbc
            'foo'
        """
        if lower_bound is None and upper_bound is None:
            raise ValueError("At least one of 'upper_bound' or 'lower_bound' must be set.")

        cdef int i
        cdef double c
        cdef CoinPackedVector* row
        row = new_CoinPackedVector();


        for i,c in coefficients:
            row.insert(i, c)

        self.si.addRow (row[0],
                        lower_bound if lower_bound is not None else -self.si.getInfinity(),
                        upper_bound if upper_bound is not None else +self.si.getInfinity())
        if name is not None:
            self.row_names.append(name)
        else:
            self.row_names.append("")

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

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.add_variables(5)                               # optional - cbc
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)       # optional - cbc
            sage: p.row(0)                                     # optional - cbc
            ([0, 1, 2, 3, 4], [0.0, 1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                              # optional - cbc
            (2.0, 2.0)
        """

        cdef list indices = []
        cdef list values = []
        cdef int * c_indices
        cdef int i
        cdef double * c_values
        cdef CoinPackedMatrix * M = <CoinPackedMatrix *> self.si.getMatrixByRow()
        cdef CoinShallowPackedVector V = <CoinShallowPackedVector> M.getVector(index)
        cdef int n = V.getNumElements()

        c_indices = <int*> V.getIndices()
        c_values = <double*> V.getElements()

        for 0<= i < n:
            indices.append(c_indices[i])
            values.append(c_values[i])

        return (indices, values)

    cpdef row_bounds(self, int i):
        r"""
        Returns the bounds of a specific constraint.

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the constraint is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.add_variables(5)                               # optional - cbc
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)       # optional - cbc
            sage: p.row(0)                                     # optional - cbc
            ([0, 1, 2, 3, 4], [0.0, 1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                              # optional - cbc
            (2.0, 2.0)
        """

        cdef double * ub
        cdef double * lb

        ub = <double*> self.si.getRowUpper()
        lb = <double*> self.si.getRowLower()

        return (lb[i] if lb[i] != - self.si.getInfinity() else None,
                ub[i] if ub[i] != + self.si.getInfinity() else None)

    cpdef col_bounds(self, int i):
        r"""
        Returns the bounds of a specific variable.

        INPUT:

        - ``index`` (integer) -- the variable's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the variable is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.add_variable()                                 # optional - cbc
            0
            sage: p.col_bounds(0)                              # optional - cbc
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                         # optional - cbc
            sage: p.col_bounds(0)                              # optional - cbc
            (0.0, 5.0)
        """

        cdef double * ub
        cdef double * lb

        ub = <double*> self.si.getColUpper()
        lb = <double*> self.si.getColLower()

        return (lb[i] if lb[i] != - self.si.getInfinity() else None,
                ub[i] if ub[i] != + self.si.getInfinity() else None)

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

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.ncols()                                       # optional - cbc
            0
            sage: p.nrows()                                       # optional - cbc
            0
            sage: p.add_linear_constraints(5, 0, None)                      # optional - cbc
            sage: p.add_col(range(5), range(5))                    # optional - cbc
            sage: p.nrows()                                       # optional - cbc
            5
        """

        cdef int n = len(indices)
        cdef int * c_indices = <int*>check_malloc(n*sizeof(int))
        cdef double * c_values  = <double*>check_malloc(n*sizeof(double))
        cdef int i

        for 0<= i< n:
            c_indices[i] = indices[i]
            c_values[i] = coeffs[i]

        self.si.addCol (n, c_indices, c_values, 0, self.si.getInfinity(), 0)

        self.col_names.append("")


    cpdef int solve(self) except -1:
        r"""
        Solves the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")    # optional - cbc
            sage: p.add_linear_constraints(5, 0, None)       # optional - cbc
            sage: p.add_col(range(5), [1,2,3,4,5])  # optional - cbc
            sage: p.solve()                         # optional - cbc
            0

        TESTS::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")    # optional - cbc
            sage: p.add_variable()                  # optional - cbc
            0
            sage: p.add_linear_constraint([(0, 1)], None, 4) # optional - cbc
            sage: p.add_linear_constraint([(0, 1)], 6, None) # optional - cbc
            sage: p.objective_coefficient(0,1)        # optional - cbc
            sage: p.solve()                         # optional - cbc
            Traceback (most recent call last):
            ...
            MIPSolverException: ...
        """

        # set up the model
        cdef OsiSolverInterface * si = self.si

        cdef CbcModel * model
        cdef int old_logLevel = self.model.logLevel()

        model = new CbcModel(si[0])
        del self.model
        self.model = model
        
        #we immediately commit to the new model so that the user has access
        #to it even when something goes wrong.

        model.setLogLevel(old_logLevel)

        # multithreading
        import multiprocessing
        model.setNumberThreads(multiprocessing.cpu_count())

        model.branchAndBound()

        if model.solver().isAbandoned():
            raise MIPSolverException("CBC : The solver has abandoned!")

        elif model.solver().isProvenPrimalInfeasible() or model.solver().isProvenDualInfeasible():
            raise MIPSolverException("CBC : The problem or its dual has been proven infeasible!")

        elif (model.solver().isPrimalObjectiveLimitReached() or model.solver().isDualObjectiveLimitReached()):
            raise MIPSolverException("CBC : The objective limit has been reached for the problem or its dual!")

        elif model.solver().isIterationLimitReached():
            raise MIPSolverException("CBC : The iteration limit has been reached!")

        elif not model.solver().isProvenOptimal():
            raise MIPSolverException("CBC : Unknown error")

        return 0

    cpdef get_objective_value(self):
        r"""
        Returns the value of the objective function.

        .. NOTE::

           Has no meaning unless ``solve`` or ``set_basis_status`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.add_variables(2)                               # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3)          # optional - cbc
            sage: p.set_objective([2, 5])                          # optional - cbc
            sage: p.solve()                                        # optional - cbc
            0
            sage: p.get_objective_value()                          # optional - cbc
            7.5
            sage: p.get_variable_value(0)                          # optional - cbc
            0.0
            sage: p.get_variable_value(1)                          # optional - cbc
            1.5
        """
        return self.model.solver().getObjValue() + self.obj_constant_term

    cpdef get_variable_value(self, int variable):
        r"""
        Returns the value of a variable given by the solver.

        .. NOTE::

           Has no meaning unless ``solve`` or ``set_basis_status`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin") # optional - cbc
            sage: p.add_variables(2)                              # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3)          # optional - cbc
            sage: p.set_objective([2, 5])                         # optional - cbc
            sage: p.solve()                                       # optional - cbc
            0
            sage: p.get_objective_value()                         # optional - cbc
            7.5
            sage: p.get_variable_value(0)                         # optional - cbc
            0.0
            sage: p.get_variable_value(1)                         # optional - cbc
            1.5
            sage: p = MixedIntegerLinearProgram("Coin")  # optional - cbc
            sage: x = p.new_variable(nonnegative=True)   # optional - cbc
            sage: p.set_min(x[0], 0.0)            # optional - cbc
            sage: p.get_values(x)                 # optional - cbc
            {0: 0.0}
        """

        cdef double * solution
        cdef double v
        solution = <double*> self.model.solver().getColSolution()
        if solution == NULL:
            v = 0.0
        else:
            v = solution[variable]
        if self.is_variable_continuous(variable):
            return v
        else:
            return round(v)

    cpdef int ncols(self):
        r"""
        Returns the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.ncols()                                       # optional - cbc
            0
            sage: p.add_variables(2)                               # optional - cbc
            1
            sage: p.ncols()                                       # optional - cbc
            2
        """

        return self.si.getNumCols()

    cpdef int nrows(self):
        r"""
        Returns the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin") # optional - cbc
            sage: p.nrows()                                      # optional - cbc
            0
            sage: p.add_linear_constraints(2, 2, None)                     # optional - cbc
            sage: p.nrows()                                      # optional - cbc
            2
        """
        return self.si.getNumRows()


    cpdef bint is_variable_binary(self, int index):
        r"""
        Tests whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.ncols()                                       # optional - cbc
            0
            sage: p.add_variable()                                 # optional - cbc
            0
            sage: p.set_variable_type(0,0)                         # optional - cbc
            sage: p.is_variable_binary(0)                          # optional - cbc
            True

        """

        return (0 == self.si.isContinuous(index) and
                self.variable_lower_bound(index) == 0 and
                self.variable_upper_bound(index) == 1)

    cpdef bint is_variable_integer(self, int index):
        r"""
        Tests whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.ncols()                                       # optional - cbc
            0
            sage: p.add_variable()                                 # optional - cbc
            0
            sage: p.set_variable_type(0,1)                         # optional - cbc
            sage: p.is_variable_integer(0)                         # optional - cbc
            True
        """
        return (0 == self.si.isContinuous(index) and
                (self.variable_lower_bound(index) != 0 or
                self.variable_upper_bound(index) != 1))

    cpdef bint is_variable_continuous(self, int index):
        r"""
        Tests whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.ncols()                                       # optional - cbc
            0
            sage: p.add_variable()                                 # optional - cbc
            0
            sage: p.is_variable_continuous(0)                      # optional - cbc
            True
            sage: p.set_variable_type(0,1)                         # optional - cbc
            sage: p.is_variable_continuous(0)                      # optional - cbc
            False

        """
        return 1 == self.si.isContinuous(index)


    cpdef bint is_maximization(self):
        r"""
        Tests whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin") # optional - cbc
            sage: p.is_maximization()                             # optional - cbc
            True
            sage: p.set_sense(-1)                             # optional - cbc
            sage: p.is_maximization()                             # optional - cbc
            False
        """

        return self.si.getObjSense() == -1

    cpdef variable_upper_bound(self, int index, value = False):
        r"""
        Returns or defines the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.add_variable()                                 # optional - cbc
            0
            sage: p.col_bounds(0)                              # optional - cbc
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                         # optional - cbc
            sage: p.col_bounds(0)                              # optional - cbc
            (0.0, 5.0)

        TESTS:

        :trac:`14581`::

            sage: P = MixedIntegerLinearProgram(solver="Coin") # optional - cbc
            sage: v = P.new_variable(nonnegative=True)         # optional - cbc
            sage: x = v["x"]                                   # optional - cbc
            sage: P.set_max(x, 0)                              # optional - cbc
            sage: P.get_max(x)                                 # optional - cbc
            0.0

        """
        cdef double * ub

        if value is False:
            ub = <double*> self.si.getColUpper()
            return ub[index] if ub[index] != + self.si.getInfinity() else None
        else:
            self.si.setColUpper(index, value if value is not None else +self.si.getInfinity())

    cpdef variable_lower_bound(self, int index, value = False):
        r"""
        Returns or defines the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.add_variable()                                 # optional - cbc
            0
            sage: p.col_bounds(0)                              # optional - cbc
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)                         # optional - cbc
            sage: p.col_bounds(0)                              # optional - cbc
            (5.0, None)

        TESTS:

        :trac:`14581`::

            sage: P = MixedIntegerLinearProgram(solver="Coin") # optional - cbc
            sage: v = P.new_variable(nonnegative=True)         # optional - cbc
            sage: x = v["x"]                                   # optional - cbc
            sage: P.set_min(x, 5)                              # optional - cbc
            sage: P.set_min(x, 0)                              # optional - cbc
            sage: P.get_min(x)                                 # optional - cbc
            0.0
        """
        cdef double * lb

        if value is False:
            lb = <double*> self.si.getColLower()
            return lb[index] if lb[index] != - self.si.getInfinity() else None
        else:
            self.si.setColLower(index, value if value is not None else -self.si.getInfinity())

    cpdef write_mps(self, char * filename, int modern):
        r"""
        Writes the problem to a .mps file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.add_variables(2)                               # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3)          # optional - cbc
            sage: p.set_objective([2, 5])                          # optional - cbc
            sage: p.write_mps(os.path.join(SAGE_TMP, "lp_problem.mps"), 0)       # optional - cbc
        """

        cdef char * mps = "mps"
        self.si.writeMps(filename, mps, -1 if self.is_maximization() else 1)

    cpdef write_lp(self, char * filename):
        r"""
        Writes the problem to a .lp file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")  # optional - cbc
            sage: p.add_variables(2)                               # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3)          # optional - cbc
            sage: p.set_objective([2, 5])                          # optional - cbc
            sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))       # optional - cbc
        """

        cdef char * lp = "lp"
        self.si.writeLp(filename, lp, 0.00001, 10, 5, -1 if self.is_maximization() else 1, 1)

    cpdef problem_name(self, char * name = NULL):
        r"""
        Returns or defines the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")   # optional - cbc
            sage: p.problem_name("There once was a french fry") # optional - cbc
            sage: print(p.problem_name())                       # optional - cbc
            There once was a french fry
        """
        if name == NULL:
            if self.prob_name is not None:
                return self.prob_name
            else:
                return ""
        else:
            self.prob_name = str(name)


    cpdef row_name(self, int index):
        r"""
        Returns the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")                                     # optional - cbc
            sage: p.add_linear_constraints(1, 2, None, names=['Empty constraint 1'])  # optional - cbc
            sage: print(p.row_name(0))                                                # optional - cbc
            Empty constraint 1
        """
        if self.row_names is not None:
            return self.row_names[index]
        else:
            return ""

    cpdef col_name(self, int index):
        r"""
        Returns the ``index`` th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")          # optional - cbc
            sage: p.add_variable(name='I am a variable')   # optional - cbc
            0
            sage: print(p.col_name(0))                     # optional - cbc
            I am a variable
        """
        if self.col_names is not None:
            return self.col_names[index]
        else:
            return ""

    cpdef __copy__(self):
        """
        Returns a copy of self.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "Coin")        # optional - cbc
            sage: b = p.new_variable(nonnegative=True)                  # optional - cbc
            sage: p.add_constraint(b[1] + b[2] <= 6)           # optional - cbc
            sage: p.set_objective(b[1] + b[2])                 # optional - cbc
            sage: copy(p).solve()                              # optional - cbc
            6.0
        """
        # create new backend
        cdef CoinBackend p = type(self)(maximization = (1 if self.is_maximization() else -1))

        # replace solver with copy of self's solver
        del p.si
        p.si = self.si.clone(1)
        p.row_names = copy(self.row_names)
        p.col_names = copy(self.col_names)
        p.obj_constant_term = self.obj_constant_term
        # Maybe I should copy this, not sure -- seems complicated, though
        p.prob_name = self.prob_name

        return p

    cpdef get_basis_status(self):
        """
        Retrieve status information for column and row variables.

        This method returns status as integer codes:

          * 0: free
          * 1: basic
          * 2: nonbasic at upper bound
          * 3: nonbasic at lower bound

        OUTPUT:

        - ``cstat`` -- The status of the column variables

        - ``rstat`` -- The status of the row variables

        .. NOTE::

            Logical variables associated with rows are all assumed to have +1
            coefficients, so for a <= constraint the logical will be at lower
            bound if the constraint is tight.

            Behaviour is undefined unless ``solve`` or ``set_basis_status`` 
            has been called before.


        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")                 # optional - cbc
            sage: p.add_variables(2)                              # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 2), (1, 3)], None, 6) # optional - cbc
            sage: p.add_linear_constraint([(0, 3), (1, 2)], None, 6) # optional - cbc
            sage: p.set_objective([1, 1], 7)                      # optional - cbc
            sage: p.solve()                                       # optional - cbc
            0
            sage: p.get_basis_status()                            # optional - cbc
            ([1, 1], [3, 3])

            sage: p = get_solver(solver = "Coin")                  # optional - cbc
            sage: p.add_variables(2)                               # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 2), (1, -3)], None, 6) # optional - cbc
            sage: p.add_linear_constraint([(0, 3), (1, 2)], None, 6)  # optional - cbc
            sage: p.set_objective([1, 1])                          # optional - cbc
            sage: p.solve()                                        # optional - cbc
            0
            sage: p.get_basis_status()                             # optional - cbc
            ([3, 1], [1, 3])

            sage: p = get_solver(solver = "Coin")                 # optional - cbc
            sage: p.add_variables(3)                              # optional - cbc
            2
            sage: p.add_linear_constraint(zip([0, 1, 2], [8, 6, 1]), None, 48)          # optional - cbc
            sage: p.add_linear_constraint(zip([0, 1, 2], [4, 2, 1.5]), None, 20)        # optional - cbc
            sage: p.add_linear_constraint(zip([0, 1, 2], [2, 1.5, 0.5]), None, 8)       # optional - cbc
            sage: p.set_objective([60, 30, 20])                   # optional - cbc
            sage: p.solve()                                       # optional - cbc
            0
            sage: p.get_basis_status()                            # optional - cbc
            ([1, 3, 1], [1, 3, 3])


            sage: lp = MixedIntegerLinearProgram(solver='Coin')   # optional - cbc
            sage: v = lp.new_variable(nonnegative=True)           # optional - cbc
            sage: x,y,z = v[0], v[1], v[2]                        # optional - cbc
            sage: lp.add_constraint(8*x + 6*y + z, max = 48)      # optional - cbc
            sage: lp.add_constraint(4*x + 2*y + 1.5*z, max = 20)  # optional - cbc
            sage: lp.add_constraint(2*x + 1.5*y + 0.5*z, max = 8) # optional - cbc
            sage: lp.set_objective(60*x + 30*y + 20*z)            # optional - cbc
            sage: lp_coin = lp.get_backend()                      # optional - cbc
            sage: lp_coin.solve()                                 # optional - cbc
            0
            sage: lp_coin.get_basis_status()                      # optional - cbc
            ([1, 3, 1], [1, 3, 3])

        """
        cdef int n = self.model.solver().getNumCols()
        cdef int m = self.model.solver().getNumRows()
        cdef int * c_cstat = <int *>check_malloc(n * sizeof(int))
        cdef int * c_rstat = <int *>check_malloc(m * sizeof(int))
        cdef list cstat
        cdef list rstat
        # enableSimplexInterface must be set to use getBasisStatus().
        # See projects.coin-or.org/Osi/ticket/84
        self.model.solver().enableSimplexInterface(True)
        try:
            sig_on()            # To catch SIGABRT
            self.model.solver().getBasisStatus(c_cstat, c_rstat)
            sig_off()
        except RuntimeError:    # corresponds to SIGABRT
            raise MIPSolverException('CBC : Signal sent, getBasisStatus() fails')
        else:
            cstat = [c_cstat[j] for j in range(n)]
            rstat = [c_rstat[j] for j in range(m)]
            return (cstat, rstat)
        finally:
            sig_free(c_cstat)
            sig_free(c_rstat)

    cpdef int set_basis_status(self, list cstat, list rstat) except -1:
        """
        Set the status of column and row variables
        and update the basis factorization and solution.

        This method returns status as integer codes:

        INPUT:

        - ``cstat`` -- The status of the column variables

        - ``rstat`` -- The status of the row variables

        .. NOTE::

            Status information should be coded as:

              * 0: free
              * 1: basic
              * 2: nonbasic at upper bound
              * 3: nonbasic at lower bound

            Logical variables associated with rows are all assumed to have +1
            coefficients, so for a <= constraint the logical will be at lower
            bound if the constraint is tight.

        OUTPUT:

        Returns 0 if all goes well, 1 if something goes wrong.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver

            sage: p = get_solver(solver = "Coin")                  # optional - cbc
            sage: p.add_variables(2)                               # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 2), (1, -3)], None, 6) # optional - cbc
            sage: p.add_linear_constraint([(0, 3), (1, 2)], None, 6)  # optional - cbc
            sage: p.set_objective([1, 1])                          # optional - cbc

            sage: p.set_basis_status([3, 3], [1, 1])               # optional - cbc
            0
            sage: p.get_objective_value()                          # optional - cbc
            0.0
            sage: p.set_basis_status([1, 3], [1, 3])               # optional - cbc
            0
            sage: p.get_objective_value()                          # optional - cbc
            2.0
            sage: p.set_basis_status([3, 1], [1, 3])               # optional - cbc
            0
            sage: p.get_objective_value()                          # optional - cbc
            3.0
            sage: p.get_basis_status()                             # optional - cbc
            ([3, 1], [1, 3])

            sage: p = get_solver(solver = "Coin")                 # optional - cbc
            sage: p.add_variables(3)                              # optional - cbc
            2
            sage: p.add_linear_constraint(zip([0, 1, 2], [8, 6, 1]), None, 48)          # optional - cbc
            sage: p.add_linear_constraint(zip([0, 1, 2], [4, 2, 1.5]), None, 20)        # optional - cbc
            sage: p.add_linear_constraint(zip([0, 1, 2], [2, 1.5, 0.5]), None, 8)       # optional - cbc
            sage: p.set_objective([60, 30, 20])                   # optional - cbc
            sage: p.set_basis_status([3, 3, 3], [1, 1, 1])        # optional - cbc
            0
            sage: p.get_objective_value()                         # optional - cbc
            0.0
            sage: p.set_basis_status([1, 3, 3], [1, 1, 3])        # optional - cbc
            0
            sage: p.get_objective_value()                         # optional - cbc
            240.0
            sage: p.get_basis_status()                            # optional - cbc
            ([1, 3, 3], [1, 1, 3])
            sage: p.set_basis_status([1, 3, 1], [1, 3, 2])        # optional - cbc
            0
            sage: p.get_basis_status()                            # optional - cbc
            ([1, 3, 1], [1, 3, 3])
            sage: p.get_objective_value()                         # optional - cbc
            280.0
        """
        cdef int n = len(cstat)
        cdef int m = len(rstat)
        cdef int * c_cstat
        cdef int * c_rstat
        cdef int result

        # set up the model
        cdef OsiSolverInterface * si = self.si

        cdef CbcModel * model
        cdef int old_logLevel = self.model.logLevel()

        model = new CbcModel(si[0])
        del self.model
        self.model = model
        
        #we immediately commit to the new model so that the user has access
        #to it even when something goes wrong.

        model.setLogLevel(old_logLevel)

        # multithreading
        import multiprocessing
        model.setNumberThreads(multiprocessing.cpu_count())
        
        if n != self.model.solver().getNumCols() or m != self.model.solver().getNumRows():
            raise ValueError("Must provide the status of every column and row variables")
        c_cstat = <int *>check_malloc(n * sizeof(int))
        c_rstat = <int *>check_malloc(m * sizeof(int))
        for i in range(n):
            c_cstat[i] = cstat[i]
        for i in range(m):
            c_rstat[i] = rstat[i]
        # enableSimplexInterface must be set to use getBasisStatus().
        # See projects.coin-or.org/Osi/ticket/84
        self.model.solver().enableSimplexInterface(True)
        try:
            sig_on()            # To catch SIGABRT
            result = self.model.solver().setBasisStatus(c_cstat, c_rstat)
            self.model.solver().setIntParam(OsiMaxNumIteration, 0)
            self.model.solver().resolve()
            sig_off()
        except RuntimeError:    # corresponds to SIGABRT
            raise MIPSolverException('CBC : Signal sent, setBasisStatus() fails')
        else:
            return result
        finally:
            sig_free(c_cstat)
            sig_free(c_rstat)

    cpdef get_binva_row(self, int i):
        """
        Return the i-th row of the tableau and the slacks.

        .. NOTE::

           Has no meaning unless ``solve`` or ``set_basis_status`` 
           has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver

            sage: p = get_solver(solver = "Coin")                  # optional - cbc
            sage: p.add_variables(2)                               # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 2), (1, -3)], None, 6) # optional - cbc
            sage: p.add_linear_constraint([(0, 3), (1, 2)], None, 6)  # optional - cbc
            sage: p.set_objective([1, 1])                          # optional - cbc

            sage: p.set_basis_status([3, 3], [1, 1])               # optional - cbc
            0
            sage: p.get_binva_row(0)                               # optional - cbc
            ([2.0, -3.0], [1.0, 0.0])
            sage: p.get_binva_row(1)                               # optional - cbc
            ([3.0, 2.0], [0.0, 1.0])

            sage: p.set_basis_status([1, 3], [1, 3])               # optional - cbc
            0
            sage: p.get_binva_row(0)                               # optional - cbc
            ([0.0, -4.333333333333333], [1.0, -0.6666666666666666])
            sage: p.get_binva_row(1)                               # optional - cbc
            ([1.0, 0.6666666666666666], [0.0, 0.3333333333333333])

            sage: p.set_basis_status([3, 1], [1, 3])               # optional - cbc
            0
            sage: p.get_binva_row(0)                               # optional - cbc
            ([6.5, 0.0], [1.0, 1.5])
            sage: p.get_binva_row(1)                               # optional - cbc
            ([1.5, 1.0], [0.0, 0.5])

        """
        cdef int n = self.model.solver().getNumCols()
        cdef int m = self.model.solver().getNumRows()
        if i < 0 or i >= m:
            raise ValueError("i = %s. The i-th row of the tableau doesn't exist" % i)
        
        cdef double * c_slack = <double *>check_malloc(m * sizeof(double))
        cdef double * c_z = <double *>check_malloc(n * sizeof(double))
        cdef list slack
        cdef list ithrow
        # enableSimplexInterface must be set to use getBasisStatus().
        # See projects.coin-or.org/Osi/ticket/84
        self.model.solver().enableSimplexInterface(True)
        try:
            sig_on()            # To catch SIGABRT
            self.model.solver().getBInvARow(i, c_z, c_slack)
            sig_off()
        except RuntimeError:    # corresponds to SIGABRT
            raise MIPSolverException('CBC : Signal sent, getBinvARow() fails')
        else:
            slack = [c_slack[j] for j in range(m)]
            ithrow = [c_z[j] for j in range(n)]
            return (ithrow, slack)
        finally:
            sig_free(c_slack)
            sig_free(c_z)

    cpdef get_binva_col(self, int j):
        """
        Return the j-th column of the tableau.
       
        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver

            sage: p = get_solver(solver = "Coin")                  # optional - cbc
            sage: p.add_variables(2)                               # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 2), (1, -3)], None, 6) # optional - cbc
            sage: p.add_linear_constraint([(0, 3), (1, 2)], None, 6)  # optional - cbc
            sage: p.set_objective([1, 1])                          # optional - cbc

            sage: p.set_basis_status([3, 3], [1, 1])               # optional - cbc
            0
            sage: p.get_binva_col(0)                               # optional - cbc
            [2.0, 3.0]
            sage: p.get_binva_col(1)                               # optional - cbc
            [-3.0, 2.0]

            sage: p.set_basis_status([1, 3], [1, 3])               # optional - cbc
            0
            sage: p.get_binva_col(0)                               # optional - cbc
            [-0.0, 1.0]
            sage: p.get_binva_col(1)                               # optional - cbc
            [-4.333333333333333, 0.6666666666666666]

            sage: p.set_basis_status([3, 1], [1, 3])               # optional - cbc
            0
            sage: p.get_binva_col(0)                               # optional - cbc
            [6.5, 1.5]
            sage: p.get_binva_col(1)                               # optional - cbc
            [-0.0, 1.0]
        """
        cdef int n = self.model.solver().getNumCols()
        cdef int m = self.model.solver().getNumRows()
        if j < 0 or j >= n + m:
            # it seems that when n <= j < m+n,
            # getBInvACol(j) is getBinvCol(j-n) 
            raise ValueError("j = %s. The j-th column of the tableau doesn't exist" % j)

        cdef double * c_vec = <double *>check_malloc(m * sizeof(double))
        cdef list jthcol
        # enableSimplexInterface must be set to use getBasisStatus().
        # See projects.coin-or.org/Osi/ticket/84
        self.model.solver().enableSimplexInterface(True)
        try:
            sig_on()            # To catch SIGABRT
            self.model.solver().getBInvACol(j, c_vec)
            sig_off()
        except RuntimeError:    # corresponds to SIGABRT
            raise MIPSolverException('CBC : Signal sent, getBinvACol() fails')
        else:
            jthcol = [c_vec[i] for i in range(m)]
            return jthcol
        finally:
            sig_free(c_vec)

    cpdef get_basics(self):
        r"""
        Returns indices of basic variables.

        The order of indices match the order of elements in the vectors returned 
        by get_binva_col() and the order of rows in get_binva_row(). 

        .. NOTE::

            Has no meaning unless ``solve`` or ``set_basis_status`` 
            has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")                  # optional - cbc
            sage: p.add_variables(2)                               # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 2), (1, -3)], None, 6) # optional - cbc
            sage: p.add_linear_constraint([(0, 3), (1, 2)], None, 6)  # optional - cbc
            sage: p.set_objective([1, 1])                          # optional - cbc
            sage: p.solve()                                        # optional - cbc
            0
            sage: p.get_basics()                                   # optional - cbc
            [2, 1]
        """
        cdef int m = self.model.solver().getNumRows()
        cdef int * c_indices = <int *>check_malloc(m * sizeof(int))
        cdef list indices 
        self.model.solver().enableSimplexInterface(True)
        try:
            sig_on()            # To catch SIGABRT
            self.model.solver().getBasics(c_indices)
            sig_off()
        except RuntimeError:    # corresponds to SIGABRT
            raise MIPSolverException('CBC : Signal sent, getBasics() fails')
        else:
            indices = [c_indices[j] for j in range(m)]
            return indices 
        finally:
            sig_free(c_indices)

    cpdef get_row_price(self):
        r"""
        Returns dual variable values.

        .. NOTE::

            Has no meaning unless ``solve`` or ``set_basis_status`` 
            has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")                  # optional - cbc
            sage: p.add_variables(2)                               # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 2), (1, -3)], None, 6) # optional - cbc
            sage: p.add_linear_constraint([(0, 3), (1, 2)], None, 6)  # optional - cbc
            sage: p.set_objective([1, 1])                          # optional - cbc
            sage: p.solve()                                        # optional - cbc
            0
            sage: p.get_row_price()                                # optional - cbc
            [0.0, -0.5]
        """
        cdef int m = self.model.solver().getNumRows()
        cdef list price
        cdef double * c_price
        self.model.solver().enableSimplexInterface(True)
        try:
            sig_on()            # To catch SIGABRT
            c_price = <double*>self.model.solver().getRowPrice()
            sig_off()
        except RuntimeError:    # corresponds to SIGABRT
            raise MIPSolverException('CBC : Signal sent, getRowPrice() fails')
        else:
            price = [c_price[j] for j in range(m)]
            return price

    cpdef get_reduced_cost(self):
        r"""
        Returns reduced costs.

        .. NOTE::

            Has no meaning unless ``solve`` or ``set_basis_status`` 
            has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Coin")                  # optional - cbc
            sage: p.add_variables(2)                               # optional - cbc
            1
            sage: p.add_linear_constraint([(0, 2), (1, -3)], None, 6) # optional - cbc
            sage: p.add_linear_constraint([(0, 3), (1, 2)], None, 6)  # optional - cbc
            sage: p.set_objective([1, 1])                          # optional - cbc
            sage: p.solve()                                        # optional - cbc
            0
            sage: p.get_reduced_cost()                             # optional - cbc
            [0.5, 0.0]
        """
        cdef int n = self.model.solver().getNumCols()
        cdef list cost
        cdef double * c_cost
        self.model.solver().enableSimplexInterface(True)
        try:
            sig_on()            # To catch SIGABRT
            c_cost = <double*>self.model.solver().getReducedCost()
            sig_off()
        except RuntimeError:    # corresponds to SIGABRT
            raise MIPSolverException('CBC : Signal sent, getReducedCost() fails')
        else:
            cost = [c_cost[i] for i in range(n)]
            return cost
