r"""
CVXOPT Backend


AUTHORS:

- Ingolfur Edvardsson (2014-05)        : initial implementation

"""
#*****************************************************************************
#       Copyright (C) 2014 Ingolfur Edvardsson <ingolfured@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.numerical.mip import MIPSolverException
from .generic_backend cimport GenericBackend
from copy import copy


cdef class CVXOPTBackend(GenericBackend):
    """
    MIP Backend that uses the CVXOPT solver.

    There is no support for integer variables.

    EXAMPLES::

        sage: p = MixedIntegerLinearProgram(solver="CVXOPT")        # optional - cvxopt

    TESTS:

    :trac:`20332`::

        sage: p                                                     # optional - cvxopt
        Mixed Integer Program (no objective, 0 variables, 0 constraints)
    """

    cdef list objective_function #c_matrix
    cdef list G_matrix
    cdef str prob_name
    cdef bint is_maximize

    cdef list row_lower_bound
    cdef list row_upper_bound
    cdef list col_lower_bound
    cdef list col_upper_bound

    cdef list row_name_var
    cdef list col_name_var
    cdef dict answer
    cdef dict param

    def __cinit__(self, maximization = True):
        """
        Cython constructor

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt

        """

        self.objective_function = [] #c_matrix in the example for cvxopt
        self.G_matrix = []
        self.prob_name = ''
        self.obj_constant_term = 0
        self.is_maximize = 1

        self.row_lower_bound = []
        self.row_upper_bound = []
        self.col_lower_bound = []
        self.col_upper_bound = []

        self.row_name_var = []
        self.col_name_var = []

        self.param = {"show_progress":False,
                      "maxiters":100,
                      "abstol":1e-7,
                      "reltol":1e-6,
                      "feastol":1e-7,
                      "refinement":0 }
        self.answer = {}
        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

    cpdef __copy__(self):
        # Added a second inequality to this doctest,
        # because otherwise CVXOPT complains: ValueError: Rank(A) < p or Rank([G; A]) < n
        """
        Returns a copy of self.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "CVXOPT")  # optional - cvxopt
            sage: b = p.new_variable()                              # optional - cvxopt
            sage: p.add_constraint(b[1] + b[2] <= 6)                # optional - cvxopt
            sage: p.add_constraint(b[2] <= 5)                       # optional - cvxopt
            sage: p.set_objective(b[1] + b[2])                      # optional - cvxopt
            sage: cp = copy(p.get_backend())                        # optional - cvxopt
            sage: cp.solve()                                        # optional - cvxopt
            0
            sage: cp.get_objective_value()                          # optional - cvxopt
            6.0
        """
        cdef CVXOPTBackend cp = type(self)()
        cp.objective_function = self.objective_function[:]
        cp.G_matrix = [row[:] for row in self.G_matrix]
        cp.prob_name = self.prob_name
        cp.obj_constant_term = self.obj_constant_term
        cp.is_maximize = self.is_maximize

        cp.row_lower_bound = self.row_lower_bound[:]
        cp.row_upper_bound = self.row_upper_bound[:]
        cp.col_lower_bound = self.col_lower_bound[:]
        cp.col_upper_bound = self.col_upper_bound[:]

        cp.row_name_var = self.row_name_var[:]
        cp.col_name_var = self.col_name_var[:]

        cp.param = copy(self.param)
        return cp

    cpdef int add_variable(self, lower_bound=0.0, upper_bound=None, binary=False, continuous=True, integer=False, obj=None, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.
        Variable types are always continuous, and thus the parameters
        ``binary``, ``integer``, and ``continuous`` have no effect.

        INPUT:

        - ``lower_bound`` - the lower bound of the variable (default: 0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is continuous (default: ``True``).

        - ``integer`` - ``True`` if the variable is integer (default: ``False``).

        - ``obj`` - (optional) coefficient of this variable in the objective function (default: 0.0)

        - ``name`` - an optional name for the newly added variable (default: ``None``).

        OUTPUT: The index of the newly created variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.ncols()                                         # optional - cvxopt
            0
            sage: p.add_variable()                                  # optional - cvxopt
            0
            sage: p.ncols()                                         # optional - cvxopt
            1
            sage: p.add_variable()                                  # optional - cvxopt
            1
            sage: p.add_variable(lower_bound=-2.0)                  # optional - cvxopt
            2
            sage: p.add_variable(continuous=True)                   # optional - cvxopt
            3
            sage: p.add_variable(name='x',obj=1.0)                  # optional - cvxopt
            4
            sage: p.col_name(3)                                     # optional - cvxopt
            'x_3'
            sage: p.col_name(4)                                     # optional - cvxopt
            'x'
            sage: p.objective_coefficient(4)                        # optional - cvxopt
            1.00000000000000

        TESTS::

            sage: p.add_variable(integer=True)                      # optional - cvxopt
            Traceback (most recent call last):
            ...
            RuntimeError: CVXOPT only supports continuous variables
            sage: p.add_variable(binary=True)                       # optional - cvxopt
            Traceback (most recent call last):
            ...
            RuntimeError: CVXOPT only supports continuous variables
        """
        if obj is None:
            obj = 0.0
        if binary or integer:
            raise RuntimeError("CVXOPT only supports continuous variables")
        self.G_matrix.append([0 for i in range(self.nrows())])
        self.col_lower_bound.append(lower_bound)
        self.col_upper_bound.append(upper_bound)
        self.objective_function.append(obj)
        self.col_name_var.append(name)
        return len(self.objective_function) - 1

    cpdef set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "cvxopt")                 # optional - cvxopt
            sage: p.add_variables(5)                                # optional - cvxopt
            4
            sage: p.set_variable_type(3, -1)                        # optional - cvxopt
            sage: p.set_variable_type(3, -2)                        # optional - cvxopt
            Traceback (most recent call last):
            ...
            ValueError: ...
        """
        if vtype != -1:
            raise ValueError('This backend does not handle integer variables ! Read the doc !')

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.is_maximization()                               # optional - cvxopt
            True
            sage: p.set_sense(-1)                                   # optional - cvxopt
            sage: p.is_maximization()                               # optional - cvxopt
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

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.add_variable()                                  # optional - cvxopt
            0
            sage: p.objective_coefficient(0)                        # optional - cvxopt
            0.0
            sage: p.objective_coefficient(0,2)                      # optional - cvxopt
            sage: p.objective_coefficient(0)                        # optional - cvxopt
            2.0
        """
        if coeff is not None:
            self.objective_function[variable] = float(coeff)
        else:
            return self.objective_function[variable]

    cpdef set_objective(self, list coeff, d = 0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        - ``d`` (double) -- the constant term in the linear function (set to `0` by default)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver    # optional - cvxopt
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.add_variables(5)                                # optional - cvxopt
            4
            sage: p.set_objective([1, 1, 2, 1, 3])                  # optional - cvxopt
            sage: [p.objective_coefficient(x) for x in range(5)]    # optional - cvxopt
            [1, 1, 2, 1, 3]
        """
        for i in range(len(coeff)):
            self.objective_function[i] = coeff[i]
        obj_constant_term = d

    cpdef set_verbosity(self, int level):
        """
        Does not apply for the cvxopt solver
        """
        pass

    cpdef add_col(self, indices, coeffs):
        """
        Add a column.

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

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.ncols()                                         # optional - cvxopt
            0
            sage: p.nrows()                                         # optional - cvxopt
            0
            sage: p.add_linear_constraints(5, 0, None)              # optional - cvxopt
            sage: p.add_col(range(5), range(5))                     # optional - cvxopt
            sage: p.nrows()                                         # optional - cvxopt
            5
        """
        column = []
        for _ in indices:
            column.append(0.0)

        for idx, ind in enumerate(indices):
            column[ind] = coeffs[idx]

        self.G_matrix.append(column)

        self.col_lower_bound.append(None)
        self.col_upper_bound.append(None)
        self.objective_function.append(0)
        self.col_name_var.append(None)

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

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.add_variables(5)                                # optional - cvxopt
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2.0, 2.0)                # optional - cvxopt
            sage: p.row(0)                                          # optional - cvxopt
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)                                   # optional - cvxopt
            (2.00000000000000, 2.00000000000000)
            sage: p.add_linear_constraint(zip(range(5), range(5)), 1.0, 1.0, name='foo')    # optional - cvxopt
            sage: p.row_name(-1)                                    # optional - cvxopt
            'foo'
        """
        coefficients = list(coefficients)
        for c in coefficients:
            while c[0] > len(self.G_matrix) - 1:
                self.add_variable()
        for i in range(len(self.G_matrix)):
            self.G_matrix[i].append(0.0)
        for c in coefficients:
            self.G_matrix[c[0]][-1] = c[1]

        self.row_lower_bound.append(lower_bound)
        self.row_upper_bound.append(upper_bound)
        self.row_name_var.append(name)

    cpdef int solve(self) except -1:
        """
        Solve the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution cannot be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver = "cvxopt", maximization=False)  # optional - cvxopt
            sage: x=p.new_variable(nonnegative=True)                # optional - cvxopt
            sage: p.set_objective(-4*x[0] - 5*x[1])                 # optional - cvxopt
            sage: p.add_constraint(2*x[0] + x[1] <= 3)              # optional - cvxopt
            sage: p.add_constraint(2*x[1] + x[0] <= 3)              # optional - cvxopt
            sage: N(p.solve(), digits=2)                            # optional - cvxopt
            -9.0
            sage: p = MixedIntegerLinearProgram(solver = "cvxopt", maximization=False)  # optional - cvxopt
            sage: x=p.new_variable(nonnegative=True)                # optional - cvxopt
            sage: p.set_objective(x[0] + 2*x[1])                    # optional - cvxopt
            sage: p.add_constraint(-5*x[0] + x[1]  <=   7)          # optional - cvxopt
            sage: p.add_constraint(-5*x[0] + x[1]  >=   7)          # optional - cvxopt
            sage: p.add_constraint(x[0] + x[1] >= 26  )             # optional - cvxopt
            sage: p.add_constraint( x[0] >= 3)                      # optional - cvxopt
            sage: p.add_constraint( x[1] >= 4)                      # optional - cvxopt
            sage: N(p.solve(),digits=4)                             # optional - cvxopt
            48.83
            sage: p = MixedIntegerLinearProgram(solver = "cvxopt")  # optional - cvxopt
            sage: x=p.new_variable(nonnegative=True)                # optional - cvxopt
            sage: p.set_objective(x[0] + x[1] + 3*x[2])             # optional - cvxopt
            sage: p.solver_parameter("show_progress",True)          # optional - cvxopt
            sage: p.add_constraint(x[0] + 2*x[1] <= 4)              # optional - cvxopt
            sage: p.add_constraint(5*x[2] - x[1] <= 8)              # optional - cvxopt
            sage: N(p.solve(), digits=2)                            # optional - cvxopt
                     pcost       dcost       gap    pres   dres   k/t
                 ...
                8.8
            sage: #CVXOPT gives different  values for variables compared to the other solvers.
            sage: c = MixedIntegerLinearProgram(solver = "cvxopt")  # optional - cvxopt
            sage: p = MixedIntegerLinearProgram(solver = "ppl")     # optional - cvxopt
            sage: g = MixedIntegerLinearProgram()                   # optional - cvxopt
            sage: xc=c.new_variable(nonnegative=True)               # optional - cvxopt
            sage: xp=p.new_variable(nonnegative=True)               # optional - cvxopt
            sage: xg=g.new_variable(nonnegative=True)               # optional - cvxopt
            sage: c.set_objective(xc[2])                            # optional - cvxopt
            sage: p.set_objective(xp[2])                            # optional - cvxopt
            sage: g.set_objective(xg[2])                            # optional - cvxopt
            sage: #we create a cube for all three solvers
            sage: c.add_constraint(xc[0] <= 100)                    # optional - cvxopt
            sage: c.add_constraint(xc[1] <= 100)                    # optional - cvxopt
            sage: c.add_constraint(xc[2] <= 100)                    # optional - cvxopt
            sage: p.add_constraint(xp[0] <= 100)                    # optional - cvxopt
            sage: p.add_constraint(xp[1] <= 100)                    # optional - cvxopt
            sage: p.add_constraint(xp[2] <= 100)                    # optional - cvxopt
            sage: g.add_constraint(xg[0] <= 100)                    # optional - cvxopt
            sage: g.add_constraint(xg[1] <= 100)                    # optional - cvxopt
            sage: g.add_constraint(xg[2] <= 100)                    # optional - cvxopt
            sage: N(c.solve(),digits=4)                             # optional - cvxopt
            100.0
            sage: N(c.get_values(xc[0]),digits=3)                   # optional - cvxopt
            50.0
            sage: N(c.get_values(xc[1]),digits=3)                   # optional - cvxopt
            50.0
            sage: N(c.get_values(xc[2]),digits=4)                   # optional - cvxopt
            100.0
            sage: N(p.solve(),digits=4)                             # optional - cvxopt
            100.0
            sage: N(p.get_values(xp[0]),2)                          # optional - cvxopt
            0.00
            sage: N(p.get_values(xp[1]),2)                          # optional - cvxopt
            0.00
            sage: N(p.get_values(xp[2]),digits=4)                   # optional - cvxopt
            100.0
            sage: N(g.solve(),digits=4)                             # optional - cvxopt
            100.0
            sage: N(g.get_values(xg[0]),2)                          # optional - cvxopt
            0.00
            sage: N(g.get_values(xg[1]),2)                          # optional - cvxopt
            0.00
            sage: N(g.get_values(xg[2]),digits=4)                   # optional - cvxopt
            100.0
        """
        from cvxopt import matrix, solvers
        h = []

        #for the equation bounds
        for eq_index in range(self.nrows()):
            h.append(self.row_upper_bound[eq_index])
            #upper bound is already in G
            if self.row_lower_bound[eq_index] is not None:
                h.append(-1 * self.row_lower_bound[eq_index])
                for cindex in range(len(self.G_matrix)):
                    if cindex == eq_index:
                        self.G_matrix[cindex].append(-1) # after multiplying the eq by -1
                    else:
                        self.G_matrix[cindex].append(0)



        #for the upper bounds (if there are any)
        for i in range(len(self.col_upper_bound)):
            if self.col_upper_bound[i] is not None:
                h.append(self.col_upper_bound[i])
                for cindex in range(len(self.G_matrix)):
                    if cindex == i:
                        self.G_matrix[cindex].append(1)
                    else:
                        self.G_matrix[cindex].append(0)
            if self.col_lower_bound[i] is not None:
                h.append(self.col_lower_bound[i])
                for cindex in range(len(self.G_matrix)):
                    if cindex == i:
                        self.G_matrix[cindex].append(-1) # after multiplying the eq by -1
                    else:
                        self.G_matrix[cindex].append(0)

        G = []
        for col in self.G_matrix:
            tempcol = []
            for i in range(len(col)):
                tempcol.append( float(col[i]) )
            G.append(tempcol)

        G = matrix(G)

        #cvxopt minimizes on default
        if self.is_maximize:
            c = [-1 * float(e) for e in self.objective_function]
        else:
            c = [float(e) for e in self.objective_function]
        c = matrix(c)

        h = [float(e) for e in h]
        h = matrix(h)

        #solvers comes from the cvxopt library
        for k,v in self.param.iteritems():
            solvers.options[k] = v
        self.answer = solvers.lp(c,G,h)

        #possible outcomes
        if self.answer['status'] == 'optimized':
            pass
        elif self.answer['status'] == 'primal infeasible':
            raise MIPSolverException("CVXOPT: primal infeasible")
        elif self.answer['status'] == 'dual infeasible':
            raise MIPSolverException("CVXOPT: dual infeasible")
        elif self.answer['status'] == 'unknown':
            raise MIPSolverException("CVXOPT: Terminated early due to numerical difficulties or because the maximum number of iterations was reached.")
        return 0


    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "cvxopt")                 # optional - cvxopt
            sage: p.add_variables(2)                                # optional - cvxopt
            1
            sage: p.add_linear_constraint([(0,1), (1,2)], None, 3)  # optional - cvxopt
            sage: p.set_objective([2, 5])                           # optional - cvxopt
            sage: p.solve()                                         # optional - cvxopt
            0
            sage: N(p.get_objective_value(),4)                      # optional - cvxopt
            7.5
            sage: N(p.get_variable_value(0),4)                      # optional - cvxopt
            3.6e-7
            sage: N(p.get_variable_value(1),4)                      # optional - cvxopt
            1.5
        """
        sum = self.obj_constant_term
        i = 0
        for v in self.objective_function:
            sum += v * float(self.answer['x'][i])
            i+=1
        return sum

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.add_variables(2)                                # optional - cvxopt
            1
            sage: p.add_linear_constraint([(0,1), (1, 2)], None, 3) # optional - cvxopt
            sage: p.set_objective([2, 5])                           # optional - cvxopt
            sage: p.solve()                                         # optional - cvxopt
            0
            sage: N(p.get_objective_value(),4)                      # optional - cvxopt
            7.5
            sage: N(p.get_variable_value(0),4)                      # optional - cvxopt
            3.6e-7
            sage: N(p.get_variable_value(1),4)                      # optional - cvxopt
            1.5
        """
        return self.answer['x'][variable]

    cpdef int ncols(self):
        """
        Return the number of columns/variables.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.ncols()                                         # optional - cvxopt
            0
            sage: p.add_variables(2)                                # optional - cvxopt
            1
            sage: p.ncols()                                         # optional - cvxopt
            2
        """

        return len(self.objective_function)

    cpdef int nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.nrows()                                         # optional - cvxopt
            0
            sage: p.add_variables(5)                                # optional - cvxopt
            4
            sage: p.add_linear_constraints(2, 2.0, None)            # optional - cvxopt
            sage: p.nrows()                                         # optional - cvxopt
            2
        """
        return len(self.row_upper_bound)


    cpdef bint is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.is_maximization()                               # optional - cvxopt
            True
            sage: p.set_sense(-1)                                   # optional - cvxopt
            sage: p.is_maximization()                               # optional - cvxopt
            False
        """
        if self.is_maximize == 1:
            return 1
        else:
            return 0

    cpdef problem_name(self, name=None):
        """
        Return or define the problem's name

        INPUT:

        - ``name`` (``str``) -- the problem's name. When set to
          ``None`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.problem_name()                                  # optional - cvxopt
            ''
            sage: p.problem_name("There once was a french fry")     # optional - cvxopt
            sage: print(p.problem_name())                           # optional - cvxopt
            There once was a french fry
        """
        if name is None:
            return self.prob_name
        self.prob_name = name


    cpdef row(self, int i):
        """
        Return a row

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_linear_constraint`` method.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.add_variables(5)                                # optional - cvxopt
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)  # optional - cvxopt
            sage: p.row(0)                                          # optional - cvxopt
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)                                   # optional - cvxopt
            (2, 2)
        """
        coeff = []
        idx = []
        index = 0
        for col in self.G_matrix:
            if col[i] != 0:
                idx.append(index)
                coeff.append(col[i])
            index += 1
        return (idx, coeff)


    cpdef row_bounds(self, int index):
        """
        Return the bounds of a specific constraint.

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the constraint is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.add_variables(5)                                # optional - cvxopt
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)  # optional - cvxopt
            sage: p.row(0)                                          # optional - cvxopt
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)                                   # optional - cvxopt
            (2, 2)
        """
        return (self.row_lower_bound[index], self.row_upper_bound[index])

    cpdef col_bounds(self, int index):
        """
        Return the bounds of a specific variable.

        INPUT:

        - ``index`` (integer) -- the variable's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the variable is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.add_variable()                                  # optional - cvxopt
            0
            sage: p.col_bounds(0)                                   # optional - cvxopt
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                      # optional - cvxopt
            sage: p.col_bounds(0)                                   # optional - cvxopt
            (0.0, 5)
        """
        return (self.col_lower_bound[index], self.col_upper_bound[index])

    cpdef bint is_variable_binary(self, int index):
        """
        Test whether the given variable is of binary type.
        CVXOPT does not allow integer variables, so this is a bit moot.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.ncols()                                         # optional - cvxopt
            0
            sage: p.add_variable()                                  # optional - cvxopt
            0
            sage: p.set_variable_type(0,0)                          # optional - cvxopt
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.is_variable_binary(0)                           # optional - cvxopt
            False

        """
        return False

    cpdef bint is_variable_integer(self, int index):
        """
        Test whether the given variable is of integer type.
        CVXOPT does not allow integer variables, so this is a bit moot.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.ncols()                                         # optional - cvxopt
            0
            sage: p.add_variable()                                  # optional - cvxopt
            0
            sage: p.set_variable_type(0,-1)                         # optional - cvxopt
            sage: p.set_variable_type(0,1)                          # optional - cvxopt
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.is_variable_integer(0)                          # optional - cvxopt
            False
        """
        return False

    cpdef bint is_variable_continuous(self, int index):
        """
        Test whether the given variable is of continuous/real type.
        CVXOPT does not allow integer variables, so this is a bit moot.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.ncols()                                         # optional - cvxopt
            0
            sage: p.add_variable()                                  # optional - cvxopt
            0
            sage: p.is_variable_continuous(0)                       # optional - cvxopt
            True
            sage: p.set_variable_type(0,1)                          # optional - cvxopt
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.is_variable_continuous(0)                       # optional - cvxopt
            True

        """
        return True

    cpdef row_name(self, int index):
        """
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.add_linear_constraints(1, 2, None, names=["Empty constraint 1"])  # optional - cvxopt
            sage: p.row_name(0)                                     # optional - cvxopt
            'Empty constraint 1'
        """
        if self.row_name_var[index] is not None:
            return self.row_name_var[index]
        return "constraint_" + repr(index)

    cpdef col_name(self, int index):
        """
        Return the ``index`` th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        - ``name`` (``char *``) -- its name. When set to ``NULL``
          (default), the method returns the current name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.add_variable(name="I am a variable")            # optional - cvxopt
            0
            sage: p.col_name(0)                                     # optional - cvxopt
            'I am a variable'
        """
        if self.col_name_var[index] is not None:
            return self.col_name_var[index]
        return "x_" + repr(index)

    cpdef variable_upper_bound(self, int index, value = False):
        """
        Return or define the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.add_variable()                                  # optional - cvxopt
            0
            sage: p.col_bounds(0)                                   # optional - cvxopt
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                      # optional - cvxopt
            sage: p.col_bounds(0)                                   # optional - cvxopt
            (0.0, 5)
        """
        if value is not False:
            self.col_upper_bound[index] = value
        else:
            return self.col_upper_bound[index]

    cpdef variable_lower_bound(self, int index, value = False):
        """
        Return or define the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")                 # optional - cvxopt
            sage: p.add_variable()                                  # optional - cvxopt
            0
            sage: p.col_bounds(0)                                   # optional - cvxopt
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)                      # optional - cvxopt
            sage: p.col_bounds(0)                                   # optional - cvxopt
            (5, None)
        """
        if value is not False:
            self.col_lower_bound[index] = value
        else:
            return self.col_lower_bound[index]

    cpdef solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter

        INPUT:

        - ``name`` (string) -- the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value.

        .. NOTE::

           The list of available parameters is available at
           :meth:`~sage.numerical.mip.MixedIntegerLinearProgram.solver_parameter`.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="CVXOPT")                   # optional - cvxopt
            sage: p.solver_parameter("show_progress")               # optional - cvxopt
            False
            sage: p.solver_parameter("show_progress", True)         # optional - cvxopt
            sage: p.solver_parameter("show_progress")               # optional - cvxopt
            True
        """
        if value is None:
            return self.param[name]
        else:
            self.param[name] = value
