r"""
CVXOPT SDP Backend


AUTHORS:

- Ingolfur Edvardsson (2014-05)        : initial implementation

- Dima Pasechnik      (2015-12)        : minor fixes

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
from __future__ import print_function

from sage.numerical.sdp import SDPSolverException
from sage.matrix.all import Matrix
from cvxopt import solvers
from .generic_sdp_backend cimport GenericSDPBackend


cdef class CVXOPTSDPBackend(GenericSDPBackend):
    cdef list objective_function #c_matrix
    cdef list coeffs_matrix
    cdef bint is_maximize

    cdef list row_name_var
    cdef list col_name_var
    cdef dict answer
    cdef dict param
    cdef str name

    def __cinit__(self, maximization = True):
        """
        Cython constructor

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")

        """

        self.objective_function = [] #c_matrix in the example for cvxopt
        self.name = ""
        self.coeffs_matrix = []
        self.obj_constant_term = 0
        self.matrices_dim = {}
        self.is_maximize = 1

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

    def get_matrix(self):
        """
        Get a block of a matrix coefficient

        EXAMPLE::

            sage: p = SemidefiniteProgram(solver="cvxopt")
            sage: x = p.new_variable()
            sage: a1 = matrix([[1, 2.], [2., 3.]])
            sage: a2 = matrix([[3, 4.], [4., 5.]])
            sage: p.add_constraint(a1*x[0] + a2*x[1] <= a1)
            sage: b = p.get_backend()
            sage: b.get_matrix()[0][0]
            (
                [-1.0 -2.0]
            -1, [-2.0 -3.0]
            )

        """
        return self.coeffs_matrix

    cpdef int add_variable(self, obj=0.0,  name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column of matrices to the matrix. By default,
        the variable is both positive and real.

        INPUT:

        - ``obj`` - (optional) coefficient of this variable in the objective function (default: 0.0)

        - ``name`` - an optional name for the newly added variable (default: ``None``).

        OUTPUT: The index of the newly created variable

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable()
            1
            sage: p.add_variable(name='x',obj=1.0)
            2
            sage: p.col_name(2)
            'x'
            sage: p.objective_coefficient(2)
            1.00000000000000
        """
        i = 0
        for row in self.coeffs_matrix:
            if self.matrices_dim.get(i) != None:
                row.append( Matrix.zero(self.matrices_dim[i], self.matrices_dim[i]) )
            else:
                row.append(0)
            i+=1
        self.col_name_var.append(name)
        self.objective_function.append(obj)
        return len(self.objective_function) - 1


    cpdef int add_variables(self, int n, names=None) except -1:
        """
        Add ``n`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        INPUT:

        - ``n`` - the number of new variables (must be > 0)

        - ``names`` - optional list of names (default: ``None``)

        OUTPUT: The index of the variable created last.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.ncols()
            0
            sage: p.add_variables(5)
            4
            sage: p.ncols()
            5
            sage: p.add_variables(2, names=['a','b'])
            6
        """
        for i in range(n):
            self.add_variable()
        return len(self.objective_function) - 1;


    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
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

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0.0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2.0
        """
        if coeff is not None:
            self.objective_function[variable] = float(coeff);
        else:
            return self.objective_function[variable]

    cpdef set_objective(self, list coeff, d=0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        - ``d`` (double) -- the constant term in the linear function (set to `0` by default)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: map(lambda x :p.objective_coefficient(x), range(5))
            [1, 1, 2, 1, 3]
        """
        for i in range(len(coeff)):
            self.objective_function[i] = coeff[i];
        obj_constant_term = d;



    cpdef add_linear_constraint(self, coefficients, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` an iterable with ``(c,v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (matrix).
          The pairs come sorted by indices. If c is -1 it
          represents the constant coefficient.

        - ``name`` - an optional name for this row (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint(  [(0, matrix([[33., -9.], [-9., 26.]])) , (1,  matrix([[-7., -11.] ,[ -11., 3.]]) )])
            sage: p.row(0)
            ([0, 1],
             [
            [ 33.0000000000000 -9.00000000000000]
            [-9.00000000000000  26.0000000000000],
            <BLANKLINE>
            [-7.00000000000000 -11.0000000000000]
            [-11.0000000000000  3.00000000000000]
            ])
            sage: p.add_linear_constraint(  [(0, matrix([[33., -9.], [-9., 26.]])) , (1,  matrix([[-7., -11.] ,[ -11., 3.]]) )],name="fun")
            sage: p.row_name(-1)
            'fun'
        """
        from sage.matrix.matrix import is_Matrix
        for t in coefficients:
            m = t[1]
            if not is_Matrix(m):
                raise ValueError("The coefficients must be matrices")
            if not m.is_square():
                raise ValueError("The matrix has to be a square")
            if self.matrices_dim.get(self.nrows()) != None and m.dimensions()[0] != self.matrices_dim.get(self.nrows()):
                raise ValueError("The matrces have to be of the same dimension")
        self.coeffs_matrix.append(coefficients)
        self.matrices_dim[self.nrows()] = m.dimensions()[0] #
        self.row_name_var.append(name)

    cpdef add_linear_constraints(self, int number, names=None):
        """
        Add constraints.

        INPUT:

        - ``number`` (integer) -- the number of constraints to add.

        - ``names`` - an optional list of names (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraints(5)
            sage: p.row(4)
            ([], [])
        """
        for i in range(number):
            self.add_linear_constraint(zip(range(self.ncols()+1),[Matrix.zero(1,1) for i in range(self.ncols()+1)]),
                                       name=None if names is None else names[i])



    cpdef int solve(self) except -1:
        """
        Solve the problem.

        .. NOTE::

            This method raises ``SDPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: p = SemidefiniteProgram(solver = "cvxopt", maximization=False)
            sage: x = p.new_variable()
            sage: p.set_objective(x[0] - x[1] + x[2])
            sage: a1 = matrix([[-7., -11.], [-11., 3.]])
            sage: a2 = matrix([[7., -18.], [-18., 8.]])
            sage: a3 = matrix([[-2., -8.], [-8., 1.]])
            sage: a4 = matrix([[33., -9.], [-9., 26.]])
            sage: b1 = matrix([[-21., -11., 0.], [-11., 10., 8.], [0.,   8., 5.]])
            sage: b2 = matrix([[0.,  10.,  16.], [10., -10., -10.], [16., -10., 3.]])
            sage: b3 = matrix([[-5.,   2., -17.], [2.,  -6.,   8.], [-17.,  8., 6.]])
            sage: b4 = matrix([[14., 9., 40.], [9., 91., 10.], [40., 10., 15.]])
            sage: p.add_constraint(a1*x[0] + a3*x[2] <= a4)
            sage: p.add_constraint(b1*x[0] + b2*x[1] + b3*x[2] <= b4)
            sage: round(p.solve(), 3)
            -3.225
            sage: p = SemidefiniteProgram(solver = "cvxopt", maximization=False)
            sage: x = p.new_variable()
            sage: p.set_objective(x[0] - x[1] + x[2])
            sage: a1 = matrix([[-7., -11.], [-11., 3.]])
            sage: a2 = matrix([[7., -18.], [-18., 8.]])
            sage: a3 = matrix([[-2., -8.], [-8., 1.]])
            sage: a4 = matrix([[33., -9.], [-9., 26.]])
            sage: b1 = matrix([[-21., -11., 0.], [-11., 10., 8.], [0.,   8., 5.]])
            sage: b2 = matrix([[0.,  10.,  16.], [10., -10., -10.], [16., -10., 3.]])
            sage: b3 = matrix([[-5.,   2., -17.], [2.,  -6.,   8.], [-17.,  8., 6.]])
            sage: b4 = matrix([[14., 9., 40.], [9., 91., 10.], [40., 10., 15.]])
            sage: p.add_constraint(a1*x[0] + a2*x[1] + a3*x[2] <= a4)
            sage: p.add_constraint(b1*x[0] + b2*x[1] + b3*x[2] <= b4)
            sage: round(p.solve(), 3)
            -3.154

        """
        from cvxopt import matrix as c_matrix, solvers
        from sage.rings.all import RDF
        G_matrix = []
        h_matrix = []
        debug_g = []
        debug_h = []
        debug_c = []

        #cvxopt minimizes on default
        if self.is_maximize:
            c = [-1 * float(e) for e in self.objective_function]
        else:
            c = [float(e) for e in self.objective_function]
        debug_c = (c)
        c = c_matrix(c)

        row_index = -1
        for row in self.coeffs_matrix:
            row_index += 1
            row.sort()
            G_temp = []
            add_null = [True for i in range(self.ncols())]
            for i,m in row:
                if i == -1:
                    h_temp = []
                    for row in m.rows():
                        row_temp = []
                        for e in row:
                            row_temp.append(-1*float(e))
                        h_temp.append(row_temp)
                    h_matrix += [c_matrix(h_temp)]
                    debug_h += [h_temp]
                else:
                    add_null[i] = False
                    m = [float(e) for e in m.list()]
                    G_temp.append(m)
            for j in range(self.ncols()):
                if add_null[j]:
                    G_temp.insert(j,[float(0) for t in range(self.matrices_dim[row_index]**2)])
            G_matrix += [c_matrix(G_temp)]
            debug_g += [(G_temp)]
        #raise Exception("G_matrix " + str(debug_g) + "\nh_matrix: " + str(debug_h) + "\nc_matrix: " + str(debug_c))

        #solvers comes from the cvxopt library
        for k,v in self.param.iteritems():
            solvers.options[k] = v

        self.answer = solvers.sdp(c,Gs=G_matrix,hs=h_matrix)

        #possible outcomes
        if self.answer['status'] == 'optimized':
            pass
        elif self.answer['status'] == 'primal infeasible':
            raise SDPSolverException("CVXOPT: primal infeasible")
        elif self.answer['status'] == 'dual infeasible':
            raise SDPSolverException("CVXOPT: dual infeasible")
        elif self.answer['status'] == 'unknown':
            raise SDPSolverException("CVXOPT: Terminated early due to numerical difficulties or because the maximum number of iterations was reached.")
        return 0


    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: p = SemidefiniteProgram(solver = "cvxopt", maximization=False)
            sage: x = p.new_variable()
            sage: p.set_objective(x[0] - x[1] + x[2])
            sage: a1 = matrix([[-7., -11.], [-11., 3.]])
            sage: a2 = matrix([[7., -18.], [-18., 8.]])
            sage: a3 = matrix([[-2., -8.], [-8., 1.]])
            sage: a4 = matrix([[33., -9.], [-9., 26.]])
            sage: b1 = matrix([[-21., -11., 0.], [-11., 10., 8.], [0.,   8., 5.]])
            sage: b2 = matrix([[0.,  10.,  16.], [10., -10., -10.], [16., -10., 3.]])
            sage: b3 = matrix([[-5.,   2., -17.], [2.,  -6.,   8.], [-17.,  8., 6.]])
            sage: b4 = matrix([[14., 9., 40.], [9., 91., 10.], [40., 10., 15.]])
            sage: p.add_constraint(a1*x[0] + a2*x[1] + a3*x[2] <= a4)
            sage: p.add_constraint(b1*x[0] + b2*x[1] + b3*x[2] <= b4)
            sage: round(p.solve(),3)
            -3.154
            sage: round(p.get_backend().get_objective_value(),3)
            -3.154
        """
        sum = self.obj_constant_term
        i = 0
        for v in self.objective_function:
            sum += v * float(self.answer['x'][i])
            i+=1
        return sum

    cpdef _get_answer(self):
        """
        return the complete output dict of the solver

        Mainly for testing purposes

        TESTS::

            sage: p = SemidefiniteProgram(maximization = False, solver='cvxopt')
            sage: x = p.new_variable()
            sage: p.set_objective(x[0] - x[1])
            sage: a1 = matrix([[1, 2.], [2., 3.]])
            sage: a2 = matrix([[3, 4.], [4., 5.]])
            sage: a3 = matrix([[5, 6.], [6., 7.]])
            sage: b1 = matrix([[1, 1.], [1., 1.]])
            sage: b2 = matrix([[2, 2.], [2., 2.]])
            sage: b3 = matrix([[3, 3.], [3., 3.]])
            sage: p.add_constraint(a1*x[0] + a2*x[1] <= a3)
            sage: p.add_constraint(b1*x[0] + b2*x[1] <= b3)
            sage: p.solve();                                     # tol 1e-08
            -3.0
            sage: p.get_backend()._get_answer()
            {...}
        """
        return self.answer

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "cvxopt")
            sage: p = SemidefiniteProgram(solver = "cvxopt", maximization=False)
            sage: x = p.new_variable()
            sage: p.set_objective(x[0] - x[1] + x[2])
            sage: a1 = matrix([[-7., -11.], [-11., 3.]])
            sage: a2 = matrix([[7., -18.], [-18., 8.]])
            sage: a3 = matrix([[-2., -8.], [-8., 1.]])
            sage: a4 = matrix([[33., -9.], [-9., 26.]])
            sage: b1 = matrix([[-21., -11., 0.], [-11., 10., 8.], [0.,   8., 5.]])
            sage: b2 = matrix([[0.,  10.,  16.], [10., -10., -10.], [16., -10., 3.]])
            sage: b3 = matrix([[-5.,   2., -17.], [2.,  -6.,   8.], [-17.,  8., 6.]])
            sage: b4 = matrix([[14., 9., 40.], [9., 91., 10.], [40., 10., 15.]])
            sage: p.add_constraint(a1*x[0] + a2*x[1] + a3*x[2] <= a4)
            sage: p.add_constraint(b1*x[0] + b2*x[1] + b3*x[2] <= b4)
            sage: round(p.solve(),3)
            -3.154
            sage: round(p.get_backend().get_variable_value(0),3)
            -0.368
            sage: round(p.get_backend().get_variable_value(1),3)
            1.898
            sage: round(p.get_backend().get_variable_value(2),3)
            -0.888

        """
        return self.answer['x'][variable]


    cpdef int ncols(self):
        """
        Return the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.ncols()
            0
            sage: p.add_variables(2)
            1
            sage: p.ncols()
            2
        """

        return len(self.objective_function)

    cpdef int nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.nrows()
            0
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraints(2)
            sage: p.nrows()
            2
        """
        return len(self.matrices_dim)


    cpdef bint is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
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

    cpdef problem_name(self, char * name = NULL):
        """
        Return or define the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.problem_name("There once was a french fry")
            sage: print(p.problem_name())
            There once was a french fry
        """
        if name == NULL:
            return self.name
        self.name = str(<bytes>name)


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

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(  [(0, matrix([[33., -9.], [-9., 26.]])) , (1,  matrix([[-7., -11.] ,[ -11., 3.]]) )])
            sage: p.row(0)
            ([0, 1],
             [
            [ 33.0000000000000 -9.00000000000000]
            [-9.00000000000000  26.0000000000000],
            <BLANKLINE>
            [-7.00000000000000 -11.0000000000000]
            [-11.0000000000000  3.00000000000000]
            ])
        """
        indices = []
        matrices = []
        for index,m in self.coeffs_matrix[i]:
            if m != Matrix.zero(self.matrices_dim[i],self.matrices_dim[i]):
                indices.append(index)
                matrices.append(m)
        return (indices, matrices)

    cpdef dual_variable(self, int i, sparse=False):
        """
        The `i`-th dual variable

        Available after self.solve() is called, otherwise the result is undefined

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        The matrix of the `i`-th dual variable

        EXAMPLE::

            sage: p = SemidefiniteProgram(maximization = False, solver='cvxopt')
            sage: x = p.new_variable()
            sage: p.set_objective(x[0] - x[1])
            sage: a1 = matrix([[1, 2.], [2., 3.]])
            sage: a2 = matrix([[3, 4.], [4., 5.]])
            sage: a3 = matrix([[5, 6.], [6., 7.]])
            sage: b1 = matrix([[1, 1.], [1., 1.]])
            sage: b2 = matrix([[2, 2.], [2., 2.]])
            sage: b3 = matrix([[3, 3.], [3., 3.]])
            sage: p.add_constraint(a1*x[0] + a2*x[1] <= a3)
            sage: p.add_constraint(b1*x[0] + b2*x[1] <= b3)
            sage: p.solve()                                     # tol 1e-08
            -3.0
            sage: B=p.get_backend()
            sage: x=p.get_values(x).values()
            sage: -(a3*B.dual_variable(0)).trace()-(b3*B.dual_variable(1)).trace()  # tol 1e-07
            -3.0
            sage: g = sum((B.slack(j)*B.dual_variable(j)).trace() for j in range(2)); g  # tol 1.5e-08
            0.0


        TESTS::

            sage: B.dual_variable(7)
            ...
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: abs(g - B._get_answer()['gap'])   # tol 1e-22
            0.0

        """
        cdef int n
        n = self.answer['zs'][i].size[0]
        assert(n == self.answer['zs'][i].size[1]) # must be square matrix
        return Matrix(n, n, list(self.answer['zs'][i]), sparse=sparse)

    cpdef slack(self, int i, sparse=False):
        """
        Slack of the `i`-th constraint

        Available after self.solve() is called, otherwise the result is undefined

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        The matrix of the slack of the `i`-th constraint

        EXAMPLE::

            sage: p = SemidefiniteProgram(maximization = False, solver='cvxopt')
            sage: x = p.new_variable()
            sage: p.set_objective(x[0] - x[1])
            sage: a1 = matrix([[1, 2.], [2., 3.]])
            sage: a2 = matrix([[3, 4.], [4., 5.]])
            sage: a3 = matrix([[5, 6.], [6., 7.]])
            sage: b1 = matrix([[1, 1.], [1., 1.]])
            sage: b2 = matrix([[2, 2.], [2., 2.]])
            sage: b3 = matrix([[3, 3.], [3., 3.]])
            sage: p.add_constraint(a1*x[0] + a2*x[1] <= a3)
            sage: p.add_constraint(b1*x[0] + b2*x[1] <= b3)
            sage: p.solve()                         # tol 1e-08
            -3.0
            sage: B = p.get_backend()
            sage: B1 = B.slack(1); B1               # tol 1e-08
            [0.0 0.0]
            [0.0 0.0]
            sage: B1.is_positive_definite()
            True
            sage: x = p.get_values(x).values()
            sage: x[0]*b1 + x[1]*b2 - b3 + B1       # tol 1e-09
            [0.0 0.0]
            [0.0 0.0]

        TESTS::

            sage: B.slack(7)
            ...
            Traceback (most recent call last):
            ...
            IndexError: list index out of range

        """
        cdef int n
        n = self.answer['ss'][i].size[0]
        assert(n == self.answer['ss'][i].size[1]) # must be square matrix
        return Matrix(n, n, list(self.answer['ss'][i]), sparse=sparse)

    cpdef row_name(self, int index):
        """
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_linear_constraints(1, names="A")
            sage: p.row_name(0)
            'A'

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

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_variable(name="I am a variable")
            0
            sage: p.col_name(0)
            'I am a variable'
        """
        if self.col_name_var[index] is not None:
            return self.col_name_var[index]
        return "x_" + repr(index)



    cpdef solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter

        INPUT:

        - ``name`` (string) -- the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value.

        .. NOTE::

           The list of available parameters is available at
           :meth:`~sage.numerical.sdp.SemidefiniteProgram.solver_parameter`.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.solver_parameter("show_progress")
            False
            sage: p.solver_parameter("show_progress", True)
            sage: p.solver_parameter("show_progress")
            True
        """
        if value == None:
            return self.param[name]
        else:
            self.param[name] = value

