# sage_setup: distribution = sagemath-polyhedra
# sage.doctest: needs cvxopt
r"""
CVXOPT SDP Backend


AUTHORS:

- Ingolfur Edvardsson (2014-05): initial implementation

- Dima Pasechnik      (2015-12): minor fixes
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

from sage.numerical.sdp import SDPSolverException
from sage.matrix.constructor import Matrix
from sage.numerical.backends.matrix_sdp_backend cimport MatrixSDPBackend


cdef class CVXOPTSDPBackend(MatrixSDPBackend):

    cdef dict answer
    cdef dict param

    def __init__(self, maximization=True, base_ring=None):
        """
        Cython constructor.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver='CVXOPT')
        """

        from sage.rings.real_double import RDF
        if base_ring is None:
            base_ring = RDF
        if base_ring is not RDF:
            raise ValueError("only base_ring=RDF is supported")
        MatrixSDPBackend.__init__(self, maximization, base_ring=base_ring)

        self.param = {"show_progress":False,
                      "maxiters":100,
                      "abstol":1e-7,
                      "reltol":1e-6,
                      "feastol":1e-7,
                      "refinement":1 }
        self.answer = {}

    cpdef int solve(self) except -1:
        """
        Solve the problem.

        .. NOTE::

            This method raises :class:`SDPSolverException` exceptions when
            the solution cannot be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLES::

            sage: p = SemidefiniteProgram(solver='cvxopt', maximization=False)
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
            sage: N(p.solve(), digits=4)
            -3.225
            sage: p = SemidefiniteProgram(solver='cvxopt', maximization=False)
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
            sage: N(p.solve(), digits=4)
            -3.154
        """
        from cvxopt import matrix as c_matrix, solvers
        from sage.rings.real_double import RDF
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

        EXAMPLES::

            sage: p = SemidefiniteProgram(solver='cvxopt', maximization=False)
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
            sage: N(p.solve(), digits=4)
            -3.154
            sage: N(p.get_backend().get_objective_value(), digits=4)
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
        Return the complete output dict of the solver.

        Mainly for testing purposes

        TESTS::

            sage: p = SemidefiniteProgram(maximization=False, solver='cvxopt')
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

        EXAMPLES::

            sage: p = SemidefiniteProgram(solver='cvxopt', maximization=False)
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
            sage: N(p.solve(), digits=4)
            -3.154
            sage: N(p.get_backend().get_variable_value(0), digits=3)
            -0.368
            sage: N(p.get_backend().get_variable_value(1), digits=4)
            1.898
            sage: N(p.get_backend().get_variable_value(2), digits=3)
            -0.888
        """
        return self.answer['x'][variable]

    cpdef dual_variable(self, int i, sparse=False):
        """
        The `i`-th dual variable.

        Available after ``self.solve()`` is called, otherwise the result is
        undefined.

        - ``index`` -- integer; the constraint's id

        OUTPUT: the matrix of the `i`-th dual variable

        EXAMPLES::

            sage: p = SemidefiniteProgram(maximization=False, solver='cvxopt')
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
            sage: B = p.get_backend()
            sage: x = p.get_values(x).values()
            sage: -(a3*B.dual_variable(0)).trace() - (b3*B.dual_variable(1)).trace()     # tol 1e-07
            -3.0
            sage: g = sum((B.slack(j)*B.dual_variable(j)).trace() for j in range(2)); g  # tol 1.5e-08
            0.0


        TESTS::

            sage: B.dual_variable(7)
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
        Slack of the `i`-th constraint.

        Available after ``self.solve()`` is called, otherwise the result is
        undefined.

        - ``index`` -- integer; the constraint's id

        OUTPUT: the matrix of the slack of the `i`-th constraint

        EXAMPLES::

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
            sage: x = sorted(p.get_values(x).values())
            sage: x[0]*b1 + x[1]*b2 - b3 + B1       # tol 1e-09
            [0.0 0.0]
            [0.0 0.0]

        TESTS::

            sage: B.slack(7)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        cdef int n
        n = self.answer['ss'][i].size[0]
        assert(n == self.answer['ss'][i].size[1]) # must be square matrix
        return Matrix(n, n, list(self.answer['ss'][i]), sparse=sparse)

    cpdef solver_parameter(self, name, value=None):
        """
        Return or define a solver parameter.

        INPUT:

        - ``name`` -- string; the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value

        .. NOTE::

           The list of available parameters is available at
           :meth:`~sage.numerical.sdp.SemidefiniteProgram.solver_parameter`.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver='CVXOPT')
            sage: p.solver_parameter("show_progress")
            False
            sage: p.solver_parameter("show_progress", True)
            sage: p.solver_parameter("show_progress")
            True
        """
        if value is None:
            return self.param[name]
        else:
            self.param[name] = value
