r"""
Matrix Backend for SDP solvers

It stores the SDP data in Sage matrices.  It allows users to specify a base ring
and can store either floating-point SDPs or exact SDPs with rational or algebraic data.

The class does not provide a solver method.  It can be used as a base class for
other classes implementing solvers.
"""

#*****************************************************************************
#       Copyright (C) 2014 Ingolfur Edvardsson <ingolfured@gmail.com>
#       Copyright (C) 2020 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import matrix
from sage.numerical.backends.generic_sdp_backend cimport GenericSDPBackend

cdef class MatrixSDPBackend(GenericSDPBackend):

    def __init__(self, maximization=True, base_ring=None):
        """
        Cython constructor.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
        """
        self.objective_function = []
        self.name = ""
        self.coeffs_matrix = []
        self.obj_constant_term = 0
        self.matrices_dim = {}
        self.is_maximize = 1

        self.row_name_var = []
        self.col_name_var = []

        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ
        self._base_ring = base_ring

    cpdef base_ring(self):
        """
        The base ring.

        TESTS::

            sage: from sage.numerical.backends.matrix_sdp_backend import MatrixSDPBackend
            sage: MatrixSDPBackend(base_ring=QQ).base_ring()
            Rational Field
        """
        return self._base_ring

    def get_matrix(self):
        """
        Get a block of a matrix coefficient.

        EXAMPLES::

            sage: p = SemidefiniteProgram(solver='cvxopt')
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

        - ``obj`` -- (optional) coefficient of this variable in the objective
          function (default: 0.0)

        - ``name`` -- an optional name for the newly added variable (default:
          ``None``)

        OUTPUT: the index of the newly created variable

        EXAMPLES::

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
            if self.matrices_dim.get(i) is not None:
                row.append( matrix.zero(self.matrices_dim[i], self.matrices_dim[i]) )
            else:
                row.append(0)
            i += 1
        self.col_name_var.append(name)
        self.objective_function.append(obj)
        return len(self.objective_function) - 1

    cpdef int add_variables(self, int n, names=None) except -1:
        """
        Add ``n`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        INPUT:

        - ``n`` -- the number of new variables (must be > 0)

        - ``names`` -- list of names (default: ``None``)

        OUTPUT: the index of the variable created last

        EXAMPLES::

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
        return len(self.objective_function) - 1

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` -- integer:

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLES::

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

        - ``variable`` -- integer; the variable's id

        - ``coeff`` -- double; its coefficient

        EXAMPLES::

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
            self.objective_function[variable] = float(coeff)
        else:
            return self.objective_function[variable]

    cpdef set_objective(self, list coeff, d=0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- list of real values, whose i-th element is the
          coefficient of the i-th variable in the objective function

        - ``d`` -- double; the constant term in the linear function (set to `0`
          by default)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1, 1, 2, 1, 3]
        """
        for i in range(len(coeff)):
            self.objective_function[i] = coeff[i]
        obj_constant_term = d

    cpdef add_linear_constraint(self, coefficients, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` an iterable with ``(c,v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (matrix).
          The pairs come sorted by indices. If c is -1 it
          represents the constant coefficient.

        - ``name`` -- an optional name for this row (default: ``None``)

        EXAMPLES::

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
            sage: p.add_linear_constraint(  [(0, matrix([[33., -9.], [-9., 26.]])) , (1,  matrix([[-7., -11.] ,[ -11., 3.]]) )],name='fun')
            sage: p.row_name(-1)
            'fun'
        """
        coefficients = list(coefficients)
        from sage.structure.element import Matrix
        for t in coefficients:
            m = t[1]
            if not isinstance(m, Matrix):
                raise ValueError("The coefficients must be matrices")
            if not m.is_square():
                raise ValueError("The matrix has to be a square")
            if self.matrices_dim.get(self.nrows()) is not None and m.dimensions()[0] != self.matrices_dim.get(self.nrows()):
                raise ValueError("the matrices have to be of the same dimension")
        self.coeffs_matrix.append(coefficients)
        self.matrices_dim[self.nrows()] = m.dimensions()[0] #
        self.row_name_var.append(name)

    cpdef add_linear_constraints(self, int number, names=None):
        """
        Add constraints.

        INPUT:

        - ``number`` -- integer; the number of constraints to add

        - ``names`` -- an optional list of names (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraints(5)
            sage: p.row(4)
            ([], [])
        """
        for i in range(number):
            self.add_linear_constraint(zip(range(self.ncols()+1),[matrix.zero(1,1) for i in range(self.ncols()+1)]),
                                       name=None if names is None else names[i])

    cpdef int ncols(self) noexcept:
        """
        Return the number of columns/variables.

        EXAMPLES::

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

    cpdef int nrows(self) noexcept:
        """
        Return the number of rows/constraints.

        EXAMPLES::

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

    cpdef bint is_maximization(self) noexcept:
        """
        Test whether the problem is a maximization

        EXAMPLES::

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

    cpdef problem_name(self, name=None):
        """
        Return or define the problem's name.

        INPUT:

        - ``name`` -- string; the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
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
            if m != matrix.zero(self.matrices_dim[i],self.matrices_dim[i]):
                indices.append(index)
                matrices.append(m)
        return (indices, matrices)

    cpdef row_name(self, int index):
        """
        Return the ``index``-th row name.

        INPUT:

        - ``index`` -- integer; the row's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver = "CVXOPT")
            sage: p.add_linear_constraints(1, names='A')
            sage: p.row_name(0)
            'A'
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
