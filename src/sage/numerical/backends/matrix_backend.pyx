r"""
Matrix Backend

It stores the problem data in Sage matrices. It allows users to specify a base ring.

The class does not provide a solver method. It can be used as a base class for
other classes implementing solvers.
"""

# *****************************************************************************
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#                     2023 Alexandre Maranhão
#                     2023 Zhongling Xu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.matrix.constructor import Matrix as matrix
from copy import copy


cdef class MatrixBackend(GenericBackend):
    r"""
    MIP Backend that stores problem data in Sage matrices.

    There is no solver; calling the :meth:`solve` method raises a
    runtime error.

    EXAMPLES::

        sage: p = MixedIntegerLinearProgram(solver="Matrix")

    TESTS:

    :trac:`20332`::

        sage: p
        Mixed Integer Program (no objective, 0 variables, 0 constraints)
    """

    def __cinit__(self, maximization=True, base_ring=None, sparse=None, implementation=None, **kwds):
        r"""
        Cython constructor

        INPUT:

        - ``sparse`` -- passed on to :func:`matrix` constructor

        - ``implementation`` -- passed on to :func:`matrix` constructor

        - ``**kwds`` -- ignored

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
        """

        self.prob_name = ''
        self.obj_constant_term = 0
        self.is_maximize = 1

        self._init_base_ring(base_ring=base_ring)

        kwds = {}
        if implementation is not None:
            kwds['implementation'] = implementation
        if sparse is not None:
            kwds['sparse'] = sparse

        # Instead of vectors, we use 1-row matrices so that
        # an "implementation" keyword can be passed consistently
        self.objective_coefficients = matrix(self._base_ring, 1, [], **kwds)
        self.Matrix = matrix(self._base_ring, [], **kwds)
        self.row_lower_bound = matrix(self._base_ring, 1, [], **kwds)
        self.row_lower_bound_indicator = []
        self.row_upper_bound = matrix(self._base_ring, 1, [], **kwds)
        self.row_upper_bound_indicator = []
        self.col_lower_bound = matrix(self._base_ring, 1, [], **kwds)
        self.col_lower_bound_indicator = []
        self.col_upper_bound = matrix(self._base_ring, 1, [], **kwds)
        self.col_upper_bound_indicator = []

        self.row_name_var = []
        self.col_name_var = []

        self.is_continuous = []
        self.is_binary = []
        self.is_integer = []
        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

    def _init_base_ring(self, base_ring=None):
        r"""
        Handle a ``base_ring`` parameter passed to the constructor.

        This method can be overridden to implement a different defaulting behavior.

        This implementation use ``QQ`` as the default base ring.

        EXAMPLES::

            sage: from sage.numerical.backends.matrix_backend import MatrixBackend
            sage: MatrixBackend(base_ring=None).base_ring()  # indirect doctest
            Rational Field
        """
        if base_ring is None:
            from sage.rings.rational_field import QQ
            base_ring = QQ

        self._base_ring = base_ring

    cpdef base_ring(self):
        """
        The base ring

        TESTS::

            sage: from sage.numerical.backends.matrix_backend import MatrixBackend
            sage: MatrixBackend(base_ring=QQ).base_ring()
            Rational Field
            sage: MatrixBackend(base_ring=AA).base_ring()
            Algebraic Real Field
            sage: MatrixBackend(base_ring=RDF).base_ring()
            Real Double Field
            sage: MatrixBackend(base_ring=RealField(100)).base_ring()
            Real Field with 100 bits of precision
        """
        return self._base_ring

    cpdef __copy__(self):
        """
        Returns a copy of ``self``.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver="Matrix")
            sage: b = p.new_variable()
            sage: p.add_constraint(b[1] + b[2] <= 6)
            sage: p.set_objective(b[1] + b[2])
            sage: mat = copy(p.get_backend())
            sage: mat.ncols()
            2
        """
        cdef MatrixBackend mat = type(self)(base_ring=self.base_ring())
        mat.Matrix = copy(self.Matrix)
        mat.row_lower_bound = copy(self.row_lower_bound)
        mat.row_upper_bound = copy(self.row_upper_bound)
        mat.col_lower_bound = copy(self.col_lower_bound)
        mat.col_upper_bound = copy(self.col_upper_bound)
        mat.objective_coefficients = copy(self.objective_coefficients)
        mat.obj_constant_term = self.obj_constant_term
        return mat

    cpdef int add_variable(self, lower_bound=0, upper_bound=None,
                           binary=False, continuous=True, integer=False,
                           obj=None, name=None, coefficients=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both nonnegative and real.

        INPUT:

        - ``lower_bound`` -- the lower bound of the variable (default: 0)

        - ``upper_bound`` -- the upper bound of the variable (default: ``None``)

        - ``binary`` -- ``True`` if the variable is binary (default: ``False``)

        - ``continuous`` -- ``True`` if the variable is continuous (default: ``True``)

        - ``integer`` -- ``True`` if the variable is integer (default: ``False``)

        - ``obj`` -- (optional, default: 0) coefficient of this variable in the objective function

        - ``name`` -- an optional name for the newly added variable (default: ``None``)

        OUTPUT: The index of the newly created variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable()
            1
            sage: p.add_variable(lower_bound=-2.0)      # converted into base ring (QQ)
            2
            sage: p.add_variable(continuous=True)
            3
            sage: p.add_variable(name='x', obj=1.0)     # converted into base ring (QQ)
            4
            sage: p.col_name(3)
            'x_3'
            sage: p.col_name(4)
            'x'
            sage: p.objective_coefficient(4)
            1
            sage: p.objective_coefficient(4).parent()
            Rational Field

            sage: p = get_solver(solver="Matrix", base_ring=AA)
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable()
            1
            sage: p.add_variable(lower_bound=-2)        # converted into base ring (AA)
            2
            sage: p.add_variable(continuous=True)
            3
            sage: p.add_variable(name='x', obj=1)       # converted into base ring (AA)
            4
            sage: p.col_name(3)
            'x_3'
            sage: p.col_name(4)
            'x'
            sage: p.objective_coefficient(4)
            1
            sage: p.objective_coefficient(4).parent()
            Algebraic Real Field

            sage: p = get_solver(solver="Matrix", base_ring=RDF)
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable()
            1
            sage: p.add_variable(lower_bound=-2.0)      # converted into base ring (RDF)
            2
            sage: p.add_variable(continuous=True)
            3
            sage: p.add_variable(name='x', obj=1.0)     # converted into base ring (RDF)
            4
            sage: p.col_name(3)
            'x_3'
            sage: p.col_name(4)
            'x'
            sage: p.objective_coefficient(4)
            1.0
            sage: p.objective_coefficient(4).parent()
            Real Double Field

            sage: p = get_solver(solver="Matrix", base_ring=RealField(100))
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable()
            1
            sage: p.add_variable(lower_bound=-2.0)      # converted into base ring (RealField(100))
            2
            sage: p.add_variable(continuous=True)
            3
            sage: p.add_variable(name='x', obj=1.0)     # converted into base ring (RealField(100))
            4
            sage: p.col_name(3)
            'x_3'
            sage: p.col_name(4)
            'x'
            sage: p.objective_coefficient(4)
            1.0000000000000000000000000000
            sage: p.objective_coefficient(4).parent()
            Real Field with 100 bits of precision

        """
        if obj is None:
            obj = self._base_ring.zero()

        column = matrix(self._base_ring, self.nrows(), 1)
        if coefficients:
            for ind, coeff in coefficients:
                column[ind] = coeff

        self.Matrix = self.Matrix.augment(column)

        if lower_bound is None:
            self.col_lower_bound_indicator.append(False)
            self.col_lower_bound = self.col_lower_bound.augment(matrix([self._base_ring.zero()]))
        else:
            self.col_lower_bound_indicator.append(True)
            self.col_lower_bound = self.col_lower_bound.augment(matrix([lower_bound]))

        if upper_bound is None:
            self.col_upper_bound_indicator.append(False)
            self.col_upper_bound = self.col_upper_bound.augment(matrix([self._base_ring.zero()]))
        else:
            self.col_upper_bound_indicator.append(True)
            self.col_upper_bound = self.col_upper_bound.augment(matrix([upper_bound]))

        self.objective_coefficients = self.objective_coefficients.augment(matrix([obj]))
        self.col_name_var.append(name)
        self.is_binary.append(binary)
        self.is_continuous.append(continuous)
        self.is_integer.append(integer)

        return self.objective_coefficients.dimensions()[1] - 1

    cpdef set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable

        INPUT:

        - ``variable`` -- (integer) the variable's id

        - ``vtype`` -- (integer) one of
          *  `1`  Integer
          *  `0`  Binary
          *  `-1`  Continuous

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)
            True
            sage: p.set_variable_type(1,1)
            Traceback (most recent call last):
            ...
            IndexError: ...
        """
        if not 0 <= variable < self.ncols():
            raise IndexError("Variable index out of bounds")
        self.is_binary[variable] = (vtype == 0)
        self.is_integer[variable] = (vtype == 1)
        self.is_continuous[variable] = (vtype == -1)

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Matrix")
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
        Set or get the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` -- (integer) the variable's id

        - ``coeff`` -- its coefficient

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2
            sage: p = get_solver(solver="Matrix", base_ring=AA)
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2
            sage: p.objective_coefficient(0).parent()
            Algebraic Real Field
            sage: p = get_solver(solver="Matrix", base_ring=RDF)
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0.0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2.0
            sage: p.objective_coefficient(0).parent()
            Real Double Field
            sage: p = get_solver(solver="Matrix", base_ring=RealField(100))
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0.00000000000000000000000000000
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2.0000000000000000000000000000
            sage: p.objective_coefficient(0).parent()
            Real Field with 100 bits of precision

        """

        if coeff is not None:
            self.objective_coefficients[0, variable] = coeff
        else:
            return self.objective_coefficients[0, variable]

    cpdef set_objective(self, list coeff, d=0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose `i`-th element is the
          coefficient of the `i`-th variable in the objective function.

        - ``d`` -- the constant term in the linear function (set to `0` by default)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1, 1, 2, 1, 3]

            sage: p = get_solver(solver="Matrix", base_ring=AA)
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1, 1, 2, 1, 3]
            sage: [p.objective_coefficient(x).parent() for x in range(5)]
            [Algebraic Real Field,
             Algebraic Real Field,
             Algebraic Real Field,
             Algebraic Real Field,
             Algebraic Real Field]

            sage: p = get_solver(solver="Matrix", base_ring=RDF)
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1.0, 1.0, 2.0, 1.0, 3.0]
            sage: [p.objective_coefficient(x).parent() for x in range(5)]
            [Real Double Field,
             Real Double Field,
             Real Double Field,
             Real Double Field,
             Real Double Field]

            sage: p = get_solver(solver="Matrix", base_ring=RealField(100))
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1.0000000000000000000000000000,
             1.0000000000000000000000000000,
             2.0000000000000000000000000000,
             1.0000000000000000000000000000,
             3.0000000000000000000000000000]
            sage: [p.objective_coefficient(x).parent() for x in range(5)]
            [Real Field with 100 bits of precision,
             Real Field with 100 bits of precision,
             Real Field with 100 bits of precision,
             Real Field with 100 bits of precision,
             Real Field with 100 bits of precision]

        """
        for i in range(len(coeff)):
            self.objective_coefficients[0, i] = coeff[i]
        self.obj_constant_term = d

    cpdef set_verbosity(self, int level):
        """
        Does not apply
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
            sage: p = get_solver(solver="Matrix")
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5

            sage: p = get_solver(solver="Matrix", base_ring=AA)
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5

            sage: p = get_solver(solver="Matrix", base_ring=RDF)
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5

            sage: p = get_solver(solver="Matrix", base_ring=RealField(100))
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5
        """
        self.add_variable(coefficients=zip(indices, coeffs))

    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` -- an iterable with ``(c,v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (real value)

        - ``lower_bound`` -- a lower bound, either a real value or ``None``

        - ``upper_bound`` -- an upper bound, either a real value or ``None``

        - ``name`` -- an optional name for this row (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2.0, 2.0)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)
            (2, 2)
            sage: p.add_linear_constraint(((1,7), (3,10)), 2.0, 2.0)
            sage: p.row(1)
            ([1, 3], [7, 10])
            sage: p.add_linear_constraint(zip(range(5), range(5)), 1.0, 1.0, name='foo')
            sage: p.row_name(-1)
            'foo'
        """
        coefficients = list(coefficients)
        for c in coefficients:
            while c[0] > self.ncols() - 1:
                self.add_variable()

        self.Matrix = self.Matrix.stack(matrix(self.base_ring(), 1, self.ncols()))

        for c in coefficients:
            self.Matrix[-1, c[0]] = c[1]

        if lower_bound is None:
            self.row_lower_bound_indicator.append(False)
            self.row_lower_bound = self.row_lower_bound.augment(matrix([self._base_ring.zero()]))
        else:
            self.row_lower_bound_indicator.append(True)
            self.row_lower_bound = self.row_lower_bound.augment(matrix([lower_bound]))

        if upper_bound is None:
            self.row_upper_bound_indicator.append(False)
            self.row_upper_bound = self.row_upper_bound.augment(matrix([self._base_ring.zero()]))
        else:
            self.row_upper_bound_indicator.append(True)
            self.row_upper_bound = self.row_upper_bound.augment(matrix([upper_bound]))

        self.row_name_var.append(name)

    cpdef int ncols(self):
        """
        Return the number of columns/variables.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.ncols()
            0
            sage: p.add_variables(2)
            1
            sage: p.ncols()
            2
        """
        return self.objective_coefficients.dimensions()[1]

    cpdef int nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.nrows()
            0
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraints(2, 2.0, None)
            sage: p.nrows()
            2
        """
        return self.row_upper_bound.dimensions()[1]

    cpdef bint is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        return self.is_maximize == 1

    cpdef problem_name(self, name=None):
        """
        Return or define the problem's name

        INPUT:

        - ``name`` (``str``) -- the problem's name. When set to
          ``None`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.problem_name()
            ''
            sage: p.problem_name("There once was a french fry")
            sage: print(p.problem_name())
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
        :meth:`add_linear_constraint` method.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)
            (2, 2)

            sage: p = get_solver(solver="Matrix", base_ring=AA)
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: [i.parent() for i in p.row(0)[1]]
            [Algebraic Real Field,
             Algebraic Real Field,
             Algebraic Real Field,
             Algebraic Real Field]
            sage: p.row_bounds(0)
            (2, 2)
            sage: [i.parent() for i in p.row_bounds(0)]
            [Algebraic Real Field, Algebraic Real Field]

            sage: p = get_solver(solver="Matrix", base_ring=RDF)
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: [i.parent() for i in p.row(0)[1]]
            [Real Double Field, Real Double Field, Real Double Field, Real Double Field]
            sage: p.row_bounds(0)
            (2.0, 2.0)
            sage: [i.parent() for i in p.row_bounds(0)]
            [Real Double Field, Real Double Field]

            sage: p = get_solver(solver="Matrix", base_ring=RealField(100))
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4],
            [1.0000000000000000000000000000,
             2.0000000000000000000000000000,
             3.0000000000000000000000000000,
             4.0000000000000000000000000000])
            sage: [i.parent() for i in p.row(0)[1]]
            [Real Field with 100 bits of precision,
             Real Field with 100 bits of precision,
             Real Field with 100 bits of precision,
             Real Field with 100 bits of precision]
            sage: p.row_bounds(0)
            (2.0000000000000000000000000000, 2.0000000000000000000000000000)
            sage: [i.parent() for i in p.row_bounds(0)]
            [Real Field with 100 bits of precision, Real Field with 100 bits of precision]

        """
        coeff = []
        idx = []
        index = 0

        for num in self.Matrix[i]:
            if num != 0:
                idx.append(index)
                coeff.append(num)
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
            sage: p = get_solver(solver="Matrix")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)
            (2, 2)
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, None)
            sage: p.row_bounds(1)
            (2, None)
        """
        if self.row_lower_bound_indicator[index]:
            lower = self.row_lower_bound[0, index]
        else:
            lower = None

        if self.row_upper_bound_indicator[index]:
            upper = self.row_upper_bound[0, index]
        else:
            upper = None

        return (lower, upper)

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
            sage: p = get_solver(solver="Matrix")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0, 5)

            sage: p = get_solver(solver="Matrix", base_ring=AA)
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0, None)
            sage: p.col_bounds(0)[0].parent()
            Algebraic Real Field
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0, 5)
            sage: [i.parent() for i in p.col_bounds(0)]
            [Algebraic Real Field, Algebraic Real Field]

            sage: p = get_solver(solver="Matrix", base_ring=RDF)
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.col_bounds(0)[0].parent()
            Real Double Field
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0.0, 5.0)
            sage: [i.parent() for i in p.col_bounds(0)]
            [Real Double Field, Real Double Field]

            sage: p = get_solver(solver="Matrix", base_ring=RealField(100))
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.00000000000000000000000000000, None)
            sage: p.col_bounds(0)[0].parent()
            Real Field with 100 bits of precision
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0.00000000000000000000000000000, 5.0000000000000000000000000000)
            sage: [i.parent() for i in p.col_bounds(0)]
            [Real Field with 100 bits of precision, Real Field with 100 bits of precision]

        """
        if self.col_lower_bound_indicator[index]:
            lower = self.col_lower_bound[0, index]
        else:
            lower = None

        if self.col_upper_bound_indicator[index]:
            upper = self.col_upper_bound[0, index]
        else:
            upper = None

        return (lower, upper)

    cpdef bint is_variable_binary(self, int index):
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0, 0)
            sage: p.is_variable_binary(0)
            True
            sage: p.is_variable_binary(1)
            False
        """
        if not 0 <= index < self.ncols():
            # This is how the other backends behave, and this method is
            # unable to raise a python exception as currently defined.
            return False

        return self.is_binary[index]

    cpdef bint is_variable_integer(self, int index):
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0, -1)
            sage: p.set_variable_type(0, 1)
            sage: p.is_variable_integer(0)
            True
            sage: p.is_variable_integer(1)
            False
        """
        if not 0 <= index < self.ncols():
            # This is how the other backends behave, and this method is
            # unable to raise a python exception as currently defined.
            return False

        return self.is_integer[index]

    cpdef bint is_variable_continuous(self, int index):
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_continuous(0)
            True
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_continuous(0)
            False
            sage: p.is_variable_binary(0)
            False
            sage: p.is_variable_integer(0)
            True
            sage: p.is_variable_continuous(1)
            False
        """
        if not 0 <= index < self.ncols():
            # This is how the other backends behave, and this method is
            # unable to raise a python exception as currently defined.
            return False

        return self.is_continuous[index]

    cpdef row_name(self, int index):
        """
        Return the ``index``-th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.add_linear_constraints(1, 2, None, names=["Empty constraint 1"])
            sage: p.row_name(0)
            'Empty constraint 1'
        """
        if self.row_name_var[index] is not None:
            return self.row_name_var[index]
        return "constraint_" + repr(index)

    cpdef col_name(self, int index):
        """
        Return the ``index``-th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        - ``name`` (``char *``) -- its name. When set to ``NULL``
          (default), the method returns the current name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
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
        Return or define the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0, 5)
        """
        if value is None:
            self.col_upper_bound_indicator[index] = False
        elif value is not False:
            self.col_upper_bound_indicator[index] = True
            self.col_upper_bound[0, index] = value
        else:
            return self.col_upper_bound[0, index]

    cpdef variable_lower_bound(self, int index, value=False):
        """
        Return or define the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver="Matrix")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0, None)
            sage: p.variable_lower_bound(0, 5)
            sage: p.col_bounds(0)
            (5, None)

            sage: p = get_solver(solver="Matrix", base_ring=AA)
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0, None)
            sage: p.col_bounds(0)[0].parent()
            Algebraic Real Field
            sage: p.variable_lower_bound(0, 5)
            sage: p.col_bounds(0)
            (5, None)
            sage: p.col_bounds(0)[0].parent()
            Algebraic Real Field

            sage: p = get_solver(solver="Matrix", base_ring=RDF)
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.col_bounds(0)[0].parent()
            Real Double Field
            sage: p.variable_lower_bound(0, 5)
            sage: p.col_bounds(0)
            (5.0, None)
            sage: p.col_bounds(0)[0].parent()
            Real Double Field

            sage: p = get_solver(solver="Matrix", base_ring=RealField(100))
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.00000000000000000000000000000, None)
            sage: p.col_bounds(0)[0].parent()
            Real Field with 100 bits of precision
            sage: p.variable_lower_bound(0, 5)
            sage: p.col_bounds(0)
            (5.0000000000000000000000000000, None)
            sage: p.col_bounds(0)[0].parent()
            Real Field with 100 bits of precision
        """
        if value is None:
            self.col_lower_bound_indicator[index] = False
        elif value is not False:
            self.col_lower_bound_indicator[index] = True
            self.col_lower_bound[0, index] = value
        else:
            return self.col_lower_bound[0, index]
