r"""
Matrix Backend

It stores the problem data in Sage matrices. It allows users to specify a base ring.

The class does not provide a solver method. It can be used as a base class for
other classes implementing solvers.

"""

#*****************************************************************************
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import Matrix as matrix
from .generic_backend cimport GenericBackend
import numpy

cdef class MatrixBackend(GenericBackend):
    """
    Cython constructor

    There is no support for integer variables.

    EXAMPLES::

        sage: p = MixedIntegerLinearProgram(solver="Matrix")

    TESTS:

    :trac:`20332`::

        sage: p                                                     
        Mixed Integer Program (no objective, 0 variables, 0 constraints)
    """
    def __cinit__(self, maximization = True, base_ring=None, numpy_implementation = False, **kwds):
        """
        Cython constructor

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Matrix")                 

        """

        #Sage Matrix and Vector instead of Python lists

        self.prob_name = ''
        self.obj_constant_term = 0
        self.is_maximize = 1

        #if base_ring is None:
        #    from sage.rings.rational_field import QQ
        #    base_ring = QQ
        self._init_base_ring(base_ring=base_ring)

        if numpy_implementation:
            kwds = dict(implementation = "numpy")
        else:
            kwds = {}

        self.objective_function = matrix(self._base_ring, 1, [], **kwds)
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

        #self._base_ring = base_ring
    
    def _init_base_ring(self, base_ring=None):
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


    cpdef int add_variable(self, lower_bound=0, upper_bound=None, binary=False, continuous=True, integer=False, obj=None, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.ncols()                                         
            0
            sage: p.add_variable()                                  
            0
            sage: p.ncols()                                         
            1
            sage: p.add_variable()                                  
            1
            sage: p.add_variable(lower_bound=-2.0)                  
            2
            sage: p.add_variable(continuous=True)                   
            3
            sage: p.add_variable(name='x',obj=1.0)                  
            4
            sage: p.col_name(3)                                     
            'x_3'
            sage: p.col_name(4)                                     
            'x'
            sage: p.objective_coefficient(4)                        
            1
            sage: p.objective_coefficient(4).parent()
            Rational Field
            sage: p = get_solver(solver = "Matrix", base_ring = AA)                 
            sage: p.ncols()                                         
            0
            sage: p.add_variable()                                  
            0
            sage: p.ncols()                                         
            1
            sage: p.add_variable()                                  
            1
            sage: p.add_variable(lower_bound=-2)                  
            2
            sage: p.add_variable(continuous=True)                   
            3
            sage: p.add_variable(name='x',obj=1)                  
            4
            sage: p.col_name(3)                                     
            'x_3'
            sage: p.col_name(4)                                     
            'x'
            sage: p.objective_coefficient(4)                        
            1
            sage: p.objective_coefficient(4).parent()
            Algebraic Real Field
            sage: p = get_solver(solver = "Matrix", base_ring = RDF)                 
            sage: p.ncols()                                         
            0
            sage: p.add_variable()                                  
            0
            sage: p.ncols()                                         
            1
            sage: p.add_variable()                                  
            1
            sage: p.add_variable(lower_bound=-2.0)                  
            2
            sage: p.add_variable(continuous=True)                   
            3
            sage: p.add_variable(name='x',obj=1.0)                  
            4
            sage: p.col_name(3)                                     
            'x_3'
            sage: p.col_name(4)                                     
            'x'
            sage: p.objective_coefficient(4)                        
            1.0
            sage: p.objective_coefficient(4).parent()
            Real Double Field
            sage: p = get_solver(solver = "Matrix", base_ring = RealField(100))
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable()
            1
            sage: p.add_variable(lower_bound=-2.0)
            2
            sage: p.add_variable(continuous=True)
            3
            sage: p.add_variable(name='x',obj=1.0)
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

        if self.nrows() == 0:
            pass
        else:
            self.Matrix = self.Matrix.augment(matrix([self._base_ring.zero() for i in range(self.nrows())]))
        
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

        self.objective_function = self.objective_function.augment(matrix([obj]))
        self.col_name_var.append(name)
        self.is_binary.append(binary)
        self.is_continuous.append(continuous)
        self.is_integer.append(integer)

        return self.objective_function.dimensions()[1] - 1

    cpdef set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``vtype`` (integer) :

            *  1  Integer
            *  0  Binary
            *  -1  Continuous

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Matrix")   
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
        if variable >= len(self.is_continuous):
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
        Set or get the coefficient of a variable in the objective
        function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variable()                                  
            0
            sage: p.objective_coefficient(0)                        
            0
            sage: p.objective_coefficient(0,2)                      
            sage: p.objective_coefficient(0)                        
            2
            sage: p = get_solver(solver = "Matrix", base_ring=AA)
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2
            sage: p.objective_coefficient(0).parent()
            Algebraic Real Field
            sage: p = get_solver(solver = "Matrix", base_ring=RDF)
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0.0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2.0
            sage: p.objective_coefficient(0).parent()
            Real Double Field
            sage: p = get_solver(solver = "Matrix", base_ring=RealField(100))
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
            self.objective_function[0, variable] = coeff
        else:
            return self.objective_function[0, variable]

    cpdef set_objective(self, list coeff, d = 0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        - ``d`` (double) -- the constant term in the linear function (set to `0` by default)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver    
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variables(5)                                
            4
            sage: p.set_objective([1, 1, 2, 1, 3])                  
            sage: [p.objective_coefficient(x) for x in range(5)]    
            [1, 1, 2, 1, 3]
            sage: p = get_solver(solver = "Matrix", base_ring = AA)
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
            sage: p = get_solver(solver = "Matrix", base_ring = RDF)
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
            sage: p = get_solver(solver = "Matrix", base_ring = RealField(100))
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
            self.objective_function[0, i] = coeff[i]
        obj_constant_term = d

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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.ncols()                                         
            0
            sage: p.nrows()                                         
            0
            sage: p.add_linear_constraints(5, 0, None)              
            sage: p.add_col(range(5), range(5))                     
            sage: p.nrows()                                         
            5
            sage: p = get_solver(solver = "Matrix", base_ring=AA)
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5
            sage: p = get_solver(solver = "Matrix", base_ring=RDF)
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5
            sage: p = get_solver(solver = "Matrix", base_ring=RealField(100))
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5
        """
        column = []
        for _ in indices:
            column.append(self._base_ring.zero())

        for idx, ind in enumerate(indices):
            column[ind] = coeffs[idx]

        if self.ncols() == 0:
            self.Matrix = matrix(len(column), column)
        else:
            self.Matrix = self.Matrix.augment(matrix(len(column), column))

        self.col_lower_bound_indicator.append(None)
        if self.col_lower_bound.dimensions()[1] == 0:
            self.col_lower_bound = matrix(1, [self._base_ring.zero()])
        else:
            self.col_lower_bound = self.col_lower_bound.augment(matrix([self._base_ring.zero()]))

        self.col_upper_bound_indicator.append(None)
        if self.col_upper_bound.dimensions()[1] == 0:
            self.col_upper_bound = matrix(1, [self._base_ring.zero()])
        else:
            self.col_upper_bound = self.col_upper_bound.augment(matrix([self._base_ring.zero()]))
            
        self.objective_function = self.objective_function.augment(matrix([self._base_ring.zero()]))
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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variables(5)                                
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2.0, 2.0)                
            sage: p.row(0)                                          
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)                                   
            (2, 2)
            sage: p.add_linear_constraint(zip(range(5), range(5)), 1.0, 1.0, name='foo')    
            sage: p.row_name(-1)                                    
            'foo'
        """
        coefficients = list(coefficients)
        for c in coefficients:
            while c[0] > self.ncols() - 1:
                self.add_variable()
                     
        if self.Matrix.dimensions()[0] == 0:
            self.Matrix = matrix(1, [self._base_ring.zero() for i in range(len(coefficients))])
        else:
            self.Matrix = self.Matrix.stack(matrix([self._base_ring.zero() for i in range(self.Matrix.dimensions()[1])]))

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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.ncols()                                         
            0
            sage: p.add_variables(2)                                
            1
            sage: p.ncols()                                         
            2
        """

        return self.objective_function.dimensions()[1]

    cpdef int nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Matrix")                 
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
            sage: p = get_solver(solver = "Matrix")                 
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
        Return or define the problem's name

        INPUT:

        - ``name`` (``str``) -- the problem's name. When set to
          ``None`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Matrix")                 
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
        ``add_linear_constraint`` method.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variables(5)                                
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)  
            sage: p.row(0)                                          
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)                                   
            (2, 2)
            sage: p = get_solver(solver = "Matrix", base_ring = AA)
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
            sage: p = get_solver(solver = "Matrix", base_ring = RDF)
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: [i.parent() for i in p.row(0)[1]]
            [Real Double Field, Real Double Field, Real Double Field, Real Double Field]
            sage:  p.row_bounds(0)
            (2.0, 2.0)
            sage: [i.parent() for i in p.row_bounds(0)]
            [Real Double Field, Real Double Field]
            sage: p = get_solver(solver = "Matrix", base_ring = RealField(100))
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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variables(5)                                
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)  
            sage: p.row(0)                                          
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)                                   
            (2, 2)
        """
        if self.row_lower_bound_indicator[index] == True:
            lower = self.row_lower_bound[0, index]
        else:
            lower = None
        
        if self.row_upper_bound_indicator[index] == True:
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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variable()                                  
            0
            sage: p.col_bounds(0)                                   
            (0, None)
            sage: p.variable_upper_bound(0, 5)                      
            sage: p.col_bounds(0)                                   
            (0, 5)
            sage: p = get_solver(solver = "Matrix", base_ring = AA)
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
            sage: p = get_solver(solver = "Matrix", base_ring = RDF)
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
            sage: p = get_solver(solver = "Matrix", base_ring = RealField(100))
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
        if self.col_lower_bound_indicator[index] == True:
            lower = self.col_lower_bound[0, index]
        else:
            lower = None
        
        if self.col_upper_bound_indicator[index] == True:
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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.ncols()                                         
            0
            sage: p.add_variable()                                  
            0
            sage: p.set_variable_type(0,0)
            sage: p.is_variable_binary(0)                           
            True
            sage: p.is_variable_binary(1)    
            False
        """
        if index < 0 or index >= len(self.is_binary):
            # This is how the other backends behave, and this method is
            # unable to raise a python exception as currently defined.
            return False

        #if index >= len(self.is_binary):
        #    raise IndexError("Variable index out of bounds")
            
        return self.is_binary[index]

    cpdef bint is_variable_integer(self, int index):
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.ncols()                                         
            0
            sage: p.add_variable()                                  
            0
            sage: p.set_variable_type(0,-1)                         
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)                          
            True
            sage: p.is_variable_integer(1)        
            False
        """

        if index < 0 or index >= len(self.is_integer):
            # This is how the other backends behave, and this method is
            # unable to raise a python exception as currently defined.
            return False

        #if index >= len(self.is_integer):
        #    raise IndexError("Variable index out of bounds")
        
        
        return self.is_integer[index]

    cpdef bint is_variable_continuous(self, int index):
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Matrix")                 
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

        if index < 0 or index >= len(self.is_continuous):
            # This is how the other backends behave, and this method is
            # unable to raise a python exception as currently defined.
            return False

        #if index >= len(self.is_continuous):
        #    raise IndexError("Variable index out of bounds")

        return self.is_continuous[index]

    cpdef row_name(self, int index):
        """
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_linear_constraints(1, 2, None, names=["Empty constraint 1"])  
            sage: p.row_name(0)                                     
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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variable(name="I am a variable")            
            0
            sage: p.col_name(0)                                     
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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variable()                                  
            0
            sage: p.col_bounds(0)                                   
            (0, None)
            sage: p.variable_upper_bound(0, 5)                      
            sage: p.col_bounds(0)                                   
            (0, 5)
        """
        if value is not False:
            self.col_upper_bound_indicator[index] = True
            self.col_upper_bound[0, index] = value
        else:
            return self.col_upper_bound[0, index]

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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variable()                                  
            0
            sage: p.col_bounds(0)                                   
            (0, None)
            sage: p.variable_lower_bound(0, 5)                      
            sage: p.col_bounds(0)                                   
            (5, None)
            sage: p = get_solver(solver = "Matrix", base_ring = AA)
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
            sage: p = get_solver(solver = "Matrix", base_ring = RDF)
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
            sage: p = get_solver(solver = "Matrix", base_ring = RealField(100))
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
        if value is not False:
            self.col_lower_bound_indicator[index] = True
            self.col_lower_bound[0, index] = value
        else:
            return self.col_lower_bound[0, index]