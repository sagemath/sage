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

from sage.matrix.constructor import Matrix
from sage.matrix.constructor import Vector
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
    def __cinit__(self, maximization = True, base_ring=None):
        """
        Cython constructor

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Matrix")                 

"""
        
        #Sage Matrix and Vector instead of Python lists
        self.objective_function = Vector(QQ, [])
        self.G_matrix = Matrix(QQ, [])

        self.prob_name = ''
        self.obj_constant_term = 0
        self.is_maximize = 1

        self.row_lower_bound = Vector(QQ, [])
        self.row_upper_bound = Vector(QQ, [])
        self.col_lower_bound = Vector(QQ, [])
        self.col_upper_bound = Vector(QQ, [])

        self.row_name_var = []
        self.col_name_var = []

        self.is_continuous = []
        self.is_binary = []
        self.is_integer = []
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
    The base ring
    
    TESTS::
    
            sage: from sage.numerical.backends.matrix_backend import MatrixBackend
            sage: MatrixBackend(base_ring=QQ).base_ring()
            Rational Field
        """
        return self._base_ring


    cpdef int add_variable(self, lower_bound=0.0, upper_bound=None, binary=False, continuous=True, integer=False, obj=None, name=None) except -1:
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
            1.00000000000000

        """
        if obj is None:
            obj = 0.0
        self.G_matrix.augment(Vector([0 for i in range(self.nrows())]))
        self.col_lower_bound.augment(lower_bound)
        self.col_upper_bound.augment(upper_bound)
        self.objective_function.augment(obj)
        self.col_name_var.augment(name)
        self.is_binary.augment(binary)
        self.is_continuous.augment(continuous)
        self.is_integer.augment(integer)
        return len(self.objective_function) - 1

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
            0.0
            sage: p.objective_coefficient(0,2)                      
            sage: p.objective_coefficient(0)                        
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

            sage: from sage.numerical.backends.generic_backend import get_solver    
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variables(5)                                
            4
            sage: p.set_objective([1, 1, 2, 1, 3])                  
            sage: [p.objective_coefficient(x) for x in range(5)]    
            [1, 1, 2, 1, 3]
        """
        for i in range(len(coeff)):
            self.objective_function[i] = coeff[i]
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
        """
        column = []
        for _ in indices:
            column.augment(Vector([0.0]))

        for idx, ind in enumerate(indices):
            column[ind] = coeffs[idx]

        self.G_matrix.augment(column)

        self.col_lower_bound.augment(None)
        self.col_upper_bound.augment(None)
        self.objective_function.augment(0)
        self.col_name_var.augment(None)

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
            (2.00000000000000, 2.00000000000000)
            sage: p.add_linear_constraint(zip(range(5), range(5)), 1.0, 1.0, name='foo')    
            sage: p.row_name(-1)                                    
            'foo'
        """
        coefficients = list(coefficients)
        for c in coefficients:
            while c[0] > len(self.G_matrix) - 1:
                self.add_variable()
        for i in range(len(self.G_matrix)):
            self.G_matrix[i].augment(Vector([0.0]))
        for c in coefficients:
            self.G_matrix[c[0]][-1] = c[1]

        self.row_lower_bound.augment(lower_bound)
        self.row_upper_bound.augment(upper_bound)
        self.row_name_var.augment(name)



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

        return len(self.objective_function)

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
        return len(self.row_upper_bound)


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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variables(5)                                
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)  
            sage: p.row(0)                                          
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)                                   
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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variable()                                  
            0
            sage: p.col_bounds(0)                                   
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                      
            sage: p.col_bounds(0)                                   
            (0.0, 5)
        """
        return (self.col_lower_bound[index], self.col_upper_bound[index])

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
            IndexError: ...
            ...

        """
        if index >= len(self.is_binary):
            raise IndexError("Variable index out of bounds")
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
            IndexError: ...
            ...
        """
        if index >= len(self.is_integer):
            raise IndexError("Variable index out of bounds")
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
            IndexError: ...
            ...

        """
        if index >= len(self.is_continuous):
            raise IndexError("Variable index out of bounds")
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
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                      
            sage: p.col_bounds(0)                                   
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
            sage: p = get_solver(solver = "Matrix")                 
            sage: p.add_variable()                                  
            0
            sage: p.col_bounds(0)                                   
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)                      
            sage: p.col_bounds(0)                                   
            (5, None)
        """
        if value is not False:
            self.col_lower_bound[index] = value
        else:
            return self.col_lower_bound[index]

cdef class NumpyMatrixBackend(GenericBackend):

