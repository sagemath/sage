# sage.doctest: optional - highspy
"""
HiGHS Backend

AUTHORS:

- SageMath Developers (2025): initial implementation

This backend uses the HiGHS optimization solver, which supports Linear Programming (LP),
Quadratic Programming (QP), and Mixed Integer Programming (MIP).

HiGHS is available under the MIT License.
"""

# ****************************************************************************
#       Copyright (C) 2025 SageMath Developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.numerical.mip import MIPSolverException
from copy import copy


cdef class HiGHSBackend(GenericBackend):
    """
    MIP Backend that uses the HiGHS solver.

    HiGHS is a high-performance solver for large-scale LP, QP, and MIP.
    """

    def __cinit__(self, maximization=True):
        """
        Constructor.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='HiGHS')
        """
        try:
            import highspy
        except ImportError:
            raise ImportError("HiGHS is not available. Please install the 'highspy' package.")
        
        self.highs_model = highspy.Highs()
        # Suppress HiGHS output messages
        self.highs_model.setOptionValue("log_to_console", False)
        self.prob_name = ""
        self.col_name_var = {}
        self.row_name_var = {}
        self.numcols = 0
        self.numrows = 0
        self.obj_constant_term = 0.0
        
        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)
    
    cpdef int add_variable(self, lower_bound=0.0, upper_bound=None, binary=False, 
                          continuous=False, integer=False, obj=0.0, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive, real and the coefficient in the
        objective function is 0.0.

        INPUT:

        - ``lower_bound`` -- the lower bound of the variable (default: 0)

        - ``upper_bound`` -- the upper bound of the variable (default: ``None``)

        - ``binary`` -- ``True`` if the variable is binary (default: ``False``)

        - ``continuous`` -- ``True`` if the variable is continuous (default: ``True``)

        - ``integer`` -- ``True`` if the variable is integral (default: ``False``)

        - ``obj`` -- (optional) coefficient of this variable in the objective function (default: 0.0)

        - ``name`` -- an optional name for the newly added variable (default: ``None``)

        OUTPUT: the index of the newly created variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "HiGHS")
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
            sage: p.add_variable(name='x', obj=1.0)
            3
            sage: p.col_name(3)                    
            'x'
            sage: p.objective_coefficient(3)       
            1.0
        """
        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")
        
        # Determine variable type
        import highspy
        
        if binary:
            var_type = highspy.HighsVarType.kInteger
            if lower_bound is None or lower_bound < 0:
                lower_bound = 0.0
            if upper_bound is None or upper_bound > 1:
                upper_bound = 1.0
        elif integer:
            var_type = highspy.HighsVarType.kInteger
        else:
            var_type = highspy.HighsVarType.kContinuous
        
        # Handle bounds
        if lower_bound is None:
            lower_bound = -highspy.kHighsInf
        if upper_bound is None:
            upper_bound = highspy.kHighsInf
        
        # Add column to model
        col_idx = self.numcols
        self.highs_model.addVar(float(lower_bound), float(upper_bound))
        
        # Set objective coefficient
        if obj is not None and obj != 0.0:
            self.highs_model.changeColCost(col_idx, float(obj))
        
        # Set variable type
        if var_type == highspy.HighsVarType.kInteger:
            self.highs_model.changeColIntegrality(col_idx, highspy.HighsVarType.kInteger)
        
        # Set name if provided
        if name is not None:
            self.col_name_var[col_idx] = name
        
        self.numcols += 1
        return col_idx
    
    cpdef int add_variables(self, int number, lower_bound=0.0, upper_bound=None,
                           binary=False, continuous=False, integer=False, obj=0.0, 
                           names=None) except -1:
        """
        Add ``number`` new variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive, real and their coefficient in
        the objective function is 0.0.

        INPUT:

        - ``n`` -- the number of new variables (must be > 0)

        - ``lower_bound`` -- the lower bound of the variable (default: 0)

        - ``upper_bound`` -- the upper bound of the variable (default: ``None``)

        - ``binary`` -- ``True`` if the variable is binary (default: ``False``)

        - ``continuous`` -- ``True`` if the variable is binary (default: ``True``)

        - ``integer`` -- ``True`` if the variable is binary (default: ``False``)

        - ``obj`` -- coefficient of all variables in the objective function (default: 0.0)

        - ``names`` -- list of names (default: ``None``)

        OUTPUT: the index of the variable created last

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "HiGHS")
            sage: p.ncols()                       
            0
            sage: p.add_variables(5)              
            4
            sage: p.ncols()                       
            5
            sage: p.add_variables(2, lower_bound=-2.0, integer=True, obj=42.0, names=['a','b'])
            6

        TESTS:

        Check that arguments are used::

            sage: p.col_bounds(5)  # tol 1e-8   
            (-2.0, None)
            sage: p.is_variable_integer(5)       
            True
            sage: p.col_name(5)                  
            'a'
            sage: p.objective_coefficient(5)     
            42.0
        """
        cdef int i
        cdef int start_idx = self.numcols
        
        for i in range(number):
            name = None
            if names is not None and i < len(names):
                name = names[i]
            self.add_variable(lower_bound=lower_bound, upper_bound=upper_bound,
                            binary=binary, continuous=continuous, integer=integer,
                            obj=obj, name=name)
        
        return self.numcols - 1
    
    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` -- integer:

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "HiGHS")
            sage: p.is_maximization()             
            True
            sage: p.set_sense(-1)                 
            sage: p.is_maximization()             
            False
        """
        import highspy
        if sense == 1:
            self.highs_model.changeObjectiveSense(highspy.ObjSense.kMaximize)
        else:
            self.highs_model.changeObjectiveSense(highspy.ObjSense.kMinimize)
    
    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective function.

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``coeff`` -- double; its coefficient or ``None`` for
          reading (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "HiGHS")
            sage: p.add_variable()                
            0
            sage: p.objective_coefficient(0)      
            0.0
            sage: p.objective_coefficient(0,2)    
            sage: p.objective_coefficient(0)      
            2.0
        """
        if coeff is None:
            lp = self.highs_model.getLp()
            return lp.col_cost_[variable]
        else:
            self.highs_model.changeColCost(variable, float(coeff))
    
    cpdef problem_name(self, name=None):
        """
        Return or define the problem's name.

        INPUT:

        - ``name`` -- string; the problem's name. When set to
          ``None`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "HiGHS")
            sage: p.problem_name("There once was a french fry")
            sage: print(p.problem_name())         
            There once was a french fry
        """
        if name is None:
            return self.prob_name
        else:
            self.prob_name = str(name)
    
    cpdef set_objective(self, list coeff, d=0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- list of real values, whose i-th element is the
          coefficient of the i-th variable in the objective function
        - ``d`` -- constant term in objective function (default: 0.0)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(5)            
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
        """
        cdef int i
        for i in range(len(coeff)):
            if i < self.numcols:
                self.objective_coefficient(i, coeff[i])
        
        self.obj_constant_term = d
    
    cpdef set_verbosity(self, int level):
        """
        Set the log (verbosity) level.

        INPUT:

        - ``level`` -- integer; from 0 (no verbosity) to 3

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.set_verbosity(0)            
        """
        # HiGHS uses different verbosity levels
        # 0 = no output, 1 = minimal, 2 = detailed, 3 = verbose
        self.highs_model.setOptionValue("log_to_console", level > 0)
        if level > 0:
            self.highs_model.setOptionValue("log_dev_level", min(level - 1, 2))
    
    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` -- an iterable of pairs ``(i, v)`` where ``i`` is a
          variable index and ``v`` is a value
        - ``lower_bound`` -- a lower bound, either a real value or ``None``
        - ``upper_bound`` -- an upper bound, either a real value or ``None``
        - ``name`` -- optional name for this constraint

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(5)            
            4
            sage: p.add_linear_constraint([(0, 1), (1, 1)], None, 2.0)
        """
        import highspy
        
        if lower_bound is None:
            lower_bound = -highspy.kHighsInf
        if upper_bound is None:
            upper_bound = highspy.kHighsInf
        
        # Build constraint
        indices = []
        values = []
        for (i, v) in coefficients:
            if i < 0 or i >= self.ncols():
                raise ValueError(f"invalid variable index {i}")
            indices.append(int(i))
            values.append(float(v))
        
        # Add row to model using numpy arrays
        import numpy as np
        row_idx = self.numrows
        num_nz = len(indices)
        if num_nz > 0:
            index_array = np.array(indices, dtype=np.int32)
            value_array = np.array(values, dtype=np.double)
            self.highs_model.addRow(float(lower_bound), float(upper_bound), num_nz, index_array, value_array)
        else:
            # Empty constraint
            self.highs_model.addRow(float(lower_bound), float(upper_bound), 0, 0, 0)
        
        if name is not None:
            self.row_name_var[name] = row_idx
        
        self.numrows += 1
    
    cpdef int solve(self) except -1:
        """
        Solve the problem.

        Sage uses HiGHS's implementation of the branch-and-cut
        algorithm to solve mixed-integer linear programs. HiGHS
        automatically selects the most appropriate algorithm based
        on the problem type.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution cannot be computed for any reason (none
            exists, or the solver was not able to find it, etc...)

        EXAMPLES::

            sage: lp = MixedIntegerLinearProgram(solver = 'HiGHS', maximization = False)
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

        TESTS::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "HiGHS")
            sage: p.add_variables(2)              
            1
            sage: p.add_linear_constraint([(0, 1), (1, 1)], None, 2.0)
            sage: p.set_objective([1, 1])         
            sage: p.solve()                       
            0
            sage: p.objective_coefficient(0,1)    
            sage: p.solve()                       
            0
        """
        import highspy
        
        status = self.highs_model.run()
        
        # Check if solution was successful
        model_status = self.highs_model.getModelStatus()
        
        if model_status == highspy.HighsModelStatus.kOptimal:
            return 0  # Success
        elif model_status == highspy.HighsModelStatus.kInfeasible:
            raise MIPSolverException("HiGHS: Problem is infeasible")
        elif model_status == highspy.HighsModelStatus.kUnbounded:
            raise MIPSolverException("HiGHS: Problem is unbounded")
        else:
            raise MIPSolverException(f"HiGHS: Solver failed with status {model_status}")
    
    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(2)            
            1
            sage: p.add_linear_constraint([(0, 1), (1, 1)], None, 2.0)
            sage: p.set_objective([1, 1])       
            sage: p.solve()                     
            0
            sage: p.get_objective_value()       
            2.0
        """
        info = self.highs_model.getInfo()
        return info.objective_function_value + self.obj_constant_term
    
    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable in the current solution.

        INPUT:

        - ``variable`` -- the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(2)            
            1
            sage: p.add_linear_constraint([(0, 1), (1, 1)], None, 2.0)
            sage: p.set_objective([1, 1])       
            sage: p.solve()                     
            0
            sage: p.get_variable_value(0)       
            2.0
            sage: p.get_variable_value(1)       
            0.0
        """
        solution = self.highs_model.getSolution()
        return solution.col_value[variable]
    
    cpdef int ncols(self) noexcept:
        """
        Return the number of columns/variables.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.ncols()                     
            0
            sage: p.add_variables(2)            
            1
            sage: p.ncols()                     
            2
        """
        return self.numcols
    
    cpdef int nrows(self) noexcept:
        """
        Return the number of rows/constraints.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.nrows()                     
            0
            sage: p.add_variables(2)            
            1
            sage: p.add_linear_constraint([(0, 1), (1, 1)], None, 2.0)
            sage: p.nrows()                     
            1
        """
        return self.numrows
    
    cpdef bint is_maximization(self) noexcept:
        """
        Test whether the problem is a maximization.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.is_maximization()           
            True
        """
        import highspy
        lp = self.highs_model.getLp()
        return lp.sense_ == highspy.ObjSense.kMaximize
    
    cpdef __copy__(self):
        """
        Return a copy of ``self``.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(2)            
            1
            sage: q = copy(p)                   
            sage: q.ncols()                     
            2
        """
        import highspy
        import tempfile
        import os
        
        cdef HiGHSBackend p = type(self)(maximization=self.is_maximization())
        
        # Copy by writing and reading through a temporary MPS file
        # since Highs objects cannot be directly copied
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mps', delete=False) as f:
            temp_file = f.name
        
        try:
            self.highs_model.writeModel(temp_file)
            p.highs_model.readModel(temp_file)
            p.highs_model.setOptionValue("log_to_console", False)
        finally:
            if os.path.exists(temp_file):
                os.unlink(temp_file)
        
        p.prob_name = self.prob_name
        p.col_name_var = copy(self.col_name_var)
        p.row_name_var = copy(self.row_name_var)
        p.numcols = self.numcols
        p.numrows = self.numrows
        p.obj_constant_term = self.obj_constant_term
        return p
    
    cpdef variable_upper_bound(self, int index, value=None):
        """
        Set or get the upper bound of a variable.

        INPUT:

        - ``index`` -- the variable's id
        - ``value`` -- real value or ``None`` to get the current value

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variable()              
            0
            sage: p.variable_upper_bound(0)     
            sage: p.variable_upper_bound(0, 10.0)
            sage: p.variable_upper_bound(0)     
            10.0
        """
        import highspy
        
        if value is None:
            lp = self.highs_model.getLp()
            ub = lp.col_upper_[index]
            if ub >= highspy.kHighsInf:
                return None
            return ub
        else:
            lb = self.variable_lower_bound(index)
            if lb is None:
                lb = -highspy.kHighsInf
            if value is None:
                value = highspy.kHighsInf
            else:
                value = float(value)
            self.highs_model.changeColBounds(index, lb, value)
    
    cpdef variable_lower_bound(self, int index, value=None):
        """
        Set or get the lower bound of a variable.

        INPUT:

        - ``index`` -- the variable's id
        - ``value`` -- real value or ``None`` to get the current value

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variable()              
            0
            sage: p.variable_lower_bound(0)     
            0.0
            sage: p.variable_lower_bound(0, -10.0)
            sage: p.variable_lower_bound(0)     
            -10.0
        """
        import highspy
        
        if value is None:
            lp = self.highs_model.getLp()
            lb = lp.col_lower_[index]
            if lb <= -highspy.kHighsInf:
                return None
            return lb
        else:
            ub = self.variable_upper_bound(index)
            if ub is None:
                ub = highspy.kHighsInf
            if value is None:
                value = -highspy.kHighsInf
            else:
                value = float(value)
            self.highs_model.changeColBounds(index, value, ub)
    
    cpdef col_name(self, int index):
        """
        Return the ``index``-th column name.

        INPUT:

        - ``index`` -- integer; the column's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variable(name='x')
            0
            sage: p.col_name(0)
            'x'
        """
        if index < 0 or index >= self.ncols():
            raise ValueError(f"invalid column index {index}")
        return self.col_name_var.get(index, f"x_{index}")
    
    cpdef col_bounds(self, int index):
        """
        Return the bounds of a specific variable.

        INPUT:

        - ``index`` -- integer; the variable's id

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the variable is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0.0, 5.0)
        """
        import highspy
        if index < 0 or index >= self.ncols():
            raise ValueError(f"invalid column index {index}")
        
        lp = self.highs_model.getLp()
        lower = lp.col_lower_[index]
        upper = lp.col_upper_[index]
        
        # HiGHS uses very large negative/positive values for infinity
        # Check with a threshold rather than exact comparison
        if lower <= -1e20:
            lower = None
        if upper >= 1e20:
            upper = None
            
        return (lower, upper)
    
    cpdef bint is_variable_integer(self, int index) noexcept:
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variable()
            0
            sage: p.is_variable_integer(0)
            False
            sage: p.add_variable(integer=True)
            1
            sage: p.is_variable_integer(1)
            True
        """
        import highspy
        if index < 0 or index >= self.ncols():
            return False
        lp = self.highs_model.getLp()
        if index >= len(lp.integrality_):
            return False
        return lp.integrality_[index] == highspy.HighsVarType.kInteger
    
    cpdef bint is_variable_binary(self, int index) noexcept:
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variable()
            0
            sage: p.is_variable_binary(0)
            False
            sage: p.add_variable(binary=True)
            1
            sage: p.is_variable_binary(1)
            True
        """
        import highspy
        if index < 0 or index >= self.ncols():
            return False
        lp = self.highs_model.getLp()
        if index >= len(lp.integrality_):
            return False
        if lp.integrality_[index] == highspy.HighsVarType.kInteger:
            lower = lp.col_lower_[index]
            upper = lp.col_upper_[index]
            return lower == 0.0 and upper == 1.0
        return False
    
    cpdef bint is_variable_continuous(self, int index) noexcept:
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` -- integer; the variable's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variable()
            0
            sage: p.is_variable_continuous(0)
            True
        """
        import highspy
        if index < 0 or index >= self.ncols():
            return False
        lp = self.highs_model.getLp()
        if index >= len(lp.integrality_):
            return True  # Default to continuous
        return lp.integrality_[index] == highspy.HighsVarType.kContinuous
    
    cpdef set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable.

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``vtype`` -- integer:

            * `1` = Integer
            * `0` = Binary
            * `-1` = Real (Continuous)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variable()
            0
            sage: p.is_variable_continuous(0)
            True
            sage: p.set_variable_type(0, 1)
            sage: p.is_variable_integer(0)
            True
        """
        import highspy
        if vtype == 1:
            # Integer
            self.highs_model.changeColIntegrality(variable, highspy.HighsVarType.kInteger)
        elif vtype == 0:
            # Binary
            self.highs_model.changeColIntegrality(variable, highspy.HighsVarType.kInteger)
            self.highs_model.changeColBounds(variable, 0.0, 1.0)
        elif vtype == -1:
            # Real (Continuous)
            self.highs_model.changeColIntegrality(variable, highspy.HighsVarType.kContinuous)
        else:
            raise ValueError(f"Unknown variable type {vtype}")
    
    cpdef write_lp(self, filename):
        """
        Write the problem to a .lp file.

        INPUT:

        - ``filename`` -- string; the file name

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(2)            
            1
            sage: p.add_linear_constraint([(0, 1), (1, 1)], None, 2.0)
            sage: import tempfile                                     
            sage: with tempfile.NamedTemporaryFile(suffix='.lp') as f:
            ....:     p.write_lp(f.name)        
        """
        self.highs_model.writeModel(str(filename))
    
    cpdef write_mps(self, filename, int modern):
        """
        Write the problem to a .mps file.

        INPUT:

        - ``filename`` -- string; the file name
        - ``modern`` -- integer; whether to use modern MPS format (ignored for HiGHS)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(2)            
            1
            sage: p.add_linear_constraint([(0, 1), (1, 1)], None, 2.0)
            sage: import tempfile                                     
            sage: with tempfile.NamedTemporaryFile(suffix='.mps') as f:
            ....:     p.write_mps(f.name, 1)    
        """
        self.highs_model.writeModel(str(filename))
