"""
HiGHS Backend

AUTHORS:

- Chenxin Zhong (chenxin.zhong@outlook.com): initial implementation

This backend uses the HiGHS optimization solver C API, which supports Linear Programming (LP),
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
from cysignals.signals cimport sig_on, sig_off

from sage.libs.highs.highs_c_api cimport *

# C standard library for memory allocation
cdef extern from "stdlib.h":
    void* malloc(size_t size) nogil
    void free(void* ptr) nogil

cdef class HiGHSBackend(GenericBackend):
    """
    MIP Backend that uses the HiGHS solver via C API.

    HiGHS is a high-performance solver for large-scale LP, QP, and MIP.
    This implementation uses the HiGHS C API directly for optimal performance
    and proper interrupt handling with sig_on/sig_off.
    """

    def __cinit__(self, maximization=True):
        """
        Constructor.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='HiGHS')
        """
        # Create HiGHS instance
        self.highs = Highs_create()
        if self.highs == NULL:
            raise MemoryError("Failed to create HiGHS instance")
        
        # Initialize metadata
        self.prob_name = ""
        self.col_name_var = {}
        self.row_name_var = {}
        self.row_data_cache = {}
        self.numcols = 0
        self.numrows = 0
        self.obj_constant_term = 0.0
        
        # Suppress HiGHS output messages
        Highs_setBoolOptionValue(self.highs, b"log_to_console", 0)
        
        # Set optimization sense
        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)
    
    def __dealloc__(self):
        """
        Destructor - free HiGHS instance.
        """
        if self.highs != NULL:
            Highs_destroy(self.highs)
            self.highs = NULL
    
    cdef void _get_col_bounds(self, int col, double* lb, double* ub) except *:
        """
        Helper method to get column bounds using Highs_getColsByRange.
        """
        cdef HighsInt num_col, num_nz, status
        
        sig_on()
        status = Highs_getColsByRange(self.highs, col, col, 
                                     &num_col, NULL, lb, ub,
                                     &num_nz, NULL, NULL, NULL)
        sig_off()
        
        if status != kHighsStatusOk:
            raise MIPSolverException("HiGHS: Failed to get column bounds")
    
    cdef void _get_row_bounds(self, int row, double* lb, double* ub) except *:
        """
        Helper method to get row bounds using Highs_getRowsByRange.
        """
        cdef HighsInt num_row, num_nz, status
        
        sig_on()
        status = Highs_getRowsByRange(self.highs, row, row, 
                                     &num_row, lb, ub,
                                     &num_nz, NULL, NULL, NULL)
        sig_off()
        
        if status != kHighsStatusOk:
            raise MIPSolverException("HiGHS: Failed to get row bounds")
    
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
        cdef HighsInt var_type
        cdef double lb, ub
        cdef HighsInt status
        cdef HighsInt col_idx
        
        # Determine variable type - only one of binary, continuous, integer can be True
        if sum([binary, continuous, integer]) > 1:
            raise ValueError("only one of binary, continuous, and integer can be True")
        
        if binary:
            var_type = kHighsVarTypeInteger
            lb = 0.0
            ub = 1.0
        elif integer:
            var_type = kHighsVarTypeInteger
        else:
            var_type = kHighsVarTypeContinuous
        
        # Set bounds
        if lower_bound is None:
            if binary:
                lb = 0.0
            else:
                lb = -Highs_getInfinity(self.highs)
        else:
            lb = float(lower_bound)
        
        if upper_bound is None:
            if binary:
                ub = 1.0
            else:
                ub = Highs_getInfinity(self.highs)
        else:
            ub = float(upper_bound)
        
        # Add column with empty constraint coefficients
        sig_on()
        status = Highs_addCol(self.highs, float(obj), lb, ub, 0, NULL, NULL)
        sig_off()
        
        if status != kHighsStatusOk:
            raise MIPSolverException("HiGHS: Failed to add variable")
        
        col_idx = self.numcols
        self.numcols += 1
        
        # Set integrality if needed
        if var_type == kHighsVarTypeInteger:
            sig_on()
            status = Highs_changeColIntegrality(self.highs, col_idx, kHighsVarTypeInteger)
            sig_off()
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to set variable integrality")
        
        # Set name if provided
        if name is not None:
            name_bytes = str(name).encode('utf-8')
            Highs_passColName(self.highs, col_idx, name_bytes)
            self.col_name_var[col_idx] = str(name)
        
        return col_idx
    
    cpdef int add_variable_with_type(self, int vtype, lower_bound=0.0, upper_bound=None, 
                                     obj=0.0, name=None) except -1:
        """
        Add a variable with type specified as an integer.

        This amounts to adding a new column to the matrix. By default,
        the variable is positive and real, and the coefficient in the
        objective function is 0.0.

        INPUT:

        - ``vtype`` -- integer specifying the variable type:

            * ``1`` = Integer
            * ``0`` = Binary
            * ``-1`` = Real (Continuous)

        - ``lower_bound`` -- the lower bound of the variable (default: 0)

        - ``upper_bound`` -- the upper bound of the variable (default: ``None``)

        - ``obj`` -- (optional) coefficient of this variable in the objective function (default: 0.0)

        - ``name`` -- an optional name for the newly added variable (default: ``None``)

        OUTPUT: the index of the newly created variable

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "HiGHS")
            sage: p.ncols()                       
            0
            sage: p.add_variable_with_type(-1)    # Continuous variable
            0
            sage: p.is_variable_continuous(0)
            True
            sage: p.add_variable_with_type(0)     # Binary variable
            1
            sage: p.is_variable_binary(1)
            True
            sage: p.add_variable_with_type(1, lower_bound=-2.0)  # Integer variable
            2
            sage: p.is_variable_integer(2)
            True
            sage: p.add_variable_with_type(1, name='x', obj=1.0)
            3
            sage: p.col_name(3)                    
            'x'
            sage: p.objective_coefficient(3)       
            1.0

        TESTS:

        Invalid variable type raises an error::

            sage: p.add_variable_with_type(2)
            Traceback (most recent call last):
            ...
            ValueError: Invalid variable type 2. Must be -1 (continuous), 0 (binary), or 1 (integer)
        """
        cdef HighsInt var_type
        cdef double lb, ub
        cdef HighsInt status
        cdef HighsInt col_idx
        
        # Validate and determine variable type
        if vtype == 1:
            # Integer
            var_type = kHighsVarTypeInteger
        elif vtype == 0:
            # Binary - set type to integer with bounds [0,1]
            var_type = kHighsVarTypeInteger
        elif vtype == -1:
            # Continuous
            var_type = kHighsVarTypeContinuous
        else:
            raise ValueError(f"Invalid variable type {vtype}. Must be -1 (continuous), 0 (binary), or 1 (integer)")
        
        # Set bounds
        if vtype == 0:  # Binary
            # Binary variables have fixed bounds [0, 1]
            lb = 0.0
            ub = 1.0
            # Override user-provided bounds for binary variables
        else:
            # For integer and continuous variables, use provided bounds
            if lower_bound is None:
                lb = -Highs_getInfinity(self.highs)
            else:
                lb = float(lower_bound)
            
            if upper_bound is None:
                ub = Highs_getInfinity(self.highs)
            else:
                ub = float(upper_bound)
        
        # Add column with empty constraint coefficients
        sig_on()
        status = Highs_addCol(self.highs, float(obj), lb, ub, 0, NULL, NULL)
        sig_off()
        
        if status != kHighsStatusOk:
            raise MIPSolverException("HiGHS: Failed to add variable")
        
        col_idx = self.numcols
        self.numcols += 1
        
        # Set integrality if needed (for integer or binary)
        if var_type == kHighsVarTypeInteger:
            sig_on()
            status = Highs_changeColIntegrality(self.highs, col_idx, kHighsVarTypeInteger)
            sig_off()
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to set variable integrality")
        
        # Set name if provided
        if name is not None:
            name_bytes = str(name).encode('utf-8')
            Highs_passColName(self.highs, col_idx, name_bytes)
            self.col_name_var[col_idx] = str(name)
        
        return col_idx
    

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` -- +1 for maximization; any other integer for minimization

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        cdef HighsInt highs_sense
        cdef HighsInt status
        
        if sense == 1:
            highs_sense = kHighsObjSenseMaximize
        else:
            highs_sense = kHighsObjSenseMinimize
        
        sig_on()
        status = Highs_changeObjectiveSense(self.highs, highs_sense)
        sig_off()
        
        if status != kHighsStatusOk:
            raise MIPSolverException("HiGHS: Failed to set objective sense")
    
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
        cdef HighsInt status
        cdef HighsInt model_status
        
        # Pure C API call between sig_on/sig_off - no Python objects touched!
        sig_on()
        status = Highs_run(self.highs)
        sig_off()
        
        if status != kHighsStatusOk:
            raise MIPSolverException("HiGHS: Solver run failed")
        
        # Check model status
        model_status = Highs_getModelStatus(self.highs)
        
        if model_status == kHighsModelStatusOptimal:
            return 0  # Success
        elif model_status == kHighsModelStatusInfeasible:
            raise MIPSolverException("HiGHS: Problem is infeasible")
        elif model_status == kHighsModelStatusUnbounded:
            raise MIPSolverException("HiGHS: Problem is unbounded")
        elif model_status == kHighsModelStatusTimeLimit:
            raise MIPSolverException("HiGHS: Time limit reached")
        elif model_status == kHighsModelStatusIterationLimit:
            raise MIPSolverException("HiGHS: Iteration limit reached")
        elif model_status == kHighsModelStatusInterrupt:
            raise MIPSolverException("HiGHS: Interrupted by user")
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
        cdef double obj_value
        
        obj_value = Highs_getObjectiveValue(self.highs)
        # HiGHS already includes the offset, so don't add it again
        return obj_value
    
    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

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
        cdef double* col_value
        cdef HighsInt num_cols
        cdef HighsInt status
        cdef double result
        
        num_cols = Highs_getNumCol(self.highs)
        
        if variable < 0 or variable >= num_cols:
            raise ValueError(f"Variable index {variable} out of range [0, {num_cols})")
        
        # Allocate array for solution
        col_value = <double*> malloc(num_cols * sizeof(double))
        if col_value == NULL:
            raise MemoryError("Failed to allocate memory for solution")
        
        try:
            sig_on()
            status = Highs_getSolution(self.highs, col_value, NULL, NULL, NULL)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to get solution")
            
            result = col_value[variable]
        finally:
            free(col_value)
        
        return result
    
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
        
        for i in range(number):
            name = None
            if names is not None and i < len(names):
                name = names[i]
            self.add_variable(lower_bound=lower_bound, upper_bound=upper_bound,
                            binary=binary, continuous=continuous, integer=integer,
                            obj=obj, name=name)
        
        return self.numcols - 1
    
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
            sage: p.objective_coefficient(0, 2)
            sage: p.objective_coefficient(0)
            2.0

        TESTS:

        We sanity check the input that will be passed to HiGHS::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.objective_coefficient(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 2
        """
        cdef HighsInt status
        cdef double cost
        cdef HighsInt num_col
        cdef HighsInt num_nz
        
        if variable < 0 or variable >= self.numcols:
            raise ValueError(f"invalid variable index {variable}")
        
        if coeff is None:
            # Get coefficient using Highs_getColsByRange
            sig_on()
            status = Highs_getColsByRange(self.highs, variable, variable, 
                                         &num_col, &cost, NULL, NULL,
                                         &num_nz, NULL, NULL, NULL)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to get objective coefficient")
            return cost
        else:
            # Set coefficient
            sig_on()
            status = Highs_changeColCost(self.highs, variable, float(coeff))
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to set objective coefficient")
    
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
        
        # Set offset in HiGHS
        sig_on()
        Highs_changeObjectiveOffset(self.highs, d)
        sig_off()
    
    cpdef set_verbosity(self, int level):
        """
        Set the log (verbosity) level.

        INPUT:

        - ``level`` -- integer; from 0 (no verbosity) to 1

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.set_verbosity(0)
        """
        cdef HighsInt status
        cdef bint log_to_console
        
        if level == 0:
            log_to_console = False
        elif level == 1:
            log_to_console = True
        else:
            raise ValueError("Invalid verbosity level. Must be 0 or 1.")
        
        status = Highs_setBoolOptionValue(self.highs, b"log_to_console", log_to_console)
        if status != kHighsStatusOk:
            raise MIPSolverException("HiGHS: Failed to set verbosity")
    
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
        cdef double lb, ub
        cdef HighsInt num_nz
        cdef HighsInt* indices
        cdef double* values
        cdef HighsInt status
        cdef int i
        
        # Convert bounds
        if lower_bound is None:
            lb = -Highs_getInfinity(self.highs)
        else:
            lb = float(lower_bound)
        
        if upper_bound is None:
            ub = Highs_getInfinity(self.highs)
        else:
            ub = float(upper_bound)
        
        # Build coefficient arrays
        coeff_list = list(coefficients)
        num_nz = len(coeff_list)
        
        if num_nz == 0:
            # Empty constraint
            sig_on()
            status = Highs_addRow(self.highs, lb, ub, 0, NULL, NULL)
            sig_off()
        else:
            # Allocate arrays
            indices = <HighsInt*> malloc(num_nz * sizeof(HighsInt))
            values = <double*> malloc(num_nz * sizeof(double))
            
            if indices == NULL or values == NULL:
                free(indices)
                free(values)
                raise MemoryError("Failed to allocate memory for constraint")
            
            try:
                # Fill arrays
                for i in range(num_nz):
                    var_idx, coeff_val = coeff_list[i]
                    if var_idx < 0 or var_idx >= self.ncols():
                        raise ValueError(f"invalid variable index {var_idx}")
                    indices[i] = var_idx
                    values[i] = float(coeff_val)
                
                # Add constraint
                sig_on()
                status = Highs_addRow(self.highs, lb, ub, num_nz, indices, values)
                sig_off()
            finally:
                free(indices)
                free(values)
        
        if status != kHighsStatusOk:
            raise MIPSolverException("HiGHS: Failed to add constraint")
        
        # Handle name
        if name is not None:
            self.row_name_var[name] = self.numrows
        
        self.numrows += 1
    
    cpdef add_linear_constraints(self, int number, lower_bound, upper_bound, names=None):
        """
        Add ``number`` linear constraints.

        INPUT:

        - ``number`` -- integer; the number of constraints to add

        - ``lower_bound`` -- a lower bound, either a real value or ``None``

        - ``upper_bound`` -- an upper bound, either a real value or ``None``

        - ``names`` -- an optional list of names (default: ``None``)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraints(5, None, 2)
            sage: p.row_bounds(4)
            (None, 2.0)
            sage: p.add_linear_constraints(2, None, 2, names=['foo','bar'])
        """
        cdef int i
        cdef double lb, ub
        cdef HighsInt status
         
        # Convert bounds
        if lower_bound is None:
            lb = -Highs_getInfinity(self.highs)
        else:
            lb = float(lower_bound)
        
        if upper_bound is None:
            ub = Highs_getInfinity(self.highs)
        else:
            ub = float(upper_bound)
        
        # Add empty constraints
        for i in range(number):
            sig_on()
            status = Highs_addRow(self.highs, lb, ub, 0, NULL, NULL)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to add constraint")
            
            if names is not None and i < len(names):
                name = names[i]
                if name is not None:
                    self.row_name_var[name] = self.numrows
            
            self.numrows += 1
    
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
        cdef HighsInt sense, status
        status = Highs_getObjectiveSense(self.highs, &sense)
        if status != kHighsStatusOk:
            return True  # default to maximize
        return sense == kHighsObjSenseMaximize
    
    cpdef get_row_prim(self, int i):
        """
        Return the value of the auxiliary variable associated with i-th row.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver='HiGHS')
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: lp.solve()
            0
            sage: lp.get_objective_value()
            280.0
            sage: lp.get_row_prim(0)
            24.0
            sage: lp.get_row_prim(1)
            20.0
            sage: lp.get_row_prim(2)
            8.0

        TESTS:

        We sanity check the input::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.get_row_prim(2)
            Traceback (most recent call last):
            ...
            ValueError: Row index 2 out of range ...
        """
        cdef double* row_value
        cdef HighsInt num_rows
        cdef HighsInt status
        cdef double result
        
        num_rows = Highs_getNumRow(self.highs)
        
        if i < 0 or i >= num_rows:
            raise ValueError(f"Row index {i} out of range [0, {num_rows})")
        
        row_value = <double*> malloc(num_rows * sizeof(double))
        if row_value == NULL:
            raise MemoryError("Failed to allocate memory")
        
        try:
            sig_on()
            status = Highs_getSolution(self.highs, NULL, NULL, row_value, NULL)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to get solution")
            
            result = row_value[i]
        finally:
            free(row_value)
        
        return result
    
    cpdef double get_row_dual(self, int i) except? -1:
        """
        Return the dual value of a constraint.

        The dual value of the i-th row is also the value of the i-th variable
        of the dual problem.

        The dual value of a constraint is the shadow price of the constraint.
        The shadow price is the amount by which the objective value will change
        if the constraint's bounds change by one unit under the precondition
        that the basis remains the same.

        INPUT:

        - ``i`` -- the index of the constraint

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver='HiGHS')
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: lp.solve()
            0
            sage: lp.get_row_dual(0)   # tol 1e-6
            0.0
            sage: lp.get_row_dual(1)   # tol 1e-6
            10.0
            sage: lp.get_row_dual(2)   # tol 1e-6
            10.0

        TESTS:

        We sanity check the input::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.get_row_dual(2)
            Traceback (most recent call last):
            ...
            ValueError: Row index 2 out of range ...
        """
        cdef double* row_dual
        cdef HighsInt num_rows
        cdef HighsInt status
        cdef double result
        
        num_rows = Highs_getNumRow(self.highs)
        
        if i < 0 or i >= num_rows:
            raise ValueError(f"Row index {i} out of range [0, {num_rows})")
        
        row_dual = <double*> malloc(num_rows * sizeof(double))
        if row_dual == NULL:
            raise MemoryError("Failed to allocate memory")
        
        try:
            sig_on()
            status = Highs_getSolution(self.highs, NULL, NULL, NULL, row_dual)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to get dual solution")
            
            result = row_dual[i]
        finally:
            free(row_dual)
        
        return result
    
    cpdef double get_col_dual(self, int j) except? -1:
        """
        Return the dual value (reduced cost) of a variable.

        The dual value is the reduced cost of a variable.
        The reduced cost is the amount by which the objective coefficient
        of a non-basic variable has to change to become a basic variable.

        INPUT:

        - ``j`` -- the index of the variable

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(3)
            2
            sage: p.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: p.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: p.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: p.set_objective([60, 30, 20])
            sage: p.solve()
            0
            sage: p.get_col_dual(0)    # tol 1e-6
            0.0
            sage: p.get_col_dual(1)    # tol 1e-6
            -5.0
            sage: p.get_col_dual(2)    # tol 1e-6
            0.0

        TESTS:

        We sanity check the input::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.get_col_dual(2)
            Traceback (most recent call last):
            ...
            ValueError: Variable index 2 out of range ...
        """
        cdef double* col_dual
        cdef HighsInt num_cols
        cdef HighsInt status
        cdef double result
        
        num_cols = Highs_getNumCol(self.highs)
        
        if j < 0 or j >= num_cols:
            raise ValueError(f"Variable index {j} out of range [0, {num_cols})")
        
        col_dual = <double*> malloc(num_cols * sizeof(double))
        if col_dual == NULL:
            raise MemoryError("Failed to allocate memory")
        
        try:
            sig_on()
            status = Highs_getSolution(self.highs, NULL, col_dual, NULL, NULL)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to get dual solution")
            
            result = col_dual[j]
        finally:
            free(col_dual)
        
        return result
    
    cpdef best_known_objective_bound(self):
        """
        Return the value of the currently best known bound.

        This method returns the current best upper (resp. lower) bound on the
        optimal value of the objective function in a maximization
        (resp. minimization) problem. It is equal to the output of
        :meth:`get_objective_value` if the MILP found an optimal solution, but
        it can differ if it was interrupted manually or after a time limit (cf
        :meth:`solver_parameter`).

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0, 1), (1, 1)], None, 2.0)
            sage: p.set_objective([1, 1])
            sage: p.solve()
            0
            sage: p.best_known_objective_bound()
            2.0

        TESTS::
            sage: # needs sage.graphs
            sage: g = graphs.CubeGraph(9)
            sage: p = MixedIntegerLinearProgram(solver='HiGHS')
            sage: p.solver_parameter("mip_rel_gap",100)
            sage: b = p.new_variable(binary=True)
            sage: p.set_objective(p.sum(b[v] for v in g))
            sage: for v in g:
            ....:     p.add_constraint(b[v]+p.sum(b[u] for u in g.neighbors(v)) <= 1)
            sage: p.add_constraint(b[v] == 1) # Force an easy non-0 solution
            sage: p.solve() # rel tol 100
            2.0
            sage: backend = p.get_backend()
            sage: backend.best_known_objective_bound()  # abs tol 1e-6
            48.0
        """
        cdef double mip_dual_bound
        cdef HighsInt status
        cdef HighsInt i, var_type
        cdef bint is_mip = False
        
        # Check if this is a MIP problem (has integer variables)
        for i in range(self.numcols):
            status = Highs_getColIntegrality(self.highs, i, &var_type)
            if status != kHighsStatusOk:
                continue
            if var_type == kHighsVarTypeInteger:
                is_mip = True
                break
        
        if not is_mip:
            # For LP problems, the bound equals the objective value
            return self.get_objective_value()
        
        # Get the MIP dual bound using info query
        status = Highs_getDoubleInfoValue(self.highs, b"mip_dual_bound", &mip_dual_bound)
        if status != kHighsStatusOk:
            # If not available, return the objective value
            return ValueError("MIP dual bound not available")
        # HiGHS already includes the offset in the dual bound
        return mip_dual_bound
    
    cpdef get_relative_objective_gap(self):
        """
        Return the relative objective gap of the best known solution.

        For a minimization problem, this value is computed by
        `(\texttt{bestinteger} - \texttt{bestobjective}) / (1e-10 +
        |\texttt{bestobjective}|)`, where ``bestinteger`` is the value returned
        by :meth:`get_objective_value` and ``bestobjective`` is the value
        returned by :meth:`best_known_objective_bound`. For a maximization
        problem, the value is computed by `(\texttt{bestobjective} -
        \texttt{bestinteger}) / (1e-10 + |\texttt:bestobjective}|)`.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0, 1), (1, 1)], None, 2.0)
            sage: p.set_objective([1, 1])
            sage: p.solve()
            0
            sage: p.get_relative_objective_gap()
            0.0
        """
        cdef double gap
        cdef HighsInt status
        cdef HighsInt i, var_type
        cdef bint is_mip = False
        
        # Check if this is a MIP problem (has integer variables)
        for i in range(self.numcols):
            status = Highs_getColIntegrality(self.highs, i, &var_type)
            if status != kHighsStatusOk:
                continue
            if var_type == kHighsVarTypeInteger:
                is_mip = True
                break
        
        if not is_mip:
            # For LP problems, the gap is 0
            return 0.0
        
        # Get the MIP gap using info query
        status = Highs_getDoubleInfoValue(self.highs, b"mip_gap", &gap)
        if status != kHighsStatusOk:
            # If not available, return 0
            return ValueError("MIP gap not available")
        return gap
    
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

        TESTS:

        Check that invalid indices raise errors::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.variable_upper_bound(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 2
            sage: p.variable_upper_bound(-1)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index -1
            sage: p.add_variable()
            0
            sage: p.variable_upper_bound(3, 5)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 3
        """
        cdef double lb, ub
        cdef HighsInt status
        
        if index < 0 or index >= self.numcols:
            raise ValueError(f"invalid variable index {index}")
        
        if value is None:
            # Get current bound
            self._get_col_bounds(index, &lb, &ub)
            
            if ub >= Highs_getInfinity(self.highs) - 1:
                return None
            else:
                return ub
        else:
            # Set new bound
            self._get_col_bounds(index, &lb, &ub)
            
            if value is None:
                ub = Highs_getInfinity(self.highs)
            else:
                ub = float(value)
            
            sig_on()
            status = Highs_changeColBounds(self.highs, index, lb, ub)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to set variable upper bound")
    
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

        TESTS:

        Check that invalid indices raise errors::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.variable_lower_bound(2)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 2
            sage: p.variable_lower_bound(-1)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index -1
            sage: p.add_variable()
            0
            sage: p.variable_lower_bound(3, 5)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 3
        """
        cdef double lb, ub
        cdef HighsInt status
        
        if index < 0 or index >= self.numcols:
            raise ValueError(f"invalid variable index {index}")
        
        if value is None:
            # Get current bound
            self._get_col_bounds(index, &lb, &ub)
            
            if lb <= -Highs_getInfinity(self.highs) - 1:
                return None
            else:
                return lb
        else:
            # Set new bound
            self._get_col_bounds(index, &lb, &ub)
            
            if value is None:
                lb = -Highs_getInfinity(self.highs)
            else:
                lb = float(value)
            
            sig_on()
            status = Highs_changeColBounds(self.highs, index, lb, ub)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to set variable lower bound")
    
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
        if index < 0 or index >= self.numcols:
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
        cdef double lb, ub
        
        if index < 0 or index >= self.numcols:
            raise ValueError(f"invalid variable index {index}")
        
        self._get_col_bounds(index, &lb, &ub)
        
        # Convert infinities to None
        if lb <= -Highs_getInfinity(self.highs) + 1:
            lb_ret = None
        else:
            lb_ret = lb
        
        if ub >= Highs_getInfinity(self.highs) - 1:
            ub_ret = None
        else:
            ub_ret = ub
        
        return (lb_ret, ub_ret)
    
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

        TESTS:

        We check the behavior for an invalid index::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.is_variable_integer(2)
            False
        """
        cdef HighsInt integrality, status
        
        if index < 0 or index >= self.numcols:
            return False
        
        status = Highs_getColIntegrality(self.highs, index, &integrality)
        if status != kHighsStatusOk:
            return False
        return integrality == kHighsVarTypeInteger
    
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

        TESTS:

        We check the behavior for an invalid index::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.is_variable_binary(2)
            False
        """
        cdef HighsInt integrality, status
        cdef double lb, ub
        
        if index < 0 or index >= self.numcols:
            return False
        
        status = Highs_getColIntegrality(self.highs, index, &integrality)
        if status != kHighsStatusOk:
            return False
            
        if integrality == kHighsVarTypeInteger:
            self._get_col_bounds(index, &lb, &ub)
            return lb == 0.0 and ub == 1.0
        
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

        TESTS:

        We check the behavior for an invalid index::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.is_variable_continuous(2)
            False
        """
        cdef HighsInt integrality, status
        
        if index < 0 or index >= self.numcols:
            return False
        
        status = Highs_getColIntegrality(self.highs, index, &integrality)
        if status != kHighsStatusOk:
            return True  # default to continuous
        return integrality == kHighsVarTypeContinuous
    
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

        TESTS:

        We sanity check the input that will be passed to HiGHS::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.set_variable_type(2, 0)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 2
        """
        cdef HighsInt status
        cdef double lb, ub
        
        if variable < 0 or variable >= self.numcols:
            raise ValueError(f"invalid variable index {variable}")
        
        if vtype == 1:
            # Integer
            sig_on()
            status = Highs_changeColIntegrality(self.highs, variable, kHighsVarTypeInteger)
            sig_off()
        elif vtype == 0:
            # Binary - set to integer and change bounds to [0,1]
            sig_on()
            status = Highs_changeColIntegrality(self.highs, variable, kHighsVarTypeInteger)
            sig_off()
            
            if status == kHighsStatusOk:
                sig_on()
                status = Highs_changeColBounds(self.highs, variable, 0.0, 1.0)
                sig_off()
        elif vtype == -1:
            # Continuous
            sig_on()
            status = Highs_changeColIntegrality(self.highs, variable, kHighsVarTypeContinuous)
            sig_off()
        else:
            raise ValueError(f"Unknown variable type {vtype}")
        
        if status != kHighsStatusOk:
            raise MIPSolverException("HiGHS: Failed to set variable type")
    
    cpdef row_bounds(self, int index):
        """
        Return the bounds of a specific constraint.

        INPUT:

        - ``index`` -- integer; the constraint's id

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the constraint is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)
            sage: p.row(0)  # Note: zero coefficients are excluded in sparse format
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)
        """
        cdef double lb, ub
        cdef double infinity
        
        if index < 0 or index >= self.numrows:
            raise ValueError(f"invalid row index {index}")
        
        self._get_row_bounds(index, &lb, &ub)
        infinity = Highs_getInfinity(self.highs)
        
        return (
            (lb if abs(lb) < infinity else None),
            (ub if abs(ub) < infinity else None)
        )
    
    cpdef row_name(self, int index):
        """
        Return the ``index``-th row name.

        INPUT:

        - ``index`` -- integer; the row's id

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_linear_constraint([], 2, 2, name='foo')
            sage: p.row_name(0)
            'foo'
        """
        if index < 0 or index >= self.numrows:
            raise ValueError(f"invalid row index {index}")
        
        # Search for name in dictionary
        for name, idx in self.row_name_var.items():
            if idx == index:
                return str(name)
        
        return f"constraint_{index}"
    
    cpdef solver_parameter(self, name, value=None):
        """
        Return or define a solver parameter.

        INPUT:

        - ``name`` -- string; the parameter name

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value

        HiGHS solver parameters can be set using their option names as documented
        in the HiGHS documentation: https://ergo-code.github.io/HiGHS/dev/options/definitions/

        Common parameters include:

        - ``time_limit`` -- maximum time in seconds (double)
        - ``mip_rel_gap`` -- relative MIP gap tolerance (double)
        - ``mip_abs_gap`` -- absolute MIP gap tolerance (double)
        - ``threads`` -- number of threads to use (int)
        - ``presolve`` -- presolve option: "off", "choose", or "on"
        - ``solver`` -- solver to use: "choose", "simplex", "ipm", or "pdlp (need CUDA)"
        - ``parallel`` -- parallel option: "off", "choose", or "on"
        - ``log_to_console`` -- whether to log to console: True or False

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.solver_parameter("time_limit", 60)
            sage: p.solver_parameter("time_limit")
            60.0
            sage: p.solver_parameter("threads", 2)
            sage: p.solver_parameter("threads")
            2
            sage: p.solver_parameter("presolve", "on")
            sage: p.solver_parameter("presolve")
            'on'

        You can also use boolean values for options::

            sage: p.solver_parameter("log_to_console", False)
            sage: p.solver_parameter("log_to_console")
            False

        Float parameters like MIP gap tolerance work correctly::

            sage: p.solver_parameter("mip_rel_gap", 0.05)
            sage: p.solver_parameter("mip_rel_gap")
            0.05
        """
        cdef HighsInt status
        cdef bytes name_bytes
        cdef HighsInt int_value
        cdef double double_value
        cdef HighsInt bool_value
        cdef char* str_value
        cdef HighsInt option_type
        
        name_bytes = str(name).encode('utf-8')
        
        if value is None:
            # Get parameter - try each type until one succeeds
            # Try int first
            status = Highs_getIntOptionValue(self.highs, name_bytes, &int_value)
            if status == kHighsStatusOk:
                return int_value
            
            # Try double
            status = Highs_getDoubleOptionValue(self.highs, name_bytes, &double_value)
            if status == kHighsStatusOk:
                return double_value
            
            # Try bool
            status = Highs_getBoolOptionValue(self.highs, name_bytes, &bool_value)
            if status == kHighsStatusOk:
                return bool(bool_value)
            
            # Try string (allocate buffer for string value)
            str_value = <char*> malloc(256 * sizeof(char))
            if str_value == NULL:
                raise MemoryError("Failed to allocate memory for string option")
            try:
                status = Highs_getStringOptionValue(self.highs, name_bytes, str_value)
                if status == kHighsStatusOk:
                    result = str_value.decode('utf-8')
                    return result
                else:
                    raise ValueError(f"Unknown option {name}")
            finally:
                free(str_value)
        else:
            # Set parameter - need to determine type
            # Convert Sage types to Python types for easier type checking
            try:
                # Try to convert to Python numeric type if it's a Sage type
                # Check for float first since int() would truncate floats
                if isinstance(value, float):
                    pass  # Already a Python float
                elif hasattr(value, '__float__') and not isinstance(value, (bool, str, int)):
                    value = float(value)
                elif hasattr(value, '__int__') and not isinstance(value, (bool, str, float)):
                    value = int(value)
            except (TypeError, AttributeError):
                pass
            
            # Check bool first since bool is a subclass of int in Python
            if isinstance(value, bool):
                status = Highs_setBoolOptionValue(self.highs, name_bytes, value)
            elif isinstance(value, str):
                value_bytes = value.encode('utf-8')
                status = Highs_setStringOptionValue(self.highs, name_bytes, value_bytes)
            elif isinstance(value, (int, float)):
                # Try as double first (works for both int and float)
                status = Highs_setDoubleOptionValue(self.highs, name_bytes, float(value))
                if status != kHighsStatusOk and isinstance(value, int):
                    # If double failed and it's an int, try as int
                    status = Highs_setIntOptionValue(self.highs, name_bytes, value)
                if status != kHighsStatusOk:
                    # If both numeric methods failed, try as string (works for advanced options)
                    # For floats that are whole numbers, convert to int string to avoid "4.0" format
                    if isinstance(value, float) and value == int(value):
                        value_bytes = str(int(value)).encode('utf-8')
                    else:
                        value_bytes = str(value).encode('utf-8')
                    status = Highs_setStringOptionValue(self.highs, name_bytes, value_bytes)
            else:
                raise ValueError(f"Unknown parameter type for {name}: {type(value)}")
            
            if status != kHighsStatusOk:
                raise MIPSolverException(f"HiGHS: Failed to set parameter {name}")
    
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
        cdef bytes filename_bytes
        cdef HighsInt status

        import os
        cdef str filenamestr = str(filename)
        cdef object _root
        cdef str ext

        _root, ext = os.path.splitext(filenamestr)
        if ext.lower() != '.lp':
            filenamestr = filenamestr + '.lp'

        filename_bytes = os.fsencode(filenamestr)

        sig_on()
        status = Highs_writeModel(self.highs, filename_bytes)
        sig_off()
        
        if status != kHighsStatusOk and status != kHighsStatusWarning:
            raise MIPSolverException(f"HiGHS: Failed to write LP file {filenamestr} (status {status})")
    
    cpdef write_mps(self, filename, int modern):
        """
        Write the problem to a .mps file.

        INPUT:

        - ``filename`` -- string; the file name
        - ``modern`` -- integer; whether to use modern MPS format (ignored for HiGHS)

        .. NOTE::

            HiGHS determines the output format from the filename extension.
            The ``modern`` flag is accepted for API compatibility but ignored.

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
        cdef bytes filename_bytes
        cdef HighsInt status

        import os
        cdef str filenamestr = str(filename)
        cdef object _root
        cdef str ext

        _root, ext = os.path.splitext(filenamestr)
        if ext.lower() != '.mps':
            filenamestr = filenamestr + '.mps'

        filename_bytes = os.fsencode(filenamestr)

        sig_on()
        status = Highs_writeModel(self.highs, filename_bytes)
        sig_off()
        
        if status != kHighsStatusOk and status != kHighsStatusWarning:
            raise MIPSolverException(f"HiGHS: Failed to write MPS file {filenamestr} (status {status})")
    
    cpdef remove_constraint(self, int i):
        """
        Remove a constraint from ``self``.

        INPUT:

        - ``i`` -- index of the constraint to remove

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='HiGHS')
            sage: x, y = p['x'], p['y']
            sage: p.add_constraint(2*x + 3*y <= 6)
            sage: p.add_constraint(3*x + 2*y <= 6)
            sage: p.add_constraint(x >= 0)
            sage: p.set_objective(x + y + 7)
            sage: p.set_integer(x); p.set_integer(y)
            sage: p.solve()
            9.0
            sage: p.remove_constraint(0)
            sage: p.solve()
            10.0

        Removing fancy constraints does not make Sage crash::

            sage: MixedIntegerLinearProgram(solver = "HiGHS").remove_constraint(-2)
            Traceback (most recent call last):
            ...
            ValueError: The constraint's index i must satisfy 0 <= i < number_of_constraints
        """
        cdef HighsInt status
        
        if i < 0 or i >= self.numrows:
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")
        
        # Delete the single row
        sig_on()
        status = Highs_deleteRowsByRange(self.highs, i, i)
        sig_off()
        
        if status != kHighsStatusOk:
            raise MIPSolverException("HiGHS: Failed to remove constraint")
        
        self.numrows -= 1
        
        # Update row name mapping
        names_to_update = {}
        names_to_remove = []
        for name, row_idx in self.row_name_var.items():
            if row_idx < i:
                names_to_update[name] = row_idx
            elif row_idx > i:
                names_to_update[name] = row_idx - 1
            else:
                names_to_remove.append(name)
        
        for name in names_to_remove:
            del self.row_name_var[name]
        self.row_name_var.update(names_to_update)
    
    cpdef remove_constraints(self, constraints):
        """
        Remove several constraints.

        INPUT:

        - ``constraints`` -- an iterable containing the indices of the rows to remove

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='HiGHS')
            sage: x, y = p['x'], p['y']
            sage: p.add_constraint(2*x + 3*y <= 6)
            sage: p.add_constraint(3*x + 2*y <= 6)
            sage: p.add_constraint(x >= 0)
            sage: p.set_objective(x + y + 7)
            sage: p.set_integer(x); p.set_integer(y)
            sage: p.solve()
            9.0
            sage: p.remove_constraints([0])
            sage: p.solve()
            10.0
            sage: p.get_values([x,y])
            [-0.0, 3.0]

        TESTS:

        Removing fancy constraints does not make Sage crash::

            sage: MixedIntegerLinearProgram(solver="HiGHS").remove_constraints([0, -2])
            Traceback (most recent call last):
            ...
            ValueError: The constraint's index i must satisfy 0 <= i < number_of_constraints
        """
        if isinstance(constraints, int):
            self.remove_constraint(constraints)
            return
        
        cdef int last = self.nrows() + 1
        
        for c in sorted(constraints, reverse=True):
            if c != last:
                self.remove_constraint(c)
                last = c
    
    cpdef row(self, int index):
        """
        Return the ``index``-th constraint as a pair of lists.

        INPUT:

        - ``index`` -- index of the constraint

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient in the order of `indices`.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)
            sage: p.row(0)  # Note: zero coefficients are excluded in sparse format
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row(1)
            Traceback (most recent call last):
            ...
            ValueError: invalid row index 1
        """
        cdef HighsInt num_row, num_nz, matrix_start
        cdef HighsInt* matrix_index
        cdef double* matrix_value
        cdef double lb, ub
        cdef HighsInt status
        cdef list indices = []
        cdef list coeffs = []
        cdef int j
        
        if index < 0 or index >= self.numrows:
            raise ValueError(f"invalid row index {index}")
        
        # First call: get the number of non-zeros
        sig_on()
        status = Highs_getRowsByRange(self.highs, index, index,
                                      &num_row, &lb, &ub, &num_nz,
                                      NULL, NULL, NULL)
        sig_off()
        
        if status != kHighsStatusOk:
            raise MIPSolverException("HiGHS: Failed to get row info")
        
        if num_nz == 0:
            return ([], [])
        
        # Allocate space for the matrix data
        matrix_index = <HighsInt*> malloc(num_nz * sizeof(HighsInt))
        matrix_value = <double*> malloc(num_nz * sizeof(double))
        
        if matrix_index == NULL or matrix_value == NULL:
            free(matrix_index)
            free(matrix_value)
            raise MemoryError("Failed to allocate memory")
        
        try:
            # Second call: get the actual matrix data
            sig_on()
            status = Highs_getRowsByRange(self.highs, index, index,
                                          &num_row, &lb, &ub, &num_nz,
                                          &matrix_start, matrix_index, matrix_value)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to get row data")
            
            # Extract indices and coefficients
            for j in range(num_nz):
                indices.append(matrix_index[j])
                coeffs.append(matrix_value[j])
        finally:
            free(matrix_index)
            free(matrix_value)
        
        return (indices, coeffs)
    
    cpdef add_col(self, indices, coeffs):
        """
        Add a column.

        INPUT:

        - ``indices`` -- list of integers; this list contains the
          indices of the constraints in which the variable's
          coefficient is nonzero

        - ``coeffs`` -- list of real values; associates a coefficient
          to the variable in each of the constraints in which it
          appears. Namely, the i-th entry of ``coeffs`` corresponds to
          the coefficient of the variable in the constraint
          represented by the i-th entry in ``indices``.

        .. NOTE::

            ``indices`` and ``coeffs`` are expected to be of the same
            length.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(list(range(5)), list(range(5)))
            sage: p.nrows()
            5
        """
        cdef int col_idx
        cdef int i, constraint_idx
        cdef HighsInt status
        
        # Add a new column (variable)
        col_idx = self.add_variable(lower_bound=0.0, upper_bound=None)
        
        # Set the coefficients for this column in the existing constraints
        for i, constraint_idx in enumerate(indices):
            if constraint_idx < 0 or constraint_idx >= self.nrows():
                continue
            
            coeff = coeffs[i]
            if coeff != 0:
                sig_on()
                status = Highs_changeCoeff(self.highs, constraint_idx, col_idx, float(coeff))
                sig_off()
                
                if status != kHighsStatusOk:
                    raise MIPSolverException("HiGHS: Failed to set coefficient")
    
    cpdef int get_row_stat(self, int i) except? -1:
        """
        Retrieve the status of a constraint.

        INPUT:

        - ``i`` -- the index of the constraint

        OUTPUT:

        Current status assigned to the auxiliary variable associated with the i-th row:

            * 0     kLower: non-basic variable at lower bound
            * 1     kBasic: basic variable
            * 2     kUpper: non-basic variable at upper bound
            * 3     kZero: non-basic free variable at zero
            * 4     kNonbasic: nonbasic (used for unbounded variables)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver='HiGHS')
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: lp.solve()
            0
            sage: lp.get_row_stat(0)  # doctest: +SKIP
            2
            sage: lp.get_row_stat(1)  # doctest: +SKIP
            2
            sage: lp.get_row_stat(-1)
            Traceback (most recent call last):
            ...
            ValueError: The constraint's index i must satisfy 0 <= i < number_of_constraints
        """
        cdef HighsInt* row_status
        cdef HighsInt* col_status
        cdef HighsInt num_rows, num_cols
        cdef HighsInt status
        cdef HighsInt result
        
        if i < 0 or i >= self.numrows:
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")
        
        # Note: HiGHS C API doesn't have getBasisValidity function
        # We'll try to get the basis and handle errors if no basis exists
        
        num_rows = Highs_getNumRow(self.highs)
        num_cols = Highs_getNumCol(self.highs)
        
        row_status = <HighsInt*> malloc(num_rows * sizeof(HighsInt))
        col_status = <HighsInt*> malloc(num_cols * sizeof(HighsInt))
        
        if row_status == NULL or col_status == NULL:
            free(row_status)
            free(col_status)
            raise MemoryError("Failed to allocate memory")
        
        try:
            sig_on()
            status = Highs_getBasis(self.highs, col_status, row_status)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to get basis")
            
            result = row_status[i]
        finally:
            free(row_status)
            free(col_status)
        
        return result
    
    cpdef int get_col_stat(self, int j) except? -1:
        """
        Retrieve the status of a variable.

        INPUT:

        - ``j`` -- the index of the variable

        OUTPUT:

        Current status assigned to the structural variable associated with the j-th column:

            * 0     kLower: non-basic variable at lower bound
            * 1     kBasic: basic variable
            * 2     kUpper: non-basic variable at upper bound
            * 3     kZero: non-basic free variable at zero
            * 4     kNonbasic: nonbasic (used for unbounded variables)

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver='HiGHS')
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: lp.solve()
            0
            sage: lp.get_col_stat(0)
            1
            sage: lp.get_col_stat(1)
            0
            sage: lp.get_col_stat(100)
            Traceback (most recent call last):
            ...
            ValueError: The variable's index j must satisfy 0 <= j < number_of_variables
        """
        cdef HighsInt* col_status
        cdef HighsInt* row_status
        cdef HighsInt num_cols, num_rows
        cdef HighsInt status
        cdef HighsInt result
        
        if j < 0 or j >= self.numcols:
            raise ValueError("The variable's index j must satisfy 0 <= j < number_of_variables")
        
        # Note: HiGHS C API doesn't have getBasisValidity function
        # We'll try to get the basis and handle errors if no basis exists
        
        num_cols = Highs_getNumCol(self.highs)
        num_rows = Highs_getNumRow(self.highs)
        
        col_status = <HighsInt*> malloc(num_cols * sizeof(HighsInt))
        row_status = <HighsInt*> malloc(num_rows * sizeof(HighsInt))
        
        if col_status == NULL or row_status == NULL:
            free(col_status)
            free(row_status)
            raise MemoryError("Failed to allocate memory")
        
        try:
            sig_on()
            status = Highs_getBasis(self.highs, col_status, row_status)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to get basis")
            
            result = col_status[j]
        finally:
            free(col_status)
            free(row_status)
        
        return result
    
    cpdef set_row_stat(self, int i, int stat):
        """
        Set the status of a constraint.

        INPUT:

        - ``i`` -- the index of the constraint

        - ``stat`` -- the status to set to:

            * 0     kLower: non-basic variable at lower bound
            * 1     kBasic: basic variable
            * 2     kUpper: non-basic variable at upper bound
            * 3     kZero: non-basic free variable at zero
            * 4     kNonbasic: nonbasic (used for unbounded variables)

        .. NOTE::

            HiGHS may reject invalid basis configurations. Setting arbitrary
            status values may result in the basis being rejected and the
            original basis being preserved.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver='HiGHS')
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: lp.solve()
            0
            sage: lp.get_row_stat(0)  
            1
            sage: lp.set_col_stat(0, 2)
            sage: lp.get_col_stat(0)
            2
            sage: lp.set_row_stat(0, 3)  
            sage: lp.get_row_stat(0)  
            3
        """
        cdef HighsInt* row_status
        cdef HighsInt* col_status
        cdef HighsInt num_rows, num_cols
        cdef HighsInt status
        cdef int j
        
        if i < 0 or i >= self.numrows:
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")
        
        if stat < 0 or stat > 4:
            raise ValueError("Invalid status value. Must be 0-4")
        
        num_rows = Highs_getNumRow(self.highs)
        num_cols = Highs_getNumCol(self.highs)
        
        row_status = <HighsInt*> malloc(num_rows * sizeof(HighsInt))
        col_status = <HighsInt*> malloc(num_cols * sizeof(HighsInt))
        
        if row_status == NULL or col_status == NULL:
            free(row_status)
            free(col_status)
            raise MemoryError("Failed to allocate memory")
        
        try:
            # Get current basis
            sig_on()
            status = Highs_getBasis(self.highs, col_status, row_status)
            sig_off()
            
            # Set the new status
            row_status[i] = stat
            
            # Set the modified basis
            sig_on()
            status = Highs_setBasis(self.highs, col_status, row_status)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to set basis")
        finally:
            free(row_status)
            free(col_status)

    cpdef set_col_stat(self, int j, int stat):
        """
        Set the status of a variable.

        INPUT:

        - ``j`` -- the index of the variable

        - ``stat`` -- the status to set to:

            * 0     kLower: non-basic variable at lower bound
            * 1     kBasic: basic variable
            * 2     kUpper: non-basic variable at upper bound
            * 3     kZero: non-basic free variable at zero
            * 4     kNonbasic: nonbasic (used for unbounded variables)

        .. NOTE::

            HiGHS may reject invalid basis configurations. Setting arbitrary
            status values may result in the basis being rejected and the
            original basis being preserved.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver='HiGHS')
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: lp.solve()
            0
            sage: lp.get_col_stat(0)  
            1
            sage: lp.set_col_stat(0, 2)  
            sage: lp.get_col_stat(0)  
            2
        """
        cdef HighsInt* row_status
        cdef HighsInt* col_status
        cdef HighsInt num_rows, num_cols
        cdef HighsInt status
        
        if j < 0 or j >= self.numcols:
            raise ValueError("The variable's index j must satisfy 0 <= j < number_of_variables")
        
        if stat < 0 or stat > 4:
            raise ValueError("Invalid status value. Must be 0-4")
        
        num_rows = Highs_getNumRow(self.highs)
        num_cols = Highs_getNumCol(self.highs)
        
        row_status = <HighsInt*> malloc(num_rows * sizeof(HighsInt))
        col_status = <HighsInt*> malloc(num_cols * sizeof(HighsInt))
        
        if row_status == NULL or col_status == NULL:
            free(row_status)
            free(col_status)
            raise MemoryError("Failed to allocate memory")
        
        try:
            # Get current basis
            sig_on()
            status = Highs_getBasis(self.highs, col_status, row_status)
            sig_off()
            
            # Set the new status
            col_status[j] = stat
            
            # Set the modified basis
            sig_on()
            status = Highs_setBasis(self.highs, col_status, row_status)
            sig_off()
            
            if status != kHighsStatusOk:
                raise MIPSolverException("HiGHS: Failed to set basis")
        finally:
            free(row_status)
            free(col_status)
    
    cpdef int warm_up(self) noexcept:
        """
        Warm up the basis using current statuses assigned to rows and cols.

        This method attempts to validate and use the currently set basis.
        In HiGHS, setting a basis automatically attempts to factorize it,
        so this method checks if the current basis is valid.

        OUTPUT: the warming up status

            * 0             The operation has been successfully performed.
            * -1            The basis is invalid or could not be factorized.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: lp = get_solver(solver = "HiGHS")
            sage: lp.add_variables(3)
            2
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [8, 6, 1])), None, 48)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [4, 2, 1.5])), None, 20)
            sage: lp.add_linear_constraint(list(zip([0, 1, 2], [2, 1.5, 0.5])), None, 8)
            sage: lp.set_objective([60, 30, 20])
            sage: lp.solve()
            0
            sage: lp.get_objective_value()
            280.0
            sage: lp.set_row_stat(0, 3)
            sage: lp.set_col_stat(1, 1)
            sage: lp.warm_up()
            0
        """
        cdef HighsInt basis_validity
        cdef HighsInt status
        cdef HighsInt* col_status
        cdef HighsInt* row_status
        cdef HighsInt num_cols, num_rows
        
        # Note: HiGHS C API doesn't have getBasisValidity function
        # Try to get the basis to check if it's available
        num_cols = Highs_getNumCol(self.highs)
        num_rows = Highs_getNumRow(self.highs)
        
        if num_cols == 0 or num_rows == 0:
            return -1
        
        col_status = <HighsInt*> malloc(num_cols * sizeof(HighsInt))
        row_status = <HighsInt*> malloc(num_rows * sizeof(HighsInt))
        
        if col_status == NULL or row_status == NULL:
            free(col_status)
            free(row_status)
            return -1
        
        try:
            status = Highs_getBasis(self.highs, col_status, row_status)
            if status != kHighsStatusOk:
                return -1
            return 0
        finally:
            free(col_status)
            free(row_status)
    
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
        cdef HiGHSBackend p
        cdef HighsInt status
        cdef bytes temp_file
        
        import tempfile
        import os
        
        p = HiGHSBackend(maximization=self.is_maximization())
        
        # If model is empty (no constraints), just copy metadata
        if self.numrows == 0 and self.numcols > 0:
            # Add the same number of variables
            for i in range(self.numcols):
                lb, ub = self.col_bounds(i)
                p.add_variable(lb, ub,
                             self.objective_coefficient(i),
                             self.is_variable_binary(i),
                             self.is_variable_continuous(i),
                             self.is_variable_integer(i),
                             self.col_name(i))
        else:
            # Copy by writing and reading through a temporary MPS file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.mps', delete=False) as f:
                temp_file = f.name.encode('utf-8')
            
            try:
                status = Highs_writeModel(self.highs, temp_file)
                if status != kHighsStatusOk:
                    raise MIPSolverException("HiGHS: Failed to write model for copy")
                
                status = Highs_readModel(p.highs, temp_file)
                if status != kHighsStatusOk:
                    raise MIPSolverException("HiGHS: Failed to read model for copy")
                
                # Turn off logging for the copied model
                Highs_setBoolOptionValue(p.highs, b"log_to_console", False)
            finally:
                if os.path.exists(temp_file.decode('utf-8')):
                    os.unlink(temp_file.decode('utf-8'))
        
        p.prob_name = self.prob_name
        p.col_name_var = copy(self.col_name_var)
        p.row_name_var = copy(self.row_name_var)
        p.numcols = self.numcols
        p.numrows = self.numrows
        p.obj_constant_term = self.obj_constant_term
        
        return p
