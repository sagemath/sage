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
        self.row_data_cache = {}
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

        - ``level`` -- integer; from 0 (no verbosity) to 1

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.set_verbosity(0)            
        """
        # HiGHS uses different verbosity levels
        if int(level) == 0:
            self.highs_model.setOptionValue("log_to_console", False)
        elif int(level) == 1:
            self.highs_model.setOptionValue("log_to_console", True)
        else:
            raise ValueError("HiGHS supports verbosity levels 0 and 1 only.")
    
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
        
        # Cache the row data for later retrieval
        self.row_data_cache[row_idx] = (indices, values)
        
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
        if lower_bound is None and upper_bound is None:
            raise ValueError("At least one of 'upper_bound' or 'lower_bound' must be set.")
        
        import highspy
        import numpy as np
        
        # Prepare bounds
        lb = float(lower_bound) if lower_bound is not None else -highspy.kHighsInf
        ub = float(upper_bound) if upper_bound is not None else highspy.kHighsInf
        
        # Add empty constraints with the specified bounds
        for i in range(number):
            # Add empty row (no coefficients)
            self.highs_model.addRow(lb, ub, 0, 0, 0)
            
            # Store name if provided
            if names is not None and i < len(names):
                name = names[i]
                if name is not None:
                    self.row_name_var[name] = self.numrows
            
            self.numrows += 1
    
    cpdef remove_constraint(self, int i):
        r"""
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
        if i < 0 or i >= self.nrows():
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")

        import numpy as np
        # HiGHS deleteRows expects num_rows and an array of row indices
        rows_to_delete = np.array([i], dtype=np.int32)
        self.highs_model.deleteRows(1, rows_to_delete)
        
        # Update the number of rows
        self.numrows -= 1
        
        # Update the row_data_cache: shift all rows after i down by one
        new_cache = {}
        for row_idx, data in self.row_data_cache.items():
            if row_idx < i:
                new_cache[row_idx] = data
            elif row_idx > i:
                new_cache[row_idx - 1] = data
            # row_idx == i is removed, don't add to new_cache
        self.row_data_cache = new_cache
        
        # Update the row_name_var mapping
        names_to_update = {}
        names_to_remove = []
        for name, row_idx in self.row_name_var.items():
            if row_idx < i:
                names_to_update[name] = row_idx
            elif row_idx > i:
                names_to_update[name] = row_idx - 1
            else:  # row_idx == i
                names_to_remove.append(name)
        
        for name in names_to_remove:
            del self.row_name_var[name]
        self.row_name_var.update(names_to_update)

    cpdef remove_constraints(self, constraints):
        r"""
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
            sage: backend.best_known_objective_bound()
            48.0
        """
        import highspy
        info = self.highs_model.getInfo()
        lp = self.highs_model.getLp()
        
        # Check if this is a MIP (has integer variables)
        has_integer_vars = False
        if hasattr(lp, 'integrality_') and len(lp.integrality_) > 0:
            for var_type in lp.integrality_:
                if var_type == highspy.HighsVarType.kInteger:
                    has_integer_vars = True
                    break
        
        if has_integer_vars:
            # For MIP problems, use mip_dual_bound
            return info.mip_dual_bound + self.obj_constant_term
        else:
            # For LP problems, the bound equals the objective value at optimum
            return info.objective_function_value + self.obj_constant_term
    
    cpdef get_relative_objective_gap(self):
        r"""
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
        import highspy
        info = self.highs_model.getInfo()
        lp = self.highs_model.getLp()
        
        # Check if this is a MIP (has integer variables)
        has_integer_vars = False
        if hasattr(lp, 'integrality_') and len(lp.integrality_) > 0:
            for var_type in lp.integrality_:
                if var_type == highspy.HighsVarType.kInteger:
                    has_integer_vars = True
                    break
        
        if has_integer_vars:
            # For MIP problems, HiGHS provides mip_gap
            return info.mip_gap
        else:
            # For LP problems, the gap is 0 at optimum
            return 0.0
    
    cpdef get_row_prim(self, int i):
        r"""
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
            ValueError: invalid row index 2
        """
        if i < 0 or i >= self.nrows():
            raise ValueError(f"invalid row index {i}")
        
        solution = self.highs_model.getSolution()
        return solution.row_value[i]
    
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
        p.row_data_cache = copy(self.row_data_cache)
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
        if index < 0 or index >= self.ncols():
            raise ValueError(f"invalid variable index {index}")
        
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
        if index < 0 or index >= self.ncols():
            raise ValueError(f"invalid variable index {index}")
        
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

        TESTS:

        We check the behavior for an invalid index::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.is_variable_integer(2)
            False
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

        TESTS:

        We check the behavior for an invalid index::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.is_variable_binary(2)
            False
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

        TESTS:

        We check the behavior for an invalid index::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.is_variable_continuous(2)
            False
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

        TESTS:

        We sanity check the input that will be passed to HiGHS::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.set_variable_type(2, 0)
            Traceback (most recent call last):
            ...
            ValueError: invalid variable index 2
        """
        if variable < 0 or variable >= self.ncols():
            raise ValueError(f"invalid variable index {variable}")
        
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
        # Add a new column (variable)
        cdef int col_idx = self.add_variable(lower_bound=0.0, upper_bound=None)
        
        # Set the coefficients for this column in the existing constraints
        # HiGHS uses a different API - we need to update the constraint matrix
        # For each constraint index in indices, we need to add the coefficient
        # This is done by getting the constraint and updating it
        import highspy
        
        # Get current model info
        lp = self.highs_model.getLp()
        
        # For each constraint where this variable appears
        for i, constraint_idx in enumerate(indices):
            if constraint_idx < 0 or constraint_idx >= self.nrows():
                continue
            
            coeff = coeffs[i]
            if coeff != 0:
                # Add this coefficient to the constraint matrix
                # HiGHS stores the matrix, so we need to add an entry
                self.highs_model.changeCoeff(constraint_idx, col_idx, float(coeff))
                
                # Update the cached row data
                if constraint_idx in self.row_data_cache:
                    cached_indices, cached_coeffs = self.row_data_cache[constraint_idx]
                    # Add the new column to the cached data
                    cached_indices = list(cached_indices) + [col_idx]
                    cached_coeffs = list(cached_coeffs) + [float(coeff)]
                    self.row_data_cache[constraint_idx] = (cached_indices, cached_coeffs)
                else:
                    # Create new cache entry
                    self.row_data_cache[constraint_idx] = ([col_idx], [float(coeff)])
    
    cpdef row(self, int index):
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

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='HiGHS')
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(list(zip(range(5), range(5))), 2, 2)
            sage: p.row(0)
            ([0, 1, 2, 3, 4], [0.0, 1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)
        """
        if index < 0 or index >= self.nrows():
            raise ValueError(f"invalid row index {index}")
        
        # Return cached row data if available
        if index in self.row_data_cache:
            return self.row_data_cache[index]
        
        # Otherwise return empty (constraint was added without coefficients)
        return ([], [])
    
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
            sage: p.row(0)
            ([0, 1, 2, 3, 4], [0.0, 1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)
        """
        if index < 0 or index >= self.nrows():
            raise ValueError(f"invalid row index {index}")
        
        import highspy
        
        # Get the LP model
        lp = self.highs_model.getLp()
        
        # Get bounds for the constraint
        lower_bound = lp.row_lower_[index]
        upper_bound = lp.row_upper_[index]
        
        # HiGHS uses infinity constants
        inf = float('inf')
        
        # Convert infinities to None
        lower = lower_bound if lower_bound != -inf and lower_bound > -1e20 else None
        upper = upper_bound if upper_bound != inf and upper_bound < 1e20 else None
        
        return (lower, upper)
    
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
        if index < 0 or index >= self.nrows():
            raise ValueError(f"invalid row index {index}")
        
        # Look up the name in our stored mapping
        for name, idx in self.row_name_var.items():
            if idx == index:
                return name
        
        return None
    
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
        - ``solver`` -- solver to use: "choose", "simplex", "ipm", or "pdlp"
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
            sage: p.solver_parameter("threads", 4)
            sage: p.solver_parameter("threads")
            4
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
        if value is None:
            # Get the current value - getOptionValue returns (status, value)
            status, current_value = self.highs_model.getOptionValue(name)
            return current_value
        else:
            # Set the value
            # Convert Sage types to native Python types to avoid issues with highspy
            # highspy doesn't recognize Sage's Integer/RealLiteral and may misinterpret them
            if isinstance(value, (int, float, bool, str)):
                # Already a native Python type - use as-is
                py_value = value
            elif hasattr(value, 'is_integral') and value.is_integral():
                # Sage Integer or other integral type
                py_value = int(value)
            elif hasattr(value, '__float__'):
                # Sage RealLiteral or other float-like type
                py_value = float(value)
            else:
                # Other type - try to pass as-is
                py_value = value
            self.highs_model.setOptionValue(name, py_value)
    
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
        if i < 0 or i >= self.nrows():
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")
        
        cdef object basis = self.highs_model.getBasis()
        if not basis.valid:
            raise ValueError("No valid basis available. Solve the problem first.")
        
        return int(basis.row_status[i])

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
        if j < 0 or j >= self.ncols():
            raise ValueError("The variable's index j must satisfy 0 <= j < number_of_variables")
        
        cdef object basis = self.highs_model.getBasis()
        if not basis.valid:
            raise ValueError("No valid basis available. Solve the problem first.")
        
        return int(basis.col_status[j])

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
        if i < 0 or i >= self.nrows():
            raise ValueError("The constraint's index i must satisfy 0 <= i < number_of_constraints")
        
        if stat < 0 or stat > 4:
            raise ValueError("Invalid status value. Must be 0-4 (kLower, kBasic, kUpper, kZero, kNonbasic)")
        
        import highspy
        
        cdef object old_basis = self.highs_model.getBasis()
        
        # Create a fresh HighsBasis object to avoid issues with modifying retrieved basis
        cdef object new_basis = highspy.HighsBasis()
        new_basis.valid = True
        
        # Build new status lists - must create the complete list first before assigning
        # because modifying individual elements after assignment doesn't work with highspy
        cdef list new_row_status
        cdef list new_col_status
        
        if old_basis.valid:
            # Copy existing statuses and modify the specific row
            new_col_status = [s for s in old_basis.col_status]
            new_row_status = [highspy.HighsBasisStatus(s) if j != i else highspy.HighsBasisStatus(stat)
                             for j, s in enumerate(old_basis.row_status)]
        else:
            # If no valid basis exists, create default values
            new_col_status = [highspy.HighsBasisStatus.kLower for _ in range(self.ncols())]
            new_row_status = [highspy.HighsBasisStatus.kBasic if j != i else highspy.HighsBasisStatus(stat)
                             for j in range(self.nrows())]
        
        # Assign the complete lists to the new basis
        new_basis.col_status = new_col_status
        new_basis.row_status = new_row_status
        
        # Set the new basis
        self.highs_model.setBasis(new_basis)

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
        if j < 0 or j >= self.ncols():
            raise ValueError("The variable's index j must satisfy 0 <= j < number_of_variables")
        
        if stat < 0 or stat > 4:
            raise ValueError("Invalid status value. Must be 0-4 (kLower, kBasic, kUpper, kZero, kNonbasic)")
        
        import highspy
        
        cdef object old_basis = self.highs_model.getBasis()
        
        # Create a fresh HighsBasis object to avoid issues with modifying retrieved basis
        cdef object new_basis = highspy.HighsBasis()
        new_basis.valid = True
        
        # Build new status lists - must create the complete list first before assigning
        # because modifying individual elements after assignment doesn't work with highspy
        cdef list new_col_status
        cdef list new_row_status
        
        if old_basis.valid:
            # Copy existing statuses and modify the specific column
            new_col_status = [highspy.HighsBasisStatus(s) if k != j else highspy.HighsBasisStatus(stat)
                             for k, s in enumerate(old_basis.col_status)]
            new_row_status = [s for s in old_basis.row_status]
        else:
            # If no valid basis exists, create default values
            new_col_status = [highspy.HighsBasisStatus.kLower if k != j else highspy.HighsBasisStatus(stat)
                             for k in range(self.ncols())]
            new_row_status = [highspy.HighsBasisStatus.kBasic for _ in range(self.nrows())]
        
        # Assign the complete lists to the new basis
        new_basis.col_status = new_col_status
        new_basis.row_status = new_row_status
        
        # Set the new basis
        self.highs_model.setBasis(new_basis)
    
    cpdef int warm_up(self) noexcept:
        r"""
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
        import highspy
        
        # Get the current basis
        cdef object basis = self.highs_model.getBasis()
        
        # Check if the basis is valid
        if not basis.valid:
            return -1
        
        # In HiGHS, setting a basis already validates it
        # So if we have a valid basis, we're good
        return 0
    
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
