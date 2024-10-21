# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.geometry.polyhedron
r"""
Generic Backend for SDP solvers

This class only lists the methods that should be defined by any
interface with a SDP Solver. All these methods immediately raise
:exc:`NotImplementedError` exceptions when called, and are obviously
meant to be replaced by the solver-specific method. This file can also
be used as a template to create a new interface : one would only need
to replace the occurrences of ``"Nonexistent_SDP_solver"`` by the
solver's name, and replace ``GenericSDPBackend`` by
``SolverName(GenericSDPBackend)`` so that the new solver extends this
class.

AUTHORS:

- Ingolfur Edvardsson (2014-07): initial implementation
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

cdef class GenericSDPBackend:

    cpdef base_ring(self):
        """
        The base ring.

        TESTS::

            sage: from sage.numerical.backends.generic_sdp_backend import GenericSDPBackend
            sage: GenericSDPBackend().base_ring()
            Real Double Field
        """
        from sage.rings.real_double import RDF
        return RDF

    cpdef zero(self):
        """
        Zero of the base ring.

        TESTS::

            sage: from sage.numerical.backends.generic_sdp_backend import GenericSDPBackend
            sage: GenericSDPBackend().zero()
            0.0
        """
        return self.base_ring().zero()

    cpdef int add_variable(self, obj=0.0, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

        INPUT:

        - ``obj`` -- (optional) coefficient of this variable in the objective
          function (default: 0.0)

        - ``name`` -- an optional name for the newly added variable (default:
          ``None``)

        OUTPUT: the index of the newly created variable

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable(name='x', obj=1.0)
            3
            sage: p.col_name(3)
            'x'
            sage: p.objective_coefficient(3)
            1.0
        """
        raise NotImplementedError()

    cpdef int add_variables(self, int n, names=None) except -1:
        """
        Add ``n`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        INPUT:

        - ``n`` -- the number of new variables (must be > 0)

        - ``obj`` -- coefficient of all variables in the objective function (default: 0.0)

        - ``names`` -- list of names (default: ``None``)

        OUTPUT: the index of the variable created last

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.ncols()
            0
            sage: p.add_variables(5)
            4
            sage: p.ncols()
            5
            sage: p.add_variables(2, lower_bound=-2.0, integer=True, names=['a','b'])
            6
        """
        raise NotImplementedError()

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` -- integer:

          * `+1` => Maximization
          * `-1` => Minimization

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        raise NotImplementedError()

    cpdef  objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective
        function

        INPUT:

        - ``variable`` -- integer; the variable's id

        - ``coeff`` -- double; its coefficient

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.add_variable()
            1
            sage: p.objective_coefficient(0)
            0.0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2.0
        """
        raise NotImplementedError()

    cpdef  set_objective(self, list coeff, d=0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- list of real values, whose i-th element is the
          coefficient of the i-th variable in the objective function

        - ``d`` -- double; the constant term in the linear function (set to `0`
          by default)

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.add_variables(5)
            5
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: [p.objective_coefficient(x) for x in range(5)]
            [1.0, 1.0, 2.0, 1.0, 3.0]

        Constants in the objective function are respected.
        """
        raise NotImplementedError()

    cpdef add_linear_constraint(self, coefficients, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` an iterable with ``(c,v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (real
          value).

        - ``lower_bound`` -- a lower bound, either a real value or ``None``

        - ``upper_bound`` -- an upper bound, either a real value or ``None``

        - ``name`` -- an optional name for this row (default: ``None``)

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2.0, 2.0)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])                   # optional - Nonexistent_LP_solver
            sage: p.row_bounds(0)
            (2.0, 2.0)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1.0, 1.0, name='foo')
            sage: p.row_name(-1)
            "foo"
        """
        raise NotImplementedError()

    cpdef add_linear_constraints(self, int number, names=None):
        """
        Add constraints.

        INPUT:

        - ``number`` -- integer; the number of constraints to add

        - ``lower_bound`` -- a lower bound, either a real value or ``None``

        - ``upper_bound`` -- an upper bound, either a real value or ``None``

        - ``names`` -- an optional list of names (default: ``None``)

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.add_variables(5)
            5
            sage: p.add_linear_constraints(5, None, 2)
            sage: p.row(4)
            ([], [])
            sage: p.row_bounds(4)
            (None, 2.0)
        """
        raise NotImplementedError()

    cpdef int solve(self) except -1:
        """
        Solve the problem.

        .. NOTE::

            This method raises :class:`SDPSolverException` exceptions when
            the solution cannot be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.solve()
            0
            sage: p.objective_coefficient(0,1)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            SDPSolverException: ...
        """
        raise NotImplementedError()

    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.add_variables(2)
            2
            sage: p.add_linear_constraint([(0,1), (1,2)], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            7.5
            sage: p.get_variable_value(0)
            0.0
            sage: p.get_variable_value(1)
            1.5
        """
        raise NotImplementedError()

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.add_variables(2)
            2
            sage: p.add_linear_constraint([(0,1), (1, 2)], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            7.5
            sage: p.get_variable_value(0)
            0.0
            sage: p.get_variable_value(1)
            1.5
        """

        raise NotImplementedError()

    cpdef int ncols(self) noexcept:
        """
        Return the number of columns/variables.

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.ncols()
            0
            sage: p.add_variables(2)
            2
            sage: p.ncols()
            2
        """

        raise NotImplementedError()

    cpdef int nrows(self) noexcept:
        """
        Return the number of rows/constraints.

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(2, 2.0, None)
            sage: p.nrows()
            2
        """
        raise NotImplementedError()

    cpdef bint is_maximization(self) noexcept:
        """
        Test whether the problem is a maximization

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        raise NotImplementedError()

    cpdef problem_name(self, name=None):
        """
        Return or define the problem's name.

        INPUT:

        - ``name`` -- string; the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.problem_name("There once was a french fry")
            sage: print(p.problem_name())
            There once was a french fry
        """

        raise NotImplementedError()

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

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.add_variables(5)
            5
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)
        """
        raise NotImplementedError()

    cpdef row_name(self, int index):
        """
        Return the ``index``-th row name.

        INPUT:

        - ``index`` -- integer; the row's id

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.add_linear_constraints(1, 2, None, name="Empty constraint 1")
            sage: p.row_name(0)
            'Empty constraint 1'
        """
        raise NotImplementedError()

    cpdef col_name(self, int index):
        """
        Return the ``index``-th col name.

        INPUT:

        - ``index`` -- integer; the col's id

        - ``name`` -- (``char *``) its name; when set to ``NULL``
          (default), the method returns the current name

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.add_variable(name="I am a variable")
            1
            sage: p.col_name(0)
            'I am a variable'
        """
        raise NotImplementedError()

    cpdef dual_variable(self, int i, sparse=False):
        """
        The `i`-th dual variable.

        Available after ``self.solve()`` is called, otherwise the result is undefined

        - ``index`` -- integer; the constraint's id

        OUTPUT: the matrix of the `i`-th dual variable

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: p = SemidefiniteProgram(maximization=False, solver="Nonexistent_LP_solver")
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
            sage: p.solve()
            -3.0
            sage: B = p.get_backend()
            sage: x = p.get_values(x).values()
            sage: -(a3*B.dual_variable(0)).trace()-(b3*B.dual_variable(1)).trace()
            -3.0
            sage: g = sum((B.slack(j)*B.dual_variable(j)).trace() for j in range(2)); g
            0.0

        TESTS::

            sage: B.dual_variable(7)  # optional - Nonexistent_LP_solver
            ...
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: abs(g - B._get_answer()['gap'])  # optional - Nonexistent_LP_solver # tol 1e-22
            0.0
        """
        raise NotImplementedError()

    cpdef slack(self, int i, sparse=False):
        """
        Slack of the `i`-th constraint.

        Available after ``self.solve()`` is called, otherwise the result is
        undefined.

        - ``index`` -- integer; the constraint's id

        OUTPUT: the matrix of the slack of the `i`-th constraint

        EXAMPLES::

            sage: # optional - nonexistent_lp_solver
            sage: p = SemidefiniteProgram(maximization=False, solver="Nonexistent_LP_solver")
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
            sage: p.solve()
            -3.0
            sage: B = p.get_backend()
            sage: B1 = B.slack(1); B1
            [0.0 0.0]
            [0.0 0.0]
            sage: B1.is_positive_definite()
            True
            sage: x = p.get_values(x).values()
            sage: x[0]*b1 + x[1]*b2 - b3 + B1
            [0.0 0.0]
            [0.0 0.0]

        TESTS::

            sage: B.slack(7)  # optional - Nonexistent_LP_solver
            ...
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        raise NotImplementedError()

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

            sage: # optional - nonexistent_lp_solver
            sage: from sage.numerical.backends.generic_sdp_backend import get_solver
            sage: p = get_solver(solver="Nonexistent_LP_solver")
            sage: p.solver_parameter("timelimit")
            sage: p.solver_parameter("timelimit", 60)
            sage: p.solver_parameter("timelimit")
        """
        raise NotImplementedError()


default_solver = None


def default_sdp_solver(solver=None):
    """
    Return/set the default SDP solver used by Sage.

    INPUT:

    - ``solver`` -- one of the following:

      - the string ``'CVXOPT'``, to make the use of the CVXOPT solver
        (see the `CVXOPT <http://cvxopt.org/>`_ web site) the default;

      - a subclass of
        :class:`sage.numerical.backends.generic_sdp_backend.GenericSDPBackend`,
        to make it the default; or

      - ``None`` -- (default) in which case the current default solver
        (a string or a class) is returned

    OUTPUT:

    This function returns the current default solver (a string or a
    class) if ``solver = None`` (default). Otherwise, it sets the
    default solver to the one given. If this solver does not exist, or
    is not available, a :exc:`ValueError` exception is raised.

    EXAMPLES::

        sage: former_solver = default_sdp_solver()
        sage: default_sdp_solver("Cvxopt")
        sage: default_sdp_solver()
        'Cvxopt'
        sage: default_sdp_solver("Yeahhhhhhhhhhh")
        Traceback (most recent call last):
        ...
        ValueError: 'solver' should be set to ...
        sage: default_sdp_solver(former_solver)
        sage: from sage.numerical.backends.generic_sdp_backend import GenericSDPBackend
        sage: class my_sdp_solver(GenericSDPBackend): pass
        sage: default_sdp_solver(my_sdp_solver)
        sage: default_sdp_solver() is my_sdp_solver
        True
    """
    global default_solver

    if solver is None:

        if default_solver is not None:
            return default_solver

        else:
            for s in ["Cvxopt"]:
                try:
                    default_sdp_solver(s)
                    return s
                except ValueError:
                    pass

            from warnings import warn
            warn("default_sdp_solver set to 'Matrix' (MatrixSDPBackend), which can construct but not solve problems. Install cvxopt for actual solver functionality")
            default_sdp_solver("Matrix")

    if callable(solver):
        default_solver = solver
        return

    solver = solver.capitalize()

    if solver == "Cvxopt":
        try:
            from sage.numerical.backends.cvxopt_sdp_backend import CVXOPTSDPBackend
            default_solver = solver
        except ImportError:
            raise ValueError("CVXOPT is not available. Please refer to the documentation to install it.")
    elif solver == "Matrix":
        default_solver = solver

    else:
        raise ValueError("'solver' should be set to 'CVXOPT', 'Matrix', a class, or None.")


cpdef GenericSDPBackend get_solver(solver=None, base_ring=None):
    """
    Return a solver according to the given preferences.

    INPUT:

    - ``solver`` -- one of the following:

      - the string ``'CVXOPT'``, designating the use of the CVXOPT solver
        (see the `CVXOPT <http://cvxopt.org/>`_ web site);

      - a subclass of
        :class:`sage.numerical.backends.generic_sdp_backend.GenericSDPBackend`

      - ``None`` -- (default) in which case the default solver is used (see
        :func:`default_sdp_solver`)

    .. SEEALSO::

        - :func:`default_sdp_solver` -- returns/sets the default SDP solver

    EXAMPLES::

        sage: from sage.numerical.backends.generic_sdp_backend import get_solver
        sage: p = get_solver()

    Passing a class::

        sage: from sage.numerical.backends.generic_sdp_backend import GenericSDPBackend
        sage: class MockSDPBackend(GenericSDPBackend):
        ....:     def solve(self):
        ....:         raise RuntimeError("SDP is too slow")
        sage: P = SemidefiniteProgram(solver=MockSDPBackend)
        sage: P.solve()
        Traceback (most recent call last):
        ...
        RuntimeError: SDP is too slow
    """
    if solver is None:
        solver = default_sdp_solver()

    if callable(solver):
        return solver(base_ring=base_ring)

    solver = solver.capitalize()

    if solver == "Cvxopt":
        from sage.numerical.backends.cvxopt_sdp_backend import CVXOPTSDPBackend
        return CVXOPTSDPBackend(base_ring=base_ring)
    elif solver == "Matrix":
        from sage.numerical.backends.matrix_sdp_backend import MatrixSDPBackend
        return MatrixSDPBackend(base_ring=base_ring)
    else:
        raise ValueError("'solver' should be set to 'CVXOPT', 'Matrix', a class, or None (in which case the default one is used).")
