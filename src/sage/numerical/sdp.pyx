r"""
Semidefinite Programming

A semidefinite program (:wikipedia:`SDP <Semidefinite_programming>`)
is an optimization problem (:wikipedia:`Optimization_(mathematics)>`)
of the following form

.. MATH::

    \min \sum_{i,j=1}^n C_{ij}X_{ij} & \qquad \text{(Dual problem)}\\
    \text{Subject to:} & \sum_{i,j=1}^n A_{ijk}X_{ij} = b_k, \qquad k=1\dots m\\
    &X \succeq 0

where the `X_{ij}`, `1 \leq i,j \leq n` are `n^2` variables satisfying the symmetry
conditions `x_{ij} = x_{ji}` for all `i,j`, the `C_{ij}=C_{ji}`, `A_{ijk}=A_{kji}` and `b_k`
are real coefficients, and `X` is positive semidefinite, i.e., all the eigenvalues of `X` are nonnegative.
The closely related dual problem of this one is the following, where we denote by
`A_k` the matrix `(A_{kij})` and by `C` the matrix `(C_{ij})`,

.. MATH::

    \max \sum_k b_k x_k & \qquad \text{(Primal problem)}\\
    \text{Subject to:} & \sum_k x_k A_k \preceq C.

Here `(x_1,...,x_m)` is a vector of scalar variables.
A wide variety of problems in optimization can be formulated in one of these two standard
forms. Then, solvers are able to calculate an approximation to a solution.
Here we refer to the latter problem as primal, and to the former problem as dual.
The optimal value of the dual is always at least the
optimal value of the primal, and usually (although not always) they are equal.

For instance, suppose you want to maximize `x_1 - x_0` subject to

.. MATH::

 \left( \begin{array}{cc} 1 & 2  \\ 2 & 3  \end{array} \right) x_0 +
 \left( \begin{array}{cc} 3 & 4  \\ 4 & 5  \end{array} \right) x_1 \preceq
 \left( \begin{array}{cc} 5 & 6  \\ 6 & 7  \end{array} \right),\quad
 \left( \begin{array}{cc} 1 & 1  \\ 1 & 1  \end{array} \right) x_0 +
 \left( \begin{array}{cc} 2 & 2  \\ 2 & 2  \end{array} \right) x_1 \preceq
 \left( \begin{array}{cc} 3 & 3  \\ 3 & 3  \end{array} \right),
 \quad x_0\geq 0, x_1\geq 0.

An SDP can give you an answer to the problem above. Here is how it's done:

  #. You have to create an instance of :class:`SemidefiniteProgram`.
  #. Create a dictionary `x` of integer variables via :meth:`~SemidefiniteProgram.new_variable`,
     for example doing ``x = p.new_variable()`` if ``p`` is the name of the SDP instance.
  #. Add those two matrix inequalities as inequality constraints via
     :meth:`~SemidefiniteProgram.add_constraint`.
  #. Add another matrix inequality to specify nonnegativity of `x`.
  #. Specify the objective function via :meth:`~SemidefiniteProgram.set_objective`.
     In our case it is  `x_1 - x_0`. If it
     is a pure constraint satisfaction problem, specify it as ``None``.
  #. To check if everything is set up correctly, you can print the problem via
     :meth:`show <sage.numerical.sdp.SemidefiniteProgram.show>`.
  #. :meth:`Solve <sage.numerical.sdp.SemidefiniteProgram.solve>` it and print the solution.

The following example shows all these steps::

    sage: p = SemidefiniteProgram()
    sage: x = p.new_variable()
    sage: p.set_objective(x[1] - x[0])
    sage: a1 = matrix([[1, 2.], [2., 3.]])
    sage: a2 = matrix([[3, 4.], [4., 5.]])
    sage: a3 = matrix([[5, 6.], [6., 7.]])
    sage: b1 = matrix([[1, 1.], [1., 1.]])
    sage: b2 = matrix([[2, 2.], [2., 2.]])
    sage: b3 = matrix([[3, 3.], [3., 3.]])
    sage: c1 = matrix([[1.0, 0],[0,0]],sparse=True)
    sage: c2 = matrix([[0.0, 0],[0,1]],sparse=True)
    sage: p.add_constraint(a1*x[0] + a2*x[1] <= a3)
    sage: p.add_constraint(b1*x[0] + b2*x[1] <= b3)
    sage: p.add_constraint(c1*x[0] + c2*x[1] >= matrix.zero(2,2,sparse=True))

    sage: # needs cvxopt
    sage: p.solver_parameter("show_progress", True)
    sage: opt = p.solve()
        pcost       dcost       gap    pres   dres   k/t
    0: ...
    ...
    Optimal solution found.
    sage: print('Objective Value: {}'.format(N(opt,3)))
    Objective Value: 1.0
    sage: [N(x, 3) for x in sorted(p.get_values(x).values())]
    [3.0e-8, 1.0]
    sage: p.show()
    Maximization:
      x_0 - x_1
    Constraints:
      constraint_0: [3.0 4.0][4.0 5.0]x_0 + [1.0 2.0][2.0 3.0]x_1 <=  [5.0 6.0][6.0 7.0]
      constraint_1: [2.0 2.0][2.0 2.0]x_0 + [1.0 1.0][1.0 1.0]x_1 <=  [3.0 3.0][3.0 3.0]
      constraint_2: [ 0.0  0.0][ 0.0 -1.0]x_0 + [-1.0  0.0][ 0.0  0.0]x_1 <=  [0 0][0 0]
    Variables:
       x_0,  x_1

Most solvers, e.g. the default Sage SDP solver CVXOPT, solve simultaneously the pair
of primal and dual problems. Thus we can get the optimizer `X` of the dual problem
as follows, as diagonal blocks, one per each constraint, via :meth:`~SemidefiniteProgram.dual_variable`.
E.g.::

    sage: p.dual_variable(1)  # rel tol 2e-03                                           # needs cvxopt
    [ 85555.0 -85555.0]
    [-85555.0  85555.0]

We can see that the optimal value of the dual is equal (up to numerical noise) to `opt`.::

    sage: opt - ((p.dual_variable(0)*a3).trace()  # tol 8e-08                           # needs cvxopt
    ....:        + (p.dual_variable(1)*b3).trace())
    0.0

Dual variable blocks at optimality are orthogonal to "slack variables", that is,
matrices `C-\sum_k x_k A_k`, cf. (Primal problem) above, available via
:meth:`~SemidefiniteProgram.slack`. E.g.::

    sage: (p.slack(0)*p.dual_variable(0)).trace()  # tol 2e-07                          # needs cvxopt
    0.0


More interesting example, the :func:`Lovasz theta <sage.graphs.lovasz_theta.lovasz_theta>` of the 7-gon::

    sage: # needs sage.graphs
    sage: c = graphs.CycleGraph(7)
    sage: c2 = c.distance_graph(2).adjacency_matrix()
    sage: c3 = c.distance_graph(3).adjacency_matrix()
    sage: p.<y> = SemidefiniteProgram()
    sage: p.add_constraint((1/7)*matrix.identity(7)>=-y[0]*c2-y[1]*c3)
    sage: p.set_objective(y[0]*(c2**2).trace()+y[1]*(c3**2).trace())
    sage: x = p.solve(); x + 1                                                          # needs cvxopt
    3.31766...

Unlike in the previous example, the slack variable is very far from 0::

    sage: p.slack(0).trace()   # tol 1e-14                                              # needs cvxopt sage.graphs
    1.0

The default CVXOPT backend computes with the Real Double Field, for example::

    sage: # needs cvxopt
    sage: p = SemidefiniteProgram(solver='cvxopt')
    sage: p.base_ring()
    Real Double Field
    sage: x = p.new_variable()
    sage: 0.5 + 3/2*x[1]
    0.5 + 1.5*x_0

For representing an SDP with exact data, there is another backend::

    sage: from sage.numerical.backends.matrix_sdp_backend import MatrixSDPBackend
    sage: p = SemidefiniteProgram(solver=MatrixSDPBackend, base_ring=QQ)
    sage: p.base_ring()
    Rational Field
    sage: x = p.new_variable()
    sage: 1/2 + 3/2 * x[1]
    1/2 + 3/2*x_0


Linear Variables and Expressions
--------------------------------


To make your code more readable, you can construct
:class:`SDPVariable` objects that can be arbitrarily named and
indexed. Internally, this is then translated back to the `x_i`
variables. For example::

    sage: sdp.<a,b> = SemidefiniteProgram()
    sage: a
    SDPVariable
    sage: 5 + a[1] + 2*b[3]
    5 + x_0 + 2*x_1

Indices can be any object, not necessarily integers. Multi-indices are
also allowed::

    sage: a[4, 'string', QQ]
    x_2
    sage: a[4, 'string', QQ] - 7*b[2]
    x_2 - 7*x_3
    sage: sdp.show()
    Maximization:
    <BLANKLINE>
    Constraints:
    Variables:
      a[1],  b[3],  a[(4, 'string', Rational Field)],  b[2]

Index of functions and methods
------------------------------

Below are listed the methods of :class:`SemidefiniteProgram`. This module
also implements the :class:`SDPSolverException` exception, as well as the
:class:`SDPVariable` class.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~SemidefiniteProgram.add_constraint`            | Add a constraint to the ``SemidefiniteProgram``
    :meth:`~SemidefiniteProgram.base_ring`                 | Return the base ring
    :meth:`~SemidefiniteProgram.dual_variable`             | Return optimal dual variable block
    :meth:`~SemidefiniteProgram.get_backend`               | Return the backend instance used
    :meth:`~SemidefiniteProgram.get_values`                | Return values found by the previous call to ``solve()``
    :meth:`~SemidefiniteProgram.linear_constraints_parent` | Return the parent for all linear constraints
    :meth:`~SemidefiniteProgram.linear_function`           | Construct a new linear function
    :meth:`~SemidefiniteProgram.linear_functions_parent`   | Return the parent for all linear functions
    :meth:`~SemidefiniteProgram.new_variable`              | Return an instance of ``SDPVariable`` associated to the ``SemidefiniteProgram``
    :meth:`~SemidefiniteProgram.number_of_constraints`     | Return the number of constraints assigned so far
    :meth:`~SemidefiniteProgram.number_of_variables`       | Return the number of variables used so far
    :meth:`~SemidefiniteProgram.set_objective`             | Set the objective of the ``SemidefiniteProgram``
    :meth:`~SemidefiniteProgram.set_problem_name`          | Set the name of the ``SemidefiniteProgram``
    :meth:`~SemidefiniteProgram.slack`                     | Return the slack variable block at the optimum
    :meth:`~SemidefiniteProgram.show`                      | Display the ``SemidefiniteProgram`` in a human-readable way
    :meth:`~SemidefiniteProgram.solve`                     | Solve the ``SemidefiniteProgram``
    :meth:`~SemidefiniteProgram.solver_parameter`          | Return or define a solver parameter
    :meth:`~SemidefiniteProgram.sum`                       | Efficiently compute the sum of a sequence of LinearFunction elements

AUTHORS:

- Ingolfur Edvardsson (2014/08): added extension for exact computation

- Dima Pasechnik      (2014-)    : supervision, minor fixes, duality
"""
# ****************************************************************************
#       Copyright (C) 2014 Ingolfur Edvardsson <ingolfured@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.numerical.linear_functions import LinearFunction, LinearConstraint
from sage.matrix.constructor import matrix
from sage.structure.element import Matrix


cdef class SemidefiniteProgram(SageObject):
    r"""
    The ``SemidefiniteProgram`` class is the link between Sage, semidefinite
    programming (SDP) and semidefinite programming solvers.

    A Semidefinite Programming (SDP) consists of variables, linear
    constraints on these variables, and an objective function which is to be
    maximised or minimised under these constraints.

    See the :wikipedia:`Semidefinite_programming` for further information on semidefinite
    programming, and the :mod:`SDP module <sage.numerical.sdp>` for its use in
    Sage.

    INPUT:

    - ``solver`` -- selects a solver:

      - CVXOPT (``solver="CVXOPT"``). See the `CVXOPT <http://www.cvxopt.org/>`_
        website.

      - If ``solver=None`` (default), the default solver is used (see
        :func:`default_sdp_solver`)

    - ``maximization``

      - When set to ``True`` (default), the ``SemidefiniteProgram``
        is defined as a maximization.

      - When set to ``False``, the ``SemidefiniteProgram`` is
        defined as a minimization.


    .. SEEALSO::

        - :func:`default_sdp_solver` -- returns/sets the default SDP solver

    EXAMPLES:

    Computation of a basic Semidefinite Program::

         sage: p = SemidefiniteProgram(maximization=False)
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
         sage: N(p.solve(), 2)                                                          # needs cvxopt
         -3.0
    """

    def __init__(self, solver=None, maximization=True,
                 names=tuple(), **kwds):
        r"""
        Constructor for the ``SemidefiniteProgram`` class.

        INPUT:

        - ``solver`` -- the following solvers should be available through this class:

          - CVXOPT (``solver="CVXOPT"``). See the `CVXOPT <http://www.cvxopt.org/>`_
              web site.

          - If ``solver=None`` (default), the default solver is used (see
           ``default_sdp_solver`` method.

        - ``maximization``

          - When set to ``True`` (default), the ``SemidefiniteProgram``
            is defined as a maximization.
          - When set to ``False``, the ``SemidefiniteProgram`` is
            defined as a minimization.

        - ``names`` -- list/tuple/iterable of string. Default names of
          the SDP variables. Used to enable the ``sdp.<x> =
          SemidefiniteProgram()`` syntax.

        - other keyword arguments are passed to the solver.

        .. SEEALSO::

        - :meth:`default_sdp_solver` -- returns/Sets the default SDP solver

        EXAMPLES::

            sage: p = SemidefiniteProgram(maximization=True)
        """
        self._first_variable_names = list(names)
        from sage.numerical.backends.generic_sdp_backend import get_solver
        self._backend = get_solver(solver=solver, **kwds)
        if not maximization:
            self._backend.set_sense(-1)

        # Associates an index to the variables
        self._variables = {}

    def linear_functions_parent(self):
        """
        Return the parent for all linear functions.

        EXAMPLES::

             sage: p = SemidefiniteProgram()
             sage: p.linear_functions_parent()
             Linear functions over Real Double Field
        """
        if self._linear_functions_parent is None:
            base_ring = self._backend.base_ring()
            from sage.numerical.linear_functions import LinearFunctionsParent
            self._linear_functions_parent = LinearFunctionsParent(base_ring)
        return self._linear_functions_parent

    def linear_constraints_parent(self):
        """
        Return the parent for all linear constraints.

        See :mod:`~sage.numerical.linear_functions` for more
        details.

        EXAMPLES::

             sage: p = SemidefiniteProgram()
             sage: p.linear_constraints_parent()
             Linear constraints over Real Double Field
        """
        if self._linear_constraints_parent is None:
            from sage.numerical.linear_functions import LinearConstraintsParent
            LF = self.linear_functions_parent()
            self._linear_constraints_parent = LinearConstraintsParent(LF)
        return self._linear_constraints_parent

    def __call__(self, x):
        """
        Construct a new linear function.

        EXAMPLES::

             sage: p = SemidefiniteProgram()
             sage: p.linear_function({0:1})
             x_0
        """
        parent = self.linear_functions_parent()
        return parent(x)

    linear_function = __call__

    def _repr_(self):
        r"""
        Return a short description of the ``SemidefiniteProgram``.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: x = p.new_variable()
            sage: a1 = matrix([[1, 2.], [2., 3.]])
            sage: a2 = matrix([[3, 4.], [4., 5.]])
            sage: a3 = matrix([[5, 6.], [6., 7.]])
            sage: b1 = matrix([[1, 1.], [1., 1.]])
            sage: b2 = matrix([[2, 2.], [2., 2.]])
            sage: b3 = matrix([[3, 3.], [3., 3.]])
            sage: p.add_constraint(a1*x[0] + a2*x[1] <= a3)
            sage: p.add_constraint(b1*x[0] + b2*x[1] <= b3)
            sage: p.add_constraint(b1*x[0] + b2*x[1] <= a1)
            sage: print(p)
            Semidefinite Program ( maximization, 2 variables, 3 constraints )
        """
        cdef GenericSDPBackend b = self._backend

        return ("Semidefinite Program " +

                ( "\"" +self._backend.problem_name()+ "\""
                  if (str(self._backend.problem_name()) != "") else "")+

                " ( " + ("maximization" if b.is_maximization() else "minimization" ) +

                ", " + str(b.ncols()) + " variables, " +
                str(b.nrows()) + " constraints )")

    def __getitem__(self, v):
        r"""
        Return the symbolic variable corresponding to the key
        from a default dictionary.

        It returns the element asked, and otherwise creates it.
        If necessary, it also creates the default dictionary.

        This method lets the user define LinearProgram without having to
        define independent dictionaries when it is not necessary for him.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: p.set_objective(p['x'] + p['z'])
            sage: p['x']
            x_0
        """

        try:
            return self._default_sdpvariable[v]
        except TypeError:
            self._default_sdpvariable = self.new_variable()
            return self._default_sdpvariable[v]

    def base_ring(self):
        """
        Return the base ring.

        OUTPUT: a ring. The coefficients that the chosen solver supports

        EXAMPLES::

            sage: p = SemidefiniteProgram(solver='cvxopt')
            sage: p.base_ring()
            Real Double Field
        """
        return self._backend.base_ring()

    def set_problem_name(self,name):
        r"""
        Set the name of the ``SemidefiniteProgram``.

        INPUT:

        - ``name`` -- string representing the name of the
          ``SemidefiniteProgram``

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: p.set_problem_name("Test program")
            sage: p
            Semidefinite Program "Test program" ( maximization, 0 variables, 0 constraints )
        """
        self._backend.problem_name(name)

    def new_variable(self, name=""):
        r"""
        Return an instance of :class:`SDPVariable` associated
        to the current instance of :class:`SemidefiniteProgram`.

        A new variable ``x`` is defined by::

            sage: p = SemidefiniteProgram()
            sage: x = p.new_variable()

        It behaves exactly as an usual dictionary would. It can use any key
        argument you may like, as ``x[5]`` or ``x["b"]``, and has methods
        ``items()`` and ``keys()``.

        INPUT:

        - ``dim`` -- integer; defines the dimension of the dictionary
          If ``x`` has dimension `2`, its fields will be of the form
          ``x[key1][key2]``. Deprecated.

        - ``name`` -- string; associates a name to the variable

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: x = p.new_variable()
            sage: a1 = matrix([[1, 2.], [2., 3.]])
            sage: p.add_constraint(a1*x[0] + a1*x[3] <= 0)
            sage: p.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              constraint_0: [1.0 2.0][2.0 3.0]x_0 + [1.0 2.0][2.0 3.0]x_1 <=  [0 0][0 0]
            Variables:
               x_0,  x_1
        """
        if not name and self._first_variable_names:
            name = self._first_variable_names.pop(0)

        return sdp_variable_parent(self,
                      name=name)

    def _first_ngens(self, n):
        """
        Construct the first `n` SDPVariables.

        This method is used for the generator syntax (see below). You
        probably should not use it for anything else.

        INPUT:

        - ``n`` -- integer; the number of variables to construct

        OUTPUT:

        A tuple of not necessarily positive :class:`SDPVariable`
        instances.

        EXAMPLES::

            sage: sdp.<a,b> = SemidefiniteProgram()
            sage: a[0] + b[2]
            x_0 + x_1
            sage: sdp.show()
            Maximization:
            <BLANKLINE>
            Constraints:
            Variables:
              a[0],  b[2]
        """
        return tuple(self.new_variable() for i in range(n))

    def gen(self, i):
        """
        Return the linear variable `x_i`.

        EXAMPLES::

            sage: sdp = SemidefiniteProgram()
            sage: sdp.gen(0)
            x_0
            sage: [sdp.gen(i) for i in range(10)]
            [x_0, x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9]
        """
        return self.linear_functions_parent().gen(i)

    cpdef int number_of_constraints(self) noexcept:
        r"""
        Return the number of constraints assigned so far.

        EXAMPLES::

            sage: p = SemidefiniteProgram(solver = "cvxopt")
            sage: x = p.new_variable()
            sage: a1 = matrix([[1, 2.], [2., 3.]])
            sage: a2 = matrix([[3, 4.], [4., 5.]])
            sage: a3 = matrix([[5, 6.], [6., 7.]])
            sage: b1 = matrix([[1, 1.], [1., 1.]])
            sage: b2 = matrix([[2, 2.], [2., 2.]])
            sage: b3 = matrix([[3, 3.], [3., 3.]])
            sage: p.add_constraint(a1*x[0] + a2*x[1] <= a3)
            sage: p.add_constraint(b1*x[0] + b2*x[1] <= b3)
            sage: p.add_constraint(b1*x[0] + a2*x[1] <= b3)
            sage: p.number_of_constraints()
            3
        """
        return self._backend.nrows()

    cpdef int number_of_variables(self) noexcept:
        r"""
        Return the number of variables used so far.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: a = matrix([[1, 2.], [2., 3.]])
            sage: p.add_constraint(a*p[0] - a*p[2] <=  2*a*p[4]  )
            sage: p.number_of_variables()
            3
        """
        return self._backend.ncols()

    def show(self):
        r"""
        Display the :class:`SemidefiniteProgram` in a human-readable way.

        EXAMPLES:

        When constraints and variables have names ::

              sage: p = SemidefiniteProgram()
              sage: x = p.new_variable(name='hihi')
              sage: a1 = matrix([[1,2],[2,3]])
              sage: a2 = matrix([[2,3],[3,4]])
              sage: a3 = matrix([[3,4],[4,5]])
              sage: p.set_objective(x[0] - x[1])
              sage: p.add_constraint(a1*x[0] + a2*x[1]<= a3)
              sage: p.show()
              Maximization:
                hihi[0] - hihi[1]
              Constraints:
                constraint_0: [1.0 2.0][2.0 3.0]hihi[0] + [2.0 3.0][3.0 4.0]hihi[1] <=  [3.0 4.0][4.0 5.0]
              Variables:
                 hihi[0],  hihi[1]
        """
        cdef int i, j
        cdef GenericSDPBackend b = self._backend

        # inv_variables associates a SDPVariable object to an id
        inv_variables = {}
        for (v, id) in self._variables.iteritems():
            inv_variables[id] = v

        # varid_name associates variables id to names
        varid_name = {}
        for 0<= i < b.ncols():
            s = b.col_name(i)
            varid_name[i] = s if s else 'x_'+str(i)

        ##### Sense and objective function
        print("Maximization:" if b.is_maximization() else "Minimization:")
        print(" ", end=" ")
        first = True
        for 0<= i< b.ncols():
            c = b.objective_coefficient(i)
            if c == 0:
                continue
            print((("+ " if (not first and c>0) else "") +
                   ("" if c == 1 else ("- " if c == -1 else str(c)+" ")) + varid_name[i]
                   ), end=" ")
            first = False
        if b.obj_constant_term > self._backend.zero():
            print("+ {}".format(b.obj_constant_term))
        elif b.obj_constant_term < self._backend.zero():
            print("- {}".format(-b.obj_constant_term))
        print("\n")

        ##### Constraints
        print("Constraints:")
        for 0<= i < b.nrows():
            indices, values = b.row(i)
            print(" ", end=" ")
            # Constraint's name
            if b.row_name(i):
                print(b.row_name(i)+":", end=" ")
            first = True
            l = sorted(zip(indices,values))
            l.reverse()
            if l[-1][0] == -1:
                last_i,last_value = l.pop()
            else:
                last_value = matrix.zero( l[0][1].dimensions()[0],l[0][1].dimensions()[1]  )
            l.reverse()
            for j, c in l:
                if c == 0:
                    continue
                print((("+ " if (not first) else "") +
                        ( str(repr(c).replace('\n',"")  )  )+varid_name[j]),
                      end=" ")
                first = False
            print(("<= "), end=" ")
            print(repr(-last_value).replace('\n',""))

        ##### Variables
        print("Variables:")
        print("  ", end=" ")
        for 0<= i < b.ncols()-1:
            print(str(varid_name[i]) + ", ", end=" ")
        print(str(varid_name[b.ncols()-1]))

    def get_values(self, *lists):
        r"""
        Return values found by the previous call to ``solve()``.

        INPUT:

        - Any instance of :class:`SDPVariable` (or one of its elements),
          or lists of them.

        OUTPUT:

        - Each instance of :class:`SDPVariable` is replaced by a dictionary
          containing the numerical values found for each
          corresponding variable in the instance.
        - Each element of an instance of a :class:`SDPVariable` is replaced
          by its corresponding numerical value.

        EXAMPLES::

            sage: p = SemidefiniteProgram(solver = "cvxopt", maximization = False)
            sage: x = p.new_variable()
            sage: p.set_objective(x[3] - x[5])
            sage: a1 = matrix([[1, 2.], [2., 3.]])
            sage: a2 = matrix([[3, 4.], [4., 5.]])
            sage: a3 = matrix([[5, 6.], [6., 7.]])
            sage: b1 = matrix([[1, 1.], [1., 1.]])
            sage: b2 = matrix([[2, 2.], [2., 2.]])
            sage: b3 = matrix([[3, 3.], [3., 3.]])
            sage: p.add_constraint(a1*x[3] + a2*x[5] <= a3)
            sage: p.add_constraint(b1*x[3] + b2*x[5] <= b3)
            sage: N(p.solve(), 3)                                                       # needs cvxopt
            -3.0

        To return  the optimal value of ``x[3]``::

            sage: N(p.get_values(x[3]),3)                                               # needs cvxopt
            -1.0

        To get a dictionary identical to ``x`` containing optimal
        values for the corresponding variables ::

            sage: x_sol = p.get_values(x)                                               # needs cvxopt
            sage: sorted(x_sol)                                                         # needs cvxopt
            [3, 5]
        """
        val = []
        for l in lists:
            if isinstance(l, SDPVariable):
                    c = {}
                    for k, v in l.items():
                        c[k] = self._backend.get_variable_value(self._variables[v])
                    val.append(c)
            elif isinstance(l, list):
                if len(l) == 1:
                    val.append([self.get_values(l[0])])
                else:
                    c = []
                    [c.append(self.get_values(ll)) for ll in l]
                    val.append(c)
            elif l in self._variables:
                #val.append(self._values[l])
                val.append(self._backend.get_variable_value(self._variables[l]))

        if len(lists) == 1:
            return val[0]
        else:
            return val

    def set_objective(self, obj):
        r"""
        Set the objective of the :class:`SemidefiniteProgram`.

        INPUT:

        - ``obj`` -- a semidefinite function to be optimized
          (can also be set to ``None`` or ``0`` when just
          looking for a feasible solution)

        EXAMPLES:

        Let's solve the following semidefinite program:

        .. MATH::

            \begin{aligned}
                \text{maximize} &\qquad x + 5y  \qquad \\
                \text{subject to} &\qquad \left( \begin{array}{cc} 1 & 2  \\ 2 & 3  \end{array} \right) x +
                \left( \begin{array}{cc} 1 & 1  \\ 1 & 1  \end{array} \right) y \preceq
                \left( \begin{array}{cc} 1 & -1  \\ -1 & 1  \end{array} \right)
            \end{aligned}

        This SDP can be solved as follows::

            sage: p = SemidefiniteProgram(maximization=True)
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + 5*x[2])
            sage: a1 = matrix([[1,2],[2,3]])
            sage: a2 = matrix([[1,1],[1,1]])
            sage: a3 = matrix([[1,-1],[-1,1]])
            sage: p.add_constraint(a1*x[1] + a2*x[2] <= a3)
            sage: N(p.solve(), digits=3)                                                # needs cvxopt
            16.2
            sage: p.set_objective(None)
            sage: _ = p.solve()                                                         # needs cvxopt
        """
        cdef list values = []

        # If the objective is None, or a constant, we want to remember
        # that the objective function has been defined ( the user did not
        # forget it ). In some SDO problems, you just want a feasible solution
        # and do not care about any function being optimal.
        cdef int i

        if obj is not None:
            f = obj.dict()
        else:
            f = {-1 : 0}
        d = f.pop(-1,self._backend.zero())

        for i in range(self._backend.ncols()):
            values.append(f.get(i,self._backend.zero()))
        self._backend.set_objective(values, d)

    def add_constraint(self, linear_function, name=None):
        r"""
        Add a constraint to the ``SemidefiniteProgram``.

        INPUT:

        - ``linear_function`` -- two different types of arguments are possible:

          - A linear function. In this case, arguments ``min`` or ``max``
            have to be specified.
          - A linear constraint of the form ``A <= B``, ``A >= B``,
            ``A <= B <= C``, ``A >= B >= C`` or ``A == B``. In this
            case, arguments ``min`` and ``max`` will be ignored.

        - ``name`` -- a name for the constraint

        EXAMPLES:

        Let's solve the following semidefinite program:

        .. MATH::

            \begin{aligned}
                \text{maximize} &\qquad x + 5y  \qquad \\
                \text{subject to} &\qquad \left( \begin{array}{cc} 1 & 2  \\ 2 & 3  \end{array} \right) x +
                \left( \begin{array}{cc} 1 & 1  \\ 1 & 1  \end{array} \right) y \preceq
                \left( \begin{array}{cc} 1 & -1  \\ -1 & 1  \end{array} \right)
            \end{aligned}

        This SDP can be solved as follows::

            sage: p = SemidefiniteProgram(maximization=True)
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + 5*x[2])
            sage: a1 = matrix([[1,2],[2,3]])
            sage: a2 = matrix([[1,1],[1,1]])
            sage: a3 = matrix([[1,-1],[-1,1]])
            sage: p.add_constraint(a1*x[1] + a2*x[2] <= a3)
            sage: N(p.solve(), digits=3)                                                # needs cvxopt
            16.2

        One can also define double-bounds or equality using the symbol
        ``>=`` or ``==``::

            sage: p = SemidefiniteProgram(maximization=True)
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + 5*x[2])
            sage: a1 = matrix([[1,2],[2,3]])
            sage: a2 = matrix([[1,1],[1,1]])
            sage: a3 = matrix([[1,-1],[-1,1]])
            sage: p.add_constraint(a3 >= a1*x[1] + a2*x[2])
            sage: N(p.solve(), digits=3)                                                # needs cvxopt
            16.2

        TESTS:

        Complex constraints::

            sage: p = SemidefiniteProgram()
            sage: b = p.new_variable()
            sage: a1 = matrix([[1,2],[2,3]])
            sage: a2 = matrix([[1,-2],[-2,4]])
            sage: p.add_constraint(a1*b[8] - a1*b[15] <= a2*b[8])
            sage: p.show()
            Maximization:
            <BLANKLINE>
            Constraints:
                constraint_0: [ 0.0  4.0][ 4.0 -1.0]x_0 + [-1.0 -2.0][-2.0 -3.0]x_1 <=  [0 0][0 0]
            Variables:
              x_0, x_1

        Empty constraint::

            sage: p = SemidefiniteProgram()
            sage: p.add_constraint(sum([]))
        """
        if linear_function is 0:
            return

        from sage.numerical.linear_tensor_constraints import LinearTensorConstraint
        from sage.numerical.linear_tensor import LinearTensor

        if isinstance(linear_function, LinearTensorConstraint) or isinstance(linear_function, LinearConstraint):
            c = linear_function
            if c.is_equation():
                self.add_constraint(c.lhs()-c.rhs(), name=name)
                self.add_constraint(-c.lhs()+c.rhs(), name=name)
            else:
                self.add_constraint(c.lhs()-c.rhs(), name=name)

        elif isinstance(linear_function, LinearFunction) or isinstance(linear_function, LinearTensor):
            l = sorted(linear_function.dict().items())
            self._backend.add_linear_constraint(l, name)

        else:
            raise ValueError('argument must be a linear function or constraint, got '+str(linear_function))

    def solve(self, objective_only=False):
        r"""
        Solve the :class:`SemidefiniteProgram`.

        INPUT:

        - ``objective_only`` -- boolean:

          - when set to ``True``, only the objective function is returned
          - when set to ``False`` (default), the optimal numerical values
            are stored (takes computational time)

        OUTPUT: the optimal value taken by the objective function

        TESTS:

        The SDP from the header of this module::

            sage: p = SemidefiniteProgram(maximization=False)
            sage: x = p.new_variable()
            sage: p.set_objective(x[0] - x[1])
            sage: a1 = matrix([[1, 2.], [2., 3.]])
            sage: a2 = matrix([[3, 4.], [4., 2.]])
            sage: a3 = matrix([[5, 6.], [6., 7.]])
            sage: b1 = matrix([[1, 1.], [1., 1.]])
            sage: b2 = matrix([[2, 2.], [2., 1.]])
            sage: b3 = matrix([[3, 3.], [3., 3.]])
            sage: p.add_constraint(a1*x[0] + a2*x[1] <= a3)
            sage: p.add_constraint(b1*x[0] + b2*x[1] <= b3)

            sage: # needs cvxopt
            sage: N(p.solve(), 4)
            -11.
            sage: x = p.get_values(x)
            sage: N(x[0],4)
            -8.0
            sage: N(x[1],4)
            3.0
        """
        self._backend.solve()
        return self._backend.get_objective_value()

    cpdef dual_variable(self, int i, sparse=False):
        """
        The `i`-th dual variable.

        Available after ``self.solve()`` is called, otherwise the result is
        undefined.

        INPUT:

        - ``index`` -- integer; the constraint's id

        OUTPUT: the matrix of the `i`-th dual variable

        EXAMPLES:

        Dual objective value is the same as the primal one::

            sage: p = SemidefiniteProgram(maximization=False)
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

            sage: # needs cvxopt
            sage: p.solve()                                                              # tol 1e-08
            -3.0
            sage: x = p.get_values(x).values()
            sage: -(a3*p.dual_variable(0)).trace() - (b3*p.dual_variable(1)).trace()     # tol 1e-07
            -3.0

        Dual variable is orthogonal to the slack ::

            sage: # needs cvxopt
            sage: g = sum((p.slack(j)*p.dual_variable(j)).trace() for j in range(2)); g  # tol 1.2e-08
            0.0

        TESTS::

            sage: p.dual_variable(7)                                                    # needs cvxopt
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self._backend.dual_variable(i, sparse=sparse)

    cpdef slack(self, int i, sparse=False):
        """
        Slack of the `i`-th constraint.

        Available after ``self.solve()`` is called, otherwise the result is
        undefined.

        INPUT:

        - ``index`` -- integer; the constraint's id

        OUTPUT: the matrix of the slack of the `i`-th constraint

        EXAMPLES::

            sage: p = SemidefiniteProgram(maximization = False)
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

            sage: # needs cvxopt
            sage: p.solve()                         # tol 1e-08
            -3.0
            sage: B1 = p.slack(1); B1               # tol 1e-08
            [0.0 0.0]
            [0.0 0.0]
            sage: B1.is_positive_definite()
            True
            sage: x = sorted(p.get_values(x).values())
            sage: x[0]*b1 + x[1]*b2 - b3 + B1       # tol 1e-09
            [0.0 0.0]
            [0.0 0.0]

        TESTS::

            sage: p.slack(7)                                                            # needs cvxopt
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self._backend.slack(i, sparse=sparse)

    def solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter.

        The solver parameters are by essence solver-specific, which
        means their meaning heavily depends on the solver used.

        (If you do not know which solver you are using, then you are
        using CVXOPT).

        INPUT:

        - ``name`` -- string; the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value

        EXAMPLES::

            sage: # needs cvxopt
            sage: p.<x> = SemidefiniteProgram(solver='cvxopt',
            ....:                             maximization=False)
            sage: p.solver_parameter("show_progress", True)
            sage: p.solver_parameter("show_progress")
            True
            sage: p.set_objective(x[0] - x[1])
            sage: a1 = matrix([[1, 2.], [2., 3.]])
            sage: a2 = matrix([[3, 4.], [4., 2.]])
            sage: a3 = matrix([[5, 6.], [6., 7.]])
            sage: b1 = matrix([[1, 1.], [1., 1.]])
            sage: b2 = matrix([[2, 2.], [2., 1.]])
            sage: b3 = matrix([[3, 3.], [3., 3.]])
            sage: p.add_constraint(a1*x[0] + a2*x[1] <= a3)
            sage: p.add_constraint(b1*x[0] + b2*x[1] <= b3)
            sage: N(p.solve(), 4)
                 pcost       dcost       gap    pres   dres   k/t
             0:  1...
            ...
            Optimal solution found.
            -11.
        """
        if value is None:
            return self._backend.solver_parameter(name)
        else:
            self._backend.solver_parameter(name, value)

    cpdef sum(self, L):
        r"""
        Efficiently compute the sum of a sequence of
        :class:`~sage.numerical.linear_functions.LinearFunction` elements.

        INPUT:

        - ``L`` -- list of
          :class:`~sage.numerical.linear_functions.LinearFunction` instances

        .. NOTE::

            The use of the regular ``sum`` function is not recommended
            as it is much less efficient than this one.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable()

        The following command::

            sage: s = p.sum(v[i] for i in range(90))

        is much more efficient than::

            sage: s = sum(v[i] for i in range(90))
        """
        d = {}
        for v in L:
            for id, coeff in v.iteritems():
                d[id] = coeff + d.get(id, 0)
        return self.linear_functions_parent()(d)

    def get_backend(self):
        r"""
        Return the backend instance used.

        This might be useful when access to additional functions provided by
        the backend is needed.

        EXAMPLES:

        This example prints a matrix coefficient::

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
        return self._backend


class SDPSolverException(RuntimeError):
    r"""
    Exception raised when the solver fails.

    ``SDPSolverException`` is the exception raised when the solver fails.

    EXAMPLES::

        sage: from sage.numerical.sdp import SDPSolverException
        sage: SDPSolverException("Error")
        SDPSolverException('Error'...)

    TESTS:

    No solution::

        sage: # needs cvxopt
        sage: p = SemidefiniteProgram(solver='cvxopt')
        sage: x = p.new_variable()
        sage: p.set_objective(x[0])
        sage: a = matrix([[1,2],[2,4]])
        sage: b = matrix([[1,9],[9,4]])
        sage: p.add_constraint(a*x[0] == b)
        sage: p.solve()
        Traceback (most recent call last):
        ...
        SDPSolverException: ...

    The value of the exception::

        sage: from sage.numerical.sdp import SDPSolverException
        sage: e = SDPSolverException("Error")
        sage: print(e)
        Error
    """
    pass


cdef class SDPVariable(Element):
    r"""
    ``SDPVariable`` is a variable used by the class
    ``SemidefiniteProgram``.

    .. warning::

        You should not instantiate this class directly. Instead, use
        :meth:`SemidefiniteProgram.new_variable`.
    """

    def __init__(self, parent, sdp, name):
        r"""
        Constructor for ``SDPVariable``.

        INPUT:

        - ``parent`` -- :class:`SDPVariableParent`; the parent of the
          SDP variable

        - ``sdp`` -- :class:`SemidefiniteProgram`; the
          underlying linear program

        - ``name`` -- a name for the ``SDPVariable``

        - ``lower_bound``, ``upper_bound`` -- lower bound and upper
          bound on the variable. Set to ``None`` to indicate that the
          variable is unbounded.

        For more informations, see the method
        ``SemidefiniteProgram.new_variable``.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: p.new_variable()
            SDPVariable
        """
        super().__init__(parent)
        self._dict = {}
        self._p = sdp
        self._name = name

    def __getitem__(self, i):
        r"""
        Return the symbolic variable corresponding to the key.

        Returns the element asked, otherwise creates it.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[0] + v[1])
            sage: v[0]
            x_0
        """
        cdef int j
        if i in self._dict:
            return self._dict[i]
        zero = self._p._backend.zero()
        name = self._name + "[" + str(i) + "]" if self._name else None
        j = self._p._backend.add_variable( obj=zero, name=name)
        v = self._p.linear_function({j : 1})
        self._p._variables[v] = j
        self._dict[i] = v
        return v

    def _repr_(self):
        r"""
        Return a representation of ``self``.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable()
            sage: v
            SDPVariable
        """
        return "SDPVariable"

    def keys(self):
        r"""
        Return the keys already defined in the dictionary.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[0] + v[1])
            sage: sorted(v.keys())
            [0, 1]
        """
        return self._dict.keys()

    def items(self):
        r"""
        Return the pairs (keys,value) contained in the dictionary.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[0] + v[1])
            sage: sorted(v.items())
            [(0, x_0), (1, x_1)]
        """
        return self._dict.items()

    def values(self):
        r"""
        Return the symbolic variables associated to the current dictionary.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[0] + v[1])
            sage: sorted(v.values(), key=str)
            [x_0, x_1]
        """
        return self._dict.values()

    cdef _matrix_rmul_impl(self, m):
        """
        Implement the action of a matrix multiplying from the right.
        """
        result = dict()
        for i, row in enumerate(m.rows()):
            x = self[i]
            x_index, = x.dict().keys()
            result[x_index] = row
        from sage.modules.free_module import FreeModule
        V = FreeModule(self._p.base_ring(), m.ncols())
        T = self._p.linear_functions_parent().tensor(V)
        return T(result)

    cdef _matrix_lmul_impl(self, m):
        """
        Implement the action of a matrix multiplying from the left.
        """
        result = dict()
        for i, col in enumerate(m.columns()):
            x = self[i]
            x_index, = x.dict().keys()
            result[x_index] = col
        from sage.modules.free_module import FreeModule
        V = FreeModule(self._p.base_ring(), m.nrows())
        T = self._p.linear_functions_parent().tensor(V)
        return T(result)

    cpdef _acted_upon_(self, mat, bint self_on_left):
        """
        Act with matrices on SDPVariables.

        EXAMPLES::

            sage: p = SemidefiniteProgram()
            sage: v = p.new_variable()
            sage: m = matrix([[1,2], [3,4]])
            sage: v * m
            (1.0, 2.0)*x_0 + (3.0, 4.0)*x_1
            sage: m * v
            (1.0, 3.0)*x_0 + (2.0, 4.0)*x_1
        """
        from sage.structure.element import Matrix
        if isinstance(mat, Matrix):
            return self._matrix_rmul_impl(mat) if self_on_left else self._matrix_lmul_impl(mat)


cdef class SDPVariableParent(Parent):
    """
    Parent for :class:`SDPVariable`.

    .. warning::

        This class is for internal use. You should not instantiate it
        yourself. Use :meth:`SemidefiniteProgram.new_variable`
        to generate sdp variables.
    """

    Element = SDPVariable

    def _repr_(self):
        r"""
        Return representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: sdp.<v> = SemidefiniteProgram()
            sage: v.parent()
            Parent of SDPVariables
        """
        return 'Parent of SDPVariables'

    def _an_element_(self):
        """
        Construct a SDP variable.

        OUTPUT:

        This is required for the coercion framework. We raise a
        :exc:`TypeError` to abort search for any coercion to another
        parent for binary operations. The only interesting operations
        involving :class:`SDPVariable` elements are actions by
        matrices.

        EXAMPLES::

            sage: sdp.<x> = SemidefiniteProgram()
            sage: parent = x.parent()
            sage: parent.an_element()    # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: disallow coercion
        """
        raise TypeError('disallow coercion')

    def _element_constructor_(self, sdp, name=""):
        """
        The Element constructor.

        INPUT/OUTPUT: see :meth:`SDPVariable.__init__`

        EXAMPLES::

            sage: sdp = SemidefiniteProgram()
            sage: sdp.new_variable()    # indirect doctest
            SDPVariable
        """
        return self.element_class(self, sdp, name)


sdp_variable_parent = SDPVariableParent()
