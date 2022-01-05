r"""
Solve SAT problems Integer Linear Programming

The class defined here is a :class:`~sage.sat.solvers.satsolver.SatSolver` that
solves its instance using :class:`MixedIntegerLinearProgram`. Its performance
can be expected to be slower than when using
:class:`~sage.sat.solvers.cryptominisat.cryptominisat.CryptoMiniSat`.
"""
from .satsolver import SatSolver
from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

class SatLP(SatSolver):
    def __init__(self, solver=None, verbose=0, *, integrality_tolerance=1e-3):
        r"""
        Initializes the instance

        INPUT:

        - ``solver`` -- (default: ``None``) Specify a Mixed Integer Linear Programming
          (MILP) solver to be used. If set to ``None``, the default one is used. For
          more information on MILP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of verbosity
          of the LP solver. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- parameter for use with MILP solvers over an
          inexact base ring; see :meth:`MixedIntegerLinearProgram.get_values`.

        EXAMPLES::

            sage: S=SAT(solver="LP"); S
            an ILP-based SAT Solver
        """
        SatSolver.__init__(self)
        self._LP = MixedIntegerLinearProgram(solver=solver)
        self._LP_verbose = verbose
        self._vars = self._LP.new_variable(binary=True)
        self._integrality_tolerance = integrality_tolerance

    def var(self):
        """
        Return a *new* variable.

        EXAMPLES::

            sage: S=SAT(solver="LP"); S
            an ILP-based SAT Solver
            sage: S.var()
            1
        """
        nvars = n = self._LP.number_of_variables()
        while nvars==self._LP.number_of_variables():
            n += 1
            self._vars[n] # creates the variable if needed
        return n

    def nvars(self):
        """
        Return the number of variables.

        EXAMPLES::

            sage: S=SAT(solver="LP"); S
            an ILP-based SAT Solver
            sage: S.var()
            1
            sage: S.var()
            2
            sage: S.nvars()
            2
        """
        return self._LP.number_of_variables()

    def add_clause(self, lits):
        """
        Add a new clause to set of clauses.

        INPUT:

        - ``lits`` - a tuple of integers != 0

        .. note::

            If any element ``e`` in ``lits`` has ``abs(e)`` greater
            than the number of variables generated so far, then new
            variables are created automatically.

        EXAMPLES::

            sage: S=SAT(solver="LP"); S
            an ILP-based SAT Solver
            sage: for u,v in graphs.CycleGraph(6).edges(labels=False):
            ....:     u,v = u+1,v+1
            ....:     S.add_clause((u,v))
            ....:     S.add_clause((-u,-v))
        """
        if 0 in lits:
            raise ValueError("0 should not appear in the clause: {}".format(lits))
        p = self._LP
        p.add_constraint(p.sum(self._vars[x] if x>0 else 1-self._vars[-x] for x in lits)
                         >=1)

    def __call__(self):
        """
        Solve this instance.

        OUTPUT:

        - If this instance is SAT: A tuple of length ``nvars()+1``
          where the ``i``-th entry holds an assignment for the
          ``i``-th variables (the ``0``-th entry is always ``None``).

        - If this instance is UNSAT: ``False``

        EXAMPLES::

            sage: def is_bipartite_SAT(G):
            ....:     S=SAT(solver="LP"); S
            ....:     for u,v in G.edges(labels=False):
            ....:         u,v = u+1,v+1
            ....:         S.add_clause((u,v))
            ....:         S.add_clause((-u,-v))
            ....:     return S
            sage: S = is_bipartite_SAT(graphs.CycleGraph(6))
            sage: S() # random
            [None, True, False, True, False, True, False]
            sage: True in S()
            True
            sage: S = is_bipartite_SAT(graphs.CycleGraph(7))
            sage: S()
            False
        """
        try:
            self._LP.solve(log=self._LP_verbose)
        except MIPSolverException:
            return False

        b = self._LP.get_values(self._vars, convert=bool, tolerance=self._integrality_tolerance)
        n = max(b)
        return [None]+[b.get(i, False) for i in range(1,n+1)]

    def __repr__(self):
        """
        TESTS::

            sage: S=SAT(solver="LP"); S
            an ILP-based SAT Solver
        """
        return "an ILP-based SAT Solver"
