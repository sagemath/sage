"""
Abstract SAT Solver

All SAT solvers must inherit from this class.

.. note::

    Our SAT solver interfaces are 1-based, i.e., literals start at
    1. This is consistent with the popular DIMACS format for SAT
    solving but not with Pythion's 0-based convention. However, this
    also allows to construct clauses using simple integers.

AUTHORS:

- Martin Albrecht (2012): first version
"""

cdef class SatSolver:
    def __cinit__(self, *args, **kwds):
        """
        Constuct a new SATSolver.

        EXAMPLE::

            sage: from sage.sat.solvers.satsolver import SatSolver
            sage: solver = SatSolver()
        """
        pass

    def var(self, decision=None):
        """
        Return a *new* variable.

        INPUT:

        - ``decision`` - is this variable a deicison variable?

        EXAMPLE::

            sage: from sage.sat.solvers.satsolver import SatSolver
            sage: solver = SatSolver()
            sage: solver.var()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def nvars(self):
        """
        Return the number of variables.

        EXAMPLE::

            sage: from sage.sat.solvers.satsolver import SatSolver
            sage: solver = SatSolver()
            sage: solver.nvars()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def add_clause(self, lits):
        """
        Add a new clause to set of clauses.

        INPUT:

        - ``lits`` - a tuple of integers != 0

        .. note::

            If any element ``e`` in ``lits`` has ``abs(e)`` greater
            than the number of variables generated so far, then new
            variables are created automatically.

        EXAMPLE::

            sage: from sage.sat.solvers.satsolver import SatSolver
            sage: solver = SatSolver()
            sage: solver.add_clause( (1, -2 , 3) )
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __call__(self, assumptions=None):
        """
        Solve this instance.

        INPUT:

        - ``assumptions`` - assumed variable assignments (default: ``None``)

        OUTPUT:

        - If this instance is SAT: A tuple of length ``nvars()+1``
          where the ``i``-th entry holds an assignment for the
          ``i``-th variables (the ``0``-th entry is always ``None``).

        - If this instance is UNSAT: ``False``

        - If the solver was interrupted before deciding satisfiability
          ``None``.

        EXAMPLE::

            sage: from sage.sat.solvers.satsolver import SatSolver
            sage: solver = SatSolver()
            sage: solver()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def conflict_clause(self):
        """
        Return conflict clause if this instance is UNSAT and the last
        call used assumptions.

        EXAMPLE::

            sage: from sage.sat.solvers.satsolver import SatSolver
            sage: solver = SatSolver()
            sage: solver.conflict_clause()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def learnt_clauses(self, unitary_only=False):
        """
        Return learnt clauses.

        INPUT:

        - ``unitary_only`` - return only unitary learnt clauses (default: ``False``)

        EXAMPLE::

            sage: from sage.sat.solvers.satsolver import SatSolver
            sage: solver = SatSolver()
            sage: solver.learnt_clauses()
            Traceback (most recent call last):
            ...
            NotImplementedError

            sage: solver.learnt_clauses(unitary_only=True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __repr__(self):
        """
        TESTS::

        sage: from sage.sat.solvers.satsolver import SatSolver
        sage: solver = SatSolver()
        sage: solver
        a generic SAT solver (don't use me, inherit from me)
        """
        return "a generic SAT solver (don't use me, inherit from me)"

    def clauses(self, filename=None):
        """
        Return original clauses.

        INPUT:

        - ``filename'' - if not ``None`` clauses are written to ``filename`` in
          DIMACS format (default: ``None``)

        OUTPUT:

            If ``filename`` is ``None`` then a list of ``lits, is_xor, rhs``
            tuples is returned, where ``lits`` is a tuple of literals,
            ``is_xor`` is always ``False`` and ``rhs`` is always ``None``.

            If ``filename`` points to a writable file, then the list of original
            clauses is written to that file in DIMACS format.


        EXAMPLE::

            sage: from sage.sat.solvers.satsolver import SatSolver
            sage: solver = SatSolver()
            sage: solver.clauses()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __getattr__(self, name):
        """
        EXAMPLE::

            sage: from sage.sat.solvers.satsolver import SatSolver
            sage: solver = SatSolver()
            sage: solver.gens() # __getattr__ points this to clauses
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if name == "gens":
            return self.clauses
        else:
            raise AttributeError("'%s' object has no attribute '%s'"%(self.__class__.__name__,name))

    def trait_names(self):
        """
        Allow alias to appear in tab completion.

        EXAMPLE::

            sage: from sage.sat.solvers.satsolver import SatSolver
            sage: solver = SatSolver()
            sage: solver.trait_names()
            ['gens']
        """
        return ["gens"]

