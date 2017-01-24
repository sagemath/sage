"""
Abstract SAT Solver

All SAT solvers must inherit from this class.

.. NOTE::

    Our SAT solver interfaces are 1-based, i.e., literals start at
    1. This is consistent with the popular DIMACS format for SAT
    solving but not with Pythion's 0-based convention. However, this
    also allows to construct clauses using simple integers.

AUTHORS:

- Martin Albrecht (2012): first version
"""
from sage.misc.package import PackageNotFoundError

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

        - ``decision`` - is this variable a decision variable?

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

        .. NOTE::

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

    def read(self, filename):
        r"""
        Reads DIMAC files.

        Reads in DIMAC formatted lines (lazily) from a
        file or file object and adds the corresponding
        clauses into this solver instance. Note that the
        DIMACS format is not well specified, see
        http://people.sc.fsu.edu/~jburkardt/data/cnf/cnf.html,
        http://www.satcompetition.org/2009/format-benchmarks2009.html,
        and http://elis.dvo.ru/~lab_11/glpk-doc/cnfsat.pdf.
        The differences were summarized in the discussion on
        the ticket :trac:`16924`. This method assumes the following
        DIMACS format

        - Any line starting with "c" is a comment
        - Any line starting with "p" is a header
        - Any variable 1-n can be used
        - Every line containing a clause must end with a "0"

        INPUT:

        - ``filename`` - The name of a file as a string or a file object

        EXAMPLE::

            sage: from six import StringIO # for python 2/3 support
            sage: file_object = StringIO("c A sample .cnf file.\np cnf 3 2\n1 -3 0\n2 3 -1 0 ")
            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: solver = DIMACS()
            sage: solver.read(file_object)
            sage: solver.clauses()
            [((1, -3), False, None), ((2, 3, -1), False, None)]
        """
        if isinstance(filename,str):
            file_object = open(filename,"r")
        else:
            file_object = filename
        for line in file_object:
            if line.startswith("c"):
                continue # comment
            if line.startswith("p"):
                continue # header
            line = line.split(" ")
            clause = map(int,[e for e in line if e])
            clause = clause[:-1] # strip trailing zero
            self.add_clause(clause)


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

def SAT(solver=None):
    r"""
    Return a :class:`SatSolver` instance.

    Through this class, one can define and solve `SAT
    <https://en.wikipedia.org/wiki/Boolean_satisfiability_problem>`__ problems.

    INPUT:

    - ``solver`` (string) -- select a solver. Admissible values are:

        - ``"cryptominisat"`` -- note that the cryptominisat package must be
          installed.

        - ``"LP"`` -- use :class:`~sage.sat.solvers.sat_lp.SatLP` to solve the
          SAT instance.

        - ``None`` (default) -- use CryptoMiniSat if available, and a LP solver
          otherwise.

    EXAMPLE::

        sage: SAT(solver="LP")
        an ILP-based SAT Solver

    TESTS::

        sage: SAT(solver="Wouhouuuuuu")
        Traceback (most recent call last):
        ...
        ValueError: Solver 'Wouhouuuuuu' is not available

    Forcing CryptoMiniSat::

        sage: SAT(solver="cryptominisat") # optional - cryptominisat
        CryptoMiniSat
        #vars:       0, #lits:       0, #clauses:       0, #learnt:       0, #assigns:       0

    """
    if solver is None:
        try:
            from sage.sat.solvers.cryptominisat.cryptominisat import CryptoMiniSat
            solver = "cryptominisat"
        except ImportError:
            solver = "LP"

    if solver == 'cryptominisat':
        try:
            from sage.sat.solvers.cryptominisat.cryptominisat import CryptoMiniSat
        except ImportError:
            raise PackageNotFoundError("cryptominisat")
        return CryptoMiniSat()
    elif solver == "LP":
        from sat_lp import SatLP
        return SatLP()
    else:
        raise ValueError("Solver '{}' is not available".format(solver))
