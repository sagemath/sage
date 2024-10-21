# sage_setup: distribution = sagemath-combinat
r"""
PicoSAT Solver

This solver relies on the ``pycosat`` Python bindings to ``PicoSAT``.

The ``pycosat`` package should be installed on your Sage installation.

AUTHORS:

- Thierry Monteil (2018): initial version.
"""

# ****************************************************************************
#       Copyright (C) 2018 Thierry Monteil <sage!lma.metelu.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from .satsolver import SatSolver

from sage.misc.lazy_import import lazy_import
from sage.features.sat import Pycosat

lazy_import('pycosat', ['solve'], feature=Pycosat())


class PicoSAT(SatSolver):
    r"""
    PicoSAT Solver.

    INPUT:

    - ``verbosity`` -- integer between 0 and 2 (default: 0)

    - ``prop_limit`` -- integer (default: 0); the propagation limit

    EXAMPLES::

        sage: from sage.sat.solvers.picosat import PicoSAT
        sage: solver = PicoSAT()                           # optional - pycosat
    """
    def __init__(self, verbosity=0, prop_limit=0):
        r"""
        Construct a new PicoSAT instance.

        See the documentation class for the description of inputs.

        EXAMPLES::

            sage: from sage.sat.solvers.picosat import PicoSAT
            sage: solver = PicoSAT()                       # optional - pycosat
        """
        self._verbosity = int(verbosity)
        if prop_limit is None:
            self._prop_limit = 0
        else:
            self._prop_limit = int(prop_limit)
        self._solve = solve
        self._nvars = 0
        self._clauses = []

    def var(self, decision=None):
        r"""
        Return a *new* variable.

        INPUT:

        - ``decision`` -- ignored; accepted for compatibility with other solvers

        EXAMPLES::

            sage: from sage.sat.solvers.picosat import PicoSAT
            sage: solver = PicoSAT()                       # optional - pycosat
            sage: solver.var()                             # optional - pycosat
            1

            sage: solver.add_clause((-1,2,-4))             # optional - pycosat
            sage: solver.var()                             # optional - pycosat
            5
        """
        self._nvars += 1
        return self._nvars

    def nvars(self):
        r"""
        Return the number of variables.

        Note that for compatibility with DIMACS convention, the number
        of variables corresponds to the maximal index of the variables used.

        EXAMPLES::

            sage: from sage.sat.solvers.picosat import PicoSAT
            sage: solver = PicoSAT()                       # optional - pycosat
            sage: solver.nvars()                           # optional - pycosat
            0

        If a variable with intermediate index is not used, it is still
        considered as a variable::

            sage: solver.add_clause((1,-2,4))              # optional - pycosat
            sage: solver.nvars()                           # optional - pycosat
            4
        """
        return self._nvars

    def add_clause(self, lits):
        r"""
        Add a new clause to set of clauses.

        INPUT:

        - ``lits`` -- tuple of nonzero integers

        .. NOTE::

            If any element ``e`` in ``lits`` has ``abs(e)`` greater
            than the number of variables generated so far, then new
            variables are created automatically.

        EXAMPLES::

            sage: from sage.sat.solvers.picosat import PicoSAT
            sage: solver = PicoSAT()                       # optional - pycosat
            sage: solver.add_clause((1, -2 , 3))           # optional - pycosat
        """
        if 0 in lits:
            raise ValueError("0 should not appear in the clause: {}".format(lits))
        # pycosat does not handle Sage integers
        lits = [int(i) for i in lits]
        self._nvars = max(self._nvars, max(abs(i) for i in lits))
        self._clauses.append(lits)

    def __call__(self, assumptions=None):
        r"""
        Solve this instance.

        OUTPUT:

        - If this instance is SAT: A tuple of length ``nvars() + 1``,
          where the ``i``-th entry holds an assignment for the
          ``i``-th variables (the ``0``-th entry is always ``None``).

        - If this instance is UNSAT: ``False``.

        EXAMPLES::

            sage: # optional - pycosat
            sage: from sage.sat.solvers.picosat import PicoSAT
            sage: solver = PicoSAT()
            sage: solver.add_clause((1,2))
            sage: solver.add_clause((-1,2))
            sage: solver.add_clause((-1,-2))
            sage: solver()
            (None, False, True)

            sage: solver.add_clause((1,-2))                # optional - pycosat
            sage: solver()                                 # optional - pycosat
            False
        """
        #import pycosat
        #self._solve = pycosat.solve
        sol = self._solve(self._clauses, verbose=self._verbosity,
                          prop_limit=self._prop_limit, vars=self._nvars)
        # sol = pycosat.solve(self._clauses)
        if sol == 'UNSAT':
            return False
        else:
            return (None,) + tuple([s > 0 for s in sol])

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.sat.solvers.picosat import PicoSAT
            sage: solver = PicoSAT()                       # optional - pycosat
            sage: solver                                   # optional - pycosat
            PicoSAT solver: 0 variables, 0 clauses.
        """
        return "PicoSAT solver: {} variables, {} clauses.".format(self.nvars(), len(self.clauses()))

    def clauses(self, filename=None):
        r"""
        Return original clauses.

        INPUT:

        - ``filename`` -- (optional) if given, clauses are written to
          ``filename`` in DIMACS format

        OUTPUT:

        If ``filename`` is ``None`` then a list of ``lits`` is returned,
        where ``lits`` is a list of literals.

        If ``filename`` points to a writable file, then the list of original
        clauses is written to that file in DIMACS format.

        EXAMPLES::

            sage: from sage.sat.solvers.picosat import PicoSAT
            sage: solver = PicoSAT()                       # optional - pycosat
            sage: solver.add_clause((1,2,3,4,5,6,7,8,-9))  # optional - pycosat
            sage: solver.clauses()                         # optional - pycosat
            [[1, 2, 3, 4, 5, 6, 7, 8, -9]]

        DIMACS format output::

            sage: # optional - pycosat
            sage: from sage.sat.solvers.picosat import PicoSAT
            sage: solver = PicoSAT()
            sage: solver.add_clause((1, 2, 4))
            sage: solver.add_clause((1, 2, -4))
            sage: fn = tmp_filename()
            sage: solver.clauses(fn)
            sage: print(open(fn).read())
            p cnf 4 2
            1 2 4 0
            1 2 -4 0
            <BLANKLINE>
        """
        if filename is None:
            return self._clauses
        from sage.sat.solvers.dimacs import DIMACS
        DIMACS.render_dimacs(self._clauses, filename, self.nvars())
