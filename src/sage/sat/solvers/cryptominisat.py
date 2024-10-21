# sage_setup: distribution = sagemath-combinat
r"""
CryptoMiniSat Solver

This solver relies on Python bindings provided by upstream cryptominisat.

AUTHORS:

- Thierry Monteil (2017): complete rewrite, using upstream Python bindings,
  works with cryptominisat 5.
- Martin Albrecht (2012): first version, as a cython interface, works with
  cryptominisat 2.
"""

# ****************************************************************************
#       Copyright (C) 2017 Thierry Monteil <sage!lma.metelu.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from .satsolver import SatSolver

from sage.misc.lazy_import import lazy_import
from sage.features.sat import Pycryptosat

lazy_import('pycryptosat', ['Solver'], feature=Pycryptosat())


class CryptoMiniSat(SatSolver):
    r"""
    CryptoMiniSat Solver.

    INPUT:

    - ``verbosity`` -- integer between 0 and 15 (default: 0)

    - ``confl_limit`` -- integer (default: ``None``); abort after this many
      conflicts. If set to ``None``, never aborts.

    - ``threads`` -- integer (default: ``None``); the number of thread to
      use. If set to ``None``, the number of threads used corresponds to the
      number of cpus.

    EXAMPLES::

        sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat
        sage: solver = CryptoMiniSat()                                  # optional - pycryptosat
    """
    def __init__(self, verbosity=0, confl_limit=None, threads=None):
        r"""
        Construct a new CryptoMiniSat instance.

        See the documentation class for the description of inputs.

        EXAMPLES::

            sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat
            sage: solver = CryptoMiniSat(threads=1)                     # optional - pycryptosat
        """
        if threads is None:
            from sage.parallel.ncpus import ncpus
            threads = ncpus()
        if confl_limit is None:
            from sys import maxsize
            confl_limit = maxsize
        self._solver = Solver(verbose=int(verbosity), confl_limit=int(confl_limit), threads=int(threads))
        self._nvars = 0
        self._clauses = []

    def var(self, decision=None):
        r"""
        Return a *new* variable.

        INPUT:

        - ``decision`` -- accepted for compatibility with other solvers; ignored

        EXAMPLES::

            sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat
            sage: solver = CryptoMiniSat()                                  # optional - pycryptosat
            sage: solver.var()                                              # optional - pycryptosat
            1

            sage: solver.add_clause((-1,2,-4))                              # optional - pycryptosat
            sage: solver.var()                                              # optional - pycryptosat
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

            sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat
            sage: solver = CryptoMiniSat()                                  # optional - pycryptosat
            sage: solver.nvars()                                            # optional - pycryptosat
            0

        If a variable with intermediate index is not used, it is still
        considered as a variable::

            sage: solver.add_clause((1,-2,4))                               # optional - pycryptosat
            sage: solver.nvars()                                            # optional - pycryptosat
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

            sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat
            sage: solver = CryptoMiniSat()                                  # optional - pycryptosat
            sage: solver.add_clause((1, -2 , 3))                            # optional - pycryptosat
        """
        if 0 in lits:
            raise ValueError("0 should not appear in the clause: {}".format(lits))
        # cryptominisat does not handle Sage integers
        lits = tuple(int(i) for i in lits)
        self._nvars = max(self._nvars, max(abs(i) for i in lits))
        self._solver.add_clause(lits)
        self._clauses.append((lits, False, None))

    def add_xor_clause(self, lits, rhs=True):
        r"""
        Add a new XOR clause to set of clauses.

        INPUT:

        - ``lits`` -- tuple of positive integers

        - ``rhs`` -- boolean (default: ``True``); whether this XOR clause should
          be evaluated to ``True`` or ``False``

        EXAMPLES::

            sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat
            sage: solver = CryptoMiniSat()                                  # optional - pycryptosat
            sage: solver.add_xor_clause((1, 2 , 3), False)                  # optional - pycryptosat
        """
        if 0 in lits:
            raise ValueError("0 should not appear in the clause: {}".format(lits))
        # cryptominisat does not handle Sage integers
        lits = tuple(int(i) for i in lits)
        self._nvars = max(self._nvars, max(abs(i) for i in lits))
        self._solver.add_xor_clause(lits, rhs)
        self._clauses.append((lits, True, rhs))

    def __call__(self, assumptions=None):
        r"""
        Solve this instance.

        OUTPUT:

        - If this instance is SAT: A tuple of length ``nvars()+1``
          where the ``i``-th entry holds an assignment for the
          ``i``-th variables (the ``0``-th entry is always ``None``).

        - If this instance is UNSAT: ``False``.

        EXAMPLES::

            sage: # optional - pycryptosat
            sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat
            sage: solver = CryptoMiniSat()
            sage: solver.add_clause((1,2))
            sage: solver.add_clause((-1,2))
            sage: solver.add_clause((-1,-2))
            sage: solver()
            (None, False, True)

            sage: solver.add_clause((1,-2))                                 # optional - pycryptosat
            sage: solver()                                                  # optional - pycryptosat
            False
        """
        satisfiable, assignments = self._solver.solve()
        if satisfiable:
            return assignments
        else:
            return False

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat
            sage: solver = CryptoMiniSat()                                  # optional - pycryptosat
            sage: solver                                                    # optional - pycryptosat
            CryptoMiniSat solver: 0 variables, 0 clauses.
        """
        return "CryptoMiniSat solver: {} variables, {} clauses.".format(self.nvars(), len(self.clauses()))

    def clauses(self, filename=None):
        r"""
        Return original clauses.

        INPUT:

        - ``filename`` -- if not ``None`` clauses are written to ``filename`` in
          DIMACS format (default: ``None``)

        OUTPUT:

            If ``filename`` is ``None`` then a list of ``lits, is_xor, rhs``
            tuples is returned, where ``lits`` is a tuple of literals,
            ``is_xor`` is always ``False`` and ``rhs`` is always ``None``.

            If ``filename`` points to a writable file, then the list of original
            clauses is written to that file in DIMACS format.

        EXAMPLES::

            sage: # optional - pycryptosat
            sage: from sage.sat.solvers import CryptoMiniSat
            sage: solver = CryptoMiniSat()
            sage: solver.add_clause((1,2,3,4,5,6,7,8,-9))
            sage: solver.add_xor_clause((1,2,3,4,5,6,7,8,9), rhs=True)
            sage: solver.clauses()
            [((1, 2, 3, 4, 5, 6, 7, 8, -9), False, None),
            ((1, 2, 3, 4, 5, 6, 7, 8, 9), True, True)]

        DIMACS format output::

            sage: # optional - pycryptosat
            sage: from sage.sat.solvers import CryptoMiniSat
            sage: solver = CryptoMiniSat()
            sage: solver.add_clause((1, 2, 4))
            sage: solver.add_clause((1, 2, -4))
            sage: fn = tmp_filename()
            sage: solver.clauses(fn)
            sage: print(open(fn).read())
            p cnf 4 2
            1 2 4 0
            1 2 -4 0
            <BLANKLINE>

        Note that in cryptominisat, the DIMACS standard format is augmented with
        the following extension: having an ``x`` in front of a line makes that
        line an XOR clause::

            sage: solver.add_xor_clause((1,2,3), rhs=True)      # optional - pycryptosat
            sage: solver.clauses(fn)                            # optional - pycryptosat
            sage: print(open(fn).read())                        # optional - pycryptosat
            p cnf 4 3
            1 2 4 0
            1 2 -4 0
            x1 2 3 0
            <BLANKLINE>

        Note that inverting an xor-clause is equivalent to inverting one of the
        variables::

            sage: solver.add_xor_clause((1,2,5),rhs=False)      # optional - pycryptosat
            sage: solver.clauses(fn)                            # optional - pycryptosat
            sage: print(open(fn).read())                        # optional - pycryptosat
            p cnf 5 4
            1 2 4 0
            1 2 -4 0
            x1 2 3 0
            x1 2 -5 0
            <BLANKLINE>
        """
        if filename is None:
            return self._clauses
        from sage.sat.solvers.dimacs import DIMACS
        DIMACS.render_dimacs(self._clauses, filename, self.nvars())
