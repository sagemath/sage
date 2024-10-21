# sage_setup: distribution = sagemath-combinat
"""
SAT-Solvers via DIMACS Files

Sage supports calling SAT solvers using the popular DIMACS format. This module implements
infrastructure to make it easy to add new such interfaces and some example interfaces.

Currently, interfaces to **RSat** and **Glucose** are included by default.

.. NOTE::

    Our SAT solver interfaces are 1-based, i.e., literals start at 1. This is consistent with the
    popular DIMACS format for SAT solving but not with Pythion's 0-based convention. However, this
    also allows to construct clauses using simple integers.

AUTHORS:

- Martin Albrecht (2012): first version
- Sébastien Labbé (2018): adding Glucose SAT solver
- Sébastien Labbé (2023): adding Kissat SAT solver

Classes and Methods
-------------------
"""
##############################################################################
#  Copyright (C) 2012 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
##############################################################################

import os
import sys
import subprocess
import shlex

from sage.sat.solvers.satsolver import SatSolver
from sage.misc.temporary_file import tmp_filename
from time import sleep


class DIMACS(SatSolver):
    """
    Generic DIMACS Solver.

    .. NOTE::

        Usually, users won't have to use this class directly but some
        class which inherits from this class.

    .. automethod:: __init__
    .. automethod:: __call__
    """

    command = ""

    def __init__(self, command=None, filename=None, verbosity=0, **kwds):
        """
        Construct a new generic DIMACS solver.

        INPUT:

        - ``command`` -- a named format string with the command to
          run. The string must contain {input} and may contain
          {output} if the solvers writes the solution to an output
          file. For example "sat-solver {input}" is a valid
          command. If ``None`` then the class variable ``command`` is
          used. (default: ``None``)

        - ``filename`` -- a filename to write clauses to in DIMACS
          format, must be writable. If ``None`` a temporary filename
          is chosen automatically. (default: ``None``)

        - ``verbosity`` -- a verbosity level, where zero means silent
          and anything else means verbose output. (default: ``0``)

        - ``**kwds`` -- accepted for compatibility with other solvers; ignored

        TESTS::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: DIMACS()
            DIMACS Solver: ''
        """
        self._headname_file_created_during_init = False
        if filename is None:
            filename = tmp_filename()
            self._headname_file_created_during_init = True

        self._headname = filename
        self._verbosity = verbosity

        if command is not None:
            self._command = command
        else:
            self._command = self.__class__.command

        self._tail = open(tmp_filename(), 'w')
        self._var = 0
        self._lit = 0

    def __repr__(self):
        """
        TESTS::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: DIMACS(command="iliketurtles {input}")
            DIMACS Solver: 'iliketurtles {input}'
        """
        return "DIMACS Solver: '%s'" % (self._command)

    def __del__(self):
        """
        TESTS::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: d = DIMACS(command="iliketurtles {input}")
            sage: del d

        We check that files created during initialization are properly
        deleted (:issue:`38328`)::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: d = DIMACS(command="iliketurtles {input}")
            sage: filename = d._headname
            sage: os.path.exists(filename)
            True
            sage: del d
            sage: os.path.exists(filename)
            False

        ::

            sage: fn = tmp_filename()
            sage: d = DIMACS(filename=fn)
            sage: del d
        """
        if not self._tail.closed:
            self._tail.close()
        if os.path.exists(self._tail.name):
            os.unlink(self._tail.name)
        if self._headname_file_created_during_init and os.path.exists(self._headname):
            os.unlink(self._headname)

    def var(self, decision=None):
        """
        Return a *new* variable.

        INPUT:

        - ``decision`` -- accepted for compatibility with other solvers; ignored

        EXAMPLES::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: solver = DIMACS()
            sage: solver.var()
            1
        """
        self._var += 1
        return self._var

    def nvars(self):
        """
        Return the number of variables.

        EXAMPLES::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: solver = DIMACS()
            sage: solver.var()
            1
            sage: solver.var(decision=True)
            2
            sage: solver.nvars()
            2
        """
        return self._var

    def add_clause(self, lits):
        """
        Add a new clause to set of clauses.

        INPUT:

        - ``lits`` -- tuple of nonzero integers

        .. NOTE::

            If any element ``e`` in ``lits`` has ``abs(e)`` greater
            than the number of variables generated so far, then new
            variables are created automatically.

        EXAMPLES::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: solver = DIMACS()
            sage: solver.var()
            1
            sage: solver.var(decision=True)
            2
            sage: solver.add_clause( (1, -2 , 3) )
            sage: solver
            DIMACS Solver: ''
        """
        l = []
        for lit in lits:
            lit = int(lit)
            while abs(lit) > self.nvars():
                self.var()
            l.append(str(lit))
        l.append("0\n")
        self._tail.write(" ".join(l) )
        self._lit += 1

    def write(self, filename=None):
        """
        Write DIMACS file.

        INPUT:

        - ``filename`` -- if ``None`` default filename specified at initialization is used for
          writing to (default: ``None``)

        EXAMPLES::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS(filename=fn)
            sage: solver.add_clause( (1, -2 , 3) )
            sage: _ = solver.write()
            sage: for line in open(fn).readlines():
            ....:     print(line)
            p cnf 3 1
            1 -2 3 0

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS()
            sage: solver.add_clause( (1, -2 , 3) )
            sage: _ = solver.write(fn)
            sage: for line in open(fn).readlines():
            ....:      print(line)
            p cnf 3 1
            1 -2 3 0
        """
        headname = self._headname if filename is None else filename
        head = open(headname, "w")
        head.truncate(0)
        head.write("p cnf %d %d\n" % (self._var,self._lit))
        head.close()

        tail = self._tail
        tail.close()

        head = open(headname,"a")
        tail = open(self._tail.name)
        head.write(tail.read())
        tail.close()
        head.close()

        self._tail = open(self._tail.name,"a")
        return headname

    def clauses(self, filename=None):
        """
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

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS()
            sage: solver.add_clause( (1, 2, 3) )
            sage: solver.clauses()
            [((1, 2, 3), False, None)]

            sage: solver.add_clause( (1, 2, -3) )
            sage: solver.clauses(fn)
            sage: print(open(fn).read())
            p cnf 3 2
            1 2 3 0
            1 2 -3 0
            <BLANKLINE>
        """
        if filename is not None:
            self.write(filename)
        else:
            tail = self._tail
            tail.close()
            tail = open(self._tail.name)

            clauses = []
            for line in tail.readlines():
                if line.startswith("p") or line.startswith("c"):
                    continue
                clause = []
                for lit in line.split(" "):
                    lit = int(lit)
                    if lit == 0:
                        break
                    clause.append(lit)
                clauses.append( ( tuple(clause), False, None ) )
            tail.close()
            self._tail = open(self._tail.name, "a")
            return clauses

    @staticmethod
    def render_dimacs(clauses, filename, nlits):
        """
        Produce DIMACS file ``filename`` from ``clauses``.

        INPUT:

        - ``clauses`` -- list of clauses, either in simple format as a list of
          literals or in extended format for CryptoMiniSat: a tuple of literals,
          ``is_xor`` and ``rhs``.

        - ``filename`` -- the file to write to

        - ``nlits -- the number of literals appearing in ``clauses``

        EXAMPLES::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS()
            sage: solver.add_clause( (1, 2, -3) )
            sage: DIMACS.render_dimacs(solver.clauses(), fn, solver.nvars())
            sage: print(open(fn).read())
            p cnf 3 1
            1 2 -3 0
            <BLANKLINE>

        This is equivalent to::

            sage: solver.clauses(fn)
            sage: print(open(fn).read())
            p cnf 3 1
            1 2 -3 0
            <BLANKLINE>

        This function also accepts a "simple" format::

            sage: DIMACS.render_dimacs([ (1,2), (1,2,-3) ], fn, 3)
            sage: print(open(fn).read())
            p cnf 3 2
            1 2 0
            1 2 -3 0
            <BLANKLINE>
        """
        fh = open(filename, "w")
        fh.write("p cnf %d %d\n" % (nlits,len(clauses)))
        for clause in clauses:
            if len(clause) == 3 and clause[1] in (True, False) and clause[2] in (True,False,None):
                lits, is_xor, rhs = clause
            else:
                lits, is_xor, rhs = clause, False, None

            if is_xor:
                closing = lits[-1] if rhs else -lits[-1]
                fh.write("x" + " ".join(map(str, lits[:-1])) + " %d 0\n" % closing)
            else:
                fh.write(" ".join(map(str, lits)) + " 0\n")
        fh.close()

    def _run(self):
        r"""
        Run 'command' and collect output.

        TESTS:

        This class is not meant to be called directly::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS(filename=fn)
            sage: solver.add_clause( (1, -2 , 3) )
            sage: solver._run()
            Traceback (most recent call last):
            ...
            ValueError: no SAT solver command selected

        It is used by subclasses::

            sage: from sage.sat.solvers import Glucose
            sage: solver = Glucose()
            sage: solver.add_clause( (1, 2, 3) )
            sage: solver.add_clause( (-1,) )
            sage: solver.add_clause( (-2,) )
            sage: solver._run()                       # optional - glucose
            sage: solver._output                      # optional - glucose
            [...
             's SATISFIABLE\n',
             'v -1 -2 3 0\n']
        """
        from sage.misc.verbose import get_verbose

        self.write()
        output_filename = None
        self._output = []

        command = self._command.strip()

        if not command:
            raise ValueError("no SAT solver command selected")

        if "{output}" in command:
            output_filename = tmp_filename()
        command = command.format(input=self._headname, output=output_filename)

        args = shlex.split(command)

        try:
            process = subprocess.Popen(args, stdout=subprocess.PIPE)
        except OSError:
            raise OSError("Could run '%s', perhaps you need to add your SAT solver to $PATH?" % (" ".join(args)))

        try:
            while process.poll() is None:
                for line in iter(process.stdout.readline, b''):
                    if get_verbose() or self._verbosity:
                        print(line)
                        sys.stdout.flush()
                    self._output.append(line.decode('utf-8'))
                sleep(0.1)
            if output_filename:
                self._output.extend(open(output_filename).readlines())
        except BaseException:
            process.kill()
            raise

    def __call__(self, assumptions=None):
        """
        Solve this instance and return the parsed output.

        INPUT:

        - ``assumptions`` -- ignored, accepted for compatibility with
          other solvers (default: ``None``)

        OUTPUT:

        - If this instance is SAT: A tuple of length ``nvars()+1``
          where the ``i``-th entry holds an assignment for the
          ``i``-th variables (the ``0``-th entry is always ``None``).

        - If this instance is UNSAT: ``False``

        EXAMPLES:

        When the problem is SAT::

            sage: from sage.sat.solvers import RSat
            sage: solver = RSat()
            sage: solver.add_clause( (1, 2, 3) )
            sage: solver.add_clause( (-1,) )
            sage: solver.add_clause( (-2,) )
            sage: solver()                            # optional - rsat
            (None, False, False, True)

        When the problem is UNSAT::

            sage: solver = RSat()
            sage: solver.add_clause((1,2))
            sage: solver.add_clause((-1,2))
            sage: solver.add_clause((1,-2))
            sage: solver.add_clause((-1,-2))
            sage: solver()                            # optional - rsat
            False

        With Glucose::

            sage: from sage.sat.solvers.dimacs import Glucose
            sage: solver = Glucose()
            sage: solver.add_clause((1,2))
            sage: solver.add_clause((-1,2))
            sage: solver.add_clause((1,-2))
            sage: solver()                           # optional - glucose
            (None, True, True)
            sage: solver.add_clause((-1,-2))
            sage: solver()                           # optional - glucose
            False

        With GlucoseSyrup::

            sage: from sage.sat.solvers.dimacs import GlucoseSyrup
            sage: solver = GlucoseSyrup()
            sage: solver.add_clause((1,2))
            sage: solver.add_clause((-1,2))
            sage: solver.add_clause((1,-2))
            sage: solver()                          # optional - glucose
            (None, True, True)
            sage: solver.add_clause((-1,-2))
            sage: solver()                          # optional - glucose
            False

        TESTS::

            sage: from sage.sat.boolean_polynomials import solve as solve_sat
            sage: sr = mq.SR(1, 1, 1, 4, gf2=True, polybori=True)                       # needs sage.rings.finite_rings sage.rings.polynomial.pbori
            sage: while True:  # workaround (see :issue:`31891`)                         # needs sage.rings.finite_rings sage.rings.polynomial.pbori
            ....:     try:
            ....:         F, s = sr.polynomial_system()
            ....:         break
            ....:     except ZeroDivisionError:
            ....:         pass
            sage: solve_sat(F, solver=sage.sat.solvers.RSat)    # optional - rsat, needs sage.rings.finite_rings sage.rings.polynomial.pbori
        """
        if assumptions is not None:
            raise NotImplementedError("Assumptions are not supported for DIMACS based solvers.")

        self._run()

        v_lines = []
        for line in self._output:
            if line.startswith("c"):
                continue
            if line.startswith("s"):
                if "UNSAT" in line:
                    return False
            if line.startswith("v"):
                v_lines.append(line[1:].strip())

        if v_lines:
            L = " ".join(v_lines).split(" ")
            assert L[-1] == "0", "last digit of solution line must be zero (not {})".format(L[-1])
            return (None,) + tuple(int(e) > 0 for e in L[:-1])
        else:
            raise ValueError("When parsing the output(={}), no line starts with letter v or s".format(self._output))


class RSat(DIMACS):
    """
    An instance of the RSat solver.

    For information on RSat see: http://reasoning.cs.ucla.edu/rsat/

    EXAMPLES::

        sage: from sage.sat.solvers import RSat
        sage: solver = RSat()
        sage: solver
        DIMACS Solver: 'rsat {input} -v -s'

    When the problem is SAT::

        sage: from sage.sat.solvers import RSat
        sage: solver = RSat()
        sage: solver.add_clause( (1, 2, 3) )
        sage: solver.add_clause( (-1,) )
        sage: solver.add_clause( (-2,) )
        sage: solver()                            # optional - rsat
        (None, False, False, True)

    When the problem is UNSAT::

        sage: solver = RSat()
        sage: solver.add_clause((1,2))
        sage: solver.add_clause((-1,2))
        sage: solver.add_clause((1,-2))
        sage: solver.add_clause((-1,-2))
        sage: solver()                            # optional - rsat
        False
    """
    command = "rsat {input} -v -s"


class Glucose(DIMACS):
    """
    An instance of the Glucose solver.

    For information on Glucose see: http://www.labri.fr/perso/lsimon/glucose/

    EXAMPLES::

        sage: from sage.sat.solvers import Glucose
        sage: solver = Glucose()
        sage: solver
        DIMACS Solver: 'glucose -verb=0 -model {input}'

    When the problem is SAT::

        sage: from sage.sat.solvers import Glucose
        sage: solver1 = Glucose()
        sage: solver1.add_clause( (1, 2, 3) )
        sage: solver1.add_clause( (-1,) )
        sage: solver1.add_clause( (-2,) )
        sage: solver1()                            # optional - glucose
        (None, False, False, True)

    When the problem is UNSAT::

        sage: solver2 = Glucose()
        sage: solver2.add_clause((1,2))
        sage: solver2.add_clause((-1,2))
        sage: solver2.add_clause((1,-2))
        sage: solver2.add_clause((-1,-2))
        sage: solver2()                            # optional - glucose
        False

    With one hundred variables::

        sage: solver3 = Glucose()
        sage: solver3.add_clause( (1, 2, 100) )
        sage: solver3.add_clause( (-1,) )
        sage: solver3.add_clause( (-2,) )
        sage: solver3()                            # optional - glucose
        (None, False, False, ..., True)

    TESTS::

        sage: print(''.join(solver1._output))      # optional - glucose
        c...
        s SATISFIABLE
        v -1 -2 3 0

    ::

        sage: print(''.join(solver2._output))      # optional - glucose
        c...
        s UNSATISFIABLE

    Glucose gives large solution on one single line::

        sage: print(''.join(solver3._output))      # optional - glucose
        c...
        s SATISFIABLE
        v -1 -2 ... 100 0
    """
    command = "glucose -verb=0 -model {input}"


class GlucoseSyrup(DIMACS):
    """
    An instance of the Glucose-syrup parallel solver.

    For information on Glucose see: http://www.labri.fr/perso/lsimon/glucose/

    EXAMPLES::

        sage: from sage.sat.solvers import GlucoseSyrup
        sage: solver = GlucoseSyrup()
        sage: solver
        DIMACS Solver: 'glucose-syrup -model -verb=0 {input}'

    When the problem is SAT::

        sage: solver1 = GlucoseSyrup()
        sage: solver1.add_clause( (1, 2, 3) )
        sage: solver1.add_clause( (-1,) )
        sage: solver1.add_clause( (-2,) )
        sage: solver1()                            # optional - glucose
        (None, False, False, True)

    When the problem is UNSAT::

        sage: solver2 = GlucoseSyrup()
        sage: solver2.add_clause((1,2))
        sage: solver2.add_clause((-1,2))
        sage: solver2.add_clause((1,-2))
        sage: solver2.add_clause((-1,-2))
        sage: solver2()                            # optional - glucose
        False

    With one hundred variables::

        sage: solver3 = GlucoseSyrup()
        sage: solver3.add_clause( (1, 2, 100) )
        sage: solver3.add_clause( (-1,) )
        sage: solver3.add_clause( (-2,) )
        sage: solver3()                            # optional - glucose
        (None, False, False, ..., True)

    TESTS::

        sage: print(''.join(solver1._output))      # optional - glucose
        c...
        s SATISFIABLE
        v -1 -2 3 0

    ::

        sage: print(''.join(solver2._output))      # optional - glucose
        c...
        s UNSATISFIABLE

    GlucoseSyrup gives large solution on one single line::

        sage: print(''.join(solver3._output))      # optional - glucose
        c...
        s SATISFIABLE
        v -1 -2 ... 100 0
    """
    command = "glucose-syrup -model -verb=0 {input}"


class Kissat(DIMACS):
    """
    An instance of the Kissat SAT solver.

    For information on Kissat see: http://fmv.jku.at/kissat/

    EXAMPLES::

        sage: from sage.sat.solvers import Kissat
        sage: solver = Kissat()
        sage: solver
        DIMACS Solver: 'kissat -q {input}'

    When the problem is SAT::

        sage: solver1 = Kissat()
        sage: solver1.add_clause( (1, 2, 3) )
        sage: solver1.add_clause( (-1,) )
        sage: solver1.add_clause( (-2,) )
        sage: solver1()                           # optional - kissat
        (None, False, False, True)

    When the problem is UNSAT::

        sage: solver2 = Kissat()
        sage: solver2.add_clause((1,2))
        sage: solver2.add_clause((-1,2))
        sage: solver2.add_clause((1,-2))
        sage: solver2.add_clause((-1,-2))
        sage: solver2()                           # optional - kissat
        False

    With one hundred variables::

        sage: solver3 = Kissat()
        sage: solver3.add_clause( (1, 2, 100) )
        sage: solver3.add_clause( (-1,) )
        sage: solver3.add_clause( (-2,) )
        sage: solver3()                           # optional - kissat
        (None, False, False, ..., True)

    TESTS::

        sage: print(''.join(solver1._output))     # optional - kissat
        s SATISFIABLE
        v -1 -2 3 0

    ::

        sage: print(''.join(solver2._output))     # optional - kissat
        s UNSATISFIABLE

    Here the output contains many lines starting with letter "v"::

        sage: print(''.join(solver3._output))     # optional - kissat
        s SATISFIABLE
        v -1 -2 ...
        v ...
        v ...
        v ... 100 0
    """

    command = "kissat -q {input}"
