Sat
===

Sage supports solving clauses in Conjunctive Normal Form (see :wikipedia:`Conjunctive_normal_form`),
i.e., SAT solving, via an interface inspired by the usual DIMACS format used in SAT solving
[SG09]_. For example, to express that::

   x1 OR x2 OR (NOT x3)

should be true, we write::

   (1, 2, -3)

.. WARNING::

    Variable indices **must** start at one.

Solvers
-------

Any SAT solver supporting the DIMACS input format is easily interfaced using the
:class:`sage.sat.solvers.dimacs.DIMACS` blueprint. Sage ships with pre-written interfaces for *RSat*
[RS]_ and *Glucose* [GL]_. Furthermore, Sage provides a C++ interface to the *CryptoMiniSat* [CMS]_ SAT
solver which can be used interchangably with DIMACS-based solvers, but also provides advanced
features. For this, the optional CryptoMiniSat package must be installed, this can be accomplished
by typing::

    sage: install_package('cryptominisat') # not tested

and by running ``sage -b`` from the shell afterwards to build Sage's CryptoMiniSat extension module.

Since by default Sage does not include any SAT solver, we demonstrate key features by instantiating
a fake DIMACS-based solver. We start with a trivial example::

    (x1 OR x2 OR x3) AND (x1 OR x2 OR (NOT x3))

In Sage's notation::

    sage: from sage.sat.solvers.dimacs import DIMACS
    sage: solver = DIMACS(command="sat-solver")
    sage: solver.add_clause( ( 1,  2,  3) )
    sage: solver.add_clause( ( 1,  2, -3) )

.. NOTE::

    :meth:`sage.sat.solvers.dimacs.DIMACS.add_clause` creates new variables when necessary. In
    particular, it creates *all* variables up to the given index. Hence, adding a literal involving
    the variable 1000 creates up to 1000 internal variables.

DIMACS-base solvers can also be used to write DIMACS files::

    sage: from sage.sat.solvers.dimacs import DIMACS
    sage: fn = tmp_filename()
    sage: solver = DIMACS(filename=fn)
    sage: solver.add_clause( ( 1,  2,  3) )
    sage: solver.add_clause( ( 1,  2, -3) )
    sage: _ = solver.write()
    sage: for line in open(fn).readlines():
    ....:    print line,
    p cnf 3 2
    1 2 3 0
    1 2 -3 0

Alternatively, there is :meth:`sage.sat.solvers.dimacs.DIMACS.clauses`::

    sage: from sage.sat.solvers.dimacs import DIMACS
    sage: fn = tmp_filename()
    sage: solver = DIMACS()
    sage: solver.add_clause( ( 1,  2,  3) )
    sage: solver.add_clause( ( 1,  2, -3) )
    sage: solver.clauses(fn)
    sage: for line in open(fn).readlines():
    ....:    print line,
    p cnf 3 2
    1 2 3 0
    1 2 -3 0

These files can then be passed external SAT solvers.

We demonstrate solving using CryptoMiniSat::

    sage: from sage.sat.solvers import CryptoMiniSat # optional - cryptominisat
    sage: cms = CryptoMiniSat()                      # optional - cryptominisat
    sage: cms.add_clause((1,2,-3))                   # optional - cryptominisat
    sage: cms()                                      # optional - cryptominisat
    (None, True, True, False)

Details on Specific Solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   sage/sat/solvers/satsolver
   sage/sat/solvers/dimacs
.. optional - cryptominisat
.. sage/sat/solvers/cryptominisat/cryptominisat
.. sage/sat/solvers/cryptominisat/solverconf

Converters
----------

Sage supports conversion from Boolean polynomials (also known as Algebraic Normal Form) to
Conjunctive Normal Form::

    sage: B.<a,b,c> = BooleanPolynomialRing()
    sage: from sage.sat.converters.polybori import CNFEncoder
    sage: from sage.sat.solvers.dimacs import DIMACS
    sage: fn = tmp_filename()
    sage: solver = DIMACS(filename=fn)
    sage: e = CNFEncoder(solver, B)
    sage: e.clauses_sparse(a*b + a + 1)
    sage: _ = solver.write()
    sage: print open(fn).read()
    p cnf 3 2
    1 0
    -2 0
    <BLANKLINE>

Details on Specific Converterts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   sage/sat/converters/polybori

Highlevel Interfaces
--------------------

Sage provides various highlevel functions which make working with Boolean polynomials easier. We
construct a very small-scale AES system of equations and pass it to a SAT solver::

    sage: sr = mq.SR(1,1,1,4,gf2=True,polybori=True)
    sage: F,s = sr.polynomial_system()
    sage: from sage.sat.boolean_polynomials import solve as solve_sat # optional - cryptominisat
    sage: s = solve_sat(F)                                            # optional - cryptominisat
    sage: F.subs(s[0])                                                # optional - cryptominisat
    Polynomial Sequence with 36 Polynomials in 0 Variables

Details on Specific Highlevel Interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   sage/sat/boolean_polynomials


REFERENCES:

.. [RS] http://reasoning.cs.ucla.edu/rsat/

.. [GL] http://www.lri.fr/~simon/?page=glucose

.. [CMS] http://www.msoos.org/cryptominisat2/

.. [SG09] http://www.satcompetition.org/2009/format-benchmarks2009.html

.. include:: ../footer.txt
