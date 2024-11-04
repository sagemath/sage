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

By default, Sage solves SAT instances as an Integer Linear Program (see
:mod:`sage.numerical.mip`), but any SAT solver supporting the DIMACS input
format is easily interfaced using the :class:`sage.sat.solvers.dimacs.DIMACS`
blueprint. Sage ships with pre-written interfaces for *RSat* [RS]_ and *Glucose*
[GL]_. Furthermore, Sage provides an interface to the *CryptoMiniSat* [CMS]_ SAT
solver which can be used interchangeably with DIMACS-based solvers. For this last
solver, the optional CryptoMiniSat package must be installed, this can be
accomplished by typing the following in the shell::

    sage -i cryptominisat sagelib

We now show how to solve a simple SAT problem. ::

    (x1 OR x2 OR x3) AND (x1 OR x2 OR (NOT x3))

In Sage's notation::

    sage: solver = SAT()
    sage: solver.add_clause( ( 1,  2,  3) )
    sage: solver.add_clause( ( 1,  2, -3) )
    sage: solver()       # random
    (None, True, True, False)

.. NOTE::

    :meth:`~sage.sat.solvers.dimacs.DIMACS.add_clause` creates new variables
    when necessary. When using CryptoMiniSat, it creates *all* variables up to
    the given index. Hence, adding a literal involving the variable 1000 creates
    up to 1000 internal variables.

DIMACS-base solvers can also be used to write DIMACS files::

    sage: from sage.sat.solvers.dimacs import DIMACS
    sage: fn = tmp_filename()
    sage: solver = DIMACS(filename=fn)
    sage: solver.add_clause( ( 1,  2,  3) )
    sage: solver.add_clause( ( 1,  2, -3) )
    sage: _ = solver.write()
    sage: for line in open(fn).readlines():
    ....:    print(line)
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
    ....:    print(line)
    p cnf 3 2
    1 2 3 0
    1 2 -3 0

These files can then be passed external SAT solvers.

Details on Specific Solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1

   sage/sat/solvers/satsolver
   sage/sat/solvers/dimacs
   sage/sat/solvers/picosat
   sage/sat/solvers/sat_lp
   sage/sat/solvers/cryptominisat

Converters
----------

Sage supports conversion from Boolean polynomials (also known as Algebraic Normal Form) to
Conjunctive Normal Form::

    sage: # needs sage.rings.polynomial.pbori
    sage: B.<a,b,c> = BooleanPolynomialRing()
    sage: from sage.sat.converters.polybori import CNFEncoder
    sage: from sage.sat.solvers.dimacs import DIMACS
    sage: fn = tmp_filename()
    sage: solver = DIMACS(filename=fn)
    sage: e = CNFEncoder(solver, B)
    sage: e.clauses_sparse(a*b + a + 1)
    sage: _ = solver.write()
    sage: print(open(fn).read())
    p cnf 3 2
    -2 0
    1 0
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

    sage: # needs sage.rings.polynomial.pbori
    sage: sr = mq.SR(1,1,1,4,gf2=True,polybori=True)
    sage: while True:
    ....:     try:
    ....:         F,s = sr.polynomial_system()
    ....:         break
    ....:     except ZeroDivisionError:
    ....:         pass
    sage: from sage.sat.boolean_polynomials import solve as solve_sat # optional - pycryptosat
    sage: s = solve_sat(F)                                            # optional - pycryptosat
    sage: F.subs(s[0])                                                # optional - pycryptosat
    Polynomial Sequence with 36 Polynomials in 0 Variables

Details on Specific Highlevel Interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   sage/sat/boolean_polynomials


REFERENCES:

.. [RS] http://reasoning.cs.ucla.edu/rsat/

.. [GL] http://www.lri.fr/~simon/?page=glucose

.. [CMS] http://www.msoos.org

.. [SG09] http://www.satcompetition.org/2009/format-benchmarks2009.html

.. include:: ../footer.txt
