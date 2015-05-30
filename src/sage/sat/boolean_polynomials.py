"""
SAT Functions for Boolean Polynomials

These highlevel functions support solving and learning from Boolean polynomial systems. In this
context, "learning" means the construction of new polynomials in the ideal spanned by the original
polynomials.

AUTHOR:

- Martin Albrecht (2012): initial version

Functions
^^^^^^^^^
"""
##############################################################################
#  Copyright (C) 2012 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.sat.solvers import SatSolver
from sage.sat.converters import ANF2CNFConverter

from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence


def solve(F, converter=None, solver=None, n=1, target_variables=None, **kwds):
    """
    Solve system of Boolean polynomials ``F`` by solving the
    SAT-problem -- produced by ``converter`` -- using ``solver``.

    INPUT:

    - ``F`` - a sequence of Boolean polynomials

    - ``n`` - number of solutions to return. If ``n`` is +infinity
      then all solutions are returned. If ``n <infinity`` then ``n``
      solutions are returned if ``F`` has at least ``n``
      solutions. Otherwise, all solutions of ``F`` are
      returned. (default: ``1``)

    - ``converter`` - an ANF to CNF converter class or object.  If
      ``converter`` is ``None`` then
      :class:`sage.sat.converters.polybori.CNFEncoder` is used to
      construct a new converter. (default: ``None``)

    - ``solver`` - a SAT-solver class or object. If ``solver`` is
      ``None`` then :class:`sage.sat.solvers.cryptominisat.CryptoMiniSat`
      is used to construct a new converter.  (default: ``None``)

    - ``target_variables`` - a list of variables. The elements of the list are
      used to exclude a particular combination of variable assignments of a
      solution from any further solution. Furthermore ``target_variables``
      denotes which variable-value pairs appear in the solutions. If
      ``target_variables`` is ``None`` all variables appearing in the
      polynomials of ``F`` are used to construct exclusion clauses.
      (default: ``None``)

    - ``**kwds`` - parameters can be passed to the converter and the
       solver by prefixing them with ``c_`` and ``s_`` respectively. For
       example, to increase CryptoMiniSat's verbosity level, pass
       ``s_verbosity=1``.

    OUTPUT:

        A list of dictionaries, each of which contains a variable
        assignment solving ``F``.

    EXAMPLE:

    We construct a very small-scale AES system of equations::

        sage: sr = mq.SR(1,1,1,4,gf2=True,polybori=True)
        sage: F,s = sr.polynomial_system()

    and pass it to a SAT solver::

        sage: from sage.sat.boolean_polynomials import solve as solve_sat # optional - cryptominisat
        sage: s = solve_sat(F)                                            # optional - cryptominisat
        sage: F.subs(s[0])                                                # optional - cryptominisat
        Polynomial Sequence with 36 Polynomials in 0 Variables

    This time we pass a few options through to the converter and the solver::

        sage: s = solve_sat(F, s_verbosity=1, c_max_vars_sparse=4, c_cutting_number=8) # optional - cryptominisat
        c Flit...
        ...
        sage: F.subs(s[0])                                                             # optional - cryptominisat
        Polynomial Sequence with 36 Polynomials in 0 Variables

    We construct a very simple system with three solutions and ask for a specific number of solutions::

        sage: B.<a,b> = BooleanPolynomialRing() # optional - cryptominisat
        sage: f = a*b                           # optional - cryptominisat
        sage: l = solve_sat([f],n=1)            # optional - cryptominisat
        sage: len(l) == 1, f.subs(l[0])         # optional - cryptominisat
        (True, 0)

        sage: l = sorted(solve_sat([a*b],n=2))        # optional - cryptominisat
        sage: len(l) == 2, f.subs(l[0]), f.subs(l[1]) # optional - cryptominisat
        (True, 0, 0)

        sage: sorted(solve_sat([a*b],n=3))         # optional - cryptominisat
        [{b: 0, a: 0}, {b: 0, a: 1}, {b: 1, a: 0}]
        sage: sorted(solve_sat([a*b],n=4))         # optional - cryptominisat
        [{b: 0, a: 0}, {b: 0, a: 1}, {b: 1, a: 0}]
        sage: sorted(solve_sat([a*b],n=infinity))  # optional - cryptominisat
        [{b: 0, a: 0}, {b: 0, a: 1}, {b: 1, a: 0}]

    In the next example we see how the ``target_variables`` parameter works::

        sage: from sage.sat.boolean_polynomials import solve as solve_sat # optional - cryptominisat
        sage: R.<a,b,c,d> = BooleanPolynomialRing()                       # optional - cryptominisat
        sage: F = [a+b,a+c+d]                                             # optional - cryptominisat

    First the normal use case::

        sage: sorted(solve_sat(F,n=infinity))                             # optional - cryptominisat
        [{d: 0, c: 0, b: 0, a: 0},
         {d: 0, c: 1, b: 1, a: 1},
         {d: 1, c: 0, b: 1, a: 1},
         {d: 1, c: 1, b: 0, a: 0}]

    Now we are only interested in the solutions of the variables a and b::

        sage: solve_sat(F,n=infinity,target_variables=[a,b])              # optional - cryptominisat
        [{b: 0, a: 0}, {b: 1, a: 1}]

    .. NOTE::

       Although supported, passing converter and solver objects
       instead of classes is discouraged because these objects are
       stateful.
    """
    assert(n>0)

    try:
        len(F)
    except AttributeError:
        F = F.gens()
        len(F)

    P = next(iter(F)).parent()
    K = P.base_ring()

    if target_variables is None:
        target_variables = PolynomialSequence(F).variables()
    else:
        target_variables = PolynomialSequence(target_variables).variables()
        assert(set(target_variables).issubset(set(P.gens())))

    # instantiate the SAT solver

    if solver is None:
        from sage.sat.solvers.cryptominisat import CryptoMiniSat as solver

    if not isinstance(solver, SatSolver):
        solver_kwds = {}
        for k, v in kwds.iteritems():
            if k.startswith("s_"):
                solver_kwds[k[2:]] = v

        solver = solver(**solver_kwds)

    # instantiate the ANF to CNF converter

    if converter is None:
        from sage.sat.converters.polybori import CNFEncoder as converter

    if not isinstance(converter, ANF2CNFConverter):
        converter_kwds = {}
        for k, v in kwds.iteritems():
            if k.startswith("c_"):
                converter_kwds[k[2:]] = v

        converter = converter(solver, P, **converter_kwds)

    phi = converter(F)
    rho = dict((phi[i], i) for i in range(len(phi)))

    S = []

    while True:
        s = solver()

        if s:
            S.append(dict((x, K(s[rho[x]])) for x in target_variables))

            if n is not None and len(S) == n:
                break

            exclude_solution = tuple(-rho[x] if s[rho[x]] else rho[x] for x in target_variables)
            solver.add_clause(exclude_solution)

        else:
            try:
                learnt = solver.learnt_clauses(unitary_only=True)
                if learnt:
                    S.append(dict((phi[abs(i)-1], K(i<0)) for i in learnt))
                else:
                    S.append(s)
                    break
            except (AttributeError, NotImplementedError):
                # solver does not support recovering learnt clauses
                S.append(s)
                break

    if len(S) == 1:
        if S[0] is False:
            return False
        if S[0] is None:
            return None
    elif S[-1] is False:
            return S[0:-1]
    return S


def learn(F, converter=None, solver=None, max_learnt_length=3, interreduction=False, **kwds):
    """
    Learn new polynomials by running SAT-solver ``solver`` on
    SAT-instance produced by ``converter`` from ``F``.

    INPUT:

    - ``F`` - a sequence of Boolean polynomials

    - ``converter`` - an ANF to CNF converter class or object.  If ``converter`` is ``None`` then
      :class:`sage.sat.converters.polybori.CNFEncoder` is used to construct a new
      converter. (default: ``None``)

    - ``solver`` - a SAT-solver class or object. If ``solver`` is ``None`` then
      :class:`sage.sat.solvers.cryptominisat.CryptoMiniSat` is used to construct a new converter.
      (default: ``None``)

    - ``max_learnt_length`` - only clauses of length <= ``max_length_learnt`` are considered and
      converted to polynomials. (default: ``3``)

    - ``interreduction`` - inter-reduce the resulting polynomials (default: ``False``)

    .. NOTE::

       More parameters can be passed to the converter and the solver by prefixing them with ``c_`` and
       ``s_`` respectively. For example, to increase CryptoMiniSat's verbosity level, pass
       ``s_verbosity=1``.

    OUTPUT:

        A sequence of Boolean polynomials.

    EXAMPLE::

       sage: from sage.sat.boolean_polynomials import learn as learn_sat # optional - cryptominisat

    We construct a simple system and solve it::

       sage: set_random_seed(2300)                      # optional - cryptominisat
       sage: sr = mq.SR(1,2,2,4,gf2=True,polybori=True) # optional - cryptominisat
       sage: F,s = sr.polynomial_system()               # optional - cryptominisat
       sage: H = learn_sat(F)                           # optional - cryptominisat
       sage: H[-1]                                      # optional - cryptominisat
       k033 + 1

    We construct a slightly larger equation system and recover some
    equations after 20 restarts::

       sage: set_random_seed(2303)                        # optional - cryptominisat
       sage: sr = mq.SR(1,4,4,4,gf2=True,polybori=True)   # optional - cryptominisat
       sage: F,s = sr.polynomial_system()                 # optional - cryptominisat
       sage: from sage.sat.boolean_polynomials import learn as learn_sat # optional - cryptominisat
       sage: H = learn_sat(F, s_maxrestarts=20, interreduction=True)     # optional - cryptominisat
       sage: H[-1]                                        # optional - cryptominisat, output random
       k001200*s031*x011201 + k001200*x011201

    .. NOTE::

       This function is meant to be called with some parameter such
       that the SAT-solver is interrupted. For CryptoMiniSat this is
       max_restarts, so pass 'c_max_restarts' to limit the number of
       restarts CryptoMiniSat will attempt. If no such parameter is
       passed, then this function behaves essentially like
       :func:`solve` except that this function does not support
       ``n>1``.

    TESTS:

    We test that :trac:`17351` is fixed, by checking that the following doctest does not raise an
    error::

        sage: P.<a,b,c> = BooleanPolynomialRing()
        sage: F = [a*c + a + b*c + c + 1,  a*b + a*c + a + c + 1,  a*b + a*c + a + b*c + 1]
        sage: from sage.sat.boolean_polynomials import learn as learn_sat # optional - cryptominisat
        sage: learn_sat(F, s_maxrestarts=0, interreduction=True)      # optional - cryptominisat
        []
    """
    try:
        len(F)
    except AttributeError:
        F = F.gens()
        len(F)

    P = next(iter(F)).parent()
    K = P.base_ring()

    # instantiate the SAT solver

    if solver is None:
        from sage.sat.solvers.cryptominisat import CryptoMiniSat as solver

    solver_kwds = {}
    for k, v in kwds.iteritems():
        if k.startswith("s_"):
            solver_kwds[k[2:]] = v

    solver = solver(**solver_kwds)

    # instantiate the ANF to CNF converter

    if converter is None:
        from sage.sat.converters.polybori import CNFEncoder as converter

    converter_kwds = {}
    for k, v in kwds.iteritems():
        if k.startswith("c_"):
            converter_kwds[k[2:]] = v

    converter = converter(solver, P, **converter_kwds)

    phi = converter(F)
    rho = dict((phi[i], i) for i in range(len(phi)))

    s = solver()

    if s:
        learnt = [x + K(s[rho[x]]) for x in P.gens()]
    else:
        learnt = []
        for c in solver.learnt_clauses():
            if len(c) <= max_learnt_length:
                try:
                    learnt.append(converter.to_polynomial(c))
                except (ValueError, NotImplementedError, AttributeError):
                    # the solver might have learnt clauses that contain CNF
                    # variables which have no correspondence to variables in our
                    # polynomial ring (XOR chaining variables for example)
                    pass

    learnt = PolynomialSequence(P, learnt)

    if interreduction:
        learnt = learnt.ideal().interreduced_basis()
    return learnt
