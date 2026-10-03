r"""
Regression tests for certified BSD computations

These tests check the evidence, not agreement with the old proof engine.

Analytic assumptions cannot populate the certified cache::

    sage: from unittest.mock import patch
    sage: from sage.schemes.elliptic_curves import bsd_certificates as cert
    sage: from sage.schemes.elliptic_curves import BSD as bsd
    sage: from sage.schemes.elliptic_curves.sha_tate import Sha
    sage: E = EllipticCurve('11a1')
    sage: with patch.object(Sha, 'an', return_value=9):
    ....:     assumed = cert._AnalyticData(E, False)
    ....:     assert assumed.sha(E) == 9
    ....:     proven = cert._AnalyticData(E, True)
    ....:     assert proven.sha(E) == 1
    sage: assumed.certificates[E.a_invariants()]['certified']
    False
    sage: proven.certificates[E.a_invariants()]['certified']
    True
    sage: E.prove_BSD(proof=True, rigorous=False, return_BSD=True).assumptions
    ['analytic_rank() and sha().an() are assumed correct']
    sage: E.prove_BSD(proof=True, return_BSD=True).assumptions
    []

An unsuccessful derivative test is inconclusive at every precision::

    sage: E = EllipticCurve('37a1')
    sage: with patch.object(cert, '_derivative_ball', return_value=RBF(0, rad=1)) as f:
    ....:     a = cert._AnalyticData(E, True)
    ....:     assert a.rank is None
    ....:     assert [call.args[1] for call in f.call_args_list] == [64, 128, 256, 512]
    sage: derivative = cert._derivative_ball
    sage: def initially_ambiguous(E, prec):
    ....:     return RBF(0, rad=1) if prec == 64 else derivative(E, prec)
    sage: with patch.object(cert, '_derivative_ball', side_effect=initially_ambiguous):
    ....:     a = cert._AnalyticData(E, True)
    ....:     assert a.rank == 1 and 128 in a._derivatives

Refine an ambiguous reconstruction rather than rounding its midpoint::

    sage: unique = cert._unique_integer
    sage: with patch.object(cert, '_unique_integer', side_effect=[None, 4]) as f:
    ....:     a = cert._AnalyticData(E, True)
    ....:     assert a.sha(E) == 1
    ....:     assert f.call_count == 2 and 128 in a._derivatives
    sage: with patch.object(cert, '_unique_integer', return_value=3):
    ....:     cert._AnalyticData(E, True).sha(E)
    Traceback (most recent call last):
    ...
    ArithmeticError: Gross-Zagier reconstruction is not a positive square

Reconstruct a nontrivial rank-one analytic order.  This example needs more
than the default PARI stack limit on some installations; restore both stack
settings afterwards::

    sage: from sage.libs.pari import pari
    sage: stack = pari.stacksize(), pari.stacksizemax()  # long time
    sage: try:  # long time
    ....:     pari.allocatemem(2^28, max(stack[1], 2^32), silent=True)
    ....:     F = EllipticCurve([1,1,1,508,-2551])
    ....:     a = cert._AnalyticData(F, True)
    ....:     assert a.rank == 1 and a.sha(F) == 4
    ....:     witness = a.certificates[F.a_invariants()]
    ....:     assert witness['index_square'] == 256 and witness['lattice'] == 1/64
    ....: finally:
    ....:     pari.allocatemem(*stack, silent=True)

A supplied point need not already be a generator.  Full saturation proves
the basis used by both height computations::

    sage: from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field
    sage: P = E(0,0)
    sage: with patch.object(EllipticCurve_rational_field, 'gens', return_value=[3*P]):
    ....:     G = cert._AnalyticData(E, True).gens(E)
    sage: G[0] in [P, -P]
    True

The original proof flag still returns immediately, before any analytic work::

    sage: with patch.object(cert, '_AnalyticData', side_effect=AssertionError('called')):
    ....:     assert E.prove_BSD(proof=False, return_BSD=True) == []
    ....:     with proof.WithProof('elliptic_curve', False):
    ....:         assert E.prove_BSD() == []

Each selected descent backend runs.  Native two-isogeny descent records its
inapplicability when rational 2-torsion is absent; PARI then supplements it::

    sage: for backend in ['mwrank', 'pari', 'sage']:
    ....:     B = EllipticCurve('14a1').prove_BSD(proof=True, two_desc=backend, return_BSD=True)
    ....:     assert backend in B.descent and B.primes == []
    sage: B = E.prove_BSD(proof=True, two_desc='sage', return_BSD=True)
    sage: 'inconclusive' in B.descent['sage'] and 'pari' in B.descent
    True
    sage: B = bsd.BSD_data()
    sage: F = EllipticCurve('66b3')
    sage: bsd._two_primary(F, cert._AnalyticData(F, True), 'pari', B)
    True
    sage: B.descent['pari']['sha2_dimension'], B.descent['pari']['cassels_tate_dimension']
    (2, 2)
    sage: B = bsd.BSD_data()
    sage: bsd._two_primary(F, cert._AnalyticData(F, True), 'mwrank', B)
    True
    sage: sorted(B.descent)
    ['mwrank', 'pari']
    sage: F = EllipticCurve('210e7')
    sage: B = bsd.BSD_data()
    sage: bsd._two_primary(F, cert._AnalyticData(F, True), 'pari', B)
    False
    sage: B.bounds[2]
    (4, +Infinity)
    sage: a = cert._AnalyticData(E, True)
    sage: with patch.object(a, 'sha', return_value=QQ(4)):
    ....:     bsd._two_primary(E, a, 'pari', bsd.BSD_data())
    Traceback (most recent call last):
    ...
    ArithmeticError: Cassels-Tate pairing contradicts the analytic 2-part

Skinner's theorem records an actual auxiliary multiplicative prime::

    sage: E = EllipticCurve('5389a1')
    sage: B = E.prove_BSD(proof=True, return_BSD=True)
    sage: B.primes, B.sha_an, B.proof[3]['reference']
    ([], 9, 'Skinner2016, Theorem C')
    sage: q = B.proof[3]['witnesses']['q']
    sage: q != 3 and E.has_multiplicative_reduction(q) and E.discriminant().valuation(q) % 3 != 0
    True

BCS is dispatched and its rank and reduction restrictions are enforced::

    sage: E = EllipticCurve('37a1')
    sage: jetchev = bsd._jetchev_exceptions
    sage: def extra_exception(E, a):
    ....:     primes, witness = jetchev(E, a)
    ....:     return primes | {ZZ(5)}, witness
    sage: with patch.object(bsd, '_jetchev_exceptions', side_effect=extra_exception):
    ....:     B = E.prove_BSD(proof=True, return_BSD=True)
    sage: B.proof[5]['reference']
    'BCS2025, Corollary 1.3.1'
    sage: E = EllipticCurve('50b1')
    sage: a = cert._AnalyticData(E, True)
    sage: bsd._bcs(E, 5, a) is None and bsd._cgs(E, 5, a) is None
    True
    sage: bsd._castella(E, 5, a) is None
    True
    sage: E = EllipticCurve('389a1')
    sage: bsd._bcs(E, 5, cert._AnalyticData(E, False)) is None
    True

Rational torsion and a trivial dual local character both obstruct CGS::

    sage: for E in EllipticCurve('91b1').isogeny_class():
    ....:     assert bsd._cgs(E, 3, cert._AnalyticData(E, True)) is None
    sage: E = EllipticCurve('50b1')
    sage: bsd._cgs(E, 3, cert._AnalyticData(E, True))['reference']
    'CGS2025, Theorem D'
    sage: E = EllipticCurve('30a1')
    sage: E.discriminant().valuation(3) % 3 == 0
    True
    sage: bsd._multiplicative_local_torsion(E, 3)
    True
    sage: E = EllipticCurve('82a1')
    sage: bsd._castella(E, 41, cert._AnalyticData(E, True))['reference']
    "Castella2018err, Theorem A'"

The Iwasawa fallback compares its bound with the analytic valuation and
does not accept a reducibility/surjectivity guard with the logic reversed::

    sage: E = EllipticCurve('11a3')
    sage: bsd._iwasawa(E, 5, cert._AnalyticData(E, True))['witnesses']['bound']
    0
    sage: a = cert._AnalyticData(E, False)
    sage: with patch.object(a, 'sha', return_value=QQ(25)):
    ....:     assert bsd._iwasawa(E, 5, a) is None
    sage: E = EllipticCurve('37a1')
    sage: rho = E.galois_representation()
    sage: with patch.object(type(rho), 'is_surjective', return_value=False):
    ....:     with patch.object(type(rho), 'is_reducible', return_value=False):
    ....:         assert bsd._iwasawa(E, 5, cert._AnalyticData(E, True)) is None
    sage: E = EllipticCurve('53a1')
    sage: cert._padic_sha_bound(E, 5, cert._AnalyticData(E, True))  # long time
    0

Mod-3 surjectivity alone is not used as a Tate-module certificate::

    sage: E = EllipticCurve('121a1')
    sage: E.galois_representation().is_surjective(3)
    True
    sage: bsd._tate_module_surjectivity(E, 3) is None
    True
    sage: B = E.prove_BSD(proof=True, return_BSD=True)
    sage: 3 in B.proof['outside']['witnesses']['exceptional_primes']
    True

Deterministic choice and ordering, CM, and nonminimal models::

    sage: E = EllipticCurve('32a1')
    sage: B = E.prove_BSD(proof=True, return_BSD=True)
    sage: B.primes, B.proof['outside']['reference']
    ([], 'Rub1991, main theorem')
    sage: E = EllipticCurve('11a1')
    sage: B = E.prove_BSD(proof=True, return_BSD=True)
    sage: B.curve.a_invariants() == min(C.a_invariants() for C in E.isogeny_class())
    True
    sage: E.change_weierstrass_model([1/2,0,0,0]).prove_BSD(proof=True) == B.primes
    True
    sage: E = EllipticCurve('50b1')
    sage: E.prove_BSD(proof=True)
    [5]
    sage: EllipticCurve('389a1').prove_BSD(proof=True) == Primes()
    True
    sage: B = EllipticCurve([-25,0]).prove_BSD(proof=True, return_BSD=True)
    sage: B.rank, B.primes, 'Kob2013' in B.proof['outside']['reference']
    (1, [5], True)

A curve beyond the mini database (conductor 9999) needs no database label::

    sage: E = EllipticCurve([0,0,1,-3,5])
    sage: E.conductor()
    10179
    sage: B = E.prove_BSD(proof=True, return_BSD=True)  # long time
    sage: B.primes, B.certification['analytic_sha']['certified']  # long time
    ([3, 29], True)

An unavailable calculation gives an inconclusive result, and never []::

    sage: with patch.object(cert._AnalyticData, 'sha', side_effect=NotImplementedError('test')):
    ....:     B = EllipticCurve('11a1').prove_BSD(proof=True, return_BSD=True)
    sage: B.primes == Primes() and 'unavailable' in B.unresolved['all']
    True
    sage: with patch.object(bsd, 'mwrank_two_descent_work', side_effect=NotImplementedError('test')):
    ....:     B = EllipticCurve('11a1').prove_BSD(proof=True, return_BSD=True)
    sage: B.primes, B.descent['mwrank']['inconclusive'], 'pari' in B.descent
    ([], 'test', True)

Verbosity zero stays silent, including on the conditional path::

    sage: import contextlib, io
    sage: output = io.StringIO()
    sage: with contextlib.redirect_stdout(output):
    ....:     _ = EllipticCurve('11a1').prove_BSD(proof=True, rigorous=False)
    sage: output.getvalue()
    ''

Unpublished results are not silently enabled; the witnesses for discussion
remain unresolved at 3::

    sage: EllipticCurve('91b1').prove_BSD(proof=True)
    [3]
    sage: EllipticCurve('2534f1').prove_BSD(proof=True)
    [3]
"""

# Distributed under the terms of the GNU General Public License (GPL).
