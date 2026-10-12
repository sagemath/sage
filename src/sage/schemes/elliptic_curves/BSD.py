r"""
Birch and Swinnerton-Dyer formulas

The proof engine first establishes a finite set of exceptional primes, then
checks explicit published criteria at each of them.  Its analytic inputs are
certified by default; see :func:`prove_BSD` for the conditional mode.

This replaces the former use of disputed cohomological arguments and
conductor-based tables by individually recorded theorem applications.

AUTHORS:

- Robert L. Miller, William Stein, and Christian Wuthrich: earlier implementation
- Christian Wuthrich (2026): theorem-based rewrite
- Barinder Banwait (2026): integration, certification, and Skinner's criterion

The rewrite is based on Christian Wuthrich's ``new_bsd_prove.py``, revision
``fae09f41a2596593ed0e0ac4720b5a86b52cb9e0``.
"""

# Distributed under the terms of the GNU General Public License (GPL).

from cypari2.handle_error import PariError

from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.sets.primes import Primes


class BSD_data:
    """
    Information and evidence collected while attempting to prove BSD.

    ``proof[p]`` records a successful theorem application.  ``proof['outside']``
    records the theorem establishing the finite exceptional set.  ``unresolved``
    records reasons for inconclusive results.  ``certification`` and
    ``assumptions`` distinguish proved analytic inputs from assumed ones.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import BSD_data
        sage: D = BSD_data()
        sage: D.Sha is None
        True
        sage: D.curve = EllipticCurve('11a1')
        sage: D.update()
        sage: D.N, D.sha_an is None
        (11, True)

    Updating the curve does not silently evaluate a numerical analytic order.
    The proof engine supplies that value together with its provenance.
    """
    def __init__(self):
        self.curve = None
        self.original_curve = None
        self.two_tor_rk = None
        self.Sha = None
        self.sha_an = None
        self.N = None
        self.rank = None
        self.gens = None
        self.bounds = {}
        self.primes = None
        self.heegner_indexes = {}
        self.heegner_index_upper_bound = {}
        self.N_factorization = None
        self.proof = {}
        self.rigorous = True
        self.assumptions = []
        self.certification = {'analytic_rank': {'certified': False},
                              'analytic_sha': {'certified': False}}
        self.unresolved = {}
        self.descent = {}

    def update(self):
        """
        Update exact arithmetic properties from ``curve``.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.BSD import BSD_data
            sage: D = BSD_data()
            sage: D.curve = EllipticCurve('14a1')
            sage: D.update()
            sage: D.N, D.two_tor_rk
            (14, 1)
        """
        self.two_tor_rk = self.curve.two_torsion_rank()
        self.Sha = self.curve.sha()
        self.N = self.curve.conductor()
        self.N_factorization = self.N.factor()


def mwrank_two_descent_work(E, two_tor_rk) -> tuple:
    """
    Prepare the output from mwrank two-descent.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``two_tor_rk`` -- its two-torsion rank

    OUTPUT:

    - a lower bound on the rank

    - an upper bound on the rank

    - a lower bound on the rank of Sha[2]

    - an upper bound on the rank of Sha[2]

    - a list of the generators found

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import mwrank_two_descent_work
        sage: E = EllipticCurve('14a')
        sage: mwrank_two_descent_work(E, E.two_torsion_rank())
        (0, 0, 0, 0, [])
        sage: E = EllipticCurve('37a')
        sage: mwrank_two_descent_work(E, E.two_torsion_rank())
        (1, 1, 0, 0, [(0 : -1 : 1)])
    """
    MWRC = E.mwrank_curve()
    rank_upper_bd = MWRC.rank_bound()
    rank_lower_bd = MWRC.rank()
    gens = [E(P) for P in MWRC.gens()]
    sha2_lower_bd = MWRC.selmer_rank() - two_tor_rk - rank_upper_bd
    sha2_upper_bd = MWRC.selmer_rank() - two_tor_rk - rank_lower_bd
    return rank_lower_bd, rank_upper_bd, sha2_lower_bd, sha2_upper_bd, gens


def pari_two_descent_work(E) -> tuple:
    r"""
    Prepare the output from pari by two-isogeny.

    INPUT:

    - ``E`` -- an elliptic curve

    OUTPUT: a tuple of 5 elements with the first 4 being integers

    - a lower bound on the rank

    - an upper bound on the rank

    - a lower bound on the rank of Sha[2]

    - an upper bound on the rank of Sha[2]

    - a list of the generators found

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import pari_two_descent_work
        sage: E = EllipticCurve('14a')
        sage: pari_two_descent_work(E)
        (0, 0, 0, 0, [])
        sage: E = EllipticCurve('37a')
        sage: pari_two_descent_work(E) # random, up to sign
        (1, 1, 0, 0, [(0 : -1 : 1)])
        sage: E = EllipticCurve('210e7')
        sage: pari_two_descent_work(E)
        (0, 2, 0, 2, [])
        sage: E = EllipticCurve('66b3')
        sage: pari_two_descent_work(E)
        (0, 0, 2, 2, [])
    """
    ep = E.pari_curve()
    lower, rank_upper_bd, s, pts = ep.ellrank()
    gens = sorted([E.point([QQ(x[0]), QQ(x[1])], check=True) for x in pts])
    gens = E.saturation(gens)[0]
    # this is explained in the pari-gp documentation:
    # s is the dimension of Sha[2]/2Sha[4],
    # which is a lower bound for dim Sha[2]
    # dim Sha[2] = dim Sel2 - rank E(Q) - dim tors
    # rank_upper_bd = dim Sel_2 - dim tors - s
    sha_upper_bd = rank_upper_bd - lower + s
    return ZZ(lower), ZZ(rank_upper_bd), ZZ(s), ZZ(sha_upper_bd), gens


def native_two_isogeny_descent_work(E, two_tor_rk) -> tuple:
    """
    Prepare the output from two-descent by two-isogeny.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``two_tor_rk`` -- its two-torsion rank

    OUTPUT:

    - a lower bound on the rank

    - an upper bound on the rank

    - a lower bound on the rank of Sha[2]

    - an upper bound on the rank of Sha[2]

    - a list of the generators found
      (currently ``None``, since we do not store them)

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import native_two_isogeny_descent_work
        sage: E = EllipticCurve('14a')
        sage: native_two_isogeny_descent_work(E, E.two_torsion_rank())
        (0, 0, 0, 0, None)
        sage: E = EllipticCurve('65a')
        sage: native_two_isogeny_descent_work(E, E.two_torsion_rank())
        (1, 1, 0, 0, None)
    """
    from sage.schemes.elliptic_curves.descent_two_isogeny import (
        two_descent_by_two_isogeny,
    )
    result_two_descent = [ZZ(n) for n in two_descent_by_two_isogeny(E)]
    # safety check that all numbers in the result are powers of two
    if not all(n.is_power_of(2) for n in result_two_descent):
        raise RuntimeError("not a power of 2 in two-descent")

    e1, e2, e1p, e2p = (n.valuation(2) for n in result_two_descent)
    rank_lower_bd = e1 + e1p - 2
    rank_upper_bd = e2 + e2p - 2
    sha_upper_bd = e2 + e2p - e1 - e1p
    gens = None  # right now, we are not keeping track of them
    return rank_lower_bd, rank_upper_bd, 0, sha_upper_bd, gens


def _evidence(reference, **witnesses):
    """
    Construct a theorem record.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _evidence
        sage: _evidence('test', prime=3)
        {'reference': 'test', 'witnesses': {'prime': 3}}
    """
    return {'reference': reference, 'witnesses': witnesses}


def _ramified_multiplicative_primes(E, p, nonsplit=False):
    r"""
    Find primes `q\ne p` of multiplicative reduction with ramified `E[p]`.

    The Tate-curve criterion is `p\nmid v_q(\Delta_{\min})`.  It applies
    equally to nonsplit multiplicative reduction, after an unramified twist.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _ramified_multiplicative_primes
        sage: _ramified_multiplicative_primes(EllipticCurve('5389a1'), 3)
        [17, 317]
        sage: 3 in _ramified_multiplicative_primes(EllipticCurve('681b1'), 3)
        False
    """
    E = E.minimal_model()
    return [q for q in E.conductor().prime_divisors()
            if q != p and E.has_multiplicative_reduction(q)
            and (not nonsplit or E.has_nonsplit_multiplicative_reduction(q))
            and E.discriminant().valuation(q) % p != 0]


def _kernel_has_local_point(phi, p):
    r"""
    Test whether the kernel of an odd-degree isogeny has a nonzero `\QQ_p`-point.

    For each irreducible kernel-polynomial factor, construct the field of the
    x-coordinate and then of the y-coordinate.  A `\QQ_p`-embedding exists
    exactly when there is a prime above `p` with ramification and residue
    degrees both one.  Thus this test uses exact number-field arithmetic,
    not fixed-precision root approximations.  Degrees are at most `p-1`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _kernel_has_local_point
        sage: E = EllipticCurve('11a1')
        sage: phis = E.isogenies_prime_degree(5)
        sage: any(_kernel_has_local_point(phi, 5) for phi in phis)
        True
    """
    from sage.rings.number_field.number_field import NumberField
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    E = phi.domain()
    b2, b4, b6, _ = E.b_invariants()
    for f, _ in phi.kernel_polynomial().factor():
        if f.degree() == 1:
            K, x = QQ, -f[0]/f[1]
        else:
            K = NumberField(f, 'a')
            x = K.gen()
            if not any(q.ramification_index() == q.residue_class_degree() == 1
                       for q in K.primes_above(p)):
                continue
        y2 = 4*x**3 + b2*x**2 + 2*b4*x + b6
        if y2.is_square():
            return True
        t = PolynomialRing(K, 't').gen()
        if K is QQ:
            L = NumberField(t**2 - y2, 'b')
        else:
            L = K.extension(t**2 - y2, 'b').absolute_field('z')
        if any(q.ramification_index() == q.residue_class_degree() == 1
               for q in L.primes_above(p)):
            return True
    return False


def _tate_module_surjectivity(E, p):
    r"""
    Certify surjectivity on the full Tate module, or return ``None``.

    For `p\geq5`, mod-`p` surjectivity suffices ([SW2013]_, Proposition 7.2).
    At 3 use a multiplicative prime `q\ne3` with `3\nmid v_q(\Delta)`.
    Tate uniformization gives the entire unipotent subgroup `U(\ZZ_3)`
    in the image of inertia at `q`.  A conjugate moving its fixed line
    modulo 3 gives a second primitive line.  In the basis of these two
    lines the two groups are the upper and lower elementary matrices;
    they generate `\mathrm{SL}_2(\ZZ_3)`.  The cyclotomic determinant is
    surjective.  This verifies [Kat2004]_, condition (12.5.2), without
    assuming a mod-3 lifting assertion.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _tate_module_surjectivity
        sage: _tate_module_surjectivity(EllipticCurve('37a1'), 3)
        {'multiplicative_inertia_prime': 37, 'surjective_mod_p': True}
        sage: _tate_module_surjectivity(EllipticCurve('11a1'), 5) is None
        True
    """
    if p < 3 or not E.galois_representation().is_surjective(p):
        return None
    if p >= 5:
        return {'surjective_mod_p': True, 'lifting': 'SW2013, Proposition 7.2'}
    qs = _ramified_multiplicative_primes(E, p)
    if qs:
        return {'surjective_mod_p': True, 'multiplicative_inertia_prime': qs[0]}


def _multiplicative_local_torsion(E, p):
    r"""
    Test for nonzero `p`-torsion over `\QQ_p` at an odd multiplicative prime.

    In the split case, write the Tate parameter as `q=p^v u`.  It is a
    `p`-th power precisely when `p\mid v` and `u^{p-1}=1\pmod {p^2}`.
    In the nonsplit case the reduction and component groups have order
    coprime to `p`, and the formal group has no `p`-torsion.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _multiplicative_local_torsion
        sage: _multiplicative_local_torsion(EllipticCurve('11a1'), 11)
        False
        sage: _multiplicative_local_torsion(EllipticCurve('37a1'), 37)
        False
        sage: _multiplicative_local_torsion(EllipticCurve('11a1'), 5)
        Traceback (most recent call last):
        ...
        ValueError: an odd prime of multiplicative reduction is required
    """
    if p <= 2 or not E.has_multiplicative_reduction(p):
        raise ValueError("an odd prime of multiplicative reduction is required")
    if E.has_nonsplit_multiplicative_reduction(p):
        return False
    v = E.minimal_model().discriminant().valuation(p)
    if v % p:
        return False
    q = E.tate_curve(p).parameter(prec=v + 8)
    u = q / p**v
    if u.precision_absolute() < 2:
        raise ArithmeticError("insufficient precision in the Tate parameter")
    return (u**(p-1) - 1).valuation() >= 2


def _skinner(E, p, analytic):
    r"""
    Apply [Skinner2016]_, Theorem C (analytic rank zero).

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _skinner
        sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
        sage: E = EllipticCurve('5389a1')
        sage: _skinner(E, 3, _AnalyticData(E, False))['witnesses']['q']
        17
        sage: E = EllipticCurve('2534f1')
        sage: _skinner(E, 3, _AnalyticData(E, False)) is None
        True
    """
    if p <= 2 or analytic.rank != 0:
        return None
    if not (E.has_multiplicative_reduction(p)
            or (E.has_good_reduction(p) and E.ap(p) % p != 0)):
        return None
    if not E.galois_representation().is_irreducible(p):
        return None
    qs = _ramified_multiplicative_primes(E, p)
    if qs:
        return _evidence('Skinner2016, Theorem C', q=qs[0])


def _bcs(E, p, analytic):
    """
    Apply [BCS2025]_, Corollary 1.3.1, using surjectivity to imply (im).

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _bcs
        sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
        sage: E = EllipticCurve('37a1')
        sage: _bcs(E, 5, _AnalyticData(E, False))['reference']
        'BCS2025, Corollary 1.3.1'
        sage: _bcs(E, 3, _AnalyticData(E, False)) is None
        True
    """
    if (p > 3 and analytic.rank in (0, 1) and not E.has_cm()
            and E.has_good_reduction(p) and E.ap(p) % p != 0
            and E.galois_representation().is_surjective(p)):
        return _evidence('BCS2025, Corollary 1.3.1', surjective=True)


def _cgs(E, p, analytic):
    r"""
    Apply [CGS2025]_, Theorem D, with its local character hypotheses.

    The kernel characters of an isogeny and its dual are `\phi` and
    `\omega\phi^{-1}`.  Neither may be trivial over `\QQ_p`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _cgs
        sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
        sage: E = EllipticCurve('91b1')
        sage: _cgs(E, 3, _AnalyticData(E, False)) is None
        True
        sage: E = EllipticCurve('50b1')
        sage: _cgs(E, 5, _AnalyticData(E, False)) is None
        True
    """
    if p <= 2 or analytic.rank not in (0, 1) or not E.has_good_reduction(p):
        return None
    if not E.galois_representation().is_reducible(p):
        return None
    for phi in E.isogenies_prime_degree(p):
        if not _kernel_has_local_point(phi, p) and not _kernel_has_local_point(phi.dual(), p):
            return _evidence('CGS2025, Theorem D',
                             isogeny_codomain=phi.codomain().a_invariants(),
                             local_characters_nontrivial=True)


def _castella(E, p, analytic):
    """
    Apply [Castella2018err]_, Theorem A', including the restriction p > 3.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _castella
        sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
        sage: E = EllipticCurve('123a1')
        sage: _castella(E, 3, _AnalyticData(E, False)) is None
        True
    """
    if p <= 3 or analytic.rank != 1 or not E.has_multiplicative_reduction(p):
        return None
    if not E.galois_representation().is_irreducible(p):
        return None
    qs = _ramified_multiplicative_primes(E, p, nonsplit=True)
    if qs and not _multiplicative_local_torsion(E, p):
        return _evidence("Castella2018err, Theorem A'", q=qs[0], local_torsion=False)


def _stein_wuthrich(E, p, analytic):
    """
    Apply Theorems 8.1 and 9.1 of [SW2013]_ when the analytic valuation is zero.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _stein_wuthrich
        sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
        sage: E = EllipticCurve('11a1')
        sage: _stein_wuthrich(E, 11, _AnalyticData(E, True))['reference']
        'SW2013, Theorem 8.1'
    """
    if p <= 2 or analytic.sha(E).valuation(p) != 0:
        return None
    tate_module = _tate_module_surjectivity(E, p)
    if tate_module is None:
        return None
    if analytic.rank == 0 and not E.has_additive_reduction(p):
        return _evidence('SW2013, Theorem 8.1', tate_module=tate_module)
    if analytic.rank == 1 and E.has_good_reduction(p) and E.ap(p) % p != 0:
        P = analytic.gens(E)[0]
        # A zero to finite p-adic precision is inconclusive, never a proof
        # that the height vanishes.
        h = E.padic_height(p)(P)
        if h != 0:
            return _evidence('SW2013, Theorem 9.1', padic_height=h, tate_module=tate_module)


def _iwasawa(E, p, analytic):
    """
    Use the upper-bound algorithm of [SW2013]_ and [Wuthrich2014]_.

    Both the upper bound and the analytic valuation must be zero.  The
    reducible case uses the corrected integrality theorem of [Wuthrich2014]_,
    not the earlier argument with a gap.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _iwasawa
        sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
        sage: E = EllipticCurve('11a3')
        sage: _iwasawa(E, 5, _AnalyticData(E, True))['witnesses']['bound']
        0
        sage: E = EllipticCurve('91b1')
        sage: _iwasawa(E, 3, _AnalyticData(E, True)) is None
        True
    """
    if p <= 2 or analytic.rank not in (0, 1) or E.has_additive_reduction(p):
        return None
    if analytic.sha(E).valuation(p) != 0:
        return None
    rho = E.galois_representation()
    reducible = rho.is_reducible(p)
    tate_module = None if reducible else _tate_module_surjectivity(E, p)
    if not reducible and tate_module is None:
        return None
    if E.has_nonsplit_multiplicative_reduction(p) and (analytic.rank or p == 3):
        return None
    if p == 3 and (analytic.rank or (E.has_good_reduction(p) and E.ap(p) % p == 0)):
        return None
    from sage.schemes.elliptic_curves.bsd_certificates import _padic_sha_bound
    bound = _padic_sha_bound(E, p, analytic)
    if bound == 0:
        return _evidence('SW2013, upper-bound algorithm; Wuthrich2014, Theorem 3',
                         bound=bound, reducible=reducible, tate_module=tate_module)


def _two_primary(E, analytic, backend, result):
    """
    Compare descent information with the analytic 2-adic valuation.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import BSD_data, _two_primary
        sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
        sage: E = EllipticCurve('66b3')
        sage: d = BSD_data()
        sage: _two_primary(E, _AnalyticData(E, True), 'pari', d)
        True
        sage: d.descent['pari']['sha2_dimension']
        2
    """
    r, v = analytic.rank, analytic.sha(E).valuation(2)
    if backend == 'sage' and E.two_torsion_rank() == 0:
        result.descent['sage'] = {'inconclusive': 'two-isogeny descent requires rational 2-torsion'}
    elif backend != 'pari':
        if backend == 'mwrank':
            lo, hi, slo, shi, pts = mwrank_two_descent_work(E, E.two_torsion_rank())
        else:
            lo, hi, slo, shi, pts = native_two_isogeny_descent_work(E, E.two_torsion_rank())
        result.descent[backend] = {'rank_bounds': (lo, hi),
                                   'sha2_bounds': (slo, shi), 'points': pts}
        if not lo <= r <= hi:
            raise ArithmeticError('descent contradicts the analytic rank')
        if shi == 0:
            if v != 0:
                raise ArithmeticError('trivial Sha[2] contradicts the analytic 2-part')
            result.bounds[2] = (0, 0)
            result.proof[2] = _evidence('2-descent', backend=backend)
            return True
    lo, hi, s, pts = E.pari_curve().ellrank()
    lo, hi, s = ZZ(lo), ZZ(hi), ZZ(s)
    d = hi + s - r
    result.descent['pari'] = {'rank_bounds': (lo, hi), 'cassels_tate_dimension': s,
                              'sha2_dimension': d, 'points': pts.sage()}
    if not (lo <= r <= hi and 0 <= s <= d):
        raise ArithmeticError('PARI descent contradicts the analytic rank')
    result.bounds[2] = (2*d-s, infinity if d > s else d)
    if d == s:
        if v != d:
            raise ArithmeticError('Cassels-Tate pairing contradicts the analytic 2-part')
        result.proof[2] = _evidence('PARI 2-descent and Cassels-Tate pairing',
                                    sha2_dimension=d, cassels_tate_dimension=s)
        return True
    result.unresolved[2] = '2-descent does not determine the full 2-primary order'
    return False


def _jetchev_exceptions(E, analytic):
    r"""
    Establish a finite exceptional set using a specified Heegner point.

    Let `J=4\hat h(y_D)/R` be the index square from Gross--Zagier.
    For odd `p`, `E(\QQ)` has index prime to `p` in the free part of
    `E(K)`.  If `E[p]` is surjective, neither group has `p`-torsion.
    Hence `v_p(J)=2m_0` in [Jet2008]_, Corollary 1.5.  At good primes
    that corollary bounds the order of Sha over `K` by
    `p^{v_p(J)-2\max_q v_p(c_q)}`.  Restriction injects the odd-primary
    Sha over `\QQ` into Sha over `K`.

    In particular, away from the support of `J/\operatorname{lcm}(c_q)^2`
    this proves triviality.  Include bad primes, nonsurjective primes,
    and primes dividing the analytic order separately.  We conservatively
    exclude primes dividing the Manin constant as well.  The hypotheses
    of that corollary do not require `p\nmid D`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import _jetchev_exceptions
        sage: from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
        sage: E = EllipticCurve('37a1')
        sage: primes, evidence = _jetchev_exceptions(E, _AnalyticData(E, True))
        sage: sorted(primes)
        [2, 37]
        sage: evidence['index_square']
        4
    """
    from sage.arith.functions import lcm

    h = analytic.heegner_data(E)
    J = analytic.sha(E) / h['lattice']
    if J not in ZZ or J <= 0 or not ZZ(J).is_square():
        raise ArithmeticError('analytic inputs contradict the Heegner index square')
    J = ZZ(J)
    quotient = J / lcm(E.tamagawa_numbers())**2
    primes = set(E.galois_representation().non_surjective())
    for n in (E.conductor(), h['manin_constant'],
              quotient.numerator(), quotient.denominator(),
              analytic.sha(E).numerator(), analytic.sha(E).denominator()):
        primes.update(n.prime_divisors())
    return primes, dict(h, index_square=J, tamagawa_corrected_bound=quotient,
                        certified=analytic.rigorous)


def prove_BSD(E, verbosity=0, two_desc='mwrank', proof=None, secs_hi=5,
              return_BSD=False, *, rigorous=True):
    r"""
    Attempt to prove the Birch and Swinnerton-Dyer formula for ``E``.

    Return a sorted list of primes at which this algorithm has not proved
    the formula.  An empty list proves the formula when ``proof`` and
    ``rigorous`` are both true.  This is not a list of all primes at which
    the formula is unknown in the mathematical literature.

    INPUT:

    - ``E`` -- an elliptic curve over the rational numbers
    - ``verbosity`` -- integer (default: 0); 0 prints nothing, 1 prints
      theorem applications, and 2 also describes unresolved primes
    - ``two_desc`` -- ``'mwrank'`` (default), ``'pari'``, or ``'sage'``;
      initial descent backend; PARI supplements an inconclusive result
    - ``proof`` -- boolean or ``None`` (default); ``None`` uses the global
      elliptic-curve proof flag; ``False`` immediately returns ``[]``,
      assuming BSD rather than attempting a proof
    - ``secs_hi`` -- retained for compatibility; no longer used, and does
      not limit the analytic certification computations
    - ``return_BSD`` -- boolean (default: ``False``); return :class:`BSD_data`
      containing the primes, certificates, assumptions, and proof evidence
    - ``rigorous`` -- boolean (default: ``True``); certify analytic inputs.
      If ``False``, assume the values from ``analytic_rank()`` and
      ``sha().an()`` and make the result conditional on these values

    OUTPUT:

    A list of primes, or :class:`~sage.sets.primes.Primes` when no finite
    exceptional set is established; alternatively, a :class:`BSD_data`.

    ALGORITHM:

    Certify rank zero or one and the rational analytic order of Sha using
    :mod:`sage.schemes.elliptic_curves.bsd_certificates`.  Use isogeny
    invariance, 2-descent, and the Cassels-Tate pairing.  For non-CM curves,
    [Kat2004]_ (Theorem 14.5) or [Jet2008]_ (Corollary 1.5) establishes a
    finite exceptional set.  Test published criteria of [SW2013]_,
    [Skinner2016]_, [Castella2018err]_, [CGS2025]_, and [BCS2025]_, followed
    by the Iwasawa upper-bound algorithm with [Wuthrich2014]_.  For CM
    curves use [Rub1991]_, [PR1987]_, and [Kob2013]_.

    Exact certification, especially the Manin constant and saturation, can
    be expensive.  ``rigorous=False`` bypasses the analytic certification;
    it does not bypass descent or theorem hypotheses.  Unpublished BSTW,
    Keller--Yin, and CCSS criteria are not used.

    EXAMPLES::

        sage: EllipticCurve('11a1').prove_BSD(proof=True)
        []
        sage: EllipticCurve('91b1').prove_BSD(proof=True)
        [3]
        sage: EllipticCurve('389a1').prove_BSD(proof=True)
        Set of all prime numbers: 2, 3, 5, 7, ...
        sage: EllipticCurve('5389a1').prove_BSD(proof=True)  # long time
        []

    Inspect the evidence and the conditional mode::

        sage: E = EllipticCurve('11a1')
        sage: B = E.prove_BSD(proof=True, return_BSD=True)
        sage: B.primes, B.rigorous, B.assumptions
        ([], True, [])
        sage: B.certification['analytic_sha']['certified']
        True
        sage: B = E.prove_BSD(proof=True, rigorous=False, return_BSD=True)
        sage: B.primes, B.rigorous
        ([], False)
        sage: B.assumptions
        ['analytic_rank() and sha().an() are assumed correct']
        sage: E.prove_BSD(proof=False)
        []

    TESTS:

    No assumption is reused by a subsequent rigorous call::

        sage: E.prove_BSD(proof=True, rigorous=True, return_BSD=True).assumptions
        []
        sage: for backend in ['mwrank', 'pari', 'sage']:
        ....:     assert EllipticCurve('14a1').prove_BSD(proof=True, two_desc=backend) == []
        sage: E.prove_BSD(proof=True, two_desc='unknown')
        Traceback (most recent call last):
        ...
        ValueError: two_desc must be 'mwrank', 'pari', or 'sage'
        sage: F = E.change_weierstrass_model([1/2, 0, 0, 0])
        sage: F.prove_BSD(proof=True) == E.prove_BSD(proof=True)
        True
    """
    from sage.schemes.elliptic_curves.bsd_certificates import _AnalyticData
    from sage.structure.proof.proof import get_flag

    if not get_flag(proof, 'elliptic_curve'):
        return []
    if two_desc not in ('mwrank', 'pari', 'sage'):
        raise ValueError("two_desc must be 'mwrank', 'pari', or 'sage'")
    result = BSD_data()
    result.original_curve = E
    result.rigorous = bool(rigorous)
    if not rigorous:
        result.assumptions = ['analytic_rank() and sha().an() are assumed correct']
        if verbosity:
            print('Conditional on analytic_rank() and sha().an() being correct.')
    E = E.minimal_model()
    try:
        analytic = _AnalyticData(E, rigorous)
    except (NotImplementedError, PariError, MemoryError) as error:
        result.curve = E
        result.update()
        result.primes = Primes()
        result.unresolved['all'] = 'analytic rank certification unavailable: %s' % error
        if verbosity:
            print(result.unresolved['all'])
        return result if return_BSD else result.primes
    result.certification['analytic_rank'] = analytic.rank_certificate
    result.rank = analytic.rank
    if analytic.rank not in (0, 1):
        result.curve = E
        result.update()
        result.primes = Primes()
        result.unresolved['all'] = 'analytic rank zero or one has not been certified' if rigorous else 'assumed analytic rank exceeds one'
        if verbosity:
            print(result.unresolved['all'])
        return result if return_BSD else result.primes

    curves = [C.minimal_model() for C in E.isogeny_class().curves]
    try:
        selected = min(curves, key=lambda C: (analytic.sha(C), C.a_invariants()))
    except (NotImplementedError, PariError, MemoryError) as error:
        result.curve = E
        result.update()
        result.primes = Primes()
        result.certification['isogeny_class'] = analytic.certificates
        result.unresolved['all'] = 'analytic Sha certification unavailable: %s' % error
        if verbosity:
            print(result.unresolved['all'])
        return result if return_BSD else result.primes
    result.curve = selected
    result.update()
    result.sha_an = analytic.sha(selected)
    result.certification['analytic_sha'] = analytic.certificates[selected.a_invariants()]
    result.certification['isogeny_class'] = analytic.certificates
    result.gens = analytic.gens(selected)
    if verbosity and selected != E:
        print('Using the isogenous minimal curve with coefficients %s.' % (selected.a_invariants(),))

    remaining = set()
    try:
        two_proved = _two_primary(selected, analytic, two_desc, result)
    except (NotImplementedError, PariError, MemoryError) as error:
        result.descent[two_desc] = {'inconclusive': str(error)}
        two_proved = False
        if two_desc != 'pari':
            try:
                two_proved = _two_primary(selected, analytic, 'pari', result)
            except (NotImplementedError, PariError, MemoryError) as pari_error:
                result.descent['pari'] = {'inconclusive': str(pari_error)}
        if not two_proved:
            result.unresolved[2] = 'descent computation unavailable; see descent diagnostics'
    if not two_proved:
        remaining.add(ZZ(2))
    elif verbosity:
        print('BSD at 2 by %s.' % result.proof[2]['reference'])
    witnesses = {}
    if selected.has_cm():
        maximal_order = min((C for C in curves
                             if C.cm_discriminant().is_fundamental_discriminant()),
                            key=lambda C: C.a_invariants())
        witnesses['maximal_order_curve'] = maximal_order.a_invariants()
        if analytic.rank == 0:
            exceptions = {ZZ(3)} if maximal_order.j_invariant() == 0 else set()
            reference = 'Rub1991, main theorem'
        else:
            exceptions = set(selected.conductor().prime_divisors()) - {ZZ(2)}
            reference = 'PR1987; Kob2013, CM rank-one BSD corollary'
    else:
        exceptions = set(selected.galois_representation().non_surjective())
        exceptions.update(result.sha_an.numerator().prime_divisors())
        exceptions.update(result.sha_an.denominator().prime_divisors())
        if analytic.rank == 0:
            exceptions.update(selected.j_invariant().denominator().prime_divisors())
            tate_module = _tate_module_surjectivity(selected, ZZ(3))
            if tate_module is None:
                exceptions.add(ZZ(3))
            else:
                witnesses['tate_module_at_3'] = tate_module
            reference = 'Kat2004, Theorem 14.5'
        else:
            try:
                exceptions, witnesses = _jetchev_exceptions(selected, analytic)
            except (NotImplementedError, PariError, MemoryError) as error:
                result.primes = Primes()
                result.unresolved['all'] = 'Heegner bound unavailable: %s' % error
                if verbosity:
                    print(result.unresolved['all'])
                return result if return_BSD else result.primes
            reference = 'Jet2008, Corollary 1.5'
        exceptions.discard(ZZ(2))
    result.proof['outside'] = _evidence(reference, **witnesses,
                                        exceptional_primes=sorted(exceptions | {ZZ(2)}))
    if verbosity:
        print('BSD outside %s by %s.' % (sorted(exceptions | {ZZ(2)}), reference))
    criteria = (_skinner, _bcs, _cgs, _castella, _stein_wuthrich, _iwasawa)
    for p in sorted(exceptions):
        failures = []
        for criterion in criteria:
            try:
                evidence = criterion(selected, p, analytic)
            except (NotImplementedError, PariError, MemoryError) as error:
                failures.append('%s: %s' % (criterion.__name__, error))
                continue
            if evidence is not None:
                evidence['certified_analytic_inputs'] = analytic.rigorous
                evidence['witnesses'].update(analytic_rank=analytic.rank,
                                             analytic_valuation=result.sha_an.valuation(p))
                result.proof[p] = evidence
                result.bounds[p] = (result.sha_an.valuation(p), result.sha_an.valuation(p))
                if verbosity:
                    print('BSD at %s by %s.' % (p, evidence['reference']))
                break
        else:
            remaining.add(p)
            if selected.has_good_reduction(p):
                reduction = 'good ordinary' if selected.ap(p) % p else 'good supersingular'
            elif selected.has_multiplicative_reduction(p):
                reduction = 'split multiplicative' if selected.ap(p) == 1 else 'nonsplit multiplicative'
            else:
                reduction = 'additive'
            rho = selected.galois_representation()
            if rho.is_reducible(p):
                representation = 'reducible'
                if selected.has_good_reduction(p):
                    failures.append('CGS local-character hypothesis fails')
            else:
                representation = 'surjective mod p' if rho.is_surjective(p) else 'surjectivity not established'
            failures.append('no implemented criterion verified (rank %s, %s, %s, analytic valuation %s)'
                            % (analytic.rank, reduction, representation, result.sha_an.valuation(p)))
            result.unresolved[p] = '; '.join(failures)
    result.primes = sorted(remaining)
    if verbosity > 1:
        for p in result.primes:
            print('Unresolved at %s: %s.' % (p, result.unresolved[p]))
    return result if return_BSD else result.primes
