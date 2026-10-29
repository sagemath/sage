r"""
Certified bounds and enumeration for rational integral points

The elliptic logarithm method [ST1994b]_ and David's final lower bound
[Dav1995]_, Theorem 2.1, give a finite coefficient bound. Arb balls
enclose the analytic quantities, exact Gram--Schmidt lengths certify
lattice reduction, and interval enumeration has rational fallbacks.

Write `H(P)=\lim_{n\to\infty}h(x(2^nP))/4^n`. Silverman's comparison
[Sil1990]_, Theorem 1.1, bounds `|H(P)-h(x(P))|`. Exact doubling and
polarization enclose the height matrix; interval LDL decomposition
certifies `H(\sum m_iP_i+T)\geq\lambda M^2`, where
`M=\max|m_i|` and `\lambda>0`. For the short coordinate `s=x+b_2/12`,
the cubic root bound separates a bounded region from one where the
real elliptic integral gives `|2u(P)|\leq K\exp(-\lambda M^2/2)`.

The logarithms of `2P_i` lie on the identity component, and
`2u(P)=\sum m_i u(2P_i)+b\omega`, where the denominator of `b` divides
the torsion exponent `e`. Apply David's final theorem to `e\,2u(P)`,
including the period as a logarithm of the origin. A common parameter
bounds curve heights, analytic terms and the projective canonical
heights `6H(P_i)` of `2P_i`. The enlarged integer constant and a
polynomial majorant give a finite initial bound; its logarithmic
derivative excludes every larger coefficient. See [ST1999]_ for the
distinction between the preliminary and final bounds.

For reduction, round normalized logarithms at an integer scale,
enclosing the rounding errors. Every nonzero lattice vector has
squared length at least the minimum exact squared Gram--Schmidt
length. Use the resulting lower bound for the linear form only after
certifying positivity, increasing scale and precision as needed.
In lattice dimensions at most eight, [FP1985]_ enumeration computes
the exact shortest-vector length. The final search intersects the
coefficient box with a certified height ellipsoid and includes all
torsion translates and both signs.

AUTHORS:

- Manvir Jaswal (2026): certified bounds and enumeration

EXAMPLES::

    sage: from sage.schemes.elliptic_curves.integral_points import integral_points
    sage: [P[0] for P in integral_points(EllipticCurve([1,-1,1,-11,16]))]
    [-1, 2, 20829482]
    sage: integral_points(EllipticCurve([0,0,0,0,-365]))
    []
"""

# *****************************************************************************
#       Copyright (C) 2026 Manvir Jaswal
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.
# *****************************************************************************

import warnings
from contextlib import contextmanager
from math import isqrt

from sage.libs.pari import pari
from sage.matrix.constructor import matrix
from sage.misc.decorators import sage_wraps
from sage.rings.complex_arb import ComplexBallField
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.qqbar import QQbar
from sage.rings.rational_field import QQ
from sage.rings.real_arb import RealBallField
from sage.rings.real_double import RDF
from sage.structure.proof.proof import WithProof, get_flag


@contextmanager
def _proved_arithmetic():
    """Require proved factorizations and restore the caller's settings."""
    previous = pari.default('factor_proven')
    try:
        pari.default('factor_proven', 1)
        with WithProof('arithmetic', True), WithProof('elliptic_curve', True):
            yield
    finally:
        pari.default('factor_proven', previous)


def _with_proved_arithmetic(func):
    """Run ``func`` with proved arithmetic, keeping its signature."""
    @sage_wraps(func)
    def wrapper(*args, **kwds):
        with _proved_arithmetic():
            return func(*args, **kwds)
    return wrapper


def _upper(x):
    """Return a rational upper endpoint of a finite real ball."""
    if not x.is_finite():
        raise ArithmeticError("a finite enclosure could not be computed")
    return x.upper().exact_rational()


def _ceil_upper(x):
    """Round a real ball upwards, including its error."""
    return _upper(x).ceil()


@contextmanager
def _pari_seed(seed):
    """
    Fix PARI's random seed and restore the caller's state afterwards.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _pari_seed
        sage: pari.setrand(5); state = pari.getrand()
        sage: with _pari_seed(1):
        ....:     a = pari.random()
        sage: pari.getrand() == state
        True
        sage: with _pari_seed(1):
        ....:     a == pari.random()
        True
    """
    state = pari.getrand()
    pari.setrand(seed)
    try:
        yield
    finally:
        pari.setrand(state)


def _pari_descent(E, certify=True):
    """Initialize PARI's descent, optionally certifying its cubic field."""
    data = E.pari_curve().ellrankinit()
    if data.type() != 't_VEC' or len(data) != 3 or not 1 <= len(data[2]) <= 3:
        raise RuntimeError("unsupported PARI descent initialization")
    if certify:
        _certify_pari_descent(data)
    return data


def _certify_pari_descent(data):
    """Certify cubic class and unit groups before using PARI's rank bound."""
    # Reducible descent uses linear or quadratic nf structures, not bnf.
    fields = data[2]
    if len(fields) == 1 and fields[0].bnfcertify() != 1:
        raise RuntimeError("the descent number field could not be certified")


def _pari_rank_data(E, data, effort):
    """Return PARI's rank estimate and candidate points for descent data."""
    _lower, upper, _sha, points = data.ellrank(effort)
    return ZZ(upper), [E([QQ(P[0]), QQ(P[1])]) for P in points]


def _mwrank_search(E, first_limit=6, second_limit=4, saturation=-1):
    r"""
    Return eclib's rank upper bound and candidate points as integer triples.

    Each call uses a fresh eclib curve; ``(20, 12, 1000)`` gives the
    search and saturation limits used by :meth:`~sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field.gens`.
    """
    from sage.libs.eclib.interface import mwrank_EllipticCurve
    curve = mwrank_EllipticCurve(E.ainvs())
    curve.two_descent(verbose=False, first_limit=first_limit,
                      second_limit=second_limit)
    curve.saturate(bound=saturation)
    return (int(curve.rank_bound()),
            [[int(c) for c in P] for P in curve.gens()])


def _saturation(E, points, rank=None):
    r"""
    Return subgroup rank, saturation status, index, unresolved primes and points.

    Use a fresh eclib curve, saturating at all primes; if ``rank`` is
    supplied, first check that the subgroup has that rank.
    """
    from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil
    group = mwrank_MordellWeil(mwrank_EllipticCurve(E.ainvs()), verbose=False)
    group.process(points, saturation_bound=0)
    if rank is None:
        rank = int(group.rank())
    if group.rank() != rank:
        return int(group.rank()), False, 0, '', []
    try:
        ok, index, unsaturated = group.saturate(max_prime=-1)
    except (ArithmeticError, RuntimeError):
        ok, index, unsaturated = False, 0, 'unknown primes'
    return (rank, bool(ok), int(index), str(unsaturated),
            [[int(c) for c in P] for P in group.points()])


def _rank_zero_from_lseries(E, max_terms=2000000):
    r"""
    Prove rank zero by enclosing a nonzero central `L`-value.

    Use the series in [GJPST2009]_, Section 2.2, with Hasse's tail bound,
    and modularity and Kolyvagin's theorem [Kol1989]_ for the conclusion.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _rank_zero_from_lseries
        sage: _rank_zero_from_lseries(EllipticCurve([0,0,0,0,-1221]))
        True
        sage: _rank_zero_from_lseries(EllipticCurve('37a1'))
        False

    A tail that dominates the possible partial sum needs no coefficients::

        sage: from unittest.mock import patch
        sage: F = EllipticCurve('389a1')
        sage: with patch.object(F, 'anlist', side_effect=AssertionError('unneeded')):
        ....:     _rank_zero_from_lseries(F, max_terms=1)
        False
    """
    if E.root_number() != 1:
        return False
    R = RealBallField(128)
    conductor = R(E.conductor())
    q = (-2*R.pi()/conductor.sqrt()).exp()
    if not 1-q > 0:
        return False
    n = min(ZZ(max_terms), max(ZZ(64), _ceil_upper(
        conductor.sqrt()*(conductor.log()+20)/(2*R.pi()))))
    # The partial sum has absolute value at most 4*q*(1-q**n)/(1-q).
    # If the tail dominates this, its enclosure must contain zero.
    if 2*q**n >= 1:
        return False
    coefficients = E.anlist(n)
    polynomial = PolynomialRing(R, 't')(coefficients[1:]).integral()
    # Hasse and multiplicativity give |a_n| <= d(n)*sqrt(n) <= 2*n.
    # Thus the omitted part of 2*sum a_n*q^n/n has this absolute bound.
    tail = 4*q**(n+1)/(1-q)
    value = (2*polynomial(q)).add_error(tail)
    return 0 not in value


def _validate_mw_base(E, mw_base):
    """Validate supplied nontorsion points, or return ``None`` for ``'auto'``."""
    if isinstance(mw_base, str) and mw_base == 'auto':
        return None
    if not isinstance(mw_base, (list, tuple)):
        raise TypeError("mw_base must be a list or 'auto'")
    from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
    points = list(mw_base)
    if not all(isinstance(P, EllipticCurvePoint_field) and P.curve() is E
               for P in points):
        raise ValueError("points are not on the correct curve")
    if any(not P or P.has_finite_order() for P in points):
        raise ValueError("mw_base contains a torsion point")
    return points


@_with_proved_arithmetic
def certified_basis(E, mw_base='auto', verbose=False):
    r"""
    Return independent generators with certified full rank and saturation.

    Try cheap certificates and candidate searches before the database,
    deferring cubic class-group proofs for fields of discriminant above
    `10^{16}` until eclib has been tried.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import certified_basis
        sage: E = EllipticCurve('37a1'); P = E(0,0)
        sage: len(certified_basis(E, [P]))
        1
        sage: B = certified_basis(E, [2*P])
        sage: B[0] in [P, -P]
        True
        sage: certified_basis(E, [])
        Traceback (most recent call last):
        ...
        ValueError: mw_base does not have the certified rank
        sage: certified_basis(E, [0])
        Traceback (most recent call last):
        ...
        ValueError: points are not on the correct curve
        sage: F = EllipticCurve('389a1'); Q = F(0, 0)
        sage: certified_basis(F, [Q, 2*Q])
        Traceback (most recent call last):
        ...
        ValueError: mw_base is not independent

    An unproved PARI estimate is never used as a rank certificate::

        sage: from unittest.mock import Mock, patch
        sage: H = EllipticCurve([0,0,0,0,-10000]); Q = H(25,75)
        sage: U = H.isomorphism_to(H.minimal_model())(Q)
        sage: x, y = U.xy(); d = lcm(x.denominator(), y.denominator())
        sage: coordinates = [ZZ(d*x), ZZ(d*y), ZZ(d)]
        sage: target = 'sage.schemes.elliptic_curves.integral_points.'
        sage: with (patch.object(H, 'database_attributes', side_effect=LookupError()),
        ....:       patch.object(H, 'gens_certain', return_value=False),
        ....:       patch(target + '_rank_zero_from_lseries', return_value=False),
        ....:       patch(target + '_pari_rank_data', return_value=(ZZ(1), [Q])),
        ....:       patch(target + '_certify_pari_descent', side_effect=RuntimeError('unproved')),
        ....:       patch(target + '_mwrank_search', return_value=(2, [coordinates]))):
        ....:     certified_basis(H)
        Traceback (most recent call last):
        ...
        RuntimeError: unable to certify a Mordell-Weil basis

    Large fields try eclib before the class-group certificate::

        sage: field = Mock(); field.idealnorm.return_value = ZZ(10)^17
        sage: with (patch.object(H, 'database_attributes', side_effect=LookupError()),
        ....:       patch.object(H, 'gens_certain', return_value=False),
        ....:       patch(target + '_rank_zero_from_lseries', return_value=False),
        ....:       patch(target + '_pari_descent', return_value={2: [field]}),
        ....:       patch(target + '_pari_rank_data', return_value=(ZZ(1), [Q])),
        ....:       patch(target + '_certify_pari_descent', side_effect=AssertionError('unneeded')),
        ....:       patch(target + '_mwrank_search', return_value=(1, [coordinates])) as search):
        ....:     B = certified_basis(H)
        sage: len(B), search.call_args.args[1:]
        (1, (20, 12, 1000))

    A cheap certificate avoids opening the database::

        sage: F = EllipticCurve([1, -1, 0, -79, 289])
        sage: with patch.object(F, 'database_attributes', side_effect=AssertionError('unneeded')):
        ....:     len(certified_basis(F))
        4

    The eligible Heegner route still precedes an eclib search::

        sage: F = EllipticCurve([0, -1, 0, 3107360, 14371118908])
        sage: with (patch.object(F, 'database_attributes', side_effect=AssertionError('unneeded')),
        ....:       patch(target + '_mwrank_search', side_effect=AssertionError('unneeded'))):
        ....:     len(certified_basis(F))
        1

    Uncertified supplied points report how to continue::

        sage: with (patch.object(E, 'database_attributes', side_effect=LookupError()),
        ....:       patch(target + '_pari_descent', side_effect=RuntimeError()),
        ....:       patch(target + '_mwrank_search', side_effect=RuntimeError())):
        ....:     E.integral_points([P], proof=True)
        Traceback (most recent call last):
        ...
        RuntimeError: unable to certify a Mordell-Weil rank; pass proof=False to use the given points without a proof
    """
    points = _validate_mw_base(E, mw_base)
    automatic = points is None
    minimal = E.minimal_model()
    to_E = minimal.isomorphism_to(E)
    # Descent gives an upper bound, certified once independent points
    # meet it. Defer the exact database rank until cheap searches fail.
    rank = None
    exact_rank = False

    def update(upper, found):
        nonlocal rank, points
        if upper is not None and not exact_rank:
            rank = upper if rank is None else min(rank, upper)
        if automatic and len(found) > len(points):
            points = found

    def certified():
        return rank is not None and len(points) == rank

    if automatic:
        points = []
        # Generators cached by gens(proof=True) are candidates only.
        try:
            if E.gens_certain():
                points = list(E.gens(proof=True))
        except RuntimeError:
            pass
    lseries_tried = False
    if rank is None and not points and E.conductor() <= 10**8:
        lseries_tried = True
        if _rank_zero_from_lseries(E):
            return []
    descent = None
    certificate_state = None
    pari_upper = None
    pari_certified = False
    pari_failed = False
    large_field = False
    full_search_tried = False

    def certify_pari(force=False):
        nonlocal pari_certified, pari_failed
        if (certified() or exact_rank or descent is None or pari_upper is None
                or pari_failed or len(points) < pari_upper):
            return
        if large_field and not force:
            return
        try:
            with _pari_seed(certificate_state):
                _certify_pari_descent(descent)
            pari_certified = True
            update(pari_upper, [])
        except (ArithmeticError, RuntimeError):
            pari_failed = True

    # Candidate points need no class-group certificate. PARI's upper
    # bound is used only after its cubic class and unit groups are proved.
    if not certified():
        try:
            with _pari_seed(1):
                descent = _pari_descent(E, certify=False)
                certificate_state = pari.getrand()
                fields = descent[2]
                pari_certified = len(fields) != 1
                if automatic and not exact_rank and len(fields) == 1:
                    field = fields[0]
                    large_field = field.idealnorm(field.nf_get_diff()) > 10**16
                pari_upper, found = _pari_rank_data(E, descent, 0)
                update(pari_upper if pari_certified else None, found)
        except (ArithmeticError, RuntimeError):
            pass
    certify_pari()
    if not lseries_tried and not points and not exact_rank and not certified():
        if _rank_zero_from_lseries(E):
            rank = ZZ(0)
            if automatic:
                points = []
    heegner_tried = False

    def heegner():
        nonlocal points, heegner_tried
        expected_rank = rank if rank is not None else pari_upper
        if (automatic and expected_rank == 1 and not points and not heegner_tried
                and E.conductor() <= 10**7):
            heegner_tried = True
            try:
                coordinates = E.pari_curve().ellheegner()
                P = E([QQ(coordinates[0]), QQ(coordinates[1])])
                if P and not P.has_finite_order():
                    points = [P]
            except (ArithmeticError, RuntimeError, TypeError, ValueError, IndexError):
                pass

    heegner()
    certify_pari()
    # Retain the source order for large fields: their full eclib search
    # precedes the class-group certificate and the short eclib search.
    if large_field and not certified():
        full_search_tried = True
        try:
            upper, found = _mwrank_search(minimal, 20, 12, 1000)
            update(ZZ(upper), [to_E(minimal(P)) for P in found])
        except (ArithmeticError, RuntimeError):
            pass
        certify_pari(force=True)
    if not certified():
        try:
            upper, found = _mwrank_search(minimal)
            update(ZZ(upper), [to_E(minimal(P)) for P in found])
        except (ArithmeticError, RuntimeError):
            pass
    certify_pari(force=full_search_tried)
    if not certified():
        try:
            attributes = E.database_attributes()
            rank = ZZ(attributes['rank'])
            exact_rank = True
            if automatic and 'gens' in attributes:
                points = [to_E(minimal(P)) for P in attributes['gens']]
        except LookupError:
            pass
    if rank == 0:
        if not automatic and points:
            raise ValueError("mw_base does not have the certified rank")
        return []
    if exact_rank and not automatic and len(points) != rank:
        raise ValueError("mw_base does not have the certified rank")
    if exact_rank:
        # A database rank of one can make the original Heegner route
        # eligible even when the preceding descent bound did not.
        heegner()
    if automatic and descent is not None and not certified():
        try:
            with _pari_seed(1):
                upper, found = _pari_rank_data(E, descent, 1)
            pari_upper = upper if pari_upper is None else min(pari_upper, upper)
            update(pari_upper if pari_certified else None, found)
        except (ArithmeticError, RuntimeError):
            pass
    certify_pari(force=full_search_tried)
    if automatic and not certified() and not full_search_tried:
        try:
            upper, found = _mwrank_search(minimal, 20, 12, 1000)
            update(ZZ(upper), [to_E(minimal(P)) for P in found])
        except (ArithmeticError, RuntimeError):
            pass
    certify_pari(force=True)
    proof_hint = "" if automatic else "; pass proof=False to use the given points without a proof"
    if rank is None:
        raise RuntimeError("unable to certify a Mordell-Weil rank" + proof_hint)
    if len(points) != rank:
        if automatic or (not exact_rank and len(points) < rank):
            raise RuntimeError("unable to certify a Mordell-Weil basis" + proof_hint)
        raise ValueError("mw_base does not have the certified rank")
    if not rank:
        return []

    phi = E.isomorphism_to(minimal)
    projective = []
    for P in points:
        if not P or P.has_finite_order():
            raise ValueError("mw_base contains a torsion point")
        x, y = phi(P).xy()
        d = x.denominator().lcm(y.denominator())
        projective.append([ZZ(d*x), ZZ(d*y), ZZ(d)])
    try:
        independent, ok, index, unsaturated, saturated = _saturation(
            minimal, projective, int(rank))
    except RuntimeError as error:
        raise RuntimeError("saturation could not be certified") from error
    if independent != rank:
        raise ValueError("mw_base is not independent")
    if not ok:
        raise RuntimeError("saturation could not be certified at %s" % unsaturated)
    basis = [to_E(minimal(P)) for P in saturated]
    return _lll_basis(E, basis)


def _lll_basis(E, basis):
    """Reduce a basis, checking the integral change of basis exactly."""
    rank = len(basis)
    if rank > 1:
        try:
            reduced, U = E.lll_reduce(basis)
            if (abs(U.det()) == 1 and all(u in ZZ for u in U.list())
                    and all(sum((U[j,i]*basis[j] for j in range(rank)), E(0)) == reduced[i]
                            for i in range(rank))):
                basis = reduced
        except (ArithmeticError, ValueError):
            pass
    return basis


def _unproved_gens(E):
    """Return Sage's unproved generators without caching a search on ``E``."""
    from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field
    fresh = EllipticCurve_rational_field(list(E.ainvs()))
    result = []
    for P in fresh.gens(proof=False):
        x, y = P.xy()
        d = x.denominator().lcm(y.denominator())
        result.append([int(d*x), int(d*y), int(d)])
    return result


@_with_proved_arithmetic
def _uncertified_basis(E, mw_base='auto', verbose=False):
    r"""
    Return candidate generators, saturated when possible.

    Warn once if PARI's lower rank estimate or descent parity requires
    more independent points, or if full saturation is not certified;
    the parity test assumes finite Tate--Shafarevich group.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _uncertified_basis
        sage: E = EllipticCurve('37a1'); P = E(0, 0)
        sage: B = _uncertified_basis(E, [3*P])
        sage: B[0] in [P, -P]
        True

    A shortfall gives one warning and retains the subgroup found::

        sage: import warnings
        sage: F = EllipticCurve('389a1')
        sage: with warnings.catch_warnings(record=True) as caught:
        ....:     warnings.simplefilter('always')
        ....:     B = _uncertified_basis(F, [F(0, 0)])
        sage: len(B), len(caught), 'requires at least 2' in str(caught[0].message)
        (1, 1, True)

    Failed saturation also retains supplied points with one warning::

        sage: from unittest.mock import patch
        sage: target = 'sage.schemes.elliptic_curves.integral_points._saturation'
        sage: with warnings.catch_warnings(record=True) as caught:
        ....:     warnings.simplefilter('always')
        ....:     with patch(target, side_effect=RuntimeError('failed')):
        ....:         B = _uncertified_basis(E, [P])
        sage: B == [P], len(caught), 'saturation' in str(caught[0].message)
        (True, 1, True)
    """
    points = _validate_mw_base(E, mw_base)
    automatic = points is None
    minimal = E.minimal_model()
    to_E = minimal.isomorphism_to(E)
    sources = [] if automatic else [points]
    lower = ZZ(0)
    upper = parity = None
    search_failed = False
    if automatic:
        try:
            E.gens_certain()
            sources.append(list(E.gens(proof=False)))
        except RuntimeError:
            pass
        try:
            attributes = E.database_attributes()
            if 'gens' in attributes:
                sources.append([to_E(minimal(P)) for P in attributes['gens']])
        except LookupError:
            pass
    try:
        with _pari_seed(1):
            lower, upper, _sha, found = E.pari_curve().ellrank(1)
        lower, upper = ZZ(lower), ZZ(upper)
        # PARI's bound is C - T - s, with s even. Under finite Sha,
        # dim Sha[2] is even as well, giving this rank parity.
        parity = upper % 2
        if automatic:
            sources.append([E([QQ(P[0]), QQ(P[1])]) for P in found])
    except (ArithmeticError, RuntimeError):
        search_failed = True
    if automatic:
        if upper is None or max(map(len, sources), default=0) < upper:
            try:
                found = _unproved_gens(minimal)
                sources.append([to_E(minimal(P)) for P in found])
            except (ArithmeticError, RuntimeError):
                search_failed = True
        points = list(dict.fromkeys(P for source in sources for P in source
                                    if P and not P.has_finite_order()))
    ok = True
    basis = points
    if points:
        phi = E.isomorphism_to(minimal)
        projective = []
        for P in points:
            x, y = phi(P).xy()
            d = x.denominator().lcm(y.denominator())
            projective.append([ZZ(d*x), ZZ(d*y), ZZ(d)])
        try:
            rank, ok, _index, _unsaturated, saturated = _saturation(minimal, projective)
            if not automatic and rank != len(points):
                raise ValueError("mw_base is not independent")
            basis = [to_E(minimal(P)) for P in saturated]
        except (ArithmeticError, RuntimeError):
            # Individual search sources give independent lists. The height
            # bound also checks positive definiteness before enumeration.
            basis = max(sources, key=len)
            ok = False
    reasons = []
    expected = lower
    if parity is not None and len(basis) % 2 != parity:
        expected = max(expected, len(basis) + 1)
    if len(basis) < expected:
        reasons.append("found %s independent points, but PARI's rank estimate requires at least %s"
                       % (len(basis), expected))
    if not ok:
        reasons.append("saturation could not be certified")
    if search_failed and upper is None:
        reasons.append("generator searches did not determine the rank")
    if reasons:
        warnings.warn("Integral-point list may be incomplete: %s." % '; '.join(reasons),
                      stacklevel=3)
    return _lll_basis(E, basis)


def _height_difference(E, R):
    r"""Return an integer bounding the difference of full and naive heights."""
    j = E.j_invariant()
    b = E.b2()/12
    h_j = R(max(abs(j.numerator()), j.denominator())).log()
    # Silverman's bound, enlarged to avoid its decimal constants.
    value = (h_j/12 + R(abs(E.discriminant())).log()/6
             + R(max(QQ(1), abs(j))).log()/6
             + R(max(QQ(1), abs(b))).log() + R(2).log() + 3)
    return _ceil_upper(value)


def _positive_definite(H):
    """Verify positive definiteness by interval LDL decomposition."""
    n = len(H)
    L = [[0 for j in range(n)] for i in range(n)]
    diagonal = []
    for j in range(n):
        d = H[j][j] - sum(L[j][k]**2 * diagonal[k] for k in range(j))
        if not d > 0:
            return False
        diagonal.append(d)
        for i in range(j+1, n):
            L[i][j] = (H[i][j] - sum(L[i][k]*L[j][k]*diagonal[k]
                                     for k in range(j))) / d
    return True


def _height_matrix_bound(E, basis, B, R, include_matrix=False):
    r"""Enclose canonical heights and certify a positive lower eigenvalue."""
    r = len(basis)
    points = list(basis)
    pairs = [(i, j) for i in range(r) for j in range(i)]
    points.extend(basis[i] + basis[j] for i, j in pairs)
    multiples = [8*P for P in points]
    for n in range(3, 11):
        scale = ZZ(4)**n
        heights = []
        for P in multiples:
            if not P:
                heights.append(R(0))
            else:
                x = P[0]
                naive = R(max(abs(x.numerator()), x.denominator())).log()
                heights.append((naive/scale).add_error(R(B)/scale))
        H = [[R(0) for j in range(r)] for i in range(r)]
        for i in range(r):
            H[i][i] = heights[i]
        for h, (i, j) in zip(heights[r:], pairs):
            H[i][j] = H[j][i] = (h - heights[i] - heights[j])/2
        if r == 1:
            lower = H[0][0].lower().exact_rational()
            if lower > 0 and lower > 3*H[0][0].mid()/4:
                return (lower, heights[:r], H) if include_matrix else (lower, heights[:r])
        elif _positive_definite(H):
            mid = matrix(RDF, [[float(x.mid()) for x in row] for row in H])
            estimate = min(x.real() for x in mid.eigenvalues())
            uncertainty = max(sum(_upper(abs(x-R(x.mid()))) for x in row)
                              for row in H)
            if uncertainty > QQ(estimate)/32:
                multiples = [2*P for P in multiples]
                continue
            candidate = QQ(estimate)*QQ(7)/8
            for unused in range(16):
                if candidate > 0:
                    shifted = [[H[i][j] - (candidate if i == j else 0)
                                for j in range(r)] for i in range(r)]
                    if _positive_definite(shifted):
                        return (candidate, heights[:r], H) if include_matrix else (candidate, heights[:r])
                candidate /= 2
        multiples = [2*P for P in multiples]
    raise ArithmeticError("a positive height bound could not be certified")


class _Logarithms:
    """Enclose real logarithms and periods using Carlson's integral."""

    def __init__(self, E, basis):
        self.E = E
        self.basis = basis
        f = PolynomialRing(QQ, 'x')([-E.c6()/864, -E.c4()/48, 0, 1])
        self.roots = f.roots(QQbar, multiplicities=False)
        self.real = sorted(z for z in self.roots if z.imag() == 0)
        self.root = self.real[-1]

    def enclosures(self, precision):
        E = self.E
        C = ComplexBallField(precision)
        R = RealBallField(precision)
        e1 = C(self.root)
        e2, e3 = [C(z) for z in self.roots if z != self.root]
        w = (2*C(0).elliptic_rf(e1-e2, e1-e3)).real()
        if E.discriminant() > 0:
            a, b, c = self.real
            v = (2*C(0).elliptic_rf(C(b-a), C(c-a))).real()
        else:
            v = C(0).elliptic_rf(e2-e1, e3-e1).real()
        if not (w > 0 and v > 0):
            raise ArithmeticError("period enclosures are not positive")
        logarithms = []
        for P in self.basis:
            Q = 2*P
            s = Q[0] + E.b2()/12
            if not R(s) - e1.real() > 0:
                raise ArithmeticError("the real logarithm is not separated from a branch point")
            u = C(C(s)-e1).elliptic_rf(C(s)-e2, C(s)-e3).real()
            # The derivative of wp is negative on (0, w/2).
            if Q[1] + (E.a1()*Q[0] + E.a3())/2 > 0:
                u = -u
            logarithms.append(u/w)
        root_bound = max(_upper(abs(C(z))) for z in self.roots)
        return logarithms, w, w*v, root_bound


def _shortest_vector_squared(basis, orthogonal=None):
    r"""
    Return the exact squared minimum of an integer row lattice of rank at most eight.

    Use Fincke--Pohst enumeration [FP1985]_ with rational Gram--Schmidt
    coefficients and outward integer bounds; no floating-point test prunes a vector.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import _shortest_vector_squared
        sage: _shortest_vector_squared(matrix(ZZ, [[2, 1], [1, 2]]))
        2
        sage: _shortest_vector_squared(matrix(ZZ, [[3]]))
        9
        sage: _shortest_vector_squared(identity_matrix(ZZ, 8))
        1
    """
    dimension = basis.nrows()
    if not 1 <= dimension <= 8:
        raise ValueError("the lattice rank must be between one and eight")
    if orthogonal is None:
        orthogonal = basis.gram_schmidt()[0]
    stars = orthogonal.rows()
    lengths = [row.dot_product(row) for row in stars]
    if any(d <= 0 for d in lengths):
        raise ArithmeticError("the lattice basis is not independent")
    rows = basis.rows()
    mu = [[rows[j].dot_product(stars[i])/lengths[i] for i in range(j)]
          for j in range(dimension)]
    best = min(row.dot_product(row) for row in rows)
    coefficients = [ZZ(0)] * dimension

    def visit(i, partial, nonzero):
        nonlocal best
        if i < 0:
            if nonzero:
                best = min(best, partial)
            return
        center = sum(mu[j][i]*coefficients[j] for j in range(i+1, dimension))
        radius = ZZ(isqrt(int(((best-partial)/lengths[i]).ceil())) + 1)
        first, last = (-center-radius).ceil(), (-center+radius).floor()
        if not nonzero:
            first = max(first, ZZ(0))
        for n in range(first, last+1):
            total = partial + lengths[i]*(n+center)**2
            if total < best:
                coefficients[i] = ZZ(n)
                visit(i-1, total, nonzero or n != 0)

    visit(dimension-1, QQ(0), False)
    return ZZ(best)


def coefficient_bound(E, basis, verbose=False, with_height=False):
    r"""
    Return a certified bound for free coefficients of integral points.

    The subgroup includes all torsion translates. With ``with_height=True``,
    also return the enclosed height matrix and a rational height upper bound.

    TESTS:

    The point with x-coordinate 20829482 is a seventh multiple::

        sage: from sage.schemes.elliptic_curves.integral_points import certified_basis, coefficient_bound
        sage: E = EllipticCurve([1,-1,1,-11,16])
        sage: P = E(2,-2)
        sage: (7*P)[0]
        20829482
        sage: coefficient_bound(E, certified_basis(E, [P])) >= 7
        True
    """
    r = len(basis)
    if not r:
        return (ZZ(0), [], QQ(0)) if with_height else ZZ(0)
    R = RealBallField(128)
    B = _height_difference(E, R)
    eigenvalue, heights, height_matrix = _height_matrix_bound(E, basis, B, R, include_matrix=True)
    geometry = _Logarithms(E, basis)
    precision = 128
    while True:
        try:
            logs, w, area, root_bound = geometry.enclosures(precision)
            break
        except ArithmeticError:
            precision *= 2
            if precision > 65536:
                raise ArithmeticError("real logarithms could not be enclosed")
    R = w.parent()
    # Outside this interval, the elementary elliptic-integral bound applies.
    shift = abs(E.b2()/12)
    z = max(2*root_bound, shift, QQ(1))
    small_height = R(max(QQ(1), z+shift)).log()+B
    small_bound = _ceil_upper((small_height/eigenvalue).sqrt())
    height_limit = None
    log_upper = (R(4)*R(2).sqrt()).log() + (R(B)+R(2).log())/2

    # David's theorem, k = r+1 including the real period.  The integer
    # constant below exceeds the constant in the final Theorem 2.1 by
    # more than 10^6. This also absorbs log(e), e <= 12, when dividing
    # the integer-coefficient form e*2*u(P) by the torsion exponent e.
    g2, g3 = E.c4()/12, E.c6()/216
    d = g2.denominator().lcm(g3.denominator())
    h = max(ZZ(1), _ceil_upper(R(max(d, abs(d*g2), abs(d*g3))).log()),
            _ceil_upper(R(max(abs(E.j_invariant().numerator()),
                             E.j_invariant().denominator())).log()))
    A = max(ZZ(1), h, _ceil_upper(3*R.pi()*w**2/area),
            *(_ceil_upper(8*height) for height in heights))
    torsion_exponent = E.torsion_subgroup().exponent()
    k = r+1
    constant = (ZZ(3)*10**(6*k+12)*4**(2*(k+1)**2)
                * (k+2)**(2*k*k+13*k+30))
    offset = A + 3*h + _ceil_upper(R(torsion_exponent*(r+2)).log())
    power = ZZ(k+1)
    while True:
        s = power*R(10).log()
        beta = s + offset
        rhs = constant*R(A)**k*(beta+1)*(beta+h+1)**(k+1)+log_upper
        if R(eigenvalue)*R(10)**(2*power)/2 > rhs:
            break
        power *= 2
    N = max(ZZ(10)**power, small_bound)
    if verbose:
        print("Height difference bound:", B)
        print("Initial coefficient bound:", N)

    for unused in range(32):
        if N <= small_bound:
            break
        # Start near the coefficient scale when the exact minimum is
        # available. A failed positivity test only enlarges the scale.
        scale = ((torsion_exponent*N)**(r+1) if r+1 <= 8
                 else (10*torsion_exponent*(r+1)*N)**(r+1))
        for attempt in range(32):
            needed = scale.nbits()+64
            if needed > precision:
                precision = needed
                logs, w, area, root_bound = geometry.enclosures(precision)
            R = w.parent()
            rounded = [ZZ((scale*x).mid().round()) for x in logs]
            errors = [_upper(abs(scale*x-a)) for x, a in zip(logs, rounded)]
            lattice = matrix(ZZ, r+1, r+1, 1)
            for i, a in enumerate(rounded):
                lattice[i, r] = a
            lattice[r, r] = scale
            reduced = lattice.LLL()
            gs = reduced.gram_schmidt()[0]
            shortest_squared = (_shortest_vector_squared(reduced, gs) if r+1 <= 8
                                else min(row.dot_product(row) for row in gs.rows()))
            gap = QQ(shortest_squared)/torsion_exponent**2 - r*N**2
            if gap > 0:
                lower = (R(gap).sqrt()-N*sum(errors))/scale
                if lower > 0:
                    # Re-enclose the analytic upper constant at this precision.
                    log_K = ((R(4)*R(2).sqrt()).log()
                             + (R(B)+R(2).log())/2-w.log())
                    analytic_height = 2*(log_K-lower.log())
                    limit = max(_upper(small_height), _upper(analytic_height))
                    height_limit = limit if height_limit is None else min(height_limit, limit)
                    square = analytic_height/eigenvalue
                    if square < 0:
                        new_bound = small_bound
                    else:
                        new_bound = max(small_bound, _ceil_upper(square.sqrt()))
                    break
            scale *= 16
        else:
            raise ArithmeticError("lattice reduction could not certify a bound")
        if new_bound >= N:
            break
        N = ZZ(new_bound)
    else:
        raise ArithmeticError("coefficient reduction did not terminate")
    if verbose:
        print("Coefficient bound:", N)
    if with_height:
        if height_limit is None:
            height_limit = sum(_upper(abs(x))*N**2 for row in height_matrix for x in row)
        margin = QQ(1)/2**32
        height_limit = height_limit*(1+margin)+margin
        grid = ZZ(2)**64
        return ZZ(N), height_matrix, QQ((height_limit*grid).ceil())/grid
    return ZZ(N)


def _height_ellipsoid(height_matrix, height_limit):
    r"""Return an exact LDL ellipsoid containing the certified height ellipsoid."""
    r = len(height_matrix)
    grid = ZZ(2)**64
    mid = [[QQ((x.mid().exact_rational()*grid).round())/grid for x in row]
           for row in height_matrix]
    error = max(sum(_upper(abs(x-x.parent()(m))) for x, m in zip(row, centers))
                for row, centers in zip(height_matrix, mid))
    error = QQ((error*grid).ceil())/grid
    # A symmetric error matrix has operator norm at most its maximum
    # absolute row sum. Thus mid-error*I is below the true height matrix.
    A = [[mid[i][j] - (error if i == j else 0) for j in range(r)] for i in range(r)]
    L = [[QQ(0) for j in range(r)] for i in range(r)]
    diagonal = []
    for j in range(r):
        d = A[j][j] - sum(L[j][k]**2*diagonal[k] for k in range(j))
        if d <= 0:
            return None
        diagonal.append(d)
        for i in range(j+1, r):
            L[i][j] = (A[i][j] - sum(L[i][k]*L[j][k]*diagonal[k] for k in range(j)))/d
    mu = [[L[r-1-k][r-1-j] for k in range(j)] for j in range(r)]
    return list(reversed(diagonal)), mu, QQ(height_limit)


def bounded_integral_x_coordinates(E, basis, N, height_data=None):
    r"""
    Enumerate integral x-coordinates with coefficients bounded by ``N``.

    With ``height_data=(matrix, bound)``, intersect the box with its certified
    height ellipsoid. Include all torsion translates, with exact fallbacks.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import bounded_integral_x_coordinates
        sage: E = EllipticCurve('389a1')
        sage: sorted(bounded_integral_x_coordinates(E, E.gens(), 4))
        [-2, -1, 0, 1, 3, 4, 6, 39, 133, 188]
        sage: R = RealBallField(128)
        sage: H = [[R(1), R(0)], [R(0), R(1)]]
        sage: sorted(bounded_integral_x_coordinates(E, [E(0, 0), E(-1, 1)], 4, (H, 1)))
        [-1, 0]
        sage: E = EllipticCurve('11a1')
        sage: sorted(bounded_integral_x_coordinates(E, [], 0))
        [5, 16]
        sage: bounded_integral_x_coordinates(E, [], -1)
        Traceback (most recent call last):
        ...
        ValueError: the coefficient bound must be nonnegative

    Wide and non-finite enclosures fall back to exact group operations::

        sage: from unittest.mock import patch
        sage: E = EllipticCurve([0,0,0,-25,0]); P = E(-4,6)
        sage: target = 'sage.schemes.elliptic_curves.integral_points.RealBallField'
        sage: R = RealBallField(128)
        sage: with patch(target, return_value=RealBallField(2)):
        ....:     sorted(bounded_integral_x_coordinates(E, [P,P], 2))
        [-5, -4, 0, 5, 45]
        sage: with patch(target, return_value=lambda q: R('nan')):
        ....:     sorted(bounded_integral_x_coordinates(E, [P,P], 2))
        [-5, -4, 0, 5, 45]
    """
    N = ZZ(N)
    if N < 0:
        raise ValueError("the coefficient bound must be nonnegative")
    torsion = E.torsion_points()
    zero = E(0)
    xs = set()

    def use(P):
        if P and P[0].is_integral() and P[1].is_integral():
            xs.add(ZZ(P[0]))

    r = len(basis)
    ellipsoid = _height_ellipsoid(*height_data) if height_data is not None and r else None
    if ellipsoid is not None:
        diagonal, mu, height_limit = ellipsoid
        basis = list(reversed(basis))
        if r == 1:
            N = min(N, ZZ(isqrt(int((height_limit/diagonal[0]).floor()))))
    if r <= 1:
        P = zero
        for n in range(N+1 if r else 1):
            for T in torsion:
                use(P+T)
            if r:
                P += basis[0]
        return xs

    R = RealBallField(128)
    a1, a2, a3 = map(R, E.a_invariants()[:3])
    has_a1, has_a2 = bool(E.a1()), bool(E.a2())

    def enclose(P):
        return (R(P[0]), R(P[1])) if P else None

    def add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        denominator = x2-x1
        if denominator.contains_zero():
            raise ArithmeticError("ambiguous ball addition")
        slope = (y2-y1)/denominator
        x3 = slope**2-x1-x2
        if has_a1:
            x3 += a1*slope
        if has_a2:
            x3 -= a2
        y3 = -y1-a1*x3-a3+slope*(x1-x3)
        return x3, y3

    multiples = []
    for P in basis:
        row = {0: zero}
        Q = zero
        for n in range(1, N+1):
            Q += P
            row[n], row[-n] = Q, -Q
        multiples.append({n: (Q, enclose(Q)) for n, Q in row.items()})
    # Absorb torsion into the final summand. At the last level only its
    # x-coordinate is needed, so no final y-coordinate enclosure is formed.
    last = {}
    for n, (Q, ball) in multiples[-1].items():
        translates = [Q+T if T else Q for T in torsion]
        last[n] = [(P, enclose(P)) for P in translates]

    def choices(j, negative, coefficients, partial):
        last = N if negative else ZZ(0)
        if ellipsoid is None:
            for n in range(-N, last+1):
                yield n, partial
            return
        center = sum(a*n for a, n in zip(mu[j], coefficients))
        radius = ZZ(isqrt(int(((height_limit-partial)/diagonal[j]).ceil())) + 1)
        first = max(-N, (-center-radius).ceil())
        last = min(last, (-center+radius).floor())
        for n in range(first, last+1):
            total = partial + diagonal[j]*(n+center)**2
            if total <= height_limit:
                yield n, total

    def visit(j, enclosure, terms, negative, coefficients, partial):
        if j < r-1:
            for n, total in choices(j, negative, coefficients, partial):
                Q, ball = multiples[j][n]
                new_terms = terms + (Q,)
                try:
                    new_ball = add(enclosure, ball)
                except ArithmeticError:
                    new_ball = enclose(sum(new_terms, zero))
                visit(j+1, new_ball, new_terms, negative or n < 0,
                      coefficients + (n,), total)
            return
        exact = None
        for n, total in choices(j, negative, coefficients, partial):
            for Q, ball in last[n]:
                if enclosure is None:
                    if ball is None:
                        continue
                    x = ball[0]
                elif ball is None:
                    x = enclosure[0]
                else:
                    x1, y1 = enclosure
                    x2, y2 = ball
                    denominator = x2-x1
                    if denominator.contains_zero():
                        x = None
                    else:
                        slope = (y2-y1)/denominator
                        x = slope**2-x1-x2
                        if has_a1:
                            x += a1*slope
                        if has_a2:
                            x -= a2
                # An ambiguous or non-finite enclosure must be checked
                # exactly. Cache the exact prefix across final summands.
                if x is None or not x.is_finite() or x.contains_integer():
                    if exact is None:
                        exact = sum(terms, zero)
                    use(exact+Q)

    visit(0, None, (), False, (), QQ(0))
    return xs


def integral_points(E, mw_base='auto', both_signs=False, verbose=False, proof=None):
    r"""
    Compute integral points, proving completeness when ``proof`` requires it.

    See :meth:`~sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field.integral_points`
    for the arguments and completeness condition.

    TESTS::

        sage: from sage.schemes.elliptic_curves.integral_points import integral_points
        sage: E = EllipticCurve([0,0,0,-16,0])
        sage: [P[0] for P in integral_points(E)]
        [-4, 0, 4]
        sage: points = integral_points(E, both_signs=True)
        sage: all(-P in points for P in points)
        True

    Curves constructed directly from the class are supported::

        sage: from sage.schemes.elliptic_curves.ell_rational_field import EllipticCurve_rational_field
        sage: E = EllipticCurve_rational_field([0,1,1,-2,0])
        sage: [P[0] for P in integral_points(E)]
        [-2, -1, 0, 1, 3, 4, 6, 39, 133, 188]

    Automatic generator searches do not require a prior call to ``gens``::

        sage: E = EllipticCurve([0,0,0,0,-5134])
        sage: integral_points(E)
        []
        sage: generators = E.gens(proof=True)
        sage: integral_points(E)
        []
    """
    if not E.is_integral():
        raise ValueError("integral_points() can only be called on an integral model")
    proof = get_flag(proof, "elliptic_curve")
    with _proved_arithmetic():
        if proof:
            basis = certified_basis(E, mw_base, verbose)
        else:
            basis = _uncertified_basis(E, mw_base, verbose)
        if verbose:
            print("Mordell-Weil basis:", basis)
        N, height_matrix, height_limit = coefficient_bound(E, basis, verbose, with_height=True)
        xs = bounded_integral_x_coordinates(E, basis, N, (height_matrix, height_limit))
    if both_signs:
        result = [P for x in xs for P in E.lift_x(x, all=True)]
    else:
        result = [E.lift_x(x) for x in xs]
    if verbose:
        print("Found %s integral points." % len(result))
    return sorted(result)
