import pytest


def _polynomial_ring():
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ

    return PolynomialRing(QQ, 'x')


def test_equivalent_approximations_have_the_same_limit():
    from sage.rings.infinity import infinity
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.gauss_valuation import GaussValuation
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = _polynomial_ring()
    x = R.gen()
    vK = QQ.valuation(2)
    gauss = GaussValuation(R, vK)

    a = gauss.augmentation(x, 1)
    b = gauss.augmentation(x + 2, 1)
    assert a is not b
    assert a >= b >= a
    assert LimitValuation(a, x) is LimitValuation(b, x)

    G = x**2 + 1
    approximant = vK.mac_lane_approximants(
        G, require_incomparability=True)[0]
    limit = LimitValuation(gauss, G)
    assert limit is LimitValuation(approximant, G)
    assert limit is LimitValuation(
        approximant.augmentation(G, infinity), G)


def test_later_approximation_of_an_infinite_chain_has_the_same_limit():
    from sage.rings.function_field.constructor import FunctionField
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.limit_valuation import LimitValuation

    K = FunctionField(QQ, 't')
    t = K.gen()
    R = PolynomialRing(K, 'y')
    y = R.gen()
    G = y**2 - t
    vK = K.valuation(t)
    approximant = vK.mac_lane_approximants(
        G, require_incomparability=True)[0]
    later = approximant.mac_lane_step(G)[0]

    assert LimitValuation(approximant, G) is LimitValuation(later, G)


def test_composite_polynomial_has_a_canonical_support():
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = _polynomial_ring()
    x = R.gen()
    vK = QQ.valuation(2)
    G = (x**2 + 7) * (x**2 + 9)
    approximants = vK.mac_lane_approximants(
        G, require_incomparability=True)
    limits = [LimitValuation(v, G) for v in approximants]

    for v, limit in zip(approximants, limits):
        support, = [factor.monic() for factor, _ in G.factor()
                    if not v.is_equivalence_unit(factor)]
        assert limit._G == support
        assert limit is LimitValuation(v, support)

    assert all(limits[i] is not limits[j]
               for i in range(len(limits))
               for j in range(i))


def test_canonical_support_uses_incomparable_approximants():
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.gauss_valuation import GaussValuation
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = _polynomial_ring()
    x = R.gen()
    vK = QQ.valuation(2)
    base = GaussValuation(R, vK).augmentation(x - 3, 3)
    support = x**2 + 7
    G = (x + 1) * support
    approximants = vK.mac_lane_approximants(
        support, require_incomparability=True)

    limit = LimitValuation(base, G)
    assert limit._G == support
    assert limit._initial_approximation is vK.mac_lane_approximant(
        support, base, approximants=approximants)
    limits = [LimitValuation(approximant, support)
              for approximant in approximants]
    for i, left in enumerate(limits):
        for j, right in enumerate(limits):
            assert (left >= right) == (i == j)


def test_ambiguous_approximation_is_rejected():
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.gauss_valuation import GaussValuation
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = _polynomial_ring()
    x = R.gen()

    gauss = GaussValuation(R, QQ.valuation(2))
    G = (x**2 + 7) * (x**2 + 9)
    with pytest.raises(ValueError, match="single out one irreducible factor"):
        LimitValuation(gauss, G)

    gauss = GaussValuation(R, QQ.valuation(5))
    with pytest.raises(ValueError, match="does not approximate a unique extension"):
        LimitValuation(gauss, x**2 + 1)


def test_unchecked_canonical_key_uses_the_same_cache_entry():
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = _polynomial_ring()
    x = R.gen()
    vK = QQ.valuation(2)

    G = x**3 + 2
    approximant = vK.mac_lane_approximants(
        G, require_incomparability=True)[0]
    unchecked = LimitValuation(approximant, G, check=False)
    assert LimitValuation(approximant, G) is unchecked


def test_unchecked_composite_polynomial_is_refined_on_demand():
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = _polynomial_ring()
    x = R.gen()
    vK = QQ.valuation(2)
    G = (x**2 + 7) * (x**2 + 9)
    approximant = vK.mac_lane_approximants(
        G, require_incomparability=True)[0]

    unchecked = LimitValuation(approximant, G, check=False)
    current = LimitValuation(approximant, G)
    assert unchecked._G == G
    assert unchecked is not current
    assert unchecked >= current >= unchecked
    assert unchecked._G == current._G


def test_legacy_factory_key_is_refined_on_demand():
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = _polynomial_ring()
    x = R.gen()
    vK = QQ.valuation(2)
    G = (x**2 + 7) * (x**2 + 9)
    approximant = vK.mac_lane_approximants(
        G, require_incomparability=True)[0]
    current = LimitValuation(approximant, G)
    legacy = LimitValuation.get_object((0,), (approximant, 2*G), {})

    assert legacy._G == G
    assert legacy(approximant.phi()) == current(approximant.phi())
    assert legacy(current._G) == current(current._G)
    assert legacy._G == current._G


def test_legacy_keys_with_shared_support_compare_by_branch():
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = _polynomial_ring()
    x = R.gen()
    support = x**2 + 7
    left_polynomial = support * (x**2 + 9)
    right_polynomial = support * x
    approximants = QQ.valuation(2).mac_lane_approximants(
        left_polynomial, require_incomparability=True)
    approximants = [v for v in approximants
                    if not v.is_equivalence_unit(support)]

    limits = []
    for i, approximant in enumerate(approximants):
        left = LimitValuation.get_object(
            (1, i, 0), (approximant, left_polynomial), {})
        right = LimitValuation.get_object(
            (1, i, 1), (approximant, right_polynomial), {})
        assert left >= right >= left
        assert left._G == right._G == support
        limits.append(left)

    assert len(limits) == 2
    assert not (limits[0] >= limits[1])
    assert not (limits[1] >= limits[0])


def test_bounded_mac_lane_step_can_need_the_full_principal_part():
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.gauss_valuation import GaussValuation
    from sage.rings.valuation.inductive_valuation import (
        EquivalenceDecompositionTooSmall,
    )

    R = _polynomial_ring()
    x = R.gen()
    valuation = GaussValuation(R, QQ.valuation(2))
    G = x**2 - 2*x - 1
    options = {
        'assume_squarefree': True,
        'assume_equivalence_irreducible': True,
        'check': False,
        'report_degree_bounds_and_caches': True,
    }

    with pytest.raises(EquivalenceDecompositionTooSmall):
        valuation.mac_lane_step(G, principal_part_bound=1, **options)
    steps = valuation.mac_lane_step(
        G, principal_part_bound=None, **options)
    assert len(steps) == 1
    assert steps[0][0].mu() == QQ(1) / 2


def test_invalid_polynomial_is_rejected():
    from sage.rings.integer_ring import ZZ
    from sage.rings.padics.factory import Qp
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.gauss_valuation import GaussValuation
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = _polynomial_ring()
    x = R.gen()
    gauss = GaussValuation(R, QQ.valuation(2))
    with pytest.raises(ValueError, match="squarefree"):
        LimitValuation(gauss, x**2)

    S = PolynomialRing(ZZ, 'y')
    y = S.gen()
    gauss = GaussValuation(S, ZZ.valuation(2))
    assert LimitValuation(gauss, -y) is LimitValuation(gauss, y)
    with pytest.raises(ValueError, match="leading coefficient"):
        LimitValuation(gauss, 2*y**2 + y + 2)

    K = Qp(5)
    T = PolynomialRing(K, 'z')
    gauss = GaussValuation(T, K.valuation())
    with pytest.raises(NotImplementedError, match="inexact rings"):
        LimitValuation(gauss, T.gen())


@pytest.mark.parametrize('degree', [1, 2])
@pytest.mark.parametrize('operation', ['reduce', 'upper_bound'])
@pytest.mark.parametrize('denominator', [1, 2])
def test_fresh_limit_operations_on_support(degree, operation, denominator):
    from sage.rings.infinity import infinity
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.limit_valuation import LimitValuation

    # Separate factory keys keep other operations from refining the object
    # before the operation under test encounters its support.
    R = PolynomialRing(QQ, f'support_{degree}_{operation}_{denominator}')
    x = R.gen()
    G = x - 1 if degree == 1 else x**2 + 1
    approximant, = QQ.valuation(2).mac_lane_approximants(
        G, require_incomparability=True)
    limit = LimitValuation(approximant.augmentation(G, infinity), G)

    result = getattr(limit, operation)(G / denominator)
    assert result == (0 if operation == 'reduce' else infinity)


@pytest.mark.parametrize('branch', [0, 1])
@pytest.mark.parametrize('operation', ['reduce', 'upper_bound'])
def test_support_operations_refine_unchecked_composite(branch, operation):
    from sage.rings.infinity import infinity
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = PolynomialRing(QQ, f'composite_support_{branch}_{operation}')
    x = R.gen()
    support = x**2 + 7
    G = x * support
    approximants = QQ.valuation(2).mac_lane_approximants(
        G, require_incomparability=True)
    approximants = [a for a in approximants
                   if not a.is_equivalence_unit(support)]
    assert len(approximants) == 2
    limit = LimitValuation(approximants[branch], G, check=False)

    # This quadratic splits over QQ_2: no finite approximation can be
    # augmented by assigning infinity to the whole support polynomial.
    result = getattr(limit, operation)(support / 2**10)
    assert result == (0 if operation == 'reduce' else infinity)
    assert limit._G == support


def test_residue_ring_reuses_checks_until_approximation_or_support_changes(monkeypatch):
    from sage.rings.infinity import infinity
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = PolynomialRing(QQ, 'residue_cache')
    x = R.gen()
    support = x**2 + 7
    G = x * support
    approximant = next(a for a in QQ.valuation(2).mac_lane_approximants(
        G, require_incomparability=True)
        if not a.is_equivalence_unit(support))
    limit = LimitValuation(approximant, G, check=False)

    checks = []
    approximation_class = type(approximant)
    is_key = approximation_class.is_key

    def record_key_check(self, polynomial, *args, **kwds):
        checks.append((self, polynomial))
        return is_key(self, polynomial, *args, **kwds)

    monkeypatch.setattr(approximation_class, 'is_key', record_key_check)
    residue = limit.residue_ring()
    assert checks
    count = len(checks)
    for _ in range(3):
        assert limit.residue_ring() is residue
        assert limit.reduce(R.one()) == residue.one()
        assert limit.lift(residue.one()) == R.one()
    assert len(checks) == count

    # Discovering the support can change G without changing the approximation.
    assert limit(support) is infinity
    assert limit._G == support
    assert limit._approximation is approximant
    checks.clear()
    assert limit.residue_ring() is residue
    assert (approximant, support) in checks

    # A later Mac Lane step also invalidates the check for the old state.
    limit._improve_approximation()
    assert limit._approximation is not approximant
    checks.clear()
    assert limit.residue_ring() is residue
    assert (limit._approximation, support) in checks
    count = len(checks)
    assert limit.residue_ring() is residue
    assert len(checks) == count


def test_nontrivial_residue_field_round_trip_after_refinement():
    from sage.rings.infinity import infinity
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.gauss_valuation import GaussValuation
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = PolynomialRing(QQ, 'residue_extension')
    x = R.gen()
    G = (x**2 + x + 1)**2 + 2
    limit = LimitValuation(GaussValuation(R, QQ.valuation(2)), G)
    residue = limit.residue_ring()
    assert residue.cardinality() == 4
    assert limit._approximation.mu() is not infinity
    for element in residue:
        assert limit.reduce(limit.lift(element)) == element

    limit._improve_approximation()
    assert limit._approximation.mu() is infinity
    residue = limit.residue_ring()
    assert residue.cardinality() == 4
    for element in residue:
        assert limit.reduce(limit.lift(element)) == element


@pytest.mark.parametrize('final_approximation', [False, True])
@pytest.mark.parametrize('operation', ['residue_ring', 'reduce_support',
                                      'reduce_one', 'lift_one'])
def test_unsupported_integer_residue_operations(final_approximation, operation):
    from sage.rings.infinity import infinity
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.gauss_valuation import GaussValuation
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = PolynomialRing(ZZ, f'integer_residue_{final_approximation}_{operation}')
    x = R.gen()
    G = x**2 + 1
    approximation = GaussValuation(R, ZZ.valuation(2)).augmentation(x + 1, QQ(1)/2)
    if final_approximation:
        approximation = approximation.augmentation(G, infinity)
    limit = LimitValuation(approximation, G, check=not final_approximation)
    method, args = {
        'residue_ring': (limit.residue_ring, ()),
        'reduce_support': (limit.reduce, (G,)),
        'reduce_one': (limit.reduce, (R.one(),)),
        'lift_one': (limit.lift, (1,)),
    }[operation]
    with pytest.raises(NotImplementedError,
                       match=r'residue polynomial|polynomial rings over fields'):
        method(*args)


def test_integer_residual_polynomial_requires_normalizer():
    from sage.rings.infinity import infinity
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.gauss_valuation import GaussValuation

    R = PolynomialRing(ZZ, 'integer_residual_polynomial')
    x = R.gen()
    approximation = GaussValuation(R, ZZ.valuation(2)).augmentation(x + 1, QQ(1)/2)
    final_approximation = approximation.augmentation(x**2 + 1, infinity)
    with pytest.raises(NotImplementedError, match='residue polynomial'):
        final_approximation.psi()


@pytest.mark.parametrize('degree', [1, 2])
def test_supported_integer_residue_operations(degree):
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.valuation.gauss_valuation import GaussValuation
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = PolynomialRing(ZZ, f'supported_integer_residue_{degree}')
    x = R.gen()
    G = x - 1 if degree == 1 else x**2 + x + 1
    limit = LimitValuation(GaussValuation(R, ZZ.valuation(2)), G)
    residue = limit.residue_ring()
    assert residue.cardinality() == 2**degree
    assert limit.reduce(G) == residue.zero()
    for element in residue:
        lifted = limit.lift(element)
        assert lifted.parent() is R
        assert limit.reduce(lifted) == element


def test_rational_residue_operations_with_negative_normalizer():
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ
    from sage.rings.valuation.gauss_valuation import GaussValuation
    from sage.rings.valuation.limit_valuation import LimitValuation

    R = PolynomialRing(QQ, 'rational_negative_normalizer')
    x = R.gen()
    G = x**2 + 1
    approximation = GaussValuation(R, QQ.valuation(2)).augmentation(x + 1, QQ(1)/2)
    limit = LimitValuation(approximation, G)
    residue = limit.residue_ring()
    assert residue.cardinality() == 2
    assert limit.reduce(G) == residue.zero()
    for element in residue:
        assert limit.reduce(limit.lift(element)) == element


def test_integer_equivalence_unit_rejects_negative_value():
    from sage.rings.integer_ring import ZZ
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.valuation.gauss_valuation import GaussValuation

    R = PolynomialRing(ZZ, 'integer_equivalence_unit')
    gauss = GaussValuation(R, ZZ.valuation(2))
    with pytest.raises(ValueError, match='Semigroup'):
        gauss.equivalence_unit(-1)
