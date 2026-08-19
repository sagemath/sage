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
    assert a >= b and b >= a
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
    assert unchecked >= current and current >= unchecked
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
        assert left >= right and right >= left
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
