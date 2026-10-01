import pytest


@pytest.fixture
def w():
    from sage.rings.function_field.constructor import FunctionField
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.integer import Integer
    from sage.rings.rational_field import QQ

    K = FunctionField(QQ, 'x')
    x = K.gen()
    R = PolynomialRing(K, 'y')
    y = R.gen()
    L = K.extension(y**2 - x)
    v = K.valuation(Integer(1))

    return v.extensions(L)


@pytest.mark.long
@pytest.mark.parametrize("idx", [0,1])
def test_finite_extension_from_limit_valuation_w(w, idx):
    r"""
    Run the ``TestSuite()`` for two examples given in the
    ``FiniteExtensionFromLimitValuation`` documentation.
    """
    from sage.misc.sage_unittest import TestSuite

    # fewer max_runs, these are kind of slow
    TestSuite(w[idx]).run(verbose=True,
                          raise_on_failure=True,
                          max_runs=512)


def test_finite_extension_preserves_an_infinite_approximant():
    from sage.rings.function_field.constructor import FunctionField
    from sage.rings.integer import Integer
    from sage.rings.infinity import infinity
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.rational_field import QQ

    K = FunctionField(QQ, 'x')
    x = K.gen()
    R = PolynomialRing(K, 'y')
    y = R.gen()
    L = K.extension(y**2 - x)
    valuation = K.valuation(Integer(2)).extension(L)
    limit = valuation._base_valuation

    assert limit._initial_approximation.mu() is infinity
    assert limit._initial_approximation in valuation._approximants
    assert repr(valuation) == '(x - 2)-adic valuation'
