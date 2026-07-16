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
def test_finite_extension_from_limit_valuation_w(w, idx, run_test_suite):
    r"""
    Run the ``TestSuite()`` for two examples given in the
    ``FiniteExtensionFromLimitValuation`` documentation.
    """
    # fewer max_runs, these are kind of slow
    run_test_suite(w[idx], verbose=True, max_runs=512)
