import pytest


def test_finite_extension_from_limit_valuation():
    r"""
    Run the ``TestSuite()`` for two examples given in the
    ``FiniteExtensionFromLimitValuation`` documentation.
    """
    from sage.misc.sage_unittest import TestSuite
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
    w = v.extensions(L)

    for ext in w:
        # fewer max_runs, these are kind of slow
        TestSuite(ext).run(verbose=True,
                           raise_on_failure=True,
                           max_runs=512)
