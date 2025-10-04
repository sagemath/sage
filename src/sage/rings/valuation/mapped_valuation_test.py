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


@pytest.fixture
def w0(w):
    return w[0]


@pytest.fixture
def w1(w):
    return w[1]


# Use strings for the fixture names here, and then later convert them
# to the actual fixture objects using request.getfixturevalue(). This
# is a workaround for being unable to pass fixtures directly as
# parameters:
#
#   https://github.com/pytest-dev/pytest/issues/349
#
extensions = ( "w0", "w1" )


@pytest.mark.parametrize("e", extensions)
def test_finite_extension_from_limit_valuation(e, request):
    r"""
    Run the ``TestSuite()`` for two examples given in the
    ``FiniteExtensionFromLimitValuation`` documentation.
    """
    from sage.misc.sage_unittest import TestSuite

    # Convert the string to a real fixture
    e = request.getfixturevalue(e)

    # fewer max_runs, these are kind of slow
    TestSuite(e).run(verbose=True,
                     raise_on_failure=True,
                     max_runs=512)
