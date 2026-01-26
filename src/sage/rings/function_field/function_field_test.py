import pytest

from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.function_field.constructor import FunctionField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ


@pytest.fixture
def F():
    return FunctionField(QQ, 'x')


@pytest.fixture
def J():
    return FunctionField(FiniteField(5), 'x')


@pytest.fixture
def K():
    return FunctionField(FiniteField(5**2, 'a'), 'x')


@pytest.fixture
def L(F):
    x = F.gen()
    Y = PolynomialRing(F, 'Y').gen()
    return F.extension(Y**2 + Y + x + 1 / x, 'y')


@pytest.fixture
def M(K, R, S):
    x = K.gen()
    y = R.gen()
    L = K.extension(y**3 - (x**3 + 2 * x * y + 1 / x))
    t = S.gen()
    return L.extension(t**2 - x * y)


@pytest.fixture
def N(K):
    return FunctionField(K, 'u')


@pytest.fixture
def O(L):
    return L.maximal_order()


@pytest.fixture
def R(K):
    return PolynomialRing(K, 'y')


@pytest.fixture
def S(K, R):
    x = K.gen()
    y = R.gen()
    L = K.extension(y**3 - (x**3 + 2 * x * y + 1 / x))
    return PolynomialRing(L, 't')


@pytest.fixture
def T(F):
    return PolynomialRing(F, 'Y')


# Use strings for the fixture names here, and then later convert them
# to the actual fixture objects using request.getfixturevalue(). This
# is a workaround for being unable to pass fixtures directly as
# parameters:
#
#   https://github.com/pytest-dev/pytest/issues/349
#
pairs = [("J", None),
         ("K", 16),
         ("L", 2),
         ("M", 1),
         ("N", 1),
         ("O", None),
         ("T", None),
         ("S", 8)]


@pytest.mark.parametrize("ff,max_runs", pairs)
def test_function_field_testsuite(ff, max_runs, request) -> None:
    r"""
    Run the TestSuite() on some function fields that are
    constructed in the documentation. They are slow, random, and not
    intended for end users. All of this makes them more appropriate to
    be run separately, here, than in the doctests.

    INPUT:

    The inputs are essentially all fixtures.

    - ``ff`` -- string; a function field fixture name
    - ``max_runs`` -- integer; the maximum number of times to
      repeat the test suite
    - ``request`` -- fixture; a pytest built-in

    """
    # The sage.misc.sage_unittest.TestSuite import is local to avoid
    # pytest warnings.
    from sage.misc.sage_unittest import TestSuite

    # Convert the fixture name (string) to an actual object using the
    # built-in "request" fixture.
    ff = request.getfixturevalue(ff)

    # Pass max_runs only if it's not None; otherwise use the default
    run_args = {"verbose": True, "raise_on_failure": True}
    if max_runs:
        run_args["max_runs"] = max_runs

    TestSuite(ff).run(**run_args)
