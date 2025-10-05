import pytest
from sage.rings.padics.factory import ZpLC, ZpLF, QpLC, QpLF

# ZpLC, ZpLF, QpLC, and QpLF all raise FutureWarnings
from warnings import catch_warnings, filterwarnings
filterwarnings("ignore", category=FutureWarning)


@pytest.fixture
def R1():
    return ZpLC(2)


@pytest.fixture
def R2():
    return ZpLF(2)


@pytest.fixture
def R3():
    return QpLC(2)


@pytest.fixture
def R4():
    return QpLF(2)


# Use strings for the fixture names here, and then later convert them
# to the actual fixture objects using request.getfixturevalue(). This
# is a workaround for being unable to pass fixtures directly as
# parameters:
#
#   https://github.com/pytest-dev/pytest/issues/349
#
elements = ( "R1", "R2", "R3", "R4" )


@pytest.mark.parametrize("e", elements)
def test_padic_lattice_element(e, request):
    r"""
    Run the ``TestSuite()`` for some examples that previously
    lived in the TESTS:: block of the padic_lattice_element module.
    """
    from sage.misc.sage_unittest import TestSuite

    # Convert the string to a real fixture
    e = request.getfixturevalue(e)

    # Only do a few runs, _test_matrix_smith() in particular is slow.
    TestSuite(e).run(verbose=True,
                     raise_on_failure=True,
                     skip="_test_teichmuller",
                     max_runs=8)
