import pytest
from sage.rings.padics.factory import ZpLC, ZpLF, QpLC, QpLF


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
elements = ("R1", "R2", "R3", "R4")


# ZpLC, ZpLF, QpLC, and QpLF all raise FutureWarnings; ignore them just for
# these tests (rather than mutating the warning filters process-wide).
@pytest.mark.filterwarnings("ignore::FutureWarning")
@pytest.mark.long
@pytest.mark.parametrize("e", elements)
def test_padic_lattice_element(e, request, run_test_suite):
    r"""
    Run the ``TestSuite()`` for some examples that previously
    lived in the TESTS:: block of the padic_lattice_element module.
    """
    # Convert the string to a real fixture
    e = request.getfixturevalue(e)

    # Only do a few runs, _test_matrix_smith() in particular is slow.
    run_test_suite(e, verbose=True, skip="_test_teichmuller", max_runs=8)
