import pytest
import warnings


def test_padic_lattice_element():
    r"""
    Run the ``TestSuite()`` for some examples that previously
    lived in the TESTS:: block of the padic_lattice_element module.
    """
    from sage.misc.sage_unittest import TestSuite
    from sage.rings.padics.factory import ZpLC, ZpLF, QpLC, QpLF

    with warnings.catch_warnings(category=FutureWarning):
        warnings.filterwarnings("ignore", category=FutureWarning)
        # These all raise FutureWarnings
        R1 = ZpLC(2)
        R2 = ZpLF(2)
        R3 = QpLC(2)
        R4 = QpLF(2)

    for R in (R1, R2, R3, R4):
        # Only do a few runs, _test_matrix_smith() in particular is
        # sloooooow.
        TestSuite(R).run(verbose=True,
                         raise_on_failure=True,
                         skip="_test_teichmuller",
                         max_runs=8)
