import pytest


def test_positive_integer_semigroup():
    r"""
    Run the ``TestSuite()`` for ``PositiveIntegerSemigroup``
    (this can take quite a long time).
    """
    from sage.misc.sage_unittest import TestSuite
    from sage.combinat.backtrack import PositiveIntegerSemigroup
    PP = PositiveIntegerSemigroup()

    # fewer max_runs since these are kind of slow
    TestSuite(PP).run(verbose=True,
                      raise_on_failure=True,
                      max_runs=256)
