def test_positive_integer_semigroup(run_test_suite):
    r"""
    Run the ``TestSuite()`` for ``PositiveIntegerSemigroup``
    (this can take quite a long time).
    """
    from sage.combinat.backtrack import PositiveIntegerSemigroup
    PP = PositiveIntegerSemigroup()

    # fewer max_runs since these are kind of slow
    run_test_suite(PP, verbose=True, max_runs=256)
