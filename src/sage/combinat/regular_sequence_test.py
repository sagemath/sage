import pytest

@pytest.mark.longlong
def test_regular_sequence_ring_testsuite():
    from itertools import islice
    from sage.combinat.regular_sequence import RegularSequenceRing
    from sage.misc.sage_unittest import TestSuite
    from sage.rings.integer_ring import ZZ

    Seq2 = RegularSequenceRing(2, ZZ)
    elts = tuple(islice(Seq2.some_elements(), 4))
    TestSuite(Seq2).run(elements=elts,
                        verbose=True,
                        raise_on_failure=True,
                        max_runs=16)
