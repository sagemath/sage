from sage.tests.memcheck.verify_no_leak import verify_no_leak


def test_sqrt_sqrt_2() -> None:
    from sage.misc.functional import sqrt
    T2 = sqrt(2)

    def sqrt_T2() -> None:
        sqrt(T2)

    verify_no_leak(sqrt_T2)
