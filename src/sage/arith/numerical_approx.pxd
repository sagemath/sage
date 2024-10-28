# sage_setup: distribution = sagemath-objects
# sage.doctest: needs sage.rings.real_mpfr

cpdef inline long digits_to_bits(d) except -1:
    """
    EXAMPLES::

        sage: from sage.arith.numerical_approx import digits_to_bits
        sage: digits_to_bits(None)
        53
        sage: digits_to_bits(15)
        54
        sage: digits_to_bits(-1)
        Traceback (most recent call last):
        ...
        ValueError: number of digits must be positive

    TESTS::

        sage: digits_to_bits("10")
        Traceback (most recent call last):
        ...
        TypeError: must be real number, not str
    """
    if d is None:
        return 53
    cdef double x = d
    if x <= 0:
        raise ValueError("number of digits must be positive")

    # Add 1 because of the way how we display real numbers by default
    x += 1
    cdef double LOG_TEN_TWO = 3.321928094887362
    # Add 1 to round up
    return <long>(x * LOG_TEN_TWO) + 1
