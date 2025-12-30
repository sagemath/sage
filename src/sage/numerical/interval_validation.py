from sage.rings.real_mpfr import RR
from sage.rings.real_interval_field import RealIntervalField

RIF = RealIntervalField(53)

def validate_numeric_result(f, x):
    """
    Validate a numerical computation using interval arithmetic.

    INPUT:
        - f: a callable accepting a real argument
        - x: a real number

    OUTPUT:
        - (value, interval, reliable)

    EXAMPLES::

        sage: from sage.numerical.interval_validation import validate_numeric_result
        sage: f = lambda t: t^2 - 2
        sage: val, interval, ok = validate_numeric_result(f, 1.41421356)
        sage: ok
        True
    """
    val = f(RR(x))
    ival = f(RIF(x))
    return val, ival, val in ival

