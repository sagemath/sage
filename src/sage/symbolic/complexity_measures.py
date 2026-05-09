"""
Complexity Measures

Some measures of symbolic expression complexity. Each complexity
measure is expected to take a symbolic expression as an argument, and
return a number.
"""


def string_length(expr):
    """
    Return the length of ``expr`` after converting it to a string.

    INPUT:

    - ``expr`` -- the expression whose complexity we want to measure

    OUTPUT: a real number representing the complexity of ``expr``

    RATIONALE:

    If the expression is longer on-screen, then a human would probably
    consider it more complex.

    EXAMPLES:

    This expression has three characters, ``x``, ``^``, and ``2``::

        sage: from sage.symbolic.complexity_measures import string_length
        sage: f = x^2
        sage: string_length(f)
        3
    """
    return len(str(expr))
