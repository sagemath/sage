# sage_setup: distribution = sagemath-repl
##          Sage Doctest File         ##
#**************************************#
#*    Generated from PreTeXt source   *#
#*    on 2017-08-24T11:43:34-07:00    *#
#*                                    *#
#*   http://mathbook.pugetsound.edu   *#
#*                                    *#
#**************************************#
##
"""
Please contact Rob Beezer (beezer@ups.edu) with
any test failures here that need to be changed
as a result of changes accepted into Sage.  You
may edit/change this file in any sensible way, so
that development work may procede.  Your changes
may later be replaced by the authors of "Abstract
Algebra: Theory and Applications" when the text is
updated, and a replacement of this file is proposed
for review.
"""
##
## To execute doctests in these files, run
##   $ $SAGE_ROOT/sage -t <directory-of-these-files>
## or
##   $ $SAGE_ROOT/sage -t <a-single-file>
##
## Replace -t by "-tp n" for parallel testing,
##   "-tp 0" will use a sensible number of threads
##
## See: http://www.sagemath.org/doc/developer/doctesting.html
##   or run  $ $SAGE_ROOT/sage --advanced  for brief help
##
## Generated at 2017-08-24T11:43:34-07:00
## From "Abstract Algebra"
## At commit 26d3cac0b4047f4b8d6f737542be455606e2c4b4
##
## Section 2.6 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: r = 14 % 3
    sage: r
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: q = (14 - r)/3
    sage: q
    4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 14
    sage: b = 3
    sage: a.quo_rem(b)
    (4, 2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (20 % 5) == 0
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (17 % 4) == 0
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c = 5
    sage: c.divides(20)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: d = 4
    sage: d.divides(17)
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: gcd(2776, 2452)
    4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 31049
    sage: b = 2105
    sage: gcd(a, b) == 1
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 3563
    sage: b = 2947
    sage: gcd(a, b) == 1
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: xgcd(633,331)
    (1, -137, 262)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 633
    sage: b = 331
    sage: extended = xgcd(a, b)
    sage: g = extended[0]
    sage: r = extended[1]
    sage: s = extended[2]
    sage: g == r*a + s*b
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 117371
    sage: a.is_prime()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b = 14547073
    sage: b.is_prime()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b == 1597 * 9109
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = random_prime(10^21, proof=True)
    sage: a   # random
    424729101793542195193

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a.is_prime()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: prime_range(500, 550)
    [503, 509, 521, 523, 541, 547]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 2600
    sage: a.factor()
    2^3 * 5^2 * 13

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 2600
    sage: factored = a.factor()
    sage: first_term = factored[0]
    sage: first_term
    (2, 3)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: second_term = factored[1]
    sage: second_term
    (5, 2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: third_term = factored[2]
    sage: third_term
    (13, 1)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: first_prime = first_term[0]
    sage: first_prime
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: first_exponent = first_term[1]
    sage: first_exponent
    3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: list(factored)
    [(2, 3), (5, 2), (13, 1)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: len(factored)
    3
"""
