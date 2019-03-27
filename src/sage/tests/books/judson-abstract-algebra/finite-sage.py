##      -*-   coding: utf-8   -*-     ##
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
## Section 22.6 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.<a> = GF(7^15); F
    Finite Field in a of size 7^15

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.polynomial()
    a^15 + 5*a^6 + 6*a^5 + 6*a^4 + 4*a^3 + a^2 + 2*a + 4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a^15 + 5*a^6 + 6*a^5 + 6*a^4 + 4*a^3 + a^2 + 2*a + 4
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: conway_polynomial(7, 15)
    x^15 + 5*x^6 + 6*x^5 + 6*x^4 + 4*x^3 + x^2 + 2*x + 4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: y = polygen(Integers(7), 'y')
    sage: P = y.parent()
    sage: p = P([4, 5, 2, 6, 3, 3, 6, 2, 1, 1, 2, 5, 6, 3, 5, 1]); p
    y^15 + 5*y^14 + 3*y^13 + 6*y^12 + 5*y^11 + 2*y^10 + y^9 +
    y^8 + 2*y^7 + 6*y^6 + 3*y^5 + 3*y^4 + 6*y^3 + 2*y^2 + 5*y + 4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p.is_irreducible()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: T.<b> = GF(7^15, modulus=p); T
    Finite Field in b of size 7^15

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.<a> = GF(5^4)
    sage: a^458
    3*a^3 + 2*a^2 + a + 3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (3*a^3 + 2*a^2 + a + 3).log(a)
    458

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: exponent = (3*a^3 + 2*a^2 + a + 3).log(2*a^3 + 4*a^2 + 4*a)
    sage: exponent
    211

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (2*a^3 + 4*a^2 + 4*a)^exponent == 3*a^3 + 2*a^2 + a + 3
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (3*a^3 + 2*a^2 + a + 3).log(a^2 + 4*a + 4)
    Traceback (most recent call last):
    ...
    ValueError: No discrete log of 3*a^3 + 2*a^2 + a + 3 found
    to base a^2 + 4*a + 4

"""
