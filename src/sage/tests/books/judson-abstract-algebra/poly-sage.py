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
## Section 17.6 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: R.<x> = Integers(8)[]; R
    Univariate Polynomial Ring in x over Ring of integers modulo 8

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S.<y> = ZZ[]; S
    Univariate Polynomial Ring in y over Integer Ring

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: T.<z> = QQ[]; T
    Univariate Polynomial Ring in z over Rational Field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: R.is_finite()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: R.is_integral_domain()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S.is_integral_domain()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: T.is_field()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: R.characteristic()
    8

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: T.characteristic()
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: y in S
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: x in S
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: q = (3/2) + (5/4)*z^2
    sage: q in T
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: 3 in S
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: r = 3
    sage: r.parent()
    Integer Ring

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: s = 3*y^0
    sage: s.parent()
    Univariate Polynomial Ring in y over Integer Ring

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p = 3 + 5*x + 2*x^2
    sage: p.parent()
    Univariate Polynomial Ring in x over Ring of integers modulo 8

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p(1)
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [p(t) for t in Integers(8)]
    [3, 2, 5, 4, 7, 6, 1, 0]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: q = 4*x^2+4*x
    sage: [q(t) for t in Integers(8)]
    [0, 0, 0, 0, 0, 0, 0, 0]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: M.<s, t> = QQ[]; M
    Multivariate Polynomial Ring in s, t over Rational Field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: R.<x> = QQ[]
    sage: p = 1/4*x^4 - x^3 + x^2 - x - 1/2
    sage: p.is_irreducible()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p.factor()
    (1/4) * (x^4 - 4*x^3 + 4*x^2 - 4*x - 2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: q = 2*x^5 + 5/2*x^4 + 3/4*x^3 - 25/24*x^2 - x - 1/2
    sage: q.is_irreducible()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: q.factor()
    (2) * (x^2 + 3/2*x + 3/4) * (x^3 - 1/4*x^2 - 1/3)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.<a> = FiniteField(5^2)
    sage: S.<y> = F[]
    sage: p = 2*y^5 + 2*y^4 + 4*y^3 + 2*y^2 + 3*y + 1
    sage: p.is_irreducible()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p.factor()
    (2) * (y^5 + y^4 + 2*y^3 + y^2 + 4*y + 3)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: q = 3*y^4+2*y^3-y+4; q.factor()
    (3) * (y^2 + (a + 4)*y + 2*a + 3) * (y^2 + 4*a*y + 3*a)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: r = y^4+2*y^3+3*y^2+4; r.factor()
    (y + 4) * (y^3 + 3*y^2 + y + 1)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: s = 3*y^4+2*y^3-y+3; s.factor()
    (3) * (y + 1) * (y + 3) * (y + 2*a + 4) * (y + 3*a + 1)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.modulus()
    x^2 + 4*x + 2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [t for t in F if r(t)==0]
    [1]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [t for t in F if s(t)==0]
    [2, 3*a + 1, 4, 2*a + 4]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: W.<w> = QQ[]
    sage: p = 16*w^5 - 9*w^4 +3*w^2 + 6*w -21
    sage: p.is_irreducible()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: W.<w> = QQ[]
    sage: r = -w^5 + 5*w^4 - 4*w^3 + 14*w^2 - 67*w + 17
    sage: s = 3*w^5 - 14*w^4 + 12*w^3 - 6*w^2 + w
    sage: S = W.ideal(r, s)
    sage: S
    Principal ideal (w^2 - 4*w + 1) of
    Univariate Polynomial Ring in w over Rational Field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (w^2)*r + (3*w-6)*s in S
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F = Integers(7)
    sage: R.<x> = F[]
    sage: p = x^5+ x + 4
    sage: p.is_irreducible()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: id = R.ideal(p)
    sage: Q = R.quotient(id); Q
    Univariate Quotient Polynomial Ring in xbar over
    Ring of integers modulo 7 with modulus x^5 + x + 4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.is_field()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.order() == 7^5
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.gen(0)
    xbar

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.<t> = R.quotient(id); Q
    Univariate Quotient Polynomial Ring in t over
    Ring of integers modulo 7 with modulus x^5 + x + 4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: t^5 + t + 4
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: t^5 == -(t+4)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: t^5
    6*t + 3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (3*t^3 + t + 5)*(t^2 + 4*t + 2)
    5*t^4 + 2*t^2 + 5*t + 5

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 3*t^4 - 6*t^3 + 3*t^2 + 5*t + 2
    sage: ainv = a^-1; ainv
    6*t^4 + 5*t^2 + 4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a*ainv
    1
"""
