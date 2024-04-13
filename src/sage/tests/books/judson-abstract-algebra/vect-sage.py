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
## Section 20.6 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: V = QQ^4; V
    Vector space of dimension 4 over Rational Field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.<a> = FiniteField(3^4)
    sage: W = F^5; W
    Vector space of dimension 5 over Finite Field in a of size 3^4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: v = vector(QQ, [1, 1/2, 1/3, 1/4]); v
    (1, 1/2, 1/3, 1/4)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: v in V
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: w = vector(F, [1, a^2, a^4, a^6, a^8]); w
    (1, a^2, a^3 + 1, a^3 + a^2 + a + 1, a^2 + a + 2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: w in W
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: u = vector(QQ, [ 1, 2,  3, 4,   5,  6])
    sage: v = vector(QQ, [-1, 2, -4, 8, -16, 32])
    sage: 3*u - 2*v
    (5, 2, 17, -4, 47, -46)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: w = vector(F, [1, a^2, a^4, a^6,  a^8])
    sage: x = vector(F, [1,   a, 2*a,   a,    1])
    sage: y = vector(F, [1, a^3, a^6, a^9, a^12])
    sage: a^25*w + a^43*x + a^66*y
    (a^3 + a^2 + a + 2, a^2 + 2*a, 2*a^3 + a^2 + 2, 2*a^3 + a^2 + a,
     a^3 + 2*a^2 + a + 2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: u = vector(QQ, [1, -1, 3])
    sage: v = vector(QQ, [2, 1, -1])
    sage: w = vector(QQ, [3, 0, 2])
    sage: S = (QQ^3).subspace([u, v, w]); S
    Vector space of degree 3 and dimension 2 over Rational Field
    Basis matrix:
    [   1    0  2/3]
    [   0    1 -7/3]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: 3*u - 6*v + (1/2)*w in S
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: vector(QQ, [4, -1, -2]) in S
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S.basis()
    [
    (1, 0, 2/3),
    (0, 1, -7/3)
    ]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: u = vector(QQ, [1, -1,  3])
    sage: v = vector(QQ, [2,  1, -1])
    sage: w = vector(QQ, [3,  0,  2])
    sage: u + v == w
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S1 = (QQ^3).subspace([u, v, w])
    sage: S2 = (QQ^3).subspace([u-v, v-w, w-u])
    sage: S1 == S2
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: u = vector(QQ, [1, -1,  3,  4])
    sage: v = vector(QQ, [2,  1, -1, -2])
    sage: S = (QQ^4).subspace([u, v, 2*u + 3*v, -u + 2*v])
    sage: S.dimension()
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.<a> = FiniteField(3^4)
    sage: u = vector(F, [a^i for i in range(0,  7, 1)])
    sage: v = vector(F, [a^i for i in range(0, 14, 2)])
    sage: w = vector(F, [a^i for i in range(0, 21, 3)])
    sage: S = (F^7).subspace([u, v, w])
    sage: S.dimension()
    3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S = (F^7).subspace([u, v, a^3*u + a^11*v])
    sage: S.dimension()
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.<a> = FiniteField(7^6)
    sage: u = 2*a^5 + 6*a^4 + 2*a^3 + 3*a^2 + 2*a + 3
    sage: v = 4*a^5 + 4*a^4 + 4*a^3 + 6*a^2 + 5*a + 6
    sage: u + v
    6*a^5 + 3*a^4 + 6*a^3 + 2*a^2 + 2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: 4*u
    a^5 + 3*a^4 + a^3 + 5*a^2 + a + 5

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: 2*u + 5*v
    3*a^5 + 4*a^4 + 3*a^3 + a^2 + a + 1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: V = F.vector_space(map=False); V
    Vector space of dimension 6 over Finite Field of size 7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: R = V.base_ring(); R
    Finite Field of size 7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: R == FiniteField(7)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: V.dimension()
    6

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: x = V(u); x
    (3, 2, 3, 2, 6, 2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: y = V(v); y
    (6, 5, 6, 4, 4, 4)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: V(u + v) == V(u) + V(v)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: two = R(2)
    sage: V(two*u) == two*V(u)
    True

"""
