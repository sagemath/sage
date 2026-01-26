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
## Section 21.6 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: M.<a> = QQ[sqrt(2)+sqrt(3)]; M
    Number Field in a with defining polynomial x^4 - 10*x^2 + 1 with a = 3.146264369941973?

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.<y> = QQ[]
    sage: p = y^3 - 1/4*y^2 - 1/16*y + 1/4
    sage: p.is_irreducible()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: N.<b> = NumberField(p, 'b'); N
    Number Field in b with
    defining polynomial y^3 - 1/4*y^2 - 1/16*y + 1/4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: z = polygen(QQ, 'z')
    sage: q = z^3 - 1/4*z^2 - 1/16*z + 1/4
    sage: q.parent()
    Univariate Polynomial Ring in z over Rational Field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: P.<c> = NumberField(q, 'c'); P
    Number Field in c with
    defining polynomial z^3 - 1/4*z^2 - 1/16*z + 1/4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: M.polynomial()
    x^4 - 10*x^2 + 1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: N.polynomial()
    y^3 - 1/4*y^2 - 1/16*y + 1/4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: element = -b^2 + 1/3*b + 4
    sage: element.parent()
    Number Field in b with
    defining polynomial y^3 - 1/4*y^2 - 1/16*y + 1/4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: r = element.minpoly('t'); r
    t^3 - 571/48*t^2 + 108389/2304*t - 13345/216

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: r.parent()
    Univariate Polynomial Ring in t over Rational Field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: r.subs(t=element)
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A.<a,b> = QQ[sqrt(2), sqrt(3)]
    sage: A
    Number Field in sqrt2 with defining polynomial x^2 - 2 over
    its base field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: B = A.base_field(); B
    Number Field in sqrt3 with defining polynomial x^2 - 3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A.is_relative()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: B.is_relative()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C.<c> = A.absolute_field()
    sage: C
    Number Field in c with defining polynomial x^4 - 10*x^2 + 1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromC, toC = C.structure()
    sage: fromC(c)
    sqrt2 - sqrt3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: toC(a)
    1/2*c^3 - 9/2*c

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: toC(b)
    1/2*c^3 - 11/2*c

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: B.degree()
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A.absolute_degree()
    4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A.relative_degree()
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: x = polygen(QQ, 'x')
    sage: p = x^4 + x^2 - 1
    sage: p.parent()
    Univariate Polynomial Ring in x over Rational Field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p.is_irreducible()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: M.<a> = NumberField(p, 'a')
    sage: y = polygen(M, 'y')
    sage: p = p.subs(x = y)
    sage: p
    y^4 + y^2 - 1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p.parent()
    Univariate Polynomial Ring in y over Number Field in a with
    defining polynomial x^4 + x^2 - 1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p.factor()
    (y - a) * (y + a) * (y^2 + a^2 + 1)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a^2 + 1 in QQ
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: q = y^2 + a^2 + 1
    sage: N.<b> = NumberField(q, 'b')
    sage: R.<z> = N[]
    sage: s = R(p)
    sage: s
    z^4 + z^2 - 1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: s.parent()
    Univariate Polynomial Ring in z over Number Field in b with
    defining polynomial y^2 + a^2 + 1 over its base field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: s.factor()
    (z + b) * (z + a) * (z - a) * (z - b)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a in N, b in N
    (True, True)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: P.<c> = N.absolute_field()
    sage: w = polygen(P, 'w')
    sage: p = w^4 + w^2- 1
    sage: p.factor()
    (w - 7/18966*c^7 + 110/9483*c^5 +  923/9483*c^3 +  3001/6322*c) *
    (w - 7/37932*c^7 +  55/9483*c^5 + 923/18966*c^3 - 3321/12644*c) *
    (w + 7/37932*c^7 -  55/9483*c^5 - 923/18966*c^3 + 3321/12644*c) *
    (w + 7/18966*c^7 - 110/9483*c^5 -  923/9483*c^3 -  3001/6322*c)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromP, toP = P.structure()
    sage: fromP(7/18966*c^7 - 110/9483*c^5 - 923/9483*c^3 - 3001/6322*c)
    -b

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: M.degree()
    4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: N.relative_degree()
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: P.degree()
    8

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: M.degree()*N.relative_degree() == P.degree()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: x = polygen(QQ, 'x')
    sage: p = x^4 + x^2 - 1
    sage: r = p.roots(ring=QQbar); r
    [(-0.7861513777574233?,  1), (0.7861513777574233?,  1),
     (-1.272019649514069?*I, 1), (1.272019649514069?*I, 1)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: r1 = r[0][0]; r1
    -0.7861513777574233?

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: r1.as_number_field_element()
    (Number Field in a with defining polynomial y^4 - y^2 - 1,
     a^3 - a,
     Ring morphism:
       From: Number Field in a with defining polynomial y^4 - y^2 - 1
       To:   Algebraic Real Field
       Defn: a |--> -1.272019649514069?)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: r1^4 + r1^2 - 1
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: N, rexact, homomorphism = r1.as_number_field_element()
    sage: (rexact)^4 + rexact^2 - 1
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: log(QQ[cos(pi/9)].degree(), 2) in ZZ
    False
"""
