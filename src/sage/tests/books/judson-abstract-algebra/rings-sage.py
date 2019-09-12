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
## Section 16.9 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F = QuadraticField(7)
    sage: F
    Number Field in a with defining polynomial x^2 - 7 with a = 2.645751311064591?

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: root = F.gen(0)
    sage: root^2
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: root
    a

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (2*root)^3
    56*a

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.<b> = QuadraticField(7)
    sage: F
    Number Field in b with defining polynomial x^2 - 7 with b = 2.645751311064591?

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b^2
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (2*b)^3
    56*b

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C.<t> = CyclotomicField(8)
    sage: C.random_element()   # random
    -2/11*t^2 + t - 1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z7 = Integers(7)
    sage: Z9 = Integers(9)
    sage: Q = QuadraticField(-11)
    sage: F.<a> = FiniteField(3^2)
    sage: P.<x> = Z7[]
    sage: S.<f,g,h> = QuaternionAlgebra(-7, 3)
    

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: QQ.is_exact()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: RR.is_exact()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z7.is_finite()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z7.is_finite()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z7.is_integral_domain()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z9.is_integral_domain()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z9.is_field()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.is_field()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.is_field()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.is_commutative()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S.is_commutative()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z7.characteristic()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z9.characteristic()
    9

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.characteristic()
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F.characteristic()
    3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: P.characteristic()
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S.characteristic()
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b = Z9.zero(); b
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b.parent()
    Ring of integers modulo 9

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c = Q.zero(); c
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c.parent()
    Number Field in a with defining polynomial x^2 + 11 with a = 3.316624790355400?*I

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b == c
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: d = Z9.one(); d
    1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: d.parent()
    Ring of integers modulo 9

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: e = Q.one(); e
    1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: e.parent()
    Number Field in a with defining polynomial x^2 + 11 with a = 3.316624790355400?*I

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: d == e
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: QQ.is_subring(Q)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: QQ.is_subring(S)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: QQ.is_subring(F)
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: three = Z9(3)
    sage: three.is_unit()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: three*three
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: four = Z9(4)
    sage: four.is_unit()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: g = four^-1; g
    7

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: four*g
    1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: I1 = ZZ.ideal(4)
    sage: I2 = 4*ZZ
    sage: I3 = (-4)*ZZ
    sage: I1 == I2
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: I2 == I3
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q = ZZ.quotient(I1); Q
    Ring of integers modulo 4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q == Integers(4)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z7 = Integers(7)
    sage: P.<y> = Z7[]
    sage: M = P.ideal(y^2+4)
    sage: Q = P.quotient(M)
    sage: Q
    Univariate Quotient Polynomial Ring in ybar over
    Ring of integers modulo 7 with modulus y^2 + 4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.random_element()   # random
    2*ybar + 6

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.order()
    49

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.is_field()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.<t> = P.quotient(M); Q
    Univariate Quotient Polynomial Ring in t over
    Ring of integers modulo 7 with modulus y^2 + 4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.random_element()   # random
    4*t + 6

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z7 = Integers(7)
    sage: P.<y> = Z7[]
    sage: M = P.ideal(y^2+3)
    sage: Q.<t> = P.quotient(M)
    sage: Q
    Univariate Quotient Polynomial Ring in t over
    Ring of integers modulo 7 with modulus y^2 + 3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.random_element()   # random
    3*t + 1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.order()
    49

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.is_field()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z7 = Integers(7)
    sage: P.<y> = Z7[]
    sage: M = P.ideal(y^2+4)
    sage: N = P.ideal(y^2+3)
    sage: M.is_maximal()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: N.is_maximal()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: M.is_prime()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: N.is_prime()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H = Hom(ZZ, QQ)
    sage: phi = H([1])
    sage: phi
    Ring morphism:
      From: Integer Ring
      To:   Rational Field
      Defn: 1 |--> 1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: phi.parent()
    Set of Homomorphisms from Integer Ring to Rational Field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 3; a
    3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a.parent()
    Integer Ring

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b = phi(3); b
    3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b.parent()
    Rational Field

"""
