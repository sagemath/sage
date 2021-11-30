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
## Section 23.6 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: x = polygen(QQ, 'x')
    sage: N.<a> = NumberField(x^4 - 2); N
    Number Field in a with defining polynomial x^4 - 2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: L.<b> = N.galois_closure(); L
    Number Field in b with defining polynomial x^8 + 28*x^4 + 2500

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: L.degree()
    8

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: y = polygen(L, 'y')
    sage: (y^4 - 2).factor()
    (y - 1/120*b^5 -  19/60*b) *
    (y - 1/240*b^5 + 41/120*b) *
    (y + 1/240*b^5 - 41/120*b) *
    (y + 1/120*b^5 +  19/60*b)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = L.galois_group(); G
    Galois group 8T4 ([4]2) with order 8 of x^8 + 28*x^4 + 2500

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G.is_abelian()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G.order()
    8

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G.list()
    [(), (1,2,8,7)(3,4,6,5),
    (1,3)(2,5)(4,7)(6,8), (1,4)(2,3)(5,8)(6,7),
    (1,5)(2,6)(3,7)(4,8), (1,6)(2,4)(3,8)(5,7),
    (1,7,8,2)(3,5,6,4), (1,8)(2,7)(3,6)(4,5)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G.is_isomorphic(DihedralGroup(4))
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: x = polygen(QQ, 'x')
    sage: p = x^4 - 2
    sage: N.<a> = NumberField(p); N
    Number Field in a with defining polynomial x^4 - 2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: y = polygen(N, 'y')
    sage: p = p.subs(x=y)
    sage: p.factor()
    (y - a) * (y + a) * (y^2 + a^2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: M.<b> = NumberField(y^2 + a^2); M
    Number Field in b with defining polynomial y^2 + a^2 over
    its base field

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: z = polygen(M, 'z')
    sage: (z^4 - 2).factor()
    (z - b) * (z - a) * (z + a) * (z + b)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: L.<c> = M.absolute_field(); L
    Number Field in c with defining polynomial x^8 + 28*x^4 + 2500

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL, toL = L.structure()

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: roots = p.roots(ring=L, multiplicities=False); roots
    [1/120*c^5 +  19/60*c,
     1/240*c^5 - 41/120*c,
    -1/240*c^5 + 41/120*c,
    -1/120*c^5 -  19/60*c]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [fromL(r) for r in roots]
    [b, a, -a, -b]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = End(L); G
    Automorphism group of Number Field in c with
    defining polynomial x^8 + 28*x^4 + 2500

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [tau(1) for tau in G]
    [1, 1, 1, 1, 1, 1, 1, 1]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Sequence([[fromL(tau(r)) for r in roots] for tau in G], cr=True)
    [
    [b, a, -a, -b],
    [-b, -a, a, b],
    [a, -b, b, -a],
    [b, -a, a, -b],
    [-a, -b, b, a],
    [a, b, -b, -a],
    [-b, a, -a, b],
    [-a, b, -b, a]
    ]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S4 = SymmetricGroup(4)
    sage: elements = [S4([1, 2, 3, 4]),
    ....:             S4([4, 3, 2, 1]),
    ....:             S4([2, 4, 1, 3]),
    ....:             S4([1, 3, 2, 4]),
    ....:             S4([3, 4, 1, 2]),
    ....:             S4([2, 1, 4, 3]),
    ....:             S4([4, 2, 3, 1]),
    ....:             S4([3, 1, 4, 2])]
    sage: elements
    [(), (1,4)(2,3), (1,2,4,3), (2,3), (1,3)(2,4),
     (1,2)(3,4), (1,4), (1,3,4,2)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: P = S4.subgroup(elements)
    sage: P.is_isomorphic(DihedralGroup(4))
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: basis = L.power_basis(); basis
    [1, c, c^2, c^3, c^4, c^5, c^6, c^7]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: tau = G[3]
    sage: z = 4 + 5*c+ 6*c^3-7*c^6
    sage: tz = tau(4 + 5*c+ 6*c^3-7*c^6); tz
    11/250*c^7 - 98/25*c^6 + 1/12*c^5 + 779/125*c^3 +
    6006/25*c^2 - 11/6*c + 4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: tz.vector()
    (4, -11/6, 6006/25, 779/125, 0, 1/12, -98/25, 11/250)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: tau_matrix = column_matrix([tau(be).vector() for be in basis])
    sage: tau_matrix
    [  1       0       0        0  -28       0        0          0]
    [  0  -11/30       0        0    0  779/15        0          0]
    [  0       0  -14/25        0    0       0  -858/25          0]
    [  0       0       0  779/750    0       0        0  -4031/375]
    [  0       0       0        0   -1       0        0          0]
    [  0    1/60       0        0    0   11/30        0          0]
    [  0       0   -1/50        0    0       0    14/25          0]
    [  0       0       0  11/1500    0       0        0   -779/750]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: tau_matrix*z.vector()
    (4, -11/6, 6006/25, 779/125, 0, 1/12, -98/25, 11/250)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: tau_matrix*(z.vector()) == (tau(z)).vector()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: K = (tau_matrix-identity_matrix(8)).right_kernel(); K
    Vector space of degree 8 and dimension 4 over Rational Field
    Basis matrix:
    [    1     0     0     0     0     0     0     0]
    [    0     1     0     0     0  1/38     0     0]
    [    0     0     1     0     0     0 -1/22     0]
    [    0     0     0     1     0     0     0 1/278]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(1)
    1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(c + (1/38)*c^5)
    60/19*b

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(c^2 - (1/22)*c^6)
    150/11*a^2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(c^3 + (1/278)*c^7)
    1500/139*a^2*b

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a^2 + b^2
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sg = P.subgroups()
    sage: [H.gens() for H in sg]
    [[()],
     [(1,4)(2,3)],
     [(2,3)],
     [(1,4)],
     [(1,2)(3,4)],
     [(1,3)(2,4)],
     [(2,3), (1,4)(2,3)],
     [(1,2,4,3), (1,4)(2,3)],
     [(1,2)(3,4), (1,4)(2,3)],
     [(2,3), (1,2,4,3), (1,4)(2,3)]]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [H.order() for H in sg]
    [1, 2, 2, 2, 2, 2, 4, 4, 4, 8]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: tau = G[4]
    sage: tau_matrix = column_matrix([tau(be).vector() for be in basis])
    sage: (tau_matrix-identity_matrix(8)).right_kernel()
    Vector space of degree 8 and dimension 4 over Rational Field
    Basis matrix:
    [     1      0      0      0      0      0      0      0]
    [     0      1      0      0      0  1/158      0      0]
    [     0      0      1      0      0      0   1/78      0]
    [     0      0      0      1      0      0      0 13/614]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(tau(1))
    1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(tau(c+(1/158)*c^5))
    120/79*b - 120/79*a

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(tau(c^2+(1/78)*c^6))
    -200/39*a*b

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(tau(c^3+(13/614)*c^7))
    3000/307*a^2*b + 3000/307*a^3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (b-a)^2
    -2*a*b

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: (b-a)^3
    2*a^2*b + 2*a^3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: subinfo = L.subfield((79/120)*(c+(1/158)*c^5)); subinfo
    (Number Field in c0 with defining polynomial x^4 + 8 with c0 = 1/240*c^5 + 79/120*c,
     Ring morphism:
       From: Number Field in c0 with defining polynomial x^4 + 8 with c0 = 1/240*c^5 + 79/120*c
       To:   Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
       Defn: c0 |--> 1/240*c^5 + 79/120*c)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: V = QQ^8
    sage: for tau in [G[0], G[1], G[3], G[6]]:
    ....:   tau_matrix = column_matrix([tau(be).vector() for be in basis])
    ....:   K = (tau_matrix-identity_matrix(8)).right_kernel()
    ....:   V = V.intersection(K)
    sage: V
    Vector space of degree 8 and dimension 2 over Rational Field
    Basis matrix:
    [    1     0     0     0     0     0     0     0]
    [    0     0     1     0     0     0 -1/22     0]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(tau(c^2 - (1/22)*c^6))
    150/11*a^2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F, mapping = L.subfield((11/150)*(c^2 - (1/22)*c^6))
    sage: F
    Number Field in c0 with defining polynomial x^2 - 2 with c0 = -1/300*c^6 + 11/150*c^2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: V = QQ^8
    sage: for tau in [G[0], G[1], G[2], G[7]]:
    ....:   tau_matrix = column_matrix([tau(be).vector() for be in basis])
    ....:   K = (tau_matrix-identity_matrix(8)).right_kernel()
    ....:   V = V.intersection(K)
    sage: V
    Vector space of degree 8 and dimension 2 over Rational Field
    Basis matrix:
    [1 0 0 0 0 0 0 0]
    [0 0 0 0 1 0 0 0]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(tau(c^4))
    -24*a^3*b - 14

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: F, mapping = L.subfield((c^4+14)/-48)
    sage: F
    Number Field in c0 with defining polynomial x^2 + 1 with c0 = -1/48*c^4 - 7/24

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: L.subfields()
    [
    (Number Field in c0 with defining polynomial x,
     Ring morphism:
       From: Number Field in c0 with defining polynomial x
       To:   Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
       Defn: 0 |--> 0,
     None),
    (Number Field in c1 with defining polynomial x^2 + 112*x + 40000,
     Ring morphism:
       From: Number Field in c1 with defining polynomial x^2 + 112*x + 40000
       To:   Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
       Defn: c1 |--> 4*c^4,
     None),
    (Number Field in c2 with defining polynomial x^2 + 512,
     Ring morphism:
       From: Number Field in c2 with defining polynomial x^2 + 512
       To:   Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
       Defn: c2 |--> 1/25*c^6 + 78/25*c^2,
     None),
    (Number Field in c3 with defining polynomial x^2 - 288,
     Ring morphism:
       From: Number Field in c3 with defining polynomial x^2 - 288
       To:   Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
       Defn: c3 |--> -1/25*c^6 + 22/25*c^2,
     None),
    (Number Field in c4 with defining polynomial x^4 + 112*x^2 + 40000,
     Ring morphism:
       From: Number Field in c4 with defining polynomial x^4 + 112*x^2 + 40000
       To:   Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
       Defn: c4 |--> 2*c^2,
     None),
    (Number Field in c5 with defining polynomial x^4 + 8,
    Ring morphism:
      From: Number Field in c5 with defining polynomial x^4 + 8
      To:   Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
      Defn: c5 |--> -1/80*c^5 + 1/40*c,
      None),
    (Number Field in c6 with defining polynomial x^4 + 648,
     Ring morphism:
       From: Number Field in c6 with defining polynomial x^4 + 648
       To:   Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
       Defn: c6 |--> 1/80*c^5 + 79/40*c,
     None),
    (Number Field in c7 with defining polynomial x^4 - 512,
     Ring morphism:
       From: Number Field in c7 with defining polynomial x^4 - 512
       To:   Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
       Defn: c7 |--> -1/60*c^5 + 41/30*c,
     None),
    (Number Field in c8 with defining polynomial x^4 - 32,
     Ring morphism:
       From: Number Field in c8 with defining polynomial x^4 - 32
       To:   Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
       Defn: c8 |--> 1/60*c^5 + 19/30*c,
     None),
    (Number Field in c9 with defining polynomial x^8 + 28*x^4 + 2500,
     Ring morphism:
       From: Number Field in c9 with defining polynomial x^8 + 28*x^4 + 2500
       To:   Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
       Defn: c9 |--> c,
     Ring morphism:
       From: Number Field in c with defining polynomial x^8 + 28*x^4 + 2500
       To:   Number Field in c9 with defining polynomial x^8 + 28*x^4 + 2500
       Defn: c |--> c9)
    ]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: tau = G[6]
    sage: tau_matrix = column_matrix([tau(be).vector() for be in basis])
    sage: (tau_matrix-identity_matrix(8)).right_kernel()
    Vector space of degree 8 and dimension 4 over Rational Field
    Basis matrix:
    [    1     0     0     0     0     0     0     0]
    [    0     1     0     0     0 -1/82     0     0]
    [    0     0     1     0     0     0 -1/22     0]
    [    0     0     0     1     0     0     0 11/58]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(tau(1))
    1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(tau(c+(-1/82)*c^5))
    -120/41*a

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(tau(c^2+(-1/22)*c^6))
    150/11*a^2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: fromL(tau(c^3+(11/58)*c^7))
    3000/29*a^3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: sg[2].is_normal(P)
    False

"""
