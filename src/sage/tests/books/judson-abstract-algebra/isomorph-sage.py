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
## Section 9.4 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: m = 12
    sage: n = 7
    sage: gcd(m, n)
    1

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = CyclicPermutationGroup(m)
    sage: H = CyclicPermutationGroup(n)
    sage: dp = direct_product_permgroups([G, H])
    sage: K = CyclicPermutationGroup(m*n)
    sage: K.is_isomorphic(dp)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: m = 15
    sage: n = 21
    sage: gcd(m, n)
    3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = CyclicPermutationGroup(m)
    sage: H = CyclicPermutationGroup(n)
    sage: dp = direct_product_permgroups([G, H])
    sage: K = CyclicPermutationGroup(m*n)
    sage: K.is_isomorphic(dp)
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: m = 6
    sage: n = 5
    sage: r = 7
    sage: G = CyclicPermutationGroup(m)
    sage: H = CyclicPermutationGroup(n)
    sage: L = CyclicPermutationGroup(r)
    sage: dp = direct_product_permgroups([G, H, L])
    sage: K = CyclicPermutationGroup(m*n*r)
    sage: K.is_isomorphic(dp)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [CyclicPermutationGroup(p) for p in [1, 2, 3, 5, 7, 11, 13]]
    [Cyclic group of order 1 as a permutation group,
     Cyclic group of order 2 as a permutation group,
     Cyclic group of order 3 as a permutation group,
     Cyclic group of order 5 as a permutation group,
     Cyclic group of order 7 as a permutation group,
     Cyclic group of order 11 as a permutation group,
     Cyclic group of order 13 as a permutation group]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = CyclicPermutationGroup(4)
    sage: H = KleinFourGroup()
    sage: T1 = CyclicPermutationGroup(2)
    sage: T2 = CyclicPermutationGroup(2)
    sage: K = direct_product_permgroups([T1, T2])
    sage: G.is_isomorphic(H)
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.is_isomorphic(K)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = CyclicPermutationGroup(6)
    sage: H = SymmetricGroup(3)
    sage: G.is_isomorphic(H)
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = DihedralGroup(3)
    sage: H = SymmetricGroup(3)
    sage: G.is_isomorphic(H)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z36 = Integers(36)
    sage: U = [x for x in Z36 if gcd(ZZ(x), 36) == 1]
    sage: U
    [1, 5, 7, 11, 13, 17, 19, 23, 25, 29, 31, 35]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [x.multiplicative_order() for x in U]
    [1, 6, 6, 6, 3, 2, 2, 6, 3, 6, 6, 2]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = U[1]
    sage: A = [a^i for i in srange(6)]
    sage: A
    [1, 5, 25, 17, 13, 29]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b = U[11]
    sage: B = [b^i for i in srange(2)]
    sage: B
    [1, 35]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [x for x in A if x in B]
    [1]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: all(x*y == y*x for x in A for y in B)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: T = sorted(x*y for x in A for y in B)
    sage: T
    [1, 5, 7, 11, 13, 17, 19, 23, 25, 29, 31, 35]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: T == U
    True

"""
