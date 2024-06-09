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
## Section 10.4 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = DihedralGroup(8)
    sage: quarter_turn = G('(1,3,5,7)(2,4,6,8)')
    sage: S = G.subgroup([quarter_turn])
    sage: C = G.cosets(S)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p = C[1][0]*C[3][0]
    sage: [i for i in srange(len(C)) if p in C[i]]
    [2]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p = C[1][2]*C[3][1]
    sage: [i for i in srange(len(C)) if p in C[i]]
    [2]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: def coset_product(i, j, C):
    ....:    p = C[i][0]*C[j][0]
    ....:    c = [k for k in srange(len(C)) if p in C[k]]
    ....:    return c[0]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: coset_product(1, 3, C)
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = AlternatingGroup(4)
    sage: face_turn = G("(1,2,3)")
    sage: S = G.subgroup([face_turn])
    sage: C = G.cosets(S)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p = C[1][0]*C[3][0]
    sage: [i for i in srange(len(C)) if p in C[i]]
    [0]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: p = C[1][0]*C[3][1]
    sage: [i for i in srange(len(C)) if p in C[i]]
    [2]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = DihedralGroup(8)
    sage: quarter_turn = G('(1,3,5,7)(2,4,6,8)')
    sage: S = G.subgroup([quarter_turn])
    sage: S.is_normal(G)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = AlternatingGroup(4)
    sage: face_turn = G("(1,2,3)")
    sage: S = G.subgroup([face_turn])
    sage: S.is_normal(G)
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = AlternatingGroup(5)
    sage: G.is_simple()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = DihedralGroup(8)
    sage: quarter_turn = G('(1,3,5,7)(2,4,6,8)')
    sage: S = G.subgroup([quarter_turn])
    sage: Q = G.quotient(S)
    sage: Q
    Permutation Group with generators [(1,2)(3,4), (1,3)(2,4)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.is_isomorphic(KleinFourGroup())
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = DihedralGroup(8)
    sage: N = G.normal_subgroups()
    sage: l=[H.order() for H in N]; l.sort(); l
    [1, 2, 4, 8, 8, 8, 16]
"""
