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
## Generated at 2017-08-24T11:43:34-07:00
## From "Abstract Algebra"
## At commit 26d3cac0b4047f4b8d6f737542be455606e2c4b4
##
## Exercises 11.6 Sage Exercises
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = CyclicPermutationGroup(3)
    sage: H = DihedralGroup(4)
    sage: results = G.direct_product(H)
    sage: phi = results[2]
    sage: H.gens()
    ((1,2,3,4), (1,4)(2,3))

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = H.gen(0); a
    (1,2,3,4)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: phi(a)
    (4,5,6,7)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = DihedralGroup(20)
    sage: l=[H.order() for H in G.normal_subgroups()]; l.sort(); l
    [1, 2, 4, 5, 10, 20, 20, 20, 40]
"""
