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
## Section 13.6 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S = SymmetricGroup(3)
    sage: D = DihedralGroup(3)
    sage: S.is_isomorphic(D)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C3 = CyclicPermutationGroup(3)
    sage: C5 = CyclicPermutationGroup(5)
    sage: DP = direct_product_permgroups([C3, C5])
    sage: C  = CyclicPermutationGroup(15)
    sage: DP.is_isomorphic(C)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q  = QuaternionGroup()
    sage: DI = DiCyclicGroup(2)
    sage: Q.is_isomorphic(DI)
    True

"""
