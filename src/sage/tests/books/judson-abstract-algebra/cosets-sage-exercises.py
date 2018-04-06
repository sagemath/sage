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
## Exercises 6.6 Sage Exercises
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 5
    sage: b = 10
    sage: c = 6
    sage: d = 13
    sage: a.divides(b)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: not (b in [c,d])
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a.divides(b) and not (b in [c,d])
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 20
    sage: b = 6
    sage: a.mod(b)
    2

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: prime_range(50, 100)
    [53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: all([True, True, True, True])
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: all([True, True, False, True])
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [a/b for a in srange(9) for b in srange(1,a) if gcd(a,b)==1]
    [2, 3, 3/2, 4, 4/3, 5, 5/2, 5/3, 5/4, 6, 6/5,
     7, 7/2, 7/3, 7/4, 7/5, 7/6, 8, 8/3, 8/5, 8/7]

"""
