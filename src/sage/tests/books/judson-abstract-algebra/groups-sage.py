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
## Section 3.7 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z8 = Integers(8)
    sage: Z8
    Ring of integers modulo 8

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z8.list()
    [0, 1, 2, 3, 4, 5, 6, 7]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = Z8.an_element(); a
    0

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a.parent()
    Ring of integers modulo 8

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 6
    sage: a
    6

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a.parent()
    Integer Ring

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b = 7
    sage: c = a + b; c
    13

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: d = Z8(6)
    sage: d
    6

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: d.parent()
    Ring of integers modulo 8

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: e = Z8(7)
    sage: f = d+e; f
    5

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: g = Z8(85); g
    5

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: f == g
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Z8.addition_table(names='elements')
    +  0 1 2 3 4 5 6 7
     +----------------
    0| 0 1 2 3 4 5 6 7
    1| 1 2 3 4 5 6 7 0
    2| 2 3 4 5 6 7 0 1
    3| 3 4 5 6 7 0 1 2
    4| 4 5 6 7 0 1 2 3
    5| 5 6 7 0 1 2 3 4
    6| 6 7 0 1 2 3 4 5
    7| 7 0 1 2 3 4 5 6

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: triangle = SymmetricGroup(3)
    sage: rho2 = triangle([3,1,2])
    sage: rho2
    (1,3,2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [rho2(x) for x in triangle.domain()]
    [3, 1, 2]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [[a(x) for x in triangle.domain()] for a in triangle]
    [[1, 2, 3], [3, 1, 2], [2, 3, 1], [1, 3, 2], [3, 2, 1], [2, 1, 3]]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: mu1 = triangle([1,3,2])
    sage: mu2 = triangle([3,2,1])
    sage: mu3 = triangle([2,1,3])
    sage: rho1 = triangle([2,3,1])
    sage: product = rho1*mu1
    sage: product == mu2
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [product(x) for x in triangle.domain()]
    [3, 2, 1]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: rho1*mu1 == mu1*rho1
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: mu1*rho1 == mu3
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: triangle.cayley_table()
    *  a b c d e f
     +------------
    a| a b c d e f
    b| b a d c f e
    c| c e a f b d
    d| d f b e a c
    e| e c f a d b
    f| f d e b c a

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: triangle.cayley_table(names='elements')
          *       ()   (2,3)   (1,2) (1,2,3) (1,3,2)   (1,3)
           +------------------------------------------------
         ()|      ()   (2,3)   (1,2) (1,2,3) (1,3,2)   (1,3)
      (2,3)|   (2,3)      () (1,2,3)   (1,2)   (1,3) (1,3,2)
      (1,2)|   (1,2) (1,3,2)      ()   (1,3)   (2,3) (1,2,3)
    (1,2,3)| (1,2,3)   (1,3)   (2,3) (1,3,2)      ()   (1,2)
    (1,3,2)| (1,3,2)   (1,2)   (1,3)      () (1,2,3)   (2,3)
      (1,3)|   (1,3) (1,2,3) (1,3,2)   (2,3)   (1,2)      ()

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: triangle.cayley_table(names=['id','u1','u3','r1','r2','u2'])
     *  id u1 u3 r1 r2 u2
      +------------------
    id| id u1 u3 r1 r2 u2
    u1| u1 id r1 u3 u2 r2
    u3| u3 r2 id u2 u1 r1
    r1| r1 u2 u1 r2 id u3
    r2| r2 u3 u2 id r1 u1
    u2| u2 r1 r2 u1 u3 id

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q = QuaternionGroup()
    sage: [[a(x) for x in Q.domain()] for a in Q]
    [[1, 2, 3, 4, 5, 6, 7, 8],
     [3, 4, 1, 2, 7, 8, 5, 6],
     [4, 1, 2, 3, 8, 5, 6, 7],
     [2, 3, 4, 1, 6, 7, 8, 5],
     [7, 6, 5, 8, 1, 4, 3, 2],
     [5, 8, 7, 6, 3, 2, 1, 4],
     [8, 7, 6, 5, 2, 1, 4, 3],
     [6, 5, 8, 7, 4, 3, 2, 1]]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.cayley_table()
    *  a b c d e f g h
     +----------------
    a| a b c d e f g h
    b| b c d a h e f g
    c| c d a b g h e f
    d| d a b c f g h e
    e| e f g h c d a b
    f| f g h e b c d a
    g| g h e f a b c d
    h| h e f g d a b c

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: id = Q.identity()
    sage: [id(x) for x in Q.domain()]
    [1, 2, 3, 4, 5, 6, 7, 8]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: minus_one = Q([3, 4, 1, 2, 7, 8, 5, 6])
    sage: minus_one*minus_one == Q.identity()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.is_finite()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.order()
    8

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.is_abelian()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S8 = SymmetricGroup(8)
    sage: a = S8.random_element()
    sage: [a(x) for x in S8.domain()]     # random
    [5, 2, 6, 4, 1, 8, 3, 7]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S8.order()
    40320

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q.is_subgroup(S8)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H = [CC(1), CC(-1), CC(I), CC(-I)]
    sage: CC.multiplication_table(elements=H,
    ....:                         names=['1', '-1', 'i', '-i'])
    *   1 -1  i -i
      +------------
     1|  1 -1  i -i
    -1| -1  1 -i  i
     i|  i -i -1  1
    -i| -i  i  1 -1
"""
