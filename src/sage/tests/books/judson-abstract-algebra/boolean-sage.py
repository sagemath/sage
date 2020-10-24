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
## Section 19.7 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: X = (24).divisors()
    sage: X
    [1, 2, 3, 4, 6, 8, 12, 24]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: R = [(a,b) for a in X for b in X if a.divides(b)]; R
    [(1, 1), (1, 2), (1, 3), (1, 4), (1, 6), (1, 8), (1, 12), (1, 24),
     (2, 2), (2, 4), (2, 6), (2, 8), (2, 12), (2, 24), (3, 3), (3, 6),
     (3, 12), (3, 24), (4, 4), (4, 8), (4, 12), (4, 24), (6, 6),
     (6, 12), (6, 24), (8, 8), (8, 24), (12, 12), (12, 24), (24, 24)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: D = Poset([X, R])
    sage: D.plot()   # not tested
    

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: divisible = lambda x, y: x.divides(y)
    sage: L = Poset([X, divisible])
    sage: L == D
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: L.plot()   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Q = Posets.PentagonPoset()
    sage: Q.plot()   # not tested
    

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: S = Posets.BooleanLattice(4)
    sage: S.plot()   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: T = Posets.RandomPoset(20,0.05)
    sage: T.plot()   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: D.is_lequal(4, 8)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: D.is_lequal(4, 4)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: D.is_less_than(4, 8)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: D.is_less_than(4, 4)
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: D.is_lequal(6, 8)
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: D.is_lequal(8, 6)
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: D.is_gequal(8, 4)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: D.is_greater_than(4, 8)
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: X = range(20)
    sage: C = [[18, 7],  [9, 11], [9, 10], [11, 8], [6, 10],
    ....:      [10, 2],   [0, 2],  [2, 1],  [1, 8], [8, 12],
    ....:      [8, 3],  [3, 15], [15, 7], [7, 16],  [7, 4],
    ....:      [16, 17], [16, 13], [4, 19], [4, 14], [14, 5]]
    sage: P = Poset([X, C])
    sage: P.plot()   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: P.minimal_elements()
    [18, 9, 6, 0]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: P.maximal_elements()
    [5, 19, 13, 17, 12]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: P.level_sets()
    [[18, 9, 6, 0], [11, 10], [2], [1], [8], [3, 12],
     [15], [7], [16, 4], [13, 17, 14, 19], [5]]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: linear = P.linear_extension(); linear
    [18, 9, 11, 6, 10, 0, 2, 1, 8, 3, 15, 
     7, 4, 14, 5, 19, 16, 13, 17, 12]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: level = P.level_sets()
    sage: bottomhalf = sum([level[i] for i in range(5)], [])
    sage: B = P.subposet(bottomhalf)
    sage: B.plot()   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Pdual = P.dual()
    sage: Pdual.plot()   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: Ddual = D.dual()
    sage: Ddual.plot()   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: P = Posets.AntichainPoset(8)
    sage: P.is_lattice()
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: LatticePoset(P)
    Traceback (most recent call last):
    ...
    ValueError: not a meet-semilattice: no bottom element

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: CP = Posets.IntegerCompositions(5)
    sage: C = LatticePoset(CP)
    sage: C.plot()   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: par = C.an_element().parent()
    sage: a = par([1, 1, 1, 2])
    sage: b = par([2, 1, 1, 1])
    sage: a, b
    ([1, 1, 1, 2], [2, 1, 1, 1])

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C.meet(a, b)
    [2, 1, 2]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: c = par([1, 4])
    sage: d = par([2, 3])
    sage: c, d
    ([1, 4], [2, 3])

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C.join(c, d)
    [1, 1, 3]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C.is_distributive()
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C.top()
    [1, 1, 1, 1, 1]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C.bottom()
    [5]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: comp = C.complements()
    sage: comp[par([1, 1, 1, 2])]
    [[4, 1]]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: comp   # not tested
    {[1, 1, 1, 1, 1]: [[5]],
     [1, 1, 1, 2]: [[4, 1]],
     [1, 1, 2, 1]: [[3, 2]],
     [1, 1, 3]: [[3, 1, 1]],
     [1, 2, 1, 1]: [[2, 3]],
     [1, 2, 2]: [[2, 2, 1]],
     [1, 3, 1]: [[2, 1, 2]],
     [1, 4]: [[2, 1, 1, 1]],
     [2, 1, 1, 1]: [[1, 4]],
     [2, 1, 2]: [[1, 3, 1]],
     [2, 2, 1]: [[1, 2, 2]],
     [2, 3]: [[1, 2, 1, 1]],
     [3, 1, 1]: [[1, 1, 3]],
     [3, 2]: [[1, 1, 2, 1]],
     [4, 1]: [[1, 1, 1, 2]],
     [5]: [[1, 1, 1, 1, 1]]}

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: [len(e[1]) for e in comp.items()]
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C.is_complemented()
    True

"""
