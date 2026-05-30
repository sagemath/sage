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
## Section 1.5 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::


~~~~~~~~~~~~~~~~~~~~~~ ::


~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: a = 10
    sage: b = 6
    sage: b = b - 10
    sage: a = a + 20
    sage: a
    30

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b = b + 50

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: b + 20
    66

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A = matrix([[3, 1], [5,2]]); A
    [3 1]
    [5 2]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: print(A); print(); print(A.inverse())
    [3 1]
    [5 2]
    <BLANKLINE>
    [ 2 -1]
    [-5  3]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A.   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A.inverse?   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: A.inverse??   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: B = matrix([[2, 20], [5, 50]])
    sage: B.inverse()
    Traceback (most recent call last):
    ...
    ZeroDivisionError: matrix must be nonsingular

~~~~~~~~~~~~~~~~~~~~~~ ::


~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: zoo = ['snake', 'parrot', 'elephant', 'baboon', 'beetle']
    sage: zoo
    ['snake', 'parrot', 'elephant', 'baboon', 'beetle']

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: zoo[2]
    'elephant'

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: zoo.append('ostrich'); zoo
    ['snake', 'parrot', 'elephant', 'baboon', 'beetle', 'ostrich']

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: zoo.remove('parrot')
    sage: zoo
    ['snake', 'elephant', 'baboon', 'beetle', 'ostrich']

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: mammals = zoo[1:3]
    sage: mammals
    ['elephant', 'baboon']

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: newzoo = sorted(zoo)
    sage: newzoo
    ['baboon', 'beetle', 'elephant', 'ostrich', 'snake']

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: zoo.sort()
    sage: zoo
    ['baboon', 'beetle', 'elephant', 'ostrich', 'snake']

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: plurality_zoo = [animal+'s' for animal in zoo]
    sage: plurality_zoo
    ['baboons', 'beetles', 'elephants', 'ostrichs', 'snakes']

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: dozen = srange(12); dozen
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: teens = srange(13, 20); teens
    [13, 14, 15, 16, 17, 18, 19]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: decades = srange(1900, 2000, 10); decades
    [1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990]
"""
