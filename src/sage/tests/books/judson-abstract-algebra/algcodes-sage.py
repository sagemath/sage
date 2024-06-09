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
## Section 8.8 Sage
##
r"""
~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: codes.   # not tested

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H = codes.HammingCode(GF(2), 3); H
    [7, 4] Hamming Code over GF(2)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.dimension()
    4

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.list()
    [(0, 0, 0, 0, 0, 0, 0), (1, 0, 0, 0, 0, 1, 1), (0, 1, 0, 0, 1, 0, 1),
     (1, 1, 0, 0, 1, 1, 0), (0, 0, 1, 0, 1, 1, 0), (1, 0, 1, 0, 1, 0, 1),
     (0, 1, 1, 0, 0, 1, 1), (1, 1, 1, 0, 0, 0, 0), (0, 0, 0, 1, 1, 1, 1),
     (1, 0, 0, 1, 1, 0, 0), (0, 1, 0, 1, 0, 1, 0), (1, 1, 0, 1, 0, 0, 1),
     (0, 0, 1, 1, 0, 0, 1), (1, 0, 1, 1, 0, 1, 0), (0, 1, 1, 1, 1, 0, 0),
     (1, 1, 1, 1, 1, 1, 1)]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.minimum_distance()
    3

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C = H.parity_check_matrix(); C
    [1 0 1 0 1 0 1]
    [0 1 1 0 0 1 1]
    [0 0 0 1 1 1 1]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: G = H.generator_matrix(); G
    [1 0 0 0 0 1 1]
    [0 1 0 0 1 0 1]
    [0 0 1 0 1 1 0]
    [0 0 0 1 1 1 1]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C*G.transpose() == zero_matrix(3, 4)
    True

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.systematic_generator_matrix()
    [1 0 0 0 0 1 1]
    [0 1 0 0 1 0 1]
    [0 0 1 0 1 1 0]
    [0 0 0 1 1 1 1]

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: r = vector(GF(2), [1, 1, 1, 1, 0, 0, 1]); r
    (1, 1, 1, 1, 0, 0, 1)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: C*r
    (1, 1, 0)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.decode_to_code(r)
    (1, 1, 0, 1, 0, 0, 1)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: message = vector(GF(2), [1, 1, 0, 1, 0, 0, 1])
    sage: errors = vector(GF(2), [0, 0, 1, 0, 1, 1, 0])
    sage: received = message + errors
    sage: received
    (1, 1, 1, 1, 1, 1, 1)

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.decode_to_code(received) == message
    False

~~~~~~~~~~~~~~~~~~~~~~ ::

    sage: H.decode_to_code(received) == received
    True
"""
