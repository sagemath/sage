"""
Optimized computing of relation matrices in certain cases
"""

#############################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL) v2+.
#  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
#############################################################################

from sage.misc.verbose import verbose
from sage.rings.rational cimport Rational


def sparse_2term_quotient_only_pm1(rels, n):
    r"""
    Perform Sparse Gauss elimination on a matrix all of whose columns
    have at most 2 nonzero entries with relations all 1 or -1.

    This algorithm is more subtle than just "identify symbols in pairs",
    since complicated relations can cause generators to equal 0.

    .. NOTE::

        Note the condition on the s,t coefficients in the relations
        being 1 or -1 for this optimized function.  There is a more
        general function in relation_matrix.py, which is much, much
        slower.

    INPUT:

    - ``rels`` -- iterable made of pairs ((i,s), (j,t)). The pair
      represents the relation `s x_i + t x_j = 0`, where the `i, j` must
      be Python int's, and the `s, t` must all be 1 or -1.

    - ``n`` -- integer; the `x_i` are `x_0, \ldots, x_{n-1}`

    OUTPUT:

    - ``mod`` -- list such that mod[i] = (j,s), which means that x_i
      is equivalent to s*x_j, where the x_j are a basis for the
      quotient.

    The output depends on the order of the input.

    EXAMPLES::

        sage: from sage.modular.modsym.relation_matrix_pyx import sparse_2term_quotient_only_pm1
        sage: rels = [((0,1), (1,-1)), ((1,1), (3,1)), ((2,1),(3,1)), ((4,1),(5,-1))]
        sage: n = 6
        sage: sparse_2term_quotient_only_pm1(rels, n)
        [(3, -1), (3, -1), (3, -1), (3, 1), (5, 1), (5, 1)]
    """
    n = int(n)

    tm = verbose("Starting optimized integer sparse 2-term quotient...")

    cdef int c0, c1, i, die
    cdef list free = list(range(n))
    cdef list coef = [1] * n
    cdef list related_to_me = [[] for i in range(n)]

    for v0, v1 in rels:
        c0 = coef[v0[0]] * v0[1]
        c1 = coef[v1[0]] * v1[1]

        # Mod out by the following relation:
        #
        #    c0*free[v0[0]] + c1*free[v1[0]] = 0.
        #
        die = -1
        if c0 == 0 and c1 == 0:
            pass
        elif c0 == 0 and c1 != 0:  # free[v1[0]] --> 0
            die = free[v1[0]]
        elif c1 == 0 and c0 != 0:
            die = free[v0[0]]
        elif free[v0[0]] == free[v1[0]]:
            if c0 + c1 != 0:
                # all xi equal to free[v0[0]] must now equal to zero.
                die = free[v0[0]]
        else:  # x1 = -c1/c0 * x2.
            x = free[v0[0]]
            free[x] = free[v1[0]]
            if c0 != 1 and c0 != -1:
                raise ValueError("coefficients must all be -1 or 1.")
            coef[x] = -c1 * c0
            for i in related_to_me[x]:
                free[i] = free[x]
                coef[i] *= coef[x]
                related_to_me[free[v1[0]]].append(i)
            related_to_me[free[v1[0]]].append(x)
        if die != -1:
            for i in related_to_me[die]:
                free[i] = 0
                coef[i] = 0
            free[die] = 0
            coef[die] = 0

    # Special casing the rationals leads to a huge speedup,
    # actually.  (All the code above is slower than just this line
    # without this special case.)
    mod = [(fi, Rational(ci)) for fi, ci in zip(free, coef)]

    verbose("finished", tm)
    return mod
