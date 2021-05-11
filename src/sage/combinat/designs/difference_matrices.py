r"""
Difference Matrices

This module gathers code related to difference matrices. One can build those
objects (or know if they can be built) with :func:`difference_matrix`::

    sage: G,DM = designs.difference_matrix(9,5,1)

Functions
---------
"""

from sage.misc.unknown import Unknown
from sage.misc.cachefunc import cached_function
from sage.categories.sets_cat import EmptySetError
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.arith.all import is_prime_power, divisors
from .designs_pyx import is_difference_matrix
from .database import DM as DM_constructions

@cached_function
def find_product_decomposition(g, k, lmbda=1):
    r"""
    Try to find a product decomposition construction for difference matrices.

    INPUT:

    - ``g,k,lmbda`` -- integers, parameters of the difference matrix

    OUTPUT:

    A pair of pairs ``(g1,lmbda),(g2,lmbda2)`` if Sage knows how to build
    `(g1,k,lmbda1)` and `(g2,k,lmbda2)` difference matrices and ``False``
    otherwise.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_matrices import find_product_decomposition
        sage: find_product_decomposition(77,6)
        ((7, 1), (11, 1))
        sage: find_product_decomposition(616,7)
        ((7, 1), (88, 1))
        sage: find_product_decomposition(24,10)
        False
    """
    for lmbda1 in divisors(lmbda):
        lmbda2 = lmbda//lmbda1

        # To avoid infinite loop:
        # if lmbda1 == lmbda, then g1 should not be g
        # if lmbda2 == lmbda, then g2 should not be g
        if lmbda1 == lmbda:
            if lmbda2 == lmbda:
                div = divisors(g)[1:-1]
            else:
                div = divisors(g)[:-1]
        else:
            if lmbda2 == lmbda:
                div = divisors(g)[1:]
            else:
                div = divisors(g)

        for g1 in div:
            g2 = g//g1
            if g1 > g2:
                break
            if (difference_matrix(g1,k,lmbda1,existence=True) is True and
                difference_matrix(g2,k,lmbda2,existence=True) is True):
                return (g1,lmbda1),(g2,lmbda2)

    return False

def difference_matrix_product(k, M1, G1, lmbda1, M2, G2, lmbda2, check=True):
    r"""
    Return the product of the ``(G1,k,lmbda1)`` and ``(G2,k,lmbda2)`` difference
    matrices ``M1`` and ``M2``.

    The result is a `(G1 \times G2, k, \lambda_1 \lambda_2)`-difference matrix.

    INPUT:

    - ``k,lmbda1,lmbda2`` -- positive integer

    - ``G1, G2`` -- groups

    - ``M1, M2`` -- ``(G1,k,lmbda1)`` and ``(G,k,lmbda2)`` difference
      matrices

    - ``check`` (boolean) -- if ``True`` (default), the output is checked before
      being returned.

    EXAMPLES::

        sage: from sage.combinat.designs.difference_matrices import (
        ....:     difference_matrix_product,
        ....:     is_difference_matrix)
        sage: G1,M1 = designs.difference_matrix(11,6)
        sage: G2,M2 = designs.difference_matrix(7,6)
        sage: G,M = difference_matrix_product(6,M1,G1,1,M2,G2,1)
        sage: G1
        Finite Field of size 11
        sage: G2
        Finite Field of size 7
        sage: G
        The Cartesian product of (Finite Field of size 11, Finite Field of size 7)
        sage: is_difference_matrix(M,G,6,1)
        True
    """
    g1 = G1.cardinality()
    g2 = G2.cardinality()
    g = g1*g2
    lmbda = lmbda1*lmbda2
    from sage.categories.cartesian_product import cartesian_product
    G = cartesian_product([G1,G2])

    M = [[G((M1[j1][i],M2[j2][i])) for i in range(k)] for j1 in range(lmbda1*g1) for j2 in range(lmbda2*g2)]

    if check and not is_difference_matrix(M,G,k,lmbda,True):
        raise RuntimeError("In the product construction, Sage built something which is not a ({},{},{})-DM!".format(g,k,lmbda))

    return G,M

def difference_matrix(g,k,lmbda=1,existence=False,check=True):
    r"""
    Return a `(g,k,\lambda)`-difference matrix

    A matrix `M` is a `(g,k,\lambda)`-difference matrix if it has size `\lambda
    g\times k`, its entries belong to the group `G` of cardinality `g`, and
    for any two rows `R,R'` of `M` and `x\in G` there are exactly `\lambda`
    values `i` such that `R_i-R'_i=x`.

    INPUT:

    - ``k`` -- (integer) number of columns. If ``k=None`` it is set to the
      largest value available.

    - ``g`` -- (integer) cardinality of the group `G`

    - ``lmbda`` -- (integer; default: 1) -- number of times each element of `G`
      appears as a difference.

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to
      ``True`` by default.

    - ``existence`` (boolean) -- instead of building the design, return:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that the design does not exist.

      .. NOTE::

          When ``k=None`` and ``existence=True`` the function returns an
          integer, i.e. the largest `k` such that we can build a
          `(g,k,\lambda)`-DM.

    EXAMPLES::

        sage: G,M = designs.difference_matrix(25,10); G
        Finite Field in x of size 5^2
        sage: designs.difference_matrix(993,None,existence=1)
        32

    Here we print for each `g` the maximum possible `k` for which Sage knows
    how to build a `(g,k,1)`-difference matrix::

        sage: for g in range(2,30):
        ....:     k_max = designs.difference_matrix(g=g,k=None,existence=True)
        ....:     print("{:2} {}".format(g, k_max))
        ....:     _ = designs.difference_matrix(g,k_max)
         2 2
         3 3
         4 4
         5 5
         6 2
         7 7
         8 8
         9 9
        10 2
        11 11
        12 6
        13 13
        14 2
        15 3
        16 16
        17 17
        18 2
        19 19
        20 4
        21 6
        22 2
        23 23
        24 8
        25 25
        26 2
        27 27
        28 6
        29 29

    TESTS::

        sage: designs.difference_matrix(10,12,1,existence=True)
        False
        sage: designs.difference_matrix(10,12,1)
        Traceback (most recent call last):
        ...
        EmptySetError: No (10,12,1)-Difference Matrix exists as k(=12)>g(=10)
        sage: designs.difference_matrix(10,9,1,existence=True)
        Unknown
        sage: designs.difference_matrix(10,9,1)
        Traceback (most recent call last):
        ...
        NotImplementedError: I don't know how to build a (10,9,1)-Difference Matrix!
    """

    if lmbda == 1 and k is not None and k>g:
        if existence:
            return False
        raise EmptySetError("No ({},{},{})-Difference Matrix exists as k(={})>g(={})".format(g,k,lmbda,k,g))

    # Prime powers
    elif lmbda == 1 and is_prime_power(g):
        if k is None:
            if existence:
                return g
            else:
                k = g
        elif existence:
            return True
        F       = FiniteField(g,'x')
        F_set   = list(F)
        F_k_set = F_set[:k]

        G = F
        M = [[x*y for y in F_k_set] for x in F_set]

    # Treat the case k=None
    # (find the max k such that there exists a DM)
    elif k is None:
        i = 2
        while difference_matrix(g=g,k=i,lmbda=lmbda,existence=True) is True:
            i += 1
        return i-1

    # From the database
    elif (g,lmbda) in DM_constructions and DM_constructions[g,lmbda][0]>=k:
        if existence:
            return True
        _,f = DM_constructions[g,lmbda]
        G, M = f()
        M = [R[:k] for R in M]

    # Product construction
    elif find_product_decomposition(g,k,lmbda):
        if existence:
            return True
        (g1,lmbda1),(g2,lmbda2) = find_product_decomposition(g,k,lmbda)
        G1,M1 = difference_matrix(g1,k,lmbda1)
        G2,M2 = difference_matrix(g2,k,lmbda2)
        G,M = difference_matrix_product(k,M1,G1,lmbda1,M2,G2,lmbda2,check=False)

    else:
        if existence:
            return Unknown
        raise NotImplementedError("I don't know how to build a ({},{},{})-Difference Matrix!".format(g,k,lmbda))

    if check and not is_difference_matrix(M,G,k,lmbda,1):
        raise RuntimeError("Sage built something which is not a ({},{},{})-DM!".format(g,k,lmbda))

    return G,M
