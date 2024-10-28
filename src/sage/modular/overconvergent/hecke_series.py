# sage.doctest: needs sage.libs.flint sage.libs.pari
r"""
Atkin/Hecke series for overconvergent modular forms

This file contains a function :func:`~hecke_series` to compute the
characteristic series `P(t)` modulo `p^m` of the Atkin/Hecke operator `U_p`
upon the space of `p`-adic overconvergent modular forms of level `\Gamma_0(N)`.
The input weight ``k`` can also be a list ``klist`` of weights which must all
be congruent modulo `(p-1)`.

Two optional parameters ``modformsring`` and ``weightbound`` can be specified,
and in most cases for levels `N > 1` they can be used to obtain the output more
quickly. When `m \le k-1` the output `P(t)` is also equal modulo `p^m` to the
reverse characteristic polynomial of the Atkin operator `U_p` on the space of
classical modular forms of weight `k` and level `\Gamma_0(Np)`. In addition,
provided `m \le (k-2)/2` the output `P(t)` is equal modulo `p^m` to the reverse
characteristic polynomial of the Hecke operator `T_p` on the space of classical
modular forms of weight `k` and level `\Gamma_0(N)`. The function is based upon
the main algorithm in [Lau2011]_, and has linear running time in the logarithm of
the weight `k`.

AUTHORS:

- Alan G.B. Lauder (2011-11-10): original implementation.
- David Loeffler (2011-12): minor optimizations in review stage.

EXAMPLES:

The characteristic series of the `U_{11}` operator modulo `11^{10}` on the space of 11-adic overconvergent
modular forms of level 1 and weight 10000::

    sage: hecke_series(11, 1, 10000, 10)
    10009319650*x^4 + 25618839103*x^3 + 6126165716*x^2 + 10120524732*x + 1

The characteristic series of the U_5 operator modulo 5^5 on the space of 5-adic overconvergent
modular forms of level 3 and weight 1000::

    sage: hecke_series(5, 3, 1000, 5)
    1875*x^6 + 1250*x^5 + 1200*x^4 + 1385*x^3 + 1131*x^2 + 2533*x + 1

The characteristic series of the `U_7` operator modulo `7^5` on the space of 7-adic overconvergent
modular forms of level 5 and weight 1000. Here the optional parameter ``modformsring`` is set to ``True``::

    sage: hecke_series(7, 5, 1000, 5, modformsring=True)  # long time (21s on sage.math, 2012)
    12005*x^7 + 10633*x^6 + 6321*x^5 + 6216*x^4 + 5412*x^3 + 4927*x^2 + 4906*x + 1

The characteristic series of the `U_{13}` operator modulo `13^5` on the space of 13-adic overconvergent
modular forms of level 2 and weight 10000. Here the optional parameter ``weightbound`` is set to 4::

    sage: hecke_series(13, 2, 10000, 5, weightbound=4)  # long time (17s on sage.math, 2012)
    325156*x^5 + 109681*x^4 + 188617*x^3 + 220858*x^2 + 269566*x + 1

A list containing the characteristic series of the `U_{23}` operator modulo `23^{10}` on the spaces of
23-adic overconvergent modular forms of level 1 and weights 1000 and 1022, respectively.

::

    sage: hecke_series(23, 1, [1000, 1022], 10)
    [7204610645852*x^6 + 2117949463923*x^5 + 24152587827773*x^4 + 31270783576528*x^3 + 30336366679797*x^2
    + 29197235447073*x + 1, 32737396672905*x^4 + 36141830902187*x^3 + 16514246534976*x^2 + 38886059530878*x + 1]
"""

# ****************************************************************************
#       Copyright (C) 2011 Alan Lauder <lauder@maths.ox.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import valuation
from sage.functions.other import floor, ceil
from sage.matrix.constructor import matrix, random_matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.functional import dimension, transpose, charpoly
from sage.misc.timing import cputime
from sage.misc.verbose import verbose
from sage.modular.dims import dimension_modular_forms
from sage.modular.modform.all import ModularForms, ModularFormsRing, delta_qexp, eisenstein_series_qexp
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ

# AUXILIARY CODE: SPACES OF MODULAR FORMS AND LINEAR ALGEBRA

def compute_G(p, F):
    r"""
    Given a power series `F \in R[[q]]^\times`, for some ring `R`, and an
    integer `p`, compute the quotient

    .. MATH::

        \frac{F(q)}{F(q^p)}.

    Used by :func:`level1_UpGj` and by :func:`higher_level_UpGj`, with `F` equal
    to the Eisenstein series `E_{p-1}`.

    INPUT:

    - ``p`` -- integer
    - ``F`` -- power series (with invertible constant term)

    OUTPUT:

    the power series `F(q) / F(q^p)`, to the same precision as `F`

    EXAMPLES::

        sage: E = sage.modular.overconvergent.hecke_series.eisenstein_series_qexp(2, 12, Zmod(9),normalization='constant')
        sage: sage.modular.overconvergent.hecke_series.compute_G(3, E)
        1 + 3*q + 3*q^4 + 6*q^7 + O(q^12)
    """
    Fp = (F.truncate_powerseries(ceil(F.prec() / ZZ(p)))).V(p)
    return F / Fp


def low_weight_bases(N, p, m, NN, weightbound):
    r"""
    Return a list of integral bases of modular forms of level `N` and (even)
    weight at most ``weightbound``, as `q`-expansions modulo `(p^m,q^{NN})`.

    These forms are obtained by reduction mod `p^m` from an integral basis in
    Hermite normal form (so they are not necessarily in reduced row echelon
    form mod `p^m`, but they are not far off).

    INPUT:

    - ``N`` -- positive integer (level)
    - ``p`` -- prime
    - ``m``, ``NN`` -- positive integers
    - ``weightbound`` -- (even) positive integer

    OUTPUT: list of lists of `q`-expansions modulo `(p^m,q^{NN})`

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import low_weight_bases
        sage: low_weight_bases(2, 5, 3, 5, 6)
        [[1 + 24*q + 24*q^2 + 96*q^3 + 24*q^4 + O(q^5)],
        [1 + 115*q^2 + 35*q^4 + O(q^5), q + 8*q^2 + 28*q^3 + 64*q^4 + O(q^5)],
        [1 + 121*q^2 + 118*q^4 + O(q^5), q + 32*q^2 + 119*q^3 + 24*q^4 + O(q^5)]]
    """
    generators = []

    for k in range(2,weightbound + 2, 2):
        b = ModularForms(N, k, base_ring=Zmod(p ** m)).q_expansion_basis(prec=NN)
        generators.append(list(b))
    return generators

def random_low_weight_bases(N,p,m,NN,weightbound):
    r"""
    Return list of random integral bases of modular forms of level `N` and
    (even) weight at most weightbound with coefficients reduced modulo
    `(p^m,q^{NN})`.

    INPUT:

    - ``N`` -- positive integer (level)
    - ``p`` -- prime
    - ``m``, ``NN`` -- positive integers
    - ``weightbound`` -- (even) positive integer

    OUTPUT: list of lists of `q`-expansions modulo `(p^m,q^{NN})`

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import random_low_weight_bases
        sage: S = random_low_weight_bases(3, 7, 2, 5, 6); S # random
        [[4 + 48*q + 46*q^2 + 48*q^3 + 42*q^4 + O(q^5)],
        [3 + 5*q + 45*q^2 + 22*q^3 + 22*q^4 + O(q^5),
        1 + 3*q + 27*q^2 + 27*q^3 + 23*q^4 + O(q^5)],
        [2*q + 4*q^2 + 16*q^3 + 48*q^4 + O(q^5),
        2 + 6*q + q^2 + 3*q^3 + 43*q^4 + O(q^5),
        1 + 2*q + 6*q^2 + 14*q^3 + 4*q^4 + O(q^5)]]
        sage: S[0][0].parent()
        Power Series Ring in q over Ring of integers modulo 49
        sage: S[0][0].prec()
        5
    """
    LWB = low_weight_bases(N,p,m,NN,weightbound)
    # this is "approximately" row reduced (it's the mod p^n reduction of a
    # matrix over ZZ in Hermite form)
    RandomLWB = []
    for i in range(len(LWB)):
        n = len(LWB[i])
        c = random_matrix(Zmod(p ** m), n)
        while c.det() % p == 0:
            c = random_matrix(Zmod(p ** m), n)
        RandomLWB.append([ sum([c[j, k] * LWB[i][k] for k in range(n)]) for j in range(n) ])

    return RandomLWB

def low_weight_generators(N,p,m,NN):
    r"""
    Return a list of lists of modular forms, and an even natural number.

    The first output is a list of lists of modular forms reduced modulo
    `(p^m,q^{NN})` which generate the `(\ZZ / p^m \ZZ)`-algebra of mod `p^m`
    modular forms of weight at most 8, and the second output is the largest
    weight among the forms in the generating set.

    We (Alan Lauder and David Loeffler, the author and reviewer of this patch)
    conjecture that forms of weight at most 8 are always sufficient to generate
    the algebra of mod `p^m` modular forms of all weights. (We believe 6 to be
    sufficient, and we can prove that 4 is sufficient when there are no
    elliptic points, but using weights up to 8 acts as a consistency check.)

    INPUT:

    - ``N`` -- positive integer (level)
    - ``p`` -- prime
    - ``m``, ``NN`` -- positive integers

    OUTPUT: a tuple consisting of:

    - a list of lists of `q`-expansions modulo `(p^m,q^{NN})`,
    - an even natural number (twice the length of the list).

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import low_weight_generators
        sage: low_weight_generators(3, 7, 3, 10)
        ([[1 + 12*q + 36*q^2 + 12*q^3 + 84*q^4 + 72*q^5 + 36*q^6 + 96*q^7 + 180*q^8 + 12*q^9 + O(q^10)],
        [1 + 240*q^3 + 102*q^6 + 203*q^9 + O(q^10)],
        [1 + 182*q^3 + 175*q^6 + 161*q^9 + O(q^10)]], 6)
        sage: low_weight_generators(11, 5, 3, 10)
        ([[1 + 12*q^2 + 12*q^3 + 12*q^4 + 12*q^5 + 24*q^6 + 24*q^7 + 36*q^8 + 36*q^9 + O(q^10),
        q + 123*q^2 + 124*q^3 + 2*q^4 + q^5 + 2*q^6 + 123*q^7 + 123*q^9 + O(q^10)],
        [q + 116*q^4 + 115*q^5 + 102*q^6 + 121*q^7 + 96*q^8 + 106*q^9 + O(q^10)]], 4)
    """
    M = ModularFormsRing(N, base_ring=Zmod(p))

    b = M.gen_forms(maxweight=8)

    weightbound = max([f.weight() for f in b])
    generators = []

    for k in range(2, weightbound + 2, 2):
        generators.append([f.qexp(NN).change_ring(Zmod(p ** m)) for f in b if f.weight() == k])

    return generators, weightbound

def random_solution(B,K):
    r"""
    Return a random solution in nonnegative integers to the equation `a_1 + 2
    a_2 + 3 a_3 + ... + B a_B = K`, using a greedy algorithm.

    Note that this is *much* faster than using
    ``WeightedIntegerVectors.random_element()``.

    INPUT:

    - ``B``, ``K`` -- nonnegative integers

    OUTPUT: list

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import random_solution
        sage: s = random_solution(5, 10)
        sage: sum(s[i] * (i + 1) for i in range(5))
        10
        sage: S = set()
        sage: while len(S) != 30:
        ....:     s = random_solution(5, 10)
        ....:     assert sum(s[i] * (i + 1) for i in range(5)) == 10
        ....:     S.add(tuple(s))
    """
    a = []
    for i in range(B,1,-1):
        ai = ZZ.random_element((K // i) + 1)
        a.append(ai)
        K = K - ai*i
    a.append(K)
    a.reverse()

    return a


# AUXILIARY CODE: ECHELON FORM

def ech_form(A, p):
    r"""
    Return echelon form of matrix ``A`` over the ring of integers modulo
    `p^m`, for some prime `p` and `m \ge 1`.

    .. TODO::

        This should be moved to :mod:`sage.matrix.matrix_modn_dense` at some
        point.

    INPUT:

    - ``A`` -- matrix over ``Zmod(p^m)`` for some m
    - ``p`` -- prime p

    OUTPUT: matrix over ``Zmod(p^m)``

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import ech_form
        sage: A = MatrixSpace(Zmod(5 ** 3), 3)([1, 2, 3, 4, 5, 6, 7, 8, 9])
        sage: ech_form(A, 5)
        [1 2 3]
        [0 1 2]
        [0 0 0]
    """
    S = A[0, 0].parent()
    a = A.nrows()
    b = A.ncols()

    k = 0 # position pivoting row will be swapped to
    for j in range(b):
        if k < a:
            pivj = k # find new pivot
            for i in range(k + 1, a):
                if valuation(A[i, j], p) < valuation(A[pivj, j], p):
                    pivj = i
            if valuation(A[pivj, j], p) < +Infinity: # else column already reduced
                A.swap_rows(pivj, k)
                A.set_row_to_multiple_of_row(k, k, S(ZZ(A[k, j])/(p ** valuation(A[k, j], p))) ** (-1))
                for i in range(k + 1, a):
                    A.add_multiple_of_row(i, k, S(-ZZ(A[i, j]) / ZZ(A[k, j])))
                k = k + 1

    return A


# *** COMPLEMENTARY SPACES FOR LEVEL N > 1 ***

def random_new_basis_modp(N, p, k, LWBModp, TotalBasisModp, elldash, bound):
    r"""
    Return a list of lists of lists ``[j, a]`` encoding a choice of basis for
    the `i`-th complementary space `W_i`, as explained in the documentation for the
    function :func:`complementary_spaces_modp`.

    INPUT:

    - ``N`` -- positive integer at least 2 and not divisible by `p` (level)
    - ``p`` -- prime at least 5
    - ``k`` -- nonnegative integer
    - ``LWBModp`` -- list of list of `q`-expansions modulo `(p,q^\text{elldash})`
    - ``TotalBasisModp`` -- matrix over `\mathrm{GF}(p)`
    - ``elldash`` -- positive integer
    - ``bound`` -- positive even integer (twice the length of the list ``LWBModp``)

    OUTPUT: list of lists of lists ``[j, a]``

    .. NOTE::

        As well as having a non-trivial return value, this function also
        modifies the input matrix ``TotalBasisModp``.

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import random_low_weight_bases, complementary_spaces_modp
        sage: LWB = random_low_weight_bases(2, 5, 2, 4, 6)
        sage: LWBModp = [ [f.change_ring(GF(5)) for f in x] for x in LWB]
        sage: complementary_spaces_modp(2, 5, 2, 3, 4, LWBModp, 4) # random, indirect doctest
        [[[[0, 0]]], [[[0, 0], [1, 1]]], [[[0, 0], [1, 0], [1, 1]]], [[[0, 0], [1, 0], [1, 1], [1, 1]]]]
    """

    R = LWBModp[0][0].parent()

    # Case k0 + i(p-1) = 0 + 0(p-1) = 0

    if k == 0:
        TotalBasisModp[0,0] = 1
        return [[]]

    # Case k = k0 + i(p-1) > 0

    di = dimension_modular_forms(N, k)
    diminus1 = dimension_modular_forms(N, k - (p - 1))
    mi = di - diminus1

    NewBasisCode = []
    rk = diminus1
    for i in range(1, mi + 1):
        while (rk < diminus1 + i):
            # take random product of basis elements
            exps = random_solution(bound // 2, k // 2)
            TotalBasisi = R(1)
            TotalBasisiCode = []
            for j in range(len(exps)):
                for l in range(exps[j]):
                    a = ZZ.random_element(len(LWBModp[j]))
                    TotalBasisi = TotalBasisi * LWBModp[j][a]
                    TotalBasisiCode.append([j, a])
            TotalBasisModp[rk] = [TotalBasisi[j] for j in range(elldash)]
            TotalBasisModp.echelonize()
            rk = TotalBasisModp.rank()
        NewBasisCode.append(TotalBasisiCode) # this choice increased the rank

    return NewBasisCode

def complementary_spaces_modp(N, p, k0, n, elldash, LWBModp, bound):
    r"""
    Return a list of lists of lists of lists ``[j, a]``. The pairs ``[j, a]``
    encode the choice of the `a`-th element in the `j`-th list of the input
    ``LWBModp``, i.e., the `a`-th element in a particular basis modulo
    `(p,q^\text{elldash})` for the space of modular forms of level
    `\Gamma_0(N)` and weight `2(j+1)`. The list ``[[j_1, a_1], ...,[j_r, a_r]]``
    then encodes the product of the r modular forms associated to each
    ``[j_i, a_i]``; this has weight `k + (p-1)i` for some `0 \le i \le n`; here
    the `i` is such that this *list of lists* occurs in the `i`-th list of the
    output. The `i`-th list of the output thus encodes a choice of basis for the
    complementary space `W_i` which occurs in Step 2 of Algorithm 2 in [Lau2011]_.
    The idea is that one searches for this space `W_i` first modulo
    `(p,q^\text{elldash})` and then, having found the correct products of
    generating forms, one can reconstruct these spaces modulo
    `(p^\text{mdash},q^\text{elldashp})` using the output of this function.
    (This idea is based upon a suggestion of John Voight.)

    INPUT:

    - ``N`` -- positive integer at least 2 and not divisible by `p` (level)
    - ``p`` -- prime at least 5
    - ``k0`` -- integer in range 0 to `p-1`
    - ``n``, ``elldash`` -- positive integers
    - ``LWBModp`` -- list of lists of `q`-expansions over `GF(p)`
    - ``bound`` -- positive even integer (twice the length of the list ``LWBModp``)

    OUTPUT: list of list of list of lists

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import random_low_weight_bases, complementary_spaces_modp
        sage: LWB = random_low_weight_bases(2, 5, 2, 4, 6)
        sage: LWBModp = [[f.change_ring(Zmod(5)) for f in x] for x in LWB]
        sage: complementary_spaces_modp(2, 5, 0, 3, 4, LWBModp, 6) # random, indirect doctest
        [[[]], [[[0, 0], [0, 0]]], [[[0, 0], [2, 1]]], [[[0, 0], [0, 0], [0, 0], [2, 1]]]]
    """
    CompSpacesCode = []
    ell = dimension_modular_forms(N, k0 + n * (p - 1))
    TotalBasisModp = matrix(GF(p), ell, elldash)  # zero matrix

    for i in range(n + 1):
        NewBasisCodemi = random_new_basis_modp(N, p, k0 + i * (p - 1), LWBModp, TotalBasisModp, elldash, bound)
        # TotalBasisModp is passed by reference and updated in function
        CompSpacesCode.append(NewBasisCodemi)

    return CompSpacesCode

def complementary_spaces(N, p, k0, n, mdash, elldashp, elldash, modformsring, bound):
    r"""
    Return a list ``Ws``, each element in which is a list ``Wi`` of
    `q`-expansions modulo `(p^\text{mdash},q^\text{elldashp})`. The list ``Wi`` is
    a basis for a choice of complementary space in level `\Gamma_0(N)` and
    weight `k` to the image of weight `k - (p-1)` forms under multiplication by
    the Eisenstein series `E_{p-1}`.

    The lists ``Wi`` play the same role as `W_i` in Step 2 of Algorithm 2 in
    [Lau2011]_. (The parameters ``k0, n, mdash, elldash, elldashp = elldash * p`` are
    defined as in Step 1 of that algorithm when this function is used in
    :func:`hecke_series`.) However, the complementary spaces are computed in a
    different manner, combining a suggestion of David Loeffler with one of John
    Voight. That is, one builds these spaces recursively using random products
    of forms in low weight, first searching for suitable products modulo
    `(p,q^\text{elldash})`, and then later reconstructing only the required
    products to the full precision modulo `(p^\text{mdash},q^\text{elldashp})`. The
    forms in low weight are chosen from either bases of all forms up to weight
    ``bound`` or from a (tentative) generating set for the ring of all modular
    forms, according to whether ``modformsring`` is ``False`` or ``True``.

    INPUT:

    - ``N`` -- positive integer at least 2 and not divisible by p (level)
    - ``p`` -- prime at least 5
    - ``k0`` -- integer in range 0 to `p - 1`
    - ``n``, ``mdash``, ``elldashp``, ``elldash`` -- positive integers
    - ``modformsring`` -- boolean
    - ``bound`` -- positive (even) integer (ignored if ``modformsring`` is True)

    OUTPUT:

    - list of lists of `q`-expansions modulo
      `(p^\text{mdash},q^\text{elldashp})`.

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import complementary_spaces
        sage: complementary_spaces(2, 5, 0, 3, 2, 5, 4, True, 6) # random
        [[1],
        [1 + 23*q + 24*q^2 + 19*q^3 + 7*q^4 + O(q^5)],
        [1 + 21*q + 2*q^2 + 17*q^3 + 14*q^4 + O(q^5)],
        [1 + 19*q + 9*q^2 + 11*q^3 + 9*q^4 + O(q^5)]]
        sage: complementary_spaces(2, 5, 0, 3, 2, 5, 4, False, 6) # random
        [[1],
        [3 + 4*q + 2*q^2 + 12*q^3 + 11*q^4 + O(q^5)],
        [2 + 2*q + 14*q^2 + 19*q^3 + 18*q^4 + O(q^5)],
        [6 + 8*q + 10*q^2 + 23*q^3 + 4*q^4 + O(q^5)]]
    """
    if not modformsring:
        LWB = random_low_weight_bases(N, p, mdash, elldashp, bound)
    else:
        LWB, bound = low_weight_generators(N, p, mdash, elldashp)

    LWBModp = [ [ f.change_ring(GF(p)).truncate_powerseries(elldash) for f in x] for x in LWB]

    CompSpacesCode = complementary_spaces_modp(N, p, k0, n, elldash, LWBModp, bound)

    Ws = []
    Epm1 = eisenstein_series_qexp(p - 1, prec=elldashp, K=Zmod(p**mdash), normalization='constant')
    for i in range(n + 1):
        CompSpacesCodemi = CompSpacesCode[i]
        Wi = []
        for k in range(len(CompSpacesCodemi)):
            CompSpacesCodemik = CompSpacesCodemi[k]
            Wik = Epm1.parent()(1)
            for j in range(len(CompSpacesCodemik)):
                l = CompSpacesCodemik[j][0]
                index = CompSpacesCodemik[j][1]
                Wik = Wik*LWB[l][index]
            Wi.append(Wik)
        Ws.append(Wi)

    return Ws

# AUXILIARY CODE: KATZ EXPANSIONS

def higher_level_katz_exp(p, N, k0, m, mdash, elldash, elldashp, modformsring, bound):
    r"""
    Return a matrix `e` of size ``ell x elldashp`` over the integers modulo
    `p^\text{mdash}`, and the Eisenstein series `E_{p-1} = 1 + .\dots \bmod
    (p^\text{mdash},q^\text{elldashp})`. The matrix `e` contains the coefficients
    of the elements `e_{i,s}` in the Katz expansions basis in Step 3 of
    Algorithm 2 in [Lau2011]_ when one takes as input to that algorithm
    `p, N, m` and `k` and define ``k0``, ``mdash``, ``n``, ``elldash``,
    ``elldashp = ell * dashp`` as in Step 1.

    INPUT:

    - ``p`` -- prime at least 5
    - ``N`` -- positive integer at least 2 and not divisible by `p` (level)
    - ``k0`` -- integer in range 0 to `p-1`
    - ``m``, ``mdash, ``elldash``, ``elldashp`` -- positive integers
    - ``modformsring`` -- boolean
    - ``bound`` -- positive (even) integer

    OUTPUT: matrix and `q`-expansion

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import higher_level_katz_exp
        sage: e, Ep1 = higher_level_katz_exp(5, 2, 0, 1, 2, 4, 20, True, 6)
        sage: e
        [ 1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
        [ 0  1 18 23 19  6  9  9 17  7  3 17 12  8 22  8 11 19  1  5]
        [ 0  0  1 11 20 16  0  8  4  0 18 15 24  6 15 23  5 18  7 15]
        [ 0  0  0  1  4 16 23 13  6  5 23  5  2 16  4 18 10 23  5 15]
        sage: Ep1
        1 + 15*q + 10*q^2 + 20*q^3 + 20*q^4 + 15*q^5 + 5*q^6 + 10*q^7 +
        5*q^9 + 10*q^10 + 5*q^11 + 10*q^12 + 20*q^13 + 15*q^14 + 20*q^15 + 15*q^16 +
        10*q^17 + 20*q^18 + O(q^20)
    """
    ordr = 1 / (p + 1)
    S = Zmod(p ** mdash)
    Ep1 = eisenstein_series_qexp(p - 1, prec=elldashp, K=S, normalization='constant')

    n = floor(((p + 1) / (p - 1)) * (m + 1))
    Wjs = complementary_spaces(N, p, k0, n, mdash, elldashp, elldash, modformsring, bound)

    Basis = []
    for j in range(n + 1):
        Wj = Wjs[j]
        dimj = len(Wj)
        Ep1minusj = Ep1 ** (-j)
        for i in range(dimj):
            wji = Wj[i]
            b = p ** floor(ordr * j) * wji * Ep1minusj
            Basis.append(b)

    # extract basis as a matrix

    ell = len(Basis)
    M = matrix(S, ell, elldashp)
    for i in range(ell):
        for j in range(elldashp):
            M[i, j] = Basis[i][j]

    ech_form(M, p) # put it into echelon form

    return M, Ep1

def compute_elldash(p, N, k0, n):
    r"""
    Return the "Sturm bound" for the space of modular forms of level
    `\Gamma_0(N)` and weight `k_0 + n(p-1)`.

    .. SEEALSO::

        :meth:`~sage.modular.modform.space.ModularFormsSpace.sturm_bound`

    INPUT:

    - ``p`` -- prime
    - ``N`` -- positive integer (level)
    - ``k0``, ``n`` -- nonnegative integers not both zero

    OUTPUT: positive integer

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import compute_elldash
        sage: compute_elldash(11, 5, 4, 10)
        53
    """

    return ModularForms(N, k0 + n * (p - 1)).sturm_bound()

# *** DEGREE BOUND ON HECKE SERIES ***

def hecke_series_degree_bound(p, N, k, m):
    r"""
    Return the ``Wan bound`` on the degree of the characteristic series of the
    Atkin operator on `p`-adic overconvergent modular forms of level
    `\Gamma_0(N)` and weight `k` when reduced modulo `p^m`.

    This bound depends only upon `p, k \pmod{p-1}`, and `N`. It uses Lemma 3.1 in
    [Wan1998]_.

    INPUT:

    - ``p`` -- prime at least 5
    - ``N`` -- positive integer not divisible by `p`
    - ``k`` -- even integer
    - ``m`` -- positive integer

    OUTPUT: nonnegative integer

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import hecke_series_degree_bound
        sage: hecke_series_degree_bound(13,11,100,5)
        39
    """
    k0 = k % (p - 1)
    ds = [dimension_modular_forms(N, k0)]
    ms = [ds[0]]
    sum = 0
    u = 1

    ord = 0
    while ord < m:
        ds.append(dimension_modular_forms(N, k0 + u * (p - 1)))
        ms.append(ds[u] - ds[u - 1])
        sum = sum + u * ms[u]
        ord = floor(((p - 1) / (p + 1)) * sum - ds[u])
        u = u + 1

    return (ds[u - 1] - 1)

# *** MAIN FUNCTION FOR LEVEL > 1 ***

# Returns matrix A modulo p^m from Step 6 of Algorithm 2.

def higher_level_UpGj(p, N, klist, m, modformsring, bound, extra_data=False):
    r"""
    Return a list ``[A_k]`` of square matrices over ``IntegerRing(p^m)``
    parameterised by the weights `k` in ``klist``.

    The matrix `A_k` is the finite square matrix which occurs on input
    `p, k, N` and `m` in Step 6 of Algorithm 2 in [Lau2011]_.

    Notational change from paper: In Step 1 following Wan we defined
    `j` by `k = k_0 + j(p-1)` with `0 \le k_0 < p-1`. Here we replace `j` by
    ``kdiv`` so that we may use `j` as a column index for matrices.)

    INPUT:

    - ``p`` -- prime at least 5
    - ``N`` -- integer at least 2 and not divisible by `p` (level)
    - ``klist`` -- list of integers all congruent modulo `(p-1)` (the weights)
    - ``m`` -- positive integer
    - ``modformsring`` -- boolean
    - ``bound`` -- (even) positive integer
    - ``extra_data`` -- boolean (default: ``False``)

    OUTPUT:

    - list of square matrices. If ``extra_data`` is ``True``, return also
      extra intermediate data, namely the matrix `E` in [Lau2011]_ and
      the integers ``elldash`` and ``mdash``.

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import higher_level_UpGj
        sage: A = Matrix([
        ....:     [1,  0,  0,  0,  0,  0],
        ....:     [0,  1,  0,  0,  0,  0],
        ....:     [0,  7,  0,  0,  0,  0],
        ....:     [0,  5, 10, 20,  0,  0],
        ....:     [0,  7, 20,  0, 20,  0],
        ....:     [0,  1, 24,  0, 20,  0]])
        sage: B = Matrix([
        ....:     [1,  0,  0,  0,  0,  0],
        ....:     [0,  1,  0,  0,  0,  0],
        ....:     [0,  7,  0,  0,  0,  0],
        ....:     [0, 19,  0, 20,  0,  0],
        ....:     [0,  7, 20,  0, 20,  0],
        ....:     [0,  1, 24,  0, 20,  0]])
        sage: C = higher_level_UpGj(5, 3, [4], 2, True, 6)
        sage: len(C)
        1
        sage: C[0] in (A, B)
        True
        sage: len(higher_level_UpGj(5, 3, [4], 2, True, 6, extra_data=True))
        4
    """
    t = cputime()
    # Step 1

    k0 = klist[0] % (p - 1)
    n = floor(((p + 1) / (p - 1)) * (m + 1))
    elldash = compute_elldash(p, N, k0, n)
    elldashp = elldash * p
    mdash = m + ceil(n / (p + 1))

    verbose("done step 1", t)
    t = cputime()
    # Steps 2 and 3

    e, Ep1 = higher_level_katz_exp(p, N, k0, m, mdash, elldash, elldashp,
                                   modformsring, bound)
    ell = dimension(transpose(e)[0].parent())
    S = e[0, 0].parent()

    verbose("done steps 2+3", t)
    t = cputime()
    # Step 4

    R = Ep1.parent()
    G = compute_G(p, Ep1)
    Alist = []

    verbose("done step 4a", t)
    t = cputime()
    for k in klist:
        k = ZZ(k) # convert to sage integer
        kdiv = k // (p - 1)
        Gkdiv = G ** kdiv

        T = matrix(S, ell, elldash)
        for i in range(ell):
            ei = R(e[i].list())
            Gkdivei = Gkdiv * ei  # act by G^kdiv
            for j in range(elldash):
                T[i, j] = Gkdivei[p * j]

        verbose("done steps 4b and 5", t)
        t = cputime()

        # Step 6: solve T = AE using fact E is upper triangular.
        # Warning: assumes that T = AE (rather than pT = AE) has
        # a solution over Z/(p^mdash). This has always been the case in
        # examples computed by the author, see Note 3.1.

        A = matrix(S, ell, ell)
        verbose("solving a square matrix problem of dimension %s" % ell)
        verbose("elldash is %s" % elldash)

        for i in range(ell):
            Ti = T[i]
            for j in range(ell):
                ej = Ti.parent()([e[j][l] for l in range(elldash)])
                ejleadpos = ej.nonzero_positions()[0]
                lj = ZZ(ej[ejleadpos])
                A[i, j] = S(ZZ(Ti[j]) / lj)
                Ti = Ti - A[i, j] * ej

        Alist.append(MatrixSpace(Zmod(p ** m), ell, ell)(A))
        verbose("done step 6", t)

    if extra_data:
        return Alist, e, elldash, mdash
    else:
        return Alist


#  *** LEVEL 1 CODE ***

def compute_Wi(k, p, h, hj, E4, E6):
    r"""
    This function computes a list `W_i` of `q`-expansions, together with an
    auxiliary quantity `h^j` (see below) which is to be used on the next
    call of this function. (The precision is that of input `q`-expansions.)

    The list `W_i` is a certain subset of a basis of the modular forms of
    weight `k` and level 1. Suppose `(a, b)` is the pair of nonnegative
    integers with `4a + 6b = k` and `a` minimal among such pairs. Then this
    space has a basis given by

    .. MATH::

        \{ \Delta^j E_6^{b - 2j} E_4^a : 0 \le j < d\}

    where `d` is the dimension.

    What this function returns is the subset of the above basis corresponding
    to `e \le j < d` where `e` is the dimension of the space of modular forms
    of weight `k - (p-1)`. This set is a basis for the complement of the image
    of the weight `k - (p-1)` forms under multiplication by `E_{p-1}`.

    This function is used repeatedly in the construction of the Katz expansion
    basis. Hence considerable care is taken to reuse steps in the computation
    wherever possible: we keep track of powers of the form `h = \Delta /
    E_6^2`.

    INPUT:

    - ``k`` -- nonnegative integer
    - ``p`` -- prime at least 5
    - ``h`` -- `q`-expansion of `h` (to some finite precision)
    - ``hj`` -- `q`-expansion of `h^j` where `j` is the dimension of the space of
      modular forms of level 1 and weight `k - (p-1)` (to same finite
      precision)
    - ``E4`` -- `q`-expansion of `E_4` (to same finite precision)
    - ``E6`` -- `q`-expansion of `E_6` (to same finite precision)

    The Eisenstein series `q`-expansions should be normalized to have constant
    term 1.

    OUTPUT:

    - list of `q`-expansions (to same finite precision), and `q`-expansion.

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import compute_Wi
        sage: p = 17
        sage: prec = 10
        sage: k = 24
        sage: S = Zmod(17^3)
        sage: E4 = eisenstein_series_qexp(4, prec, K=S, normalization='constant')
        sage: E6 = eisenstein_series_qexp(6, prec, K=S, normalization='constant')
        sage: h = delta_qexp(prec, K=S) / E6^2
        sage: from sage.modular.dims import dimension_modular_forms
        sage: j = dimension_modular_forms(1, k - (p - 1))
        sage: hj = h ** j
        sage: c = compute_Wi(k, p, h, hj, E4, E6); c
        ([q + 3881*q^2 + 4459*q^3 + 4665*q^4 + 2966*q^5 + 1902*q^6 + 1350*q^7 + 3836*q^8 + 1752*q^9 + O(q^10), q^2 + 4865*q^3 + 1080*q^4 + 4612*q^5 + 1343*q^6 + 1689*q^7 + 3876*q^8 + 1381*q^9 + O(q^10)], q^3 + 2952*q^4 + 1278*q^5 + 3225*q^6 + 1286*q^7 + 589*q^8 + 122*q^9 + O(q^10))
        sage: c == ([delta_qexp(10) * E6^2, delta_qexp(10)^2], h**3)
        True
    """

    # Define a and b
    a = k % 3
    b = (k // 2) % 2

    # Compute dimensions required for Miller basis
    d = dimension_modular_forms(1, k) - 1
    e = dimension_modular_forms(1, k - (p - 1)) - 1

    # This next line is a bit of a bottleneck, particularly when m is large but
    # p is small. It would be good to reuse values calculated on the previous
    # call here somehow.
    r = E6 ** (2 * d + b) * E4 ** a

    prec = E4.prec() # everything gets truncated to this precision

    # Construct basis for Wi
    Wi = []
    for j in range(e + 1,d + 1):
        # compute aj = delta^j*E6^(2*(d-j) + b)*E4^a
        verbose("k = %s, computing Delta^%s E6^%s E4^%s" % (k, j, 2 * (d - j) + b, a), level=2)
        aj = (hj * r).truncate_powerseries(prec)
        hj = (hj * h).truncate_powerseries(prec)
        Wi.append(aj)

    return Wi, hj

def katz_expansions(k0, p, ellp, mdash, n):
    r"""
    Return a list `e` of `q`-expansions, and the Eisenstein series `E_{p-1} = 1 +
    \dots`, all modulo `(p^\text{mdash},q^\text{ellp})`. The list `e` contains
    the elements `e_{i,s}` in the Katz expansions basis in Step 3 of Algorithm
    1 in [Lau2011]_ when one takes as input to that algorithm `p,m` and `k` and define
    ``k0``, ``mdash``, ``n``, ``ellp = ell * p`` as in Step 1.

    INPUT:

    - ``k0`` -- integer in range 0 to `p - 1`
    - ``p`` -- prime at least 5
    - ``ellp``, ``mdash``, ``n`` -- positive integers

    OUTPUT:

    - list of `q`-expansions and the Eisenstein series `E_{p-1}` modulo
      `(p^\text{mdash},q^\text{ellp})`.

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import katz_expansions
        sage: katz_expansions(0, 5, 10, 3, 4)
        ([1 + O(q^10), q + 6*q^2 + 27*q^3 + 98*q^4 + 65*q^5 + 37*q^6 + 81*q^7 + 85*q^8 + 62*q^9 + O(q^10)],
        1 + 115*q + 35*q^2 + 95*q^3 + 20*q^4 + 115*q^5 + 105*q^6 + 60*q^7 + 25*q^8 + 55*q^9 + O(q^10))
    """
    S = Zmod(p ** mdash)

    Ep1 = eisenstein_series_qexp(p - 1, ellp, K=S, normalization='constant')
    E4 = eisenstein_series_qexp(4, ellp, K=S, normalization='constant')
    E6 = eisenstein_series_qexp(6, ellp, K=S, normalization='constant')

    delta = delta_qexp(ellp, K=S)
    h = delta / E6 ** 2
    hj = delta.parent()(1)
    e = []

    # We compute negative powers of E_(p-1) successively (this saves a great
    # deal of time). The effect is that Ep1mi = Ep1 ** (-i).
    Ep1m1 = ~Ep1
    Ep1mi = 1
    for i in range(0,n + 1):
        Wi,hj = compute_Wi(k0 + i * (p - 1), p, h, hj, E4, E6)
        for bis in Wi:
            eis = p ** floor(i / (p + 1)) * Ep1mi * bis
            e.append(eis)
        Ep1mi = Ep1mi * Ep1m1

    return e, Ep1

# *** MAIN FUNCTION FOR LEVEL 1 ***

def level1_UpGj(p, klist, m, extra_data=False):
    r"""
    Return a list `[A_k]` of square matrices over ``IntegerRing(p^m)``
    parameterised by the weights `k` in ``klist``.

    The matrix `A_k` is the finite square matrix which occurs on input
    `p, k` and `m` in Step 6 of Algorithm 1 in [Lau2011]_.

    Notational change from paper: In Step 1 following Wan we defined
    `j` by `k = k_0 + j(p-1)` with `0 \le k_0 < p-1`. Here we replace `j` by
    ``kdiv`` so that we may use `j` as a column index for matrices.

    INPUT:

    - ``p`` -- prime at least 5
    - ``klist`` -- list of integers congruent modulo `(p-1)` (the weights)
    - ``m`` -- positive integer
    - ``extra_data`` -- boolean (default: ``False``)

    OUTPUT:

    - list of square matrices. If ``extra_data`` is ``True``, return also
      extra intermediate data, namely the matrix `E` in [Lau2011]_ and
      the integers ``elldash`` and ``mdash``.

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import level1_UpGj
        sage: level1_UpGj(7, [100], 5)
        [
        [    1   980  4802     0     0]
        [    0 13727 14406     0     0]
        [    0 13440  7203     0     0]
        [    0  1995  4802     0     0]
        [    0  9212 14406     0     0]
        ]
        sage: len(level1_UpGj(7, [100], 5, extra_data=True))
        4
    """
    # Step 1
    t = cputime()

    k0 = klist[0] % (p - 1)
    n = floor(((p + 1) / (p - 1)) * (m + 1))
    ell = dimension_modular_forms(1, k0 + n * (p - 1))
    ellp = ell * p
    mdash = m + ceil(n / (p + 1))

    verbose("done step 1", t)
    t = cputime()
    # Steps 2 and 3

    e, Ep1 = katz_expansions(k0, p, ellp, mdash, n)

    verbose("done steps 2+3", t)
    t = cputime()
    # Step 4

    G = compute_G(p, Ep1)
    Alist = []

    verbose("done step 4a", t)
    t = cputime()
    for k in klist:
        k = ZZ(k) # convert to sage integer
        kdiv = k // (p - 1)
        Gkdiv = G ** kdiv
        u = []
        for i in range(0, ell):
            ei = e[i]
            ui = Gkdiv * ei
            u.append(ui)

        verbose("done step 4b", t)
        t = cputime()
        # Step 5 and computation of T in Step 6

        S = e[0][0].parent()
        T = matrix(S, ell, ell)

        for i in range(0, ell):
            for j in range(0, ell):
                T[i, j] = u[i][p * j]

        verbose("done step 5", t)
        t = cputime()
        # Step 6: solve T = AE using fact E is upper triangular.
        # Warning: assumes that T = AE (rather than pT = AE) has
        # a solution over Z/(p^mdash). This has always been the case in
        # examples computed by the author, see Note 3.1.

        A = matrix(S, ell, ell)
        verbose("solving a square matrix problem of dimension %s" % ell, t)

        for i in range(0, ell):
            Ti = T[i]
            for j in range(0,ell):
                ej = Ti.parent()([e[j][l] for l in range(0, ell)])
                lj = ZZ(ej[j])
                A[i, j] = S(ZZ(Ti[j]) / lj)
                Ti = Ti - A[i, j] * ej

        Alist.append(MatrixSpace(Zmod(p ** m), ell, ell)(A))
        verbose("done step 6", t)

    if extra_data:
        return Alist, e, ell, mdash
    else:
        return Alist

# *** CODE FOR GENERAL LEVEL ***

def is_valid_weight_list(klist, p):
    r"""
    This function checks that ``klist`` is a nonempty list of integers all of
    which are congruent modulo `(p-1)`. Otherwise, it will raise a
    :exc:`ValueError`.

    INPUT:

    - ``klist`` -- list of integers
    - ``p`` -- prime

    EXAMPLES::

        sage: from sage.modular.overconvergent.hecke_series import is_valid_weight_list
        sage: is_valid_weight_list([10, 20, 30],11)
        sage: is_valid_weight_list([-3, 1], 5)
        sage: is_valid_weight_list([], 3)
        Traceback (most recent call last):
        ...
        ValueError: List of weights must be non-empty
        sage: is_valid_weight_list([-3, 2], 5)
        Traceback (most recent call last):
        ...
        ValueError: List of weights must be all congruent modulo p-1 = 4, but given list contains -3 and 2 which are not congruent
    """
    if len(klist) == 0:
        raise ValueError("List of weights must be non-empty")
    k0 = klist[0] % (p - 1)
    for i in range(1,len(klist)):
        if (klist[i] % (p-1)) != k0:
            raise ValueError("List of weights must be all congruent modulo p-1 = %s, but given list contains %s and %s which are not congruent" % (p - 1, klist[0], klist[i]))

def hecke_series(p, N, klist, m, modformsring=False, weightbound=6):
    r"""
    Return the characteristic series modulo `p^m` of the Atkin operator `U_p`
    acting upon the space of `p`-adic overconvergent modular forms of level
    `\Gamma_0(N)` and weight ``klist``.

    The input ``klist`` may also be a list of weights congruent modulo `(p-1)`,
    in which case the output is the corresponding list of characteristic series
    for each `k` in ``klist``; this is faster than performing the computation
    separately for each `k`, since intermediate steps in the computation may be
    reused.

    If ``modformsring`` is ``True``, then for `N > 1` the algorithm computes at one
    step ``ModularFormsRing(N).generators()``. This will often be faster but
    the algorithm will default to ``modformsring=False`` if the generators
    found are not `p`-adically integral. Note that ``modformsring`` is ignored
    for `N = 1` and the ring structure of modular forms is *always* used in
    this case.

    When ``modformsring`` is ``False`` and `N > 1`, ``weightbound`` is a bound set on
    the weight of generators for a certain subspace of modular forms. The
    algorithm will often be faster if ``weightbound=4``, but it may fail to
    terminate for certain exceptional small values of `N`, when this bound is
    too small.

    The algorithm is based upon that described in [Lau2011]_.

    INPUT:

    - ``p`` -- a prime greater than or equal to 5
    - ``N`` -- positive integer not divisible by `p`
    - ``klist`` -- either a list of integers congruent modulo `(p-1)`, or a single integer
    - ``m`` -- positive integer
    - ``modformsring`` -- boolean (default: ``False``); ignored if `N = 1`
    - ``weightbound`` -- a positive even integer (default: 6). Ignored
      if `N = 1` or ``modformsring`` is ``True``

    OUTPUT: either a list of polynomials or a single polynomial over the integers modulo `p^m`

    EXAMPLES::

        sage: hecke_series(5, 7, 10000, 5, modformsring=True) # long time (3.4s)
        250*x^6 + 1825*x^5 + 2500*x^4 + 2184*x^3 + 1458*x^2 + 1157*x + 1
        sage: hecke_series(7, 3, 10000, 3, weightbound=4)
        196*x^4 + 294*x^3 + 197*x^2 + 341*x + 1
        sage: hecke_series(19, 1, [10000, 10018], 5)
        [1694173*x^4 + 2442526*x^3 + 1367943*x^2 + 1923654*x + 1,
        130321*x^4 + 958816*x^3 + 2278233*x^2 + 1584827*x + 1]

    Check that silly weights are handled correctly::

        sage: hecke_series(5, 7, [2, 3], 5)
        Traceback (most recent call last):
        ...
        ValueError: List of weights must be all congruent modulo p-1 = 4, but given list contains 2 and 3 which are not congruent
        sage: hecke_series(5, 7, [3], 5)
        [1]
        sage: hecke_series(5, 7, 3, 5)
        1
    """
    # convert to sage integers
    p = ZZ(p)
    N = ZZ(N)
    m = ZZ(m)
    weightbound = ZZ(weightbound)

    oneweight = False
    # convert single weight to list
    if ((isinstance(klist, int)) or (isinstance(klist, Integer))):
        klist = [klist]
        oneweight = True # input is single weight

    # algorithm may finish with false output unless:
    is_valid_weight_list(klist, p)
    if not p.is_prime():
        raise ValueError("p (=%s) is not prime" % p)
    if p < 5:
        raise ValueError("p = 2 and p = 3 not supported")
    if not N % p:
        raise ValueError("Level (=%s) should be prime to p (=%s)" % (N, p))

    # return all 1 list for odd weights
    if klist[0] % 2 == 1:
        if oneweight:
            return 1
        else:
            return [1 for i in range(len(klist))]

    if N == 1:
        Alist = level1_UpGj(p, klist, m)
    else:
        Alist = higher_level_UpGj(p, N, klist, m, modformsring, weightbound)

    Plist = []
    for A in Alist:
        P = charpoly(A).reverse()
        Plist.append(P)

    if oneweight:
        return Plist[0]
    else:
        return Plist
