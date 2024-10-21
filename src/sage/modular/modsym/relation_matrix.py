# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.flint
"""
Relation matrices for ambient modular symbols spaces

This file contains functions that are used by the various ambient modular
symbols classes to compute presentations of spaces in terms of generators and
relations, using the standard methods based on Manin symbols.
"""
# ****************************************************************************
#       Sage: Open Source Mathematical Software
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.matrix_space import MatrixSpace
from sage.misc.search import search
from sage.misc.verbose import verbose
from sage.modular.modsym.manin_symbol_list import ManinSymbolList
from sage.rings.rational_field import RationalField
from sage.categories.rings import Rings


SPARSE = True

# S = [0,-1; 1,0]
# T = [0,-1; 1,-1],
# T^2 = [-1, 1, -1, 0]
# I = [-1,0; 0,1]


######################################################################
# The following four functions are used to compute the quotient
# modulo the S, I, and T relations more efficiently that the generic
# code in the relation_matrix file:
#    modS_relations -- compute the S relations.
#    modI_quotient --  compute the I relations.
#    T_relation_matrix -- matrix whose echelon form gives
#                                 the quotient by 3-term T relations.
#    gens_to_basis_matrix -- compute echelon form of 3-term
#                                    relation matrix, and read off each
#                                    generator in terms of basis.
# These four functions are orchestrated in the function
#    compute_presentation
# which is defined below.  See the comment at the beginning
# of that function for an overall description of the algorithm.
######################################################################


def modS_relations(syms):
    r"""
    Compute quotient of Manin symbols by the S relations.

    Here S is the 2x2 matrix [0, -1; 1, 0].

    INPUT:

    - ``syms`` -- :class:`ManinSymbolList`

    OUTPUT:

    - ``rels`` -- set of pairs of pairs (j, s), where if
      mod[i] = (j,s), then x_i = s\*x_j (mod S relations)

    EXAMPLES::

        sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma0
        sage: from sage.modular.modsym.relation_matrix import modS_relations

    ::

        sage: syms = ManinSymbolList_gamma0(2, 4); syms
        Manin Symbol List of weight 4 for Gamma0(2)
        sage: modS_relations(syms)
        {((0, 1), (7, 1)),
         ((1, 1), (6, 1)),
         ((2, 1), (8, 1)),
         ((3, -1), (4, 1)),
         ((3, 1), (4, -1)),
         ((5, -1), (5, 1))}

    ::

        sage: syms = ManinSymbolList_gamma0(7, 2); syms
        Manin Symbol List of weight 2 for Gamma0(7)
        sage: modS_relations(syms)
        {((0, 1), (1, 1)), ((2, 1), (7, 1)), ((3, 1), (4, 1)), ((5, 1), (6, 1))}

    Next we do an example with Gamma1::

        sage: from sage.modular.modsym.manin_symbol_list import ManinSymbolList_gamma1
        sage: syms = ManinSymbolList_gamma1(3,2); syms
        Manin Symbol List of weight 2 for Gamma1(3)
        sage: modS_relations(syms)
        {((0, 1), (2, 1)),
         ((0, 1), (5, 1)),
         ((1, 1), (2, 1)),
         ((1, 1), (5, 1)),
         ((3, 1), (4, 1)),
         ((3, 1), (6, 1)),
         ((4, 1), (7, 1)),
         ((6, 1), (7, 1))}
    """
    if not isinstance(syms, ManinSymbolList):
        raise TypeError("syms must be a ManinSymbolList")
    tm = verbose()
    # We will fill in this set with the relations x_i + s*x_j = 0,
    # where the notation is as in _sparse_2term_quotient.
    rels = set()
    for i in range(len(syms)):
        j, s = syms.apply_S(i)
        assert j != -1
        if i < j:
            rels.add(((i, 1), (j, s)))
        else:
            rels.add(((j, s), (i, 1)))
    verbose("finished creating S relations", tm)
    return rels


def modI_relations(syms, sign):
    r"""
    Compute quotient of Manin symbols by the I relations.

    INPUT:

    - ``syms`` -- :class:`ManinSymbolList`

    - ``sign`` -- integer (either -1, 0, or 1)

    OUTPUT:

    - ``rels`` -- set of pairs of pairs (j, s), where if
      mod[i] = (j,s), then x_i = s\*x_j (mod S relations)

    EXAMPLES::

        sage: L = sage.modular.modsym.manin_symbol_list.ManinSymbolList_gamma1(4, 3)
        sage: sage.modular.modsym.relation_matrix.modI_relations(L, 1)
        {((0, 1), (0, -1)),
         ((1, 1), (1, -1)),
         ((2, 1), (8, -1)),
         ((3, 1), (9, -1)),
         ((4, 1), (10, -1)),
         ((5, 1), (11, -1)),
         ((6, 1), (6, -1)),
         ((7, 1), (7, -1)),
         ((8, 1), (2, -1)),
         ((9, 1), (3, -1)),
         ((10, 1), (4, -1)),
         ((11, 1), (5, -1)),
         ((12, 1), (12, 1)),
         ((13, 1), (13, 1)),
         ((14, 1), (20, 1)),
         ((15, 1), (21, 1)),
         ((16, 1), (22, 1)),
         ((17, 1), (23, 1)),
         ((18, 1), (18, 1)),
         ((19, 1), (19, 1)),
         ((20, 1), (14, 1)),
         ((21, 1), (15, 1)),
         ((22, 1), (16, 1)),
         ((23, 1), (17, 1))}

    .. warning::

       We quotient by the involution eta((u,v)) = (-u,v), which has
       the opposite sign as the involution in Merel's Springer LNM
       1585 paper! Thus our +1 eigenspace is his -1 eigenspace,
       etc. We do this for consistency with MAGMA.
    """
    tm = verbose()
    # We will fill in this set with the relations x_i - sign*s*x_j = 0,
    # where the notation is as in _sparse_2term_quotient.
    rels = set()
    for i in range(len(syms)):
        j, s = syms.apply_I(i)
        assert j != -1
        rels.add(((i, 1), (j, -sign * s)))
    verbose("finished creating I relations", tm)
    return rels


def T_relation_matrix_wtk_g0(syms, mod, field, sparse):
    r"""
    Compute a matrix whose echelon form gives the quotient by 3-term T
    relations. Despite the name, this is used for all modular symbols spaces
    (including those with character and those for `\Gamma_1` and `\Gamma_H`
    groups), not just `\Gamma_0`.

    INPUT:

    - ``syms`` -- :class:`ManinSymbolList`

    - ``mod`` -- list that gives quotient modulo some two-term relations, i.e.,
      the S relations, and if sign is nonzero, the I relations

    - ``field`` -- ``base_ring``

    - ``sparse`` -- boolean; whether to use sparse rather than dense
      linear algebra

    OUTPUT: a sparse matrix whose rows correspond to the reduction of
    the `T` relations modulo the `S` and `I` relations.

    EXAMPLES::

        sage: from sage.modular.modsym.relation_matrix import sparse_2term_quotient, T_relation_matrix_wtk_g0, modS_relations
        sage: L = sage.modular.modsym.manin_symbol_list.ManinSymbolList_gamma_h(GammaH(36, [17,19]), 2)
        sage: modS = sparse_2term_quotient(modS_relations(L), 216, QQ)
        sage: T_relation_matrix_wtk_g0(L, modS, QQ, False)
        72 x 216 dense matrix over Rational Field (use the '.str()' method to see the entries)
        sage: T_relation_matrix_wtk_g0(L, modS, GF(17), True)
        72 x 216 sparse matrix over Finite Field of size 17 (use the '.str()' method to see the entries)
    """
    tm = verbose()
    row = 0
    entries = {}
    already_seen = set()
    w = syms.weight()
    for i in range(len(syms)):
        if i in already_seen:
            continue
        iT_plus_iTT = syms.apply_T(i) + syms.apply_TT(i)
        j0, s0 = mod[i]
        v = {j0: s0}
        for j, s in iT_plus_iTT:
            if w == 2:
                already_seen.add(j)
            j0, s0 = mod[j]
            s0 = s * s0
            if j0 in v:
                v[j0] += s0
            else:
                v[j0] = s0
        for j0 in v:
            entries[(row, j0)] = v[j0]
        row += 1

    MAT = MatrixSpace(field, row, len(syms), sparse=True)
    R = MAT(entries)
    if not sparse:
        R = R.dense_matrix()
    verbose("finished (number of rows=%s)" % row, tm)
    return R


def gens_to_basis_matrix(syms, relation_matrix, mod, field, sparse):
    """
    Compute echelon form of 3-term relation matrix, and read off each
    generator in terms of basis.

    INPUT:

    - ``syms`` -- :class:`ManinSymbolList`

    - ``relation_matrix`` -- as output by
      ``__compute_T_relation_matrix(self, mod)``

    - ``mod`` -- quotient of modular symbols modulo the
      2-term S (and possibly I) relations

    - ``field`` -- base field

    - ``sparse`` -- boolean; whether or not matrix should be sparse

    OUTPUT:

    ``matrix`` -- a matrix whose `i`-th row expresses the Manin symbol
    generators in terms of a basis of Manin symbols (modulo the S, (possibly
    I,) and T rels). Note that the entries of the matrix need not be integers.

    - ``list`` -- integers `i`, such that the Manin symbols `x_i` are a basis

    EXAMPLES::

        sage: from sage.modular.modsym.relation_matrix import sparse_2term_quotient, T_relation_matrix_wtk_g0, gens_to_basis_matrix, modS_relations
        sage: L = sage.modular.modsym.manin_symbol_list.ManinSymbolList_gamma1(4, 3)
        sage: modS = sparse_2term_quotient(modS_relations(L), 24, GF(3))
        sage: gens_to_basis_matrix(L, T_relation_matrix_wtk_g0(L, modS, GF(3), 24), modS, GF(3), True)
        (24 x 2 sparse matrix over Finite Field of size 3, [13, 23])
    """
    from sage.structure.element import Matrix
    if not isinstance(relation_matrix, Matrix):
        raise TypeError("relation_matrix must be a matrix")
    if not isinstance(mod, list):
        raise TypeError("mod must be a list")

    verbose(str(relation_matrix.parent()))

    try:
        h = relation_matrix.height()
    except AttributeError:
        h = 9999999
    tm = verbose("putting relation matrix in echelon form (height = %s)" % h)
    if h < 10:
        A = relation_matrix.echelon_form(algorithm='multimodular',
                                         height_guess=1)
    else:
        A = relation_matrix.echelon_form()
    A.set_immutable()
    tm = verbose('finished echelon', tm)

    tm = verbose("Now creating gens --> basis mapping")

    basis_set = set(A.nonpivots())
    pivots = A.pivots()

    basis_mod2 = {j for j, c in mod if c != 0}

    basis_set = basis_set.intersection(basis_mod2)
    basis = sorted(basis_set)

    ONE = field(1)

    verbose("done doing setup", tm)

    tm = verbose("now forming quotient matrix")
    M = MatrixSpace(field, len(syms), len(basis), sparse=sparse)

    B = M(0)
    cols_index = {basis[i]: i for i in range(len(basis))}

    for i in basis_mod2:
        t, l = search(basis, i)
        if t:
            B[i, l] = ONE
        else:
            _, r = search(pivots, i)    # so pivots[r] = i
            # Set row i to -(row r of A), but where we only take
            # the non-pivot columns of A:
            B._set_row_to_negative_of_row_of_A_using_subset_of_columns(i, A, r, basis, cols_index)

    verbose("done making quotient matrix", tm)

    # The following is very fast (over Q at least).
    tm = verbose('now filling in the rest of the matrix')
    k = 0
    for i in range(len(mod)):
        j, s = mod[i]
        if j != i and s != 0:   # ignored in the above matrix
            k += 1
            B.set_row_to_multiple_of_row(i, j, s)
    verbose("set %s rows" % k)
    tm = verbose("time to fill in rest of matrix", tm)

    return B, basis


def compute_presentation(syms, sign, field, sparse=None):
    r"""
    Compute the presentation for self, as a quotient of Manin symbols
    modulo relations.

    INPUT:

    - ``syms`` -- :class:`ManinSymbolList`

    - ``sign`` -- integer (-1, 0, 1)

    - ``field`` -- a field

    OUTPUT:

    - sparse matrix whose rows give each generator
      in terms of a basis for the quotient

    - list of integers that give the basis for the
      quotient

    - mod: list where mod[i]=(j,s) means that x_i
      = s\*x_j modulo the 2-term S (and possibly I) relations.


    ALGORITHM:

    #. Let `S = [0,-1; 1,0], T = [0,-1; 1,-1]`, and
       `I = [-1,0; 0,1]`.

    #. Let `x_0,\ldots, x_{n-1}` by a list of all
       non-equivalent Manin symbols.

    #. Form quotient by 2-term S and (possibly) I relations.

    #. Create a sparse matrix `A` with `m` columns,
       whose rows encode the relations

       .. MATH::

                          [x_i] + [x_i T] + [x_i T^2] = 0.


       There are about n such rows. The number of nonzero entries per row
       is at most 3\*(k-1). Note that we must include rows for *all* i,
       since even if `[x_i] = [x_j]`, it need not be the case
       that `[x_i T] = [x_j T]`, since `S` and
       `T` do not commute. However, in many cases we have an a
       priori formula for the dimension of the quotient by all these
       relations, so we can omit many relations and just check that there
       are enough at the end--if there aren't, we add in more.

    #. Compute the reduced row echelon form of `A` using sparse
       Gaussian elimination.

    #. Use what we've done above to read off a sparse matrix R that
       uniquely expresses each of the n Manin symbols in terms of a subset
       of Manin symbols, modulo the relations. This subset of Manin
       symbols is a basis for the quotient by the relations.

    EXAMPLES::

        sage: L = sage.modular.modsym.manin_symbol_list.ManinSymbolList_gamma0(8,2)
        sage: sage.modular.modsym.relation_matrix.compute_presentation(L, 1, GF(9,'a'), True)
        (
        [2 0 0]
        [1 0 0]
        [0 0 0]
        [0 2 0]
        [0 0 0]
        [0 0 2]
        [0 0 0]
        [0 2 0]
        [0 0 0]
        [0 1 0]
        [0 1 0]
        [0 0 1], [1, 9, 11], [(1, 2), (1, 1), (0, 0), (9, 2), (0, 0), (11, 2), (0, 0), (9, 2), (0, 0), (9, 1), (9, 1), (11, 1)]
        )
    """
    if sparse is None:
        if syms.weight() >= 6:
            sparse = False
        else:
            sparse = True
    R, mod = relation_matrix_wtk_g0(syms, sign, field, sparse)
    B, basis = gens_to_basis_matrix(syms, R, mod, field, sparse)
    return B, basis, mod


def relation_matrix_wtk_g0(syms, sign, field, sparse):
    r"""
    Compute the matrix of relations. Despite the name, this is used for all
    spaces (not just for Gamma0). For a description of the algorithm, see the
    docstring for ``compute_presentation``.

    INPUT:

    - ``syms`` -- :class:`ManinSymbolList`

    - ``sign`` -- integer (0, 1 or -1)

    - ``field`` -- the base field (non-field base rings not supported at present)

    - ``sparse`` -- boolean; whether to use sparse arithmetic

    Note that ManinSymbolList objects already have a specific weight, so there
    is no need for an extra ``weight`` parameter.

    OUTPUT: a pair (R, mod) where

    - R is a matrix as output by ``T_relation_matrix_wtk_g0``

    - mod is a set of 2-term relations as output by ``sparse_2term_quotient``

    EXAMPLES::

        sage: L = sage.modular.modsym.manin_symbol_list.ManinSymbolList_gamma0(8,2)
        sage: A = sage.modular.modsym.relation_matrix.relation_matrix_wtk_g0(L, 0, GF(2), True); A
        (
        [0 0 0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 1 1 1 0]
        [0 0 0 0 0 0 1 0 0 1 1 0]
        [0 0 0 0 0 0 1 0 0 0 0 0], [(1, 1), (1, 1), (8, 1), (10, 1), (6, 1), (11, 1), (6, 1), (9, 1), (8, 1), (9, 1), (10, 1), (11, 1)]
        )
        sage: A[0].is_sparse()
        True
    """
    rels = modS_relations(syms)
    if sign != 0:
        # Let rels = rels union I relations.
        rels.update(modI_relations(syms, sign))

    rels = sorted(rels)
    # required for stability of doctests with python3

    if syms._apply_S_only_0pm1() and isinstance(field, RationalField):
        from . import relation_matrix_pyx
        mod = relation_matrix_pyx.sparse_2term_quotient_only_pm1(rels, len(syms))
    else:
        mod = sparse_2term_quotient(rels, len(syms), field)

    R = T_relation_matrix_wtk_g0(syms, mod, field, sparse)
    return R, mod


def sparse_2term_quotient(rels, n, F):
    r"""
    Perform Sparse Gauss elimination on a matrix all of whose columns
    have at most 2 nonzero entries. We use an obvious algorithm, which
    runs fast enough. (Typically making the list of relations takes
    more time than computing this quotient.) This algorithm is more
    subtle than just "identify symbols in pairs", since complicated
    relations can cause generators to surprisingly equal 0.

    INPUT:

    - ``rels`` -- iterable made of pairs ((i,s), (j,t)). The pair
      represents the relation `s x_i + t x_j = 0`, where the `i, j` must
      be Python int's.

    - ``n`` -- integer, the `x_i` are `x_0, \ldots, x_{n-1}`

    - ``F`` -- base field

    OUTPUT:

    ``mod`` -- list such that ``mod[i] = (j,s)``, which means that `x_i` is
    equivalent to `s x_j`, where the `x_j` are a basis for the quotient.

    EXAMPLES: We quotient out by the relations

    .. MATH::

        3*x0 - x1 = 0,\qquad  x1 + x3 = 0,\qquad   x2 + x3 = 0,\qquad  x4 - x5 = 0

    to get::

        sage: rels = [((int(0),3), (int(1),-1)), ((int(1),1), (int(3),1)), ((int(2),1),(int(3),1)), ((int(4),1),(int(5),-1))]
        sage: n = 6
        sage: from sage.modular.modsym.relation_matrix import sparse_2term_quotient
        sage: sparse_2term_quotient(rels, n, QQ)
        [(3, -1/3), (3, -1), (3, -1), (3, 1), (5, 1), (5, 1)]
    """
    n = int(n)
    if F not in Rings():
        raise TypeError("F must be a ring")

    tm = verbose("Starting sparse 2-term quotient...")
    free = list(range(n))
    ONE = F.one()
    ZERO = F.zero()
    coef = [ONE for i in range(n)]
    related_to_me = [[] for i in range(n)]
    for v0, v1 in sorted(rels):
        c0 = coef[v0[0]] * F(v0[1])
        c1 = coef[v1[0]] * F(v1[1])

        # Mod out by the following relation:
        #
        #    c0*free[v0[0]] + c1*free[v1[0]] = 0.
        #
        die = None
        if c0 == ZERO and c1 == ZERO:
            pass
        elif c0 == ZERO and c1 != ZERO:  # free[v1[0]] --> 0
            die = free[v1[0]]
        elif c1 == ZERO and c0 != ZERO:
            die = free[v0[0]]
        elif free[v0[0]] == free[v1[0]]:
            if c0 + c1 != 0:
                # all xi equal to free[v0[0]] must now equal to zero.
                die = free[v0[0]]
        else:  # x1 = -c1/c0 * x2.
            x = free[v0[0]]
            free[x] = free[v1[0]]
            coef[x] = -c1 / c0
            for i in related_to_me[x]:
                free[i] = free[x]
                coef[i] *= coef[x]
                related_to_me[free[v1[0]]].append(i)
            related_to_me[free[v1[0]]].append(x)
        if die is not None:
            for i in related_to_me[die]:
                free[i] = 0
                coef[i] = ZERO
            free[die] = 0
            coef[die] = ZERO

    mod = [(free[i], coef[i]) for i in range(len(free))]
    verbose("finished", tm)
    return mod
