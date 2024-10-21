# sage_setup: distribution = sagemath-gap
# sage.doctest: needs sage.libs.pari
r"""
Canonical forms and automorphisms for linear codes over finite fields

We implemented the algorithm described in [Feu2009]_ which computes, a unique
code (canonical form) in the equivalence class of a given
linear code `C \leq \GF{q}^n`. Furthermore, this algorithm will return the
automorphism group of `C`, too. You will find more details about the algorithm
in the documentation of the class
:class:`~sage.coding.codecan.autgroup_can_label.LinearCodeAutGroupCanLabel`.

The equivalence of codes is modeled as a group action by the group
`G = {\GF{q}^*}^n \rtimes (Aut(\GF{q}) \times S_n)` on the set
of subspaces of `\GF{q}^n` . The group `G` will be called the
semimonomial group of degree `n`.

The algorithm is started by initializing the class
:class:`~sage.coding.codecan.autgroup_can_label.LinearCodeAutGroupCanLabel`.
When the object gets available, all computations are already finished and
you can access the relevant data using the member functions:

* :meth:`~sage.coding.codecan.autgroup_can_label.LinearCodeAutGroupCanLabel.get_canonical_form`

* :meth:`~sage.coding.codecan.autgroup_can_label.LinearCodeAutGroupCanLabel.get_transporter`

* :meth:`~sage.coding.codecan.autgroup_can_label.LinearCodeAutGroupCanLabel.get_autom_gens`

People do also use some weaker notions of equivalence, namely
**permutational** equivalence and monomial equivalence (**linear** isometries).
These can be seen as the subgroups `S_n` and `{\GF{q}^*}^n \rtimes S_n` of `G`.
If you are interested in one of these notions, you can just pass
the optional parameter ``algorithm_type``.

A second optional parameter ``P`` allows you to restrict the
group of permutations `S_n` to a subgroup which respects the coloring given
by ``P``.

AUTHORS:

- Thomas Feulner (2012-11-15): initial version

EXAMPLES::

    sage: from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
    sage: C = codes.HammingCode(GF(3), 3).dual_code()
    sage: P = LinearCodeAutGroupCanLabel(C)
    sage: P.get_canonical_form().generator_matrix()
    [1 0 0 0 0 1 1 1 1 1 1 1 1]
    [0 1 0 1 1 0 0 1 1 2 2 1 2]
    [0 0 1 1 2 1 2 1 2 1 2 0 0]
    sage: LinearCode(P.get_transporter()*C.generator_matrix()) == P.get_canonical_form()
    True
    sage: A = P.get_autom_gens()
    sage: all(LinearCode(a*C.generator_matrix()) == C for a in A)
    True
    sage: P.get_autom_order() == GL(3, GF(3)).order()
    True

If the dimension of the dual code is smaller, we will work on this code::

    sage: C2 = codes.HammingCode(GF(3), 3)
    sage: P2 = LinearCodeAutGroupCanLabel(C2)
    sage: P2.get_canonical_form().parity_check_matrix() == P.get_canonical_form().generator_matrix()
    True

There is a specialization of this algorithm to pass a coloring on the
coordinates. This is just a list of lists, telling the algorithm which
columns do share the same coloring::

    sage: C = codes.HammingCode(GF(4, 'a'), 3).dual_code()
    sage: P = LinearCodeAutGroupCanLabel(C, P=[ [0], [1], list(range(2, C.length())) ])
    sage: P.get_autom_order()
    864
    sage: A = [a.get_perm() for a in P.get_autom_gens()]
    sage: H = SymmetricGroup(21).subgroup(A)
    sage: H.orbits()
    ((1,),
     (2,),
     (3, 5, 4),
     (6, 19, 16, 21, 9, 10, 15, 8, 20, 11, 14, 13, 7, 12, 18, 17))

We can also restrict the group action to linear isometries::

    sage: P = LinearCodeAutGroupCanLabel(C, algorithm_type='linear')
    sage: P.get_autom_order() == GL(3, GF(4, 'a')).order()
    True

and to the action of the symmetric group only::

    sage: P = LinearCodeAutGroupCanLabel(C, algorithm_type='permutational')
    sage: P.get_autom_order() == C.permutation_automorphism_group().order()
    True
"""
# ****************************************************************************
#       Copyright (C) 2012 Thomas Feulner <thomas.feulner@uni-bayreuth.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.coding.codecan.codecan import PartitionRefinementLinearCode
from sage.combinat.permutation import Permutation
from sage.arith.misc import factorial


def _cyclic_shift(n, p):
    r"""
    If ``p`` is a list of pairwise distinct coordinates in ``range(n)``,
    then this function returns the cyclic shift of
    the coordinates contained in ``p``, when acting on a vector of length `n`.

    Note that the domain of a ``Permutation`` is ``range(1, n+1)``.

    EXAMPLES::

        sage: from sage.coding.codecan.autgroup_can_label import _cyclic_shift
        sage: p = _cyclic_shift(10, [2,7,4,1]); p
        [1, 3, 8, 4, 2, 6, 7, 5, 9, 10]

    We prove that the action is as expected::

        sage: t = list(range(10))
        sage: p.action(t)
        [0, 2, 7, 3, 1, 5, 6, 4, 8, 9]
    """
    x = list(range(1, n + 1))
    for i in range(1, len(p)):
        x[p[i - 1]] = p[i] + 1
    x[p[len(p) - 1]] = p[0] + 1
    return Permutation(x)


class LinearCodeAutGroupCanLabel:
    r"""
    Canonical representatives and automorphism group computation for linear
    codes over finite fields.

    There are several notions of equivalence for linear codes:
    Let `C`, `D` be linear codes of length `n` and dimension `k`.
    The codes `C` and `D` are said to be

    - permutational equivalent, if there is some permutation `\pi \in S_n`
      such that `(c_{\pi(0)}, \ldots, c_{\pi(n-1)}) \in D` for all `c \in C`.

    - linear equivalent, if there is some permutation `\pi \in S_n` and a
      vector `\phi \in {\GF{q}^*}^n` of units of length `n` such that
      `(c_{\pi(0)} \phi_0^{-1}, \ldots, c_{\pi(n-1)} \phi_{n-1}^{-1}) \in D`
      for all `c \in C`.

    - semilinear equivalent, if there is some permutation `\pi \in S_n`, a
      vector `\phi` of units of length `n` and a field automorphism `\alpha`
      such that
      `(\alpha(c_{\pi(0)}) \phi_0^{-1}, \ldots, \alpha( c_{\pi(n-1)}) \phi_{n-1}^{-1} ) \in D`
      for all `c \in C`.

    These are group actions. This class provides an algorithm that will compute
    a unique representative `D` in the orbit of the given linear code `C`.
    Furthermore, the group element `g` with `g * C = D` and the automorphism
    group of `C` will be computed as well.

    There is also the possibility to restrict the permutational part of this
    action to a Young subgroup of `S_n`. This could be achieved by passing a
    partition `P` (as a list of lists) of the set `\{0, \ldots, n-1\}`. This is
    an option which is also available in the computation of a canonical form of
    a graph, see :meth:`sage.graphs.generic_graph.GenericGraph.canonical_label`.

    EXAMPLES::

        sage: from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
        sage: C = codes.HammingCode(GF(3), 3).dual_code()
        sage: P = LinearCodeAutGroupCanLabel(C)
        sage: P.get_canonical_form().generator_matrix()
        [1 0 0 0 0 1 1 1 1 1 1 1 1]
        [0 1 0 1 1 0 0 1 1 2 2 1 2]
        [0 0 1 1 2 1 2 1 2 1 2 0 0]
        sage: LinearCode(P.get_transporter()*C.generator_matrix()) == P.get_canonical_form()
        True
        sage: a = P.get_autom_gens()[0]
        sage: (a*C.generator_matrix()).echelon_form() == C.generator_matrix().echelon_form()
        True
        sage: P.get_autom_order() == GL(3, GF(3)).order()
        True
    """

    def __init__(self, C, P=None, algorithm_type='semilinear'):
        """
        See :class:`LinearCodeAutGroupCanLabel`.

        INPUT:

        - ``C`` -- a linear code

        - ``P`` -- (optional) a coloring of the coordinates i.e. a partition
          (list of disjoint lists) of [0 , ..., C.length()-1 ]

        - ``algorithm_type`` -- (optional) defines the acting group, either

            * ``permutational``

            * ``linear``

            * ``semilinear``

        EXAMPLES::

            sage: from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
            sage: C = codes.HammingCode(GF(2), 3).dual_code()
            sage: P = LinearCodeAutGroupCanLabel(C)
            sage: P.get_canonical_form().generator_matrix()
            [1 0 0 0 1 1 1]
            [0 1 0 1 0 1 1]
            [0 0 1 1 1 1 0]
            sage: P2 = LinearCodeAutGroupCanLabel(C, P=[[0,3,5],[1,2,4,6]],
            ....:      algorithm_type='permutational')
            sage: P2.get_canonical_form().generator_matrix()
            [1 1 1 0 0 0 1]
            [0 1 0 1 1 0 1]
            [0 0 1 0 1 1 1]
        """
        from sage.groups.semimonomial_transformations.semimonomial_transformation_group import SemimonomialTransformationGroup
        from sage.coding.linear_code import LinearCode, AbstractLinearCode

        if not isinstance(C, AbstractLinearCode):
            raise TypeError("%s is not a linear code"%C)

        self.C = C
        mat = C.generator_matrix()
        F = mat.base_ring()
        S = SemimonomialTransformationGroup(F, mat.ncols())

        if P is None:
            P = [list(range(mat.ncols()))]

        pos2P = [-1] * mat.ncols()
        for i in range(len(P)):
            P[i].sort(reverse=True)
            for x in P[i]:
                pos2P[x] = i

        col_list = mat.columns()
        nz = [i for i in range(mat.ncols()) if not col_list[i].is_zero()]
        z = [(pos2P[i], i) for i in range(mat.ncols()) if col_list[i].is_zero()]
        z.sort()
        z = [i for (p, i) in z]

        normalization_factors = [ F.one() ] * mat.ncols()
        if algorithm_type == "permutational":
            for c in col_list:
                c.set_immutable()
        else:
            for x in nz:
                normalization_factors[x] = col_list[x][(col_list[x].support())[0]]
                col_list[x] = normalization_factors[x] ** (-1) * col_list[x]
                col_list[x].set_immutable()

        normalization = S(v=normalization_factors)
        normalization_inverse = normalization ** (-1)
        col_set = list({col_list[y] for y in nz })
        col2pos = []
        col2P = []
        for c in col_set:
            X = [(pos2P[y], y) for y in range(mat.ncols()) if col_list[y] == c ]
            X.sort()
            col2pos.append([b for (a, b) in X ])
            col2P.append([a for (a, b) in X ])

        zipped = sorted(zip(col2P, col_set, col2pos))

        col2P = [qty for (qty, c, pos) in zipped]
        col_set = [c for (qty, c, pos) in zipped]
        col2pos = [pos for (qty, c, pos) in zipped]
        P_refined = []
        p = [0]
        act_qty = col2P[0]
        for i in range(1, len(col_set)):
            if act_qty == col2P[i]:
                p.append(i)
            else:
                P_refined.append(p)
                p = [i]
                act_qty = col2P[i]
        P_refined.append(p)
        # now we can start the main algorithm:
        from sage.matrix.constructor import matrix

        if len(col_set) >= 2 * mat.nrows():
            # the dimension of the code is smaller or equal than
            # the dimension of the dual code
            # in this case we work with the code itself.
            pr = PartitionRefinementLinearCode(len(col_set),
                matrix(col_set).transpose(), P=P_refined, algorithm_type=algorithm_type)

            # this command allows you some advanced debugging
            # it prints the backtrack tree -> must be activated when installing
            # pr._latex_view(title='MyTitle') #this will provide you some visual representation of what is going on

            can_transp = pr.get_transporter()
            can_col_set = pr.get_canonical_form().columns()
            self._PGammaL_autom_gens = self._compute_PGammaL_automs(pr.get_autom_gens(),
                normalization, normalization_inverse, col2pos)
            self._PGammaL_autom_size = pr.get_autom_order_permutation()
            self._PGammaL_autom_size *= pr.get_autom_order_inner_stabilizer()
            self._full_autom_order = self._PGammaL_autom_size
            self._PGammaL_autom_size /= (len(self.C.base_ring()) - 1)
        elif mat.nrows() == len(col_set):
            # it could happen (because of the recursive call, see `else`) that
            # the canonical form is the identity matrix
            n = len(col_set)
            from sage.modules.free_module import VectorSpace
            can_col_set = VectorSpace(F, n).gens()
            A = []
            self._full_autom_order = 1
            S_short = SemimonomialTransformationGroup(F, n)
            can_transp = S_short.one()
            if algorithm_type != "permutational":
                # linear or semilinear
                if algorithm_type == "semilinear":
                    A.append(S_short(autom=F.hom([F.gen() ** F.characteristic()])))
                    self._full_autom_order *= F.degree()
                for p in P_refined:
                    m = [F.one()] * n
                    m[p[0]] = F.gen()
                    A.append(S_short(v=m))
                self._full_autom_order *= (len(F) - 1) ** n
            for p in P_refined:
                if len(p) > 1:
                    A.append(S_short(perm=_cyclic_shift(n, p[:2])))
                    if len(p) > 2:
                        # cyclic shift of all elements
                        A.append(S_short(perm=_cyclic_shift(n, p)))
                    self._full_autom_order *= factorial(len(p))
            self._PGammaL_autom_size = self._full_autom_order / (len(F) - 1)
            self._PGammaL_autom_gens = self._compute_PGammaL_automs(A,
                normalization, normalization_inverse, col2pos)
        else:
            # use the dual code for the computations
            # this might have zero columns or multiple columns, hence
            # we call this algorithm again.
            short_dual_code = LinearCode(matrix(col_set).transpose()).dual_code()
            agcl = LinearCodeAutGroupCanLabel(short_dual_code,
                P=P_refined, algorithm_type=algorithm_type)
            can_transp = agcl.get_transporter()
            can_transp.invert_v()
            can_col_set = agcl.get_canonical_form().parity_check_matrix().columns()
            A = agcl.get_autom_gens()
            for a in A:
                a.invert_v()
            self._PGammaL_autom_gens = self._compute_PGammaL_automs(A,
                normalization, normalization_inverse, col2pos)
            self._PGammaL_autom_size = agcl.get_PGammaL_order()
            self._full_autom_order = agcl.get_autom_order()

        count = 0
        block_ptr = []
        canonical_form = matrix(F, mat.ncols(), mat.nrows())

        perm = [-1] * mat.ncols()
        mult = [F.one()] * mat.ncols()

        for i in range(len(can_col_set)):
            img = can_transp.get_perm()(i + 1)
            for j in col2pos[img - 1]:
                pos = P[ pos2P[j] ].pop()
                canonical_form[ pos ] = can_col_set[i]
                mult[pos] = can_transp.get_v()[i]
                perm[pos] = j + 1

        it = iter(z)
        for p in P:
            while p:
                pos = p.pop()
                perm[pos] = next(it) + 1

        self._canonical_form = LinearCode(canonical_form.transpose())
        self._transporter = S(perm=Permutation(perm), v=mult, autom=can_transp.get_autom()) * normalization
        self._trivial_autom_gens, a = self._compute_trivial_automs(normalization,
            normalization_inverse, z, [pos2P[x] for x in z], zero_column_case=True)
        self._full_autom_order *= a

        for i in range(len(col2P)):
            if len(col2P[i]) > 1:
                A, a = self._compute_trivial_automs(normalization,
                    normalization_inverse, col2pos[i], col2P[i])
                self._full_autom_order *= a
                self._trivial_autom_gens += A

    @staticmethod
    def _compute_trivial_automs(normalization, normalization_inverse, col2pos, col2P, zero_column_case=False):
        """
        In order to call the algorithm described in
        :class:`PartitionRefinementLinearCode` we remove
        zero columns and multiple occurrences of the same column (up to
        normalization).

        This function computes a set of generators for the allowed permutations
        of such a block of equal columns (depending on the partition).
        If ``zero_column_case==True`` we also add a generator for the
        multiplication by units. Furthermore, we return the order of this group.

        INPUT:

        - ``normalization`` -- an element in the semimonomial transformation group `S`

        - ``normalization_inverse`` -- the inverse of ``normalization``

        - ``col2pos`` -- list of disjoint indices in ``range(n)``

        - ``col2P`` -- an increasing list of integers, with
          ``len(col2P) == len(col2pos)`` with ``col2P[i] == col2P[j]`` if and
          only if the indices ``col2pos[i]`` and ``col2pos[j]`` are in the
          same partition

        - ``zero_column_case`` -- boolean; set to ``True`` iff we are dealing
          with the zero column

        OUTPUT:

        - list of generators in `S`

        - the order of this group

        EXAMPLES::

            sage: from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
            sage: S = SemimonomialTransformationGroup(GF(3), 10)
            sage: s = S.one()
            sage: col2pos = [1,4,2,3,5]
            sage: col2P = [1,1,3,3,3]
            sage: LinearCodeAutGroupCanLabel._compute_trivial_automs(s, s, col2pos, col2P, True)
            ([((1, 2, 1, 1, 1, 1, 1, 1, 1, 1); (), Ring endomorphism of Finite Field of size 3
              Defn: 1 |--> 1), ((1, 1, 1, 1, 1, 1, 1, 1, 1, 1); (2,5), Ring endomorphism of Finite Field of size 3
              Defn: 1 |--> 1), ((1, 1, 2, 1, 1, 1, 1, 1, 1, 1); (), Ring endomorphism of Finite Field of size 3
              Defn: 1 |--> 1), ((1, 1, 1, 1, 1, 1, 1, 1, 1, 1); (3,4), Ring endomorphism of Finite Field of size 3
              Defn: 1 |--> 1), ((1, 1, 1, 1, 1, 1, 1, 1, 1, 1); (3,4,6), Ring endomorphism of Finite Field of size 3
              Defn: 1 |--> 1)], 384)
        """
        S = normalization.parent()
        F = S.base_ring()
        n = S.degree()

        beg = 0
        if zero_column_case:
            aut_order = (len(F) - 1) ** len(col2P)
        else:
            aut_order = 1
        A = []
        while beg < len(col2P):
            if zero_column_case:
                mult = [F.one()] * n
                mult[col2pos[beg]] = F.primitive_element()
                A.append(S(v=mult))
            P_indx = col2P[beg]
            j = beg + 1
            while j < len(col2P) and col2P[j] == P_indx:
                j += 1

            if j - beg > 1:
                aut_order *= factorial(j - beg)
                # we append a transposition of the first two elements
                A.append(normalization_inverse *
                    S(perm=_cyclic_shift(n, col2pos[beg:beg + 2])) * normalization)
                if j - beg > 2:
                    # we append a cycle on all elements
                    A.append(normalization_inverse *
                        S(perm=_cyclic_shift(n, col2pos[beg:j])) * normalization)
            beg = j
        return A, aut_order

    @staticmethod
    def _compute_PGammaL_automs(gens, normalization, normalization_inverse, col2pos):
        """
        In order to call the algorithm described in
        :class:`sage.coding.codecan.codecan.PartitionRefinementLinearCode` we removed
        zero columns and multiple occurrences of the same column (up to normalization).
        This function lifts the generators ``gens`` which were returned to their full length.

        INPUT:

        - ``gens`` -- list of semimonomial transformation group elements of length `m`

        - ``normalization`` -- a semimonomial transformation of length `n`

        - ``normalization_inverse`` -- the inverse of ``normalization``

        - ``col2pos`` -- a partition of ``range(n)`` into `m` disjoint sets,
          given as a list of lists. The elements `g` in ``gens`` are only
          allowed to permute entries of ``col2pos`` of equal length.

        OUTPUT:

        - a list of semimonomial transformations containing
          ``normalization`` which are the lifts of the elements in ``gens``

        EXAMPLES::

            sage: from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
            sage: S = SemimonomialTransformationGroup(GF(3), 10)
            sage: T = SemimonomialTransformationGroup(GF(3), 5)
            sage: s = S.one()
            sage: col2pos = [[1,4,2,3,5], [0],[6],[8],[7,9]]
            sage: gens = [T(perm=Permutation([1,3,4,2,5]), v=[2,1,1,1,1])]
            sage: LinearCodeAutGroupCanLabel._compute_PGammaL_automs(gens, s, s, col2pos)
            [((1, 2, 2, 2, 2, 2, 1, 1, 1, 1); (1,7,9), Ring endomorphism of Finite Field of size 3
              Defn: 1 |--> 1)]
        """
        S = normalization.parent()
        n = S.degree()
        A = []
        for g in gens:
            perm = list(range(1, n + 1))
            mult = [S.base_ring().one()] * n
            short_perm = g.get_perm()
            short_mult = g.get_v()
            for i in range(len(col2pos)):
                c = col2pos[i]
                img_iter = iter(col2pos[short_perm(i + 1) - 1])
                for x in c:
                    perm[x] = next(img_iter) + 1
                    mult[x] = short_mult[i]
            A.append(normalization_inverse * S(perm=perm, v=mult, autom=g.get_autom()) * normalization)
        return A

    def get_canonical_form(self):
        """
        Return the canonical orbit representative we computed.

        EXAMPLES::

            sage: from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
            sage: C = codes.HammingCode(GF(3), 3).dual_code()
            sage: CF1 = LinearCodeAutGroupCanLabel(C).get_canonical_form()
            sage: s = SemimonomialTransformationGroup(GF(3), C.length()).an_element()
            sage: C2 = LinearCode(s*C.generator_matrix())
            sage: CF2 = LinearCodeAutGroupCanLabel(C2).get_canonical_form()
            sage: CF1 == CF2
            True
        """
        return self._canonical_form

    def get_transporter(self):
        """
        Return the element which maps the code to its canonical form.

        EXAMPLES::

            sage: from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
            sage: C = codes.HammingCode(GF(2), 3).dual_code()
            sage: P = LinearCodeAutGroupCanLabel(C)
            sage: g = P.get_transporter()
            sage: D = P.get_canonical_form()
            sage: (g*C.generator_matrix()).echelon_form() == D.generator_matrix().echelon_form()
            True
        """
        return self._transporter

    def get_autom_gens(self):
        """
        Return a generating set for the automorphism group of the code.

        EXAMPLES::

            sage: from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
            sage: C = codes.HammingCode(GF(2), 3).dual_code()
            sage: A = LinearCodeAutGroupCanLabel(C).get_autom_gens()
            sage: Gamma = C.generator_matrix().echelon_form()
            sage: all((g*Gamma).echelon_form() == Gamma for g in A)
            True
        """
        return self._PGammaL_autom_gens + self._trivial_autom_gens

    def get_autom_order(self):
        """
        Return the size of the automorphism group of the code.

        EXAMPLES::

            sage: from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
            sage: C = codes.HammingCode(GF(2), 3).dual_code()
            sage: LinearCodeAutGroupCanLabel(C).get_autom_order()
            168
        """
        return self._full_autom_order

    def get_PGammaL_gens(self):
        r"""
        Return the set of generators translated to the group `P\Gamma L(k,q)`.

        There is a geometric point of view of code equivalence. A linear
        code is identified with the multiset of points in the finite projective
        geometry `PG(k-1, q)`. The equivalence of codes translates to the
        natural action of `P\Gamma L(k,q)`. Therefore, we may interpret the
        group as a subgroup of `P\Gamma L(k,q)` as well.

        EXAMPLES::

            sage: from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
            sage: C = codes.HammingCode(GF(4, 'a'), 3).dual_code()
            sage: A = LinearCodeAutGroupCanLabel(C).get_PGammaL_gens()
            sage: Gamma = C.generator_matrix()
            sage: N = [ x.monic() for x in Gamma.columns() ]
            sage: all((g[0]*n.apply_map(g[1])).monic() in N for n in N for g in A)
            True
        """
        Gamma = self.C.generator_matrix()
        res = []
        for a in self._PGammaL_autom_gens:
            B = Gamma.solve_left(a * Gamma, check=True)
            res.append([B, a.get_autom()])

        return res

    def get_PGammaL_order(self):
        r"""
        Return the size of the automorphism group as a subgroup of
        `P\Gamma L(k,q)`.

        There is a geometric point of view of code equivalence. A linear
        code is identified with the multiset of points in the finite projective
        geometry `PG(k-1, q)`. The equivalence of codes translates to the
        natural action of `P\Gamma L(k,q)`. Therefore, we may interpret the
        group as a subgroup of `P\Gamma L(k,q)` as well.

        EXAMPLES::

            sage: from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
            sage: C = codes.HammingCode(GF(4, 'a'), 3).dual_code()
            sage: LinearCodeAutGroupCanLabel(C).get_PGammaL_order() == GL(3, GF(4, 'a')).order()*2/3
            True
        """
        return self._PGammaL_autom_size
