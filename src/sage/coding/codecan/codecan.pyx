# sage.doctest: needs sage.libs.pari
r"""
Canonical forms and automorphism group computation for linear codes over finite fields

We implemented the algorithm described in [Feu2009]_ which computes the unique
semilinearly isometric code (canonical form) in the equivalence class of a given
linear code ``C``. Furthermore, this algorithm will return the automorphism
group of ``C``, too.

The algorithm should be started via a further class
:class:`~sage.coding.codecan.autgroup_can_label.LinearCodeAutGroupCanLabel`.
This class removes duplicated columns (up to multiplications
by units) and zero columns. Hence, we can suppose that the input for the algorithm
developed here is a set of points in `PG(k-1, q)`.

The implementation is based on the class
:class:`sage.groups.perm_gps.partn_ref2.refinement_generic.PartitionRefinement_generic`.
See the description of this algorithm in
:mod:`sage.groups.perm_gps.partn_ref2.refinement_generic`.
In the language given there, we have to implement the group action of
`G = (GL(k,q) \times {\GF{q}^*}^n ) \rtimes Aut(\GF{q})` on the set `X =
(\GF{q}^k)^n` of `k \times n` matrices over `\GF{q}` (with the above
restrictions).

The derived class here implements the stabilizers
`G_{\Pi^{(I)}(x)}` of the projections `\Pi^{(I)}(x)` of `x` to
the coordinates specified in the sequence `I`. Furthermore, we implement
the inner minimization, i.e. the computation of a canonical form of
the projection `\Pi^{(I)}(x)` under the action of `G_{\Pi^{(I^{(i-1)})}(x)}` .
Finally, we provide suitable homomorphisms of group actions for the refinements
and methods to compute the applied group elements in `G \rtimes S_n`.

The algorithm also uses Jeffrey Leon's idea of maintaining an
invariant set of codewords which is computed in the beginning, see
:meth:`~sage.coding.codecan.codecan.PartitionRefinementLinearCode._init_point_hyperplane_incidence`.
An example for such a set is the set of all codewords of weight `\leq w` for
some uniquely defined `w`. In our case, we interpret the codewords as a set of
hyperplanes (via the corresponding information word) and compute invariants of
the bipartite, colored derived subgraph of the point-hyperplane incidence graph,
see :meth:`PartitionRefinementLinearCode._point_refine` and
:meth:`PartitionRefinementLinearCode._hyp_refine`.

Since we are interested in subspaces (linear codes) instead of matrices, our
group elements returned in
:meth:`PartitionRefinementLinearCode.get_transporter` and
:meth:`PartitionRefinementLinearCode.get_autom_gens`
will be elements in the group
`({\GF{q}^*}^n  \rtimes Aut(\GF{q})) \rtimes S_n =
({\GF{q}^*}^n  \rtimes (Aut(\GF{q}) \times S_n)`.

AUTHORS:

- Thomas Feulner (2012-11-15): initial version

REFERENCES:

- [Feu2009]

EXAMPLES:

Get the canonical form of the Simplex code::

    sage: from sage.coding.codecan.codecan import PartitionRefinementLinearCode
    sage: mat = codes.HammingCode(GF(3), 3).dual_code().generator_matrix()
    sage: P = PartitionRefinementLinearCode(mat.ncols(), mat)
    sage: cf = P.get_canonical_form(); cf
    [1 0 0 0 0 1 1 1 1 1 1 1 1]
    [0 1 0 1 1 0 0 1 1 2 2 1 2]
    [0 0 1 1 2 1 2 1 2 1 2 0 0]

The transporter element is a group element which maps the input
to its canonical form::

    sage: cf.echelon_form() == (P.get_transporter() * mat).echelon_form()
    True

The automorphism group of the input, i.e. the stabilizer under this group action,
is returned by generators::

    sage: P.get_autom_order_permutation() == GL(3, GF(3)).order()/(len(GF(3))-1)
    True
    sage: A = P.get_autom_gens()
    sage: all((a*mat).echelon_form() == mat.echelon_form() for a in A)
    True
"""

#*******************************************************************************
#       Copyright (C) 2012 Thomas Feulner <thomas.feulner@uni-bayreuth.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*******************************************************************************
from itertools import repeat
from copy import copy
from cysignals.memory cimport check_allocarray, sig_free

from sage.rings.integer cimport Integer
cimport sage.groups.perm_gps.partn_ref2.refinement_generic
from sage.modules.finite_submodule_iter cimport FiniteFieldsubspace_projPoint_iterator as FFSS_projPoint
from sage.groups.perm_gps.partn_ref.data_structures cimport *
from sage.data_structures.bitset_base cimport *


cdef class InnerGroup:
    r"""
    This class implements the stabilizers `G_{\Pi^{(I)}(x)}` described in
    :mod:`sage.groups.perm_gps.partn_ref2.refinement_generic` with
    `G = (GL(k,q) \times \GF{q}^n ) \rtimes Aut(\GF{q})`.

    Those stabilizers can be stored as triples:

    - ``rank`` -- integer in `\{0, \ldots, k\}`
    - ``row_partition`` -- a partition of `\{0, \ldots, k-1\}` with
      discrete cells for all integers `i` `\geq` ``rank``
    - ``frob_pow`` -- integer `s` in `\{0, \ldots, r-1\}` if `q = p^r`

    The group `G_{\Pi^{(I)}(x)}` contains all elements `(A, \varphi, \alpha) \in G`,
    where

    - `A` is a `2 \times 2` blockmatrix, whose upper left matrix
      is a `k \times k` diagonal matrix whose entries `A_{i,i}` are constant
      on the cells of the partition ``row_partition``.
      The lower left matrix is zero.
      And the right part is arbitrary.
    - The support of the columns given by `i \in I` intersect exactly one
      cell of the partition. The entry `\varphi_i` is equal to the entries
      of the corresponding diagonal entry of `A`.
    - `\alpha` is a power of `\tau^s`, where `\tau` denotes the
      Frobenius automorphism of the finite field `\GF{q}` and `s` = ``frob_pow``.

    See [Feu2009]_ for more details.
    """
    def __cinit__(self, k=0, algorithm='semilinear', **kwds):
        r"""
        See :class:`sage.coding.codecan.codecan.InnerGroup`.

        INPUT:

        - ``k`` -- integer; gives the dimension of the matrix component
        - ``algorithm`` -- either

            * "semilinear" --  full group
            * "linear" -- no field automorphisms, i.e. `G = (GL(k,q) \times \GF{q}^n )`
            * "permutational -- no field automorphisms and no column multiplications
              i.e. `G = GL(k,q)`

        - ``transporter`` -- (optional) set to an element of the group
          :class:`sage.groups.semimonomial_transformations.semimonomial_transformation_group.SemimonomialTransformationGroup`
          if you would like to modify this element simultaneously

        EXAMPLES::

            sage: from sage.coding.codecan.codecan import InnerGroup
            sage: IG = InnerGroup(10)
            sage: IG = InnerGroup(10, "linear")
            sage: IG = InnerGroup(10, "permutational")

        ::

            sage: S = SemimonomialTransformationGroup(GF(4, 'a'), 8)
            sage: IG = InnerGroup(3, transporter=S.an_element())
        """

        self.rank = 0
        if k > 0:
            self.row_partition = OP_new(k)

        if algorithm == "permutational":
            self.frob_pow = 0
            self.permutational_only = 1
            for i in range(1, k):
                OP_join(self.row_partition, 0, i)
        else:
            self.permutational_only = 0
            if algorithm == "semilinear":
                self.frob_pow = 1
            elif algorithm == "linear":
                self.frob_pow = 0

        self.compute_transporter = False
        if "transporter" in kwds:
            self.transporter = kwds["transporter"]
            self.compute_transporter = True

    def __dealloc__(self):
        r"""
        Deallocate ``self``.
        """
        OP_dealloc(self.row_partition)

    cdef int get_rep(self, int pos) noexcept:
        """
        Get the index of the cell of ``self.row_partition`` containing ``pos``.
        """
        return OP_find(self.row_partition, pos)

    cdef bint has_semilinear_action(self) noexcept:
        """
        Return ``True`` iff the field automorphism group component of ``self``
        is non-trivial.
        """
        return (self.frob_pow > 0)

    cdef int join_rows(self, int rep1, int rep2) noexcept:
        """
        Join the cells with unique representatives
        ``rep1`` and ``rep2`` of ``self.row_partition``.
        Return the index of the join.
        """
        OP_join(self.row_partition, rep1, rep2)
        return self.get_rep(rep1)

    cdef void copy_from(self, InnerGroup other) noexcept:
        """
        Copy the group ``other`` to ``self``.
        """
        self.rank = other.rank
        self.frob_pow = other.frob_pow
        self.permutational_only = other.permutational_only
        OP_copy_from_to(other.row_partition, self.row_partition)

    cdef minimize_by_row_mult(self, FreeModuleElement w):
        r"""
        We suppose `v \in \GF{q}^k` and the entries `v_i = 0` for all
        ``i >= self.rank``.
        We compute the smallest element ``w`` in the orbit of ``v`` under the
        group action of the matrix components of all elements in ``self``.

        We return ``d, w``, where ``d`` is a dictionary mapping
        ``self.row_partition`` (accessed via their unique representatives)
        to its necessary multiplication. Non-occurring cells are multiplicated
        by 1.
        """
        cdef FreeModuleElement v = w.__copy__()
        cdef dict d = dict()
        if self.permutational_only:
            return d, v

        cdef list nz_pos = v.nonzero_positions()
        for r in nz_pos:
            r_rep = self.get_rep(r)
            if r_rep not in d:
                d[r_rep] = v[r] ** (-1)
                v[r] = 1
            else:
                v[r] *= d[r_rep]
        return d, v

    cdef minimize_matrix_col(self, object m, int pos, list fixed_minimized_cols,
                             bint *group_changed):
        r"""
        Minimize the column at position ``pos`` of the matrix ``m`` by the
        action of ``self``. ``m`` should have no zero column. ``self`` is set to
        the stabilizer of this column.

        We return ``group_changed, mm`` where ``group_changed`` is a boolean
        indicating if ``self`` got changed and ``mm`` is the modification of
        ``m``.

        In the array ``fixed_minimized_cols`` we store, those
        columns of ``m`` which are known to be invariant under ``self``.
        """
        group_changed[0] = False
        cdef SemimonomialTransformation my_trans
        cdef FreeModuleElement act_col = m.column(pos)
        cdef int pivot = -1
        cdef list nz_pos = act_col.nonzero_positions()
        cdef int applied_frob, i, col, row, first_nz_rep

        F = m.base_ring()

        for i in nz_pos:
            if i >= self.rank:
                pivot = i
                break

        if pivot == -1:
            if self.permutational_only:
                return m
            # this column is linearly dependent on those already fixed
            first_nz = nz_pos.pop(0)
            first_nz_rep = self.get_rep(first_nz)
            factor = m[first_nz, pos] ** (-1)
            m.rescale_col(pos, factor)

            if self.compute_transporter:
                n = self.transporter.parent().degree()
                v = (F.one(),)*(pos) + (factor**(-1), ) + (F.one(),)*(n-pos-1)
                my_trans = self.transporter.parent()(v=v)

            d, _ = self.minimize_by_row_mult(factor * act_col)
            d.pop(first_nz_rep)
            if len(d):  # there is at least one more multiplication
                group_changed[0] = True
                for i in range(self.rank):
                    factor = d.get(self.get_rep(i))
                    if factor and not factor.is_zero():
                        m.rescale_row(i, factor)
                for i in d:
                    first_nz_rep = self.join_rows(first_nz_rep, i)
                # rescale the already fixed part by column multiplications
                for col in fixed_minimized_cols:
                    col_nz = m.column(col).nonzero_positions()
                    if col_nz:
                        row = col_nz[0]
                        if self.compute_transporter:
                            my_trans.v = (my_trans.v[:col] + (m[row, col],) +
                                          my_trans.v[col+1:])
                        m.rescale_col(col, m[row, col] ** (-1))
            if self.has_semilinear_action():
                applied_frob = 0
                self.minimize_by_frobenius(m[nz_pos].column(pos), &applied_frob, &self.frob_pow)
                f = F.hom([F.gen() ** (F.characteristic() ** applied_frob)])
                m = m.apply_map(f)  # this would change the reference!

                if self.compute_transporter:
                    my_trans.v = tuple([my_trans.v[i].frobenius(applied_frob)
                                        for i in range(len(my_trans.v))])
                    my_trans.alpha = f

            if self.compute_transporter:
                self.transporter = my_trans * self.transporter
        else:
            # this column is linearly independent on those already fixed,
            # map it to the self._rank-th unit vector
            group_changed[0] = True
            self.gaussian_elimination(m, pos, pivot, nz_pos)
            self.rank += 1
        return m

    cdef void gaussian_elimination(self, object m, int pos, int pivot, list nz_pos) noexcept:
        r"""
        Minimize the column at position ``pos`` of the matrix ``m`` by the
        action of ``self``. We know that there is some nonzero entry of this
        column at ``pivot >= self.rank``. All nonzero entries are stored in
        the list ``nz_pos``.

        ``self`` is not modified by this function, but ``m`` is.
        """
        nz_pos.remove(pivot)
        m.rescale_row(pivot, m[pivot, pos] ** (-1))

        for r in nz_pos:
            m.add_multiple_of_row(r, pivot, -m[r, pos])  # Gaussian elimination
        if pivot != self.rank:
            m.swap_rows(self.rank, pivot)

    cdef InnerGroup _new_c(self):
        r"""
        Make a new copy of ``self``.
        """
        cdef InnerGroup res = InnerGroup()
        res.frob_pow = self.frob_pow
        res.rank = self.rank
        res.row_partition = OP_copy(self.row_partition)
        res.permutational_only = self.permutational_only
        return res

    cdef SemimonomialTransformation get_transporter(self):
        r"""
        Return the group element we have applied. Should only be called if
        we passed an element in
        :meth:`sage.coding.codecan.codecan.InnerGroup.__cinit__`.
        """
        return self.transporter

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from sage.coding.codecan.codecan import InnerGroup
            sage: InnerGroup(10)
            Subgroup of (GL(k,q) times \GF{q}^n ) rtimes Aut(\GF{q}) with rank = 0,
            frobenius power = 1 and partition = 0 -> 0 1 -> 1 2 -> 2 3 -> 3 4 -> 4 5 -> 5
            6 -> 6 7 -> 7 8 -> 8 9 -> 9
        """
        return r"Subgroup of (GL(k,q) times \GF{q}^n ) rtimes Aut(\GF{q}) " + \
            "with rank = %s, frobenius power = %s and partition =%s" % (self.rank,
            self.frob_pow, OP_string(self.row_partition))

    cdef void minimize_by_frobenius(self, object v, int *applied_frob, int *stab_pow) noexcept:
        r"""
        Minimize the vector ``v \in \GF{q}^k`` by the
        action of the field automorphism component of ``self``.
        ``self`` and ``v`` are not modified by this function.

        Let `\tau` denote the Frobenius automorphism of ``\GF{q}``. Then
        ``applied_frob``-th power of `\tau` will give us the minimal element.
        The ``stab_pow``-th power of `\tau` will generate the stabilizer of `v`.
        """
        stab_pow[0] = self.frob_pow
        applied_frob[0] = 0
        cdef int loc_frob, min_pow = 0
        for el in v:
            x = el.frobenius(applied_frob[0])
            y = x  # the elements in the cyclic(!) orbit
            m = x  # a candidate for the minimal element

            loc_frob = 0
            min_pow = 0

            while True:
                loc_frob += 1
                y = y.frobenius(stab_pow[0])
                if y == x:
                    break
                if y < m:
                    m = y
                    min_pow = loc_frob

            # now x.frobenius(stab_pow*loc_frob) == x
            applied_frob[0] += min_pow * stab_pow[0]
            stab_pow[0] *= loc_frob
            if stab_pow[0] == el.parent().degree():
                stab_pow[0] = 0
                break  # for

    cpdef int get_frob_pow(self) noexcept:
        r"""
        Return the power of the Frobenius automorphism which generates
        the corresponding component of ``self``.

        EXAMPLES::

            sage: from sage.coding.codecan.codecan import InnerGroup
            sage: I = InnerGroup(10)
            sage: I.get_frob_pow()
            1
        """
        return self.frob_pow

    cpdef column_blocks(self, mat):
        r"""
        Let ``mat`` be a matrix which is stabilized by ``self`` having no zero
        columns. We know that for each column of ``mat`` there is a uniquely
        defined cell in ``self.row_partition`` having a nontrivial intersection
        with the support of this particular column.

        This function returns a partition (as list of lists) of the columns
        indices according to the partition of the rows given by ``self``.

        EXAMPLES::

            sage: from sage.coding.codecan.codecan import InnerGroup
            sage: I = InnerGroup(3)
            sage: mat = Matrix(GF(3), [[0,1,0],[1,0,0], [0,0,1]])
            sage: I.column_blocks(mat)
            [[1], [0], [2]]
        """
        if self.row_partition.num_cells == 1:
            return [list(range(mat.ncols()))]

        r = [[] for _ in repeat(None, mat.ncols())]
        cols = iter(mat.columns())
        for i in range(mat.ncols()):
            # there should be no zero columns by assumption!
            m = OP_find(self.row_partition, next(cols).nonzero_positions()[0])
            r[m].append(i)
        return [x for x in r if x]


cdef class PartitionRefinementLinearCode(PartitionRefinement_generic):
    """
    See :mod:`sage.coding.codecan.codecan`.

    EXAMPLES::

        sage: from sage.coding.codecan.codecan import PartitionRefinementLinearCode
        sage: mat = codes.HammingCode(GF(3), 3).dual_code().generator_matrix()
        sage: P = PartitionRefinementLinearCode(mat.ncols(), mat)
        sage: cf = P.get_canonical_form(); cf
        [1 0 0 0 0 1 1 1 1 1 1 1 1]
        [0 1 0 1 1 0 0 1 1 2 2 1 2]
        [0 0 1 1 2 1 2 1 2 1 2 0 0]

    ::

        sage: cf.echelon_form() == (P.get_transporter() * mat).echelon_form()
        True

    ::

        sage: P.get_autom_order_permutation() == GL(3, GF(3)).order()/(len(GF(3))-1)
        True
        sage: A = P.get_autom_gens()
        sage: all((a*mat).echelon_form() == mat.echelon_form() for a in A)
        True
    """
    def __cinit__(self):
        r"""
        Initialization. See :meth:`__init__`.

        TESTS::

            sage: from sage.coding.codecan.codecan import PartitionRefinementLinearCode
            sage: C = PartitionRefinementLinearCode.__new__(PartitionRefinementLinearCode, 0)
        """
        self._hyp2points = NULL
        self._points2hyp = NULL
        self._hyp_part = NULL
        self._hyp_refine_vals_scratch = NULL
        self._nr_of_supp_refine_calls = 0
        self._nr_of_point_refine_calls = 0
        self._stored_states = dict()

    def __init__(self, n, generator_matrix, P=None, algorithm_type='semilinear'):
        r"""
        Initialization, we immediately start the algorithm
        (see :mod:`sage.coding.codecan.codecan`)
        to compute the canonical form and automorphism group of the linear code
        generated by ``generator_matrix``.

        INPUT:

        - ``n`` -- integer
        - ``generator_matrix`` -- a `k \times n` matrix over `\GF{q}` of full row rank,
          i.e. `k<n` and without zero columns
        - partition --  (optional) a partition (as list of lists) of the set
          `\{0, \ldots, n-1\}` which restricts the action of the permutational
          part of the group to the stabilizer of this partition
        - algorithm_type -- (optional) use one of the following options

          * "semilinear" --  full group
          * "linear" -- no field automorphisms, i.e. `G = (GL(k,q) \times \GF{q}^n )`
          * "permutational - no field automorphisms and no column multiplications
            i.e. `G = GL(k,q)`

        EXAMPLES::

            sage: from sage.coding.codecan.codecan import PartitionRefinementLinearCode
            sage: mat = codes.HammingCode(GF(3), 3).dual_code().generator_matrix()
            sage: P = PartitionRefinementLinearCode(mat.ncols(), mat)

        ::

            sage: from sage.coding.codecan.codecan import PartitionRefinementLinearCode
            sage: mat = codes.HammingCode(GF(3), 3).dual_code().generator_matrix()
            sage: P = PartitionRefinementLinearCode(mat.ncols(), mat)
        """
        self._k = generator_matrix.nrows()
        self._q = len(generator_matrix.base_ring())
        self._matrix = copy(generator_matrix)
        self._root_matrix = generator_matrix
        self._supp_refine_vals = _BestValStore(n)
        self._point_refine_vals = _BestValStore(n)
        # self._hyp_refine_vals will initialized after
        # we computed the set of codewords

        self._run(P, algorithm_type)

    def __dealloc__(self):
        r"""
        Deallocate ``self``.
        """
        cdef int i
        if self._points2hyp is not NULL:
            for i in range(self._n):
                bitset_free(self._points2hyp[i])
            sig_free(self._points2hyp)

        if self._points2hyp is not NULL:
            for i in range(self._hyp_part.degree):
                bitset_free(self._hyp2points[i])
            sig_free(self._hyp2points)

        PS_dealloc(self._hyp_part)
        sig_free(self._hyp_refine_vals_scratch)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.coding.codecan.codecan import PartitionRefinementLinearCode
            sage: mat = codes.HammingCode(GF(3), 3).dual_code().generator_matrix()
            sage: PartitionRefinementLinearCode(mat.ncols(), mat)
            Canonical form algorithm for linear code generated by
            [1 0 1 1 0 1 0 1 1 1 0 1 1]
            [0 1 1 2 0 0 1 1 2 0 1 1 2]
            [0 0 0 0 1 1 1 1 1 2 2 2 2]
        """
        return "Canonical form algorithm for linear code generated" + \
            " by\n%s" % (self._root_matrix)

    def _run(self, P, algorithm_type):
        """
        Start the main algorithm, this method is called in :meth:`init`.
        See this method for the description of the input.

        EXAMPLES::

            sage: from sage.coding.codecan.codecan import PartitionRefinementLinearCode
            sage: mat = codes.HammingCode(GF(3), 3).dual_code().generator_matrix()
            sage: P = PartitionRefinementLinearCode(mat.ncols(), mat) #indirect doctest
            sage: P.get_canonical_form()
            [1 0 0 0 0 1 1 1 1 1 1 1 1]
            [0 1 0 1 1 0 0 1 1 2 2 1 2]
            [0 0 1 1 2 1 2 1 2 1 2 0 0]
        """
        self._init_point_hyperplane_incidence()
        F = self._matrix.base_ring()
        if F.order() == 2:
            algorithm_type = "permutational"
        elif self._matrix.base_ring().is_prime_field() and algorithm_type != "permutational":
            algorithm_type = "linear"
        self._inner_group = InnerGroup(self._k, algorithm_type)

        self._init_partition_stack(P)
        self._init_point_hyperplane_incidence()
        self._start_Sn_backtrack()  # start the main computation

        # up to now, we just computed the permutational part of the group action
        # compute the other components of the transporter
        from sage.combinat.permutation import Permutation
        from sage.groups.semimonomial_transformations.semimonomial_transformation_group import SemimonomialTransformationGroup
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup

        S = SemimonomialTransformationGroup(self._matrix.base_ring(), self._n)
        S_n = SymmetricGroup(self._n)

        self._transporter = S(perm= S_n(self._to_best.sage()))
        self._transporter, self._best_candidate, remaining_inner_group = self._compute_group_element(self._transporter, algorithm_type)

        # compute the other components of the automorphism group generators
        self._autom_group_generators = []
        transp_inv = self._transporter ** (-1)

        for a in self._known_automorphisms.small_generating_set():
            x = S(perm=self._transporter.get_perm() * Permutation(S_n(a)))
            x, _, _ = self._compute_group_element(x, algorithm_type)
            self._autom_group_generators.append(transp_inv * x)

        if algorithm_type == "permutational":
            self._inner_group_stabilizer_order = 1
        else:
            P = remaining_inner_group.column_blocks(self._best_candidate)
            for p in P:
                x = S(v=[ F.primitive_element() if i in p else F.one()  for i in range(self._n) ])
                self._autom_group_generators.append(transp_inv * x * self._transporter)
            self._inner_group_stabilizer_order = (len(F) - 1) ** len(P)

        if remaining_inner_group.get_frob_pow() > 0:
            x = S(autom=F.hom([F.primitive_element() ** (remaining_inner_group.get_frob_pow() * F.characteristic())]))
            self._autom_group_generators.append(transp_inv * x * self._transporter)
            self._inner_group_stabilizer_order *= Integer(F.degree() / remaining_inner_group.get_frob_pow())

    cdef _compute_group_element(self, SemimonomialTransformation trans, str algorithm_type):
        """
        Apply ``trans`` to ``self._root_matrix`` and minimize this matrix
        column by column under the inner minimization. The action is
        simultaneously applied to ``trans``.

        The output of this function is a triple containing, the modified
        group element ``trans``, the minimized matrix and the stabilizer of this
        matrix under the inner group.
        """
        cdef InnerGroup inner_group = InnerGroup(self._k, algorithm_type, transporter=trans)
        cdef bint group_changed = False
        cdef int i
        cdef list fixed_pos = []
        mat = trans * self._root_matrix

        for i in range(self._n):
            mat = inner_group.minimize_matrix_col(mat, i, fixed_pos,
                                                  &group_changed)
            fixed_pos.append(i)

        trans = inner_group.get_transporter()
        return trans, mat, inner_group

    def get_canonical_form(self):
        r"""
        Return the canonical form for this matrix.

        EXAMPLES::

            sage: from sage.coding.codecan.codecan import PartitionRefinementLinearCode
            sage: mat = codes.HammingCode(GF(3), 3).dual_code().generator_matrix()
            sage: P1 = PartitionRefinementLinearCode(mat.ncols(), mat)
            sage: CF1 = P1.get_canonical_form()
            sage: s = SemimonomialTransformationGroup(GF(3), mat.ncols()).an_element()
            sage: P2 = PartitionRefinementLinearCode(mat.ncols(), s*mat)
            sage: CF1 == P2.get_canonical_form()
            True
        """
        return self._best_candidate

    def get_transporter(self):
        """
        Return the transporter element, mapping the initial matrix to its
        canonical form.

        EXAMPLES::

            sage: from sage.coding.codecan.codecan import PartitionRefinementLinearCode
            sage: mat = codes.HammingCode(GF(3), 3).dual_code().generator_matrix()
            sage: P = PartitionRefinementLinearCode(mat.ncols(), mat)
            sage: CF = P.get_canonical_form()
            sage: t = P.get_transporter()
            sage: (t*mat).echelon_form() == CF.echelon_form()
            True
        """
        return self._transporter

    def get_autom_gens(self):
        """
        Return generators of the automorphism group of the initial matrix.

        EXAMPLES::

            sage: from sage.coding.codecan.codecan import PartitionRefinementLinearCode
            sage: mat = codes.HammingCode(GF(3), 3).dual_code().generator_matrix()
            sage: P = PartitionRefinementLinearCode(mat.ncols(), mat)
            sage: A = P.get_autom_gens()
            sage: all((a*mat).echelon_form() == mat.echelon_form() for a in A)
            True
        """
        return self._autom_group_generators

    def get_autom_order_inner_stabilizer(self):
        """
        Return the order of the stabilizer of the initial matrix under
        the action of the inner group `G`.

        EXAMPLES::

            sage: from sage.coding.codecan.codecan import PartitionRefinementLinearCode
            sage: mat = codes.HammingCode(GF(3), 3).dual_code().generator_matrix()
            sage: P = PartitionRefinementLinearCode(mat.ncols(), mat)
            sage: P.get_autom_order_inner_stabilizer()
            2
            sage: mat2 = Matrix(GF(4, 'a'), [[1,0,1], [0,1,1]])
            sage: P2 = PartitionRefinementLinearCode(mat2.ncols(), mat2)
            sage: P2.get_autom_order_inner_stabilizer()
            6
        """
        return self._inner_group_stabilizer_order

    cdef _init_point_hyperplane_incidence(self):
        r"""
        Compute a set of codewords `W` of `C` (generated by self) which is compatible
        with the group action, i.e. if we start with some other code `(g,\pi)C`
        the result should be `(g,\pi)W`.

        The set `W` will consist of all normalized codewords up to some weight
        `w`, where `w` is the smallest integer such that `W` spans the linear code `C`.

        This set is then transformed to an incidence matrix ``self._points2hyp``
        of the point-hyperplane graph (points correspond to rows, hyperplanes to
        columns). The hyperplanes correspond to the information
        words. For performance reasons, we also store the transpose
        ``self._hyp2points`` of ``self._points2hyp``.

        This graph will be later used in the refinement procedures.
        """
        cdef FFSS_projPoint iter = FFSS_projPoint(self._matrix)
        cdef mp_bitcnt_t i,j

        ambient_space = (self._matrix.base_ring()) ** (self._n)
        weights2size = [0] * (self.len() + 1)
        W = [[] for _ in repeat(None, self.len() + 1)]
        span = [ambient_space.zero_subspace()] * (self.len() + 1)
        min_weight = self.len()
        max_weight = self.len()

        while True:  # compute an invariant set of (normalized) codewords which span the subspace
            try:
                cw = next(iter)
            except StopIteration:
                break
            w = cw.hamming_weight()
            if min_weight > w:
                min_weight = w
            if w <= max_weight:
                X = ambient_space.subspace([cw])
                for i in range(w, max_weight):
                    old_dim = span[i].dimension()
                    span[i] += X
                    if span[i].dimension() == old_dim:
                        break  # this will also be the case for all others
                    if old_dim + 1 == self._k:
                        # the codewords of weight <= max_weight span the code
                        max_weight = i
                        break
                W[w].append(cw)

        flat_W = sum(W[min_weight: max_weight + 1], [])
        cdef int __hyp2points_size = len(flat_W)
        self._hyp_part = PS_new(__hyp2points_size, 1)
        s = -1
        for x in W[min_weight: max_weight]:
            s += len(x)
            if s >= 0:
                self._hyp_part.levels[s] = 0

        self._points2hyp = <bitset_t*>check_calloc(self._n, sizeof(bitset_t))
        for i in range(self._n):
            bitset_init(self._points2hyp[i], self._hyp_part.degree)

        self._hyp2points = <bitset_t*>check_calloc(self._hyp_part.degree, sizeof(bitset_t))
        for i in range(self._hyp_part.degree):
            bitset_init(self._hyp2points[i], self._n)
            for j in flat_W[i].support():
                bitset_add(self._hyp2points[i], j)
                bitset_add(self._points2hyp[j], i)

        self._hyp_refine_vals_scratch = <long*>check_allocarray(
                self._hyp_part.degree, sizeof(long))

        self._hyp_refine_vals = _BestValStore(self._hyp_part.degree)

    cdef bint _minimization_allowed_on_col(self, int pos) noexcept:
        r"""
        Decide if we are allowed to perform the inner minimization on position
        ``pos`` which is supposed to be a singleton. For linear codes over finite
        fields, we can always return ``True``.
        """
        return True

    cdef bint _inner_min_(self, int pos, bint *inner_group_changed) noexcept:
        r"""
        Minimize the node by the action of the inner group on the ``pos``-th position.
        Sets ``inner_group_changed`` to ``True`` if and only if the inner group
        has changed.

        INPUT:

        - ``pos`` -- a position in  ``range(self.n)``

        OUTPUT:

        - ``True`` if and only if the actual node compares less or equal
          to the candidate for the canonical form.
        """
        self._matrix = self._inner_group.minimize_matrix_col(self._matrix, pos,
                            self._fixed_minimized, inner_group_changed)

        # finally compare the new column with the best candidate
        if self._is_candidate_initialized:
            A = self._matrix.column(pos)
            B = self._best_candidate.column(
                self._inner_min_order_best[len(self._fixed_minimized)])
            if B < A:
                return False
            if A < B:
                # the next leaf will become the next candidate
                self._is_candidate_initialized = False
        return True

    cdef bint _refine(self, bint *part_changed,
                      bint inner_group_changed, bint first_step) noexcept:
        """
        Refine the partition ``self.part``. Set  ``part_changed`` to ``True``
        if and only if ``self.part`` was refined.

        OUTPUT: ``False`` -- only if the actual node compares larger than the
        candidate for the canonical form
        """
        part_changed[0] = False
        cdef bint res, hyp_part_changed = not first_step
        cdef bint n_partition_changed = first_step
        cdef bint n_partition_changed_copy = n_partition_changed

        while hyp_part_changed or n_partition_changed:
            inner_group_changed = False
            res = self._inner_min_refine(&inner_group_changed, &n_partition_changed)
            if not res:
                return False

            part_changed[0] |= n_partition_changed
            n_partition_changed = n_partition_changed_copy
            n_partition_changed_copy = True

            if n_partition_changed:
                if PS_is_discrete(self._part):
                    return True
                if inner_group_changed:
                    continue

            while hyp_part_changed or n_partition_changed:
                if n_partition_changed:
                    res = self._hyp_refine(&hyp_part_changed)
                    if not res:
                        return False
                    n_partition_changed = False
                else:
                    res = self._point_refine(&inner_group_changed, &n_partition_changed)
                    if not res:
                        return False
                    part_changed[0] |= n_partition_changed
                    hyp_part_changed = False
                    if inner_group_changed:
                        break  # perform the inner_min_refine first!
                    if n_partition_changed and PS_is_discrete(self._part):
                        return True
        return True

    cdef bint _inner_min_refine(self, bint *inner_stab_changed, bint *changed_partition) noexcept:
        """
        Refine the partition ``self.part`` by computing the orbit (respectively
        the hash of a canonical form) of each column vector under the inner group.

        New fixed points of ``self.part`` get refined by the inner group. If this
        leads to a smaller group then we set ``inner_stab_changed`` to ``True``.
        ``changed_partition``  is set to ``True`` if and only if ``self.part``
        was refined.

        OUTPUT:

        - ``False`` only if the image under this homomorphism of group actions
          compares larger than the image of the candidate for the canonical form.
        """
        cdef int i, j, res, stab_pow, apply_pow

        if self._inner_group.rank < 2:
            return True

        lower = iter(self._matrix[:self._inner_group.rank].columns())
        upper = iter(self._matrix[self._inner_group.rank:].columns())

        for i in range(self._n):
            l = next(lower)
            u = next(upper)

            if u.is_zero() and i not in self._fixed_minimized:
                # minimize by self._inner_group as in _inner_min:
                _, l = self._inner_group.minimize_by_row_mult(l)

                if self._inner_group.has_semilinear_action():
                    stab_pow = self._inner_group.frob_pow
                    apply_pow = 0
                    self._inner_group.minimize_by_frobenius(l, &apply_pow, &stab_pow)

                    F = self._matrix.base_ring()
                    f = F.hom([F.gen() ** (F.characteristic() ** apply_pow)])
                    l = l.apply_map(f)
                res = 0
                for r in iter(l):
                    res *= self._q
                    res += hash(r)
                self._refine_vals_scratch[i] = res
            else:
                self._refine_vals_scratch[i] = -1

        # provide some space to store the result (if not already exists)
        cdef long * best_vals = self._supp_refine_vals.get_row(self._nr_of_supp_refine_calls)
        self._nr_of_supp_refine_calls += 1
        return self._one_refinement(best_vals, 0, self._n, inner_stab_changed,
                                    changed_partition, "supp_refine")

    cdef bint _point_refine(self, bint *inner_stab_changed, bint *changed_partition) noexcept:
        """
        Refine the partition ``self.part`` by counting
        (colored) neighbours in the point-hyperplane graph.

        New fixed points of ``self.part`` get refined by the inner group. If this
        leads to a smaller group then we set ``inner_stab_changed`` to ``True``.
        ``changed_partition``  is set to ``True`` if and only if ``self.part``
        was refined.

        OUTPUT:

        - ``False`` only if the image under this homomorphism of group actions
          compares larger than the image of the candidate for the canonical form.
        """

        self._part.depth += 1
        PS_clear(self._part)

        cdef bitset_t *nonsingletons = NULL
        cdef bitset_t scratch
        bitset_init(scratch, self._hyp_part.degree)
        cdef int nr_cells = PS_all_new_cells(self._hyp_part, & nonsingletons)

        for i in range(self._n):
            res = [0] * nr_cells
            for j in range(nr_cells):
                bitset_and(scratch, self._points2hyp[i], nonsingletons[j])
                res[j] = bitset_hamming_weight(scratch)
            self._refine_vals_scratch[i] = hash(tuple(res))

        for j in range(nr_cells):
            bitset_free(nonsingletons[j])
        sig_free(nonsingletons)
        bitset_free(scratch)

        # provide some space to store the result (if not already exists)
        cdef long * best_vals = self._point_refine_vals.get_row(self._nr_of_point_refine_calls)
        self._nr_of_point_refine_calls += 1
        cdef bint ret_val = self._one_refinement(best_vals, 0, self._n,
            inner_stab_changed, changed_partition, "point_refine")

        if not changed_partition[0]:
            self._part.depth -= 1
        return ret_val

    cdef bint _hyp_refine(self, bint *changed_partition) noexcept:
        """
        Refine the partition of the hyperplanes by counting
        (colored) neighbours in the point-hyperplane graph.

        ``changed_partition``  is set to ``True`` if and only if the partition
        was refined.

        OUTPUT:

        - ``False`` only if the image under this homomorphism of group actions
          compares larger than the image of the candidate for the canonical form.
        """
        self._hyp_part.depth += 1
        PS_clear(self._hyp_part)
        cdef bitset_t *nonsingletons = NULL
        cdef bitset_t scratch
        bitset_init(scratch, self._part.degree)
        cdef int nr_cells = PS_all_new_cells(self._part, & nonsingletons)

        for i in range(self._hyp_part.degree):
            res = [0] * nr_cells
            for j in range(nr_cells):
                bitset_and(scratch, self._hyp2points[i], nonsingletons[j])
                res[j] = bitset_hamming_weight(scratch)
            self._hyp_refine_vals_scratch[i] = hash(tuple(res))

        for j in range(nr_cells):
            bitset_free(nonsingletons[j])
        sig_free(nonsingletons)
        bitset_free(scratch)

        # provide some space to store the result (if not already exists)
        cdef long * best_vals = self._hyp_refine_vals.get_row(self._nr_of_hyp_refine_calls)
        self._nr_of_hyp_refine_calls += 1

        cdef tuple ret_val = PS_refinement(self._hyp_part,
            self._hyp_refine_vals_scratch, best_vals, 0, self._hyp_part.degree,
            &self._is_candidate_initialized, changed_partition)

        if not changed_partition[0]:
            self._hyp_part.depth -= 1
        return ret_val[0]

    cdef tuple _store_state_(self):
        r"""
        Store the current state of the node to a tuple, such that it can be
        restored by :meth:`_restore_state_`.
        """
        return (self._inner_group._new_c(), self._nr_of_supp_refine_calls,
                self._nr_of_point_refine_calls, self._nr_of_hyp_refine_calls,
                self._hyp_part.depth)

    cdef void _restore_state_(self, tuple act_state) noexcept:
        r"""
        The inverse of :meth:`_store_state_`.
        """
        self._inner_group.copy_from(act_state[0])
        self._nr_of_supp_refine_calls = act_state[1]
        self._nr_of_point_refine_calls = act_state[2]
        self._nr_of_hyp_refine_calls = act_state[3]
        self._hyp_part.depth = act_state[4]

    cdef void _store_best_(self) noexcept:
        """
        Store this node as the actual best candidate for the canonical form.
        """
        self._best_candidate = copy(self._matrix)

    cdef void _latex_act_node(self, str comment='', int printlvl=0) noexcept:
        """
        Print the actual status as latex (tikz) commands to
        ``self._latex_debug_string``. Only needed if one wants to visualize
        the algorithm using latex. This is very helpful for debugging purposes.
        """
        if not BACKTRACK_WITHLATEX_DEBUG:
            return

        self._latex_debug_string += "\\node{\\begin{tabular}{"
        cdef int i, j, last = -1
        divide_sign = ""
        for i from 0 <= i < self._part.degree:
            if self._part.levels[i] <= self._part.depth:
                self._latex_debug_string += divide_sign + "*{" + str(i - last) + "}{c}"
                last = i
                divide_sign = "|"
        self._latex_debug_string += "}"

        # Print the applied permutation. We do highlight the fixed positions.
        for i from 0 <= i < (self._n):
            if self._part.entries[i] in self._fixed_minimized:
                self._latex_debug_string += "\\color{red}"

            self._latex_debug_string += str(self._part.entries[i])
            if i == self._n - 1:
                self._latex_debug_string += " \\\\\\hline\n"
            else:
                self._latex_debug_string += " & "

        permuted_matrix = self._matrix.matrix_from_columns([self._part.entries[i] for i in range(self._n) ])

        # Now we will finally print the matrix.
        for i from 0 <= i < self._k:
            for j from 0 <= j < (self._n - 1):
                self._latex_debug_string += "$" + permuted_matrix[i, j]._latex_() + "$ & "
            self._latex_debug_string += "$" + permuted_matrix[i, self._n - 1]._latex_() + "$ \\\\\n"

        if comment != "":
            self._latex_debug_string += "\\multicolumn{" + str(self._n) + "}{c}{" + comment + "}\n"
        self._latex_debug_string += "\\end{tabular}};\n"
