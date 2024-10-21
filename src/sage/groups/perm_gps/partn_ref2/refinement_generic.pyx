# sage_setup: distribution = sagemath-gap
r"""
Automorphism groups and canonical labels

For details see section 3 of [Feu2013]_.

Definitions
###########

Let `G` be a group which acts on a finite set `X` and let `\mathcal{L}(G)`
denote the set of subgroups of `G`. We say that
a mapping `Can_G: X \rightarrow X \times G \times \mathcal{L}(G),
x \mapsto \left( CF_G(x), T_G(x), G_x \right)` with

    - `CF_G(gx) = CF_G(x)` for all `x \in X` and all `g\in G` (canonical form)
    - `CF_G(x) = T_G(x) x` for all `x \in X` (transporter element)
    - `G_x = \{g \in G \mid gx=x\}` (stabilizer)

is a canonization algorithm.

Let `H` be another group acting on a set `Y`. A pair `(\theta, \omega)` with
a homomorphism `\theta: G \rightarrow H` and `\omega: X \rightarrow Y`
is called a homomorphism of group actions if `\omega(gx) =
\theta(g)\omega(x)` for all `g \in G`, for all `x \in X`. In the case that
`G=H` and `\theta = id_G`, we call `\omega` a `G`-homomorphism.

A standard
partition is a sequence `C = (C_1, ..., C_k)` of subsets of `[n] = \{0, \ldots,
n-1\}` where the `C_i` are
disjoint, consecutive intervals and their union is equal to `[n]`. Each element
is
called a cell and the elements lying in a cell of cardinality 1 are called
fixed points. Let `I_C` be a sequence of all fixed points of `C` (the exact
ordering is important, but will be determined by the algorithm).
The stabilizer of `(S_n)_C` is a standard Young subgroup of `S_n`.


Requirements
############

In the following we want to present an efficient algorithm for the canonization
problem in the case that the group action is of type `G \rtimes S_n` on
`X^n`. For this goal, we suppose that the group action of `G = G \rtimes \{id\}`
on `X^n` has the following properties:

    -- The action of `G \rtimes \{id\}` on `X^n` is given by a direct product of
       (maybe different) group actions of `G` on `X`.

    -- For every `i \in \{0, \ldots, n-1\}` let `\Pi^{(i)}: X^n \rightarrow X`
       be the projection to the `i`-th coordinate. The function `\Pi^{(i)}`
       is an invariant under the action of the subgroup
       `(S_n)_i = \{\pi \in S_n \mid \pi(i) = i\}`.

Let `(i_0, \ldots, i_{k-1}) = I \subseteq \{0, \ldots, n-1\}` be an
injective sequence. We define
`\Pi^{(I)}(x) := (\Pi^{(i_0)}(x), \ldots, \Pi^{(i_{k-1})}(x))`.
We call an element `x \in X` `I`-semicanonical if and only if
`\Pi^{(I)}(x) \leq \Pi^{(I)}(gx)` for all `g \in G`.

For each `x \in X^n` there is obviously an `I`-semicanonical element in the orbit `Gx`.
Suppose that `x` is already `I`-semicanonical
and `J` is an injective sequence which extends
`I` by one further coordinate `i`, then we call the process of computing
a `J`-semicanonical representative the inner minimization at position `i`.

This could be described by a group action of the stabilizer `G_{\Pi^{(I)}(x)}`.
on `x_i`.


Partitions and refinements
##########################

Now the idea for the canonization algorithm for the action of `G \rtimes
S_n` on `X^n` can be formulated as follows.
Let us suppose that we want
to canonize `x\in X^n`. We build a
backtrack tree in the following manner:

 - each node gets represented by a quadruple `(C, I_C, y, G_{\Pi^{(I_C)}(y))})`
   where `C` is an ordered partition, `I_C` a subsequence of its fixed
   points, `y \in X^n` is `I_C`-semicanonical
   and `G_{\Pi^{(I_C)}(y))}` is the remaining group of the inner action
 - the root node is on level `0` and gets represented by `(([n]), (), x,  G)`
 - nodes with `| I_C | = n`, i.e. all coordinates are fixed, will be leafs
   and we have `G_{\Pi^{(I_C)}(y))} = G_y`
 - the children of some other node `(C, I_C, y, G_{\Pi^{(I_C)}(y))})` are
   constructed by
    * picking one cell `C_j` with `| C_j |` (the choice has to be invariant in
      some sense)
    * separating the point `\min(C_j)` from its cell `C_j` which leads to the partition
      `D` and the sequence `I_D`
    * building the `| C_j |` successors `(D, I_D, z_k,  G_{\Pi^{(I_D)}(z_k))})`
      by applying the permutation
      `\sigma_k := (\min(C_j), k)` for all `k \in C_j` to `y` and computing an
      `I_D`-semicanonical representative `z_k` of `\sigma_k y`
 - if the projection `\Pi^{(I_D)}(z_k)` is not optimal then we stop the
   backtracking in this node (i.e. prune the subtree below this node)
 - the tree is traversed in a depth-first manner
 - the smallest reached leaf node is defined to be the canonical form
   (the defined ordering takes all comparisons in predecessors also into account)
 - equal leaf nodes correspond to automorphism of `x`, they could be used to
   define further pruning methods.

We are able to speed up the computation by making use of refinements
(i.e. `(S_n)_C`-homomorphisms). Suppose we constructed a node
`(C, I_C, y, G_{\Pi^{(I_C)}(y))})`. We define
`Y := \{z \in X^n \mid \Pi^{(I_C)}(z) = \Pi^{(I_C)}(y)\}`
and we search for functions `\omega: Y \rightarrow \ZZ^n` which are
constant on the orbits of `G_{\Pi^{(I_C)}(y)}` and compatible with the action
of `(S_n)_C` on both sides. Then we compute a permutation `\pi \in (S_n)_C`
which cell-wisely sorts the entries of `\omega(y)`. Afterwards we are allowed
to reduce ourselves to the stabilizer `(S_n)_D` of `\omega(\pi y)`. If this
stabilizer contains further fixed points, they are appended (in some fixed
order) to the sequence `I_C` in order to define `I_D`. Furthermore the element
`y` gets updated by an `I_D`-semicanonical representative of this orbit. These
refinements could also be applied iteratively.

As said before, the backtrack search tree is traversed in a depth-first search
manner, hence we maintain a candidate for the result. Each newly constructed
node is compared to this candidate. If it
 - compares larger, then we can prune the
   subtree rooted in this node
 - compares smaller, then we replace the candidate by the actual node
   (more precisely by the next leaf node we will construct, for this we
   use the boolean flag ``_is_candidate_initialized``)
 - compares equal we just continue.


Implementation Details
######################

The class :class:`~sage.groups.perm_gps.partn_ref2.PartitionRefinement_generic`
provides a framework for such a backtracking.
It maintains the partition `C` and the sequence `I_C`.
Instead of permuting the elements `x \in X`
during the construction of nodes and in the refinements,
we use a ``PartitionStack``,
see :mod:`sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label`.
Hence, `C = (C_1, ..., C_k)` and `I_C` will be maintained implicitly. If `\pi \in S_n` is the permutation
applied to this node, then we will store `\pi(I_C)` and in the partition stack
we will store `(\pi(C_1), ..., \pi(C_k))`.

The class further decides when we have to
call the inner minimization and which subtrees could be pruned by the use
of automorphisms, see :class:`sage.groups.perm_gps.partn_ref2.LabelledBranching` for more details.

Derived classes
###############

Derived classes have to implement the following methods, see the
method descriptions for more information:
 - Cython Functions:
     - bint _inner_min_(self, int pos, bint * inner_group_changed)
     - bint _refine(self, bint * part_changed)
     - tuple _store_state_(self)
     - void _restore_state_(self, tuple act_state)
     - void _store_best_(self)
     - void _latex_act_node(self, str comment=None) (to use the implemented debugging method
       via printing the backtrack tree using latex)
     - bint _minimization_allowed_on_col(self, int pos)

 - Python functions:
     - get_canonical_form(self)
     - get_transporter(self)
     - get_autom_gens(self)

AUTHORS:

- Thomas Feulner (2012-11-15): initial version

REFERENCES:

- [Feu2013]_
"""

#*****************************************************************************
#       Copyright (C) 2012 Thomas Feulner <thomas.feulner@uni-bayreuth.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy

from cysignals.memory cimport sig_malloc, sig_realloc, sig_free

from sage.groups.perm_gps.partn_ref.data_structures cimport *
from sage.data_structures.bitset_base cimport *


cdef tuple PS_refinement(PartitionStack * part, long *refine_vals, long *best,
                         int begin, int end,
                         bint * cand_initialized, bint *changed_partition):
    """
    Refine the partition stack by the values given by ``refine_vals``.
    We also compare our actual refinement result with the vector ``best`` in the
    case that ``cand_initialized`` is set to ``True``.

    If the current (canonical) vector is smaller than ``best``, then we will
    replace the values stored in ``best`` and set ``cand_initialized``
    to ``False``.

    The function will set ``changed_partition`` to ``True`` if and only if
    some cell was refined.

    The function will return the boolean value ``False`` if and only if
    the (canonical) vector is larger than ``best``. The second entry of the
    tuple will contain all new singleton cells, which appeared in the
    refinement of ``refine_vals``.
    """
    global global_refine_vals_array

    changed_partition[0] = False
    cdef int i = begin, loc_begin = begin, j, last_val, act_val
    cdef list newly_fixed = []

    while i < end:
        if part.levels[i] <= part.depth:
            # [loc_begin, ..., i] is a block of the old partition
            if i > loc_begin:
                global_refine_vals_array = refine_vals
                qsort(part.entries + loc_begin, (i+1)-loc_begin, sizeof(int), my_comp_func)
            last_val = refine_vals[ part.entries[loc_begin] ]
            j = loc_begin
            while True:
                act_val = refine_vals[ part.entries[j] ]
                if act_val != last_val:
                    if j == 1 or part.levels[j - 2] <= part.depth:
                        # a block of size 1
                        newly_fixed.append(part.entries[j - 1])

                    part.levels[j - 1] = part.depth  # a new block
                    last_val = act_val
                    changed_partition[0] = True

                # compare with the candidate
                if cand_initialized[0]:
                    if act_val < best[j]:
                        cand_initialized[0] = False
                        best[j] = act_val
                    if act_val > best[j]:
                        return (False, [])
                else:
                    best[j] = act_val

                if j == i:  # check for end of while
                    if loc_begin != i and part.levels[i - 1] <= part.depth:
                        # this is a singleton
                        newly_fixed.append(part.entries[i])
                    break

                j += 1

            loc_begin = i + 1
        i += 1
    return (True, newly_fixed)


cdef class _BestValStore:
    r"""
    This class implements a dynamic array of integer vectors of length `n`.
    """
    def __cinit__(self, int n):
        r"""
        Initialization.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import _BestValStore
            sage: B = _BestValStore(5)
        """
        self.default_data_length = n
        self.storage_length = 0
        self.values = NULL

    def __dealloc__(self):
        """
        Dealloc.
        """
        sig_free(self.values)

    cdef long * get_row(self, int i) noexcept:
        r"""
        Return the `i`-th row.

        If this row is not yet initialized (i >= self.storage_length), then
        we extend the array and fill the new positions
        with the all zero vector.
        """
        if i >= self.storage_length:
            self.values = < long *> sig_realloc(self.values, (i + 1)* self.default_data_length * sizeof(long))
            if self.values is NULL:
                raise MemoryError('resizing _BestValStore')
            self.storage_length = i + 1
        return self.values+(i*self.default_data_length)


cdef class LabelledBranching:
    r"""
    This class implements complete labelled branchings.

    To each subgroup of `S_n` we can uniquely assign a directed forest
    on `n` vertices, where the
    edges `(i,j)` fulfill `i<j` and some further conditions on the edge labels
    (which we do not want to state).
    This graph is called a complete labelled branching.

    The edges `(i,j)` will be stored in a vector ``father`` with
    ``father_j = -1`` if `j` is a root of a tree, and
    ``father_j = i`` if `i` is the predecessor of `j`, i.e. is an `(i,j)` is an edge.

    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import LabelledBranching
        sage: L = LabelledBranching(3)
        sage: L.add_gen(libgap.eval('(1,2,3)'))
        sage: L.get_order()
        3
        sage: L.small_generating_set()
        [(1,2,3)]
    """

    def __cinit__(self, n):
        r"""
        Initialization by the identity subgroup.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import LabelledBranching
            sage: L = LabelledBranching(3)
        """
        from sage.libs.gap.libgap import libgap
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup

        self.n = n
        self.group = libgap.eval("Group(())")
        self.ClosureGroup = libgap.eval("ClosureGroup")
        self.father = < int *> sig_malloc(n * sizeof(int))
        self.sym_gp = SymmetricGroup(self.n)
        if self.father is NULL:
            raise MemoryError('allocating LabelledBranching')
        self.act_perm = < int *> sig_malloc(n * sizeof(int))
        if self.act_perm is NULL:
            sig_free(self.father)
            raise MemoryError('allocating LabelledBranching')
        cdef int i
        for i from 0 <= i < self.n:
            self.father[i] = -1

    def __dealloc__(self):
        """
        Dealloc.
        """
        sig_free(self.father)
        sig_free(self.act_perm)

    cpdef add_gen(self, GapElement_Permutation gen):
        r"""
        Add a further generator to the group and
        update the complete labeled branching.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import LabelledBranching
            sage: L = LabelledBranching(3)
            sage: L.add_gen(libgap.eval('(1,2,3)'))
        """
        from sage.libs.gap.libgap import libgap

        self.group = self.ClosureGroup(self.group, gen)
        libgap.StabChain(self.group)
        cdef GapElement h = self.group
        cdef int i, i_1, j, f

        if libgap.IsTrivial(h).sage():
            return
        for i in range(self.n):
            self.father[i] = -1

        while True:
            i = libgap.SmallestMovedPoint(h).sage()
            i_1 = i - 1
            f = self.father[i_1]
            for j in libgap.Orbit(h, i).sage():
                self.father[j - 1] = i_1
            self.father[i_1] = f
            h = libgap.Stabilizer(h, i)
            if libgap.IsTrivial(h).sage():
                break

    cdef bint has_empty_intersection(self, PartitionStack * part) noexcept:
        r"""
        Pruning by automorphisms.

        The partition stack ``part`` encodes a right coset `G\pi` where
        `G <= S_n` is a standard Young subgroup. On the other hand
        this complete labeled branching of the group `A` naturally defines
        a transversal of left coset representatives `T` of `S_n/A`.

        This function returns true if and only if `T \cap G\pi` is empty, see
        [Feu13]_.
        """
        cdef int i, j, k
        for i from 0 <= i < self.n:
            self.act_perm[part.entries[i]] = i

        for i from 0 <= i < self.n:
            j = self.father[i]
            if j != -1:
                j = self.act_perm[j]
                k = self.act_perm[i]

                while k < j:
                    if part.levels[k] <= part.depth:
                        return 1  # True
                    k += 1
        return 0  # False

    def small_generating_set(self):
        r"""
        Return a small set of generators of the group stored by ``self``.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import LabelledBranching
            sage: L = LabelledBranching(3)
            sage: L.small_generating_set()
            [()]
            sage: L.add_gen(libgap.eval('(1,2,3)'))
            sage: L.small_generating_set()
            [(1,2,3)]
        """
        from sage.libs.gap.libgap import libgap
        return libgap.SmallGeneratingSet(self.group).sage(parent=self.sym_gp)

    def get_order(self):
        r"""
        Return the order of the group stored by ``self``.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import LabelledBranching
            sage: L = LabelledBranching(3)
            sage: L.get_order()
            1
            sage: L.add_gen(libgap.eval('(1,2,3)'))
            sage: L.get_order()
            3
        """
        from sage.libs.gap.libgap import libgap

        return libgap.Order(self.group).sage()

cdef class PartitionRefinement_generic:
    r"""
    Implement the partition and refinement framework for
    group actions `G \rtimes S_n` on `X^n` as described in
    :mod:`sage.groups.perm_gps.partn_ref2.refinement_generic`.
    """

    def __cinit__(self, n, *args, **kwds):
        """
        Init.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import PartitionRefinement_generic
            sage: P = PartitionRefinement_generic(5)
        """
        self._n = n  # int
        self._is_candidate_initialized = False

        self._known_automorphisms = LabelledBranching(n)
        self._refine_vals_scratch = <long *> sig_malloc(n * sizeof(long))
        if self._refine_vals_scratch is NULL:
            raise MemoryError('allocating PartitionRefinement_generic')

        self._latex_debug_string = ""
        self._fixed_not_minimized = []
        self._fixed_minimized = []
        self._allowance_best = []
        self._nr_of_inner_min_unmin_calls=0

    def __dealloc__(self):
        """
        Dealloc.
        """
        PS_dealloc(self._part)
        sig_free(self._refine_vals_scratch)

    #####################################################################
    # The following functions have to be implemented by derived classes
    #####################################################################
    cdef bint _inner_min_(self, int pos, bint * inner_group_changed) noexcept:
        """
        Minimize the node by the action of the inner group on the `i`-th position.

        INPUT:

        - ``pos`` -- a position in ``range(self.n)``
        - ``inner_group_changed`` -- will be set to ``True`` if `G_y` got smaller

        OUTPUT: ``True`` if and only if the actual node compares less or equal
        to the candidate for the canonical form
        """
        raise NotImplementedError

    cdef bint _refine(self, bint *part_changed, bint inner_group_changed, bint first_step) noexcept:
        r"""
        Refine the partition ``self._part``.

        Set  ``part_changed`` to ``True`` if and only if the refinement let
        to a smaller subgroup of `S_n`. This function also has to take
        care on ``self._is_candidate_initialized``.

        OUTPUT: ``False`` only if the actual node compares larger than the
        candidate for the canonical form
        """
        raise NotImplementedError

    cdef tuple _store_state_(self):
        r"""
        Store the current state of the node to a tuple, such that it can be
        restored by :meth:`_restore_state_`.
        """
        raise NotImplementedError

    cdef void _restore_state_(self, tuple act_state) noexcept:
        r"""
        The inverse of :meth:`_store_state_`.
        """
        raise NotImplementedError

    cdef void _store_best_(self) noexcept:
        r"""
        Store this leaf node as the actual best candidate for the canonical form.
        """
        raise NotImplementedError

    cdef bint _minimization_allowed_on_col(self, int pos) noexcept:
        r"""
        Decide if we are allowed to perform the inner minimization on position
        ``pos`` which is supposed to be a singleton.
        """
        raise NotImplementedError

    def get_canonical_form(self):
        r"""
        Return the canonical form we have computed.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import PartitionRefinement_generic
            sage: P = PartitionRefinement_generic(5)
            sage: P.get_canonical_form()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def get_transporter(self):
        r"""
        Return the transporter element we have computed.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import PartitionRefinement_generic
            sage: P = PartitionRefinement_generic(5)
            sage: P.get_transporter()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def get_autom_gens(self):
        r"""
        Return a list of generators we have computed.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import PartitionRefinement_generic
            sage: P = PartitionRefinement_generic(5)
            sage: P.get_autom_gens()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def get_autom_order_permutation(self):
        r"""
        Return the order of the automorphism group we have computed.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import PartitionRefinement_generic
            sage: P = PartitionRefinement_generic(5)
            sage: P.get_autom_order_permutation()
            1
        """
        return self._known_automorphisms.get_order()

    ###############################################
    # THE BACKTRACK ALGORITHM:
    ###############################################
    cdef void _init_partition_stack(self, list partition) noexcept:
        r"""
        Initialize the partition stack.

        If the list ``partition`` is not ``None``, it must
        represent a coloring of the coordinates `\{0, \ldots, n-1\}`, i.e.
        a list of disjoint lists.
        """
        cdef bint inner_group_changed = False
        if partition is None:
            self._part = PS_new(self._n, 1)
        else:
            self._part = PS_from_list(partition)
            # minimize singletons which already exist in this partition
            self._fixed_not_minimized = PS_singletons(self._part)
            self._inner_min_unminimized(&inner_group_changed)
        if self._part is NULL:
            sig_free(self._refine_vals_scratch)
            raise MemoryError('initializing the partition stack')

    cdef void _start_Sn_backtrack(self) noexcept:
        r"""
        Start the partition refinement algorithm.
        """
        self._init_latex()
        self._backtrack(True)
        self._finish_latex()

    cdef void _backtrack(self, bint first_step=False) noexcept:
        r"""
        Backtracking with pruning.

         - We first check if the subtree below this node contains any
           valuable information otherwise we return. This test is based on
           the subgroup of already known automorphisms.
         - Perform a refinement. If this node is a leaf, we do all necessary updates.
         - Else, we split the points of a fixed cell and backtrack recursively.
           This leads to a depth-first-search construction of the search tree.
        """
        self._latex_new_lvl()
        cdef bint n_partition_changed = False, inner_group_changed=False

        if self._cut_by_known_automs():
            self._latex_act_node("autcut")
            self._latex_finish_lvl()
            return

        if not self._inner_min_unminimized(&inner_group_changed):
            self._latex_act_node("inner-min")
            self._latex_finish_lvl()
            return

        if not self._refine(&n_partition_changed, inner_group_changed, first_step):
            self._latex_finish_lvl()
            return

        if PS_is_discrete(self._part):
            self._leaf_computations()
            self._latex_finish_lvl()
            return

        # store some important states of the current node
        cdef int old_partition_depth = self._part.depth
        cdef int old_fixed_minimized_len = len(self._fixed_minimized)
        cdef int old_nr_of_inner_min_unmin_calls = self._nr_of_inner_min_unmin_calls
        cdef list old_fixed_not_minimized = copy(self._fixed_not_minimized)
        cdef tuple act_state = self._store_state_()

        # determine a cell C of the partition stack, which we will individualize
        cdef bitset_t b
        bitset_init(b, self._n)
        PS_move_all_mins_to_front(self._part)
        cdef int second_pos
        cdef int smallest = PS_first_smallest_PR(self._part, b, &second_pos, self)
        if second_pos != -1:
            self._fixed_not_minimized.append(second_pos)
        cdef int pos = smallest
        self._latex_act_node()

        # recursively individualize the elements of the cell
        while True:  # run over all elements in this cell
            self._part.depth = old_partition_depth + 1
            PS_clear(self._part)
            PS_split_point(self._part, pos)  # <- individualization!!!
            self._fixed_not_minimized.append(pos)
            self._backtrack()  # backtracking

            #restore the old state and continue backtracking
            self._part.depth = old_partition_depth
            self._fixed_minimized = self._fixed_minimized[:old_fixed_minimized_len]
            self._fixed_not_minimized = copy(old_fixed_not_minimized)
            self._nr_of_inner_min_unmin_calls = old_nr_of_inner_min_unmin_calls
            self._restore_state_(act_state)  # we fixed at least one position, undo this change

            if second_pos != -1:
                # special case when block has size 2
                # update for second run on while, i.e. pos == second_pos
                self._fixed_not_minimized.append(smallest)

            pos = bitset_next(b, pos + 1)  # find the next candidate to individualize
            if pos == -1:
                # all candidates have been individualized!
                break

        bitset_free(b)
        self._latex_finish_lvl()

    cdef bint _inner_min_unminimized(self, bint *inner_group_changed) noexcept:
        r"""
        Compute a subsequence `J = (j_0, \ldots, j_{l-1})`
        of fixed coordinates, which were not yet used in
        the inner minimization process (we call this sequence of used coordinates
        `I`) and compute the `(I+J)`-semicanonical representative.
        """
        cdef bint loc_inner_group_changed, reset_allowance_best = False
        cdef int i, j, best_end, best_ind

        if self._is_candidate_initialized:
            allowance_best = self._allowance_best[self._nr_of_inner_min_unmin_calls]
            best_end = len(allowance_best)
            best_ind = 0
        else:
            allowance_best = []
            reset_allowance_best = True
        inner_group_changed[0] = False
        i = len(self._fixed_not_minimized)-1
        while i>=0:
            pos = self._fixed_not_minimized[i]
            if self._minimization_allowed_on_col(pos):
                for j from 0 <= j < self._n:
                    if self._part.entries[j] == pos:
                        my_final_pos = j
                        break

                if self._is_candidate_initialized:
                    if best_ind == best_end:
                        self._is_candidate_initialized = False
                        reset_allowance_best = True
                    else:
                        if allowance_best[best_ind] < my_final_pos:
                            return False
                        if allowance_best[best_ind] > my_final_pos:
                            self._is_candidate_initialized = False
                            allowance_best = allowance_best[:best_ind]
                            reset_allowance_best = True
                    best_ind += 1 # to avoid to come to this point again

                if not self._is_candidate_initialized:
                    allowance_best.append(my_final_pos)

                if not self._inner_min_(pos, &loc_inner_group_changed):
                    return False

                if not reset_allowance_best and not self._is_candidate_initialized:
                    allowance_best = allowance_best[:best_ind]
                    reset_allowance_best = True

                self._fixed_not_minimized.remove(pos)
                self._fixed_minimized.append(pos)
                inner_group_changed[0] |= loc_inner_group_changed
                if loc_inner_group_changed:
                    i = len(self._fixed_not_minimized)
            i -= 1

        if self._is_candidate_initialized:
            if best_ind < best_end:
                return False
        else:
            self._allowance_best = self._allowance_best[:self._nr_of_inner_min_unmin_calls]
            self._allowance_best.append(allowance_best)
        self._nr_of_inner_min_unmin_calls += 1
        return True

    cdef bint _one_refinement(self, long * best, int begin, int end, bint *inner_group_changed,
                             bint *changed_partition, str refine_name) noexcept:
        r"""
        Let ``self._refine_vals_scratch`` contain the result of the `(S_n)_C`-homomorphism
        `\omega` and ``best`` be the result of `\omega` applied to the candidate.

        This method computes performs the refinements, i.e. updates the partition stack
        and furthermore compares the result with ``best`` and decides
         - if we could prune the subtree -> return value is set to False
         - if the inner group has changed -> sets ``inner_group_changed`` to True
         - if the partition has changed -> sets ``changed_partition`` to True

         The string ``refine_name`` is only necessary for printing within the
         latex output (if activated).
        """
        cdef tuple res = PS_refinement(self._part, self._refine_vals_scratch, best, begin, end,
                          &self._is_candidate_initialized, changed_partition)
        inner_group_changed[0] = False
        if res[0]:
            if changed_partition and self._cut_by_known_automs():
                self._latex_act_node("autcut in " + refine_name)
                return False

            if res[1]:
                self._fixed_not_minimized += res[1]
                if self._inner_min_unminimized(inner_group_changed):
                    return True
                else:
                    self._latex_act_node("inner-min in " + refine_name)
                    return False
            return True
        self._latex_act_node(refine_name)
        return False

    cdef int len(self) noexcept:
        r"""
        Return the degree of the acting symmetric group.
        """
        return self._n

    cdef void _leaf_computations(self) noexcept:
        r"""
        All necessary computations which have to be performed in a leaf.

        There are to possibilities depending on the flag
        ``self._is_candidate_initialized``:
        - ``True`` -- There is another leaf equal to this node. This defines an automorphism
          which we add to the group of known automorphisms stored in
          ``self._known_automorphisms``.
        - ``False`` -- We have shown, that the current leaf is smaller than all
          other leaf nodes already visited. We set it to be the
          candidate for the canonical form.
        """
        from sage.libs.gap.libgap import libgap
        cdef int i

        if self._is_candidate_initialized:
            transp = libgap.PermList([self._part.entries[i] + 1 for i in range(self._n)])
            # transp is the inverse of the actual transporter element (seen as action from left!)
            self._known_automorphisms.add_gen(self._to_best_inverse * transp)  # first apply _to_best_inverse then transp
            self._latex_act_node("AUTOM")
        else:
            self._is_candidate_initialized = True
            self._inner_min_order_best = self._fixed_minimized
            self._to_best = libgap.PermList([self._part.entries[i] + 1 for i in range(self._n)])
            self._to_best_inverse = self._to_best.Inverse()

            # call the method from the derived class
            self._store_best_()
            self._latex_act_node("NEW")

    cdef bint _cut_by_known_automs(self) noexcept:
        r"""
        Return ``True`` if the subtree below this node does not contain
        any *interesting* information.
        """
        return self._is_candidate_initialized and self._known_automorphisms.has_empty_intersection(self._part)

    ###########################################################################
    # These functions are used to produce some latex output:
    # it writes the actual node
    # to a 'tikz node' in latex. You just have to provide '_latex_act_node'
    # in derived classes and to set
    # BACKTRACK_WITHLATEX_DEBUG = 1 in the setup.py script when
    # building this module!
    ###########################################################################
    cdef void _latex_act_node(self, str comment='', int printlvl=0) noexcept:
        r"""
        Append the actual node as a string of latex-commands to
        ``self._latex_debug_string``.
        """
        raise NotImplementedError  # must be implemented by derived classes

    def _latex_view(self, title=None):
        r"""
        Evaluate the latex commands written to ``self._latex_debug_string``
        and view the result.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref2.refinement_generic import PartitionRefinement_generic
            sage: P = PartitionRefinement_generic(5)
            sage: P._latex_view()
            sorry, no debug output was written. Set BACKTRACK_WITHLATEX_DEBUG to True if interested in this information
        """
        if BACKTRACK_WITHLATEX_DEBUG:
            from sage.misc.latex import latex, view
            latex.extra_preamble('')
            latex.add_to_preamble("\\usepackage{tikz}")
            latex.add_to_preamble("\\usepackage{tikz-qtree}")
            view(self._latex_debug_string, engine="pdflatex", title=title)
        else:
            print("sorry, no debug output was written. " +
                  "Set BACKTRACK_WITHLATEX_DEBUG to True if interested in this information")

    cdef void _init_latex(self) noexcept:
        r"""
        Add some initial commands to the string
        ``self._latex_debug_string``.
        """
        if BACKTRACK_WITHLATEX_DEBUG:
            from sage.misc.latex import LatexExpr
            self._latex_debug_string = LatexExpr(
                "\\resizebox{\\textwidth}{!}{\n" +
                "\\begin{tikzpicture}\n" +
                "\\tikzset{level distance=3cm, edge from parent/.style=" +
                "{draw, edge from parent path={(\\tikzparentnode.south) -- (\\tikzchildnode.north)}}}\n" +
                "\\Tree")
            self._latex_debug_string += "[."
            self._latex_act_node()

    cdef void _finish_latex(self) noexcept:
        r"""
        Add some final commands to the string
        ``self._latex_debug_string``.
        """
        if BACKTRACK_WITHLATEX_DEBUG:
            self._latex_debug_string += "]\n"
            self._latex_debug_string += "\\end{tikzpicture} }\n"

    cdef void _latex_new_lvl(self) noexcept:
        r"""
        Add some commands to the string ``self._latex_debug_string``,
        in the case that the depth was increased during the backtracking.
        """
        if BACKTRACK_WITHLATEX_DEBUG:
            self._latex_debug_string += "[."

    cdef void _latex_finish_lvl(self) noexcept:
        r"""
        Add some commands to the string ``self._latex_debug_string``,
        in the case that we return to a node in the backtracking.
        """
        if BACKTRACK_WITHLATEX_DEBUG:
            self._latex_debug_string += "]\n"


cdef int PS_first_smallest_PR(PartitionStack *PS, bitset_t b, int *second_pos=NULL,
                              PartitionRefinement_generic partn_ref_alg=None) noexcept:
    """
    Find the first occurrence of the smallest cell of size greater than one,
    which is admissible (checked by the function ``test_allowance``).
    Its entries are stored to b and its minimum element is returned.

    This generalizes :func:`sage.groups.perm_gps.partn_ref.data_structures.PS_first_smallest`.
    """
    cdef int i = 0, j = 0, location = 0, n = PS.degree
    bitset_zero(b)
    while True:
        if PS.levels[i] <= PS.depth:
            if i != j and n > i - j + 1 and (partn_ref_alg is None or
                                partn_ref_alg._minimization_allowed_on_col(PS.entries[j])):
                n = i - j + 1
                location = j
            j = i + 1
        if PS.levels[i] == -1:
            break
        i += 1
    # location now points to the beginning of the first, smallest,
    # nontrivial cell
    i = location
    while True:
        bitset_flip(b, PS.entries[i])
        if PS.levels[i] <= PS.depth:
            break
        i += 1

    if second_pos != NULL:
        if n == 2:
            second_pos[0] = PS.entries[location + 1]
        else:
            second_pos[0] = -1

    return PS.entries[location]
