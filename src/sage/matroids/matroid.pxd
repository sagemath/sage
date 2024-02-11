from sage.structure.sage_object cimport SageObject

cdef class Matroid(SageObject):
    cdef public _SageObject__custom_name
    cdef public _cached_info
    cdef int _stored_full_rank
    cdef int _stored_size

    # virtual methods
    cpdef groundset(self) noexcept
    cpdef _rank(self, X) noexcept

    # internal methods, assuming verified input
    cpdef _max_independent(self, X) noexcept
    cpdef _circuit(self, X) noexcept
    cpdef _fundamental_circuit(self, B, e) noexcept
    cpdef _closure(self, X) noexcept
    cpdef _corank(self, X) noexcept
    cpdef _max_coindependent(self, X) noexcept
    cpdef _cocircuit(self, X) noexcept
    cpdef _fundamental_cocircuit(self, B, e) noexcept
    cpdef _coclosure(self, X) noexcept
    cpdef _augment(self, X, Y) noexcept

    cpdef _is_independent(self, X) noexcept
    cpdef _is_basis(self, X) noexcept
    cpdef _is_circuit(self, X) noexcept
    cpdef _is_closed(self, X) noexcept
    cpdef _is_coindependent(self, X) noexcept
    cpdef _is_cobasis(self, X) noexcept
    cpdef _is_cocircuit(self, X) noexcept
    cpdef _is_coclosed(self, X) noexcept

    cpdef _minor(self, contractions, deletions) noexcept
    cpdef _has_minor(self, N, bint certificate=*) noexcept
    cpdef _line_length(self, F) noexcept
    cpdef _extension(self, element, hyperplanes) noexcept

    cdef inline _subset_internal(self, X) noexcept:
        """
        Convert ``X`` to a ``frozenset`` and check that it is a subset
        of the groundset.

        See ``_subset`` for the corresponding Python method.
        """
        S = frozenset(X)
        if not self.groundset().issuperset(S):
            raise ValueError(f"{X!r} is not a subset of the groundset")
        return S

    cdef inline __subset_all(self, X) noexcept:
        """
        If ``X`` is ``None``, return the groundset.

        Otherwise, do like ``_subset``:
        convert ``X`` to a ``frozenset`` and check that it is a subset
        of the groundset.

        See ``_subset_all`` for the corresponding Python method.
        """
        if X is None:
            return self.groundset()
        S = frozenset(X)
        if not self.groundset().issuperset(S):
            raise ValueError(f"{X!r} is not a subset of the groundset")
        return S

    # ** user-facing methods **
    cpdef size(self) noexcept

    # matroid oracle
    cpdef rank(self, X=*) noexcept
    cpdef full_rank(self) noexcept
    cpdef basis(self) noexcept
    cpdef max_independent(self, X) noexcept
    cpdef circuit(self, X=*) noexcept
    cpdef fundamental_circuit(self, B, e) noexcept
    cpdef closure(self, X) noexcept
    cpdef k_closure(self, X, k) noexcept

    cpdef augment(self, X, Y=*) noexcept

    cpdef corank(self, X=*) noexcept
    cpdef full_corank(self) noexcept
    cpdef cobasis(self) noexcept
    cpdef max_coindependent(self, X) noexcept
    cpdef cocircuit(self, X=*) noexcept
    cpdef fundamental_cocircuit(self, B, e) noexcept
    cpdef coclosure(self, X) noexcept

    cpdef loops(self) noexcept
    cpdef is_independent(self, X) noexcept
    cpdef is_dependent(self, X) noexcept
    cpdef is_basis(self, X) noexcept
    cpdef is_circuit(self, X) noexcept
    cpdef is_closed(self, X) noexcept
    cpdef is_subset_k_closed(self, X, int k) noexcept

    cpdef coloops(self) noexcept
    cpdef is_coindependent(self, X) noexcept
    cpdef is_codependent(self, X) noexcept
    cpdef is_cobasis(self, X) noexcept
    cpdef is_cocircuit(self, X) noexcept
    cpdef is_coclosed(self, X) noexcept

    # verification
    cpdef is_valid(self) noexcept

    # enumeration
    cpdef circuits(self) noexcept
    cpdef nonspanning_circuits(self) noexcept
    cpdef cocircuits(self) noexcept
    cpdef noncospanning_cocircuits(self) noexcept
    cpdef circuit_closures(self) noexcept
    cpdef nonspanning_circuit_closures(self) noexcept
    cpdef bases(self) noexcept
    cpdef independent_sets(self) noexcept
    cpdef independent_r_sets(self, long r) noexcept
    cpdef nonbases(self) noexcept
    cpdef dependent_r_sets(self, long r) noexcept
    cpdef _extend_flags(self, flags) noexcept
    cpdef _flags(self, r) noexcept
    cpdef flats(self, r) noexcept
    cpdef coflats(self, r) noexcept
    cpdef hyperplanes(self) noexcept
    cpdef f_vector(self) noexcept
    cpdef broken_circuits(self, ordering=*) noexcept
    cpdef no_broken_circuits_sets(self, ordering=*) noexcept

    # isomorphism
    cpdef is_isomorphic(self, other, certificate=*) noexcept
    cpdef _is_isomorphic(self, other, certificate=*) noexcept
    cpdef isomorphism(self, other) noexcept
    cpdef _isomorphism(self, other) noexcept
    cpdef equals(self, other) noexcept
    cpdef is_isomorphism(self, other, morphism) noexcept
    cpdef _is_isomorphism(self, other, morphism) noexcept

    # minors, dual, truncation
    cpdef minor(self, contractions=*, deletions=*) noexcept
    cpdef contract(self, X) noexcept
    cpdef delete(self, X) noexcept
    cpdef _backslash_(self, X) noexcept
    cpdef dual(self) noexcept
    cpdef truncation(self) noexcept
    cpdef has_minor(self, N, bint certificate=*) noexcept
    cpdef has_line_minor(self, k, hyperlines=*, certificate=*) noexcept
    cpdef _has_line_minor(self, k, hyperlines, certificate=*) noexcept

    # extension
    cpdef extension(self, element=*, subsets=*) noexcept
    cpdef coextension(self, element=*, subsets=*) noexcept
    cpdef modular_cut(self, subsets) noexcept
    cpdef linear_subclasses(self, line_length=*, subsets=*) noexcept
    cpdef extensions(self, element=*, line_length=*, subsets=*) noexcept

    # connectivity
    cpdef simplify(self) noexcept
    cpdef cosimplify(self) noexcept
    cpdef is_simple(self) noexcept
    cpdef is_cosimple(self) noexcept
    cpdef components(self) noexcept
    cpdef is_connected(self, certificate=*) noexcept
    cpdef connectivity(self, S, T=*) noexcept
    cpdef _connectivity(self, S, T) noexcept
    cpdef is_kconnected(self, k, certificate=*) noexcept
    cpdef link(self, S, T) noexcept
    cpdef _link(self, S, T) noexcept
    cpdef _is_3connected_shifting(self, certificate=*) noexcept
    cpdef _is_4connected_shifting(self, certificate=*) noexcept
    cpdef _shifting_all(self, X, P_rows, P_cols, Q_rows, Q_cols, m) noexcept
    cpdef _shifting(self, X, X_1, Y_2, X_2, Y_1, m) noexcept
    cpdef is_3connected(self, certificate=*, algorithm=*) noexcept
    cpdef is_4connected(self, certificate=*, algorithm=*) noexcept
    cpdef _is_3connected_CE(self, certificate=*) noexcept
    cpdef _is_3connected_BC(self, certificate=*) noexcept
    cpdef _is_3connected_BC_recursion(self, basis, fund_cocircuits) noexcept
    cpdef is_paving(self) noexcept
    cpdef is_sparse_paving(self) noexcept
    cpdef girth(self) noexcept

    # representability
    cpdef _local_binary_matroid(self, basis=*) noexcept
    cpdef binary_matroid(self, randomized_tests=*, verify=*) noexcept
    cpdef is_binary(self, randomized_tests=*) noexcept
    cpdef _local_ternary_matroid(self, basis=*) noexcept
    cpdef ternary_matroid(self, randomized_tests=*, verify=*) noexcept
    cpdef is_ternary(self, randomized_tests=*) noexcept

    # matroid k-closed
    cpdef is_k_closed(self, int k) noexcept

    # matroid chordality
    cpdef _is_circuit_chordal(self, frozenset C, bint certificate=*) noexcept
    cpdef is_circuit_chordal(self, C, bint certificate=*) noexcept
    cpdef is_chordal(self, k1=*, k2=*, bint certificate=*) noexcept
    cpdef chordality(self) noexcept

    # optimization
    cpdef max_weight_independent(self, X=*, weights=*) noexcept
    cpdef max_weight_coindependent(self, X=*, weights=*) noexcept
    cpdef is_max_weight_independent_generic(self, X=*, weights=*) noexcept
    cpdef is_max_weight_coindependent_generic(self, X=*, weights=*) noexcept
    cpdef intersection(self, other, weights=*) noexcept
    cpdef _intersection(self, other, weights) noexcept
    cpdef _intersection_augmentation(self, other, weights, Y) noexcept
    cpdef intersection_unweighted(self, other) noexcept
    cpdef _intersection_unweighted(self, other) noexcept
    cpdef _intersection_augmentation_unweighted(self, other, Y) noexcept
    cpdef partition(self) noexcept

    # invariants
    cpdef _internal(self, B) noexcept
    cpdef _external(self, B) noexcept
    cpdef tutte_polynomial(self, x=*, y=*) noexcept
    cpdef flat_cover(self, solver=*, verbose=*, integrality_tolerance=*) noexcept

    # misc
    cpdef automorphism_group(self) noexcept
    cpdef bergman_complex(self) noexcept
    cpdef augmented_bergman_complex(self) noexcept

    # visualization
    cpdef plot(self,B=*,lineorders=*,pos_method=*,pos_dict=*,save_pos=*) noexcept
    cpdef show(self,B=*,lineorders=*,pos_method=*,pos_dict=*,save_pos=*,lims=*) noexcept
    cpdef _fix_positions(self,pos_dict=*,lineorders=*) noexcept
