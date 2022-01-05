from sage.structure.sage_object cimport SageObject

cdef class Matroid(SageObject):
    cdef public __custom_name
    cdef public _custom_name
    cdef public _cached_info
    cdef int _stored_full_rank
    cdef int _stored_size

    # virtual methods
    cpdef groundset(self)
    cpdef _rank(self, X)

    # internal methods, assuming verified input
    cpdef _max_independent(self, X)
    cpdef _circuit(self, X)
    cpdef _fundamental_circuit(self, B, e)
    cpdef _closure(self, X)
    cpdef _corank(self, X)
    cpdef _max_coindependent(self, X)
    cpdef _cocircuit(self, X)
    cpdef _fundamental_cocircuit(self, B, e)
    cpdef _coclosure(self, X)
    cpdef _augment(self, X, Y)

    cpdef _is_independent(self, X)
    cpdef _is_basis(self, X)
    cpdef _is_circuit(self, X)
    cpdef _is_closed(self, X)
    cpdef _is_coindependent(self, X)
    cpdef _is_cobasis(self, X)
    cpdef _is_cocircuit(self, X)
    cpdef _is_coclosed(self, X)

    cpdef _minor(self, contractions, deletions)
    cpdef _has_minor(self, N, bint certificate=*)
    cpdef _line_length(self, F)
    cpdef _extension(self, element, hyperplanes)

    cdef inline __subset(self, X):
        """
        Convert ``X`` to a ``frozenset`` and check that it is a subset
        of the groundset.

        See ``_subset`` for the corresponding Python method.
        """
        S = frozenset(X)
        if not self.groundset().issuperset(S):
            raise ValueError(f"{X!r} is not a subset of the groundset")
        return S

    cdef inline __subset_all(self, X):
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
    cpdef size(self)

    # matroid oracle
    cpdef rank(self, X=*)
    cpdef full_rank(self)
    cpdef basis(self)
    cpdef max_independent(self, X)
    cpdef circuit(self, X=*)
    cpdef fundamental_circuit(self, B, e)
    cpdef closure(self, X)
    cpdef k_closure(self, X, k)

    cpdef augment(self, X, Y=*)

    cpdef corank(self, X=*)
    cpdef full_corank(self)
    cpdef cobasis(self)
    cpdef max_coindependent(self, X)
    cpdef cocircuit(self, X=*)
    cpdef fundamental_cocircuit(self, B, e)
    cpdef coclosure(self, X)

    cpdef loops(self)
    cpdef is_independent(self, X)
    cpdef is_dependent(self, X)
    cpdef is_basis(self, X)
    cpdef is_circuit(self, X)
    cpdef is_closed(self, X)
    cpdef is_subset_k_closed(self, X, int k)

    cpdef coloops(self)
    cpdef is_coindependent(self, X)
    cpdef is_codependent(self, X)
    cpdef is_cobasis(self, X)
    cpdef is_cocircuit(self, X)
    cpdef is_coclosed(self, X)

    # verification
    cpdef is_valid(self)

    # enumeration
    cpdef circuits(self)
    cpdef nonspanning_circuits(self)
    cpdef cocircuits(self)
    cpdef noncospanning_cocircuits(self)
    cpdef circuit_closures(self)
    cpdef nonspanning_circuit_closures(self)
    cpdef bases(self)
    cpdef independent_sets(self)
    cpdef independent_r_sets(self, long r)
    cpdef nonbases(self)
    cpdef dependent_r_sets(self, long r)
    cpdef _extend_flags(self, flags)
    cpdef _flags(self, r)
    cpdef flats(self, r)
    cpdef coflats(self, r)
    cpdef hyperplanes(self)
    cpdef f_vector(self)
    cpdef broken_circuits(self, ordering=*)
    cpdef no_broken_circuits_sets(self, ordering=*)

    # isomorphism
    cpdef is_isomorphic(self, other, certificate=*)
    cpdef _is_isomorphic(self, other, certificate=*)
    cpdef isomorphism(self, other)
    cpdef _isomorphism(self, other)
    cpdef equals(self, other)
    cpdef is_isomorphism(self, other, morphism)
    cpdef _is_isomorphism(self, other, morphism)

    # minors, dual, truncation
    cpdef minor(self, contractions=*, deletions=*)
    cpdef contract(self, X)
    cpdef delete(self, X)
    cpdef _backslash_(self, X)
    cpdef dual(self)
    cpdef truncation(self)
    cpdef has_minor(self, N, bint certificate=*)
    cpdef has_line_minor(self, k, hyperlines=*, certificate=*)
    cpdef _has_line_minor(self, k, hyperlines, certificate=*)

    # extension
    cpdef extension(self, element=*, subsets=*)
    cpdef coextension(self, element=*, subsets=*)
    cpdef modular_cut(self, subsets)
    cpdef linear_subclasses(self, line_length=*, subsets=*)
    cpdef extensions(self, element=*, line_length=*, subsets=*)

    # connectivity
    cpdef simplify(self)
    cpdef cosimplify(self)
    cpdef is_simple(self)
    cpdef is_cosimple(self)
    cpdef components(self)
    cpdef is_connected(self, certificate=*)
    cpdef connectivity(self, S, T=*)
    cpdef _connectivity(self, S, T)
    cpdef is_kconnected(self, k, certificate=*)
    cpdef link(self, S, T)
    cpdef _link(self, S, T)
    cpdef _is_3connected_shifting(self, certificate=*)
    cpdef _is_4connected_shifting(self, certificate=*)
    cpdef _shifting_all(self, X, P_rows, P_cols, Q_rows, Q_cols, m)
    cpdef _shifting(self, X, X_1, Y_2, X_2, Y_1, m)
    cpdef is_3connected(self, certificate=*, algorithm=*)
    cpdef is_4connected(self, certificate=*, algorithm=*)
    cpdef _is_3connected_CE(self, certificate=*)
    cpdef _is_3connected_BC(self, certificate=*)
    cpdef _is_3connected_BC_recursion(self, basis, fund_cocircuits)

    # representability
    cpdef _local_binary_matroid(self, basis=*)
    cpdef binary_matroid(self, randomized_tests=*, verify=*)
    cpdef is_binary(self, randomized_tests=*)
    cpdef _local_ternary_matroid(self, basis=*)
    cpdef ternary_matroid(self, randomized_tests=*, verify=*)
    cpdef is_ternary(self, randomized_tests=*)

    # matroid k-closed
    cpdef is_k_closed(self, int k)

    # matroid chordality
    cpdef _is_circuit_chordal(self, frozenset C, bint certificate=*)
    cpdef is_circuit_chordal(self, C, bint certificate=*)
    cpdef is_chordal(self, k1=*, k2=*, bint certificate=*)
    cpdef chordality(self)

    # optimization
    cpdef max_weight_independent(self, X=*, weights=*)
    cpdef max_weight_coindependent(self, X=*, weights=*)
    cpdef is_max_weight_independent_generic(self, X=*, weights=*)
    cpdef is_max_weight_coindependent_generic(self, X=*, weights=*)
    cpdef intersection(self, other, weights=*)
    cpdef _intersection(self, other, weights)
    cpdef _intersection_augmentation(self, other, weights, Y)
    cpdef intersection_unweighted(self, other)
    cpdef _intersection_unweighted(self, other)
    cpdef _intersection_augmentation_unweighted(self, other, Y)
    cpdef partition(self)

    # invariants
    cpdef _internal(self, B)
    cpdef _external(self, B)
    cpdef tutte_polynomial(self, x=*, y=*)
    cpdef flat_cover(self, solver=*, verbose=*, integrality_tolerance=*)

    # visualization
    cpdef plot(self,B=*,lineorders=*,pos_method=*,pos_dict=*,save_pos=*)
    cpdef show(self,B=*,lineorders=*,pos_method=*,pos_dict=*,save_pos=*,lims=*)
    cpdef _fix_positions(self,pos_dict=*,lineorders=*)
