from sage.structure.sage_object cimport SageObject
from sage.matroids.set_system cimport SetSystem

cdef class Matroid(SageObject):
    cdef public _SageObject__custom_name
    cdef public _cached_info
    cdef int _stored_full_rank
    cdef int _stored_size

    # virtual methods
    cpdef frozenset groundset(self)
    cpdef int _rank(self, frozenset X) except? -1

    # internal methods, assuming verified input
    cpdef frozenset _max_independent(self, frozenset X)
    cpdef frozenset _circuit(self, frozenset X)
    cpdef frozenset _fundamental_circuit(self, frozenset B, e)
    cpdef frozenset _closure(self, frozenset X)
    cpdef int _corank(self, frozenset X) noexcept
    cpdef frozenset _max_coindependent(self, frozenset X)
    cpdef frozenset _cocircuit(self, frozenset X)
    cpdef frozenset _fundamental_cocircuit(self, frozenset B, e)
    cpdef frozenset _coclosure(self, frozenset X)
    cpdef frozenset _augment(self, frozenset X, frozenset Y)

    cpdef bint _is_independent(self, frozenset X) noexcept
    cpdef bint _is_basis(self, frozenset X) noexcept
    cpdef bint _is_circuit(self, frozenset X) noexcept
    cpdef bint _is_closed(self, frozenset X) noexcept
    cpdef bint _is_coindependent(self, frozenset X) noexcept
    cpdef bint _is_cobasis(self, frozenset X) noexcept
    cpdef bint _is_cocircuit(self, frozenset X) noexcept
    cpdef bint _is_coclosed(self, frozenset X) noexcept

    cpdef _minor(self, contractions, deletions)
    cpdef _has_minor(self, N, bint certificate=*)
    cpdef _line_length(self, F)
    cpdef _extension(self, element, hyperplanes)

    cdef inline _subset_internal(self, X):
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
    cpdef bint is_valid(self) noexcept

    # enumeration
    cpdef SetSystem circuits(self, k=*)
    cpdef SetSystem nonspanning_circuits(self)
    cpdef SetSystem cocircuits(self)
    cpdef SetSystem noncospanning_cocircuits(self)
    cpdef dict circuit_closures(self)
    cpdef dict nonspanning_circuit_closures(self)
    cpdef SetSystem bases(self)
    cpdef SetSystem independent_sets(self, long k=*)
    cdef SetSystem _independent_sets(self)
    cpdef SetSystem nonbases(self)
    cpdef SetSystem dependent_sets(self, long k)
    cpdef list _extend_flags(self, list flags)
    cpdef list _flags(self, long k)
    cpdef SetSystem flats(self, long k)
    cpdef SetSystem coflats(self, long k)
    cpdef SetSystem hyperplanes(self)
    cpdef list f_vector(self)
    cpdef list whitney_numbers(self)
    cpdef list whitney_numbers2(self)
    cpdef SetSystem broken_circuits(self, ordering=*)
    cpdef SetSystem no_broken_circuits_sets(self, ordering=*)

    # polytopes
    cpdef matroid_polytope(self)
    cpdef independence_matroid_polytope(self)

    # isomorphism
    cpdef is_isomorphic(self, other, certificate=*)
    cpdef _is_isomorphic(self, other, certificate=*)
    cpdef isomorphism(self, other)
    cpdef _isomorphism(self, other)
    cpdef equals(self, other)
    cpdef is_isomorphism(self, other, morphism)
    cpdef _is_isomorphism(self, other, morphism)
    cpdef _relabel_map(self, mapping)

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
    cpdef coextensions(self, element=*, coline_length=*, subsets=*)

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
    cpdef bint is_paving(self) noexcept
    cpdef bint is_sparse_paving(self) noexcept
    cpdef girth(self)

    # representability
    cpdef _local_binary_matroid(self, basis=*)
    cpdef binary_matroid(self, randomized_tests=*, verify=*)
    cpdef is_binary(self, randomized_tests=*)
    cpdef _local_ternary_matroid(self, basis=*)
    cpdef ternary_matroid(self, randomized_tests=*, verify=*)
    cpdef is_ternary(self, randomized_tests=*)
    cpdef bint is_regular(self) noexcept
    cpdef bint is_graphic(self) noexcept

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
    cpdef characteristic_polynomial(self, la=*)
    cpdef flat_cover(self, solver=*, verbose=*, integrality_tolerance=*)

    # misc
    cpdef automorphism_group(self)
    cpdef bergman_complex(self)
    cpdef augmented_bergman_complex(self)
    cpdef broken_circuit_complex(self, ordering=*)

    # visualization
    cpdef plot(self, B=*, lineorders=*, pos_method=*, pos_dict=*, save_pos=*)
    cpdef show(self, B=*, lineorders=*, pos_method=*, pos_dict=*, save_pos=*, lims=*)
    cpdef _fix_positions(self, pos_dict=*, lineorders=*)

    # construction
    cpdef direct_sum(self, matroids)
