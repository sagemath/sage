from sage.data_structures.bitset cimport bitset_t

from sage.matroids.matroid cimport Matroid
from sage.matroids.basis_exchange_matroid cimport BasisExchangeMatroid
from sage.matroids.lean_matrix cimport LeanMatrix, GenericMatrix, BinaryMatrix, TernaryMatrix, QuaternaryMatrix


cdef class LinearMatroid(BasisExchangeMatroid):
    cdef LeanMatrix _A, _representation
    cdef long *_prow
    cdef object _zero, _one

    cpdef _forget(self) noexcept
    cpdef base_ring(self) noexcept
    cpdef characteristic(self) noexcept

    cdef list _setup_internal_representation(self, matrix, reduced_matrix, ring, keep_initial_representation) noexcept
    cdef _exchange_value_internal(self, long x, long y) noexcept

    cpdef representation(self, B=*, reduced=*, labels=*, order=*, lift_map=*) noexcept
    cpdef _current_rows_cols(self, B=*) noexcept
    cpdef representation_vectors(self) noexcept
    cpdef LeanMatrix _basic_representation(self, B=*) noexcept
    cpdef LeanMatrix _reduced_representation(self, B=*) noexcept

    cpdef bint _is_field_isomorphism(self, LinearMatroid other, morphism) noexcept
    cpdef is_field_equivalent(self, other) noexcept
    cpdef is_field_isomorphism(self, other, morphism) noexcept
    # cpdef is_field_isomorphic(self, other)  # TODO: currently only works as ``def``
    cpdef _fast_isom_test(self, other) noexcept

    cpdef _minor(self, contractions, deletions) noexcept
    cpdef dual(self) noexcept
    cpdef has_line_minor(self, k, hyperlines=*, certificate=*) noexcept
    cpdef has_field_minor(self, N) noexcept

    cpdef _exchange_value(self, e, f) noexcept
    cpdef fundamental_cycle(self, B, e) noexcept
    cpdef fundamental_cocycle(self, B, e) noexcept

    cpdef _line_ratios(self, F) noexcept
    cpdef _line_length(self, F) noexcept

    cpdef _line_cross_ratios(self, F) noexcept
    cpdef cross_ratios(self, hyperlines=*) noexcept
    cpdef cross_ratio(self, F, a, b, c, d) noexcept
    cpdef _line_cross_ratio_test(self, F, x, fundamentals) noexcept
    cpdef _cross_ratio_test(self, x, fundamentals, hyperlines=*) noexcept

    cpdef linear_extension(self, element, chain=*, col=*) noexcept
    cpdef linear_coextension(self, element, cochain=*, row=*) noexcept
    cpdef _linear_extensions(self, element, chains) noexcept
    cpdef _linear_coextensions(self, element, cochains) noexcept
    cdef _extend_chains(self, C, f, fundamentals=*) noexcept
    cpdef _linear_extension_chains(self, F, fundamentals=*) noexcept
    cpdef linear_extension_chains(self, F=*, simple=*, fundamentals=*) noexcept
    cpdef linear_coextension_cochains(self, F=*, cosimple=*, fundamentals=*) noexcept
    cpdef linear_extensions(self, element=*, F=*, simple=*, fundamentals=*) noexcept
    cpdef linear_coextensions(self, element=*, F=*, cosimple=*, fundamentals=*) noexcept

    cpdef _is_3connected_shifting(self, certificate=*) noexcept
    cpdef _is_4connected_shifting(self, certificate=*) noexcept

    cpdef is_valid(self) noexcept

cdef class BinaryMatroid(LinearMatroid):
    cdef tuple _b_invariant, _b_partition
    cdef BinaryMatrix _b_projection, _eq_part

    cpdef base_ring(self) noexcept
    cpdef characteristic(self) noexcept

    cpdef _current_rows_cols(self, B=*) noexcept
    cpdef LeanMatrix _basic_representation(self, B=*) noexcept
    cpdef LeanMatrix _reduced_representation(self, B=*) noexcept

    cdef  __fundamental_cocircuit(self, bitset_t, long x) noexcept

    cpdef _is_isomorphic(self, other, certificate=*) noexcept

    cpdef _minor(self, contractions, deletions) noexcept

    cpdef _make_invariant(self) noexcept
    cpdef _invariant(self) noexcept
    cpdef bicycle_dimension(self) noexcept
    cpdef brown_invariant(self) noexcept
    cpdef _principal_tripartition(self) noexcept
    cpdef BinaryMatrix _projection(self) noexcept
    cpdef BinaryMatrix _projection_partition(self) noexcept
    cpdef _fast_isom_test(self, other) noexcept

    cpdef is_graphic(self) noexcept
    cpdef is_valid(self) noexcept


cdef class TernaryMatroid(LinearMatroid):
    cdef object _two
    cdef tuple _t_invariant, _t_partition
    cdef TernaryMatrix _t_projection

    cpdef base_ring(self) noexcept
    cpdef characteristic(self) noexcept

    cpdef _current_rows_cols(self, B=*) noexcept
    cpdef LeanMatrix _basic_representation(self, B=*) noexcept
    cpdef LeanMatrix _reduced_representation(self, B=*) noexcept

    cdef  __fundamental_cocircuit(self, bitset_t, long x) noexcept

    cpdef _is_isomorphic(self, other, certificate=*) noexcept

    cpdef _minor(self, contractions, deletions) noexcept

    cpdef _make_invariant(self) noexcept
    cpdef _invariant(self) noexcept
    cpdef bicycle_dimension(self) noexcept
    cpdef character(self) noexcept
    cpdef _principal_quadripartition(self) noexcept
    cpdef TernaryMatrix _projection(self) noexcept
    cpdef _fast_isom_test(self, other) noexcept

    cpdef is_valid(self) noexcept

cdef class QuaternaryMatroid(LinearMatroid):
    cdef object _x_zero, _x_one
    cdef tuple _q_invariant, _q_partition
    cdef QuaternaryMatrix _q_projection

    cpdef base_ring(self) noexcept
    cpdef characteristic(self) noexcept

    cpdef _current_rows_cols(self, B=*) noexcept
    cpdef LeanMatrix _basic_representation(self, B=*) noexcept
    cpdef LeanMatrix _reduced_representation(self, B=*) noexcept

    cdef  __fundamental_cocircuit(self, bitset_t, long x) noexcept

    cpdef _is_isomorphic(self, other, certificate=*) noexcept

    cpdef _minor(self, contractions, deletions) noexcept

    cpdef _make_invariant(self) noexcept
    cpdef _invariant(self) noexcept
    cpdef bicycle_dimension(self) noexcept
    cpdef _principal_tripartition(self) noexcept
    cpdef _fast_isom_test(self, other) noexcept

    cpdef is_valid(self) noexcept

cdef class RegularMatroid(LinearMatroid):
    cdef _bases_count, _r_invariant
    cdef _r_projection, _r_hypergraph
    cdef _hypergraph_vertex_partition, _hypergraph_tuples

    cpdef base_ring(self) noexcept
    cpdef characteristic(self) noexcept

    cpdef _is_isomorphic(self, other, certificate=*) noexcept

    cpdef _invariant(self) noexcept
    cpdef _fast_isom_test(self, other) noexcept

    cpdef bases_count(self) noexcept
    cpdef _projection(self) noexcept
    cpdef _hypergraph(self) noexcept
    cdef _hypertest(self, other) noexcept
    cpdef has_line_minor(self, k, hyperlines=*, certificate=*) noexcept
    cpdef _linear_extension_chains(self, F, fundamentals=*) noexcept

    cpdef is_graphic(self) noexcept
    cpdef is_valid(self) noexcept
