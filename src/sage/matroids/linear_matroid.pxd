from sage.misc.bitset cimport bitset_t

from matroid cimport Matroid
from basis_exchange_matroid cimport BasisExchangeMatroid
from sage.matrix.matrix2 cimport Matrix
from lean_matrix cimport LeanMatrix, GenericMatrix, BinaryMatrix, TernaryMatrix, QuaternaryMatrix

cdef inline gauss_jordan_reduce(LeanMatrix A, columns)
cdef inline characteristic(LeanMatrix A)

cdef class LinearMatroid(BasisExchangeMatroid):
    cdef LeanMatrix _A, _representation
    cdef long *_prow
    cdef object _zero, _one

    cpdef _forget(self)
    cpdef base_ring(self)
    cpdef characteristic(self)

    cdef list _setup_internal_representation(self, matrix, reduced_matrix, ring, keep_initial_representation)
    cdef  bint __is_exchange_pair(self, long x, long y)
    cdef  bint __exchange(self, long x, long y)
    cdef  __exchange_value(self, long x, long y)

    cpdef representation(self, B=*, reduced=*, labels=*, order=*)
    cpdef _current_rows_cols(self, B=*)
    cpdef representation_vectors(self)
    cpdef LeanMatrix _basic_representation(self, B=*)
    cpdef LeanMatrix _reduced_representation(self, B=*)

    cpdef bint _is_field_isomorphism(self, LinearMatroid other, morphism)
    cpdef is_field_equivalent(self, other)
    cpdef is_field_isomorphism(self, other, morphism)
    # cpdef is_field_isomorphic(self, other)  # TODO: currently only works as ``def``
    cpdef _fast_isom_test(self, other)

    cpdef _minor(self, contractions, deletions)
    cpdef dual(self)
    cpdef has_line_minor(self, k, hyperlines=*)
    cpdef has_field_minor(self, N)

    cpdef _exchange_value(self, e, f)
    cpdef fundamental_cycle(self, B, e)
    cpdef fundamental_cocycle(self, B, e)

    cpdef _line_ratios(self, F)
    cpdef _line_length(self, F)

    cpdef _line_cross_ratios(self, F)
    cpdef cross_ratios(self, hyperlines=*)
    cpdef cross_ratio(self, F, a, b, c, d)
    cpdef _line_cross_ratio_test(self, F, x, fundamentals)
    cpdef _cross_ratio_test(self, x, fundamentals, hyperlines=*)

    cpdef linear_extension(self, element, chain=*, col=*)
    cpdef linear_coextension(self, element, cochain=*, row=*)
    cpdef _linear_extensions(self, element, chains)
    cpdef _linear_coextensions(self, element, cochains)
    cdef _extend_chains(self, C, f, fundamentals=*)
    cpdef _linear_extension_chains(self, F, fundamentals=*)
    cpdef linear_extension_chains(self, F=*, simple=*, fundamentals=*)
    cpdef linear_coextension_cochains(self, F=*, cosimple=*, fundamentals=*)
    cpdef linear_extensions(self, element=*, F=*, simple=*, fundamentals=*)
    cpdef linear_coextensions(self, element=*, F=*, cosimple=*, fundamentals=*)

    cpdef is_valid(self)

cdef class BinaryMatroid(LinearMatroid):
    cdef tuple _b_invariant, _b_partition
    cdef BinaryMatrix _b_projection, _eq_part

    cpdef base_ring(self)
    cpdef characteristic(self)

    cdef  bint __is_exchange_pair(self, long x, long y)
    cdef  bint __exchange(self, long x, long y)
    cdef  __exchange_value(self, long x, long y)

    cpdef _current_rows_cols(self, B=*)
    cpdef LeanMatrix _basic_representation(self, B=*)
    cpdef LeanMatrix _reduced_representation(self, B=*)

    cdef  __fundamental_cocircuit(self, bitset_t, long x)

    cpdef _is_isomorphic(self, other)

    cpdef _minor(self, contractions, deletions)

    cpdef _make_invariant(self)
    cpdef _invariant(self)
    cpdef bicycle_dimension(self)
    cpdef brown_invariant(self)
    cpdef _principal_tripartition(self)
    cpdef BinaryMatrix _projection(self)
    cpdef BinaryMatrix _projection_partition(self)
    cpdef _fast_isom_test(self, other)

    cpdef is_graphic(self)
    cpdef is_valid(self)


cdef class TernaryMatroid(LinearMatroid):
    cdef object _two
    cdef tuple _t_invariant, _t_partition
    cdef TernaryMatrix _t_projection

    cpdef base_ring(self)
    cpdef characteristic(self)

    cdef  bint __is_exchange_pair(self, long x, long y)
    cdef  bint __exchange(self, long x, long y)
    cdef  __exchange_value(self, long x, long y)

    cpdef _current_rows_cols(self, B=*)
    cpdef LeanMatrix _basic_representation(self, B=*)
    cpdef LeanMatrix _reduced_representation(self, B=*)

    cdef  __fundamental_cocircuit(self, bitset_t, long x)

    cpdef _is_isomorphic(self, other)

    cpdef _minor(self, contractions, deletions)

    cpdef _make_invariant(self)
    cpdef _invariant(self)
    cpdef bicycle_dimension(self)
    cpdef character(self)
    cpdef _principal_quadripartition(self)
    cpdef TernaryMatrix _projection(self)
    cpdef _fast_isom_test(self, other)

    cpdef is_valid(self)

cdef class QuaternaryMatroid(LinearMatroid):
    cdef object _x_zero, _x_one
    cdef tuple _q_invariant, _q_partition
    cdef QuaternaryMatrix _q_projection

    cpdef base_ring(self)
    cpdef characteristic(self)

    cdef  bint __is_exchange_pair(self, long x, long y)
    cdef  bint __exchange(self, long x, long y)
    cdef  __exchange_value(self, long x, long y)

    cpdef _current_rows_cols(self, B=*)
    cpdef LeanMatrix _basic_representation(self, B=*)
    cpdef LeanMatrix _reduced_representation(self, B=*)

    cdef  __fundamental_cocircuit(self, bitset_t, long x)

    cpdef _is_isomorphic(self, other)

    cpdef _minor(self, contractions, deletions)

    cpdef _make_invariant(self)
    cpdef _invariant(self)
    cpdef bicycle_dimension(self)
    cpdef _principal_tripartition(self)
    cpdef _fast_isom_test(self, other)

    cpdef is_valid(self)

cdef class RegularMatroid(LinearMatroid):
    cdef _bases_count, _r_invariant
    cdef _r_projection, _r_hypergraph
    cdef _hypergraph_vertex_partition, _hypergraph_tuples

    cpdef base_ring(self)
    cpdef characteristic(self)

    cdef  bint __is_exchange_pair(self, long x, long y)
    cdef  bint __exchange(self, long x, long y)
    cdef  __exchange_value(self, long x, long y)

    cpdef _is_isomorphic(self, other)

    cpdef _invariant(self)
    cpdef _fast_isom_test(self, other)

    cpdef bases_count(self)
    cpdef _projection(self)
    cpdef _hypergraph(self)
    cdef _hypertest(self, other)
    cpdef has_line_minor(self, k, hyperlines=*)
    cpdef _linear_extension_chains(self, F, fundamentals=*)

    cpdef is_graphic(self)
    cpdef is_valid(self)
