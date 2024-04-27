from sage.data_structures.bitset cimport bitset_t
from sage.libs.gmp.types cimport mpq_t

cdef class LeanMatrix:
    cdef long _nrows
    cdef long _ncols

    cdef LeanMatrix copy(self)   # Deprecated Sage matrix operation
    cdef int resize(self, long k) except -1   # Not a Sage matrix operation
    cdef LeanMatrix stack(self, LeanMatrix M)
    cdef LeanMatrix augment(self, LeanMatrix M)
    cdef LeanMatrix prepend_identity(self)   # Not a Sage matrix operation

    cpdef long ncols(self) except -1
    cpdef long nrows(self) except -1
    cpdef base_ring(self)
    cpdef characteristic(self)   # Not a Sage matrix operation

    cdef get_unsafe(self, long r, long c)
    cdef int set_unsafe(self, long r, long c, x) except -1
    cdef bint is_nonzero(self, long r, long c) except -2  # Not a Sage matrix operation

    cdef int add_multiple_of_row_c(self, long x, long y, s, bint col_start) except -1
    cdef int swap_rows_c(self, long x, long y) except -1
    cdef int rescale_row_c(self, long x, s, bint col_start) except -1
    cdef int rescale_column_c(self, long y, s, bint start_row) except -1
    cdef int pivot(self, long x, long y) except -1  # Not a Sage matrix operation
    cdef list gauss_jordan_reduce(self, columns)   # Not a Sage matrix operation

    cdef list nonzero_positions_in_row(self, long r)

    cdef LeanMatrix transpose(self)
    cdef LeanMatrix _matrix_times_matrix_(self, LeanMatrix other)
    cdef LeanMatrix matrix_from_rows_and_columns(self, rows, columns)

    cdef shifting_all(self, P_rows, P_cols, Q_rows, Q_cols, int m)
    cdef shifting(self, U_1, V_2, U_2, V_1, z2, z1, int m)

cdef class GenericMatrix(LeanMatrix):
    cdef _base_ring, _characteristic
    cdef list _entries
    cdef _zero
    cdef _one

    cdef inline row_inner_product(self, long i, long j)   # Not a Sage matrix operation

cdef class BinaryMatrix(LeanMatrix):
    cdef bitset_t* _M
    cdef bitset_t _temp

    cdef inline long row_len(self, long i) except -1   # Not a Sage matrix operation
    cdef inline bint row_inner_product(self, long i, long j) noexcept   # Not a Sage matrix operation

    cdef inline bint get(self, long x, long y) noexcept   # Not a Sage matrix operation
    cdef inline void set(self, long x, long y) noexcept   # Not a Sage matrix operation

    cdef inline list row_sum(self, object L)   # Not a Sage matrix operation
    cdef inline list row_union(self, object L)   # Not a Sage matrix operation

    cdef LeanMatrix matrix_from_rows_and_columns(self, rows, columns)
    cdef matrix_from_rows_and_columns_reordered(self, rows, columns)

    cdef list _character(self, bitset_t x)
    cdef BinaryMatrix _distinguish_by(self, BinaryMatrix P)
    cdef BinaryMatrix _splice_by(self, BinaryMatrix P)
    cdef BinaryMatrix _isolate(self, long r)
    cdef BinaryMatrix equitable_partition(self, BinaryMatrix P=*)   # Not a Sage matrix operation
    cdef bint is_isomorphic(self, BinaryMatrix other, BinaryMatrix s_eq=*, BinaryMatrix o_eq=*) except -2   # Not a Sage matrix operation


cdef class TernaryMatrix(LeanMatrix):
    cdef bitset_t *_M0    # _M0[i] = support of row i
    cdef bitset_t *_M1    # _M1[i] = negative support of row i
    cdef bitset_t _s, _t, _u    # registers

    cdef inline long get(self, long r, long c) noexcept   # Not a Sage matrix operation
    cdef inline int set(self, long r, long c, x) except -1   # Not a Sage matrix operation

    cdef bint _is_negative(self, long r, long c) noexcept

    cdef inline long row_len(self, long i) noexcept   # Not a Sage matrix operation
    cdef inline long row_inner_product(self, long i, long j) noexcept   # Not a Sage matrix operation
    cdef void row_subs(self, long x, long y) noexcept   # Not a Sage matrix operation
    cdef void _row_negate(self, long x) noexcept
    cdef matrix_from_rows_and_columns_reordered(self, rows, columns)

cdef class QuaternaryMatrix(LeanMatrix):
    cdef bitset_t *_M0    # _M0[i] = 1-support of row i
    cdef bitset_t *_M1    # _M1[i] = x- support of row i
    cdef bitset_t _s, _t, _u    # registers
    cdef object _gf4, _zero, _one, _x_zero, _x_one

    cdef inline get(self, long r, long c)   # Not a Sage matrix operation
    cdef inline int set(self, long r, long c, x) except -1   # Not a Sage matrix operation

    cdef inline long row_len(self, long i) except -1   # Not a Sage matrix operation
    cdef inline row_inner_product(self, long i, long j)   # Not a Sage matrix operation
    cdef inline int _row_div(self, long x, object s) except -1
    cdef matrix_from_rows_and_columns_reordered(self, rows, columns)
    cdef void conjugate(self) noexcept   # Not a Sage matrix operation


cdef class PlusMinusOneMatrix(LeanMatrix):
    cdef int* _entries

    cdef inline int get(self, long r, long c) noexcept   # Not a Sage matrix operation
    cdef inline void set(self, long r, long c, int x) noexcept   # Not a Sage matrix operation

    cdef inline long row_len(self, long i) except -1   # Not a Sage matrix operation
    cdef inline row_inner_product(self, long i, long j)   # Not a Sage matrix operation

cdef class RationalMatrix(LeanMatrix):
    cdef mpq_t* _entries

    cdef inline long index(self, long r, long c) noexcept   # Not a Sage matrix operation
    cdef inline void set(self, long r, long c, mpq_t x) noexcept   # Not a Sage matrix operation

    cdef inline long row_len(self, long i) except -1   # Not a Sage matrix operation
    cdef inline row_inner_product(self, long i, long j)   # Not a Sage matrix operation

    cdef int add_multiple_of_row_mpq(self, long x, long y, mpq_t s, bint col_start) except -1
    cdef int rescale_row_mpq(self, long x, mpq_t s, bint col_start) except -1
    cdef int rescale_column_mpq(self, long y, mpq_t s, bint start_row) except -1

cpdef GenericMatrix generic_identity(n, ring)
