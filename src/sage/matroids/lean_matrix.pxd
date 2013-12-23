from sage.misc.bitset cimport bitset_t

cdef class LeanMatrix:
    cdef long _nrows
    cdef long _ncols

    cdef LeanMatrix copy(self)   # Deprecated Sage matrix operation
    cdef void resize(self, long k)   # Not a Sage matrix operation
    cdef LeanMatrix stack(self, LeanMatrix M)
    cdef LeanMatrix augment(self, LeanMatrix M)
    cdef LeanMatrix prepend_identity(self)   # Not a Sage matrix operation

    cpdef long ncols(self)
    cpdef long nrows(self)
    cpdef base_ring(self)
    cpdef characteristic(self)   # Not a Sage matrix operation

    cdef get_unsafe(self, long r, long c)
    cdef void set_unsafe(self, long r, long c, x)
    cdef bint is_nonzero(self, long r, long c)   # Not a Sage matrix operation

    cdef void add_multiple_of_row_c(self, long x, long y, s, bint col_start)
    cdef void swap_rows_c(self, long x, long y)
    cdef void rescale_row_c(self, long x, s, bint col_start)
    cdef void rescale_column_c(self, long y, s, bint start_row)
    cdef void pivot(self, long x, long y)   # Not a Sage matrix operation
    cdef list gauss_jordan_reduce(self, columns)   # Not a Sage matrix operation

    cdef list nonzero_positions_in_row(self, long r)

    cdef LeanMatrix transpose(self)
    cdef LeanMatrix _matrix_times_matrix_(self, LeanMatrix other)
    cdef LeanMatrix matrix_from_rows_and_columns(self, rows, columns)

cdef class GenericMatrix(LeanMatrix):
    cdef _base_ring, _characteristic
    cdef list _entries
    cdef _zero
    cdef _one

    cdef LeanMatrix copy(self)   # Deprecated Sage matrix operation
    cdef void resize(self, long k)   # Not a Sage matrix operation
    cdef LeanMatrix stack(self, LeanMatrix M)
    cdef LeanMatrix augment(self, LeanMatrix M)
    cdef LeanMatrix prepend_identity(self)   # Not a Sage matrix operation

    cpdef long ncols(self)
    cpdef long nrows(self)
    cpdef base_ring(self)
    cpdef characteristic(self)   # Not a Sage matrix operation

    cdef get_unsafe(self, long r, long c)
    cdef void set_unsafe(self, long r, long c, x)

    cdef void swap_rows_c(self, long x, long y)

    cdef LeanMatrix transpose(self)
    cdef inline row_inner_product(self, long i, long j)   # Not a Sage matrix operation
    cdef LeanMatrix _matrix_times_matrix_(self, LeanMatrix other)

cdef class BinaryMatrix(LeanMatrix):
    cdef bitset_t* _M
    cdef bitset_t _temp

    cdef LeanMatrix copy(self)   # Deprecated Sage matrix operation
    cdef void resize(self, long k)   # Not a Sage matrix operation
    cdef LeanMatrix stack(self, LeanMatrix M)
    cdef LeanMatrix augment(self, LeanMatrix M)
    cdef LeanMatrix prepend_identity(self)   # Not a Sage matrix operation

    cpdef long ncols(self)
    cpdef long nrows(self)
    cpdef base_ring(self)
    cpdef characteristic(self)   # Not a Sage matrix operation

    cdef get_unsafe(self, long r, long c)
    cdef void set_unsafe(self, long r, long c, x)
    cdef bint is_nonzero(self, long r, long c)   # Not a Sage matrix operation

    cdef void add_multiple_of_row_c(self, long i, long j, s, bint col_start)
    cdef void swap_rows_c(self, long i, long j)
    cdef void pivot(self, long x, long y)   # Not a Sage matrix operation

    cdef list nonzero_positions_in_row(self, long i)

    cdef LeanMatrix transpose(self)
    cdef LeanMatrix _matrix_times_matrix_(self, LeanMatrix other)

    cdef inline long row_len(self, long i)   # Not a Sage matrix operation
    cdef inline bint row_inner_product(self, long i, long j)   # Not a Sage matrix operation

    cdef inline bint get(self, long x, long y)   # Not a Sage matrix operation
    cdef inline void set(self, long x, long y)   # Not a Sage matrix operation

    cdef inline list row_sum(self, object L)   # Not a Sage matrix operation
    cdef inline list row_union(self, object L)   # Not a Sage matrix operation

    cdef LeanMatrix matrix_from_rows_and_columns(self, rows, columns)

    cdef list _character(self, bitset_t x)
    cdef BinaryMatrix _distinguish_by(self, BinaryMatrix P)
    cdef BinaryMatrix _splice_by(self, BinaryMatrix P)
    cdef BinaryMatrix _isolate(self, long r)
    cdef BinaryMatrix equitable_partition(self, BinaryMatrix P=*)   # Not a Sage matrix operation
    cdef bint is_isomorphic(self, BinaryMatrix other, BinaryMatrix s_eq=*, BinaryMatrix o_eq=*)   # Not a Sage matrix operation


cdef class TernaryMatrix(LeanMatrix):
    cdef bitset_t *_M0, *_M1    # _M0[i] = support of row i, _M1[i] = negative support of row i
    cdef bitset_t _s, _t, _u    # registers

    cdef get_unsafe(self, long r, long c)
    cdef void set_unsafe(self, long r, long c, x)

    cdef LeanMatrix copy(self)   # Deprecated Sage matrix operation
    cdef void resize(self, long k)   # Not a Sage matrix operation
    cdef LeanMatrix stack(self, LeanMatrix M)
    cdef LeanMatrix augment(self, LeanMatrix M)
    cdef LeanMatrix prepend_identity(self)   # Not a Sage matrix operation

    cpdef long ncols(self)
    cpdef long nrows(self)
    cpdef base_ring(self)
    cpdef characteristic(self)   # Not a Sage matrix operation

    cdef inline long get(self, long r, long c)   # Not a Sage matrix operation
    cdef inline void set(self, long r, long c, x)   # Not a Sage matrix operation

    cdef bint is_nonzero(self, long r, long c)   # Not a Sage matrix operation
    cdef bint _is_negative(self, long r, long c)

    cdef inline long row_len(self, long i)   # Not a Sage matrix operation
    cdef inline long row_inner_product(self, long i, long j)   # Not a Sage matrix operation
    cdef void add_multiple_of_row_c(self, long x, long y, s, bint col_start)
    cdef void row_subs(self, long x, long y)   # Not a Sage matrix operation
    cdef void _row_negate(self, long x)
    cdef void swap_rows_c(self, long x, long y)
    cdef void pivot(self, long x, long y)   # Not a Sage matrix operation

    cdef list nonzero_positions_in_row(self, long r)
    cdef LeanMatrix transpose(self)
    cdef LeanMatrix _matrix_times_matrix_(self, LeanMatrix other)

cdef class QuaternaryMatrix(LeanMatrix):
    cdef bitset_t *_M0, *_M1    # _M0[i] = 1-support of row i, _M1[i] = x- support of row i
    cdef bitset_t _s, _t, _u    # registers
    cdef object _gf4, _zero, _one, _x_zero, _x_one

    cdef inline get(self, long r, long c)   # Not a Sage matrix operation
    cdef inline void set(self, long r, long c, x)   # Not a Sage matrix operation

    cdef get_unsafe(self, long r, long c)
    cdef void set_unsafe(self, long r, long c, x)
    cdef bint is_nonzero(self, long r, long c)   # Not a Sage matrix operation

    cdef LeanMatrix copy(self)   # Deprecated Sage matrix operation
    cdef void resize(self, long k)   # Not a Sage matrix operation
    cdef LeanMatrix stack(self, LeanMatrix M)
    cdef LeanMatrix augment(self, LeanMatrix M)
    cdef LeanMatrix prepend_identity(self)   # Not a Sage matrix operation

    cpdef base_ring(self)
    cpdef characteristic(self)   # Not a Sage matrix operation

    cdef inline long row_len(self, long i)   # Not a Sage matrix operation
    cdef inline row_inner_product(self, long i, long j)   # Not a Sage matrix operation
    cdef void add_multiple_of_row_c(self, long x, long y, s, bint col_start)
    cdef void swap_rows_c(self, long x, long y)
    cdef inline void _row_div(self, long x, object s)
    cdef void pivot(self, long x, long y)   # Not a Sage matrix operation

    cdef list nonzero_positions_in_row(self, long r)
    cdef LeanMatrix transpose(self)
    cdef void conjugate(self)   # Not a Sage matrix operation
    cdef LeanMatrix _matrix_times_matrix_(self, LeanMatrix other)

cdef class IntegerMatrix(LeanMatrix):
    cdef int* _entries

    cdef inline get(self, long r, long c)   # Not a Sage matrix operation
    cdef inline void set(self, long r, long c, int x)   # Not a Sage matrix operation

    cdef get_unsafe(self, long r, long c)
    cdef void set_unsafe(self, long r, long c, x)
    cdef bint is_nonzero(self, long r, long c)   # Not a Sage matrix operation

    cdef LeanMatrix copy(self)   # Deprecated Sage matrix operation
    cdef void resize(self, long k)   # Not a Sage matrix operation
    cdef LeanMatrix stack(self, LeanMatrix M)
    cdef LeanMatrix augment(self, LeanMatrix M)
    cdef LeanMatrix prepend_identity(self)   # Not a Sage matrix operation

    cpdef base_ring(self)
    cpdef characteristic(self)   # Not a Sage matrix operation

    cdef inline long row_len(self, long i)   # Not a Sage matrix operation
    cdef inline row_inner_product(self, long i, long j)   # Not a Sage matrix operation
    cdef void add_multiple_of_row_c(self, long x, long y, s, bint col_start)
    cdef void swap_rows_c(self, long x, long y)
    cdef void rescale_row_c(self, long x, s, bint col_start)
    cdef void rescale_column_c(self, long y, s, bint start_row)
    cdef void pivot(self, long x, long y)   # Not a Sage matrix operation

    cdef list nonzero_positions_in_row(self, long r)
    cdef LeanMatrix transpose(self)
    cdef LeanMatrix _matrix_times_matrix_(self, LeanMatrix other)

    cdef list gauss_jordan_reduce(self, columns)   # Not a Sage matrix operation

cpdef GenericMatrix generic_identity(n, ring)
