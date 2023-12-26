from sage.data_structures.bitset cimport bitset_t
from sage.matroids.matroid cimport Matroid
from sage.matroids.basis_exchange_matroid cimport BasisExchangeMatroid
from sage.matroids.set_system cimport SetSystem

cdef class BasisMatroid(BasisExchangeMatroid):
    cdef bitset_t _bb
    cdef bitset_t _b
    cdef SetSystem _nonbases
    cdef _bases_invariant_var
    cdef SetSystem _bases_partition_var
    cdef _bases_invariant2_var
    cdef SetSystem _bases_partition2_var
    cdef _bases_invariant3_var
    cdef SetSystem _bases_partition3_var

    cdef reset_current_basis(self) noexcept

    cpdef _is_basis(self, X) noexcept

    cpdef bases_count(self) noexcept
    cpdef bases(self) noexcept
    cpdef nonbases(self) noexcept

    cpdef truncation(self) noexcept
    cpdef _extension(self, e, H) noexcept
    cpdef _with_coloop(self, e) noexcept
    cpdef relabel(self, l) noexcept

    cpdef _bases_invariant(self) noexcept
    cpdef _bases_partition(self) noexcept
    cpdef _bases_invariant2(self) noexcept
    cpdef _bases_partition2(self) noexcept
    cpdef _bases_invariant3(self) noexcept
    cpdef _bases_partition3(self) noexcept
    cdef _reset_invariants(self) noexcept
    cpdef  bint is_distinguished(self, e) noexcept
    cpdef _is_relaxation(self, M, morphism) noexcept
    cpdef _is_isomorphism(self, M, morphism) noexcept
    cpdef _isomorphism(self, other) noexcept
    cpdef _is_isomorphic(self, other, certificate=*) noexcept


cdef  binom_init(long n, long k) noexcept
cdef  long set_to_index(bitset_t S) noexcept
cdef  index_to_set(bitset_t, long, long, long) noexcept
