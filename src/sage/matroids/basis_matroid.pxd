from sage.data_structures.bitset cimport bitset_t
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

    cdef reset_current_basis(self)

    cpdef bint _is_basis(self, frozenset X)

    cpdef bases_count(self)
    cpdef SetSystem bases(self)
    cpdef SetSystem nonbases(self)

    cpdef truncation(self)
    cpdef _extension(self, e, H)
    cpdef _with_coloop(self, e)
    # cpdef relabel(self, mapping)

    cpdef _bases_invariant(self)
    cpdef _bases_partition(self)
    cpdef _bases_invariant2(self)
    cpdef _bases_partition2(self)
    cpdef _bases_invariant3(self)
    cpdef _bases_partition3(self)
    cdef _reset_invariants(self)
    cpdef  bint is_distinguished(self, e) noexcept
    cpdef _is_relaxation(self, M, morphism)
    cpdef _is_isomorphism(self, M, morphism)
    cpdef _isomorphism(self, other)
    cpdef _is_isomorphic(self, other, certificate=*)

cdef  binom_init(long n, long k)
cdef  long set_to_index(bitset_t S) noexcept
cdef  index_to_set(bitset_t, long, long, long)
