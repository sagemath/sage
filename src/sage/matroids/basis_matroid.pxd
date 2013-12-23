from sage.misc.bitset cimport bitset_t
from matroid cimport Matroid
from basis_exchange_matroid cimport BasisExchangeMatroid
from set_system cimport SetSystem

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

    cdef  bint __is_exchange_pair(self, long x, long y)
    cdef reset_current_basis(self)

    cpdef _is_basis(self, X)

    cpdef bases_count(self)
    cpdef bases(self)
    cpdef nonbases(self)

    cpdef truncation(self)
    cpdef _extension(self, e, H)
    cpdef _with_coloop(self, e)
    cpdef relabel(self, l)

    cpdef _bases_invariant(self)
    cpdef _bases_partition(self)
    cpdef _bases_invariant2(self)
    cpdef _bases_partition2(self)
    cpdef _bases_invariant3(self)
    cpdef _bases_partition3(self)
    cdef _reset_invariants(self)
    cpdef  bint is_distinguished(self, e)
    cpdef _is_relaxation(self, M, morphism)
    cpdef _is_isomorphism(self, M, morphism)
    cpdef _is_isomorphic(self, other)


cdef  binom_init(long n, long k)
cdef  long set_to_index(bitset_t S)
cdef  index_to_set(bitset_t, long, long, long)
