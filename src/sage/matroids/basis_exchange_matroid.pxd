from sage.data_structures.bitset cimport *
from sage.data_structures.bitset_base cimport bitset_t, bitset_s

from sage.matroids.matroid cimport Matroid
from sage.matroids.set_system cimport SetSystem

cdef class BasisExchangeMatroid(Matroid):
    cdef long _groundset_size, _matroid_rank, _bitset_size
    cdef bitset_t _current_basis, _inside, _outside, _input, _input2, _output, _temp
    cdef tuple _E
    cdef dict _idx
    cdef frozenset _groundset

    cdef _bcount
    cdef _weak_invariant_var, _strong_invariant_var, _heuristic_invariant_var
    cdef SetSystem _weak_partition_var, _strong_partition_var, _heuristic_partition_var

    cdef _relabel(self, l) noexcept

    cdef _pack(self, bitset_t, X) noexcept
    cdef __unpack(self, bitset_t) noexcept
    cdef bint _is_exchange_pair(self, long x, long y) except -1
    cdef int _exchange(self, long x, long y) except -1
    cdef int _move(self, bitset_t X, bitset_t Y) except -1
    cdef __fundamental_cocircuit(self, bitset_t, long x) noexcept
    cdef __fundamental_circuit(self, bitset_t, long y) noexcept

    cdef __max_independent(self, bitset_t, bitset_t) noexcept
    cdef __circuit(self, bitset_t, bitset_t) noexcept
    cdef __closure(self, bitset_t, bitset_t) noexcept
    cdef __max_coindependent(self, bitset_t, bitset_t) noexcept
    cdef __cocircuit(self, bitset_t, bitset_t) noexcept
    cdef _coclosure_internal(self, bitset_t, bitset_t) noexcept

    cdef __augment(self, bitset_t, bitset_t, bitset_t) noexcept
    cdef bint __is_independent(self, bitset_t F) except -1
    cdef __move_current_basis(self, bitset_t, bitset_t) noexcept

    cdef bint _set_current_basis(self, F) noexcept

    cpdef groundset(self) noexcept
    cpdef groundset_list(self) noexcept
    cpdef full_rank(self) noexcept
    cpdef full_corank(self) noexcept

    cpdef basis(self) noexcept
    cpdef _move_current_basis(self, X, Y) noexcept

    cpdef _max_independent(self, F) noexcept
    cpdef _rank(self, F) noexcept
    cpdef _circuit(self, F) noexcept
    cpdef _fundamental_circuit(self, B, e) noexcept
    cpdef _closure(self, F) noexcept

    cpdef _max_coindependent(self, F) noexcept
    cpdef _corank(self, F) noexcept
    cpdef _cocircuit(self, F) noexcept
    cpdef _fundamental_cocircuit(self, B, e) noexcept
    cpdef _coclosure(self, F) noexcept

    cpdef _augment(self, X, Y) noexcept
    cpdef _is_independent(self, F) noexcept

    cpdef f_vector(self) noexcept
    cdef  _f_vector_rec(self, object f_vec, bitset_t* flats, bitset_t* todo, long elt, long rnk) noexcept
    cpdef flats(self, R) noexcept
    cdef  _flats_rec(self, SetSystem Rflats, long R, bitset_t* flats, bitset_t* todo, long elt, long rnk) noexcept
    cpdef coflats(self, R) noexcept
    cdef  _coflats_rec(self, SetSystem Rcoflats, long R, bitset_t* coflats, bitset_t* todo, long elt, long cornk) noexcept
    cdef _flat_element_inv(self, long k) noexcept
    cdef  _flat_element_inv_rec(self, object f_inc, long R, bitset_t* flats, bitset_t* todo, long elt, long i) noexcept

    cpdef bases_count(self) noexcept
    cpdef independent_r_sets(self, long r) noexcept
    cpdef bases(self) noexcept
    cpdef dependent_r_sets(self, long r) noexcept
    cpdef nonbases(self) noexcept

    cpdef nonspanning_circuits(self) noexcept
    cpdef cocircuits(self) noexcept
    cpdef circuits(self) noexcept

    cpdef _characteristic_setsystem(self) noexcept
    cpdef _weak_invariant(self) noexcept
    cpdef _weak_partition(self) noexcept
    cpdef _strong_invariant(self) noexcept
    cpdef _strong_partition(self) noexcept
    cpdef _heuristic_invariant(self) noexcept
    cpdef _heuristic_partition(self) noexcept
    cdef _flush(self) noexcept

    cpdef _equitable_partition(self, P=*) noexcept
    cpdef _is_isomorphic(self, other, certificate=*) noexcept
    cpdef _isomorphism(self, other) noexcept
    cpdef _is_isomorphism(self, other, morphism) noexcept
    cdef bint __is_isomorphism(self, BasisExchangeMatroid other, morphism) noexcept

    cpdef is_valid(self) noexcept

cdef bint nxksrd(bitset_s *b, long n, long k, bint succ) noexcept
