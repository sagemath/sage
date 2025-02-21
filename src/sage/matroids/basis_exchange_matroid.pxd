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

    cdef _relabel(self, mapping)

    cdef _pack(self, bitset_t, X)
    cdef __unpack(self, bitset_t)
    cdef bint _is_exchange_pair(self, long x, long y) except -1
    cdef int _exchange(self, long x, long y) except -1
    cdef int _move(self, bitset_t X, bitset_t Y) except -1
    cdef __fundamental_cocircuit(self, bitset_t, long x)
    cdef __fundamental_circuit(self, bitset_t, long y)

    cdef __max_independent(self, bitset_t, bitset_t)
    cdef __circuit(self, bitset_t, bitset_t)
    cdef __closure(self, bitset_t, bitset_t)
    cdef __max_coindependent(self, bitset_t, bitset_t)
    cdef __cocircuit(self, bitset_t, bitset_t)
    cdef _coclosure_internal(self, bitset_t, bitset_t)

    cdef __augment(self, bitset_t, bitset_t, bitset_t)
    cdef bint __is_independent(self, bitset_t F) except -1
    cdef __move_current_basis(self, bitset_t, bitset_t)

    cdef bint _set_current_basis(self, F) noexcept

    cpdef frozenset groundset(self)
    cpdef list groundset_list(self)
    cpdef full_rank(self)
    cpdef full_corank(self)

    cpdef basis(self)
    cpdef _move_current_basis(self, X, Y)

    cpdef frozenset _max_independent(self, frozenset F)
    cpdef int _rank(self, frozenset F) except? -1
    cpdef frozenset _circuit(self, frozenset F)
    cpdef frozenset _fundamental_circuit(self, frozenset B, e)
    cpdef frozenset _closure(self, frozenset F)

    cpdef frozenset _max_coindependent(self, frozenset F)
    cpdef int _corank(self, frozenset F) noexcept
    cpdef frozenset _cocircuit(self, frozenset F)
    cpdef frozenset _fundamental_cocircuit(self, frozenset B, e)
    cpdef frozenset _coclosure(self, frozenset F)

    cpdef frozenset _augment(self, frozenset X, frozenset Y)
    cpdef bint _is_independent(self, frozenset F) noexcept

    cpdef list whitney_numbers2(self)
    cdef  _whitney_numbers2_rec(self, object f_vec, bitset_t* flats, bitset_t* todo, long elt, long rnk)
    cdef  _flats_rec(self, SetSystem Rflats, long R, bitset_t* flats, bitset_t* todo, long elt, long rnk)
    cdef  _coflats_rec(self, SetSystem Rcoflats, long R, bitset_t* coflats, bitset_t* todo, long elt, long cornk)
    cdef _flat_element_inv(self, long k)
    cdef  _flat_element_inv_rec(self, object f_inc, long R, bitset_t* flats, bitset_t* todo, long elt, long i)

    cpdef bases_count(self)
    cpdef SetSystem independent_sets(self, long k=*)
    cpdef SetSystem dependent_sets(self, long k)

    cpdef SetSystem nonspanning_circuits(self)
    cpdef SetSystem cocircuits(self)
    cpdef SetSystem circuits(self, k=*)

    cpdef _characteristic_setsystem(self)
    cpdef _weak_invariant(self)
    cpdef _weak_partition(self)
    cpdef _strong_invariant(self)
    cpdef _strong_partition(self)
    cpdef _heuristic_invariant(self)
    cpdef _heuristic_partition(self)
    cdef _flush(self)

    cpdef _equitable_partition(self, P=*)
    cpdef _is_isomorphic(self, other, certificate=*)
    cpdef _isomorphism(self, other)
    cpdef _is_isomorphism(self, other, morphism)
    cdef bint __is_isomorphism(self, BasisExchangeMatroid other, morphism) noexcept

    cpdef is_valid(self, certificate=*)

cdef bint nxksrd(bitset_s *b, long n, long k, bint succ) noexcept
