from sage.misc.bitset cimport *

DEF BINT_EXCEPT = -2**31 - 1  # def. repeated in .pyx file

from matroid cimport Matroid
from set_system cimport SetSystem

cdef class BasisExchangeMatroid(Matroid):
    cdef long _groundset_size, _matroid_rank, _bitset_size
    cdef bitset_t _current_basis, _inside, _outside, _input, _input2, _output, _temp
    cdef list _E
    cdef dict _idx
    cdef frozenset _groundset

    cdef _bcount
    cdef _weak_invariant_var, _strong_invariant_var, _heuristic_invariant_var
    cdef SetSystem _weak_partition_var, _strong_partition_var, _heuristic_partition_var

    cdef __relabel(self, l)

    cdef  __pack(self, bitset_t, X)
    cdef  __unpack(self, bitset_t)
    cdef  bint __is_exchange_pair(self, long x, long y) except BINT_EXCEPT
    cdef  bint __exchange(self, long x, long y) except BINT_EXCEPT
    cdef  __move(self, bitset_t X, bitset_t Y)
    cdef  __fundamental_cocircuit(self, bitset_t, long x)
    cdef  __fundamental_circuit(self, bitset_t, long y)

    cdef  __max_independent(self, bitset_t, bitset_t)
    cdef  __circuit(self, bitset_t, bitset_t)
    cdef  __closure(self, bitset_t, bitset_t)
    cdef  __max_coindependent(self, bitset_t, bitset_t)
    cdef  __cocircuit(self, bitset_t, bitset_t)
    cdef  __coclosure(self, bitset_t, bitset_t)

    cdef  __augment(self, bitset_t, bitset_t, bitset_t)
    cdef  bint __is_independent(self, bitset_t F)
    cdef  __move_current_basis(self, bitset_t, bitset_t)

    cdef bint _set_current_basis(self, F)

    cpdef groundset(self)
    cpdef groundset_list(self)
    cpdef full_rank(self)
    cpdef full_corank(self)

    cpdef basis(self)
    cpdef _move_current_basis(self, X, Y)

    cpdef _max_independent(self, F)
    cpdef _rank(self, F)
    cpdef _circuit(self, F)
    cpdef _fundamental_circuit(self, B, e)
    cpdef _closure(self, F)

    cpdef _max_coindependent(self, F)
    cpdef _corank(self, F)
    cpdef _cocircuit(self, F)
    cpdef _fundamental_cocircuit(self, B, e)
    cpdef _coclosure(self, F)

    cpdef _augment(self, X, Y)
    cpdef _is_independent(self, F)

    cpdef f_vector(self)
    cdef  _f_vector_rec(self, object f_vec, bitset_t* flats, bitset_t* todo, long elt, long rnk)
    cpdef flats(self, R)
    cdef  _flats_rec(self, SetSystem Rflats, long R, bitset_t* flats, bitset_t* todo, long elt, long rnk)
    cpdef coflats(self, R)
    cdef  _coflats_rec(self, SetSystem Rcoflats, long R, bitset_t* coflats, bitset_t* todo, long elt, long cornk)
    cdef _flat_element_inv(self, long k)
    cdef  _flat_element_inv_rec(self, object f_inc, long R, bitset_t* flats, bitset_t* todo, long elt, long i)

    cpdef bases_count(self)
    cpdef independent_r_sets(self, long r)
    cpdef bases(self)
    cpdef dependent_r_sets(self, long r)
    cpdef nonbases(self)

    cpdef nonspanning_circuits(self)
    cpdef cocircuits(self)
    cpdef circuits(self)

    cpdef _characteristic_setsystem(self)
    cpdef _weak_invariant(self)
    cpdef _weak_partition(self)
    cpdef _strong_invariant(self)
    cpdef _strong_partition(self)
    cpdef _heuristic_invariant(self)
    cpdef _heuristic_partition(self)
    cdef _flush(self)

    cpdef _equitable_partition(self, P=*)
    cpdef _is_isomorphic(self, other)
    cpdef _is_isomorphism(self, other, morphism)
    cdef bint __is_isomorphism(self, BasisExchangeMatroid other, morphism)

    cpdef is_valid(self)

cdef bint nxksrd(bitset_s *b, long n, long k, bint succ)
