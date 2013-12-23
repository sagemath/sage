from sage.misc.bitset cimport bitset_t

cdef class SetSystem:
    cdef long _groundset_size, _bitset_size
    cdef list _groundset
    cdef dict _idx
    cdef bitset_t* _subsets
    cdef long _len, _capacity
    cdef bitset_t _temp

    cdef copy(self)
    cdef _relabel(self, l)
    cpdef _complements(self)

    cdef inline resize(self, k=*)
    cdef inline _append(self, bitset_t X)
    cdef inline append(self, X)
    cdef inline _subset(self, long k)
    cdef subset(self, k)
    cpdef _get_groundset(self)

    cdef list _incidence_count(self, E)
    cdef SetSystem _groundset_partition(self, SetSystem P, list cnt)
    cdef long subset_characteristic(self, SetSystem P, long e)
    cdef subsets_partition(self, SetSystem P=*, E=*)
    cdef _distinguish(self, v)

    cdef initial_partition(self, SetSystem P=*, E=*)
    cpdef _equitable_partition(self, SetSystem P=*, EP=*)
    cpdef _heuristic_partition(self, SetSystem P=*, EP=*)
    cpdef _isomorphism(self, SetSystem other, SetSystem SP=*, SetSystem OP=*)
    cpdef _equivalence(self, is_equiv, SetSystem other, SetSystem SP=*, SetSystem OP=*)

cdef class SetSystemIterator:
    cdef SetSystem _H
    cdef long _pointer, _len
