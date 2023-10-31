from sage.data_structures.bitset cimport bitset_t

cdef class SetSystem:
    cdef long _groundset_size, _bitset_size
    cdef tuple _groundset
    cdef dict _idx
    cdef bitset_t* _subsets
    cdef long _len, _capacity
    cdef bitset_t _temp

    cdef copy(self) noexcept
    cdef _relabel(self, l) noexcept
    cpdef _complements(self) noexcept

    cdef resize(self, k=*) noexcept
    cdef _append(self, bitset_t X) noexcept
    cdef append(self, X) noexcept
    cdef _subset(self, long k) noexcept
    cdef subset(self, k) noexcept
    cpdef _get_groundset(self) noexcept

    cdef list _incidence_count(self, E) noexcept
    cdef SetSystem _groundset_partition(self, SetSystem P, list cnt) noexcept
    cdef long subset_characteristic(self, SetSystem P, long e) noexcept
    cdef subsets_partition(self, SetSystem P=*, E=*) noexcept
    cdef _distinguish(self, Py_ssize_t v) noexcept
    cpdef is_connected(self) noexcept

    cdef initial_partition(self, SetSystem P=*, E=*) noexcept
    cpdef _equitable_partition(self, SetSystem P=*, EP=*) noexcept
    cpdef _heuristic_partition(self, SetSystem P=*, EP=*) noexcept
    cpdef _isomorphism(self, SetSystem other, SetSystem SP=*, SetSystem OP=*) noexcept
    cpdef _equivalence(self, is_equiv, SetSystem other, SetSystem SP=*, SetSystem OP=*) noexcept

cdef class SetSystemIterator:
    cdef SetSystem _H
    cdef long _pointer, _len
