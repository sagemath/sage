from sage.structure.element cimport AlgebraElement, RingElement

cdef class PowerSeries(AlgebraElement):
    cdef char _is_gen
    cdef _prec
    cdef common_prec_c(self, PowerSeries other)
    #_prec(self, RingElement right_r)

    # UNSAFE, only call from an inplace operator
    # may return a new element if not possible to modify inplace
    cdef _inplace_truncate(self, long prec)
