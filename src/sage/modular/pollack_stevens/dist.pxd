from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport ModuleElement
from sage.categories.action cimport Action
from sage.rings.padics.pow_computer cimport PowComputer_class


cdef class Dist(ModuleElement):
    cpdef normalize(self, include_zeroth_moment=*)
    cdef long ordp
    cpdef long _ord_p(self) noexcept
    cdef long _relprec(self) noexcept
    cdef _unscaled_moment(self, long i)

cdef class Dist_vector(Dist):
    cdef public _moments
    cdef Dist_vector _new_c(self)
    cdef Dist_vector _addsub(self, Dist_vector right, bint negate)
    cpdef _add_(self, other)


cdef class WeightKAction(Action):
    cdef public _k
    cdef public _character
    cdef public _adjuster
    cdef public _p
    cdef public _Np
    cdef public _actmat
    cdef public _maxprecs
    cdef public _symk
    cdef public _dettwist
    cdef public _Sigma0

    cpdef acting_matrix(self, g, M)
    cpdef _compute_acting_matrix(self, g, M)

cdef class WeightKAction_vector(WeightKAction):
    pass
