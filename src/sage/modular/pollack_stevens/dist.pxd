from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport ModuleElement
from sage.categories.action cimport Action
from sage.rings.padics.pow_computer cimport PowComputer_class
#from sage.libs.flint.ulong_extras cimport *

#cdef extern from "../../../ext/multi_modular.h":
#    ctypedef unsigned long mod_int
#    mod_int MOD_INT_MAX



cdef class Dist(ModuleElement):
    cpdef normalize(self)
    cdef long ordp
    cpdef long _ord_p(self)
    cdef long _relprec(self)
    cdef _unscaled_moment(self, long i)

cdef class Dist_vector(Dist):
    cdef public _moments
    cdef Dist_vector _new_c(self)
    cdef Dist_vector _addsub(self, Dist_vector right, bint negate)

# cdef class Dist_simple(Dist):
#     cdef public _moments
#     cdef Dist_simple _new_c(self)
#     cdef Dist_simple _addsub(self, Dist_simple right, bint negate)

#cdef class Dist2(Dist): # only works on 64-bit....
#    cdef long[60] moments
#    cdef int prec
#    cdef public PowComputer_long prime_pow
#    cdef Dist2 _new_c(self)

# cdef class Dist_long(Dist):
#     cdef long[60] _moments # 38 once 2 is special-cased
#     cdef int relprec
#     cdef public PowComputer_class prime_pow
#     cdef int quasi_normalize(self) except -1
#     cdef Dist_long _new_c(self)
#     cdef Dist_long _addsub(self, Dist_long right, bint negate)

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

cdef class WeightKAction_simple(WeightKAction):
    pass

cdef class SimpleMat(SageObject):
    cdef long* _mat
    cdef long M
    cdef bint _inited

cdef class WeightKAction_long(WeightKAction):
    pass

cdef class iScale(Action):
    pass
