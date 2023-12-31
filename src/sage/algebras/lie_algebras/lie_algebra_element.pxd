from sage.structure.element cimport Element
from sage.structure.element_wrapper cimport ElementWrapper
from sage.structure.sage_object cimport SageObject
from sage.modules.with_basis.indexed_element cimport IndexedFreeModuleElement

cdef class LieAlgebraElement(IndexedFreeModuleElement):
    cpdef lift(self) noexcept

cdef class LieAlgebraElementWrapper(ElementWrapper):
    cpdef _add_(self, right) noexcept
    cpdef _sub_(self, right) noexcept

cdef class LieAlgebraMatrixWrapper(LieAlgebraElementWrapper):
    pass

cdef class LieSubalgebraElementWrapper(LieAlgebraElementWrapper):
    cdef dict _monomial_coefficients
    cpdef dict monomial_coefficients(self, bint copy=*) noexcept

cdef class StructureCoefficientsElement(LieAlgebraMatrixWrapper):
    cpdef bracket(self, right) noexcept
    cpdef _bracket_(self, right) noexcept
    cpdef to_vector(self, bint sparse=*) noexcept
    cpdef dict monomial_coefficients(self, bint copy=*) noexcept
    # cpdef lift(self)

cdef class UntwistedAffineLieAlgebraElement(Element):
    cdef dict _t_dict
    cdef _c_coeff
    cdef _d_coeff
    cdef long _hash

    cpdef _add_(self, other) noexcept
    cpdef _sub_(self, other) noexcept
    cpdef _neg_(self) noexcept

    cpdef dict t_dict(self) noexcept
    cpdef c_coefficient(self) noexcept
    cpdef d_coefficient(self) noexcept

    cpdef bracket(self, y) noexcept
    cpdef _bracket_(self, y) noexcept
    cpdef canonical_derivation(self) noexcept
    cpdef monomial_coefficients(self, bint copy=*) noexcept

cdef class LieObject(SageObject):
    cdef tuple _word
    cdef public tuple _index_word
    cpdef tuple to_word(self) noexcept

cdef class LieGenerator(LieObject):
    cdef public str _name
    cpdef lift(self, dict UEA_gens_dict)

cdef class LieBracket(LieObject):
    cdef public LieObject _left
    cdef public LieObject _right
    cdef long _hash

    cpdef lift(self, dict UEA_gens_dict) noexcept

cdef class GradedLieBracket(LieBracket):
    cdef public _grade

cdef class LyndonBracket(GradedLieBracket):
    pass
