from sage.structure.element cimport Element
from sage.structure.element_wrapper cimport ElementWrapper
from sage.structure.sage_object cimport SageObject
from sage.modules.with_basis.indexed_element cimport IndexedFreeModuleElement

cdef class LieAlgebraElement(IndexedFreeModuleElement):
    cpdef lift(self)

cdef class LieAlgebraElementWrapper(ElementWrapper):
    cpdef _add_(self, right)
    cpdef _sub_(self, right)

cdef class LieAlgebraMatrixWrapper(LieAlgebraElementWrapper):
    pass

cdef class LieSubalgebraElementWrapper(LieAlgebraElementWrapper):
    cdef dict _monomial_coefficients
    cpdef dict monomial_coefficients(self, bint copy=*)

cdef class StructureCoefficientsElement(LieAlgebraMatrixWrapper):
    cpdef bracket(self, right)
    cpdef _bracket_(self, right)
    cpdef _vector_(self, bint sparse=*, order=*)
    cpdef to_vector(self, bint sparse=*, order=*)
    cpdef dict monomial_coefficients(self, bint copy=*)
    # cpdef lift(self)

cdef class UntwistedAffineLieAlgebraElement(Element):
    cdef dict _t_dict
    cdef _c_coeff
    cdef _d_coeff
    cdef long _hash

    cpdef _add_(self, other)
    cpdef _sub_(self, other)
    cpdef _neg_(self)

    cpdef dict t_dict(self)
    cpdef c_coefficient(self)
    cpdef d_coefficient(self)

    cpdef bracket(self, y)
    cpdef _bracket_(self, y)
    cpdef canonical_derivation(self)
    cpdef monomial_coefficients(self, bint copy=*)

cdef class LieObject(SageObject):
    cdef tuple _word
    cdef public tuple _index_word
    cpdef tuple to_word(self)

cdef class LieGenerator(LieObject):
    cdef public str _name
    cpdef lift(self, dict UEA_gens_dict)

cdef class LieBracket(LieObject):
    cdef public LieObject _left
    cdef public LieObject _right
    cdef long _hash

    cpdef lift(self, dict UEA_gens_dict)

cdef class GradedLieBracket(LieBracket):
    cdef public _grade

cdef class LyndonBracket(GradedLieBracket):
    pass
