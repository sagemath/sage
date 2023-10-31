from sage.structure.parent cimport Parent, Parent_richcmp_element_without_coercion
from sage.structure.element cimport ModuleElement, RingElement, Element

cpdef is_LinearFunction(x) noexcept

cdef class LinearFunctionOrConstraint(ModuleElement):
    pass

cdef class LinearFunctionsParent_class(Parent):
    cpdef _element_constructor_(self, x) noexcept
    cpdef _coerce_map_from_(self, R) noexcept
    cdef public _multiplication_symbol

cdef class LinearFunction(LinearFunctionOrConstraint):
    cdef dict _f
    cpdef _add_(self, other) noexcept
    cpdef iteritems(self) noexcept
    cpdef _acted_upon_(self, x, bint self_on_left) noexcept
    cpdef is_zero(self) noexcept
    cpdef equals(LinearFunction left, LinearFunction right) noexcept

cdef class LinearConstraintsParent_class(Parent):
    cdef LinearFunctionsParent_class _LF
    cpdef _element_constructor_(self, left, right=?, equality=?) noexcept
    cpdef _coerce_map_from_(self, R) noexcept

cdef class LinearConstraint(LinearFunctionOrConstraint):
    cdef bint equality
    cdef list constraints
    cpdef equals(LinearConstraint left, LinearConstraint right) noexcept
