
from sage.structure.element cimport Element

cdef class ElementWrapper(Element):
    cdef public object value

    cpdef _richcmp_(left, right, int op) noexcept
    cpdef bint _lt_by_value(self, other) noexcept

cdef class ElementWrapperCheckWrappedClass(ElementWrapper):
    pass

