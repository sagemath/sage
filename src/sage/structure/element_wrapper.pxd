# sage_setup: distribution = sagemath-objects

from sage.structure.element cimport Element

cdef class ElementWrapper(Element):
    cdef public object value

    cpdef _richcmp_(left, right, int op)
    cpdef bint _lt_by_value(self, other) noexcept

cdef class ElementWrapperCheckWrappedClass(ElementWrapper):
    pass

