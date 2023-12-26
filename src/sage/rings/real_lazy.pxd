from sage.rings.ring cimport Field
from sage.structure.element cimport RingElement, ModuleElement, Element, FieldElement

cdef class LazyField(Field):
    cpdef interval_field(self, prec=*) noexcept

cdef class LazyFieldElement(FieldElement):
    cpdef _add_(self, other) noexcept
    cpdef _mul_(self, other) noexcept
    cdef LazyFieldElement _new_wrapper(self, value) noexcept
    cdef LazyFieldElement _new_binop(self, LazyFieldElement left, LazyFieldElement right, op) noexcept
    cdef LazyFieldElement _new_unop(self, LazyFieldElement arg, op) noexcept
    cpdef eval(self, R) noexcept
    cpdef int depth(self) noexcept

cdef class LazyWrapper(LazyFieldElement):
    cdef readonly _value

cdef class LazyBinop(LazyFieldElement):
    cdef readonly LazyFieldElement _left
    cdef readonly LazyFieldElement _right
    cdef readonly _op

cdef class LazyUnop(LazyFieldElement):
    cdef readonly LazyFieldElement _arg
    cdef readonly _op

cdef class LazyNamedUnop(LazyUnop):
    cdef readonly _extra_args

