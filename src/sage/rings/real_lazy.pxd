# sage_setup: distribution = sagemath-categories
from sage.rings.ring cimport Field
from sage.structure.element cimport RingElement, ModuleElement, Element, FieldElement

cdef class LazyField(Field):
    cpdef interval_field(self, prec=*)

cdef class LazyFieldElement(FieldElement):
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cdef LazyFieldElement _new_wrapper(self, value)
    cdef LazyFieldElement _new_binop(self, LazyFieldElement left, LazyFieldElement right, op)
    cdef LazyFieldElement _new_unop(self, LazyFieldElement arg, op)
    cpdef eval(self, R)
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

