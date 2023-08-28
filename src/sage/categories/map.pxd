# sage_setup: distribution = sagemath-objects
from sage.structure.parent cimport Parent
from sage.structure.element cimport Element

cdef class Map(Element):
    cdef object __weakref__

    cdef public int _coerce_cost
    # a rough measure of the cost of using this morphism in the coercion system.
    # 10 by default, 100 if a DefaultCoercionMorphism, 10000 if inexact.

    cdef _update_slots(self, dict)
    cdef dict _extra_slots(self)

    # these methods require x is an element of domain, and returns an element with parent codomain
    cpdef Element _call_(self, x)
    cpdef Element _call_with_args(self, x, args=*, kwds=*)

    cdef public domain    # will be either a weakref or a constant map
    cdef public codomain  # will be a constant map
    cdef Parent _codomain # for accessing the codomain directly
    cdef object _category_for  # category in which this is a morphism

    cdef public _repr_type_str
    cdef public bint _is_coercion

    cpdef _pow_int(self, n)


cdef class Section(Map):
    cdef Map _inverse

cdef class FormalCompositeMap(Map):
    cdef __list
