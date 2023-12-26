from sage.structure.element cimport Element
from sage.categories.morphism cimport Morphism
from sage.categories.map cimport Map
from sage.categories.functor cimport Functor

cdef class Action(Functor):
    cdef readonly G
    cdef readonly op
    cdef readonly bint _is_left
    cdef US
    cdef underlying_set(self) noexcept

    cdef _act_convert(self, g, x) noexcept
    cpdef _act_(self, g, x) noexcept


cdef class InverseAction(Action):
    cdef Action _action
    cdef Map S_precomposition

cdef class PrecomposedAction(Action):
    cdef Action _action
    cdef Map G_precomposition
    cdef Map S_precomposition

cdef class ActionEndomorphism(Morphism):
    cdef Action _action
    cdef Element _g
