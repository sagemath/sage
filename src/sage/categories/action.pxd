from sage.structure.element cimport Element
from .morphism cimport Morphism
from .map cimport Map
from .functor cimport Functor

cdef class Action(Functor):
    cdef G
    cdef US
    cdef bint _is_left
    cdef op
    cdef underlying_set(self)
    cpdef _call_(self, a, b)


cdef class InverseAction(Action):
    cdef Action _action
    cdef Map S_precomposition

cdef class PrecomposedAction(Action):
    cdef Action _action
    cdef Map left_precomposition
    cdef Map right_precomposition

cdef class ActionEndomorphism(Morphism):
    cdef Action _action
    cdef Element _g
