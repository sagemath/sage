from sage.structure.element cimport Element
from sage.categories.morphism cimport Morphism
from sage.categories.map cimport Map
from sage.categories.functor cimport Functor

cdef class Action(Functor):
    cdef _G            # strong reference (default)
    cdef _G_weakref    # weak reference (set by _make_actor_weak)
    cdef bint _pinned  # if True, _make_actor_weak is a no-op
    cdef readonly op
    cdef readonly bint _is_left
    cdef US
    cdef underlying_set(self)
    cdef _actor(self)

    cpdef _pin_actor(self)
    cpdef _make_actor_weak(self)

    cdef _act_convert(self, g, x)
    cpdef _act_(self, g, x)


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
