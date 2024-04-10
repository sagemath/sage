# sage_setup: distribution = sagemath-objects

from sage.structure.parent cimport Parent
from sage.structure.coerce_dict cimport TripleDict

cpdef py_scalar_parent(py_type)
cpdef py_scalar_to_element(py)
cpdef bint parent_is_integers(P) except -1
cpdef bint is_numpy_type(t) noexcept
cpdef bint is_mpmath_type(t) noexcept


cdef class CoercionModel:
    # This MUST be a mapping to tuples, where each
    # tuple contains at least two elements that are either
    # None or of type Morphism.
    cdef readonly TripleDict _coercion_maps

    # This MUST be a mapping to actions.
    cdef readonly TripleDict _action_maps

    cpdef canonical_coercion(self, x, y)
    cpdef bin_op(self, x, y, op)
    cpdef richcmp(self, x, y, int op)

    cpdef coercion_maps(self, R, S)
    cpdef discover_coercion(self, R, S)
    cpdef verify_coercion_maps(self, R, S, homs, bint fix=*)
    cpdef verify_action(self, action, R, S, op, bint fix=*)

    cpdef get_action(self, R, S, op=*, r=*, s=*)
    cpdef discover_action(self, R, S, op, r=*, s=*)

    cdef bint _record_exceptions
    cpdef _record_exception(self)
    cdef readonly list _exception_stack
    cdef bint _exceptions_cleared

    cdef TripleDict _division_parents
    cpdef analyse(self, xp, yp, op=*)
    cpdef division_parent(self, Parent P)


# Unique global coercion_model instance
cdef CoercionModel coercion_model
