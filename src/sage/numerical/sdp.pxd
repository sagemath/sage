from sage.structure.sage_object cimport SageObject
from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.numerical.backends.generic_sdp_backend cimport GenericSDPBackend
cdef class SDPVariable


cdef class SemidefiniteProgram(SageObject):
    cdef GenericSDPBackend _backend
    cdef list _first_variable_names
    cdef list _sdpvariables
    cdef SDPVariable _default_sdpvariable
    cdef dict _variables
    cdef object _linear_functions_parent
    cdef object _linear_constraints_parent
    cpdef int number_of_constraints(self) noexcept
    cpdef int number_of_variables(self) noexcept
    cdef list _constraints
    cpdef sum(self, L) noexcept
    cpdef dual_variable(self, int i, sparse=*) noexcept
    cpdef slack(self, int i, sparse=*) noexcept


cdef class SDPVariable(Element):
    cdef SemidefiniteProgram _p
    cdef dict _dict
    cdef str _name
    cdef bint _hasname
    cdef _matrix_rmul_impl(self, m) noexcept
    cdef _matrix_lmul_impl(self, m) noexcept
    cpdef _acted_upon_(self, mat, bint self_on_left) noexcept


cdef class SDPVariableParent(Parent):
    pass


cdef SDPVariableParent sdp_variable_parent
