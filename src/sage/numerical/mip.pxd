cdef extern from *:
    ctypedef double* const_double_ptr "const double*"
    cdef int BINARY = 1
    cdef int REAL = -1
    cdef int INTEGER = 0

from sage.structure.sage_object cimport SageObject
from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.numerical.backends.generic_backend cimport GenericBackend
cdef class MIPVariable


cdef class MixedIntegerLinearProgram(SageObject):
    cdef GenericBackend _backend
    cdef list _first_variable_names
    cdef list _mipvariables
    cdef MIPVariable _default_mipvariable
    cdef dict _variables
    cdef int __BINARY
    cdef int __REAL
    cdef int __INTEGER
    cdef object _linear_functions_parent
    cdef object _linear_constraints_parent
    cpdef int number_of_constraints(self)
    cpdef int number_of_variables(self)
    cdef int _check_redundant
    cdef list _constraints
    cpdef sum(self, L)


cdef class MIPVariable(Element):
    cdef MixedIntegerLinearProgram _p
    cdef int _dim
    cdef dict _dict
    cdef int _vtype
    cdef str _name
    cdef bint _hasname
    cdef object _lower_bound
    cdef object _upper_bound
    cdef _matrix_rmul_impl(self, m)
    cdef _matrix_lmul_impl(self, m)
    cpdef _acted_upon_(self, mat, bint self_on_left)


cdef class MIPVariableParent(Parent):
    pass

    
cdef MIPVariableParent mip_variable_parent

