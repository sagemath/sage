from sage.categories.map cimport Map
from sage.structure.parent cimport Parent


cdef class RingExtension_generic(Parent):
    cdef _type
    cdef _backend
    cdef _defining_morphism
    cdef _backend_defining_morphism
    cdef dict _print_options
    cdef bint _import_methods
    cdef bint _is_backend_exposed
    # For division
    cdef RingExtension_generic _fraction_field
    cdef type _fraction_field_type

    cpdef is_defined_over(self, base)
    cpdef Parent _check_base(self, Parent base)
    cpdef _degree_over(self, Parent base)
    cpdef _is_finite_over(self, Parent base)
    cpdef _is_free_over(self, Parent base)
    cdef Map _defining_morphism_fraction_field(self, bint extend_base)


cdef class RingExtensionFractionField(RingExtension_generic):
    cdef _ring


cdef class RingExtensionWithBasis(RingExtension_generic):
    cdef _basis
    cdef _basis_names
    cdef _basis_latex_names

    cpdef _basis_over(self, Parent base)
    # cpdef _free_module(self, Parent base, bint map)


cdef class RingExtensionWithGen(RingExtensionWithBasis):
    cdef _gen
    cdef _name
    cdef public object _latex_names
