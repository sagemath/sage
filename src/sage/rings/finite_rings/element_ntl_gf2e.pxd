from sage.libs.ntl.types cimport GF2E_c, GF2EContext_c
from sage.rings.finite_rings.finite_field_base cimport FiniteField
from sage.rings.finite_rings.element_base cimport FinitePolyExtElement
from sage.structure.sage_object cimport SageObject

cdef class FiniteField_ntl_gf2eElement(FinitePolyExtElement)

cdef class Cache_ntl_gf2e(SageObject):
    cdef GF2EContext_c F
    cdef FiniteField _parent
    cdef public FiniteField_ntl_gf2eElement _zero_element
    cdef public FiniteField_ntl_gf2eElement _one_element
    cdef public FiniteField_ntl_gf2eElement _gen
    cdef FiniteField_ntl_gf2eElement _new(self)
    cpdef FiniteField_ntl_gf2eElement fetch_int(self, number)

cdef class FiniteField_ntl_gf2eElement(FinitePolyExtElement):
    cdef GF2E_c x
    cdef Cache_ntl_gf2e _cache
    cdef FiniteField_ntl_gf2eElement _new(FiniteField_ntl_gf2eElement self)
