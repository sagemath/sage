# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport GF2E_c, GF2EContext_c
from sage.rings.finite_rings.finite_field_base cimport FiniteField
from sage.rings.finite_rings.element_base cimport FinitePolyExtElement, Cache_base
from sage.rings.integer cimport Integer

cdef class FiniteField_ntl_gf2eElement(FinitePolyExtElement)

cdef class Cache_ntl_gf2e(Cache_base):
    cdef GF2EContext_c F
    cdef FiniteField _parent
    cdef public FiniteField_ntl_gf2eElement _zero_element
    cdef public FiniteField_ntl_gf2eElement _one_element
    cdef public FiniteField_ntl_gf2eElement _gen
    cdef Integer _order
    cdef Integer _degree
    cdef FiniteField_ntl_gf2eElement _new(self)

cdef class FiniteField_ntl_gf2eElement(FinitePolyExtElement):
    cdef GF2E_c x
    cdef Cache_ntl_gf2e _cache
    cdef FiniteField_ntl_gf2eElement _new(FiniteField_ntl_gf2eElement self)
