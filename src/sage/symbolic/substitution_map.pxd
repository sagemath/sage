from ginac cimport *

from sage.symbolic.expression cimport Expression
from sage.structure.sage_object cimport SageObject

cdef class SubstitutionMap(SageObject):
    cdef GExMap _gmapobj
    cpdef Expression apply_to(self, Expression expr, unsigned options)
    
cdef SubstitutionMap new_SubstitutionMap_from_GExMap(const GExMap& smap)

