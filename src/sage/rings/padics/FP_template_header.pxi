# sage_setup: distribution = sagemath-pari
"""
This file provides the declaration for ``FPElement`` and the morphisms
to and from the integers and rationals.

It is included in the .pxd files associated to gluing files, such as
padic_floating_point_element.pxd.

AUTHORS:

- David Roe (2012-03-01) -- initial version
"""

#*****************************************************************************
#       Copyright (C) 2012 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# includes the header for pAdicTemplateElement
include "padic_template_element_header.pxi"

from sage.categories.morphism cimport Morphism
from sage.rings.morphism cimport RingHomomorphism, RingMap

cdef class FPElement(pAdicTemplateElement):
    cdef celement unit
    cdef long ordp

    cdef FPElement _new_c(self)
    cdef int _normalize(self) except -1
    cdef int _set_infinity(self) except -1
    cpdef val_unit(self, p=*)

cdef class pAdicCoercion_ZZ_FP(RingHomomorphism):
    cdef FPElement _zero
    cdef RingMap _section
cdef class pAdicConvert_FP_ZZ(RingMap):
    pass
cdef class pAdicCoercion_QQ_FP(RingHomomorphism):
    cdef FPElement _zero
    cdef RingMap _section
cdef class pAdicConvert_FP_QQ(RingMap):
    pass
cdef class pAdicConvert_QQ_FP(Morphism):
    cdef FPElement _zero
    cdef RingMap _section
cdef class pAdicCoercion_FP_frac_field(RingHomomorphism):
    cdef FPElement _zero
    cdef Morphism _section
cdef class pAdicConvert_FP_frac_field(Morphism):
    cdef FPElement _zero
