"""
This file provides the declaration for ``CRElement`` and the morphisms
to and from the integers and rationals.

It is included in the .pxd files associated to gluing files, such as
padic_capped_relative_element.pxd.

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
from sage.rings.morphism cimport RingHomomorphism_coercion, RingMap

cdef class CRElement(pAdicTemplateElement):
    cdef celement unit
    cdef long ordp
    cdef long relprec

    cdef CRElement _new_c(self)
    cdef int _normalize(self) except -1

cdef class pAdicCoercion_ZZ_CR(RingHomomorphism_coercion):
    cdef CRElement _zero
    cdef RingMap _section
cdef class pAdicConvert_CR_ZZ(RingMap):
    pass
cdef class pAdicCoercion_QQ_CR(RingHomomorphism_coercion):
    cdef CRElement _zero
    cdef RingMap _section
cdef class pAdicConvert_CR_QQ(RingMap):
    pass
cdef class pAdicConvert_QQ_CR(Morphism):
    cdef CRElement _zero
    cdef RingMap _section
