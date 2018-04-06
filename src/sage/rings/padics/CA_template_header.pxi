"""
This file provides the declaration for ``CAElement`` and the morphisms to
the integers and rationals.

It is included in the .pxd files associated to gluing files, such as
padic_capped_absolute_element.pxd.

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

cdef class CAElement(pAdicTemplateElement):
    cdef celement value
    cdef long absprec

    cdef CAElement _new_c(self)

cdef class pAdicCoercion_ZZ_CA(RingHomomorphism):
    cdef CAElement _zero
    cdef RingMap _section
cdef class pAdicConvert_CA_ZZ(RingMap):
    pass
cdef class pAdicConvert_QQ_CA(Morphism):
    cdef CAElement _zero
    cdef RingMap _section
# There should also be a pAdicConvert_CA_QQ for extension rings....
cdef class pAdicCoercion_CA_frac_field(RingHomomorphism):
    cdef CRElement _zero
    cdef Morphism _section
cdef class pAdicConvert_CA_frac_field(Morphism):
    cdef CAElement _zero
