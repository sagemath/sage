"""
This file provides the declaration for the pAdicTemplateElement class,
which collects common functionality for the differen p-adic template
classes.

It is included in CR_template_header.pxi, CA_template_header.pxi and
FM_template_header.pxi.  Each of these are then included in the .pxd
files associated to gluing files, such as
padic_capped_relative_element.pxd.

.. NOTE::

    Since this file is included into others at compile time, there
    will be multiple pAdicTemplateElement classes: there is no common
    one to cimport.

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

from sage.structure.element cimport ModuleElement, RingElement
from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.rings.padics.pow_computer cimport PowComputer_class

cdef class pAdicTemplateElement(pAdicGenericElement):
    cdef PowComputer_class prime_pow
    cdef int _set(self, x, long val, long xprec, absprec, relprec) except -1
    cdef pAdicTemplateElement _lshift_c(self, long shift)
    cdef pAdicTemplateElement _rshift_c(self, long shift)
    #cpdef RingElement _floordiv_c_impl(self, RingElement right)
    cdef int check_preccap(self) except -1
    cdef pAdicTemplateElement lift_to_precision_c(self, long absprec)
    cpdef pAdicTemplateElement unit_part(self)
