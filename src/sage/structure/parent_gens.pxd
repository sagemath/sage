"""
Parent objects with generators
"""

###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

cimport parent_base

cdef class ParentWithGens(parent_base.ParentWithBase):
    cdef public object _gens
    cdef public object _latex_names
    cdef public object _list


cdef class ParentWithMultiplicativeAbelianGens(ParentWithGens):
    cdef public object _generator_orders

cdef class ParentWithAdditiveAbelianGens(ParentWithGens):
    cdef public object _generator_orders
