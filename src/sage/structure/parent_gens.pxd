# sage_setup: distribution = sagemath-objects
"""
Parent objects with generators
"""

#*****************************************************************************
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent_base cimport ParentWithBase


cdef class ParentWithGens(ParentWithBase):
    cdef public object _gens
    cdef public object _latex_names
    cdef public object _list
