# sage_setup: distribution = sagemath-objects
"""
Mutability -- Pyrex Implementation
"""

##########################################################################
#
#   Sage: Open Source Mathematical Software
#
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
##########################################################################

cdef class Mutability:
    cdef public bint _is_immutable
    cpdef _require_mutable(self)
    cpdef _require_immutable(self)
    cpdef bint is_immutable(self) noexcept
    cpdef bint is_mutable(self) noexcept
