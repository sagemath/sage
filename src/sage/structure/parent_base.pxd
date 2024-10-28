# sage_setup: distribution = sagemath-objects
###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
###############################################################################

from sage.structure.parent_old cimport Parent as Parent_old

cdef class ParentWithBase(Parent_old):
    pass
