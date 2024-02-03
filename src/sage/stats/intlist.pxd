# sage_setup: distribution = sagemath-modules
#############################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL) v2+.
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

cdef class IntList:
    cdef int* _values
    cdef Py_ssize_t _length

    cpdef int prod(self) noexcept
    cpdef int sum(self) noexcept
