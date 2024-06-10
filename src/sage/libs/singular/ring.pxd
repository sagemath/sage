"""
Wrapper for Singular's Rings

AUTHOR:

- Martin Albrecht (2009-07): initial implementation
"""
#*****************************************************************************
#       Copyright (C) 2009 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.singular.decl cimport ring


# To work with singular rings, you need to balance singular_ring_new with
# singular_ring_delete or singular_ring_reference with
# singular_ring_delete. That is, either use one of the two patterns:
#
# cdef class myclass_new():
#     cdef ring* myring;
#     def __cinit__():
#         self.myring = singular_ring_new(...)
#     def __dealloc__():
#         singular_ring_delete(self.myring)
#
# cdef class myclass_reference():
#     cdef ring* refring;
#     def __cinit__(ring* some_ring):
#         self.refring = singular_ring_reference(some_ring)
#     def __dealloc__():
#         singular_ring_delete(self.refring)
#
# You must not refer to Python/Cython classes in the Cython
# destructor, the following is INVALID:
#
# cdef class myclass_invalid():
#     cdef Parent parent;
#     def __cinit__(Parent p):
#         self.parent = p
#     def __dealloc__():
#         do_something_with(self.parent.ring)   # segfault


# create a new singular ring
cdef ring *singular_ring_new(base_ring, n, names, term_order) except NULL

# reference an existing ring
cdef ring *singular_ring_reference(ring *existing_ring) except NULL

# carefully delete a ring once its refcount is zero
cdef void singular_ring_delete(ring *doomed) noexcept

# Used internally for reference counting
cdef wrap_ring(ring* R)
