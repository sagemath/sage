"""
Example of a class wrapping an STL vector

EXAMPLES::

    sage: from sage.tests.stl_vector import stl_int_vector
    sage: v = stl_int_vector()
    sage: v
    A vector of integers
    vector<int>:
     data[0] = 123
     data[1] = 456

AUTHORS:

- Volker Braun (2012-01-18): initial version
"""

#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/interrupt.pxi'

from sage.structure.sage_object cimport SageObject
from sage.rings.integer cimport Integer
from sage.libs.gmp.mpz cimport mpz_add_ui
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef class stl_int_vector(SageObject):
    """
    Example class wrapping an STL vector

    EXAMPLES::

        sage: from sage.tests.stl_vector import stl_int_vector
        sage: v = stl_int_vector()
    """

    cdef vector[int] *data
    cdef string *name

    def __cinit__(self):
        """
        The Cython constructor.

        EXAMPLES::

            sage: from sage.tests.stl_vector import stl_int_vector
            sage: v = stl_int_vector()   # indirect doctest
            sage: TestSuite(v)
            Test suite for A vector of integers
            vector<int>:
             data[0] = 123
             data[1] = 456
        """
        self.data = new vector[int]()
        self.name = new string(<char*>"A vector of integers\n")
        self.data.push_back(123)
        self.data.push_back(456)

    def __dealloc__(self):
        """
        The Cython destructor.

        EXAMPLES::

            sage: from sage.tests.stl_vector import stl_int_vector
            sage: v = stl_int_vector()   # indirect doctest
        """
        del self.data
        del self.name

    def __getitem__(self, int i):
        """
        Return the ``i``-th element.

        EXAMPLES::

            sage: from sage.tests.stl_vector import stl_int_vector
            sage: v = stl_int_vector()
            sage: v[1]
            456
        """
        assert i>=0 and i<self.data.size()
        return self.data.at(i)

    def __repr__(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: from sage.tests.stl_vector import stl_int_vector
            sage: v = stl_int_vector()
            sage: v
            A vector of integers
            vector<int>:
             data[0] = 123
             data[1] = 456
        """
        s = self.name.c_str()
        s += 'vector<int>:\n'
        for i in range(0,self.data.size()):
            s += ' data['+str(i)+'] = '+str(self.data.at(i))+'\n'
        return s.strip()

    cpdef sum(self):
        """
        Add the elements.

        EXAMPLES::

            sage: from sage.tests.stl_vector import stl_int_vector
            sage: v = stl_int_vector()
            sage: v.sum()
            579
        """
        cdef Integer accumulator = Integer(0)
        cdef vector[int].iterator i = self.data.begin()
        sig_on()
        while i != self.data.end():
            mpz_add_ui(accumulator.value, accumulator.value, <int>(i[0]))
            i += 1
        sig_off()
        return accumulator

    def __cmp__(lhs, stl_int_vector rhs):
        """
        Compare with ``other``.

        EXAMPLES::

            sage: from sage.tests.stl_vector import stl_int_vector
            sage: u = stl_int_vector()
            sage: v = stl_int_vector()
            sage: cmp(u,v)
            0
        """
        cdef int c = cmp(lhs.data.size(), rhs.data.size())
        if c != 0:
            return c
        cdef vector[int].iterator lhs_iter = lhs.data.begin()
        cdef vector[int].iterator rhs_iter = rhs.data.begin()
        sig_on()
        try:
            while lhs_iter != lhs.data.end():
                c = cmp(<int>(lhs_iter[0]), <int>(rhs_iter[0]))
                if c != 0:
                    return c
                lhs_iter += 1
                rhs_iter += 1
        finally:
            sig_off()
        return 0

