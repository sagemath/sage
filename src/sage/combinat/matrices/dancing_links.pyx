"""
Dancing Links internal pyx code
"""
#*****************************************************************************
#       Copyright (C) 2008 Carlo Hamalainen <carlo.hamalainen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#clang c++

import sys

from cpython.list cimport *
from cpython.int cimport *
from cpython.ref cimport *

cdef extern from "dancing_links_c.h":
    ctypedef struct vector_int "std::vector<int>":
        void (* push_back)(int elem)
        void clear()
        int at(size_t loc)
        int size()

    ctypedef struct vector_vector_int "std::vector<vector<int> >":
        void (* push_back)(vector_int elem)

    ctypedef struct dancing_links:
        vector_int solution
        void add_rows(vector_vector_int rows)
        int search()
        void freemem()

cdef extern from "ccobject.h":
    dancing_links* dancing_links_construct "Construct<dancing_links>"(void *mem)
    void dancing_links_destruct "Destruct<dancing_links>"(dancing_links *mem)


from sage.rings.integer cimport Integer

cdef class dancing_linksWrapper:
    cdef dancing_links x
    cdef object rows

    def __init__(self, rows):
        """
        Initialize our wrapper (self.x) as an actual C++ object.

        We must pass a list of rows at start up. There are no methods
        for resetting the list of rows, so this class acts as a one-time
        executor of the C++ code.

        TESTS:
            sage: rows = [[0,1,2], [1, 2]]
            sage: x = make_dlxwrapper(dumps(rows))
            sage: loads(x.__reduce__()[1][0])
            [[0, 1, 2], [1, 2]]

        """
        pass

    # Note that the parameters to __cinit__, if any, must be identical to __init__
    # This is due to the way Python constructs class instance.
    def __cinit__(self, rows):
        self.rows = PyList_New(len(rows))
        dancing_links_construct(&self.x)
        if rows:
            self.add_rows(rows)

    def __dealloc__(self):
        self.x.freemem()
        dancing_links_destruct(&self.x)

    def __str__(self):
        """
        The string representation of this wrapper is just the list of
        rows as supplied at startup.

        TESTS:
            sage: rows = [[0,1,2]]
            sage: print make_dlxwrapper(dumps(rows)).__str__()
            [[0, 1, 2]]
        """

        return self.rows.__str__()

    def __reduce__(self):
        """
        This is used when pickling.

        TESTS:
            sage: rows = [[0,1,2]]
            sage: x = make_dlxwrapper(dumps(rows))
            sage: loads(x.__reduce__()[1][0])
            [[0, 1, 2]]
        """
        # A comment from sage/rings/integer.pyx:

        # This single line below took me HOURS to figure out.
        # It is the *trick* needed to pickle Cython extension types.
        # The trick is that you must put a pure Python function
        # as the first argument, and that function must return
        # the result of unpickling with the argument in the second
        # tuple as input. All kinds of problems happen
        # if we don't do this.
        #return sage.rings.integer.make_integer, (self.str(32),)

        import sage.combinat.matrices.dancing_links
        from sage.all import dumps
        return sage.combinat.matrices.dancing_links.make_dlxwrapper, (dumps(self.rows),)

    def __richcmp__(dancing_linksWrapper left, dancing_linksWrapper right, int op):
        """
        Two dancing_linksWrapper objects are equal if they were
        initialised using the same row list.

        TESTS:
            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: X = dlx_solver(rows)
            sage: Z = dlx_solver(rows)
            sage: rows += [[2]]
            sage: Y = dlx_solver(rows)
            sage: X == Z
            1
            sage: X == Y
            0
        """

        cdef int equal
        equal = left.rows == right.rows

        if op == 2: # ==
            return equal
        elif op == 3: # !=
            return not equal
        else:
            return NotImplemented

    def dumps(self):
        """
        TESTS:
            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: X = dlx_solver(rows)
            sage: X == loads(dumps(X))
            1
            sage: rows += [[2]]
            sage: Y = dlx_solver(rows)
            sage: Y == loads(dumps(X))
            0
        """
        return self.rows.dumps()

    def add_rows(self, rows):
        """
        Initialize our instance of dancing_links with the given rows.
        This is for internal use by dlx_solver.

        This doctest tests add_rows vicariously!

        TESTS:
            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows+= [[0,2]]
            sage: rows+= [[1]]
            sage: rows+= [[3]]
            sage: x = dlx_solver(rows)
            sage: print x.search()
            1

        The following example would crash in Sage's debug version
        from :trac:`13864` prior to the fix from :trac:`13822`::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: x = dlx_solver([])          # indirect doctest
            sage: x.get_solution()
            []

        """
        if not rows:
            return

        cdef vector_int v
        cdef vector_vector_int vv

        cdef int i = 0

        for row in rows:
            v.clear()

            Py_INCREF(row);
            PyList_SET_ITEM(self.rows, i, row)
            i += 1

            for x in row:
                v.push_back(x)

            vv.push_back(v)

        self.x.add_rows(vv)

    def get_solution(self):
        """
        After calling search(), we can extract a solution
        from the instance variable self.x.solution, a C++ vector<int>
        listing the rows that make up the current solution.

        TESTS:
            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows+= [[0,2]]
            sage: rows+= [[1]]
            sage: rows+= [[3]]
            sage: x = dlx_solver(rows)
            sage: print x.search()
            1
            sage: print x.get_solution()
            [3, 0]
        """

        s = []
        for i in range(self.x.solution.size()):
            s.append(self.x.solution.at(i))

        return s

    def search(self):
        """
        TESTS:
            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows+= [[0,2]]
            sage: rows+= [[1]]
            sage: rows+= [[3]]
            sage: x = dlx_solver(rows)
            sage: print x.search()
            1
            sage: print x.get_solution()
            [3, 0]
        """

        x = self.x.search()

        return x

def dlx_solver(rows):
    """
    Internal-use wrapper for the dancing links C++ code.

    EXAMPLES:
        sage: from sage.combinat.matrices.dancing_links import dlx_solver
        sage: rows = [[0,1,2]]
        sage: rows+= [[0,2]]
        sage: rows+= [[1]]
        sage: rows+= [[3]]
        sage: x = dlx_solver(rows)
        sage: print x.search()
        1
        sage: print x.get_solution()
        [3, 0]
        sage: print x.search()
        1
        sage: print x.get_solution()
        [3, 1, 2]
        sage: print x.search()
        0
    """

    cdef dancing_linksWrapper dlw

    dlw = dancing_linksWrapper(rows)
    #dlw.add_rows(rows)

    return dlw


def make_dlxwrapper(s):
    """
    Create a dlx wrapper from a Python *string* s.
    This is used in unpickling. We expect s to be dumps(rows) where
    rows is the list of rows used to instantiate the object.

    TESTS:
        sage: rows = [[0,1,2]]
        sage: x = make_dlxwrapper(dumps(rows))
        sage: print x.__str__()
        [[0, 1, 2]]
    """

    from sage.all import loads

    cdef dancing_linksWrapper dlw
    dlw = dancing_linksWrapper(loads(s))
    return dlw

