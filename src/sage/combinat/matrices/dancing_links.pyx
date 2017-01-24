# distutils: language = c++
"""
Dancing Links internal pyx code
"""

#*****************************************************************************
#       Copyright (C) 2008 Carlo Hamalainen <carlo.hamalainen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

include "cysignals/signals.pxi"
from cpython.object cimport PyObject_RichCompare

from libcpp.vector cimport vector

cdef extern from "dancing_links_c.h":
    ctypedef struct dancing_links:
        vector[int] solution
        int number_of_columns()
        void add_rows(vector[vector[int]] rows)
        int search()
        void freemem()

cdef extern from "ccobject.h":
    dancing_links* dancing_links_construct "Construct<dancing_links>"(void *mem)
    void dancing_links_destruct "Destruct<dancing_links>"(dancing_links *mem)

cdef class dancing_linksWrapper:
    r"""
    A simple class that implements dancing links.

    The main methods to list the solutions are :meth:`search` and
    :meth:`get_solution`. You can also use :meth:`number_of_solutions` to count
    them.

    This class simply wraps a C++ implementation of Carlo Hamalainen.
    """
    cdef dancing_links _x
    cdef list _rows

    def __init__(self, rows):
        """
        Initialize our wrapper (self._x) as an actual C++ object.

        We must pass a list of rows at start up. There are no methods
        for resetting the list of rows, so this class acts as a one-time
        executor of the C++ code.

        TESTS::

            sage: rows = [[0,1,2], [1, 2]]
            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: x = dlx_solver(rows)
            sage: x
            Dancing links solver for 3 columns and 2 rows
            sage: type(x)
            <type 'sage.combinat.matrices.dancing_links.dancing_linksWrapper'>
        """
        self._init_rows(rows)

    def __cinit__(self):
        dancing_links_construct(&self._x)

    def __dealloc__(self):
        self._x.freemem()
        dancing_links_destruct(&self._x)

    def __repr__(self):
        """
        The string representation of this wrapper is just the list of
        rows as supplied at startup.

        TESTS::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [1,2], [0]]
            sage: dlx_solver(rows)
            Dancing links solver for 3 columns and 3 rows
        """
        return "Dancing links solver for {} columns and {} rows".format(
                self._x.number_of_columns(),
                len(self._rows))

    def rows(self):
        r"""
        Return the list of rows.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [1,2], [0]]
            sage: x = dlx_solver(rows)
            sage: x.rows()
            [[0, 1, 2], [1, 2], [0]]
        """
        return self._rows

    def __reduce__(self):
        """
        This is used when pickling.

        TESTS::

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
        return type(self), (self._rows,)

    def __richcmp__(dancing_linksWrapper left, dancing_linksWrapper right, int op):
        """
        Two dancing_linksWrapper objects are equal if they were
        initialised using the same row list.

        TESTS::

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
        return PyObject_RichCompare(left._rows, right._rows, op)

    def _init_rows(self, rows):
        """
        Initialize our instance of dancing_links with the given rows.

        This is for internal use by dlx_solver only.

        TESTS:

        This doctest tests ``_init_rows`` vicariously! ::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows+= [[0,2]]
            sage: rows+= [[1]]
            sage: rows+= [[3]]
            sage: x = dlx_solver(rows)
            sage: print(x.search())
            1

        The following example would crash in Sage's debug version
        from :trac:`13864` prior to the fix from :trac:`13882`::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: x = dlx_solver([])          # indirect doctest
            sage: x.get_solution()
            []

        """
        cdef vector[int] v
        cdef vector[vector[int]] vv

        self._rows = [row for row in rows]

        for row in self._rows:
            v.clear()

            for x in row:
                v.push_back(x)

            vv.push_back(v)

        sig_on()
        self._x.add_rows(vv)
        sig_off()

    def get_solution(self):
        """
        Return the current solution.

        After a new solution is found using the method :meth:`search` this
        method return the rows that make up the current solution.

        TESTS::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows+= [[0,2]]
            sage: rows+= [[1]]
            sage: rows+= [[3]]
            sage: x = dlx_solver(rows)
            sage: print(x.search())
            1
            sage: print(x.get_solution())
            [3, 0]
        """
        cdef size_t i
        cdef list s = []
        for i in range(self._x.solution.size()):
            s.append(self._x.solution.at(i))

        return s

    def search(self):
        """
        Search for a new solution.

        Return ``1`` if a new solution is found and ``0`` otherwise. To recover
        the solution, use the method :meth:`get_solution`.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows+= [[0,2]]
            sage: rows+= [[1]]
            sage: rows+= [[3]]
            sage: x = dlx_solver(rows)
            sage: print(x.search())
            1
            sage: print(x.get_solution())
            [3, 0]

        TESTS:

        Test that :trac:`11814` is fixed::

            sage: dlx_solver([]).search()
            0
            sage: dlx_solver([[]]).search()
            0

        If search is called once too often, it keeps returning 0::

            sage: x = dlx_solver([[0]])
            sage: x.search()
            1
            sage: x.search()
            0
            sage: x.search()
            0
        """
        sig_on()
        x = self._x.search()
        sig_off()
        return x

    def split(self, column):
        r"""
        Return a dict of rows solving independent cases.

        For each ``i``-th row containing a ``1`` in the ``column``, the
        dict associates the list of rows where each other row containing a
        ``1`` in the same ``column`` is replaced by the empty list.

        This is used for parallel computations.

        INPUT:

        - ``column`` -- integer, the column used to split the problem

        OUTPUT:

            dict where keys are row numbers and values are lists of rows

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: d
            Dancing links solver for 6 columns and 6 rows
            sage: sorted(d.solutions_iterator())
            [[0, 1], [2, 3], [4, 5]]

        After the split each subproblem has the same number of columns and
        rows as the original one::

            sage: D = d.split(0)
            sage: D
            {0: [[0, 1, 2], [3, 4, 5], [], [2, 3, 4, 5], [], [1, 2, 3, 4, 5]],
             2: [[], [3, 4, 5], [0, 1], [2, 3, 4, 5], [], [1, 2, 3, 4, 5]],
             4: [[], [3, 4, 5], [], [2, 3, 4, 5], [0], [1, 2, 3, 4, 5]]}

        The (disjoint) union of the solutions of the subproblems is equal to the
        set of solutions shown above::

            sage: for a in D.values(): list(dlx_solver(a).solutions_iterator())
            [[0, 1]]
            [[2, 3]]
            [[4, 5]]
        """
        from copy import deepcopy
        indices = [i for (i,row) in enumerate(self._rows) if column in row]
        D = dict()
        for i in indices:
            rows = deepcopy(self._rows)
            for j in indices:
                if j != i:
                    rows[j] = []
            D[i] = rows
        return D

    def solutions_iterator(self):
        r"""
        Return an iterator of the solutions.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: list(d.solutions_iterator())
            [[0, 1], [2, 3], [4, 5]]
        """
        while self.search():
            yield self.get_solution()

    def _number_of_solutions_iterator(self, ncpus=1, column=0):
        r"""
        Return an iterator over the number of solutions using each row
        containing a ``1`` in the given ``column``.

        INPUT:

        - ``ncpus`` -- integer (default: ``1``), maximal number of
          subprocesses to use at the same time
        - ``column`` -- integer (default: ``0``), the column used to split
          the problem

        OUTPUT:

            iterator of tuples (row number, number of solutions)

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: sorted(d._number_of_solutions_iterator(ncpus=2, column=3))
            [(1, 1), (3, 1), (5, 1)]

        ::

            sage: S = Subsets(range(5))
            sage: rows = map(list, S)
            sage: d = dlx_solver(rows)
            sage: d.number_of_solutions()
            52
            sage: sum(b for a,b in d._number_of_solutions_iterator(ncpus=2, column=3))
            52
        """
        D = self.split(column)
        from sage.parallel.decorate import parallel

        @parallel(ncpus=ncpus)
        def nb_sol(i):
            return dlx_solver(D[i]).number_of_solutions()
        K = sorted(D)
        for ((args, kwds), val) in nb_sol(K):
            yield args[0], val

    def number_of_solutions(self, ncpus=1, column=0):
        r"""
        Return the number of distinct solutions.

        INPUT:

        - ``ncpus`` -- integer (default: ``1``), maximal number of
          subprocesses to use at the same time. If `ncpus>1` the dancing
          links problem is split into independent subproblems to
          allow parallel computation.
        - ``column`` -- integer (default: ``0``), the column used to split
          the problem (ignored if ``ncpus`` is ``1``)

        OUTPUT:

        integer

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows += [[0,2]]
            sage: rows += [[1]]
            sage: rows += [[3]]
            sage: x = dlx_solver(rows)
            sage: x.number_of_solutions()
            2

        ::

            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: x = dlx_solver(rows)
            sage: x.number_of_solutions(ncpus=2, column=3)
            3

        TESTS::

            sage: dlx_solver([]).number_of_solutions()
            0
        """
        cdef int N = 0
        if ncpus == 1:
            while self.search():
                N += 1
            return N
        else:
            it = self._number_of_solutions_iterator(ncpus, column)
            return sum(val for (k,val) in it)

def dlx_solver(rows):
    """
    Internal-use wrapper for the dancing links C++ code.

    EXAMPLES::

        sage: from sage.combinat.matrices.dancing_links import dlx_solver
        sage: rows = [[0,1,2]]
        sage: rows+= [[0,2]]
        sage: rows+= [[1]]
        sage: rows+= [[3]]
        sage: x = dlx_solver(rows)
        sage: print(x.search())
        1
        sage: print(x.get_solution())
        [3, 0]
        sage: print(x.search())
        1
        sage: print(x.get_solution())
        [3, 1, 2]
        sage: print(x.search())
        0
    """
    return dancing_linksWrapper(rows)


def make_dlxwrapper(s):
    """
    Create a dlx wrapper from a Python *string* s.

    This was historically used in unpickling and is kept for backwards
    compatibility. We expect s to be ``dumps(rows)`` where rows is the
    list of rows used to instantiate the object.

    TESTS::

        sage: from sage.combinat.matrices.dancing_links import make_dlxwrapper
        sage: rows = [[0,1,2]]
        sage: x = make_dlxwrapper(dumps(rows))
        sage: print(x.__str__())
        Dancing links solver for 3 columns and 1 rows
    """
    from sage.all import loads
    return dancing_linksWrapper(loads(s))
