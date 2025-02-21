# distutils: language = c++
"""
Dancing Links internal pyx code

EXAMPLES::

    sage: from sage.combinat.matrices.dancing_links import dlx_solver
    sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
    sage: x = dlx_solver(rows)
    sage: x
    Dancing links solver for 6 columns and 6 rows

The number of solutions::

    sage: x.number_of_solutions()
    3

Iterate over the solutions::

    sage: sorted(map(sorted, x.solutions_iterator()))
    [[0, 1], [2, 3], [4, 5]]

All solutions (computed in parallel)::

    sage: sorted(map(sorted, x.all_solutions()))
    [[0, 1], [2, 3], [4, 5]]

Return the first solution found when the computation is done in parallel::

    sage: sorted(x.one_solution(ncpus=2)) # random
    [0, 1]

Find all solutions using some specific rows::

    sage: x_using_row_2 = x.restrict([2])
    sage: x_using_row_2
    Dancing links solver for 7 columns and 6 rows
    sage: sorted(map(sorted, x_using_row_2.solutions_iterator()))
    [[2, 3]]

The two basic methods that are wrapped in this class are ``search`` which
returns ``1`` if a solution is found or ``0`` otherwise and ``get_solution``
which return the current solution::

    sage: x = dlx_solver(rows)
    sage: x.search()
    1
    sage: x.get_solution()
    [0, 1]
    sage: x.search()
    1
    sage: x.get_solution()
    [2, 3]
    sage: x.search()
    1
    sage: x.get_solution()
    [4, 5]
    sage: x.search()
    0

There is also a method ``reinitialize`` to reinitialize the algorithm::

    sage: x.reinitialize()
    sage: x.search()
    1
    sage: x.get_solution()
    [0, 1]
"""
# ****************************************************************************
#       Copyright (C) 2008 Carlo Hamalainen <carlo.hamalainen@gmail.com>
#       Copyright (C) 2015-2018 Sébastien Labbé <slabqc@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport PyObject_RichCompare
from libcpp.vector cimport vector
from cysignals.signals cimport sig_on, sig_off

cdef extern from "dancing_links_c.h":
    cdef cppclass dancing_links:
        dancing_links()
        vector[int] solution
        int number_of_columns()
        void add_rows(vector[vector[int]] rows)
        int search_is_started()
        int search()

from sage.misc.cachefunc import cached_method

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

        We must pass a list of rows at start up.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows+= [[0,2]]
            sage: rows+= [[1]]
            sage: rows+= [[3]]
            sage: x = dlx_solver(rows)
            sage: x
            Dancing links solver for 4 columns and 4 rows
            sage: x.search()
            1
            sage: x.get_solution()
            [3, 0]

        ::

            sage: rows = [[0,1,2], [1, 2]]
            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: x = dlx_solver(rows)
            sage: x
            Dancing links solver for 3 columns and 2 rows
            sage: type(x)
            <class 'sage.combinat.matrices.dancing_links.dancing_linksWrapper'>

        TESTS:

        The following example would crash in Sage's debug version
        from :issue:`13864` prior to the fix from :issue:`13882`::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: x = dlx_solver([])
            sage: x.get_solution()
            []
        """
        self._rows = [row for row in rows]
        self._initialize()

    def _initialize(self):
        r"""
        Initialization of the search algorithm.

        This adds the rows to the instance of dancing_links. This method is
        used by ``__init__`` and ``reinitialize`` methods and should not be
        used directly.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: x = dlx_solver(rows)         # indirect doctest
            sage: x.get_solution() if x.search() else None
            [0, 1]
            sage: x.get_solution() if x.search() else None
            [2, 3]

        Reinitialization of the algorithm::

            sage: x.reinitialize()             # indirect doctest
            sage: x.get_solution() if x.search() else None
            [0, 1]
        """
        cdef vector[int] v
        cdef vector[vector[int]] vv

        for row in self._rows:
            v.clear()
            for x in row:
                v.push_back(x)
            vv.push_back(v)

        sig_on()
        self._x.add_rows(vv)
        sig_off()

    def reinitialize(self):
        r"""
        Reinitialization of the search algorithm.

        This recreates an empty ``dancing_links`` object and adds the rows to
        the instance of ``dancing_links.``

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: x = dlx_solver(rows)
            sage: x.get_solution() if x.search() else None
            [0, 1]
            sage: x.get_solution() if x.search() else None
            [2, 3]

        Reinitialization of the algorithm::

            sage: x.reinitialize()
            sage: x.get_solution() if x.search() else None
            [0, 1]
            sage: x.get_solution() if x.search() else None
            [2, 3]
            sage: x.get_solution() if x.search() else None
            [4, 5]
            sage: x.get_solution() if x.search() else None

        Reinitialization works after solutions are exhausted::

            sage: x.reinitialize()
            sage: x.get_solution() if x.search() else None
            [0, 1]
            sage: x.get_solution() if x.search() else None
            [2, 3]
            sage: x.get_solution() if x.search() else None
            [4, 5]
            sage: x.get_solution() if x.search() else None
        """
        sig_on()
        self._x = dancing_links()
        sig_off()

        self._initialize()

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
                self.ncols(), self.nrows())

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

    def ncols(self):
        """
        Return the number of columns.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [1,2], [0], [3,4,5]]
            sage: dlx = dlx_solver(rows)
            sage: dlx.ncols()
            6
        """
        return self._x.number_of_columns()

    def nrows(self):
        """
        Return the number of rows.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [1,2], [0], [3,4,5]]
            sage: dlx = dlx_solver(rows)
            sage: dlx.nrows()
            4
        """
        return len(self._rows)

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

        Test that :issue:`11814` is fixed::

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

    def restrict(self, indices):
        r"""
        Return a dancing links solver solving the subcase which uses some
        given rows.

        For every row that is wanted in the solution, we add a new column
        to the row to make sure it is in the solution.

        INPUT:

        - ``indices`` -- list; row indices to be found in the solution

        OUTPUT: dancing links solver

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: d
            Dancing links solver for 6 columns and 6 rows
            sage: sorted(map(sorted, d.solutions_iterator()))
            [[0, 1], [2, 3], [4, 5]]

        To impose that the `0`-th row is part of the solution, the rows of the new
        problem are::

            sage: d_using_0 = d.restrict([0])
            sage: d_using_0
            Dancing links solver for 7 columns and 6 rows
            sage: d_using_0.rows()
            [[0, 1, 2, 6], [3, 4, 5], [0, 1], [2, 3, 4, 5], [0], [1, 2, 3, 4, 5]]

        After restriction the subproblem has one more columns and the same
        number of rows as the original one::

            sage: d.restrict([1]).rows()
            [[0, 1, 2], [3, 4, 5, 6], [0, 1], [2, 3, 4, 5], [0], [1, 2, 3, 4, 5]]
            sage: d.restrict([2]).rows()
            [[0, 1, 2], [3, 4, 5], [0, 1, 6], [2, 3, 4, 5], [0], [1, 2, 3, 4, 5]]

        This method allows to find solutions where the `0`-th row is part of a
        solution::

            sage: sorted(map(sorted, d.restrict([0]).solutions_iterator()))
            [[0, 1]]

        Some other examples::

            sage: sorted(map(sorted, d.restrict([1]).solutions_iterator()))
            [[0, 1]]
            sage: sorted(map(sorted, d.restrict([2]).solutions_iterator()))
            [[2, 3]]
            sage: sorted(map(sorted, d.restrict([3]).solutions_iterator()))
            [[2, 3]]

        Here there are no solutions using both 0th and 3rd row::

            sage: list(d.restrict([0,3]).solutions_iterator())
            []

        TESTS::

            sage: d.restrict([]).rows()
            [[0, 1, 2], [3, 4, 5], [0, 1], [2, 3, 4, 5], [0], [1, 2, 3, 4, 5]]
        """
        from copy import copy
        rows = copy(self._rows)
        ncols = self.ncols()
        for i, row_index in enumerate(indices):
            # in the line below we want the creation of a new list
            rows[row_index] = rows[row_index] + [ncols+i]
        return dlx_solver(rows)

    def split(self, column):
        r"""
        Return a dict of independent solvers.

        For each ``i``-th row containing a ``1`` in the ``column``, the
        dict associates the solver giving all solution using the ``i``-th
        row.

        This is used for parallel computations.

        INPUT:

        - ``column`` -- integer; the column used to split the problem into
          independent subproblems

        OUTPUT: dict where keys are row numbers and values are dlx solvers

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: d
            Dancing links solver for 6 columns and 6 rows
            sage: sorted(map(sorted, d.solutions_iterator()))
            [[0, 1], [2, 3], [4, 5]]

        After the split each subproblem has one more column and the same
        number of rows as the original problem::

            sage: D = d.split(0)
            sage: D
            {0: Dancing links solver for 7 columns and 6 rows,
             2: Dancing links solver for 7 columns and 6 rows,
             4: Dancing links solver for 7 columns and 6 rows}

        The (disjoint) union of the solutions of the subproblems is equal to the
        set of solutions shown above::

            sage: for x in D.values(): sorted(map(sorted, x.solutions_iterator()))
            [[0, 1]]
            [[2, 3]]
            [[4, 5]]

        TESTS::

            sage: d.split(6)
            Traceback (most recent call last):
            ...
            ValueError: column(=6) must be in range(ncols) where ncols=6

        This use to take a lot of time and memory. Not anymore since
        :issue:`24315`::

            sage: S = Subsets(range(11))
            sage: rows = map(list, S)
            sage: dlx = dlx_solver(rows)
            sage: dlx
            Dancing links solver for 11 columns and 2048 rows
            sage: d = dlx.split(0)
            sage: d[1]
            Dancing links solver for 12 columns and 2048 rows
        """
        if not 0 <= column < self.ncols():
            raise ValueError("column(={}) must be in range(ncols) "
                             "where ncols={}".format(column, self.ncols()))
        indices = (i for i, row in enumerate(self._rows) if column in row)
        return {i: self.restrict([i]) for i in indices}

    def solutions_iterator(self):
        r"""
        Return an iterator of the solutions.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: sorted(map(sorted, d.solutions_iterator()))
            [[0, 1], [2, 3], [4, 5]]

        TESTS:

        The algorithm is automatically reinitialized if needed, for example
        when iterating the solutions a second time (:issue:`25125`)::

            sage: sorted(map(sorted, d.solutions_iterator()))
            [[0, 1], [2, 3], [4, 5]]
        """
        if self._x.search_is_started():
            self.reinitialize()
        while self.search():
            yield self.get_solution()

    def one_solution(self, ncpus=None, column=None):
        r"""
        Return the first solution found.

        This method allows parallel computations which might be useful for
        some kind of problems when it is very hard just to find one
        solution.

        INPUT:

        - ``ncpus`` -- integer (default: ``None``); maximal number of
          subprocesses to use at the same time. If ``None``, it detects the
          number of effective CPUs in the system using
          :func:`sage.parallel.ncpus.ncpus()`.
          If ``ncpus=1``, the first solution is searched serially.
        - ``column`` -- integer (default: ``None``); the column used to split
          the problem (see :meth:`restrict`). If ``None``, a random column
          is chosen. This argument is ignored if ``ncpus=1``.

        OUTPUT: list of rows or ``None`` if no solution is found

        .. NOTE::

            For some case, increasing the number of cpus makes it
            faster. For other instances, ``ncpus=1`` is faster. It all
            depends on problem which is considered.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: solutions = [[0,1], [2,3], [4,5]]
            sage: sorted(d.one_solution()) in solutions
            True

        The number of CPUs can be specified as input::

            sage: sorted(d.one_solution(ncpus=2)) in solutions
            True

        The column used to split the problem for parallel computations can
        be given::

            sage: sorted(d.one_solution(ncpus=2, column=4)) in solutions
            True

        When no solution is found::

            sage: rows = [[0,1,2], [2,3,4,5], [0,1,2,3]]
            sage: d = dlx_solver(rows)
            sage: d.one_solution() is None
            True

        TESTS::

            sage: [d.one_solution(column=i) for i in range(6)]
            [None, None, None, None, None, None]

        The preprocess needed to start the parallel computation is not so
        big (less than 50ms in the example below)::

            sage: S = Subsets(range(11))
            sage: rows = list(map(list, S))
            sage: dlx = dlx_solver(rows)
            sage: dlx
            Dancing links solver for 11 columns and 2048 rows
            sage: solution = dlx.one_solution()
            sage: subsets = [set(rows[i]) for i in solution]

        We make sure the solution is an exact cover::

            sage: set.union(*subsets)
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
            sage: from itertools import combinations
            sage: any(p.intersection(q) for p,q in combinations(subsets, 2))
            False
        """
        if ncpus == 1:
            return self.get_solution() if self.search() else None

        if column is None:
            from random import randrange
            column = randrange(self.ncols())

        if not 0 <= column < self.ncols():
            raise ValueError("column(={}) must be in range(ncols) "
                             "where ncols={}".format(column, self.ncols()))

        from sage.parallel.decorate import parallel

        @parallel(ncpus=ncpus)
        def first_solution(i):
            dlx = self.restrict([i])
            if dlx.search():
                return dlx.get_solution()
            else:
                return None

        indices = [i for (i, row) in enumerate(self._rows) if column in row]
        for (args_kwds, val) in first_solution(indices):
            if val is not None:
                return val

    def all_solutions(self, ncpus=None, column=None):
        r"""
        Return all solutions found after splitting the problem to allow
        parallel computation.

        INPUT:

        - ``ncpus`` -- integer (default: ``None``); maximal number of
          subprocesses to use at the same time. If ``None``, it detects the
          number of effective CPUs in the system using
          :func:`sage.parallel.ncpus.ncpus()`.
        - ``column`` -- integer (default: ``None``); the column used to split
          the problem, if ``None`` a random column is chosen

        OUTPUT: list of solutions

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: S = d.all_solutions()
            sage: sorted(sorted(s) for s in S)
            [[0, 1], [2, 3], [4, 5]]

        The number of CPUs can be specified as input::

            sage: S = Subsets(range(4))
            sage: rows = map(list, S)
            sage: dlx = dlx_solver(rows)
            sage: dlx
            Dancing links solver for 4 columns and 16 rows
            sage: dlx.number_of_solutions()
            15
            sage: sorted(sorted(s) for s in dlx.all_solutions(ncpus=2))
            [[1, 2, 3, 4],
             [1, 2, 10],
             [1, 3, 9],
             [1, 4, 8],
             [1, 14],
             [2, 3, 7],
             [2, 4, 6],
             [2, 13],
             [3, 4, 5],
             [3, 12],
             [4, 11],
             [5, 10],
             [6, 9],
             [7, 8],
             [15]]

        If ``ncpus=1``, the computation is not done in parallel::

            sage: sorted(sorted(s) for s in dlx.all_solutions(ncpus=1))
            [[1, 2, 3, 4],
             [1, 2, 10],
             [1, 3, 9],
             [1, 4, 8],
             [1, 14],
             [2, 3, 7],
             [2, 4, 6],
             [2, 13],
             [3, 4, 5],
             [3, 12],
             [4, 11],
             [5, 10],
             [6, 9],
             [7, 8],
             [15]]

        TESTS:

        When no solution is found::

            sage: rows = [[0,1,2], [2,3,4,5], [0,1,2,3]]
            sage: d = dlx_solver(rows)
            sage: d.all_solutions()
            []

        ::

            sage: [d.all_solutions(column=i) for i in range(6)]
            [[], [], [], [], [], []]
        """
        if ncpus == 1:
            if self._x.search_is_started():
                self.reinitialize()
            L = []
            while self.search():
                L.append(self.get_solution())
            return L

        if column is None:
            from random import randrange
            column = randrange(self.ncols())

        if not 0 <= column < self.ncols():
            raise ValueError("column(={}) must be in range(ncols) "
                             "where ncols={}".format(column, self.ncols()))

        from sage.parallel.decorate import parallel

        @parallel(ncpus=ncpus)
        def all_solutions(i):
            dlx = self.restrict([i])
            L = []
            while dlx.search():
                L.append(dlx.get_solution())
            return L

        indices = [i for i, row in enumerate(self._rows) if column in row]
        L = []
        for (args_kwds, val) in all_solutions(indices):
            L.extend(val)
        return L

    def number_of_solutions(self, ncpus=None, column=None):
        r"""
        Return the number of distinct solutions.

        INPUT:

        - ``ncpus`` -- integer (default: ``None``); maximal number of
          subprocesses to use at the same time. If ``ncpus>1`` the dancing
          links problem is split into independent subproblems to allow
          parallel computation. If ``None``, it detects the number of
          effective CPUs in the system using
          :func:`sage.parallel.ncpus.ncpus()`.
        - ``column`` -- integer (default: ``None``); the column used to split
          the problem, if ``None`` a random column is chosen (this argument
          is ignored if ``ncpus`` is ``1``)

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2]]
            sage: rows += [[0,2]]
            sage: rows += [[1]]
            sage: rows += [[3]]
            sage: x = dlx_solver(rows)
            sage: x.number_of_solutions()
            2

        The number of CPUs can be specified as input::

            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: x = dlx_solver(rows)
            sage: x.number_of_solutions(ncpus=2, column=3)
            3

        ::

            sage: S = Subsets(range(5))
            sage: rows = map(list, S)
            sage: d = dlx_solver(rows)
            sage: d.number_of_solutions()
            52

        TESTS:

        The algorithm is automatically reinitialized if needed, for example
        when counting the number of solutions a second time (:issue:`25125`)::

            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: x = dlx_solver(rows)
            sage: x.number_of_solutions(ncpus=1)
            3
            sage: x.number_of_solutions(ncpus=1)
            3

        Works with empty rows::

            sage: dlx_solver([]).number_of_solutions(ncpus=None)
            0
            sage: dlx_solver([]).number_of_solutions(ncpus=1)
            0
        """
        cdef int N = 0
        if ncpus == 1:
            if self._x.search_is_started():
                self.reinitialize()
            while self.search():
                N += 1
            return N

        if self.ncols() == 0:
            return 0

        if column is None:
            from random import randrange
            column = randrange(self.ncols())

        if not 0 <= column < self.ncols():
            raise ValueError("column(={}) must be in range(ncols) "
                             "where ncols={}".format(column, self.ncols()))

        from sage.parallel.decorate import parallel

        @parallel(ncpus=ncpus)
        def nb_sol(i):
            dlx = self.restrict([i])
            N = 0
            while dlx.search():
                N += 1
            return N

        indices = [i for i, row in enumerate(self._rows) if column in row]
        return sum(val for (args_kwds, val) in nb_sol(indices))

    @cached_method
    def to_sat_solver(self, solver=None):
        r"""
        Return the SAT solver solving an equivalent problem.

        Note that row index `i` in the dancing links solver corresponds to
        the boolean variable index `ì+1` for the SAT solver to avoid
        the variable index `0`.

        See also :mod:`sage.sat.solvers.satsolver`.

        INPUT:

        - ``solver`` -- string or ``None`` (default: ``None``),
          possible values include ``'picosat'``, ``'cryptominisat'``,
          ``'LP'``, ``'glucose'``, ``'glucose-syrup'``.

        OUTPUT: SAT solver instance

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [0,2], [1], [3]]
            sage: x = dlx_solver(rows)
            sage: s = x.to_sat_solver()                                                 # needs sage.sat

        Using some optional SAT solvers::

            sage: x.to_sat_solver('cryptominisat')      # optional - pycryptosat        # needs sage.sat
            CryptoMiniSat solver: 4 variables, 7 clauses.
        """
        from sage.sat.solvers.satsolver import SAT
        s = SAT(solver)

        # Note that row number i is associated to SAT variable i+1 to
        # avoid a variable zero
        columns = [[] for _ in range(self.ncols())]
        for i, row in enumerate(self.rows(), start=1):
            for a in row:
                columns[a].append(i)

        # At least one 1 in each column
        for clause in columns:
            s.add_clause(clause)

        # At most one 1 in each column
        import itertools
        for clause in columns:
            for p, q in itertools.combinations(clause, 2):
                sub_clause = [-p, -q]
                s.add_clause(sub_clause)

        return s

    def one_solution_using_sat_solver(self, solver=None):
        r"""
        Return a solution found using a SAT solver.

        INPUT:

        - ``solver`` -- string or ``None`` (default: ``None``),
          possible values include ``'picosat'``, ``'cryptominisat'``,
          ``'LP'``, ``'glucose'``, ``'glucose-syrup'``.

        OUTPUT: list of rows or ``None`` if no solution is found

        .. NOTE::

            When comparing the time taken by method ``one_solution``,
            have in mind that ``one_solution_using_sat_solver`` first
            creates the SAT solver instance from the dancing links
            solver. This copy of data may take many seconds depending on
            the size of the problem.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: solutions = [[0,1], [2,3], [4,5]]
            sage: d.one_solution_using_sat_solver() in solutions                        # needs sage.sat
            True

        Using optional solvers::

            sage: s = d.one_solution_using_sat_solver('glucose')                # optional - glucose, needs sage.sat
            sage: s in solutions                                                # optional - glucose, needs sage.sat
            True

        When no solution is found::

            sage: rows = [[0,1,2], [2,3,4,5], [0,1,2,3]]
            sage: d = dlx_solver(rows)
            sage: d.one_solution_using_sat_solver() is None                             # needs sage.sat
            True
        """
        sat_solver = self.to_sat_solver(solver)
        solution = sat_solver()
        if not solution:
            return None
        return [key for key, val in enumerate(solution, start=-1) if val]

    @cached_method
    def to_milp(self, solver=None):
        r"""
        Return the mixed integer linear program (MILP) representing an
        equivalent problem.

        See also :mod:`sage.numerical.mip.MixedIntegerLinearProgram`.

        INPUT:

        - ``solver`` -- string or ``None`` (default: ``None``); possible
          values include ``'GLPK'``, ``'GLPK/exact'``, ``'Coin'``,
          ``'CPLEX'``, ``'Gurobi'``, ``'CVXOPT'``, ``'PPL'``,
          ``'InteractiveLP'``

        OUTPUT:

        - MixedIntegerLinearProgram instance
        - MIPVariable with binary components

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [0,2], [1], [3]]
            sage: d = dlx_solver(rows)
            sage: p,x = d.to_milp()                                                     # needs sage.numerical.mip
            sage: p                                                                     # needs sage.numerical.mip
            Boolean Program (no objective, 4 variables, ... constraints)
            sage: x                                                                     # needs sage.numerical.mip
            MIPVariable with 4 binary components

        In the reduction, the boolean variable `x_i` is ``True`` if and only if
        the `i`-th row is in the solution::

            sage: p.show()                                                              # needs sage.numerical.mip
            Maximization:
            <BLANKLINE>
            <BLANKLINE>
            Constraints:...
              one 1 in 0-th column: 1.0 <= x_0 + x_1 <= 1.0
              one 1 in 1-th column: 1.0 <= x_0 + x_2 <= 1.0
              one 1 in 2-th column: 1.0 <= x_0 + x_1 <= 1.0
              one 1 in 3-th column: 1.0 <= x_3 <= 1.0
            Variables:
              x_0 is a boolean variable (min=0.0, max=1.0)
              x_1 is a boolean variable (min=0.0, max=1.0)
              x_2 is a boolean variable (min=0.0, max=1.0)
              x_3 is a boolean variable (min=0.0, max=1.0)

        Using some optional MILP solvers::

            sage: d.to_milp('gurobi')           # optional - gurobi sage_numerical_backends_gurobi, needs sage.numerical.mip
            (Boolean Program (no objective, 4 variables, 4 constraints),
             MIPVariable with 4 binary components)
        """
        from sage.numerical.mip import MixedIntegerLinearProgram
        p = MixedIntegerLinearProgram(solver=solver)

        # x[i] == True iff i-th dlx row is in the solution
        x = p.new_variable(binary=True, indices=range(self.nrows()))

        # Construction of the columns (transpose of the rows)
        columns = [[] for _ in range(self.ncols())]
        for i, row in enumerate(self.rows()):
            for a in row:
                columns[a].append(i)

        # Constraints: exactly one 1 in each column
        for j, column in enumerate(columns):
            S = p.sum(x[a] for a in column)
            name = "one 1 in {}-th column".format(j)
            p.add_constraint(S == 1, name=name)

        return p, x

    def one_solution_using_milp_solver(self, solver=None, integrality_tolerance=1e-3):
        r"""
        Return a solution found using a MILP solver.

        INPUT:

        - ``solver`` -- string or ``None`` (default: ``None``); possible
          values include ``'GLPK'``, ``'GLPK/exact'``, ``'Coin'``,
          ``'CPLEX'``, ``'Gurobi'``, ``'CVXOPT'``, ``'PPL'``,
          ``'InteractiveLP'``

        OUTPUT: list of rows or ``None`` if no solution is found

        .. NOTE::

            When comparing the time taken by method ``one_solution``, have in
            mind that ``one_solution_using_milp_solver`` first creates (and
            caches) the MILP solver instance from the dancing links solver.
            This copy of data may take many seconds depending on the size
            of the problem.

        EXAMPLES::

            sage: from sage.combinat.matrices.dancing_links import dlx_solver
            sage: rows = [[0,1,2], [3,4,5], [0,1], [2,3,4,5], [0], [1,2,3,4,5]]
            sage: d = dlx_solver(rows)
            sage: solutions = [[0,1], [2,3], [4,5]]
            sage: d.one_solution_using_milp_solver() in solutions                       # needs sage.numerical.mip
            True

        Using optional solvers::

            sage: # optional - gurobi sage_numerical_backends_gurobi, needs sage.numerical.mip
            sage: s = d.one_solution_using_milp_solver('gurobi')
            sage: s in solutions
            True

        When no solution is found::

            sage: rows = [[0,1,2], [2,3,4,5], [0,1,2,3]]
            sage: d = dlx_solver(rows)
            sage: d.one_solution_using_milp_solver() is None                            # needs sage.numerical.mip
            True
        """
        from sage.numerical.mip import MIPSolverException
        p, x = self.to_milp(solver)
        try:
            p.solve()
        except MIPSolverException:
            return None

        soln = p.get_values(x, convert=bool, tolerance=integrality_tolerance)
        return sorted(key for key in soln if soln[key])


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
    from sage.misc.persist import loads
    return dancing_linksWrapper(loads(s))
