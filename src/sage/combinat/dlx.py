"""
Exact Cover Problem via Dancing Links
"""
# dlx.py
# Copyright (c) 2006,2008 Antti Ajanki <antti.ajanki@iki.fi>

# Heavily Modified Feb 2008, Tom Boothby
#  * Added a function which takes a Sage matrix and solves
#    the exact cover problem for it.
#  * Recursive search is now iterative
#  * Removed callback functionality
#  * Revamped the class to be a pythonic generator; new usage:
#        for cover in DLXMatrix(ones, initialsolution):
#            blah(cover)
#  * DOCUMENTATION AND TESTS GALORE HOLYCRAP 100% COVERAGE!

# DLXMatrix class is used to store and solve an exact cover problem
# with help of Dancing Links [1] technique by Donald Knuth. The
# problem can be stated as an attempt to leave out rows of a 0/1
# matrix until remaining matrix has exactly one 1 in each column.
# Knuth proposes a fast solution technique based on clever trick with
# double linked list.
#
# This implementation will return row numbers of the solution instead
# column indexes (like in the Knuth's paper).
#
# [1] Donald E Knuth, Dancing links, preprint, available at
# http://www-cs-faculty.stanford.edu/~knuth/preprints.html

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

ROOTNODE = 0

# Node's attributes
# LEFT, RIGHT, UP and DOWN are identifiers of corresponding neighbors
# INDEX is the (row) index of the node (None in the header nodes)
# COLUMN/COUNT is column index on a regular node and count of nodes on
# a header node
LEFT = 0
RIGHT = 1
UP = 2
DOWN = 3
COLUMN = 4
INDEX = 5
COUNT = 5


class DLXMatrix:
    def __init__(self, ones, initialsolution=None):
        """
        Solve the Exact Cover problem by using the Dancing Links algorithm
        described by Knuth.

        Consider a matrix M with entries of 0 and 1, and compute a subset
        of the rows of this matrix which sum to the vector of all 1s.

        The dancing links algorithm works particularly well for sparse
        matrices, so the input is a list of lists of the form: (note the
        1-index!)::

          [
           [1, [i_11,i_12,...,i_1r]]
           ...
           [m,[i_m1,i_m2,...,i_ms]]
          ]

        where M[j][i_jk] = 1.

        The first example below corresponds to the matrix::

           1110
           1010
           0100
           0001

        which is exactly covered by::

           1110
           0001

        and

        ::

           1010
           0100
           0001

        EXAMPLES::

            sage: from sage.combinat.dlx import *
            sage: ones = [[1,[1,2,3]]]
            sage: ones+= [[2,[1,3]]]
            sage: ones+= [[3,[2]]]
            sage: ones+= [[4,[4]]]
            sage: DLXM = DLXMatrix(ones,[4])
            sage: for C in DLXM:
            ....:      print(C)
            [4, 1]
            [4, 2, 3]

        .. NOTE::

            The 0 entry is reserved internally for headers in the
            sparse representation, so rows and columns begin their
            indexing with 1.  Apologies for any heartache this
            causes. Blame the original author, or fix it yourself.
        """
        if initialsolution is None:
            initialsolution = []
        self._cursolution = []
        self._nodes = [[0, 0, None, None, None, None]]
        self._constructmatrix(ones, initialsolution)
        self._level = 0
        self._stack = [(None, None)]

    def __eq__(self, other):
        r"""
        Return ``True`` if every attribute of ``other`` matches the attribute
        of ``self``.

        INPUT:

        - ``other`` -- a DLX matrix

        EXAMPLES::

            sage: from sage.combinat.dlx import *
            sage: M = DLXMatrix([[1,[1]]])
            sage: M == loads(dumps(M))
            True
        """
        if not isinstance(other, DLXMatrix):
            return False
        return self.__dict__ == other.__dict__

    def __iter__(self):
        """
        Return ``self``.

        TESTS::

            sage: from sage.combinat.dlx import *
            sage: M = DLXMatrix([[1,[1]]])
            sage: M.__iter__() is M
            True
        """

        return self

    def _walknodes(self, firstnode, direction):
        """
        Generator for iterating over all nodes in given ``direction`` (not
        including ``firstnode``).

        TESTS::

            sage: from sage.combinat.dlx import *
            sage: ones = [[1,[1,2,3]]]
            sage: ones+= [[2,[1,3]]]
            sage: ones+= [[3,[2]]]
            sage: ones+= [[4,[4]]]
            sage: DLX = DLXMatrix(ones,[4])
            sage: count = 0
            sage: for c in DLX._walknodes(ROOTNODE,RIGHT):
            ....:     count += DLX._nodes[c][COUNT]
            ....:     for d in DLX._walknodes(c,DOWN):
            ....:         count -= 1
            sage: count
            0
        """
        nodetable = self._nodes
        n = nodetable[firstnode][direction]
        while n != firstnode:
            yield n
            n = nodetable[n][direction]

    def _constructmatrix(self, ones, initialsolution=None):
        """
        Construct the (sparse) DLX matrix based on list ``'ones'``.

        The first component in the list elements is row index (which
        will be returned by solve) and the second component is list of
        column indexes of ones in given row.

        'initialsolution' is list of row indexes that are required to be
        part of the solution. They will be removed from the matrix.

        .. NOTE::

            Rows and cols are 1-indexed ; the zero index is reserved for
            the root node and column heads.

        TESTS::

            sage: from sage.combinat.dlx import *
            sage: ones = [[1,[1,2,3]]]
            sage: ones+= [[2,[1,3]]]
            sage: ones+= [[3,[2]]]
            sage: ones+= [[4,[4]]]
            sage: DLX = DLXMatrix([[1,[1]]])
            sage: DLX._constructmatrix(ones,[4])
            sage: c = DLX._nodes[ROOTNODE][RIGHT]
            sage: fullcount = 0
            sage: while c != ROOTNODE:
            ....:     fullcount += DLX._nodes[c][COUNT]
            ....:     d = DLX._nodes[c][DOWN]
            ....:     while d != c:
            ....:         bad = DLX._nodes[DLX._nodes[d][DOWN]][UP] != d
            ....:         bad|= DLX._nodes[DLX._nodes[d][UP]][DOWN] != d
            ....:         bad|= DLX._nodes[DLX._nodes[d][LEFT]][RIGHT] != d
            ....:         bad|= DLX._nodes[DLX._nodes[d][RIGHT]][LEFT] != d
            ....:         if bad:
            ....:             raise RuntimeError("linked list inconsistent")
            ....:         d = DLX._nodes[d][DOWN]
            ....:     c = DLX._nodes[c][RIGHT]
            sage: fullcount
            6
        """
        if initialsolution is None:
            initialsolution = []
        self._cursolution = []
        # LEFT, RIGHT, UP, DOWN, COLUMN, INDEX/COUNT
        self._nodes = [[ROOTNODE, ROOTNODE, None, None, None, None]]

        # optimization: local variables are faster
        nodetable = self._nodes
        ones.sort()
        pruneNodes = []
        headers = [ROOTNODE]  # indexes of header nodes for faster access
        for r in ones:
            curRow = r[0]  # row index
            columns = r[1]  # column indexes
            if not (columns):
                continue
            columns.sort()

            # Do we need more headers?
            while len(headers) <= columns[-1]:
                lastheader = headers[-1]
                newind = len(nodetable)
                nodetable.append([lastheader, ROOTNODE, newind, newind, None, 0])
                nodetable[ROOTNODE][LEFT] = newind
                nodetable[lastheader][RIGHT] = newind
                headers.append(newind)

            # Add new nodes to indexes newind..newind+len(columns)-1
            # LEFT and RIGHT links can be calculated when created the
            # node, only UP, DOWN and COUNT have to be updated when
            # adding new nodes
            newind = len(nodetable)
            for i, c in enumerate(columns):
                h = headers[c]
                l = newind + ((i - 1) % len(columns))
                r = newind + ((i + 1) % len(columns))
                nodetable.append([l, r, nodetable[h][UP], h, h, curRow])
                nodetable[nodetable[h][UP]][DOWN] = newind + i
                nodetable[h][UP] = newind + i
                nodetable[h][COUNT] += 1

            if curRow in initialsolution:
                pruneNodes.append(newind)

        # Remove columns that are required to be in the solution
        for n in pruneNodes:
            self._cursolution += [nodetable[n][INDEX]]
            self._covercolumn(nodetable[n][COLUMN])
            for j in self._walknodes(n, RIGHT):
                self._covercolumn(nodetable[j][COLUMN])

    def _covercolumn(self, c):
        """
        Perform the column covering operation, as described by Knuth's
        pseudocode::

           cover(c):
                i <- D[c]
                while i != c:
                    j <- R[i]
                    while j != i
                        D[U[j]] <- D[j]
                        U[D[j]] <- U[j]
                        N[C[j]] <- N[C[j]] - 1
                        j <- R[j]
                    i <- D[i]

        This is undone by the uncover operation.

        TESTS::

            sage: from sage.combinat.dlx import *
            sage: M = DLXMatrix([[1,[1,3]],[2,[1,2]],[3,[2]]])
            sage: one = M._nodes[ROOTNODE][RIGHT]
            sage: M._covercolumn(one)
            sage: two = M._nodes[ROOTNODE][RIGHT]
            sage: three = M._nodes[two][RIGHT]
            sage: M._nodes[three][RIGHT] == ROOTNODE
            True
            sage: M._nodes[two][COUNT]
            1
            sage: M._nodes[three][COUNT]
            0
        """
        nodetable = self._nodes
        nodetable[nodetable[c][RIGHT]][LEFT] = nodetable[c][LEFT]
        nodetable[nodetable[c][LEFT]][RIGHT] = nodetable[c][RIGHT]
        for i in self._walknodes(c, DOWN):
            for j in self._walknodes(i, RIGHT):
                nodetable[nodetable[j][DOWN]][UP] = nodetable[j][UP]
                nodetable[nodetable[j][UP]][DOWN] = nodetable[j][DOWN]
                nodetable[nodetable[j][COLUMN]][COUNT] -= 1

    def _uncovercolumn(self, c):
        """
        Perform the column uncovering operation, as described by Knuth's
        pseudocode::

            uncover(c):
                i <- U[c]
                while i != c:
                    j <- L[i]
                    while j != i
                        U[j] <- U[D[j]]
                        D[j] <- D[U[j]]
                        N[C[j]] <- N[C[j]] + 1
                        j <- L[j]
                    i <- U[i]

        This undoes by the cover operation since everything is done in the
        reverse order.

        TESTS::

            sage: from sage.combinat.dlx import *
            sage: M = DLXMatrix([[1,[1,3]],[2,[1,2]],[3,[2]]])
            sage: one = M._nodes[ROOTNODE][RIGHT]
            sage: M._covercolumn(one)
            sage: two = M._nodes[ROOTNODE][RIGHT]
            sage: M._uncovercolumn(one)
            sage: M._nodes[two][LEFT] == one
            True
            sage: M._nodes[two][COUNT]
            2
        """
        nodetable = self._nodes
        for i in self._walknodes(c, UP):
            for j in self._walknodes(i, LEFT):
                nodetable[nodetable[j][DOWN]][UP] = j
                nodetable[nodetable[j][UP]][DOWN] = j
                nodetable[nodetable[j][COLUMN]][COUNT] += 1
        nodetable[nodetable[c][RIGHT]][LEFT] = c
        nodetable[nodetable[c][LEFT]][RIGHT] = c

    def __next__(self):
        """
        Search for the first solution we can find, and return it.

        Knuth describes the Dancing Links algorithm recursively, though
        actually implementing it as a recursive algorithm is permissible
        only for highly restricted problems. (for example, the original
        author implemented this for Sudoku, and it works beautifully
        there)

        What follows is an iterative description of DLX::

            stack <- [(NULL)]
            level <- 0
            while level >= 0:
                cur <- stack[level]
                if cur = NULL:
                    if R[h] = h:
                        level <- level - 1
                        yield solution
                    else:
                        cover(best_column)
                        stack[level] = best_column
                else if D[cur] != C[cur]:
                    if cur != C[cur]:
                        delete solution[level]
                        for j in L[cur], L[L[cur]], ... , while j != cur:
                            uncover(C[j])
                    cur <- D[cur]
                    solution[level] <- cur
                    stack[level] <- cur
                    for j in R[cur], R[R[cur]], ... , while j != cur:
                        cover(C[j])
                    level <- level + 1
                    stack[level] <- (NULL)
                else:
                    if C[cur] != cur:
                        delete solution[level]
                        for j in L[cur], L[L[cur]], ... , while j != cur:
                            uncover(C[j])
                    uncover(cur)
                    level <- level - 1

        TESTS::

            sage: from sage.combinat.dlx import *
            sage: M = DLXMatrix([[1,[1,2]],[2,[2,3]],[3,[1,3]]])
            sage: while 1:
            ....:     try:
            ....:         C = next(M)
            ....:     except StopIteration:
            ....:         print("StopIteration")
            ....:         break
            ....:     print(C)
            StopIteration
            sage: M = DLXMatrix([[1,[1,2]],[2,[2,3]],[3,[3]]])
            sage: for C in M:
            ....:       print(C)
            [1, 3]
            sage: M = DLXMatrix([[1,[1]],[2,[2,3]],[3,[2]],[4,[3]]])
            sage: for C in M:
            ....:       print(C)
            [1, 2]
            [1, 3, 4]
        """
        nodetable = self._nodes  # optimization: local variables are faster

        while self._level >= 0:
            c, r = self._stack[self._level]
            if c is None:
                if nodetable[ROOTNODE][RIGHT] == ROOTNODE:
                    self._level -= 1
                    return self._cursolution
                else:
                    c = nodetable[ROOTNODE][RIGHT]
                    maxcount = nodetable[nodetable[ROOTNODE][RIGHT]][COUNT]
                    for j in self._walknodes(ROOTNODE, RIGHT):
                        if nodetable[j][COUNT] < maxcount:
                            c = j
                            maxcount = nodetable[j][COUNT]
                    self._covercolumn(c)
                    self._stack[self._level] = (c, c)
            elif nodetable[r][DOWN] != c:
                if c != r:
                    self._cursolution = self._cursolution[:-1]
                    for j in self._walknodes(r, LEFT):
                        self._uncovercolumn(nodetable[j][COLUMN])
                r = nodetable[r][DOWN]
                self._cursolution += [nodetable[r][INDEX]]
                for j in self._walknodes(r, RIGHT):
                    self._covercolumn(nodetable[j][COLUMN])
                self._stack[self._level] = (c, r)
                self._level += 1
                if len(self._stack) == self._level:
                    self._stack.append((None, None))
                else:
                    self._stack[self._level] = (None, None)
            else:
                if c != r:
                    self._cursolution = self._cursolution[:-1]
                    for j in self._walknodes(r, LEFT):
                        self._uncovercolumn(nodetable[j][COLUMN])
                self._uncovercolumn(c)
                self._level -= 1

        raise StopIteration

    next = __next__


def AllExactCovers(M):
    """
    Use A. Ajanki's DLXMatrix class to solve the exact cover
    problem on the matrix M (treated as a dense binary matrix).

    EXAMPLES::

        sage: # needs sage.modules
        sage: M = Matrix([[1,1,0],[1,0,1],[0,1,1]])  # no exact covers
        sage: for cover in AllExactCovers(M):
        ....:     print(cover)
        sage: M = Matrix([[1,1,0],[1,0,1],[0,0,1],[0,1,0]]) # two exact covers
        sage: for cover in AllExactCovers(M):
        ....:     print(cover)
        [(1, 1, 0), (0, 0, 1)]
        [(1, 0, 1), (0, 1, 0)]
    """
    ones = []
    r = 1    # damn 1-indexing
    for R in M.rows():
        row = [i for i, Ri in enumerate(R, start=1) if Ri]
        ones.append([r, row])
        r += 1
    for s in DLXMatrix(ones):
        yield [M.row(i - 1) for i in s]  # damn 1-indexing


def OneExactCover(M):
    """
    Use A. Ajanki's DLXMatrix class to solve the exact cover
    problem on the matrix M (treated as a dense binary matrix).

    EXAMPLES::

        sage: M = Matrix([[1,1,0],[1,0,1],[0,1,1]])  # no exact covers                  # needs sage.modules
        sage: OneExactCover(M)                                                          # needs sage.modules

        sage: M = Matrix([[1,1,0],[1,0,1],[0,0,1],[0,1,0]])  # two exact covers         # needs sage.modules
        sage: OneExactCover(M)                                                          # needs sage.modules
        [(1, 1, 0), (0, 0, 1)]
    """
    for s in AllExactCovers(M):
        return s
