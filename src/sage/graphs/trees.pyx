r"""
Generation of trees

This is an implementation of the algorithm for generating trees with `n`
vertices (up to isomorphism) in constant time per tree described in
[WROM1986]_.

AUTHORS:

- Ryan Dingman (2009-04-16): initial version
"""

from libc.limits cimport INT_MAX
from cysignals.memory cimport check_allocarray, sig_free

# from networkx import MultiGraph

from sage.graphs.graph import Graph
from sage.graphs.base.sparse_graph cimport SparseGraph
from sage.graphs.base.sparse_graph cimport SparseGraphBackend

cdef class TreeIterator:
    r"""
    This class iterates over all trees with n vertices (up to isomorphism).

    EXAMPLES::

        sage: from sage.graphs.trees import TreeIterator
        sage: def check_trees(n):
        ....:     trees = []
        ....:     for t in TreeIterator(n):
        ....:         if not t.is_tree():
        ....:             return False
        ....:         if t.num_verts() != n:
        ....:             return False
        ....:         if t.num_edges() != n - 1:
        ....:             return False
        ....:         for tree in trees:
        ....:             if tree.is_isomorphic(t):
        ....:                 return False
        ....:         trees.append(t)
        ....:     return True
        sage: check_trees(10)
        True

    ::

        sage: from sage.graphs.trees import TreeIterator
        sage: count = 0
        sage: for t in TreeIterator(15):
        ....:     count += 1
        sage: count
        7741
    """

    def __init__(self, int vertices):
        r"""
        Initialize an iterator over all trees with `n` vertices.

        EXAMPLES::

            sage: from sage.graphs.trees import TreeIterator
            sage: t = TreeIterator(100) # indirect doctest
            sage: print(t)
            Iterator over all trees with 100 vertices
        """
        self.vertices = vertices
        self.l = NULL
        self.current_level_sequence = NULL
        self.first_time = 1

    def __dealloc__(self):
        r"""
        EXAMPLES::

            sage: from sage.graphs.trees import TreeIterator
            sage: t = TreeIterator(100)
            sage: t = None # indirect doctest
        """
        sig_free(self.l)
        sig_free(self.current_level_sequence)

    def __str__(self):
        r"""
        EXAMPLES::

            sage: from sage.graphs.trees import TreeIterator
            sage: t = TreeIterator(100)
            sage: print(t)  # indirect doctest
            Iterator over all trees with 100 vertices
        """
        return "Iterator over all trees with %s vertices" % (self.vertices)

    def __iter__(self):
        r"""
        Return an iterator over all the trees with `n` vertices.

        EXAMPLES::

            sage: from sage.graphs.trees import TreeIterator
            sage: t = TreeIterator(4)
            sage: list(iter(t))
            [Graph on 4 vertices, Graph on 4 vertices]
        """
        return self

    def __next__(self):
        r"""
        Return the next tree with `n` vertices.

        EXAMPLES::

            sage: from sage.graphs.trees import TreeIterator
            sage: T = TreeIterator(5)
            sage: [t for t in T] # indirect doctest
            [Graph on 5 vertices, Graph on 5 vertices, Graph on 5 vertices]


        TESTS:

        This used to be broken for trees with no vertices
        and was fixed in :issue:`13719` ::

            sage: from sage.graphs.trees import TreeIterator
            sage: T = TreeIterator(0)
            sage: [t for t in T] # indirect doctest
            [Graph on 0 vertices]
        """

        if not self.first_time and not self.q:
            raise StopIteration

        if self.first_time == 1:
            self.first_time = 0
            if self.vertices:
                self.l = <int *>check_allocarray(self.vertices, sizeof(int))
                self.current_level_sequence = <int *>check_allocarray(self.vertices, sizeof(int))

                self.generate_first_level_sequence()
            else:
                self.q = 0
        else:
            self.generate_next_level_sequence()

        cdef int i
        cdef int vertex1
        cdef int vertex2
        cdef object G

        G = Graph(self.vertices, sparse=True)
        cdef SparseGraph SG = (<SparseGraphBackend?> G._backend)._cg

        for i in range(2, self.vertices + 1):
            vertex1 = i - 1
            vertex2 = self.current_level_sequence[i - 1] - 1
            SG.add_arc_unsafe(vertex1, vertex2)

        return G

    cdef int generate_first_level_sequence(self) noexcept:
        r"""
        Generates the level sequence representing the first tree with `n` vertices
        """
        cdef int i
        cdef int k

        k = (self.vertices / 2) + 1

        if self.vertices == 4:
            self.p = 3
        else:
            self.p = self.vertices
        self.q = self.vertices - 1
        self.h1 = k
        self.h2 = self.vertices
        if self.vertices % 2:
            self.c = INT_MAX  # oo
        else:
            self.c = self.vertices + 1

        self.r = k

        for i in range(1, k + 1):
            self.l[i - 1] = i
        for i in range(k + 1, self.vertices + 1):
            self.l[i - 1] = i - k + 1
        for i in range(self.vertices):
            self.current_level_sequence[i] = i
        if self.vertices > 2:
            self.current_level_sequence[k] = 1
        if self.vertices <= 3:
            self.q = 0

        return 0

    cdef int generate_next_level_sequence(self) noexcept:
        r"""
        Generates the level sequence representing the next tree with `n` vertices
        """
        cdef int i
        cdef int fixit = 0

        cdef int needr = 0
        cdef int needc = 0
        cdef int needh2 = 0

        cdef int n = self.vertices
        cdef int p = self.p
        cdef int q = self.q
        cdef int h1 = self.h1
        cdef int h2 = self.h2
        cdef int c = self.c
        cdef int r = self.r
        cdef int *l = self.l
        cdef int *w = self.current_level_sequence

        if c == n + 1 or p == h2 and (l[h1 - 1] == l[h2 - 1] + 1 and n - h2 > r - h1 or l[h1 - 1] == l[h2 - 1] and n - h2 + 1 < r - h1):
            if l[r - 1] > 3:
                p = r
                q = w[r - 1]
                if h1 == r:
                    h1 = h1 - 1
                fixit = 1
            else:
                p = r
                r = r - 1
                q = 2

        if p <= h1:
            h1 = p - 1
        if p <= r:
            needr = 1
        elif p <= h2:
            needh2 = 1
        elif l[h2 - 1] == l[h1 - 1] - 1 and n - h2 == r - h1:
            if p <= c:
                needc = 1
        else:
            c = INT_MAX

        cdef int oldp = p
        cdef int delta = q - p
        cdef int oldlq = l[q - 1]
        cdef int oldwq = w[q - 1]
        p = INT_MAX

        for i in range(oldp, n + 1):
            l[i - 1] = l[i - 1 + delta]
            if l[i - 1] == 2:
                w[i - 1] = 1
            else:
                p = i
                if l[i - 1] == oldlq:
                    q = oldwq
                else:
                    q = w[i - 1 + delta] - delta
                w[i - 1] = q
            if needr == 1 and l[i - 1] == 2:
                needr = 0
                needh2 = 1
                r = i - 1
            if needh2 == 1 and l[i - 1] <= l[i - 2] and i > r + 1:
                needh2 = 0
                h2 = i - 1
                if l[h2 - 1] == l[h1 - 1] - 1 and n - h2 == r - h1:
                    needc = 1
                else:
                    c = INT_MAX
            if needc == 1:
                if l[i - 1] != l[h1 - h2 + i - 1] - 1:
                    needc = 0
                    c = i
                else:
                    c = i + 1

        if fixit == 1:
            r = n - h1 + 1
            for i in range(r + 1, n + 1):
                l[i - 1] = i - r + 1
                w[i - 1] = i - 1
            w[r] = 1
            h2 = n
            p = n
            q = p - 1
            c = INT_MAX
        else:
            if p == INT_MAX:
                if l[oldp - 2] != 2:
                    p = oldp - 1
                else:
                    p = oldp - 2
                q = w[p - 1]
            if needh2 == 1:
                h2 = n
                if l[h2 - 1] == l[h1 - 1] - 1 and h1 == r:
                    c = n + 1
                else:
                    c = INT_MAX

        self.p = p
        self.q = q
        self.h1 = h1
        self.h2 = h2
        self.c = c
        self.r = r
        self.l = l
        self.current_level_sequence = w

        return 0
