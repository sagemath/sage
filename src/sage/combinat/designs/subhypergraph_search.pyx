r"""
Hypergraph isomorphic copy search

This module implements a code for the following problem:

    **INPUT:** two hypergraphs `H_1, H_2`

    **OUTPUT:** a copy of `H_2` in `H_1`

It is also possible to enumerate all such copies, and to require that such
copies be induced copies. More formally:

    A copy of `H_2` in `H_1` is an injection `f:V(H_2)\mapsto V(H_1)` such that
    for any set `S_2\in E(H_2)` we have `f(S_2)\in E(H_1)`.

    It is an *induced* copy if no other set of `E(H_1)` is contained in
    `f(V(H_2))`, i.e. `|E(H_2)|=\{S:S\in E(H_1)\text{ and }S\subseteq
    f(V(H_2))\}`.

The functions implemented here lists all such injections. In particular, the
number of copies of `H` in itself is equal to `|Aut(H)|`.

The feature is available through
:meth:`IncidenceStructure.isomorphic_substructures_iterator`.

Implementation
--------------

A hypergraph is stored as a list of edges, each of which is a "dense" bitset
over `|V(H_1)|` points. In particular, two sets of distinct cardinalities
require the same memory space. A hypergraph is a C struct with the following
fields:

* ``n``, ``m`` -- (``int``) number of points and edges

* ``limbs`` -- (``int``) number of 64-bits blocks per set

* ``set_space`` -- (``uint64_t *``) address of the memory used to store the
  sets

* ``sets`` -- (``uint64_t **``) ``sets[i]`` points toward the ``limbs``
  blocks encoding set `i`. Note also that ``sets[i][limbs]`` is equal to the
  cardinality of ``set[i]``, so that ``sets`` has length
  ``m*(limbs+1)*sizeof(uint64_t)``.

* ``names`` -- (``int *``) associates an integer 'name' to each of the ``n``
  points

The operations used on this data structure are:

* ``void permute(hypergraph * h, int n1, int n2)`` -- exchanges points `n1`
  and `n2` in the data structure. Note that their names are also exchanged
  so that we still know which is which.

* ``int induced_hypergraph(hypergraph * h1, int n, hypergraph * tmp1)`` --
  stores in ``tmp1`` the hypergraph induced by the first `n` points,
  i.e. all sets `S` such that `S\subseteq \{0,...,n-1\}`. The function
  returns the number of such sets.

* ``void trace_hypergraph64(hypergraph * h, int n, hypergraph * tmp)`` -- stores
  in ``tmp1`` the trace of `h` on the first `n` points, i.e. all sets of the
  form `S\cap \{0, \ldots, n-1\}`.

Algorithm
---------

We try all possible assignments of a representant `r_i\in H_1` for every `i\in
H_2`. When we have picked a representant for the first `n<` points
`\{0, \ldots, n-1\}\subsetneq V(H_2)`, we check that:

* The hypergraph induced by the (ordered) list `0, \ldots, n-1` in `H_2` is
  equal to the one induced by `r_0, \ldots, r_{n-1}` in `H_1`.

* If `S\subseteq \{0,...,n-1\}` is contained in `c` sets of size `k` in
  `H_2`, then `\{r_i:i\in S\}` is contained in `\geq c` sets of size `k` in
  `H_1`. This is done by comparing the trace of the hypergraphs while
  remembering the original size of each set.

As we very often need to build the hypergraph obtained by the trace of the first
`n` points (for all possible `n`), those hypergraphs are cached. The hypergraphs
induced by the same points are handled similarly.

Limitations
-----------

**Number of points** For efficiency reason the implementation assumes that `H_2`
has `\leq 64` points. Making this work for larger values means that calls to
``qsort`` have to be replaced by calls to ``qsort_r`` (i.e. to sort the edges
you need to know the number of limbs per edge) and that induces a big slowdown
for small cases (~50% when this code was implemented). Also, 64 points for `H_2`
is already very very big considering the problem at hand. Even `|V(H_1)|> 64`
seems too much.

**Vertex ordering** The order of vertices in `H_2` has a huge influence on the
performance of the algorithm. If no set of `H_2` contains more that one of the
first `k<n` points, then almost all partial assignments of representants are
possible for the first `k` points (though the degree of the vertices is taken
into account). For this reason it is best to pick an ordering such that the
first vertices are contained in as many sets as possible together. A heuristic
is implemented at
:meth:`~sage.combinat.designs.subhypergraph_search.SubHypergraphSearch.relabel_heuristic`.

AUTHORS:

- Nathann Cohen (November 2014, written in various airports between Nice and
  Chennai).

Methods
-------
"""
# ****************************************************************************
#       Copyright (C) 2014 Nathann Cohen <nathann.cohen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from libc.stdlib cimport qsort
from libc.stdint cimport uint64_t
from cysignals.memory cimport sig_malloc, sig_calloc, sig_free

ctypedef struct hypergraph:
    int n
    int m
    int limbs
    uint64_t ** sets
    uint64_t * set_space
    int * names

cdef inline int bs_get(uint64_t * bitset, int index) noexcept:
    r"""
    Return a bit of a bitset
    """
    return (bitset[index//64]>>(index%64))&1

cdef inline void bs_set(uint64_t * bitset, int index, int bit) noexcept:
    r"""
    Set a bit of a bitset.

    "bit" *MUST* be equal to either 0 or to 1. The code does not involve any
    "if".
    """
    bitset[index//64] &= ~((<uint64_t> 1)<<index%64)
    bitset[index//64] |= (<uint64_t> bit)<<index%64

cdef inline int bs_issubset64(uint64_t * b1, uint64_t b2, int limbs) noexcept:
    r"""
    Test whether bitset ``b1`` (on ``limbs`` blocks) is a subset of b2 (one block).

    It implies in particular that all last `limbs-1` blocks of ``b1`` are equal
    to zero.
    """
    # Checking that the whole field is zero can be done by computing b1[0] |
    # b1[1] | b[2] | ..., then calling "if" on that variable. May be faster, but
    # I have no real situation on which to test it.
    cdef int i
    for i in range(1,limbs):
        if b1[i]:
            return 0
    return (b1[0]&(~b2)) == 0

cdef void h_free(hypergraph h) noexcept:
    r"""
    Free the hypergraph
    """
    sig_free(h.names)
    sig_free(h.set_space)
    sig_free(h.sets)
    h.names = NULL
    h.set_space = NULL
    h.sets = NULL

cdef hypergraph h_init(int n, list H) noexcept:
    r"""
    Build a C hypergraph from a list `H` of sets on `\{0,...,n-1\}`.
    """
    cdef int x,i
    cdef hypergraph h
    h.n = n
    h.m = len(H)
    h.limbs = (n+63) // 64  # =ceil(n/64)
    h.names = <int *> sig_malloc(sizeof(int)*n)
    h.sets = <uint64_t **> sig_malloc(h.m*sizeof(uint64_t *))
    h.set_space = <uint64_t *> sig_calloc(h.m*(h.limbs+1),sizeof(uint64_t))

    # Consistency check
    for S in H:
        for x in S:
            if x < 0 or x >= n:
                h.n = -1

    if (h.names == NULL or
        h.sets == NULL or
        h.set_space == NULL or
        h.n == -1):
        h.n = -1
        return h

    for i in range(n):
        h.names[i] = i
    for i in range(h.m):
        h.sets[i] = h.set_space+i*(h.limbs+1)

    for i,S in enumerate(H):
        for x in S:
            bs_set(h.sets[i],x,1)
        h.sets[i][h.limbs] = len(S)

    return h

cdef inline void permute(hypergraph * h, int n1, int n2) noexcept:
    r"""
    Permutes two points of h inplace.

    This is only a data structure change, as h still represents the same
    hypergraph. In particular, the names of `n1` and `n2` are also exchanged.
    """
    if n1==n2:
        return

    h.names[n1],h.names[n2] = h.names[n2],h.names[n1]

    cdef int i,b1,b2
    for i in range(h.m):
        b1 = bs_get(h.sets[i],n1)
        b2 = bs_get(h.sets[i],n2)
        bs_set(h.sets[i],n1,b2)
        bs_set(h.sets[i],n2,b1)

cdef induced_hypergraph(hypergraph * h, int n, hypergraph * tmp):
    r"""
    Fills tmp with the hypergraph induced by points {0,...,n-1} in h.

    Assumes `n<=64`.
    """
    cdef int i,num_sets
    num_sets = 0
    cdef uint64_t current_set = ((<uint64_t> 1)<<(n+1))-1
    for i in range(h.m):
        if bs_issubset64(h.sets[i],current_set,h.limbs):
            tmp.sets[num_sets][0] = h.sets[i][0]
            num_sets += 1
    tmp.m = num_sets
    tmp.n = n
    tmp.limbs = 1


cdef void trace_hypergraph64(hypergraph * h, int n, hypergraph * tmp) noexcept:
    r"""
    Store in `tmp` the trace of the sets on `\{0, \ldots, n-1\}` in h1.

    Note that the size of the sets are kept as they are, i.e. the size of a set
    stored in tmp is what it was in h. This is useful information we use to cut
    the exploration.

    Assumes `n<=64`.
    """
    cdef int i
    cdef uint64_t current_set = ((<uint64_t> 1)<<(n+1))-1
    for i in range(h.m):
        tmp.sets[i][0] = h.sets[i][0]&current_set
        tmp.sets[i][1] = h.sets[i][h.limbs]

    tmp.limbs = 1

cdef int is_subhypergraph_admissible(hypergraph h1, hypergraph * h2_trace, int n, hypergraph tmp1) noexcept:
    r"""
    If there are `c` sets of size `k` containing `S\subseteq \{0,...,n-1\}` in
    `h2`, then there must be `>=c` sets of size `k` containing `S` in h1. This
    function checks that this property hold.

    ``h2_trace`` is expected to point toward the trace of ``h2`` on its first
    `n` points. Besides, its sets must be sorted with respect to
    ``cmp_128_bits``.

    Assumes `n<=64`.
    """
    trace_hypergraph64(&h1,n,&tmp1)
    qsort(tmp1.sets,h1.m,sizeof(uint64_t *),cmp_128_bits)

    cdef int i1,i2
    i1 = -1
    for i2 in range(h2_trace.m):
        i1 += 1
        while (i1<h1.m and
               (tmp1.sets[i1][0] < h2_trace.sets[i2][0] or
                tmp1.sets[i1][1] < h2_trace.sets[i2][1])):
            i1 += 1
        if (i1>=h1.m or
            (tmp1.sets[i1][0] > h2_trace.sets[i2][0] or
             tmp1.sets[i1][1] > h2_trace.sets[i2][1])):
            return 0

    return 1

cdef int cmp_128_bits(const void * a, const void * b) noexcept nogil:
    r"""
    Lexicographic order on 128-bits words
    """
    cdef uint64_t * p1 = (<uint64_t **> a)[0]
    cdef uint64_t * p2 = (<uint64_t **> b)[0]
    if p1[0] > p2[0]:
        return 1
    elif p1[0] == p2[0]:
        return 1 if p1[1] > p2[1] else -1
    else:
        return -1

cdef int is_induced_admissible64(hypergraph h1, hypergraph * h2_induced, int n, hypergraph tmp1) noexcept:
    r"""
    Test if the hypergraph induced in h1 by 0,...,n-1 is equal to the hypergraph
    induced in h2 by 0,...,n-1.

    ``h2_induced`` is expected to be a pointer toward the hypergraph induced by
    the first n points of ``h2``. Its sets must be sorted according to
    ``cmp_128_bits``.

    Assumes `n<=64`.
    """
    # *MUST* cache the info of tmp2 which only depends on n !!!!
    # *MUST* cache the sets' size
    induced_hypergraph(&h1,n,&tmp1)

    if tmp1.m!=h2_induced.m:
        return 0

    qsort(tmp1.sets,tmp1.m,sizeof(uint64_t *),cmp_128_bits)

    cdef int i
    for i in range(tmp1.m):
        if tmp1.sets[i][0] != h2_induced.sets[i][0]:
            return 0

    return 1

cdef class SubHypergraphSearch:

    cdef hypergraph h1,h2,tmp1,tmp2
    cdef list points1,points2
    cdef int induced
    cdef int * step
    cdef hypergraph * h2_traces
    cdef hypergraph * h2_induced

    def __cinit__(self, H1, H2, induced):
        r"""
        See the documentation's class.

        EXAMPLES::

            sage: from sage.combinat.designs.subhypergraph_search import SubHypergraphSearch
            sage: g1 = IncidenceStructure(graphs.PetersenGraph().edges(sort=True, labels=False))
            sage: g2 = IncidenceStructure(graphs.CycleGraph(5).edges(sort=True, labels=False))
            sage: S = SubHypergraphSearch(g1,g2,0)
            sage: sum(1 for _ in S)
            120
        """
        self.points1 = H1._points
        self.points2 = H2._points
        self.induced = induced
        cdef int n1 = H1.num_points()
        cdef int n2 = H2.num_points()

        if n2>64:
            raise RuntimeError("H2 has {}>64 points".format(n2))

        self.h1   = h_init(n1,H1._blocks)
        self.h2   = h_init(n2,H2._blocks)
        self.tmp1 = h_init(n1,H1._blocks) # No actual need to fill them,
        self.tmp2 = h_init(n2,H2._blocks) # only allocate the memory

        self.step = <int *> sig_malloc((n2+1)*sizeof(int))

        # all possible traces/induced subgraphs for h2
        #
        # (calloc sets all internal pointers to NULL)
        self.h2_traces  = <hypergraph *> sig_calloc(n2+1,sizeof(hypergraph))
        self.h2_induced = <hypergraph *> sig_calloc(n2+1,sizeof(hypergraph))

        if (self.h1.n   == -1 or
            self.h2.n   == -1 or
            self.tmp1.n == -1 or
            self.tmp2.n == -1 or
            self.h2_traces  == NULL or
            self.h2_induced == NULL or
            self.step       == NULL):
                raise MemoryError # also calls __dealloc__

        self.relabel_heuristic()

        # cache all n2+1 traces of h2, we will need them often.
        cdef int i
        for i in range(n2+1):
            self.h2_traces[i] = h_init(n2,H2._blocks)
            if self.h2_traces[i].n == -1:
                raise MemoryError
            trace_hypergraph64(&self.h2,i,&self.h2_traces[i])
            qsort(self.h2_traces[i].sets,self.h2.m,sizeof(uint64_t *),cmp_128_bits)

        # cache all n2+1 induced subhypergraphs of h2, we will need them often.
        for i in range(n2+1):
            self.h2_induced[i] = h_init(n2,H2._blocks)
            if self.h2_induced[i].n == -1:
                raise MemoryError
            induced_hypergraph(&self.h2,i,&self.h2_induced[i])
            qsort(self.h2_induced[i].sets,self.h2_induced[i].m,sizeof(uint64_t *),cmp_128_bits)

    def __dealloc__(self):
        r"""
        Free the ressources
        """
        # If self.h2 was not allocated, then self.h2.n = -1
        cdef int i
        for i in range(self.h2.n+1):
            h_free(self.h2_traces[i])
            h_free(self.h2_induced[i])
        h_free(self.h1)
        h_free(self.h2)
        h_free(self.tmp1)
        h_free(self.tmp2)
        sig_free(self.step)
        sig_free(self.h2_traces)
        sig_free(self.h2_induced)

    def relabel_heuristic(self):
        r"""
        Relabel `H_2` in order to make the algorithm faster.

        Objective: we try to pick an ordering `p_1,...,p_k` of the points of
        `H_2` that maximizes the number of sets involving the first points in
        the ordering. One way to formalize the problems indicates that it may be
        NP-Hard (generalizes the max clique problem for graphs) so we do not try
        to solve it exactly: we just need a sufficiently good heuristic.

        Assuming that the first points are `p_1,...,p_k`, we determine `p_{k+1}`
        as the point `x` such that the number of sets `S` with `x\in S` and
        `S\cap \{p_1,...,p_k\}\neq \emptyset` is maximal. In case of ties, we
        take a point with maximum degree.

        This function is called when an instance of :class:`SubHypergraphSearch`
        is created.

        EXAMPLES::

            sage: d = designs.projective_plane(3)                                       # needs sage.schemes
            sage: d.isomorphic_substructures_iterator(d).relabel_heuristic()            # needs sage.schemes
        """
        cdef hypergraph h2 = self.h2
        cdef int x,y,i
        cdef list degree = [0]*h2.n
        cdef list degree_with_ordered_points = [0]*h2.n
        cdef uint64_t current_set = 0

        # pre-compute the degree of each point
        for x in range(h2.n):
            y = 0
            for i in range(h2.m):
                y += bs_get(h2.sets[i],x)
            degree[x] = y

        # Compute the actual ordering
        for x in range(h2.n):
            y = 0
            for i in range(h2.m):
                if h2.sets[i][0] & current_set:
                    y += bs_get(h2.sets[i],x)
            degree_with_ordered_points[x] = y
            next_point = max((degree_with_ordered_points[i],degree[i],i)
                             for i in range(x,h2.n))[2]
            permute(&h2,x,next_point)
            current_set += 1<<x

    def __iter__(self):
        r"""
        Iterate over all copies of h2 in h1.

        EXAMPLES:

        How many distinct `C_5` in Petersen's graph ? ::

            sage: P = graphs.PetersenGraph()
            sage: C = graphs.CycleGraph(5)
            sage: IP = IncidenceStructure(P.edges(sort=True, labels=False))
            sage: IC = IncidenceStructure(C.edges(sort=True, labels=False))
            sage: sum(1 for _ in IP.isomorphic_substructures_iterator(IC))
            120
        """
        cdef hypergraph h1 = self.h1
        cdef hypergraph h2 = self.h2
        cdef hypergraph tmp1 = self.tmp1
        cdef int * step = self.step

        if h1.n<h2.n or h1.m<h2.m:
            return

        self.step[0] = -1

        # n is the current point of H2 whose representant must be decided
        cdef int n = 0

        while n>=0:
            step[n] += 1

            # If we went through all permutations beginning with
            # h1.names[0],...,h1.names[n-1].
            if n+step[n]==h1.n:
                step[n] = -1
                n -= 1
                if n>=0:
                    permute(&h1,n,n+step[n])
                continue

            permute(&h1,n,n+step[n])

            # Filter the assignments. If this 'if' is replaced by 'if True:', the
            # code enumerates all n1!/(n1-n2)! injections from [n2] into [n1].
            if ((not self.induced or is_induced_admissible64(h1,&self.h2_induced[n],n,tmp1)) and
                is_subhypergraph_admissible(h1,&self.h2_traces[n],n,tmp1)):
                n += 1
                step[n] = -1
                if n == h2.n:
                    yield {self.points2[h2.names[i]]:self.points1[h1.names[i]] for i in range(h2.n)}
                    n -= 1
                else:
                    continue

            permute(&h1,n,n+step[n])
