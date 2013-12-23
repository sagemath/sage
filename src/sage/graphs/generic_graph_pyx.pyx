"""
GenericGraph Cython functions

AUTHORS:
    -- Robert L. Miller   (2007-02-13): initial version
    -- Robert W. Bradshaw (2007-03-31): fast spring layout algorithms
    -- Nathann Cohen                  : exhaustive search
"""

#*****************************************************************************
#           Copyright (C) 2007 Robert L. Miller <rlmillster@gmail.com>
#                         2007 Robert W. Bradshaw <robertwb@math.washington.edu>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/interrupt.pxi"
include 'sage/ext/cdefs.pxi'
include 'sage/ext/stdsage.pxi'

# import from Python standard library
from sage.misc.prandom import random

# import from third-party library
from sage.graphs.base.sparse_graph cimport SparseGraph


cdef class GenericGraph_pyx(SageObject):
    pass

def spring_layout_fast_split(G, **options):
    """
    Graphs each component of G separately, placing them adjacent to
    each other. This is done because on a disconnected graph, the
    spring layout will push components further and further from each
    other without bound, resulting in very tight clumps for each
    component.

    NOTE:
        If the axis are scaled to fit the plot in a square, the
        horizontal distance may end up being "squished" due to
        the several adjacent components.

    EXAMPLES:

        sage: G = graphs.DodecahedralGraph()
        sage: for i in range(10): G.add_cycle(range(100*i, 100*i+3))
        sage: from sage.graphs.generic_graph_pyx import spring_layout_fast_split
        sage: spring_layout_fast_split(G)
        {0: [0.452..., 0.247...], ..., 502: [25.7..., 0.505...]}

    AUTHOR:
        Robert Bradshaw
    """
    Gs = G.connected_components_subgraphs()
    pos = {}
    left = 0
    buffer = 1/sqrt(len(G))
    for g in Gs:
        cur_pos = spring_layout_fast(g, **options)
        xmin = min([x[0] for x in cur_pos.values()])
        xmax = max([x[0] for x in cur_pos.values()])
        if len(g) > 1:
            buffer = (xmax - xmin)/sqrt(len(g))
        for v, loc in cur_pos.iteritems():
            loc[0] += left - xmin + buffer
            pos[v] = loc
        left += xmax - xmin + buffer
    return pos

def spring_layout_fast(G, iterations=50, int dim=2, vpos=None, bint rescale=True, bint height=False, by_component = False, **options):
    """
    Spring force model layout

    This function primarily acts as a wrapper around run_spring,
    converting to and from raw c types.

    This kind of speed cannot be achieved by naive Cythonification of the
    function alone, especially if we require a function call (let alone
    an object creation) every time we want to add a pair of doubles.

    INPUT:

     - ``by_component`` - a boolean

    EXAMPLES::

        sage: G = graphs.DodecahedralGraph()
        sage: for i in range(10): G.add_cycle(range(100*i, 100*i+3))
        sage: from sage.graphs.generic_graph_pyx import spring_layout_fast
        sage: spring_layout_fast(G)
        {0: [-0.0733..., 0.157...], ..., 502: [-0.551..., 0.682...]}

    With ``split=True``, each component of G is layed out separately,
    placing them adjacent to each other. This is done because on a
    disconnected graph, the spring layout will push components further
    and further from each other without bound, resulting in very tight
    clumps for each component.

    NOTE:

        If the axis are scaled to fit the plot in a square, the
        horizontal distance may end up being "squished" due to
        the several adjacent components.

        sage: G = graphs.DodecahedralGraph()
        sage: for i in range(10): G.add_cycle(range(100*i, 100*i+3))
        sage: from sage.graphs.generic_graph_pyx import spring_layout_fast
        sage: spring_layout_fast(G, by_component = True)
        {0: [2.12..., -0.321...], ..., 502: [26.0..., -0.812...]}
    """

    if by_component:
        return spring_layout_fast_split(G, iterations=iterations, dim = dim,
                                        vpos = vpos, rescale = rescale, height = height,
                                        **options)

    G = G.to_undirected()
    vlist = G.vertices() # this defines a consistent order

    cdef int i, j, x
    cdef int n = G.order()
    if n == 0:
        return {}

    cdef double* pos = <double*>sage_malloc(n * dim * sizeof(double))
    if pos is NULL:
            raise MemoryError, "error allocating scratch space for spring layout"

    # convert or create the starting positions as a flat list of doubles
    if vpos is None:  # set the initial positions randomly in 1x1 box
        for i from 0 <= i < n*dim:
            pos[i] = random()
    else:
        for i from 0 <= i < n:
            loc = vpos[vlist[i]]
            for x from 0 <= x <dim:
                pos[i*dim + x] = loc[x]


    # here we construct a lexicographically ordered list of all edges
    # where elist[2*i], elist[2*i+1] represents the i-th edge
    cdef int* elist = <int*>sage_malloc( (2 * len(G.edges()) + 2) * sizeof(int)  )
    if elist is NULL:
        sage_free(pos)
        raise MemoryError, "error allocating scratch space for spring layout"

    cdef int cur_edge = 0

    for i from 0 <= i < n:
        for j from i < j < n:
            if G.has_edge(vlist[i], vlist[j]):
                elist[cur_edge] = i
                elist[cur_edge+1] = j
                cur_edge += 2

    # finish the list with -1, -1 which never gets matched
    # but does get compared against when looking for the "next" edge
    elist[cur_edge] = -1
    elist[cur_edge+1] = -1

    run_spring(iterations, dim, pos, elist, n, height)

    # recenter
    cdef double* cen
    cdef double r, r2, max_r2 = 0
    if rescale:
        cen = <double *>sage_malloc(sizeof(double) * dim)
        if cen is NULL:
            sage_free(elist)
            sage_free(pos)
            raise MemoryError, "error allocating scratch space for spring layout"
        for x from 0 <= x < dim: cen[x] = 0
        for i from 0 <= i < n:
            for x from 0 <= x < dim:
                cen[x] += pos[i*dim + x]
        for x from 0 <= x < dim: cen[x] /= n
        for i from 0 <= i < n:
            r2 = 0
            for x from 0 <= x < dim:
                pos[i*dim + x] -= cen[x]
                r2 += pos[i*dim + x] * pos[i*dim + x]
            if r2 > max_r2:
                max_r2 = r2
        r = 1 if max_r2 == 0 else sqrt(max_r2)
        for i from 0 <= i < n:
            for x from 0 <= x < dim:
                pos[i*dim + x] /= r
        sage_free(cen)

    # put the data back into a position dictionary
    vpos = {}
    for i from 0 <= i < n:
        vpos[vlist[i]] = [pos[i*dim+x] for x from 0 <= x < dim]

    sage_free(pos)
    sage_free(elist)

    return vpos


cdef run_spring(int iterations, int dim, double* pos, int* edges, int n, bint height):
    """
    Find a locally optimal layout for this graph, according to the
    constraints that neighboring nodes want to be a fixed distance
    from each other, and non-neighboring nodes always repel.

    This is not a true physical model of mutually-repulsive particles
    with springs, rather it is more a model of such things traveling,
    without any inertia, through an (ever thickening) fluid.

    TODO: The inertial model could be incorporated (with F=ma)
    TODO: Are the hard-coded constants here optimal?

    INPUT:
        iterations -- number of steps to take
        dim        -- number of dimensions of freedom
        pos        -- already initialized initial positions
                      Each vertex is stored as [dim] consecutive doubles.
                      These doubles are then placed consecutively in the array.
                      For example, if dim=3, we would have
                      pos = [x_1, y_1, z_1, x_2, y_2, z_2, ... , x_n, y_n, z_n]
        edges      -- List of edges, sorted lexicographically by the first
                      (smallest) vertex, terminated by -1, -1.
                      The first two entries represent the first edge, and so on.
        n          -- number of vertices in the graph
        height     -- if True, do not update the last coordinate ever

    OUTPUT:
        Modifies contents of pos.

    AUTHOR:
        Robert Bradshaw
    """

    cdef int cur_iter, cur_edge
    cdef int i, j, x

    # k -- the equilibrium distance between two adjacent nodes
    cdef double t = 1, dt = t/(1e-20 + iterations), k = sqrt(1.0/n)
    cdef double square_dist, force, scale
    cdef double* disp_i
    cdef double* disp_j
    cdef double* delta

    cdef double* disp = <double*>sage_malloc((n+1) * dim * sizeof(double))
    if disp is NULL:
            raise MemoryError, "error allocating scratch space for spring layout"
    delta = &disp[n*dim]

    if height:
        update_dim = dim-1
    else:
        update_dim = dim

    sig_on()

    for cur_iter from 0 <= cur_iter < iterations:
      cur_edge = 1 # offset by one for fast checking against 2nd element first
      # zero out the disp vectors
      memset(disp, 0, n * dim * sizeof(double))
      for i from 0 <= i < n:
          disp_i = disp + (i*dim)
          for j from i < j < n:
              disp_j = disp + (j*dim)

              for x from 0 <= x < dim:
                  delta[x] = pos[i*dim+x] - pos[j*dim+x]

              square_dist = delta[0] * delta[0]
              for x from 1 <= x < dim:
                  square_dist += delta[x] * delta[x]

              if square_dist < 0.01:
                  square_dist = 0.01

              # they repel according to the (capped) inverse square law
              force = k*k/square_dist

              # and if they are neighbors, attract according Hooke's law
              if edges[cur_edge] == j and edges[cur_edge-1] == i:
                  force -= sqrt(square_dist)/k
                  cur_edge += 2

              # add this factor into each of the involved points
              for x from 0 <= x < dim:
                  disp_i[x] += delta[x] * force
                  disp_j[x] -= delta[x] * force

      # now update the positions
      for i from 0 <= i < n:
          disp_i = disp + (i*dim)

          square_dist = disp_i[0] * disp_i[0]
          for x from 1 <= x < dim:
              square_dist += disp_i[x] * disp_i[x]

          scale = t / (1 if square_dist < 0.01 else sqrt(square_dist))

          for x from 0 <= x < update_dim:
              pos[i*dim+x] += disp_i[x] * scale

      t -= dt

    sig_off()

    sage_free(disp)

def binary(n, length=None):
    """
    A quick python int to binary string conversion.

    EXAMPLE:
        sage: sage.graphs.generic_graph_pyx.binary(389)
        '110000101'
        sage: Integer(389).binary()
        '110000101'
        sage: sage.graphs.generic_graph_pyx.binary(2007)
        '11111010111'
    """
    cdef mpz_t i
    mpz_init(i)
    mpz_set_ui(i,n)
    cdef char* s=mpz_get_str(NULL, 2, i)
    t=str(s)
    sage_free(s)
    mpz_clear(i)
    return t

def R(x):
    """
    A helper function for the graph6 format. Described in [McK]

    EXAMPLE:
        sage: from sage.graphs.generic_graph_pyx import R
        sage: R('110111010110110010111000001100000001000000001')
        'vUqwK@?G'

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)
    """
    # pad on the right to make a multiple of 6
    x += '0' * ( (6 - (len(x) % 6)) % 6)

    # split into groups of 6, and convert numbers to decimal, adding 63
    six_bits = ''
    cdef int i
    for i from 0 <= i < len(x)/6:
        six_bits += chr( int( x[6*i:6*(i+1)], 2) + 63 )
    return six_bits

def N(n):
    """
    A helper function for the graph6 format. Described in [McK]

    EXAMPLE:
        sage: from sage.graphs.generic_graph_pyx import N
        sage: N(13)
        'L'
        sage: N(136)
        '~?AG'

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)
    """
    if n < 63:
        return chr(n + 63)
    else:
        # get 18-bit rep of n
        n = binary(n)
        n = '0'*(18-len(n)) + n
        return chr(126) + R(n)

def N_inverse(s):
    """
    A helper function for the graph6 format. Described in [McK]

    EXAMPLE:
        sage: from sage.graphs.generic_graph_pyx import N_inverse
        sage: N_inverse('~??~?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O??????G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?CC?A?CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@EG???O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAACd?_C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?aC_G@A??O??_?A?????O@Z?_@M????GQ@_G@?C?')
        (63, '?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O??????G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?CC?A?CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@EG???O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAACd?_C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?aC_G@A??O??_?A?????O@Z?_@M????GQ@_G@?C?')
        sage: N_inverse('_???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG????I?J??Q??O?_@@??@??????')
        (32, '???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG????I?J??Q??O?_@@??@??????')

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)
    """
    if s[0] == chr(126): # first four bytes are N
        a = binary(ord(s[1]) - 63).zfill(6)
        b = binary(ord(s[2]) - 63).zfill(6)
        c = binary(ord(s[3]) - 63).zfill(6)
        n = int(a + b + c,2)
        s = s[4:]
    else: # only first byte is N
        o = ord(s[0])
        if o > 126 or o < 63:
            raise RuntimeError("The string seems corrupt: valid characters are \n" + ''.join([chr(i) for i in xrange(63,127)]))
        n = o - 63
        s = s[1:]
    return n, s

def R_inverse(s, n):
    """
    A helper function for the graph6 format. Described in [McK]

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)

    EXAMPLE:
        sage: from sage.graphs.generic_graph_pyx import R_inverse
        sage: R_inverse('?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O??????G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?CC?A?CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@EG???O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAACd?_C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?aC_G@A??O??_?A?????O@Z?_@M????GQ@_G@?C?', 63)
        '0000000000000000000000000000001000000000010000000001000010000000000000000000110000000000000000010100000010000000000001000000000010000000000...10000000000000000000000000000000010000000001011011000000100000000001001110000000000000000000000000001000010010000001100000001000000001000000000100000000'
        sage: R_inverse('???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG????I?J??Q??O?_@@??@??????', 32)
        '0000000000000000000001000000000000010000100000100000001000000000000000100000000100000...010000000000000100010000001000000000000000000000000000001010000000001011000000000000010010000000000000010000000000100000000001000001000000000000000001000000000000000000000000000000000000'

    """
    l = []
    cdef int i
    for i from 0 <= i < len(s):
        o = ord(s[i])
        if o > 126 or o < 63:
            raise RuntimeError("The string seems corrupt: valid characters are \n" + ''.join([chr(i) for i in xrange(63,127)]))
        a = binary(o-63)
        l.append( '0'*(6-len(a)) + a )
    m = "".join(l)
    return m

def D_inverse(s, n):
    """
    A helper function for the dig6 format.

    EXAMPLE:
        sage: from sage.graphs.generic_graph_pyx import D_inverse
        sage: D_inverse('?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O??????G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?CC?A?CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@EG???O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAACd?_C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?aC_G@A??O??_?A?????O@Z?_@M????GQ@_G@?C?', 63)
        '0000000000000000000000000000001000000000010000000001000010000000000000000000110000000000000000010100000010000000000001000000000010000000000...10000000000000000000000000000000010000000001011011000000100000000001001110000000000000000000000000001000010010000001100000001000000001000000000100000000'
        sage: D_inverse('???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG????I?J??Q??O?_@@??@??????', 32)
        '0000000000000000000001000000000000010000100000100000001000000000000000100000000100000...010000000000000100010000001000000000000000000000000000001010000000001011000000000000010010000000000000010000000000100000000001000001000000000000000001000000000000000000000000000000000000'

    """
    l = []
    cdef int i
    for i from 0 <= i < len(s):
        o = ord(s[i])
        if o > 126 or o < 63:
            raise RuntimeError("The string seems corrupt: valid characters are \n" + ''.join([chr(i) for i in xrange(63,127)]))
        a = binary(o-63)
        l.append( '0'*(6-len(a)) + a )
    m = "".join(l)
    return m[:n*n]

# Exhaustive search in graphs

cdef class SubgraphSearch:
    r"""
    This class implements methods to exhaustively search for labelled
    copies of a graph `H` in a larger graph `G`.

    It is possible to look for induced subgraphs instead, and to
    iterate or count the number of their occurrences.

    ALGORITHM:

    The algorithm is a brute-force search.  Let `V(H) =
    \{h_1,\dots,h_k\}`.  It first tries to find in `G` a possible
    representant of `h_1`, then a representant of `h_2` compatible
    with `h_1`, then a representant of `h_3` compatible with the first
    two, etc.

    This way, most of the time we need to test far less than `k!
    \binom{|V(G)|}{k}` subsets, and hope this brute-force technique
    can sometimes be useful.
    """
    def __init__(self, G, H, induced = False):
        r"""
        Constructor

        This constructor only checks there is no inconsistency in the
        input : `G` and `H` are both graphs or both digraphs and that `H`
        has order at least 2.

        EXAMPLE::

            sage: g = graphs.PetersenGraph()
            sage: g.subgraph_search(graphs.CycleGraph(5))
            Subgraph of (Petersen graph): Graph on 5 vertices

        TESTS:

        Test proper initialization and deallocation, see :trac:`14067`.
        We intentionally only create the class without doing any
        computations with it::

            sage: from sage.graphs.generic_graph_pyx import SubgraphSearch
            sage: SubgraphSearch(Graph(5), Graph(1))
            Traceback (most recent call last):
            ...
            ValueError: Searched graph should have at least 2 vertices.
            sage: SubgraphSearch(Graph(5), Graph(2))
            <sage.graphs.generic_graph_pyx.SubgraphSearch ...>
        """
        if H.order() <= 1:
            raise ValueError("Searched graph should have at least 2 vertices.")

        if sum([G.is_directed(), H.is_directed()]) == 1:
            raise ValueError("One graph can not be directed while the other is not.")

        self._initialization()

    def __iter__(self):
        r"""
        Returns an iterator over all the labeleld subgraphs of `G`
        isomorphic to `H`.

        EXAMPLE:

        Iterating through all the `P_3` of `P_5`::

            sage: from sage.graphs.generic_graph_pyx import SubgraphSearch
            sage: g = graphs.PathGraph(5)
            sage: h = graphs.PathGraph(3)
            sage: S = SubgraphSearch(g, h)
            sage: for p in S:
            ...      print p
            [0, 1, 2]
            [1, 2, 3]
            [2, 1, 0]
            [2, 3, 4]
            [3, 2, 1]
            [4, 3, 2]
        """
        self._initialization()
        return self

    def cardinality(self):
        r"""
        Returns the number of labelled subgraphs of `G` isomorphic to
        `H`.

        .. NOTE::

           This method counts the subgraphs by enumerating them all !
           Hence it probably is not a good idea to count their number
           before enumerating them :-)

        EXAMPLE:

        Counting the number of labelled `P_3` in `P_5`::

            sage: from sage.graphs.generic_graph_pyx import SubgraphSearch
            sage: g = graphs.PathGraph(5)
            sage: h = graphs.PathGraph(3)
            sage: S = SubgraphSearch(g, h)
            sage: S.cardinality()
            6
        """
        if self.nh > self.ng:
            return 0

        self._initialization()
        cdef int i

        i=0
        for _ in self:
            i+=1

        return i

    def _initialization(self):
        r"""
        Initialization of the variables.

        Once the memory allocation is done in :meth:`__cinit__`,
        several variables need to be set to a default value. As this
        operation needs to be performed before any call to
        :meth:`__iter__` or to :meth:`cardinality`, it is cleaner to
        create a dedicated method.

        EXAMPLE:

        Finding two times the first occurrence through the
        re-initialization of the instance ::

            sage: from sage.graphs.generic_graph_pyx import SubgraphSearch
            sage: g = graphs.PathGraph(5)
            sage: h = graphs.PathGraph(3)
            sage: S = SubgraphSearch(g, h)
            sage: S.__next__()
            [0, 1, 2]
            sage: S._initialization()
            sage: S.__next__()
            [0, 1, 2]
        """

        memset(self.busy, 0, self.ng * sizeof(int))
        # 0 is the first vertex we use, so it is at first busy
        self.busy[0] = 1
        # stack -- list of the vertices which are part of the partial copy of H
        # in G.
        #
        # stack[i] -- the integer corresponding to the vertex of G representing
        # the i-th vertex of H.
        #
        # stack[i] = -1 means that i is not represented
        # ... yet!

        self.stack[0] = 0
        self.stack[1] = -1

        # Number of representants we have already found. Set to 1 as vertex 0
        # is already part of the partial copy of H in G.
        self.active = 1

    def __cinit__(self, G, H, induced = False):
        r"""
        Cython constructor

        This method initializes all the C values.

        EXAMPLE::

            sage: g = graphs.PetersenGraph()
            sage: g.subgraph_search(graphs.CycleGraph(5))
            Subgraph of (Petersen graph): Graph on 5 vertices
        """

        # Storing the number of vertices
        self.ng = G.order()
        self.nh = H.order()

        # Storing the list of vertices
        self.g_vertices = G.vertices()

        # Are the graphs directed (in __init__(), we check
        # whether both are of the same type)
        self.directed = G.is_directed()

        cdef int i, j, k

        self.tmp_array = <int *>sage_malloc(self.ng * sizeof(int))
        if self.tmp_array == NULL:
            raise MemoryError()

        # Should we look for induced subgraphs ?
        if induced:
            self.is_admissible = vectors_equal
        else:
            self.is_admissible = vectors_inferior

        # static copies of the two graphs for more efficient operations
        self.g = DenseGraph(self.ng)
        self.h = DenseGraph(self.nh)

        # copying the adjacency relations in both G and H
        i = 0
        for row in G.adjacency_matrix():
            j = 0
            for k in row:
                if k:
                    self.g.add_arc(i, j)
                j += 1
            i += 1
        i = 0
        for row in H.adjacency_matrix():
            j = 0
            for k in row:
                if k:
                    self.h.add_arc(i, j)
                j += 1
            i += 1

        # A vertex is said to be busy if it is already part of the partial copy
        # of H in G.
        self.busy = <int *>sage_malloc(self.ng * sizeof(int))
        self.stack = <int *>sage_malloc(self.nh * sizeof(int))

        # vertices is equal to range(nh), as an int *variable
        self.vertices = <int *>sage_malloc(self.nh * sizeof(int))
        for 0 <= i < self.nh:
            self.vertices[i] = i

        # line_h_out[i] represents the adjacency sequence of vertex i
        # in h relative to vertices 0, 1, ..., i-1
        self.line_h_out = <int **>sage_malloc(self.nh * sizeof(int *))
        for 0 <= i < self.nh:
            self.line_h_out[i] = <int *> sage_malloc(self.nh * sizeof(int *))
            if self.line_h_out[i] is NULL:
                raise MemoryError()
            self.h.adjacency_sequence_out(i, self.vertices, i, self.line_h_out[i])

        # Similarly in the opposite direction (only useful if the
        # graphs are directed)
        if self.directed:
            self.line_h_in = <int **>sage_malloc(self.nh * sizeof(int *))
            for 0 <= i < self.nh:
                self.line_h_in[i] = <int *> sage_malloc(self.nh * sizeof(int *))
                if self.line_h_in[i] is NULL:
                    raise MemoryError()

                self.h.adjacency_sequence_in(i, self.vertices, i, self.line_h_in[i])

    def __next__(self):
        r"""
        Returns the next isomorphic subgraph if any, and raises a
        ``StopIteration`` otherwise.

        EXAMPLE::

            sage: from sage.graphs.generic_graph_pyx import SubgraphSearch
            sage: g = graphs.PathGraph(5)
            sage: h = graphs.PathGraph(3)
            sage: S = SubgraphSearch(g, h)
            sage: S.__next__()
            [0, 1, 2]
        """
        sig_on()
        cdef bint is_admissible
        cdef int * tmp_array = self.tmp_array

        # as long as there is a non-void partial copy of H in G
        while self.active >= 0:
            # If we are here and found nothing yet, we try the next possible
            # vertex as a representant of the active i-th vertex of H.
            self.i = self.stack[self.active] + 1
            # Looking for a vertex that is not busy and compatible with the
            # partial copy we have of H.
            while self.i < self.ng:
                if self.busy[self.i]:
                    self.i += 1
                else:
                    # Testing whether the vertex we picked is a
                    # correct extension by checking the edges from the
                    # vertices already selected to self.i satisfy the
                    # constraints
                    self.g.adjacency_sequence_out(self.active, self.stack, self.i, tmp_array)
                    is_admissible = self.is_admissible(self.active, tmp_array, self.line_h_out[self.active])

                    # If G and H are digraphs, we also need to ensure
                    # the edges going in the opposite direction
                    # satisfy the constraints
                    if is_admissible and self.directed:
                        self.g.adjacency_sequence_in(self.active, self.stack, self.i, tmp_array)
                        is_admissible = is_admissible and self.is_admissible(self.active, tmp_array, self.line_h_in[self.active])

                    if is_admissible:
                        break
                    else:
                        self.i += 1

            # If we have found a good representant of H's i-th vertex in G
            if self.i < self.ng:

                # updating the last vertex of the stack
                if self.stack[self.active] != -1:
                    self.busy[self.stack[self.active]] = 0
                self.stack[self.active] = self.i

                # We have found our copy !!!
                if self.active == self.nh-1:
                    sig_off()
                    return [self.g_vertices[self.stack[l]] for l in xrange(self.nh)]

                # We are still missing several vertices ...
                else:
                    self.busy[self.stack[self.active]] = 1
                    self.active += 1

                    # we begin the search of the next vertex at 0
                    self.stack[self.active] = -1

            # If we found no representant for the i-th vertex, it
            # means that we cannot extend the current copy of H so we
            # update the status of stack[active] and prepare to change
            # the previous vertex.

            else:
                if self.stack[self.active] != -1:
                    self.busy[self.stack[self.active]] = 0
                self.stack[self.active] = -1
                self.active -= 1

        sig_off()
        raise StopIteration

    def __dealloc__(self):
        r"""
        Freeing the allocated memory.
        """

        # Free the memory
        sage_free(self.busy)
        sage_free(self.stack)
        sage_free(self.vertices)
        for 0 <= i < self.nh:
            sage_free(self.line_h_out[i])
        sage_free(self.line_h_out)

        if self.directed:
            for 0 <= i < self.nh:
                sage_free(self.line_h_in[i])
            sage_free(self.line_h_in)

        if self.tmp_array != NULL:
            sage_free(self.tmp_array)

cdef inline bint vectors_equal(int n, int *a, int *b):
    r"""
    Tests whether the two given vectors are equal. Two integer vectors
    `a = (a_1, a_2, \dots, a_n)` and `b = (b_1, b_2, \dots, b_n)` are equal
    iff `a_i = b_i` for all `i = 1, 2, \dots, n`. See the function
    ``_test_vectors_equal_inferior()`` for unit tests.

    INPUT:

    - ``n`` -- positive integer; length of the vectors.

    - ``a``, ``b`` -- two vectors of integers.

    OUTPUT:

    - ``True`` if ``a`` and ``b`` are the same vector; ``False`` otherwise.
    """
    cdef int i = 0
    for 0 <= i < n:
        if a[i] != b[i]:
            return False
    return True

cdef inline bint vectors_inferior(int n, int *a, int *b):
    r"""
    Tests whether the second vector of integers is inferior to the first. Let
    `u = (u_1, u_2, \dots, u_k)` and `v = (v_1, v_2, \dots, v_k)` be two
    integer vectors of equal length. Then `u` is said to be less than
    (or inferior to) `v` if `u_i \leq v_i` for all `i = 1, 2, \dots, k`. See
    the function ``_test_vectors_equal_inferior()`` for unit tests. Given two
    equal integer vectors `u` and `v`, `u` is inferior to `v` and vice versa.
    We could also define two vectors `a` and `b` to be equal if `a` is
    inferior to `b` and `b` is inferior to `a`.

    INPUT:

    - ``n`` -- positive integer; length of the vectors.

    - ``a``, ``b`` -- two vectors of integers.

    OUTPUT:

    - ``True`` if ``b`` is inferior to (or less than) ``a``; ``False``
      otherwise.
    """
    cdef int i = 0
    for 0 <= i < n:
        if a[i] < b[i]:
            return False
    return True

##############################
# Further tests. Unit tests for methods, functions, classes defined with cdef.
##############################

def _test_vectors_equal_inferior():
    """
    Unit testing the function ``vectors_equal()``. No output means that no
    errors were found in the random tests.

    TESTS::

        sage: from sage.graphs.generic_graph_pyx import _test_vectors_equal_inferior
        sage: _test_vectors_equal_inferior()
    """
    from sage.misc.prandom import randint
    n = randint(500, 10**3)
    cdef int *u = <int *>sage_malloc(n * sizeof(int))
    cdef int *v = <int *>sage_malloc(n * sizeof(int))
    cdef int i
    # equal vectors: u = v
    for 0 <= i < n:
        u[i] = randint(-10**6, 10**6)
        v[i] = u[i]
    try:
        assert vectors_equal(n, u, v)
        assert vectors_equal(n, v, u)
        # Since u and v are equal vectors, then u is inferior to v and v is
        # inferior to u. One could also define u and v as being equal if
        # u is inferior to v and vice versa.
        assert vectors_inferior(n, u, v)
        assert vectors_inferior(n, v, u)
    except AssertionError:
        sage_free(u)
        sage_free(v)
        raise AssertionError("Vectors u and v should be equal.")
    # Different vectors: u != v because we have u_j > v_j for some j. Thus,
    # u_i = v_i for 0 <= i < j and u_j > v_j. For j < k < n - 2, we could have:
    # (1) u_k = v_k,
    # (2) u_k < v_k, or
    # (3) u_k > v_k.
    # And finally, u_{n-1} < v_{n-1}.
    cdef int j = randint(1, n//2)
    cdef int k
    for 0 <= i < j:
        u[i] = randint(-10**6, 10**6)
        v[i] = u[i]
    u[j] = randint(-10**6, 10**6)
    v[j] = u[j] - randint(1, 10**6)
    for j < k < n:
        u[k] = randint(-10**6, 10**6)
        v[k] = randint(-10**6, 10**6)
    u[n - 1] = v[n - 1] - randint(1, 10**6)
    try:
        assert not vectors_equal(n, u, v)
        assert not vectors_equal(n, v, u)
        # u is not inferior to v because at least u_j > v_j
        assert u[j] > v[j]
        assert not vectors_inferior(n, v, u)
        # v is not inferior to u because at least v_{n-1} > u_{n-1}
        assert v[n - 1] > u[n - 1]
        assert not vectors_inferior(n, u, v)
    except AssertionError:
        sage_free(u)
        sage_free(v)
        raise AssertionError("".join([
                    "Vectors u and v should not be equal. ",
                    "u should not be inferior to v, and vice versa."]))
    # Different vectors: u != v because we have u_j < v_j for some j. Thus,
    # u_i = v_i for 0 <= i < j and u_j < v_j. For j < k < n - 2, we could have:
    # (1) u_k = v_k,
    # (2) u_k < v_k, or
    # (3) u_k > v_k.
    # And finally, u_{n-1} > v_{n-1}.
    j = randint(1, n//2)
    for 0 <= i < j:
        u[i] = randint(-10**6, 10**6)
        v[i] = u[i]
    u[j] = randint(-10**6, 10**6)
    v[j] = u[j] + randint(1, 10**6)
    for j < k < n:
        u[k] = randint(-10**6, 10**6)
        v[k] = randint(-10**6, 10**6)
    u[n - 1] = v[n - 1] + randint(1, 10**6)
    try:
        assert not vectors_equal(n, u, v)
        assert not vectors_equal(n, v, u)
        # u is not inferior to v because at least u_{n-1} > v_{n-1}
        assert u[n - 1] > v[n - 1]
        assert not vectors_inferior(n, v, u)
        # v is not inferior to u because at least u_j < v_j
        assert u[j] < v[j]
        assert not vectors_inferior(n, u, v)
    except AssertionError:
        sage_free(u)
        sage_free(v)
        raise AssertionError("".join([
                    "Vectors u and v should not be equal. ",
                    "u should not be inferior to v, and vice versa."]))
    # different vectors u != v
    # What's the probability of two random vectors being equal?
    for 0 <= i < n:
        u[i] = randint(-10**6, 10**6)
        v[i] = randint(-10**6, 10**6)
    try:
        assert not vectors_equal(n, u, v)
        assert not vectors_equal(n, v, u)
    except AssertionError:
        sage_free(u)
        sage_free(v)
        raise AssertionError("Vectors u and v should not be equal.")
    # u is inferior to v, but v is not inferior to u
    for 0 <= i < n:
        v[i] = randint(-10**6, 10**6)
        u[i] = randint(-10**6, 10**6)
        while u[i] > v[i]:
            u[i] = randint(-10**6, 10**6)
    try:
        assert not vectors_equal(n, u, v)
        assert not vectors_equal(n, v, u)
        assert vectors_inferior(n, v, u)
        assert not vectors_inferior(n, u, v)
    except AssertionError:
        raise AssertionError(
            "u should be inferior to v, but v is not inferior to u.")
    finally:
        sage_free(u)
        sage_free(v)

cpdef tuple find_hamiltonian( G, long max_iter=100000, long reset_bound=30000, long backtrack_bound=1000, find_path=False ):
    r"""
    Randomized backtracking for finding hamiltonian cycles and paths.

    ALGORITHM:

    A path ``P`` is maintained during the execution of the algorithm. Initially
    the path will contain an edge of the graph. Every 10 iterations the path
    is reversed. Every ``reset_bound`` iterations the path will be cleared
    and the procedure is restarted. Every ``backtrack_bound`` steps we discard
    the last five vertices and continue with the procedure. The total number
    of steps in the algorithm is controlled by ``max_iter``. If a hamiltonian
    cycle or hamiltonian path is found it is returned. If the number of steps reaches
    ``max_iter`` then a longest path is returned. See OUTPUT for more details.


    INPUT:

    - ``G`` - Graph.

    - ``max_iter`` - Maximum number of iterations.

    - ``reset_bound`` - Number of iterations before restarting the
       procedure.

    - ``backtrack_bound`` - Number of iterations to elapse before
       discarding the last 5 vertices of the path.

    - ``find_path`` - If set to ``True``, will search a hamiltonian
       path. If ``False``, will search for a hamiltonian
       cycle. Default value is ``False``.

    OUTPUT:

    A pair ``(B,P)``, where ``B`` is a Boolean and ``P`` is a list of vertices.

        * If ``B`` is ``True`` and ``find_path`` is ``False``, ``P``
          represents a hamiltonian cycle.

        * If ``B`` is ``True`` and ``find_path`` is ``True``, ``P``
          represents a hamiltonian path.

        * If ``B`` is false, then ``P`` represents the longest path
          found during the execution of the algorithm.

    .. WARNING::

        May loop endlessly when run on a graph with vertices of degree
        1.

    EXAMPLES:

    First we try the algorithm in the Dodecahedral graph, which is
    hamiltonian, so we are able to find a hamiltonian cycle and a
    hamiltonian path ::

        sage: from sage.graphs.generic_graph_pyx import find_hamiltonian as fh
        sage: G=graphs.DodecahedralGraph()
        sage: fh(G)
        (True, [9, 10, 0, 19, 3, 2, 1, 8, 7, 6, 5, 4, 17, 18, 11, 12, 16, 15, 14, 13])
        sage: fh(G,find_path=True)
        (True, [8, 9, 10, 11, 18, 17, 4, 3, 19, 0, 1, 2, 6, 7, 14, 13, 12, 16, 15, 5])

    Another test, now in the Moebius-Kantor graph which is also
    hamiltonian, as in our previous example, we are able to find a
    hamiltonian cycle and path ::

        sage: G=graphs.MoebiusKantorGraph()
        sage: fh(G)
        (True, [5, 4, 3, 2, 10, 15, 12, 9, 1, 0, 7, 6, 14, 11, 8, 13])
        sage: fh(G,find_path=True)
        (True, [4, 5, 6, 7, 15, 12, 9, 1, 0, 8, 13, 10, 2, 3, 11, 14])

    Now, we try the algorithm on a non hamiltonian graph, the Petersen
    graph.  This graph is known to be hypohamiltonian, so a
    hamiltonian path can be found ::

        sage: G=graphs.PetersenGraph()
        sage: fh(G)
        (False, [7, 9, 4, 3, 2, 1, 0, 5, 8, 6])
        sage: fh(G,find_path=True)
        (True, [3, 8, 6, 1, 2, 7, 9, 4, 0, 5])

    We now show the algorithm working on another known hypohamiltonian
    graph, the generalized Petersen graph with parameters 11 and 2 ::

        sage: G=graphs.GeneralizedPetersenGraph(11,2)
        sage: fh(G)
        (False, [13, 11, 0, 10, 9, 20, 18, 16, 14, 3, 2, 1, 12, 21, 19, 8, 7, 6, 17, 15, 4, 5])
        sage: fh(G,find_path=True)
        (True, [7, 18, 20, 9, 8, 19, 17, 6, 5, 16, 14, 3, 4, 15, 13, 11, 0, 10, 21, 12, 1, 2])

    Finally, an example on a graph which does not have a hamiltonian
    path ::

        sage: G=graphs.HyperStarGraph(5,2)
        sage: fh(G,find_path=False)
        (False, ['00011', '10001', '01001', '11000', '01010', '10010', '00110', '10100', '01100'])
        sage: fh(G,find_path=True)
        (False, ['00101', '10001', '01001', '11000', '01010', '10010', '00110', '10100', '01100'])
    """

    from sage.misc.prandom import randint
    cdef int n = G.order()
    cdef int m = G.num_edges()

    #Initialize the path.
    cdef int *path = <int *>sage_malloc(n * sizeof(int))
    memset(path, -1, n * sizeof(int))

    #Initialize the membership array
    cdef bint *member = <bint *>sage_malloc(n * sizeof(int))
    memset(member, 0, n * sizeof(int))

    # static copy of the graph for more efficient operations
    cdef SparseGraph g = SparseGraph(n)
    # copying the adjacency relations in G
    cdef int i
    cdef int j
    i = 0
    for row in G.adjacency_matrix():
        j = 0
        for k in row:
            if k:
                g.add_arc(i, j)
            j += 1
        i += 1
    # Cache copy of the vertices
    cdef list vertices = g.verts()

    # A list to store the available vertices at each step
    cdef list available_vertices=[]

    #We now work towards picking a random edge
    #  First we pick a random vertex u
    cdef int x = randint( 0, n-1 )
    cdef int u = vertices[x]
    #  Then we pick at random a neighbor of u
    x = randint( 0, len(g.out_neighbors( u ))-1 )
    cdef int v = g.out_neighbors( u )[x]
    # This will be the first edge in the path
    cdef int length=2
    path[ 0 ] = u
    path[ 1 ] = v
    member[ u ] = True
    member[ v ] = True

    #Initialize all the variables neccesary to start iterating
    cdef bint done = False
    cdef long counter = 0
    cdef long bigcount = 0
    cdef int longest = length

    #Initialize a path to contain the longest path
    cdef int *longest_path = <int *>sage_malloc(n * sizeof(int))
    memset(longest_path, -1, n * sizeof(int))
    i = 0
    for 0 <= i < length:
        longest_path[ i ] = path[ i ]

    #Initialize a temporary path for flipping
    cdef int *temp_path = <int *>sage_malloc(n * sizeof(int))
    memset(temp_path, -1, n * sizeof(int))

    cdef bint longer = False
    cdef bint good = True

    while not done:
        counter = counter + 1
        if counter%10 == 0:
            #Reverse the path

            i=0
            for 0<= i < length/2:
                t=path[ i ]
                path[ i ] = path[ length - i - 1]
                path[ length -i -1 ] = t

        if counter > reset_bound:
            bigcount = bigcount + 1
            counter = 1

            #Time to reset the procedure
            for 0 <= i < n:
                member[ i ]=False
            #  First we pick a random vertex u
            x = randint( 0, n-1 )
            u = vertices[x]
            #  Then we pick at random a neighbor of u
            degree = len(g.out_neighbors( u ))
            x = randint( 0, degree-1 )
            v = g.out_neighbors( u )[x]
            #  This will be the first edge in the path
            length=2
            path[ 0 ] = u
            path[ 1 ] = v
            member[ u ] = True
            member[ v ] = True

        if counter%backtrack_bound == 0:
            for 0 <= i < 5:
                member[ path[length - i - 1] ] = False
            length = length - 5
        longer = False

        available_vertices = []
        for u in g.out_neighbors( path[ length-1 ] ):
            if not member[ u ]:
                available_vertices.append( u )

        n_available=len( available_vertices )
        if  n_available > 0:
            longer = True
            x=randint( 0, n_available-1 )
            path[ length ] = available_vertices[ x ]
            length = length + 1
            member [ available_vertices[ x ] ] = True

        if not longer and length > longest:

            for 0 <= i < length:
                longest_path[ i ] = path[ i ]

            longest = length
        if not longer:

            memset(temp_path, -1, n * sizeof(int))
            degree = len(g.out_neighbors( path[ length-1 ] ))
            while True:
                x = randint( 0, degree-1 )
                u = g.out_neighbors(path[length - 1])[ x ]
                if u != path[length - 2]:
                    break

            flag = False
            i=0
            j=0
            for 0 <= i < length:
                if i > length-j-1:
                    break
                if flag:
                    t=path[ i ]
                    path[ i ] = path[ length - j - 1]
                    path[ length - j - 1 ] = t
                    j=j+1
                if path[ i ] == u:
                    flag = True
        if length == n:
            if find_path:
                done=True
            else:
                done = g.has_arc( path[n-1], path[0] )

        if bigcount*reset_bound > max_iter:
            verts=G.vertices()
            output=[ verts[ longest_path[i] ] for i from 0<= i < longest ]
            sage_free( member )
            sage_free( path )
            sage_free( longest_path )
            sage_free( temp_path )
            return (False, output)
    # #
    # # Output test
    # #

    # Test adjacencies
    for 0 <=i < n-1:
        u = path[i]
        v = path[i + 1]
        #Graph is simple, so both arcs are present
        if not g.has_arc( u, v ):
            good = False
            break
    if good == False:
        raise RuntimeError( 'Vertices %d and %d are consecutive in the cycle but are not ajacent.'%(u,v) )
    if not find_path and not g.has_arc( path[0], path[n-1] ):
        raise RuntimeError( 'Vertices %d and %d are not ajacent.'%(path[0],path[n-1]) )
    for 0 <= u < n:
        member[ u ]=False

    for 0 <= u < n:
        if member[ u ]:
            good = False
            break
        member[ u ] = True
    if good == False:
        raise RuntimeError( 'Vertex %d appears twice in the cycle.'%(u) )
    verts=G.vertices()
    output=[ verts[path[i]] for i from 0<= i < length ]
    sage_free( member )
    sage_free( path )
    sage_free( longest_path )
    sage_free( temp_path )

    return (True,output)

