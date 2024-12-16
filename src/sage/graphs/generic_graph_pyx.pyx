"""
GenericGraph Cython functions

AUTHORS:

- Robert L. Miller   (2007-02-13): initial version
- Robert W. Bradshaw (2007-03-31): fast spring layout algorithms
- Nathann Cohen                  : exhaustive search
"""

# ****************************************************************************
#       Copyright (C) 2007 Robert L. Miller <rlmillster@gmail.com>
#                     2007 Robert W. Bradshaw <robertwb@math.washington.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport check_allocarray, check_calloc, sig_free
from cysignals.signals cimport sig_on, sig_off

import cython

from sage.data_structures.binary_matrix cimport *
from libc.math cimport sqrt, fabs
from libc.string cimport memset
from memory_allocator cimport MemoryAllocator

from sage.cpython.string cimport char_to_str
from sage.libs.gmp.mpz cimport *
from sage.misc.prandom import random
from sage.graphs.base.static_sparse_graph cimport short_digraph
from sage.graphs.base.static_sparse_graph cimport init_short_digraph
from sage.graphs.base.static_sparse_graph cimport init_reverse
from sage.graphs.base.static_sparse_graph cimport free_short_digraph
from sage.graphs.base.static_sparse_graph cimport out_degree, has_edge


cdef class GenericGraph_pyx(SageObject):
    pass


def layout_split(layout_function, G, **options):
    """
    Graph each component of ``G`` separately with ``layout_function``,
    placing them adjacent to each other.

    This is done because several layout methods need the input graph to
    be connected. For instance, on a disconnected graph, the spring
    layout will push components further and further from each other
    without bound, resulting in very tight clumps for each component.

    .. NOTE::

        If the axis are scaled to fit the plot in a square, the
        horizontal distance may end up being "squished" due to
        the several adjacent components.

    EXAMPLES::

        sage: G = graphs.DodecahedralGraph()
        sage: for i in range(10): G.add_cycle(list(range(100*i, 100*i+3)))
        sage: from sage.graphs.generic_graph_pyx import layout_split, spring_layout_fast
        sage: D = layout_split(spring_layout_fast, G); D  # random
        {0: [0.77..., 0.06...],
         ...
         902: [3.13..., 0.22...]}

    AUTHOR:

    Robert Bradshaw
    """
    from copy import copy
    Gs = G.connected_components_subgraphs()
    pos = {}
    left = 0
    buffer = 1/sqrt(len(G))

    on_embedding = options.get('on_embedding', None)
    forest_roots = options.get('forest_roots', None)
    try:
        forest_roots = list(forest_roots) if forest_roots else None
    except TypeError:
        raise TypeError('forest_roots should be an iterable of vertices')

    if forest_roots or on_embedding:
        options = copy(options)
        options.pop('forest_roots', None)
        options.pop('on_embedding', None)

    for g in Gs:
        if on_embedding:
            # Restrict ``on_embedding`` to ``g``
            embedding_g = {v: on_embedding[v] for v in g}
            cur_pos = layout_function(g, on_embedding=embedding_g, **options)
        elif forest_roots:
            # Find a root for ``g`` (if any)
            tree_root = next((v for v in forest_roots if v in g), None)
            cur_pos = layout_function(g, tree_root=tree_root, **options)
        else:
            cur_pos = layout_function(g, **options)

        xmin = min(x[0] for x in cur_pos.values())
        xmax = max(x[0] for x in cur_pos.values())
        if len(g) > 1:
            buffer = max(1, (xmax - xmin)/sqrt(len(g)))
        for v, loc in cur_pos.items():
            loc[0] += left - xmin + buffer
            pos[v] = loc
        left += xmax - xmin + buffer

    if options.get('set_embedding', None):
        embedding = dict()
        for g in Gs:
            embedding.update(g.get_embedding())
        G.set_embedding(embedding)
    return pos


def spring_layout_fast(G, iterations=50, int dim=2, vpos=None, bint rescale=True,
                       bint height=False, by_component=False, **options):
    """
    Spring force model layout.

    This function primarily acts as a wrapper around :func:`run_spring`,
    converting to and from raw C types.

    This kind of speed cannot be achieved by naive Cythonification of the
    function alone, especially if we require a function call (let alone
    an object creation) every time we want to add a pair of doubles.

    INPUT:

    - ``by_component`` -- boolean

    EXAMPLES::

        sage: G = graphs.DodecahedralGraph()
        sage: for i in range(10): G.add_cycle(list(range(100*i, 100*i+3)))
        sage: from sage.graphs.generic_graph_pyx import spring_layout_fast
        sage: pos = spring_layout_fast(G)
        sage: pos[0]  # random
        [0.00..., 0.03...]
        sage: sorted(pos.keys()) == sorted(G)
        True

    With ``split=True``, each component of G is laid out separately,
    placing them adjacent to each other. This is done because on a
    disconnected graph, the spring layout will push components further
    and further from each other without bound, resulting in very tight
    clumps for each component.

    If the axis are scaled to fit the plot in a square, the
    horizontal distance may end up being "squished" due to
    the several adjacent components. ::

        sage: G = graphs.DodecahedralGraph()
        sage: for i in range(10): G.add_cycle(list(range(100*i, 100*i+3)))
        sage: from sage.graphs.generic_graph_pyx import spring_layout_fast
        sage: pos = spring_layout_fast(G, by_component = True)
        sage: pos[0]  # random
        [2.21..., -0.00...]
        sage: len(pos) == G.order()
        True
    """
    if by_component:
        return layout_split(spring_layout_fast, G, iterations=iterations,
                            dim=dim, vpos=vpos, rescale=rescale,
                            height=height, **options)
    elif not G:
        return {}

    G = G.to_undirected()
    cdef list vlist = list(G)  # this defines a consistent order

    cdef int i, j, x
    cdef int n = G.order()

    cdef double* pos = NULL  # position of each vertex (for dim=2: x1,y1,x2,y2,...)
    cdef int* elist = NULL   # lexicographically ordered list of edges (u1,v1,u2,v2,...)
    cdef double* cen = NULL  # array of 'dim' doubles
    try:
        elist = <int*>check_allocarray(2*G.size() + 2, sizeof(int))
        pos = <double*>check_allocarray(n*dim, sizeof(double))
        cen = <double*>check_calloc(dim, sizeof(double))
    except MemoryError:
        sig_free(pos)
        sig_free(elist)
        sig_free(cen)
        raise

    # Initialize the starting positions
    if vpos is None:
        for i in range(n*dim):
            pos[i] = random()  # random in 1x1 box
    else:
        for i in range(n):
            loc = vpos[vlist[i]]
            for x in range(dim):
                pos[i*dim + x] = loc[x]

    # Lexicographically ordered list of edges
    cdef int cur_edge = 0

    for i in range(n):
        for j in range(i+1, n):
            if G.has_edge(vlist[i], vlist[j]):
                elist[cur_edge] = i
                elist[cur_edge+1] = j
                cur_edge += 2

    # finish the list with -1, -1 which never gets matched
    # but does get compared against when looking for the "next" edge
    elist[cur_edge] = -1
    elist[cur_edge+1] = -1

    if dim == 2:
        run_spring(<int> iterations, <D_TWO> NULL, <double*> pos, <int*>elist, <int> n, <int> G.size(), <bint> height)
    elif dim == 3:
        run_spring(<int> iterations, <D_THREE> NULL, <double*> pos, <int*>elist, <int> n, <int> G.size(), <bint> height)
    else:
        raise ValueError("'dim' must be equal to 2 or 3")

    # recenter
    cdef double r, r2, max_r2 = 0
    if rescale:
        for i in range(n):
            for x in range(dim):
                cen[x] += pos[i*dim + x]
        for x in range(dim):
            cen[x] /= n
        for i in range(n):
            r2 = 0
            for x in range(dim):
                pos[i*dim + x] -= cen[x]
                r2 += pos[i*dim + x] * pos[i*dim + x]
            if r2 > max_r2:
                max_r2 = r2
        r = 1 if max_r2 == 0 else sqrt(max_r2)
        for i in range(n):
            for x in range(dim):
                pos[i*dim + x] /= r

    # put the data back into a position dictionary
    vpos = {}
    for i in range(n):
        vpos[vlist[i]] = [pos[i*dim + x] for x in range(dim)]

    sig_free(pos)
    sig_free(elist)
    sig_free(cen)

    return vpos


@cython.cdivision(True)
cdef run_spring(int iterations, dimension_t _dim, double* pos, int* edges, int n, int m, bint height):
    r"""
    Find a locally optimal layout for this graph, according to the
    constraints that neighboring nodes want to be a fixed distance
    from each other, and non-neighboring nodes always repel.

    This is not a true physical model of mutually-repulsive particles
    with springs, rather it is more a model of such things traveling,
    without any inertia, through an (ever thickening) fluid.

    TODO: The inertial model could be incorporated (with F=ma)
    TODO: Are the hard-coded constants here optimal?

    INPUT:

    - ``iterations`` -- number of steps to take
    - ``_dim`` -- number of dimensions of freedom. Provide a value of type `D_TWO`
      for 2 dimensions, or type `D_THREE` for three dimensions. The actual
      value does not matter: only its type is important.
    - ``pos`` -- already initialized initial positions. Each vertex is stored as
      [dim] consecutive doubles. These doubles are then placed consecutively
      in the array. For example, if dim=3, we would have
      pos = [x_1, y_1, z_1, x_2, y_2, z_2, ... , x_n, y_n, z_n]
    - ``edges`` -- List of edges, sorted lexicographically by the first (smallest)
      vertex, terminated by -1, -1. The first two entries represent the first
      edge, and so on.
    - ``n`` -- number of vertices in the graph
    - ``height`` -- if ``True``, do not update the last coordinate ever

    OUTPUT: modifies contents of pos

    AUTHOR:

    Robert Bradshaw
    """
    cdef int dim
    cdef int cur_iter, cur_edge
    cdef int i, j, x

    if dimension_t is D_TWO:
        dim = 2
    else:
        dim = 3

    # k -- the equilibrium distance between two adjacent nodes
    cdef double t = 1, dt = t/(1e-20 + iterations), k = sqrt(1.0/n)
    cdef double square_dist, dist, force, scale
    cdef double* disp_i
    cdef double* disp_j
    cdef double delta[3]
    cdef double d_tmp
    cdef double xx, yy, zz

    cdef double* disp = <double*>check_allocarray(n, dim * sizeof(double))

    if height:
        update_dim = dim - 1
    else:
        update_dim = dim

    sig_on()

    for cur_iter in range(iterations):
        cur_edge = 1  # offset by one for fast checking against 2nd element first
        # zero out the disp vectors
        memset(disp, 0, n * dim * sizeof(double))
        for i in range(n):
            disp_i = disp + (i*dim)
            for j in range(i + 1, n):
                disp_j = disp + (j*dim)

                for x in range(dim):
                    delta[x] = pos[i*dim + x] - pos[j*dim + x]

                xx = delta[0] * delta[0]
                yy = delta[1] * delta[1]
                if dim == 2:
                    square_dist = xx+yy
                else:
                    zz = delta[2] * delta[2]
                    square_dist = xx+yy+zz

                if square_dist < 0.0001:
                    square_dist = 0.0001

                # they repel according to the (capped) inverse square law
                force = (k*k)/square_dist

                # and if they are neighbors, attract according Hooke's law
                if edges[cur_edge] == j and edges[cur_edge - 1] == i:
                    if dim == 2:
                        dist = sqrt_approx(delta[0], delta[1], xx, yy)
                    else:
                        dist = sqrt(square_dist)
                    force -= dist/k
                    cur_edge += 2

                # add this factor into each of the involved points
                for x in range(dim):
                    d_tmp = delta[x] * force
                    disp_i[x] += d_tmp
                    disp_j[x] -= d_tmp

        # now update the positions
        for i in range(n):
            disp_i = disp + (i*dim)

            square_dist = disp_i[0] * disp_i[0]
            for x in range(1, dim):
                square_dist += disp_i[x] * disp_i[x]

            if square_dist < 0.0001:
                scale = 1
            else:
                scale = t/sqrt(square_dist)

            for x in range(update_dim):
                pos[i*dim+x] += disp_i[x] * scale

        t -= dt

    sig_off()
    sig_free(disp)


@cython.cdivision(True)
cdef inline double sqrt_approx(double x, double y, double xx, double yy) noexcept:
    r"""
    Approximation of `\sqrt(x^2+y^2)`.

    Assuming that `x > y > 0`, it is a taylor expansion at `x^2`. To see how
    'bad' the approximation is::

        sage: def dist(x, y):
        ....:    x = abs(x)
        ....:    y = abs(y)
        ....:    return max(x,y) + min(x,y)**2/(2*max(x,y))

        sage: polar_plot([1, lambda x: dist(cos(x), sin(x))], (0, 2*math.pi))           # needs sage.plot sage.symbolic
        Graphics object consisting of 2 graphics primitives
    """
    if xx < yy:
        x, y = y, x
        xx, yy = yy, xx

    x = fabs(x)

    return x + yy/(2*x)


def int_to_binary_string(n):
    """
    A quick python int to binary string conversion.

    INPUT:

    - ``n`` -- integer

    EXAMPLES::

        sage: sage.graphs.generic_graph_pyx.int_to_binary_string(389)
        '110000101'
        sage: Integer(389).binary()
        '110000101'
        sage: sage.graphs.generic_graph_pyx.int_to_binary_string(2007)
        '11111010111'
    """
    cdef mpz_t i
    cdef char* s
    mpz_init(i)
    mpz_set_ui(i, n)
    s = mpz_get_str(NULL, 2, i)
    t = char_to_str(s)
    sig_free(s)
    mpz_clear(i)
    return t


def binary_string_to_graph6(x):
    r"""
    Transform a binary string into its graph6 representation.

    This helper function is named `R` in [McK2015]_.

    INPUT:

    - ``x`` -- a binary string

    EXAMPLES::

        sage: from sage.graphs.generic_graph_pyx import binary_string_to_graph6
        sage: binary_string_to_graph6('110111010110110010111000001100000001000000001')
        'vUqwK@?G'
    """
    # The length of x must be a multiple of 6. We extend it with 0s.
    x += '0' * ((6 - (len(x) % 6)) % 6)

    # Split into groups of 6, and convert numbers to decimal, adding 63
    six_bits = ''
    cdef int i
    for i in range(len(x)/6):
        six_bits += chr(int(x[6*i:6*(i + 1)], 2) + 63)
    return six_bits


def small_integer_to_graph6(n):
    r"""
    Encode a small integer (i.e. a number of vertices) as a graph6 string.

    This helper function is named `N` [McK2015]_.

    INPUT:

    - ``n`` -- integer

    EXAMPLES::

        sage: from sage.graphs.generic_graph_pyx import small_integer_to_graph6
        sage: small_integer_to_graph6(13)
        'L'
        sage: small_integer_to_graph6(136)
        '~?AG'
    """
    if n < 63:
        return chr(n + 63)
    # get 18-bit rep of n
    n = int_to_binary_string(n)
    n = '0'*(18 - len(n)) + n
    return chr(126) + binary_string_to_graph6(n)


def length_and_string_from_graph6(s):
    r"""
    Return a pair ``(length, graph6_string)`` from a graph6 string of unknown length.

    This helper function is the inverse of `N` from [McK2015]_.

    INPUT:

    - ``s`` -- a graph6 string describing a binary vector (and encoding its
      length)

    EXAMPLES::

        sage: from sage.graphs.generic_graph_pyx import length_and_string_from_graph6
        sage: g6 = '~??~?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O?'
        sage: g6 += '?????G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?C'
        sage: g6 += 'C?A?CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@E'
        sage: g6 += 'G???O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAA'
        sage: g6 += 'Cd?_C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?'
        sage: g6 += 'aC_G@A??O??_?A?????O@Z?_@M????GQ@_G@?C?'
        sage: length_and_string_from_graph6(g6)
        (63, '?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O??????G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?CC?A?CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@EG???O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAACd?_C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?aC_G@A??O??_?A?????O@Z?_@M????GQ@_G@?C?')
        sage: g6 = '_???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG?'
        sage: g6 += '???I?J??Q??O?_@@??@??????'
        sage: length_and_string_from_graph6(g6)
        (32, '???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG????I?J??Q??O?_@@??@??????')
    """
    if s[0] == chr(126):  # first four bytes are N
        a = int_to_binary_string(ord(s[1]) - 63).zfill(6)
        b = int_to_binary_string(ord(s[2]) - 63).zfill(6)
        c = int_to_binary_string(ord(s[3]) - 63).zfill(6)
        n = int(a + b + c, 2)
        s = s[4:]
    else:  # only first byte is N
        o = ord(s[0])
        if o > 126 or o < 63:
            raise RuntimeError("the string seems corrupt: valid characters are \n" + ''.join(chr(i) for i in range(63, 127)))
        n = o - 63
        s = s[1:]
    return n, s


def binary_string_from_graph6(s, n):
    r"""
    Decode a binary string from its graph6 representation.

    This helper function is the inverse of `R` from [McK2015]_.

    INPUT:

    - ``s`` -- a graph6 string

    - ``n`` -- the length of the binary string encoded by ``s``

    EXAMPLES::

        sage: from sage.graphs.generic_graph_pyx import binary_string_from_graph6
        sage: g6 = '?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O?????'
        sage: g6 += '?G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?CC?A?'
        sage: g6 += 'CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@EG???'
        sage: g6 += 'O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAACd?_'
        sage: g6 += 'C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?aC_G'
        sage: g6 += '@A??O??_?A?????O@Z?_@M????GQ@_G@?C?'
        sage: binary_string_from_graph6(g6, 63)
        '0000000000000000000000000000001000000000010000000001000010000000000000000000110000000000000000010100000010000000000001000000000010000000000...10000000000000000000000000000000010000000001011011000000100000000001001110000000000000000000000000001000010010000001100000001000000001000000000100000000'
        sage: g6 = '???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG??'
        sage: g6 += '??I?J??Q??O?_@@??@??????'
        sage: binary_string_from_graph6(g6, 32)
        '0000000000000000000001000000000000010000100000100000001000000000000000100000000100000...010000000000000100010000001000000000000000000000000000001010000000001011000000000000010010000000000000010000000000100000000001000001000000000000000001000000000000000000000000000000000000'
    """
    cdef list l = []
    cdef int i
    for i in range(len(s)):
        o = ord(s[i])
        if o > 126 or o < 63:
            raise RuntimeError("the string seems corrupt: valid characters are \n" + ''.join(chr(i) for i in range(63, 127)))
        a = int_to_binary_string(o - 63)
        l.append('0'*(6 - len(a)) + a)
    m = "".join(l)
    return m


def binary_string_from_dig6(s, n):
    """
    A helper function for the dig6 format.

    INPUT:

    - ``s`` -- a graph6 string

    - ``n`` -- the length of the binary string encoded by ``s``

    EXAMPLES::

        sage: from sage.graphs.generic_graph_pyx import binary_string_from_dig6
        sage: d6 = '?????_@?CG??B??@OG?C?G???GO??W@a???CO???OACC?OA?P@G??O?????'
        sage: d6 += '?G??C????c?G?CC?_?@???C_??_?C????PO?C_??AA?OOAHCA___?CC?A?'
        sage: d6 += 'CAOGO??????A??G?GR?C?_o`???g???A_C?OG??O?G_IA????_QO@EG???'
        sage: d6 += 'O??C?_?C@?G???@?_??AC?AO?a???O?????A?_Dw?H???__O@AAOAACd?_'
        sage: d6 += 'C??G?G@??GO?_???O@?_O??W??@P???AG??B?????G??GG???A??@?aC_G'
        sage: d6 += '@A??O??_?A?????O@Z?_@M????GQ@_G@?C?'
        sage: binary_string_from_dig6(d6, 63)
        '0000000000000000000000000000001000000000010000000001000010000000000000000000110000000000000000010100000010000000000001000000000010000000000...10000000000000000000000000000000010000000001011011000000100000000001001110000000000000000000000000001000010010000001100000001000000001000000000100000000'
        sage: d6 = '???C?@AA?_?A?O?C??S??O?q_?P?CHD??@?C?GC???C??GG?C_??O?COG??'
        sage: d6 += '??I?J??Q??O?_@@??@??????'
        sage: binary_string_from_dig6(d6, 32)
        '0000000000000000000001000000000000010000100000100000001000000000000000100000000100000...010000000000000100010000001000000000000000000000000000001010000000001011000000000000010010000000000000010000000000100000000001000001000000000000000001000000000000000000000000000000000000'
    """
    cdef list l = []
    cdef int i
    for i in range(len(s)):
        o = ord(s[i])
        if o > 126 or o < 63:
            raise RuntimeError("the string seems corrupt: valid characters are \n" + ''.join(chr(i) for i in range(63, 127)))
        a = int_to_binary_string(o - 63)
        l.append('0'*(6 - len(a)) + a)
    m = "".join(l)
    return m[:n*n]


# Exhaustive search in graphs

cdef class SubgraphSearch:
    r"""
    This class implements methods to exhaustively search for
    copies of a graph `H` in a larger graph `G`.

    It is possible to look for induced subgraphs instead, and to
    iterate or count the number of their occurrences.

    ALGORITHM:

    The algorithm is a brute-force search.  Let `V(H) =
    \{h_1,\dots,h_k\}`.  It first tries to find in `G` a possible
    representative of `h_1`, then a representative of `h_2` compatible
    with `h_1`, then a representative of `h_3` compatible with the first
    two, etc.

    This way, most of the time we need to test far less than `k!
    \binom{|V(G)|}{k}` subsets, and hope this brute-force technique
    can sometimes be useful.

    .. NOTE::

        This algorithm does not take vertex/edge labels into account.
    """
    def __init__(self, G, H, induced=False):
        r"""
        Constructor.

        This constructor only checks there is no inconsistency in the
        input : `G` and `H` are both graphs or both digraphs and that `H`
        has order at least 2.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: g.subgraph_search(graphs.CycleGraph(5))                               # needs sage.modules
            Subgraph of (Petersen graph): Graph on 5 vertices

        TESTS:

        Test proper initialization and deallocation, see :issue:`14067`.
        We intentionally only create the class without doing any
        computations with it::

            sage: from sage.graphs.generic_graph_pyx import SubgraphSearch
            sage: SubgraphSearch(Graph(5), Graph(1))                                    # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: searched graph should have at least 2 vertices
            sage: SubgraphSearch(Graph(5), Graph(2))                                    # needs sage.modules
            <sage.graphs.generic_graph_pyx.SubgraphSearch ...>
        """
        if H.order() <= 1:
            raise ValueError("searched graph should have at least 2 vertices")

        if G.is_directed() != H.is_directed():
            raise ValueError("one graph cannot be directed while the other is not")

        G._scream_if_not_simple(allow_loops=True)
        H._scream_if_not_simple(allow_loops=True)

        self._initialization()

    def __iter__(self):
        r"""
        Return an iterator over all the labeled subgraphs of `G`
        isomorphic to `H`.

        EXAMPLES:

        Iterating through all the `P_3` of `P_5`::

            sage: from sage.graphs.generic_graph_pyx import SubgraphSearch
            sage: g = graphs.PathGraph(5)
            sage: h = graphs.PathGraph(3)
            sage: S = SubgraphSearch(g, h)                                              # needs sage.modules
            sage: for p in S:                                                           # needs sage.modules
            ....:     print(p)
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
        Return the number of labelled subgraphs of `G` isomorphic to
        `H`.

        .. NOTE::

           This method counts the subgraphs by enumerating them all !
           Hence it probably is not a good idea to count their number
           before enumerating them :-)

        EXAMPLES:

        Counting the number of labelled `P_3` in `P_5`::

            sage: from sage.graphs.generic_graph_pyx import SubgraphSearch
            sage: g = graphs.PathGraph(5)
            sage: h = graphs.PathGraph(3)
            sage: S = SubgraphSearch(g, h)                                              # needs sage.modules
            sage: S.cardinality()                                                       # needs sage.modules
            6

        Check that the method is working even when vertices or edges are of
        incomparable types (see :issue:`35904`)::

            sage: from sage.graphs.generic_graph_pyx import SubgraphSearch
            sage: G = Graph()
            sage: G.add_cycle(['A', 1, 2, 3, ('a', 1)])
            sage: H = Graph()
            sage: H.add_path("xyz")
            sage: S = SubgraphSearch(G, H)                                              # needs sage.modules
            sage: S.cardinality()                                                       # needs sage.modules
            10
        """
        if self.nh > self.ng:
            return 0

        self._initialization()
        cdef int i

        i = 0
        for _ in self:
            i += 1

        from sage.rings.integer import Integer
        return Integer(i)

    def _initialization(self):
        r"""
        Initialization of the variables.

        Once the memory allocation is done in :meth:`__cinit__`,
        several variables need to be set to a default value. As this
        operation needs to be performed before any call to
        :meth:`__iter__` or to :meth:`cardinality`, it is cleaner to
        create a dedicated method.

        EXAMPLES:

        Finding two times the first occurrence through the
        re-initialization of the instance ::

            sage: from sage.graphs.generic_graph_pyx import SubgraphSearch
            sage: g = graphs.PathGraph(5)
            sage: h = graphs.PathGraph(3)
            sage: S = SubgraphSearch(g, h)                                              # needs sage.modules
            sage: S.__next__()                                                          # needs sage.modules
            [0, 1, 2]
            sage: S._initialization()                                                   # needs sage.modules
            sage: S.__next__()                                                          # needs sage.modules
            [0, 1, 2]

        TESTS:

        Check that :issue:`21828` is fixed::

            sage: Poset().is_incomparable_chain_free(1,1)   # indirect doctest          # needs sage.modules
            True
        """
        cdef int i

        if self.ng > 0:
            # 0 is the first vertex we use, so it is at first busy
            self.busy[0] = 1
            for i in range(1, self.ng):
                self.busy[i] = 0
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

        # Number of representatives we have already found. Set to 1 as vertex 0
        # is already part of the partial copy of H in G.
        self.active = 1

    def __cinit__(self, G, H, induced=False):
        r"""
        Cython constructor.

        This method initializes all the C values.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: g.subgraph_search(graphs.CycleGraph(5))                               # needs sage.modules
            Subgraph of (Petersen graph): Graph on 5 vertices
        """
        self.mem = MemoryAllocator()

        # Storing the number of vertices
        self.ng = G.order()
        self.nh = H.order()

        # Storing the list of vertices
        self.g_vertices = list(G)
        cdef list h_vertices = list(H)

        # Are the graphs directed (in __init__(), we check
        # whether both are of the same type)
        self.directed = G.is_directed()

        cdef int i, j

        # A vertex is said to be busy if it is already part of the partial copy
        # of H in G.
        self.busy = <int *> self.mem.allocarray(self.ng, sizeof(int))
        self.tmp_array = <int *> self.mem.allocarray(self.ng, sizeof(int))
        self.stack = <int *> self.mem.allocarray(self.nh, sizeof(int))
        self.vertices = <int *> self.mem.allocarray(self.nh, sizeof(int))
        self.line_h_out = <int **> self.mem.allocarray(self.nh, sizeof(int *))
        self.line_h_in = <int **> self.mem.allocarray(self.nh, sizeof(int *)) if self.directed else NULL

        self.line_h_out[0] = <int *> self.mem.allocarray(self.nh*self.nh,
                                                         sizeof(int))
        if self.directed:
            self.line_h_in[0] = <int *> self.mem.allocarray(self.nh*self.nh,
                                                            sizeof(int))

        # Should we look for induced subgraphs ?
        if induced:
            self.is_admissible = vectors_equal
        else:
            self.is_admissible = vectors_inferior

        # static copies of the two graphs for more efficient operations
        self.g = DenseGraph(self.ng)
        self.h = DenseGraph(self.nh)

        # copying the adjacency relations in both G and H
        cdef dict vertex_to_int = {v: i for i, v in enumerate(self.g_vertices)}
        cdef bint undirected = not G.is_directed()
        for u, v in G.edge_iterator(labels=False):
            i = vertex_to_int[u]
            j = vertex_to_int[v]
            self.g.add_arc(i, j)
            if undirected:
                self.g.add_arc(j, i)

        vertex_to_int = {v: i for i, v in enumerate(h_vertices)}
        for u, v in H.edge_iterator(labels=False):
            i = vertex_to_int[u]
            j = vertex_to_int[v]
            self.h.add_arc(i, j)
            if undirected:
                self.h.add_arc(j, i)

        # vertices is equal to range(nh), as an int *variable
        for i in range(self.nh):
            self.vertices[i] = i

        # line_h_out[i] represents the adjacency sequence of vertex i
        # in h relative to vertices 0, 1, ..., i-1
        for i in range(self.nh):
            self.line_h_out[i] = self.line_h_out[0] + i*self.nh
            self.h.adjacency_sequence_out(i, self.vertices, i, self.line_h_out[i])

        # Similarly in the opposite direction (only useful if the
        # graphs are directed)
        if self.directed:
            for i in range(self.nh):
                self.line_h_in[i] = self.line_h_in[0] + i*self.nh
                self.h.adjacency_sequence_in(i, self.vertices, i, self.line_h_in[i])

    def __next__(self):
        r"""
        Return the next isomorphic subgraph if any, and raises a
        ``StopIteration`` otherwise.

        EXAMPLES::

            sage: from sage.graphs.generic_graph_pyx import SubgraphSearch
            sage: g = graphs.PathGraph(5)
            sage: h = graphs.PathGraph(3)
            sage: S = SubgraphSearch(g, h)                                              # needs sage.modules
            sage: S.__next__()                                                          # needs sage.modules
            [0, 1, 2]
        """
        if not self.ng:
            raise StopIteration
        sig_on()
        cdef bint is_admissible
        cdef int * tmp_array = self.tmp_array

        # as long as there is a non-void partial copy of H in G
        while self.active >= 0:
            # If we are here and found nothing yet, we try the next possible
            # vertex as a representative of the active i-th vertex of H.
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

            # If we have found a good representative of H's i-th vertex in G
            if self.i < self.ng:

                # updating the last vertex of the stack
                if self.stack[self.active] != -1:
                    self.busy[self.stack[self.active]] = 0
                self.stack[self.active] = self.i

                # We have found our copy !!!
                if self.active == self.nh-1:
                    sig_off()
                    return [self.g_vertices[self.stack[l]]
                            for l in range(self.nh)]

                # We are still missing several vertices ...
                else:
                    self.busy[self.stack[self.active]] = 1
                    self.active += 1

                    # we begin the search of the next vertex at 0
                    self.stack[self.active] = -1

            # If we found no representative for the i-th vertex, it
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

cdef inline bint vectors_equal(int n, int *a, int *b) noexcept:
    r"""
    Test whether the two given vectors are equal. Two integer vectors
    `a = (a_1, a_2, \dots, a_n)` and `b = (b_1, b_2, \dots, b_n)` are equal
    iff `a_i = b_i` for all `i = 1, 2, \dots, n`. See the function
    ``_test_vectors_equal_inferior()`` for unit tests.

    INPUT:

    - ``n`` -- positive integer; length of the vectors

    - ``a``, ``b`` -- two vectors of integers

    OUTPUT: ``True`` if ``a`` and ``b`` are the same vector; ``False`` otherwise
    """
    cdef int i
    for i in range(n):
        if a[i] != b[i]:
            return False
    return True

cdef inline bint vectors_inferior(int n, int *a, int *b) noexcept:
    r"""
    Test whether the second vector of integers is inferior to the first. Let
    `u = (u_1, u_2, \dots, u_k)` and `v = (v_1, v_2, \dots, v_k)` be two
    integer vectors of equal length. Then `u` is said to be less than
    (or inferior to) `v` if `u_i \leq v_i` for all `i = 1, 2, \dots, k`. See
    the function ``_test_vectors_equal_inferior()`` for unit tests. Given two
    equal integer vectors `u` and `v`, `u` is inferior to `v` and vice versa.
    We could also define two vectors `a` and `b` to be equal if `a` is
    inferior to `b` and `b` is inferior to `a`.

    INPUT:

    - ``n`` -- positive integer; length of the vectors

    - ``a``, ``b`` -- two vectors of integers

    OUTPUT:

    - ``True`` if ``b`` is inferior to (or less than) ``a``; ``False``
      otherwise.
    """
    cdef int i
    for i in range(n):
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
    cdef int *u = <int *>check_allocarray(n, sizeof(int))
    cdef int *v = <int *>check_allocarray(n, sizeof(int))
    cdef int i
    # equal vectors: u = v
    for i in range(n):
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
        sig_free(u)
        sig_free(v)
        raise AssertionError("Vectors u and v should be equal.")
    # Different vectors: u != v because we have u_j > v_j for some j. Thus,
    # u_i = v_i for 0 <= i < j and u_j > v_j. For j < k < n - 2, we could have:
    # (1) u_k = v_k,
    # (2) u_k < v_k, or
    # (3) u_k > v_k.
    # And finally, u_{n-1} < v_{n-1}.
    cdef int j = randint(1, n//2)
    cdef int k
    for i in range(j):
        u[i] = randint(-10**6, 10**6)
        v[i] = u[i]
    u[j] = randint(-10**6, 10**6)
    v[j] = u[j] - randint(1, 10**6)
    for k in range(j + 1, n):
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
        sig_free(u)
        sig_free(v)
        raise AssertionError("Vectors u and v should not be equal. "
                             "u should not be inferior to v, and vice versa.")
    # Different vectors: u != v because we have u_j < v_j for some j. Thus,
    # u_i = v_i for 0 <= i < j and u_j < v_j. For j < k < n - 2, we could have:
    # (1) u_k = v_k,
    # (2) u_k < v_k, or
    # (3) u_k > v_k.
    # And finally, u_{n-1} > v_{n-1}.
    j = randint(1, n//2)
    for i in range(j):
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
        sig_free(u)
        sig_free(v)
        raise AssertionError("Vectors u and v should not be equal. "
                             "u should not be inferior to v, and vice versa.")
    # different vectors u != v
    # What's the probability of two random vectors being equal?
    for i in range(n):
        u[i] = randint(-10**6, 10**6)
        v[i] = randint(-10**6, 10**6)
    try:
        assert not vectors_equal(n, u, v)
        assert not vectors_equal(n, v, u)
    except AssertionError:
        sig_free(u)
        sig_free(v)
        raise AssertionError("Vectors u and v should not be equal.")
    # u is inferior to v, but v is not inferior to u
    for i in range(n):
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
        raise AssertionError("u should be inferior to v, but v is not inferior to u.")
    finally:
        sig_free(u)
        sig_free(v)


cpdef tuple find_hamiltonian(G, long max_iter=100000, long reset_bound=30000,
                             long backtrack_bound=1000, find_path=False):
    r"""
    Randomized backtracking for finding Hamiltonian cycles and paths.

    ALGORITHM:

    A path ``P`` is maintained during the execution of the algorithm.
    Initially the path will contain an edge of the graph. Every 10
    iterations the path is reversed. Every ``reset_bound`` iterations
    the path will be cleared and the procedure is restarted. Every
    ``backtrack_bound`` steps we discard the last five vertices and
    continue with the procedure. The total number of steps in the
    algorithm is controlled by ``max_iter``. If a Hamiltonian cycle or
    Hamiltonian path is found it is returned. If the number of steps
    reaches ``max_iter`` then a longest path is returned. See OUTPUT
    for more details.

    INPUT:

    - ``G`` -- graph

    - ``max_iter`` -- maximum number of iterations

    - ``reset_bound`` -- number of iterations before restarting the
      procedure

    - ``backtrack_bound`` -- number of iterations to elapse before
      discarding the last 5 vertices of the path

    - ``find_path`` -- boolean (default: ``False``); if set to ``True``, will
      search a Hamiltonian path. If ``False``, will search for a Hamiltonian
      cycle.

    OUTPUT:

    A pair ``(B, P)``, where ``B`` is a Boolean and ``P`` is a list
    of vertices.

    * If ``B`` is ``True`` and ``find_path`` is ``False``, ``P``
      represents a Hamiltonian cycle.

    * If ``B`` is ``True`` and ``find_path`` is ``True``, ``P``
      represents a Hamiltonian path.

    * If ``B`` is ``False``, then ``P`` represents the longest path
      found during the execution of the algorithm.

    .. WARNING::

        May loop endlessly when run on a graph with vertices of degree 1.

    EXAMPLES:

    For demonstration purposes we fix a random seed::

        sage: set_random_seed(0)

    First we try the algorithm in the Dodecahedral graph, which is
    Hamiltonian, so we are able to find a Hamiltonian cycle and a
    Hamiltonian path::

        sage: from sage.graphs.generic_graph_pyx import find_hamiltonian as fh
        sage: G=graphs.DodecahedralGraph()
        sage: fh(G)
        (True, [12, 11, 10, 9, 13, 14, 15, 5, 4, 3, 2, 6, 7, 8, 1, 0, 19, 18, 17, 16])
        sage: fh(G,find_path=True)
        (True, [10, 0, 19, 3, 4, 5, 15, 16, 17, 18, 11, 12, 13, 9, 8, 1, 2, 6, 7, 14])

    Another test, now in the Möbius-Kantor graph which is also
    Hamiltonian, as in our previous example, we are able to find a
    Hamiltonian cycle and path::

        sage: G=graphs.MoebiusKantorGraph()
        sage: fh(G)
        (True, [15, 10, 2, 3, 4, 5, 13, 8, 11, 14, 6, 7, 0, 1, 9, 12])
        sage: fh(G,find_path=True)
        (True, [10, 15, 7, 6, 5, 4, 12, 9, 14, 11, 3, 2, 1, 0, 8, 13])

    Now, we try the algorithm on a non Hamiltonian graph, the Petersen
    graph.  This graph is known to be hypohamiltonian, so a
    Hamiltonian path can be found::

        sage: G=graphs.PetersenGraph()
        sage: fh(G)
        (False, [9, 4, 0, 1, 6, 8, 5, 7, 2, 3])
        sage: fh(G,find_path=True)
        (True, [7, 2, 1, 0, 5, 8, 6, 9, 4, 3])

    We now show the algorithm working on another known hypohamiltonian
    graph, the generalized Petersen graph with parameters 11 and 2::

        sage: G=graphs.GeneralizedPetersenGraph(11,2)
        sage: fh(G)
        (False, [7, 8, 9, 10, 0, 1, 2, 3, 14, 12, 21, 19, 17, 6, 5, 4, 15, 13, 11, 20, 18, 16])
        sage: fh(G,find_path=True)
        (True, [2, 1, 12, 21, 10, 0, 11, 13, 15, 17, 19, 8, 7, 6, 5, 4, 3, 14, 16, 18, 20, 9])

    Finally, an example on a graph which does not have a Hamiltonian
    path::

        sage: G = graphs.HyperStarGraph(5, 2)
        sage: G.order()
        10
        sage: b, P = fh(G,find_path=False)
        sage: b, len(P)
        (False, 9)
        sage: b, P = fh(G,find_path=True)
        sage: b, len(P)
        (False, 9)

    The method can also be used for directed graphs::

        sage: G = DiGraph([(0, 1), (1, 2), (2, 3)])
        sage: fh(G)
        (False, [0, 1, 2, 3])
        sage: G = G.reverse()
        sage: fh(G)
        (False, [3, 2, 1, 0])
        sage: G = DiGraph()
        sage: G.add_cycle([0, 1, 2, 3, 4, 5])
        sage: b, P = fh(G)
        sage: b, len(P)
        (True, 6)

    TESTS:

    :issue:`10206` -- Hamiltonian cycle in small (di)graphs::

        sage: for n in range(3):
        ....:     for G in graphs(n):
        ....:         print('order {} and size {}: {}'.format(G.order(),G.size(),fh(G, find_path=False)))
        order 0 and size 0: (False, [])
        order 1 and size 0: (False, [0])
        order 2 and size 0: (False, [0])
        order 2 and size 1: (False, [0, 1])
        sage: for n in range(3):
        ....:     for G in digraphs(n):
        ....:         print('order {} and size {}: {}'.format(G.order(),G.size(),fh(G, find_path=False)))
        order 0 and size 0: (False, [])
        order 1 and size 0: (False, [0])
        order 2 and size 0: (False, [0])
        order 2 and size 1: (False, [0, 1])
        order 2 and size 2: (False, [0, 1])

    :issue:`10206` -- Hamiltonian path in small (di)graphs::

        sage: for n in range(3):
        ....:     for G in graphs(n):
        ....:         print('order {} and size {}: {}'.format(G.order(),G.size(),fh(G, find_path=True)))
        order 0 and size 0: (False, [])
        order 1 and size 0: (False, [0])
        order 2 and size 0: (False, [0])
        order 2 and size 1: (True, [0, 1])
        sage: for n in range(3):
        ....:     for G in digraphs(n):
        ....:         print('order {} and size {}: {}'.format(G.order(),G.size(),fh(G, find_path=True)))
        order 0 and size 0: (False, [])
        order 1 and size 0: (False, [0])
        order 2 and size 0: (False, [0])
        order 2 and size 1: (True, [0, 1])
        order 2 and size 2: (True, [0, 1])

    :issue:`10206` -- disconnected graphs::

        sage: G = graphs.CompleteGraph(4) + Graph(1)
        sage: fh(G, find_path=False)
        (False, [0, 1, 2, 3])
        sage: fh(G, find_path=True)
        (False, [0, 1, 2, 3])

    Check that the method is robust to incomparable vertices::

        sage: G = Graph([(1, 'a'), ('a', 2), (2, 3), (3, 1)])
        sage: b, C = fh(G, find_path=False)
        sage: b, len(C)
        (True, 4)
    """
    G._scream_if_not_simple()

    from sage.misc.prandom import randint
    cdef int n = G.order()

    # Easy cases
    if n < 2:
        return False, list(G)

    # To clean the output when find_path is None or a number
    find_path = (find_path > 0)

    if G.is_clique(induced=False):
        # We have a hamiltonian path since n >= 2, but we have a hamiltonian
        # cycle only if n >= 3
        return find_path or n >= 3, list(G)

    cdef list best_path, p
    if not G.is_connected():
        # The (Di)Graph has no hamiltonian path or cycle. We search for the
        # longest path in its connected components.
        best_path = []
        for H in G.connected_components_subgraphs():
            if H.order() <= len(best_path):
                continue
            _, p = find_hamiltonian(H, max_iter=max_iter, reset_bound=reset_bound,
                                    backtrack_bound=backtrack_bound, find_path=True)
            if len(p) > len(best_path):
                best_path = p
        return False, best_path

    # Misc variables used below
    cdef int i, j
    cdef bint directed = G.is_directed()

    # Initialize the path.
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *path = <int *>mem.allocarray(n, sizeof(int))

    # Initialize the membership array
    cdef bint *member = <bint *>mem.allocarray(n, sizeof(int))
    memset(member, 0, n * sizeof(int))

    # static copy of the graph for more efficient operations
    cdef list int_to_vertex = list(G)
    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_vertex)
    cdef short_digraph rev_sd
    cdef bint reverse = False
    if directed:
        init_reverse(rev_sd, sd)

    # A list to store the available vertices at each step
    cdef list available_vertices = []

    # We now work towards picking a random edge
    # First we pick a random vertex u of (out-)degree at least one
    cdef int u = randint(0, n - 1)
    while not out_degree(sd, u):
        u = randint(0, n - 1)
    # Then we pick at random a neighbor of u
    cdef int x = randint(0, out_degree(sd, u) - 1)
    cdef int v = sd.neighbors[u][x]
    # This will be the first edge in the path
    cdef int length = 2
    path[0] = u
    path[1] = v
    member[u] = True
    member[v] = True

    # Initialize all the variables necessary to start iterating
    cdef bint done = False
    cdef long counter = 0
    cdef long bigcount = 0
    cdef int longest = length

    # Initialize a path to contain the longest path
    cdef int *longest_path = <int *>mem.allocarray(n, sizeof(int))
    for i in range(length):
        longest_path[i] = path[i]

    cdef bint longer = False
    cdef bint longest_reversed = False
    cdef bint flag

    while not done:
        counter = counter + 1
        if counter % 10 == 0:
            # Reverse the path
            for i in range(length//2):
                t = path[i]
                path[i] = path[length - i - 1]
                path[length - i - 1] = t

            if directed:
                # We now work on the reverse graph
                reverse = not reverse

        if counter > reset_bound:
            bigcount = bigcount + 1
            counter = 1

            # Time to reset the procedure
            memset(member, 0, n * sizeof(int))
            if directed and reverse:
                # We restore the original orientation
                reverse = False

            # First we pick a random vertex u of (out-)degree at least one
            u = randint(0, n - 1)
            while not out_degree(sd, u):
                u = randint(0, n - 1)
            # Then we pick at random a neighbor of u
            x = randint(0, out_degree(sd, u) - 1)
            v = sd.neighbors[u][x]
            # This will be the first edge in the path
            length = 2
            path[0] = u
            path[1] = v
            member[u] = True
            member[v] = True

        if length > 5 and counter % backtrack_bound == 0:
            for i in range(5):
                member[path[length - i - 1]] = False
            length = length - 5
        longer = False

        # We search for a possible extension of the path
        available_vertices = []
        u = path[length - 1]
        if directed and reverse:
            for i in range(out_degree(rev_sd, u)):
                v = rev_sd.neighbors[u][i]
                if not member[v]:
                    available_vertices.append(v)
        else:
            for i in range(out_degree(sd, u)):
                v = sd.neighbors[u][i]
                if not member[v]:
                    available_vertices.append(v)

        if available_vertices:
            longer = True
            x = randint(0, len(available_vertices) - 1)
            v = available_vertices[x]
            path[length] = v
            length = length + 1
            member[v] = True

        if not longer and length > longest:
            # Store the current best solution
            for i in range(length):
                longest_path[i] = path[i]

            longest = length
            longest_reversed = reverse

        if not directed and not longer and out_degree(sd, path[length - 1]) > 1:
            # We revert a cycle to change the extremity of the path
            degree = out_degree(sd, path[length - 1])
            while True:
                x = randint(0, degree - 1)
                u = sd.neighbors[path[length - 1]][x]
                if u != path[length - 2]:
                    break

            flag = False
            j = 0
            for i in range(length):
                if i > length - j - 1:
                    break
                if flag:
                    t = path[i]
                    path[i] = path[length - j - 1]
                    path[length - j - 1] = t
                    j += 1
                if path[i] == u:
                    flag = True

        if length == n:
            if find_path:
                done = True
            elif directed and reverse:
                done = has_edge(rev_sd, path[0], path[n - 1]) != NULL
            else:
                done = has_edge(sd, path[n - 1], path[0]) != NULL

        if bigcount * reset_bound > max_iter:
            output = [int_to_vertex[longest_path[i]] for i in range(longest)]
            free_short_digraph(sd)
            if directed:
                free_short_digraph(rev_sd)
                if longest_reversed:
                    return (False, output[::-1])
            return (False, output)
    # #
    # # Output test
    # #

    if directed and reverse:
        # We revert the path to work on sd
        for i in range(length//2):
            t = path[i]
            path[i] = path[length - i - 1]
            path[length - i - 1] = t

    # Test adjacencies
    cdef bint good = True
    for i in range(n - 1):
        u = path[i]
        v = path[i + 1]
        if has_edge(sd, u, v) == NULL:
            good = False
            break
    if good is False:
        raise RuntimeError(f"vertices {int_to_vertex[u]} and {int_to_vertex[v]}"
                           " are consecutive in the cycle but are not adjacent")
    if not find_path and has_edge(sd, path[n - 1], path[0]) == NULL:
        raise RuntimeError(f"vertices {int_to_vertex[path[n - 1]]} and "
                           f"{int_to_vertex[path[0]]} are not adjacent")

    output = [int_to_vertex[path[i]] for i in range(length)]
    free_short_digraph(sd)
    if directed:
        free_short_digraph(rev_sd)

    return (True, output)


def transitive_reduction_acyclic(G):
    r"""
    Return the transitive reduction of an acyclic digraph.

    INPUT:

    - ``G`` -- an acyclic digraph

    EXAMPLES::

        sage: from sage.graphs.generic_graph_pyx import transitive_reduction_acyclic
        sage: G = posets.BooleanLattice(4).hasse_diagram()
        sage: G == transitive_reduction_acyclic(G.transitive_closure())
        True
    """
    cdef int  n = G.order()
    cdef dict v_to_int = {vv: i for i, vv in enumerate(G)}
    cdef int  u, v, i

    cdef list linear_extension

    is_acyclic, linear_extension = G.is_directed_acyclic(certificate=True)
    if not is_acyclic:
        raise ValueError("The graph is not directed acyclic")

    linear_extension.reverse()

    cdef binary_matrix_t closure

    # Build the transitive closure of G
    #
    # A point is reachable from u if it is one of its neighbours, or if it is
    # reachable from one of its neighbours.
    binary_matrix_init(closure, n, n)
    for uu in linear_extension:
        u = v_to_int[uu]
        for vv in G.neighbors_out(uu):
            v = v_to_int[vv]
            binary_matrix_set1(closure, u, v)
            bitset_or(closure.rows[u], closure.rows[u], closure.rows[v])

    # Build the transitive reduction of G
    #
    # An edge uv belongs to the transitive reduction of G if no outneighbor of u
    # can reach v (except v itself, of course).
    linear_extension.reverse()
    cdef list useful_edges = []
    for uu in linear_extension:
        u = v_to_int[uu]
        for vv in G.neighbors_out(uu):
            v = v_to_int[vv]
            bitset_difference(closure.rows[u], closure.rows[u], closure.rows[v])
        for vv in G.neighbors_out(uu):
            v = v_to_int[vv]
            if binary_matrix_get(closure, u, v):
                useful_edges.append((uu, vv))

    from sage.graphs.digraph import DiGraph
    reduced = DiGraph()
    reduced.add_edges(useful_edges)
    reduced.add_vertices(linear_extension)

    binary_matrix_free(closure)

    return reduced
