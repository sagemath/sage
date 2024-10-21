# sage_setup: distribution = sagemath-graphs
r"""
Chessboard graphs

The methods defined here appear in :mod:`sage.graphs.graph_generators`.

- :meth:`BishopGraph <sage.graphs.graph_generators.GraphGenerators.BishopGraph>`
- :meth:`KingGraph <sage.graphs.graph_generators.GraphGenerators.KingGraph>`
- :meth:`KnightGraph <sage.graphs.graph_generators.GraphGenerators.KnightGraph>`
- :meth:`QueenGraph <sage.graphs.graph_generators.GraphGenerators.QueenGraph>`
- :meth:`RookGraph <sage.graphs.graph_generators.GraphGenerators.RookGraph>`

AUTHORS:

- David Coudert    (2012)
"""

# ****************************************************************************
#           Copyright (C) 2012 David Coudert <david.coudert@inria.fr>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
# ****************************************************************************

from itertools import product
from itertools import combinations
from sage.graphs.graph import Graph


def ChessboardGraphGenerator(dim_list, rook=True, rook_radius=None,
                             bishop=True, bishop_radius=None,
                             knight=True, knight_x=1, knight_y=2,
                             relabel=False):
    r"""
    Return a Graph built on a `d`-dimensional chessboard with prescribed
    dimensions and interconnections.

    This function allows to generate many kinds of graphs corresponding to legal
    movements on a `d`-dimensional chessboard: Queen Graph, King Graph, Knight
    Graphs, Bishop Graph, and many generalizations. It also allows to avoid
    redundant code.

    INPUT:

    - ``dim_list`` -- iterable (list, set, dict); provides the dimensions
      `n_1, n_2, \ldots, n_d`, with `n_i \geq 1`, of the chessboard

    - ``rook`` -- boolean (default: ``True``); indicates whether the chess piece
      is able to move as a rook, that is at any distance along a dimension

    - ``rook_radius`` -- integer (default: ``None``); restriction on the
      rook-like movements to distance at most ``rook_radius``

    - ``bishop`` -- boolean (default: ``True``); indicates whether the chess
      piece is able to move like a bishop, that is along diagonals

    - ``bishop_radius`` -- integer (default: ``None``); restriction on the
      bishop-like movements to distance at most ``bishop_radius``

    - ``knight`` -- boolean (default: ``True``); indicating whether the chess
      piece is able to move like a knight

    - ``knight_x`` -- integer (default: `1`); indicates the number on steps
      the chess piece moves in one dimension when moving like a knight

    - ``knight_y`` -- integer (default: `2`); indicates the number on steps
      the chess piece moves in the second dimension when moving like a knight

    - ``relabel`` -- boolean (default: ``False``); indicates whether the
      vertices must be relabeled as integers

    OUTPUT:

    - A Graph build on a `d`-dimensional chessboard with prescribed dimensions,
      and with edges according given parameters.

    - A string encoding the dimensions. This is mainly useful for providing
      names to graphs.

    EXAMPLES:

    A `(2,2)`-King Graph is isomorphic to the complete graph on 4 vertices::

        sage: G, _ = graphs.ChessboardGraphGenerator( [2,2] )
        sage: G.is_isomorphic( graphs.CompleteGraph(4) )
        True

    A Rook's Graph in 2 dimensions is isomorphic to the Cartesian product of 2
    complete graphs::

        sage: G, _ = graphs.ChessboardGraphGenerator([3,4], rook=True, rook_radius=None, bishop=False, knight=False)
        sage: H = (graphs.CompleteGraph(3)).cartesian_product(graphs.CompleteGraph(4))
        sage: G.is_isomorphic(H)
        True

    TESTS:

    Giving dimensions less than 2::

        sage: graphs.ChessboardGraphGenerator([0, 2])
        Traceback (most recent call last):
        ...
        ValueError: the dimensions must be positive integers larger than 1

    Giving non integer dimensions::

        sage: graphs.ChessboardGraphGenerator([4.5, 2])
        Traceback (most recent call last):
        ...
        ValueError: the dimensions must be positive integers larger than 1

    Giving too few dimensions::

        sage: graphs.ChessboardGraphGenerator([2])
        Traceback (most recent call last):
        ...
        ValueError: the chessboard must have at least 2 dimensions

    Giving a non-iterable object as first parameter::

        sage: graphs.ChessboardGraphGenerator(2, 3)
        Traceback (most recent call last):
        ...
        TypeError: the first parameter must be an iterable object

    Giving too small rook radius::

        sage: graphs.ChessboardGraphGenerator([2, 3], rook=True, rook_radius=0)
        Traceback (most recent call last):
        ...
        ValueError: the rook_radius must be either None or have an integer value >= 1

    Giving wrong values for knights movements::

        sage: graphs.ChessboardGraphGenerator([2, 3], rook=False, bishop=False, knight=True, knight_x=1, knight_y=-1)
        Traceback (most recent call last):
        ...
        ValueError: the knight_x and knight_y values must be integers of value >= 1
    """
    from sage.rings.integer_ring import ZZ

    # We decode the dimensions of the chessboard
    try:
        dim = list(dim_list)
        nb_dim = len(dim)
    except TypeError:
        raise TypeError('the first parameter must be an iterable object')
    if nb_dim < 2:
        raise ValueError('the chessboard must have at least 2 dimensions')
    if any(a not in ZZ or a < 1 for a in dim):
        raise ValueError('the dimensions must be positive integers larger than 1')
    dimstr = str(tuple(dim))

    # We check the radius toward neighbors
    if rook:
        if rook_radius is None:
            rook_radius = max(dim)
        elif rook_radius not in ZZ or rook_radius < 1:
            raise ValueError('the rook_radius must be either None or have an integer value >= 1')
    if bishop:
        if bishop_radius is None:
            bishop_radius = max(dim)
        elif bishop_radius not in ZZ or bishop_radius < 1:
            raise ValueError('the bishop_radius must be either None or have an integer value >= 1')
    if knight and (knight_x not in ZZ or knight_y not in ZZ or knight_x < 1 or knight_y < 1):
        raise ValueError('the knight_x and knight_y values must be integers of value >= 1')

    # We build the set of vertices of the d-dimensional chessboard
    V = [list(x) for x in list(product(*[range(_) for _ in dim]))]

    G = Graph()
    for u in V:
        uu = tuple(u)
        G.add_vertex(uu)

        if rook:
            # We add edges to vertices we can reach when moving in one dimension
            for d in range(nb_dim):
                v = u[:]
                for k in range(v[d] + 1, min(dim[d], v[d] + 1 + rook_radius)):
                    v[d] = k
                    G.add_edge(uu, tuple(v))

        if bishop or knight:
            # We add edges to vertices we can reach when moving in two dimensions
            for dx, dy in combinations(range(nb_dim), 2):
                n = dim[dx]
                m = dim[dy]
                v = u[:]
                i = u[dx]
                j = u[dy]

                if bishop:
                    # Diagonal
                    for k in range(1, min(n - i, m - j, bishop_radius + 1)):
                        v[dx] = i + k
                        v[dy] = j + k
                        G.add_edge(uu, tuple(v))

                    # Anti-diagonal
                    for k in range(min(i, m - j - 1, bishop_radius)):
                        v[dx] = i - k - 1
                        v[dy] = j + k + 1
                        G.add_edge(uu, tuple(v))

                if knight:
                    # Moving knight_x in one dimension and knight_y in another
                    # dimension
                    if i + knight_y < n:
                        if j + knight_x < m:
                            v[dx] = i + knight_y
                            v[dy] = j + knight_x
                            G.add_edge(uu, tuple(v))
                        if j - knight_x >= 0:
                            v[dx] = i + knight_y
                            v[dy] = j - knight_x
                            G.add_edge(uu, tuple(v))
                    if j + knight_y < m:
                        if i + knight_x < n:
                            v[dx] = i + knight_x
                            v[dy] = j + knight_y
                            G.add_edge(uu, tuple(v))
                        if i - knight_x >= 0:
                            v[dx] = i - knight_x
                            v[dy] = j + knight_y
                            G.add_edge(uu, tuple(v))

    if relabel:
        G.relabel(inplace=True)
    return G, dimstr


def QueenGraph(dim_list, radius=None, relabel=False):
    r"""
    Return the `d`-dimensional Queen Graph with prescribed dimensions.

    The 2-dimensional Queen Graph of parameters `n` and `m` is a graph with `nm`
    vertices in which each vertex represents a square in an `n \times m`
    chessboard, and each edge corresponds to a legal move by a queen.

    The `d`-dimensional Queen Graph with `d >= 2` has for vertex set the cells
    of a `d`-dimensional grid with prescribed dimensions, and each edge
    corresponds to a legal move by a queen in either one or two dimensions.

    All 2-dimensional Queen Graphs are Hamiltonian and biconnected. The
    chromatic number of a `(n,n)`-Queen Graph is at least `n`, and it is exactly
    `n` when `n\equiv 1,5 \bmod{6}`.

    INPUT:

    - ``dim_list`` -- iterable (list, set, dict); provides the dimensions
      `n_1, n_2, \ldots, n_d`, with `n_i \geq 1`, of the chessboard

    - ``radius`` -- integer (default: ``None``); by setting the radius to a
      positive integer, one may reduce the visibility of the queen to at most
      ``radius`` steps. When radius is 1, the resulting graph is a King Graph.

    - ``relabel`` -- boolean (default: ``False``); indicates whether the
      vertices must be relabeled as integers

    EXAMPLES:

    The `(2,2)`-Queen Graph is isomorphic to the complete graph on 4 vertices::

        sage: G = graphs.QueenGraph([2, 2])
        sage: G.is_isomorphic(graphs.CompleteGraph(4))
        True

    The Queen Graph with radius 1 is isomorphic to the King Graph::

        sage: G = graphs.QueenGraph([4, 5], radius=1)
        sage: H = graphs.KingGraph([5, 4])
        sage: G.is_isomorphic(H)
        True

    Also ``True`` in higher dimensions::

        sage: G = graphs.QueenGraph([3, 4, 5], radius=1)
        sage: H = graphs.KingGraph([5, 3, 4])
        sage: G.is_isomorphic(H)
        True

    The Queen Graph can be obtained from the Rook Graph and the Bishop Graph::

        sage: for d in range(3,12):   # long time
        ....:     for r in range(1,d+1):
        ....:         G = graphs.QueenGraph([d,d],radius=r)
        ....:         H = graphs.RookGraph([d,d],radius=r)
        ....:         B = graphs.BishopGraph([d,d],radius=r)
        ....:         H.add_edges(B.edges(sort=False))
        ....:         if not G.is_isomorphic(H):
        ....:             print("that's not good!")
    """
    G, dimstr = ChessboardGraphGenerator(dim_list,
                                         rook=True, rook_radius=radius,
                                         bishop=True, bishop_radius=radius,
                                         knight=False,
                                         relabel=relabel)
    if radius is None:
        G.name(dimstr + "-Queen Graph")
    else:
        G.name(dimstr + "-Queen Graph with radius {}".format(radius))
    return G


def KingGraph(dim_list, radius=None, relabel=False):
    r"""
    Return the `d`-dimensional King Graph with prescribed dimensions.

    The 2-dimensional King Graph of parameters `n` and `m` is a graph with `nm`
    vertices in which each vertex represents a square in an `n \times m`
    chessboard, and each edge corresponds to a legal move by a king.

    The d-dimensional King Graph with `d >= 2` has for vertex set the cells of a
    d-dimensional grid with prescribed dimensions, and each edge corresponds to
    a legal move by a king in either one or two dimensions.

    All 2-dimensional King Graphs are Hamiltonian, biconnected, and have
    chromatic number 4 as soon as both dimensions are larger or equal to 2.

    INPUT:

    - ``dim_list`` -- iterable (list, set, dict); provides the dimensions
      `n_1, n_2, \ldots, n_d`, with `n_i \geq 1`, of the chessboard

    - ``radius`` -- integer (default: ``None``); by setting the radius to a
      positive integer, one may increase the power of the king to at least
      ``radius`` steps. When the radius equals the higher size of the
      dimensions, the resulting graph is a Queen Graph.

    - ``relabel`` -- boolean (default: ``False``); indicates whether the
      vertices must be relabeled as integers

    EXAMPLES:

    The `(2,2)`-King Graph is isomorphic to the complete graph on 4 vertices::

        sage: G = graphs.QueenGraph( [2, 2] )
        sage: G.is_isomorphic( graphs.CompleteGraph(4) )
        True

    The King Graph with large enough radius is isomorphic to a Queen Graph::

        sage: G = graphs.KingGraph( [5, 4], radius=5 )
        sage: H = graphs.QueenGraph( [4, 5] )
        sage: G.is_isomorphic( H )
        True

    Also ``True`` in higher dimensions::

        sage: G = graphs.KingGraph( [2, 5, 4], radius=5 )
        sage: H = graphs.QueenGraph( [4, 5, 2] )
        sage: G.is_isomorphic( H )
        True
    """
    rook_radius = 1 if radius is None else radius
    bishop_radius = 1 if radius is None else radius
    G, dimstr = ChessboardGraphGenerator(dim_list,
                                         rook=True, rook_radius=rook_radius,
                                         bishop=True, bishop_radius=bishop_radius,
                                         knight=False,
                                         relabel=relabel)
    if radius is None:
        G.name(dimstr + "-King Graph")
    else:
        G.name(dimstr + "-King Graph with radius {}".format(radius))
    return G


def KnightGraph(dim_list, one=1, two=2, relabel=False):
    r"""
    Return the d-dimensional Knight Graph with prescribed dimensions.

    The 2-dimensional Knight Graph of parameters `n` and `m` is a graph with
    `nm` vertices in which each vertex represents a square in an `n \times m`
    chessboard, and each edge corresponds to a legal move by a knight.

    The d-dimensional Knight Graph with `d >= 2` has for vertex set the cells of
    a d-dimensional grid with prescribed dimensions, and each edge corresponds
    to a legal move by a knight in any pairs of dimensions.

    The `(n,n)`-Knight Graph is Hamiltonian for even `n > 4`.

    INPUT:

    - ``dim_list`` -- iterable (list, set, dict); provides the dimensions
      `n_1, n_2, \ldots, n_d`, with `n_i \geq 1`, of the chessboard

    - ``one`` -- integer (default: `1`); indicates the number of steps in the
      first dimension

    - ``two`` -- integer (default: `2`); indicates the number of steps in the
      second dimension

    - ``relabel`` -- boolean (default: ``False``); indicates whether the
      vertices must be relabeled as integers

    EXAMPLES:

    The `(3,3)`-Knight Graph has an isolated vertex::

        sage: G = graphs.KnightGraph( [3, 3] )
        sage: G.degree( (1,1) )
        0

    The `(3,3)`-Knight Graph minus vertex (1,1) is a cycle of order 8::

        sage: G = graphs.KnightGraph( [3, 3] )
        sage: G.delete_vertex( (1,1) )
        sage: G.is_isomorphic( graphs.CycleGraph(8) )
        True

    The `(6,6)`-Knight Graph is Hamiltonian::

        sage: G = graphs.KnightGraph( [6, 6] )
        sage: G.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
    """
    G, dimstr = ChessboardGraphGenerator(dim_list,
                                         rook=False, bishop=False,
                                         knight=True, knight_x=one, knight_y=two,
                                         relabel=relabel)
    if one + two == 3:
        G.name(dimstr + "-Knight Graph")
    else:
        G.name(dimstr + "-Knight Graph with edges at distance ({}, {})".format(one, two))
    return G


def RookGraph(dim_list, radius=None, relabel=False):
    r"""
    Return the `d`-dimensional Rook's Graph with prescribed dimensions.

    The 2-dimensional Rook's Graph of parameters `n` and `m` is a graph with
    `nm` vertices in which each vertex represents a square in an `n \times m`
    chessboard, and each edge corresponds to a legal move by a rook.

    The `d`-dimensional Rook Graph with `d >= 2` has for vertex set the cells of
    a `d`-dimensional grid with prescribed dimensions, and each edge corresponds
    to a legal move by a rook in any of the dimensions.

    The Rook's Graph for an `n\times m` chessboard may also be defined as the
    Cartesian product of two complete graphs `K_n \square K_m`.

    INPUT:

    - ``dim_list`` -- iterable (list, set, dict); provides the dimensions
      `n_1, n_2, \ldots, n_d`, with `n_i \geq 1`, of the chessboard

    - ``radius`` -- integer (default: ``None``); by setting the radius to a
      positive integer, one may decrease the power of the rook to at most
      ``radius`` steps. When the radius is 1, the resulting graph is a
      `d`-dimensional grid.

    - ``relabel`` -- boolean (default: ``False``); indicates whether the
      vertices must be relabeled as integers

    EXAMPLES:

    The `(n,m)`-Rook's Graph is isomorphic to the Cartesian product of two
    complete graphs::

        sage: G = graphs.RookGraph( [3, 4] )
        sage: H = ( graphs.CompleteGraph(3) ).cartesian_product( graphs.CompleteGraph(4) )
        sage: G.is_isomorphic( H )
        True

    When the radius is 1, the Rook's Graph is a grid::

        sage: G = graphs.RookGraph( [3, 3, 4], radius=1 )
        sage: H = graphs.GridGraph( [3, 4, 3] )
        sage: G.is_isomorphic( H )
        True
    """
    G, dimstr = ChessboardGraphGenerator(dim_list,
                                         rook=True, rook_radius=radius,
                                         bishop=False, knight=False,
                                         relabel=relabel)
    if radius is None:
        G.name(dimstr + "-Rook Graph")
    else:
        G.name(dimstr + "-Rook Graph with radius {}".format(radius))
    return G


def BishopGraph(dim_list, radius=None, relabel=False):
    r"""
    Return the `d`-dimensional Bishop Graph with prescribed dimensions.

    The 2-dimensional Bishop Graph of parameters `n` and `m` is a graph with
    `nm` vertices in which each vertex represents a square in an `n \times m`
    chessboard, and each edge corresponds to a legal move by a bishop.

    The `d`-dimensional Bishop Graph with `d >= 2` has for vertex set the cells
    of a `d`-dimensional grid with prescribed dimensions, and each edge
    corresponds to a legal move by a bishop in any pairs of dimensions.

    The Bishop Graph is not connected.

    INPUT:

    - ``dim_list`` -- iterable (list, set, dict); provides the dimensions
      `n_1, n_2, \ldots, n_d`, with `n_i \geq 1`, of the chessboard

    - ``radius`` -- integer (default: ``None``); by setting the radius to a
      positive integer, one may decrease the power of the bishop to at most
      ``radius`` steps.

    - ``relabel`` -- boolean (default: ``False``); indicates whether the
      vertices must be relabeled as integers

    EXAMPLES:

    The (n,m)-Bishop Graph is not connected::

        sage: G = graphs.BishopGraph( [3, 4] )
        sage: G.is_connected()
        False

    The Bishop Graph can be obtained from Knight Graphs::

        sage: for d in range(3,12):   # long time
        ....:     H = Graph()
        ....:     for r in range(1,d+1):
        ....:         B = graphs.BishopGraph([d,d],radius=r)
        ....:         H.add_edges( graphs.KnightGraph([d,d],one=r,two=r).edges(sort=False) )
        ....:         if not B.is_isomorphic(H):
        ....:            print("that's not good!")
    """
    G, dimstr = ChessboardGraphGenerator(dim_list,
                                         rook=False, knight=False,
                                         bishop=True, bishop_radius=radius,
                                         relabel=relabel)
    if radius is None:
        G.name(dimstr + "-Bishop Graph")
    else:
        G.name(dimstr + "-Bishop Graph with radius {}".format(radius))
    return G
