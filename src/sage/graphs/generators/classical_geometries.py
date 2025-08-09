# sage.doctest: optional - sage.modules
r"""
Families of graphs derived from classical geometries over finite fields

These include graphs of polar spaces, affine polar graphs, graphs
related to Hermitean unitals, graphs on nonisotropic points, etc.

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""

# ****************************************************************************
#           Copyright (C) 2015 Sagemath project
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.graphs.graph import Graph
from sage.arith.misc import is_prime_power
from sage.rings.finite_rings.finite_field_constructor import FiniteField


def SymplecticPolarGraph(d, q, algorithm=None):
    r"""
    Return the Symplectic Polar Graph `Sp(d,q)`.

    The Symplectic Polar Graph `Sp(d,q)` is built from a projective space of
    dimension `d-1` over a field `F_q`, and a symplectic form `f`. Two vertices
    `u,v` are made adjacent if `f(u,v)=0`.

    See the page `on symplectic graphs on Andries Brouwer's website
    <https://www.win.tue.nl/~aeb/graphs/Sp.html>`_.

    INPUT:

    - ``d``, ``q`` -- integers; note that only even values of `d` are accepted
      by the function

    - ``algorithm`` -- string (default: ``None``); if set to ``'gap'``, then the
      computation is carried via GAP library interface, computing totally
      singular subspaces, which is faster for `q>3`.  Otherwise it is done
      directly.

    EXAMPLES:

    Computation of the spectrum of `Sp(6,2)`::

        sage: g = graphs.SymplecticPolarGraph(6, 2)
        sage: g.is_strongly_regular(parameters=True)
        (63, 30, 13, 15)
        sage: set(g.spectrum()) == {-5, 3, 30}                                          # needs sage.rings.number_field
        True

    The parameters of `Sp(4,q)` are the same as of `O(5,q)`, but they are
    not isomorphic if `q` is odd::

        sage: G = graphs.SymplecticPolarGraph(4, 3)
        sage: G.is_strongly_regular(parameters=True)
        (40, 12, 2, 4)

        sage: # needs sage.libs.gap
        sage: O = graphs.OrthogonalPolarGraph(5, 3)
        sage: O.is_strongly_regular(parameters=True)
        (40, 12, 2, 4)
        sage: O.is_isomorphic(G)
        False
        sage: S = graphs.SymplecticPolarGraph(6, 4, algorithm='gap')    # not tested (long time)
        sage: S.is_strongly_regular(parameters=True)                    # not tested (long time)
        (1365, 340, 83, 85)

    TESTS::

        sage: graphs.SymplecticPolarGraph(4,4,algorithm='gap').is_strongly_regular(parameters=True)                     # needs sage.libs.gap
        (85, 20, 3, 5)
        sage: graphs.SymplecticPolarGraph(4,4).is_strongly_regular(parameters=True)     # needs sage.libs.pari
        (85, 20, 3, 5)
        sage: graphs.SymplecticPolarGraph(4,4,algorithm='blah')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm!
    """
    if d < 1 or d % 2:
        raise ValueError("d must be even and greater than 2")

    if algorithm == "gap":     # faster for larger (q>3)  fields
        from sage.libs.gap.libgap import libgap
        G = _polar_graph(d, q, libgap.SymplecticGroup(d, q))

    elif algorithm is None:    # faster for small (q<4) fields
        from sage.modules.free_module import VectorSpace
        from sage.schemes.projective.projective_space import ProjectiveSpace
        from sage.matrix.constructor import identity_matrix, block_matrix, zero_matrix

        F = FiniteField(q, "x")
        M = block_matrix(F, 2, 2,
                         [zero_matrix(F, d/2),
                          identity_matrix(F, d/2),
                          -identity_matrix(F, d/2),
                          zero_matrix(F, d/2)])

        V = VectorSpace(F, d)
        PV = list(ProjectiveSpace(d - 1, F))
        G = Graph([[tuple(_) for _ in PV], lambda x, y: V(x)*(M*V(y)) == 0], loops=False)

    else:
        raise ValueError("unknown algorithm!")

    G.name("Symplectic Polar Graph Sp({},{})".format(d, q))
    G.relabel()
    return G


def AffineOrthogonalPolarGraph(d, q, sign='+'):
    r"""
    Return the affine polar graph `VO^+(d,q),VO^-(d,q)` or `VO(d,q)`.

    Affine Polar graphs are built from a `d`-dimensional vector space over
    `F_q`, and a quadratic form which is hyperbolic, elliptic or parabolic
    according to the value of ``sign``.

    Note that `VO^+(d,q),VO^-(d,q)` are strongly regular graphs, while `VO(d,q)`
    is not.

    For more information on Affine Polar graphs, see `Affine Polar Graphs page
    of Andries Brouwer's website <https://www.win.tue.nl/~aeb/graphs/VO.html>`_.

    INPUT:

    - ``d`` -- integer; ``d`` must be even if ``sign is not None``, and odd
      otherwise

    - ``q`` -- integer; a power of a prime number, as `F_q` must exist

    - ``sign`` -- string (default: ``'+'``); must be equal to ``'+'``, ``'-'``,
      or ``None`` to compute (respectively) `VO^+(d,q),VO^-(d,q)` or
      `VO(d,q)`

    .. NOTE::

        The graph `VO^\epsilon(d,q)` is the graph induced by the
        non-neighbors of a vertex in an :meth:`Orthogonal Polar Graph
        <OrthogonalPolarGraph>` `O^\epsilon(d+2,q)`.

    EXAMPLES:

    The :meth:`Brouwer-Haemers graph <BrouwerHaemersGraph>` is isomorphic to
    `VO^-(4,3)`::

        sage: g = graphs.AffineOrthogonalPolarGraph(4,3,"-")                            # needs sage.libs.gap
        sage: g.is_isomorphic(graphs.BrouwerHaemersGraph())                             # needs sage.libs.gap
        True

    Some examples from `Brouwer's table or strongly regular graphs
    <https://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`_::

        sage: # needs sage.libs.gap
        sage: g = graphs.AffineOrthogonalPolarGraph(6,2,"-"); g
        Affine Polar Graph VO^-(6,2): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 27, 10, 12)
        sage: g = graphs.AffineOrthogonalPolarGraph(6,2,"+"); g
        Affine Polar Graph VO^+(6,2): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 35, 18, 20)

    When ``sign is None``::

        sage: # needs sage.libs.gap
        sage: g = graphs.AffineOrthogonalPolarGraph(5,2,None); g
        Affine Polar Graph VO^-(5,2): Graph on 32 vertices
        sage: g.is_strongly_regular(parameters=True)
        False
        sage: g.is_regular()
        True
        sage: g.is_vertex_transitive()
        True
    """
    if sign in ["+", "-"]:
        s = 1 if sign == "+" else -1
        if d % 2:
            raise ValueError("d must be even when sign is not None")
    else:
        if not d % 2:
            raise ValueError("d must be odd when sign is None")
        s = 0

    from sage.modules.free_module import VectorSpace
    from sage.matrix.constructor import Matrix
    from sage.libs.gap.libgap import libgap
    from itertools import combinations

    M = Matrix(libgap.InvariantQuadraticForm(libgap.GeneralOrthogonalGroup(s, d, q))['matrix'])
    F = libgap.GF(q).sage()
    V = list(VectorSpace(F, d))

    G = Graph()
    G.add_vertices([tuple(_) for _ in V])
    for x, y in combinations(V, 2):
        if not (x - y)*M*(x - y):
            G.add_edge(tuple(x), tuple(y))

    G.name("Affine Polar Graph VO^" + str('+' if s == 1 else '-') + "(" + str(d) + "," + str(q) + ")")
    G.relabel()
    return G


def _orthogonal_polar_graph(m, q, sign='+', point_type=[0]):
    r"""
    A helper function to build ``OrthogonalPolarGraph`` and ``NO2,3,5`` graphs.

    See the `page of
    Andries Brouwer's website <https://www.win.tue.nl/~aeb/graphs/srghub.html>`_.

    INPUT:

    - ``m``, ``q`` -- integers; `q` must be a prime power

    - ``sign`` -- string (default: ``'+'``); must be ``'+'`` or ``'-'`` if `m`
      is even, ``'+'`` (default) otherwise

    - ``point_type`` -- list of elements from `F_q`

    EXAMPLES:

    Petersen graph::

        sage: from sage.graphs.generators.classical_geometries import _orthogonal_polar_graph
        sage: g = _orthogonal_polar_graph(3,5,point_type=[2,3])                         # needs sage.libs.gap
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.gap
        (10, 3, 0, 1)

    A locally Petersen graph (a.k.a. Doro graph, a.k.a. Hall graph)::

        sage: g = _orthogonal_polar_graph(4,5,'-',point_type=[2,3])                     # needs sage.libs.gap
        sage: g.is_distance_regular(parameters=True)                                    # needs sage.libs.gap
        ([10, 6, 4, None], [None, 1, 2, 5])

    Various big and slow to build graphs:

    `NO^+(7,3)`::

        sage: g = _orthogonal_polar_graph(7,3,point_type=[1])    # not tested (long time)
        sage: g.is_strongly_regular(parameters=True)             # not tested (long time)
        (378, 117, 36, 36)

    `NO^-(7,3)`::

        sage: g = _orthogonal_polar_graph(7,3,point_type=[-1])   # not tested (long time)
        sage: g.is_strongly_regular(parameters=True)             # not tested (long time)
        (351, 126, 45, 45)

    `NO^+(6,3)`::

        sage: g = _orthogonal_polar_graph(6,3,point_type=[1])                           # needs sage.libs.gap
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.gap
        (117, 36, 15, 9)

    `NO^-(6,3)`::

        sage: g = _orthogonal_polar_graph(6,3,'-',point_type=[1])                       # needs sage.libs.gap
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.gap
        (126, 45, 12, 18)

    `NO^{-,\perp}(5,5)`::

        sage: g = _orthogonal_polar_graph(5,5,point_type=[2,3])         # long time, needs sage.libs.gap
        sage: g.is_strongly_regular(parameters=True)                    # long time, needs sage.libs.gap
        (300, 65, 10, 15)

    `NO^{+,\perp}(5,5)`::

        sage: g = _orthogonal_polar_graph(5,5,point_type=[1,-1]) # not tested (long time)
        sage: g.is_strongly_regular(parameters=True)             # not tested (long time)
        (325, 60, 15, 10)

    TESTS::

        sage: # needs sage.libs.gap
        sage: g = _orthogonal_polar_graph(5,3,point_type=[-1])
        sage: g.is_strongly_regular(parameters=True)
        (45, 12, 3, 3)
        sage: g = _orthogonal_polar_graph(5,3,point_type=[1])
        sage: g.is_strongly_regular(parameters=True)
        (36, 15, 6, 6)
    """
    from sage.schemes.projective.projective_space import ProjectiveSpace
    from sage.modules.free_module_element import free_module_element as vector
    from sage.matrix.constructor import Matrix
    from sage.libs.gap.libgap import libgap

    if not m % 2:
        if sign != "+" and sign != "-":
            raise ValueError("sign must be equal to either '-' or '+' when "
                             "m is even")
    else:
        if sign != "" and sign != "+":
            raise ValueError("sign must be equal to either '' or '+' when "
                             "m is odd")
        sign = ""

    e = {'+': 1,
         '-': -1,
         '': 0}[sign]

    M = Matrix(libgap.InvariantQuadraticForm(libgap.GeneralOrthogonalGroup(e, m, q))['matrix'])
    Fq = libgap.GF(q).sage()
    PG = [vector(s) for s in ProjectiveSpace(m - 1, Fq)]

    for v in PG:
        v.set_immutable()

    def F(x):
        return x*M*x

    if not q % 2:
        def P(x, y):
            return F(x - y)
    else:
        def P(x, y):
            return x*M*y + y*M*x

    V = [x for x in PG if F(x) in point_type]

    G = Graph([V, lambda x, y: P(x, y) == 0], loops=False)

    G.relabel()
    return G


def OrthogonalPolarGraph(m, q, sign='+'):
    r"""
    Return the Orthogonal Polar Graph `O^{\epsilon}(m,q)`.

    For more information on Orthogonal Polar graphs, see the `page of Andries
    Brouwer's website <https://www.win.tue.nl/~aeb/graphs/srghub.html>`_.

    INPUT:

    - ``m``, ``q`` -- integers; `q` must be a prime power

    - ``sign`` -- string (default: ``'+'``); must be ``'+'`` or ``'-'`` if `m`
      is even, ``'+'`` (default) otherwise

    EXAMPLES::

        sage: # needs sage.libs.gap
        sage: G = graphs.OrthogonalPolarGraph(6,3,"+"); G
        Orthogonal Polar Graph O^+(6, 3): Graph on 130 vertices
        sage: G.is_strongly_regular(parameters=True)
        (130, 48, 20, 16)
        sage: G = graphs.OrthogonalPolarGraph(6,3,"-"); G
        Orthogonal Polar Graph O^-(6, 3): Graph on 112 vertices
        sage: G.is_strongly_regular(parameters=True)
        (112, 30, 2, 10)
        sage: G = graphs.OrthogonalPolarGraph(5,3); G
        Orthogonal Polar Graph O(5, 3): Graph on 40 vertices
        sage: G.is_strongly_regular(parameters=True)
        (40, 12, 2, 4)
        sage: G = graphs.OrthogonalPolarGraph(8,2,"+"); G
        Orthogonal Polar Graph O^+(8, 2): Graph on 135 vertices
        sage: G.is_strongly_regular(parameters=True)
        (135, 70, 37, 35)
        sage: G = graphs.OrthogonalPolarGraph(8,2,"-"); G
        Orthogonal Polar Graph O^-(8, 2): Graph on 119 vertices
        sage: G.is_strongly_regular(parameters=True)
        (119, 54, 21, 27)

    TESTS::

        sage: G = graphs.OrthogonalPolarGraph(4,3,"")                                   # needs sage.libs.gap
        Traceback (most recent call last):
        ...
        ValueError: sign must be equal to either '-' or '+' when m is even
        sage: G = graphs.OrthogonalPolarGraph(5,3,"-")                                  # needs sage.libs.gap
        Traceback (most recent call last):
        ...
        ValueError: sign must be equal to either '' or '+' when m is odd
    """
    G = _orthogonal_polar_graph(m, q, sign=sign)
    if m % 2:
        sign = ""
    G.name("Orthogonal Polar Graph O" + ("^" + sign if sign else "") + str((m, q)))
    return G


def NonisotropicOrthogonalPolarGraph(m, q, sign='+', perp=None):
    r"""
    Return the Graph `NO^{\epsilon,\perp}_{m}(q)`.

    Let the vectorspace of dimension `m` over `F_q` be endowed with a
    nondegenerate quadratic form `F`, of type ``sign`` for `m` even.

    * `m` even: assume further that `q=2` or `3`. Returns the graph of the
      points (in the underlying projective space) `x` satisfying `F(x)=1`, with
      adjacency given by orthogonality w.r.t. `F`. Parameter ``perp`` is
      ignored.

    * `m` odd: if ``perp`` is not ``None``, then we assume that `q=5` and return
      the graph of the points `x` satisfying `F(x)=\pm 1` if ``sign="+"``,
      respectively `F(x) \in \{2,3\}` if ``sign="-"``, with adjacency given by
      orthogonality w.r.t. `F` (cf. Sect 7.D of [BL1984]_). Otherwise return
      the graph of nongenerate hyperplanes of type ``sign``, adjacent whenever
      the intersection is degenerate (cf. Sect. 7.C of [BL1984]_).
      Note that for `q=2` one will get a complete graph.

    For more information, see Sect. 9.9 of [BH2012]_ and [BL1984]_. Note that
    the `page of Andries Brouwer's website
    <https://www.win.tue.nl/~aeb/graphs/srghub.html>`_ uses different notation.

    INPUT:

    - ``m`` -- integer;  half the dimension of the underlying vectorspace

    - ``q`` -- a power of a prime number, the size of the underlying field

    - ``sign`` -- string (default: ``'+'``); must be either ``'+'`` or ``'-'``

    EXAMPLES:

    `NO^-(4,2)` is isomorphic to Petersen graph::

        sage: g = graphs.NonisotropicOrthogonalPolarGraph(4,2,'-'); g                   # needs sage.libs.gap
        NO^-(4, 2): Graph on 10 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.gap
        (10, 3, 0, 1)

    `NO^-(6,2)` and `NO^+(6,2)`::

        sage: # needs sage.libs.gap
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(6,2,'-')
        sage: g.is_strongly_regular(parameters=True)
        (36, 15, 6, 6)
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(6,2,'+'); g
        NO^+(6, 2): Graph on 28 vertices
        sage: g.is_strongly_regular(parameters=True)
        (28, 15, 6, 10)

    `NO^+(8,2)`::

        sage: g = graphs.NonisotropicOrthogonalPolarGraph(8,2,'+')                      # needs sage.libs.gap
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.gap
        (120, 63, 30, 36)

    Wilbrink's graphs for `q=5`::

        sage: # needs sage.libs.gap
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(5,5,perp=1)
        sage: g.is_strongly_regular(parameters=True)    # long time
        (325, 60, 15, 10)
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(5,5,'-',perp=1)
        sage: g.is_strongly_regular(parameters=True)    # long time
        (300, 65, 10, 15)

    Wilbrink's graphs::

        sage: # needs sage.libs.gap
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(5,4,'+')
        sage: g.is_strongly_regular(parameters=True)
        (136, 75, 42, 40)
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(5,4,'-')
        sage: g.is_strongly_regular(parameters=True)
        (120, 51, 18, 24)
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(7,4,'+'); g        # not tested (long time)
        NO^+(7, 4): Graph on 2080 vertices
        sage: g.is_strongly_regular(parameters=True)                         # not tested (long time)
        (2080, 1071, 558, 544)

    TESTS::

        sage: # needs sage.libs.gap
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(4,2); g
        NO^+(4, 2): Graph on 6 vertices
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(4,3,'-')
        sage: g.is_strongly_regular(parameters=True)
        (15, 6, 1, 3)
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(3,5,'-',perp=1); g
        NO^-,perp(3, 5): Graph on 10 vertices
        sage: g.is_strongly_regular(parameters=True)
        (10, 3, 0, 1)

        sage: # long time, needs sage.libs.gap
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(6,3,'+')
        sage: g.is_strongly_regular(parameters=True)
        (117, 36, 15, 9)
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(6,3,'-'); g
        NO^-(6, 3): Graph on 126 vertices
        sage: g.is_strongly_regular(parameters=True)
        (126, 45, 12, 18)
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(5,5,'-')
        sage: g.is_strongly_regular(parameters=True)
        (300, 104, 28, 40)
        sage: g = graphs.NonisotropicOrthogonalPolarGraph(5,5,'+')
        sage: g.is_strongly_regular(parameters=True)
        (325, 144, 68, 60)

        sage: g = graphs.NonisotropicOrthogonalPolarGraph(6,4,'+')
        Traceback (most recent call last):
        ...
        ValueError: for m even q must be 2 or 3
    """
    p, k = is_prime_power(q, get_data=True)
    if not k:
        raise ValueError('q must be a prime power')
    dec = ''
    if not m % 2:
        if q in [2, 3]:
            G = _orthogonal_polar_graph(m, q, sign=sign, point_type=[1])
        else:
            raise ValueError("for m even q must be 2 or 3")
    elif perp is not None:
        if q == 5:
            pt = [-1, 1] if sign == '+' else [2, 3] if sign == '-' else []
            G = _orthogonal_polar_graph(m, q, point_type=pt)
            dec = ",perp"
        else:
            raise ValueError("for perp not None q must be 5")
    else:
        if sign not in ['+', '-']:
            raise ValueError("sign must be '+' or '-'")
        from sage.libs.gap.libgap import libgap
        g0 = libgap.GeneralOrthogonalGroup(m, q)
        g = libgap.Group(libgap.List(g0.GeneratorsOfGroup(), libgap.TransposedMat))
        F = libgap.GF(q)  # F_q
        W = libgap.FullRowSpace(F, m)  # F_q^m
        e = 1 if sign == '+' else -1
        n = (m - 1)/2
        # we build (q^n(q^n+e)/2, (q^n-e)(q^(n-1)+e), 2(q^(2n-2)-1)+eq^(n-1)(q-1),
        #                                          2q^(n-1)(q^(n-1)+e))-srg
        # **use** v and k to select appropriate orbit and orbital
        nvert = (q**n)*(q**n + e)/2  # v
        deg = (q**n - e)*(q**(n - 1) + e)   # k
        S = [libgap.Elements(libgap.Basis(x))[0]
             for x in libgap.Elements(libgap.Subspaces(W, 1))]
        (V,) = (x for x in libgap.Orbits(g, S, libgap.OnLines)
                if len(x) == nvert)
        gp = libgap.Action(g, V, libgap.OnLines)  # make a permutation group
        h = libgap.Stabilizer(gp, 1)
        (Vh,) = (x for x in libgap.Orbits(h, libgap.Orbit(gp, 1))
                 if len(x) == deg)
        Vh = Vh[0]
        L = libgap.Orbit(gp, [1, Vh], libgap.OnSets)
        G = Graph()
        G.add_edges(L)
    G.name("NO^" + sign + dec + str((m, q)))
    return G


def _polar_graph(m, q, g, intersection_size=None):
    r"""
    The helper function to build graphs `(D)U(m,q)` and `(D)Sp(m,q)`.

    Building a graph on an orbit of a group `g` of `m\times m` matrices over
    `GF(q)` on the points (or subspaces of dimension ``m//2``) isotropic
    w.r.t. the form `F` left invariant by the group `g`.

    The only constraint is that the first ``m//2`` elements of the standard
    basis must generate a totally isotropic w.r.t. `F` subspace; this is the
    case with these groups coming from GAP; namely, `F` has the anti-diagonal
    all-1 matrix.

    INPUT:

    - ``m`` -- the dimension of the underlying vector space

    - ``q`` -- the size of the field

    - ``g`` -- the group acting

    - ``intersection_size`` -- (default: ``None``) if ``None``, build the graph
      on the isotropic points, with adjacency being orthogonality w.r.t. `F`.
      Otherwise, build the graph on the maximal totally isotropic subspaces,
      with adjacency specified by ``intersection_size`` being as given.

    TESTS::

        sage: from sage.graphs.generators.classical_geometries import _polar_graph
        sage: _polar_graph(4, 4, libgap.GeneralUnitaryGroup(4, 2))                                  # needs sage.libs.gap
        Graph on 45 vertices
        sage: _polar_graph(4, 4, libgap.GeneralUnitaryGroup(4, 2), intersection_size=1)             # needs sage.libs.gap
        Graph on 27 vertices
    """
    from sage.libs.gap.libgap import libgap
    from itertools import combinations
    W = libgap.FullRowSpace(libgap.GF(q), m)   # F_q^m
    B = libgap.Elements(libgap.Basis(W))       # the standard basis of W
    V = libgap.Orbit(g, B[0], libgap.OnLines)  # orbit on isotropic points
    gp = libgap.Action(g, V, libgap.OnLines)   # make a permutation group
    s = libgap.Subspace(W, [B[i] for i in range(m//2)])  # a totally isotropic subspace
    # and the points there
    sp = [libgap.Elements(libgap.Basis(x))[0] for x in libgap.Elements(s.Subspaces(1))]
    h = libgap.Set([libgap.Position(V, x) for x in sp])  # indices of the points in s
    L = [tuple(o) for o in libgap.Orbit(gp, h, libgap.OnSets)]  # orbit on these subspaces
    if intersection_size is None:
        G = Graph()
        for x in L:  # every pair of points in the subspace is adjacent to each other in G
            G.add_edges(combinations(x, 2))
        return G
    else:
        return Graph([L, lambda i, j: libgap.Size(libgap.Intersection(i, j)) == intersection_size],
                     loops=False)


def UnitaryPolarGraph(m, q, algorithm='gap'):
    r"""
    Return the Unitary Polar Graph `U(m,q)`.

    For more information on Unitary Polar graphs, see the `page of Andries
    Brouwer's website <https://www.win.tue.nl/~aeb/graphs/srghub.html>`_.

    INPUT:

    - ``m``, ``q`` -- integers; `q` must be a prime power

    - ``algorithm`` -- string (default: ``'gap'``); if set to ``'gap'`` then the
      computation is carried via GAP library interface, computing totally
      singular subspaces, which is faster for large examples (especially with
      `q>2`). Otherwise it is done directly.

    EXAMPLES::

        sage: # needs sage.libs.gap
        sage: G = graphs.UnitaryPolarGraph(4,2); G
        Unitary Polar Graph U(4, 2); GQ(4, 2): Graph on 45 vertices
        sage: G.is_strongly_regular(parameters=True)
        (45, 12, 3, 3)
        sage: graphs.UnitaryPolarGraph(5,2).is_strongly_regular(parameters=True)
        (165, 36, 3, 9)
        sage: graphs.UnitaryPolarGraph(6,2)     # not tested (long time)
        Unitary Polar Graph U(6, 2): Graph on 693 vertices

    TESTS::

        sage: graphs.UnitaryPolarGraph(4,3, algorithm='gap').is_strongly_regular(parameters=True)   # needs sage.libs.gap
        (280, 36, 8, 4)
        sage: graphs.UnitaryPolarGraph(4,3).is_strongly_regular(parameters=True)                    # needs sage.libs.gap
        (280, 36, 8, 4)
        sage: graphs.UnitaryPolarGraph(4,3, algorithm='foo')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm!
    """
    if algorithm == "gap":
        from sage.libs.gap.libgap import libgap
        G = _polar_graph(m, q**2, libgap.GeneralUnitaryGroup(m, q))

    elif algorithm is None:  # slow on large examples
        from sage.schemes.projective.projective_space import ProjectiveSpace
        from sage.modules.free_module_element import free_module_element as vector
        Fq = FiniteField(q**2, 'a')
        PG = map(vector, ProjectiveSpace(m - 1, Fq))

        for v in PG:
            v.set_immutable()

        def P(x, y):
            return sum(x[j] * y[m - 1 - j] ** q for j in range(m)) == 0

        V = [x for x in PG if P(x, x)]
        # bottleneck is here, of course
        G = Graph([V, lambda x, y: P(x, y)], loops=False)
    else:
        raise ValueError("unknown algorithm!")

    G.relabel()
    G.name("Unitary Polar Graph U" + str((m, q)))
    if m == 4:
        G.name(G.name() + '; GQ' + str((q**2, q)))
    if m == 5:
        G.name(G.name() + '; GQ' + str((q**2, q**3)))
    return G


def NonisotropicUnitaryPolarGraph(m, q):
    r"""
    Return the Graph `NU(m,q)`.

    This returns the graph on nonisotropic, with respect to a nondegenerate
    Hermitean form, points of the `(m-1)`-dimensional projective space over
    `F_q`, with points adjacent whenever they lie on a tangent (to the set of
    isotropic points) line.  For more information, see Sect. 9.9 of [BH2012]_
    and series C14 in [Hub1975]_.

    INPUT:

    - ``m``, ``q`` -- integers; `q` must be a prime power

    EXAMPLES::

        sage: g = graphs.NonisotropicUnitaryPolarGraph(5,2); g                          # needs sage.libs.gap
        NU(5, 2): Graph on 176 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.gap
        (176, 135, 102, 108)

    TESTS::

        sage: graphs.NonisotropicUnitaryPolarGraph(4,2).is_strongly_regular(parameters=True)        # needs sage.libs.gap
        (40, 27, 18, 18)
        sage: graphs.NonisotropicUnitaryPolarGraph(4,3).is_strongly_regular(parameters=True)  # long time, needs sage.libs.gap
        (540, 224, 88, 96)
        sage: graphs.NonisotropicUnitaryPolarGraph(6,6)
        Traceback (most recent call last):
        ...
        ValueError: q must be a prime power
    """
    p, k = is_prime_power(q, get_data=True)
    if not k:
        raise ValueError('q must be a prime power')
    from sage.libs.gap.libgap import libgap
    from itertools import combinations
    F = libgap.GF(q**2)  # F_{q^2}
    W = libgap.FullRowSpace(F, m)  # F_{q^2}^m
    B = libgap.Elements(libgap.Basis(W))  # the standard basis of W
    if m % 2:
        point = B[(m - 1)/2]
    else:
        if p == 2:
            point = B[m/2] + F.PrimitiveRoot()*B[(m - 2)/2]
        else:
            point = B[(m - 2)/2] + B[m/2]
    g = libgap.GeneralUnitaryGroup(m, q)
    V = libgap.Orbit(g, point, libgap.OnLines)  # orbit on nonisotropic points
    gp = libgap.Action(g, V, libgap.OnLines)  # make a permutation group

    s = libgap.Subspace(W, [point, point + B[0]])  # a tangent line on point

    # and the points there
    sp = [libgap.Elements(libgap.Basis(x))[0] for x in libgap.Elements(s.Subspaces(1))]
    h = libgap.Set([libgap.Position(V, x)
                    for x in libgap.Intersection(V, sp)])  # indices
    L = libgap.Orbit(gp, h, libgap.OnSets)  # orbit on the tangent lines
    G = Graph()
    for x in L:  # every pair of points in the subspace is adjacent to each other in G
        G.add_edges(combinations(x, 2))
    G.relabel()
    G.name("NU" + str((m, q)))
    return G


def UnitaryDualPolarGraph(m, q):
    r"""
    Return the Dual Unitary Polar Graph `U(m,q)`.

    For more information on Unitary Dual Polar graphs, see [BCN1989]_ and
    Sect. 2.3.1 of [Coh1981]_.

    INPUT:

    - ``m``, ``q`` -- integers; `q` must be a prime power

    EXAMPLES:

    The point graph of a generalized quadrangle (see
    :wikipedia:`Generalized_quadrangle`, [PT2009]_) of order (8,4)::

        sage: G = graphs.UnitaryDualPolarGraph(5,2); G  # long time                     # needs sage.libs.gap
        Unitary Dual Polar Graph DU(5, 2); GQ(8, 4): Graph on 297 vertices
        sage: G.is_strongly_regular(parameters=True)    # long time                     # needs sage.libs.gap
        (297, 40, 7, 5)

    Another way to get the  generalized quadrangle of order (2,4)::

        sage: G = graphs.UnitaryDualPolarGraph(4,2); G                                  # needs sage.libs.gap
        Unitary Dual Polar Graph DU(4, 2); GQ(2, 4): Graph on 27 vertices
        sage: G.is_isomorphic(graphs.OrthogonalPolarGraph(6,2,'-'))                     # needs sage.libs.gap
        True

    A bigger graph::

        sage: G = graphs.UnitaryDualPolarGraph(6,2); G   # not tested (long time)
        Unitary Dual Polar Graph DU(6, 2): Graph on 891 vertices
        sage: G.is_distance_regular(parameters=True)     # not tested (long time)
        ([42, 40, 32, None], [None, 1, 5, 21])

    TESTS::

        sage: graphs.UnitaryDualPolarGraph(6,6)                                         # needs sage.libs.gap
        Traceback (most recent call last):
        ...
        GAPError: Error, <subfield> must be a prime or a finite field
    """
    from sage.libs.gap.libgap import libgap
    G = _polar_graph(m, q**2, libgap.GeneralUnitaryGroup(m, q),
                     intersection_size=int((q**(2*(m//2 - 1)) - 1)/(q**2 - 1)))
    G.relabel()
    G.name("Unitary Dual Polar Graph DU" + str((m, q)))
    if m == 4:
        G.name(G.name() + '; GQ' + str((q, q**2)))
    if m == 5:
        G.name(G.name() + '; GQ' + str((q**3, q**2)))
    return G


def SymplecticDualPolarGraph(m, q):
    r"""
    Return the Symplectic Dual Polar Graph `DSp(m,q)`.

    For more information on Symplectic Dual Polar graphs, see [BCN1989]_ and
    Sect. 2.3.1 of [Coh1981]_.

    INPUT:

    - ``m``, ``q`` -- integers; `q` must be a prime power, and `m` must be even

    EXAMPLES::

        sage: G = graphs.SymplecticDualPolarGraph(6,3); G       # not tested (long time)
        Symplectic Dual Polar Graph DSp(6, 3): Graph on 1120 vertices
        sage: G.is_distance_regular(parameters=True)            # not tested (long time)
        ([39, 36, 27, None], [None, 1, 4, 13])

    TESTS::

        sage: G = graphs.SymplecticDualPolarGraph(6,2); G                               # needs sage.libs.gap
        Symplectic Dual Polar Graph DSp(6, 2): Graph on 135 vertices
        sage: G.is_distance_regular(parameters=True)                                    # needs sage.libs.gap
        ([14, 12, 8, None], [None, 1, 3, 7])
        sage: graphs.SymplecticDualPolarGraph(6,6)                                      # needs sage.libs.gap
        Traceback (most recent call last):
        ...
        GAPError: Error, <subfield> must be a prime or a finite field
    """
    from sage.libs.gap.libgap import libgap
    G = _polar_graph(m, q, libgap.SymplecticGroup(m, q),
                     intersection_size=int((q**(m/2 - 1) - 1)/(q - 1)))

    G.relabel()
    G.name("Symplectic Dual Polar Graph DSp" + str((m, q)))
    if m == 4:
        G.name(G.name() + '; GQ' + str((q, q)))
    return G


def TaylorTwographDescendantSRG(q, clique_partition=False):
    r"""
    Return the descendant graph of the Taylor's two-graph for `U_3(q)`, `q` odd.

    This is a strongly regular graph with parameters
    `(v,k,\lambda,\mu)=(q^3, (q^2+1)(q-1)/2, (q-1)^3/4-1, (q^2+1)(q-1)/4)`
    obtained as a two-graph descendant of the
    :func:`Taylor's two-graph <sage.combinat.designs.twographs.taylor_twograph>` `T`.
    This graph admits a partition into cliques of size `q`, which are useful in
    :func:`~sage.graphs.graph_generators.GraphGenerators.TaylorTwographSRG`,
    a strongly regular graph on `q^3+1` vertices in the
    Seidel switching class of `T`, for which we need `(q^2+1)/2` cliques.
    The cliques are the `q^2` lines on `v_0` of the projective plane containing
    the unital for `U_3(q)`, and intersecting the unital (i.e. the vertices of
    the graph and the point we remove) in `q+1` points. This is all taken from
    §7E of [BL1984]_.

    INPUT:

    - ``q`` -- a power of an odd prime number

    - ``clique_partition`` -- boolean (default: ``False``); when set to
      ``True``, return `q^2-1` cliques of size `q` with empty pairwise
      intersection. (Removing all of them leaves a clique, too), and the point
      removed from the unital.

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: g = graphs.TaylorTwographDescendantSRG(3); g
        Taylor two-graph descendant SRG: Graph on 27 vertices
        sage: g.is_strongly_regular(parameters=True)
        (27, 10, 1, 5)
        sage: from sage.combinat.designs.twographs import taylor_twograph
        sage: T = taylor_twograph(3)                            # long time
        sage: g.is_isomorphic(T.descendant(T.ground_set()[1]))  # long time
        True
        sage: g = graphs.TaylorTwographDescendantSRG(5)         # not tested (long time)
        sage: g.is_strongly_regular(parameters=True)            # not tested (long time)
        (125, 52, 15, 26)

    TESTS::

        sage: # needs sage.rings.finite_rings
        sage: g,l,_ = graphs.TaylorTwographDescendantSRG(3, clique_partition=True)
        sage: all(g.is_clique(x) for x in l)
        True
        sage: graphs.TaylorTwographDescendantSRG(4)
        Traceback (most recent call last):
        ...
        ValueError: q must be an odd prime power
        sage: graphs.TaylorTwographDescendantSRG(6)
        Traceback (most recent call last):
        ...
        ValueError: q must be an odd prime power
    """
    p, k = is_prime_power(q, get_data=True)
    if not k or p == 2:
        raise ValueError('q must be an odd prime power')
    from sage.schemes.projective.projective_space import ProjectiveSpace
    from sage.rings.finite_rings.integer_mod import mod
    Fq = FiniteField(q**2, 'a')
    PG = map(tuple, ProjectiveSpace(2, Fq))

    def S(x, y):
        return sum(x[j] * y[2 - j] ** q for j in range(3))

    V = [x for x in PG if S(x, x) == 0]  # the points of the unital
    v0 = V[0]
    V.remove(v0)
    if mod(q, 4) == 1:
        G = Graph([V, lambda y, z: not (S(v0, y)*S(y, z)*S(z, v0)).is_square()], loops=False)
    else:
        G = Graph([V, lambda y, z: (S(v0, y)*S(y, z)*S(z, v0)).is_square()], loops=False)
    G.name("Taylor two-graph descendant SRG")
    if clique_partition:
        lines = [[t for t in V if t[0] + z * t[1] == 0]
                 for z in Fq if z]
        return (G, lines, v0)
    else:
        return G


def TaylorTwographSRG(q):
    r"""
    Return a strongly regular graph from the Taylor's two-graph for `U_3(q)`,
    `q` odd

    This is a strongly regular graph with parameters
    `(v,k,\lambda,\mu)=(q^3+1, q(q^2+1)/2, (q^2+3)(q-1)/4, (q^2+1)(q+1)/4)`
    in the Seidel switching class of
    :func:`Taylor two-graph <sage.combinat.designs.twographs.taylor_twograph>`.
    Details are in §7E of [BL1984]_.

    INPUT:

    - ``q`` -- a power of an odd prime number

    .. SEEALSO::

        * :meth:`~sage.graphs.graph_generators.GraphGenerators.TaylorTwographDescendantSRG`

    EXAMPLES::

        sage: t = graphs.TaylorTwographSRG(3); t                                        # needs sage.rings.finite_rings
        Taylor two-graph SRG: Graph on 28 vertices
        sage: t.is_strongly_regular(parameters=True)                                    # needs sage.rings.finite_rings
        (28, 15, 6, 10)
    """
    G, l, v0 = TaylorTwographDescendantSRG(q, clique_partition=True)
    G.add_vertex(v0)
    G.seidel_switching(sum(l[:(q**2 + 1)/2], []))
    G.name("Taylor two-graph SRG")
    return G


def AhrensSzekeresGeneralizedQuadrangleGraph(q, dual=False):
    r"""
    Return the collinearity graph of the generalized quadrangle `AS(q)`, or of
    its dual

    Let `q` be an odd prime power.  `AS(q)` is a generalized quadrangle
    (:wikipedia:`Generalized_quadrangle`) of
    order `(q-1,q+1)`, see 3.1.5 in [PT2009]_. Its points are elements
    of `F_q^3`, and lines are sets of size `q` of the form

    * `\{ (\sigma, a, b) \mid \sigma\in F_q \}`
    * `\{ (a, \sigma, b) \mid \sigma\in F_q \}`
    * `\{ (c \sigma^2 - b \sigma + a, -2 c \sigma + b, \sigma) \mid \sigma\in F_q \}`,

    where `a`, `b`, `c` are arbitrary elements of `F_q`.

    INPUT:

    - ``q`` -- a power of an odd prime number

    - ``dual`` -- boolean (default: ``False``); whether to return the
      collinearity graph of `AS(q)` or of the dual `AS(q)` (when ``True``)

    EXAMPLES::

        sage: g = graphs.AhrensSzekeresGeneralizedQuadrangleGraph(5); g
        AS(5); GQ(4, 6): Graph on 125 vertices
        sage: g.is_strongly_regular(parameters=True)
        (125, 28, 3, 7)
        sage: g = graphs.AhrensSzekeresGeneralizedQuadrangleGraph(5, dual=True); g
        AS(5)*; GQ(6, 4): Graph on 175 vertices
        sage: g.is_strongly_regular(parameters=True)
        (175, 30, 5, 5)
    """
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    p, k = is_prime_power(q, get_data=True)
    if not k or p == 2:
        raise ValueError('q must be an odd prime power')
    F = FiniteField(q, 'a')
    L = []
    for a in F:
        for b in F:
            L.append(tuple((s, a, b) for s in F))
            L.append(tuple((a, s, b) for s in F))
            for c in F:
                L.append(tuple((c*s**2 - b*s + a, -2*c*s + b, s) for s in F))
    if dual:
        G = IncidenceStructure(L).intersection_graph()
        G.name('AS(' + str(q) + ')*; GQ' + str((q + 1, q - 1)))
    else:
        G = IncidenceStructure(L).dual().intersection_graph()
        G.name('AS(' + str(q) + '); GQ' + str((q - 1, q + 1)))
    return G


def T2starGeneralizedQuadrangleGraph(q, dual=False, hyperoval=None, field=None, check_hyperoval=True):
    r"""
    Return the collinearity graph of the generalized quadrangle `T_2^*(q)`, or
    of its dual

    Let `q=2^k` and `\Theta=PG(3,q)`.  `T_2^*(q)` is a generalized quadrangle
    (:wikipedia:`Generalized_quadrangle`)
    of order `(q-1,q+1)`, see 3.1.3 in [PT2009]_. Fix a plane `\Pi \subset
    \Theta` and a
    `hyperoval <http://en.wikipedia.org/wiki/Oval_(projective_plane)#Even_q>`__
    `O \subset \Pi`. The points of `T_2^*(q):=T_2^*(O)` are the points of
    `\Theta` outside `\Pi`, and the lines are the lines of `\Theta` outside
    `\Pi` that meet `\Pi` in a point of `O`.

    INPUT:

    - ``q`` -- a power of two

    - ``dual`` -- boolean (default: ``False``); whether to return the graph of
      `T_2^*(O)` or of the dual `T_2^*(O)` (when ``True``)

    - ``hyperoval`` -- a hyperoval (i.e. a complete 2-arc; a set of points in
      the plane meeting every line in 0 or 2 points) in the plane of points with
      0th coordinate 0 in `PG(3,q)` over the field ``field``. Each point of
      ``hyperoval`` must be a length 4 vector over ``field`` with 1st non-0
      coordinate equal to 1. By default, ``hyperoval`` and ``field`` are not
      specified, and constructed on the fly. In particular, ``hyperoval`` we
      build is the classical one, i.e. a conic with the point of intersection of
      its tangent lines.

    - ``field`` -- an instance of a finite field of order `q`, must be provided
      if ``hyperoval`` is provided

    - ``check_hyperoval`` -- boolean (default: ``True``); whether to check
      ``hyperoval`` for correctness or not

    EXAMPLES:

    using the built-in construction::

        sage: # needs sage.combinat sage.rings.finite_rings
        sage: g = graphs.T2starGeneralizedQuadrangleGraph(4); g
        T2*(O,4); GQ(3, 5): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 18, 2, 6)
        sage: g = graphs.T2starGeneralizedQuadrangleGraph(4, dual=True); g
        T2*(O,4)*; GQ(5, 3): Graph on 96 vertices
        sage: g.is_strongly_regular(parameters=True)
        (96, 20, 4, 4)

    supplying your own hyperoval::

        sage: # needs sage.combinat sage.rings.finite_rings
        sage: F = GF(4,'b')
        sage: O = [vector(F,(0,0,0,1)),vector(F,(0,0,1,0))] + [vector(F, (0,1,x^2,x))
        ....:                                                  for x in F]
        sage: g = graphs.T2starGeneralizedQuadrangleGraph(4, hyperoval=O, field=F); g
        T2*(O,4); GQ(3, 5): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 18, 2, 6)

    TESTS::

        sage: # needs sage.combinat sage.rings.finite_rings
        sage: F = GF(4,'b')  # repeating a point...
        sage: O = [vector(F,(0,1,0,0)),vector(F,(0,0,1,0))]+[vector(F, (0,1,x^2,x)) for x in F]
        sage: graphs.T2starGeneralizedQuadrangleGraph(4, hyperoval=O, field=F)
        Traceback (most recent call last):
        ...
        RuntimeError: incorrect hyperoval size
        sage: O = [vector(F,(0,1,1,0)),vector(F,(0,0,1,0))]+[vector(F, (0,1,x^2,x)) for x in F]
        sage: graphs.T2starGeneralizedQuadrangleGraph(4, hyperoval=O, field=F)
        Traceback (most recent call last):
        ...
        RuntimeError: incorrect hyperoval
    """
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    from sage.combinat.designs.block_design import ProjectiveGeometryDesign as PG

    p, k = is_prime_power(q, get_data=True)
    if not k or p != 2:
        raise ValueError('q must be a power of 2')
    if field is None:
        F = FiniteField(q, 'a')
    else:
        F = field

    Theta = PG(3, 1, F, point_coordinates=1)
    Pi = set(x for x in Theta.ground_set() if x[0] == F.zero())
    if hyperoval is None:
        HO = set(x for x in Pi
                 if (x[1] + x[2] * x[3] == 0) or
                    (x[1] == 1 and x[2] == x[3] == 0))
    else:
        for v in hyperoval:
            v.set_immutable()

        HO = set(hyperoval)
        if check_hyperoval:
            if len(HO) != q + 2:
                raise RuntimeError("incorrect hyperoval size")
            for L in Theta.blocks():
                if set(L).issubset(Pi):
                    if len(HO.intersection(L)) not in [0, 2]:
                        raise RuntimeError("incorrect hyperoval")

    L = [[y for y in z if y not in HO]
         for z in [x for x in Theta.blocks() if len(HO.intersection(x)) == 1]]

    if dual:
        G = IncidenceStructure(L).intersection_graph()
        G.name('T2*(O,' + str(q) + ')*; GQ' + str((q + 1, q - 1)))
    else:
        G = IncidenceStructure(L).dual().intersection_graph()
        G.name('T2*(O,' + str(q) + '); GQ' + str((q - 1, q + 1)))
    return G


def HaemersGraph(q, hyperoval=None, hyperoval_matching=None, field=None, check_hyperoval=True):
    r"""
    Return the Haemers graph obtained from `T_2^*(q)^*`.

    Let `q` be a power of 2. In Sect. 8.A of [BL1984]_ one finds a construction
    of a strongly regular graph with parameters `(q^2(q+2),q^2+q-1,q-2,q)` from
    the graph of `T_2^*(q)^*`, constructed by
    :func:`~sage.graphs.graph_generators.GraphGenerators.T2starGeneralizedQuadrangleGraph`,
    by redefining adjacencies in the way specified by an arbitrary
    ``hyperoval_matching`` of the points (i.e. partitioning into size two parts)
    of ``hyperoval`` defining `T_2^*(q)^*`.

    While [BL1984]_ gives the construction in geometric terms, it can be
    formulated, and is implemented, in graph-theoretic ones, of re-adjusting the
    edges. Namely, `G=T_2^*(q)^*` has a partition into `q+2` independent sets
    `I_k` of size `q^2` each. Each vertex in `I_j` is adjacent to `q` vertices
    from `I_k`. Each `I_k` is paired to some `I_{k'}`, according to
    ``hyperoval_matching``. One adds edges `(s,t)` for `s,t \in I_k` whenever
    `s` and `t` are adjacent to some `u \in I_{k'}`, and removes all the edges
    between `I_k` and `I_{k'}`.

    INPUT:

    - ``q`` -- a power of two

    - ``hyperoval_matching`` -- if ``None`` (default), pair each `i`-th point of
      ``hyperoval`` with `(i+1)`-th. Otherwise, specifies the pairing
      in the format `((i_1,i'_1),(i_2,i'_2),...)`.

    - ``hyperoval`` -- a hyperoval defining `T_2^*(q)^*`. If ``None`` (default),
      the classical hyperoval obtained from a conic is used. See the
      documentation of
      :func:`~sage.graphs.graph_generators.GraphGenerators.T2starGeneralizedQuadrangleGraph`,
      for more information.

    - ``field`` -- an instance of a finite field of order `q`, must be provided
      if ``hyperoval`` is provided

    - ``check_hyperoval`` -- boolean (default: ``True``); whether to check
      ``hyperoval`` for correctness or not

    EXAMPLES:

    using the built-in constructions::

        sage: # needs sage.combinat sage.rings.finite_rings
        sage: g = graphs.HaemersGraph(4); g
        Haemers(4): Graph on 96 vertices
        sage: g.is_strongly_regular(parameters=True)
        (96, 19, 2, 4)

    supplying your own hyperoval_matching::

        sage: # needs sage.combinat sage.rings.finite_rings
        sage: g = graphs.HaemersGraph(4, hyperoval_matching=((0,5),(1,4),(2,3))); g
        Haemers(4): Graph on 96 vertices
        sage: g.is_strongly_regular(parameters=True)
        (96, 19, 2, 4)

    TESTS::

        sage: # needs sage.combinat sage.rings.finite_rings
        sage: F = GF(4,'b')  # repeating a point...
        sage: O = [vector(F,(0,1,0,0)),vector(F,(0,0,1,0))]+[vector(F, (0,1,x^2,x)) for x in F]
        sage: graphs.HaemersGraph(4, hyperoval=O, field=F)
        Traceback (most recent call last):
        ...
        RuntimeError: incorrect hyperoval size
        sage: O = [vector(F,(0,1,1,0)),vector(F,(0,0,1,0))]+[vector(F, (0,1,x^2,x)) for x in F]
        sage: graphs.HaemersGraph(4, hyperoval=O, field=F)
        Traceback (most recent call last):
        ...
        RuntimeError: incorrect hyperoval

        sage: g = graphs.HaemersGraph(8); g             # not tested (long time)        # needs sage.rings.finite_rings
        Haemers(8): Graph on 640 vertices
        sage: g.is_strongly_regular(parameters=True)    # not tested (long time)        # needs sage.rings.finite_rings
        (640, 71, 6, 8)
    """
    from sage.modules.free_module_element import free_module_element as vector
    from sage.rings.finite_rings.finite_field_constructor import GF
    from itertools import combinations

    p, k = is_prime_power(q, get_data=True)
    if not k or p != 2:
        raise ValueError('q must be a power of 2')

    if hyperoval_matching is None:
        hyperoval_matching = [(2 * K + 1, 2 * K) for K in range(1 + q // 2)]
    if field is None:
        F = GF(q, 'a')
    else:
        F = field

    # for q=8, 95% of CPU time taken by this function is spent in the following call
    G = T2starGeneralizedQuadrangleGraph(q, field=F, dual=True,
                                         hyperoval=hyperoval,
                                         check_hyperoval=check_hyperoval)

    def normalize(v):  # make sure the 1st non-0 coordinate is 1.
        d = next(x for x in v if x != F.zero())
        return vector([x / d for x in v])

    # build the partition into independent sets
    P = [tuple(normalize(v[0] - v[1])) for v in G.vertices(sort=True)]
    Pi_to_int = {Pi: i for i, Pi in enumerate(set(P))}
    I_ks = {x: [] for x in range(q + 2)}  # the partition into I_k's
    for i, Pi in enumerate(P):
        I_ks[Pi_to_int[tuple(Pi)]].append(i)

    # perform the adjustment of the edges, as described.
    G.relabel(range(G.order()))
    cliques = []
    for i, j in hyperoval_matching:
        Pij = set(I_ks[i] + I_ks[j])
        for v in Pij:
            cliques.append(Pij.intersection(G.neighbors(v)))
        G.delete_edges(G.edge_boundary(I_ks[i], I_ks[j]))  # edges on (I_i,I_j)
    G.add_edges(e for c in cliques for e in combinations(c, 2))
    G.name('Haemers(' + str(q) + ')')
    return G


def CossidentePenttilaGraph(q):
    r"""
    Return the Cossidente-Penttila
    `((q^3+1)(q+1)/2,(q^2+1)(q-1)/2,(q-3)/2,(q-1)^2/2)`-strongly regular graph

    For each odd prime power `q`, one can partition the points of the
    `O_6^-(q)`-generalized quadrangle `GQ(q,q^2)` into two parts, so that on any
    of them the induced subgraph of the point graph of the GQ has parameters as
    above [CP2005]_.

    Directly following the construction in [CP2005]_ is not efficient, as one
    then needs to construct the dual `GQ(q^2,q)`. Thus we describe here a more
    efficient approach that we came up with, following a suggestion by
    T.Penttila. Namely, this partition is invariant under the subgroup
    `H=\Omega_3(q^2)<O_6^-(q)`. We build the appropriate `H`, which leaves the
    form `B(X,Y,Z)=XY+Z^2` invariant, and pick up two orbits of `H` on the
    `F_q`-points. One them is `B`-isotropic, and we take the representative
    `(1:0:0)`. The other one corresponds to the points of `PG(2,q^2)` that have
    all the lines on them either missing the conic specified by `B`, or
    intersecting the conic in two points. We take `(1:1:e)` as the
    representative. It suffices to pick `e` so that `e^2+1` is not a square in
    `F_{q^2}`. Indeed, The conic can be viewed as the union of `\{(0:1:0)\}` and
    `\{(1:-t^2:t) | t \in F_{q^2}\}`.  The coefficients of a generic line on
    `(1:1:e)` are `[1:-1-eb:b]`, for `-1\neq eb`.  Thus, to make sure the
    intersection with the conic is always even, we need that the discriminant of
    `1+(1+eb)t^2+tb=0` never vanishes, and this is if and only if `e^2+1` is not
    a square. Further, we need to adjust `B`, by multiplying it by appropriately
    chosen `\nu`, so that `(1:1:e)` becomes isotropic under the relative trace
    norm `\nu B(X,Y,Z)+(\nu B(X,Y,Z))^q`. The latter is used then to define the
    graph.

    INPUT:

    - ``q`` -- an odd prime power

    EXAMPLES:

    For `q=3` one gets Sims-Gewirtz graph. ::

        sage: G = graphs.CossidentePenttilaGraph(3)     # optional - gap_package_grape
        sage: G.is_strongly_regular(parameters=True)    # optional - gap_package_grape
        (56, 10, 0, 2)

    For `q>3` one gets new graphs. ::

        sage: G = graphs.CossidentePenttilaGraph(5)     # optional - gap_package_grape
        sage: G.is_strongly_regular(parameters=True)    # optional - gap_package_grape
        (378, 52, 1, 8)

    TESTS::

        sage: G = graphs.CossidentePenttilaGraph(7)     # optional - gap_package_grape, long time
        sage: G.is_strongly_regular(parameters=True)    # optional - gap_package_grape, long time
        (1376, 150, 2, 18)
        sage: graphs.CossidentePenttilaGraph(2)
        Traceback (most recent call last):
        ...
        ValueError: q(=2) must be an odd prime power
    """
    p, k = is_prime_power(q, get_data=True)
    if not k or p == 2:
        raise ValueError('q(={}) must be an odd prime power'.format(q))

    from sage.features.gap import GapPackage
    GapPackage("grape", spkg='gap_packages').require()

    from sage.libs.gap.libgap import libgap
    adj_list = libgap.function_factory("""function(q)
        local z, e, so, G, nu, G1, G0, B, T, s, O1, O2, x, sqo;
        LoadPackage("grape");
        G0:=SO(3,q^2);
        so:=GeneratorsOfGroup(G0);
        G1:=Group(Comm(so[1],so[2]),Comm(so[1],so[3]),Comm(so[2],so[3]));
        B:=InvariantBilinearForm(G0).matrix;
        z:=Z(q^2); e:=z; sqo:=(q^2-1)/2;
        if IsInt(sqo/Order(e^2+z^0)) then
            e:=z^First([2..q^2-2], x-> not IsInt(sqo/Order(z^(2*x)+z^0)));
        fi;
        nu:=z^First([0..q^2-2], x->z^x*(e^2+z^0)+(z^x*(e^2+z^0))^q=0*z);
        T:=function(x)
            local r;
            r:=nu*x*B*x;
            return r+r^q;
        end;
        s:=Group([Z(q)*IdentityMat(3,GF(q))]);
        O1:=Orbit(G1, Set(Orbit(s,z^0*[1,0,0])), OnSets);
        O2:=Orbit(G1, Set(Orbit(s,z^0*[1,1,e])), OnSets);
        G:=Graph(G1,Concatenation(O1,O2),OnSets,
            function(x,y) return x<>y and 0*z=T(x[1]+y[1]); end);
        return List([1..OrderGraph(G)],x->Adjacency(G,x));
        end;""")

    adj = adj_list(q)  # for each vertex, we get the list of vertices it is adjacent to
    G = Graph(((i, int(j - 1))
               for i, ni in enumerate(adj) for j in ni),
              format='list_of_edges', multiedges=False)
    G.name('CossidentePenttila(' + str(q) + ')')
    return G


def Nowhere0WordsTwoWeightCodeGraph(q, hyperoval=None, field=None, check_hyperoval=True):
    r"""
    Return the subgraph of nowhere 0 words from two-weight code of projective
    plane hyperoval.

    Let `q=2^k` and `\Pi=PG(2,q)`.  Fix a
    `hyperoval <http://en.wikipedia.org/wiki/Oval_(projective_plane)#Even_q>`__
    `O \subset \Pi`. Let `V=F_q^3` and `C` the two-weight 3-dimensional linear
    code over `F_q` with words `c(v)` obtained from `v\in V` by computing

    .. MATH::

        c(v)=(\langle v,o_1 \rangle,...,\langle v,o_{q+2} \rangle), o_j \in O.

    `C` contains `q(q-1)^2/2` words without 0 entries. The subgraph of the
    strongly regular graph of `C` induced on the latter words is also strongly
    regular, assuming `q>4`. This is a construction due to A.E.Brouwer
    [Bro2016]_, and leads to graphs with parameters also given by a construction
    in [HHL2009]_.  According to [Bro2016]_, these two constructions are likely
    to produce isomorphic graphs.

    INPUT:

    - ``q`` -- a power of two

    - ``hyperoval`` -- a hyperoval (i.e. a complete 2-arc; a set of points in
      the plane meeting every line in 0 or 2 points) in `PG(2,q)` over the field
      ``field``.  Each point of ``hyperoval`` must be a length 3 vector over
      ``field`` with 1st non-0 coordinate equal to 1. By default, ``hyperoval``
      and ``field`` are not specified, and constructed on the fly. In
      particular, ``hyperoval`` we build is the classical one, i.e. a conic with
      the point of intersection of its tangent lines.

    - ``field`` -- an instance of a finite field of order `q`; must be provided
      if ``hyperoval`` is provided

    - ``check_hyperoval`` -- boolean (default: ``True``); whether to check
      ``hyperoval`` for correctness or not

    .. SEEALSO::

        - :func:`~sage.graphs.strongly_regular_db.is_nowhere0_twoweight`

    EXAMPLES:

    using the built-in construction::

        sage: # needs sage.combinat sage.rings.finite_rings
        sage: g = graphs.Nowhere0WordsTwoWeightCodeGraph(8); g
        Nowhere0WordsTwoWeightCodeGraph(8): Graph on 196 vertices
        sage: g.is_strongly_regular(parameters=True)
        (196, 60, 14, 20)
        sage: g = graphs.Nowhere0WordsTwoWeightCodeGraph(16)  # not tested (long time)
        sage: g.is_strongly_regular(parameters=True)          # not tested (long time)
        (1800, 728, 268, 312)

    supplying your own hyperoval::

        sage: # needs sage.combinat sage.rings.finite_rings
        sage: F = GF(8)
        sage: O = [vector(F,(0,0,1)),vector(F,(0,1,0))] + [vector(F, (1,x^2,x))
        ....:                                              for x in F]
        sage: g = graphs.Nowhere0WordsTwoWeightCodeGraph(8,hyperoval=O,field=F); g
        Nowhere0WordsTwoWeightCodeGraph(8): Graph on 196 vertices
        sage: g.is_strongly_regular(parameters=True)
        (196, 60, 14, 20)

    TESTS::

        sage: # needs sage.combinat sage.rings.finite_rings
        sage: F = GF(8)  # repeating a point...
        sage: O = [vector(F,(1,0,0)),vector(F,(0,1,0))]+[vector(F, (1,x^2,x)) for x in F]
        sage: graphs.Nowhere0WordsTwoWeightCodeGraph(8,hyperoval=O,field=F)
        Traceback (most recent call last):
        ...
        RuntimeError: incorrect hyperoval size
        sage: O = [vector(F,(1,1,0)),vector(F,(0,1,0))]+[vector(F, (1,x^2,x)) for x in F]
        sage: graphs.Nowhere0WordsTwoWeightCodeGraph(8,hyperoval=O,field=F)
        Traceback (most recent call last):
        ...
        RuntimeError: incorrect hyperoval
    """
    from sage.combinat.designs.block_design import ProjectiveGeometryDesign as PG
    from sage.matrix.constructor import matrix

    p, k = is_prime_power(q, get_data=True)
    if not k or p != 2:
        raise ValueError('q must be a power of 2')
    if k < 3:
        raise ValueError('q must be a at least 8')
    if field is None:
        F = FiniteField(q, 'a')
    else:
        F = field

    Theta = PG(2, 1, F, point_coordinates=1)
    Pi = Theta.ground_set()
    if hyperoval is None:
        hyperoval = [x for x in Pi
                     if (x[0] + x[1] * x[2] == 0) or
                        (x[0] == 1 and x[1] == x[2] == 0)]
    else:
        for v in hyperoval:
            v.set_immutable()

        if check_hyperoval:
            HO = set(hyperoval)
            if len(HO) != q + 2:
                raise RuntimeError("incorrect hyperoval size")
            for L in Theta.blocks():
                if set(L).issubset(Pi):
                    if len(HO.intersection(L)) not in [0, 2]:
                        raise RuntimeError("incorrect hyperoval")
    M = matrix(hyperoval)
    F_0 = F.zero()
    C = [p for p in [M*x for x in F**3] if F_0 not in p]

    for x in C:
        x.set_immutable()
    G = Graph([C, lambda x, y: F.zero() not in x + y])
    G.name('Nowhere0WordsTwoWeightCodeGraph(' + str(q) + ')')
    G.relabel()
    return G


def OrthogonalDualPolarGraph(e, d, q):
    r"""
    Return the dual polar graph on `GO^e(n,q)` of diameter `d`.

    The value of `n` is determined by `d` and `e`.

    The graph is distance-regular with classical parameters `(d, q, 0, q^e)`.

    INPUT:

    - ``e`` -- integer; type of the orthogonal polar space to consider;
      must be `-1, 0` or  `1`

    - ``d`` -- integer; diameter of the graph

    - ``q`` -- integer; prime power; order of the finite field over which to
      build the polar space

    EXAMPLES::

        sage: # needs sage.libs.gap
        sage: G = graphs.OrthogonalDualPolarGraph(1,3,2)
        sage: G.is_distance_regular(True)
        ([7, 6, 4, None], [None, 1, 3, 7])
        sage: G = graphs.OrthogonalDualPolarGraph(0,3,3)        # long time
        sage: G.is_distance_regular(True)                       # long time
        ([39, 36, 27, None], [None, 1, 4, 13])
        sage: G.order()                                         # long time
        1120

    REFERENCES:

    See [BCN1989]_ pp. 274-279 or [VDKT2016]_ p. 22.

    TESTS::

        sage: # needs sage.libs.gap
        sage: G = graphs.OrthogonalDualPolarGraph(0,3,2)
        sage: G.is_distance_regular(True)
        ([14, 12, 8, None], [None, 1, 3, 7])
        sage: G = graphs.OrthogonalDualPolarGraph(-1,3,2)       # long time
        sage: G.is_distance_regular(True)                       # long time
        ([28, 24, 16, None], [None, 1, 3, 7])
        sage: G = graphs.OrthogonalDualPolarGraph(1,3,4)
        sage: G.is_distance_regular(True)
        ([21, 20, 16, None], [None, 1, 5, 21])
        sage: G = graphs.OrthogonalDualPolarGraph(1,4,2)
        sage: G.is_distance_regular(True)
        ([15, 14, 12, 8, None], [None, 1, 3, 7, 15])
    """
    from sage.libs.gap.libgap import libgap
    from sage.matrix.constructor import Matrix
    from sage.modules.free_module import VectorSpace
    from sage.rings.finite_rings.finite_field_constructor import GF
    import itertools

    def hashable(v):
        v.set_immutable()
        return v

    if e not in {0, 1, -1}:
        raise ValueError("e must by 0, +1 or -1")

    m = 2*d + 1 - e

    group = libgap.GeneralOrthogonalGroup(e, m, q)
    M = Matrix(libgap.InvariantQuadraticForm(group)["matrix"])
    # Q(x) = xMx is our quadratic form

    # we need to find a totally isotropic subspace of dimension d
    K = M.kernel()
    isotropicBasis = list(K.basis())

    # extend K to a maximal isotropic subspace
    if K.dimension() < d:
        V = VectorSpace(GF(q), m)

        # get all projective points not in K
        candidates = set(map(hashable, [P.basis()[0] for P in V.subspaces(1)]))
        hashableK = map(hashable, [P.basis()[0] for P in K.subspaces(1)])
        candidates = candidates.difference(hashableK)

        nonZeroScalars = [x for x in GF(q) if not x.is_zero()]
        while K.dimension() < d:
            found = False
            while not found:
                v = candidates.pop()
                if v*M*v == 0:
                    # found another isotropic point
                    # check if we can add it to K
                    found = True
                    for w in isotropicBasis:
                        if w*M*v + v*M*w != 0:
                            found = False
                            break
            # here we found a valid point
            isotropicBasis.append(v)

            # remove new points of K
            newVectors = map(hashable,
                             [k + s*v for k in K for s in nonZeroScalars])
            candidates.difference(newVectors)
            K = V.span(isotropicBasis)

        # here K is a totally isotropic subspace of dimension d
        isotropicBasis = list(K.basis())

    W = libgap.FullRowSpace(libgap.GF(q), m)
    isoS = libgap.Subspace(W, isotropicBasis)  # gap version of K

    allIsoPoints = libgap.Orbit(group, isotropicBasis[0], libgap.OnLines)
    permutation = libgap.Action(group, allIsoPoints, libgap.OnLines)
    # this is the permutation group generated by GO^e(n,q) acting on
    # projective isotropic points

    # convert isoS into a list of ints representing the projective points
    isoSPoints = [libgap.Elements(libgap.Basis(x))[0]
                  for x in libgap.Elements(isoS.Subspaces(1))]
    isoSPointsInt = libgap.Set([libgap.Position(allIsoPoints, x)
                                for x in isoSPoints])

    # all isotropic subspaces of dimension d
    allIsoSubspaces = libgap.Orbit(permutation, isoSPointsInt, libgap.OnSets)

    # number of projective points in a (d-1)-subspace
    intersection_size = (q**(d-1) - 1) // (q-1)

    edges = []
    n = len(allIsoSubspaces)
    for i, j in itertools.combinations(range(n), 2):
        if libgap.Size(libgap.Intersection(allIsoSubspaces[i],
                                           allIsoSubspaces[j])) \
                                           == intersection_size:
            edges.append((i, j))

    G = Graph(edges, format='list_of_edges')
    G.name("Dual Polar Graph on Orthogonal group (%d, %d, %d)" % (e, m, q))
    return G
