r"""
Some useful functions for the matroid class.

For direct access to the methods :meth:`newlabel`, :meth:`setprint` and
:meth:`get_nonisomorphic_matroids`, type::

    sage: from sage.matroids.advanced import *

See also :mod:`sage.matroids.advanced`.

AUTHORS:

- Stefan van Zwam (2011-06-24): initial version
"""
# ****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections.abc import Iterable
from sage.matrix.constructor import Matrix
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.all import SageObject
from operator import itemgetter


def setprint(X):
    """
    Print nested data structures nicely.

    Python's data structures ``set`` and ``frozenset`` do not print nicely.
    This function can be used as replacement for ``print`` to overcome this.
    For direct access to ``setprint``, run::

        sage: from sage.matroids.advanced import *

    .. NOTE::

        This function will be redundant when Sage moves to Python 3, since the
        default ``print`` will suffice then.

    INPUT:

    - ``X`` -- any Python object

    OUTPUT: none; however, the function prints a nice representation of ``X``

    EXAMPLES:

    Output looks much better::

        sage: from sage.matroids.advanced import setprint
        sage: L = [{1, 2, 3}, {1, 2, 4}, {2, 3, 4}, {4, 1, 3}]
        sage: print(L)
        [{1, 2, 3}, {1, 2, 4}, {2, 3, 4}, {1, 3, 4}]
        sage: setprint(L)
        [{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}]

    Note that for iterables, the effect can be undesirable::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.Fano().delete('efg')
        sage: M.bases()
        SetSystem of 3 sets over 4 elements
        sage: setprint(M.bases())
        [{'a', 'b', 'c'}, {'a', 'b', 'd'}, {'a', 'c', 'd'}]

    An exception was made for subclasses of SageObject::

        sage: from sage.matroids.advanced import setprint
        sage: G = graphs.PetersenGraph()                                                # needs sage.graphs
        sage: list(G)                                                                   # needs sage.graphs
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: setprint(G)                                                               # needs sage.graphs
        Petersen graph: Graph on 10 vertices
    """
    print(setprint_s(X, toplevel=True))


def setprint_s(X, toplevel=False):
    """
    Create the string for use by ``setprint()``.

    INPUT:

    - ``X`` -- any Python object
    - ``toplevel`` -- boolean (default: ``False``); indicates whether this is a
      recursion or not

    OUTPUT:

    A string representation of the object, with nice notation for sets and
    frozensets.

    EXAMPLES::

        sage: from sage.matroids.utilities import setprint_s
        sage: L = [{1, 2, 3}, {1, 2, 4}, {2, 3, 4}, {4, 1, 3}]
        sage: setprint_s(L)
        '[{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}]'

    The ``toplevel`` argument only affects strings, to mimic ``print``'s
    behavior::

        sage: X = 'abcd'
        sage: setprint_s(X)
        "'abcd'"
        sage: setprint_s(X, toplevel=True)
        'abcd'
    """
    if isinstance(X, (frozenset, set)):
        return '{' + ', '.join(sorted(setprint_s(x) for x in X)) + '}'
    elif isinstance(X, dict):
        return '{' + ', '.join(sorted(setprint_s(key) + ': ' + setprint_s(val)
                                      for key, val in X.items())) + '}'
    elif isinstance(X, str):
        if toplevel:
            return X
        else:
            return "'" + X + "'"
    elif isinstance(X, Iterable) and not isinstance(X, SageObject):
        return '[' + ', '.join(sorted(setprint_s(x) for x in X)) + ']'
    else:
        return repr(X)


def newlabel(groundset):
    r"""
    Create a new element label different from the labels in ``groundset``.

    INPUT:

    - ``groundset`` -- set of objects

    OUTPUT: string not in the set ``groundset``

    For direct access to ``newlabel``, run::

        sage: from sage.matroids.advanced import *

    ALGORITHM:

    #. Create a set of all one-character alphanumeric strings.
    #. Remove the string representation of each existing element from this
       list.
    #. If the list is nonempty, return any element.
    #. Otherwise, return the concatenation of the strings of each existing
       element, preceded by 'e'.

    EXAMPLES::

        sage: from sage.matroids.advanced import newlabel
        sage: S = set(['a', 42, 'b'])
        sage: newlabel(S) in S
        False

        sage: S = set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789')
        sage: t = newlabel(S)
        sage: len(t)
        63
        sage: t[0]
        'e'
    """
    char_list = set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789')
    char_list.difference_update([str(e) for e in groundset])
    try:
        s = char_list.pop()
    except KeyError:
        s = 'e' + ''.join([str(e) for e in groundset])
    return s


def sanitize_contractions_deletions(matroid, contractions, deletions):
    r"""
    Return a fixed version of sets ``contractions`` and ``deletions``.

    INPUT:

    - ``matroid`` -- a :class:`Matroid <sage.matroids.matroid.Matroid>`
      instance
    - ``contractions`` -- a subset of the groundset
    - ``deletions`` -- a subset of the groundset

    OUTPUT: an independent set ``C`` and a coindependent set ``D`` such that

        ``matroid / contractions \ deletions == matroid / C \ D``

    Raise an error if either is not a subset of the groundset of ``matroid``
    or if they are not disjoint.

    This function is used by the
    :meth:`Matroid.minor() <sage.matroids.matroid.Matroid.minor>` method.

    EXAMPLES::

        sage: from sage.matroids.utilities import setprint
        sage: from sage.matroids.utilities import sanitize_contractions_deletions
        sage: M = matroids.catalog.Fano()
        sage: setprint(sanitize_contractions_deletions(M, 'abc', 'defg'))
        [{'a', 'b', 'c'}, {'d', 'e', 'f', 'g'}]
        sage: setprint(sanitize_contractions_deletions(M, 'defg', 'abc'))
        [{'a', 'b', 'c', 'f'}, {'d', 'e', 'g'}]
        sage: setprint(sanitize_contractions_deletions(M, [1, 2, 3], 'efg'))
        Traceback (most recent call last):
        ...
        ValueError: [1, 2, 3] is not a subset of the groundset
        sage: setprint(sanitize_contractions_deletions(M, 'efg', [1, 2, 3]))
        Traceback (most recent call last):
        ...
        ValueError: [1, 2, 3] is not a subset of the groundset
        sage: setprint(sanitize_contractions_deletions(M, 'ade', 'efg'))
        Traceback (most recent call last):
        ...
        ValueError: contraction and deletion sets are not disjoint.
    """
    if not contractions:
        contractions = frozenset()
    else:
        contractions = matroid._subset(contractions)

    if not deletions:
        deletions = frozenset()
    else:
        deletions = matroid._subset(deletions)

    if not contractions.isdisjoint(deletions):
        raise ValueError("contraction and deletion sets are not disjoint.")

    conset = matroid._max_independent(contractions)
    delset = matroid._max_coindependent(deletions)

    return conset.union(deletions.difference(delset)), delset.union(contractions.difference(conset))


def make_regular_matroid_from_matroid(matroid):
    r"""
    Attempt to construct a regular representation of a matroid.

    INPUT:

    - ``matroid`` -- matroid

    OUTPUT:

    Return a `(0, 1, -1)`-matrix over the integers such that, if the input is
    a regular matroid, then the output is a totally unimodular matrix
    representing that matroid.

    EXAMPLES::

        sage: from sage.matroids.utilities import make_regular_matroid_from_matroid
        sage: make_regular_matroid_from_matroid(                                        # needs sage.graphs
        ....:               matroids.CompleteGraphic(6)).is_isomorphic(
        ....:                                     matroids.CompleteGraphic(6))
        True
    """
    import sage.matroids.linear_matroid
    M = matroid
    if isinstance(M, sage.matroids.linear_matroid.RegularMatroid):
        return M

    from sage.graphs.bipartite_graph import BipartiteGraph

    rk = M.full_rank()
    # First create a reduced 0-1 matrix
    B = list(M.basis())
    NB = list(M.groundset().difference(B))
    dB = {}
    i = 0
    for e in B:
        dB[e] = i
        i += 1
    dNB = {}
    i = 0
    for e in NB:
        dNB[e] = i
        i += 1
    A = Matrix(ZZ, len(B), len(NB), 0)
    G = BipartiteGraph(A.transpose())  # Sage's BipartiteGraph uses the column set as first color class. This is an edgeless graph.
    for e in NB:
        C = M.circuit(B + [e])
        for f in C.difference([e]):
            A[dB[f], dNB[e]] = 1
    # Change some entries from -1 to 1
    entries = list(BipartiteGraph(A.transpose()).edges(labels=False, sort=False))
    while entries:
        L = [G.shortest_path(u, v) for u, v in entries]
        mindex, minval = min(enumerate(L), key=lambda x: len(x[1]))

        # if minval = 0, there is an edge not spanned by the current subgraph. Its entry is free to be scaled any way.
        if len(minval) > 0:  # DUBIOUS !!
            # Check the subdeterminant
            S = frozenset(L[mindex])
            rows = []
            cols = []
            for i in S:
                if i < rk:
                    rows.append(i)
                else:
                    cols.append(i - rk)
            if A[rows, cols].det() != 0:
                A[entries[mindex][0], entries[mindex][1] - rk] = -1
        G.add_edge(entries[mindex][0], entries[mindex][1])
        entries.pop(mindex)
    return sage.matroids.linear_matroid.RegularMatroid(groundset=B + NB, reduced_matrix=A)


def get_nonisomorphic_matroids(MSet):
    """
    Return non-isomorphic members of the matroids in set ``MSet``.

    For direct access to ``get_nonisomorphic_matroids``, run::

        sage: from sage.matroids.advanced import *

    INPUT:

    - ``MSet`` -- an iterable whose members are matroids

    OUTPUT:

    A list containing one representative of each isomorphism class of
    members of ``MSet``.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: L = matroids.Uniform(3, 5).extensions()
        sage: len(list(L))
        32
        sage: len(get_nonisomorphic_matroids(L))
        5
    """
    OutSet = []
    for M in MSet:
        seen = False
        for N in OutSet:
            if N.is_isomorphic(M):
                seen = True
                break
        if not seen:
            OutSet.append(M)
    return OutSet


def spanning_forest(M):
    r"""
    Return a list of edges of a spanning forest of the bipartite
    graph defined by `M`

    INPUT:

    - ``M`` -- a matrix defining a bipartite graph G. The vertices are the
      rows and columns, if `M[i,j]` is nonzero, then there is an edge
      between row `i` and column `j`.

    OUTPUT:

    A list of tuples `(r_i,c_i)` representing edges between row `r_i` and
    column `c_i`.

    EXAMPLES::

        sage: from sage.matroids.utilities import spanning_forest
        sage: len(spanning_forest(matrix([[1,1,1],[1,1,1],[1,1,1]])))                   # needs sage.graphs
        5
        sage: len(spanning_forest(matrix([[0,0,1],[0,1,0],[0,1,0]])))                   # needs sage.graphs
        3
    """
    from sage.graphs.graph import Graph
    from sage.graphs.spanning_tree import kruskal

    # Given a matrix, produce a spanning tree
    G = Graph()
    m = M.ncols()
    for (x, y) in M.dict():
        G.add_edge(x + m, y)
    T = []
    # find spanning tree in each component
    for component in G.connected_components_subgraphs():
        spanning_tree = kruskal(component)
        for (x, y, z) in spanning_tree:
            if x < m:
                t = x
                x = y
                y = t
            T.append((x - m, y))
    return T


def spanning_stars(M):
    r"""
    Return the edges of a connected subgraph that is a union of
    all edges incident some subset of vertices.

    INPUT:

    - ``M`` -- a matrix defining a bipartite graph G. The vertices are the
      rows and columns, if `M[i,j]` is nonzero, then there is an edge
      between row i and column 0.

    OUTPUT:

    A list of tuples `(row,column)` in a spanning forest of the bipartite graph
    defined by ``M``.

    EXAMPLES::

        sage: from sage.matroids.utilities import spanning_stars
        sage: edges = spanning_stars(matrix([[1,1,1],[1,1,1],[1,1,1]]))                 # needs sage.graphs
        sage: Graph([(x+3, y) for x,y in edges]).is_connected()                         # needs sage.graphs
        True
    """
    from sage.graphs.graph import Graph

    G = Graph()
    m = M.ncols()
    for x, y in M.dict():
        G.add_edge(x + m, y)

    delta = (M.nrows() + m)**0.5
    # remove low degree vertices
    H = []
    # candidate vertices
    V_0 = set([])
    d = 0
    while G.order():
        x, d = min(G.degree_iterator(labels=True), key=itemgetter(1))
        if d < delta:
            V_0.add(x)
            H.extend(G.edges_incident(x, False))
            G.delete_vertex(x)
        else:
            break

    # min degree is at least sqrt(n)
    # greedily remove vertices
    G2 = G.copy()
    # set of picked vertices
    V_1 = set([])
    while G2.order():
        # choose vertex with maximum degree in G2
        x, d = max(G2.degree_iterator(labels=True), key=itemgetter(1))
        V_1.add(x)
        G2.delete_vertices(G2.neighbors(x))
        G2.delete_vertex(x)

    # G2 is a graph of all edges incident to V_1
    G2 = Graph()
    for v in V_1:
        for u in G.neighbors(v):
            G2.add_edge(u, v)

    V = V_0 | V_1
    # compute a spanning tree
    T = spanning_forest(M)
    for x, y in T:
        if x not in V and y not in V:
            V.add(v)

    for v in V:
        if G.has_vertex(v):  # some vertices are not in G
            H.extend(G.edges_incident(v, False))

    # T contain all edges in some spanning tree
    T = []
    for x, y in H:
        if x < m:
            t = x
            x = y
            y = t
        T.append((x - m, y))
    return T

# Partial fields and lifting


def lift_cross_ratios(A, lift_map=None):
    r"""
    Return a matrix which arises from the given matrix by lifting cross ratios.

    INPUT:

    - ``A`` -- a matrix over a ring ``source_ring``
    - ``lift_map`` -- a Python dictionary, mapping each cross ratio of ``A`` to
      some element of a target ring, and such that
      ``lift_map[source_ring(1)] = target_ring(1)``

    OUTPUT: ``Z`` -- a matrix over the ring ``target_ring``

    The intended use of this method is to create a (reduced) matrix representation of a
    matroid ``M`` over a ring ``target_ring``, given a (reduced) matrix representation of
    ``A`` of ``M`` over a ring ``source_ring`` and a map ``lift_map`` from ``source_ring``
    to ``target_ring``.

    This method will create a unique candidate representation ``Z``, but will not verify
    if ``Z`` is indeed a representation of ``M``. However, this is guaranteed if the
    conditions of the lift theorem (see [PvZ2010]_) hold for the lift map in combination with
    the matrix ``A``.

    For a lift map `f` and a matrix `A` these conditions are as follows. First of all
    `f: S \rightarrow T`, where `S` is a set of invertible elements of the source ring and
    `T` is a set of invertible elements of the target ring. The matrix `A` has entries
    from the source ring, and each cross ratio of `A` is contained in `S`. Moreover:

    - `1 \in S`, `1 \in T`;
    - for all `x \in S`: `f(x) = 1` if and only if `x = 1`;
    - for all `x, y \in S`: if `x + y = 0` then `f(x) + f(y) = 0`;
    - for all `x, y \in S`: if `x + y = 1` then `f(x) + f(y) = 1`;
    - for all `x, y, z \in S`: if  `xy = z` then `f(x)f(y) = f(z)`.

    Any ring homomorphism `h: P \rightarrow R` induces a lift map from the set of units `S` of
    `P` to the set of units `T` of `R`. There exist lift maps which do not arise in
    this manner. Several such maps can be created by the function
    :meth:`lift_map() <sage.matroids.utilities.lift_map>`.

    .. SEEALSO::

        :meth:`lift_map() <sage.matroids.utilities.lift_map>`

    EXAMPLES::

        sage: # needs sage.graphs
        sage: from sage.matroids.advanced import lift_cross_ratios, lift_map, LinearMatroid
        sage: R = GF(7)
        sage: to_sixth_root_of_unity = lift_map('sru')                                  # needs sage.rings.number_field
        sage: A = Matrix(R, [[1, 0, 6, 1, 2],[6, 1, 0, 0, 1],[0, 6, 3, 6, 0]])
        sage: A
        [1 0 6 1 2]
        [6 1 0 0 1]
        [0 6 3 6 0]
        sage: Z = lift_cross_ratios(A, to_sixth_root_of_unity)                          # needs sage.rings.finite_rings sage.rings.number_field
        sage: Z                                                                         # needs sage.rings.finite_rings sage.rings.number_field
        [     1      0      1      1      1]
        [-z + 1      1      0      0      1]
        [     0     -1      1 -z + 1      0]
        sage: M = LinearMatroid(reduced_matrix=A)
        sage: sorted(M.cross_ratios())
        [3, 5]
        sage: N = LinearMatroid(reduced_matrix=Z)                                       # needs sage.rings.finite_rings sage.rings.number_field
        sage: sorted(N.cross_ratios())                                                  # needs sage.rings.finite_rings sage.rings.number_field
        [-z + 1, z]
        sage: M.is_isomorphism(N, {e:e for e in M.groundset()})                         # needs sage.rings.finite_rings sage.rings.number_field
        True
    """
    from sage.graphs.graph import Graph

    for s, t in lift_map.items():
        source_ring = s.parent()
        target_ring = t.parent()
        break
    plus_one1 = source_ring(1)
    minus_one1 = source_ring(-1)
    plus_one2 = target_ring(1)
    minus_one2 = target_ring(-1)

    G = Graph([((r, 0), (c, 1), (r, c)) for r, c in A.nonzero_positions()])
    # write the entries of (a scaled version of) A as products of cross ratios of A
    T = Graph()
    for C in G.connected_components_subgraphs():
        T.add_edges(C.min_spanning_tree())
    # - fix a tree of the support graph G to units (= empty dict, product of 0 terms)
    F = {entry: dict() for entry in T.edge_labels()}
    W = set(G.edge_iterator()) - set(T.edge_iterator())
    H = G.subgraph(edges=T.edge_iterator())
    while W:
        # - find an edge in W to process, closing a circuit in H which is induced in G
        edge = W.pop()
        path = H.shortest_path(edge[0], edge[1])
        path_s = set(path)
        retry = True
        while retry:
            retry = False
            for edge2 in W:
                if edge2[0] in path_s and edge2[1] in path_s:
                    W.add(edge)
                    edge = edge2
                    W.remove(edge)
                    path = H.shortest_path(edge[0], edge[1])
                    path_s = set(path)
                    retry = True
                    break
        entry = edge[2]
        entries = []
        for i in range(len(path) - 1):
            v = path[i]
            w = path[i + 1]
            if v[1] == 0:
                entries.append((v[0], w[0]))
            else:
                entries.append((w[0], v[0]))
        # - compute the cross ratio `cr` of this whirl
        cr = source_ring(A[entry])
        div = True
        for entry2 in entries:
            if div:
                cr /= A[entry2]
            else:
                cr *= A[entry2]
            div = not div

        monomial = {}
        if len(path) % 4 == 0:
            if not cr == plus_one1:
                monomial[cr] = 1
        else:
            cr = -cr
            if cr != plus_one1:
                monomial[cr] = 1
            if minus_one1 in monomial:
                monomial[minus_one1] = monomial[minus_one1] + 1
            else:
                monomial[minus_one1] = 1

        if cr != plus_one1 and cr not in lift_map:
            raise ValueError("Input matrix has a cross ratio " + str(cr) + ", which is not in the lift_map")
        # - write the entry as a product of cross ratios of A
        div = True
        for entry2 in entries:
            if div:
                for cr, degree in F[entry2].items():
                    if cr in monomial:
                        monomial[cr] = monomial[cr] + degree
                    else:
                        monomial[cr] = degree
            else:
                for cr, degree in F[entry2].items():
                    if cr in monomial:
                        monomial[cr] = monomial[cr] - degree
                    else:
                        monomial[cr] = -degree
            div = not div
        F[entry] = monomial
        # - current edge is done, can be used in next iteration
        H.add_edge(edge)

    # compute each entry of Z as the product of lifted cross ratios
    Z = Matrix(target_ring, A.nrows(), A.ncols())
    for entry, monomial in F.items():
        Z[entry] = plus_one2
        for cr, degree in monomial.items():
            if cr == minus_one1:
                Z[entry] = Z[entry] * (minus_one2**degree)
            else:
                Z[entry] = Z[entry] * (lift_map[cr]**degree)

    return Z


def lift_map(target):
    r"""
    Create a lift map, to be used for lifting the cross ratios of a matroid
    representation.

    .. SEEALSO::

        :meth:`lift_cross_ratios() <sage.matroids.utilities.lift_cross_ratios>`

    INPUT:

    - ``target`` -- string describing the target (partial) field

    OUTPUT: dictionary

    Depending on the value of ``target``, the following lift maps will be created:

    - "reg": a lift map from `\GF3` to the regular partial field `(\ZZ, <-1>)`.

    - "sru": a lift map from `\GF7` to the
      sixth-root-of-unity partial field `(\QQ(z), <z>)`, where `z` is a sixth root
      of unity. The map sends 3 to `z`.

    - "dyadic": a lift map from `\GF{11}` to the dyadic partial field `(\QQ, <-1, 2>)`.

    - "gm": a lift map from `\GF{19}` to the golden mean partial field
      `(\QQ(t), <-1,t>)`, where `t` is a root of `t^2-t-1`. The map sends `5` to `t`.

    The example below shows that the latter map satisfies three necessary conditions stated in
    :meth:`lift_cross_ratios() <sage.matroids.utilities.lift_cross_ratios>`

    EXAMPLES::

        sage: from sage.matroids.utilities import lift_map
        sage: lm = lift_map('gm')                                                       # needs sage.rings.finite_rings sage.rings.number_field
        sage: for x in lm:                                                              # needs sage.rings.finite_rings sage.rings.number_field
        ....:     if (x == 1) is not (lm[x] == 1):
        ....:         print('not a proper lift map')
        ....:     for y in lm:
        ....:         if (x+y == 0) and not (lm[x]+lm[y] == 0):
        ....:             print('not a proper lift map')
        ....:         if (x+y == 1) and not (lm[x]+lm[y] == 1):
        ....:             print('not a proper lift map')
        ....:         for z in lm:
        ....:             if (x*y==z) and not (lm[x]*lm[y]==lm[z]):
        ....:                 print('not a proper lift map')
    """
    from sage.rings.finite_rings.finite_field_constructor import GF

    if target == "reg":
        R = GF(3)
        return {R(1): ZZ(1)}

    if target == "sru":
        from sage.rings.number_field.number_field import NumberField

        R = GF(7)
        z = ZZ['z'].gen()
        S = NumberField(z * z - z + 1, 'z')
        z = S(z)
        return {R.one(): S.one(), R(3): z, R(3)**(-1): z**5}

    if target == "dyadic":
        R = GF(11)
        return {R(1): QQ(1), R(-1): QQ(-1), R(2): QQ(2), R(6): QQ((1, 2))}

    if target == "gm":
        from sage.rings.number_field.number_field import NumberField

        R = GF(19)
        t = QQ['t'].gen()
        G = NumberField(t * t - t - 1, 't')
        return {R(1): G(1), R(5): G(t),
                R(1) / R(5): G(1) / G(t), R(-5): G(-t),
                R(-5)**(-1): G(-t)**(-1), R(5)**2: G(t)**2,
                R(5)**(-2): G(t)**(-2)}

    raise NotImplementedError(target)


def split_vertex(G, u, v=None, edges=None):
    """
    Split a vertex in a graph.

    If an edge is inserted between the vertices after splitting, this
    corresponds to a graphic coextension of a matroid.

    INPUT:

    - ``G`` -- a SageMath :class:`Graph`
    - ``u`` -- a vertex in ``G``
    - ``v`` -- (optional) the name of the new vertex after the splitting. If
      ``v`` is specified and already in the graph, it must be an isolated vertex.
    - ``edges`` -- (optional) iterable container of edges on ``u`` that
      move to ``v`` after the splitting. If ``None``, ``v`` will be an isolated
      vertex. The edge labels must be specified.

    EXAMPLES::

        sage: # needs sage.graphs
        sage: from sage.matroids.utilities import split_vertex
        sage: G = graphs.BullGraph()
        sage: split_vertex(G, u=1, v=55, edges=[(1, 3)])
        Traceback (most recent call last):
        ...
        ValueError: the edges are not all incident with u
        sage: split_vertex(G, u=1, v=55, edges=[(1, 3, None)])
        sage: list(G.edges(sort=True))
        [(0, 1, None), (0, 2, None), (1, 2, None), (2, 4, None), (3, 55, None)]
    """
    if v is None:
        v = G.add_vertex()
    elif v not in G:
        G.add_vertex(v)
    elif G.degree(v):
        raise ValueError("v must be a new vertex or an isolated vertex")
    if edges is None:
        edges = []

    edges_on_u = G.edges_incident(u)

    for e in edges:
        if e not in edges_on_u:
            # if e is a loop, put it on u and v
            # otherwise raise an error
            if e[0] == e[1]:
                G.add_edge(u, v, e[2])
                G.delete_edge(e)
            else:
                raise ValueError("the edges are not all incident with u")

        elif e[0] == u:
            G.add_edge(v, e[1], e[2])
        elif e[1] == u:
            G.add_edge(e[0], v, e[2])
        G.delete_edge(e)

    # This modifies the graph without needing to return anything


def cmp_elements_key(x):
    """
    A helper function to compare elements which may be integers or strings.

    EXAMPLES::

        sage: from sage.matroids.utilities import cmp_elements_key
        sage: l = ['a', 'b', 1, 3, 2, 10, 111, 100, 'c', 'aa']
        sage: sorted(l, key=cmp_elements_key)
        [1, 2, 3, 10, 100, 111, 'a', 'aa', 'b', 'c']
    """
    return (isinstance(x, str), x)
