"""
Independence Poset

Software implementation of independence posets. Includes functions to compute maximal orthogonal pairs, tight orthogonal pairs, rowmotion, labelling, and generating Galois graphs of lattices.

REFERENCES:

- [TW2017]_

- [TW2018]_

AUTHORS:

- Nathan Williams (2022-03-24): 1.0
- Deep Desai (2022-03-24): 1.0
- Ariba Tahsin (2022-03-24): 1.0
- Ryan Lofton (2022-03-24): 1.0
- Ramesh Kanakala (2022-03-24): 1.0
- Savitha Venkatesh (2022-03-24): 1.0

"""

# ****************************************************************************
#       Copyright (C) 2013 Nathan Williams <Nathan.Williams1@utdallas.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.graphs.digraph import DiGraph
from sage.combinat.posets.posets import Poset
from sage.sets.set import Set
from sage.combinat.subset import Subsets

##############################################################
# code for maximal orthogonal pairs
##############################################################

def _inc(a, b):
    '''
    INPUT:

    ``a`` -- a pair of sets

    ``b`` -- a pair of sets

    OUTPUT:

    compares two orthogonal pairs by the inclusion relation on the first sets and reverse inclusion on the second sets.
    '''
    if a[0].issubset(b[0]) and a[1].issubset(b[1]):
        return True
    else:
        return False


def _inc1(a, b):
    '''
    INPUT:

    ``a`` -- a pair of sets

    ``b`` -- a pair of sets

    OUTPUT:
    
    compares an orthogonal pair by inclusion of first sets.
    '''
    if a[0].issubset(b[0]):
        return True
    else:
        return False


def orthogonal_pairs(G):
    '''
    INPUT:

    ``G`` -- an acyclic directed graph

    OUTPUT:

    the set of all orthogonal pairs of `G`, defined as a pair of sets `(X, Y)` such that there is no edge from a node in `X` to a node in `Y`

    .. NOTE:

    This is a slow, brute force implementation.

    EXAMPLES:

    The orthogonal pairs for a small graph ::

        sage: from sage.combinat.independence_posets import orthogonal_pairs
        sage: G = Graph([[1, 2], [(1, 2)]])
        sage: G.show()
        sage: orthogonal_pairs(G)
        [({}, {}), ({}, {1}), ({}, {2}), ({}, {1, 2}), ({1}, {}), ({2}, {}), ({2}, {1}), ({1, 2}, {})]

    The orthogonal pairs for a complete graph ::

        sage: G = graphs.CompleteGraph(3)
        sage: G.show()
        sage: orthogonal_pairs(G)
        [({}, {}), ({}, {0}), ({}, {1}), ({}, {2}), ({}, {0, 1}), ({}, {0, 2}), ({}, {1, 2}), ({}, {0, 1, 2}), ({0}, {}), ({1}, {}), ({1}, {0}), ({2}, {}), ({2}, {0}), ({2}, {1}), ({2}, {0, 1}), ({0, 1}, {}), ({0, 2}, {}), ({1, 2}, {}), ({1, 2}, {0}), ({0, 1, 2}, {})]

    '''
    sets = []
    edges = [(i[0], i[1]) for i in G.edges()]
    for s in Subsets(G):
        for t in Subsets(Set(G).difference(s)):
            no = False
            for i in s:
                for j in t:
                    if (i, j) in edges:
                        no = True
                        break
                if no is True:
                    break
            if no is False:
                sets.append((s, t))
    return sets


def poset_of_orthogonal_pairs(G):
    '''
    INPUT:

    ``G`` -- an acyclic directed graph

    OUTPUT:

    a poset of orthogonal pairs of g with the inclusion relations on both elements

    EXAMPLES:

    The poset of orthogonal pairs for a small graph ::

        sage: from sage.combinat.independence_posets import poset_of_orthogonal_pairs
        sage: G = Graph([[1, 2], [(1, 2)]])
        sage: poset_of_orthogonal_pairs(G)
        Finite poset containing 8 elements

    The poset of orthogonal pairs for a completely connected graph ::

        sage: G = graphs.CompleteGraph(3)
        sage: poset_of_orthogonal_pairs(G)
        Finite poset containing 20 elements
    '''
    return Poset([orthogonal_pairs(G), _inc])


def mops(G):
    '''
    INPUT:

    ``G`` -- an acyclic directed graph

    OUTPUT:

    pairs of sets `(X, Y)` containing the maximal elements in a poset of orthogonal pairs with respect to the condition that no edges run from `X` to `Y`.

    EXAMPLES:

    The maximal orthogonal pairs for a connected two-vertex graph ::

        sage: from sage.combinat.independence_posets import mops
        sage: G = Graph([[1, 2], [(1, 2)]])
        sage: mops(G)
        [({}, {1, 2}), ({2}, {1}), ({1, 2}, {})]

    The maximal orthogonal pairs for a completely connected 3 vertex graph ::

        sage: G = graphs.CompleteGraph(3)
        sage: mops(G)
        [({2}, {0, 1}), ({}, {0, 1, 2}), ({1, 2}, {0}), ({0, 1, 2}, {})]
    '''
    return poset_of_orthogonal_pairs(G).maximal_elements()


def poset_of_mops(G):
    '''
    INPUT:

    ``G`` -- an acyclic directed graph

    OUTPUT:

    a poset of the maximal orthogonal pairs of ``G`` with the inclusion relation for only the first element

    EXAMPLES:

    The poset of maximal orthogonal pairs for a connected two-vertex graph ::

        sage: from sage.combinat.independence_posets import poset_of_mops
        sage: G = Graph([[1, 2], [(1, 2)]])
        sage: poset_of_mops(G)
        Finite poset containing 3 elements

    The poset of maximal orthogonal pairs for a completely connected 3 vertex graph ::

        sage: G = graphs.CompleteGraph(3)
        sage: poset_of_mops(G)
        Finite poset containing 4 elements
    '''
    return Poset([mops(G), _inc1])


def complete_mop(G, S):
    '''
    INPUT:

    - ``G`` -- an acyclic directed graph

    - ``S`` -- a subset of the graph ``G`` of type ``Set``

    OUTPUT:

    a set containing the maximal elements in a poset of orthogonal pairs `(X, Y)` such that `X \\subset S`

    EXAMPLES:
    
    Example of the ``complete_mop`` ::

        sage: from sage.combinat.independence_posets import complete_mop
        sage: G = DiGraph([[1..5], lambda i, j: i < j])
        sage: S = Set([2..4])
        sage: complete_mop(G, S)
        [({2, 3, 4}, {1})]

    .. NOTE:

    This is a slow, brute force implementation.
    '''
    sets = []
    edges = [(i[0], i[1]) for i in G.edges()]
    for t in Subsets(Set(G).difference(S)):
        no = False
        for i in S:
            for j in t:
                if (i, j) in edges:
                    no = True
                    break
            if no is True:
                break
        if no is False:
            sets.append((S, t))
    P = Poset([sets, _inc])
    return P.maximal_elements()


##############################################################
# code for left modular lattices
##############################################################
def is_left_modular(L, H=None, verbose=False):
    '''
    INPUT:

    ``L`` -- left-modular lattice

    ``H`` -- subset of elements; ``H`` is taken as ``L`` if none is given

    ``verbose`` -- indicates whether to give a list of failures; false by default

    OUTPUT:

    if ``verbose == True``, outputs a list of tuples `(y, x, z)` which fail left-modularity.  if ``verbose == False``, outputs ``False`` if any one `x \\in H` fails to be left-modular and ``True`` otherwise.

    ALGORITHM:

    Given a lattice `L` and a subset of elements `H`, an element `x \\in H` is left-modular if for every `y,z \\in L, y \\leq z` the equality `(y \\vee x) \\wedge z = y \\vee (x \\wedge z)`. 
    
    EXAMPLES:

    Testing a lattice that isn't left-modular ::

        sage: from sage.combinat.independence_posets import is_left_modular
        sage: L = LatticePoset(([1,2,3,4,5],[[1,2],[1,3],[3,4],[4,5],[2,5]]))
        sage: is_left_modular(L)
        False

    Testing a left-modular lattice ::

        sage: L = LatticePoset(([1,2,3,4,5,6],[[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: is_left_modular(L)
        True
    '''
    if H == None:
        H = L
    out = []
    for x in H:
        for z in L:
            for y in L.principal_lower_set(z):
                if (L.join(y, L.meet(x, z)) != L.meet(L.join(y, x), z)):
                    if verbose is False:
                        return False
                    else:
                        out += [(y, x, z)]
    if verbose is False:
        return True
    else:
        return out


def join_irr_label(L, c, j):
    '''
    INPUT:

    ``L`` -- a left-modular lattice

    ``c`` -- a left-modular chain in ``L``

    ``j`` -- join-irreducible element

    OUTPUT:

    returns the label of the join-irreducible element ``j``, defined as the minimum index `i` such that `j < x_i`, where `x_i \\in c`
    
    EXAMPLES ::
    
        sage: from sage.combinat.independence_posets import join_irr_label
        sage: L = LatticePoset(([1,2,3,4,5,6],[[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,3,4,5]
        sage: join_irr_label(L, chain, 5)
        3
    '''
    for i in range(len(c)):
        if L.le(j, c[i]):
            return i


def meet_irr_label(L, c, m):
    '''
    INPUT:

    ``L`` -- a left-modular lattice

    ``c`` -- a left-modular chain in ``L``

    ``m`` -- meet-irreducible element

    OUTPUT:

    returns the label of the meet-irreducible element ``m``, defined as the maximum index `i` such that `m > x_{i-1}`, where `x_{i-1} \\in c`
    
    EXAMPLES ::

        sage: from sage.combinat.independence_posets import meet_irr_label
        sage: L = LatticePoset(([1,2,3,4,5,6],[[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,3,4,5]
        sage: meet_irr_label(L, chain, 5)
        4
    '''
    for i in range(len(c) - 1, -1, -1):
        if L.ge(m, c[i]):
            return i + 1


def downward_labels(L, c, z):
    '''
    INPUT:

    ``L`` -- a trim lattice L

    ``c`` -- a left-modular chain c

    ``z`` -- an element z in the lattice ``L``

    OUTPUT:

    the list of downward labels for z, defined as `\\min\\{\\beta_J(j) | j \\in J, y \\vee j = z\\}` where `J` is the set of all join-irreducibles in `L` and `\\beta_J(j)` is the join-irreducible label of `j` for all `y \\in L, y \\leq z`
    
    EXAMPLES ::

        sage: from sage.combinat.independence_posets import downward_labels
        sage: L = LatticePoset(([1,2,3,4,5,6],[[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: downward_labels(L, chain, 4)
        [2, 1]
        sage: L = LatticePoset(([1,2,3,4,5,6],[[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: downward_labels(L, chain, 5)
        [2, 3]
    '''
    J = L.join_irreducibles()
    JL = dict([(j, join_irr_label(L, c, j)) for j in J])
    cov = L.lower_covers(z)
    return [min([JL[j] for j in J if L.join([y, j]) == z]) for y in cov]


def upward_labels(L, c, z):
    '''
    INPUT:

    ``L`` -- a trim lattice L

    ``c`` -- a left-modular chain c

    ``z`` -- an element z in the lattice ``L``

    OUTPUT:

    the list of upward labels for ``z``, defined as `\\max\\{\\beta_M(m) | m \\in M, z \\wedge m = y\\}` where `M` is the set of all meet-irreducibles in `L` and `\\beta_M(m)` is the meet-irreducible label of `m` for all `y \\in L, y \\geq z`
    
    EXAMPLES ::
    
        sage: from sage.combinat.independence_posets import upward_labels
        sage: L = LatticePoset(([1,2,3,4,5,6],[[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: upward_labels(L, chain, 1)
        [1, 2]
        sage: L = LatticePoset(([1,2,3,4,5,6],[[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: upward_labels(L, chain, 2)
        [3, 2]
    '''
    M = L.meet_irreducibles()
    ML = dict([(m, meet_irr_label(L, c, m)) for m in M])
    cov = L.upper_covers(z)
    return [max([ML[m] for m in M if L.meet([y, m]) == z]) for y in cov]


def cover_label(L, c, edge):
    '''
    INPUT:

    ``L`` -- a trim lattice L

    ``c`` -- a left-modular chain c

    ``edge`` -- a covering relation ``[x,y]`` where `x \\leq y`

    OUTPUT:

    the label for the covering relation ``edge`` in ``L``
    
    EXAMPLES ::

        sage: from sage.combinat.independence_posets import cover_label
        sage: L = LatticePoset(([1,2,3,4,5,6],[[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: cover_label(L, chain, [2, 4])
        2
        sage: L = LatticePoset(([1,2,3,4,5,6],[[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: chain = [1,2,4,5]
        sage: cover_label(L, chain, [2, 6])
        3
    '''
    x, y = edge
    J = L.join_irreducibles()
    JL = dict([(j, join_irr_label(L, c, j)) for j in J])
    return min([JL[j] for j in J if L.join([x, j]) == y])


##############################################################
# trim lattices
##############################################################
def is_trim(L):
    '''
    INPUT:

    - ``L`` -- a lattice

    OUTPUT:
    returns a tuple ``(True, maximal_chains)`` if ``L`` is a trim lattice, where ``maximal_chains`` is a list of maximal chains in ``L``; returns ``(False, False)`` if ``L`` is not a trim lattice. ``L`` is trim, if there exists a maximal chain of length `n+1` containing only left-modular elements, and there are exactly `n` join-irreducible elements and meet-irreducible elements.

    EXAMPLES:

    Testing a trim lattice ::

        sage: from sage.combinat.independence_posets import is_trim
        sage: L = LatticePoset(([1,2,3,4,5,6],[[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5],[2,4]]))
        sage: is_trim(L)
        (True, [[1, 2, 6, 5], [1, 2, 4, 5], [1, 3, 4, 5]])

    Testing a lattice which is not trim ::

        sage: L = LatticePoset(([1,2,3,4,5,6],[[1,2],[1,3],[3,4],[4,5],[2,5],[2,6],[6,5]]))
        sage: is_trim(L)
        (False, False)

    NOTE:
    computes all maximal length chains of ``L`` which can take some time.
    '''
    J = L.join_irreducibles()
    M = L.meet_irreducibles()
    C = L.maximal_chains()
    n = max(map(len, C))
    #print n,len(J),len(M)
    maximal_chains = [c for c in C if len(c) == n]
    #print is_left_modular(mC[0])
    if len(J) == len(M) and len(J) == n - 1 \
            and is_left_modular(L, maximal_chains[0]):
        return (True, maximal_chains)
    else:
        return (False, False)


##############################################################
# tight orthogonal pairs (tops)
##############################################################
def complete_top(G, u):
    '''
    INPUT:

    - ``G`` -- a finite acyclic directed graph

    - ``u`` -- an independent set in ``G``

    OUTPUT:
    returns a tuple (``m``, ``u``) independent sets in ``G`` which are orthogonal (there is no edge from ``m`` to ``u`` in ``G``) and tight (if an element in ``m`` is increased (removed and replaced by a larger element) or any element of ``u`` is decreased, or a new element is added into either ``m`` or ``u`` then the result will no longer be an orthogonal pair)

    EXAMPLES ::
    
        sage: from sage.combinat.independence_posets import complete_top
        sage: G = DiGraph([[1, 2, 3, 4], [[4,3], [3,2], [2,1]]])
        sage: complete_top(G, Set([2, 4]))
        ({1}, {2, 4})
        sage: complete_top(G, Set([1]))
        ({4}, {1})
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9], [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: complete_top(G, Set([1,4]))
        ({8, 3, 7}, {1, 4})
        sage: complete_top(G, Set([5, 6, 9]))
        ({1, 4}, {9, 5, 6})
    '''
    P = Poset(G)
    m = Set([])
    for k in P.linear_extension():
        ngk = Set([i[0] for i in G.incoming_edges(k)])
        nlk = Set([i[1] for i in G.outgoing_edges(k)])
        if ngk.intersection(m) == Set([]) \
                and nlk.intersection(u) == Set([]) and k not in u:
            m = m.union(Set([k]))
    return (m, u)


def top_to_mop(G, du):
    '''
    INPUT:

    - ``G`` -- a finite acyclic directed graph

    - ``du`` -- a tight orthogonal pair (``d``, ``u``)

    OUTPUT:
    a maximal orthogonal pair corresponding to the top (``d``, ``u``) given

    EXAMPLES ::

        sage: from sage.combinat.independence_posets import complete_top, top_to_mop
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9], [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([4, 7]))
        sage: top_to_mop(G, du)
        ({1, 3, 6, 8, 9}, {4, 7})
    '''
    X = du[0]
    Y = du[1]
    for b in Set(G).difference(Y):
        if Set(G.neighbors_out(b)).intersection(Y) == Set([]):
            X = X.union(Set([b]))
    for a in Set(G).difference(X):
        if Set(G.neighbors_in(a)).intersection(X) == Set([]):
            Y = Y.union(Set([a]))
    return (X, Y)


def maximal_top(G):
    '''
    INPUT:

    - ``G`` -- a finite acyclic directed graph

    OUTPUT:
    a tuple (``m``, ``u``) where ``m`` is the set which is tight and orthogonal to an empty ``u``.

    EXAMPLES ::

        sage: from sage.combinat.independence_posets import maximal_top
        sage: G = DiGraph([[1, 2, 3, 4], [[4,3], [3,2], [2,1]]])
        sage: maximal_top(G)
        ({2, 4}, {})
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9], [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: maximal_top(G)
        ({1, 4, 5, 6, 9}, {})
    '''
    return complete_top(G, Set([]))


def minimal_top(G):
    '''
    INPUT:

    - ``G`` -- a finite acyclic directed graph

    OUTPUT:
    a tuple (``m``, ``u``) where ``u`` is the set which is tight and orthogonal to an empty ``m``.

    EXAMPLES ::
    
        sage: from sage.combinat.independence_posets import minimal_top
        sage: G = DiGraph([[1, 2, 3, 4], [[4,3], [3,2], [2,1]]])
        sage: minimal_top(G)
        ({}, {1, 3})
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9], [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: minimal_top(G)
        ({}, {1, 4, 5, 6, 9})
    '''
    P = Poset(G)
    m = Set([])
    for k in reversed(P.linear_extension()):
        if Set([i[1] for i in G.outgoing_edges(k)]).intersection(m) == Set([]):
            m = m.union(Set([k]))
    return (Set([]), m)


def flip(G, du, j):
    '''
    INPUT:

    - ``G`` -- a finite acyclic directed graph

    - ``du`` -- a top (D, U)

    - ``j`` -- an element `j`

    OUTPUT:
    
    Performs a flip if `j \\in U` or `j \\in D`. Returns a new top ``(d, u)`` after the flip.

    EXAMPLES:

    Flipping on an element in U ::

        sage: from sage.combinat.independence_posets import complete_top, flip
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9], [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip(G, du, 4)
        ({8, 3, 4}, {2})

    Flipping on an element in D ::

        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9], [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip(G, du, 3)
        ({9, 5, 6}, {3, 4})

    Flipping on an element in neither D nor U ::

        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9], [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip(G, du, 5)
        ({8, 3, 7}, {1, 4})
    '''
    P = Poset(G)
    a = du[0]
    b = du[1]
    if j not in b and j not in a:
        return du
    elif j in b:
        a2 = Set([k for k in a if not (P.ge(k, j))] + [j])
        b2 = Set([k for k in b if not (P.le(k, j))])
    elif j in a:
        a2 = Set([k for k in a if not (P.ge(k, j))])
        b2 = Set([k for k in b if not (P.le(k, j))] + [j])
    for k in P.linear_extension():
        if not (P.le(k, j)) \
                and Set([i[1] for i in G.outgoing_edges(k)]).intersection(b) == Set([]) \
                and Set([i[0] for i in G.incoming_edges(k)]).intersection(a2) == Set([]) \
                and (k not in b2):
            a2 = a2.union(Set([k]))
    for k in reversed(P.linear_extension()):
        if not (P.ge(k, j)) \
                and Set([i[0] for i in G.incoming_edges(k)]).intersection(a) == Set([]) \
                and Set([i[1] for i in G.outgoing_edges(k)]).intersection(b2) == Set([]) \
                and (k not in a2):
            b2 = b2.union(Set([k]))
    return (a2, b2)


def flip_up(G, du, j):
    '''
    INPUT:

    - ``G`` -- a finite acyclic directed graph

    - ``du`` -- a top (D, U)

    - ``j`` -- an element j in U

    OUTPUT:
    
    conducts a flip only if `j \\in U`. Returns a new top ``(d, u)`` after the flip

    EXAMPLES:

    Flipping on an element in U ::

        sage: from sage.combinat.independence_posets import complete_top, flip_up
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9], [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip_up(G, du, 4)
        ({8, 3, 4}, {2})

    Flipping on an element in D ::

        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9], [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip_up(G, du, 3)
        ({8, 3, 7}, {1, 4})

    Flipping on an element in neither D nor U ::

        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9], [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1,4]))
        sage: flip_up(G, du, 5)
        ({8, 3, 7}, {1, 4})
    '''
    if j not in du[1]:
        return du
    elif j in du[1]:
        return flip(G, du, j)


def flips(G, du, s):
    '''
    INPUT:

    - ``G`` -- a directed acyclic graph

    - ``du`` -- a top (D, U)

    - ``s`` -- a sequence of elements `j` in G

    OUTPUT:
    
    returns a new top ``(d, u)`` after flipping on each element in ``s`` (i.e `\\displaystyle{ \\left( \\prod_{j \\in s} \\text{flip}_j \\right)(D, U) }`)

    EXAMPLES ::
    
        sage: from sage.combinat.independence_posets import complete_top, flip, flips
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6, 7, 8, 9], [[1,2], [1,3], [2,4], [2,5], [3,5], [3,6], [4,7], [5,7], [5,8], [6,8], [7,9], [8,9]]])
        sage: du = complete_top(G, Set([1, 4]))
        sage: a = flip(G, du, 4)
        sage: b = flip(G, a, 2)
        sage: flip(G, b, 3)
        ({2, 6, 7}, {3})
        sage: flips(G, du, [4, 2, 3])
        ({2, 6, 7}, {3})

    '''
    if not type(s) is list:
        s_iter = [s]
    else:
        s_iter = s
    dv = tuple(du)
    for i in s_iter:
        dv = flip(G, dv, i)
    return dv


def make_flip_poset(G):
    '''
    INPUT:

    - ``G`` -- a finite directed acyclic graph

    OUTPUT:
    
    returns the independence poset corresponding to the graph ``G``.

    EXAMPLES ::
    
        sage: from sage.combinat.independence_posets import make_flip_poset
        sage: G = DiGraph([[1, 2, 3, 4], [[4, 3], [3, 2], [2, 1]]])
        sage: make_flip_poset(G)
        Finite poset containing 8 elements

    The below graph generates a Tamari lattice with 14 elements ::

        sage: G = DiGraph([[1, 2, 3, 4, 5, 6], [[6, 5], [6, 4], [5, 4], [5, 3], [5, 2], [4, 2], [4, 1], [3, 2], [2, 1]]])
        sage: make_flip_poset(G) # this generates the Tamari lattice with 14 elements
        Finite poset containing 14 elements
    '''
    edges = []
    els = Set([])
    last = []
    new = [minimal_top(G)]
    while new != []:
        last = new
        els = els.union(Set(last))
        new = []
        for a in last:
            for i in a[1]:
                b = flip(G, a, i)
                edges.append([a, b])
                if b not in new:
                    new = new + [b]
    return Poset((list(els), edges))

def rowmotion(L, x):
    '''
    INPUT:

    ``L`` -- trim lattice

    ``x`` -- element x to apply rowmotion on

    OUTPUT:
    
    Rowmotion (with respect with a particular labeling on the lattice `L`) on an element x is the unique element `y` such that the set of labels of upper covers `s > y` equals the set of labels of lower covers `s < x`
    
    EXAMPLES:

    Example of rowmotion on an independent poset ::

        sage: from sage.combinat.independence_posets import make_flip_poset, rowmotion
        sage: G = DiGraph([[1, 2, 3, 4], [[4, 3], [3, 2], [2, 1]]])
        sage: L = make_flip_poset(G)
        sage: x = L.list()[3]
        sage: rowmotion(L, x)
        ({1, 3}, {4})
    '''
    lower = dict([(l, l[0]) for l in L])
    upper = dict([(l[1], l) for l in L])
    return upper[lower[x]]

def toggle(G, x):
    '''
    INPUT:

    ``G`` -- a finite directed acyclic graph

    ``x`` -- an element in ``G``

    OUTPUT:
    the directed acyclic graph given by reversing all edges incident to `x` in G
    
    EXAMPLES:
    
    Example of a usage of toggle. We take a graph, ``G`` and compute the independent poset on ``G`` and ``tog_3(G)`` ::

        sage: from sage.combinat.independence_posets import toggle, make_flip_poset
        sage: G = DiGraph([[1, 2, 3], [[3, 2], [2, 1]]])
        sage: G_3 = toggle(G, 3)
        sage: P_3 = make_flip_poset(G_3)
        sage: P_3
        Finite poset containing 5 elements
    '''
    H = DiGraph(G)
    for i in H.edges():
        if x in i:
            H.reverse_edge(i)
    return (H)

def toggles(G, s):
    '''
    INPUT:

    - ``G`` -- a finite directed acyclic graph

    - ``s`` -- a sequence of elements

    OUTPUT:
    
    toggles all elements `x` in graph `G` (i.e. reverses all edges incident to `x`) as they appear in the sequence `s` starting at index 0.

    EXAMPLES:

    Example of a usage of toggles. We take a graph, `G` and compute the independent poset on `G` and `tog_2(tog_4(G))` two different ways: once manually using ``toggle`` and once using ``toggles`` ::

        sage: from sage.combinat.independence_posets import toggle, toggles, make_flip_poset
        sage: G = DiGraph([[1, 2, 3, 4, 5, 6], [[6, 5], [6, 4], [5, 4], [5, 3], [5, 2], [4, 2], [4, 1], [3, 2], [2, 1]]])
        sage: P = make_flip_poset(G)
        sage: G_2 = toggle(G, 2)
        sage: G_24_1 = toggle(G_2, 4)
        sage: G_24_2 = toggles(G, [2, 4])
        sage: G_24_1 == G_24_2
        True

    ''toggles'' will also accept a singular element for ``s``; in which case it will be the same thing as a ``toggle`` ::

        sage: G = DiGraph([[1, 2, 3], [[3, 2], [2, 1]]])
        sage: P_1 = toggles(G, 2)
        sage: P_2 = toggle(G, 2)
        sage: P_1 == P_2
        True
    '''
    if not type(s) is list:
        s_iter = [s]
    else:
        s_iter = s
    H = DiGraph(G)
    for i in s_iter:
        H = DiGraph(toggle(H, i))
        H.show()
    return(H)