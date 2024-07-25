# sage.doctest: needs sage.graphs
"""
Helper functions related to Gauss codes of knots

This contains in particular a function that takes a signed Gauss code
and builds the list of local orientations at every crossing. This involves
making choices, but allows to build a knot from just a signed Gauss code.

The following convention is used for signed Gauss codes: a negative
number corresponds to an over-crossing.

REFERENCES:

- https://math.stackexchange.com/questions/1233249/recovering-knot-crossing-orientations-from-a-gauss-code-or-dowker-notation

.. [Kauf1999] Louis H. Kauffman, *Virtual Knot Theory*, Europ. J. Combinatorics
   (1999) 20, 663-691 ; http://homepages.math.uic.edu/~kauffman/VKT.pdf
"""
# ****************************************************************************
#       Copyright (C) 2019 Frédéric Chapoton <fchapoton2@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.graphs.graph import Graph


def dowker_to_gauss(code):
    """
    Convert from Dowker-Thistlethwaite code to signed Gauss code.

    EXAMPLES::

        sage: from sage.knots.gauss_code import dowker_to_gauss
        sage: dowker_to_gauss([6,-12,2,8,-4,-10])
        [-3, 1, 6, -2, -1, 3, -4, 4, 2, -5, 5, -6]
        sage: dowker_to_gauss([-4,-6,-2])
        [2, -1, 3, -2, 1, -3]

    TESTS::

        sage: dowker_to_gauss([])
        []
    """
    n = len(code)
    # vertices are numbered by half of their even label
    signes = {abs(j): (1 if j > 0 else -1) for j in code}
    gauss = []
    for i in range(1, 2 * n + 1):
        if i % 2:
            letter = code[(i - 1) // 2] // 2
            gauss.append(-letter)
        else:
            gauss.append(signes[abs(i)] * i // 2)
    return gauss


def recover_orientations(gauss):
    """
    Create diagrammatic information from signed Gauss code.

    This method is an auxiliary method, used for two different
    goals. The first goal is to create a knot from the signed Gauss
    code. This requires choosing at every crossing a local
    orientation, in a compatible way. The second goal is to picture
    the knot diagram using unicode. This goes through the creation of
    a noncrossing matching diagram.

    INPUT:

    - signed Gauss code

    OUTPUT:

    - word, a permutation of the (unsigned) Gauss code
    - list of positive (upwards) couplings in this word
    - list of negative (downwards) couplings in this word
    - list of signs, local orientations at each crossing

    ALGORITHM:

    The algorithm comes from Section 3 of [Kauf1999]_.

    EXAMPLES::

        sage: from sage.knots.gauss_code import recover_orientations
        sage: G = [4,1,5,2,1,3,6,4,2,5,3,6]
        sage: recover_orientations(G)
        ([4, 2, 1, 4, 6, 3, 1, 5, 2, 5, 3, 6],
         [(1, 8), (2, 6)],
         [(0, 3), (4, 11), (5, 10), (7, 9)],
         [1, -1, 1, 1, -1, 1])

        sage: G = [1,2,3,1,2,3]
        sage: recover_orientations(G)
        ([1, 3, 2, 1, 2, 3], [(0, 3)], [(1, 5), (2, 4)], [1, 1, 1])

        sage: G = [1, 2, 4, 5, 8, 1, 3, 4, 6, 7, 2, 3, 5, 6, 7, 8]
        sage: recover_orientations(G)
        ([1, 8, 7, 6, 5, 4, 6, 7, 2, 4, 3, 2, 1, 3, 5, 8],
         [(0, 12), (2, 7), (3, 6), (8, 11)],
         [(1, 15), (4, 14), (5, 9), (10, 13)],
         [1, -1, 1, -1, 1, -1, -1, 1])

    TESTS::

        sage: recover_orientations([])
        ([], [], [], [])
    """
    signs_overunder = [1 if letter > 0 else -1 for letter in gauss]
    gauss = [abs(x) for x in gauss]
    n = len(gauss) // 2

    # first step, perform subword reversals
    changed = list(gauss)
    for i in range(1, n + 1):
        id0 = changed.index(i)
        start = changed[:id0 + 1]
        after0 = changed[id0 + 1:]
        id1 = after0.index(i)
        middle = list(reversed(after0[:id1]))
        end = after0[id1:]
        changed = start + middle + end

    # second step, build the noncrossing planar diagram
    pairings = [[] for _ in range(n)]
    for pos, letter in enumerate(changed):
        pairings[letter - 1].append(pos)
    pairings = [tuple(ab) for ab in pairings]

    def cross(ab, cd):
        a, b = ab
        c, d = cd
        if c < a:
            a, b, c, d = c, d, a, b
        return a < c < b < d

    colors = Graph([pairings, cross]).coloring()
    if colors is None:  # empty graph
        positive = []
        negative = []
    elif len(colors) == 1:
        positive = colors[0]
        negative = []
    else:
        positive, negative = sorted(colors)
    positive.sort()
    negative.sort()

    # third step, determine the local orientations
    # this is done by following the knot

    signs_local = {d: 1 for d in range(1, n + 1)}
    for x, _ in negative:
        signs_local[changed[x]] *= -1

    jump_from = {}
    for x, y in positive + negative:
        jump_from[x] = y
        jump_from[y] = x

    direction = 1
    position = 0
    for k in range(2 * n):
        letter = changed[position]
        new_position = jump_from[position]
        signe = signs_local[letter]
        if new_position < position:
            signe = -signe
        if direction == 1:
            signe *= signs_overunder[k]
        signs_local[letter] = signe
        direction = -direction
        position = new_position + direction

    signs_final = [signs_local[d] for d in range(1, n + 1)]
    return changed, positive, negative, signs_final


def rectangular_diagram(gauss):
    """
    Return a rectangular diagram and crossing coordinates.

    INPUT:

    - signed Gauss code

    OUTPUT:

    - graph whose vertices are the corners of the knot diagram

    - positions of the horizontal and vertical crossings

    EXAMPLES::

        sage: from sage.knots.gauss_code import rectangular_diagram
        sage: G = [4,1,5,2,1,3,6,4,2,5,3,6]
        sage: rectangular_diagram(G)
        (Graph on 18 vertices, ([(1, 3), (3, 7), (4, 6), (6, 2)],
         [(4, 3), (6, 5)]))

        sage: G = [1,2,3,1,2,3]
        sage: rectangular_diagram(G)
        (Graph on 10 vertices, ([(1, 3), (2, 2)], [(2, 1)]))

    TESTS::

        sage: rectangular_diagram([])
        (Graph on 4 vertices, ([], []))
    """
    if not gauss:
        G = Graph()
        verts = [(0, 0), (1, 0), (1, 1), (0, 1)]
        for i in range(4):
            G.add_edge(verts[i], verts[(i + 1) % 4])
        G.set_pos({v: v for v in verts})
        return G, ([], [])

    changed, positive, negative, ori = recover_orientations(gauss)

    n = len(changed)
    G = Graph()
    horizontal = []
    vertical = []
    for a, b in positive:
        G.add_edge((a, a), (b, a))
        G.add_edge((b, a), (b, b))
        G.add_edge((a + 1, a + 1), (b + 1, a + 1))
        G.add_edge((b + 1, a + 1), (b + 1, b + 1))
        letter = changed[a]
        if b != a + 1:
            if ori[letter - 1] == -1:
                horizontal.append((b, a + 1))
            else:
                vertical.append((b, a + 1))
    for a, b in negative:
        G.add_edge((a, a), (a, b))
        G.add_edge((a, b), (b, b))
        G.add_edge((a + 1, a + 1), (a + 1, b + 1))
        G.add_edge((a + 1, b + 1), (b + 1, b + 1))
        letter = changed[a]
        if a + 1 != b:
            if ori[letter - 1] == 1:
                horizontal.append((a + 1, b))
            else:
                vertical.append((a + 1, b))

    # closing the loop on either side
    total = positive + negative
    debut = next(iter(b for a, b in total if a == 0))
    fin = next(iter(n - 1 - a for a, b in total if b == n - 1))
    if debut > fin:
        G.add_edge((0, 0), (0, n))
        G.add_edge((0, n), (n, n))
    else:
        G.add_edge((0, 0), (n, 0))
        G.add_edge((n, 0), (n, n))

    # remove useless vertices and double edges
    for a, b in list(G):
        (x, y), (xx, yy) = sorted(G.neighbors((a, b)))
        if x == xx == a:
            G.delete_vertex((a, b))
            G.add_edge((a, y), (a, yy))
        elif y == yy == b:
            G.delete_vertex((a, b))
            G.add_edge((x, b), (xx, b))

    # renumber lines and columns to remove empty ones
    lignes = sorted(set(a for a, _ in G))
    colonnes = sorted(set(b for _, b in G))
    d_lig = {a: i for i, a in enumerate(lignes)}
    d_col = {b: i for i, b in enumerate(colonnes)}

    horizontal = sorted((d_lig[a], d_col[b]) for a, b in horizontal)
    vertical = sorted((d_lig[a], d_col[b]) for a, b in vertical)

    def standard(ab):
        a, b = ab
        return (d_lig[a], d_col[b])

    G.relabel(standard)

    # setting the positions
    positions = {ab: ab for ab in G}
    G.set_pos(positions)
    return G, (horizontal, vertical)
