# cython: binding=True
"""
Chromatic Polynomial

AUTHORS:

- Gordon Royle - original C implementation
- Robert Miller - transplant

REFERENCE:

    Ronald C Read, An improved method for computing the chromatic polynomials of
    sparse graphs.
"""

#*****************************************************************************
#       Copyright (C) 2008 Robert Miller
#       Copyright (C) 2008 Gordon Royle
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_check

from sage.libs.gmp.mpz cimport *
from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer
from sage.ext.memory_allocator cimport MemoryAllocator
from sage.misc.all import prod


def chromatic_polynomial(G, return_tree_basis=False):
    """
    Compute the chromatic polynomial of the graph G.

    The algorithm used is a recursive one, based on the following observations
    of Read:

        - The chromatic polynomial of a tree on n vertices is x(x-1)^(n-1).

        - If e is an edge of G, G' is the result of deleting the edge e, and G''
          is the result of contracting e, then the chromatic polynomial of G is
          equal to that of G' minus that of G''.

    EXAMPLES::

        sage: graphs.CycleGraph(4).chromatic_polynomial()
        x^4 - 4*x^3 + 6*x^2 - 3*x
        sage: graphs.CycleGraph(3).chromatic_polynomial()
        x^3 - 3*x^2 + 2*x
        sage: graphs.CubeGraph(3).chromatic_polynomial()
        x^8 - 12*x^7 + 66*x^6 - 214*x^5 + 441*x^4 - 572*x^3 + 423*x^2 - 133*x
        sage: graphs.PetersenGraph().chromatic_polynomial()
        x^10 - 15*x^9 + 105*x^8 - 455*x^7 + 1353*x^6 - 2861*x^5 + 4275*x^4 - 4305*x^3 + 2606*x^2 - 704*x
        sage: graphs.CompleteBipartiteGraph(3,3).chromatic_polynomial()
        x^6 - 9*x^5 + 36*x^4 - 75*x^3 + 78*x^2 - 31*x
        sage: for i in range(2,7):
        ....:     graphs.CompleteGraph(i).chromatic_polynomial().factor()
        (x - 1) * x
        (x - 2) * (x - 1) * x
        (x - 3) * (x - 2) * (x - 1) * x
        (x - 4) * (x - 3) * (x - 2) * (x - 1) * x
        (x - 5) * (x - 4) * (x - 3) * (x - 2) * (x - 1) * x
        sage: graphs.CycleGraph(5).chromatic_polynomial().factor()
        (x - 2) * (x - 1) * x * (x^2 - 2*x + 2)
        sage: graphs.OctahedralGraph().chromatic_polynomial().factor()
        (x - 2) * (x - 1) * x * (x^3 - 9*x^2 + 29*x - 32)
        sage: graphs.WheelGraph(5).chromatic_polynomial().factor()
        (x - 2) * (x - 1) * x * (x^2 - 5*x + 7)
        sage: graphs.WheelGraph(6).chromatic_polynomial().factor()
        (x - 3) * (x - 2) * (x - 1) * x * (x^2 - 4*x + 5)
        sage: C(x)=graphs.LCFGraph(24, [12,7,-7], 8).chromatic_polynomial()  # long time (6s on sage.math, 2011)
        sage: C(2)  # long time
        0

    By definition, the chromatic number of a graph G is the least integer k such that
    the chromatic polynomial of G is strictly positive at k::

        sage: G = graphs.PetersenGraph()
        sage: P = G.chromatic_polynomial()
        sage: min(i for i in range(11) if P(i) > 0) == G.chromatic_number()
        True

        sage: G = graphs.RandomGNP(10,0.7)
        sage: P = G.chromatic_polynomial()
        sage: min(i for i in range(11) if P(i) > 0) == G.chromatic_number()
        True

    TESTS:

    Check that :trac:`21502` is solved::

        sage: graphs.EmptyGraph().chromatic_polynomial()
        1

    Check that :trac:`27966` is solved::

        sage: Graph([[1, 1]], multiedges=True, loops=True).chromatic_polynomial()
        0
    """
    if not G:
        R = ZZ['x']
        return R.one()
    if G.has_loops():
        R = ZZ['x']
        return R.zero()
    if not G.is_connected():
        return prod([chromatic_polynomial(g) for g in G.connected_components_subgraphs()])
    R = ZZ['x']
    x = R.gen()
    if G.is_tree():
        return x * (x - 1) ** (G.num_verts() - 1)

    cdef int nverts, nedges, i, j, u, v, top, bot, num_chords, next_v
    cdef int *queue
    cdef int *chords1
    cdef int *chords2
    cdef int *bfs_reorder
    cdef int *parent
    cdef mpz_t m, coeff
    cdef mpz_t *tot
    cdef mpz_t *coeffs
    G = G.relabel(inplace=False)
    G.remove_multiple_edges()
    G.remove_loops()
    nverts = G.num_verts()
    nedges = G.num_edges()

    cdef MemoryAllocator mem = MemoryAllocator()
    queue       = <int *>   mem.allocarray(nverts, sizeof(int))
    chords1     = <int *>   mem.allocarray((nedges - nverts + 1), sizeof(int))
    chords2     = <int *>   mem.allocarray((nedges - nverts + 1), sizeof(int))
    parent      = <int *>   mem.allocarray(nverts, sizeof(int))
    bfs_reorder = <int *>   mem.allocarray(nverts, sizeof(int))
    tot         = <mpz_t *> mem.allocarray((nverts+1), sizeof(mpz_t))
    coeffs      = <mpz_t *> mem.allocarray((nverts+1), sizeof(mpz_t))
    num_chords = 0

    # Breadth first search from 0:
    bfs_reorder[0] = 0
    mpz_init(tot[0]) # sets to 0
    for i from 0 < i < nverts:
        bfs_reorder[i] = -1
        mpz_init(tot[i]) # sets to 0
    mpz_init(tot[nverts]) # sets to 0
    queue[0] = 0
    top = 1
    bot = 0
    next_v = 1
    while top > bot:
        v = queue[bot]
        bot += 1
        for u in G.neighbor_iterator(v):
            if bfs_reorder[u] == -1: # if u is not yet in tree
                bfs_reorder[u] = next_v
                next_v += 1
                queue[top] = u
                top += 1
                parent[bfs_reorder[u]] = bfs_reorder[v]
            else:
                if bfs_reorder[u] > bfs_reorder[v]:
                    chords1[num_chords] = bfs_reorder[u]
                    chords2[num_chords] = bfs_reorder[v]
                else:
                    continue
                i = num_chords
                num_chords += 1
                # bubble sort the chords
                while i > 0:
                    if chords1[i-1] > chords1[i]:
                        break
                    if chords1[i-1] == chords1[i] and chords2[i-1] > chords2[i]:
                        break
                    j = chords1[i-1]
                    chords1[i-1] = chords1[i]
                    chords1[i] = j
                    j = chords2[i-1]
                    chords2[i-1] = chords2[i]
                    chords2[i] = j
                    i -= 1
    try:
        contract_and_count(chords1, chords2, num_chords, nverts, tot, parent)
    except BaseException:
        for i in range(nverts):
            mpz_clear(tot[i])
        raise
    for i from 0 <= i <= nverts:
        mpz_init(coeffs[i]) # also sets them to 0
    mpz_init(coeff)
    mpz_init_set_si(m, -1)
    # start with the zero polynomial: f(x) = 0
    for i from nverts >= i > 0:
        if not mpz_sgn(tot[i]):
            continue
        mpz_neg(m, m)

        # do this:
        # f += tot[i]*m*x*(x-1)**(i-1)
        mpz_addmul(coeffs[i], m, tot[i])
        mpz_set_si(coeff, 1)
        for j from 1 <= j < i:
            # an iterative method for binomial coefficients...
            mpz_mul_si(coeff, coeff, j-i)
            mpz_divexact_ui(coeff, coeff, j)
            # coeffs[i-j] += tot[i]*m*coeff
            mpz_mul(coeff, coeff, m)
            mpz_addmul(coeffs[i-j], coeff, tot[i])
            mpz_mul(coeff, coeff, m)
    coeffs_ZZ = []
    cdef Integer c_ZZ
    for i from 0 <= i <= nverts:
        c_ZZ = Integer(0)
        mpz_set(c_ZZ.value, coeffs[i])
        coeffs_ZZ.append(c_ZZ)
    f = R(coeffs_ZZ)

    for i from 0 <= i <= nverts:
        mpz_clear(tot[i])
        mpz_clear(coeffs[i])

    mpz_clear(coeff)
    mpz_clear(m)

    return f


cdef int contract_and_count(int *chords1, int *chords2, int num_chords, int nverts,
                         mpz_t *tot, int *parent) except -1:
    if num_chords == 0:
        mpz_add_ui(tot[nverts], tot[nverts], 1)
        return 0
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *new_chords1 = <int *> mem.allocarray(num_chords, sizeof(int))
    cdef int *new_chords2 = <int *> mem.allocarray(num_chords, sizeof(int))
    cdef int *ins_list1   = <int *> mem.allocarray(num_chords, sizeof(int))
    cdef int *ins_list2   = <int *> mem.allocarray(num_chords, sizeof(int))
    cdef int i, j, k, x1, xj, z, num, insnum, parent_checked
    for i in range(num_chords):
        sig_check()

        # contract chord i, and recurse
        z = chords1[i]
        x1 = chords2[i]
        j = i + 1
        insnum = 0
        parent_checked = 0
        while j < num_chords and chords1[j] == z:
            xj = chords2[j]
            if parent[z] > xj:
                parent_checked = 1
                # now try adding {x1, parent[z]} to the list
                if not parent[x1] == parent[z]:
                    if x1 > parent[z]:
                        ins_list1[insnum] = x1
                        ins_list2[insnum] = parent[z]
                    else:
                        ins_list1[insnum] = parent[z]
                        ins_list2[insnum] = x1
                    insnum += 1
            if not parent[x1] == xj: # then {x1, xj} isn't already a tree edge
                ins_list1[insnum] = x1
                ins_list2[insnum] = xj
                insnum += 1
            j += 1
        if not parent_checked:
            if not parent[x1] == parent[z]:
                if x1 > parent[z]:
                    ins_list1[insnum] = x1
                    ins_list2[insnum] = parent[z]
                else:
                    ins_list1[insnum] = parent[z]
                    ins_list2[insnum] = x1
                insnum += 1

        # now merge new_chords and ins_list
        num = 0
        k = 0
        while k < insnum and j < num_chords:
            if chords1[j] > ins_list1[k] or \
              (chords1[j] == ins_list1[k] and chords2[j] > ins_list2[k]):
                new_chords1[num] = chords1[j]
                new_chords2[num] = chords2[j]
                num += 1
                j += 1
            elif chords1[j] < ins_list1[k] or \
              (chords1[j] == ins_list1[k] and chords2[j] < ins_list2[k]):
                new_chords1[num] = ins_list1[k]
                new_chords2[num] = ins_list2[k]
                num += 1
                k += 1
            else:
                new_chords1[num] = chords1[j]
                new_chords2[num] = chords2[j]
                num += 1
                j += 1
                k += 1
        if j == num_chords:
            while k < insnum:
                new_chords1[num] = ins_list1[k]
                new_chords2[num] = ins_list2[k]
                num += 1
                k += 1
        elif k == insnum:
            while j < num_chords:
                new_chords1[num] = chords1[j]
                new_chords2[num] = chords2[j]
                num += 1
                j += 1
        contract_and_count(new_chords1, new_chords2, num, nverts - 1, tot, parent)
    mpz_add_ui(tot[nverts], tot[nverts], 1)
