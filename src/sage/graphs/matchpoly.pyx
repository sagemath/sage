# cython: binding=True
"""
Matching polynomial

This module contains the following methods:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`matching_polynomial` | Computes the matching polynomial of a given graph
    :meth:`complete_poly` | Compute the matching polynomial of the complete graph on `n` vertices.

AUTHORS:

- Robert Miller, Tom Boothby - original implementation

REFERENCE:

[God1993]_


Methods
-------
"""

# ****************************************************************************
#       Copyright (C) 2010 Robert Miller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport check_allocarray, sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer

from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *

x = polygen(ZZ, 'x')


def matching_polynomial(G, complement=True, name=None):
    r"""
    Computes the matching polynomial of the graph `G`.

    If `p(G, k)` denotes the number of `k`-matchings (matchings with `k` edges)
    in `G`, then the matching polynomial is defined as [God1993]_:

    .. MATH::

        \mu(x)=\sum_{k \geq 0} (-1)^k p(G,k) x^{n-2k}

    INPUT:

    - ``complement`` - (default: ``True``) whether to use Godsil's duality
      theorem to compute the matching polynomial from that of the graphs
      complement (see ALGORITHM).

    - ``name`` - optional string for the variable name in the polynomial

    .. NOTE::

        The ``complement`` option uses matching polynomials of complete graphs,
        which are cached. So if you are crazy enough to try computing the
        matching polynomial on a graph with millions of vertices, you might not
        want to use this option, since it will end up caching millions of
        polynomials of degree in the millions.

    ALGORITHM:

    The algorithm used is a recursive one, based on the following observation
    [God1993]_:

    - If `e` is an edge of `G`, `G'` is the result of deleting the edge `e`, and
      `G''` is the result of deleting each vertex in `e`, then the matching
      polynomial of `G` is equal to that of `G'` minus that of `G''`.

      (the algorithm actually computes the *signless* matching polynomial, for
      which the recursion is the same when one replaces the subtraction by an
      addition. It is then converted into the matching polynomial and returned)

    Depending on the value of ``complement``, Godsil's duality theorem
    [God1993]_ can also be used to compute `\mu(x)` :

    .. MATH::

        \mu(\overline{G}, x) = \sum_{k \geq 0} p(G,k) \mu( K_{n-2k}, x)


    Where `\overline{G}` is the complement of `G`, and `K_n` the complete graph
    on `n` vertices.

    EXAMPLES::

        sage: g = graphs.PetersenGraph()
        sage: g.matching_polynomial()
        x^10 - 15*x^8 + 75*x^6 - 145*x^4 + 90*x^2 - 6
        sage: g.matching_polynomial(complement=False)
        x^10 - 15*x^8 + 75*x^6 - 145*x^4 + 90*x^2 - 6
        sage: g.matching_polynomial(name='tom')
        tom^10 - 15*tom^8 + 75*tom^6 - 145*tom^4 + 90*tom^2 - 6
        sage: g = Graph()
        sage: L = [graphs.RandomGNP(8, .3) for i in range(1, 6)]
        sage: prod([h.matching_polynomial() for h in L]) == sum(L, g).matching_polynomial()  # long time (up to 10s on sage.math, 2011)
        True

    ::

        sage: for i in range(1, 12):  # long time (10s on sage.math, 2011)
        ....:     for t in graphs.trees(i):
        ....:         if t.matching_polynomial() != t.characteristic_polynomial():
        ....:             raise RuntimeError('bug for a tree A of size {0}'.format(i))
        ....:         c = t.complement()
        ....:         if c.matching_polynomial(complement=False) != c.matching_polynomial():
        ....:             raise RuntimeError('bug for a tree B of size {0}'.format(i))

    ::

        sage: from sage.graphs.matchpoly import matching_polynomial
        sage: matching_polynomial(graphs.CompleteGraph(0))
        1
        sage: matching_polynomial(graphs.CompleteGraph(1))
        x
        sage: matching_polynomial(graphs.CompleteGraph(2))
        x^2 - 1
        sage: matching_polynomial(graphs.CompleteGraph(3))
        x^3 - 3*x
        sage: matching_polynomial(graphs.CompleteGraph(4))
        x^4 - 6*x^2 + 3
        sage: matching_polynomial(graphs.CompleteGraph(5))
        x^5 - 10*x^3 + 15*x
        sage: matching_polynomial(graphs.CompleteGraph(6))
        x^6 - 15*x^4 + 45*x^2 - 15
        sage: matching_polynomial(graphs.CompleteGraph(7))
        x^7 - 21*x^5 + 105*x^3 - 105*x
        sage: matching_polynomial(graphs.CompleteGraph(8))
        x^8 - 28*x^6 + 210*x^4 - 420*x^2 + 105
        sage: matching_polynomial(graphs.CompleteGraph(9))
        x^9 - 36*x^7 + 378*x^5 - 1260*x^3 + 945*x
        sage: matching_polynomial(graphs.CompleteGraph(10))
        x^10 - 45*x^8 + 630*x^6 - 3150*x^4 + 4725*x^2 - 945
        sage: matching_polynomial(graphs.CompleteGraph(11))
        x^11 - 55*x^9 + 990*x^7 - 6930*x^5 + 17325*x^3 - 10395*x
        sage: matching_polynomial(graphs.CompleteGraph(12))
        x^12 - 66*x^10 + 1485*x^8 - 13860*x^6 + 51975*x^4 - 62370*x^2 + 10395
        sage: matching_polynomial(graphs.CompleteGraph(13))
        x^13 - 78*x^11 + 2145*x^9 - 25740*x^7 + 135135*x^5 - 270270*x^3 + 135135*x

    ::

        sage: G = Graph({0:[1,2], 1:[2]})
        sage: matching_polynomial(G)
        x^3 - 3*x
        sage: G = Graph({0:[1,2]})
        sage: matching_polynomial(G)
        x^3 - 2*x
        sage: G = Graph({0:[1], 2:[]})
        sage: matching_polynomial(G)
        x^3 - x
        sage: G = Graph({0:[], 1:[], 2:[]})
        sage: matching_polynomial(G)
        x^3

    ::

        sage: matching_polynomial(graphs.CompleteGraph(0), complement=False)
        1
        sage: matching_polynomial(graphs.CompleteGraph(1), complement=False)
        x
        sage: matching_polynomial(graphs.CompleteGraph(2), complement=False)
        x^2 - 1
        sage: matching_polynomial(graphs.CompleteGraph(3), complement=False)
        x^3 - 3*x
        sage: matching_polynomial(graphs.CompleteGraph(4), complement=False)
        x^4 - 6*x^2 + 3
        sage: matching_polynomial(graphs.CompleteGraph(5), complement=False)
        x^5 - 10*x^3 + 15*x
        sage: matching_polynomial(graphs.CompleteGraph(6), complement=False)
        x^6 - 15*x^4 + 45*x^2 - 15
        sage: matching_polynomial(graphs.CompleteGraph(7), complement=False)
        x^7 - 21*x^5 + 105*x^3 - 105*x
        sage: matching_polynomial(graphs.CompleteGraph(8), complement=False)
        x^8 - 28*x^6 + 210*x^4 - 420*x^2 + 105
        sage: matching_polynomial(graphs.CompleteGraph(9), complement=False)
        x^9 - 36*x^7 + 378*x^5 - 1260*x^3 + 945*x
        sage: matching_polynomial(graphs.CompleteGraph(10), complement=False)
        x^10 - 45*x^8 + 630*x^6 - 3150*x^4 + 4725*x^2 - 945
        sage: matching_polynomial(graphs.CompleteGraph(11), complement=False)
        x^11 - 55*x^9 + 990*x^7 - 6930*x^5 + 17325*x^3 - 10395*x
        sage: matching_polynomial(graphs.CompleteGraph(12), complement=False)
        x^12 - 66*x^10 + 1485*x^8 - 13860*x^6 + 51975*x^4 - 62370*x^2 + 10395
        sage: matching_polynomial(graphs.CompleteGraph(13), complement=False)
        x^13 - 78*x^11 + 2145*x^9 - 25740*x^7 + 135135*x^5 - 270270*x^3 + 135135*x

    TESTS:

    Non-integer labels should work, (:trac:`15545`)::

        sage: G = Graph(10)
        sage: G.add_vertex((0,1))
        sage: G.add_vertex('X')
        sage: G.matching_polynomial()
        x^12
    """
    if G.has_multiple_edges():
        raise NotImplementedError

    cdef int i, j, d
    cdef fmpz_poly_t pol
    cdef nverts = G.num_verts()

    # Using Godsil's duality theorem when the graph is dense

    if complement and G.density() > 0.5:  # this cutoff could probably be tuned
        f_comp = matching_polynomial(G.complement()).list()
        f = x.parent().zero()
        for i from 0 <= i <= nverts / 2:  # implicit floor
            j = nverts - 2 * i
            f += complete_poly(j) * f_comp[j] * (-1)**i
        return f

    cdef int nedges = G.num_edges()

    # Relabelling the vertices of the graph as [0...n-1] so that they are sorted
    # in increasing order of degree

    cdef list L = []
    for v, d in G.degree_iterator(labels=True):
        L.append((d, v))
    L.sort(key=lambda pair: pair[0])
    G = G.relabel(perm={L[i][1]: i for i in range(nverts)}, inplace=False)
    G.allow_loops(False)

    # Initialization of pol, edges* variables.

    # The edges_mem table is of size (2 * nedges * nedges), and is to be read as
    # nedges blocks of size (2 * nedges). These blocks of size (2 * nedges) are
    # themselves to be read as two blocks of length nedges, addressed as edges1
    # and edges2.
    #
    # Only the first block of size (2 * nedges) is here filled. The function
    # delete_and_add will need the rest of the memory.

    fmpz_poly_init(pol)  # sets to zero
    cdef int** edges = <int**>check_allocarray(2 * nedges, sizeof(int*))
    cdef int* edges_mem = <int*>check_allocarray(2 * nedges * nedges, sizeof(int))

    for i in range(2 * nedges):
        edges[i] = edges_mem + i * nedges

    cdef int* edges1 = edges_mem           # edges[0]
    cdef int* edges2 = edges_mem + nedges  # edges[1]

    cdef int cur = 0
    for i, j in sorted(map(sorted, G.edge_iterator(labels=False))):
        edges1[cur] = i
        edges2[cur] = j
        cur += 1

    # Computing the signless matching polynomial

    sig_on()
    delete_and_add(edges, nverts, nedges, nverts, 0, pol)
    sig_off()

    # Building the actual matching polynomial

    cdef list coeffs_ZZ = []
    cdef Integer c_ZZ
    for i in range(nverts + 1):
        c_ZZ = Integer(0)
        fmpz_poly_get_coeff_mpz(c_ZZ.value, pol, i)
        coeffs_ZZ.append(c_ZZ * (-1)**((nverts - i) // 2))

    f = x.parent()(coeffs_ZZ)
    sig_free(edges_mem)
    sig_free(edges)
    fmpz_poly_clear(pol)
    if name is not None:
        return f.change_variable_name(name)
    return f


# The following is a cache of complete graph matching polynomials.

cdef list complete_matching_polys = [x.parent().one(), x]


def complete_poly(n):
    """
    Compute the matching polynomial of the complete graph on `n` vertices.

    INPUT:

    - ``n`` -- order of the complete graph

    .. TODO::

        This code could probably be made more efficient by using FLINT
        polynomials and being written in Cython, using an array of
        fmpz_poly_t pointers or something...  Right now just about the
        whole complement optimization is written in Python, and could
        be easily sped up.

    EXAMPLES::

        sage: from sage.graphs.matchpoly import complete_poly
        sage: f = complete_poly(10)
        sage: f
        x^10 - 45*x^8 + 630*x^6 - 3150*x^4 + 4725*x^2 - 945
        sage: f = complete_poly(20)
        sage: f[8]
        1309458150
        sage: f = complete_poly(1000)
        sage: len(str(f))
        406824

    TESTS:

    Checking the numerical results up to 20::

        sage: from sage.functions.orthogonal_polys import hermite                       # needs sage.symbolic
        sage: p = lambda n: 2^(-n/2)*hermite(n, x/sqrt(2))                              # needs sage.symbolic
        sage: all(p(i) == complete_poly(i) for i in range(2, 20))                       # needs sage.symbolic
        True
    """
    # global complete_matching_polys # if we do eventually make it a C array...
    cdef int lcmp = len(complete_matching_polys)

    if n < lcmp:
        return complete_matching_polys[n]
    lcmp -= 2
    a = complete_matching_polys[lcmp]
    b = complete_matching_polys[lcmp + 1]
    while lcmp + 2 <= n:
        a, b, lcmp = b, x * b - (lcmp + 1) * a, lcmp + 1
        complete_matching_polys.append(b)
    return b


cdef void delete_and_add(int **edges, int nverts, int nedges, int totverts, int depth, fmpz_poly_t pol):
    """
    Add matching polynomial to pol via recursion.

    NOTE : at the end of this function, pol represents the *SIGNLESS*
    matching polynomial.
    """
    cdef fmpz* coeff

    if nverts == 3:
        coeff = fmpz_poly_get_coeff_ptr(pol, 3)
        if coeff is NULL:
            fmpz_poly_set_coeff_ui(pol, 3, 1)
        else:
            fmpz_add_ui(coeff, coeff, 1)
        coeff = fmpz_poly_get_coeff_ptr(pol, 1)
        if coeff is NULL:
            fmpz_poly_set_coeff_ui(pol, 1, nedges)
        else:
            fmpz_add_ui(coeff, coeff, nedges)
        return

    if not nedges:
        coeff = fmpz_poly_get_coeff_ptr(pol, nverts)
        if coeff is NULL:
            fmpz_poly_set_coeff_ui(pol, nverts, 1)
        else:
            fmpz_add_ui(coeff, coeff, 1)
        return

    cdef int* edges1 = edges[2 * depth]
    cdef int* edges2 = edges[2 * depth + 1]
    cdef int* new_edges1 = edges[2 * depth + 2]
    cdef int* new_edges2 = edges[2 * depth + 3]

    nedges -= 1

    # The last edge is (edge1, edge2)
    cdef int edge1 = edges1[nedges]
    cdef int edge2 = edges2[nedges]
    cdef int new_nedges = 0
    cdef int i, new_edge1, new_edge2

    # The new edges are all the edges that are not incident with (edge1, edge2)

    for i in range(nedges):
        if edge1 == edges1[i]:
            # Since the rest of the edges are incident to edge1 (the edges are
            # sorted by increasing order of their first component), we can stop
            break

        if edge1 != edges2[i] and edge2 != edges2[i]:
            # since edge2 > edge1 > edges1[i], we only need to check edges2[i]
            new_edges1[new_nedges] = edges1[i]
            new_edges2[new_nedges] = edges2[i]
            new_nedges += 1

    delete_and_add(edges, nverts - 2, new_nedges, totverts, depth + 1, pol)
    delete_and_add(edges, nverts, nedges, totverts, depth, pol)
