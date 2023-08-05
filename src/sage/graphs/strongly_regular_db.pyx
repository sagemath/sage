# -*- coding: utf-8 -*-
r"""
Database of strongly regular graphs

This module manages a database associating to a set of four integers
`(v,k,\lambda,\mu)` a strongly regular graphs with these parameters, when one
exists.

Using Andries Brouwer's `database of strongly regular graphs
<https://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__, it can also return
non-existence results. Note that some constructions are missing, and that some
strongly regular graphs that exist in the database cannot be automatically built
by Sage. Help us if you know any.
An outline of the implementation can be found in [CP2016]_.

.. NOTE::

    Any missing/incorrect information in the database must be reported to
    `Andries E. Brouwer <https://www.win.tue.nl/~aeb/>`__ directly, in order to
    have a unique and updated source of information.

REFERENCES:

[BL1984]_

Functions
---------
"""

import json
import os

from libc.math cimport sqrt, floor
from libc.stdint cimport uint_fast32_t

from sage.arith.misc import divisors, is_prime_power, is_square
from sage.categories.sets_cat import EmptySetError
from sage.graphs.graph import Graph
from sage.misc.cachefunc import cached_function
from sage.misc.lazy_import import LazyImport
from sage.misc.unknown import Unknown
from sage.rings.sum_of_squares cimport two_squares_c

orthogonal_array = LazyImport('sage.combinat.designs.orthogonal_arrays', 'orthogonal_array')
balanced_incomplete_block_design = LazyImport('sage.combinat.designs.bibd', 'balanced_incomplete_block_design')
GF = LazyImport('sage.rings.finite_rings.finite_field_constructor', 'GF')
Matrix = LazyImport('sage.matrix.constructor', 'Matrix')
LinearCode = LazyImport('sage.coding.linear_code', 'LinearCode')

cdef dict _brouwer_database = None
_small_srg_database = None


@cached_function
def is_paley(int v, int k, int l, int mu):
    r"""
    Test whether some Paley graph is `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_paley
        sage: t = is_paley(13,6,2,3); t
        (..., 13)
        sage: g = t[0](*t[1:]); g                                                       # needs sage.rings.finite_rings
        Paley graph with parameter 13: Graph on 13 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.rings.finite_rings
        (13, 6, 2, 3)
        sage: t = is_paley(5,5,5,5); t
    """
    if (v % 4 == 1 and is_prime_power(v) and
            k == (v - 1)//2 and
            l == (v - 5)//4 and
            mu == (v - 1)//4):
        from sage.graphs.generators.families import PaleyGraph
        return (PaleyGraph, v)


@cached_function
def is_mathon_PC_srg(int v, int k, int l, int mu):
    r"""
    Test whether some Mathon's Pseudocyclic s.r.g. is `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    .. TODO::

        The current implementation only gives a subset of all possible graphs that can be
        obtained using this construction. A  full implementation should rely on a database
        of conference matrices (or, equivalently, on a database of s.r.g.'s with parameters
        `(4t+1,2t,t-1,t)`. Currently we make an extra assumption that `4t+1` is a prime power.
        The first case where we miss a construction is `t=11`, where we could (recursively)
        use the graph for `t=1` to construct a graph on 83205 vertices.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_mathon_PC_srg
        sage: t = is_mathon_PC_srg(45,22,10,11); t                                      # needs sage.libs.pari
        (..., 1)
        sage: g = t[0](*t[1:]); g                                                       # needs sage.libs.pari
        Mathon's PC SRG on 45 vertices: Graph on 45 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.pari
        (45, 22, 10, 11)

    TESTS::

        sage: t = is_mathon_PC_srg(5,5,5,5); t                                          # needs sage.libs.pari
        sage: mu = 1895  # t=5 case -- the construction cannot work                     # needs sage.libs.pari
        sage: t = is_mathon_PC_srg(4*mu+1,2*mu,mu-1,mu); t                              # needs sage.libs.pari
    """
    cdef int t
    if (v % 4 == 1 and
            k == (v - 1)//2 and
            l == (v - 5)//4 and
            mu == (v - 1)//4):
        from sage.rings.integer_ring import ZZ
        K = ZZ['x']
        x = K.gen()
        rpoly = (w for w in (x*(4*x*(4*x - 1) - 1) - mu).roots() if w[0] > 0)
        try:
            t = next(rpoly)[0]
            if (is_prime_power(4*t - 1) and
                    is_prime_power(4*t + 1)):  # extra assumption in TODO!
                from sage.graphs.generators.families import \
                                    MathonPseudocyclicStronglyRegularGraph
                return (MathonPseudocyclicStronglyRegularGraph, t)
        except StopIteration:
            pass


@cached_function
def is_muzychuk_S6(int v, int k, int l, int mu):
    r"""
    Test whether some Muzychuk S6 graph is (v, k, l, mu)-strongly regular.

    Tests whether a :func:`~sage.graphs.graph_generators.GraphGenerators.MuzychukS6Graph`
    has parameters (v, k, l, mu).

    INPUT:

    - ``v, k, l, mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the required graph if it exists,
    and ``None`` otherwise.

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: from sage.graphs.strongly_regular_db import is_muzychuk_S6
        sage: t = is_muzychuk_S6(378, 116, 34, 36)
        sage: G = t[0](*t[1:]); G
        Muzychuk S6 graph with parameters (3,3): Graph on 378 vertices
        sage: G.is_strongly_regular(parameters=True)
        (378, 116, 34, 36)
        sage: t = is_muzychuk_S6(5, 5, 5, 5); t
    """
    cdef int n, d
    from sage.rings.integer_ring import ZZ
    n_list = [n for n in range(l - 1) if ZZ(n).is_prime_power()]
    for n in n_list:
        d = 2
        while n**d * ((n**d - 1)//(n - 1) + 1) <= v:
            if (v == n**d * ((n**d - 1)//(n - 1) + 1) and
                    k == n**(d - 1)*(n**d - 1)//(n - 1) - 1 and
                    l == mu - 2 and
                    mu == n**(d - 1) * (n**(d - 1) - 1)//(n - 1)):
                from sage.graphs.generators.families import MuzychukS6Graph
                return (MuzychukS6Graph, n, d)
            d += 1


@cached_function
def is_orthogonal_array_block_graph(int v, int k, int l, int mu):
    r"""
    Test whether some (pseudo)Orthogonal Array graph is `(v,k,\lambda,\mu)`-strongly regular.

    We know how to construct graphs with parameters of an Orthogonal Array (`OA(m,n)`),
    also known as Latin squares graphs `L_m(n)`, in several cases where no orthogonal
    array is known, or even in some cases for which they are known not to exist.

    Such graphs are usually called pseudo-Latin squares graphs. Namely, Sage
    can construct a graph with parameters of an `OA(m,n)`-graph whenever there
    exists a skew-Hadamard matrix of order `n+1`, and `m=(n+1)/2` or
    `m=(n-1)/2`. The construction in the former case is due to Goethals-Seidel
    [BL1984]_, and in the latter case due to Pasechnik [Pas1992]_.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: # needs sage.combinat sage.modules
        sage: from sage.graphs.strongly_regular_db import is_orthogonal_array_block_graph
        sage: t = is_orthogonal_array_block_graph(64, 35, 18, 20); t
        (..., 5, 8)
        sage: g = t[0](*t[1:]); g
        OA(5,8): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 35, 18, 20)
        sage: t = is_orthogonal_array_block_graph(225,98,43,42); t
        (..., 4)
        sage: g = t[0](*t[1:]); g
        Pasechnik Graph_4: Graph on 225 vertices
        sage: g.is_strongly_regular(parameters=True)
        (225, 98, 43, 42)
        sage: t = is_orthogonal_array_block_graph(225,112,55,56); t
        (..., 4)
        sage: g = t[0](*t[1:]); g
        skewhad^2_4: Graph on 225 vertices
        sage: g.is_strongly_regular(parameters=True)
        (225, 112, 55, 56)

        sage: t = is_orthogonal_array_block_graph(5,5,5,5); t                           # needs sage.combinat sage.modules
    """
    # notations from
    # https://www.win.tue.nl/~aeb/graphs/OA.html
    from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix
    try:
        m, n = latin_squares_graph_parameters(v, k, l, mu)
    except Exception:
        return
    if orthogonal_array(m, n, existence=True) is True:
        from sage.graphs.generators.intersection import OrthogonalArrayBlockGraph
        return (lambda m, n: OrthogonalArrayBlockGraph(m, n), m, n)

    elif n > 2 and skew_hadamard_matrix(n+1, existence=True) is True:
        if m == (n + 1)/2:
            from sage.graphs.generators.families import SquaredSkewHadamardMatrixGraph as G
        elif m == (n - 1)//2:
            from sage.graphs.generators.families import PasechnikGraph as G
        else:
            return
        return (G, (n+1)//4)


@cached_function
def is_johnson(int v, int k, int l, int mu):
    r"""
    Test whether some Johnson graph is `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_johnson
        sage: t = is_johnson(10,6,3,4); t
        (..., 5)
        sage: g = t[0](*t[1:]); g
        Johnson graph with parameters 5,2: Graph on 10 vertices
        sage: g.is_strongly_regular(parameters=True)
        (10, 6, 3, 4)

        sage: t = is_johnson(5,5,5,5); t
    """
    # Using notations of https://www.win.tue.nl/~aeb/graphs/Johnson.html
    #
    # J(n,m) has parameters v = m(m – 1)/2, k = 2(m – 2), λ = m – 2, μ = 4.
    m = l + 2
    if (mu == 4 and
            k == 2*(m - 2) and
            v == m*(m - 1)//2):
        from sage.graphs.generators.families import JohnsonGraph
        return (lambda m: JohnsonGraph(m, 2), m)


@cached_function
def is_steiner(int v, int k, int l, int mu):
    r"""
    Test whether some Steiner graph is `(v,k,\lambda,\mu)`-strongly regular.

    A Steiner graph is the intersection graph of a Steiner set system. For more
    information, see https://www.win.tue.nl/~aeb/graphs/S.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_steiner
        sage: t = is_steiner(26,15,8,9); t
        (..., 13, 3)
        sage: g = t[0](*t[1:]); g
        Intersection Graph: Graph on 26 vertices
        sage: g.is_strongly_regular(parameters=True)
        (26, 15, 8, 9)

        sage: t = is_steiner(5,5,5,5); t
    """
    # Using notations from https://www.win.tue.nl/~aeb/graphs/S.html
    #
    # The block graph of a Steiner 2-design S(2,m,n) has parameters:
    # v = n(n-1)/m(m-1), k = m(n-m)/(m-1), λ = (m-1)^2 + (n-1)/(m–1)–2, μ = m^2.
    if mu <= 1 or not is_square(mu):
        return
    m = int(sqrt(mu))
    n = (k*(m - 1))//m + m

    if (v == (n*(n - 1))/(m*(m - 1)) and
            k == m*(n - m)/(m - 1) and
            l == (m - 1)**2 + (n - 1)/(m - 1) - 2 and
            balanced_incomplete_block_design(n, m, existence=True) is True):
        from sage.graphs.generators.intersection import IntersectionGraph
        return (lambda n, m: IntersectionGraph([frozenset(b) for b in balanced_incomplete_block_design(n, m)]), n, m)


@cached_function
def is_affine_polar(int v, int k, int l, int mu):
    r"""
    Test whether some Affine Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see https://www.win.tue.nl/~aeb/graphs/VO.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_affine_polar
        sage: t = is_affine_polar(81,32,13,12); t                                       # needs sage.rings.finite_rings
        (..., 4, 3)
        sage: g = t[0](*t[1:]); g                                                       # needs sage.rings.finite_rings
        Affine Polar Graph VO^+(4,3): Graph on 81 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.rings.finite_rings
        (81, 32, 13, 12)

        sage: t = is_affine_polar(5,5,5,5); t
    """
    # Using notations from https://www.win.tue.nl/~aeb/graphs/VO.html
    #
    # VO+(2e,q) has parameters: v = q^(2e), k = (q^(e−1) + 1)(q^e − 1), λ =
    # q(q^(e−2) + 1)(q^(e−1) − 1) + q − 2, μ = q^(e−1)(q^(e−1) + 1)
    #
    # VO−(2e,q) has parameters v = q^(2e), k = (q^(e−1) - 1)(q^e + 1), λ =
    # q(q^(e−2) - 1)(q^(e−1) + 1) + q − 2, μ = q^(e−1)(q^(e−1) - 1)
    if not is_square(v) or not is_prime_power(v):
        return
    prime, power = is_prime_power(v, get_data=True)
    if power % 2:
        return
    for e in divisors(power/2):
        q = prime**(power//(2*e))
        assert v == q**(2*e)
        if (k == (q**(e - 1) + 1)*(q**e - 1) and
                l == q*(q**(e - 2) + 1)*(q**(e - 1) - 1) + q - 2 and
                mu == q**(e - 1)*(q**(e - 1) + 1)):
            from sage.graphs.generators.classical_geometries import AffineOrthogonalPolarGraph
            return (lambda d, q: AffineOrthogonalPolarGraph(d, q, sign='+'), 2*e, q)
        if (k == (q**(e - 1) - 1)*(q**e + 1) and
                l == q*(q**(e - 2) - 1)*(q**(e - 1) + 1) + q - 2 and
                mu == q**(e - 1)*(q**(e - 1) - 1)):
            from sage.graphs.generators.classical_geometries import AffineOrthogonalPolarGraph
            return (lambda d, q: AffineOrthogonalPolarGraph(d, q, sign='-'), 2*e, q)


@cached_function
def is_orthogonal_polar(int v, int k, int l, int mu):
    r"""
    Test whether some Orthogonal Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see https://www.win.tue.nl/~aeb/graphs/srghub.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_orthogonal_polar
        sage: t = is_orthogonal_polar(85, 20, 3, 5); t
        (<function OrthogonalPolarGraph at ...>, 5, 4, '')
        sage: g = t[0](*t[1:]); g                                                       # needs sage.rings.finite_rings
        Orthogonal Polar Graph O(5, 4): Graph on 85 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.rings.finite_rings
        (85, 20, 3, 5)

        sage: t = is_orthogonal_polar(5,5,5,5); t                                       # needs sage.rings.finite_rings

    TESTS:

    All of ``O(2m+1,q)``, ``O^+(2m,q)`` and ``O^-(2m,q)`` appear::

        sage: is_orthogonal_polar(85, 20, 3, 5)
        (<function OrthogonalPolarGraph at ...>, 5, 4, '')
        sage: is_orthogonal_polar(119,54,21,27)
        (<function OrthogonalPolarGraph at ...>, 8, 2, '-')
        sage: is_orthogonal_polar(130,48,20,16)                                         # needs sage.rings.finite_rings
        (<function OrthogonalPolarGraph at ...>, 6, 3, '+')

    """
    r, s = eigenvalues(v, k, l, mu)
    if r is None:
        return
    q_pow_m_minus_one = -s-1 if abs(s) > r else r+1

    if is_prime_power(q_pow_m_minus_one):
        prime, power = is_prime_power(q_pow_m_minus_one, get_data=True)
        for d in divisors(power):
            q = prime**d
            m = (power//d) + 1

            # O(2m+1,q)
            if (v == (q**(2*m) - 1)//(q - 1) and
                    k == q*(q**(2*m - 2) - 1)//(q - 1) and
                    l == q**2*(q**(2*m - 4) - 1)//(q - 1) + q - 1 and
                    mu == (q**(2*m - 2) - 1)//(q - 1)):
                from sage.graphs.generators.classical_geometries import OrthogonalPolarGraph
                return (OrthogonalPolarGraph, 2*m+1, q, "")

            # O^+(2m,q)
            if (v == (q**(2*m - 1) - 1)//(q - 1) + q**(m - 1) and
                    k == q*(q**(2*m - 3) - 1)//(q - 1) + q**(m - 1) and
                    k == q**(2*m - 3) + l + 1 and
                    mu == k//q):
                from sage.graphs.generators.classical_geometries import OrthogonalPolarGraph
                return (OrthogonalPolarGraph, 2*m, q, "+")

            # O^+(2m+1,q)
            if (v == (q**(2*m - 1) - 1)//(q - 1) - q**(m - 1) and
                    k == q*(q**(2*m - 3) - 1)//(q - 1) - q**(m - 1) and
                    k == q**(2*m - 3) + l + 1 and
                    mu == k//q):
                from sage.graphs.generators.classical_geometries import OrthogonalPolarGraph
                return (OrthogonalPolarGraph, 2*m, q, "-")


@cached_function
def is_goethals_seidel(int v, int k, int l, int mu):
    r"""
    Test whether some
    :func:`~sage.graphs.graph_generators.GraphGenerators.GoethalsSeidelGraph` graph is
    `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_goethals_seidel
        sage: t = is_goethals_seidel(28, 15, 6, 10); t                                  # needs sage.combinat sage.modules
        [<function GoethalsSeidelGraph at ...>, 3, 3]
        sage: g = t[0](*t[1:]); g                                                       # needs sage.combinat sage.modules
        Graph on 28 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.combinat sage.modules
        (28, 15, 6, 10)

        sage: t = is_goethals_seidel(256, 135, 70, 72); t                               # needs sage.combinat sage.modules
        [<function GoethalsSeidelGraph at ...>, 2, 15]
        sage: g = t[0](*t[1:]); g                                                       # needs sage.combinat sage.modules
        Graph on 256 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.combinat sage.modules
        (256, 135, 70, 72)

        sage: t = is_goethals_seidel(5,5,5,5); t                                        # needs sage.combinat sage.modules

    TESTS::

        sage: for p in [(16, 9, 4, 6), (28, 15, 6, 10),                                 # needs sage.combinat sage.modules
        ....:           (64, 35, 18, 20), (120, 63, 30, 36),
        ....:           (144, 77, 40, 42), (256, 135, 70, 72), (400, 209, 108, 110),
        ....:           (496, 255, 126, 136), (540, 275, 130, 150), (576, 299, 154, 156),
        ....:           (780, 399, 198, 210), (784, 405, 208, 210), (976, 495, 238, 264)]:
        ....:     print(is_goethals_seidel(*p))
        [<function GoethalsSeidelGraph at ...>, 2, 3]
        [<function GoethalsSeidelGraph at ...>, 3, 3]
        [<function GoethalsSeidelGraph at ...>, 2, 7]
        [<function GoethalsSeidelGraph at ...>, 3, 7]
        [<function GoethalsSeidelGraph at ...>, 2, 11]
        [<function GoethalsSeidelGraph at ...>, 2, 15]
        [<function GoethalsSeidelGraph at ...>, 2, 19]
        [<function GoethalsSeidelGraph at ...>, 3, 15]
        [<function GoethalsSeidelGraph at ...>, 5, 11]
        [<function GoethalsSeidelGraph at ...>, 2, 23]
        [<function GoethalsSeidelGraph at ...>, 3, 19]
        [<function GoethalsSeidelGraph at ...>, 2, 27]
        [<function GoethalsSeidelGraph at ...>, 5, 15]
    """
    from sage.combinat.designs.bibd import balanced_incomplete_block_design
    from sage.combinat.matrices.hadamard_matrix import hadamard_matrix

    # here we guess the parameters v_bibd,k_bibd and r_bibd of the block design
    #
    # - the number of vertices v is equal to v_bibd*(r_bibd+1)
    # - the degree k of the graph is equal to k=(v+r_bibd-1)/2

    r_bibd = k - (v - 1 - k)
    v_bibd = v//(r_bibd + 1)
    k_bibd = (v_bibd - 1)//r_bibd + 1 if r_bibd > 0 else -1

    if (v == v_bibd*(r_bibd + 1) and
            2*k == v + r_bibd - 1 and
            4*l == -2*v + 6*k - v_bibd - k_bibd and
            hadamard_matrix(r_bibd + 1, existence=True) is True and
            balanced_incomplete_block_design(v_bibd, k_bibd, existence=True) is True):
        from sage.graphs.generators.families import GoethalsSeidelGraph
        return [GoethalsSeidelGraph, k_bibd, r_bibd]


@cached_function
def is_NOodd(int v, int k, int l, int mu):
    r"""
    Test whether some NO^e(2n+1,q) graph is `(v,k,\lambda,\mu)`-strongly regular.

    Here `q>2`, for in the case `q=2` this graph is complete. For more
    information, see
    :func:`sage.graphs.graph_generators.GraphGenerators.NonisotropicOrthogonalPolarGraph`
    and Sect. 7.C of [BL1984]_.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_NOodd
        sage: t = is_NOodd(120, 51, 18, 24); t                                          # needs sage.libs.pari
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 4, '-')
        sage: g = t[0](*t[1:]); g                                                       # needs sage.libs.pari
        NO^-(5, 4): Graph on 120 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.pari
        (120, 51, 18, 24)

    TESTS:

    All of ``NO^+(2m+1,q)`` and ``NO^-(2m+1,q)`` appear::

        sage: # needs sage.libs.pari
        sage: t = is_NOodd(120, 51, 18, 24); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 4, '-')
        sage: t = is_NOodd(136, 75, 42, 40); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 4, '+')
        sage: t = is_NOodd(378, 260, 178, 180); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 7, 3, '+')
        sage: t = is_NOodd(45, 32, 22, 24); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 3, '+')
        sage: t = is_NOodd(351, 224, 142, 144); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 7, 3, '-')
        sage: t = is_NOodd(325, 144, 68, 60); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 5, '+')
        sage: t = is_NOodd(300, 104, 28, 40); t
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 5, '-')
        sage: t = is_NOodd(5,5,5,5); t
    """
    cdef int n, q
    r, s = eigenvalues(v, k, l, mu)  # -eq^(n-1)-1 and eq^(n-1)(q-2)-1; q=3 is special case
    if r is None:
        return
    r += 1
    s += 1
    if abs(r) > abs(s):
        (r, s) = (s, r)  # r=-eq^(n-1) s= eq^(n-1)(q-2)
    q = 2 - s//r
    p, t = is_prime_power(q, get_data=True)
    pp, kk = is_prime_power(abs(r), get_data=True)
    if p == pp and t:
        n = kk//t + 1
        e = 1 if v == (q**n)*(q**n + 1)//2 else -1
        if (v == (q**n)*(q**n + e)//2 and
                k == (q**n - e)*(q**(n - 1) + e) and
                l == 2*(q**(2*n - 2) - 1) + e*q**(n - 1)*(q - 1) and
                mu == 2*q**(n - 1)*(q**(n - 1) + e)):
            from sage.graphs.generators.classical_geometries import NonisotropicOrthogonalPolarGraph
            return (NonisotropicOrthogonalPolarGraph, 2*n + 1, q, '+' if e == 1 else '-')


@cached_function
def is_NOperp_F5(int v, int k, int l, int mu):
    r"""
    Test whether some NO^e,perp(2n+1,5) graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see
    :func:`sage.graphs.graph_generators.GraphGenerators.NonisotropicOrthogonalPolarGraph`
    and Sect. 7.D of [BL1984]_.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_NOperp_F5
        sage: t = is_NOperp_F5(10, 3, 0, 1); t                                          # needs sage.libs.pari
        (<function NonisotropicOrthogonalPolarGraph at ...>, 3, 5, '-', 1)
        sage: g = t[0](*t[1:]); g                                                       # needs sage.libs.pari
        NO^-,perp(3, 5): Graph on 10 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.pari
        (10, 3, 0, 1)

    TESTS:

    All of ``NO^+,perp(2m+1,5)`` and ``NO^-,perp(2m+1,5)`` appear::

        sage: t = is_NOperp_F5(325, 60, 15, 10); t                                      # needs sage.libs.pari
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 5, '+', 1)
        sage: t = is_NOperp_F5(300, 65, 10, 15); t                                      # needs sage.libs.pari
        (<function NonisotropicOrthogonalPolarGraph at ...>, 5, 5, '-', 1)
        sage: t = is_NOperp_F5(5,5,5,5); t                                              # needs sage.libs.pari
    """
    cdef int n
    r, s = eigenvalues(v, k, l, mu)  # 2*e*5**(n-1), -e*5**(n-1); note exceptional case n=1
    if r is None:
        return
    if abs(r) < abs(s):
        (r, s) = (s, r)
    e = 1 if s < 0 else -1
    p, n = is_prime_power(abs(s), get_data=True)
    if (5 == p and n) or (abs(r) == 2 and abs(s) == 1):
        n += 1
        if (v == (5**n)*(5**n + e)//2 and
                k == (5**n - e)*5**(n - 1)//2 and
                l == 5**(n - 1)*(5**(n - 1) + e)//2 and
                mu == 5**(n - 1)*(5**(n - 1) - e)//2):
            from sage.graphs.generators.classical_geometries import NonisotropicOrthogonalPolarGraph
            return (NonisotropicOrthogonalPolarGraph, 2*n + 1, 5, '+' if e == 1 else '-', 1)


@cached_function
def is_NO_F2(int v, int k, int l, int mu):
    r"""
    Test whether some NO^e,perp(2n,2) graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see
    :func:`sage.graphs.graph_generators.GraphGenerators.NonisotropicOrthogonalPolarGraph`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_NO_F2
        sage: t = is_NO_F2(10, 3, 0, 1); t                                              # needs sage.libs.pari
        (<function NonisotropicOrthogonalPolarGraph at ...>, 4, 2, '-')
        sage: g = t[0](*t[1:]); g                                                       # needs sage.libs.pari
        NO^-(4, 2): Graph on 10 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.pari
        (10, 3, 0, 1)

    TESTS:

    All of ``NO^+(2m,2)`` and ``NO^-(2m,2)`` appear::

        sage: t = is_NO_F2(36, 15, 6, 6); t                                             # needs sage.libs.pari
        (<function NonisotropicOrthogonalPolarGraph at ...>, 6, 2, '-')
        sage: t = is_NO_F2(28, 15, 6, 10); t                                            # needs sage.libs.pari
        (<function NonisotropicOrthogonalPolarGraph at ...>, 6, 2, '+')
        sage: t = is_NO_F2(5,5,5,5); t                                                  # needs sage.libs.pari
    """
    cdef int n, e, p
    p, n = is_prime_power(k+1, get_data=True)  # k+1==2**(2*n-2)
    if 2 == p and n and not n % 2:
        n = (n+2)//2
        e = (2**(2*n-1)-v)//2**(n-1)
        if (abs(e) == 1 and
                v == 2**(2*n - 1) - e*2**(n - 1) and
                k == 2**(2*n - 2) - 1 and
                l == 2**(2*n - 3) - 2 and
                mu == 2**(2*n - 3) + e*2**(n - 2)):
            from sage.graphs.generators.classical_geometries import NonisotropicOrthogonalPolarGraph
            return (NonisotropicOrthogonalPolarGraph, 2*n, 2, '+' if e == 1 else '-')


@cached_function
def is_NO_F3(int v, int k, int l, int mu):
    r"""
    Test whether some NO^e,perp(2n,3) graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see
    :func:`sage.graphs.graph_generators.GraphGenerators.NonisotropicOrthogonalPolarGraph`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_NO_F3
        sage: t = is_NO_F3(15, 6, 1, 3); t                                              # needs sage.libs.pari
        (<function NonisotropicOrthogonalPolarGraph at ...>, 4, 3, '-')
        sage: g = t[0](*t[1:]); g                                                       # needs sage.libs.pari
        NO^-(4, 3): Graph on 15 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.pari
        (15, 6, 1, 3)

    TESTS:

    All of ``NO^+(2m,3)`` and ``NO^-(2m,3)`` appear::

        sage: t = is_NO_F3(126, 45, 12, 18); t                                          # needs sage.libs.pari
        (<function NonisotropicOrthogonalPolarGraph at ...>, 6, 3, '-')
        sage: t = is_NO_F3(117, 36, 15, 9); t                                           # needs sage.libs.pari
        (<function NonisotropicOrthogonalPolarGraph at ...>, 6, 3, '+')
        sage: t = is_NO_F3(5,5,5,5); t                                                  # needs sage.libs.pari
    """
    cdef int n, e, p
    r, s = eigenvalues(v, k, l, mu)  # e*3**(n-1), -e*3**(n-2)
    if r is None:
        return
    if abs(r) < abs(s):
        (r, s) = (s, r)
    e = 1 if r > 0 else -1
    p, n = is_prime_power(abs(r), get_data=True)
    if 3 == p and n:
        n += 1
        if (v == 3**(n - 1)*(3**n - e)//2 and
                k == 3**(n - 1)*(3**(n - 1) - e)//2 and
                l == 3**(n - 2)*(3**(n - 1) + e)//2 and
                mu == 3**(n - 1)*(3**(n - 2) - e)//2):
            from sage.graphs.generators.classical_geometries import NonisotropicOrthogonalPolarGraph
            return (NonisotropicOrthogonalPolarGraph, 2*n, 3, '+' if e == 1 else '-')


@cached_function
def is_NU(int v, int k, int l, int mu):
    r"""
    Test whether some NU(n,q)-graph, is `(v,k,\lambda,\mu)`-strongly regular.

    Note that n>2; for n=2 there is no s.r.g. For more information, see
    :func:`sage.graphs.graph_generators.GraphGenerators.NonisotropicUnitaryPolarGraph`
    and series C14 in [Hub1975]_.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_NU
        sage: t = is_NU(40, 27, 18, 18); t                                              # needs sage.libs.pari
        (<function NonisotropicUnitaryPolarGraph at ...>, 4, 2)
        sage: g = t[0](*t[1:]); g                                                       # needs sage.libs.pari
        NU(4, 2): Graph on 40 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.pari
        (40, 27, 18, 18)

    TESTS::

        sage: # needs sage.libs.pari
        sage: t = is_NU(176, 135, 102, 108); t
        (<function NonisotropicUnitaryPolarGraph at ...>, 5, 2)
        sage: t = is_NU(540, 224, 88, 96); t
        (<function NonisotropicUnitaryPolarGraph at ...>, 4, 3)
        sage: t = is_NU(208, 75, 30, 25); t
        (<function NonisotropicUnitaryPolarGraph at ...>, 3, 4)
        sage: t = is_NU(5,5,5,5); t
    """
    cdef int n, q, e               # special cases: n=3 or q=2
    r, s = eigenvalues(v, k, l, mu)  # r,s = eq^{n-2} - 1, -e(q^2-q-1)q^{n-3} - 1, e=(-1)^n
    if r is None:
        return
    r += 1
    s += 1
    if abs(r) > abs(s):
        (r, s) = (s, r)
    p, t = is_prime_power(abs(r), get_data=True)
    if p == 2:  # it can be that q=2, then we'd have r>s now
        pp, kk = is_prime_power(abs(s), get_data=True)
        if pp == 2 and kk > 0:
            (r, s) = (s, r)
            p, t = is_prime_power(abs(r), get_data=True)
    if r == 1:
        return
    kr = k//(r-1)  # eq^{n-1}+1
    e = 1 if kr > 0 else -1
    q = (kr-1)//r
    pp, kk = is_prime_power(q, get_data=True)
    if p == pp and kk:
        n = t//kk + 2
        if (v == q**(n - 1)*(q**n - e)//(q + 1) and
                k == (q**(n - 1) + e)*(q**(n - 2) - e) and
                l == q**(2*n - 5)*(q + 1) - e*q**(n - 2)*(q - 1) - 2 and
                mu == q**(n - 3)*(q + 1)*(q**(n - 2) - e)):
            from sage.graphs.generators.classical_geometries import NonisotropicUnitaryPolarGraph
            return (NonisotropicUnitaryPolarGraph, n, q)


@cached_function
def is_haemers(int v, int k, int l, int mu):
    r"""
    Test whether some HaemersGraph graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see
    :func:`~sage.graphs.graph_generators.GraphGenerators.HaemersGraph`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_haemers
        sage: t = is_haemers(96, 19, 2, 4); t                                           # needs sage.libs.pari
        (<function HaemersGraph at ...>, 4)
        sage: g = t[0](*t[1:]); g                                                       # needs sage.libs.pari
        Haemers(4): Graph on 96 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.pari
        (96, 19, 2, 4)

    TESTS::

        sage: t = is_haemers(5,5,5,5); t                                                # needs sage.libs.pari
    """
    cdef int q, n, p
    p, n = is_prime_power(mu, get_data=True)
    q = mu
    if 2 == p and n:
        if (v == q**2*(q + 2) and
                k == q*(q + 1) - 1 and
                l == q - 2):
            from sage.graphs.generators.classical_geometries import HaemersGraph
            return (HaemersGraph, q)


@cached_function
def is_cossidente_penttila(int v, int k, int l, int mu):
    r"""
    Test whether some CossidentePenttilaGraph graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see
    :func:`~sage.graphs.graph_generators.GraphGenerators.CossidentePenttilaGraph`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_cossidente_penttila
        sage: t =  is_cossidente_penttila(378, 52, 1, 8); t                             # needs sage.libs.pari
        (<function CossidentePenttilaGraph at ...>, 5)
        sage: g = t[0](*t[1:]); g                       # optional - gap_package_design, needs sage.libs.pari
        CossidentePenttila(5): Graph on 378 vertices
        sage: g.is_strongly_regular(parameters=True)    # optional - gap_package_design, needs sage.libs.pari
        (378, 52, 1, 8)

    TESTS::

        sage: t =  is_cossidente_penttila(56,10,0,2); t                                 # needs sage.libs.pari
        (<function CossidentePenttilaGraph at ...>, 3)
        sage: t =  is_cossidente_penttila(1376,150,2,18); t                             # needs sage.libs.pari
        (<function CossidentePenttilaGraph at ...>, 7)
        sage: t = is_cossidente_penttila(5,5,5,5); t                                    # needs sage.libs.pari
    """
    cdef int q, n, p
    q = 2*l + 3
    p, n = is_prime_power(q, get_data=True)
    if 2 < p and n:
        if (v == (q**3 + 1)*(q + 1)//2 and
                k == (q**2 + 1)*(q - 1)//2 and
                mu == (q - 1)**2//2):
            from sage.graphs.generators.classical_geometries import CossidentePenttilaGraph
            return (CossidentePenttilaGraph, q)


@cached_function
def is_complete_multipartite(int v, int k, int l, int mu):
    r"""
    Test whether some complete multipartite graph is `(v,k,\lambda,\mu)`-strongly regular.

    Any complete multipartite graph with parts of the same size is strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_complete_multipartite
        sage: t = is_complete_multipartite(12,8,4,8); t
        (<cyfunction is_complete_multipartite.<locals>.CompleteMultipartiteSRG at ...>,
         3,
         4)
        sage: g = t[0](*t[1:]); g
        Multipartite Graph with set sizes [4, 4, 4]: Graph on 12 vertices
        sage: g.is_strongly_regular(parameters=True)
        (12, 8, 4, 8)

    TESTS::

        sage: t = is_complete_multipartite(5,5,5,5); t
        sage: t = is_complete_multipartite(11,8,4,8); t
        sage: t = is_complete_multipartite(20,16,12,16)
        sage: g = t[0](*t[1:]); g
        Multipartite Graph with set sizes [4, 4, 4, 4, 4]: Graph on 20 vertices
        sage: g.is_strongly_regular(parameters=True)
        (20, 16, 12, 16)
    """
    if v > k:
        r = v//(v - k)  # number of parts (of size v-k each)
        if l == (v - k)*(r - 2) and k == mu and v == r*(v - k):
            from sage.graphs.generators.basic import CompleteMultipartiteGraph

            def CompleteMultipartiteSRG(nparts, partsize):
                return CompleteMultipartiteGraph([partsize] * nparts)
            return (CompleteMultipartiteSRG, r, v - k)


@cached_function
def is_polhill(int v, int k, int l, int mu):
    r"""
    Test whether some graph from [Pol2009]_ is `(1024,k,\lambda,\mu)`-strongly
    regular.

    .. NOTE::

        This function does not actually explore *all* strongly regular graphs
        produced in [Pol2009]_, but only those on 1024 vertices.

        John Polhill offered his help if we attempt to write a code to guess,
        given `(v,k,\lambda,\mu)`, which of his construction must be applied to
        find the graph.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if the
    parameters match, and ``None`` otherwise.

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: from sage.graphs.strongly_regular_db import is_polhill
        sage: t = is_polhill(1024, 231,  38,  56); t
        [<cyfunction is_polhill.<locals>.<lambda> at ...>]
        sage: g = t[0](*t[1:]); g               # not tested (too long)
        Graph on 1024 vertices
        sage: g.is_strongly_regular(parameters=True)    # not tested (too long)
        (1024, 231, 38, 56)
        sage: t = is_polhill(1024, 264,  56,  72); t
        [<cyfunction is_polhill.<locals>.<lambda> at ...>]
        sage: t = is_polhill(1024, 297,  76,  90); t
        [<cyfunction is_polhill.<locals>.<lambda> at ...>]
        sage: t = is_polhill(1024, 330,  98, 110); t
        [<cyfunction is_polhill.<locals>.<lambda> at ...>]
        sage: t = is_polhill(1024, 462, 206, 210); t
        [<cyfunction is_polhill.<locals>.<lambda> at ...>]
    """
    if (v, k, l, mu) not in [(1024, 231,  38,  56),
                             (1024, 264,  56,  72),
                             (1024, 297,  76,  90),
                             (1024, 330,  98, 110),
                             (1024, 462, 206, 210)]:
        return

    from itertools import product
    from sage.categories.cartesian_product import cartesian_product
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    from copy import copy

    def additive_cayley(vertices):
        g = Graph()
        g.add_vertices(vertices[0].parent())
        edges = [(x, x + vv)
                 for vv in set(vertices)
                 for x in g]
        g.add_edges(edges)
        g.relabel()
        return g

    # D is a Partial Difference Set of (Z4)^2, see section 2.
    G = cartesian_product([IntegerModRing(4), IntegerModRing(4)])
    D = [[(2, 0), (0, 1), (0, 3), (1, 1), (3, 3)],
         [(1, 0), (3, 0), (0, 2), (1, 3), (3, 1)],
         [(1, 2), (3, 2), (2, 1), (2, 3), (2, 2)]]
    D = [[G(e) for e in x] for x in D]

    # The K_i are hyperplanes partitionning the nonzero elements of
    # GF(2^s)^2. See section 6.
    s = 3
    G1 = GF(2**s,'x')
    Gp = cartesian_product([G1, G1])
    K = [Gp((x, 1)) for x in G1] + [Gp((1, 0))]
    K = [[x for x in Gp if x[0]*uu + x[1]*vv == 0] for (uu, vv) in K]

    # We now define the P_{i,j}. see section 6.

    P = {}
    P[0,1] = list(xrange((-1) + 1                  , 2**(s-2)+1))
    P[1,1] = list(xrange((-1) + 2**(s-2)+2         , 2**(s-1)+1))
    P[2,1] = list(xrange((-1) + 2**(s-1)+2         , 2**(s-1)+2**(s-2)+1))
    P[3,1] = list(xrange((-1) + 2**(s-1)+2**(s-2)+2, 2**(s)+1))

    P[0,2] = list(xrange((-1) + 2**(s-2)+2         , 2**(s-1)+2))
    P[1,2] = list(xrange((-1) + 2**(s-1)+3         , 2**(s-1)+2**(s-2)+2))
    P[2,2] = list(xrange((-1) + 2**(s-1)+2**(s-2)+3, 2**(s)+1)) + [0]
    P[3,2] = list(xrange((-1) + 2                  , 2**(s-2)+1))

    P[0,3] = list(xrange((-1) + 2**(s-1)+3         , 2**(s-1)+2**(s-2)+3))
    P[1,3] = list(xrange((-1) + 2**(s-1)+2**(s-2)+4, 2**(s)+1)) + [0,1]
    P[2,3] = list(xrange((-1) + 3                  , 2**(s-2)+2))
    P[3,3] = list(xrange((-1) + 2**(s-2)+3         , 2**(s-1)+2))

    P[0,4] = list(xrange((-1) + 2**(s-1)+2**(s-2)+4, 2**(s)+1))
    P[1,4] = list(xrange((-1) + 3                  , 2**(s-2)+1)) + [2**(s-1)+1,2**(s-1)+2**(s-2)+2]
    P[2,4] = list(xrange((-1) + 2**(s-2)+3         , 2**(s-1)+1)) + [2**(s-1)+2**(s-2)+1,1]
    P[3,4] = list(xrange((-1) + 2**(s-1)+3         , 2**(s-1)+2**(s-2)+1)) + [2**(s-2)+1,0]

    R = {x: copy(P[x]) for x in P}

    for x in P:
        P[x] = [K[i] for i in P[x]]
        P[x] = set(sum(P[x], [])).difference([Gp((0, 0))])

    P[1, 4].add(Gp((0, 0)))
    P[2, 4].add(Gp((0, 0)))
    P[3, 4].add(Gp((0, 0)))

    # We now define the R_{i,j}. see *end* of section 6.

    R[0,3] = list(xrange((-1) + 2**(s-1)+3         , 2**(s-1)+2**(s-2)+2))
    R[1,3] = list(xrange((-1) + 2**(s-1)+2**(s-2)+4, 2**(s)+1)) + [0,1,2**(s-1)+2**(s-2)+2]
    R[0,4] = list(xrange((-1) + 2**(s-1)+2**(s-2)+4, 2**(s)+1)) + [2**(s-1)+2**(s-2)+2]
    R[1,4] = list(xrange((-1) + 3                  , 2**(s-2)+1)) + [2**(s-1)+1]

    for x in R:
        R[x] = [K[i] for i in R[x]]
        R[x] = set(sum(R[x], [])).difference([Gp((0, 0))])

    R[1, 3].add(Gp((0, 0)))
    R[2, 4].add(Gp((0, 0)))
    R[3, 4].add(Gp((0, 0)))

    # Dabcd = Da, Db, Dc, Dd (cf. p273)
    # D1234 = D1, D2, D3, D4 (cf. p276)
    Dabcd = []
    D1234 = []

    Gprod = cartesian_product([G, Gp])
    for DD,PQ in [(Dabcd, P), (D1234, R)]:
        for i in range(1, 5):
            Dtmp = [product([G.zero()], PQ[0, i]),
                    product(D[0], PQ[1, i]),
                    product(D[1], PQ[2, i]),
                    product(D[2], PQ[3, i])]
            Dtmp = map(set, Dtmp)
            Dtmp = [Gprod(e) for e in sum(map(list, Dtmp), [])]
            DD.append(Dtmp)

    # Now that we have the data, we can return the graphs.
    if k == 231:
        return [lambda: additive_cayley(Dabcd[0])]
    if k == 264:
        return [lambda: additive_cayley(D1234[2])]
    if k == 297:
        return [lambda: additive_cayley(D1234[0] + D1234[1] + D1234[2]).complement()]
    if k == 330:
        return [lambda: additive_cayley(Dabcd[0] + Dabcd[1] + Dabcd[2]).complement()]
    if k == 462:
        return [lambda: additive_cayley(Dabcd[0] + Dabcd[1])]


def is_RSHCD(int v, int k, int l, int mu):
    r"""
    Test whether some RSHCD graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see :func:`SRG_from_RSHCD`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_RSHCD
        sage: t = is_RSHCD(64,27,10,12); t                                              # needs sage.combinat sage.modules
        [<built-in function SRG_from_RSHCD>, 64, 27, 10, 12]
        sage: g = t[0](*t[1:]); g                                                       # needs sage.combinat sage.modules
        Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.combinat sage.modules
        (64, 27, 10, 12)
    """
    if SRG_from_RSHCD(v, k, l, mu, existence=True) is True:
        return [SRG_from_RSHCD, v, k, l, mu]


def SRG_from_RSHCD(v, k, l, mu, existence=False, check=True):
    r"""
    Return a `(v,k,l,mu)`-strongly regular graph from a RSHCD

    This construction appears in 8.D of [BL1984]_. For more information, see
    :func:`~sage.combinat.matrices.hadamard_matrix.regular_symmetric_hadamard_matrix_with_constant_diagonal`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    - ``existence`` (boolean) -- whether to return a graph or to test if Sage
      can build such a graph.

    - ``check`` (boolean) -- whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to ``True``
      by default.

    EXAMPLES:

    some graphs ::

        sage: from sage.graphs.strongly_regular_db import SRG_from_RSHCD
        sage: SRG_from_RSHCD(784, 0, 14, 38, existence=True)                            # needs sage.combinat sage.modules
        False
        sage: SRG_from_RSHCD(784, 377, 180, 182, existence=True)                        # needs sage.combinat sage.modules
        True
        sage: SRG_from_RSHCD(144, 65, 28, 30)                                           # needs sage.combinat sage.modules
        Graph on 144 vertices

    an example with vertex-transitive automorphism group, found during the
    implementation of the case `v=324` ::

        sage: # long time, needs sage.combinat sage.modules
        sage: G = SRG_from_RSHCD(324,152,70,72)
        sage: a = G.automorphism_group()
        sage: a.order()
        2592
        sage: len(a.orbits())
        1

    TESTS::

        sage: SRG_from_RSHCD(784, 0, 14, 38)                                            # needs sage.combinat sage.modules
        Traceback (most recent call last):
        ...
        ValueError: I do not know how to build a (784, 0, 14, 38)-SRG from a RSHCD
    """
    from sage.combinat.matrices.hadamard_matrix import regular_symmetric_hadamard_matrix_with_constant_diagonal

    def sgn(x):
        return 1 if x >= 0 else -1
    n = v
    a = (n-4*mu)//2
    e = 2*k - n + 1 + a
    t = abs(a//2)

    if (e**2 == 1 and
            k == (n-1-a+e)/2 and
            l == (n-2*a)/4 - (1-e) and
            mu== (n-2*a)/4 and
            regular_symmetric_hadamard_matrix_with_constant_diagonal(n, sgn(a)*e, existence=True) is True):
        if existence:
            return True
        from sage.matrix.constructor import identity_matrix as I
        from sage.matrix.constructor import ones_matrix as J

        H = regular_symmetric_hadamard_matrix_with_constant_diagonal(n, sgn(a)*e)
        if list(H.column(0)[1:]).count(1) == k:
            H = -H
        G = Graph((J(n) - I(n) - H + H[0, 0]*I(n)) / 2,
                  loops=False, multiedges=False, format="adjacency_matrix")
        if check:
            assert G.is_strongly_regular(parameters=True) == (v, k, l, mu)
        return G

    if existence:
        return False
    raise ValueError("I do not know how to build a {}-SRG from a RSHCD".format((v, k, l, mu)))


@cached_function
def is_unitary_polar(int v, int k, int l, int mu):
    r"""
    Test whether some Unitary Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see https://www.win.tue.nl/~aeb/graphs/srghub.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_unitary_polar
        sage: t = is_unitary_polar(45, 12, 3, 3); t                                     # needs sage.libs.pari
        (<function UnitaryPolarGraph at ...>, 4, 2)
        sage: g = t[0](*t[1:]); g                                                       # needs sage.libs.pari
        Unitary Polar Graph U(4, 2); GQ(4, 2): Graph on 45 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.pari
        (45, 12, 3, 3)

        sage: t = is_unitary_polar(5,5,5,5); t                                          # needs sage.libs.pari

    TESTS:

    All the ``U(n,q)`` appear::

        sage: # needs sage.libs.pari
        sage: t = is_unitary_polar(45, 12, 3, 3); t
        (<function UnitaryPolarGraph at ...>, 4, 2)
        sage: t = is_unitary_polar(165, 36, 3, 9); t
        (<function UnitaryPolarGraph at ...>, 5, 2)
        sage: t = is_unitary_polar(693, 180, 51, 45); t
        (<function UnitaryPolarGraph at ...>, 6, 2)
        sage: t = is_unitary_polar(1105, 80, 15, 5); t
        (<function UnitaryPolarGraph at ...>, 4, 4)
    """
    r, s = eigenvalues(v, k, l, mu)
    if r is None:
        return
    q = k//mu
    if q*mu != k or q < 2:
        return
    p, t = is_prime_power(q, get_data=True)
    if p**t != q or t % 2:
        return
    # at this point we know that we should have U(n,q) for some n and q=p^t, t even
    if r > 0:
        q_pow_d_minus_one = r+1
    else:
        q_pow_d_minus_one = -s-1
    ppp, ttt = is_prime_power(q_pow_d_minus_one, get_data=True)
    d = ttt//t + 1
    if ppp != p or (d-1)*t != ttt:
        return
    t //= 2
    # U(2d+1,q); write q^(1/2) as p^t
    if (v == (q**d - 1)*((q**d)*p**t + 1)//(q - 1) and
            k == q*(q**(d-1) - 1)*((q**d)//(p**t) + 1)//(q - 1) and
            l == q*q*(q**(d-2)-1)*((q**(d-1))//(p**t) + 1)//(q - 1) + q - 1):
        from sage.graphs.generators.classical_geometries import UnitaryPolarGraph
        return (UnitaryPolarGraph, 2*d+1, p**t)

    # U(2d,q);
    if (v == (q**d - 1)*((q**d)//(p**t) + 1)//(q - 1) and
            k == q*(q**(d-1) - 1)*((q**(d-1))//(p**t) + 1)//(q - 1) and
            l == q*q*(q**(d-2)-1)*((q**(d-2))//(p**t) + 1)//(q - 1) + q - 1):
        from sage.graphs.generators.classical_geometries import UnitaryPolarGraph
        return (UnitaryPolarGraph, 2*d, p**t)


@cached_function
def is_unitary_dual_polar(int v, int k, int l, int mu):
    r"""
    Test whether some Unitary Dual Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    This must be the U_5(q) on totally isotropic lines.
    For more information, see https://www.win.tue.nl/~aeb/graphs/srghub.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: from sage.graphs.strongly_regular_db import is_unitary_dual_polar
        sage: t = is_unitary_dual_polar(297, 40, 7, 5); t
        (<function UnitaryDualPolarGraph at ...>, 5, 2)
        sage: g = t[0](*t[1:]); g
        Unitary Dual Polar Graph DU(5, 2); GQ(8, 4): Graph on 297 vertices
        sage: g.is_strongly_regular(parameters=True)
        (297, 40, 7, 5)
        sage: t = is_unitary_dual_polar(5,5,5,5); t

    TESTS::

        sage: is_unitary_dual_polar(6832, 270, 26, 10)                                  # needs sage.libs.pari
        (<function UnitaryDualPolarGraph at ...>, 5, 3)
    """
    r, s = eigenvalues(v, k, l, mu)
    if r is None:
        return
    q = mu - 1
    if q < 2:
        return
    p, t = is_prime_power(q, get_data=True)
    if p**t != q or t % 2:
        return
    if (r < 0 and q != -r - 1) or (s < 0 and q != -s - 1):
        return
    t //= 2
    # we have correct mu, negative eigenvalue, and q=p^(2t)
    if (v == (q**2*p**t + 1)*(q*p**t + 1) and
            k == q*p**t*(q + 1) and
            l == k - 1 - q**2*p**t):
        from sage.graphs.generators.classical_geometries import UnitaryDualPolarGraph
        return (UnitaryDualPolarGraph, 5, p**t)


@cached_function
def is_GQqmqp(int v, int k, int l, int mu):
    r"""
    Test whether some `GQ(q-1,q+1)` or `GQ(q+1,q-1)`-graph is `(v,k,\lambda,\mu)`-srg.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: from sage.graphs.strongly_regular_db import is_GQqmqp
        sage: t = is_GQqmqp(27,10,1,5); t
        (<function AhrensSzekeresGeneralizedQuadrangleGraph at ...>, 3, False)
        sage: g = t[0](*t[1:]); g
        AS(3); GQ(2, 4): Graph on 27 vertices
        sage: t = is_GQqmqp(45,12,3,3); t
        (<function AhrensSzekeresGeneralizedQuadrangleGraph at ...>, 3, True)
        sage: g = t[0](*t[1:]); g
        AS(3)*; GQ(4, 2): Graph on 45 vertices
        sage: g.is_strongly_regular(parameters=True)
        (45, 12, 3, 3)
        sage: t = is_GQqmqp(16,6,2,2); t
        (<function T2starGeneralizedQuadrangleGraph at ...>, 2, True)
        sage: g = t[0](*t[1:]); g
        T2*(O,2)*; GQ(3, 1): Graph on 16 vertices
        sage: g.is_strongly_regular(parameters=True)
        (16, 6, 2, 2)
        sage: t = is_GQqmqp(64,18,2,6); t
        (<function T2starGeneralizedQuadrangleGraph at ...>, 4, False)
        sage: g = t[0](*t[1:]); g
        T2*(O,4); GQ(3, 5): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 18, 2, 6)

    TESTS::

        sage: # needs sage.libs.pari
        sage: (S,T) = (127,129)
        sage: t = is_GQqmqp((S+1)*(S*T+1), S*(T+1), S-1, T+1); t
        (<function T2starGeneralizedQuadrangleGraph at ...>, 128, False)
        sage: (S,T) = (129,127)
        sage: t = is_GQqmqp((S+1)*(S*T+1), S*(T+1), S-1, T+1); t
        (<function T2starGeneralizedQuadrangleGraph at ...>, 128, True)
        sage: (S,T) = (124,126)
        sage: t = is_GQqmqp((S+1)*(S*T+1), S*(T+1), S-1, T+1); t
        (<function AhrensSzekeresGeneralizedQuadrangleGraph at ...>, 125, False)
        sage: (S,T) = (126,124)
        sage: t = is_GQqmqp((S+1)*(S*T+1), S*(T+1), S-1, T+1); t
        (<function AhrensSzekeresGeneralizedQuadrangleGraph at ...>, 125, True)
        sage: t = is_GQqmqp(5,5,5,5); t
    """
    # do we have GQ(s,t)? we must have mu=t+1, s=l+1,
    # v=(s+1)(st+1), k=s(t+1)
    S = l + 1
    T = mu - 1
    q = (S+T)//2
    p, w = is_prime_power(q, get_data=True)
    if (v == (S+1)*(S*T+1) and
            k == S*(T+1) and
            q == p**w and
            (S+T)//2 == q):
        if p % 2 == 0:
            from sage.graphs.generators.classical_geometries\
                    import T2starGeneralizedQuadrangleGraph as F
        else:
            from sage.graphs.generators.classical_geometries\
                    import AhrensSzekeresGeneralizedQuadrangleGraph as F
        if (S, T) == (q-1, q+1):
            return (F, q, False)
        elif (S, T) == (q+1, q-1):
            return (F, q, True)


@cached_function
def is_twograph_descendant_of_srg(int v, int k0, int l, int mu):
    r"""
    Test whether some descendant graph of a s.r.g. is `(v,k_0,\lambda,\mu)`-s.r.g.

    We check whether there can exist `(v+1,k,\lambda^*,\mu^*)`-s.r.g. `G` so
    that ``self`` is a descendant graph of the regular two-graph specified
    by `G`.
    Specifically, we must have that `v+1=2(2k-\lambda^*-\mu^*)`, and
    `k_0=2(k-\mu^*)`, `\lambda=k+\lambda^*-2\mu^*`, `\mu=k-\mu^*`, which give 2
    independent linear conditions, say `k-\mu^*=\mu` and
    `\lambda^*-\mu^*=\lambda-\mu`.  Further, there is a quadratic relation
    `2 k^2-(v+1+4 \mu) k+ 2 v \mu=0`.

    If we can construct such `G` then we return a function to build a
    `(v,k_0,\lambda,\mu)`-s.r.g.  For more information,
    see 10.3 in https://www.win.tue.nl/~aeb/2WF02/spectra.pdf

    INPUT:

    - ``v,k0,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists and is known, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_twograph_descendant_of_srg
        sage: t = is_twograph_descendant_of_srg(27, 10, 1, 5); t                        # needs sage.rings.finite_rings
        (<cyfunction is_twograph_descendant_of_srg.<locals>.la at...
        sage: g = t[0](*t[1:]); g                                                       # needs sage.rings.finite_rings
        descendant of complement(Johnson graph with parameters 8,2) at {0, 1}: Graph on 27 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.rings.finite_rings
        (27, 10, 1, 5)
        sage: t = is_twograph_descendant_of_srg(5,5,5,5); t

    TESTS::

        sage: graphs.strongly_regular_graph(279, 150, 85, 75, existence=True)           # needs sage.combinat
        True
        sage: graphs.strongly_regular_graph(279, 150, 85, 75).is_strongly_regular(parameters=True)  # optional - gap_package_design internet
        (279, 150, 85, 75)
    """
    cdef int b, k, s
    if k0 != 2*mu or not v % 2:
        return
    b = v+1+4*mu
    D = sqrt(b**2-16*v*mu)
    if int(D) == D:
        for kf in [(-D+b)//4, (D+b)//4]:
            k = int(kf)
            if (k == kf and
                    strongly_regular_graph(v+1, k, l - 2*mu + k, k - mu,  existence=True) is True):
                try:
                    g = strongly_regular_graph_lazy(v+1, k, l - 2*mu + k)  # Sage might not know how to build g

                    def la(*gr):
                        from sage.combinat.designs.twographs import twograph_descendant
                        gg = g[0](*gr)
                        if (gg.name() is None) or (gg.name() == ''):
                            gg = Graph(gg, name=str((v+1, k, l - 2*mu + k, k - mu))+"-strongly regular graph")
                        return twograph_descendant(gg, next(gg.vertex_iterator()),
                                                   name=True)
                    return (la, *g[1:])
                except RuntimeError:
                    pass
    return


@cached_function
def is_taylor_twograph_srg(int v, int k, int l, int mu):
    r"""
    Test whether some Taylor two-graph SRG is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see §7E of [BL1984]_.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph
    :func:`TaylorTwographSRG
    <sage.graphs.graph_generators.GraphGenerators.TaylorTwographSRG>` if the
    parameters match, and ``None`` otherwise.

    EXAMPLES::

        sage: # needs sage.libs.pari
        sage: from sage.graphs.strongly_regular_db import is_taylor_twograph_srg
        sage: t = is_taylor_twograph_srg(28, 15, 6, 10); t
        (<function TaylorTwographSRG at ...>, 3)
        sage: g = t[0](*t[1:]); g
        Taylor two-graph SRG: Graph on 28 vertices
        sage: g.is_strongly_regular(parameters=True)
        (28, 15, 6, 10)
        sage: t = is_taylor_twograph_srg(5,5,5,5); t

    TESTS::

        sage: is_taylor_twograph_srg(730, 369, 168, 205)                                # needs sage.libs.pari
        (<function TaylorTwographSRG at ...>, 9)

    """
    r, s = eigenvalues(v, k, l, mu)
    if r is None:
        return
    p, t = is_prime_power(v-1, get_data=True)
    if p**t+1 != v or t % 3 != 0 or p % 2 == 0:
        return
    q = p**(t//3)
    if (k, l, mu) == (q*(q**2+1)//2, (q**2+3)*(q-1)//4, (q**2+1)*(q+1)//4):
        from sage.graphs.generators.classical_geometries import TaylorTwographSRG
        return (TaylorTwographSRG, q)
    return


def is_switch_skewhad(int v, int k, int l, int mu):
    r"""
    Test whether some ``switch skewhad^2+*`` is `(v,k,\lambda,\mu)`-strongly regular.

    The ``switch skewhad^2+*`` graphs appear on `Andries Brouwer's database
    <https://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__ and are built by
    adding an isolated vertex to the complement of
    :func:`~sage.graphs.graph_generators.GraphGenerators.SquaredSkewHadamardMatrixGraph`,
    and a :meth:`Seidel switching <Graph.seidel_switching>` a set of disjoint
    `n`-cocliques.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if the
    parameters match, and ``None`` otherwise.

    EXAMPLES::

        sage: graphs.strongly_regular_graph(226, 105, 48, 49)                           # needs sage.combinat sage.modules
        switch skewhad^2+*_4: Graph on 226 vertices

    TESTS::

        sage: from sage.graphs.strongly_regular_db import is_switch_skewhad
        sage: t = is_switch_skewhad(5,5,5,5); t                                         # needs sage.combinat sage.modules

    """
    from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix
    from sage.graphs.generators.families import SwitchedSquaredSkewHadamardMatrixGraph
    cdef int n
    r, s = eigenvalues(v, k, l, mu)
    if r is None:
        return
    if r < s:
        r, s = s, r
    n = -s // 2
    if (int(r) == 2*n-1 and
            v == (4*n-1)**2 + 1 and
            k == (4*n-1)*(2*n-1) and
            skew_hadamard_matrix(4*n, existence=True) is True):
        return (SwitchedSquaredSkewHadamardMatrixGraph, n)


def is_switch_OA_srg(int v, int k, int l, int mu):
    r"""
    Test whether some *switch* `OA(k,n)+*` is `(v,k,\lambda,\mu)`-strongly regular.

    The "switch* `OA(k,n)+*` graphs appear on `Andries Brouwer's database
    <https://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__ and are built by
    adding an isolated vertex to a
    :meth:`~sage.graphs.graph_generators.GraphGenerators.OrthogonalArrayBlockGraph`,
    and a :meth:`Seidel switching <Graph.seidel_switching>` a set of disjoint
    `n`-cocliques.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if the
    parameters match, and ``None`` otherwise.

    EXAMPLES::

        sage: graphs.strongly_regular_graph(170, 78, 35, 36)  # indirect doctest        # needs sage.combinat sage.modules
        Graph on 170 vertices

    TESTS::

        sage: from sage.graphs.strongly_regular_db import is_switch_OA_srg
        sage: t = is_switch_OA_srg(5,5,5,5); t
        sage: t = is_switch_OA_srg(170, 78, 35, 36)                                     # needs sage.schemes
        sage: t[0](*t[1:]).is_strongly_regular(parameters=True)                         # needs sage.schemes
        (170, 78, 35, 36)
        sage: t = is_switch_OA_srg(290, 136,  63,  64)                                  # needs sage.schemes
        sage: t[0](*t[1:]).is_strongly_regular(parameters=True)                         # needs sage.schemes
        (290, 136, 63, 64)
        sage: is_switch_OA_srg(626, 300, 143, 144)                                      # needs sage.schemes
        (<cyfunction is_switch_OA_srg.<locals>.switch_OA_srg at ..., 12, 25)
        sage: is_switch_OA_srg(842, 406, 195, 196)                                      # needs sage.schemes
        (<cyfunction is_switch_OA_srg.<locals>.switch_OA_srg at ..., 14, 29)
    """
    cdef int n_2_p_1 = v
    cdef int n = <int> floor(sqrt(n_2_p_1 - 1))

    if n*n != n_2_p_1 - 1:  # is it a square?
        return None

    cdef int c = k//n
    if (k % n or l != c*c-1 or k != 1+(c-1)*(c+1)+(n-c)*(n-c-1) or
            orthogonal_array(c+1, n, existence=True, resolvable=True) is not True):
        return None

    def switch_OA_srg(c, n):
        OA = [tuple(x) for x in orthogonal_array(c+1, n, resolvable=True)]
        g = Graph([OA, lambda x, y: any(xx == yy for xx, yy in zip(x, y))],
                  loops=False)
        g.add_vertex(0)
        g.seidel_switching(OA[:c*n])
        return g

    return (switch_OA_srg, c, n)


def is_nowhere0_twoweight(int v, int k, int l, int mu):
    r"""
    Test whether some graph of nowhere 0 words is `(v,k,\lambda,\mu)`-strongly regular.

    Test whether a :meth:`~sage.graphs.graph_generators.GraphGenerators.Nowhere0WordsTwoWeightCodeGraph`
    is `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if the
    parameters match, and ``None`` otherwise.

    EXAMPLES::

        sage: graphs.strongly_regular_graph(196, 60, 14, 20)                            # needs sage.combinat sage.modules
        Nowhere0WordsTwoWeightCodeGraph(8): Graph on 196 vertices

    TESTS::

        sage: from sage.graphs.strongly_regular_db import is_nowhere0_twoweight
        sage: t = is_nowhere0_twoweight(1800, 728, 268, 312); t                         # needs sage.libs.pari
        (<function Nowhere0WordsTwoWeightCodeGraph at ...>, 16)
        sage: t = is_nowhere0_twoweight(5,5,5,5); t                                     # needs sage.libs.pari

    """
    from sage.graphs.generators.classical_geometries import Nowhere0WordsTwoWeightCodeGraph
    cdef int q
    r, s = eigenvalues(v, k, l, mu)
    if r is None:
        return
    if r < s:
        r, s = s, r
    q = r*2
    if (q > 4 and is_prime_power(q) and not r % 2 and
            v == r*(q-1)**2 and
            4*k == q*(q-2)*(q-3) and
            8*mu == q*(q-3)*(q-4)):
        return (Nowhere0WordsTwoWeightCodeGraph, q)


cdef eigenvalues(int v, int k, int l, int mu):
    r"""
    Return the eigenvalues of a (v,k,l,mu)-strongly regular graph.

    If the set of parameters is not feasible, or if they correspond to a
    conference graph, the function returns ``(None,None)``. Otherwise
    it returns the pair [r,s] of eigenvalues, satisfying r>s.

    INPUT:

    - ``v,k,l,mu`` (integers)

    """
    # See 1.3.1 of [Distance-regular graphs]
    b = (mu-l)
    c = (mu-k)
    D = b**2-4*c
    if not is_square(D):
        return [None, None]
    return [(-b+sqrt(D))/2.0,
            (-b-sqrt(D))/2.0]


def eigenmatrix(int v, int k, int l, int mu):
    r"""
    Return the first eigenmatrix of a `(v,k,l,mu)`-strongly regular graph.

    The adjacency matrix `A` of an s.r.g. commutes with the adjacency matrix
    `A'=J-A-I` of its complement (here `J` is all-1 matrix, and `I` the identity
    matrix).  Thus, they can be simultaneously diagonalized and so `A` and `A'`
    share eigenspaces.

    The eigenvalues of `J` are `v` with multiplicity 1, and 0 with multiplicity
    `v-1`. Thus the eigenvalue of `A'` corresponding to the 1-dimension
    `k`-eigenspace of `A` is `v-k-1`.  Respectively, the eigenvalues of `A'`
    corresponding to `t`-eigenspace of `A`, with `t` unequal to `k`, equals
    `-t-1`. The 1st eigenmatrix `P` of the C-algebra `C[A]` generated by `A`
    encodes this eigenvalue information in its three columns;
    the 2nd (resp. 3rd)
    column contains distinct eigenvalues of `A` (resp. of `A'`), and the 1st
    column contains the corresponding eigenvalues of `I`. The matrix `vP^{-1}`
    is called the 2nd eigenvalue matrix of `C[A]`.

    The most interesting feature of `vP^{-1}` is that it is the 1st eigenmatrix
    of the dual of `C[A]` if the dual is generated by the adjacency matrix of a
    strongly regular graph. See [BH2012]_ and [BI1984]_ for details.

    If the set of parameters is not feasible, or if they correspond to a
    conference graph, the function returns ``None``. Its output is stable, assuming
    that the eigenvalues r,s used satisfy r>s; this holds for the current
    implementation of eigenvalues().

    INPUT:

    - ``v,k,l,mu`` (integers)

    EXAMPLES:

    Petersen's graph's C-algebra does not have a dual coming from an s.r.g.::

        sage: from sage.graphs.strongly_regular_db import eigenmatrix
        sage: P = eigenmatrix(10,3,0,1); P                                              # needs sage.modules
        [ 1  3  6]
        [ 1  1 -2]
        [ 1 -2  1]
        sage: 10*P^-1                                                                   # needs sage.modules
        [   1    5    4]
        [   1  5/3 -8/3]
        [   1 -5/3  2/3]

    The line graph of `K_{3,3}` is self-dual::

        sage: P = eigenmatrix(9,4,1,2); P                                               # needs sage.modules
        [ 1  4  4]
        [ 1  1 -2]
        [ 1 -2  1]
        sage: 9*P^-1                                                                    # needs sage.modules
        [ 1  4  4]
        [ 1  1 -2]
        [ 1 -2  1]

    A strongly regular graph with a non-isomorphic dual coming from another
    strongly regular graph::

        sage: # needs sage.modules
        sage: graphs.strongly_regular_graph(243,220,199,200, existence=True)            # needs sage.combinat
        True
        sage: graphs.strongly_regular_graph(243,110,37,60, existence=True)              # needs sage.combinat
        True
        sage: P = eigenmatrix(243,220,199,200); P
        [  1 220  22]
        [  1   4  -5]
        [  1  -5   4]
        sage: 243*P^-1
        [  1 110 132]
        [  1   2  -3]
        [  1 -25  24]
        sage: 243*P^-1==eigenmatrix(243,110,37,60)
        True

    TESTS::

        sage: eigenmatrix(5,5,5,-5)
    """
    from sage.rings.integer_ring import ZZ
    r, s = eigenvalues(v, k, l, mu)
    if r is not None:
        return Matrix(ZZ, [[1, k, v-k-1], [1, r, -r-1], [1, s, -s-1]])


cpdef latin_squares_graph_parameters(int v, int k, int l,int mu):
    r"""
    Check whether (v,k,l,mu)-strongly regular graph has parameters of an `L_g(n)` s.r.g.

    Also known as pseudo-OA(n,g) case, i.e. s.r.g. with parameters of an OA(n,g)-graph.
    Return g and n, if they exist. See Sect. 9.1 of [BH2012]_ for details.

    INPUT:

    - ``v,k,l,mu`` -- (integers) parameters of the graph

    OUTPUT:

    - ``(g, n)`` -- parameters of an `L_g(n)` graph, or `None`

    TESTS::

        sage: from sage.graphs.strongly_regular_db import latin_squares_graph_parameters
        sage: latin_squares_graph_parameters(9,4,1,2)
        (2, 3)
        sage: latin_squares_graph_parameters(5,4,1,2)
    """
    cdef int g, n
    r, s = eigenvalues(v, k, l, mu)
    if r is None:
        return
    if r < s:
        r, s = s, r
    g = -s
    n = r+g
    if v == n**2 and k == g*(n-1) and l == (g-1)*(g-2)+n-2 and mu == g*(g-1):
        return g, n
    return


def _H_3_cayley_graph(L):
    r"""
    return the `L`-Cayley graph of the group `H_3` from Prop. 12 in [JK2003]_.

    INPUT:

    - the list of words for the generating set in the format ["abc",...,"xyz"] for
      a,b,...,z being integers between 0 and 4.

    TESTS::

        sage: from sage.graphs.strongly_regular_db import _H_3_cayley_graph
        sage: _H_3_cayley_graph(["100","110","130","140","200","230","240","300"])      # needs sage.groups
        Graph on 100 vertices
    """
    from sage.groups.free_group import FreeGroup
    from sage.groups.finitely_presented import FinitelyPresentedGroup
    G = FreeGroup('x,y,z')
    x, y, z = G.gens()
    rels = (x**5, y**5, z**4, x*y*x**(-1)*y**(-1), z*x*z**(-1)*x**(-2), z*y*z**(-1)*y**(-2))
    G = FinitelyPresentedGroup(G, rels)
    x, y, z = G.gens()
    H = G.as_permutation_group()
    L = [[int(u) for u in x] for x in L]
    x, y, z = (H.gen(0), H.gen(1), H.gen(2))
    L = [H(x**xx*y**yy*z**zz) for xx, yy, zz in L]
    return Graph(H.cayley_graph(generators=L, simple=True))


def SRG_100_44_18_20():
    r"""
    Return a `(100, 44, 18, 20)`-strongly regular graph.

    This graph is built as a Cayley graph, using the construction for `\Delta_1`
    with group `H_3` presented in Table 8.1 of [JK2003]_

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_100_44_18_20
        sage: G = SRG_100_44_18_20()                    # long time
        sage: G.is_strongly_regular(parameters=True)    # long time
        (100, 44, 18, 20)
    """
    L = ['100', '110', '130', '140', '200', '230', '240', '300', '310', '320',
         '400', '410', '420', '440', '041', '111', '221', '231', '241', '321',
         '331', '401', '421', '441', '002', '042', '112', '122', '142', '212',
         '232', '242', '322', '342', '033', '113', '143', '223', '303', '333',
         '343', '413', '433', '443']
    return _H_3_cayley_graph(L)


def SRG_100_45_20_20():
    r"""
    Return a `(100, 45, 20, 20)`-strongly regular graph.

    This graph is built as a Cayley graph, using the construction for `\Gamma_3`
    with group `H_3` presented in Table 8.1 of [JK2003]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_100_45_20_20
        sage: G = SRG_100_45_20_20()                    # long time
        sage: G.is_strongly_regular(parameters=True)    # long time
        (100, 45, 20, 20)
    """
    L = ['120', '140', '200', '210', '201', '401', '411', '321', '002', '012',
         '022', '042', '303', '403', '013', '413', '240', '031', '102', '323',
         '300', '231', '132', '133', '310', '141', '142', '233', '340', '241',
         '202', '333', '410', '341', '222', '433', '430', '441', '242', '302',
         '312', '322', '332', '442', '143']
    return _H_3_cayley_graph(L)


def SRG_105_32_4_12():
    r"""
    Return a `(105, 32, 4, 12)`-strongly regular graph.

    The vertices are the flags of the projective plane of order 4. Two flags
    `(a,A)` and `(b,B)` are adjacent if the point `a` is on the line `B` or
    the point `b` is  on the line `A`, and `a \neq b`, `A \neq B`. See
    Theorem 2.7 in [GS1970]_, and [Coo2006]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_105_32_4_12
        sage: G = SRG_105_32_4_12(); G                                                  # needs sage.rings.finite_rings
        Aut L(3,4) on flags: Graph on 105 vertices
        sage: G.is_strongly_regular(parameters=True)                                    # needs sage.rings.finite_rings
        (105, 32, 4, 12)
    """
    from sage.combinat.designs.block_design import ProjectiveGeometryDesign
    P = ProjectiveGeometryDesign(2, 1, GF(4, 'a'))
    IG = P.incidence_graph().line_graph()
    a = IG.automorphism_group()
    h = a.stabilizer(a.domain()[0])
    o = next(x for x in h.orbits() if len(x) == 32)[0]
    e = a.orbit((a.domain()[0], o), action="OnSets")
    G = Graph()
    G.add_edges(e)
    G.name('Aut L(3,4) on flags')
    return G


def SRG_120_77_52_44():
    r"""
    Return a `(120,77,52,44)`-strongly regular graph.

    To build this graph, we first build a `2-(21,7,12)` design, by removing two
    points from the :func:`~sage.combinat.designs.block_design.WittDesign` on 23
    points. We then build the intersection graph of blocks with intersection
    size 3.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_120_77_52_44
        sage: G = SRG_120_77_52_44()                    # optional - gap_package_design
        sage: G.is_strongly_regular(parameters=True)    # optional - gap_package_design
        (120, 77, 52, 44)
    """
    from sage.combinat.designs.block_design import WittDesign
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    W = WittDesign(23)
    H = IncidenceStructure([x for x in W if 22 not in x and 21 not in x])
    g = H.intersection_graph(3)
    g.name('PG(2,2)s in PG(2,4)')
    return g


def SRG_144_39_6_12():
    r"""
    Return a `(144,39,6,12)`-strongly regular graph.

    This graph is obtained as an orbit of length 2808 on sets of cardinality 2
    (among 2 such orbits) of the group `PGL_3(3)` acting on the (right) cosets of
    a subgroup of order 39.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_144_39_6_12
        sage: G = SRG_144_39_6_12()                                                     # needs sage.libs.gap
        sage: G.is_strongly_regular(parameters=True)                                    # needs sage.libs.gap
        (144, 39, 6, 12)
    """
    from sage.libs.gap.libgap import libgap
    g = libgap.ProjectiveGeneralLinearGroup(3, 3)
    ns = g.Normalizer(g.SylowSubgroup(13))
    G = g.Action(g.RightCosets(ns), libgap.OnRight)
    H = G.Stabilizer(1)
    for o in H.Orbits():
        if len(o) != 39:
            continue
        h = Graph()
        h.add_edges(G.Orbit([1, o[0]], libgap.OnSets))
        if h.is_strongly_regular():
            h.relabel()
            h.name('PGL_3(3) on cosets of 13:3')
            return h


def SRG_176_49_12_14():
    r"""
    Return a `(176,49,12,14)`-strongly regular graph.

    This graph is built from the symmetric Higman-Sims design. In
    [Bro1982]_, it is explained that there exists an involution
    `\sigma` exchanging the points and blocks of the Higman-Sims design, such
    that each point is mapped on a block that contains it (i.e. `\sigma` is a
    'polarity with all universal points'). The graph is then built by making two
    vertices `u,v` adjacent whenever `v\in \sigma(u)`.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_176_49_12_14
        sage: G = SRG_176_49_12_14()                    # long time, optional - gap_package_design
        sage: G.is_strongly_regular(parameters=True)    # long time, optional - gap_package_design
        (176, 49, 12, 14)
    """
    from sage.combinat.designs.database import HigmanSimsDesign
    d = HigmanSimsDesign()
    g = d.incidence_graph(labels=True)
    ag = g.automorphism_group().conjugacy_classes_representatives()

    # Looking for an involution that maps a point of the design to one of the
    # blocks that contains it. It is called a polarity with only absolute
    # points.
    for aut in ag:
        try:
            0 in aut(0)
        except TypeError:
            continue
        if (aut.order() == 2 and
                all(i in aut(i) for i in d.ground_set())):
            g = Graph()
            g.add_edges(((u, v) for u in d.ground_set() for v in aut(u)), loops=False)
            g.name('Higman symmetric 2-design')
            return g


def SRG_176_105_68_54():
    r"""
    Return a `(176, 105, 68, 54)`-strongly regular graph.

    To build this graph, we first build a `2-(22,7,16)` design, by removing one
    point from the :func:`~sage.combinat.designs.block_design.WittDesign` on 23
    points. We then build the intersection graph of blocks with intersection
    size 3. Known as S.7 in [Hub1975]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_176_105_68_54
        sage: G = SRG_176_105_68_54()                   # optional - gap_package_design
        sage: G.is_strongly_regular(parameters=True)    # optional - gap_package_design
        (176, 105, 68, 54)
    """
    from sage.combinat.designs.block_design import WittDesign
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    W = WittDesign(23)
    H = IncidenceStructure([x for x in W if 22 not in x])
    g = H.intersection_graph(3)
    g.name('Witt 3-(22,7,4)')
    return g


def SRG_210_99_48_45():
    r"""
    Return a strongly regular graph with parameters `(210, 99, 48, 45)`

    This graph is from Example 4.2 in [KPRWZ2010]_. One considers the action of
    the symmetric group `S_7` on the 210 digraphs isomorphic to the
    disjoint union of `K_1` and the circulant 6-vertex digraph
    ``digraphs.Circulant(6,[1,4])``. It has 16 orbitals; the package [FK1991]_
    found a megring of them, explicitly described in [KPRWZ2010]_, resulting in
    this graph.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_210_99_48_45
        sage: g = SRG_210_99_48_45()                                                    # needs sage.libs.gap
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.libs.gap
        (210, 99, 48, 45)
    """
    from sage.libs.gap.libgap import libgap
    from sage.combinat.permutation import Permutation

    def ekg(g0):
        # return arcs of the Cayley digraph of <g> on {g,g^4}
        g = Permutation(g0)
        return libgap.Set([(x, g(x)) for x in range(1, 8)] +
                          [(x, g(g(g(g(x))))) for x in range(1, 8)])

    kd = list(map(ekg,
                  [(7, 1, 2, 3, 4, 5), (7, 1, 3, 4, 5, 6),
                   (7, 3, 4, 5, 6, 2), (7, 1, 4, 3, 5, 6),
                   (7, 3, 1, 4, 5, 6), (7, 2, 4, 3, 5, 6),
                   (7, 3, 2, 4, 5, 1), (7, 2, 4, 3, 5, 1)]))
    s = libgap.SymmetricGroup(7)
    O = s.Orbit(kd[0], libgap.OnSetsTuples)
    sa = s.Action(O, libgap.OnSetsTuples)
    G = Graph()
    for g in kd[1:]:
        G.add_edges(libgap.Orbit(sa, [libgap.Position(O, kd[0]),
                                      libgap.Position(O, g)], libgap.OnSets))
    G.name('merging of S_7 on Circulant(6,[1,4])s')
    return G


def SRG_243_110_37_60():
    r"""
    Return a `(243, 110, 37, 60)`-strongly regular graph.

    Consider the orthogonal complement of the
    :func:`~sage.coding.code_constructions.TernaryGolayCode`, which has 243
    words. On them we define a graph, in which two words are adjacent
    whenever their Hamming distance is 9. This construction appears in
    [GS1975]_.

    .. NOTE::

        A strongly regular graph with the same parameters is also obtained from
        the database of 2-weight codes.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_243_110_37_60
        sage: G = SRG_243_110_37_60()                                                   # needs sage.modules sage.rings.finite_rings
        sage: G.is_strongly_regular(parameters=True)                                    # needs sage.modules sage.rings.finite_rings
        (243, 110, 37, 60)
    """
    from sage.coding.golay_code import GolayCode
    M = GolayCode(GF(3), False).generator_matrix()
    V = list(M.right_kernel())
    g = Graph([list(xrange(len(V))), lambda x, y: (V[x] - V[y]).hamming_weight() == 9])
    g.name('Ternary Golay code')
    return g


def SRG_253_140_87_65():
    r"""
    Return a `(253, 140, 87, 65)`-strongly regular graph.

    To build this graph, we first build the
    :func:`~sage.combinat.designs.block_design.WittDesign` on 23 points which is
    a `2-(23,7,21)` design. We then build the intersection graph of blocks with
    intersection size 3. Known as S.6 in [Hub1975]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_253_140_87_65
        sage: G = SRG_253_140_87_65()                   # optional - gap_package_design
        sage: G.is_strongly_regular(parameters=True)    # optional - gap_package_design
        (253, 140, 87, 65)
    """
    from sage.combinat.designs.block_design import WittDesign
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    W = WittDesign(23)
    g = W.intersection_graph(3)
    g.name('Witt 4-(23,7,1)')
    return g


def SRG_196_91_42_42():
    r"""
    Return a `(196,91,42,42)`-strongly regular graph.

    This strongly regular graph is built following the construction provided in
    Corollary 8.2.27 of [IS2006]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_196_91_42_42
        sage: G = SRG_196_91_42_42()
        sage: G.is_strongly_regular(parameters=True)
        (196, 91, 42, 42)
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    from sage.graphs.generators.intersection import IntersectionGraph
    k = 7
    R = IntegerModRing(91)
    A = list(map(R, [0, 10, 27, 28, 31, 43, 50]))
    B = list(map(R, [0, 11, 20, 25, 49, 55, 57]))
    H = list(map(R, [13 * i for i in range(k)]))
    U = list(map(frozenset, [[x + z for x in A] for z in R]))
    V = list(map(frozenset, [[x + z for x in B] for z in R]))
    W = list(map(frozenset, [[x + z for x in H] for z in R]))
    G = IntersectionGraph(U + V + W)

    G.seidel_switching(U)

    G.add_edges((-1, x) for x in U)
    G.relabel(perm={u: i for i, u in enumerate(G)})
    G.name('RSHCD+')
    return G


def SRG_220_84_38_28():
    r"""
    Return a `(220, 84, 38, 28)`-strongly regular graph.

    This graph is obtained from the
    :meth:`~IncidenceStructure.intersection_graph` of a
    :func:`~sage.combinat.designs.database.BIBD_45_9_8`. This construction
    appears in VII.11.2 from [DesignHandbook]_

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_220_84_38_28
        sage: g=SRG_220_84_38_28()
        sage: g.is_strongly_regular(parameters=True)
        (220, 84, 38, 28)
    """
    from sage.combinat.designs.database import BIBD_45_9_8
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    G = IncidenceStructure(BIBD_45_9_8()).intersection_graph(3)
    G.relabel()
    G.name('Tonchev: quasisymmetric 2-(45,9,8)')
    return G


def SRG_276_140_58_84():
    r"""
    Return a `(276, 140, 58, 84)`-strongly regular graph.

    The graph is built from
    :meth:`~sage.graphs.graph_generators.GraphGenerators.McLaughlinGraph`, with
    an added isolated vertex. We then perform a
    :meth:`~Graph.seidel_switching` on a set of 28 disjoint 5-cliques, which
    exist by cf. [HT1996]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_276_140_58_84
        sage: g = SRG_276_140_58_84()                   # long time, optional - gap_package_design
        sage: g.is_strongly_regular(parameters=True)    # long time, optional - gap_package_design
        (276, 140, 58, 84)
    """
    from sage.graphs.generators.smallgraphs import McLaughlinGraph
    g = McLaughlinGraph()
    C = [[ 0,  72,  87, 131, 136], [ 1,  35,  61, 102, 168], [ 2,  32,  97, 125, 197], [ 3,  22,  96, 103, 202],
         [ 4,  46,  74, 158, 229], [ 5,  83,  93, 242, 261], [ 6,  26,  81, 147, 176], [ 7,  42,  63, 119, 263],
         [ 8,  49,  64, 165, 227], [ 9,  70,  85, 208, 273], [10,  73,  92, 230, 268], [11,  54,  95, 184, 269],
         [12,  55,  62, 185, 205], [13,  51,  65, 162, 254], [14,  78,  88, 231, 274], [15,  40,  59, 117, 252],
         [16,  24,  71, 137, 171], [17,  39,  43, 132, 163], [18,  57,  79, 175, 271], [19,  68,  80, 217, 244],
         [20,  75,  98, 239, 267], [21,  33,  56, 113, 240], [23, 127, 152, 164, 172], [25, 101, 128, 183, 264],
         [27, 129, 154, 160, 201], [28, 126, 144, 161, 228], [29, 100, 133, 204, 266], [30, 108, 146, 200, 219]]
    g.add_vertex(-1)
    g.seidel_switching(sum(C, []))
    g.relabel()
    g.name('Haemers-Tonchev')
    return g


def SRG_280_135_70_60():
    r"""
    Return a strongly regular graph with parameters `(280, 135, 70, 60)`.

    This graph is built from the action of `J_2` on the cosets of a `3.PGL(2,9)`-subgroup.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_280_135_70_60
        sage: g=SRG_280_135_70_60()                     # long time, optional - internet
        sage: g.is_strongly_regular(parameters=True)    # long time, optional - internet
        (280, 135, 70, 60)
    """
    from sage.libs.gap.libgap import libgap
    from sage.graphs.graph import Graph

    libgap.load_package("AtlasRep")
    # A representation of J2 acting on a 3.PGL(2,9) it contains.
    J2 = libgap.AtlasGroup("J2", libgap.NrMovedPoints, 280)
    edges = J2.Orbit([1, 2], libgap.OnSets)
    g = Graph()
    g.add_edges(edges)
    g.relabel()
    g.name('J_2 on cosets of 3.PGL(2,9)')
    return g


def SRG_280_117_44_52():
    r"""
    Return a strongly regular graph with parameters `(280, 117, 44, 52)`.

    This graph is built according to a very pretty construction of Mathon and
    Rosa [MR1985]_:

        The vertices of the graph `G` are all partitions of a set of 9 elements
        into `\{\{a,b,c\},\{d,e,f\},\{g,h,i\}\}`. The cross-intersection of two
        such partitions `P=\{P_1,P_2,P_3\}` and `P'=\{P'_1,P'_2,P'_3\}` being
        defined as `\{P_i \cap P'_j: 1\leq i,j\leq 3\}`, two vertices of `G` are
        set to be adjacent if the cross-intersection of their respective
        partitions does not contain exactly 7 nonempty sets.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_280_117_44_52
        sage: g=SRG_280_117_44_52()
        sage: g.is_strongly_regular(parameters=True)
        (280, 117, 44, 52)
    """
    from sage.graphs.hypergraph_generators import hypergraphs

    # V is the set of partitions {{a,b,c},{d,e,f},{g,h,i}} of {0,...,8}
    H = hypergraphs.CompleteUniform(9, 3)
    g = H.intersection_graph()
    V = g.complement().cliques_maximal()
    V = [frozenset(u) for u in V]

    # G is the graph defined on V in which two vertices are adjacent when they
    # corresponding partitions cross-intersect on 7 nonempty sets
    G = Graph([V, lambda x, y:
               sum(any(xxx in yy for xxx in xx) for xx in x for yy in y) != 7],
              loops=False)
    G.name('Mathon-Rosa')
    return G


def strongly_regular_from_two_weight_code(L):
    r"""
    Return a strongly regular graph from a two-weight code.

    A code is said to be a *two-weight* code the weight of its nonzero codewords
    (i.e. their number of nonzero coordinates) can only be one of two integer
    values `w_1,w_2`. It is said to be *projective* if the minimum weight of the
    dual code is `\geq 3`. A strongly regular graph can be built from a
    two-weight projective code with weights `w_1,w_2` (assuming `w_1<w_2`) by
    adding an edge between any two codewords whose difference has weight
    `w_1`. For more information, see [LS1981]_ or [Del1972]_.

    INPUT:

    - ``L`` -- a two-weight linear code, or its generating matrix.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import strongly_regular_from_two_weight_code
        sage: x = ("100022021001111",
        ....:      "010011211122000",
        ....:      "001021112100011",
        ....:      "000110120222220")
        sage: M = Matrix(GF(3),[list(l) for l in x])                                    # needs sage.modules sage.rings.finite_rings
        sage: G = strongly_regular_from_two_weight_code(LinearCode(M))                  # needs sage.modules sage.rings.finite_rings
        sage: G.is_strongly_regular(parameters=True)                                    # needs sage.modules sage.rings.finite_rings
        (81, 50, 31, 30)
    """
    from sage.structure.element import is_Matrix
    if is_Matrix(L):
        L = LinearCode(L)
    V = [tuple(l) for l in L]
    w1, w2 = sorted(set(sum(map(bool, x)) for x in V).difference([0]))
    G = Graph([V, lambda u, v: sum(uu!=vv for uu, vv in zip(u, v)) == w1])
    G.relabel()
    G.name('two-weight code: '+str(L))
    return G


def SRG_416_100_36_20():
    r"""
    Return a `(416,100,36,20)`-strongly regular graph.

    This graph is obtained as an orbit on sets of cardinality 2
    (among 2 that exists) of the group `G_2(4)`.
    This graph is isomorphic to the subgraph of the from :meth:`Suzuki Graph
    <sage.graphs.graph_generators.GraphGenerators.SuzukiGraph>` induced on
    the neighbors of a vertex. Known as S.14 in [Hub1975]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_416_100_36_20
        sage: g = SRG_416_100_36_20()                   # long time, optional - internet
        sage: g.is_strongly_regular(parameters=True)    # long time, optional - internet
        (416, 100, 36, 20)
    """
    from sage.libs.gap.libgap import libgap
    libgap.load_package("AtlasRep")
    g=libgap.AtlasGroup("G2(4)", libgap.NrMovedPoints, 416)
    h = Graph()
    h.add_edges(g.Orbit([1, 5],libgap.OnSets))
    h.relabel()
    h.name('G_2(4) on cosets of HS')
    return h


def SRG_560_208_72_80():
    r"""
    Return a `(560,208,72,80)`-strongly regular graph

    This graph is obtained as the union of 4 orbits of sets of cardinality 2
    (among the 13 that exist) of the group `Sz(8)`.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_560_208_72_80
        sage: g = SRG_560_208_72_80()                   # not tested (~2s)
        sage: g.is_strongly_regular(parameters=True)    # not tested (~2s)
        (560, 208, 72, 80)
    """
    from sage.libs.gap.libgap import libgap
    libgap.load_package("AtlasRep")
    g=libgap.AtlasGroup("Sz8", libgap.NrMovedPoints, 560)

    h = Graph()
    h.add_edges(g.Orbit([1, 2],libgap.OnSets))
    h.add_edges(g.Orbit([1, 4],libgap.OnSets))
    h.add_edges(g.Orbit([1, 8],libgap.OnSets))
    h.add_edges(g.Orbit([1, 27],libgap.OnSets))
    h.relabel()
    h.name('Sz(8)-graph')
    return h


def strongly_regular_from_two_intersection_set(M):
    r"""
    Return a strongly regular graph from a 2-intersection set.

    A set of points in the projective geometry `PG(k,q)` is said to be a
    2-intersection set if it intersects every hyperplane in either `h_1` or
    `h_2` points, where `h_1,h_2\in \\NN`.

    From a 2-intersection set `S` can be defined a strongly-regular graph in the
    following way:

    - Place the points of `S` on a hyperplane `H` in `PG(k+1,q)`

    - Define the graph `G` on all points of `PG(k+1,q)\backslash H`

    - Make two points of `V(G)=PG(k+1,q)\backslash H` adjacent if the line going
      through them intersects `S`

    For more information, see e.g. [CD2013]_ where this explanation has been
    taken from.

    INPUT:

    - `M` -- a `|S| \times k` matrix with entries in `F_q` representing the points of
      the 2-intersection set. We assume that the first non-zero entry of each row is
      equal to `1`, that is, they give points in homogeneous coordinates.

    The implementation does not check that `S` is actually a 2-intersection set.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import strongly_regular_from_two_intersection_set
        sage: S = Matrix([(0,0,1),(0,1,0)] + [(1,x^2,x) for x in GF(4,'b')])            # needs sage.modules sage.rings.finite_rings
        sage: g = strongly_regular_from_two_intersection_set(S); g                      # needs sage.modules sage.rings.finite_rings
        two-intersection set in PG(3,4): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.modules sage.rings.finite_rings
        (64, 18, 2, 6)
    """
    from itertools import product
    from sage.rings.rational_field import QQ
    K = M.base_ring()
    k = M.ncols()
    g = Graph()

    M = [list(p) for p in M]

    # For every point in F_q^{k+1} not on the hyperplane of M
    for u in [tuple(x) for x in product(K,repeat=k)]:
        # For every v point of M
        for v in M:
            # u is adjacent with all vertices on a uv line.
            g.add_edges([[u, tuple([u[i] + qq*v[i] for i in range(k)])]
                         for qq in K if not qq==K.zero()])
    g.relabel()
    e = QQ((1,k))
    qq = g.num_verts()**e
    g.name('two-intersection set in PG('+str(k)+','+str(qq)+')')
    return g


def SRG_120_63_30_36():
    r"""
    Return a `(120,63,30,36)`-strongly regular graph

    It is the distance-2 graph of :meth:`JohnsonGraph(10,3)
    <sage.graphs.graph_generators.GraphGenerators.JohnsonGraph>`.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_120_63_30_36
        sage: G =  SRG_120_63_30_36()
        sage: G.is_strongly_regular(parameters=True)
        (120, 63, 30, 36)
    """
    from sage.graphs.generators.families import JohnsonGraph
    return JohnsonGraph(10, 3).distance_graph([2])


def SRG_126_25_8_4():
    r"""
    Return a `(126,25,8,4)`-strongly regular graph

    It is the distance-(1 or 4) graph of :meth:`JohnsonGraph(9,4)
    <sage.graphs.graph_generators.GraphGenerators.JohnsonGraph>`.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_126_25_8_4
        sage: G =  SRG_126_25_8_4()
        sage: G.is_strongly_regular(parameters=True)
        (126, 25, 8, 4)
    """
    from sage.graphs.generators.families import JohnsonGraph
    return JohnsonGraph(9, 4).distance_graph([1, 4])


def SRG_175_72_20_36():
    r"""
    Return a `(175,72,20,36)`-strongly regular graph

    This graph is obtained from the line graph of
    :meth:`~sage.graphs.graph_generators.GraphGenerators.HoffmanSingletonGraph`. Setting
    two vertices to be adjacent if their distance in the line graph is exactly
    2 yields the graph. For more information, see 10.B.(iv) in [BL1984]_ and
    https://www.win.tue.nl/~aeb/graphs/McL.html.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_175_72_20_36
        sage: G = SRG_175_72_20_36()
        sage: G.is_strongly_regular(parameters=True)
        (175, 72, 20, 36)
    """
    from sage.graphs.generators.smallgraphs import HoffmanSingletonGraph
    return HoffmanSingletonGraph().line_graph().distance_graph([2])


def SRG_176_90_38_54():
    r"""
    Return a `(176,90,38,54)`-strongly regular graph

    This graph is obtained from
    :func:`~sage.graphs.strongly_regular_db.SRG_175_72_20_36`
    by attaching a isolated vertex and doing Seidel switching
    with respect to disjoint union of 18 maximum cliques, following
    a construction by W.Haemers given in Sect.10.B.(vi) of [BL1984]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_176_90_38_54
        sage: G = SRG_176_90_38_54(); G
        a Seidel switching of Distance graph for distance 2 in : Graph on 176 vertices
        sage: G.is_strongly_regular(parameters=True)
        (176, 90, 38, 54)
    """
    g = SRG_175_72_20_36()
    g.relabel(range(175))
    # c=filter(lambda x: len(x)==5, g.cliques_maximal())
    # r=flatten(Hypergraph(c).packing()[:18]) # takes 3s, so we put the answer here
    r = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
         20, 21, 22, 23, 24, 25, 28, 29, 32, 38, 39, 41, 42, 43, 47, 49, 50, 51,
         52, 53, 55, 57, 61, 63, 65, 67, 69, 72, 75, 77, 79, 81, 84, 87, 88, 89,
         92, 95, 96, 97, 99, 101, 102, 104, 105, 107, 112, 114, 117, 118, 123,
         125, 129, 132, 139, 140, 141, 144, 146, 147, 153, 154, 162, 165, 166,
         167, 170, 172, 173, 174]
    g.add_vertex()
    g.seidel_switching(r)
    g.name('a Seidel switching of ' + g.name())
    return g


def SRG_630_85_20_10():
    r"""
    Return a `(630,85,20,10)`-strongly regular graph

    This graph is the line graph of `pg(5,18,2)`; its point graph is
    :func:`~sage.graphs.strongly_regular_db.SRG_175_72_20_36`.
    One selects a subset of 630 maximum cliques in the latter following
    a construction by W.Haemers given in Sect.10.B.(v) of [BL1984]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_630_85_20_10
        sage: G = SRG_630_85_20_10()                    # long time
        sage: G.is_strongly_regular(parameters=True)    # long time
        (630, 85, 20, 10)
    """
    from sage.graphs.generators.intersection import IntersectionGraph
    from sage.graphs.generators.smallgraphs import HoffmanSingletonGraph
    hs = HoffmanSingletonGraph()
    P = list(range(5)) + list(range(30, 35))  # a Petersen in hs
    mc = [0, 1, 5, 6, 12, 13, 16, 17, 22, 23, 29, 33, 39, 42, 47]
    assert(hs.subgraph(mc).is_regular(k=0))  # a maximum coclique
    assert(hs.subgraph(P).is_regular(k=3))
    h = hs.automorphism_group().stabilizer(mc, action="OnSets")
    l = h.orbit(tuple((x[0], x[1]) for x in hs.subgraph(P).matching()),
                "OnSetsSets")
    return IntersectionGraph(l)


def SRG_126_50_13_24():
    r"""
    Return a `(126,50,13,24)`-strongly regular graph

    This graph is a subgraph of
    :meth:`~sage.graphs.strongly_regular_db.SRG_175_72_20_36`.
    This construction, due to Goethals, is given in §10B.(vii) of [BL1984]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_126_50_13_24
        sage: G = SRG_126_50_13_24(); G
        Goethals graph: Graph on 126 vertices
        sage: G.is_strongly_regular(parameters=True)
        (126, 50, 13, 24)
    """
    from sage.graphs.strongly_regular_db import SRG_175_72_20_36
    from sage.graphs.generators.smallgraphs import HoffmanSingletonGraph
    hs = HoffmanSingletonGraph()
    s = set(hs.vertices(sort=False)).difference(hs.neighbors(0) + [0])
    g = SRG_175_72_20_36().subgraph(hs.edge_boundary(s, s))
    g.name('Goethals graph')
    return g


def SRG_1288_792_476_504():
    r"""
    Return a `(1288, 792, 476, 504)`-strongly regular graph.

    This graph is built on the words of weight 12 in the
    :func:`~sage.coding.code_constructions.BinaryGolayCode`. Two of them are
    then made adjacent if their symmetric difference has weight 12 (cf
    [BE1992]_).

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_1288_792_476_504
        sage: G = SRG_1288_792_476_504()                # long time                     # needs sage.rings.finite_rings
        sage: G.is_strongly_regular(parameters=True)    # long time                     # needs sage.rings.finite_rings
        (1288, 792, 476, 504)
    """
    from sage.coding.golay_code import GolayCode
    C = GolayCode(GF(2), False)
    C = [[i for i,v in enumerate(c) if v]
         for c in C]
    C = [s for s in C if len(s) == 12]
    G = Graph([[frozenset(c) for c in C],
               lambda x, y: len(x.symmetric_difference(y)) == 12])
    G.relabel()
    G.name('binary Golay code')
    return G


cdef bint seems_feasible(int v, int k, int l, int mu):
    r"""
    Check if the set of parameters seems feasible.

    INPUT:

    - ``v,k,l,mu`` (integers)

    TESTS:

    :trac:`32306` is fixed::

        sage: from sage.graphs.strongly_regular_db import strongly_regular_graph
        sage: strongly_regular_graph(16384, 8256, 4160, 4160, existence=True)           # needs sage.combinat sage.modules
        True
    """
    cdef uint_fast32_t tmp[2]

    if (v < 0 or k <= 0 or l < 0 or mu < 0 or
            k >= v - 1 or l >= k or mu > k or
            v - 2*k + mu - 2 < 0 or  # lambda of complement graph >=0
            v - 2*k + l < 0 or       # μ of complement graph >= 0
            mu*(v - k - 1) != k*(k - l - 1)):
        return False

    if mu == k:  # complete multipartite graph
        r = v//(v-k)  # number of parts (of size v-k each)
        return (l == (v-k)*(r-2) and v == r*(v-k))

    if mu == 0:  # the complement of a complete multipartite graph
        r = v//(k+1)  # number of parts (of size k+1 each)
        return (l == k-1 and v == r*(k+1))

    # Conference graphs. Only possible if 'v' is a sum of two squares (3.A of
    # [BL1984]_)
    if (v-1)*(mu-l)-2*k == 0:
        return two_squares_c(v, tmp)

    rr, ss = eigenvalues(v, k, l, mu)
    if rr is None:
        return False
    r, s = rr, ss

    # p.87 of [BL1984]_
    # "Integrality condition"
    if ((s+1)*(k-s)*k) % (mu*(s-r)) or ((r+1)*(k-r)*k) % (mu*(s-r)):
        return False

    # Theorem 21.3 of [WilsonACourse] or
    # 3.B of [BL1984]_
    # (Krein conditions)
    if (r+1)*(k+r+2*r*s) > (k+r)*(s+1)**2 or (s+1)*(k+s+2*r*s) > (k+s)*(r+1)**2:
        return False

    # multiplicity of eigenvalues 'r,s' (f=lambda_r, g=lambda_s)
    #
    # They are integers (checked by the 'integrality condition').
    f = -k*(s+1)*(k-s)//(mu*(r-s))
    g = k*(r+1)*(k-r)//(mu*(r-s))
    if 1 + f + g != v:  # the only other eigenvalue, k, has multiplicity 1
        return False

    # 3.C of [BL1984]_
    # (Absolute bound)
    if 2*v > f*(f+3) or 2*v > g*(g+3):
        return False

    # 3.D of [BL1984]_
    # (Claw bound)
    if (mu != s**2 and
            mu != s*(s+1) and
            2*(r+1) > s*(s+1)*(mu+1)):
        return False

    # 3.E of [BL1984]_
    # (the Case μ=1)
    if mu == 1:
        if k % (l+1) or (v*k) % ((l+1)*(l+2)):
            return False

    # 3.F of [BL1984]_
    # (the Case μ=2)
    if mu == 2 and 2*k < l*(l + 3) and k % (l + 1):
        return False

    return True


def strongly_regular_graph(int v, int k, int l, int mu=-1, bint existence=False, bint check=True):
    r"""
    Return a `(v,k,\lambda,\mu)`-strongly regular graph.

    This function relies partly on Andries Brouwer's `database of strongly
    regular graphs <https://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__. See
    the documentation of :mod:`sage.graphs.strongly_regular_db` for more
    information.

    INPUT:

    - ``v,k,l,mu`` (integers) -- note that ``mu``, if unspecified, is
      automatically determined from ``v,k,l``.

    - ``existence`` (boolean;``False``) -- instead of building the graph,
      return:

        - ``True`` -- meaning that a `(v,k,\lambda,\mu)`-strongly regular graph
          exists.

        - ``Unknown`` -- meaning that Sage does not know if such a strongly
          regular graph exists (see :mod:`sage.misc.unknown`).

        - ``False`` -- meaning that no such strongly regular graph exists.

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to
      ``True`` by default.

    EXAMPLES:

    Petersen's graph from its set of parameters::

        sage: graphs.strongly_regular_graph(10,3,0,1,existence=True)
        True
        sage: graphs.strongly_regular_graph(10,3,0,1)
        complement(Johnson graph with parameters 5,2): Graph on 10 vertices

    Now without specifying `\mu`::

        sage: graphs.strongly_regular_graph(10,3,0)
        complement(Johnson graph with parameters 5,2): Graph on 10 vertices

    An obviously infeasible set of parameters::

        sage: graphs.strongly_regular_graph(5,5,5,5,existence=True)
        False
        sage: graphs.strongly_regular_graph(5,5,5,5)
        Traceback (most recent call last):
        ...
        ValueError: There exists no (5, 5, 5, 5)-strongly regular graph

    An set of parameters proved in a paper to be infeasible::

        sage: graphs.strongly_regular_graph(324,57,0,12,existence=True)                 # needs sage.combinat sage.modules
        False
        sage: graphs.strongly_regular_graph(324,57,0,12)                                # needs sage.combinat sage.modules
        Traceback (most recent call last):
        ...
        EmptySetError: Andries Brouwer's database reports that no (324, 57, 0,
        12)-strongly regular graph exists. Comments: <a
        href="srgtabrefs.html#GavrilyukMakhnev05">Gavrilyuk & Makhnev</a> ...

    A set of parameters unknown to be realizable in Andries Brouwer's database::

        sage: graphs.strongly_regular_graph(324,95,22,30,existence=True)                # needs sage.combinat
        Unknown
        sage: graphs.strongly_regular_graph(324,95,22,30)                               # needs sage.combinat
        Traceback (most recent call last):
        ...
        RuntimeError: Andries Brouwer's database reports that no
        (324, 95, 22, 30)-strongly regular graph is known to exist.
        Comments:

    A large unknown set of parameters (not in Andries Brouwer's database)::

        sage: graphs.strongly_regular_graph(1394,175,0,25,existence=True)               # needs sage.combinat
        Unknown
        sage: graphs.strongly_regular_graph(1394,175,0,25)                              # needs sage.combinat
        Traceback (most recent call last):
        ...
        RuntimeError: Sage cannot figure out if a (1394, 175, 0, 25)-strongly
        regular graph exists.

    Test the Claw bound (see 3.D of [BL1984]_)::

        sage: graphs.strongly_regular_graph(2058,242,91,20,existence=True)
        False

    TESTS:

    Check that :trac:`26513` is fixed::

        sage: graphs.strongly_regular_graph(539, 288, 162, 144)                         # needs sage.combinat
        descendant of (540, 264, 138, 120)-strongly regular graph at ... 539 vertices
        sage: graphs.strongly_regular_graph(539, 250, 105, 125)                         # needs sage.combinat
        descendant of (540, 275, 130, 150)-strongly regular graph at ... 539 vertices
        sage: graphs.strongly_regular_graph(209, 100, 45, 50)                           # needs sage.libs.pari
        descendant of complement(merging of S_7 on Circulant(6,[1,4])s) at ... 209 vertices


    Check that all of our constructions are correct - you will need gap_packages spkg installed::

        sage: from sage.graphs.strongly_regular_db import apparently_feasible_parameters
        sage: for p in sorted(apparently_feasible_parameters(1300)):   # not tested, optional gap_package_design
        ....:     if graphs.strongly_regular_graph(*p,existence=True) is True:
        ....:         try:
        ....:             _ = graphs.strongly_regular_graph(*p)
        ....:             print(p, "built successfully")
        ....:         except RuntimeError as e:
        ....:             if 'Brouwer' not in str(e):
        ....:                 raise

    `\mu=0` behaves correctly (:trac:`19712`)::

        sage: graphs.strongly_regular_graph(10,2,1)
        Traceback (most recent call last):
        ...
        ValueError: There exists no (10, 2, 1, 0)-strongly regular graph
        sage: graphs.strongly_regular_graph(12,3,2)
        complement(Multipartite Graph with set sizes [4, 4, 4]): Graph on 12 vertices
        sage: graphs.strongly_regular_graph(6,3,0)
        Multipartite Graph with set sizes [3, 3]: Graph on 6 vertices
    """
    if mu == -1:
        mu = k*(k - l - 1)//(v - k - 1)
    g = strongly_regular_graph_lazy(v, k, l, mu=mu, existence=existence)
    if existence is True:
        return g
    G = g[0](*g[1:])
    if check and (v, k, l, mu) != G.is_strongly_regular(parameters=True):
        params = (v, k, l, mu)
        raise RuntimeError(f"Sage built an incorrect {params}-SRG.")
    return G


def strongly_regular_graph_lazy(int v, int k, int l, int mu=-1, bint existence=False):
    r"""
    return a promise to build an `(v,k,l,mu)`-srg

    Return a promise to build an `(v,k,l,mu)`-srg as a tuple `t`, with `t[0]` a
    function to evaluate on `*t[1:]`.

    Input as in :func:`~sage.graphs.strongly_regular_graphs_db.strongly_regular_graph`,
    although without `check`.

    TESTS::

        sage: from sage.graphs.strongly_regular_db import strongly_regular_graph_lazy
        sage: g,p=strongly_regular_graph_lazy(10,6,3); g,p
        (<cyfunction is_johnson.<locals>.<lambda> at ...>, 5)
        sage: g(p)
        Johnson graph with parameters 5,2: Graph on 10 vertices
        sage: g,p=strongly_regular_graph_lazy(10,3,0,1); g,p
        (<cyfunction strongly_regular_graph_lazy.<locals>.<lambda> at...>,
         (5,))
        sage: g(p)
        complement(Johnson graph with parameters 5,2): Graph on 10 vertices
        sage: g,p=strongly_regular_graph_lazy(12,3,2); g,p
        (<cyfunction strongly_regular_graph_lazy.<locals>.<lambda> at...>,
         (3, 4))
        sage: g(p)
        complement(Multipartite Graph with set sizes [4, 4, 4]): Graph on 12 vertices
        sage: g = strongly_regular_graph_lazy(539,250,105); g                           # needs sage.combinat sage.modules
        (<cyfunction is_twograph_descendant_of_srg.<locals>.la at...>,
         5,
         11)
        sage: g[0](*g[1:])                                                              # needs sage.combinat sage.modules
        descendant of (540, 275, 130, 150)-strongly regular graph at 0: Graph on 539 vertices
    """
    load_brouwer_database()
    if mu == -1:
        mu = k*(k - l - 1)//(v - k - 1)

    params = (v, k, l, mu)
    params_complement = (v, v - k - 1, v - 2*k + mu - 2, v - 2*k + l)

    if not seems_feasible(v, k, l, mu):
        if existence:
            return False
        raise ValueError(f"There exists no {params}-strongly regular graph")

    if _small_srg_database is None:
        _build_small_srg_database()

    if params in _small_srg_database:
        val = _small_srg_database[params]
        return True if existence else (val[0], *val[1:])
    if params_complement in _small_srg_database:
        val = _small_srg_database[params_complement]
        return True if existence else (lambda *t: val[0](*t).complement(), *val[1:])

    test_functions = [is_complete_multipartite,  # must be 1st, to prevent 0-divisions
                      is_paley, is_johnson,
                      is_orthogonal_array_block_graph,
                      is_steiner, is_affine_polar,
                      is_goethals_seidel,
                      is_orthogonal_polar,
                      is_NOodd, is_NOperp_F5, is_NO_F2, is_NO_F3, is_NU,
                      is_unitary_polar, is_unitary_dual_polar, is_GQqmqp,
                      is_RSHCD,
                      is_twograph_descendant_of_srg,
                      is_taylor_twograph_srg,
                      is_switch_OA_srg,
                      is_polhill,
                      is_haemers,
                      is_cossidente_penttila,
                      is_mathon_PC_srg,
                      is_muzychuk_S6,
                      is_nowhere0_twoweight,
                      is_switch_skewhad]

    # Going through all test functions, for the set of parameters and its
    # complement.
    for f in test_functions:
        if f(*params):
            if existence:
                return True
            ans = f(*params)
            return (ans[0], *ans[1:])
        if f(*params_complement):
            if existence:
                return True
            ans = f(*params_complement)
            return (lambda t: ans[0](*t).complement(), ans[1:])

    # From now on, we have no idea how to build the graph.
    #
    # We try to return the most appropriate error message.

    global _brouwer_database
    brouwer_data = _brouwer_database.get(params, None)

    if brouwer_data is not None:
        comments = brouwer_data['comments']
        if brouwer_data['status'] == 'impossible':
            if existence:
                return False
            raise EmptySetError(
                f"Andries Brouwer's database reports that no "
                f"{params}-strongly regular graph exists. Comments: {comments}")

        if brouwer_data['status'] == 'open':
            if existence:
                return Unknown
            raise RuntimeError(
                f"Andries Brouwer's database reports that no "
                f"{params}-strongly regular graph is known to exist.\n"
                f"Comments: {comments}")

        if brouwer_data['status'] == 'exists':
            if existence:
                return True
            raise RuntimeError(
                f"Andries Brouwer's database claims that such a "
                f"{params}-strongly regular graph exists, but Sage does not "
                f"know how to build it. If *you* do, please get in touch "
                f"with us on sage-devel!\n"
                f"Comments: {comments}")
    if existence:
        return Unknown
    raise RuntimeError(
        f"Sage cannot figure out if a {params}-strongly "
        f"regular graph exists.")


def apparently_feasible_parameters(int n):
    r"""
    Return a list of a priori feasible parameters `(v,k,\lambda,\mu)`, with `0<\mu<k`.

    Note that some of those that it returns may also be infeasible for more
    involved reasons. The condition `0<\mu<k` makes sure we skip trivial cases of
    complete multipartite graphs and their complements.

    INPUT:

    - ``n`` (integer) -- return all a-priori feasible tuples `(v,k,\lambda,\mu)`
      for `v<n`

    EXAMPLES:

    All sets of parameters with `v<20` which pass basic arithmetic tests are
    feasible::

        sage: from sage.graphs.strongly_regular_db import apparently_feasible_parameters
        sage: small_feasible = apparently_feasible_parameters(20); small_feasible
        {(5, 2, 0, 1),
         (9, 4, 1, 2),
         (10, 3, 0, 1),
         (10, 6, 3, 4),
         (13, 6, 2, 3),
         (15, 6, 1, 3),
         (15, 8, 4, 4),
         (16, 5, 0, 2),
         (16, 6, 2, 2),
         (16, 9, 4, 6),
         (16, 10, 6, 6),
         (17, 8, 3, 4)}
        sage: all(graphs.strongly_regular_graph(*x,existence=True) is True              # needs sage.libs.pari
        ....:     for x in small_feasible)
        True

    But that becomes wrong for `v<60` (because of the non-existence of a
    `(49,16,3,6)`-strongly regular graph)::

        sage: small_feasible = apparently_feasible_parameters(60)
        sage: all(graphs.strongly_regular_graph(*x,existence=True) is True              # needs sage.libs.pari
        ....:     for x in small_feasible)
        False
    """
    cdef int v, k, l, mu
    feasible = set()
    for v in range(n):
        for k in range(1, v - 1):
            for l in range(k - 1):
                mu = k*(k - l - 1)//(v - k - 1)
                if mu > 0 and mu < k and seems_feasible(v, k, l, mu):
                    feasible.add((v, k, l, mu))
    return feasible


def _build_small_srg_database():
    r"""
    Build the database of small strongly regular graphs.

    This data is stored in the module-level variable ``_small_srg_database``.
    We use formulas from Cor.3.7 of [CK1986]_ to compute parameters of the
    graph of the projective 2-intersection set associated with a 2-weight code `C`,
    and the usual theory of duality in association schemes to compute the
    parameters of the graph of words of `C`. Another relevant reference is
    Sect.9.8.3 of [BH2012]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import _build_small_srg_database
        sage: _build_small_srg_database()                                               # needs sage.modules sage.rings.finite_rings

    TESTS:

    Make sure that all two-weight codes yield the strongly regular graphs we
    expect::

        sage: graphs.strongly_regular_graph(81, 50, 31, 30)                             # needs sage.libs.pari
        complement(two-intersection set in PG(4,3)): Graph on 81 vertices
        sage: graphs.strongly_regular_graph(243, 220, 199, 200)                 # long time, needs sage.rings.finite_rings
        two-weight code: [55, 5] linear code over GF(3): Graph on 243 vertices
        sage: graphs.strongly_regular_graph(256, 153, 92, 90)                           # needs sage.combinat
        complement(two-intersection set in PG(4,4)): Graph on 256 vertices
        sage: graphs.strongly_regular_graph(256, 170, 114, 110)                         # needs sage.combinat
        complement(two-intersection set in PG(8,2)): Graph on 256 vertices
        sage: graphs.strongly_regular_graph(256, 187, 138, 132)                         # needs sage.combinat
        complement(two-intersection set in PG(8,2)): Graph on 256 vertices
        sage: graphs.strongly_regular_graph(512, 73, 12, 10)                    # not tested (too long), needs sage.rings.finite_rings
        two-weight code: [219, 9] linear code over GF(2): Graph on 512 vertices
        sage: graphs.strongly_regular_graph(512, 219, 106, 84)                  # long time
        two-intersection set in PG(9,2): Graph on 512 vertices
        sage: graphs.strongly_regular_graph(512, 315, 202, 180)                 # not tested (too long), needs sage.rings.finite_rings
        two-weight code: [70, 9] linear code over GF(2): Graph on 512 vertices
        sage: graphs.strongly_regular_graph(625, 364, 213, 210)                 # long time
        complement(two-intersection set in PG(4,5)): Graph on 625 vertices
        sage: graphs.strongly_regular_graph(625, 416, 279, 272)                 # long time
        complement(two-intersection set in PG(4,5)): Graph on 625 vertices
        sage: graphs.strongly_regular_graph(625, 468, 353, 342)                 # long time
        complement(two-intersection set in PG(4,5)): Graph on 625 vertices
        sage: graphs.strongly_regular_graph(729, 336, 153,156)                  # not tested (too long)
        two-intersection set in PG(6,3): Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 420, 243, 240)                 # not tested (too long)
        complement(two-intersection set in PG(6,3)): Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 448, 277, 272)                 # not tested (too long)
        complement(two-intersection set in PG(6,3)): Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 476, 313, 306)                 # not tested (too long)
        complement(two-intersection set in PG(6,3)): Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 532, 391, 380)                 # not tested (too long)
        complement(two-intersection set in PG(6,3)): Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 560, 433, 420)                 # not tested (too long)
        complement(two-intersection set in PG(6,3)): Graph on 729 vertices
        Graph on 729 vertices
        sage: graphs.strongly_regular_graph(729, 616, 523, 506)                 # not tested (too long)
        complement(two-intersection set in PG(6,3)): Graph on 729 vertices
        sage: graphs.strongly_regular_graph(1024, 363, 122, 132)                # not tested (too long)
        two-intersection set in PG(5,4): Graph on 1024 vertices
        sage: graphs.strongly_regular_graph(1024, 396, 148, 156)                # not tested (too long)
        two-intersection set in PG(5,4): Graph on 1024 vertices
        sage: graphs.strongly_regular_graph(1024, 429, 176, 182)                # not tested (too long)
        two-intersection set in PG(5,4): Graph on 1024 vertices
        sage: graphs.strongly_regular_graph(1024, 825, 668, 650)                # not tested (too long)
        complement(two-intersection set in PG(10,2)): Graph on 1024 vertices
    """
    from sage.graphs.generators.smallgraphs import McLaughlinGraph
    from sage.graphs.generators.smallgraphs import CameronGraph
    from sage.graphs.generators.smallgraphs import GritsenkoGraph
    from sage.graphs.generators.smallgraphs import M22Graph
    from sage.graphs.generators.smallgraphs import SimsGewirtzGraph
    from sage.graphs.generators.smallgraphs import HoffmanSingletonGraph
    from sage.graphs.generators.smallgraphs import HigmanSimsGraph
    from sage.graphs.generators.smallgraphs import IoninKharaghani765Graph
    from sage.graphs.generators.smallgraphs import JankoKharaghaniGraph
    from sage.graphs.generators.smallgraphs import LocalMcLaughlinGraph
    from sage.graphs.generators.smallgraphs import SuzukiGraph
    from sage.graphs.generators.smallgraphs import MathonStronglyRegularGraph
    from sage.graphs.generators.smallgraphs import U42Graph216
    from sage.graphs.generators.smallgraphs import U42Graph540

    global _small_srg_database
    _small_srg_database = {
        ( 36,  14,  4,  6): [Graph, ('c~rLDEOcKTPO`U`HOIj@MWFLQFAaRIT`HIWqPsQQJ'
                                     'DXGLqYM@gRLAWLdkEW@RQYQIErcgesClhKefC_ygS'
                                     'GkZ`OyHETdK[?lWStCapVgKK')],
        ( 50,   7,  0,  1): [HoffmanSingletonGraph],
        ( 56,  10,  0,  2): [SimsGewirtzGraph],
        ( 65,  32,  15, 16): [GritsenkoGraph],
        ( 77,  16,   0,  4): [M22Graph],
        (100,  22,   0,  6): [HigmanSimsGraph],
        (100,  44,  18, 20): [SRG_100_44_18_20],
        (100,  45,  20, 20): [SRG_100_45_20_20],
        (105,  32,   4, 12): [SRG_105_32_4_12],
        (120,  63,  30, 36): [SRG_120_63_30_36],
        (120,  77,  52, 44): [SRG_120_77_52_44],
        (126,  25,   8,  4): [SRG_126_25_8_4],
        (126,  50,  13, 24): [SRG_126_50_13_24],
        (144,  39,   6, 12): [SRG_144_39_6_12],
        (162,  56,  10, 24): [LocalMcLaughlinGraph],
        (175,  72,  20, 36): [SRG_175_72_20_36],
        (176,  49,  12, 14): [SRG_176_49_12_14],
        (176,  90,  38, 54): [SRG_176_90_38_54],
        (176, 105,  68, 54): [SRG_176_105_68_54],
        (196,  91,  42, 42): [SRG_196_91_42_42],
        (210,  99,  48, 45): [SRG_210_99_48_45],
        (216,  40,   4,  8): [U42Graph216],
        (220,  84,  38, 28): [SRG_220_84_38_28],
        (231,  30,   9,  3): [CameronGraph],
        (243, 110,  37, 60): [SRG_243_110_37_60],
        (253, 140,  87, 65): [SRG_253_140_87_65],
        (275, 112,  30, 56): [McLaughlinGraph],
        (276, 140,  58, 84): [SRG_276_140_58_84],
        (280, 117, 44,  52): [SRG_280_117_44_52],
        (280, 135,  70, 60): [SRG_280_135_70_60],
        (416, 100,  36, 20): [SRG_416_100_36_20],
        (540, 187,  58, 68): [U42Graph540],
        (560, 208,  72, 80): [SRG_560_208_72_80],
        (630,  85,  20, 10): [SRG_630_85_20_10],
        (765, 192,  48, 48): [IoninKharaghani765Graph],
        (784, 243,  82, 72): [MathonStronglyRegularGraph, 0],
        (784, 270, 98, 90):  [MathonStronglyRegularGraph, 1],
        (784, 297, 116, 110):[MathonStronglyRegularGraph, 2],
        (936, 375, 150,150): [JankoKharaghaniGraph, 936],
        (1288,792, 476,504): [SRG_1288_792_476_504],
        (1782,416, 100, 96): [SuzukiGraph],
        (1800,1029,588,588): [JankoKharaghaniGraph, 1800],
    }

    # Turns the known two-weight codes into SRG constructors
    #
    cdef int n, q, k, w1, w2, K, N, l, m, K_O, l_O, m_O
    import sage.coding.two_weight_db
    from sage.matrix.constructor import matrix
    from sage.rings.integer_ring import ZZ
    cinv = matrix(ZZ, [[1, 0, 0], [0, 0, 1], [0, 1, 0]])
    for code in sage.coding.two_weight_db.data:
        n, q, k, w1, w2 = code['n'], code['K'].cardinality(), code['k'], code['w1'], code['w2']
        N = q**k
        K_O = n*(q - 1)
        l_O = K_O**2 + 3*K_O - q*(w1 + w2) - K_O*q*(w1 + w2) + w1*w2*q**2
        m_O = (w1*w2*q**2)//N

        em = eigenmatrix(N, K_O, l_O, m_O)  # 1st eigenmatrix
        assert((em is not None) and (em.det() != 0))
        emi = N*em.inverse()                # 2nd eigenmatrix
        # 1st and 2nd eigenmatrices equal up to renumbering graphs?
        selfdual = em == cinv*emi*cinv
        _small_srg_database[N, K_O, l_O, m_O] = \
            [lambda x: strongly_regular_from_two_intersection_set(x.transpose()), code['M']]
        if not selfdual:  # we can build two graphs (not complements to each other!)
            K, s, r = emi[0, 1], emi[1, 1], emi[2, 1]  # by Thm 5.7 in [CK1986]_.
            l = K + r*s + r + s
            m = K + r*s
            _small_srg_database[N, K, l, m] = [strongly_regular_from_two_weight_code, code['M']]


cdef load_brouwer_database():
    r"""
    Loads Andries Brouwer's database into _brouwer_database.
    """
    global _brouwer_database
    if _brouwer_database is not None:
        return

    from sage.env import GRAPHS_DATA_DIR
    filename = os.path.join(GRAPHS_DATA_DIR, 'brouwer_srg_database.json')
    with open(filename) as fobj:
        database = json.load(fobj)

    _brouwer_database = {
        (v, k, l, mu): {'status': status, 'comments': comments}
        for (v, k, l, mu, status, comments) in database
    }


def _check_database():
    r"""
    Checks the coherence of Andries Brouwer's database with Sage.

    The function also outputs some statistics on the database.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import _check_database
        sage: _check_database() # long time
        Sage cannot build a (512  133  24   38  ) that exists. Comment ...
        ...
        In Andries Brouwer's database:
        - 462 impossible entries
        - 2911 undecided entries
        - 1165 realizable entries (Sage misses ... of them)

    """
    global _brouwer_database
    load_brouwer_database()

    # Check that all parameters detected as infeasible are actually infeasible
    # in Brouwer's database, for a test that was implemented.
    for params in set(_brouwer_database).difference(apparently_feasible_parameters(1301)):
        if _brouwer_database[params]['status'] != "impossible":
            raise RuntimeError("Brouwer's db does not seem to know that {} in unfeasible".format(params))
        comment = _brouwer_database[params]['comments']
        if ('Krein'    in comment or
            'Absolute' in comment or
            'Conf'     in comment or
            'mu=1'     in comment or
            '&mu;=2'   in comment):
            continue
        raise RuntimeError("We detected that {} was unfeasible, but maybe we should not have".format(params))

    # We empty the global database, to be sure that strongly_regular_graph does
    # not use its data to answer.
    _brouwer_database, saved_database = {}, _brouwer_database

    cdef int missed = 0
    for params, dic in sorted(saved_database.items()):
        sage_answer = strongly_regular_graph(*params, existence=True)
        if dic['status'] == 'open':
            if sage_answer is True:
                print("Sage can build a {}, Brouwer's database cannot".format(params))
            assert sage_answer is not False
        elif dic['status'] == 'exists':
            if sage_answer is not True:
                print(("Sage cannot build a ({:<4} {:<4} {:<4} {:<4}) that exists. " +
                       "Comment from Brouwer's database: ").format(*params)
                      + dic['comments'])
                missed += 1
            assert sage_answer is not False
        elif dic['status'] == 'impossible':
            assert sage_answer is not True
        else:
            assert False  # must not happen

    status = [x['status'] for x in saved_database.values()]
    print("\nIn Andries Brouwer's database:")
    print("- {} impossible entries".format(status.count('impossible')))
    print("- {} undecided entries".format(status.count('open')))
    print("- {} realizable entries (Sage misses {} of them)".format(status.count('exists'), missed))

    # Reassign its value to the global database
    _brouwer_database = saved_database
