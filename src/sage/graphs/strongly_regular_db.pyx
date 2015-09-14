# -*- coding: utf-8 -*-
r"""
Database of strongly regular graphs

This module manages a database associating to a set of four integers
`(v,k,\lambda,\mu)` a strongly regular graphs with these parameters, when one
exists.

Using Andries Brouwer's `database of strongly regular graphs
<http://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__, it can also return
non-existence results. Note that some constructions are missing, and that some
strongly regular graphs that exist in the database cannot be automatically built
by Sage. Help us if you know any.

.. NOTE::

    Any missing/incorrect information in the database must be reported to
    `Andries E. Brouwer <http://www.win.tue.nl/~aeb/>`__ directly, in order to
    have a unique and updated source of information.

REFERENCES:

.. [BvL84] A. Brouwer, J van Lint,
   Strongly regular graphs and partial geometries,
   Enumeration and design,
   (Waterloo, Ont., 1982) (1984): 85-122.
   http://oai.cwi.nl/oai/asset/1817/1817A.pdf

Functions
---------
"""
from sage.categories.sets_cat import EmptySetError
from sage.misc.unknown import Unknown
from sage.rings.arith import is_square
from sage.rings.arith import is_prime_power
from sage.misc.cachefunc import cached_function
from sage.combinat.designs.orthogonal_arrays import orthogonal_array
from sage.combinat.designs.bibd import balanced_incomplete_block_design
from sage.graphs.generators.smallgraphs import McLaughlinGraph
from sage.graphs.generators.smallgraphs import CameronGraph
from sage.graphs.generators.smallgraphs import M22Graph
from sage.graphs.generators.smallgraphs import SimsGewirtzGraph
from sage.graphs.generators.smallgraphs import HoffmanSingletonGraph
from sage.graphs.generators.smallgraphs import SchlaefliGraph
from sage.graphs.generators.smallgraphs import HigmanSimsGraph
from sage.graphs.generators.smallgraphs import LocalMcLaughlinGraph
from sage.graphs.graph import Graph
from libc.math cimport sqrt
from sage.matrix.constructor import Matrix
from sage.rings.finite_rings.constructor import FiniteField as GF
from sage.coding.linear_code import LinearCode
from sage.rings.sum_of_squares cimport two_squares_c
from libc.stdint cimport uint_fast32_t

cdef dict _brouwer_database = None

@cached_function
def is_paley(int v,int k,int l,int mu):
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
        sage: g = t[0](*t[1:]); g
        Paley graph with parameter 13: Graph on 13 vertices
        sage: g.is_strongly_regular(parameters=True)
        (13, 6, 2, 3)
        sage: t = is_paley(5,5,5,5); t
    """
    if (v%4 == 1 and is_prime_power(v) and
        k   == (v-1)/2 and
        l   == (v-5)/4 and
        mu  == (v-1)/4):
        from sage.graphs.generators.families import PaleyGraph
        return (lambda q : PaleyGraph(q),v)

@cached_function
def is_orthogonal_array_block_graph(int v,int k,int l,int mu):
    r"""
    Test whether some Orthogonal Array graph is `(v,k,\lambda,\mu)`-strongly regular.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_orthogonal_array_block_graph
        sage: t = is_orthogonal_array_block_graph(64, 35, 18, 20); t
        (..., 5, 8)
        sage: g = t[0](*t[1:]); g
        OA(5,8): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 35, 18, 20)

        sage: t = is_orthogonal_array_block_graph(5,5,5,5); t
    """
    # notations from
    # http://www.win.tue.nl/~aeb/graphs/OA.html
    if not is_square(v):
        return
    n = int(sqrt(v))
    if k % (n-1):
        return
    m = k//(n-1)
    if (l  != (m-1)*(m-2)+n-2 or
        mu != m*(m-1)):
        return
    if orthogonal_array(m,n,existence=True):
        from sage.graphs.generators.intersection import OrthogonalArrayBlockGraph
        return (lambda m,n : OrthogonalArrayBlockGraph(m, n), m,n)

@cached_function
def is_johnson(int v,int k,int l,int mu):
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
    # Using notations of http://www.win.tue.nl/~aeb/graphs/Johnson.html
    #
    # J(n,m) has parameters v = m(m – 1)/2, k = 2(m – 2), λ = m – 2, μ = 4.
    m = l + 2
    if (mu == 4 and
        k  == 2*(m-2) and
        v  == m*(m-1)/2):
        from sage.graphs.generators.families import JohnsonGraph
        return (lambda m: JohnsonGraph(m,2), m)

@cached_function
def is_steiner(int v,int k,int l,int mu):
    r"""
    Test whether some Steiner graph is `(v,k,\lambda,\mu)`-strongly regular.

    A Steiner graph is the intersection graph of a Steiner set system. For more
    information, see http://www.win.tue.nl/~aeb/graphs/S.html.

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
    # Using notations from http://www.win.tue.nl/~aeb/graphs/S.html
    #
    # The block graph of a Steiner 2-design S(2,m,n) has parameters:
    # v = n(n-1)/m(m-1), k = m(n-m)/(m-1), λ = (m-1)^2 + (n-1)/(m–1)–2, μ = m^2.
    if mu <= 1 or not is_square(mu):
        return
    m = int(sqrt(mu))
    n = (k*(m-1))//m+m
    if (v == (n*(n-1))/(m*(m-1)) and
        k == m*(n-m)/(m-1) and
        l == (m-1)**2 + (n-1)/(m-1)-2 and
        balanced_incomplete_block_design(n,m,existence=True)):
        from sage.graphs.generators.intersection import IntersectionGraph
        return (lambda n,m: IntersectionGraph(map(frozenset,balanced_incomplete_block_design(n,m))),n,m)

@cached_function
def is_affine_polar(int v,int k,int l,int mu):
    r"""
    Test whether some Affine Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see http://www.win.tue.nl/~aeb/graphs/VO.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_affine_polar
        sage: t = is_affine_polar(81,32,13,12); t
        (..., 4, 3)
        sage: g = t[0](*t[1:]); g
        Affine Polar Graph VO^+(4,3): Graph on 81 vertices
        sage: g.is_strongly_regular(parameters=True)
        (81, 32, 13, 12)

        sage: t = is_affine_polar(5,5,5,5); t
    """
    from sage.rings.arith import divisors
    # Using notations from http://www.win.tue.nl/~aeb/graphs/VO.html
    #
    # VO+(2e,q) has parameters: v = q^(2e), k = (q^(e−1) + 1)(q^e − 1), λ =
    # q(q^(e−2) + 1)(q^(e−1) − 1) + q − 2, μ = q^(e−1)(q^(e−1) + 1)
    #
    # VO−(2e,q) has parameters v = q^(2e), k = (q^(e−1) - 1)(q^e + 1), λ =
    # q(q^(e−2) - 1)(q^(e−1) + 1) + q − 2, μ = q^(e−1)(q^(e−1) - 1)
    if (not is_square(v) or
        not is_prime_power(v)):
        return
    prime,power = is_prime_power(v,get_data=True)
    if power%2:
        return
    for e in divisors(power/2):
        q = prime**(power//(2*e))
        assert v == q**(2*e)
        if (k == (q**(e-1) + 1)*(q**e-1) and
            l == q*(q**(e-2) + 1)*(q**(e-1)-1)+q-2 and
            mu== q**(e-1)*(q**(e-1) + 1)):
            from sage.graphs.generators.families import AffineOrthogonalPolarGraph
            return (lambda d,q : AffineOrthogonalPolarGraph(d,q,sign='+'),2*e,q)
        if (k == (q**(e-1) - 1)*(q**e+1) and
            l == q*(q**(e-2)- 1)*(q**(e-1)+1)+q-2 and
            mu== q**(e-1)*(q**(e-1) - 1)):
            from sage.graphs.generators.families import AffineOrthogonalPolarGraph
            return (lambda d,q : AffineOrthogonalPolarGraph(d,q,sign='-'),2*e,q)

@cached_function
def is_orthogonal_polar(int v,int k,int l,int mu):
    r"""
    Test whether some Orthogonal Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see http://www.win.tue.nl/~aeb/graphs/srghub.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_orthogonal_polar
        sage: t = is_orthogonal_polar(85, 20, 3, 5); t
        (<function OrthogonalPolarGraph at ...>, 5, 4, '')
        sage: g = t[0](*t[1:]); g
        Orthogonal Polar Graph O(5, 4): Graph on 85 vertices
        sage: g.is_strongly_regular(parameters=True)
        (85, 20, 3, 5)

        sage: t = is_orthogonal_polar(5,5,5,5); t

    TESTS:

    All of ``O(2m+1,q)``, ``O^+(2m,q)`` and ``O^-(2m,q)`` appear::

        sage: is_orthogonal_polar(85, 20, 3, 5)
        (<function OrthogonalPolarGraph at ...>, 5, 4, '')
        sage: is_orthogonal_polar(119,54,21,27)
        (<function OrthogonalPolarGraph at ...>, 8, 2, '-')
        sage: is_orthogonal_polar(130,48,20,16)
        (<function OrthogonalPolarGraph at ...>, 6, 3, '+')

    """
    from sage.rings.arith import divisors
    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return
    q_pow_m_minus_one = -s-1 if abs(s) > r else r+1

    if is_prime_power(q_pow_m_minus_one):
        prime,power = is_prime_power(q_pow_m_minus_one,get_data=True)
        for d in divisors(power):
            q = prime**d
            m = (power//d)+1

            # O(2m+1,q)
            if (v == (q**(2*m)-1)/(q-1)              and
                k == q*(q**(2*m-2)-1)/(q-1)          and
                l == q**2*(q**(2*m-4)-1)/(q-1) + q-1 and
                mu== (q**(2*m-2)-1)/(q-1)):
                from sage.graphs.generators.families import OrthogonalPolarGraph
                return (OrthogonalPolarGraph, 2*m+1, q, "")

            # O^+(2m,q)
            if (v ==   (q**(2*m-1)-1)/(q-1) + q**(m-1)   and
                k == q*(q**(2*m-3)-1)/(q-1) + q**(m-1) and
                k == q**(2*m-3) + l + 1                  and
                mu== k/q):
                from sage.graphs.generators.families import OrthogonalPolarGraph
                return (OrthogonalPolarGraph, 2*m, q, "+")

            # O^+(2m+1,q)
            if (v ==   (q**(2*m-1)-1)/(q-1) - q**(m-1)   and
                k == q*(q**(2*m-3)-1)/(q-1) - q**(m-1) and
                k == q**(2*m-3) + l + 1                  and
                mu== k/q):
                from sage.graphs.generators.families import OrthogonalPolarGraph
                return (OrthogonalPolarGraph, 2*m, q, "-")

def is_RSHCD(int v,int k,int l,int mu):
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
        sage: t = is_RSHCD(64,27,10,12); t
        [<built-in function SRG_from_RSHCD>, 64, 27, 10, 12]
        sage: g = t[0](*t[1:]); g
        Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 27, 10, 12)

    """
    if SRG_from_RSHCD(v,k,l,mu,existence=True):
        return [SRG_from_RSHCD,v,k,l,mu]

def SRG_from_RSHCD(v,k,l,mu, existence=False,check=True):
    r"""
    Return a `(v,k,l,mu)`-strongly regular graph from a RSHCD

    This construction appears in 8.D of [BvL84]_. For more information, see
    :func:`~sage.combinat.matrices.hadamard_matrix.regular_symmetric_hadamard_matrix_with_constant_diagonal`.

    INPUT:

    - ``v,k,l,mu`` (integers)

    - ``existence`` (boolean) -- whether to return a graph or to test if Sage
      can build such a graph.

    - ``check`` (boolean) -- whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to ``True``
      by default.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_from_RSHCD
        sage: SRG_from_RSHCD(784, 0, 14, 38, existence=True)
        False
        sage: SRG_from_RSHCD(784, 377, 180, 182, existence=True)
        True
        sage: SRG_from_RSHCD(144, 65, 28, 30)
        Graph on 144 vertices

    TESTS::

        sage: SRG_from_RSHCD(784, 0, 14, 38)
        Traceback (most recent call last):
        ...
        ValueError: I do not know how to build a (784, 0, 14, 38)-SRG from a RSHCD

    """
    from sage.combinat.matrices.hadamard_matrix import regular_symmetric_hadamard_matrix_with_constant_diagonal
    sgn = lambda x: 1 if x>=0 else -1
    n = v
    a = (n-4*mu)//2
    e = 2*k - n + 1 + a
    t = abs(a//2)

    if (e**2 == 1              and
        k == (n-1-a+e)/2       and
        l == (n-2*a)/4 - (1-e) and
        mu== (n-2*a)/4         and
        regular_symmetric_hadamard_matrix_with_constant_diagonal(n,sgn(a)*e,existence=True)):
        if existence:
            return True
        from sage.matrix.constructor import identity_matrix as I
        from sage.matrix.constructor import ones_matrix     as J

        H = regular_symmetric_hadamard_matrix_with_constant_diagonal(n,sgn(a)*e)
        if list(H.column(0)[1:]).count(1) == k:
            H = -H
        G = Graph((J(n)-I(n)-H+H[0,0]*I(n))/2,loops=False,multiedges=False,format="adjacency_matrix")
        if check:
            assert G.is_strongly_regular(parameters=True) == (v,k,l,mu)
        return G

    if existence:
        return False
    raise ValueError("I do not know how to build a {}-SRG from a RSHCD".format((v,k,l,mu)))

@cached_function
def is_unitary_polar(int v,int k,int l,int mu):
    r"""
    Test whether some Unitary Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see http://www.win.tue.nl/~aeb/graphs/srghub.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_unitary_polar
        sage: t = is_unitary_polar(45, 12, 3, 3); t
        (<function UnitaryPolarGraph at ...>, 4, 2)
        sage: g = t[0](*t[1:]); g
        Unitary Polar Graph U(4, 2); GQ(4, 2): Graph on 45 vertices
        sage: g.is_strongly_regular(parameters=True)
        (45, 12, 3, 3)

        sage: t = is_unitary_polar(5,5,5,5); t

    TESTS:

    All the ``U(n,q)`` appear::

        sage: t = is_unitary_polar(45, 12, 3, 3); t
        (<function UnitaryPolarGraph at ...>, 4, 2)
        sage: t = is_unitary_polar(165, 36, 3, 9); t
        (<function UnitaryPolarGraph at ...>, 5, 2)
        sage: t = is_unitary_polar(693, 180, 51, 45); t
        (<function UnitaryPolarGraph at ...>, 6, 2)
        sage: t = is_unitary_polar(1105, 80, 15, 5); t
        (<function UnitaryPolarGraph at ...>, 4, 4)
    """
    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return
    q = k/mu
    if q*mu != k or q < 2:
        return
    p,t = is_prime_power(q, get_data=True)
    if p**t != q or t % 2 != 0:
        return
    # at this point we know that we should have U(n,q) for some n and q=p^t, t even
    if r > 0:
        q_pow_d_minus_one = r+1
    else:
        q_pow_d_minus_one = -s-1
    ppp,ttt = is_prime_power(q_pow_d_minus_one, get_data=True)
    d = ttt/t + 1
    if ppp != p or (d-1)*t != ttt:
        return
    t /= 2
    # U(2d+1,q); write q^(1/2) as p^t
    if (v == (q**d - 1)*((q**d)*p**t + 1)/(q - 1)               and
        k == q*(q**(d-1) - 1)*((q**d)/(p**t) + 1)/(q - 1)       and
        l == q*q*(q**(d-2)-1)*((q**(d-1))/(p**t) + 1)/(q - 1) + q - 1):
        from sage.graphs.generators.families import UnitaryPolarGraph
        return (UnitaryPolarGraph, 2*d+1, p**t)

    # U(2d,q);
    if (v == (q**d - 1)*((q**d)/(p**t) + 1)/(q - 1)             and
        k == q*(q**(d-1) - 1)*((q**(d-1))/(p**t) + 1)/(q - 1)   and
        l == q*q*(q**(d-2)-1)*((q**(d-2))/(p**t) + 1)/(q - 1) + q - 1):
        from sage.graphs.generators.families import UnitaryPolarGraph
        return (UnitaryPolarGraph, 2*d, p**t)

@cached_function
def is_unitary_dual_polar(int v,int k,int l,int mu):
    r"""
    Test whether some Unitary Dual Polar graph is `(v,k,\lambda,\mu)`-strongly regular.

    This must be the U_5(q) on totally isotropic lines.
    For more information, see http://www.win.tue.nl/~aeb/graphs/srghub.html.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_unitary_dual_polar
        sage: t = is_unitary_dual_polar(297, 40, 7, 5); t
        (<function UnitaryDualPolarGraph at ...>, 5, 2)
        sage: g = t[0](*t[1:]); g
        Unitary Dual Polar Graph DU(5, 2); GQ(8, 4): Graph on 297 vertices
        sage: g.is_strongly_regular(parameters=True)
        (297, 40, 7, 5)
        sage: t = is_unitary_dual_polar(5,5,5,5); t

    TESTS::

        sage: is_unitary_dual_polar(6832, 270, 26, 10)
        (<function UnitaryDualPolarGraph at ...>, 5, 3)
    """
    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return
    q = mu - 1
    if q < 2:
        return
    p,t = is_prime_power(q, get_data=True)
    if p**t != q or t % 2 != 0:
        return
    if (r < 0 and q != -r - 1) or (s < 0 and q != -s - 1):
       return
    t /= 2
    # we have correct mu, negative eigenvalue, and q=p^(2t)
    if (v == (q**2*p**t + 1)*(q*p**t + 1)  and
        k == q*p**t*(q + 1)                and
        l == k - 1 - q**2*p**t):
        from sage.graphs.generators.families import UnitaryDualPolarGraph
        return (UnitaryDualPolarGraph, 5, p**t)

@cached_function
def is_twograph_descendant_of_srg(int v, int k0, int l, int mu):
    r"""
    Test whether some descendant graph of an s.r.g. is `(v,k_0,\lambda,\mu)`-s.r.g.

    We check whether there can exist `(v+1,k,\lambda^*,\mu^*)`-s.r.g. `G` so
    that ``self`` is a descendant graph of the regular two-graph specified
    by `G`.
    Specifically, we must have that `v+1=2(2k-\lambda^*-\mu^*)`, and
    `k_0=2(k-\mu^*)`, `\lambda=k+\lambda^*-2\mu^*`, `\mu=k-\mu^*`, which give 2
    independent linear conditions, say `k-\mu^*=\mu` and
    `\lambda^*-\mu^*=\lambda-\mu`.  Further, there is a quadratic relation
    `2 k^2-(v+1+4 \mu) k+ 2 v \mu=0`.

    If we can contruct such `G` then we return a function to build a
    `(v,k_0,\lambda,\mu)`-s.r.g.  For more information,
    see 10.3 in http://www.win.tue.nl/~aeb/2WF02/spectra.pdf

    INPUT:

    - ``v,k0,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph if one
    exists and is known, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_twograph_descendant_of_srg
        sage: t = is_twograph_descendant_of_srg(27, 10, 1, 5); t
        (<cyfunction is_twograph_descendant_of_srg.<locals>.la at...
        sage: g = t[0](*t[1:]); g
        descendant of complement(Johnson graph with parameters 8,2) at {5, 7}: Graph on 27 vertices
        sage: g.is_strongly_regular(parameters=True)
        (27, 10, 1, 5)
        sage: t = is_twograph_descendant_of_srg(5,5,5,5); t

    TESTS::

        sage: graphs.strongly_regular_graph(279, 150, 85, 75, existence=True)
        True
        sage: graphs.strongly_regular_graph(279, 150, 85, 75).is_strongly_regular(parameters=True) # optional - gap_packages internet
        (279, 150, 85, 75)
    """
    cdef int b, k, s
    if k0 != 2*mu or v % 2 == 0:
        return
    b = v+1+4*mu
    D = sqrt(b**2-16*v*mu)
    if int(D)==D:
        for kf in [(-D+b)/4, (D+b)/4]:
            k = int(kf)
            if k == kf and \
                strongly_regular_graph(v+1, k, l - 2*mu + k , k - mu,  existence=True):
                def la(vv):
                    from sage.combinat.designs.twographs import twograph_descendant
                    g = strongly_regular_graph(vv, k, l - 2*mu + k)
                    return twograph_descendant(g, g.vertex_iterator().next(), name=True)
                return(la, v+1)
    return

@cached_function
def is_taylor_twograph_srg(int v,int k,int l,int mu):
    r"""
    Test whether some Taylor two-graph SRG is `(v,k,\lambda,\mu)`-strongly regular.

    For more information, see §7E of [BvL84]_.

    INPUT:

    - ``v,k,l,mu`` (integers)

    OUTPUT:

    A tuple ``t`` such that ``t[0](*t[1:])`` builds the requested graph
    :func:`TaylorTwographSRG <sage.graphs.generators.families.TaylorTwographSRG>`
    if the parameters match, and ``None`` otherwise.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import is_taylor_twograph_srg
        sage: t = is_taylor_twograph_srg(28, 15, 6, 10); t
        (<function TaylorTwographSRG at ...>, 3)
        sage: g = t[0](*t[1:]); g
        Taylor two-graph SRG: Graph on 28 vertices
        sage: g.is_strongly_regular(parameters=True)
        (28, 15, 6, 10)
        sage: t = is_taylor_twograph_srg(5,5,5,5); t

    TESTS::

        sage: is_taylor_twograph_srg(730, 369, 168, 205)
        (<function TaylorTwographSRG at ...>, 9)

    """
    r,s = eigenvalues(v,k,l,mu)
    if r is None:
        return
    p,t = is_prime_power(v-1, get_data=True)
    if p**t+1 != v or t % 3 != 0 or p % 2 == 0:
        return
    q = p**(t//3)
    if (k, l, mu) == (q*(q**2+1)/2, (q**2+3)*(q-1)/4, (q**2+1)*(q+1)/4):
        from sage.graphs.generators.families import TaylorTwographSRG
        return (TaylorTwographSRG, q)
    return

cdef eigenvalues(int v,int k,int l,int mu):
    r"""
    Return the eigenvalues of a (v,k,l,mu)-strongly regular graph.

    If the set of parameters is not feasible, or if they correspond to a
    conference graph, the function returns ``(None,None)``.

    INPUT:

    - ``v,k,l,mu`` (integers)

    """
    # See 1.3.1 of [Distance-regular graphs]
    b = (mu-l)
    c = (mu-k)
    D = b**2-4*c
    if not is_square(D):
        return [None,None]
    return [(-b+sqrt(D))/2.0,
            (-b-sqrt(D))/2.0]

def _H_3_cayley_graph(L):
    r"""
    return the `L`-Cayley graph of the group `H_3` from Prop. 12 in [JK03]_.

    INPUT:

    - the list of words for the generating set in the format ["abc",...,"xyz"] for
      a,b,...,z being integers between 0 and 4.

    TESTS::

        sage: from sage.graphs.strongly_regular_db import _H_3_cayley_graph
        sage: _H_3_cayley_graph(["100","110","130","140","200","230","240","300"])
        Graph on 100 vertices
    """
    from sage.groups.free_group import FreeGroup
    from sage.groups.finitely_presented import FinitelyPresentedGroup
    G = FreeGroup('x,y,z')
    x,y,z = G.gens()
    rels = (x**5,y**5,z**4,x*y*x**(-1)*y**(-1),z*x*z**(-1)*x**(-2),z*y*z**(-1)*y**(-2))
    G = FinitelyPresentedGroup(G,rels)
    x,y,z = G.gens()
    H = G.as_permutation_group()
    L = map(lambda x:map(int,x),L)
    x,y,z=(H.gen(0),H.gen(1),H.gen(2))
    L = [H(x**xx*y**yy*z**zz) for xx,yy,zz in L]
    return Graph(H.cayley_graph(generators=L, simple=True))

def SRG_100_44_18_20():
    r"""
    Return a `(100, 44, 18, 20)`-strongly regular graph.

    This graph is built as a Cayley graph, using the construction for `\Delta_1`
    with group `H_3` presented in Table 8.1 of [JK03]_

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_100_44_18_20
        sage: G = SRG_100_44_18_20()                 # long time
        sage: G.is_strongly_regular(parameters=True) # long time
        (100, 44, 18, 20)

    REFERENCES:

    .. [JK03] L. K. Jørgensen, M. Klin, M.,
      Switching of edges in strongly regular graphs.
      I. A family of partial difference sets on 100 vertices,
      Electronic Journal of Combinatorics 10(1), 2003.
    """
    return _H_3_cayley_graph(["100","110","130","140","200","230","240","300",
             "310","320","400","410","420","440","041","111","221","231","241",
             "321","331","401","421","441","002","042","112","122","142","212",
             "232","242","322","342","033","113","143","223","303","333","343",
             "413","433","443"])

def SRG_100_45_20_20():
    r"""
    Return a `(100, 45, 20, 20)`-strongly regular graph.

    This graph is built as a Cayley graph, using the construction for `\Gamma_3`
    with group `H_3` presented in Table 8.1 of [JK03]_.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_100_45_20_20
        sage: G = SRG_100_45_20_20()              # long time
        sage: G.is_strongly_regular(parameters=True) # long time
        (100, 45, 20, 20)
    """
    return _H_3_cayley_graph(["120","140","200","210","201","401","411","321",
             "002","012","022","042","303","403","013","413","240","031","102",
             "323","300","231","132","133","310","141","142","233","340","241",
             "202","333","410","341","222","433","430","441","242","302","312",
             "322","332","442","143"])

def SRG_120_77_52_44():
    r"""
    Return a `(120,77,52,44)`-strongly regular graph.

    To build this graph, we first build a `2-(21,7,12)` design, by removing two
    points from the :func:`~sage.combinat.designs.block_design.WittDesign` on 23
    points. We then build the intersection graph of blocks with intersection
    size 3.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_120_77_52_44
        sage: G = SRG_120_77_52_44()                 # optional - gap_packages
        sage: G.is_strongly_regular(parameters=True) # optional - gap_packages
        (120, 77, 52, 44)
    """
    from sage.combinat.designs.block_design import WittDesign
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    W = WittDesign(23)
    H = IncidenceStructure([x for x in W if 22 not in x and 21 not in x])
    return H.intersection_graph(3)

def SRG_176_105_68_54():
    r"""
    Return a `(176, 105, 68, 54)`-strongly regular graph.

    To build this graph, we first build a `2-(22,7,16)` design, by removing one
    point from the :func:`~sage.combinat.designs.block_design.WittDesign` on 23
    points. We then build the intersection graph of blocks with intersection
    size 3.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_176_105_68_54
        sage: G = SRG_176_105_68_54()                # optional - gap_packages
        sage: G.is_strongly_regular(parameters=True) # optional - gap_packages
        (176, 105, 68, 54)
    """
    from sage.combinat.designs.block_design import WittDesign
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    W = WittDesign(23)
    H = IncidenceStructure([x for x in W if 22 not in x])
    return H.intersection_graph(3)

def SRG_253_140_87_65():
    r"""
    Return a `(253, 140, 87, 65)`-strongly regular graph.

    To build this graph, we first build the
    :func:`~sage.combinat.designs.block_design.WittDesign` on 23 points which is
    a `2-(23,7,21)` design. We then build the intersection graph of blocks with
    intersection size 3.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_253_140_87_65
        sage: G = SRG_253_140_87_65()                # optional - gap_packages
        sage: G.is_strongly_regular(parameters=True) # optional - gap_packages
        (253, 140, 87, 65)
    """
    from sage.combinat.designs.block_design import WittDesign
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    W = WittDesign(23)
    return W.intersection_graph(3)

def SRG_196_91_42_42():
    r"""
    Return a `(196,91,42,42)`-strongly regular graph.

    This strongly regular graph is built following the construction provided in
    Corollary 8.2.27 of [IS06]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_196_91_42_42
        sage: G = SRG_196_91_42_42()
        sage: G.is_strongly_regular(parameters=True)
        (196, 91, 42, 42)

    REFERENCE:

    .. [IS06] Y.J. Ionin, S. Shrikhande,
      Combinatorics of symmetric designs.
      Cambridge University Press, 2006.
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    from sage.graphs.generators.intersection import IntersectionGraph
    k = 7
    G = IntegerModRing(91)
    A = map(G,{0, 10, 27, 28, 31, 43, 50})
    B = map(G,{0, 11, 20, 25, 49, 55, 57})
    H = map(G,[13*i for i in range(k)])
    U = map(frozenset,[[x+z for x in A] for z in G])
    V = map(frozenset,[[x+z for x in B] for z in G])
    W = map(frozenset,[[x+z for x in H] for z in G])
    G = IntersectionGraph(U+V+W)

    G.seidel_switching(U)

    G.add_edges((-1,x) for x in U)
    G.relabel()
    return G

def SRG_280_135_70_60():
    r"""
    Return a strongly regular graph with parameters (280, 135, 70, 60).

    This graph is built from the action of `J_2` on a `3.PGL(2,9)` subgroup it
    contains.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_280_135_70_60
        sage: g=SRG_280_135_70_60()                  # long time # optional - gap_packages
        sage: g.is_strongly_regular(parameters=True) # long time # optional - gap_packages
        (280, 135, 70, 60)
    """
    from sage.interfaces.gap import gap
    from sage.groups.perm_gps.permgroup import PermutationGroup
    from sage.graphs.graph import Graph

    gap.load_package("AtlasRep")

    # A representation of J2 acting on a 3.PGL(2,9) it contains.
    J2    = PermutationGroup(gap('AtlasGenerators("J2",2).generators'))
    edges = J2.orbit((1,2),"OnSets")
    g     = Graph()
    g.add_edges(edges)
    g.relabel()
    return g

def strongly_regular_from_two_weight_code(L):
    r"""
    Return a strongly regular graph from a two-weight code.

    A code is said to be a *two-weight* code the weight of its nonzero codewords
    (i.e. their number of nonzero coordinates) can only be one of two integer
    values `w_1,w_2`. It is said to be *projective* if the minimum weight of the
    dual code is `\geq 3`. A strongly regular graph can be built from a
    two-weight projective code with weights `w_1,w_2` (assuming `w_1<w_2`) by
    adding an edge between any two codewords whose difference has weight
    `w_1`. For more information, see [vLintSchrijver81]_ or [Delsarte72]_.

    INPUT:

    - ``L`` -- a two-weight linear code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import strongly_regular_from_two_weight_code
        sage: x=("100022021001111",
        ....:    "010011211122000",
        ....:    "001021112100011",
        ....:    "000110120222220")
        sage: M = Matrix(GF(3),[list(l) for l in x])
        sage: G = strongly_regular_from_two_weight_code(LinearCode(M))
        sage: G.is_strongly_regular(parameters=True)
        (81, 50, 31, 30)

    REFERENCES:

    .. [vLintSchrijver81] J. H. van Lint, and A. Schrijver (1981),
      Construction of strongly regular graphs, two-weight codes and
      partial geometries by finite fields,
      Combinatorica, 1(1), 63-73.

    .. [Delsarte72] Ph. Delsarte,
      Weights of linear codes and strongly regular normed spaces,
      Discrete Mathematics (1972), Volume 3, Issue 1, Pages 47-64,
      http://dx.doi.org/10.1016/0012-365X(72)90024-6.

    """
    V = map(tuple,list(L))
    w1, w2 = sorted(set(sum(map(bool,x)) for x in V).difference([0]))
    G = Graph([V,lambda u,v: sum(uu!=vv for uu,vv in zip(u,v)) == w1])
    G.relabel()
    return G

def SRG_256_187_138_132():
    r"""
    Return a `(256, 187, 138, 132)`-strongly regular graph.

    This graph is built from a projective binary `[68,8]` code with weights `32,
    40`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_256_187_138_132
        sage: G = SRG_256_187_138_132()
        sage: G.is_strongly_regular(parameters=True)
        (256, 187, 138, 132)
    """
    x=("10000000100111100110000001101000100111000011100101011010111111010110",
       "01000000010011110011000000110100010011100001110010101101011111101011",
       "00100000001001111101100000011010001001110000111001010110101111110101",
       "00010000100011011100110001100101100011111011111001100001101000101100",
       "00001000110110001100011001011010011110111110011001111010001011000000",
       "00000100111100100000001101000101101000011100101001110111111010110110",
       "00000010011110010000000110100010111100001110010100101011111101011011",
       "00000001001111001100000011010001011110000111001010010101111110101101")
    M = Matrix(GF(2),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_729_532_391_380():
    r"""
    Return a `(729, 532, 391, 380)`-strongly regular graph.

    This graph is built from a projective ternary `[98,6]` code with weights
    `63, 72`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_729_532_391_380
        sage: G = SRG_729_532_391_380()               # long time
        sage: G.is_strongly_regular(parameters=True)  # long time
        (729, 532, 391, 380)
    """
    x=("10000021022112121121110122000110112002010011100120022110120200120111220220122120012012100201110210",
       "01000020121020200200211101202121120002211002210100021021202220112122012212101102010210010221221201",
       "00100021001211011111111202120022221002201111021101021212210122101020121111002000210000101222202000",
       "00010022122200222202201212211112001102200112202202121201211212010210202001222120000002110021000110",
       "00001021201002010011020210221221012112200012020011201200111021021102212120211102012002011201210221",
       "00000120112212122122202110022202210010200022002120112200101002202221111102110100210212001022201202")
    M = Matrix(GF(3),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_729_560_433_420():
    r"""
    Return a `(729, 560, 433, 420)`-strongly regular graph.

    This graph is built from a projective ternary `[84,6]` code with weights
    `54, 63`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_729_560_433_420
        sage: G = SRG_729_560_433_420()               # long time
        sage: G.is_strongly_regular(parameters=True) # long time
        (729, 560, 433, 420)
    """
    x=("100000210221121211211212100002020022102220010202100220112211111022012202220001210020",
       "010000201210202002002200010022222022012112111222010212120102222221210102112001001022",
       "001000210012110111111202101021212221000101021021021211021221000111100202101200010122",
       "000100221222002222022202010121111210202200012001222011212000211200122202100120211002",
       "000010212010020100110002001011101112122110211102212121200111102212021122100010201120",
       "000001201122121221222212000110100102011101201012001102201222221110211011100001200102")
    M = Matrix(GF(3),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_729_616_523_506():
    r"""
    Return a `(729, 616, 523, 506)`-strongly regular graph.

    This graph is built from a projective ternary `[56,6]` code with weights
    `36, 45`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_729_616_523_506
        sage: G = SRG_729_616_523_506()              # not tested (3s)
        sage: G.is_strongly_regular(parameters=True) # not tested (3s)
        (729, 616, 523, 506)
    """
    x=("10000021022112022210202200202122221120200112100200111102",
       "01000020121020221101202120220001110202220110010222122212",
       "00100021001211211020022112222122002210122100101222020020",
       "00010022122200010012221111121001121211212002110020010101",
       "00001021201002220211121011010222000111021002011201112112",
       "00000120112212111201011001002111121101002212001022222010")
    M = Matrix(GF(3),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_625_364_213_210():
    r"""
    Return a `(625, 364, 213, 210)`-strongly regular graph.

    This graph is built from a projective 5-ary `[88,5]` code with weights `64,
    72`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_625_364_213_210
        sage: G = SRG_625_364_213_210()              # long time
        sage: G.is_strongly_regular(parameters=True) # long time
        (625, 364, 213, 210)
    """
    x=("10004323434444234221223441130101034431234004441141003110400203240",
       "01003023101220331314013121123212111200011403221341101031340421204",
       "00104120244011212302124203142422240001230144213220111213034240310",
       "00012321211123213343321143204040211243210011144140014401003023101")
    M = Matrix(GF(5),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_625_416_279_272():
    r"""
    Return a `(625, 416, 279, 272)`-strongly regular graph.

    This graph is built from a projective 5-ary `[52,4]` code with weights `40,
    45`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_625_416_279_272
        sage: G = SRG_625_416_279_272()               # long time
        sage: G.is_strongly_regular(parameters=True) # long time
        (625, 416, 279, 272)
    """
    x=("1000432343444423422122344123113041011022221414310431",
       "0100302310122033131401312133032331123141114414001300",
       "0010412024401121230212420301411224123332332300210011",
       "0001232121112321334332114324420140440343341412401244")
    M = Matrix(GF(5),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_243_220_199_200():
    r"""
    Return a `(243, 220, 199, 200)`-strongly regular graph.

    This graph is built from a projective ternary `[55,5]` code with weights
    `36, 45`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_243_220_199_200
        sage: G = SRG_243_220_199_200()
        sage: G.is_strongly_regular(parameters=True)
        (243, 220, 199, 200)
    """
    x=("1000010122200120121002211022111101011212112022022020002",
       "0100011101120102100102202121022211112000020211221222002",
       "0010021021222220122011212220021121100021220002100102201",
       "0001012221012012100200102211110211121211201002202000222",
       "0000101222101201210020110221111020112121120120220200022")
    M = Matrix(GF(3),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_729_476_313_306():
    r"""
    Return a `(729, 476, 313, 306)`-strongly regular graph.

    This graph is built from a projective ternary `[126,6]` code with weights
    `81, 90`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_729_476_313_306
        sage: G = SRG_729_476_313_306()               # not tested (5s)
        sage: G.is_strongly_regular(parameters=True) # not tested (5s)
        (729, 476, 313, 306)
    """
    x=("100000210221121211211101220021210000100011020200201101121021122102020111100122122221120200110001010222000021110110011211110210",
       "010000201210202002002111012020001001110012222220221211200120201212222102210100001110202220121001110211200120221121012001221201",
       "001000210012110111111112021220210102211012212122200222212000112220212011021102122002210122122101120210120100102212112112202000",
       "000100221222002222022012122120201012021112211112111120010221100121011012202201001121211212002211120210012201120021222121000110",
       "000010212010020100110202102200200102002122111011112210121010202111121212020010222000111021000222122210001011222102100121210221",
       "000001201122121221222021100221200012000220101001022022100122112010102222002122111121101002200020221110000122202000221222201202")
    M = Matrix(GF(3),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_729_420_243_240():
    r"""
    Return a `(729, 420, 243, 240)`-strongly regular graph.

    This graph is built from a projective ternary `[154,6]` code with weights
    `99, 108`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_729_420_243_240
        sage: G = SRG_729_420_243_240()               # not tested (5s)
        sage: G.is_strongly_regular(parameters=True) # not tested (5s)
        (729, 420, 243, 240)
    """
    x=("10000021022112121121110122002121000010001102020020110112102112202221021020201"+
       "20202212102220222022222110122210022201211222111110211101121002011102101111002",
       "01000020121020200200211101202000100111001222222022121120012020122110122122221"+
       "02222102012112111221111021101101021121002001022221202211100102212212010222102",
       "00100021001211011111111202122021010221101221212220022221200011221102002202120"+
       "20121121000101000111020212200020121210011112210001001022001012222020000100212",
       "00010022122200222202201212212020101202111221111211112001022110001001221210110"+
       "12211020202200222000021101010212001022212020002112011200021100210001100121020",
       "00001021201002010011020210220020010200212211101111221012101020222021111111212"+
       "11120012122110211222201220220201222200102101111020112221020112012102211120101",
       "00000120112212122122202110022120001200022010100102202210012211211120100101022"+
       "01011212011101110111112202111200111021221112222211222020120010222012022220012")
    M = Matrix(GF(3),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_1024_825_668_650():
    r"""
    Return a `(1024, 825, 668, 650)`-strongly regular graph.

    This graph is built from a projective binary `[198,10]` code with weights
    `96, 112`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_1024_825_668_650
        sage: G = SRG_1024_825_668_650()               # not tested (13s)
        sage: G.is_strongly_regular(parameters=True)   # not tested (13s)
        (1024, 825, 668, 650)
    """
    x=("1000000000111111101010000100011001111101101010010010011110001101111001101100111000111101010101110011"+
       "11110010100111101001001100111101011110111100101101110100111100011111011011100111110010100110110000",
       "0100000000010110011100100010101010101111011010001001010110011010101011011101000110000001101101010110"+
       "10110111110101000000011011001100010111110001001011011100111100100000110001011001110110011101011000",
       "0010000000011100111110111011000011010100100011110000001100011011101111001010001100110110000001111000"+
       "11000000101011010111110101000111110010011011101110000010110100000011100010011111100100111101010010",
       "0001000000001111100010000000100101010001110111100010010010010111000100101100010001001110111101110100"+
       "10010101101100110011010011101100110100100011011101100000110011110011111000000010110101011111101111",
       "0000100000110010010000010110000111010011010101000010110100101010011011000011001100001110011011110001"+
       "11101000010000111101101100111100001011010010111011100101101001111000100011000010110111111111011100",
       "0000010000110100111001111011010000101110001011100010010010010111100101011001011011100110101110100001"+
       "01101010110010100011000101111100100001110111001001001001001100001101110110000110101010011010101101",
       "0000001000011011110010110100010010001100000011001000011101000110001101001000110110010101011011001111"+
       "01111111010011111010100110011001110001001000001110000110111011010000011101001110111001011011001011",
       "0000000100111001101011110010111100100001010100100110001100100110010101111001100101101001000101011000"+
       "10001001111101011101001001010111010011011101010011010000101010011001010110011110010000011011111001",
       "0000000010101011010101010101011100111101111110100011011001001010111101100111010110100101100110101100"+
       "00000001100011110110010101100001000000010100001101111011111000110001100101101010000001110101011100",
       "0000000001101100111101011000010000000011010100000110101010011010100111100001000011010011011101110111"+
       "01110111011110101100100100110110011100001001000001010011010010010111110011101011101001101101011010")
    M = Matrix(GF(2),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_512_73_12_10():
    r"""
    Return a `(512, 73, 12, 10)`-strongly regular graph.

    This graph is built from a projective binary `[219,9]` code with weights
    `96, 112`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_512_73_12_10
        sage: G = SRG_512_73_12_10()                 # not tested (3s)
        sage: G.is_strongly_regular(parameters=True) # not tested (3s)
        (512, 73, 12, 10)
    """
    x=("10100011001110100010100010010000100100001011110010001010011000000001101011110011001001000010011110111011111"+
       "0001001010110110110111001100111100011011101000000110101110001010100011110011111111110111010100101011000101111111",
       "01100010110101110100001000010110001010010010011000111101111001011101000011101011100111111110001100000111010"+
       "1101001000110001111011001100101011110101011110010001101011110000100000101101100010110100001111001100110011001111",
       "00010010001001011011001110011101111110000000101110101000110110011001110101011011101011011011000010010011111"+
       "1110110100111111000000110011101101000000001010000000011000111111100101100001110011110001110011110110100111100001",
       "00001000100010101110101110011100010101110011010110000001111111100111010000101110001010100100000001011010111"+
       "1001001000000011000011001100100100111010000000001010111001001100100101011110001100110001000000111001100100100111",
       "00000101010100010101101110011101001000101110000000000111101100011000000001110100000001011010101001111110110"+
       "0010110111100111000000110011110110101101110000001111100001010001100101100001110011110001101101000000000000100001",
       "00000000000000000000010000011101011100100010000110110100101011001011001100000001011000101010100111000111101"+
       "0011100011011011011111100010011100010111101001011001001101100010011010001011010110001110100001001111110010100100",
       "00000000000000000000000001011010110110101111010110101001001001000101010000000000001011000011000010100100110"+
       "0000110000111101100010000111111111101101001010110000111111101110101011010010010001011101110011111001100100101110",
       "00000000000000000000000000110111101011110010101110000110010010100010001010000000010100011000101000010011000"+
       "0110000111100110100001001011111111111010110000001010111111110011110110001100100010101011101101110110011000110110",
       "00000000000000000000000000000000000000000000000001111111111111111111111110000001111111111111111111111111111"+
       "1111111100000000000011111111111111000000111111111111111111000000000000111111111111000000000000000000111111000110")
    M = Matrix(GF(2),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_512_219_106_84():
    r"""
    Return a `(512, 219, 106, 84)`-strongly regular graph.

    This graph is built from a projective binary `[73,9]` code with weights `32,
    40`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_512_219_106_84
        sage: G = SRG_512_219_106_84()
        sage: G.is_strongly_regular(parameters=True)
        (512, 219, 106, 84)
    """
    x=("1010010100000010100000101010001100110101101101000010110010100100111011101",
       "0110000110000101101111001101000100111111101011011101110010110001100111100",
       "0001010000000001111111011010100101001111011010101100001010000001110100001",
       "0000100100000001111111100111000011110011110101000001010110000001011010001",
       "0000001010000001111110111100011000111100101110010010101100000001101001001",
       "0000000001000111001010110010011001101001011010110110011001010111100010010",
       "0000000000100100011000100100111100001100101111010001011011111000110011110",
       "0000000000010111001100101011111110101010000000000100111110000001111111100",
       "0000000000001011100001000011011010110001110101101100001100101110101110110")
    M = Matrix(GF(2),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_256_153_92_90():
    r"""
    Return a `(256, 153, 92, 90)`-strongly regular graph.

    This graph is built from a projective 4-ary `[34,4]` code with weights `24,
    28`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_256_153_92_90
        sage: G = SRG_256_153_92_90()
        sage: G.is_strongly_regular(parameters=True)
        (256, 153, 92, 90)
    """
    K = GF(4,conway=True, prefix='x')
    F = K.gens()[0]
    J = F*F
    x = [[1,0,0,0,1,F,F,J,1,0,F,F,0,1,J,F,F,J,J,J,F,F,J,J,J,1,J,F,1,0,1,F,J,1],
         [0,1,0,0,F,F,1,J,1,1,J,1,F,F,0,0,1,0,F,F,0,1,J,F,F,1,0,0,0,1,F,F,J,1],
         [0,0,1,0,1,0,0,F,F,1,J,1,1,J,1,F,F,F,J,1,0,F,F,0,1,J,F,F,1,0,0,0,1,F],
         [0,0,0,1,F,F,J,1,0,F,F,0,1,J,F,F,1,J,J,F,F,J,J,J,1,J,F,1,0,1,F,J,1,J]]
    M = Matrix(K,[map(K,l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

def SRG_256_170_114_110():
    r"""
    Return a `(256, 170, 114, 110)`-strongly regular graph.

    This graph is built from a projective binary `[85,8]` code with weights `40,
    48`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_256_170_114_110
        sage: G = SRG_256_170_114_110()
        sage: G.is_strongly_regular(parameters=True)
        (256, 170, 114, 110)
    """
    x=("1000000010011101010001000011100111000111111010110001101101000110010011001101011100001",
       "0100000011010011111001100010010100100100000111101001011011100101011010101011110010001",
       "0010000011110100101101110010101101010101111001000101000000110100111110011000100101001",
       "0001000011100111000111111010110001101101000110010011001101011100001100000001001110101",
       "0000100011101110110010111110111111110001011001111000001011101000010101001101111011011",
       "0000010011101010001000011100111000111111010110001101101000110010011001101011100001100",
       "0000001001110101000100001110011100011111101011000110110100011001001100110101110000110",
       "0000000100111010100010000111001110001111110101100011011010001100100110011010111000011")
    M = Matrix(GF(2),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

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
    return JohnsonGraph(10,3).distance_graph([2])

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
    return JohnsonGraph(9,4).distance_graph([1,4])

def SRG_175_72_20_36():
    r"""
    Return a `(175,72,20,36)`-strongly regular graph

    This graph is obtained from the line graph of
    :meth:`~sage.graphs.graph_generators.GraphGenerators.HoffmanSingletonGraph`. Setting
    two vertices to be adjacent if their distance in the line graph is exactly
    two yields the strongly regular graph. For more information, see
    http://www.win.tue.nl/~aeb/graphs/McL.html.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_175_72_20_36
        sage: G = SRG_175_72_20_36()
        sage: G.is_strongly_regular(parameters=True)
        (175, 72, 20, 36)
    """
    return HoffmanSingletonGraph().line_graph().distance_graph([2])

def SRG_126_50_13_24():
    r"""
    Return a `(126,50,13,24)`-strongly regular graph

    This graph is a subgraph of
    :meth:`~sage.graphs.strongly_regular_db.SRG_175_72_20_36`.
    This construction, due to Goethals, is given in §10B.(vii) of [BvL84]_.

    EXAMPLES::

        sage: from sage.graphs.strongly_regular_db import SRG_126_50_13_24
        sage: G = SRG_126_50_13_24()
        sage: G.is_strongly_regular(parameters=True)
        (126, 50, 13, 24)
    """
    from sage.graphs.generators.smallgraphs import HoffmanSingletonGraph
    from sage.graphs.strongly_regular_db import SRG_175_72_20_36
    hs = HoffmanSingletonGraph()
    s = set(hs.vertices()).difference(hs.neighbors(0)+[0])
    return SRG_175_72_20_36().subgraph(hs.edge_boundary(s,s))

def SRG_81_50_31_30():
    r"""
    Return a `(81, 50, 31, 30)`-strongly regular graph.

    This graph is built from a projective ternary `[15,4]` code with weights `9,
    12`, obtained from Eric Chen's `database of two-weight codes
    <http://moodle.tec.hkr.se/~chen/research/2-weight-codes/search.php>`__.

    .. SEEALSO::

        :func:`strongly_regular_from_two_weight_code` -- build a strongly regular graph from
        a two-weight code.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import SRG_81_50_31_30
        sage: G = SRG_81_50_31_30()
        sage: G.is_strongly_regular(parameters=True)
        (81, 50, 31, 30)
    """
    x=("100022021001111",
       "010011211122000",
       "001021112100011",
       "000110120222220")
    M = Matrix(GF(3),[list(l) for l in x])
    return strongly_regular_from_two_weight_code(LinearCode(M))

cdef bint seems_feasible(int v, int k, int l, int mu):
    r"""
    Tests is the set of parameters seems feasible

    INPUT:

    - ``v,k,l,mu`` (integers)
    """
    cdef int r,s,f,g
    cdef uint_fast32_t tmp[2]

    if (v<0 or k<=0 or l<0 or mu<0 or
        k>=v-1 or l>=k or mu>=k or
        v-2*k+mu-2 < 0 or # lambda of complement graph >=0
        v-2*k+l    < 0 or # μ of complement graph >=0
        mu*(v-k-1) != k*(k-l-1)):
        return False

    # Conference graphs. Only possible if 'v' is a sum of two squares (3.A of
    # [BvL84]
    if (v-1)*(mu-l)-2*k == 0:
        return two_squares_c(v,tmp)

    rr,ss = eigenvalues(v,k,l,mu)
    if rr is None:
        return False
    r,s = rr,ss

    # 1.3.1 of [Distance-regular graphs]
    # "Integrality condition"
    if ((s+1)*(k-s)*k) % (mu*(s-r)):
        return False

    # Theorem 21.3 of [WilsonACourse] or
    # 3.B of [BvL84]
    # (Krein conditions)
    if ((r+1)*(k+r+2*r*s) > (k+r)*(s+1)**2 or
        (s+1)*(k+s+2*r*s) > (k+s)*(r+1)**2):
        return False

    # multiplicity of eigenvalues 'r,s' (f=lambda_r, g=lambda_s)
    #
    # They are integers (checked by the 'integrality condition').
    f = -k*(s+1)*(k-s)/((k+r*s)*(r-s))
    g =  k*(r+1)*(k-r)/((k+r*s)*(r-s))

    # 3.C of [BvL84]
    # (Absolute bound)
    if (2*v > f*(f+3) or
        2*v > g*(g+3)):
        return False

    # 3.D of [BvL84]
    # (Claw bound)
    if (mu != s**2    and
        mu != s*(s+1) and
        2*(r+1) > s*(s+1)*(mu+1)):
        return False

    # 3.E of [BvL84]
    # (the Case μ=1)
    if mu == 1:
        if (   k  % (l+1) or
            (v*k) % ((l+1)*(l+2))):
            return False

    # 3.F of [BvL84]
    # (the Case μ=2)
    if mu == 2 and 2*k < l*(l+3) and k%(l+1):
        return False

    return True

def strongly_regular_graph(int v,int k,int l,int mu=-1,bint existence=False,bint check=True):
    r"""
    Return a `(v,k,\lambda,\mu)`-strongly regular graph.

    This function relies partly on Andries Brouwer's `database of strongly
    regular graphs <http://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__. See
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

        sage: graphs.strongly_regular_graph(324,57,0,12,existence=True)
        False
        sage: graphs.strongly_regular_graph(324,57,0,12)
        Traceback (most recent call last):
        ...
        EmptySetError: Andries Brouwer's database reports that no (324, 57, 0,
        12)-strongly regular graph exists. Comments: <a
        href="srgtabrefs.html#GavrilyukMakhnev05">Gavrilyuk & Makhnev</a> and <a
        href="srgtabrefs.html#KaskiOstergard07">Kaski & stergrd</a>

    A set of parameters unknown to be realizable in Andries Brouwer's database::

        sage: graphs.strongly_regular_graph(324,95,22,30,existence=True)
        Unknown
        sage: graphs.strongly_regular_graph(324,95,22,30)
        Traceback (most recent call last):
        ...
        RuntimeError: Andries Brouwer's database reports that no
        (324,95,22,30)-strongly regular graph is known to exist.
        Comments:

    A realizable set of parameters that Sage cannot realize (help us!)::

        sage: graphs.strongly_regular_graph(1288, 495, 206, existence=True)
        True
        sage: graphs.strongly_regular_graph(1288, 495, 206)
        Traceback (most recent call last):
        ...
        RuntimeError: Andries Brouwer's database claims that such a (1288,495,206,180)-strongly
        regular graph exists, but Sage does not know how to build it.
        ...

    A large unknown set of parameters (not in Andries Brouwer's database)::

        sage: graphs.strongly_regular_graph(1394,175,0,25,existence=True)
        Unknown
        sage: graphs.strongly_regular_graph(1394,175,0,25)
        Traceback (most recent call last):
        ...
        RuntimeError: Sage cannot figure out if a (1394,175,0,25)-strongly regular graph exists.

    Test the Claw bound (see 3.D of [BvL84]_)::

        sage: graphs.strongly_regular_graph(2058,242,91,20,existence=True)
        False
    """
    load_brouwer_database()
    if mu == -1:
        mu = k*(k-l-1)//(v-k-1)

    params = (v,k,l,mu)
    params_complement = (v,v-k-1,v-2*k+mu-2,v-2*k+l)

    if not seems_feasible(v,k,l,mu):
        if existence:
            return False
        raise ValueError("There exists no "+str(params)+"-strongly regular graph")

    def check_srg(G):
        if check and (v,k,l,mu) != G.is_strongly_regular(parameters=True):
            raise RuntimeError("Sage built an incorrect {}-SRG.".format((v,k,l,mu)))
        return G

    constructions = {
        ( 27,  16, 10,  8): [SchlaefliGraph],
        ( 36,  14,  4,  6): [Graph,('c~rLDEOcKTPO`U`HOIj@MWFLQFAaRIT`HIWqPsQQJ'+
          'DXGLqYM@gRLAWLdkEW@RQYQIErcgesClhKefC_ygSGkZ`OyHETdK[?lWStCapVgKK')],
        ( 40,  12,  2,  4): [Graph,('g}iS[A@_S@OA_BWQIGaPCQE@CcQGcQECXAgaOdS@a'+
          'CWEEAOIBH_HW?scb?f@GMBGGhIPGaQoh?q_bD_pGPq_WI`T_DBU?R_dECsSARGgogBO'+
          '{_IPBKZ?DI@Wgt_E?MPo{_?')],
        ( 45,  12,  3,  3): [Graph,('l~}CKMF_C?oB_FPCGaICQOaH@DQAHQ@Ch?aJHAQ@G'+
          'P_CQAIGcAJGO`IcGOY`@IGaGHGaKSCDI?gGDgGcE_@OQAg@PCSO_hOa`GIDADAD@XCI'+
          'ASDKB?oKOo@_SHCc?SGcGd@A`B?bOOHGQH?ROQOW`?XOPa@C_hcGo`CGJK')],
        ( 50,   7,  0,  1): [HoffmanSingletonGraph],
        ( 56,  10,  0,  2): [SimsGewirtzGraph],
        ( 64,  18,  2,  6): [Graph,('~?@?~aK[A@_[?O@_B_?O?K?B_?A??K??YQQPHGcQQ'+
          'CaPIOHAX?POhAPIC`GcgSAHDE?PCiC@BCcDADIG_QCocS@AST?OOceGG@QGcKcdCbCB'+
          'gIEHAScIDDOy?DAWaEg@IQO?maHPOhAW_dBCX?s@HOpKD@@GpOpHO?bCbHGOaGgpWQQ'+
          '?PDDDw@A_CSRIS_P?GeGpg`@?EOcaJGccbDC_dLAc_pHOe@`ocEGgo@sRo?WRAbAcPc'+
          '?iCiHEKBO_hOiOWpOSGSTBQCUAW_DDIWOqHBO?gghw_?`kOAXH?\\Ds@@@CpIDKOpc@'+
          'OCoeIS_YOgGATGaqAhKGA?cqDOwQKGc?')],
        ( 77,  16,   0,  4): [M22Graph],
        ( 81,  50,  31, 30): [SRG_81_50_31_30],
        (100,  22,   0,  6): [HigmanSimsGraph],
        (100,  44,  18, 20): [SRG_100_44_18_20],
        (100,  45,  20, 20): [SRG_100_45_20_20],
        (120,  63,  30, 36): [SRG_120_63_30_36],
        (120,  77,  52, 44): [SRG_120_77_52_44],
        (126,  25,   8,  4): [SRG_126_25_8_4],
        (126,  50,  13, 24): [SRG_126_50_13_24],
        (162,  56,  10, 24): [LocalMcLaughlinGraph],
        (175,  72,  20, 36): [SRG_175_72_20_36],
        (176, 105,  68, 54): [SRG_176_105_68_54],
        (196,  91,  42, 42): [SRG_196_91_42_42],
        (231,  30,   9,  3): [CameronGraph],
        (243, 220, 199,200): [SRG_243_220_199_200],
        (253, 140,  87, 65): [SRG_253_140_87_65],
        (256, 170, 114,110): [SRG_256_170_114_110],
        (256, 187, 138,132): [SRG_256_187_138_132],
        (256, 153,  92, 90): [SRG_256_153_92_90],
        (275, 112,  30, 56): [McLaughlinGraph],
        (280, 135,  70, 60): [SRG_280_135_70_60],
        (512, 219, 106, 84): [SRG_512_219_106_84],
        (512,  73,  12, 10): [SRG_512_73_12_10],
        (625, 416, 279,272): [SRG_625_416_279_272],
        (625, 364, 213,210): [SRG_625_364_213_210],
        (729, 616, 523,506): [SRG_729_616_523_506],
        (729, 420, 243,240): [SRG_729_420_243_240],
        (729, 560, 433,420): [SRG_729_560_433_420],
        (729, 476, 313,306): [SRG_729_476_313_306],
        (729, 532, 391,380): [SRG_729_532_391_380],
        (1024,825, 668,650): [SRG_1024_825_668_650],
    }

    if params in constructions:
        val = constructions[params]
        return True if existence else check_srg(val[0](*val[1:]))
    if params_complement in constructions:
        val = constructions[params_complement]
        return True if existence else check_srg(val[0](*val[1:]).complement())

    test_functions = [is_paley, is_johnson,
                      is_orthogonal_array_block_graph,
                      is_steiner, is_affine_polar,
                      is_orthogonal_polar,
                      is_unitary_polar,
                      is_unitary_dual_polar,
                      is_RSHCD,
                      is_twograph_descendant_of_srg,
                      is_taylor_twograph_srg]

    # Going through all test functions, for the set of parameters and its
    # complement.
    for f in test_functions:
        if f(*params):
            if existence:
                return True
            ans = f(*params)
            return check_srg(ans[0](*ans[1:]))
        if f(*params_complement):
            if existence:
                return True
            ans = f(*params_complement)
            return check_srg(ans[0](*ans[1:]).complement())

    # From now on, we have no idea how to build the graph.
    #
    # We try to return the most appropriate error message.

    global _brouwer_database
    brouwer_data = _brouwer_database.get(params,None)

    if brouwer_data is not None:
        if brouwer_data['status'] == 'impossible':
            if existence:
                return False
            raise EmptySetError("Andries Brouwer's database reports that no "+
                                str((v,k,l,mu))+"-strongly regular graph exists. "+
                                "Comments: "+brouwer_data['comments'].encode('ascii','ignore'))

        if brouwer_data['status'] == 'open':
            if existence:
                return Unknown
            raise RuntimeError(("Andries Brouwer's database reports that no "+
                                "({},{},{},{})-strongly regular graph is known "+
                                "to exist.\nComments: ").format(v,k,l,mu)
                               +brouwer_data['comments'].encode('ascii','ignore'))

        if brouwer_data['status'] == 'exists':
            if existence:
                return True
            raise RuntimeError(("Andries Brouwer's database claims that such a "+
                                "({},{},{},{})-strongly regular graph exists, but "+
                                "Sage does not know how to build it. If *you* do, "+
                                "please get in touch with us on sage-devel!\n"+
                                "Comments: ").format(v,k,l,mu)
                               +brouwer_data['comments'].encode('ascii','ignore'))
    if existence:
        return Unknown
    raise RuntimeError(("Sage cannot figure out if a ({},{},{},{})-strongly "+
                        "regular graph exists.").format(v,k,l,mu))

def apparently_feasible_parameters(int n):
    r"""
    Return a list of parameters `(v,k,\lambda,\mu)` which are a-priori feasible.

    Note that some of those that it returns may also be infeasible for more
    involved reasons.

    INPUT:

    - ``n`` (integer) -- return all a-priori feasible tuples `(v,k,\lambda,\mu)`
      for `v<n`

    EXAMPLE:

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
        sage: all(graphs.strongly_regular_graph(*x,existence=True) for x in small_feasible)
        True

    But that becomes wrong for `v<60` (because of the non-existence of a
    `(49,16,3,6)`-strongly regular graph)::

        sage: small_feasible = apparently_feasible_parameters(60)
        sage: all(graphs.strongly_regular_graph(*x,existence=True) for x in small_feasible)
        False

    """
    cdef int v,k,l,mu
    feasible = set()
    for v in range(n):
        for k in range(1,v-1):
            for l in range(k-1):
                mu = k*(k-l-1)//(v-k-1)
                if seems_feasible(v,k,l,mu):
                    feasible.add((v,k,l,mu))
    return feasible

cdef load_brouwer_database():
    r"""
    Loads Andries Brouwer's database into _brouwer_database.
    """
    global _brouwer_database
    if _brouwer_database is not None:
        return
    import json

    from sage.env import SAGE_SHARE
    with open(SAGE_SHARE+"/graphs/brouwer_srg_database.json",'r') as datafile:
        _brouwer_database = {(v,k,l,mu):{'status':status,'comments':comments}
                             for (v,k,l,mu,status,comments) in json.load(datafile)}

def _check_database():
    r"""
    Checks the coherence of Andries Brouwer's database with Sage.

    The function also outputs some statistics on the database.

    EXAMPLE::

        sage: from sage.graphs.strongly_regular_db import _check_database
        sage: _check_database() # long time
        Sage cannot build a (45   22   10   11  ) that exists. Comment from Brouwer's database: <a href="srgtabrefs.html#Mathon78">Mathon</a>; 2-graph*
        ...
        In Andries Brouwer's database:
        - 448 impossible entries
        - 2950 undecided entries
        - 1140 realizable entries (Sage misses ... of them)

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
    for params,dic in sorted(saved_database.items()):
        sage_answer = strongly_regular_graph(*params,existence=True)
        if dic['status'] == 'open':
            if sage_answer:
                print "Sage can build a {}, Brouwer's database cannot".format(params)
            assert sage_answer is not False
        elif dic['status'] == 'exists':
            if sage_answer is not True:
                print (("Sage cannot build a ({:<4} {:<4} {:<4} {:<4}) that exists. "+
                       "Comment from Brouwer's database: ").format(*params)
                       +dic['comments'].encode('ascii','ignore'))
                missed += 1
            assert sage_answer is not False
        elif dic['status'] == 'impossible':
            assert sage_answer is not True
        else:
            assert False # must not happen

    status = [x['status'] for x in saved_database.values()]
    print "\nIn Andries Brouwer's database:"
    print "- {} impossible entries".format(status.count('impossible'))
    print "- {} undecided entries".format(status.count('open'))
    print "- {} realizable entries (Sage misses {} of them)".format(status.count('exists'),missed)

    # Reassign its value to the global database
    _brouwer_database = saved_database
