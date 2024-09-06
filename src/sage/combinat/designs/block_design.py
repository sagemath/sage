# sage.doctest: needs sage.rings.finite_rings
r"""
Block designs

A *block design* is a set together with a family of subsets (repeated subsets
are allowed) whose members are chosen to satisfy some set of properties that are
deemed useful for a particular application. See :wikipedia:`Block_design`.

REFERENCES:

- Block design from wikipedia: :wikipedia:`Block_design`

- What is a block design?,
  http://designtheory.org/library/extrep/extrep-1.1-html/node4.html (in 'The
  External Representation of Block Designs' by Peter J. Cameron, Peter
  Dobcsanyi, John P. Morgan, Leonard H. Soicher)

.. [Hu57] Daniel R. Hughes, "A class of non-Desarguesian projective planes",
   The Canadian Journal of Mathematics (1957), http://cms.math.ca/cjm/v9/p378

.. [We07] Charles Weibel, "Survey of Non-Desarguesian planes" (2007), notices of
   the AMS, vol. 54 num. 10, pages 1294--1303

AUTHORS:

- Quentin Honoré (2015): construction of Hughes plane :issue:`18527`

- Vincent Delecroix (2014): rewrite the part on projective planes :issue:`16281`

- Peter Dobcsanyi and David Joyner (2007-2008)

  This is a significantly modified form of the module block_design.py (version
  0.6) written by Peter Dobcsanyi peter@designtheory.org. Thanks go to Robert
  Miller for lots of good design suggestions.

.. TODO::

    Implement more finite non-Desarguesian plane as in [We07]_ and
    :wikipedia:`Non-Desarguesian_plane`.

Functions and methods
---------------------
"""

#*****************************************************************************
#       Copyright (C) 2007 Peter Dobcsanyi <peter@designtheory.org>
#       Copyright (C) 2007 David Joyner <wdjoyner@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************
from sage.arith.misc import binomial, integer_floor, is_prime_power
from sage.categories.sets_cat import EmptySetError
from sage.misc.lazy_import import lazy_import
from sage.misc.unknown import Unknown
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ

from .incidence_structures import IncidenceStructure

lazy_import('sage.libs.gap.libgap', 'libgap')
lazy_import('sage.matrix.matrix_space', 'MatrixSpace')
lazy_import('sage.modules.free_module', 'VectorSpace')
lazy_import('sage.rings.finite_rings.finite_field_constructor', 'FiniteField')


BlockDesign = IncidenceStructure

###  utility functions  -------------------------------------------------------


def tdesign_params(t, v, k, L):
    """
    Return the design's parameters: `(t, v, b, r , k, L)`. Note that `t` must be
    given.

    EXAMPLES::

        sage: BD = BlockDesign(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
        sage: from sage.combinat.designs.block_design import tdesign_params
        sage: tdesign_params(2,7,3,1)
        (2, 7, 7, 3, 3, 1)
    """
    x = binomial(v, t)
    y = binomial(k, t)
    b = divmod(L * x, y)[0]
    x = binomial(v-1, t-1)
    y = binomial(k-1, t-1)
    r = integer_floor(L * x/y)
    return (t, v, b, r, k, L)


def are_hyperplanes_in_projective_geometry_parameters(v, k, lmbda, return_parameters=False):
    r"""
    Return ``True`` if the parameters ``(v,k,lmbda)`` are the one of hyperplanes in
    a (finite Desarguesian) projective space.

    In other words, test whether there exists a prime power ``q`` and an integer
    ``d`` greater than two such that:

    - `v = (q^{d+1}-1)/(q-1) = q^d + q^{d-1} + ... + 1`
    - `k = (q^d - 1)/(q-1) = q^{d-1} + q^{d-2} + ... + 1`
    - `lmbda = (q^{d-1}-1)/(q-1) = q^{d-2} + q^{d-3} + ... + 1`

    If it exists, such a pair ``(q,d)`` is unique.

    INPUT:

    - ``v``, ``k``, ``lmbda`` -- integers

    OUTPUT:

    - a boolean or, if ``return_parameters`` is set to ``True`` a pair
      ``(True, (q,d))`` or ``(False, (None,None))``.

    EXAMPLES::

        sage: from sage.combinat.designs.block_design import are_hyperplanes_in_projective_geometry_parameters
        sage: are_hyperplanes_in_projective_geometry_parameters(40, 13, 4)
        True
        sage: are_hyperplanes_in_projective_geometry_parameters(40, 13, 4,
        ....:                                                   return_parameters=True)
        (True, (3, 3))
        sage: PG = designs.ProjectiveGeometryDesign(3, 2, GF(3))                        # needs sage.combinat
        sage: PG.is_t_design(return_parameters=True)                                    # needs sage.combinat
        (True, (2, 40, 13, 4))

        sage: are_hyperplanes_in_projective_geometry_parameters(15, 3, 1)
        False
        sage: are_hyperplanes_in_projective_geometry_parameters(15, 3, 1,
        ....:                                                   return_parameters=True)
        (False, (None, None))

    TESTS::

        sage: sgp = lambda q,d: ((q**(d+1)-1)//(q-1), (q**d-1)//(q-1), (q**(d-1)-1)//(q-1))
        sage: for q in [3,4,5,7,8,9,11]:
        ....:     for d in [2,3,4,5]:
        ....:         v,k,l = sgp(q,d)
        ....:         assert are_hyperplanes_in_projective_geometry_parameters(v,k,l,True) == (True, (q,d))
        ....:         assert are_hyperplanes_in_projective_geometry_parameters(v+1,k,l) is False
        ....:         assert are_hyperplanes_in_projective_geometry_parameters(v-1,k,l) is False
        ....:         assert are_hyperplanes_in_projective_geometry_parameters(v,k+1,l) is False
        ....:         assert are_hyperplanes_in_projective_geometry_parameters(v,k-1,l) is False
        ....:         assert are_hyperplanes_in_projective_geometry_parameters(v,k,l+1) is False
        ....:         assert are_hyperplanes_in_projective_geometry_parameters(v,k,l-1) is False
    """
    import sage.arith.all as arith

    q1 = Integer(v - k)
    q2 = Integer(k - lmbda)

    if (lmbda <= 0 or q1 < 4 or q2 < 2 or
        not q1.is_prime_power() or
        not q2.is_prime_power()):
        return (False,(None,None)) if return_parameters else False

    p1,e1 = q1.factor()[0]
    p2,e2 = q2.factor()[0]

    k = arith.gcd(e1,e2)
    d = e1//k
    q = p1**k
    if e2//k != d-1 or lmbda != (q**(d-1)-1)//(q-1):
        return (False,(None,None)) if return_parameters else False

    return (True, (q,d)) if return_parameters else True


def ProjectiveGeometryDesign(n, d, F, algorithm=None, point_coordinates=True, check=True):
    r"""
    Return a projective geometry design.

    The projective geometry design `PG_d(n,q)` has for points the lines of
    `\GF{q}^{n+1}`, and for blocks the `(d+1)`-dimensional subspaces of
    `\GF{q}^{n+1}`, each of which contains `\frac {|\GF{q}|^{d+1}-1} {|\GF{q}|-1}` lines.
    It is a `2`-design with parameters

    .. MATH::

        v = \binom{n+1}{1}_q,\ k = \binom{d+1}{1}_q,\ \lambda =
        \binom{n-1}{d-1}_q

    where the `q`-binomial coefficient `\binom{m}{r}_q` is defined by

    .. MATH::

        \binom{m}{r}_q = \frac{(q^m - 1)(q^{m-1} - 1) \cdots (q^{m-r+1}-1)}
              {(q^r-1)(q^{r-1}-1)\cdots (q-1)}

    .. SEEALSO::

        :func:`AffineGeometryDesign`

    INPUT:

    - ``n`` -- the projective dimension

    - ``d`` -- the dimension of the subspaces which make up the blocks

    - ``F`` -- a finite field or a prime power

    - ``algorithm`` -- set to ``None`` by default, which results in using Sage's
      own implementation. In order to use GAP's implementation instead (i.e. its
      ``PGPointFlatBlockDesign`` function) set ``algorithm="gap"``. Note that
      GAP's "design" package must be available in this case, and that it can be
      installed with the ``gap_packages`` spkg.

    - ``point_coordinates`` -- ``True`` by default. Ignored and assumed to be ``False`` if
      ``algorithm="gap"``. If ``True``, the ground set is indexed by coordinates
      in `\GF{q}^{n+1}`.  Otherwise the ground set is indexed by integers.

    - ``check`` -- boolean (default: ``True``); whether to check the output

    EXAMPLES:

    The set of `d`-dimensional subspaces in a `n`-dimensional projective space
    forms `2`-designs (or balanced incomplete block designs)::

        sage: PG = designs.ProjectiveGeometryDesign(4, 2, GF(2)); PG                    # needs sage.combinat
        Incidence structure with 31 points and 155 blocks
        sage: PG.is_t_design(return_parameters=True)                                    # needs sage.combinat
        (True, (2, 31, 7, 7))

        sage: PG = designs.ProjectiveGeometryDesign(3, 1, GF(4))                        # needs sage.combinat
        sage: PG.is_t_design(return_parameters=True)                                    # needs sage.combinat
        (True, (2, 85, 5, 1))

    Check with ``F`` being a prime power::

        sage: PG = designs.ProjectiveGeometryDesign(3, 2, 4); PG
        Incidence structure with 85 points and 85 blocks

    Use coordinates::

        sage: PG = designs.ProjectiveGeometryDesign(2, 1, GF(3))
        sage: PG.blocks()[0]
        [(1, 0, 0), (1, 0, 1), (1, 0, 2), (0, 0, 1)]

    Use indexing by integers::

        sage: PG = designs.ProjectiveGeometryDesign(2, 1, GF(3), point_coordinates=0)
        sage: PG.blocks()[0]
        [0, 1, 2, 12]

    Check that the constructor using gap also works::

        sage: BD = designs.ProjectiveGeometryDesign(2, 1, GF(2), algorithm='gap')  # optional - gap_package_design
        sage: BD.is_t_design(return_parameters=True)                               # optional - gap_package_design
        (True, (2, 7, 3, 1))
    """
    try:
        q = int(F)
    except TypeError:
        q = F.cardinality()
    else:
        from sage.rings.finite_rings.finite_field_constructor import GF
        F = GF(q)

    if algorithm is None:
        from sage.matrix.echelon_matrix import reduced_echelon_matrix_iterator

        points = {p:i for i,p in enumerate(reduced_echelon_matrix_iterator(F,1,n+1,copy=True,set_immutable=True))}
        blocks = []
        for m1 in reduced_echelon_matrix_iterator(F,d+1,n+1,copy=False):
            b = []
            for m2 in reduced_echelon_matrix_iterator(F,1,d+1,copy=False):
                m = m2*m1
                m.echelonize()
                m.set_immutable()
                b.append(points[m])
            blocks.append(b)
        B = BlockDesign(len(points), blocks, name='ProjectiveGeometryDesign', check=check)
        if point_coordinates:
            B.relabel({i:p[0] for p,i in points.items()})

    elif algorithm == "gap":   # Requires GAP's Design
        libgap.load_package("design")
        D = libgap.PGPointFlatBlockDesign(n, F.order(), d)
        v = D['v'].sage()
        gblcks = D['blocks'].sage()
        gB = []
        for b in gblcks:
            gB.append([x - 1 for x in b])
        B = BlockDesign(v, gB, name='ProjectiveGeometryDesign', check=check)

    if check:
        from sage.combinat.q_analogues import q_binomial
        q = F.cardinality()
        if not B.is_t_design(t=2, v=q_binomial(n+1,1,q),
                                  k=q_binomial(d+1,1,q),
                                  l=q_binomial(n-1, d-1, q)):
            raise RuntimeError("error in ProjectiveGeometryDesign "
                    "construction. Please e-mail sage-devel@googlegroups.com")
    return B


def DesarguesianProjectivePlaneDesign(n, point_coordinates=True, check=True):
    r"""
    Return the Desarguesian projective plane of order ``n`` as a 2-design.

    The Desarguesian projective plane of order `n` can also be defined as the
    projective plane over a field of order `n`. For more information, have a
    look at :wikipedia:`Projective_plane`.

    INPUT:

    - ``n`` -- integer which must be a power of a prime number

    - ``point_coordinates`` -- boolean (default: ``True``); whether to label the
      points with their homogeneous coordinates (default) or with integers

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.

    .. SEEALSO::

        :func:`ProjectiveGeometryDesign`

    EXAMPLES::

        sage: designs.DesarguesianProjectivePlaneDesign(2)
        (7,3,1)-Balanced Incomplete Block Design
        sage: designs.DesarguesianProjectivePlaneDesign(3)
        (13,4,1)-Balanced Incomplete Block Design
        sage: designs.DesarguesianProjectivePlaneDesign(4)
        (21,5,1)-Balanced Incomplete Block Design
        sage: designs.DesarguesianProjectivePlaneDesign(5)
        (31,6,1)-Balanced Incomplete Block Design
        sage: designs.DesarguesianProjectivePlaneDesign(6)
        Traceback (most recent call last):
        ...
        ValueError: the order of a finite field must be a prime power
    """
    K = FiniteField(n, 'a')
    n2 = n**2
    relabel = {x:i for i,x in enumerate(K)}
    Kiter = relabel  # it is much faster to iterate through a dict than through
                     # the finite field K

    # we decompose the (equivalence class) of points [x:y:z] of the projective
    # plane into an affine plane, an affine line and a point. At the same time,
    # we relabel the points with the integers from 0 to n^2 + n as follows:
    # - the affine plane is the set of points [x:y:1] (i.e. the third coordinate
    #   is nonzero) and gets relabeled from 0 to n^2-1
    affine_plane = lambda x,y: relabel[x] + n * relabel[y]

    # - the affine line is the set of points [x:1:0] (i.e. the third coordinate is
    #   zero but not the second one) and gets relabeled from n^2 to n^2 + n - 1
    line_infinity = lambda x: n2 + relabel[x]

    # - the point is [1:0:0] and gets relabeled n^2 + n
    point_infinity = n2 + n

    blcks = []

    # the n^2 lines of the form "x = sy + az"
    for s in Kiter:
        for a in Kiter:
            # points in the affine plane
            blcks.append([affine_plane(s*y+a, y) for y in Kiter])
            # point at infinity
            blcks[-1].append(line_infinity(s))

    # the n horizontals of the form "y = az"
    for a in Kiter:
        # points in the affine plane
        blcks.append([affine_plane(x,a) for x in Kiter])
        # point at infinity
        blcks[-1].append(point_infinity)

    # the line at infinity "z = 0"
    blcks.append(range(n2,n2+n+1))
    if check:
        from .designs_pyx import is_projective_plane
        if not is_projective_plane(blcks):
            raise RuntimeError('There is a problem in the function DesarguesianProjectivePlane')
    from .bibd import BalancedIncompleteBlockDesign
    B = BalancedIncompleteBlockDesign(n2+n+1, blcks, check=check)

    if point_coordinates:
        zero = K.zero()
        one = K.one()
        d = {affine_plane(x,y): (x,y,one)
             for x in Kiter
             for y in Kiter}
        d.update({line_infinity(x): (x,one,zero)
                  for x in Kiter})
        d[n2+n] = (one,zero,zero)
        B.relabel(d)

    return B


def q3_minus_one_matrix(K):
    r"""
    Return a companion matrix in `GL(3, K)` whose multiplicative order is `q^3 - 1`.

    This function is used in :func:`HughesPlane`.

    EXAMPLES::

        sage: from sage.combinat.designs.block_design import q3_minus_one_matrix
        sage: m = q3_minus_one_matrix(GF(3))
        sage: m.multiplicative_order() == 3**3 - 1
        True

        sage: m = q3_minus_one_matrix(GF(4, 'a'))
        sage: m.multiplicative_order() == 4**3 - 1
        True

        sage: m = q3_minus_one_matrix(GF(5))
        sage: m.multiplicative_order() == 5**3 - 1
        True

        sage: m = q3_minus_one_matrix(GF(9, 'a'))
        sage: m.multiplicative_order() == 9**3 - 1
        True
    """
    q = K.cardinality()
    M = MatrixSpace(K, 3)

    if q.is_prime():
        from sage.rings.finite_rings.conway_polynomials import conway_polynomial
        try:
            a,b,c,_ = conway_polynomial(q, 3)
        except RuntimeError:  # the polynomial is not in the database
            pass
        else:
            return M([0,0,-a,1,0,-b,0,1,-c])

    m = M()
    m[1,0] = m[2,1] = K.one()
    while True:
        m[0,2] = K._random_nonzero_element()
        m[1,2] = K.random_element()
        m[2,2] = K.random_element()
        if m.multiplicative_order() == q**3 - 1:
            return m


def normalize_hughes_plane_point(p, q):
    r"""
    Return the normalized form of point ``p`` as a 3-tuple.

    In the Hughes projective plane over the finite field `K`, all triples `(xk,
    yk, zk)` with `k \in K` represent the same point (where the multiplication
    is over the nearfield built from `K`). This function chooses a canonical
    representative among them.

    This function is used in :func:`HughesPlane`.

    INPUT:

    - ``p`` -- point with the coordinates `(x,y,z)` (a list, a vector, a tuple...)

    - ``q`` -- cardinality of the underlying finite field

    EXAMPLES::

        sage: from sage.combinat.designs.block_design import normalize_hughes_plane_point
        sage: K = FiniteField(9,'x')
        sage: x = K.gen()
        sage: normalize_hughes_plane_point((x, x + 1, x), 9)
        (1, x, 1)
        sage: normalize_hughes_plane_point(vector((x,x,x)), 9)
        (1, 1, 1)
        sage: zero = K.zero()
        sage: normalize_hughes_plane_point((2*x + 2, zero, zero), 9)
        (1, 0, 0)
        sage: one = K.one()
        sage: normalize_hughes_plane_point((2*x, one, zero), 9)
        (2*x, 1, 0)
    """
    for i in [2,1,0]:
        if p[i].is_one():
            return tuple(p)
        elif not p[i].is_zero():
            k = ~p[i]
            if k.is_square():
                return (p[0] * k,p[1] * k,p[2] * k)
            else:
                return ((p[0] * k)**q,(p[1]*k)**q,(p[2]*k)**q)


def HughesPlane(q2, check=True):
    r"""
    Return the Hughes projective plane of order ``q2``.

    Let `q` be an odd prime, the Hughes plane of order `q^2` is a finite
    projective plane of order `q^2` introduced by D. Hughes in [Hu57]_. Its
    construction is as follows.

    Let `K = GF(q^2)` be a finite field with `q^2` elements and `F = GF(q)
    \subset K` be its unique subfield with `q` elements. We define a twisted
    multiplication on `K` as

    .. MATH::

        x \circ y =
        \begin{cases}
        x\ y & \text{if y is a square in K}\\
        x^q\ y & \text{otherwise}
        \end{cases}

    The points of the Hughes plane are the triples `(x, y, z)` of points in `K^3
    \backslash \{0,0,0\}` up to the equivalence relation `(x,y,z) \sim (x \circ
    k, y \circ k, z \circ k)` where `k \in K`.

    For `a = 1` or `a \in (K \backslash F)` we define a block `L(a)` as the set of
    triples `(x,y,z)` so that `x + a \circ y + z = 0`. The rest of the blocks
    are obtained by letting act the group `GL(3, F)` by its standard action.

    For more information, see :wikipedia:`Hughes_plane` and [We07].

    .. SEEALSO::

        :func:`DesarguesianProjectivePlaneDesign` to build the Desarguesian
        projective planes

    INPUT:

    - ``q2`` -- an even power of an odd prime number

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.

    EXAMPLES::

        sage: H = designs.HughesPlane(9); H
        (91,10,1)-Balanced Incomplete Block Design

    We prove in the following computations that the Desarguesian plane ``H`` is
    not Desarguesian. Let us consider the two triangles `(0,1,10)` and `(57, 70,
    59)`. We show that the intersection points `D_{0,1} \cap D_{57,70}`,
    `D_{1,10} \cap D_{70,59}` and `D_{10,0} \cap D_{59,57}` are on the same line
    while `D_{0,70}`, `D_{1,59}` and `D_{10,57}` are not concurrent::

        sage: blocks = H.blocks()
        sage: line = lambda p,q: next(b for b in blocks if p in b and q in b)

        sage: b_0_1 = line(0, 1)
        sage: b_1_10 = line(1, 10)
        sage: b_10_0 = line(10, 0)
        sage: b_57_70 = line(57, 70)
        sage: b_70_59 = line(70, 59)
        sage: b_59_57 = line(59, 57)

        sage: set(b_0_1).intersection(b_57_70)
        {2}
        sage: set(b_1_10).intersection(b_70_59)
        {73}
        sage: set(b_10_0).intersection(b_59_57)
        {60}

        sage: line(2, 73) == line(73, 60)
        True

        sage: b_0_57 = line(0, 57)
        sage: b_1_70 = line(1, 70)
        sage: b_10_59 = line(10, 59)

        sage: p = set(b_0_57).intersection(b_1_70)
        sage: q = set(b_1_70).intersection(b_10_59)
        sage: p == q
        False

    TESTS:

    Some wrong input::

        sage: designs.HughesPlane(5)
        Traceback (most recent call last):
        ...
        EmptySetError: No Hughes plane of non-square order exists.

        sage: designs.HughesPlane(16)
        Traceback (most recent call last):
        ...
        EmptySetError: No Hughes plane of even order exists.

    Check that it works for non-prime `q`::

        sage: designs.HughesPlane(3**4)    # not tested - 10 secs
        (6643,82,1)-Balanced Incomplete Block Design
    """
    if not q2.is_square():
        raise EmptySetError("No Hughes plane of non-square order exists.")
    if q2 % 2 == 0:
        raise EmptySetError("No Hughes plane of even order exists.")
    q = q2.sqrt()
    K = FiniteField(q2, prefix='x')
    F = FiniteField(q, prefix='y')
    A = q3_minus_one_matrix(F)
    A = A.change_ring(K)
    m = K.list()
    V = VectorSpace(K, 3)
    zero = K.zero()
    one = K.one()
    points = [(x, y, one) for x in m for y in m] + \
             [(x, one, zero) for x in m] + \
             [(one, zero, zero)]
    relabel = {tuple(p):i for i,p in enumerate(points)}
    blcks = []
    for a in m:
        if a not in F or a == 1:
            # build L(a)
            aa = ~a
            l = []
            l.append(V((-a, one, zero)))
            for x in m:
                y = - aa * (x+one)
                if not y.is_square():
                    y *= aa**(q-1)
                l.append(V((x, y, one)))
            # compute the orbit of L(a)
            blcks.append([relabel[normalize_hughes_plane_point(p,q)] for p in l])
            for i in range(q2 + q):
                l = [A*j for j in l]
                blcks.append([relabel[normalize_hughes_plane_point(p,q)] for p in l])
    from .bibd import BalancedIncompleteBlockDesign
    return BalancedIncompleteBlockDesign(q2**2+q2+1, blcks, check=check)


def projective_plane_to_OA(pplane, pt=None, check=True):
    r"""
    Return the orthogonal array built from the projective plane ``pplane``.

    The orthogonal array `OA(n+1,n,2)` is obtained from the projective plane
    ``pplane`` by removing the point ``pt`` and the `n+1` lines that pass
    through it`. These `n+1` lines form the `n+1` groups while the remaining
    `n^2+n` lines form the transversals.

    INPUT:

    - ``pplane`` -- a projective plane as a 2-design

    - ``pt`` -- a point in the projective plane ``pplane``; if it is not provided,
      then it is set to `n^2 + n`

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.

    EXAMPLES::

        sage: from sage.combinat.designs.block_design import projective_plane_to_OA
        sage: p2 = designs.DesarguesianProjectivePlaneDesign(2, point_coordinates=False)
        sage: projective_plane_to_OA(p2)
        [[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0]]
        sage: p3 = designs.DesarguesianProjectivePlaneDesign(3, point_coordinates=False)
        sage: projective_plane_to_OA(p3)
        [[0, 0, 0, 0],
         [0, 1, 2, 1],
         [0, 2, 1, 2],
         [1, 0, 2, 2],
         [1, 1, 1, 0],
         [1, 2, 0, 1],
         [2, 0, 1, 1],
         [2, 1, 0, 2],
         [2, 2, 2, 0]]

        sage: pp = designs.DesarguesianProjectivePlaneDesign(16, point_coordinates=False)
        sage: _ = projective_plane_to_OA(pp, pt=0)
        sage: _ = projective_plane_to_OA(pp, pt=3)
        sage: _ = projective_plane_to_OA(pp, pt=7)
    """
    from .bibd import _relabel_bibd
    pplane = pplane.blocks()
    n = len(pplane[0]) - 1

    if pt is None:
        pt = n**2+n

    assert len(pplane) == n**2+n+1, "pplane is not a projective plane"
    assert all(len(B) == n+1 for B in pplane), "pplane is not a projective plane"

    pplane = _relabel_bibd(pplane,n**2+n+1,p=n**2+n)
    OA = [[x % n for x in sorted(X)] for X in pplane if n**2+n not in X]

    assert len(OA) == n**2, "pplane is not a projective plane"

    if check:
        from .designs_pyx import is_orthogonal_array
        is_orthogonal_array(OA,n+1,n,2)

    return OA


def projective_plane(n, check=True, existence=False):
    r"""
    Return a projective plane of order ``n`` as a 2-design.

    A finite projective plane is a 2-design with `n^2+n+1` lines (or blocks) and
    `n^2+n+1` points. For more information on finite projective planes, see the
    :wikipedia:`Projective_plane#Finite_projective_planes`.

    If no construction is possible, then the function raises a :exc:`EmptySetError`,
    whereas if no construction is available, the function raises a
    :exc:`NotImplementedError`.

    INPUT:

    - ``n`` -- the finite projective plane's order

    EXAMPLES::

        sage: # needs sage.schemes
        sage: designs.projective_plane(2)
        (7,3,1)-Balanced Incomplete Block Design
        sage: designs.projective_plane(3)
        (13,4,1)-Balanced Incomplete Block Design
        sage: designs.projective_plane(4)
        (21,5,1)-Balanced Incomplete Block Design
        sage: designs.projective_plane(5)
        (31,6,1)-Balanced Incomplete Block Design
        sage: designs.projective_plane(6)
        Traceback (most recent call last):
        ...
        EmptySetError: By the Bruck-Ryser theorem, no projective plane of order 6 exists.
        sage: designs.projective_plane(10)
        Traceback (most recent call last):
        ...
        EmptySetError: No projective plane of order 10 exists by
        C. Lam, L. Thiel and S. Swiercz "The nonexistence of finite
        projective planes of order 10" (1989), Canad. J. Math.
        sage: designs.projective_plane(12)
        Traceback (most recent call last):
        ...
        NotImplementedError: If such a projective plane exists,
        we do not know how to build it.
        sage: designs.projective_plane(14)
        Traceback (most recent call last):
        ...
        EmptySetError: By the Bruck-Ryser theorem, no projective plane of order 14 exists.

    TESTS::

        sage: # needs sage.schemes
        sage: designs.projective_plane(2197, existence=True)
        True
        sage: designs.projective_plane(6, existence=True)
        False
        sage: designs.projective_plane(10, existence=True)
        False
        sage: designs.projective_plane(12, existence=True)
        Unknown
    """
    from sage.combinat.designs.bibd import BruckRyserChowla_check

    if n <= 1:
        if existence:
            return False
        raise EmptySetError("There is no projective plane of order <= 1")

    if n == 10:
        if existence:
            return False
        ref = ("C. Lam, L. Thiel and S. Swiercz \"The nonexistence of finite "
               "projective planes of order 10\" (1989), Canad. J. Math.")
        raise EmptySetError("No projective plane of order 10 exists by %s" % ref)

    if BruckRyserChowla_check(n*n+n+1, n+1, 1) is False:
        if existence:
            return False
        raise EmptySetError("By the Bruck-Ryser theorem, no projective"
                         " plane of order {} exists.".format(n))

    if not is_prime_power(n):
        if existence:
            return Unknown
        raise NotImplementedError("If such a projective plane exists, we do "
                                  "not know how to build it.")

    if existence:
        return True
    else:
        return DesarguesianProjectivePlaneDesign(n, point_coordinates=False, check=check)


def AffineGeometryDesign(n, d, F, point_coordinates=True, check=True):
    r"""
    Return an affine geometry design.

    The affine geometry design `AG_d(n,q)` is the 2-design whose blocks are the
    `d`-vector subspaces in `\GF{q}^n`. It has parameters

    .. MATH::

        v = q^n,\ k = q^d,\ \lambda = \binom{n-1}{d-1}_q

    where the `q`-binomial coefficient `\binom{m}{r}_q` is defined by

    .. MATH::

        \binom{m}{r}_q = \frac{(q^m - 1)(q^{m-1} - 1) \cdots (q^{m-r+1}-1)}
              {(q^r-1)(q^{r-1}-1)\cdots (q-1)}

    .. SEEALSO::

        :func:`ProjectiveGeometryDesign`

    INPUT:

    - ``n`` -- integer; the Euclidean dimension. The number of points of the
      design is `v=|\GF{q}^n|`

    - ``d`` -- integer; the dimension of the (affine) subspaces of `\GF{q}^n`
      which make up the blocks

    - ``F`` -- a finite field or a prime power

    - ``point_coordinates`` -- boolean (default: ``True``); whether we use
      coordinates in `\GF{q}^n` or plain integers for the points of the design

    - ``check`` -- boolean (default: ``True``); whether to check the output

    EXAMPLES::

        sage: # needs sage.combinat
        sage: BD = designs.AffineGeometryDesign(3, 1, GF(2))
        sage: BD.is_t_design(return_parameters=True)
        (True, (2, 8, 2, 1))
        sage: BD = designs.AffineGeometryDesign(3, 2, GF(4))
        sage: BD.is_t_design(return_parameters=True)
        (True, (2, 64, 16, 5))
        sage: BD = designs.AffineGeometryDesign(4, 2, GF(3))
        sage: BD.is_t_design(return_parameters=True)
        (True, (2, 81, 9, 13))

    With ``F`` an integer instead of a finite field::

        sage: BD = designs.AffineGeometryDesign(3, 2, 4)
        sage: BD.is_t_design(return_parameters=True)
        (True, (2, 64, 16, 5))

    Testing the option ``point_coordinates``::

        sage: designs.AffineGeometryDesign(3, 1, GF(2),
        ....:                              point_coordinates=True).blocks()[0]
        [(0, 0, 0), (0, 0, 1)]
        sage: designs.AffineGeometryDesign(3, 1, GF(2),
        ....:                              point_coordinates=False).blocks()[0]
        [0, 1]
    """
    try:
        q = int(F)
    except TypeError:
        q = F.cardinality()
    else:
        from sage.rings.finite_rings.finite_field_constructor import GF
        F = GF(q)

    n = int(n)
    d = int(d)

    from itertools import islice
    from sage.combinat.q_analogues import q_binomial
    from sage.matrix.echelon_matrix import reduced_echelon_matrix_iterator

    points = {p:i for i,p in enumerate(reduced_echelon_matrix_iterator(F,1,n+1,copy=True,set_immutable=True)) if p[0,0]}

    blocks = []
    l1 = int(q_binomial(n+1, d+1, q) - q_binomial(n, d+1, q))
    l2 = q**d
    for m1 in islice(reduced_echelon_matrix_iterator(F,d+1,n+1,copy=False),
                     int(l1)):
        b = []
        for m2 in islice(reduced_echelon_matrix_iterator(F,1,d+1,copy=False),
                         int(l2)):
            m = m2*m1
            m.echelonize()
            m.set_immutable()
            b.append(points[m])
        blocks.append(b)

    B = BlockDesign(len(points), blocks, name='AffineGeometryDesign', check=check)

    if point_coordinates:
        rd = {i: p[0][1:] for p, i in points.items()}
        for v in rd.values():
            v.set_immutable()
        B.relabel(rd)

    if check:
        if not B.is_t_design(t=2, v=q**n, k=q**d, l=q_binomial(n-1, d-1, q)):
            raise RuntimeError("error in AffineGeometryDesign "
                    "construction. Please e-mail sage-devel@googlegroups.com")
    return B


def CremonaRichmondConfiguration():
    r"""
    Return the Cremona-Richmond configuration.

    The Cremona-Richmond configuration is a set system whose incidence graph
    is equal to the
    :meth:`~sage.graphs.graph_generators.GraphGenerators.TutteCoxeterGraph`. It
    is a generalized quadrangle of parameters `(2,2)`.

    For more information, see the
    :wikipedia:`Cremona-Richmond_configuration`.

    EXAMPLES::

        sage: H = designs.CremonaRichmondConfiguration(); H                             # needs networkx
        Incidence structure with 15 points and 15 blocks
        sage: g = graphs.TutteCoxeterGraph()                                            # needs networkx
        sage: H.incidence_graph().is_isomorphic(g)                                      # needs networkx
        True
    """
    from sage.graphs.generators.smallgraphs import TutteCoxeterGraph
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    g = TutteCoxeterGraph()
    H = IncidenceStructure([g.neighbors(v)
                            for v in g.bipartite_sets()[0]])
    H.relabel()
    return H


def WittDesign(n):
    """
    INPUT:

    - ``n`` -- integer in `9,10,11,12,21,22,23,24`

    Wraps GAP Design's WittDesign. If ``n=24`` then this function returns the
    large Witt design `W_{24}`, the unique (up to isomorphism) `5-(24,8,1)`
    design. If ``n=12`` then this function returns the small Witt design
    `W_{12}`, the unique (up to isomorphism) `5-(12,6,1)` design.  The other
    values of `n` return a block design derived from these.

    .. NOTE::

        Requires GAP's Design package (included in the gap_packages Sage spkg).

    EXAMPLES::

        sage: # optional - gap_package_design
        sage: BD = designs.WittDesign(9)
        sage: BD.is_t_design(return_parameters=True)
        (True, (2, 9, 3, 1))
        sage: BD
        Incidence structure with 9 points and 12 blocks
        sage: print(BD)
        Incidence structure with 9 points and 12 blocks
    """
    libgap.load_package("design")
    B = libgap.WittDesign(n)
    v = B['v'].sage()
    gB = [[x - 1 for x in b] for b in B['blocks'].sage()]
    return BlockDesign(v, gB, name='WittDesign', check=True)


def HadamardDesign(n):
    """
    As described in Section 1, p. 10, in [CvL]_. The input n must have the
    property that there is a Hadamard matrix of order `n+1` (and that a
    construction of that Hadamard matrix has been implemented...).

    EXAMPLES::

        sage: # needs sage.combinat sage.modules
        sage: designs.HadamardDesign(7)
        Incidence structure with 7 points and 7 blocks
        sage: print(designs.HadamardDesign(7))
        Incidence structure with 7 points and 7 blocks

    For example, the Hadamard 2-design with `n = 11` is a design whose parameters are `2-(11, 5, 2)`.
    We verify that `NJ = 5J` for this design. ::

        sage: # needs sage.combinat sage.modules
        sage: D = designs.HadamardDesign(11); N = D.incidence_matrix()
        sage: J = matrix(ZZ, 11, 11, [1]*11*11); N*J
        [5 5 5 5 5 5 5 5 5 5 5]
        [5 5 5 5 5 5 5 5 5 5 5]
        [5 5 5 5 5 5 5 5 5 5 5]
        [5 5 5 5 5 5 5 5 5 5 5]
        [5 5 5 5 5 5 5 5 5 5 5]
        [5 5 5 5 5 5 5 5 5 5 5]
        [5 5 5 5 5 5 5 5 5 5 5]
        [5 5 5 5 5 5 5 5 5 5 5]
        [5 5 5 5 5 5 5 5 5 5 5]
        [5 5 5 5 5 5 5 5 5 5 5]
        [5 5 5 5 5 5 5 5 5 5 5]

    REFERENCES:

    - [CvL] P. Cameron, J. H. van Lint, Designs, graphs, codes and
      their links, London Math. Soc., 1991.
    """
    from sage.combinat.matrices.hadamard_matrix import hadamard_matrix
    from sage.matrix.constructor import matrix
    H = hadamard_matrix(n+1) #assumed to be normalised.
    H1 = H.matrix_from_columns(range(1,n+1))
    H2 = H1.matrix_from_rows(range(1,n+1))
    J = matrix(ZZ,n,n,[1]*n*n)
    MS = J.parent()
    A = MS((H2+J)/2) # convert -1's to 0's; coerce entries to ZZ
    # A is the incidence matrix of the block design
    return IncidenceStructure(incidence_matrix=A,name='HadamardDesign')


def Hadamard3Design(n):
    r"""
    Return the Hadamard 3-design with parameters `3-(n, \frac n 2, \frac n 4 - 1)`.

    This is the unique extension of the Hadamard `2`-design (see
    :meth:`HadamardDesign`).  We implement the description from pp. 12 in
    [CvL]_.

    INPUT:

    - ``n`` -- integer; a multiple of 4 such that `n>4`

    EXAMPLES::

        sage: # needs sage.combinat sage.modules
        sage: designs.Hadamard3Design(12)
        Incidence structure with 12 points and 22 blocks

    We verify that any two blocks of the Hadamard `3`-design `3-(8, 4, 1)`
    design meet in `0` or `2` points. More generally, it is true that any two
    blocks of a Hadamard `3`-design meet in `0` or `\frac{n}{4}` points (for `n
    > 4`).

    ::

        sage: # needs sage.combinat sage.modules
        sage: D = designs.Hadamard3Design(8)
        sage: N = D.incidence_matrix()
        sage: N.transpose()*N
        [4 2 2 2 2 2 2 2 2 2 2 2 2 0]
        [2 4 2 2 2 2 2 2 2 2 2 2 0 2]
        [2 2 4 2 2 2 2 2 2 2 2 0 2 2]
        [2 2 2 4 2 2 2 2 2 2 0 2 2 2]
        [2 2 2 2 4 2 2 2 2 0 2 2 2 2]
        [2 2 2 2 2 4 2 2 0 2 2 2 2 2]
        [2 2 2 2 2 2 4 0 2 2 2 2 2 2]
        [2 2 2 2 2 2 0 4 2 2 2 2 2 2]
        [2 2 2 2 2 0 2 2 4 2 2 2 2 2]
        [2 2 2 2 0 2 2 2 2 4 2 2 2 2]
        [2 2 2 0 2 2 2 2 2 2 4 2 2 2]
        [2 2 0 2 2 2 2 2 2 2 2 4 2 2]
        [2 0 2 2 2 2 2 2 2 2 2 2 4 2]
        [0 2 2 2 2 2 2 2 2 2 2 2 2 4]


    REFERENCES:

    .. [CvL] \P. Cameron, J. H. van Lint, Designs, graphs, codes and
      their links, London Math. Soc., 1991.
    """
    if n == 1 or n == 4:
        raise ValueError("The Hadamard design with n = %s does not extend to a three design." % n)
    from sage.combinat.matrices.hadamard_matrix import hadamard_matrix
    from sage.matrix.constructor import matrix, block_matrix
    H = hadamard_matrix(n) #assumed to be normalised.
    H1 = H.matrix_from_columns(range(1, n))
    J = matrix(ZZ, n, n-1, [1]*(n-1)*n)
    A1 = (H1+J)/2
    A2 = (J-H1)/2
    A = block_matrix(1, 2, [A1, A2]) #the incidence matrix of the design.
    return IncidenceStructure(incidence_matrix=A, name='HadamardThreeDesign')
