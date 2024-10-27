# sage.doctest: needs sage.libs.gap
r"""
Database of generalised quadrangles with spread

This module implements some construction of generalised quadrangles
with spread.

EXAMPLES::

    sage: GQ, S = designs.generalised_quadrangle_with_spread(4, 16, check=False)
    sage: GQ
    Incidence structure with 325 points and 1105 blocks
    sage: GQ2, O = designs.generalised_quadrangle_hermitian_with_ovoid(4)
    sage: GQ2
    Incidence structure with 1105 points and 325 blocks
    sage: GQ3 = GQ.dual()
    sage: set(GQ3._points) == set(GQ2._points)
    True
    sage: GQ2.is_isomorphic(GQ3) # long time
    True

REFERENCES:

- [PT2009]_

- [TP1994]_

- :wikipedia:`Generalized_quadrangle`

AUTHORS:

- Ivo Maffei (2020-07-26): initial version
"""

# ****************************************************************************
#       Copyright (C) 2020 Ivo Maffei <ivomaffei@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


def generalised_quadrangle_with_spread(const int s, const int t,
                                       existence=False, check=True):
    r"""
    Construct a generalised quadrangle GQ of order `(s,t)` with a spread S.

    INPUT:

    - ``s``, ``t`` -- integers; order of the generalised quadrangle

    - ``existence`` -- boolean;

    - ``check`` -- boolean; if ``True``, then Sage checks that the object built
      is correct. (default: ``True``)

    OUTPUT:

    A pair `(GQ, S)` where `GQ` is a :class:`IncidenceStructure` representing
    the generalised quadrangle and `S` is a list of blocks of `GQ` representing
    the spread of `GQ`.

    EXAMPLES::

        sage: t = designs.generalised_quadrangle_with_spread(3, 9)
        sage: t[0]
        Incidence structure with 112 points and 280 blocks
        sage: designs.generalised_quadrangle_with_spread(5, 25, existence=True)
        True
        sage: (designs.generalised_quadrangle_with_spread(4, 16, check=False))[0]
        Incidence structure with 325 points and 1105 blocks
        sage: designs.generalised_quadrangle_with_spread(0, 2, existence=True)
        False

    REFERENCES:

    For more on generalised quadrangles and their spread see [PT2009]_ or
    [TP1994]_.

    TESTS::

        sage: GQ, S = designs.generalised_quadrangle_with_spread(2, 4)
        sage: GQ
        Incidence structure with 27 points and 45 blocks
        sage: designs.generalised_quadrangle_with_spread(3, 4)
        Traceback (most recent call last):
        ...
        RuntimeError: Sage can't build a GQ of order (3, 4) with a spread
        sage: designs.generalised_quadrangle_with_spread(3, 4, existence=True)
        Unknown
    """
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    from sage.misc.unknown import Unknown
    from sage.arith.misc import is_prime_power

    if s < 1 or t < 1:
        if existence:
            return False
        raise RuntimeError(f"No GQ of order ({s}, {t}) exists")

    if s == 1 and t == 1:  # we have a square
        if existence:
            return True
        D = IncidenceStructure([[0, 1], [1, 2], [2, 3], [3, 0]])
        return (D, [[0, 1], [2, 3]])

    if is_prime_power(s) and t == s * s:
        if existence:
            return True
        (GQ, S) = dual_GQ_ovoid(*generalised_quadrangle_hermitian_with_ovoid(s))
        if check:
            if not is_GQ_with_spread(GQ, S, s=s, t=t):
                raise RuntimeError("Sage built a wrong GQ with spread")
        return (GQ, S)

    if existence:
        return Unknown
    raise RuntimeError(
        f"Sage can't build a GQ of order ({s}, {t}) with a spread")


def is_GQ_with_spread(GQ, S, s=None, t=None):
    r"""
    Check if GQ is a generalised quadrangle of order `(s,t)` and
    check that S is a spread of GQ

    INPUT:

    - ``GQ`` -- IncidenceStructure; the incidence structure that is supposed to
      be a generalised quadrangle

    - ``S`` -- iterable; the spread of ``GQ`` as an
      iterable of the blocks of ``GQ``

    - ``s``, ``t`` -- integers (optional); if `(s,t)` are given, then we check that
      ``GQ`` has order `(s,t)`

    EXAMPLES::

        sage: from sage.combinat.designs.gen_quadrangles_with_spread import *
        sage: t = generalised_quadrangle_hermitian_with_ovoid(3)
        sage: is_GQ_with_spread(*t)
        Traceback (most recent call last):
        ...
        TypeError: 'int' object is not iterable
        sage: t = dual_GQ_ovoid(*t)
        sage: is_GQ_with_spread(*t)
        True
        sage: is_GQ_with_spread(*t, s=3)
        True

    TESTS::

       sage: from sage.combinat.designs.gen_quadrangles_with_spread import *
       sage: t = generalised_quadrangle_hermitian_with_ovoid(2)
       sage: t = dual_GQ_ovoid(*t)
       sage: is_GQ_with_spread(*t, s=2, t=4)
       True
       sage: is_GQ_with_spread(*t, s=2)
       True
       sage: is_GQ_with_spread(*t, s=3)
       False
    """
    res = GQ.is_generalized_quadrangle(parameters=True)
    if res is False \
       or (s is not None and s != res[0]) \
       or (t is not None and t != res[1]):
        return False

    # check spread
    return GQ.is_spread(S)


def dual_GQ_ovoid(GQ, O):
    r"""
    Compute the dual incidence structure of GQ
    and return the image of `O` under the dual map

    INPUT:

    - ``GQ`` -- IncidenceStructure; the generalised quadrangle we want
      the dual of

    - ``O`` -- iterable; the iterable of blocks we want to compute the dual

    OUTPUT: a pair ``(D, S)`` where ``D`` is the dual of ``GQ`` and
    ``S`` is the dual of ``O``

    EXAMPLES::

        sage: from sage.combinat.designs.gen_quadrangles_with_spread import \
        ....: dual_GQ_ovoid
        sage: t = designs.generalised_quadrangle_hermitian_with_ovoid(3)
        sage: t[0].is_generalized_quadrangle(parameters=True)
        (9, 3)
        sage: t = dual_GQ_ovoid(*t)
        sage: t[0].is_generalized_quadrangle(parameters=True)
        (3, 9)
        sage: all([x in t[0] for x in t[1]])
        True


    TESTS::

        sage: from sage.combinat.designs.gen_quadrangles_with_spread import *
        sage: t = designs.generalised_quadrangle_hermitian_with_ovoid(2)
        sage: t = dual_GQ_ovoid(*t)
        sage: t[0].is_generalized_quadrangle(parameters=True)
        (2, 4)
        sage: is_GQ_with_spread(*t)
        True
    """
    from sage.combinat.designs.incidence_structures import IncidenceStructure

    # GQ.ground_set()[i] becomes newBlocks[i]
    # GQ.blocks()[i] becomes i
    newBlocks = [[] for _ in range(GQ.num_points())]
    pointsToInt = {p: i for i, p in enumerate(GQ.ground_set())}

    for i, b in enumerate(GQ.blocks()):
        for p in b:
            newBlocks[pointsToInt[p]].append(i)

    S = [newBlocks[pointsToInt[p]] for p in O]

    D = IncidenceStructure(newBlocks)
    return (D, S)


def generalised_quadrangle_hermitian_with_ovoid(const int q):
    r"""
    Construct the generalised quadrangle `H(3,q^2)` with an ovoid.

    The GQ has order `(q^2,q)`.

    INPUT:

    - ``q`` -- integer; a prime power

    OUTPUT:

    A pair ``(D, O)`` where ``D`` is an IncidenceStructure representing the
    generalised quadrangle and ``O`` is a list of points of ``D`` which
    constitute an ovoid of ``D``.

    EXAMPLES::

        sage: t = designs.generalised_quadrangle_hermitian_with_ovoid(4)
        sage: t[0]
        Incidence structure with 1105 points and 325 blocks
        sage: len(t[1])
        65
        sage: G = t[0].intersection_graph([1])  # line graph
        sage: G.is_strongly_regular(True)
        (325, 68, 3, 17)
        sage: set(t[0].block_sizes())
        {17}

    REFERENCES:

    For more on `H(3,q^2)` and the construction implemented here see [PT2009]_
    or [TP1994]_.

    TESTS::

        sage: from sage.combinat.designs.gen_quadrangles_with_spread import \
              is_GQ_with_spread, dual_GQ_ovoid
        sage: t = designs.generalised_quadrangle_hermitian_with_ovoid(3)
        sage: t = dual_GQ_ovoid(*t)
        sage: is_GQ_with_spread(*t, s=3, t=9)
        True
        sage: t = dual_GQ_ovoid(*(
        ....: designs.generalised_quadrangle_hermitian_with_ovoid(2)))
        sage: t[0]
        Incidence structure with 27 points and 45 blocks
        sage: len(t[1])
        9
    """
    from sage.libs.gap.libgap import libgap
    from sage.combinat.designs.incidence_structures import IncidenceStructure
    from sage.arith.misc import is_prime_power

    GU = libgap.GU(4, q)
    Fq = libgap.GF(q * q)
    zero = libgap.Zero(Fq)
    one = libgap.One(Fq)
    V = libgap.FullRowSpace(Fq, 4)

    e1 = [one, zero, zero, zero]  # isotropic point

    points = libgap.Orbit(GU, e1, libgap.OnLines)  # all isotropic points
    pointInt = {tuple(x): i for i, x in enumerate(points)}

    GUp = libgap.Action(GU, points, libgap.OnLines)

    e2 = [zero, one, zero, zero]  # another isotropic point
    line = V.Subspace([e1, e2])  # totally isotropic line
    lineAsPoints = [tuple(libgap.Elements(libgap.Basis(b))[0])
                    for b in libgap.Elements(line.Subspaces(1))]
    line = libgap.Set([pointInt[p] + 1 for p in lineAsPoints])

    lines = libgap.Orbit(GUp, line, libgap.OnSets)  # all isotropic lines
    lines = [[int(x)-1 for x in b] for b in lines]  # convert to int
    # lines defines the GQ H(3, q^2)

    # to find an ovoid, we embed H(3,q^2) in H(4,q^2)
    # the embedding is (a,b,c,d) -> (a,b,0,c,d)
    # then we find a point in the latter and not in the former
    # this point will be collinear (in H(3,q^2)) to all points in an ovoid
    if q % 2 == 1:
        p, _ = is_prime_power(q, get_data=True)
        a = (p - 1) // 2
        aGap = zero
        for i in range(a):
            aGap += one
        p = [zero, one, one, aGap, zero]
    else:
        a = libgap.PrimitiveRoot(Fq)**(q-1)
        p = [zero, one, a+one, a, zero]

    J = [[0, 0, 0, 0, 1], [0, 0, 0, 1, 0], [0, 0, 1, 0, 0],
         [0, 1, 0, 0, 0], [1, 0, 0, 0, 0]]
    J = libgap(J)  # matrix of the invariant form of GU(5,q)

    # p' is collinear to p iff p'Jp^q = 0
    # note that p'Jp^q = bx^q + c where p' = (a,b,0,c,d) and p = (0,1,y,x,0)
    ovoid = []
    xq = p[3]**q
    for p2, i in pointInt.items():
        if p2[1]*xq + p2[2] == zero:  # p collinear to p2
            ovoid.append(i)

    D = IncidenceStructure(lines)
    return (D, ovoid)
