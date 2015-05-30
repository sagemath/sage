"""
Block designs

A *block design* is a set together with a family of subsets (repeated subsets
are allowed) whose members are chosen to satisfy some set of properties that are
deemed useful for a particular application. See :wikipedia:`Block_design`.

REFERENCES:

.. [1] Block design from wikipedia,
  :wikipedia:`Block_design`

.. [2] What is a block design?,
  http://designtheory.org/library/extrep/extrep-1.1-html/node4.html (in 'The
  External Representation of Block Designs' by Peter J. Cameron, Peter
  Dobcsanyi, John P. Morgan, Leonard H. Soicher)

.. [We07] Charles Weibel, "Survey of Non-Desarguesian planes" (2007), notices of
   the AMS, vol. 54 num. 10, pages 1294--1303

AUTHORS:

- Vincent Delecroix (2014): rewrite the part on projective planes :trac:`16281`

- Peter Dobcsanyi and David Joyner (2007-2008)

  This is a significantly modified form of the module block_design.py (version
  0.6) written by Peter Dobcsanyi peter@designtheory.org. Thanks go to Robert
  Miller for lots of good design suggestions.

.. TODO::

    Implement finite non-Desarguesian plane as in [We07]_ and
    :wikipedia:`Non-Desarguesian_plane`.

Functions and methods
---------------------
"""
#***************************************************************************
#                              Copyright (C) 2007                          #
#                                                                          #
#                Peter Dobcsanyi       and         David Joyner            #
#           <peter@designtheory.org>          <wdjoyner@gmail.com>         #
#                                                                          #
#                                                                          #
#    Distributed under the terms of the GNU General Public License (GPL)   #
#    as published by the Free Software Foundation; either version 2 of     #
#    the License, or (at your option) any later version.                   #
#                    http://www.gnu.org/licenses/                          #
#***************************************************************************

from sage.modules.free_module import VectorSpace
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.arith import binomial, integer_floor
from incidence_structures import IncidenceStructure
from sage.misc.decorators import rename_keyword
from sage.rings.finite_rings.constructor import FiniteField
from sage.categories.sets_cat import EmptySetError
from sage.misc.unknown import Unknown

BlockDesign = IncidenceStructure

###  utility functions  -------------------------------------------------------

def tdesign_params(t, v, k, L):
    """
    Return the design's parameters: `(t, v, b, r , k, L)`. Note that `t` must be
    given.

    EXAMPLES::

        sage: BD = designs.BlockDesign(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
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

    - ``v,k,lmbda`` (integers)

    OUTPUT:

    - a boolean or, if ``return_parameters`` is set to ``True`` a pair
      ``(True, (q,d))`` or ``(False, (None,None))``.

    EXAMPLES::

        sage: from sage.combinat.designs.block_design import are_hyperplanes_in_projective_geometry_parameters
        sage: are_hyperplanes_in_projective_geometry_parameters(40,13,4)
        True
        sage: are_hyperplanes_in_projective_geometry_parameters(40,13,4,return_parameters=True)
        (True, (3, 3))
        sage: PG = designs.ProjectiveGeometryDesign(3,2,GF(3))
        sage: PG.is_t_design(return_parameters=True)
        (True, (2, 40, 13, 4))

        sage: are_hyperplanes_in_projective_geometry_parameters(15,3,1)
        False
        sage: are_hyperplanes_in_projective_geometry_parameters(15,3,1,return_parameters=True)
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
    import sage.rings.arith as arith

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

def ProjectiveGeometryDesign(n, d, F, algorithm=None, check=True):
    """
    Return a projective geometry design.

    A projective geometry design of parameters `n,d,F` has for points the lines
    of `F^{n+1}`, and for blocks the `d+1`-dimensional subspaces of `F^{n+1}`,
    each of which contains `\\frac {|F|^{d+1}-1} {|F|-1}` lines.

    INPUT:

    - ``n`` is the projective dimension

    - ``d`` is the dimension of the subspaces of `P = PPn(F)` which
      make up the blocks.

    - ``F`` is a finite field.

    - ``algorithm`` -- set to ``None`` by default, which results in using Sage's
      own implementation. In order to use GAP's implementation instead (i.e. its
      ``PGPointFlatBlockDesign`` function) set ``algorithm="gap"``. Note that
      GAP's "design" package must be available in this case, and that it can be
      installed with the ``gap_packages`` spkg.

    EXAMPLES:

    The set of `d`-dimensional subspaces in a `n`-dimensional projective space
    forms `2`-designs (or balanced incomplete block designs)::

        sage: PG = designs.ProjectiveGeometryDesign(4,2,GF(2))
        sage: PG
        Incidence structure with 31 points and 155 blocks
        sage: PG.is_t_design(return_parameters=True)
        (True, (2, 31, 7, 7))

        sage: PG = designs.ProjectiveGeometryDesign(3,1,GF(4,'z'))
        sage: PG.is_t_design(return_parameters=True)
        (True, (2, 85, 5, 1))

    Check that the constructor using gap also works::

        sage: BD = designs.ProjectiveGeometryDesign(2, 1, GF(2), algorithm="gap") # optional - gap_packages (design package)
        sage: BD.is_t_design(return_parameters=True)                              # optional - gap_packages (design package)
        (True, (2, 7, 3, 1))
    """
    if algorithm is None:
        from sage.matrix.echelon_matrix import reduced_echelon_matrix_iterator
        from copy import copy

        points = {}
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
        return BlockDesign(len(points), blocks, name="ProjectiveGeometryDesign", check=check)
    if algorithm == "gap":   # Requires GAP's Design
        from sage.interfaces.gap import gap
        gap.load_package("design")
        gap.eval("D := PGPointFlatBlockDesign( %s, %s, %d )"%(n,F.order(),d))
        v = eval(gap.eval("D.v"))
        gblcks = eval(gap.eval("D.blocks"))
        gB = []
        for b in gblcks:
            gB.append([x-1 for x in b])
        return BlockDesign(v, gB, name="ProjectiveGeometryDesign", check=check)

def DesarguesianProjectivePlaneDesign(n, check=True):
    r"""
    Return the Desarguesian projective plane of order ``n`` as a 2-design.

    The Desarguesian projective plane of order `n` can also be defined as the
    projective plane over a field of order `n`. For more information, have a
    look at :wikipedia:`Projective_plane`.

    INPUT:

    - ``n`` -- an integer which must be a power of a prime number

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to
      ``True`` by default.

    .. SEEALSO::

        :func:`ProjectiveGeometryDesign`

    EXAMPLES::

        sage: designs.DesarguesianProjectivePlaneDesign(2)
        Incidence structure with 7 points and 7 blocks
        sage: designs.DesarguesianProjectivePlaneDesign(3)
        Incidence structure with 13 points and 13 blocks
        sage: designs.DesarguesianProjectivePlaneDesign(4)
        Incidence structure with 21 points and 21 blocks
        sage: designs.DesarguesianProjectivePlaneDesign(5)
        Incidence structure with 31 points and 31 blocks
        sage: designs.DesarguesianProjectivePlaneDesign(6)
        Traceback (most recent call last):
        ...
        ValueError: the order of a finite field must be a prime power
    """
    K = FiniteField(n, 'x')
    n2 = n**2
    relabel = {x:i for i,x in enumerate(K)}
    Kiter = relabel  # it is much faster to iterate throug a dict than through
                     # the finite field K

    # we decompose the (equivalence class) of points [x:y:z] of the projective
    # plane into an affine plane, an affine line and a point. At the same time,
    # we relabel the points with the integers from 0 to n^2 + n as follows:
    # - the affine plane is the set of points [x:y:1] (i.e. the third coordinate
    #   is non-zero) and gets relabeled from 0 to n^2-1
    affine_plane   = lambda x,y: relabel[x] + n * relabel[y]

    # - the affine line is the set of points [x:1:0] (i.e. the third coordinate is
    #   zero but not the second one) and gets relabeld from n^2 to n^2 + n - 1
    line_infinity  = lambda x: n2 + relabel[x]

    # - the point is [1:0:0] and gets relabeld n^2 + n
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

    return BlockDesign(n2+n+1, blcks, name="Desarguesian projective plane of order %d"%n, check=check)

def projective_plane_to_OA(pplane, pt=None, check=True):
    r"""
    Return the orthogonal array built from the projective plane ``pplane``.

    The orthogonal array `OA(n+1,n,2)` is obtained from the projective plane
    ``pplane`` by removing the point ``pt`` and the `n+1` lines that pass
    through it`. These `n+1` lines form the `n+1` groups while the remaining
    `n^2+n` lines form the transversals.

    INPUT:

    - ``pplane`` - a projective plane as a 2-design

    - ``pt`` - a point in the projective plane ``pplane``. If it is not provided
      then it is set to `n^2 + n`.

    - ``check`` -- (boolean) Whether to check that output is correct before
      returning it. As this is expected to be useless (but we are cautious
      guys), you may want to disable it whenever you want speed. Set to
      ``True`` by default.

    EXAMPLES::

        sage: from sage.combinat.designs.block_design import projective_plane_to_OA
        sage: p2 = designs.DesarguesianProjectivePlaneDesign(2)
        sage: projective_plane_to_OA(p2)
        [[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0]]
        sage: p3 = designs.DesarguesianProjectivePlaneDesign(3)
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

        sage: pp = designs.DesarguesianProjectivePlaneDesign(16)
        sage: _ = projective_plane_to_OA(pp, pt=0)
        sage: _ = projective_plane_to_OA(pp, pt=3)
        sage: _ = projective_plane_to_OA(pp, pt=7)
    """
    from bibd import _relabel_bibd
    pplane = pplane.blocks()
    n = len(pplane[0]) - 1

    if pt is None:
        pt = n**2+n

    assert len(pplane) == n**2+n+1, "pplane is not a projective plane"
    assert all(len(B) == n+1 for B in pplane), "pplane is not a projective plane"

    pplane = _relabel_bibd(pplane,n**2+n+1,p=n**2+n)
    OA = [[x%n for x in sorted(X)] for X in pplane if not n**2+n in X]

    assert len(OA) == n**2, "pplane is not a projective plane"

    if check:
        from designs_pyx import is_orthogonal_array
        is_orthogonal_array(OA,n+1,n,2)

    return OA

def projective_plane(n, check=True, existence=False):
    r"""
    Return a projective plane of order ``n`` as a 2-design.

    A finite projective plane is a 2-design with `n^2+n+1` lines (or blocks) and
    `n^2+n+1` points. For more information on finite projective planes, see the
    :wikipedia:`Projective_plane#Finite_projective_planes`.

    If no construction is possible, then the function raises a ``EmptySetError``
    whereas if no construction is available the function raises a
    ``NotImplementedError``.

    INPUT:

    - ``n`` -- the finite projective plane's order

    EXAMPLES::

        sage: designs.projective_plane(2)
        Incidence structure with 7 points and 7 blocks
        sage: designs.projective_plane(3)
        Incidence structure with 13 points and 13 blocks
        sage: designs.projective_plane(4)
        Incidence structure with 21 points and 21 blocks
        sage: designs.projective_plane(5)
        Incidence structure with 31 points and 31 blocks
        sage: designs.projective_plane(6)
        Traceback (most recent call last):
        ...
        EmptySetError: By the Bruck-Ryser theorem, no projective plane of order 6 exists.
        sage: designs.projective_plane(10)
        Traceback (most recent call last):
        ...
        EmptySetError: No projective plane of order 10 exists by C. Lam, L. Thiel and S. Swiercz "The nonexistence of finite projective planes of order 10" (1989), Canad. J. Math.
        sage: designs.projective_plane(12)
        Traceback (most recent call last):
        ...
        NotImplementedError: If such a projective plane exists, we do not know how to build it.
        sage: designs.projective_plane(14)
        Traceback (most recent call last):
        ...
        EmptySetError: By the Bruck-Ryser theorem, no projective plane of order 14 exists.

    TESTS::

        sage: designs.projective_plane(2197, existence=True)
        True
        sage: designs.projective_plane(6, existence=True)
        False
        sage: designs.projective_plane(10, existence=True)
        False
        sage: designs.projective_plane(12, existence=True)
        Unknown
    """
    from sage.rings.arith import is_prime_power
    from sage.rings.sum_of_squares import is_sum_of_two_squares_pyx

    if n <= 1:
        if existence:
            return False
        raise EmptySetError("There is no projective plane of order <= 1")

    if n == 10:
        if existence:
            return False
        ref = ("C. Lam, L. Thiel and S. Swiercz \"The nonexistence of finite "
               "projective planes of order 10\" (1989), Canad. J. Math.")
        raise EmptySetError("No projective plane of order 10 exists by %s"%ref)

    if (n%4) in [1,2] and not is_sum_of_two_squares_pyx(n):
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
        return DesarguesianProjectivePlaneDesign(n, check=check)

def AffineGeometryDesign(n, d, F):
    r"""
    Return an Affine Geometry Design.

    INPUT:

    - `n` (integer) -- the Euclidean dimension. The number of points is
      `v=|F^n|`.

    - `d` (integer) -- the dimension of the (affine) subspaces of `P = GF(q)^n`
      which make up the blocks.

    - `F` -- a Finite Field (i.e. ``FiniteField(17)``), or a prime power
      (i.e. an integer)

    `AG_{n,d} (F)`, as it is sometimes denoted, is a `2` - `(v, k, \lambda)`
    design of points and `d`- flats (cosets of dimension `n`) in the affine
    geometry `AG_n (F)`, where

    .. math::

             v = q^n,\  k = q^d ,
             \lambda =\frac{(q^{n-1}-1) \cdots (q^{n+1-d}-1)}{(q^{n-1}-1) \cdots (q-1)}.

    Wraps some functions used in GAP Design's ``PGPointFlatBlockDesign``.  Does
    *not* require GAP's Design package.

    EXAMPLES::

        sage: BD = designs.AffineGeometryDesign(3, 1, GF(2))
        sage: BD.is_t_design(return_parameters=True)
        (True, (2, 8, 2, 1))
        sage: BD = designs.AffineGeometryDesign(3, 2, GF(2))
        sage: BD.is_t_design(return_parameters=True)
        (True, (3, 8, 4, 1))

    With an integer instead of a Finite Field::

        sage: BD = designs.AffineGeometryDesign(3, 2, 4)
        sage: BD.is_t_design(return_parameters=True)
        (True, (2, 64, 16, 5))
    """
    try:
        q = int(F)
    except TypeError:
        q = F.order()

    from sage.interfaces.gap import gap, GapElement
    from sage.sets.set import Set
    gap.eval("V:=GaloisField(%s)^%s"%(q,n))
    gap.eval("points:=AsSet(V)")
    gap.eval("Subs:=AsSet(Subspaces(V,%s));"%d)
    gap.eval("CP:=Cartesian(points,Subs)")
    flats = gap.eval("flats:=List(CP,x->Sum(x))") # affine spaces
    gblcks = eval(gap.eval("Set(List(flats,f->Filtered([1..Length(points)],i->points[i] in f)));"))
    v = q**n
    gB = []
    for b in gblcks:
       gB.append([x-1 for x in b])
    return BlockDesign(v, gB, name="AffineGeometryDesign")

def WittDesign(n):
    """
    INPUT:

    - ``n`` is in `9,10,11,12,21,22,23,24`.

    Wraps GAP Design's WittDesign. If ``n=24`` then this function returns the
    large Witt design `W_{24}`, the unique (up to isomorphism) `5-(24,8,1)`
    design. If ``n=12`` then this function returns the small Witt design
    `W_{12}`, the unique (up to isomorphism) `5-(12,6,1)` design.  The other
    values of `n` return a block design derived from these.

    .. NOTE:

        Requires GAP's Design package (included in the gap_packages Sage spkg).

    EXAMPLES::

        sage: BD = designs.WittDesign(9)             # optional - gap_packages (design package)
        sage: BD.is_t_design(return_parameters=True) # optional - gap_packages (design package)
        (True, (2, 9, 3, 1))
        sage: BD                             # optional - gap_packages (design package)
        Incidence structure with 9 points and 12 blocks
        sage: print BD                       # optional - gap_packages (design package)
        WittDesign<points=[0, 1, 2, 3, 4, 5, 6, 7, 8], blocks=[[0, 1, 7], [0, 2, 5], [0, 3, 4], [0, 6, 8], [1, 2, 6], [1, 3, 5], [1, 4, 8], [2, 3, 8], [2, 4, 7], [3, 6, 7], [4, 5, 6], [5, 7, 8]]>
    """
    from sage.interfaces.gap import gap, GapElement
    gap.load_package("design")
    gap.eval("B:=WittDesign(%s)"%n)
    v = eval(gap.eval("B.v"))
    gblcks = eval(gap.eval("B.blocks"))
    gB = []
    for b in gblcks:
       gB.append([x-1 for x in b])
    return BlockDesign(v, gB, name="WittDesign", check=True)

def HadamardDesign(n):
    """
    As described in Section 1, p. 10, in [CvL]. The input n must have the
    property that there is a Hadamard matrix of order `n+1` (and that a
    construction of that Hadamard matrix has been implemented...).

    EXAMPLES::

        sage: designs.HadamardDesign(7)
        Incidence structure with 7 points and 7 blocks
        sage: print designs.HadamardDesign(7)
        HadamardDesign<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>

    For example, the Hadamard 2-design with `n = 11` is a design whose parameters are 2-(11, 5, 2).
    We verify that `NJ = 5J` for this design. ::
     
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
    return IncidenceStructure(incidence_matrix=A,name="HadamardDesign")

def Hadamard3Design(n):
    """
    Return the Hadamard 3-design with parameters `3-(n, \\frac n 2, \\frac n 4 - 1)`.

    This is the unique extension of the Hadamard `2`-design (see
    :meth:`HadamardDesign`).  We implement the description from pp. 12 in
    [CvL]_.

    INPUT:

    - ``n`` (integer) -- a multiple of 4 such that `n>4`.

    EXAMPLES::

        sage: designs.Hadamard3Design(12)
        Incidence structure with 12 points and 22 blocks

    We verify that any two blocks of the Hadamard `3`-design `3-(8, 4, 1)`
    design meet in `0` or `2` points. More generally, it is true that any two
    blocks of a Hadamard `3`-design meet in `0` or `\\frac{n}{4}` points (for `n
    > 4`).

    ::

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

    .. [CvL] P. Cameron, J. H. van Lint, Designs, graphs, codes and
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
    return IncidenceStructure(incidence_matrix=A, name="HadamardThreeDesign")
