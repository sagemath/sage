# sage.doctest: needs sage.rings.finite_rings
r"""
Orthogonal arrays (build recursive constructions)

This module implements several constructions of
:mod:`Orthogonal Arrays<sage.combinat.designs.orthogonal_arrays>`.
As their input can be complex, they all have a counterpart in the
:mod:`~sage.combinat.designs.orthogonal_arrays_find_recursive` module
that automatically computes it.

All these constructions are automatically queried when the
:func:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array` function is
called.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_3` | Return an `OA(k,nm+i)`.
    :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_4` | Return a `OA(k,nm+rs)`.
    :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_5` | Return an `OA(k,nm+r+s+t)`.
    :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_3_6` | Return a `OA(k,nm+i)`.
    :func:`~sage.combinat.designs.orthogonal_arrays_build_recursive.construction_q_x` | Return an `OA(k,(q-1)*(q-x)+x+2)` using the `q-x` construction.
    :func:`OA_and_oval` | Return a `OA(q+1,q)` whose blocks contains `\leq 2` zeroes in the last `q` columns.
    :func:`thwart_lemma_3_5` | Return an `OA(k,nm+a+b+c+d)`.
    :func:`thwart_lemma_4_1` | Return an `OA(k,nm+4(n-2))`.
    :func:`three_factor_product` | Return an `OA(k+1,n_1n_2n_3)`.
    :func:`brouwer_separable_design` | Return a `OA(k,t(q^2+q+1)+x)` using Brouwer's result on separable designs.

Functions
---------
"""
from itertools import repeat
from .orthogonal_arrays import orthogonal_array, wilson_construction, is_orthogonal_array


def construction_3_3(k,n,m,i,explain_construction=False):
    r"""
    Return an `OA(k,nm+i)`.

    This is Wilson's construction with `i` truncated columns of size 1 and such
    that a block `B_0` of the incomplete OA intersects all truncated columns. As
    a consequence, all other blocks intersect only `0` or `1` of the last `i`
    columns. This allow to consider the block `B_0` only up to its first `k`
    coordinates and then use a `OA(k,i)` instead of a `OA(k,m+i) - i.OA(k,1)`.

    This is construction 3.3 from [AC07]_.

    INPUT:

    - ``k``, ``n``, ``m``, ``i`` -- integers such that the following designs are
      available: `OA(k,n)`, `OA(k,m)`, `OA(k,m+1)`, `OA(k,r)`

    - ``explain_construction`` -- boolean; return a string describing
      the construction

    .. SEEALSO::

        :func:`~sage.combinat.designs.orthogonal_arrays_find_recursive.find_construction_3_3`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_construction_3_3
        sage: from sage.combinat.designs.orthogonal_arrays_build_recursive import construction_3_3
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: k = 11; n = 177
        sage: is_orthogonal_array(construction_3_3(*find_construction_3_3(k,n)[1]),k,n,2)           # needs sage.schemes
        True

        sage: print(designs.orthogonal_arrays.explain_construction(9,91))
        Construction 3.3 with n=11,m=8,i=3 from:
           Julian R. Abel, Nicholas Cavenagh
           Concerning eight mutually orthogonal latin squares,
           Vol. 15, n.3, pp. 255-261,
           Journal of Combinatorial Designs, 2007
    """
    from .orthogonal_arrays import wilson_construction, OA_relabel, incomplete_orthogonal_array
    if explain_construction:
        return (("Construction 3.3 with n={},m={},i={} from:\n"
                 "  Julian R. Abel, Nicholas Cavenagh\n" +
                 "  Concerning eight mutually orthogonal latin squares,\n" +
                 "  Vol. 15, n.3, pp. 255-261,\n" +
                 "  Journal of Combinatorial Designs, 2007").format(n,m,i))

    # Builds an OA(k+i,n) containing a block [0]*(k+i)
    OA = incomplete_orthogonal_array(k+i,n,(1,))
    OA = [[(x+1) % n for x in B] for B in OA]

    # Truncated version
    OA = [B[:k]+[0 if x == 0 else None for x in B[k:]] for B in OA]

    OA = wilson_construction(OA,k,n,m,[1]*i,check=False)[:-i]
    matrix = [list(range(m)) + list(range(n*m, n*m+i))] * k
    OA.extend(OA_relabel(orthogonal_array(k,m+i),k,m+i,matrix=matrix))
    assert is_orthogonal_array(OA,k,n*m+i)
    return OA


def construction_3_4(k,n,m,r,s,explain_construction=False):
    r"""
    Return a `OA(k,nm+rs)`.

    This is Wilson's construction applied to a truncated `OA(k+r+1,n)` with `r`
    columns of size `1` and one column of size `s`.

    The unique elements of the `r` truncated columns are picked so that a block
    `B_0` contains them all.

    - If there exists an `OA(k,m+r+1)` the column of size `s` is truncated in
      order to intersect `B_0`.

    - Otherwise, if there exists an `OA(k,m+r)`, the last column must not
      intersect `B_0`

    This is construction 3.4 from [AC07]_.

    INPUT:

    - ``k``, ``n``, ``m``, ``r``, ``s`` -- integers; we assume that `s<n` and
      `1\leq r,s`

      The following designs must be available: `OA(k,n)`, `OA(k,m)`,
      `OA(k,m+1)`, `OA(k,m+2)`, `OA(k,s)`. Additionally, it requires either a
      `OA(k,m+r)` or a `OA(k,m+r+1)`.

    - ``explain_construction`` -- boolean; return a string describing
      the construction

    .. SEEALSO::

        :func:`~sage.combinat.designs.orthogonal_arrays_find_recursive.find_construction_3_4`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_construction_3_4
        sage: from sage.combinat.designs.orthogonal_arrays_build_recursive import construction_3_4
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: k = 8; n = 196
        sage: is_orthogonal_array(construction_3_4(*find_construction_3_4(k,n)[1]),k,n,2)           # needs sage.schemes
        True

        sage: print(designs.orthogonal_arrays.explain_construction(8,164))
        Construction 3.4 with n=23,m=7,r=2,s=1 from:
           Julian R. Abel, Nicholas Cavenagh
           Concerning eight mutually orthogonal latin squares,
           Vol. 15, n.3, pp. 255-261,
           Journal of Combinatorial Designs, 2007
    """
    if explain_construction:
        return ("Construction 3.4 with n={},m={},r={},s={} from:\n" +
                 "  Julian R. Abel, Nicholas Cavenagh\n" +
                 "  Concerning eight mutually orthogonal latin squares,\n" +
                 "  Vol. 15, n.3, pp. 255-261,\n" +
                 "  Journal of Combinatorial Designs, 2007").format(n,m,r,s)

    from .orthogonal_arrays import wilson_construction, OA_relabel
    assert s < n
    master_design = orthogonal_array(k+r+1,n)

    # Defines the first k+r columns of the matrix of labels
    matrix = [list(range(n))] * k + [[None]*n]*(r) + [[None]*n]
    B0 = master_design[0]
    for i in range(k,k+r):
        matrix[i][B0[i]] = 0

    # Last column
    if orthogonal_array(k, m+r  ,existence=True):
        last_group = [x for x in range(s+1) if x != B0[-1]][:s]
    elif orthogonal_array(k,m+r+1,existence=True):
        last_group = [x for x in range(s+1) if x != B0[-1]][:s-1] + [B0[-1]]
    else:
        raise RuntimeError

    for i, x in enumerate(last_group):
        matrix[-1][x] = i

    OA = OA_relabel(master_design,k+r+1,n, matrix=matrix)
    OA = wilson_construction(OA,k,n,m,[1]*r+[s],check=False)
    return OA


def construction_3_5(k,n,m,r,s,t,explain_construction=False):
    r"""
    Return an `OA(k,nm+r+s+t)`.

    This is exactly Wilson's construction with three truncated groups
    except we make sure that all blocks have size `>k`, so we don't
    need a `OA(k,m+0)` but only `OA(k,m+1)`, `OA(k,m+2)` ,`OA(k,m+3)`.

    This is construction 3.5 from [AC07]_.

    INPUT:

    - ``k``, ``n``, ``m`` -- integers

    - ``r``, ``s``, ``t`` -- integers; sizes of the three truncated groups,
      such that `r\leq s` and `(q-r-1)(q-s) \geq (q-s-1)*(q-r)`

    - ``explain_construction`` -- boolean; return a string describing
      the construction

    The following designs must be available : `OA(k,n)`, `OA(k,r)`, `OA(k,s)`,
    `OA(k,t)`, `OA(k,m+1)`, `OA(k,m+2)`, `OA(k,m+3)`.

    .. SEEALSO::

        :func:`~sage.combinat.designs.orthogonal_arrays_find_recursive.find_construction_3_5`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_construction_3_5
        sage: from sage.combinat.designs.orthogonal_arrays_build_recursive import construction_3_5
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: k=8;n=111
        sage: is_orthogonal_array(construction_3_5(*find_construction_3_5(k,n)[1]),k,n,2)           # needs sage.schemes
        True

        sage: print(designs.orthogonal_arrays.explain_construction(8,90))
        Construction 3.5 with n=11,m=6,r=8,s=8,t=8 from:
           Julian R. Abel, Nicholas Cavenagh
           Concerning eight mutually orthogonal latin squares,
           Vol. 15, n.3, pp. 255-261,
           Journal of Combinatorial Designs, 2007
    """
    from .orthogonal_arrays import wilson_construction, OA_relabel
    assert r <= s
    q = n
    assert (q-r-1)*(q-s) >= (q-s-1)*(q-r)

    if explain_construction:
        return (("Construction 3.5 with n={},m={},r={},s={},t={} from:\n"
                 "  Julian R. Abel, Nicholas Cavenagh\n" +
                 "  Concerning eight mutually orthogonal latin squares,\n" +
                 "  Vol. 15, n.3, pp. 255-261,\n" +
                 "  Journal of Combinatorial Designs, 2007").format(n,m,r,s,t))

    master_design = orthogonal_array(k+3,q)

    # group k+1 has cardinality r
    # group k+2 has cardinality s
    # group k+3 has cardinality t

    # Taking q-s blocks going through 0 in the last block
    blocks_crossing_0 = [B[-3:] for B in master_design if B[-1] == 0][:q-s]

    # defining the undeleted points of the groups k+1,k+2
    group_k_1 = [x[0] for x in blocks_crossing_0]
    group_k_1 = [x for x in range(q) if x not in group_k_1][:r]

    group_k_2 = [x[1] for x in blocks_crossing_0]
    group_k_2 = [x for x in range(q) if x not in group_k_2][:s]

    # All blocks that have a deleted point in groups k+1 and k+2 MUST contain a
    # point in group k+3
    group_k_3 = [B[-1] for B in master_design if B[-3] not in group_k_1 and B[-2] not in group_k_2]
    group_k_3 = list(set(group_k_3))
    assert len(group_k_3) <= t
    group_k_3.extend([x for x in range(q) if x not in group_k_3])
    group_k_3 = group_k_3[:t]

    # Relabelling the OA
    r1 = [None]*q
    r2 = [None]*q
    r3 = [None]*q
    for i, x in enumerate(group_k_1):
        r1[x] = i
    for i, x in enumerate(group_k_2):
        r2[x] = i
    for i, x in enumerate(group_k_3):
        r3[x] = i

    OA = OA_relabel(master_design, k+3,q, matrix=[list(range(q))]*k+[r1,r2,r3])
    OA = wilson_construction(OA,k,q,m,[r,s,t], check=False)
    return OA


def construction_3_6(k,n,m,i,explain_construction=False):
    r"""
    Return a `OA(k,nm+i)`.

    This is Wilson's construction with `r` columns of order `1`, in which each
    block intersects at most two truncated columns. Such a design exists when
    `n` is a prime power and is returned by :func:`OA_and_oval`.

    INPUT:

    - ``k``, ``n``, ``m``, ``i`` -- integers; `n` must be a prime power. The
      following designs must be available: `OA(k+r,q)`, `OA(k,m)`, `OA(k,m+1)`,
      `OA(k,m+2)`

    - ``explain_construction`` -- boolean; return a string describing
      the construction

    This is construction 3.6 from [AC07]_.

    .. SEEALSO::

        - :func:`~sage.combinat.designs.orthogonal_arrays_find_recursive.find_construction_3_6`

        - :func:`OA_and_oval`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_find_recursive import find_construction_3_6
        sage: from sage.combinat.designs.orthogonal_arrays_build_recursive import construction_3_6
        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: k=8;n=95
        sage: is_orthogonal_array(construction_3_6(*find_construction_3_6(k,n)[1]),k,n,2)           # needs sage.schemes
        True

        sage: print(designs.orthogonal_arrays.explain_construction(10,756))
        Construction 3.6 with n=16,m=47,i=4 from:
           Julian R. Abel, Nicholas Cavenagh
           Concerning eight mutually orthogonal latin squares,
           Vol. 15, n.3, pp. 255-261,
           Journal of Combinatorial Designs, 2007
    """
    if explain_construction:
        return (("Construction 3.6 with n={},m={},i={} from:\n"
                 "  Julian R. Abel, Nicholas Cavenagh\n" +
                 "  Concerning eight mutually orthogonal latin squares,\n" +
                 "  Vol. 15, n.3, pp. 255-261,\n" +
                 "  Journal of Combinatorial Designs, 2007").format(n,m,i))

    from .orthogonal_arrays import wilson_construction
    OA = OA_and_oval(n)
    OA = [B[:k+i] for B in OA]
    OA = [B[:k] + [x if x == 0 else None for x in B[k:]] for B in OA]
    OA = wilson_construction(OA,k,n,m,[1]*i)
    assert is_orthogonal_array(OA,k,n*m+i)
    return OA


def OA_and_oval(q, *, solver=None, integrality_tolerance=1e-3):
    r"""
    Return a `OA(q+1,q)` whose blocks contains `\leq 2` zeroes in the last `q`
    columns.

    This `OA` is build from a projective plane of order `q`, in which there
    exists an oval `O` of size `q+1` (i.e. a set of `q+1` points no three of
    which are [colinear/contained in a common set of the projective plane]).

    Removing an element `x\in O` and all sets that contain it, we obtain a
    `TD(q+1,q)` in which `O` intersects all columns except one. As `O` is an
    oval, no block of the `TD` intersects it more than twice.

    INPUT:

    - ``q`` -- a prime power

    - ``solver`` -- (default: ``None``) specify a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``integrality_tolerance`` -- parameter for use with MILP solvers over an
      inexact base ring; see :meth:`MixedIntegerLinearProgram.get_values`

    .. NOTE::

        This function is called by :func:`construction_3_6`, an implementation
        of Construction 3.6 from [AC07]_.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_build_recursive import OA_and_oval
        sage: _ = OA_and_oval
    """
    from sage.arith.misc import is_prime_power
    from sage.combinat.designs.block_design import projective_plane
    from .orthogonal_arrays import OA_relabel

    assert is_prime_power(q)
    B = projective_plane(q, check=False)

    # We compute the oval with a linear program
    from sage.numerical.mip import MixedIntegerLinearProgram
    p = MixedIntegerLinearProgram(solver=solver)
    b = p.new_variable(binary=True)
    V = B.ground_set()
    p.add_constraint(p.sum([b[i] for i in V]) == q+1)
    for bl in B:
        p.add_constraint(p.sum([b[i] for i in bl]) <= 2)
    p.solve()
    b = p.get_values(b, convert=bool, tolerance=integrality_tolerance)
    oval = [x for x,i in b.items() if i]
    assert len(oval) == q+1

    # We remove one element from the oval
    x = oval.pop()
    oval.sort()

    # We build the TD by relabelling the point set, and removing those which
    # contain x.
    r = {}

    # Make sure that the first set containing x in B is the one
    # which contains no other oval point
    B = sorted(B, key=lambda b: any(xx in oval for xx in b))

    BB = []
    for b in B:
        if x in b:
            for xx in b:
                if xx == x:
                    continue
                r[xx] = len(r)
        else:
            BB.append(b)

    assert len(r) == (q+1)*q # all points except x have an image
    assert len(set(r.values())) == len(r) # the images are different

    # Relabelling/sorting the blocks and the oval
    BB = [[r[xx] for xx in b] for b in BB]
    oval = [r[xx] for xx in oval]

    for b in BB:
        b.sort()
    oval.sort()

    # Turning the TD into an OA
    BB = [[xx % q for xx in b] for b in BB]
    oval = [xx % q for xx in oval]
    assert len(oval) == q

    # We relabel the "oval" as relabelled as [0,...,0]
    OA = OA_relabel(BB+([[0]+oval]),q+1,q,blocks=[[0]+oval])
    OA = [[(x+1) % q for x in B] for B in OA]
    OA.remove([0]*(q+1))

    assert all(sum([xx == 0 for xx in b[1:]]) <= 2 for b in OA)
    return OA


def construction_q_x(k, q, x, check=True, explain_construction=False):
    r"""
    Return an `OA(k,(q-1)*(q-x)+x+2)` using the `q-x` construction.

    Let `v=(q-1)*(q-x)+x+2`. If there exists a projective plane of order `q`
    (e.g. when `q` is a prime power) and `0<x<q` then there exists a
    `(v-1,\{q-x-1,q-x+1\})`-GDD of type `(q-1)^{q-x}(x+1)^1` (see [Greig99]_ or
    Theorem 2.50, section IV.2.3 of [DesignHandbook]_). By adding to the ground
    set one point contained in all groups of the GDD, one obtains a
    `(v,\{q-x-1,q-x+1,q,x+2\})`-PBD with exactly one set of size `x+2`.

    Thus, assuming that we have the following:

    - `OA(k,q-x-1)-(q-x-1).OA(k,1)`
    - `OA(k,q-x+1)-(q-x+1).OA(k,1)`
    - `OA(k,q)-q.OA(k,1)`
    - `OA(k,x+2)`

    Then we can build from the PBD an `OA(k,v)`.

    Construction of the PBD (shared by Julian R. Abel):

        Start with a resolvable `(q^2,q,1)`-BIBD and put the points into a `q\times q`
        array so that rows form a parallel class and columns form another.

        Now delete:

        - All `x(q-1)` points from the first `x` columns and not in the first
          row

        - All `q-x` points in the last `q-x` columns AND the first row.

        Then add a point `p_1` to the blocks that are rows. Add a second point
        `p_2` to the `q-x` blocks that are columns of size `q-1`, plus the first
        row of size `x+1`.

    INPUT:

    - ``k``, ``q``, ``x`` -- integers such that `0<x<q` and such that Sage can
      build:

        - A projective plane of order `q`
        - `OA(k,q-x-1)-(q-x-1).OA(k,1)`
        - `OA(k,q-x+1)-(q-x+1).OA(k,1)`
        - `OA(k,q)-q.OA(k,1)`
        - `OA(k,x+2)`

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.

    - ``explain_construction`` -- boolean; return a string describing
      the construction

    .. SEEALSO::

        - :func:`~sage.combinat.designs.orthogonal_arrays_find_recursive.find_q_x`
        - :func:`~sage.combinat.designs.block_design.projective_plane`
        - :func:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array`
        - :func:`~sage.combinat.designs.orthogonal_arrays.OA_from_PBD`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_build_recursive import construction_q_x
        sage: _ = construction_q_x(9,16,6)                                                          # needs sage.schemes

        sage: print(designs.orthogonal_arrays.explain_construction(9,158))
        (q-x)-construction with q=16,x=6 from:
           Malcolm Greig,
           Designs from projective planes and PBD bases,
           vol. 7, num. 5, pp. 341--374,
           Journal of Combinatorial Designs, 1999

    REFERENCES:

    .. [Greig99] Designs from projective planes and PBD bases
      Malcolm Greig
      Journal of Combinatorial Designs
      vol. 7, num. 5, pp. 341--374
      1999
    """
    from sage.combinat.designs.orthogonal_arrays import OA_from_PBD
    from sage.combinat.designs.orthogonal_arrays import incomplete_orthogonal_array

    if explain_construction:
        return ("(q-x)-construction with q={},x={} from:\n" +
                "   Malcolm Greig,\n" +
                "   Designs from projective planes and PBD bases,\n" +
                "   vol. 7, num. 5, pp. 341--374,\n" +
                "   Journal of Combinatorial Designs, 1999").format(q, x)

    n = (q-1)*(q-x)+x+2

    # We obtain the qxq matrix from a OA(q,q)-q.OA(1,q). We will need to add
    # blocks corresponding to the rows/columns
    OA = incomplete_orthogonal_array(q,q,(1,)*q)
    TD = [[i*q+xx for i, xx in enumerate(B)] for B in OA]

    # Add rows, extended with p1 and p2
    p1 = q**2
    p2 = p1 + 1
    TD.extend([ii*q + i for ii in range(q)] + [p1] for i in range(1, q))
    TD.append([ii*q for ii in range(q)] + [p1, p2])

    # Add Columns. We do not add some columns which would have size 1 after we
    # delete points.
    #
    # TD.extend([range(i*q,(i+1)*q) for i in range(x)])
    TD.extend([list(range(i*q,(i+1)*q))+[p2] for i in range(x,q)])

    points_to_delete = set([i*q+j for i in range(x) for j in range(1,q)]+[i*q for i in range(x,q)])
    points_to_keep = set(range(q**2+2))-points_to_delete
    relabel = {i:j for j,i in enumerate(points_to_keep)}

    # PBD is a (n,[q,q-x-1,q-x+1,x+2])-PBD
    PBD = [[relabel[xx] for xx in B if xx not in points_to_delete] for B in TD]

    # Taking the unique block of size x+2
    assert list(map(len,PBD)).count(x+2) == 1
    for B in PBD:
        if len(B) == x+2:
            break

    # We call OA_from_PBD without the block of size x+2 as there may not exist a
    # OA(k,x+2)-(x+2).OA(k,1)
    PBD.remove(B)
    OA = OA_from_PBD(k,(q-1)*(q-x)+x+2,PBD,check=False)

    # Filling the hole
    for xx in B:
        OA.remove([xx]*k)

    for BB in orthogonal_array(k, x+2):
        OA.append([B[x] for x in BB])

    if check:
        assert is_orthogonal_array(OA,k,n,2)

    return OA


def thwart_lemma_3_5(k,n,m,a,b,c,d=0,complement=False,explain_construction=False):
    r"""
    Return an `OA(k,nm+a+b+c+d)`.

    *(When `d=0`)*

    According to [Thwarts]_ when `n` is a prime power and `a+b+c\leq n+1`, one
    can build an `OA(k+3,n)` with three truncated columns of sizes `a,b,c` in
    such a way that all blocks have size `\leq k+2`.

    (in order to build a `OA(k,nm+a+b+c)` the following designs must also exist:
    `OA(k,a)`, `OA(k,b)`, `OA(k,c)`, `OA(k,m+0)`, `OA(k,m+1)`, `OA(k,m+2)`)

    Considering the complement of each truncated column, it is also possible to
    build an `OA(k+3,n)` with three truncated columns of sizes `a,b,c` in such a
    way that all blocks have size `>k` whenever `(n-a)+(n-b)+(n-c)\leq n+1`.

    (in order to build a `OA(k,nm+a+b+c)` the following designs must also exist:
    `OA(k,a)`, `OA(k,b)`, `OA(k,c)`, `OA(k,m+1)`, `OA(k,m+2)`, `OA(k,m+3)`)

    Here is the proof of Lemma 3.5 from [Thwarts]_ enriched with explanations
    from Julian R. Abel:

        For any prime power `n` one can build `k-1` MOLS by associating to every
        nonzero `x\in \mathbb F_n` the latin square:

        .. MATH::

            M_x(i,j) = i+x*j \text{ where }i,j\in \mathbb F_n

        In particular `M_1(i,j)=i+j`, whose `n` columns and lines are indexed by
        the elements of `\mathbb F_n`. If we order the elements of `\mathbb F_n`
        as `0,1,...,n-1,x+0,...,x+n-1,x^2+0,...` and reorder the columns
        and lines of `M_1` accordingly, the top-left `a\times b` squares
        contains at most `a+b-1` distinct symbols.

    *(When* `d\neq 0` *)*

    If there exists an `OA(k+3,n)` with three truncated columns of sizes `a,b,c`
    in such a way that all blocks have size `\leq k+2`, by truncating
    arbitrarily another column to size `d` one obtains an `OA` with 4 truncated
    columns whose blocks miss at least one value. Thus, following the proof
    again one can build an `OA(k+4)` with four truncated columns of sizes
    `a,b,c,d` with blocks of size `\leq k+3`.

    (in order to build a `OA(k,nm+a+b+c+d)` the following designs must also
    exist: `OA(k,a)`, `OA(k,b)`, `OA(k,c)`, `OA(k,d)`, `OA(k,m+0)`, `OA(k,m+1)`,
    `OA(k,m+2)`, `OA(k,m+3)`)

    As before, this also shows that one can build an `OA(k+4,n)` with four
    truncated columns of sizes `a,b,c,d` in such a way that all blocks have size
    `>k` whenever `(n-a)+(n-b)+(n-c)\leq n+1`

    (in order to build a `OA(k,nm+a+b+c+d)` the following designs must also
    exist: `OA(k,n-a)`, `OA(k,n-b)`, `OA(k,n-c)`, `OA(k,d)`, `OA(k,m+1)`,
    `OA(k,m+2)`, `OA(k,m+3)`, `OA(k,m+4)`)

    INPUT:

    - ``k``, ``n``, ``m``, ``a``, ``b``, ``c``, ``d`` -- integers which must
      satisfy the constraints above. In particular, `a+b+c\leq n+1` must hold
      By default, `d=0`.

    - ``complement`` -- boolean; whether to complement the sets, i.e. follow
      the `n-a,n-b,n-c` variant described above

    - ``explain_construction`` -- boolean; return a string describing
      the construction

    .. SEEALSO::

        - :func:`~sage.combinat.designs.orthogonal_arrays_find_recursive.find_thwart_lemma_3_5`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_build_recursive import thwart_lemma_3_5
        sage: from sage.combinat.designs.designs_pyx import is_orthogonal_array
        sage: OA = thwart_lemma_3_5(6,23,7,5,7,8)                                                   # needs sage.schemes
        sage: is_orthogonal_array(OA,6,23*7+5+7+8,2)                                                # needs sage.schemes
        True

        sage: print(designs.orthogonal_arrays.explain_construction(10,408))                         # needs sage.schemes
        Lemma 4.1 with n=13,m=28 from:
           Charles J.Colbourn, Jeffrey H. Dinitz, Mieczyslaw Wojtas,
           Thwarts in transversal designs,
           Designs, Codes and Cryptography 5, no. 3 (1995): 189-197.

    With sets of parameters from [Thwarts]_::

        sage: l = [
        ....:    [11, 27, 78, 16, 17, 25, 0],
        ....:    [12, 19, 208, 11, 13, 16, 0],
        ....:    [12, 19, 208, 13, 13, 16, 0],
        ....:    [10, 13, 78, 9, 9, 13, 1],
        ....:    [10, 13, 79, 9, 9, 13, 1]]
        sage: for k,n,m,a,b,c,d in l:                                       # not tested -- too long
        ....:     OA = thwart_lemma_3_5(k,n,m,a,b,c,d,complement=True)
        ....:     assert is_orthogonal_array(OA,k,n*m+a+b+c+d,verbose=True)

        sage: print(designs.orthogonal_arrays.explain_construction(10,1046))
        Lemma 3.5 with n=13,m=79,a=9,b=1,c=0,d=9 from:
           Charles J.Colbourn, Jeffrey H. Dinitz, Mieczyslaw Wojtas,
           Thwarts in transversal designs,
           Designs, Codes and Cryptography 5, no. 3 (1995): 189-197.

    REFERENCE:

    .. [Thwarts] Thwarts in transversal designs
      Charles J.Colbourn, Jeffrey H. Dinitz, Mieczyslaw Wojtas.
      Designs, Codes and Cryptography 5, no. 3 (1995): 189-197.
    """
    from sage.arith.misc import is_prime_power
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF

    if complement:
        a,b,c = n-a,n-b,n-c

    if explain_construction:
        return ("Lemma 3.5 with n={},m={},a={},b={},c={},d={} from:\n" +
                "   Charles J.Colbourn, Jeffrey H. Dinitz, Mieczyslaw Wojtas,\n" +
                "   Thwarts in transversal designs,\n" +
                "   Designs, Codes and Cryptography 5, no. 3 (1995): 189-197.").format(n,m,a,b,c,d)

    assert is_prime_power(n), "n(={}) must be a prime power".format(n)
    assert a <= n and b <= n and c <= n and d <= n, "a,b,c,d (={},{},{},{}) must be <=n(={})".format(a,b,c,d,n)
    assert a+b+c <= n+1, "{}={}+{}+{}=a+b+c>n+1={}+1 violates the assumptions".format(a+b+c,a,b,c,n)
    assert k+3+bool(d) <= n+1, "There exists no OA({},{}).".format(k+3+bool(d),n)
    G = GF(n,prefix='x')
    G_set = sorted(G) # sorted by lexicographic order, G[1] = 1
    assert G_set[0] == G.zero() and G_set[1] == G.one(), "problem with the ordering of {}".format(G)
    G_to_int = {v:i for i,v in enumerate(G_set)}

    # Builds an OA(n+1,n) whose last n-1 columns are
    #
    # \forall x \in G and x!=0, C_x(i,j) = i+x*j
    #
    # (only the necessary columns are built)
    OA = [[G_to_int[i+x*j] for i in G_set for j in G_set] for x in G_set[1:k+2+bool(d)]]
    # Adding the first two trivial columns
    OA.insert(0,[j for i in range(n) for j in range(n)])
    OA.insert(0,[i for i in range(n) for j in range(n)])
    OA = sorted(zip(*OA))

    # Moves the first three columns to the end
    OA = [list(B[3:]+B[:3]) for B in OA]

    # Set of values in the axb square
    third_complement = set(B[-1] for B in OA if B[-3] < a and B[-2] < b)

    assert n - len(third_complement) >= c

    # The keepers
    first_set = list(range(a))
    second_set = list(range(b))
    third_set = [x for x in range(n) if x not in third_complement][:c]

    last_sets = [first_set, second_set, third_set]

    if complement:
        last_sets = [set(range(n)).difference(s) for s in last_sets]

    sizes = [len(_) for _ in last_sets]
    last_sets_dict = [{v:i for i,v in enumerate(s)} for s in last_sets]

    # Truncating the OA
    for i,D in enumerate(last_sets_dict):
        kk = len(OA[0])-3+i
        for R in OA:
            R[kk] = D[R[kk]] if R[kk] in D else None

    if d:
        for R in OA:
            if R[-4] >= d:
                R[-4] = None
        sizes.insert(0,d)

    return wilson_construction(OA,k,n,m,sizes, check=False)


def thwart_lemma_4_1(k,n,m,explain_construction=False):
    r"""
    Return an `OA(k,nm+4(n-2))`.

    Implements Lemma 4.1 from [Thwarts]_.

        If `n\equiv 0,1\pmod{3}` is a prime power, then there exists a truncated
        `OA(n+1,n)` whose last four columns have size `n-2` and intersect every
        block on `1,3` or `4` values. Consequently, if there exists an
        `OA(k,m+1)`, `OA(k,m+3)`, `OA(k,m+4)` and a `OA(k,n-2)` then there
        exists an `OA(k,nm+4(n-2)`

        Proof: form the transversal design by removing one point of the
        `AG(2,3)` (Affine Geometry) contained in the Desarguesian Projective
        Plane `PG(2,n)`.

    The affine geometry on 9 points contained in the projective geometry
    `PG(2,n)` is given explicitly in [OS64]_ (Thanks to Julian R. Abel for
    finding the reference!).

    INPUT:

    - ``k``, ``n``, ``m`` -- integers

    - ``explain_construction`` -- boolean; return a string describing
      the construction

    .. SEEALSO::

        - :func:`~sage.combinat.designs.orthogonal_arrays_find_recursive.find_thwart_lemma_4_1`

    EXAMPLES::

        sage: print(designs.orthogonal_arrays.explain_construction(10,408))                         # needs sage.schemes
        Lemma 4.1 with n=13,m=28 from:
           Charles J.Colbourn, Jeffrey H. Dinitz, Mieczyslaw Wojtas,
           Thwarts in transversal designs,
           Designs, Codes and Cryptography 5, no. 3 (1995): 189-197.


    REFERENCES:

    .. [OS64] Finite projective planes with affine subplanes,
      T. G. Ostrom and F. A. Sherk.
      Canad. Math. Bull vol7 num.4 (1964)
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    from sage.arith.misc import is_prime_power
    from .block_design import DesarguesianProjectivePlaneDesign
    from itertools import chain

    if explain_construction:
        return ("Lemma 4.1 with n={},m={} from:\n" +
                "   Charles J.Colbourn, Jeffrey H. Dinitz, Mieczyslaw Wojtas,\n" +
                "   Thwarts in transversal designs,\n" +
                "   Designs, Codes and Cryptography 5, no. 3 (1995): 189-197.").format(n,m)

    assert is_prime_power(n), "n(={}) must be a prime power"
    assert k+4 <= n+1

    q = n
    K = FiniteField(q, 'x')
    relabel = {x: i for i, x in enumerate(K)}
    PG = DesarguesianProjectivePlaneDesign(q, check=False,
                                           point_coordinates=False).blocks()

    if q % 3 == 0:
        t = K.one()
    elif q % 3 == 1:
        t = K.multiplicative_generator()**((q - 1)//3)
    else:
        raise ValueError("q(={}) must be congruent to 0 or 1 mod 3".format(q))

    # The projective plane is labelled with integer coordinates. This code
    # relabels to integers the following points (given by homogeneous
    # coordinates in the projective space):
    #
    # - (1+t,t,1+t), (1,1,1), (1+t,t,t), (1,1,2), (0,0,1), (1,0,1), (0,1,1+t),
    #   (0,1,1), (1,0,-t)
    points = [(1+t,t,1+t), (1,1,1), (1+t,t,t), (1,1,2), (0,0,1), (1,0,1), (0,1,1+t), (0,1,1), (1,0,-t)]
    points = [[K(_) for _ in t] for t in points] # triples of K^3
    AG_2_3 = []
    for x,y,z in points:
        if z != 0:
            x,y,z = x/z,y/z,z/z
            AG_2_3.append(relabel[x]+n*relabel[y])
        elif y != 0:
            x,y,z = x/y,y/y,z
            AG_2_3.append(q**2+relabel[x])
        else:
            AG_2_3.append(q**2+q)

    AG_2_3 = set(AG_2_3)

    # All blocks of PG should intersect 'AG_2_3' on !=2 AG_2_3.
    assert all(len(AG_2_3.intersection(B)) != 2 for B in PG)

    p = list(AG_2_3)[0]
    # We now build a TD from the PG by removing p, in such a way that the last
    # two elements of the last 4 columns are elements of AG_2_3
    blocks = []
    columns = []
    for B in PG:
        if p not in B:
            blocks.append(B)
        else:
            B.remove(p)
            columns.append(B)

    # The columns containing elements from the AG are the last ones, and those
    # elements should be the last two
    columns.sort(key=lambda x:len(AG_2_3.intersection(x)))
    for i in range(4):
        columns[-i-1].sort(key=lambda x: int(x in AG_2_3))

    relabel = {v:i for i,v in enumerate(chain(columns))}

    TD = [sorted(relabel[x] for x in B) for B in blocks]

    # We build the OA, removing unnecessary columns
    OA = [[x % q for x in B[-k-4:]] for B in TD]
    for B in OA:
        for i in range(4):
            if B[k+i] >= n-2:
                B[k+i] = None

    return wilson_construction(OA,k,n,m,[n-2,]*4,check=False)


def three_factor_product(k,n1,n2,n3,check=False,explain_construction=False):
    r"""
    Return an `OA(k+1,n_1n_2n_3)`.

    The three factor product construction from [DukesLing14]_ does the following:

        If `n_1\leq n_2\leq n_3` are such that there exists an
        `OA(k,n_1)`, `OA(k+1,n_2)` and `OA(k+1,n_3)`, then there exists a
        `OA(k+1,n_1n_2n_3)`.

    It works with a modified product of orthogonal arrays ([Rees93]_, [Rees00]_)
    which keeps track of parallel classes in the `OA` (the definition is given
    for transversal designs).

        A subset of blocks in an `TD(k,n)` is called a `c`-parallel class if
        every point is covered exactly `c` times. A 1-parallel class is a
        parallel class.

    The modified product:

        If there exists an `OA(k,n_1)`, and if there exists an `OA(k,n_2)` whose
        blocks are partitionned into `s` `n_1`-parallel classes and `n_2-sn_1`
        parallel classes, then there exists an `OA(k,n_1n_2)` whose blocks can
        be partitionned into `sn_1^2` parallel classes and
        `(n_1n_2-sn_1^2)/n_1=n_2-sn_1` `n_1`-parallel classes.

        Proof:

        - The product of the blocks of a parallel class with an `OA(k,n_1)`
          yields an `n_1`-parallel class of an `OA(k,n_1n_2)`.

        - The product of the blocks of a `n_1`-parallel class of `OA(k,n_2)`
          with an `OA(k,n_1)` can be done in such a way that it yields `n_1n_2`
          parallel classes of `OA(k,n_1n_2)`. Those classes cover exactly the
          pairs that would have been covered with the usual product.

          This can be achieved by simple cyclic permutations. Let us build the
          product of the `n_1`-parallel class `\mathcal P\subseteq OA(k,n_2)`
          with `OA(k,n_1)`: when computing the product of `P\in\mathcal P` with
          `B^1\in OA(k,n_1)` the `i`-th coordinate should not be `(B^1_i,P_i)`
          but `(B^1_i+r,P_i)` (the sum is mod `n_1`) where `r` is the number of
          blocks of `\mathcal P` we have already processed whose `i`-th
          coordinate is equal to `P_i` (note that `r< n_1` as `\mathcal P` is
          `n_1`-parallel).

    With these tools, one can obtain the designs promised by the three factors
    construction applied to `k,n_1,n_2,n_3` (thanks to Julian R. Abel's help):

        1) Let `s` be the largest integer `\leq n_3/n_1`. Apply the product
           construction to `OA(k,n_1)` and a resolvable `OA(k,n_3)` whose blocks
           are partitionned into `s` `n_1`-parallel classes and `n_3-sn_1`
           parallel classes. It results in a `OA(k,n_1n_3)` partitionned into
           `sn_1^2` parallel classes plus `(n_1n_3-sn_1^2)/n_1=n_3-sn_1`
           `n_1`-parallel classes.

        2) Add `n_3-n_1` parallel classes to every `n_1`-parallel class to turn
           them into `n_3`-parallel classes. Apply the product construction to
           this partitionned `OA(k,n_1n_3)` with a resolvable `OA(k,n_2)`.

        3) As `OA(k,n_2)` is resolvable, the `n_2`-parallel classes of
           `OA(k,n_1n_2n_3)` are actually the union of `n_2` parallel classes,
           thus the `OA(k,n_1n_2n_3)` is resolvable and can be turned into an
           `OA(k+1,n_1n_2n_3)`

    INPUT:

    - ``k``, ``n1``, ``n2``, ``n3`` -- integers

    - ``check`` -- boolean; whether to check that everything is going smoothly
      while the design is being built. It is disabled by default, as the
      constructor of orthogonal arrays checks the final design anyway.

    - ``explain_construction`` -- boolean; return a string describing
      the construction

    .. SEEALSO::

        - :func:`~sage.combinat.designs.orthogonal_arrays_find_recursive.find_three_factor_product`

    EXAMPLES::

        sage: # needs sage.schemes
        sage: from sage.combinat.designs.designs_pyx import is_orthogonal_array
        sage: from sage.combinat.designs.orthogonal_arrays_build_recursive import three_factor_product
        sage: OA = three_factor_product(4,4,4,4)
        sage: is_orthogonal_array(OA,5,64)
        True
        sage: OA = three_factor_product(4,3,4,5)
        sage: is_orthogonal_array(OA,5,60)
        True
        sage: OA = three_factor_product(5,4,5,7)
        sage: is_orthogonal_array(OA,6,140)
        True
        sage: OA = three_factor_product(9,8,9,9) # long time
        sage: is_orthogonal_array(OA,10,8*9*9)   # long time
        True
        sage: print(designs.orthogonal_arrays.explain_construction(10,648))
        Three-factor product with n=8.9.9 from:
           Peter J. Dukes, Alan C.H. Ling,
           A three-factor product construction for mutually orthogonal latin squares,
           https://arxiv.org/abs/1401.1466

    REFERENCE:

    .. [DukesLing14] A three-factor product construction for mutually orthogonal latin squares,
      Peter J. Dukes, Alan C.H. Ling,
      :arxiv:`1401.1466`

    .. [Rees00] Truncated Transversal Designs: A New Lower Bound on the Number of Idempotent MOLS of Side,
      Rolf S. Rees,
      Journal of Combinatorial Theory, Series A 90.2 (2000): 257-266.

    .. [Rees93] Two new direct product-type constructions for resolvable group-divisible designs,
      Rolf S. Rees,
      Journal of Combinatorial Designs 1.1 (1993): 15-26.
    """
    assert n1 <= n2 and n2 <= n3

    if explain_construction:
        return ("Three-factor product with n={}.{}.{} from:\n" +
                "   Peter J. Dukes, Alan C.H. Ling,\n" +
                "   A three-factor product construction for mutually orthogonal latin squares,\n" +
                "   https://arxiv.org/abs/1401.1466").format(n1, n2, n3)

    def assert_c_partition(classs,k,n,c):
        r"""
        Makes sure that ``classs`` contains blocks `B` of size `k` such that the list of
        ``B[i]`` covers `[n]` exactly `c` times for every index `i`.
        """
        c = int(c)
        assert all(len(B) == k for B in classs), "A block has length {}!=k(={})".format(len(B),k)
        assert len(classs) == n*c, "not the right number of blocks"
        for p in zip(*classs):
            assert all(x == i//c for i,x in enumerate(sorted(p))), "A class is not c(={})-parallel".format(c)

    def product_with_parallel_classes(OA1,k,g1,g2,g1_parall,parall,check=True):
        r"""
        Return the product of two OA while keeping track of parallel classes.

        INPUT:

        - ``OA1`` -- (an `OA(k,g_1)`

        - ``k``, ``g1``, ``g2`` -- integers

        - ``g1_parall`` -- list of `g_1`-parallel classes

        - ``parall`` -- list of parallel classes

        .. NOTE::

            The list ``g1_parall+parall`` should be an `OA(k,g_2)`

        OUTPUT:

        Two lists of classes ``g1_parall`` and ``parallel`` which are respectively
        `g_1`-parallel and parallel classes such that ``g1_parall+parallel`` is an
        ``OA(k,g1*g2)``.
        """
        if check:
            for classs in g1_parall:
                assert_c_partition(classs,k,g2,g1)
            for classs in parall:
                assert_c_partition(classs,k,g2,1)

        # New parallel classes, built from a g1-parallel class with shifted copies
        # of OA1

        new_parallel_classes = []
        for classs2 in g1_parall:

            # Keep track of how many times we saw each point of [k]x[g2]
            count = [[0]*g2 for _ in range(k)]

            copies_of_OA1 = []
            for B2 in classs2:
                copy_of_OA1 = []

                shift = [count[i][x2] for i,x2 in enumerate(B2)]
                assert max(shift) < g1

                for B1 in OA1:
                    copy_of_OA1.append([x2*g1+(x1+sh) % g1 for sh,x1,x2 in zip(shift,B1,B2)])

                copies_of_OA1.append(copy_of_OA1)

                # Update the counts
                for i,x2 in enumerate(B2):
                    count[i][x2] += 1

            new_parallel_classes.extend([list(_) for _ in zip(*copies_of_OA1)])

        # New g1-parallel classes, each one built from the product of a parallel
        # class with a OA1

        new_g1_parallel_classes = []
        for classs2 in parall:
            disjoint_copies_of_OA1 = []
            for B2 in classs2:
                for B1 in OA1:
                    disjoint_copies_of_OA1.append([x2*g1+x1 for x1,x2 in zip(B1,B2)])
            new_g1_parallel_classes.append(disjoint_copies_of_OA1)

        # Check our stuff before we return it
        if check:
            for classs in new_g1_parallel_classes:
                assert_c_partition(classs, k, g2 * g1, g1)
            for classs in new_parallel_classes:
                assert_c_partition(classs, k, g2 * g1, 1)

        return new_g1_parallel_classes, new_parallel_classes

    # The three factors product construction begins !
    #
    # OA1 and resolvable OA2 and OA3
    OA1 = orthogonal_array(k,n1)
    OA3 = sorted(orthogonal_array(k+1,n3))
    OA3 = [B[1:] for B in OA3]
    OA2 = orthogonal_array(k+1,n2)
    OA2.sort()
    OA2 = [B[1:] for B in OA2]

    # We split OA3 into as many n1-parallel classes as possible, i.e. n3//n1 classes of size n1*n3
    OA3_n1_parall = [OA3[i:i+n1*n3] for i in range(0,(n3-n1)*n3,n1*n3)]

    # Leftover blocks become parallel classes. We must split them into slices of
    # length n3
    OA3_parall = [OA3[i:i+n3] for i in range(len(OA3_n1_parall)*n1*n3, len(OA3), n3)]

    # First product: OA1 and OA3
    n1_parall, parall = product_with_parallel_classes(OA1,k,n1,n3,OA3_n1_parall,OA3_parall,check=check)

    if check:
        OA_13 = [block for classs in parall+n1_parall for block in classs]
        assert is_orthogonal_array(OA_13,k,n1*n3,2,1)

    # Add parallel classes to turn the n1-parall classes into n2-parallel classes
    for classs in n1_parall:
        for i in range(n2-n1):
            classs.extend(parall.pop())

    n2_parall = n1_parall
    del n1_parall

    # We compute the product of OA2 with our decomposition of OA1xOA2 into
    # n2-parallel classes and parallel classes
    n2_parall, parall = product_with_parallel_classes(OA2,k,n2,n1*n3,n2_parall,parall,check=check)
    for n2_classs in n2_parall:
        for i in range(n2):
            partition = [B for j in range(n1*n3) for B in n2_classs[j*n2**2+i*n2:j*n2**2+(i+1)*n2]]
            parall.append(partition)

    # That's what we fought for: this design is resolvable, so let's add a last
    # column to them
    for i,classs in enumerate(parall):
        for B in classs:
            B.append(i)

    OA = [block for classs in parall for block in classs]

    if check:
        assert is_orthogonal_array(OA,k+1,n1*n2*n3,2,1)

    return OA


def _reorder_matrix(matrix):
    r"""
    Return a matrix which is obtained from ``matrix`` by permutation of each row
    in which each column contain every symbol exactly once.

    The input must be a `N \times k` matrix with entries in `\{0,\ldots,N-1\}`
    such that:

    - the symbols on each row are distinct (and hence can be identified with
      subsets of `\{0,\ldots,N-1\}`),
    - each symbol appear exactly `k` times.

    The problem is equivalent to an edge coloring of a bipartite graph. This
    function is used by :func:`brouwer_separable_design`.

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays_build_recursive import _reorder_matrix
        sage: N = 4; k = 3
        sage: M = [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]
        sage: M2 = _reorder_matrix(M)
        sage: all(set(M2[i][0] for i in range(N)) == set(range(N)) for i in range(k))
        True

        sage: M = [list(range(10))] * 10
        sage: N = k = 10
        sage: M2 = _reorder_matrix(M)
        sage: all(set(M2[i][0] for i in range(N)) == set(range(N)) for i in range(k))
        True
    """
    from sage.graphs.graph import Graph

    N = len(matrix)
    k = len(matrix[0])

    g = Graph()
    g.add_edges((x,N+i) for i,S in enumerate(matrix) for x in S)
    matrix = []
    for _ in range(k):
        matching = g.matching(algorithm='LP')
        col = [0]*N
        for x,i,_ in matching:
            if i < N:
                x,i = i,x
            col[i-N] = x
        matrix.append(col)
        g.delete_edges(matching)

    return list(zip(*matrix))


def brouwer_separable_design(k,t,q,x,check=False,verbose=False,explain_construction=False):
    r"""
    Return a `OA(k,t(q^2+q+1)+x)` using Brouwer's result on separable designs.

    This method is an implementation of Brouwer's construction presented in
    [Brouwer80]_. It consists in a systematic application of the usual
    transformation from PBD to OA, applied to a specific PBD.

    **Baer subplanes**

    When `q` is a prime power, the projective plane `PG(2,q^2)` can be
    partitionned into subplanes `PG(2,q)` (called Baer subplanes), giving
    `PG(2,q^2)=B_1\cup \dots\cup B_{q^2-q+1}`. As a result, every line of the
    `PG(2,q^2)` intersects one of the subplane on `q+1` points and all others on
    `1` point.

    The `OA` are built by considering `B_1\cup\dots\cup B_t`, for a total of
    `t(q^2+q+1)` points (to which `x` new points are then added). The blocks of
    this subdesign belong to two categories:

    * The blocks of size `t`: they come from the lines which intersect a
      `B_i` on `q+1` points for some `i>t`. The blocks of size `t` can be partitionned
      into `q^2-q+t-1` parallel classes according to their associated subplane `B_i`
      with `i>t`.

    * The blocks of size `q+t`: those blocks form a symmetric design, as every
      point is incident with `q+t` of them.

    **Constructions**

    In the following, we write `N=t(q^2+q+1)+x`. The code is also heavily
    commented, and will clear any doubt.

    * i) `x=0`: in that case we build a resolvable `OA(k-1,N)` that will then be
      completed into an `OA(k,N)`.

        * *Sets of size* `t`)

          We take the product of each parallel class with the parallel classes
          of a resolvable `OA(k-1,t)-t.OA(k-1,t)`, yielding new parallel
          classes.

        * *Sets of size* `q+t`)

          A `N \times (q+t)` array is built whose rows are the sets of size
          `q+t` such that every value appears once per column. For each block of
          a `OA(k-1,q+t)-(q+t).OA(k-1,t)`, the product with the rows of the
          matrix yields a parallel class.

    * ii) `x=q+t`

        * *Sets of size* `t`)

          Each set of size `t` gives a `OA(k,t)-t.OA(k,1)`, except if there is
          only one parallel class in which case a `OA(k,t)` is sufficient.

        * *Sets of size* `q+t`)

          A `(N-x) \times (q+t)` array `M` is built whose `N-x` rows are the
          sets of size `q+t` such that every value appears once per column. For
          each of the new `x=q+t` points `p_1,\dots,p_{q+t}` we build a matrix
          `M_i` obtained from `M` by adding a column equal to `(p_i,p_i,p_i\dots
          )`. We add to the OA the product of all rows of the `M_i` with the
          block of the `x=q+t` parallel classes of a resolvable
          `OA(k,t+q+1)-(t+q+1).OA(k,1)`.

        * *Set of size* `x`) An `OA(k,x)`

    * iii) `x = q^2-q+1-t`

        * *Sets of size* `t`)

          All blocks of the `i`-th parallel class are extended with the `i`-th
          new point. The blocks are then replaced by a `OA(k,t+1)-(t+1).OA(k,1)`
          or, if there is only one parallel class (i.e. `x=1`) by a
          `OA(k,t+1)-OA(k,1)`.

        * *Set of size* `q+t`)

          They are replaced by `OA(k,q+t)-(q+t).OA(k,1)`.

        * *Set of size* `x`) An `OA(k,x)`

    * iv) `x = q^2+1`

        * *Sets of size* `t`)

          All blocks of the `i`-th parallel class are extended with the `i`-th
          new point (the other `x-q-t` new points are not touched at this
          step). The blocks are then replaced by a `OA(k,t+1)-(t+1).OA(k,1)` or,
          if there is only one parallel class (i.e. `x=1`) by a
          `OA(k,t+1)-OA(k,1)`.

        * *Sets of size* `q+t`) Same as for ii)

        * *Set of size* `x`) An `OA(k,x)`

    * v) `0<x<q^2-q+1-t`

        * *Sets of size* `t`)

          The blocks of the first `x` parallel class are extended with the `x`
          new points, and replaced with `OA(k.t+1)-(t+1).OA(k,1)` or, if `x=1`,
          by `OA(k.t+1)-.OA(k,1)`

          The blocks of the other parallel classes are replaced by
          `OA(k,t)-t.OA(k,t)` or, if there is only one class left, by
          `OA(k,t)-OA(k,t)`

        * *Sets of size* `q+t`)

          They are replaced with `OA(k,q+t)-(q+t).OA(k,1)`.

        * *Set of size* `x`) An `OA(k,x)`

    * vi) `t+q<x<q^2+1`

        * *Sets of size* `t`) Same as in v) with an `x` equal to `x-q+t`.

        * *Sets of size* `t`) Same as in vii)

        * *Set of size* `x`) An `OA(k,x)`

    INPUT:

    - ``k``, ``t``, ``q``, ``x`` -- integers

    - ``check`` -- boolean (default: ``False``); whether to check that output
      is correct before returning it

    - ``verbose`` -- boolean; whether to print some information on the
      construction and parameters being used

    - ``explain_construction`` -- boolean; return a string describing
      the construction

    .. SEEALSO::

        - :func:`~sage.combinat.designs.orthogonal_arrays_find_recursive.find_brouwer_separable_design`

    REFERENCES:

    .. [Brouwer80] A Series of Separable Designs with Application to Pairwise Orthogonal Latin Squares,
      Andries E. Brouwer,
      Vol. 1, n. 1, pp. 39-41,
      European Journal of Combinatorics, 1980
      http://www.sciencedirect.com/science/article/pii/S0195669880800199

    EXAMPLES:

    Test all possible cases::

        sage: from sage.combinat.designs.orthogonal_arrays_build_recursive import brouwer_separable_design
        sage: k,q,t=4,4,3; _=brouwer_separable_design(k,q,t,0,verbose=True)
        Case i) with k=4,q=3,t=4,x=0
        sage: k,q,t=3,3,3; _=brouwer_separable_design(k,t,q,t+q,verbose=True,check=True)
        Case ii) with k=3,q=3,t=3,x=6,e3=1
        sage: k,q,t=3,3,6; _=brouwer_separable_design(k,t,q,t+q,verbose=True,check=True)
        Case ii) with k=3,q=3,t=6,x=9,e3=0
        sage: k,q,t=3,3,6; _=brouwer_separable_design(k,t,q,q**2-q+1-t,verbose=True,check=True)
        Case iii) with k=3,q=3,t=6,x=1,e2=0
        sage: k,q,t=3,4,6; _=brouwer_separable_design(k,t,q,q**2-q+1-t,verbose=True,check=True)
        Case iii) with k=3,q=4,t=6,x=7,e2=1
        sage: k,q,t=3,4,6; _=brouwer_separable_design(k,t,q,q**2+1,verbose=True,check=True)
        Case iv) with k=3,q=4,t=6,x=17,e4=1
        sage: k,q,t=3,2,2; _=brouwer_separable_design(k,t,q,q**2+1,verbose=True,check=True)
        Case iv) with k=3,q=2,t=2,x=5,e4=0
        sage: k,q,t=3,4,7; _=brouwer_separable_design(k,t,q,3,verbose=True,check=True)
        Case v) with k=3,q=4,t=7,x=3,e1=1,e2=1
        sage: k,q,t=3,4,7; _=brouwer_separable_design(k,t,q,1,verbose=True,check=True)
        Case v) with k=3,q=4,t=7,x=1,e1=1,e2=0
        sage: k,q,t=3,4,7; _=brouwer_separable_design(k,t,q,q**2-q-t,verbose=True,check=True)
        Case v) with k=3,q=4,t=7,x=5,e1=0,e2=1
        sage: k,q,t=5,4,7; _=brouwer_separable_design(k,t,q,t+q+3,verbose=True,check=True)
        Case vi) with k=5,q=4,t=7,x=14,e3=1,e4=1
        sage: k,q,t=5,4,8; _=brouwer_separable_design(k,t,q,t+q+1,verbose=True,check=True)
        Case vi) with k=5,q=4,t=8,x=13,e3=1,e4=0
        sage: k,q,t=5,4,8; _=brouwer_separable_design(k,t,q,q**2,verbose=True,check=True)
        Case vi) with k=5,q=4,t=8,x=16,e3=0,e4=1

        sage: print(designs.orthogonal_arrays.explain_construction(10,189))
        Brouwer's separable design construction with t=9,q=4,x=0 from:
           Andries E. Brouwer,
           A series of separable designs with application to pairwise orthogonal Latin squares
           Vol. 1, n. 1, pp. 39-41,
           European Journal of Combinatorics, 1980
    """
    from sage.combinat.designs.orthogonal_arrays import OA_from_PBD
    from .difference_family import difference_family
    from .orthogonal_arrays import incomplete_orthogonal_array
    from sage.arith.misc import is_prime_power

    if explain_construction:
        return ("Brouwer's separable design construction with t={},q={},x={} from:\n" +
                "    Andries E. Brouwer,\n" +
                "    A series of separable designs with application to pairwise orthogonal Latin squares\n" +
                "    Vol. 1, n. 1, pp. 39-41,\n" +
                "    European Journal of Combinatorics, 1980").format(t,q,x)

    ###########################################################
    # Part 1: compute the separable PBD on t(q^2+q+1) points. #
    ###########################################################

    assert t < q**2-q+1
    assert x >= 0
    assert is_prime_power(q)
    N2 = q**4+q**2+1
    N1 = q**2 + q + 1

    # A projective plane on (q^2-q+1)*(q^2+q+1)=q^4+q^2+1 points
    B = difference_family(N2,q**2+1,1)[1][0]
    BIBD = [[(xx+i) % N2 for xx in B] for i in range(N2)]

    # Each congruence class mod q^2-q+1 yields a Baer subplane. Let's check that:
    m = q**2-q+1
    for i in range(m):
        for B in BIBD:
            assert sum((xx % m) == i for xx in B) in [1,q+1], sum((xx % m) == i for xx in B)

    # We are only interested by the points of the first t Baer subplanes (each
    # has size q**2+q+1). Note that each block of the projective plane:
    #
    # - Intersects one Baer plane on q+1 points.
    # - Intersects all other Baer planes on 1 point.
    #
    # When the design its truncated to its first t Baer subplanes, all blocks
    # now have size t or t+q, and cover t(q^2+q+1) points.
    #
    # 1) The blocks of size t can be partitionned into q**2-q+1-t parallel
    #    classes, according to the Baer plane in which they contain q+1
    #    elements.
    #
    # 2) The blocks of size q+t are a symmetric design

    blocks_of_size_q_plus_t = []
    partition_of_blocks_of_size_t = [[] for _ in repeat(None, m - t)]

    relabel = {i+j*m: N1*i+j for i in range(t) for j in range(N1)}

    for B in BIBD:
        # Find the Baer subplane which B intersects on more than 1 point
        B_mod = sorted(xx % m for xx in B)
        while B_mod.pop(0) != B_mod[0]:
            pass
        plane = B_mod[0]
        if plane < t:
            blocks_of_size_q_plus_t.append([relabel[xx] for xx in B if xx % m < t])
        else:
            partition_of_blocks_of_size_t[plane-t].append([relabel[xx] for xx in B if xx % m < t])

    ###########################################################################
    # Separable design built !
    # ------------------------
    #
    # At this point we have a PBD on t*(q**2+q+1) points. Its blocks are
    # split into:
    #
    # - partition_of_blocks_of_size_t : contains all blocks of size t split into
    #                                   q^2-q+t-1 parallel classes.
    #
    # - blocks_of_size_q_plus_t : contains all t*(q**2+q+1)blocks of size q+t,
    #                             covering the same number of points: it is a
    #                             symmetric design.
    ###########################################################################

    ##############################################
    # Part 2: Build an OA on t(q^2+q+1)+x points #
    ##############################################

    e1 = int(x != q**2-q-t)
    e2 = int(x != 1)
    e3 = int(x != q**2)
    e4 = int(x != t+q+1)
    N = t*N1+x

    # i)
    if x == 0:

        if verbose:
            print("Case i) with k={},q={},t={},x={}".format(k, q, t, x))

        # 1) We build a resolvable OA(k-1,t)-t.OA(k-1,1).
        #    With it, from every parallel class with blocks of size t we build a
        #    parallel class of a resolvable OA(k-1,N)

        rOA_N_classes = []

        # A resolvable OA(k-1,t)-t.OA(k-1,1)
        OA_t = incomplete_orthogonal_array(k-1,t,[1]*t,resolvable=True)
        OA_t_classes = [OA_t[i*t:(i+1)*t] for i in range(t-1)]

        # We can now build (t-1)(q^2-q+1-t) parallel classes of the resolvable
        # OA(k-1,N)
        for PBD_parallel_class in partition_of_blocks_of_size_t:
            for OA_class in OA_t_classes:
                rOA_N_classes.append([[B[x] for x in BB]
                                      for BB in OA_class
                                      for B in PBD_parallel_class])

        # 2) We build a Nx(q+t) matrix such that:
        #
        #    a) Each row is a set of size q+t of the PBD
        #    b) an element appears exactly once per column.
        #
        # (This is equivalent to an edge coloring of the (bipartite) incidence
        # graph of points and sets)

        block_of_size_q_plus_t = _reorder_matrix(blocks_of_size_q_plus_t)

        # 3) We now create blocks of an OA(k-1,N) as the product of
        #    a) A set of size q+t (i.e. a row of the matrix)
        #    b) An OA(k-1,q+t)-(q+t).OA(k-1,1)
        #
        # Thanks to the ordering of the points in each set of size q+t, the
        # product of a block B of the incomplete OA with all blocks of size q+t
        # yields a parallel class of an OA(k-1,N)
        OA = incomplete_orthogonal_array(k-1,q+t,[1]*(q+t))
        for B in OA:
            rOA_N_classes.append([[R[x] for x in B] for R in block_of_size_q_plus_t])

        # 4) A last parallel class with blocks [0,0,...], [1,1,...],...
        rOA_N_classes.append([[i]*(k-1) for i in range(N)])

        # 5) We now build the OA(k,N) from the N parallel classes of our resolvable OA(k-1,N)
        OA = [B for classs in rOA_N_classes for B in classs]
        for i,B in enumerate(OA):
            B.append(i//N)

    # ii)
    elif (x == t+q and
          orthogonal_array(k+e3,  t  ,existence=True) and
          orthogonal_array( k  , t+q ,existence=True) and
          orthogonal_array( k+1,t+q+1,existence=True)):

        if verbose:
            print("Case ii) with k={},q={},t={},x={},e3={}".format(k,q,t,x,e3))

        # The sets of size t:
        #
        # This is the usual OA_from_PBD replacement. If there is only one class
        # an OA(k,t) can be used instead of an OA(k+1,t)

        if x == q**2:
            assert e3 == 0, "equivalent to x==q^2"
            assert len(partition_of_blocks_of_size_t) == 1, "also equivalent to exactly one partition into sets of size t"
            OA = [[B[xx] for xx in R] for R in orthogonal_array(k,t) for B in partition_of_blocks_of_size_t[0]]
        else:
            OA = OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t,[]),check=False)[:-N]
            OA.extend([i]*k for i in range(N-x))

        # The sets of size q+t:
        #
        # We build an OA(k,t+q+1)-(t+q+1).OA(k,t+q+1) and the reordered
        # matrix. We then compute the product of every parallel class of the OA
        # (x classes in total) with the rows of the ordered matrix (extended
        # with one of the new x points).

        # Resolvable OA(k,t+q+1)-(t+q+1).OA(k,t+q+1)
        OA_tq1 = incomplete_orthogonal_array(k,t+q+1,[1]*(t+q+1),resolvable=True)
        OA_tq1_classes = [OA_tq1[i*(t+q+1):(i+1)*(t+q+1)] for i in range(t+q)]

        blocks_of_size_q_plus_t = _reorder_matrix(blocks_of_size_q_plus_t)

        for i,classs in enumerate(OA_tq1_classes):
            OA.extend([R[xx] if xx < t+q else N-i-1 for xx in B] for R in blocks_of_size_q_plus_t for B in classs)

        # The set of size x
        OA.extend([N-1-xx for xx in R] for R in orthogonal_array(k,x))

    # iii)
    elif (x == q**2-q+1-t and
          orthogonal_array( k  ,  x  ,existence=True) and # d0
          orthogonal_array(k+e2, t+1 ,existence=True) and # d2-e2
          orthogonal_array(k+1 , t+q ,existence=True)):   # d3-e1
        if verbose:
            print("Case iii) with k={},q={},t={},x={},e2={}".format(k,q,t,x,e2))

        OA = []

        # Each of the x partition into blocks of size t is extended with one of
        # the new x points.

        if x == 1:
            assert e2 == 0, "equivalent to x=1"
            # There is one partition into blocks of size t, which we extend with
            # the new vertex. The OA on t+1 points does not have to be resolvable.

            OA.extend([B[xx] if xx < t else N-1 for xx in R]
                       for R in incomplete_orthogonal_array(k,t+1,[1])
                       for B in partition_of_blocks_of_size_t[0])

        else:
            assert e2 == 1, "equivalent to x!=1"
            # Extending the x partitions into blocks of size t with
            # each of the new x points.

            for i,partition in enumerate(partition_of_blocks_of_size_t):
                for B in partition:
                    B.append(N-i-1)
            OA = OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t,[]),check=False)[:-x]

        # The blocks of size q+t are covered with a resolvable OA(k,q+t)
        OA.extend(OA_from_PBD(k,N,blocks_of_size_q_plus_t,check=False)[:-N])

        # The set of size x
        OA.extend([N-xx-1 for xx in B] for B in orthogonal_array(k,x))

    # iv)
    elif (x == q**2 + 1 and
          orthogonal_array(k, x, existence=True) and            # d0
          orthogonal_array(k + e4, t + 1, existence=True) and   # d2 - e4
          orthogonal_array(k + 1, t + q + 1, existence=True)):  # d4 - 1

        if verbose:
            print(f"Case iv) with k={k},q={q},t={t},x={x},e4={e4}")

        # Sets of size t:
        #
        # All partitions of t-sets are extended with as many new points

        if e4 == 0:
            # Only one partition into t-sets. The OA(k,t+1) needs not be resolvable
            OA = [[B[xx] if xx < t else N-x for xx in R]
                  for R in incomplete_orthogonal_array(k,t+1,[1])
                  for B in partition_of_blocks_of_size_t[0]]
        else:
            for i,classs in enumerate(partition_of_blocks_of_size_t):
                for B in classs:
                    B.append(N-x+i)
            OA = OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t,[]),check=False)[:-x]

        # The sets of size q+t:
        #
        # We build an OA(k,t+q+1)-(t+q+1).OA(k,t+q+1) and the reordered
        # matrix. We then compute the product of every parallel class of the OA
        # (q+t classes in total) with the rows of the ordered matrix (extended
        # with the last q+t new points).

        # Resolvable OA(k,t+q+1)-(t+q+1).OA(k,t+q+1)
        OA_tq1 = incomplete_orthogonal_array(k,t+q+1,[1]*(t+q+1),resolvable=True)
        OA_tq1_classes = [OA_tq1[i*(t+q+1):(i+1)*(t+q+1)] for i in range(t+q)]

        blocks_of_size_q_plus_t = _reorder_matrix(blocks_of_size_q_plus_t)

        for i,classs in enumerate(OA_tq1_classes):
            OA.extend([R[xx] if xx < t+q else N-i-1 for xx in B] for R in blocks_of_size_q_plus_t for B in classs)

        # Set of size x
        OA_k_x = orthogonal_array(k,x)
        OA.extend([N-i-1 for i in R] for R in OA_k_x)

    # v)
    elif (0 < x and x < q**2-q+1-t and (e1 or e2) and # The result is wrong when e1=e2=0
          orthogonal_array(k   ,x  ,existence=True) and # d0
          orthogonal_array(k+e1,t  ,existence=True) and # d1-e1
          orthogonal_array(k+e2,t+1,existence=True) and # d2-e2
          orthogonal_array(k+1,t+q,existence=True)):   # d3-1
        if verbose:
            print("Case v) with k={},q={},t={},x={},e1={},e2={}".format(k,q,t,x,e1,e2))

        OA = []

        # Sets of size t+1
        #
        # We extend x partitions into blocks of size t with the new x elements
        if e2:
            assert x != 1, "equivalent to e2==1"
            for i,classs in enumerate(partition_of_blocks_of_size_t[:x]):
                for B in classs:
                    B.append(N-1-i)
            OA.extend(OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t[:x],[]),check=False)[:-N])

        else:
            assert x == 1, "equivalent to e2==0"
            # Only one class, the OA(k,t+1) need not be resolvable.
            OA.extend([B[xx] if xx < t else N-1 for xx in R]
                       for R in incomplete_orthogonal_array(k,t+1,[1])
                       for B in partition_of_blocks_of_size_t[0])

        # Sets of size t
        if e1:
            assert x != q**2-q-t, "equivalent to e1=1"
            OA.extend(OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t[x:],[]),check=False)[:-N])
        else:
            assert x == q**2-q-t, "equivalent to e1=0"
            # Only one class. The OA(k,t) needs not be resolvable
            OA.extend([B[xx] for xx in R] for R in orthogonal_array(k,t) for B in partition_of_blocks_of_size_t[-1])

        if e1 and e2:
            OA.extend([i]*k for i in range(N-x))

        if e1 == 0 and e2 == 0:
            raise RuntimeError("Brouwer's construction does not work for case v) with e2=e1=0")

        # Sets of size q+t
        OA.extend(OA_from_PBD(k,N,blocks_of_size_q_plus_t,check=False)[:-N])

        # Set of size x
        OA.extend([N-i-1 for i in R] for R in orthogonal_array(k,x))

    # vi)
    elif (t+q < x and x < q**2+1 and (e3 or e4) and # The result is wrong when e3=e4=0
          orthogonal_array(k   ,x   ,existence=True) and # d0
          orthogonal_array(k+e3,t   ,existence=True) and # d1-e3
          orthogonal_array(k+e4,t+1 ,existence=True) and # d2-e4
          orthogonal_array(k+1,t+q+1,existence=True)):    # d4-1
        if verbose:
            print("Case vi) with k={},q={},t={},x={},e3={},e4={}".format(k,q,t,x,e3,e4))

        OA = []

        # Sets of size t+1
        #
        # All x-(q+t) parallel classes with blocks of size t are extended with
        # x-(q+t) of the new points.
        if e4:
            assert x != q+t+1, "equivalent to e4=1"
            for i,classs in enumerate(partition_of_blocks_of_size_t[:x-(q+t)]):
                for B in classs:
                    B.append(N-x+i)
            OA.extend(OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t[:x-(q+t)],[]),check=False)[:-N])
        else:
            assert x == q+t+1, "equivalent to e4=0"
            # Only one class. The OA(k,t+1) needs not be resolvable.
            OA.extend([B[xx] if xx < t else N-x for xx in R]
                       for R in incomplete_orthogonal_array(k,t+1,[1])
                       for B in partition_of_blocks_of_size_t[0])

        # Sets of size t
        if e3:
            assert x != q**2, "equivalent to e3=1"
            OA.extend(OA_from_PBD(k,N,sum(partition_of_blocks_of_size_t[x-(q+t):],[]),check=False)[:-N])
        else:
            assert x == q**2, "equivalent to e3=0"
            # Only one class. The OA(k,t) needs not be resolvable.
            OA.extend([B[xx] for xx in R]
                       for R in orthogonal_array(k,t)
                       for B in partition_of_blocks_of_size_t[-1])

        if e3 and e4:
            OA.extend([i]*k for i in range(N-x))

        elif e3 == 0 and e4 == 0:
            raise RuntimeError("Brouwer's construction does not work for case v) with e3=e4=0")

        # The sets of size q+t:
        #
        # We build an OA(k,t+q+1)-(t+q+1).OA(k,t+q+1) and the reordered
        # matrix. We then compute the product of every parallel class of the OA
        # (q+t classes in total) with the rows of the ordered matrix (extended
        # with the last q+t new points).

        # Resolvable OA(k,t+q+1)-(t+q+1).OA(k,t+q+1)
        OA_tq1 = incomplete_orthogonal_array(k,t+q+1,[1]*(t+q+1),resolvable=True)
        OA_tq1_classes = [OA_tq1[i*(t+q+1):(i+1)*(t+q+1)] for i in range(t+q)]

        blocks_of_size_q_plus_t = _reorder_matrix(blocks_of_size_q_plus_t)

        for i,classs in enumerate(OA_tq1_classes):
            OA.extend([R[xx] if xx < t+q else N-i-1 for xx in B]
                       for R in blocks_of_size_q_plus_t
                       for B in classs)

        # Set of size x
        OA.extend([N-xx-1 for xx in B] for B in orthogonal_array(k,x))

    else:
        raise ValueError("this input is not handled by Brouwer's result")

    if check:
        assert is_orthogonal_array(OA,k,N,2,1)
    return OA
