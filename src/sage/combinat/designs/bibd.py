r"""
Balanced Incomplete Block Designs (BIBD)

This module gathers everything related to Balanced Incomplete Block Designs. One can build a
BIBD (or check that it can be built) with :func:`balanced_incomplete_block_design`::

    sage: BIBD = designs.balanced_incomplete_block_design(7,3,1)                        # needs sage.schemes

In particular, Sage can build a `(v,k,1)`-BIBD when one exists for all `k\leq
5`. The following functions are available:


.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`balanced_incomplete_block_design` | Return a BIBD of parameters `v,k,\lambda`.
    :func:`BIBD_from_TD` | Return a BIBD through TD-based constructions.
    :func:`BIBD_from_difference_family` | Return the BIBD associated to the difference family ``D`` on the group ``G``.
    :func:`BIBD_from_PBD` | Return a `(v,k,1)`-BIBD from a `(r,K)`-PBD where `r=(v-1)/(k-1)`.
    :func:`PBD_from_TD` | Return a `(kt,\{k,t\})`-PBD if `u=0` and a `(kt+u,\{k,k+1,t,u\})`-PBD otherwise.
    :func:`steiner_triple_system` | Return a Steiner Triple System.
    :func:`v_5_1_BIBD` | Return a `(v,5,1)`-BIBD.
    :func:`v_4_1_BIBD` | Return a `(v,4,1)`-BIBD.
    :func:`PBD_4_5_8_9_12` | Return a `(v,\{4,5,8,9,12\})`-PBD on `v` elements.
    :func:`BIBD_5q_5_for_q_prime_power` | Return a `(5q,5,1)`-BIBD with `q\equiv 1\pmod 4` a prime power.


**Construction of BIBD when** `k=4`

Decompositions of `K_v` into `K_4` (i.e. `(v,4,1)`-BIBD) are built following
Douglas Stinson's construction as presented in [Stinson2004]_ page 167. It is
based upon the construction of `(v\{4,5,8,9,12\})`-PBD (see the doc of
:func:`PBD_4_5_8_9_12`), knowing that a `(v\{4,5,8,9,12\})`-PBD on `v` points
can always be transformed into a `((k-1)v+1,4,1)`-BIBD, which covers all
possible cases of `(v,4,1)`-BIBD.

**Construction of BIBD when** `k=5`

Decompositions of `K_v` into `K_4` (i.e. `(v,4,1)`-BIBD) are built following
Clayton Smith's construction [ClaytonSmith]_.

.. [ClaytonSmith] On the existence of `(v,5,1)`-BIBD.
  http://www.argilo.net/files/bibd.pdf
  Clayton Smith


Functions
---------
"""

from sage.arith.misc import binomial, is_prime_power, is_square
from sage.categories.sets_cat import EmptySetError
from sage.misc.lazy_import import lazy_import
from sage.misc.unknown import Unknown

from .design_catalog import transversal_design  # type:ignore
from .designs_pyx import is_pairwise_balanced_design
from .group_divisible_designs import GroupDivisibleDesign

lazy_import('sage.schemes.plane_conics.constructor', 'Conic')


def biplane(n, existence=False):
    r"""
    Return a biplane of order `n`.

    A biplane of order `n` is a symmetric `(1+\frac {(n+1)(n+2)} {2}, n+2, 2)`-BIBD.
    A symmetric (or square) `(v,k,\lambda)`-BIBD is a `(v,k,\lambda)`-BIBD with `v` blocks.

    INPUT:

    - ``n`` -- integer; order of the biplane

    - ``existence`` -- boolean; instead of building the design, return:

      - ``True`` -- meaning that Sage knows how to build the design

      - ``Unknown`` -- meaning that Sage does not know how to build the
        design, but that the design may exist (see :mod:`sage.misc.unknown`)

      - ``False`` -- meaning that the design does not exist

    .. SEEALSO::

        * :func:`balanced_incomplete_block_design`

    EXAMPLES::

        sage: designs.biplane(4)                                                        # needs sage.rings.finite_rings
        (16,6,2)-Balanced Incomplete Block Design
        sage: designs.biplane(7, existence=True)                                        # needs sage.schemes
        True
        sage: designs.biplane(11)                                                       # needs sage.schemes
        (79,13,2)-Balanced Incomplete Block Design

    TESTS::

        sage: designs.biplane(9)                                                        # needs sage.libs.gap
        (56,11,2)-Balanced Incomplete Block Design

    Check all known biplanes::

        sage: [n for n in [0,1,2,3,4,7,9,11]                                            # needs sage.schemes
        ....:  if designs.biplane(n, existence=True) is True]
        [0, 1, 2, 3, 4, 7, 9, 11]
    """
    k = n+2
    v = (k*(k-1))//2 + 1
    return balanced_incomplete_block_design(v, k, lambd=2, existence=existence)


def balanced_incomplete_block_design(v, k, lambd=1, existence=False, use_LJCR=False):
    r"""
    Return a BIBD of parameters `v,k, \lambda`.

    A Balanced Incomplete Block Design of parameters `v,k,\lambda` is a collection
    `\mathcal C` of `k`-subsets of `V=\{0,\dots,v-1\}` such that for any two
    distinct elements `x,y\in V` there are `\lambda` elements `S\in \mathcal C`
    such that `x,y\in S`.

    For more information on BIBD, see the
    :wikipedia:`corresponding Wikipedia entry <Block_design#Definition_of_a_BIBD_.28or_2-design.29>`.

    INPUT:

    - ``v``, ``k``, ``lambd`` -- integers

    - ``existence`` -- boolean; instead of building the design, return:

        - ``True`` -- meaning that Sage knows how to build the design

        - ``Unknown`` -- meaning that Sage does not know how to build the
          design, but that the design may exist (see :mod:`sage.misc.unknown`)

        - ``False`` -- meaning that the design does not exist

    - ``use_LJCR`` -- boolean; whether to query the La Jolla Covering
      Repository for the design when Sage does not know how to build it (see
      :func:`~sage.combinat.designs.covering_design.best_known_covering_design_www`). This
      requires internet.

    .. SEEALSO::

        * :func:`steiner_triple_system`
        * :func:`v_4_1_BIBD`
        * :func:`v_5_1_BIBD`

    .. TODO::

        Implement other constructions from the Handbook of Combinatorial
        Designs.

    EXAMPLES::

        sage: designs.balanced_incomplete_block_design(7, 3, 1).blocks()                # needs sage.schemes
        [[0, 1, 3], [0, 2, 4], [0, 5, 6], [1, 2, 6], [1, 4, 5], [2, 3, 5], [3, 4, 6]]
        sage: B = designs.balanced_incomplete_block_design(66, 6, 1,         # optional - internet
        ....:                                              use_LJCR=True)
        sage: B                                                              # optional - internet
        (66,6,1)-Balanced Incomplete Block Design
        sage: B.blocks()                                                     # optional - internet
        [[0, 1, 2, 3, 4, 65], [0, 5, 22, 32, 38, 58], [0, 6, 21, 30, 43, 48], ...
        sage: designs.balanced_incomplete_block_design(216, 6, 1)
        Traceback (most recent call last):
        ...
        NotImplementedError: I don't know how to build a (216,6,1)-BIBD!

    TESTS::

        sage: designs.balanced_incomplete_block_design(85,5,existence=True)
        True
        sage: _ = designs.balanced_incomplete_block_design(85,5)                        # needs sage.libs.pari

    A BIBD from a Finite Projective Plane::

        sage: _ = designs.balanced_incomplete_block_design(21,5)                        # needs sage.schemes

    Some trivial BIBD::

        sage: designs.balanced_incomplete_block_design(10,10)
        (10,10,1)-Balanced Incomplete Block Design
        sage: designs.balanced_incomplete_block_design(1,10)
        (1,0,1)-Balanced Incomplete Block Design

    Existence of BIBD with `k=3,4,5`::

        sage: [v for v in range(50) if designs.balanced_incomplete_block_design(v,3,existence=True)]                    # needs sage.schemes
        [1, 3, 7, 9, 13, 15, 19, 21, 25, 27, 31, 33, 37, 39, 43, 45, 49]
        sage: [v for v in range(100) if designs.balanced_incomplete_block_design(v,4,existence=True)]                   # needs sage.schemes
        [1, 4, 13, 16, 25, 28, 37, 40, 49, 52, 61, 64, 73, 76, 85, 88, 97]
        sage: [v for v in range(150) if designs.balanced_incomplete_block_design(v,5,existence=True)]                   # needs sage.schemes
        [1, 5, 21, 25, 41, 45, 61, 65, 81, 85, 101, 105, 121, 125, 141, 145]

    For `k > 5` there are currently very few constructions::

        sage: [v for v in range(300) if designs.balanced_incomplete_block_design(v,6,existence=True) is True]           # needs sage.schemes
        [1, 6, 31, 66, 76, 91, 96, 106, 111, 121, 126, 136, 141, 151, 156, 171, 181, 186, 196, 201, 211, 241, 271]
        sage: [v for v in range(300) if designs.balanced_incomplete_block_design(v,6,existence=True) is Unknown]        # needs sage.schemes
        [51, 61, 81, 166, 216, 226, 231, 246, 256, 261, 276, 286, 291]

    Here are some constructions with `k \geq 7` and `v` a prime power::

        sage: # needs sage.libs.pari
        sage: designs.balanced_incomplete_block_design(169,7)
        (169,7,1)-Balanced Incomplete Block Design
        sage: designs.balanced_incomplete_block_design(617,8)
        (617,8,1)-Balanced Incomplete Block Design
        sage: designs.balanced_incomplete_block_design(433,9)
        (433,9,1)-Balanced Incomplete Block Design
        sage: designs.balanced_incomplete_block_design(1171,10)
        (1171,10,1)-Balanced Incomplete Block Design

    And we know some nonexistence results::

        sage: designs.balanced_incomplete_block_design(21,6,existence=True)
        False

    Some BIBDs with `\lambda \neq 1`::

        sage: designs.balanced_incomplete_block_design(176, 50, 14, existence=True)
        True
        sage: designs.balanced_incomplete_block_design(64,28,12)                        # needs sage.libs.pari
        (64,28,12)-Balanced Incomplete Block Design
        sage: designs.balanced_incomplete_block_design(37,9,8)                          # needs sage.libs.pari
        (37,9,8)-Balanced Incomplete Block Design
        sage: designs.balanced_incomplete_block_design(15,7,3)                          # needs sage.schemes
        (15,7,3)-Balanced Incomplete Block Design

    Some BIBDs from the recursive construction ::

        sage: designs.balanced_incomplete_block_design(76,16,4)                         # needs sage.libs.pari
        (76,16,4)-Balanced Incomplete Block Design
        sage: designs.balanced_incomplete_block_design(10,4,2)                          # needs sage.libs.pari
        (10,4,2)-Balanced Incomplete Block Design
        sage: designs.balanced_incomplete_block_design(50,25,24)                        # needs sage.schemes
        (50,25,24)-Balanced Incomplete Block Design
        sage: designs.balanced_incomplete_block_design(29,15,15)                        # needs sage.libs.pari sage.schemes
        (29,15,15)-Balanced Incomplete Block Design
    """
    # Trivial BIBD
    if v == 1:
        if existence:
            return True
        return BIBD(v, [], check=False)

    if k == v:
        if existence:
            return True
        return BIBD(v, [list(range(v)) for _ in range(lambd)],lambd=lambd, check=False, copy=False)

    # Non-existence of BIBD
    if (v < k or
        k < 2 or
        (lambd*(v-1)) % (k-1) != 0 or
        (lambd*v*(v-1)) % (k*(k-1)) != 0 or
        # From the Handbook of combinatorial designs:
        #
        # With lambda>1 other exceptions are
        # (15,5,2),(21,6,2),(22,7,2),(22,8,4).
        (k == 6 and v in [36,46]) or
        (k == 7 and v == 43) or
        # Fisher's inequality
        (lambd*v*(v-1))/(k*(k-1)) < v):
        if existence:
            return False
        raise EmptySetError("There exists no ({},{},{})-BIBD".format(v, k, lambd))

    # Non-existence by BRC Theorem
    if BruckRyserChowla_check(v, k, lambd) is False:
        if existence:
            return False
        raise EmptySetError("There exists no ({},{},{})-BIBD by Bruck-Ryser-Chowla Theorem".format(v,k,lambd))

    if k == 2:
        if existence:
            return True
        return BIBD(v, [[x, y] for _ in range(lambd) for x in range(v) for y in range(x+1, v) if x != y], lambd=lambd, check=False, copy=True)
    if k == 3 and lambd == 1:
        if existence:
            return v % 6 == 1 or v % 6 == 3
        return steiner_triple_system(v)
    if k == 4 and lambd == 1:
        if existence:
            return v % 12 == 1 or v % 12 == 4
        return BIBD(v, v_4_1_BIBD(v), copy=False)
    if k == 5 and lambd == 1:
        if existence:
            return v % 20 == 1 or v % 20 == 5
        return BIBD(v, v_5_1_BIBD(v), copy=False)

    from .difference_family import difference_family
    from .database import BIBD_constructions

    if (v, k, lambd) in BIBD_constructions:
        if existence:
            return True
        return BIBD(v,BIBD_constructions[(v, k, lambd)](), lambd=lambd, copy=False)
    if lambd == 1 and BIBD_from_arc_in_desarguesian_projective_plane(v, k, existence=True):
        if existence:
            return True
        B = BIBD_from_arc_in_desarguesian_projective_plane(v, k)
        return BIBD(v, B, copy=False)
    if lambd == 1 and BIBD_from_TD(v, k, existence=True) is True:
        if existence:
            return True
        return BIBD(v, BIBD_from_TD(v, k), copy=False)
    if lambd == 1 and v == (k-1)**2+k and is_prime_power(k-1):
        if existence:
            return True
        from .block_design import projective_plane
        return BIBD(v, projective_plane(k-1),copy=False)
    if difference_family(v, k, l=lambd, existence=True) is True:
        if existence:
            return True
        G, D = difference_family(v, k, l=lambd)
        return BIBD(v, BIBD_from_difference_family(G, D, check=False), lambd=lambd, copy=False)
    if lambd == 1 and use_LJCR:
        from .covering_design import best_known_covering_design_www
        values_in_db = False
        try:
            B = best_known_covering_design_www(v, k, 2)
            values_in_db = True
        except ValueError:
            # the parameters are not in the LJCR database
            pass

        if values_in_db:
            # Is it a BIBD or just a good covering?
            expected_n_of_blocks = binomial(v, 2) // binomial(k, 2)
            if B.low_bd() > expected_n_of_blocks:
                if existence:
                    return False
                raise EmptySetError(f"there exists no ({v},{k},{lambd})-BIBD")
            B = B.incidence_structure()
            if B.num_blocks() == expected_n_of_blocks:
                if existence:
                    return True
                else:
                    return BIBD(B.ground_set(), B.blocks(), k=k, lambd=1, copy=False)

    if ( (k+lambd)*(k+lambd-1) == lambd*(v+k+lambd-1) and
         balanced_incomplete_block_design(v+k+lambd, k+lambd, lambd, existence=True) is True):
        # By removing a block and all points of that block from the
        # symmetric (v+k+lambd, k+lambd, lambd) BIBD
        # we get a (v, k, lambd) BIBD
        if existence:
            return True

        D = balanced_incomplete_block_design(v+k+lambd, k+lambd, lambd)
        Br = D.blocks()[0]  # block to remove
        blocks = D.blocks()[1:]

        blocks = [set(B).difference(Br) for B in blocks]
        points = set(D.ground_set()).difference(Br)

        return BalancedIncompleteBlockDesign(points, blocks, k=k, lambd=lambd, copy=False)

    if existence:
        return Unknown
    else:
        raise NotImplementedError("I don't know how to build a ({},{},{})-BIBD!".format(v, k, lambd))

def BruckRyserChowla_check(v, k, lambd):
    r"""
    Check whether the parameters passed satisfy the Bruck-Ryser-Chowla theorem.

    For more information on the theorem, see the
    :wikipedia:`corresponding Wikipedia entry <Bruck–Ryser–Chowla_theorem>`.

    INPUT:

    - ``v``, ``k``, ``lambd`` -- integers; parameters to check

    OUTPUT: ``True`` -- the parameters satisfy the theorem

    - ``False`` -- the theorem fails for the given parameters

    - ``Unknown`` -- the preconditions of the theorem are not met

    EXAMPLES:

        sage: from sage.combinat.designs.bibd import BruckRyserChowla_check
        sage: BruckRyserChowla_check(22,7,2)
        False

    Nonexistence of projective planes of order 6 and 14

        sage: from sage.combinat.designs.bibd import BruckRyserChowla_check
        sage: BruckRyserChowla_check(43,7,1)                                            # needs sage.schemes
        False
        sage: BruckRyserChowla_check(211,15,1)                                          # needs sage.schemes
        False

    Existence of symmetric BIBDs with parameters `(79,13,2)` and `(56,11,2)`

        sage: from sage.combinat.designs.bibd import BruckRyserChowla_check
        sage: BruckRyserChowla_check(79,13,2)                                           # needs sage.schemes
        True
        sage: BruckRyserChowla_check(56,11,2)
        True

    TESTS:

    Test some non-symmetric parameters::

        sage: from sage.combinat.designs.bibd import BruckRyserChowla_check
        sage: BruckRyserChowla_check(89,11,3)
        Unknown
        sage: BruckRyserChowla_check(25,23,2)
        Unknown

    Clearly wrong parameters satisfying the theorem::

        sage: from sage.combinat.designs.bibd import BruckRyserChowla_check
        sage: BruckRyserChowla_check(13,25,50)                                          # needs sage.schemes
        True
    """
    from sage.rings.rational_field import QQ

    # design is not symmetric
    if k*(k-1) != lambd*(v-1):
        return Unknown

    if v % 2 == 0:
        return is_square(k-lambd)

    g = 1 if v % 4 == 1 else -1
    C = Conic(QQ, [1, lambd - k, -g * lambd])

    (flag, sol) = C.has_rational_point(point=True)

    return flag

def steiner_triple_system(n):
    r"""
    Return a Steiner Triple System.

    A Steiner Triple System (STS) of a set `\{0,...,n-1\}`
    is a family `S` of 3-sets such that for any `i \not = j`
    there exists exactly one set of `S` in which they are
    both contained.

    It can alternatively be thought of as a factorization of
    the complete graph `K_n` with triangles.

    A Steiner Triple System of a `n`-set exists if and only if
    `n \equiv 1 \pmod 6` or `n \equiv 3 \pmod 6`, in which case
    one can be found through Bose's and Skolem's constructions,
    respectively [AndHonk97]_.

    INPUT:

    - ``n`` -- return a Steiner Triple System of `\{0,...,n-1\}`

    EXAMPLES:

    A Steiner Triple System on `9` elements ::

        sage: sts = designs.steiner_triple_system(9)
        sage: sts
        (9,3,1)-Balanced Incomplete Block Design
        sage: list(sts)
        [[0, 1, 5], [0, 2, 4], [0, 3, 6], [0, 7, 8], [1, 2, 3],
         [1, 4, 7], [1, 6, 8], [2, 5, 8], [2, 6, 7], [3, 4, 8],
         [3, 5, 7], [4, 5, 6]]

    As any pair of vertices is covered once, its parameters are ::

        sage: sts.is_t_design(return_parameters=True)
        (True, (2, 9, 3, 1))

    An exception is raised for invalid values of ``n`` ::

        sage: designs.steiner_triple_system(10)
        Traceback (most recent call last):
        ...
        EmptySetError: Steiner triple systems only exist for n = 1 mod 6 or n = 3 mod 6

    REFERENCE:

    .. [AndHonk97] A short course in Combinatorial Designs,
      Ian Anderson, Iiro Honkala,
      Internet Editions, Spring 1997,
      http://www.utu.fi/~honkala/designs.ps
    """

    name = "Steiner Triple System on "+str(n)+" elements"

    if n % 6 == 3:
        t = (n-3) // 6
        Z = list(range(2 * t + 1))

        T = lambda x_y : x_y[0] + (2*t+1)*x_y[1]

        sts = [[(i,0),(i,1),(i,2)] for i in Z] + \
            [[(i,k),(j,k),(((t+1)*(i+j)) % (2*t+1),(k+1) % 3)] for k in range(3) for i in Z for j in Z if i != j]

    elif n % 6 == 1:

        t = (n-1) // 6
        N = list(range(2 * t))
        T = lambda x_y : x_y[0]+x_y[1]*t*2 if x_y != (-1,-1) else n-1

        L1 = lambda i,j : (i+j) % ((n-1)//3)
        L = lambda i,j : L1(i,j)//2 if L1(i,j) % 2 == 0 else t+(L1(i,j)-1)//2

        sts = [[(i,0),(i,1),(i,2)] for i in range(t)] + \
            [[(-1,-1),(i,k),(i-t,(k+1) % 3)] for i in range(t,2*t) for k in [0,1,2]] + \
            [[(i,k),(j,k),(L(i,j),(k+1) % 3)] for k in [0,1,2] for i in N for j in N if i < j]

    else:
        raise EmptySetError("Steiner triple systems only exist for n = 1 mod 6 or n = 3 mod 6")

    # apply T and remove duplicates
    sts = set(frozenset(T(xx) for xx in x) for x in sts)

    return BIBD(n, sts, name=name,check=False)


def BIBD_from_TD(v,k,existence=False):
    r"""
    Return a BIBD through TD-based constructions.

    INPUT:

    - ``v``, ``k`` -- integers; computes a `(v,k,1)`-BIBD

    - ``existence`` -- boolean; instead of building the design, return:

      - ``True`` -- meaning that Sage knows how to build the design

      - ``Unknown`` -- meaning that Sage does not know how to build the
        design, but that the design may exist (see :mod:`sage.misc.unknown`)

      - ``False`` -- meaning that the design does not exist

    This method implements three constructions:

    - If there exists a `TD(k,v)` and a `(v,k,1)`-BIBD then there exists a
      `(kv,k,1)`-BIBD.

      The BIBD is obtained from all blocks of the `TD`, and from the blocks of
      the `(v,k,1)`-BIBDs defined over the `k` groups of the `TD`.

    - If there exists a `TD(k,v)` and a `(v+1,k,1)`-BIBD then there exists a
      `(kv+1,k,1)`-BIBD.

      The BIBD is obtained from all blocks of the `TD`, and from the blocks of
      the `(v+1,k,1)`-BIBDs defined over the sets `V_1\cup \infty,\dots,V_k\cup
      \infty` where the `V_1,\dots,V_k` are the groups of the TD.

    - If there exists a `TD(k,v)` and a `(v+k,k,1)`-BIBD then there exists a
      `(kv+k,k,1)`-BIBD.

      The BIBD is obtained from all blocks of the `TD`, and from the blocks of
      the `(v+k,k,1)`-BIBDs defined over the sets `V_1\cup
      \{\infty_1,\dots,\infty_k\},\dots,V_k\cup \{\infty_1,\dots,\infty_k\}`
      where the `V_1,\dots,V_k` are the groups of the TD. By making sure that
      all copies of the `(v+k,k,1)`-BIBD contain the block
      `\{\infty_1,\dots,\infty_k\}`, the result is also a BIBD.

    These constructions can be found in
    `<http://www.argilo.net/files/bibd.pdf>`_.

    EXAMPLES:

    First construction::

        sage: from sage.combinat.designs.bibd import BIBD_from_TD
        sage: BIBD_from_TD(25,5,existence=True)                                         # needs sage.schemes
        True
        sage: _ = BlockDesign(25,BIBD_from_TD(25,5))                                    # needs sage.schemes

    Second construction::

        sage: from sage.combinat.designs.bibd import BIBD_from_TD
        sage: BIBD_from_TD(21,5,existence=True)                                         # needs sage.schemes
        True
        sage: _ = BlockDesign(21,BIBD_from_TD(21,5))                                    # needs sage.schemes

    Third construction::

        sage: from sage.combinat.designs.bibd import BIBD_from_TD
        sage: BIBD_from_TD(85,5,existence=True)                                         # needs sage.schemes
        True
        sage: _ = BlockDesign(85,BIBD_from_TD(85,5))                                    # needs sage.schemes

    No idea::

        sage: from sage.combinat.designs.bibd import BIBD_from_TD
        sage: BIBD_from_TD(20,5,existence=True)
        Unknown
        sage: BIBD_from_TD(20,5)
        Traceback (most recent call last):
        ...
        NotImplementedError: I do not know how to build a (20,5,1)-BIBD!
    """
    # First construction
    if (v % k == 0 and
        balanced_incomplete_block_design(v//k, k, existence=True) is True and
        transversal_design(k, v//k, existence=True) is True):

        if existence:
            return True

        v = v//k
        BIBDvk = balanced_incomplete_block_design(v,k)._blocks
        TDkv = transversal_design(k,v,check=False)

        BIBD = TDkv._blocks
        for i in range(k):
            BIBD.extend([[x+i*v for x in B] for B in BIBDvk])

    # Second construction
    elif ((v-1) % k == 0 and
        balanced_incomplete_block_design((v-1)//k+1,k,existence=True) is True and
        transversal_design(k,(v-1)//k,existence=True)) is True:

        if existence:
            return True

        v = (v-1)//k
        BIBDv1k = balanced_incomplete_block_design(v+1,k)._blocks
        TDkv = transversal_design(k,v,check=False)._blocks

        inf = v*k
        BIBD = TDkv
        for i in range(k):
            BIBD.extend([[inf if x == v else x+i*v for x in B] for B in BIBDv1k])

    # Third construction
    elif ((v-k) % k == 0 and
        balanced_incomplete_block_design((v-k)//k+k,k,existence=True) is True
        and transversal_design(k,(v-k)//k,existence=True) is True):
        if existence:
            return True

        v = (v-k)//k
        BIBDvpkk = balanced_incomplete_block_design(v+k,k)
        TDkv = transversal_design(k,v,check=False)._blocks
        inf = v*k
        BIBD = TDkv

        # makes sure that [v,...,v+k-1] is a block of BIBDvpkk. Then, we remove it.
        BIBDvpkk = _relabel_bibd(BIBDvpkk,v+k)
        BIBDvpkk = [B for B in BIBDvpkk if min(B) < v]

        for i in range(k):
            BIBD.extend([[(x-v)+inf if x >= v else x+i*v for x in B] for B in BIBDvpkk])

        BIBD.append(list(range(k * v, v * k + k)))

    # No idea ...
    else:
        if existence:
            return Unknown
        else:
            raise NotImplementedError("I do not know how to build a ({},{},1)-BIBD!".format(v,k))

    return BIBD


def BIBD_from_difference_family(G, D, lambd=None, check=True):
    r"""
    Return the BIBD associated to the difference family ``D`` on the group ``G``.

    Let `G` be a group. A `(G,k,\lambda)`-*difference family* is a family `B =
    \{B_1,B_2,\ldots,B_b\}` of `k`-subsets of `G` such that for each element of
    `G \backslash \{0\}` there exists exactly `\lambda` pairs of elements
    `(x,y)`, `x` and `y` belonging to the same block, such that `x - y = g` (or
    x y^{-1} = g` in multiplicative notation).

    If `\{B_1, B_2, \ldots, B_b\}` is a `(G,k,\lambda)`-difference family then
    its set of translates `\{B_i \cdot g; i \in \{1,\ldots,b\}, g \in G\}` is a
    `(v,k,\lambda)`-BIBD where `v` is the cardinality of `G`.

    INPUT:

    - ``G`` -- a finite additive Abelian group

    - ``D`` -- a difference family on ``G`` (short blocks are allowed)

    - ``lambd`` -- the `\lambda` parameter (optional, only used if ``check`` is
      ``True``)

    - ``check`` -- boolean (default: ``True``); whether or not we check the output

    EXAMPLES::

        sage: G = Zmod(21)
        sage: D = [[0,1,4,14,16]]
        sage: sorted(G(x-y) for x in D[0] for y in D[0] if x != y)
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

        sage: from sage.combinat.designs.bibd import BIBD_from_difference_family
        sage: BIBD_from_difference_family(G, D)
        [[0, 1, 4, 14, 16],
         [1, 2, 5, 15, 17],
         [2, 3, 6, 16, 18],
         [3, 4, 7, 17, 19],
         [4, 5, 8, 18, 20],
         [5, 6, 9, 19, 0],
         [6, 7, 10, 20, 1],
         [7, 8, 11, 0, 2],
         [8, 9, 12, 1, 3],
         [9, 10, 13, 2, 4],
         [10, 11, 14, 3, 5],
         [11, 12, 15, 4, 6],
         [12, 13, 16, 5, 7],
         [13, 14, 17, 6, 8],
         [14, 15, 18, 7, 9],
         [15, 16, 19, 8, 10],
         [16, 17, 20, 9, 11],
         [17, 18, 0, 10, 12],
         [18, 19, 1, 11, 13],
         [19, 20, 2, 12, 14],
         [20, 0, 3, 13, 15]]
    """
    from .difference_family import group_law, block_stabilizer
    identity, mul, inv = group_law(G)
    bibd = []
    Gset = set(G)
    p_to_i = {g: i for i, g in enumerate(Gset)}
    for b in D:
        b = [G(w) for w in b]
        S = block_stabilizer(G, b)
        GG = Gset.copy()
        while GG:
            g = GG.pop()
            if S:
                GG.difference_update(mul(s,g) for s in S)
            bibd.append([p_to_i[mul(i,g)] for i in b])

    if check:
        if lambd is None:
            k = len(bibd[0])
            v = G.cardinality()
            lambd = (len(bibd) * k * (k-1)) // (v * (v-1))
        assert is_pairwise_balanced_design(bibd, G.cardinality(), [len(D[0])], lambd=lambd)

    return bibd

################
# (v,4,1)-BIBD #
################

def v_4_1_BIBD(v, check=True):
    r"""
    Return a `(v,4,1)`-BIBD.

    A `(v,4,1)`-BIBD is an edge-decomposition of the complete graph `K_v` into
    copies of `K_4`. For more information, see
    :func:`balanced_incomplete_block_design`. It exists if and only if `v\equiv 1,4
    \pmod {12}`.

    See page 167 of [Stinson2004]_ for the construction details.

    .. SEEALSO::

        * :func:`balanced_incomplete_block_design`

    INPUT:

    - ``v`` -- integer; number of points

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import v_4_1_BIBD  # long time
        sage: for n in range(13,100):                            # long time
        ....:    if n%12 in [1,4]:
        ....:       _ = v_4_1_BIBD(n, check = True)

    TESTS:

    Check that the `(25,4)` and `(37,4)`-difference family are available::

        sage: assert designs.difference_family(25,4,existence=True)
        sage: _ = designs.difference_family(25,4)
        sage: assert designs.difference_family(37,4,existence=True)
        sage: _ = designs.difference_family(37,4)

    Check some larger `(v,4,1)`-BIBD (see :issue:`17557`)::

        sage: for v in range(400):                                      # long time
        ....:     if v%12 in [1,4]:
        ....:         _ = designs.balanced_incomplete_block_design(v,4)
    """
    k = 4
    if v == 0:
        return []
    if v <= 12 or v % 12 not in [1,4]:
        raise EmptySetError("A K_4-decomposition of K_v exists iif v=2,4 mod 12, v>12 or v==0")

    # Step 1. Base cases.
    if v == 13:
        # note: this construction can also be obtained from difference_family
        from .block_design import projective_plane
        return projective_plane(3)._blocks
    if v == 16:
        from .block_design import AffineGeometryDesign
        from sage.rings.finite_rings.finite_field_constructor import FiniteField
        return AffineGeometryDesign(2,1,FiniteField(4,'x'))._blocks
    if v == 25 or v == 37:
        from .difference_family import difference_family
        G,D = difference_family(v,4)
        return BIBD_from_difference_family(G,D,check=False)
    if v == 28:
        return [[0, 1, 23, 26], [0, 2, 10, 11], [0, 3, 16, 18], [0, 4, 15, 20],
                [0, 5, 8, 9], [0, 6, 22, 25], [0, 7, 14, 21], [0, 12, 17, 27],
                [0, 13, 19, 24], [1, 2, 24, 27], [1, 3, 11, 12], [1, 4, 17, 19],
                [1, 5, 14, 16], [1, 6, 9, 10], [1, 7, 20, 25], [1, 8, 15, 22],
                [1, 13, 18, 21], [2, 3, 21, 25], [2, 4, 12, 13], [2, 5, 18, 20],
                [2, 6, 15, 17], [2, 7, 19, 22], [2, 8, 14, 26], [2, 9, 16, 23],
                [3, 4, 22, 26], [3, 5, 7, 13], [3, 6, 14, 19], [3, 8, 20, 23],
                [3, 9, 15, 27], [3, 10, 17, 24], [4, 5, 23, 27], [4, 6, 7, 8],
                [4, 9, 14, 24], [4, 10, 16, 21], [4, 11, 18, 25], [5, 6, 21, 24],
                [5, 10, 15, 25], [5, 11, 17, 22], [5, 12, 19, 26], [6, 11, 16, 26],
                [6, 12, 18, 23], [6, 13, 20, 27], [7, 9, 17, 18], [7, 10, 26, 27],
                [7, 11, 23, 24], [7, 12, 15, 16], [8, 10, 18, 19], [8, 11, 21, 27],
                [8, 12, 24, 25], [8, 13, 16, 17], [9, 11, 19, 20], [9, 12, 21, 22],
                [9, 13, 25, 26], [10, 12, 14, 20], [10, 13, 22, 23], [11, 13, 14, 15],
                [14, 17, 23, 25], [14, 18, 22, 27], [15, 18, 24, 26], [15, 19, 21, 23],
                [16, 19, 25, 27], [16, 20, 22, 24], [17, 20, 21, 26]]

    # Step 2 : this is function PBD_4_5_8_9_12
    PBD = PBD_4_5_8_9_12((v-1)//(k-1),check=False)

    # Step 3 : Theorem 7.20
    bibd = BIBD_from_PBD(PBD,v,k,check=False)

    if check:
        assert is_pairwise_balanced_design(bibd,v,[k])

    return bibd


def BIBD_from_PBD(PBD, v, k, check=True, base_cases=None):
    r"""
    Return a `(v,k,1)`-BIBD from a `(r,K)`-PBD where `r=(v-1)/(k-1)`.

    This is Theorem 7.20 from [Stinson2004]_.

    INPUT:

    - ``v``, ``k`` -- integers

    - ``PBD`` -- a PBD on `r=(v-1)/(k-1)` points, such that for any block of
      ``PBD`` of size `s` there must exist a `((k-1)s+1,k,1)`-BIBD

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.

    - ``base_cases`` -- caching system, for internal use

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import PBD_4_5_8_9_12
        sage: from sage.combinat.designs.bibd import BIBD_from_PBD
        sage: from sage.combinat.designs.bibd import is_pairwise_balanced_design
        sage: PBD = PBD_4_5_8_9_12(17)                                                  # needs sage.schemes
        sage: bibd = is_pairwise_balanced_design(BIBD_from_PBD(PBD,52,4),52,[4])        # needs sage.schemes
    """
    if base_cases is None:
        base_cases = {}
    r = (v-1) // (k-1)
    bibd = []
    for X in PBD:
        n = len(X)
        N = (k-1)*n+1
        if (n,k) not in base_cases:
            base_cases[n,k] = _relabel_bibd(balanced_incomplete_block_design(N,k), N)

        for XX in base_cases[n,k]:
            if N-1 in XX:
                continue
            bibd.append([X[x//(k-1)] + (x % (k-1))*r for x in XX])

    for x in range(r):
        bibd.append([x+i*r for i in range(k-1)]+[v-1])

    if check:
        assert is_pairwise_balanced_design(bibd,v,[k])

    return bibd

def _relabel_bibd(B,n,p=None):
    r"""
    Relabel the BIBD on `n` points and blocks of size k such that
    `\{0,...,k-2,n-1\},\{k-1,...,2k-3,n-1\},...,\{n-k,...,n-2,n-1\}` are blocks
    of the BIBD.

    INPUT:

    - ``B`` -- list of blocks

    - ``n`` -- integer; number of points

    - ``p`` -- (optional) the point that will be labeled with `n-1`

    EXAMPLES::

        sage: designs.balanced_incomplete_block_design(40,4).blocks()  # indirect doctest           # needs sage.schemes
        [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 10],
         [0, 5, 7, 11], [0, 13, 26, 39], [0, 14, 25, 28],
         [0, 15, 27, 38], [0, 16, 22, 32], [0, 17, 23, 34],
        ...
    """
    if p is None:
        p = n-1
    found = 0
    last = n-1
    d = {}
    for X in B:
        if last in X:
            for x in X:
                if x == last:
                    continue
                d[x] = found
                found += 1
            if found == n-1:
                break
    d[p] = n-1
    return [[d[x] for x in X] for X in B]


def PBD_4_5_8_9_12(v, check=True):
    r"""
    Return a `(v,\{4,5,8,9,12\})`-PBD on `v` elements.

    A `(v,\{4,5,8,9,12\})`-PBD exists if and only if `v\equiv 0,1 \pmod 4`. The
    construction implemented here appears page 168 in [Stinson2004]_.

    INPUT:

    - ``v`` -- integer congruent to `0` or `1` modulo `4`

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.

    EXAMPLES::

        sage: designs.balanced_incomplete_block_design(40,4).blocks()  # indirect doctest           # needs sage.schemes
        [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 10],
         [0, 5, 7, 11], [0, 13, 26, 39], [0, 14, 25, 28],
         [0, 15, 27, 38], [0, 16, 22, 32], [0, 17, 23, 34],
        ...

    Check that :issue:`16476` is fixed::

        sage: from sage.combinat.designs.bibd import PBD_4_5_8_9_12
        sage: for v in (0,1,4,5,8,9,12,13,16,17,20,21,24,25):                           # needs sage.schemes
        ....:     _ = PBD_4_5_8_9_12(v)
    """
    if v % 4 not in [0, 1]:
        raise ValueError
    if v <= 1:
        PBD = []
    elif v <= 12:
        PBD = [list(range(v))]
    elif v == 13 or v == 28:
        PBD = v_4_1_BIBD(v, check=False)
    elif v == 29:
        TD47 = transversal_design(4,7)._blocks
        four_more_sets = [[28]+[i*7+j for j in range(7)] for i in range(4)]
        PBD = TD47 + four_more_sets
    elif v == 41:
        TD59 = transversal_design(5,9)
        PBD = ([[x for x in X if x < 41] for X in TD59]
                + [[i*9+j for j in range(9)] for i in range(4)]
                + [[36,37,38,39,40]])
    elif v == 44:
        TD59 = transversal_design(5,9)
        PBD = ([[x for x in X if x < 44] for X in TD59]
                + [[i*9+j for j in range(9)] for i in range(4)]
                + [[36,37,38,39,40,41,42,43]])
    elif v == 45:
        TD59 = transversal_design(5,9)._blocks
        PBD = (TD59+[[i*9+j for j in range(9)] for i in range(5)])
    elif v == 48:
        TD4_12 = transversal_design(4,12)._blocks
        PBD = (TD4_12+[[i*12+j for j in range(12)] for i in range(4)])
    elif v == 49:
        # Lemma 7.16 : A (49,{4,13})-PBD
        TD4_12 = transversal_design(4,12)._blocks

        # Replacing the block of size 13 with a BIBD
        BIBD_13_4 = v_4_1_BIBD(13)
        for i in range(4):
            for B in BIBD_13_4:
                TD4_12.append([i*12+x if x != 12 else 48
                               for x in B])

        PBD = TD4_12
    else:
        t,u = _get_t_u(v)
        TD = transversal_design(5,t)
        TD = [[x for x in X if x < 4*t+u] for X in TD]
        for B in [list(range(t*i,t*(i+1))) for i in range(4)]:
            TD.extend(_PBD_4_5_8_9_12_closure([B]))

        if u > 1:
            TD.extend(_PBD_4_5_8_9_12_closure([list(range(4*t,4*t+u))]))

        PBD = TD

    if check:
        assert is_pairwise_balanced_design(PBD,v,[4,5,8,9,12])

    return PBD

def _PBD_4_5_8_9_12_closure(B):
    r"""
    Makes sure all blocks of `B` have size in `\{4,5,8,9,12\}`.

    This is a helper function for :func:`PBD_4_5_8_9_12`. Given that
    `\{4,5,8,9,12\}` is PBD-closed, any block of size not in `\{4,5,8,9,12\}`
    can be decomposed further.

    EXAMPLES::

        sage: designs.balanced_incomplete_block_design(40,4).blocks()  # indirect doctest           # needs sage.schemes
        [[0, 1, 2, 12], [0, 3, 6, 9], [0, 4, 8, 10],
         [0, 5, 7, 11], [0, 13, 26, 39], [0, 14, 25, 28],
         [0, 15, 27, 38], [0, 16, 22, 32], [0, 17, 23, 34],
        ...
    """
    BB = []
    for X in B:
        if len(X) not in [4,5,8,9,12]:
            PBD = PBD_4_5_8_9_12(len(X), check=False)
            X = [[X[i] for i in XX] for XX in PBD]
            BB.extend(X)
        else:
            BB.append(X)
    return BB


table_7_1 = {
    0:{'t':-4,'u':16,'s':2},
    1:{'t':-4,'u':17,'s':2},
    4:{'t':1,'u':0,'s':1},
    5:{'t':1,'u':1,'s':1},
    8:{'t':1,'u':4,'s':1},
    9:{'t':1,'u':5,'s':1},
    12:{'t':1,'u':8,'s':1},
    13:{'t':1,'u':9,'s':1},
    16:{'t':4,'u':0,'s':0},
    17:{'t':4,'u':1,'s':0},
    20:{'t':5,'u':0,'s':0},
    21:{'t':5,'u':1,'s':0},
    24:{'t':5,'u':4,'s':0},
    25:{'t':5,'u':5,'s':0},
    28:{'t':5,'u':8,'s':1},
    29:{'t':5,'u':9,'s':1},
    32:{'t':8,'u':0,'s':0},
    33:{'t':8,'u':1,'s':0},
    36:{'t':8,'u':4,'s':0},
    37:{'t':8,'u':5,'s':0},
    40:{'t':8,'u':8,'s':0},
    41:{'t':8,'u':9,'s':1},
    44:{'t':8,'u':12,'s':1},
    45:{'t':8,'u':13,'s':1},
    }


def _get_t_u(v):
    r"""
    Return the parameters of table 7.1 from [Stinson2004]_.

    INPUT:

    - ``v`` -- integer

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import _get_t_u
        sage: _get_t_u(20)
        (5, 0)
    """
    # Table 7.1
    v = int(v)
    global table_7_1
    d = table_7_1[v % 48]
    s = v//48
    if s < d['s']:
        raise RuntimeError("This should not have happened.")
    t = 12*s+d['t']
    u = d['u']
    return t,u

################
# (v,5,1)-BIBD #
################


def v_5_1_BIBD(v, check=True):
    r"""
    Return a `(v,5,1)`-BIBD.

    This method follows the construction from [ClaytonSmith]_.

    INPUT:

    - ``v`` -- integer

    .. SEEALSO::

        * :func:`balanced_incomplete_block_design`

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import v_5_1_BIBD
        sage: i = 0
        sage: while i<200:                                                              # needs sage.libs.pari sage.schemes
        ....:    i += 20
        ....:    _ = v_5_1_BIBD(i+1)
        ....:    _ = v_5_1_BIBD(i+5)

    TESTS:

    Check that the needed difference families are there::

        sage: for v in [21,41,61,81,141,161,281]:                                       # needs sage.libs.pari
        ....:     assert designs.difference_family(v,5,existence=True)
        ....:     _ = designs.difference_family(v,5)
    """
    v = int(v)

    assert (v > 1)
    assert (v % 20 == 5 or v % 20 == 1)  # note: equivalent to (v-1)%4 == 0 and (v*(v-1))%20 == 0

    # Lemma 27
    if v % 5 == 0 and (v//5) % 4 == 1 and is_prime_power(v//5):
        bibd = BIBD_5q_5_for_q_prime_power(v//5)
    # Lemma 28
    elif v in [21,41,61,81,141,161,281]:
        from .difference_family import difference_family
        G,D = difference_family(v,5)
        bibd = BIBD_from_difference_family(G, D, check=False)
    # Lemma 29
    elif v == 165:
        bibd = BIBD_from_PBD(v_5_1_BIBD(41,check=False),165,5,check=False)
    elif v == 181:
        bibd = BIBD_from_PBD(v_5_1_BIBD(45,check=False),181,5,check=False)
    elif v in (201,285,301,401,421,425):
        # Call directly the BIBD_from_TD function
        # note: there are (201,5,1) and (421,5)-difference families that can be
        # obtained from the general constructor
        bibd = BIBD_from_TD(v,5)
    # Theorem 31.2
    elif (v-1)//4 in [80, 81, 85, 86, 90, 91, 95, 96, 110, 111, 115, 116, 120, 121, 250, 251, 255, 256, 260, 261, 265, 266, 270, 271]:
        r = (v-1)//4
        if r <= 96:
            k,t,u = 5, 16, r-80
        elif r <= 121:
            k,t,u = 10, 11, r-110
        else:
            k,t,u = 10, 25, r-250
        bibd = BIBD_from_PBD(PBD_from_TD(k,t,u),v,5,check=False)

    else:
        r,s,t,u = _get_r_s_t_u(v)
        bibd = BIBD_from_PBD(PBD_from_TD(5,t,u),v,5,check=False)

    if check:
        assert is_pairwise_balanced_design(bibd,v,[5])

    return bibd

def _get_r_s_t_u(v):
    r"""
    Implement the table from [ClaytonSmith]_.

    Return the parameters ``r,s,t,u`` associated with an integer ``v``.

    INPUT:

    - ``v`` -- integer

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import _get_r_s_t_u
        sage: _get_r_s_t_u(25)
        (6, 0, 1, 1)
    """
    r = int((v-1)/4)
    s = r//150
    x = r % 150

    if x == 0:
        t,u = 30*s-5,  25
    elif x == 1:
        t,u = 30*s-5,  26
    elif x <= 21:
        t,u = 30*s+1,  x-5
    elif x == 25:
        t,u = 30*s+5,  0
    elif x == 26:
        t,u = 30*s+5,  1
    elif x == 30:
        t,u = 30*s+5,  5
    elif x <= 51:
        t,u = 30*s+5,  x-25
    elif x <= 121:
        t,u = 30*s+11, x-55
    elif x <= 146:
        t,u = 30*s+25, x-125

    return r,s,t,u


def PBD_from_TD(k,t,u):
    r"""
    Return a `(kt,\{k,t\})`-PBD if `u=0` and a `(kt+u,\{k,k+1,t,u\})`-PBD otherwise.

    This is theorem 23 from [ClaytonSmith]_. The PBD is obtained from the blocks
    a truncated `TD(k+1,t)`, to which are added the blocks corresponding to the
    groups of the TD. When `u=0`, a `TD(k,t)` is used instead.

    INPUT:

    - ``k``, ``t``, ``u`` -- integers such that `0\leq u \leq t`

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import PBD_from_TD
        sage: from sage.combinat.designs.bibd import is_pairwise_balanced_design
        sage: PBD = PBD_from_TD(2,2,1); PBD
        [[0, 2, 4], [0, 3], [1, 2], [1, 3, 4], [0, 1], [2, 3]]
        sage: is_pairwise_balanced_design(PBD,2*2+1,[2,3])
        True
    """
    from .orthogonal_arrays import transversal_design
    TD = transversal_design(k+bool(u),t, check=False)
    TD = [[x for x in X if x < k*t+u] for X in TD]
    for i in range(k):
        TD.append(list(range(t*i,t*i+t)))
    if u >= 2:
        TD.append(list(range(k*t,k*t+u)))
    return TD

def BIBD_5q_5_for_q_prime_power(q):
    r"""
    Return a `(5q,5,1)`-BIBD with `q\equiv 1\pmod 4` a prime power.

    See Theorem 24 [ClaytonSmith]_.

    INPUT:

    - ``q`` -- integer; a prime power such that `q\equiv 1\pmod 4`

    EXAMPLES::

        sage: from sage.combinat.designs.bibd import BIBD_5q_5_for_q_prime_power
        sage: for q in [25, 45, 65, 85, 125, 145, 185, 205, 305, 405, 605]: # long time
        ....:     _ = BIBD_5q_5_for_q_prime_power(q/5)
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField

    if q % 4 != 1 or not is_prime_power(q):
        raise ValueError("q is not a prime power or q%4!=1.")

    d = (q-1)//4
    B = []
    F = FiniteField(q,'x')
    a = F.primitive_element()
    L = {b:i for i,b in enumerate(F)}
    for b in L:
        B.append([i*q + L[b] for i in range(5)])
        for i in range(5):
            for j in range(d):
                B.append([        i*q + L[b          ],
                          ((i+1) % 5)*q + L[ a**j+b    ],
                          ((i+1) % 5)*q + L[-a**j+b    ],
                          ((i+4) % 5)*q + L[ a**(j+d)+b],
                          ((i+4) % 5)*q + L[-a**(j+d)+b],
                          ])

    return B


def BIBD_from_arc_in_desarguesian_projective_plane(n,k,existence=False):
    r"""
    Return a `(n,k,1)`-BIBD from a maximal arc in a projective plane.

    This function implements a construction from Denniston [Denniston69]_, who
    describes a maximal :meth:`arc
    <sage.combinat.designs.bibd.BalancedIncompleteBlockDesign.arc>` in a
    :func:`Desarguesian Projective Plane
    <sage.combinat.designs.block_design.DesarguesianProjectivePlaneDesign>` of
    order `2^k`. From two powers of two `n,q` with `n<q`, it produces a
    `((n-1)(q+1)+1,n,1)`-BIBD.

    INPUT:

    - ``n``, ``k`` -- integers; must be powers of two (among other restrictions)

    - ``existence`` -- boolean; whether to return the BIBD obtained through
      this construction (default), or to merely indicate with a boolean return
      value whether this method *can* build the requested BIBD.

    EXAMPLES:

    A `(232,8,1)`-BIBD::

        sage: from sage.combinat.designs.bibd import BIBD_from_arc_in_desarguesian_projective_plane
        sage: from sage.combinat.designs.bibd import BalancedIncompleteBlockDesign
        sage: D = BIBD_from_arc_in_desarguesian_projective_plane(232,8)                 # needs sage.libs.gap sage.modules sage.rings.finite_rings
        sage: BalancedIncompleteBlockDesign(232,D)                                      # needs sage.libs.gap sage.modules sage.rings.finite_rings
        (232,8,1)-Balanced Incomplete Block Design

    A `(120,8,1)`-BIBD::

        sage: D = BIBD_from_arc_in_desarguesian_projective_plane(120,8)                 # needs sage.libs.gap sage.modules sage.rings.finite_rings
        sage: BalancedIncompleteBlockDesign(120,D)                                      # needs sage.libs.gap sage.modules sage.rings.finite_rings
        (120,8,1)-Balanced Incomplete Block Design

    Other parameters::

        sage: all(BIBD_from_arc_in_desarguesian_projective_plane(n,k,existence=True)
        ....:     for n,k in
        ....:       [(120, 8), (232, 8), (456, 8), (904, 8), (496, 16),
        ....:        (976, 16), (1936, 16), (2016, 32), (4000, 32), (8128, 64)])
        True

    Of course, not all can be built this way::

        sage: BIBD_from_arc_in_desarguesian_projective_plane(7,3,existence=True)
        False
        sage: BIBD_from_arc_in_desarguesian_projective_plane(7,3)
        Traceback (most recent call last):
        ...
        ValueError: This function cannot produce a (7,3,1)-BIBD

    REFERENCE:

    .. [Denniston69] \R. H. F. Denniston,
       Some maximal arcs in finite projective planes.
       Journal of Combinatorial Theory 6, no. 3 (1969): 317-319.
       :doi:`10.1016/S0021-9800(69)80095-5`
    """
    q = (n-1)//(k-1)-1
    if (k % 2                 or
        q % 2                 or
        q <= k                or
        n != (k-1)*(q+1)+1    or
        not is_prime_power(k) or
        not is_prime_power(q)):
        if existence:
            return False
        raise ValueError("This function cannot produce a ({},{},1)-BIBD".format(n,k))

    if existence:
        return True

    n = k

    # From now on, the code assumes the notations of [Denniston69] for n,q, so
    # that the BIBD returned by the method will have the requested parameters.

    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
    from sage.libs.gap.libgap import libgap
    from sage.matrix.constructor import Matrix

    K = GF(q,'a')
    one = K.one()

    # An irreducible quadratic form over K[X,Y]
    GO = libgap.GeneralOrthogonalGroup(-1,2,q)
    M = libgap.InvariantQuadraticForm(GO)['matrix']
    M = Matrix(M)
    M = M.change_ring(K)
    Q = lambda xx,yy : M[0,0]*xx**2+(M[0,1]+M[1,0])*xx*yy+M[1,1]*yy**2

    # Here, the additive subgroup H (of order n) of K mentioned in
    # [Denniston69] is the set of all elements of K of degree < log_n
    # (seeing elements of K as polynomials in 'a')

    K_iter = list(K) # faster iterations
    log_n = is_prime_power(n,get_data=True)[1]
    C = [(x,y,one)
         for x in K_iter
         for y in K_iter
         if Q(x,y).polynomial().degree() < log_n]

    from sage.combinat.designs.block_design import DesarguesianProjectivePlaneDesign
    return DesarguesianProjectivePlaneDesign(q).trace(C)._blocks

class PairwiseBalancedDesign(GroupDivisibleDesign):
    r"""
    Pairwise Balanced Design (PBD).

    A Pairwise Balanced Design, or `(v,K,\lambda)`-PBD, is a collection
    `\mathcal B` of blocks defined on a set `X` of size `v`, such that any block
    pair of points `p_1,p_2\in X` occurs in exactly `\lambda` blocks of
    `\mathcal B`. Besides, for every block `B\in \mathcal B` we must have
    `|B|\in K`.

    INPUT:

    - ``points`` -- the underlying set; if ``points`` is an integer `v`, then
      the set is considered to be `\{0, ..., v-1\}`

    - ``blocks`` -- collection of blocks

    - ``K`` -- list of integers of which the sizes of the blocks must be
      elements; set to ``None`` (automatic guess) by default

    - ``lambd`` -- integer; value of `\lambda`, set to `1` by default

    - ``check`` -- boolean; whether to check that the design is a `PBD` with
      the right parameters

    - ``copy`` -- (use with caution) if set to ``False`` then ``blocks`` must be
      a list of lists of integers. The list will not be copied but will be
      modified in place (each block is sorted, and the whole list is
      sorted). Your ``blocks`` object will become the instance's internal data.
    """
    def __init__(self, points, blocks, K=None, lambd=1, check=True, copy=True,**kwds):
        r"""
        Constructor.

        EXAMPLES::

            sage: designs.balanced_incomplete_block_design(13,3) # indirect doctest
            (13,3,1)-Balanced Incomplete Block Design
        """
        try:
            i = int(points)
        except TypeError:
            pass
        else:
            points = list(range(i))

        GroupDivisibleDesign.__init__(self,
                                      points,
                                      [[x] for x in points],
                                      blocks,
                                      K=K,
                                      lambd=lambd,
                                      check=check,
                                      copy=copy,
                                      **kwds)

    def __repr__(self):
        r"""
        Return a string describing the PBD.

        EXAMPLES::

            sage: designs.balanced_incomplete_block_design(13,3) # indirect doctest
            (13,3,1)-Balanced Incomplete Block Design
        """
        bsizes = list(frozenset(self.block_sizes()))
        return "Pairwise Balanced Design on {} points with sets of sizes in {}".format(self.num_points(), bsizes)


class BalancedIncompleteBlockDesign(PairwiseBalancedDesign):
    r"""
    Balanced Incomplete Block Design (BIBD).

    INPUT:

    - ``points`` -- the underlying set. If ``points`` is an integer `v`, then
      the set is considered to be `\{0, ..., v-1\}`

    - ``blocks`` -- collection of blocks

    - ``k`` -- integer; size of the blocks. Set to ``None`` (automatic guess)
      by default

    - ``lambd`` -- integer; value of `\lambda`, set to `1` by default

    - ``check`` -- boolean; whether to check that the design is a `PBD` with
      the right parameters

    - ``copy`` -- (use with caution) if set to ``False`` then ``blocks`` must be
      a list of lists of integers. The list will not be copied but will be
      modified in place (each block is sorted, and the whole list is
      sorted). Your ``blocks`` object will become the instance's internal data.

    EXAMPLES::

        sage: b=designs.balanced_incomplete_block_design(9,3); b
        (9,3,1)-Balanced Incomplete Block Design
    """
    def __init__(self, points, blocks, k=None, lambd=1, check=True, copy=True,**kwds):
        r"""
        Constructor.

        EXAMPLES::

            sage: b=designs.balanced_incomplete_block_design(9,3); b
            (9,3,1)-Balanced Incomplete Block Design
        """
        PairwiseBalancedDesign.__init__(self,
                                        points,
                                        blocks,
                                        K=[k] if k is not None else None,
                                        lambd=lambd,
                                        check=check,
                                        copy=copy,
                                        **kwds)

    def __repr__(self):
        r"""
        A string to describe ``self``.

        EXAMPLES::

            sage: b=designs.balanced_incomplete_block_design(9,3); b
            (9,3,1)-Balanced Incomplete Block Design
        """
        v = self.num_points()
        k = len(self._blocks[0]) if self._blocks else 0
        l = self._lambd
        return "({},{},{})-Balanced Incomplete Block Design".format(v,k,l)

    def arc(self, s=2, solver=None, verbose=0, *, integrality_tolerance=1e-3):
        r"""
        Return the ``s``-arc with maximum cardinality.

        A `s`-arc is a subset of points in a BIBD that intersects each block on
        at most `s` points. It is one possible generalization of independent set
        for graphs.

        A simple counting shows that the cardinality of a `s`-arc is at most
        `(s-1) * r + 1` where `r` is the number of blocks incident to any point.
        A `s`-arc in a BIBD with cardinality `(s-1) * r + 1` is called maximal
        and is characterized by the following property: it is not empty and each
        block either contains `0` or `s` points of this arc. Equivalently, the
        trace of the BIBD on these points is again a BIBD (with block size `s`).

        For more informations, see :wikipedia:`Arc_(projective_geometry)`.

        INPUT:

        - ``s`` -- (default: `2`) the maximum number of points from the arc
          in each block

        - ``solver`` -- (default: ``None``) specify a Mixed Integer Linear
          Programming (MILP) solver to be used. If set to ``None``, the default
          one is used. For more information on MILP solvers and which default
          solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: 0); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- parameter for use with MILP solvers over
          an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        EXAMPLES::

            sage: # needs sage.schemes
            sage: B = designs.balanced_incomplete_block_design(21, 5)
            sage: a2 = B.arc(); a2  # random
            [5, 9, 10, 12, 15, 20]
            sage: len(a2)
            6
            sage: a4 = B.arc(4); a4  # random
            [0, 1, 2, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20]
            sage: len(a4)
            16

        The `2`-arc and `4`-arc above are maximal. One can check that they
        intersect the blocks in either 0 or `s` points. Or equivalently that the
        traces are again BIBD::

            sage: r = (21-1)//(5-1)
            sage: 1 + r*1
            6
            sage: 1 + r*3
            16

            sage: B.trace(a2).is_t_design(2, return_parameters=True)                    # needs sage.schemes
            (True, (2, 6, 2, 1))
            sage: B.trace(a4).is_t_design(2, return_parameters=True)                    # needs sage.schemes
            (True, (2, 16, 4, 1))

        Some other examples which are not maximal::

            sage: # needs sage.numerical.mip
            sage: B = designs.balanced_incomplete_block_design(25, 4)
            sage: a2 = B.arc(2)
            sage: r = (25-1)//(4-1)
            sage: len(a2), 1 + r
            (8, 9)
            sage: sa2 = set(a2)
            sage: set(len(sa2.intersection(b)) for b in B.blocks())
            {0, 1, 2}
            sage: B.trace(a2).is_t_design(2)
            False

            sage: # needs sage.numerical.mip
            sage: a3 = B.arc(3)
            sage: len(a3), 1 + 2*r
            (15, 17)
            sage: sa3 = set(a3)
            sage: set(len(sa3.intersection(b)) for b in B.blocks()) == set([0,3])
            False
            sage: B.trace(a3).is_t_design(3)
            False

        TESTS:

        Test consistency with relabeling::

            sage: b = designs.balanced_incomplete_block_design(7,3)                     # needs sage.schemes
            sage: b.relabel(list("abcdefg"))                                            # needs sage.schemes
            sage: set(b.arc()).issubset(b.ground_set())                                 # needs sage.schemes
            True
        """
        s = int(s)

        # trivial cases
        if s <= 0:
            return []
        elif s >= max(self.block_sizes()):
            return self._points[:]

        # integer linear program
        from sage.numerical.mip import MixedIntegerLinearProgram

        p = MixedIntegerLinearProgram(solver=solver)
        b = p.new_variable(binary=True)
        p.set_objective(p.sum(b[i] for i in range(len(self._points))))
        for i in self._blocks:
            p.add_constraint(p.sum(b[k] for k in i) <= s)
        p.solve(log=verbose)

        values = p.get_values(b, convert=bool, tolerance=integrality_tolerance)
        return [self._points[i] for (i,j) in values.items() if j]


BIBD = BalancedIncompleteBlockDesign
