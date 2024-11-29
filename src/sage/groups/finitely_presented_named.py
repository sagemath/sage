r"""
Named Finitely Presented Groups

Construct groups of small order and "named" groups as quotients of free groups. These
groups are available through tab completion by typing ``groups.presentation.<tab>``
or by importing the required methods. Tab completion is made available through
Sage's :ref:`group catalog <sage.groups.groups_catalog>`. Some examples are engineered
from entries in [TW1980]_.

Groups available as finite presentations:

- Alternating group, `A_n` of order `n!/2` --
  :func:`groups.presentation.Alternating <sage.groups.finitely_presented_named.AlternatingPresentation>`

- the `n`-fruit Cactus group, a standard notation for which is `J_n` --
  :func:`groups.presentation.Cactus <sage.groups.finitely_presented_named.CactusPresentation>`

- Cyclic group, `C_n` of order `n` --
  :func:`groups.presentation.Cyclic <sage.groups.finitely_presented_named.CyclicPresentation>`

- Dicyclic group, nonabelian groups of order `4n` with a unique element of
  order 2 --
  :func:`groups.presentation.DiCyclic <sage.groups.finitely_presented_named.DiCyclicPresentation>`

- Dihedral group, `D_n` of order `2n` --
  :func:`groups.presentation.Dihedral <sage.groups.finitely_presented_named.DihedralPresentation>`

- Finitely generated abelian group, `\ZZ_{n_1} \times \ZZ_{n_2} \times \cdots \times \ZZ_{n_k}` --
  :func:`groups.presentation.FGAbelian <sage.groups.finitely_presented_named.FinitelyGeneratedAbelianPresentation>`

- Finitely generated Heisenberg group --
  :func:`groups.presentation.Heisenberg <sage.groups.finitely_presented_named.FinitelyGeneratedHeisenbergPresentation>`

- Klein four group, `C_2 \times C_2` --
  :func:`groups.presentation.KleinFour <sage.groups.finitely_presented_named.KleinFourPresentation>`

- Quaternion group of order 8 --
  :func:`groups.presentation.Quaternion <sage.groups.finitely_presented_named.QuaternionPresentation>`

- Symmetric group, `S_n` of order `n!` --
  :func:`groups.presentation.Symmetric <sage.groups.finitely_presented_named.SymmetricPresentation>`

AUTHORS:

- Davis Shurbert (2013-06-21): initial version

EXAMPLES::

    sage: groups.presentation.Cyclic(4)
    Finitely presented group < a | a^4 >

You can also import the desired functions::

    sage: from sage.groups.finitely_presented_named import CyclicPresentation
    sage: CyclicPresentation(4)
    Finitely presented group < a | a^4 >
"""
# ****************************************************************************
#       Copyright (C) 2013 Davis Shurbert <davis.sprout@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer import Integer
from sage.groups.free_group import FreeGroup
from sage.groups.finitely_presented import FinitelyPresentedGroup
from sage.libs.gap.libgap import libgap
from sage.matrix.constructor import diagonal_matrix
from sage.modules.fg_pid.fgp_module import FGP_Module
from sage.rings.integer_ring import ZZ


def CyclicPresentation(n):
    r"""
    Build cyclic group of order `n` as a finitely presented group.

    INPUT:

    - ``n`` -- the order of the cyclic presentation to be returned

    OUTPUT: the cyclic group of order `n` as finite presentation

    EXAMPLES::

        sage: groups.presentation.Cyclic(10)
        Finitely presented group < a | a^10 >
        sage: n = 8; C = groups.presentation.Cyclic(n)
        sage: C.as_permutation_group().is_isomorphic(CyclicPermutationGroup(n))
        True

    TESTS::

        sage: groups.presentation.Cyclic(0)
        Traceback (most recent call last):
        ...
        ValueError: finitely presented group order must be positive
    """
    n = Integer(n)
    if n < 1:
        raise ValueError('finitely presented group order must be positive')
    F = FreeGroup('a')
    rls = F([1])**n,
    return FinitelyPresentedGroup(F, rls)


def FinitelyGeneratedAbelianPresentation(int_list):
    r"""
    Return canonical presentation of finitely generated abelian group.

    INPUT:

    - ``int_list`` -- list of integers defining the group to be returned, the defining list
      is reduced to the invariants of the input list before generating the corresponding
      group

    OUTPUT:

    Finitely generated abelian group, `\ZZ_{n_1} \times \ZZ_{n_2} \times \cdots \times \ZZ_{n_k}`
    as a finite presentation, where `n_i` forms the invariants of the input list.

    EXAMPLES::

        sage: groups.presentation.FGAbelian([2,2])
        Finitely presented group < a, b | a^2, b^2, a^-1*b^-1*a*b >
        sage: groups.presentation.FGAbelian([2,3])
        Finitely presented group < a | a^6 >
        sage: groups.presentation.FGAbelian([2,4])
        Finitely presented group < a, b | a^2, b^4, a^-1*b^-1*a*b >

    You can create free abelian groups::

        sage: groups.presentation.FGAbelian([0])
        Finitely presented group < a |  >
        sage: groups.presentation.FGAbelian([0,0])
        Finitely presented group < a, b | a^-1*b^-1*a*b >
        sage: groups.presentation.FGAbelian([0,0,0])
        Finitely presented group < a, b, c | a^-1*b^-1*a*b, a^-1*c^-1*a*c, b^-1*c^-1*b*c >

    And various infinite abelian groups::

        sage: groups.presentation.FGAbelian([0,2])
        Finitely presented group < a, b | a^2, a^-1*b^-1*a*b >
        sage: groups.presentation.FGAbelian([0,2,2])
        Finitely presented group < a, b, c | a^2, b^2, a^-1*b^-1*a*b, a^-1*c^-1*a*c, b^-1*c^-1*b*c >

    Outputs are reduced to minimal generators and relations::

        sage: groups.presentation.FGAbelian([3,5,2,7,3])
        Finitely presented group < a, b | a^3, b^210, a^-1*b^-1*a*b >
        sage: groups.presentation.FGAbelian([3,210])
        Finitely presented group < a, b | a^3, b^210, a^-1*b^-1*a*b >

    The trivial group is an acceptable output::

        sage: groups.presentation.FGAbelian([])
        Finitely presented group <  |  >
        sage: groups.presentation.FGAbelian([1])
        Finitely presented group <  |  >
        sage: groups.presentation.FGAbelian([1,1,1,1,1,1,1,1,1,1])
        Finitely presented group <  |  >

    Input list must consist of positive integers::

        sage: groups.presentation.FGAbelian([2,6,3,9,-4])
        Traceback (most recent call last):
        ...
        ValueError: input list must contain nonnegative entries
        sage: groups.presentation.FGAbelian([2,'a',4])
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 'a' to an integer

    TESTS::

        sage: ag = groups.presentation.FGAbelian([2,2])
        sage: ag.as_permutation_group().is_isomorphic(groups.permutation.KleinFour())
        True
        sage: G = groups.presentation.FGAbelian([2,4,8])
        sage: C2 = CyclicPermutationGroup(2)
        sage: C4 = CyclicPermutationGroup(4)
        sage: C8 = CyclicPermutationGroup(8)
        sage: gg = (C2.direct_product(C4)[0]).direct_product(C8)[0]
        sage: gg.is_isomorphic(G.as_permutation_group())
        True
        sage: all(groups.presentation.FGAbelian([i]).as_permutation_group().is_isomorphic(groups.presentation.Cyclic(i).as_permutation_group()) for i in [2..35])
        True
    """
    from sage.groups.free_group import _lexi_gen
    check_ls = [Integer(x) for x in int_list if Integer(x) >= 0]
    if len(check_ls) != len(int_list):
        raise ValueError('input list must contain nonnegative entries')

    col_sp = diagonal_matrix(int_list).column_space()
    invariants = FGP_Module(ZZ**(len(int_list)), col_sp).invariants()
    name_gen = _lexi_gen()
    F = FreeGroup([next(name_gen) for i in invariants])
    ret_rls = [F([i+1])**invariants[i] for i in range(len(invariants)) if invariants[i] != 0]

    # Build commutator relations
    gen_pairs = [[F.gen(i),F.gen(j)] for i in range(F.ngens()-1) for j in range(i+1,F.ngens())]
    ret_rls = ret_rls + [x[0]**(-1)*x[1]**(-1)*x[0]*x[1] for x in gen_pairs]
    return FinitelyPresentedGroup(F, tuple(ret_rls))


def FinitelyGeneratedHeisenbergPresentation(n=1, p=0):
    r"""
    Return a finite presentation of the Heisenberg group.

    The Heisenberg group is the group of `(n+2) \times (n+2)` matrices
    over a ring `R` with diagonal elements equal to 1, first row and
    last column possibly nonzero, and all the other entries equal to zero.

    INPUT:

    - ``n`` -- the degree of the Heisenberg group

    - ``p`` -- (optional) a prime number, where we construct the
      Heisenberg group over the finite field `\ZZ/p\ZZ`

    OUTPUT: finitely generated Heisenberg group over the finite field
    of order ``p`` or over the integers

    .. SEEALSO::

        :class:`~sage.groups.matrix_gps.heisenberg.HeisenbergGroup`

    EXAMPLES::

        sage: H = groups.presentation.Heisenberg(); H
        Finitely presented group < x1, y1, z |
         x1*y1*x1^-1*y1^-1*z^-1, z*x1*z^-1*x1^-1, z*y1*z^-1*y1^-1 >
        sage: H.order()
        +Infinity
        sage: r1, r2, r3 = H.relations()
        sage: A = matrix([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
        sage: B = matrix([[1, 0, 0], [0, 1, 1], [0, 0, 1]])
        sage: C = matrix([[1, 0, 1], [0, 1, 0], [0, 0, 1]])
        sage: r1(A, B, C)
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: r2(A, B, C)
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: r3(A, B, C)
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: p = 3
        sage: Hp = groups.presentation.Heisenberg(p=3)
        sage: Hp.order() == p**3
        True
        sage: Hnp = groups.presentation.Heisenberg(n=2, p=3)
        sage: len(Hnp.relations())
        13

    REFERENCES:

    - :wikipedia:`Heisenberg_group`
    """
    n = Integer(n)
    if n < 1:
        raise ValueError('n must be a positive integer')

    # generators' names are x1, .., xn, y1, .., yn, z
    vx = ['x' + str(i) for i in range(1,n+1)]
    vy = ['y' + str(i) for i in range(1,n+1)]
    str_generators = ', '.join(vx + vy + ['z'])

    F = FreeGroup(str_generators)
    x = F.gens()[0:n] # list of generators x1, x2, ..., xn
    y = F.gens()[n:2*n] # list of generators x1, x2, ..., xn
    z = F.gen(n*2)

    def commutator(a, b):
        return a * b * a**-1 * b**-1
    # First set of relations: [xi, yi] = z
    r1 = [commutator(x[i], y[i]) * z**-1 for i in range(n)]
    # Second set of relations: [z, xi] = 1
    r2 = [commutator(z, x[i]) for i in range(n)]
    # Third set of relations: [z, yi] = 1
    r3 = [commutator(z, y[i]) for i in range(n)]
    # Fourth set of relations: [xi, yi] = 1 for i != j
    r4 = [commutator(x[i], y[j]) for i in range(n) for j in range(n) if i != j]
    rls = r1 + r2 + r3 + r4

    from sage.sets.primes import Primes
    if p not in Primes() and p != 0:
        raise ValueError("p must be 0 or a prime number")
    if p > 0:
        rls += [w**p for w in F.gens()]
    return FinitelyPresentedGroup(F, tuple(rls))


def DihedralPresentation(n):
    r"""
    Build the Dihedral group of order `2n` as a finitely presented group.

    INPUT:

    - ``n`` -- the size of the set that `D_n` is acting on

    OUTPUT: Dihedral group of order `2n`

    EXAMPLES::

        sage: D = groups.presentation.Dihedral(7); D
        Finitely presented group < a, b | a^7, b^2, (a*b)^2 >
        sage: D.as_permutation_group().is_isomorphic(DihedralGroup(7))
        True

    TESTS::

        sage: n = 9
        sage: D = groups.presentation.Dihedral(n)
        sage: D.ngens() == 2
        True
        sage: groups.presentation.Dihedral(0)
        Traceback (most recent call last):
        ...
        ValueError: finitely presented group order must be positive
    """
    n = Integer( n )
    if n < 1:
        raise ValueError('finitely presented group order must be positive')
    F = FreeGroup([ 'a', 'b' ])
    rls = F([1])**n, F([2])**2, (F([1])*F([2]))**2
    return FinitelyPresentedGroup( F, rls )


def DiCyclicPresentation(n):
    r"""
    Build the dicyclic group of order `4n`, for `n \geq 2`, as a finitely
    presented group.

    INPUT:

    - ``n`` -- positive integer, 2 or greater, determining the order of
      the group (`4n`)

    OUTPUT: the dicyclic group of order `4n` is defined by the presentation

    .. MATH::

        \langle a, x \mid a^{2n}=1, x^{2}=a^{n}, x^{-1}ax=a^{-1} \rangle

    .. NOTE::

        This group is also available as a permutation group via
        :class:`groups.permutation.DiCyclic <sage.groups.perm_gps.permgroup_named.DiCyclicGroup>`.

    EXAMPLES::

        sage: D = groups.presentation.DiCyclic(9); D
        Finitely presented group < a, b | a^18, b^2*a^-9, b^-1*a*b*a >
        sage: D.as_permutation_group().is_isomorphic(groups.permutation.DiCyclic(9))
        True

    TESTS::

        sage: Q = groups.presentation.DiCyclic(2)
        sage: Q.as_permutation_group().is_isomorphic(QuaternionGroup())
        True
        sage: for i in [5, 8, 12, 32]:
        ....:     A = groups.presentation.DiCyclic(i).as_permutation_group()
        ....:     B = groups.permutation.DiCyclic(i)
        ....:     assert A.is_isomorphic(B)
        sage: groups.presentation.DiCyclic(1)
        Traceback (most recent call last):
        ...
        ValueError: input integer must be greater than 1
    """
    n = Integer(n)
    if n < 2:
        raise ValueError('input integer must be greater than 1')

    F = FreeGroup(['a','b'])
    rls = F([1])**(2*n), F([2,2])*F([-1])**n, F([-2,1,2,1])
    return FinitelyPresentedGroup(F, rls)


def SymmetricPresentation(n):
    r"""
    Build the Symmetric group of order `n!` as a finitely presented group.

    INPUT:

    - ``n`` -- the size of the underlying set of arbitrary symbols being acted
      on by the Symmetric group of order `n!`

    OUTPUT:

    Symmetric group as a finite presentation, implementation uses GAP to find an
    isomorphism from a permutation representation to a finitely presented group
    representation. Due to this fact, the exact output presentation may not be
    the same for every method call on a constant ``n``.

    EXAMPLES::

        sage: S4 = groups.presentation.Symmetric(4)
        sage: S4.as_permutation_group().is_isomorphic(SymmetricGroup(4))
        True

    TESTS::

        sage: S = [groups.presentation.Symmetric(i) for i in range(1,4)]; S[0].order()
        1
        sage: S[1].order(), S[2].as_permutation_group().is_isomorphic(DihedralGroup(3))
        (2, True)
        sage: S5 = groups.presentation.Symmetric(5)
        sage: perm_S5 = S5.as_permutation_group(); perm_S5.is_isomorphic(SymmetricGroup(5))
        True
        sage: groups.presentation.Symmetric(8).order()
        40320
    """
    from sage.groups.perm_gps.permgroup_named import SymmetricGroup
    from sage.groups.free_group import _lexi_gen

    n = Integer(n)
    if n <= 1:
        return FinitelyPresentedGroup(FreeGroup(()), ())

    perm_rep = SymmetricGroup(n)
    GAP_fp_rep = libgap.Image(libgap.IsomorphismFpGroupByGenerators(perm_rep, perm_rep.gens()))
    image_gens = GAP_fp_rep.FreeGeneratorsOfFpGroup()
    name_itr = _lexi_gen()  # Python generator object for variable names
    F = FreeGroup([next(name_itr) for x in perm_rep.gens()])
    ret_rls = tuple([F(rel_word.TietzeWordAbstractWord(image_gens).sage())
                     for rel_word in GAP_fp_rep.RelatorsOfFpGroup()])
    return FinitelyPresentedGroup(F, ret_rls)


def QuaternionPresentation():
    r"""
    Build the Quaternion group of order 8 as a finitely presented group.

    OUTPUT: Quaternion group as a finite presentation

    EXAMPLES::

        sage: Q = groups.presentation.Quaternion(); Q
        Finitely presented group < a, b | a^4, b^2*a^-2, a*b*a*b^-1 >
        sage: Q.as_permutation_group().is_isomorphic(QuaternionGroup())
        True

    TESTS::

        sage: Q = groups.presentation.Quaternion()
        sage: Q.order(), Q.is_abelian()
        (8, False)
        sage: Q.is_isomorphic(groups.presentation.DiCyclic(2))
        #I  Forcing finiteness test
        True
    """
    F = FreeGroup(['a','b'])
    rls = F([1])**4, F([2,2,-1,-1]), F([1,2,1,-2])
    return FinitelyPresentedGroup(F, rls)


def AlternatingPresentation(n):
    r"""
    Build the Alternating group of order `n!/2` as a finitely presented group.

    INPUT:

    - ``n`` -- the size of the underlying set of arbitrary symbols being acted
      on by the Alternating group of order `n!/2`

    OUTPUT:

    Alternating group as a finite presentation, implementation uses GAP to find an
    isomorphism from a permutation representation to a finitely presented group
    representation. Due to this fact, the exact output presentation may not be
    the same for every method call on a constant ``n``.

    EXAMPLES::

        sage: A6 = groups.presentation.Alternating(6)
        sage: A6.as_permutation_group().is_isomorphic(AlternatingGroup(6)), A6.order()
        (True, 360)

    TESTS:

    Even permutation tests::

        sage: A1 = groups.presentation.Alternating(1); A2 = groups.presentation.Alternating(2)
        sage: A1.is_isomorphic(A2), A1.order()
        (True, 1)
        sage: A3 = groups.presentation.Alternating(3); A3.order(), A3.as_permutation_group().is_cyclic()
        (3, True)
        sage: A8 = groups.presentation.Alternating(8); A8.order()
        20160
    """
    from sage.groups.perm_gps.permgroup_named import AlternatingGroup
    from sage.groups.free_group import _lexi_gen

    n = Integer(n)
    if n <= 2:
        return FinitelyPresentedGroup(FreeGroup(()), ())

    perm_rep = AlternatingGroup(n)
    GAP_fp_rep = libgap.Image(libgap.IsomorphismFpGroupByGenerators(perm_rep, perm_rep.gens()))
    image_gens = GAP_fp_rep.FreeGeneratorsOfFpGroup()
    name_itr = _lexi_gen()  # Python generator object for variable names
    F = FreeGroup([next(name_itr) for x in perm_rep.gens()])
    ret_rls = tuple([F(rel_word.TietzeWordAbstractWord(image_gens).sage())
                     for rel_word in GAP_fp_rep.RelatorsOfFpGroup()])
    return FinitelyPresentedGroup(F, ret_rls)


def KleinFourPresentation():
    r"""
    Build the Klein group of order `4` as a finitely presented group.

    OUTPUT: Klein four group (`C_2 \times C_2`) as a finitely presented group

    EXAMPLES::

        sage: K = groups.presentation.KleinFour(); K
        Finitely presented group < a, b | a^2, b^2, a^-1*b^-1*a*b >
    """
    F = FreeGroup(['a','b'])
    rls = F([1])**2, F([2])**2, F([-1])*F([-2])*F([1])*F([2])
    return FinitelyPresentedGroup(F, rls)


def BinaryDihedralPresentation(n):
    r"""
    Build a binary dihedral group of order `4n` as a finitely presented group.

    The binary dihedral group `BD_n` has the following presentation
    (note that there is a typo in [Sun2010]_):

    .. MATH::

        BD_n = \langle x, y, z | x^2 = y^2 = z^n = x y z \rangle.

    INPUT:

    - ``n`` -- the value `n`

    OUTPUT: the binary dihedral group of order `4n` as finite presentation

    EXAMPLES::

        sage: groups.presentation.BinaryDihedral(9)
        Finitely presented group < x, y, z | x^-2*y^2, x^-2*z^9, x^-1*y*z >

    TESTS::

        sage: for n in range(3, 9):                                                     # needs sage.modules
        ....:     P = groups.presentation.BinaryDihedral(n)
        ....:     M = groups.matrix.BinaryDihedral(n)
        ....:     assert P.is_isomorphic(M)
        #I  Forcing finiteness test
        #I  Forcing finiteness test
        #I  Forcing finiteness test
        #I  Forcing finiteness test
        #I  Forcing finiteness test
        #I  Forcing finiteness test
    """
    F = FreeGroup('x,y,z')
    x,y,z = F.gens()
    rls = (x**-2 * y**2, x**-2 * z**n, x**-2 * x*y*z)
    return FinitelyPresentedGroup(F, rls)


def CactusPresentation(n):
    r"""
    Build the `n`-fruit cactus group as a finitely presented group.

    OUTPUT: Cactus group `J_n` as a finitely presented group

    EXAMPLES::

        sage: J3 = groups.presentation.Cactus(3); J3                                    # needs sage.graphs
        Finitely presented group < s12, s13, s23 |
         s12^2, s13^2, s23^2, s13*s12*s13^-1*s23^-1, s13*s23*s13^-1*s12^-1 >
    """
    from sage.groups.cactus_group import CactusGroup
    G = CactusGroup(n)
    F = FreeGroup(G.variable_names())
    gens = F.gens()
    rls = [g**2 for g in gens]
    Gg = G.group_generators()
    K = Gg.keys()
    for i,key in enumerate(K):
        for j,key2 in enumerate(K):
            if i == j:
                continue
            x,y = (Gg[key] * Gg[key2])._data
            if key == x and key2 == y:
                continue
            elt = gens[i] * gens[j] * ~gens[K.index(y)] * ~gens[K.index(x)]
            rls.append(elt)
    return FinitelyPresentedGroup(F, tuple(rls))
