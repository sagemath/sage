r"""
Documentation for the matroids in the catalog

This module contains implementations for many of the functions accessible
through :mod:`matroids. <sage.matroids.matroids_catalog>` and
:mod:`matroids.named_matroids. <sage.matroids.matroids_catalog>`
(type those lines in Sage and hit ``tab`` for a list).

The docstrings include educational information about each named matroid with
the hopes that this class can be used as a reference. However, for a more
comprehensive list of properties we refer to the appendix of [Oxl2011]_.

.. TODO::

    Add optional argument ``groundset`` to each method so users can customize
    the groundset of the matroid. We probably want some means of relabeling to
    accomplish that.

    Add option to specify the field for represented matroids.

AUTHORS:

- Michael Welsh, Stefan van Zwam (2013-04-01): initial version
- Giorgos Mousa, Andreas Triantafyllos (2023-12-08): reorganization

Functions
=========
"""
# ****************************************************************************
#       Copyright (C) 2013 Michael Welsh <michael@welsh.co.nz >
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.constructor import Matrix
from sage.matroids.constructor import Matroid
from sage.matroids.circuit_closures_matroid import CircuitClosuresMatroid
from sage.matroids.linear_matroid import (
    BinaryMatroid,
    TernaryMatroid,
    QuaternaryMatroid,
)
from sage.rings.finite_rings.finite_field_constructor import GF


def NonVamos():
    """
    Return the non-Vamos matroid.

    The non-Vamos matroid, or `V_8^+` is an 8-element matroid of rank 4. It is
    a tightening of the Vamos matroid. It is representable over some field.
    See [Oxl2011]_, p. 72, 84.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.NonVamos(); M
        NonVamos: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
             {'c', 'd', 'e', 'f'}, {'c', 'd', 'g', 'h'}, {'e', 'f', 'g', 'h'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: setprint(M.nonbases())
        [{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
         {'c', 'd', 'e', 'f'}, {'c', 'd', 'g', 'h'}, {'e', 'f', 'g', 'h'}]
        sage: M.is_dependent(['c', 'd', 'g', 'h'])
        True
        sage: M.is_valid() # long time
        True
    """
    E = "abcdefgh"
    CC = {3: ["abcd", "abef", "cdef", "abgh", "cdgh", "efgh"], 4: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("NonVamos: " + repr(M))
    return M


def NotP8():
    """
    Return the matroid ``NotP8``.

    This is a matroid that is not `P_8`, found on page 512 of [Oxl1992]_ (the
    first edition).

    EXAMPLES::

        sage: M = matroids.named_matroids.P8()
        sage: N = matroids.named_matroids.NotP8()
        sage: M.is_isomorphic(N)
        False
        sage: M.is_valid()
        True
    """
    A = Matrix(
        GF(3),
        [
            [1, 0, 0, 0, 0, 1, 1, -1],
            [0, 1, 0, 0, 1, 0, 1, 1],
            [0, 0, 1, 0, 1, 1, 0, 1],
            [0, 0, 0, 1, -1, 1, 1, 1],
        ],
    )
    M = TernaryMatroid(A, "abcdefgh")
    M.rename("NotP8: " + repr(M))
    return M


def AG23minus():
    """
    Return the ternary affine plane minus a point.

    This is a sixth-roots-of-unity matroid, and an excluded minor for the
    class of near-regular matroids.
    See [Oxl2011]_, p. 653.

    EXAMPLES::

        sage: M = matroids.named_matroids.AG23minus()
        sage: M.is_valid()
        True

    """
    E = "abcdefgh"
    CC = {2: ["abc", "ceh", "fgh", "adf", "aeg", "cdg", "bdh", "bef"], 3: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("AG23minus: " + repr(M))
    return M


def P9():
    """
    Return the matroid `P_9`.

    This is the matroid referred to as `P_9` by Oxley in his paper "The binary
    matroids with no 4-wheel minor"

    EXAMPLES::

        sage: M = matroids.named_matroids.P9()
        sage: M
        P9: Binary matroid of rank 4 on 9 elements, type (1, 1)
        sage: M.is_valid()
        True
    """
    A = Matrix(
        GF(2),
        [
            [1, 0, 0, 0, 1, 0, 0, 1, 1],
            [0, 1, 0, 0, 1, 1, 0, 0, 1],
            [0, 0, 1, 0, 0, 1, 1, 0, 1],
            [0, 0, 0, 1, 0, 0, 1, 1, 0],
        ],
    )
    M = BinaryMatroid(A, "abcdefghi")
    M.rename("P9: " + repr(M))
    return M


def R9A():
    """
    Return the matroid `R_9^A`.

    The matroid `R_9^A` is not representable over any field, yet none of the
    cross-ratios in its Tuttegroup equal 1. It is one of the 4 matroids on at
    most 9 elements with this property, the others being `{R_9^A}^*`, `R_9^B`
    and `{R_9^B}^*`.

    EXAMPLES::

        sage: M = matroids.named_matroids.R9A()
        sage: M.is_valid() # long time
        True

    """
    E = "abcdefghi"
    CC = {
        3: [
            "abde", "bcdf", "aceg", "abch", "aefh", "adgh", "acdi",
            "abfi", "defi", "begi", "bdhi", "cehi", "fghi",
        ],
        4: [E],
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("R9A: " + repr(M))
    return M


def R9B():
    """
    Return the matroid `R_9^B`.

    The matroid `R_9^B` is not representable over any field, yet none of the
    cross-ratios in its Tuttegroup equal 1.
    It is one of the 4 matroids on at most 9 elements with this property, the
    others being `{R_9^B}^*`, `R_9^A` and `{R_9^A}^*`.

    EXAMPLES::

        sage: M = matroids.named_matroids.R9B()
        sage: M.is_valid() # long time
        True

    """
    E = "abcdefghi"
    CC = {
        3: [
            "abde", "bcdf", "aceg", "abch", "befh", "cdgh", "bcei",
            "adfi", "abgi", "degi", "bdhi", "aehi", "fghi",
        ],
        4: [E],
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("R9B: " + repr(M))
    return M


def Block_9_4():
    """
    Return the paving matroid whose non-spanning circuits form the blocks of a
    `2-(9, 4, 3)` design.

    EXAMPLES::

        sage: M = matroids.named_matroids.Block_9_4()
        sage: M.is_valid() # long time
        True
        sage: BD = BlockDesign(M.groundset(), M.nonspanning_circuits())                 # needs sage.graphs
        sage: BD.is_t_design(return_parameters=True)                                    # needs sage.graphs
        (True, (2, 9, 4, 3))
    """
    E = "abcdefghi"
    CC = {
        3: [
            "abcd", "acef", "bdef", "cdeg", "abfg", "adeh",
            "bcfh", "acgh", "begh", "dfgh", "abei", "cdfi",
            "bcgi", "adgi", "efgi", "bdhi", "cehi", "afhi",
        ],
        4: [E],
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("Block(9, 4): " + repr(M))
    return M


def TicTacToe():
    """
    Return the TicTacToe matroid.

    The dual of the TicTacToe matroid is not algebraic; it is unknown whether
    the TicTacToe matroid itself is algebraic. See [Hoc]_.

    EXAMPLES::

        sage: M = matroids.named_matroids.TicTacToe()
        sage: M.is_valid() # long time
        True

    """
    E = "abcdefghi"
    CC = {
        4: ["abcdg", "adefg", "abceh", "abcfi",
            "cdefi", "adghi", "beghi", "cfghi"],
        5: [E],
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("TicTacToe: " + repr(M))
    return M


def N1():
    r"""
    Return the matroid `N_1`, represented over `\GF{3}`.

    `N_1` is an excluded minor for the dyadic matroids. See [Oxl2011]_, p. 554.

    EXAMPLES::

        sage: M = matroids.named_matroids.N1()
        sage: M.is_field_isomorphic(M.dual())
        True
        sage: M.is_valid()
        True

    """
    A = Matrix(
        GF(3),
        [
            [1, 0, 0, 0, 0, 2, 0, 0, 1, 1],
            [0, 1, 0, 0, 0, 1, 2, 0, 0, 1],
            [0, 0, 1, 0, 0, 0, 1, 2, 0, 1],
            [0, 0, 0, 1, 0, 0, 0, 1, 2, 2],
            [0, 0, 0, 0, 1, 1, 1, 1, 2, 0],
        ],
    )
    M = TernaryMatroid(A, "abcdefghij")
    M.rename("N1: " + repr(M))
    return M


def Block_10_5():
    """
    Return the paving matroid whose non-spanning circuits form the blocks of a
    `3-(10, 5, 3)` design.

    EXAMPLES::

        sage: M = matroids.named_matroids.Block_10_5()
        sage: M.is_valid() # long time
        True
        sage: BD = BlockDesign(M.groundset(), M.nonspanning_circuits())                 # needs sage.graphs
        sage: BD.is_t_design(return_parameters=True)                                    # needs sage.graphs
        (True, (3, 10, 5, 3))
    """
    E = "abcdefghij"
    CC = {
        4: [
            "abcde", "acdfg", "bdefg", "bcdfh", "abefh", "abcgh", "adegh",
            "cefgh", "bcefi", "adefi", "bcdgi", "acegi", "abfgi", "abdhi",
            "cdehi", "acfhi", "beghi", "dfghi", "abdfj", "acefj", "abegj",
            "cdegj", "bcfgj", "acdhj", "bcehj", "defhj", "bdghj", "afghj",
            "abcij", "bdeij", "cdfij", "adgij", "efgij", "aehij", "bfhij",
            "cghij",
        ],
        5: [E],
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("Block(10, 5): " + repr(M))
    return M


def Q10():
    r"""
    Return the matroid `Q_{10}`, represented over `\GF{4}`.

    `Q_{10}` is a 10-element, rank-5, self-dual matroid. It is representable
    over `\GF{3}` and `\GF{4}`, and hence is a sixth-roots-of-unity matroid.
    `Q_{10}` is a splitter for the class of sixth-root-of-unity matroids.

    EXAMPLES::

        sage: M = matroids.named_matroids.Q10()                                         # needs sage.rings.finite_rings
        sage: M.is_isomorphic(M.dual())                                                 # needs sage.rings.finite_rings
        True
        sage: M.is_valid()                                                              # needs sage.rings.finite_rings
        True

    Check the splitter property. By Seymour's Theorem, and using self-duality,
    we only need to check that all 3-connected single-element extensions have
    an excluded minor for sixth-roots-of-unity. The only excluded minors that
    are quaternary are `U_{2, 5}, U_{3, 5}, F_7, F_7^*`. As it happens, it
    suffices to check for `U_{2, 5}`:

        sage: S = matroids.named_matroids.Q10().linear_extensions(simple=True)          # needs sage.rings.finite_rings
        sage: [M for M in S if not M.has_line_minor(5)]         # long time, needs sage.rings.finite_rings
        []
    """
    F = GF(4, "x")
    x = F.gens()[0]
    A = Matrix(
        F,
        [
            [1, 0, 0, 0, 0, 1, x, 0, 0, x + 1],
            [0, 1, 0, 0, 0, x + 1, 1, x, 0, 0],
            [0, 0, 1, 0, 0, 0, x + 1, 1, x, 0],
            [0, 0, 0, 1, 0, 0, 0, x + 1, 1, x],
            [0, 0, 0, 0, 1, x, 0, 0, x + 1, 1],
        ],
    )
    M = QuaternaryMatroid(A, "abcdefghij")
    M.rename("Q10: " + repr(M))
    return M


def BetsyRoss():
    """
    Return the Betsy Ross matroid, represented by circuit closures.

    An extremal golden-mean matroid.
    That is, if `M` is simple, rank 3, has the Betsy Ross matroid as a
    restriction and is a Golden Mean matroid, then `M` is the Betsy Ross
    matroid.

    EXAMPLES::

        sage: M = matroids.named_matroids.BetsyRoss()
        sage: len(M.circuit_closures()[2])
        10
        sage: M.is_valid() # long time
        True
    """
    E = "abcdefghijk"
    CC = {
        2: ["acfg", "bdgh", "cehi", "befj", "adij",
            "dfk", "egk", "ahk", "bik", "cjk"],
        3: [E],
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("BetsyRoss: " + repr(M))
    return M


def N2():
    r"""
    Return the matroid `N_2`, represented over `\GF{3}`.

    `N_2` is an excluded minor for the dyadic matroids. See [Oxl2011]_, p. 554.

    EXAMPLES::

        sage: M = matroids.named_matroids.N2()
        sage: M.is_field_isomorphic(M.dual())
        True
        sage: M.is_valid()
        True

    """
    A = Matrix(
        GF(3),
        [
            [1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1, 1],
            [0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 0, 1],
            [0, 0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 1],
            [0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 1, 0],
            [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1],
            [0, 0, 0, 0, 0, 1, 1, 2, 2, 1, 0, 1],
        ],
    )
    M = TernaryMatroid(A, "abcdefghijkl")
    M.rename("N2: " + repr(M))
    return M


def ExtendedTernaryGolayCode():
    """
    Return the matroid of the extended ternary Golay code.

    See
    :class:`GolayCode <sage.coding.golay_code.GolayCode>`

    EXAMPLES::

        sage: M = matroids.named_matroids.ExtendedTernaryGolayCode()
        sage: C = LinearCode(M.representation())
        sage: C.is_permutation_equivalent(codes.GolayCode(GF(3)))       # long time, needs sage.rings.finite_rings
        True
        sage: M.is_valid()
        True
    """
    A = Matrix(
        GF(3),
        [
            [1, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 0],
            [0, 1, 0, 0, 0, 0, 1, 1, 2, 1, 0, 2],
            [0, 0, 1, 0, 0, 0, 1, 2, 1, 0, 1, 2],
            [0, 0, 0, 1, 0, 0, 1, 2, 0, 1, 2, 1],
            [0, 0, 0, 0, 1, 0, 1, 0, 2, 2, 1, 1],
            [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1],
        ],
    )
    M = TernaryMatroid(A, "abcdefghijkl")
    M.rename("Extended Ternary Golay Code: " + repr(M))
    return M


def D16():  # A.K.A. the Carolyn Chun Matroid
    """
    Return the matroid `D_{16}`.

    Let `M` be a 4-connected binary matroid and `N` an internally 4-connected
    proper minor of `M` with at least 7 elements. Then some element of `M` can
    be deleted or contracted preserving an `N`-minor, unless `M` is `D_{16}`.
    See [CMO2012]_.

    EXAMPLES::

        sage: M = matroids.named_matroids.D16()
        sage: M
        D16: Binary matroid of rank 8 on 16 elements, type (0, 0)
        sage: M.is_valid()
        True

    """
    A = Matrix(
        GF(2),
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1],
            [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0],
        ],
    )
    M = BinaryMatroid(A, "abcdefghijklmnop")
    M.rename("D16: " + repr(M))
    return M


def Terrahawk():  # A.K.A. the Dillon Mayhew Matroid
    """
    Return the Terrahawk matroid.

    The Terrahawk is a binary matroid that is a sporadic exception in a chain
    theorem for internally 4-connected binary matroids. See [CMO2011]_.

    EXAMPLES::

        sage: M = matroids.named_matroids.Terrahawk()
        sage: M
        Terrahawk: Binary matroid of rank 8 on 16 elements, type (0, 4)
        sage: M.is_valid()
        True

    """
    A = Matrix(
        GF(2),
        [
            [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0],
        ],
    )
    M = BinaryMatroid(A, "abcdefghijklmnop")
    M.rename("Terrahawk: " + repr(M))
    return M


def ExtendedBinaryGolayCode():
    """
    Return the matroid of the extended binary Golay code.

    See
    :class:`GolayCode <sage.coding.golay_code.GolayCode>`
    documentation for more on this code.

    EXAMPLES::

        sage: M = matroids.named_matroids.ExtendedBinaryGolayCode()
        sage: C = LinearCode(M.representation())
        sage: C.is_permutation_equivalent(codes.GolayCode(GF(2)))       # long time, needs sage.rings.finite_rings
        True
        sage: M.is_valid()
        True
    """
    A = Matrix(
        GF(2),
        [
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0],
            [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1],
            [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0],
            [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0],
            [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0],
            [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1],
            [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1],
            [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        ],
    )
    M = BinaryMatroid(A, "abcdefghijklmnopqrstuvwx")
    M.rename("Extended Binary Golay Code: " + repr(M))
    return M


def CompleteGraphic(n):
    """
    Return the cycle matroid of the complete graph on `n` vertices.

    INPUT:

    - ``n`` -- an integer, the number of vertices of the underlying complete
      graph.

    OUTPUT:

    The graphic matroid associated with the `n`-vertex complete graph.
    This matroid has rank `n - 1`.

    EXAMPLES::

        sage: # needs sage.graphs
        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.CompleteGraphic(5); M
        M(K5): Graphic matroid of rank 4 on 10 elements
        sage: M.has_minor(matroids.Uniform(2, 4))
        False
        sage: simplify(M.contract(randrange(0,
        ....:                 10))).is_isomorphic(matroids.CompleteGraphic(4))
        True
        sage: setprint(M.closure([0, 2, 4, 5]))
        {0, 1, 2, 4, 5, 7}
        sage: M.is_valid()
        True
    """
    from sage.graphs.graph_generators import graphs

    M = Matroid(
        groundset=list(range((n * (n - 1)) // 2)),
        graph=graphs.CompleteGraph(n)
    )
    M.rename("M(K" + str(n) + "): " + repr(M))
    return M
