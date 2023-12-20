r"""
Collection of Oxley's matroids

This module contains implementations of Oxley's interesting matroids,
accessible through :mod:`matroids.catalog. <sage.matroids.catalog>` (type
those lines in Sage and hit :kbd:`Tab` for a list).

The docstrings include educational information about each named matroid with
the hopes that this class can be used as a reference. However, for a more
comprehensive list of properties we refer to the appendix of [Oxl2011]_.

AUTHORS:

- Michael Welsh, Stefan van Zwam (2013-04-01): initial version
- Giorgos Mousa, Andreas Triantafyllos (2023-12-08): more matroids

.. TODO::

    Add option to specify the field for represented matroids.

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
    RegularMatroid,
    BinaryMatroid,
    TernaryMatroid,
    QuaternaryMatroid
)
from sage.matroids.database.various_matroids import CompleteGraphic
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.schemes.projective.projective_space import ProjectiveSpace


# The order is the same as in Oxley.


def U24():
    r"""
    The uniform matroid of rank `2` on `4` elements.

    The `4`-point line; isomorphic to `\mathcal{W}^2` , the rank-`2` whirl.
    The unique excluded minor for the class of binary matroids.

    See [Oxl2011]_, p. 639.

    EXAMPLES::

        sage: M = matroids.Uniform(2, 4)
        sage: N = matroids.catalog.U24()
        sage: M.is_isomorphic(N)
        True
        sage: W2 = matroids.Whirl(2)
        sage: W2.is_isomorphic(N)
        True
        sage: N.is_3connected()
        True

    """
    M = Uniform(2, 4)
    M.rename("U(2, 4): " + repr(M))
    return M


def U25():
    """
    The uniform matroid of rank `2` on `5` elements.

    `U_{2,5}` is the `5`-point line. Dual to `U_{3,5}`.

    See [Oxl2011]_, p. 640.

    EXAMPLES::

        sage: U25 = matroids.catalog.U25()
        sage: U35 = matroids.catalog.U35()
        sage: U25.is_isomorphic(U35.dual())
        True

    """
    M = Uniform(2, 5)
    M.rename("U(2, 5): " + repr(M))
    return M


def U35():
    """
    The uniform matroid of rank `3` on `5` elements.

    `U_{3,5}` is five points freely placed in the plane. Dual to `U_{2,5}`.

    See [Oxl2011]_, p. 640.

    EXAMPLES::

        sage: U35 = matroids.catalog.U35()
        sage: U25 = matroids.catalog.U25()
        sage: U35.is_isomorphic(U25.dual())
        True

    """
    M = Uniform(3, 5)
    M.rename("U(3, 5): " + repr(M))
    return M


def K4():
    """
    The graphic matroid of the complete graph `K_4`.

    Isomorphic to `M(W_3)`, the rank-`3` wheel, and to the tipless binary
    `3`-spike.

    See [Oxl2011]_, p. 640.

    EXAMPLES::

        sage: M = matroids.catalog.K4()
        sage: W3 = matroids.Wheel(3)
        sage: M.is_isomorphic(W3)
        True
        sage: Z = matroids.Z(3, False)
        sage: M.is_isomorphic(Z)
        True

    """
    M = CompleteGraphic(4)
    return M


def Whirl3():
    """
    The rank-`3` whirl.

    The unique relaxation of `M(K_4)`. Self-dual but not identically self-dual.

    See [Oxl2011]_, p. 641.

    EXAMPLES::

        sage: W = matroids.catalog.Whirl3()
        sage: W.equals(W.dual())
        False
        sage: W.is_isomorphic(W.dual())
        True

    """
    M = Whirl(3)
    M.rename("Whirl(3): " + repr(M))
    return M


def Q6():
    """
    Return the matroid `Q_6`, represented over `GF(4)`.

    The matroid `Q_6` is a 6-element matroid of rank-3.
    It is representable over a field if and only if that field has at least
    four elements. It is the unique relaxation of the rank-3 whirl.
    See [Oxl2011]_, p. 641.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.Q6(); M                                       # needs sage.rings.finite_rings
        Q6: Quaternary matroid of rank 3 on 6 elements
        sage: setprint(M.hyperplanes())                                                 # needs sage.rings.finite_rings
        [{'a', 'b', 'd'}, {'a', 'c'}, {'a', 'e'}, {'a', 'f'}, {'b', 'c', 'e'},
         {'b', 'f'}, {'c', 'd'}, {'c', 'f'}, {'d', 'e'}, {'d', 'f'},
         {'e', 'f'}]
        sage: M.nonspanning_circuits() == M.noncospanning_cocircuits()                  # needs sage.rings.finite_rings
        False
    """
    F = GF(4, "x")
    x = F.gens()[0]
    A = Matrix(F, [[1, 0, 0, 1, 0, 1], [0, 1, 0, 1, 1, x], [0, 0, 1, 0, 1, 1]])
    M = QuaternaryMatroid(A, "abcdef")
    M.rename("Q6: " + repr(M))
    return M


def P6():
    """
    Return the matroid `P_6`, represented as circuit closures.

    The matroid `P_6` is a 6-element matroid of rank-3.
    It is representable over a field if and only if that field has at least
    five elements.
    It is the unique relaxation of `Q_6`.
    It is an excluded minor for the class of quaternary matroids.
    See [Oxl2011]_, p. 641.

    EXAMPLES::

        sage: M = matroids.catalog.P6(); M
        P6: Matroid of rank 3 on 6 elements with circuit-closures
        {2: {{'a', 'b', 'c'}}, 3: {{'a', 'b', 'c', 'd', 'e', 'f'}}}
        sage: len(set(M.nonspanning_circuits()).difference(M.nonbases())) == 0
        True
        sage: Matroid(matrix=random_matrix(GF(4, 'a'), ncols=5,                         # needs sage.rings.finite_rings
        ....:                                          nrows=5)).has_minor(M)
        False
        sage: M.is_valid()
        True

    """
    E = "abcdef"
    CC = {2: ["abc"], 3: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("P6: " + repr(M))
    return M


def U36():
    """
    The uniform matroid of rank `3` on `6` elements.

    Six points freely placed in the plane; the tipless free 3-spike.
    Identically self-dual.

    See [Oxl2011]_, p. 642.

    EXAMPLES::

        sage: U36 = matroids.catalog.U36()
        sage: Z = matroids.Spike(3, False)
        sage: U36.is_isomorphic(Z)
        True
        sage: U36.equals(U36.dual())
        True

    """
    M = Uniform(3, 6)
    M.rename("U(3, 6): " + repr(M))
    return M


def R6():
    """
    Return the matroid `R_6`, represented over `GF(3)`.

    The matroid `R_6` is a 6-element matroid of rank-3.
    It is representable over a field if and only if that field has at least
    three elements.
    It is isomorphic to the 2-sum of two copies of `U_{2, 4}`.
    See [Oxl2011]_, p. 642.

    EXAMPLES::

        sage: M = matroids.catalog.R6(); M
        R6: Ternary matroid of rank 3 on 6 elements, type 2+
        sage: M.equals(M.dual())
        True
        sage: M.is_connected()
        True
        sage: M.is_3connected()                                                         # needs sage.graphs
        False
    """
    A = Matrix(GF(3), [[1, 0, 0, 1, 1, 1], [0, 1, 0, 1, 2, 1], [0, 0, 1, 1, 0, 2]])
    M = TernaryMatroid(A, "abcdef")
    M.rename("R6: " + repr(M))
    return M


def Fano():
    r"""
    Return the Fano matroid, represented over `GF(2)`.

    The Fano matroid, or Fano plane, or `F_7`, is a 7-element matroid of
    rank-3.
    It is representable over a field if and only if that field has
    characteristic two.
    It is also the projective plane of order two, i.e. `\mathrm{PG}(2, 2)`.
    `F_7` is isomorphic to the unique `S(2, 3, 7)` and to the unique binary
    3-spike.

    See [Oxl2011]_, p. 643.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.Fano(); M
        Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
        sage: setprint(sorted(M.nonspanning_circuits()))
        [{'a', 'b', 'f'}, {'a', 'c', 'e'}, {'a', 'd', 'g'}, {'b', 'c', 'd'},
         {'b', 'e', 'g'}, {'c', 'f', 'g'}, {'d', 'e', 'f'}]
        sage: M.delete(M.groundset_list()[randrange(0,                                  # needs sage.graphs
        ....:                  7)]).is_isomorphic(matroids.CompleteGraphic(4))
        True
        sage: M.is_isomorphic(matroids.Z(3))
        True

    """
    A = Matrix(
        GF(2),
        [[1, 0, 0, 0, 1, 1, 1], [0, 1, 0, 1, 0, 1, 1], [0, 0, 1, 1, 1, 0, 1]]
    )
    M = BinaryMatroid(A, "abcdefg")
    M.rename("Fano: " + repr(M))
    return M


def FanoDual():
    """
    Return the dual of the Fano matroid.

    `F_7^*` is a `7`-element matroid of rank-`3`.

    See [Oxl2011]_, p. 643.

    EXAMPLES::

        sage: F7 = matroids.catalog.Fano()
        sage: F7D = matroids.catalog.FanoDual()
        sage: F7.is_isomorphic(F7D.dual())
        True

    """
    M = Fano().dual()
    M.rename("FanoDual: " + repr(M))
    return M


def NonFano():
    """
    Return the non-Fano matroid, represented over `GF(3)`.

    The non-Fano matroid, or `F_7^-`, is a 7-element matroid of rank-3.
    It is representable over a field if and only if that field has
    characteristic other than two.
    It is the unique relaxation of `F_7`. See [Oxl2011]_, p. 643-4.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.NonFano(); M
        NonFano: Ternary matroid of rank 3 on 7 elements, type 0-
        sage: setprint(M.nonbases())
        [{'a', 'b', 'f'}, {'a', 'c', 'e'}, {'a', 'd', 'g'}, {'b', 'c', 'd'},
         {'b', 'e', 'g'}, {'c', 'f', 'g'}]
        sage: M.delete('f').is_isomorphic(matroids.CompleteGraphic(4))                  # needs sage.graphs
        True
        sage: M.delete('g').is_isomorphic(matroids.CompleteGraphic(4))                  # needs sage.graphs
        False

    """
    A = Matrix(
        GF(3),
        [[1, 0, 0, 0, 1, 1, 1], [0, 1, 0, 1, 0, 1, 1], [0, 0, 1, 1, 1, 0, 1]]
    )
    M = TernaryMatroid(A, "abcdefg")
    M.rename("NonFano: " + repr(M))
    return M


def NonFanoDual():
    r"""
    Return the dual of the non-Fano matroid.

    `(F_7^-)^*` is a 7-element matroid of rank-3. Every single-element
    contraction of `(F_7^-)^*` is isomorphic to `M(K_4)` or `\mathcal{W}^3`.

    See [Oxl2011]_, p. 643-4.

    EXAMPLES::

        sage: M = matroids.catalog.NonFanoDual()
        sage: sorted(M.groundset())
        ['a', 'b', 'c', 'd', 'e', 'f', 'g']

    """
    M = NonFano().dual()
    M.rename("NonFanoDual: " + repr(M))
    return M


def O7():
    """
    Return the matroid `O_7`, represented over `GF(3)`.

    The matroid `O_7` is a 7-element matroid of rank-3.
    It is representable over a field if and only if that field has at least
    three elements.
    It is obtained by freely adding a point to any line of `M(K_4)`.
    See [Oxl2011]_, p. 644

    EXAMPLES::

        sage: M = matroids.catalog.O7(); M
        O7: Ternary matroid of rank 3 on 7 elements, type 0+
        sage: M.delete('e').is_isomorphic(matroids.CompleteGraphic(4))                  # needs sage.graphs
        True
        sage: M.tutte_polynomial()
        y^4 + x^3 + x*y^2 + 3*y^3 + 4*x^2 + 5*x*y + 5*y^2 + 4*x + 4*y

    """
    A = Matrix(
        GF(3),
        [[1, 0, 0, 1, 1, 1, 1], [0, 1, 0, 0, 1, 2, 2], [0, 0, 1, 1, 0, 1, 0]]
    )
    M = TernaryMatroid(A, "abcdefg")
    M.rename("O7: " + repr(M))
    return M


def P7():
    """
    Return the matroid `P_7`, represented over `GF(3)`.

    The matroid `P_7` is a 7-element matroid of rank-3.
    It is representable over a field if and only if that field has at least
    3 elements.
    It is one of two ternary 3-spikes, with the other being `F_7^-`.
    See [Oxl2011]_, p. 644.

    EXAMPLES::

        sage: M = matroids.catalog.P7(); M
        P7: Ternary matroid of rank 3 on 7 elements, type 1+
        sage: M.f_vector()
        [1, 7, 11, 1]
        sage: M.has_minor(matroids.CompleteGraphic(4))                                  # needs sage.graphs
        False
        sage: M.is_valid()
        True

    """
    A = Matrix(
        GF(3),
        [[1, 0, 0, 2, 1, 1, 0], [0, 1, 0, 1, 1, 0, 1], [0, 0, 1, 1, 0, 1, 1]]
    )
    M = TernaryMatroid(A, "abcdefg")
    M.rename("P7: " + repr(M))
    return M


def AG32():
    """
    Return the matroid `AG(3, 2)`.

    The binary affine cube. Isomorphic to the unique tipless binary 4-spike.

    EXAMPLES::

        sage: M = matroids.catalog.AG32()
        sage: M.is_valid()
        True
        sage: M.is_isomorphic(matroids.Z(4, False))
        True

    """
    M = AG(3, 2)
    M.rename("AG(3, 2): " + repr(M))
    return M


def AG32prime():
    """
    Return the matroid `AG(3, 2)'`, represented as circuit closures.

    The matroid `AG(3, 2)'` is a 8-element matroid of rank-4.
    It is a smallest non-representable matroid.
    It is the unique relaxation of `AG(3, 2)`. See [Oxl2011]_, p. 646.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.AG32prime(); M
        AG(3, 2)': Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
             {'a', 'c', 'd', 'f'}, {'a', 'c', 'e', 'g'}, {'a', 'd', 'g', 'h'},
             {'a', 'e', 'f', 'h'}, {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'},
             {'b', 'e', 'g', 'h'}, {'c', 'd', 'e', 'h'}, {'c', 'f', 'g', 'h'},
             {'d', 'e', 'f', 'g'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: M.contract('c').is_isomorphic(matroids.catalog.Fano())
        True
        sage: setprint(M.noncospanning_cocircuits())
        [{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
         {'a', 'c', 'd', 'f'}, {'a', 'd', 'g', 'h'}, {'a', 'e', 'f', 'h'},
         {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'}, {'b', 'd', 'f', 'h'},
         {'b', 'e', 'g', 'h'}, {'c', 'd', 'e', 'h'}, {'c', 'f', 'g', 'h'},
         {'d', 'e', 'f', 'g'}]
        sage: M.is_valid()                      # long time, needs sage.rings.finite_rings
        True

    """
    E = "abcdefgh"
    CC = {
        3: [
            "abfg", "bcdg", "defg", "cdeh", "aefh", "abch", "abed",
            "cfgh", "bcef", "adgh", "acdf", "begh", "aceg",
        ],
        4: [E],
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("AG(3, 2)': " + repr(M))
    return M


def R8():
    """
    Return the matroid `R_8`, represented over `GF(3)`.

    The matroid `R_8` is a 8-element matroid of rank-4.
    It is representable over a field if and only if the characteristic of that
    field is not two.
    It is the real affine cube. See [Oxl2011]_, p. 646.

    EXAMPLES::

        sage: M = matroids.catalog.R8(); M
        R8: Ternary matroid of rank 4 on 8 elements, type 0+
        sage: M.contract(M.groundset_list()[randrange(0,
        ....:            8)]).is_isomorphic(matroids.catalog.NonFano())
        True
        sage: M.equals(M.dual())
        True
        sage: M.has_minor(matroids.catalog.Fano())
        False
    """
    A = Matrix(
        GF(3),
        [
            [1, 0, 0, 0, 2, 1, 1, 1],
            [0, 1, 0, 0, 1, 2, 1, 1],
            [0, 0, 1, 0, 1, 1, 2, 1],
            [0, 0, 0, 1, 1, 1, 1, 2],
        ],
    )
    M = TernaryMatroid(A, "abcdefgh")
    M.rename("R8: " + repr(M))
    return M


def F8():
    """
    Return the matroid `F_8`, represented as circuit closures.

    The matroid `F_8` is a 8-element matroid of rank-4.
    It is a smallest non-representable matroid. See [Oxl2011]_, p. 647.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.catalog.F8(); M
        F8: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
             {'a', 'c', 'd', 'f'}, {'a', 'c', 'e', 'g'}, {'a', 'd', 'g', 'h'},
             {'a', 'e', 'f', 'h'}, {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'},
             {'c', 'd', 'e', 'h'}, {'c', 'f', 'g', 'h'}, {'d', 'e', 'f', 'g'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: D = get_nonisomorphic_matroids([M.contract(i)
        ....:                                         for i in M.groundset()])
        sage: len(D)
        3
        sage: [N.is_isomorphic(matroids.catalog.Fano()) for N in D]
        [...True...]
        sage: [N.is_isomorphic(matroids.catalog.NonFano()) for N in D]
        [...True...]
        sage: M.is_valid()                      # long time, needs sage.rings.finite_rings
        True

    """
    E = "abcdefgh"
    CC = {
        3: [
            "abfg", "bcdg", "defg", "cdeh",
            "aefh", "abch", "abed", "cfgh",
            "bcef", "adgh", "acdf", "aceg",
        ],
        4: [E],
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("F8: " + repr(M))
    return M


def Q8():
    """
    Return the matroid `Q_8`, represented as circuit closures.

    The matroid `Q_8` is a 8-element matroid of rank-4.
    It is a smallest non-representable matroid. See [Oxl2011]_, p. 647.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.Q8(); M
        Q8: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
             {'a', 'c', 'd', 'f'}, {'a', 'd', 'g', 'h'}, {'a', 'e', 'f', 'h'},
             {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'}, {'c', 'd', 'e', 'h'},
             {'c', 'f', 'g', 'h'}, {'d', 'e', 'f', 'g'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: setprint(M.flats(3))
        [{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
         {'a', 'c', 'd', 'f'}, {'a', 'c', 'e'}, {'a', 'c', 'g'},
         {'a', 'd', 'g', 'h'}, {'a', 'e', 'f', 'h'}, {'a', 'e', 'g'},
         {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'}, {'b', 'd', 'f'},
         {'b', 'd', 'h'}, {'b', 'e', 'g'}, {'b', 'e', 'h'}, {'b', 'f', 'h'},
         {'b', 'g', 'h'}, {'c', 'd', 'e', 'h'}, {'c', 'e', 'g'},
         {'c', 'f', 'g', 'h'}, {'d', 'e', 'f', 'g'}, {'d', 'f', 'h'},
         {'e', 'g', 'h'}]
        sage: M.is_valid() # long time
        True

    """
    E = "abcdefgh"
    CC = {
        3: [
            "abfg", "bcdg", "defg", "cdeh", "aefh", "abch",
            "abed", "cfgh", "bcef", "adgh", "acdf",
        ],
        4: [E],
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("Q8: " + repr(M))
    return M


def L8():
    """
    Return the matroid `L_8`, represented as circuit closures.

    The matroid `L_8` is a 8-element matroid of rank-4.
    It is representable over all fields with at least five elements.
    It is a cube, yet it is not a tipless spike. See [Oxl2011]_, p. 648.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.L8(); M
        L8: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'h'}, {'a', 'b', 'f', 'g'}, {'a', 'c', 'e', 'g'},
             {'a', 'e', 'f', 'h'}, {'b', 'c', 'd', 'g'}, {'b', 'd', 'f', 'h'},
             {'c', 'd', 'e', 'h'}, {'d', 'e', 'f', 'g'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: M.equals(M.dual())
        True
        sage: M.is_valid() # long time
        True

    """
    E = "abcdefgh"
    CC = {3: ["abfg", "bcdg", "defg", "cdeh", "aefh", "abch", "aceg", "bdfh"],
          4: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("L8: " + repr(M))
    return M


def S8():
    """
    Return the matroid `S_8`, represented over `GF(2)`.

    The matroid `S_8` is a 8-element matroid of rank-4.
    It is representable over a field if and only if that field has
    characteristic two.
    It is the unique deletion of a non-tip element from the binary 4-spike.
    See [Oxl2011]_, p. 648.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.catalog.S8(); M
        S8: Binary matroid of rank 4 on 8 elements, type (2, 0)
        sage: M.contract('d').is_isomorphic(matroids.catalog.Fano())
        True
        sage: M.delete('d').is_isomorphic(
        ....:                           matroids.catalog.Fano().dual())
        False
        sage: M.is_graphic()
        False
        sage: D = get_nonisomorphic_matroids(
        ....:       list(matroids.catalog.Fano().linear_coextensions(
        ....:                                                 cosimple=True)))
        sage: len(D)
        2
        sage: [N.is_isomorphic(M) for N in D]
        [...True...]

    """
    A = Matrix(
        GF(2),
        [
            [1, 0, 0, 0, 0, 1, 1, 1],
            [0, 1, 0, 0, 1, 0, 1, 1],
            [0, 0, 1, 0, 1, 1, 0, 1],
            [0, 0, 0, 1, 1, 1, 1, 1],
        ],
    )
    M = BinaryMatroid(A, "abcdefgh")
    M.rename("S8: " + repr(M))
    return M


def Vamos():
    """
    Return the Vamos matroid, represented as circuit closures.

    The Vamos matroid, or Vamos cube, or `V_8` is a 8-element matroid of
    rank-4.
    It violates Ingleton's condition for representability over a division
    ring.
    It is not algebraic. See [Oxl2011]_, p. 649.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.Vamos(); M
        Vamos: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
             {'c', 'd', 'e', 'f'}, {'e', 'f', 'g', 'h'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: setprint(M.nonbases())
        [{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
         {'c', 'd', 'e', 'f'}, {'e', 'f', 'g', 'h'}]
        sage: M.is_dependent(['c', 'd', 'g', 'h'])
        False
        sage: M.is_valid() # long time
        True

    """
    E = "abcdefgh"
    CC = {3: ["abcd", "abef", "cdef", "abgh", "efgh"], 4: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("Vamos: " + repr(M))
    return M


def T8():
    """
    Return the matroid `T_8`, represented over `GF(3)`.

    The matroid `T_8` is a 8-element matroid of rank-4.
    It is representable over a field if and only if that field has
    characteristic three.
    It is an excluded minor for the dyadic matroids. See [Oxl2011]_, p. 649.

    EXAMPLES::

        sage: M = matroids.catalog.T8(); M
        T8: Ternary matroid of rank 4 on 8 elements, type 0-
        sage: M.truncation().is_isomorphic(matroids.Uniform(3, 8))
        True
        sage: M.contract('e').is_isomorphic(matroids.catalog.P7())
        True
        sage: M.has_minor(matroids.Uniform(3, 8))
        False

    """
    A = Matrix(
        GF(3),
        [
            [1, 0, 0, 0, 0, 1, 1, 1],
            [0, 1, 0, 0, 1, 0, 1, 1],
            [0, 0, 1, 0, 1, 1, 0, 1],
            [0, 0, 0, 1, 1, 1, 1, 0],
        ],
    )
    M = TernaryMatroid(A, "abcdefgh")
    M.rename("T8: " + repr(M))
    return M


def J():
    """
    Return the matroid `J`, represented over `GF(3)`.

    The matroid `J` is a 8-element matroid of rank-4.
    It is representable over a field if and only if that field has at least
    three elements. See [Oxl2011]_, p. 650.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.J(); M
        J: Ternary matroid of rank 4 on 8 elements, type 0-
        sage: setprint(M.truncation().nonbases())
        [{'a', 'b', 'f'}, {'a', 'c', 'g'}, {'a', 'd', 'h'}]
        sage: M.is_isomorphic(M.dual())
        True
        sage: M.has_minor(matroids.CompleteGraphic(4))                                  # needs sage.graphs
        False
        sage: M.is_valid()
        True

    """
    A = Matrix(
        GF(3),
        [
            [1, 0, 0, 0, 0, 1, 1, 1],
            [0, 1, 0, 0, 1, 1, 0, 0],
            [0, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 1, 1, 0, 0, 1],
        ],
    )
    M = TernaryMatroid(A, "abcdefgh")
    M.rename("J: " + repr(M))
    return M


def P8():
    """
    Return the matroid `P_8`, represented over `GF(3)`.

    The matroid `P_8` is a 8-element matroid of rank-4.
    It is uniquely representable over all fields of characteristic other than
    two.
    It is an excluded minor for all fields of characteristic two with four or
    more elements. See [Oxl2011]_, p. 650.

    EXAMPLES::

        sage: M = matroids.catalog.P8(); M
        P8: Ternary matroid of rank 4 on 8 elements, type 2+
        sage: M.is_isomorphic(M.dual())
        True
        sage: Matroid(matrix=random_matrix(GF(4, 'a'), ncols=5,                         # needs sage.rings.finite_rings
        ....:                              nrows=5)).has_minor(M)
        False
        sage: M.bicycle_dimension()
        2

    """
    A = Matrix(
        GF(3),
        [
            [1, 0, 0, 0, 2, 1, 1, 0],
            [0, 1, 0, 0, 1, 1, 0, 1],
            [0, 0, 1, 0, 1, 0, 1, 1],
            [0, 0, 0, 1, 0, 1, 1, 2],
        ],
    )
    M = TernaryMatroid(A, "abcdefgh")
    M.rename("P8: " + repr(M))
    return M


def P8pp():
    """
    Return the matroid `P_8^=`, represented as circuit closures.

    The matroid `P_8^=` is a 8-element matroid of rank-4.
    It can be obtained from `P_8` by relaxing the unique pair of disjoint
    circuit-hyperplanes.
    It is an excluded minor for `GF(4)`-representability.
    It is representable over all fields with at least five elements.
    See [Oxl2011]_, p. 651.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.catalog.P8pp(); M
        P8'': Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'f', 'h'}, {'a', 'c', 'e', 'f'}, {'a', 'c', 'g', 'h'},
             {'a', 'd', 'e', 'g'}, {'b', 'c', 'e', 'g'}, {'b', 'd', 'e', 'h'},
             {'b', 'd', 'f', 'g'}, {'c', 'd', 'f', 'h'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: M.is_isomorphic(M.dual())
        True
        sage: len(get_nonisomorphic_matroids([M.contract(i)
        ....:                                        for i in M.groundset()]))
        1
        sage: M.is_valid() # long time
        True

    """
    E = "abcdefgh"
    CC = {3: ["abfh", "bceg", "cdfh", "adeg", "acef", "bdfg", "acgh", "bdeh"],
          4: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("P8'': " + repr(M))
    return M


def Wheel4():
    """
    Return the rank-4 wheel.

    Self-dual but not identically self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.Wheel4()
        sage: M.is_valid()
        True
        sage: M.equals(M.dual())
        False
        sage: M.is_isomorphic(M.dual())
        True

    """
    M = Wheel(4)
    M.rename("Wheel(4): " + repr(M))
    return M


def Whirl4():
    """
    Return the rank-4 whirl.

    Self-dual but not identically self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.Whirl4()
        sage: M.is_valid()
        True
        sage: M.equals(M.dual())
        False
        sage: M.is_isomorphic(M.dual())
        True

    """
    M = Wheel(4)
    M.rename("Whirl(4): " + repr(M))
    return M


def K33dual():
    """
    Return the matroid `M^*(K_{3, 3})`, represented over the regular partial
    field.

    The matroid `M^*(K_{3, 3})` is a 9-element matroid of rank-4.
    It is an excluded minor for the class of graphic matroids.
    It is the graft matroid of the 4-wheel with every vertex except the hub
    being coloured. See [Oxl2011]_, p. 652-3.

    EXAMPLES::

        sage: M = matroids.catalog.K33dual(); M                                  # needs sage.graphs
        M*(K3, 3): Regular matroid of rank 4 on 9 elements with 81 bases
        sage: any(N.is_3connected()                                                     # needs sage.graphs
        ....:     for N in M.linear_extensions(simple=True))
        False
        sage: M.is_valid()                      # long time, needs sage.graphs
        True

    """
    from sage.graphs.graph_generators import graphs

    E = "abcdefghi"
    G = graphs.CompleteBipartiteGraph(3, 3)
    M = Matroid(groundset=E, graph=G, regular=True)
    M = M.dual()
    M.rename("M*(K3, 3): " + repr(M))
    return M


def K33():
    r"""
    Return the graphic matroid `M(K_{3,3})`.

    `M(K_{3,3})` is an excluded minor for the class of cographic matroids.

    See [Oxl2011]_, p. 652-3.

    EXAMPLES::

        sage: M = matroids.catalog.K33()
        sage: M.is_valid()
        True

    """
    from sage.graphs.graph_generators import graphs

    E = "abcdefghi"
    G = graphs.CompleteBipartiteGraph(3, 3)
    M = Matroid(groundset=E, graph=G, regular=True)
    M.rename("M(K3, 3): " + repr(M))
    return M


def AG23():
    """
    Return the matroid `AG(2, 3)`.

    The ternary affine plane.

    See [Oxl2011]_, p. 653.

    EXAMPLES::

        sage: M = matroids.catalog.AG23()
        sage: M.is_valid()
        True

    """
    M = AG(2, 3)
    M.rename("AG(2, 3): " + repr(M))
    return M


def TernaryDowling3():
    """
    Return the matroid `Q_3(GF(3)^\times)`, represented over `GF(3)`.

    The matroid `Q_3(GF(3)^\times)` is a 9-element matroid of rank-3.
    It is the rank-3 ternary Dowling geometry.
    It is representable over a field if and only if that field does not have
    characteristic two. See [Oxl2011]_, p. 654.

    EXAMPLES::

        sage: M = matroids.catalog.TernaryDowling3(); M
        Q3(GF(3)x): Ternary matroid of rank 3 on 9 elements, type 0-
        sage: len(list(M.linear_subclasses()))
        72
        sage: M.fundamental_cycle('abc', 'd')
        {'a': 2, 'b': 1, 'd': 1}

    """
    A = Matrix(
        GF(3),
        [
            [1, 0, 0, 1, 1, 0, 0, 1, 1],
            [0, 1, 0, 2, 1, 1, 1, 0, 0],
            [0, 0, 1, 0, 0, 2, 1, 2, 1],
        ],
    )
    M = TernaryMatroid(A, "abcdefghi")
    M.rename("Q3(GF(3)x): " + repr(M))
    return M


def Pappus():
    """
    Return the Pappus matroid.

    The Pappus matroid is a 9-element matroid of rank-3.
    It is representable over a field if and only if that field either has 4
    elements or more than 7 elements.
    It is an excluded minor for the class of GF(5)-representable matroids.
    See [Oxl2011]_, p. 655.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.Pappus(); M
        Pappus: Matroid of rank 3 on 9 elements with circuit-closures
        {2: {{'a', 'b', 'c'}, {'a', 'e', 'i'}, {'a', 'f', 'h'},
             {'b', 'd', 'i'}, {'b', 'f', 'g'}, {'c', 'd', 'h'},
             {'c', 'e', 'g'}, {'d', 'e', 'f'}, {'g', 'h', 'i'}},
         3: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'}}}
        sage: setprint(M.nonspanning_circuits())
        [{'a', 'b', 'c'}, {'a', 'e', 'i'}, {'a', 'f', 'h'}, {'b', 'd', 'i'},
         {'b', 'f', 'g'}, {'c', 'd', 'h'}, {'c', 'e', 'g'}, {'d', 'e', 'f'},
         {'g', 'h', 'i'}]
        sage: M.is_dependent(['d', 'e', 'f'])
        True
        sage: M.is_valid() # long time
        True

    """
    E = "abcdefghi"
    CC = {2: ["abc", "def", "ceg", "bfg", "cdh", "afh", "bdi", "aei", "ghi"],
          3: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("Pappus: " + repr(M))
    return M


def NonPappus():
    """
    Return the non-Pappus matroid.

    The non-Pappus matroid is a 9-element matroid of rank-3.
    It is not representable over any commutative field.
    It is the unique relaxation of the Pappus matroid. See [Oxl2011]_, p. 655.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.NonPappus(); M
        NonPappus: Matroid of rank 3 on 9 elements with circuit-closures
        {2: {{'a', 'b', 'c'}, {'a', 'e', 'i'}, {'a', 'f', 'h'},
             {'b', 'd', 'i'}, {'b', 'f', 'g'}, {'c', 'd', 'h'},
             {'c', 'e', 'g'}, {'g', 'h', 'i'}},
         3: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'}}}
        sage: setprint(M.nonspanning_circuits())
        [{'a', 'b', 'c'}, {'a', 'e', 'i'}, {'a', 'f', 'h'}, {'b', 'd', 'i'},
         {'b', 'f', 'g'}, {'c', 'd', 'h'}, {'c', 'e', 'g'}, {'g', 'h', 'i'}]
        sage: M.is_dependent(['d', 'e', 'f'])
        False
        sage: M.is_valid() # long time
        True

    """
    E = "abcdefghi"
    CC = {2: ["abc", "ceg", "bfg", "cdh", "afh", "bdi", "aei", "ghi"],
          3: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("NonPappus: " + repr(M))
    return M


def K5():
    """
    Return the graphic matroid `M(K_5)`.

    `M(K_5)` is an excluded minor for the class of cographic matroids.

    See [Oxl2011]_, p. 656.

    EXAMPLES::

        sage: M = matroids.catalog.K5()
        sage: M.is_valid()
        True

    """
    M = CompleteGraphic(5)
    return M


def K5dual():
    """
    Return the matroid `M^*(K_5)`.

    `M^*(K_5)` is an excluded minor for the class of graphic matroids.

    See [Oxl2011]_, p. 656.

    EXAMPLES::

        sage: M = matroids.catalog.K5dual()
        sage: M.is_3connected()
        True

    """
    M = CompleteGraphic(5).dual()
    M.rename("M*(K5): " + repr(M))
    return M


def R10():
    """
    Return the matroid `R_{10}`, represented over the regular partial field.

    The matroid `R_{10}` is a 10-element regular matroid of rank-5.
    It is the unique splitter for the class of regular matroids.
    It is the graft matroid of `K_{3, 3}` in which every vertex is coloured.
    See [Oxl2011]_, p. 656.

    EXAMPLES::

        sage: M = matroids.catalog.R10(); M
        R10: Regular matroid of rank 5 on 10 elements with 162 bases
        sage: cct = []
        sage: for i in M.circuits():
        ....:      cct.append(len(i))
        sage: Set(cct)
        {4, 6}
        sage: M.equals(M.dual())
        False
        sage: M.is_isomorphic(M.dual())                                                 # needs sage.graphs
        True
        sage: M.is_valid()
        True

    Check the splitter property::

        sage: matroids.catalog.R10().linear_extensions(simple=True)
        []

    """
    A = Matrix(
        ZZ,
        [
            [1, 0, 0, 0, 0, -1, 1, 0, 0, 1],
            [0, 1, 0, 0, 0, 1, -1, 1, 0, 0],
            [0, 0, 1, 0, 0, 0, 1, -1, 1, 0],
            [0, 0, 0, 1, 0, 0, 0, 1, -1, 1],
            [0, 0, 0, 0, 1, 1, 0, 0, 1, -1],
        ],
    )
    M = RegularMatroid(A, "abcdefghij")
    M.rename("R10: " + repr(M))
    return M


# NonDesargues


def R12():
    """
    Return the matroid `R_{12}`, represented over the regular partial field.

    The matroid `R_{12}` is a 12-element regular matroid of rank-6.
    It induces a 3-separation in its 3-connected majors within the class of
    regular matroids.
    An excluded minor for the class of graphic or cographic matroids.
    See [Oxl2011]_, p. 657.

    EXAMPLES::

        sage: M = matroids.catalog.R12(); M
        R12: Regular matroid of rank 6 on 12 elements with 441 bases
        sage: M.equals(M.dual())
        False
        sage: M.is_isomorphic(M.dual())                                                 # needs sage.graphs
        True
        sage: M.is_valid() # long time
        True

    """
    A = Matrix(
        ZZ,
        [
            [1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0],
            [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
            [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1],
            [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, -1, -1],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -1, -1],
        ],
    )
    M = RegularMatroid(A, "abcdefghijkl")
    M.rename("R12: " + repr(M))
    return M


def ExtendedTernaryGolayCode():
    """
    Return the matroid of the extended ternary Golay code.

    This is the unique Steiner system `S(5, 6, 12)`.

    See [Oxl2011]_, p. 658, and
    :class:`GolayCode <sage.coding.golay_code.GolayCode>`

    EXAMPLES::

        sage: M = matroids.catalog.ExtendedTernaryGolayCode()
        sage: C = LinearCode(M.representation())
        sage: C.is_permutation_equivalent(codes.GolayCode(GF(3)))       # long time, needs sage.rings.finite_rings
        True
        sage: M.is_valid()
        True

    `S(5, 6, 12)` is an identically self-dual matroid::

        sage: M.equals(M.dual())
        True

    Every contraction of three elements is isomorphic to `AG(2, 3)`; every
    contraction of two elements and deletion of two elements is isomorphic to
    `P8`::

        sage: import random, itertools
        sage: C = list(itertools.combinations(M.groundset(), 3))
        sage: elements = random.choice(C)
        sage: N = M.contract(elements)
        sage: N.is_isomorphic(matroids.catalog.AG23())
        True
        sage: C = list(itertools.combinations(M.groundset(), 4))
        sage: elements = list(random.choice(C))
        sage: random.shuffle(elements)
        sage: N = M.contract(elements[:2]).delete(elements[2:4])
        sage: N.is_isomorphic(matroids.catalog.P8())
        True

    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 0],
        [0, 1, 0, 0, 0, 0, 1, 1, 2, 1, 0, 2],
        [0, 0, 1, 0, 0, 0, 1, 2, 1, 0, 1, 2],
        [0, 0, 0, 1, 0, 0, 1, 2, 0, 1, 2, 1],
        [0, 0, 0, 0, 1, 0, 1, 0, 2, 2, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1]
    ])
    M = TernaryMatroid(A, "abcdefghijkl")
    M.rename('Extended Ternary Golay Code: ' + repr(M))
    return M


def T12():
    """
    Return the matroid `T_{12}`.

    The edges of the Petersen graph can be labeled by the 4-circuits of
    `T_{12}` so that two edges are adjacent if and only if the corresponding
    4-circuits overlap in exactly two elements.
    Relaxing a circuit-hyperplane yields an excluded minor for the class of
    matroids that are either binary or ternary. See [Oxl2011]_, p. 658.

    EXAMPLES::

        sage: M = matroids.catalog.T12()
        sage: M
        T12: Binary matroid of rank 6 on 12 elements, type (2, None)
        sage: M.is_valid()
        True

    """
    A = Matrix(
        GF(2),
        [
            [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
            [0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1],
            [0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0],
            [0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1],
        ],
    )
    M = BinaryMatroid(A, "abcdefghijkl")
    M.rename("T12: " + repr(M))
    return M


def PG23():
    """
    Return the matroid `PG23`.

    The second smallest projective plane. Not graphic, not cographic, not
    regular, not near-regular.

    See [Oxl2011]_, p. 659.

    EXAMPLES::

        sage: M = matroids.catalog.PG23()
        sage: M.is_3connected()
        True

    """
    M = PG(2, 3)
    M.rename("PG(2, 3): " + repr(M))
    return M


def Wheel(n, field=None, ring=None):
    r"""
    Return the rank-`n` wheel.

    INPUT:

    - ``n`` -- a positive integer. The rank of the desired matroid.
    - ``ring`` -- any ring. If provided, output will be a linear matroid
      over the ring or field ``ring``. If the ring is `\ZZ`, then output
      will be a regular matroid.
    - ``field`` -- any field. Same as ``ring``, but only fields are allowed.

    OUTPUT:

    The rank-`n` wheel matroid, represented as a regular matroid.

    See [Oxl2011]_, p. 659.

    EXAMPLES::

        sage: M = matroids.Wheel(5); M
        Wheel(5): Regular matroid of rank 5 on 10 elements with 121 bases
        sage: M.tutte_polynomial()
        x^5 + y^5 + 5*x^4 + 5*x^3*y + 5*x^2*y^2 + 5*x*y^3 + 5*y^4 + 10*x^3 +
        15*x^2*y + 15*x*y^2 + 10*y^3 + 10*x^2 + 16*x*y + 10*y^2 + 4*x + 4*y
        sage: M.is_valid()
        True
        sage: M = matroids.Wheel(3)
        sage: M.is_isomorphic(matroids.CompleteGraphic(4))                              # needs sage.graphs
        True
        sage: M.is_isomorphic(matroids.Wheel(3, field=GF(3)))
        True
        sage: M = matroids.Wheel(3, field=GF(3)); M
        Wheel(3): Ternary matroid of rank 3 on 6 elements, type 0+
    """
    base_ring = ZZ
    if field is not None and ring is not None:
        raise ValueError("only one of ring and field can be specified.")
    if field is not None:
        base_ring = field
        try:
            if not base_ring.is_field():
                raise TypeError("specified ``field`` is not a field.")
        except AttributeError:
            raise TypeError("specified ``field`` is not a field.")
    if ring is not None:
        base_ring = ring
    A = Matrix(base_ring, n, 2 * n, sparse=True)
    for i in range(n):
        A[i, i] = 1
        A[i, n + i] = 1
        if i != 0:
            A[i, i + n - 1] = -1
        else:
            A[i, 2 * n - 1] = -1
    if base_ring is ZZ:
        M = RegularMatroid(A)
    else:
        M = Matroid(A)
    M.rename("Wheel(" + str(n) + "): " + repr(M))
    return M


def Whirl(n):
    """
    Return the rank-`n` whirl.

    INPUT:

    - ``n`` -- a positive integer. The rank of the desired matroid.

    OUTPUT:

    The rank-`n` whirl matroid, represented as a ternary matroid.

    The whirl is the unique relaxation of the wheel. See [Oxl2011]_, p. 659.

    EXAMPLES::

        sage: M = matroids.Whirl(5); M
        Whirl(5): Ternary matroid of rank 5 on 10 elements, type 0-
        sage: M.is_valid()
        True
        sage: M.tutte_polynomial()
        x^5 + y^5 + 5*x^4 + 5*x^3*y + 5*x^2*y^2 + 5*x*y^3 + 5*y^4 + 10*x^3 +
        15*x^2*y + 15*x*y^2 + 10*y^3 + 10*x^2 + 15*x*y + 10*y^2 + 5*x + 5*y
        sage: M.is_isomorphic(matroids.Wheel(5))
        False
        sage: M = matroids.Whirl(3)
        sage: M.is_isomorphic(matroids.CompleteGraphic(4))                              # needs sage.graphs
        False

    .. TODO::

        Optional arguments ``ring`` and ``x``, such that the resulting matroid
        is represented over ``ring`` by a reduced matrix like
        ``[-1  0  x]``
        ``[ 1 -1  0]``
        ``[ 0  1 -1]``
    """
    A = Matrix(GF(3), n, 2 * n, sparse=True)
    for i in range(n):
        A[i, i] = 1
        A[i, n + i] = 1
        if i != 0:
            A[i, i + n - 1] = -1
        else:
            A[i, 2 * n - 1] = 1
    M = TernaryMatroid(A)
    M.rename("Whirl(" + str(n) + "): " + repr(M))
    return M


def Uniform(r, n):
    """
    Return the uniform matroid of rank `r` on `n` elements.

    INPUT:

    - ``r`` -- a nonnegative integer. The rank of the uniform matroid.
    - ``n`` -- a nonnegative integer. The number of elements of the uniform
      matroid.

    OUTPUT:

    The uniform matroid `U_{r,n}`.

    All subsets of size `r` or less are independent; all larger subsets are
    dependent. Representable when the field is sufficiently large. The precise
    bound is the subject of the MDS conjecture from coding theory.
    See [Oxl2011]_, p. 660.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.Uniform(2, 5); M
        U(2, 5): Matroid of rank 2 on 5 elements with circuit-closures
        {2: {{0, 1, 2, 3, 4}}}
        sage: M.dual().is_isomorphic(matroids.Uniform(3, 5))
        True
        sage: setprint(M.hyperplanes())
        [{0}, {1}, {2}, {3}, {4}]
        sage: M.has_line_minor(6)
        False
        sage: M.is_valid()
        True

    Check that bug :trac:`15292` was fixed::

        sage: M = matroids.Uniform(4,4)
        sage: len(M.circuit_closures())
        0
    """
    E = list(range(n))
    if r < n:
        CC = {r: [E]}
    else:
        CC = {}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename("U(" + str(r) + ", " + str(n) + "): " + repr(M))
    return M


def PG(n, q, x=None):
    """
    Return the projective geometry of dimension ``n`` over the finite field
    of order ``q``.

    INPUT:

    - ``n`` -- a positive integer. The dimension of the projective space. This
      is one less than the rank of the resulting matroid.
    - ``q`` -- a positive integer that is a prime power. The order of the
      finite field.
    - ``x`` -- (default: ``None``) a string. The name of the generator of a
      non-prime field, used for non-prime fields. If not supplied, ``'x'`` is
      used.

    OUTPUT:

    A linear matroid whose elements are the points of `PG(n, q)`.

    EXAMPLES::

        sage: M = matroids.PG(2, 2)
        sage: M.is_isomorphic(matroids.catalog.Fano())
        True
        sage: matroids.PG(5, 4, 'z').size() == (4^6 - 1) / (4 - 1)                      # needs sage.rings.finite_rings
        True
        sage: M = matroids.PG(4, 7); M
        PG(4, 7): Linear matroid of rank 5 on 2801 elements represented over
        the Finite Field of size 7
    """
    if x is None:
        x = "x"
    F = GF(q, x)
    P = ProjectiveSpace(n, F)
    A = Matrix(F, [list(p) for p in list(P)]).transpose()
    M = Matroid(A)
    M.rename("PG(" + str(n) + ", " + str(q) + "): " + repr(M))
    return M


def AG(n, q, x=None):
    r"""
    Return the affine geometry of dimension ``n`` over the finite field of
    order ``q``.

    INPUT:

    - ``n`` -- a positive integer. The dimension of the projective space. This
      is one less than the rank of the resulting matroid.
    - ``q`` -- a positive integer that is a prime power. The order of the
      finite field.
    - ``x`` -- (default: ``None``) a string. The name of the generator of a
      non-prime field, used for non-prime fields. If not supplied, ``'x'`` is
      used.

    OUTPUT:

    A linear matroid whose elements are the points of `AG(n, q)`.

    The affine geometry can be obtained from the projective geometry by
    removing a hyperplane.

    EXAMPLES::

        sage: M = matroids.AG(2, 3).delete(8)
        sage: M.is_isomorphic(matroids.catalog.AG23minus())
        True
        sage: matroids.AG(5, 4, 'z').size() == ((4 ^ 6 - 1) / (4 - 1) -                 # needs sage.rings.finite_rings
        ....:                                             (4 ^ 5 - 1)/(4 - 1))
        True
        sage: M = matroids.AG(4, 2); M
        AG(4, 2): Binary matroid of rank 5 on 16 elements, type (5, 0)

    """
    if x is None:
        x = "x"
    F = GF(q, x)
    P = ProjectiveSpace(n, F)
    A = Matrix(
        F, [list(p) for p in list(P) if not list(p)[0] == 0]
    ).transpose()
    M = Matroid(A)
    M.rename("AG(" + str(n) + ", " + str(q) + "): " + repr(M))
    return M


def Z(r, t=True):
    r"""
    Return the unique rank-`r` binary spike.

    Defined for all `r \ge 3`.

    See p. 661-2 of [Oxl2011]_.

    INPUT:

    - ``r`` -- an integer (`r \ge 3`); the rank of the spike
    - ``t`` -- a Boolean (default: ``True``); whether the spike is tipped

    OUTPUT:

    a matroid; the unique rank-`r` binary spike (tipped or tipless)

    EXAMPLES::

        sage: import random

    It holds that `Z_3 \setminus e \cong M(K4)`, for all `e`::

        sage: Z3 = matroids.Z(3)
        sage: E = sorted(Z3.groundset())
        sage: e = random.choice(E)
        sage: Z3.delete(e).is_isomorphic(matroids.catalog.K4())
        True

    `Z_3 \cong F_7`::

        sage: F7 = matroids.catalog.Fano()
        sage: Z3.is_isomorphic(F7)
        True

    `Z_4 \setminus t \cong AG(3, 2)`::

        sage: Z4 = matroids.Z(4, False)
        sage: Z4.is_isomorphic(matroids.catalog.AG32())
        True

    The tipless binary spike is self-dual; it is identically self-dual if and
    only if r is even::

        sage: r = random.choice(range(3, 8))
        sage: Z = matroids.Z(r, False)
        sage: Z.is_isomorphic(Z.dual())
        True
        sage: Z.equals(Z.dual()) != (r % 2 == 1)  # XOR
        True

    """
    from sage.matrix.special import identity_matrix, ones_matrix
    Id = Matrix(GF(2), identity_matrix(r))
    J = Matrix(GF(2), ones_matrix(r))
    tip = Matrix(GF(2), ones_matrix(r, 1))
    A = Id.augment(J-Id).augment(tip)

    M = Matroid(A)
    # X = ["x"+str(i) for i in range(1, r+1)]
    # Y = ["y"+str(i) for i in range(1, r+1)]
    if t:
        # M = M.relabel(X+Y+["t"])
        M.rename("Z_" + str(r) + ": " + repr(M))
    else:
        M = M.delete(2*r)
        # M = M.relabel(X+Y)
        M.rename("Z_" + str(r) + "\\t: " + repr(M))
    return M


def Spike(r, t=True, C3=[]):
    r"""
    Return a rank-r spike.

    Defined for all `r \ge 3`; a rank-r spike with tip `t` and legs `L_1,
    L_2, \ldots, L_r`. Deleting `t` gives a tipless rank-`r` spike.

    The groundset is `E = \{t, x_1, x_2, \ldots, x_r, y_1, y_2, \ldots,
    y_r\}` with `r(E) = r`.

    Let `L_i = \{t, x_i , y_i\}`. The non-spanning circuits are `\{L_1 , L_2 ,
    \ldots, L_r\}`, all sets of the form `(L_i \cup L_j) \setminus t` for `1
    \le i < j \le r`, and some (possibly empty) collection `C_3` of sets of
    the form `\{z_1, z_2, \ldots, z_r\}` where `z_i \in \{x_i, y_i\}` for all
    `i`, and no two members of `C_3` have more than `r-2` common elements.

    See p. 662 of [Oxl2011]_.

    INPUT:

    - ``r`` -- an integer (`r \ge 3`); the rank of the spike
    - ``t`` -- a boolean (default: ``True``); whether the spike is tipped
    - ``C3`` -- a list (default: ``[]``); a list of extra nonspanning circuits.
      The default (i.e. the empty list) results in a free `r`-spike

    OUTPUT:

    a matroid; a rank-`r` spike (tipped or tipless)

    EXAMPLES::

        sage: M = matroids.Spike(3, False)
        sage: M.is_isomorphic(matroids.Uniform(3, 6))
        True
        sage: import random
        sage: r = random.choice(range(3, 8))
        sage: M = matroids.Spike(r)
        sage: M.is_3connected()
        True

    Each of `F_7`, `F_7^-`, and `P_7`, is a 3-spike. After inspection of the
    nonspanning circuits of these matroids, it becomes clear that they indeed
    constitute tipped 3-spikes. This can be verified by using an appropriate
    choice of extra circuits in `C_3`::

        sage: M = matroids.Spike(3, C3=[['x1', 'x2', 'y3'],
        ....:                           ['x1', 'x3', 'y2'],
        ....:                           ['x2', 'x3', 'y1'],
        ....:                           ['y1', 'y2', 'y3']])
        sage: M.is_isomorphic(matroids.catalog.Fano())
        True
        sage: M = matroids.Spike(3, C3=[['x1', 'x2', 'x3'],
        ....:                           ['x1', 'y2', 'y3'],
        ....:                           ['x2', 'y1', 'y3']])
        sage: M.is_isomorphic(matroids.catalog.NonFano())
        True
        sage: M = matroids.Spike(3, C3=[['x1', 'x2', 'y3'],
        ....:                           ['x3', 'y1', 'y2']])
        sage: M.is_isomorphic(matroids.catalog.P7())
        True

    Deleting any element gives a self-dual matroid. The tipless free spike
    (i.e., when `C_3` is empty) is identically self-dual::

        sage: M = matroids.Spike(6)
        sage: e = random.choice(list(M.groundset()))
        sage: Minor = M.delete(e)
        sage: Minor.is_isomorphic(Minor.dual())
        True
        sage: r = random.choice(range(3, 8))
        sage: M = matroids.Spike(r, False)
        sage: M.equals(M.dual())
        True

    """
    if not (r >= 3):
        raise ValueError("The r-spike is defined for r >= 3.")

    E = ["t"]
    X, Y = [], []
    for i in range(1, r + 1):
        X.append("x" + str(i))
        Y.append("y" + str(i))
    E += X
    E += Y

    for S in C3:
        for xy in S:
            if xy not in X+Y:
                raise ValueError(
                    "The sets in C3 must contain elements x_i and y_i only"
                )
        for T in C3:
            if S != T and len(set(S).intersection(set(T))) > r - 2:
                raise ValueError(
                    "Every pair of sets in C3 must not have more than r - 2 "
                    + "common elements"
                )

    NSC = []  # nonspanning_circuits
    NSC += C3
    for i in range(1, len(X)+1):
        NSC += [["t", "x"+str(i), "y"+str(i)]]
        for j in range(i+1, len(Y)+1):
            NSC += [["x"+str(i), "y"+str(i), "x"+str(j), "y"+str(j)]]

    import itertools
    B = []  # bases
    for b in itertools.combinations(E, r):
        flag = True
        for C in NSC:
            if set(b) >= set(C):
                flag = False
                break
        if flag:
            B += [list(b)]

    M = Matroid(groundset=E, bases=B)
    free = "Free " if C3 == [] else ""
    tip = "" if t else "\\t"
    M = M if t else M.delete('t')
    M.rename(free + str(r) + "-spike" + tip + ": " + repr(M))
    return M


# Q_r(A)


def Theta(n):
    r"""
    Return the matroid `\Theta_n`.

    Defined for all `n \ge 2`. `\Theta_2 \cong U_{1,2} \bigoplus U_{1,2}` and
    `\Theta_3 \cong M(K_4)`.

    See [Oxl2011]_, p. 663-4.

    INPUT:

    - ``n`` -- an integer (`n \ge 2`); the rank of the matroid

    OUTPUT:

    a matroid (`\Theta_n`)

    EXAMPLES::

        sage: M = matroids.Theta(2)
        sage: U12 = matroids.Uniform(1, 2)
        sage: U = U12.direct_sum(U12)
        sage: M.is_isomorphic(U)
        True
        sage: M = matroids.Theta(3)
        sage: M.is_isomorphic(matroids.catalog.K4())
        True

    `\Theta_n` is self-dual; identically self-dual if and only if `n = 2`::

        sage: M = matroids.Theta(2)
        sage: M.equals(M.dual())
        True
        sage: import random
        sage: n = random.choice(range(3, 10))
        sage: M = matroids.Theta(n)
        sage: M.equals(M.dual())
        False
        sage: M.is_isomorphic(M.dual())
        True

    """
    X = ["x"+str(i) for i in range(1, n+1)]
    Y = ["y"+str(i) for i in range(1, n+1)]
    E = X + Y

    import itertools
    C = []
    C += list(itertools.combinations(X, 3))
    for i in range(1, n+1):
        Yi = [Y[j] for j in range(len(Y)) if j != i-1]
        C += [Yi + ["x"+str(i)]]

    for u in range(1, n+1):
        for s in range(1, n+1):
            for t in range(1, n+1):
                if u != s and u != t and s != t:
                    Yu = [Y[i] for i in range(len(Y)) if i != u-1]
                    C += [Yu + ["x"+str(s)] + ["x"+str(t)]]

    M = Matroid(groundset=E, circuits=C)
    M.rename("Theta_" + str(n) + ": " + repr(M))
    return M


def Psi(r):
    r"""
    Return the matroid `\Psi_r`.

    The rank-`r` free swirl; defined for all `r \ge 3`.

    See [Oxl2011]_, p. 664.

    INPUT:

    - ``r`` -- an integer (`r \ge 3`); the rank of the matroid

    OUTPUT:

    a matroid (`\Psi_r`)

    EXAMPLES:

    The matroid `\Psi_r` is `3`-connected but, for all `r \ge 4`, not
    `4`-connected::

        sage: M = matroids.Psi(3)
        sage: M.is_4connected()
        True
        sage: import random
        sage: r = random.choice(range(4, 10))
        sage: M = matroids.Psi(r)
        sage: M.is_4connected()
        False

    `\Psi_3 \cong U_{3, 6}`::

        sage: M = matroids.Psi(3)
        sage: M.is_isomorphic(matroids.catalog.U36())
        True

    It is identically self-dual::

        sage: r = random.choice(range(3, 10))
        sage: M = matroids.Psi(r)
        sage: M.equals(M.dual())
        True

    """
    A = ["a"+str(i) for i in range(0, r)]
    B = ["b"+str(i) for i in range(0, r)]
    E = A + B

    def generate_binary_strings(bit_count):
        binary_strings = []

        def genbin(n, bs=""):
            if len(bs) == n:
                binary_strings.append(bs)
            else:
                genbin(n, bs + "a")
                genbin(n, bs + "b")

        genbin(bit_count)
        return binary_strings

    NSC = []  # nonspanning circuits
    for i in range(0, r):
        for k in range(1, r-2):
            I0 = ["a"+str(i), "b"+str(i)]
            IK = ["a"+str((i+k) % r), "b"+str((i+k) % r)]
            for AB in generate_binary_strings(k-1):
                C = []
                C += I0 + IK
                j = 1
                for z in AB:
                    C += [z+str((i+j) % r)]
                    j += 1
                NSC += [C]

    import itertools
    B = []  # bases
    for b in itertools.combinations(E, r):
        flag = True
        for C in NSC:
            if set(b) >= set(C):
                flag = False
                break
        if flag:
            B += [list(b)]

    M = Matroid(groundset=E, bases=B)
    M.rename("Psi_" + str(r) + ": " + repr(M))
    return M
