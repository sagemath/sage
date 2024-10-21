# sage_setup: distribution = sagemath-modules
r"""
Database of matroids

This module contains the implementation and documentation for all matroids in
the database, accessible through :mod:`matroids. <sage.matroids.>` and
:mod:`matroids.catalog. <sage.matroids.catalog>` (type those lines followed by
:kbd:`Tab` for a list).

AUTHORS:

- Michael Welsh, Stefan van Zwam (2013-04-01): initial version
- Giorgos Mousa, Andreas Triantafyllos (2023-12-08): more matroids

REFERENCES:

For more information on the matroids that belong to Oxley's matroid collection,
as well as for most parametrized matroids, see [Oxl2011]_.

For more information on the matroids that belong to Brettell's matroid
collection, see the associated publications, [Bre2023]_ and [BP2023]_.

.. NOTE::

    The grouping of the matroids into collections can be viewed in the
    :mod:`catalog of matroids <sage.matroids.matroids_catalog>`.
"""

# ****************************************************************************
#       Copyright (C) 2013 Michael Welsh <michael@welsh.co.nz>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#       Copyright (C) 2023 Giorgos Mousa <gmousa@proton.me>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.constructor import Matrix
from sage.matroids.constructor import Matroid
from sage.matroids.linear_matroid import (
    RegularMatroid,
    BinaryMatroid,
    TernaryMatroid,
    QuaternaryMatroid
)
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.schemes.projective.projective_space import ProjectiveSpace


# **************************** #
#  Oxley's matroid collection  #
#                              #
# **************************** #


# The order is the same as in Oxley.


def U24(groundset='abcd'):
    r"""
    The uniform matroid of rank `2` on `4` elements.

    The `4`-point line; isomorphic to `\mathcal{W}^2` , the rank-`2` whirl.
    The unique excluded minor for the class of binary matroids.

    EXAMPLES::

        sage: M = matroids.catalog.U24(); M
        U(2, 4): Matroid of rank 2 on 4 elements with circuit-closures
        {2: {{'a', 'b', 'c', 'd'}}}
        sage: N = matroids.Uniform(2, 4)
        sage: M.is_isomorphic(N)
        True
        sage: M.automorphism_group().structure_description()
        'S4'

    `U_{2,4}` is isomorphic to `\mathcal{W}^2`::

        sage: W2 = matroids.Whirl(2)
        sage: W2.is_isomorphic(M)
        True

    identically self-dual::

        sage: M.equals(M.dual())
        True

    and `3`-connected::

        sage: M.is_3connected()
        True

    REFERENCES:

    [Oxl2011]_, p. 639.

    TESTS::

        sage: M = matroids.catalog.U24(range(4))
        sage: sorted(M.groundset())
        [0, 1, 2, 3]
        sage: M.is_isomorphic(matroids.catalog.U24())
        True
    """
    M = Uniform(2, 4)
    M = _rename_and_relabel(M, "U(2, 4)", groundset)
    return M


def U25(groundset='abcde'):
    """
    The uniform matroid of rank `2` on `5` elements.

    `U_{2,5}` is the `5`-point line. Dual to `U_{3,5}`.

    EXAMPLES::

        sage: U25 = matroids.catalog.U25(); U25
        U(2, 5): Matroid of rank 2 on 5 elements with circuit-closures
        {2: {{'a', 'b', 'c', 'd', 'e'}}}
        sage: U25.is_graphic() or U25.is_regular()
        False
        sage: U35 = matroids.catalog.U35()
        sage: U25.is_isomorphic(U35.dual())
        True

    REFERENCES:

    [Oxl2011]_, p. 640.

    TESTS::

        sage: M = matroids.catalog.U25(range(5))
        sage: sorted(M.groundset())
        [0, 1, 2, 3, 4]
    """
    M = Uniform(2, 5)
    M = _rename_and_relabel(M, "U(2, 5)", groundset)
    return M


def U35(groundset='abcde'):
    """
    The uniform matroid of rank `3` on `5` elements.

    `U_{3,5}` is five points freely placed in the plane. Dual to `U_{2,5}`.

    EXAMPLES::

        sage: U35 = matroids.catalog.U35(); U35
        U(3, 5): Matroid of rank 3 on 5 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'd', 'e'}}}
        sage: U35.is_graphic() or U35.is_regular()
        False
        sage: U25 = matroids.catalog.U25()
        sage: U35.is_isomorphic(U25.dual())
        True

    REFERENCES:

    [Oxl2011]_, p. 640.

    TESTS::

        sage: M = matroids.catalog.U35(range(5))
        sage: sorted(M.groundset())
        [0, 1, 2, 3, 4]
    """
    M = Uniform(3, 5)
    M = _rename_and_relabel(M, "U(3, 5)", groundset)
    return M


def K4(groundset='abcdef'):
    r"""
    The graphic matroid of the complete graph `K_4`.

    EXAMPLES::

        sage: M = matroids.catalog.K4(); M
        M(K4): Graphic matroid of rank 3 on 6 elements
        sage: M.is_graphic()
        True

    `M(K_4)` is isomorphic to `M(\mathcal{W}_3)`, the rank-`3` wheel::

        sage: W3 = matroids.Wheel(3)
        sage: M.is_isomorphic(W3)
        True

    and to the tipless binary `3`-spike::

        sage: Z = matroids.Z(3, False)
        sage: M.is_isomorphic(Z)
        True

    It has a transitive automorphism group::

        sage: M.automorphism_group().is_transitive()
        True

    REFERENCES:

    [Oxl2011]_, p. 640.

    TESTS::

        sage: M = matroids.catalog.K4(range(6))
        sage: sorted(M.groundset())
        [0, 1, 2, 3, 4, 5]
    """
    M = CompleteGraphic(4)
    M = _rename_and_relabel(M, "M(K4)", groundset)
    return M


def Whirl3(groundset='abcdef'):
    r"""
    The rank-`3` whirl.

    The unique relaxation of `M(K_4)`. Self-dual but not identically self-dual.

    EXAMPLES::

        sage: W = matroids.catalog.Whirl3(); W
        Whirl(3): Ternary matroid of rank 3 on 6 elements, type 0-
        sage: W.equals(W.dual())
        False
        sage: W.is_isomorphic(W.dual())
        True
        sage: W.automorphism_group().is_transitive()
        False

    For all elements `e`, neither `\mathcal{W}_3 \setminus \{e\}` nor `\mathcal{W}_3 / \{e\}`
    is `3`-connected::

        sage: import random
        sage: e = random.choice(list(W.groundset()))
        sage: W.delete(e).is_3connected()
        False
        sage: W.contract(e).is_3connected()
        False

    REFERENCES:

    [Oxl2011]_, p. 641.

    TESTS::

        sage: W = matroids.catalog.Whirl3(range(6))
        sage: sorted(W.groundset())
        [0, 1, 2, 3, 4, 5]
    """
    M = Whirl(3)
    M = _rename_and_relabel(M, "Whirl(3)", groundset)
    return M


def Q6(groundset='abcdef'):
    """
    Return the matroid `Q_6`, represented over `GF(4)`.

    The matroid `Q_6` is a `6`-element matroid of rank-`3`.
    It is representable over a field if and only if that field has at least
    four elements. It is the unique relaxation of the rank-`3` whirl.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.Q6(); M
        Q6: Quaternary matroid of rank 3 on 6 elements
        sage: setprint(M.hyperplanes())
        [{'a', 'b', 'd'}, {'a', 'c'}, {'a', 'e'}, {'a', 'f'}, {'b', 'c', 'e'},
         {'b', 'f'}, {'c', 'd'}, {'c', 'f'}, {'d', 'e'}, {'d', 'f'},
         {'e', 'f'}]
        sage: M.nonspanning_circuits() == M.noncospanning_cocircuits()
        False
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 641.
    """
    F = GF(4, 'x')
    x = F.gens()[0]
    A = Matrix(F, [[1, 0, 0, 1, 0, 1], [0, 1, 0, 1, 1, x], [0, 0, 1, 0, 1, 1]])
    M = QuaternaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "Q6")
    return M


def P6(groundset=None):
    """
    Return the matroid `P_6`, represented as circuit closures.

    The matroid `P_6` is a `6`-element matroid of rank-`3`. It is
    representable over a field if and only if that field has at least five
    elements. It is the unique relaxation of `Q_6`. It is an excluded minor
    for the class of quaternary matroids.

    EXAMPLES::

        sage: M = matroids.catalog.P6(); M
        P6: Matroid of rank 3 on 6 elements with circuit-closures
        {2: {{'a', 'b', 'c'}}, 3: {{'a', 'b', 'c', 'd', 'e', 'f'}}}
        sage: len(set(M.nonspanning_circuits()).difference(M.nonbases())) == 0
        True
        sage: Matroid(matrix=random_matrix(GF(4, 'a'), ncols=5, nrows=5)).has_minor(M)
        False
        sage: M.is_valid()
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 641-2.
    """
    CC = {2: ['abc'], 3: ['abcdef']}
    M = Matroid(circuit_closures=CC)
    M = _rename_and_relabel(M, "P6", groundset)
    return M


def U36(groundset='abcdef'):
    """
    The uniform matroid of rank `3` on `6` elements.

    Six points freely placed in the plane; the tipless free `3`-spike.
    Identically self-dual.

    EXAMPLES::

        sage: U36 = matroids.catalog.U36(); U36
        U(3, 6): Matroid of rank 3 on 6 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'd', 'e', 'f'}}}
        sage: Z = matroids.Spike(3, False)
        sage: U36.is_isomorphic(Z)
        True
        sage: U36.equals(U36.dual())
        True
        sage: U36.automorphism_group().structure_description()
        'S6'

    REFERENCES:

    [Oxl2011]_, p. 642.

    TESTS::

        sage: M = matroids.catalog.U36(range(6))
        sage: sorted(M.groundset())
        [0, 1, 2, 3, 4, 5]
        sage: M.is_isomorphic(matroids.catalog.U36())
        True
    """
    M = Uniform(3, 6)
    M = _rename_and_relabel(M, "U(3, 6)", groundset)
    return M


def R6(groundset='abcdef'):
    """
    Return the matroid `R_6`, represented over `GF(3)`.

    The matroid `R_6` is a `6`-element matroid of rank-`3`. It is
    representable over a field if and only if that field has at least three
    elements. It is isomorphic to the `2`-sum of two copies of `U_{2, 4}`.

    EXAMPLES::

        sage: M = matroids.catalog.R6(); M
        R6: Ternary matroid of rank 3 on 6 elements, type 2+
        sage: M.equals(M.dual())
        True
        sage: M.is_connected()
        True
        sage: M.is_3connected()
        False
        sage: M.automorphism_group().is_transitive()
        True

    REFERENCES:

    [Oxl2011]_, p. 642.
    """
    A = Matrix(
        GF(3), [[1, 0, 0, 1, 1, 1], [0, 1, 0, 1, 2, 1], [0, 0, 1, 1, 0, 2]]
    )
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "R6")
    return M


def Fano(groundset='abcdefg'):
    r"""
    Return the Fano matroid, represented over `GF(2)`.

    The Fano matroid, or Fano plane, or `F_7`, is a `7`-element matroid of
    rank-`3`. It is representable over a field if and only if that field has
    characteristic two.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.Fano(); M
        Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
        sage: M.automorphism_group().is_transitive()
        True
        sage: M.automorphism_group().structure_description()
        'PSL(3,2)'

    Every single-element deletion of `F_7` is isomorphic to `M(K_4)`::

        sage: setprint(sorted(M.nonspanning_circuits()))
        [{'a', 'b', 'f'}, {'a', 'c', 'e'}, {'a', 'd', 'g'}, {'b', 'c', 'd'},
         {'b', 'e', 'g'}, {'c', 'f', 'g'}, {'d', 'e', 'f'}]
        sage: M.delete(M.groundset_list()[randrange(0,
        ....:                  7)]).is_isomorphic(matroids.CompleteGraphic(4))
        True

    It is also the projective plane of order two, i.e. `PG(2, 2)`::

        sage: M.is_isomorphic(matroids.PG(2, 2))
        True

    `F_7` is isomorphic to the unique binary `3`-spike::

        sage: M.is_isomorphic(matroids.Z(3))
        True

    REFERENCES:

    [Oxl2011]_, p. 643.
    """
    A = Matrix(
        GF(2),
        [[1, 0, 0, 0, 1, 1, 1], [0, 1, 0, 1, 0, 1, 1], [0, 0, 1, 1, 1, 0, 1]]
    )
    M = BinaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "Fano")
    return M


def FanoDual(groundset='abcdefg'):
    """
    Return the dual of the Fano matroid.

    `F_7^*` is a `7`-element matroid of rank-`3`.

    EXAMPLES::

        sage: F7 = matroids.catalog.Fano()
        sage: F7D = matroids.catalog.FanoDual(); F7D
        F7*: Binary matroid of rank 4 on 7 elements, type (3, 7)
        sage: F7.is_isomorphic(F7D.dual())
        True
        sage: F7D.automorphism_group().is_transitive()
        True
        sage: F7D.automorphism_group().structure_description()
        'PSL(3,2)'

    Every single-element deletion of `F_7^*` is isomorphic to `M(K_{2, 3})`::

        sage: K2_3 = Matroid(graphs.CompleteBipartiteGraph(2, 3))
        sage: import random
        sage: e = random.choice(list(F7D.groundset()))
        sage: F7D.delete(e).is_isomorphic(K2_3)
        True

    REFERENCES:

    [Oxl2011]_, p. 643.

    TESTS::

        sage: F7D = matroids.catalog.FanoDual(range(7))
        sage: sorted(F7D.groundset())
        [0, 1, 2, 3, 4, 5, 6]
    """
    M = Fano().dual()
    M = _rename_and_relabel(M, "F7*", groundset)
    return M


def NonFano(groundset='abcdefg'):
    """
    Return the non-Fano matroid, represented over `GF(3)`.

    The non-Fano matroid, or `F_7^-`, is a `7`-element matroid of rank-`3`. It
    is representable over a field if and only if that field has characteristic
    other than two. It is the unique relaxation of `F_7`.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.NonFano(); M
        NonFano: Ternary matroid of rank 3 on 7 elements, type 0-
        sage: setprint(M.nonbases())
        [{'a', 'b', 'f'}, {'a', 'c', 'e'}, {'a', 'd', 'g'},
         {'b', 'c', 'd'}, {'b', 'e', 'g'}, {'c', 'f', 'g'}]
        sage: M.delete('f').is_isomorphic(matroids.CompleteGraphic(4))
        True
        sage: M.delete('g').is_isomorphic(matroids.CompleteGraphic(4))
        False

    REFERENCES:

    [Oxl2011]_, p. 643-4.
    """
    A = Matrix(
        GF(3),
        [[1, 0, 0, 0, 1, 1, 1], [0, 1, 0, 1, 0, 1, 1], [0, 0, 1, 1, 1, 0, 1]]
    )
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "NonFano")
    return M


def NonFanoDual(groundset='abcdefg'):
    r"""
    Return the dual of the non-Fano matroid.

    `(F_7^-)^*` is a `7`-element matroid of rank-`3`. Every single-element
    contraction of `(F_7 )^*` is isomorphic to `M(K_4)` or `\mathcal{W}^3`.

    EXAMPLES::

        sage: M = matroids.catalog.NonFanoDual(); M
        NonFano*: Ternary matroid of rank 4 on 7 elements, type 0-
        sage: sorted(M.groundset())
        ['a', 'b', 'c', 'd', 'e', 'f', 'g']

    Every single-element contraction of `(F_7^-)^*` is isomorphic to `M(K_4)` or
    `\mathcal{W}^3`::

        sage: import random
        sage: e = random.choice(list(M.groundset()))
        sage: N = M.contract(e)
        sage: K4 = matroids.catalog.K4()
        sage: W3 = matroids.catalog.Whirl3()
        sage: N.is_isomorphic(K4) or N.is_isomorphic(W3)
        True

    REFERENCES:

    [Oxl2011]_, p. 643-4.

    TESTS::

        sage: M = matroids.catalog.NonFanoDual(range(7))
        sage: sorted(M.groundset())
        [0, 1, 2, 3, 4, 5, 6]
        sage: M.is_isomorphic(matroids.catalog.NonFanoDual())
        True
    """
    M = NonFano().dual()
    M = _rename_and_relabel(M, "NonFano*", groundset)
    return M


def O7(groundset='abcdefg'):
    """
    Return the matroid `O_7`, represented over `GF(3)`.

    The matroid `O_7` is a `7`-element matroid of rank-`3`. It is representable
    over a field if and only if that field has at least three elements. It is
    obtained by freely adding a point to any line of `M(K_4)`.

    EXAMPLES::

        sage: M = matroids.catalog.O7(); M
        O7: Ternary matroid of rank 3 on 7 elements, type 0+
        sage: M.delete('e').is_isomorphic(matroids.CompleteGraphic(4))
        True
        sage: M.tutte_polynomial()
        y^4 + x^3 + x*y^2 + 3*y^3 + 4*x^2 + 5*x*y + 5*y^2 + 4*x + 4*y

    REFERENCES:

    [Oxl2011]_, p. 644.
    """
    A = Matrix(
        GF(3),
        [[1, 0, 0, 1, 1, 1, 1], [0, 1, 0, 0, 1, 2, 2], [0, 0, 1, 1, 0, 1, 0]]
    )
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "O7")
    return M


def P7(groundset='abcdefg'):
    """
    Return the matroid `P_7`, represented over `GF(3)`.

    The matroid `P_7` is a `7`-element matroid of rank-`3`. It is
    representable over a field if and only if that field has at least three
    elements. It is one of two ternary `3`-spikes, with the other being
    `F_7^-`.

    EXAMPLES::

        sage: M = matroids.catalog.P7(); M
        P7: Ternary matroid of rank 3 on 7 elements, type 1+
        sage: M.whitney_numbers2()
        [1, 7, 11, 1]
        sage: M.has_minor(matroids.CompleteGraphic(4))
        False
        sage: M.is_valid()
        True

    REFERENCES:

    [Oxl2011]_, p. 644-5.
    """
    A = Matrix(
        GF(3),
        [[1, 0, 0, 2, 1, 1, 0], [0, 1, 0, 1, 1, 0, 1], [0, 0, 1, 1, 0, 1, 1]]
    )
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "P7")
    return M


def AG32(groundset='abcdefgh'):
    """
    Return the matroid `AG(3, 2)`.

    The binary affine cube.

    EXAMPLES::

        sage: M = matroids.catalog.AG32(); M
        AG(3, 2): Binary matroid of rank 4 on 8 elements, type (4, 0)
        sage: M.is_valid() and M.is_3connected()
        True

    `AG(3, 2)` is isomorphic to the unique tipless binary `4`-spike::

        sage: M.is_isomorphic(matroids.Z(4, False))
        True

    and it is identically self-dual::

        sage: M.equals(M.dual())
        True

    Every single-element deletion is isomorphic to `F_7^*` and every single-element
    contraction is isomorphic to `F_7`::

        sage: F7 = matroids.catalog.Fano()
        sage: F7D = matroids.catalog.FanoDual()
        sage: import random
        sage: e = random.choice(list(M.groundset()))
        sage: M.delete(e).is_isomorphic(F7D)
        True
        sage: M.contract(e).is_isomorphic(F7)
        True

    REFERENCES:

    [Oxl2011]_, p. 645.
    """
    M = AG(3, 2)
    M = _rename_and_relabel(M, "AG(3, 2)", groundset)
    return M


def AG32prime(groundset=None):
    """
    Return the matroid `AG(3, 2)'`, represented as circuit closures.

    The matroid `AG(3, 2)'` is an `8`-element matroid of rank-`4`. It is a
    smallest non-representable matroid. It is the unique relaxation of
    `AG(3, 2)`.

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
        sage: setprint(M.noncospanning_cocircuits())
        [{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
         {'a', 'c', 'd', 'f'}, {'a', 'd', 'g', 'h'}, {'a', 'e', 'f', 'h'},
         {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'}, {'b', 'd', 'f', 'h'},
         {'b', 'e', 'g', 'h'}, {'c', 'd', 'e', 'h'}, {'c', 'f', 'g', 'h'},
         {'d', 'e', 'f', 'g'}]
        sage: M.is_valid()
        True
        sage: M.automorphism_group().is_transitive()
        False

    Self-dual but not identically self-dual::

        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True

    Every single-element deletion is isomorphic to `F_7^*` or `(F_7^-)^*` and every
    single-element contraction is isomorphic to `F_7` or `F_7^-`::

        sage: F7 = matroids.catalog.Fano()
        sage: F7D = matroids.catalog.FanoDual()
        sage: F7m = matroids.catalog.NonFano()
        sage: F7mD = matroids.catalog.NonFanoDual()
        sage: import random
        sage: e = random.choice(list(M.groundset()))
        sage: M.delete(e).is_isomorphic(F7D) or M.delete(e).is_isomorphic(F7mD)
        True
        sage: Me = M.contract(e)
        sage: Me.is_isomorphic(F7) or Me.is_isomorphic(F7m)
        True

    REFERENCES:

    [Oxl2011]_, p. 646.
    """
    CC = {
        3: [
            'abfg', 'bcdg', 'defg', 'cdeh', 'aefh', 'abch', 'abed',
            'cfgh', 'bcef', 'adgh', 'acdf', 'begh', 'aceg',
        ],
        4: ['abcdefgh'],
    }
    M = Matroid(circuit_closures=CC)
    M = _rename_and_relabel(M, "AG(3, 2)'", groundset)
    return M


def R8(groundset='abcdefgh'):
    """
    Return the matroid `R_8`, represented over `GF(3)`.

    The matroid `R_8` is an `8`-element matroid of rank-`4`. It is
    representable over a field if and only if the characteristic of that field
    is not two. It is the real affine cube.

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
        sage: M.automorphism_group().is_transitive()
        True

    Every single-element deletion is isomorphic to (F_7^-)^* and every
    single-element contraction is isomorphic to F_7^-::

        sage: F7m = matroids.catalog.NonFano()
        sage: F7mD = matroids.catalog.NonFanoDual()
        sage: import random
        sage: e = random.choice(list(M.groundset()))
        sage: M.delete(e).is_isomorphic(F7mD)
        True
        sage: M.contract(e).is_isomorphic(F7m)
        True

    REFERENCES:

    [Oxl2011]_, p. 646.
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
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "R8")
    return M


def F8(groundset=None):
    """
    Return the matroid `F_8`, represented as circuit closures.

    The matroid `F_8` is an `8`-element matroid of rank-`4`. It is a smallest
    non-representable matroid.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.catalog.F8(); M
        F8: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
             {'a', 'c', 'd', 'f'}, {'a', 'c', 'e', 'g'}, {'a', 'd', 'g', 'h'},
             {'a', 'e', 'f', 'h'}, {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'},
             {'c', 'd', 'e', 'h'}, {'c', 'f', 'g', 'h'}, {'d', 'e', 'f', 'g'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: D = get_nonisomorphic_matroids([M.contract(i) for i in M.groundset()])
        sage: len(D)
        3
        sage: [N.is_isomorphic(matroids.catalog.Fano()) for N in D]
        [...True...]
        sage: [N.is_isomorphic(matroids.catalog.NonFano()) for N in D]
        [...True...]
        sage: M.is_valid() and M.is_paving()
        True
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 647.
    """
    CC = {
        3: [
            'abfg', 'bcdg', 'defg', 'cdeh',
            'aefh', 'abch', 'abed', 'cfgh',
            'bcef', 'adgh', 'acdf', 'aceg',
        ],
        4: ['abcdefgh'],
    }
    M = Matroid(circuit_closures=CC)
    M = _rename_and_relabel(M, "F8", groundset)
    return M


def Q8(groundset=None):
    """
    Return the matroid `Q_8`, represented as circuit closures.

    The matroid `Q_8` is an `8`-element matroid of rank-`4`. It is a smallest
    non-representable matroid.

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
        sage: M.is_valid()  # long time
        True
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 647.
    """
    CC = {
        3: [
            'abfg', 'bcdg', 'defg', 'cdeh', 'aefh', 'abch',
            'abed', 'cfgh', 'bcef', 'adgh', 'acdf',
        ],
        4: ['abcdefgh'],
    }
    M = Matroid(circuit_closures=CC)
    M = _rename_and_relabel(M, "Q8", groundset)
    return M


def L8(groundset=None):
    """
    Return the matroid `L_8`, represented as circuit closures.

    The matroid `L_8` is an `8`-element matroid of rank-`4`. It is
    representable over all fields with at least five elements. It is a cube,
    yet it is not a tipless spike.

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
        sage: M.is_valid()  # long time
        True
        sage: M.automorphism_group().is_transitive()
        True

    Every single-element contraction is isomorphic to the free extension of
    `M(K_4)`::

        sage: K4 = matroids.catalog.K4(range(6))
        sage: Bext = [list(b) for b in K4.bases()] + [list(I)+[6] for I in
        ....:                                         K4.independent_sets(2)]
        sage: K4ext = Matroid(bases=Bext)
        sage: import random
        sage: e = random.choice(list(M.groundset()))
        sage: M.contract(e).is_isomorphic(K4ext)
        True

    REFERENCES:

    [Oxl2011]_, p. 648.
    """
    CC = {3: ['abfg', 'bcdg', 'defg', 'cdeh', 'aefh', 'abch', 'aceg', 'bdfh'],
          4: ['abcdefgh']}
    M = Matroid(circuit_closures=CC)
    M = _rename_and_relabel(M, "L8", groundset)
    return M


def S8(groundset='abcdefgh'):
    """
    Return the matroid `S_8`, represented over `GF(2)`.

    The matroid `S_8` is an `8`-element matroid of rank-`4`. It is
    representable over a field if and only if that field has characteristic
    two. It is the unique deletion of a non-tip element from the binary
    `4`-spike.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.catalog.S8(); M
        S8: Binary matroid of rank 4 on 8 elements, type (2, 0)
        sage: M.contract('d').is_isomorphic(matroids.catalog.Fano())
        True
        sage: M.delete('d').is_isomorphic(matroids.catalog.FanoDual())
        False
        sage: M.delete('h').is_isomorphic(matroids.catalog.FanoDual())
        True
        sage: M.is_graphic()
        False
        sage: D = get_nonisomorphic_matroids(
        ....:       list(matroids.catalog.Fano().linear_coextensions(cosimple=True)))
        sage: len(D)
        2
        sage: [N.is_isomorphic(M) for N in D]
        [...True...]
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 648.
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
    M = BinaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "S8")
    return M


def Vamos(groundset=None):
    r"""
    Return the `V\acute{a}mos` matroid, represented as circuit closures.

    The `V\acute{a}mos` matroid, or `V\acute{a}mos` cube, or `V_8` is an
    `8`-element matroid of rank-`4`. It violates Ingleton's condition for
    representability over a division ring. It is not algebraic.

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
        sage: M.is_valid() and M.is_paving()  # long time
        True
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 649.
    """
    CC = {3: ['abcd', 'abef', 'cdef', 'abgh', 'efgh'], 4: ['abcdefgh']}
    M = Matroid(circuit_closures=CC)
    M = _rename_and_relabel(M, "Vamos", groundset)
    return M


def T8(groundset='abcdefgh'):
    """
    Return the matroid `T_8`, represented over `GF(3)`.

    The matroid `T_8` is an `8`-element matroid of rank-`4`. It is
    representable over a field if and only if that field has characteristic
    three. It is an excluded minor for the dyadic matroids.

    EXAMPLES::

        sage: M = matroids.catalog.T8(); M
        T8: Ternary matroid of rank 4 on 8 elements, type 0-
        sage: M.truncation().is_isomorphic(matroids.Uniform(3, 8))
        True
        sage: M.contract('e').is_isomorphic(matroids.catalog.P7())
        True
        sage: M.has_minor(matroids.Uniform(3, 8))
        False
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 649.
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
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "T8")
    return M


def J(groundset='abcdefgh'):
    """
    Return the matroid `J`, represented over `GF(3)`.

    The matroid `J` is an `8`-element matroid of rank-`4`. It is representable
    over a field if and only if that field has at least three elements.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.J(); M
        J: Ternary matroid of rank 4 on 8 elements, type 0-
        sage: setprint(M.truncation().nonbases())
        [{'a', 'b', 'f'}, {'a', 'c', 'g'}, {'a', 'd', 'h'}]
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.has_minor(matroids.CompleteGraphic(4))
        False
        sage: M.is_valid()
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 650.
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
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "J")
    return M


def P8(groundset='abcdefgh'):
    """
    Return the matroid `P_8`, represented over `GF(3)`.

    The matroid `P_8` is an `8`-element matroid of rank-`4`. It is uniquely
    representable over all fields of characteristic other than two. It is an
    excluded minor for all fields of characteristic two with four or more
    elements.

    EXAMPLES::

        sage: M = matroids.catalog.P8(); M
        P8: Ternary matroid of rank 4 on 8 elements, type 2+
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: Matroid(matrix=random_matrix(GF(4, 'a'), ncols=5, nrows=5)).has_minor(M)
        False
        sage: M.bicycle_dimension()
        2

    REFERENCES:

    [Oxl2011]_, p. 650-1.
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
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "P8")
    return M


def P8pp(groundset=None):
    """
    Return the matroid `P_8^=`, represented as circuit closures.

    The matroid `P_8^=` is an `8`-element matroid of rank-`4`. It can be
    obtained from `P_8` by relaxing the unique pair of disjoint
    circuit-hyperplanes. It is an excluded minor for `GF(4)`-representability.
    It is representable over all fields with at least five elements.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.catalog.P8pp(); M
        P8'': Matroid of rank 4 on 8 elements with 8 nonspanning circuits
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: len(get_nonisomorphic_matroids([M.contract(i) for i in M.groundset()]))
        1
        sage: M.is_valid() and M.is_paving()
        True

    REFERENCES:

    [Oxl2011]_, p. 651.
    """
    NSC = ['abfh', 'bceg', 'cdfh', 'adeg', 'acef', 'bdfg', 'acgh', 'bdeh']
    M = Matroid(rank=4, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, "P8''", groundset)
    return M


def Wheel4(groundset='abcdefgh'):
    """
    Return the rank-`4` wheel.

    A regular, graphic, and cographic matroid. Self-dual but not identically
    self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.Wheel4(); M
        Wheel(4): Regular matroid of rank 4 on 8 elements with 45 bases
        sage: M.is_valid() and M.is_graphic() and M.dual().is_graphic()
        True
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 651-2.
    """
    M = Wheel(4)
    M = _rename_and_relabel(M, "Wheel(4)", groundset)
    return M


def Whirl4(groundset='abcdefgh'):
    """
    Return the rank-`4` whirl.

    A matroid which is not graphic, not cographic, and not regular. Self-dual
    but not identically self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.Whirl4(); M
        Whirl(4): Ternary matroid of rank 4 on 8 elements, type 0+
        sage: M.is_valid()
        True
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 652.
    """
    M = Whirl(4)
    M = _rename_and_relabel(M, "Whirl(4)", groundset)
    return M


def K33dual(groundset='abcdefghi'):
    """
    Return the matroid `M*(K_{3, 3})`, represented over the regular partial
    field.

    The matroid `M*(K_{3, 3})` is a `9`-element matroid of rank-`4`. It is an
    excluded minor for the class of graphic matroids. It is the graft matroid
    of the `4`-wheel with every vertex except the hub being coloured.

    EXAMPLES::

        sage: # needs sage.graphs
        sage: M = matroids.catalog.K33dual(); M
        M*(K3, 3): Regular matroid of rank 4 on 9 elements with 81 bases
        sage: any(N.is_3connected() for N in M.linear_extensions(simple=True))
        False
        sage: M.is_valid()
        True
        sage: M.automorphism_group().is_transitive()
        True

    REFERENCES:

    [Oxl2011]_, p. 652-3.
    """
    from sage.graphs.graph_generators import graphs

    G = graphs.CompleteBipartiteGraph(3, 3)
    M = Matroid(groundset=groundset, graph=G, regular=True)
    M = M.dual()
    M = _rename_and_relabel(M, "M*(K3, 3)")
    return M


def K33(groundset='abcdefghi'):
    r"""
    Return the graphic matroid `M(K_{3,3})`.

    `M(K_{3,3})` is an excluded minor for the class of cographic matroids.

    EXAMPLES::

        sage: # needs sage.graphs
        sage: M = matroids.catalog.K33(); M
        M(K3, 3): Regular matroid of rank 5 on 9 elements with 81 bases
        sage: M.is_valid()
        True
        sage: G1 = M.automorphism_group()
        sage: G2 = matroids.catalog.K33dual().automorphism_group()
        sage: G1.is_isomorphic(G2)
        True

    REFERENCES:

    [Oxl2011]_, p. 652-3.
    """
    from sage.graphs.graph_generators import graphs

    G = graphs.CompleteBipartiteGraph(3, 3)
    M = Matroid(groundset=groundset, graph=G, regular=True)
    M = _rename_and_relabel(M, "M(K3, 3)")
    return M


def AG23(groundset='abcdefghi'):
    """
    Return the matroid `AG(2, 3)`.

    The ternary affine plane. A matroid which is not graphic, not cographic,
    and not regular.

    EXAMPLES::

        sage: M = matroids.catalog.AG23(); M
        AG(2, 3): Ternary matroid of rank 3 on 9 elements, type 3+
        sage: M.is_valid() and M.is_3connected() and M.is_ternary()
        True
        sage: M.has_minor(matroids.catalog.K4())
        False
        sage: import random
        sage: e = random.choice(list(M.groundset()))
        sage: M.delete(e).is_isomorphic(matroids.catalog.AG23minus())
        True
        sage: M.automorphism_group().is_transitive()
        True

    REFERENCES:

    [Oxl2011]_, p. 653.
    """
    M = AG(2, 3)
    M = _rename_and_relabel(M, "AG(2, 3)", groundset)
    return M


def TernaryDowling3(groundset='abcdefghi'):
    r"""
    Return the matroid `Q_3(GF(3)^\times)`, represented over `GF(3)`.

    The matroid `Q_3(GF(3)^\times)` is a `9`-element matroid of rank-`3`. It
    is the rank-`3` ternary Dowling geometry. It is representable over a field
    if and only if that field does not have characteristic two.

    EXAMPLES::

        sage: M = matroids.catalog.TernaryDowling3(); M
        Q3(GF(3)x): Ternary matroid of rank 3 on 9 elements, type 0-
        sage: len(list(M.linear_subclasses()))
        72
        sage: M.fundamental_cycle('abc', 'd')
        {'a': 2, 'b': 1, 'd': 1}
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 654.
    """
    A = Matrix(
        GF(3),
        [
            [1, 0, 0, 1, 1, 0, 0, 1, 1],
            [0, 1, 0, 2, 1, 1, 1, 0, 0],
            [0, 0, 1, 0, 0, 2, 1, 2, 1],
        ],
    )
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "Q3(GF(3)x)")
    return M


def R9(groundset=None):
    """
    Return the matroid `R_9`.

    The ternary Reid geometry. The only `9`-element rank-`3` simple ternary
    matroids are `R_9`, `Q_3(GF(3)^\times)`, and `AG(2, 3)`. It is not graphic,
    not cographic, and not regular.

    EXAMPLES::

        sage: M = matroids.catalog.R9(); M
        R9: Matroid of rank 3 on 9 elements with 15 nonspanning circuits
        sage: M.is_valid()
        True
        sage: len(M.nonspanning_circuits())
        15
        sage: M.is_simple() and M.is_ternary()
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 654.
    """
    NSC = ['abc', 'abd', 'acd', 'aef', 'agh', 'bcd', 'bfh', 'bgi',
           'ceg', 'cfi', 'deh', 'dei', 'dfg', 'dhi', 'ehi']
    M = Matroid(rank=3, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, "R9", groundset)
    return M


def Pappus(groundset=None):
    """
    Return the Pappus matroid.

    The Pappus matroid is a `9`-element matroid of rank-`3`. It is
    representable over a field if and only if that field either has 4 elements
    or more than 7 elements. It is an excluded minor for the class of
    GF(5)-representable matroids.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.Pappus(); M
        Pappus: Matroid of rank 3 on 9 elements with 9 nonspanning circuits
        sage: setprint(M.nonspanning_circuits())
        [{'a', 'b', 'c'}, {'a', 'e', 'i'}, {'a', 'f', 'h'}, {'b', 'd', 'i'},
         {'b', 'f', 'g'}, {'c', 'd', 'h'}, {'c', 'e', 'g'}, {'d', 'e', 'f'},
         {'g', 'h', 'i'}]
        sage: M.is_dependent(['d', 'e', 'f'])
        True
        sage: M.is_valid()
        True
        sage: M.automorphism_group().is_transitive()
        True

    REFERENCES:

    [Oxl2011]_, p. 655.
    """
    NSC = ['abc', 'def', 'ceg', 'bfg', 'cdh', 'afh', 'bdi', 'aei', 'ghi']
    M = Matroid(rank=3, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, "Pappus", groundset)
    return M


def NonPappus(groundset=None):
    """
    Return the non-Pappus matroid.

    The non-Pappus matroid is a `9`-element matroid of rank-`3`. It is not
    representable over any commutative field. It is the unique relaxation of
    the Pappus matroid.

    EXAMPLES::

        sage: M = matroids.catalog.NonPappus(); M
        NonPappus: Matroid of rank 3 on 9 elements with 8 nonspanning circuits
        sage: NSC = set([('a', 'b', 'c'), ('a', 'e', 'i'), ('a', 'f', 'h'),
        ....:            ('b', 'd', 'i'), ('b', 'f', 'g'), ('c', 'd', 'h'),
        ....:            ('c', 'e', 'g'), ('g', 'h', 'i')])
        sage: NSC == set(tuple(sorted(C)) for C in M.nonspanning_circuits())
        True
        sage: M.is_dependent(['d', 'e', 'f'])
        False
        sage: M.is_valid() and M.is_paving()
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 655.
    """
    NSC = ['abc', 'ceg', 'bfg', 'cdh', 'afh', 'bdi', 'aei', 'ghi']
    M = Matroid(rank=3, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, "NonPappus", groundset)
    return M


def K5(groundset='abcdefghij'):
    """
    Return the graphic matroid `M(K_5)`.

    `M(K_5)` is an excluded minor for the class of cographic matroids. It is
    the `3`-dimensional Desargues conﬁguration.

    EXAMPLES::

        sage: M = matroids.catalog.K5(); M
        M(K5): Graphic matroid of rank 4 on 10 elements
        sage: M.is_valid()
        True
        sage: M.automorphism_group().is_transitive()
        True

    REFERENCES:

    [Oxl2011]_, p. 656.
    """
    M = CompleteGraphic(5)
    M = _rename_and_relabel(M, "M(K5)", groundset)
    return M


def K5dual(groundset='abcdefghij'):
    """
    Return the matroid `M^*(K_5)`.

    `M^*(K_5)` is an excluded minor for the class of graphic matroids.

    EXAMPLES::

        sage: M = matroids.catalog.K5dual(); M
        M*(K5): Dual of 'Graphic matroid of rank 4 on 10 elements'
        sage: M.is_3connected()
        True
        sage: G1 = M.automorphism_group()
        sage: G2 = matroids.catalog.K5().automorphism_group()
        sage: G1.is_isomorphic(G2)
        True

    REFERENCES:

    [Oxl2011]_, p. 656.
    """
    M = CompleteGraphic(5).dual()
    M = _rename_and_relabel(M, "M*(K5)", groundset)
    return M


def R10(groundset='abcdefghij'):
    """
    Return the matroid `R_{10}`, represented over the regular partial field.

    The NonDesargues matroid is a `10`-element matroid of rank-5. It is the
    unique splitter for the class of regular matroids. It is the graft matroid
    of `K_{3, 3}` in which every vertex is coloured.

    EXAMPLES::

        sage: M = matroids.catalog.R10(); M
        R10: Regular matroid of rank 5 on 10 elements with 162 bases
        sage: cct = []
        sage: for i in M.circuits():
        ....:     cct.append(len(i))
        sage: Set(cct)
        {4, 6}
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.is_valid()
        True
        sage: M.automorphism_group().is_transitive()
        True

    Every single-element deletion is isomorphic to `M(K_{3, 3})`, and every
    single-element contraction is isomorphic to `M^*(K_{3, 3})`::

        sage: K33 = matroids.catalog.K33()
        sage: K33D = matroids.catalog.K33dual()
        sage: import random
        sage: e = random.choice(list(M.groundset()))
        sage: M.delete(e).is_isomorphic(K33)
        True
        sage: M.contract(e).is_isomorphic(K33D)
        True

    Check the splitter property::

        sage: matroids.catalog.R10().linear_extensions(simple=True)
        []

    REFERENCES:

    [Oxl2011]_, p. 656-7.
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
    M = RegularMatroid(A, groundset)
    M = _rename_and_relabel(M, "R10")
    return M


def NonDesargues(groundset=None):
    """
    Return the NonDesargues matroid.

    The NonDesargues matroid is a `10`-element matroid of rank-`3`. It is not
    representable over any division ring. It is not graphic, not cographic, and
    not regular.

    EXAMPLES::

        sage: M = matroids.catalog.NonDesargues(); M
        NonDesargues: Matroid of rank 3 on 10 elements with 9 nonspanning circuits
        sage: M.is_valid()
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 657.
    """
    NSC = ['acj', 'aef', 'bce', 'bfj', 'bgi', 'chi', 'dfg', 'dij', 'egh']
    M = Matroid(rank=3, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, "NonDesargues", groundset)
    return M


def R12(groundset='abcdefghijkl'):
    """
    Return the matroid `R_{12}`, represented over the regular partial field.

    The matroid `R_{12}` is a `12`-element regular matroid of rank-`6`. It
    induces a `3`-separation in its `3`-connected majors within the class of
    regular matroids. An excluded minor for the class of graphic or cographic
    matroids.

    EXAMPLES::

        sage: M = matroids.catalog.R12(); M
        R12: Regular matroid of rank 6 on 12 elements with 441 bases
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.is_valid()
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 657.
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
    M = RegularMatroid(A, groundset)
    M = _rename_and_relabel(M, "R12")
    return M


def ExtendedTernaryGolayCode(groundset='abcdefghijkl'):
    """
    Return the matroid of the extended ternary Golay code.

    This is the unique Steiner system `S(5, 6, 12)`.

    EXAMPLES::

        sage: M = matroids.catalog.ExtendedTernaryGolayCode(); M
        Extended Ternary Golay Code: Ternary matroid of rank 6 on 12 elements,
        type 6+
        sage: C = LinearCode(M.representation())
        sage: C.is_permutation_equivalent(codes.GolayCode(GF(3)))
        True
        sage: M.is_valid()
        True

    The automorphism group is the `5`-transitive Mathieu group `M12`:

        sage: # long time
        sage: G = M.automorphism_group()
        sage: G.is_transitive()
        True
        sage: G.structure_description()
        'M12'

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

    .. SEEALSO::

        :class:`GolayCode <sage.coding.golay_code.GolayCode>`

    REFERENCES:

    [Oxl2011]_, p. 658.
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 0],
        [0, 1, 0, 0, 0, 0, 1, 1, 2, 1, 0, 2],
        [0, 0, 1, 0, 0, 0, 1, 2, 1, 0, 1, 2],
        [0, 0, 0, 1, 0, 0, 1, 2, 0, 1, 2, 1],
        [0, 0, 0, 0, 1, 0, 1, 0, 2, 2, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1]
    ])
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "Extended Ternary Golay Code")
    return M


def T12(groundset='abcdefghijkl'):
    """
    Return the matroid `T_{12}`.

    The edges of the Petersen graph can be labeled by the `4`-circuits of
    `T_{12}` so that two edges are adjacent if and only if the corresponding
    `4`-circuits overlap in exactly two elements. Relaxing a
    circuit-hyperplane yields an excluded minor for the class of matroids that
    are either binary or ternary.

    EXAMPLES::

        sage: M = matroids.catalog.T12(); M
        T12: Binary matroid of rank 6 on 12 elements, type (2, None)
        sage: M.is_valid()
        True
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.automorphism_group().is_transitive()
        True

    REFERENCES:

    [Oxl2011]_, p. 658-9.
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
    M = BinaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "T12")
    return M


def PG23(groundset=None):
    """
    Return the matroid `PG23`.

    The second smallest projective plane. Not graphic, not cographic, not
    regular, not near-regular.

    EXAMPLES::

        sage: M = matroids.catalog.PG23(); M
        PG(2, 3): Ternary matroid of rank 3 on 13 elements, type 3+
        sage: M.is_3connected()
        True
        sage: M.automorphism_group().is_transitive()
        True

    REFERENCES:

    [Oxl2011]_, p. 659.
    """
    M = PG(2, 3)
    M = _rename_and_relabel(M, groundset=groundset)
    return M


def Wheel(r, field=None, ring=None, groundset=None):
    r"""
    Return the rank-`r` wheel.

    INPUT:

    - ``r`` -- positive integer; the rank of the matroid
    - ``ring`` -- any ring; if provided, output will be a linear matroid
      over the ring or field ``ring``. If the ring is `\ZZ`, then output
      will be a regular matroid.
    - ``field`` -- any field; same as ``ring``, but only fields are allowed
    - ``groundset`` -- string (optional); the groundset of the matroid

    OUTPUT: the rank-`r` wheel matroid, represented as a regular matroid

    EXAMPLES::

        sage: M = matroids.Wheel(5); M
        Wheel(5): Regular matroid of rank 5 on 10 elements with 121 bases
        sage: M.tutte_polynomial()
        x^5 + y^5 + 5*x^4 + 5*x^3*y + 5*x^2*y^2 + 5*x*y^3 + 5*y^4 + 10*x^3 +
        15*x^2*y + 15*x*y^2 + 10*y^3 + 10*x^2 + 16*x*y + 10*y^2 + 4*x + 4*y
        sage: M.is_valid()
        True
        sage: M = matroids.Wheel(3)
        sage: M.is_isomorphic(matroids.CompleteGraphic(4))
        True
        sage: M.is_isomorphic(matroids.Wheel(3, field=GF(3)))
        True
        sage: M = matroids.Wheel(3, field=GF(3)); M
        Wheel(3): Ternary matroid of rank 3 on 6 elements, type 0+

    For `r \ge 2`, the wheel is self-dual but not identically self-dual, and
    for `r \ge 4` it has a non-transitive automorphism group::

        sage: import random
        sage: r = random.choice(range(4, 8))
        sage: M = matroids.Wheel(r)
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 659-60.
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
    A = Matrix(base_ring, r, 2 * r, sparse=True)
    for i in range(r):
        A[i, i] = 1
        A[i, r + i] = 1
        if i != 0:
            A[i, i + r - 1] = -1
        else:
            A[i, 2 * r - 1] = -1
    if base_ring is ZZ:
        M = RegularMatroid(A)
    else:
        M = Matroid(A)
    M = _rename_and_relabel(M, f'Wheel({r})', groundset)
    return M


def Whirl(r, groundset=None):
    r"""
    Return the rank-`r` whirl.

    The whirl is the unique relaxation of the wheel.

    INPUT:

    - ``r`` -- positive integer; the rank of the matroid
    - ``groundset`` -- string (optional); the groundset of the matroid

    OUTPUT: the rank-`r` whirl matroid, represented as a ternary matroid

    EXAMPLES::

        sage: M = matroids.Whirl(5); M
        Whirl(5): Ternary matroid of rank 5 on 10 elements, type 0-
        sage: M.is_valid()
        True
        sage: M.tutte_polynomial()
        x^5 + y^5 + 5*x^4 + 5*x^3*y + 5*x^2*y^2 + 5*x*y^3 + 5*y^4 + 10*x^3 + 15*x^2*y +
         15*x*y^2 + 10*y^3 + 10*x^2 + 15*x*y + 10*y^2 + 5*x + 5*y
        sage: M.is_isomorphic(matroids.Wheel(5))
        False
        sage: M = matroids.Whirl(3)
        sage: M.is_isomorphic(matroids.CompleteGraphic(4))
        False

    For `r \ge 3`, the whirl is self-dual but not identically self-dual::

        sage: import random
        sage: r = random.choice(range(3, 10))
        sage: M = matroids.Whirl(r)
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True

    Except for `\mathcal{W}^2`, which is isomorphic to `U_{2, 4}`, these
    matroids have non-transitive automorphism groups::

        sage: r = random.choice(range(3, 8))
        sage: M = matroids.Whirl(r)
        sage: M.automorphism_group().is_transitive()
        False

    .. TODO::

        Optional arguments ``ring`` and ``x``, such that the resulting matroid
        is represented over ``ring`` by a reduced matrix like
        ``[-1  0  x]``
        ``[ 1 -1  0]``
        ``[ 0  1 -1]``

    REFERENCES:

    [Oxl2011]_, p. 659-60.
    """
    A = Matrix(GF(3), r, 2 * r, sparse=True)
    for i in range(r):
        A[i, i] = 1
        A[i, r + i] = 1
        if i != 0:
            A[i, i + r - 1] = -1
        else:
            A[i, 2 * r - 1] = 1
    M = TernaryMatroid(A)
    M = _rename_and_relabel(M, f'Whirl({r})', groundset)
    return M


def Uniform(r, n, groundset=None):
    """
    Return the uniform matroid of rank `r` on `n` elements.

    All subsets of size `r` or less are independent; all larger subsets are
    dependent. Representable when the field is sufficiently large. The precise
    bound is the subject of the MDS conjecture from coding theory.

    INPUT:

    - ``r`` -- nonnegative integer; the rank of the uniform matroid
    - ``n`` -- nonnegative integer; the number of elements of the uniform
      matroid
    - ``groundset`` -- string (optional); the groundset of the matroid

    OUTPUT: the uniform matroid `U_{r,n}`

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

    Check that bug :issue:`15292` was fixed::

        sage: M = matroids.Uniform(4,4)
        sage: len(M.circuit_closures())
        0

    REFERENCES:

    [Oxl2011]_, p. 660.
    """
    E = range(n)
    if r < n:
        CC = {r: [E]}
    else:
        CC = {}
    M = Matroid(groundset=E, circuit_closures=CC)
    M = _rename_and_relabel(M, f'U({r}, {n})', groundset)
    return M


def PG(n, q, x=None, groundset=None):
    """
    Return the projective geometry of dimension `n` over the finite field
    of order `q`.

    INPUT:

    - ``n`` -- positive integer; the dimension of the projective space. This
      is one less than the rank of the resulting matroid.
    - ``q`` -- positive integer that is a prime power; the order of the
      finite field
    - ``x`` -- string (default: ``None``); the name of the generator of a
      non-prime field, used for non-prime fields. If not supplied, ``'x'`` is
      used.
    - ``groundset`` -- string (optional); the groundset of the matroid

    OUTPUT: a linear matroid whose elements are the points of `PG(n, q)`

    EXAMPLES::

        sage: M = matroids.PG(2, 2)
        sage: M.is_isomorphic(matroids.catalog.Fano())
        True
        sage: matroids.PG(5, 4, 'z').size() == (4^6 - 1) / (4 - 1)
        True
        sage: M = matroids.PG(4, 7); M
        PG(4, 7): Linear matroid of rank 5 on 2801 elements represented over the Finite Field
         of size 7

    REFERENCES:

    [Oxl2011]_, p. 660.
    """
    if x is None:
        x = 'x'
    F = GF(q, x)
    P = ProjectiveSpace(n, F)
    A = Matrix(F, [list(p) for p in list(P)]).transpose()
    M = Matroid(A)
    M = _rename_and_relabel(M, f'PG({n}, {q})', groundset)
    return M


def AG(n, q, x=None, groundset=None):
    r"""
    Return the affine geometry of dimension ``n`` over the finite field of
    order ``q``.

    The affine geometry can be obtained from the projective geometry by
    removing a hyperplane.

    INPUT:

    - ``n`` -- positive integer; the dimension of the projective space. This
      is one less than the rank of the resulting matroid.
    - ``q`` -- positive integer that is a prime power; the order of the
      finite field
    - ``x`` -- string (default: ``None``); the name of the generator of a
      non-prime field, used for non-prime fields. If not supplied, ``'x'`` is
      used.
    - ``groundset`` -- string (optional); the groundset of the matroid

    OUTPUT: a linear matroid whose elements are the points of `AG(n, q)`

    EXAMPLES::

        sage: M = matroids.AG(2, 3).delete(8)
        sage: M.is_isomorphic(matroids.catalog.AG23minus())
        True
        sage: matroids.AG(5, 4, 'z').size() == ((4 ^ 6 - 1) / (4 - 1) - (4 ^ 5 - 1)/(4 - 1))
        True
        sage: M = matroids.AG(4, 2); M
        AG(4, 2): Binary matroid of rank 5 on 16 elements, type (5, 0)

    REFERENCES:

    [Oxl2011]_, p. 661.
    """
    if x is None:
        x = 'x'
    F = GF(q, x)
    P = ProjectiveSpace(n, F)
    A = Matrix(
        F, [list(p) for p in list(P) if not list(p)[0] == 0]
    ).transpose()
    M = Matroid(A)
    M = _rename_and_relabel(M, f'AG({n}, {q})', groundset)
    return M


def Z(r, t=True, groundset=None):
    r"""
    Return the unique rank-`r` binary spike.

    Defined for all `r \ge 3`.

    INPUT:

    - ``r`` -- integer (`r \ge 3`); the rank of the spike
    - ``t`` -- boolean (default: ``True``); whether the spike is tipped
    - ``groundset`` -- string (optional); the groundset of the matroid

    OUTPUT: matroid; the unique rank-`r` binary spike (tipped or tipless)

    EXAMPLES::

        sage: matroids.Z(8)
        Z_8: Binary matroid of rank 8 on 17 elements, type (7, 1)
        sage: matroids.Z(9)
        Z_9: Binary matroid of rank 9 on 19 elements, type (9, None)
        sage: matroids.Z(20, False)
        Z_20\t: Binary matroid of rank 20 on 40 elements, type (20, 0)

    It holds that `Z_3 \setminus e \cong M(K4)`, for all `e`::

        sage: import random
        sage: Z3 = matroids.Z(3)
        sage: E = sorted(Z3.groundset()); E
        ['t', 'x1', 'x2', 'x3', 'y1', 'y2', 'y3']
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

    and `Z_4 \setminus e \cong S_8`, for all `e \neq t`::

        sage: Z4 = matroids.Z(4)
        sage: E = sorted(Z4.groundset())
        sage: E.remove('t')
        sage: e = random.choice(E)
        sage: S8 = matroids.catalog.S8()
        sage: Z4.delete(e).is_isomorphic(S8)
        True
        sage: Z4.delete('t').is_isomorphic(S8)
        False

    The tipless binary spike is self-dual; it is identically self-dual if and
    only if r is even. It also has a transitive automorphism group::

        sage: r = random.choice(range(3, 8))
        sage: Z = matroids.Z(r, False)
        sage: Z.is_isomorphic(Z.dual())
        True
        sage: Z.equals(Z.dual()) != (r % 2 == 1)  # XOR
        True
        sage: Z.automorphism_group().is_transitive()  # long time
        True

    REFERENCES:

    [Oxl2011]_, p. 661-2.
    """
    from sage.matrix.special import identity_matrix, ones_matrix
    Id = Matrix(GF(2), identity_matrix(r))
    J = Matrix(GF(2), ones_matrix(r))
    tip = Matrix(GF(2), ones_matrix(r, 1))
    A = Id.augment(J-Id).augment(tip)

    M = Matroid(A)
    X = [f'x{i}' for i in range(1, r + 1)]
    Y = [f'y{i}' for i in range(1, r + 1)]
    if t:
        M = M.relabel(X + Y + ['t'])
        M.rename(f'Z_{r}: ' + repr(M))
    else:
        M = M.delete(2 * r)
        M = M.relabel(X + Y)
        M.rename(f'Z_{r}\\t: ' + repr(M))
    M = _rename_and_relabel(M, groundset=groundset)
    return M


def Spike(r, t=True, C3=[], groundset=None):
    r"""
    Return a rank-`r` spike.

    Defined for all `r \ge 3`; a rank-`r` spike with tip `t` and legs `L_1,
    L_2, \ldots, L_r`, where `L_i = \{t, x_i , y_i\}`. Deleting `t` gives a
    tipless rank-`r` spike.

    The groundset is `E = \{t, x_1, x_2, \ldots, x_r, y_1, y_2, \ldots,
    y_r\}` with `r(E) = r`.

    The nonspanning circuits are `\{L_1, L_2, \ldots, L_r\}`, all sets of the
    form `(L_i \cup L_j) \setminus t` for `1 \le i < j \le r`, and some
    (possibly empty) collection `C_3` of sets of the form `\{z_1, z_2, \ldots,
    z_r\}` where `z_i \in \{x_i, y_i\}` for all `i`, and no two members of
    `C_3` have more than `r-2` common elements.

    INPUT:

    - ``r`` -- integer (`r \ge 3`); the rank of the spike
    - ``t`` -- boolean (default: ``True``); whether the spike is tipped
    - ``C3`` -- list (default: ``[]``); a list of extra nonspanning circuits.
      The default (i.e. the empty list) results in a free `r`-spike.
    - ``groundset`` -- string (optional); the groundset of the matroid

    OUTPUT: matroid; a rank-`r` spike (tipped or tipless)

    EXAMPLES::

        sage: M = matroids.Spike(3, False); M
        Free 3-spike\t: M \ {'t'}, where M is Matroid of rank 3 on 7 elements with 3
         nonspanning circuits
        sage: M.is_isomorphic(matroids.Uniform(3, 6))
        True
        sage: len(matroids.Spike(8).bases())
        4864
        sage: import random
        sage: r = random.choice(range(3, 20))
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

    REFERENCES:

    [Oxl2011]_, p. 662.
    """
    if not (r >= 3):
        raise ValueError("The r-spike is defined for r >= 3.")

    E = ['t']
    X, Y = [], []
    for i in range(1, r + 1):
        X.append(f'x{i}')
        Y.append(f'y{i}')
    E += X
    E += Y

    if C3 == [] and r > 3:
        # free spike (can be defined fast through circuit closures)
        lines = [['t', f'x{i}', f'y{i}'] for i in range(1, r + 1)]
        planes = [['t', f'x{i}', f'y{i}', f'x{j}', f'y{j}']
                  for i in range(1, r + 1) for j in range(i + 1, r + 1)]
        CC = {2: lines, 3: planes, r: [E]}
        M = Matroid(circuit_closures=CC)
    else:
        for S in C3:
            for xy in S:
                if xy not in X+Y:
                    raise ValueError(
                        "The sets in C3 must contain elements xi and yi only."
                    )
            for T in C3:
                if S != T and len(set(S).intersection(set(T))) > r - 2:
                    raise ValueError(
                        "Every pair of sets in C3 must not have more than "
                        + "r - 2 common elements."
                    )

        NSC = []  # nonspanning_circuits
        NSC += C3
        for i in range(1, r + 1):
            NSC += [['t', f'x{i}', f'y{i}']]
            for j in range(i + 1, r + 1):
                NSC += [[f'x{i}', f'y{i}', f'x{j}', f'y{j}']]

        M = Matroid(rank=r, nonspanning_circuits=NSC)

    free = "Free " if C3 == [] else ""
    tip = "" if t else "\\t"
    M = M if t else M.delete('t')
    M = _rename_and_relabel(M, f'{free}{r}-spike{tip}', groundset)
    return M


# Q_r(A)


def Theta(n, groundset=None):
    r"""
    Return the matroid `\Theta_n`.

    Defined for all `n \ge 2`. `\Theta_2 \cong U_{1,2} \bigoplus U_{1,2}` and
    `\Theta_3 \cong M(K_4)`.

    INPUT:

    - ``n`` -- integer (`n \ge 2`); the rank of the matroid
    - ``groundset`` -- string (optional); the groundset of the matroid

    OUTPUT: matroid (`\Theta_n`)

    EXAMPLES::

        sage: matroids.Theta(30)
        Theta_30: Matroid of rank 30 on 60 elements with 16270 circuits
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
        sage: n = random.choice(range(3, 7))
        sage: M = matroids.Theta(n)
        sage: M.is_isomorphic(M.dual()) and not M.equals(M.dual())
        True

    For `n \le 3`, its automorphism group is transitive, while for `n \ge 4` it is not::

        sage: n = random.choice(range(4, 8))
        sage: M = matroids.Theta(2 + n % 2)
        sage: M.automorphism_group().is_transitive()
        True
        sage: M = matroids.Theta(n)
        sage: M.automorphism_group().is_transitive()
        False

    REFERENCES:

    [Oxl2011]_, p. 663-4.
    """
    X = [f'x{i}' for i in range(n)]
    Y = [f'y{i}' for i in range(n)]

    import itertools
    C = []
    C += list(itertools.combinations(X, 3))
    for i in range(n):
        Yi = [Y[j] for j in range(len(Y)) if j != i]
        C += [Yi + [f'x{i}']]

    for u in range(n):
        for s in range(n):
            for t in range(s + 1, n):
                if u != s and u != t and s != t:
                    Yu = [Y[i] for i in range(len(Y)) if i != u]
                    C += [Yu + [f'x{s}'] + [f'x{t}']]

    M = Matroid(circuits=C)
    M = _rename_and_relabel(M, f'Theta_{n}', groundset)
    return M


def Psi(r, groundset=None):
    r"""
    Return the matroid `\Psi_r`.

    The rank-`r` free swirl; defined for all `r \ge 3`.

    INPUT:

    - ``r`` -- integer (`r \ge 3`); the rank of the matroid
    - ``groundset`` -- string (optional); the groundset of the matroid

    OUTPUT: matroid (`\Psi_r`)

    EXAMPLES::

        sage: matroids.Psi(7)
        Psi_7: Matroid of rank 7 on 14 elements with 105 nonspanning circuits

    The matroid `\Psi_r` is `3`-connected but, for all `r \ge 4`, not `4`-connected::

        sage: M = matroids.Psi(3)
        sage: M.is_4connected()
        True
        sage: import random
        sage: r = random.choice(range(4, 8))
        sage: M = matroids.Psi(r)
        sage: M.is_4connected()
        False

    `\Psi_3 \cong U_{3, 6}`::

        sage: M = matroids.Psi(3)
        sage: M.is_isomorphic(matroids.catalog.U36())
        True

    It is identically self-dual with a transitive automorphism group::

        sage: M = matroids.Psi(r)
        sage: M.equals(M.dual())
        True
        sage: M.automorphism_group().is_transitive()  # long time
        True

    REFERENCES:

    [Oxl2011]_, p. 664.
    """
    A = [f'a{i}' for i in range(0, r)]
    B = [f'b{i}' for i in range(0, r)]
    E = A + B

    def generate_binary_strings(bit_count):
        binary_strings = []

        def genbin(n, bs=""):
            if len(bs) == n:
                binary_strings.append(bs)
            else:
                genbin(n, bs + 'a')
                genbin(n, bs + 'b')

        genbin(bit_count)
        return binary_strings

    NSC = []  # nonspanning circuits
    for i in range(0, r):
        for k in range(1, r - 2):
            I0 = [f'a{i}', f'b{i}']
            IK = [f'a{(i+k) % r}', f'b{(i+k) % r}']
            for AB in generate_binary_strings(k - 1):
                C = []
                C += I0 + IK
                j = 1
                for z in AB:
                    C += [f'{z}{(i+j) % r}']
                    j += 1
                NSC += [C]

    M = Matroid(groundset=E, rank=r, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, f'Psi_{r}', groundset)
    return M


# ******************************* #
#  Brettell's matroid collection  #
#                                 #
# ******************************* #


# 7 elements:


def RelaxedNonFano(groundset=None):
    """
    Return the relaxed NonFano matroid.

    An excluded minor for `2`-regular matroids. UPF is `K_2`.

    EXAMPLES::

        sage: M = matroids.catalog.RelaxedNonFano(); M
        F7=: Quaternary matroid of rank 3 on 7 elements
        sage: M.is_valid()
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(GF4, [[1, 1, 0, 1], [1, 0, 1, 1], [0, 1, w, 1]])
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "F7=", groundset)
    return M


def TippedFree3spike(groundset=None):
    """
    Return the tipped free `3`-spike.

    Unique 3-connected extension of
    :func:`U36 <sage.matroids.database_matroids.U36>`. Stabilizer for `K_2`.

    EXAMPLES::

        sage: M = matroids.catalog.TippedFree3spike(); M
        Tipped rank-3 free spike: Quaternary matroid of rank 3 on 7 elements
        sage: M.has_minor(matroids.Uniform(3,6))
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(GF4, [[1, 1, 1, 1], [1, w + 1, 0, w], [1, 0, w + 1, w]])
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[0, 3, 5, 1, 4, 6, 2]
    )
    M = _rename_and_relabel(M, "Tipped rank-3 free spike", groundset)
    return M


# 8 elements:


def AG23minusDY(groundset=None):
    r"""
    Return the matroid `AG23minusDY`.

    The matroid obtained from a `AG(2, 3)\setminus e` by a single `\delta-Y`
    exchange on a triangle. An excluded minor for near-regular matroids. UPF
    is `S`.

    EXAMPLES::

        sage: M = matroids.catalog.AG23minusDY(); M
        Delta-Y of AG(2,3)\e: Ternary matroid of rank 4 on 8 elements, type 0-
        sage: M.is_valid()
        True
    """
    A = Matrix(GF(3), [[1, 1, 1, 1], [1, 0, 1, 2], [2, 0, 1, 2], [2, 1, 1, 0]])
    M = TernaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "Delta-Y of AG(2,3)\\e", groundset)
    return M


def TQ8(groundset=None):
    """
    Return the matroid `TQ8`.

    An excluded minor for `2`-regular matroids. UPF is `K_2`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.TQ8(); M
        TQ8: Quaternary matroid of rank 4 on 8 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4, [[0, w, 1, 1], [1, 0, w, w + 1], [1, w, 0, w], [1, w + 1, 1, 0]]
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 7, 5, 3, 8, 6, 4, 2]
    )
    M = _rename_and_relabel(M, "TQ8", groundset)
    return M


def P8p(groundset=None):
    """
    Return the matroid `P8^-`.

    `P8^-` is obtained by relaxing one of the disjoint circuit-hyperplanes of
    :func:`P8 <sage.matroids.database_matroids.P8>`. An excluded minor for
    `2`-regular matroids. UPF is `K_2`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.P8p(); M
        P8-: Quaternary matroid of rank 4 on 8 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4, [[1, 1, 1, w], [1, w + 1, 1, 0], [1, 0, w, w], [0, 1, 1, 1]]
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=['a', 'c', 'b', 'f', 'd', 'e', 'g', 'h']
    )
    M = _rename_and_relabel(M, "P8-", groundset)
    return M


def KP8(groundset=None):
    """
    Return the matroid `KP8`.

    An excluded minor for `K_2`-representable matroids. UPF is `G`. Self-dual.
    Uniquely `GF(5)`-representable. (An excluded minor for `H_2`-representable
    matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.KP8(); M
        KP8: Quaternary matroid of rank 4 on 8 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4, [[0, 1, 1, 1], [1, 0, w, w], [1, 1, 1, 1 + w], [1, 1, 1 + w, 0]]
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 4, 3, 5, 6, 7, 0, 2]
    )
    M = _rename_and_relabel(M, "KP8", groundset)
    return M


def Sp8(groundset=None):
    """
    Return the matroid `Sp8`.

    An excluded minor for `G`- and `K_2`-representable matroids. UPF is
    `U_1^{(2)}`. Self-dual. (An excluded minor for `H_2`- and
    `GF(5)`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.Sp8(); M
        Sp8: Quaternary matroid of rank 4 on 8 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4, [[1, 1, w + 1, 0], [1, 1, 0, w + 1], [1, 0, w, w], [0, 1, 1, 1]]
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 2, 3, 5, 4, 6, 7, 8]
    )
    M = _rename_and_relabel(M, "Sp8", groundset)
    return M


def Sp8pp(groundset=None):
    """
    Return the matroid `Sp8=`.

    An excluded minor for `G`- and `K_2`-representable matroids. UPF is
    `(GF(2)(a,b),<a,b,a+1,b+1,ab+a+b>)`. Self-dual. (An excluded minor for
    `H_2`- and `GF(5)`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.Sp8pp(); M
        Sp8=: Quaternary matroid of rank 4 on 8 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(GF4, [[1, w, 1, 0], [1, 1, 1, 1], [w, 0, 1, w], [0, w, 1, 1]])
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 5, 6, 7, 2, 3, 4, 8]
    )
    M = _rename_and_relabel(M, "Sp8=", groundset)
    return M


def LP8(groundset=None):
    """
    Return the matroid `LP8`.

    An excluded minor for `G`- and `K_2`-representable matroids. Self-dual.
    UPF is `W`. (Also an excluded minor for `H_2`- and `GF(5)`-representable
    matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.LP8(); M
        LP8: Quaternary matroid of rank 4 on 8 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4, [[1, 1, 1, 1], [w + 1, w, 0, 1], [1, 0, w + 1, 1], [0, w, w, 1]]
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=['a', 'b', 'd', 'e', 'c', 'f', 'g', 'h']
    )
    M = _rename_and_relabel(M, "LP8", groundset)
    return M


def WQ8(groundset=None):
    r"""
    Return the matroid `WQ8`.

    An excluded minor for `G`, `K_2`, `H_4`, and `GF(5)`-representable
    matroids. Self-dual. UPF is `(Z[\zeta,a], <\zeta,a-\zeta>)` where `\zeta`
    is solution to `x^2-x+1 = 0` and `a` is an indeterminate.

    EXAMPLES::

        sage: M = matroids.catalog.WQ8(); M
        WQ8: Quaternary matroid of rank 4 on 8 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4, [[1, 0, 1, w + 1], [1, 1, 1, 1], [w, 1, 1, 0], [0, w, 1, 1]]
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[0, 1, 3, 4, 2, 5, 6, 7]
    )
    M = _rename_and_relabel(M, "WQ8", groundset)
    return M


# 9 elements:


def BB9(groundset=None):
    """
    Return the matroid `BB9`.

    An excluded minor for `K_2`-representable matroids, and a restriction of
    the Betsy Ross matroid. The UPF is `G`. Uniquely `GF(5)`-representable.
    (An excluded minor for `H_2`-representable matroids.)

    EXAMPLES::

        sage: BB = matroids.catalog.BB9(); BB
        BB9: Quaternary matroid of rank 3 on 9 elements
        sage: BR = matroids.catalog.BetsyRoss()
        sage: from itertools import combinations
        sage: pairs = combinations(sorted(BR.groundset()), 2)
        sage: for pair in pairs:
        ....:     if BR.delete(pair).is_isomorphic(BB):
        ....:         print(pair)
        ('a', 'h')
        ('b', 'i')
        ('c', 'j')
        ('d', 'f')
        ('e', 'g')
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4, [[1, 0, 1, 1, 1, 1],
              [0, 1, w, 1, 0, w],
              [w + 1, 1, w + 1, 1, w, 0]]
    )
    M = QuaternaryMatroid(
        reduced_matrix=A,
        groundset=['i', 'b', 'd', 'j', 'h', 'f', 'c', 'a', 'k']
    )
    M = _rename_and_relabel(M, "BB9", groundset)
    return M


def TQ9(groundset=None):
    """
    Return the matroid `TQ9`.

    An excluded minor for `K_2`-representable matroids, and a single-element
    extension of :func:`TQ8 <sage.matroids.database_matroids.TQ8>`.
    The UPF is `G`. Uniquely `GF(5)`-representable. (An excluded minor for
    `H_2`-representable matroids.)

    EXAMPLES::

        sage: TQ8 = matroids.catalog.TQ8()
        sage: TQ9 = matroids.catalog.TQ9(); TQ9
        TQ9: Quaternary matroid of rank 4 on 9 elements
        sage: for M in TQ8.extensions():
        ....:     if M.is_isomorphic(TQ9):
        ....:         print(True)
        ....:         break
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4, [[1, 0, w, 1, 1],
              [w + 1, 0, 0, w, 1],
              [1, w, 0, 0, w + 1],
              [1, 1, 1, 1, 0]],
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 4, 6, 0, 2, 5, 3, 7, 8]
    )
    M = _rename_and_relabel(M, "TQ9", groundset)
    return M


def TQ9p(groundset=None):
    """
    Return the matroid `TQ9^-`.

    An excluded minor for `G`- and `K_2`-representable matroids, and a
    single-element extension of
    :func:`TQ8 <sage.matroids.database_matroids.TQ8>`. UPF is
    `U_1^{(2)}`. (An excluded minor for `H_2`- and `GF(5)`-representable
    matroids.)

    EXAMPLES::

        sage: TQ8 = matroids.catalog.TQ8()
        sage: TQ9p = matroids.catalog.TQ9p(); TQ9p
        TQ9': Quaternary matroid of rank 4 on 9 elements
        sage: for M in TQ8.extensions():
        ....:     if M.is_isomorphic(TQ9p):
        ....:         print(True)
        ....:         break
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 1, 1, w + 1, 1],
            [w + 1, w, w + 1, w + 1, w + 1],
            [w, 0, w, 1, w + 1],
            [0, 1, w + 1, w + 1, 0],
        ],
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 4, 7, 8, 0, 6, 5, 2, 3]
    )
    M = _rename_and_relabel(M, "TQ9'", groundset)
    return M


def M8591(groundset=None):
    r"""
    Return the matroid `M8591`.

    An excluded minor for `K_2`-representable matroids. A `Y-\delta` exchange
    on the unique triad gives :func:`A9 <sage.matroids.database_matroids.A9>`.
    The UPF is `P_4`.

    EXAMPLES::

        sage: M = matroids.catalog.M8591(); M
        M8591: Quaternary matroid of rank 4 on 9 elements
        sage: M.is_valid()
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4, [[1, 1, 0, w, 1],
              [0, 1, 1, w, w + 1],
              [1, 0, w, w, 1],
              [0, 0, 1, 1, 0]]
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "M8591", groundset)
    return M


def PP9(groundset=None):
    """
    Return the matroid `PP9`.

    An excluded minor for `K_2`-representable matroids. A single-element
    extension of `P8^-`. The UPF is `P_4`. Has a
    :func:`P8p <sage.matroids.database_matroids.P8p>`-minor (delete `z`).
    Uniquely `GF(5)`-representable. (An excluded minor for `H_2`-representable
    matroids.)

    EXAMPLES::

        sage: P8p = matroids.catalog.P8p()
        sage: PP9 = matroids.catalog.PP9(); PP9
        PP9: Quaternary matroid of rank 4 on 9 elements
        sage: for M in P8p.extensions():
        ....:     if M.is_isomorphic(PP9):
        ....:         print(True)
        ....:         break
        True
        sage: M = PP9.delete('z')
        sage: M.is_isomorphic(P8p)
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4, [[1, 1, 1, w, w],
              [1, 1 + w, 1, 0, w],
              [1, 0, w, w, w],
              [0, 1, 1, 1, 1]]
    )
    M = QuaternaryMatroid(
        reduced_matrix=A,
        groundset=['a', 'c', 'b', 'f', 'd', 'e', 'g', 'h', 'z']
    )
    M = _rename_and_relabel(M, "PP9", groundset)
    return M


def BB9gDY(groundset=None):
    r"""
    Return the matroid `BB9gDY`.

    An excluded minor for `K_2`-representable matroids. The UPF is `G`. In a
    `DY^*`-equivalence class of 4 matroids, one of which can be obtained from
    :func:`BB9 <sage.matroids.database_matroids.BB9>` by a segment-cosegment
    exchange on `\{a,d,i,j\}`. Uniquely `GF(5)`-representable. (An excluded
    minor for `H_2`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.BB9gDY(); M
        Segment cosegment exchange on BB9: Quaternary matroid of rank 5 on 9 elements
        sage: M.is_valid()
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w, w, w + 1, 1],
            [w, 1, 0, 0],
            [w + 1, 1, 0, 0],
            [1, w, w + 1, w],
            [w, 0, 1, 1],
        ],
    )
    # M9573
    M = QuaternaryMatroid(
        reduced_matrix=A,
        groundset=['c', 'd', 'i', 'f', 'h', 'a', 'j', 'k', 'b']
    )
    M = _rename_and_relabel(M, "Segment cosegment exchange on BB9", groundset)
    return M


def A9(groundset=None):
    """
    Return the matroid `A9`.

    An excluded minor for `K_2`-representable matroids. The UPF is `P_4`.
    Uniquely `GF(5)`-representable. (An excluded minor for `H_2`-representable
    matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.A9(); M
        A9: Quaternary matroid of rank 3 on 9 elements
        sage: M.is_valid()
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4, [[w + 1, 1, w, w, w, w],
              [0, 1, 1, w + 1, 0, w],
              [w, 0, 1, w + 1, w, 1]]
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[6, 5, 4, 1, 2, 3, 7, 8, 0]
    )
    M = _rename_and_relabel(M, "A9", groundset)
    return M


def FN9(groundset=None):
    """
    Return the matroid `FN9`.

    An excluded minor for `G`- and `K_2`-representable matroids. In a
    `DY^*`-equivalence class of `10` matroids. UPF is `U_1^{(2)}`. (An
    excluded minor for `H_2`- and `GF(5)`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.FN9(); M
        FN9: Quaternary matroid of rank 3 on 9 elements
        sage: M.is_valid()
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w + 1, w, w + 1, w, 1, 0],
            [1, w + 1, 0, 1, w + 1, 1],
            [w + 1, w + 1, w, w + 1, 1, 1],
        ],
    )
    # M3209
    M = QuaternaryMatroid(
        reduced_matrix=A,
        groundset=['b0', 'a', 'y', 'z', 'x', "c0", 'b', 'c', 'a0']
    )
    M = _rename_and_relabel(M, "FN9", groundset)
    return M


def FX9(groundset=None):
    """
    Return the matroid `FX9`.

    An excluded minor for `G`- and `K_2`-representable matroids. UPF is
    `(Q(a,b), <-1,a,b,a-1,b-1,a-b,a+b,a+b-2,a+b-2ab>)`. (An excluded minor for
    `H_2`- and `GF(5)`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.FX9(); M
        FX9: Quaternary matroid of rank 4 on 9 elements
        sage: M.is_valid()
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, w + 1, 0, w, 1],
            [1, w, 1, w + 1, 1],
            [w + 1, w + 1, w, w + 1, w + 1],
            [w, w, w + 1, w + 1, 1],
        ],
    )
    # M48806
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FX9", groundset)
    return M


def KR9(groundset=None):
    """
    Return the matroid `KR9`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). In a `DY`-equivalence class of `4`
    matroids. Has a :func:`KP8 <sage.matroids.database_matroids.KP8>`-minor
    (delete `8`). UPF is `GF(4)`.

    EXAMPLES::

        sage: KR9 = matroids.catalog.KR9(); KR9
        KR9: Quaternary matroid of rank 4 on 9 elements
        sage: KP8 = matroids.catalog.KP8()
        sage: KP8.is_isomorphic(KR9.delete(8))
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w + 1, w, w + 1, w, 1],
            [0, 1, w + 1, w, 1],
            [w, w + 1, 1, 1, 1],
            [w + 1, w + 1, w, w + 1, 0],
        ],
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[2, 4, 0, 6, 1, 5, 3, 7, 8]
    )
    M = _rename_and_relabel(M, "KR9", groundset)
    return M


def KQ9(groundset=None):
    """
    Return the matroid `KQ9`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Has a
    :func:`TQ8 <sage.matroids.database_matroids.TQ8>`-minor` (delete `6`) and a
    :func:`KP8 <sage.matroids.database_matroids.KP8>`-minor (delete `8`). UPF
    is `GF(4)`.

    EXAMPLES::

        sage: KQ9 = matroids.catalog.KQ9(); KQ9
        KQ9: Quaternary matroid of rank 4 on 9 elements
        sage: TQ8 = matroids.catalog.TQ8()
        sage: TQ8.is_isomorphic(KQ9.delete(6))
        True
        sage: KP8 = matroids.catalog.KP8()
        sage: KP8.is_isomorphic(KQ9.delete(8))
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w + 1, w, w + 1, 1, w + 1],
            [1, 0, w + 1, w + 1, 1],
            [0, 1, w, w + 1, 1],
            [1, 1, w + 1, 0, w + 1],
        ],
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[5, 0, 4, 3, 2, 6, 8, 7, 1]
    )
    M = _rename_and_relabel(M, "KQ9", groundset)
    return M


# 10 elements:


def UG10(groundset=None):
    """
    Return the matroid `UG10`.

    An excluded minor for `K_2`- and `P_4`-representable matroids. Self-dual.
    An excluded minor for `H_3`- and `H_2`-representable matroids. Uniquely
    `GF(5)`-representable. Although not `P_4`-representable, it is
    `O`-representable, and hence is representable over all fields of size at
    least four.

    EXAMPLES::

        sage: M = matroids.catalog.UG10(); M
        UG10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 0, 1, w, w + 1],
            [1, 1, 1, 1, 1],
            [1, 0, w, 1, w + 1],
            [1, w + 1, w, 1, 0],
            [1, 1, 1, 0, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "UG10", groundset)
    return M


def FF10(groundset=None):
    """
    Return the matroid `FF10`.

    An excluded minor for `K_2`-representable matroids. UPF is `P_4`.
    Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FF10(); M
        FF10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 0, 1, 1, 1],
            [1, 0, 0, 1, 1],
            [1 + w, 1, 0, 0, 1],
            [1, 1, w, 0, 1],
            [1 + w, 1 + w, w, w, 1 + w],
        ],
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    )
    M = _rename_and_relabel(M, "FF10", groundset)
    return M


def GP10(groundset=None):
    """
    Return the matroid `GP10`.

    An excluded minor for `K_2`-representable matroids. UPF is `G`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.GP10(); M
        GP10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w + 1, w, 0, 1, 1],
            [w, w, 0, 0, 1],
            [0, 0, w + 1, w, 1],
            [1, 0, w, 0, 1],
            [1, 1, 1, 1, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "GP10", groundset)
    return M


def FZ10(groundset=None):
    """
    Return the matroid `FZ10`.

    An excluded minor for `K_2`- and `G`-representable matroids (and `H_2`-
    and `GF(5)`-representable matroids). UPF is `W`. Not self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FZ10(); M
        FZ10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 0, 1, w, 1],
            [1, 0, w, w + 1, w + 1],
            [1, 1, w, w, w],
            [0, 1, 0, w, w],
            [1, 1, w + 1, 1, w],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FZ10", groundset)
    return M


def UQ10(groundset=None):
    """
    Return the matroid `UQ10`.

    An excluded minor for `K_2`- and `G`-representable matroids (and `H_2`-
    and `GF(5)`-representable matroids). Self-dual. UPF is
    `(Q(a,b), <-1,a,b,a-1,b-1,a-b,a+b,a+1,ab+b-1,ab-b+1>)`.

    EXAMPLES::

        sage: M = matroids.catalog.UQ10(); M
        UQ10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 0, 1, w, 1],
            [1, 1, 1, 1, 0],
            [1, w + 1, 0, 1, 0],
            [w + 1, w, 0, 0, 1],
            [1, 1, 1, w + 1, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "UQ10", groundset)
    return M


def FP10(groundset=None):
    """
    Return the matroid `FP10`.

    An excluded minor for `K_2`- and `G`-representable matroids (and `H_2`-
    and `GF(5)`-representable matroids). UPF is `U_1^{(2)}`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FP10(); M
        FP10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 1, 1, 1, 1],
            [0, 1, 1, 1, 0],
            [w, w, w + 1, 0, 0],
            [0, 1, 1, 0, w],
            [1, w, 1, w, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FP10", groundset)
    return M


def TQ10(groundset=None):
    """
    Return the matroid `TQ10`.

    An excluded minor for `K_2`-representable matroids. UPF is `G`. Self-dual.
    Has :func:`TQ8 <sage.matroids.database_matroids.TQ8>` as a minor (delete
    'd' and contract 'c').

    EXAMPLES::

        sage: M = matroids.catalog.TQ10(); M
        TQ10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
        sage: N = M.delete('d').contract('c')
        sage: N.is_isomorphic(matroids.catalog.TQ8())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, w, 0, w + 1, 1],
            [w + 1, w + 1, w + 1, w, 1],
            [w + 1, w, 1, 0, 1],
            [1, 1, 0, w + 1, 0],
            [w + 1, 0, w, w + 1, 0],
        ],
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 6, 8, 'c', 3, 7, 'd', 2, 5, 4]
    )
    M = _rename_and_relabel(M, "TQ10", groundset)
    return M


def FY10(groundset=None):
    """
    Return the matroid `FY10`.

    An excluded minor for `P_4`-representable matroids. UPF is `G`. Not
    self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FY10(); M
        FY10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 1, 1, 0, 0],
            [0, 1, w + 1, 1, 0],
            [1, 1, 1, w + 1, 1],
            [0, 0, w, 1, 1],
            [1, w + 1, 1, 1, w],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FY10", groundset)
    return M


def PP10(groundset=None):
    """
    Return the matroid `PP10`.

    An excluded minor for `P_4`-representable matroids. UPF is `U_1^{(2)}`.
    Has a :func:`TQ8 <sage.matroids.database_matroids.TQ8>`-minor (e.g. delete
    'a' and contract 'e') and a
    :func:`PP9 <sage.matroids.database_matroids.PP9>` (and hence
    :func:`P8p <sage.matroids.database_matroids.P8p>`) minor (contract 'x').

    EXAMPLES::

        sage: PP10 = matroids.catalog.PP10(); PP10
        PP10: Quaternary matroid of rank 5 on 10 elements
        sage: M = PP10.delete('a').contract('e')
        sage: M.is_isomorphic(matroids.catalog.TQ8())
        True
        sage: M = PP10.contract('x')
        sage: M.is_isomorphic(matroids.catalog.PP9())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w + 1, 0, w + 1, 0, w],
            [w, w, 1, w, 1],
            [w + 1, w + 1, 0, w + 1, 1],
            [1, 0, 1, w + 1, 1],
            [w, 1, w + 1, w, w + 1],
        ],
    )
    M = QuaternaryMatroid(
        reduced_matrix=A,
        groundset=['z', 'f', 'c', 'g', 'e', 'b', 'a', 'h', 'd', 'x']
    )
    M = _rename_and_relabel(M, "PP10", groundset)
    return M


def FU10(groundset=None):
    """
    Return the matroid `FU10`.

    An excluded minor for `P_4`-representable matroids. UPF is `G`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FU10(); M
        FU10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 1, w, 1, 1],
            [w, w + 1, 1, 1, w],
            [1, 1, 0, w + 1, 0],
            [1, 1, 0, w, w + 1],
            [w + 1, w, w + 1, w, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FU10", groundset)
    return M


def D10(groundset=None):
    """
    Return the matroid `D10`.

    An excluded minor for `P_4`-representable matroids. UPF is `G`. Has a
    :func:`TQ8 <sage.matroids.database_matroids.TQ8>`-minor. In a
    `DY^*`-equivalence class of `13` matroids.

    EXAMPLES::

        sage: M = matroids.catalog.D10(); M
        D10: Quaternary matroid of rank 4 on 10 elements
        sage: M.has_minor(matroids.catalog.TQ8())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w, 1, w, 1, w + 1, w],
            [w, 0, w + 1, w + 1, w, w],
            [w + 1, 0, 0, w + 1, w + 1, w + 1],
            [w + 1, 1, 0, 1, w, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "D10", groundset)
    return M


def UK10(groundset=None):
    """
    Return the matroid `UK10`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Not self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.UK10(); M
        UK10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, w, w + 1, w, w + 1],
            [1, w, w, w + 1, w],
            [1, w + 1, w + 1, 0, w],
            [1, 1, w, 0, w],
            [w + 1, 0, 0, 1, w],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "UK10", groundset)
    return M


def PK10(groundset=None):
    """
    Return the matroid `PK10`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Not self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.PK10(); M
        PK10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 1, w, 0, w + 1],
            [0, 1, 1, 1, 1],
            [w + 1, w, w, w + 1, w],
            [0, 1, w, w + 1, 1],
            [w + 1, w + 1, 0, 0, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "PK10", groundset)
    return M


def GK10(groundset=None):
    """
    Return the matroid `GK10`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Not self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.GK10(); M
        GK10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 1, 0, 1, 1],
            [1, 0, 0, w, w],
            [w, w + 1, 1, 1, 0],
            [w, w + 1, w + 1, w + 1, w + 1],
            [1, w, w, w + 1, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "GK10", groundset)
    return M


def FT10(groundset=None):
    """
    Return the matroid `FT10`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.FT10(); M
        FT10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w, 0, w, w + 1, w + 1],
            [0, 1, w + 1, w, 1],
            [w, 1, 0, w + 1, w + 1],
            [w, 1, 0, 0, w],
            [0, 1, w, 0, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FT10", groundset)
    return M


def TK10(groundset=None):
    """
    Return the matroid `TK10`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.TK10(); M
        TK10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 1, 0, 0, w],
            [0, w, 0, w + 1, 1],
            [w + 1, 0, w + 1, w, w + 1],
            [w, 0, 1, 0, w],
            [0, w, 1, w, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "TK10", groundset)
    return M


def KT10(groundset=None):
    """
    Return the matroid `KT10`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.KT10(); M
        KT10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 1, 1, w, 0],
            [0, 1, 0, w + 1, 1],
            [w + 1, 0, 1, w + 1, w + 1],
            [w, w + 1, w + 1, 1, w],
            [w + 1, w + 1, w + 1, 1, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "KT10", groundset)
    return M


def TU10(groundset=None):
    """
    Return the matroid `TU10`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.TU10(); M
        TU10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, w + 1, 1, 0, w],
            [w, 0, 1, w + 1, w + 1],
            [w, w, w + 1, w + 1, 1],
            [w + 1, 0, w + 1, 1, w + 1],
            [w + 1, w + 1, 1, w, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "TU10", groundset)
    return M


def UT10(groundset=None):
    """
    Return the matroid `UT10`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Self-dual. UPF is `I`.

    EXAMPLES::

        sage: M = matroids.catalog.UT10(); M
        UT10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w, w + 1, 0, w + 1, 0],
            [w + 1, w + 1, 1, w, 1],
            [1, 0, 1, w + 1, w + 1],
            [w + 1, w + 1, 1, 1, w],
            [1, w + 1, 1, w, w + 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "UT10", groundset)
    return M


def FK10(groundset=None):
    """
    Return the matroid `FK10`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.FK10(); M
        FK10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 0, 1, 1, w],
            [w, 1, w, 1, 1],
            [1, 1, 0, 0, w + 1],
            [w + 1, 0, w + 1, 0, 1],
            [w + 1, w + 1, 1, 1, w + 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FK10", groundset)
    return M


def KF10(groundset=None):
    """
    Return the matroid `KF10`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.KF10(); M
        KF10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w + 1, w + 1, 1, w, 1],
            [0, w + 1, w, 1, 0],
            [0, w, 0, 1, w + 1],
            [w + 1, 0, w, w + 1, 1],
            [w + 1, 1, 0, w + 1, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "KF10", groundset)
    return M


# 11 elements:


def FA11(groundset=None):
    """
    Return the matroid `FA11`.

    An excluded minor for `P_4`-representable matroids. UPF is `PT`. In a
    `DY^*`-equivalence class of `6` matroids. Has an
    :func:`FF10 <sage.matroids.database_matroids.FF10>`-minor (delete `10`).

    EXAMPLES::

        sage: FA11 = matroids.catalog.FA11(); FA11
        FA11: Quaternary matroid of rank 5 on 11 elements
        sage: FF10 = matroids.catalog.FF10()
        sage: FF10.is_isomorphic(FA11.delete(10))
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w, 0, w, w, 1, 0],
            [1, w, 0, w, w, 1],
            [0, w, 1, w + 1, 0, w],
            [0, w, 0, w + 1, 0, 1],
            [w + 1, w + 1, w + 1, 0, w, 0],
        ],
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 3, 4, 2, 8, 7, 9, 0, 5, 10, 6]
    )
    M = _rename_and_relabel(M, "FA11", groundset)
    return M


# 12 elements:


def FR12(groundset=None):
    """
    Return the matroid `FR12`.

    An excluded minor for `K_2`-representable matroids. UPF is `P_4`.
    Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FR12(); M
        FR12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 1, 1, 1, 0, 1],
            [1, 1, 1, 0, 1, 1],
            [1, 1, 0, 0, 0, 1],
            [1, 0, 0, 0, w + 1, 1],
            [0, 1, 0, w + 1, 0, 1],
            [1, 1, 1, 1, 1, w + 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FR12", groundset)
    return M


def GP12(groundset=None):
    """
    Return the matroid `GP12`.

    An excluded minor for `K_2`-representable matroids. UPF is `G`. Not
    self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.GP12(); M
        GP12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1 + w, 1 + w, w, 1 + w, 0, 0],
            [1, 0, w, 1 + w, 0, 0],
            [1, 0, 1, 0, 1 + w, 1 + w],
            [1, 1, 0, 1, w, w],
            [1, 0, 1, 1, 0, 1 + w],
            [1, 1, 1, 1, 1, 1 + w],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "GP12", groundset)
    return M


def FQ12(groundset=None):
    """
    Return the matroid `FQ12`.

    An excluded minor for `P_4`-representable matroids. UPF is `PT`. Has` a
    :func:`PP9 <sage.matroids.database_matroids.PP9>`-minor (contract `4` and
    `7`, delete `6`) and
    :func:`FF10 <sage.matroids.database_matroids.FF10>`-minor (contract 'c'
    and delete 'd').

    EXAMPLES::

        sage: FQ12 = matroids.catalog.FQ12(); FQ12
        FQ12: Quaternary matroid of rank 6 on 12 elements
        sage: PP9 = matroids.catalog.PP9()
        sage: PP9.is_isomorphic(FQ12.contract([4,7]).delete(6))
        True
        sage: FF10 = matroids.catalog.FF10()
        sage: FF10.is_isomorphic(FQ12.contract('c').delete('d'))
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 0, w, w, 1, 0],
            [0, 0, w + 1, w + 1, 1, 1],
            [1, 1, w, 1, 1, 1],
            [w, 0, 1, 1, 0, 0],
            [w, w, w + 1, w + 1, 1, 1],
            [0, 1, 1, w, w + 1, 1],
        ],
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[7, 4, 5, 9, 2, 1, 0, 6, 'd', 'c', 8, 3]
    )
    M = _rename_and_relabel(M, "FQ12", groundset)
    return M


def FF12(groundset=None):
    """
    Return the matroid `FF12`.

    An excluded minor for `P_4`-representable matroids. Self-dual. UPF is
    `(Q(a,b),<-1,a,b,a-2,a-1,a+1,b-1,ab-a+b,ab-a-b,ab-a-2b>)`. Has an
    :func:`FF10 <sage.matroids.database_matroids.FF10>`-minor (contract 'c'
    and delete 'd').

    EXAMPLES::

        sage: M = matroids.catalog.FF12(); M
        FF12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        True
        sage: FF10 = matroids.catalog.FF10()
        sage: FF10.is_isomorphic(M.contract('c').delete('d'))
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 1, 0, 1, 0, 0],
            [1, 1, w, 1, 0, 1],
            [w, w + 1, 1, 0, 1, 1],
            [1, 1, w, 1, w + 1, 0],
            [1, 1, 0, 0, w + 1, 0],
            [1, w + 1, 1, 0, w + 1, w],
        ],
    )
    M = QuaternaryMatroid(
        reduced_matrix=A, groundset=[0, 4, 'c', 3, 5, 'd', 8, 9, 2, 7, 1, 6]
    )
    M = _rename_and_relabel(M, "FF12", groundset)
    return M


def FZ12(groundset=None):
    """
    Return the matroid `FZ12`.

    An excluded minor for `K_2`- and `G`-representable matroids (and `H_2`-
    and `GF(5)`-representable matroids). UPF is `W`. Not self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FZ12(); M
        FZ12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, w + 1, 0, w + 1, 0, w + 1],
            [w + 1, w + 1, w, 0, w, 1],
            [w, 1, 0, w, 0, w],
            [w + 1, 1, 1, w + 1, 1, w],
            [w, w, 1, 0, 0, 0],
            [w + 1, 1, 0, w, 1, w],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FZ12", groundset)
    return M


def UQ12(groundset=None):
    """
    Return the matroid `UQ12`.

    An excluded minor for `K_2` and `G`-representable matroids (and `H2` and
    `GF(5)`-representable matroids). UPF is `P_{pappus}`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.UQ12(); M
        UQ12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 0, 0, w + 1, 1, 0],
            [w + 1, w, w + 1, 1, 1, 1],
            [w + 1, w, w + 1, w + 1, w, 1],
            [1, 0, 0, w, 1, w + 1],
            [1, 0, 1, w, 1, w],
            [1, 1, 1, w, 1, w + 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "UQ12", groundset)
    return M


def FP12(groundset=None):
    """
    Return the matroid `FP12`.

    An excluded minor for `K_2`- and `G`-representable matroids (and `H_2`-
    and `GF(5)`-representable matroids). UPF is `W`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FP12(); M
        FP12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, w + 1, 1, 0, 1, w],
            [0, w + 1, 1, 0, w + 1, 0],
            [w + 1, 1, w, 0, w + 1, w],
            [w, 1, w, 1, w, w + 1],
            [w + 1, 0, w + 1, w, w + 1, 0],
            [w, 0, w + 1, w + 1, w, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FP12", groundset)
    return M


def FS12(groundset=None):
    """
    Return the matroid `FS12`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Rank `5`. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.FS12(); M
        FS12: Quaternary matroid of rank 5 on 12 elements
        sage: M.rank()
        5
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 1, 1, 1, 1, 0, 1],
            [0, 1, w + 1, w, w + 1, 1, 1],
            [w, w + 1, 0, 0, w + 1, 0, 0],
            [0, 0, w + 1, 0, 0, 1, w],
            [1, 0, w, w, w, 1, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FS12", groundset)
    return M


def UK12(groundset=None):
    """
    Return the matroid `UK12`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Self-dual. UPF is `I`.

    EXAMPLES::

        sage: M = matroids.catalog.UK12(); M
        UK12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, w + 1, 1, 0, 0, w],
            [w, 1, w, w + 1, w, w + 1],
            [1, w + 1, w + 1, 0, 0, 0],
            [0, w, 1, 1, w + 1, 1],
            [0, w, 1, 1, 0, w + 1],
            [w, w, w + 1, w, w, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "UK12", groundset)
    return M


def UA12(groundset=None):
    """
    Return the matroid `UA12`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Not self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.UA12(); M
        UA12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, w + 1, w, w + 1, w + 1, 0],
            [1, 1, 1, w + 1, w, 1],
            [1, w + 1, 0, 1, w + 1, 1],
            [0, 0, 0, w + 1, w, 1],
            [0, w + 1, 0, w, w + 1, 1],
            [1, w, w + 1, w + 1, w, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "UA12", groundset)
    return M


def AK12(groundset=None):
    """
    Return the matroid `AK12`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Not self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.AK12(); M
        AK12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, w, 0, 0, w + 1, w],
            [w + 1, w, 0, w, w + 1, w + 1],
            [0, w, w + 1, 0, w, 0],
            [1, 0, w, 0, 0, w + 1],
            [0, 1, 0, 1, 1, 0],
            [1, 0, 1, 1, 0, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "AK12", groundset)
    return M


def FK12(groundset=None):
    """
    Return the matroid `FK12`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.FK12(); M
        FK12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w + 1, w, w, 0, w, w + 1],
            [w, 1, 1, w + 1, w + 1, w],
            [1, w + 1, w, 0, w + 1, 0],
            [w, 1, 1, w, w + 1, w + 1],
            [w + 1, w + 1, w, 0, w, 0],
            [1, w, w, w, 1, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FK12", groundset)
    return M


def KB12(groundset=None):
    """
    Return the matroid `KB12`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.KB12(); M
        KB12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 0, w, 0, 0, 1],
            [1, w, 0, w, 1, w + 1],
            [1, 1, w + 1, 0, 1, w + 1],
            [w, w, 1, w, 0, 0],
            [1, 1, 1, 0, 0, 0],
            [w, w, 1, w, 1, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "KB12", groundset)
    return M


def AF12(groundset=None):
    """
    Return the matroid `AF12`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.AF12(); M
        AF12: Quaternary matroid of rank 6 on 12 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 0, 0, 0, 1, 1],
            [0, 1, 0, w, w, w + 1],
            [0, 1, 0, 0, w + 1, w + 1],
            [w, 1, 1, 0, 0, 1],
            [0, 1, w + 1, w, w, 1],
            [1, 1, 1, w + 1, w + 1, w],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "AF12", groundset)
    return M


def NestOfTwistedCubes(groundset=None):
    r"""
    Return the NestOfTwistedCubes matroid.

    A matroid with no `U(2,4)`-detachable pairs (only `\{e_i,f_i\}` pairs are
    detachable).

    EXAMPLES::

        sage: M = matroids.catalog.NestOfTwistedCubes(); M
        NestOfTwistedCubes: Matroid of rank 6 on 12 elements with 57 circuits
        sage: M.is_3connected()
        True
    """
    # utility function
    def complement(groundset, subset):
        return list(set(groundset).difference(subset))

    gs = ["e1", "e2", "e3", "e4", "e5", "e6",
          "f1", "f2", "f3", "f4", "f5", "f6"]
    M = Matroid(
        groundset=gs,
        circuit_closures={
            3: [
                ["e1", "e2", "f3", "f4"],
                ["e3", "e4", "f5", "f6"],
                ["e5", "e6", "f1", "f2"],
                ["e4", "e5", "f1", "f3"],
                ["e1", "e3", "f6", "f2"],
                ["e6", "e2", "f4", "f5"],
                ["e2", "e5", "f3", "f6"],
                ["e1", "e4", "f2", "f5"],
                ["e3", "e6", "f1", "f4"],
                ["e5", "e1", "f4", "f6"],
                ["e4", "e6", "f2", "f3"],
                ["e2", "e3", "f5", "f1"],
                ["e2", "e4", "f6", "f1"],
                ["e6", "e1", "f3", "f5"],
                ["e3", "e5", "f2", "f4"],
            ],
            5: [
                complement(gs, ["e1", "e2", "f5", "f6"]),
                complement(gs, ["e3", "e4", "f1", "f2"]),
                complement(gs, ["e5", "e6", "f3", "f4"]),
                complement(gs, ["e4", "e5", "f6", "f2"]),
                complement(gs, ["e1", "e3", "f4", "f5"]),
                complement(gs, ["e6", "e2", "f1", "f3"]),
                complement(gs, ["e2", "e5", "f1", "f4"]),
                complement(gs, ["e1", "e4", "f3", "f6"]),
                complement(gs, ["e3", "e6", "f2", "f5"]),
                complement(gs, ["e5", "e1", "f2", "f3"]),
                complement(gs, ["e4", "e6", "f5", "f1"]),
                complement(gs, ["e2", "e3", "f4", "f6"]),
                complement(gs, ["e2", "e4", "f3", "f5"]),
                complement(gs, ["e6", "e1", "f2", "f4"]),
                complement(gs, ["e3", "e5", "f6", "f1"]),
            ],
            6: [gs],
        },
    )
    M = Matroid(circuits=list(M.circuits()))
    M = _rename_and_relabel(M, "NestOfTwistedCubes", groundset)
    return M


# 13 elements:


def XY13(groundset=None):
    """
    Return the matroid `XY13`.

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids). UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.XY13(); M
        XY13: Quaternary matroid of rank 6 on 13 elements
        sage: M.is_3connected()
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 0, 1, 1, 0, 1, 1],
            [w, 1, w, w + 1, 1, 1, w + 1],
            [0, 0, w + 1, 1, 1, w, 1],
            [0, w, 1, 1, w + 1, 1, 1],
            [w, w + 1, w, w + 1, 1, 1, w],
            [1, 0, 0, 1, 0, w, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "XY13", groundset)
    return M


# 14 elements:


def N3(groundset=None):
    """
    Return the matroid `N3`.

    An excluded minor for dyadic matroids (and `GF(5)`-representable matroids).
    UPF is `GF(3)`. `4`- (but not `5`-) connected. Self-dual.

    EXAMPLES::

        sage: N3 = matroids.catalog.N3(); N3
        N3: Ternary matroid of rank 7 on 14 elements, type 0+
        sage: N3.is_isomorphic(N3.dual())
        True
        sage: N3.is_kconnected(4)
        True
        sage: N3.is_kconnected(5)
        False
    """
    A = Matrix(
        GF(3),
        [
            [2, 0, 0, 2, 1, 1, 2],
            [1, 2, 0, 0, 2, 0, 2],
            [0, 1, 2, 1, 0, 0, 2],
            [0, 2, 2, 0, 0, 0, 2],
            [1, 0, 0, 0, 1, 0, 2],
            [2, 0, 0, 1, 2, 1, 1],
            [1, 1, 1, 2, 2, 2, 0],
        ],
    )
    M = TernaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "N3", groundset)
    return M


def N3pp(groundset=None):
    """
    Return the matroid `N3pp`.

    An excluded minor for `K_2`-representable matroids. Self-dual.
    Obtained by relaxing the two complementary circuit-hyperplanes of
    :func:`N4 <sage.matroids.database_matroids.N4>`. Not
    `P_4`-representable, but `O`-representable, and hence representable
    over all fields of size at least four.

    EXAMPLES::

        sage: M = matroids.catalog.N3pp(); M
        N3=: Quaternary matroid of rank 7 on 14 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 1, 0, 1, 0, w, w],
            [1, w, 1, 0, 1, 0, 1],
            [1, 0, 0, 0, w, w, w],
            [1, 0, 0, 0, 1, 0, 1],
            [1, 1, 1, 0, 1, 0, 0],
            [1, w, 1, w, 1, 1, 1],
            [1, 1, 0, 1, 0, 0, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "N3=", groundset)
    return M


def UP14(groundset=None):
    """
    Return the matroid `UP14`.

    An excluded minor for `K_2`-representable matroids. Has disjoint
    circuit-hyperplanes. UPF is `W`. Not self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.UP14(); M
        UP14: Quaternary matroid of rank 7 on 14 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, w, w, 1, w, 0, 0],
            [w, w, w, 1, w, w, 0],
            [w, w, 1, w, 1, 0, 1],
            [1, 1, 0, 1, 0, 0, 0],
            [w, w, 1, w, w, w, 1],
            [0, w, 1, 1, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "UP14", groundset)
    return M


def VP14(groundset=None):
    """
    Return the matroid `VP14`.

    An excluded minor for `K_2`-representable matroids. Has disjoint
    circuit-hyperplanes. UPF is `W`. Not self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.VP14(); M
        VP14: Quaternary matroid of rank 7 on 14 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [w, 0, 1, 1, w, 1, 1],
            [w, w, 1, 1, w, 1, w],
            [1, 1, 0, 1, 1, w, 0],
            [1, 1, 0, 1, 0, 0, 0],
            [0, 0, 0, w, 0, 1, 1],
            [1, 0, 0, w, 1, 0, 0],
            [1, 1, 1, 1, w, 1, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "VP14", groundset)
    return M


def FV14(groundset=None):
    """
    Return the matroid `FV14`.

    An excluded minor for `P_4`-representable matroids. Not self-dual. UPF is
    `PT`.

    EXAMPLES::

        sage: M = matroids.catalog.FV14(); M
        FV14: Quaternary matroid of rank 7 on 14 elements
        sage: M.is_isomorphic(M.dual())
        False
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 1, 0, 1, 1, 1, 1],
            [0, 1, 0, 1, 0, 1, 0],
            [1, 0, 0, 1, 1, 1, 1],
            [0, 1 - w, 1, 1, 1 - w, 1 - w, 1],
            [0, 0, 1, 0, 1 - w, -w, 1],
            [1 - w, 0, -1, -w, 0, 0, -w],
            [1, 0, 0, 1, 0, 0, 1],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FV14", groundset)
    return M


def OW14(groundset=None):
    """
    Return the matroid `OW14`.

    An excluded minor for `P_4`-representable matroids. Self-dual. UPF is
    `Orthrus`.

    EXAMPLES::

        sage: M = matroids.catalog.OW14(); M
        OW14: Quaternary matroid of rank 7 on 14 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [0, 1, 1, w, 1, 0, w],
            [0, w + 1, 1, w, w + 1, 0, w + 1],
            [0, w + 1, 1, w, 0, w + 1, w + 1],
            [1, 1, w + 1, w + 1, 0, 1, 1],
            [1, 0, 0, w, 0, 0, w],
            [1, 0, 1, 0, 1, 0, 0],
            [0, 1, 0, w, 0, 1, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "OW14", groundset)
    return M


def FM14(groundset=None):
    """
    Return the matroid `FM14`.

    An excluded minor for `P_4`-representable matroids. Self-dual. UPF is `PT`.

    EXAMPLES::

        sage: M = matroids.catalog.FM14(); M
        FM14: Quaternary matroid of rank 7 on 14 elements
        sage: M.is_isomorphic(M.dual())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 0, w, 1, 1, w, 0],
            [0, 1, 0, 0, 1, 1, 1],
            [0, 1, 0, 0, w, w, 1],
            [w, 1, 0, w, 0, 1, 1],
            [0, 1, 1, 0, 0, 0, 1],
            [1, w, w, 0, 0, 0, 0],
            [0, w, 0, 1, 0, w, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FM14", groundset)
    return M


# 15 elements:


def FA15(groundset=None):
    """
    Return the matroid `FA15`.

    An excluded minor for `O`-representable matroids. UPF is `PT`. In a
    `DY^*`-equivalence class of `6` matroids. Has an
    :func:`SQ14 <sage.matroids.database_matroids.N3pp>`-minor.

    EXAMPLES::

        sage: M = matroids.catalog.FA15(); M
        FA15: Quaternary matroid of rank 7 on 15 elements
        sage: M.has_minor(matroids.catalog.N3pp())
        True
    """
    GF4 = GF(4, 'w')
    w = GF4('w')
    A = Matrix(
        GF4,
        [
            [1, 1, 1, 1, 1, 1, 1, 1],
            [1, 0, 1, 1, 1, 1, 0, 1],
            [0, 0, 1, 0, 0, 1, 0, 1],
            [0, 0, w, w, 0, 1, 1, 1],
            [1, 1, 0, 0, 0, 0, w, w],
            [w, 0, 1, 1, 1, 1, 0, 0],
            [w, w, 1, 1, 0, 0, 0, 0],
        ],
    )
    M = QuaternaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "FA15", groundset)
    return M


# 16 elements:


def N4(groundset=None):
    """
    Return the matroid `N4`.

    An excluded minor for dyadic matroids (and `GF(5)`-representable matroids).
    UPF is `GF(3)`. `4`- (but not `5`-) connected. Self-dual.

    EXAMPLES::

        sage: N4 = matroids.catalog.N4(); N4
        N4: Ternary matroid of rank 8 on 16 elements, type 0+
        sage: N4.is_isomorphic(N4.dual())
        True
        sage: N4.is_kconnected(4)
        True
        sage: N4.is_kconnected(5)
        False
    """
    A = Matrix(
        GF(3),
        [
            [2, 0, 2, 1, 2, 1, 0, 0],
            [2, 1, 0, 2, 0, 2, 2, 0],
            [2, 0, 2, 0, 0, 1, 0, 0],
            [2, 0, 2, 2, 2, 2, 2, 2],
            [0, 1, 1, 1, 1, 1, 1, 1],
            [2, 1, 0, 0, 0, 2, 0, 0],
            [2, 0, 0, 2, 2, 2, 2, 0],
            [1, 0, 1, 2, 1, 2, 1, 1],
        ],
    )
    M = TernaryMatroid(reduced_matrix=A)
    M = _rename_and_relabel(M, "N4", groundset)
    return M


# ******************************** #
#  Collection of various matroids  #
#                                  #
# ******************************** #


def NonVamos(groundset=None):
    r"""
    Return the non-`V\acute{a}mos` matroid.

    The non-`V\acute{a}mos` matroid, or `V_8^+` is an `8`-element matroid of
    rank `4`. It is a tightening of the `V\acute{a}mos` matroid. It is
    representable over some field.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.catalog.NonVamos(); M
        NonVamos: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
             {'c', 'd', 'e', 'f'}, {'c', 'd', 'g', 'h'}, {'e', 'f', 'g', 'h'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: setprint(M.nonbases())
        [{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
         {'c', 'd', 'e', 'f'}, {'c', 'd', 'g', 'h'}, {'e', 'f', 'g', 'h'}]
        sage: M.is_dependent(['c', 'd', 'g', 'h'])
        True
        sage: M.is_valid()  # long time
        True

    REFERENCES:

    [Oxl2011]_, p. 72, 84.
    """
    CC = {
        3: ['abcd', 'abef', 'cdef', 'abgh', 'cdgh', 'efgh'],
        4: ['abcdefgh']
    }
    M = Matroid(circuit_closures=CC)
    M = _rename_and_relabel(M, "NonVamos", groundset)
    return M


def NotP8(groundset='abcdefgh'):
    """
    Return the matroid ``NotP8``.

    EXAMPLES::

        sage: M = matroids.catalog.P8()
        sage: N = matroids.catalog.NotP8(); N
        NotP8: Ternary matroid of rank 4 on 8 elements, type 0-
        sage: M.is_isomorphic(N)
        False
        sage: M.is_valid()
        True

    REFERENCES:

    [Oxl1992]_, p.512 (the first edition).
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 0, 1, 1, -1],
        [0, 1, 0, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 0, 1],
        [0, 0, 0, 1, -1, 1, 1, 1]
    ])
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "NotP8")
    return M


def AG23minus(groundset=None):
    """
    Return the ternary affine plane minus a point.

    This is a sixth-roots-of-unity matroid, and an excluded minor for the
    class of near-regular matroids.

    EXAMPLES::

        sage: M = matroids.catalog.AG23minus(); M
        AG23minus: Matroid of rank 3 on 8 elements with circuit-closures
        {2: {{'a', 'b', 'c'}, {'a', 'd', 'f'}, {'a', 'e', 'g'},
        {'b', 'd', 'h'}, {'b', 'e', 'f'}, {'c', 'd', 'g'},
        {'c', 'e', 'h'}, {'f', 'g', 'h'}},
        3: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: M.is_valid()
        True

    REFERENCES:

    [Oxl2011]_, p. 653.
    """
    CC = {2: ['abc', 'ceh', 'fgh', 'adf', 'aeg', 'cdg', 'bdh', 'bef'],
          3: ['abcdefgh']}
    M = Matroid(circuit_closures=CC)
    M = _rename_and_relabel(M, "AG23minus", groundset)
    return M


def P9(groundset='abcdefghi'):
    """
    Return the matroid `P_9`.

    EXAMPLES::

        sage: M = matroids.catalog.P9(); M
        P9: Binary matroid of rank 4 on 9 elements, type (1, 1)
        sage: M.is_valid()
        True

    REFERENCES:

    This is the matroid referred to as `P_9` by Oxley in his paper "The binary
    matroids with no 4-wheel minor", [Oxl1987]_.
    """
    A = Matrix(GF(2), [
        [1, 0, 0, 0, 1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1, 1, 0, 0, 1],
        [0, 0, 1, 0, 0, 1, 1, 0, 1],
        [0, 0, 0, 1, 0, 0, 1, 1, 0]
    ])
    M = BinaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "P9")
    return M


def R9A(groundset=None):
    """
    Return the matroid `R_9^A`.

    The matroid `R_9^A` is not representable over any field, yet none of the
    cross-ratios in its Tuttegroup equal 1. It is one of the 4 matroids on at
    most 9 elements with this property, the others being `{R_9^A}^*`, `R_9^B`
    and `{R_9^B}^*`.

    EXAMPLES::

        sage: M = matroids.catalog.R9A(); M
        R9A: Matroid of rank 4 on 9 elements with 13 nonspanning circuits
        sage: M.is_valid()
        True
    """
    NSC = ['abch', 'abde', 'abfi', 'acdi', 'aceg', 'adgh', 'aefh', 'bcdf',
           'bdhi', 'begi', 'cehi', 'defi', 'fghi']
    M = Matroid(rank=4, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, "R9A", groundset)
    return M


def R9B(groundset=None):
    """
    Return the matroid `R_9^B`.

    The matroid `R_9^B` is not representable over any field, yet none of the
    cross-ratios in its Tuttegroup equal 1.
    It is one of the 4 matroids on at most 9 elements with this property, the
    others being `{R_9^B}^*`, `R_9^A` and `{R_9^A}^*`.

    EXAMPLES::

        sage: M = matroids.catalog.R9B(); M
        R9B: Matroid of rank 4 on 9 elements with 13 nonspanning circuits
        sage: M.is_valid() and M.is_paving()
        True
    """
    NSC = ['abde', 'bcdf', 'aceg', 'abch', 'befh', 'cdgh', 'bcei', 'adfi',
           'abgi', 'degi', 'bdhi', 'aehi', 'fghi']
    M = Matroid(rank=4, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, "R9B", groundset)
    return M


def Block_9_4(groundset=None):
    """
    Return the paving matroid whose nonspanning circuits form the blocks of a
    `2-(9, 4, 3)` design.

    EXAMPLES::

        sage: M = matroids.catalog.Block_9_4(); M
        Block(9, 4): Matroid of rank 4 on 9 elements with 18 nonspanning
        circuits
        sage: M.is_valid() and M.is_paving()
        True
        sage: BD = BlockDesign(M.groundset(), list(M.nonspanning_circuits()))
        sage: BD.is_t_design(return_parameters=True)
        (True, (2, 9, 4, 3))
    """
    NSC = ['abcd', 'acef', 'bdef', 'cdeg', 'abfg', 'adeh', 'bcfh', 'acgh',
           'begh', 'dfgh', 'abei', 'cdfi', 'bcgi', 'adgi', 'efgi', 'bdhi',
           'cehi', 'afhi']
    M = Matroid(rank=4, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, "Block(9, 4)", groundset)
    return M


def TicTacToe(groundset=None):
    """
    Return the TicTacToe matroid.

    The dual of the TicTacToe matroid is not algebraic; it is unknown whether
    the TicTacToe matroid itself is algebraic.

    EXAMPLES::

        sage: M = matroids.catalog.TicTacToe(); M
        TicTacToe: Matroid of rank 5 on 9 elements with 8 nonspanning circuits
        sage: M.is_valid() and M.is_paving()
        True

    REFERENCES:

    [Hoc]_
    """
    NSC = ['abcdg', 'adefg', 'abceh', 'abcfi', 'cdefi', 'adghi', 'beghi',
           'cfghi']
    M = Matroid(rank=5, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, "TicTacToe", groundset)
    return M


def N1(groundset='abcdefghij'):
    r"""
    Return the matroid `N_1`, represented over `\GF{3}`.

    `N_1` is an excluded minor for the dyadic matroids.

    EXAMPLES::

        sage: M = matroids.catalog.N1(); M
        N1: Ternary matroid of rank 5 on 10 elements, type 0+
        sage: M.is_field_isomorphic(M.dual())
        True
        sage: M.is_valid()
        True

    REFERENCES:

    [Oxl2011]_, p. 554.
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 0, 2, 0, 0, 1, 1],
        [0, 1, 0, 0, 0, 1, 2, 0, 0, 1],
        [0, 0, 1, 0, 0, 0, 1, 2, 0, 1],
        [0, 0, 0, 1, 0, 0, 0, 1, 2, 2],
        [0, 0, 0, 0, 1, 1, 1, 1, 2, 0]
    ])
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "N1")
    return M


def Block_10_5(groundset=None):
    """
    Return the paving matroid whose nonspanning circuits form the blocks of a
    `3-(10, 5, 3)` design.

    EXAMPLES::

        sage: M = matroids.catalog.Block_10_5(); M
        Block(10, 5): Matroid of rank 5 on 10 elements with 36 nonspanning
        circuits
        sage: M.is_valid() and M.is_paving()
        True
        sage: BD = BlockDesign(M.groundset(), list(M.nonspanning_circuits()))
        sage: BD.is_t_design(return_parameters=True)
        (True, (3, 10, 5, 3))
    """
    NSC = ['abcde', 'acdfg', 'bdefg', 'bcdfh', 'abefh', 'abcgh', 'adegh',
           'cefgh', 'bcefi', 'adefi', 'bcdgi', 'acegi', 'abfgi', 'abdhi',
           'cdehi', 'acfhi', 'beghi', 'dfghi', 'abdfj', 'acefj', 'abegj',
           'cdegj', 'bcfgj', 'acdhj', 'bcehj', 'defhj', 'bdghj', 'afghj',
           'abcij', 'bdeij', 'cdfij', 'adgij', 'efgij', 'aehij', 'bfhij',
           'cghij']
    M = Matroid(rank=5, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, "Block(10, 5)", groundset)
    return M


def Q10(groundset='abcdefghij'):
    r"""
    Return the matroid `Q_{10}`, represented over `\GF{4}`.

    `Q_{10}` is a `10`-element, rank-`5`, self-dual matroid. It is
    representable over `\GF{3}` and `\GF{4}`, and hence is a
    sixth-roots-of-unity matroid. `Q_{10}` is a splitter for the class of
    sixth-root-of-unity matroids.

    EXAMPLES::

        sage: M = matroids.catalog.Q10(); M
        Q10: Quaternary matroid of rank 5 on 10 elements
        sage: M.is_isomorphic(M.dual())
        True
        sage: M.is_valid()
        True

    Check the splitter property. By Seymour's Theorem, and using self-duality,
    we only need to check that all 3-connected single-element extensions have
    an excluded minor for sixth-roots-of-unity. The only excluded minors that
    are quaternary are `U_{2, 5}, U_{3, 5}, F_7, F_7^*`. As it happens, it
    suffices to check for `U_{2, 5}`::

        sage: S = matroids.catalog.Q10().linear_extensions(simple=True)
        sage: [M for M in S if not M.has_line_minor(5)]
        []
    """
    F = GF(4, 'x')
    x = F.gens()[0]
    A = Matrix(F, [
        [1, 0, 0, 0, 0, 1, x, 0, 0, x + 1],
        [0, 1, 0, 0, 0, x + 1, 1, x, 0, 0],
        [0, 0, 1, 0, 0, 0, x + 1, 1, x, 0],
        [0, 0, 0, 1, 0, 0, 0, x + 1, 1, x],
        [0, 0, 0, 0, 1, x, 0, 0, x + 1, 1]
    ])
    M = QuaternaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "Q10")
    return M


def BetsyRoss(groundset=None):
    """
    Return the Betsy Ross matroid, represented by circuit closures.

    An extremal golden-mean matroid. That is, if `M` is simple, rank `3`, has
    the Betsy Ross matroid as a restriction and is a Golden Mean matroid, then
    `M` is the Betsy Ross matroid.

    EXAMPLES::

        sage: M = matroids.catalog.BetsyRoss(); M
        BetsyRoss: Matroid of rank 3 on 11 elements with 25 nonspanning
        circuits
        sage: len(M.circuit_closures()[2])
        10
        sage: M.is_valid()
        True
    """
    NSC = ['acf', 'acg', 'adi', 'adj', 'afg', 'ahk', 'aij', 'bdg', 'bdh',
           'bef', 'bej', 'bfj', 'bgh', 'bik', 'ceh', 'cei', 'cfg', 'chi',
           'cjk', 'dfk', 'dgh', 'dij', 'efj', 'egk', 'ehi']
    M = Matroid(rank=3, nonspanning_circuits=NSC)
    M = _rename_and_relabel(M, "BetsyRoss", groundset)
    return M


def N2(groundset='abcdefghijkl'):
    r"""
    Return the matroid `N_2`, represented over `\GF{3}`.

    `N_2` is an excluded minor for the dyadic matroids.

    EXAMPLES::

        sage: M = matroids.catalog.N2(); M
        N2: Ternary matroid of rank 6 on 12 elements, type 0+
        sage: M.is_field_isomorphic(M.dual())
        True
        sage: M.is_valid()
        True

    REFERENCES:

    [Oxl2011]_, p. 554.
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1, 1],
        [0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 0, 1],
        [0, 0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 1],
        [0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 1, 0],
        [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0, 1, 1, 2, 2, 1, 0, 1]
    ])
    M = TernaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "N2")
    return M


def D16(groundset='abcdefghijklmnop'):  # A.K.A. the Carolyn Chun Matroid
    """
    Return the matroid `D_{16}`.

    Let `M` be a `4`-connected binary matroid and `N` an internally
    `4`-connected proper minor of `M` with at least 7 elements. Then some
    element of `M` can be deleted or contracted preserving an `N`-minor,
    unless `M` is `D_{16}`.

    EXAMPLES::

        sage: M = matroids.catalog.D16(); M
        D16: Binary matroid of rank 8 on 16 elements, type (0, 0)
        sage: M.is_valid()
        True

    REFERENCES:

    [CMO2012]_
    """
    A = Matrix(GF(2), [
        [1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1],
        [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0]
    ])
    M = BinaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "D16")
    return M


def Terrahawk(groundset='abcdefghijklmnop'):  # aka the Dillon Mayhew Matroid
    """
    Return the Terrahawk matroid.

    The Terrahawk is a binary matroid that is a sporadic exception in a chain
    theorem for internally `4`-connected binary matroids.

    EXAMPLES::

        sage: M = matroids.catalog.Terrahawk(); M
        Terrahawk: Binary matroid of rank 8 on 16 elements, type (0, 4)
        sage: M.is_valid()
        True

    REFERENCES:

    [CMO2011]_
    """
    A = Matrix(GF(2), [
        [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0]
    ])
    M = BinaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "Terrahawk")
    return M


def ExtendedBinaryGolayCode(groundset='abcdefghijklmnopqrstuvwx'):
    """
    Return the matroid of the extended binary Golay code.

    EXAMPLES::

        sage: M = matroids.catalog.ExtendedBinaryGolayCode(); M
        Extended Binary Golay Code: Binary matroid of rank 12 on 24 elements, type (12, 0)
        sage: C = LinearCode(M.representation())
        sage: C.is_permutation_equivalent(codes.GolayCode(GF(2)))
        True
        sage: M.is_valid()
        True

    .. SEEALSO::

        :class:`GolayCode <sage.coding.golay_code.GolayCode>`
    """
    A = Matrix(GF(2), [
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1],
        [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0],
        [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0],
        [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0],
        [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1],
        [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
         0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
         0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
         0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
         0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    ])
    M = BinaryMatroid(A, groundset)
    M = _rename_and_relabel(M, "Extended Binary Golay Code")
    return M


def CompleteGraphic(n, groundset=None):
    """
    Return the cycle matroid of the complete graph on `n` vertices.

    INPUT:

    - ``n`` -- integer; the number of vertices of the underlying complete
      graph

    OUTPUT: the graphic matroid associated with the `n`-vertex complete graph.
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
    M = _rename_and_relabel(M, f'M(K{n})', groundset)
    return M


# helper function


def _rename_and_relabel(M, name=None, groundset=None):
    """
    Return a renamed and relabeled matroid.

    This is a helper function for easily renaming and relabeling matroids upon
    definition in the context of the database of matroids.

    INPUT:

    - ``M`` -- matroid
    - ``name`` -- string (optional)
    - ``groundset`` -- string (optional)

    OUTPUT: matroid
    """
    if groundset is not None:
        if len(groundset) != len(M.groundset()):
            raise ValueError(
                "the groundset should be of size %s (%s given)" %
                (len(M.groundset()), len(groundset))
            )
        M = M.relabel(dict(zip(M.groundset(), groundset)))

    if name is not None:
        M.rename(name+": " + repr(M))

    return M
