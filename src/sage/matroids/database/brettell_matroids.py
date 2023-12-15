r"""
Documentation for the matroids in the catalog

This module contains implementations of Brettell's interesting matroids,
accessible through :mod:`matroids.catalog. <sage.matroids.catalog>` (type
and hit ``tab`` for a list).

AUTHORS:

- Nick Brettell (2023-02-25, date of associated publication): initial version
                (retrieved online on 2023-12-06)
- Giorgos Mousa (2023-12-08): import to SageMath and add examples

Functions
=========

"""
from sage.matrix.constructor import Matrix
from sage.matroids.circuit_closures_matroid import CircuitClosuresMatroid
from sage.matroids.linear_matroid import TernaryMatroid, QuaternaryMatroid
from sage.matroids.database.oxley_matroids import Uniform
from sage.rings.finite_rings.finite_field_constructor import GF


# 7 elements:


def RelaxedNonFano():
    """
    Return the relaxed NonFano matroid

    An excluded minor for `2`-regular matroids. UPF is `K_2`.

    EXAMPLES::

        sage: M = matroids.catalog.RelaxedNonFano()
        sage: M.is_valid()
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(gf4, [[1, 1, 0, 1], [1, 0, 1, 1], [0, 1, w, 1]])
    RelaxedNonFano = QuaternaryMatroid(reduced_matrix=A)
    RelaxedNonFano.rename("F7=: " + repr(RelaxedNonFano))
    return RelaxedNonFano


def TippedFree3spike():
    """
    Return the tipped free `3`-spike

    Unique 3-connected extension of
    :func:`U36 <sage.matroids.database.oxley_matroids.U36>`. Stabilizer for
    `K_2`.

    EXAMPLES::

        sage: M = matroids.catalog.TippedFree3spike()
        sage: M.has_minor(matroids.Uniform(3,6))
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(gf4, [[1, 1, 1, 1], [1, w + 1, 0, w], [1, 0, w + 1, w]])
    spike = QuaternaryMatroid(
        reduced_matrix=A, groundset=[0, 3, 5, 1, 4, 6, 2]
    )
    spike.rename("Tipped rank-3 free spike: " + repr(spike))
    return spike


# 8 elements:


def AG23minusDY():
    r"""
    Return the matroid `AG23minusDY`

    The matroid obtained from a `AG(2, 3)\setminus e` by a single `\delta-Y`
    exchange on a triangle. An excluded minor for near-regular matroids. UPF
    is `S`.

    EXAMPLES::

        sage: M = matroids.catalog.AG23minusDY()
        sage: M.is_valid()
        True

    """
    A = Matrix(GF(3), [[1, 1, 1, 1], [1, 0, 1, 2], [2, 0, 1, 2], [2, 1, 1, 0]])
    dy = TernaryMatroid(reduced_matrix=A)
    dy.rename("Delta-Y of AG(2,3)\\e: " + repr(dy))
    return dy


def TQ8():
    """
    Return the matroid `TQ8`

    An excluded minor for `2`-regular matroids.
    UPF is `K_2`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.TQ8()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4, [[0, w, 1, 1], [1, 0, w, w + 1], [1, w, 0, w], [1, w + 1, 1, 0]]
    )
    tp8 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 7, 5, 3, 8, 6, 4, 2]
    )
    tp8.rename("TQ8: " + repr(tp8))
    return tp8


def P8p():
    """
    Return the matroid `P8^-`

    `P8^-` is obtained by relaxing one of the disjoint circuit-hyperplanes of
    :func:`P8 <sage.matroids.database.oxley_matroids.P8>`. An excluded minor
    for `2`-regular matroids. UPF is `K_2`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.P8p()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4, [[1, 1, 1, w], [1, w + 1, 1, 0], [1, 0, w, w], [0, 1, 1, 1]]
    )
    p8p = QuaternaryMatroid(
        reduced_matrix=A, groundset=["a", "c", "b", "f", "d", "e", "g", "h"]
    )
    p8p.rename("P8-: " + repr(p8p))
    return p8p


def KP8():
    """
    Return the matroid `KP8`

    An excluded minor for `K_2`-representable matroids.
    UPF is `G`. Self-dual. Uniquely `GF(5)`-representable.
    (An excluded minor for `H_2`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.KP8()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4, [[0, 1, 1, 1], [1, 0, w, w], [1, 1, 1, 1 + w], [1, 1, 1 + w, 0]]
    )
    kp8 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 4, 3, 5, 6, 7, 0, 2]
    )
    kp8.rename("KP8: " + repr(kp8))
    return kp8


def Sp8():
    """
    Return the matroid `Sp8`

    An excluded minor for `G`- and `K_2`-representable matroids.
    UPF is `U_1^{(2)}`. Self-dual.
    (An excluded minor for `H_2`- and `GF(5)`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.Sp8()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4, [[1, 1, w + 1, 0], [1, 1, 0, w + 1], [1, 0, w, w], [0, 1, 1, 1]]
    )
    sp = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 2, 3, 5, 4, 6, 7, 8]
    )
    sp.rename("Sp8: " + repr(sp))
    return sp


def Sp8pp():
    """
    Return the matroid `Sp8=`

    An excluded minor for `G`- and `K_2`-representable matroids.
    UPF is `(GF(2)(a,b),<a,b,a+1,b+1,ab+a+b>)`. Self-dual.
    (An excluded minor for `H_2`- and `GF(5)`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.Sp8pp()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(gf4, [[1, w, 1, 0], [1, 1, 1, 1], [w, 0, 1, w], [0, w, 1, 1]])
    sp = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 5, 6, 7, 2, 3, 4, 8]
    )
    sp.rename("Sp8=: " + repr(sp))
    return sp


def LP8():
    """
    Return the matroid `LP8`

    An excluded minor for `G`- and `K_2`-representable matroids.
    Self-dual. UPF is `W`.
    (Also an excluded minor for `H_2`- and `GF(5)`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.LP8()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4, [[1, 1, 1, 1], [w + 1, w, 0, 1], [1, 0, w + 1, 1], [0, w, w, 1]]
    )
    lp8 = QuaternaryMatroid(
        reduced_matrix=A, groundset=["a", "b", "d", "e", "c", "f", "g", "h"]
    )
    lp8.rename("LP8: " + repr(lp8))
    return lp8


def WQ8():
    r"""
    Return the matroid `WQ8`

    An excluded minor for `G`, `K_2`, `H_4`, and `GF(5)`-representable
    matroids. Self-dual. UPF is `(Z[\zeta,a], <\zeta,a-\zeta>)` where `\zeta`
    is solution to `x^2-x+1 = 0` and `a` is an indeterminate.

    EXAMPLES::

        sage: M = matroids.catalog.WQ8()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4, [[1, 0, 1, w + 1], [1, 1, 1, 1], [w, 1, 1, 0], [0, w, 1, 1]]
    )
    wq8 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[0, 1, 3, 4, 2, 5, 6, 7]
    )
    wq8.rename("WQ8: " + repr(wq8))
    return wq8


# 9 elements:


def BB9():
    """
    Return the matroid `BB9`

    An excluded minor for `K_2`-representable matroids, and a restriction of
    the Betsy Ross matroid. The UPF is `G`. Uniquely `GF(5)`-representable.
    (An excluded minor for `H_2`-representable matroids.)

    EXAMPLES::

        sage: BB = matroids.catalog.BB9()
        sage: BR = matroids.catalog.BetsyRoss()
        sage: for M in BB.extensions(): # long time
        ....:     for N in M.extensions():
        ....:         if N.is_isomorphic(BR):
        ....:             print(True)
        True
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [[1, 0, 1, 1, 1, 1], [0, 1, w, 1, 0, w], [w + 1, 1, w + 1, 1, w, 0]]
    )
    bb9 = QuaternaryMatroid(
        reduced_matrix=A,
        groundset=["i", "b", "d", "j", "h", "f", "c", "a", "k"]
    )
    bb9.rename("BB9: " + repr(bb9))
    return bb9


def TQ9():
    """
    Return the matroid `TQ9`

    An excluded minor for `K_2`-representable matroids, and
    a single-element extension of
    :func:`TQ8 <sage.matroids.database.brettell_matroids.TQ8>`.
    The UPF is `G`. Uniquely `GF(5)`-representable.
    (An excluded minor for `H_2`-representable matroids.)

    EXAMPLES::

        sage: TQ8 = matroids.catalog.TQ8()
        sage: TQ9 = matroids.catalog.TQ9()
        sage: for M in TQ8.extensions():
        ....:     if M.is_isomorphic(TQ9):
        ....:         print(True)
        ....:         break
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4, [[1, 0, w, 1, 1],
              [w + 1, 0, 0, w, 1],
              [1, w, 0, 0, w + 1],
              [1, 1, 1, 1, 0]],
    )
    tq9 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 4, 6, 0, 2, 5, 3, 7, 8]
    )
    tq9.rename("TQ9: " + repr(tq9))
    return tq9


def TQ9p():
    """
    Return the matroid `TQ9^-`

    An excluded minor for `G`- and `K_2`-representable matroids, and
    a single-element extension of
    :func:`TQ8 <sage.matroids.database.brettell_matroids.TQ8>`. UPF is
    `U_1^{(2)}`. (An excluded minor for `H_2`- and
    `GF(5)`-representable matroids.)

    EXAMPLES::

        sage: TQ8 = matroids.catalog.TQ8()
        sage: TQ9p = matroids.catalog.TQ9p()
        sage: for M in TQ8.extensions():
        ....:     if M.is_isomorphic(TQ9p):
        ....:         print(True)
        ....:         break
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 1, 1, w + 1, 1],
            [w + 1, w, w + 1, w + 1, w + 1],
            [w, 0, w, 1, w + 1],
            [0, 1, w + 1, w + 1, 0],
        ],
    )
    tq9p = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 4, 7, 8, 0, 6, 5, 2, 3]
    )
    tq9p.rename("TQ9': " + repr(tq9p))
    return tq9p


def M8591():
    r"""
    Return the matroid `M8591`

    An excluded minor for `K_2`-representable matroids.
    A `Y-\delta` exchange on the unique triad gives
    :func:`A9 <sage.matroids.database.brettell_matroids.A9>`. The UPF is `P_4`.

    EXAMPLES::

        sage: M = matroids.catalog.M8591()
        sage: M.is_valid()
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4, [[1, 1, 0, w, 1],
              [0, 1, 1, w, w + 1],
              [1, 0, w, w, 1],
              [0, 0, 1, 1, 0]]
    )
    m8591 = QuaternaryMatroid(reduced_matrix=A)
    m8591.rename("M8591: " + repr(m8591))
    return m8591


def PP9():
    """
    Return the matroid `PP9`

    An excluded minor for `K_2`-representable matroids. A single-element
    extension of `P8^-`. The UPF is `P_4`. Has a
    :func:`P8p <sage.matroids.database.brettell_matroids.P8p>`-minor (delete
    `z`). Uniquely `GF(5)`-representable. (An excluded minor for
    `H_2`-representable matroids.)

    EXAMPLES::

        sage: P8p = matroids.catalog.P8p()
        sage: PP9 = matroids.catalog.PP9()
        sage: for M in P8p.extensions():
        ....:     if M.is_isomorphic(PP9):
        ....:         print(True)
        ....:         break
        True
        sage: M = PP9.delete('z')
        sage: M.is_isomorphic(P8p)
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4, [[1, 1, 1, w, w],
              [1, 1 + w, 1, 0, w],
              [1, 0, w, w, w],
              [0, 1, 1, 1, 1]]
    )
    pp9 = QuaternaryMatroid(
        reduced_matrix=A,
        groundset=["a", "c", "b", "f", "d", "e", "g", "h", "z"]
    )
    pp9.rename("PP9: " + repr(pp9))
    return pp9


def BB9gDY():
    r"""
    Return the matroid `BB9gDY`

    An excluded minor for `K_2`-representable matroids. The UPF is `G`. In a
    `DY^*`-equivalence class of 4 matroids, one of which can be obtained from
    :func:`BB9 <sage.matroids.database.brettell_matroids.BB9>` by a
    segment-cosegment exchange on `\{a,d,i,j\}`. Uniquely
    `GF(5)`-representable. (An excluded minor for `H_2`-representable
    matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.BB9gDY()
        sage: M.is_valid()
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w, w, w + 1, 1],
            [w, 1, 0, 0],
            [w + 1, 1, 0, 0],
            [1, w, w + 1, w],
            [w, 0, 1, 1],
        ],
    )
    m9573 = QuaternaryMatroid(
        reduced_matrix=A, groundset=["c", "d", "i", "f", "h",
                                     "a", "j", "k", "b"]
    )
    m9573.rename("Segment cosegment exchange on BB9: " + repr(m9573))
    return m9573


def A9():
    """
    Return the matroid `A9`

    An excluded minor for `K_2`-representable matroids.
    The UPF is `P_4`. Uniquely `GF(5)`-representable.
    (An excluded minor for `H_2`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.A9()
        sage: M.is_valid()
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4, [[w + 1, 1, w, w, w, w],
              [0, 1, 1, w + 1, 0, w],
              [w, 0, 1, w + 1, w, 1]]
    )
    a9 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[6, 5, 4, 1, 2, 3, 7, 8, 0]
    )
    a9.rename("A9: " + repr(a9))
    return a9


def FN9():
    """
    Return the matroid `FN9`

    An excluded minor for `G`- and `K_2`-representable matroids.
    In a `DY^*`-equivalence class of `10` matroids. UPF is `U_1^{(2)}`.
    (An excluded minor for `H_2`- and `GF(5)`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.FN9()
        sage: M.is_valid()
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w + 1, w, w + 1, w, 1, 0],
            [1, w + 1, 0, 1, w + 1, 1],
            [w + 1, w + 1, w, w + 1, 1, 1],
        ],
    )
    m3209 = QuaternaryMatroid(
        reduced_matrix=A, groundset=["b0", "a", "y", "z", "x",
                                     "c0", "b", "c", "a0"]
    )
    m3209.rename("FN9: " + repr(m3209))
    return m3209


def FX9():
    """
    Return the matroid `FX9`

    An excluded minor for `G`- and `K_2`-representable matroids.
    UPF is `(Q(a,b), <-1,a,b,a-1,b-1,a-b,a+b,a+b-2,a+b-2ab>)`.
    (An excluded minor for `H_2`- and `GF(5)`-representable matroids.)

    EXAMPLES::

        sage: M = matroids.catalog.FX9()
        sage: M.is_valid()
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, w + 1, 0, w, 1],
            [1, w, 1, w + 1, 1],
            [w + 1, w + 1, w, w + 1, w + 1],
            [w, w, w + 1, w + 1, 1],
        ],
    )
    m48806 = QuaternaryMatroid(reduced_matrix=A)
    m48806.rename("FX9: " + repr(m48806))
    return m48806


def KR9():
    """
    Return the matroid `KR9`

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids.) In a `DY`-equivalence class of `4`
    matroids. Has a
    :func:`KP8 <sage.matroids.database.brettell_matroids.KP8>`-minor (delete
    `8`). UPF is `GF(4)`.

    EXAMPLES::

        sage: KR9 = matroids.catalog.KR9()
        sage: KP8 = matroids.catalog.KP8()
        sage: KP8.is_isomorphic(KR9.delete(8))
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w + 1, w, w + 1, w, 1],
            [0, 1, w + 1, w, 1],
            [w, w + 1, 1, 1, 1],
            [w + 1, w + 1, w, w + 1, 0],
        ],
    )
    kr9 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[2, 4, 0, 6, 1, 5, 3, 7, 8]
    )
    kr9.rename("KR9: " + repr(kr9))
    return kr9


def KQ9():
    """
    Return the matroid `KQ9`

    An excluded minor for `G`-representable matroids (and
    `GF(5)`-representable matroids.) Has a
    :func:`TQ8 <sage.matroids.database.brettell_matroids.TQ8>`-minor`
    (delete `6`) and a
    :func:`KP8 <sage.matroids.database.brettell_matroids.KP8>`-minor
    (delete `8`). UPF is `GF(4)`.

    EXAMPLES::

        sage: KQ9 = matroids.catalog.KQ9()
        sage: TQ8 = matroids.catalog.TQ8()
        sage: TQ8.is_isomorphic(KQ9.delete(6))
        True
        sage: KP8 = matroids.catalog.KP8()
        sage: KP8.is_isomorphic(KQ9.delete(8))
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w + 1, w, w + 1, 1, w + 1],
            [1, 0, w + 1, w + 1, 1],
            [0, 1, w, w + 1, 1],
            [1, 1, w + 1, 0, w + 1],
        ],
    )
    kq9 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[5, 0, 4, 3, 2, 6, 8, 7, 1]
    )
    kq9.rename("KQ9: " + repr(kq9))
    return kq9


# 10 elements:


def UG10():
    """
    Return the matroid `UG10`

    An excluded minor for `K_2`- and `P_4`-representable matroids. Self-dual.
    An excluded minor for `H_3`- and `H_2`-representable matroids.
    Uniquely `GF(5)`-representable.
    Although not `P_4`-representable, it is `O`-representable,
    and hence is representable over all fields of size at least four.

    EXAMPLES::

        sage: M = matroids.catalog.UG10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 0, 1, w, w + 1],
            [1, 1, 1, 1, 1],
            [1, 0, w, 1, w + 1],
            [1, w + 1, w, 1, 0],
            [1, 1, 1, 0, 0],
        ],
    )
    ug10 = QuaternaryMatroid(reduced_matrix=A)
    ug10.rename("UG10: " + repr(ug10))
    return ug10


def FF10():
    """
    Return the matroid `FF10`

    An excluded minor for `K_2`-representable matroids.
    UPF is `P_4`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FF10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, 0, 1, 1, 1],
            [1, 0, 0, 1, 1],
            [1 + w, 1, 0, 0, 1],
            [1, 1, w, 0, 1],
            [1 + w, 1 + w, w, w, 1 + w],
        ],
    )
    ff10 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    )
    ff10.rename("FF10: " + repr(ff10))
    return ff10


def GP10():
    """
    Return the matroid `GP10`

    An excluded minor for `K_2`-representable matroids.
    UPF is `G`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.GP10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w + 1, w, 0, 1, 1],
            [w, w, 0, 0, 1],
            [0, 0, w + 1, w, 1],
            [1, 0, w, 0, 1],
            [1, 1, 1, 1, 0],
        ],
    )
    gp10 = QuaternaryMatroid(reduced_matrix=A)
    gp10.rename("GP10: " + repr(gp10))
    return gp10


def FZ10():
    """
    Return the matroid `FZ10`

    An excluded minor for `K_2`- and `G`-representable matroids
    (and `H_2`- and `GF(5)`-representable matroids).
    UPF is `W`. Not self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FZ10()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, 0, 1, w, 1],
            [1, 0, w, w + 1, w + 1],
            [1, 1, w, w, w],
            [0, 1, 0, w, w],
            [1, 1, w + 1, 1, w],
        ],
    )
    fz10 = QuaternaryMatroid(reduced_matrix=A)
    fz10.rename("FZ10: " + repr(fz10))
    return fz10


def UQ10():
    """
    Return the matroid `UQ10`

    An excluded minor for `K_2`- and `G`-representable matroids
    (and `H_2`- and `GF(5)`-representable matroids).
    Self-dual. UPF is `(Q(a,b), <-1,a,b,a-1,b-1,a-b,a+b,a+1,ab+b-1,ab-b+1>)`.

    EXAMPLES::

        sage: M = matroids.catalog.UQ10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, 0, 1, w, 1],
            [1, 1, 1, 1, 0],
            [1, w + 1, 0, 1, 0],
            [w + 1, w, 0, 0, 1],
            [1, 1, 1, w + 1, 1],
        ],
    )
    uq10 = QuaternaryMatroid(reduced_matrix=A)
    uq10.rename("UQ10: " + repr(uq10))
    return uq10


def FP10():
    """
    Return the matroid `FP10`

    An excluded minor for `K_2`- and `G`-representable matroids
    (and `H_2`- and `GF(5)`-representable matroids).
    UPF is `U_1^{(2)}`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FP10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 1, 1, 1, 1],
            [0, 1, 1, 1, 0],
            [w, w, w + 1, 0, 0],
            [0, 1, 1, 0, w],
            [1, w, 1, w, 1],
        ],
    )
    fp10 = QuaternaryMatroid(reduced_matrix=A)
    fp10.rename("FP10: " + repr(fp10))
    return fp10


def TQ10():
    """
    Return the matroid `TQ10`

    An excluded minor for `K_2`-representable matroids. UPF is `G`. Self-dual.
    Has :func:`TQ8 <sage.matroids.database.brettell_matroids.TQ8>` as a minor
    (delete 'd' and contract 'c').

    EXAMPLES::

        sage: M = matroids.catalog.TQ10()
        sage: M.is_isomorphic(M.dual())
        True
        sage: N = M.delete('d').contract('c')
        sage: N.is_isomorphic(matroids.catalog.TQ8())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, w, 0, w + 1, 1],
            [w + 1, w + 1, w + 1, w, 1],
            [w + 1, w, 1, 0, 1],
            [1, 1, 0, w + 1, 0],
            [w + 1, 0, w, w + 1, 0],
        ],
    )
    tp10 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 6, 8, "c", 3, 7, "d", 2, 5, 4]
    )
    tp10.rename("TQ10: " + repr(tp10))
    return tp10


def FY10():
    """
    Return the matroid `FY10`

    An excluded minor for `P_4`-representable matroids. UPF is `G`. Not
    self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FY10()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 1, 1, 0, 0],
            [0, 1, w + 1, 1, 0],
            [1, 1, 1, w + 1, 1],
            [0, 0, w, 1, 1],
            [1, w + 1, 1, 1, w],
        ],
    )
    fy10 = QuaternaryMatroid(reduced_matrix=A)
    fy10.rename("FY10: " + repr(fy10))
    return fy10


def PP10():
    """
    Return the matroid `PP10`

    An excluded minor for `P_4`-representable matroids. UPF is `U_1^{(2)}`.
    Has a :func:`TQ8 <sage.matroids.database.brettell_matroids.TQ8>`-minor
    (e.g. delete 'a' and contract 'e') and a
    :func:`PP9 <sage.matroids.database.brettell_matroids.PP9>` (and hence
    :func:`P8p <sage.matroids.database.brettell_matroids.P8p>`) minor
    (contract 'x').

    EXAMPLES::

        sage: PP10 = matroids.catalog.PP10()
        sage: M = PP10.delete('a').contract('e')
        sage: M.is_isomorphic(matroids.catalog.TQ8())
        True
        sage: M = PP10.contract('x')
        sage: M.is_isomorphic(matroids.catalog.PP9())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w + 1, 0, w + 1, 0, w],
            [w, w, 1, w, 1],
            [w + 1, w + 1, 0, w + 1, 1],
            [1, 0, 1, w + 1, 1],
            [w, 1, w + 1, w, w + 1],
        ],
    )
    pp10 = QuaternaryMatroid(
        reduced_matrix=A, groundset=["z", "f", "c", "g", "e",
                                     "b", "a", "h", "d", "x"]
    )
    pp10.rename("PP10: " + repr(pp10))
    return pp10


def FU10():
    """
    Return the matroid `FU10`

    An excluded minor for `P_4`-representable matroids. UPF is `G`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FU10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, 1, w, 1, 1],
            [w, w + 1, 1, 1, w],
            [1, 1, 0, w + 1, 0],
            [1, 1, 0, w, w + 1],
            [w + 1, w, w + 1, w, 0],
        ],
    )
    fu10 = QuaternaryMatroid(reduced_matrix=A)
    fu10.rename("FU10: " + repr(fu10))
    return fu10


def D10():
    """
    Return the matroid `D10`

    An excluded minor for `P_4`-representable matroids.
    UPF is `G`. Has a
    :func:`TQ8 <sage.matroids.database.brettell_matroids.TQ8>`-minor.
    In a `DY^*`-equivalence class of `13` matroids.

    EXAMPLES::

        sage: M = matroids.catalog.D10()
        sage: M.has_minor(matroids.catalog.TQ8())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w, 1, w, 1, w + 1, w],
            [w, 0, w + 1, w + 1, w, w],
            [w + 1, 0, 0, w + 1, w + 1, w + 1],
            [w + 1, 1, 0, 1, w, 0],
        ],
    )
    d10 = QuaternaryMatroid(reduced_matrix=A)
    d10.rename("D10: " + repr(d10))
    return d10


def UK10():
    """
    Return the matroid `UK10`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Not self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.UK10()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, w, w + 1, w, w + 1],
            [1, w, w, w + 1, w],
            [1, w + 1, w + 1, 0, w],
            [1, 1, w, 0, w],
            [w + 1, 0, 0, 1, w],
        ],
    )
    uk10 = QuaternaryMatroid(reduced_matrix=A)
    uk10.rename("UK10: " + repr(uk10))
    return uk10


def PK10():
    """
    Return the matroid `PK10`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Not self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.PK10()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 1, w, 0, w + 1],
            [0, 1, 1, 1, 1],
            [w + 1, w, w, w + 1, w],
            [0, 1, w, w + 1, 1],
            [w + 1, w + 1, 0, 0, 1],
        ],
    )
    pk10 = QuaternaryMatroid(reduced_matrix=A)
    pk10.rename("PK10: " + repr(pk10))
    return pk10


def GK10():
    """
    Return the matroid `GK10`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Not self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.GK10()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 1, 0, 1, 1],
            [1, 0, 0, w, w],
            [w, w + 1, 1, 1, 0],
            [w, w + 1, w + 1, w + 1, w + 1],
            [1, w, w, w + 1, 1],
        ],
    )
    gk10 = QuaternaryMatroid(reduced_matrix=A)
    gk10.rename("GK10: " + repr(gk10))
    return gk10


def FT10():
    """
    Return the matroid `FT10`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.FT10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w, 0, w, w + 1, w + 1],
            [0, 1, w + 1, w, 1],
            [w, 1, 0, w + 1, w + 1],
            [w, 1, 0, 0, w],
            [0, 1, w, 0, 1],
        ],
    )
    ft10 = QuaternaryMatroid(reduced_matrix=A)
    ft10.rename("FT10: " + repr(ft10))
    return ft10


def TK10():
    """
    Return the matroid `TK10`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.TK10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 1, 0, 0, w],
            [0, w, 0, w + 1, 1],
            [w + 1, 0, w + 1, w, w + 1],
            [w, 0, 1, 0, w],
            [0, w, 1, w, 0],
        ],
    )
    tk10 = QuaternaryMatroid(reduced_matrix=A)
    tk10.rename("TK10: " + repr(tk10))
    return tk10


def KT10():
    """
    Return the matroid `KT10`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.KT10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, 1, 1, w, 0],
            [0, 1, 0, w + 1, 1],
            [w + 1, 0, 1, w + 1, w + 1],
            [w, w + 1, w + 1, 1, w],
            [w + 1, w + 1, w + 1, 1, 1],
        ],
    )
    kt10 = QuaternaryMatroid(reduced_matrix=A)
    kt10.rename("KT10: " + repr(kt10))
    return kt10


def TU10():
    """
    Return the matroid `TU10`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.TU10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, w + 1, 1, 0, w],
            [w, 0, 1, w + 1, w + 1],
            [w, w, w + 1, w + 1, 1],
            [w + 1, 0, w + 1, 1, w + 1],
            [w + 1, w + 1, 1, w, 1],
        ],
    )
    tu10 = QuaternaryMatroid(reduced_matrix=A)
    tu10.rename("TU10: " + repr(tu10))
    return tu10


def UT10():
    """
    Return the matroid `UT10`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Self-dual. UPF is `I`.

    EXAMPLES::

        sage: M = matroids.catalog.UT10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w, w + 1, 0, w + 1, 0],
            [w + 1, w + 1, 1, w, 1],
            [1, 0, 1, w + 1, w + 1],
            [w + 1, w + 1, 1, 1, w],
            [1, w + 1, 1, w, w + 1],
        ],
    )
    ut10 = QuaternaryMatroid(reduced_matrix=A)
    ut10.rename("UT10: " + repr(ut10))
    return ut10


def FK10():
    """
    Return the matroid `FK10`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.FK10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, 0, 1, 1, w],
            [w, 1, w, 1, 1],
            [1, 1, 0, 0, w + 1],
            [w + 1, 0, w + 1, 0, 1],
            [w + 1, w + 1, 1, 1, w + 1],
        ],
    )
    fk10 = QuaternaryMatroid(reduced_matrix=A)
    fk10.rename("FK10: " + repr(fk10))
    return fk10


def KF10():
    """
    Return the matroid `KF10`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.KF10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w + 1, w + 1, 1, w, 1],
            [0, w + 1, w, 1, 0],
            [0, w, 0, 1, w + 1],
            [w + 1, 0, w, w + 1, 1],
            [w + 1, 1, 0, w + 1, 0],
        ],
    )
    kf10 = QuaternaryMatroid(reduced_matrix=A)
    kf10.rename("KF10: " + repr(kf10))
    return kf10


# 11 elements:


def FA11():
    """
    Return the matroid `FA11`

    An excluded minor for `P_4`-representable matroids. UPF is `PT`. In a
    `DY^*`-equivalence class of `6` matroids. Has an
    :func:`FF10 <sage.matroids.database.brettell_matroids.FF10>`-minor (delete
    `10`).

    EXAMPLES::

        sage: FA11 = matroids.catalog.FA11()
        sage: FF10 = matroids.catalog.FF10()
        sage: FF10.is_isomorphic(FA11.delete(10))
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w, 0, w, w, 1, 0],
            [1, w, 0, w, w, 1],
            [0, w, 1, w + 1, 0, w],
            [0, w, 0, w + 1, 0, 1],
            [w + 1, w + 1, w + 1, 0, w, 0],
        ],
    )
    fa11 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[1, 3, 4, 2, 8, 7, 9, 0, 5, 10, 6]
    )
    fa11.rename("FA11: " + repr(fa11))
    return fa11


# 12 elements:


def FR12():
    """
    Return the matroid `FR12`

    An excluded minor for `K_2`-representable matroids.
    UPF is `P_4`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FR12()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, 1, 1, 1, 0, 1],
            [1, 1, 1, 0, 1, 1],
            [1, 1, 0, 0, 0, 1],
            [1, 0, 0, 0, w + 1, 1],
            [0, 1, 0, w + 1, 0, 1],
            [1, 1, 1, 1, 1, w + 1],
        ],
    )
    fr12 = QuaternaryMatroid(reduced_matrix=A)
    fr12.rename("FR12: " + repr(fr12))
    return fr12


def GP12():
    """
    Return the matroid `GP12`

    An excluded minor for `K_2`-representable matroids.
    UPF is `G`. Not self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.GP12()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1 + w, 1 + w, w, 1 + w, 0, 0],
            [1, 0, w, 1 + w, 0, 0],
            [1, 0, 1, 0, 1 + w, 1 + w],
            [1, 1, 0, 1, w, w],
            [1, 0, 1, 1, 0, 1 + w],
            [1, 1, 1, 1, 1, 1 + w],
        ],
    )
    gp12 = QuaternaryMatroid(reduced_matrix=A)
    gp12.rename("GP12: " + repr(gp12))
    return gp12


def FQ12():
    """
    Return the matroid `FQ12`

    An excluded minor for `P_4`-representable matroids. UPF is `PT`. Has` a
    :func:`PP9 <sage.matroids.database.brettell_matroids.PP9>`-minor (contract
    `4` and `7`, delete `6`) and
    :func:`FF10 <sage.matroids.database.brettell_matroids.FF10>`-minor
    (contract 'c' and delete 'd').

    EXAMPLES::

        sage: FQ12 = matroids.catalog.FQ12()
        sage: PP9 = matroids.catalog.PP9()
        sage: PP9.is_isomorphic(FQ12.contract([4,7]).delete(6))
        True
        sage: FF10 = matroids.catalog.FF10()
        sage: FF10.is_isomorphic(FQ12.contract('c').delete('d'))
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, 0, w, w, 1, 0],
            [0, 0, w + 1, w + 1, 1, 1],
            [1, 1, w, 1, 1, 1],
            [w, 0, 1, 1, 0, 0],
            [w, w, w + 1, w + 1, 1, 1],
            [0, 1, 1, w, w + 1, 1],
        ],
    )
    fq12 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[7, 4, 5, 9, 2, 1, 0, 6, "d", "c", 8, 3]
    )
    fq12.rename("FQ12: " + repr(fq12))
    return fq12


def FF12():
    """
    Return the matroid `FF12`

    An excluded minor for `P_4`-representable matroids. Self-dual. UPF is
    `(Q(a,b),<-1,a,b,a-2,a-1,a+1,b-1,ab-a+b,ab-a-b,ab-a-2b>)`. Has an
    :func:`FF10 <sage.matroids.database.brettell_matroids.FF10>`-minor
    (contract 'c' and delete 'd').

    EXAMPLES::

        sage: M = matroids.catalog.FF12()
        sage: M.is_isomorphic(M.dual())
        True
        sage: FF10 = matroids.catalog.FF10()
        sage: FF10.is_isomorphic(M.contract('c').delete('d'))
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 1, 0, 1, 0, 0],
            [1, 1, w, 1, 0, 1],
            [w, w + 1, 1, 0, 1, 1],
            [1, 1, w, 1, w + 1, 0],
            [1, 1, 0, 0, w + 1, 0],
            [1, w + 1, 1, 0, w + 1, w],
        ],
    )
    ff12 = QuaternaryMatroid(
        reduced_matrix=A, groundset=[0, 4, "c", 3, 5, "d", 8, 9, 2, 7, 1, 6]
    )
    ff12.rename("FF12: " + repr(ff12))
    return ff12


def FZ12():
    """
    Return the matroid `FZ12`

    An excluded minor for `K_2`- and `G`-representable matroids
    (and `H_2`- and `GF(5)`-representable matroids).
    UPF is `W`. Not self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FZ12()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, w + 1, 0, w + 1, 0, w + 1],
            [w + 1, w + 1, w, 0, w, 1],
            [w, 1, 0, w, 0, w],
            [w + 1, 1, 1, w + 1, 1, w],
            [w, w, 1, 0, 0, 0],
            [w + 1, 1, 0, w, 1, w],
        ],
    )
    fz12 = QuaternaryMatroid(reduced_matrix=A)
    fz12.rename("FZ12: " + repr(fz12))
    return fz12


def UQ12():
    """
    Return the matroid `UQ12`

    An excluded minor for `K_2` and `G`-representable matroids
    (and `H2` and `GF(5)`-representable matroids).
    UPF is `P_{pappus}`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.UQ12()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 0, 0, w + 1, 1, 0],
            [w + 1, w, w + 1, 1, 1, 1],
            [w + 1, w, w + 1, w + 1, w, 1],
            [1, 0, 0, w, 1, w + 1],
            [1, 0, 1, w, 1, w],
            [1, 1, 1, w, 1, w + 1],
        ],
    )
    uq12 = QuaternaryMatroid(reduced_matrix=A)
    uq12.rename("UQ12: " + repr(uq12))
    return uq12


def FP12():
    """
    Return the matroid `FP12`

    An excluded minor for `K_2`- and `G`-representable matroids
    (and `H_2`- and `GF(5)`-representable matroids).
    UPF is `W`. Self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.FP12()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, w + 1, 1, 0, 1, w],
            [0, w + 1, 1, 0, w + 1, 0],
            [w + 1, 1, w, 0, w + 1, w],
            [w, 1, w, 1, w, w + 1],
            [w + 1, 0, w + 1, w, w + 1, 0],
            [w, 0, w + 1, w + 1, w, 0],
        ],
    )
    fp12 = QuaternaryMatroid(reduced_matrix=A)
    fp12.rename("FP12: " + repr(fp12))
    return fp12


def FS12():
    """
    Return the matroid `FS12`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Rank `5`. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.FS12()
        sage: M.rank()
        5
    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 1, 1, 1, 1, 0, 1],
            [0, 1, w + 1, w, w + 1, 1, 1],
            [w, w + 1, 0, 0, w + 1, 0, 0],
            [0, 0, w + 1, 0, 0, 1, w],
            [1, 0, w, w, w, 1, 0],
        ],
    )
    fs12 = QuaternaryMatroid(reduced_matrix=A)
    fs12.rename("FS12: " + repr(fs12))
    return fs12


def UK12():
    """
    Return the matroid `UK12`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Self-dual. UPF is I.

    EXAMPLES::

        sage: M = matroids.catalog.UK12()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, w + 1, 1, 0, 0, w],
            [w, 1, w, w + 1, w, w + 1],
            [1, w + 1, w + 1, 0, 0, 0],
            [0, w, 1, 1, w + 1, 1],
            [0, w, 1, 1, 0, w + 1],
            [w, w, w + 1, w, w, 1],
        ],
    )
    uk12 = QuaternaryMatroid(reduced_matrix=A)
    uk12.rename("UK12: " + repr(uk12))
    return uk12


def UA12():
    """
    Return the matroid `UA12`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Not self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.UA12()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, w + 1, w, w + 1, w + 1, 0],
            [1, 1, 1, w + 1, w, 1],
            [1, w + 1, 0, 1, w + 1, 1],
            [0, 0, 0, w + 1, w, 1],
            [0, w + 1, 0, w, w + 1, 1],
            [1, w, w + 1, w + 1, w, 1],
        ],
    )
    ua12 = QuaternaryMatroid(reduced_matrix=A)
    ua12.rename("UA12: " + repr(ua12))
    return ua12


def AK12():
    """
    Return the matroid `AK12`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids).
    Not self-dual. UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.AK12()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, w, 0, 0, w + 1, w],
            [w + 1, w, 0, w, w + 1, w + 1],
            [0, w, w + 1, 0, w, 0],
            [1, 0, w, 0, 0, w + 1],
            [0, 1, 0, 1, 1, 0],
            [1, 0, 1, 1, 0, 1],
        ],
    )
    ak12 = QuaternaryMatroid(reduced_matrix=A)
    ak12.rename("AK12: " + repr(ak12))
    return ak12


def FK12():
    """
    Return the matroid `FK12`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids). Self-dual.
    UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.UT10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [w + 1, w, w, 0, w, w + 1],
            [w, 1, 1, w + 1, w + 1, w],
            [1, w + 1, w, 0, w + 1, 0],
            [w, 1, 1, w, w + 1, w + 1],
            [w + 1, w + 1, w, 0, w, 0],
            [1, w, w, w, 1, 1],
        ],
    )
    fk12 = QuaternaryMatroid(reduced_matrix=A)
    fk12.rename("FK12: " + repr(fk12))
    return fk12


def KB12():
    """
    Return the matroid `KB12`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids). Self-dual.
    UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.UT10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 0, w, 0, 0, 1],
            [1, w, 0, w, 1, w + 1],
            [1, 1, w + 1, 0, 1, w + 1],
            [w, w, 1, w, 0, 0],
            [1, 1, 1, 0, 0, 0],
            [w, w, 1, w, 1, 1],
        ],
    )
    kb12 = QuaternaryMatroid(reduced_matrix=A)
    kb12.rename("KB12: " + repr(kb12))
    return kb12


def AF12():
    """
    Return the matroid `AF12`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids). Self-dual.
    UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.UT10()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [1, 0, 0, 0, 1, 1],
            [0, 1, 0, w, w, w + 1],
            [0, 1, 0, 0, w + 1, w + 1],
            [w, 1, 1, 0, 0, 1],
            [0, 1, w + 1, w, w, 1],
            [1, 1, 1, w + 1, w + 1, w],
        ],
    )
    af12 = QuaternaryMatroid(reduced_matrix=A)
    af12.rename("AF12: " + repr(af12))
    return af12


# utility function
def complement(groundset, subset):
    return list(set(groundset).difference(subset))


def NestOfTwistedCubes():
    r"""
    Return the NestOfTwistedCubes matroid

    A matroid with no `U(2,4)`-detachable pairs (only `\{e_i,f_i\}` pairs are
    detachable).

    EXAMPLES::

        sage: M = matroids.catalog.NestOfTwistedCubes()
        sage: M.is_valid() # long time
        True

    """
    gs = ["e1", "e2", "e3", "e4", "e5", "e6",
          "f1", "f2", "f3", "f4", "f5", "f6"]
    ntc = CircuitClosuresMatroid(
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
    return ntc


# 13 elements:


def XY13():
    """
    Return the matroid `XY13`

    An excluded minor for `G`-representable matroids
    (and `GF(5)`-representable matroids). UPF is `GF(4)`.

    EXAMPLES::

        sage: M = matroids.catalog.XY13()
        sage: M.is_3connected()
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
        [
            [0, 0, 1, 1, 0, 1, 1],
            [w, 1, w, w + 1, 1, 1, w + 1],
            [0, 0, w + 1, 1, 1, w, 1],
            [0, w, 1, 1, w + 1, 1, 1],
            [w, w + 1, w, w + 1, 1, 1, w],
            [1, 0, 0, 1, 0, w, 0],
        ],
    )
    xy13 = QuaternaryMatroid(reduced_matrix=A)
    xy13.rename("XY13: " + repr(xy13))
    return xy13


# 14 elements:


def N3():
    """
    Return the matroid `N3`

    An excluded minor for dyadic matroids (and `GF(5)`-representable matroids).
    UPF is `GF(3)`. `4`- (but not `5`-) connected. Self-dual.

    EXAMPLES::

        sage: N3 = matroids.catalog.N3()
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
    n3 = TernaryMatroid(reduced_matrix=A)
    n3.rename("N3: " + repr(n3))
    return n3


def N3pp():
    """
    Return the matroid `N3pp`

    An excluded minor for `K_2`-representable matroids. Self-dual.
    Obtained by relaxing the two complementary circuit-hyperplanes of
    :func:`N4 <sage.matroids.database.brettell_matroids.N4>`. Not
    `P_4`-representable, but `O`-representable, and hence representable
    over all fields of size at least four.

    EXAMPLES::

        sage: M = matroids.catalog.N3pp()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
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
    sq14 = QuaternaryMatroid(reduced_matrix=A)
    sq14.rename("N3=: " + repr(sq14))
    return sq14


def UP14():
    """
    Return the matroid `UP14`

    An excluded minor for `K_2`-representable matroids.
    Has disjoint circuit-hyperplanes.
    UPF is `W`. Not self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.UP14()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
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
    up14 = QuaternaryMatroid(reduced_matrix=A)
    up14.rename("UP14: " + repr(up14))
    return up14


def VP14():
    """
    Return the matroid `VP14`

    An excluded minor for `K_2`-representable matroids.
    Has disjoint circuit-hyperplanes.
    UPF is `W`. Not self-dual.

    EXAMPLES::

        sage: M = matroids.catalog.VP14()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
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
    vp14 = QuaternaryMatroid(reduced_matrix=A)
    vp14.rename("VP14: " + repr(vp14))
    return vp14


def FV14():
    """
    Return the matroid `FV14`

    An excluded minor for `P_4`-representable matroids.
    Not self-dual. UPF is `PT`.

    EXAMPLES::

        sage: M = matroids.catalog.FV14()
        sage: M.is_isomorphic(M.dual())
        False

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
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
    fv14 = QuaternaryMatroid(reduced_matrix=A)
    fv14.rename("FV14: " + repr(fv14))
    return fv14


def OW14():
    """
    Return the matroid `OW14`

    An excluded minor for `P_4`-representable matroids.
    Self-dual. UPF is `Orthrus`.

    EXAMPLES::

        sage: M = matroids.catalog.OW14()
        sage: M.is_isomorphic(M.dual())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
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
    ow14 = QuaternaryMatroid(reduced_matrix=A)
    ow14.rename("OW14: " + repr(ow14))
    return ow14


def FM14():
    """
    Return the matroid `FM14`

    An excluded minor for `P_4`-representable matroids.
    Self-dual. UPF is `PT`.

    EXAMPLES::

        sage: M = matroids.catalog.FM14()
        sage: M.is_isomorphic(M.dual())
        True

    """

    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
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
    fm14 = QuaternaryMatroid(reduced_matrix=A)
    fm14.rename("FM14: " + repr(fm14))
    return fm14


# 15 elements:


def FA15():
    """
    Return the matroid `FA15`

    An excluded minor for `O`-representable matroids. UPF is `PT`.
    In a `DY^*`-equivalence class of `6` matroids. Has an
    :func:`SQ14 <sage.matroids.database.brettell_matroids.N3pp>`-minor.

    EXAMPLES::

        sage: M = matroids.catalog.FA15()
        sage: M.has_minor(matroids.catalog.N3pp())
        True

    """
    gf4 = GF(4, "w")
    w = gf4("w")
    A = Matrix(
        gf4,
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
    fa15 = QuaternaryMatroid(reduced_matrix=A)
    fa15.rename("FA15: " + repr(fa15))
    return fa15


# 16 elements:


def N4():
    """
    Return the matroid `N4`

    An excluded minor for dyadic matroids (and `GF(5)`-representable matroids).
    UPF is `GF(3)`. `4`- (but not `5`-) connected. Self-dual.

    EXAMPLES::

        sage: N4 = matroids.catalog.N4()
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
    n4 = TernaryMatroid(reduced_matrix=A)
    n4.rename("N4: " + repr(n4))
    return n4


# 2r elements:


def FreeSpike(r, t=True):
    r"""
    Return the rank-`r` free spike

    For the tipless free spike (`t` = ``False``), when `r = 3`, it is
    isomorphic to :func:`U36 <sage.matroids.database.oxley_matroids.U36>`;
    when `r = 4`, it is the unique tightening of the
    :func:`Vamos matroid <sage.matroids.database.oxley_matroids.Vamos>`.

    EXAMPLES::

        sage: M = matroids.FreeSpike(3, False)
        sage: M.is_isomorphic(matroids.Uniform(3, 6))
        True
        sage: M = matroids.FreeSpike(8)
        sage: M.is_3connected()
        True

    """
    if t:  # tipped free spike
        if r == 3:
            return TippedFree3spike()
        elif r < 3:
            raise ValueError("Tipped free spike must have rank at least 3.")

        E = range(2 * r + 1)
        tris = [[0, 2 * i + 1, 2 * i + 2] for i in range(r)]
        planes = [
            [2 * i + 1, 2 * i + 2, 2 * j + 1, 2 * j + 2]
            for i in range(r)
            for j in range(i + 1, r)
        ]
        CC = {2: tris, 3: planes, r: [E]}
        spike = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
        spike.rename("Tipped rank-" + str(r) + " free spike: " + repr(spike))
        return spike
    else:  # tipless free spike
        if r == 3:
            return Uniform(3, 6)
        elif r < 3:
            raise ValueError("Tipless free spike must have rank at least 3.")

        E = range(1, 2 * r + 1)
        circs = [
            [2 * i + 1, 2 * i + 2, 2 * j + 1, 2 * j + 2]
            for i in range(r)
            for j in range(i + 1, r)
        ]
        CC = {3: circs, r: [E]}
        spike = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
        spike.rename(
            "(Tipless) rank-" + str(r) + " free spike: " + repr(spike)
        )
        return spike
