r"""
Collections of matroids

This module contains functions that access the collections of matroids in the
database. Each of these functions returns an iterator over the nonparametrized
matroids from the corresponding collection. These functions can be viewed by
typing ``matroids.`` + :kbd:`Tab`.

AUTHORS:

- Giorgos Mousa (2023-12-08): initial version
"""

# ****************************************************************************
#       Copyright (C) 2023 Giorgos Mousa <gmousa@proton.me>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


def OxleyMatroids():
    """
    Return an iterator over Oxley's matroid collection.

    EXAMPLES::

        sage: OM = list(matroids.OxleyMatroids()); len(OM)
        44
        sage: import random
        sage: M = random.choice(OM)
        sage: M.is_valid()
        True

    .. SEEALSO::

        :mod:`Matroid catalog <sage.matroids.matroids_catalog>`, under
        ``Oxley's matroid collection``.

    REFERENCES:

    These matroids are the nonparametrized matroids that appear in the
    Appendix ``Some Interesting Matroids`` in [Oxl2011]_ (p. 639-64).
    """
    from sage.matroids.database_matroids import (
        U24, U25, U35, K4, Whirl3, Q6, P6, U36, R6,
        Fano, FanoDual, NonFano, NonFanoDual, O7, P7,
        AG32, AG32prime, R8, F8, Q8, L8, S8, Vamos, T8, J, P8, P8pp,
        Wheel4, Whirl4,
        K33dual, K33, AG23, TernaryDowling3, R9, Pappus, NonPappus,
        K5, K5dual, R10, NonDesargues,
        R12, ExtendedTernaryGolayCode, T12,
        PG23
    )

    lst = [U24,  # 4
           U25, U35,  # 5
           K4, Whirl3, Q6, P6, U36, R6,  # 6
           Fano, FanoDual, NonFano, NonFanoDual, O7, P7,  # 7
           AG32, AG32prime,
           R8, F8, Q8, L8, S8,
           Vamos, T8, J, P8, P8pp,
           Wheel4, Whirl4,  # 8
           K33dual, K33, AG23, TernaryDowling3, R9, Pappus, NonPappus,  # 9
           K5, K5dual, R10, NonDesargues,  # 10
           R12, ExtendedTernaryGolayCode, T12,  # 12
           PG23]  # 13
    for M in lst:
        yield M()


def BrettellMatroids():
    """
    Return an iterator over Brettell's matroid collection.

    EXAMPLES::

        sage: BM = list(matroids.BrettellMatroids()); len(BM)
        68
        sage: import random
        sage: M = random.choice(BM)
        sage: M.is_valid()
        True

    .. SEEALSO::

        :mod:`Matroid catalog <sage.matroids.matroids_catalog>`, under
        ``Brettell's matroid collection``.
    """
    from sage.matroids.database_matroids import (
        RelaxedNonFano, TippedFree3spike,
        AG23minusDY, TQ8, P8p, KP8, Sp8, Sp8pp, LP8, WQ8,
        BB9, TQ9, TQ9p, M8591, PP9, BB9gDY, A9, FN9, FX9, KR9, KQ9,
        UG10, FF10, GP10, FZ10, UQ10, FP10, TQ10, FY10, PP10, FU10, D10, UK10,
        PK10, GK10, FT10, TK10, KT10, TU10, UT10, FK10, KF10,
        FA11,
        FR12, GP12, FQ12, FF12, FZ12, UQ12, FP12, FS12, UK12, UA12, AK12,
        FK12, KB12, AF12, NestOfTwistedCubes,
        XY13,
        N3, N3pp, UP14, VP14, FV14, OW14, FM14,
        FA15,
        N4
    )

    lst = [RelaxedNonFano, TippedFree3spike,  # 7
           AG23minusDY, TQ8, P8p, KP8, Sp8, Sp8pp, LP8, WQ8,  # 8
           BB9, TQ9, TQ9p, M8591, PP9, BB9gDY, A9, FN9, FX9, KR9, KQ9,  # 9
           UG10, FF10, GP10, FZ10, UQ10, FP10, TQ10, FY10, PP10, FU10, D10,
           UK10, PK10, GK10, FT10, TK10, KT10, TU10, UT10, FK10, KF10,  # 10
           FA11,  # 11
           FR12, GP12, FQ12, FF12, FZ12, UQ12, FP12, FS12, UK12, UA12, AK12,
           FK12, KB12, AF12, NestOfTwistedCubes,  # 12
           XY13,  # 13
           N3, N3pp, UP14, VP14, FV14, OW14, FM14,  # 14
           FA15,  # 15
           N4]  # 16
    for M in lst:
        yield M()


def VariousMatroids():
    """
    Return an iterator over various other named matroids.

    EXAMPLES::

        sage: VM = list(matroids.VariousMatroids()); len(VM)
        16
        sage: import random
        sage: M = random.choice(VM)
        sage: M.is_valid()
        True

    .. SEEALSO::

        :mod:`Matroid catalog <sage.matroids.matroids_catalog>`, under
        ``Collection of various matroids``.
    """
    from sage.matroids.database_matroids import (
        NonVamos, NotP8, AG23minus,
        P9, R9A, R9B, Block_9_4, TicTacToe,
        N1, Block_10_5, Q10,
        BetsyRoss,
        N2,
        D16, Terrahawk,
        ExtendedBinaryGolayCode
    )

    lst = [NonVamos, NotP8, AG23minus,  # 8
           P9, R9A, R9B, Block_9_4, TicTacToe,  # 9
           N1, Block_10_5, Q10,  # 10
           BetsyRoss,  # 11
           N2,  # 12
           D16, Terrahawk,  # 16
           ExtendedBinaryGolayCode]  # 24
    for M in lst:
        yield M()
