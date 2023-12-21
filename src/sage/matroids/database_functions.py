r"""
Collections of matroids

This module contains functions that access the collections of matroids in the
database. Each of these functions returns a complete list of the
nonparametrized matroids from the corresponding collection. These functions
can be viewed by typing ``matroids.`` + :kbd:`Tab`.

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
    Return the list of Oxley's matroids.

    EXAMPLES::

        sage: import random
        sage: OM = matroids.OxleyMatroids(); len(OM)
        42
        sage: M = random.choice(OM)
        sage: M.is_valid() # long time
        True

    .. SEEALSO::

        :mod:`matroids.database_matroids <sage.matroids.database_matroids>`

    REFERENCES:

    These matroids are the nonparametrized matroids that appear in the
    Appendix ``Some Interesting Matroids`` in [Oxl2011]_ (p. 639-64).

    """
    Matroids = []
    from sage.matroids.database_matroids import (
        U24, U25, U35, K4, Whirl3, Q6, P6, U36, R6,
        Fano, FanoDual, NonFano, NonFanoDual, O7, P7,
        AG32, AG32prime, R8, F8, Q8, L8, S8, Vamos, T8, J, P8, P8pp,
        Wheel4, Whirl4,
        K33dual, K33, AG23, TernaryDowling3,  # R9
        Pappus, NonPappus,
        K5, K5dual, R10,  # NonDesargues,
        R12, ExtendedTernaryGolayCode, T12,
        PG23
    )

    lst = {
        4: [U24],
        5: [U25, U35],
        6: [K4, Whirl3, Q6, P6, U36, R6],
        7: [Fano, FanoDual, NonFano, NonFanoDual, O7, P7],
        8: [
            AG32, AG32prime,
            R8, F8, Q8, L8, S8,
            Vamos, T8, J, P8, P8pp,
            Wheel4, Whirl4
        ],
        9: [K33dual, K33, AG23, TernaryDowling3, Pappus, NonPappus],
        10: [K5, K5dual, R10],
        12: [R12, ExtendedTernaryGolayCode, T12],
        13: [PG23],
    }
    for i in lst:
        for M in lst[i]:
            Matroids.append(M())
    return Matroids


def BrettellMatroids():
    """
    Return the list of Brettell's matroids.

    EXAMPLES::

        sage: BM = matroids.BrettellMatroids(); len(BM)
        68
        sage: import random
        sage: M = random.choice(BM)
        sage: M.is_valid() # long time
        True

    .. SEEALSO::

        :mod:`matroids.database_matroids <sage.matroids.database_matroids>`

    """
    Matroids = []
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

    lst = {
        7: [RelaxedNonFano, TippedFree3spike],
        8: [AG23minusDY, TQ8, P8p, KP8, Sp8, Sp8pp, LP8, WQ8],
        9: [BB9, TQ9, TQ9p, M8591, PP9, BB9gDY, A9, FN9, FX9, KR9, KQ9],
        10: [
            UG10, FF10, GP10, FZ10, UQ10, FP10, TQ10, FY10, PP10, FU10, D10,
            UK10, PK10, GK10, FT10, TK10, KT10, TU10, UT10, FK10, KF10
        ],
        11: [FA11],
        12: [
            FR12, GP12, FQ12, FF12, FZ12, UQ12, FP12, FS12, UK12, UA12, AK12,
            FK12, KB12, AF12, NestOfTwistedCubes
        ],
        13: [XY13],
        14: [N3, N3pp, UP14, VP14, FV14, OW14, FM14],
        15: [FA15],
        16: [N4],
    }
    for i in lst:
        for M in lst[i]:
            Matroids.append(M())
    return Matroids


def VariousMatroids():
    """
    Return a list of various other named matroids.

    EXAMPLES::

        sage: import random
        sage: VM = matroids.VariousMatroids(); len(VM)
        16
        sage: M = random.choice(VM)
        sage: M.is_valid() # long time
        True

    .. SEEALSO::

        :mod:`matroids.database_matroids <sage.matroids.database_matroids>`

    """
    Matroids = []
    from sage.matroids.database_matroids import (
        NonVamos, NotP8, AG23minus,
        P9, R9A, R9B, Block_9_4, TicTacToe,
        N1, Block_10_5, Q10,
        BetsyRoss,
        N2,
        D16, Terrahawk,
        ExtendedBinaryGolayCode
    )

    lst = {
        8: [NonVamos, NotP8, AG23minus],
        9: [P9, R9A, R9B, Block_9_4, TicTacToe],
        10: [N1, Block_10_5, Q10],
        11: [BetsyRoss],
        12: [N2],
        16: [D16, Terrahawk],
        24: [ExtendedBinaryGolayCode],
    }
    for i in lst:
        for M in lst[i]:
            Matroids.append(M())
    return Matroids
