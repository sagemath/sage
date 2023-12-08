r"""
This module contains driver functions to easily access the collections of matroids in the database.
These functions can be viewed by typing ``matroids.`` and hitting the ``tab`` button.

AUTHORS:
- Giorgos Mousa (2023-12-08): initial version

Functions
=========
"""
# **********************************************************************
#       Copyright (C) 2023 Giorgos Mousa <gmousa@proton.me>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# **********************************************************************


def OxleyMatroids():
    """
    Return an iterator of (nonparameterized) matroids as listed in the Appendix ``Some Interesting Matroids`` in [Oxl2011]_.

    EXAMPLES::

        sage: for M in matroids.OxleyMatroids(): # long time
        ....:     assert M.is_valid()
    """
    all = []
    from sage.matroids.database.oxley_matroids import U24, U25, U35, K4, Whirl3, Q6, P6, U36, R6, Fano, FanoDual, NonFano,NonFanoDual, O7, P7, AG32, AG32prime, R8, F8, Q8, L8, S8, Vamos, T8, J, P8, P8pp, Wheel4, Whirl4, K33dual, K33, AG23, TernaryDowling3, Pappus, NonPappus, K5, K5dual, R10, R12, T12, PG23
    lst = { 4 : [ U24 ],
            5 : [ U25, U35 ],
            6 : [ K4, Whirl3, Q6, P6, U36, R6 ],
            7 : [ Fano, FanoDual, NonFano,NonFanoDual, O7, P7 ],
            8 : [ AG32, AG32prime, R8, F8, Q8, L8, S8, Vamos, T8, J, P8, P8pp, Wheel4, Whirl4 ],
            9 : [ K33dual, K33, AG23, TernaryDowling3, Pappus, NonPappus ],
            10 : [ K5, K5dual, R10 ],
            # NonDesargues,
            12: [ R12,
            # S_5_6_12,
            T12 ],
            13: [ PG23 ]}
    for i in lst:
        for M in lst[i]:
            all.append(M())
    return iter(all)


def BrettellMatroids():
    """
    Return an iterator of interesting matroids as listed in [].

        EXAMPLES::

        sage: for M in matroids.OxleyMatroids(): # long time
        ....:     assert M.is_valid()
    """
    all = []
    from sage.matroids.database.brettell_matroids import RelaxedNonFano, TippedFree3spike, AG23minusDY, TQ8, P8p, KP8, Sp8, Sp8pp, LP8, WQ8, BB9, TQ9, TQ9p, M8591, PP9, BB9gDY, A9, FN9, FX9, KR9, KQ9, UG10, FF10, GP10, FZ10, UQ10, FP10, TQ10, FY10, PP10, FU10, D10, UK10, PK10, GK10, FT10, TK10, KT10, TU10, UT10, FK10, KF10, FA11, FR12, GP12, FQ12, FF12, FZ12, UQ12, FP12, FS12, UK12, UA12, AK12, FK12, KB12, AF12, NestOfTwistedCubes, XY13, N3, N3pp, UP14, VP14, FV14, OW14, FM14, FA15, N4
    lst = { 7 : [ RelaxedNonFano, TippedFree3spike ],
            8 : [ AG23minusDY, TQ8, P8p, KP8, Sp8, Sp8pp, LP8, WQ8 ],
            9 : [ BB9, TQ9, TQ9p, M8591, PP9, BB9gDY, A9, FN9, FX9, KR9, KQ9 ],
            10 : [ UG10, FF10, GP10, FZ10, UQ10, FP10, TQ10, FY10, PP10, FU10, D10, UK10, PK10, GK10, FT10, TK10, KT10, TU10, UT10, FK10, KF10 ],
            11 : [ FA11 ],
            12 : [ FR12, GP12, FQ12, FF12, FZ12, UQ12, FP12, FS12, UK12, UA12, AK12, FK12, KB12, AF12, NestOfTwistedCubes ],
            13 : [ XY13 ],
            14 : [ N3, N3pp, UP14, VP14, FV14, OW14, FM14 ],
            15 : [ FA15 ],
            16 : [ N4 ]}
    for i in lst:
        for M in lst[i]:
            all.append(M())
    return iter(all)


def VariousMatroids():
    """
    Return an iterator of various other named matroids.

    EXAMPLES::

        sage: for M in matroids.VariousMatroids(): # long time
        ....:     assert M.is_valid()
    """
    all = []
    from sage.matroids.database.various_matroids import NonVamos, TicTacToe, Q10, N1, N2, BetsyRoss, Block_9_4, Block_10_5, ExtendedBinaryGolayCode, ExtendedTernaryGolayCode, AG23minus, NotP8, D16, Terrahawk, R9A, R9B, P9
    lst = { 8 : [ NonVamos, NotP8, AG23minus ],
            9 : [ P9, R9A, R9B, Block_9_4, TicTacToe ],
            10 : [ N1, Block_10_5, Q10 ],
            11 : [ BetsyRoss ],
            12 : [ N2, ExtendedTernaryGolayCode ],
            16 : [ D16, Terrahawk ],
            24 : [ ExtendedBinaryGolayCode ] }
    for i in lst:
        for M in lst[i]:
            all.append(M())
    return iter(all)