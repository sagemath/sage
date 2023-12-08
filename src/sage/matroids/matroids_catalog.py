r"""
Catalog of matroids

A module containing constructors for several common matroids.

A list of all matroids in this module is available via tab completion.
Let ``<tab>`` indicate pressing the :kbd:`Tab` key. So begin by typing
``matroids.<tab>`` to see the various constructions available. Many special
matroids can be accessed from the submenu ``matroids.catalog.<tab>``.

To create a custom matroid using a variety of inputs, see the function
:func:`Matroid() <sage.matroids.constructor.Matroid>`.

- Parametrized matroid constructors (``matroids.<tab>``)

    - :func:`matroids.AG <sage.matroids.database.oxley_matroids.AG>`
    - :func:`matroids.CompleteGraphic <sage.matroids.database.various_matroids.CompleteGraphic>`
    - :func:`matroids.FreeSpike <sage.matroids.database.brettell_matroids.FreeSpike>`
    - :func:`matroids.CompleteGraphic <sage.matroids.database.various_matroids.CompleteGraphic>`
    - :func:`matroids.PG <sage.matroids.database.oxley_matroids.PG>`
    - :func:`matroids.TippedFreeSpike <sage.matroids.database.brettell_matroids.TippedFreeSpike>`
    - :func:`matroids.Uniform <sage.matroids.database.oxley_matroids.Uniform>`
    - :func:`matroids.Wheel <sage.matroids.database.oxley_matroids.Wheel>`
    - :func:`matroids.Whirl <sage.matroids.database.oxley_matroids.Whirl>`


- List of collections of matroids (``matroids.<tab>``)

    - :func:`matroids.oxley_matroids <sage.matroids.database.oxley_matroids>`
    - :func:`matroids.brettell_matroids <sage.matroids.database.brettell_matroids>`
    - :func:`matroids.various_matroids <sage.matroids.database.various_matroids>`


- List of matroids in the catalog (``matroids.catalog.<tab>``)

    - :func:`matroids.catalog.U24 <sage.matroids.database.oxley_matroids.U24>`
    - :func:`matroids.catalog.U25 <sage.matroids.database.oxley_matroids.U25>`
    - :func:`matroids.catalog.U35 <sage.matroids.database.oxley_matroids.U35>`
    - :func:`matroids.catalog.K4 <sage.matroids.database.oxley_matroids.K4>`
    - :func:`matroids.catalog.Whirl3 <sage.matroids.database.oxley_matroids.Whirl3>`
    - :func:`matroids.catalog.Q6 <sage.matroids.database.oxley_matroids.Q6>`
    - :func:`matroids.catalog.P6 <sage.matroids.database.oxley_matroids.P6>`
    - :func:`matroids.catalog.U36 <sage.matroids.database.oxley_matroids.U36>`
    - :func:`matroids.catalog.R6 <sage.matroids.database.oxley_matroids.R6>`
    - :func:`matroids.catalog.Fano <sage.matroids.database.oxley_matroids.Fano>`
    - :func:`matroids.catalog.FanoDual <sage.matroids.database.oxley_matroids.FanoDual>`
    - :func:`matroids.catalog.NonFano <sage.matroids.database.oxley_matroids.NonFano>`
    - :func:`matroids.catalog.NonFanoDual <sage.matroids.database.oxley_matroids.NonFanoDual>`
    - :func:`matroids.catalog.O7 <sage.matroids.database.oxley_matroids.O7>`
    - :func:`matroids.catalog.P7 <sage.matroids.database.oxley_matroids.P7>`
    - :func:`matroids.catalog.AG32 <sage.matroids.database.oxley_matroids.AG32>`
    - :func:`matroids.catalog.AG32prime <sage.matroids.database.oxley_matroids.AG32prime>`
    - :func:`matroids.catalog.R8 <sage.matroids.database.oxley_matroids.R8>`
    - :func:`matroids.catalog.F8 <sage.matroids.database.oxley_matroids.F8>`
    - :func:`matroids.catalog.Q8 <sage.matroids.database.oxley_matroids.Q8>`
    - :func:`matroids.catalog.L8 <sage.matroids.database.oxley_matroids.L8>`
    - :func:`matroids.catalog.S8 <sage.matroids.database.oxley_matroids.S8>`
    - :func:`matroids.catalog.Vamos <sage.matroids.database.oxley_matroids.Vamos>`
    - :func:`matroids.catalog.T8 <sage.matroids.database.oxley_matroids.T8>`
    - :func:`matroids.catalog.J <sage.matroids.database.oxley_matroids.J>`
    - :func:`matroids.catalog.P8 <sage.matroids.database.oxley_matroids.P8>`
    - :func:`matroids.catalog.P8pp <sage.matroids.database.oxley_matroids.P8pp>`
    - :func:`matroids.catalog.Wheel4 <sage.matroids.database.oxley_matroids.Wheel4>`
    - :func:`matroids.catalog.Whirl4 <sage.matroids.database.oxley_matroids.Whirl4>`
    - :func:`matroids.catalog.K33dual <sage.matroids.database.oxley_matroids.K33dual>`
    - :func:`matroids.catalog.K33 <sage.matroids.database.oxley_matroids.K33>`
    - :func:`matroids.catalog.AG23 <sage.matroids.database.oxley_matroids.AG23>`
    - :func:`matroids.catalog.TernaryDowling3 <sage.matroids.database.oxley_matroids.TernaryDowling3>`
    - :func:`matroids.catalog.Pappus <sage.matroids.database.oxley_matroids.Pappus>`
    - :func:`matroids.catalog.NonPappus <sage.matroids.database.oxley_matroids.NonPappus>`
    - :func:`matroids.catalog.K5 <sage.matroids.database.oxley_matroids.K5>`
    - :func:`matroids.catalog.K5dual <sage.matroids.database.oxley_matroids.K5dual>`
    - :func:`matroids.catalog.R10 <sage.matroids.database.oxley_matroids.R10>`
    - :func:`matroids.catalog.R12 <sage.matroids.database.oxley_matroids.R12>`
    - :func:`matroids.catalog.T12 <sage.matroids.database.oxley_matroids.T12>`
    - :func:`matroids.catalog.PG23 <sage.matroids.database.oxley_matroids.PG23>`

    - :func:`matroids.catalog.relaxedNonFano <sage.matroids.database.brettell_matroids.relaxedNonFano>`
    - :func:`matroids.catalog.tippedFree3spike <sage.matroids.database.brettell_matroids.tippedFree3spike>`
    - :func:`matroids.catalog.AG23minusDY <sage.matroids.database.brettell_matroids.AG23minusDY>`
    - :func:`matroids.catalog.TQ8 <sage.matroids.database.brettell_matroids.TQ8>`
    - :func:`matroids.catalog.P8p <sage.matroids.database.brettell_matroids.P8p>`
    - :func:`matroids.catalog.KP8 <sage.matroids.database.brettell_matroids.KP8>`
    - :func:`matroids.catalog.Sp8 <sage.matroids.database.brettell_matroids.Sp8>`
    - :func:`matroids.catalog.Sp8pp <sage.matroids.database.brettell_matroids.Sp8pp>`
    - :func:`matroids.catalog.LP8 <sage.matroids.database.brettell_matroids.LP8>`
    - :func:`matroids.catalog.WQ8 <sage.matroids.database.brettell_matroids.WQ8>`
    - :func:`matroids.catalog.BB9 <sage.matroids.database.brettell_matroids.BB9>`
    - :func:`matroids.catalog.TQ9 <sage.matroids.database.brettell_matroids.TQ9>`
    - :func:`matroids.catalog.TQ9p <sage.matroids.database.brettell_matroids.TQ9p>`
    - :func:`matroids.catalog.M8591 <sage.matroids.database.brettell_matroids.M8591>`
    - :func:`matroids.catalog.PP9 <sage.matroids.database.brettell_matroids.PP9>`
    - :func:`matroids.catalog.BB9gDY <sage.matroids.database.brettell_matroids.BB9gDY>`
    - :func:`matroids.catalog.A9 <sage.matroids.database.brettell_matroids.A9>`
    - :func:`matroids.catalog.FN9 <sage.matroids.database.brettell_matroids.FN9>`
    - :func:`matroids.catalog.FX9 <sage.matroids.database.brettell_matroids.FX9>`
    - :func:`matroids.catalog.KR9 <sage.matroids.database.brettell_matroids.KR9>`
    - :func:`matroids.catalog.KQ9 <sage.matroids.database.brettell_matroids.KQ9>`
    - :func:`matroids.catalog.UG10 <sage.matroids.database.brettell_matroids.UG10>`
    - :func:`matroids.catalog.FF10 <sage.matroids.database.brettell_matroids.FF10>`
    - :func:`matroids.catalog.GP10 <sage.matroids.database.brettell_matroids.GP10>`
    - :func:`matroids.catalog.FZ10 <sage.matroids.database.brettell_matroids.FZ10>`
    - :func:`matroids.catalog.UQ10 <sage.matroids.database.brettell_matroids.UQ10>`
    - :func:`matroids.catalog.FP10 <sage.matroids.database.brettell_matroids.FP10>`
    - :func:`matroids.catalog.TQ10 <sage.matroids.database.brettell_matroids.TQ10>`
    - :func:`matroids.catalog.FY10 <sage.matroids.database.brettell_matroids.FY10>`
    - :func:`matroids.catalog.PP10 <sage.matroids.database.brettell_matroids.PP10>`
    - :func:`matroids.catalog.FU10 <sage.matroids.database.brettell_matroids.FU10>`
    - :func:`matroids.catalog.D10 <sage.matroids.database.brettell_matroids.D10>`
    - :func:`matroids.catalog.UK10 <sage.matroids.database.brettell_matroids.UK10>`
    - :func:`matroids.catalog.PK10 <sage.matroids.database.brettell_matroids.PK10>`
    - :func:`matroids.catalog.GK10 <sage.matroids.database.brettell_matroids.GK10>`
    - :func:`matroids.catalog.FT10 <sage.matroids.database.brettell_matroids.FT10>`
    - :func:`matroids.catalog.TK10 <sage.matroids.database.brettell_matroids.TK10>`
    - :func:`matroids.catalog.KT10 <sage.matroids.database.brettell_matroids.KT10>`
    - :func:`matroids.catalog.TU10 <sage.matroids.database.brettell_matroids.TU10>`
    - :func:`matroids.catalog.UT10 <sage.matroids.database.brettell_matroids.UT10>`
    - :func:`matroids.catalog.FK10 <sage.matroids.database.brettell_matroids.FK10>`
    - :func:`matroids.catalog.KF10 <sage.matroids.database.brettell_matroids.KF10>`
    - :func:`matroids.catalog.FA11 <sage.matroids.database.brettell_matroids.FA11>`
    - :func:`matroids.catalog.FR12 <sage.matroids.database.brettell_matroids.FR12>`
    - :func:`matroids.catalog.GP12 <sage.matroids.database.brettell_matroids.GP12>`
    - :func:`matroids.catalog.FQ12 <sage.matroids.database.brettell_matroids.FQ12>`
    - :func:`matroids.catalog.FF12 <sage.matroids.database.brettell_matroids.FF12>`
    - :func:`matroids.catalog.FZ12 <sage.matroids.database.brettell_matroids.FZ12>`
    - :func:`matroids.catalog.UQ12 <sage.matroids.database.brettell_matroids.UQ12>`
    - :func:`matroids.catalog.FP12 <sage.matroids.database.brettell_matroids.FP12>`
    - :func:`matroids.catalog.FS12 <sage.matroids.database.brettell_matroids.FS12>`
    - :func:`matroids.catalog.UK12 <sage.matroids.database.brettell_matroids.UK12>`
    - :func:`matroids.catalog.UA12 <sage.matroids.database.brettell_matroids.UA12>`
    - :func:`matroids.catalog.AK12 <sage.matroids.database.brettell_matroids.AK12>`
    - :func:`matroids.catalog.FK12 <sage.matroids.database.brettell_matroids.FK12>`
    - :func:`matroids.catalog.KB12 <sage.matroids.database.brettell_matroids.KB12>`
    - :func:`matroids.catalog.AF12 <sage.matroids.database.brettell_matroids.AF12>`
    - :func:`matroids.catalog.nestOfTwistedCubes <sage.matroids.database.brettell_matroids.nestOfTwistedCubes>`
    - :func:`matroids.catalog.XY13 <sage.matroids.database.brettell_matroids.XY13>`
    - :func:`matroids.catalog.N3 <sage.matroids.database.brettell_matroids.N3>`
    - :func:`matroids.catalog.N3pp <sage.matroids.database.brettell_matroids.N3pp>`
    - :func:`matroids.catalog.UP14 <sage.matroids.database.brettell_matroids.UP14>`
    - :func:`matroids.catalog.VP14 <sage.matroids.database.brettell_matroids.VP14>`
    - :func:`matroids.catalog.FV14 <sage.matroids.database.brettell_matroids.FV14>`
    - :func:`matroids.catalog.OW14 <sage.matroids.database.brettell_matroids.OW14>`
    - :func:`matroids.catalog.FM14 <sage.matroids.database.brettell_matroids.FM14>`
    - :func:`matroids.catalog.FA15 <sage.matroids.database.brettell_matroids.FA15>`
    - :func:`matroids.catalog.N4 <sage.matroids.database.brettell_matroids.N4>`

    - :func:`matroids.catalog.NonVamos <sage.matroids.database.various_matroids.NonVamos>`
    - :func:`matroids.catalog.TicTacToe <sage.matroids.database.various_matroids.TicTacToe>`
    - :func:`matroids.catalog.Q10 <sage.matroids.database.various_matroids.Q10>`
    - :func:`matroids.catalog.N1 <sage.matroids.database.various_matroids.N1>`
    - :func:`matroids.catalog.N2 <sage.matroids.database.various_matroids.N2>`
    - :func:`matroids.catalog.BetsyRoss <sage.matroids.database.various_matroids.BetsyRoss>`
    - :func:`matroids.catalog.Block_9_4 <sage.matroids.database.various_matroids.Block_9_4>`
    - :func:`matroids.catalog.Block_10_5 <sage.matroids.database.various_matroids.Block_10_5>`
    - :func:`matroids.catalog.ExtendedBinaryGolayCode <sage.matroids.database.various_matroids.ExtendedBinaryGolayCode>`
    - :func:`matroids.catalog.ExtendedTernaryGolayCode <sage.matroids.database.various_matroids.ExtendedTernaryGolayCode>`
    - :func:`matroids.catalog.AG23minus <sage.matroids.database.various_matroids.AG23minus>`
    - :func:`matroids.catalog.NotP8 <sage.matroids.database.various_matroids.NotP8>`
    - :func:`matroids.catalog.D16 <sage.matroids.database.various_matroids.D16>`
    - :func:`matroids.catalog.Terrahawk <sage.matroids.database.various_matroids.Terrahawk>`
    - :func:`matroids.catalog.R9A <sage.matroids.database.various_matroids.R9A>`
    - :func:`matroids.catalog.R9B <sage.matroids.database.various_matroids.R9B>`
    - :func:`matroids.catalog.P9 <sage.matroids.database.various_matroids.P9>`
"""

# Do not add code to this file, only imports.
# Workaround for help in the notebook (needs parentheses in this file)

# user-accessible:
from sage.misc.lazy_import import lazy_import
lazy_import('sage.matroids.database.driver_functions', ('BrettellMatroids', 'OxleyMatroids', 'VariousMatroids'))
lazy_import('sage.matroids.database.oxley_matroids', ('Wheel', 'Whirl', 'Uniform', 'PG', 'AG'))
lazy_import('sage.matroids.database.brettell_matroids', 'FreeSpike')
lazy_import('sage.matroids.database.various_matroids', 'CompleteGraphic')
lazy_import('sage.matroids', 'catalog')
lazy_import('sage.matroids', 'named_matroids')
