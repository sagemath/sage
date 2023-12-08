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

    - :func:`matroids.AG <sage.matroids.database.OxleyMatroids.AG>`
    - :func:`matroids.CompleteGraphic <sage.matroids.database.VariousMatroids.CompleteGraphic>`
    - :func:`matroids.FreeSpike <sage.matroids.database.BrettellMatroids.FreeSpike>`
    - :func:`matroids.CompleteGraphic <sage.matroids.database.VariousMatroids.CompleteGraphic>`
    - :func:`matroids.PG <sage.matroids.database.OxleyMatroids.PG>`
    - :func:`matroids.TippedFreeSpike <sage.matroids.database.BrettellMatroids.TippedFreeSpike>`
    - :func:`matroids.Uniform <sage.matroids.database.OxleyMatroids.Uniform>`
    - :func:`matroids.Wheel <sage.matroids.database.OxleyMatroids.Wheel>`
    - :func:`matroids.Whirl <sage.matroids.database.OxleyMatroids.Whirl>`


- List of collections of matroids (``matroids.<tab>``)

    - :func:`matroids.OxleyMatroids <sage.matroids.database.OxleyMatroids>`
    - :func:`matroids.BrettellMatroids <sage.matroids.database.BrettellMatroids>`
    - :func:`matroids.VariousMatroids <sage.matroids.database.VariousMatroids>`


- List of matroids in the catalog (``matroids.catalog.<tab>``)

    - :func:`matroids.catalog.U24 <sage.matroids.database.OxleyMatroids.U24>`
    - :func:`matroids.catalog.U25 <sage.matroids.database.OxleyMatroids.U25>`
    - :func:`matroids.catalog.U35 <sage.matroids.database.OxleyMatroids.U35>`
    - :func:`matroids.catalog.K4 <sage.matroids.database.OxleyMatroids.K4>`
    - :func:`matroids.catalog.Whirl3 <sage.matroids.database.OxleyMatroids.Whirl3>`
    - :func:`matroids.catalog.Q6 <sage.matroids.database.OxleyMatroids.Q6>`
    - :func:`matroids.catalog.P6 <sage.matroids.database.OxleyMatroids.P6>`
    - :func:`matroids.catalog.U36 <sage.matroids.database.OxleyMatroids.U36>`
    - :func:`matroids.catalog.R6 <sage.matroids.database.OxleyMatroids.R6>`
    - :func:`matroids.catalog.Fano <sage.matroids.database.OxleyMatroids.Fano>`
    - :func:`matroids.catalog.FanoDual <sage.matroids.database.OxleyMatroids.FanoDual>`
    - :func:`matroids.catalog.NonFano <sage.matroids.database.OxleyMatroids.NonFano>`
    - :func:`matroids.catalog.NonFanoDual <sage.matroids.database.OxleyMatroids.NonFanoDual>`
    - :func:`matroids.catalog.O7 <sage.matroids.database.OxleyMatroids.O7>`
    - :func:`matroids.catalog.P7 <sage.matroids.database.OxleyMatroids.P7>`
    - :func:`matroids.catalog.AG32 <sage.matroids.database.OxleyMatroids.AG32>`
    - :func:`matroids.catalog.AG32prime <sage.matroids.database.OxleyMatroids.AG32prime>`
    - :func:`matroids.catalog.R8 <sage.matroids.database.OxleyMatroids.R8>`
    - :func:`matroids.catalog.F8 <sage.matroids.database.OxleyMatroids.F8>`
    - :func:`matroids.catalog.Q8 <sage.matroids.database.OxleyMatroids.Q8>`
    - :func:`matroids.catalog.L8 <sage.matroids.database.OxleyMatroids.L8>`
    - :func:`matroids.catalog.S8 <sage.matroids.database.OxleyMatroids.S8>`
    - :func:`matroids.catalog.Vamos <sage.matroids.database.OxleyMatroids.Vamos>`
    - :func:`matroids.catalog.T8 <sage.matroids.database.OxleyMatroids.T8>`
    - :func:`matroids.catalog.J <sage.matroids.database.OxleyMatroids.J>`
    - :func:`matroids.catalog.P8 <sage.matroids.database.OxleyMatroids.P8>`
    - :func:`matroids.catalog.P8pp <sage.matroids.database.OxleyMatroids.P8pp>`
    - :func:`matroids.catalog.Wheel4 <sage.matroids.database.OxleyMatroids.Wheel4>`
    - :func:`matroids.catalog.Whirl4 <sage.matroids.database.OxleyMatroids.Whirl4>`
    - :func:`matroids.catalog.K33dual <sage.matroids.database.OxleyMatroids.K33dual>`
    - :func:`matroids.catalog.K33 <sage.matroids.database.OxleyMatroids.K33>`
    - :func:`matroids.catalog.AG23 <sage.matroids.database.OxleyMatroids.AG23>`
    - :func:`matroids.catalog.TernaryDowling3 <sage.matroids.database.OxleyMatroids.TernaryDowling3>`
    - :func:`matroids.catalog.Pappus <sage.matroids.database.OxleyMatroids.Pappus>`
    - :func:`matroids.catalog.NonPappus <sage.matroids.database.OxleyMatroids.NonPappus>`
    - :func:`matroids.catalog.K5 <sage.matroids.database.OxleyMatroids.K5>`
    - :func:`matroids.catalog.K5dual <sage.matroids.database.OxleyMatroids.K5dual>`
    - :func:`matroids.catalog.R10 <sage.matroids.database.OxleyMatroids.R10>`
    - :func:`matroids.catalog.R12 <sage.matroids.database.OxleyMatroids.R12>`
    - :func:`matroids.catalog.T12 <sage.matroids.database.OxleyMatroids.T12>`
    - :func:`matroids.catalog.PG23 <sage.matroids.database.OxleyMatroids.PG23>`

    - :func:`matroids.catalog.relaxedNonFano <sage.matroids.database.BrettellMatroids.relaxedNonFano>`
    - :func:`matroids.catalog.tippedFree3spike <sage.matroids.database.BrettellMatroids.tippedFree3spike>`
    - :func:`matroids.catalog.AG23minusDY <sage.matroids.database.BrettellMatroids.AG23minusDY>`
    - :func:`matroids.catalog.TQ8 <sage.matroids.database.BrettellMatroids.TQ8>`
    - :func:`matroids.catalog.P8p <sage.matroids.database.BrettellMatroids.P8p>`
    - :func:`matroids.catalog.KP8 <sage.matroids.database.BrettellMatroids.KP8>`
    - :func:`matroids.catalog.Sp8 <sage.matroids.database.BrettellMatroids.Sp8>`
    - :func:`matroids.catalog.Sp8pp <sage.matroids.database.BrettellMatroids.Sp8pp>`
    - :func:`matroids.catalog.LP8 <sage.matroids.database.BrettellMatroids.LP8>`
    - :func:`matroids.catalog.WQ8 <sage.matroids.database.BrettellMatroids.WQ8>`
    - :func:`matroids.catalog.BB9 <sage.matroids.database.BrettellMatroids.BB9>`
    - :func:`matroids.catalog.TQ9 <sage.matroids.database.BrettellMatroids.TQ9>`
    - :func:`matroids.catalog.TQ9p <sage.matroids.database.BrettellMatroids.TQ9p>`
    - :func:`matroids.catalog.M8591 <sage.matroids.database.BrettellMatroids.M8591>`
    - :func:`matroids.catalog.PP9 <sage.matroids.database.BrettellMatroids.PP9>`
    - :func:`matroids.catalog.BB9gDY <sage.matroids.database.BrettellMatroids.BB9gDY>`
    - :func:`matroids.catalog.A9 <sage.matroids.database.BrettellMatroids.A9>`
    - :func:`matroids.catalog.FN9 <sage.matroids.database.BrettellMatroids.FN9>`
    - :func:`matroids.catalog.FX9 <sage.matroids.database.BrettellMatroids.FX9>`
    - :func:`matroids.catalog.KR9 <sage.matroids.database.BrettellMatroids.KR9>`
    - :func:`matroids.catalog.KQ9 <sage.matroids.database.BrettellMatroids.KQ9>`
    - :func:`matroids.catalog.UG10 <sage.matroids.database.BrettellMatroids.UG10>`
    - :func:`matroids.catalog.FF10 <sage.matroids.database.BrettellMatroids.FF10>`
    - :func:`matroids.catalog.GP10 <sage.matroids.database.BrettellMatroids.GP10>`
    - :func:`matroids.catalog.FZ10 <sage.matroids.database.BrettellMatroids.FZ10>`
    - :func:`matroids.catalog.UQ10 <sage.matroids.database.BrettellMatroids.UQ10>`
    - :func:`matroids.catalog.FP10 <sage.matroids.database.BrettellMatroids.FP10>`
    - :func:`matroids.catalog.TQ10 <sage.matroids.database.BrettellMatroids.TQ10>`
    - :func:`matroids.catalog.FY10 <sage.matroids.database.BrettellMatroids.FY10>`
    - :func:`matroids.catalog.PP10 <sage.matroids.database.BrettellMatroids.PP10>`
    - :func:`matroids.catalog.FU10 <sage.matroids.database.BrettellMatroids.FU10>`
    - :func:`matroids.catalog.D10 <sage.matroids.database.BrettellMatroids.D10>`
    - :func:`matroids.catalog.UK10 <sage.matroids.database.BrettellMatroids.UK10>`
    - :func:`matroids.catalog.PK10 <sage.matroids.database.BrettellMatroids.PK10>`
    - :func:`matroids.catalog.GK10 <sage.matroids.database.BrettellMatroids.GK10>`
    - :func:`matroids.catalog.FT10 <sage.matroids.database.BrettellMatroids.FT10>`
    - :func:`matroids.catalog.TK10 <sage.matroids.database.BrettellMatroids.TK10>`
    - :func:`matroids.catalog.KT10 <sage.matroids.database.BrettellMatroids.KT10>`
    - :func:`matroids.catalog.TU10 <sage.matroids.database.BrettellMatroids.TU10>`
    - :func:`matroids.catalog.UT10 <sage.matroids.database.BrettellMatroids.UT10>`
    - :func:`matroids.catalog.FK10 <sage.matroids.database.BrettellMatroids.FK10>`
    - :func:`matroids.catalog.KF10 <sage.matroids.database.BrettellMatroids.KF10>`
    - :func:`matroids.catalog.FA11 <sage.matroids.database.BrettellMatroids.FA11>`
    - :func:`matroids.catalog.FR12 <sage.matroids.database.BrettellMatroids.FR12>`
    - :func:`matroids.catalog.GP12 <sage.matroids.database.BrettellMatroids.GP12>`
    - :func:`matroids.catalog.FQ12 <sage.matroids.database.BrettellMatroids.FQ12>`
    - :func:`matroids.catalog.FF12 <sage.matroids.database.BrettellMatroids.FF12>`
    - :func:`matroids.catalog.FZ12 <sage.matroids.database.BrettellMatroids.FZ12>`
    - :func:`matroids.catalog.UQ12 <sage.matroids.database.BrettellMatroids.UQ12>`
    - :func:`matroids.catalog.FP12 <sage.matroids.database.BrettellMatroids.FP12>`
    - :func:`matroids.catalog.FS12 <sage.matroids.database.BrettellMatroids.FS12>`
    - :func:`matroids.catalog.UK12 <sage.matroids.database.BrettellMatroids.UK12>`
    - :func:`matroids.catalog.UA12 <sage.matroids.database.BrettellMatroids.UA12>`
    - :func:`matroids.catalog.AK12 <sage.matroids.database.BrettellMatroids.AK12>`
    - :func:`matroids.catalog.FK12 <sage.matroids.database.BrettellMatroids.FK12>`
    - :func:`matroids.catalog.KB12 <sage.matroids.database.BrettellMatroids.KB12>`
    - :func:`matroids.catalog.AF12 <sage.matroids.database.BrettellMatroids.AF12>`
    - :func:`matroids.catalog.nestOfTwistedCubes <sage.matroids.database.BrettellMatroids.nestOfTwistedCubes>`
    - :func:`matroids.catalog.XY13 <sage.matroids.database.BrettellMatroids.XY13>`
    - :func:`matroids.catalog.N3 <sage.matroids.database.BrettellMatroids.N3>`
    - :func:`matroids.catalog.N3pp <sage.matroids.database.BrettellMatroids.N3pp>`
    - :func:`matroids.catalog.UP14 <sage.matroids.database.BrettellMatroids.UP14>`
    - :func:`matroids.catalog.VP14 <sage.matroids.database.BrettellMatroids.VP14>`
    - :func:`matroids.catalog.FV14 <sage.matroids.database.BrettellMatroids.FV14>`
    - :func:`matroids.catalog.OW14 <sage.matroids.database.BrettellMatroids.OW14>`
    - :func:`matroids.catalog.FM14 <sage.matroids.database.BrettellMatroids.FM14>`
    - :func:`matroids.catalog.FA15 <sage.matroids.database.BrettellMatroids.FA15>`
    - :func:`matroids.catalog.N4 <sage.matroids.database.BrettellMatroids.N4>`

    - :func:`matroids.catalog.NonVamos <sage.matroids.database.VariousMatroids.NonVamos>`
    - :func:`matroids.catalog.TicTacToe <sage.matroids.database.VariousMatroids.TicTacToe>`
    - :func:`matroids.catalog.Q10 <sage.matroids.database.VariousMatroids.Q10>`
    - :func:`matroids.catalog.N1 <sage.matroids.database.VariousMatroids.N1>`
    - :func:`matroids.catalog.N2 <sage.matroids.database.VariousMatroids.N2>`
    - :func:`matroids.catalog.BetsyRoss <sage.matroids.database.VariousMatroids.BetsyRoss>`
    - :func:`matroids.catalog.Block_9_4 <sage.matroids.database.VariousMatroids.Block_9_4>`
    - :func:`matroids.catalog.Block_10_5 <sage.matroids.database.VariousMatroids.Block_10_5>`
    - :func:`matroids.catalog.ExtendedBinaryGolayCode <sage.matroids.database.VariousMatroids.ExtendedBinaryGolayCode>`
    - :func:`matroids.catalog.ExtendedTernaryGolayCode <sage.matroids.database.VariousMatroids.ExtendedTernaryGolayCode>`
    - :func:`matroids.catalog.AG23minus <sage.matroids.database.VariousMatroids.AG23minus>`
    - :func:`matroids.catalog.NotP8 <sage.matroids.database.VariousMatroids.NotP8>`
    - :func:`matroids.catalog.D16 <sage.matroids.database.VariousMatroids.D16>`
    - :func:`matroids.catalog.Terrahawk <sage.matroids.database.VariousMatroids.Terrahawk>`
    - :func:`matroids.catalog.R9A <sage.matroids.database.VariousMatroids.R9A>`
    - :func:`matroids.catalog.R9B <sage.matroids.database.VariousMatroids.R9B>`
    - :func:`matroids.catalog.P9 <sage.matroids.database.VariousMatroids.P9>`
"""

# Do not add code to this file, only imports.
# Workaround for help in the notebook (needs parentheses in this file)

# user-accessible:
from sage.misc.lazy_import import lazy_import
lazy_import('sage.matroids.database.driver_functions', ('AllMatroids', 'BrettellMatroids', 'OxleyMatroids', 'VariousMatroids'))
lazy_import('sage.matroids.database.OxleyMatroids', ('Wheel', 'Whirl', 'Uniform', 'PG', 'AG'))
lazy_import('sage.matroids.database.BrettellMatroids', 'FreeSpike')
lazy_import('sage.matroids.database.VariousMatroids', 'CompleteGraphic')
lazy_import('sage.matroids', 'catalog')
lazy_import('sage.matroids', 'named_matroids')
