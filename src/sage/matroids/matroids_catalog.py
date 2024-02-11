r"""
Catalog of matroids

A list of parametrized and individual matroids. The individual matroids are
grouped into collections.

All listed matroids are available via tab completion. Simply type ``matroids.``
+ :kbd:`Tab` to see the various constructions available. For individual
matroids type ``matroids.catalog.`` + :kbd:`Tab`.

.. NOTE::

    To create a custom matroid using a variety of inputs, see the function
    :func:`Matroid() <sage.matroids.constructor.Matroid>`.

- Parametrized matroids (``matroids.`` + :kbd:`Tab`)
    - :func:`matroids.AG <sage.matroids.database_matroids.AG>`
    - :func:`matroids.CompleteGraphic <sage.matroids.database_matroids.CompleteGraphic>`
    - :func:`matroids.PG <sage.matroids.database_matroids.PG>`
    - :func:`matroids.Psi <sage.matroids.database_matroids.Psi>`
    - :func:`matroids.Spike <sage.matroids.database_matroids.Spike>`
    - :func:`matroids.Theta <sage.matroids.database_matroids.Theta>`
    - :func:`matroids.Uniform <sage.matroids.database_matroids.Uniform>`
    - :func:`matroids.Wheel <sage.matroids.database_matroids.Wheel>`
    - :func:`matroids.Whirl <sage.matroids.database_matroids.Whirl>`
    - :func:`matroids.Z <sage.matroids.database_matroids.Z>`

- Oxley's matroid collection (``matroids.catalog.`` + :kbd:`Tab`)
    - :func:`matroids.catalog.U24 <sage.matroids.database_matroids.U24>`
    - :func:`matroids.catalog.U25 <sage.matroids.database_matroids.U25>`
    - :func:`matroids.catalog.U35 <sage.matroids.database_matroids.U35>`
    - :func:`matroids.catalog.K4 <sage.matroids.database_matroids.K4>`
    - :func:`matroids.catalog.Whirl3 <sage.matroids.database_matroids.Whirl3>`
    - :func:`matroids.catalog.Q6 <sage.matroids.database_matroids.Q6>`
    - :func:`matroids.catalog.P6 <sage.matroids.database_matroids.P6>`
    - :func:`matroids.catalog.U36 <sage.matroids.database_matroids.U36>`
    - :func:`matroids.catalog.R6 <sage.matroids.database_matroids.R6>`
    - :func:`matroids.catalog.Fano <sage.matroids.database_matroids.Fano>`
    - :func:`matroids.catalog.FanoDual <sage.matroids.database_matroids.FanoDual>`
    - :func:`matroids.catalog.NonFano <sage.matroids.database_matroids.NonFano>`
    - :func:`matroids.catalog.NonFanoDual <sage.matroids.database_matroids.NonFanoDual>`
    - :func:`matroids.catalog.O7 <sage.matroids.database_matroids.O7>`
    - :func:`matroids.catalog.P7 <sage.matroids.database_matroids.P7>`
    - :func:`matroids.catalog.AG32 <sage.matroids.database_matroids.AG32>`
    - :func:`matroids.catalog.AG32prime <sage.matroids.database_matroids.AG32prime>`
    - :func:`matroids.catalog.R8 <sage.matroids.database_matroids.R8>`
    - :func:`matroids.catalog.F8 <sage.matroids.database_matroids.F8>`
    - :func:`matroids.catalog.Q8 <sage.matroids.database_matroids.Q8>`
    - :func:`matroids.catalog.L8 <sage.matroids.database_matroids.L8>`
    - :func:`matroids.catalog.S8 <sage.matroids.database_matroids.S8>`
    - :func:`matroids.catalog.Vamos <sage.matroids.database_matroids.Vamos>`
    - :func:`matroids.catalog.T8 <sage.matroids.database_matroids.T8>`
    - :func:`matroids.catalog.J <sage.matroids.database_matroids.J>`
    - :func:`matroids.catalog.P8 <sage.matroids.database_matroids.P8>`
    - :func:`matroids.catalog.P8pp <sage.matroids.database_matroids.P8pp>`
    - :func:`matroids.catalog.Wheel4 <sage.matroids.database_matroids.Wheel4>`
    - :func:`matroids.catalog.Whirl4 <sage.matroids.database_matroids.Whirl4>`
    - :func:`matroids.catalog.K33dual <sage.matroids.database_matroids.K33dual>`
    - :func:`matroids.catalog.K33 <sage.matroids.database_matroids.K33>`
    - :func:`matroids.catalog.AG23 <sage.matroids.database_matroids.AG23>`
    - :func:`matroids.catalog.TernaryDowling3 <sage.matroids.database_matroids.TernaryDowling3>`
    - :func:`matroids.catalog.R9 <sage.matroids.database_matroids.R9>`
    - :func:`matroids.catalog.Pappus <sage.matroids.database_matroids.Pappus>`
    - :func:`matroids.catalog.NonPappus <sage.matroids.database_matroids.NonPappus>`
    - :func:`matroids.catalog.K5 <sage.matroids.database_matroids.K5>`
    - :func:`matroids.catalog.K5dual <sage.matroids.database_matroids.K5dual>`
    - :func:`matroids.catalog.R10 <sage.matroids.database_matroids.R10>`
    - :func:`matroids.catalog.NonDesargues <sage.matroids.database_matroids.NonDesargues>`
    - :func:`matroids.catalog.R12 <sage.matroids.database_matroids.R12>`
    - :func:`matroids.catalog.ExtendedTernaryGolayCode <sage.matroids.database_matroids.ExtendedTernaryGolayCode>`
    - :func:`matroids.catalog.T12 <sage.matroids.database_matroids.T12>`
    - :func:`matroids.catalog.PG23 <sage.matroids.database_matroids.PG23>`

- Brettell's matroid collection (``matroids.catalog.`` + :kbd:`Tab`)
    - :func:`matroids.catalog.RelaxedNonFano <sage.matroids.database_matroids.RelaxedNonFano>`
    - :func:`matroids.catalog.AG23minusDY <sage.matroids.database_matroids.AG23minusDY>`
    - :func:`matroids.catalog.TQ8 <sage.matroids.database_matroids.TQ8>`
    - :func:`matroids.catalog.P8p <sage.matroids.database_matroids.P8p>`
    - :func:`matroids.catalog.KP8 <sage.matroids.database_matroids.KP8>`
    - :func:`matroids.catalog.Sp8 <sage.matroids.database_matroids.Sp8>`
    - :func:`matroids.catalog.Sp8pp <sage.matroids.database_matroids.Sp8pp>`
    - :func:`matroids.catalog.LP8 <sage.matroids.database_matroids.LP8>`
    - :func:`matroids.catalog.WQ8 <sage.matroids.database_matroids.WQ8>`
    - :func:`matroids.catalog.BB9 <sage.matroids.database_matroids.BB9>`
    - :func:`matroids.catalog.TQ9 <sage.matroids.database_matroids.TQ9>`
    - :func:`matroids.catalog.TQ9p <sage.matroids.database_matroids.TQ9p>`
    - :func:`matroids.catalog.M8591 <sage.matroids.database_matroids.M8591>`
    - :func:`matroids.catalog.PP9 <sage.matroids.database_matroids.PP9>`
    - :func:`matroids.catalog.BB9gDY <sage.matroids.database_matroids.BB9gDY>`
    - :func:`matroids.catalog.A9 <sage.matroids.database_matroids.A9>`
    - :func:`matroids.catalog.FN9 <sage.matroids.database_matroids.FN9>`
    - :func:`matroids.catalog.FX9 <sage.matroids.database_matroids.FX9>`
    - :func:`matroids.catalog.KR9 <sage.matroids.database_matroids.KR9>`
    - :func:`matroids.catalog.KQ9 <sage.matroids.database_matroids.KQ9>`
    - :func:`matroids.catalog.UG10 <sage.matroids.database_matroids.UG10>`
    - :func:`matroids.catalog.FF10 <sage.matroids.database_matroids.FF10>`
    - :func:`matroids.catalog.GP10 <sage.matroids.database_matroids.GP10>`
    - :func:`matroids.catalog.FZ10 <sage.matroids.database_matroids.FZ10>`
    - :func:`matroids.catalog.UQ10 <sage.matroids.database_matroids.UQ10>`
    - :func:`matroids.catalog.FP10 <sage.matroids.database_matroids.FP10>`
    - :func:`matroids.catalog.TQ10 <sage.matroids.database_matroids.TQ10>`
    - :func:`matroids.catalog.FY10 <sage.matroids.database_matroids.FY10>`
    - :func:`matroids.catalog.PP10 <sage.matroids.database_matroids.PP10>`
    - :func:`matroids.catalog.FU10 <sage.matroids.database_matroids.FU10>`
    - :func:`matroids.catalog.D10 <sage.matroids.database_matroids.D10>`
    - :func:`matroids.catalog.UK10 <sage.matroids.database_matroids.UK10>`
    - :func:`matroids.catalog.PK10 <sage.matroids.database_matroids.PK10>`
    - :func:`matroids.catalog.GK10 <sage.matroids.database_matroids.GK10>`
    - :func:`matroids.catalog.FT10 <sage.matroids.database_matroids.FT10>`
    - :func:`matroids.catalog.TK10 <sage.matroids.database_matroids.TK10>`
    - :func:`matroids.catalog.KT10 <sage.matroids.database_matroids.KT10>`
    - :func:`matroids.catalog.TU10 <sage.matroids.database_matroids.TU10>`
    - :func:`matroids.catalog.UT10 <sage.matroids.database_matroids.UT10>`
    - :func:`matroids.catalog.FK10 <sage.matroids.database_matroids.FK10>`
    - :func:`matroids.catalog.KF10 <sage.matroids.database_matroids.KF10>`
    - :func:`matroids.catalog.FA11 <sage.matroids.database_matroids.FA11>`
    - :func:`matroids.catalog.FR12 <sage.matroids.database_matroids.FR12>`
    - :func:`matroids.catalog.GP12 <sage.matroids.database_matroids.GP12>`
    - :func:`matroids.catalog.FQ12 <sage.matroids.database_matroids.FQ12>`
    - :func:`matroids.catalog.FF12 <sage.matroids.database_matroids.FF12>`
    - :func:`matroids.catalog.FZ12 <sage.matroids.database_matroids.FZ12>`
    - :func:`matroids.catalog.UQ12 <sage.matroids.database_matroids.UQ12>`
    - :func:`matroids.catalog.FP12 <sage.matroids.database_matroids.FP12>`
    - :func:`matroids.catalog.FS12 <sage.matroids.database_matroids.FS12>`
    - :func:`matroids.catalog.UK12 <sage.matroids.database_matroids.UK12>`
    - :func:`matroids.catalog.UA12 <sage.matroids.database_matroids.UA12>`
    - :func:`matroids.catalog.AK12 <sage.matroids.database_matroids.AK12>`
    - :func:`matroids.catalog.FK12 <sage.matroids.database_matroids.FK12>`
    - :func:`matroids.catalog.KB12 <sage.matroids.database_matroids.KB12>`
    - :func:`matroids.catalog.AF12 <sage.matroids.database_matroids.AF12>`
    - :func:`matroids.catalog.NestOfTwistedCubes <sage.matroids.database_matroids.NestOfTwistedCubes>`
    - :func:`matroids.catalog.XY13 <sage.matroids.database_matroids.XY13>`
    - :func:`matroids.catalog.N3 <sage.matroids.database_matroids.N3>`
    - :func:`matroids.catalog.N3pp <sage.matroids.database_matroids.N3pp>`
    - :func:`matroids.catalog.UP14 <sage.matroids.database_matroids.UP14>`
    - :func:`matroids.catalog.VP14 <sage.matroids.database_matroids.VP14>`
    - :func:`matroids.catalog.FV14 <sage.matroids.database_matroids.FV14>`
    - :func:`matroids.catalog.OW14 <sage.matroids.database_matroids.OW14>`
    - :func:`matroids.catalog.FM14 <sage.matroids.database_matroids.FM14>`
    - :func:`matroids.catalog.FA15 <sage.matroids.database_matroids.FA15>`
    - :func:`matroids.catalog.N4 <sage.matroids.database_matroids.N4>`

- Collection of various matroids (``matroids.catalog.`` + :kbd:`Tab`)
    - :func:`matroids.catalog.NonVamos <sage.matroids.database_matroids.NonVamos>`
    - :func:`matroids.catalog.TicTacToe <sage.matroids.database_matroids.TicTacToe>`
    - :func:`matroids.catalog.Q10 <sage.matroids.database_matroids.Q10>`
    - :func:`matroids.catalog.N1 <sage.matroids.database_matroids.N1>`
    - :func:`matroids.catalog.N2 <sage.matroids.database_matroids.N2>`
    - :func:`matroids.catalog.BetsyRoss <sage.matroids.database_matroids.BetsyRoss>`
    - :func:`matroids.catalog.Block_9_4 <sage.matroids.database_matroids.Block_9_4>`
    - :func:`matroids.catalog.Block_10_5 <sage.matroids.database_matroids.Block_10_5>`
    - :func:`matroids.catalog.ExtendedBinaryGolayCode <sage.matroids.database_matroids.ExtendedBinaryGolayCode>`
    - :func:`matroids.catalog.AG23minus <sage.matroids.database_matroids.AG23minus>`
    - :func:`matroids.catalog.NotP8 <sage.matroids.database_matroids.NotP8>`
    - :func:`matroids.catalog.D16 <sage.matroids.database_matroids.D16>`
    - :func:`matroids.catalog.Terrahawk <sage.matroids.database_matroids.Terrahawk>`
    - :func:`matroids.catalog.R9A <sage.matroids.database_matroids.R9A>`
    - :func:`matroids.catalog.R9B <sage.matroids.database_matroids.R9B>`
    - :func:`matroids.catalog.P9 <sage.matroids.database_matroids.P9>`
"""

# Do not add code to this file, only imports.
# Workaround for help in the notebook (needs parentheses in this file)
from sage.misc.lazy_import import lazy_import

# user-accessible:
lazy_import(
    "sage.matroids.database_collections",
    ("AllMatroids", "BrettellMatroids", "OxleyMatroids", "VariousMatroids"),
)
lazy_import(
    "sage.matroids.database_matroids",
    (
        "AG",
        "CompleteGraphic",
        "PG",
        "Psi",
        "Spike",
        "Theta",
        "Uniform",
        "Wheel",
        "Whirl",
        "Z",
    ),
)
lazy_import("sage.matroids", "catalog")
lazy_import("sage.matroids", "named_matroids")

del lazy_import  # delete so it doesn't appear under tab completion
