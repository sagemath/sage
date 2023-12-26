r"""
Catalog of matroids

A module containing constructors for several common matroids.

A list of all matroids in this module is available via tab completion.
Let ``<tab>`` indicate pressing the :kbd:`Tab` key. So begin by typing
``matroids.<tab>`` to see the various constructions available. Many special
matroids can be accessed from the submenu ``matroids.named_matroids.<tab>``.

To create a custom matroid using a variety of inputs, see the function
:func:`Matroid() <sage.matroids.constructor.Matroid>`.

- Parametrized matroid constructors
    - :func:`matroids.AG <sage.matroids.catalog.AG>`
    - :func:`matroids.CompleteGraphic <sage.matroids.catalog.CompleteGraphic>`
    - :func:`matroids.PG <sage.matroids.catalog.PG>`
    - :func:`matroids.Uniform <sage.matroids.catalog.Uniform>`
    - :func:`matroids.Wheel <sage.matroids.catalog.Wheel>`
    - :func:`matroids.Whirl <sage.matroids.catalog.Whirl>`

- Named matroids (``matroids.named_matroids.<tab>``)
    - :func:`matroids.named_matroids.AG23minus <sage.matroids.catalog.AG23minus>`
    - :func:`matroids.named_matroids.AG32prime <sage.matroids.catalog.AG32prime>`
    - :func:`matroids.named_matroids.BetsyRoss <sage.matroids.catalog.BetsyRoss>`
    - :func:`matroids.named_matroids.Block_9_4 <sage.matroids.catalog.Block_9_4>`
    - :func:`matroids.named_matroids.Block_10_5 <sage.matroids.catalog.Block_10_5>`
    - :func:`matroids.named_matroids.D16 <sage.matroids.catalog.D16>`
    - :func:`matroids.named_matroids.ExtendedBinaryGolayCode <sage.matroids.catalog.ExtendedBinaryGolayCode>`
    - :func:`matroids.named_matroids.ExtendedTernaryGolayCode <sage.matroids.catalog.ExtendedTernaryGolayCode>`
    - :func:`matroids.named_matroids.F8 <sage.matroids.catalog.F8>`
    - :func:`matroids.named_matroids.Fano <sage.matroids.catalog.Fano>`
    - :func:`matroids.named_matroids.J <sage.matroids.catalog.J>`
    - :func:`matroids.named_matroids.K33dual <sage.matroids.catalog.K33dual>`
    - :func:`matroids.named_matroids.L8 <sage.matroids.catalog.L8>`
    - :func:`matroids.named_matroids.N1 <sage.matroids.catalog.N1>`
    - :func:`matroids.named_matroids.N2 <sage.matroids.catalog.N2>`
    - :func:`matroids.named_matroids.NonFano <sage.matroids.catalog.NonFano>`
    - :func:`matroids.named_matroids.NonPappus <sage.matroids.catalog.NonPappus>`
    - :func:`matroids.named_matroids.NonVamos <sage.matroids.catalog.NonVamos>`
    - :func:`matroids.named_matroids.NotP8 <sage.matroids.catalog.NotP8>`
    - :func:`matroids.named_matroids.O7 <sage.matroids.catalog.O7>`
    - :func:`matroids.named_matroids.P6 <sage.matroids.catalog.P6>`
    - :func:`matroids.named_matroids.P7 <sage.matroids.catalog.P7>`
    - :func:`matroids.named_matroids.P8 <sage.matroids.catalog.P8>`
    - :func:`matroids.named_matroids.P8pp <sage.matroids.catalog.P8pp>`
    - :func:`matroids.named_matroids.P9 <sage.matroids.catalog.P9>`
    - :func:`matroids.named_matroids.Pappus <sage.matroids.catalog.Pappus>`
    - :func:`matroids.named_matroids.Q6 <sage.matroids.catalog.Q6>`
    - :func:`matroids.named_matroids.Q8 <sage.matroids.catalog.Q8>`
    - :func:`matroids.named_matroids.Q10 <sage.matroids.catalog.Q10>`
    - :func:`matroids.named_matroids.R6 <sage.matroids.catalog.R6>`
    - :func:`matroids.named_matroids.R8 <sage.matroids.catalog.R8>`
    - :func:`matroids.named_matroids.R9A <sage.matroids.catalog.R9A>`
    - :func:`matroids.named_matroids.R9B <sage.matroids.catalog.R9B>`
    - :func:`matroids.named_matroids.R10 <sage.matroids.catalog.R10>`
    - :func:`matroids.named_matroids.R12 <sage.matroids.catalog.R12>`
    - :func:`matroids.named_matroids.S8 <sage.matroids.catalog.S8>`
    - :func:`matroids.named_matroids.T8 <sage.matroids.catalog.T8>`
    - :func:`matroids.named_matroids.T12 <sage.matroids.catalog.T12>`
    - :func:`matroids.named_matroids.TernaryDowling3 <sage.matroids.catalog.TernaryDowling3>`
    - :func:`matroids.named_matroids.Terrahawk <sage.matroids.catalog.Terrahawk>`
    - :func:`matroids.named_matroids.TicTacToe <sage.matroids.catalog.TicTacToe>`
    - :func:`matroids.named_matroids.Vamos <sage.matroids.catalog.Vamos>`

"""

# Do not add code to this file, only imports.
# Workaround for help in the notebook (needs parentheses in this file)

# user-accessible:
from sage.matroids.catalog import AG, CompleteGraphic, PG, Uniform, Wheel, Whirl
from sage.matroids import named_matroids
