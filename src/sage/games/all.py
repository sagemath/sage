# sage_setup: distribution = sagemath-combinat

from sage.misc.lazy_import import lazy_import

lazy_import('sage.games.sudoku', ['Sudoku', 'sudoku'])
lazy_import('sage.games.hexad', ['Minimog'])
del lazy_import
