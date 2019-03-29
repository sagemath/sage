"""
Test for deprecations of imports into global namespace::

    sage: backtrack_all
    doctest:warning...:
    DeprecationWarning:
    Importing backtrack_all from here is deprecated. If you need to use it, please import it directly from sage.games.sudoku_backtrack
    See https://trac.sagemath.org/27066 for details.
    ...
"""
from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import

lazy_import("sage.games.sudoku_backtrack", 'backtrack_all', deprecation=27066)

from .sudoku import Sudoku, sudoku
from .hexad import Minimog

del absolute_import
