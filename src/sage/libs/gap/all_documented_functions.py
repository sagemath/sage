"""nodoctest
All Documented GAP Functions

This Python module contains all documented GAP functions, they can be
thought of as the official API of GAP.

EXAMPLES::

    sage: from sage.libs.gap.all_documented_functions import *
    sage: DihedralGroup(8)
    <pc group of size 8 with 3 generators>
    sage: GeneratorsOfGroup(_)
    [ f1, f2, f3 ]
    sage: List(_, Order)
    [ 2, 4, 2 ]
"""

from sage.libs.gap.libgap import libgap
from sage.libs.gap.assigned_names import FUNCTIONS as _FUNCTIONS



for _f in _FUNCTIONS:
    globals()[_f] = libgap.function_factory(_f)
