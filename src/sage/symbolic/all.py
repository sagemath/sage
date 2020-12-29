from sage.misc.lazy_import import lazy_import

import sage.libs.pynac.pynac # initialize pynac before .ring
lazy_import("sage.libs.pynac.pynac", "I", deprecation=(18036,
        "import I from sage.symbolic.constants for the imaginary unit viewed as an element of SR, or from sage.rings.imaginary_unit for the element of ZZ[i]"))
lazy_import("sage.libs.pynac.pynac", "I", as_="i", deprecation=(18036,
        "import I from sage.symbolic.constants for the imaginary unit viewed as an element of SR, or from sage.rings.imaginary_unit for the element of ZZ[i]"))

from .ring import SR
from .constants import (pi, e, NaN, golden_ratio, log2, euler_gamma, catalan,
                       khinchin, twinprime, mertens, glaisher)
from .expression import Expression, solve_diophantine, hold
from .callable import CallableSymbolicExpressionRing

from sage.symbolic.relation import solve, solve_mod, solve_ineq
from sage.symbolic.assumptions import assume, forget, assumptions, assuming

from .units import units

π = pi
