# sage_setup: distribution = sagemath-symbolics
from sage.symbolic.ring import SR
from sage.symbolic.constants import (pi, e, NaN, golden_ratio, log2, euler_gamma, catalan,
                                     khinchin, twinprime, mertens, glaisher)
from sage.symbolic.expression import Expression, solve_diophantine, hold
from sage.symbolic.callable import CallableSymbolicExpressionRing

from sage.symbolic.relation import solve, solve_mod, solve_ineq
from sage.symbolic.assumptions import assume, forget, assumptions, assuming

from sage.symbolic.units import units

π = pi

from sage.symbolic.operators import D
