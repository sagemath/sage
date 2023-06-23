# sage_setup: distribution = sagemath-mpmath
import sage.libs.mpmath

# Patch mpmath to use Cythonized functions
from . import utils as _utils

# Also import internal functions
from sage.libs.mpmath._vendor.mpmath.libmp import *

# Main namespace
from sage.libs.mpmath import *

# Utilities
from .utils import call, mpmath_to_sage, sage_to_mpmath

# Use mpmath internal functions for constants, to avoid unnecessary overhead
_constants_funcs = {
  'glaisher': glaisher_fixed,
  'khinchin': khinchin_fixed,
  'twinprime': twinprime_fixed,
  'mertens': mertens_fixed
}

def eval_constant(name, ring):
    prec = ring.precision() + 20
    return ring(_constants_funcs[name](prec)) >> prec
