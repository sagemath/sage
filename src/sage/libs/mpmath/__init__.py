# sage_setup: distribution = sagemath-mpmath

from ._vendor.mpmath import *
from ._vendor.mpmath.libmp.backend import BACKEND as _BACKEND

assert _BACKEND == 'sage'

del _BACKEND
