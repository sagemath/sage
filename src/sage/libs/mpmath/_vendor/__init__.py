
from .mpmath import *
from .mpmath.libmp.backend import BACKEND as _BACKEND

assert _BACKEND == 'sage'

del _BACKEND
