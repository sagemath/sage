# sage_setup: distribution = sagemath-flint

from .rings.all__sagemath_flint import *

from sage.rings.qqbar import _init_qqbar
_init_qqbar()
