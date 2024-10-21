# sage_setup: distribution = sagemath-modules

from sage.rings.all__sagemath_categories import *

from sage.rings.function_field.all__sagemath_modules import *
from sage.rings.polynomial.all__sagemath_modules import *

# Real numbers
from sage.rings.real_mpfr import (RealField, RR,
                                  create_RealNumber as RealNumber)   # this is used by the preparser to wrap real literals -- very important.
Reals = RealField

# Complex numbers

from sage.rings.complex_mpfr import ComplexField
from sage.rings.complex_mpfr import create_ComplexNumber as ComplexNumber
Complexes = ComplexField

from sage.rings.complex_double import ComplexDoubleField, ComplexDoubleElement, CDF

from sage.rings.cc import CC

from sage.rings.complex_mpc import MPComplexField

# invariant theory
from sage.rings.invariants.all import *
