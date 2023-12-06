"""
Rings
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.lazy_import import lazy_import

from sage.rings.all__sagemath_categories import *

# Finite fields
from sage.rings.finite_rings.all import *

from sage.rings.all__sagemath_combinat import *
from sage.rings.all__sagemath_flint import *
from sage.rings.all__sagemath_gap import *
from sage.rings.all__sagemath_modules import *

try:
    from sage.rings.all__sagemath_symbolics import *
except ImportError:
    pass

# Function field
from sage.rings.function_field.all import *

# Polynomial Rings and Polynomial Quotient Rings
from sage.rings.polynomial.all import *

# c-finite sequences
from sage.rings.cfinite_sequence import CFiniteSequence, CFiniteSequences

from sage.rings.fast_arith import prime_range

# asymptotic ring
# from sage.rings.asymptotic.all import *
lazy_import('sage.rings.asymptotic.asymptotic_ring', 'AsymptoticRing')
lazy_import('sage.rings.asymptotic.asymptotic_expansion_generators',
            'asymptotic_expansions')

# Register classes in numbers abc
from sage.rings import numbers_abc
del lazy_import
