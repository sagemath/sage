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

# Following will go to all__sagemath_categories.py in #36566

# Quotient
from sage.rings.quotient_ring import QuotientRing

# Infinities
from sage.rings.infinity import infinity, Infinity, InfinityRing, unsigned_infinity, UnsignedInfinityRing
oo = infinity

# Rational integers.
from sage.rings.integer_ring import IntegerRing, ZZ, crt_basis
from sage.rings.integer import Integer

# Rational numbers
from sage.rings.rational import Rational
from sage.rings.rational_field import RationalField, QQ
Rationals = RationalField

# Integers modulo n.
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing, Zmod
from sage.rings.finite_rings.integer_mod import IntegerMod, Mod, mod
Integers = IntegerModRing

# up to here (#36566)

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

# Following will go to all__sagemath_categories.py in #36566

# Semirings
from sage.rings.semirings.all import *

# Double precision floating point numbers
from sage.rings.real_double import RealDoubleField, RDF, RealDoubleElement

# Lazy reals
from sage.rings.real_lazy import RealLazyField, RLF, ComplexLazyField, CLF

# up to here (#36566)

# Polynomial Rings and Polynomial Quotient Rings
from sage.rings.polynomial.all import *

# Following will go to all__sagemath_categories.py in #36566

# Power series rings
from sage.rings.power_series_ring import PowerSeriesRing

# Laurent series ring in one variable
from sage.rings.laurent_series_ring import LaurentSeriesRing

# Puiseux series ring
from sage.rings.puiseux_series_ring import PuiseuxSeriesRing

# Big-oh notation
from sage.rings.big_oh import O

# Fraction field
from sage.rings.fraction_field import FractionField
Frac = FractionField

# Localization
from sage.rings.localization import Localization

# up to here (#36566)

# c-finite sequences
from sage.rings.cfinite_sequence import CFiniteSequence, CFiniteSequences

from sage.rings.fast_arith import prime_range

# Following will go to all__sagemath_categories.py in #36566

# continued fractions
from sage.rings.continued_fraction import (continued_fraction,
                                           continued_fraction_list)

# up to here (#36566)

# asymptotic ring
# from sage.rings.asymptotic.all import *
lazy_import('sage.rings.asymptotic.asymptotic_ring', 'AsymptoticRing')
lazy_import('sage.rings.asymptotic.asymptotic_expansion_generators',
            'asymptotic_expansions')

# Register classes in numbers abc
from sage.rings import numbers_abc
del lazy_import
