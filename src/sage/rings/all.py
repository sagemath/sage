"""
Rings

Tests for deprecations of imports in global namespace from :trac:`33602`::

    sage: PowerSeries
    doctest:warning...:
    DeprecationWarning:
    Importing PowerSeries from here is deprecated;
    please use "from sage.rings.power_series_ring_element import PowerSeries" instead.
    See https://github.com/sagemath/sage/issues/33602 for details.
    ...
    sage: PuiseuxSeries
    doctest:warning...:
    DeprecationWarning:
    Importing PuiseuxSeries from here is deprecated;
    please use "from sage.rings.puiseux_series_ring_element import PuiseuxSeries" instead.
    See https://github.com/sagemath/sage/issues/33602 for details.
    ...
    sage: LaurentSeries
    doctest:warning...:
    DeprecationWarning:
    Importing LaurentSeries from here is deprecated;
    please use "from sage.rings.laurent_series_ring_element import LaurentSeries" instead.
    See https://github.com/sagemath/sage/issues/33602 for details.
    ...
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

from .all__sagemath_polyhedra import *

# Finite fields
from .finite_rings.all import *

# Number field
from .number_field.all import *

# Finite residue fields
from .finite_rings.residue_field import ResidueField

# p-adic field
from .padics.all import *
from .padics.padic_printing import _printer_defaults as padic_printing

# valuations
from .valuation.all import *

# Real numbers
from .real_mpfr import (RealField, RR,
                       create_RealNumber as RealNumber)   # this is used by the preparser to wrap real literals -- very important.
Reals = RealField

from .real_lazy import RealLazyField, RLF, ComplexLazyField, CLF

from sage.rings.real_arb import RealBallField, RBF

# Polynomial Rings and Polynomial Quotient Rings
from .polynomial.all import *


# Algebraic numbers
from .qqbar import (AlgebraicRealField, AA,
                   AlgebraicReal,
                   AlgebraicField, QQbar,
                   AlgebraicNumber,
                   number_field_elements_from_algebraics)
from .universal_cyclotomic_field import UniversalCyclotomicField, E

# Intervals
from .real_mpfi import (RealIntervalField,
                       RIF,
                       RealInterval)

# Complex numbers
from .complex_mpfr import ComplexField
from .complex_mpfr import create_ComplexNumber as ComplexNumber
Complexes = ComplexField
from .complex_interval_field import ComplexIntervalField
from .complex_interval import (create_ComplexIntervalFieldElement as ComplexIntervalFieldElement)

from .complex_double import ComplexDoubleField, ComplexDoubleElement, CDF

from .complex_mpc import MPComplexField

from sage.rings.complex_arb import ComplexBallField, CBF

lazy_import("sage.rings.imaginary_unit", "I")

# Tate algebras
from .tate_algebra import TateAlgebra

# Puiseux series ring
from .puiseux_series_ring import PuiseuxSeriesRing
lazy_import('sage.rings.puiseux_series_ring_element', 'PuiseuxSeries', deprecation=33602)

# Pseudo-ring of PARI objects.
from .pari_ring import PariRing, Pari

# c-finite sequences
from .cfinite_sequence import CFiniteSequence, CFiniteSequences

from .bernoulli_mod_p import bernoulli_mod_p, bernoulli_mod_p_single

from .monomials import monomials

from .cc import CC
from .cif import CIF

# invariant theory
from .invariants.all import *

from .fast_arith import prime_range

# continued fractions
from sage.rings.continued_fraction import (continued_fraction,
                                           continued_fraction_list)

# asymptotic ring
from .asymptotic.all import *

# Register classes in numbers abc
from . import numbers_abc
