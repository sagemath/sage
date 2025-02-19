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

# Ring base classes
from sage.rings.ring import (Ring, Field, CommutativeRing, IntegralDomain,
                             PrincipalIdealDomain)

lazy_import("sage.rings.ring", "DedekindDomain")

# Ring element base classes
from sage.structure.element import (CommutativeAlgebraElement,
        RingElement, CommutativeRingElement, IntegralDomainElement,
        DedekindDomainElement, PrincipalIdealDomainElement,
        EuclideanDomainElement, FieldElement)

# Ideals
from sage.rings.ideal import Ideal
ideal = Ideal

# Quotient
from sage.rings.quotient_ring import QuotientRing

# Infinities
from sage.rings.infinity import infinity, Infinity, InfinityRing, unsigned_infinity, UnsignedInfinityRing

# Rational integers.
from sage.rings.integer_ring import IntegerRing, ZZ, crt_basis
from sage.rings.integer import Integer

# Rational numbers
from sage.rings.rational_field import RationalField, QQ
from sage.rings.rational import Rational
Rationals = RationalField

# Integers modulo n.
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing, Zmod
from sage.rings.finite_rings.integer_mod import IntegerMod, Mod, mod
Integers = IntegerModRing

# Finite fields
from sage.rings.finite_rings.all import *

# Number field
from sage.rings.number_field.all import *

# Function field
from sage.rings.function_field.all import *

# Finite residue fields
from sage.rings.finite_rings.residue_field import ResidueField

# p-adic field
from sage.rings.padics.all import *
from sage.rings.padics.padic_printing import _printer_defaults as padic_printing

# valuations
from sage.rings.valuation.all import *

# Semirings
from sage.rings.semirings.all import *

# Real numbers
from sage.rings.real_mpfr import (RealField, RR,
                       create_RealNumber as RealNumber)   # this is used by the preparser to wrap real literals -- very important.
Reals = RealField

from sage.rings.real_double import RealDoubleField, RDF, RealDoubleElement

from sage.rings.real_lazy import RealLazyField, RLF, ComplexLazyField, CLF

from sage.rings.real_arb import RealBallField, RBF

# Polynomial Rings and Polynomial Quotient Rings
from sage.rings.polynomial.all import *


# Algebraic numbers
from sage.rings.qqbar import (AlgebraicRealField, AA,
                   AlgebraicReal,
                   AlgebraicField, QQbar,
                   AlgebraicNumber,
                   number_field_elements_from_algebraics)
from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField, E

# Intervals
from sage.rings.real_mpfi import (RealIntervalField,
                       RIF,
                       RealInterval)

# Complex numbers
from sage.rings.complex_mpfr import ComplexField
from sage.rings.complex_mpfr import create_ComplexNumber as ComplexNumber
Complexes = ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.complex_interval import (create_ComplexIntervalFieldElement as ComplexIntervalFieldElement)

from sage.rings.complex_double import ComplexDoubleField, ComplexDoubleElement, CDF

from sage.rings.complex_mpc import MPComplexField

from sage.rings.complex_arb import ComplexBallField, CBF

lazy_import("sage.rings.imaginary_unit", "I")

# Power series rings
from sage.rings.power_series_ring import PowerSeriesRing

# Laurent series ring in one variable
from sage.rings.laurent_series_ring import LaurentSeriesRing

# Lazy Laurent series ring
lazy_import('sage.rings.lazy_series_ring', ['LazyLaurentSeriesRing', 'LazyPowerSeriesRing',
                                            'LazySymmetricFunctions', 'LazyDirichletSeriesRing'])

# Tate algebras
from sage.rings.tate_algebra import TateAlgebra

# Puiseux series ring
from sage.rings.puiseux_series_ring import PuiseuxSeriesRing

# Pseudo-ring of PARI objects.
from sage.rings.pari_ring import PariRing, Pari

# Big-oh notation
from sage.rings.big_oh import O

# Fraction field
from sage.rings.fraction_field import FractionField
Frac = FractionField

# Localization
from sage.rings.localization import Localization

# c-finite sequences
from sage.rings.cfinite_sequence import CFiniteSequence, CFiniteSequences

from sage.rings.bernoulli_mod_p import bernoulli_mod_p, bernoulli_mod_p_single

from sage.rings.monomials import monomials

from sage.rings.cc import CC
from sage.rings.cif import CIF

# invariant theory
from sage.rings.invariants.all import *

from sage.rings.fast_arith import prime_range

# continued fractions
from sage.rings.continued_fraction import (continued_fraction,
                                           continued_fraction_list)

# asymptotic ring
from sage.rings.asymptotic.all import *

# Register classes in numbers abc
from sage.rings import numbers_abc
