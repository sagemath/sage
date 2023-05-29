# Ring base classes
from .ring import (Ring, Field, CommutativeRing, IntegralDomain,
    DedekindDomain, PrincipalIdealDomain, EuclideanDomain)

# Ring element base classes
from sage.structure.element import (CommutativeAlgebraElement,
        RingElement, CommutativeRingElement, IntegralDomainElement,
        DedekindDomainElement, PrincipalIdealDomainElement,
        EuclideanDomainElement, FieldElement)

# Rational numbers
from .rational import Rational
from .rational_field import RationalField, QQ
Rationals = RationalField

# Rational integers.
from .integer_ring import IntegerRing, ZZ, crt_basis
from .integer import Integer

# Integers modulo n.
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing, Zmod
from sage.rings.finite_rings.integer_mod import IntegerMod, Mod, mod
Integers = IntegerModRing

# Infinities
from .infinity import infinity, Infinity, InfinityRing, unsigned_infinity, UnsignedInfinityRing
oo = infinity

# Quotient
from .quotient_ring import QuotientRing

# Localization
from .localization import Localization

# Fraction field
from .fraction_field import FractionField
Frac = FractionField

# Function field
from .function_field.all__sagemath_categories import *

# Double precision floating point numbers
from .real_double import RealDoubleField, RDF, RealDoubleElement

# Ideals
from .ideal import Ideal
ideal = Ideal

# Semirings
from .semirings.all import *

from .finite_rings.all__sagemath_categories import *
from .number_field.all__sagemath_categories import *
from .padics.all__sagemath_categories import *
from .polynomial.all__sagemath_categories import *

# Power series rings
from .power_series_ring import PowerSeriesRing

from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.power_series_ring_element', 'PowerSeries', deprecation=33602)

# Laurent series ring in one variable
from .laurent_series_ring import LaurentSeriesRing
lazy_import('sage.rings.laurent_series_ring_element', 'LaurentSeries', deprecation=33602)

# Puiseux series ring
from .puiseux_series_ring import PuiseuxSeriesRing
lazy_import('sage.rings.puiseux_series_ring_element', 'PuiseuxSeries', deprecation=33602)

del lazy_import

# Big-oh notation
from .big_oh import O

# continued fractions
from sage.rings.continued_fraction import (continued_fraction,
                                           continued_fraction_list)

# Lazy reals
from .real_lazy import RealLazyField, RLF, ComplexLazyField, CLF

# Preliminary version of real numbers for doctesting without sage.rings.real_mpfr.
# sage.rings.all redefines it.
RealNumber = RR = RDF             # used by the preparser to wrap real literals
