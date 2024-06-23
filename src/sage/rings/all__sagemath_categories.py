# sage_setup: distribution = sagemath-categories

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

# Integers modulo n.
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing, Zmod
from sage.rings.finite_rings.integer_mod import IntegerMod, Mod, mod
Integers = IntegerModRing

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

# Function field
from sage.rings.function_field.all__sagemath_categories import *

# Double precision floating point numbers
from sage.rings.real_double import RealDoubleField, RDF, RealDoubleElement

# Lazy reals
from sage.rings.real_lazy import RealLazyField, RLF, ComplexLazyField, CLF

# Ideals
from sage.rings.ideal import Ideal
ideal = Ideal

# Semirings
from sage.rings.semirings.all import *

from sage.rings.finite_rings.all__sagemath_categories import *
from sage.rings.function_field.all__sagemath_categories import *
from sage.rings.number_field.all__sagemath_categories import *
from sage.rings.padics.all__sagemath_categories import *
from sage.rings.polynomial.all__sagemath_categories import *

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

# continued fractions
from sage.rings.continued_fraction import (continued_fraction,
                                           continued_fraction_list)

# Preliminary version of real numbers for doctesting without sage.rings.real_mpfr.
# sage.rings.all redefines it.
RealNumber = RR = RDF             # used by the preparser to wrap real literals

del lazy_import
