from .all__sagemath_categories import *

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

# Infinities
from .infinity import infinity, Infinity, InfinityRing, unsigned_infinity, UnsignedInfinityRing

# Quotient
from .quotient_ring import QuotientRing

# Localization
from .localization import Localization

# Double precision floating point numbers
from .real_double import RealDoubleField, RDF, RealDoubleElement


# Preliminary version of real numbers for doctesting without sage.rings.real_mpfr.
# sage.rings.all redefines it.

RealNumber = RDF             # used by the preparser to wrap real literals
