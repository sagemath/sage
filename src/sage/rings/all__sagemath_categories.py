# Ring base classes
from sage.rings.ring import (Ring, Field, CommutativeRing, IntegralDomain,
                             DedekindDomain, PrincipalIdealDomain, EuclideanDomain)

# Ring element base classes
from sage.structure.element import (CommutativeAlgebraElement,
                                    RingElement, CommutativeRingElement, IntegralDomainElement,
                                    DedekindDomainElement, PrincipalIdealDomainElement,
                                    EuclideanDomainElement, FieldElement)

# Ideals
from sage.rings.ideal import Ideal
ideal = Ideal

# To be added in #36566:

# from sage.rings.finite_rings.all__sagemath_categories import *
# from sage.rings.function_field.all__sagemath_categories import *
# from sage.rings.number_field.all__sagemath_categories import *
# from sage.rings.padics.all__sagemath_categories import *
# from sage.rings.polynomial.all__sagemath_categories import *
