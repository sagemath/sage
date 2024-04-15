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

# Ideals
from sage.rings.ideal import Ideal
ideal = Ideal

# To be added in #36566:

# from sage.rings.finite_rings.all__sagemath_categories import *
# from sage.rings.function_field.all__sagemath_categories import *
# from sage.rings.number_field.all__sagemath_categories import *
# from sage.rings.padics.all__sagemath_categories import *
# from sage.rings.polynomial.all__sagemath_categories import *

del lazy_import
