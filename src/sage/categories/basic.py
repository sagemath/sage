# sage_setup: distribution = sagemath-objects
r"""
A subset of sage.categories.all with just the basic categories needed
for sage startup (i.e. to define ZZ, QQ, ...).
"""
# *****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************

from sage.categories.objects import Objects
from sage.categories.sets_cat import Sets, EmptySetError
from sage.categories.posets import Posets

# For backward compatibility; will be deprecated at some point
PartiallyOrderedSets = Posets
OrderedSets = Posets

from sage.categories.additive_magmas import AdditiveMagmas
from sage.categories.commutative_additive_semigroups import (
    CommutativeAdditiveSemigroups,
)
from sage.categories.commutative_additive_monoids import CommutativeAdditiveMonoids
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups

from sage.categories.magmas import Magmas
from sage.categories.semigroups import Semigroups
from sage.categories.monoids import Monoids
from sage.categories.groups import Groups
from sage.categories.partially_ordered_monoids import PartiallyOrderedMonoids

# For backward compatibility; might be deprecated at some point
OrderedMonoids = PartiallyOrderedMonoids

from sage.categories.rngs import Rngs
from sage.categories.semirings import Semirings
from sage.categories.rings import Rings
from sage.categories.domains import Domains
from sage.categories.division_rings import DivisionRings

from sage.categories.commutative_rings import CommutativeRings
from sage.categories.integral_domains import IntegralDomains
from sage.categories.gcd_domains import GcdDomains
from sage.categories.dedekind_domains import DedekindDomains
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.categories.euclidean_domains import EuclideanDomains
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains

from sage.categories.fields import Fields
from sage.categories.quotient_fields import QuotientFields
from sage.categories.finite_fields import FiniteFields
from sage.categories.discrete_valuation import (
    DiscreteValuationRings,
    DiscreteValuationFields,
)
from sage.categories.complete_discrete_valuation import (
    CompleteDiscreteValuationRings,
    CompleteDiscreteValuationFields,
)
