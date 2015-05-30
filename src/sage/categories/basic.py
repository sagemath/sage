r"""
A subset of sage.categories.all with just the basic categories needed
for sage startup (i.e. to define ZZ, QQ, ...).
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from objects import Objects
from sets_cat import Sets, EmptySetError
from posets import Posets
# For backward compatibility; will be deprecated at some point
PartiallyOrderedSets = Posets
OrderedSets = Posets

from additive_magmas import AdditiveMagmas
from commutative_additive_semigroups import CommutativeAdditiveSemigroups
from commutative_additive_monoids import CommutativeAdditiveMonoids
from commutative_additive_groups import CommutativeAdditiveGroups

from magmas import Magmas
from semigroups import Semigroups
from monoids import Monoids
from groups import Groups
from partially_ordered_monoids import PartiallyOrderedMonoids
# For backward compatibility; might be deprecated at some point
OrderedMonoids = PartiallyOrderedMonoids

from rngs import Rngs
from semirings import Semirings
from rings import Rings
from domains import Domains
from division_rings import DivisionRings

from commutative_rings import CommutativeRings
from integral_domains import IntegralDomains
from gcd_domains import GcdDomains
from principal_ideal_domains import PrincipalIdealDomains
from euclidean_domains import EuclideanDomains
from unique_factorization_domains import UniqueFactorizationDomains
from complete_discrete_valuation import CompleteDiscreteValuationRings

from fields import Fields
from quotient_fields import QuotientFields
from finite_fields import FiniteFields
from discrete_valuation import DiscreteValuationRings, DiscreteValuationFields
from complete_discrete_valuation import CompleteDiscreteValuationRings, CompleteDiscreteValuationFields
