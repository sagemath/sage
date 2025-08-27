# sage_setup: distribution = sagemath-categories
r"""
Sage categories quickref

- ``sage.categories.primer?``                      a primer on Elements, Parents, and Categories
- ``sage.categories.tutorial?``                    a tutorial on Elements, Parents, and Categories
- ``Category?``                                    technical background on categories
- ``Sets()``, ``Semigroups()``, ``Algebras(QQ)``   some categories
- ``SemiGroups().example()??``                     sample implementation of a semigroup
- ``Hom(A, B)``, ``End(A, Algebras())``            homomorphisms sets
- ``tensor``, ``cartesian_product``                functorial constructions

Module layout:

- :mod:`sage.categories.basic`                the basic categories
- :mod:`sage.categories.all`                  all categories
- :mod:`sage.categories.semigroups`           the ``Semigroups()`` category
- :mod:`sage.categories.examples.semigroups`  the example of ``Semigroups()``
- :mod:`sage.categories.homset`               morphisms, ...
- :mod:`sage.categories.map`
- :mod:`sage.categories.morphism`
- :mod:`sage.categories.functors`
- :mod:`sage.categories.cartesian_product`    functorial constructions
- :mod:`sage.categories.tensor`
- :mod:`sage.categories.dual`
"""
# install the docstring of this module to the containing package
from sage.misc.namespace_package import install_doc
install_doc(__package__, __doc__)

from sage.categories import primer

from sage.misc.lazy_import import lazy_import

from sage.categories.all__sagemath_objects import *

from sage.categories.basic import *

from sage.categories.chain_complexes import ChainComplexes, HomologyFunctor

from sage.categories.simplicial_complexes import SimplicialComplexes

from sage.categories.tensor import tensor
from sage.categories.signed_tensor import tensor_signed

from sage.categories.g_sets import GSets
from sage.categories.pointed_sets import PointedSets

from sage.categories.sets_with_grading import SetsWithGrading

from sage.categories.groupoid import Groupoid
from sage.categories.permutation_groups import PermutationGroups

# enumerated sets
from sage.categories.finite_sets import FiniteSets
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets

# posets
from sage.categories.posets import Posets
from sage.categories.finite_posets import FinitePosets
from sage.categories.lattice_posets import LatticePosets
from sage.categories.finite_lattice_posets import FiniteLatticePosets

# finite groups/...
from sage.categories.finite_semigroups import FiniteSemigroups
from sage.categories.finite_monoids import FiniteMonoids
from sage.categories.finite_groups import FiniteGroups
from sage.categories.finite_permutation_groups import FinitePermutationGroups

# fields
from sage.categories.number_fields import NumberFields
from sage.categories.function_fields import FunctionFields

# modules
from sage.categories.left_modules import LeftModules
from sage.categories.right_modules import RightModules
from sage.categories.bimodules import Bimodules

from sage.categories.modules import Modules
RingModules = Modules
from sage.categories.vector_spaces import VectorSpaces

# (Hopf) algebra structures
from sage.categories.algebras import Algebras
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.coalgebras import Coalgebras
from sage.categories.bialgebras import Bialgebras
from sage.categories.hopf_algebras import HopfAlgebras
from sage.categories.lie_algebras import LieAlgebras

# operads
from .set_operads import SetOperads
from .operads import Operads

# specific algebras
from sage.categories.monoid_algebras import MonoidAlgebras
from sage.categories.group_algebras import GroupAlgebras
from sage.categories.matrix_algebras import MatrixAlgebras

# ideals
from sage.categories.ring_ideals import RingIdeals
Ideals = RingIdeals
from sage.categories.commutative_ring_ideals import CommutativeRingIdeals
from sage.categories.algebra_modules import AlgebraModules
from sage.categories.algebra_ideals import AlgebraIdeals
from sage.categories.commutative_algebra_ideals import CommutativeAlgebraIdeals

# schemes and varieties
from sage.categories.modular_abelian_varieties import ModularAbelianVarieties
from sage.categories.schemes import Schemes, AbelianVarieties, Jacobians

# * with basis
from sage.categories.modules_with_basis import ModulesWithBasis
FreeModules = ModulesWithBasis
from sage.categories.hecke_modules import HeckeModules
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.coalgebras_with_basis import CoalgebrasWithBasis
from sage.categories.bialgebras_with_basis import BialgebrasWithBasis
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from .operads_with_basis import OperadsWithBasis

# finite dimensional * with basis
from sage.categories.finite_dimensional_modules_with_basis import FiniteDimensionalModulesWithBasis
from sage.categories.finite_dimensional_algebras_with_basis import FiniteDimensionalAlgebrasWithBasis
from sage.categories.finite_dimensional_coalgebras_with_basis import FiniteDimensionalCoalgebrasWithBasis
from sage.categories.finite_dimensional_bialgebras_with_basis import FiniteDimensionalBialgebrasWithBasis
from sage.categories.finite_dimensional_hopf_algebras_with_basis import FiniteDimensionalHopfAlgebrasWithBasis

# graded *
from sage.categories.graded_modules import GradedModules
from sage.categories.graded_algebras import GradedAlgebras
from sage.categories.graded_coalgebras import GradedCoalgebras
from sage.categories.graded_bialgebras import GradedBialgebras
from sage.categories.graded_hopf_algebras import GradedHopfAlgebras

# graded * with basis
from sage.categories.graded_modules_with_basis import GradedModulesWithBasis
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.graded_coalgebras_with_basis import GradedCoalgebrasWithBasis
from sage.categories.graded_bialgebras_with_basis import GradedBialgebrasWithBasis
from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis

# Coxeter groups
from sage.categories.coxeter_groups import CoxeterGroups
lazy_import('sage.categories.finite_coxeter_groups', 'FiniteCoxeterGroups')
from sage.categories.weyl_groups import WeylGroups
from sage.categories.finite_weyl_groups import FiniteWeylGroups
from sage.categories.affine_weyl_groups import AffineWeylGroups

# crystal bases
from sage.categories.crystals import Crystals
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.classical_crystals import ClassicalCrystals

# polyhedra
lazy_import('sage.categories.polyhedra', 'PolyhedralSets')

# lie conformal algebras
lazy_import('sage.categories.lie_conformal_algebras', 'LieConformalAlgebras')
del lazy_import
del install_doc
