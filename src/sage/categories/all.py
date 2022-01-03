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

from . import primer

from sage.misc.lazy_import import lazy_import

from .category import Category

from .category_types import Elements

from .chain_complexes import ChainComplexes, HomologyFunctor

from .simplicial_complexes import SimplicialComplexes

from .tensor import tensor
from .signed_tensor import tensor_signed
from .cartesian_product import cartesian_product

from .functor  import (ForgetfulFunctor,
                      IdentityFunctor)

from .homset   import (Hom, hom,
                      End, end,
                      Homset, HomsetWithBase)

from .morphism import Morphism

from .basic import *

from .realizations import Realizations

from .g_sets import GSets
from .pointed_sets import PointedSets

from .sets_with_partial_maps import SetsWithPartialMaps
from .sets_with_grading import SetsWithGrading

from .groupoid import Groupoid
from .permutation_groups import PermutationGroups

# enumerated sets
from .finite_sets import FiniteSets
from .enumerated_sets import EnumeratedSets
from .finite_enumerated_sets import FiniteEnumeratedSets
from .infinite_enumerated_sets import InfiniteEnumeratedSets

# posets
from .posets import Posets
from .finite_posets import FinitePosets
from .lattice_posets import LatticePosets
from .finite_lattice_posets import FiniteLatticePosets

# finite groups/...
from .finite_semigroups import FiniteSemigroups
from .finite_monoids import FiniteMonoids
from .finite_groups import FiniteGroups
from .finite_permutation_groups import FinitePermutationGroups

# fields
from .number_fields import NumberFields
from .function_fields import FunctionFields

# modules
from .left_modules import LeftModules
from .right_modules import RightModules
from .bimodules import Bimodules

from .modules import Modules
RingModules = Modules
from .vector_spaces import VectorSpaces

# (hopf) algebra structures
from .algebras import Algebras
from .commutative_algebras import CommutativeAlgebras
from .coalgebras import Coalgebras
from .bialgebras import Bialgebras
from .hopf_algebras import HopfAlgebras
from .lie_algebras import LieAlgebras

# specific algebras
from .monoid_algebras import MonoidAlgebras
from .group_algebras import GroupAlgebras
from .matrix_algebras import MatrixAlgebras

# ideals
from .ring_ideals import RingIdeals
Ideals = RingIdeals
from .commutative_ring_ideals import CommutativeRingIdeals
from .algebra_modules import AlgebraModules
from .algebra_ideals import AlgebraIdeals
from .commutative_algebra_ideals import CommutativeAlgebraIdeals

# schemes and varieties
from .modular_abelian_varieties import ModularAbelianVarieties
from .schemes import Schemes

# * with basis
from .modules_with_basis import ModulesWithBasis
FreeModules = ModulesWithBasis
from .hecke_modules            import HeckeModules
from .algebras_with_basis      import AlgebrasWithBasis
from .coalgebras_with_basis    import CoalgebrasWithBasis
from .bialgebras_with_basis    import BialgebrasWithBasis
from .hopf_algebras_with_basis import HopfAlgebrasWithBasis

# finite dimensional * with basis
from .finite_dimensional_modules_with_basis       import FiniteDimensionalModulesWithBasis
from .finite_dimensional_algebras_with_basis      import FiniteDimensionalAlgebrasWithBasis
from .finite_dimensional_coalgebras_with_basis    import FiniteDimensionalCoalgebrasWithBasis
from .finite_dimensional_bialgebras_with_basis    import FiniteDimensionalBialgebrasWithBasis
from .finite_dimensional_hopf_algebras_with_basis import FiniteDimensionalHopfAlgebrasWithBasis

# graded *
from .graded_modules       import GradedModules
from .graded_algebras      import GradedAlgebras
from .graded_coalgebras    import GradedCoalgebras
from .graded_bialgebras    import GradedBialgebras
from .graded_hopf_algebras import GradedHopfAlgebras

# graded * with basis
from .graded_modules_with_basis       import GradedModulesWithBasis
from .graded_algebras_with_basis      import GradedAlgebrasWithBasis
from .graded_coalgebras_with_basis    import GradedCoalgebrasWithBasis
from .graded_bialgebras_with_basis    import GradedBialgebrasWithBasis
from .graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis

# Coxeter groups
from .coxeter_groups import CoxeterGroups
lazy_import('sage.categories.finite_coxeter_groups', 'FiniteCoxeterGroups')
from .weyl_groups import WeylGroups
from .finite_weyl_groups import FiniteWeylGroups
from .affine_weyl_groups import AffineWeylGroups

# crystal bases
from .crystals import Crystals
from .highest_weight_crystals import HighestWeightCrystals
from .regular_crystals import RegularCrystals
from .finite_crystals import FiniteCrystals
from .classical_crystals import ClassicalCrystals

# polyhedra
lazy_import('sage.categories.polyhedra', 'PolyhedralSets')

# lie conformal algebras
lazy_import('sage.categories.lie_conformal_algebras', 'LieConformalAlgebras')
