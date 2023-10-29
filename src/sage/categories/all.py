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


from sage.categories.all__sagemath_objects import *
from sage.categories.basic import *

# enumerated sets
# posets
# finite groups/...
# fields
# modules
from sage.categories.modules import Modules
from sage.misc.lazy_import import lazy_import

RingModules = Modules

# (hopf) algebra structures

# specific algebras

# ideals
from sage.categories.ring_ideals import RingIdeals

Ideals = RingIdeals

# schemes and varieties

# * with basis
from sage.categories.modules_with_basis import ModulesWithBasis

FreeModules = ModulesWithBasis

# finite dimensional * with basis

# graded *

# graded * with basis

# Coxeter groups
lazy_import('sage.categories.finite_coxeter_groups', 'FiniteCoxeterGroups')

# crystal bases

# polyhedra
lazy_import('sage.categories.polyhedra', 'PolyhedralSets')

# lie conformal algebras
lazy_import('sage.categories.lie_conformal_algebras', 'LieConformalAlgebras')
