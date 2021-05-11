r"""
Root system features that are imported by default in the interpreter namespace
"""

from sage.misc.lazy_import import lazy_import

from .cartan_type import CartanType
from .dynkin_diagram import DynkinDiagram
from .cartan_matrix import CartanMatrix
from .coxeter_matrix import CoxeterMatrix
from .coxeter_type import CoxeterType
from .root_system import RootSystem, WeylDim
lazy_import('sage.combinat.root_system.weyl_group', ['WeylGroup',
                                                     'WeylGroupElement'])
lazy_import('sage.combinat.root_system.reflection_group_real',
            'ReflectionGroup')
lazy_import('sage.combinat.root_system.extended_affine_weyl_group',
            'ExtendedAffineWeylGroup')
lazy_import('sage.combinat.root_system.coxeter_group', 'CoxeterGroup')
lazy_import('sage.combinat.root_system.weyl_characters', ['WeylCharacterRing',
                                                          'WeightRing'])
lazy_import('sage.combinat.root_system.fusion_ring', ['FusionRing'])
from .branching_rules import BranchingRule, branching_rule_from_plethysm, branching_rule

lazy_import('sage.combinat.root_system.non_symmetric_macdonald_polynomials', 'NonSymmetricMacdonaldPolynomials')
lazy_import('sage.combinat.root_system.integrable_representations', 'IntegrableRepresentation')
