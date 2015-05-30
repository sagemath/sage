r"""
Root system features that are imported by default in the interpreter namespace
"""
from sage.misc.lazy_import import lazy_import
lazy_import('sage.combinat.root_system.associahedron', 'Associahedron')

from cartan_type import CartanType
from dynkin_diagram import DynkinDiagram
from cartan_matrix import CartanMatrix, cartan_matrix
from coxeter_matrix import coxeter_matrix
from root_system import RootSystem, WeylDim
from weyl_group import WeylGroup, WeylGroupElement
from coxeter_group import CoxeterGroup
from weyl_characters import WeylCharacterRing, WeightRing
from branching_rules import BranchingRule, branching_rule_from_plethysm, branching_rule
lazy_import('sage.combinat.root_system.non_symmetric_macdonald_polynomials', 'NonSymmetricMacdonaldPolynomials')
lazy_import('sage.combinat.root_system.integrable_representations', 'IntegrableRepresentation')

