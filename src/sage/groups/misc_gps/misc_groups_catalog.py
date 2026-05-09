r"""
Type ``groups.misc.<tab>`` to access examples
of various groups not listed elsewhere.
"""

# groups imported here will be available
# via  groups.misc.<tab>
#
# Do not use this file for code
#
# If you import a new group, then add an
# entry to the list in the module-level
# docstring of groups/groups_catalog.py

from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup as AdditiveAbelian
from sage.groups.abelian_gps.abelian_group import AbelianGroup as MultiplicativeAbelian
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
from sage.groups.free_group import FreeGroup as Free
from sage.groups.artin import ArtinGroup as Artin
from sage.groups.braid import BraidGroup as Braid
from sage.groups.semimonomial_transformations.semimonomial_transformation_group import SemimonomialTransformationGroup as SemimonomialTransformation
from sage.combinat.root_system.coxeter_group import CoxeterGroup
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.combinat.colored_permutations import ShephardToddFamilyGroup as ShephardToddFamily
from sage.groups.raag import RightAngledArtinGroup as RightAngledArtin
from sage.combinat.root_system.reflection_group_real import ReflectionGroup
from sage.groups.cactus_group import CactusGroup as Cactus
from sage.groups.cactus_group import PureCactusGroup as PureCactus
