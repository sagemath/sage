# sage_setup: distribution = sagemath-gap
from sage.groups.perm_gps.permgroup_named import (SymmetricGroup, AlternatingGroup,
                                                  DihedralGroup, SplitMetacyclicGroup,
                                                  SemidihedralGroup, CyclicPermutationGroup,
                                                  DiCyclicGroup, TransitiveGroup,
                                                  PGL, PSL, PSp, PSU, PGU,
                                                  MathieuGroup, KleinFourGroup, QuaternionGroup,
                                                  PrimitiveGroup, PrimitiveGroups,
                                                  SuzukiGroup, TransitiveGroups,
                                                  GeneralDihedralGroup, SmallPermutationGroup)

from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic, PermutationGroup_subgroup, direct_product_permgroups

from sage.groups.perm_gps.constructor import PermutationGroupElement

from sage.groups.perm_gps.permgroup_morphism import (PermutationGroupMorphism as PermutationGroupMap,
                                                     PermutationGroupMorphism_im_gens,
                                                     PermutationGroupMorphism_id)
PermutationGroupMorphism = PermutationGroupMorphism_im_gens

from sage.groups.perm_gps.cubegroup import CubeGroup, RubiksCube
