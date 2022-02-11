r"""
Catalog of permutation groups

Type ``groups.permutation.<tab>`` to access examples
of groups implemented as permutation groups.
"""

# groups imported here will be available
# via  groups.permutation.<tab>
#
# Do not use this file for code
#
# If you import a new group, then add an
# entry to the list in the module-level
# docstring of groups/groups_catalog.py

from .permgroup_named import KleinFourGroup as KleinFour
from .permgroup_named import QuaternionGroup as Quaternion
from .permgroup_named import SymmetricGroup as Symmetric
from .permgroup_named import AlternatingGroup as Alternating
from .permgroup_named import CyclicPermutationGroup as Cyclic
from .permgroup_named import ComplexReflectionGroup as ComplexReflection
from .permgroup_named import DihedralGroup as Dihedral
from .permgroup_named import DiCyclicGroup as DiCyclic
from .permgroup_named import MathieuGroup as Mathieu
from .permgroup_named import JankoGroup as Janko
from .permgroup_named import SuzukiSporadicGroup as SuzukiSporadic
from .permgroup_named import SuzukiGroup as Suzuki
from .permgroup_named import (PGL, PSL, PSp,PSU,PGU,)
from .permgroup_named import TransitiveGroup as Transitive
from .cubegroup import CubeGroup as RubiksCube
