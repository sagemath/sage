# sage_setup: distribution = sagemath-modules
r"""
Library of Interesting Groups

Type ``groups.matrix.<tab>`` to access examples
of groups implemented as permutation groups.
"""

# groups imported here will be available
# via  groups.matrix.<tab>
#
# Do not use this file for code
#
# If you import a new group, then add an
# entry to the list in the module-level
# docstring of groups/groups_catalog.py

from .all__sagemath_modules import GL, SL, Sp, SU, GU, SO, GO
from .all__sagemath_modules import QuaternionMatrixGroupGF3 as QuaternionGF3

from sage.misc.lazy_import import lazy_import

lazy_import('sage.groups.matrix_gps.binary_dihedral', 'BinaryDihedralGroup', as_='BinaryDihedral')
lazy_import('sage.groups.matrix_gps.heisenberg', 'HeisenbergGroup', as_='Heisenberg')

del lazy_import
