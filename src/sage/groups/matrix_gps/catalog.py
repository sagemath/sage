r"""
Library of Interesting Groups

Type ``groups.matrix.<tab>`` to access examples
of groups implemented as permutation groups.

TESTS::

    sage: 'absolute_import' in dir(groups.matrix)
    False
"""
from __future__ import absolute_import

# groups imported here will be available
# via  groups.matrix.<tab>
#
# Do not use this file for code
#
# If you import a new group, then add an
# entry to the list in the module-level
# docstring of groups/groups_catalog.py

from .all import GL, SL, Sp, SU, GU, SO, GO
from .all import QuaternionMatrixGroupGF3 as QuaternionGF3
from sage.groups.matrix_gps.binary_dihedral import BinaryDihedralGroup as BinaryDihedral
from sage.groups.matrix_gps.heisenberg import HeisenbergGroup as Heisenberg

# We don't want this to appear tab completion
del absolute_import
