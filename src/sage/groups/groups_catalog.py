# sage_setup: distribution = sagemath-categories
r"""
Examples of Groups

The ``groups`` object may be used to access examples of various groups.
Using tab-completion on this object is an easy way to discover and quickly
create the groups that are available (as listed here).

Let ``<tab>`` indicate pressing the :kbd:`Tab` key.  So begin by typing
``groups.<tab>`` to the see primary divisions, followed by (for example)
``groups.matrix.<tab>`` to access various groups implemented as sets of matrices.

- Permutation Groups  (``groups.permutation.<tab>``)

  - :class:`groups.permutation.Symmetric <sage.groups.perm_gps.permgroup_named.SymmetricGroup>`
  - :class:`groups.permutation.Alternating <sage.groups.perm_gps.permgroup_named.AlternatingGroup>`
  - :class:`groups.permutation.KleinFour <sage.groups.perm_gps.permgroup_named.KleinFourGroup>`
  - :class:`groups.permutation.Quaternion <sage.groups.perm_gps.permgroup_named.QuaternionGroup>`
  - :class:`groups.permutation.Cyclic <sage.groups.perm_gps.permgroup_named.CyclicPermutationGroup>`
  - :class:`groups.permutation.ComplexReflection <sage.groups.perm_gps.permgroup_named.ComplexReflectionGroup>`
  - :class:`groups.permutation.Dihedral <sage.groups.perm_gps.permgroup_named.DihedralGroup>`
  - :class:`groups.permutation.DiCyclic <sage.groups.perm_gps.permgroup_named.DiCyclicGroup>`
  - :class:`groups.permutation.Mathieu <sage.groups.perm_gps.permgroup_named.MathieuGroup>`
  - :class:`groups.permutation.Suzuki <sage.groups.perm_gps.permgroup_named.SuzukiGroup>`
  - :class:`groups.permutation.PGL <sage.groups.perm_gps.permgroup_named.PGL>`
  - :class:`groups.permutation.PSL <sage.groups.perm_gps.permgroup_named.PSL>`
  - :class:`groups.permutation.PSp <sage.groups.perm_gps.permgroup_named.PSp>`
  - :class:`groups.permutation.PSU <sage.groups.perm_gps.permgroup_named.PSU>`
  - :class:`groups.permutation.PGU <sage.groups.perm_gps.permgroup_named.PGU>`
  - :class:`groups.permutation.Transitive <sage.groups.perm_gps.permgroup_named.TransitiveGroup>`
  - :class:`groups.permutation.RubiksCube <sage.groups.perm_gps.cubegroup.CubeGroup>`

- Matrix Groups (``groups.matrix.<tab>``)

  - :func:`groups.matrix.QuaternionGF3 <sage.groups.matrix_gps.finitely_generated.QuaternionMatrixGroupGF3>`
  - :func:`groups.matrix.GL <sage.groups.matrix_gps.linear.GL>`
  - :func:`groups.matrix.SL <sage.groups.matrix_gps.linear.SL>`
  - :func:`groups.matrix.Sp <sage.groups.matrix_gps.symplectic.Sp>`
  - :func:`groups.matrix.GU <sage.groups.matrix_gps.unitary.GU>`
  - :func:`groups.matrix.SU <sage.groups.matrix_gps.unitary.SU>`
  - :func:`groups.matrix.GO <sage.groups.matrix_gps.orthogonal.GO>`
  - :func:`groups.matrix.SO <sage.groups.matrix_gps.orthogonal.SO>`

- Finitely Presented Groups (``groups.presentation.<tab>``)

  - :func:`groups.presentation.Alternating <sage.groups.finitely_presented_named.AlternatingPresentation>`
  - :func:`groups.presentation.Cactus <sage.groups.finitely_presented_named.CactusPresentation>`
  - :func:`groups.presentation.Cyclic <sage.groups.finitely_presented_named.CyclicPresentation>`
  - :func:`groups.presentation.Dihedral <sage.groups.finitely_presented_named.DihedralPresentation>`
  - :func:`groups.presentation.DiCyclic <sage.groups.finitely_presented_named.DiCyclicPresentation>`
  - :func:`groups.presentation.FGAbelian <sage.groups.finitely_presented_named.FinitelyGeneratedAbelianPresentation>`
  - :func:`groups.presentation.KleinFour <sage.groups.finitely_presented_named.KleinFourPresentation>`
  - :func:`groups.presentation.Quaternion <sage.groups.finitely_presented_named.QuaternionPresentation>`
  - :func:`groups.presentation.Symmetric <sage.groups.finitely_presented_named.SymmetricPresentation>`

- Affine Groups (``groups.affine.<tab>``)

  - :func:`groups.affine.Affine <sage.groups.affine_gps.affine_group.AffineGroup>`
  - :func:`groups.affine.Euclidean <sage.groups.affine_gps.euclidean_group.EuclideanGroup>`

- Lie Groups  (``groups.lie.<tab>``)

  - :class:`groups.lie.Nilpotent <sage.groups.lie_gps.nilpotent_lie_group.NilpotentLieGroup>`

- Miscellaneous Groups (``groups.misc.<tab>``)

  - Coxeter, reflection and related groups

    - :func:`groups.misc.Braid <sage.groups.braid.BraidGroup>`
    - :func:`groups.misc.Cactus <sage.groups.cactus_group.CactusGroup>`
    - :func:`groups.misc.CoxeterGroup <sage.combinat.root_system.coxeter_group.CoxeterGroup>`
    - :func:`groups.misc.PureCactus <sage.groups.cactus_group.PureCactusGroup>`
    - :func:`groups.misc.ReflectionGroup <sage.combinat.root_system.reflection_group_real.ReflectionGroup>`
    - :class:`groups.misc.RightAngledArtin <sage.groups.raag.RightAngledArtinGroup>`
    - :func:`groups.misc.WeylGroup <sage.combinat.root_system.weyl_group.WeylGroup>`

  - other miscellaneous groups

    - :func:`groups.misc.AdditiveAbelian <sage.groups.additive_abelian.additive_abelian_group.AdditiveAbelianGroup>`
    - :class:`groups.misc.AdditiveCyclic <sage.rings.finite_rings.integer_mod_ring.IntegerModFactory>`
    - :func:`groups.misc.Free <sage.groups.free_group.FreeGroup>`
    - :func:`groups.misc.SemimonomialTransformation <sage.groups.semimonomial_transformations.semimonomial_transformation_group.SemimonomialTransformationGroup>`
"""

# Implementation notes:
#
#   With this  groups_catalog.py  module imported as
#   "groups" in all.py then  groups.<tab>  is available
#
#   "catalog" modules are made available
#   as groups.matrix, etc by imports below
#
#   Do not use this file for code
#
#   Please keep this top-level clean, use
#   groups.misc for one-off examples
#
#   Candidates for new primary divisions:
#       groups.sporadic - 26 sporadic groups
#       groups.misc - one-off stuff (implemented, but empty)
#       groups.presentation - free groups with relations
#       groups.symmetries - permutation groups of regular solids, or similar

from sage.misc.lazy_import import lazy_import

lazy_import('sage.groups.matrix_gps', 'catalog', as_='matrix')
lazy_import('sage.groups.perm_gps', 'permutation_groups_catalog', as_='permutation')
lazy_import('sage.groups.misc_gps', 'misc_groups_catalog', as_='misc')
lazy_import('sage.groups.affine_gps', 'catalog', as_='affine')
lazy_import('sage.groups', 'finitely_presented_catalog', as_='presentation')
lazy_import('sage.groups.lie_gps', 'catalog', as_='lie')

del lazy_import
