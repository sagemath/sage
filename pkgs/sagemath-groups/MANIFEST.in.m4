dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

prune .tox
exclude *.m4
include requirements.txt

graft sage/groups
# exclude what is in sagemath-objects
exclude sage/groups/group.*
exclude sage/groups/old.*
# exclude what is in sagemath-categories
exclude sage/groups/generic.p*
exclude sage/groups/groups_catalog.p*
# exclude what is in sagemath-modules
prune sage/groups/abelian_gps
include sage/groups/abelian_gps/all.p*
exclude sage/groups/galois_group.p*
prune sage/groups/additive_abelian
prune sage/groups/affine_gps
prune sage/groups/matrix_gps
# exclude what is in sagemath-gap
prune sage/groups/perm_gps
exclude sage/groups/*gap*.p*
exclude sage/groups/galois_group_perm.p*

global-exclude all__sagemath_*.py
global-include all__sagemath_groups.py

global-exclude *.py[co]
global-exclude *.so
global-exclude *.bak
