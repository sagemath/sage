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
# exclude what is in sagemath-modules
prune sage/groups/additive_abelian
prune sage/groups/affine_gps
prune sage/groups/matrix_gps
# exclude what is in sagemath-gap
prune sage/groups/perm_gps
exclude sage/groups/*gap*.p*

global-include all__sagemath_groups.py

global-exclude *.py[co]
global-exclude *.so
global-exclude *.bak
