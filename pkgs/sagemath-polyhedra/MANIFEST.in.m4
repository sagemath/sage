dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.

dnl Include all from sagemath-categories (via m4 include)
include(`../sagemath_categories/src/MANIFEST.in')

# Extra in sagemath-polyhedra:

include sage/rings/integer*.*
include sage/rings/rational*.*
#include sage/rings/infinity.*

graft sage/data_structures
exclude sage/data_structures/bounded_integer_sequences.*   # depends on flint

graft sage/geometry
prune sage/geometry/hyperbolic_space
prune sage/geometry/riemannian_manifolds
