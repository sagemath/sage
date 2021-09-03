dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.

dnl Include all from sagemath-categories (via m4 include)
include(`../sagemath_categories/src/MANIFEST.in')

# Extra in sagemath-polyhedra:

include sage/rings/integer*.*
include sage/rings/rational*.*
#include sage/rings/infinity.*
include sage/arith/rational_reconstruction.*
include sage/misc/allocator.*
include sage/ext/mod_int.*

include sage/rings/finite_rings/__init__.py
include sage/rings/finite_rings/stdint.*
include sage/rings/finite_rings/integer_mod_limits.h
include sage/rings/finite_rings/integer_mod.pxd   # .pyx depends on pari

graft sage/modules
exclude sage/modules/vector_*double*.*  # depends on numpy
exclude sage/modules/vector_mod2*.*     # depends on m4ri
exclude sage/modules/vector_*symbol*.*  # --> sagemath-symbolics

graft sage/matrix
exclude sage/matrix/misc.*  # until refactored
exclude sage/matrix/args.pyx  # until refactored
exclude sage/matrix/matrix_gap.*
exclude sage/matrix/matrix_*ball*.*     # depends on arb
exclude sage/matrix/matrix_*double*.*   # depends on numpy
exclude sage/matrix/matrix_*cyclo*.*    # depends on ntl
exclude sage/matrix/matrix_*gap*.*      # depends on gap
exclude sage/matrix/matrix_gf2*.*       # depends on m4ri, m4rie
exclude sage/matrix/matrix_gfpn*.*      # depends on meataxe
exclude sage/matrix/matrix_integer_*.*  # depends on flint, pari, iml, linbox
exclude sage/matrix/matrix_mod2*.*      # depends on m4ri
exclude sage/matrix/matrix_modn*.*      # depends on linbox or flint
exclude sage/matrix/matrix_mpolynom*.*  # depends on singular
exclude sage/matrix/matrix_rational_*.* # depends on flint, pari
exclude sage/matrix/matrix_symbolic_*.* # --> sagemath-symbolics

graft sage/data_structures
exclude sage/data_structures/bounded_integer_sequences.*   # depends on flint

graft sage/geometry
prune sage/geometry/hyperbolic_space
prune sage/geometry/riemannian_manifolds
exclude sage/geometry/integral_points.pyx  # depends on matrix_integer_dense
