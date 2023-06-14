dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.

dnl Include all from sagemath-categories (via m4 include)
include(`../sagelib/src/MANIFEST.in')

exclude *.m4
include requirements.txt

global-include all__sagemath_standard_no_symbolics.py

prune sage/symbolic
prune sage/manifolds
prune sage/lfunctions
prune sage/geometry/riemannian_manifolds
prune sage/geometry/hyperbolic_space
prune sage/dynamics/complex_dynamics
prune sage/groups/lie_gps
prune sage/modular/modform_hecketriangle
prune sage/rings/asymptotic

prune sage/plot

exclude sage/modules/vector_*symbol*.*
exclude sage/matrix/matrix_symbolic_*.*
exclude sage/groups/misc_gps/argument_groups.*
exclude sage/groups/misc_gps/imaginary_groups.*

prune sage/libs/giac
prune sage/libs/pynac
exclude sage/libs/ecl.p*

exclude sage/interfaces/fricas.p*
exclude sage/interfaces/giac.p*
exclude sage/interfaces/magma*.p*
exclude sage/interfaces/maple*.p*
exclude sage/interfaces/mathematica.p*
exclude sage/interfaces/mathics.p*
exclude sage/interfaces/maxima*.p*
exclude sage/interfaces/qepcad.p*
exclude sage/interfaces/sympy*.p*
exclude sage/interfaces/tides.p*


# Exclude what is included in other distros
prune sage/algebras/finite_dimensional_algebras
prune sage/arith
prune sage/calculus
prune sage/categories
prune sage/coding/guruswami_sudan
prune sage/combinat/crystals
prune sage/combinat/designs
prune sage/combinat/integer_lists
prune sage/combinat/rigged_configurations
prune sage/combinat/sf
prune sage/combinat/species
prune sage/combinat/words
prune sage/crypto
prune sage/doctest
prune sage/dynamics/arithmetic_dynamics
prune sage/features
prune sage/functions
prune sage/game_theory
prune sage/games
prune sage/geometry/hyperplane_arrangement
prune sage/geometry/polyhedron
prune sage/geometry/triangulation
prune sage/groups/additive_abelian
prune sage/libs/gsl
prune sage/libs/lrcalc
prune sage/libs/mpfr
prune sage/libs/mpmath
prune sage/libs/symmetrica
prune sage/matroids
prune sage/monoids
prune sage/parallel
prune sage/repl
prune sage/rings/semirings
prune sage/sandpiles
prune sage/sat
prune sage/schemes/affine
prune sage/schemes/product_projective
prune sage/schemes/projective
prune sage/sets
prune sage/stats/distributions
prune sage/structure
prune sage/symbolic
prune sage/tensor
prune sage/topology

include sage/calculus/all__sagemath_standard_no_symbolics.py
include sage/calculus/integration.p*

include sage/*/all.py
include sage/*/*/all.py
include sage/*/*/*/all.py
