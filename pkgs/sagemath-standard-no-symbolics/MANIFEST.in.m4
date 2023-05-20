dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.

dnl Include all from sagemath-categories (via m4 include)
include(`../sagelib/src/MANIFEST.in')

exclude *.m4
include requirements.txt

global-include all__sagemath_standard_no_symbolics.py

prune sage/symbolic
include sage/symbolic/all__sagemath_standard_no_symbolics.py
include sage/symbolic/symbols.p*
include sage/symbolic/function.p*

prune sage/manifolds
prune sage/lfunctions
prune sage/geometry/riemannian_manifolds
prune sage/geometry/hyperbolic_space
prune sage/dynamics/complex_dynamics
prune sage/groups/lie_gps
prune sage/modular/modform_hecketriangle
prune sage/rings/asymptotic

prune sage/plot
include sage/plot/all.p*
include sage/plot/colors.p*                     # needed by sage.graphs even when not plotting

exclude sage/calculus/all.*
exclude sage/calculus/calculus.*
exclude sage/calculus/desolvers.*
exclude sage/calculus/predefined.*
exclude sage/calculus/tests.*
exclude sage/calculus/var.*
prune sage/calculus/transforms

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
