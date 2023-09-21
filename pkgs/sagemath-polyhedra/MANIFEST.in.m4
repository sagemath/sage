dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

exclude *.m4
include requirements.txt

include sage/interfaces/latte.p*
include sage/interfaces/four_ti_2.p*

graft sage/geometry
exclude sage/geometry/abc.p*                            # in sagemath-categories
exclude sage/geometry/toric_lattice*.p*                 # in sagemath-modules
exclude sage/geometry/all.py
prune sage/geometry/hyperbolic_space
prune sage/geometry/riemannian_manifolds
exclude sage/geometry/ribbon_graph.p*                    # depends on sage.groups.perm_gps
exclude sage/geometry/integral_points_integer_dense.pyx  # depends on matrix_integer_dense (flint)

graft sage/game_theory

graft sage/numerical
exclude sage/numerical/backends/generic*backend.p*      # sagemath-categories
exclude sage/numerical/gauss_legendre.p*                # sagemath-modules
exclude sage/numerical/optimize.p*                      # sagemath-modules
exclude sage/numerical/backends/glpk*.p*                # sagemath-glpk

include sage/rings/groebner_fan.p*

graft sage/schemes/toric

global-exclude all__sagemath_*.py
global-include all__sagemath_polyhedra.py

global-exclude __pycache__
global-exclude *.py[co]
global-exclude *.bak
global-exclude *.so
global-exclude *~
prune .tox
prune build
prune dist
