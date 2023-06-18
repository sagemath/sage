dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

prune .tox
exclude *.m4
include requirements.txt

global-include all__sagemath_polyhedra.py

include sage/interfaces/latte.p*
include sage/interfaces/four_ti_2.p*

graft sage/geometry
exclude sage/geometry/abc.py                            # in sagemath-categories
exclude sage/geometry/toric_lattice*.p*                 # in sagemath-modules
exclude sage/geometry/all.py
prune sage/geometry/hyperbolic_space
prune sage/geometry/riemannian_manifolds
exclude sage/geometry/ribbon_graph.p*                    # depends on sage.groups.perm_gps
exclude sage/geometry/integral_points_integer_dense.pyx  # depends on matrix_integer_dense

graft sage/game_theory

graft sage/numerical
exclude sage/numerical/gauss_legendre.p*                # sagemath-modules
exclude sage/numerical/optimize.p*                      # sagemath-modules
exclude sage/numerical/backends/glpk*.p*                # sagemath-glpk

global-exclude all__sagemath_categories.py
global-exclude all__sagemath_modules.py

global-exclude *.py[co]
global-exclude *.so
global-exclude *.bak
