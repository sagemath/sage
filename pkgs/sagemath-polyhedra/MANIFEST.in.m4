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

include sage/numerical/mip.p*
include sage/numerical/interactive_simplex_method.p*
include sage/numerical/linear_*.p*
include sage/numerical/sdp.p*
include sage/numerical/backends/cvxopt_*.p*
include sage/numerical/backends/cvxpy_*.p*
include sage/numerical/backends/generic_*.p*
include sage/numerical/backends/interactivelp_*.p*
include sage/numerical/backends/logging_*.p*
include sage/numerical/backends/matrix_*.p*
include sage/numerical/backends/ppl_*.p*
include sage/numerical/backends/scip_*.p*


global-exclude *.py[co]
global-exclude *.so
global-exclude *.bak
