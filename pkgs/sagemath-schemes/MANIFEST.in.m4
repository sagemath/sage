dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

prune .tox
exclude *.m4
include requirements.txt

global-include all__sagemath_schemes.py

graft sage/schemes
graft sage/modular
graft sage/dynamics/arithmetic_dynamics
# Included in sagemath-categories
prune sage/schemes/affine
prune sage/schemes/projective
prune sage/schemes/generic
# in sagemath-polyhedra
prune sage/schemes/toric
# Has compile-time dependencies on sagemath-ntl
exclude sage/schemes/hyperelliptic_curves/hypellfrob.pyx
# Has compile-time dependencies on flint
exclude sage/modular/modform/eis_series_cython.p*
exclude sage/modular/modsym/apply.p*
exclude sage/modular/modsym/heilbronn.p*
exclude sage/modular/pollack_stevens/dist.p*
exclude sage/schemes/elliptic_curves/descent_two_isogeny.p*
#exclude sage/modular/arithgroup/arithgroup_element.pyx          # via Matrix_integer_dense
#exclude sage/modular/arithgroup/congroup.pyx                    # via Matrix_integer_dense


global-exclude *.py[co]
global-exclude *.so
global-exclude *.bak
