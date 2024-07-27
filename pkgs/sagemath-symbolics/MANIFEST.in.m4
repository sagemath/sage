dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

prune .tox
exclude *.m4
include requirements.txt

graft sage/symbolic
# exclude what is included in sagemath-categories
exclude sage/symbolic/symbols.p*
exclude sage/symbolic/function.p*

include sage/calculus/all.p*
include sage/calculus/calculus.p*
include sage/calculus/desolvers.p*
include sage/calculus/predefined.p*
include sage/calculus/tests.p*
include sage/calculus/var.p*

graft sage/manifolds

graft sage/geometry/riemannian_manifolds
graft sage/geometry/hyperbolic_space
graft sage/dynamics/complex_dynamics

include sage/modules/vector_*symbol*.p*
include sage/matrix/matrix_symbolic_*.p*

#graft sage/libs/giac
graft sage/libs/pynac
include sage/libs/ecl.p*
include sage/libs/eclsig.h

include sage/interfaces/fricas.p*
include sage/interfaces/giac.p*
include sage/interfaces/magma*.p*
include sage/interfaces/maple*.p*
include sage/interfaces/mathematica.p*
include sage/interfaces/mathics.p*
include sage/interfaces/maxima*.p*
include sage/interfaces/qepcad.p*
include sage/interfaces/sympy*.p*
include sage/interfaces/tides.p*

include sage/interfaces/sage-maxima.lisp

graft sage/ext_data/kenzo
graft sage/ext_data/magma

global-exclude all__sagemath_*.py
global-include all__sagemath_symbolics.py

global-exclude __pycache__
global-exclude *.py[co]
global-exclude *.bak
global-exclude *.so
global-exclude *~
prune .tox
prune build
prune dist
