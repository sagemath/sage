dnl MANIFEST.in is generated from this file by SAGE_ROOT/bootstrap via m4.
prune sage

include VERSION.txt

prune .tox
exclude *.m4
include requirements.txt

global-include all__sagemath_symbolics.py

graft sage/symbolic
# exclude what is included in sagemath-categories
exclude sage/symbolic/symbols.p*
exclude sage/symbolic/function.p*

graft sage/calculus
# exclude what is included in sagemath-categories
exclude sage/calculus/functional.p*

graft sage/manifolds

graft sage/geometry/riemannian_manifolds
graft sage/geometry/hyperbolic_space
graft sage/dynamics/complex_dynamics

#include sage/modules/vector_*symbol*.p*
#include sage/matrix/matrix_symbolic_*.p*

graft sage/libs/giac
graft sage/libs/gsl
graft sage/libs/pynac
include sage/libs/ecl.p*
include sage/libs/eclsig.h

graft sage/interfaces
# include sage/interfaces/fricas.p*
# include sage/interfaces/giac.p*
# include sage/interfaces/magma*.p*
# include sage/interfaces/maple*.p*
# include sage/interfaces/mathematica.p*
# include sage/interfaces/mathics.p*
# include sage/interfaces/maxima*.p*
# include sage/interfaces/qepcad.p*
# include sage/interfaces/sympy*.p*
# include sage/interfaces/tides.p*
exclude sage/interfaces/four_ti_2.p*
exclude sage/interfaces/kenzo.py
exclude sage/interfaces/gap.py
exclude sage/interfaces/tachyon.py
exclude sage/interfaces/singular.py

# Temporary:
graft sage/ext/interpreters
include sage/ext/fast*.p*
graft sage/libs/mpfr
graft sage/libs/mpc
graft sage/libs/mpmath
include sage/rings/real_mpfr.p*
include sage/rings/real_field.p*
include sage/rings/cc.p*
include sage/rings/complex_double.p*
include sage/rings/complex_field.p*
include sage/rings/complex_mpfr.p*
include sage/rings/complex_conversion.p*
include sage/rings/polynomial/*_mpfr_*.p*


# Including singular is too tricky; singular.pyx pulls in pari, ntl, flint, givaro
#graft sage/libs/singular
#include sage/rings/polynomial/*_libsingular.p*


global-exclude *.py[co]
global-exclude *.so
global-exclude *.bak
