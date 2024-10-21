# sage_setup: distribution = sagemath-symbolics
from sage.misc.lazy_import import lazy_import

lazy_import('sage.geometry.hyperbolic_space.hyperbolic_interface', 'HyperbolicPlane')
del lazy_import
