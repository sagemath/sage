from sage.misc.lazy_import import lazy_import

from .generic_ds import DynamicalSystem
from .affine_ds import DynamicalSystem_affine
from .projective_ds import DynamicalSystem_projective
from .product_projective_ds import DynamicalSystem_product_projective
from .berkovich_ds import DynamicalSystem_Berkovich
lazy_import('sage.dynamics.arithmetic_dynamics.wehlerK3', 'WehlerK3Surface')
lazy_import('sage.dynamics.arithmetic_dynamics.wehlerK3', 'random_WehlerK3Surface')
