# sage_setup: distribution = sagemath-modules
from sage.groups.all__sagemath_categories import *

from sage.groups.additive_abelian.all import *
from sage.groups.abelian_gps.all__sagemath_modules import *
from sage.groups.matrix_gps.all__sagemath_modules import *

from sage.misc.lazy_import import lazy_import

lazy_import('sage.groups.affine_gps.affine_group', 'AffineGroup')
lazy_import('sage.groups.affine_gps.euclidean_group', 'EuclideanGroup')
del lazy_import
