# sage_setup: distribution = sagemath-symbolics
from sage.misc.lazy_import import lazy_import

lazy_import('sage.libs.giac.giac', 'libgiac')
del lazy_import
