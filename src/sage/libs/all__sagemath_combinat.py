# sage_setup: distribution = sagemath-combinat
from sage.misc.lazy_import import lazy_import

lazy_import('sage.libs.symmetrica', 'all', as_='symmetrica')

del lazy_import
