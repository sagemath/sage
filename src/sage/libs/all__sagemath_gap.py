# sage_setup: distribution = sagemath-gap
from sage.misc.lazy_import import lazy_import
lazy_import('sage.libs.gap.libgap', 'libgap')
del lazy_import
