# sage_setup: distribution = sagemath-gap

from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.universal_cyclotomic_field', ['UniversalCyclotomicField', 'E'])

del lazy_import
