# sage_setup: distribution = sagemath-categories
from .constructor import FunctionField

from sage.misc.lazy_import import lazy_import

lazy_import("sage.rings.function_field.drinfeld_modules.drinfeld_module", "DrinfeldModule")
del lazy_import
