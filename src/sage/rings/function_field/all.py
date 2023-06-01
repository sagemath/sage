from .all__sagemath_modules import *

from sage.misc.lazy_import import lazy_import

lazy_import("sage.rings.function_field.drinfeld_modules.drinfeld_module", "DrinfeldModule")
