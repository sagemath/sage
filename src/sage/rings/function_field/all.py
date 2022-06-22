from .constructor import FunctionField

from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.function_field.drinfeld_module', 'DrinfeldModule')
lazy_import('sage.rings.function_field.drinfeld_module', 'DrinfeldModuleAction')
