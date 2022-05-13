from .constructor import FunctionField

from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.function_field.finite_drinfeld_module', 'FiniteDrinfeldModule')
lazy_import('sage.rings.function_field.finite_drinfeld_module',
        'FiniteDrinfeldModule_rank_two')
lazy_import('sage.rings.function_field.finite_drinfeld_module', 'FiniteDrinfeldModuleAction')
