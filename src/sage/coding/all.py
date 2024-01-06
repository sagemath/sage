# sage_setup: distribution = sagemath-modules

from sage.misc.lazy_import import lazy_import

lazy_import("sage.coding.code_constructions", ["permutation_action",
                                               "walsh_matrix"])

lazy_import("sage.coding.linear_code", "LinearCode")

# Functions removed from the global namespace

lazy_import('sage.coding', 'codes_catalog', 'codes')
lazy_import('sage.coding', 'channels_catalog', 'channels')

del lazy_import
