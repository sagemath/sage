from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import as _lazy_import

_lazy_import("sage.coding.code_constructions", ["permutation_action",
            "walsh_matrix"])

_lazy_import("sage.coding.linear_code", [
            "LinearCode",
            "LinearCodeFromVectorSpace",
            "self_orthogonal_binary_codes"])

# Functions removed from the global namespace
_lazy_import('sage.coding.databases','best_linear_code_in_guava', "best_known_linear_code",
    deprecation=(21165, "best_known_linear_code has moved to sage.coding.databases.best_linear_code_in_guava"))
_lazy_import('sage.coding.databases','best_linear_code_in_guava', "best_known_linear_code_www",
    deprecation=(21165, "best_known_linear_code_www has moved to sage.coding.databases.best_linear_code_in_guava"))
_lazy_import('sage.coding.databases','bounds_on_minimum_distance_in_guava', "bounds_minimum_distance",
    deprecation=(21165, "bounds_minimum_distance has moved to sage.coding.databases.bounds_on_minimum_distance_in_guava"))
_lazy_import('sage.coding.databases','self_orthogonal_binary_codes', "self_orthogonal_binary_codes",
    deprecation=(21165, "self_orthogonal_binary_codes has moved to sage.coding.databases.self_orthogonal_binary_codes"))
_lazy_import('sage.coding.databases','self_dual_binary_codes', "self_dual_codes_binary",
    deprecation=(21165, "self_dual_codes_binary has moved to sage.coding.databases.self_dual_binary_codes"))

_lazy_import('sage.coding.delsarte_bounds','krawtchouk', "Krawtchouk",
    deprecation=(20908, "Krawtchouk will be removed from the global namespace. Please use codes.bounds.krawtchouk instead."))
_lazy_import('sage.coding.delsarte_bounds','krawtchouk', "Kravchuk",
    deprecation=(20908, "Kravchuk will be removed from the global namespace. Please use codes.bounds.krawtchouk instead."))

_lazy_import('sage.coding.delsarte_bounds',
    ["delsarte_bound_hamming_space", "delsarte_bound_additive_hamming_space"],
    deprecation=(20908, "This function will soon be removed from the global namespace. "
                "Please call it using codes.bounds.... instead"))

_lazy_import('sage.coding', 'codes_catalog', 'codes')
_lazy_import('sage.coding', 'channels_catalog', 'channels')
