from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import

lazy_import("sage.coding.code_constructions", ["permutation_action",
            "walsh_matrix"])

from sage.misc.superseded import deprecated_callable_import

deprecated_callable_import(19315,
            "sage.coding.code_bounds",
            globals(),
            locals(),
            ["codesize_upper_bound",
            "dimension_upper_bound",
            "volume_hamming",
            "gilbert_lower_bound",
            "plotkin_upper_bound",
            "griesmer_upper_bound",
            "elias_upper_bound",
            "hamming_upper_bound",
            "singleton_upper_bound",
            "gv_info_rate",
            "entropy",
            "gv_bound_asymp",
            "hamming_bound_asymp",
            "singleton_bound_asymp",
            "plotkin_bound_asymp",
            "elias_bound_asymp",
            "mrrw1_bound_asymp"],
            ("This method soon will not be available in that way."
            "Please call codes.bounds.%(name)s instead"))

lazy_import("sage.coding.linear_code", [
            "LinearCode",
            "LinearCodeFromVectorSpace",
            "self_orthogonal_binary_codes"])

# Functions removed from the global namespace
lazy_import('sage.coding.databases','best_linear_code_in_guava', "best_known_linear_code",
    deprecation=(21165, "best_known_linear_code has moved to sage.coding.databases.best_linear_code_in_guava"))
lazy_import('sage.coding.databases','best_linear_code_in_guava', "best_known_linear_code_www",
    deprecation=(21165, "best_known_linear_code_www has moved to sage.coding.databases.best_linear_code_in_guava"))
lazy_import('sage.coding.databases','bounds_on_minimum_distance_in_guava', "bounds_minimum_distance",
    deprecation=(21165, "bounds_minimum_distance has moved to sage.coding.databases.bounds_on_minimum_distance_in_guava"))
lazy_import('sage.coding.databases','self_orthogonal_binary_codes', "self_orthogonal_binary_codes",
    deprecation=(21165, "self_orthogonal_binary_codes has moved to sage.coding.databases.self_orthogonal_binary_codes"))
lazy_import('sage.coding.databases','self_dual_binary_codes', "self_dual_codes_binary",
    deprecation=(21165, "self_dual_codes_binary has moved to sage.coding.databases.self_dual_binary_codes"))

lazy_import("sage.coding.delsarte_bounds", [
            "Krawtchouk",
            "Kravchuk",
            "delsarte_bound_hamming_space",
            "delsarte_bound_additive_hamming_space"])

lazy_import('sage.coding', 'codes_catalog', 'codes')
lazy_import('sage.coding', 'channels_catalog', 'channels')
