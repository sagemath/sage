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

lazy_import("sage.coding.linear_code", ["LinearCode",
            "LinearCodeFromVectorSpace",
            "best_known_linear_code",
            "best_known_linear_code_www",
            "bounds_minimum_distance",
            "self_orthogonal_binary_codes"])

lazy_import("sage.coding.delsarte_bounds", ["Krawtchouk", "Kravchuk",
            "delsarte_bound_hamming_space",
            "delsarte_bound_additive_hamming_space"])

from .sd_codes import self_dual_codes_binary
lazy_import('sage.coding', 'codes_catalog', 'codes')
lazy_import('sage.coding', 'channels_catalog', 'channels')
