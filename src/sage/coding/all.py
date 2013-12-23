from sage.misc.lazy_import import lazy_import

from ag_code import ag_code

from code_constructions import (permutation_action,
                   walsh_matrix,cyclotomic_cosets)

from sage.misc.superseded import deprecated_callable_import
deprecated_callable_import(15445,
                           'sage.coding.code_constructions',
                           globals(),
                           locals(),
                           ["BinaryGolayCode",
                            "BCHCode",
                            "CyclicCode",
                            "CyclicCodeFromGeneratingPolynomial",
                            "CyclicCodeFromCheckPolynomial",
                            "DuadicCodeEvenPair",
                            "DuadicCodeOddPair",
                            "ExtendedBinaryGolayCode",
                            "ExtendedQuadraticResidueCode",
                            "ExtendedTernaryGolayCode",
                            "HammingCode",
                            "LinearCodeFromCheckMatrix",
                            "QuadraticResidueCode",
                            "QuadraticResidueCodeEvenPair",
                            "QuadraticResidueCodeOddPair",
                            "RandomLinearCode",
                            "ReedSolomonCode",
                            "TernaryGolayCode",
                            "ToricCode",
                            "TrivialCode",
                            "WalshCode"],
                           ("This method soon will not be available in that "
                            "way anymore. To use it, you can now call it by "
                            "typing codes.%(name)s"))

deprecated_callable_import(15445,
                           'sage.coding.guava',
                           globals(),
                           locals(),
                           ["BinaryReedMullerCode",
                            "QuasiQuadraticResidueCode",
                            "RandomLinearCodeGuava"],
                            ("This method soon will not be available in that "
                            "way anymore. To use it, you can now call it by "
                            "typing codes.%(name)s"))

from code_bounds import (codesize_upper_bound,
                         dimension_upper_bound,
                         volume_hamming,
                         gilbert_lower_bound,
                         plotkin_upper_bound,
                         griesmer_upper_bound,
                         elias_upper_bound,
                         hamming_upper_bound,
                         singleton_upper_bound,
                         gv_info_rate,
                         entropy,
                         gv_bound_asymp,
                         hamming_bound_asymp,
                         singleton_bound_asymp,
                         plotkin_bound_asymp,
                         elias_bound_asymp,
                         mrrw1_bound_asymp)

from linear_code import (LinearCode, LinearCodeFromVectorSpace,
                         hamming_weight,
                         best_known_linear_code,
                         best_known_linear_code_www,
                         bounds_minimum_distance,
                         self_orthogonal_binary_codes)

from sd_codes import self_dual_codes_binary

lazy_import("sage.coding.delsarte_bounds",
    ["Krawtchouk", "delsarte_bound_hamming_space", "delsarte_bound_additive_hamming_space"])

lazy_import('sage.coding', 'codes_catalog', 'codes')
