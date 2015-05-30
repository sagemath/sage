r"""
Features that are imported by default in the interpreter namespace
"""
from rigged_configurations import RiggedConfigurations

# Deprecated classes
from rigged_configurations import HighestWeightRiggedConfigurations
from tensor_product_kr_tableaux import HighestWeightTensorProductOfKirillovReshetikhinTableaux

# Deprecations from global namespace
from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.rigged_configurations.tensor_product_kr_tableaux',
            'TensorProductOfKirillovReshetikhinTableaux',
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.TensorProductOfKirillovReshetikhinTableaux instead"))

lazy_import('sage.combinat.rigged_configurations.kr_tableaux',
            'KirillovReshetikhinTableaux',
            deprecation=(15882, "this is being removed from the global namespace. Use crystals.KirillovReshetikhin with model='KR' instead"))

