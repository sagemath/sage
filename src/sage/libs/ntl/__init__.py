from sage.features import FeatureNotPresentError
from sage.libs.ntl.error import setup_NTL_error_callback

try:
    setup_NTL_error_callback()
except FeatureNotPresentError:
    pass
