
from sage.features import FeatureNotPresentError
from sage.rings.number_field.number_field import GaussianField

try:
    I = GaussianField().gen()
except FeatureNotPresentError:
    I = None # Needs NTL
