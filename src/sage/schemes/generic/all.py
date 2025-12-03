# code exports
from sage.misc.lazy_import import lazy_import

from sage.schemes.generic.spec import Spec
from sage.schemes.generic.hypersurface import ProjectiveHypersurface, AffineHypersurface
lazy_import('sage.schemes.generic.igusa_top_zeta', 'ZetaFunctions')
