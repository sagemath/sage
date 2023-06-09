from .all__sagemath_categories import *

try:  # extra
    from sage.all__sagemath_graphs import *
except ImportError:
    pass

try:  # extra
    from sage.all__sagemath_modules import *
except ImportError:
    pass

from sage.algebras.all__sagemath_combinat import *
from sage.combinat.all__sagemath_combinat import *
from sage.rings.all__sagemath_combinat import *
from sage.monoids.all import *
from sage.games.all import *
from sage.sat.all import *
