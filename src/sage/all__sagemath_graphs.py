try:  # extra
    from sage.all__sagemath_polyhedra import *
except ImportError:
    pass

try:  # extra
    from sage.all__sagemath_modules import *
except ImportError:
    pass

from .all__sagemath_categories import *

from sage.graphs.all     import *

from sage.topology.all   import *

from sage.combinat.all__sagemath_graphs import *

from sage.sandpiles.all import *
