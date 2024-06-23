# sage_setup: distribution = sagemath-categories
from sage.all__sagemath_objects import *

try:
    # For doctesting
    from sage.all__sagemath_repl import *
except ImportError:
    pass

from sage.categories.all import *
from sage.rings.all__sagemath_categories import *
from sage.sets.all import *
from sage.combinat.all__sagemath_categories import *
from sage.arith.all import *
from sage.data_structures.all import *
from sage.ext.all__sagemath_categories import *
from sage.groups.all__sagemath_categories import *
from sage.interfaces.all import *
from sage.misc.all__sagemath_categories import *
from sage.typeset.all import *
from sage.schemes.all__sagemath_categories import *

from sage.calculus.all__sagemath_categories import *
from sage.functions.all import *

from sage.parallel.all import *
