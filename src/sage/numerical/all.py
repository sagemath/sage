from optimize import (find_fit,
                      find_local_maximum,
                      find_local_minimum,
                      find_root,
                      linear_program,
                      minimize,
                      minimize_constrained)
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.numerical.backends.generic_backend import default_mip_solver

from sage.misc.lazy_import import lazy_import
lazy_import("sage.numerical.interactive_simplex_method",
            ["InteractiveLPProblem", "InteractiveLPProblemStandardForm"])

from sage.misc.superseded import deprecated_callable_import

# Deprecation (#17867) -- two lines at the end of
# sage.numerical.interactive_simplex_method should be removed too.
deprecated_callable_import(17867,
                           "sage.numerical.interactive_simplex_method",
                           globals(),
                           locals(),
                           ["LPProblem"],
                           ("This class meant for **educational purposes only** has been renamed to InteractiveLPProblem"))

deprecated_callable_import(17867,
                           "sage.numerical.interactive_simplex_method",
                           globals(),
                           locals(),
                           ["LPProblemStandardForm"],
                           ("This class meant for **educational purposes only** has been renamed to InteractiveLPProblemStandardForm"))
