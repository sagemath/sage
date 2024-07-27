# sage_setup: distribution = sagemath-modules
from sage.misc.lazy_import import lazy_import
lazy_import("sage.numerical.optimize",
            ["find_fit", "find_local_maximum", "find_local_minimum",
             "find_root", "minimize", "minimize_constrained"])
del lazy_import
