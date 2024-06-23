from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute

# from sage.misc.all__sagemath_objects import *
from sage.misc.all__sagemath_environment import *
from sage.misc.all__sagemath_categories import *
from sage.misc.all__sagemath_repl import *

from sage.misc.misc import (BackslashOperator,
                            exists, forall, is_iterator,
                            random_sublist,
                            pad_zeros,
                            SAGE_DB,
                            newton_method_sizes, compose,
                            nest)

from sage.misc.remote_file import get_remote_file

lazy_import('sage.misc.pager', 'pager')

from sage.misc.classgraph import class_graph

from sage.misc.mathml import mathml

lazy_import("sage.misc.cython", "cython_lambda")
lazy_import("sage.misc.cython", "cython_compile", "cython")

from sage.misc.func_persist import func_persist

lazy_import('sage.repl.interpreter', 'logstr', deprecation=34259)
