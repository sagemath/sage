# sage_setup: distribution = sagemath-modules
# All of sage.misc.all except for development tools, session management,
# and deprecated functionality

from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute

from sage.misc.all__sagemath_categories import *

from sage.misc.misc import (BackslashOperator,               # Depends on sage.env -- can lower to sagemath-objects after splitting this module
                            exists, forall, is_iterator,
                            random_sublist,
                            pad_zeros,
                            SAGE_DB,
                            newton_method_sizes, compose,
                            nest)

from sage.misc.temporary_file import tmp_dir, tmp_filename  # Depends on sage.env

from sage.misc.mathml import mathml

from sage.misc.func_persist import func_persist
