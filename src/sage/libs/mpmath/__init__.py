import os
import sys

import sage


os.environ['MPMATH_SAGE'] = '1'

# https://github.com/mpmath/mpmath/blob/master/mpmath/libmp/backend.py tries to load sage.all
if hasattr(sage, 'all'):
    sage_toplevel = None
else:
    import sage.all__sagemath_categories as sage_toplevel
    sys.modules['sage.all'] = sage.all = sage_toplevel
try:
    import mpmath.libmp.backend
finally:
    if sage_toplevel:
        del sage.all
        del sys.modules['sage.all']

assert mpmath.libmp.backend.BACKEND == 'sage'
