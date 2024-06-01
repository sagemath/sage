# sage_setup: distribution = sagemath-repl

# Set up warning filters before importing Sage stuff

import sys
import warnings

# This is a Python debug build (--with-pydebug)
__with_pydebug = hasattr(sys, 'gettotalrefcount')
if __with_pydebug:
    # a debug build does not install the default warning filters. Sadly, this breaks doctests so we
    # have to re-add them:
    warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
    warnings.filterwarnings('ignore', category=ImportWarning)
    warnings.filterwarnings('ignore', category=ResourceWarning)
else:
    deprecationWarning = ('ignore', None, DeprecationWarning, None, 0)
    if deprecationWarning in warnings.filters:
        warnings.filters.remove(deprecationWarning)

# Ignore all deprecations from IPython etc.
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        module='(IPython|ipykernel|jupyter_client|jupyter_core|nbformat|notebook|ipywidgets|storemagic|jedi)')

# scipy 1.18 introduced reprecation warnings on a number of things they are moving to
# numpy, e.g. DeprecationWarning: scipy.array is deprecated
#             and will be removed in SciPy 2.0.0, use numpy.array instead
# This affects networkx 2.2 up and including 2.4 (cf. :issue:29766)
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        module='(scipy|networkx)')

# However, be sure to keep OUR deprecation warnings
warnings.filterwarnings('default', category=DeprecationWarning,
                        message=r'[\s\S]*See https?://trac\.sagemath\.org/[0-9]* for details.')

# Ignore Python 3.9 deprecation warnings
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        module='ast')

# Ignore packaging 20.5 deprecation warnings
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        module='(.*[.]_vendor[.])?packaging')

# Ignore a few warnings triggered by pythran 0.12.1
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message='\n\n  `numpy.distutils` is deprecated since NumPy 1.23.0',
                        module='pythran.dist')
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message='pkg_resources is deprecated as an API|'
                        'Deprecated call to `pkg_resources.declare_namespace(.*)`',
                        module='pkg_resources|setuptools.sandbox')
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message='msvccompiler is deprecated and slated to be removed',
                        module='distutils.msvccompiler')

warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message='The distutils(.sysconfig module| package) is deprecated',
                        module='Cython|distutils|numpy|sage.env|sage.features')

# triggered by cython 0.29.32
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message="'cgi' is deprecated and slated for removal in Python 3.13",
                        module='Cython')

# triggered by pyparsing 2.4.7
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message="module 'sre_constants' is deprecated",
                        module='pyparsing')

# importlib.resources.path and ...read_binary are deprecated in python 3.11,
# but the replacement importlib.resources.files needs python 3.9
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message=r'(path|read_binary) is deprecated\. Use files\(\) instead\.',
                        module='sage.repl.rich_output.output_(graphics|graphics3d|video)')

# triggered by sphinx
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message="'imghdr' is deprecated and slated for removal in Python 3.13",
                        module='sphinx.util.images')

# triggered by docutils 0.19 on Python 3.11
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message=r"Use setlocale\(\), getencoding\(\) and getlocale\(\) instead",
                        module='docutils.io')

# triggered by dateutil 2.8.2 and sphinx 7.0.1 on Python 3.12
# see: https://github.com/dateutil/dateutil/pull/1285
# see: https://github.com/sphinx-doc/sphinx/pull/11468
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message=r"datetime.datetime.utcfromtimestamp\(\) is deprecated",
                        module='dateutil.tz.tz|sphinx.(builders.gettext|util.i18n)')

# triggered on Python 3.12
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message=r"This process.* is multi-threaded, "
                                r"use of .*\(\) may lead to deadlocks in the child.")

# pickling of itertools is deprecated in Python 3.12
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message=r"Pickle, copy, and deepcopy support will be "
                                r"removed from itertools in Python 3.14.")

# triggered in Python 3.9 on Redhat-based distributions
# https://github.com/sagemath/sage/issues/37863
# https://github.com/networkx/networkx/issues/7101
warnings.filterwarnings('ignore', category=RuntimeWarning,
                        message="networkx backend defined more than once: nx-loopback")

from sage.all__sagemath_objects import *
from sage.all__sagemath_environment import *

from sage.doctest.all import *
from sage.repl.all import *
from sage.misc.all__sagemath_repl import *

# For doctesting. These are overwritten later

Integer = int
RealNumber = float
