# sage_setup: distribution = sagemath-repl

from sage.misc.all__sagemath_objects import *

from sage.misc.sage_eval import sage_eval, sageobj

from sage.misc.sage_input import sage_input

from sage.misc.banner import version

lazy_import('sage.misc.sagedoc', ['browse_sage_doc',
                                  'search_src', 'search_def', 'search_doc',
                                  'tutorial', 'reference', 'manual', 'developer',
                                  'constructions', 'help'])

lazy_import('pydoc', 'help', 'python_help')
