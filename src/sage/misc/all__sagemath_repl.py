from .sage_eval import sage_eval, sageobj

from .sage_input import sage_input

from .banner import version, banner

lazy_import('sage.misc.sagedoc', ['browse_sage_doc',
        'search_src', 'search_def', 'search_doc',
        'tutorial', 'reference', 'manual', 'developer',
        'constructions', 'help'])

lazy_import('pydoc', 'help', 'python_help')
