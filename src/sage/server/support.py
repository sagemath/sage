"""
Support for Notebook Introspection and Setup

AUTHORS:

- William Stein (much of this code is from IPython).

- Nick Alexander

TESTS:

Test deprecation::

    sage: from sage.server.support import syseval
    sage: syseval(gap, "2+3")  # long time
    doctest:...: DeprecationWarning: 
    Importing syseval from here is deprecated. If you need to use it, please import it directly from sagenb.misc.support
    See http://trac.sagemath.org/2891 for details.
    '5'
"""

EMBEDDED_MODE = False

from sage.misc.lazy_import import lazy_import
lazy_import('sagenb.misc.support',
        ('sage_globals', 'globals_at_init',
        'global_names_at_init', 'init', 'setup_systems', 'help',
        'get_rightmost_identifier', 'completions', 'docstring',
        'source_code', 'syseval'),
    deprecation=2891)

lazy_import('sage.misc.cython', ('cython_import', 'cython_import_all'),
    deprecation=9552)
