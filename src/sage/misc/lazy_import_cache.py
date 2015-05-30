"""
Lazy import cache

This is a pure Python file with no dependencies so it can be used in setup.py.
"""

import os

from ..env import SAGE_SRC, DOT_SAGE

def get_cache_file():
    """
    Returns a per-branch file for caching names of lazily imported modules.

    EXAMPLES::

        sage: from sage.misc.lazy_import_cache import get_cache_file
        sage: get_cache_file()
        '...-lazy_import_cache.pickle'
        sage: get_cache_file().startswith(DOT_SAGE)
        True
        sage: 'src' in get_cache_file()
        True

    It shouldn't matter whether DOT_SAGE ends with a slash::

        sage: OLD = DOT_SAGE
        sage: sage.misc.lazy_import_cache.DOT_SAGE = '/tmp'
        sage: get_cache_file().startswith('/tmp/')
        True
        sage: sage.misc.lazy_import_cache.DOT_SAGE = OLD
    """
    mangled = os.path.realpath(SAGE_SRC).replace(os.sep, '_')
    return os.path.join(DOT_SAGE, 'cache',
                        "%s-lazy_import_cache.pickle" % mangled)
