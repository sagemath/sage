"""
Lazy import cache

This is a pure Python file with no dependencies so it can be used in setup.py.
"""
import os
import hashlib

from ..env import SAGE_LIB, DOT_SAGE


def get_cache_file():
    """
    Return the canonical filename for caching names of lazily imported
    modules.

    EXAMPLES::

        sage: from sage.misc.lazy_import_cache import get_cache_file
        sage: get_cache_file()
        '...-lazy_import_cache.pickle'
        sage: get_cache_file().startswith(DOT_SAGE)
        True
        sage: 'cache' in get_cache_file()
        True

    It should not matter whether DOT_SAGE ends with a slash::

        sage: OLD = DOT_SAGE
        sage: sage.misc.lazy_import_cache.DOT_SAGE = '/tmp'
        sage: get_cache_file().startswith('/tmp/')
        True
        sage: sage.misc.lazy_import_cache.DOT_SAGE = OLD
    """
    mangled = hashlib.sha256(os.path.realpath(SAGE_LIB).encode('utf-8')).hexdigest()
    return os.path.join(DOT_SAGE, 'cache',
                        "%s-lazy_import_cache.pickle" % mangled)
