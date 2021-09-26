# sage.cpython is an ordinary package, not a namespace package.

# This package is imported very early, which is why workarounds/monkey-patching
# are done in this file.

# Make sure that the correct zlib library is loaded. This is needed
# to prevent the system zlib to be loaded instead of the Sage one.
# See https://trac.sagemath.org/ticket/23122
import zlib as _zlib
del _zlib

# Monkey-patch ExtensionFileLoader to allow IPython to find the sources
# of Cython files. See https://trac.sagemath.org/ticket/24681
from importlib.machinery import ExtensionFileLoader as _ExtensionFileLoader
del _ExtensionFileLoader.get_source
del _ExtensionFileLoader

# Work around a Cygwin-specific bug caused by sqlite3; see
# https://trac.sagemath.org/ticket/30157 and the docstring for
# fix_for_ticket_30157
# Here we monkey-patch the sqlite3 module to ensure the fix is
# applied the very first time a connection is made to a sqlite3
# database
import sys as _sys
if _sys.platform == 'cygwin':
    def _patch_sqlite3():
        try:
            from sage.misc.sage_ostools import fix_for_ticket_30157
        except ImportError:
            # The module might not have been re-built yet; don't worry about it
            # then
            return

        import sqlite3
        import functools
        orig_sqlite3_connect = sqlite3.connect

        @functools.wraps(orig_sqlite3_connect)
        def connect(*args, **kwargs):
            if fix_for_ticket_30157():
                raise RuntimeError(
                    'patch for Trac ticket #30157 failed; please report this '
                    'bug to https://trac.sagemath.org')

            # Undo the monkey-patch
            try:
                return orig_sqlite3_connect(*args, **kwargs)
            finally:
                sqlite3.connect = orig_sqlite3_connect

        sqlite3.connect = connect

    _patch_sqlite3()
    del _patch_sqlite3
    del _sys
