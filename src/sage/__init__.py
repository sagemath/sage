__all__ = ['all']

# Set sage.__version__ to the current version number. This is analogous
# to many other Python packages.
from sage.version import version as __version__

import sys
# Make sure that the correct zlib library is loaded. This is needed
# to prevent the system zlib to be loaded instead of the Sage one.
# See https://trac.sagemath.org/ticket/23122
import zlib

# IPython calls this when starting up
def load_ipython_extension(*args):
    import sage.repl.ipython_extension
    sage.repl.ipython_extension.load_ipython_extension(*args)


# Monkey-patch inspect.isfunction() to support Cython functions.
def isfunction(obj):
    """
    Check whether something is a function.

    We assume that anything which has a genuine ``__code__``
    attribute (not using ``__getattr__`` overrides) is a function.
    This is meant to support Cython functions.

    EXAMPLES::

        sage: from inspect import isfunction
        sage: def f(): pass
        sage: isfunction(f)
        True
        sage: isfunction(lambda x:x)
        True
        sage: from sage.categories.coercion_methods import _mul_parent
        sage: isfunction(_mul_parent)
        True
        sage: isfunction(Integer.digits)     # unbound method
        False
        sage: isfunction(Integer(1).digits)  # bound method
        False

    Verify that ipywidgets can correctly determine signatures of Cython
    functions::

        sage: from ipywidgets.widgets.interaction import signature
        sage: from sage.dynamics.complex_dynamics.mandel_julia_helper import fast_mandelbrot_plot
        sage: signature(fast_mandelbrot_plot)  # random
        <IPython.utils._signatures.Signature object at 0x7f3ec8274e10>
    """
    # We use type(obj) instead of just obj to avoid __getattr__().
    # Some types, like methods, will return the __code__ of the
    # underlying function in __getattr__() but we don't want to
    # detect those as functions.
    return hasattr(type(obj), "__code__")

import inspect
inspect.isfunction = isfunction


# Monkey-patch ExtensionFileLoader to allow IPython to find the sources
# of Cython files. See https://trac.sagemath.org/ticket/24681
try:
    from importlib.machinery import ExtensionFileLoader
except ImportError:
    pass  # Python 2
else:
    del ExtensionFileLoader.get_source


# Work around a Cygwin-specific bug caused by sqlite3; see
# https://trac.sagemath.org/ticket/30157 and the docstring for
# fix_for_ticket_30157
# Here we monkey-patch the sqlite3 module to ensure the fix is
# applied the very first time a connection is made to a sqlite3
# database
if sys.platform == 'cygwin':
    def patch_sqlite3():
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

    patch_sqlite3()
    del patch_sqlite3
