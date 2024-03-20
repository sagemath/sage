# sage_setup: distribution = sagemath-objects
# sage.cpython is an ordinary package, not a namespace package.

# This package is imported very early, which is why workarounds/monkey-patching
# are done in this file.

# Make sure that the correct zlib library is loaded. This is needed
# to prevent the system zlib to be loaded instead of the Sage one.
# See https://github.com/sagemath/sage/issues/23122
import zlib as _zlib
del _zlib

# Monkey-patch ExtensionFileLoader to allow IPython to find the sources
# of Cython files. See https://github.com/sagemath/sage/issues/24681
from importlib.machinery import ExtensionFileLoader as _ExtensionFileLoader
if hasattr(_ExtensionFileLoader, 'get_source'):
    del _ExtensionFileLoader.get_source
del _ExtensionFileLoader
