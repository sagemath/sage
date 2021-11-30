"""
Utility functions for namespace packages in Sage
"""
from importlib import import_module

def install_doc(package, doc):
    """
    Install the docstring ``doc`` to the package.

    TESTS:

        sage: from sage.misc.namespace_package import install_doc
        sage: install_doc('sage', 'hello')
        sage: from inspect import getdoc
        sage: getdoc(sage)
        'hello'
    """
    pkg = import_module(package)
    pkg.__doc__ = doc         # enable sage.package?
    pkg.getdoc = lambda: doc  # enable help(sage.package)
