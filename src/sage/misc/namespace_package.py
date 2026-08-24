"""
Utility functions for namespace packages in Sage
"""
from collections.abc import Mapping
from importlib import import_module


def install_doc(package: str, doc: str) -> None:
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


def install_dict(package: str, dic: Mapping[str, object]) -> None:
    """
    Install ``dic`` to the ``__dict__`` of the package.

    TESTS:

        sage: from sage.misc.namespace_package import install_dict
        sage: install_dict('sage', {'greeting': 'hello'})
        sage: sage.greeting
        'hello'
    """
    pkg = import_module(package)
    pkg.__dict__.update(dic)
