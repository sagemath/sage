r"""
Enumeration of all defined features
"""

def all_features():
    r"""
    Return an iterable of all features.

    EXAMPLES::

        sage: from sage.features.all import all_features
        sage: sorted(all_features(), key=lambda f: f.name)  # random
        [...Feature('sage.combinat')...]
    """
    import pkgutil
    import importlib
    import sage.features
    # Following https://packaging.python.org/guides/creating-and-discovering-plugins/#using-namespace-packages
    for finder, name, ispkg in pkgutil.iter_modules(sage.features.__path__, sage.features.__name__ + "."):
        module = importlib.import_module(name)
        try:
            af = module.all_features
        except AttributeError:
            pass
        else:
            if af != all_features:
                yield from af()
