# sage_setup: distribution = sagemath-environment
r"""
Enumeration of all defined features
"""

# *****************************************************************************
#       Copyright (C) 2021-2023 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import itertools


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


def module_feature(module_name):
    r"""
    Find a top-level :class:`Feature` that provides the Python module of the given ``module_name``.

    Only features known to :func:`all_features` are considered.

    INPUT:

    - ``module_name`` -- string

    OUTPUT: a :class:`Feature` or ``None``

    EXAMPLES::

        sage: from sage.features.all import module_feature
        sage: module_feature('sage.combinat.tableau')                                   # needs sage.combinat
        Feature('sage.combinat')
        sage: module_feature('sage.combinat.posets.poset')                              # needs sage.graphs
        Feature('sage.graphs')
        sage: module_feature('sage.schemes.toric.variety')                              # needs sage.geometry.polyhedron
        Feature('sage.geometry.polyhedron')
        sage: module_feature('scipy')                                                   # needs scipy
        Feature('scipy')
        sage: print(module_feature('sage.structure.element'))
        None
        sage: print(module_feature('sage.does_not_exist'))
        None
    """
    longest_prefix = ''
    longest_prefix_feature = None
    for feature in all_features():
        for joined in itertools.chain([feature], feature.joined_features()):
            if joined.name == module_name:
                return feature
            if (joined.name + '.').startswith(longest_prefix):
                if (module_name + '.').startswith(joined.name + '.'):
                    longest_prefix = feature.name + '.'
                    longest_prefix_feature = feature
    return longest_prefix_feature


def name_feature(name, toplevel=None):
    r"""
    Find a top-level :class:`Feature` that provides the top-level ``name``.

    Only features known to :func:`all_features` are considered.

    INPUT:

    - ``name`` -- string

    - ``toplevel`` -- a module or other namespace

    OUTPUT: a :class:`Feature` or ``None``

    EXAMPLES::

        sage: from sage.features.all import name_feature
        sage: name_feature('QuadraticField')                                            # needs sage.rings.number_field
        Feature('sage.rings.number_field')
        sage: name_feature('line')                                                      # needs sage.plot
        Feature('sage.plot')
        sage: print(name_feature('ZZ'))
        None
        sage: print(name_feature('does_not_exist'))
        None
    """
    if toplevel is None:
        try:
            import sage.all as toplevel
        except ImportError:
            return None
    try:
        obj = getattr(toplevel, name)
    except AttributeError:
        return None

    from sage.misc.sageinspect import find_object_modules

    for module, names in find_object_modules(obj).items():
        if name in names and (feature := module_feature(module)):
            return feature

    return None
