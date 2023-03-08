r"""
Join features
"""

# ****************************************************************************
#       Copyright (C) 2021-2022 Matthias Koeppe
#                     2021-2022 Kwankyu Lee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import Feature, FeatureTestResult


class JoinFeature(Feature):
    r"""
    Join of several :class:`~sage.features.Feature` instances.

    This creates a new feature as the union of the given features. Typically
    these are executables of an SPKG. For an example, see
    :class:`~sage.features.rubiks.Rubiks`.

    Furthermore, this can be the union of a single feature. This is used to map
    the given feature to a more convenient name to be used in ``optional`` tags
    of doctests. Thus you can equip a feature such as a
    :class:`~sage.features.PythonModule` with a tag name that differs from the
    systematic tag name. As an example for this use case, see
    :class:`~sage.features.meataxe.Meataxe`.

    EXAMPLES::

        sage: from sage.features import Executable
        sage: from sage.features.join_feature import JoinFeature
        sage: F = JoinFeature("shell-boolean",
        ....:                 (Executable('shell-true', 'true'),
        ....:                  Executable('shell-false', 'false')))
        sage: F.is_present()
        FeatureTestResult('shell-boolean', True)
        sage: F = JoinFeature("asdfghjkl",
        ....:                 (Executable('shell-true', 'true'),
        ....:                  Executable('xxyyyy', 'xxyyyy-does-not-exist')))
        sage: F.is_present()
        FeatureTestResult('xxyyyy', False)
    """
    def __init__(self, name, features, spkg=None, url=None, description=None):
        """
        TESTS:

        The empty join feature is present::

            sage: from sage.features.join_feature import JoinFeature
            sage: JoinFeature("empty", ()).is_present()
            FeatureTestResult('empty', True)
        """
        if spkg is None:
            spkgs = set(f.spkg for f in features if f.spkg)
            if len(spkgs) > 1:
                raise ValueError('given features have more than one spkg; provide spkg argument')
            elif len(spkgs) == 1:
                spkg = next(iter(spkgs))
        if url is None:
            urls = set(f.url for f in features if f.url)
            if len(urls) > 1:
                raise ValueError('given features have more than one url; provide url argument')
            elif len(urls) == 1:
                url = next(iter(urls))
        super().__init__(name, spkg=spkg, url=url, description=description)
        self._features = features

    def _is_present(self):
        r"""
        Test for the presence of the join feature.

        EXAMPLES::

            sage: from sage.features.latte import Latte
            sage: Latte()._is_present()  # optional - latte_int
            FeatureTestResult('latte_int', True)
        """
        for f in self._features:
            test = f._is_present()
            if not test:
                return test
        return FeatureTestResult(self, True)

    def is_functional(self):
        r"""
        Test whether the join feature is functional.

        This method is deprecated. Use :meth:`Feature.is_present` instead.

        EXAMPLES::

            sage: from sage.features.latte import Latte
            sage: Latte().is_functional()  # optional - latte_int
            doctest:warning...
            DeprecationWarning: method JoinFeature.is_functional; use is_present instead
            See https://github.com/sagemath/sage/issues/33114 for details.
            FeatureTestResult('latte_int', True)
        """
        try:
            from sage.misc.superseded import deprecation
        except ImportError:
            # The import can fail because sage.misc.superseded is provided by
            # the distribution sagemath-objects, which is not an
            # install-requires of the distribution sagemath-environment.
            pass
        else:
            deprecation(33114, 'method JoinFeature.is_functional; use is_present instead')
        for f in self._features:
            test = f.is_functional()
            if not test:
                return test
        return FeatureTestResult(self, True)
