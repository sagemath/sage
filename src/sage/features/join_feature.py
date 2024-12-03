# sage_setup: distribution = sagemath-environment
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

    def __init__(self, name, features, spkg=None, url=None, description=None, type=None,
                 **kwds):
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
        if type is None:
            if any(f._spkg_type() == 'experimental' for f in features):
                type = 'experimental'
            elif any(f._spkg_type() == 'optional' for f in features):
                type = 'optional'
            else:
                type = 'standard'

        super().__init__(name, spkg=spkg, url=url, description=description, type=type, **kwds)
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

    def hide(self):
        r"""
        Hide this feature and all its joined features.

        EXAMPLES::

            sage: from sage.features.sagemath import sage__groups
            sage: f = sage__groups()
            sage: f.hide()
            sage: f._features[0].is_present()
            FeatureTestResult('sage.groups.perm_gps.permgroup', False)

            sage: f.require()
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: sage.groups is not available.
            Feature `sage.groups` is hidden.
            Use method `unhide` to make it available again.
        """
        for f in self._features:
            f.hide()
        super().hide()

    def unhide(self):
        r"""
        Revert what :meth:`hide` did.

        EXAMPLES::

            sage: from sage.features.sagemath import sage__groups
            sage: f = sage__groups()
            sage: f.hide()
            sage: f.is_present()
            FeatureTestResult('sage.groups', False)
            sage: f._features[0].is_present()
            FeatureTestResult('sage.groups.perm_gps.permgroup', False)

            sage: f.unhide()
            sage: f.is_present()    # optional sage.groups
            FeatureTestResult('sage.groups', True)
            sage: f._features[0].is_present() # optional sage.groups
            FeatureTestResult('sage.groups.perm_gps.permgroup', True)
        """
        for f in self._features:
            f.unhide()
        super().unhide()
