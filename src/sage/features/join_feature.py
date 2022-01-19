r"""
Join features
"""

from . import Feature, FeatureTestResult


class JoinFeature(Feature):
    r"""
    Join of several :class:`~sage.features.Feature` instances.

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

        EXAMPLES::

            sage: from sage.features.latte import Latte
            sage: Latte().is_functional()  # optional - latte_int
            FeatureTestResult('latte_int', True)
        """
        for f in self._features:
            test = f.is_functional()
            if not test:
                return test
        return FeatureTestResult(self, True)
