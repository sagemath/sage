r"""
Feature for testing if the Internet is available
"""

from . import Feature, FeatureTestResult


class Internet(Feature):
    r"""
    A :class:`~sage.features.Feature` describing if Internet is available.

    Failure of connecting to the site "https://www.sagemath.org" within a second
    is regarded as internet being not available.

    EXAMPLES::

        sage: from sage.features.internet import Internet
        sage: Internet()
        Feature('internet')
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.internet import Internet
            sage: Internet() is Internet()
            True
        """
        Feature.__init__(self, 'internet')

    def _is_present(self):
        r"""
        Test whether Internet is available.

        EXAMPLES::

            sage: from sage.features.internet import Internet
            sage: Internet()._is_present()  # random, optional - internet
            FeatureTestResult('internet', True)
        """
        import urllib.error
        from urllib.request import Request, urlopen
        from ssl import SSLContext

        req = Request("https://www.sagemath.org", headers={"User-Agent": "sage-doctest"})
        try:
            urlopen(req, timeout=1, context=SSLContext())
            return FeatureTestResult(self, True)
        except urllib.error.URLError:
            return FeatureTestResult(self, False)


def all_features():
    return [Internet()]
