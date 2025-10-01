r"""
Features for testing the presence of ``flatter``
"""

from . import Executable


class flatter(Executable):
    """
    A :class:`~sage.features.Feature` describing the presence of ``flatter``.

    EXAMPLES::

        sage: from sage.features.flatter import flatter
        sage: flatter().is_present()  # optional - flatter
        FeatureTestResult('flatter', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.flatter import flatter
            sage: isinstance(flatter(), flatter)
            True
        """
        Executable.__init__(self, "flatter", executable="flatter")


def all_features():
    return [flatter()]
