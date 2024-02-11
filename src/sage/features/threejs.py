import os

from . import StaticFile


class Threejs(StaticFile):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of
    threejs-sage in a few standard locations.

    EXAMPLES::

        sage: from sage.features.threejs import Threejs
        sage: bool(Threejs().is_present())  # needs threejs
        True
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.threejs import Threejs
            sage: isinstance(Threejs(), Threejs)
            True
        """
        from sage.env import SAGE_SHARE, THREEJS_DIR

        version = self.required_version()

        threejs_search_path = THREEJS_DIR or (
            os.path.join(SAGE_SHARE, "jupyter", "nbextensions", "threejs-sage"),
            os.path.join(SAGE_SHARE, "sagemath", "threejs-sage"),
            os.path.join(SAGE_SHARE, "sage", "threejs"),
            os.path.join(SAGE_SHARE, "threejs-sage")
            )

        StaticFile.__init__(
            self, name="threejs",
            filename=os.path.join(version, "three.min.js"),
            spkg="threejs",
            type="standard",
            search_path=threejs_search_path,
            description="JavaScript library to display 3D graphics")

    def required_version(self):
        """
        Return the version of threejs that Sage requires.

        EXAMPLES::

            sage: from sage.features.threejs import Threejs
            sage: Threejs().required_version()
            'r...'
        """
        from sage.env import SAGE_EXTCODE

        filename = os.path.join(SAGE_EXTCODE, 'threejs', 'threejs-version.txt')

        with open(filename) as f:
            return f.read().strip()


def all_features():
    return [Threejs()]
