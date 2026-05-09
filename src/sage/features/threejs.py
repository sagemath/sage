from pathlib import Path

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
        from sage.env import sage_data_paths, THREEJS_DIR

        threejs_search_path = THREEJS_DIR or list(sage_data_paths('jupyter/nbextensions/threejs-sage')) + list(sage_data_paths('threejs-sage')) + list(sage_data_paths('threejs'))

        try:
            version = self.required_version()
            filename = Path(version) / "three.min.js"
        except FileNotFoundError:
            filename = 'unknown'

        StaticFile.__init__(
            self, name='threejs',
            filename=filename,
            spkg='threejs',
            type='standard',
            search_path=threejs_search_path,
            description="JavaScript library to display 3D graphics")

    def required_version(self):
        """
        Return the version of threejs that Sage requires.

        Defining what version is required is delegated to the distribution package
        that provides the file ``threejs-version.txt`` in :mod:`sage.ext_data.threejs`.

        If the file is not provided, :exc:`FileNotFoundError` is raised.

        EXAMPLES::

            sage: from sage.features.threejs import Threejs
            sage: Threejs().required_version()
            'r...'
        """
        from sage.env import SAGE_EXTCODE

        filename = Path(SAGE_EXTCODE) / 'threejs' / 'threejs-version.txt'

        with open(filename) as f:
            return f.read().strip()


def all_features():
    return [Threejs()]
