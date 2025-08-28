# sage_setup: distribution = sagemath-environment
import os

from . import StaticFile


class JmolDataJar(StaticFile):
    r"""
    A :class:`~sage.features.Feature` which describes the presence of
    JmolData.jar in a few standard locations.

    EXAMPLES::

        sage: from sage.features.jmol import JmolDataJar
        sage: bool(JmolDataJar().is_present())  # needs jmol
        True
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.jmol import JmolDataJar
            sage: isinstance(JmolDataJar(), JmolDataJar)
            True
        """
        from sage.env import SAGE_SHARE, JMOL_DIR

        jmol_search_path = JMOL_DIR or (
                os.path.join(SAGE_SHARE, "sagemath", "jmol"),
                os.path.join(SAGE_SHARE, "jmol")
                )

        StaticFile.__init__(
            self, name='jmol',
            filename='JmolData.jar',
            search_path=jmol_search_path,
            spkg='jmol',
            type='optional',
            description="Java viewer for chemical structures in 3D")


def all_features():
    return [JmolDataJar()]
