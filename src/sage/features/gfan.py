r"""
Features for testing the presence of ``gfan``
"""

from . import Executable


class GfanExecutable(Executable):
    r"""
    A :class:`~sage.features.Feature` for the gfan executables.
    """
    def __init__(self, cmd=None):
        r"""
        TESTS::

            sage: from sage.features.gfan import GfanExecutable
            sage: isinstance(GfanExecutable('groebnercone'), GfanExecutable)
            True
        """
        if cmd is None:
            name = "gfan"
        else:
            name = f"gfan_{cmd}"
        Executable.__init__(self, name, executable=name, spkg="gfan")


def all_features():
    return [GfanExecutable()]
