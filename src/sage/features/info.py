# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of ``info``, from GNU Info
"""

from . import Executable


class Info(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :ref:`info <spkg_info>`.

    EXAMPLES::

        sage: from sage.features.info import Info
        sage: Info()
        Feature('info')
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.info import Info
            sage: isinstance(Info(), Info)
            True
        """
        Executable.__init__(self, 'info', executable='info',
                            spkg='info', type='standard')


def all_features():
    return [Info()]
