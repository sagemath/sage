import os

from interpreters.internal import rebuild


def autogen_all():
    """
    Regenerate the automatically generated files of the Sage library.

    Return a list of sub-packages that should be appended to the list
    of packages built/installed by setup.py.
    """
    from sage.env import SAGE_SRC
    rebuild(os.path.join(SAGE_SRC, "sage", "ext", "interpreters"))

    return ['sage.ext.interpreters']
