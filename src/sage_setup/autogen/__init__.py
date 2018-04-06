import os


def autogen_all():
    """
    Regenerate the automatically generated files of the Sage library.

    Return a list of sub-packages that should be appended to the list
    of packages built/installed by setup.py.
    """

    from . import interpreters
    interpreters.rebuild(os.path.join("sage", "ext", "interpreters"))

    return ['sage.ext.interpreters']
