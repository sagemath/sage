import os
from sage.env import SAGE_SRC

def autogen_all():
    """
    Regenerate the automatically generated files of the Sage library.
    """
    import pari
    pari.rebuild()

    import interpreters
    interpreters.rebuild(os.path.join(SAGE_SRC, "sage", "ext", "interpreters"))
