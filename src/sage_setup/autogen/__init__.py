import os


def autogen_all():
    """
    Regenerate the automatically generated files of the Sage library.
    """
    from sage_setup.autogen import pari
    pari.rebuild()

    from sage_setup.autogen import interpreters
    interpreters.rebuild(os.path.join("sage", "ext", "interpreters"))
