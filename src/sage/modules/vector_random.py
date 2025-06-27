from sage.modules.free_module import FreeModule

def vector_random(R, n, *args, **kwargs):
    """
    Return a random vector over ring R of length n.

    Mirrors matrix.random for API symmetry.
    """
    return FreeModule(R, n).random_element(*args, **kwargs)
