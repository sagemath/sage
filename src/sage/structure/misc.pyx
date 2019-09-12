"""
Miscellaneous utilities
"""

from sage.misc.superseded import deprecation
deprecation(27099, "the module sage.structure.misc is deprecated")

def is_extension_type(cls):
    """
    Return whether or not instances of ``cls`` have a ``__dict__``.

    This is deprecated as there should not be any use case for it.

    INPUT:

    - ``cls`` -- a class

    EXAMPLES::

        sage: from sage.structure.misc import is_extension_type
        doctest:...: DeprecationWarning: the module sage.structure.misc is deprecated
        See https://trac.sagemath.org/27099 for details.
        sage: is_extension_type(int)
        True
        sage: is_extension_type(list)
        True
        sage: is_extension_type(ZZ.__class__)
        True
        sage: is_extension_type(QQ.__class__)
        False
    """
    # Robert B claims that this should be robust
    try:
        return cls.__dictoffset__ == 0
    except AttributeError:
        pass
    return False
