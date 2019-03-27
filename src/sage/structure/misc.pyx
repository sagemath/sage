"""
Miscellaneous utilities
"""

def is_extension_type(cls):
    """
    INPUT:

    - cls: a class

    Tests whether cls is an extension type (int, list, cython compiled classes, ...)

    EXAMPLES::

        sage: from sage.structure.parent import is_extension_type
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
