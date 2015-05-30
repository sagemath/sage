"""
Deprecated two by two matrices over the integers.

See :trac:`17824` for more informations.
"""

def MatrixSpace_ZZ_2x2():
    """
    Return the space of 2x2 integer matrices.

    See :trac:`17824` for more informations.

    EXAMPLES::

        sage: from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
        sage: M = MatrixSpace_ZZ_2x2()
        doctest:...: DeprecationWarning: MatrixSpace_ZZ_2x2 is deprecated.
        Please use MatrixSpace(ZZ,2) instead See http://trac.sagemath.org/17824
        for details.
        sage: M
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: M is MatrixSpace_ZZ_2x2()
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(17824, 'MatrixSpace_ZZ_2x2 is deprecated. Please use MatrixSpace(ZZ,2) instead')

    from sage.matrix.matrix_space import MatrixSpace
    from sage.rings.integer_ring import ZZ
    return MatrixSpace(ZZ,2)
