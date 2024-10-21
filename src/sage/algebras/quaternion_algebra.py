# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.modules

############################################################
# Backwards compatible unpickling
############################################################

def unpickle_QuaternionAlgebra_v0(*key):
    """
    The `0`-th version of pickling for quaternion algebras.

    EXAMPLES::

        sage: t = (QQ, -5, -19, ('i', 'j', 'k'))
        sage: import sage.algebras.quaternion_algebra
        sage: sage.algebras.quaternion_algebra.unpickle_QuaternionAlgebra_v0(*t)
        Quaternion Algebra (-5, -19) with base ring Rational Field
    """
    from .quatalg.quaternion_algebra import QuaternionAlgebra
    return QuaternionAlgebra(*key)
