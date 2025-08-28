# sage.doctest: needs sage.modules

#######################################################################
# Backward compatible unpickle functions
#######################################################################

from .quatalg.quaternion_algebra_element import (QuaternionAlgebraElement_generic,
                                                QuaternionAlgebraElement_rational_field,
                                                QuaternionAlgebraElement_number_field)


def unpickle_QuaternionAlgebraElement_generic_v0(*args):
    """
    EXAMPLES::

        sage: K.<X> = QQ[]
        sage: Q.<i,j,k> = QuaternionAlgebra(Frac(K), -5,-19); z = 2/3 + i*X - X^2*j + X^3*k
        sage: f, t = z.__reduce__()
        sage: import sage.algebras.quaternion_algebra_element
        sage: sage.algebras.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_generic_v0(*t)
        2/3 + X*i + (-X^2)*j + X^3*k
        sage: sage.algebras.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_generic_v0(*t) == z
        True
    """
    return QuaternionAlgebraElement_generic(*args)


def unpickle_QuaternionAlgebraElement_rational_field_v0(*args):
    """
    EXAMPLES::

        sage: Q.<i,j,k> = QuaternionAlgebra(-5,-19); a = 2/3 + i*5/7 - j*2/5 +19/2
        sage: f, t = a.__reduce__()
        sage: import sage.algebras.quaternion_algebra_element
        sage: sage.algebras.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_rational_field_v0(*t)
        61/6 + 5/7*i - 2/5*j
    """
    return QuaternionAlgebraElement_rational_field(*args)


def unpickle_QuaternionAlgebraElement_number_field_v0(*args):
    """
    EXAMPLES::

        sage: # needs sage.symbolic
        sage: K.<a> = QQ[2^(1/3)]; Q.<i,j,k> = QuaternionAlgebra(K, -3, a); z = i + j
        sage: f, t = z.__reduce__()
        sage: import sage.algebras.quaternion_algebra_element
        sage: sage.algebras.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_number_field_v0(*t)
        i + j
        sage: sage.algebras.quaternion_algebra_element.unpickle_QuaternionAlgebraElement_number_field_v0(*t) == z
        True
    """
    return QuaternionAlgebraElement_number_field(*args)
