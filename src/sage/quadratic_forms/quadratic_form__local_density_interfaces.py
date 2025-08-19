"""
Local Density Interfaces
"""
# // This is needed in the filter for primitivity...
# #include "../max-min.h"

from sage.arith.misc import valuation
from sage.rings.rational_field import QQ


def local_density(self, p, m):
    """
    Return the local density.

    .. NOTE::

        This screens for imprimitive forms, and puts the quadratic
        form in local normal form, which is a *requirement* of the
        routines performing the computations!

    INPUT:

    - ``p`` -- a prime number > 0
    - ``m`` -- integer

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])   # NOTE: This is already in local normal form for *all* primes p!
        sage: Q.local_density(p=2, m=1)
        1
        sage: Q.local_density(p=3, m=1)
        8/9
        sage: Q.local_density(p=5, m=1)
        24/25
        sage: Q.local_density(p=7, m=1)
        48/49
        sage: Q.local_density(p=11, m=1)
        120/121
    """
    n = self.dim()
    if n == 0:
        raise TypeError("we do not currently handle 0-dim'l forms")

    # Find the local normal form and p-scale of Q     --  Note: This uses the valuation ordering of local_normal_form.
    #                                                     TO DO:  Write a separate p-scale and p-norm routines!
    Q_local = self.local_normal_form(p)
    if n == 1:
        p_valuation = valuation(Q_local[0, 0], p)
    else:
        p_valuation = min(valuation(Q_local[0, 0], p),
                          valuation(Q_local[0, 1], p))

    # If m is less p-divisible than the matrix, return zero
    if ((m != 0) and (valuation(m, p) < p_valuation)):   # Note: The (m != 0) condition protects taking the valuation of zero.
        return QQ(0)

    # If the form is imprimitive, rescale it and call the local density routine
    p_adjustment = QQ(1) / p**p_valuation
    m_prim = QQ(m) / p**p_valuation
    Q_prim = Q_local.scale_by_factor(p_adjustment)

    # Return the densities for the reduced problem
    return Q_prim.local_density_congruence(p, m_prim)


def local_primitive_density(self, p, m):
    """
    Return the local primitive density -- should be called by the user.

    NOTE: This screens for imprimitive forms, and puts the
    quadratic form in local normal form, which is a *requirement* of
    the routines performing the computations!

    INPUT:

    - ``p`` -- a prime number > 0
    - ``m`` -- integer

    OUTPUT: a rational number

    EXAMPLES::

        sage: Q = QuadraticForm(ZZ, 4, range(10))
        sage: Q[0,0] = 5
        sage: Q[1,1] = 10
        sage: Q[2,2] = 15
        sage: Q[3,3] = 20
        sage: Q
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 5 1 2 3 ]
        [ * 10 5 6 ]
        [ * * 15 8 ]
        [ * * * 20 ]
        sage: Q.theta_series(20)                                                        # needs sage.libs.pari
        1 + 2*q^5 + 2*q^10 + 2*q^14 + 2*q^15 + 2*q^16 + 2*q^18 + O(q^20)
        sage: Q.local_normal_form(2)
        Quadratic form in 4 variables over Integer Ring with coefficients:
        [ 0 1 0 0 ]
        [ * 0 0 0 ]
        [ * * 0 1 ]
        [ * * * 0 ]

        sage: Q.local_primitive_density(2, 1)
        3/4
        sage: Q.local_primitive_density(5, 1)
        24/25

        sage: Q.local_primitive_density(2, 5)
        3/4
        sage: Q.local_density(2, 5)
        3/4
    """
    n = self.dim()
    if n == 0:
        raise TypeError("we do not currently handle 0-dim'l forms")

    # Find the local normal form and p-scale of Q     --  Note: This uses the valuation ordering of local_normal_form.
    #                                                     TO DO:  Write a separate p-scale and p-norm routines!
    Q_local = self.local_normal_form(p)
    if n == 1:
        p_valuation = valuation(Q_local[0, 0], p)
    else:
        p_valuation = min(valuation(Q_local[0, 0], p),
                          valuation(Q_local[0, 1], p))

    # If m is less p-divisible than the matrix, return zero
    if m != 0 and valuation(m, p) < p_valuation:   # Note: The (m != 0) condition protects taking the valuation of zero.
        return QQ.zero()

    # If the form is imprimitive, rescale it and call the local density routine
    p_adjustment = QQ.one() / p**p_valuation
    m_prim = QQ(m) / p**p_valuation
    Q_prim = Q_local.scale_by_factor(p_adjustment)

    # Return the densities for the reduced problem
    return Q_prim.local_primitive_density_congruence(p, m_prim)
