# sage_setup: distribution = sagemath-pari
# sage.doctest: needs sage.libs.pari sage.modules
"""
Siegel Products
"""

# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import (bernoulli,
                             fundamental_discriminant,
                             kronecker as kronecker_symbol,
                             prime_divisors)
from sage.misc.functional import sqrt
from sage.misc.verbose import verbose
from sage.quadratic_forms.special_values import QuadraticBernoulliNumber
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

# /*!  \brief Computes the product of all local densities for comparison with independently computed Eisenstein coefficients.
# *
# *  \todo We fixed the generic factors to compensate for using the matrix of 2Q, but we need to document this better! =)
# */

# /////////////////////////////////////////////////////////////////
# ///
# /////////////////////////////////////////////////////////////////

# mpq_class Matrix_mpz::siegel_product(mpz_class u) const {


def siegel_product(self, u):
    """
    Compute the infinite product of local densities of the quadratic
    form for the number `u`.

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1,1])
        sage: Q.theta_series(11)
        1 + 8*q + 24*q^2 + 32*q^3 + 24*q^4 + 48*q^5 + 96*q^6 + 64*q^7 + 24*q^8 + 104*q^9 + 144*q^10 + O(q^11)

        sage: Q.siegel_product(1)
        8
        sage: Q.siegel_product(2)      # This one is wrong -- expect 24, and the higher powers of 2 don't work... =(
        24
        sage: Q.siegel_product(3)
        32
        sage: Q.siegel_product(5)
        48
        sage: Q.siegel_product(6)
        96
        sage: Q.siegel_product(7)
        64
        sage: Q.siegel_product(9)
        104

        sage: Q.local_density(2,1)
        1
        sage: M = 4; len([v  for v in mrange([M,M,M,M])  if Q(v) % M == 1]) / M^3
        1
        sage: M = 16; len([v  for v in mrange([M,M,M,M])  if Q(v) % M == 1]) / M^3  # long time (2s on sage.math, 2014)
        1

        sage: Q.local_density(2,2)
        3/2
        sage: M = 4; len([v  for v in mrange([M,M,M,M])  if Q(v) % M == 2]) / M^3
        3/2
        sage: M = 16; len([v  for v in mrange([M,M,M,M])  if Q(v) % M == 2]) / M^3  # long time (2s on sage.math, 2014)
        3/2

    TESTS::

        sage: [1] + [Q.siegel_product(ZZ(a))  for a in range(1,11)] == Q.theta_series(11).list()  # long time (2s on sage.math, 2014)
        True
    """
    # Protect u (since it fails often if it's an just an int!)
    u = ZZ(u)

    n = self.dim()
    d = self.det()       # ??? Warning: This is a factor of 2^n larger than it should be!

    # DIAGNOSTIC
    verbose("n = " + str(n))
    verbose("d = " + str(d))
    verbose("In siegel_product:  d = " + str(d) + "\n")

    # Product of "bad" places to omit
    S = 2 * d * u

    # DIAGNOSTIC
    verbose("siegel_product Break 1. \n")
    verbose(" u = " + str(u) + "\n")

    # Make the odd generic factors
    if n % 2:
        m = (n - 1) // 2
        d1 = fundamental_discriminant(((-1)**m) * 2*d * u)     # Replaced d by 2d here to compensate for the determinant
        f = abs(d1)
        # gaining an odd power of 2 by using the matrix of 2Q instead
        # of the matrix of Q.
        #  --> Old d1 = CoreDiscriminant((mpz_class(-1)^m) * d * u);

        # Make the ratio of factorials factor: [(2m)! / m!] * prod_{i=1}^m (2*i-1)
        factor1 = 1
        for i in range(1, m+1):
            factor1 *= 2*i - 1
        for i in range(m+1, 2*m + 1):
            factor1 *= i

        genericfactor = factor1 * ((u / f) ** m) \
            * QQ(sqrt((2 ** n) * f) / (u * d)) \
            * abs(QuadraticBernoulliNumber(m, d1) / bernoulli(2*m))

    # DIAGNOSTIC
    verbose("siegel_product Break 2. \n")

    # Make the even generic factor
    if ((n % 2) == 0):
        m = n // 2
        d1 = fundamental_discriminant(((-1)**m) * d)
        f = abs(d1)

        # DIAGNOSTIC
        # cout << " mpz_class(-1)^m = " << (mpz_class(-1)^m) << " and d = " << d << endl;
        # cout << " f = " << f << " and d1 = " << d1 << endl;

        genericfactor = m / QQ(sqrt(f*d)) \
            * ((u/2) ** (m-1)) * (f ** m) \
            / abs(QuadraticBernoulliNumber(m, d1)) \
            * (2 ** m)                                               # This last factor compensates for using the matrix of 2*Q

    # return genericfactor

    # Omit the generic factors in S and compute them separately
    omit = 1
    include = 1

    S_divisors = prime_divisors(S)

    # DIAGNOSTIC
    # cout << "\n S is " << S << endl;
    # cout << " The Prime divisors of S are :";
    # PrintV(S_divisors);

    for p in S_divisors:
        Q_normal = self.local_normal_form(p)

        # DIAGNOSTIC
        verbose(" p = " + str(p) + " and its Kronecker symbol (d1/p) = (" + str(d1) + "/" + str(p) + ") is " + str(kronecker_symbol(d1, p)) + "\n")

        omit *= 1 / (1 - (kronecker_symbol(d1, p) / (p**m)))

        # DIAGNOSTIC
        verbose(" omit = " + str(omit) + "\n")
        verbose(" Q_normal is \n" + str(Q_normal) + "\n")
        verbose(" Q_normal = \n" + str(Q_normal))
        verbose(" p = " + str(p) + "\n")
        verbose(" u = " + str(u) + "\n")
        verbose(" include = " + str(include) + "\n")

        include *= Q_normal.local_density(p, u)

        # DIAGNOSTIC
        # cout << " Including the p = " << p << " factor: " << local_density(Q_normal, p, u) << endl;

        # DIAGNOSTIC
        verbose("    ---  Exiting loop \n")

    # // ****************  Important *******************
    # // Additional fix (only included for n=4) to deal
    # // with the power of 2 introduced at the real place
    # // by working with Q instead of 2*Q.  This needs to
    # // be done for all other n as well...
    # /*
    # if (n==4)
    #  genericfactor = 4 * genericfactor;
    # */

    # DIAGNOSTIC
    # cout << endl;
    # cout << " generic factor = " << genericfactor << endl;
    # cout << " omit = " << omit << endl;
    # cout << " include = " << include << endl;
    # cout << endl;

    # DIAGNOSTIC
    # //  cout << "siegel_product Break 3. " << endl;

    # Return the final factor (and divide by 2 if n=2)
    if n == 2:
        return genericfactor * omit * include / 2
    return genericfactor * omit * include
