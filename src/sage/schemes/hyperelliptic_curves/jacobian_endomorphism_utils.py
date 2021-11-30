# -*- coding: utf-8 -*-
r"""
Some functions regarding geometric endomorphism rings of Jacobians of
hyperelliptic curves.

There are currently two main functions in this module, associated to the
Jacobian `J` of a genus 2 curve over the rational numbers `\QQ`:

- ``get_is_geom_field``
- ``is_geom_trivial_when_field``

``get_is_geom_field`` determines whether the geometric endomorphism
algebra \textup{End}_{\overline{\QQ}}(J) \otimes \QQ` is a field.

``is_geom_trivial_when_field`` determines whether the geometric endomorphism ring
`\textup{End}_{\overline{\QQ}}(J)` is equal to the integer ring `\ZZ`, assuming
the geometric endomorphism algebra is a field. If this is the case, one says
that `J` is ``generic``.

These functions are used in `jacobian_generic.py` in the following two
analogous methods associated to the Jacobian `J` of a genus 2 curve over the
rational numbers `\QQ`:

- ``geometric_endomorphism_algebra_is_field(J)``
- ``geometric_endomorphism_ring_is_ZZ(J)``

Both of these are important attributes of the Jacobian J.

The algorithms in the functions in this module are implementations of
Algorithms 4.10 and 4.15 from Lombardo's paper [Lom2019]_.

The following examples have been verified with the corresponding LMFDB entries.

EXAMPLES:

Here is an example of a generic Jacobian; the LMFDB label of the curve is
249.a.249.1::

    sage: R.<x> = QQ[]
    sage: f = x^6 + 2*x^3 + 4*x^2 + 4*x + 1
    sage: C = HyperellipticCurve(f)
    sage: A = C.jacobian()
    sage: A.geometric_endomorphism_ring_is_ZZ()
    True

Here is an example of a Jacobian whose endomorphism algebra is a field but not
the rational number field; the LMFDB label of the curve is 529.a.529.1::

    sage: f = x^6 - 4*x^5 + 2*x^4 + 2*x^3 + x^2 + 2*x + 1
    sage: C = HyperellipticCurve(f)
    sage: A = C.jacobian()
    sage: A.geometric_endomorphism_algebra_is_field()
    True
    sage: A.geometric_endomorphism_ring_is_ZZ()
    False

Here is an example of a Jacobian whose endomorphism algebra is not a field;
the LMFDB label of the curve is 169.a.169.1::

    sage: f = x^6 + 4*x^5 + 6*x^4 + 2*x^3 + x^2 + 2*x + 1
    sage: C = HyperellipticCurve(f)
    sage: A = C.jacobian()
    sage: A.geometric_endomorphism_algebra_is_field()
    False

.. WARNING:

    There is a very small chance that the algorithms return ``False`` for the
    two methods described above when in fact one or both of them are ``True``.
    In this case, as explained in the discussion immediately preceding
    Algorithm 4.15 of [Lom2019]_, this can be established by increasing the
    optional `B` parameter of the two methods. Mathematically, the algorithms
    give the correct answer only in the limit as `B \to \infty`, although in
    practice `B = 200` was sufficient to correctly verify every single entry
    in the LMFDB.

AUTHORS:

- Barinder S. Banwait and Davide Lombardo (2021-06-09): initial version

"""

# ****************************************************************************
#       Copyright (C) 2021 Barinder S. Banwait and Davide Lombardo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.all import QQ, ZZ, PolynomialRing, FiniteField, NumberField
from sage.misc.lazy_import import lazy_import
from sage.rings.fast_arith import prime_range
from sage.arith.all import gcd
lazy_import('sage.interfaces.genus2reduction', ['genus2reduction', 'Genus2reduction'])


def satisfies_coefficient_condition(g, p):
    """
    This is the coefficient condition in the definition of Omega_K'
    on page 912 of the published version of paper.

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves.jacobian_endomorphism_utils import satisfies_coefficient_condition
        sage: R.<x> = ZZ[]
        sage: f = x^4 + x^3 + 17*x^2 + 5*x
        sage: satisfies_coefficient_condition(f,17)
        False
        sage: f = x^4 + x^3 + 17*x^2 + 23*x + 23^2
        sage: satisfies_coefficient_condition(f,23)
        True
    """

    if g[0] != p**2:
        return False
    if g[3]*p != g[1]:
        return False
    if g[2]%p == 0:
        return False
    return True


def get_is_geom_field(f, C, bad_primes, B=200):
    r"""
    Determine whether the geometric endomorphism algebra is a field.

    This is Algorithm 4.10 in [Lom2019]_. The computation done here
    may allow one to immediately conclude that the geometric endomorphism
    ring is trivial (i.e. the integer ring); this information is output
    in a second boolean to avoid unnecessary subsequent computation.

    An additional optimisation comes from Part (2) of Theorem 4.8 in
    [Lom2019]_, from which we can conclude that the endomorphism ring
    is geometrically trivial, and from Proposition 4.7 in loc. cit. from
    which we can rule out potential QM.

    INPUT:

    - ``f`` -- a polynomial defining the hyperelliptic curve.

    - ``C`` -- the hyperelliptic curve.

    - ``bad_primes`` -- the list of odd primes of bad reduction.

    - ``B`` -- (default: 200) the bound which appears in the statement of
      the algorithm from [Lom2019]_

    OUTPUT:

    Pair of booleans (bool1, bool2). `bool1` indicates if the
    geometric endomorphism algebra is a field; `bool2` indicates if the
    geometric endomorphism algebra is the field of rational numbers.

    WARNING:

    There is a very small chance that this algorithm return ``False`` when in
    fact it is ``True``. In this case, as explained in the discussion
    immediately preceding Algorithm 4.15 of [Lom2019]_, this can be established
    by increasing the optional `B` parameter. Mathematically, this algorithm
    gives the correct answer only in the limit as `B \to \infty`, although in
    practice `B = 200` was sufficient to correctly verify every single entry
    in the LMFDB. However, strictly speaking, a ``False`` returned by this
    function is not provably ``False``.

    EXAMPLES:

    This is LMFDB curve 940693.a.960693.1::

        sage: from sage.schemes.hyperelliptic_curves.jacobian_endomorphism_utils import get_is_geom_field
        sage: R.<x> = QQ[]
        sage: f = 4*x^6 - 12*x^5 + 20*x^3 - 8*x^2 - 4*x + 1
        sage: C = HyperellipticCurve(f)
        sage: get_is_geom_field(f,C,[13,269])
        (False, False)

    This is LMFDB curve 3125.a.3125.1::

        sage: f = 4*x^5 + 1
        sage: C = HyperellipticCurve(f)
        sage: get_is_geom_field(f,C,[5])
        (True, False)

    This is LMFDB curve 277.a.277.2::

        sage: f = 4*x^6 - 36*x^4 + 56*x^3 - 76*x^2 + 44*x - 23
        sage: C = HyperellipticCurve(f)
        sage: get_is_geom_field(f,C,[277])
        (True, True)
    """

    if C.has_odd_degree_model():
        C_odd = C.odd_degree_model()
        f_odd, h_odd = C_odd.hyperelliptic_polynomials()
        # if f was odd to begin with, then f_odd = f
        assert f_odd.degree() == 5

        if (4*f_odd + h_odd**2).degree() == 5:
            f_new = 4*f_odd + h_odd**2
            if f_new.is_irreducible():
                # i.e. the Jacobian is geometrically simple
                f_disc_odd_prime_exponents = [v for _,v in f_new.discriminant().prime_to_S_part([ZZ(2)]).factor()]
                if 1 in f_disc_odd_prime_exponents:
                    return (True, True)  # Theorem 4.8 (2)
                # At this point we are in the situation of Algorithm 4.10
                # Step 1, so either the geometric endomorphism algebra is a
                # field or it is a quaternion algebra. This latter case implies
                # that the Jacobian is the square of an elliptic curve modulo
                # a prime p where f is also irreducible (which exist by
                # Chebotarev density). This contradicts Prop 4.7, hence we can
                # conclude as follows.
                return (True, False)

    if f.is_irreducible():
        assert f.degree() == 6  # else we should already have exited by now
        G = f.galois_group()
        if G.order() in [360, 720]:
            return (True, True)  # Algorithm 4.10 Step 2

    R = PolynomialRing(ZZ,2,"xv")
    x,v = R.gens()
    T = PolynomialRing(QQ,'v')
    g = v - x**12

    for p in prime_range(3,B):
        if p not in bad_primes:
            fp = C.change_ring(FiniteField(p)).frobenius_polynomial()

            # This defines the polynomial f_v**[12] from the paper
            fp12 = T(R(fp).resultant(g))

            if fp12.is_irreducible():
                # i.e. the Jacobian is geometrically simple
                f_disc_odd_prime_exponents = [v for _,v in f.discriminant().prime_to_S_part([ZZ(2)]).factor()]
                if 1 in f_disc_odd_prime_exponents:
                    return (True, True)  # Theorem 4.8 (2)
                return (True, False) # Algorithm 4.10 Step 3 plus Prop 4.7 as above
    return (False, False)


def is_geom_trivial_when_field(C, bad_primes, B=200):
    r"""
    Determine if the geometric endomorphism ring is trivial assuming the
    geometric endomorphism algebra is a field.

    This is Algorithm 4.15 in [Lom2019]_.

    INPUT:

    - ``C`` -- the hyperelliptic curve.

    - ``bad_primes`` -- the list of odd primes of bad reduction.

    - ``B`` -- (default: 200) the bound which appears in the statement of
      the algorithm from [Lom2019]_

    OUTPUT:

    Boolean indicating whether or not the geometric endomorphism
    algebra is the field of rational numbers.

    WARNING:

    There is a very small chance that this algorithm returns ``False`` when in
    fact it is ``True``. In this case, as explained in the discussion
    immediately preceding Algorithm 4.15 of [Lom2019]_, this can be established
    by increasing the optional `B` parameter. Mathematically, this algorithm
    gives the correct answer only in the limit as `B \to \infty`, although in
    practice `B = 200` was sufficient to correctly verify every single entry
    in the LMFDB. However, strictly speaking, a ``False`` returned by this
    function is not provably ``False``.

    EXAMPLES:

    This is LMFDB curve 461.a.461.2::

        sage: from sage.schemes.hyperelliptic_curves.jacobian_endomorphism_utils import is_geom_trivial_when_field
        sage: R.<x> = QQ[]
        sage: f = 4*x^5 - 4*x^4 - 156*x^3 + 40*x^2 + 1088*x - 1223
        sage: C = HyperellipticCurve(f)
        sage: is_geom_trivial_when_field(C,[461])
        True

    This is LMFDB curve 4489.a.4489.1::

        sage: f = x^6 + 4*x^5 + 2*x^4 + 2*x^3 + x^2 - 2*x + 1
        sage: C = HyperellipticCurve(f)
        sage: is_geom_trivial_when_field(C,[67])
        False
    """

    running_gcd = 0
    R = PolynomialRing(ZZ,2,"xv")
    x,v = R.gens()
    T = PolynomialRing(QQ,'v')
    g = v - x**4

    for p in prime_range(3,B):
        if p not in bad_primes:
            Cp = C.change_ring(FiniteField(p))
            fp = Cp.frobenius_polynomial()
            if satisfies_coefficient_condition(fp, p):
                # This defines the polynomial f_v**[4] from the paper
                fp4 = T(R(fp).resultant(g))
                if fp4.is_irreducible():
                    running_gcd = gcd(running_gcd, NumberField(fp,'a').discriminant())
                    if running_gcd <= 24:
                        return True
    return False
