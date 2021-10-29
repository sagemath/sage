r"""
Reconstruction of Algebraic Forms

This module reconstructs algebraic forms from the values of their
invariants. Given a set of (classical) invariants, it returns a
form that attains this values as invariants (up to scaling).

AUTHORS:

- Jesper Noordsij (2018-06): initial version
"""

# ****************************************************************************
#     Copyright (C) 2018 Jesper Noordsij <jesper.noordsij@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


def binary_quadratic_coefficients_from_invariants(discriminant, invariant_choice='default'):
    """
    Reconstruct a binary quadratic from the value of its discriminant.

    INPUT:

    - ``discriminant`` -- The value of the discriminant of the
      binary quadratic.

    - ``invariant_choice`` -- The type of invariants provided. The accepted
      options are ``'discriminant'`` and ``'default'``, which are the same. No
      other options are implemented.

    OUTPUT:

    A set of coefficients of a binary quadratic, whose discriminant
    is equal to the given ``discriminant`` up to a scaling.

    EXAMPLES::

        sage: from sage.rings.invariants.reconstruction import binary_quadratic_coefficients_from_invariants
        sage: quadratic = invariant_theory.binary_form_from_invariants(2, [24]) # indirect doctest
        sage: quadratic
        Binary quadratic with coefficients (1, -6, 0)
        sage: quadratic.discriminant()
        24
        sage: binary_quadratic_coefficients_from_invariants(0)
        (1, 0, 0)
    """
    if invariant_choice not in ['default', 'discriminant']:
        raise ValueError('unknown choice of invariants {} for a binary '
                         'quadratic'.format(invariant_choice))
    if discriminant == 0:
        return (1, 0, 0)
    else:
        try:
            return (1, 0, -discriminant/4)
        except ZeroDivisionError:
            return (0, 1, 0)


def binary_cubic_coefficients_from_invariants(discriminant, invariant_choice='default'):
    """
    Reconstruct a binary cubic from the value of its discriminant.

    INPUT:

    - ``discriminant`` -- The value of the discriminant of the
      binary cubic.

    - ``invariant_choice`` -- The type of invariants provided. The accepted
      options are ``'discriminant'`` and ``'default'``, which are the same. No
      other options are implemented.

    OUTPUT:

    A set of coefficients of a binary cubic, whose discriminant
    is equal to the given ``discriminant`` up to a scaling.

    EXAMPLES::

        sage: from sage.rings.invariants.reconstruction import binary_cubic_coefficients_from_invariants
        sage: coeffs = binary_cubic_coefficients_from_invariants(1)
        sage: coeffs
        (0, 1, -1, 0)
        sage: R.<x> = QQ[]
        sage: R(coeffs).discriminant()
        1

    The two non-equivalent cubics `x^3` and `x^2*z` with discriminant 0 can't
    be distinguished based on their discriminant, hence an error is raised::

        sage: binary_cubic_coefficients_from_invariants(0)
        Traceback (most recent call last):
        ...
        ValueError: no unique reconstruction possible for binary cubics with a double root
    """
    if invariant_choice not in ['default', 'discriminant']:
        raise ValueError('unknown choice of invariants {} for a binary cubic'
                         .format(invariant_choice))
    if discriminant == 0:
        raise ValueError('no unique reconstruction possible for binary '
                         'cubics with a double root')
    else:
        return (0, 1, -1, 0)


def binary_quintic_coefficients_from_invariants(invariants, K=None, invariant_choice='default', scaling='none'):
    r"""
    Reconstruct a binary quintic from the values of its (Clebsch) invariants.

    INPUT:

    - ``invariants`` -- A list or tuple of values of the three or four
      invariants. The default option requires the Clebsch invariants `A`, `B`,
      `C` and `R` of the binary quintic.

    - ``K`` -- The field over which the quintic is defined.

    - ``invariant_choice`` -- The type of invariants provided. The accepted
      options are ``'clebsch'`` and ``'default'``, which are the same. No
      other options are implemented.

    - ``scaling`` -- How the coefficients should be scaled. The accepted
      values are ``'none'`` for no scaling, ``'normalized'`` to scale in such
      a way that the resulting coefficients are independent of the scaling of
      the input invariants and ``'coprime'`` which scales the input invariants
      by dividing them by their gcd.

    OUTPUT:

    A set of coefficients of a binary quintic, whose invariants are equal to
    the given ``invariants`` up to a scaling.

    EXAMPLES:

    First we check the general case, where the invariant `M` is non-zero::

        sage: R.<x0, x1> = QQ[]
        sage: p = 3*x1^5 + 6*x1^4*x0 + 3*x1^3*x0^2 + 4*x1^2*x0^3 - 5*x1*x0^4 + 4*x0^5
        sage: quintic = invariant_theory.binary_quintic(p, x0, x1)
        sage: invs = quintic.clebsch_invariants(as_tuple=True)
        sage: reconstructed = invariant_theory.binary_form_from_invariants(5, invs, variables=quintic.variables()) # indirect doctest
        sage: reconstructed
        Binary quintic with coefficients (9592267437341790539005557/244140625000000,
        2149296928207625556323004064707/610351562500000000,
        11149651890347700974453304786783/76293945312500000,
        122650775751894638395648891202734239/47683715820312500000,
        323996630945706528474286334593218447/11920928955078125000,
        1504506503644608395841632538558481466127/14901161193847656250000)

    We can see that the invariants of the reconstructed form match the ones of
    the original form by scaling the invariants `B` and `C`::

        sage: scale = invs[0]/reconstructed.A_invariant()
        sage: invs[1] == reconstructed.B_invariant()*scale^2
        True
        sage: invs[2] == reconstructed.C_invariant()*scale^3
        True

    If we compare the form obtained by this reconstruction to the one found by
    letting the covariants `\alpha` and `\beta` be the coordinates of the form,
    we find the forms are the same up to a power of the determinant of `\alpha`
    and `\beta`::

        sage: alpha = quintic.alpha_covariant()
        sage: beta = quintic.beta_covariant()
        sage: g = matrix([[alpha(x0=1,x1=0),alpha(x0=0,x1=1)],[beta(x0=1,x1=0),beta(x0=0,x1=1)]])^-1
        sage: transformed = tuple([g.determinant()^-5*x for x in quintic.transformed(g).coeffs()])
        sage: transformed == reconstructed.coeffs()
        True

    This can also be seen by computing the `\alpha` covariant of the obtained
    form::

        sage: reconstructed.alpha_covariant().coefficient(x1)
        0
        sage: reconstructed.alpha_covariant().coefficient(x0) != 0
        True

    If the invariant `M` vanishes, then the coefficients are computed in a
    different way::

        sage: [A,B,C] = [3,1,2]
        sage: M = 2*A*B - 3*C
        sage: M
        0
        sage: from sage.rings.invariants.reconstruction import binary_quintic_coefficients_from_invariants
        sage: reconstructed = binary_quintic_coefficients_from_invariants([A,B,C])
        sage: reconstructed
        (-66741943359375/2097152,
         -125141143798828125/134217728,
         0,
         52793920040130615234375/34359738368,
         19797720015048980712890625/1099511627776,
         -4454487003386020660400390625/17592186044416)
        sage: newform = sum([ reconstructed[i]*x0^i*x1^(5-i) for i in range(6) ])
        sage: newquintic = invariant_theory.binary_quintic(newform, x0, x1)
        sage: scale = 3/newquintic.A_invariant()
        sage: [3, newquintic.B_invariant()*scale^2, newquintic.C_invariant()*scale^3]
        [3, 1, 2]

    Several special cases::

        sage: quintic = invariant_theory.binary_quintic(x0^5 - x1^5, x0, x1)
        sage: invs = quintic.clebsch_invariants(as_tuple=True)
        sage: binary_quintic_coefficients_from_invariants(invs)
        (1, 0, 0, 0, 0, 1)
        sage: quintic = invariant_theory.binary_quintic(x0*x1*(x0^3-x1^3), x0, x1)
        sage: invs = quintic.clebsch_invariants(as_tuple=True)
        sage: binary_quintic_coefficients_from_invariants(invs)
        (0, 1, 0, 0, 1, 0)
        sage: quintic = invariant_theory.binary_quintic(x0^5 + 10*x0^3*x1^2 - 15*x0*x1^4, x0, x1)
        sage: invs = quintic.clebsch_invariants(as_tuple=True)
        sage: binary_quintic_coefficients_from_invariants(invs)
        (1, 0, 10, 0, -15, 0)
        sage: quintic = invariant_theory.binary_quintic(x0^2*(x0^3 + x1^3), x0, x1)
        sage: invs = quintic.clebsch_invariants(as_tuple=True)
        sage: binary_quintic_coefficients_from_invariants(invs)
        (1, 0, 0, 1, 0, 0)
        sage: quintic = invariant_theory.binary_quintic(x0*(x0^4 + x1^4), x0, x1)
        sage: invs = quintic.clebsch_invariants(as_tuple=True)
        sage: binary_quintic_coefficients_from_invariants(invs)
        (1, 0, 0, 0, 1, 0)

    For fields of characteristic 2, 3 or 5, there is no reconstruction
    implemented. This is part of :trac:`26786`.::

        sage: binary_quintic_coefficients_from_invariants([3,1,2], K=GF(5))
        Traceback (most recent call last):
        ...
        NotImplementedError: no reconstruction of binary quintics implemented for fields of characteristic 2, 3 or 5

    TESTS::

        sage: from sage.rings.invariants.reconstruction import binary_quintic_coefficients_from_invariants
        sage: binary_quintic_coefficients_from_invariants([1,2,3], scaling='unknown')
        Traceback (most recent call last):
        ...
        ValueError: unknown scaling option 'unknown'
    """
    if invariant_choice not in ['default', 'clebsch']:
        raise ValueError('unknown choice of invariants {} for a binary quintic'
                         .format(invariant_choice))
    if scaling not in ['none', 'normalized', 'coprime']:
        raise ValueError("unknown scaling option '%s'" % scaling)
    if scaling == 'coprime':
        if len(invariants) == 3:
            invariants = _reduce_invariants(invariants, [1,2,3])
        elif len(invariants) == 4:
            invariants = _reduce_invariants(invariants, [2,4,6,9])
    A, B, C = invariants[0:3]
    if K is None:
        from sage.rings.fraction_field import FractionField
        K = FractionField(A.parent())
    if K.characteristic() in [2, 3, 5]:
        raise NotImplementedError('no reconstruction of binary quintics '
                          'implemented for fields of characteristic 2, 3 or 5')
    M = 2*A*B - 3*C
    N = K(2)**-1 * (A*C-B**2)
    R2 = -K(2)**-1 * (A*N**2-2*B*M*N+C*M**2)
    scale = [1,1,1,1,1,1]
    from sage.functions.all import binomial
    from sage.misc.functional import sqrt
    if len(invariants) == 3:
        if R2.is_square():
            R = sqrt(R2)
        else:
            # if R2 is not a square, we scale the invariants in a suitable way
            # so that the 'new' R2 is a square
            [A, B, C] = [R2*A, R2**2*B, R2**3*C]
            [M, N] = [R2**3*M, R2**4*N]
            R = R2**5
    elif len(invariants) == 4:
        if invariants[3]**2 != R2:
            raise ValueError('provided invariants do not satisfy the syzygy '
                             'for Clebsch invariants of a binary quintic')
        R = invariants[3]
    else:
        raise ValueError('incorrect number of invariants provided, this '
                         'method requires 3 or 4 invariants')
    if M == 0:
        if N == 0:
            if A == 0:
                raise ValueError('no unique reconstruction possible for '
                                 'quintics with a treefold linear factor')
            else:
                if B == 0:
                    return (1,0,0,0,0,1)
                else:
                    return (0,1,0,0,1,0)
        else:
            # case corresponding to using alpha and gamma as coordinates
            if A == 0:
                return (1,0,0,0,1,0)
            else:
                if scaling == 'normalized':
                    # scaling z by (R/A**3)
                    scale = [ (-N)**-5*A**6*(R/A**3)**i for i in range(6) ]
                D = -N
                Delta = C
                a = [0]
                a.append((2*K(3)**-1*A**2-B)*N*B*K(2)**-1 - N**2*K(2)**-1)
                B0 = 2*K(3)**-1*A*R
                B1 = A*N*B*K(3)**-1
                C0 = 2*K(3)**-1*R
                C1 = B*N
    else:
        # case corresponding to using alpha and beta as coordinates
        if R == 0:
            if A == 0:
                return (1,0,10,0,-15,0)
            elif scaling == 'normalized':
                # scaling x by A and z by sqrt(A)
                scale = [ (-M)**(-5)*sqrt(A)**(12+i) for i in range(6) ]
        else:
            if A == 0:
                if B == 0:
                    return (1,0,0,1,0,0)
                elif scaling == 'normalized':
                    # scaling y by R/B**2
                    scale = [ (-M)**(-3)*(R/B**2)**i for i in range(6) ]
            elif scaling == 'normalized':
                # scaling y by R/A**4
                scale = [ (-M)**(-3)*(R/A**4)**i for i in range(6) ]
        D = -M
        Delta = A
        a = [0]
        a.append((2*K(3)**-1*A**2-B)*(N*A-M*B)*K(2)**-1 \
                    - M*(N*K(2)**-1-M*A*K(3)**-1))
        B0 = R
        B1 = K(2)**-1*(N*A-M*B)
        C0 = 0
        C1 = -M
    a[0] = (2*K(3)**-1*A**2-B)*R
    a.append(-D*B0 - K(2)**-1*Delta*a[0])
    a.append(-D*B1 - K(2)**-1*Delta*a[1])
    a.append(D**2*C0 + D*Delta*B0 + K(4)**-1*Delta**2*a[0])
    a.append(D**2*C1 + D*Delta*B1 + K(4)**-1*Delta**2*a[1])
    coeffs = tuple([K((-1)**i*binomial(5,i)*scale[5-i]*a[i]) for i in range(6)])
    if scaling == 'coprime':
        from sage.arith.misc import gcd
        return tuple([coeffs[i]/gcd(coeffs) for i in range(6)])
    else:
        return coeffs


######################################################################

def _reduce_invariants(invariants, weights):
    """
    Reduce a list of invariants of given weights.

    Reduces the given invariants over a field whose rings of integers is a gcd
    domain, in such a way that their greatest common weighted divisor is a unit
    in this ring.

    INPUT:

    - ``invariants`` -- The values of the invariants.

    - ``weights`` -- The respective weights of the invariants.

    OUTPUT:

    A list of invariants that is equivalent to the input and
    has no common weighted divisor.

    EXAMPLES::

        sage: from sage.rings.invariants.reconstruction import _reduce_invariants
        sage: invariants = [6/5, 12, 16]
        sage: weights = [1, 2, 3]
        sage: _reduce_invariants(invariants, weights)
        [3, 75, 250]
    """
    from sage.rings.integer_ring import ZZ
    factors = [dict(I.factor()) for I in invariants]
    scalar = ZZ(1)
    n = len(weights)
    from sage.arith.misc import gcd
    for prime in gcd(invariants).factor():
        p = prime[0]
        for D in factors:
            if p not in D:
                D[p] = 0
        scalar = scalar*p**min([factors[i][p]//weights[i] for i in range(n)])
    return [invariants[i]*scalar**-weights[i] for i in range(n)]
