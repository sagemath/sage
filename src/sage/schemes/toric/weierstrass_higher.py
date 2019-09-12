r"""
Weierstrass for elliptic curves in higher codimension

The :mod:`~sage.schemes.toric.weierstrass` module lets you transform a
genus-one curve, given as a hypersurface in a toric surface, into
Weierstrass form. The purpose of this module is to extend this to
higher codimension subschemes of toric varieties. In general, this is
an unsolved problem. However, for certain special cases this is known.

The simplest codimension-two case is the complete intersection of two
quadratic equations in `\mathbb{P}^3` ::

    sage: R.<w,x,y,z> = QQ[]
    sage: quadratic1 = w^2+x^2+y^2
    sage: quadratic2 = z^2 + w*x
    sage: WeierstrassForm([quadratic1, quadratic2])
    (-1/4, 0)

Hence, the Weierstrass form of this complete intersection is $Y^2 =
X^3 - \frac{1}{4} X Z^4$.
"""

#*****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import PolynomialRing
from sage.rings.invariants.invariant_theory import invariant_theory
from sage.schemes.toric.weierstrass import _check_homogeneity


######################################################################
def WeierstrassForm2(polynomial, variables=None, transformation=False):
    r"""
    Helper function for :func:`~sage.schemes.toric.weierstrass.WeierstrassForm`

    Currently, only the case of the complete intersection of two
    quadratic equations in `\mathbb{P}^3` is supported.

    INPUT / OUTPUT:

    See :func:`~sage.schemes.toric.weierstrass.WeierstrassForm`

    TESTS::

        sage: from sage.schemes.toric.weierstrass_higher import WeierstrassForm2
        sage: R.<w,x,y,z> = QQ[]
        sage: quadratic1 = w^2+x^2+y^2
        sage: quadratic2 = z^2 + w*x
        sage: WeierstrassForm2([quadratic1, quadratic2])
        (-1/4, 0)
    """
    if transformation:
        return WeierstrassMap_P3(*polynomial, variables=variables)
    else:
        return WeierstrassForm_P3(*polynomial, variables=variables)


######################################################################
#
#  Weierstrass form of complete intersection of two quadratics  in P^3
#
######################################################################
def _check_polynomials_P3(quadratic1, quadratic2, variables):
    """
    Check that the polynomial is weighted homogeneous in standard variables.

    INPUT:

    - ``quadratic1``, ``quadratic2`` -- two quadratic polynomials in 4
      homogeneous or 3 inhomogeneous variables.

    - ``variables`` -- the variables or ``None`` (default).

    OUTPUT:

    This function returns ``variables``, potentially guessed from the
    polynomial ring. A ``ValueError`` is raised if the polynomial is
    not homogeneous.

    EXAMPLES::

        sage: from sage.schemes.toric.weierstrass_higher import _check_polynomials_P3
        sage: R.<w,x,y,z> = QQ[]
        sage: quadratic = w^2+x^2+y^2+z^2
        sage: _check_polynomials_P3(w^2, quadratic, [w,x,y,z])
        (w, x, y, z)
        sage: _check_polynomials_P3(w^2, quadratic, None)
        (w, x, y, z)
        sage: _check_polynomials_P3(z^2, quadratic.subs(w=0), None)
        (x, y, z, None)
        sage: R.<w,x,y,z,t> = QQ[]
        sage: quadratic = w^2+x^2+y^2+z^2 + t*(x*y+y*z+z*w+w*x)
        sage: _check_polynomials_P3(w^2, quadratic, [w,x,y,z])
        (w, x, y, z)
        sage: _check_polynomials_P3(w^2, quadratic, [w,x,y,t])
        Traceback (most recent call last):
        ...
        ValueError: The polynomial is not homogeneous with weights (1, 1, 1, 1)
    """
    if quadratic1.parent() is not quadratic2.parent():
        raise ValueError('The two quadratics must be in the same polynomial ring.')
    if variables is None:
        variables = quadratic1.variables() + quadratic2.variables()
        variables = sorted(set(variables), reverse=True)
    if len(variables) == 4:
        w, x, y, z = variables
        _check_homogeneity(quadratic1, [w, x, y, z], (1, 1, 1, 1), 2)
        _check_homogeneity(quadratic2, [w, x, y, z], (1, 1, 1, 1), 2)
    elif len(variables) == 3:
        w, x, y = variables
        z = None
    else:
        raise ValueError('Need three or four variables, got '+str(variables))
    return (w, x, y, z)



######################################################################
def _biquadratic_syzygy_quartic(quadratic1, quadratic2, variables=None):
    r"""
    Helper function for the Weierstrass form of a biquadratic in $`\mathbb{P}^3$

    The invariants and covariants of a quaternary biquadratic satisfy
    the relation
    :meth:`sage.rings.invariant_theory.TwoQuaternaryQuadratics.syzygy`,
    which is (modulo the two quadratic equations) of the form $J^2 =
    p_4(T, T')$ where

    * $J$, $T$, $T'$ are the covariants of the biquadratic.

    * $p_4$ is some quartic polynomial whose coefficients are
      invariants of the biquadratic.

    INPUT:

    See :func:`WeierstrassForm_P3`

    OUTPUT:

    A triple consisting of

    - The quaternary biquadratic as an algebraic form
      :class:`~sage.rings.invariant_theory.TwoQuaternaryQuadratics`

    - The binary quartic $p_4$ as a
      :class:`~sage.rings.invariant_theory.BinaryQuartic`

    - The dictionary of variable substitutions from the variables of
      the quartic to the variables of the biquadratic.

    EXAMPLES::

        sage: from sage.schemes.toric.weierstrass_higher import _biquadratic_syzygy_quartic
        sage: R.<w,x,y,z> = QQ[]
        sage: _biquadratic_syzygy_quartic(w^2+x^2+y^2, z^2)
        (Joint quaternary quadratic with coefficients (1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
         and quaternary quadratic with coefficients (0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
         Binary quartic with coefficients (0, 0, 0, -1, 0), {aux...})
    """
    w, x, y, z = _check_polynomials_P3(quadratic1, quadratic2, variables)
    biquadratic = invariant_theory.quaternary_biquadratic(quadratic1, quadratic2, [w, x, y, z])

    # construct auxiliary polynomial ring to work with the rhs of the syzygy
    R = biquadratic.ring()
    n = R.ngens()
    R_aux = PolynomialRing(R.base_ring(), n+2, 'aux')
    to_aux = dict()
    from_aux = dict()
    for var, var_aux in zip(R.gens(), R_aux.gens()[0:n]):
        to_aux[var] = var_aux
        from_aux[var_aux] = var
    T, T_prime = R_aux.gens()[n:]
    from_aux[T] = biquadratic.T_covariant()
    from_aux[T_prime] = biquadratic.T_prime_covariant()

    # Syzygy is J^2 = syz_rhs + (terms that vanish on the biquadratic) with
    # J = biquadratic.J_covariant()
    syz_rhs = T**4 * biquadratic.Delta_invariant().subs(to_aux) \
        - T**3*T_prime * biquadratic.Theta_invariant().subs(to_aux) \
        + T**2*T_prime**2 * biquadratic.Phi_invariant().subs(to_aux) \
        - T*T_prime**3 * biquadratic.Theta_prime_invariant().subs(to_aux) \
        + T_prime**4 * biquadratic.Delta_prime_invariant().subs(to_aux)
    quartic = invariant_theory.binary_quartic(syz_rhs, [T, T_prime])
    return (biquadratic, quartic, from_aux)


######################################################################
def WeierstrassForm_P3(quadratic1, quadratic2, variables=None):
    r"""
    Bring a complete intersection of two quadratics into Weierstrass form.

    Input/output is the same as
    :func:`sage.schemes.toric.weierstrass.WeierstrassForm`, except
    that the two input polynomials must be quadratic polynomials in
    `\mathbb{P}^3`.

    EXAMPLES::

        sage: from sage.schemes.toric.weierstrass_higher import WeierstrassForm_P3
        sage: R.<w,x,y,z> = QQ[]
        sage: quadratic1 = w^2+x^2+y^2
        sage: quadratic2 = z^2 + w*x
        sage: WeierstrassForm_P3(quadratic1, quadratic2)
        (-1/4, 0)

    TESTS::

        sage: R.<w,x,y,z,a0,a1,a2,a3,b0,b1,b2,b3,b4,b5> = QQ[]
        sage: p1 = w^2 + x^2 + y^2 + z^2
        sage: p2 = a0*w^2 + a1*x^2 + a2*y^2 + a3*z^2
        sage: p2 += b0*x*y + b1*x*z + b2*x*w + b3*y*z + b4*y*w + b5*z*w
        sage: a, b = WeierstrassForm_P3(p1, p2, [w,x,y,z])
        sage: a.total_degree(), len(a.coefficients())
        (4, 107)
        sage: b.total_degree(), len(b.coefficients())
        (6, 648)
    """
    biquadratic, quartic, from_aux = \
        _biquadratic_syzygy_quartic(quadratic1, quadratic2, variables=variables)
    a = quartic.EisensteinD().subs(from_aux)
    b = quartic.EisensteinE().subs(from_aux)
    return (-4*a, 16*b)


######################################################################
def WeierstrassMap_P3(quadratic1, quadratic2, variables=None):
    r"""
    Bring a complete intersection of two quadratics into Weierstrass form.

    Input/output is the same as
    :func:`sage.schemes.toric.weierstrass.WeierstrassForm`, except
    that the two input polynomials must be quadratic polynomials in
    `\mathbb{P}^3`.

    EXAMPLES::

        sage: from sage.schemes.toric.weierstrass_higher import \
        ....:     WeierstrassMap_P3, WeierstrassForm_P3
        sage: R.<w,x,y,z> = QQ[]
        sage: quadratic1 = w^2+x^2+y^2
        sage: quadratic2 = z^2 + w*x
        sage: X, Y, Z = WeierstrassMap_P3(quadratic1, quadratic2)
        sage: X
        1/1024*w^8 + 3/256*w^6*x^2 + 19/512*w^4*x^4 + 3/256*w^2*x^6 + 1/1024*x^8
        sage: Y
        1/32768*w^12 - 7/16384*w^10*x^2 - 145/32768*w^8*x^4 - 49/8192*w^6*x^6
        - 145/32768*w^4*x^8 - 7/16384*w^2*x^10 + 1/32768*x^12
        sage: Z
        -1/8*w^2*y*z + 1/8*x^2*y*z

        sage: a, b = WeierstrassForm_P3(quadratic1, quadratic2);  a, b
        (-1/4, 0)

        sage: ideal = R.ideal(quadratic1, quadratic2)
        sage: (-Y^2 + X^3 + a*X*Z^4 + b*Z^6).reduce(ideal)
        0

    TESTS::

        sage: R.<w,x,y,z,a0,a1,a2,a3> = GF(101)[]
        sage: p1 = w^2 + x^2 + y^2 + z^2
        sage: p2 = a0*w^2 + a1*x^2 + a2*y^2 + a3*z^2
        sage: X, Y, Z = WeierstrassMap_P3(p1, p2, [w,x,y,z])
        sage: X.total_degree(), len(X.coefficients())
        (22, 4164)
        sage: Y.total_degree(), len(Y.coefficients())
        (33, 26912)
        sage: Z.total_degree(), len(Z.coefficients())
        (10, 24)
        sage: Z
        w*x*y*z*a0^3*a1^2*a2 - w*x*y*z*a0^2*a1^3*a2 - w*x*y*z*a0^3*a1*a2^2
        + w*x*y*z*a0*a1^3*a2^2 + w*x*y*z*a0^2*a1*a2^3 - w*x*y*z*a0*a1^2*a2^3
        - w*x*y*z*a0^3*a1^2*a3 + w*x*y*z*a0^2*a1^3*a3 + w*x*y*z*a0^3*a2^2*a3
        - w*x*y*z*a1^3*a2^2*a3 - w*x*y*z*a0^2*a2^3*a3 + w*x*y*z*a1^2*a2^3*a3
        + w*x*y*z*a0^3*a1*a3^2 - w*x*y*z*a0*a1^3*a3^2 - w*x*y*z*a0^3*a2*a3^2
        + w*x*y*z*a1^3*a2*a3^2 + w*x*y*z*a0*a2^3*a3^2 - w*x*y*z*a1*a2^3*a3^2
        - w*x*y*z*a0^2*a1*a3^3 + w*x*y*z*a0*a1^2*a3^3 + w*x*y*z*a0^2*a2*a3^3
        - w*x*y*z*a1^2*a2*a3^3 - w*x*y*z*a0*a2^2*a3^3 + w*x*y*z*a1*a2^2*a3^3
    """
    biquadratic, quartic, from_aux = \
        _biquadratic_syzygy_quartic(quadratic1, quadratic2, variables=variables)
    J = biquadratic.J_covariant()
    g = quartic.g_covariant().subs(from_aux)
    h = quartic.h_covariant().subs(from_aux)
    return (4*g, 4*h, J)

