# sage.doctest: needs sage.geometry.polyhedron sage.graphs
r"""
Map to the Weierstrass form of a toric elliptic curve

There are 16 reflexive polygons in 2-d. Each defines a toric Fano
variety, which (since it is 2-d) has a unique crepant resolution to a smooth
toric surface. An anticanonical hypersurface defines a genus one curve
`C` in this ambient space, with Jacobian elliptic curve `J(C)` which
can be defined by the Weierstrass model `y^2 = x^3 + f x + g`. The
coefficients `f` and `g` can be computed with the
:mod:`~sage.schemes.toric.weierstrass` module. The purpose of this
model is to give an explicit rational map `C \to J(C)`. This is an
`n^2`-cover, where `n` is the minimal multi-section of `C`.

Since it is technically often easier to deal with polynomials than
with fractions, we return the rational map in terms of homogeneous
coordinates. That is, the ambient space for the Weierstrass model is
the weighted projective space `\mathbb{P}^2[2,3,1]` with homogeneous
coordinates `[X:Y:Z] = [\lambda^2 X, \lambda^3 Y, \lambda Z]`. The
homogenized Weierstrass equation is

.. MATH::

    Y^2 = X^3 + f X Z^4 + g Z^6

EXAMPLES::

    sage: R.<x,y> = QQ[]
    sage: cubic = x^3 + y^3 + 1
    sage: f, g = WeierstrassForm(cubic);  (f,g)
    (0, -27/4)

That is, this hypersurface `C \in \mathbb{P}^2` has a Weierstrass
equation `Y^2 = X^3 + 0 \cdot X Z^4 - \frac{27}{4} Z^6` where
`[X:Y:Z]` are projective coordinates on `\mathbb{P}^2[2,3,1]`. The
form of the map `C\to J(C)` is::

    sage: X,Y,Z = WeierstrassForm(cubic, transformation=True);  (X,Y,Z)
    (-x^3*y^3 - x^3 - y^3,
     1/2*x^6*y^3 - 1/2*x^3*y^6 - 1/2*x^6 + 1/2*y^6 + 1/2*x^3 - 1/2*y^3,
     x*y)

Note that plugging in `[X:Y:Z]` to the Weierstrass equation is a
complicated polynomial, but contains the hypersurface equation as a
factor::

    sage: -Y^2 + X^3 + f*X*Z^4 + g*Z^6
    -1/4*x^12*y^6 - 1/2*x^9*y^9 - 1/4*x^6*y^12 + 1/2*x^12*y^3
    - 7/2*x^9*y^6 - 7/2*x^6*y^9 + 1/2*x^3*y^12 - 1/4*x^12 - 7/2*x^9*y^3
    - 45/4*x^6*y^6 - 7/2*x^3*y^9 - 1/4*y^12 - 1/2*x^9 - 7/2*x^6*y^3
    - 7/2*x^3*y^6 - 1/2*y^9 - 1/4*x^6 + 1/2*x^3*y^3 - 1/4*y^6
    sage: cubic.divides(-Y^2 + X^3 + f*X*Z^4 + g*Z^6)
    True

If you prefer you can also use homogeneous coordinates for `C \in
\mathbb{P}^2` ::

    sage: R.<x,y,z> = QQ[]
    sage: cubic = x^3 + y^3 + z^3
    sage: f, g = WeierstrassForm(cubic);  (f,g)
    (0, -27/4)
    sage: X,Y,Z = WeierstrassForm(cubic, transformation=True)
    sage: cubic.divides(-Y^2 + X^3 + f*X*Z^4 + g*Z^6)
    True

The 16 toric surfaces corresponding to the 16 reflexive polygons can
all be blown down to `\mathbb{P}^2`, `\mathbb{P}^1\times\mathbb{P}^1`,
or `\mathbb{P}^{2}[1,1,2]`. Their (and hence in all 16 cases)
anticanonical hypersurface can equally be brought into Weierstrass
form. For example, here is an anticanonical hypersurface in
`\mathbb{P}^{2}[1,1,2]` ::

    sage: P2_112 = toric_varieties.P2_112()
    sage: C = P2_112.anticanonical_hypersurface(coefficients=[1]*4);  C
    Closed subscheme of 2-d CPR-Fano toric variety
    covered by 3 affine patches defined by:
      z0^4 + z2^4 + z0*z1*z2 + z1^2
    sage: eq = C.defining_polynomials()[0]
    sage: f, g = WeierstrassForm(eq)
    sage: X,Y,Z = WeierstrassForm(eq, transformation=True)
    sage: (-Y^2 + X^3 + f*X*Z^4 + g*Z^6).reduce(C.defining_ideal())
    0

Finally, you sometimes have to manually specify the variables to
use. This is either because the equation is degenerate or because it
contains additional variables that you want to treat as coefficients::

    sage: R.<a, x,y,z> = QQ[]
    sage: cubic = x^3 + y^3 + z^3 + a*x*y*z
    sage: f, g = WeierstrassForm(cubic, variables=[x,y,z])
    sage: X,Y,Z = WeierstrassForm(cubic, variables=[x,y,z], transformation=True)
    sage: cubic.divides(-Y^2 + X^3 + f*X*Z^4 + g*Z^6)
    True

REFERENCES:

- [AKMMMP2002]_
"""

########################################################################
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
########################################################################

from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.rings.invariants.invariant_theory import invariant_theory
from sage.schemes.toric.weierstrass import (
    _partial_discriminant,
    _check_polynomial_P2,
    _check_polynomial_P1xP1,
    _check_polynomial_P2_112,
)


######################################################################

def WeierstrassMap(polynomial, variables=None):
    r"""
    Return the Weierstrass form of an anticanonical hypersurface.

    You should use
    :meth:`sage.schemes.toric.weierstrass.WeierstrassForm` with
    ``transformation=True`` to get the transformation. This function
    is only for internal use.

    INPUT:

    - ``polynomial`` -- a polynomial. The toric hypersurface
      equation. Can be either a cubic, a biquadric, or the
      hypersurface in `\mathbb{P}^2[1,1,2]`. The equation need not be
      in any standard form, only its Newton polyhedron is used.

    - ``variables`` -- list of variables of the parent polynomial
      ring or ``None`` (default). In the latter case, all variables
      are taken to be polynomial ring variables. If a subset of
      polynomial ring variables are given, the Weierstrass form is
      determined over the function field generated by the remaining
      variables.

    OUTPUT:

    A triple `(X,Y,Z)` of polynomials defining a rational map of the
    toric hypersurface to its Weierstrass form in
    `\mathbb{P}^2[2,3,1]`. That is, the triple satisfies

    .. MATH::

        Y^2 = X^3 + f X Z^4 + g Z^6

    when restricted to the toric hypersurface.

    EXAMPLES::

        sage: R.<x,y,z> = QQ[]
        sage: cubic = x^3 + y^3 + z^3
        sage: X,Y,Z = WeierstrassForm(cubic, transformation=True);  (X,Y,Z)
        (-x^3*y^3 - x^3*z^3 - y^3*z^3,
         1/2*x^6*y^3 - 1/2*x^3*y^6 - 1/2*x^6*z^3 + 1/2*y^6*z^3
             + 1/2*x^3*z^6 - 1/2*y^3*z^6,
         x*y*z)
        sage: f, g = WeierstrassForm(cubic);  (f,g)
        (0, -27/4)
        sage: cubic.divides(-Y^2 + X^3 + f*X*Z^4 + g*Z^6)
        True

    Only the affine span of the Newton polytope of the polynomial
    matters. For example::

        sage: WeierstrassForm(cubic.subs(z=1), transformation=True)
        (-x^3*y^3 - x^3 - y^3,
         1/2*x^6*y^3 - 1/2*x^3*y^6 - 1/2*x^6
             + 1/2*y^6 + 1/2*x^3 - 1/2*y^3,
         x*y)
        sage: WeierstrassForm(x * cubic, transformation=True)
        (-x^3*y^3 - x^3*z^3 - y^3*z^3,
         1/2*x^6*y^3 - 1/2*x^3*y^6 - 1/2*x^6*z^3 + 1/2*y^6*z^3
             + 1/2*x^3*z^6 - 1/2*y^3*z^6,
         x*y*z)

    This allows you to work with either homogeneous or inhomogeneous
    variables. For example, here is the del Pezzo surface of degree 8::

        sage: dP8 = toric_varieties.dP8()
        sage: dP8.inject_variables()
        Defining t, x, y, z
        sage: WeierstrassForm(x*y^2 + y^2*z + t^2*x^3 + t^2*z^3, transformation=True)
        (-1/27*t^4*x^6 - 2/27*t^4*x^5*z - 5/27*t^4*x^4*z^2
             - 8/27*t^4*x^3*z^3 - 5/27*t^4*x^2*z^4 - 2/27*t^4*x*z^5
             - 1/27*t^4*z^6 - 4/81*t^2*x^4*y^2 - 4/81*t^2*x^3*y^2*z
             - 4/81*t^2*x*y^2*z^3 - 4/81*t^2*y^2*z^4 - 2/81*x^2*y^4
             - 4/81*x*y^4*z - 2/81*y^4*z^2,
         0,
         1/3*t^2*x^2*z + 1/3*t^2*x*z^2 - 1/9*x*y^2 - 1/9*y^2*z)
        sage: WeierstrassForm(x*y^2 + y^2 + x^3 + 1, transformation=True)
        (-1/27*x^6 - 4/81*x^4*y^2 - 2/81*x^2*y^4 - 2/27*x^5
             - 4/81*x^3*y^2 - 4/81*x*y^4 - 5/27*x^4 - 2/81*y^4 - 8/27*x^3
             - 4/81*x*y^2 - 5/27*x^2 - 4/81*y^2 - 2/27*x - 1/27,
         0,
         -1/9*x*y^2 + 1/3*x^2 - 1/9*y^2 + 1/3*x)

    By specifying only certain variables we can compute the
    Weierstrass form over the function field generated by the
    remaining variables. For example, here is a cubic over `\QQ[a]` ::

        sage: R.<a, x,y,z> = QQ[]
        sage: cubic = x^3 + a*y^3 + a^2*z^3
        sage: WeierstrassForm(cubic, variables=[x,y,z], transformation=True)
        (-a^9*y^3*z^3 - a^8*x^3*z^3 - a^7*x^3*y^3,
         -1/2*a^14*y^3*z^6 + 1/2*a^13*y^6*z^3 + 1/2*a^13*x^3*z^6
             - 1/2*a^11*x^3*y^6 - 1/2*a^11*x^6*z^3 + 1/2*a^10*x^6*y^3,
         a^3*x*y*z)

    TESTS::

        sage: for P in ReflexivePolytopes(2):
        ....:     S = ToricVariety(FaceFan(P))
        ....:     p = sum( (-S.K()).sections_monomials() )
        ....:     f, g = WeierstrassForm(p)
        ....:     X,Y,Z = WeierstrassForm(p, transformation=True)
        ....:     assert p.divides(-Y^2 + X^3 + f*X*Z^4 + g*Z^6)
    """
    if variables is None:
        variables = polynomial.variables()
    # switch to suitable inhomogeneous coordinates
    from sage.geometry.polyhedron.ppl_lattice_polygon import (
        polar_P2_polytope, polar_P1xP1_polytope, polar_P2_112_polytope)
    from sage.schemes.toric.weierstrass import Newton_polygon_embedded
    newton_polytope, polynomial_aff, variables_aff = \
        Newton_polygon_embedded(polynomial, variables)
    polygon = newton_polytope.embed_in_reflexive_polytope('polytope')
    # Compute the map in inhomogeneous coordinates
    if polygon is polar_P2_polytope():
        X, Y, Z = WeierstrassMap_P2(polynomial_aff, variables_aff)
    elif polygon is polar_P1xP1_polytope():
        X, Y, Z = WeierstrassMap_P1xP1(polynomial_aff, variables_aff)
    elif polygon is polar_P2_112_polytope():
        X, Y, Z = WeierstrassMap_P2_112(polynomial_aff, variables_aff)
    else:
        assert False, 'Newton polytope is not contained in a reflexive polygon'
    # homogenize again
    R = polynomial.parent()
    x = R.gens().index(variables_aff[0])
    y = R.gens().index(variables_aff[1])
    hom = newton_polytope.embed_in_reflexive_polytope('hom')

    def homogenize(inhomog, degree):
        e = tuple(hom._A * vector(ZZ, [inhomog[x], inhomog[y]]) + degree * hom._b)
        result = list(inhomog)
        for i, var in enumerate(variables):
            result[R.gens().index(var)] = e[i]
        result = vector(ZZ, result)
        result.set_immutable()
        return result
    X_dict = {homogenize(e, 2): v for e, v in X.monomial_coefficients().items()}
    Y_dict = {homogenize(e, 3): v for e, v in Y.monomial_coefficients().items()}
    Z_dict = {homogenize(e, 1): v for e, v in Z.monomial_coefficients().items()}
    # shift to nonnegative exponents if necessary
    min_deg = [0] * R.ngens()
    for var in variables:
        i = R.gens().index(var)
        min_X = min([e[i] for e in X_dict]) if X_dict else 0
        min_Y = min([e[i] for e in Y_dict]) if Y_dict else 0
        min_Z = min([e[i] for e in Z_dict]) if Z_dict else 0
        min_deg[i] = min(min_X / 2, min_Y / 3, min_Z)
    min_deg = vector(min_deg)
    X_dict = {tuple(e - 2 * min_deg): v for e, v in X_dict.items()}
    Y_dict = {tuple(e - 3 * min_deg): v for e, v in Y_dict.items()}
    Z_dict = {tuple(e - 1 * min_deg): v for e, v in Z_dict.items()}
    return (R(X_dict), R(Y_dict), R(Z_dict))


######################################################################
#
#  Weierstrass form of cubic in P^2
#
######################################################################

def WeierstrassMap_P2(polynomial, variables=None):
    r"""
    Map a cubic to its Weierstrass form.

    Input/output is the same as :func:`WeierstrassMap`, except that
    the input polynomial must be a cubic in `\mathbb{P}^2`,

    .. MATH::

        \begin{split}
          p(x,y) =&\;
          a_{30} x^{3} + a_{21} x^{2} y + a_{12} x y^{2} +
          a_{03} y^{3} + a_{20} x^{2} +
          \\ &\;
          a_{11} x y +
          a_{02} y^{2} + a_{10} x + a_{01} y + a_{00}
        \end{split}

    EXAMPLES::

        sage: from sage.schemes.toric.weierstrass import WeierstrassForm_P2
        sage: from sage.schemes.toric.weierstrass_covering import WeierstrassMap_P2
        sage: R.<x,y,z> = QQ[]
        sage: equation =  x^3 + y^3 + z^3 + x*y*z
        sage: f, g = WeierstrassForm_P2(equation)
        sage: X,Y,Z = WeierstrassMap_P2(equation)
        sage: equation.divides(-Y^2 + X^3 + f*X*Z^4 + g*Z^6)
        True

        sage: from sage.schemes.toric.weierstrass import WeierstrassForm_P2
        sage: from sage.schemes.toric.weierstrass_covering import WeierstrassMap_P2
        sage: R.<x,y> = QQ[]
        sage: equation =  x^3 + y^3 + 1
        sage: f, g = WeierstrassForm_P2(equation)
        sage: X,Y,Z = WeierstrassMap_P2(equation)
        sage: equation.divides(-Y^2 + X^3 + f*X*Z^4 + g*Z^6)
        True
    """
    x, y, z = _check_polynomial_P2(polynomial, variables)
    cubic = invariant_theory.ternary_cubic(polynomial, x, y, z)
    H = cubic.Hessian()
    Theta = cubic.Theta_covariant()
    J = cubic.J_covariant()
    F = polynomial.parent().base_ring()
    return (Theta, J / F(2), H)


######################################################################
#
#  Weierstrass form of biquadric in P1 x P1
#
######################################################################

def WeierstrassMap_P1xP1(polynomial, variables=None):
    r"""
    Map an anticanonical hypersurface in
    `\mathbb{P}^1 \times \mathbb{P}^1` into Weierstrass form.

    Input/output is the same as :func:`WeierstrassMap`, except that
    the input polynomial must be a standard anticanonical hypersurface
    in the toric surface `\mathbb{P}^1 \times \mathbb{P}^1`:

    EXAMPLES::

        sage: from sage.schemes.toric.weierstrass_covering import WeierstrassMap_P1xP1
        sage: from sage.schemes.toric.weierstrass import WeierstrassForm_P1xP1
        sage: R.<x0,x1,y0,y1,a> = QQ[]
        sage: biquadric = (x0^2*y0^2 + x1^2*y0^2 + x0^2*y1^2 + x1^2*y1^2 +
        ....:     a * x0*x1*y0*y1*5)
        sage: f, g = WeierstrassForm_P1xP1(biquadric, [x0, x1, y0, y1]);  (f,g)
        (-625/48*a^4 + 25/3*a^2 - 16/3, 15625/864*a^6 - 625/36*a^4 - 100/9*a^2 + 128/27)
        sage: X, Y, Z = WeierstrassMap_P1xP1(biquadric, [x0, x1, y0, y1])
        sage: (-Y^2 + X^3 + f*X*Z^4 + g*Z^6).reduce(R.ideal(biquadric))                 # needs sage.libs.singular
        0

        sage: R = PolynomialRing(QQ, 'x,y,s,t', order='lex')
        sage: R.inject_variables()
        Defining x, y, s, t
        sage: equation = (s^2*(x^2+2*x*y+3*y^2) + s*t*(4*x^2+5*x*y+6*y^2)
        ....:             + t^2*(7*x^2+8*x*y+9*y^2))
        sage: X, Y, Z = WeierstrassMap_P1xP1(equation, [x,y,s,t])
        sage: f, g = WeierstrassForm_P1xP1(equation, variables=[x,y,s,t])
        sage: (-Y^2 + X^3 + f*X*Z^4 + g*Z^6).reduce(R.ideal(equation))                  # needs sage.libs.singular
        0

        sage: R = PolynomialRing(QQ, 'x,s', order='lex')
        sage: R.inject_variables()
        Defining x, s
        sage: equation = s^2*(x^2+2*x+3) + s*(4*x^2+5*x+6) + (7*x^2+8*x+9)
        sage: X, Y, Z = WeierstrassMap_P1xP1(equation)
        sage: f, g = WeierstrassForm_P1xP1(equation)
        sage: (-Y^2 + X^3 + f*X*Z^4 + g*Z^6).reduce(R.ideal(equation))                  # needs sage.libs.singular
        0
    """
    x, y, s, t = _check_polynomial_P1xP1(polynomial, variables)
    a00 = polynomial.coefficient({s: 2})
    V = polynomial.coefficient({s: 1})
    U = - _partial_discriminant(polynomial, s, t) / 4
    Q = invariant_theory.binary_quartic(U, x, y)
    g = Q.g_covariant()
    h = Q.h_covariant()
    if t is None:
        t = 1
    return (4 * g * t**2, 4 * h * t**3, (a00 * s + V / 2))


######################################################################
#
#  Weierstrass form of anticanonical hypersurface in WP2[1,1,2]
#
######################################################################

def WeierstrassMap_P2_112(polynomial, variables=None):
    r"""
    Map an anticanonical hypersurface in `\mathbb{P}^2[1,1,2]` into Weierstrass form.

    Input/output is the same as :func:`WeierstrassMap`, except that
    the input polynomial must be a standard anticanonical hypersurface
    in weighted projective space `\mathbb{P}^2[1,1,2]`:

    .. MATH::

        \begin{split}
          p(x,y) =&\;
          a_{40} x^4 +
          a_{30} x^3 +
          a_{21} x^2 y +
          a_{20} x^2 +
          \\ &\;
          a_{11} x y +
          a_{02} y^2 +
          a_{10} x +
          a_{01} y +
          a_{00}
        \end{split}

    EXAMPLES::

        sage: from sage.schemes.toric.weierstrass_covering import WeierstrassMap_P2_112
        sage: from sage.schemes.toric.weierstrass import WeierstrassForm_P2_112
        sage: R = PolynomialRing(QQ, 'x,y,a0,a1,a2,a3,a4', order='lex')
        sage: R.inject_variables()
        Defining x, y, a0, a1, a2, a3, a4
        sage: equation = y^2 + a0*x^4 + 4*a1*x^3 + 6*a2*x^2 + 4*a3*x + a4
        sage: X, Y, Z = WeierstrassMap_P2_112(equation, [x,y])
        sage: f, g = WeierstrassForm_P2_112(equation, variables=[x,y])
        sage: (-Y^2 + X^3 + f*X*Z^4 + g*Z^6).reduce(R.ideal(equation))                  # needs sage.libs.singular
        0

    Another example, this time in homogeneous coordinates::

        sage: fan = Fan(rays=[(1,0),(0,1),(-1,-2),(0,-1)],
        ....:           cones=[[0,1],[1,2],[2,3],[3,0]])
        sage: P112.<x,y,z,t> = ToricVariety(fan)
        sage: (-P112.K()).sections_monomials()
        (z^4*t^2, x*z^3*t^2, x^2*z^2*t^2, x^3*z*t^2,
         x^4*t^2, y*z^2*t, x*y*z*t, x^2*y*t, y^2)
        sage: C_eqn = sum(_)
        sage: C = P112.subscheme(C_eqn)
        sage: WeierstrassForm_P2_112(C_eqn, [x,y,z,t])
        (-97/48, 17/864)
        sage: X, Y, Z = WeierstrassMap_P2_112(C_eqn, [x,y,z,t])
        sage: (-Y^2 + X^3 - 97/48*X*Z^4 + 17/864*Z^6).reduce(C.defining_ideal())        # needs sage.libs.singular
        0
    """
    x, y, z, t = _check_polynomial_P2_112(polynomial, variables)
    a00 = polynomial.coefficient({y: 2})
    V = polynomial.coefficient({y: 1})
    U = - _partial_discriminant(polynomial, y, t) / 4
    Q = invariant_theory.binary_quartic(U, x, z)
    g = Q.g_covariant()
    h = Q.h_covariant()
    if t is None:
        t = 1
    return (4 * g * t**2, 4 * h * t**3, (a00 * y + V / 2))
