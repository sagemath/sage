"""
Some Extras
"""
##
## Some extra routines to make the QuadraticForm class more useful.
##

from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.quadratic_forms.quadratic_form import QuadraticForm



def BezoutianQuadraticForm(f, g):
    r"""
    Compute the Bezoutian of two polynomials defined over a common base ring.  This is defined by

    .. MATH::

        {\rm Bez}(f, g) := \frac{f(x) g(y) - f(y) g(x)}{y - x}

    and has size defined by the maximum of the degrees of `f` and `g`.

    INPUT:

    - `f`, `g` -- polynomials in `R[x]`, for some ring `R`

    OUTPUT:

    a quadratic form over `R`

    EXAMPLES::

        sage: R = PolynomialRing(ZZ, 'x')
        sage: f = R([1,2,3])
        sage: g = R([2,5])
        sage: Q = BezoutianQuadraticForm(f, g) ; Q
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 1 -12 ]
        [ * -15 ]

    AUTHORS:

    - Fernando Rodriguez-Villegas, Jonathan Hanke -- added on 11/9/2008

    """
    ## Check that f and g are polynomials with a common base ring
    if not is_Polynomial(f) or not is_Polynomial(g):
        raise TypeError("Oops!  One of your inputs is not a polynomial. =(")
    if f.base_ring() != g.base_ring():                   ## TO DO:  Change this to allow coercion!
        raise TypeError("Oops!  These polynomials are not defined over the same coefficient ring.")

    ## Initialize the quadratic form
    R = f.base_ring()
    P = PolynomialRing(R, ['x','y'])
    a, b = P.gens()
    n = max(f.degree(), g.degree())
    Q = QuadraticForm(R, n)

    ## Set the coefficients of Bezoutian
    bez_poly = (f(a) * g(b) - f(b) * g(a)) // (b - a)    ## Truncated (exact) division here
    for i in range(n):
        for j in range(i, n):
            if i == j:
                Q[i,j] = bez_poly.coefficient({a:i,b:j})
            else:
                Q[i,j] = bez_poly.coefficient({a:i,b:j}) * 2

    return Q


def HyperbolicPlane_quadratic_form(R, r=1):
    """
    Constructs the direct sum of `r` copies of the quadratic form `xy`
    representing a hyperbolic plane defined over the base ring `R`.

    INPUT:

    - `R`: a ring
    - `n` (integer, default 1) number of copies

    EXAMPLES::

        sage: HyperbolicPlane_quadratic_form(ZZ)
        Quadratic form in 2 variables over Integer Ring with coefficients:
        [ 0 1 ]
        [ * 0 ]

    """
    r = ZZ(r)
    ## Check that the multiplicity is a natural number
    if r < 1:
        raise TypeError("The multiplicity r must be a natural number.")

    H = QuadraticForm(R, 2, [0, 1, 0])
    return sum([H  for i in range(r-1)], H)
