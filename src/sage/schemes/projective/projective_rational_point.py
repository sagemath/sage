r"""
Enumeration of rational points on projective schemes

Naive algorithms for enumerating rational points over `\QQ` or finite fields
over for general schemes.

.. WARNING::

    Incorrect results and infinite loops may occur if using a wrong function.
    (For instance using an affine function for a projective scheme or a finite
    field function for a scheme defined over an infinite field.)

EXAMPLES:

Projective, over `\QQ`::

    sage: from sage.schemes.projective.projective_rational_point import enum_projective_rational_field
    sage: P.<X,Y,Z> = ProjectiveSpace(2,QQ)
    sage: C = P.subscheme([X+Y-Z])
    sage: enum_projective_rational_field(C,3)
    [(-2 : 3 : 1), (-1 : 1 : 0), (-1 : 2 : 1), (-1/2 : 3/2 : 1),
     (0 : 1 : 1), (1/3 : 2/3 : 1), (1/2 : 1/2 : 1), (2/3 : 1/3 : 1),
     (1 : 0 : 1), (3/2 : -1/2 : 1), (2 : -1 : 1), (3 : -2 : 1)]


Projective over a finite field::

    sage: from sage.schemes.projective.projective_rational_point import enum_projective_finite_field
    sage: E = EllipticCurve('72').change_ring(GF(19))
    sage: enum_projective_finite_field(E)
    [(0 : 1 : 0), (1 : 0 : 1), (3 : 0 : 1), (4 : 9 : 1), (4 : 10 : 1),
     (6 : 6 : 1), (6 : 13 : 1), (7 : 6 : 1), (7 : 13 : 1), (9 : 4 : 1),
     (9 : 15 : 1), (12 : 8 : 1), (12 : 11 : 1), (13 : 8 : 1), (13 : 11 : 1),
     (14 : 3 : 1), (14 : 16 : 1), (15 : 0 : 1), (16 : 9 : 1), (16 : 10 : 1),
     (17 : 7 : 1), (17 : 12 : 1), (18 : 9 : 1), (18 : 10 : 1)]


AUTHORS:

- David R. Kohel <kohel@maths.usyd.edu.au>: original version.

- John Cremona and Charlie Turner <charlotteturner@gmail.com> (06-2010):
  improvements to clarity and documentation.
"""


#*******************************************************************************
#  Copyright (C) 2010 William Stein, David Kohel, John Cremona, Charlie Turner
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************


from sage.rings.arith import gcd
from sage.rings.all import ZZ
from sage.misc.all import srange, cartesian_product_iterator
from sage.schemes.generic.scheme import is_Scheme

def enum_projective_rational_field(X,B):
    r"""
    Enumerates projective, rational points on scheme ``X`` of height up to
    bound ``B``.

    INPUT:

    - ``X`` -  a scheme or set of abstract rational points of a scheme;
    - ``B`` -  a positive integer bound.

    OUTPUT:

    - a list containing the projective points of ``X`` of height up to ``B``,
      sorted.

    EXAMPLES::

        sage: P.<X,Y,Z> = ProjectiveSpace(2,QQ)
        sage: C = P.subscheme([X+Y-Z])
        sage: from sage.schemes.projective.projective_rational_point import enum_projective_rational_field
        sage: enum_projective_rational_field(C(QQ),6)
        [(-5 : 6 : 1), (-4 : 5 : 1), (-3 : 4 : 1), (-2 : 3 : 1),
         (-3/2 : 5/2 : 1), (-1 : 1 : 0), (-1 : 2 : 1), (-2/3 : 5/3 : 1),
         (-1/2 : 3/2 : 1), (-1/3 : 4/3 : 1), (-1/4 : 5/4 : 1),
         (-1/5 : 6/5 : 1), (0 : 1 : 1), (1/6 : 5/6 : 1), (1/5 : 4/5 : 1),
         (1/4 : 3/4 : 1), (1/3 : 2/3 : 1), (2/5 : 3/5 : 1), (1/2 : 1/2 : 1),
         (3/5 : 2/5 : 1), (2/3 : 1/3 : 1), (3/4 : 1/4 : 1), (4/5 : 1/5 : 1),
         (5/6 : 1/6 : 1), (1 : 0 : 1), (6/5 : -1/5 : 1), (5/4 : -1/4 : 1),
         (4/3 : -1/3 : 1), (3/2 : -1/2 : 1), (5/3 : -2/3 : 1), (2 : -1 : 1),
         (5/2 : -3/2 : 1), (3 : -2 : 1), (4 : -3 : 1), (5 : -4 : 1),
         (6 : -5 : 1)]
        sage: enum_projective_rational_field(C,6) == enum_projective_rational_field(C(QQ),6)
        True

    ::

        sage: P3.<W,X,Y,Z> = ProjectiveSpace(3,QQ)
        sage: enum_projective_rational_field(P3,1)
        [(-1 : -1 : -1 : 1), (-1 : -1 : 0 : 1), (-1 : -1 : 1 : 0), (-1 : -1 : 1 : 1),
        (-1 : 0 : -1 : 1), (-1 : 0 : 0 : 1), (-1 : 0 : 1 : 0), (-1 : 0 : 1 : 1),
        (-1 : 1 : -1 : 1), (-1 : 1 : 0 : 0), (-1 : 1 : 0 : 1), (-1 : 1 : 1 : 0),
        (-1 : 1 : 1 : 1), (0 : -1 : -1 : 1), (0 : -1 : 0 : 1), (0 : -1 : 1 : 0),
        (0 : -1 : 1 : 1), (0 : 0 : -1 : 1), (0 : 0 : 0 : 1), (0 : 0 : 1 : 0),
        (0 : 0 : 1 : 1), (0 : 1 : -1 : 1), (0 : 1 : 0 : 0), (0 : 1 : 0 : 1),
        (0 : 1 : 1 : 0), (0 : 1 : 1 : 1), (1 : -1 : -1 : 1), (1 : -1 : 0 : 1),
        (1 : -1 : 1 : 0), (1 : -1 : 1 : 1), (1 : 0 : -1 : 1), (1 : 0 : 0 : 0),
        (1 : 0 : 0 : 1), (1 : 0 : 1 : 0), (1 : 0 : 1 : 1), (1 : 1 : -1 : 1),
        (1 : 1 : 0 : 0), (1 : 1 : 0 : 1), (1 : 1 : 1 : 0), (1 : 1 : 1 : 1)]

    ALGORITHM:

    We just check all possible projective points in correct dimension
    of projective space to see if they lie on ``X``.

    AUTHORS:

    - John Cremona and Charlie Turner (06-2010)
    """
    from sage.schemes.projective.projective_space import is_ProjectiveSpace
    if(is_Scheme(X)):
        if (not is_ProjectiveSpace(X.ambient_space())):
            raise TypeError("Ambient space must be projective space over the rational field")
        X = X(X.base_ring())
    else:
        if (not is_ProjectiveSpace(X.codomain().ambient_space())):
            raise TypeError("Codomain must be projective space over the rational field")

    n = X.codomain().ambient_space().ngens()
    zero = (0,) * n
    pts = []
    for c in cartesian_product_iterator([srange(-B,B+1) for _ in range(n)]):
        if gcd(c) == 1 and c > zero:
            try:
                pts.append(X(c))
            except TypeError:
                pass
    pts.sort()
    return pts


def enum_projective_number_field(X,B, prec=53):
    """
    Enumerates projective points on scheme ``X`` defined over a number field. Simply checks all of the
    points of absolute height of at most ``B`` and adds those that are on the scheme to the list.

    INPUT:

    - ``X`` - a scheme defined over a number field

    - ``B`` - a real number

    - ``prec`` - the precision to use for computing the elements of bounded height of number fields

    OUTPUT:

     - a list containing the projective points of ``X`` of absolute height up to ``B``,
       sorted.

    .. WARNING::

       In the current implementation, the output of the [Doyle-Krumm] algorithm
       for elements of bounded height cannot be guaranteed to be correct due to
       the necessity of floating point computations. In some cases, the default
       53-bit precision is considerably lower than would be required for the
       algorithm to generate correct output.

    EXAMPLES::

        sage: from sage.schemes.projective.projective_rational_point import enum_projective_number_field
        sage: u = QQ['u'].0
        sage: K = NumberField(u^3 - 5,'v')
        sage: P.<x,y,z> = ProjectiveSpace(K, 2)
        sage: X = P.subscheme([x - y])
        sage: enum_projective_number_field(X(K), 5^(1/3), prec=2^10)
        [(0 : 0 : 1), (-1 : -1 : 1), (1 : 1 : 1), (-1/5*v^2 : -1/5*v^2 : 1), (-v : -v : 1),
        (1/5*v^2 : 1/5*v^2 : 1), (v : v : 1), (1 : 1 : 0)]

    ::

        sage: u = QQ['u'].0
        sage: K = NumberField(u^2 + 3, 'v')
        sage: A.<x,y> = ProjectiveSpace(K,1)
        sage: X=A.subscheme(x-y)
        sage: from sage.schemes.projective.projective_rational_point import enum_projective_number_field
        sage: enum_projective_number_field(X, 2)
        [(1 : 1)]
    """
    from sage.schemes.projective.projective_space import is_ProjectiveSpace
    if(is_Scheme(X)):
        if (not is_ProjectiveSpace(X.ambient_space())):
            raise TypeError("Ambient space must be projective space over a number field")
        X = X(X.base_ring())
    else:
        if (not is_ProjectiveSpace(X.codomain().ambient_space())):
            raise TypeError("Codomain must be projective space over a number field")

    R = X.codomain().ambient_space()

    pts = []

    for P in R.points_of_bounded_height(B, prec):
        try:
            pts.append(X(P))
        except TypeError:
            pass
    pts.sort()
    return pts


def enum_projective_finite_field(X):
    """
    Enumerates projective points on scheme ``X`` defined over a finite field.

    INPUT:

    - ``X`` -  a scheme defined over a finite field or a set of abstract
      rational points of such a scheme.

    OUTPUT:

    - a list containing the projective points of ``X`` over the finite field,
      sorted.

    EXAMPLES::

        sage: F = GF(53)
        sage: P.<X,Y,Z> = ProjectiveSpace(2,F)
        sage: from sage.schemes.projective.projective_rational_point import enum_projective_finite_field
        sage: len(enum_projective_finite_field(P(F)))
        2863
        sage: 53^2+53+1
        2863

    ::

        sage: F = GF(9,'a')
        sage: P.<X,Y,Z> = ProjectiveSpace(2,F)
        sage: C = Curve(X^3-Y^3+Z^2*Y)
        sage: enum_projective_finite_field(C(F))
        [(0 : 0 : 1), (0 : 1 : 1), (0 : 2 : 1), (1 : 1 : 0), (a + 1 : 2*a : 1),
        (a + 1 : 2*a + 1 : 1), (a + 1 : 2*a + 2 : 1), (2*a + 2 : a : 1),
        (2*a + 2 : a + 1 : 1), (2*a + 2 : a + 2 : 1)]

    ::

        sage: F = GF(5)
        sage: P2F.<X,Y,Z> = ProjectiveSpace(2,F)
        sage: enum_projective_finite_field(P2F)
        [(0 : 0 : 1), (0 : 1 : 0), (0 : 1 : 1), (0 : 2 : 1), (0 : 3 : 1), (0 : 4 : 1),
        (1 : 0 : 0), (1 : 0 : 1), (1 : 1 : 0), (1 : 1 : 1), (1 : 2 : 1), (1 : 3 : 1),
        (1 : 4 : 1), (2 : 0 : 1), (2 : 1 : 0), (2 : 1 : 1), (2 : 2 : 1), (2 : 3 : 1),
        (2 : 4 : 1), (3 : 0 : 1), (3 : 1 : 0), (3 : 1 : 1), (3 : 2 : 1), (3 : 3 : 1),
        (3 : 4 : 1), (4 : 0 : 1), (4 : 1 : 0), (4 : 1 : 1), (4 : 2 : 1), (4 : 3 : 1),
        (4 : 4 : 1)]

    ALGORITHM:

    Checks all points in projective space to see if they lie on X.

    .. WARNING::

        If ``X`` is defined over an infinite field, this code will not finish!

    AUTHORS:

    - John Cremona and Charlie Turner (06-2010).
    """
    from sage.schemes.projective.projective_space import is_ProjectiveSpace
    if(is_Scheme(X)):
        if (not is_ProjectiveSpace(X.ambient_space())):
            raise TypeError("Ambient space must be projective space over a finite")
        X = X(X.base_ring())
    else:
        if (not is_ProjectiveSpace(X.codomain().ambient_space())):
            raise TypeError("Codomain must be projective space over a finite field")

    n = X.codomain().ambient_space().ngens()-1
    F = X.value_ring()
    pts = []
    for k in range(n+1):
        for c in cartesian_product_iterator([F for _ in range(k)]):
            try:
                pts.append(X(list(c)+[1]+[0]*(n-k)))
            except TypeError:
                pass
    pts.sort()
    return pts

