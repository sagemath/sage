r"""
Enumeration of rational points on affine schemes

Naive algorithms for enumerating rational points over `\QQ` or finite fields
over for general schemes.

.. WARNING::

    Incorrect results and infinite loops may occur if using a wrong function.
    (For instance using an affine function for a projective scheme or a finite
    field function for a scheme defined over an infinite field.)

EXAMPLES:

Affine, over `\QQ`::

    sage: from sage.schemes.affine.affine_rational_point import enum_affine_rational_field
    sage: A.<x,y,z> = AffineSpace(3,QQ)
    sage: S = A.subscheme([2*x-3*y])
    sage: enum_affine_rational_field(S,2)
    [(0, 0, -2), (0, 0, -1), (0, 0, -1/2), (0, 0, 0),
     (0, 0, 1/2), (0, 0, 1), (0, 0, 2)]

Affine over a finite field::

    sage: from sage.schemes.affine.affine_rational_point import enum_affine_finite_field
    sage: A.<w,x,y,z> = AffineSpace(4,GF(2))
    sage: enum_affine_finite_field(A(GF(2)))
    [(0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0), (0, 0, 1, 1), (0, 1, 0, 0),
     (0, 1, 0, 1), (0, 1, 1, 0), (0, 1, 1, 1), (1, 0, 0, 0), (1, 0, 0, 1),
     (1, 0, 1, 0), (1, 0, 1, 1), (1, 1, 0, 0), (1, 1, 0, 1), (1, 1, 1, 0),
     (1, 1, 1, 1)]

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


from sage.rings.all import ZZ
from sage.misc.all import srange, cartesian_product_iterator
from sage.schemes.generic.scheme import is_Scheme


def enum_affine_rational_field(X,B):
    """
    Enumerates affine rational points on scheme ``X`` (defined over `\QQ`) up
    to bound ``B``.

    INPUT:

    - ``X`` -  a scheme or set of abstract rational points of a scheme;
    - ``B`` -  a positive integer bound.

    OUTPUT:

    - a list containing the affine points of ``X`` of height up to ``B``,
      sorted.

    EXAMPLES::

        sage: A.<x,y,z> = AffineSpace(3,QQ)
        sage: from sage.schemes.affine.affine_rational_point import enum_affine_rational_field
        sage: enum_affine_rational_field(A(QQ),1)
        [(-1, -1, -1), (-1, -1, 0), (-1, -1, 1), (-1, 0, -1), (-1, 0, 0), (-1, 0, 1),
        (-1, 1, -1), (-1, 1, 0), (-1, 1, 1), (0, -1, -1), (0, -1, 0), (0, -1, 1),
        (0, 0, -1), (0, 0, 0), (0, 0, 1), (0, 1, -1), (0, 1, 0), (0, 1, 1), (1, -1, -1),
        (1, -1, 0), (1, -1, 1), (1, 0, -1), (1, 0, 0), (1, 0, 1), (1, 1, -1), (1, 1, 0),
        (1, 1, 1)]

    ::

        sage: A.<w,x,y,z> = AffineSpace(4,QQ)
        sage: S = A.subscheme([x^2-y*z+3,w^3+z+y^2])
        sage: enum_affine_rational_field(S(QQ),2)
        []
        sage: enum_affine_rational_field(S(QQ),3)
        [(-2, 0, -3, -1)]

    ::

        sage: A.<x,y> = AffineSpace(2,QQ)
        sage: C = Curve(x^2+y-x)
        sage: enum_affine_rational_field(C,10)
        [(-2, -6), (-1, -2), (0, 0), (1, 0), (2, -2), (3, -6)]


    AUTHORS:

    - David R. Kohel <kohel@maths.usyd.edu.au>: original version.

    - Charlie Turner (06-2010): small adjustments.
    """
    from sage.schemes.affine.affine_space import is_AffineSpace
    if(is_Scheme(X)):
        if (not is_AffineSpace(X.ambient_space())):
            raise TypeError("Ambient space must be affine space over the rational field")
        X = X(X.base_ring())
    else:
        if (not is_AffineSpace(X.codomain().ambient_space())):
            raise TypeError("Codomain must be affine space over the rational field")

    n = X.codomain().ambient_space().ngens()
    if X.value_ring() is ZZ:
        Q = [ 1 ]
    else: # rational field
        Q = range(1, B + 1)
    R = [ 0 ] + [ s*k for k in range(1, B+1) for s in [1, -1] ]
    pts = []
    P = [0] * n
    m = ZZ(0)
    try:
        pts.append(X(P))
    except TypeError:
        pass
    iters = [ iter(R) for _ in range(n) ]
    for it in iters:
        next(it)
    i = 0
    while i < n:
        try:
            a = ZZ(next(iters[i]))
        except StopIteration:
            iters[i] = iter(R) # reset
            P[i] = next(iters[i]) # reset P[i] to 0 and increment
            i += 1
            continue
        m = m.gcd(a)
        P[i] = a
        for b in Q:
            if m.gcd(b) == 1:
                try:
                    pts.append(X([ num/b for num in P ]))
                except TypeError:
                    pass
        i = 0
        m = ZZ(0)
    pts.sort()
    return pts


def enum_affine_number_field(X,B):
    """
    Enumerates affine points on scheme ``X`` defined over a number field. Simply checks all of the
    points of absolute height up to ``B`` and adds those that are on the scheme to the list.

    INPUT:

    - ``X`` - a scheme defined over a number field

    - ``B`` - a real number

    OUTPUT:

     - a list containing the affine points of ``X`` of absolute height up to ``B``,
       sorted.

    EXAMPLES::

        sage: from sage.schemes.affine.affine_rational_point import enum_affine_number_field
        sage: u = QQ['u'].0
        sage: K = NumberField(u^2 + 2,'v')
        sage: A.<x,y,z> = AffineSpace(K, 3)
        sage: X = A.subscheme([y^2 - x])
        sage: enum_affine_number_field(X(K),4)
        [(0, 0, -1), (0, 0, -v), (0, 0, -1/2*v), (0, 0, 0), (0, 0, 1/2*v), (0, 0, v), (0, 0, 1),
        (1, -1, -1), (1, -1, -v), (1, -1, -1/2*v), (1, -1, 0), (1, -1, 1/2*v), (1, -1, v), (1, -1, 1),
        (1, 1, -1), (1, 1, -v), (1, 1, -1/2*v), (1, 1, 0), (1, 1, 1/2*v), (1, 1, v), (1, 1, 1)]

    ::

        sage: u = QQ['u'].0
        sage: K = NumberField(u^2 + 3,'v')
        sage: A.<x,y> = AffineSpace(K,2)
        sage: X=A.subscheme(x-y)
        sage: from sage.schemes.affine.affine_rational_point import enum_affine_number_field
        sage: enum_affine_number_field(X,3)
        [(-1, -1), (-1/2*v - 1/2, -1/2*v - 1/2), (1/2*v - 1/2, 1/2*v - 1/2), (0, 0), (-1/2*v + 1/2, -1/2*v + 1/2),
        (1/2*v + 1/2, 1/2*v + 1/2), (1, 1)]
    """
    from sage.schemes.affine.affine_space import is_AffineSpace
    if(is_Scheme(X)):
        if (not is_AffineSpace(X.ambient_space())):
            raise TypeError("Ambient space must be affine space over a number field")
        X = X(X.base_ring())
    else:
        if (not is_AffineSpace(X.codomain().ambient_space())):
            raise TypeError("Codomain must be affine space over a number field")

    R = X.codomain().ambient_space()

    pts = []
    for P in R.points_of_bounded_height(B):
        try:
            pts.append(X(P))
        except TypeError:
            pass
    pts.sort()
    return pts


def enum_affine_finite_field(X):
    r"""
    Enumerates affine points on scheme ``X`` defined over a finite field.

    INPUT:

    - ``X`` -  a scheme defined over a finite field or a set of abstract
      rational points of such a scheme.

    OUTPUT:

    - a list containing the affine points of ``X`` over the finite field,
      sorted.

    EXAMPLES::

        sage: F = GF(7)
        sage: A.<w,x,y,z> = AffineSpace(4,F)
        sage: C = A.subscheme([w^2+x+4,y*z*x-6,z*y+w*x])
        sage: from sage.schemes.affine.affine_rational_point import enum_affine_finite_field
        sage: enum_affine_finite_field(C(F))
        []
        sage: C = A.subscheme([w^2+x+4,y*z*x-6])
        sage: enum_affine_finite_field(C(F))
        [(0, 3, 1, 2), (0, 3, 2, 1), (0, 3, 3, 3), (0, 3, 4, 4), (0, 3, 5, 6),
        (0, 3, 6, 5), (1, 2, 1, 3), (1, 2, 2, 5), (1, 2, 3, 1), (1, 2, 4, 6),
        (1, 2, 5, 2), (1, 2, 6, 4), (2, 6, 1, 1), (2, 6, 2, 4), (2, 6, 3, 5),
        (2, 6, 4, 2), (2, 6, 5, 3), (2, 6, 6, 6), (3, 1, 1, 6), (3, 1, 2, 3),
        (3, 1, 3, 2), (3, 1, 4, 5), (3, 1, 5, 4), (3, 1, 6, 1), (4, 1, 1, 6),
        (4, 1, 2, 3), (4, 1, 3, 2), (4, 1, 4, 5), (4, 1, 5, 4), (4, 1, 6, 1),
        (5, 6, 1, 1), (5, 6, 2, 4), (5, 6, 3, 5), (5, 6, 4, 2), (5, 6, 5, 3),
        (5, 6, 6, 6), (6, 2, 1, 3), (6, 2, 2, 5), (6, 2, 3, 1), (6, 2, 4, 6),
        (6, 2, 5, 2), (6, 2, 6, 4)]

    ::

        sage: A.<x,y,z> = AffineSpace(3,GF(3))
        sage: S = A.subscheme(x+y)
        sage: enum_affine_finite_field(S)
        [(0, 0, 0), (0, 0, 1), (0, 0, 2), (1, 2, 0), (1, 2, 1), (1, 2, 2),
        (2, 1, 0), (2, 1, 1), (2, 1, 2)]

    ALGORITHM:

    Checks all points in affine space to see if they lie on X.

    .. WARNING::

        If ``X`` is defined over an infinite field, this code will not finish!

    AUTHORS:

    - John Cremona and Charlie Turner (06-2010)
    """
    from sage.schemes.affine.affine_space import is_AffineSpace
    if(is_Scheme(X)):
        if (not is_AffineSpace(X.ambient_space())):
            raise TypeError("Ambient space must be affine space over a finite field")
        X = X(X.base_ring())
    else:
        if (not is_AffineSpace(X.codomain().ambient_space())):
            raise TypeError("Codomain must be affine space over a finite field")

    n = X.codomain().ambient_space().ngens()
    F = X.value_ring()
    pts = []
    for c in cartesian_product_iterator([F]*n):
        try:
            pts.append(X(c))
        except Exception:
            pass
    pts.sort()
    return pts
