r"""
Enumeration of rational points on product projective schemes

Naive algorithms for enumerating rational points over `\QQ`, number fields or
finite fields over general schemes.

.. WARNING::

    Incorrect results and infinite loops may occur if using a wrong function.
    (For instance using an affine function for a product projective scheme
    or a finite field function for a scheme defined over an infinite field.)

EXAMPLES:

Product Projective, over `\QQ`::

    sage: PP.<x,y,z> = ProductProjectiveSpaces([1,0], QQ)
    sage: from sage.schemes.product_projective.rational_point import \
            enum_product_projective_rational_field
    sage: enum_product_projective_rational_field(PP,3)
    [(-3 : 1 , 1), (-2 : 1 , 1), (-3/2 : 1 , 1),
     (-1 : 1 , 1), (-2/3 : 1 , 1), (-1/2 : 1 , 1),
     (-1/3 : 1 , 1), (0 : 1 , 1), (1/3 : 1 , 1),
     (1/2 : 1 , 1), (2/3 : 1 , 1), (1 : 0 , 1),
     (1 : 1 , 1), (3/2 : 1 , 1), (2 : 1 , 1),
     (3 : 1 , 1)]

Product projective over finite field::

    sage: P1.<x,y,a,b> = ProductProjectiveSpaces([1,1], GF(7))
    sage: X = P1.subscheme([2*x+3*y])
    sage: from sage.schemes.product_projective.rational_point import \
            enum_product_projective_finite_field
    sage: enum_product_projective_finite_field(X)
    [(2 : 1 , 0 : 1), (2 : 1 , 1 : 0), (2 : 1 , 1 : 1),
     (2 : 1 , 2 : 1), (2 : 1 , 3 : 1), (2 : 1 , 4 : 1),
     (2 : 1 , 5 : 1), (2 : 1 , 6 : 1)]

AUTHORS:

- Volker Braun and Ben Hutz (2014): initial version

- Raghukul Raman (2018): code cleanup and added support for rational fields

"""

# ****************************************************************************
# Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#                    Ben Hutz <bn4941@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.schemes.generic.scheme import is_Scheme
from sage.schemes.product_projective.space import is_ProductProjectiveSpaces


def enum_product_projective_rational_field(X, B):
    r"""
    Enumerate projective, rational points on scheme ``X`` of height up to
    bound ``B``.

    INPUT:

    - ``X`` -- a scheme or set of abstract rational points of a scheme

    - ``B`` -- a positive integer bound

    OUTPUT:

    - a list containing the product projective points of ``X`` of height up
      to ``B``, sorted.

    EXAMPLES::

        sage: PP.<x0,x1,x2,x3,x4> = ProductProjectiveSpaces([1, 2], QQ)
        sage: from sage.schemes.product_projective.rational_point import \
                enum_product_projective_rational_field
        sage: enum_product_projective_rational_field(PP,1)
        [(-1 : 1 , -1 : -1 : 1), (-1 : 1 , -1 : 0 : 1), (-1 : 1 , -1 : 1 : 0),
         (-1 : 1 , -1 : 1 : 1), (-1 : 1 , 0 : -1 : 1), (-1 : 1 , 0 : 0 : 1),
         (-1 : 1 , 0 : 1 : 0), (-1 : 1 , 0 : 1 : 1), (-1 : 1 , 1 : -1 : 1),
         (-1 : 1 , 1 : 0 : 0), (-1 : 1 , 1 : 0 : 1), (-1 : 1 , 1 : 1 : 0),
         (-1 : 1 , 1 : 1 : 1), (0 : 1 , -1 : -1 : 1), (0 : 1 , -1 : 0 : 1),
         (0 : 1 , -1 : 1 : 0), (0 : 1 , -1 : 1 : 1), (0 : 1 , 0 : -1 : 1),
         (0 : 1 , 0 : 0 : 1), (0 : 1 , 0 : 1 : 0), (0 : 1 , 0 : 1 : 1),
         (0 : 1 , 1 : -1 : 1), (0 : 1 , 1 : 0 : 0), (0 : 1 , 1 : 0 : 1),
         (0 : 1 , 1 : 1 : 0), (0 : 1 , 1 : 1 : 1), (1 : 0 , -1 : -1 : 1),
         (1 : 0 , -1 : 0 : 1), (1 : 0 , -1 : 1 : 0), (1 : 0 , -1 : 1 : 1),
         (1 : 0 , 0 : -1 : 1), (1 : 0 , 0 : 0 : 1), (1 : 0 , 0 : 1 : 0),
         (1 : 0 , 0 : 1 : 1), (1 : 0 , 1 : -1 : 1), (1 : 0 , 1 : 0 : 0),
         (1 : 0 , 1 : 0 : 1), (1 : 0 , 1 : 1 : 0), (1 : 0 , 1 : 1 : 1),
         (1 : 1 , -1 : -1 : 1), (1 : 1 , -1 : 0 : 1), (1 : 1 , -1 : 1 : 0),
         (1 : 1 , -1 : 1 : 1), (1 : 1 , 0 : -1 : 1), (1 : 1 , 0 : 0 : 1),
         (1 : 1 , 0 : 1 : 0), (1 : 1 , 0 : 1 : 1), (1 : 1 , 1 : -1 : 1),
         (1 : 1 , 1 : 0 : 0), (1 : 1 , 1 : 0 : 1), (1 : 1 , 1 : 1 : 0),
         (1 : 1 , 1 : 1 : 1)]

    ::

        sage: PP.<x,y,z,u,v> = ProductProjectiveSpaces([2,1], QQ)
        sage: X = PP.subscheme([x^2 + x*y + y*z, u*u-v*u])
        sage: from sage.schemes.product_projective.rational_point import \
                enum_product_projective_rational_field
        sage: enum_product_projective_rational_field(X,4)
        [(-2 : 4 : 1 , 0 : 1), (-2 : 4 : 1 , 1 : 1), (-1 : 1 : 0 , 0 : 1),
         (-1 : 1 : 0 , 1 : 1), (-2/3 : -4/3 : 1 , 0 : 1), (-2/3 : -4/3 : 1 , 1 : 1),
         (-1/2 : -1/2 : 1 , 0 : 1), (-1/2 : -1/2 : 1 , 1 : 1),
         (0 : 0 : 1 , 0 : 1), (0 : 0 : 1 , 1 : 1), (0 : 1 : 0 , 0 : 1),
         (0 : 1 : 0 , 1 : 1), (1 : -1/2 : 1 , 0 : 1), (1 : -1/2 : 1 , 1 : 1)]
    """
    if(is_Scheme(X)):
        if (not is_ProductProjectiveSpaces(X.ambient_space())):
            raise TypeError("ambient space must be product of projective space over the rational field")
        X = X(X.base_ring())
    else:
        if (not is_ProductProjectiveSpaces(X.codomain().ambient_space())):
            raise TypeError("codomain must be product of projective space over the rational field")

    R = X.codomain().ambient_space()
    m = R.num_components()
    iters = [ R[i].points_of_bounded_height(bound=B) for i in range(m) ]
    dim = [R[i].dimension_relative() + 1 for i in range(m)]
    
    dim_prefix = [0, dim[0]] # prefixes dim list
    for i in range(1, len(dim)):
        dim_prefix.append(dim_prefix[i] + dim[i])

    pts = []
    P = []
    for i in range(m):
        pt = next(iters[i])
        for j in range(dim[i]):
            P.append(pt[j]) # initial value of P

    try: # add the initial point
        pts.append(X(P))
    except TypeError:
        pass

    i = 0
    while i < m:
        try:
            pt = next(iters[i])
            for j in range(dim[i]):
                P[dim_prefix[i] + j] = pt[j]
            try:
                pts.append(X(P))
            except TypeError:
                pass
            i = 0
        except StopIteration:
            iters[i] = R[i].points_of_bounded_height(bound=B)
            pt = next(iters[i]) # reset
            for j in range(dim[i]):
                P[dim_prefix[i] + j] = pt[j]
            i += 1
    pts.sort()

    return pts

def enum_product_projective_number_field(X, **kwds):
    r"""
    Enumerates product projective points on scheme ``X`` defined over a number field.

    Simply checks all of the points of absolute height of at most ``B``
    and adds those that are on the scheme to the list.

    This algorithm computes 2 lists: L containing elements x in `K` such that
    H_k(x) <= B, and a list L' containing elements x in `K` that, due to
    floating point issues,
    may be slightly larger then the bound. This can be controlled
    by lowering the tolerance.

    ALGORITHM:

    This is an implementation of the revised algorithm (Algorithm 4) in
    [DK2013]_. Algorithm 5 is used for imaginary quadratic fields.
    
    INPUT:

    kwds:

    - ``bound`` - a real number

    - ``tolerance`` - a rational number in (0,1] used in doyle-krumm algorithm-4

    - ``precision`` - the precision to use for computing the elements of bounded height of number fields.

    OUTPUT:

    - a list containing the product projective points of ``X`` of
      absolute height up to ``B``, sorted.

    EXAMPLES::

        sage: u = QQ['u'].0
        sage: K = NumberField(u^2 + 2, 'v')
        sage: PP.<x,y,z,w> = ProductProjectiveSpaces([1, 1], K)
        sage: X = PP.subscheme([x^2 + 2*y^2])
        sage: from sage.schemes.product_projective.rational_point import \
                enum_product_projective_number_field
        sage: enum_product_projective_number_field(X, bound=1.5)
        [(-v : 1 , -1 : 1), (-v : 1 , -v : 1), (-v : 1 , -1/2*v : 1),
         (-v : 1 , 0 : 1), (-v : 1 , 1/2*v : 1), (-v : 1 , v : 1),
         (-v : 1 , 1 : 0), (-v : 1 , 1 : 1), (v : 1 , -1 : 1),
         (v : 1 , -v : 1), (v : 1 , -1/2*v : 1), (v : 1 , 0 : 1),
         (v : 1 , 1/2*v : 1), (v : 1 , v : 1), (v : 1 , 1 : 0),
         (v : 1 , 1 : 1)]
    """
    B = kwds.pop('bound')
    tol = kwds.pop('tolerance', 1e-2)
    prec = kwds.pop('precision', 53)

    if(is_Scheme(X)):
        if (not is_ProductProjectiveSpaces(X.ambient_space())):
            raise TypeError("ambient space must be product of projective space over the rational field")
        X = X(X.base_ring())
    else:
        if (not is_ProductProjectiveSpaces(X.codomain().ambient_space())):
            raise TypeError("codomain must be product of projective space over the rational field")

    R = X.codomain().ambient_space()

    pts = []

    for P in R.points_of_bounded_height(bound=B, tolerance=tol, precision=prec):
        try:
            pts.append(X(P))
        except TypeError:
            pass
    pts.sort()
    return pts

def enum_product_projective_finite_field(X):
    r"""
    Enumerates projective points on scheme ``X`` defined over a finite field.

    INPUT:

    - ``X`` -  a scheme defined over a finite field or a set of abstract
      rational points of such a scheme.

    OUTPUT:

    - a list containing the projective points of ``X`` over the finite field,
      sorted.

    EXAMPLES::

        sage: PP.<x,y,z,w> = ProductProjectiveSpaces([1, 1], GF(3))
        sage: from sage.schemes.product_projective.rational_point import \
                enum_product_projective_finite_field
        sage: enum_product_projective_finite_field(PP)
        [(0 : 1 , 0 : 1), (0 : 1 , 1 : 0), (0 : 1 , 1 : 1),
         (0 : 1 , 2 : 1), (1 : 0 , 0 : 1), (1 : 0 , 1 : 0),
         (1 : 0 , 1 : 1), (1 : 0 , 2 : 1), (1 : 1 , 0 : 1),
         (1 : 1 , 1 : 0), (1 : 1 , 1 : 1), (1 : 1 , 2 : 1),
         (2 : 1 , 0 : 1), (2 : 1 , 1 : 0), (2 : 1 , 1 : 1),
         (2 : 1 , 2 : 1)]

    ::

        sage: PP.<x0,x1,x2,x3> = ProductProjectiveSpaces([1, 1], GF(17))
        sage: X = PP.subscheme([x0^2 + 2*x1^2])
        sage: from sage.schemes.product_projective.rational_point import \
                enum_product_projective_finite_field
        sage: len(enum_product_projective_finite_field(X))
        36
    """
    if(is_Scheme(X)):
        if (not is_ProductProjectiveSpaces(X.ambient_space())):
            raise TypeError("ambient space must be product of projective space over the rational field")
        X = X(X.base_ring())
    else:
        if (not is_ProductProjectiveSpaces(X.codomain().ambient_space())):
            raise TypeError("codomain must be product of projective space over the rational field")

    R = X.codomain().ambient_space()
    pts = []

    for P in R.rational_points():
        try:
            pts.append(X(P))
        except TypeError:
            pass
    pts.sort()

    return pts
