r"""
Products of projective spaces

This class builds on the projective space class and its point and morphism
classes.

Products of projective spaces of varying dimension are convenient ambient
spaces for complete intersections.

Group actions on them, and the interplay with representation theory, provide
many interesting examples of algebraic varieties.

EXAMPLES:

We construct products projective spaces of various dimensions over the same ring::

    sage: P1 = ProjectiveSpace(ZZ, 1, 'x')
    sage: P2 = ProjectiveSpace(ZZ, 2, 'y')
    sage: ProductProjectiveSpaces([P1, P2])
    Product of projective spaces P^1 x P^2 over Integer Ring

We can also construct the product by specifying the dimensions and the base ring::

    sage: ProductProjectiveSpaces([1, 2, 3], QQ, 'z')
    Product of projective spaces P^1 x P^2 x P^3 over Rational Field

    sage: P2xP2 = ProductProjectiveSpaces([2, 2], QQ, names=['x', 'y'])
    sage: P2xP2.coordinate_ring().inject_variables()
    Defining x0, x1, x2, y0, y1, y2
"""

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2014 Ben Hutz <bn4941@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.categories.fields import Fields
from sage.categories.commutative_rings import CommutativeRings
from sage.rings.polynomial.polydict import ETuple
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
from sage.schemes.generic.ambient_space import AmbientSpace
from sage.schemes.projective.projective_space import ProjectiveSpace, ProjectiveSpace_ring
from sage.schemes.product_projective.homset import (SchemeHomset_points_product_projective_spaces_ring,
                                                    SchemeHomset_points_product_projective_spaces_field)
from sage.schemes.product_projective.point import (ProductProjectiveSpaces_point_ring,
                                                   ProductProjectiveSpaces_point_field,
                                                   ProductProjectiveSpaces_point_finite_field)
from sage.schemes.product_projective.morphism import ProductProjectiveSpaces_morphism_ring
from sage.schemes.product_projective.subscheme import AlgebraicScheme_subscheme_product_projective


def is_ProductProjectiveSpaces(x):
    r"""
    Return ``True`` if ``x`` is a product of projective spaces.

    This is an ambient space defined by `\mathbb{P}^n_R \times \cdots \times \mathbb{P}^m_R`,
    where `R` is a ring and `n,\ldots, m\geq 0` are integers.

    OUTPUT: boolean

    EXAMPLES::

        sage: is_ProductProjectiveSpaces(ProjectiveSpace(5, names='x'))
        doctest:warning...
        DeprecationWarning: The function is_ProductProjectiveSpaces is deprecated; use 'isinstance(..., ProductProjectiveSpaces_ring)' instead.
        See https://github.com/sagemath/sage/issues/38022 for details.
        False
        sage: is_ProductProjectiveSpaces(ProductProjectiveSpaces([1, 2, 3], ZZ, 'x'))
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38022, "The function is_ProductProjectiveSpaces is deprecated; use 'isinstance(..., ProductProjectiveSpaces_ring)' instead.")
    return isinstance(x, ProductProjectiveSpaces_ring)


def ProductProjectiveSpaces(n, R=None, names='x'):
    r"""
    Return the Cartesian product of projective spaces.

    The input ``n`` is either a list of projective space over the same base
    ring or the list of dimensions, ``R`` the base ring, and ``names`` the
    variable names.

    INPUT:

    - ``n`` -- list of integers or a list of projective spaces

    - ``R`` -- a ring

    - ``names`` -- string or list of strings

    EXAMPLES::

        sage: P1 = ProjectiveSpace(QQ, 2, 'x')
        sage: P2 = ProjectiveSpace(QQ, 3, 'y')
        sage: ProductProjectiveSpaces([P1, P2])
        Product of projective spaces P^2 x P^3 over Rational Field

    ::

        sage: ProductProjectiveSpaces([2, 2], GF(7), 'y')
        Product of projective spaces P^2 x P^2 over Finite Field of size 7

    ::

        sage: P1 = ProjectiveSpace(ZZ, 2, 'x')
        sage: P2 = ProjectiveSpace(QQ, 3, 'y')
        sage: ProductProjectiveSpaces([P1, P2])
        Traceback (most recent call last):
        ...
        AttributeError: components must be over the same base ring
    """
    if isinstance(R, (list, tuple)):
        n, R = R, n
    if not isinstance(n, (tuple, list)):
        raise TypeError("must be a list of dimensions")
    if R is None:
        R = QQ  # default is the rationals
    if isinstance(n[0], ProjectiveSpace_ring):
        # this should be a list of projective spaces
        names = []
        N = []
        R = None
        for PS in n:
            if not isinstance(PS,ProjectiveSpace_ring):
                raise TypeError("must be a list of projective spaces or (dimensions, base ring, names)")
            if R is None:
                R = PS.base_ring()
            elif R != PS.base_ring():
                raise AttributeError("components must be over the same base ring")
            N.append(PS.dimension_relative())
            names += PS.variable_names()
        if isinstance(R, FiniteField):
            X = ProductProjectiveSpaces_finite_field(N, R, names)
        elif R in Fields():
            X = ProductProjectiveSpaces_field(N, R, names)
        else:
            X = ProductProjectiveSpaces_ring(N, R, names)
        X._components = n
    else:
        if not isinstance(n,(list, tuple)):
            raise ValueError("need list or tuple of dimensions")
        if R not in CommutativeRings():
            raise ValueError("must be a commutative ring")
        from sage.structure.category_object import normalize_names
        n_vars = sum(d + 1 for d in n)
        if isinstance(names, str):
            names = normalize_names(n_vars, names)
        else:
            name_list = list(names)
            if len(name_list) == len(n):
                names = []
                for name, dim in zip(name_list, n):
                    names += normalize_names(dim + 1, name)
            else:
                n_vars = sum(1+d for d in n)
                names = normalize_names(n_vars, name_list)
        if isinstance(R, FiniteField):
            X = ProductProjectiveSpaces_finite_field(n, R, names)
        elif R in Fields():
            X = ProductProjectiveSpaces_field(n, R, names)
        else:
            X = ProductProjectiveSpaces_ring(n, R, names)
    return X


class ProductProjectiveSpaces_ring(AmbientSpace):
    r"""
    Cartesian product of projective spaces `\mathbb{P}^{n_1} \times \cdots \times \mathbb{P}^{n_r}`.

    EXAMPLES::

        sage: P.<x0,x1,x2,x3,x4> = ProductProjectiveSpaces([1, 2], QQ); P
        Product of projective spaces P^1 x P^2 over Rational Field
        sage: P.coordinate_ring()
        Multivariate Polynomial Ring in x0, x1, x2, x3, x4 over Rational Field
        sage: P[0]
        Projective Space of dimension 1 over Rational Field
        sage: P[1]
        Projective Space of dimension 2 over Rational Field
        sage: Q = P(6, 3, 2, 2, 2); Q
        (2 : 1 , 1 : 1 : 1)
        sage: Q[0]
        (2 : 1)
        sage: H = Hom(P,P)
        sage: f = H([x0^2*x3, x2*x1^2, x2^2, 2*x3^2, x4^2])
        sage: f(Q)
        (4 : 1 , 1 : 2 : 1)
    """
    def __init__(self, N, R=QQ, names=None):
        r"""
        The Python constructor.

        INPUT:

        - ``N`` -- list or tuple of positive integers

        - ``R`` -- a ring

        - ``names`` -- tuple or list of strings; this must either be a single
          variable name or the complete list of variables

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([2, 2], QQ)
            sage: T
            Product of projective spaces P^2 x P^2 over Rational Field
            sage: T.coordinate_ring()
            Multivariate Polynomial Ring in x, y, z, u, v, w over Rational Field
            sage: T[1].coordinate_ring()
            Multivariate Polynomial Ring in u, v, w over Rational Field

        ::

            sage: ProductProjectiveSpaces([1,1,1],ZZ, ['x', 'y', 'z', 'u', 'v', 'w'])
            Product of projective spaces P^1 x P^1 x P^1 over Integer Ring

        ::

            sage: T = ProductProjectiveSpaces([1, 1], QQ, 'z')
            sage: T.coordinate_ring()
            Multivariate Polynomial Ring in z0, z1, z2, z3 over Rational Field
        """
        assert isinstance(N, (tuple, list))
        N = [Integer(n) for n in N]
        assert R in CommutativeRings()
        if len(N) < 2:
            raise ValueError("must be at least two components for a product")
        AmbientSpace.__init__(self, sum(N), R)
        self._dims = N
        start = 0
        self._components = []
        for i in range(len(N)):
            self._components.append(ProjectiveSpace(N[i],R,names[start:start+N[i]+1]))
            start += N[i]+1
        # Note that the coordinate ring should really be the tensor product of
        # the component coordinate rings. But we just deal with them as
        # multihomogeneous polynomial rings.
        self._coordinate_ring = PolynomialRing(R,sum(N) + len(N),names)
        self._assign_names(names)

    def _repr_(self):
        r"""
        Return a string representation of this space.

        EXAMPLES::

            sage: ProductProjectiveSpaces([1, 1, 1], ZZ, ['x', 'y', 'z', 'u', 'v', 'w'])
            Product of projective spaces P^1 x P^1 x P^1 over Integer Ring
        """
        return ''.join(['Product of projective spaces ',
                        ' x '.join('P^{}'.format(d) for d in self._dims),
                        ' over ', str(self.base_ring())])

    def _repr_generic_point(self, v=None):
        """
        Return a string representation of the generic point
        on this product space.

        If ``v`` is ``None``, the representation of the generic point of
        the product space is returned.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 2, 1], QQ, 'x')
            sage: T._repr_generic_point()
            '(x0 : x1 , x2 : x3 : x4 , x5 : x6)'
        """
        if v is None:
            v = list(self.gens())
        else:
            v = list(v)
        splitv = self._factors(v)
        return '(%s)' % (" , ".join((" : ".join(str(t) for t in P))
                                    for P in splitv))

    def _latex_(self):
        r"""
        Return a LaTeX representation of this product space.

        EXAMPLES::

            sage: latex(ProductProjectiveSpaces([1, 2, 3], ZZ, 'x'))
            {\mathbf P}_{\Bold{Z}}^{1} \times {\mathbf P}_{\Bold{Z}}^{2} \times
             {\mathbf P}_{\Bold{Z}}^{3}
        """
        return " \\times ".join(PS._latex_() for PS in self)

    def _latex_generic_point(self, v=None):
        """
        Return a LaTeX representation of the generic point
        on this product space.

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 2, 3], ZZ, 'x')
            sage: T._latex_generic_point()
            '\\left(x_{0} : x_{1} , x_{2} : x_{3} : x_{4} , x_{5} : x_{6} : x_{7} : x_{8}\\right)'
        """
        if v is None:
            v = list(self.gens())
        else:
            v = list(v)
        splitv = self._factors(v)
        return '\\left(%s\\right)' % (" , ".join((" : ".join(t._latex_() for t in P)) for P in splitv))

    def __getitem__(self, i):
        r"""
        Return the `i`-th component of the product.

        INPUT:

        - ``i`` -- positive integer

        OUTPUT: a projective space

        EXAMPLES::

            sage: T.<a,x,y,z,u,v,w> = ProductProjectiveSpaces([3, 2], QQ)
            sage: T[0]
            Projective Space of dimension 3 over Rational Field
        """
        return self._components[i]

    def __eq__(self, right):
        r"""
        Check equality of two products of projective spaces.

        INPUT:

        - ``right`` -- a product of projective spaces

        OUTPUT: boolean

        EXAMPLES::

            sage: S.<a,x,y,z,u,v,w> = ProductProjectiveSpaces([3, 2], QQ)
            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([2, 2], QQ)
            sage: S == T
            False
        """
        if not isinstance(right, ProductProjectiveSpaces_ring):
            return False
        else:
            return self._components == right._components

    def __ne__(self, other):
        """
        Check non-equality of two products of projective spaces.

        INPUT:

        - ``other`` -- a product of projective spaces

        OUTPUT: boolean

        EXAMPLES::

            sage: S.<a,x,y,z,u,v,w> = ProductProjectiveSpaces([3, 2], QQ)
            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([2, 2], QQ)
            sage: S != T
            True
        """
        return not (self == other)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: S.<a,x,y,z,u,v,w> = ProductProjectiveSpaces([3, 2], QQ)
            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([2, 2], QQ)
            sage: U.<a,x,y,z,u,v,w> = ProductProjectiveSpaces([3, 2], QQ)
            sage: hash(S) == hash(T)
            False
            sage: hash(S) == hash(U)
            True
        """
        return hash(tuple(self._components))

    def __pow__(self, m):
        """
        Return the Cartesian power of this space.

        INPUT:

        - ``m`` -- integer

        OUTPUT: product of projective spaces

        EXAMPLES::

            sage: P1 = ProductProjectiveSpaces([2, 1], QQ, 'x')
            sage: P1^3
            Product of projective spaces P^2 x P^1 x P^2 x P^1 x P^2 x P^1
             over Rational Field

        As you see, custom variable names are not preserved by power operator,
        since there is no natural way to make new ones in general.
        """
        mm = int(m)
        if mm != m:
            raise ValueError("m must be an integer")
        return ProductProjectiveSpaces(self.dimension_relative_components()*mm, self.base_ring())

    def __mul__(self, right):
        r"""
        Create the product of projective spaces.

        INPUT:

        - ``right`` -- a projective space, product of projective spaces, or subscheme

        OUTPUT: a product of projective spaces or subscheme

        EXAMPLES::

            sage: S.<t,x,y,z,u,v,w> = ProductProjectiveSpaces([3, 2], QQ)
            sage: T.<a,b> = ProjectiveSpace(QQ, 1)
            sage: S*T
            Product of projective spaces P^3 x P^2 x P^1 over Rational Field

        ::

            sage: S = ProductProjectiveSpaces([3, 2], QQ, 'x')
            sage: T = ProductProjectiveSpaces([2, 2], QQ, 'y')
            sage: S*T
            Product of projective spaces P^3 x P^2 x P^2 x P^2 over Rational Field

        ::

            sage: S = ProductProjectiveSpaces([1, 2, 1], ZZ, 't')
            sage: T = ProductProjectiveSpaces([2, 2], ZZ, 'x')
            sage: T.inject_variables()
            Defining x0, x1, x2, x3, x4, x5
            sage: X = T.subscheme([x0*x4 - x1*x3])
            sage: S*X
            Closed subscheme of Product of projective spaces P^1 x P^2 x P^1 x P^2 x P^2 over Integer Ring defined by:
              -x1*x3 + x0*x4

        ::

            sage: S = ProductProjectiveSpaces([3, 2], QQ,'x')
            sage: T = AffineSpace(2, QQ, 'y')
            sage: S*T
            Traceback (most recent call last):
            ...
            TypeError: Affine Space of dimension 2 over Rational Field must be a projective space,
            product of projective spaces, or subscheme
        """
        if self is right:
            return self.__pow__(2)
        if isinstance(right, ProductProjectiveSpaces_ring):
            return ProductProjectiveSpaces(self.components() + right.components())
        elif isinstance(right, ProjectiveSpace_ring):
            return ProductProjectiveSpaces(self.components() + [right])
        elif isinstance(right, AlgebraicScheme_subscheme):
            AS = self*right.ambient_space()
            CR = AS.coordinate_ring()
            n = self.ambient_space().coordinate_ring().ngens()

            phi = self.ambient_space().coordinate_ring().hom(list(CR.gens()[:n]), CR)
            psi = right.ambient_space().coordinate_ring().hom(list(CR.gens()[n:]), CR)
            return AS.subscheme([phi(t) for t in self.defining_polynomials()] + [psi(t) for t in right.defining_polynomials()])
        else:
            raise TypeError('%s must be a projective space, product of projective spaces, or subscheme' % right)

    def components(self):
        r"""
        Return the components of this product of projective spaces.

        OUTPUT: list of projective spaces

        EXAMPLES::

            sage: P.<x,y,z,u,v> = ProductProjectiveSpaces(QQ, [2, 1])
            sage: P.components()
            [Projective Space of dimension 2 over Rational Field,
            Projective Space of dimension 1 over Rational Field]
        """
        return self._components

    def dimension_relative(self):
        r"""
        Return the relative dimension of the product of projective spaces.

        OUTPUT: a positive integer

        EXAMPLES::

            sage: T.<a,x,y,z,u,v,w> = ProductProjectiveSpaces([3, 2], QQ)
            sage: T.dimension_relative()
            5
        """
        return sum(self._dims)

    def dimension_absolute(self):
        r"""
        Return the absolute dimension of the product of projective spaces.

        OUTPUT: a positive integer

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([2, 2], GF(17))
            sage: T.dimension_absolute()
            4
            sage: T.dimension()
            4
        """
        base = self.base_scheme()
        if base.is_noetherian():
            return sum([self[i].dimension_relative() + base.dimension() for i in range(self.num_components())])
        raise NotImplementedError("cannot compute the dimension of this scheme.")

    dimension = dimension_absolute

    def dimension_relative_components(self):
        r"""
        Return the relative dimension of the product of projective spaces.

        OUTPUT: list of positive integers

        EXAMPLES::

            sage: T.<a,x,y,z,u,v,w> = ProductProjectiveSpaces([3, 2], QQ)
            sage: T.dimension_relative_components()
            [3, 2]
        """
        return self._dims

    def dimension_absolute_components(self):
        r"""
        Return the absolute dimension of the product of projective spaces.

        OUTPUT: list of positive integers

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([2, 2], GF(17))
            sage: T.dimension_absolute_components()
            [2, 2]
            sage: T.dimension_components()
            [2, 2]
        """
        base = self.base_scheme()
        if base.is_noetherian():
            return [self[i].dimension_relative() + base.dimension() for i in range(self.num_components())]
        raise NotImplementedError("cannot compute the dimension of this scheme.")

    dimension_components = dimension_absolute_components

    def num_components(self):
        r"""
        Return the number of components of this space.

        OUTPUT: integer

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1, 1], GF(5), 'x')
            sage: T.num_components()
            3
        """
        return len(self._components)

    def ngens(self):
        r"""
        Return the number of generators of this space.

        This is the number of variables in the coordinate ring of the
        projective space.

        OUTPUT: integer

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1, 1], GF(5), 'x')
            sage: T.ngens()
            6
        """
        return sum([P.ngens() for P in self._components])

    def _factors(self, v):
        r"""
        Return the tuple/list ``v`` split into the components of this space.

        INPUT:

        - ``v`` -- list or tuple

        OUTPUT: list of lists

        EXAMPLES::

            sage: T = ProductProjectiveSpaces([1, 1, 1], QQ, 'x')
            sage: T._factors([1, 2, 3, 4, 5, 6])
            [[1, 2], [3, 4], [5, 6]]
        """
        if not isinstance(v, (list, tuple, ETuple)):
            raise TypeError("%s, must be a list or tuple" % v)
        if len(v) != self.ngens():
            raise ValueError("%s must have %s elements" % (v, self.ngens()))
        index = 0
        splitv = []
        dims = self._dims
        for i in range(len(dims)):
            splitv.append(v[index:index+dims[i]+1])
            index += dims[i]+1
        return splitv

    def _degree(self, polynomial):
        r"""
        Return the homogeneous degrees.

        INPUT:

        - ``polynomial`` -- a polynomial in the coordinate_ring

        OUTPUT: a tuple of integers, one for each projective space component. A
        :exc:`ValueError` is raised if the polynomial is not multihomogeneous.

        EXAMPLES::

            sage: P1xP1.<x,y,s,t> = ProductProjectiveSpaces([1, 1], QQ)
            sage: P1xP1._degree(x^2*t + y^2*s)
            [2, 1]
            sage: P1xP1._degree(x + s)
            Traceback (most recent call last):
            ...
            ValueError: polynomial is not multi-homogeneous
        """
        E = polynomial.exponents()
        if not E:
            return []
        d = [sum(t) for t in self._factors(E[0])]
        for k in range(len(E)):
            if d != [sum(t) for t in self._factors(E[k])]:
                raise ValueError("polynomial is not multi-homogeneous")
        return d

    def _point(self, *args, **kwds):
        """
        Construct a point.

        For internal use only. See :mod:`morphism` for details.

        EXAMPLES::

            sage: P = ProductProjectiveSpaces([1, 1], QQ, 'z')
            sage: point_homset = P._point_homset(Spec(QQ), P)
            sage: P._point(point_homset, [2, 2, 1, 1])
            (1 : 1 , 1 : 1)
        """
        return ProductProjectiveSpaces_point_ring(*args, **kwds)

    def _morphism(self, *args, **kwds):
        """
        Construct a morphism.

        For internal use only. See :mod:`morphism` for details.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1, 1], QQ)
            sage: P._morphism(P.Hom(P), [x, y, z, w])
            Scheme endomorphism of Product of projective spaces P^1 x P^1 over Rational Field
              Defn: Defined by sending (x : y , z : w) to
                    (x : y , z : w).
        """
        return ProductProjectiveSpaces_morphism_ring(*args, **kwds)

    def _point_homset(self, *args, **kwds):
        """
        Construct a point Hom-set.

        For internal use only. See :mod:`morphism` for details.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1,1], ZZ)
            sage: P._point_homset(Spec(ZZ), P)
            Set of rational points of Product of projective spaces P^1 x P^1 over
            Integer Ring
        """
        return SchemeHomset_points_product_projective_spaces_ring(*args, **kwds)

    def _validate(self, polynomials):
        r"""
        If ``polynomials`` is a tuple of valid polynomial functions on this space,
        return ``polynomials``, otherwise raise a :exc:`TypeError`.

        Since this is a product of projective spaces, the polynomials must be multi-homogeneous.

        INPUT:

        - ``polynomials`` -- tuple of polynomials in the coordinate ring of
          this projective space

        OUTPUT: tuple of polynomials in the coordinate ring of this space

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: T._validate([x^2*u, y^2*w, z^2*u, w^2, u^2])
            [x^2*u, y^2*w, z^2*u, w^2, u^2]

        ::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: T._validate([x^2+w^2, y^2*w, z^2*u, w^2, u^2])
            Traceback (most recent call last):
            ...
            ValueError: polynomial is not multi-homogeneous

        ::

            sage: R.<t> = PolynomialRing(GF(5))
            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: T._validate([t, t, t, w^2, u^2])
            Traceback (most recent call last):
            ...
            TypeError: polynomials (=[t, t, t, w^2, u^2]) must be elements of Multivariate
            Polynomial Ring in x, y, z, w, u over Rational Field
        """
        if not isinstance(polynomials, (list, tuple)):
            raise TypeError('the argument polynomials=%s must be a list or tuple' % polynomials)
        #check in the coordinate ring
        source_ring = self.coordinate_ring()
        try:
            polynomials = [source_ring(poly) for poly in polynomials]
        except TypeError:
            raise TypeError("polynomials (=%s) must be elements of %s" % (polynomials,source_ring))
        for f in polynomials:
            self._degree(f)  # raises a ValueError if not multi-homogeneous
        return polynomials

    def _check_satisfies_equations(self, v):
        """
        Return ``True`` if ``v`` defines a point on the scheme this space;
        raise a :exc:`TypeError` otherwise.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: T._check_satisfies_equations([0, 1, 1, 1, 1])
            True

        ::

            sage: R.<t> = PolynomialRing(GF(7))
            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], R)
            sage: T._check_satisfies_equations([1 + t, 1, 0, 0, 1])
            True

        ::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], ZZ)
            sage: T._check_satisfies_equations([1, 1, 1, 0, 0])
            Traceback (most recent call last):
            ...
            TypeError: the zero vector is not a point in projective space

        ::

            sage: T.<x,y,z,w> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: T._check_satisfies_equations([1, 1, 1, 0, 0])
            Traceback (most recent call last):
            ...
            TypeError: the list v=[1, 1, 1, 0, 0] must have 4 components

        ::

            sage: T.<x,y,z,w> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: T._check_satisfies_equations([1, 1/2, 1, 0])
            Traceback (most recent call last):
            ...
            TypeError: the components of v=[1, 1/2, 1, 0] must be elements of Integer Ring
        """
        if not isinstance(v, (list, tuple)):
            raise TypeError('the argument v=%s must be a list or tuple' % v)
        n = self.ngens()
        if not len(v) == n:
            raise TypeError('the list v=%s must have %s components' % (v, n))
        R = self.base_ring()
        try:
            n = [R(w) for w in v]
        except TypeError:
            raise TypeError('the components of v=%s must be elements of %s' % (v, R))
        #check if any of the component points are 0
        N = self._dims
        start = 0
        for i in range(len(N)):
            if v[start:start + N[i]+1] == [R(0)]*(N[i]+1):
                raise TypeError('the zero vector is not a point in projective space')
            start += N[i]+1
        return True

    def _an_element_(self):
        r"""
        Return a (preferably typical) element of this space.

        This is used both for illustration and testing purposes.

        OUTPUT: a point in the this projective space

        EXAMPLES::

            sage: ProductProjectiveSpaces([1, 2, 3], ZZ).an_element()
            (7 : 1 , 7 : 6 : 1 , 7 : 6 : 5 : 1)
            sage: ProductProjectiveSpaces([3, 2, 1], PolynomialRing(ZZ, 'y')).an_element()
            (7*y : 6*y : 5*y : 1 , 7*y : 6*y : 1 , 7*y : 1)
        """
        v = [R.an_element() for R in self._components]
        return self(v)

    def subscheme(self, X):
        r"""
        Return the closed subscheme defined by ``X``.

        INPUT:

        - ``X`` -- list or tuple of equations

        OUTPUT: :class:`AlgebraicScheme_subscheme_projective_cartesian_product`

        EXAMPLES::

            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1, 1], GF(5))
            sage: X = P.subscheme([x - y, z - w]); X
            Closed subscheme of Product of projective spaces P^1 x P^1
             over Finite Field of size 5 defined by:
              x - y,
              z - w
            sage: X.defining_polynomials()
            [x - y, z - w]
            sage: I = X.defining_ideal(); I
            Ideal (x - y, z - w) of Multivariate Polynomial Ring in x, y, z, w
             over Finite Field of size 5
            sage: X.dimension()                                                         # needs sage.libs.singular
            0
            sage: X.base_ring()
            Finite Field of size 5
            sage: X.base_scheme()
            Spectrum of Finite Field of size 5
            sage: X.structure_morphism()
            Scheme morphism:
              From: Closed subscheme of Product of projective spaces P^1 x P^1
                    over Finite Field of size 5 defined by: x - y, z - w
              To:   Spectrum of Finite Field of size 5
              Defn: Structure map
        """
        return AlgebraicScheme_subscheme_product_projective(self, X)

    def change_ring(self, R):
        r"""
        Return a product of projective spaces over a ring ``R`` and otherwise
        the same as this projective space.

        INPUT:

        - ``R`` -- commutative ring or morphism

        OUTPUT: product of projective spaces over ``R``

        .. NOTE::

            There is no need to have any relation between ``R`` and the base ring
            of this space, if you want to have such a relation, use
            ``self.base_extend(R)`` instead.

        EXAMPLES::

            sage: T.<x,y,z,u,v,w> = ProductProjectiveSpaces([2, 2], QQ)
            sage: T.change_ring(GF(17))
            Product of projective spaces P^2 x P^2 over Finite Field of size 17
        """
        new_components = [P.change_ring(R) for P in self._components]
        return ProductProjectiveSpaces(new_components)

    def affine_patch(self, I, return_embedding=False):
        r"""
        Return the `I`-th affine patch of this projective space product
        where ``I`` is a multi-index.

        INPUT:

        - ``I`` -- list or tuple of positive integers

        - ``return_embedding`` -- boolean; if ``True`` the projective embedding
          is also returned

        OUTPUT:

        - An affine space.

        - An embedding into a product of projective spaces (optional).

        EXAMPLES::

            sage: PP = ProductProjectiveSpaces([2, 2, 2], ZZ, 'x')
            sage: phi = PP.affine_patch([0, 1, 2], True)
            sage: phi.domain()
            Affine Space of dimension 6 over Integer Ring
            sage: phi
            Scheme morphism:
              From: Affine Space of dimension 6 over Integer Ring
              To:   Product of projective spaces P^2 x P^2 x P^2 over Integer Ring
              Defn: Defined on coordinates by sending (x0, x1, x2, x3, x4, x5) to
                    (1 : x0 : x1 , x2 : 1 : x3 , x4 : x5 : 1)
        """
        if not isinstance(I, (list, tuple)):
            raise TypeError('the argument I=%s must be a list or tuple of positive integers' % I)
        PP = self.ambient_space()
        N = PP._dims
        if len(I) != len(N):
            raise ValueError('the argument I=%s must have %s entries' % (I,len(N)))
        I = tuple([int(i) for i in I])   # implicit type checking
        for i in range(len(I)):
            if I[i] < 0 or I[i] > N[i]:
                raise ValueError("argument i (= %s) must be between 0 and %s." % (I[i], N[i]))
        try:
            if return_embedding:
                return self.__affine_patches[I][1]
            else:
                return self.__affine_patches[I][0]
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        from sage.schemes.affine.affine_space import AffineSpace
        AA = AffineSpace(PP.base_ring(), sum(N), 'x')
        v = list(AA.gens())
        index = 0
        for i in range(len(I)):
            v.insert(index+I[i], 1)
            index += N[i]+1
        phi = AA.hom(v, self)
        self.__affine_patches.update({I:(AA, phi)})
        if return_embedding:
            return phi
        else:
            return AA

    @cached_method
    def segre_embedding(self, PP=None, var='u'):
        r"""
        Return the Segre embedding of this space into the appropriate
        projective space.

        INPUT:

        - ``PP`` -- (default: ``None``) ambient image projective space;
            this is constructed if it is not given

        - ``var`` -- string (default: ``'u'``); variable name of the image
          projective space

        OUTPUT: Hom; from this space to the appropriate subscheme of projective
        space

        .. TODO::

            Cartesian products with more than two components.

        EXAMPLES::

            sage: X.<y0,y1,y2,y3,y4,y5> = ProductProjectiveSpaces(ZZ, [2, 2])
            sage: phi = X.segre_embedding(); phi                                        # needs sage.libs.singular
            Scheme morphism:
              From: Product of projective spaces P^2 x P^2 over Integer Ring
              To:   Closed subscheme of Projective Space of dimension 8 over Integer Ring
                    defined by:
                      -u5*u7 + u4*u8,       -u5*u6 + u3*u8,       -u4*u6 + u3*u7,
                      -u2*u7 + u1*u8,       -u2*u4 + u1*u5,       -u2*u6 + u0*u8,
                      -u1*u6 + u0*u7,       -u2*u3 + u0*u5,       -u1*u3 + u0*u4
              Defn: Defined by sending (y0 : y1 : y2 , y3 : y4 : y5) to
                    (y0*y3 : y0*y4 : y0*y5 : y1*y3 : y1*y4 : y1*y5 : y2*y3 : y2*y4 : y2*y5).

            ::

            sage: T = ProductProjectiveSpaces([1, 2], CC, 'z')                          # needs sage.rings.real_mpfr
            sage: T.segre_embedding()                                                   # needs sage.libs.singular sage.rings.real_mpfr
            Scheme morphism:
              From: Product of projective spaces P^1 x P^2
                    over Complex Field with 53 bits of precision
              To:   Closed subscheme of Projective Space of dimension 5
                    over Complex Field with 53 bits of precision defined by:
                      -u2*u4 + u1*u5,       -u2*u3 + u0*u5,       -u1*u3 + u0*u4
              Defn: Defined by sending (z0 : z1 , z2 : z3 : z4) to
                    (z0*z2 : z0*z3 : z0*z4 : z1*z2 : z1*z3 : z1*z4).

            ::

            sage: T = ProductProjectiveSpaces([1, 2, 1], QQ, 'z')
            sage: T.segre_embedding()                                                   # needs sage.libs.singular
            Scheme morphism:
              From: Product of projective spaces P^1 x P^2 x P^1 over Rational Field
              To:   Closed subscheme of Projective Space of dimension 11
                    over Rational Field defined by:
                      -u9*u10 + u8*u11,     -u7*u10 + u6*u11,     -u7*u8 + u6*u9,
                      -u5*u10 + u4*u11,     -u5*u8 + u4*u9,       -u5*u6 + u4*u7,
                      -u5*u9 + u3*u11,      -u5*u8 + u3*u10,      -u5*u8 + u2*u11,
                      -u4*u8 + u2*u10,      -u3*u8 + u2*u9,       -u3*u6 + u2*u7,
                      -u3*u4 + u2*u5,       -u5*u7 + u1*u11,      -u5*u6 + u1*u10,
                      -u3*u7 + u1*u9,       -u3*u6 + u1*u8,       -u5*u6 + u0*u11,
                      -u4*u6 + u0*u10,      -u3*u6 + u0*u9,       -u2*u6 + u0*u8,
                      -u1*u6 + u0*u7,       -u1*u4 + u0*u5,       -u1*u2 + u0*u3
              Defn: Defined by sending (z0 : z1 , z2 : z3 : z4 , z5 : z6) to
                    (z0*z2*z5 : z0*z2*z6 : z0*z3*z5 : z0*z3*z6 : z0*z4*z5 : z0*z4*z6
                     : z1*z2*z5 : z1*z2*z6 : z1*z3*z5 : z1*z3*z6 : z1*z4*z5 : z1*z4*z6).
        """
        N = self._dims
        M = prod([n+1 for n in N]) - 1
        CR = self.coordinate_ring()

        vars = list(self.coordinate_ring().variable_names()) + [var + str(i) for i in range(M+1)]
        R = PolynomialRing(self.base_ring(), self.ngens()+M+1, vars, order='lex')

        #set-up the elimination for the segre embedding
        mapping = []
        k = self.ngens()
        index = self.num_components()*[0]
        for count in range(M + 1):
            mapping.append(R.gen(k+count)-prod([CR(self[i].gen(index[i])) for i in range(len(index))]))
            for i in range(len(index)-1, -1, -1):
                if index[i] == N[i]:
                    index[i] = 0
                else:
                    index[i] += 1
                    break #only increment once

        #change the defining ideal of the subscheme into the variables
        I = R.ideal(list(self.defining_polynomials()) + mapping)
        J = I.groebner_basis()
        s = set(R.gens()[:self.ngens()])
        n = len(J)-1
        L = []
        while s.isdisjoint(J[n].variables()):
            L.append(J[n])
            n = n-1

        #create new subscheme
        if PP is None:
            PS = ProjectiveSpace(self.base_ring(), M, R.variable_names()[self.ngens():])
            Y = PS.subscheme(L)
        else:
            if PP.dimension_relative() != M:
                raise ValueError("projective Space %s must be dimension %s") % (PP, M)
            S = PP.coordinate_ring()
            psi = R.hom([0]*k + list(S.gens()), S)
            L = [psi(l) for l in L]
            Y = PP.subscheme(L)

        #create embedding for points
        mapping = []
        index = self.num_components()*[0]
        for count in range(M + 1):
            mapping.append(prod([CR(self[i].gen(index[i])) for i in range(len(index))]))
            for i in range(len(index)-1, -1, -1):
                if index[i] == N[i]:
                    index[i] = 0
                else:
                    index[i] += 1
                    break #only increment once
        phi = self.hom(mapping, Y)

        return phi

class ProductProjectiveSpaces_field(ProductProjectiveSpaces_ring):
    def _point(self, *args, **kwds):
        """
        Construct a point.

        For internal use only. See :mod:`morphism` for details.

        EXAMPLES::

            sage: u = QQ['u'].0
            sage: K = NumberField(u^2 - 2, 'v')                                         # needs sage.rings.number_field
            sage: P = ProductProjectiveSpaces([1, 2], K, 'x')                           # needs sage.rings.number_field
            sage: P([1, 3, u, 1, 1])                                                    # needs sage.rings.number_field
            (1/3 : 1 , v : 1 : 1)
        """
        return ProductProjectiveSpaces_point_field(*args, **kwds)

    def _point_homset(self, *args, **kwds):
        """
        Construct a point Hom-set.

        For internal use only. See :mod:`morphism` for details.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1, 1], GF(5))
            sage: P._point_homset(Spec(GF(5)), P)
            Set of rational points of Product of projective spaces P^1 x P^1
             over Finite Field of size 5
        """
        return SchemeHomset_points_product_projective_spaces_field(*args, **kwds)

    def points_of_bounded_height(self, **kwds):
        r"""
        Return an iterator of the points in this product of projective spaces
        with the absolute heights of the components of at most the given bound.

        Bound check is strict for the rational field. Requires the base field
        of this space to be a number field. Uses the Doyle-Krumm algorithm 4
        (algorithm 5 for imaginary quadratic) for computing algebraic numbers
        up to a given height [DK2013]_.

        The algorithm requires floating point arithmetic, so the user is
        allowed to specify the precision for such calculations.
        Additionally, due to floating point issues, points
        slightly larger than the bound may be returned. This can be controlled
        by lowering the tolerance.


        INPUT:

        - ``bound`` -- a real number

        - ``tolerance`` -- a rational number in (0,1] used in Doyle-Krumm
          algorithm-4

        - ``precision`` -- the precision to use for computing the elements of
          bounded height of number fields

        OUTPUT: an iterator of points in this space

        EXAMPLES::

            sage: PP = ProductProjectiveSpaces(QQ, [1, 2])
            sage: sorted(list(PP.points_of_bounded_height(bound=1)))
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

            sage: u = QQ['u'].0
            sage: P = ProductProjectiveSpaces([1, 1], NumberField(u^2 - 2, 'v'))        # needs sage.rings.number_field
            sage: sorted(list(P.points_of_bounded_height(bound=1.5)))                   # needs sage.rings.number_field
            [(-v : 1 , -v : 1), (-v : 1 , -1 : 1), (-v : 1 , -1/2*v : 1), (-v : 1 , 0 : 1),
             (-v : 1 , 1/2*v : 1), (-v : 1 , 1 : 0), (-v : 1 , 1 : 1), (-v : 1 , v : 1),
             (-1 : 1 , -v : 1), (-1 : 1 , -1 : 1), (-1 : 1 , -1/2*v : 1), (-1 : 1 , 0 : 1),
             (-1 : 1 , 1/2*v : 1), (-1 : 1 , 1 : 0), (-1 : 1 , 1 : 1), (-1 : 1 , v : 1),
             (-1/2*v : 1 , -v : 1), (-1/2*v : 1 , -1 : 1), (-1/2*v : 1 , -1/2*v : 1),
             (-1/2*v : 1 , 0 : 1), (-1/2*v : 1 , 1/2*v : 1), (-1/2*v : 1 , 1 : 0),
             (-1/2*v : 1 , 1 : 1), (-1/2*v : 1 , v : 1), (0 : 1 , -v : 1), (0 : 1 , -1 : 1),
             (0 : 1 , -1/2*v : 1), (0 : 1 , 0 : 1), (0 : 1 , 1/2*v : 1), (0 : 1 , 1 : 0),
             (0 : 1 , 1 : 1), (0 : 1 , v : 1), (1/2*v : 1 , -v : 1), (1/2*v : 1 , -1 : 1),
             (1/2*v : 1 , -1/2*v : 1), (1/2*v : 1 , 0 : 1), (1/2*v : 1 , 1/2*v : 1),
             (1/2*v : 1 , 1 : 0), (1/2*v : 1 , 1 : 1), (1/2*v : 1 , v : 1), (1 : 0 , -v : 1),
             (1 : 0 , -1 : 1), (1 : 0 , -1/2*v : 1), (1 : 0 , 0 : 1), (1 : 0 , 1/2*v : 1),
             (1 : 0 , 1 : 0), (1 : 0 , 1 : 1), (1 : 0 , v : 1), (1 : 1 , -v : 1),
             (1 : 1 , -1 : 1), (1 : 1 , -1/2*v : 1), (1 : 1 , 0 : 1), (1 : 1 , 1/2*v : 1),
             (1 : 1 , 1 : 0), (1 : 1 , 1 : 1), (1 : 1 , v : 1), (v : 1 , -v : 1),
             (v : 1 , -1 : 1), (v : 1 , -1/2*v : 1), (v : 1 , 0 : 1), (v : 1 , 1/2*v : 1),
             (v : 1 , 1 : 0), (v : 1 , 1 : 1), (v : 1 , v : 1)]
        """
        B = kwds.pop('bound')
        tol = kwds.pop('tolerance', 1e-2)
        prec = kwds.pop('precision', 53)
        m = self.num_components()
        iters = [ self[i].points_of_bounded_height(bound=B, tolerance=tol, precision=prec) for i in range(m) ]
        dim = [self[i].dimension_relative() + 1 for i in range(m)]

        dim_prefix = [0, dim[0]] # prefixes dim list
        for i in range(1, len(dim)):
            dim_prefix.append(dim_prefix[i] + dim[i])

        P = []
        for i in range(m):
            pt = next(iters[i])
            P.extend(pt[j] for j in range(dim[i]))
            # initial value of P
        yield self(P)

        i = 0
        while i < m:
            try:
                pt = next(iters[i])
                for j in range(dim[i]):
                    P[dim_prefix[i] + j] = pt[j]
                yield self(P)
                i = 0
            except StopIteration:
                iters[i] = self[i].points_of_bounded_height(bound=B, tolerance=tol, precision=prec)
                pt = next(iters[i]) # reset
                for j in range(dim[i]):
                    P[dim_prefix[i] + j] = pt[j]
                i += 1

class ProductProjectiveSpaces_finite_field(ProductProjectiveSpaces_field):
    def _point(self, *args, **kwds):
        r"""
        Construct a point.

        For internal use only. See :mod:`morphism` for details.

        EXAMPLES::

            sage: P = ProductProjectiveSpaces([1, 2], GF(11))
            sage: P([3, 7, 4, 5, 9])
            (2 : 1 , 9 : 3 : 1)
        """
        return ProductProjectiveSpaces_point_finite_field(*args, **kwds)

    def __iter__(self):
        r"""
        Return an iterator over the elements of this product of projective spaces.

        EXAMPLES::

            sage: P = ProductProjectiveSpaces([2, 1], GF(3))
            sage: [x for x in P]
            [(0 : 0 : 1 , 0 : 1),
             (0 : 1 : 1 , 0 : 1),
             (0 : 2 : 1 , 0 : 1),
             ...
             (1 : 1 : 0 , 1 : 0),
             (2 : 1 : 0 , 1 : 0),
             (1 : 0 : 0 , 1 : 0)]
        """
        iters = [iter(T) for T in self._components]
        L = [next(x) for x in iters]  # put at zero
        yield self(L)
        j = 0
        while j < self.num_components():
            try:
                L[j] = next(iters[j])
                yield self(L)
                j = 0
            except StopIteration:
                iters[j] = iter(self[j])  # reset
                L[j] = next(iters[j])  # put at zero
                j += 1

    def rational_points(self, F=None):
        r"""
        Return the list of `F`-rational points on this product of projective spaces,
        where `F` is a given finite field, or the base ring of this space.

        EXAMPLES::

            sage: P = ProductProjectiveSpaces([1, 1], GF(5))
            sage: P.rational_points()
            [(0 : 1 , 0 : 1), (1 : 1 , 0 : 1), (2 : 1 , 0 : 1), (3 : 1 , 0 : 1), (4 : 1 , 0 : 1), (1 : 0 , 0 : 1),
             (0 : 1 , 1 : 1), (1 : 1 , 1 : 1), (2 : 1 , 1 : 1), (3 : 1 , 1 : 1), (4 : 1 , 1 : 1), (1 : 0 , 1 : 1),
             (0 : 1 , 2 : 1), (1 : 1 , 2 : 1), (2 : 1 , 2 : 1), (3 : 1 , 2 : 1), (4 : 1 , 2 : 1), (1 : 0 , 2 : 1),
             (0 : 1 , 3 : 1), (1 : 1 , 3 : 1), (2 : 1 , 3 : 1), (3 : 1 , 3 : 1), (4 : 1 , 3 : 1), (1 : 0 , 3 : 1),
             (0 : 1 , 4 : 1), (1 : 1 , 4 : 1), (2 : 1 , 4 : 1), (3 : 1 , 4 : 1), (4 : 1 , 4 : 1), (1 : 0 , 4 : 1),
             (0 : 1 , 1 : 0), (1 : 1 , 1 : 0), (2 : 1 , 1 : 0), (3 : 1 , 1 : 0), (4 : 1 , 1 : 0), (1 : 0 , 1 : 0)]

        ::

            sage: P = ProductProjectiveSpaces([1, 1], GF(2))
            sage: sorted(P.rational_points(GF(2^2, 'a')), key=str)                      # needs sage.rings.finite_rings
            [(0 : 1 , 0 : 1), (0 : 1 , 1 : 0), (0 : 1 , 1 : 1), (0 : 1 , a + 1 : 1), (0 : 1 , a : 1),
             (1 : 0 , 0 : 1), (1 : 0 , 1 : 0), (1 : 0 , 1 : 1), (1 : 0 , a + 1 : 1), (1 : 0 , a : 1),
             (1 : 1 , 0 : 1), (1 : 1 , 1 : 0), (1 : 1 , 1 : 1), (1 : 1 , a + 1 : 1), (1 : 1 , a : 1),
             (a + 1 : 1 , 0 : 1), (a + 1 : 1 , 1 : 0), (a + 1 : 1 , 1 : 1), (a + 1 : 1 , a + 1 : 1), (a + 1 : 1 , a : 1),
             (a : 1 , 0 : 1), (a : 1 , 1 : 0), (a : 1 , 1 : 1), (a : 1 , a + 1 : 1), (a : 1 , a : 1)]
        """
        if F is None:
            return list(self)
        elif not isinstance(F, FiniteField):
            raise TypeError("second argument (= %s) must be a finite field" % F)
        return list(self.base_extend(F))
