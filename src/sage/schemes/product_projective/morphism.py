r"""
Polynomial morphisms for products of projective spaces

This class builds on the projective space class and its point and morphism classes.

EXAMPLES::

    sage: P1xP1.<x,y,u,v> = ProductProjectiveSpaces(QQ, [1, 1])
    sage: H = End(P1xP1)
    sage: H([x^2*u, y^2*v, x*v^2, y*u^2])
    Scheme endomorphism of Product of projective spaces P^1 x P^1 over Rational Field
      Defn: Defined by sending (x : y , u : v) to
            (x^2*u : y^2*v , x*v^2 : y*u^2).
"""
# ****************************************************************************
# Copyright (C) 2014 Ben Hutz <bn4941@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
# https://www.gnu.org/licenses/
# ****************************************************************************
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.categories.fields import Fields
from sage.categories.number_fields import NumberFields
from sage.rings.number_field.order import is_NumberFieldOrder
from sage.rings.fraction_field import FractionField
from sage.rings.qqbar import QQbar
_Fields = Fields()


class ProductProjectiveSpaces_morphism_ring(SchemeMorphism_polynomial):
    r"""
    The class of morphisms on products of projective spaces.

    The components are projective space morphisms.

    EXAMPLES::

        sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
        sage: H = T.Hom(T)
        sage: H([x^2, y^2, z^2, w^2, u^2])
        Scheme endomorphism of Product of projective spaces P^2 x P^1 over Rational Field
          Defn: Defined by sending (x : y : z , w : u) to
                (x^2 : y^2 : z^2 , w^2 : u^2).
    """

    def __init__(self, parent, polys, check=True):
        r"""
        The Python constructor.

        INPUT:

        - ``parent`` -- Hom-set.

        - ``polys`` -- anything that defines a point in the class.

        - ``check`` -- Boolean. Whether or not to perform input checks.
          (Default:`` True``)

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: H = T.Hom(T)
            sage: H([x^2*u, y^2*w, z^2*u, w^2, u^2])
            Scheme endomorphism of Product of projective spaces P^2 x P^1 over Rational Field
              Defn: Defined by sending (x : y : z , w : u) to
                    (x^2*u : y^2*w : z^2*u , w^2 : u^2).

        ::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: H = T.Hom(T)
            sage: H([x^2*u, y^2*w, z^2*u, w^2, u*z])
            Traceback (most recent call last):
            ...
            TypeError: polys (=[x^2*u, y^2*w, z^2*u, w^2, z*u]) must be
            multi-homogeneous of the same degrees (by component)

        ::

            sage: R.<s,t> = PolynomialRing(QQ)
            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], QQ)
            sage: P.<u,v,w,s,t,r> = ProductProjectiveSpaces([3, 1], QQ)
            sage: H = Hom(Z,P)
            sage: f = H([a^2,b^2,a^2,a*b,a*x,b*z]); f
            Scheme morphism:
              From: Product of projective spaces P^1 x P^2 over Rational Field
              To:   Product of projective spaces P^3 x P^1 over Rational Field
              Defn: Defined by sending (a : b , x : y : z) to
                    (a^2 : b^2 : a^2 : a*b , a*x : b*z).

        ::

            sage: Z.<a,b,c,x,y,z> = ProductProjectiveSpaces([1, 3], QQ)
            sage: P.<u,v,w,s,t,r> = ProductProjectiveSpaces([2, 2], QQ)
            sage: H = Hom(Z,P)
            sage: f = H([a^2,b^2,c^2,x^2,y^2,z^2])
            Traceback (most recent call last):
            ...
            TypeError: polys (=[a^2, b^2, c^2, x^2, y^2, z^2]) must be
            multi-homogeneous of the same degrees (by component)
        """
        if check:
            #check multi-homogeneous
            #if self is a subscheme, we may need the lift of the polynomials
            try:
                polys[0].exponents()
            except AttributeError:
                polys = [f.lift() for f in polys]

            target = parent.codomain().ambient_space()
            dom = parent.domain().ambient_space()
            from sage.schemes.product_projective.space import is_ProductProjectiveSpaces
            if is_ProductProjectiveSpaces(target):
                splitpolys = target._factors(polys)
                for m in range(len(splitpolys)):
                    d = dom._degree(splitpolys[m][0])
                    if not all(d == dom._degree(f) for f in splitpolys[m]):
                        raise  TypeError("polys (=%s) must be multi-homogeneous of the same degrees (by component)"%polys)
            else:
                #we are mapping into some other kind of space
                target._validate(polys)

        SchemeMorphism_polynomial.__init__(self, parent, polys, check)

    def __getitem__(self, i):
        r"""
        Return the ``i``-th coordinate polynomial.

        INPUT:

        - ``i`` -- integer.

        OUTPUT:

        The (multi)-homogeneous polynomial that is the ``i``-th coordinate.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: H = T.Hom(T)
            sage: F = H([x^2*u, y^2*w, z^2*u, w^2, u^2])
            sage: F[2]
            z^2*u
        """
        return self._polys[i]

    def _repr_defn(self):
        r"""
        Return a string representation of this morphism.

        OUTPUT: String.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProductProjectiveSpaces([1, 1], QQ)
            sage: H = Hom(P,P)
            sage: f = H([x^2, y^2, z, w])
            sage: f._repr_defn()
            'Defined by sending (x : y , z : w) to \n(x^2 : y^2 , z : w).'
        """
        s  = 'Defined by sending '
        s += self.domain().ambient_space()._repr_generic_point()
        s += ' to \n'
        s += self.codomain().ambient_space()._repr_generic_point(self._polys)
        s += '.'
        return s

    def __call__(self, P, check=True):
        r"""
        Make morphisms of products of projective spaces callable.

        INPUT:

        - ``P`` -- a point in the domain.

        - ``check`` -- Boolean - whether or not to perform the input checks
          on the image point (Default: ``True``).

        OUTPUT: The image point in the codomain.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: H = T.Hom(T)
            sage: F = H([x^2*u, y^2*w, z^2*u, w^2, u^2])
            sage: F(T([2, 1, 3, 0, 1]))
            (4/9 : 0 : 1 , 0 : 1)

        ::

            sage: PP.<x,y,z,u,v,w> = ProductProjectiveSpaces(QQ, [1, 1, 1])
            sage: HP = End(PP)
            sage: f = HP([v*x^2,w*y^2,z^2,u^2,v^2,w^2])
            sage: Q = PP([0,1,1,1,1,1])
            sage: f(Q)
            (0 : 1 , 1 : 1 , 1 : 1)

        ::

            sage: PP.<t0,t1,t2,t3,t4> = ProductProjectiveSpaces([2,1], ZZ)
            sage: Q = PP([1,1,1,2,1])
            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1,2], ZZ)
            sage: H = End(Z)
            sage: f = H([a^3, b^3+a*b^2, x^2, y^2-z^2, z*y])
            sage: f(Q)
            Traceback (most recent call last):
            ...
            TypeError: (1 : 1 : 1 , 2 : 1) fails to convert into the map's domain
            Product of projective spaces P^1 x P^2 over Integer Ring, but a
            `pushforward` method is not properly implemented
            sage: f([1,1,1,2,1])
            (1 : 2 , 1 : 3 : 2)

        ::

            sage: PP.<x,y,u,v> = ProductProjectiveSpaces(ZZ, [1, 1])
            sage: HP = End(PP)
            sage: g = HP([x^2, y^2, u^2, v^2])
            sage: g([0, 0, 0, 0],check=False)
            (0 : 0 , 0 : 0)
        """
        from sage.schemes.product_projective.point import ProductProjectiveSpaces_point_ring
        if check:
            if not isinstance(P, ProductProjectiveSpaces_point_ring):
                try:
                    P = self.domain()(P)
                except (TypeError, NotImplementedError):
                    raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(P, self.domain()))
            elif self.domain()!= P.codomain():
                raise TypeError("%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(P, self.domain()))

        A = self.codomain()
        Q = list(P)
        newP = [f(Q) for f in self.defining_polynomials()]
        return A.point(newP, check)

    def __eq__(self, right):
        """
        Tests the equality of two product projective morphisms.

        INPUT:

        - ``right`` - a map on product of projective space.

        OUTPUT:

        - Boolean - True if ``self`` and ``right`` define the same product projective
          map. False otherwise.

        EXAMPLES::

            sage: P1.<x1,x2,x3,x4> = ProductProjectiveSpaces([1, 1], QQ)
            sage: P2.<y1,y2,y3,y4> = ProductProjectiveSpaces([1, 1], CC)
            sage: H1 = End(P1); H2 = End(P2)
            sage: f = H1([x1*x2, x2^2, x3*x4, x4^2])
            sage: g = H2([y1*y2, y2^2, y3*y4, y4^2])
            sage: f == g
            False

        ::

            sage: P.<x,y,u,v> = ProductProjectiveSpaces([1, 1], QQ)
            sage: H = Hom(P, P)
            sage: f = H([x^2, y^2, u^2, v^2])
            sage: g = H([x^2, x*y, u*v, u^2])
            sage: f == g
            False

        ::

            sage: PP.<x,y,z,a,b> = ProductProjectiveSpaces([2,1], ZZ)
            sage: H = End(PP)
            sage: f = H([x^2*y*z, x*y^2*z, x*y*z^2, a^2, b^2])
            sage: g = H([x, y, z, a^3, a*b^2])
            sage: f == g
            True
        """
        if not isinstance(right, SchemeMorphism_polynomial):
            return False
        if self.parent() != right.parent():
            return False
        PP = self.parent().codomain()

        n = PP.num_components()
        dim = [ P.ngens() for P in PP ]
        dim_prefix = [0,dim[0]]

        for i in range(1,n):
            dim_prefix.append(dim_prefix[i] + dim[i])

        # compare ratio of coordinates for each projective component
        for m in range(n):
            l = dim_prefix[m]
            r = dim_prefix[m] + dim[m]
            for i in range(l,r):
                for j in range(i+1,r):
                    if self[i]*right[j] != self[j]*right[i]:
                        return False
        return True

    def __ne__(self, right):
        """
        Tests the inequality of two prduct projective morphisms.

        INPUT:

        - ``right`` -- a map on product of projective space.

        OUTPUT:

        - Boolean -- True if ``self`` and ``right`` define different product
          projective maps. False otherwise.

        EXAMPLES::

            sage: PP.<a,b,x,y,z> = ProductProjectiveSpaces([1,2], ZZ)
            sage: E = End(PP)
            sage: f = E([a^3, a*b^2, x*y, y*z, z*x])
            sage: g = E([a*b, a^2, x^2, y^2, z^2])
            sage: f != g
            True
            sage: f != f
            False
        """
        if not isinstance(right, SchemeMorphism_polynomial):
            return True
        if self.parent() != right.parent():
            return True
        PP = self.parent().codomain()

        n = PP.num_components()
        dim = [ P.ngens() for P in PP ]
        dim_prefix = [0,dim[0]]

        for i in range(1,n):
            dim_prefix.append(dim_prefix[i] + dim[i])

        for m in range(n):
            l = dim_prefix[m]
            r = dim_prefix[m] + dim[m]
            for i in range(l,r):
                for j in range(i+1,r):
                    if self[i]*right[j] != self[j]*right[i]:
                        return True
        return False

    def is_morphism(self):
        r"""
        Returns ``True`` if this mapping is a morphism of products of projective spaces.

        For each component space of the codomain of this mapping we consider the subscheme of
        the domain of this map generated by the corresponding coordinates of the map.
        This map is a morphism if and only if each of these subschemes has no points.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1, 2], ZZ)
            sage: H = End(Z)
            sage: f = H([a^2, b^2, x*z-y*z, x^2-y^2, z^2])
            sage: f.is_morphism()
            False

        ::

            sage: P.<x,y,z,u,v,w>=ProductProjectiveSpaces([2, 2], QQ)
            sage: H = End(P)
            sage: f = H([u, v, w, u^2, v^2, w^2])
            sage: f.is_morphism()
            True

        ::

            sage: P.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: Q.<a,b,c,d,e> = ProductProjectiveSpaces([1, 2], QQ)
            sage: H = Hom(P, Q)
            sage: f = H([x^2, y^2, u^3, w^3, u^3])
            sage: f.is_morphism()
            False
        """
        m = 0
        T = self.domain().ambient_space()
        S = self.codomain().ambient_space()

        if T.base_ring().is_field():
            f = self
        else:
            f = self.change_ring(T.base_ring().fraction_field())
            T = T.change_ring(T.base_ring().fraction_field())

        for i in range(S.num_components()):
            t = S[i].dimension_relative() + 1
            X = T.subscheme(list(f)[m : m + t])
            if X.dimension() > -1:
                return False
            m = m + t

        return True

    def as_dynamical_system(self):
        """
        Return this endomorphism as a :class:`DynamicalSystem_producte_projective`.

        OUTPUT:

        - :class:`DynamicalSystem_produce_projective`

        EXAMPLES::

            sage: Z.<a,b,x,y,z> = ProductProjectiveSpaces([1 , 2], ZZ)
            sage: H = End(Z)
            sage: f = H([a^3, b^3, x^2, y^2, z^2])
            sage: type(f.as_dynamical_system())
            <class 'sage.dynamics.arithmetic_dynamics.product_projective_ds.DynamicalSystem_product_projective'>
        """
        if not self.is_endomorphism():
            raise TypeError("must be an endomorphism")
        from sage.dynamics.arithmetic_dynamics.product_projective_ds import DynamicalSystem_product_projective
        return DynamicalSystem_product_projective(list(self), self.domain())

    def global_height(self, prec=None):
        r"""
        Returns the maximum of the absolute logarithmic heights of the coefficients
        in any of the coordinate functions of this map.

        INPUT:

        - ``prec`` -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        - a real number.

        .. TODO::

            Add functionality for `\QQbar`, implement function to convert
            the map defined over `\QQbar` to map over a number field.

        EXAMPLES::

            sage: P1xP1.<x,y,u,v> = ProductProjectiveSpaces([1, 1], ZZ)
            sage: H = End(P1xP1)
            sage: f = H([x^2*u, 3*y^2*v, 5*x*v^2, y*u^2])
            sage: f.global_height()
            1.60943791243410

        ::

            sage: u = QQ['u'].0
            sage: R = NumberField(u^2 - 2, 'v')
            sage: PP.<x,y,a,b> = ProductProjectiveSpaces([1, 1], R)
            sage: H = End(PP)
            sage: O = R.maximal_order()
            sage: g = H([3*O(u)*x^2, 13*x*y, 7*a*y, 5*b*x + O(u)*a*y])
            sage: g.global_height()
            2.56494935746154
        """
        K = self.domain().base_ring()
        if K in NumberFields() or is_NumberFieldOrder(K):
            H = 0
            for i in range(self.domain().ambient_space().ngens()):
                C = self[i].coefficients()
                h = max(c.global_height(prec=prec) for c in C)
                H = max(H, h)
            return H
        elif K == QQbar:
            raise NotImplementedError("not implemented for QQbar")
        else:
            raise TypeError("Must be over a Numberfield or a Numberfield Order or QQbar")

    def local_height(self, v, prec=None):
        r"""
        Returns the maximum of the local height of the coefficients in any
        of the coordinate functions of this map.

        INPUT:

        - ``v`` -- a prime or prime ideal of the base ring.

        - ``prec`` -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        - a real number.

        EXAMPLES::

            sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
            sage: H = T.Hom(T)
            sage: f = H([4*x^2+3/100*y^2, 8/210*x*y, 1/10000*z^2, 20*w^2, 1/384*u*w])
            sage: f.local_height(2)
            4.85203026391962

        ::

            sage: R.<z> = PolynomialRing(QQ)
            sage: K.<w> = NumberField(z^2-5)
            sage: P.<x,y,a,b> = ProductProjectiveSpaces([1, 1], K)
            sage: H = Hom(P,P)
            sage: f = H([2*x^2 + w/3*y^2, 1/w*y^2, a^2, 6*b^2 + 1/9*a*b])
            sage: f.local_height(K.ideal(3))
            2.19722457733622
        """
        K = FractionField(self.domain().base_ring())
        if K not in NumberFields():
            raise TypeError("must be over a number field or a number field order")

        H = 0
        for i in range(self.domain().ambient_space().ngens()):
            C = self[i].coefficients()
            h = max(K(c).local_height(v, prec) for c in C)
            H = max(H, h)
        return H
