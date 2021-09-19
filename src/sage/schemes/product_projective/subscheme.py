r"""
Subschemes of products of projective spaces

AUTHORS:

- Ben Hutz (2014): subschemes of Cartesian products of projective space
"""

#*****************************************************************************
# Copyright (C) 2014 Ben Hutz <bn4941@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
# http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.misc_c import prod
from sage.misc.cachefunc import cached_method
from sage.rings.fraction_field import FractionField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.schemes.affine.affine_space import AffineSpace
from sage.schemes.projective.projective_subscheme import AlgebraicScheme_subscheme_projective
from sage.schemes.projective.projective_space import ProjectiveSpace

class AlgebraicScheme_subscheme_product_projective(AlgebraicScheme_subscheme_projective):
    r"""
    Construct an algebraic subscheme of a product of projective spaces.

    .. WARNING::

        You should not create objects of this class directly. The
        preferred method to construct such subschemes is to use
        :meth:`~sage.schemes.product_projective.subscheme`
        method of :class:`Product of Projective Spaces
        <sage.schemes.product_projective.space.ProductProjectiveSpaces_ring>`.

    INPUT:

    - ``A`` -- ambient :class:`Product of Projective Spaces
      <sage.schemes.product_projective.space.ProductProjectiveSpaces_ring>`.

    - ``polynomials`` -- single polynomial, ideal or iterable of
      defining multi-homogeneous polynomials.

    EXAMPLES::

        sage: P.<x, y, u, v> = ProductProjectiveSpaces([1,1], QQ)
        sage: P.subscheme([u*x^2-v*y*x])
        Closed subscheme of Product of projective spaces P^1 x P^1 over Rational
        Field defined by:
          x^2*u - x*y*v

    TESTS::

        sage: from sage.schemes.product_projective.subscheme \
              import AlgebraicScheme_subscheme_product_projective
        sage: AlgebraicScheme_subscheme_product_projective(P, [u*x^2-v*y*x])
        Closed subscheme of Product of projective spaces P^1 x P^1 over Rational
        Field defined by:
          x^2*u - x*y*v
    """

    @cached_method
    def segre_embedding(self, PP=None):
        r"""
        Return the Segre embedding of this subscheme into the appropriate projective
        space.

        INPUT:

        - ``PP`` -- (default: ``None``) ambient image projective space;
          this is constructed if it is not given.

        OUTPUT:

        Hom from this subscheme to the appropriate subscheme of projective space

        EXAMPLES::

            sage: X.<x,y,z,w,u,v> = ProductProjectiveSpaces([2,2], QQ)
            sage: P = ProjectiveSpace(QQ,8,'t')
            sage: L = (-w - v)*x + (-w*y - u*z)
            sage: Q = (-u*w - v^2)*x^2 + ((-w^2 - u*w + (-u*v - u^2))*y + (-w^2 - u*v)*z)*x + \
            ((-w^2 - u*w - u^2)*y^2 + (-u*w - v^2)*z*y + (-w^2 + (-v - u)*w)*z^2)
            sage: W = X.subscheme([L,Q])
            sage: phi = W.segre_embedding(P)
            sage: phi.codomain().ambient_space() == P
            True

        ::

            sage: PP.<x,y,u,v,s,t> = ProductProjectiveSpaces([1,1,1], CC)
            sage: PP.subscheme([]).segre_embedding()
            Scheme morphism:
              From: Closed subscheme of Product of projective spaces P^1 x P^1 x P^1
            over Complex Field with 53 bits of precision defined by:
              (no polynomials)
              To:   Closed subscheme of Projective Space of dimension 7 over Complex
            Field with 53 bits of precision defined by:
              -u5*u6 + u4*u7,
              -u3*u6 + u2*u7,
              -u3*u4 + u2*u5,
              -u3*u5 + u1*u7,
              -u3*u4 + u1*u6,
              -u3*u4 + u0*u7,
              -u2*u4 + u0*u6,
              -u1*u4 + u0*u5,
              -u1*u2 + u0*u3
              Defn: Defined by sending (x : y , u : v , s : t) to
                    (x*u*s : x*u*t : x*v*s : x*v*t : y*u*s : y*u*t : y*v*s : y*v*t).

        ::

            sage: PP.<x,y,z,u,v,s,t> = ProductProjectiveSpaces([2,1,1], ZZ)
            sage: PP.subscheme([x^3, u-v, s^2-t^2]).segre_embedding()
            Scheme morphism:
              From: Closed subscheme of Product of projective spaces P^2 x P^1 x P^1
            over Integer Ring defined by:
              x^3,
              u - v,
              s^2 - t^2
              To:   Closed subscheme of Projective Space of dimension 11 over
            Integer Ring defined by:
              u10^2 - u11^2,
              u9 - u11,
              u8 - u10,
              -u7*u10 + u6*u11,
              u6*u10 - u7*u11,
              u6^2 - u7^2,
              u5 - u7,
              u4 - u6,
              u3^3,
              -u3*u10 + u2*u11,
              u2*u10 - u3*u11,
              -u3*u6 + u2*u7,
              u2*u6 - u3*u7,
              u2*u3^2,
              u2^2 - u3^2,
              u1 - u3,
              u0 - u2
              Defn: Defined by sending (x : y : z , u : v , s : t) to
                    (x*u*s : x*u*t : x*v*s : x*v*t : y*u*s : y*u*t : y*v*s : y*v*t :
            z*u*s : z*u*t : z*v*s : z*v*t).
        """
        AS = self.ambient_space()
        CR = AS.coordinate_ring()
        N = AS.dimension_relative_components()
        M = prod([n+1 for n in N]) - 1

        vars = list(AS.coordinate_ring().variable_names()) + ['u' + str(i) for i in range(M+1)]
        R = PolynomialRing(AS.base_ring(), AS.ngens()+M+1, vars, order='lex')

        #set-up the elimination for the segre embedding
        mapping = []
        k = AS.ngens()
        index = AS.num_components()*[0]
        for count in range(M + 1):
            mapping.append(R.gen(k+count)-prod([CR(AS[i].gen(index[i])) for i in range(len(index))]))
            for i in range(len(index)-1, -1, -1):
                if index[i] == N[i]:
                    index[i] = 0
                else:
                    index[i] += 1
                    break #only increment once

        #change the defining ideal of the subscheme into the variables
        I = R.ideal(list(self.defining_polynomials()) + mapping)
        J = I.groebner_basis()
        s = set(R.gens()[:AS.ngens()])
        n = len(J)-1
        L = []
        while s.isdisjoint(J[n].variables()):
            L.append(J[n])
            n = n-1

        #create new subscheme
        if PP is None:
            PS = ProjectiveSpace(self.base_ring(), M, R.gens()[AS.ngens():])
            Y = PS.subscheme(L)
        else:
            if PP.dimension_relative() != M:
                raise ValueError("projective space %s must be dimension %s")%(PP, M)
            S = PP.coordinate_ring()
            psi = R.hom([0]*k + list(S.gens()), S)
            L = [psi(l) for l in L]
            Y = PP.subscheme(L)

        #create embedding for points
        mapping = []
        index = AS.num_components()*[0]
        for count in range(M + 1):
            mapping.append(prod([CR(AS[i].gen(index[i])) for i in range(len(index))]))
            for i in range(len(index)-1, -1, -1):
                if index[i] == N[i]:
                    index[i] = 0
                else:
                    index[i] += 1
                    break #only increment once
        phi = self.hom(mapping, Y)

        return phi

    def dimension(self):
        r"""
        Return the dimension of the algebraic subscheme.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: X.<x,y,z,w,u,v> = ProductProjectiveSpaces([2,2],QQ)
            sage: L = (-w - v)*x + (-w*y - u*z)
            sage: Q = (-u*w - v^2)*x^2 + ((-w^2 - u*w + (-u*v - u^2))*y + (-w^2 - u*v)*z)*x + \
            ((-w^2 - u*w - u^2)*y^2 + (-u*w - v^2)*z*y + (-w^2 + (-v - u)*w)*z^2)
            sage: W = X.subscheme([L,Q])
            sage: W.dimension()
            2

        ::

            sage: PP.<x,y,z,u,v,s,t> = ProductProjectiveSpaces([2,1,1], QQ)
            sage: X = PP.subscheme([x^3, x^5+y^5, z^6, x*u-v*y, s^2-t^2])
            sage: X.dimension()
            -1

        ::

            sage: PP = ProductProjectiveSpaces([2,1,3], CC, 't')
            sage: PP.subscheme([]).dimension()
            6

        ::

            sage: PP = ProductProjectiveSpaces([1,3,1], ZZ, 't')
            sage: PP.subscheme([]).dimension()
            5

        ::

            sage: PP.<x,y,u,v,s,t> = ProductProjectiveSpaces([1,1,1], CC)
            sage: X = PP.subscheme([x^2-y^2, u-v, s^2-t^2])
            sage: X.dimension()
            0
        """
        try:
            return self.__dimension
        except AttributeError:
            try:
                #move to field to compute radical
                X = self.change_ring(FractionField(self.base_ring()))
                PP = X.ambient_space()
                I = X.defining_ideal().radical()
                #check if the irrelevant ideal of any component is in the radical
                if any(all(t in I for t in PS.gens())
                       for PS in PP.components()):
                    self.__dimension = -1
                else:
                    self.__dimension = I.dimension() - PP.num_components()
            except TypeError:  #cannot compute radical for this base ring
                phi = self.segre_embedding()
                self.__dimension = phi.codomain().defining_ideal().dimension() - 1
            return self.__dimension

    def is_smooth(self, point=None):
        r"""
        Test whether the algebraic subscheme is smooth.

        EXAMPLES::

            sage: X.<x,y,z,w,u,v> = ProductProjectiveSpaces([2,2],QQ)
            sage: L = (-w - v)*x + (-w*y - u*z)
            sage: Q = (-u*w - v^2)*x^2 + ((-w^2 - u*w + (-u*v - u^2))*y + (-w^2 - u*v)*z)*x + \
            ((-w^2 - u*w - u^2)*y^2 + (-u*w - v^2)*z*y + (-w^2 + (-v - u)*w)*z^2)
            sage: W = X.subscheme([L,Q])
            sage: W.is_smooth()
            Traceback (most recent call last):
            ...
            NotImplementedError: Not Implemented
        """
        raise NotImplementedError("Not Implemented")

    def affine_patch(self, I, return_embedding = False):
        r"""
        Return the `I^{th}` affine patch of this projective scheme
        where 'I' is a multi-index.

        INPUT:

        - ``I`` -- a list or tuple of positive integers

        - ``return_embedding`` -- Boolean, if true the projective embedding is also returned

        OUTPUT:

        - An affine algebraic scheme

        - An embedding into a product of projective space (optional)

        EXAMPLES::

            sage: PP.<x,y,z,w,u,v> = ProductProjectiveSpaces([3,1],QQ)
            sage: W = PP.subscheme([y^2*z-x^3,z^2-w^2,u^3-v^3])
            sage: W.affine_patch([0,1],True)
            (Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
              x0^2*x1 - 1,
              x1^2 - x2^2,
              x3^3 - 1, Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 4 over Rational Field defined by:
              x0^2*x1 - 1,
              x1^2 - x2^2,
              x3^3 - 1
              To:   Closed subscheme of Product of projective spaces P^3 x P^1 over Rational Field defined by:
              -x^3 + y^2*z,
              z^2 - w^2,
              u^3 - v^3
              Defn: Defined on coordinates by sending (x0, x1, x2, x3) to
                    (1 : x0 : x1 : x2 , x3 : 1))
        """
        if not isinstance(I, (list, tuple)):
            raise TypeError('The argument I=%s must be a list or tuple of positive integers' % I)
        PP = self.ambient_space()
        N = PP.dimension_relative_components()
        if len(I) != len(N):
            raise ValueError('The argument I=%s must have %s entries'%(I,len(N)))
        I = tuple([int(i) for i in I])   # implicit type checking
        for i in range(len(I)):
            if I[i] < 0 or I[i] > N[i]:
                raise ValueError("Argument i (= %s) must be between 0 and %s."%(I[i], N[i]))
        #see if we've already created this affine patch
        try:
            if return_embedding:
                return self.__affine_patches[I]
            else:
                return self.__affine_patches[I][0]
        except AttributeError:
            self.__affine_patches = {}
        except KeyError:
            pass
        AA = AffineSpace(PP.base_ring(),sum(N),'x')
        v = list(AA.gens())
        # create the projective embedding
        index = 0
        for i in range(len(I)):
            v.insert(index+I[i],1)
            index += N[i]+1
        phi = AA.hom(v,self)
        #find the image of the subscheme
        polys = self.defining_polynomials()
        xi = phi.defining_polynomials()
        U = AA.subscheme([ f(xi) for f in polys ])
        phi = U.hom(v,self)
        self.__affine_patches.update({I:(U,phi)})
        if return_embedding:
            return U,phi
        else:
            return U

    def intersection_multiplicity(self, X, P):
        r"""
        Return the intersection multiplicity of this subscheme and the subscheme ``X`` at the point ``P``.

        This uses the intersection_multiplicity function for affine subschemes on affine patches of this subscheme
        and ``X`` that contain ``P``.

        INPUT:

        - ``X`` -- subscheme in the same ambient space as this subscheme.

        - ``P`` -- a point in the intersection of this subscheme with ``X``.

        OUTPUT: An integer.

        EXAMPLES:

        Multiplicity of a fixed point of the map `z^2 + \frac{1}{4}`::

            sage: PP.<x,y,u,v> = ProductProjectiveSpaces(QQ, [1,1])
            sage: G = PP.subscheme([(x^2 + 1/4*y^2)*v - y^2*u])
            sage: D = PP.subscheme([x*v - y*u])
            sage: sorted(G.intersection(D).rational_points())
            [(1/2 : 1 , 1/2 : 1), (1 : 0 , 1 : 0)]
            sage: Q = PP([1/2,1,1/2,1])
            sage: G.intersection_multiplicity(D, Q)
            2

        ::

            sage: F.<a> = GF(4)
            sage: PP.<x,y,z,u,v,w> = ProductProjectiveSpaces(F, [2,2])
            sage: X = PP.subscheme([z^5 + 3*x*y^4 + 8*y^5, u^2 - v^2])
            sage: Y = PP.subscheme([x^6 + z^6, w*z - v*y])
            sage: Q = PP([a,a+1,1,a,a,1])
            sage: X.intersection_multiplicity(Y, Q)
            16

        ::

            sage: PP.<x,y,z,u,v,w> = ProductProjectiveSpaces(QQ, [2,2])
            sage: X = PP.subscheme([x^2*u^3 + y*z*u*v^2, x - y])
            sage: Y = PP.subscheme([u^3 - w^3, x*v - y*w, z^3*w^2 - y^3*u*v])
            sage: Q = PP([0,0,1,0,1,0])
            sage: X.intersection_multiplicity(Y, Q)
            Traceback (most recent call last):
            ...
            TypeError: the intersection of this subscheme and (=Closed subscheme of Affine Space of dimension 4
            over Rational Field defined by: x2^3 - x3^3, -x1*x3 + x0, -x1^3*x2 + x3^2) must be proper and finite
        """
        PP = self.ambient_space()
        try:
            PP(P)
        except TypeError:
            raise TypeError("(=%s) must be a point in the ambient space of this subscheme and (=%s)"%(P,X))
        # find an affine chart of the ambient space of this subscheme that contains P
        indices = []
        aff_pt = []
        for i in range(PP.num_components()):
            Q = P[i]
            j = 0
            while Q[j] == 0:
                j = j + 1
            indices.append(j)
            T = list(Q)
            t = T.pop(j)
            aff_pt.extend([1/t*T[k] for k in range(PP.components()[i].dimension_relative())])
        X1 = self.affine_patch(indices)
        X2 = X.affine_patch(indices)
        return X1.intersection_multiplicity(X2, X1.ambient_space()(aff_pt))

    def multiplicity(self, P):
        r"""
        Return the multiplicity of ``P`` on this subscheme.

        This is computed as the multiplicity of the corresponding point on an affine patch of this subscheme
        that contains ``P``. This subscheme must be defined over a field. An error is returned if ``P``
        not a point on this subscheme.

        INPUT:

        - ``P`` -- a point on this subscheme.

        OUTPUT: an integer.

        EXAMPLES::

            sage: PP.<x,y,z,w> = ProductProjectiveSpaces(QQ, [1,1])
            sage: X = PP.subscheme([x^4*z^3 - y^4*w^3])
            sage: Q1 = PP([1,1,1,1])
            sage: X.multiplicity(Q1)
            1
            sage: Q2 = PP([0,1,1,0])
            sage: X.multiplicity(Q2)
            3

        ::

            sage: PP.<x,y,z,w,u> = ProductProjectiveSpaces(GF(11), [1,2])
            sage: X = PP.subscheme([x^7*u - y^7*z, u^6*x^2 - w^3*z^3*x*y - w^6*y^2])
            sage: Q1 = PP([1,0,10,1,0])
            sage: X.multiplicity(Q1)
            1
            sage: Q2 = PP([1,0,1,0,0])
            sage: X.multiplicity(Q2)
            4
        """
        PP = self.ambient_space()
        try:
            PP(P)
        except TypeError:
            raise TypeError("(={}) must be a point in the ambient space of this "
                            "subscheme and (={})".format(P, self))
        # find an affine chart of the ambient space of this subscheme that contains P
        indices = []
        aff_pt = []
        for i in range(PP.num_components()):
            Q = P[i]
            j = 0
            while Q[j] == 0:
                j = j + 1
            indices.append(j)
            T = list(Q)
            t = T.pop(j)
            aff_pt.extend([1/t*T[k] for k in range(PP.components()[i].dimension_relative())])
        X = self.affine_patch(indices)
        return X.multiplicity(X.ambient_space()(aff_pt))
