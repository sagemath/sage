"""
Morphisms of free modules

AUTHORS:
    - William Stein: initial version

    - Miguel Marco (2010-06-19): added eigenvalues, eigenvectors and minpoly functions


TESTS::

    sage: V = ZZ^2; f = V.hom([V.1, -2*V.0])
    sage: loads(dumps(f))
    Free module morphism defined by the matrix
    [ 0  1]
    [-2  0]
    Domain:   Ambient free module of rank 2 over the principal ideal domain ...
    Codomain: Ambient free module of rank 2 over the principal ideal domain ...
    sage: loads(dumps(f)) == f
    True
"""

####################################################################################
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
####################################################################################

# A matrix morphism is a morphism that is defined by multiplication by a
# matrix.  Elements of domain must either have a method "vector()" that
# returns a vector that the defining matrix can hit from the left, or
# be coercible into vector space of appropriate dimension.

import sage.modules.free_module as free_module

from sage.categories.finite_dimensional_modules_with_basis import FiniteDimensionalModulesWithBasis
from sage.categories.morphism import Morphism
from sage.modules import free_module_homspace, matrix_morphism
from sage.structure.richcmp import rich_to_bool, richcmp


def is_FreeModuleMorphism(x):
    """
    This function is deprecated.

    EXAMPLES::

        sage: V = ZZ^2; f = V.hom([V.1, -2*V.0])
        sage: sage.modules.free_module_morphism.is_FreeModuleMorphism(f)
        doctest:warning...
        DeprecationWarning: is_FreeModuleMorphism is deprecated;
        use isinstance(..., FreeModuleMorphism) or categories instead
        See https://github.com/sagemath/sage/issues/37731 for details.
        True
        sage: sage.modules.free_module_morphism.is_FreeModuleMorphism(0)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(37731,
                "is_FreeModuleMorphism is deprecated; "
                "use isinstance(..., FreeModuleMorphism) or categories instead")
    return isinstance(x, FreeModuleMorphism)


class FreeModuleMorphism(matrix_morphism.MatrixMorphism, FiniteDimensionalModulesWithBasis.Homsets.Endset.ElementMethods):

    def __init__(self, parent, A, side='left'):
        """
        INPUT:

            - ``parent`` -- a homspace in a (sub) category of free modules

            - ``A`` -- matrix

            - ``side`` -- side of the vectors acted on by the matrix  (default: ``'left'``)

        EXAMPLES::

            sage: V = ZZ^3; W = span([[1,2,3], [-1,2,8]], ZZ)
            sage: phi = V.hom(matrix(ZZ, 3, [1..9]))
            sage: type(phi)
            <class 'sage.modules.free_module_morphism.FreeModuleMorphism'>
        """
        if not isinstance(parent, free_module_homspace.FreeModuleHomspace):
            raise TypeError("parent (=%s) must be a free module hom space" % parent)
        if isinstance(A, matrix_morphism.MatrixMorphism):
            A = A.matrix()
        A = parent._matrix_space(side)(A)
        matrix_morphism.MatrixMorphism.__init__(self, parent, A, side=side)

    def pushforward(self, x):
        """
        Compute the image of a sub-module of the domain.

        EXAMPLES::

            sage: V = QQ^3; W = span([[1,2,3], [-1,2,5/3]], QQ)
            sage: phi = V.hom(matrix(QQ, 3, [1..9]))
            sage: phi.rank()
            2
            sage: phi(V)   #indirect doctest
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]

        We compute the image of a submodule of a ZZ-module embedded in
        a rational vector space::

            sage: V = QQ^3; W = V.span_of_basis([[2,2,3], [-1,2,5/3]], ZZ)
            sage: phi = W.hom([W.0, W.0 - W.1]); phi
            Free module morphism defined by the matrix
            [ 1  0]
            [ 1 -1]...
            sage: phi(span([2*W.1], ZZ))
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [  6   0 8/3]
            sage: phi(2*W.1)
            (6, 0, 8/3)
        """
        if isinstance(x, free_module.FreeModule_generic):
            V = self.domain().submodule(x)
            return self.restrict_domain(V).image()
        raise TypeError("`pushforward` is only defined for submodules")

    def _repr_(self):
        r"""
        Return string representation of this morphism of free modules.

        EXAMPLES::

            sage: V = ZZ^3; W = span([[1,2,3], [-1,2,8]], ZZ)
            sage: phi = V.hom(matrix(ZZ, 3, [1..9]))
            sage: phi._repr_()
            'Free module morphism defined by the matrix\n[1 2 3]\n[4 5 6]\n[7 8 9]\nDomain: Ambient free module of rank 3 over the principal ideal domain Integer Ring\nCodomain: Ambient free module of rank 3 over the principal ideal domain Integer Ring'

            sage: V = ZZ^6
            sage: W = ZZ^4
            sage: m = matrix(QQ, [[1, 0, 0 ,0], [0]*4, [0]*4, [0]*4, [0]*4, [0]*4])
            sage: phi = V.hom(m, W)
            sage: rho = phi.restrict_codomain(W.span([W.0]))
            sage: rho
            Free module morphism defined by the matrix
            [1]
            [0]
            [0]
            [0]
            [0]
            [0]
            Domain: Ambient free module of rank 6 over the principal ideal domain Integer Ring
            Codomain: Free module of degree 4 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1 0 0 0]

            sage: V = QQ^40
            sage: m = matrix(QQ, 40, 40, 1600)
            sage: phi = V.hom(m, V)
            sage: phi
            Vector space morphism represented by the matrix:
            40 x 40 dense matrix over Rational Field
            Domain:   Vector space of dimension 40 over Rational Field
            Codomain: Vector space of dimension 40 over Rational Field

        The representation displays which side of the vectors the matrix is acting::

            sage: V = ZZ^3
            sage: h = V.hom([V.1, V.2, V.0]); h
            Free module morphism defined by the matrix
            [0 1 0]
            [0 0 1]
            [1 0 0]
            Domain:   Ambient free module of rank 3 over the principal ideal domain Integer Ring
            Codomain: Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: h2 = V.hom([V.1, V.2, V.0], side='right'); h2
            Free module morphism defined as left-multiplication by the matrix
            [0 0 1]
            [1 0 0]
            [0 1 0]
            Domain:   Ambient free module of rank 3 over the principal ideal domain Integer Ring
            Codomain: Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        r = "Free module morphism defined {}by the matrix\n{!r}\nDomain: {}\nCodomain: {}"
        act = ""
        if self.side() == "right":
            act = "as left-multiplication "
        return r.format(act, self.matrix(), self.domain(), self.codomain())

    def change_ring(self, R):
        """
        Change the ring over which this morphism is defined.

        This changes the ring of the domain, codomain, and underlying matrix.

        EXAMPLES::

            sage: V0 = span([[0,0,1],[0,2,0]], ZZ); V1 = span([[1/2,0],[0,2]], ZZ)
            sage: W = span([[1,0],[0,6]], ZZ)
            sage: h = V0.hom([-3*V1.0 - 3*V1.1, -3*V1.0 - 3*V1.1])
            sage: h.base_ring()
            Integer Ring
            sage: h
            Free module morphism defined by the matrix
            [-3 -3]
            [-3 -3]...
            sage: h.change_ring(QQ).base_ring()
            Rational Field
            sage: f = h.change_ring(QQ); f
            Vector space morphism represented by the matrix:
            [-3 -3]
            [-3 -3]
            Domain:   Vector space of degree 3 and dimension 2 over Rational Field
                      Basis matrix:
                      [0 1 0]
                      [0 0 1]
            Codomain: Vector space of degree 2 and dimension 2 over Rational Field
                      Basis matrix:
                      [1 0]
                      [0 1]
            sage: f = h.change_ring(GF(7)); f
            Vector space morphism represented by the matrix:
            [4 4]
            [4 4]
            Domain:   Vector space of degree 3 and dimension 2 over Finite Field of size 7
                      Basis matrix:
                      [0 1 0]
                      [0 0 1]
            Codomain: Vector space of degree 2 and dimension 2 over Finite Field of size 7
                      Basis matrix:
                      [1 0]
                      [0 1]
        """
        D = self.domain().change_ring(R)
        C = self.codomain().change_ring(R)
        A = self.matrix().change_ring(R)
        return D.hom(A, C, side=self.side())

    def inverse_image(self, V):
        """
        Given a submodule V of the codomain of self, return the
        inverse image of V under self, i.e., the biggest submodule of
        the domain of ``self`` that maps into V.

        EXAMPLES:

        We test computing inverse images over a field::

            sage: V = QQ^3; W = span([[1,2,3], [-1,2,5/3]], QQ)
            sage: phi = V.hom(matrix(QQ, 3, [1..9]))
            sage: phi.rank()
            2
            sage: I = phi.inverse_image(W); I
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0    0]
            [   0    1 -1/2]
            sage: phi(I.0) in W
            True
            sage: phi(I.1) in W
            True
            sage: W = phi.image()
            sage: phi.inverse_image(W) == V
            True

        We test computing inverse images between two spaces embedded in different
        ambient spaces.::

            sage: V0 = span([[0,0,1],[0,2,0]],ZZ); V1 = span([[1/2,0],[0,2]],ZZ)
            sage: W = span([[1,0],[0,6]],ZZ)
            sage: h = V0.hom([-3*V1.0 - 3*V1.1, -3*V1.0 - 3*V1.1])
            sage: h.inverse_image(W)
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [0 2 1]
            [0 0 2]
            sage: h(h.inverse_image(W)).is_submodule(W)
            True
            sage: h(h.inverse_image(W)).index_in(W)
            +Infinity
            sage: h(h.inverse_image(W))
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 3 12]


        We test computing inverse images over the integers::

            sage: V = QQ^3; W = V.span_of_basis([[2,2,3],[-1,2,5/3]], ZZ)
            sage: phi = W.hom([W.0, W.0 - W.1])
            sage: Z = W.span([2*W.1]); Z
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [    2    -4 -10/3]
            sage: Y = phi.inverse_image(Z); Y
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [  6   0 8/3]
            sage: phi(Y) == Z
            True

        We test that :issue:`24590` is resolved::

            sage: A = FreeQuadraticModule(ZZ, 1, matrix([2]))
            sage: f = A.Hom(A).an_element()
            sage: f.inverse_image(A)
            Free module of degree 1 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1]

        We test that it respects the ``side``::

            sage: V = ZZ^2
            sage: m = matrix(2, [1, 1, 0, 1])
            sage: h = V.hom(m, side='right')
            sage: h
            Free module morphism defined as left-multiplication by the matrix
            [1 1]
            [0 1]...
            sage: SV = V.span([V.0])
            sage: h.inverse_image(SV)
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [1 0]
            sage: V.hom(m).inverse_image(SV)
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 1 -1]
        """
        if self.rank() == 0:
            # Special case -- if this is the 0 map, then the only possibility
            # for the inverse image is that it is the whole domain.
            return self.domain()

        R = self.base_ring()
        if self.side() == "left":
            A = self.matrix()
        else:
            A = self.matrix().transpose()

        # Replace the module V that we are going to pullback by a
        # submodule that is contained in the image of self, since our
        # plan is to lift all generators of V.
        V = self.image().intersection(V)
        # Write V in terms of the basis for the codomain.
        V = self.codomain().coordinate_module(V)
        B = V.basis_matrix()

        # Compute the kernel, which is contained in the inverse image.
        K = self.kernel()

        if R.is_field():
            # By solving, find lifts of each of the basis elements of V.
            # Each row of C gives a linear combination of the basis for the domain
            # that maps to one of the basis elements V.
            C = A.solve_left(B)

        else:
            if not hasattr(A, 'hermite_form'):
                raise NotImplementedError("base ring (%s) must have hermite_form algorithm in order to compute inverse image" % R)

            # 1. Compute H such that U*A = H = hnf(A) without zero
            # rows. What this "does" is find a basis for the image of
            # A and explicitly represents each element in this basis
            # as the image of some element of the domain (the rows of
            # U give these elements of the domain).
            H, U = A.hermite_form(transformation=True, include_zero_rows=False)

            # 2. Next we find the unique solution to the equation
            #    Y*H = B.  This writes each basis element of V in
            # terms of our image basis found in the previous step.
            Y = H.solve_left(B)

            # 3. Multiply Y by U then takes those same linear combinations
            # from step 2 above and lifts them to coefficients that define
            # linear combinations of the basis for the domain.
            C = Y*U

        # Finally take the linear combinations of the basis for the
        # domain defined by C. Together with the kernel K, this spans
        # the inverse image of V.
        dom = self.domain()
        if not dom.is_ambient():
            C = C * dom.basis_matrix()
        L = dom.submodule(C.rows())
        return K + L

    def lift(self, x):
        r"""
        Given an element of the image, return an element of the codomain that maps onto it.

        Note that ``lift`` and ``preimage_representative`` are
        equivalent names for this method, with the latter suggesting
        that the return value is a coset representative of the domain
        modulo the kernel of the morphism.

        EXAMPLES::

            sage: X = QQ**2
            sage: V = X.span([[2, 0], [0, 8]], ZZ)
            sage: W = (QQ**1).span([[1/12]], ZZ)
            sage: f = V.hom([W([1/3]), W([1/2])], W)
            sage: l=f.lift([1/3]); l  # random
            (8, -16)
            sage: f(l)
            (1/3)
            sage: f(f.lift([1/2]))
            (1/2)
            sage: f(f.lift([1/6]))
            (1/6)
            sage: f.lift([1/12])
            Traceback (most recent call last):
            ...
            ValueError: element is not in the image
            sage: f.lift([1/24])
            Traceback (most recent call last):
            ...
            TypeError: element [1/24] is not in free module

        This works for vector spaces, too::

            sage: V = VectorSpace(GF(3), 2)
            sage: W = VectorSpace(GF(3), 3)
            sage: f = V.hom([W.1, W.1 - W.0])
            sage: f.lift(W.1)
            (1, 0)
            sage: f.lift(W.2)
            Traceback (most recent call last):
            ...
            ValueError: element is not in the image
            sage: w = W((17, -2, 0))
            sage: f(f.lift(w)) == w
            True

        This example illustrates the use of the ``preimage_representative``
        as an equivalent name for this method.  ::

            sage: V = ZZ^3
            sage: W = ZZ^2
            sage: w = vector(ZZ, [1,2])
            sage: f = V.hom([w, w, w], W)
            sage: f.preimage_representative(vector(ZZ, [10, 20]))
            (0, 0, 10)

        ::

            sage: V = QQ^2; m = matrix(2, [1, 1, 0, 1])
            sage: V.hom(m, side='right').lift(V.0 + V.1)
            (0, 1)
            sage: V.hom(m).lift(V.0 + V.1)
            (1, 0)
        """
        from .free_module_element import vector
        x = self.codomain()(x)
        if self.side() == "right":
            A = self.matrix().transpose()
        else:
            A = self.matrix()
        R = self.base_ring()
        if R.is_field():
            try:
                C = A.solve_left(x)
            except ValueError:
                raise ValueError("element is not in the image")
        else:
            # see inverse_image for similar code but with comments
            if not hasattr(A, 'hermite_form'):
                raise NotImplementedError("base ring (%s) must have hermite_form algorithm in order to compute inverse image" % R)
            H, U = A.hermite_form(transformation=True, include_zero_rows=False)
            Y = H.solve_left(vector(self.codomain().coordinates(x)))
            C = Y*U
        try:
            t = self.domain().linear_combination_of_basis(C)
        except TypeError:
            raise ValueError("element is not in the image")
        assert self(t) == x
        return t

    preimage_representative = lift


class BaseIsomorphism1D(Morphism):
    """
    An isomorphism between a ring and a free rank-1 module over the ring.

    EXAMPLES::

        sage: R.<x,y> = QQ[]
        sage: V, from_V, to_V = R.free_module(R)
        sage: from_V
        Isomorphism morphism:
          From: Ambient free module of rank 1 over the integral domain
                Multivariate Polynomial Ring in x, y over Rational Field
          To:   Multivariate Polynomial Ring in x, y over Rational Field
    """
    def _repr_type(self):
        r"""
        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: V, from_V, to_V = R.free_module(R)
            sage: from_V._repr_type()
            'Isomorphism'
        """
        return "Isomorphism"

    def is_injective(self):
        r"""
        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: V, from_V, to_V = R.free_module(R)
            sage: from_V.is_injective()
            True
        """
        return True

    def is_surjective(self):
        r"""
        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: V, from_V, to_V = R.free_module(R)
            sage: from_V.is_surjective()
            True
        """
        return True

    def _richcmp_(self, other, op):
        r"""
        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: V, fr, to = R.free_module(R)
            sage: fr == loads(dumps(fr))
            True
        """
        if isinstance(other, BaseIsomorphism1D):
            return richcmp(self._basis, other._basis, op)
        else:
            return rich_to_bool(op, 1)


class BaseIsomorphism1D_to_FM(BaseIsomorphism1D):
    """
    An isomorphism from a ring to its 1-dimensional free module.

    INPUT:

    - ``parent`` -- the homset
    - ``basis`` -- (default: 1) an invertible element of the ring

    EXAMPLES::

        sage: R = Zmod(8)
        sage: V, from_V, to_V = R.free_module(R)
        sage: v = to_V(2); v
        (2)
        sage: from_V(v)
        2
        sage: W, from_W, to_W = R.free_module(R, basis=3)
        sage: W is V
        True
        sage: w = to_W(2); w
        (6)
        sage: from_W(w)
        2

    The basis vector has to be a unit so that the map is an isomorphism::

        sage: W, from_W, to_W = R.free_module(R, basis=4)
        Traceback (most recent call last):
        ...
        ValueError: basis element must be a unit
    """
    def __init__(self, parent, basis=None):
        """
        TESTS::

            sage: R = Zmod(8)
            sage: W, from_W, to_W = R.free_module(R, basis=3)
            sage: TestSuite(to_W).run()
        """
        Morphism.__init__(self, parent)
        self._basis = basis

    def _call_(self, x):
        """
        TESTS::

            sage: R = Zmod(8)
            sage: W, from_W, to_W = R.free_module(R, basis=3)
            sage: to_W(6) # indirect doctest
            (2)
        """
        if self._basis is not None:
            x *= self._basis
        return self.codomain()([x])


class BaseIsomorphism1D_from_FM(BaseIsomorphism1D):
    """
    An isomorphism to a ring from its 1-dimensional free module.

    INPUT:

    - ``parent`` -- the homset
    - ``basis`` -- (default: 1) an invertible element of the ring

    EXAMPLES::

        sage: R.<x> = QQ[[]]
        sage: V, from_V, to_V = R.free_module(R)
        sage: v = to_V(1+x); v
        (1 + x)
        sage: from_V(v)
        1 + x
        sage: W, from_W, to_W = R.free_module(R, basis=(1-x))
        sage: W is V
        True
        sage: w = to_W(1+x); w
        (1 - x^2)
        sage: from_W(w)
        1 + x + O(x^20)

    The basis vector has to be a unit so that the map is an isomorphism::

        sage: W, from_W, to_W = R.free_module(R, basis=x)
        Traceback (most recent call last):
        ...
        ValueError: basis element must be a unit
    """
    def __init__(self, parent, basis=None):
        """
        TESTS::

            sage: R.<x> = QQ[[]]
            sage: W, from_W, to_W = R.free_module(R, basis=(1-x))
            sage: TestSuite(from_W).run(skip='_test_nonzero_equal')
        """
        Morphism.__init__(self, parent)
        self._basis = basis

    def _call_(self, x):
        """
        TESTS::

            sage: R.<x> = QQ[[]]
            sage: W, from_W, to_W = R.free_module(R, basis=(1-x))
            sage: w = to_W(1+x); w
            (1 - x^2)
            sage: from_W(w)
            1 + x + O(x^20)
        """
        if self._basis is None:
            return x[0]
        return self.codomain()(x[0] / self._basis)
