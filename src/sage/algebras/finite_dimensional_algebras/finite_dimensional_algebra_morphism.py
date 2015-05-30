"""
Morphisms Between Finite Algebras
"""

#*****************************************************************************
#  Copyright (C) 2011 Johan Bosman <johan.g.bosman@gmail.com>
#  Copyright (C) 2011, 2013 Peter Bruin <peter.bruin@math.uzh.ch>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.homset import Hom
from sage.rings.morphism import RingHomomorphism_im_gens
from sage.rings.homset import RingHomset_generic
from sage.matrix.constructor import matrix
from sage.matrix.matrix import is_Matrix

class FiniteDimensionalAlgebraMorphism(RingHomomorphism_im_gens):
    """
    Create a morphism between two :class:`finite-dimensional algebras <FiniteDimensionalAlgebra>`.

    INPUT:

    - ``parent`` -- the parent homset

    - ``f`` -- matrix of the underlying `k`-linear map

    - ``unitary`` -- boolean (default: ``True``); if ``True`` and ``check``
      is also ``True``, raise a ``ValueError`` unless ``A`` and ``B`` are
      unitary and ``f`` respects unit elements

    - ``check`` -- boolean (default: ``True``); check whether the given
      `k`-linear map really defines a (not necessarily unitary)
      `k`-algebra homomorphism

    The algebras ``A`` and ``B`` must be defined over the same base field.

    EXAMPLES::

        sage: from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_morphism import FiniteDimensionalAlgebraMorphism
        sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
        sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([1])])
        sage: H = Hom(A, B)
        sage: f = H(Matrix([[1], [0]]))
        sage: f.domain() is A
        True
        sage: f.codomain() is B
        True
        sage: f(A.basis()[0])
        e
        sage: f(A.basis()[1])
        0

    .. TODO:: An example illustrating unitary flag.
    """
    def __init__(self, parent, f, check=True, unitary=True):
        """
        TESTS::

            sage: from sage.algebras.finite_dimensional_algebras.finite_dimensional_algebra_morphism import FiniteDimensionalAlgebraMorphism
            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([1])])
            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: H = Hom(A, B)
            sage: phi = FiniteDimensionalAlgebraMorphism(H, Matrix([[1, 0]]))
            sage: TestSuite(phi).run(skip="_test_category")
        """
        A = parent.domain()
        B = parent.codomain()

        RingHomomorphism_im_gens.__init__(self, parent=parent, im_gens=f.rows(), check=check)
        self._matrix = f

        if unitary and check and (not A.is_unitary()
                                  or not B.is_unitary()
                                  or self(A.one()) != B.one()):
            raise ValueError("homomorphism does not respect unit elements")

    def _repr_(self):
        """
        TESTS::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.maximal_ideal()
            sage: q = A.quotient_map(I)
            sage: q._repr_()
            'Morphism from Finite-dimensional algebra of degree 2 over Rational Field to Finite-dimensional algebra of degree 1 over Rational Field given by matrix\n[1]\n[0]'
        """
        return "Morphism from {} to {} given by matrix\n{}".format(
                self.domain(), self.codomain(), self._matrix)

    def __call__(self, x):
        """
        TESTS::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.maximal_ideal()
            sage: q = A.quotient_map(I)
            sage: q(0) == 0 and q(1) == 1
            True
        """
        x = self.domain()(x)
        B = self.codomain()
        return B.element_class(B, x.vector() * self._matrix)

    def __eq__(self, other):
        """
        Check equality.

        TESTS::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([1])])
            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: H = Hom(A, B)
            sage: phi = H(Matrix([[1, 0]]))
            sage: psi = H(Matrix([[1, 0]]))
            sage: phi == psi
            True
            sage: phi == H.zero()
            False
        """
        return (isinstance(other, FiniteDimensionalAlgebraMorphism)
                and self.parent() == other.parent()
                and self._matrix == other._matrix)

    def __ne__(self, other):
        """
        Check not equals.

        TESTS::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([1])])
            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: H = Hom(A, B)
            sage: phi = H(Matrix([[1, 0]]))
            sage: psi = H(Matrix([[1, 0]]))
            sage: phi != psi
            False
            sage: phi != H.zero()
            True
        """
        return not self.__eq__(other)

    def matrix(self):
        """
        Return the matrix of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([1])])
            sage: M = Matrix([[1], [0]])
            sage: H = Hom(A, B)
            sage: f = H(M)
            sage: f.matrix() == M
            True
        """
        return self._matrix

    def inverse_image(self, I):
        """
        Return the inverse image of ``I`` under ``self``.

        INPUT:

        - ``I`` -- ``FiniteDimensionalAlgebraIdeal``, an ideal of ``self.codomain()``

        OUTPUT:

        -- ``FiniteDimensionalAlgebraIdeal``, the inverse image of `I` under ``self``.

        EXAMPLE::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.maximal_ideal()
            sage: q = A.quotient_map(I)
            sage: B = q.codomain()
            sage: q.inverse_image(B.zero_ideal()) == I
            True
        """
        coker_I = I.basis_matrix().transpose().kernel().basis_matrix().transpose()
        return self.domain().ideal((self._matrix * coker_I).kernel().basis_matrix(), given_by_matrix=True)

class FiniteDimensionalAlgebraHomset(RingHomset_generic):
    """
    Set of morphisms between two finite-dimensional algebras.
    """
    @cached_method
    def zero(self):
        """
        Construct the zero morphism of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([1])])
            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: H = Hom(A, B)
            sage: H.zero()
            Morphism from Finite-dimensional algebra of degree 1 over Rational Field to
             Finite-dimensional algebra of degree 2 over Rational Field given by matrix
            [0 0]
        """
        from sage.matrix.constructor import matrix
        return FiniteDimensionalAlgebraMorphism(self, matrix.zero(self.domain().ngens(),
                                                       self.codomain().ngens()),
                                     False, False)

    def __call__(self, f, check=True, unitary=True):
        """
        Construct a homomorphism.

        .. TODO::

            Implement taking generator images and converting them to a matrix.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(QQ, [Matrix([1])])
            sage: B = FiniteDimensionalAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: H = Hom(A, B)
            sage: H(Matrix([[1, 0]]))
            Morphism from Finite-dimensional algebra of degree 1 over Rational Field to
             Finite-dimensional algebra of degree 2 over Rational Field given by matrix
            [1 0]
        """
        if isinstance(f, FiniteDimensionalAlgebraMorphism):
            if f.parent() is self:
                return f
            if f.parent() == self:
                return FiniteDimensionalAlgebraMorphism(self, f._matrix, check, unitary)
        elif is_Matrix(f):
            return FiniteDimensionalAlgebraMorphism(self, f, check, unitary)
        try:
            from sage.matrix.constructor import Matrix
            return FiniteDimensionalAlgebraMorphism(self, Matrix(f), check, unitary)
        except Exception:
            return RingHomset_generic.__call__(self, f, check)

