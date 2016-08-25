"""
Ideals of Finite Algebras
"""
from __future__ import absolute_import

#*****************************************************************************
#  Copyright (C) 2011 Johan Bosman <johan.g.bosman@gmail.com>
#  Copyright (C) 2011, 2013 Peter Bruin <peter.bruin@math.uzh.ch>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .finite_dimensional_algebra_element import FiniteDimensionalAlgebraElement

from sage.matrix.constructor import Matrix
from sage.matrix.matrix import is_Matrix
from sage.rings.ideal import Ideal_generic
from sage.structure.element import parent
from sage.structure.sage_object import SageObject

from sage.misc.cachefunc import cached_method
from functools import reduce


class FiniteDimensionalAlgebraIdeal(Ideal_generic):
    """
    An ideal of a :class:`FiniteDimensionalAlgebra`.

    INPUT:

    - ``A`` -- a finite-dimensional algebra
    - ``gens`` -- the generators of this ideal
    - ``given_by_matrix`` -- (default: ``False``) whether the basis matrix is
      given by ``gens``

    EXAMPLES::

        sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
        sage: A.ideal(A([0,1]))
        Ideal (e1) of Finite-dimensional algebra of degree 2 over Finite Field of size 3
    """
    def __init__(self, A, gens=None, given_by_matrix=False):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.ideal(A([0,1]))
            sage: TestSuite(I).run(skip="_test_category") # Currently ideals are not using the category framework
        """
        k = A.base_ring()
        n = A.degree()
        if given_by_matrix:
            self._basis_matrix = gens
            gens = gens.rows()
        elif gens is None:
            self._basis_matrix = Matrix(k, 0, n)
        elif isinstance(gens, (list, tuple)):
            B = [FiniteDimensionalAlgebraIdeal(A, x).basis_matrix() for x in gens]
            B = reduce(lambda x, y: x.stack(y), B, Matrix(k, 0, n))
            self._basis_matrix = B.echelon_form().image().basis_matrix()
        elif is_Matrix(gens):
            gens = FiniteDimensionalAlgebraElement(A, gens)
        elif isinstance(gens, FiniteDimensionalAlgebraElement):
            gens = gens.vector()
            B = Matrix([gens * b for b in A.table()])
            self._basis_matrix = B.echelon_form().image().basis_matrix()
        Ideal_generic.__init__(self, A, gens)

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.ideal(A([1,1]))
            sage: J = A.ideal(A([0,1]))
            sage: I == J
            False
            sage: I == I
            True
            sage: I == I + J
            True

            sage: A2 = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: A is A2
            True
            sage: A == A2
            True
            sage: I2 = A.ideal(A([1,1]))
            sage: I == I2
            True
        """
        if self is other:
            return True
        if not isinstance(other, FiniteDimensionalAlgebraIdeal):
            return False
        if self.ring() != other.ring():
            return False
        return self.basis_matrix() == other.basis_matrix()

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.ideal(A([1,1]))
            sage: J = A.ideal(A([0,1]))
            sage: I != J
            True
            sage: I != I
            False
            sage: I != I + J
            False
        """
        return not self == other

    def __contains__(self, elt):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: J = A.ideal(A([0,1]))
            sage: A([0,1]) in J
            True
            sage: A([1,0]) in J
            False
        """
        if self.ring() is not parent(elt):
            return False
        return elt.vector() in self.vector_space()

    def __le__(self, other):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.ideal(A([1,1]))
            sage: J = A.ideal(A([0,1]))
            sage: I <= J
            False
            sage: I <= I
            True
            sage: I <= I + J
            True
        """
        if self is other:
            return True
        if self.ring() is not other.ring():
            return False
        return self.vector_space().is_subspace(other.vector_space())

    def __lt__(self, other):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.ideal(A([1,1]))
            sage: J = A.ideal(A([0,1]))
            sage: I < J
            False
            sage: I < I
            False
            sage: I < I + J
            False
        """
        return self != other and self <= other

    def __ge__(self, other):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.ideal(A([1,1]))
            sage: J = A.ideal(A([0,1]))
            sage: I >= J
            True
            sage: I >= I
            True
            sage: I >= I + J
            True
        """
        return other <= self

    def __gt__(self, other):
        """
        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.ideal(A([1,1]))
            sage: J = A.ideal(A([0,1]))
            sage: I > J
            True
            sage: I > I
            False
            sage: I > I + J
            False
        """
        return other < self

    def basis_matrix(self):
        """
        Return the echelonized matrix whose rows form a basis of ``self``.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.ideal(A([1,1]))
            sage: I.basis_matrix()
            [1 0]
            [0 1]
        """
        return self._basis_matrix

    @cached_method
    def vector_space(self):
        """
        Return ``self`` as a vector space.

        EXAMPLES::

            sage: A = FiniteDimensionalAlgebra(GF(3), [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.ideal(A([1,1]))
            sage: I.vector_space()
            Vector space of degree 2 and dimension 2 over Finite Field of size 3
            Basis matrix:
            [1 0]
            [0 1]
        """
        return self.basis_matrix().image()

