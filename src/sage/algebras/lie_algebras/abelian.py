"""
Abelian Lie Algebras

AUTHORS:

- Travis Scrimshaw (2016-06-07): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.indexed_generators import (IndexedGenerators,
                                               standardize_names_index_set)
from sage.categories.lie_algebras import LieAlgebras
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.algebras.lie_algebras.lie_algebra import InfinitelyGeneratedLieAlgebra
from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.infinity import infinity
from sage.sets.family import Family


class AbelianLieAlgebra(LieAlgebraWithStructureCoefficients):
    r"""
    An abelian Lie algebra.

    A Lie algebra `\mathfrak{g}` is abelian if `[x, y] = 0` for all
    `x, y \in \mathfrak{g}`.

    EXAMPLES::

        sage: L.<x, y> = LieAlgebra(QQ, abelian=True)
        sage: L.bracket(x, y)
        0
    """
    @staticmethod
    def __classcall_private__(cls, R, names=None, index_set=None, category=None, **kwds):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: L1 = LieAlgebra(QQ, 'x,y', {})
            sage: L2.<x, y> = LieAlgebra(QQ, abelian=True)
            sage: L1 is L2
            True
        """
        names, index_set = standardize_names_index_set(names, index_set)
        if index_set.cardinality() == infinity:
            return InfiniteDimensionalAbelianLieAlgebra(R, index_set, **kwds)
        return super().__classcall__(cls, R, names, index_set, category=category, **kwds)

    def __init__(self, R, names, index_set, category, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
            sage: TestSuite(L).run()
        """
        cat = LieAlgebras(R).FiniteDimensional().WithBasis().Nilpotent()
        category = cat.or_subcategory(category)
        LieAlgebraWithStructureCoefficients.__init__(self, R, Family({}), names,
                                                     index_set, category, **kwds)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, 3, 'x', abelian=True)
            Abelian Lie algebra on 3 generators (x0, x1, x2) over Rational Field
        """
        gens = self.lie_algebra_generators()
        if gens.cardinality() == 1:
            return "Abelian Lie algebra on generator {} over {}".format(tuple(gens)[0], self.base_ring())
        return "Abelian Lie algebra on {} generators {} over {}".format(
            gens.cardinality(), tuple(gens), self.base_ring())

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
            sage: L._construct_UEA()
            Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        """
        return PolynomialRing(self.base_ring(), self.variable_names())

    def is_abelian(self):
        """
        Return ``True`` since ``self`` is an abelian Lie algebra.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
            sage: L.is_abelian()
            True
        """
        return True

    # abelian => nilpotent => solvable
    is_nilpotent = is_solvable = is_abelian

    class Element(LieAlgebraWithStructureCoefficients.Element):
        def _bracket_(self, y):
            """
            Return the Lie bracket ``[self, y]``.

            EXAMPLES::

                sage: L.<x, y> = LieAlgebra(QQ, abelian=True)
                sage: L.bracket(x, y)
                0
            """
            return self.parent().zero()


class InfiniteDimensionalAbelianLieAlgebra(InfinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    An infinite dimensional abelian Lie algebra.

    A Lie algebra `\mathfrak{g}` is abelian if `[x, y] = 0` for all
    `x, y \in \mathfrak{g}`.
    """
    def __init__(self, R, index_set, prefix='L', **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, index_set=ZZ, abelian=True)
            sage: TestSuite(L).run()
        """
        cat = LieAlgebras(R).WithBasis()
        InfinitelyGeneratedLieAlgebra.__init__(self, R, category=cat)
        IndexedGenerators.__init__(self, index_set, prefix=prefix, **kwds)

    def dimension(self):
        r"""
        Return the dimension of ``self``, which is `\infty`.

        EXAMPLES::

            sage: L = lie_algebras.abelian(QQ, index_set=ZZ)
            sage: L.dimension()
            +Infinity
        """
        return infinity

    def is_abelian(self):
        """
        Return ``True`` since ``self`` is an abelian Lie algebra.

        EXAMPLES::

            sage: L = lie_algebras.abelian(QQ, index_set=ZZ)
            sage: L.is_abelian()
            True
        """
        return True

    # abelian => nilpotent => solvable
    is_nilpotent = is_solvable = is_abelian

    # For compatibility with CombinatorialFreeModuleElement
    _repr_term = IndexedGenerators._repr_generator
    _latex_term = IndexedGenerators._latex_generator

    class Element(LieAlgebraElement):
        def _bracket_(self, other):
            """
            Return the Lie bracket ``[self, y]``.

            EXAMPLES::

                sage: L = lie_algebras.abelian(QQ, index_set=ZZ)
                sage: B = L.basis()
                sage: l1 = B[1]
                sage: l5 = B[5]
                sage: l1.bracket(l5)
                0
            """
            return self.parent().zero()
