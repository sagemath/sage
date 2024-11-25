# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.combinat sage.modules
r"""
Examples of a Lie algebra with basis
"""
#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.categories.lie_algebras import LieAlgebras
from sage.categories.algebras import Algebras
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.combinat.free_module import CombinatorialFreeModule


class AbelianLieAlgebra(CombinatorialFreeModule):
    r"""
    An example of a Lie algebra: the abelian Lie algebra.

    This class illustrates a minimal implementation of a Lie algebra with
    a distinguished basis.
    """
    def __init__(self, R, gens):
        """
        EXAMPLES::

            sage: L = LieAlgebras(QQ).WithBasis().example()
            sage: TestSuite(L).run()
        """
        cat = LieAlgebras(R).WithBasis()
        CombinatorialFreeModule.__init__(self, R, gens, category=cat)

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).WithBasis().example()
            sage: L._construct_UEA()
            Polynomial algebra with generators indexed by Partitions over Rational Field
        """
        return IndexedPolynomialRing(self.base_ring(), self._indices)

    def _repr_(self):
        """
        EXAMPLES::

            sage: LieAlgebras(QQ).WithBasis().example()
            An example of a Lie algebra: the abelian Lie algebra on the
             generators indexed by Partitions over Rational Field
        """
        return "An example of a Lie algebra: the abelian Lie algebra on the" \
               " generators indexed by {} over {}".format(
                        self.basis().keys(), self.base_ring())

    def lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a Lie algebra.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).WithBasis().example()
            sage: L.lie_algebra_generators()
            Lazy family (Term map from Partitions to
             An example of a Lie algebra: the abelian Lie algebra on the
             generators indexed by Partitions over Rational
             Field(i))_{i in Partitions}
        """
        return self.basis()

    def bracket_on_basis(self, x, y):
        """
        Return the Lie bracket on basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).WithBasis().example()
            sage: L.bracket_on_basis(Partition([4,1]), Partition([2,2,1]))
            0
        """
        return self.zero()

    class Element(CombinatorialFreeModule.Element):
        def lift(self):
            """
            Return the lift of ``self`` to the universal enveloping algebra.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).WithBasis().example()
                sage: elt = L.an_element()
                sage: elt.lift()
                3*P[F[2]] + 2*P[F[1]] + 2*P[F[]]
            """
            UEA = self.parent().universal_enveloping_algebra()
            I = UEA._indices
            return UEA.sum_of_terms((I.gen(t), c) for t, c in self)


Example = AbelianLieAlgebra

##############


class IndexedPolynomialRing(CombinatorialFreeModule):
    """
    Polynomial ring whose generators are indexed by an arbitrary set.

    .. TODO::

        Currently this is just used as the universal enveloping algebra
        for the example of the abelian Lie algebra. This should be
        factored out into a more complete class.
    """
    def __init__(self, R, indices, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).WithBasis().example()
            sage: UEA = L.universal_enveloping_algebra()
            sage: TestSuite(UEA).run()
        """
        if 'category' not in kwds:
            kwds['category'] = Algebras(R).WithBasis()
        if 'prefix' not in kwds:
            kwds['prefix'] = 'P'
        # This is a workaround until IndexedFree(Abelian)Monoid elements compare properly
        kwds['sorting_key'] = lambda x: x.to_word_list()
        kwds['sorting_reverse'] = True
        M = IndexedFreeAbelianMonoid(indices, bracket='')
        CombinatorialFreeModule.__init__(self, R, M, **kwds)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).WithBasis().example()
            sage: L.universal_enveloping_algebra()
            Polynomial algebra with generators indexed by Partitions over Rational Field
        """
        return "Polynomial algebra with generators indexed by {} over {}".format(
            self._indices._indices, self.base_ring())

    def one_basis(self):
        """
        Return the index of element `1`.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).WithBasis().example()
            sage: UEA = L.universal_enveloping_algebra()
            sage: UEA.one_basis()
            1
            sage: UEA.one_basis().parent()
            Free abelian monoid indexed by Partitions
        """
        return self._indices.one()

    def product_on_basis(self, x, y):
        """
        Return the product of the monomials indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).WithBasis().example()
            sage: UEA = L.universal_enveloping_algebra()
            sage: I = UEA._indices
            sage: UEA.product_on_basis(I.an_element(), I.an_element())
            P[F[]^4*F[1]^4*F[2]^6]
        """
        return self.monomial(x*y)

    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).WithBasis().example()
            sage: UEA = L.universal_enveloping_algebra()
            sage: UEA.algebra_generators()
            Lazy family (algebra generator map(i))_{i in Partitions}
        """
        I = self._indices
        return Family(I._indices, lambda x: self.monomial(I.gen(x)),
                      name="algebra generator map")
