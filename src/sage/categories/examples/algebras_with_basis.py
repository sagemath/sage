r"""
Examples of algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.categories.all import AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.words import Words

class FreeAlgebra(CombinatorialFreeModule):
    r"""
    An example of an algebra with basis: the free algebra

    This class illustrates a minimal implementation of an algebra with basis.
    """

    def __init__(self, R, alphabet = ("a", "b", "c")):
        """
        EXAMPLES::

            sage: A = AlgebrasWithBasis(QQ).example(); A
            An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
            sage: TestSuite(A).run()

        """
        self._alphabet = alphabet
        CombinatorialFreeModule.__init__(self, R, Words(alphabet, infinite=False), category = AlgebrasWithBasis(R))

    def _repr_(self):
        """
        EXAMPLES::

            sage: AlgebrasWithBasis(QQ).example() # indirect doctest
            An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
        """
        return "An example of an algebra with basis: the free algebra on the generators %s over %s"%(self._alphabet, self.base_ring())

    @cached_method
    def one_basis(self):
        """
        Returns the empty word, which index the one of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::r

            sage: A = AlgebrasWithBasis(QQ).example()
            sage: A.one_basis()
            word:
            sage: A.one()
            B[word: ]
        """
        return self.basis().keys()([])

    def product_on_basis(self, w1, w2):
        r"""
        Product of basis elements, as per
        :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis`.

        EXAMPLES::

            sage: A = AlgebrasWithBasis(QQ).example()
            sage: Words = A.basis().keys()
            sage: A.product_on_basis(Words("acb"), Words("cba"))
            B[word: acbcba]
            sage: (a,b,c) = A.algebra_generators()
            sage: a * (1-b)^2 * c
            B[word: abbc] - 2*B[word: abc] + B[word: ac]
        """
        return self.basis()[w1 + w2]

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra, as per :meth:`~.magmatic_algebras.MagmaticAlgebras.ParentMethods.algebra_generators`.

        EXAMPLES::

            sage: A = AlgebrasWithBasis(QQ).example(); A
            An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
            sage: A.algebra_generators()
            Family (B[word: a], B[word: b], B[word: c])
        """
        Words = self.basis().keys()
        return Family( [self.monomial(Words(a)) for a in self._alphabet] )
        # FIXME: use this once the keys argument of FiniteFamily will be honoured
        # for the specifying the order of the elements in the family
        #return Family(self._alphabet, lambda a: self.term(self.basis().keys()(a)))

Example = FreeAlgebra
