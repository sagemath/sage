r"""
Examples of Hopf algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.categories.all import HopfAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.all import tensor

class MyGroupAlgebra(CombinatorialFreeModule):
    r"""
    An example of a Hopf algebra with basis: the group algebra of a group

    This class illustrates a minimal implementation of a Hopf algebra with basis.
    """

    def __init__(self, R, G):
        """
        EXAMPLES::

            sage: from sage.categories.examples.hopf_algebras_with_basis import MyGroupAlgebra
            sage: A = MyGroupAlgebra(QQ, DihedralGroup(6))
            sage: A.category()
            Category of finite dimensional hopf algebras with basis over Rational Field
            sage: TestSuite(A).run()
        """
        self._group = G
        CombinatorialFreeModule.__init__(self, R, G, category = HopfAlgebrasWithBasis(R))

    def _repr_(self):
        """
        EXAMPLES::

            sage: HopfAlgebrasWithBasis(QQ).example() # indirect doctest
            An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
        """
        return "An example of Hopf algebra with basis: the group algebra of the %s over %s"%(self._group, self.base_ring())

    @cached_method
    def one_basis(self):
        """
        Returns the one of the group, which index the one of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = HopfAlgebrasWithBasis(QQ).example()
            sage: A.one_basis()
            ()
            sage: A.one()
            B[()]
        """
        return self._group.one()

    def product_on_basis(self, g1, g2):
        r"""
        Product, on basis elements, as per
        :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis`.

        The product of two basis elements is induced by the product of
        the corresponding elements of the group.

        EXAMPLES::

            sage: A = HopfAlgebrasWithBasis(QQ).example()
            sage: (a, b) = A._group.gens()
            sage: a*b
            (1,2)
            sage: A.product_on_basis(a, b)
            B[(1,2)]
        """
        return self.basis()[g1 * g2]

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra, as per :meth:`~.magmatic_algebras.MagmaticAlgebras.ParentMethods.algebra_generators`.

        They correspond to the generators of the group.

        EXAMPLES::

            sage: A = HopfAlgebrasWithBasis(QQ).example(); A
            An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
            sage: A.algebra_generators()
            Finite family {(1,2,3): B[(1,2,3)], (1,3): B[(1,3)]}
        """
        return Family(self._group.gens(), self.monomial)

    def coproduct_on_basis(self, g):
        r"""
        Coproduct, on basis elements, as per :meth:`HopfAlgebrasWithBasis.ParentMethods.coproduct_on_basis`.

        The basis elements are group like: `\Delta(g) = g \otimes g`.

        EXAMPLES::

            sage: A = HopfAlgebrasWithBasis(QQ).example()
            sage: (a, b) = A._group.gens()
            sage: A.coproduct_on_basis(a)
            B[(1,2,3)] # B[(1,2,3)]
        """
        g = self.monomial(g)
        return tensor([g, g])

    def counit_on_basis(self, g):
        r"""
        Counit, on basis elements, as per :meth:`HopfAlgebrasWithBasis.ParentMethods.counit_on_basis`.

        The counit on the basis elements is 1.

        EXAMPLES::

            sage: A = HopfAlgebrasWithBasis(QQ).example()
            sage: (a, b) = A._group.gens()
            sage: A.counit_on_basis(a)
            1
        """
        return self.base_ring().one()

    def antipode_on_basis(self, g):
        r"""
        Antipode, on basis elements, as per :meth:`HopfAlgebrasWithBasis.ParentMethods.antipode_on_basis`.

        It is given, on basis elements, by `\nu(g) = g^{-1}`

        EXAMPLES::

            sage: A = HopfAlgebrasWithBasis(QQ).example()
            sage: (a, b) = A._group.gens()
            sage: A.antipode_on_basis(a)
            B[(1,3,2)]
        """
        return self.monomial(~g)
