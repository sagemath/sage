# sage_setup: distribution = sagemath-categories
r"""
Group Algebras

This module implements the category of group algebras for arbitrary
groups over arbitrary commutative rings. For details, see
:mod:`sage.categories.algebra_functor`.

AUTHOR:

- David Loeffler (2008-08-24): initial version
- Martin Raum (2009-08): update to use new coercion model -- see
  :issue:`6670`.
- John Palmieri (2011-07): more updates to coercion, categories, etc.,
  group algebras constructed using CombinatorialFreeModule -- see
  :issue:`6670`.
- Nicolas M. Thiéry (2010-2017), Travis Scrimshaw (2017):
  generalization to a covariant functorial construction for
  monoid algebras, and beyond -- see e.g. :issue:`18700`.
"""

#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2017 Nicolas M. Thiéry <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.algebra_functor import AlgebrasCategory


class GroupAlgebras(AlgebrasCategory):
    r"""
    The category of group algebras over a given base ring.

    EXAMPLES::

        sage: C = Groups().Algebras(ZZ); C
        Category of group algebras over Integer Ring
        sage: C.super_categories()
        [Category of Hopf algebras with basis over Integer Ring,
         Category of monoid algebras over Integer Ring]

    We can also construct this category with::

        sage: C is GroupAlgebras(ZZ)
        True

    Here is how to create the group algebra of a group `G`::

        sage: G = DihedralGroup(5)                                                      # needs sage.groups
        sage: QG = G.algebra(QQ); QG                                                    # needs sage.groups sage.modules
        Algebra of
         Dihedral group of order 10 as a permutation group over Rational Field

    and an example of computation::

        sage: g = G.an_element(); g                                                     # needs sage.groups sage.modules
        (1,4)(2,3)
        sage: (QG.term(g) + 1)**3                                                       # needs sage.groups sage.modules
        4*() + 4*(1,4)(2,3)

    .. TODO::

        - Check which methods would be better located in
          ``Monoid.Algebras`` or ``Groups.Finite.Algebras``.

    TESTS::

        sage: # needs sage.groups sage.modules
        sage: A = GroupAlgebras(QQ).example(GL(3, GF(11)))
        sage: A.one_basis()
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: A = SymmetricGroupAlgebra(QQ, 4)                                          # needs sage.combinat
        sage: x = Permutation([4,3,2,1])
        sage: A.product_on_basis(x, x)                                                  # needs sage.combinat
        [1, 2, 3, 4]

        sage: C = GroupAlgebras(ZZ)
        sage: TestSuite(C).run()
    """
    def extra_super_categories(self):
        """
        Implement the fact that the algebra of a group is a Hopf
        algebra.

        EXAMPLES::

            sage: C = Groups().Algebras(QQ)
            sage: C.extra_super_categories()
            [Category of Hopf algebras over Rational Field]
            sage: sorted(C.super_categories(), key=str)
            [Category of Hopf algebras with basis over Rational Field,
             Category of monoid algebras over Rational Field]
        """
        from sage.categories.hopf_algebras import HopfAlgebras
        return [HopfAlgebras(self.base_ring())]

    def example(self, G=None):
        """
        Return an example of group algebra.

        EXAMPLES::

            sage: GroupAlgebras(QQ['x']).example()                                      # needs sage.groups sage.modules
            Algebra of Dihedral group of order 8 as a permutation group
             over Univariate Polynomial Ring in x over Rational Field

        An other group can be specified as optional argument::

            sage: GroupAlgebras(QQ).example(AlternatingGroup(4))                        # needs sage.groups sage.modules
            Algebra of
             Alternating group of order 4!/2 as a permutation group over Rational Field
        """
        from sage.groups.perm_gps.permgroup_named import DihedralGroup
        if G is None:
            G = DihedralGroup(4)
        return G.algebra(self.base_ring())

    class ParentMethods:
        def __init_extra__(self):
            """
            Enable coercion from the defining group.

            EXAMPLES::

                sage: # needs sage.combinat sage.groups sage.modules
                sage: A = GroupAlgebra(SymmetricGroup(4), QQ)
                sage: B = GroupAlgebra(SymmetricGroup(3), ZZ)
                sage: A.has_coerce_map_from(B)
                True
                sage: B.has_coerce_map_from(A)
                False
                sage: A.has_coerce_map_from(ZZ)
                True
                sage: A.has_coerce_map_from(CC)
                False
                sage: A.has_coerce_map_from(SymmetricGroup(5))
                False
                sage: A.has_coerce_map_from(SymmetricGroup(2))
                True
            """
            if not self.base_ring().has_coerce_map_from(self.group()):
                ## some matrix groups assume that coercion is only valid to
                ## other matrix groups. This is a workaround
                ## call _element_constructor_ to coerce group elements
                #try:
                self._populate_coercion_lists_(coerce_list=[self.group()])

        def _latex_(self):
            r"""
            Latex string of ``self``.

            EXAMPLES::

                sage: A = GroupAlgebra(KleinFourGroup(), ZZ)                            # needs sage.groups sage.modules
                sage: latex(A)  # indirect doctest                                      # needs sage.groups sage.modules
                \Bold{Z}[\langle (3,4), (1,2) \rangle]
            """
            from sage.misc.latex import latex
            return "%s[%s]" % (latex(self.base_ring()), latex(self.group()))

        def group(self):
            r"""
            Return the underlying group of the group algebra.

            EXAMPLES::

                sage: GroupAlgebras(QQ).example(GL(3, GF(11))).group()                  # needs sage.groups sage.modules
                General Linear Group of degree 3 over Finite Field of size 11
                sage: SymmetricGroup(10).algebra(QQ).group()                            # needs sage.combinat sage.groups sage.modules
                Symmetric group of order 10! as a permutation group
            """
            return self.basis().keys()

        @cached_method
        def center_basis(self):
            r"""
            Return a basis of the center of the group algebra.

            The canonical basis of the center of the group algebra
            is the family `(f_\sigma)_{\sigma\in C}`, where `C` is
            any collection of representatives of the conjugacy
            classes of the group, and `f_\sigma` is the sum of the
            elements in the conjugacy class of `\sigma`.

            OUTPUT: tuple of elements of ``self``

            .. WARNING::

                - This method requires the underlying group to
                  have a method ``conjugacy_classes``
                  (every permutation group has one, thanks GAP!).

            EXAMPLES::

                sage: SymmetricGroup(3).algebra(QQ).center_basis()                      # needs sage.combinat sage.groups sage.modules
                ((), (2,3) + (1,2) + (1,3), (1,2,3) + (1,3,2))

            .. SEEALSO::

                - :meth:`Groups.Algebras.ElementMethods.central_form`
                - :meth:`Monoids.Algebras.ElementMethods.is_central`
            """
            return tuple([self.sum_of_monomials(conj) for conj in
                          self.basis().keys().conjugacy_classes()])

        # Hopf algebra structure

        def coproduct_on_basis(self, g):
            r"""
            Return the coproduct of the element ``g`` of the basis.

            Each basis element ``g`` is group-like. This method is
            used to compute the coproduct of any element.

            EXAMPLES::

                sage: # needs sage.groups sage.modules
                sage: A = CyclicPermutationGroup(6).algebra(ZZ); A
                Algebra of
                 Cyclic group of order 6 as a permutation group over Integer Ring
                sage: g = CyclicPermutationGroup(6).an_element(); g
                (1,2,3,4,5,6)
                sage: A.coproduct_on_basis(g)
                (1,2,3,4,5,6) # (1,2,3,4,5,6)
                sage: a = A.an_element(); a
                () + 3*(1,2,3,4,5,6) + 3*(1,3,5)(2,4,6)
                sage: a.coproduct()
                () # () + 3*(1,2,3,4,5,6) # (1,2,3,4,5,6) + 3*(1,3,5)(2,4,6) # (1,3,5)(2,4,6)
            """
            from sage.categories.tensor import tensor
            g = self.term(g)
            return tensor([g, g])

        def antipode_on_basis(self,g):
            r"""
            Return the antipode of the element ``g`` of the basis.

            Each basis element ``g`` is group-like, and so has
            antipode `g^{-1}`. This method is used to compute the
            antipode of any element.

            EXAMPLES::

                sage: # needs sage.groups sage.modules
                sage: A = CyclicPermutationGroup(6).algebra(ZZ); A
                Algebra of
                 Cyclic group of order 6 as a permutation group over Integer Ring
                sage: g = CyclicPermutationGroup(6).an_element(); g
                (1,2,3,4,5,6)
                sage: A.antipode_on_basis(g)
                (1,6,5,4,3,2)
                sage: a = A.an_element(); a
                () + 3*(1,2,3,4,5,6) + 3*(1,3,5)(2,4,6)
                sage: a.antipode()
                () + 3*(1,5,3)(2,6,4) + 3*(1,6,5,4,3,2)
            """
            return self.term(~g)

        def counit_on_basis(self,g):
            r"""
            Return the counit of the element ``g`` of the basis.

            Each basis element ``g`` is group-like, and so has
            counit `1`. This method is used to compute the
            counit of any element.

            EXAMPLES::

                sage: A = CyclicPermutationGroup(6).algebra(ZZ); A                      # needs sage.groups sage.modules
                Algebra of
                 Cyclic group of order 6 as a permutation group over Integer Ring
                sage: g = CyclicPermutationGroup(6).an_element(); g                     # needs sage.groups sage.modules
                (1,2,3,4,5,6)
                sage: A.counit_on_basis(g)                                              # needs sage.groups sage.modules
                1
            """
            return self.base_ring().one()

        def counit(self,x):
            r"""
            Return the counit of the element ``x`` of the group
            algebra.

            This is the sum of all coefficients of ``x`` with respect
            to the standard basis of the group algebra.

            EXAMPLES::

                sage: A = CyclicPermutationGroup(6).algebra(ZZ); A                      # needs sage.groups sage.modules
                Algebra of
                 Cyclic group of order 6 as a permutation group over Integer Ring
                sage: a = A.an_element(); a                                             # needs sage.groups sage.modules
                () + 3*(1,2,3,4,5,6) + 3*(1,3,5)(2,4,6)
                sage: a.counit()                                                        # needs sage.groups sage.modules
                7
            """
            return self.base_ring().sum(x.coefficients())

        def is_integral_domain(self, proof=True):
            r"""
            Return ``True`` if ``self`` is an integral domain.

            This is false unless ``self.base_ring()`` is an integral
            domain, and even then it is false unless ``self.group()``
            has no nontrivial elements of finite order. I don't know
            if this condition suffices, but it obviously does if the
            group is abelian and finitely generated.

            EXAMPLES::

                sage: # needs sage.groups sage.modules
                sage: S2 = SymmetricGroup(2)
                sage: GroupAlgebra(S2).is_integral_domain()                             # needs sage.combinat
                False
                sage: S1 = SymmetricGroup(1)
                sage: GroupAlgebra(S1).is_integral_domain()                             # needs sage.combinat
                True
                sage: GroupAlgebra(S1, IntegerModRing(4)).is_integral_domain()          # needs sage.combinat
                False
                sage: GroupAlgebra(AbelianGroup(1)).is_integral_domain()
                True
                sage: GroupAlgebra(AbelianGroup(2, [0,2])).is_integral_domain()
                False
                sage: GroupAlgebra(GL(2, ZZ)).is_integral_domain()      # not implemented
                False
            """
            from sage.sets.set import Set
            ans = False
            try:
                if self.base_ring().is_integral_domain():
                    if self.group().is_finite():
                        if self.group().order() > 1:
                            ans = False
                        else:
                            ans = True
                    else:
                        if self.group().is_abelian():
                            invs = self.group().invariants()
                            if Set(invs) != Set([0]):
                                ans = False
                            else:
                                ans = True
                        else:
                            raise NotImplementedError
                else:
                    ans = False
            except (AttributeError, NotImplementedError):
                if proof:
                    raise NotImplementedError("cannot determine whether self is an integral domain")

            return ans

        # I haven't written is_noetherian(), because I don't know when group
        # algebras are Noetherian, and I haven't written is_prime_field(), because
        # I don't know if that means "is canonically isomorphic to a prime field"
        # or "is identical to a prime field".

    class ElementMethods:

        def central_form(self):
            r"""
            Return ``self`` expressed in the canonical basis of the center
            of the group algebra.

            INPUT:

            - ``self`` -- an element of the center of the group algebra

            OUTPUT:

            - A formal linear combination of the conjugacy class
              representatives representing its coordinates in the
              canonical basis of the center. See
              :meth:`Groups.Algebras.ParentMethods.center_basis` for
              details.

            .. WARNING::

                - This method requires the underlying group to
                  have a method ``conjugacy_classes_representatives``
                  (every permutation group has one, thanks GAP!).
                - This method does not check that the element is
                  indeed central. Use the method
                  :meth:`Monoids.Algebras.ElementMethods.is_central`
                  for this purpose.
                - This function has a complexity linear in the
                  number of conjugacy classes of the group. One
                  could easily implement a function whose
                  complexity is linear in the size of the support
                  of ``self``.

            EXAMPLES::

                sage: # needs sage.combinat sage.groups sage.modules
                sage: QS3 = SymmetricGroup(3).algebra(QQ)
                sage: A = QS3([2,3,1]) + QS3([3,1,2])
                sage: A.central_form()
                B[(1,2,3)]
                sage: QS4 = SymmetricGroup(4).algebra(QQ)
                sage: B = sum(len(s.cycle_type()) * QS4(s) for s in Permutations(4))
                sage: B.central_form()
                4*B[()] + 3*B[(1,2)] + 2*B[(1,2)(3,4)] + 2*B[(1,2,3)] + B[(1,2,3,4)]

            The following test fails due to a bug involving combinatorial free modules and
            the coercion system (see :issue:`28544`)::

                sage: # needs sage.groups sage.modules
                sage: G = PermutationGroup([[(1,2,3),(4,5)], [(3,4)]])
                sage: QG = GroupAlgebras(QQ).example(G)
                sage: s = sum(QG.basis())
                sage: s.central_form()          # not tested
                B[()] + B[(4,5)] + B[(3,4,5)] + B[(2,3)(4,5)]
                + B[(2,3,4,5)] + B[(1,2)(3,4,5)] + B[(1,2,3,4,5)]

            .. SEEALSO::

                - :meth:`Groups.Algebras.ParentMethods.center_basis`
                - :meth:`Monoids.Algebras.ElementMethods.is_central`
            """
            from sage.combinat.free_module import CombinatorialFreeModule
            conj_classes_reps = self.parent().basis().keys().conjugacy_classes_representatives()
            Z = CombinatorialFreeModule(self.base_ring(), conj_classes_reps)
            return sum(self[i] * Z.basis()[i] for i in Z.basis().keys())
