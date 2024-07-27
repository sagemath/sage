# sage_setup: distribution = sagemath-groups
r"""
Kernel Subgroups

The kernel of a homomorphism implemented as a subgroup.

AUTHORS:

- Travis Scrimshaw (1-2023): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2023 Travis Scrimshaw <tcscrims at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.groups import Groups
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class KernelSubgroup(UniqueRepresentation, Parent):
    r"""
    The kernel (normal) subgroup.

    Let `\phi : G \to H` be a group homomorphism. The kernel
    `K = \{\phi(g) = 1 | g \in G\}` is a normal subgroup of `G`.
    """
    def __init__(self, morphism):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: S2 = SymmetricGroup(2)
            sage: S3 = SymmetricGroup(3)
            sage: H = Hom(S3, S2)
            sage: phi = H(S2.__call__)
            sage: from sage.groups.kernel_subgroup import KernelSubgroup
            sage: K = KernelSubgroup(phi)
            sage: TestSuite(K).run()
        """
        self._morphism = morphism
        cat = Groups().Subobjects()
        base_cat = morphism.domain().category()
        if base_cat in Groups().Finite():
            cat = cat.Finite()
        elif base_cat in Groups().Enumerated():
            cat = cat.Enumerated()
        Parent.__init__(self, category=cat)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: S2 = SymmetricGroup(2)
            sage: S3 = SymmetricGroup(3)
            sage: H = Hom(S3, S2)
            sage: phi = H(S2.__call__)
            sage: from sage.groups.kernel_subgroup import KernelSubgroup
            sage: KernelSubgroup(phi)
            Kernel subgroup defined by Generic morphism:
              From: Symmetric group of order 3! as a permutation group
              To:   Symmetric group of order 2! as a permutation group
        """
        return "Kernel subgroup defined by {}".format(self._morphism)

    def gens(self) -> tuple:
        r"""
        Return the generators of ``self``.

        EXAMPLES::

            sage: S2 = SymmetricGroup(2)
            sage: S3 = SymmetricGroup(3)
            sage: H = Hom(S3, S2)
            sage: phi = H(S2.__call__)
            sage: from sage.groups.kernel_subgroup import KernelSubgroup
            sage: K = KernelSubgroup(phi)
            sage: K.gens()
            ((),)
        """
        if self.ambient() in Groups().Finite():
            return tuple(self)
        raise NotImplementedError("only implemented for finite groups")

    def defining_morphism(self):
        r"""
        Return the defining morphism of ``self``.

        EXAMPLES::

            sage: PJ3 = groups.misc.PureCactus(3)                                       # needs sage.rings.number_field
            sage: PJ3.defining_morphism()                                               # needs sage.rings.number_field
            Conversion via _from_cactus_group_element map:
              From: Cactus Group with 3 fruit
              To:   Symmetric group of order 3! as a permutation group
        """
        return self._morphism

    @cached_method
    def ambient(self):
        r"""
        Return the ambient group of ``self``.

        EXAMPLES::

            sage: PJ3 = groups.misc.PureCactus(3)                                       # needs sage.rings.number_field
            sage: PJ3.ambient()                                                         # needs sage.rings.number_field
            Cactus Group with 3 fruit
        """
        return self._morphism.domain()

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: PJ3 = groups.misc.PureCactus(3)                                       # needs sage.rings.number_field
            sage: PJ3.an_element()                                                      # needs sage.rings.number_field
            1
        """
        return self.element_class(self, self.ambient().one())

    def lift(self, x):
        r"""
        Lift ``x`` to the ambient group of ``self``.

        EXAMPLES::

            sage: PJ3 = groups.misc.PureCactus(3)                                       # needs sage.rings.number_field
            sage: PJ3.lift(PJ3.an_element()).parent()                                   # needs sage.rings.number_field
            Cactus Group with 3 fruit
        """
        return x.value

    def retract(self, x):
        r"""
        Convert ``x`` to an element of ``self``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: J3 = groups.misc.Cactus(3)
            sage: s12,s13,s23 = J3.group_generators()
            sage: PJ3 = groups.misc.PureCactus(3)
            sage: elt = PJ3.retract(s23*s12*s23*s13); elt
            s[2,3]*s[1,2]*s[2,3]*s[1,3]
            sage: elt.parent() is PJ3
            True
        """
        return self._element_constructor_(x)

    def _element_constructor_(self, x):
        r"""
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: J3 = groups.misc.Cactus(3)
            sage: s12,s13,s23 = J3.group_generators()
            sage: PJ3 = groups.misc.PureCactus(3)
            sage: elt = PJ3(s23*s12*s23*s13)
            sage: elt.parent() is PJ3
            True
        """
        if self._morphism(x) != self._morphism.codomain().one():
            raise ValueError("{} is not in the kernel of {}".format(x, self._morphism))
        return self.element_class(self, x)

    def __iter__(self):
        r"""
        Iterate through ``self``.

        EXAMPLES::

            sage: S2 = SymmetricGroup(2)
            sage: S3 = SymmetricGroup(3)
            sage: H = Hom(S3, S2)
            sage: phi = H(S2.__call__)
            sage: from sage.groups.kernel_subgroup import KernelSubgroup
            sage: K = KernelSubgroup(phi)
            sage: list(K)
            [()]
        """
        for g in self.ambient():
            try:
                yield self(g)
            except ValueError:
                pass

    class Element(ElementWrapper):
        def _mul_(self, other):
            r"""
            Multiply ``self`` and ``other``.

            EXAMPLES::

                sage: # needs sage.rings.number_field
                sage: J3 = groups.misc.Cactus(3)
                sage: s12,s13,s23 = J3.group_generators()
                sage: PJ3 = groups.misc.PureCactus(3)
                sage: elt = PJ3(s23*s12*s23*s13)
                sage: elt * elt
                s[2,3]*s[1,2]*s[2,3]*s[1,2]*s[2,3]*s[1,2]
            """
            return type(self)(self.parent(), self.value * other.value)

        def __invert__(self):
            r"""
            Return the inverse of ``self``.

            EXAMPLES::

                sage: # needs sage.rings.number_field
                sage: J3 = groups.misc.Cactus(3)
                sage: s12,s13,s23 = J3.group_generators()
                sage: PJ3 = groups.misc.PureCactus(3)
                sage: elt = PJ3(s23*s12*s23*s13)
                sage: ~elt
                s[1,2]*s[2,3]*s[1,2]*s[1,3]
            """
            return type(self)(self.parent(), ~self.value)
