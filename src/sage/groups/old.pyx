# sage_setup: distribution = sagemath-objects
"""
Deprecated base class for groups
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

doc = """
Deprecated base class for all groups
"""
import sage.rings.integer_ring
from sage.misc.superseded import deprecation


cdef class Group(sage.structure.parent.Parent):
    """
    Generic group class.
    """
    def __init__(self, category=None):
        """

        TESTS::

            sage: from sage.groups.old import Group
            sage: G = Group()
            doctest:warning...:
            DeprecationWarning: do not use the old Group class
            See https://github.com/sagemath/sage/issues/37449 for details.
            sage: G.category()
            Category of groups
            sage: G = Group(category = Groups()) # todo: do the same test with some subcategory of Groups when there will exist one
            sage: G.category()
            Category of groups
            sage: G = Group(category = CommutativeAdditiveGroups())
            Traceback (most recent call last):
            ...
            AssertionError: Category of commutative additive groups is not a subcategory of Category of groups

         Check for :issue:`8119`::

            sage: # needs sage.groups
            sage: G = SymmetricGroup(2)
            sage: h = hash(G)
            sage: G.rename('S2')
            sage: h == hash(G)
            True
        """
        deprecation(37449, 'do not use the old Group class')
        from sage.categories.basic import Groups
        if category is None:
            category = Groups()
        else:
            assert category.is_subcategory(Groups()), "%s is not a subcategory of %s" % (category, Groups())

        sage.structure.parent.Parent.__init__(self,
                base=sage.rings.integer_ring.ZZ, category=category)

    def __contains__(self, x):
        r"""
        Return ``True`` if coercion of ``x`` into ``self`` is defined.

        EXAMPLES::

            sage: from sage.groups.old import Group
            sage: G = Group()
            sage: 4 in G               #indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot construct elements of <sage.groups.old.Group object at ...>
        """
        try:
            self(x)
        except TypeError:
            return False
        return True

    def is_abelian(self):
        """
        Return ``True`` if this group is abelian.

        EXAMPLES::

            sage: from sage.groups.old import Group
            sage: G = Group()
            sage: G.is_abelian()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_commutative(self):
        r"""
        Return ``True`` if this group is commutative. This is an alias for
        is_abelian, largely to make groups work well with the Factorization
        class.

        (Note for developers: Derived classes should override is_abelian, not
        is_commutative.)

        EXAMPLES::

            sage: SL(2, 7).is_commutative()                                             # needs sage.libs.gap sage.modules
            False
        """
        return self.is_abelian()

    def order(self):
        """
        Return the number of elements of this group, which is either a
        positive integer or infinity.

        EXAMPLES::

            sage: from sage.groups.old import Group
            sage: G = Group()
            sage: G.order()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_multiplicative(self):
        r"""
        Return ``True`` if the group operation is given by \* (rather than +).

        Override for additive groups.

        EXAMPLES::

            sage: from sage.groups.old import Group
            sage: G = Group()
            sage: G.is_multiplicative()
            True
        """
        return True

    def random_element(self, bound=None):
        """
        Return a random element of this group.

        EXAMPLES::

            sage: from sage.groups.old import Group
            sage: G = Group()
            sage: G.random_element()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def quotient(self, H, **kwds):
        """
        Return the quotient of this group by the normal subgroup `H`.

        EXAMPLES::

            sage: from sage.groups.old import Group
            sage: G = Group()
            sage: G.quotient(G)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

cdef class AbelianGroup(Group):
    """
    Generic abelian group.
    """
    def is_abelian(self):
        """
        Return True.

        EXAMPLES::

            sage: from sage.groups.old import AbelianGroup
            sage: G = AbelianGroup()
            sage: G.is_abelian()
            True
        """
        return True

cdef class FiniteGroup(Group):
    """
    Generic finite group.
    """
    def is_finite(self):
        """
        Return True.

        EXAMPLES::

            sage: from sage.groups.old import FiniteGroup
            sage: G = FiniteGroup()
            sage: G.is_finite()
            True
        """
        return True

cdef class AlgebraicGroup(Group):
    """
    Generic algebraic group.
    """
