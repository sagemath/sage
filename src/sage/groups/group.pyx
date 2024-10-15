"""
Base class for groups
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
from sage.misc.superseded import deprecation
from sage.rings.infinity import infinity
from sage.structure.parent cimport Parent


def is_Group(x):
    """
    Return whether ``x`` is a group object.

    INPUT:

    - ``x`` -- anything

    OUTPUT: boolean

    EXAMPLES::

        sage: from sage.groups.group import is_Group
        sage: is_Group("a string")
        doctest:warning...DeprecationWarning: use instead G in Groups()
        See https://github.com/sagemath/sage/issues/37449 for details.
        False
        sage: F.<a,b> = FreeGroup()                                                     # needs sage.groups
        sage: is_Group(F)                                                               # needs sage.groups
        True
    """
    deprecation(37449, 'use instead G in Groups()')
    return isinstance(x, Group)


cdef class Group(Parent):
    """
    Base class for all groups.

    TESTS::

        sage: from sage.groups.group import Group
        sage: G = Group()
        sage: TestSuite(G).run(skip = ["_test_an_element",\
                                       "_test_associativity",\
                                       "_test_elements",\
                                       "_test_elements_eq_reflexive",\
                                       "_test_elements_eq_symmetric",\
                                       "_test_elements_eq_transitive",\
                                       "_test_elements_neq",\
                                       "_test_inverse",\
                                       "_test_one",\
                                       "_test_pickling",\
                                       "_test_prod",\
                                       "_test_some_elements"])

    Generic groups have very little functionality::

        sage: 4 in G
        Traceback (most recent call last):
        ...
        NotImplementedError: cannot construct elements of <sage.groups.group.Group object at ...>
    """
    def __init__(self, base=None, category=None):
        """
        The Python constructor.

        TESTS::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.category()
            Category of groups
            sage: G = Group(category=Groups()) # todo: do the same test with some subcategory of Groups when there will exist one
            sage: G.category()
            Category of groups
            sage: G = Group(category=CommutativeAdditiveGroups())
            Traceback (most recent call last):
            ...
            ValueError: (Category of commutative additive groups,) is not a subcategory of Category of groups
            sage: G._repr_option('element_is_atomic')
            False

        Check for :issue:`8119`::

            sage: # needs sage.groups
            sage: G = SymmetricGroup(2)
            sage: h = hash(G)
            sage: G.rename('S2')
            sage: h == hash(G)
            True
        """
        from sage.categories.groups import Groups
        if category is None:
            category = Groups()
        else:
            if not isinstance(category, tuple):
                category = (category,)
            if not any(cat.is_subcategory(Groups()) for cat in category):
                raise ValueError("%s is not a subcategory of %s"%(category, Groups()))
        Parent.__init__(self, base=base, category=category)

    def is_abelian(self):
        """
        Test whether this group is abelian.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.is_abelian()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_commutative(self):
        r"""
        Test whether this group is commutative.

        This is an alias for is_abelian, largely to make groups work
        well with the Factorization class.

        (Note for developers: Derived classes should override is_abelian, not
        is_commutative.)

        EXAMPLES::

            sage: SL(2, 7).is_commutative()                                             # needs sage.libs.gap sage.modules sage.rings.finite_rings
            False
        """
        return self.is_abelian()

    def order(self):
        """
        Return the number of elements of this group.

        This is either a positive integer or infinity.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.order()
            Traceback (most recent call last):
            ...
            NotImplementedError

        TESTS::

            sage: H = SL(2, QQ)                                                         # needs sage.modules
            sage: H.order()                                                             # needs sage.modules
            +Infinity
        """
        try:
            return self.cardinality()
        except AttributeError:
            raise NotImplementedError

    def is_finite(self):
        """
        Return ``True`` if this group is finite.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.is_finite()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return self.order() != infinity

    def is_trivial(self):
        r"""
        Return ``True`` if this group is the trivial group.

        A group is trivial, if it consists only of the identity
        element.

        .. WARNING::

            It is in principle undecidable whether a group is
            trivial, for example, if the group is given by a finite
            presentation.  Thus, this method may not terminate.

        EXAMPLES::

            sage: groups.presentation.Cyclic(1).is_trivial()                            # needs sage.groups
            True

            sage: # needs sage.groups
            sage: G.<a,b> = FreeGroup('a, b')
            sage: H = G / (a^2, b^3, a*b*~a*~b)
            sage: H.is_trivial()
            False

        A non-trivial presentation of the trivial group::

            sage: # needs sage.groups
            sage: F.<a,b> = FreeGroup()
            sage: J = F / ((~a)*b*a*(~b)^2, (~b)*a*b*(~a)^2)
            sage: J.is_trivial()
            True
        """
        return self.order() == 1

    def is_multiplicative(self):
        r"""
        Return ``True`` if the group operation is given by ``*`` (rather than ``+``).

        Override for additive groups.

        EXAMPLES::

            sage: from sage.groups.group import Group
            sage: G = Group()
            sage: G.is_multiplicative()
            True
        """
        return True

    def _an_element_(self):
        """
        Return an element.

        OUTPUT: an element of the group

        EXAMPLES::

            sage: G = AbelianGroup([2,3,4,5])                                           # needs sage.modules
            sage: G.an_element()                                                        # needs sage.modules
            f0*f1*f2*f3
        """
        return self.prod(self.gens())

    def quotient(self, H, **kwds):
        """
        Return the quotient of this group by the normal subgroup `H`.

        EXAMPLES::

            sage: from sage.groups.group import Group
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

            sage: from sage.groups.group import AbelianGroup
            sage: G = AbelianGroup()
            sage: G.is_abelian()
            True
        """
        return True

cdef class FiniteGroup(Group):
    """
    Generic finite group.
    """

    def __init__(self, base=None, category=None):
        """
        The Python constructor.

        TESTS::

            sage: from sage.groups.group import FiniteGroup
            sage: G = FiniteGroup()
            sage: G.category()
            Category of finite groups
        """
        from sage.categories.finite_groups import FiniteGroups
        if category is None:
            category = FiniteGroups()
        else:
            if not isinstance(category, tuple):
                category = (category,)
            if not any(cat.is_subcategory(FiniteGroups()) for cat in category):
                raise ValueError("%s is not a subcategory of %s" % (category, FiniteGroups()))
        Parent.__init__(self, base=base, category=category)

    def is_finite(self):
        """
        Return ``True``.

        EXAMPLES::

            sage: from sage.groups.group import FiniteGroup
            sage: G = FiniteGroup()
            sage: G.is_finite()
            True
        """
        return True
