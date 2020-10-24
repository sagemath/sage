"""
Indexed Free Groups

Free groups and free abelian groups implemented using an indexed set of
generators.

AUTHORS:

- Travis Scrimshaw (2013-10-16): Initial version
"""

##############################################################################
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
##############################################################################

from sage.categories.groups import Groups
from sage.categories.poor_man_map import PoorManMap
from sage.groups.group import Group, AbelianGroup
from sage.monoids.indexed_free_monoid import (IndexedMonoid,
        IndexedFreeMonoidElement, IndexedFreeAbelianMonoidElement)
from sage.misc.cachefunc import cached_method
import sage.data_structures.blas_dict as blas
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.sets.family import Family


class IndexedGroup(IndexedMonoid):
    """
    Base class for free (abelian) groups whose generators are indexed
    by a set.

    TESTS:

    We check finite properties::

        sage: G = Groups().free(index_set=ZZ)
        sage: G.is_finite()
        False
        sage: G = Groups().free(index_set='abc')
        sage: G.is_finite()
        False
        sage: G = Groups().free(index_set=[])
        sage: G.is_finite()
        True

    ::

        sage: G = Groups().Commutative().free(index_set=ZZ)
        sage: G.is_finite()
        False
        sage: G = Groups().Commutative().free(index_set='abc')
        sage: G.is_finite()
        False
        sage: G = Groups().Commutative().free(index_set=[])
        sage: G.is_finite()
        True
    """
    def order(self):
        r"""
        Return the number of elements of ``self``, which is `\infty` unless
        this is the trivial group.

        EXAMPLES::

            sage: G = Groups().free(index_set=ZZ)
            sage: G.order()
            +Infinity
            sage: G = Groups().Commutative().free(index_set='abc')
            sage: G.order()
            +Infinity
            sage: G = Groups().Commutative().free(index_set=[])
            sage: G.order()
            1
        """
        return self.cardinality()

    def rank(self):
        """
        Return the rank of ``self``.

        This is the number of generators of ``self``.

        EXAMPLES::

            sage: G = Groups().free(index_set=ZZ)
            sage: G.rank()
            +Infinity
            sage: G = Groups().free(index_set='abc')
            sage: G.rank()
            3
            sage: G = Groups().free(index_set=[])
            sage: G.rank()
            0

        ::

            sage: G = Groups().Commutative().free(index_set=ZZ)
            sage: G.rank()
            +Infinity
            sage: G = Groups().Commutative().free(index_set='abc')
            sage: G.rank()
            3
            sage: G = Groups().Commutative().free(index_set=[])
            sage: G.rank()
            0
        """
        return self.group_generators().cardinality()

    @cached_method
    def group_generators(self):
        """
        Return the group generators of ``self``.

        EXAMPLES::

            sage: G = Groups.free(index_set=ZZ)
            sage: G.group_generators()
            Lazy family (Generator map from Integer Ring to
             Free group indexed by Integer Ring(i))_{i in Integer Ring}
            sage: G = Groups().free(index_set='abcde')
            sage: sorted(G.group_generators())
            [F['a'], F['b'], F['c'], F['d'], F['e']]
        """
        if self._indices.cardinality() == infinity:
            gen = PoorManMap(self.gen, domain=self._indices, codomain=self, name="Generator map")
            return Family(self._indices, gen)
        return Family(self._indices, self.gen)

    gens = group_generators

class IndexedFreeGroup(IndexedGroup, Group):
    """
    An indexed free group.

    EXAMPLES::

        sage: G = Groups().free(index_set=ZZ)
        sage: G
        Free group indexed by Integer Ring
        sage: G = Groups().free(index_set='abcde')
        sage: G
        Free group indexed by {'a', 'b', 'c', 'd', 'e'}
    """
    def __init__(self, indices, prefix, category=None, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: G = Groups().free(index_set=ZZ)
            sage: TestSuite(G).run()
            sage: G = Groups().free(index_set='abc')
            sage: TestSuite(G).run()
        """
        category = Groups().or_subcategory(category)
        IndexedGroup.__init__(self, indices, prefix, category, **kwds)

    def _repr_(self):
        """
        Return a string representation of ``self``

        TESTS::

            sage: Groups().free(index_set=ZZ)       # indirect doctest
            Free group indexed by Integer Ring
        """
        return 'Free group indexed by {}'.format(self._indices)

    @cached_method
    def one(self):
        """
        Return the identity element of ``self``.

        EXAMPLES::

            sage: G = Groups().free(ZZ)
            sage: G.one()
            1
        """
        return self.element_class(self, ())

    def gen(self, x):
        """
        The generator indexed by ``x`` of ``self``.

        EXAMPLES::

            sage: G = Groups().free(index_set=ZZ)
            sage: G.gen(0)
            F[0]
            sage: G.gen(2)
            F[2]
        """
        if x not in self._indices:
            raise IndexError("{} is not in the index set".format(x))
        try:
            return self.element_class(self, ((self._indices(x),1),))
        except TypeError: # Backup (if it is a string)
            return self.element_class(self, ((x,1),))

    class Element(IndexedFreeMonoidElement):
        def __len__(self):
            """
            Return the length of ``self``.

            EXAMPLES::

                sage: G = Groups().free(index_set=ZZ)
                sage: a,b,c,d,e = [G.gen(i) for i in range(5)]
                sage: elt = a*c^-3*b^-2*a
                sage: elt.length()
                7
                sage: len(elt)
                7

                sage: G = Groups().free(index_set=ZZ)
                sage: a,b,c,d,e = [G.gen(i) for i in range(5)]
                sage: elt = a*c^-3*b^-2*a
                sage: elt.length()
                7
                sage: len(elt)
                7
            """
            return sum(abs(exp) for gen,exp in self._monomial)

        length = __len__

        def _mul_(self, other):
            """
            Multiply ``self`` by ``other``.

            EXAMPLES::

                sage: G = Groups().free(index_set=ZZ)
                sage: a,b,c,d,e = [G.gen(i) for i in range(5)]
                sage: a*b^2*e*d
                F[0]*F[1]^2*F[4]*F[3]
                sage: (a*b^2*d^2) * (d^-4*b*e)
                F[0]*F[1]^2*F[3]^-2*F[1]*F[4]
                sage: (a*b^-2*d^2) * (d^-2*b^2*a^-1)
                1
            """
            if not self._monomial:
                return other
            if not other._monomial:
                return self

            ret = list(self._monomial)
            rhs = list(other._monomial)
            while ret and rhs and ret[-1][0] == rhs[0][0]:
                rhs[0] = (rhs[0][0], rhs[0][1] + ret.pop()[1])
                if rhs[0][1] == 0:
                    rhs.pop(0)
            ret += rhs
            return self.__class__(self.parent(), tuple(ret))

        def __invert__(self):
            """
            Return the inverse of ``self``.

            EXAMPLES::

                sage: G = Groups().free(index_set=ZZ)
                sage: a,b,c,d,e = [G.gen(i) for i in range(5)]
                sage: x = a*b^2*e^-1*d; ~x
                F[3]^-1*F[4]*F[1]^-2*F[0]^-1
                sage: x * ~x
                1
            """
            return self.__class__(self.parent(),
                   tuple((x[0], -x[1]) for x in reversed(self._monomial)))

        def to_word_list(self):
            """
            Return ``self`` as a word represented as a list whose entries
            are the pairs ``(i, s)`` where ``i`` is the index and ``s`` is
            the sign.

            EXAMPLES::

                sage: G = Groups().free(index_set=ZZ)
                sage: a,b,c,d,e = [G.gen(i) for i in range(5)]
                sage: x = a*b^2*e*a^-1
                sage: x.to_word_list()
                [(0, 1), (1, 1), (1, 1), (4, 1), (0, -1)]
            """
            sign = lambda x: 1 if x > 0 else -1 # It is never 0
            return [ (k, sign(e)) for k,e in self._sorted_items()
                     for dummy in range(abs(e))]

class IndexedFreeAbelianGroup(IndexedGroup, AbelianGroup):
    """
    An indexed free abelian group.

    EXAMPLES::

        sage: G = Groups().Commutative().free(index_set=ZZ)
        sage: G
        Free abelian group indexed by Integer Ring
        sage: G = Groups().Commutative().free(index_set='abcde')
        sage: G
        Free abelian group indexed by {'a', 'b', 'c', 'd', 'e'}
    """
    def __init__(self, indices, prefix, category=None, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: G = Groups().Commutative().free(index_set=ZZ)
            sage: TestSuite(G).run()
            sage: G = Groups().Commutative().free(index_set='abc')
            sage: TestSuite(G).run()
        """
        category = Groups().or_subcategory(category)
        IndexedGroup.__init__(self, indices, prefix, category, **kwds)

    def _repr_(self):
        """
        TESTS::

            sage: Groups.Commutative().free(index_set=ZZ)
            Free abelian group indexed by Integer Ring
        """
        return 'Free abelian group indexed by {}'.format(self._indices)

    def _element_constructor_(self, x=None):
        """
        Create an element of ``self`` from ``x``.

        EXAMPLES::

            sage: G = Groups().Commutative().free(index_set=ZZ)
            sage: G(G.gen(2))
            F[2]
            sage: G([[1, 3], [-2, 12]])
            F[-2]^12*F[1]^3
            sage: G({1: 3, -2: 12})
            F[-2]^12*F[1]^3
            sage: G(-5)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert -5, use gen() instead

        TESTS::

            sage: G([(1, 3), (1, -5)])
            F[1]^-2

            sage: G([(42, 0)])
            1
            sage: G([(42, 3), (42, -3)])
            1
            sage: G({42: 0})
            1
        """
        if isinstance(x, (list, tuple)):
            d = dict()
            for k, v in x:
                if k in d:
                    d[k] += v
                else:
                    d[k] = v
            x = d
        if isinstance(x, dict):
            x = {k: v for k, v in x.items() if v != 0}
        return IndexedGroup._element_constructor_(self, x)

    @cached_method
    def one(self):
        """
        Return the identity element of ``self``.

        EXAMPLES::

            sage: G = Groups().Commutative().free(index_set=ZZ)
            sage: G.one()
            1
        """
        return self.element_class(self, {})

    def gen(self, x):
        """
        The generator indexed by ``x`` of ``self``.

        EXAMPLES::

            sage: G = Groups().Commutative().free(index_set=ZZ)
            sage: G.gen(0)
            F[0]
            sage: G.gen(2)
            F[2]
        """
        if x not in self._indices:
            raise IndexError("{} is not in the index set".format(x))
        try:
            return self.element_class(self, {self._indices(x):1})
        except TypeError: # Backup (if it is a string)
            return self.element_class(self, {x:1})

    class Element(IndexedFreeAbelianMonoidElement, IndexedFreeGroup.Element):
        def _mul_(self, other):
            """
            Multiply ``self`` by ``other``.

            EXAMPLES::

                sage: G = Groups().Commutative().free(index_set=ZZ)
                sage: a,b,c,d,e = [G.gen(i) for i in range(5)]
                sage: a*b^2*e^-1*d
                F[0]*F[1]^2*F[3]*F[4]^-1
                sage: (a*b^2*d^2) * (d^-4*b^-2*e)
                F[0]*F[3]^-2*F[4]
                sage: (a*b^-2*d^2) * (d^-2*b^2*a^-1)
                1
            """
            return self.__class__(self.parent(),
                                  blas.add(self._monomial, other._monomial))

        def __invert__(self):
            """
            Return the inverse of ``self``.

            EXAMPLES::

                sage: G = Groups().Commutative().free(index_set=ZZ)
                sage: a,b,c,d,e = [G.gen(i) for i in range(5)]
                sage: x = a*b^2*e^-1*d; ~x
                F[0]^-1*F[1]^-2*F[3]^-1*F[4]
                sage: x * ~x
                1
            """
            return self ** -1

        def __floordiv__(self, a):
            """
            Return the division of ``self`` by ``a``.

            EXAMPLES::

                sage: G = Groups().Commutative().free(index_set=ZZ)
                sage: a,b,c,d,e = [G.gen(i) for i in range(5)]
                sage: elt = a*b*c^3*d^2; elt
                F[0]*F[1]*F[2]^3*F[3]^2
                sage: elt // a
                F[1]*F[2]^3*F[3]^2
                sage: elt // c
                F[0]*F[1]*F[2]^2*F[3]^2
                sage: elt // (a*b*d^2)
                F[2]^3
                sage: elt // a^4
                F[0]^-3*F[1]*F[2]^3*F[3]^2
            """
            return self * ~a

        def __pow__(self, n):
            """
            Raise ``self`` to the power of ``n``.

            EXAMPLES::

                sage: G = Groups().Commutative().free(index_set=ZZ)
                sage: a,b,c,d,e = [G.gen(i) for i in range(5)]
                sage: x = a*b^2*e^-1*d; x
                F[0]*F[1]^2*F[3]*F[4]^-1
                sage: x^3
                F[0]^3*F[1]^6*F[3]^3*F[4]^-3
                sage: x^0
                1
                sage: x^-3
                F[0]^-3*F[1]^-6*F[3]^-3*F[4]^3
            """
            if not isinstance(n, (int, Integer)):
                raise TypeError("Argument n (= {}) must be an integer".format(n))
            if n == 1:
                return self
            if n == 0:
                return self.parent().one()
            return self.__class__(self.parent(), {k:v*n for k,v in self._monomial.items()})

