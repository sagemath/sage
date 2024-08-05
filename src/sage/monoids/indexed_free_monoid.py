"""
Indexed Monoids

AUTHORS:

- Travis Scrimshaw (2013-10-15)
"""

# ****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import MonoidElement
from sage.structure.indexed_generators import IndexedGenerators, parse_indices_names
from sage.structure.richcmp import op_EQ, op_NE, richcmp, rich_to_bool
import sage.data_structures.blas_dict as blas

from sage.categories.monoids import Monoids
from sage.categories.poor_man_map import PoorManMap
from sage.categories.sets_cat import Sets
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.sets.family import Family


class IndexedMonoidElement(MonoidElement):
    """
    An element of an indexed monoid.

    This is an abstract class which uses the (abstract) method
    :meth:`_sorted_items` for all of its functions. So to implement an
    element of an indexed monoid, one just needs to implement
    :meth:`_sorted_items`, which returns a list of pairs ``(i, p)`` where
    ``i`` is the index and ``p`` is the corresponding power, sorted in some
    order. For example, in the free monoid there is no such choice, but for
    the free abelian monoid, one could want lex order or have the highest
    powers first.

    Indexed monoid elements are ordered lexicographically with respect to
    the result of :meth:`_sorted_items` (which for abelian free monoids is
    influenced by the order on the indexing set).
    """
    def __init__(self, F, x):
        """
        Create the element ``x`` of an indexed free abelian monoid ``F``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: F.gen(1)
            F[1]
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: x = a^2 * b^3 * a^2 * b^4; x
            F[0]^4*F[1]^7
            sage: TestSuite(x).run()

            sage: F = FreeMonoid(index_set=tuple('abcde'))
            sage: a,b,c,d,e = F.gens()
            sage: a in F
            True
            sage: a*b in F
            True
            sage: TestSuite(a*d^2*e*c*a).run()
        """
        MonoidElement.__init__(self, F)
        self._monomial = x

    @abstract_method
    def _sorted_items(self):
        """
        Return the sorted items (i.e factors) of ``self``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: x = a*b^2*e*d
            sage: x._sorted_items()
            ((0, 1), (1, 2), (4, 1), (3, 1))

        .. SEEALSO::

            :meth:`_repr_`, :meth:`_latex_`, :meth:`print_options`
        """

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: a*b^2*e*d
            F[0]*F[1]^2*F[3]*F[4]
        """
        if not self._monomial:
            return '1'

        monomial = self._sorted_items()
        P = self.parent()

        scalar_mult = P._print_options['scalar_mult']

        exp = lambda v: f'^{v}' if v != 1 else ''
        return scalar_mult.join(P._repr_generator(g) + exp(v) for g,v in monomial)

    def _ascii_art_(self):
        r"""
        Return an ASCII art representation of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: ascii_art(a*e*d)
            F *F *F
             0  3  4
            sage: ascii_art(a*b^2*e*d)
                2
            F *F *F *F
             0  1  3  4
        """
        from sage.typeset.ascii_art import AsciiArt, ascii_art, empty_ascii_art

        if not self._monomial:
            return AsciiArt(["1"])

        monomial = self._sorted_items()
        P = self.parent()
        scalar_mult = P._print_options['scalar_mult']

        if all(x[1] == 1 for x in monomial):
            ascii_art_gen = lambda m: P._ascii_art_generator(m[0])
        else:
            pref = AsciiArt([P.prefix()])

            def ascii_art_gen(m):
                if m[1] != 1:
                    r = (AsciiArt([" " * len(pref)]) + ascii_art(m[1]))
                else:
                    r = empty_ascii_art
                r = r * P._ascii_art_generator(m[0])
                r._baseline = r._h - 2
                return r
        b = ascii_art_gen(monomial[0])
        for x in monomial[1:]:
            b = b + AsciiArt([scalar_mult]) + ascii_art_gen(x)
        return b

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: latex(a*b^2*e*d)
            F_{0} F_{1}^{2} F_{3} F_{4}
        """
        if not self._monomial:
            return '1'

        monomial = self._sorted_items()
        P = self.parent()

        scalar_mult = P._print_options['latex_scalar_mult']
        if scalar_mult is None:
            scalar_mult = P._print_options['scalar_mult']
            if scalar_mult == "*":
                scalar_mult = " "

        exp = lambda v: f'^{{{v}}}' if v != 1 else ''
        return scalar_mult.join(P._latex_generator(g) + exp(v) for g,v in monomial)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: list(b*a*c^3*b)
            [(F[1], 1), (F[0], 1), (F[2], 3), (F[1], 1)]

        ::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: list(b*c^3*a)
            [(F[0], 1), (F[1], 1), (F[2], 3)]
        """
        return ((self.parent().gen(index), exp) for (index,exp) in self._sorted_items())

    def _richcmp_(self, other, op):
        r"""
        Comparisons.

        TESTS::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: a == a
            True
            sage: a*e == a*e
            True
            sage: a*b*c^3*b*d == (a*b*c)*(c^2*b*d)
            True
            sage: a != b
            True
            sage: a*b != b*a
            True
            sage: a*b*c^3*b*d != (a*b*c)*(c^2*b*d)
            False

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: a < b
            True
            sage: a*b < b*a
            True
            sage: a*b < a*a
            False
            sage: a^2*b < a*b*b
            True
            sage: b > a
            True
            sage: a*b > b*a
            False
            sage: a*b > a*a
            True
            sage: a*b <= b*a
            True
            sage: a*b <= b*a
            True

            sage: FA = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [FA.gen(i) for i in range(5)]
            sage: a == a
            True
            sage: a*e == e*a
            True
            sage: a*b*c^3*b*d == a*d*(b^2*c^2)*c
            True
            sage: a != b
            True
            sage: a*b != a*a
            True
            sage: a*b*c^3*b*d != a*d*(b^2*c^2)*c
            False
        """
        if self._monomial == other._monomial:
            # Equal
            return rich_to_bool(op, 0)
        if op == op_EQ or op == op_NE:
            # Not equal
            return rich_to_bool(op, 1)
        return richcmp(self.to_word_list(), other.to_word_list(), op)

    def support(self):
        """
        Return a list of the objects indexing ``self`` with
        nonzero exponents.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: (b*a*c^3*b).support()
            [0, 1, 2]

        ::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: (a*c^3).support()
            [0, 2]
        """
        supp = {key for key, exp in self._sorted_items() if exp != 0}
        try:
            return sorted(supp, key=print_options['sorting_key'],
                          reverse=print_options['sorting_reverse'])
        except Exception:  # Sorting the output is a plus, but if we can't, no big deal
            return list(supp)

    def leading_support(self):
        """
        Return the support of the leading generator of ``self``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: (b*a*c^3*a).leading_support()
            1

        ::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: (b*c^3*a).leading_support()
            0
        """
        if not self:
            return None
        return self._sorted_items()[0][0]

    def trailing_support(self):
        """
        Return the support of the trailing generator of ``self``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: (b*a*c^3*a).trailing_support()
            0

        ::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: (b*c^3*a).trailing_support()
            2
        """
        if not self:
            return None
        return self._sorted_items()[-1][0]

    def to_word_list(self):
        """
        Return ``self`` as a word represented as a list whose entries
        are indices of ``self``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: (b*a*c^3*a).to_word_list()
            [1, 0, 2, 2, 2, 0]

        ::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: (b*c^3*a).to_word_list()
            [0, 1, 2, 2, 2]
        """
        return [k for k,e in self._sorted_items() for dummy in range(e)]

    def is_one(self) -> bool:
        """
        Return if ``self`` is the identity element.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: (b*a*c^3*a).is_one()
            False
            sage: F.one().is_one()
            True

        ::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: (b*c^3*a).is_one()
            False
            sage: F.one().is_one()
            True
        """
        return not self._monomial


class IndexedFreeMonoidElement(IndexedMonoidElement):
    """
    An element of an indexed free abelian monoid.
    """
    def __init__(self, F, x):
        """
        Create the element ``x`` of an indexed free abelian monoid ``F``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=tuple('abcde'))
            sage: x = F( [(1, 2), (0, 1), (3, 2), (0, 1)] )
            sage: y = F( ((1, 2), (0, 1), [3, 2], [0, 1]) )
            sage: z = F( reversed([(0, 1), (3, 2), (0, 1), (1, 2)]) )
            sage: x == y and y == z
            True
            sage: TestSuite(x).run()
        """
        IndexedMonoidElement.__init__(self, F, tuple(map(tuple, x)))

    def __hash__(self):
        r"""
        TESTS::

            sage: F = FreeMonoid(index_set=tuple('abcde'))
            sage: hash(F ([(1,2),(0,1)]) ) == hash(((1, 2), (0, 1)))
            True
            sage: hash(F ([(0,2),(1,1)]) ) == hash(((0, 2), (1, 1)))
            True
        """
        return hash(self._monomial)

    def _sorted_items(self):
        """
        Return the sorted items (i.e factors) of ``self``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: x = a*b^2*e*d
            sage: x._sorted_items()
            ((0, 1), (1, 2), (4, 1), (3, 1))
            sage: F.print_options(sorting_reverse=True)
            sage: x._sorted_items()
            ((0, 1), (1, 2), (4, 1), (3, 1))
            sage: F.print_options(sorting_reverse=False) # reset to original state

        .. SEEALSO::

            :meth:`_repr_`, :meth:`_latex_`, :meth:`print_options`
        """
        return self._monomial

    def _mul_(self, other):
        """
        Multiply ``self`` by ``other``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: a*b^2*e*d
            F[0]*F[1]^2*F[4]*F[3]
            sage: (a*b^2*d^2) * (d^4*b*e)
            F[0]*F[1]^2*F[3]^6*F[1]*F[4]
        """
        if not self._monomial:
            return other
        if not other._monomial:
            return self

        ret = list(self._monomial)
        rhs = list(other._monomial)
        if ret[-1][0] == rhs[0][0]:
            rhs[0] = (rhs[0][0], rhs[0][1] + ret.pop()[1])
        ret += rhs
        return self.__class__(self.parent(), tuple(ret))

    def __len__(self):
        """
        Return the length of ``self``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: elt = a*c^3*b^2*a
            sage: elt.length()
            7
            sage: len(elt)
            7
        """
        return sum(exp for gen,exp in self._monomial)

    length = __len__


class IndexedFreeAbelianMonoidElement(IndexedMonoidElement):
    """
    An element of an indexed free abelian monoid.
    """
    def __init__(self, F, x):
        """
        Create the element ``x`` of an indexed free abelian monoid ``F``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: x = F([(0, 1), (2, 2), (-1, 2)])
            sage: y = F({0:1, 2:2, -1:2})
            sage: z = F(reversed([(0, 1), (2, 2), (-1, 2)]))
            sage: x == y and y == z
            True
            sage: TestSuite(x).run()
        """
        IndexedMonoidElement.__init__(self, F, dict(x))

    def _sorted_items(self):
        """
        Return the sorted items (i.e factors) of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: x = a*b^2*e*d
            sage: x._sorted_items()
            [(0, 1), (1, 2), (3, 1), (4, 1)]
            sage: F.print_options(sorting_reverse=True)
            sage: x._sorted_items()
            [(4, 1), (3, 1), (1, 2), (0, 1)]
            sage: F.print_options(sorting_reverse=False) # reset to original state

        .. SEEALSO::

            :meth:`_repr_`, :meth:`_latex_`, :meth:`print_options`
        """
        print_options = self.parent().print_options()
        v = list(self._monomial.items())
        try:
            v.sort(key=print_options['sorting_key'],
                   reverse=print_options['sorting_reverse'])
        except Exception:  # Sorting the output is a plus, but if we can't, no big deal
            pass
        return v

    def __hash__(self):
        r"""
        TESTS::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: H1 = hash( F([(0,1), (2,2)]) )
            sage: H2 = hash( F([(2,1)]) )
            sage: H1 == H2
            False
        """
        return hash(frozenset(self._monomial.items()))

    def _mul_(self, other):
        """
        Multiply ``self`` by ``other``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: a*b^2*e*d
            F[0]*F[1]^2*F[3]*F[4]
        """
        return self.__class__(self.parent(),
                              blas.add(self._monomial, other._monomial))

    def __pow__(self, n):
        """
        Raise ``self`` to the power of ``n``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: x = a*b^2*e*d; x
            F[0]*F[1]^2*F[3]*F[4]
            sage: x^3
            F[0]^3*F[1]^6*F[3]^3*F[4]^3
            sage: x^0
            1
        """
        if not isinstance(n, (int, Integer)):
            raise TypeError(f"Argument n (= {n}) must be an integer")
        if n < 0:
            raise ValueError(f"Argument n (= {n}) must be positive")
        if n == 1:
            return self
        if n == 0:
            return self.parent().one()
        return self.__class__(self.parent(), {k:v*n for k,v in self._monomial.items()})

    def __floordiv__(self, elt):
        """
        Cancel the element ``elt`` out of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: elt = a*b*c^3*d^2; elt
            F[0]*F[1]*F[2]^3*F[3]^2
            sage: elt // a
            F[1]*F[2]^3*F[3]^2
            sage: elt // c
            F[0]*F[1]*F[2]^2*F[3]^2
            sage: elt // (a*b*d^2)
            F[2]^3
            sage: elt // a^4
            Traceback (most recent call last):
            ...
            ValueError: invalid cancellation
            sage: elt // e^4
            Traceback (most recent call last):
            ...
            ValueError: invalid cancellation
        """
        d = copy(self._monomial)
        for k, v in elt._monomial.items():
            if k not in d:
                raise ValueError("invalid cancellation")
            diff = d[k] - v
            if diff < 0:
                raise ValueError("invalid cancellation")
            elif diff == 0:
                del d[k]
            else:
                d[k] = diff
        return self.__class__(self.parent(), d)

    def divides(self, m) -> bool:
        r"""
        Return whether ``self`` divides ``m``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: elt = a*b*c^3*d^2
            sage: a.divides(elt)
            True
            sage: c.divides(elt)
            True
            sage: (a*b*d^2).divides(elt)
            True
            sage: (a^4).divides(elt)
            False
            sage: e.divides(elt)
            False
        """
        other = m._monomial
        return all(k in other and v <= other[k] for k, v in self._monomial.items())

    def __len__(self):
        """
        Return the length of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: elt = a*c^3*b^2*a
            sage: elt.length()
            7
            sage: len(elt)
            7
        """
        return sum(self._monomial.values())

    length = __len__

    def dict(self):
        """
        Return ``self`` as a dictionary.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: a,b,c,d,e = [F.gen(i) for i in range(5)]
            sage: (a*c^3).dict()
            {0: 1, 2: 3}
        """
        return copy(self._monomial)

class IndexedMonoid(Parent, IndexedGenerators, UniqueRepresentation):
    """
    Base class for monoids with an indexed set of generators.

    INPUT:

    - ``indices`` -- the indices for the generators

    For the optional arguments that control the printing, see
    :class:`~sage.structure.indexed_generators.IndexedGenerators`.
    """
    @staticmethod
    def __classcall__(cls, indices, prefix=None, names=None, **kwds):
        """
        TESTS::

            sage: F = FreeAbelianMonoid(index_set=['a','b','c'])
            sage: G = FreeAbelianMonoid(index_set=('a','b','c'))
            sage: H = FreeAbelianMonoid(index_set=tuple('abc'))
            sage: F is G and F is H
            True

            sage: F = FreeAbelianMonoid(index_set=['a','b','c'], latex_bracket=['LEFT', 'RIGHT'])
            sage: F.print_options()['latex_bracket']
            ('LEFT', 'RIGHT')
            sage: F is G
            False
            sage: Groups.Commutative.free()
            Traceback (most recent call last):
            ...
            ValueError: either the indices or names must be given
        """
        names, indices, prefix = parse_indices_names(names, indices, prefix, kwds)
        if prefix is None:
            prefix = "F"

        # bracket or latex_bracket might be lists, so convert
        # them to tuples so that they're hashable.
        bracket = kwds.get('bracket', None)
        if isinstance(bracket, list):
            kwds['bracket'] = tuple(bracket)
        latex_bracket = kwds.get('latex_bracket', None)
        if isinstance(latex_bracket, list):
            kwds['latex_bracket'] = tuple(latex_bracket)

        return super().__classcall__(cls, indices, prefix,
                                     names=names, **kwds)

    def __init__(self, indices, prefix, category=None, names=None, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: TestSuite(F).run()
            sage: F = FreeMonoid(index_set=tuple('abcde'))
            sage: TestSuite(F).run()

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: TestSuite(F).run()
            sage: F = FreeAbelianMonoid(index_set=tuple('abcde'))
            sage: TestSuite(F).run()
        """
        self._indices = indices
        category = Monoids().or_subcategory(category)
        if indices.cardinality() == 0:
            category = category.Finite()
        else:
            category = category.Infinite()
        if indices in Sets().Finite():
            category = category.FinitelyGeneratedAsMagma()
        Parent.__init__(self, names=names, category=category)

        # ignore the optional 'key' since it only affects CachedRepresentation
        kwds.pop('key', None)
        IndexedGenerators.__init__(self, indices, prefix, **kwds)

    def _first_ngens(self, n):
        """
        Used by the preparser for ``F.<x> = ...``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: F._first_ngens(3)
            (F[0], F[1], F[-1])
        """
        it = iter(self._indices)
        return tuple(self.gen(next(it)) for i in range(n))

    def _element_constructor_(self, x=None):
        """
        Create an element of this abelian monoid from ``x``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: F(F.gen(2))
            F[2]
            sage: F([[1, 3], [-2, 12]])
            F[-2]^12*F[1]^3
            sage: F(-5)
            F[-5]
        """
        if x is None:
            return self.one()
        if x in self._indices:
            return self.gen(x)
        return self.element_class(self, x)

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: G = FreeAbelianMonoid(index_set=ZZ)
            sage: G.an_element()
            F[-1]^3*F[0]*F[1]^3
            sage: G = FreeMonoid(index_set=tuple('ab'))
            sage: G.an_element()
            F['a']^2*F['b']^2
        """
        x = self.one()
        I = self._indices
        try:
            x *= self.gen(I.an_element())
        except Exception:
            pass
        try:
            g = iter(self._indices)
            for c in range(1,4):
                x *= self.gen(next(g)) ** c
        except Exception:
            pass
        return x

    def cardinality(self):
        r"""
        Return the cardinality of ``self``, which is `\infty` unless this is
        the trivial monoid.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: F.cardinality()
            +Infinity
            sage: F = FreeMonoid(index_set=())
            sage: F.cardinality()
            1

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: F.cardinality()
            +Infinity
            sage: F = FreeAbelianMonoid(index_set=())
            sage: F.cardinality()
            1
        """
        if self._indices.cardinality() == 0:
            return ZZ.one()
        return infinity

    @cached_method
    def monoid_generators(self):
        """
        Return the monoid generators of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: F.monoid_generators()
            Lazy family (Generator map from Integer Ring to
             Free abelian monoid indexed by Integer Ring(i))_{i in Integer Ring}
            sage: F = FreeAbelianMonoid(index_set=tuple('abcde'))
            sage: sorted(F.monoid_generators())
            [F['a'], F['b'], F['c'], F['d'], F['e']]
        """
        if self._indices.cardinality() == infinity:
            gen = PoorManMap(self.gen, domain=self._indices, codomain=self, name="Generator map")
            return Family(self._indices, gen)
        return Family(self._indices, self.gen)

    gens = monoid_generators

class IndexedFreeMonoid(IndexedMonoid):
    """
    Free monoid with an indexed set of generators.

    INPUT:

    - ``indices`` -- the indices for the generators

    For the optional arguments that control the printing, see
    :class:`~sage.structure.indexed_generators.IndexedGenerators`.

    EXAMPLES::

        sage: F = FreeMonoid(index_set=ZZ)
        sage: F.gen(15)^3 * F.gen(2) * F.gen(15)
        F[15]^3*F[2]*F[15]
        sage: F.gen(1)
        F[1]

    Now we examine some of the printing options::

        sage: F = FreeMonoid(index_set=ZZ, prefix='X', bracket=['|','>'])
        sage: F.gen(2) * F.gen(12)
        X|2>*X|12>
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: FreeMonoid(index_set=ZZ)
            Free monoid indexed by Integer Ring
        """
        return f"Free monoid indexed by {self._indices}"

    Element = IndexedFreeMonoidElement

    @cached_method
    def one(self):
        """
        Return the identity element of ``self``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: F.one()
            1
        """
        return self.element_class(self, ())

    def gen(self, x):
        """
        The generator indexed by ``x`` of ``self``.

        EXAMPLES::

            sage: F = FreeMonoid(index_set=ZZ)
            sage: F.gen(0)
            F[0]
            sage: F.gen(2)
            F[2]

        TESTS::

            sage: F = FreeMonoid(index_set=[1,2])
            sage: F.gen(2)
            F[2]
            sage: F.gen(0)
            Traceback (most recent call last):
            ...
            IndexError: 0 is not in the index set
        """
        if x not in self._indices:
            raise IndexError(f"{x} is not in the index set")
        try:
            return self.element_class(self, ((self._indices(x), ZZ.one()),))
        except (ValueError, TypeError, NotImplementedError): # Backup (e.g., if it is a string)
            return self.element_class(self, ((x, ZZ.one()),))

class IndexedFreeAbelianMonoid(IndexedMonoid):
    """
    Free abelian monoid with an indexed set of generators.

    INPUT:

    - ``indices`` -- the indices for the generators

    For the optional arguments that control the printing, see
    :class:`~sage.structure.indexed_generators.IndexedGenerators`.

    EXAMPLES::

        sage: F = FreeAbelianMonoid(index_set=ZZ)
        sage: F.gen(15)^3 * F.gen(2) * F.gen(15)
        F[2]*F[15]^4
        sage: F.gen(1)
        F[1]

    Now we examine some of the printing options::

        sage: F = FreeAbelianMonoid(index_set=Partitions(), prefix='A', bracket=False, scalar_mult='%')
        sage: F.gen([3,1,1]) * F.gen([2,2])
        A[2, 2]%A[3, 1, 1]

    .. TODO::

        Implement a subclass when the index sets is finite that utilizes
        vectors or the polydict monomials with the index order fixed.
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: FreeAbelianMonoid(index_set=ZZ)
            Free abelian monoid indexed by Integer Ring
        """
        return f"Free abelian monoid indexed by {self._indices}"

    def _element_constructor_(self, x=None):
        """
        Create an element of ``self`` from ``x``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: F(F.gen(2))
            F[2]
            sage: F([[1, 3], [-2, 12]])
            F[-2]^12*F[1]^3
            sage: F({1:3, -2: 12})
            F[-2]^12*F[1]^3

        TESTS::

            sage: F([(1, 3), (1, 2)])
            F[1]^5

            sage: F([(42, 0)])
            1
            sage: F({42: 0})
            1
        """
        if isinstance(x, (list, tuple)):
            d = {}
            for k, v in x:
                if k in d:
                    d[k] += v
                else:
                    d[k] = v
            x = d
        if isinstance(x, dict):
            x = {k: v for k, v in x.items() if v != 0}
        return IndexedMonoid._element_constructor_(self, x)

    Element = IndexedFreeAbelianMonoidElement

    @cached_method
    def one(self):
        """
        Return the identity element of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: F.one()
            1
        """
        return self.element_class(self, {})

    def gen(self, x):
        """
        The generator indexed by ``x`` of ``self``.

        EXAMPLES::

            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: F.gen(0)
            F[0]
            sage: F.gen(2)
            F[2]

        TESTS::

            sage: F = FreeAbelianMonoid(index_set=[1,2])
            sage: F.gen(2)
            F[2]
            sage: F.gen(0)
            Traceback (most recent call last):
            ...
            IndexError: 0 is not in the index set

            sage: F = lie_algebras.VirasoroAlgebra(QQ).pbw_basis().indices(); F
            Free abelian monoid indexed by Disjoint union of Family ({'c'}, Integer Ring)
            sage: F.gen('c')
            PBW['c']
        """
        if x not in self._indices:
            raise IndexError(f"{x} is not in the index set")
        try:
            return self.element_class(self, {self._indices(x): ZZ.one()})
        except (ValueError, TypeError, NotImplementedError):  # Backup (e.g., if it is a string)
            return self.element_class(self, {x: ZZ.one()})
