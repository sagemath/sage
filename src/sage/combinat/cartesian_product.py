r"""
Cartesian Products
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from inspect import isgenerator
import sage.misc.prandom as rnd
from sage.misc.mrange import xmrange_iter, _is_finite, _len
from .combinat import CombinatorialClass
from .ranker import unrank
from sage.rings.infinity import infinity


def CartesianProduct(*iters):
    """
    This is deprecated. Use :obj:`cartesian_product` instead.

    EXAMPLES::

        sage: cp = CartesianProduct([1,2], [3,4]); cp
        doctest:...: DeprecationWarning: CartesianProduct is deprecated. Use
        cartesian_product instead
        See http://trac.sagemath.org/18411 for details.
        The Cartesian product of ({1, 2}, {3, 4})
        sage: cp.list()
        [(1, 3), (1, 4), (2, 3), (2, 4)]

    Note that you must not use a generator-type object that is
    returned by a function (using "yield"). They cannot be copied or
    rewound (you cannot jump back to the beginning), but this is
    necessary to construct the Cartesian product::

        sage: def a(n): yield 1*n; yield 2*n
        sage: def b(): yield 'a'; yield 'b'
        sage: CartesianProduct(a(3), b()).list()
        Traceback (most recent call last):
        ...
        ValueError: generators are not allowed, see the
        documentation (type "CartesianProduct?") for a workaround

    The usage of iterable is also deprecated, so the following will no longer be
    supported::

        sage: from sage.combinat.misc import IterableFunctionCall
        sage: C = CartesianProduct(IterableFunctionCall(a, 3), IterableFunctionCall(b))
        doctest:...: DeprecationWarning: Usage of IterableFunctionCall in
        CartesianProduct is deprecated. You can use EnumeratedSetFromIterator
        (in sage.sets.set_from_iterator) instead.
        See http://trac.sagemath.org/18411 for details.
        sage: list(C)
        doctest:...: UserWarning: Sage is not able to determine whether the
        factors of this Cartesian product are finite. The lexicographic ordering
        might not go through all elements.
        [(3, 'a'), (3, 'b'), (6, 'a'), (6, 'b')]

    You might use
    :class:`~sage.sets.set_from_iterator.EnumeratedSetFromIterator` for that
    purpose.::

        sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
        sage: A = EnumeratedSetFromIterator(a, (3,), category=FiniteEnumeratedSets())
        sage: B = EnumeratedSetFromIterator(b, category=FiniteEnumeratedSets())
        sage: C = cartesian_product([A, B])
        sage: C.list()
        [(3, 'a'), (3, 'b'), (6, 'a'), (6, 'b')]
    """
    if any(isgenerator(i) for i in iters):
        raise ValueError('generators are not allowed, see the documentation '+
                         '(type "CartesianProduct?") for a workaround')

    from sage.misc.superseded import deprecation
    deprecation(18411, "CartesianProduct is deprecated. Use cartesian_product instead")

    from sage.combinat.misc import IterableFunctionCall
    from sage.sets.set_from_iterator import EnumeratedSetFromIterator
    deprecate_ifc = False
    iiters = []
    for a in iters:
        if isinstance(a, IterableFunctionCall):
            deprecate_ifc = True
            iiters.append(EnumeratedSetFromIterator(a.f, a.args, a.kwargs))
        else:
            iiters.append(a)
    iters = tuple(iiters)

    if deprecate_ifc:
        deprecation(18411, """Usage of IterableFunctionCall in CartesianProduct is deprecated. You can use EnumeratedSetFromIterator (in sage.sets.set_from_iterator) instead.""")

    from sage.categories.cartesian_product import cartesian_product
    return cartesian_product(iters)

class CartesianProduct_iters(CombinatorialClass):
    r"""
    Cartesian product of finite sets.

    This class will soon be deprecated (see :trac:`18411` and :trac:`19195`).
    One should instead use the functorial construction
    :class:`cartesian_product <sage.categories.cartesian_product.CartesianProductFunctor>`.
    The main differences in behavior are:

    - construction: ``CartesianProduct`` takes as many argument as
      there are factors whereas ``cartesian_product`` takes a single
      list (or iterable) of factors;

    - representation of elements: elements are represented by plain
      Python list for ``CartesianProduct`` versus a custom element
      class for ``cartesian_product``;

    - membership testing: because of the above, plain Python lists are
      not considered as elements of a ``cartesian_product``.

    All of these is illustrated in the examples below.

    EXAMPLES::

        sage: F1 = ['a', 'b']
        sage: F2 = [1, 2, 3, 4]
        sage: F3 = Permutations(3)
        sage: from sage.combinat.cartesian_product import CartesianProduct_iters
        sage: C = CartesianProduct_iters(F1, F2, F3)
        sage: c = cartesian_product([F1, F2, F3])

        sage: type(C.an_element())
        <type 'list'>
        sage: type(c.an_element())
        <class 'sage.sets.cartesian_product.CartesianProduct_with_category.element_class'>

        sage: l = ['a', 1, Permutation([3,2,1])]
        sage: l in C
        True
        sage: l in c
        False
        sage: elt = c(l)
        sage: elt
        ('a', 1, [3, 2, 1])
        sage: elt in c
        True
        sage: elt.parent() is c
        True
    """
    def __init__(self, *iters):
        """
        TESTS::

            sage: import sage.combinat.cartesian_product as cartesian_product
            sage: cp = cartesian_product.CartesianProduct_iters([1,2],[3,4]); cp
            Cartesian product of [1, 2], [3, 4]
            sage: loads(dumps(cp)) == cp
            True
        """
        self.iters = iters
        self._mrange = xmrange_iter(iters)
        CombinatorialClass.__init__(self)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: cp = CartesianProduct_iters([1,2],[3,4])
            sage: [1,3] in cp
            True
            sage: [1,2] in cp
            False
            sage: [1, 3, 1] in cp
            False

        Note that it differs with the behavior of Cartesian products::

            sage: cp = cartesian_product([[1,2], [3,4]])
            sage: [1,3] in cp
            False
        """
        try:
            return len(x) == len(self.iters) and all(x[i] in self.iters[i] for i in range(len(self.iters)))
        except (TypeError, IndexError):
            return False

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: CartesianProduct_iters(range(2), range(3))
            Cartesian product of [0, 1], [0, 1, 2]
        """
        return "Cartesian product of " + ", ".join(map(str, self.iters))

    def cardinality(self):
        """
        Returns the number of elements in the Cartesian product of
        everything in \*iters.

        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: CartesianProduct_iters(range(2), range(3)).cardinality()
            6
            sage: CartesianProduct_iters(range(2), xrange(3)).cardinality()
            6
            sage: CartesianProduct_iters(range(2), xrange(3), xrange(4)).cardinality()
            24

        This works correctly for infinite objects::

            sage: CartesianProduct_iters(ZZ, QQ).cardinality()
            +Infinity
            sage: CartesianProduct_iters(ZZ, []).cardinality()
            0
        """
        return self._mrange.cardinality()

    def __len__(self):
        """
        Return the number of elements of the Cartesian product.

        OUTPUT:

        An ``int``, the number of elements in the Cartesian product. If the
        number of elements is infinite or does not fit into a python ``int``, a
        ``TypeError`` is raised.

        .. SEEALSO::

            :meth:`cardinality`

        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: C = CartesianProduct_iters(xrange(3), xrange(4))
            sage: len(C)
            12
            sage: C = CartesianProduct_iters(ZZ, QQ)
            sage: len(C)
            Traceback (most recent call last):
            ...
            TypeError: cardinality does not fit into a Python int.
            sage: C = CartesianProduct_iters(ZZ, [])
            sage: len(C)
            0
        """
        return len(self._mrange)

    def list(self):
        """
        Returns

        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: CartesianProduct_iters(range(3), range(3)).list()
            [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
            sage: CartesianProduct_iters('dog', 'cat').list()
            [['d', 'c'],
             ['d', 'a'],
             ['d', 't'],
             ['o', 'c'],
             ['o', 'a'],
             ['o', 't'],
             ['g', 'c'],
             ['g', 'a'],
             ['g', 't']]
        """
        return [e for e in self]


    def __iter__(self):
        """
        An iterator for the elements in the Cartesian product of the
        iterables \*iters.

        From Recipe 19.9 in the Python Cookbook by Alex Martelli and David
        Ascher.

        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: [e for e in CartesianProduct_iters(range(3), range(3))]
            [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
            sage: [e for e in CartesianProduct_iters('dog', 'cat')]
            [['d', 'c'],
             ['d', 'a'],
             ['d', 't'],
             ['o', 'c'],
             ['o', 'a'],
             ['o', 't'],
             ['g', 'c'],
             ['g', 'a'],
             ['g', 't']]
        """
        return iter(self._mrange)

    def is_finite(self):
        """
        The Cartesian product is finite if all of its inputs are
        finite, or if any input is empty.

        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: CartesianProduct_iters(ZZ, []).is_finite()
            True
            sage: CartesianProduct_iters(4,4).is_finite()
            Traceback (most recent call last):
            ...
            ValueError: Unable to determine whether this product is finite
        """
        finites = [_is_finite(L, fallback=None) for L in self.iters]
        if any(f is None for f in finites):
            raise ValueError("Unable to determine whether this product is finite")
        if all(f is True for f in finites):
            return True
        lens = [_len(L) for L in self.iters]
        if any(l == 0 for l in lens):
            return True
        return False

    def unrank(self, x):
        """
        For finite Cartesian products, we can reduce unrank to the
        constituent iterators.

        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: C = CartesianProduct_iters(xrange(1000), xrange(1000), xrange(1000))
            sage: C[238792368]
            [238, 792, 368]

        Check for :trac:`15919`::

            sage: FF = IntegerModRing(29)
            sage: C = CartesianProduct_iters(FF, FF, FF)
            sage: C.unrank(0)
            [0, 0, 0]
        """
        try:
            lens = [_len(it) for it in self.iters]
        except (TypeError, AttributeError):
            return CartesianProduct_iters.unrank(self, x)
        positions = []
        for n in lens:
            if n is infinity:
                return CartesianProduct_iters.unrank(self, x)
            if n == 0:
                raise IndexError("Cartesian Product is empty")
            positions.append(x % n)
            x = x // n
        if x != 0:
            raise IndexError("x larger than the size of the Cartesian Product")
        positions.reverse()
        return [unrank(L, i) for L,i in zip(self.iters, positions)]

    def random_element(self):
        """
        Returns a random element from the Cartesian product of \*iters.

        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: CartesianProduct_iters('dog', 'cat').random_element()
            ['d', 'a']
        """
        return [rnd.choice(_) for _ in self.iters]
