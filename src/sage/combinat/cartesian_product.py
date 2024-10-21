# sage_setup: distribution = sagemath-modules
r"""
Cartesian Products
"""
# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.enumerated_sets import EnumeratedSets
from sage.sets.set_from_iterator import EnumeratedSetFromIterator

import sage.misc.prandom as rnd
from sage.misc.mrange import xmrange_iter, _is_finite, _len
from .ranker import unrank
from sage.rings.infinity import infinity


class CartesianProduct_iters(EnumeratedSetFromIterator):
    r"""
    Cartesian product of finite sets.

    This class will soon be deprecated (see :issue:`18411` and :issue:`19195`).
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
        <class 'list'>
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

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: cp = CartesianProduct_iters([1,2],[3,4]); cp
            Cartesian product of [1, 2], [3, 4]
            sage: loads(dumps(cp)) == cp
            True
            sage: TestSuite(cp).run(skip='_test_an_element')

        Check that :issue:`24558` is fixed::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: from sage.sets.set_from_iterator import EnumeratedSetFromIterator
            sage: I = EnumeratedSetFromIterator(Integers)
            sage: CartesianProduct_iters(I, I)
            Cartesian product of {0, 1, -1, 2, -2, ...}, {0, 1, -1, 2, -2, ...}
        """
        self.iters = iters
        self._mrange = xmrange_iter(iters)
        category = EnumeratedSets()
        try:
            category = category.Finite() if self.is_finite() else category.Infinite()
        except ValueError:  # Unable to determine if it is finite or not
            pass

        def iterfunc():
            # we can not use self.__iterate__ directly because
            # that leads to an infinite recursion in __eq__
            return self.__iterate__()
        name = "Cartesian product of " + ", ".join(map(str, self.iters))
        EnumeratedSetFromIterator.__init__(self, iterfunc,
                                           name=name,
                                           category=category,
                                           cache=False)

    def __hash__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: cp = CartesianProduct_iters((1,2), (3,4))
            sage: hash(cp) == CartesianProduct_iters((1,2), (3,4))
        """
        return hash(tuple(self.iters))

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

    def __reduce__(self):
        r"""
        Support for pickle.

        TESTS::

            sage: cp = cartesian_product([[1,2],range(9)])
            sage: loads(dumps(cp)) == cp
            True
        """
        return (self.__class__, (self.iters))

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: CartesianProduct_iters(list(range(2)), list(range(3)))
            Cartesian product of [0, 1], [0, 1, 2]
        """
        return "Cartesian product of " + ", ".join(map(str, self.iters))

    def cardinality(self):
        r"""
        Return the number of elements in the Cartesian product of
        everything in \*iters.

        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: CartesianProduct_iters(range(2), range(3)).cardinality()
            6
            sage: CartesianProduct_iters(range(2), range(3)).cardinality()
            6
            sage: CartesianProduct_iters(range(2), range(3), range(4)).cardinality()
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
        :exc:`TypeError` is raised.

        .. SEEALSO::

            :meth:`cardinality`

        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: C = CartesianProduct_iters(range(3), range(4))
            sage: len(C)
            12
            sage: C = CartesianProduct_iters(ZZ, QQ)
            sage: len(C)
            Traceback (most recent call last):
            ...
            TypeError: cardinality does not fit into a Python int
            sage: C = CartesianProduct_iters(ZZ, [])
            sage: len(C)
            0
        """
        return len(self._mrange)

    def list(self):
        """
        Return.

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

    def __iterate__(self):
        r"""
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
            ValueError: unable to determine whether this product is finite
        """
        finites = [_is_finite(L, fallback=None) for L in self.iters]
        if any(f is None for f in finites):
            raise ValueError("unable to determine whether this product is finite")
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
            sage: C = CartesianProduct_iters(range(1000), range(1000), range(1000))
            sage: C[238792368]
            [238, 792, 368]

        Check for :issue:`15919`::

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
        return [unrank(L, i) for L, i in zip(self.iters, positions)]

    def random_element(self):
        r"""
        Return a random element from the Cartesian product of \*iters.

        EXAMPLES::

            sage: from sage.combinat.cartesian_product import CartesianProduct_iters
            sage: c = CartesianProduct_iters('dog', 'cat').random_element()
            sage: c in CartesianProduct_iters('dog', 'cat')
            True
        """
        return [rnd.choice(w) for w in self.iters]
