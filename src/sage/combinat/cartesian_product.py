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

from inspect import isgenerator
import sage.misc.prandom as rnd
from sage.misc.mrange import xmrange_iter, _is_finite, _len
from combinat import CombinatorialClass
from ranker import unrank
from sage.rings.infinity import infinity

def CartesianProduct(*iters):
    """
    Returns the combinatorial class of the Cartesian product of
    \*iters.

    EXAMPLES::

        sage: cp = CartesianProduct([1,2], [3,4]); cp
        Cartesian product of [1, 2], [3, 4]
        sage: cp.list()
        [[1, 3], [1, 4], [2, 3], [2, 4]]

    Note that you must not use a generator-type object that is
    returned by a function (using "yield"). They cannot be copied or
    rewound (you cannot jump back to the beginning), but this is
    necessary to construct the cartesian product::

        sage: def a(n): yield 1*n; yield 2*n
        sage: def b(): yield 'a'; yield 'b'
        sage: CartesianProduct(a(3), b()).list()
        Traceback (most recent call last):
        ...
        ValueError: generators are not allowed, see the
        documentation (type "CartesianProduct?") for a workaround

    You either create a list of all values or you use
    :class:`sage.combinat.misc.IterableFunctionCall` to make a
    (copy-able) iterator::

        sage: from sage.combinat.misc import IterableFunctionCall
        sage: CartesianProduct(IterableFunctionCall(a, 3), IterableFunctionCall(b)).list()
        [[3, 'a'], [3, 'b'], [6, 'a'], [6, 'b']]

    See the documentation for
    :class:`~sage.combinat.misc.IterableFunctionCall` for more
    information.
    """
    if any(isgenerator(i) for i in iters):
        raise ValueError('generators are not allowed, see the documentation '+
                         '(type "CartesianProduct?") for a workaround')
    return CartesianProduct_iters(*iters)

class CartesianProduct_iters(CombinatorialClass):
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

            sage: cp = CartesianProduct([1,2],[3,4])
            sage: [1,3] in cp
            True
            sage: [1,2] in cp
            False
            sage: [1, 3, 1] in cp
            False
        """
        try:
            return len(x) == len(self.iters) and all(x[i] in self.iters[i] for i in range(len(self.iters)))
        except (TypeError, IndexError):
            return False

    def __repr__(self):
        """
        EXAMPLES::

            sage: CartesianProduct(range(2), range(3))
            Cartesian product of [0, 1], [0, 1, 2]
        """
        return "Cartesian product of " + ", ".join(map(str, self.iters))

    def cardinality(self):
        """
        Returns the number of elements in the cartesian product of
        everything in \*iters.

        EXAMPLES::

            sage: CartesianProduct(range(2), range(3)).cardinality()
            6
            sage: CartesianProduct(range(2), xrange(3)).cardinality()
            6
            sage: CartesianProduct(range(2), xrange(3), xrange(4)).cardinality()
            24

        This works correctly for infinite objects::

            sage: CartesianProduct(ZZ, QQ).cardinality()
            +Infinity
            sage: CartesianProduct(ZZ, []).cardinality()
            0
        """
        return self._mrange.cardinality()

    def __len__(self):
        """
        Return the number of elements of the cartesian product.

        OUTPUT:

        An ``int``, the number of elements in the cartesian product. If the
        number of elements is infinite or does not fit into a python ``int``, a
        ``TypeError`` is raised.

        .. SEEALSO::

            :meth:`cardinality`

        EXAMPLES::

            sage: C = CartesianProduct(xrange(3), xrange(4))
            sage: len(C)
            12
            sage: C = CartesianProduct(ZZ, QQ)
            sage: len(C)
            Traceback (most recent call last):
            ...
            TypeError: cardinality does not fit into a Python int.
            sage: C = CartesianProduct(ZZ, [])
            sage: len(C)
            0
        """
        return self._mrange.__len__()

    def list(self):
        """
        Returns

        EXAMPLES::

            sage: CartesianProduct(range(3), range(3)).list()
            [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
            sage: CartesianProduct('dog', 'cat').list()
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
        An iterator for the elements in the cartesian product of the
        iterables \*iters.

        From Recipe 19.9 in the Python Cookbook by Alex Martelli and David
        Ascher.

        EXAMPLES::

            sage: [e for e in CartesianProduct(range(3), range(3))]
            [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
            sage: [e for e in CartesianProduct('dog', 'cat')]
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
        return self._mrange.__iter__()

    def is_finite(self):
        """
        The cartesian product is finite if all of its inputs are
        finite, or if any input is empty.

        EXAMPLES::

            sage: CartesianProduct(ZZ, []).is_finite()
            True
            sage: CartesianProduct(4,4).is_finite()
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
        For finite cartesian products, we can reduce unrank to the
        constituent iterators.

        EXAMPLES::

            sage: C = CartesianProduct(xrange(1000), xrange(1000), xrange(1000))
            sage: C[238792368]
            [238, 792, 368]

        Check for :trac:`15919`::

            sage: FF = IntegerModRing(29)
            sage: C = CartesianProduct(FF, FF, FF)
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
        Returns a random element from the cartesian product of \*iters.

        EXAMPLES::

            sage: CartesianProduct('dog', 'cat').random_element()
            ['d', 'a']
        """
        return [rnd.choice(_) for _ in self.iters]
