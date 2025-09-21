r"""
Miscellaneous
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

from sage.misc.misc_c import prod


class DoublyLinkedList:
    """
    A doubly linked list class that provides constant time hiding and
    unhiding of entries.

    Note that this list's indexing is 1-based.

    EXAMPLES::

        sage: dll = sage.combinat.misc.DoublyLinkedList([1,2,3]); dll
        Doubly linked list of [1, 2, 3]: [1, 2, 3]
        sage: dll.hide(1); dll
        Doubly linked list of [1, 2, 3]: [2, 3]
        sage: dll.unhide(1); dll
        Doubly linked list of [1, 2, 3]: [1, 2, 3]
        sage: dll.hide(2); dll
        Doubly linked list of [1, 2, 3]: [1, 3]
        sage: dll.unhide(2); dll
        Doubly linked list of [1, 2, 3]: [1, 2, 3]
    """

    def __init__(self, l):
        """
        TESTS::

            sage: dll = sage.combinat.misc.DoublyLinkedList([1,2,3])
            sage: dll == loads(dumps(dll))
            True
        """
        n = len(l)
        self.l = l
        self.next_value = {}
        self.next_value['begin'] = l[0]
        self.next_value[l[n-1]] = 'end'
        for i in range(n-1):
            self.next_value[l[i]] = l[i+1]

        self.prev_value = {}
        self.prev_value['end'] = l[-1]
        self.prev_value[l[0]] = 'begin'
        for i in range(1,n):
            self.prev_value[l[i]] = l[i-1]

    def __eq__(self, other):
        """
        TESTS::

            sage: dll = sage.combinat.misc.DoublyLinkedList([1,2,3])
            sage: dll2 = sage.combinat.misc.DoublyLinkedList([1,2,3])
            sage: dll == dll2
            True
            sage: dll.hide(1)
            sage: dll == dll2
            False
        """
        return (isinstance(other, DoublyLinkedList) and
            self.l == other.l and
            self.next_value == other.next_value and
            self.prev_value == other.prev_value)

    def __ne__(self, other):
        """
        TESTS::

            sage: dll = sage.combinat.misc.DoublyLinkedList([1,2,3])
            sage: dll2 = sage.combinat.misc.DoublyLinkedList([1,2,3])
            sage: dll != dll2
            False
            sage: dll.hide(1)
            sage: dll != dll2
            True
        """
        return not (self == other)

    def __repr__(self):
        """
        TESTS::

            sage: repr(sage.combinat.misc.DoublyLinkedList([1,2,3]))
            'Doubly linked list of [1, 2, 3]: [1, 2, 3]'
        """
        return "Doubly linked list of %s: %s" % (self.l, list(self))

    def __iter__(self):
        """
        TESTS::

            sage: dll = sage.combinat.misc.DoublyLinkedList([1,2,3])
            sage: list(dll)
            [1, 2, 3]
        """
        j = self.next_value['begin']
        while j != 'end':
            yield j
            j = self.next_value[j]

    def hide(self, i):
        """
        TESTS::

            sage: dll = sage.combinat.misc.DoublyLinkedList([1,2,3])
            sage: dll.hide(1)
            sage: list(dll)
            [2, 3]
        """
        self.next_value[self.prev_value[i]] = self.next_value[i]
        self.prev_value[self.next_value[i]] = self.prev_value[i]

    def unhide(self, i):
        """
        TESTS::

            sage: dll = sage.combinat.misc.DoublyLinkedList([1,2,3])
            sage: dll.hide(1); dll.unhide(1)
            sage: list(dll)
            [1, 2, 3]
        """
        self.next_value[self.prev_value[i]] = i
        self.prev_value[self.next_value[i]] = i

    def head(self):
        """
        TESTS::

            sage: dll = sage.combinat.misc.DoublyLinkedList([1,2,3])
            sage: dll.head()
            1
            sage: dll.hide(1)
            sage: dll.head()
            2
        """
        return self.next_value['begin']

    def next(self, j):
        """
        TESTS::

            sage: dll = sage.combinat.misc.DoublyLinkedList([1,2,3])
            sage: dll.next(1)
            2
            sage: dll.hide(2)
            sage: dll.next(1)
            3
        """
        return self.next_value[j]

    def prev(self, j):
        """
        TESTS::

            sage: dll = sage.combinat.misc.DoublyLinkedList([1,2,3])
            sage: dll.prev(3)
            2
            sage: dll.hide(2)
            sage: dll.prev(3)
            1
        """
        return self.prev_value[j]


def _monomial_exponent_to_lower_factorial(me, x):
    r"""
    Convert a tuple of exponents to the monomial obtained by replacing
    each me[i] with `x_i*(x_i - 1)*\cdots*(x_i - a_i + 1)`

    EXAMPLES::

        sage: from sage.combinat.misc import _monomial_exponent_to_lower_factorial
        sage: R.<x,y,z> = QQ[]
        sage: a = R.gens()
        sage: _monomial_exponent_to_lower_factorial(([1,0,0]),a)
        x
        sage: _monomial_exponent_to_lower_factorial(([2,0,0]),a)
        x^2 - x
        sage: _monomial_exponent_to_lower_factorial(([0,2,0]),a)
        y^2 - y
        sage: _monomial_exponent_to_lower_factorial(([1,1,0]),a)
        x*y
        sage: _monomial_exponent_to_lower_factorial(([1,1,2]),a)
        x*y*z^2 - x*y*z
        sage: _monomial_exponent_to_lower_factorial(([2,2,2]),a)
        x^2*y^2*z^2 - x^2*y^2*z - x^2*y*z^2 - x*y^2*z^2 + x^2*y*z + x*y^2*z + x*y*z^2 - x*y*z
    """
    return prod(x[i] - j for i, mei in enumerate(me) for j in range(mei))


def umbral_operation(poly):
    r"""
    Return the umbral operation `\downarrow` applied to poly.

    The umbral operation replaces each instance of
    `x_i^{a_i}` with
    `x_i*(x_i - 1)*\cdots*(x_i - a_i + 1)`.

    EXAMPLES::

        sage: P = PolynomialRing(QQ, 2, 'x')
        sage: x = P.gens()
        sage: from sage.combinat.misc import umbral_operation
        sage: umbral_operation(x[0]^3) == x[0]*(x[0]-1)*(x[0]-2)
        True
        sage: umbral_operation(x[0]*x[1])
        x0*x1
        sage: umbral_operation(x[0]+x[1])
        x0 + x1
        sage: umbral_operation(x[0]^2*x[1]^2) == x[0]*(x[0]-1)*x[1]*(x[1]-1)
        True
    """
    x = poly.parent().gens()
    exponents = poly.exponents()
    coefficients = poly.coefficients()
    length = len(exponents)
    return sum(coefficients[i]*_monomial_exponent_to_lower_factorial(exponents[i], x) for i in range(length))


class IterableFunctionCall:
    """
    This class wraps functions with a yield statement (generators) by
    an object that can be iterated over. For example,

    EXAMPLES::

        sage: def f(): yield 'a'; yield 'b'

    This does not work::

        sage: for z in f: print(z)
        Traceback (most recent call last):
        ...
        TypeError: 'function' object is not iterable

    Use IterableFunctionCall if you want something like the above to
    work::

        sage: from sage.combinat.misc import IterableFunctionCall
        sage: g = IterableFunctionCall(f)
        sage: for z in g: print(z)
        a
        b

    If your function takes arguments, just put them after the function
    name. You needn't enclose them in a tuple or anything, just put them
    there::

        sage: def f(n, m): yield 'a' * n; yield 'b' * m; yield 'foo'
        sage: g = IterableFunctionCall(f, 2, 3)
        sage: for z in g: print(z)
        aa
        bbb
        foo
    """

    def __init__(self, f, *args, **kwargs):
        """
        EXAMPLES::

            sage: from sage.combinat.misc import IterableFunctionCall
            sage: IterableFunctionCall(iter, [1,2,3])
            Iterable function call <built-in function iter> with args=([1, 2, 3],) and kwargs={}
        """
        self.f = f
        self.args = args
        self.kwargs = kwargs

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.misc import IterableFunctionCall
            sage: list(iter(IterableFunctionCall(iter, [1,2,3])))
            [1, 2, 3]
        """
        return self.f(*self.args, **self.kwargs)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.misc import IterableFunctionCall
            sage: repr(IterableFunctionCall(iter, [1,2,3]))
            'Iterable function call <built-in function iter> with args=([1, 2, 3],) and kwargs={}'
        """
        return "Iterable function call %s with args=%s and kwargs=%s" % (self.f, self.args, self.kwargs)


def check_integer_list_constraints(l, **kwargs):
    """
    EXAMPLES::

        sage: from sage.combinat.misc import check_integer_list_constraints
        sage: cilc = check_integer_list_constraints
        sage: l = [[2,1,3],[1,2],[3,3],[4,1,1]]
        sage: cilc(l, min_part=2)
        [[3, 3]]
        sage: cilc(l, max_part=2)
        [[1, 2]]
        sage: cilc(l, length=2)
        [[1, 2], [3, 3]]
        sage: cilc(l, max_length=2)
        [[1, 2], [3, 3]]
        sage: cilc(l, min_length=3)
        [[2, 1, 3], [4, 1, 1]]
        sage: cilc(l, max_slope=0)
        [[3, 3], [4, 1, 1]]
        sage: cilc(l, min_slope=1)
        [[1, 2]]
        sage: cilc(l, outer=[2,2])
        [[1, 2]]
        sage: cilc(l, inner=[2,2])
        [[3, 3]]

    ::

        sage: cilc([1,2,3], length=3, singleton=True)
        [1, 2, 3]
        sage: cilc([1,2,3], length=2, singleton=True) is None
        True
    """
    if 'singleton' in kwargs and kwargs['singleton']:
        singleton = True
        result = [l]
        n = sum(l)
        del kwargs['singleton']
    else:
        singleton = False
        if l:
            n = sum(l[0])
            result = l
        else:
            return []

    min_part = kwargs.get('min_part', None)
    max_part = kwargs.get('max_part', None)

    min_length = kwargs.get('min_length', None)
    max_length = kwargs.get('max_length', None)

    min_slope = kwargs.get('min_slope', None)
    max_slope = kwargs.get('max_slope', None)

    length = kwargs.get('length', None)

    inner = kwargs.get('inner', None)
    outer = kwargs.get('outer', None)

    # Preprocess the constraints
    if outer is not None:
        max_length = len(outer)
        for i in range(max_length):
            if outer[i] == "inf":
                outer[i] = n+1
    if inner is not None:
        min_length = len(inner)

    if length is not None:
        max_length = length
        min_length = length

    filters = {}
    filters['length'] = lambda x: len(x) == length
    filters['min_part'] = lambda x: min(x) >= min_part
    filters['max_part'] = lambda x: max(x) <= max_part
    filters['min_length'] = lambda x: len(x) >= min_length
    filters['max_length'] = lambda x: len(x) <= max_length
    filters['min_slope'] = lambda x: min((x[i + 1] - x[i] for i in range(len(x) - 1)), default=min_slope + 1) >= min_slope
    filters['max_slope'] = lambda x: max((x[i + 1] - x[i] for i in range(len(x) - 1)), default=max_slope - 1) <= max_slope
    filters['outer'] = lambda x: len(outer) >= len(x) and min(outer[i] - x[i] for i in range(len(x))) >= 0
    filters['inner'] = lambda x: len(x) >= len(inner) and max(inner[i] - x[i] for i in range(len(inner))) <= 0

    for key in kwargs:
        result = [x for x in result if filters[key](x)]

    if singleton:
        try:
            return result[0]
        except IndexError:
            return None
    else:
        return result
