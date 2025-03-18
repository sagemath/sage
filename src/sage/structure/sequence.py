# sage_setup: distribution = sagemath-objects
r"""
Finite Homogeneous Sequences

A mutable sequence of elements with a common guaranteed category,
which can be set immutable.

Sequence derives from list, so has all the functionality of lists and
can be used wherever lists are used.  When a sequence is created
without explicitly given the common universe of the elements, the
constructor coerces the first and second element to some
*canonical* common parent, if possible, then the second and
third, etc.  If this is possible, it then coerces everything into the
canonical parent at the end.  (Note that canonical coercion is very
restrictive.)  The sequence then has a function ``universe()``
which returns either the common canonical parent (if the coercion
succeeded), or the category of all objects (Objects()).  So if you
have a list `v` and type::

    sage: v = [1, 2/3, 5]
    sage: w = Sequence(v)
    sage: w.universe()
    Rational Field

then since ``w.universe()`` is `\QQ`, you're guaranteed that all
elements of `w` are rationals::

    sage: v[0].parent()
    Integer Ring
    sage: w[0].parent()
    Rational Field

If you do assignment to `w` this property of being rationals is guaranteed
to be preserved::

    sage: w[0] = 2
    sage: w[0].parent()
    Rational Field
    sage: w[0] = 'hi'
    Traceback (most recent call last):
    ...
    TypeError: unable to convert 'hi' to a rational

However, if you do ``w = Sequence(v)`` and the resulting universe
is ``Objects()``, the elements are not guaranteed to have any
special parent.  This is what should happen, e.g., with finite field
elements of different characteristics::

    sage: v = Sequence([GF(3)(1), GF(7)(1)])
    sage: v.universe()
    Category of objects

You can make a list immutable with ``v.freeze()``.  Assignment is
never again allowed on an immutable list.

Creation of a sequence involves making a copy of the input list, and
substantial coercions.  It can be greatly sped up by explicitly
specifying the universe of the sequence::

    sage: v = Sequence(range(10000), universe=ZZ)
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.persist import register_unpickle_override
import sage.structure.sage_object
import sage.structure.coerce


def Sequence(x, universe=None, check=True, immutable=False, cr=False, cr_str=None, use_sage_types=False):
    r"""
    A mutable list of elements with a common guaranteed universe,
    which can be set immutable.

    A universe is either an object that supports coercion (e.g., a
    parent), or a category.

    INPUT:

    - ``x`` -- list or tuple instance

    - ``universe`` -- (default: ``None``) the universe of elements; if ``None``
      determined using canonical coercions and the entire list of elements.
      If list is empty, is category ``Objects()`` of all objects.

    - ``check`` -- boolean (default: ``True``); whether to coerce the elements of x
      into the universe

    - ``immutable`` -- boolean (default: ``True``); whether or not this sequence is
      immutable

    - ``cr`` -- boolean (default: ``False``); if ``True``, then print a carriage return
      after each comma when calling ``repr()`` on this sequence (see note below)

    - ``cr_str`` -- boolean (default: same as ``cr``); if ``True``, then print a carriage return
      after each comma when calling ``str()`` on this sequence (see note below)

    - ``use_sage_types`` -- boolean (default: ``False``); if ``True``, coerce the
      built-in Python numerical types int, float, complex to the corresponding
      Sage types (this makes functions like ``vector()`` more flexible)

    OUTPUT: a sequence

    EXAMPLES::

        sage: v = Sequence(range(10))
        sage: v.universe()
        <class 'int'>
        sage: v
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    We can request that the built-in Python numerical types be coerced
    to Sage objects::

        sage: v = Sequence(range(10), use_sage_types=True)
        sage: v.universe()
        Integer Ring
        sage: v
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    You can also use seq for "Sequence", which is identical to using
    Sequence::

        sage: v = seq([1,2,1/1]); v
        [1, 2, 1]
        sage: v.universe()
        Rational Field

    Note that assignment coerces if possible,::

        sage: v = Sequence(range(10), ZZ)
        sage: a = QQ(5)
        sage: v[3] = a
        sage: parent(v[3])
        Integer Ring
        sage: parent(a)
        Rational Field
        sage: v[3] = 2/3
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    Sequences can be used absolutely anywhere lists or tuples can be used::

        sage: isinstance(v, list)
        True

    Sequence can be immutable, so entries can't be changed::

        sage: v = Sequence([1,2,3], immutable=True)
        sage: v.is_immutable()
        True
        sage: v[0] = 5
        Traceback (most recent call last):
        ...
        ValueError: object is immutable; please change a copy instead.

    Only immutable sequences are hashable (unlike Python lists),
    though the hashing is potentially slow, since it first involves
    conversion of the sequence to a tuple, and returning the hash of
    that.::

        sage: v = Sequence(range(10), ZZ, immutable=True)
        sage: hash(v) == hash(tuple(range(10)))
        True


    If you really know what you are doing, you can circumvent the type
    checking (for an efficiency gain)::

        sage: list.__setitem__(v, int(1), 2/3)        # bad circumvention
        sage: v
        [0, 2/3, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: list.__setitem__(v, int(1), int(2))     # not so bad circumvention

    You can make a sequence with a new universe from an old sequence.::

        sage: w = Sequence(v, QQ)
        sage: w
        [0, 2, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: w.universe()
        Rational Field
        sage: w[1] = 2/3
        sage: w
        [0, 2/3, 2, 3, 4, 5, 6, 7, 8, 9]

    The default universe for any sequence, if no compatible parent structure
    can be found, is the universe of all Sage objects.

    This example illustrates how every element of a list is taken into account
    when constructing a sequence.::

        sage: v = Sequence([1, 7, 6, GF(5)(3)]); v
        [1, 2, 1, 3]
        sage: v.universe()
        Finite Field of size 5

    .. NOTE::

        ``cr`` and ``cr_str`` is not recommended (because IPython's pretty printer is used);
        nevertheless it is kept for backwards compatibility.

        By default ``Sequence`` are printed using IPython's pretty printer,
        so ``cr`` and ``cr_str`` are not taken into account at all::

            sage: Sequence([1, 2, 3], cr=False)
            [1, 2, 3]
            sage: Sequence([1, 2, 3], cr=True)
            [1, 2, 3]

        Nevertheless, before the pretty printer exists, ``repr()`` is used.
        Now ``cr`` and ``cr_str`` still affects the behavior of ``repr()`` and ``str()``::

            sage: repr(Sequence([1, 2, 3], cr=False))
            '[1, 2, 3]'
            sage: repr(Sequence([1, 2, 3], cr=True))
            '[\n1,\n2,\n3\n]'

        In any case, this behavior should probably not be relied upon.

    TESTS::

        sage: Sequence(["a"], universe=ZZ)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert a to an element of Integer Ring

    Here are some tests for ``cr`` and ``cr_str``, even though they shouldn't be used.
    ``cr_str`` can be weird in this case, but we keep the current implementation
    (the feature is not recommended anyway so it doesn't make much difference)::

        sage: str(Sequence([1, 2, 3], cr=True, cr_str=True))
        '[\n1,\n2,\n3\n]'
        sage: str(Sequence([1, 2, 3], cr=True, cr_str=False))
        '[\n1,\n2,\n3\n]'

    In the opposite case, ``cr_str`` works fine::

        sage: str(Sequence([1, 2, 3], cr=False, cr_str=False))
        '[1, 2, 3]'
        sage: str(Sequence([1, 2, 3], cr=False, cr_str=True))
        '[\n1,\n2,\n3\n]'
        sage: repr(Sequence([1, 2, 3], cr=False, cr_str=False))
        '[1, 2, 3]'
        sage: repr(Sequence([1, 2, 3], cr=False, cr_str=True))
        '[1, 2, 3]'
    """
    if universe is None:
        if isinstance(x, Sequence_generic):
            universe = x.universe()
            x = list(x)
        else:
            try:
                from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
            except ImportError:
                pass
            else:
                if isinstance(x, MPolynomialIdeal):
                    universe = x.ring()
                    x = x.gens()

    if universe is None:
        orig_x = x
        x = list(x)  # make a copy even if x is a list, we're going to change it

        if len(x) == 0:
            from sage.categories.objects import Objects
            universe = Objects()
        else:
            import sage.structure.element
            if use_sage_types:
                # convert any Python built-in numerical types to Sage objects
                x = [sage.structure.coerce.py_scalar_to_element(e) for e in x]
            # start the pairwise coercion
            for i in range(len(x) - 1):
                try:
                    x[i], x[i+1] = sage.structure.element.canonical_coercion(x[i], x[i+1])
                except TypeError:
                    from sage.categories.objects import Objects
                    universe = Objects()
                    x = list(orig_x)
                    check = False  # no point
                    break
            if universe is None:   # no type errors raised.
                universe = sage.structure.element.parent(x[len(x)-1])

    try:
        from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
        from sage.rings.polynomial.pbori.pbori import BooleanMonomialMonoid
        from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
        from sage.rings.quotient_ring import QuotientRing_nc
    except ImportError:
        pass
    else:
        if isinstance(universe, (MPolynomialRing_base, BooleanMonomialMonoid)) or (isinstance(universe, QuotientRing_nc) and isinstance(universe.cover_ring(), MPolynomialRing_base)):
            return PolynomialSequence(x, universe, immutable=immutable, cr=cr, cr_str=cr_str)

    return Sequence_generic(x, universe, check, immutable, cr, cr_str, use_sage_types)


class Sequence_generic(sage.structure.sage_object.SageObject, list):
    """
    A mutable list of elements with a common guaranteed universe,
    which can be set immutable.

    A universe is either an object that supports coercion (e.g., a parent),
    or a category.

    INPUT:

    - ``x`` -- list or tuple instance

    - ``universe`` -- (default: ``None``) the universe of elements; if ``None``
      determined using canonical coercions and the entire list of elements.
      If list is empty, is category ``Objects()`` of all objects.

    - ``check`` -- boolean (default: ``True``); whether to coerce the elements of x
      into the universe

    - ``immutable`` -- boolean (default: ``True``); whether or not this sequence is
      immutable

    - ``cr``, ``cr_str`` -- see :func:`Sequence`

    - ``use_sage_types`` -- boolean (default: ``False``); if ``True``, coerce the
      built-in Python numerical types int, float, complex to the corresponding
      Sage types (this makes functions like ``vector()`` more flexible)

    OUTPUT: a sequence

    EXAMPLES::

        sage: v = Sequence(range(10))
        sage: v.universe()
        <class 'int'>
        sage: v
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    We can request that the built-in Python numerical types be coerced
    to Sage objects::

        sage: v = Sequence(range(10), use_sage_types=True)
        sage: v.universe()
        Integer Ring
        sage: v
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    You can also use seq for "Sequence", which is identical to using Sequence::

        sage: v = seq([1,2,1/1]); v
        [1, 2, 1]
        sage: v.universe()
        Rational Field

    Note that assignment coerces if possible,

    ::

        sage: v = Sequence(range(10), ZZ)
        sage: a = QQ(5)
        sage: v[3] = a
        sage: parent(v[3])
        Integer Ring
        sage: parent(a)
        Rational Field
        sage: v[3] = 2/3
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    Sequences can be used absolutely anywhere lists or tuples can be used::

        sage: isinstance(v, list)
        True

    Sequence can be immutable, so entries can't be changed::

        sage: v = Sequence([1,2,3], immutable=True)
        sage: v.is_immutable()
        True
        sage: v[0] = 5
        Traceback (most recent call last):
        ...
        ValueError: object is immutable; please change a copy instead.

    Only immutable sequences are hashable (unlike Python lists),
    though the hashing is potentially slow, since it first involves
    conversion of the sequence to a tuple, and returning the hash of
    that.

    ::

        sage: v = Sequence(range(10), ZZ, immutable=True)
        sage: hash(v) == hash(tuple(range(10)))
        True


    If you really know what you are doing, you can circumvent the type
    checking (for an efficiency gain)::

        sage: list.__setitem__(v, int(1), 2/3)        # bad circumvention
        sage: v
        [0, 2/3, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: list.__setitem__(v, int(1), int(2))     # not so bad circumvention

    You can make a sequence with a new universe from an old sequence.

    ::

        sage: w = Sequence(v, QQ)
        sage: w
        [0, 2, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: w.universe()
        Rational Field
        sage: w[1] = 2/3
        sage: w
        [0, 2/3, 2, 3, 4, 5, 6, 7, 8, 9]

    The default universe for any sequence, if no compatible parent structure
    can be found, is the universe of all Sage objects.

    This example illustrates how every element of a list is taken into account
    when constructing a sequence.

    ::

        sage: v = Sequence([1, 7, 6, GF(5)(3)]); v
        [1, 2, 1, 3]
        sage: v.universe()
        Finite Field of size 5
    """
    def __init__(self, x, universe=None, check=True, immutable=False,
                 cr=False, cr_str=None, use_sage_types=False):
        r"""
        Create a sequence.

        EXAMPLES::

            sage: Sequence([1..5])
            [1, 2, 3, 4, 5]
            sage: a = Sequence([1..3], universe=QQ, check=False, immutable=True, cr=True, cr_str=False, use_sage_types=True)
            sage: a
            [1, 2, 3]
            sage: a = Sequence([1..5], universe=QQ, check=False, immutable=True, cr_str=True, use_sage_types=True)
            sage: a
            [1, 2, 3, 4, 5]
            sage: a._Sequence_generic__cr_str
            True
            sage: a.__str__()
            '[\n1,\n2,\n3,\n4,\n5\n]'
        """
        self.__hash = None

        self.__cr = cr
        if cr_str is None:
            self.__cr_str = cr
        else:
            self.__cr_str = cr_str

        if isinstance(x, Sequence_generic):
            if universe is None or universe == x.__universe:
                list.__init__(self, x)
                self.__universe = x.__universe
                self._is_immutable = immutable
                return

        self.__universe = universe
        if check:
            x = list(x)
            for i in range(len(x)):
                try:
                    x[i] = universe(x[i])
                except TypeError:
                    raise TypeError("unable to convert {} to an element of {}"
                                    .format(x[i], universe))
        list.__init__(self, x)
        self._is_immutable = immutable

    def reverse(self):
        """
        Reverse the elements of self, in place.

        EXAMPLES::

            sage: B = Sequence([1,2,3])
            sage: B.reverse(); B
            [3, 2, 1]
        """
        self._require_mutable()
        list.reverse(self)

    def __setitem__(self, n, value):
        """
        EXAMPLES::

            sage: a = Sequence([1..5])
            sage: a[2] = 19
            sage: a
            [1, 2, 19, 4, 5]
            sage: a[2] = 'hello'
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'hello' to an integer
            sage: a[2] = '5'
            sage: a
            [1, 2, 5, 4, 5]
            sage: v = Sequence([1,2,3,4], immutable=True)
            sage: v[1:3] = [5,7]
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: v = Sequence([1,2,3,4])
            sage: v[1:3] = [5, 3/1]
            sage: v
            [1, 5, 3, 4]
            sage: type(v[2])
            <class 'sage.rings.integer.Integer'>
        """
        self._require_mutable()
        if isinstance(n, slice):
            y = [self.__universe(x) for x in value]
        else:
            y = self.__universe(value)
        list.__setitem__(self, n, y)
        self.__hash = None

    def __getitem__(self, n):
        """
        EXAMPLES::

            sage: v = Sequence([1,2,3,4], immutable=True)
            sage: w = v[2:]
            sage: w
            [3, 4]
            sage: type(w)
            <class 'sage.structure.sequence.Sequence_generic'>
            sage: w[0] = 5; w
            [5, 4]
            sage: v
            [1, 2, 3, 4]
        """
        if isinstance(n, slice):
            return Sequence(list.__getitem__(self, n),
                            universe=self.__universe,
                            check=False,
                            immutable=False,
                            cr=self.__cr)

        return list.__getitem__(self, n)

    # We have to define the *slice functions as long as Sage uses Python 2.*
    # otherwise the inherited *slice functions from list are called
    def __getslice__(self, i, j):
        return self.__getitem__(slice(i, j))

    def __setslice__(self, i, j, value):
        return self.__setitem__(slice(i, j), value)

    def append(self, x):
        """
        EXAMPLES::

            sage: v = Sequence([1,2,3,4], immutable=True)
            sage: v.append(34)
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: v = Sequence([1/3,2,3,4])
            sage: v.append(4)
            sage: type(v[4])
            <class 'sage.rings.rational.Rational'>
        """
        self._require_mutable()
        y = self.__universe(x)
        list.append(self, y)

    def extend(self, iterable):
        """
        Extend list by appending elements from the iterable.

        EXAMPLES::

            sage: B = Sequence([1,2,3])
            sage: B.extend(range(4))
            sage: B
            [1, 2, 3, 0, 1, 2, 3]
        """
        self._require_mutable()
        v = [self.__universe(x) for x in iterable]
        list.extend(self, v)

    def insert(self, index, object):
        """
        Insert object before index.

        EXAMPLES::

            sage: B = Sequence([1,2,3])
            sage: B.insert(10, 5)
            sage: B
            [1, 2, 3, 5]
        """
        self._require_mutable()
        list.insert(self, index, self.__universe(object))

    def pop(self, index=-1):
        """
        Remove and return item at index ``index`` (default: last).

        EXAMPLES::

            sage: B = Sequence([1,2,3])
            sage: B.pop(1)
            2
            sage: B
            [1, 3]
        """
        self._require_mutable()
        return list.pop(self, index)

    def remove(self, value):
        """
        Remove first occurrence of ``value``.

        EXAMPLES::

            sage: B = Sequence([1,2,3])
            sage: B.remove(2)
            sage: B
            [1, 3]
        """
        self._require_mutable()
        list.remove(self, value)

    def sort(self, key=None, reverse=False):
        """
        Sort this list *IN PLACE*.

        INPUT:

        - ``key`` -- see Python ``list sort``

        - ``reverse`` -- see Python ``list sort``

        EXAMPLES::

            sage: B = Sequence([3,2,1/5])
            sage: B.sort()
            sage: B
            [1/5, 2, 3]
            sage: B.sort(reverse=True); B
            [3, 2, 1/5]
        """
        self._require_mutable()
        list.sort(self, key=key, reverse=reverse)

    def __hash__(self):
        """
        EXAMPLES::

            sage: a = Sequence([1..5])
            sage: a.__hash__()
            Traceback (most recent call last):
            ...
            ValueError: mutable sequences are unhashable
            sage: a[0] = 10
            sage: a.set_immutable()
            sage: a.__hash__() == hash(a) == hash(tuple(a))
            True
        """
        if not self._is_immutable:
            raise ValueError("mutable sequences are unhashable")
        if self.__hash is None:
            self.__hash = hash(tuple(self))
        return self.__hash

    def _repr_(self):
        """
        Return a string representation of this sequence.
        Typically, :meth:`_repr_pretty_` is used instead of this method.

        EXAMPLES::

            sage: Sequence([1,2/3,-2/5])._repr_()
            '[1, 2/3, -2/5]'
            sage: print(Sequence([1,2/3,-2/5], cr=True)._repr_())
            [
            1,
            2/3,
            -2/5
            ]
        """
        if self.__cr:
            return '[\n' + ',\n'.join(repr(x) for x in self) + '\n]'
        else:
            return list.__repr__(self)

    def _repr_pretty_(self, p, cycle):
        """
        For pretty printing in the Sage command prompt.

        Since ``Sequence`` inherits from ``list``, we just use IPython's built-in
        ``list`` pretty printer.
        When :issue:`36801` is fixed, this function will be redundant.

        EXAMPLES::

            sage: Sequence([1,2/3,-2/5])  # indirect doctest
            [1, 2/3, -2/5]
        """
        p.pretty(list(self))

    def _latex_(self):
        r"""
        TESTS::

            sage: t= Sequence([sqrt(x), exp(x), x^(x-1)], universe=SR); t               # needs sage.symbolic
            [sqrt(x), e^x, x^(x - 1)]
            sage: t._latex_()                                                           # needs sage.symbolic
            '\\left[\\sqrt{x}, e^{x}, x^{x - 1}\\right]'
            sage: latex(t)                                                              # needs sage.symbolic
            \left[\sqrt{x}, e^{x}, x^{x - 1}\right]
        """
        from sage.misc.latex import list_function as list_latex_function
        return list_latex_function(self)

    def __str__(self):
        r"""
        EXAMPLES::

            sage: s = Sequence([1,2,3], cr=False)
            sage: str(s)
            '[1, 2, 3]'
            sage: repr(s)
            '[1, 2, 3]'
            sage: print(s)
            [1, 2, 3]
            sage: s = Sequence([1,2,3], cr=True)
            sage: str(s)
            '[\n1,\n2,\n3\n]'
        """
        if self.__cr_str:
            return '[\n' + ',\n'.join(str(x) for x in self) + '\n]'
        else:
            return list.__str__(self)

    def universe(self):
        """
        Return the universe of the sequence.

        EXAMPLES::

            sage: Sequence([1, 2/3, -2/5]).universe()
            Rational Field
            sage: Sequence([1, 2/3, '-2/5']).universe()
            Category of objects
        """
        return self.__universe

    def _require_mutable(self):
        """
        EXAMPLES::

            sage: a = Sequence([1,2/3,'-2/5'])
            sage: a._require_mutable()
            sage: a.set_immutable()
            sage: a._require_mutable()
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        if self._is_immutable:
            raise ValueError("object is immutable; please change a copy instead.")

    def set_immutable(self):
        """
        Make this object immutable, so it can never again be changed.

        EXAMPLES::

            sage: v = Sequence([1, 2, 3, 4/5])
            sage: v[0] = 5
            sage: v
            [5, 2, 3, 4/5]
            sage: v.set_immutable()
            sage: v[3] = 7
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        self._is_immutable = True

    def is_immutable(self):
        """
        Return ``True`` if this object is immutable (can not be changed)
        and ``False`` if it is not.

        To make this object immutable use :meth:`set_immutable`.

        EXAMPLES::

            sage: v = Sequence([1, 2, 3, 4/5])
            sage: v[0] = 5
            sage: v
            [5, 2, 3, 4/5]
            sage: v.is_immutable()
            False
            sage: v.set_immutable()
            sage: v.is_immutable()
            True
        """
        try:
            return self._is_immutable
        except AttributeError:
            return False

    def is_mutable(self):
        """
        EXAMPLES::

            sage: a = Sequence([1, 2/3, -2/5])
            sage: a.is_mutable()
            True
            sage: a[0] = 100
            sage: type(a[0])
            <class 'sage.rings.rational.Rational'>
            sage: a.set_immutable()
            sage: a[0] = 50
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: a.is_mutable()
            False
        """
        try:
            return not self._is_immutable
        except AttributeError:
            return True

    def __reduce__(self):
        """
        Implement pickling for sequences.

        TESTS::

            sage: v = Sequence([1..5])
            sage: w = loads(dumps(v))
            sage: v == w
            True
            sage: w.is_mutable()
            True
            sage: v.set_immutable()
            sage: w = loads(dumps(v))
            sage: w.is_mutable()
            False
        """
        args = (list(self), self.__universe, False,
                self._is_immutable, self.__cr_str)
        return type(self), args

    def __copy__(self):
        """
        Return a copy of this sequence.

        EXAMPLES::

            sage: s = seq(range(10))
            sage: t = copy(s)
            sage: t == s
            True
            sage: t.is_immutable() == s.is_immutable()
            True
            sage: t.is_mutable() == s.is_mutable()
            True
        """
        return Sequence(self, universe=self.__universe,
                        check=False,
                        immutable=self._is_immutable,
                        cr=self.__cr_str)

    def __getattr__(self, name):
        """
        Strictly for unpickling old 'Sequences'.

        INPUT:

        - ``name`` -- some string

        TESTS::

            sage: S = Sequence([])
            sage: del S._Sequence_generic__universe
            sage: S.universe()
            Traceback (most recent call last):
            ...
            AttributeError: 'Sequence_generic' object has no attribute '_Sequence_generic__universe'...
            sage: S._Sequence__universe = 'foobar'
            sage: S.universe()
            'foobar'

        We test that :issue:`13998` is fixed::

            sage: S = Sequence([])
            sage: S.set_immutable()
            sage: del S._Sequence_generic__hash
            sage: hash(S)
            Traceback (most recent call last):
            ...
            AttributeError: 'Sequence_generic' object has no attribute '_Sequence_generic__hash'...
            sage: S._Sequence__hash = int(34)
            sage: hash(S)
            34
        """
        if name == "_Sequence_generic__cr" and hasattr(self, "_Sequence__cr"):
            self.__cr = self._Sequence__cr
            return self.__cr
        if name == "_Sequence_generic__cr_str" and hasattr(self, "_Sequence__cr_str"):
            self.__cr_str = self._Sequence__cr_str
            return self.__cr_str
        if name == "_Sequence_generic__immutable" and hasattr(self, "_Sequence__immutable"):
            self.__immutable = self._Sequence__immutable
            return self.__immutable
        if name == "_Sequence_generic__universe" and hasattr(self, "_Sequence__universe"):
            self.__universe = self._Sequence__universe
            return self.__universe
        if name == "_Sequence_generic__hash" and hasattr(self, "_Sequence__hash"):
            self.__hash = self._Sequence__hash
            return self.__hash

        raise AttributeError("'Sequence_generic' object has no attribute '%s'" % name)


seq = Sequence

register_unpickle_override('sage.structure.sequence', 'Sequence', Sequence_generic)
