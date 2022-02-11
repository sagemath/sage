"""
Finite Enumerated Sets
"""
# ****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
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
# *****************************************************************************

from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_cat import EmptySetError
from sage.rings.integer import Integer

#################################################################
class FiniteEnumeratedSet(UniqueRepresentation, Parent):
    """
    A class for finite enumerated set.

    Returns the finite enumerated set with elements in ``elements``
    where ``element`` is any (finite) iterable object.

    The main purpose is to provide a variant of ``list`` or ``tuple``,
    which is a parent with an interface consistent with
    ``EnumeratedSets`` and has unique representation.
    The list of the elements is expanded in memory.


    EXAMPLES::

        sage: S = FiniteEnumeratedSet([1, 2, 3])
        sage: S
        {1, 2, 3}
        sage: S.list()
        [1, 2, 3]
        sage: S.cardinality()
        3
        sage: S.random_element()  # random
        1
        sage: S.first()
        1
        sage: S.category()
        Category of facade finite enumerated sets
        sage: TestSuite(S).run()

    Note that being an enumerated set, the result depends on the order::

        sage: S1 = FiniteEnumeratedSet((1, 2, 3))
        sage: S1
        {1, 2, 3}
        sage: S1.list()
        [1, 2, 3]
        sage: S1 == S
        True
        sage: S2 = FiniteEnumeratedSet((2, 1, 3))
        sage: S2 == S
        False

    As an abuse, repeated entries in ``elements`` are allowed to model
    multisets::

        sage: S1 = FiniteEnumeratedSet((1, 2, 1, 2, 2, 3))
        sage: S1
        {1, 2, 1, 2, 2, 3}

    Finally, the elements are not aware of their parent::

        sage: S.first().parent()
        Integer Ring
    """

    @staticmethod
    def __classcall__(cls, iterable):
        """
        Standard trick to expand the iterable upon input, and
        guarantees unique representation, independently of the type of
        the iterable. See ``UniqueRepresentation``.

        TESTS::

            sage: S1 = FiniteEnumeratedSet([1, 2, 3])
            sage: S2 = FiniteEnumeratedSet((1, 2, 3))
            sage: S3 = FiniteEnumeratedSet((x for x in range(1,4)))
            sage: S1 is S2
            True
            sage: S2 is S3
            True
        """
        return super(FiniteEnumeratedSet, cls).__classcall__(
                cls,
                tuple(iterable))

    def __init__(self, elements):
        """
        TESTS::

            sage: TestSuite(FiniteEnumeratedSet([1,2,3])).run()
            sage: TestSuite(FiniteEnumeratedSet([])).run()
        """
        self._elements = elements
        Parent.__init__(self, facade = True, category = FiniteEnumeratedSets())

    def __bool__(self):
        r"""
        Conversion to boolean.

        EXAMPLES::

            sage: bool(FiniteEnumeratedSet('abc'))
            True
            sage: bool(FiniteEnumeratedSet([]))
            False
        """
        return bool(self._elements)

    __nonzero__ = __bool__

    def _repr_(self):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: repr(S)
            '{1, 2, 3}'
            sage: S = FiniteEnumeratedSet(['1','2','3'])
            sage: repr(S)
            "{'1', '2', '3'}"
            sage: S = FiniteEnumeratedSet([1])
            sage: repr(S)
            '{1}'
        """
        return '{' + ', '.join(repr(e) for e in self._elements) + '}'

    def __contains__(self, x):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: 1 in S
            True
            sage: 2 in S
            True
            sage: 4 in S
            False
            sage: ZZ in S
            False

            sage: S.is_parent_of(2)
            True
            sage: S.is_parent_of(4)
            False
        """
        return x in self._elements

    is_parent_of = __contains__

    def __iter__(self):
        r"""
        Iterator over the element of self.

        EXAMPLES::

            sage: for i in FiniteEnumeratedSet([1,2,3]): print(i)
            1
            2
            3
        """
        return iter(self._elements)

    def list(self):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: S.list()
            [1, 2, 3]
        """
        return list(self._elements)

    def an_element(self):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: S.an_element()
            1
        """
        if not self._elements:
            raise EmptySetError
        return self._elements[0]

    def first(self):
        r"""
        Return the first element of the enumeration or raise an EmptySetError if
        the set is empty.

        EXAMPLES::

            sage: S = FiniteEnumeratedSet('abc')
            sage: S.first()
            'a'
        """
        if not self._elements:
            raise EmptySetError
        return self._elements[0]

    def last(self):
        r"""
        Returns the last element of the iteration or raise an EmptySetError if
        the set is empty.

        EXAMPLES::

            sage: S = FiniteEnumeratedSet([0,'a',1.23, 'd'])
            sage: S.last()
            'd'
        """
        if not self._elements:
            raise EmptySetError
        return self._elements[-1]

    def random_element(self):
        r"""
        Return a random element.

        EXAMPLES::

            sage: S = FiniteEnumeratedSet('abc')
            sage: S.random_element()   # random
            'b'

        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: S.random_element() in S
            True
        """
        if not self._elements:
            raise EmptySetError
        from sage.misc.prandom import choice
        return choice(self._elements)

    def cardinality(self):
        """
        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: S.cardinality()
            3
        """
        return Integer(len(self._elements))

    def rank(self, x):
        """
        Returns the index of ``x`` in this finite enumerated set.

        EXAMPLES::

            sage: S = FiniteEnumeratedSet(['a','b','c'])
            sage: S.index('b')
            1
        """
        return self._elements.index(x)

    index = rank

    def unrank(self,i):
        r"""
        Return the element at position ``i``.

        EXAMPLES::

            sage: S = FiniteEnumeratedSet([1,'a',-51])
            sage: S[0], S[1], S[2]
            (1, 'a', -51)
            sage: S[3]
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
            sage: S[-1], S[-2], S[-3]
            (-51, 'a', 1)
            sage: S[-4]
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
        """
        return self._elements[i]

    def __call__(self, el):
        """
        Coerce or convert ``el`` into an element of ``self``.

        INPUT:

        - ``el`` -- some object

        As :meth:`Parent.__call__`, this tries to convert or coerce
        ``el`` into an element of ``self`` depending on the parent of
        ``el``. If no such conversion or coercion is available, this
        calls :meth:`_element_constructor_`.

        :meth:`Parent.__call__` enforces that
        :meth:`_element_constructor_` return an :class:`Element` (more
        precisely, it calls :meth:`_element_constructor_` through a
        :class:`sage.structure.coerce_maps.DefaultConvertMap_unique`,
        and any :class:`sage.categories.map.Map` requires its results
        to be instances of :class:`Element`).

        Since :class:`FiniteEnumeratedSets` is often a facade over
        plain Python objects, :trac:`16280` introduced this method
        which works around this limitation by calling directly
        :meth:`_element_constructor_` whenever ``el`` is not an
        :class:`Element`. Otherwise :meth:`Parent.__call__` is called
        as usual.

        .. WARNING::

            This workaround prevents conversions or coercions from
            facade parents over plain Python objects into ``self``.

        If the :meth:`Parent.__call__` fails, then we try
        :meth:`_element_constructor_` directly as the element returned
        may not be a subclass of :class:`Element`, which is currently
        not supported (see :trac:`19553`).

        EXAMPLES::

            sage: F = FiniteEnumeratedSet([1, 2, 'a', 'b'])
            sage: F(1)
            1
            sage: F('a')
            'a'

        We check that conversions are properly honored for usual
        parents; this is not the case for facade parents over plain
        Python objects::

            sage: F = FiniteEnumeratedSet([1, 2, 3, 'a', 'aa'])
            sage: phi = Hom(ZZ, F, Sets())(lambda i: i+i)
            sage: phi(1)
            2
            sage: phi.register_as_conversion()

            sage: from sage.sets.pythonclass import Set_PythonType_class
            sage: psi = Hom(Set_PythonType_class(str), F, Sets())(lambda s: ZZ(len(s)))
            sage: psi.register_as_conversion()
            sage: psi('a')
            1
            sage: F(1)
            2
            sage: F('a')
            'a'

        Check that :trac:`19554` is fixed::

            sage: S = FiniteEnumeratedSet(range(5))
            sage: S(1)
            1
            sage: type(S(1))
            <class 'int'>
        """
        if not isinstance(el, Element):
            return self._element_constructor_(el)
        try:
            return Parent.__call__(self, el)
        except TypeError:
            return self._element_constructor_(el)

    def _element_constructor_(self, el):
        """
        Return ``el``.

        INPUT:

        - ``el`` -- an element of ``self``

        If ``el`` is not an element of ``self``, a :class:`ValueError`
        is raised.

        TESTS::

            sage: S = FiniteEnumeratedSet([1,2,3])
            sage: one,two,three = sorted(S)
            sage: S(1) is one
            True
            sage: x = 3
            sage: x is three
            False
            sage: S(x) is three
            True
            sage: S(0)
            Traceback (most recent call last):
            ...
            ValueError: 0 not in {1, 2, 3}
        """
        try:
            return self._elements[self.rank(el)]
        except (ValueError,KeyError):
            raise ValueError("%s not in %s"%(el, self))
