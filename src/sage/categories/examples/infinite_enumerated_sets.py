# sage_setup: distribution = sagemath-categories
"""
Examples of infinite enumerated sets
"""
# ****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent import Parent
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer import Integer


class NonNegativeIntegers(UniqueRepresentation, Parent):
    r"""
    An example of infinite enumerated set: the nonnegative integers.

    This class provides a minimal implementation of an infinite enumerated set.

    EXAMPLES::

        sage: NN = InfiniteEnumeratedSets().example()
        sage: NN
        An example of an infinite enumerated set: the nonnegative integers
        sage: NN.cardinality()
        +Infinity
        sage: NN.list()
        Traceback (most recent call last):
        ...
        NotImplementedError: cannot list an infinite set
        sage: NN.element_class
        <class 'sage.rings.integer.Integer'>
        sage: it = iter(NN)
        sage: [next(it), next(it), next(it), next(it), next(it)]
        [0, 1, 2, 3, 4]
        sage: x = next(it); type(x)
        <class 'sage.rings.integer.Integer'>
        sage: x.parent()
        Integer Ring
        sage: x+3
        8
        sage: NN(15)
        15
        sage: NN.first()
        0

    This checks that the different methods of `NN` return consistent
    results::

        sage: TestSuite(NN).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_construction() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_nonzero_equal() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
    """

    def __init__(self):
        """
        TESTS::

            sage: NN = InfiniteEnumeratedSets().example()
            sage: NN
            An example of an infinite enumerated set: the nonnegative integers
            sage: NN.category()
            Category of infinite enumerated sets
            sage: TestSuite(NN).run()
        """
        Parent.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: InfiniteEnumeratedSets().example() # indirect doctest
            An example of an infinite enumerated set: the nonnegative integers
        """
        return "An example of an infinite enumerated set: the nonnegative integers"

    def __contains__(self, elt):
        """
        EXAMPLES::

            sage: NN = InfiniteEnumeratedSets().example()
            sage: 1 in NN
            True
            sage: -1 in NN
            False
        """
        return Integer(elt) >= Integer(0)

    def __iter__(self):
        """
        EXAMPLES::

            sage: NN = InfiniteEnumeratedSets().example()
            sage: g = iter(NN)
            sage: next(g), next(g), next(g), next(g)
            (0, 1, 2, 3)
        """
        i = Integer(0)
        while True:
            yield self._element_constructor_(i)
            i += 1

    def __call__(self, elt):
        """
        EXAMPLES::

            sage: NN = InfiniteEnumeratedSets().example()
            sage: NN(3)         # indirect doctest
            3
            sage: NN(3).parent()
            Integer Ring
            sage: NN(-1)
            Traceback (most recent call last):
            ...
            ValueError: Value -1 is not a nonnegative integer.
        """
        if elt in self:
            return self._element_constructor_(elt)
        raise ValueError("Value %s is not a nonnegative integer." % (elt))

    def an_element(self):
        """
        EXAMPLES::

            sage: InfiniteEnumeratedSets().example().an_element()
            42
        """
        return self._element_constructor_(Integer(42))

    def next(self, o):
        """
        EXAMPLES::

            sage: NN = InfiniteEnumeratedSets().example()
            sage: NN.next(3)
            4
        """
        return self._element_constructor_(o+1)

    def _element_constructor_(self, i):
        """
        The default implementation of _element_constructor_ assumes
        that the constructor of the element class takes the parent as
        parameter. This is not the case for ``Integer``, so we need to
        provide an implementation.


        TESTS::

            sage: NN = InfiniteEnumeratedSets().example()
            sage: x = NN(42); x
            42
            sage: type(x)
            <class 'sage.rings.integer.Integer'>
            sage: x.parent()
            Integer Ring
        """
        return self.element_class(i)

    Element = Integer


Example = NonNegativeIntegers
