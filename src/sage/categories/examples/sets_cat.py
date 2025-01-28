# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.libs.pari
"""
Examples of sets
"""
# ****************************************************************************
#  Copyright (C) 2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.sets_cat import Sets
from sage.rings.integer import Integer, IntegerWrapper
from sage.rings.integer_ring import IntegerRing
from sage.arith.misc import is_prime
from sage.structure.unique_representation import UniqueRepresentation


class PrimeNumbers(UniqueRepresentation, Parent):
    r"""
    An example of parent in the category of sets: the set of prime numbers.

    The elements are represented as plain integers in `\ZZ` (facade
    implementation).

    This is a minimal implementations. For more advanced examples of
    implementations, see also::

        sage: P = Sets().example("facade")
        sage: P = Sets().example("inherits")
        sage: P = Sets().example("wrapper")

    EXAMPLES::

        sage: P = Sets().example()
        sage: P(12)
        Traceback (most recent call last):
        ...
        AssertionError: 12 is not a prime number
        sage: a = P.an_element()
        sage: a.parent()
        Integer Ring
        sage: x = P(13); x
        13
        sage: type(x)
        <class 'sage.rings.integer.Integer'>
        sage: x.parent()
        Integer Ring
        sage: 13 in P
        True
        sage: 12 in P
        False
        sage: y = x+1; y
        14
        sage: type(y)
        <class 'sage.rings.integer.Integer'>

        sage: TestSuite(P).run(verbose=True)
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
        running ._test_eq() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
    """
    def __init__(self):
        """
        TESTS::

            sage: from sage.categories.examples.sets_cat import PrimeNumbers
            sage: P = PrimeNumbers()
            sage: P.category()
            Category of facade sets
            sage: P is Sets().example()
            True
        """
        Parent.__init__(self, facade=IntegerRing(), category=Sets())

    def _repr_(self):
        """
        TESTS::

            sage: Sets().example()       # indirect doctest
            Set of prime numbers (basic implementation)
        """
        return "Set of prime numbers (basic implementation)"

    def an_element(self):
        """
        Implement :meth:`Sets.ParentMethods.an_element`.

        TESTS::

            sage: P = Sets().example()
            sage: x = P.an_element(); x
            47
            sage: x.parent()
            Integer Ring
        """
        return self(47) # if speed is needed, call: self.element_class(47)

    def __contains__(self, p):
        """
        TESTS::

            sage: P = Sets().example()
            sage: 13 in P
            True
            sage: 12 in P
            False
        """
        return isinstance(p, Integer) and p.is_prime()

    def _element_constructor_(self, e):
        """
        TESTS::

            sage: P = Sets().example()
            sage: P._element_constructor_(13) == 13
            True
            sage: P._element_constructor_(13).parent()
            Integer Ring
            sage: P._element_constructor_(14)
            Traceback (most recent call last):
            ...
            AssertionError: 14 is not a prime number
        """
        p = self.element_class(e)
        assert is_prime(p), "%s is not a prime number" % (p)
        return p

    element_class = Integer


from sage.misc.abstract_method import abstract_method


class PrimeNumbers_Abstract(UniqueRepresentation, Parent):
    """
    This class shows how to write a parent while keeping the choice of the
    datastructure for the children open. Different class with fixed
    datastructure will then be constructed by inheriting from
    :class:`PrimeNumbers_Abstract`.

    This is used by::

        sage: P = Sets().example("facade")
        sage: P = Sets().example("inherits")
        sage: P = Sets().example("wrapper")
    """
    def __init__(self):
        """
        TESTS::

            sage: P = Sets().example("inherits")
        """
        Parent.__init__(self, category=Sets())

    def _repr_(self):
        """
        TESTS::

            sage: Sets().example("inherits")     # indirect doctest
            Set of prime numbers
        """
        return "Set of prime numbers"

    def an_element(self):
        """
        Implement :meth:`Sets.ParentMethods.an_element`.

        TESTS::

            sage: P = Sets().example("inherits")
            sage: x = P.an_element(); x
            47
            sage: x.parent()
            Set of prime numbers
        """
        return self._from_integer_(47)

    def _element_constructor_(self, i):
        """
        Construct an element of ``self`` from an integer, testing that
        this integer is indeed prime.

        EXAMPLES::

            sage: P = Sets().example("inherits")
            sage: P(13)       # indirect doctest
            13
            sage: P(42)
            Traceback (most recent call last):
            ...
            ValueError: 42 is not a prime number
        """
        if i in self:
            return self._from_integer_(i)
        else:
            raise ValueError("%s is not a prime number" % (i))

    @abstract_method
    def _from_integer_(self, i):
        """
        Fast construction of an element of ``self`` from an integer.

        No prime checking is performed. To be defined.

        EXAMPLES::

            sage: P = Sets().example("inherits")
            sage: P._from_integer_(13)
            13
            sage: P._from_integer_(42)       # Do not do that at home kids!
            42
            sage: P(42)
            Traceback (most recent call last):
            ...
            ValueError: 42 is not a prime number
        """

    def next(self, i):
        """
        Return the next prime number.

        EXAMPLES::

            sage: P = Sets().example("inherits")
            sage: x = P.next(P.an_element()); x
            53
            sage: x.parent()
            Set of prime numbers
        """
        assert i in self
        return self._from_integer_((Integer(i) + 1).next_prime())

    def some_elements(self):
        """
        Return some prime numbers.

        EXAMPLES::

            sage: P = Sets().example("inherits")
            sage: P.some_elements()
            [47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
        """
        x = self.an_element()
        res = [x]
        for i in range(10):
            x = self.next(x)
            res.append(x)
        return res

    class Element(Element):
        def is_prime(self):
            """
            Return whether ``self`` is a prime number.

            EXAMPLES::

                sage: P = Sets().example("inherits")
                sage: x = P.an_element()
                sage: P.an_element().is_prime()
                True
            """
            return True

        def next(self):
            """
            Return the next prime number.

            EXAMPLES::

                sage: P = Sets().example("inherits")
                sage: p = P.an_element(); p
                47
                sage: p.next()
                53

            .. NOTE::

                This method is not meant to implement the protocol iterator,
                and thus not subject to Python 2 vs Python 3 incompatibilities.
            """
            return self.parent().next(self)


#*************************************************************************#
class PrimeNumbers_Inherits(PrimeNumbers_Abstract):
    """
    An example of parent in the category of sets: the set of prime numbers.
    In this implementation, the element are stored as object of a new class
    which inherits from the class Integer (technically :class:`IntegerWrapper`).

    EXAMPLES::

        sage: P = Sets().example("inherits")
        sage: P
        Set of prime numbers
        sage: P(12)
        Traceback (most recent call last):
        ...
        ValueError: 12 is not a prime number
        sage: a = P.an_element()
        sage: a.parent()
        Set of prime numbers
        sage: x = P(13); x
        13
        sage: x.is_prime()
        True
        sage: type(x)
        <class 'sage.categories.examples.sets_cat.PrimeNumbers_Inherits_with_category.element_class'>
        sage: x.parent()
        Set of prime numbers
        sage: P(13) in P
        True
        sage: y = x+1; y
        14
        sage: type(y)
        <class 'sage.rings.integer.Integer'>
        sage: y.parent()
        Integer Ring
        sage: type(P(13)+P(17))
        <class 'sage.rings.integer.Integer'>
        sage: type(P(2)+P(3))
        <class 'sage.rings.integer.Integer'>

        sage: z = P.next(x); z
        17
        sage: type(z)
        <class 'sage.categories.examples.sets_cat.PrimeNumbers_Inherits_with_category.element_class'>
        sage: z.parent()
        Set of prime numbers

        sage: TestSuite(P).run(verbose=True)
        running ._test_an_element() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_construction() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass

    See also::

        sage: P = Sets().example("facade")
        sage: P = Sets().example("inherits")
        sage: P = Sets().example("wrapper")
    """

    def __init__(self):
        """
        TESTS::

            sage: P = Sets().example("inherits")
            sage: type(P(13)+P(17))
            <class 'sage.rings.integer.Integer'>
            sage: type(P(2)+P(3))
            <class 'sage.rings.integer.Integer'>
        """
        super().__init__()
        self._populate_coercion_lists_(embedding=IntegerRing())

    def __contains__(self, p):
        """
        TESTS::

            sage: P = Sets().example("inherits")
            sage: 13 in P, P(13) in P
            (True, True)
            sage: 12 in P
            False
        """
        return (isinstance(p, self.element_class) and p.parent() is self
                or isinstance(p, Integer) and p.is_prime())

    def _from_integer_(self, p):
        """
        TESTS::

            sage: P = Sets().example("inherits")
            sage: P._from_integer_(13)
            13
            sage: P._from_integer_(42)            # Don't do that at home kids!
            42
        """
        return self.element_class(self, p)

    class Element(IntegerWrapper, PrimeNumbers_Abstract.Element):
        def __init__(self, parent, p):
            """
            TESTS::

                sage: P = Sets().example("inherits")
                sage: P(12)
                Traceback (most recent call last):
                ...
                ValueError: 12 is not a prime number
                sage: x = P(13); type(x)
                <class 'sage.categories.examples.sets_cat.PrimeNumbers_Inherits_with_category.element_class'>
                sage: x.parent() is P
                True
            """
            IntegerWrapper.__init__(self, parent, p)


#*************************************************************************#
class PrimeNumbers_Wrapper(PrimeNumbers_Abstract):
    """
    An example of parent in the category of sets: the set of prime numbers.

    In this second alternative implementation, the prime integer are stored as
    a attribute of a sage object by inheriting from :class:`ElementWrapper`.  In
    this case we need to ensure conversion and coercion from this parent and
    its element to ``ZZ`` and ``Integer``.

    EXAMPLES::

        sage: P = Sets().example("wrapper")
        sage: P(12)
        Traceback (most recent call last):
        ...
        ValueError: 12 is not a prime number
        sage: a = P.an_element()
        sage: a.parent()
        Set of prime numbers (wrapper implementation)
        sage: x = P(13); x
        13
        sage: type(x)
        <class 'sage.categories.examples.sets_cat.PrimeNumbers_Wrapper_with_category.element_class'>
        sage: x.parent()
        Set of prime numbers (wrapper implementation)
        sage: 13 in P
        True
        sage: 12 in P
        False
        sage: y = x+1; y
        14
        sage: type(y)
        <class 'sage.rings.integer.Integer'>

        sage: z = P.next(x); z
        17
        sage: type(z)
        <class 'sage.categories.examples.sets_cat.PrimeNumbers_Wrapper_with_category.element_class'>
        sage: z.parent()
        Set of prime numbers (wrapper implementation)

    TESTS::

        sage: TestSuite(P).run()
    """
    def __init__(self):
        """
        TESTS::

            sage: P = Sets().example("wrapper")
            sage: P.category()
            Category of sets
            sage: P(13) == 13
            True
            sage: ZZ(P(13)) == 13
            True
            sage: P(13) + 1 == 14
            True
        """
        Parent.__init__(self, category=Sets())
        from sage.rings.integer_ring import IntegerRing
        from sage.categories.homset import Hom
        self.mor = Hom(self, IntegerRing())(lambda z: z.value)
        self._populate_coercion_lists_(embedding=self.mor)

    def _repr_(self):
        """
        TESTS::

            sage: Sets().example("wrapper")       # indirect doctest
            Set of prime numbers (wrapper implementation)
        """
        return "Set of prime numbers (wrapper implementation)"

    def __contains__(self, p):
        """
        TESTS::

            sage: P = Sets().example("wrapper")
            sage: 13 in P
            True
            sage: 12 in P
            False
        """
        return (isinstance(p, self.element_class) and p.parent() == self or
                isinstance(p, Integer) and p.is_prime())

    def _from_integer_(self, e):
        """
        TESTS::

            sage: P = Sets().example("wrapper")
            sage: P._from_integer_(13).parent()
            Set of prime numbers (wrapper implementation)
            sage: P._from_integer_(14)            # Don't do that at home kids!
            14
            sage: P._element_constructor_(14)
            Traceback (most recent call last):
            ...
            ValueError: 14 is not a prime number
        """
        return self.element_class(self, Integer(e))

    from sage.structure.element_wrapper import ElementWrapper

    class Element (ElementWrapper, PrimeNumbers_Abstract.Element):
        def _integer_(self, IntRing):
            """
            Convert to an integer.

            TESTS::

                sage: P = Sets().example("wrapper")
                sage: x = P.an_element()
                sage: Integer(x)           # indirect doctest
                47
            """
            return IntRing(self.value)


#*************************************************************************#
class PrimeNumbers_Facade(PrimeNumbers_Abstract):
    r"""
    An example of parent in the category of sets: the set of prime numbers.

    In this alternative implementation, the elements are represented
    as plain integers in `\ZZ` (facade implementation).

    EXAMPLES::

        sage: P = Sets().example("facade")
        sage: P(12)
        Traceback (most recent call last):
        ...
        ValueError: 12 is not a prime number
        sage: a = P.an_element()
        sage: a.parent()
        Integer Ring
        sage: x = P(13); x
        13
        sage: type(x)
        <class 'sage.rings.integer.Integer'>
        sage: x.parent()
        Integer Ring
        sage: 13 in P
        True
        sage: 12 in P
        False
        sage: y = x+1; y
        14
        sage: type(y)
        <class 'sage.rings.integer.Integer'>

        sage: z = P.next(x); z
        17
        sage: type(z)
        <class 'sage.rings.integer.Integer'>
        sage: z.parent()
        Integer Ring

    The disadvantage of this implementation is that the elements do not know
    that they are prime, so that prime testing is slow::

        sage: pf = Sets().example("facade").an_element()
        sage: timeit("pf.is_prime()") #    random
        625 loops, best of 3: 4.1 us per loop

    compared to the other implementations where prime testing is only done if
    needed during the construction of the element, and later on the elements
    "know" that they are prime::

        sage: pw = Sets().example("wrapper").an_element()
        sage: timeit("pw.is_prime()")    # random
        625 loops, best of 3: 859 ns per loop

        sage: pi = Sets().example("inherits").an_element()
        sage: timeit("pw.is_prime()")    # random
        625 loops, best of 3: 854 ns per loop

    Note also that the ``next`` method for the elements does not exist::

        sage: pf.next()
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no attribute 'next'...

    unlike in the other implementations::

        sage: pw.next()
        53
        sage: pi.next()
        53

    TESTS::

        sage: TestSuite(P).run(verbose = True)
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
        running ._test_eq() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
    """

    def __init__(self):
        """
        TESTS::

            sage: P = Sets().example("inherits")
        """
        Parent.__init__(self, facade=IntegerRing(), category=Sets())

    def _repr_(self):
        """
        TESTS::

            sage: Sets().example("facade")       # indirect doctest
            Set of prime numbers (facade implementation)
        """
        return "Set of prime numbers (facade implementation)"

    def __contains__(self, p):
        """
        TESTS::

            sage: P = Sets().example("facade")
            sage: 13 in P
            True
            sage: 12 in P
            False
        """
        return isinstance(p, Integer) and p.is_prime()

    def _from_integer_(self, e):
        """
        TESTS::

            sage: P = Sets().example("facade")
            sage: P._from_integer_(13).parent()
            Integer Ring
            sage: P._from_integer_(14)            # Don't do that at home kids!
            14
            sage: P._element_constructor_(14)
            Traceback (most recent call last):
            ...
            ValueError: 14 is not a prime number
        """
        return self.element_class(e)

    element_class = Integer
