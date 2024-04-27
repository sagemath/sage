r"""
Example of a set with grading
"""
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_with_grading import SetsWithGrading

from sage.rings.integer_ring import IntegerRing
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet


class NonNegativeIntegers(UniqueRepresentation, Parent):
    r"""
    Non negative integers graded by themselves.

    EXAMPLES::

        sage: E = SetsWithGrading().example(); E
        Non negative integers
        sage: E in Sets().Infinite()
        True
        sage: E.graded_component(0)
        {0}
        sage: E.graded_component(100)
        {100}
    """
    def __init__(self):
        r"""
        TESTS::

            sage: TestSuite(SetsWithGrading().example()).run()
        """
        Parent.__init__(self, category=SetsWithGrading().Infinite(),
                        facade=IntegerRing())

    def an_element(self):
        r"""
        Return 0.

        EXAMPLES::

            sage: SetsWithGrading().example().an_element()
            0
        """
        return 0

    def _repr_(self):
        r"""
        TESTS::

            sage: SetsWithGrading().example() # indirect example
            Non negative integers
        """
        return "Non negative integers"

    def graded_component(self, grade):
        r"""
        Return the component with grade ``grade``.

        EXAMPLES::

            sage: N = SetsWithGrading().example()
            sage: N.graded_component(65)
            {65}
        """
        return FiniteEnumeratedSet([grade])

    def grading(self, elt):
        r"""
        Return the grade of ``elt``.

        EXAMPLES::

            sage: N = SetsWithGrading().example()
            sage: N.grading(10)
            10
        """
        return elt

    def generating_series(self, var='z'):
        r"""
        Return `1 / (1-z)`.

        EXAMPLES::

            sage: N = SetsWithGrading().example(); N
            Non negative integers
            sage: f = N.generating_series(); f
            1/(-z + 1)
            sage: LaurentSeriesRing(ZZ,'z')(f)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + z^7 + z^8 + z^9 + z^10 + z^11 + z^12 + z^13 + z^14 + z^15 + z^16 + z^17 + z^18 + z^19 + O(z^20)
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.integer import Integer
        R = PolynomialRing(IntegerRing(), var)
        z = R.gen()
        return Integer(1) / (Integer(1) - z)


Example = NonNegativeIntegers
