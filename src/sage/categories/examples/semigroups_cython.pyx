# sage_setup: distribution = sagemath-categories
r"""
Examples of semigroups in cython
"""

from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.categories.category import Category
from sage.categories.semigroups import Semigroups
from sage.categories.examples.semigroups import LeftZeroSemigroup as LeftZeroSemigroupPython
from cpython.object cimport PyObject_RichCompare


class IdempotentSemigroups(Category):
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import IdempotentSemigroups
            sage: IdempotentSemigroups().super_categories()
            [Category of semigroups]
        """
        return [Semigroups()]

    class ElementMethods:
        def _pow_int(self, i):
            """
            EXAMPLES::

                sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
                sage: S = LeftZeroSemigroup()
                sage: S(2)._pow_int(3)
                2
            """
            assert i > 0
            return self

        def is_idempotent(self):
            """
            EXAMPLES::

                sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
                sage: S = LeftZeroSemigroup()
                sage: S(2).is_idempotent()
                True
            """
            return True


cdef class LeftZeroSemigroupElement(Element):
    cdef object _value

    def __init__(self, parent, value):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: x = S(3)
            sage: TestSuite(x).run()
        """
        Element.__init__(self, parent=parent)
        self._value = value

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S(3)                    # indirect doctest
            3
        """
        return repr(self._value)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: x = S(3)
            sage: x.__reduce__()
            (<class 'sage.categories.examples.semigroups_cython.LeftZeroSemigroupElement'>,
             (An example of a semigroup: the left zero semigroup, 3))
        """
        return LeftZeroSemigroupElement, (self._parent, self._value)

    cpdef _richcmp_(self, other, int op):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S(3) == S(3)
            True
            sage: S(3) == S(2)
            False
        """
        left = (<LeftZeroSemigroupElement>self)._value
        right = (<LeftZeroSemigroupElement>other)._value
        return PyObject_RichCompare(left, right, op)

    cpdef _mul_(self, other):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S(2) * S(3)
            2
            sage: S(2)._mul_(S(3))
            2
        """
        return self.parent().product(self, other)

    def __pow__(self, i, dummy):
        """
        EXAMPLES::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S(2)^3
            2
        """
        return self._pow_int(i)


class LeftZeroSemigroup(LeftZeroSemigroupPython):
    r"""
    An example of semigroup.

    This class illustrates a minimal implementation of a semi-group
    where the element class is an extension type, and still gets code
    from the category. The category itself must be a Python class
    though.

    This is purely a proof of concept. The code obviously needs refactorisation!

    Comments:

    - one cannot play ugly class surgery tricks (as with _mul_parent).
      available operations should really be declared to the coercion model!

    EXAMPLES::

        sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
        sage: S = LeftZeroSemigroup(); S
        An example of a semigroup: the left zero semigroup

    This is the semigroup which contains all sort of objects::

        sage: S.some_elements()
        [3, 42, 'a', 3.4, 'raton laveur']

    with product rule given by `a \times b = a` for all `a,b`. ::

        sage: S('hello') * S('world')
        'hello'

        sage: S(3)*S(1)*S(2)
        3

        sage: S(3)^12312321312321
        3

        sage: TestSuite(S).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
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

    That's really the only method which is obtained from the category ... ::

        sage: S(42).is_idempotent
        <bound method IdempotentSemigroups.ElementMethods.is_idempotent of 42>
        sage: S(42).is_idempotent()
        True

        sage: S(42)._pow_int
        <bound method IdempotentSemigroups.ElementMethods._pow_int of 42>
        sage: S(42)^10
        42

        sage: S(42).is_idempotent
        <bound method IdempotentSemigroups.ElementMethods.is_idempotent of 42>
        sage: S(42).is_idempotent()
        True
    """

    def __init__(self):
        """
        TESTS::

            sage: from sage.categories.examples.semigroups_cython import LeftZeroSemigroup
            sage: S = LeftZeroSemigroup()
            sage: S.category()
            Category of idempotent semigroups
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category=IdempotentSemigroups())

    Element = LeftZeroSemigroupElement
