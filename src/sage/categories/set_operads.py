r"""
Set Operads
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from sage.categories.sets_cat import Sets
from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from sage.categories.category_singleton import Category_singleton


class SetOperads(Category_singleton):
    """
    The category of set operads.

    EXAMPLES::

        sage: SetOperads()
        Category of set operads
        sage: SetOperads().super_categories()
        [Category of sets]

    TESTS::

        sage: C = SetOperads()
        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        """
        Return the super-categories of ``self``.

        EXAMPLES::

            sage: SetOperads().super_categories()
            [Category of sets]
        """
        return [Sets()]

    def example(self):
        """
        Return an example of set operad.

        Here, the associative operad.

        EXAMPLES::

            sage: SetOperads().example()
            An example of a set operad: the Associative operad
        """
        from sage.categories.examples.set_operads import Example
        return Example()

    class ParentMethods:

        @abstract_method(optional=True)
        def composition(self, left, right, index):
            """
            Return the composition of ``left`` with ``right`` at position
            ``index``.

            EXAMPLES::

                sage: P = SetOperads().example()
                sage: x = P.one('i')
                sage: y = P.one('j')
                sage: P.composition(x, y, 'i')
                'j'
            """

        @abstract_method(optional=True)
        def composition_with_numbers(self):
            """
            This is a variant of composition, where one assumes that
            the objects are labelled by integers from `1` to `n`. The
            result is labelled in the same way.

            EXAMPLES::

                sage: A = SetOperads().example()
                sage: A.composition_with_numbers(A('4321'),A('123'),1)
                '654123'
            """

        def global_composition_with_numbers(self, left, list_right):
            r"""
            Return the global composition of ``left`` with a list of elements.

            The elements are supposed to be labelled by consecutive integers.

            EXAMPLES::

                sage: A = SetOperads().example()
                sage: A.global_composition_with_numbers(A('21'),[A('123'),A('312')])
                '645123'
            """
            if self.composition_with_numbers is not NotImplemented:
                if left.degree() != len(list_right):
                    raise ValueError("the degree of x is not equal to the length of list_right")
                res = left
                for i in range(left.degree(), 0, -1):
                    res = self.composition_with_numbers(res,
                                                        list_right[i - 1], i)
                return res
            else:
                return NotImplemented

        @abstract_method(optional=True)
        def operad_morphism(self, arg, codomain):
            """
            Return the image of ``arg`` by a morphism from ``self`` to
            ``codomain``.

            EXAMPLES::

                sage: A = SetOperads().example()
                sage: e = A('a')
                sage: A.operad_morphism(e, A) == e
                True
                sage: x = A('abc')
                sage: A.operad_morphism(x, A) == x
                True
            """

        @abstract_method(optional=True)
        def one(self, letter):
            """
            Return the one of the operad.

            INPUT:

            ``letter`` -- the chosen labelling.

            EXAMPLES::

                sage: A = SetOperads().example()
                sage: A.one('x')
                'x'
            """

        def cardinality(self):
            """
            Return the cardinality.

            This is usually infinity.

            EXAMPLES::

                sage: A = SetOperads().example()
                sage: A.cardinality()
                +Infinity
            """
            from sage.rings.infinity import Infinity
            return Infinity

    class ElementMethods:

        def compose(self, other, index):
            """
            Return the composition of ``self`` with ``other`` at position
            ``index``.

            EXAMPLES::

                sage: A = SetOperads().example()
                sage: x = A('fou')
                sage: y = A('pi')
                sage: x.compose(y,'u')
                'fopi'
            """
            return self.parent().composition(self, other, index)

        def compose_with_numbers(self, other, index):
            """
            Return the composition of ``self`` with ``other`` at position
            ``index``.

            The elements are supposed to be labelled by consecutive integers.

            EXAMPLES::

                sage: A = SetOperads().example()
                sage: x = A('123')
                sage: y = A('12')
                sage: x.compose_with_numbers(y,3)
                '1234'
            """
            return self.parent().composition_with_numbers(self, other, index)

        @abstract_method(optional=True)
        def degree(self):
            """
            Return the degree of an element.

            EXAMPLES::

                sage: A = SetOperads().example()
                sage: y = A('12')
                sage: y.degree()
                2
            """
            pass

        @abstract_method(optional=True)
        def map_labels(self, f):
            """
            Apply the function `f` to the labels of an element.

            EXAMPLES::

                sage: A = SetOperads().example()
                sage: y = A('R2D2')
                sage: y.map_labels(lambda u: u.lower())
                'r2d2'
            """
            pass
