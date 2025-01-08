r"""
Examples of set operads
"""
# ****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.categories.all import SetOperads
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent


class AssociativeOperad(UniqueRepresentation, Parent):
    r"""
    An example of a set operad: the Associative operad

    This class illustrates a minimal implementation of a set operad.
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: A = SetOperads().example(); A
            An example of a set operad: the Associative operad
            sage: TestSuite(A).run()
        """
        Parent.__init__(self, category=SetOperads())

    def _repr_(self):
        """
        Return the string representation.

        EXAMPLES::

            sage: SetOperads().example()  # indirect doctest
            An example of a set operad: the Associative operad
        """
        return "An example of a set operad: the Associative operad"

    @cached_method
    def one(self, letter):
        """
        Return the word of length one, which index the one of this operad.

        EXAMPLES::

            sage: A = SetOperads().example()
            sage: A.one("a")
            'a'
        """
        return self(letter)

    def an_element(self):
        """
        Return a word.

        EXAMPLES::

            sage: A = SetOperads().example()
            sage: A.an_element()
            'abcd'
        """
        return self('abcd')

    def some_elements(self):
        """
        Return a few words.

        EXAMPLES::

            sage: A = SetOperads().example()
            sage: A.some_elements()
            ['abcd', 'aa', '', 'baba']
        """
        return [self(u) for u in ['abcd', 'aa', '', 'baba']]

    def composition(self, x, y, i):
        """
        Insert a word y at position i in a word x and return a word.

        This is the composition of the set-theoretic Associative operad.

        EXAMPLES::

            sage: A = SetOperads().example()
            sage: A.composition(A("acb"), A("de"),"c")
            'adeb'
        """
        pos = x.value.index(i)
        return self(x.value[:pos] + y.value + x.value[pos + 1:])

    def composition_with_numbers(self, x, y, i):
        """
        Insert a word y at position i in a word x and return a word.

        This is the composition of the set-theoretic Associative operad.

        INPUT:

        - `i` -- an integer

        EXAMPLES::

            sage: A = SetOperads().example()
            sage: A.composition_with_numbers(A("123"), A("21"),1)
            '2134'
        """
        pos = x.value.index(str(i))
        n = len(y.value)
        k = int(x.value[pos])

        def shift(z):
            if z < k:
                return str(z)
            else:
                return str(z + n - 1)
        newx = ''.join(shift(int(v)) for v in x.value)
        newy = ''.join(str(int(v) + k - 1) for v in y.value)
        return self(newx[:pos] + newy + newx[pos + 1:])

    def assoc_product(self, a, b):
        """
        Return the associative product.

        Here, this is just concatenation.

        EXAMPLES::

            sage: A = SetOperads().example()
            sage: A.assoc_product(A("123"), A("45"))
            '12345'
        """
        return self(a.value + b.value)

    def operad_morphism(self, t, cod):
        """
        Define a morphism from the Associative set operad to the target operad.

        The target operad has to possess a method called ``assoc_product``.

        The argument t (an element of A) should not have repeated labels.

        EXAMPLES::

            sage: A = SetOperads().example()
            sage: e = A('a')
            sage: A.operad_morphism(e, A) == e
            True
            sage: x = A('abc')
            sage: A.operad_morphism(x, A) == x
            True
        """
        try:
            targetAssocProduct = cod.assoc_product
        except AttributeError:
            raise
        if len(t.value) == 1:
            return cod.one(t.value[0])
        else:
            return targetAssocProduct(cod.one(t.value[0]),
                                      self.operad_morphism(self(t.value[1:]),
                                                           cod))

    class Element(ElementWrapper):
        wrapped_class = str

        def degree(self):
            """
            Return the degree of an element.

            EXAMPLES::

                sage: A = SetOperads().example()
                sage: y = A('12')
                sage: y.degree()
                2
            """
            return len(self.value)

        def map_labels(self, f):
            """
            Apply the function `f` to the labels of an element.

            EXAMPLES::

                sage: A = SetOperads().example()
                sage: y = A('Z3')
                sage: y.map_labels(lambda v:v.lower())
                'z3'
            """
            return self.parent()(''.join(str(f(v)) for v in self.value))


Example = AssociativeOperad
