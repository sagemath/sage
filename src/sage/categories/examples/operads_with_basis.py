r"""
Examples of operads with basis
"""
# ****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.cachefunc import cached_method
from sage.categories.operads_with_basis import OperadsWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.words import Words


class AssociativeOperad(CombinatorialFreeModule):
    r"""
    An example of an operad with basis: the Associative operad

    This class illustrates a minimal implementation of an operad with basis.
    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = OperadsWithBasis(QQ).example(); A
            An example of an operad with basis: the Associative operad over Rational Field
            sage: TestSuite(A).run()
        """
        CombinatorialFreeModule.__init__(self, R, Words(),
                                         category=OperadsWithBasis(R))

    def _repr_(self):
        """
        EXAMPLES::

            sage: OperadsWithBasis(QQ).example() # indirect doctest
            An example of an operad with basis: the Associative operad
            over Rational Field
        """
        return "An example of an operad with basis: the Associative operad over %s" % (self.base_ring())

    @cached_method
    def one_basis(self, letter):
        """
        Return the word of length one, which index the one of this operad,
        as per :meth:`OperadsWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = OperadsWithBasis(QQ).example()
            sage: A.one_basis("a")
            word: a
        """
        return self.basis().keys()([letter])

    def degree_on_basis(self, t):
        """
        Return the degree of a word `t` in the Associative operad.

        This is the length of the word.

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.degree_on_basis(m)
            4
        """
        return t.length()

    def map_labels(self, t, f):
        """
        Map the function `f` on the word `t`.

        INPUT:

        - t -- the index of a basis element

        - f -- a map that can be applied to the labels of t

        EXAMPLES::

            sage: A = OperadsWithBasis(QQ).example()
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.map_labels(m,lambda u:u)
            word: 4321
        """
        return self.basis().keys()([f(u) for u in t])

    def labelling_on_basis(self, t):
        """
        Put canonical labels on a word in the Associative operad.

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.labelling_on_basis(m)
            B[word: 1234]
        """
        B = self.basis()
        return B[B.keys()([1 + i for i in range(t.length())])]

    def unlabelling_on_basis(self, t):
        """
        Remove the labels of a tree in the Associative operad.

        EXAMPLES::

            sage: A = AssociativeOperad(QQ)
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.unlabelling_on_basis(m)
            B[word: 1111]
        """
        B = self.basis()
        return B[B.keys()([1 for i in range(t.length())])]

    def grafts(self, x, y, i):
        """
        Insert a word y at position i in a word x and return a word

        This is the composition of the set-theoretic Associative operad.

        EXAMPLES::

            sage: A = OperadsWithBasis(QQ).example()
            sage: Words = A.basis().keys()
            sage: A.grafts(Words("acb"), Words("de"),"c")
            word: adeb
        """
        if x[0] == i:
            return y + x[1:]
        return x[:1] + self.grafts(x[1:], y, i)

    def composition_on_basis(self, x, y, i):
        """
        Composition of basis elements, as per :meth:`OperadsWithBasis.ParentMethods.composition_on_basis`.

        EXAMPLES::

            sage: A = OperadsWithBasis(QQ).example()
            sage: Words = A.basis().keys()
            sage: A.composition_on_basis(Words("acb"), Words("de"),"c")
            B[word: adeb]
        """
        if i not in x:
            raise ValueError("the composition index is not present")
        return self.basis()[self.grafts(x, y, i)]


Example = AssociativeOperad
