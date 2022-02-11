"""
Abelian group elements

AUTHORS:

- David Joyner (2006-02); based on free_abelian_monoid_element.py, written by David Kohel.

- David Joyner (2006-05); bug fix in order

- David Joyner (2006-08); bug fix+new method in pow for negatives+fixed corresponding examples.

- David Joyner (2009-02): Fixed bug in order.

- Volker Braun (2012-11) port to new Parent base. Use tuples for immutables.


EXAMPLES:

Recall an example from abelian groups::

    sage: F = AbelianGroup(5,[4,5,5,7,8],names = list("abcde"))
    sage: (a,b,c,d,e) = F.gens()
    sage: x = a*b^2*e*d^20*e^12
    sage: x
    a*b^2*d^6*e^5
    sage: x = a^10*b^12*c^13*d^20*e^12
    sage: x
    a^2*b^2*c^3*d^6*e^4
    sage: y = a^13*b^19*c^23*d^27*e^72
    sage: y
    a*b^4*c^3*d^6
    sage: x*y
    a^3*b*c*d^5*e^4
    sage: x.list()
    [2, 2, 3, 6, 4]
"""

###########################################################################
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Copyright (C) 2006 David Joyner  <wdjoyner@gmail.com>
#  Copyright (C) 2012 Volker Braun  <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
###########################################################################

from sage.groups.abelian_gps.element_base import AbelianGroupElementBase


def is_AbelianGroupElement(x):
    """
    Return true if x is an abelian group element, i.e., an element of
    type AbelianGroupElement.

    EXAMPLES: Though the integer 3 is in the integers, and the integers
    have an abelian group structure, 3 is not an AbelianGroupElement::

        sage: from sage.groups.abelian_gps.abelian_group_element import is_AbelianGroupElement
        sage: is_AbelianGroupElement(3)
        False
        sage: F = AbelianGroup(5, [3,4,5,8,7], 'abcde')
        sage: is_AbelianGroupElement(F.0)
        True
    """
    return isinstance(x, AbelianGroupElement)


class AbelianGroupElement(AbelianGroupElementBase):
    """
    Elements of an
    :class:`~sage.groups.abelian_gps.abelian_group.AbelianGroup`

    INPUT:

    - ``x`` -- list/tuple/iterable of integers (the element vector)

    - ``parent`` -- the parent
      :class:`~sage.groups.abelian_gps.abelian_group.AbelianGroup`

    EXAMPLES::

        sage: F = AbelianGroup(5, [3,4,5,8,7], 'abcde')
        sage: a, b, c, d, e = F.gens()
        sage: a^2 * b^3 * a^2 * b^-4
        a*b^3
        sage: b^-11
        b
        sage: a^-11
        a
        sage: a*b in F
        True
    """
    def as_permutation(self):
        r"""
        Return the element of the permutation group G (isomorphic to the
        abelian group A) associated to a in A.

        EXAMPLES::

            sage: G = AbelianGroup(3,[2,3,4],names="abc"); G
            Multiplicative Abelian group isomorphic to C2 x C3 x C4
            sage: a,b,c=G.gens()
            sage: Gp = G.permutation_group(); Gp
            Permutation Group with generators [(6,7,8,9), (3,4,5), (1,2)]
            sage: a.as_permutation()
            (1,2)
            sage: ap = a.as_permutation(); ap
            (1,2)
            sage: ap in Gp
            True
        """
        from sage.libs.gap.libgap import libgap
        G = self.parent()
        A = libgap.AbelianGroup(G.gens_orders())
        phi = libgap.IsomorphismPermGroup(A)
        gens = libgap.GeneratorsOfGroup(A)
        L2 = libgap.Product([geni**Li for geni, Li in zip(gens, self.list())])
        pg = libgap.Image(phi, L2)
        return G.permutation_group()(pg)

    def word_problem(self, words):
        """
        TODO - this needs a rewrite - see stuff in the matrix_grp
        directory.

        G and H are abelian groups, g in G, H is a subgroup of G generated
        by a list (words) of elements of G. If self is in H, return the
        expression for self as a word in the elements of (words).

        This function does not solve the word problem in Sage. Rather
        it pushes it over to GAP, which has optimized (non-deterministic)
        algorithms for the word problem.

        .. warning::

            Don't use E (or other GAP-reserved letters) as a generator
            name.

        EXAMPLES::

            sage: G = AbelianGroup(2,[2,3], names="xy")
            sage: x,y = G.gens()
            sage: x.word_problem([x,y])
            [[x, 1]]
            sage: y.word_problem([x,y])
            [[y, 1]]
            sage: v = (y*x).word_problem([x,y]); v #random
            [[x, 1], [y, 1]]
            sage: prod([x^i for x,i in v]) == y*x
            True
        """
        from sage.groups.abelian_gps.abelian_group import word_problem
        return word_problem(words, self)
