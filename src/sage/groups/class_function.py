r"""
Class functions of groups.

This module implements a wrapper of GAP's ClassFunction function.

NOTE: The ordering of the columns of the character table of a group
corresponds to the ordering of the list. However, in general there is
no way to canonically list (or index) the conjugacy classes of a group.
Therefore the ordering of  the columns of the character table of
a group is somewhat random.

AUTHORS:

- Franco Saliola (November 2008): initial version

- Volker Braun (October 2010): Bugfixes, exterior and symmetric power.
"""

#*****************************************************************************
#       Copyright (C) 2008 Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.richcmp import richcmp, richcmp_method
from sage.interfaces.gap import gap
from sage.rings.integer import Integer
from sage.rings.all import CyclotomicField
from sage.libs.gap.element import GapElement
from sage.libs.gap.libgap import libgap
from sage.libs.gap.element import GapElement as LibGapElement

# TODO:
#
# This module needs to be rewritten to implement the ring of class
# functions in the usual parent/element pattern. But
# http://trac.sagemath.org/14014 is already too long...


def ClassFunction(group, values):
    """
    Construct a class function.

    INPUT:

    - ``group`` -- a group.

    - ``values`` -- list/tuple/iterable of numbers. The values of the
      class function on the conjugacy classes, in that order.


    EXAMPLES::

        sage: G = CyclicPermutationGroup(4)
        sage: G.conjugacy_classes()
        [Conjugacy class of () in Cyclic group of order 4 as a permutation group,
         Conjugacy class of (1,2,3,4) in Cyclic group of order 4 as a permutation group,
         Conjugacy class of (1,3)(2,4) in Cyclic group of order 4 as a permutation group,
         Conjugacy class of (1,4,3,2) in Cyclic group of order 4 as a permutation group]
        sage: values  = [1, -1, 1, -1]
        sage: chi = ClassFunction(G, values); chi
        Character of Cyclic group of order 4 as a permutation group
    """
    try:
        return group.class_function(values)
    except AttributeError:
        pass

    if isinstance(values, LibGapElement):
        return ClassFunction_libgap(group, values)

    return ClassFunction_gap(group, values)



#####################################################################
###
### GAP Interface-based Class Function
###
### This is old code that should be deleted once we have transitioned
### everything to using the library interface to GAP.
###
#####################################################################

@richcmp_method
class ClassFunction_gap(SageObject):
    """
    A wrapper of GAP's ClassFunction function.

    .. NOTE::

        It is *not* checked whether the given values describes a character,
        since GAP does not do this.

    EXAMPLES::

        sage: G = CyclicPermutationGroup(4)
        sage: values  = [1, -1, 1, -1]
        sage: chi = ClassFunction(G, values); chi
        Character of Cyclic group of order 4 as a permutation group
        sage: loads(dumps(chi)) == chi
        True
    """

    def __init__(self, G, values):
        r"""
        Return the character of the group ``G`` with values given by the list
        values. The order of the values must correspond to the output of
        ``G.conjugacy_classes_representatives()``.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(4)
            sage: values  = [1, -1, 1, -1]
            sage: chi = ClassFunction(G, values); chi
            Character of Cyclic group of order 4 as a permutation group
        """
        self._group = G
        if isinstance(values, GapElement) and gap.IsClassFunction(values):
            self._gap_classfunction = values
        else:
            self._gap_classfunction = gap.ClassFunction(G, list(values))

        e = self._gap_classfunction.Conductor()
        self._base_ring = CyclotomicField(e)

    def _gap_init_(self):
        r"""
        Returns a string showing how to declare / initialize self in Gap.
        Stored in the \code{self._gap_string} attribute.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(4)
            sage: values  = [1, -1, 1, -1]
            sage: ClassFunction(G, values)._gap_init_()
            'ClassFunction( CharacterTable( Group( [ (1,2,3,4) ] ) ), [ 1, -1, 1, -1 ] )'
        """
        return str(self._gap_classfunction)

    def _gap_(self, *args):
        r"""
        Coerce self into a GAP element.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(4)
            sage: values  = [1, -1, 1, -1]
            sage: chi = ClassFunction(G, values);  chi
            Character of Cyclic group of order 4 as a permutation group
            sage: type(_)
            <class 'sage.groups.class_function.ClassFunction_gap'>
            sage: chi._gap_()
            ClassFunction( CharacterTable( Group( [ (1,2,3,4) ] ) ), [ 1, -1, 1, -1 ] )
            sage: type(_)
            <class 'sage.interfaces.gap.GapElement'>
        """
        return self._gap_classfunction

    def __repr__(self):
        r"""
        Return a string representation.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: values  = [1, -1, 1, 1, -1]
            sage: ClassFunction(G, values)
            Character of Symmetric group of order 4! as a permutation group
        """
        return "Character of %s" % repr(self._group)

    def __iter__(self):
        r"""
        Iterate through the values of self evaluated on the conjugacy
        classes.

        EXAMPLES::

            sage: xi = ClassFunction(SymmetricGroup(4), [1, -1, 1, 1, -1])
            sage: list(xi)
            [1, -1, 1, 1, -1]
        """
        for v in self._gap_classfunction:
            yield self._base_ring(v)

    def __richcmp__(self, other, op):
        r"""
        Rich comparison for class functions.

        Compares groups and then the values of the class function on the
        conjugacy classes.

        EXAMPLES::

            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: chi = G.character([1, 1, 1, 1, 1, 1, 1])
            sage: H = PermutationGroup([[(1,2,3),(4,5)]])
            sage: xi = H.character([1, 1, 1, 1, 1, 1])
            sage: chi == chi
            True
            sage: xi == xi
            True
            sage: xi == chi
            False
            sage: chi < xi
            False
            sage: xi < chi
            True
        """
        if isinstance(other, ClassFunction_gap):
            return richcmp((self._group, self.values()),
                           (other._group, other.values()), op)
        else:
            return NotImplemented

    def __hash__(self):
        r"""
        TESTS::

            sage: G = SymmetricGroup(5)
            sage: chi1 = ClassFunction(G,[1,1,1,1,1,1,1])
            sage: d = {chi1:'trivial'}
        """
        return hash((self._group, tuple(self)))

    def __reduce__(self):
        r"""
        Add pickle support.

        EXAMPLES::

            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: chi = ClassFunction(G, [1, 1, 1, 1, 1, 1, 1])
            sage: type(chi)
            <class 'sage.groups.class_function.ClassFunction_gap'>
            sage: loads(dumps(chi)) == chi
            True
        """
        return ClassFunction_gap, (self._group, self.values())


    def domain(self):
        r"""
        Returns the domain of the self.

        OUTPUT:

        The underlying group of the class function.

        EXAMPLES::

            sage: ClassFunction(SymmetricGroup(4), [1,-1,1,1,-1]).domain()
            Symmetric group of order 4! as a permutation group
        """
        return self._group


    def __call__(self, g):
        """
        Evaluate the character on the group element `g`.

        Return an error if `g` is not in `G`.

        EXAMPLES::

            sage: G = GL(2,7)
            sage: values = G.gap().CharacterTable().Irr()[2].List().sage()
            sage: chi = ClassFunction(G, values)
            sage: z = G([[3,0],[0,3]]); z
            [3 0]
            [0 3]
            sage: chi(z)
            zeta3
            sage: G = GL(2,3)
            sage: chi = G.irreducible_characters()[3]
            sage: g = G.conjugacy_classes_representatives()[6]
            sage: chi(g)
            zeta8^3 + zeta8

            sage: G = SymmetricGroup(3)
            sage: h = G((2,3))
            sage: triv = G.trivial_character()
            sage: triv(h)
            1
        """
        return self._base_ring(gap(g)._operation("^", self._gap_classfunction))


    def __add__(self, other):
        r"""
        Returns the sum of the characters self and other.

        INPUT:

        - ``other`` -- a :class:`ClassFunction` of the same group as
          ``self``.

        OUTPUT:

        A :class:`ClassFunction`

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: s = chi+chi
            sage: s
            Character of Symmetric group of order 4! as a permutation group
            sage: s.values()
            [6, 2, -2, 0, -2]
        """
        if not isinstance(other, ClassFunction_gap):
            raise NotImplementedError
        s = self._gap_classfunction + other._gap_classfunction
        return ClassFunction(self._group, s)


    def __sub__(self, other):
        r"""
        Returns the difference of the characters ``self`` and ``other``.

        INPUT:

        - ``other`` -- a :class:`ClassFunction` of the same group as
          ``self``.

        OUTPUT:

        A :class:`ClassFunction`

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: chi1 = ClassFunction(G, [3, 1, -1, 0, -1])
            sage: chi2 = ClassFunction(G, [1, -1, 1, 1, -1])
            sage: s = chi1 - chi2
            sage: s
            Character of Symmetric group of order 4! as a permutation group
            sage: s.values()
            [2, 2, -2, -1, 0]
        """
        if not isinstance(other, ClassFunction_gap):
            raise NotImplementedError
        s = self._gap_classfunction - other._gap_classfunction
        return ClassFunction(self._group, s)


    def __mul__(self, other):
        r"""
        Return the product of the character with ``other``.

        INPUT:

        - ``other`` -- either a number or a :class:`ClassFunction` of
          the same group as ``self``. A number can be anything that
          can be converted into GAP: integers, rational, and elements
          of certain number fields.

        OUTPUT:

        A :class:`ClassFunction`

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: chi1 = ClassFunction(G, [3, 1, -1, 0, -1])
            sage: 3*chi1
            Character of Symmetric group of order 4! as a permutation group
            sage: 3*chi1 == chi1+chi1+chi1
            True
            sage: (3*chi1).values()
            [9, 3, -3, 0, -3]


            sage: (1/2*chi1).values()
            [3/2, 1/2, -1/2, 0, -1/2]

            sage: CF3 = CyclotomicField(3)
            sage: CF3.inject_variables()
            Defining zeta3
            sage: (zeta3 * chi1).values()
            [3*zeta3, zeta3, -zeta3, 0, -zeta3]

            sage: chi2 = ClassFunction(G, [1, -1, 1, 1, -1])
            sage: p = chi1*chi2
            sage: p
            Character of Symmetric group of order 4! as a permutation group
            sage: p.values()
            [3, -1, -1, 0, 1]
        """
        if isinstance(other, ClassFunction_gap):
            p = self._gap_classfunction * other._gap_classfunction
            return ClassFunction(self._group, p)
        else:
            return ClassFunction(self._group, other * self._gap_classfunction)


    def __rmul__(self, other):
        r"""
        Return the reverse multiplication of ``self`` and ``other``.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: chi = ClassFunction(G, [3, 1, -1, 0, -1])
            sage: chi * 4   # calls chi.__mul__
            Character of Symmetric group of order 4! as a permutation group
            sage: 4 * chi   # calls chi.__rmul__
            Character of Symmetric group of order 4! as a permutation group
            sage: (4 * chi).values()
            [12, 4, -4, 0, -4]
        """
        return self * other

    def __pos__(self):
        r"""
        Return ``self``.

        OUTPUT:

        A :class:`ClassFunction`

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: +chi
            Character of Symmetric group of order 4! as a permutation group
            sage: _.values()
            [3, 1, -1, 0, -1]
            sage: chi.__pos__() == +chi
            True
        """
        return ClassFunction(self._group, self._gap_classfunction)


    def __neg__(self):
        r"""
        Return the additive inverse of ``self``.

        OUTPUT:

        A :class:`ClassFunction`

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: -chi
            Character of Symmetric group of order 4! as a permutation group
            sage: _.values()
            [-3, -1, 1, 0, 1]
            sage: chi.__neg__() == -chi
            True
        """
        return ClassFunction(self._group, -self._gap_classfunction)


    def __pow__(self, other):
        r"""
        Returns the product of self with itself other times.

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: p = chi**3
            sage: p
            Character of Symmetric group of order 4! as a permutation group
            sage: p.values()
            [27, 1, -1, 0, -1]
        """
        if not isinstance(other, (int,Integer)):
            raise NotImplementedError
        return ClassFunction(self._group, self._gap_classfunction ** other)


    def symmetric_power(self, n):
        r"""
        Returns the symmetrized product of self with itself ``n`` times.

        INPUT:

        - ``n`` -- a positive integer.

        OUTPUT:

        The ``n``-th symmetrized power of ``self`` as a
        :class:`ClassFunction`.

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: p = chi.symmetric_power(3)
            sage: p
            Character of Symmetric group of order 4! as a permutation group
            sage: p.values()
            [10, 2, -2, 1, 0]
        """
        n = Integer(n)
        tbl = gap.UnderlyingCharacterTable(self)
        return ClassFunction(self._group, gap.SymmetricParts(tbl,[self],n)[1])


    def exterior_power(self, n):
        r"""
        Returns the anti-symmetrized product of self with itself ``n`` times.

        INPUT:

        - ``n`` -- a positive integer.

        OUTPUT:

        The ``n``-th anti-symmetrized power of ``self`` as a
        :class:`ClassFunction`.

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: p = chi.exterior_power(3)   # the highest anti-symmetric power for a 3-d character
            sage: p
            Character of Symmetric group of order 4! as a permutation group
            sage: p.values()
            [1, -1, 1, 1, -1]
            sage: p == chi.determinant_character()
            True
        """
        n = Integer(n)
        tbl = gap.UnderlyingCharacterTable(self)
        return ClassFunction(self._group, gap.AntiSymmetricParts(tbl,[self],n)[1])


    def scalar_product(self, other):
        r"""
        Returns the scalar product of self with other.

        EXAMPLES::

            sage: S4 = SymmetricGroup(4)
            sage: irr = S4.irreducible_characters()
            sage: [[x.scalar_product(y) for x in irr] for y in irr]
            [[1, 0, 0, 0, 0],
             [0, 1, 0, 0, 0],
             [0, 0, 1, 0, 0],
             [0, 0, 0, 1, 0],
             [0, 0, 0, 0, 1]]
        """
        return self._gap_classfunction.ScalarProduct(other)


    def is_irreducible(self):
        r"""
        Returns True if self cannot be written as the sum of two nonzero
        characters of self.

        EXAMPLES::

            sage: S4 = SymmetricGroup(4)
            sage: irr = S4.irreducible_characters()
            sage: [x.is_irreducible() for x in irr]
            [True, True, True, True, True]
        """
        return bool(self._gap_classfunction.IsIrreducible())


    def degree(self):
        r"""
        Returns the degree of the character self.

        EXAMPLES::

            sage: S5 = SymmetricGroup(5)
            sage: irr = S5.irreducible_characters()
            sage: [x.degree() for x in irr]
            [1, 4, 5, 6, 5, 4, 1]
        """
        return Integer(self._gap_classfunction.DegreeOfCharacter())


    def irreducible_constituents(self):
        r"""
        Returns a list of the characters that appear in the decomposition
        of chi.

        EXAMPLES::

            sage: S5 = SymmetricGroup(5)
            sage: chi = ClassFunction(S5, [22, -8, 2, 1, 1, 2, -3])
            sage: irr = chi.irreducible_constituents(); irr
            (Character of Symmetric group of order 5! as a permutation group,
             Character of Symmetric group of order 5! as a permutation group)
            sage: list(map(list, irr))
            [[4, -2, 0, 1, 1, 0, -1], [5, -1, 1, -1, -1, 1, 0]]
            sage: G = GL(2,3)
            sage: chi = ClassFunction(G, [-1, -1, -1, -1, -1, -1, -1, -1])
            sage: chi.irreducible_constituents()
            (Character of General Linear Group of degree 2 over Finite Field of size 3,)
            sage: chi = ClassFunction(G, [1, 1, 1, 1, 1, 1, 1, 1])
            sage: chi.irreducible_constituents()
            (Character of General Linear Group of degree 2 over Finite Field of size 3,)
            sage: chi = ClassFunction(G, [2, 2, 2, 2, 2, 2, 2, 2])
            sage: chi.irreducible_constituents()
            (Character of General Linear Group of degree 2 over Finite Field of size 3,)
            sage: chi = ClassFunction(G, [-1, -1, -1, -1, 3, -1, -1, 1])
            sage: ic = chi.irreducible_constituents(); ic
            (Character of General Linear Group of degree 2 over Finite Field of size 3,
             Character of General Linear Group of degree 2 over Finite Field of size 3)
            sage: list(map(list, ic))
            [[2, -1, 2, -1, 2, 0, 0, 0], [3, 0, 3, 0, -1, 1, 1, -1]]
        """
        L = self._gap_classfunction.ConstituentsOfCharacter()
        return tuple(ClassFunction(self._group, list(l)) for l in L)


    def decompose(self):
        r"""
        Returns a list of the characters that appear in the decomposition
        of chi.

        EXAMPLES::

            sage: S5 = SymmetricGroup(5)
            sage: chi = ClassFunction(S5, [22, -8, 2, 1, 1, 2, -3])
            sage: chi.decompose()
            ((3, Character of Symmetric group of order 5! as a permutation group),
             (2, Character of Symmetric group of order 5! as a permutation group))
        """
        L = []
        for irr in self.irreducible_constituents():
            L.append((self.scalar_product(irr), irr))
        return tuple(L)


    def norm(self):
        r"""
        Returns the norm of self.

        EXAMPLES::

            sage: A5 = AlternatingGroup(5)
            sage: [x.norm() for x in A5.irreducible_characters()]
            [1, 1, 1, 1, 1]
        """
        return self._gap_classfunction.Norm()


    def values(self):
        r"""
        Return the list of values of self on the conjugacy classes.

        EXAMPLES::

            sage: G = GL(2,3)
            sage: [x.values() for x in G.irreducible_characters()] #random
            [[1, 1, 1, 1, 1, 1, 1, 1],
             [1, 1, 1, 1, 1, -1, -1, -1],
             [2, -1, 2, -1, 2, 0, 0, 0],
             [2, 1, -2, -1, 0, -zeta8^3 - zeta8, zeta8^3 + zeta8, 0],
             [2, 1, -2, -1, 0, zeta8^3 + zeta8, -zeta8^3 - zeta8, 0],
             [3, 0, 3, 0, -1, -1, -1, 1],
             [3, 0, 3, 0, -1, 1, 1, -1],
             [4, -1, -4, 1, 0, 0, 0, 0]]

        TESTS::

            sage: G = GL(2,3)
            sage: k = CyclotomicField(8)
            sage: zeta8 = k.gen()
            sage: v = [tuple(x.values()) for x in G.irreducible_characters()]
            sage: set(v) == set([(1, 1, 1, 1, 1, 1, 1, 1), (1, 1, 1, 1, 1, -1, -1, -1), (2, -1, 2, -1, 2, 0, 0, 0), (2, 1, -2, -1, 0, -zeta8^3 - zeta8, zeta8^3 + zeta8, 0), (2, 1, -2, -1, 0, zeta8^3 + zeta8, -zeta8^3 - zeta8, 0), (3, 0, 3, 0, -1, -1, -1, 1), (3, 0, 3, 0, -1, 1, 1, -1), (4, -1, -4, 1, 0, 0, 0, 0)])
            True
        """
        return list(self)


    def central_character(self):
        r"""
        Returns the central character of self.

        EXAMPLES::

            sage: t = SymmetricGroup(4).trivial_character()
            sage: t.central_character().values()
            [1, 6, 3, 8, 6]
        """
        return ClassFunction(self._group, self._gap_classfunction.CentralCharacter())


    def determinant_character(self):
        r"""
        Returns the determinant character of self.

        EXAMPLES::

            sage: t = ClassFunction(SymmetricGroup(4), [1, -1, 1, 1, -1])
            sage: t.determinant_character().values()
            [1, -1, 1, 1, -1]
        """
        return ClassFunction(self._group, self._gap_classfunction.DeterminantOfCharacter())


    def tensor_product(self, other):
        r"""
        EXAMPLES::

            sage: S3 = SymmetricGroup(3)
            sage: chi1, chi2, chi3 = S3.irreducible_characters()
            sage: chi1.tensor_product(chi3).values()
            [1, -1, 1]
        """
        return ClassFunction(self._group, gap.Tensored([self],[other])[1])


    def restrict(self, H):
        r"""
        Return the restricted character.

        INPUT:

        - ``H`` -- a subgroup of the underlying group of ``self``.

        OUTPUT:

        A :class:`ClassFunction` of ``H`` defined by restriction.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: chi = ClassFunction(G, [3, -3, -1, 0, 0, -1, 3]); chi
            Character of Symmetric group of order 5! as a permutation group
            sage: H = G.subgroup([(1,2,3), (1,2), (4,5)])
            sage: chi.restrict(H)
            Character of Subgroup generated by [(4,5), (1,2), (1,2,3)] of (Symmetric group of order 5! as a permutation group)
            sage: chi.restrict(H).values()
            [3, -3, -3, -1, 0, 0]
        """
        rest = self._gap_classfunction.RestrictedClassFunction(H._gap_())
        return ClassFunction(H, rest)


    def induct(self, G):
        r"""
        Return the induced character.

        INPUT:

        - ``G`` -- A supergroup of the underlying group of ``self``.

        OUTPUT:

        A :class:`ClassFunction` of ``G`` defined by
        induction. Induction is the adjoint functor to restriction,
        see :meth:`restrict`.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: H = G.subgroup([(1,2,3), (1,2), (4,5)])
            sage: xi = H.trivial_character(); xi
            Character of Subgroup generated by [(4,5), (1,2), (1,2,3)] of (Symmetric group of order 5! as a permutation group)
            sage: xi.induct(G)
            Character of Symmetric group of order 5! as a permutation group
            sage: xi.induct(G).values()
            [10, 4, 2, 1, 1, 0, 0]
        """
        rest = self._gap_classfunction.InducedClassFunction(G._gap_())
        return ClassFunction(G, rest)


    def adams_operation(self, k):
        r"""
        Return the ``k``-th Adams operation on ``self``.

        Let `G` be a finite group. The `k`-th Adams operation `\Psi^k`
        is given by

        .. MATH::

            \Psi^k(\chi)(g) = \chi(g^k).

        The Adams operations turn the representation ring of `G`
        into a `\lambda`-ring.

        EXAMPLES::

            sage: G = groups.permutation.Alternating(5)
            sage: chars = G.irreducible_characters()
            sage: [chi.adams_operation(2).values() for chi in chars]
            [[1, 1, 1, 1, 1],
             [3, 3, 0, -zeta5^3 - zeta5^2, zeta5^3 + zeta5^2 + 1],
             [3, 3, 0, zeta5^3 + zeta5^2 + 1, -zeta5^3 - zeta5^2],
             [4, 4, 1, -1, -1],
             [5, 5, -1, 0, 0]]
            sage: chars[4].adams_operation(2).decompose()
            ((1, Character of Alternating group of order 5!/2 as a permutation group),
             (-1, Character of Alternating group of order 5!/2 as a permutation group),
             (-1, Character of Alternating group of order 5!/2 as a permutation group),
             (2, Character of Alternating group of order 5!/2 as a permutation group))

        REFERENCES:

        - :wikipedia:`Adams_operation`
        """
        reprs = self._group.conjugacy_classes_representatives()
        return ClassFunction(self._group, [self(x**k) for x in reprs])






#####################################################################
###
### Class function using the GAP library
###
#####################################################################


@richcmp_method
class ClassFunction_libgap(SageObject):
    """
    A wrapper of GAP's ``ClassFunction`` function.

    .. NOTE::

        It is *not* checked whether the given values describes a character,
        since GAP does not do this.

    EXAMPLES::

        sage: G = SO(3,3)
        sage: values  = [1, -1, -1, 1, 2]
        sage: chi = ClassFunction(G, values); chi
        Character of Special Orthogonal Group of degree 3 over Finite Field of size 3
        sage: loads(dumps(chi)) == chi
        True
    """

    def __init__(self, G, values):
        r"""
        Return the character of the group ``G`` with values given by the list
        values. The order of the values must correspond to the output of
        ``G.conjugacy_classes_representatives()``.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(4)
            sage: values  = [1, -1, 1, -1]
            sage: chi = ClassFunction(G, values); chi
            Character of Cyclic group of order 4 as a permutation group
        """
        self._group = G
        if isinstance(values, LibGapElement) and values.IsClassFunction():
            self._gap_classfunction = values
        else:
            self._gap_classfunction = libgap.ClassFunction(G._libgap_(),
                                                           list(values))
        e = self._gap_classfunction.Conductor().sage()
        self._base_ring = CyclotomicField(e)


    def gap(self):
        r"""
        Return the underlying LibGAP element.

        EXAMPLES::

            sage: G = CyclicPermutationGroup(4)
            sage: values  = [1, -1, 1, -1]
            sage: chi = ClassFunction(G, values);  chi
            Character of Cyclic group of order 4 as a permutation group
            sage: type(chi)
            <class 'sage.groups.class_function.ClassFunction_gap'>
            sage: gap(chi)
            ClassFunction( CharacterTable( Group( [ (1,2,3,4) ] ) ), [ 1, -1, 1, -1 ] )
            sage: type(_)
            <class 'sage.interfaces.gap.GapElement'>
        """
        return self._gap_classfunction

    _libgap_ = _gap_ = gap


    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: values  = [1, -1, 1, 1, -1]
            sage: ClassFunction(G, values)
            Character of Symmetric group of order 4! as a permutation group
        """
        return "Character of %s" % repr(self._group)


    def __iter__(self):
        r"""
        Iterate through the values.

        A class function assigns values to each conjugacy class. This
        method iterates over the values, in the same order as the
        conjugacy classes of the group.

        EXAMPLES::

            sage: xi = ClassFunction(SymmetricGroup(4), [1, -1, 1, 1, -1])
            sage: list(xi)
            [1, -1, 1, 1, -1]
        """
        for v in self._gap_classfunction.List():
            yield v.sage(ring=self._base_ring)

    def __richcmp__(self, other, op):
        r"""
        Rich comparison for class functions.

        Compares groups and then the values of the class function on the
        conjugacy classes.

        EXAMPLES::

            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: chi = G.character([1, 1, 1, 1, 1, 1, 1])
            sage: H = PermutationGroup([[(1,2,3),(4,5)]])
            sage: xi = H.character([1, 1, 1, 1, 1, 1])
            sage: chi == chi
            True
            sage: xi == xi
            True
            sage: xi == chi
            False
            sage: chi < xi
            False
            sage: xi < chi
            True
        """
        if isinstance(other, ClassFunction_libgap):
            return richcmp((self._group, self.values()),
                           (other._group, other.values()), op)
        else:
            return NotImplemented

    def __reduce__(self):
        r"""
        Add pickle support.

        EXAMPLES::

            sage: G = GL(2,7)
            sage: values = G.gap().CharacterTable().Irr()[2].List().sage()
            sage: chi = ClassFunction(G, values)
            sage: type(chi)
            <class 'sage.groups.class_function.ClassFunction_libgap'>
            sage: loads(dumps(chi)) == chi
            True
        """
        return ClassFunction_libgap, (self._group, self.values())


    def domain(self):
        r"""
        Return the domain of ``self``.

        OUTPUT:

        The underlying group of the class function.

        EXAMPLES::

            sage: ClassFunction(SymmetricGroup(4), [1,-1,1,1,-1]).domain()
            Symmetric group of order 4! as a permutation group
        """
        return self._group


    def __call__(self, g):
        """
        Evaluate the character on the group element `g`.

        Return an error if `g` is not in `G`.

        EXAMPLES::

            sage: G = GL(2,7)
            sage: values = G.gap().CharacterTable().Irr()[2].List().sage()
            sage: chi = ClassFunction(G, values)
            sage: z = G([[3,0],[0,3]]); z
            [3 0]
            [0 3]
            sage: chi(z)
            zeta3

            sage: G = GL(2,3)
            sage: chi = G.irreducible_characters()[3]
            sage: g = G.conjugacy_classes_representatives()[6]
            sage: chi(g)
            zeta8^3 + zeta8

            sage: G = SymmetricGroup(3)
            sage: h = G((2,3))
            sage: triv = G.trivial_character()
            sage: triv(h)
            1
        """
        value = g.gap() ** self.gap()
        return value.sage(self._base_ring)


    def __add__(self, other):
        r"""
        Return the sum of the characters ``self`` and ``other``.

        INPUT:

        - ``other`` -- a :class:`ClassFunction` of the same group as
          ``self``.

        OUTPUT:

        A :class:`ClassFunction`

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: s = chi+chi
            sage: s
            Character of Symmetric group of order 4! as a permutation group
            sage: s.values()
            [6, 2, -2, 0, -2]
        """
        if not isinstance(other, ClassFunction_libgap):
            raise NotImplementedError
        s = self._gap_classfunction + other._gap_classfunction
        return ClassFunction(self._group, s)


    def __sub__(self, other):
        r"""
        Return the difference of the characters ``self`` and ``other``.

        INPUT:

        - ``other`` -- a :class:`ClassFunction` of the same group as
          ``self``.

        OUTPUT:

        A :class:`ClassFunction`

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: chi1 = ClassFunction(G, [3, 1, -1, 0, -1])
            sage: chi2 = ClassFunction(G, [1, -1, 1, 1, -1])
            sage: s = chi1 - chi2
            sage: s
            Character of Symmetric group of order 4! as a permutation group
            sage: s.values()
            [2, 2, -2, -1, 0]
        """
        if not isinstance(other, ClassFunction_libgap):
            raise NotImplementedError
        s = self._gap_classfunction - other._gap_classfunction
        return ClassFunction(self._group, s)


    def __mul__(self, other):
        r"""
        Return the product of the character with ``other``.

        INPUT:

        - ``other`` -- either a number or a :class:`ClassFunction` of
          the same group as ``self``. A number can be anything that
          can be converted into GAP: integers, rational, and elements
          of certain number fields.

        OUTPUT:

        A :class:`ClassFunction`

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: chi1 = ClassFunction(G, [3, 1, -1, 0, -1])
            sage: 3*chi1
            Character of Symmetric group of order 4! as a permutation group
            sage: 3*chi1 == chi1+chi1+chi1
            True
            sage: (3*chi1).values()
            [9, 3, -3, 0, -3]


            sage: (1/2*chi1).values()
            [3/2, 1/2, -1/2, 0, -1/2]

            sage: CF3 = CyclotomicField(3)
            sage: CF3.inject_variables()
            Defining zeta3
            sage: (zeta3 * chi1).values()
            [3*zeta3, zeta3, -zeta3, 0, -zeta3]

            sage: chi2 = ClassFunction(G, [1, -1, 1, 1, -1])
            sage: p = chi1*chi2
            sage: p
            Character of Symmetric group of order 4! as a permutation group
            sage: p.values()
            [3, -1, -1, 0, 1]
        """
        if isinstance(other, ClassFunction_libgap):
            p = self._gap_classfunction * other._gap_classfunction
            return ClassFunction(self._group, p)
        else:
            return ClassFunction(self._group, other * self._gap_classfunction)


    def __rmul__(self, other):
        r"""
        Return the reverse multiplication of ``self`` and ``other``.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: chi = ClassFunction(G, [3, 1, -1, 0, -1])
            sage: chi * 4   # calls chi.__mul__
            Character of Symmetric group of order 4! as a permutation group
            sage: 4 * chi   # calls chi.__rmul__
            Character of Symmetric group of order 4! as a permutation group
            sage: (4 * chi).values()
            [12, 4, -4, 0, -4]
        """
        return self.__mul__(other)


    def __pos__(self):
        r"""
        Return ``self``.

        OUTPUT:

        A :class:`ClassFunction`

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: +chi
            Character of Symmetric group of order 4! as a permutation group
            sage: _.values()
            [3, 1, -1, 0, -1]
            sage: chi.__pos__() == +chi
            True
        """
        return ClassFunction(self._group, self._gap_classfunction)


    def __neg__(self):
        r"""
        Return the additive inverse of ``self``.

        OUTPUT:

        A :class:`ClassFunction`

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: -chi
            Character of Symmetric group of order 4! as a permutation group
            sage: _.values()
            [-3, -1, 1, 0, 1]
            sage: chi.__neg__() == -chi
            True
        """
        return ClassFunction(self._group, -self._gap_classfunction)


    def __pow__(self, other):
        r"""
        Return the product of ``self`` with itself ``other`` times.

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: p = chi**3
            sage: p
            Character of Symmetric group of order 4! as a permutation group
            sage: p.values()
            [27, 1, -1, 0, -1]
        """
        if not isinstance(other, (int,Integer)):
            raise NotImplementedError
        return ClassFunction(self._group, self._gap_classfunction ** other)


    def symmetric_power(self, n):
        r"""
        Return the symmetrized product of ``self`` with itself ``n`` times.

        INPUT:

        - ``n`` -- a positive integer

        OUTPUT:

        The ``n``-th symmetrized power of ``self`` as a
        :class:`ClassFunction`.

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: p = chi.symmetric_power(3)
            sage: p
            Character of Symmetric group of order 4! as a permutation group
            sage: p.values()
            [10, 2, -2, 1, 0]
        """
        n = Integer(n)
        tbl = self._gap_classfunction.UnderlyingCharacterTable(self)
        return ClassFunction(self._group, tbl.SymmetricParts([self],n)[1])


    def exterior_power(self, n):
        r"""
        Return the anti-symmetrized product of ``self`` with itself ``n`` times.

        INPUT:

        - ``n`` -- a positive integer

        OUTPUT:

        The ``n``-th anti-symmetrized power of ``self`` as a
        :class:`ClassFunction`.

        EXAMPLES::

            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: p = chi.exterior_power(3)   # the highest anti-symmetric power for a 3-d character
            sage: p
            Character of Symmetric group of order 4! as a permutation group
            sage: p.values()
            [1, -1, 1, 1, -1]
            sage: p == chi.determinant_character()
            True
        """
        n = Integer(n)
        tbl = self._gap_classfunction.UnderlyingCharacterTable(self)
        return ClassFunction(self._group, tbl.AntiSymmetricParts([self],n)[1])


    def scalar_product(self, other):
        r"""
        Return the scalar product of ``self`` with ``other``.

        EXAMPLES::

            sage: S4 = SymmetricGroup(4)
            sage: irr = S4.irreducible_characters()
            sage: [[x.scalar_product(y) for x in irr] for y in irr]
            [[1, 0, 0, 0, 0],
             [0, 1, 0, 0, 0],
             [0, 0, 1, 0, 0],
             [0, 0, 0, 1, 0],
             [0, 0, 0, 0, 1]]
        """
        return self._gap_classfunction.ScalarProduct(other).sage()


    def is_irreducible(self):
        r"""
        Return ``True`` if ``self`` cannot be written as the sum of two nonzero
        characters of ``self``.

        EXAMPLES::

            sage: S4 = SymmetricGroup(4)
            sage: irr = S4.irreducible_characters()
            sage: [x.is_irreducible() for x in irr]
            [True, True, True, True, True]
        """
        return self._gap_classfunction.IsIrreducible().sage()


    def degree(self):
        r"""
        Return the degree of the character ``self``.

        EXAMPLES::

            sage: S5 = SymmetricGroup(5)
            sage: irr = S5.irreducible_characters()
            sage: [x.degree() for x in irr]
            [1, 4, 5, 6, 5, 4, 1]
        """
        return self._gap_classfunction.DegreeOfCharacter().sage()


    def irreducible_constituents(self):
        r"""
        Return a list of the characters that appear in the decomposition
        of ``self``.

        EXAMPLES::

            sage: S5 = SymmetricGroup(5)
            sage: chi = ClassFunction(S5, [22, -8, 2, 1, 1, 2, -3])
            sage: irr = chi.irreducible_constituents(); irr
            (Character of Symmetric group of order 5! as a permutation group,
             Character of Symmetric group of order 5! as a permutation group)
            sage: list(map(list, irr))
            [[4, -2, 0, 1, 1, 0, -1], [5, -1, 1, -1, -1, 1, 0]]

            sage: G = GL(2,3)
            sage: chi = ClassFunction(G, [-1, -1, -1, -1, -1, -1, -1, -1])
            sage: chi.irreducible_constituents()
            (Character of General Linear Group of degree 2 over Finite Field of size 3,)
            sage: chi = ClassFunction(G, [1, 1, 1, 1, 1, 1, 1, 1])
            sage: chi.irreducible_constituents()
            (Character of General Linear Group of degree 2 over Finite Field of size 3,)
            sage: chi = ClassFunction(G, [2, 2, 2, 2, 2, 2, 2, 2])
            sage: chi.irreducible_constituents()
            (Character of General Linear Group of degree 2 over Finite Field of size 3,)
            sage: chi = ClassFunction(G, [-1, -1, -1, -1, 3, -1, -1, 1])
            sage: ic = chi.irreducible_constituents(); ic
            (Character of General Linear Group of degree 2 over Finite Field of size 3,
             Character of General Linear Group of degree 2 over Finite Field of size 3)
            sage: list(map(list, ic))
            [[2, -1, 2, -1, 2, 0, 0, 0], [3, 0, 3, 0, -1, 1, 1, -1]]
        """
        L = self._gap_classfunction.ConstituentsOfCharacter()
        return tuple(ClassFunction_libgap(self._group, l) for l in L)


    def decompose(self):
        r"""
        Return a list of the characters that appear in the decomposition
        of ``self``.

        EXAMPLES::

            sage: S5 = SymmetricGroup(5)
            sage: chi = ClassFunction(S5, [22, -8, 2, 1, 1, 2, -3])
            sage: chi.decompose()
            ((3, Character of Symmetric group of order 5! as a permutation group),
             (2, Character of Symmetric group of order 5! as a permutation group))
        """
        L = []
        for irr in self.irreducible_constituents():
            L.append((self.scalar_product(irr), irr))
        return tuple(L)


    def norm(self):
        r"""
        Return the norm of ``self``.

        EXAMPLES::

            sage: A5 = AlternatingGroup(5)
            sage: [x.norm() for x in A5.irreducible_characters()]
            [1, 1, 1, 1, 1]
        """
        return self._gap_classfunction.Norm().sage()


    def values(self):
        r"""
        Return the list of values of self on the conjugacy classes.

        EXAMPLES::

            sage: G = GL(2,3)
            sage: [x.values() for x in G.irreducible_characters()] #random
            [[1, 1, 1, 1, 1, 1, 1, 1],
             [1, 1, 1, 1, 1, -1, -1, -1],
             [2, -1, 2, -1, 2, 0, 0, 0],
             [2, 1, -2, -1, 0, -zeta8^3 - zeta8, zeta8^3 + zeta8, 0],
             [2, 1, -2, -1, 0, zeta8^3 + zeta8, -zeta8^3 - zeta8, 0],
             [3, 0, 3, 0, -1, -1, -1, 1],
             [3, 0, 3, 0, -1, 1, 1, -1],
             [4, -1, -4, 1, 0, 0, 0, 0]]

        TESTS::

            sage: G = GL(2,3)
            sage: k = CyclotomicField(8)
            sage: zeta8 = k.gen()
            sage: v = [tuple(x.values()) for x in G.irreducible_characters()]
            sage: set(v) == set([(1, 1, 1, 1, 1, 1, 1, 1), (1, 1, 1, 1, 1, -1, -1, -1), (2, -1, 2, -1, 2, 0, 0, 0), (2, 1, -2, -1, 0, -zeta8^3 - zeta8, zeta8^3 + zeta8, 0), (2, 1, -2, -1, 0, zeta8^3 + zeta8, -zeta8^3 - zeta8, 0), (3, 0, 3, 0, -1, -1, -1, 1), (3, 0, 3, 0, -1, 1, 1, -1), (4, -1, -4, 1, 0, 0, 0, 0)])
            True
        """
        return list(self)


    def central_character(self):
        r"""
        Return the central character of ``self``.

        EXAMPLES::

            sage: t = SymmetricGroup(4).trivial_character()
            sage: t.central_character().values()
            [1, 6, 3, 8, 6]
        """
        return ClassFunction(self._group, self._gap_classfunction.CentralCharacter())


    def determinant_character(self):
        r"""
        Return the determinant character of ``self``.

        EXAMPLES::

            sage: t = ClassFunction(SymmetricGroup(4), [1, -1, 1, 1, -1])
            sage: t.determinant_character().values()
            [1, -1, 1, 1, -1]
        """
        return ClassFunction(self._group, self._gap_classfunction.DeterminantOfCharacter())


    def tensor_product(self, other):
        r"""
        Return the tensor product of ``self`` and ``other``.

        EXAMPLES::

            sage: S3 = SymmetricGroup(3)
            sage: chi1, chi2, chi3 = S3.irreducible_characters()
            sage: chi1.tensor_product(chi3).values()
            [1, -1, 1]
        """
        product = libgap.Tensored([self], [other])
        return ClassFunction(self._group, product[0])


    def restrict(self, H):
        r"""
        Return the restricted character.

        INPUT:

        - ``H`` -- a subgroup of the underlying group of ``self``.

        OUTPUT:

        A :class:`ClassFunction` of ``H`` defined by restriction.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: chi = ClassFunction(G, [3, -3, -1, 0, 0, -1, 3]); chi
            Character of Symmetric group of order 5! as a permutation group
            sage: H = G.subgroup([(1,2,3), (1,2), (4,5)])
            sage: chi.restrict(H)
            Character of Subgroup generated by [(4,5), (1,2), (1,2,3)] of (Symmetric group of order 5! as a permutation group)
            sage: chi.restrict(H).values()
            [3, -3, -3, -1, 0, 0]
        """
        try:
            gapH = H.gap()
        except AttributeError:
            from sage.libs.gap.libgap import libgap
            gapH = libgap(H)
        rest = self._gap_classfunction.RestrictedClassFunction(gapH)
        return ClassFunction(H, rest)


    def induct(self, G):
        r"""
        Return the induced character.

        INPUT:

        - ``G`` -- A supergroup of the underlying group of ``self``.

        OUTPUT:

        A :class:`ClassFunction` of ``G`` defined by
        induction. Induction is the adjoint functor to restriction,
        see :meth:`restrict`.

        EXAMPLES::

            sage: G = SymmetricGroup(5)
            sage: H = G.subgroup([(1,2,3), (1,2), (4,5)])
            sage: xi = H.trivial_character(); xi
            Character of Subgroup generated by [(4,5), (1,2), (1,2,3)] of (Symmetric group of order 5! as a permutation group)
            sage: xi.induct(G)
            Character of Symmetric group of order 5! as a permutation group
            sage: xi.induct(G).values()
            [10, 4, 2, 1, 1, 0, 0]
        """
        try:
            gapG = G.gap()
        except AttributeError:
            from sage.libs.gap.libgap import libgap
            gapG = libgap(G)
        ind = self._gap_classfunction.InducedClassFunction(gapG)
        return ClassFunction(G, ind)


    def adams_operation(self, k):
        r"""
        Return the ``k``-th Adams operation on ``self``.

        Let `G` be a finite group. The `k`-th Adams operation `\Psi^k`
        is given by

        .. MATH::

            \Psi^k(\chi)(g) = \chi(g^k).

        The Adams operations turn the representation ring of `G`
        into a `\lambda`-ring.

        EXAMPLES::

            sage: G = GL(2,3)
            sage: chars = G.irreducible_characters()
            sage: [chi.adams_operation(2).values() for chi in chars]
            [[1, 1, 1, 1, 1, 1, 1, 1],
             [1, 1, 1, 1, 1, 1, 1, 1],
             [2, -1, 2, -1, 2, 2, 2, 2],
             [2, -1, 2, -1, -2, 0, 0, 2],
             [2, -1, 2, -1, -2, 0, 0, 2],
             [3, 0, 3, 0, 3, -1, -1, 3],
             [3, 0, 3, 0, 3, -1, -1, 3],
             [4, 1, 4, 1, -4, 0, 0, 4]]
            sage: chars[5].adams_operation(3).decompose()
            ((1, Character of General Linear Group of degree 2 over Finite Field of size 3),
             (1, Character of General Linear Group of degree 2 over Finite Field of size 3),
             (-1, Character of General Linear Group of degree 2 over Finite Field of size 3),
             (1, Character of General Linear Group of degree 2 over Finite Field of size 3))

        REFERENCES:

        - :wikipedia:`Adams_operation`
        """
        reprs = self._group.conjugacy_classes_representatives()
        return ClassFunction(self._group, [self(x**k) for x in reprs])

