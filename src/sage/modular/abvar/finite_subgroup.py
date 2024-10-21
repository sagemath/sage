# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.pari
r"""
Finite subgroups of modular abelian varieties

Sage can compute with fairly general finite subgroups of modular
abelian varieties. Elements of finite order are represented by
equivalence classes of elements in `H_1(A,\QQ)`
modulo `H_1(A,\ZZ)`. A finite subgroup can be
defined by giving generators and via various other constructions.
Given a finite subgroup, one can compute generators, as well as the
structure as an abstract group. Arithmetic on subgroups is also
supported, including adding two subgroups together, checking
inclusion, etc.

TODO: Intersection, action of Hecke operators.

AUTHORS:

- William Stein (2007-03)

EXAMPLES::

    sage: J = J0(33)
    sage: C = J.cuspidal_subgroup()
    sage: C
    Finite subgroup with invariants [10, 10] over QQ of Abelian variety J0(33) of dimension 3
    sage: C.order()
    100
    sage: C.gens()
    [[(1/10, 0, 1/10, 1/10, 1/10, 3/10)], [(0, 1/5, 1/10, 0, 1/10, 9/10)], [(0, 0, 1/2, 0, 1/2, 1/2)]]
    sage: C.0 + C.1
    [(1/10, 1/5, 1/5, 1/10, 1/5, 6/5)]
    sage: 10*(C.0 + C.1)
    [(0, 0, 0, 0, 0, 0)]
    sage: G = C.subgroup([C.0 + C.1]); G
    Finite subgroup with invariants [10] over QQbar of Abelian variety J0(33) of dimension 3
    sage: G.gens()
    [[(1/10, 1/5, 1/5, 1/10, 1/5, 1/5)]]
    sage: G.order()
    10
    sage: G <= C
    True
    sage: G >= C
    False

We make a table of the order of the cuspidal subgroup for the first
few levels::

    sage: for N in range(11,40):
    ....:     print("{} {}".format(N, J0(N).cuspidal_subgroup().order()))
    11 5
    12 1
    13 1
    14 6
    15 8
    16 1
    17 4
    18 1
    19 3
    20 6
    21 8
    22 25
    23 11
    24 8
    25 1
    26 21
    27 9
    28 36
    29 7
    30 192
    31 5
    32 8
    33 100
    34 48
    35 48
    36 12
    37 3
    38 135
    39 56

TESTS::

    sage: G = J0(11).finite_subgroup([[1/3,0], [0,1/5]]); G
    Finite subgroup with invariants [15] over QQbar of Abelian variety J0(11) of dimension 1
    sage: loads(dumps(G)) == G
    True
    sage: loads(dumps(G.0)) == G.0
    True
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.rings.abc

from sage.misc.lazy_import import lazy_import
from sage.modular.abvar.torsion_point import TorsionPoint
from sage.modules.module import Module
from sage.modules.free_module import FreeModule_generic
from sage.structure.gens_py import abelian_iterator
from sage.structure.sequence import Sequence
from sage.structure.richcmp import richcmp_method, richcmp
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.functions import lcm
from sage.misc.misc_c import prod
from sage.structure.element import coercion_model


@richcmp_method
class FiniteSubgroup(Module):
    r"""
    A finite subgroup of a modular abelian variety.

    INPUT:

    - ``abvar`` -- a modular abelian variety

    - ``field_of_definition`` -- a field over which this group is defined

    EXAMPLES:

    This is an abstract base class, so there are no instances of
    this class itself::

        sage: A = J0(37)
        sage: G = A.torsion_subgroup(3); G
        Finite subgroup with invariants [3, 3, 3, 3] over QQ of Abelian variety J0(37) of dimension 2
        sage: type(G)
        <class 'sage.modular.abvar.finite_subgroup.FiniteSubgroup_lattice_with_category'>
        sage: from sage.modular.abvar.finite_subgroup import FiniteSubgroup
        sage: isinstance(G, FiniteSubgroup)
        True
    """

    Element = TorsionPoint

    def __init__(self, abvar, field_of_definition=QQ):
        """
        Initialize ``self``.

        TESTS::

            sage: A = J0(11)
            sage: G = A.torsion_subgroup(2)
            sage: TestSuite(G).run() # long time
        """
        from sage.categories.category import Category
        from sage.categories.fields import Fields
        from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
        from sage.categories.modules import Modules
        from .abvar import ModularAbelianVariety_abstract
        if field_of_definition not in Fields():
            raise TypeError("field_of_definition must be a field")
        if not isinstance(abvar, ModularAbelianVariety_abstract):
            raise TypeError("abvar must be a modular abelian variety")
        category = Category.join((Modules(ZZ), FiniteEnumeratedSets()))
        Module.__init__(self, ZZ, category=category)
        self.__abvar = abvar
        self.__field_of_definition = field_of_definition

    ################################################################
    # DERIVED CLASS MUST OVERRIDE THE lattice METHOD
    ################################################################
    def lattice(self):
        """
        Return the lattice corresponding to this subgroup in the rational
        homology of the modular Jacobian product. The elements of the
        subgroup are represented by vectors in the ambient vector space
        (the rational homology), and this returns the lattice they span.
        EXAMPLES::

            sage: J = J0(33); C = J[0].cuspidal_subgroup(); C
            Finite subgroup with invariants [5] over QQ of
             Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: C.lattice()
            Free module of degree 6 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1/5 13/5   -2 -4/5    2 -1/5]
            [   0    3   -2   -1    2    0]
        """
        raise NotImplementedError

    def _relative_basis_matrix(self):
        """
        Return matrix of this finite subgroup, but relative to the homology
        of the parent abelian variety.

        EXAMPLES::

            sage: A = J0(43)[1]; A
            Simple abelian subvariety 43b(1,43) of dimension 2 of J0(43)
            sage: C = A.cuspidal_subgroup(); C
            Finite subgroup with invariants [7] over QQ of
             Simple abelian subvariety 43b(1,43) of dimension 2 of J0(43)
            sage: C._relative_basis_matrix()
            [  1   0   0   0]
            [  0 1/7 6/7 5/7]
            [  0   0   1   0]
            [  0   0   0   1]
        """
        try:
            return self.__relative_basis_matrix
        except AttributeError:
            M = self.__abvar.lattice().coordinate_module(self.lattice()).basis_matrix()
            self.__relative_basis_matrix = M
            return M

    # General functionality
    def __richcmp__(self, other, op):
        """
        Compare ``self`` to ``other``.

        If ``other`` is not a :class:`FiniteSubgroup`, then
        ``NotImplemented`` is returned. If ``other`` is a
        :class:`FiniteSubgroup` and the ambient abelian varieties are
        not equal, then the ambient abelian varieties are compared.
        If ``other`` is a :class:`FiniteSubgroup` and the ambient
        abelian varieties are equal, then the subgroups are compared
        via their corresponding lattices.

        EXAMPLES:

        We first compare two subgroups of `J_0(37)`::

            sage: A = J0(37)
            sage: G = A.torsion_subgroup(3); G.order()
            81
            sage: H = A.cuspidal_subgroup(); H.order()
            3
            sage: H < G
            True
            sage: H.is_subgroup(G)
            True

        The ambient varieties are compared::

            sage: A[0].cuspidal_subgroup() > J0(11).cuspidal_subgroup()
            True

        Comparing subgroups sitting in different abelian varieties::

            sage: A[0].cuspidal_subgroup() < A[1].cuspidal_subgroup()
            True
        """
        if not isinstance(other, FiniteSubgroup):
            return NotImplemented
        A = self.abelian_variety()
        B = other.abelian_variety()
        if not A.in_same_ambient_variety(B):
            return richcmp(A.ambient_variety(), B.ambient_variety(), op)
        L = A.lattice() + B.lattice()
        lx = other.lattice() + L
        rx = self.lattice() + L
        # order gets reversed in passing to lattices.
        return lx._echelon_matrix_richcmp(rx, op)

    def is_subgroup(self, other):
        """
        Return ``True`` exactly if ``self`` is a subgroup of ``other``, and both are
        defined as subgroups of the same ambient abelian variety.

        EXAMPLES::

            sage: C = J0(22).cuspidal_subgroup()
            sage: H = C.subgroup([C.0])
            sage: K = C.subgroup([C.1])
            sage: H.is_subgroup(K)
            False
            sage: K.is_subgroup(H)
            False
            sage: K.is_subgroup(C)
            True
            sage: H.is_subgroup(C)
            True
        """
        # We use that self is contained in other, whether other is
        # either a finite group or an abelian variety, if and only
        # if self doesn't shrink when intersected with other.
        try:
            return self.intersection(other).order() == self.order()
        except TypeError:
            return False

    def __add__(self, other):
        """
        Return the sum of two subgroups.

        EXAMPLES::

            sage: C = J0(22).cuspidal_subgroup()
            sage: C.gens()
            [[(1/5, 1/5, 4/5, 0)], [(0, 0, 0, 1/5)]]
            sage: A = C.subgroup([C.0]); B = C.subgroup([C.1])
            sage: A + B == C
            True

        An example where the parent abelian varieties are different::

            sage: A = J0(48); A[0].cuspidal_subgroup() + A[1].cuspidal_subgroup()
            Finite subgroup with invariants [2, 4, 4] over QQ of
             Abelian subvariety of dimension 2 of J0(48)
        """
        if not isinstance(other, FiniteSubgroup):
            raise TypeError("only addition of two finite subgroups is defined")
        A = self.abelian_variety()
        B = other.abelian_variety()
        if not A.in_same_ambient_variety(B):
            raise ValueError("self and other must be in the same ambient Jacobian")
        K = coercion_model.common_parent(self.field_of_definition(), other.field_of_definition())
        lattice = self.lattice() + other.lattice()
        if A != B:
            C = A + B
            lattice += C.lattice()
            return FiniteSubgroup_lattice(C, lattice, field_of_definition=K)
        else:
            return FiniteSubgroup_lattice(self.abelian_variety(), lattice, field_of_definition=K)

    def exponent(self):
        """
        Return the exponent of this finite abelian group.

        OUTPUT: integer

        EXAMPLES::

            sage: t = J0(33).hecke_operator(7)
            sage: G = t.kernel()[0]; G
            Finite subgroup with invariants [2, 2, 2, 2, 4, 4] over QQ of
             Abelian variety J0(33) of dimension 3
            sage: G.exponent()
            4
        """
        try:
            return self.__exponent
        except AttributeError:
            e = lcm(self.invariants())
            self.__exponent = e
            return e

    def intersection(self, other):
        """
        Return the intersection of the finite subgroups ``self`` and ``other``.

        INPUT:

        - ``other`` -- a finite group

        OUTPUT: a finite group

        EXAMPLES::

            sage: E11a0, E11a1, B = J0(33)
            sage: G = E11a0.torsion_subgroup(6); H = E11a0.torsion_subgroup(9)
            sage: G.intersection(H)
            Finite subgroup with invariants [3, 3] over QQ of
             Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: W = E11a1.torsion_subgroup(15)
            sage: G.intersection(W)
            Finite subgroup with invariants [] over QQ of
             Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: E11a0.intersection(E11a1)[0]
            Finite subgroup with invariants [5] over QQ of
             Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)

        We intersect subgroups of different abelian varieties.

        ::

            sage: E11a0, E11a1, B = J0(33)
            sage: G = E11a0.torsion_subgroup(5); H = E11a1.torsion_subgroup(5)
            sage: G.intersection(H)
            Finite subgroup with invariants [5] over QQ of
             Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: E11a0.intersection(E11a1)[0]
            Finite subgroup with invariants [5] over QQ of
             Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)

        We intersect abelian varieties with subgroups::

            sage: t = J0(33).hecke_operator(7)
            sage: G = t.kernel()[0]; G
            Finite subgroup with invariants [2, 2, 2, 2, 4, 4] over QQ of
             Abelian variety J0(33) of dimension 3
            sage: A = J0(33).old_subvariety()
            sage: A.intersection(G)
            Finite subgroup with invariants [2, 2, 2, 2] over QQ of
             Abelian subvariety of dimension 2 of J0(33)
            sage: A.hecke_operator(7).kernel()[0]
            Finite subgroup with invariants [2, 2, 2, 2] over QQ of
             Abelian subvariety of dimension 2 of J0(33)
            sage: B = J0(33).new_subvariety()
            sage: B.intersection(G)
            Finite subgroup with invariants [4, 4] over QQ of
             Abelian subvariety of dimension 1 of J0(33)
            sage: B.hecke_operator(7).kernel()[0]
            Finite subgroup with invariants [4, 4] over QQ of
             Abelian subvariety of dimension 1 of J0(33)
            sage: A.intersection(B)[0]
            Finite subgroup with invariants [3, 3] over QQ of
             Abelian subvariety of dimension 2 of J0(33)
        """
        from .abvar import ModularAbelianVariety_abstract
        A = self.abelian_variety()
        if isinstance(other, ModularAbelianVariety_abstract):
            amb = other
            B = other
            M = B.lattice().scale(Integer(1)/self.exponent())
            K = coercion_model.common_parent(self.field_of_definition(), other.base_field())
        else:
            amb = A
            if not isinstance(other, FiniteSubgroup):
                raise TypeError("only intersection with a finite subgroup or "
                        "modular abelian variety is defined")
            B = other.abelian_variety()
            if A.ambient_variety() != B.ambient_variety():
                raise TypeError("finite subgroups must be in the same ambient product Jacobian")
            M = other.lattice()
            K = coercion_model.common_parent(self.field_of_definition(), other.field_of_definition())

        L = self.lattice()
        if A != B:
            # TODO: This might be way slower than what we could do if
            # we think more carefully.
            C = A + B
            L = L + C.lattice()
            M = M + C.lattice()
        W = L.intersection(M).intersection(amb.vector_space())
        return FiniteSubgroup_lattice(amb, W, field_of_definition=K)

    def __mul__(self, right):
        """
        Multiply this subgroup by the rational number ``right``.

        If ``right`` is an integer the result is a subgroup of ``self``. If
        ``right`` is a rational number `n/m`, then this group is first
        divided by `m` then multiplied by `n`.

        INPUT:

        - ``right`` -- a rational number

        OUTPUT: a subgroup

        EXAMPLES::

            sage: J = J0(37)
            sage: H = J.cuspidal_subgroup(); H.order()
            3
            sage: G = H * 3; G.order()
            1
            sage: G = H * (1/2); G.order()
            48
            sage: J.torsion_subgroup(2) + H == G
            True
            sage: G = H*(3/2); G.order()
            16
            sage: J = J0(42)
            sage: G = J.cuspidal_subgroup(); factor(G.order())
            2^8 * 3^2
            sage: (G * 3).order()
            256
            sage: (G * 0).order()
            1
            sage: (G * (1/5)).order()
            22500000000
        """
        lattice = self.lattice().scale(right)
        return FiniteSubgroup_lattice(self.abelian_variety(), lattice,
                                field_of_definition=self.field_of_definition())

    def __rmul__(self, left):
        """
        Multiply this finite subgroup on the left by an integer.

        EXAMPLES::

            sage: J = J0(42)
            sage: G = J.cuspidal_subgroup(); factor(G.order())
            2^8 * 3^2
            sage: H = G.__rmul__(2)
            sage: H.order().factor()
            2^4 * 3^2
            sage: 2*G
            Finite subgroup with invariants [6, 24] over QQ of Abelian variety J0(42) of dimension 5
        """
        return self * left

    def abelian_variety(self):
        """
        Return the abelian variety that this is a finite subgroup of.

        EXAMPLES::

            sage: J = J0(42)
            sage: G = J.rational_torsion_subgroup(); G
            Torsion subgroup of Abelian variety J0(42) of dimension 5
            sage: G.abelian_variety()
            Abelian variety J0(42) of dimension 5
        """
        return self.__abvar

    def field_of_definition(self):
        """
        Return the field over which this finite modular abelian variety
        subgroup is defined. This is a field over which this subgroup is
        defined.

        EXAMPLES::

            sage: J = J0(42)
            sage: G = J.rational_torsion_subgroup(); G
            Torsion subgroup of Abelian variety J0(42) of dimension 5
            sage: G.field_of_definition()
            Rational Field
        """
        return self.__field_of_definition

    def _repr_(self):
        """
        Return string representation of this finite subgroup.

        EXAMPLES::

            sage: J = J0(42)
            sage: G = J.torsion_subgroup(3); G._repr_()
            'Finite subgroup with invariants [3, 3, 3, 3, 3, 3, 3, 3, 3, 3] over QQ of Abelian variety J0(42) of dimension 5'
        """
        K = self.__field_of_definition
        if isinstance(K, sage.rings.abc.AlgebraicField):
            field = "QQbar"
        elif K == QQ:
            field = "QQ"
        else:
            field = str(K)
        return "Finite subgroup %sover %s of %s" % (self._invariants_repr(), field, self.__abvar)

    def _invariants_repr(self):
        """
        The string representation of the 'invariants' part of this group.

        We make this a separate function so it is possible to create finite
        subgroups that don't print their invariants, since printing them
        could be expensive.

        EXAMPLES::

            sage: J0(42).cuspidal_subgroup()._invariants_repr()
            'with invariants [2, 2, 12, 48] '
        """
        return 'with invariants %s ' % (self.invariants(), )

    def order(self):
        """
        Return the order (number of elements) of this finite subgroup.

        EXAMPLES::

            sage: J = J0(42)
            sage: C = J.cuspidal_subgroup()
            sage: C.order()
            2304
        """
        try:
            return self.__order
        except AttributeError:
            if self.__abvar.dimension() == 0:
                self.__order = ZZ(1)
                return self.__order
            o = prod(self.invariants())
            self.__order = o
            return o

    def gens(self):
        """
        Return generators for this finite subgroup.

        EXAMPLES: We list generators for several cuspidal subgroups::

            sage: J0(11).cuspidal_subgroup().gens()
            [[(0, 1/5)]]
            sage: J0(37).cuspidal_subgroup().gens()
            [[(0, 0, 0, 1/3)]]
            sage: J0(43).cuspidal_subgroup().gens()
            [[(0, 1/7, 0, 6/7, 0, 5/7)]]
            sage: J1(13).cuspidal_subgroup().gens()
            [[(1/19, 0, 9/19, 9/19)], [(0, 1/19, 0, 9/19)]]
            sage: J0(22).torsion_subgroup(6).gens()
            [[(1/6, 0, 0, 0)], [(0, 1/6, 0, 0)], [(0, 0, 1/6, 0)], [(0, 0, 0, 1/6)]]
        """
        try:
            return self.__gens
        except AttributeError:
            pass

        B = [self.element_class(self, v) for v in self.lattice().basis() if v.denominator() > 1]
        self.__gens = Sequence(B, immutable=True)
        return self.__gens

    def gen(self, n):
        r"""
        Return `n`-th generator of ``self``.

        EXAMPLES::

            sage: J = J0(23)
            sage: C = J.torsion_subgroup(3)
            sage: C.gens()
            [[(1/3, 0, 0, 0)], [(0, 1/3, 0, 0)], [(0, 0, 1/3, 0)], [(0, 0, 0, 1/3)]]
            sage: C.gen(0)
            [(1/3, 0, 0, 0)]
            sage: C.gen(3)
            [(0, 0, 0, 1/3)]
            sage: C.gen(4)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range

        Negative indices wrap around::

            sage: C.gen(-1)
            [(0, 0, 0, 1/3)]
        """
        return self.gens()[n]

    def _element_constructor_(self, x, check=True):
        r"""
        Convert `x` into this finite subgroup.

        This works when the abelian varieties that contain `x` and
        ``self`` are the same, or if `x` is convertible into the
        rational homology (viewed as an abstract `\QQ`-vector space).

        EXAMPLES: We first construct the `11`-torsion subgroup of
        `J_0(23)`::

            sage: J = J0(23)
            sage: G = J.torsion_subgroup(11)
            sage: G.invariants()
            [11, 11, 11, 11]

        We also construct the cuspidal subgroup::

            sage: C = J.cuspidal_subgroup()
            sage: C.invariants()
            [11]

            sage: G(G.0) is G.0
            True

        We convert an element from the cuspidal subgroup into the
        `11`-torsion subgroup::

            sage: z = G(C.0); z
            [(1/11, 10/11, 0, 8/11)]
            sage: z.parent() == G
            True

        We convert a list, which defines an element of the underlying
        ``full_module`` into `G`, and verify an equality::

            sage: x = G([1/11, 1/11, 0, -1/11])
            sage: x == G([1/11, 1/11, 0, 10/11])
            True

        Finally we attempt to convert some elements that shouldn't
        work, since they are not in `G`::

            sage: G(J.torsion_subgroup(3).0)
            Traceback (most recent call last):
            ...
            TypeError: element [1/3, 0, 0, 0] is not in free module

            sage: G(J0(27).cuspidal_subgroup()(0))
            Traceback (most recent call last):
            ...
            ValueError: ambient abelian varieties are different
        """
        if isinstance(x, TorsionPoint):
            if x.parent().abelian_variety() != self.abelian_variety():
                raise ValueError('ambient abelian varieties are different')
            x = x.element()
        x = self.lattice()(x, check=check)
        return self.element_class(self, x, check=False)

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is contained in this finite subgroup.

        EXAMPLES:

        We define two distinct finite subgroups of `J_0(27)`::

            sage: G1 = J0(27).rational_cusp_subgroup(); G1
            Finite subgroup with invariants [3] over QQ of Abelian variety J0(27) of dimension 1
            sage: G1.0
            [(1/3, 0)]
            sage: G2 = J0(27).cuspidal_subgroup(); G2
            Finite subgroup with invariants [3, 3] over QQ of Abelian variety J0(27) of dimension 1
            sage: G2.gens()
            [[(1/3, 0)], [(0, 1/3)]]

        Now we check whether various elements are in `G_1` and `G_2`::

            sage: G2.0 in G1
            True
            sage: G2.1 in G1
            False
            sage: G1.0 in G1
            True
            sage: G1.0 in G2
            True

        The integer `0` is in `G_1`::

            sage: 0 in G1
            True

        Elements that have a completely different ambient product Jacobian
        are never in `G`::

            sage: J0(23).cuspidal_subgroup().0 in G1
            False
            sage: J0(23).cuspidal_subgroup()(0) in G1
            False
        """
        try:
            self(x)
        except (TypeError, ValueError):
            return False
        return True

    def subgroup(self, gens):
        """
        Return the subgroup of ``self`` spanned by the given
        generators, which must all be elements of ``self``.

        EXAMPLES::

            sage: J = J0(23)
            sage: G = J.torsion_subgroup(11); G
            Finite subgroup with invariants [11, 11, 11, 11] over QQ of
             Abelian variety J0(23) of dimension 2

        We create the subgroup of the 11-torsion subgroup of `J_0(23)`
        generated by the first `11`-torsion point::

            sage: H = G.subgroup([G.0]); H
            Finite subgroup with invariants [11] over QQbar of
             Abelian variety J0(23) of dimension 2
            sage: H.invariants()
            [11]

        We can also create a subgroup from a list of objects that can
        be converted into the ambient rational homology::

            sage: H == G.subgroup([[1/11,0,0,0]])
            True
        """
        from sage.rings.qqbar import QQbar

        if not isinstance(gens, (tuple, list)):
            raise TypeError("gens must be a list or tuple")
        A = self.abelian_variety()
        lattice = A._ambient_lattice().span([self(g).element() for g in gens])
        return FiniteSubgroup_lattice(self.abelian_variety(), lattice, field_of_definition=QQbar)

    def invariants(self):
        r"""
        Return elementary invariants of this abelian group, by which we
        mean a nondecreasing (immutable) sequence of integers
        `n_i`, `1 \leq i \leq k`, with `n_i`
        dividing `n_{i+1}`, and such that this group is abstractly
        isomorphic to
        `\ZZ/n_1\ZZ \times\cdots\times \ZZ/n_k\ZZ.`

        EXAMPLES::

            sage: J = J0(38)
            sage: C = J.cuspidal_subgroup(); C
            Finite subgroup with invariants [3, 45] over QQ of
             Abelian variety J0(38) of dimension 4
            sage: v = C.invariants(); v
            [3, 45]
            sage: v[0] = 5
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: type(v[0])
            <class 'sage.rings.integer.Integer'>

        ::

            sage: C * 3
            Finite subgroup with invariants [15] over QQ of
             Abelian variety J0(38) of dimension 4

        An example involving another cuspidal subgroup::

            sage: C = J0(22).cuspidal_subgroup(); C
            Finite subgroup with invariants [5, 5] over QQ of
             Abelian variety J0(22) of dimension 2
            sage: C.lattice()
            Free module of degree 4 and rank 4 over Integer Ring
            Echelon basis matrix:
            [1/5 1/5 4/5   0]
            [  0   1   0   0]
            [  0   0   1   0]
            [  0   0   0 1/5]
            sage: C.invariants()
            [5, 5]
        """
        try:
            return self.__invariants
        except AttributeError:
            pass
        M = self.lattice().coordinate_module(self.abelian_variety().lattice())
        E = M.basis_matrix().change_ring(ZZ).elementary_divisors()
        v = [Integer(x) for x in E if x != 1]
        I = Sequence(v)
        I.sort()
        I.set_immutable()
        self.__invariants = I
        return I

    __iter__ = abelian_iterator


class FiniteSubgroup_lattice(FiniteSubgroup):
    def __init__(self, abvar, lattice, field_of_definition=None, check=True):
        """
        A finite subgroup of a modular abelian variety that is defined by a
        given lattice.

        INPUT:

        - ``abvar`` -- a modular abelian variety

        - ``lattice`` -- a lattice that contains the lattice of abvar

        - ``field_of_definition`` -- the field of definition
          of this finite group scheme

        - ``check`` -- boolean (default: ``True``); whether or not to
          check that lattice contains the abvar lattice

        EXAMPLES::

            sage: J = J0(11)
            sage: G = J.finite_subgroup([[1/3,0], [0,1/5]]); G
            Finite subgroup with invariants [15] over QQbar of
             Abelian variety J0(11) of dimension 1
        """
        if field_of_definition is None:
            from sage.rings.qqbar import QQbar as field_of_definition
        if check:
            from .abvar import ModularAbelianVariety_abstract
            if not isinstance(lattice, FreeModule_generic) or lattice.base_ring() != ZZ:
                raise TypeError("lattice must be a free module over ZZ")
            if not isinstance(abvar, ModularAbelianVariety_abstract):
                raise TypeError("abvar must be a modular abelian variety")
            if not abvar.lattice().is_submodule(lattice):
                lattice += abvar.lattice()
            if lattice.rank() != abvar.lattice().rank():
                raise ValueError("lattice must contain the lattice of abvar with finite index")
        FiniteSubgroup.__init__(self, abvar, field_of_definition)
        self.__lattice = lattice

    def lattice(self):
        r"""
        Return lattice that defines this finite subgroup.

        EXAMPLES::

            sage: J = J0(11)
            sage: G = J.finite_subgroup([[1/3,0], [0,1/5]]); G
            Finite subgroup with invariants [15] over QQbar of
             Abelian variety J0(11) of dimension 1
            sage: G.lattice()
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1/3   0]
            [  0 1/5]
        """
        return self.__lattice
