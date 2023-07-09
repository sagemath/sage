# -*- coding: utf-8 -*-
"""
Base class for modular abelian varieties

AUTHORS:

- William Stein (2007-03)

TESTS::

    sage: A = J0(33)
    sage: D = A.decomposition(); D
    [
    Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33),
    Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33),
    Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
    ]
    sage: loads(dumps(D)) == D
    True
    sage: loads(dumps(A)) == A
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

from copy import copy
from random import randrange

import sage.rings.abc

from sage.arith.functions import lcm as LCM
from sage.arith.misc import divisors, next_prime, is_prime
from sage.categories.modular_abelian_varieties import ModularAbelianVarieties
from sage.matrix.constructor import matrix
from sage.matrix.special import block_diagonal_matrix, identity_matrix
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.modular.arithgroup.congroup_gamma0 import is_Gamma0
from sage.modular.arithgroup.congroup_gamma1 import is_Gamma1
from sage.modular.arithgroup.congroup_gammaH import is_GammaH
from sage.modular.arithgroup.congroup_generic import is_CongruenceSubgroup
from sage.modular.modform.constructor import Newform
from sage.modular.modsym.modsym import ModularSymbols
from sage.modular.modsym.space import ModularSymbolsSpace
from sage.modular.quatalg.brandt import BrandtModule
from sage.modules.free_module import is_FreeModule
from sage.modules.free_module_element import vector
from sage.rings.fast_arith import prime_range
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.ring import is_Ring
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.sets.primes import Primes
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp_method, richcmp_not_equal, rich_to_bool
from sage.structure.sequence import Sequence, Sequence_generic

lazy_import('sage.databases.cremona',
            ['cremona_letter_code', 'CremonaDatabase'])


from . import homspace
from . import lseries
from .morphism import HeckeOperator, Morphism, DegeneracyMap
from .torsion_subgroup import RationalTorsionSubgroup, QQbarTorsionSubgroup
from .finite_subgroup import (FiniteSubgroup_lattice, FiniteSubgroup,
                              TorsionPoint)
from .cuspidal_subgroup import (CuspidalSubgroup, RationalCuspidalSubgroup,
                                RationalCuspSubgroup)


def is_ModularAbelianVariety(x) -> bool:
    """
    Return True if x is a modular abelian variety.

    INPUT:

    -  ``x`` - object

    EXAMPLES::

        sage: from sage.modular.abvar.abvar import is_ModularAbelianVariety
        sage: is_ModularAbelianVariety(5)
        False
        sage: is_ModularAbelianVariety(J0(37))
        True

    Returning True is a statement about the data type not whether or
    not some abelian variety is modular::

        sage: is_ModularAbelianVariety(EllipticCurve('37a'))
        False
    """
    return isinstance(x, ModularAbelianVariety_abstract)


@richcmp_method
class ModularAbelianVariety_abstract(Parent):
    def __init__(self, groups, base_field, is_simple=None, newform_level=None,
                 isogeny_number=None, number=None, check=True):
        """
        Abstract base class for modular abelian varieties.

        INPUT:


        -  ``groups`` - a tuple of congruence subgroups

        -  ``base_field`` - a field

        -  ``is_simple`` - bool; whether or not self is
           simple

        -  ``newform_level`` - if self is isogenous to a
           newform abelian variety, returns the level of that abelian variety

        -  ``isogeny_number`` - which isogeny class the
           corresponding newform is in; this corresponds to the Cremona letter
           code

        -  ``number`` - the t number of the degeneracy map that
           this abelian variety is the image under

        -  ``check`` - whether to do some type checking on the
           defining data


        EXAMPLES: One should not create an instance of this class, but we
        do so anyways here as an example::

            sage: A = sage.modular.abvar.abvar.ModularAbelianVariety_abstract((Gamma0(37),), QQ)
            sage: type(A)
            <class 'sage.modular.abvar.abvar.ModularAbelianVariety_abstract_with_category'>


        All hell breaks loose if you try to do anything with `A`::

            sage: A
            <repr(<sage.modular.abvar.abvar.ModularAbelianVariety_abstract_with_category at 0x...>) failed: NotImplementedError: BUG -- lattice method must be defined in derived class>


        All instances of this class are in the category of modular
        abelian varieties::

            sage: A.category()
            Category of modular abelian varieties over Rational Field
            sage: J0(23).category()
            Category of modular abelian varieties over Rational Field
        """
        if check:
            if not isinstance(groups, tuple):
                raise TypeError("groups must be a tuple")
            for G in groups:
                if not is_CongruenceSubgroup(G):
                    raise TypeError("each element of groups must be a congruence subgroup")
        self.__groups = groups
        if is_simple is not None:
            self.__is_simple = is_simple
        if newform_level is not None:
            self.__newform_level = newform_level
        if number is not None:
            self.__degen_t = number
        if isogeny_number is not None:
            self.__isogeny_number = isogeny_number
        if check and not is_Ring(base_field) and base_field.is_field():
            raise TypeError("base_field must be a field")
        Parent.__init__(self, base=base_field,
                        category=ModularAbelianVarieties(base_field))

    def groups(self):
        r"""
        Return an ordered tuple of the congruence subgroups that the
        ambient product Jacobian is attached to.

        Every modular abelian variety is a finite quotient of an abelian
        subvariety of a product of modular Jacobians `J_\Gamma`.
        This function returns a tuple containing the groups
        `\Gamma`.

        EXAMPLES::

            sage: A = (J0(37) * J1(13))[0]; A
            Simple abelian subvariety 13aG1(1,13) of dimension 2 of J0(37) x J1(13)
            sage: A.groups()
            (Congruence Subgroup Gamma0(37), Congruence Subgroup Gamma1(13))
        """
        return self.__groups

    def is_J0(self) -> bool:
        """
        Return whether of not self is of the form J0(N).

        OUTPUT: bool

        EXAMPLES::

            sage: J0(23).is_J0()
            True
            sage: J1(11).is_J0()
            False
            sage: (J0(23) * J1(11)).is_J0()
            False
            sage: J0(37)[0].is_J0()
            False
            sage: (J0(23) * J0(21)).is_J0()
            False
        """
        return len(self.groups()) == 1 and is_Gamma0(self.groups()[0]) \
            and self.is_ambient()

    def is_J1(self) -> bool:
        """
        Return whether of not self is of the form J1(N).

        OUTPUT: bool

        EXAMPLES::

            sage: J1(23).is_J1()
            True
            sage: J0(23).is_J1()
            False
            sage: (J1(11) * J1(13)).is_J1()
            False
            sage: (J1(11) * J0(13)).is_J1()
            False
            sage: J1(23)[0].is_J1()
            False
        """
        return len(self.groups()) == 1 and is_Gamma1(self.groups()[0]) \
            and self.is_ambient()

    ##########################################################################
    # lattice() *must* be defined by every derived class!!!!
    def lattice(self):
        r"""
        Return lattice in ambient cuspidal modular symbols product that
        defines this modular abelian variety.

        This must be defined in each derived class.

        OUTPUT: a free module over `\ZZ`

        EXAMPLES::

            sage: A = sage.modular.abvar.abvar.ModularAbelianVariety_abstract((Gamma0(37),), QQ)
            sage: A
            <repr(<sage.modular.abvar.abvar.ModularAbelianVariety_abstract_with_category at 0x...>) failed: NotImplementedError: BUG -- lattice method must be defined in derived class>
        """
        raise NotImplementedError("BUG -- lattice method must be defined in derived class")
    ###########################################################################

    def free_module(self):
        r"""
        Synonym for ``self.lattice()``.

        OUTPUT: a free module over `\ZZ`

        EXAMPLES::

            sage: J0(37).free_module()
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
            sage: J0(37)[0].free_module()
            Free module of degree 4 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1 -1  1  0]
            [ 0  0  2 -1]
        """
        return self.lattice()

    def vector_space(self):
        r"""
        Return vector space corresponding to the modular abelian variety.

        This is the lattice tensored with `\QQ`.

        EXAMPLES::

            sage: J0(37).vector_space()
            Vector space of dimension 4 over Rational Field
            sage: J0(37)[0].vector_space()
            Vector space of degree 4 and dimension 2 over Rational Field
            Basis matrix:
            [   1   -1    0  1/2]
            [   0    0    1 -1/2]
        """
        try:
            return self.__vector_space
        except AttributeError:
            self.__vector_space = self.lattice().change_ring(QQ)
            return self.__vector_space

    def base_field(self):
        r"""
        Synonym for ``self.base_ring()``.

        EXAMPLES::

            sage: J0(11).base_field()
            Rational Field
        """
        return self.base_ring()

    def base_extend(self, K):
        """
        EXAMPLES::

            sage: A = J0(37); A
            Abelian variety J0(37) of dimension 2
            sage: A.base_extend(QQbar)
            Abelian variety J0(37) over Algebraic Field of dimension 2
            sage: A.base_extend(GF(7))
            Abelian variety J0(37) over Finite Field of size 7 of dimension 2
        """
        N = self.__newform_level if hasattr(self, '__newform_level') else None
        return ModularAbelianVariety(self.groups(), self.lattice(), K, newform_level=N)

    def __contains__(self, x) -> bool:
        """
        Determine whether or not self contains x.

        EXAMPLES::

            sage: J = J0(67); G = (J[0] + J[1]).intersection(J[1] + J[2])
            sage: G[0]
            Finite subgroup with invariants [5, 10] over QQbar of Abelian subvariety of dimension 3 of J0(67)
            sage: a = G[0].0; a
            [(1/10, 1/10, 3/10, 1/2, 1, -2, -3, 33/10, 0, -1/2)]
            sage: a in J[0]
            False
            sage: a in (J[0]+J[1])
            True
            sage: a in (J[1]+J[2])
            True
            sage: C = G[1]   # abelian variety in kernel
            sage: G[0].0
            [(1/10, 1/10, 3/10, 1/2, 1, -2, -3, 33/10, 0, -1/2)]
            sage: 5*G[0].0
            [(1/2, 1/2, 3/2, 5/2, 5, -10, -15, 33/2, 0, -5/2)]
            sage: 5*G[0].0 in C
            True
        """
        if not isinstance(x, TorsionPoint):
            return False
        if x.parent().abelian_variety().groups() != self.groups():
            return False
        v = x.element()
        n = v.denominator()
        nLambda = self.ambient_variety().lattice().scale(n)
        return n * v in self.lattice() + nLambda

    def __richcmp__(self, other, op):
        """
        Compare two modular abelian varieties.

        If other is not a modular abelian variety, compares the types of
        self and other. If other is a modular abelian variety, compares the
        groups, then if those are the same, compares the newform level and
        isogeny class number and degeneracy map numbers. If those are not
        defined or matched up, compare the underlying lattices.

        EXAMPLES::

            sage: J0(37)[0] < J0(37)[1]
            True
            sage: J0(37)[0] == J0(37)[1]
            False
            sage: J0(33)[0] < J0(33)[1]
            True
            sage: J0(33)[0] >= J0(33)[1]
            False
        """
        if self is other:
            return rich_to_bool(op, 0)
        if not isinstance(other, ModularAbelianVariety_abstract):
            return NotImplemented

        lx = self.groups()
        rx = other.groups()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        try:
            lx = self.__newform_level
            rx = other.__newform_level
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)
        except AttributeError:
            pass

        try:
            lx = self.__isogeny_number
            rx = other.__isogeny_number
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)
        except AttributeError:
            pass

        try:
            lx = self.__degen_t
            rx = other.__degen_t
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)
        except AttributeError:
            pass

        # NOTE!! having the same newform level, isogeny class number,
        # and degen_t does not imply two abelian varieties are equal.
        # See the docstring for self.label.
        return self.lattice()._echelon_matrix_richcmp(other.lattice(), op)

    def __radd__(self, other):
        """
        Return other + self when other is 0. Otherwise raise a TypeError.

        EXAMPLES::

            sage: int(0) + J0(37)
            Abelian variety J0(37) of dimension 2
        """
        if other == 0:
            return self
        raise TypeError

    def _repr_(self) -> str:
        """
        Return string representation of this modular abelian variety.

        This is just the generic base class, so it's unlikely to be called
        in practice.

        EXAMPLES::

            sage: A = J0(23)
            sage: import sage.modular.abvar.abvar as abvar
            sage: abvar.ModularAbelianVariety_abstract._repr_(A)
            'Abelian variety J0(23) of dimension 2'

        ::

            sage: (J0(11) * J0(33))._repr_()
            'Abelian variety J0(11) x J0(33) of dimension 4'
        """
        field = '' if self.base_field() == QQ else ' over %s' % self.base_field()
        # if self.newform_level(none_if_not_known=True) is None:
        simple = self.is_simple(none_if_not_known=True)
        if simple and self.dimension() > 0:
            label = self.label() + ' '
        else:
            label = ''
        simple = 'Simple a' if simple else 'A'
        if self.is_ambient():
            return '%sbelian variety %s%s of dimension %s' % (simple, self._ambient_repr(), field, self.dimension())

        if self.is_subvariety_of_ambient_jacobian():
            sub = 'subvariety'
        else:
            sub = 'variety factor'
        return "%sbelian %s %sof dimension %s of %s%s" % (
            simple, sub, label, self.dimension(), self._ambient_repr(), field)

    def label(self) -> str:
        r"""
        Return the label associated to this modular abelian variety.

        The format of the label is [level][isogeny class][group](t, ambient
        level)

        If this abelian variety `B` has the above label, this
        implies only that `B` is isogenous to the newform abelian
        variety `A_f` associated to the newform with label
        [level][isogeny class][group]. The [group] is empty for
        `\Gamma_0(N)`, is G1 for `\Gamma_1(N)` and is
        GH[...] for `\Gamma_H(N)`.

        .. warning::

           The sum of `\delta_s(A_f)` for all `s\mid t`
           contains `A`, but no sum for a proper divisor of
           `t` contains `A`. It need *not* be the case
           that `B` is equal to `\delta_t(A_f)`!!!

        OUTPUT: string

        EXAMPLES::

            sage: J0(11).label()
            '11a(1,11)'
            sage: J0(11)[0].label()
            '11a(1,11)'
            sage: J0(33)[2].label()
            '33a(1,33)'
            sage: J0(22).label()
            Traceback (most recent call last):
            ...
            ValueError: self must be simple

        We illustrate that self need not equal `\delta_t(A_f)`::

            sage: J = J0(11); phi = J.degeneracy_map(33, 1) + J.degeneracy_map(33,3)
            sage: B = phi.image(); B
            Abelian subvariety of dimension 1 of J0(33)
            sage: B.decomposition()
            [
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33)
            ]
            sage: C = J.degeneracy_map(33,3).image(); C
            Abelian subvariety of dimension 1 of J0(33)
            sage: C == B
            False
        """
        degen = str(self.degen_t()).replace(' ', '')
        return '%s%s' % (self.newform_label(), degen)

    def newform(self, names=None):
        """
        Return the newform `f` such that this abelian variety is isogenous to
        the newform abelian variety `A_f`. If this abelian variety is not
        simple, raise a ``ValueError``.

        INPUT:

        - ``names`` -- (default: None) If the newform has coefficients in a
          number field, then a generator name must be specified.

        OUTPUT: A newform `f` so that self is isogenous to `A_f`.

        EXAMPLES::

            sage: J0(11).newform()
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6)

            sage: f = J0(23).newform(names='a')
            sage: AbelianVariety(f) == J0(23)
            True

            sage: J = J0(33)
            sage: [s.newform('a') for s in J.decomposition()]
            [q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6),
             q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6),
             q + q^2 - q^3 - q^4 - 2*q^5 + O(q^6)]

        The following fails since `J_0(33)` is not simple::

            sage: J0(33).newform()
            Traceback (most recent call last):
            ...
            ValueError: self must be simple
        """
        return Newform(self.newform_label(), names=names)

    def newform_decomposition(self, names=None):
        """
        Return the newforms of the simple subvarieties in the decomposition of
        self as a product of simple subvarieties, up to isogeny.

        OUTPUT:

        - an array of newforms

        EXAMPLES::

            sage: J = J1(11) * J0(23)
            sage: J.newform_decomposition('a')
            [q - 2*q^2 - q^3 + 2*q^4 + q^5 + O(q^6),
            q + a0*q^2 + (-2*a0 - 1)*q^3 + (-a0 - 1)*q^4 + 2*a0*q^5 + O(q^6)]
        """
        return [S.newform(names=names) for S in self.decomposition()]

    def newform_label(self):
        """
        Return the label [level][isogeny class][group] of the newform
        `f` such that this abelian variety is isogenous to the
        newform abelian variety `A_f`. If this abelian variety is
        not simple, raise a ValueError.

        OUTPUT: string

        EXAMPLES::

            sage: J0(11).newform_label()
            '11a'
            sage: J0(33)[2].newform_label()
            '33a'

        The following fails since `J_0(33)` is not simple::

            sage: J0(33).newform_label()
            Traceback (most recent call last):
            ...
            ValueError: self must be simple
        """
        N, G = self.newform_level()
        if is_Gamma0(G):
            group = ''
        elif is_Gamma1(G):
            group = 'G1'
        elif is_GammaH(G):
            group = 'GH%s' % (str(G._generators_for_H()).replace(' ', ''))
        return '%s%s%s' % (N, cremona_letter_code(self.isogeny_number()), group)

    def elliptic_curve(self):
        """
        Return an elliptic curve isogenous to self. If self is not dimension 1
        with rational base ring, raise a ValueError.

        The elliptic curve is found by looking it up in the CremonaDatabase.
        The CremonaDatabase contains all curves up to some large conductor.  If
        a curve is not found in the CremonaDatabase, a RuntimeError will be
        raised. In practice, only the most committed users will see this
        RuntimeError.

        OUTPUT: an elliptic curve isogenous to self.

        EXAMPLES::

            sage: J = J0(11)
            sage: J.elliptic_curve()
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field

            sage: J = J0(49)
            sage: J.elliptic_curve()
            Elliptic Curve defined by y^2 + x*y = x^3 - x^2 - 2*x - 1 over Rational Field

            sage: A = J0(37)[1]
            sage: E = A.elliptic_curve()
            sage: A.lseries()(1)
            0.725681061936153
            sage: E.lseries()(1)
            0.725681061936153

        Elliptic curves are of dimension 1. ::

            sage: J = J0(23)
            sage: J.elliptic_curve()
            Traceback (most recent call last):
            ...
            ValueError: self must be of dimension 1

        This is only implemented for curves over QQ. ::

            sage: J = J0(11).change_ring(CC)
            sage: J.elliptic_curve()
            Traceback (most recent call last):
            ...
            ValueError: base ring must be QQ
        """

        if self.dimension() > 1:
            raise ValueError("self must be of dimension 1")
        if self.base_ring() != QQ:
            raise ValueError("base ring must be QQ")

        f = self.newform('a')
        N = f.level()

        c = CremonaDatabase()
        if N > c.largest_conductor():
            raise RuntimeError("Elliptic curve not found" +
                               " in installed database")

        isogeny_classes = c.isogeny_classes(N)
        curves = [EllipticCurve(x[0][0]) for x in isogeny_classes]

        if len(curves) == 1:
            return curves[0]
        for p in Primes():
            for E in curves:
                if E.ap(p) != f.coefficient(p):
                    curves.remove(E)
                    if len(curves) == 1:
                        return curves[0]

    def _isogeny_to_newform_abelian_variety(self):
        r"""
        Return an isogeny from self to an abelian variety `A_f`
        attached to a newform. If self is not simple (so that no such
        isogeny exists), raise a ValueError.

        EXAMPLES::

            sage: J0(22)[0]._isogeny_to_newform_abelian_variety()
            Abelian variety morphism:
              From: Simple abelian subvariety 11a(1,22) of dimension 1 of J0(22)
              To:   Newform abelian subvariety 11a of dimension 1 of J0(11)
            sage: J = J0(11); phi = J.degeneracy_map(33, 1) + J.degeneracy_map(33,3)
            sage: A = phi.image()
            sage: A._isogeny_to_newform_abelian_variety().matrix()
            [-3  3]
            [ 0 -3]
        """
        try:
            return self._newform_isogeny
        except AttributeError:
            pass

        if not self.is_simple():
            raise ValueError("self is not simple")

        t, N = self.decomposition()[0].degen_t()
        A = self.ambient_variety()
        for i in range(len(self.groups())):
            g = self.groups()[i]
            if N == g.level():
                J = g.modular_abelian_variety()
                d = J.degeneracy_map(self.newform_level()[0], t)
                p = A.project_to_factor(i)
                mat = p.matrix() * d.matrix()
                if not (self.lattice().matrix() * mat).is_zero():
                    break

        from .constructor import AbelianVariety
        Af = AbelianVariety(self.newform_label())
        H = A.Hom(Af.ambient_variety())
        m = H(Morphism(H, mat))
        self._newform_isogeny = m.restrict_domain(self).restrict_codomain(Af)
        return self._newform_isogeny

    def _simple_isogeny(self, other):
        """
        Given self and other, if both are simple, and correspond to the
        same newform with the same congruence subgroup, return an isogeny.
        Otherwise, raise a ValueError.

        INPUT:


        -  ``self, other`` - modular abelian varieties


        OUTPUT: an isogeny

        EXAMPLES::

            sage: J = J0(33); J
            Abelian variety J0(33) of dimension 3
            sage: J[0]._simple_isogeny(J[1])
            Abelian variety morphism:
              From: Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
              To:   Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33)

        The following illustrates how simple isogeny is only implemented
        when the ambients are the same::

            sage: J[0]._simple_isogeny(J1(11))
            Traceback (most recent call last):
            ...
            NotImplementedError: _simple_isogeny only implemented when both abelian variety have the same ambient product Jacobian
        """
        if not is_ModularAbelianVariety(other):
            raise TypeError("other must be a modular abelian variety")

        if not self.is_simple():
            raise ValueError("self is not simple")

        if not other.is_simple():
            raise ValueError("other is not simple")

        if self.groups() != other.groups():
            # The issue here is that the stuff below probably won't make any sense at all if we don't know
            # that the two newform abelian varieties $A_f$ are identical.
            raise NotImplementedError("_simple_isogeny only implemented when both abelian variety have the same ambient product Jacobian")

        if (self.newform_level() != other.newform_level()) or \
           (self.isogeny_number() != other.isogeny_number()):
            raise ValueError("self and other do not correspond to the same newform")

        return other._isogeny_to_newform_abelian_variety().complementary_isogeny() * \
            self._isogeny_to_newform_abelian_variety()

    def _Hom_(self, B, cat=None):
        """
        INPUT:


        -  ``B`` - modular abelian varieties

        -  ``cat`` - category


        EXAMPLES::

            sage: J0(37)._Hom_(J1(37))
            Space of homomorphisms from Abelian variety J0(37) of dimension 2 to Abelian variety J1(37) of dimension 40
            sage: J0(37)._Hom_(J1(37)).homset_category()
            Category of modular abelian varieties over Rational Field
        """
        if cat is None:
            K = self.base_field()
            L = B.base_field()
            if K == L:
                F = K
            elif isinstance(K, sage.rings.abc.AlgebraicField) or isinstance(L, sage.rings.abc.AlgebraicField):
                from sage.rings.qqbar import QQbar as F
            else:
                # TODO -- improve this
                raise ValueError("please specify a category")
            cat = ModularAbelianVarieties(F)
        if self is B:
            return self.endomorphism_ring(cat)
        else:
            return homspace.Homspace(self, B, cat)

    def in_same_ambient_variety(self, other):
        """
        Return True if self and other are abelian subvarieties of the same
        ambient product Jacobian.

        EXAMPLES::

            sage: A,B,C = J0(33)
            sage: A.in_same_ambient_variety(B)
            True
            sage: A.in_same_ambient_variety(J0(11))
            False
        """
        if not is_ModularAbelianVariety(other):
            return False
        if self.groups() != other.groups():
            return False
        if not self.is_subvariety_of_ambient_jacobian() or not other.is_subvariety_of_ambient_jacobian():
            return False
        return True

    def modular_kernel(self):
        """
        Return the modular kernel of this abelian variety, which is the
        kernel of the canonical polarization of self.

        EXAMPLES::

            sage: A = AbelianVariety('33a'); A
            Newform abelian subvariety 33a of dimension 1 of J0(33)
            sage: A.modular_kernel()
            Finite subgroup with invariants [3, 3] over QQ of Newform abelian subvariety 33a of dimension 1 of J0(33)
        """
        try:
            return self.__modular_kernel
        except AttributeError:
            _, f, _ = self.dual()
            G = f.kernel()[0]
            self.__modular_kernel = G
            return G

    def modular_degree(self):
        """
        Return the modular degree of this abelian variety, which is the
        square root of the degree of the modular kernel.

        EXAMPLES::

            sage: A = AbelianVariety('37a')
            sage: A.modular_degree()
            2
        """
        n = self.modular_kernel().order()
        return ZZ(n.sqrt())

    def intersection(self, other):
        """
        Return the intersection of self and other inside a common ambient
        Jacobian product.

        When ``other`` is a modular abelian variety, the output will be a tuple
        ``(G, A)``, where ``G`` is a finite subgroup that surjects onto the
        component group and ``A`` is the identity component. So in particular,
        the intersection is the variety ``G+A``. Note that ``G`` is not chosen
        in any canonical way. When ``other`` is a finite group, the
        intersection will be returned as a finite group.

        INPUT:


        -  ``other`` - a modular abelian variety or a finite
           group


        OUTPUT: If other is a modular abelian variety:


        -  ``G`` - finite subgroup of self

        -  ``A`` - abelian variety (identity component of
           intersection)

        If other is a finite group:

        -  ``G`` - a finite group


        EXAMPLES: We intersect some abelian varieties with finite
        intersection.

        ::

            sage: J = J0(37)
            sage: J[0].intersection(J[1])
            (Finite subgroup with invariants [2, 2] over QQ of Simple abelian subvariety 37a(1,37) of dimension 1 of J0(37), Simple abelian subvariety of dimension 0 of J0(37))

        ::

            sage: D = list(J0(65)); D
            [Simple abelian subvariety 65a(1,65) of dimension 1 of J0(65), Simple abelian subvariety 65b(1,65) of dimension 2 of J0(65), Simple abelian subvariety 65c(1,65) of dimension 2 of J0(65)]
            sage: D[0].intersection(D[1])
            (Finite subgroup with invariants [2] over QQ of Simple abelian subvariety 65a(1,65) of dimension 1 of J0(65), Simple abelian subvariety of dimension 0 of J0(65))
            sage: (D[0]+D[1]).intersection(D[1]+D[2])
            (Finite subgroup with invariants [2] over QQbar of Abelian subvariety of dimension 3 of J0(65), Abelian subvariety of dimension 2 of J0(65))

        ::

            sage: J = J0(33)
            sage: J[0].intersection(J[1])
            (Finite subgroup with invariants [5] over QQ of Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33), Simple abelian subvariety of dimension 0 of J0(33))

        Next we intersect two abelian varieties with non-finite
        intersection::

            sage: J = J0(67); D = J.decomposition(); D
            [
            Simple abelian subvariety 67a(1,67) of dimension 1 of J0(67),
            Simple abelian subvariety 67b(1,67) of dimension 2 of J0(67),
            Simple abelian subvariety 67c(1,67) of dimension 2 of J0(67)
            ]
            sage: (D[0] + D[1]).intersection(D[1] + D[2])
            (Finite subgroup with invariants [5, 10] over QQbar of Abelian subvariety of dimension 3 of J0(67), Abelian subvariety of dimension 2 of J0(67))

        When the intersection is infinite, the output is ``(G, A)``, where
        ``G`` surjects onto the component group. This choice of ``G`` is not
        canonical (see :trac:`26189`). In this following example, ``B`` is a
        subvariety of ``J``::

            sage: d1 = J0(11).degeneracy_map(22, 1)
            sage: d2 = J0(11).degeneracy_map(22, 2)
            sage: B = (d1-d2).image()
            sage: J = J0(22)
            sage: J.intersection(B)
            (Finite subgroup with invariants [] over QQbar of Abelian variety J0(22) of dimension 2,
             Abelian subvariety of dimension 1 of J0(22))
            sage: G, B = B.intersection(J); G, B
            (Finite subgroup with invariants [2] over QQbar of Abelian subvariety of dimension 1 of J0(22),
             Abelian subvariety of dimension 1 of J0(22))
            sage: G.is_subgroup(B)
            True
        """
        # First check whether we are intersecting an abelian variety
        # with a finite subgroup.  If so, call the intersection method
        # for the finite group, which does know how to intersect with
        # an abelian variety.
        if isinstance(other, FiniteSubgroup):
            return other.intersection(self)

        # Now both self and other are abelian varieties.  We require
        # at least that the ambient Jacobian product is the same for
        # them.
        if not self.in_same_ambient_variety(other):
            raise TypeError("other must be an abelian variety in the same ambient space")

        # 1. Compute the abelian variety (connected) part of the intersection
        V = self.vector_space().intersection(other.vector_space())
        if V.dimension() > 0:
            # If there is a nonzero abelian variety, get the actual
            # lattice that defines it.  We intersect (=saturate) in
            # the sum of the lattices, to ensure that the intersection
            # is an abelian subvariety of both self and other (even if
            # they aren't subvarieties of the ambient Jacobian).
            lattice = V.intersection(self.lattice() + other.lattice())
            A = ModularAbelianVariety(self.groups(), lattice, self.base_field(), check=False)
        else:
            A = self.zero_subvariety()

        # 2. Compute the finite intersection group when the
        # intersection is finite, or a group that maps surjectively
        # onto the component group in general.

        # First we get basis matrices for the lattices that define
        # both abelian varieties.
        L = self.lattice().basis_matrix()
        M = other.lattice().basis_matrix()

        # Then we stack matrices and find a subset that forms a
        # basis.
        LM = L.stack(M)
        P = LM.pivot_rows()
        V = (ZZ**L.ncols()).span_of_basis([LM.row(p) for p in P])
        S = (self.lattice() + other.lattice()).saturation()
        n = self.lattice().rank()
        # Finally we project onto the L factor.
        gens = [L.linear_combination_of_rows(v.list()[:n])
                for v in V.coordinate_module(S).basis()]

        if A.dimension() > 0:
            from sage.rings.qqbar import QQbar as finitegroup_base_field
        else:
            finitegroup_base_field = self.base_field()
        G = self.finite_subgroup(gens, field_of_definition=finitegroup_base_field)

        return G, A

    def __add__(self, other):
        r"""
        Return the sum of the *images* of self and other inside the
        ambient Jacobian product.

        Here self and other must be abelian
        subvarieties of the ambient Jacobian product.

        .. WARNING::

            The sum of course only makes sense in some ambient variety,
            and by definition this function takes the sum of the images
            of both self and other in the ambient product Jacobian.

        EXAMPLES:

        We compute the sum of two abelian varieties of `J_0(33)`::

            sage: J = J0(33)
            sage: J[0] + J[1]
            Abelian subvariety of dimension 2 of J0(33)

        We sum all three and get the full `J_0(33)`::

            sage: (J[0] + J[1]) + (J[1] + J[2])
            Abelian variety J0(33) of dimension 3

        Adding to zero works::

            sage: J[0] + 0
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)

        Hence the sum command works::

            sage: sum([J[0], J[2]])
            Abelian subvariety of dimension 2 of J0(33)

        We try to add something in `J_0(33)` to something in
        `J_0(11)`; this shouldn't and doesn't work.

        ::

            sage: J[0] + J0(11)
            Traceback (most recent call last):
            ...
            TypeError: sum not defined since ambient spaces different

        We compute the diagonal image of `J_0(11)` in
        `J_0(33)`, then add the result to the new elliptic curve
        of level `33`.

        ::

            sage: A = J0(11)
            sage: B = (A.degeneracy_map(33,1) + A.degeneracy_map(33,3)).image()
            sage: B + J0(33)[2]
            Abelian subvariety of dimension 2 of J0(33)

        TESTS:

        This exposed a bug in HNF (see :trac:`4527`)::

            sage: A = J0(206).new_subvariety().decomposition()[3] ; A # long time
            Simple abelian subvariety 206d(1,206) of dimension 4 of J0(206)
            sage: B = J0(206).old_subvariety(2) ; B # long time
            Abelian subvariety of dimension 16 of J0(206)
            sage: A+B # long time
            Abelian subvariety of dimension 20 of J0(206)
        """
        if not is_ModularAbelianVariety(other):
            if other == 0:
                return self
            raise TypeError("other must be a modular abelian variety")
        if self.groups() != other.groups():
            raise ValueError("incompatible ambient Jacobians")
        L = self.vector_space() + other.vector_space()
        M = L.intersection(self._ambient_lattice())
        return ModularAbelianVariety(self.groups(), M, self.base_field(), check=False)

    def direct_product(self, other):
        """
        Compute the direct product of self and other.

        INPUT:


        -  ``self, other`` - modular abelian varieties


        OUTPUT: abelian variety

        EXAMPLES::

            sage: J0(11).direct_product(J1(13))
            Abelian variety J0(11) x J1(13) of dimension 3
            sage: A = J0(33)[0].direct_product(J0(33)[1]); A
            Abelian subvariety of dimension 2 of J0(33) x J0(33)
            sage: A.lattice()
            Free module of degree 12 and rank 4 over Integer Ring
            Echelon basis matrix:
            [ 1  1 -2  0  2 -1  0  0  0  0  0  0]
            [ 0  3 -2 -1  2  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  1  0  0  0 -1  2]
            [ 0  0  0  0  0  0  0  1 -1  1  0 -2]
        """
        return self * other

    def __pow__(self, n):
        """
        Return `n^{th}` power of self.

        INPUT:


        -  ``n`` - a nonnegative integer


        OUTPUT: an abelian variety

        EXAMPLES::

            sage: J = J0(37)
            sage: J^0
            Simple abelian subvariety of dimension 0 of J0(37)
            sage: J^1
            Abelian variety J0(37) of dimension 2
            sage: J^1 is J
            True
        """
        n = ZZ(n)
        if n < 0:
            raise ValueError("n must be nonnegative")
        if n == 0:
            return self.zero_subvariety()
        if n == 1:
            return self
        groups = self.groups() * n
        L = self.lattice().basis_matrix()
        lattice = block_diagonal_matrix([L] * n).row_module(ZZ)
        return ModularAbelianVariety(groups, lattice, self.base_field(), check=False)

    def __mul__(self, other):
        """
        Compute the direct product of self and other.

        EXAMPLES: Some modular Jacobians::

            sage: J0(11) * J0(33)
            Abelian variety J0(11) x J0(33) of dimension 4
            sage: J0(11) * J0(33) * J0(11)
            Abelian variety J0(11) x J0(33) x J0(11) of dimension 5

        We multiply some factors of `J_0(65)`::

            sage: d = J0(65).decomposition()
            sage: d[0] * d[1] * J0(11)
            Abelian subvariety of dimension 4 of J0(65) x J0(65) x J0(11)
        """
        if not is_ModularAbelianVariety(other):
            raise TypeError("other must be a modular abelian variety")
        if other.base_ring() != self.base_ring():
            raise TypeError("self and other must have the same base ring")
        groups = tuple(list(self.groups()) + list(other.groups()))
        lattice = self.lattice().direct_sum(other.lattice())
        base_field = self.base_ring()
        return ModularAbelianVariety(groups, lattice, base_field, check=False)

    def quotient(self, other, **kwds):
        """
        Compute the quotient of self and other, where other is either an
        abelian subvariety of self or a finite subgroup of self.

        INPUT:


        -  ``other`` - a finite subgroup or subvariety
        -  further named arguments, that are currently ignored.


        OUTPUT: a pair (A, phi) with phi the quotient map from self to A

        EXAMPLES: We quotient `J_0(33)` out by an abelian
        subvariety::

            sage: Q, f = J0(33).quotient(J0(33)[0])
            sage: Q
            Abelian variety factor of dimension 2 of J0(33)
            sage: f
            Abelian variety morphism:
              From: Abelian variety J0(33) of dimension 3
              To:   Abelian variety factor of dimension 2 of J0(33)

        We quotient `J_0(33)` by the cuspidal subgroup::

            sage: C = J0(33).cuspidal_subgroup()
            sage: Q, f = J0(33).quotient(C)
            sage: Q
            Abelian variety factor of dimension 3 of J0(33)
            sage: f.kernel()[0]
            Finite subgroup with invariants [10, 10] over QQ of Abelian variety J0(33) of dimension 3
            sage: C
            Finite subgroup with invariants [10, 10] over QQ of Abelian variety J0(33) of dimension 3
            sage: J0(11).direct_product(J1(13))
            Abelian variety J0(11) x J1(13) of dimension 3
        """
        return self / other

    def __truediv__(self, other):
        """
        Compute the quotient of self and other, where other is either an
        abelian subvariety of self or a finite subgroup of self.

        INPUT:


        -  ``other`` - a finite subgroup or subvariety


        EXAMPLES: Quotient out by a finite group::

            sage: J = J0(67); G = (J[0] + J[1]).intersection(J[1] + J[2])
            sage: Q, _ = J/G[0]; Q
            Abelian variety factor of dimension 5 of J0(67) over Algebraic Field
            sage: Q.base_field()
            Algebraic Field
            sage: Q.lattice()
            Free module of degree 10 and rank 10 over Integer Ring
            Echelon basis matrix:
            [1/10 1/10 3/10  1/2    0    0    0 3/10    0  1/2]
            [   0  1/5  4/5  4/5    0    0    0    0    0  3/5]
            ...

        Quotient out by an abelian subvariety::

            sage: A, B, C = J0(33)
            sage: Q, phi = J0(33)/A
            sage: Q
            Abelian variety factor of dimension 2 of J0(33)
            sage: phi.domain()
            Abelian variety J0(33) of dimension 3
            sage: phi.codomain()
            Abelian variety factor of dimension 2 of J0(33)
            sage: phi.kernel()
            (Finite subgroup with invariants [2] over QQbar of Abelian variety J0(33) of dimension 3,
             Abelian subvariety of dimension 1 of J0(33))
            sage: phi.kernel()[1] == A
            True

        The abelian variety we quotient out by must be an abelian
        subvariety.

        ::

            sage: Q = (A + B)/C; Q
            Traceback (most recent call last):
            ...
            TypeError: other must be a subgroup or abelian subvariety
        """
        if isinstance(other, FiniteSubgroup):
            if other.abelian_variety() != self:
                other = self.finite_subgroup(other)
            return self._quotient_by_finite_subgroup(other)
        elif isinstance(other, ModularAbelianVariety_abstract) and other.is_subvariety(self):
            return self._quotient_by_abelian_subvariety(other)
        else:
            raise TypeError("other must be a subgroup or abelian subvariety")

    def degeneracy_map(self, M_ls, t_ls):
        """
        Return the degeneracy map with domain self and given
        level/parameter. If self.ambient_variety() is a product of
        Jacobians (as opposed to a single Jacobian), then one can provide a
        list of new levels and parameters, corresponding to the ambient
        Jacobians in order. (See the examples below.)

        INPUT:


        -  ``M, t`` - integers level and `t`, or

        -  ``Mlist, tlist`` - if self is in a nontrivial
           product ambient Jacobian, input consists of a list of levels and
           corresponding list of `t`'s.


        OUTPUT: a degeneracy map

        EXAMPLES: We make several degeneracy maps related to
        `J_0(11)` and `J_0(33)` and compute their
        matrices.

        ::

            sage: d1 = J0(11).degeneracy_map(33, 1); d1
            Degeneracy map from Abelian variety J0(11) of dimension 1 to Abelian variety J0(33) of dimension 3 defined by [1]
            sage: d1.matrix()
            [ 0 -3  2  1 -2  0]
            [ 1 -2  0  1  0 -1]
            sage: d2 = J0(11).degeneracy_map(33, 3); d2
            Degeneracy map from Abelian variety J0(11) of dimension 1 to Abelian variety J0(33) of dimension 3 defined by [3]
            sage: d2.matrix()
            [-1  0  0  0  1 -2]
            [-1 -1  1 -1  1  0]
            sage: d3 = J0(33).degeneracy_map(11, 1); d3
            Degeneracy map from Abelian variety J0(33) of dimension 3 to Abelian variety J0(11) of dimension 1 defined by [1]

        He we verify that first mapping from level `11` to level
        `33`, then back is multiplication by `4`::

            sage: d1.matrix() * d3.matrix()
            [4 0]
            [0 4]

        We compute a more complicated degeneracy map involving nontrivial
        product ambient Jacobians; note that this is just the block direct
        sum of the two matrices at the beginning of this example::

            sage: d = (J0(11)*J0(11)).degeneracy_map([33,33], [1,3]); d
            Degeneracy map from Abelian variety J0(11) x J0(11) of dimension 2 to Abelian variety J0(33) x J0(33) of dimension 6 defined by [1, 3]
            sage: d.matrix()
            [ 0 -3  2  1 -2  0  0  0  0  0  0  0]
            [ 1 -2  0  1  0 -1  0  0  0  0  0  0]
            [ 0  0  0  0  0  0 -1  0  0  0  1 -2]
            [ 0  0  0  0  0  0 -1 -1  1 -1  1  0]
        """
        if not isinstance(M_ls, list):
            M_ls = [M_ls]
        if not isinstance(t_ls, list):
            t_ls = [t_ls]

        groups = self.groups()
        length = len(M_ls)
        if length != len(t_ls):
            raise ValueError("must have same number of Ms and ts")
        if length != len(groups):
            raise ValueError("must have same number of Ms and groups in ambient variety")

        for i in range(length):
            N = groups[i].level()
            if (M_ls[i] % N) and (N % M_ls[i]):
                raise ValueError("one level must divide the other in %s-th component" % i)
            if ((max(M_ls[i], N) // min(M_ls[i], N)) % t_ls[i]):
                raise ValueError("each t must divide the quotient of the levels")

        ls = [self.groups()[i].modular_abelian_variety().degeneracy_map(M_ls[i], t_ls[i]).matrix()
              for i in range(length)]

        new_codomain = prod([self.groups()[i]._new_group_from_level(M_ls[i]).modular_abelian_variety()
                             for i in range(length)])
        M = block_diagonal_matrix(ls, subdivide=False)

        H = self.Hom(new_codomain)
        return H(DegeneracyMap(H, M.restrict_domain(self.lattice()), t_ls))

    def _quotient_by_finite_subgroup(self, G):
        """
        Return the quotient of self by the finite subgroup `G`.
        This is used internally by the quotient and __div__ commands.

        INPUT:


        -  ``G`` - a finite subgroup of self


        OUTPUT: abelian variety - the quotient `Q` of self by `G`


        -  ``morphism`` - from self to the quotient
           `Q`


        EXAMPLES: We quotient the elliptic curve `J_0(11)` out by
        its cuspidal subgroup.

        ::

            sage: A = J0(11)
            sage: G = A.cuspidal_subgroup(); G
            Finite subgroup with invariants [5] over QQ of Abelian variety J0(11) of dimension 1
            sage: Q, f = A._quotient_by_finite_subgroup(G)
            sage: Q
            Abelian variety factor of dimension 1 of J0(11)
            sage: f
            Abelian variety morphism:
              From: Abelian variety J0(11) of dimension 1
              To:   Abelian variety factor of dimension 1 of J0(11)

        We compute the finite kernel of `f` (hence the [0]) and
        note that it equals the subgroup `G` that we quotiented out
        by::

            sage: f.kernel()[0] == G
            True

        Quotienting by the identity should return the original variety and the
        identity morphism. :trac:`6392` ::

            sage: J = J0(22)
            sage: G = J.zero_subgroup()
            sage: Q, f = J / G
            sage: Q
            Abelian variety J0(22) of dimension 2
            sage: f
            Abelian variety endomorphism of Abelian variety J0(22) of dimension 2

        """
        if G.order() == 1:
            return self, self.endomorphism_ring().identity()
        L = self.lattice() + G.lattice()
        A = ModularAbelianVariety(self.groups(), L, G.field_of_definition())
        M = L.coordinate_module(self.lattice()).basis_matrix()
        phi = self.Hom(A)(M)
        return A, phi

    def _quotient_by_abelian_subvariety(self, B):
        """
        Return the quotient of self by the abelian variety `B`.
        This is used internally by the quotient and __div__ commands.

        INPUT:


        -  ``B`` - an abelian subvariety of self


        OUTPUT:


        -  ``abelian variety`` - quotient `Q` of self
           by B

        -  ``morphism`` - from self to the quotient
           `Q`


        EXAMPLES: We compute the new quotient of `J_0(33)`.

        ::

            sage: A = J0(33); B = A.old_subvariety()
            sage: Q, f = A._quotient_by_abelian_subvariety(B)

        Note that the quotient happens to also be an abelian subvariety::

            sage: Q
            Abelian subvariety of dimension 1 of J0(33)
            sage: Q.lattice()
            Free module of degree 6 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0  0 -1  0  0]
            [ 0  0  1  0  1 -1]
            sage: f
            Abelian variety morphism:
              From: Abelian variety J0(33) of dimension 3
              To:   Abelian subvariety of dimension 1 of J0(33)

        We verify that `B` is equal to the kernel of the quotient
        map.

        ::

            sage: f.kernel()[1] == B
            True

        Next we quotient `J_0(33)` out by `Q` itself::

            sage: C, g = A._quotient_by_abelian_subvariety(Q)

        The result is not a subvariety::

            sage: C
            Abelian variety factor of dimension 2 of J0(33)
            sage: C.lattice()
            Free module of degree 6 and rank 4 over Integer Ring
            Echelon basis matrix:
            [ 1/3    0    0  2/3   -1    0]
            [   0    1    0    0   -1    1]
            [   0    0  1/3    0 -2/3  2/3]
            [   0    0    0    1   -1   -1]
        """

        # We first compute the complement of B in self to get
        # an abelian variety C also in self such that self/B
        # is isogenous to C. This is the case because the
        # projection map pi:self --> C is surjective and has
        # kernel a finite extension of the abelian variety B.
        C = B.complement(self)

        # Now that we have C we need to find some abelian variety Q
        # isogenous to C and a map self --> Q whose kernel is exactly
        # B.  We do this by computing the kernel of the map pi below,
        # which is an extension of the abelian variety B by a finite
        # group Phi of complements.  Our strategy is to enlarge the
        # lattice that defines C so that the map pi below suddenly
        # has connected kernel.

        pi = self.projection(C)
        psi = pi.factor_out_component_group()
        Q = psi.codomain()
        return Q, psi

    def projection(self, A, check=True):
        """
        Given an abelian subvariety A of self, return a projection morphism
        from self to A. Note that this morphism need not be unique.

        INPUT:


        -  ``A`` - an abelian variety


        OUTPUT: a morphism

        EXAMPLES::

            sage: a,b,c = J0(33)
            sage: pi = J0(33).projection(a); pi.matrix()
            [ 3 -2]
            [-5  5]
            [-4  1]
            [ 3 -2]
            [ 5  0]
            [ 1  1]
            sage: pi = (a+b).projection(a); pi.matrix()
            [ 0  0]
            [-3  2]
            [-4  1]
            [-1 -1]
            sage: pi = a.projection(a); pi.matrix()
            [1 0]
            [0 1]

        We project onto a factor in a product of two Jacobians::

            sage: A = J0(11)*J0(11); A
            Abelian variety J0(11) x J0(11) of dimension 2
            sage: A[0]
            Simple abelian subvariety 11a(1,11) of dimension 1 of J0(11) x J0(11)
            sage: A.projection(A[0])
            Abelian variety morphism:
              From: Abelian variety J0(11) x J0(11) of dimension 2
              To:   Simple abelian subvariety 11a(1,11) of dimension 1 of J0(11) x J0(11)
            sage: A.projection(A[0]).matrix()
            [0 0]
            [0 0]
            [1 0]
            [0 1]
            sage: A.projection(A[1]).matrix()
            [1 0]
            [0 1]
            [0 0]
            [0 0]
        """
        if check and not A.is_subvariety(self):
            raise ValueError("A must be an abelian subvariety of self")

        W = A.complement(self)
        mat = A.lattice().basis_matrix().stack(W.lattice().basis_matrix())

        # solve  X * mat = self, i.e. write each row of self in terms of the
        # rows of mat.
        X = mat.solve_left(self.lattice().basis_matrix())

        # The projection map is got from the first 2*dim(A) columns of X.
        X = X.matrix_from_columns(range(2 * A.dimension()))

        X, _ = X._clear_denom()

        return Morphism(self.Hom(A), X)

    def project_to_factor(self, n):
        """
        If self is an ambient product of Jacobians, return a projection
        from self to the nth such Jacobian.

        EXAMPLES::

            sage: J = J0(33)
            sage: J.project_to_factor(0)
            Abelian variety endomorphism of Abelian variety J0(33) of dimension 3

        ::

            sage: J = J0(33) * J0(37) * J0(11)
            sage: J.project_to_factor(2)
            Abelian variety morphism:
              From: Abelian variety J0(33) x J0(37) x J0(11) of dimension 6
              To:   Abelian variety J0(11) of dimension 1
            sage: J.project_to_factor(2).matrix()
            [0 0]
            [0 0]
            [0 0]
            [0 0]
            [0 0]
            [0 0]
            [0 0]
            [0 0]
            [0 0]
            [0 0]
            [1 0]
            [0 1]
        """
        if not self.is_ambient():
            raise ValueError("self is not ambient")
        if n >= len(self.groups()):
            raise IndexError("index (=%s) too large (max = %s)" % (n, len(self.groups())))

        G = self.groups()[n]
        A = G.modular_abelian_variety()
        index = sum([gp.modular_symbols().cuspidal_subspace().dimension()
                     for gp in self.groups()[0: n]])

        H = self.Hom(A)
        mat = H.matrix_space()(0)
        mat.set_block(index, 0, identity_matrix(2 * A.dimension()))

        return H(Morphism(H, mat))

    def is_subvariety_of_ambient_jacobian(self):
        """
        Return True if self is (presented as) a subvariety of the ambient
        product Jacobian.

        Every abelian variety in Sage is a quotient of a subvariety of an
        ambient Jacobian product by a finite subgroup.

        EXAMPLES::

            sage: J0(33).is_subvariety_of_ambient_jacobian()
            True
            sage: A = J0(33)[0]; A
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: A.is_subvariety_of_ambient_jacobian()
            True
            sage: B, phi = A / A.torsion_subgroup(2)
            sage: B
            Abelian variety factor of dimension 1 of J0(33)
            sage: phi.matrix()
            [2 0]
            [0 2]
            sage: B.is_subvariety_of_ambient_jacobian()
            False
        """
        try:
            return self.__is_sub_ambient
        except AttributeError:
            self.__is_sub_ambient = (self.lattice().denominator() == 1)
            return self.__is_sub_ambient

    def ambient_variety(self):
        """
        Return the ambient modular abelian variety that contains this
        abelian variety. The ambient variety is always a product of
        Jacobians of modular curves.

        OUTPUT: abelian variety

        EXAMPLES::

            sage: A = J0(33)[0]; A
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: A.ambient_variety()
            Abelian variety J0(33) of dimension 3
        """
        try:
            return self.__ambient_variety
        except AttributeError:
            A = ModularAbelianVariety(self.groups(), ZZ**(2 * self._ambient_dimension()),
                                      self.base_field(), check=False)
            self.__ambient_variety = A
            return A

    def ambient_morphism(self):
        """
        Return the morphism from self to the ambient variety. This is
        injective if self is natural a subvariety of the ambient product
        Jacobian.

        OUTPUT: morphism

        The output is cached.

        EXAMPLES: We compute the ambient structure morphism for an abelian
        subvariety of `J_0(33)`::

            sage: A,B,C = J0(33)
            sage: phi = A.ambient_morphism()
            sage: phi.domain()
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: phi.codomain()
            Abelian variety J0(33) of dimension 3
            sage: phi.matrix()
            [ 1  1 -2  0  2 -1]
            [ 0  3 -2 -1  2  0]

        phi is of course injective

        ::

            sage: phi.kernel()
            (Finite subgroup with invariants [] over QQ of Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33),
             Abelian subvariety of dimension 0 of J0(33))

        This is the same as the basis matrix for the lattice corresponding
        to self::

            sage: A.lattice()
            Free module of degree 6 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  1 -2  0  2 -1]
            [ 0  3 -2 -1  2  0]

        We compute a non-injective map to an ambient space::

            sage: Q,pi = J0(33)/A
            sage: phi = Q.ambient_morphism()
            sage: phi.matrix()
            [  1   4   1   9  -1  -1]
            [  0  15   0   0  30 -75]
            [  0   0   5  10  -5  15]
            [  0   0   0  15 -15  30]
            sage: phi.kernel()[0]
            Finite subgroup with invariants [5, 15, 15] over QQ of Abelian variety factor of dimension 2 of J0(33)
        """
        try:
            return self.__ambient_morphism
        except AttributeError:
            matrix, _ = self.lattice().basis_matrix()._clear_denom()
            phi = Morphism(self.Hom(self.ambient_variety()), matrix)
            self.__ambient_morphism = phi
            return phi

    def is_ambient(self) -> bool:
        """
        Return True if self equals the ambient product Jacobian.

        OUTPUT: bool

        EXAMPLES::

            sage: A,B,C = J0(33)
            sage: A.is_ambient()
            False
            sage: J0(33).is_ambient()
            True
            sage: (A+B).is_ambient()
            False
            sage: (A+B+C).is_ambient()
            True
        """
        try:
            return self.__is_ambient
        except AttributeError:
            pass
        L = self.lattice()
        self.__is_ambient = (self.lattice() == ZZ**L.degree())
        return self.__is_ambient

    def dimension(self):
        """
        Return the dimension of this abelian variety.

        EXAMPLES::

            sage: A = J0(23)
            sage: A.dimension()
            2
        """
        return self.lattice().rank() // 2

    def conductor(self):
        """
        Return the conductor of this abelian variety.

        EXAMPLES::

            sage: A = J0(23)
            sage: A.conductor().factor()
            23^2

            sage: A = J1(25)
            sage: A.conductor().factor()
            5^24

            sage: A = J0(11^2); A.decomposition()
            [
            Simple abelian subvariety 11a(1,121) of dimension 1 of J0(121),
            Simple abelian subvariety 11a(11,121) of dimension 1 of J0(121),
            Simple abelian subvariety 121a(1,121) of dimension 1 of J0(121),
            Simple abelian subvariety 121b(1,121) of dimension 1 of J0(121),
            Simple abelian subvariety 121c(1,121) of dimension 1 of J0(121),
            Simple abelian subvariety 121d(1,121) of dimension 1 of J0(121)
            ]
            sage: A.conductor().factor()
            11^10

            sage: A = J0(33)[0]; A
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: A.conductor()
            11
            sage: A.elliptic_curve().conductor()
            11
        """
        if not self.base_ring() == QQ:
            raise ValueError("base ring must be QQ")
        return prod(f.level() ** f.base_ring().degree()
                    for f in self.newform_decomposition('a'))

    def rank(self):
        """
        Return the rank of the underlying lattice of self.

        EXAMPLES::

            sage: J = J0(33)
            sage: J.rank()
            6
            sage: J[1]
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33)
            sage: (J[1] * J[1]).rank()
            4
        """
        return self.lattice().rank()

    def degree(self):
        """
        Return the degree of this abelian variety, which is the dimension
        of the ambient Jacobian product.

        EXAMPLES::

            sage: A = J0(23)
            sage: A.dimension()
            2
        """
        return self._ambient_dimension()

    def endomorphism_ring(self, category=None):
        """
        Return the endomorphism ring of self.

        OUTPUT: b = self.sturm_bound()

        EXAMPLES: We compute a few endomorphism rings::

            sage: J0(11).endomorphism_ring()
            Endomorphism ring of Abelian variety J0(11) of dimension 1
            sage: J0(37).endomorphism_ring()
            Endomorphism ring of Abelian variety J0(37) of dimension 2
            sage: J0(33)[2].endomorphism_ring()
            Endomorphism ring of Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)

        No real computation is done::

            sage: J1(123456).endomorphism_ring()
            Endomorphism ring of Abelian variety J1(123456) of dimension 423185857
        """
        try:
            return self.__endomorphism_ring
        except AttributeError:
            pass

        self.__endomorphism_ring = homspace.EndomorphismSubring(self, category=category)
        return self.__endomorphism_ring

    def sturm_bound(self):
        r"""
        Return a bound `B` such that all Hecke operators
        `T_n` for `n\leq B` generate the Hecke algebra.

        OUTPUT: integer

        EXAMPLES::

            sage: J0(11).sturm_bound()
            2
            sage: J0(33).sturm_bound()
            8
            sage: J1(17).sturm_bound()
            48
            sage: J1(123456).sturm_bound()
            1693483008
            sage: JH(37,[2,3]).sturm_bound()
            7
            sage: J1(37).sturm_bound()
            228
        """
        try:
            return self.__sturm_bound
        except AttributeError:
            B = max([G.sturm_bound(2) for G in self.groups()])
            self.__sturm_bound = B
            return B

    def is_hecke_stable(self) -> bool:
        """
        Return True if self is stable under the Hecke operators of its
        ambient Jacobian.

        OUTPUT: bool

        EXAMPLES::

            sage: J0(11).is_hecke_stable()
            True
            sage: J0(33)[2].is_hecke_stable()
            True
            sage: J0(33)[0].is_hecke_stable()
            False
            sage: (J0(33)[0] + J0(33)[1]).is_hecke_stable()
            True
        """
        try:
            return self._is_hecke_stable
        except AttributeError:
            pass

        # b = self.modular_symbols().sturm_bound()
        b = max([m.sturm_bound()
                 for m in self._ambient_modular_symbols_spaces()])
        J = self.ambient_variety()
        L = self.lattice()
        B = self.lattice().basis()

        for n in prime_range(1, b + 1):
            Tn_matrix = J.hecke_operator(n).matrix()
            for v in B:
                if v * Tn_matrix not in L:
                    self._is_hecke_stable = False
                    return False

        self._is_hecke_stable = True
        return True

    def is_subvariety(self, other) -> bool:
        """
        Return True if self is a subvariety of other as they sit in a
        common ambient modular Jacobian. In particular, this function will
        only return True if self and other have exactly the same ambient
        Jacobians.

        EXAMPLES::

            sage: J = J0(37); J
            Abelian variety J0(37) of dimension 2
            sage: A = J[0]; A
            Simple abelian subvariety 37a(1,37) of dimension 1 of J0(37)
            sage: A.is_subvariety(A)
            True
            sage: A.is_subvariety(J)
            True
        """
        if not is_ModularAbelianVariety(other):
            return False
        if self is other:
            return True
        if self.groups() != other.groups():
            return False
        L = self.lattice()
        M = other.lattice()
        # self is an abelian subvariety of other if and only if
        #   1. L is a subset of M (so the abelian subvarieties of
        #      the ambient J are equal), and
        #   2. L is relatively saturated in M, i.e., M/L is
        #      torsion free.
        if not L.is_submodule(M):
            return False
        # To determine if L is relatively saturated we compute the
        # intersection of M with (L tensor Q) and see if that equals
        # L.
        return L.change_ring(QQ).intersection(M) == L

    def change_ring(self, R):
        """
        Change the base ring of this modular abelian variety.

        EXAMPLES::

            sage: A = J0(23)
            sage: A.change_ring(QQ)
            Abelian variety J0(23) of dimension 2
        """
        return ModularAbelianVariety(self.groups(), self.lattice(), R, check=False)

    def level(self):
        """
        Return the level of this modular abelian variety, which is an
        integer N (usually minimal) such that this modular abelian variety
        is a quotient of `J_1(N)`. In the case that the ambient
        variety of self is a product of Jacobians, return the LCM of their
        levels.

        EXAMPLES::

            sage: J1(5077).level()
            5077
            sage: JH(389,[4]).level()
            389
            sage: (J0(11)*J0(17)).level()
            187
        """
        try:
            return self.__level
        except AttributeError:
            self.__level = LCM([G.level() for G in self.groups()])
            return self.__level

    def newform_level(self, none_if_not_known=False):
        """
        Write self as a product (up to isogeny) of newform abelian
        varieties `A_f`. Then this function return the least
        common multiple of the levels of the newforms `f`, along
        with the corresponding group or list of groups (the groups do not
        appear with multiplicity).

        INPUT:


        -  ``none_if_not_known`` - (default: False) if True,
           return None instead of attempting to compute the newform level, if
           it isn't already known. This None result is not cached.


        OUTPUT: integer group or list of distinct groups

        EXAMPLES::

            sage: J0(33)[0].newform_level()
            (11, Congruence Subgroup Gamma0(33))
            sage: J0(33)[0].newform_level(none_if_not_known=True)
            (11, Congruence Subgroup Gamma0(33))

        Here there are multiple groups since there are in fact multiple
        newforms::

            sage: (J0(11) * J1(13)).newform_level()
            (143, [Congruence Subgroup Gamma0(11), Congruence Subgroup Gamma1(13)])
        """
        try:
            return self.__newform_level
        except AttributeError:
            if none_if_not_known:
                return None
            level = LCM([f.level() for f in self.newform_decomposition('a')])
            groups = sorted(set([f.group() for f in
                                 self.newform_decomposition('a')]))
            if len(groups) == 1:
                groups = groups[0]
            self.__newform_level = level, groups
            return self.__newform_level

    def zero_subvariety(self):
        """
        Return the zero subvariety of self.

        EXAMPLES::

            sage: J = J0(37)
            sage: J.zero_subvariety()
            Simple abelian subvariety of dimension 0 of J0(37)
            sage: J.zero_subvariety().level()
            37
            sage: J.zero_subvariety().newform_level()
            (1, [])
        """
        try:
            return self.__zero_subvariety
        except AttributeError:
            lattice = (ZZ**(2 * self.degree())).zero_submodule()
            A = ModularAbelianVariety(self.groups(), lattice, self.base_field(),
                                      is_simple=True, check=False)
            self.__zero_subvariety = A
            return A

    # #########################################################################
    # Properties of the ambient product of Jacobians
    # #########################################################################
    def _ambient_repr(self):
        """
        OUTPUT: string

        EXAMPLES::

            sage: (J0(33)*J1(11))._ambient_repr()
            'J0(33) x J1(11)'
        """
        v = []
        for G in self.groups():
            if is_Gamma0(G):
                v.append('J0(%s)' % G.level())
            elif is_Gamma1(G):
                v.append('J1(%s)' % G.level())
            elif is_GammaH(G):
                v.append('JH(%s,%s)' % (G.level(), G._generators_for_H()))
        return ' x '.join(v)

    def _ambient_latex_repr(self):
        """
        Return Latex representation of the ambient product.

        OUTPUT: string

        EXAMPLES::

            sage: (J0(11) * J0(33))._ambient_latex_repr()
            'J_0(11) \\times J_0(33)'
        """
        v = []
        for G in self.groups():
            if is_Gamma0(G):
                v.append('J_0(%s)' % G.level())
            elif is_Gamma1(G):
                v.append('J_1(%s)' % G.level())
            elif is_GammaH(G):
                v.append('J_H(%s,%s)' % (G.level(), G._generators_for_H()))
        return ' \\times '.join(v)

    def _ambient_lattice(self):
        """
        Return free lattice of rank twice the degree of self. This is the
        lattice corresponding to the ambient product Jacobian.

        OUTPUT: lattice

        EXAMPLES: We compute the ambient lattice of a product::

            sage: (J0(33)*J1(11))._ambient_lattice()
            Ambient free module of rank 8 over the principal ideal domain Integer Ring

        We compute the ambient lattice of an abelian subvariety
        `J_0(33)`, which is the same as the lattice for the
        `J_0(33)` itself::

            sage: A = J0(33)[0]; A._ambient_lattice()
            Ambient free module of rank 6 over the principal ideal domain Integer Ring
            sage: J0(33)._ambient_lattice()
            Ambient free module of rank 6 over the principal ideal domain Integer Ring
        """
        try:
            return self.__ambient_lattice
        except AttributeError:
            self.__ambient_lattice = ZZ**(2 * self.degree())
            return self.__ambient_lattice

    def _ambient_modular_symbols_spaces(self):
        """
        Return a tuple of the ambient cuspidal modular symbols spaces that
        make up the Jacobian product that contains self.

        OUTPUT: tuple of cuspidal modular symbols spaces

        EXAMPLES::

            sage: (J0(11) * J0(33))._ambient_modular_symbols_spaces()
            (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field,
             Modular Symbols subspace of dimension 6 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field)
            sage: (J0(11) * J0(33)[0])._ambient_modular_symbols_spaces()
            (Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field,
             Modular Symbols subspace of dimension 6 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field)
        """
        if not self.is_ambient():
            return self.ambient_variety()._ambient_modular_symbols_spaces()
        try:
            return self.__ambient_modular_symbols_spaces
        except AttributeError:
            X = tuple([ModularSymbols(G).cuspidal_subspace() for G in self.groups()])
            self.__ambient_modular_symbols_spaces = X
            return X

    def _ambient_modular_symbols_abvars(self):
        """
        Return a tuple of the ambient modular symbols abelian varieties
        that make up the Jacobian product that contains self.

        OUTPUT: tuple of modular symbols abelian varieties

        EXAMPLES::

            sage: (J0(11) * J0(33))._ambient_modular_symbols_abvars()
            (Abelian variety J0(11) of dimension 1, Abelian variety J0(33) of dimension 3)
        """
        if not self.is_ambient():
            return self.ambient_variety()._ambient_modular_symbols_abvars()
        try:
            return self.__ambient_modular_symbols_abvars
        except AttributeError:
            X = tuple([ModularAbelianVariety_modsym(M) for M in self._ambient_modular_symbols_spaces()])
            self.__ambient_modular_symbols_abvars = X
            return X

    def _ambient_dimension(self):
        """
        Return the dimension of the ambient Jacobian product.

        EXAMPLES::

            sage: A = J0(37) * J1(13); A
            Abelian variety J0(37) x J1(13) of dimension 4
            sage: A._ambient_dimension()
            4
            sage: B = A[0]; B
            Simple abelian subvariety 13aG1(1,13) of dimension 2 of J0(37) x J1(13)
            sage: B._ambient_dimension()
            4

        This example is fast because it implicitly calls
        _ambient_dimension.

        ::

            sage: J0(902834082394)
            Abelian variety J0(902834082394) of dimension 113064825881
        """
        try:
            return self.__ambient_dimension
        except AttributeError:
            d = sum([G.dimension_cusp_forms(2) for G in self.groups()], Integer(0))
            self.__ambient_dimension = d
            return d

    def _ambient_hecke_matrix_on_modular_symbols(self, n):
        r"""
        Return block direct sum of the matrix of the Hecke operator
        `T_n` acting on each of the ambient modular symbols
        spaces.

        INPUT:


        -  ``n`` - an integer `\geq 1`.


        OUTPUT: a matrix

        EXAMPLES::

            sage: (J0(11) * J1(13))._ambient_hecke_matrix_on_modular_symbols(2)
            [-2  0  0  0  0  0]
            [ 0 -2  0  0  0  0]
            [ 0  0 -1 -1 -1  1]
            [ 0  0  1 -2 -1  0]
            [ 0  0  0  0 -2  1]
            [ 0  0  0  0 -1 -1]
        """
        if not self.is_ambient():
            return self.ambient_variety()._ambient_hecke_matrix_on_modular_symbols(n)
        try:
            return self.__ambient_hecke_matrix_on_modular_symbols[n]
        except AttributeError:
            self.__ambient_hecke_matrix_on_modular_symbols = {}
        except KeyError:
            pass
        M = self._ambient_modular_symbols_spaces()
        if not M:
            return matrix(QQ, 0)
        T = M[0].hecke_matrix(n)
        for i in range(1, len(M)):
            T = T.block_sum(M[i].hecke_matrix(n))
        self.__ambient_hecke_matrix_on_modular_symbols[n] = T
        return T

    def rational_torsion_order(self, proof=True):
        """
        Return the order of the rational torsion subgroup of this modular
        abelian variety.

        This function is really an alias for
        :meth:`~sage.modular.abvar.torsion_subgroup.RationalTorsionSubgroup.order`
        See the docstring there for a more in-depth reference and more
        interesting examples.

        INPUT:

        - ``proof`` -- a boolean (default: True)

        OUTPUT:

        The order of the rational torsion subgroup of this modular abelian
        variety.


        EXAMPLES::

            sage: J0(11).rational_torsion_subgroup().order()
            5
            sage: J0(11).rational_torsion_order()
            5
        """
        return self.rational_torsion_subgroup().order(proof=proof)

    def number_of_rational_points(self):
        """
        Return the number of rational points of this modular abelian variety.

        It is not always possible to compute the order of the torsion
        subgroup. The BSD conjecture is assumed to compute the algebraic rank.

        OUTPUT: a positive integer or infinity

        EXAMPLES::

            sage: J0(23).number_of_rational_points()
            11
            sage: J0(29).number_of_rational_points()
            7
            sage: J0(37).number_of_rational_points()
            +Infinity

            sage: J0(12); J0(12).number_of_rational_points()
            Abelian variety J0(12) of dimension 0
            1

            sage: J1(17).number_of_rational_points()
            584

            sage: J1(16).number_of_rational_points()
            Traceback (most recent call last):
            ...
            RuntimeError: Unable to compute order of torsion subgroup (it is in [1, 2, 4, 5, 10, 20])
        """
        # Check easy dimension zero case
        if self.dimension() == 0:
            return ZZ(1)

        positive_rank = self.lseries().vanishes_at_1()

        if positive_rank:
            return infinity
        else:
            return self.rational_torsion_subgroup().order()

    def frobenius_polynomial(self, p, var='x'):
        """
        Computes the frobenius polynomial at `p`.

        INPUT:

        - ``p`` -- prime number

        OUTPUT:

        a monic integral polynomial

        EXAMPLES::

            sage: f = Newform('39b','a')
            sage: A = AbelianVariety(f)
            sage: A.frobenius_polynomial(5)
            x^4 + 2*x^2 + 25

            sage: J = J0(23)
            sage: J.frobenius_polynomial(997)
            x^4 + 20*x^3 + 1374*x^2 + 19940*x + 994009

            sage: J = J0(33)
            sage: J.frobenius_polynomial(7)
            x^6 + 9*x^4 - 16*x^3 + 63*x^2 + 343

            sage: J = J0(19)
            sage: J.frobenius_polynomial(3, var='y')
            y^2 + 2*y + 3

            sage: J = J0(3); J
            Abelian variety J0(3) of dimension 0
            sage: J.frobenius_polynomial(11)
            1

            sage: A = J1(27)[1]; A
            Simple abelian subvariety 27bG1(1,27) of dimension 12 of J1(27)
            sage: A.frobenius_polynomial(11)
            x^24 - 3*x^23 - 15*x^22 + 126*x^21 - 201*x^20 - 1488*x^19 + 7145*x^18 - 1530*x^17 - 61974*x^16 + 202716*x^15 - 19692*x^14 - 1304451*x^13 + 4526883*x^12 - 14348961*x^11 - 2382732*x^10 + 269814996*x^9 - 907361334*x^8 - 246408030*x^7 + 12657803345*x^6 - 28996910448*x^5 - 43086135081*x^4 + 297101409066*x^3 - 389061369015*x^2 - 855935011833*x + 3138428376721

            sage: J = J1(33)
            sage: J.frobenius_polynomial(11)
            Traceback (most recent call last):
            ...
            ValueError: p must not divide the level of self
            sage: J.frobenius_polynomial(4)
            Traceback (most recent call last):
            ...
            ValueError: p must be prime
        """
        if self.dimension() == 0:
            return ZZ(1)
        if self.level() % p == 0:
            raise ValueError("p must not divide the level of self")
        if not is_prime(p):
            raise ValueError("p must be prime")
        if not self.is_simple():
            from .constructor import AbelianVariety
            decomp = [AbelianVariety(f) for f in
                      self.newform_decomposition('a')]
            return prod((s.frobenius_polynomial(p) for s in
                         decomp))
        f = self.newform('a')
        Kf = f.base_ring()
        eps = f.character()
        Qe = eps.base_ring()
        d = Kf.degree() / Qe.degree()

        ZZpoly = PolynomialRing(ZZ, var)
        z = ZZpoly.gens()[0]

        if Qe == QQ:
            Gp = self.hecke_polynomial(p)
            return ZZpoly(z**d * Gp(z + p / z))

        # relativize number fields to compute charpoly of
        # left multiplication of ap on Kf as a Qe-vector
        # space.
        Lf = Kf.relativize(Qe.gen(), 'a')
        to_Lf = Lf.structure()[1]

        # Re is just Qe with a different name
        Re = Lf.base_field()
        to_Re = Qe.hom(Re.gens()[0], Re)

        name = Kf._names[0]
        ap = to_Lf(f.modular_symbols(1).eigenvalue(p, name))
        Gp = ap.charpoly(var='x')

        S = PolynomialRing(Re, 'x')
        x = S.gens()[0]
        h = S(x**d * Gp(x=x + to_Re(eps(p)) * p / x))

        # take Qe norm
        R = PolynomialRing(QQ, ['x', 'y'])
        x, y = R.gens()
        g = Qe.defining_polynomial()
        H = sum(h[i].lift() * y**i for i in range(h.degree() + 1))
        ans = g.resultant(H)

        return ZZpoly(ans.univariate_polynomial())

    ###############################################################################
    # Rational and Integral Homology
    ###############################################################################
    def _rational_homology_space(self):
        """
        Return the rational homology of this modular abelian variety.

        EXAMPLES::

            sage: J = J0(11)
            sage: J._rational_homology_space()
            Vector space of dimension 2 over Rational Field

        The result is cached::

            sage: J._rational_homology_space() is J._rational_homology_space()
            True
        """
        try:
            return self.__rational_homology_space
        except AttributeError:
            HQ = self.rational_homology().free_module()
            self.__rational_homology_space = HQ
            return HQ

    def homology(self, base_ring=ZZ):
        """
        Return the homology of this modular abelian variety.

        .. warning::

           For efficiency reasons the basis of the integral homology
           need not be the same as the basis for the rational
           homology.

        EXAMPLES::

            sage: J0(389).homology(GF(7))
            Homology with coefficients in Finite Field of size 7 of Abelian variety J0(389) of dimension 32
            sage: J0(389).homology(QQ)
            Rational Homology of Abelian variety J0(389) of dimension 32
            sage: J0(389).homology(ZZ)
            Integral Homology of Abelian variety J0(389) of dimension 32
        """
        from . import homology
        try:
            return self._homology[base_ring]
        except AttributeError:
            self._homology = {}
        except KeyError:
            pass
        if base_ring == ZZ:
            H = homology.IntegralHomology(self)
        elif base_ring == QQ:
            H = homology.RationalHomology(self)
        else:
            H = homology.Homology_over_base(self, base_ring)
        self._homology[base_ring] = H
        return H

    def integral_homology(self):
        """
        Return the integral homology of this modular abelian variety.

        EXAMPLES::

            sage: H = J0(43).integral_homology(); H
            Integral Homology of Abelian variety J0(43) of dimension 3
            sage: H.rank()
            6
            sage: H = J1(17).integral_homology(); H
            Integral Homology of Abelian variety J1(17) of dimension 5
            sage: H.rank()
            10

        If you just ask for the rank of the homology, no serious
        calculations are done, so the following is fast::

            sage: H = J0(50000).integral_homology(); H
            Integral Homology of Abelian variety J0(50000) of dimension 7351
            sage: H.rank()
            14702

        A product::

            sage: H = (J0(11) * J1(13)).integral_homology()
            sage: H.hecke_operator(2)
            Hecke operator T_2 on Integral Homology of Abelian variety J0(11) x J1(13) of dimension 3
            sage: H.hecke_operator(2).matrix()
            [-2  0  0  0  0  0]
            [ 0 -2  0  0  0  0]
            [ 0  0 -1 -1 -1  1]
            [ 0  0  1 -2 -1  0]
            [ 0  0  0  0 -2  1]
            [ 0  0  0  0 -1 -1]
        """
        return self.homology(ZZ)

    def rational_homology(self):
        """
        Return the rational homology of this modular abelian variety.

        EXAMPLES::

            sage: H = J0(37).rational_homology(); H
            Rational Homology of Abelian variety J0(37) of dimension 2
            sage: H.rank()
            4
            sage: H.base_ring()
            Rational Field
            sage: H = J1(17).rational_homology(); H
            Rational Homology of Abelian variety J1(17) of dimension 5
            sage: H.rank()
            10
            sage: H.base_ring()
            Rational Field
        """
        return self.homology(QQ)

    ###############################################################################
    # L-series
    ###############################################################################
    def lseries(self):
        """
        Return the complex `L`-series of this modular abelian
        variety.

        EXAMPLES::

            sage: A = J0(37)
            sage: A.lseries()
            Complex L-series attached to Abelian variety J0(37) of dimension 2
        """
        try:
            return self.__lseries
        except AttributeError:
            pass
        self.__lseries = lseries.Lseries_complex(self)
        return self.__lseries

    def padic_lseries(self, p):
        """
        Return the `p`-adic `L`-series of this modular
        abelian variety.

        EXAMPLES::

            sage: A = J0(37)
            sage: A.padic_lseries(7)
            7-adic L-series attached to Abelian variety J0(37) of dimension 2
        """
        p = int(p)
        try:
            return self.__lseries_padic[p]
        except AttributeError:
            self.__lseries_padic = {}
        except KeyError:
            pass
        self.__lseries_padic[p] = lseries.Lseries_padic(self, p)
        return self.__lseries_padic[p]

    ###############################################################################
    # Hecke Operators
    ###############################################################################
    def hecke_operator(self, n):
        """
        Return the `n^{th}` Hecke operator on the modular abelian
        variety, if this makes sense [[elaborate]]. Otherwise raise a
        ValueError.

        EXAMPLES: We compute `T_2` on `J_0(37)`.

        ::

            sage: t2 = J0(37).hecke_operator(2); t2
            Hecke operator T_2 on Abelian variety J0(37) of dimension 2
            sage: t2.charpoly().factor()
            x * (x + 2)
            sage: t2.index()
            2

        Note that there is no matrix associated to Hecke operators on
        modular abelian varieties. For a matrix, instead consider, e.g.,
        the Hecke operator on integral or rational homology.

        ::

            sage: t2.action_on_homology().matrix()
            [-1  1  1 -1]
            [ 1 -1  1  0]
            [ 0  0 -2  1]
            [ 0  0  0  0]
        """
        try:
            return self._hecke_operator[n]
        except AttributeError:
            self._hecke_operator = {}
        except KeyError:
            pass
        Tn = HeckeOperator(self, n)
        self._hecke_operator[n] = Tn
        return Tn

    def hecke_polynomial(self, n, var='x'):
        r"""
        Return the characteristic polynomial of the `n^{th}` Hecke
        operator `T_n` acting on self. Raises an ArithmeticError
        if self is not Hecke equivariant.

        INPUT:


        -  ``n`` - integer `\geq 1`

        -  ``var`` - string (default: 'x'); valid variable
           name


        EXAMPLES::

            sage: J0(33).hecke_polynomial(2)
            x^3 + 3*x^2 - 4
            sage: f = J0(33).hecke_polynomial(2, 'y'); f
            y^3 + 3*y^2 - 4
            sage: f.parent()
            Univariate Polynomial Ring in y over Rational Field
            sage: J0(33)[2].hecke_polynomial(3)
            x + 1
            sage: J0(33)[0].hecke_polynomial(5)
            x - 1
            sage: J0(33)[0].hecke_polynomial(11)
            x - 1
            sage: J0(33)[0].hecke_polynomial(3)
            Traceback (most recent call last):
            ...
            ArithmeticError: subspace is not invariant under matrix
        """
        n = Integer(n)
        if n <= 0:
            raise ValueError("n must be a positive integer")
        key = (n, var)
        try:
            return self.__hecke_polynomial[key]
        except AttributeError:
            self.__hecke_polynomial = {}
        except KeyError:
            pass
        f = self._compute_hecke_polynomial(n, var=var)
        self.__hecke_polynomial[key] = f
        return f

    def _compute_hecke_polynomial(self, n, var='x'):
        """
        Return the Hecke polynomial of index `n` in terms of the
        given variable.

        INPUT:


        -  ``n`` - positive integer

        -  ``var`` - string (default: 'x')


        EXAMPLES::

            sage: A = J0(33)*J0(11)
            sage: A._compute_hecke_polynomial(2)
            x^4 + 5*x^3 + 6*x^2 - 4*x - 8
        """
        return self.hecke_operator(n).charpoly(var=var)

    def _integral_hecke_matrix(self, n):
        """
        Return the matrix of the Hecke operator `T_n` acting on
        the integral homology of this modular abelian variety, if the
        modular abelian variety is stable under `T_n`. Otherwise,
        raise an ArithmeticError.

        EXAMPLES::

            sage: A = J0(23)
            sage: t = A._integral_hecke_matrix(2); t
            [ 0  1 -1  0]
            [ 0  1 -1  1]
            [-1  2 -2  1]
            [-1  1  0 -1]
            sage: t.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Integer Ring
        """
        A = self._ambient_hecke_matrix_on_modular_symbols(n)
        return A.restrict(self.lattice())

    def _rational_hecke_matrix(self, n):
        r"""
        Return the matrix of the Hecke operator `T_n` acting on
        the rational homology `H_1(A,\QQ)` of this modular
        abelian variety, if this action is defined. Otherwise, raise an
        ArithmeticError.

        EXAMPLES::

            sage: A = J0(23)
            sage: t = A._rational_hecke_matrix(2); t
            [ 0  1 -1  0]
            [ 0  1 -1  1]
            [-1  2 -2  1]
            [-1  1  0 -1]
            sage: t.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Rational Field
        """
        return self._integral_hecke_matrix(n)

    ###############################################################################
    # Subgroups
    ###############################################################################
    def qbar_torsion_subgroup(self):
        r"""
        Return the group of all points of finite order in the algebraic
        closure of this abelian variety.

        EXAMPLES::

            sage: T = J0(33).qbar_torsion_subgroup(); T
            Group of all torsion points in QQbar on Abelian variety J0(33) of dimension 3

        The field of definition is the same as the base field of the
        abelian variety.

        ::

            sage: T.field_of_definition()
            Rational Field

        On the other hand, T is a module over `\ZZ`.

        ::

            sage: T.base_ring()
            Integer Ring
        """
        try:
            return self.__qbar_torsion_subgroup
        except AttributeError:
            G = QQbarTorsionSubgroup(self)
            self.__qbar_torsion_subgroup = G
            return G

    def rational_torsion_subgroup(self):
        """
        Return the maximal torsion subgroup of self defined over QQ.

        EXAMPLES::

            sage: J = J0(33)
            sage: A = J.new_subvariety()
            sage: A
            Abelian subvariety of dimension 1 of J0(33)
            sage: t = A.rational_torsion_subgroup(); t
            Torsion subgroup of Abelian subvariety of dimension 1 of J0(33)
            sage: t.multiple_of_order()
            4
            sage: t.divisor_of_order()
            4
            sage: t.order()
            4
            sage: t.gens()
            [[(1/2, 0, 0, -1/2, 0, 0)], [(0, 0, 1/2, 0, 1/2, -1/2)]]
        """
        try:
            return self.__rational_torsion_subgroup
        except AttributeError:
            T = RationalTorsionSubgroup(self)
            self.__rational_torsion_subgroup = T
            return T

    def cuspidal_subgroup(self):
        """
        Return the cuspidal subgroup of this modular abelian variety. This
        is the subgroup generated by rational cusps.

        EXAMPLES::

            sage: J = J0(54)
            sage: C = J.cuspidal_subgroup()
            sage: C.gens()
            [[(1/3, 0, 0, 0, 0, 1/3, 0, 2/3)], [(0, 1/3, 0, 0, 0, 2/3, 0, 1/3)], [(0, 0, 1/9, 1/9, 1/9, 1/9, 1/9, 2/9)], [(0, 0, 0, 1/3, 0, 1/3, 0, 0)], [(0, 0, 0, 0, 1/3, 1/3, 0, 1/3)], [(0, 0, 0, 0, 0, 0, 1/3, 2/3)]]
            sage: C.invariants()
            [3, 3, 3, 3, 3, 9]
            sage: J1(13).cuspidal_subgroup()
            Finite subgroup with invariants [19, 19] over QQ of Abelian variety J1(13) of dimension 2
            sage: A = J0(33)[0]
            sage: A.cuspidal_subgroup()
            Finite subgroup with invariants [5] over QQ of Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
        """
        try:
            return self._cuspidal_subgroup
        except AttributeError:
            if not self.is_subvariety_of_ambient_jacobian():
                raise ValueError("self must be a subvariety of the ambient variety")
            if self.is_ambient():
                T = self._ambient_cuspidal_subgroup(rational_only=False)
            else:
                T = self.ambient_variety().cuspidal_subgroup().intersection(self)
            self._cuspidal_subgroup = T
            return T

    def _ambient_cuspidal_subgroup(self, rational_only=False, rational_subgroup=False):
        """
        EXAMPLES::

            sage: (J1(13)*J0(11))._ambient_cuspidal_subgroup()
            Finite subgroup with invariants [19, 95] over QQ of Abelian variety J1(13) x J0(11) of dimension 3
            sage: (J0(33))._ambient_cuspidal_subgroup()
            Finite subgroup with invariants [10, 10] over QQ of Abelian variety J0(33) of dimension 3
            sage: (J0(33)*J0(33))._ambient_cuspidal_subgroup()
            Finite subgroup with invariants [10, 10, 10, 10] over QQ of Abelian variety J0(33) x J0(33) of dimension 6
        """
        n = 2 * self.degree()
        i = 0
        lattice = (ZZ**n).zero_submodule()
        if rational_subgroup:
            CS = RationalCuspidalSubgroup
        elif rational_only:
            CS = RationalCuspSubgroup
        else:
            CS = CuspidalSubgroup
        for J in self._ambient_modular_symbols_abvars():
            L = CS(J).lattice().basis_matrix()
            Z_left = matrix(QQ, L.nrows(), i)
            Z_right = matrix(QQ, L.nrows(), n - i - L.ncols())
            lattice += (Z_left.augment(L).augment(Z_right)).row_module(ZZ)
            i += L.ncols()
        return FiniteSubgroup_lattice(self, lattice, field_of_definition=self.base_field())

    def shimura_subgroup(self):
        r"""
        Return the Shimura subgroup of this modular abelian variety.

        This is the kernel of `J_0(N) \rightarrow J_1(N)` under the
        natural map.

        Here we compute the Shimura subgroup as the kernel of `J_0(N)
        \rightarrow J_0(Np)` where the map is the difference between
        the two degeneracy maps.

        EXAMPLES::

            sage: J = J0(11)
            sage: J.shimura_subgroup()
            Finite subgroup with invariants [5] over QQ of Abelian variety J0(11) of dimension 1

            sage: J = J0(17)
            sage: G = J.cuspidal_subgroup(); G
            Finite subgroup with invariants [4] over QQ of Abelian variety J0(17) of dimension 1
            sage: S = J.shimura_subgroup(); S
            Finite subgroup with invariants [4] over QQ of Abelian variety J0(17) of dimension 1
            sage: G.intersection(S)
            Finite subgroup with invariants [2] over QQ of Abelian variety J0(17) of dimension 1

            sage: J = J0(33)
            sage: A = J.decomposition()[0]
            sage: A.shimura_subgroup()
            Finite subgroup with invariants [5] over QQ of Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: J.shimura_subgroup()
            Finite subgroup with invariants [10] over QQ of Abelian variety J0(33) of dimension 3
        """
        N = self.level()
        J = self.ambient_variety()
        for p in prime_range(100):
            if N % p:
                break
        phi = J.degeneracy_map(N * p, 1)
        phip = J.degeneracy_map(N * p, p)
        SIG = (phi - phip).kernel()
        assert SIG[1].dimension() == 0, "The intersection should have dimension 0"

        return self.intersection(SIG[0])

    def rational_cusp_subgroup(self):
        r"""
        Return the subgroup of this modular abelian variety generated by
        rational cusps.

        This is a subgroup of the group of rational points in the cuspidal
        subgroup.

        .. warning::

           This is only currently implemented for
           `\Gamma_0(N)`.

        EXAMPLES::

            sage: J = J0(54)
            sage: CQ = J.rational_cusp_subgroup(); CQ
            Finite subgroup with invariants [3, 3, 9] over QQ of Abelian variety J0(54) of dimension 4
            sage: CQ.gens()
            [[(1/3, 0, 0, 1/3, 2/3, 1/3, 0, 1/3)], [(0, 0, 1/9, 1/9, 7/9, 7/9, 1/9, 8/9)], [(0, 0, 0, 0, 0, 0, 1/3, 2/3)]]
            sage: factor(CQ.order())
            3^4
            sage: CQ.invariants()
            [3, 3, 9]

        In this example the rational cuspidal subgroup and the cuspidal
        subgroup differ by a lot.

        ::

            sage: J = J0(49)
            sage: J.cuspidal_subgroup()
            Finite subgroup with invariants [2, 14] over QQ of Abelian variety J0(49) of dimension 1
            sage: J.rational_cusp_subgroup()
            Finite subgroup with invariants [2] over QQ of Abelian variety J0(49) of dimension 1

        Note that computation of the rational cusp subgroup isn't
        implemented for `\Gamma_1`.

        ::

            sage: J = J1(13)
            sage: J.cuspidal_subgroup()
            Finite subgroup with invariants [19, 19] over QQ of Abelian variety J1(13) of dimension 2
            sage: J.rational_cusp_subgroup()
            Traceback (most recent call last):
            ...
            NotImplementedError: computation of rational cusps only implemented in Gamma0 case.
        """
        try:
            return self._rational_cusp_subgroup
        except AttributeError:
            if not self.is_subvariety_of_ambient_jacobian():
                raise ValueError("self must be a subvariety of the ambient variety")
            if self.is_ambient():
                T = self._ambient_cuspidal_subgroup(rational_only=True)
            else:
                T = self.ambient_variety().rational_cusp_subgroup().intersection(self)
            self._rational_cusp_subgroup = T
            return T

    def rational_cuspidal_subgroup(self):
        r"""
        Return the rational subgroup of the cuspidal subgroup of this
        modular abelian variety.

        This is a subgroup of the group of rational points in the
        cuspidal subgroup.

        .. warning::

           This is only currently implemented for
           `\Gamma_0(N)`.

        EXAMPLES::

            sage: J = J0(54)
            sage: CQ = J.rational_cuspidal_subgroup(); CQ
            Finite subgroup with invariants [3, 3, 9] over QQ of Abelian variety J0(54) of dimension 4
            sage: CQ.gens()
            [[(1/3, 0, 0, 1/3, 2/3, 1/3, 0, 1/3)], [(0, 0, 1/9, 1/9, 7/9, 7/9, 1/9, 8/9)], [(0, 0, 0, 0, 0, 0, 1/3, 2/3)]]
            sage: factor(CQ.order())
            3^4
            sage: CQ.invariants()
            [3, 3, 9]

        In this example the rational cuspidal subgroup and the cuspidal
        subgroup differ by a lot.

        ::

            sage: J = J0(49)
            sage: J.cuspidal_subgroup()
            Finite subgroup with invariants [2, 14] over QQ of Abelian variety J0(49) of dimension 1
            sage: J.rational_cuspidal_subgroup()
            Finite subgroup with invariants [2] over QQ of Abelian variety J0(49) of dimension 1

        Note that computation of the rational cusp subgroup isn't
        implemented for `\Gamma_1`.

        ::

            sage: J = J1(13)
            sage: J.cuspidal_subgroup()
            Finite subgroup with invariants [19, 19] over QQ of Abelian variety J1(13) of dimension 2
            sage: J.rational_cuspidal_subgroup()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented when group is Gamma0
        """
        try:
            return self._rational_cuspidal_subgroup
        except AttributeError:
            if not self.is_subvariety_of_ambient_jacobian():
                raise ValueError("self must be a subvariety of the ambient variety")
            if self.is_ambient():
                T = self._ambient_cuspidal_subgroup(rational_subgroup=True)
            else:
                T = self.ambient_variety().rational_cuspidal_subgroup().intersection(self)
            self._rational_cuspidal_subgroup = T
            return T

    def zero_subgroup(self):
        """
        Return the zero subgroup of this modular abelian variety, as a
        finite group.

        EXAMPLES::

            sage: A =J0(54); G = A.zero_subgroup(); G
            Finite subgroup with invariants [] over QQ of Abelian variety J0(54) of dimension 4
            sage: G.is_subgroup(A)
            True
        """
        try:
            return self.__zero_subgroup
        except AttributeError:
            G = FiniteSubgroup_lattice(self, self.lattice(), field_of_definition=QQ)
            self.__zero_subgroup = G
            return G

    def finite_subgroup(self, X, field_of_definition=None, check=True):
        """
        Return a finite subgroup of this modular abelian variety.

        INPUT:


        -  ``X`` - list of elements of other finite subgroups
           of this modular abelian variety or elements that coerce into the
           rational homology (viewed as a rational vector space); also X could
           be a finite subgroup itself that is contained in this abelian
           variety.

        -  ``field_of_definition`` - (default: None) field
           over which this group is defined. If None try to figure out the
           best base field.


        OUTPUT: a finite subgroup of a modular abelian variety

        EXAMPLES::

            sage: J = J0(11)
            sage: J.finite_subgroup([[1/5,0], [0,1/3]])
            Finite subgroup with invariants [15] over QQbar of Abelian variety J0(11) of dimension 1

        ::

            sage: J = J0(33); C = J[0].cuspidal_subgroup(); C
            Finite subgroup with invariants [5] over QQ of Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: J.finite_subgroup([[0,0,0,0,0,1/6]])
            Finite subgroup with invariants [6] over QQbar of Abelian variety J0(33) of dimension 3
            sage: J.finite_subgroup(C)
            Finite subgroup with invariants [5] over QQ of Abelian variety J0(33) of dimension 3

        This method gives a way of changing the ambient abelian variety of a
        finite subgroup. This caused an issue in :trac:`6392` but should be
        fixed now. ::

            sage: A, B = J0(43)
            sage: G, _ = A.intersection(B)
            sage: G
            Finite subgroup with invariants [2, 2] over QQ of Simple abelian subvariety 43a(1,43) of dimension 1 of J0(43)
            sage: B.finite_subgroup(G)
            Finite subgroup with invariants [2, 2] over QQ of Simple abelian subvariety 43b(1,43) of dimension 2 of J0(43)

        """
        if isinstance(X, (list, tuple)):
            X = self._ambient_lattice().span(X)
        elif isinstance(X, FiniteSubgroup):
            if field_of_definition is None:
                field_of_definition = X.field_of_definition()
            A = X.abelian_variety()
            if A.groups() != self.groups():
                raise ValueError("ambient product Jacobians must be equal")
            if A == self:
                X = X.lattice()
            else:
                if X.is_subgroup(self):
                    X = (X.lattice() +
                         self.ambient_variety().lattice()).intersection(
                             self.vector_space())
                else:
                    raise ValueError("X must be a subgroup of self.")

        if field_of_definition is None:
            from sage.rings.qqbar import QQbar as field_of_definition

        return FiniteSubgroup_lattice(
            self, X, field_of_definition=field_of_definition, check=check)

    def torsion_subgroup(self, n):
        """
        If n is an integer, return the subgroup of points of order n.
        Return the `n`-torsion subgroup of elements of order
        dividing `n` of this modular abelian variety `A`,
        i.e., the group `A[n]`.

        EXAMPLES::

            sage: J1(13).torsion_subgroup(19)
            Finite subgroup with invariants [19, 19, 19, 19] over QQ of Abelian variety J1(13) of dimension 2

        ::

            sage: A = J0(23)
            sage: G = A.torsion_subgroup(5); G
            Finite subgroup with invariants [5, 5, 5, 5] over QQ of Abelian variety J0(23) of dimension 2
            sage: G.order()
            625
            sage: G.gens()
            [[(1/5, 0, 0, 0)], [(0, 1/5, 0, 0)], [(0, 0, 1/5, 0)], [(0, 0, 0, 1/5)]]
            sage: A = J0(23)
            sage: A.torsion_subgroup(2).order()
            16
        """
        try:
            return self.__torsion_subgroup[n]
        except KeyError:
            pass
        except AttributeError:
            self.__torsion_subgroup = {}
        lattice = self.lattice().scale(1 / Integer(n))
        H = FiniteSubgroup_lattice(self, lattice, field_of_definition=self.base_field())
        self.__torsion_subgroup[n] = H
        return H

    # #########################################################################
    # Decomposition
    # #########################################################################
    def degen_t(self, none_if_not_known=False):
        """
        If this abelian variety is obtained via decomposition then it gets
        labeled with the newform label along with some information about
        degeneracy maps. In particular, the label ends in a pair
        `(t,N)`, where `N` is the ambient level and
        `t` is an integer that divides the quotient of `N`
        by the newform level. This function returns the tuple
        `(t,N)`, or raises a ValueError if self isn't simple.

        .. note::

           It need not be the case that self is literally equal to the
           image of the newform abelian variety under the `t^{th}`
           degeneracy map. See the documentation for the label method
           for more details.

        INPUT:


        -  ``none_if_not_known`` - (default: False) - if
           True, return None instead of attempting to compute the degen map's
           `t`, if it isn't known. This None result is not cached.


        OUTPUT: a pair (integer, integer)

        EXAMPLES::

            sage: D = J0(33).decomposition(); D
            [
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33),
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33),
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
            ]
            sage: D[0].degen_t()
            (1, 33)
            sage: D[1].degen_t()
            (3, 33)
            sage: D[2].degen_t()
            (1, 33)
            sage: J0(33).degen_t()
            Traceback (most recent call last):
            ...
            ValueError: self must be simple
        """
        try:
            return self.__degen_t
        except AttributeError:
            if none_if_not_known:
                return None
            elif self.dimension() > 0 and self.is_simple():
                self.__degen_t = self.decomposition()[0].degen_t()
                return self.__degen_t
            raise ValueError("self must be simple")

    def isogeny_number(self, none_if_not_known=False):
        """
        Return the number (starting at 0) of the isogeny class of new
        simple abelian varieties that self is in. If self is not simple,
        raises a ValueError exception.

        INPUT:


        -  ``none_if_not_known`` - bool (default: False); if
           True then this function may return None instead of True of False if
           we don't already know the isogeny number of self.


        EXAMPLES: We test the none_if_not_known flag first::

            sage: J0(33).isogeny_number(none_if_not_known=True) is None
            True

        Of course, `J_0(33)` is not simple, so this function
        raises a ValueError::

            sage: J0(33).isogeny_number()
            Traceback (most recent call last):
            ...
            ValueError: self must be simple

        Each simple factor has isogeny number 1, since that's the number at
        which the factor is new.

        ::

            sage: J0(33)[1].isogeny_number()
            0
            sage: J0(33)[2].isogeny_number()
            0

        Next consider `J_0(37)` where there are two distinct
        newform factors::

            sage: J0(37)[1].isogeny_number()
            1
        """
        try:
            return self.__isogeny_number
        except AttributeError:
            if none_if_not_known:
                return None
            elif self.is_simple():
                self.__isogeny_number = self.decomposition()[0].isogeny_number()
                return self.__isogeny_number
            else:
                raise ValueError("self must be simple")

    def is_simple(self, none_if_not_known=False):
        """
        Return whether or not this modular abelian variety is simple, i.e.,
        has no proper nonzero abelian subvarieties.

        INPUT:


        -  ``none_if_not_known`` - bool (default: False); if
           True then this function may return None instead of True of False if
           we don't already know whether or not self is simple.


        EXAMPLES::

            sage: J0(5).is_simple(none_if_not_known=True) is None  # this may fail if J0(5) comes up elsewhere...
            True
            sage: J0(33).is_simple()
            False
            sage: J0(33).is_simple(none_if_not_known=True)
            False
            sage: J0(33)[1].is_simple()
            True
            sage: J1(17).is_simple()
            False
        """
        try:
            return self.__is_simple
        except AttributeError:
            if none_if_not_known:
                return None
            self.__is_simple = len(self.decomposition()) <= 1
            return self.__is_simple

    def decomposition(self, simple=True, bound=None):
        """
        Return a sequence of abelian subvarieties of self that are all
        simple, have finite intersection and sum to self.

        INPUT: simple- bool (default: True) if True, all factors are
        simple. If False, each factor returned is isogenous to a power of a
        simple and the simples in each factor are distinct.


        -  ``bound`` - int (default: None) if given, only use
           Hecke operators up to this bound when decomposing. This can give
           wrong answers, so use with caution!


        EXAMPLES::

            sage: m = ModularSymbols(11).cuspidal_submodule()
            sage: d1 = m.degeneracy_map(33,1).matrix(); d3=m.degeneracy_map(33,3).matrix()
            sage: w = ModularSymbols(33).submodule((d1 + d3).image(), check=False)
            sage: A = w.abelian_variety(); A
            Abelian subvariety of dimension 1 of J0(33)
            sage: D = A.decomposition(); D
            [
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33)
            ]
            sage: D[0] == A
            True
            sage: B = A + J0(33)[0]; B
            Abelian subvariety of dimension 2 of J0(33)
            sage: dd = B.decomposition(simple=False); dd
            [
            Abelian subvariety of dimension 2 of J0(33)
            ]
            sage: dd[0] == B
            True
            sage: dd = B.decomposition(); dd
            [
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33),
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33)
            ]
            sage: sum(dd) == B
            True

        We decompose a product of two Jacobians::

            sage: (J0(33) * J0(11)).decomposition()
            [
            Simple abelian subvariety 11a(1,11) of dimension 1 of J0(33) x J0(11),
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33) x J0(11),
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33) x J0(11),
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33) x J0(11)
            ]
        """
        try:
            return self.__decomposition[(simple, bound)]
        except KeyError:
            pass
        except AttributeError:
            self.__decomposition = {}

        if self.is_ambient():
            # Decompose each piece, then lift
            if not self.groups():
                D = []
            elif len(self.groups()) == 1:
                D = ModularAbelianVariety_modsym(ModularSymbols(self.groups()[0], sign=0).cuspidal_submodule()).decomposition(simple=simple, bound=bound)
            else:
                # Decompose each ambient modular symbols factor.
                # X = [ModularAbelianVariety_modsym(ModularSymbols(G,sign=0).cuspidal_submodule()) for G in self.groups()]
                from .abvar_ambient_jacobian import ModAbVar_ambient_jacobian_class
                X = [ModAbVar_ambient_jacobian_class(G) for G in self.groups()]
                E = [A.decomposition(simple=simple, bound=bound) for A in X]
                i = 0
                n = 2 * self.dimension()
                # Now lift each factor of the decomposition to self.
                G = self.groups()
                D = []
                K = self.base_field()
                for C in E:
                    for B in C:
                        L = B.lattice().basis_matrix()
                        if simple:
                            is_simple = True
                        else:
                            is_simple = None
                        lattice = matrix(QQ, L.nrows(), i).augment(L).augment(matrix(QQ, L.nrows(), n - i - L.ncols())).row_module(ZZ)
                        D.append(ModularAbelianVariety(G, lattice, K, is_simple=is_simple, newform_level=B.newform_level(),
                                                       isogeny_number=B.isogeny_number(none_if_not_known=True),
                                                       number=B.degen_t(none_if_not_known=True)))
                    if C:
                        i += L.ncols()
        elif not simple:
            # In this case decompose the ambient space into powers of
            # simple abelian varieties (i.e. with
            # \code{simple=False)}, and then intersect the lattice
            # corresponding to self with each of these factors.
            D = []
            L = self.lattice()
            groups = self.groups()
            K = self.base_ring()
            for X in self.ambient_variety().decomposition(simple=False):
                lattice = L.intersection(X.vector_space())
                if lattice.rank() > 0:
                    the_factor = ModularAbelianVariety(groups, lattice, K, is_simple=X.is_simple(none_if_not_known=True), newform_level=X.newform_level(), isogeny_number=X.isogeny_number(none_if_not_known=True), number=X.degen_t(none_if_not_known=True))
                    D.append(the_factor)

        else:
            # See the documentation for self._classify_ambient_factors
            # in order to understand what we're doing here.
            I_F, I_E, X = self._classify_ambient_factors(simple=simple, bound=bound)
            Z_E = [X[i] for i in I_E]
            Z_F = [X[i] for i in I_F]
            F = sum(Z_F, self.zero_subvariety())
            # Now self is isogenous to the sum of the factors in Z.
            # We use this isogeny to obtain a product decomposition of
            # self.
            if F == self:
                # The easy case -- it is already such a decomposition
                D = Z_F
            else:
                # The hard case -- now we have to pull back the
                # factorization

                # Suppose $B$ is an abelian variety and there is a
                # finite degree map $B\to J$, where $J$ is an ambient
                # Jacobian.  Suppose further that we find abelian
                # subvarieties $E$ and $F$ of $J$ such that $E + F =
                # J$, $E$ and $F$ have finite intersection, the
                # composition $B \to J \to J/E$ is an isogeny, and we
                # know an explicit decomposition of $F$.  Then we can
                # compute a decomposition of $B$ as follows.  Let
                # $L_E$ and $L_F$ be the lattices corresponding to $E$
                # and $F$ inside of $L_J$.  Compute a matrix $\Phi$
                # representing the composition $L_B \to L_J \to L_F
                # \otimes \QQ$, where the map $L_J$ to $L_F\otimes
                # \QQ$ is projection onto the second factor in the
                # decomposition of $L_J$ as $L_E + L_F$ (up to finite
                # index).  Finally, for each factor $A_i$ of $F$ with
                # lattice $L_{A_i}$, compute the saturation $S_i$ of
                # $\Phi^{-1}(L_{A_i})$.  Then the $S_i$ define a
                # decomposition of $B$.
                E = sum(Z_E, self.zero_subvariety())
                L_B = self.lattice()
                L_E = E.lattice()
                L_F = F.lattice()
                decomp_matrix = L_E.basis_matrix().stack(L_F.basis_matrix())
                # Now we compute explicitly the ZZ-linear map (over
                # QQ) from L_B that is "projection onto L_F".  This
                # means write each element of a basis for L_B in terms
                # of decomp_matrix, then take the bottom coordinates.
                X = decomp_matrix.solve_left(L_B.basis_matrix())
                # Now row of X gives each element of L_B as a linear
                # combination of the rows of decomp_matrix.  We
                # project onto L_F by taking the right-most part of
                # this matrix.
                n = X.ncols()
                proj = X.matrix_from_columns(range(n - L_F.rank(), n))
                # Now proj is the matrix of projection that goes from
                # L_B to L_F, wrt the basis of those spaces.
                section = proj**(-1)

                # Now section maps L_F to L_B (tensor QQ).  Now we
                # just take each factor of F, which corresponds to a
                # submodule of L_F, and map it over to L_B tensor QQ
                # and saturate.
                D = []
                groups = self.groups()
                K = self.base_field()
                for A in Z_F:
                    L_A = A.lattice()
                    M = L_F.coordinate_module(L_A).basis_matrix() * section
                    M, _ = M._clear_denom()
                    M = M.saturation()
                    M = M * L_B.basis_matrix()
                    lattice = M.row_module(ZZ)
                    the_factor = ModularAbelianVariety(groups, lattice, K, is_simple=True, newform_level=A.newform_level(),
                                                       isogeny_number=A.isogeny_number(), number=A.degen_t())
                    D.append(the_factor)

        ################

        if isinstance(D, Sequence_generic):
            S = D
        else:
            D.sort()
            S = Sequence(D, immutable=True, cr=True, universe=self.category())
        self.__decomposition[(simple, bound)] = S
        return S

    def _classify_ambient_factors(self, simple=True, bound=None):
        r"""
        This function implements the following algorithm, which produces
        data useful in finding a decomposition or complement of self.


        #. Suppose `A_1 + \cdots + A_n` is a simple decomposition
           of the ambient space.

        #. For each `i`, let
           `B_i = A_1 + \cdots + A_i`.

        #. For each `i`, compute the intersection `C_i` of
           `B_i` and self.

        #. For each `i`, if the dimension of `C_i` is
           bigger than `C_{i-1}` put `i` in the "in" list;
           otherwise put `i` in the "out" list.


        Then one can show that self is isogenous to the sum of the
        `A_i` with `i` in the "in" list. Moreover, the sum
        of the `A_j` with `i` in the "out" list is a
        complement of self in the ambient space.

        INPUT:


        -  ``simple`` - bool (default: True)

        -  ``bound`` - integer (default: None); if given,
           passed onto decomposition function


        OUTPUT: IN list OUT list simple (or power of simple) factors

        EXAMPLES::

            sage: d1 = J0(11).degeneracy_map(33, 1); d1
            Degeneracy map from Abelian variety J0(11) of dimension 1 to Abelian variety J0(33) of dimension 3 defined by [1]
            sage: d2 = J0(11).degeneracy_map(33, 3); d2
            Degeneracy map from Abelian variety J0(11) of dimension 1 to Abelian variety J0(33) of dimension 3 defined by [3]
            sage: A = (d1 + d2).image(); A
            Abelian subvariety of dimension 1 of J0(33)
            sage: A._classify_ambient_factors()
            ([1], [0, 2], [
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33),
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33),
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
            ])
        """
        # Decompose an arbitrary abelian variety
        amb = self.ambient_variety()
        S = self.vector_space()
        X = amb.decomposition(simple=simple, bound=bound)
        IN = []
        OUT = []
        V = 0
        last_dimension = 0
        for j in range(len(X)):
            V += X[j].vector_space()
            d = S.intersection(V).dimension()
            if d > last_dimension:
                IN.append(j)
                last_dimension = d
            else:
                OUT.append(j)
        return IN, OUT, X

    def _isogeny_to_product_of_simples(self):
        r"""
        Given an abelian variety `A`, return an isogeny
        `\phi: A \rightarrow B_1 \times \cdots \times B_n`, where
        each `B_i` is simple. Note that this isogeny is not
        unique.

        EXAMPLES::

            sage: J = J0(37) ; J.decomposition()
            [
            Simple abelian subvariety 37a(1,37) of dimension 1 of J0(37),
            Simple abelian subvariety 37b(1,37) of dimension 1 of J0(37)
            ]
            sage: phi = J._isogeny_to_product_of_simples() ; phi
            Abelian variety morphism:
              From: Abelian variety J0(37) of dimension 2
              To:   Abelian subvariety of dimension 2 of J0(37) x J0(37)
            sage: J[0].intersection(J[1]) == phi.kernel()
            True

        ::

            sage: J = J0(22) * J0(37)
            sage: J._isogeny_to_product_of_simples()
            Abelian variety morphism:
              From: Abelian variety J0(22) x J0(37) of dimension 4
              To:   Abelian subvariety of dimension 4 of J0(11) x J0(11) x J0(37) x J0(37)
        """
        try:
            return self._simple_product_isogeny
        except AttributeError:
            pass

        D = self.decomposition()
        dest = prod([d._isogeny_to_newform_abelian_variety().image() for d in D])
        A = self.ambient_variety()
        dim = sum([d.dimension() for d in D])

        proj_ls = [A.projection(factor) for factor in D]

        mat = matrix(ZZ, 2 * self.dimension(), 2 * dim)
        ind = 0

        for i in range(len(D)):
            factor = D[i]
            proj = proj_ls[i]
            mat.set_block(0, ind, proj.restrict_domain(self).matrix())
            ind += 2 * factor.dimension()

        H = self.Hom(dest)
        self._simple_product_isogeny = H(Morphism(H, mat))
        return self._simple_product_isogeny

    def _isogeny_to_product_of_powers(self):
        r"""
        Given an abelian variety `A`, return an isogeny
        `\phi: A \rightarrow B_1 \times \cdots \times B_n`, where
        each `B_i` is a power of a simple abelian variety. These
        factors will be exactly those returned by
        ``self.decomposition(simple=False)``. Note that this isogeny is not
        unique.

        EXAMPLES::

            sage: J = J0(33) ; D = J.decomposition(simple=False) ; len(D)
            2
            sage: phi = J._isogeny_to_product_of_powers() ; phi
            Abelian variety morphism:
              From: Abelian variety J0(33) of dimension 3
              To:   Abelian subvariety of dimension 3 of J0(33) x J0(33)

        ::

            sage: J = J0(22) * J0(37)
            sage: J._isogeny_to_product_of_powers()
            Abelian variety morphism:
              From: Abelian variety J0(22) x J0(37) of dimension 4
              To:   Abelian subvariety of dimension 4 of J0(22) x J0(37) x J0(22) x J0(37) x J0(22) x J0(37)
        """
        try:
            return self._simple_power_product_isogeny
        except AttributeError:
            pass

        D = self.decomposition(simple=False)
        A = self.ambient_variety()
        proj_ls = [A.projection(factor) for factor in D]
        dest = prod([phi.image() for phi in proj_ls])
        dim = sum([d.dimension() for d in D])

        mat = matrix(ZZ, 2 * self.dimension(), 2 * dim)
        ind = 0

        for i in range(len(D)):
            factor = D[i]
            proj = proj_ls[i]
            mat.set_block(0, ind, proj.restrict_domain(self).matrix())
            ind += 2 * factor.dimension()

        H = self.Hom(dest)
        self._simple_power_product_isogeny = H(Morphism(H, mat))
        return self._simple_power_product_isogeny

    def complement(self, A=None):
        """
        Return a complement of this abelian variety.

        INPUT:

        -  ``A`` - (default: None); if given, A must be an
           abelian variety that contains self, in which case the complement of
           self is taken inside A. Otherwise the complement is taken in the
           ambient product Jacobian.


        OUTPUT: abelian variety

        EXAMPLES::

            sage: a,b,c = J0(33)
            sage: (a+b).complement()
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
            sage: (a+b).complement() == c
            True
            sage: a.complement(a+b)
            Abelian subvariety of dimension 1 of J0(33)
        """
        try:
            C = self.__complement
        except AttributeError:
            pass
        if self.dimension() == 0:
            if A is None:
                C = self.ambient_variety()
            else:
                C = A
        elif A is not None and self.dimension() == A.dimension():
            if not self.is_subvariety(A):
                raise ValueError("self must be a subvariety of A")
            C = self.zero_subvariety()
        else:
            _, factors, X = self._classify_ambient_factors()
            D = [X[i] for i in factors]
            C = sum(D)
            if C:
                self.__complement = C
                if A is not None:
                    C = C.intersection(A)[1]
            else:
                C = self.zero_subvariety()
        return C

    def dual(self):
        r"""
        Return the dual of this abelian variety.

        OUTPUT:

        - dual abelian variety
        - morphism from self to dual
        - covering morphism from J to dual

        .. warning::

           This is currently only implemented when self is an abelian
           subvariety of the ambient Jacobian product, and the
           complement of self in the ambient product Jacobian share no
           common factors. A more general implementation will require
           implementing computation of the intersection pairing on
           integral homology and the resulting Weil pairing on
           torsion.

        EXAMPLES: We compute the dual of the elliptic curve newform abelian
        variety of level `33`, and find the kernel of the modular
        map, which has structure `(\ZZ/3)^2`.

        ::

            sage: A,B,C = J0(33)
            sage: C
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
            sage: Cd, f, pi = C.dual()
            sage: f.matrix()
            [3 0]
            [0 3]
            sage: f.kernel()[0]
            Finite subgroup with invariants [3, 3] over QQ of Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)

        By a theorem the modular degree must thus be `3`::

            sage: E = EllipticCurve('33a')
            sage: E.modular_degree()
            3

        Next we compute the dual of a `2`-dimensional new simple
        abelian subvariety of `J_0(43)`.

        ::

            sage: A = AbelianVariety('43b'); A
            Newform abelian subvariety 43b of dimension 2 of J0(43)
            sage: Ad, f, pi = A.dual()

        The kernel shows that the modular degree is `2`::

            sage: f.kernel()[0]
            Finite subgroup with invariants [2, 2] over QQ of Newform abelian subvariety 43b of dimension 2 of J0(43)

        Unfortunately, the dual is not implemented in general::

            sage: A = J0(22)[0]; A
            Simple abelian subvariety 11a(1,22) of dimension 1 of J0(22)
            sage: A.dual()
            Traceback (most recent call last):
            ...
            NotImplementedError: dual not implemented unless complement shares no simple factors with self.
        """
        try:
            return self.__dual
        except AttributeError:
            if not self.is_subvariety_of_ambient_jacobian():
                raise NotImplementedError("dual not implemented unless abelian variety is a subvariety of the ambient Jacobian product")
            if not self._complement_shares_no_factors_with_same_label():
                raise NotImplementedError("dual not implemented unless complement shares no simple factors with self.")
            C = self.complement()
            Q, phi = self.ambient_variety().quotient(C)
            psi = self.ambient_morphism()
            self.__dual = Q, phi * psi, phi
            return self.__dual

    def _factors_with_same_label(self, other):
        """
        Given two modular abelian varieties self and other, this function
        returns a list of simple abelian subvarieties appearing in the
        decomposition of self that have the same newform labels. Each
        simple factor with a given newform label appears at most one.

        INPUT:


        -  ``other`` - abelian variety


        OUTPUT: list of simple abelian varieties

        EXAMPLES::

            sage: D = J0(33).decomposition(); D
            [
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33),
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33),
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
            ]
            sage: D[0]._factors_with_same_label(D[1])
            [Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)]
            sage: D[0]._factors_with_same_label(D[2])
            []
            sage: (D[0]+D[1])._factors_with_same_label(D[1] + D[2])
            [Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)]

        This illustrates that the multiplicities in the returned list are
        1::

            sage: (D[0]+D[1])._factors_with_same_label(J0(33))
            [Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)]

        This illustrates that the ambient product Jacobians do not have to
        be the same::

            sage: (D[0]+D[1])._factors_with_same_label(J0(22))
            [Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)]

        This illustrates that the actual factor labels are relevant, not
        just the isogeny class.

        ::

            sage: (D[0]+D[1])._factors_with_same_label(J1(11))
            []
            sage: J1(11)[0].newform_label()
            '11aG1'
        """
        if not isinstance(other, ModularAbelianVariety_abstract):
            raise TypeError("other must be an abelian variety")
        D = self.decomposition()
        C = set([A.newform_label() for A in other.decomposition()])
        Z = []
        for X in D:
            lbl = X.newform_label()
            if lbl in C:
                Z.append(X)
                C.remove(lbl)
        Z.sort()
        return Z

    def _complement_shares_no_factors_with_same_label(self):
        """
        Return True if no simple factor of self has the same newform_label
        as any factor in a Poincaré complement of self in the ambient
        product Jacobian.

        EXAMPLES: `J_0(37)` is made up of two non-isogenous
        elliptic curves::

            sage: J0(37)[0]._complement_shares_no_factors_with_same_label()
            True

        `J_0(33)` decomposes as a product of two isogenous
        elliptic curves with a third nonisogenous curve::

            sage: D = J0(33).decomposition(); D
            [
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33),
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33),
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
            ]
            sage: D[0]._complement_shares_no_factors_with_same_label()
            False
            sage: (D[0]+D[1])._complement_shares_no_factors_with_same_label()
            True
            sage: D[2]._complement_shares_no_factors_with_same_label()
            True

        This example illustrates the relevance of the ambient product
        Jacobian.

        ::

            sage: D = (J0(11) * J0(11)).decomposition(); D
            [
            Simple abelian subvariety 11a(1,11) of dimension 1 of J0(11) x J0(11),
            Simple abelian subvariety 11a(1,11) of dimension 1 of J0(11) x J0(11)
            ]
            sage: D[0]._complement_shares_no_factors_with_same_label()
            False

        This example illustrates that it is the newform label, not the
        isogeny, class that matters::

            sage: D = (J0(11)*J1(11)).decomposition(); D
            [
            Simple abelian subvariety 11aG1(1,11) of dimension 1 of J0(11) x J1(11),
            Simple abelian subvariety 11a(1,11) of dimension 1 of J0(11) x J1(11)
            ]
            sage: D[0]._complement_shares_no_factors_with_same_label()
            True
            sage: D[0].newform_label()
            '11aG1'
            sage: D[1].newform_label()
            '11a'
        """
        try:
            return self.__complement_shares
        except AttributeError:
            t = len(self._factors_with_same_label(self.complement())) == 0
            self.__complement_shares = t
            return t

    def __getitem__(self, i):
        """
        Return the `i^{th}` decomposition factor of self
        or return the slice `i` of decompositions of self.

        EXAMPLES::

            sage: J = J0(389)
            sage: J.decomposition()
            [
            Simple abelian subvariety 389a(1,389) of dimension 1 of J0(389),
            Simple abelian subvariety 389b(1,389) of dimension 2 of J0(389),
            Simple abelian subvariety 389c(1,389) of dimension 3 of J0(389),
            Simple abelian subvariety 389d(1,389) of dimension 6 of J0(389),
            Simple abelian subvariety 389e(1,389) of dimension 20 of J0(389)
            ]
            sage: J[2]
            Simple abelian subvariety 389c(1,389) of dimension 3 of J0(389)
            sage: J[-1]
            Simple abelian subvariety 389e(1,389) of dimension 20 of J0(389)
            sage: J = J0(125); J.decomposition()
            [
            Simple abelian subvariety 125a(1,125) of dimension 2 of J0(125),
            Simple abelian subvariety 125b(1,125) of dimension 2 of J0(125),
            Simple abelian subvariety 125c(1,125) of dimension 4 of J0(125)
            ]
            sage: J[:2]
            [
            Simple abelian subvariety 125a(1,125) of dimension 2 of J0(125),
            Simple abelian subvariety 125b(1,125) of dimension 2 of J0(125)
            ]
        """
        return self.decomposition()[i]


class ModularAbelianVariety(ModularAbelianVariety_abstract):
    def __init__(self, groups, lattice=None, base_field=QQ, is_simple=None, newform_level=None,
                 isogeny_number=None, number=None, check=True):
        r"""
        Create a modular abelian variety with given level and base field.

        INPUT:


        -  ``groups`` - a tuple of congruence subgroups

        -  ``lattice`` - (default: `\ZZ^n`) a
           full lattice in `\ZZ^n`, where `n` is the
           sum of the dimensions of the spaces of cuspidal modular symbols
           corresponding to each `\Gamma \in` groups

        -  ``base_field`` - a field (default:
           `\QQ`)


        EXAMPLES::

            sage: J0(23)
            Abelian variety J0(23) of dimension 2
        """
        ModularAbelianVariety_abstract.__init__(self, groups, base_field, is_simple=is_simple, newform_level=newform_level,
                                                isogeny_number=isogeny_number, number=number, check=check)
        if lattice is None:
            lattice = ZZ**(2 * self._ambient_dimension())
        if check:
            n = self._ambient_dimension()
            if not is_FreeModule(lattice):
                raise TypeError("lattice must be a free module")
            if lattice.base_ring() != ZZ:
                raise TypeError("lattice must be over ZZ")
            if lattice.degree() != 2 * n:
                raise ValueError("lattice must have degree 2*n (=%s)" % (2 * n))
            if not lattice.saturation().is_submodule(lattice):  # potentially expensive
                raise ValueError("lattice must be full")
        self.__lattice = lattice

    def lattice(self):
        """
        Return the lattice that defines this abelian variety.

        OUTPUT:


        -  ``lattice`` - a lattice embedded in the rational
           homology of the ambient product Jacobian


        EXAMPLES::

            sage: A = (J0(11) * J0(37))[1]; A
            Simple abelian subvariety 37a(1,37) of dimension 1 of J0(11) x J0(37)
            sage: type(A)
            <class 'sage.modular.abvar.abvar.ModularAbelianVariety_with_category'>
            sage: A.lattice()
            Free module of degree 6 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 0  0  1 -1  1  0]
            [ 0  0  0  0  2 -1]
        """
        return self.__lattice


class ModularAbelianVariety_modsym_abstract(ModularAbelianVariety_abstract):
    # Anything that derives from this class must define the
    # modular_symbols method, which returns a cuspidal modular symbols
    # space over QQ.  It can have any sign.
    def _modular_symbols(self):
        """
        Return the space of modular symbols corresponding to this modular
        symbols abelian variety.

        EXAMPLES: This function is in the abstract base class, so it raises
        a NotImplementedError::

            sage: M = ModularSymbols(37).cuspidal_submodule()
            sage: A = M.abelian_variety(); A
            Abelian variety J0(37) of dimension 2
            sage: sage.modular.abvar.abvar.ModularAbelianVariety_modsym_abstract._modular_symbols(A)
            Traceback (most recent call last):
            ...
            NotImplementedError: bug -- must define this

        Of course this function isn't called in practice, so this works::

            sage: A._modular_symbols()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
        """
        raise NotImplementedError("bug -- must define this")

    def __add__(self, other):
        """
        Add two modular abelian variety factors.

        EXAMPLES::

            sage: A = J0(42); D = A.decomposition(); D
            [
            Simple abelian subvariety 14a(1,42) of dimension 1 of J0(42),
            Simple abelian subvariety 14a(3,42) of dimension 1 of J0(42),
            Simple abelian subvariety 21a(1,42) of dimension 1 of J0(42),
            Simple abelian subvariety 21a(2,42) of dimension 1 of J0(42),
            Simple abelian subvariety 42a(1,42) of dimension 1 of J0(42)
            ]
            sage: D[0] + D[1]
            Abelian subvariety of dimension 2 of J0(42)
            sage: D[1].is_subvariety(D[0] + D[1])
            True
            sage: D[0] + D[1] + D[2]
            Abelian subvariety of dimension 3 of J0(42)
            sage: D[0] + D[0]
            Abelian subvariety of dimension 1 of J0(42)
            sage: D[0] + D[0] == D[0]
            True
            sage: sum(D, D[0]) == A
            True
        """
        if not is_ModularAbelianVariety(other):
            if other == 0:
                return self
            raise TypeError("sum not defined")
        if not isinstance(other, ModularAbelianVariety_modsym_abstract):
            return ModularAbelianVariety_abstract.__add__(self, other)
        if self.groups() != other.groups():
            raise TypeError("sum not defined since ambient spaces different")
        M = self.modular_symbols() + other.modular_symbols()
        return ModularAbelianVariety_modsym(M)

    def groups(self):
        """
        Return the tuple of groups associated to the modular symbols
        abelian variety. This is always a 1-tuple.

        OUTPUT: tuple

        EXAMPLES::

            sage: A = ModularSymbols(33).cuspidal_submodule().abelian_variety(); A
            Abelian variety J0(33) of dimension 3
            sage: A.groups()
            (Congruence Subgroup Gamma0(33),)
            sage: type(A)
            <class 'sage.modular.abvar.abvar.ModularAbelianVariety_modsym_with_category'>
        """
        return (self._modular_symbols().group(), )

    def lattice(self):
        r"""
        Return the lattice defining this modular abelian variety.

        OUTPUT:

        A free `\ZZ`-module embedded in an ambient `\QQ`-vector space.

        EXAMPLES::

            sage: A = ModularSymbols(33).cuspidal_submodule()[0].abelian_variety(); A
            Abelian subvariety of dimension 1 of J0(33)
            sage: A.lattice()
            Free module of degree 6 and rank 2 over Integer Ring
            User basis matrix:
            [ 1  0  0 -1  0  0]
            [ 0  0  1  0  1 -1]
            sage: type(A)
            <class 'sage.modular.abvar.abvar.ModularAbelianVariety_modsym_with_category'>
        """
        try:
            return self.__lattice
        except AttributeError:
            M = self.modular_symbols()
            S = M.ambient_module().cuspidal_submodule()
            if M.dimension() == S.dimension():
                L = ZZ**M.dimension()
            else:
                K0 = M.integral_structure()
                K1 = S.integral_structure()
                L = K1.coordinate_module(K0)
            self.__lattice = L
            return self.__lattice

    def _set_lattice(self, lattice):
        """
        Set the lattice of this modular symbols abelian variety.

        .. warning::

           This is only for internal use. Do not use this unless you
           really really know what you're doing. That's why there is
           an underscore in this method name.

        INPUT:


        -  ``lattice`` - a lattice


        EXAMPLES: We do something evil - there's no type checking since
        this function is for internal use only::

            sage: A = ModularSymbols(33).cuspidal_submodule().abelian_variety()
            sage: A._set_lattice(5)
            sage: A.lattice()
            5
        """
        self.__lattice = lattice

    def modular_symbols(self, sign=0):
        """
        Return space of modular symbols (with given sign) associated to
        this modular abelian variety, if it can be found by cutting down
        using Hecke operators. Otherwise raise a RuntimeError exception.

        EXAMPLES::

            sage: A = J0(37)
            sage: A.modular_symbols()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
            sage: A.modular_symbols(1)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(37) of weight 2 with sign 1 over Rational Field

        More examples::

            sage: J0(11).modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: J0(11).modular_symbols(sign=1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: J0(11).modular_symbols(sign=0)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: J0(11).modular_symbols(sign=-1)
            Modular Symbols space of dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Rational Field

        Even more examples::

            sage: A = J0(33)[1]; A
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33)
            sage: A.modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field

        It is not always possible to determine the sign subspaces::

            sage: A.modular_symbols(1)
            Traceback (most recent call last):
            ...
            RuntimeError: unable to determine sign (=1) space of modular symbols

        ::

            sage: A.modular_symbols(-1)
            Traceback (most recent call last):
            ...
            RuntimeError: unable to determine sign (=-1) space of modular symbols
        """
        M = self._modular_symbols().modular_symbols_of_sign(sign)
        if (sign != 0 and M.dimension() != self.dimension()) or (sign == 0 and M.dimension() != 2 * self.dimension()):
            raise RuntimeError("unable to determine sign (=%s) space of modular symbols" % sign)
        return M

    def _compute_hecke_polynomial(self, n, var='x'):
        r"""
        Return the characteristic polynomial of the `n^{th}` Hecke
        operator on ``self``.

        .. note::

            If self has dimension d, then this is a polynomial of
            degree d. It is not of degree 2\*d, so it is the square
            root of the characteristic polynomial of the Hecke operator
            on integral or rational homology (which has degree 2\*d).

        EXAMPLES::

            sage: J0(11).hecke_polynomial(2)
            x + 2
            sage: J0(23)._compute_hecke_polynomial(2)
            x^2 + x - 1
            sage: J1(13).hecke_polynomial(2)
            x^2 + 3*x + 3
            sage: factor(J0(43).hecke_polynomial(2))
            (x + 2) * (x^2 - 2)

        The Hecke polynomial is the square root of the characteristic
        polynomial::

            sage: factor(J0(43).hecke_operator(2).charpoly())
            (x + 2) * (x^2 - 2)
        """
        b, r = self.modular_symbols().hecke_polynomial(n, var).is_square(True)
        return r

    def _integral_hecke_matrix(self, n, sign=0):
        """
        Return the action of the Hecke operator `T_n` on the
        integral homology of self.

        INPUT:


        -  ``n`` - a positive integer

        -  ``sign`` - 0, +1, or -1; if 1 or -1 act on the +1 or
           -1 quotient of the integral homology.


        EXAMPLES::

            sage: J1(13)._integral_hecke_matrix(2)     # slightly random choice of basis
            [-2  0 -1  1]
            [ 1 -1  0 -1]
            [ 1  1 -2  0]
            [ 0  1 -1 -1]
            sage: J1(13)._integral_hecke_matrix(2,sign=1)  # slightly random choice of basis
            [-1  1]
            [-1 -2]
            sage: J1(13)._integral_hecke_matrix(2,sign=-1)  # slightly random choice of basis
            [-2 -1]
            [ 1 -1]
        """
        return self.modular_symbols(sign).integral_hecke_matrix(n)

    def _rational_hecke_matrix(self, n, sign=0):
        """
        Return the action of the Hecke operator `T_n` on the
        rational homology of self.

        INPUT:


        -  ``n`` - a positive integer

        -  ``sign`` - 0, +1, or -1; if 1 or -1 act on the +1 or
           -1 quotient of the rational homology.


        EXAMPLES::

            sage: J1(13)._rational_hecke_matrix(2)    # slightly random choice of basis
            [-2  0 -1  1]
            [ 1 -1  0 -1]
            [ 1  1 -2  0]
            [ 0  1 -1 -1]
            sage: J0(43)._rational_hecke_matrix(2,sign=1)  # slightly random choice of basis
            [-2  0  1]
            [-1 -2  2]
            [-2  0  2]
        """
        return self._integral_hecke_matrix(n, sign=sign).change_ring(QQ)

    def group(self):
        """
        Return the congruence subgroup associated that this modular abelian
        variety is associated to.

        EXAMPLES::

            sage: J0(13).group()
            Congruence Subgroup Gamma0(13)
            sage: J1(997).group()
            Congruence Subgroup Gamma1(997)
            sage: JH(37,[3]).group()
            Congruence Subgroup Gamma_H(37) with H generated by [3]
            sage: J0(37)[1].groups()
            (Congruence Subgroup Gamma0(37),)
        """
        return self.modular_symbols().group()

    def is_subvariety(self, other):
        """
        Return True if self is a subvariety of other.

        EXAMPLES::

            sage: J = J0(37); J
            Abelian variety J0(37) of dimension 2
            sage: A = J[0]; A
            Simple abelian subvariety 37a(1,37) of dimension 1 of J0(37)
            sage: A.is_subvariety(J)
            True
            sage: A.is_subvariety(J0(11))
            False

        There may be a way to map `A` into `J_0(74)`, but
        `A` is not equipped with any special structure of an
        embedding.

        ::

            sage: A.is_subvariety(J0(74))
            False

        Some ambient examples::

            sage: J = J0(37)
            sage: J.is_subvariety(J)
            True
            sage: J.is_subvariety(25)
            False

        More examples::

            sage: A = J0(42); D = A.decomposition(); D
            [
            Simple abelian subvariety 14a(1,42) of dimension 1 of J0(42),
            Simple abelian subvariety 14a(3,42) of dimension 1 of J0(42),
            Simple abelian subvariety 21a(1,42) of dimension 1 of J0(42),
            Simple abelian subvariety 21a(2,42) of dimension 1 of J0(42),
            Simple abelian subvariety 42a(1,42) of dimension 1 of J0(42)
            ]
            sage: D[0].is_subvariety(A)
            True
            sage: D[1].is_subvariety(D[0] + D[1])
            True
            sage: D[2].is_subvariety(D[0] + D[1])
            False
        """
        if not is_ModularAbelianVariety(other):
            return False
        if not isinstance(other, ModularAbelianVariety_modsym_abstract):
            return ModularAbelianVariety_abstract.is_subvariety(self, other)
        return self.modular_symbols().is_submodule(other.modular_symbols())

    def is_ambient(self):
        """
        Return True if this abelian variety attached to a modular symbols
        space is attached to the cuspidal subspace of the ambient
        modular symbols space.

        OUTPUT: bool

        EXAMPLES::

            sage: A = ModularSymbols(43).cuspidal_subspace().abelian_variety(); A
            Abelian variety J0(43) of dimension 3
            sage: A.is_ambient()
            True
            sage: type(A)
            <class 'sage.modular.abvar.abvar.ModularAbelianVariety_modsym_with_category'>
            sage: A = ModularSymbols(43).cuspidal_subspace()[1].abelian_variety(); A
            Abelian subvariety of dimension 2 of J0(43)
            sage: A.is_ambient()
            False
        """
        return self.degree() == self.dimension()

    def dimension(self):
        """
        Return the dimension of this modular abelian variety.

        EXAMPLES::

            sage: J0(37)[0].dimension()
            1
            sage: J0(43)[1].dimension()
            2
            sage: J1(17)[1].dimension()
            4
        """
        try:
            return self._dimension
        except AttributeError:
            M = self._modular_symbols()
            if M.sign() == 0:
                d = M.dimension() // 2
            else:
                d = M.dimension()
            self._dimension = d
            return d

    def new_subvariety(self, p=None):
        """
        Return the new or `p`-new subvariety of self.

        INPUT:


        -  ``self`` - a modular abelian variety

        -  ``p`` - prime number or None (default); if p is a
           prime, return the p-new subvariety. Otherwise return the full new
           subvariety.


        EXAMPLES::

            sage: J0(34).new_subvariety()
            Abelian subvariety of dimension 1 of J0(34)
            sage: J0(100).new_subvariety()
            Abelian subvariety of dimension 1 of J0(100)
            sage: J1(13).new_subvariety()
            Abelian variety J1(13) of dimension 2
        """
        try:
            return self.__new_subvariety[p]
        except AttributeError:
            self.__new_subvariety = {}
        except KeyError:
            pass
        A = self.modular_symbols()
        N = A.new_submodule(p=p)
        B = ModularAbelianVariety_modsym(N)
        self.__new_subvariety[p] = B
        return B

    def old_subvariety(self, p=None):
        """
        Return the old or `p`-old abelian variety of self.

        INPUT:


        -  ``self`` - a modular abelian variety

        -  ``p`` - prime number or None (default); if p is a
           prime, return the p-old subvariety. Otherwise return the full old
           subvariety.


        EXAMPLES::

            sage: J0(33).old_subvariety()
            Abelian subvariety of dimension 2 of J0(33)
            sage: J0(100).old_subvariety()
            Abelian subvariety of dimension 6 of J0(100)
            sage: J1(13).old_subvariety()
            Abelian subvariety of dimension 0 of J1(13)
        """
        try:
            return self.__old_subvariety[p]
        except AttributeError:
            self.__old_subvariety = {}
        except KeyError:
            pass
        A = self.modular_symbols()
        N = A.old_submodule(p=p)
        B = ModularAbelianVariety_modsym(N)
        self.__old_subvariety[p] = B
        return B

    def decomposition(self, simple=True, bound=None):
        r"""
        Decompose this modular abelian variety as a product of abelian
        subvarieties, up to isogeny.

        INPUT: simple- bool (default: True) if True, all factors are
        simple. If False, each factor returned is isogenous to a power of a
        simple and the simples in each factor are distinct.


        -  ``bound`` - int (default: None) if given, only use
           Hecke operators up to this bound when decomposing. This can give
           wrong answers, so use with caution!


        EXAMPLES::

            sage: J = J0(33)
            sage: J.decomposition()
            [
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33),
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(33),
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(33)
            ]
            sage: J1(17).decomposition()
            [
            Simple abelian subvariety 17aG1(1,17) of dimension 1 of J1(17),
            Simple abelian subvariety 17bG1(1,17) of dimension 4 of J1(17)
            ]
        """
        try:
            return self.__decomposition[(simple, bound)]
        except KeyError:
            pass
        except AttributeError:
            self.__decomposition = {}
        if not self.is_ambient():
            S = ModularAbelianVariety_abstract.decomposition(self, simple=simple, bound=bound)
        else:
            A = self.modular_symbols()
            amb = A.ambient_module()
            G = amb.group()
            S = amb.cuspidal_submodule().integral_structure()
            if simple:
                M = A.level()
                D = []
                for N in reversed(divisors(M)):
                    if N > 1:
                        isogeny_number = 0
                        A = amb.modular_symbols_of_level(N).cuspidal_subspace().new_subspace()
                        if bound is None:
                            X = factor_new_space(A)
                        else:
                            X = A.decomposition(bound=bound)
                        for B in X:
                            for t in divisors(M // N):
                                D.append(ModularAbelianVariety_modsym(B.degeneracy_map(M, t).image(),
                                                                      is_simple=True, newform_level=(N, G),
                                                                      isogeny_number=isogeny_number,
                                                                      number=(t, M)))
                            isogeny_number += 1
            elif A == amb.cuspidal_submodule():
                D = [ModularAbelianVariety_modsym(B)
                     for B in A.decomposition(bound=bound)]
            else:
                D = ModularAbelianVariety_abstract.decomposition(self, simple=simple, bound=bound)
            D.sort()
            S = Sequence(D, immutable=True, cr=True, universe=self.category())
        self.__decomposition[(simple, bound)] = S
        return S


class ModularAbelianVariety_modsym(ModularAbelianVariety_modsym_abstract):

    def __init__(self, modsym, lattice=None, newform_level=None,
                 is_simple=None, isogeny_number=None, number=None, check=True):
        """
        Modular abelian variety that corresponds to a Hecke stable space of
        cuspidal modular symbols.

        EXAMPLES: We create a modular abelian variety attached to a space
        of modular symbols.

        ::

            sage: M = ModularSymbols(23).cuspidal_submodule()
            sage: A = M.abelian_variety(); A
            Abelian variety J0(23) of dimension 2
        """
        if check:
            if not isinstance(modsym, ModularSymbolsSpace):
                raise TypeError("modsym must be a modular symbols space")
            if modsym.sign() != 0:
                raise TypeError("modular symbols space must have sign 0")
            if not modsym.is_cuspidal():
                raise ValueError("modsym must be cuspidal")

        ModularAbelianVariety_abstract.__init__(self, (modsym.group(), ), modsym.base_ring(),
                             newform_level=newform_level, is_simple=is_simple,
                             isogeny_number=isogeny_number, number=number, check=check)
        if lattice is not None:
            self._set_lattice(lattice)
        self.__modsym = modsym

    def _modular_symbols(self):
        """
        Return the modular symbols space that defines this modular abelian
        variety.

        OUTPUT: space of modular symbols

        EXAMPLES::

            sage: M = ModularSymbols(37).cuspidal_submodule()
            sage: A = M.abelian_variety(); A
            Abelian variety J0(37) of dimension 2
            sage: A._modular_symbols()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
        """
        return self.__modsym

    def component_group_order(self, p):
        """
        Return the order of the component group of the special fiber
        at p of the Neron model of self.

        NOTE: For bad primes, this is only implemented when the group
        if Gamma0 and p exactly divides the level.

        NOTE: the input abelian variety must be simple

        ALGORITHM: See "Component Groups of Quotients of J0(N)" by Kohel and Stein.  That
        paper is about optimal quotients; however, section 4.1 of Conrad-Stein "Component
        Groups of Purely Toric Quotients", one sees that the component group of an optimal
        quotient is the same as the component group of its dual (which is the subvariety).

        INPUT:

        - p -- a prime number

        OUTPUT:

        - Integer

        EXAMPLES::

            sage: A = J0(37)[1]
            sage: A.component_group_order(37)
            3
            sage: A = J0(43)[1]
            sage: A.component_group_order(37)
            1
            sage: A.component_group_order(43)
            7
            sage: A = J0(23)[0]
            sage: A.component_group_order(23)
            11
        """
        if not self.is_simple():
            raise ValueError("self must be simple")
        p = Integer(p)
        if not p.is_prime():
            raise ValueError("p must be a prime integer")
        try:
            return self.__component_group[p][0]
        except AttributeError:
            self.__component_group = {}
        except KeyError:
            pass
        # Easy special case -- a prime of good reduction
        if self.level() % p:
            one = Integer(1)
            self.__component_group[p] = (one, one, one)
            return one
        # Cases that we don't know how to handle yet.
        if not is_Gamma0(self.group()):
            raise NotImplementedError("computation of component group not implemented when group isn't Gamma0")
        if self.level() % (p * p) == 0:
            raise NotImplementedError("computation of component group not implemented when p^2 divides the level")

        # Now we're on Gamma0(p*M) with gcd(p,M) = 1.
        # 1. Compute factor of Brandt module space, and put integral structure on it.

        # TODO -- in case self.level() is prime, should use
        # supersingular module instead for massive speedup...  Of
        # course, then one can just use Emertons theorem that the
        # component group order equals the torsion order, and avoid
        # all of this!
        XI = self.brandt_module(p)
        Y = XI.ambient_module()
        n = Y.dimension()

        # X_ZZ is the submodule of degree 0 divisors
        M = ZZ**n
        deg_zero = []
        for k in range(1, n):
            v = vector(ZZ, n)
            v[0] = 1
            v[k] = -1
            deg_zero.append(v)
        X_ZZ = M.span(deg_zero, ZZ)
        XI_ZZ = XI.free_module().intersection(M)

        # 2. Compute the map alpha: X --> Hom(X[I],Z) over ZZ
        # todo -- this could be done more quickly with a clever matrix multiply
        B = [XI(vv) for vv in XI_ZZ.basis()]
        mat = []
        for v in M.basis():
            w = Y(v)
            mat.append([w.monodromy_pairing(b) for b in B])
        monodromy = matrix(ZZ, mat)
        alpha = X_ZZ.basis_matrix().change_ring(ZZ) * monodromy

        # 3. Compute invariants:
        #        * Phi_X = #coker(alpha)
        #        * m_X = #(alpha(X)/alpha(X[I]))
        alphaX = alpha.row_module()
        Phi_X_invariants = alphaX.basis_matrix().change_ring(ZZ).elementary_divisors()
        Phi_X = prod(Phi_X_invariants + [Integer(1)])

        W = alphaX.span([b * monodromy for b in XI_ZZ.basis()], ZZ)
        m_X = Integer(W.index_in(alphaX))

        # 4. Compute the modular degree
        moddeg = self.modular_degree()

        # 5. Obtain the component group order using Theorem 1 of [Kohel-Stein]
        Phi = Phi_X * moddeg / m_X

        # 6. Record the answer
        self.__component_group[p] = (Phi, Phi_X_invariants, m_X)
        return Phi

    def _invariants_of_image_of_component_group_of_J0(self, p):
        """
        Return the elementary invariants of the image of the component
        group of J0(N).  The API of this function is subject to
        change, which is why it starts with an underscore.

        INPUT:

        - p -- integer

        OUTPUT:

        - list -- of elementary invariants

        EXAMPLES::

            sage: A = J0(62).new_subvariety()[1]; A
            Simple abelian subvariety 62b(1,62) of dimension 2 of J0(62)
            sage: A._invariants_of_image_of_component_group_of_J0(2)
            [1, 6]
            sage: A.component_group_order(2)
            66
        """
        self.component_group_order(p)
        return list(self.__component_group[p][1])   # make a copy

    def tamagawa_number(self, p):
        """
        Return the Tamagawa number of this abelian variety at p.

        NOTE: For bad primes, this is only implemented when the group
        if Gamma0 and p exactly divides the level and Atkin-Lehner
        acts diagonally on this abelian variety (e.g., if this variety
        is new and simple).  See the self.component_group command for
        more information.

        NOTE: the input abelian variety must be simple

        In cases where this function doesn't work, consider using the
        self.tamagawa_number_bounds functions.

        INPUT:

        - p -- a prime number

        OUTPUT:

        - Integer

        EXAMPLES::

            sage: A = J0(37)[1]
            sage: A.tamagawa_number(37)
            3
            sage: A = J0(43)[1]
            sage: A.tamagawa_number(37)
            1
            sage: A.tamagawa_number(43)
            7
            sage: A = J0(23)[0]
            sage: A.tamagawa_number(23)
            11
        """
        try:
            return self.__tamagawa_number[p]
        except AttributeError:
            self.__tamagawa_number = {}
        except KeyError:
            pass
        if not self.is_simple():
            raise ValueError("self must be simple")
        try:
            self.component_group_order(p)
        except NotImplementedError:
            raise NotImplementedError("Tamagawa number can't be determined using known algorithms, so consider using the tamagawa_number_bounds function instead")
        div, mul, mul_primes = self.tamagawa_number_bounds(p)
        if div == mul:
            cp = div
        else:
            raise NotImplementedError("the Tamagawa number at %s is a power of 2, but the exact power can't be determined using known algorithms.  Consider using the tamagawa_number_bounds function instead." % p)
        self.__tamagawa_number[p] = cp
        return cp

    def tamagawa_number_bounds(self, p):
        """
        Return a divisor and multiple of the Tamagawa number of self at `p`.

        NOTE: the input abelian variety must be simple.

        INPUT:

        - `p` -- a prime number

        OUTPUT:

        - div -- integer; divisor of Tamagawa number at `p`
        - mul -- integer; multiple of Tamagawa number at `p`
        - mul_primes -- tuple; in case mul==0, a list of all
          primes that can possibly divide the Tamagawa number at `p`

        EXAMPLES::

            sage: A = J0(63).new_subvariety()[1]; A
            Simple abelian subvariety 63b(1,63) of dimension 2 of J0(63)
            sage: A.tamagawa_number_bounds(7)
            (3, 3, ())
            sage: A.tamagawa_number_bounds(3)
            (1, 0, (2, 3, 5))
        """
        try:
            return self.__tamagawa_number_bounds[p]
        except AttributeError:
            self.__tamagawa_number_bounds = {}
        except KeyError:
            pass
        if not self.is_simple():
            raise ValueError("self must be simple")
        N = self.level()
        div = 1
        mul = 0
        mul_primes = []
        if N % p:
            div = 1
            mul = 1
        elif N.valuation(p) == 1:
            M = self.modular_symbols(sign=1)
            if is_Gamma0(M.group()):
                g = self.component_group_order(p)
                W = M.atkin_lehner_operator(p).matrix()
                if W == -1:
                    # Frob acts trivially
                    div = g
                    mul = g
                elif W == 1:
                    # Frob acts by -1
                    n = g.valuation(2)
                    if n <= 1:
                        div = 2**n
                    else:
                        phi_X_invs = self._invariants_of_image_of_component_group_of_J0(p)
                        m = max(1, len([z for z in phi_X_invs if z % 2 == 0]))
                        div = 2**m
                    mul = 2**n
                else:
                    raise NotImplementedError("Atkin-Lehner at p must act as a scalar")
        else:
            mul_primes = sorted(set([p] + [q for q in prime_range(2, 2 * self.dimension() + 2)]))
        div = Integer(div)
        mul = Integer(mul)
        mul_primes = tuple(mul_primes)
        self.__tamagawa_number_bounds[p] = (div, mul, mul_primes)
        return (div, mul, mul_primes)

    def brandt_module(self, p):
        """
        Return the Brandt module at p that corresponds to self.  This
        is the factor of the vector space on the ideal class set in an
        order of level N in the quaternion algebra ramified at p and
        infinity.

        INPUT:

        - p -- prime that exactly divides the level

        OUTPUT:

        - Brandt module space that corresponds to self.

        EXAMPLES::

            sage: J0(43)[1].brandt_module(43)
            Subspace of dimension 2 of Brandt module of dimension 4 of level 43 of weight 2 over Rational Field
            sage: J0(43)[1].brandt_module(43).basis()
            ((1, 0, -1/2, -1/2), (0, 1, -1/2, -1/2))
            sage: J0(43)[0].brandt_module(43).basis()
            ((0, 0, 1, -1),)
            sage: J0(35)[0].brandt_module(5).basis()
            ((1, 0, -1, 0),)
            sage: J0(35)[0].brandt_module(7).basis()
            ((1, -1, 1, -1),)
        """
        try:
            return self.__brandt_module[p]
        except AttributeError:
            self.__brandt_module = {}
        except KeyError:
            pass
        p = Integer(p)
        if not is_Gamma0(self.group()):
            raise NotImplementedError("Brandt module only defined on Gamma0")
        if not p.is_prime():
            raise ValueError("p must be a prime integer")
        if self.level().valuation(p) != 1:
            raise ValueError("p must exactly divide the level")
        M = self.level() / p
        V = BrandtModule(p, M)
        # now cut out version of self in B
        S = self.modular_symbols(sign=1)
        B = S.hecke_bound()
        if self.dimension() <= 3:
            q = 2
            while V.dimension() > self.dimension() and q <= B:
                f = S.hecke_polynomial(q)
                V = f(V.hecke_operator(q)).kernel()
                q = next_prime(q)
            if V.dimension() > self.dimension():
                raise RuntimeError("unable to cut out Brandt module (got dimension %s instead of %s)" % (V.dimension(), self.dimension()))
        else:
            D = V.decomposition()
            D = [A for A in D if A.dimension() == self.dimension()]
            # now figure out which element of D is isomorphic to self.
            q = 2
            while len(D) > 1 and q <= B:
                f = S.hecke_polynomial(q)
                D = [A for A in D if A.hecke_polynomial(q) == f]
                q = next_prime(q)
            if len(D) != 1:
                raise RuntimeError("unable to locate Brandt module (got %s candidates instead of 1)" % (len(D)))
            V = D[0]
        self.__brandt_module[p] = V
        return V


def sqrt_poly(f):
    """
    Return the square root of the polynomial `f`.

    .. note::

        At some point something like this should be a member of the
        polynomial class. For now this is just used internally by some
        charpoly functions above.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: f = (x-1)*(x+2)*(x^2 + 1/3*x + 5)
        sage: f
        x^4 + 4/3*x^3 + 10/3*x^2 + 13/3*x - 10
        sage: sage.modular.abvar.abvar.sqrt_poly(f^2)
        x^4 + 4/3*x^3 + 10/3*x^2 + 13/3*x - 10
        sage: sage.modular.abvar.abvar.sqrt_poly(f)
        Traceback (most recent call last):
        ...
        ValueError: f must be a perfect square
        sage: sage.modular.abvar.abvar.sqrt_poly(2*f^2)
        Traceback (most recent call last):
        ...
        ValueError: f must be monic
    """
    if not f.is_monic():
        raise ValueError("f must be monic")
    try:
        return prod([g**Integer(e / Integer(2)) for g, e in f.factor()])
    except TypeError:
        raise ValueError("f must be a perfect square")


##############################################################################
# Useful for decomposing exactly the sort of modular symbols spaces that come up here.


def random_hecke_operator(M, t=None, p=2):
    """
    Return a random Hecke operator acting on `M`, got by adding
    to `t` a random multiple of `T_p`

    INPUT:


    -  ``M`` - modular symbols space

    -  ``t`` - None or a Hecke operator

    -  ``p`` - a prime


    OUTPUT: Hecke operator prime

    EXAMPLES::

        sage: M = ModularSymbols(11).cuspidal_subspace()
        sage: t, p = sage.modular.abvar.abvar.random_hecke_operator(M)
        sage: p
        3
        sage: t, p = sage.modular.abvar.abvar.random_hecke_operator(M, t, p)
        sage: p
        5
    """
    r = 0
    while r == 0:
        r = randrange(1, p // 2 + 1) * ZZ.random_element()
    t = (0 if t is None else t) + r * M.hecke_operator(p)
    return t, next_prime(p)


def factor_new_space(M):
    """
    Given a new space `M` of modular symbols, return the
    decomposition into simple of `M` under the Hecke
    operators.

    INPUT:


    -  ``M`` - modular symbols space


    OUTPUT: list of factors

    EXAMPLES::

        sage: M = ModularSymbols(37).cuspidal_subspace()
        sage: sage.modular.abvar.abvar.factor_new_space(M)
        [
        Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field,
        Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
        ]
    """
    t = None
    p = 2
    for i in range(200):
        t, p = random_hecke_operator(M, t, p)
        f = t.charpoly()
        cube_free = True
        for _, e in f.factor():
            if e > 2:
                cube_free = False
                break
        if cube_free:
            return t.decomposition()
        t, p = random_hecke_operator(M, t, p)
    raise RuntimeError("unable to factor new space -- this should not happen")


def factor_modsym_space_new_factors(M):
    """
    Return the factorizations of all the new subspaces of `M`.

    INPUT:

    - ``M`` -- ambient modular symbols space

    OUTPUT: list of decompositions corresponding to each new space.

    EXAMPLES::

        sage: M = ModularSymbols(33)
        sage: sage.modular.abvar.abvar.factor_modsym_space_new_factors(M)
        [[
        Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
        ],
         [
        Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
        ]]
    """
    eps = M.character()
    K = eps.conductor() if eps is not None else 1
    N = [M.modular_symbols_of_level(d).cuspidal_subspace().new_subspace()
         for d in M.level().divisors() if d % K == 0 and (d == 11 or d >= 13)]
    return [factor_new_space(A) for A in N]


def simple_factorization_of_modsym_space(M, simple=True):
    """
    Return the canonical factorization of `M` into (simple) subspaces.

    INPUT:

    - ``M`` -- ambient modular symbols space

    - ``simple`` -- boolean (default: ``True``); if set to ``False``,
      isogenous simple factors are grouped together

    OUTPUT: sequence

    EXAMPLES::

        sage: M = ModularSymbols(33)
        sage: sage.modular.abvar.abvar.simple_factorization_of_modsym_space(M)
        [
        (11, 0, 1, Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field),
        (11, 0, 3, Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field),
        (33, 0, 1, Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field)
        ]
        sage: sage.modular.abvar.abvar.simple_factorization_of_modsym_space(M, simple=False)
        [
        (11, 0, None, Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field),
        (33, 0, None, Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field)
        ]

    TESTS:

    Check that :trac:`21799` is fixed::

        sage: JH(28, [15]).decomposition()
        [
        Simple abelian subvariety 14aGH[15](1,28) of dimension 1 of JH(28,[15]),
        Simple abelian subvariety 14aGH[15](2,28) of dimension 1 of JH(28,[15]),
        Simple abelian subvariety 28aGH[15](1,28) of dimension 2 of JH(28,[15])
        ]
    """
    D = []
    for G in factor_modsym_space_new_factors(M):
        if G:
            # Compute the matrices of the degeneracy maps up.
            T = (M.level() // G[0].level()).divisors()
            degen = [G[0].ambient_module().degeneracy_map(M, t).matrix() for t in T]
            # Construct a matrix with rows the basis for all the factors
            # stacked on top of each other.  We just multiply this by each
            # degeneracy matrix to get the basis for the images of the
            # factors at higher level.  By doing matrix multiplies, we
            # save time over taking images of individual factors.
            matrix = G[0].basis_matrix()
            for A in G[1:]:
                matrix = matrix.stack(A.basis_matrix())

            # Compute the actual images
            ims = [matrix * z for z in degen]

            # Construct the corresponding subspaces at higher level.
            j = 0
            for (isog, A) in enumerate(G):
                d = A.dimension()
                if simple:
                    for i in range(len(T)):
                        V = ims[i].matrix_from_rows(range(j, j + d)).row_module()
                        W = M.submodule(V, check=False)
                        D.append((A.level(), isog, T[i], W))
                else:
                    V = sum(ims[i].matrix_from_rows(range(j, j + d)).row_module() for i in range(len(T)))
                    W = M.submodule(V, check=False)
                    D.append((A.level(), isog, None, W))
                j += d
    return Sequence(D, cr=True)


def modsym_lattices(M, factors):
    """
    Append lattice information to the output of
    simple_factorization_of_modsym_space.

    INPUT:


    -  ``M`` - modular symbols spaces

    -  ``factors`` - Sequence
       (simple_factorization_of_modsym_space)


    OUTPUT: sequence with more information for each factor (the
    lattice)

    EXAMPLES::

        sage: M = ModularSymbols(33)
        sage: factors = sage.modular.abvar.abvar.simple_factorization_of_modsym_space(M, simple=False)
        sage: sage.modular.abvar.abvar.modsym_lattices(M, factors)
        [
        (11, 0, None, Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field, Free module of degree 6 and rank 4 over Integer Ring
        Echelon basis matrix:
        [ 1  0  0  0 -1  2]
        [ 0  1  0  0 -1  1]
        [ 0  0  1  0 -2  2]
        [ 0  0  0  1 -1 -1]),
        (33, 0, None, Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field, Free module of degree 6 and rank 2 over Integer Ring
        Echelon basis matrix:
        [ 1  0  0 -1  0  0]
        [ 0  0  1  0  1 -1])
        ]
    """
    # 1. Change basis of everything to the ambient integral modular symbols space
    # 2. Clear denominator.
    # 3. Echelonize/saturate each factor
    if not factors:
        return factors

    D = []
    I = M.cuspidal_submodule().integral_structure().basis_matrix()
    A = factors[0][-1].basis_matrix()
    rows = [range(A.nrows())]
    for F in factors[1:]:
        mat = F[-1].basis_matrix()
        i = rows[-1][-1] + 1
        rows.append(range(i, i + mat.nrows()))
        A = A.stack(mat)
    X = I.solve_left(A)
    X, _ = X._clear_denom()
    for i, R in enumerate(rows):
        A = X.matrix_from_rows(R)
        A = copy(A.saturation())
        A.echelonize()
        D.append(tuple(list(factors[i]) + [A.row_module()]))
    return Sequence(D, cr=True)
