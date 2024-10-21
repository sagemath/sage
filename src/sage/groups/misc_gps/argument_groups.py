# sage_setup: distribution = sagemath-groups
r"""
Groups of elements representing (complex) arguments.

This includes

- :class:`RootsOfUnityGroup` (containing all roots of unity)

- :class:`UnitCircleGroup` (representing elements on the unit circle by
  `e^{2\pi\cdot\mathit{exponent}}`)

- :class:`ArgumentByElementGroup` (whose elements are defined via
  formal arguments by `e^{I\cdot\mathrm{arg}(\mathit{element})}`.

Use the factory :class:`ArgumentGroup` for creating such a group conveniently.

.. NOTE::

    One main purpose of such groups is in an
    :mod:`asymptotic ring's <sage.rings.asymptotic.asymptotic_ring>`
    :mod:`growth group <sage.rings.asymptotic.growth_group>`
    when an element like `z^n`
    (for some constant `z`) is split into
    `\lvert z \rvert^n \cdot e^{I\cdot \mathrm{arg}(z) n}`.
    (Note that the first factor determines the growth of that product,
    the second does not influence the growth.)

AUTHORS:

- Daniel Krenn (2018)

Classes and Methods
===================
"""
# ****************************************************************************
# Copyright (C) 2018 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element import MultiplicativeGroupElement
from sage.structure.factory import UniqueFactory
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp_by_eq_and_lt
from sage.structure.unique_representation import UniqueRepresentation
import sage.rings.abc


class AbstractArgument(MultiplicativeGroupElement):
    r"""
    An element of :class:`AbstractArgumentGroup`. This abstract class
    encapsulates an element of the parent's base, i.e. it can be seen
    as a wrapper class.

    INPUT:

    - ``parent`` -- a SageMath parent

    - ``element`` -- an element of parent's base

    - ``normalize`` -- boolean (default: ``True``)
    """

    def __init__(self, parent, element, normalize=True):
        r"""
        See :class:`AbstractArgument` for more information.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup, RootsOfUnityGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=1/2)  # indirect doctest
            e^(2*pi*0.500000000000000)
            sage: C(exponent=3/2)
            e^(2*pi*0.500000000000000)

            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0)
            1
            sage: U(exponent=1)
            1
            sage: U(exponent=2/3) == U(exponent=5/3)
            True
        """
        if parent is None:
            raise ValueError('parent must be provided')
        super().__init__(parent=parent)

        try:
            element = parent.base()(element)
        except (TypeError, ValueError) as e:
            from sage.rings.asymptotic.misc import combine_exceptions
            from sage.structure.element import parent as parent_function
            raise combine_exceptions(
                ValueError(
                    '{} ({}) is not in {}'.format(element,
                                                  parent_function(element),
                                                  parent.base())),
                e)

        if normalize:
            element = self._normalize_(element)
        self._element_ = element

    @staticmethod
    def _normalize_(element):
        r"""
        Normalize the given element.

        INPUT:

        - ``element`` -- an element of the parent's base

        OUTPUT: an element

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import AbstractArgument
            sage: AbstractArgument._normalize_(3/2)
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented in concrete realizations
        """
        raise NotImplementedError('only implemented in concrete realizations')

    def __hash__(self):
        r"""
        Return a hash value of this argument.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: C = UnitCircleGroup(RR)
            sage: hash(C(exponent=1/3))  # indirect doctest, random
            42
        """
        return hash((self.parent(), self._element_))

    def _symbolic_(self, R=None):
        r"""
        Return this argument as a symbolic expression.

        INPUT:

        - ``R`` -- (a subring of) the symbolic ring or ``None``.
          The output will be an element of ``R``. If ``None``,
          then the symbolic ring is used.

        OUTPUT: a symbolic expression

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import AbstractArgument
            sage: class MyArgument(AbstractArgument):
            ....:     @staticmethod
            ....:     def _normalize_(element):
            ....:         return element
            sage: a = MyArgument(ZZ, -1)
            sage: a._symbolic_()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented in concrete realizations
        """
        raise NotImplementedError('only implemented in concrete realizations')

    _richcmp_ = richcmp_by_eq_and_lt("_eq_", "_lt_")

    def _eq_(self, other):
        r"""
        Return whether this argument equals ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0) == U(exponent=1)
            True
            sage: U(exponent=2/3) == U(exponent=5/3)
            True
            sage: U(exponent=2/3) == U(exponent=-2/3)
            False

        ::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(I) == C(I)                                                          # needs sage.symbolic
            True

        As we do not have normalization in :class:`ArgumentByElement`,
        then following, although equal, is not equal::

            sage: C(I) == C(2*I)                                                        # needs sage.symbolic
            False
        """
        return self._element_ == other._element_

    def _lt_(self, other):
        r"""
        Raise an error since points on the unit circle are not comparable.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0) < U(exponent=0)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for the roots of unity 1 and 1
            sage: U(exponent=0) < U(exponent=1/2)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for the roots of unity 1 and -1
            sage: U(exponent=0) > U(exponent=1/2)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for the roots of unity -1 and 1
        """
        raise RuntimeError("cannot decide '<' "
                           "for the roots of unity "
                           "{} and {}".format(self, other))

    def _act_on_(self, other, is_left):
        r"""
        Return the action of this argument on ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: U(-1) * 4
            -4
            sage: _.parent()
            Symbolic Ring
            sage: 4 * U(-1)
            -4
            sage: _.parent()
            Symbolic Ring

            sage: P = Permutation([1,2,3])
            sage: U(-1) * P
            Traceback (most recent call last):
            ...
            TypeError: -1 (Group of Roots of Unity) cannot
            (left-)act on [1, 2, 3] (Standard permutations)
            > *previous* TypeError: no canonical coercion from
              Standard permutations to Symbolic Ring

        ::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(SR)                                        # needs sage.symbolic
            sage: C(-1) * 4                                                             # needs sage.symbolic
            -4
            sage: _.parent()
            Symbolic Ring
            sage: 4 * C(-1)                                                             # needs sage.symbolic
            -4
            sage: _.parent()
            Symbolic Ring
        """
        from sage.symbolic.ring import SymbolicRing, SR

        P = other.parent()
        S = P if isinstance(P, SymbolicRing) else SR
        try:
            other = S.coerce(other)
        except (TypeError, ValueError) as e:
            from sage.rings.asymptotic.misc import combine_exceptions
            raise combine_exceptions(
                TypeError('{} ({}) cannot ({}-)act on '
                          '{} ({})'.format(
                              self, self.parent(),
                              'left' if is_left else 'right',
                              other, P)),
                e)
        return self._symbolic_(S) * other

    def __abs__(self):
        r"""
        Return the absolute value of this argument which equals `1`.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: abs(U(exponent=1/4))  # indirect doctest
            1
            sage: _.parent()
            Integer Ring
        """
        from sage.rings.integer_ring import ZZ
        return ZZ.one()


class AbstractArgumentGroup(UniqueRepresentation, Parent):
    r"""
    A group whose elements represent (complex) arguments.

    INPUT:

    - ``base`` -- a SageMath parent

    - ``category`` -- a category
    """

    Element = AbstractArgument

    @staticmethod
    def __classcall__(cls, base, category=None):
        r"""
        See :class:`AbstractArgumentGroup` for more information.

        TESTS:

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: UnitCircleGroup(RR).category()  # indirect doctest
            Category of commutative groups
        """
        category = cls._determine_category_(category)
        return super().__classcall__(cls, base, category)

    @staticmethod
    def _determine_category_(category):
        r"""
        Return the category of this argument group.

        INPUT:

        - ``category`` -- a category or ``None`` (in which case the output
          equals ``category``)

        OUTPUT: a category

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: UnitCircleGroup._determine_category_(None)
            Category of commutative groups
            sage: UnitCircleGroup._determine_category_(Groups())
            Category of groups
        """
        if category is None:
            from sage.categories.groups import Groups
            category = Groups().Commutative()
        return category

    def __init__(self, base, category):
        r"""
        See :class:`AbstractArgumentGroup` for more information.

        TESTS:

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: UnitCircleGroup(RR).base()  # indirect doctest
            Real Field with 53 bits of precision
        """
        super().__init__(category=category, base=base)

    def __hash__(self):
        r"""
        Return a hash value of this argument group.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: hash(UnitCircleGroup(RR))  # indirect doctest, random
            42
        """
        return hash((self.__class__, self.base()))

    def _an_element_(self):
        r"""
        Return an element of this argument group.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: UnitCircleGroup(RR).an_element()  # indirect doctest
            e^(2*pi*0.000000000000000)
        """
        return self.element_class(self, self.base().an_element())


class UnitCirclePoint(AbstractArgument):
    r"""
    An element of :class:`UnitCircleGroup`
    which is `e^{2\pi\cdot\mathit{exponent}}`.

    INPUT:

    - ``parent`` -- a SageMath parent

    - ``exponent`` -- a number (of a subset of the reals)

    - ``normalize`` -- boolean (default: ``True``)
    """

    @staticmethod
    def _normalize_(exponent):
        r"""
        Normalize the given exponent so that it is in `[0,1)`.

        INPUT:

        - ``exponent`` -- an element of the parent's base

        OUTPUT: an element

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCirclePoint
            sage: UnitCirclePoint._normalize_(3/2)
            1/2
        """
        return exponent - exponent.floor()

    @property
    def exponent(self):
        r"""
        The exponent of this point on the unit circle.

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=4/3).exponent
            0.333333333333333
        """
        return self._element_

    def _repr_(self):
        r"""
        Return a representation string of this point on the unit circle.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=1/3)
            e^(2*pi*0.333333333333333)
        """
        return 'e^(2*pi*{})'.format(self.exponent)

    def _symbolic_(self, R=None):
        r"""
        Return this point on the unit circle as a symbolic expression.

        INPUT:

        - ``R`` -- (a subring of) the symbolic ring or ``None``.
          The output will be an element of ``R``. If ``None``,
          then the symbolic ring is used.

        OUTPUT: a symbolic expression

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=1/4)._symbolic_()
            e^(0.500000000000000*I*pi)
            sage: _.parent()
            Symbolic Ring

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: U(exponent=1/4)._symbolic_()
            I
            sage: _.parent()
            Symbolic Ring
        """
        from sage.functions.log import exp
        from sage.symbolic.ring import SR

        if R is None:
            R = SR

        return exp(2 * R('pi') * R('I') * self.exponent)

    def _mul_(self, other):
        r"""
        Return the product of this point on the unit circle and ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=0.3) * C(exponent=0.4)
            e^(2*pi*0.700000000000000)
        """
        P = self.parent()
        return P.element_class(P, self.exponent + other.exponent)

    def __pow__(self, exponent):
        r"""
        Return the power of this point on the unit circle
        to the given ``exponent``.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup, RootsOfUnityGroup

            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=0.1)^2
            e^(2*pi*0.200000000000000)
            sage: _.parent()
            Unit Circle Group with Exponents in
            Real Field with 53 bits of precision modulo ZZ
            sage: C(exponent=0.1)^QQ(2/1)
            e^(2*pi*0.200000000000000)
            sage: _.parent()
            Unit Circle Group with Exponents in
            Real Field with 53 bits of precision modulo ZZ

            sage: U = RootsOfUnityGroup()
            sage: a = U(exponent=1/7); a
            zeta7
            sage: a^(7/3)
            zeta3
            sage: _.parent()
            Group of Roots of Unity

            sage: U(exponent=1/3)^(0.25)
            e^(2*pi*0.0833333333333333)
            sage: _.parent()
            Unit Circle Group with Exponents in
            Real Field with 53 bits of precision modulo ZZ

            sage: U(exponent=1/3)^x
            (1/2*I*sqrt(3) - 1/2)^x
            sage: U(exponent=1/2)^x
            (-1)^x
            sage: U(exponent=1/4)^x
            I^x
            sage: U(exponent=1/4)^SR(8)
            1
        """
        from sage.symbolic.ring import SymbolicRing

        new_exponent = self.exponent * exponent
        parent = new_exponent.parent()
        if isinstance(new_exponent.parent(), SymbolicRing):
            return self._symbolic_(parent) ** exponent
        return self.parent()._create_element_in_extension_(new_exponent)

    def __invert__(self):
        r"""
        Return the inverse of this point on the unit circle.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: C = UnitCircleGroup(RR)
            sage: ~C(exponent=0.4)
            e^(2*pi*0.600000000000000)
            sage: C(1) / C(exponent=0.4)
            e^(2*pi*0.600000000000000)
            sage: C(exponent=0) / C(exponent=0.42)
            e^(2*pi*0.580000000000000)
        """
        P = self.parent()
        return P.element_class(P, -self.exponent)

    def is_one(self):
        r"""
        Return whether this point on the unit circle is `1`.

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: C = UnitCircleGroup(QQ)
            sage: C(exponent=0).is_one()
            True
            sage: C(exponent=1/2).is_one()
            False
            sage: C(exponent=2/3).is_one()
            False
            sage: C(exponent=42).is_one()
            True
        """
        return self.exponent == 0

    def is_minus_one(self):
        r"""
        Return whether this point on the unit circle is `-1`.

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: C = UnitCircleGroup(QQ)
            sage: C(exponent=0).is_minus_one()
            False
            sage: C(exponent=1/2).is_minus_one()
            True
            sage: C(exponent=2/3).is_minus_one()
            False
        """
        from sage.rings.rational_field import QQ
        return self.exponent == QQ((1, 2))


class UnitCircleGroup(AbstractArgumentGroup):
    r"""
    A group of points on the unit circle. These points are
    represented by `e^{2\pi\cdot\mathit{exponent}}`.

    INPUT:

    - ``base`` -- a SageMath parent representing a subset of the reals

    - ``category`` -- a category

    EXAMPLES::

        sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup

        sage: R = UnitCircleGroup(RR); R
        Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
        sage: R(exponent=2.42)
        e^(2*pi*0.420000000000000)

        sage: Q = UnitCircleGroup(QQ); Q
        Unit Circle Group with Exponents in Rational Field modulo ZZ
        sage: Q(exponent=6/5)
        e^(2*pi*1/5)
    """

    Element = UnitCirclePoint

    def _repr_(self):
        r"""
        Return a representation string of this unit circle group.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: UnitCircleGroup(RR)  # indirect doctest
            Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
        """
        return 'Unit Circle Group with Exponents in {} modulo ZZ'.format(self.base())

    def _repr_short_(self):
        r"""
        Return a short representation string of this unit circle group.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup
            sage: UnitCircleGroup(RR)._repr_short_()
            'UU_RR'
        """
        from sage.rings.asymptotic.misc import parent_to_repr_short
        s = parent_to_repr_short(self.base())
        if ' ' in s:
            s = '({})'.format(s)
        return 'UU_{}'.format(s)

    def _element_constructor_(self, data, exponent=None, **kwds):
        r"""
        Construct an element out of the given data.

        INPUT:

        - ``data`` -- an object

        - ``exponent`` -- a number (of a subset of the reals) or ``None``

        - ``kwds`` -- are passed on to element

        OUTPUT: a :class:`UnitCirclePoint`

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup, RootsOfUnityGroup, SignGroup
            sage: R = UnitCircleGroup(RR)
            sage: R(exponent=1/2)
            e^(2*pi*0.500000000000000)

            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0)
            1
            sage: U(exponent=1)
            1
            sage: U(exponent=1/2)
            -1
            sage: U(exponent=1/4)
            I
            sage: U(exponent=1/3)
            zeta3

            sage: C.<z> = CyclotomicField(6)                                            # needs sage.rings.number_field
            sage: z, U(z)                                                               # needs sage.rings.number_field
            (z, zeta6)
            sage: z^2, U(z^2)                                                           # needs sage.rings.number_field
            (z - 1, zeta3)

            sage: U(ZZ(-1))
            -1
            sage: U(QQ(-1))
            -1
            sage: U(int(-1))
            -1

            sage: S = SignGroup()
            sage: R(S(1))
            e^(2*pi*0.000000000000000)
            sage: R(S(-1))
            e^(2*pi*0.500000000000000)
            sage: U(S(1))
            1
            sage: U(S(-1))
            -1

            sage: R(U(exponent=1/3))
            e^(2*pi*0.333333333333333)

            sage: U(exponent=5/2, normalize=False)
            zeta2^5
        """
        from sage.groups.generic import discrete_log
        import sage.rings.abc
        from sage.rings.asymptotic.misc import combine_exceptions
        from sage.rings.rational_field import QQ

        if exponent is None:
            if isinstance(data, int) and data == 0:
                raise ValueError('no input specified')

            elif isinstance(data, self.element_class):
                if data.parent() == self:
                    return data
                exponent = data.exponent

            elif data == 1 or data == '1':
                exponent = 0

            elif data == -1 or data == '-1':
                exponent = QQ((1, 2))

            else:
                try:
                    P = data.parent()
                except AttributeError:
                    raise TypeError('{} is not in {}'.format(data, self))

                if isinstance(P, SignGroup):
                    if data.is_one():
                        exponent = 0
                    elif data.is_minus_one():
                        exponent = QQ((1, 2))

                elif isinstance(P, UnitCircleGroup):
                    exponent = data.exponent

                elif isinstance(P, sage.rings.abc.NumberField_cyclotomic):
                    zeta = P.gen()
                    n = zeta.multiplicative_order()
                    try:
                        exponent = QQ(discrete_log(data, zeta)) / QQ(n)
                    except ValueError as e:
                        raise combine_exceptions(
                            ValueError('{} is not in {}'.format(data, self)), e)

            if exponent is None:
                raise ValueError('{} is not in {}'.format(data, self))

        elif not isinstance(data, int) or data != 0:
            raise ValueError('input is ambigous: '
                             '{} as well as exponent={} '
                             'specified'.format(data, exponent))

        return self.element_class(self, exponent, **kwds)

    def _create_element_in_extension_(self, exponent):
        r"""
        Create an element in an extension of this unit circle group which
        is chosen according to the input ``exponent``.

        INPUT:

        - ``exponent`` -- the element data

        OUTPUT: an element

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup, RootsOfUnityGroup

            sage: C = UnitCircleGroup(QQ)
            sage: C._create_element_in_extension_(2.12).parent()                        # needs sage.rings.number_field
            Unit Circle Group with Exponents in
            Real Field with 53 bits of precision modulo ZZ

            sage: U = RootsOfUnityGroup()
            sage: U._create_element_in_extension_(2.12).parent()                        # needs sage.rings.number_field
            Unit Circle Group with Exponents in
            Real Field with 53 bits of precision modulo ZZ
        """
        if exponent.parent() is self.base():
            parent = self
        else:
            parent = ArgumentGroup(exponents=exponent.parent())
        return parent(exponent=exponent)

    def _coerce_map_from_(self, R):
        r"""
        Return whether ``R`` coerces into this unit circle group.

        INPUT:

        - ``R`` -- a parent

        OUTPUT: boolean

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import UnitCircleGroup, RootsOfUnityGroup, SignGroup
            sage: R = UnitCircleGroup(RR)
            sage: Q = UnitCircleGroup(QQ)
            sage: U = RootsOfUnityGroup()
            sage: S = SignGroup()
            sage: for A in (S, U, Q, R):  # indirect doctest
            ....:     for B in (S, U, Q, R):
            ....:         print('{} has {}coerce map from {}'.format(
            ....:             A,
            ....:             '' if A.has_coerce_map_from(B) else 'no ',
            ....:             B))
            Sign Group
              has coerce map from
              Sign Group
            Sign Group
              has no coerce map from
              Group of Roots of Unity
            Sign Group
              has no coerce map from
              Unit Circle Group with Exponents in Rational Field modulo ZZ
            Sign Group
              has no coerce map from
              Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
            Group of Roots of Unity
              has coerce map from
              Sign Group
            Group of Roots of Unity
              has coerce map from
              Group of Roots of Unity
            Group of Roots of Unity
              has coerce map from
              Unit Circle Group with Exponents in Rational Field modulo ZZ
            Group of Roots of Unity
              has no coerce map from
              Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
            Unit Circle Group with Exponents in Rational Field modulo ZZ
              has coerce map from
              Sign Group
            Unit Circle Group with Exponents in Rational Field modulo ZZ
              has coerce map from
              Group of Roots of Unity
            Unit Circle Group with Exponents in Rational Field modulo ZZ
              has coerce map from
              Unit Circle Group with Exponents in Rational Field modulo ZZ
            Unit Circle Group with Exponents in Rational Field modulo ZZ
              has no coerce map from
              Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
            Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
              has coerce map from
              Sign Group
            Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
              has coerce map from
              Group of Roots of Unity
            Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
              has coerce map from
              Unit Circle Group with Exponents in Rational Field modulo ZZ
            Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
              has coerce map from
              Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
        """
        if isinstance(R, UnitCircleGroup):
            return self.base().has_coerce_map_from(R.base())
        if isinstance(R, SignGroup):
            return True


class RootOfUnity(UnitCirclePoint):
    r"""
    A root of unity (i.e. an element of :class:`RootsOfUnityGroup`)
    which is `e^{2\pi\cdot\mathit{exponent}}` for a rational ``exponent``.
    """

    def exponent_numerator(self):
        r"""
        Return the numerator of the rational quotient in `[0,1)`
        representing the exponent of this root of unity.

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: a = U(exponent=2/3); a
            zeta3^2
            sage: a.exponent_numerator()
            2
        """
        return self.exponent.numerator()

    def exponent_denominator(self):
        r"""
        Return the denominator of the rational quotient in `[0,1)`
        representing the exponent of this root of unity.

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: a = U(exponent=2/3); a
            zeta3^2
            sage: a.exponent_denominator()
            3
        """
        return self.exponent.denominator()

    def _repr_(self):
        r"""
        Return a representation string of this root of unity.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0)
            1
            sage: U(exponent=1/2)
            -1
            sage: U(exponent=1/4)
            I
            sage: U(exponent=3/4)
            -I
            sage: U(exponent=1/3)
            zeta3
            sage: U(exponent=2/3)
            zeta3^2
        """
        from sage.rings.rational_field import QQ
        if self.exponent == 0:
            return '1'
        if self.exponent == QQ((1, 2)):
            return '-1'
        if self.exponent == QQ((1, 4)):
            return 'I'
        if self.exponent == QQ((3, 4)):
            return '-I'
        num = self.exponent_numerator()
        den = self.exponent_denominator()
        zeta = 'zeta{}'.format(den)
        if num == 1:
            return zeta
        return '{}^{}'.format(zeta, num)


class RootsOfUnityGroup(UnitCircleGroup):
    r"""
    The group of all roots of unity.

    INPUT:

    - ``category`` -- a category

    This is a specialized :class:`UnitCircleGroup` with base `\QQ`.

    EXAMPLES::

        sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
        sage: U = RootsOfUnityGroup(); U
        Group of Roots of Unity
        sage: U(exponent=1/4)
        I
    """

    Element = RootOfUnity

    @staticmethod
    def __classcall__(cls, category=None):
        r"""
        See :class:`RootsOfUnityGroup` for more information.

        TESTS:

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: RootsOfUnityGroup().category()  # indirect doctest
            Category of commutative groups
        """
        category = cls._determine_category_(category)
        return super(AbstractArgumentGroup, cls).__classcall__(cls, category)

    def __init__(self, category):
        r"""
        See :class:`RootsOfUnityGroup` for more information.

        TESTS:

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: RootsOfUnityGroup().base()  # indirect doctest
            Rational Field
        """
        from sage.rings.rational_field import QQ
        return super().__init__(base=QQ, category=category)

    def _repr_(self):
        r"""
        Return a representation string of this roots of unity group.

        TESTS:

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: RootsOfUnityGroup()  # indirect doctest
            Group of Roots of Unity
        """
        return 'Group of Roots of Unity'

    def _repr_short_(self):
        r"""
        Return a short representation string of this roots of unity group.

        TESTS:

            sage: from sage.groups.misc_gps.argument_groups import RootsOfUnityGroup
            sage: RootsOfUnityGroup()._repr_short_()
            'UU'
        """
        return 'UU'


class ArgumentByElement(AbstractArgument):
    r"""
    An element of :class:`ArgumentByElementGroup`.

    INPUT:

    - ``parent`` -- a SageMath parent

    - ``element`` -- a nonzero element of the parent's base

    - ``normalize`` -- boolean (default: ``True``)
    """

    def __init__(self, parent, element, normalize=True):
        r"""
        See :class:`ArgumentByElement` for more information.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(1+2*I)  # indirect doctest                                          # needs sage.symbolic
            e^(I*arg(1.00000000000000 + 2.00000000000000*I))
        """
        super().__init__(parent, element, normalize=normalize)
        if self._element_ == 0:
            raise ValueError('{} is not allowed'.format(element))

    @staticmethod
    def _normalize_(element):
        r"""
        Normalize the given element.

        This is the identity for :class:`ArgumentByElement`.

        INPUT:

        - ``element`` -- an element of the parent's base

        OUTPUT: an element

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElement
            sage: ArgumentByElement._normalize_(3/2)
            3/2
        """
        return element

    def _repr_(self):
        r"""
        Return a representation string of this argument by element.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(2+3*I)  # indirect doctest                                          # needs sage.symbolic
            e^(I*arg(2.00000000000000 + 3.00000000000000*I))
        """
        return 'e^(I*arg({}))'.format(self._element_)

    def _symbolic_(self, R=None):
        r"""
        Return this argument by element as a symbolic expression.

        INPUT:

        - ``R`` -- (a subring of) the symbolic ring or ``None``.
          The output will be an element of ``R``. If ``None``,
          then the symbolic ring is used.

        OUTPUT: a symbolic expression

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(ZZ)
            sage: C(-2)._symbolic_()                                                    # needs sage.symbolic
            -1
            sage: _.parent()
            Symbolic Ring
        """
        from sage.functions.log import exp
        from sage.functions.other import arg
        from sage.symbolic.ring import SR

        if R is None:
            R = SR

        return exp(R('I') * arg(self._element_))

    def _mul_(self, other):
        r"""
        Return the product of this argument by element with ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(I) * C(1 + I)  # indirect doctest                                   # needs sage.symbolic
            e^(I*arg(-1.00000000000000 + 1.00000000000000*I))
        """
        P = self.parent()
        return P.element_class(P, self._element_ * other._element_)

    def __pow__(self, exponent):
        r"""
        Return the power of this argument by element
        to the given ``exponent``.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(I)^5  # indirect doctest                                            # needs sage.symbolic
            e^(I*arg(1.00000000000000*I))
            sage: _.parent()
            Unit Circle Group with Argument of Elements in
            Complex Field with 53 bits of precision
            sage: C(1+I)^3  # indirect doctest                                          # needs sage.symbolic
            e^(I*arg(-2.00000000000000 + 2.00000000000000*I))
            sage: _.parent()
            Unit Circle Group with Argument of Elements in
            Complex Field with 53 bits of precision

            sage: C = ArgumentByElementGroup(RR)
            sage: C(0.42)^CC(2.4)
            e^(I*arg(0.124680431591996))
            sage: _.parent()
            Unit Circle Group with Argument of Elements in
            Complex Field with 53 bits of precision

            sage: C = ArgumentByElementGroup(QQ)
            sage: a = C(-20)^x; a
            (-1)^x
            sage: a.parent()
            Symbolic Ring
        """
        from sage.rings.abc import SymbolicRing

        element = self._element_ ** exponent
        parent = element.parent()
        if isinstance(parent, SymbolicRing):
            return self._symbolic_(parent) ** exponent
        return self.parent()._create_element_in_extension_(element)

    def __invert__(self):
        r"""
        Return the inverse of this argument by element.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: ~C(I)  # indirect doctest                                             # needs sage.symbolic
            e^(I*arg(-1.00000000000000*I))
        """
        P = self.parent()
        return P.element_class(P, ~self._element_)


class ArgumentByElementGroup(AbstractArgumentGroup):
    r"""
    A group of (complex) arguments. The arguments are represented
    by a the formal argument of an element, i.e.,
    by `\mathrm{arg}(\mathit{element})`.

    INPUT:

    - ``base`` -- a SageMath parent representing a subset of the complex plane

    - ``category`` -- a category

    EXAMPLES::

        sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
        sage: C = ArgumentByElementGroup(CC); C
        Unit Circle Group with Argument of Elements in
        Complex Field with 53 bits of precision
        sage: C(1 + 2*I)                                                                # needs sage.symbolic
        e^(I*arg(1.00000000000000 + 2.00000000000000*I))
    """

    Element = ArgumentByElement

    def _repr_(self):
        r"""
        Return a representation string of this argument by element group.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: ArgumentByElementGroup(CC)  # indirect doctest
            Unit Circle Group with Argument of Elements in
            Complex Field with 53 bits of precision
        """
        return 'Unit Circle Group with Argument of Elements in {}'.format(self.base())

    def _repr_short_(self):
        r"""
        Return a short representation string of this  argument by element group.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: ArgumentByElementGroup(CC)._repr_short_()
            'Arg_CC'
        """
        from sage.rings.asymptotic.misc import parent_to_repr_short, repr_op
        return repr_op('Arg', '_', parent_to_repr_short(self.base()))

    def _element_constructor_(self, data, **kwds):
        r"""
        Construct an element out of the given data.

        INPUT:

        - ``data`` -- an object

        - ``kwds`` -- are passed on to element

        OUTPUT: a :class:`ArgumentByElement`

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(CC)
            sage: C(1 + 2*I)  # indirect doctest                                        # needs sage.symbolic
            e^(I*arg(1.00000000000000 + 2.00000000000000*I))
            sage: C(1)
            e^(I*arg(1.00000000000000))
            sage: C(-1)
            e^(I*arg(-1.00000000000000))
            sage: C(ZZ(-1))
            e^(I*arg(-1.00000000000000))
            sage: C(QQ(-1))
            e^(I*arg(-1.00000000000000))
            sage: C(int(-1))
            e^(I*arg(-1.00000000000000))
            sage: C('-1')
            e^(I*arg(-1.00000000000000))

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: C(S(1))
            e^(I*arg(1.00000000000000))
            sage: C(S(-1))
            e^(I*arg(-1.00000000000000))
        """
        from sage.structure.element import parent

        if isinstance(data, int) and data == 0:
            raise ValueError('no input specified')

        elif isinstance(data, self.element_class):
            if data.parent() == self:
                return data
            element = data._element_

        elif data == '1':
            element = 1

        elif data == '-1':
            element = -1

        else:
            P = parent(data)

            if isinstance(P, SignGroup):
                if data.is_one():
                    element = 1
                elif data.is_minus_one():
                    element = -1

            else:
                element = data

        return self.element_class(self, element, **kwds)

    def _create_element_in_extension_(self, element):
        r"""
        Create an element in an extension of this
        argument by element group which
        is chosen according to the input ``element``.

        INPUT:

        - ``element`` -- the element data

        OUTPUT: an element

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: C = ArgumentByElementGroup(ZZ)
            sage: C._create_element_in_extension_(2/3).parent()
            Unit Circle Group with Argument of Elements in
            Rational Field
            sage: C._create_element_in_extension_(0.23).parent()
            Unit Circle Group with Argument of Elements in
            Real Field with 53 bits of precision

            sage: C = ArgumentByElementGroup(RR)
            sage: C._create_element_in_extension_(CC(0.23)).parent()
            Unit Circle Group with Argument of Elements in
            Complex Field with 53 bits of precision
        """
        if element.parent() is self.base():
            parent = self
        else:
            parent = ArgumentByElementGroup(element.parent())
        return parent(element)

    def _coerce_map_from_(self, R):
        r"""
        Return whether ``R`` coerces into this argument by element group.

        INPUT:

        - ``R`` -- a parent

        OUTPUT: boolean

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentByElementGroup
            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: C = ArgumentByElementGroup(ZZ)
            sage: S = SignGroup()
            sage: C.has_coerce_map_from(S)  # indirect doctest
            True
        """
        if isinstance(R, SignGroup):
            return True


class Sign(AbstractArgument):
    r"""
    An element of :class:`SignGroup`.

    INPUT:

    - ``parent`` -- a SageMath parent

    - ``element`` -- a nonzero element of the parent's base

    - ``normalize`` -- boolean (default: ``True``)
    """

    def __init__(self, parent, element, normalize=True):
        r"""
        See :class:`Sign` for more information.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S.an_element()  # indirect doctest
            -1
        """
        super().__init__(parent, int(element), normalize=normalize)
        if self._element_ not in (-1, 1):
            raise ValueError('{} is not allowed '
                             '(only -1 or 1 is)'.format(element))

    @staticmethod
    def _normalize_(element):
        r"""
        Normalize the given element.

        This is the identity for :class:`Sign`.

        INPUT:

        - ``element`` -- an element of the parent's base

        OUTPUT: an element

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import Sign
            sage: Sign._normalize_(True)
            True
        """
        return element

    def _repr_(self):
        r"""
        Return a representation string of this sign.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S(1)
            1
            sage: S(-1)
            -1
        """
        return repr(self._element_)

    def _mul_(self, other):
        r"""
        Return the product of this sign with ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S(1) * S(-1)  # indirect doctest
            -1
        """
        P = self.parent()
        return P.element_class(P, self._element_ * other._element_)

    def __pow__(self, exponent):
        r"""
        Return the power of this sign to the given ``exponent``.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S(-1)^4  # indirect doctest
            1
            sage: S(-1)^3  # indirect doctest
            -1
        """
        P = self.parent()
        return P.element_class(P, self._element_ ** exponent)

    def __invert__(self):
        r"""
        Return the inverse of this sign.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: ~S(-1)  # indirect doctest
            -1
            sage: _.parent()
            Sign Group
        """
        return self

    def _act_on_(self, other, is_left):
        r"""
        Return the action of this sign on ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S(-1) * 4
            -4
            sage: _.parent()
            Integer Ring
            sage: 4 * S(-1)
            -4
            sage: _.parent()
            Integer Ring

            sage: S(-1) * ZZ(4)
            -4
            sage: _.parent()
            Integer Ring
            sage: S(-1) * int(4)
            -4
            sage: type(_)
            <class 'int'>
            sage: S(-1) * QQ(4)
            -4
            sage: _.parent()
            Rational Field
            sage: S(-1) * RR(4)
            -4.00000000000000
            sage: _.parent()
            Real Field with 53 bits of precision
            sage: S(-1) * CC(4)
            -4.00000000000000
            sage: _.parent()
            Complex Field with 53 bits of precision
            sage: S(-1) * SR.var('x')                                                   # needs sage.symbolic
            -x
            sage: _.parent()
            Symbolic Ring

        ::

            sage: P = Permutation([1,2,3])
            sage: S(-1) * P
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent
            for unary -: 'Standard permutations'
        """
        if self.is_one():
            return other
        if self.is_minus_one():
            return -other

    def is_one(self):
        r"""
        Return whether this sign is `1`.

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S(-1).is_one()
            False
            sage: S(1).is_one()
            True
        """
        return self._element_ == 1

    def is_minus_one(self):
        r"""
        Return whether this sign is `-1`.

        EXAMPLES::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S(1).is_minus_one()
            False
            sage: S(-1).is_minus_one()
            True
        """
        return self._element_ == -1


class SignGroup(AbstractArgumentGroup):
    r"""
    A group of the signs `-1` and `1`.

    INPUT:

    - ``category`` -- a category

    EXAMPLES::

        sage: from sage.groups.misc_gps.argument_groups import SignGroup
        sage: S = SignGroup(); S
        Sign Group
        sage: S(-1)
        -1
    """

    Element = Sign

    @staticmethod
    def __classcall__(cls, category=None):
        r"""
        See :class:`SignGroup` for more information.

        TESTS:

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S.category()  # indirect doctest
            Category of finite commutative groups
        """
        category = cls._determine_category_(category).Finite()
        return super(AbstractArgumentGroup, cls).__classcall__(
            cls, category)

    def __init__(self, category):
        r"""
        See :class:`SignGroup` for more information.

        TESTS:

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S.base()  # indirect doctest
            <class 'int'>
        """
        return super().__init__(base=int, category=category)

    def _repr_(self):
        r"""
        Return a representation string of the sign group.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: SignGroup()  # indirect doctest
            Sign Group
        """
        return 'Sign Group'

    def _repr_short_(self):
        r"""
        Return a short representation string of this sign group.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S._repr_short_()
            'Signs'
        """
        return 'Signs'

    def _an_element_(self):
        r"""
        Return a short representation string of this sign group.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S.an_element()
            -1
        """
        return self.element_class(self, -1)

    def _element_constructor_(self, data):
        r"""
        Construct an element out of the given data.

        INPUT:

        - ``data`` -- an object

        OUTPUT: a :class:`Sign`

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import SignGroup
            sage: S = SignGroup()
            sage: S(1)  # indirect doctest
            1
            sage: S(-1)  # indirect doctest
            -1
        """
        if isinstance(data, int) and data == 0:
            raise ValueError('no input specified')

        elif isinstance(data, self.element_class):
            if data.parent() == self:
                return data
            element = data._element_

        else:
            element = data

        return self.element_class(self, element)


class ArgumentGroupFactory(UniqueFactory):
    r"""
    A factory for creating argument groups.

    INPUT:

    - ``data`` -- an object

      The factory will analyze ``data`` and interpret it as
      ``specification`` or ``domain``.

    - ``specification`` -- string

      The following is possible:

      - ``'Signs'`` give the :class:`SignGroup`

      - ``'UU'`` give the :class:`RootsOfUnityGroup`

      - ``'UU_P'``, where ``'P'`` is
        a string representing a SageMath parent which is interpreted as
        ``exponents``

      - ``'Arg_P'``, where ``'P'`` is
        a string representing a SageMath parent which is interpreted as
        ``domain``

    - ``domain`` -- a SageMath parent representing a subset of the complex plane;
      an instance of :class:`ArgumentByElementGroup` will be created with the given
      ``domain``

    - ``exponents`` -- a SageMath parent representing a subset of the reals;
      an instance of :class`UnitCircleGroup` will be created with the given
      ``exponents``

    Exactly one of ``data``, ``specification``, ``exponents`` has to be provided.

    Further keyword parameters will be carried on to the initialization of
    the group.

    EXAMPLES::

        sage: from sage.groups.misc_gps.argument_groups import ArgumentGroup

        sage: ArgumentGroup('UU')                                                       # needs sage.rings.number_field
        Group of Roots of Unity

        sage: # needs sage.rings.number_field
        sage: ArgumentGroup(ZZ)
        Sign Group
        sage: ArgumentGroup(QQ)
        Sign Group
        sage: ArgumentGroup('UU_QQ')
        Group of Roots of Unity
        sage: ArgumentGroup(AA)
        Sign Group

        sage: ArgumentGroup(RR)                                                         # needs sage.rings.number_field
        Sign Group
        sage: ArgumentGroup('Arg_RR')                                                   # needs sage.rings.number_field
        Sign Group
        sage: ArgumentGroup(RIF)                                                        # needs sage.rings.real_interval_field
        Sign Group
        sage: ArgumentGroup(RBF)
        Sign Group

        sage: ArgumentGroup(CC)                                                         # needs sage.rings.number_field
        Unit Circle Group with Exponents in
        Real Field with 53 bits of precision modulo ZZ
        sage: ArgumentGroup('Arg_CC')                                                   # needs sage.rings.number_field
        Unit Circle Group with Exponents in
        Real Field with 53 bits of precision modulo ZZ
        sage: ArgumentGroup(CIF)
        Unit Circle Group with Exponents in
        Real Interval Field with 53 bits of precision modulo ZZ
        sage: ArgumentGroup(CBF)
        Unit Circle Group with Exponents in
        Real ball field with 53 bits of precision modulo ZZ

        sage: ArgumentGroup(CyclotomicField(3))                                         # needs sage.rings.number_field
        Unit Circle Group with Argument of Elements in
        Cyclotomic Field of order 3 and degree 2
    """
    def create_key_and_extra_args(self,
                                  data=None,
                                  specification=None,
                                  domain=None,
                                  exponents=None,
                                  **kwds):
        r"""
        Normalize the input.

        See :class:`ArgumentGroupFactory` for a description and examples.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentGroup

            sage: # needs sage.rings.number_field
            sage: ArgumentGroup(specification='UU')
            Group of Roots of Unity
            sage: ArgumentGroup('UU') is ArgumentGroup(exponents=QQ)  # indirect doctest
            True
            sage: ArgumentGroup('Arg_CC') is ArgumentGroup(exponents=RR)  # indirect doctest
            True
            sage: ArgumentGroup('Arg_CC') is ArgumentGroup(domain=CC)  # indirect doctest
            True
        """
        from sage.rings.integer_ring import ZZ
        from sage.misc.misc import exactly_one_is_true
        from sage.rings.qqbar import AA
        from sage.rings.rational_field import QQ

        if not exactly_one_is_true(
                (data is not None,
                 specification is not None,
                 domain is not None,
                 exponents is not None)):
            raise ValueError(
                'input ambigous: ' +
                ', '.join('{}={}'.format(s, v) for s, v in
                          [('data', data), ('specification', specification),
                           ('domain', domain), ('exponents', exponents)]
                          if v is not None))

        if data is not None:
            if isinstance(data, str):
                specification = data
            else:
                domain = data

        if specification is not None:
            if specification == 'UU':
                return (RootsOfUnityGroup, ()), kwds
            if specification == 'Signs':
                return (SignGroup, ()), kwds
            elif specification.startswith('UU_'):
                from sage.rings.asymptotic.misc import repr_short_to_parent
                exponents = repr_short_to_parent(specification[3:])
            elif specification.startswith('Arg_') or specification.startswith('arg_'):
                from sage.rings.asymptotic.misc import repr_short_to_parent
                domain = repr_short_to_parent(specification[4:])
            else:
                raise ValueError('unknown specification {}'.format(specification))

        if domain is not None:
            if domain in (ZZ, QQ, AA) \
               or isinstance(domain, (sage.rings.abc.RealField,
                                      sage.rings.abc.RealIntervalField,
                                      sage.rings.abc.RealBallField)):
                return (SignGroup, ()), kwds
            elif isinstance(domain, (sage.rings.abc.ComplexField,
                                     sage.rings.abc.ComplexIntervalField,
                                     sage.rings.abc.ComplexBallField)):
                return (UnitCircleGroup, (domain._real_field(),)), kwds
            else:
                return (ArgumentByElementGroup, (domain,)), kwds

        elif exponents is not None:
            if exponents == ZZ:
                return (SignGroup, ()), kwds
            elif exponents == QQ:
                return (RootsOfUnityGroup, ()), kwds
            else:
                return (UnitCircleGroup, (exponents,)), kwds

    def create_object(self, version, key, **kwds):
        r"""
        Create an object from the given arguments.

        TESTS::

            sage: from sage.groups.misc_gps.argument_groups import ArgumentGroup
            sage: ArgumentGroup('UU')  # indirect doctest                               # needs sage.rings.number_field
            Group of Roots of Unity
        """
        cls, args = key
        return cls(*args, **kwds)


ArgumentGroup = ArgumentGroupFactory('sage.groups.misc_gps.argument_groups.ArgumentGroup')
r"""
A factory for argument groups.

This is an instance of :class:`ArgumentGroupFactory` whose documentation
provides more details.
"""
