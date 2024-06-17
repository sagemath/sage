r"""
Groups of imaginary elements

.. NOTE::

    One main purpose of such groups is in an
    :mod:`asymptotic ring's <sage.rings.asymptotic.asymptotic_ring>`
    :mod:`growth group <sage.rings.asymptotic.growth_group>`
    when an element like `n^z`
    (for some constant `z`) is split into
    `n^{\Re z + I\Im z}`.
    (Note that the first summand in the exponent determines the growth,
    the second does not influence the growth.)

AUTHORS:

- Daniel Krenn (2018)

Classes and Methods
===================
"""
#*****************************************************************************
# Copyright (C) 2018 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import AdditiveGroupElement
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp_by_eq_and_lt
from sage.structure.unique_representation import UniqueRepresentation


class ImaginaryElement(AdditiveGroupElement):
    r"""
    An element of :class:`ImaginaryGroup`.

    INPUT:

    - ``parent`` -- a SageMath parent

    - ``imag`` -- an element of parent's base
    """

    def __init__(self, parent, imag):
        r"""
        See :class:`ImaginaryGroup` for more information.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=-42)  # indirect doctest
            -42*I
        """
        if parent is None:
            raise ValueError('parent must be provided')
        super().__init__(parent=parent)

        try:
            self._imag_ = parent.base()(imag)
        except (TypeError, ValueError) as e:
            from sage.rings.asymptotic.misc import combine_exceptions
            from sage.structure.element import parent as parent_function
            raise combine_exceptions(
                ValueError(
                    '{} ({}) is not in {}'.format(imag,
                                                  parent_function(imag),
                                                  parent.base())),
                e)

    def imag(self):
        r"""
        Return the imaginary part of this imaginary element.

        EXAMPLES::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(I).imag()                                                           # needs sage.symbolic
            1
            sage: imag_part(J(I))  # indirect doctest                                   # needs sage.symbolic
            1
        """
        return self._imag_

    def real(self):
        r"""
        Return the real part (`=0`) of this imaginary element.

        EXAMPLES::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(I).real()                                                           # needs sage.symbolic
            0
            sage: real_part(J(I))  # indirect doctest                                   # needs sage.symbolic
            0
        """
        return self.parent().base().zero()

    def __hash__(self):
        r"""
        Return a hash value of this imaginary element.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: hash(J(I))  # indirect doctest, random                                # needs sage.symbolic
            42
        """
        return hash((self.parent(), self._imag_))

    _richcmp_ = richcmp_by_eq_and_lt("_eq_", "_lt_")

    def _eq_(self, other):
        r"""
        Return whether this imaginary part equals ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=3) == J(imag=4)
            False
            sage: J(imag=42) == J(imag=42)
            True
        """
        return self._imag_ == other._imag_

    def _lt_(self, other):
        r"""
        Raise an error since imaginary elements are not comparable.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=1) < J(imag=1)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for imaginary elements I and I
            sage: J(imag=1) < J(imag=2)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for imaginary elements I and 2*I
            sage: J(imag=1) > J(imag=2)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for imaginary elements 2*I and I
        """
        raise RuntimeError("cannot decide '<' "
                           "for imaginary elements "
                           "{} and {}".format(self, other))

    def _repr_(self):
        r"""
        Return a representation string of this imaginary element.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=0)  # indirect doctest
            0
            sage: J(imag=1)  # indirect doctest
            I
            sage: J(imag=-1)  # indirect doctest
            -I
            sage: J(imag=42)  # indirect doctest
            42*I
            sage: J(imag=-42)  # indirect doctest
            -42*I
        """
        from sage.rings.asymptotic.misc import repr_op
        if self._imag_ == 0:
            return '0'
        if self._imag_ == 1:
            return 'I'
        if self._imag_ == -1:
            return '-I'
        return repr_op(self._imag_, '*', 'I')

    def _add_(self, other):
        r"""
        Return the sum of this imaginary element and ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=1) + J(imag=-2)  # indirect doctest
            -I
        """
        P = self.parent()
        return P.element_class(P, self._imag_ + other._imag_)

    def _sub_(self, other):
        r"""
        Return the difference of this imaginary element and ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=1) - J(imag=-2)  # indirect doctest
            3*I
        """
        P = self.parent()
        return P.element_class(P, self._imag_ - other._imag_)

    def __neg__(self):
        r"""
        Return the negative of this imaginary element.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: -J(imag=-2)  # indirect doctest
            2*I
        """
        P = self.parent()
        return P.element_class(P, -self._imag_)


class ImaginaryGroup(UniqueRepresentation, Parent):
    r"""
    A group whose elements are purely imaginary.

    INPUT:

    - ``base`` -- a SageMath parent

    - ``category`` -- a category

    EXAMPLES::

        sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
        sage: J = ImaginaryGroup(ZZ)
        sage: J(0)
        0
        sage: J(imag=100)
        100*I
        sage: J(3*I)                                                                    # needs sage.symbolic
        3*I
        sage: J(1 + 2*I)                                                                # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ValueError: 2*I + 1 is not in
        Imaginary Group over Integer Ring
        because it is not purely imaginary
    """

    Element = ImaginaryElement

    @staticmethod
    def __classcall__(cls, base, category=None):
        r"""
        See :class:`ImaginaryGroup` for more information.

        TESTS:

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)  # indirect doctest
            sage: J.category()
            Category of commutative additive groups
        """
        category = cls._determine_category_(category)
        return super().__classcall__(cls, base, category)

    @staticmethod
    def _determine_category_(category):
        r"""
        Return the category of this imaginary group.

        INPUT:

        - ``category`` -- a category or ``None`` (in which case the output
          equals ``category``)

        OUTPUT: a category

        EXAMPLES::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: from sage.categories.additive_groups import AdditiveGroups
            sage: ImaginaryGroup._determine_category_(None)
            Category of commutative additive groups
            sage: ImaginaryGroup._determine_category_(AdditiveGroups())
            Category of additive groups
        """
        if category is None:
            from sage.categories.additive_groups import AdditiveGroups
            category = AdditiveGroups().AdditiveCommutative()
        return category

    def __init__(self, base, category):
        r"""
        See :class:`ImaginaryGroup` for more information.

        TESTS:

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)  # indirect doctest
        """
        super().__init__(category=category, base=base)

    def __hash__(self):
        r"""
        Return a hash value of this imaginary group.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: hash(J)  # indirect doctest, random
            42
        """
        return hash((self.__class__, self.base()))

    def _an_element_(self):
        r"""
        Return an element of this imaginary group.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J.an_element()  # indirect doctest
            I
        """
        return self.element_class(self, self.base().an_element())

    def _repr_(self):
        r"""
        Return a representation string of this imaginary group.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J  # indirect doctest
            Imaginary Group over Integer Ring
        """
        return 'Imaginary Group over {}'.format(self.base())

    def _repr_short_(self):
        r"""
        Return a short representation string of this imaginary group.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J._repr_short_()
            'ZZ*I'
        """
        from sage.rings.asymptotic.misc import parent_to_repr_short, repr_op
        return repr_op(parent_to_repr_short(self.base()), '*', 'I')

    def _element_constructor_(self, data, imag=None):
        r"""
        Construct an element out of the given data.

        INPUT:

        - ``data`` -- an object

        - ``imag`` -- a number (of a subset of the reals) or ``None``

        OUTPUT: a :class:`ImaginaryElement`

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(0)
            0
            sage: J('0')
            0
            sage: J('-I')
            -I
            sage: J(imag=100)
            100*I

            sage: # needs sage.symbolic
            sage: J(3*I)
            3*I
            sage: J(1+2*I)
            Traceback (most recent call last):
            ...
            ValueError: 2*I + 1 is not in
            Imaginary Group over Integer Ring
            because it is not purely imaginary
            sage: i = CC(I)
            sage: J(3*i)
            3*I
            sage: J(1+2*i)
            Traceback (most recent call last):
            ...
            ValueError: 1.00000000000000 + 2.00000000000000*I is not in
            Imaginary Group over Integer Ring
            because it is not purely imaginary
            sage: i = CBF(I)
            sage: J(3*i)
            3*I
            sage: J(1+2*i)
            Traceback (most recent call last):
            ...
            ValueError: 1.000000000000000 + 2.000000000000000*I is not in
            Imaginary Group over Integer Ring
            because it is not purely imaginary
            sage: i = CIF(I)
            sage: J(3*i)
            3*I
            sage: J(1+2*i)
            Traceback (most recent call last):
            ...
            ValueError: 1 + 2*I is not in
            Imaginary Group over Integer Ring
            because it is not purely imaginary

            sage: J(x)
            Traceback (most recent call last):
            ...
            ValueError: x is not in
            Imaginary Group over Integer Ring
            because it is not purely imaginary
            sage: J(imag=1/2)
            Traceback (most recent call last):
            ...
            ValueError: 1/2 (Rational Field) is not in Integer Ring
            > *previous* TypeError: no conversion of this rational to integer
        """
        if imag is None:
            if isinstance(data, int) and data == 0:
                imag = data

            elif isinstance(data, self.element_class):
                if data.parent() == self:
                    return data
                imag = data._imag_

            elif data == 0 or data == '0':
                imag = 0

            elif data == 'I':
                imag = 1

            elif data == '-I':
                imag = -1

            else:
                try:
                    if data.real() == 0:
                        imag = data.imag()
                    else:
                        raise ValueError(
                            '{} is not in {} because it is not '
                            'purely imaginary'.format(data, self))
                except AttributeError:
                    pass

            if imag is None:
                raise ValueError('{} is not in {}'.format(data, self))

        elif not isinstance(data, int) or data != 0:
            raise ValueError('input is ambigous: '
                             '{} as well as imag={} '
                             'specified'.format(data, imag))

        return self.element_class(self, imag)
