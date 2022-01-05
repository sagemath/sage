r"""
(Asymptotic) Term Monoids

This module implements asymptotic term monoids. The elements of these
monoids are used behind the scenes when performing calculations in an
:doc:`asymptotic ring <asymptotic_ring>`.

The monoids build upon the (asymptotic) growth groups. While growth
elements only model the growth of a function as it tends towards
infinity (or tends towards another fixed point; see
:doc:`growth_group` for more details), an
asymptotic term additionally specifies its "type" and performs the
actual arithmetic operations (multiplication and partial
addition/absorption of terms).

Besides an abstract base term :class:`GenericTerm`, this module
implements the following types of terms:

- :class:`OTerm` -- `O`-terms at infinity, see
  :wikipedia:`Big_O_notation`.
- :class:`TermWithCoefficient` -- abstract base class for
  asymptotic terms with coefficients.
- :class:`ExactTerm` -- this class represents a growth element
  multiplied with some non-zero coefficient from a coefficient ring.

A characteristic property of asymptotic terms is that some terms are
able to "absorb" other terms (see
:meth:`~sage.rings.asymptotic.term_monoid.GenericTerm.absorb`). For
instance, `O(x^2)` is able to absorb `O(x)` (with result
`O(x^2)`), and `3\cdot x^5` is able to absorb `-2\cdot x^5` (with result
`x^5`). Essentially, absorption can be interpreted as the
addition of "compatible" terms (partial addition).


.. _term_absorption:

Absorption of Asymptotic Terms
==============================

A characteristic property of asymptotic terms is that some terms are
able to "absorb" other terms. This is realized with the method
:meth:`~sage.rings.asymptotic.term_monoid.GenericTerm.absorb`.

For instance, `O(x^2)` is able to absorb `O(x)` (with result
`O(x^2)`). This is because the functions bounded by linear growth
are bounded by quadratic growth as well. Another example would be
that `3x^5` is able to absorb `-2x^5` (with result `x^5`), which
simply corresponds to addition.

Essentially, absorption can be interpreted
as the addition of "compatible" terms (partial addition).

We want to show step by step which terms can be absorbed
by which other terms. We start by defining the necessary
term monoids and some terms::

    sage: from sage.rings.asymptotic.term_monoid import OTermMonoid, ExactTermMonoid
    sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
    sage: from sage.rings.asymptotic.growth_group import GrowthGroup
    sage: G = GrowthGroup('x^ZZ'); x = G.gen()
    sage: OT = OTermMonoid(TermMonoid, growth_group=G, coefficient_ring=QQ)
    sage: ET = ExactTermMonoid(TermMonoid, growth_group=G, coefficient_ring=QQ)
    sage: ot1 = OT(x); ot2 = OT(x^2)
    sage: et1 = ET(x^2, 2)

- Because of the definition of `O`-terms (see
  :wikipedia:`Big_O_notation`), :class:`OTerm` are able to absorb all
  other asymptotic terms with weaker or equal growth. In our
  implementation, this means that :class:`OTerm` is able to absorb
  other :class:`OTerm`, as well as :class:`ExactTerm`, as long as the
  growth of the other term is less than or equal to the growth of this
  element::

      sage: ot1, ot2
      (O(x), O(x^2))
      sage: ot1.can_absorb(ot2), ot2.can_absorb(ot1)
      (False, True)
      sage: et1
      2*x^2
      sage: ot1.can_absorb(et1)
      False
      sage: ot2.can_absorb(et1)
      True

  The result of this absorption always is the dominant
  (absorbing) :class:`OTerm`::

      sage: ot1.absorb(ot1)
      O(x)
      sage: ot2.absorb(ot1)
      O(x^2)
      sage: ot2.absorb(et1)
      O(x^2)

  These examples correspond to `O(x) + O(x) = O(x)`,
  `O(x^2) + O(x) = O(x^2)`, and `O(x^2) + 2x^2 = O(x^2)`.

- :class:`ExactTerm` can only absorb another
  :class:`ExactTerm` if the growth coincides with the
  growth of this element::

      sage: et1.can_absorb(ET(x^2, 5))
      True
      sage: any(et1.can_absorb(t) for t in [ot1, ot2])
      False

  As mentioned above, absorption directly corresponds
  to addition in this case::

      sage: et1.absorb(ET(x^2, 5))
      7*x^2

  When adding two exact terms, they might cancel out.
  For technical reasons, ``None`` is returned in this
  case::

      sage: ET(x^2, 5).can_absorb(ET(x^2, -5))
      True
      sage: ET(x^2, 5).absorb(ET(x^2, -5)) is None
      True

- The abstract base terms :class:`GenericTerm` and
  :class:`TermWithCoefficient` can neither absorb any
  other term, nor be absorbed by any other term.

If ``absorb`` is called on a term that cannot be absorbed, an
:python:`ArithmeticError<library/exceptions.html#exceptions.ArithmeticError>`
is raised::

    sage: ot1.absorb(ot2)
    Traceback (most recent call last):
    ...
    ArithmeticError: O(x) cannot absorb O(x^2)

This would only work the other way around::

    sage: ot2.absorb(ot1)
    O(x^2)

Comparison
==========

The comparison of asymptotic terms with `\leq` is implemented as follows:

- When comparing ``t1 <= t2``, the coercion framework first tries to
  find a common parent for both terms. If this fails, ``False`` is
  returned.

- In case the coerced terms do not have a coefficient in their
  common parent (e.g. :class:`OTerm`), the growth of the two terms
  is compared.

- Otherwise, if the coerced terms have a coefficient (e.g.
  :class:`ExactTerm`), we compare whether ``t1`` has a growth that is
  strictly weaker than the growth of ``t2``. If so, we return
  ``True``. If the terms have equal growth, then we return ``True``
  if and only if the coefficients coincide as well.

  In all other cases, we return ``False``.

Long story short: we consider terms with different coefficients that
have equal growth to be incomparable.


Various
=======

.. TODO::

    - Implementation of more term types (e.g. `L` terms,
      `\Omega` terms, `o` terms, `\Theta` terms).

AUTHORS:

- Benjamin Hackl (2015)
- Daniel Krenn (2015)
- Clemens Heuberger (2016)

ACKNOWLEDGEMENT:

- Benjamin Hackl, Clemens Heuberger and Daniel Krenn are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

- Benjamin Hackl is supported by the Google Summer of Code 2015.


Classes and Methods
===================
"""
# *****************************************************************************
# Copyright (C) 2014--2015 Benjamin Hackl <benjamin.hackl@aau.at>
#               2014--2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# https://www.gnu.org/licenses/
# *****************************************************************************

from sage.rings.big_oh import O
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.factory import UniqueFactory
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.richcmp import richcmp_by_eq_and_lt
from .misc import WithLocals


class ZeroCoefficientError(ValueError):
    pass


def absorption(left, right):
    r"""
    Let one of the two passed terms absorb the other.

    Helper function used by
    :class:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticExpansion`.

    .. NOTE::

        If neither of the terms can absorb the other, an
        :python:`ArithmeticError<library/exceptions.html#exceptions.ArithmeticError>`
        is raised.

        See the :ref:`module description <term_absorption>` for a
        detailed explanation of absorption.

    INPUT:

    - ``left`` -- an asymptotic term.

    - ``right`` -- an asymptotic term.

    OUTPUT:

    An asymptotic term or ``None``.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
        sage: from sage.rings.asymptotic.term_monoid import absorption
        sage: T = TermMonoid('O', GrowthGroup('x^ZZ'), ZZ)
        sage: absorption(T(x^2), T(x^3))
        O(x^3)
        sage: absorption(T(x^3), T(x^2))
        O(x^3)

    ::

        sage: T = TermMonoid('exact', GrowthGroup('x^ZZ'), ZZ)
        sage: absorption(T(x^2), T(x^3))
        Traceback (most recent call last):
        ...
        ArithmeticError: Absorption between x^2 and x^3 is not possible.
    """
    try:
        return left.absorb(right)
    except ArithmeticError:
        try:
            return right.absorb(left)
        except ArithmeticError:
            raise ArithmeticError('Absorption between %s and %s is not possible.' % (left, right))


def can_absorb(left, right):
    r"""
    Return whether one of the two input terms is able to absorb the
    other.

    Helper function used by
    :class:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticExpansion`.

    INPUT:

    - ``left`` -- an asymptotic term.

    - ``right`` -- an asymptotic term.

    OUTPUT:

    A boolean.

    .. NOTE::

        See the :ref:`module description <term_absorption>` for a
        detailed explanation of absorption.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
        sage: from sage.rings.asymptotic.term_monoid import can_absorb
        sage: T = TermMonoid('O', GrowthGroup('x^ZZ'), ZZ)
        sage: can_absorb(T(x^2), T(x^3))
        True
        sage: can_absorb(T(x^3), T(x^2))
        True
    """
    return left.can_absorb(right) or right.can_absorb(left)


class GenericTerm(MultiplicativeGroupElement):
    r"""
    Base class for asymptotic terms. Mainly the structure and
    several properties of asymptotic terms are handled here.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
        sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

        sage: G = GrowthGroup('x^ZZ'); x = G.gen()
        sage: T = GenericTermMonoid(TermMonoid, G, QQ)
        sage: t1 = T(x); t2 = T(x^2); (t1, t2)
        (Generic Term with growth x, Generic Term with growth x^2)
        sage: t1 * t2
        Generic Term with growth x^3
        sage: t1.can_absorb(t2)
        False
        sage: t1.absorb(t2)
        Traceback (most recent call last):
        ...
        ArithmeticError: Generic Term with growth x cannot absorb Generic Term with growth x^2
        sage: t1.can_absorb(t1)
        False
    """

    def __init__(self, parent, growth):
        r"""
        See :class:`GenericTerm` for more information.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = GenericTermMonoid(TermMonoid, G, ZZ)
            sage: T(x^2)
            Generic Term with growth x^2

        ::

            sage: from sage.rings.asymptotic.term_monoid import GenericTerm
            sage: GenericTerm(parent=None, growth=x)
            Traceback (most recent call last):
            ...
            ValueError: The parent must be provided
            sage: GenericTerm(T, GrowthGroup('y^ZZ').gen())
            Traceback (most recent call last):
            ...
            ValueError: y is not in Growth Group x^ZZ.
        """
        if parent is None:
            raise ValueError('The parent must be provided')
        self.growth = parent.growth_group(growth)

        super(GenericTerm, self).__init__(parent=parent)

    def _mul_(self, other):
        r"""
        Multiplication of this term by another.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A :class:`GenericTerm`.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element, as well as ``other``
            are from a common parent.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = GenericTermMonoid(TermMonoid, G, ZZ)
            sage: t1 = T(x); t2 = T(x^2)
            sage: t1, t2
            (Generic Term with growth x, Generic Term with growth x^2)
            sage: t1 * t2
            Generic Term with growth x^3
        """
        return self.parent()(self.growth * other.growth)

    def __invert__(self):
        r"""
        Invert this term.

        OUTPUT:

        A :class:`GenericTerm`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = GenericTermMonoid(TermMonoid, G, QQ)
            sage: ~T(x) # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Inversion of Generic Term with growth x
            not implemented (in this abstract method).

        ::

            sage: t1 = T(x); t2 = T(x^2)
            sage: t1 / t2
            Traceback (most recent call last):
            ...
            NotImplementedError: Inversion of Generic Term with growth x^2
            not implemented (in this abstract method).
        """
        raise NotImplementedError('Inversion of %s not implemented '
                                  '(in this abstract method).' % (self,))

    def __pow__(self, exponent):
        r"""
        Calculate the power of this element to the given ``exponent``.

        INPUT:

        - ``exponent`` -- an element.

        OUTPUT:

        Raise a :python:`NotImplementedError<library/exceptions.html#exceptions.NotImplementedError>`
        since it is an abstract base class.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('z^ZZ')
            sage: t = GenericTermMonoid(TermMonoid, G, ZZ).an_element(); t
            Generic Term with growth z
            sage: t^3  # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Taking powers of Generic Term with growth z
            not implemented (in this abstract method).
        """
        raise NotImplementedError('Taking powers of %s not implemented '
                                  '(in this abstract method).' % (self,))

    def _calculate_pow_test_zero_(self, exponent):
        r"""
        Helper function for :meth:`__pow__` which calculates the power of this
        element to the given ``exponent`` only if zero to this exponent is possible.

        INPUT:

        - ``exponent`` -- an element.

        OUTPUT:

        A term.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('z^ZZ')
            sage: t = GenericTermMonoid(TermMonoid, G, ZZ).an_element(); t
            Generic Term with growth z
            sage: t._calculate_pow_test_zero_(3)
            Generic Term with growth z^3
            sage: t._calculate_pow_test_zero_(-2)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot take Generic Term with growth z to exponent -2.
            > *previous* ZeroDivisionError: rational division by zero

        ::

            sage: TermMonoid('O', G, QQ)('z')._calculate_pow_test_zero_(-1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot take O(z) to exponent -1.
            > *previous* ZeroDivisionError: rational division by zero
        """
        # This assumes `0 = O(g)` for any `g` in the growth group, which
        # is valid in the case of a variable going to `\infty`.
        # Once non-standard asymptotics are supported, this has to be
        # rewritten.
        # See also #19083, comment 64, 27.

        zero = self.parent().coefficient_ring.zero()
        try:
            zero ** exponent
        except (TypeError, ValueError, ZeroDivisionError) as e:
            from .misc import combine_exceptions
            raise combine_exceptions(
                ZeroDivisionError('Cannot take %s to exponent %s.' %
                                  (self, exponent)), e)
        return self._calculate_pow_(exponent)

    def _calculate_pow_(self, exponent, new_coefficient=None):
        r"""
        Helper function for :meth:`__pow__` which calculates the power of this
        element to the given ``exponent``.

        INPUT:

        - ``exponent`` -- an element.

        - ``new_coefficient`` -- if not ``None`` this is passed on to the
          construction of the element (in particular, not taken to any power).

        OUTPUT:

        A term.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('z^ZZ')
            sage: t = GenericTermMonoid(TermMonoid, G, ZZ).an_element(); t
            Generic Term with growth z
            sage: t._calculate_pow_(3)
            Generic Term with growth z^3
            sage: t._calculate_pow_(3, new_coefficient=2)
            Traceback (most recent call last):
            ...
            ValueError: Coefficient 2 is not 1, but Generic Term Monoid z^ZZ with
            (implicit) coefficients in Integer Ring does not support coefficients.
            sage: t._calculate_pow_(-2)
            Generic Term with growth z^(-2)
            sage: t._calculate_pow_(-2, new_coefficient=2)
            Traceback (most recent call last):
            ...
            ValueError: Coefficient 2 is not 1, but Generic Term Monoid z^ZZ with
            (implicit) coefficients in Integer Ring does not support coefficients.
        """
        try:
            g = self.growth ** exponent
        except (ValueError, TypeError, ZeroDivisionError) as e:
            from .misc import combine_exceptions
            raise combine_exceptions(
                ValueError('Cannot take %s to the exponent %s.' % (self, exponent)), e)

        return self.parent()._create_element_in_extension_(g, new_coefficient)

    def can_absorb(self, other):
        r"""
        Check whether this asymptotic term is able to absorb
        the asymptotic term ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            A :class:`GenericTerm` cannot absorb any other term.

            See the :ref:`module description <term_absorption>` for a
            detailed explanation of absorption.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GenericGrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GenericGrowthGroup(ZZ)
            sage: T = GenericTermMonoid(TermMonoid, G, QQ)
            sage: g1 = G(raw_element=21); g2 = G(raw_element=42)
            sage: t1 = T(g1); t2 = T(g2)
            sage: t1.can_absorb(t2)  # indirect doctest
            False
            sage: t2.can_absorb(t1)  # indirect doctest
            False
        """
        return False

    def absorb(self, other, check=True):
        r"""
        Absorb the asymptotic term ``other`` and return the resulting
        asymptotic term.

        INPUT:

        - ``other`` -- an asymptotic term.

        - ``check`` -- a boolean. If ``check`` is ``True`` (default),
          then ``can_absorb`` is called before absorption.

        OUTPUT:

        An asymptotic term or ``None`` if a cancellation occurs. If no
        absorption can be performed, an :python:`ArithmeticError<library/exceptions.html#exceptions.ArithmeticError>`
        is raised.

        .. NOTE::

            Setting ``check`` to ``False`` is meant to be used in
            cases where the respective comparison is done externally
            (in order to avoid duplicate checking).

            For a more detailed explanation of the *absorption* of
            asymptotic terms see
            the :ref:`module description <term_absorption>`.

        EXAMPLES:

        We want to demonstrate in which cases an asymptotic term
        is able to absorb another term, as well as explain the output
        of this operation. We start by defining some parents and
        elements::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G_QQ = GrowthGroup('x^QQ'); x = G_QQ.gen()
            sage: OT = TermMonoid('O', G_QQ, coefficient_ring=ZZ)
            sage: ET = TermMonoid('exact', G_QQ, coefficient_ring=QQ)
            sage: ot1 = OT(x); ot2 = OT(x^2)
            sage: et1 = ET(x, 100); et2 = ET(x^2, 2)
            sage: et3 = ET(x^2, 1); et4 = ET(x^2, -2)

        `O`-Terms are able to absorb other `O`-terms and exact terms
        with weaker or equal growth. ::

            sage: ot1.absorb(ot1)
            O(x)
            sage: ot1.absorb(et1)
            O(x)
            sage: ot1.absorb(et1) is ot1
            True

        :class:`ExactTerm` is able to absorb another
        :class:`ExactTerm` if the terms have the same growth. In this
        case, *absorption* is nothing else than an addition of the
        respective coefficients::

            sage: et2.absorb(et3)
            3*x^2
            sage: et3.absorb(et2)
            3*x^2
            sage: et3.absorb(et4)
            -x^2

        Note that, for technical reasons, the coefficient `0` is not
        allowed, and thus ``None`` is returned if two exact terms
        cancel each other out::

            sage: et2.absorb(et4)
            sage: et4.absorb(et2) is None
            True

        TESTS:

        When disabling the ``check`` flag, absorb might produce
        wrong results::

            sage: et1.absorb(ot2, check=False)
            O(x)
        """
        from sage.structure.element import have_same_parent

        if check:
            if not self.can_absorb(other):
                raise ArithmeticError('%s cannot absorb %s' % (self, other))

        if have_same_parent(self, other):
            return self._absorb_(other)

        from sage.structure.element import get_coercion_model

        return get_coercion_model().bin_op(self, other,
                                           lambda left, right:
                                           left._absorb_(right))

    def _absorb_(self, other):
        r"""
        Let this element absorb ``other``.

        INPUT:

        - ``other`` -- an asymptotic term from the same parent as
          this element.

        OUTPUT:

        An asymptotic term or ``None``.

        .. NOTE::

            This is not implemented for abstract base classes. For
            concrete realizations see, for example, :meth:`OTerm._absorb_`
            or :meth:`ExactTerm._absorb_`.
            Override this in derived class.

        EXAMPLES:

        First, we define some asymptotic terms::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = GenericTermMonoid(TermMonoid, G, QQ)
            sage: t1 = T(x); t2 = T(x^2)

        When it comes to absorption, note that the method
        :meth:`can_absorb` (which is called before absorption takes
        place) does not allow the absorption of generic terms. Thus,
        an :python:`ArithmeticError<library/exceptions.html#exceptions.ArithmeticError>`
        is raised::

            sage: t2.absorb(t1)
            Traceback (most recent call last):
            ...
            ArithmeticError: Generic Term with growth x^2 cannot absorb Generic Term with growth x

        TESTS::

            sage: t2._absorb_(t1)
            Traceback (most recent call last):
            ...
            NotImplementedError: Not implemented in abstract base classes
        """
        raise NotImplementedError('Not implemented in abstract base classes')

    def log_term(self, base=None, locals=None):
        r"""
        Determine the logarithm of this term.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.

        - ``locals`` -- a dictionary which may contain the following keys and values:

          - ``'log'`` -- value: a function. If not used, then the usual
            :class:`log <sage.functions.log.Function_log>` is taken.

        OUTPUT:

        A tuple of terms.

        .. NOTE::

            This abstract method raises a
            :python:`NotImplementedError<library/exceptions.html#exceptions.NotImplementedError>`.
            See :class:`ExactTerm` and :class:`OTerm` for a concrete
            implementation.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: T = GenericTermMonoid(TermMonoid, GrowthGroup('x^ZZ'), QQ)
            sage: T.an_element().log_term()
            Traceback (most recent call last):
            ...
            NotImplementedError: This method is not implemented in
            this abstract base class.

        ::

            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: T = TermWithCoefficientMonoid(TermMonoid, GrowthGroup('x^ZZ'), QQ)
            sage: T.an_element().log_term()
            Traceback (most recent call last):
            ...
            NotImplementedError: This method is not implemented in
            this abstract base class.

        .. SEEALSO::

            :meth:`ExactTerm.log_term`,
            :meth:`OTerm.log_term`.
        """
        raise NotImplementedError('This method is not implemented in this '
                                  'abstract base class.')

    def _log_growth_(self, base=None, locals=None):
        r"""
        Helper function to calculate the logarithm of the growth of this element.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.

        - ``locals`` -- a dictionary which may contain the following keys and values:

          - ``'log'`` -- value: a function. If not used, then the usual
            :class:`log <sage.functions.log.Function_log>` is taken.

        OUTPUT:

        A tuple of terms.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('O', GrowthGroup('x^ZZ * log(x)^ZZ'), QQ)
            sage: T(x^2)._log_growth_()
            (O(log(x)),)
            sage: T(x^1234).log_term()  # indirect doctest
            (O(log(x)),)

        .. SEEALSO::

            :meth:`ExactTerm.log_term`,
            :meth:`OTerm.log_term`.
        """
        return tuple(self.parent()._create_element_in_extension_(g, c)
                     for g, c in
                     self.growth.log_factor(base=base,
                                            locals=locals))

    _richcmp_ = richcmp_by_eq_and_lt("_eq_", "_lt_")

    def _lt_(self, other):
        r"""
        Return whether this generic term grows less than ``other``.

        INPUT:

        - ``other`` -- an asymptotic term

        OUTPUT:

        A boolean.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element, as well as ``other``
            are from the same parent.

            Also, this method **only** compares the growth of the
            input terms!

        EXAMPLES:

        First, we define some asymptotic terms (and their parents)::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: GT = GenericTermMonoid(TermMonoid, G, QQ)
            sage: OT = TermMonoid('O', G, QQ)
            sage: ET_ZZ = TermMonoid('exact', G, ZZ)
            sage: ET_QQ = TermMonoid('exact', G, QQ)
            sage: g1 = GT(x); g2 = GT(x^2); g1, g2
            (Generic Term with growth x, Generic Term with growth x^2)
            sage: o1 = OT(x^-1); o2 = OT(x^3); o1, o2
            (O(x^(-1)), O(x^3))
            sage: t1 = ET_ZZ(x^2, 5); t2 = ET_QQ(x^3, 2/7); t1, t2
            (5*x^2, 2/7*x^3)

        In order for the comparison to work, the terms have to come from
        or coerce into the same parent.

        ::

            sage: g1 <= g2
            True

        If the elements of the common parent do not possess
        coefficients, then only the growth is compared::

            sage: o1 <= o1
            True
            sage: o1 <= o2
            True
            sage: o1 <= t1 and t1 <= o2
            True

        For terms with coefficient (like exact terms), comparison
        works similarly, with the sole exception that terms with
        equal growth are considered incomparable. Thus, `\leq`
        only holds if the coefficients are equal as well::

            sage: t1 <= t2
            True
            sage: ET_ZZ(x, -5) <= ET_ZZ(x, 42)
            False
            sage: ET_ZZ(x, 5) <= ET_ZZ(x, 5)
            True

        ::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: ET = TermMonoid('exact', G, QQ)
            sage: t1 = ET(x, 5); t2 = ET(x^2, 3); t3 = ET(x^2, 42)
            sage: t1 <= t2  # indirect doctest
            True
            sage: t2 <= t1  # indirect doctest
            False
            sage: t2 <= t3  # indirect doctest
            False
            sage: t3 <= t2  # indirect doctest
            False
            sage: t2 <= t2  # indirect doctest
            True

        TESTS::

            sage: ET(x, -2) <= ET(x, 1)
            False

        ::

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = GenericTermMonoid(TermMonoid, G, QQ)
            sage: t1 = T(x^-2); t2 = T(x^5); t1, t2
            (Generic Term with growth x^(-2), Generic Term with growth x^5)
            sage: t1 <= t2  # indirect doctest
            True
            sage: t2 <= t1  # indirect doctest
            False
        """
        return self.growth < other.growth

    def _eq_(self, other):
        r"""
        Return whether this generic term is equal to ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method gets called by the coercion framework, so it
            can be assumed that this asymptotic term is from the
            same parent as ``other``.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import (GenericTermMonoid,
            ....:      ExactTermMonoid, OTermMonoid)
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: GT = GenericTermMonoid(TermMonoid, GrowthGroup('x^ZZ'), QQ)
            sage: ET = ExactTermMonoid(TermMonoid, GrowthGroup('x^ZZ'), ZZ)
            sage: OT = OTermMonoid(TermMonoid, GrowthGroup('x^ZZ'), QQ)
            sage: g = GT.an_element(); e = ET.an_element(); o = OT.an_element()
            sage: g, e, o
            (Generic Term with growth x, x, O(x))
            sage: e == e^2  # indirect doctest
            False
            sage: e == ET(x,1)  # indirect doctest
            True
            sage: o == OT(x^2)  # indirect doctest
            False

        TESTS::

            sage: T = GenericTermMonoid(TermMonoid, GrowthGroup('x^ZZ'), QQ)
            sage: t = T.an_element()
            sage: t == t
            True

        ::

            sage: OT = OTermMonoid(TermMonoid, GrowthGroup('x^ZZ'), QQ)
            sage: t = OT.an_element(); t
            O(x)
            sage: t == OT(x)  # indirect doctest
            True
            sage: t == OT(x^2)  # indirect doctest
            False
        """
        return self.growth == other.growth

    def is_constant(self):
        r"""
        Return whether this term is an (exact) constant.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        .. NOTE::

            Only :class:`ExactTerm` with constant growth (`1`) are
            constant.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = GenericTermMonoid(TermMonoid, GrowthGroup('x^ZZ * log(x)^ZZ'), QQ)
            sage: t = T.an_element(); t
            Generic Term with growth x*log(x)
            sage: t.is_constant()
            False

        ::

            sage: T = TermMonoid('O', GrowthGroup('x^ZZ'), QQ)
            sage: T('x').is_constant()
            False
            sage: T(1).is_constant()
            False
        """
        return False

    def is_exact(self):
        r"""
        Return whether this term is an exact term.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = GenericTermMonoid(TermMonoid, GrowthGroup('x^ZZ * log(x)^ZZ'), QQ)
            sage: T.an_element().is_exact()
            False
        """
        return False

    def is_little_o_of_one(self):
        r"""
        Return whether this generic term is of order `o(1)`.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import (GenericTermMonoid,
            ....:                                                TermWithCoefficientMonoid)
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: T = GenericTermMonoid(TermMonoid, GrowthGroup('x^ZZ'), QQ)
            sage: T.an_element().is_little_o_of_one()
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot check whether Generic Term with growth x is o(1)
            in the abstract base class
            Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field.
            sage: T = TermWithCoefficientMonoid(TermMonoid, GrowthGroup('x^ZZ'), QQ)
            sage: T.an_element().is_little_o_of_one()
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot check whether Term with coefficient 1/2 and growth x
            is o(1) in the abstract base class
            Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field.
        """
        raise NotImplementedError('Cannot check whether %s is o(1) in the '
                                  'abstract base class %s.' % (self, self.parent()))

    def rpow(self, base):
        r"""
        Return the power of ``base`` to this generic term.

        INPUT:

        - ``base`` -- an element or ``'e'``.

        OUTPUT:

        A term.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: T = GenericTermMonoid(TermMonoid, GrowthGroup('x^ZZ * log(x)^ZZ'), QQ)
            sage: T.an_element().rpow('e')
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot take e to the exponent
            Generic Term with growth x*log(x) in the abstract base class
            Generic Term Monoid x^ZZ * log(x)^ZZ with (implicit) coefficients in Rational Field.
        """
        raise NotImplementedError('Cannot take %s to the exponent %s in the '
                                  'abstract base class %s.' % (base, self, self.parent()))

    def _repr_(self):
        r"""
        A representation string for this generic term.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = GenericTermMonoid(TermMonoid, G, QQ)
            sage: T(x)._repr_()
            'Generic Term with growth x'
            sage: T(x^7)._repr_()
            'Generic Term with growth x^7'
        """
        return 'Generic Term with growth ' + repr(self.growth)

    def _latex_(self):
        r"""
        A LaTeX-representation string for this generic term.

        OUTPUT:

        A NotImplementedError is raised.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = GenericTermMonoid(TermMonoid, G, QQ)
            sage: latex(T(x))
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _substitute_(self, rules):
        r"""
        Substitute the given ``rules`` in this generic term.

        INPUT:

        - ``rules`` -- a dictionary.

        OUTPUT:

        Nothing since a
        :python:`TypeError<library/exceptions.html#exceptions.TypeError>`
        is raised.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: t = GenericTermMonoid(TermMonoid, GrowthGroup('x^ZZ'), ZZ).an_element()
            sage: t._substitute_({})
            Traceback (most recent call last):
            ...
            TypeError: Cannot substitute in Generic Term with growth x in
            Generic Term Monoid x^ZZ with (implicit) coefficients in Integer Ring.
            > *previous* TypeError: Cannot substitute in the abstract base class
            Generic Term Monoid x^ZZ with (implicit) coefficients in Integer Ring.
        """
        from .misc import substitute_raise_exception
        substitute_raise_exception(self, TypeError(
            'Cannot substitute in the abstract '
            'base class %s.' % (self.parent(),)))

    def variable_names(self):
        r"""
        Return the names of the variables of this term.

        OUTPUT:

        A tuple of strings.

        EXAMPLES::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('exact', 'QQ^m * m^QQ * log(n)^ZZ', QQ)
            sage: T('4 * 2^m * m^4 * log(n)').variable_names()
            ('m', 'n')
            sage: T('4 * 2^m * m^4').variable_names()
            ('m',)
            sage: T('4 * log(n)').variable_names()
            ('n',)
            sage: T('4 * m^3').variable_names()
            ('m',)
            sage: T('4 * m^0').variable_names()
            ()
        """
        return self.growth.variable_names()

    def _factorial_(self):
        r"""
        Return the factorial of this generic term.

        OUTPUT:

        A term.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: T = GenericTermMonoid(TermMonoid, GrowthGroup('x^QQ'), QQ)
            sage: T.an_element()._factorial_()
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot build the factorial of
            Generic Term with growth x^(1/2).
        """
        raise NotImplementedError(
            'Cannot build the factorial of {}.'.format(self))

    def _singularity_analysis_(self, var, zeta, precision):
        r"""
        Perform singularity analysis on this term.

        INPUT:

        - ``var`` -- a string denoting the variable

        - ``zeta`` -- a number

        - ``precision`` -- an integer

        OUTPUT:

        An asymptotic expansion for `[z^n] f` where `n` is ``var``
        and `f` has this term as a singular expansion
        in `T=\frac{1}{1-\frac{z}{\zeta}}\to \infty` where this element
        is a monomial in `T`.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: t = GenericTermMonoid(TermMonoid, GrowthGroup('x^ZZ'), ZZ).an_element()
            sage: t._singularity_analysis_('n', 2, precision=3)
            Traceback (most recent call last):
            ...
            NotImplementedError: singularity analysis of Generic Term with growth x
            not implemented
        """
        raise NotImplementedError('singularity analysis of {} '
                                  'not implemented '.format(self))


class GenericTermMonoid(UniqueRepresentation, Parent, WithLocals):
    r"""
    Parent for generic asymptotic terms.

    INPUT:

    - ``growth_group`` -- a growth group (i.e. an instance of
      :class:`~sage.rings.asymptotic.growth_group.GenericGrowthGroup`).

    - ``coefficient_ring`` -- a ring which contains the (maybe implicit)
      coefficients of the elements.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. It has to be a subcategory
      of ``Join of Category of Monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    In this class the base
    structure for asymptotic term monoids will be handled. These
    monoids are the parents of asymptotic terms (for example, see
    :class:`GenericTerm` or :class:`OTerm`). Basically, asymptotic
    terms consist of a ``growth`` (which is an asymptotic growth
    group element, for example
    :class:`~sage.rings.asymptotic.growth_group.MonomialGrowthElement`);
    additional structure and properties are added by the classes inherited
    from :class:`GenericTermMonoid`.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
        sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

        sage: G_x = GrowthGroup('x^ZZ'); x = G_x.gen()
        sage: G_y = GrowthGroup('y^QQ'); y = G_y.gen()
        sage: T_x_ZZ = GenericTermMonoid(TermMonoid, G_x, QQ)
        sage: T_y_QQ = GenericTermMonoid(TermMonoid, G_y, QQ)
        sage: T_x_ZZ
        Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
        sage: T_y_QQ
        Generic Term Monoid y^QQ with (implicit) coefficients in Rational Field
    """

    # enable the category framework for elements
    Element = GenericTerm

    @staticmethod
    def __classcall__(cls, term_monoid_factory,
                      growth_group, coefficient_ring,
                      category=None):
        r"""
        Normalize the input in order to ensure a unique
        representation of the parent.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('x^ZZ')
            sage: T = GenericTermMonoid(TermMonoid, G, QQ)
            sage: T._underlying_class()(TermMonoid, G, QQ) is T
            True

        ::

            sage: GenericTermMonoid(TermMonoid, None, ZZ)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: No growth group specified.
            sage: GenericTermMonoid(TermMonoid, int, ZZ)  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: <... 'int'> is not a valid growth group.
            sage: GenericTermMonoid(TermMonoid, G, None)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: No coefficient ring specified.
            sage: GenericTermMonoid(TermMonoid, G, int)  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: <... 'int'> is not a valid coefficient ring.
        """
        if term_monoid_factory is None:
            raise ValueError('No term monoid factory specified.')
        if not isinstance(term_monoid_factory, TermMonoidFactory):
            raise TypeError('%s is not a valid term monoid factory.' % (term_monoid_factory,))

        if growth_group is None:
            raise ValueError('No growth group specified.')
        if not isinstance(growth_group, Parent):
            raise TypeError('%s is not a valid growth group.' % (growth_group,))

        if coefficient_ring is None:
            raise ValueError('No coefficient ring specified.')
        if not isinstance(coefficient_ring, Parent):
            raise TypeError('%s is not a valid coefficient ring.' % (coefficient_ring,))

        if category is None:
            from sage.categories.monoids import Monoids
            from sage.categories.posets import Posets
            category = Monoids() & Posets()

        return super(GenericTermMonoid, cls).__classcall__(
            cls, term_monoid_factory, growth_group, coefficient_ring, category)

    def __init__(self, term_monoid_factory, growth_group, coefficient_ring, category):
        r"""
        See :class:`GenericTermMonoid` for more information.

        EXAMPLES::

            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid, TermWithCoefficientMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G_x = GrowthGroup('x^ZZ')
            sage: T_x = GenericTermMonoid(TermMonoid, G_x, QQ); T_x
            Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
            sage: T_x.growth_group
            Growth Group x^ZZ
            sage: G_y = GrowthGroup('y^QQ')
            sage: T_y = GenericTermMonoid(TermMonoid, G_y, QQ); T_y
            Generic Term Monoid y^QQ with (implicit) coefficients in Rational Field
            sage: T_x is T_y
            False

        ::

            sage: GenericTermMonoid(TermMonoid, None, None)
            Traceback (most recent call last):
            ...
            ValueError: No growth group specified.

        ::

            sage: G = GrowthGroup('x^ZZ')
            sage: T_ZZ = TermWithCoefficientMonoid(TermMonoid, G, ZZ); T_ZZ
            Generic Term Monoid x^ZZ with (implicit) coefficients in Integer Ring
            sage: T_QQ = TermWithCoefficientMonoid(TermMonoid, G, QQ); T_QQ
            Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
            sage: T_QQ.category()
            Join of Category of monoids and Category of posets
        """
        self._term_monoid_factory_ = term_monoid_factory
        self._growth_group_ = growth_group
        self._coefficient_ring_ = coefficient_ring
        super(GenericTermMonoid, self).__init__(category=category)

    @property
    def term_monoid_factory(self):
        r"""
        The term monoid factory capable of creating this term monoid.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory
            sage: DefaultTermMonoidFactory('exact', GrowthGroup('x^ZZ'), ZZ).term_monoid_factory
            Term Monoid Factory 'sage.rings.asymptotic.term_monoid.DefaultTermMonoidFactory'

            sage: from sage.rings.asymptotic.term_monoid import TermMonoidFactory
            sage: TermMonoid = TermMonoidFactory('__main__.TermMonoid')

            sage: TermMonoid('exact', GrowthGroup('x^ZZ'), ZZ).term_monoid_factory
            Term Monoid Factory '__main__.TermMonoid'
        """
        return self._term_monoid_factory_

    def term_monoid(self, type):
        r"""
        Return the term monoid of specified ``type``.

        INPUT:

        - ``type`` -- 'O' or 'exact', or an instance of an existing
          term monoid.
          See :class:`~sage.rings.asymptotic.term_monoid.TermMonoidFactory`
          for more details.

        OUTPUT:

        A term monoid object derived from
        :class:`~sage.rings.asymptotic.term_monoid.GenericTermMonoid`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: E = TermMonoid('exact', GrowthGroup('x^ZZ'), ZZ); E
            Exact Term Monoid x^ZZ with coefficients in Integer Ring
            sage: E.term_monoid('O')
            O-Term Monoid x^ZZ with implicit coefficients in Integer Ring

        TESTS::

            sage: E.term_monoid('exact') is E
            True
        """
        TermMonoid = self.term_monoid_factory
        return TermMonoid(type,
                          growth_group=self.growth_group,
                          coefficient_ring=self.coefficient_ring)

    @property
    def growth_group(self):
        r"""
        The growth group underlying this term monoid.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: TermMonoid('exact', GrowthGroup('x^ZZ'), ZZ).growth_group
            Growth Group x^ZZ
        """
        return self._growth_group_

    @property
    def coefficient_ring(self):
        r"""
        The coefficient ring of this term monoid, i.e. the ring where
        the coefficients are from.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: GenericTermMonoid(TermMonoid, GrowthGroup('x^ZZ'), ZZ).coefficient_ring
            Integer Ring
        """
        return self._coefficient_ring_

    def change_parameter(self, growth_group=None, coefficient_ring=None):
        r"""
        Return a term monoid with a change in one or more of the
        given parameters.

        INPUT:

        - ``growth_group`` -- (default: ``None``) the new growth group.

        - ``coefficient_ring`` -- (default: ``None``) the new coefficient ring.

        OUTPUT:

        A term monoid.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: E = TermMonoid('exact', GrowthGroup('n^ZZ'), ZZ)
            sage: E.change_parameter(coefficient_ring=QQ)
            Exact Term Monoid n^ZZ with coefficients in Rational Field
            sage: E.change_parameter(growth_group=GrowthGroup('n^QQ'))
            Exact Term Monoid n^QQ with coefficients in Integer Ring

        TESTS::

            sage: E.change_parameter() is E
            True
            sage: E.change_parameter(growth_group=None) is E
            True
            sage: E.change_parameter(coefficient_ring=None) is E
            True
            sage: E.change_parameter(growth_group=None, coefficient_ring=None) is E
            True
        """
        if growth_group is None:
            growth_group = self.growth_group
        if coefficient_ring is None:
            coefficient_ring = self.coefficient_ring
        if self.growth_group is growth_group and \
                self.coefficient_ring is coefficient_ring:
            return self
        return self.term_monoid_factory(self, growth_group, coefficient_ring)

    def _repr_(self):
        r"""
        A representation string for this generic term monoid.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import (GenericGrowthGroup, GrowthGroup)
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: GenericTermMonoid(TermMonoid, GenericGrowthGroup(ZZ), QQ)._repr_()
            'Generic Term Monoid Generic(ZZ) with (implicit) coefficients in Rational Field'
            sage: GenericTermMonoid(TermMonoid, GrowthGroup('x^ZZ'), QQ)._repr_()
            'Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field'
        """
        return 'Generic Term Monoid %s with (implicit) coefficients in %s' % \
            (self.growth_group._repr_short_(), self.coefficient_ring)

    def _coerce_map_from_(self, S):
        r"""
        Return whether ``S`` coerces into this term monoid.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            Another generic term monoid ``S`` coerces into this term
            monoid if and only if both, the growth group of ``S`` coerces
            into the growth group of this term monoid and the coefficient
            ring of ``S`` coerces into the coefficient ring of this term
            monoid.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: T_ZZ = GenericTermMonoid(TermMonoid, G_ZZ, QQ); T_ZZ
            Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
            sage: G_QQ = GrowthGroup('x^QQ')
            sage: T_QQ = GenericTermMonoid(TermMonoid, G_QQ, QQ); T_QQ
            Generic Term Monoid x^QQ with (implicit) coefficients in Rational Field
            sage: T_QQ.has_coerce_map_from(T_ZZ)  # indirect doctest
            True
            sage: T_QQ_ZZ = GenericTermMonoid(TermMonoid, G_QQ, ZZ); T_QQ_ZZ
            Generic Term Monoid x^QQ with (implicit) coefficients in Integer Ring
            sage: T_QQ.has_coerce_map_from(T_QQ_ZZ)
            True
            sage: T_QQ_ZZ.has_coerce_map_from(T_QQ)
            False

        ::

            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: TC_ZZ = TermWithCoefficientMonoid(TermMonoid, G_ZZ, ZZ); TC_ZZ
            Generic Term Monoid x^ZZ with (implicit) coefficients in Integer Ring
            sage: TC_QQ = TermWithCoefficientMonoid(TermMonoid, G_QQ, QQ); TC_QQ
            Generic Term Monoid x^QQ with (implicit) coefficients in Rational Field
            sage: TC_QQ.has_coerce_map_from(TC_ZZ)  # indirect doctest
            True
            sage: TC_ZZ.has_coerce_map_from(TC_QQ)  # indirect doctest
            False
        """
        if isinstance(S, self.__class__):
            if self.growth_group.has_coerce_map_from(S.growth_group) and \
                    self.coefficient_ring.has_coerce_map_from(S.coefficient_ring):
                return True

    def _element_constructor_(self, data, coefficient=None):
        r"""
        Convert the given object to this term monoid.

        INPUT:

        - ``data`` -- a growth element or an object representing the
          element to be initialized.

        - ``coefficient`` -- (default: ``None``)
          an element of the coefficient ring.

        OUTPUT:

        An element of this term monoid.

        EXAMPLES::

            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: G_QQ = GrowthGroup('x^QQ')
            sage: T_ZZ = GenericTermMonoid(TermMonoid, G_ZZ, QQ)
            sage: T_QQ = GenericTermMonoid(TermMonoid, G_QQ, QQ)
            sage: term1 = T_ZZ(G_ZZ.gen())
            sage: term2 = T_QQ(G_QQ.gen()^2)

        In order for two terms to be compared, a coercion into
        a common parent has to be found::

            sage: term1.parent()
            Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
            sage: term2.parent()
            Generic Term Monoid x^QQ with (implicit) coefficients in Rational Field
            sage: term1 <= term2
            True

        In this case, this works because ``T_ZZ``, the parent of
        ``term1``, coerces into ``T_QQ``::

            sage: T_QQ.coerce(term1)  # coercion does not fail
            Generic Term with growth x

        The conversion of growth elements also works for the creation
        of terms::

            sage: x = SR('x'); x.parent()
            Symbolic Ring
            sage: T_ZZ(x^42)
            Generic Term with growth x^42
            sage: x = PolynomialRing(ZZ, 'x').gen(); x.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: T_ZZ(x^10)
            Generic Term with growth x^10
            sage: T_ZZ(10 * x^2)
            Traceback (most recent call last):
            ...
            ValueError: 10*x^2 is not in Generic Term Monoid x^ZZ
            with (implicit) coefficients in Rational Field.
            > *previous* ValueError: Factor 10*x^2 of 10*x^2 is neither a
            coefficient (in Rational Field) nor growth (in Growth Group x^ZZ).

        ::

            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: G = GrowthGroup('x^ZZ')
            sage: T = TermWithCoefficientMonoid(TermMonoid, G, ZZ)
            sage: t1 = T(x^2, 5); t1  # indirect doctest
            Term with coefficient 5 and growth x^2

        TESTS::

            sage: O_ZZ = TermMonoid('O', G_ZZ, QQ)
            sage: O_ZZ(x^11)
            O(x^11)

        ::

            sage: T(G.gen()^10)
            Term with coefficient 1 and growth x^10
            sage: T(G.gen()^10, coefficient=10)
            Term with coefficient 10 and growth x^10
            sage: T(x^123)
            Term with coefficient 1 and growth x^123

        ::

            sage: T(x)
            Term with coefficient 1 and growth x

        ::

            sage: G_log = GrowthGroup('log(x)^ZZ')
            sage: T_log = TermWithCoefficientMonoid(TermMonoid, G_log, ZZ)
            sage: T_log(log(x))
            Term with coefficient 1 and growth log(x)

        """
        if isinstance(data, self.element_class) and data.parent() == self:
            return data
        elif isinstance(data, TermWithCoefficient):
            return self._create_element_(data.growth, data.coefficient)
        elif isinstance(data, GenericTerm):
            return self._create_element_(data.growth, None)
        elif isinstance(data, int) and data == 0:
            raise ValueError('No input specified. Cannot continue '
                             'creating an element of %s.' % (self,))

        from .misc import combine_exceptions
        if coefficient is not None:
            try:
                data = self.growth_group(data)
            except (ValueError, TypeError) as e:
                raise combine_exceptions(
                    ValueError('Growth %s is not in %s.' % (data, self)), e)
            return self._create_element_(data, coefficient)

        try:
            growth, coefficient = self._split_growth_and_coefficient_(data)
        except ValueError as e:
            raise combine_exceptions(
                ValueError('%s is not in %s.' % (data, self)), e)

        return self._create_element_(growth, coefficient)

    def _create_element_(self, growth, coefficient):
        r"""
        Helper method which creates an element by using the ``element_class``.

        INPUT:

        - ``growth`` -- a growth element.

        - ``coefficient`` -- an element of the coefficient ring.

        OUTPUT:

        A term.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: T_ZZ = GenericTermMonoid(TermMonoid, G_ZZ, QQ)
            sage: T_ZZ(G_ZZ.gen())  # indirect doctest
            Generic Term with growth x
        """
        if coefficient is not None and coefficient != self.coefficient_ring.one():
            raise ValueError('Coefficient %s is not 1, but %s does not '
                             'support coefficients.' % (coefficient, self))
        return self.element_class(self, growth)

    def _create_element_in_extension_(self, growth, coefficient):
        r"""
        Create an element in an extension of this term monoid which
        is chosen according to the input.

        INPUT:

        - ``growth`` and ``coefficient`` -- the element data.

        OUTPUT:

        An element.

        EXAMPLES::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('z^ZZ')
            sage: T = TermMonoid('exact', G, ZZ)
            sage: T._create_element_in_extension_(G.an_element(), 3)
            3*z
            sage: T._create_element_in_extension_(G.an_element(), 3/2).parent()
            Exact Term Monoid z^ZZ with coefficients in Rational Field
        """
        if (growth.parent() is self.growth_group) and \
           (coefficient is None or coefficient.parent() is self.coefficient_ring):
            parent = self
        else:
            parent = self._underlying_class()(self.term_monoid_factory,
                                              growth.parent(),
                                              coefficient.parent()
                                              if coefficient is not None
                                              else self.coefficient_ring,
                                              category=self.category())
        return parent(growth, coefficient)

    def _split_growth_and_coefficient_(self, data):
        r"""
        Split given ``data`` into a growth element and a coefficient.

        INPUT:

        ``data`` -- an element.

        OUTPUT:

        A pair ``(growth, coefficient``).

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ')
            sage: T = TermMonoid('exact', G, QQ)
            sage: T._split_growth_and_coefficient_('2*x^3')
            (x^3, 2)

        ::

            sage: T._split_growth_and_coefficient_('2.7 * x^3')
            Traceback (most recent call last):
            ...
            ValueError: Factor 2.7 of 2.7 * x^3 is neither a coefficient
            (in Rational Field) nor growth (in Growth Group x^ZZ).

        ::

            sage: G = GrowthGroup('QQ^x * x^ZZ * log(x)^ZZ')
            sage: T = TermMonoid('exact', G, QQ)
            sage: T._split_growth_and_coefficient_('3/4 * 2^x * log(x)')
            (2^x*log(x), 3/4)
            sage: T._split_growth_and_coefficient_('3 * x^2 * 4 * log(x) * x')
            (x^3*log(x), 12)
            sage: var('x')
            x
            sage: T._split_growth_and_coefficient_(log(x)^5 * x^2 * 4)
            (x^2*log(x)^5, 4)

        ::

            sage: T = TermMonoid('exact', G, SR)
            sage: T._split_growth_and_coefficient_(log(x)^5 * x^2 * 4)
            (x^2*log(x)^5, 4)
            sage: var('y')
            y
            sage: T._split_growth_and_coefficient_(2^x * y * 4)
            (2^x, 4*y)
        """
        factors = self._get_factors_(data)

        growth_group = self.growth_group
        coefficient_ring = self.coefficient_ring

        coefficients = []
        growths = []
        for f in factors:
            try:
                growths.append(growth_group(f))
                continue
            except (ValueError, TypeError):
                pass

            try:
                coefficients.append(coefficient_ring(f))
                continue
            except (ValueError, TypeError):
                pass

            raise ValueError('Factor %s of %s is neither a coefficient (in %s) '
                             'nor growth (in %s).' %
                             (f, data, coefficient_ring, growth_group))

        from sage.misc.misc_c import prod
        growth = prod(growths) if growths else growth_group.one()
        coefficient = prod(coefficients) if coefficients else coefficient_ring.one()
        return (growth, coefficient)

    def _get_factors_(self, data):
        r"""
        Split given ``data`` into separate factors.

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        A tuple.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: T_ZZ = TermMonoid('exact', G_ZZ, QQ)
            sage: T_ZZ._get_factors_(x^2 * log(x))
            (x^2, log(x))
        """
        if isinstance(data, str):
            from .misc import split_str_by_op
            return split_str_by_op(data, '*')

        try:
            P = data.parent()
        except AttributeError:
            return (data,)

        from sage.symbolic.ring import SymbolicRing
        if isinstance(P, SymbolicRing):
            from sage.symbolic.operators import mul_vararg
            if data.operator() == mul_vararg:
                return tuple(data.operands())

        return (data,)

    def _an_element_(self):
        r"""
        Return an element of this term monoid.

        INPUT:

        Nothing.

        OUTPUT:

        An element of this term monoid.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: G = GrowthGroup('x^ZZ')
            sage: TermMonoid('O', G, QQ).an_element()  # indirect doctest
            O(x)
            sage: GenericTermMonoid(TermMonoid, G, QQ).an_element()  # indirect doctest
            Generic Term with growth x
        """
        return self(self.growth_group.an_element())

    def some_elements(self):
        r"""
        Return some elements of this term monoid.

        See :class:`TestSuite` for a typical use case.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ')
            sage: tuple(TermMonoid('O', G, QQ).some_elements())
            (O(1), O(x), O(x^(-1)), O(x^2), O(x^(-2)), O(x^3), ...)
        """
        return iter(self(g) for g in self.growth_group.some_elements())

    def le(self, left, right):
        r"""
        Return whether the term ``left`` is at most (less than or equal
        to) the term ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = GenericTermMonoid(TermMonoid, G, QQ)
            sage: t1 = T(x); t2 = T(x^2)
            sage: T.le(t1, t2)
            True
        """
        return self(left) <= self(right)


class OTerm(GenericTerm):
    r"""
    Class for an asymptotic term representing an `O`-term with
    specified growth. For the mathematical properties of `O`-terms
    see :wikipedia:`Big_O_Notation`.

    `O`-terms can *absorb* terms of weaker or equal growth.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- a growth element.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: from sage.rings.asymptotic.term_monoid import OTermMonoid
        sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

        sage: G = GrowthGroup('x^ZZ'); x = G.gen()
        sage: OT = OTermMonoid(TermMonoid, G, QQ)
        sage: t1 = OT(x^-7); t2 = OT(x^5); t3 = OT(x^42)
        sage: t1, t2, t3
        (O(x^(-7)), O(x^5), O(x^42))
        sage: t1.can_absorb(t2)
        False
        sage: t2.can_absorb(t1)
        True
        sage: t2.absorb(t1)
        O(x^5)
        sage: t1 <= t2 and t2 <= t3
        True
        sage: t3 <= t1
        False

    The conversion of growth elements also works for the creation
    of `O`-terms::

        sage: x = SR('x'); x.parent()
        Symbolic Ring
        sage: OT(x^17)
        O(x^17)
    """

    def _repr_(self, latex=False):
        r"""
        A representation string for this `O`-term.

        INPUT:

        - ``latex`` -- (default: ``False``) a boolean. If set, then
          LaTeX-output is returned.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: OT = TermMonoid('O', G, QQ)
            sage: OT(x)._repr_()
            'O(x)'
            sage: OT(x^2)._repr_()
            'O(x^2)'
            sage: OT(x^3)
            O(x^3)
        """
        if latex:
            from sage.misc.latex import latex as latex_repr
            f = latex_repr
        else:
            f = repr

        if latex:
            s = r'O\!\left({g}\right)'
        else:
            s = 'O({g})'
        return s.format(g=f(self.growth))

    def _latex_(self):
        r"""
        A LaTeX-representation string for this `O`-term.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: OT = TermMonoid('O', G, QQ)
            sage: latex(OT(x))
            O\!\left(x\right)
            sage: latex(OT(x^2))
            O\!\left(x^{2}\right)
            sage: latex(OT(x^3))
            O\!\left(x^{3}\right)
        """
        return self._repr_(latex=True)

    def __invert__(self):
        r"""
        Invert this term.

        OUTPUT:

        A :class:`ZeroDivisionError` since `O`-terms cannot be inverted.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = TermMonoid('O', G, QQ)
            sage: ~T(x) # indirect doctest
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot invert O(x).
        """
        raise ZeroDivisionError('Cannot invert %s.' % (self,))

    def __pow__(self, exponent):
        r"""
        Calculate the power of this :class:`OTerm` to the given ``exponent``.

        INPUT:

        - ``exponent`` -- an element.

        OUTPUT:

        An :class:`OTerm`.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('z^ZZ')
            sage: t = TermMonoid('O', G, ZZ).an_element(); t
            O(z)
            sage: t^3  # indirect doctest
            O(z^3)
            sage: t^(1/2)  # indirect doctest
            O(z^(1/2))
            sage: t^(-1)  # indirect doctest
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot take O(z) to exponent -1.
            > *previous* ZeroDivisionError: rational division by zero
        """
        from .misc import strip_symbolic
        return self._calculate_pow_test_zero_(strip_symbolic(exponent))

    def can_absorb(self, other):
        r"""
        Check whether this `O`-term can absorb ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            An :class:`OTerm` can absorb any other asymptotic term
            with weaker or equal growth.

            See the :ref:`module description <term_absorption>` for a
            detailed explanation of absorption.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: OT = TermMonoid('O', GrowthGroup('x^ZZ'), QQ)
            sage: t1 = OT(x^21); t2 = OT(x^42)
            sage: t1.can_absorb(t2)
            False
            sage: t2.can_absorb(t1)
            True
        """
        return other <= self

    def _absorb_(self, other):
        r"""
        Let this `O`-term absorb another `O`-term ``other``.

        INPUT:

        - ``other`` -- an asymptotic `O`-term.

        OUTPUT:

        An asymptotic `O`-term.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element and ``other``
            have the same parent.

            Also, observe that the result of a "dominant" `O`-term
            absorbing another `O`-term always is the "dominant"
            `O`-term again.

            See the :ref:`module description <term_absorption>` for a
            detailed explanation on absorption.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: OT = TermMonoid('O', G, QQ)
            sage: ot1 = OT(x); ot2 = OT(x^2)
            sage: ot1.absorb(ot1)
            O(x)
            sage: ot2.absorb(ot1)
            O(x^2)
            sage: ot1.absorb(ot2)
            Traceback (most recent call last):
            ...
            ArithmeticError: O(x) cannot absorb O(x^2)
        """
        return self

    def log_term(self, base=None, locals=None):
        r"""
        Determine the logarithm of this O-term.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.

        - ``locals`` -- a dictionary which may contain the following keys and values:

          - ``'log'`` -- value: a function. If not used, then the usual
            :class:`log <sage.functions.log.Function_log>` is taken.

        OUTPUT:

        A tuple of terms.

        .. NOTE::

            This method returns a tuple with the summands that come from
            applying the rule `\log(x\cdot y) = \log(x) + \log(y)`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('O', GrowthGroup('x^ZZ * log(x)^ZZ'), QQ)
            sage: T(x^2).log_term()
            (O(log(x)),)
            sage: T(x^1234).log_term()
            (O(log(x)),)

        ::

            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: T = TermMonoid('O', GrowthGroup('x^ZZ * log(x)^ZZ * y^ZZ * log(y)^ZZ'), QQ)
            sage: T('x * y').log_term()
            (O(log(x)), O(log(y)))

        .. SEEALSO::

            :meth:`ExactTerm.log_term`.
        """
        return self._log_growth_(base=base, locals=locals)

    def is_little_o_of_one(self):
        r"""
        Return whether this O-term is of order `o(1)`.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('O', GrowthGroup('x^ZZ'), QQ)
            sage: T(x).is_little_o_of_one()
            False
            sage: T(1).is_little_o_of_one()
            False
            sage: T(x^(-1)).is_little_o_of_one()
            True

        ::

            sage: T = TermMonoid('O', GrowthGroup('x^ZZ * y^ZZ'), QQ)
            sage: T('x * y^(-1)').is_little_o_of_one()
            False
            sage: T('x^(-1) * y').is_little_o_of_one()
            False
            sage: T('x^(-2) * y^(-3)').is_little_o_of_one()
            True

        ::

            sage: T = TermMonoid('O', GrowthGroup('x^QQ * log(x)^QQ'), QQ)
            sage: T('x * log(x)^2').is_little_o_of_one()
            False
            sage: T('x^2 * log(x)^(-1234)').is_little_o_of_one()
            False
            sage: T('x^(-1) * log(x)^4242').is_little_o_of_one()
            True
            sage: T('x^(-1/100) * log(x)^(1000/7)').is_little_o_of_one()
            True
        """
        return self.growth.is_lt_one()

    def rpow(self, base):
        r"""
        Return the power of ``base`` to this O-term.

        INPUT:

        - ``base`` -- an element or ``'e'``.

        OUTPUT:

        A term.

        .. NOTE::

            For :class:`OTerm`, the powers can only be
            constructed for exponents `O(1)` or if ``base`` is `1`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('O', GrowthGroup('x^ZZ * log(x)^ZZ'), QQ)
            sage: T(1).rpow('e')
            O(1)
            sage: T(1).rpow(2)
            O(1)

        ::

            sage: T.an_element().rpow(1)
            1
            sage: T('x^2').rpow(1)
            1

        ::

            sage: T.an_element().rpow('e')
            Traceback (most recent call last):
            ...
            ValueError: Cannot take e to the exponent O(x*log(x)) in
            O-Term Monoid x^ZZ * log(x)^ZZ with implicit coefficients in Rational Field
            sage: T('log(x)').rpow('e')
            Traceback (most recent call last):
            ...
            ValueError: Cannot take e to the exponent O(log(x)) in
            O-Term Monoid x^ZZ * log(x)^ZZ with implicit coefficients in Rational Field
        """
        if self.is_one() and base != 0:
            return self
        if base == 1:
            P = self.parent()
            return P.term_monoid_factory('exact',
                                         P.growth_group,
                                         P.coefficient_ring).one()
        raise ValueError('Cannot take %s to the exponent %s in %s' %
                         (base, self, self.parent()))

    def _substitute_(self, rules):
        r"""
        Substitute the given ``rules`` in this O-Term.

        INPUT:

        - ``rules`` -- a dictionary.

        OUTPUT:

        An object.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('O', GrowthGroup('x^ZZ'), ZZ)
            sage: t = T.an_element(); t
            O(x)
            sage: t._substitute_({'x': SR.var('z')})
            Order(z)
            sage: t._substitute_({'x': SR.var('z'), 'O': function('Oh')})
            Oh(z)
            sage: u = AsymptoticRing('x^ZZ', ZZ)('2*x'); u
            2*x
            sage: t._substitute_({'x': u})
            O(x)
            sage: T(1/x)._substitute_({'x': 0})
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot substitute in O(x^(-1)) in
            O-Term Monoid x^ZZ with implicit coefficients in Integer Ring.
            > *previous* ZeroDivisionError: Cannot substitute in x^(-1) in
            Growth Group x^ZZ.
            >> *previous* ZeroDivisionError: rational division by zero
            sage: t._substitute_({'x': 3})
            O(3)
            sage: t._substitute_({'x': 'null'})
            Traceback (most recent call last):
            ...
            TypeError: Cannot substitute in O(x) in O-Term Monoid x^ZZ with implicit coefficients in Integer Ring.
            > *previous* TypeError: Cannot substitute in x in Growth Group x^ZZ.
            >> *previous* TypeError: unsupported operand type(s) for ** or pow(): 'str' and 'int'
        """
        try:
            g = self.growth._substitute_(rules)
        except (ArithmeticError, TypeError, ValueError) as e:
            from .misc import substitute_raise_exception
            substitute_raise_exception(self, e)

        try:
            return rules['O'](g)
        except KeyError:
            pass

        try:
            P = g.parent()
        except AttributeError:
            pass
        else:
            from .asymptotic_ring import AsymptoticRing
            from sage.symbolic.ring import SymbolicRing

            if isinstance(P, AsymptoticRing):
                return g.O()

            elif isinstance(P, SymbolicRing):
                return g.Order()

        try:
            return O(g)
        except (ArithmeticError, TypeError, ValueError) as e:
            from .misc import substitute_raise_exception
            substitute_raise_exception(self, e)

    def _factorial_(self):
        r"""
        Return the factorial of this O-term if it is constant
        (i.e., has growth `1`).

        OUTPUT:

        A term.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('O', 'z^QQ', QQ)
            sage: T(1)._factorial_()
            O(1)
            sage: T('z^(3/2)')._factorial_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot build the factorial of O(z^(3/2))
            since it has growth != 1.
        """
        if not self.growth.is_one():
            raise ValueError(
                'Cannot build the factorial of {} since it has growth '
                '!= 1.'.format(self))

        return self

    def _singularity_analysis_(self, var, zeta, precision):
        r"""
        Perform singularity analysis on this O-term.

        INPUT:

        - ``var`` -- a string denoting the variable

        - ``zeta`` -- a number

        - ``precision`` -- an integer

        OUTPUT:

        An asymptotic expansion for `[z^n] f` where `n` is ``var``
        and `f` has this term as a singular expansion
        in `T=\frac{1}{1-\frac{z}{\zeta}}\to \infty` where this element
        is a monomial in `T`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('O', GrowthGroup('x^ZZ'), ZZ)
            sage: T('x^1')._singularity_analysis_('n', 2, precision=3)
            O((1/2)^n)
            sage: T('x^(-1)')._singularity_analysis_('n', 2, precision=3)
            O((1/2)^n*n^(-2))
        """
        return self.growth._singularity_analysis_(
            var=var, zeta=zeta, precision=0)


class OTermMonoid(GenericTermMonoid):
    r"""
    Parent for asymptotic big `O`-terms.

    INPUT:

    - ``growth_group`` -- a growth group.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. It has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: from sage.rings.asymptotic.term_monoid import OTermMonoid
        sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
        sage: G_x_ZZ = GrowthGroup('x^ZZ')
        sage: G_y_QQ = GrowthGroup('y^QQ')
        sage: OT_x_ZZ = OTermMonoid(TermMonoid, G_x_ZZ, QQ); OT_x_ZZ
        O-Term Monoid x^ZZ with implicit coefficients in Rational Field
        sage: OT_y_QQ = OTermMonoid(TermMonoid, G_y_QQ, QQ); OT_y_QQ
        O-Term Monoid y^QQ with implicit coefficients in Rational Field

    `O`-term monoids can also be created by using the
    :class:`term factory <TermMonoidFactory>`::

        sage: TermMonoid('O', G_x_ZZ, QQ) is OT_x_ZZ
        True
        sage: TermMonoid('O', GrowthGroup('x^QQ'), QQ)
        O-Term Monoid x^QQ with implicit coefficients in Rational Field
    """

    # enable the category framework for elements
    Element = OTerm

    def _create_element_(self, growth, coefficient):
        r"""
        Helper method which creates an element by using the ``element_class``.

        INPUT:

        - ``growth`` -- a growth element.

        - ``coefficient`` -- an element of the coefficient ring (will be
          ignored since we create an O-Term).

        OUTPUT:

        An O-term.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ')
            sage: T = TermMonoid('O', G, QQ)
            sage: T(G.gen())  # indirect doctest
            O(x)
            sage: T(G.gen(), SR.var('y'))  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot create O(x) since given coefficient y
            is not valid in O-Term Monoid x^ZZ with implicit coefficients in
            Rational Field.
            > *previous* TypeError: unable to convert y to a rational
        """
        if coefficient is not None:
            try:
                self.coefficient_ring(coefficient)
            except (TypeError, ValueError) as e:
                from .misc import combine_exceptions
                raise combine_exceptions(
                    ValueError('Cannot create O(%s) since given coefficient %s '
                               'is not valid in %s.' %
                               (growth, coefficient, self)), e)
        return self.element_class(self, growth)

    def _coerce_map_from_(self, S):
        r"""
        Return whether ``S`` coerces into this term monoid.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        ``True`` or ``None``.

        .. NOTE::

            Another term monoid ``S`` coerces into this term monoid
            if ``S`` is an instance of one of the following classes:

            - :class:`OTermMonoid`

            - :class:`ExactTermMonoid`

            Additionally, the growth group underlying ``S`` has to
            coerce into the growth group of this term monoid.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G_ZZ = GrowthGroup('x^ZZ'); x_ZZ = G_ZZ.gen()
            sage: G_QQ = GrowthGroup('x^QQ'); x_QQ = G_QQ.gen()
            sage: OT_ZZ = TermMonoid('O', G_ZZ, QQ)
            sage: OT_QQ = TermMonoid('O', G_QQ, QQ)
            sage: ET = TermMonoid('exact', G_ZZ, ZZ)

        Now, the :class:`OTermMonoid` whose growth group is over the
        integer ring has to coerce into the :class:`OTermMonoid` with
        the growth group over the rational field, and the
        :class:`ExactTermMonoid` also has to coerce in each of the
        given :class:`OTermMonoid`::

            sage: OT_QQ.has_coerce_map_from(OT_ZZ)  # indirect doctest
            True
            sage: OT_QQ.has_coerce_map_from(ET)  # indirect doctest
            True
            sage: ET.has_coerce_map_from(OT_ZZ)  # indirect doctest
            False
        """
        if isinstance(S, (ExactTermMonoid,)):
            if self.growth_group.has_coerce_map_from(S.growth_group):
                return True
        else:
            return super(OTermMonoid, self)._coerce_map_from_(S)

    def _repr_(self):
        r"""
        A representation string for this `O`-term monoid.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: TermMonoid('O', G, QQ)._repr_()
            'O-Term Monoid x^ZZ with implicit coefficients in Rational Field'
        """
        return 'O-Term Monoid %s with implicit coefficients in %s' % \
            (self.growth_group._repr_short_(), self.coefficient_ring)


class TermWithCoefficient(GenericTerm):
    r"""
    Base class for asymptotic terms possessing a coefficient. For
    example, :class:`ExactTerm` directly inherits from this class.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element of
      the parent's growth group.

    - ``coefficient`` -- an element of the parent's coefficient ring.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
        sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

        sage: G = GrowthGroup('x^ZZ'); x = G.gen()
        sage: CT_ZZ = TermWithCoefficientMonoid(TermMonoid, G, ZZ)
        sage: CT_QQ = TermWithCoefficientMonoid(TermMonoid, G, QQ)
        sage: CT_ZZ(x^2, 5)
        Term with coefficient 5 and growth x^2
        sage: CT_QQ(x^3, 3/8)
        Term with coefficient 3/8 and growth x^3
    """

    def __init__(self, parent, growth, coefficient):
        r"""
        See :class:`TermWithCoefficient` for more information.

        EXAMPLES:

        First, we define some monoids::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: CT_ZZ = TermWithCoefficientMonoid(TermMonoid, G, ZZ)
            sage: CT_QQ = TermWithCoefficientMonoid(TermMonoid, G, QQ)

        The coefficients have to be from the given coefficient ring::

            sage: t = CT_ZZ(x, 1/2)
            Traceback (most recent call last):
            ...
            ValueError: 1/2 is not a coefficient in
            Generic Term Monoid x^ZZ with (implicit) coefficients in Integer Ring.
            sage: t = CT_QQ(x, 1/2); t
            Term with coefficient 1/2 and growth x

        For technical reasons, the coefficient 0 is not allowed::

            sage: t = CT_ZZ(x^42, 0)
            Traceback (most recent call last):
            ...
            ZeroCoefficientError:  Zero coefficient 0 is not allowed in
            Generic Term Monoid x^ZZ with (implicit) coefficients in Integer Ring.

        The conversion of growth elements also works for the creation
        of terms with coefficient::

            sage: x = SR('x'); x.parent()
            Symbolic Ring
            sage: CT_ZZ(x^42, 42)
            Term with coefficient 42 and growth x^42
        """
        try:
            coefficient = parent.coefficient_ring(coefficient)
        except (ValueError, TypeError):
            raise ValueError('%s is not a coefficient in %s.' %
                             (coefficient, parent))
        if coefficient == 0:
            raise ZeroCoefficientError(
                'Zero coefficient %s is not allowed in %s.' %
                (coefficient, parent))

        self.coefficient = coefficient
        super(TermWithCoefficient, self).__init__(parent=parent, growth=growth)

    def _repr_(self):
        r"""
        A representation string for this term with coefficient.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = TermWithCoefficientMonoid(TermMonoid, G, ZZ)
            sage: T(x^2, 5)._repr_()
            'Term with coefficient 5 and growth x^2'
        """
        return 'Term with coefficient %s and growth %s' % \
               (self.coefficient, self.growth)

    def _repr_product_(self, latex=False):
        r"""
        A representation string for this term with coefficient as a product.

        INPUT:

        - ``latex`` -- (default: ``False``) a boolean. If set, then
          LaTeX-output is returned.

        OUTPUT:

        A string

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: from sage.rings.asymptotic.term_monoid import TermMonoidFactory
            sage: TermMonoid = TermMonoidFactory('__main__.TermMonoid')

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = TermWithCoefficientMonoid(TermMonoid, G, ZZ)
            sage: T(x^2, 5)._repr_product_()
            '5*x^2'
        """
        if latex:
            from sage.misc.latex import latex as latex_repr
            f = latex_repr
        else:
            f = repr

        g = f(self.growth)
        c = f(self.coefficient)

        if g == '1':
            return c
        elif c == '1':
            return '{g}'.format(g=g)
        elif c == '-1':
            return '-{g}'.format(g=g)
        elif self.coefficient._is_atomic() or (-self.coefficient)._is_atomic():
            # note that -pi/2 is not atomic, but -5 is. As subtractions are handled
            # in the asymptotic ring, we ignore such non-atomicity.
            s = '{c} {g}' if latex else '{c}*{g}'
        else:
            s = r'\left({c}\right) {g}' if latex else '({c})*{g}'
        s = s.format(c=c, g=g)

        if latex:
            import re
            s = re.sub(r'([0-9])\s+([0-9])', r'\1 \\cdot \2', s)
        return s

    def _mul_(self, other):
        r"""
        Multiplication method for asymptotic terms with coefficients.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        An asymptotic term representing the product of this element
        and ``other``.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element and ``other`` have
            the same parent.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: CT = TermWithCoefficientMonoid(TermMonoid, G, ZZ)
            sage: ET = TermMonoid('exact', G, ZZ)

        This method handles the multiplication of abstract terms
        with coefficient (i.e. :class:`TermWithCoefficient`) and
        exact terms (i.e. :class:`ExactTerm`). First, an example
        for abstract terms::

            sage: t1 = CT(x^2, 2); t2 = CT(x^3, 3)
            sage: t1 * t2
            Term with coefficient 6 and growth x^5

        And now, an example for exact terms::

            sage: t1 = ET(x^2, 2); t2 = ET(x^3, 3)
            sage: t1 * t2
            6*x^5
        """
        return self.parent()(self.growth * other.growth,
                             self.coefficient * other.coefficient)

    def _calculate_pow_(self, exponent):
        r"""
        Helper function for :meth:`~ExactTerm.__pow__` which calculates
        the power of this element to the given ``exponent``.

        INPUT:

        - ``exponent`` -- an element.

        OUTPUT:

        A term.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: G = GrowthGroup('z^ZZ')
            sage: T = TermWithCoefficientMonoid(TermMonoid, G, ZZ)
            sage: t = T('2*z'); t
            Term with coefficient 2 and growth z
            sage: t._calculate_pow_(3)
            Term with coefficient 8 and growth z^3
            sage: t._calculate_pow_(-2)
            Term with coefficient 1/4 and growth z^(-2)

        ::

            sage: T = TermWithCoefficientMonoid(TermMonoid, G, CIF)
            sage: T(G.an_element(), coefficient=CIF(RIF(-1,1), RIF(-1,1)))._calculate_pow_(I)
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot take Term with coefficient 0.? + 0.?*I and
            growth z to the exponent I in Generic Term Monoid z^ZZ with
            (implicit) coefficients in Complex Interval Field with
            53 bits of precision since its coefficient 0.? + 0.?*I
            cannot be taken to this exponent.
            > *previous* ValueError: Can't take the argument of
            interval strictly containing zero
        """
        try:
            c = self.coefficient ** exponent
        except (TypeError, ValueError, ZeroDivisionError) as e:
            from .misc import combine_exceptions
            raise combine_exceptions(
                ArithmeticError('Cannot take %s to the exponent %s in %s since its '
                                'coefficient %s cannot be taken to this exponent.' %
                                (self, exponent, self.parent(), self.coefficient)), e)
        return super(TermWithCoefficient, self)._calculate_pow_(exponent, new_coefficient=c)

    def _log_coefficient_(self, base=None, locals=None):
        r"""
        Helper function to calculate the logarithm of the coefficient of this element.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.

        - ``locals`` -- a dictionary which may contain the following keys and values:

          - ``'log'`` -- value: a function. If not used, then the usual
            :class:`log <sage.functions.log.Function_log>` is taken.

        OUTPUT:

        A tuple of terms.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('exact', GrowthGroup('x^ZZ * log(x)^ZZ'), QQ)
            sage: T(3*x^2)._log_coefficient_()
            (log(3),)
            sage: T(x^1234).log_term()  # indirect doctest
            (1234*log(x),)

        .. SEEALSO::

            :meth:`ExactTerm.log_term`,
            :meth:`OTerm.log_term`.

        TESTS::

            sage: L.<log3> = QQ[]
            sage: T = TermMonoid('exact', GrowthGroup('x^ZZ * log(x)^ZZ'), L)
            sage: T(3*x^2)._log_coefficient_()
            (log(3),)
            sage: mylog=lambda z, base: log3 if z == 3 else log(z)
            sage: T(3*x^2)._log_coefficient_(locals={'log': mylog})
            (log3,)
            sage: T(3*x^2).log_term()  # indirect doctest
            (log(3), 2*log(x))
            sage: T(3*x^2).log_term(locals={'log': mylog})  # indirect doctest
            (log3, 2*log(x))
        """
        if self.coefficient.is_one():
            return tuple()
        log = self.parent().locals(locals)['log']
        return (self.parent()._create_element_in_extension_(
            self.parent().growth_group.one(),
            log(self.coefficient, base=base)),)

    def _eq_(self, other):
        r"""
        Return whether this :class:`TermWithCoefficient` is the same as
        ``other``.

        INPUT:

        - ``other`` -- an :class:`TermWithCoefficient`

        OUTPUT:

        A boolean.

        .. NOTE::

            This method gets called by the coercion model, so it can
            be assumed that this term and ``other`` come from the
            same parent.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

            sage: T = TermWithCoefficientMonoid(TermMonoid, GrowthGroup('x^ZZ'), ZZ)
            sage: t = T.an_element(); t
            Term with coefficient 1 and growth x
            sage: t == T(x, 1)
            True
            sage: t == T(x, 2)
            False
            sage: t == T(x^2, 1)
            False
        """
        return super(TermWithCoefficient, self)._eq_(other) and \
            self.coefficient == other.coefficient


class TermWithCoefficientMonoid(GenericTermMonoid):
    r"""
    This class implements the base structure for parents of
    asymptotic terms possessing a coefficient from some coefficient
    ring. In particular, this is also the parent for
    :class:`TermWithCoefficient`.

    INPUT:

    - ``growth_group`` -- a growth group.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. It has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    - ``coefficient_ring`` -- the ring which contains the
      coefficients of the elements.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
        sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

        sage: G_ZZ = GrowthGroup('x^ZZ'); x_ZZ = G_ZZ.gen()
        sage: G_QQ = GrowthGroup('x^QQ'); x_QQ = G_QQ.gen()
        sage: TC_ZZ = TermWithCoefficientMonoid(TermMonoid, G_ZZ, QQ); TC_ZZ
        Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
        sage: TC_QQ = TermWithCoefficientMonoid(TermMonoid, G_QQ, QQ); TC_QQ
        Generic Term Monoid x^QQ with (implicit) coefficients in Rational Field
        sage: TC_ZZ == TC_QQ or TC_ZZ is TC_QQ
        False
        sage: TC_QQ.coerce_map_from(TC_ZZ)
        Coercion map:
          From: Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
          To:   Generic Term Monoid x^QQ with (implicit) coefficients in Rational Field
    """

    # enable the category framework for elements
    Element = TermWithCoefficient

    def _create_element_(self, growth, coefficient):
        r"""
        Helper method which creates an element by using the ``element_class``.

        INPUT:

        - ``growth`` -- a growth element.

        - ``coefficient`` -- an element of the coefficient ring.

        OUTPUT:

        A term.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: T_ZZ = TermMonoid('exact', G_ZZ, QQ)
            sage: T_ZZ(G_ZZ.gen(), 4/3)  # indirect doctest
            4/3*x
        """
        return self.element_class(self, growth, coefficient)

    def _an_element_(self):
        r"""
        Return an element of this term with coefficient monoid.

        INPUT:

        Nothing.

        OUTPUT:

        An element of this term monoid.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: G = GrowthGroup('x^ZZ')
            sage: TermWithCoefficientMonoid(TermMonoid, G, ZZ).an_element()  # indirect doctest
            Term with coefficient 1 and growth x
            sage: TermMonoid('exact', G, ZZ).an_element()  # indirect doctest
            x
            sage: TermMonoid('exact', G, QQ).an_element()  # indirect doctest
            1/2*x
        """
        return self(self.growth_group.an_element(),
                    self.coefficient_ring.an_element())

    def some_elements(self):
        r"""
        Return some elements of this term with coefficient monoid.

        See :class:`TestSuite` for a typical use case.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from itertools import islice
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('z^QQ')
            sage: T = TermMonoid('exact', G, ZZ)
            sage: tuple(islice(T.some_elements(), int(10)))
            (z^(1/2), z^(-1/2), -z^(1/2), z^2, -z^(-1/2), 2*z^(1/2),
             z^(-2), -z^2, 2*z^(-1/2), -2*z^(1/2))
        """
        from sage.misc.mrange import cantor_product
        return iter(self(g, c) for g, c in cantor_product(
            self.growth_group.some_elements(),
            iter(c for c in self.coefficient_ring.some_elements() if c != 0)))


class ExactTerm(TermWithCoefficient):
    r"""
    Class for asymptotic exact terms. These terms primarily consist of
    an asymptotic growth element as well as a coefficient specifying
    the growth of the asymptotic term.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element from
      ``parent.growth_group``.

    - ``coefficient`` -- an element from ``parent.coefficient_ring``.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
        sage: from sage.rings.asymptotic.term_monoid import ExactTermMonoid
        sage: G = GrowthGroup('x^ZZ'); x = G.gen()
        sage: ET = ExactTermMonoid(TermMonoid, G, QQ)

    Asymptotic exact terms may be multiplied (with the usual rules
    applying)::

        sage: ET(x^2, 3) * ET(x, -1)
        -3*x^3
        sage: ET(x^0, 4) * ET(x^5, 2)
        8*x^5

    They may also be multiplied with `O`-terms::

        sage: OT = TermMonoid('O', G, QQ)
        sage: ET(x^2, 42) * OT(x)
        O(x^3)

    Absorption for asymptotic exact terms relates to addition::

        sage: ET(x^2, 5).can_absorb(ET(x^5, 12))
        False
        sage: ET(x^2, 5).can_absorb(ET(x^2, 1))
        True
        sage: ET(x^2, 5).absorb(ET(x^2, 1))
        6*x^2

    Note that, as for technical reasons, `0` is not allowed as a
    coefficient for an asymptotic term with coefficient. Instead ``None``
    is returned if two asymptotic exact terms cancel out each other
    during absorption::

        sage: ET(x^2, 42).can_absorb(ET(x^2, -42))
        True
        sage: ET(x^2, 42).absorb(ET(x^2, -42)) is None
        True

    Exact terms can also be created by converting monomials with
    coefficient from the symbolic ring, or a suitable polynomial
    or power series ring::

        sage: x = var('x'); x.parent()
        Symbolic Ring
        sage: ET(5*x^2)
        5*x^2
    """

    def _repr_(self, latex=False):
        r"""
        A representation string for this exact term.

        INPUT:

        - ``latex`` -- (default: ``False``) a boolean. If set, then
          LaTeX-output is returned.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: ET = TermMonoid('exact', G, ZZ)
            sage: ET(x^2, 2)
            2*x^2

        TESTS::

            sage: ET(x^2, 1)
            x^2
            sage: ET(x^2, -1)
            -x^2
            sage: ET(x^0, 42)
            42

        Check that :trac:`19576` is fixed::

            sage: C.<c> = AsymptoticRing('c^ZZ', SR)
            sage: (1+pi)*c
            (pi + 1)*c
            sage: R.<a> = QQ[]
            sage: S.<n> = AsymptoticRing('n^QQ', R)
            sage: (1+a)/n
            (a + 1)*n^(-1)
        """
        return self._repr_product_(latex=latex)

    def _latex_(self):
        r"""
        A LaTeX-representation string for this exact term.

        OUTPUT:

        A string.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: ET = TermMonoid('exact', G, ZZ)
            sage: latex(ET(x^2, 2))
            2 x^{2}

        ::

            sage: latex(ET(x^2, 1))
            x^{2}
            sage: latex(ET(x^2, -1))
            -x^{2}
            sage: latex(ET(x^0, 42))
            42

        ::

            sage: C.<c> = AsymptoticRing('c^ZZ', SR)
            sage: latex((1+pi)*c)
            \left(\pi + 1\right) c
            sage: R.<a> = QQ[]
            sage: S.<n> = AsymptoticRing('n^QQ', R)
            sage: latex((1+a)/n)
            \left(a + 1\right) n^{-1}

        ::

            sage: D.<d> = AsymptoticRing('QQ^d * d^ZZ', ZZ)
            sage: latex(42*3^d)
            42 \cdot 3^{d}
        """
        return self._repr_(latex=True)

    def __invert__(self):
        r"""
        Invert this term.

        OUTPUT:

        A term.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = TermMonoid('exact', G, ZZ)
            sage: ~T(x, 2)  # indirect doctest
            1/2*x^(-1)
            sage: (~T(x, 2)).parent()
            Exact Term Monoid x^ZZ with coefficients in Rational Field
        """
        try:
            c = ~self.coefficient
        except ZeroDivisionError:
            raise ZeroDivisionError('Cannot invert %s since its coefficient %s '
                                    'cannot be inverted.' % (self, self.coefficient))
        g = ~self.growth
        return self.parent()._create_element_in_extension_(g, c)

    def __pow__(self, exponent):
        r"""
        Calculate the power of this :class:`ExactTerm` to the given ``exponent``.

        INPUT:

        - ``exponent`` -- an element.

        OUTPUT:

        An :class:`ExactTerm`.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('z^ZZ')
            sage: T = TermMonoid('exact', G, ZZ)
            sage: t = T('2*z'); t
            2*z
            sage: t^3  # indirect doctest
            8*z^3
            sage: t^(1/2)  # indirect doctest
            sqrt(2)*z^(1/2)
        """
        from .misc import strip_symbolic
        return self._calculate_pow_(strip_symbolic(exponent))

    def can_absorb(self, other):
        r"""
        Check whether this exact term can absorb ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            For :class:`ExactTerm`, absorption corresponds to
            addition. This means that an exact term can absorb
            only other exact terms with the same growth.

            See the :ref:`module description <term_absorption>` for a
            detailed explanation of absorption.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: ET = TermMonoid('exact', GrowthGroup('x^ZZ'), ZZ)
            sage: t1 = ET(x^21, 1); t2 = ET(x^21, 2); t3 = ET(x^42, 1)
            sage: t1.can_absorb(t2)
            True
            sage: t2.can_absorb(t1)
            True
            sage: t1.can_absorb(t3) or t3.can_absorb(t1)
            False
        """
        return other.is_exact() and self.growth == other.growth

    def _absorb_(self, other):
        r"""
        Let this exact term absorb another exact term ``other``.

        INPUT:

        - ``other`` -- an exact term.

        OUTPUT:

        An exact term or ``None``.

        .. NOTE::

            In the context of exact terms, absorption translates
            to addition. As the coefficient `0` is not allowed,
            ``None`` is returned instead if the terms cancel out.

            See the :ref:`module description <term_absorption>` for a
            detailed explanation on absorption.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: ET = TermMonoid('exact', G, QQ)

        Asymptotic exact terms can absorb other asymptotic exact
        terms with the same growth::

            sage: et1 = ET(x^2, 2); et2 = ET(x^2, -2)
            sage: et1.absorb(et1)
            4*x^2
            sage: et1.absorb(et2) is None
            True

        If the growth differs, an ``ArithmeticException`` is raised::

            sage: ET(x^5, 1).absorb(et1)
            Traceback (most recent call last):
            ...
            ArithmeticError: x^5 cannot absorb 2*x^2
        """
        coeff_new = self.coefficient + other.coefficient
        if coeff_new.is_zero():
            return None
        else:
            return self.parent()(self.growth, coeff_new)

    def log_term(self, base=None, locals=None):
        r"""
        Determine the logarithm of this exact term.

        INPUT:

        - ``base`` -- the base of the logarithm. If ``None``
          (default value) is used, the natural logarithm is taken.


        - ``locals`` -- a dictionary which may contain the following keys and values:

          - ``'log'`` -- value: a function. If not used, then the usual
            :class:`log <sage.functions.log.Function_log>` is taken.

        OUTPUT:

        A tuple of terms.

        .. NOTE::

            This method returns a tuple with the summands that come from
            applying the rule `\log(x\cdot y) = \log(x) + \log(y)`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('exact', GrowthGroup('x^ZZ * log(x)^ZZ'), SR)
            sage: T(3*x^2).log_term()
            (log(3), 2*log(x))
            sage: T(x^1234).log_term()
            (1234*log(x),)
            sage: T(49*x^7).log_term(base=7)
            (2, 7/log(7)*log(x))

        ::

            sage: T = TermMonoid('exact', GrowthGroup('x^ZZ * log(x)^ZZ * y^ZZ * log(y)^ZZ'), SR)
            sage: T('x * y').log_term()
            (log(x), log(y))
            sage: T('4 * x * y').log_term(base=2)
            (2, 1/log(2)*log(x), 1/log(2)*log(y))

        .. SEEALSO::

            :meth:`OTerm.log_term`.
        """
        return (self._log_coefficient_(base=base, locals=locals)
                + self._log_growth_(base=base, locals=locals))

    def is_constant(self):
        r"""
        Return whether this term is an (exact) constant.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        .. NOTE::

            Only :class:`ExactTerm` with constant growth (`1`) are
            constant.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('exact', GrowthGroup('x^ZZ * log(x)^ZZ'), QQ)
            sage: T('x * log(x)').is_constant()
            False
            sage: T('3*x').is_constant()
            False
            sage: T(1/2).is_constant()
            True
            sage: T(42).is_constant()
            True
        """
        return self.growth.is_one()

    def is_little_o_of_one(self):
        r"""
        Return whether this exact term is of order `o(1)`.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('exact', GrowthGroup('x^ZZ'), QQ)
            sage: T(x).is_little_o_of_one()
            False
            sage: T(1).is_little_o_of_one()
            False
            sage: T(x^(-1)).is_little_o_of_one()
            True

        ::

            sage: T = TermMonoid('exact', GrowthGroup('x^ZZ * y^ZZ'), QQ)
            sage: T('x * y^(-1)').is_little_o_of_one()
            False
            sage: T('x^(-1) * y').is_little_o_of_one()
            False
            sage: T('x^(-2) * y^(-3)').is_little_o_of_one()
            True

        ::

            sage: T = TermMonoid('exact', GrowthGroup('x^QQ * log(x)^QQ'), QQ)
            sage: T('x * log(x)^2').is_little_o_of_one()
            False
            sage: T('x^2 * log(x)^(-1234)').is_little_o_of_one()
            False
            sage: T('x^(-1) * log(x)^4242').is_little_o_of_one()
            True
            sage: T('x^(-1/100) * log(x)^(1000/7)').is_little_o_of_one()
            True
        """
        return self.growth.is_lt_one()

    def is_exact(self):
        r"""
        Return whether this term is an exact term.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('exact', GrowthGroup('x^ZZ * log(x)^ZZ'), QQ)
            sage: T('x * log(x)').is_exact()
            True
            sage: T('3 * x^2').is_exact()
            True

        TESTS::

            sage: T = TermMonoid('O', GrowthGroup('x^ZZ'), QQ)
            sage: T('x').is_exact()
            False
        """
        return True

    def rpow(self, base):
        r"""
        Return the power of ``base`` to this exact term.

        INPUT:

        - ``base`` -- an element or ``'e'``.

        OUTPUT:

        A term.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('exact', GrowthGroup('QQ^x * x^ZZ * log(x)^ZZ'), QQ)
            sage: T('x').rpow(2)
            2^x
            sage: T('log(x)').rpow('e')
            x
            sage: T('42*log(x)').rpow('e')
            x^42
            sage: T('3*x').rpow(2)
            8^x

        ::

            sage: T('3*x^2').rpow(2)
            Traceback (most recent call last):
            ...
            ArithmeticError: Cannot construct 2^(x^2) in
            Growth Group QQ^x * x^ZZ * log(x)^ZZ * Signs^x
            > *previous* TypeError: unsupported operand parent(s) for *:
            'Growth Group QQ^x * x^ZZ * log(x)^ZZ * Signs^x' and
            'Growth Group ZZ^(x^2)'

        ::

            sage: T = TermMonoid('exact', GrowthGroup('(QQ_+)^n * n^QQ'), SR)
            sage: n = T('n')
            sage: n.rpow(2)
            2^n
            sage: _.parent()
            Exact Term Monoid QQ^n * n^QQ with coefficients in Symbolic Ring
        """
        P = self.parent()

        if self.is_constant():
            if not hasattr(base, 'parent'):
                base = P.coefficient_ring(base)
            return P._create_element_in_extension_(
                P.growth_group.one(), base ** self.coefficient)

        elem = P._create_element_in_extension_(
            self.growth.rpow(base), P.coefficient_ring.one())
        return elem ** self.coefficient

    def _substitute_(self, rules):
        r"""
        Substitute the given ``rules`` in this exact term.

        INPUT:

        - ``rules`` -- a dictionary.

        OUTPUT:

        An object.

        TESTS::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: E = TermMonoid('exact', GrowthGroup('x^ZZ'), ZZ)
            sage: e = E.an_element(); e
            x
            sage: e._substitute_({'x': SR.var('z')})
            z
            sage: E(2/x)._substitute_({'x': 0})
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot substitute in 2*x^(-1) in
            Exact Term Monoid x^ZZ with coefficients in Integer Ring.
            > *previous* ZeroDivisionError: Cannot substitute in x^(-1) in
            Growth Group x^ZZ.
            >> *previous* ZeroDivisionError: rational division by zero
            sage: E(1/x)._substitute_({'x': ZZ})
            Traceback (most recent call last):
            ...
            ValueError: Cannot substitute in x^(-1) in
            Exact Term Monoid x^ZZ with coefficients in Integer Ring.
            > *previous* ValueError: Cannot substitute in x^(-1) in Growth Group x^ZZ.
            >> *previous* ValueError: rank (=-1) must be nonnegative
        """
        try:
            g = self.growth._substitute_(rules)
        except (ArithmeticError, TypeError, ValueError) as e:
            from .misc import substitute_raise_exception
            substitute_raise_exception(self, e)

        c = self.coefficient

        try:
            return c * g
        except (ArithmeticError, TypeError, ValueError) as e:
            from .misc import substitute_raise_exception
            substitute_raise_exception(self, e)

    def _factorial_(self):
        r"""
        Return the factorial of this exact term if it is constant
        (i.e., has growth `1`).

        OUTPUT:

        A term.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('exact', 'z^QQ', QQ)
            sage: T('4')._factorial_()
            24
            sage: T('1/2')._factorial_()
            1/2*sqrt(pi)
            sage: T('4*z^(3/2)')._factorial_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot build the factorial of 4*z^(3/2)
            since it has growth != 1.
        """
        if not self.growth.is_one():
            raise ValueError(
                'Cannot build the factorial of {} since it has growth '
                '!= 1.'.format(self))

        from sage.functions.other import factorial
        return self.parent()._create_element_in_extension_(
            self.growth, factorial(self.coefficient))

    def _singularity_analysis_(self, var, zeta, precision):
        r"""
        Perform singularity analysis on this exact term.

        INPUT:

        - ``var`` -- a string denoting the variable

        - ``zeta`` -- a number

        - ``precision`` -- an integer

        OUTPUT:

        An asymptotic expansion for `[z^n] f` where `n` is ``var``
        and `f` has this term as a singular expansion
        in `T=\frac{1}{1-\frac{z}{\zeta}}\to \infty` where this element
        is a monomial in `T`.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: T = TermMonoid('exact', GrowthGroup('x^QQ'), ZZ)
            sage: T('5*x^(1/2)')._singularity_analysis_('n', 2, precision=2)
            5/sqrt(pi)*(1/2)^n*n^(-1/2)
            - 5/8/sqrt(pi)*(1/2)^n*n^(-3/2)
            + O((1/2)^n*n^(-5/2))
            sage: T('2*x^(-1)')._singularity_analysis_('n', 2, precision=3)
            Traceback (most recent call last):
            ...
            NotImplementedOZero: got O(0)
            The error term O(0) means 0 for sufficiently large n.
        """
        return self.coefficient * self.growth._singularity_analysis_(
            var=var, zeta=zeta, precision=precision)


class ExactTermMonoid(TermWithCoefficientMonoid):
    r"""
    Parent for asymptotic exact terms, implemented in
    :class:`ExactTerm`.

    INPUT:

    - ``growth_group`` -- a growth group.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. It has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    - ``coefficient_ring`` -- the ring which contains the coefficients of
      the elements.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: from sage.rings.asymptotic.term_monoid import ExactTermMonoid
        sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid

        sage: G_ZZ = GrowthGroup('x^ZZ'); x_ZZ = G_ZZ.gen()
        sage: G_QQ = GrowthGroup('x^QQ'); x_QQ = G_QQ.gen()
        sage: ET_ZZ = ExactTermMonoid(TermMonoid, G_ZZ, ZZ); ET_ZZ
        Exact Term Monoid x^ZZ with coefficients in Integer Ring
        sage: ET_QQ = ExactTermMonoid(TermMonoid, G_QQ, QQ); ET_QQ
        Exact Term Monoid x^QQ with coefficients in Rational Field
        sage: ET_QQ.coerce_map_from(ET_ZZ)
        Coercion map:
          From: Exact Term Monoid x^ZZ with coefficients in Integer Ring
          To:   Exact Term Monoid x^QQ with coefficients in Rational Field

    Exact term monoids can also be created using the
    :class:`term factory <TermMonoidFactory>`::

        sage: TermMonoid('exact', G_ZZ, ZZ) is ET_ZZ
        True
        sage: TermMonoid('exact', GrowthGroup('x^ZZ'), QQ)
        Exact Term Monoid x^ZZ with coefficients in Rational Field
    """

    # enable the category framework for elements
    Element = ExactTerm

    def _repr_(self):
        r"""
        A representation string for this exact term monoid.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: TermMonoid('exact', G, QQ)._repr_()
            'Exact Term Monoid x^ZZ with coefficients in Rational Field'
        """
        return 'Exact Term Monoid %s with coefficients in %s' % \
               (self.growth_group._repr_short_(), self.coefficient_ring)


class TermMonoidFactory(UniqueRepresentation, UniqueFactory):
    r"""
    Factory for asymptotic term monoids. It can generate the following
    term monoids:

    - :class:`OTermMonoid`,

    - :class:`ExactTermMonoid`.

    .. NOTE::

        An instance of this factory is available as ``DefaultTermMonoidFactory``.

    INPUT:

    - ``term_monoid`` -- the kind of terms held in the new term monoid.
      Either a string ``'exact'`` or ``'O'`` (capital letter ``O``),
      or an existing instance of a term monoid.

    - ``growth_group`` -- a growth group or
      a string describing a growth group.

    - ``coefficient_ring`` -- a ring.

    - ``asymptotic_ring`` -- if specified, then ``growth_group`` and
      ``coefficient_ring`` are taken from this asymptotic ring.

    OUTPUT:

    An asymptotic term monoid.

    EXAMPLES::

        sage: from sage.rings.asymptotic.growth_group import GrowthGroup
        sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
        sage: G = GrowthGroup('x^ZZ')
        sage: TermMonoid('O', G, QQ)
        O-Term Monoid x^ZZ with implicit coefficients in Rational Field
        sage: TermMonoid('exact', G, ZZ)
        Exact Term Monoid x^ZZ with coefficients in Integer Ring

    ::

        sage: R = AsymptoticRing(growth_group=G, coefficient_ring=QQ)
        sage: TermMonoid('exact', asymptotic_ring=R)
        Exact Term Monoid x^ZZ with coefficients in Rational Field
        sage: TermMonoid('O', asymptotic_ring=R)
        O-Term Monoid x^ZZ with implicit coefficients in Rational Field

        sage: TermMonoid('exact', 'QQ^m * m^QQ * log(n)^ZZ', ZZ)
        Exact Term Monoid QQ^m * m^QQ * Signs^m * log(n)^ZZ
        with coefficients in Integer Ring

    TESTS::

        sage: TermMonoid(TermMonoid('O', G, ZZ), asymptotic_ring=R)
        O-Term Monoid x^ZZ with implicit coefficients in Rational Field
        sage: TermMonoid(TermMonoid('exact', G, ZZ), asymptotic_ring=R)
        Exact Term Monoid x^ZZ with coefficients in Rational Field
        sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
        sage: TermMonoid(GenericTermMonoid(TermMonoid, G, ZZ), asymptotic_ring=R)
        Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field

    ::

        sage: TestSuite(TermMonoid('exact', GrowthGroup('x^ZZ'), QQ)).run(verbose=True)  # long time
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
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass

    ::

        sage: TestSuite(TermMonoid('O', GrowthGroup('x^QQ'), ZZ)).run(verbose=True)  # long time
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
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
    """

    def __init__(self, name,
                 exact_term_monoid_class=None,
                 O_term_monoid_class=None):
        r"""
        See :class:`TermMonoidFactory` for more information.

        TESTS::


            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import ExactTermMonoid, OTermMonoid
            sage: from sage.rings.asymptotic.term_monoid import TermMonoidFactory

            sage: class MyExactTermMonoid(ExactTermMonoid):
            ....:     pass
            sage: class MyOTermMonoid(OTermMonoid):
            ....:     pass
            sage: MyTermMonoid = TermMonoidFactory('MyTermMonoid',
            ....:                                  exact_term_monoid_class=MyExactTermMonoid,
            ....:                                  O_term_monoid_class=MyOTermMonoid)
            sage: G = GrowthGroup('x^ZZ')
            sage: type(MyTermMonoid('exact', G, QQ))
            <class '__main__.MyExactTermMonoid_with_category'>
            sage: type(MyTermMonoid('O', G, QQ))
            <class '__main__.MyOTermMonoid_with_category'>
        """
        super(TermMonoidFactory, self).__init__(name)

        if exact_term_monoid_class is None:
            exact_term_monoid_class = ExactTermMonoid
        self.ExactTermMonoid = exact_term_monoid_class

        if O_term_monoid_class is None:
            O_term_monoid_class = OTermMonoid
        self.OTermMonoid = O_term_monoid_class

    def create_key_and_extra_args(self, term_monoid,
                                  growth_group=None, coefficient_ring=None,
                                  asymptotic_ring=None,
                                  **kwds):
        r"""
        Given the arguments and keyword, create a key that uniquely
        determines this object.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ')
            sage: TermMonoid.create_key_and_extra_args('O', G, QQ)
            ((<class 'sage.rings.asymptotic.term_monoid.OTermMonoid'>,
              Growth Group x^ZZ, Rational Field), {})
            sage: TermMonoid.create_key_and_extra_args('exact', G, ZZ)
            ((<class 'sage.rings.asymptotic.term_monoid.ExactTermMonoid'>,
              Growth Group x^ZZ, Integer Ring), {})
            sage: TermMonoid.create_key_and_extra_args('exact', G)
            Traceback (most recent call last):
            ...
            ValueError: A coefficient ring has to be specified
            to create a term monoid of type 'exact'

        TESTS::

            sage: TermMonoid.create_key_and_extra_args('icecream', G)
            Traceback (most recent call last):
            ...
            ValueError: Term specification 'icecream' has to be either
            'exact' or 'O' or an instance of an existing term.
            sage: TermMonoid.create_key_and_extra_args('O', ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Integer Ring has to be an asymptotic growth group
        """
        if isinstance(term_monoid, GenericTermMonoid):
            term_class = term_monoid._underlying_class()
        elif term_monoid == 'O':
            term_class = self.OTermMonoid
        elif term_monoid == 'exact':
            term_class = self.ExactTermMonoid
        else:
            raise ValueError("Term specification '%s' has to be either 'exact' or 'O' "
                             "or an instance of an existing term." % term_monoid)

        if asymptotic_ring is not None and \
                (growth_group is not None or coefficient_ring is not None):
            raise ValueError("Input ambiguous: asymptotic ring %s as well as "
                             "growth group %s or coefficient ring %s are given." %
                             (asymptotic_ring, growth_group, coefficient_ring))

        if asymptotic_ring is not None:
            growth_group = asymptotic_ring.growth_group
            coefficient_ring = asymptotic_ring.coefficient_ring

        from .growth_group import GenericGrowthGroup
        if not isinstance(growth_group, GenericGrowthGroup):
            if isinstance(growth_group, str):
                from .growth_group import GrowthGroup
                growth_group = GrowthGroup(growth_group)
            else:
                raise ValueError('{} has to be an asymptotic growth '
                                 'group'.format(growth_group))

        if coefficient_ring is None:
            raise ValueError("A coefficient ring has to be specified to "
                             "create a term monoid of type '%s'" % (term_monoid,))

        return (term_class, growth_group, coefficient_ring), kwds

    def create_object(self, version, key, **kwds):
        r"""
        Create a object from the given arguments.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: G = GrowthGroup('x^ZZ')
            sage: TermMonoid('O', G, QQ)  # indirect doctest
            O-Term Monoid x^ZZ with implicit coefficients in Rational Field
            sage: TermMonoid('exact', G, ZZ)  # indirect doctest
            Exact Term Monoid x^ZZ with coefficients in Integer Ring
        """
        term_class, growth_group, coefficient_ring = key
        return term_class(self, growth_group, coefficient_ring, **kwds)

    def _repr_(self):
        r"""
        Return a representation string of this term monoid factory.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import TermMonoidFactory
            sage: TermMonoidFactory('TestTermMonoidFactory')  # indirect doctest
            Term Monoid Factory 'TestTermMonoidFactory'
        """
        return "Term Monoid Factory '{}'".format(self._name)

    def _cache_key(self):
        r"""
        Return a key for caching this object.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: TermMonoid._cache_key()
            (<class 'sage.rings.asymptotic.term_monoid.TermMonoidFactory'>,
             'sage.rings.asymptotic.term_monoid.DefaultTermMonoidFactory',
             <class 'sage.rings.asymptotic.term_monoid.ExactTermMonoid'>,
             <class 'sage.rings.asymptotic.term_monoid.OTermMonoid'>)
        """
        return (TermMonoidFactory,
                self._name,
                self.ExactTermMonoid,
                self.OTermMonoid)

    def __hash__(self):
        r"""
        Return a hash of this object.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory as TermMonoid
            sage: hash(TermMonoid)  # random
            42
        """
        return hash(self._cache_key())


DefaultTermMonoidFactory = TermMonoidFactory('sage.rings.asymptotic.term_monoid.DefaultTermMonoidFactory')
r"""
A factory for asymptotic term monoids.
This is an instance of :class:`TermMonoidFactory` whose documentation
provides more details.

TESTS::

    sage: from sage.rings.asymptotic.term_monoid import DefaultTermMonoidFactory
    sage: loads(dumps(DefaultTermMonoidFactory)) is DefaultTermMonoidFactory
    True
"""
