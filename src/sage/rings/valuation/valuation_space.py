r"""
Spaces of valuations

This module provides spaces of exponential pseudo-valuations on integral
domains. It currently only provides support for such valuations if they are
discrete, i.e., their image is a discrete additive subgroup of the rational
numbers extended by `\infty`.

AUTHORS:

- Julian Rüth (2016-10-14): initial version

EXAMPLES::

    sage: QQ.valuation(2).parent()
    Discrete pseudo-valuations on Rational Field

.. NOTE::

    Note that many tests not only in this module do not create instances of
    valuations directly since this gives the wrong inheritance structure on
    the resulting objects::

        sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
        sage: from sage.rings.valuation.trivial_valuation import TrivialDiscretePseudoValuation
        sage: H = DiscretePseudoValuationSpace(QQ)
        sage: v = TrivialDiscretePseudoValuation(H)
        sage: v._test_category()
        Traceback (most recent call last):
        ...
        AssertionError: False is not true

    Instead, the valuations need to be created through the
    ``__make_element_class__`` of the containing space::

        sage: from sage.rings.valuation.trivial_valuation import TrivialDiscretePseudoValuation
        sage: v = H.__make_element_class__(TrivialDiscretePseudoValuation)(H)
        sage: v._test_category()

    The factories such as ``TrivialPseudoValuation`` provide the right
    inheritance structure::

        sage: v = valuations.TrivialPseudoValuation(QQ)
        sage: v._test_category()
"""
# ****************************************************************************
#       Copyright (C) 2016-2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.homset import Homset
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.abstract_method import abstract_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.categories.action import Action


class DiscretePseudoValuationSpace(UniqueRepresentation, Homset):
    r"""
    The space of discrete pseudo-valuations on ``domain``.

    EXAMPLES::

        sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
        sage: H = DiscretePseudoValuationSpace(QQ)
        sage: QQ.valuation(2) in H
        True

    .. NOTE::

        We do not distinguish between the space of discrete valuations and the
        space of discrete pseudo-valuations. This is entirely for practical
        reasons: We would like to model the fact that every discrete valuation
        is also a discrete pseudo-valuation. At first, it seems to be
        sufficient to make sure that the ``in`` operator works which can
        essentially be achieved by overriding ``_element_constructor_`` of
        the space of discrete pseudo-valuations to accept discrete valuations
        by just returning them. Currently, however, if one does not change the
        parent of an element in ``_element_constructor_`` to ``self``, then
        one cannot register that conversion as a coercion. Consequently, the
        operators ``<=`` and ``>=`` cannot be made to work between discrete
        valuations and discrete pseudo-valuations on the same domain (because
        the implementation only calls ``_richcmp`` if both operands have the
        same parent.) Of course, we could override ``__ge__`` and  ``__le__``
        but then we would likely run into other surprises.
        So in the end, we went for a single homspace for all discrete
        valuations (pseudo or not) as this makes the implementation much
        easier.

    .. TODO::

        The comparison problem might be fixed by :issue:`22029` or similar.

    TESTS::

        sage: TestSuite(H).run() # long time
    """
    def __init__(self, domain):
        r"""
        TESTS::

            sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
            sage: isinstance(QQ.valuation(2).parent(), DiscretePseudoValuationSpace)
            True
        """
        from .value_group import DiscreteValuationCodomain
        # A valuation is a map from an additive semigroup to an additive semigroup, however, it
        # does not preserve that structure. It is therefore only a morphism in the category of sets.
        from sage.categories.sets_cat import Sets

        UniqueRepresentation.__init__(self)
        Homset.__init__(self, domain, DiscreteValuationCodomain(), category=Sets())

        from sage.categories.domains import Domains
        if domain not in Domains():
            raise ValueError("domain must be an integral domain")

    @lazy_attribute
    def _abstract_element_class(self):
        r"""
        Return an abstract base class for all valuations in this space.

        This is used to extend every valuation with a number of generic methods
        that are independent of implementation details.

        Usually, extensions of this kind would be done by implementing an
        appropriate class ``MorphismMethods`` in the category of this homset.
        However, there is no category whose arrows are the valuations, so we
        need to move this magic down to the level of the actual homset.

        EXAMPLES::

            sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
            sage: isinstance(QQ.valuation(2), DiscretePseudoValuationSpace.ElementMethods) # indirect doctest
            True
        """
        class_name = "%s._abstract_element_class" % self.__class__.__name__
        from sage.structure.dynamic_class import dynamic_class
        return dynamic_class(class_name, (super()._abstract_element_class, self.__class__.ElementMethods))

    def _get_action_(self, S, op, self_on_left):
        r"""
        Return the ``op`` action of ``S`` on elements in this space.

        EXAMPLES::

            sage: v = QQ.valuation(2)
            sage: from operator import mul
            sage: v.parent().get_action(ZZ, mul) # indirect doctest
            Right action by Integer Ring on Discrete pseudo-valuations on Rational Field
        """
        from operator import mul
        from sage.rings.infinity import InfinityRing
        from sage.rings.rational_field import QQ
        from sage.rings.integer_ring import ZZ
        if op == mul and (S is InfinityRing or S is QQ or S is ZZ):
            return ScaleAction(S, self, not self_on_left, op)
        return None

    def _an_element_(self):
        r"""
        Return a trivial valuation in this space.

        EXAMPLES::

            sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
            sage: DiscretePseudoValuationSpace(QQ).an_element() # indirect doctest
            Trivial pseudo-valuation on Rational Field
        """
        from .trivial_valuation import TrivialPseudoValuation
        return TrivialPseudoValuation(self.domain())

    def _repr_(self):
        r"""
        Return a printable representation of this space.

        EXAMPLES::

            sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
            sage: DiscretePseudoValuationSpace(QQ) # indirect doctest
            Discrete pseudo-valuations on Rational Field
        """
        return "Discrete pseudo-valuations on %r" % (self.domain(),)

    def __contains__(self, x):
        r"""
        Return whether ``x`` is a valuation in this space.

        EXAMPLES::

            sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
            sage: H = DiscretePseudoValuationSpace(QQ)
            sage: H.an_element() in H
            True
            sage: QQ.valuation(2) in H
            True
        """
        # override the logic from Homset with the original implementation for Parent
        # which entirely relies on a proper implementation of
        # _element_constructor_ and coercion maps
        from sage.structure.parent import Parent
        return Parent.__contains__(self, x)

    def __call__(self, x):
        r"""
        Create an element in this space from ``x``.

        EXAMPLES::

            sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
            sage: H = DiscretePseudoValuationSpace(QQ)
            sage: H(QQ.valuation(2))
            2-adic valuation
        """
        # override the logic from Homset with the original implementation for Parent
        # which entirely relies on a proper implementation of
        # _element_constructor_ and coercion maps
        from sage.structure.parent import Parent
        return Parent.__call__(self, x)

    def _element_constructor_(self, x):
        r"""
        Create an element in this space from ``x``.

        EXAMPLES:

        We try to convert valuations defined on different domains by changing
        their base ring::

            sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
            sage: Z = DiscretePseudoValuationSpace(ZZ)
            sage: Q = DiscretePseudoValuationSpace(QQ)
            sage: v = ZZ.valuation(2)
            sage: v in Q
            False
            sage: Q(v) in Q
            True
            sage: Q(v) in Z
            False
            sage: Z(Q(v)) in Z
            True
        """
        if isinstance(x.parent(), DiscretePseudoValuationSpace):
            if x.domain() is not self.domain():
                try:
                    return self(x.change_domain(self.domain()))
                except NotImplementedError:
                    pass
            else:
                return x
        raise ValueError("element cannot be converted into the space of %r" % (self,))

    class ElementMethods:
        r"""
        Provides methods for discrete pseudo-valuations that are added
        automatically to valuations in this space.

        EXAMPLES:

        Here is an example of a method that is automagically added to a
        discrete valuation::

            sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
            sage: H = DiscretePseudoValuationSpace(QQ)
            sage: QQ.valuation(2).is_discrete_pseudo_valuation() # indirect doctest
            True

        The methods will be provided even if the concrete type is not created
        with ``__make_element_class__``::

            sage: from sage.rings.valuation.valuation import DiscretePseudoValuation
            sage: m = DiscretePseudoValuation(H)
            sage: m.parent() is H
            True
            sage: m.is_discrete_pseudo_valuation()
            True

        However, the category framework advises you to use inheritance::

            sage: m._test_category()
            Traceback (most recent call last):
            ...
            AssertionError: False is not true

        Using ``__make_element_class__``, makes your concrete valuation inherit
        from this class::

            sage: m = H.__make_element_class__(DiscretePseudoValuation)(H)
            sage: m._test_category()
        """
        def is_discrete_pseudo_valuation(self):
            r"""
            Return whether this valuation is a discrete pseudo-valuation.

            EXAMPLES::

                sage: QQ.valuation(2).is_discrete_pseudo_valuation()
                True
            """
            return True

        @abstract_method
        def is_discrete_valuation(self):
            r"""
            Return whether this valuation is a discrete valuation, i.e.,
            whether it is a :meth:`discrete pseudo valuation
            <is_discrete_pseudo_valuation>` that only sends zero to `\infty`.

            EXAMPLES::

                sage: QQ.valuation(2).is_discrete_valuation()
                True
            """

        def is_negative_pseudo_valuation(self):
            r"""
            Return whether this valuation is a discrete pseudo-valuation that
            does attain `-\infty`, i.e., it is non-trivial and its domain
            contains an element with valuation `\infty` that has an inverse.

            EXAMPLES::

                sage: QQ.valuation(2).is_negative_pseudo_valuation()
                False
            """
            from sage.categories.fields import Fields
            if self.is_discrete_valuation():
                return False
            elif self.domain() in Fields():
                return True
            raise NotImplementedError

        @cached_method
        def is_trivial(self):
            r"""
            Return whether this valuation is trivial, i.e., whether it is
            constant `\infty` or constant zero for everything but the zero
            element.

            Subclasses need to override this method if they do not implement
            :meth:`uniformizer`.

            EXAMPLES::

                sage: QQ.valuation(7).is_trivial()
                False
            """
            from sage.rings.infinity import infinity
            if self(self.domain().one()) is infinity:
                # the constant infinity
                return True
            if self(self.uniformizer()) != 0:
                # not constant on the nonzero elements
                return False
            return True

        @abstract_method
        def uniformizer(self):
            r"""
            Return an element in the domain which has positive valuation and
            generates the value group of this valuation.

            EXAMPLES::

                sage: QQ.valuation(11).uniformizer()
                11

            Trivial valuations have no uniformizer::

                sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
                sage: v = DiscretePseudoValuationSpace(QQ).an_element()
                sage: v.is_trivial()
                True
                sage: v.uniformizer()
                Traceback (most recent call last):
                ...
                ValueError: Trivial valuations do not define a uniformizing element
            """

        @cached_method
        def value_group(self):
            r"""
            Return the value group of this discrete pseudo-valuation, the
            discrete additive subgroup of the rational numbers which is
            generated by the valuation of the :meth:`uniformizer`.

            EXAMPLES::

                sage: QQ.valuation(2).value_group()
                Additive Abelian Group generated by 1

            A pseudo-valuation that is `\infty` everywhere, does not have a
            value group::

                sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
                sage: v = DiscretePseudoValuationSpace(QQ).an_element()
                sage: v.value_group()
                Traceback (most recent call last):
                ...
                ValueError: The trivial pseudo-valuation that is infinity everywhere does not have a value group.
            """
            from .value_group import DiscreteValueGroup
            return DiscreteValueGroup(self(self.uniformizer()))

        def value_semigroup(self):
            r"""
            Return the value semigroup of this discrete pseudo-valuation, the
            additive subsemigroup of the rational numbers which is generated by
            the valuations of the elements in the domain.

            EXAMPLES:

            Most commonly, in particular over fields, the semigroup is the
            group generated by the valuation of the uniformizer::

                sage: G = QQ.valuation(2).value_semigroup(); G
                Additive Abelian Semigroup generated by -1, 1
                sage: G in AdditiveMagmas().AdditiveAssociative().AdditiveUnital().AdditiveInverse()
                True

            If the domain is a discrete valuation ring, then the semigroup
            consists of the positive elements of the :meth:`value_group`::

                sage: Zp(2).valuation().value_semigroup()
                Additive Abelian Semigroup generated by 1

            The semigroup can have a more complicated structure when the
            uniformizer is not in the domain::

                sage: v = ZZ.valuation(2)
                sage: R.<x> = ZZ[]
                sage: w = GaussValuation(R, v)
                sage: u = w.augmentation(x, 5/3)
                sage: u.value_semigroup()
                Additive Abelian Semigroup generated by 1, 5/3
            """
            from sage.categories.fields import Fields
            if self.domain() in Fields():
                from .value_group import DiscreteValueSemigroup
                # return the semigroup generated by the elements of the group
                return DiscreteValueSemigroup([]) + self.value_group()
            raise NotImplementedError("cannot determine value semigroup of %r" % (self,))

        def element_with_valuation(self, s):
            r"""
            Return an element in the domain of this valuation with valuation
            ``s``.

            EXAMPLES::

                sage: v = ZZ.valuation(2)
                sage: v.element_with_valuation(10)
                1024
            """
            from sage.rings.integer_ring import ZZ
            from sage.rings.rational_field import QQ
            s = QQ.coerce(s)
            if s not in self.value_semigroup():
                raise ValueError("s must be in the value semigroup of this valuation but %r is not in %r" % (s, self.value_semigroup()))
            if s == 0:
                return self.domain().one()
            exp = s / self.value_group().gen()
            if exp not in ZZ:
                raise NotImplementedError("s must be a multiple of %r but %r is not" % (self.value_group().gen(), s))
            ret = self.domain()(self.uniformizer() ** ZZ(exp))
            return self.simplify(ret, error=s)

        @abstract_method
        def residue_ring(self):
            r"""
            Return the residue ring of this valuation, i.e., the elements of
            nonnegative valuation modulo the elements of positive valuation.
            EXAMPLES::

                sage: QQ.valuation(2).residue_ring()
                Finite Field of size 2
                sage: valuations.TrivialValuation(QQ).residue_ring()
                Rational Field

            Note that a residue ring always exists, even when a residue field
            may not::

                sage: valuations.TrivialPseudoValuation(QQ).residue_ring()
                Quotient of Rational Field by the ideal (1)
                sage: valuations.TrivialValuation(ZZ).residue_ring()
                Integer Ring
                sage: GaussValuation(ZZ['x'], ZZ.valuation(2)).residue_ring()
                Univariate Polynomial Ring in x over Finite Field of size 2...
            """

        def residue_field(self):
            r"""
            Return the residue field of this valuation, i.e., the field of
            fractions of the :meth:`residue_ring`, the elements of nonnegative
            valuation modulo the elements of positive valuation.

            EXAMPLES::

                sage: QQ.valuation(2).residue_field()
                Finite Field of size 2
                sage: valuations.TrivialValuation(QQ).residue_field()
                Rational Field

                sage: valuations.TrivialValuation(ZZ).residue_field()
                Rational Field
                sage: GaussValuation(ZZ['x'], ZZ.valuation(2)).residue_field()
                Rational function field in x over Finite Field of size 2
            """
            ret = self.residue_ring()
            from sage.categories.fields import Fields
            if ret in Fields():
                return ret
            from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
            if isinstance(ret, PolynomialRing_generic):
                from sage.rings.function_field.constructor import FunctionField
                return FunctionField(ret.base_ring().fraction_field(), names=(ret.variable_name(),))
            return ret.fraction_field()

        @abstract_method
        def reduce(self, x):
            r"""
            Return the image of ``x`` in the :meth:`residue_ring` of this
            valuation.

            EXAMPLES::

                sage: v = QQ.valuation(2)
                sage: v.reduce(2)
                0
                sage: v.reduce(1)
                1
                sage: v.reduce(1/3)
                1
                sage: v.reduce(1/2)
                Traceback (most recent call last):
                ...
                ValueError: reduction is only defined for elements of nonnegative valuation
            """

        @abstract_method
        def lift(self, X):
            r"""
            Return a lift of ``X`` in the domain which reduces down to ``X``
            again via :meth:`reduce`.

            EXAMPLES::

                sage: v = QQ.valuation(2)
                sage: v.lift(v.residue_ring().one())
                1
            """

        def extension(self, ring):
            r"""
            Return the unique extension of this valuation to ``ring``.

            EXAMPLES::

                sage: v = ZZ.valuation(2)
                sage: w = v.extension(QQ)
                sage: w.domain()
                Rational Field
            """
            extensions = self.extensions(ring)
            assert extensions
            if len(extensions) > 1:
                raise ValueError("there is no unique extension of %r from %r to %r" % (self, self.domain(), ring))
            return extensions[0]

        def extensions(self, ring):
            r"""
            Return the extensions of this valuation to ``ring``.

            EXAMPLES::

                sage: v = ZZ.valuation(2)
                sage: v.extensions(QQ)
                [2-adic valuation]
            """
            if ring is self.domain():
                return [self]
            raise NotImplementedError("extending %r from %r to %r not implemented" % (self, self.domain(), ring))

        def restriction(self, ring):
            r"""
            Return the restriction of this valuation to ``ring``.

            EXAMPLES::

                sage: v = QQ.valuation(2)
                sage: w = v.restriction(ZZ)
                sage: w.domain()
                Integer Ring
            """
            if ring is self.domain():
                return self
            raise NotImplementedError("restricting %r from %r to %r not implemented" % (self, self.domain(), ring))

        def change_domain(self, ring):
            r"""
            Return this valuation over ``ring``.

            Unlike :meth:`extension` or :meth:`restriction`, this might not be
            completely sane mathematically. It is essentially a conversion of
            this valuation into another space of valuations.

            EXAMPLES::

                sage: v = QQ.valuation(3)
                sage: v.change_domain(ZZ)
                3-adic valuation
            """
            if ring is self.domain():
                return self
            if self.domain().is_subring(ring):
                return self.extension(ring)
            if ring.is_subring(self.domain()):
                return self.restriction(ring)
            raise NotImplementedError("changing %r from %r to %r not implemented" % (self, self.domain(), ring))

        def scale(self, scalar):
            r"""
            Return this valuation scaled by ``scalar``.

            INPUT:

            - ``scalar`` -- a nonnegative rational number or infinity

            EXAMPLES::

                sage: v = ZZ.valuation(3)
                sage: w = v.scale(3)
                sage: w(3)
                3

            Scaling can also be done through multiplication with a scalar::

                sage: w/3 == v
                True

            Multiplication by zero produces the trivial discrete valuation::

                sage: w = 0*v
                sage: w(3)
                0
                sage: w(0)
                +Infinity

            Multiplication by infinity produces the trivial discrete
            pseudo-valuation::

                sage: w = infinity*v
                sage: w(3)
                +Infinity
                sage: w(0)
                +Infinity
            """
            from sage.rings.infinity import infinity
            if scalar is infinity:
                from .trivial_valuation import TrivialPseudoValuation
                return TrivialPseudoValuation(self.domain())
            if scalar == 0:
                from .trivial_valuation import TrivialValuation
                return TrivialValuation(self.domain())
            if scalar == 1:
                return self
            if scalar < 0:
                raise ValueError("scalar must be nonnegative")
            if self.is_trivial():
                return self

            from .scaled_valuation import ScaledValuation_generic
            if isinstance(self, ScaledValuation_generic):
                return self._base_valuation.scale(scalar * self._scale)

            from .scaled_valuation import ScaledValuation
            return ScaledValuation(self, scalar)

        def separating_element(self, others):
            r"""
            Return an element in the domain of this valuation which has
            positive valuation with respect to this valuation but negative
            valuation with respect to the valuations in ``others``.

            EXAMPLES::

                sage: v2 = QQ.valuation(2)
                sage: v3 = QQ.valuation(3)
                sage: v5 = QQ.valuation(5)
                sage: v2.separating_element([v3,v5])
                4/15
            """
            try:
                iter(others)
            except TypeError:
                raise ValueError("others must be a list of valuations")

            for other in others + [self]:
                if other.parent() is not self.parent():
                    raise ValueError("all valuations must be valuations on %r but %r is a valuation on %r" % (self.domain(), other, other.domain()))
                if not other.is_discrete_valuation():
                    raise ValueError("all valuations must be discrete valuations but %r is not" % (other,))
                if other.is_trivial():
                    raise ValueError("all valuations must be non-trivial but %r is not" % (other,))

            if len(others) == 0:
                return self.uniformizer()

            # see the proof of Lemma 6.9 in http://www1.spms.ntu.edu.sg/~frederique/antchap6.pdf
            ret = self._strictly_separating_element(others[0])
            for i in range(1, len(others)):
                # ret is an element which separates self and others[:i]
                if others[i](ret) < 0:
                    # it also separates self and others[i]
                    continue

                delta = self._strictly_separating_element(others[i])
                if others[i](ret) == 0:
                    # combining powers of ret and delta, we produce a
                    # separating element for self and others[:i+1]
                    factor = ret
                    ret = delta
                    while any(other(ret) >= 0 for other in others[:i]):
                        assert others[i](ret) < 0
                        ret *= factor
                else:  # others[i](ret) > 0
                    # construct an element which approximates a unit with respect to others[i]
                    # and has negative valuation with respect to others[:i]
                    from sage.rings.semirings.non_negative_integer_semiring import NN
                    for r in iter(NN):
                        # When we enter this loop we are essentially out of
                        # luck.  The size of the coefficients is likely going
                        # through the roof here and this is not going to
                        # terminate in reasonable time.
                        factor = (ret**r)/(1+ret**r)
                        ret = factor * delta
                        if all(other(ret) < 0 for other in others[:i+1]):
                            break
            return ret

        def _strictly_separating_element(self, other):
            r"""
            Return an element in the domain of this valuation which has
            positive valuation with respect to this valuation but negative
            valuation with respect to ``other``.

            .. NOTE::

                Overriding this method tends to be a nuisance as you need to
                handle all possible types (as in Python type) of valuations.
                This is essentially the same problem that you have when
                implementing operators such as ``+`` or ``>=``. A sufficiently
                fancy multimethod implementation could solve that here but
                there is currently nothing like that in Sage/Python.

            EXAMPLES::

                sage: v2 = QQ.valuation(2)
                sage: v3 = QQ.valuation(3)
                sage: v2._strictly_separating_element(v3)
                2/3
            """
            from sage.rings.infinity import infinity

            numerator = self._weakly_separating_element(other)
            n = self(numerator)
            nn = other(numerator)
            assert n > 0
            assert nn is not infinity
            if nn < 0:
                return numerator

            denominator = other._weakly_separating_element(self)
            d = self(denominator)
            dd = other(denominator)
            assert dd > 0
            assert d is not infinity
            if d < 0:
                # The following may fail if denominator is not
                # invertible in the domain, but we don't have a better
                # option this generically.
                return self.domain()(~denominator)

            # We need nonnegative integers a and b such that
            # a*n - b*d > 0 and a*nn - b*dd < 0
            if nn == 0:
                # the above becomes b != 0 and a/b > d/n
                b = 1
                a = (d/n + 1).floor()
            else:
                # Since n,nn,d,dd are all nonnegative this is essentially equivalent to
                # a/b > d/n and b/a > nn/dd
                # which is
                # dd/nn > a/b > d/n
                assert dd / nn > d / n
                from sage.rings.continued_fraction import continued_fraction
                ab_cf = []
                dn_cf = continued_fraction(d / n)
                ddnn_cf = continued_fraction(dd / nn)
                for i, (x,y) in enumerate(zip(dn_cf, ddnn_cf)):
                    if x == y:
                        ab_cf.append(x)
                    elif x < y:
                        if y > x+1 or len(ddnn_cf) > i+1:
                            ab_cf.append(x+1)
                        else:
                            # the expansion of dd/nn is ending, so we can't append x+1
                            ab_cf.extend([x,1,1])
                    elif y < x:
                        if x > y+1 or len(dn_cf) > i+1:
                            ab_cf.append(y+1)
                        else:
                            ab_cf.extend([y,1,1])
                ab = continued_fraction(ab_cf).value()
                a,b = ab.numerator(), ab.denominator()

            ret = self.domain()(numerator**a / denominator**b)
            assert (self(ret) > 0)
            assert (other(ret) < 0)
            return ret

        def _weakly_separating_element(self, other):
            r"""
            Return an element in the domain of this valuation which has
            positive valuation with respect to this valuation and higher
            valuation with respect to this valuation than with respect to
            ``other``.

            .. NOTE::

                Overriding this method tends to be a nuisance as you need to
                handle all possible types (as in Python type) of valuations.
                This is essentially the same problem that you have when
                implementing operators such as ``+`` or ``>=``. A sufficiently
                fancy multimethod implementation could solve that here but
                there is currently nothing like that in Sage/Python.

            EXAMPLES::

                sage: v2 = QQ.valuation(2)
                sage: v3 = QQ.valuation(3)
                sage: v2._weakly_separating_element(v3)
                2
            """
            ret = self.uniformizer()
            if self(ret) > other(ret):
                return ret
            raise NotImplementedError("weakly separating element for %r and %r" % (self, other))

        def shift(self, x, s):
            r"""
            Shift ``x`` in its expansion with respect to :meth:`uniformizer` by
            ``s`` "digits".

            For nonnegative ``s``, this just returns ``x`` multiplied by a
            power of the uniformizer `\pi`.

            For negative ``s``, it does the same but when not over a field, it
            drops coefficients in the `\pi`-adic expansion which have negative
            valuation.

            EXAMPLES::

                sage: v = ZZ.valuation(2)
                sage: v.shift(1, 10)
                1024
                sage: v.shift(11, -1)
                5

            For some rings, there is no clear `\pi`-adic expansion. In this
            case, this method performs negative shifts by iterated division by
            the uniformizer and substraction of a lift of the reduction::

                sage: R.<x> = ZZ[]
                sage: v = ZZ.valuation(2)
                sage: w = GaussValuation(R, v)
                sage: w.shift(x, 1)
                2*x
                sage: w.shift(2*x, -1)
                x
                sage: w.shift(x + 2*x^2, -1)
                x^2
            """
            from sage.rings.integer_ring import ZZ
            x = self.domain().coerce(x)
            s = self.value_group()(s)
            if s == 0:
                return x

            s = ZZ(s / self.value_group().gen())
            if s > 0:
                return x * self.uniformizer()**s
            else:  # s < 0
                if ~self.uniformizer() in self.domain():
                    return self.domain()(x / self.uniformizer()**(-s))
                else:
                    for i in range(-s):
                        if self(x) < 0:
                            raise NotImplementedError("cannot compute general shifts over non-fields which contain elements of negative valuation")
                        x -= self.lift(self.reduce(x))
                        x //= self.uniformizer()
                    return x

        def simplify(self, x, error=None, force=False):
            r"""
            Return a simplified version of ``x``.

            Produce an element which differs from ``x`` by an element of
            valuation strictly greater than the valuation of ``x`` (or strictly
            greater than ``error`` if set.)

            If ``force`` is not set, then expensive simplifications may be avoided.

            EXAMPLES::

                sage: v = ZZ.valuation(2)
                sage: v.simplify(6, force=True)
                2
                sage: v.simplify(6, error=0, force=True)
                0
            """
            x = self.domain().coerce(x)

            if error is not None and self(x) > error:
                return self.domain().zero()
            return x

        def lower_bound(self, x):
            r"""
            Return a lower bound of this valuation at ``x``.

            Use this method to get an approximation of the valuation of ``x``
            when speed is more important than accuracy.

            EXAMPLES::

                sage: v = ZZ.valuation(2)
                sage: v.lower_bound(2^10)
                10
            """
            return self(x)

        def upper_bound(self, x):
            r"""
            Return an upper bound of this valuation at ``x``.

            Use this method to get an approximation of the valuation of ``x``
            when speed is more important than accuracy.

            EXAMPLES::

                sage: v = ZZ.valuation(2)
                sage: v.upper_bound(2^10)
                10
            """
            return self(x)

        def inverse(self, x, precision):
            r"""
            Return an approximate inverse of ``x``.

            The element returned is such that the product differs from 1 by an
            element of valuation at least ``precision``.

            INPUT:

            - ``x`` -- an element in the domain of this valuation

            - ``precision`` -- a rational or infinity

            EXAMPLES::

                sage: v = ZZ.valuation(2)
                sage: x = 3
                sage: y = v.inverse(3, 2); y
                3
                sage: x*y - 1
                8

            This might not be possible for elements of positive valuation::

                sage: v.inverse(2, 2)
                Traceback (most recent call last):
                ...
                ValueError: element has no approximate inverse in this ring

            Of course this always works over fields::

                sage: v = QQ.valuation(2)
                sage: v.inverse(2, 2)
                1/2
            """
            return x.inverse_of_unit()

        def _relative_size(self, x):
            r"""
            Return an estimate on the coefficient size of ``x``.

            The number returned is an estimate on the factor between the number of
            bits used by ``x`` and the minimal number of bits used by an element
            congruent to ``x``.

            This is used by :meth:`simplify` to decide whether simplification of
            coefficients is going to lead to a significant shrinking of the
            coefficients of ``x``.

            EXAMPLES::

                sage: v = Qp(2).valuation()
                sage: v._relative_size(2)
                1

            Some valuations do not overwrite this method because simplification
            does not increase the speed of valuations, e.g., some `p`-adic
            valuations::

                sage: v._relative_size(2**20)
                1
            """
            return 1

        def _test_is_negative_pseudo_valuation(self, **options):
            r"""
            Check that :meth:`is_negative_pseudo_valuation` works correctly.

            TESTS::

                sage: v = ZZ.valuation(3)
                sage: v._test_is_negative_pseudo_valuation()
            """
            tester = self._tester(**options)

            if self.is_discrete_valuation():
                tester.assertFalse(self.is_negative_pseudo_valuation())
                return

            if not self.is_negative_pseudo_valuation():
                X = self.domain().some_elements()
                for x in tester.some_elements(X):
                    from sage.rings.infinity import infinity
                    tester.assertNotEqual(self(x), -infinity)

        def _test_bounds(self, **options):
            r"""
            Check that :meth:`lower_bound` and :meth:`upper_bound` work
            correctly.

            TESTS::

                sage: v = ZZ.valuation(3)
                sage: v._test_bounds()
            """
            tester = self._tester(**options)

            X = self.domain().some_elements()
            for x in tester.some_elements(X):
                tester.assertGreaterEqual(self.upper_bound(x), self(x))
                tester.assertLessEqual(self.lower_bound(x), self(x))

        def _test_simplify(self, **options):
            r"""
            Check that :meth:`simplify` works correctly.

            TESTS::

                sage: v = ZZ.valuation(3)
                sage: v._test_simplify()
            """
            tester = self._tester(**options)

            try:
                self.residue_ring()
                has_residue_ring = True
            except NotImplementedError:
                # over non-fields (and especially polynomial rings over
                # non-fields) computation of the residue ring is often
                # difficult and not very interesting
                from sage.categories.fields import Fields
                if self.domain() not in Fields():
                    return
                raise

            X = self.domain().some_elements()
            for x in tester.some_elements(X):
                y = self.simplify(x)
                tester.assertEqual(self(x), self(y))
                if self(x) >= 0 and has_residue_ring:
                    tester.assertEqual(self.reduce(x), self.reduce(y))

            if self.is_trivial() and not self.is_discrete_valuation():
                return

            S = self.value_group().some_elements()
            from itertools import product
            for x,s in tester.some_elements(product(X, S)):
                y = self.simplify(x, error=s)
                if self.domain().is_exact():
                    tester.assertGreaterEqual(self(x-y), s)
                elif hasattr(y, 'precision_absolute'):
                    tester.assertGreaterEqual(self(x-y), min(s, y.precision_absolute()))

        def _test_shift(self, **options):
            r"""
            Check that :meth:`shift` works correctly.

            TESTS::

                sage: v = ZZ.valuation(3)
                sage: v._test_shift()
            """
            if self.is_trivial() and not self.is_discrete_valuation():
                return

            try:
                self.residue_ring()
            except Exception:
                # it is not clear what a shift should be in this case
                return

            tester = self._tester(**options)
            X = self.domain().some_elements()
            S = self.value_group().some_elements()
            from itertools import product
            for x,s in tester.some_elements(product(X, S)):
                if self(x) < 0 and ~self.uniformizer() not in self.domain():
                    # it is not clear what a shift should be in this case
                    continue
                y = self.shift(x, s)
                if s >= 0:
                    tester.assertGreaterEqual(self(y),self(x))
                from sage.categories.fields import Fields
                if self.domain().is_exact() and self.domain() in Fields():
                    # the shift here sometimes fails if elements implement
                    # __floordiv__ incorrectly, see #23971
                    tester.assertEqual(x, self.shift(y, -s))

        def _test_scale(self, **options):
            r"""
            Check that :meth:`scale` works correctly.

            TESTS::

                sage: v = ZZ.valuation(3)
                sage: v._test_scale()
            """
            tester = self._tester(**options)

            from sage.rings.rational_field import QQ
            from sage.rings.infinity import infinity
            from .trivial_valuation import TrivialValuation, TrivialPseudoValuation

            tester.assertEqual(QQ(0)*self, TrivialValuation(self.domain()))
            tester.assertEqual(infinity*self, TrivialPseudoValuation(self.domain()))

            for s in tester.some_elements(QQ.some_elements()):
                if s < 0:
                    with tester.assertRaises(ValueError):
                        s * self
                    continue
                if s == 0:
                    continue

                scaled = s * self

                tester.assertEqual(self.is_trivial(), scaled.is_trivial())
                if not self.is_trivial():
                    tester.assertEqual(self.uniformizer(), scaled.uniformizer())
                    tester.assertEqual(scaled(self.uniformizer()), s * self(self.uniformizer()))
                unscaled = scaled / s
                tester.assertEqual(self, unscaled)

        def _test_add(self, **options):
            r"""
            Check that the (strict) triangle equality is satisfied for the
            valuation of this ring.

            TESTS::

                sage: v = ZZ.valuation(3)
                sage: v._test_add()
            """
            tester = self._tester(**options)
            S = self.domain().some_elements()
            from itertools import product
            for x, y in tester.some_elements(product(S, S)):
                tester.assertGreaterEqual(self(x + y), min(self(x), self(y)))
                if self(x) != self(y):
                    tester.assertEqual(self(x + y), min(self(x), self(y)))

        def _test_infinite_zero(self, **options):
            r"""
            Check that zero is sent to infinity.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_infinite_zero()
            """
            tester = self._tester(**options)
            from sage.rings.infinity import infinity
            tester.assertEqual(self(self.domain().zero()), infinity)

        def _test_mul(self, **options):
            r"""
            Check that multiplication translates to addition of valuations.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_mul()
            """
            from sage.rings.infinity import infinity
            from itertools import product
            tester = self._tester(**options)
            S = self.domain().some_elements()
            infis = {infinity, -infinity}

            for x, y in tester.some_elements(product(S, S)):
                if {Sx := self(x), Sy := self(y)} == infis:
                    continue
                tester.assertEqual(self(x * y), Sx + Sy)

        def _test_no_infinite_units(self, **options):
            r"""
            Check that no units are sent to infinity.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_no_infinite_units()

            As multiplication translates to addition, pseudo-valuations which
            send a unit to infinity are necessarily trivial::

                sage: from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
                sage: v = DiscretePseudoValuationSpace(QQ).an_element()
                sage: v(1)
                +Infinity
                sage: v.is_trivial()
                True
            """
            if not self.is_discrete_valuation() and self.is_trivial():
                return
            if self.is_negative_pseudo_valuation():
                return

            from sage.rings.infinity import infinity
            tester = self._tester(**options)
            for x in tester.some_elements(self.domain().some_elements()):
                if self(x) is infinity:
                    tester.assertFalse(x.is_unit())

        def _test_value_group(self, **options):
            r"""
            Check correctness of the value group.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_value_group()
            """
            from sage.rings.infinity import infinity
            tester = self._tester(**options)
            # check consistency of trivial valuations first
            if self.is_trivial():
                if self(self.domain().one()) is infinity:
                    # a trivial pseudo-valuation that sends everything to infinity
                    with tester.assertRaises(ValueError):
                        self.value_group()
                    return

            # check that all valuations are in the value group
            for x in tester.some_elements(self.domain().some_elements()):
                if self(x) is not infinity and self(x) is not -infinity:
                    tester.assertIn(self(x), self.value_group())

            if not self.is_trivial():
                # check that the uniformizer generates the value group
                tester.assertEqual(self.value_group().gen(), self(self.uniformizer()))

        def _test_value_semigroup(self, **options):
            r"""
            Check correctness of the value semigroup.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_value_semigroup()                                         # needs sage.geometry.polyhedron
            """
            tester = self._tester(**options)

            if self.is_trivial() and not self.is_discrete_valuation():
                # the trivial pseudo-valuation does not have a value semigroup
                return

            for s in tester.some_elements(self.value_semigroup().some_elements()):
                tester.assertIn(s, self.value_group())

        def _test_element_with_valuation(self, **options):
            r"""
            Check correctness of :meth:`element_with_valuation`.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_element_with_valuation()                                  # needs sage.geometry.polyhedron
            """
            tester = self._tester(**options)

            if self.is_trivial() and not self.is_discrete_valuation():
                # the trivial pseudo-valuation does not have a value semigroup
                return

            for s in tester.some_elements(self.value_semigroup().some_elements()):
                tester.assertEqual(self(self.element_with_valuation(s)), s)

        def _test_residue_ring(self, **options):
            r"""
            Check the correctness of residue rings.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_residue_ring()
            """
            tester = self._tester(**options)

            try:
                r = self.residue_ring()
            except NotImplementedError:
                # over non-fields (and especially polynomial rings over
                # non-fields) computation of the residue ring is often
                # difficult and not very interesting
                from sage.categories.fields import Fields
                if self.domain() not in Fields():
                    return
                raise

            if r.zero() == r.one():
                # residue ring is the zero rng
                tester.assertGreater(self(1), 0)
                return

            c = self.residue_ring().characteristic()
            if c != 0:
                tester.assertGreater(self(c), 0)

        def _test_reduce(self, **options):
            r"""
            Check the correctness of reductions.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_reduce()
            """
            tester = self._tester(**options)

            try:
                self.residue_ring()
            except NotImplementedError:
                # over non-fields (and especially polynomial rings over
                # non-fields) computation of the residue ring is often
                # difficult and not very interesting
                from sage.categories.fields import Fields
                if self.domain() not in Fields():
                    return
                raise

            for x in tester.some_elements(self.domain().some_elements()):
                if self(x) < 0:
                    with tester.assertRaises((ValueError, ArithmeticError)):
                        self.reduce(x)
                    continue
                if self(x) == 0:
                    y = self.reduce(x)
                    tester.assertIn(y, self.residue_ring())
                    tester.assertNotEqual(y, 0)
                    if x.is_unit() and ~x in self.domain():
                        tester.assertTrue(y.is_unit())
                        tester.assertIn(~y, self.residue_ring())
                        tester.assertEqual(~y, self.reduce(self.domain()(~x)))
                if self(x) > 0:
                    tester.assertEqual(self.reduce(x), 0)

        def _test_lift(self, **options):
            r"""
            Check the correctness of lifts.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_lift()
            """
            tester = self._tester(**options)

            try:
                self.residue_ring()
            except NotImplementedError:
                # over non-fields (and especially polynomial rings over
                # non-fields) computation of the residue ring is often
                # difficult and not very interesting
                from sage.categories.fields import Fields
                if self.domain() not in Fields():
                    return
                raise

            for X in tester.some_elements(self.residue_ring().some_elements()):
                x = self.lift(X)
                y = self.reduce(x)
                tester.assertEqual(X, y)
                if X != 0:
                    tester.assertEqual(self(x), 0)

        def _test_restriction(self, **options):
            r"""
            Check the correctness of reductions.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_restriction()
            """
            tester = self._tester(**options)

            tester.assertEqual(self.restriction(self.domain()), self)

        def _test_extension(self, **options):
            r"""
            Check the correctness of extensions.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_extension()
            """
            tester = self._tester(**options)

            tester.assertEqual(self.extension(self.domain()), self)
            tester.assertEqual(self.extensions(self.domain()), [self])

        def _test_change_domain(self, **options):
            r"""
            Check the correctness of :meth:`change_domain`.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_change_domain()
            """
            tester = self._tester(**options)

            tester.assertEqual(self.change_domain(self.domain()), self)

        def _test_no_infinite_nonzero(self, **options):
            r"""
            Check that only zero is sent to infinity.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_no_infinite_nonzero()
            """
            if not self.is_discrete_valuation():
                return

            from sage.rings.infinity import infinity
            tester = self._tester(**options)
            for x in tester.some_elements(self.domain().some_elements()):
                if self(x) is infinity:
                    tester.assertEqual(x, 0)

        def _test_residue_field(self, **options):
            r"""
            Check the correctness of residue fields.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_residue_field()
            """
            if not self.is_discrete_valuation():
                return

            tester = self._tester(**options)
            try:
                self.residue_field()
            except ValueError:
                from sage.categories.fields import Fields
                # a discrete valuation on a field has a residue field
                tester.assertNotIn(self.domain(), Fields())
                return
            except NotImplementedError:
                # over non-fields (and especially polynomial rings over
                # non-fields) computation of the residue ring is often
                # difficult and not very interesting
                from sage.categories.fields import Fields
                if self.domain() not in Fields():
                    return
                raise

            try:
                r = self.residue_ring()
            except Exception:
                # If the residue ring cannot be constructed for some reason
                # then we do not check its relation to the residue field.
                # _test_residue_ring() is responsible for checking whether the
                # residue ring should be constructible or not.
                pass
            else:
                # the residue ring must coerce into the residue field
                tester.assertTrue(self.residue_field().has_coerce_map_from(r))

            c = self.residue_field().characteristic()
            if c != 0:
                tester.assertGreater(self(c), 0)

        def _test_ge(self, **options):
            r"""
            Check the correctness of the ``>=`` operator.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_ge()
            """
            tester = self._tester(**options)

            tester.assertGreaterEqual(self, self)

            if self.is_negative_pseudo_valuation():
                return

            from .trivial_valuation import TrivialPseudoValuation, TrivialValuation
            tester.assertGreaterEqual(self, TrivialValuation(self.domain()))
            tester.assertLessEqual(self, TrivialPseudoValuation(self.domain()))

        def _test_le(self, **options):
            r"""
            Check the correctness of the ``<=`` operator.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_le()
            """
            tester = self._tester(**options)

            tester.assertGreaterEqual(self, self)

            if self.is_negative_pseudo_valuation():
                return

            from .trivial_valuation import TrivialPseudoValuation, TrivialValuation
            tester.assertLessEqual(TrivialValuation(self.domain()), self)
            tester.assertGreaterEqual(TrivialPseudoValuation(self.domain()), self)

        def _test_inverse(self, **options):
            r"""
            Check the correctness of :meth:`inverse`.

            TESTS::

                sage: v = QQ.valuation(5)
                sage: v._test_inverse()
            """
            tester = self._tester(**options)

            for x in tester.some_elements(self.domain().some_elements()):
                from sage.rings.infinity import infinity
                for prec in (0, 1, 42, infinity):
                    try:
                        y = self.inverse(x, prec)
                    except ArithmeticError:  # Inverse does not exist
                        continue
                    except ValueError:
                        if prec is not infinity:
                            tester.assertNotEqual(self(x), 0)
                        tester.assertFalse(x.is_unit())
                        continue

                    tester.assertIs(y.parent(), self.domain())
                    if self.domain().is_exact():
                        tester.assertGreaterEqual(self(x * y - 1), prec)


class ScaleAction(Action):
    r"""
    Action of integers, rationals and the infinity ring on valuations by
    scaling it.

    EXAMPLES::

        sage: v = QQ.valuation(5)
        sage: from operator import mul
        sage: v.parent().get_action(ZZ, mul, self_on_left=False)
        Left action by Integer Ring on Discrete pseudo-valuations on Rational Field
    """
    def _act_(self, s, v):
        r"""
        Let ``s`` act on ``v``.

        EXAMPLES::

            sage: v = QQ.valuation(5)
            sage: 3 * v  # indirect doctest
            3 * 5-adic valuation
        """
        return v.scale(s)
