# sage_setup: distribution = sagemath-pari
r"""
Trivial valuations

AUTHORS:

- Julian Rüth (2016-10-14): initial version

EXAMPLES::

    sage: v = valuations.TrivialValuation(QQ); v
    Trivial valuation on Rational Field
    sage: v(1)
    0
"""
# ****************************************************************************
#       Copyright (C) 2016-2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from .valuation import DiscretePseudoValuation, DiscreteValuation, InfiniteDiscretePseudoValuation
from .valuation_space import DiscretePseudoValuationSpace
from sage.structure.factory import UniqueFactory


class TrivialValuationFactory(UniqueFactory):
    r"""
    Create a trivial valuation on ``domain``.

    EXAMPLES::

        sage: v = valuations.TrivialValuation(QQ); v
        Trivial valuation on Rational Field
        sage: v(1)
        0
    """
    def __init__(self, clazz, parent, *args, **kwargs):
        r"""
        TESTS::

            sage: from sage.rings.valuation.trivial_valuation import TrivialValuationFactory
            sage: isinstance(valuations.TrivialValuation, TrivialValuationFactory)
            True
        """
        UniqueFactory.__init__(self, *args, **kwargs)
        self._class = clazz
        self._parent = parent

    def create_key(self, domain):
        r"""
        Create a key that identifies this valuation.

        EXAMPLES::

            sage: valuations.TrivialValuation(QQ) is valuations.TrivialValuation(QQ) # indirect doctest
            True
        """
        return domain,

    def create_object(self, version, key, **extra_args):
        r"""
        Create a trivial valuation from ``key``.

        EXAMPLES::

            sage: valuations.TrivialValuation(QQ) # indirect doctest
            Trivial valuation on Rational Field
        """
        domain, = key
        parent = self._parent(domain)
        return parent.__make_element_class__(self._class)(parent)


class TrivialDiscretePseudoValuation_base(DiscretePseudoValuation):
    r"""
    Base class for code shared by trivial valuations.

    EXAMPLES::

        sage: v = valuations.TrivialPseudoValuation(ZZ); v
        Trivial pseudo-valuation on Integer Ring

    TESTS::

        sage: TestSuite(v).run() # long time
    """
    def uniformizer(self):
        r"""
        Return a uniformizing element for this valuation.

        EXAMPLES::

            sage: v = valuations.TrivialPseudoValuation(ZZ)
            sage: v.uniformizer()
            Traceback (most recent call last):
            ...
            ValueError: Trivial valuations do not define a uniformizing element
        """
        raise ValueError("Trivial valuations do not define a uniformizing element")

    def is_trivial(self):
        r"""
        Return whether this valuation is trivial.

        EXAMPLES::

            sage: v = valuations.TrivialPseudoValuation(QQ)
            sage: v.is_trivial()
            True
        """
        return True

    def is_negative_pseudo_valuation(self):
        r"""
        Return whether this valuation attains the value `-\infty`.

        EXAMPLES::

            sage: v = valuations.TrivialPseudoValuation(QQ)
            sage: v.is_negative_pseudo_valuation()
            False
        """
        return False


class TrivialDiscretePseudoValuation(TrivialDiscretePseudoValuation_base, InfiniteDiscretePseudoValuation):
    r"""
    The trivial pseudo-valuation that is `\infty` everywhere.

    EXAMPLES::

        sage: v = valuations.TrivialPseudoValuation(QQ); v
        Trivial pseudo-valuation on Rational Field

    TESTS::

        sage: TestSuite(v).run() # long time
    """
    def __init__(self, parent):
        r"""
        TESTS::

            sage: from sage.rings.valuation.trivial_valuation import TrivialDiscretePseudoValuation
            sage: v = valuations.TrivialPseudoValuation(QQ)
            sage: isinstance(v, TrivialDiscretePseudoValuation)
            True
        """
        TrivialDiscretePseudoValuation_base.__init__(self, parent)
        InfiniteDiscretePseudoValuation.__init__(self, parent)

    def _call_(self, x):
        r"""
        Evaluate this valuation at ``x``.

        EXAMPLES::

            sage: v = valuations.TrivialPseudoValuation(QQ)
            sage: v(0)
            +Infinity
            sage: v(1)
            +Infinity
        """
        from sage.rings.infinity import infinity
        return infinity

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: valuations.TrivialPseudoValuation(QQ) # indirect doctest
            Trivial pseudo-valuation on Rational Field
        """
        return "Trivial pseudo-valuation on %r" % (self.domain(),)

    def value_group(self):
        r"""
        Return the value group of this valuation.

        EXAMPLES:

        A trivial discrete pseudo-valuation has no value group::

            sage: v = valuations.TrivialPseudoValuation(QQ)
            sage: v.value_group()
            Traceback (most recent call last):
            ...
            ValueError: The trivial pseudo-valuation that is infinity everywhere does not have a value group.
        """
        raise ValueError("The trivial pseudo-valuation that is infinity everywhere does not have a value group.")

    def residue_ring(self):
        r"""
        Return the residue ring of this valuation.

        EXAMPLES::

            sage: valuations.TrivialPseudoValuation(QQ).residue_ring()
            Quotient of Rational Field by the ideal (1)
        """
        return self.domain().quo(self.domain().one())

    def reduce(self, x):
        r"""
        Reduce ``x`` modulo the positive elements of this valuation.

        EXAMPLES::

            sage: v = valuations.TrivialPseudoValuation(QQ)
            sage: v.reduce(1)
            0
        """
        self.domain().coerce(x)
        return self.residue_ring().zero()

    def lift(self, X):
        r"""
        Return a lift of ``X`` to the domain of this valuation.

        EXAMPLES::

            sage: v = valuations.TrivialPseudoValuation(QQ)
            sage: v.lift(v.residue_ring().zero())
            0
        """
        self.residue_ring().coerce(X)  # ignore the output
        return self.domain().zero()

    def _ge_(self, other):
        r"""
        Return whether this valuation is bigger or equal than ``other``
        everywhere.

        EXAMPLES::

            sage: v = valuations.TrivialPseudoValuation(QQ)
            sage: w = valuations.TrivialValuation(QQ)
            sage: v >= w
            True
        """
        # the trivial discrete valuation is the biggest valuation
        return True


class TrivialDiscreteValuation(TrivialDiscretePseudoValuation_base, DiscreteValuation):
    r"""
    The trivial valuation that is zero on nonzero elements.

    EXAMPLES::

        sage: v = valuations.TrivialValuation(QQ); v
        Trivial valuation on Rational Field

    TESTS::

        sage: TestSuite(v).run() # long time
    """
    def __init__(self, parent):
        r"""
        TESTS::

            sage: from sage.rings.valuation.trivial_valuation import TrivialDiscreteValuation
            sage: v = valuations.TrivialValuation(QQ)
            sage: isinstance(v, TrivialDiscreteValuation)
            True
        """
        TrivialDiscretePseudoValuation_base.__init__(self, parent)
        DiscreteValuation.__init__(self, parent)

    def _call_(self, x):
        r"""
        Evaluate this valuation at ``x``.

        EXAMPLES::

            sage: v = valuations.TrivialValuation(QQ)
            sage: v(0)
            +Infinity
            sage: v(1)
            0
        """
        from sage.rings.infinity import infinity
        return infinity if x == 0 else self.codomain().zero()

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: valuations.TrivialValuation(QQ) # indirect doctest
            Trivial valuation on Rational Field
        """
        return "Trivial valuation on %r" % (self.domain(),)

    def value_group(self):
        r"""
        Return the value group of this valuation.

        EXAMPLES:

        A trivial discrete valuation has a trivial value group::

            sage: v = valuations.TrivialValuation(QQ)
            sage: v.value_group()
            Trivial Additive Abelian Group
        """
        from .value_group import DiscreteValueGroup
        return DiscreteValueGroup(0)

    def residue_ring(self):
        r"""
        Return the residue ring of this valuation.

        EXAMPLES::

            sage: valuations.TrivialValuation(QQ).residue_ring()
            Rational Field
        """
        return self.domain()

    def reduce(self, x):
        r"""
        Reduce ``x`` modulo the positive elements of this valuation.

        EXAMPLES::

            sage: v = valuations.TrivialValuation(QQ)
            sage: v.reduce(1)
            1
        """
        return self.domain().coerce(x)

    def lift(self, X):
        r"""
        Return a lift of ``X`` to the domain of this valuation.

        EXAMPLES::

            sage: v = valuations.TrivialValuation(QQ)
            sage: v.lift(v.residue_ring().zero())
            0
        """
        return self.residue_ring().coerce(X)

    def extensions(self, ring):
        r"""
        Return the unique extension of this valuation to ``ring``.

        EXAMPLES::

            sage: v = valuations.TrivialValuation(ZZ)
            sage: v.extensions(QQ)
            [Trivial valuation on Rational Field]
        """
        if self.domain().is_subring(ring):
            return [TrivialValuation(ring)]
        return super(DiscretePseudoValuation, self).extensions(ring)

    def _ge_(self, other):
        r"""
        Return whether this valuation is bigger or equal than ``other``
        everywhere.

        EXAMPLES::

            sage: v = valuations.TrivialPseudoValuation(QQ)
            sage: w = valuations.TrivialValuation(QQ)
            sage: w >= v
            False
        """
        # the trivial discrete valuation is the smallest valuation
        return self is other


TrivialValuation = TrivialValuationFactory(TrivialDiscreteValuation, DiscretePseudoValuationSpace, "sage.rings.valuation.trivial_valuation.TrivialValuation")
TrivialPseudoValuation = TrivialValuationFactory(TrivialDiscretePseudoValuation, DiscretePseudoValuationSpace, "sage.rings.valuation.trivial_valuation.TrivialPseudoValuation")
