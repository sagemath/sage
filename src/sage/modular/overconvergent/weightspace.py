# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.rings.padics
r"""
The space of `p`-adic weights

A `p`-adic weight is a continuous character `\ZZ_p^\times \to
\CC_p^\times`. These are the `\CC_p`-points of a rigid space over
`\QQ_p`, which is isomorphic to a disjoint union of copies (indexed by
`(\ZZ/p\ZZ)^\times`) of the open unit `p`-adic disc.

Sage supports both "classical points", which are determined by the data of a
Dirichlet character modulo `p^m` for some `m` and an integer `k` (corresponding
to the character `z \mapsto z^k \chi(z)`) and "non-classical points" which are
determined by the data of an element of `(\ZZ/p\ZZ)^\times` and
an element `w \in \CC_p` with `|w - 1| < 1`.

EXAMPLES::

    sage: W = pAdicWeightSpace(17)
    sage: W
    Space of 17-adic weight-characters
     defined over 17-adic Field with capped relative precision 20
    sage: R.<x> = QQ[]
    sage: L = Qp(17).extension(x^2 - 17, names='a'); L.rename('L')
    sage: W.base_extend(L)
    Space of 17-adic weight-characters defined over L

We create a simple element of `\mathcal{W}`: the algebraic character, `x \mapsto x^6`::

    sage: kappa = W(6)
    sage: kappa(5)
    15625
    sage: kappa(5) == 5^6
    True

A locally algebraic character, `x \mapsto x^6 \chi(x)` for `\chi` a Dirichlet
character mod `p`::

    sage: kappa2 = W(6, DirichletGroup(17, Qp(17)).0^8)
    sage: kappa2(5) == -5^6
    True
    sage: kappa2(13) == 13^6
    True

A non-locally-algebraic character, sending the generator 18 of `1 + 17
\ZZ_{17}` to 35 and acting as `\mu \mapsto \mu^4` on the group of 16th
roots of unity::

    sage: kappa3 = W(35 + O(17^20), 4, algebraic=False)
    sage: kappa3(2)
    16 + 8*17 + ... + O(17^20)

AUTHORS:

- David Loeffler (2008-9)
"""

# ****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                     2008-9 David Loeffler <d.loeffler.01@cantab.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import weakref

from sage.arith.misc import divisors
from sage.categories.sets_cat import Sets
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.modular.dirichlet import DirichletGroup, trivial_character
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.rational_field import QQ
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp

lazy_import('sage.rings.padics.factory', 'Qp')
lazy_import('sage.rings.padics.padic_generic_element', 'pAdicGenericElement')


_wscache = {}


def WeightSpace_constructor(p, base_ring=None):
    r"""
    Construct the `p`-adic weight space for the given prime p.

    A `p`-adic weight
    is a continuous character `\ZZ_p^\times \to \CC_p^\times`.
    These are the `\CC_p`-points of a rigid space over `\QQ_p`,
    which is isomorphic to a disjoint union of copies (indexed by
    `(\ZZ/p\ZZ)^\times`) of the open unit `p`-adic disc.

    Note that the "base ring" of a `p`-adic weight is the smallest ring
    containing the image of `\ZZ`; in particular, although the default base
    ring is `\QQ_p`, base ring `\QQ` will also work.

    EXAMPLES::

        sage: pAdicWeightSpace(3) # indirect doctest
        Space of 3-adic weight-characters
         defined over 3-adic Field with capped relative precision 20
        sage: pAdicWeightSpace(3, QQ)
        Space of 3-adic weight-characters defined over Rational Field
        sage: pAdicWeightSpace(10)
        Traceback (most recent call last):
        ...
        ValueError: p must be prime
    """
    if base_ring is None:
        base_ring = Qp(p)
    if (p, base_ring) in _wscache:
        m = _wscache[(p, base_ring)]()
        if m is not None:
            return m
    m = WeightSpace_class(p, base_ring)
    _wscache[(p, base_ring)] = weakref.ref(m)
    return m


class WeightSpace_class(Parent):
    r"""
    The space of `p`-adic weight-characters `\mathcal{W} = {\rm
    Hom}(\ZZ_p^\times, \CC_p^\times)`.

    This is isomorphic to a
    disjoint union of `(p-1)` open discs of radius 1 (or 2 such discs if `p =
    2`), with the parameter on the open disc corresponding to the image of `1 +
    p` (or 5 if `p = 2`)

    TESTS::

        sage: W = pAdicWeightSpace(3)
        sage: W is loads(dumps(W))
        True
    """

    def __init__(self, p, base_ring):
        r"""
        Initialisation function.

        EXAMPLES::

            sage: pAdicWeightSpace(17)
            Space of 17-adic weight-characters defined over 17-adic Field with capped relative precision 20
        """
        Parent.__init__(self, base=base_ring, category=Sets())
        p = ZZ(p)
        if not p.is_prime():
            raise ValueError("p must be prime")
        self._p = p
        self._param = Qp(p)((p == 2 and 5) or (p + 1))

    def _repr_(self) -> str:
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: pAdicWeightSpace(17)._repr_()
            'Space of 17-adic weight-characters defined over 17-adic Field with capped relative precision 20'
        """
        return "Space of %s-adic weight-characters defined over %s" % (self.prime(), self.base_ring())

    def __reduce__(self):
        r"""
        Used for pickling.

        EXAMPLES::

            sage: pAdicWeightSpace(3).__reduce__()
            (<function WeightSpace_constructor at ...>, (3, 3-adic Field with capped relative precision 20))
        """
        return (WeightSpace_constructor, (self.prime(), self.base_ring()))

    def _element_constructor_(self, arg1, arg2=None, algebraic=True):
        r"""
        Create an element of this space.

        If ``algebraic = True`` (the default), create a locally algebraic
        character. The arguments should be `(k, \chi)` with `k \in \ZZ`
        and `\chi` a Dirichlet character of `p`-power conductor defined over a
        `p`-adic field; this corresponds to the weight-character `x \mapsto x^k
        \chi(x)`. If `\chi` is omitted, it defaults to the trivial character.

        If ``algebraic = False``, create a general character. The arguments are
        now (t, w) where `t \in \ZZ/(p-1)\ZZ` and `w \in
        \CC_p` with `|w - 1| < 1`. This corresponds to the character
        `\kappa` satisfying `\kappa(\mu) = \mu^t` where `\mu` is a `(p-1)`-st
        root of unity, and `\kappa(1 + p) = w`.

        EXAMPLES::

            sage: W = pAdicWeightSpace(17)
            sage: W(4)
            4
            sage: W(4, DirichletGroup(17, Qp(17)).0)
            (4, 17, [3 + 13*17 + ... + O(17^20)])
            sage: W(1 + O(17^5), 4, algebraic = False)
            [1 + O(17^5), 4]
        """
        if isinstance(arg1, WeightCharacter):
            if arg1.parent() is self:
                return arg1
            elif arg1.parent().prime() == self.prime():
                return self._coerce_in_wtchar(arg1)
            else:
                raise TypeError("Incompatible type!")

        if algebraic:
            return AlgebraicWeight(self, arg1, arg2)
        else:
            return ArbitraryWeight(self, arg1, arg2)

    @cached_method
    def zero(self):
        """
        Return the zero of this weight space.

        EXAMPLES::

            sage: W = pAdicWeightSpace(17)
            sage: W.zero()
            0
        """
        return self(0)

    def prime(self):
        r"""
        Return the prime `p` such that this is a `p`-adic weight space.

        EXAMPLES::

            sage: pAdicWeightSpace(17).prime()
            17
        """
        return self._p

    def base_extend(self, R):
        r"""
        Extend scalars to the ring R.

        There must be a canonical coercion map from the present base ring to R.

        EXAMPLES::

            sage: W = pAdicWeightSpace(3, QQ)
            sage: W.base_extend(Qp(3))
            Space of 3-adic weight-characters
             defined over 3-adic Field with capped relative precision 20
            sage: W.base_extend(IntegerModRing(12))
            Traceback (most recent call last):
            ...
            TypeError: No coercion map from 'Rational Field'
            to 'Ring of integers modulo 12' is defined
        """
        if R.has_coerce_map_from(self.base_ring()):
            return WeightSpace_constructor(self.prime(), R)
        else:
            raise TypeError("No coercion map from '%s' to '%s' is defined" % (self.base_ring(), R))

    def _coerce_map_from_(self, other):
        r"""
        Canonical coercion of ``other`` into ``self``.

        TESTS::

            sage: W1 = pAdicWeightSpace(23, QQ)
            sage: W2 = W1.base_extend(Qp(23))
            sage: w = W1(3)
            sage: W2.coerce(w) # indirect doctest
            3
        """
        return (isinstance(other, WeightSpace_class)
                and other.prime() == self.prime()
                and self.base_ring().has_coerce_map_from(other.base_ring()))

    def _coerce_in_wtchar(self, x):
        r"""
        Convert in a weight-character whose parent is different from ``self``
        (with has the prime, but possibly different base ring).

        EXAMPLES::

            sage: W1 = pAdicWeightSpace(23, Qp(3))
            sage: W2 = pAdicWeightSpace(23, QQ)
            sage: w = W1(3)
            sage: W2._coerce_in_wtchar(w)
            3
        """
        if isinstance(x, AlgebraicWeight):
            return AlgebraicWeight(self, x.k(), x.chi().change_ring(self.base_ring()))
        else:
            return ArbitraryWeight(self, self.base_ring()(x.w()), x.teichmuller_type())


class WeightCharacter(Element):
    r"""
    Abstract base class representing an element of the `p`-adic weight space
    `Hom(\ZZ_p^\times, \CC_p^\times)`.
    """

    # This should probably derive from Morphism or even from
    # AbelianGroupMorphism; but Sage does not know about the abelian group
    # Z_p^*, so Hom(Z_p^*, C_p^*) is a bit beyond it!

    def __init__(self, parent):
        r"""
        Initialisation function.

        EXAMPLES::

            sage: pAdicWeightSpace(17)(0)
            0
        """
        Element.__init__(self, parent)
        self._p = self.parent().prime()

    def base_extend(self, R):
        r"""
        Extend scalars to the base ring R.

        The ring R must have a canonical map from the current base ring.

        EXAMPLES::

            sage: w = pAdicWeightSpace(17, QQ)(3)
            sage: w.base_extend(Qp(17))
            3
        """
        return self.parent().base_extend(R).coerce(self)

    def is_even(self) -> bool:
        r"""
        Return ``True`` if this weight-character sends -1 to +1.

        EXAMPLES::

            sage: pAdicWeightSpace(17)(0).is_even()
            True
            sage: pAdicWeightSpace(17)(11).is_even()
            False
            sage: pAdicWeightSpace(17)(1 + 17 + O(17^20), 3, False).is_even()
            False
            sage: pAdicWeightSpace(17)(1 + 17 + O(17^20), 4, False).is_even()
            True
        """
        return self(-1) != -1

    def pAdicEisensteinSeries(self, ring, prec=20):
        r"""
        Calculate the `q`-expansion of the `p`-adic Eisenstein series of given
        weight-character, normalised so the constant term is 1.

        EXAMPLES::

            sage: kappa = pAdicWeightSpace(3)(3, DirichletGroup(3,QQ).0)
            sage: kappa.pAdicEisensteinSeries(QQ[['q']], 20)
            1 - 9*q + 27*q^2 - 9*q^3 - 117*q^4 + 216*q^5 + 27*q^6 - 450*q^7 + 459*q^8
             - 9*q^9 - 648*q^10 + 1080*q^11 - 117*q^12 - 1530*q^13 + 1350*q^14 + 216*q^15
             - 1845*q^16 + 2592*q^17 + 27*q^18 - 3258*q^19 + O(q^20)
        """
        if not self.is_even():
            raise ValueError("Eisenstein series not defined for odd weight-characters")
        q = ring.gen()
        s = ring(1) + 2*self.one_over_Lvalue() * sum(sum(self(d)/d for d in divisors(n)) * q**n for n in range(1, prec))
        return s.add_bigoh(prec)

    def values_on_gens(self):
        r"""
        If `\kappa` is this character, calculate the values `(\kappa(r), t)`
        where `r` is `1 + p` (or 5 if `p = 2`) and `t` is the unique element of
        `\ZZ/(p-1)\ZZ` such that `\kappa(\mu) = \mu^t` for `\mu`
        a (p-1)st root of unity. (If `p = 2`, we take `t` to be 0 or 1
        according to whether `\kappa` is odd or even.) These two values
        uniquely determine the character `\kappa`.

        EXAMPLES::

            sage: W = pAdicWeightSpace(11); W(2).values_on_gens()
            (1 + 2*11 + 11^2 + O(11^20), 2)
            sage: W(2, DirichletGroup(11, QQ).0).values_on_gens()
            (1 + 2*11 + 11^2 + O(11^20), 7)
            sage: W(1 + 2*11 + O(11^5), 4, algebraic = False).values_on_gens()
            (1 + 2*11 + O(11^5), 4)
        """

        return ( self(self.parent()._param), self.teichmuller_type())

    def is_trivial(self) -> bool:
        r"""
        Return ``True`` if and only if this is the trivial character.

        EXAMPLES::

            sage: pAdicWeightSpace(11)(2).is_trivial()
            False
            sage: pAdicWeightSpace(11)(2, DirichletGroup(11, QQ).0).is_trivial()
            False
            sage: pAdicWeightSpace(11)(0).is_trivial()
            True
        """
        return self.values_on_gens() == (1, 0)

    def _richcmp_(self, other, op) -> bool:
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: W = pAdicWeightSpace(11)
            sage: W(2) == W(3)
            False
            sage: W(2, DirichletGroup(11, QQ).0) == W(2)
            False
            sage: W(2, DirichletGroup(11, QQ).0) == W(144 + O(11^20), 7, False)
            True
        """
        return richcmp(self.values_on_gens(), other.values_on_gens(), op)

    def Lvalue(self):
        r"""
        Return the value of the `p`-adic `L`-function of `\QQ`, which can be
        regarded as a rigid-analytic function on weight space, evaluated at
        this character.

        EXAMPLES::

            sage: W = pAdicWeightSpace(11)
            sage: sage.modular.overconvergent.weightspace.WeightCharacter(W).Lvalue()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def one_over_Lvalue(self):
        r"""
        Return the reciprocal of the `p`-adic `L`-function evaluated at this
        weight-character.

        If the weight-character is odd, then the `L`-function
        is zero, so an error will be raised.

        EXAMPLES::

            sage: pAdicWeightSpace(11)(4).one_over_Lvalue()
            -12/133
            sage: pAdicWeightSpace(11)(3, DirichletGroup(11, QQ).0).one_over_Lvalue()
            -1/6
            sage: pAdicWeightSpace(11)(3).one_over_Lvalue()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
            sage: pAdicWeightSpace(11)(0).one_over_Lvalue()
            0
            sage: type(_)
            <class 'sage.rings.integer.Integer'>
        """
        if self.is_trivial():
            return ZZ(0)
        else:
            return 1 / self.Lvalue()


class AlgebraicWeight(WeightCharacter):
    r"""
    A point in weight space corresponding to a locally algebraic character, of
    the form `x \mapsto \chi(x) x^k` where `k` is an integer and `\chi` is a
    Dirichlet character modulo `p^n` for some `n`.

    TESTS::

        sage: w = pAdicWeightSpace(23)(12, DirichletGroup(23, QQ).0) # exact
        sage: w == loads(dumps(w))
        True
        sage: w = pAdicWeightSpace(23)(12, DirichletGroup(23, Qp(23)).0) # inexact
        sage: w == loads(dumps(w))
        True
        sage: w is loads(dumps(w)) # elements are not globally unique
        False
    """

    def __init__(self, parent, k, chi=None):
        r"""
        Create a locally algebraic weight-character.

        EXAMPLES::

            sage: pAdicWeightSpace(29)(13, DirichletGroup(29, Qp(29)).0)
            (13, 29, [2 + 2*29 + ... + O(29^20)])
        """
        WeightCharacter.__init__(self, parent)
        k = ZZ(k)
        self._k = k
        if chi is None:
            chi = trivial_character(self._p, QQ)
        n = ZZ(chi.conductor())
        if n == 1:
            n = self._p
        if not n.is_power_of(self._p):
            raise ValueError("Character must have %s-power conductor" % self._p)
        self._chi = DirichletGroup(n, chi.base_ring())(chi)

    def __call__(self, x):
        r"""
        Evaluate this character at an element of `\ZZ_p^\times`.

        EXAMPLES:

        Exact answers are returned when this is possible::

            sage: kappa = pAdicWeightSpace(29)(13, DirichletGroup(29, QQ).0)
            sage: kappa(1)
            1
            sage: kappa(0)
            0
            sage: kappa(12)
            -106993205379072
            sage: kappa(-1)
            -1
            sage: kappa(13 + 4*29 + 11*29^2 + O(29^3))
            9 + 21*29 + 27*29^2 + O(29^3)

        When the character chi is defined over a `p`-adic field, the results returned are inexact::

            sage: kappa = pAdicWeightSpace(29)(13, DirichletGroup(29, Qp(29)).0^14)
            sage: kappa(1)
            1 + O(29^20)
            sage: kappa(0)
            0
            sage: kappa(12)
            17 + 11*29 + 7*29^2 + 4*29^3 + 5*29^4 + 2*29^5 + 13*29^6 + 3*29^7 + 18*29^8 + 21*29^9 + 28*29^10 + 28*29^11 + 28*29^12 + 28*29^13 + 28*29^14 + 28*29^15 + 28*29^16 + 28*29^17 + 28*29^18 + 28*29^19 + O(29^20)
            sage: kappa(12) == -106993205379072
            True
            sage: kappa(-1) == -1
            True
            sage: kappa(13 + 4*29 + 11*29^2 + O(29^3))
            9 + 21*29 + 27*29^2 + O(29^3)
        """
        if isinstance(x, pAdicGenericElement):
            if x.parent().prime() != self._p:
                raise TypeError("x must be an integer or a %s-adic integer" % self._p)
            if self._p**(x.precision_absolute()) < self._chi.conductor():
                raise PrecisionError("Precision too low")
            xint = x.lift()
        else:
            xint = x
        if (xint % self._p == 0):
            return 0
        return self._chi(xint) * x**self._k

    def k(self):
        r"""
        If this character is `x \mapsto x^k \chi(x)` for an integer `k` and a
        Dirichlet character `\chi`, return `k`.

        EXAMPLES::

            sage: kappa = pAdicWeightSpace(29)(13, DirichletGroup(29, Qp(29)).0^14)
            sage: kappa.k()
            13
        """
        return self._k

    def chi(self):
        r"""
        If this character is `x \mapsto x^k \chi(x)` for an integer `k` and a
        Dirichlet character `\chi`, return `\chi`.

        EXAMPLES::

            sage: kappa = pAdicWeightSpace(29)(13, DirichletGroup(29, Qp(29)).0^14)
            sage: kappa.chi()
            Dirichlet character modulo 29 of conductor 29
             mapping 2 |--> 28 + 28*29 + 28*29^2 + ... + O(29^20)
        """
        return self._chi

    def __hash__(self):
        r"""
        TESTS::

            sage: w = pAdicWeightSpace(23)(12, DirichletGroup(23, QQ).0)
            sage: hash(w) == hash((12, 23, (-1,)))
            True
        """
        if self._chi.is_trivial():
            return hash(self._k)
        else:
            return hash( (self._k,self._chi.modulus(),self._chi) )

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: pAdicWeightSpace(17)(2)._repr_()
            '2'
            sage: pAdicWeightSpace(17)(2, DirichletGroup(17, QQ).0)._repr_()
            '(2, 17, [-1])'
            sage: pAdicWeightSpace(17)(2, DirichletGroup(17, QQ).0^2)._repr_()
            '2'
        """
        if self._chi.is_trivial():
            return "%s" % self._k
        else:
            return "(%s, %s, %s)" % (self._k, self._chi.modulus(), self._chi._repr_short_())

    def teichmuller_type(self):
        r"""
        Return the Teichmuller type of this weight-character `\kappa`.

        This is the unique `t \in \ZZ/(p-1)\ZZ` such that `\kappa(\mu)
        = \mu^t` for `\mu` a `(p-1)`-st root of 1.

        For `p = 2` this does not make sense, but we still want the Teichmuller
        type to correspond to the index of the component of weight space in
        which `\kappa` lies, so we return 1 if `\kappa` is odd and 0 otherwise.

        EXAMPLES::

            sage: pAdicWeightSpace(11)(2, DirichletGroup(11,QQ).0).teichmuller_type()
            7
            sage: pAdicWeightSpace(29)(13, DirichletGroup(29, Qp(29)).0).teichmuller_type()
            14
            sage: pAdicWeightSpace(2)(3, DirichletGroup(4,QQ).0).teichmuller_type()
            0
        """
        # Special case p == 2
        if self._p == 2:
            if self.is_even():
                return IntegerModRing(2).zero()
            else:
                return IntegerModRing(2).one()
        m = IntegerModRing(self._p).multiplicative_generator()
        x = [y for y in IntegerModRing(self._chi.modulus()) if y == m and y**(self._p - 1) == 1]
        if len(x) != 1:
            raise ArithmeticError
        x = x[0]
        f = IntegerModRing(self._p)(self._chi(x)).log(m)
        return IntegerModRing(self._p - 1)(self._k + f)

    def Lvalue(self):
        r"""
        Return the value of the `p`-adic `L`-function of `\QQ` evaluated at
        this weight-character.

        If the character is `x \mapsto x^k \chi(x)`
        where `k > 0` and `\chi` has conductor a power of `p`, this is an
        element of the number field generated by the values of `\chi`, equal to
        the value of the complex `L`-function `L(1-k, \chi)`. If `\chi` is
        trivial, it is equal to `(1 - p^{k-1})\zeta(1-k)`.

        At present this is not implemented in any other cases, except the
        trivial character (for which the value is `\infty`).

        .. TODO::

            Implement this more generally using the Amice transform
            machinery in
            sage/schemes/elliptic_curves/padic_lseries.py, which
            should clearly be factored out into a separate class.

        EXAMPLES::

            sage: pAdicWeightSpace(7)(4).Lvalue() == (1 - 7^3)*zeta__exact(-3)
            True
            sage: pAdicWeightSpace(7)(5, DirichletGroup(7, Qp(7)).0^4).Lvalue()
            0
            sage: pAdicWeightSpace(7)(6, DirichletGroup(7, Qp(7)).0^4).Lvalue()
            1 + 2*7 + 7^2 + 3*7^3 + 3*7^5 + 4*7^6 + 2*7^7 + 5*7^8 + 2*7^9 + 3*7^10 + 6*7^11
             + 2*7^12 + 3*7^13 + 5*7^14 + 6*7^15 + 5*7^16 + 3*7^17 + 6*7^18 + O(7^19)
        """
        if self._k > 0:
            return -self._chi.bernoulli(self._k) / self._k
        if self.is_trivial():
            return Infinity
        else:
            raise NotImplementedError("Don't know how to compute value of this L-function")


class ArbitraryWeight(WeightCharacter):

    def __init__(self, parent, w, t):
        r"""
        Create the element of `p`-adic weight space in the given component
        mapping 1 + p to w.

        Here w must be an element of a `p`-adic field, with finite
        precision.

        EXAMPLES::

            sage: pAdicWeightSpace(17)(1 + 17^2 + O(17^3), 11, False)
            [1 + 17^2 + O(17^3), 11]
        """
        WeightCharacter.__init__(self, parent)

        self.t = ZZ(t) % (self._p > 2 and (self._p - 1) or 2)
                # do we store w precisely?
        if (w - 1).valuation() <= 0:
            raise ValueError("Must send generator to something nearer 1")
        self.w = w

    def _repr_(self):
        r"""
        String representation of this character.

        EXAMPLES::

            sage: pAdicWeightSpace(97)(1 + 2*97 + O(97^20), 12, False)._repr_()
            '[1 + 2*97 + O(97^20), 12]'
        """
        return "[%s, %s]" % (self.w, self.t)

    def __call__(self, x):
        r"""
        Evaluate this character at an element of `\ZZ_p^\times`.

        EXAMPLES::

            sage: kappa = pAdicWeightSpace(23)(1 + 23^2 + O(23^20), 4, False)
            sage: kappa(2)
            16 + 7*23 + 7*23^2 + 16*23^3 + 23^4 + 20*23^5 + 15*23^7 + 11*23^8 + 12*23^9 + 8*23^10 + 22*23^11 + 16*23^12 + 13*23^13 + 4*23^14 + 19*23^15 + 6*23^16 + 7*23^17 + 11*23^19 + O(23^20)
            sage: kappa(-1)
            1 + O(23^20)
            sage: kappa(23)
            0
            sage: kappa(2 + 2*23 + 11*23^2 + O(23^3))
            16 + 7*23 + O(23^3)
        """
        from sage.misc.verbose import verbose

        if not isinstance(x, pAdicGenericElement):
            x = Qp(self._p)(x)
        if x.valuation() != 0:
            return 0

        teich = x.parent().teichmuller(x, x.precision_absolute())
        xx = x / teich
        if (xx - 1).valuation() <= 0:
            raise ArithmeticError
        verbose("Normalised element is %s" % xx)

        e = xx.log() / self.parent()._param.log()
        verbose("Exponent is %s" % e)

        return teich**(self.t) * (self.w.log() * e).exp()

    def teichmuller_type(self):
        r"""
        Return the Teichmuller type of this weight-character `\kappa`.

        This is
        the unique `t \in \ZZ/(p-1)\ZZ` such that `\kappa(\mu) =
        \mu^t` for \mu a `(p-1)`-st root of 1.

        For `p = 2` this does not make sense, but we still want the Teichmuller
        type to correspond to the index of the component of weight space in
        which `\kappa` lies, so we return 1 if `\kappa` is odd and 0 otherwise.

        EXAMPLES::

            sage: pAdicWeightSpace(17)(1 + 3*17 + 2*17^2 + O(17^3), 8, False).teichmuller_type()
            8
            sage: pAdicWeightSpace(2)(1 + 2 + O(2^2), 1, False).teichmuller_type()
            1
        """
        return self.t
