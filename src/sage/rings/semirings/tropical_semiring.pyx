r"""
Tropical Semirings

AUTHORS:

- Travis Scrimshaw (2013-04-28) - Initial version
"""
# ****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element cimport Element, ModuleElement
from sage.structure.richcmp cimport rich_to_bool
from sage.categories.semirings import Semirings
from sage.categories.map cimport Map
from sage.rings.integer_ring import ZZ


cdef class TropicalSemiringElement(Element):
    r"""
    An element in the tropical semiring over an ordered additive semigroup
    `R`. Either in `R` or `\infty`. The operators `+, \cdot` are defined as
    the tropical operators `\oplus, \odot` respectively.
    """
    cdef ModuleElement _val

    cdef TropicalSemiringElement _new(self):
        """
        Return a new tropical semiring element with parent ``self``.
        """
        cdef TropicalSemiringElement x
        x = TropicalSemiringElement.__new__(TropicalSemiringElement)
        x._parent = self._parent
        x._val = self._val
        return x

    def __init__(self, parent, ModuleElement val=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: elt = T(2)
            sage: TestSuite(elt).run()
        """
        Element.__init__(self, parent)
        self._val = val

    def __reduce__(self):
        """
        Used in pickling tropical semiring elements.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: elt = T(2)
            sage: elt.__reduce__()
            (<class 'sage.rings.semirings.tropical_semiring.TropicalSemiringElement'>,
             (Tropical semiring over Rational Field, 2))
        """
        return (TropicalSemiringElement, (self.parent(), self._val))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: T(2)
            2
            sage: T.infinity()
            +infinity
            sage: T = TropicalSemiring(QQ, False)
            sage: T.infinity()
            -infinity
        """
        if self._val is None:
            if self.parent()._use_min:
                return "+infinity"
            return "-infinity"
        return repr(self._val)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: T(2)
            2
            sage: latex(T.infinity())
            \infty
            sage: T = TropicalSemiring(QQ, False)
            sage: latex(T.infinity())
            -\infty
        """
        if self._val is None:
            if self.parent()._use_min:
                return "\\infty"
            return "-\\infty"
        return repr(self._val)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: hash(T(2)) == hash(T(2))
            True
            sage: hash(T.infinity()) == hash(T.infinity())
            True
        """
        return hash(self._val)

    # Comparisons
    cpdef _richcmp_(left, right, int op):
        r"""
        Return the standard comparison of ``left`` and ``right``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: T(2) == T(2)
            True
            sage: T.infinity() == T.infinity()
            True
            sage: T(2) != T(4)
            True
            sage: T.infinity() != T.infinity()
            False
            sage: T(2) < T(4)
            True
            sage: T(2) < T.infinity()
            True
            sage: T.infinity() < T.infinity()
            False
            sage: T(2) <= T(4)
            True
            sage: T.infinity() <= T.infinity()
            True
            sage: T(2) > T(4)
            False
            sage: T(2) > T.infinity()
            False
            sage: T.infinity() > T.infinity()
            False
            sage: T(2) >= T(4)
            False
            sage: T.infinity() >= T.infinity()
            True

        Using the `\max` definition::

            sage: T = TropicalSemiring(QQ, False)
            sage: T(2) == T(2)
            True
            sage: T(2) != T(4)
            True
            sage: T(2) < T(4)
            True
            sage: T(2) > T(4)
            False
            sage: T.infinity() == T.infinity()
            True
            sage: T(4) <= T.infinity()
            False
        """
        cdef TropicalSemiringElement self, x
        self = left
        x = right

        if self._val is None:
            if x._val is None:
                return rich_to_bool(op, 0)
            if self.parent()._use_min:
                return rich_to_bool(op, 1)
            return rich_to_bool(op, -1)

        if x._val is None:
            if self.parent()._use_min:
                return rich_to_bool(op, -1)
            return rich_to_bool(op, 1)

        if self._val < x._val:
            return rich_to_bool(op, -1)
        if self._val > x._val:
            return rich_to_bool(op, 1)
        return rich_to_bool(op, 0)

    cpdef _add_(left, right):
        """
        Add ``left`` to ``right``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: T(2) + T(4)
            2
            sage: T(2) + T.infinity()
            2
            sage: T.infinity() + T(2)
            2
            sage: T = TropicalSemiring(QQ, False)
            sage: T(2) + T(4)
            4
            sage: T(2) + T.infinity()
            2
            sage: T.infinity() + T(2)
            2

        TESTS:

        Check some additions with trivial coercion of ``int(0)``::

            sage: T = TropicalSemiring(QQ)
            sage: T(1) + int(0)
            0
            sage: T(-1) + int(0)
            -1
            sage: int(0) + T(1)
            0
            sage: int(0) + T(-1)
            -1
        """
        cdef TropicalSemiringElement self, rhs
        self = left
        rhs = right
        if self._val is None:
            return rhs
        if rhs._val is None:
            return self
        cdef TropicalSemiringElement x
        x = self._new()
        if self.parent()._use_min:
            x._val = min(self._val, rhs._val)
        else:
            x._val = max(self._val, rhs._val)
        return x

    def __neg__(self):
        r"""
        Return the additive inverse, which only exists for `\infty`.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: -T.infinity()
            +infinity
            sage: -T(2)
            Traceback (most recent call last):
            ...
            ArithmeticError: cannot negate any non-infinite element
        """
        if self._val is None:
            return self
        raise ArithmeticError("cannot negate any non-infinite element")

    cpdef _mul_(left, right):
        """
        Multiply ``left`` and ``right``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: T(2) * T(4)
            6
            sage: T(2) * T.infinity()
            +infinity
            sage: T(-2) * T(1)
            -1
        """
        cdef TropicalSemiringElement self, rhs
        self = left
        rhs = right
        if self._val is None:
            return self
        if rhs._val is None:
            return rhs
        cdef TropicalSemiringElement x
        x = self._new()
        x._val = self._val + rhs._val
        return x

    cpdef _div_(left, right):
        """
        Divide ``left`` by ``right``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: T(2) / T(4)
            -2
            sage: T.infinity() / T(2)
            +infinity
        """
        cdef TropicalSemiringElement self, rhs
        self = left
        rhs = right

        if rhs._val is None:
            raise ZeroDivisionError("Tropical division by infinity")
        if self._val is None:
            return self
        cdef TropicalSemiringElement x
        x = self._new()
        x._val = self._val - rhs._val
        return x

    def __invert__(self):
        """
        Return the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: ~T(2)
            -2
        """
        if self.is_one():
            return self
        if self._val is None:
            raise ZeroDivisionError("Tropical division by infinity")
        cdef TropicalSemiringElement x
        x = self._new()
        x._val = -self._val
        return x

    def __pow__(base, exp, dummy):
        """
        Return ``self`` to ``exp``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: elt = T(2)
            sage: elt**3
            6
            sage: elt**-2
            -4
            sage: elt**(3/7)
            6/7
            sage: elt**0
            0

            sage: elt = T.infinity()
            sage: elt**0
            0
            sage: elt**(1/2)
            +infinity
            sage: elt*33
            +infinity
        """
        cdef TropicalSemiringElement self, x
        self = base
        if self._val is None:
            if exp > 0:
                return self
            elif exp == 0:
                return self.parent().multiplicative_identity()
            raise ZeroDivisionError("Tropical division by infinity")
        x = self._new()
        x._val = exp*self._val
        return x

    def multiplicative_order(self):
        """
        Return the multiplicative order of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: T.multiplicative_identity().multiplicative_order()
            1
            sage: T.additive_identity().multiplicative_order()
            +Infinity
        """
        if self.is_one():
            return ZZ.one()
        from sage.rings.infinity import infinity
        return infinity

    cpdef ModuleElement lift(self):
        """
        Return the value of ``self`` lifted to the base.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: elt = T(2)
            sage: elt.lift()
            2
            sage: elt.lift().parent() is QQ
            True
            sage: T.additive_identity().lift().parent()
            The Infinity Ring
        """
        if self._val is None:
            from sage.rings.infinity import infinity
            return infinity
        return self._val


class TropicalSemiring(Parent, UniqueRepresentation):
    r"""
    The tropical semiring.

    Given an ordered additive semigroup `R`, we define the tropical
    semiring `T = R \cup \{+\infty\}` by defining tropical addition
    and multiplication as follows:

    .. MATH::

        a \oplus b = \min(a, b), \quad \quad a \odot b = a + b.

    In particular, note that there are no (tropical) additive inverses
    (except for `\infty`), and every element in `R` has a (tropical)
    multiplicative inverse.

    There is an alternative definition where we define `T = R \cup \{-\infty\}`
    and alter tropical addition to be defined by

    .. MATH::

        a \oplus b = \max(a, b).

    To use the `\max` definition, set the argument ``use_min = False``.

    .. WARNING::

        :meth:`zero` and :meth:`one` refer to the tropical additive
        and multiplicative identities respectively. These are **not** the
        same as calling ``T(0)`` and ``T(1)`` respectively as these are **not**
        the tropical additive and multiplicative identities respectively.

        Specifically do not use ``sum(...)`` as this converts `0` to `0` as
        a tropical element, which is not the same as :meth:`zero`. Instead
        use the ``sum`` method of the tropical semiring::

            sage: T = TropicalSemiring(QQ)

            sage: sum([T(1), T(2)]) # This is wrong
            0
            sage: T.sum([T(1), T(2)]) # This is correct
            1

        Be careful about using code that has not been checked for tropical
        safety.

    INPUT:

    - ``base`` -- the base ordered additive semigroup `R`
    - ``use_min`` -- boolean (default: ``True``); if ``True``, then the semiring uses
      `a \oplus b = \min(a, b)`. Otherwise uses `a \oplus b = \max(a, b)`.

    EXAMPLES::

        sage: T = TropicalSemiring(QQ)
        sage: elt = T(2); elt
        2

    Recall that tropical addition is the minimum of two elements::

        sage: T(3) + T(5)
        3

    Tropical multiplication is the addition of two elements::

        sage: T(2) * T(3)
        5
        sage: T(0) * T(-2)
        -2

    We can also do tropical division and arbitrary tropical exponentiation::

        sage: T(2) / T(1)
        1
        sage: T(2)^(-3/7)
        -6/7

    Note that "zero" and "one" are the additive and multiplicative
    identities of the tropical semiring. In general, they are **not**
    the elements `0` and `1` of `R`, respectively, even if such elements
    exist (e.g., for `R = \ZZ`), but instead the (tropical) additive and
    multiplicative identities `+\infty` and `0` respectively::

        sage: T.zero() + T(3) == T(3)
        True
        sage: T.one() * T(3) == T(3)
        True
        sage: T.zero() == T(0)
        False
        sage: T.one() == T(1)
        False
    """
    def __init__(self, base, use_min=True):
        r"""
        Initialize ``self``.

        TESTS::

            sage: T = TropicalSemiring(QQ); T
            Tropical semiring over Rational Field
            sage: TestSuite(T).run()
        """
        self._use_min = use_min
        self._names = ('x', 'infty')
        Parent.__init__(self, base=base, category=Semirings())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: TropicalSemiring(QQ)
            Tropical semiring over Rational Field
        """
        return "Tropical semiring over {}".format(self.base())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(TropicalSemiring(QQ))
            \Bold{T}\left( \Bold{Q} \right)
        """
        return "\\Bold{T}\\left( " + self.base()._latex_() + " \\right)"

    def _coerce_map_from_(self, S):
        """
        Canonical coercion into ``self`` from ``S``.

        The only objects that canonically coerce to a tropical semiring are
        tropical semirings whose base rings have a coercion.

        EXAMPLES::

            sage: TQ = TropicalSemiring(QQ)
            sage: TQ.has_coerce_map_from(TQ)
            True
            sage: TQ.has_coerce_map_from(TropicalSemiring(ZZ))
            True

            sage: # needs sage.rings.real_mpfr
            sage: TR = TropicalSemiring(RR)
            sage: T60 = TropicalSemiring(RealField(60))
            sage: TR.has_coerce_map_from(T60)
            True
            sage: TR.has_coerce_map_from(TR)
            True
            sage: TR.has_coerce_map_from(TQ)
            True
            sage: TR.has_coerce_map_from(float)
            False
            sage: TR.has_coerce_map_from(RR)
            False
            sage: TR.has_coerce_map_from(QQ)
            False
            sage: TR.coerce_map_from(T60)(T60(2))
            2.00000000000000
            sage: TR.coerce(T60(3.4))
            3.40000000000000
            sage: TR.coerce(T60.infinity())
            +infinity
            sage: TQ.coerce(TR(3.4))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Tropical semiring over
             Real Field with 53 bits of precision to Tropical semiring over Rational Field
        """
        if isinstance(S, TropicalSemiring) and self._use_min == S._use_min \
                and self.base().has_coerce_map_from(S.base()):
            return TropicalToTropical(S, self)

    def _element_constructor_(self, val):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: T(2)
            2
        """
        if val is not None:
            val = self.base()(val)
        return self.element_class(self, val)

    Element = TropicalSemiringElement

    @cached_method
    def zero(self):
        r"""
        Return the (tropical) additive identity element `+\infty`.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: T.zero()
            +infinity
        """
        return self.element_class(self, None)

    infinity = zero
    additive_identity = zero

    @cached_method
    def one(self):
        """
        Return the (tropical) multiplicative identity element `0`.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: T.one()
            0
        """
        return self.element_class(self, self.base().zero())

    multiplicative_identity = one

    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: T.gens()
            (1, +infinity)
        """
        return (self.element_class(self, self.base().one()), self.infinity())


cdef class TropicalToTropical(Map):
    """
    Map from the tropical semiring to itself (possibly with different bases).
    Used in coercion.
    """
    cpdef TropicalSemiringElement _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.semirings.tropical_semiring import TropicalToTropical
            sage: TZ = TropicalSemiring(ZZ)
            sage: TQ = TropicalSemiring(QQ)
            sage: f = TropicalToTropical(TZ, TQ)
            sage: a = TZ(2)
            sage: f(a)
            2
            sage: f(TZ.infinity())
            +infinity
        """
        return self.codomain()((<TropicalSemiringElement>x)._val)
