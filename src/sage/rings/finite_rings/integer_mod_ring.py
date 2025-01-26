r"""
Ring `\ZZ/n\ZZ` of integers modulo `n`

EXAMPLES::

    sage: R = Integers(97)
    sage: a = R(5)
    sage: a**100000000000000000000000000000000000000000000000000000000000000
    61

This example illustrates the relation between
`\ZZ/p\ZZ` and `\GF{p}`. In
particular, there is a canonical map to `\GF{p}`, but not in
the other direction.

::

    sage: r = Integers(7)
    sage: s = GF(7)
    sage: r.has_coerce_map_from(s)
    False
    sage: s.has_coerce_map_from(r)
    True
    sage: s(1) + r(1)
    2
    sage: parent(s(1) + r(1))
    Finite Field of size 7
    sage: parent(r(1) + s(1))
    Finite Field of size 7

We list the elements of `\ZZ/3\ZZ`::

    sage: R = Integers(3)
    sage: list(R)
    [0, 1, 2]

AUTHORS:

- William Stein (initial code)

- David Joyner (2005-12-22): most examples

- Robert Bradshaw (2006-08-24): convert to SageX (Cython)

- William Stein (2007-04-29): square_roots_of_one

- Simon King (2011-04-21): allow to prescribe a category

- Simon King (2013-09): Only allow to prescribe the category of fields

- Kyle Hofmann (2024-02): New implementation of root-finding
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


import sage.misc.prandom as random

from sage.arith.misc import factor
from sage.arith.misc import primitive_root
from sage.arith.misc import CRT_basis
from sage.rings.ring import Field, CommutativeRing
from sage.misc.mrange import cartesian_product_iterator
import sage.rings.abc
from sage.rings.finite_rings import integer_mod
import sage.rings.integer as integer
import sage.rings.integer_ring as integer_ring
import sage.rings.quotient_ring as quotient_ring

try:
    from sage.libs.pari.all import pari, PariError
except ImportError:
    class PariError(Exception):
        pass

from sage.misc.cachefunc import cached_method

from sage.structure.factory import UniqueFactory
from sage.structure.richcmp import richcmp, richcmp_method

from sage.interfaces.abc import GapElement


class IntegerModFactory(UniqueFactory):
    r"""
    Return the quotient ring `\ZZ / n\ZZ`.

    INPUT:

    - ``order`` -- integer (default: 0); positive or negative
    - ``is_field`` -- boolean (default: ``False``); assert that the order is
      prime and hence the quotient ring belongs to the category of fields
    - ``category`` -- (optional) the category that the quotient ring belongs to

    .. NOTE::

        The optional argument ``is_field`` is not part of the cache key.
        Hence, this factory will create precisely one instance of `\ZZ /
        n\ZZ`.  However, if ``is_field`` is true, then a previously created
        instance of the quotient ring will be updated to be in the category of
        fields.

        **Use with care!** Erroneously putting `\ZZ / n\ZZ` into the category
        of fields may have consequences that can compromise a whole Sage
        session, so that a restart will be needed.

    EXAMPLES::

        sage: IntegerModRing(15)
        Ring of integers modulo 15
        sage: IntegerModRing(7)
        Ring of integers modulo 7
        sage: IntegerModRing(-100)
        Ring of integers modulo 100

    Note that you can also use ``Integers``, which is a
    synonym for ``IntegerModRing``.

    ::

        sage: Integers(18)
        Ring of integers modulo 18
        sage: Integers() is Integers(0) is ZZ
        True

    .. NOTE::

        Testing whether a quotient ring `\ZZ / n\ZZ` is a field can of
        course be very costly. By default, it is not tested whether `n`
        is prime or not, in contrast to
        :func:`~sage.rings.finite_rings.finite_field_constructor.GF`. If the user
        is sure that the modulus is prime and wants to avoid a primality
        test, (s)he can provide ``category=Fields()`` when constructing
        the quotient ring, and then the result will behave like a field.
        If the category is not provided during initialisation, and it is
        found out later that the ring is in fact a field, then the category
        will be changed at runtime, having the same effect as providing
        ``Fields()`` during initialisation.

    EXAMPLES::

        sage: R = IntegerModRing(5)
        sage: R.category()
        Join of Category of finite commutative rings
            and Category of subquotients of monoids
            and Category of quotients of semigroups
            and Category of finite enumerated sets
        sage: R in Fields()
        True
        sage: R.category()
        Join of Category of finite enumerated fields
            and Category of subquotients of monoids
            and Category of quotients of semigroups
        sage: S = IntegerModRing(5, is_field=True)
        sage: S is R
        True

    .. WARNING::

        If the optional argument ``is_field`` was used by mistake, there is
        currently no way to revert its impact, even though
        :meth:`IntegerModRing_generic.is_field` with the optional argument
        ``proof=True`` would return the correct answer.  So, prescribe
        ``is_field=True`` only if you know what your are doing!

    EXAMPLES::

        sage: R = IntegerModRing(33, is_field=True)
        sage: R in Fields()
        True
        sage: R.is_field()
        True

    If the optional argument `proof=True` is provided, primality is tested and
    the mistaken category assignment is reported::

        sage: R.is_field(proof=True)
        Traceback (most recent call last):
        ...
        ValueError: THIS SAGE SESSION MIGHT BE SERIOUSLY COMPROMISED!
        The order 33 is not prime, but this ring has been put
        into the category of fields. This may already have consequences
        in other parts of Sage. Either it was a mistake of the user,
        or a probabilistic primality test has failed.
        In the latter case, please inform the developers.

    However, the mistaken assignment is not automatically corrected::

        sage: R in Fields()
        True

    To avoid side-effects of this test on other tests, we clear the cache of
    the ring factory::

        sage: IntegerModRing._cache.clear()
    """
    def get_object(self, version, key, extra_args):
        out = super().get_object(version, key, extra_args)
        category = extra_args.get('category', None)
        if category is not None:
            out._refine_category_(category)
            out._factory_data[3]['category'] = category
        return out

    def create_key_and_extra_args(self, order=0, is_field=False, category=None):
        """
        An integer mod ring is specified uniquely by its order.

        EXAMPLES::

            sage: Zmod.create_key_and_extra_args(7)
            (7, {})
            sage: Zmod.create_key_and_extra_args(7, True)
            (7, {'category': Category of fields})
        """
        if is_field:
            from sage.categories.fields import Fields
            return order, {'category': Fields()}
        return order, {}

    def create_object(self, version, order, **kwds):
        """
        EXAMPLES::

            sage: R = Integers(10)
            sage: TestSuite(R).run() # indirect doctest
        """
        if isinstance(order, tuple):
            # this is for unpickling old data
            order, category = order
            kwds.setdefault('category', category)
        if order < 0:
            order = -order
        if order == 0:
            return integer_ring.IntegerRing(**kwds)
        else:
            return IntegerModRing_generic(order, **kwds)


Zmod = Integers = IntegerModRing = IntegerModFactory("IntegerModRing")


from sage.categories.commutative_rings import CommutativeRings
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.category import JoinCategory
default_category = JoinCategory((CommutativeRings(), FiniteEnumeratedSets()))
ZZ = integer_ring.IntegerRing()


def _unit_gens_primepowercase(p, r):
    r"""
    Return a list of generators for `(\ZZ/p^r\ZZ)^*` and their orders.

    EXAMPLES::

        sage: from sage.rings.finite_rings.integer_mod_ring import _unit_gens_primepowercase
        sage: _unit_gens_primepowercase(2, 3)
        [(7, 2), (5, 2)]
        sage: _unit_gens_primepowercase(17, 1)                                          # needs sage.libs.pari
        [(3, 16)]
        sage: _unit_gens_primepowercase(3, 3)                                           # needs sage.libs.pari
        [(2, 18)]
    """
    pr = p**r
    if p == 2:
        if r == 1:
            return []
        if r == 2:
            return [(integer_mod.Mod(3, 4), integer.Integer(2))]
        return [(integer_mod.Mod(-1, pr), integer.Integer(2)),
                (integer_mod.Mod(5, pr), integer.Integer(2**(r - 2)))]

    # odd prime
    return [(integer_mod.Mod(primitive_root(pr, check=False), pr),
             integer.Integer(p**(r - 1) * (p - 1)))]


@richcmp_method
class IntegerModRing_generic(quotient_ring.QuotientRing_generic, sage.rings.abc.IntegerModRing):
    """
    The ring of integers modulo `N`.

    INPUT:

    - ``order`` -- integer

    - ``category`` -- a subcategory of ``CommutativeRings()`` (the default)

    OUTPUT: the ring of integers modulo `N`

    EXAMPLES:

    First we compute with integers modulo `29`.

    ::

        sage: FF = IntegerModRing(29)
        sage: FF
        Ring of integers modulo 29
        sage: FF.category()
        Join of Category of finite commutative rings
            and Category of subquotients of monoids
            and Category of quotients of semigroups
            and Category of finite enumerated sets
        sage: FF.is_field()
        True
        sage: FF.characteristic()
        29
        sage: FF.order()
        29

        sage: # needs sage.groups
        sage: gens = FF.unit_gens()
        sage: a = gens[0]
        sage: a
        2
        sage: a.is_square()
        False
        sage: def pow(i): return a**i
        sage: [pow(i) for i in range(16)]
        [1, 2, 4, 8, 16, 3, 6, 12, 24, 19, 9, 18, 7, 14, 28, 27]
        sage: TestSuite(FF).run()

    We have seen above that an integer mod ring is, by default, not
    initialised as an object in the category of fields. However, one
    can force it to be. Moreover, testing containment in the category
    of fields my re-initialise the category of the integer mod ring::

        sage: F19 = IntegerModRing(19, is_field=True)
        sage: F19.category().is_subcategory(Fields())
        True
        sage: F23 = IntegerModRing(23)
        sage: F23.category().is_subcategory(Fields())
        False
        sage: F23 in Fields()
        True
        sage: F23.category().is_subcategory(Fields())
        True
        sage: TestSuite(F19).run()
        sage: TestSuite(F23).run()

    By :issue:`15229`, there is a unique instance of the
    integral quotient ring of a given order. Using the
    :func:`IntegerModRing` factory twice, and using
    ``is_field=True`` the second time, will update the
    category of the unique instance::

        sage: F31a = IntegerModRing(31)
        sage: F31a.category().is_subcategory(Fields())
        False
        sage: F31b = IntegerModRing(31, is_field=True)
        sage: F31a is F31b
        True
        sage: F31a.category().is_subcategory(Fields())
        True

    Next we compute with the integers modulo `16`.

    ::

        sage: Z16 = IntegerModRing(16)
        sage: Z16.category()
        Join of Category of finite commutative rings
            and Category of subquotients of monoids
            and Category of quotients of semigroups
            and Category of finite enumerated sets
        sage: Z16.is_field()
        False
        sage: Z16.order()
        16
        sage: Z16.characteristic()
        16

        sage: # needs sage.groups
        sage: gens = Z16.unit_gens()
        sage: gens
        (15, 5)
        sage: a = gens[0]
        sage: b = gens[1]
        sage: def powa(i): return a**i
        sage: def powb(i): return b**i
        sage: gp_exp = FF.unit_group_exponent()
        sage: gp_exp
        28
        sage: [powa(i) for i in range(15)]
        [1, 15, 1, 15, 1, 15, 1, 15, 1, 15, 1, 15, 1, 15, 1]
        sage: [powb(i) for i in range(15)]
        [1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9, 13, 1, 5, 9]
        sage: a.multiplicative_order()
        2
        sage: b.multiplicative_order()
        4

        sage: TestSuite(Z16).run()

    Saving and loading::

        sage: R = Integers(100000)
        sage: TestSuite(R).run()  # long time (17s on sage.math, 2011)

    Testing ideals and quotients::

        sage: Z10 = Integers(10)
        sage: I = Z10.principal_ideal(0)
        sage: Z10.quotient(I) == Z10
        True
        sage: I = Z10.principal_ideal(2)
        sage: Z10.quotient(I) == Z10
        False
        sage: I.is_prime()
        True

    ::

        sage: R = IntegerModRing(97)
        sage: a = R(5)
        sage: a**(10^62)
        61
    """
    def __init__(self, order, cache=None, category=None):
        """
        Create with the command ``IntegerModRing(order)``.

        TESTS::

            sage: FF = IntegerModRing(29)
            sage: TestSuite(FF).run()
            sage: F19 = IntegerModRing(19, is_field=True)
            sage: TestSuite(F19).run()
            sage: F23 = IntegerModRing(23)
            sage: F23 in Fields()
            True
            sage: TestSuite(F23).run()
            sage: Z16 = IntegerModRing(16)
            sage: TestSuite(Z16).run()
            sage: R = Integers(100000)
            sage: TestSuite(R).run()  # long time (17s on sage.math, 2011)

            sage: R = IntegerModRing(18)
            sage: R.is_finite()
            True
        """
        order = ZZ(order)
        if order <= 0:
            raise ZeroDivisionError("order must be positive")
        self.__order = order
        self._pyx_order = integer_mod.NativeIntStruct(order)
        global default_category
        if category is None:
            category = default_category
        else:
            # If the category is given, e.g., as Fields(), then we still
            # know that the result will also live in default_category.
            # Hence, we use the join of the default and the given category.
            category = category.join([category, default_category])
        # Give the generator a 'name' to make quotients work.  The
        # name 'x' is used because it's also used for the ring of
        # integers: see the __init__ method for IntegerRing_class in
        # sage/rings/integer_ring.pyx.
        quotient_ring.QuotientRing_generic.__init__(self, ZZ, ZZ.ideal(order),
                                                    names=('x',),
                                                    category=category)
        # We want that the ring is its own base ring.
        self._base = self
        if cache is None:
            cache = order < 500
        if cache:
            self._precompute_table()
        self._zero_element = integer_mod.IntegerMod(self, 0)
        self._one_element = integer_mod.IntegerMod(self, 1)

    def _macaulay2_init_(self, macaulay2=None):
        """
        EXAMPLES::

            sage: macaulay2(Integers(7))  # optional - macaulay2
            ZZ
            --
             7

        ::

            sage: macaulay2(Integers(10)) # optional - macaulay2
            Traceback (most recent call last):
            ...
            TypeError: Error evaluating Macaulay2 code.
            IN:...
            OUT:...error: ZZ/n not implemented yet for composite n...
        """
        return "ZZ/{}".format(self.order())

    def _axiom_init_(self):
        """
        Return a string representation of ``self`` in (Pan)Axiom.

        EXAMPLES::

            sage: Z7 = Integers(7)
            sage: Z7._axiom_init_()
            'IntegerMod(7)'

            sage: axiom(Z7)  #optional - axiom
            IntegerMod 7

            sage: fricas(Z7) #optional - fricas
            IntegerMod(7)
        """
        return 'IntegerMod({})'.format(self.order())

    _fricas_init_ = _axiom_init_

    def krull_dimension(self):
        """
        Return the Krull dimension of ``self``.

        EXAMPLES::

            sage: Integers(18).krull_dimension()
            0
        """
        return integer.Integer(0)

    def is_noetherian(self):
        """
        Check if ``self`` is a Noetherian ring.

        EXAMPLES::

            sage: Integers(8).is_noetherian()
            True
        """
        return True

    def extension(self, poly, name=None, names=None, **kwds):
        """
        Return an algebraic extension of ``self``. See
        :meth:`sage.rings.ring.CommutativeRing.extension()` for more
        information.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: Integers(8).extension(t^2 - 3)
            Univariate Quotient Polynomial Ring in t
             over Ring of integers modulo 8 with modulus t^2 + 5
        """
        if self.modulus() == 1:
            return self

        from sage.rings.ring import CommutativeRing
        return CommutativeRing.extension(self, poly, name, names, **kwds)

    @cached_method
    def is_prime_field(self):
        """
        Return ``True`` if the order is prime.

        EXAMPLES::

            sage: Zmod(7).is_prime_field()
            True
            sage: Zmod(8).is_prime_field()
            False
        """
        return self.__order.is_prime()

    def _precompute_table(self):
        """
        Compute a table of elements so that elements are unique.

        EXAMPLES::

            sage: R = Zmod(500); R._precompute_table()
            sage: R(7) + R(13) is R(3) + R(17)
            True
        """
        self._pyx_order.precompute_table(self)

    def list_of_elements_of_multiplicative_group(self):
        """
        Return a list of all invertible elements, as python ints.

        EXAMPLES::

            sage: R = Zmod(12)
            sage: L = R.list_of_elements_of_multiplicative_group(); L
            [1, 5, 7, 11]
            sage: type(L[0])
            <... 'int'>
            sage: Zmod(1).list_of_elements_of_multiplicative_group()
            [0]
        """
        import sage.rings.fast_arith as a
        if self.__order <= 46340:   # todo: don't hard code
            gcd = a.arith_int().gcd_int
        elif self.__order <= 2147483647:   # todo: don't hard code
            gcd = a.arith_llong().gcd_longlong
        else:
            raise NotImplementedError("list_of_elements_of_multiplicative_group() is not implemented for large moduli")
        N = self.__order
        # Don't use N.coprime_integers() here because we want Python ints
        return [i for i in range(N) if gcd(i, N) == 1]

    @cached_method
    def multiplicative_subgroups(self):
        r"""
        Return generators for each subgroup of `(\ZZ/N\ZZ)^*`.

        EXAMPLES::

            sage: # optional - gap_package_polycyclic, needs sage.groups
            sage: Integers(5).multiplicative_subgroups()
            ((2,), (4,), ())
            sage: Integers(15).multiplicative_subgroups()
            ((11, 7), (11, 4), (2,), (11,), (14,), (7,), (4,), ())
            sage: Integers(2).multiplicative_subgroups()
            ((),)
            sage: len(Integers(341).multiplicative_subgroups())
            80

        TESTS::

            sage: IntegerModRing(1).multiplicative_subgroups()                          # needs sage.groups
            ((),)
            sage: IntegerModRing(2).multiplicative_subgroups()                          # needs sage.groups
            ((),)
            sage: IntegerModRing(3).multiplicative_subgroups()  # optional - gap_package_polycyclic, needs sage.groups
            ((2,), ())
        """
        return tuple(tuple(g.value() for g in H.gens())
                     for H in self.unit_group().subgroups())

    def is_integral_domain(self, proof=None):
        """
        Return ``True`` if and only if the order of ``self`` is prime.

        EXAMPLES::

            sage: Integers(389).is_integral_domain()
            True
            sage: Integers(389^2).is_integral_domain()                                  # needs sage.libs.pari
            False

        TESTS:

        Check that :issue:`17453` is fixed::

            sage: R = Zmod(5)
            sage: R in IntegralDomains()
            True
        """
        return self.is_field(proof)

    def is_unique_factorization_domain(self, proof=None):
        """
        Return ``True`` if and only if the order of ``self`` is prime.

        EXAMPLES::

            sage: Integers(389).is_unique_factorization_domain()
            True
            sage: Integers(389^2).is_unique_factorization_domain()                      # needs sage.libs.pari
            False
        """
        return self.is_field(proof)

    @cached_method
    def is_field(self, proof=None):
        r"""
        Return ``True`` precisely if the order is prime.

        INPUT:

        - ``proof`` -- boolean or ``None`` (default). If ``False``, then test
          whether the category of the quotient is a subcategory of
          ``Fields()``, or do a probabilistic primality test. If ``None``, then
          test the category and then do a primality test according to the
          global arithmetic proof settings. If ``True``, do a deterministic
          primality test.

        If it is found (perhaps probabilistically) that the ring is a field,
        then the category of the ring is refined to include the category
        of fields. This may change the Python class of the ring!

        EXAMPLES::

            sage: R = IntegerModRing(18)
            sage: R.is_field()
            False
            sage: FF = IntegerModRing(17)
            sage: FF.is_field()
            True

        By :issue:`15229`, the category of the ring is refined,
        if it is found that the ring is in fact a field::

            sage: R = IntegerModRing(127)
            sage: R.category()
            Join of Category of finite commutative rings
                and Category of subquotients of monoids
                and Category of quotients of semigroups
                and Category of finite enumerated sets
            sage: R.is_field()
            True
            sage: R.category()
            Join of Category of finite enumerated fields
                and Category of subquotients of monoids
                and Category of quotients of semigroups

        It is possible to mistakenly put `\ZZ/n\ZZ` into the category of fields.
        In this case, :meth:`is_field` will return ``True`` without performing a
        primality check. However, if the optional argument ``proof=True`` is
        provided, primality is tested and the mistake is uncovered in a warning
        message::

            sage: R = IntegerModRing(21, is_field=True)
            sage: R.is_field()
            True
            sage: R.is_field(proof=True)
            Traceback (most recent call last):
            ...
            ValueError: THIS SAGE SESSION MIGHT BE SERIOUSLY COMPROMISED!
            The order 21 is not prime, but this ring has been put
            into the category of fields. This may already have consequences
            in other parts of Sage. Either it was a mistake of the user,
            or a probabilistic primality test has failed.
            In the latter case, please inform the developers.

        To avoid side-effects of this test on other tests, we clear the cache
        of the ring factory::

            sage: IntegerModRing._cache.clear()
        """
        from sage.categories.fields import Fields
        if not proof:
            if self.category().is_subcategory(Fields()):
                return True
        is_prime = self.order().is_prime(proof=proof)
        if is_prime:
            self._refine_category_(Fields())
            self._factory_data[3]['category'] = Fields()
        else:
            if self.category().is_subcategory(Fields()):
                raise ValueError("""THIS SAGE SESSION MIGHT BE SERIOUSLY COMPROMISED!
The order {} is not prime, but this ring has been put
into the category of fields. This may already have consequences
in other parts of Sage. Either it was a mistake of the user,
or a probabilistic primality test has failed.
In the latter case, please inform the developers.""".format(self.order()))
        return is_prime

    @cached_method
    def field(self):
        """
        If this ring is a field, return the corresponding field as a finite
        field, which may have extra functionality and structure. Otherwise,
        raise a :exc:`ValueError`.

        EXAMPLES::

            sage: R = Integers(7); R
            Ring of integers modulo 7
            sage: R.field()
            Finite Field of size 7
            sage: R = Integers(9)
            sage: R.field()
            Traceback (most recent call last):
            ...
            ValueError: self must be a field
        """
        try:
            return self.__field
        except AttributeError:
            if not self.is_field():
                raise ValueError("self must be a field")
            from . import finite_field_constructor
            k = finite_field_constructor.FiniteField(self.order())
            self.__field = k
            return k

    def _pseudo_fraction_field(self):
        """
        If ``self`` is composite, we may still want to do division by elements
        of ``self``.

        EXAMPLES::

            sage: Integers(15).fraction_field()
            Traceback (most recent call last):
            ...
            TypeError: self must be an integral domain.
            sage: Integers(15)._pseudo_fraction_field()
            Ring of integers modulo 15
            sage: R.<x> = Integers(15)[]
            sage: (x+5)/2
            8*x + 10

        This should be very fast::

            sage: R.<x> = Integers(next_prime(10^101)*next_prime(10^100))[]             # needs sage.libs.pari
            sage: x / R.base_ring()(2)                                                  # needs sage.libs.pari
            500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000013365000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000401*x
        """
        return self

    @cached_method
    def multiplicative_group_is_cyclic(self):
        """
        Return ``True`` if the multiplicative group of this field is cyclic.
        This is the case exactly when the order is less than 8, a power
        of an odd prime, or twice a power of an odd prime.

        EXAMPLES::

            sage: R = Integers(7); R
            Ring of integers modulo 7
            sage: R.multiplicative_group_is_cyclic()
            True
            sage: R = Integers(9)
            sage: R.multiplicative_group_is_cyclic()                                    # needs sage.libs.pari
            True
            sage: Integers(8).multiplicative_group_is_cyclic()
            False
            sage: Integers(4).multiplicative_group_is_cyclic()
            True
            sage: Integers(25*3).multiplicative_group_is_cyclic()                       # needs sage.libs.pari
            False

        We test that :issue:`5250` is fixed::

            sage: Integers(162).multiplicative_group_is_cyclic()                        # needs sage.libs.pari
            True
        """
        n = self.order()
        if n < 8:
            return True

        if n % 4 == 0:
            return False  # know n > 7, so n=4 case not a problem
        if n % 4 == 2:
            n = n // 2

        return n.is_prime_power()

    @cached_method
    def multiplicative_generator(self):
        """
        Return a generator for the multiplicative group of this ring,
        assuming the multiplicative group is cyclic.

        Use the unit_gens function to obtain generators even in the
        non-cyclic case.

        EXAMPLES::

            sage: # needs sage.groups sage.libs.pari
            sage: R = Integers(7); R
            Ring of integers modulo 7
            sage: R.multiplicative_generator()
            3
            sage: R = Integers(9)
            sage: R.multiplicative_generator()
            2
            sage: Integers(8).multiplicative_generator()
            Traceback (most recent call last):
            ...
            ValueError: multiplicative group of this ring is not cyclic
            sage: Integers(4).multiplicative_generator()
            3
            sage: Integers(25*3).multiplicative_generator()
            Traceback (most recent call last):
            ...
            ValueError: multiplicative group of this ring is not cyclic
            sage: Integers(25*3).unit_gens()
            (26, 52)
            sage: Integers(162).unit_gens()
            (83,)
        """
        try:
            return self.__mult_gen
        except AttributeError:
            if self.is_field():
                a = self(self.field().multiplicative_generator())
                self.__mult_gen = a
                return a
            if self.multiplicative_group_is_cyclic():
                v = self.unit_gens()
                if len(v) != 1:
                    raise ArithmeticError
                return v[0]

            raise ValueError("multiplicative group of this ring is not cyclic")

    def quadratic_nonresidue(self):
        """
        Return a quadratic non-residue in ``self``.

        EXAMPLES::

            sage: R = Integers(17)
            sage: R.quadratic_nonresidue()                                              # needs sage.libs.pari
            3
            sage: R(3).is_square()
            False
        """
        try:
            return self._nonresidue
        except AttributeError:
            for a in self:
                if not a.is_square():
                    self._nonresidue = a
                    return a

    def square_roots_of_one(self):
        """
        Return all square roots of 1 in self, i.e., all solutions to
        `x^2 - 1 = 0`.

        OUTPUT: the square roots of 1 in ``self`` as a tuple

        EXAMPLES::

            sage: R = Integers(2^10)
            sage: [x for x in R if x^2 == 1]
            [1, 511, 513, 1023]
            sage: R.square_roots_of_one()
            (1, 511, 513, 1023)

        ::

            sage: # needs sage.libs.pari
            sage: v = Integers(9*5).square_roots_of_one(); v
            (1, 19, 26, 44)
            sage: [x^2 for x in v]
            [1, 1, 1, 1]
            sage: v = Integers(9*5*8).square_roots_of_one(); v
            (1, 19, 71, 89, 91, 109, 161, 179, 181, 199, 251, 269, 271, 289, 341, 359)
            sage: [x^2 for x in v]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        try:
            return self.__square_roots_of_one
        except AttributeError:
            pass
        n = self.__order
        if n.is_prime_power():
            if n % 2 == 0:
                # power of 2
                if n == 2:
                    v = [self(1)]
                elif n == 4:
                    v = [self(1), self(3)]
                else:  # n >= 8
                    half_ord = n // 2
                    v = [self(1), self(-1),
                         self(half_ord - 1), self(half_ord + 1)]
            else:
                v = [self(1), self(-1)]
        else:
            # Reduce to the prime power case.
            F = self.factored_order()
            vmod = []
            moduli = []
            for p, e in F:
                k = p**e
                R = IntegerModRing(p**e)
                w = [self(x) for x in R.square_roots_of_one()]
                vmod.append(w)
                moduli.append(k)
            # Now combine in all possible ways using the CRT
            basis = CRT_basis(moduli)
            from sage.misc.mrange import cartesian_product_iterator
            v = []
            for x in cartesian_product_iterator(vmod):
                # x is a specific choice of roots modulo each prime power divisor
                a = sum([basis[i] * x[i] for i in range(len(x))])
                v.append(a)

        v.sort()
        v = tuple(v)
        self.__square_roots_of_one = v
        return v

    @cached_method
    def factored_order(self):
        """
        EXAMPLES::

            sage: R = IntegerModRing(18)
            sage: FF = IntegerModRing(17)
            sage: R.factored_order()
            2 * 3^2
            sage: FF.factored_order()
            17
        """
        return factor(self.__order, int_=(self.__order < 2**31))

    def factored_unit_order(self):
        r"""
        Return a list of :class:`Factorization` objects, each the factorization
        of the order of the units in a `\ZZ / p^n \ZZ` component of this group
        (using the Chinese Remainder Theorem).

        EXAMPLES::

            sage: R = Integers(8*9*25*17*29)
            sage: R.factored_unit_order()
            [2^2, 2 * 3, 2^2 * 5, 2^4, 2^2 * 7]
        """
        ans = []
        from sage.structure.factorization import Factorization
        for p, e in self.factored_order():
            ans.append(Factorization([(p, e - 1)]) *
                       factor(p - 1, int_=(self.__order < 2**31)))
        return ans

    def characteristic(self):
        """
        EXAMPLES::

            sage: R = IntegerModRing(18)
            sage: FF = IntegerModRing(17)
            sage: FF.characteristic()
            17
            sage: R.characteristic()
            18
        """
        return self.__order

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: Zmod(87)
            Ring of integers modulo 87
        """
        return "Ring of integers modulo {}".format(self.__order)

    def _latex_(self):
        r"""
        Latex representation.

        EXAMPLES::

            sage: latex(Zmod(87))
            \ZZ/87\ZZ
        """
        return "\\ZZ/{}\\ZZ".format(self.__order)

    def modulus(self):
        r"""
        Return the polynomial `x - 1` over this ring.

        .. NOTE::

           This function exists for consistency with the finite-field
           modulus function.

        EXAMPLES::

            sage: R = IntegerModRing(18)
            sage: R.modulus()
            x + 17
            sage: R = IntegerModRing(17)
            sage: R.modulus()
            x + 16
        """
        try:
            return self.__modulus
        except AttributeError:
            x = self['x'].gen()
            self.__modulus = x - 1
            return self.__modulus

    def order(self):
        """
        Return the order of this ring.

        EXAMPLES::

            sage: Zmod(87).order()
            87
        """
        return self.__order

    def cardinality(self):
        """
        Return the cardinality of this ring.

        EXAMPLES::

            sage: Zmod(87).cardinality()
            87
        """
        return self.order()

    def _pari_order(self):
        """
        Return the pari integer representing the order of this ring.

        EXAMPLES::

            sage: Zmod(87)._pari_order()                                                # needs sage.libs.pari
            87
        """
        try:
            return self.__pari_order
        except AttributeError:
            self.__pari_order = pari(self.order())
            return self.__pari_order

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: K2 = GF(2)
            sage: K3 = GF(3)
            sage: K8 = GF(8, 'a')
            sage: K8(5)  # indirect doctest
            1
            sage: K8('a+1')
            a + 1
            sage: K8(K2(1))
            1

        The following test refers to :issue:`6468`::

            sage: class foo_parent(Parent):
            ....:     pass
            sage: class foo(RingElement):
            ....:     def lift(self):
            ....:         raise PariError
            sage: P = foo_parent()
            sage: F = foo(P)
            sage: GF(2)(F)
            Traceback (most recent call last):
            ...
            TypeError: error coercing to finite field

        The following test refers to :issue:`8970`::

            sage: R = Zmod(13); a = R(2)
            sage: a == R(gap(a))                                                        # needs sage.libs.gap
            True

        libgap interface (:issue:`23714`)::

            sage: a = libgap.eval("Z(13)^2")                                            # needs sage.libs.gap
            sage: a.sage()                                                              # needs sage.libs.gap
            4
            sage: libgap(a.sage()) == a                                                 # needs sage.libs.gap
            True

        better syntax for libgap interface::

            sage: a = libgap.Z(13)^2                                                    # needs sage.libs.gap
            sage: libgap(a.sage()) == a                                                 # needs sage.libs.gap
            True
        """
        try:
            return integer_mod.IntegerMod(self, x)
        except (NotImplementedError, PariError):
            raise TypeError("error coercing to finite field")
        except TypeError:
            if isinstance(x, GapElement):
                from sage.libs.gap.libgap import libgap
                return libgap(x).sage()
            raise  # Continue up with the original TypeError

    def __iter__(self):
        """
        EXAMPLES::

            sage: R = IntegerModRing(3)
            sage: for i in R:
            ....:     print(i)
            0
            1
            2
            sage: L = [i for i in R]
            sage: L[0].parent()
            Ring of integers modulo 3
        """
        i = 0
        order = int(self.__order)
        while i < order:
            yield self(i)
            i = i + 1

    def _coerce_map_from_(self, S):
        r"""
        EXAMPLES::

            sage: R = Integers(15)
            sage: f = R.coerce_map_from(Integers(450)); f  # indirect doctest
            Natural morphism:
              From: Ring of integers modulo 450
              To:   Ring of integers modulo 15
            sage: f(-1)
            14
            sage: f = R.coerce_map_from(int); f
            Native morphism:
              From: Set of Python objects of class 'int'
              To:   Ring of integers modulo 15
            sage: f(-1r)
            14
            sage: f = R.coerce_map_from(ZZ); f
            Natural morphism:
              From: Integer Ring
              To:   Ring of integers modulo 15
            sage: f(-1)
            14
            sage: f = R.coerce_map_from(Integers(10)); print(f)
            None
            sage: f = R.coerce_map_from(QQ); print(f)
            None

            sage: R = IntegerModRing(17)
            sage: a = R(3)
            sage: b = R.coerce(3)
            sage: b
            3
            sage: a==b
            True

        This is allowed::

            sage: R(2/3)
            12

        But this is not, since there is no (canonical or not!) ring
        homomorphism from `\QQ` to `\GF{17}`.

        ::

            sage: R.coerce(2/3)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Ring of integers modulo 17

        We do not allow the coercion ``GF(p) -> Z/pZ``, because in case of a
        canonical isomorphism, there is a coercion map in only one
        direction, i.e., to the object in the smaller category.
        """
        if S is int:
            return integer_mod.Int_to_IntegerMod(self)
        elif S is integer_ring.ZZ:
            return integer_mod.Integer_to_IntegerMod(self)
        elif isinstance(S, IntegerModRing_generic):
            if isinstance(S, Field):
                return None
            try:
                return integer_mod.IntegerMod_to_IntegerMod(S, self)
            except TypeError:
                pass
        to_ZZ = integer_ring.ZZ._internal_coerce_map_from(S)
        if to_ZZ is not None:
            return integer_mod.Integer_to_IntegerMod(self) * to_ZZ

    def _convert_map_from_(self, other):
        """
        Conversion from `p`-adic fields.

        EXAMPLES::

            sage: Zmod(81).convert_map_from(Qp(3))                                      # needs sage.rings.padics
            Reduction morphism:
              From: 3-adic Field with capped relative precision 20
              To:   Ring of integers modulo 81
        """
        from sage.rings.padics.padic_generic import pAdicGeneric, ResidueReductionMap
        if isinstance(other, pAdicGeneric) and other.degree() == 1:
            p = other.prime()
            N = self.cardinality()
            n = N.exact_log(p)
            if p**n == N:
                return ResidueReductionMap._create_(other, self)

    def __richcmp__(self, other, op):
        """
        EXAMPLES::

            sage: Z11 = IntegerModRing(11); Z11
            Ring of integers modulo 11
            sage: Z12 = IntegerModRing(12); Z12
            Ring of integers modulo 12
            sage: Z13 = IntegerModRing(13); Z13
            Ring of integers modulo 13
            sage: Z11 == Z11, Z11 == Z12, Z11 == Z13
            (True, False, False)
            sage: F = GF(11); F
            Finite Field of size 11
            sage: Z11 == F
            False

        In :issue:`15229`, the following was implemented::

            sage: R1 = IntegerModRing(5)
            sage: R2 = IntegerModRing(5, is_field=True)
            sage: R1 is R2    # used to return False
            True
            sage: R2 == GF(5)
            False
        """
        # We want that GF(p) and IntegerModRing(p) evaluate unequal.
        # However, we cannot just compare the types, since the
        # choice of a different category also changes the type.
        # But if we go to the base class, we avoid the influence
        # of the category.
        try:
            c = bool(other.__class__.__base__ != self.__class__.__base__)
        except AttributeError:  # __base__ does not always exists
            c = bool(type(other) is not type(self))
        if c:
            return NotImplemented
        return richcmp(self.__order, other.__order, op)

    def unit_gens(self, **kwds):
        r"""
        Return generators for the unit group `(\ZZ/N\ZZ)^*`.

        We compute the list of generators using a deterministic algorithm, so
        the generators list will always be the same. For each odd prime divisor
        of `N` there will be exactly one corresponding generator; if `N` is
        even there will be 0, 1 or 2 generators according to whether 2 divides
        `N` to order 1, 2 or `\geq 3`.

        OUTPUT: a tuple containing the units of ``self``

        EXAMPLES::

            sage: R = IntegerModRing(18)
            sage: R.unit_gens()                                                         # needs sage.groups
            (11,)
            sage: R = IntegerModRing(17)
            sage: R.unit_gens()                                                         # needs sage.groups
            (3,)
            sage: IntegerModRing(next_prime(10^30)).unit_gens()                         # needs sage.groups
            (5,)

        The choice of generators is affected by the optional keyword
        ``algorithm``; this can be ``'sage'`` (default) or ``'pari'``.
        See :meth:`unit_group` for details. ::

            sage: A = Zmod(55)
            sage: A.unit_gens(algorithm='sage')                                         # needs sage.groups
            (12, 46)
            sage: A.unit_gens(algorithm='pari')                                         # needs sage.groups sage.libs.pari
            (2, 21)

        TESTS::

            sage: IntegerModRing(2).unit_gens()                                         # needs sage.groups
            ()
            sage: IntegerModRing(4).unit_gens()                                         # needs sage.groups
            (3,)
            sage: IntegerModRing(8).unit_gens()                                         # needs sage.groups
            (7, 5)
        """
        return self.unit_group(**kwds).gens_values()

    def unit_group_exponent(self):
        """
        EXAMPLES::

            sage: R = IntegerModRing(17)
            sage: R.unit_group_exponent()                                               # needs sage.groups
            16
            sage: R = IntegerModRing(18)
            sage: R.unit_group_exponent()                                               # needs sage.groups
            6
        """
        return self.unit_group().exponent()

    def unit_group_order(self):
        """
        Return the order of the unit group of this residue class ring.

        EXAMPLES::

            sage: R = Integers(500)
            sage: R.unit_group_order()                                                  # needs sage.groups
            200
        """
        return self.unit_group().order()

    @cached_method
    def unit_group(self, algorithm='sage'):
        r"""
        Return the unit group of ``self``.

        INPUT:

        - ``self`` -- the ring `\ZZ/n\ZZ` for a positive integer `n`

        - ``algorithm`` -- either ``'sage'`` (default) or ``'pari'``

        OUTPUT:

        The unit group of ``self``.  This is a finite Abelian group
        equipped with a distinguished set of generators, which is
        computed using a deterministic algorithm depending on the
        ``algorithm`` parameter.

        - If ``algorithm == 'sage'``, the generators correspond to the
          prime factors `p \mid n` (one generator for each odd `p`;
          the number of generators for `p = 2` is 0, 1 or 2 depending
          on the order to which 2 divides `n`).

        - If ``algorithm == 'pari'``, the generators are chosen such
          that their orders form a decreasing sequence with respect to
          divisibility.

        EXAMPLES:

        The output of the algorithms ``'sage'`` and ``'pari'`` can
        differ in various ways.  In the following example, the same
        cyclic factors are computed, but in a different order::

            sage: # needs sage.groups
            sage: A = Zmod(15)
            sage: G = A.unit_group(); G
            Multiplicative Abelian group isomorphic to C2 x C4
            sage: G.gens_values()
            (11, 7)
            sage: H = A.unit_group(algorithm='pari'); H                                 # needs sage.libs.pari
            Multiplicative Abelian group isomorphic to C4 x C2
            sage: H.gens_values()                                                       # needs sage.libs.pari
            (7, 11)

        Here are two examples where the cyclic factors are isomorphic,
        but are ordered differently and have different generators::

            sage: # needs sage.groups
            sage: A = Zmod(40)
            sage: G = A.unit_group(); G
            Multiplicative Abelian group isomorphic to C2 x C2 x C4
            sage: G.gens_values()
            (31, 21, 17)
            sage: H = A.unit_group(algorithm='pari'); H                                 # needs sage.libs.pari
            Multiplicative Abelian group isomorphic to C4 x C2 x C2
            sage: H.gens_values()                                                       # needs sage.libs.pari
            (17, 31, 21)

            sage: # needs sage.groups
            sage: A = Zmod(192)
            sage: G = A.unit_group(); G
            Multiplicative Abelian group isomorphic to C2 x C16 x C2
            sage: G.gens_values()
            (127, 133, 65)
            sage: H = A.unit_group(algorithm='pari'); H                                 # needs sage.libs.pari
            Multiplicative Abelian group isomorphic to C16 x C2 x C2
            sage: H.gens_values()                                                       # needs sage.libs.pari
            (133, 127, 65)

        In the following examples, the cyclic factors are not even
        isomorphic::

            sage: A = Zmod(319)
            sage: A.unit_group()                                                        # needs sage.groups
            Multiplicative Abelian group isomorphic to C10 x C28
            sage: A.unit_group(algorithm='pari')                                        # needs sage.groups sage.libs.pari
            Multiplicative Abelian group isomorphic to C140 x C2

            sage: A = Zmod(30.factorial())
            sage: A.unit_group()                                                        # needs sage.groups
            Multiplicative Abelian group isomorphic to
             C2 x C16777216 x C3188646 x C62500 x C2058 x C110 x C156 x C16 x C18 x C22 x C28
            sage: A.unit_group(algorithm='pari')                                        # needs sage.groups sage.libs.pari
            Multiplicative Abelian group isomorphic to
             C20499647385305088000000 x C55440 x C12 x C12 x C4 x C2 x C2 x C2 x C2 x C2 x C2

        TESTS:

        We test the cases where the unit group is trivial::

            sage: # needs sage.groups
            sage: A = Zmod(1)
            sage: A.unit_group()
            Trivial Abelian group
            sage: A.unit_group(algorithm='pari')                                        # needs sage.libs.pari
            Trivial Abelian group
            sage: A = Zmod(2)
            sage: A.unit_group()
            Trivial Abelian group
            sage: A.unit_group(algorithm='pari')                                        # needs sage.libs.pari
            Trivial Abelian group

            sage: Zmod(3).unit_group(algorithm='bogus')                                 # needs sage.groups
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm 'bogus' for computing the unit group
        """
        from sage.groups.abelian_gps.values import AbelianGroupWithValues
        if algorithm == 'sage':
            n = self.order()
            gens = []
            orders = []
            for p, r in self.factored_order():
                m = n // (p**r)
                for g, o in _unit_gens_primepowercase(p, r):
                    x = g.crt(integer_mod.Mod(1, m))
                    gens.append(x)
                    orders.append(o)
        elif algorithm == 'pari':
            _, orders, gens = self.order().__pari__().znstar()
            gens = [self(g) for g in gens]
            orders = [integer.Integer(o) for o in orders]
        else:
            raise ValueError('unknown algorithm %r for computing the unit group' % algorithm)
        return AbelianGroupWithValues(gens, orders, values_group=self)

    def random_element(self, bound=None):
        """
        Return a random element of this ring.

        INPUT:

        - ``bound`` -- positive integer or ``None`` (the default); if given,
          return  the coercion of an integer in the interval
          ``[-bound, bound]`` into this ring

        EXAMPLES::

            sage: R = IntegerModRing(18)
            sage: R.random_element().parent() is R
            True
            sage: found = [False]*18
            sage: while not all(found):
            ....:     found[R.random_element()] = True

        We test ``bound``-option::

            sage: R.random_element(2) in [R(16), R(17), R(0), R(1), R(2)]
            True
        """
        if bound is not None:
            return CommutativeRing.random_element(self, bound)
        a = random.randint(0, self.order() - 1)
        return self(a)

    @staticmethod
    def _lift_residue_field_root(p, e, f, fprime, root):
        """Lifts a root of f

        INPUT:

        - ``p`` -- integer, a prime number
        - ``e`` -- positive integer
        - ``f`` -- polynomial with coefficients in ``IntegerModRing(p**e)``
        - ``fprime`` -- derivative of ``f``
        - ``root`` -- Element of ``IntegerModRing(p)`` with ``f(root) = 0``

        OUTPUT: iterable of roots of ``f`` modulo ``p**e``.  Each root is an
        ``IntegerModRing(p**e)`` element.

        TESTS::

            sage: R = Zmod(2)
            sage: S.<x> = R[]
            sage: R._lift_residue_field_root(2, 1, S.zero(), S.zero(), R(0))
            (0,)

        Lifting roots of the zero polynomial::

            sage: R = Zmod(41)
            sage: S.<x> = R[]
            sage: R._lift_residue_field_root(41, 1, S.zero(), S.zero(), R(12))
            (12,)
            sage: R = Zmod(5**2)
            sage: S.<x> = R[]
            sage: R._lift_residue_field_root(5, 2, S.zero(), S.zero(), R(2))
            [2, 7, 12, 17, 22]
            sage: R = Zmod(2**3)
            sage: S.<x> = R[]
            sage: R._lift_residue_field_root(2, 3, S.zero(), S.zero(), R(1))
            [1, 5, 3, 7]

        Trivial case where ``e == 1``::

            sage: R = Zmod(41)
            sage: S.<x> = R[]
            sage: f = x^2 - 2
            sage: R._lift_residue_field_root(41, 1, f, f.derivative(), R(17))
            (17,)

            sage: R = Zmod(43)
            sage: S.<x> = R[]
            sage: f = x^43 - 3
            sage: R._lift_residue_field_root(43, 1, f, f.derivative(), R(3))
            (3,)

        Non-singular cases with one step of lifting::

            sage: R = Zmod(2**2)
            sage: S.<x> = R[]
            sage: f = x - 1
            sage: R._lift_residue_field_root(2, 2, f, f.derivative(), R(1))
            (1,)
            sage: f = x - 3
            sage: R._lift_residue_field_root(2, 2, f, f.derivative(), R(1))
            (3,)

            sage: R = Zmod(4001**2)
            sage: S.<x> = R[]
            sage: f = x^3 - 2
            sage: R._lift_residue_field_root(4001, 2, f, f.derivative(), R(3981))
            (5309307,)
            sage: f = x^3 - 3
            sage: R._lift_residue_field_root(4001, 2, f, f.derivative(), R(1091))
            (11035849,)

        Non-singular cases with multiple steps of lifting::

            sage: R = Zmod(2**10)
            sage: S.<x> = R[]
            sage: f = x + 1
            sage: R._lift_residue_field_root(2, 10, f, f.derivative(), Zmod(2)(1))
            (1023,)

            sage: R = Zmod(2**16)
            sage: S.<x> = R[]
            sage: f = x + 1
            sage: R._lift_residue_field_root(2, 16, f, f.derivative(), Zmod(2)(1))
            (65535,)

            sage: R = Zmod(7**4)
            sage: S.<x> = R[]
            sage: f = x^4 - 2
            sage: R._lift_residue_field_root(7, 4, f, f.derivative(), Zmod(7)(2))
            (121,)

        Singular cases::

            sage: R = Zmod(2**3)
            sage: S.<x> = R[]
            sage: f = x^2 - 1
            sage: R._lift_residue_field_root(2, 3, f, f.derivative(), Zmod(2)(1))
            [1, 5, 3, 7]
            sage: f = 2*x
            sage: R._lift_residue_field_root(2, 3, f, f.derivative(), Zmod(2)(0))
            [0, 4]

            sage: R = Zmod(11**2)
            sage: S.<x> = R[]
            sage: f = x^2 + 13*x + 1
            sage: R._lift_residue_field_root(11, 2, f, f.derivative(), Zmod(11)(10))
            []

            sage: R = Zmod(11**3)
            sage: S.<x> = R[]
            sage: f = x^2 + 123*x + 1
            sage: R._lift_residue_field_root(11, 3, f, f.derivative(), Zmod(11)(10))
            [10, 131, 252, 373, 494, 615, 736, 857, 978, 1099, 1220, 109, 230, 351, 472, 593, 714, 835, 956, 1077, 1198, 1319]
        """
        if e == 1:
            # Nothing to do
            return (root,)
        deriv = fprime(root)
        if deriv:
            # Unique lift, use Newton iteration
            prec = 1
            while True:
                prec = min(2*prec, e)
                Zp_prec = Zmod(p**prec)
                root = Zp_prec(root.lift())
                deriv = fprime(root)
                step = f(root) / deriv
                root -= step
                if prec >= e:
                    return (root,)
        else:
            # Non-unique lift, go one power at a time
            prec = 1
            new_power = 1
            new_mod = p
            current_roots = (root,)
            for _ in range(e - 1):
                prec += 1
                new_power = new_mod
                new_mod *= p
                new_roots = []
                Zp_prec = Zmod(new_mod)
                for rt in current_roots:
                    rt = Zp_prec(rt.lift())
                    if f(rt):
                        continue
                    new_roots.append(rt)
                    for _ in range(p - 1):
                        rt += new_power
                        new_roots.append(rt)
                current_roots = new_roots

            return current_roots

    def _roots_univariate_polynomial(self, f, ring=None, multiplicities=True, algorithm=None):
        r"""
        Return the roots of ``f`` in the ring ``ring``.

        INPUT:

        - ``f`` - a polynomial defined over this ring

        - ``ring`` - the ring to find roots in.  Control flow elsewhere
          ensures that the only cases we need to handle are ``self`` and
          ``None``.  Otherwise we raise ``NotImplementedError``.

        - ``multiplicities`` - bool (default: ``True``). If ``True``, return
          list of pairs `(r, n)`, where `r` is a root and `n` is its
          multiplicity. If ``False``, just return the unique roots, with no
          information about multiplicities. Multiplicities are only defined
          over fields, and this method raises ``NotImplementedError`` if this
          is ``True`` but the ring is not a field.

        - ``algorithm`` - ignored

        ALGORITHM:

        The algorithm is adapted from [Gou2020]_, section 4.5, and [Coh1993]_,
        section 3.5.3.  It is a combination of the Chinese Remainder Theorem
        and Hensel's lemma.  As a base case, if `N` is prime, then we find
        roots by factoring `f`.  If `N` is a prime power `p^e`, then we find
        roots modulo `p` and lift them.  Finally, for general `N`, we first
        factor the modulus `N` into prime powers, list all roots modulo those
        prime powers, and combine the roots using the Chinese Remainder
        Theorem.

        Suppose that we are trying to find roots modulo `p^e` and that `r` is
        a root of `f(x)` modulo `p`.  The easy case is when `f'(r) \not\equiv
        0 \pmod{p}`, for then Hensel's lemma implies that there is a unique
        `r_e \in \Zmod{p^e}` with `r_e \equiv r \pmod{p}`.
        Moreover, this `r_e` can be found by applying Newton's method for
        numerically approximating roots.  Each iteration of Newton's method
        doubles the precision to which the root is known.

        But if `f'(r) \equiv 0 \pmod{p}`, then this is no longer true.  In
        fact, in this case roots modulo `p^e` are not the same as `p`-adic
        roots, and finding all the latter does not guarantee that we have
        found all the former.  For example, if `f(x) = 2x` and `p = 2`, then
        there is only one `p`-adic root, namely zero.  But the solutions of
        `2x \equiv 0 \pmod{2^k}` are `0` and `2^{k-1}`; the former lifts to
        two roots modulo `2^{k+1}`, namely `0` and `2^k`, while the latter
        does not lift at all.  We handle this case by lifting one power at a
        time.  While we can no longer use Newton's method to solve for a lift,
        the Taylor series it is based on still yields constraints on the roots
        modulo `p^{k+1}`:  If `r_k` is a root of `f` modulo `p^k`, then either
        every lift of `r_k` to `\Zmod{p^{k + 1}}` is a root of `f`
        modulo `p^{k+1}` or none of them are.  Consequently we may find roots
        modulo `p^e` by lifting one power at a time.

        When `f'(r) \equiv 0 \pmod{p}`, an alternative approach is to change
        variables, factor out the root, and then factor out powers of `p`.
        This has the advantage that it will eventually reach a situation where
        the lift converges quadratically, but it is not presently implemented.
        A different form of Hensel's lemma applies once we are close enough to
        a `p`-adic root (see [Gou2020]_, problem 120), but it seems delicate
        to use it directly to find all roots modulo `p^e` (consider our
        earlier example of `f(x) = 2x`), so we do not presently attempt to
        apply Hensel's lemma in this way.

        EXAMPLES::

            sage: R.<x> = Zmod(41)[]
            sage: (x^3 + x).roots()
            [(0, 1), (32, 1), (9, 1)]
            sage: (x^3 + x).roots(multiplicities=False)
            [0, 32, 9]
            sage: (x^6 + x^5 + 9*x^4 + 20*x^3 + 3*x^2 + 18*x + 7).roots()
            [(19, 1), (20, 2), (21, 3)]
            sage: (x^6 + x^5 + 9*x^4 + 20*x^3 + 3*x^2 + 18*x + 7).roots(multiplicities=False)
            [19, 20, 21]

        We can find roots without multiplicities over a ring whose modulus is
        a prime power, even a big power:

            sage: R.<x> = Zmod(7^3)[]
            sage: (x^2 + x + 1).roots(multiplicities=False)
            [18, 324]
            sage: R.<x> = Zmod(2^50)[]
            sage: (x + 1).roots(multiplicities=False)
            [1125899906842623]

        We can also find roots without multiplicities over a ring whose modulus
        is a product of primes or prime powers:

            sage: R.<x> = Zmod(60)[]
            sage: (x^2 - 1).roots(multiplicities=False)
            [29, 41, 49, 1, 59, 11, 19, 31]

        We may also ask for roots modulo a quotient of the ring over which the
        polynomial is defined:

            sage: R.<x> = Zmod(120)[]
            sage: (x^2 - 1).roots(multiplicities=False)
            [89, 41, 49, 1, 29, 101, 109, 61, 59, 11, 19, 91, 119, 71, 79, 31]
            sage: (x^2 - 1).roots(Zmod(60), multiplicities=False)
            [29, 41, 49, 1, 59, 11, 19, 31]

        TESTS::

            sage: R.<x> = Zmod(2)[]
            sage: x.roots()
            [(0, 1)]

        Test polynomials with content:

            sage: R.<x> = Zmod(4)[]
            sage: (2*x).roots(multiplicities=False)
            [0, 2]

            sage: R.<x> = Zmod(6)[]
            sage: (3*x).roots(multiplicities=False)
            [0, 4, 2]

        Test polynomial with many roots:

            sage: R.<x> = Zmod(6)[]
            sage: f = x * (x - 1) * (x - 2) * (x - 3) * (x - 4) * (x - 5)
            sage: len(f.roots(multiplicities=False))
            6

        Test finding roots over large prime powers:

            sage: R.<x> = Zmod(2**16)[]
            sage: (x^3 + 5).roots(multiplicities=False)
            [45475]
            sage: (x^2 + 46*x + 1).roots(multiplicities=False)
            [421, 33189, 16805, 49573, 8613, 41381, 24997, 57765, 7725, 40493, 24109, 56877, 15917, 48685, 32301, 65069]

            sage: R.<x> = Zmod(3**16)[]
            sage: (x^2 + 2).roots(multiplicities=False)
            [24620738, 18425983]
            sage: (x^2 + 11*x + 1).roots(multiplicities=False)
            [633836, 14982743, 29331650, 13715060, 28063967, 42412874]
            sage: (x^3 + 8).roots(multiplicities=False)
            [14348905, 28697812, 43046719]

        Test some larger primes:

            sage: R.<x> = Zmod(41**4)[]
            sage: (x^2 + 2).roots(multiplicities=False)
            [2208905, 616856]
            sage: R.<x> = Zmod(43**4)[]
            sage: (x^2 + 3).roots(multiplicities=False)
            [3269879, 148922]

        We can't find roots with multiplicities in non-fields:

            sage: R.<x> = Zmod(6)[]
            sage: (x + 1).roots()
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding with multiplicities for this polynomial not implemented (try the multiplicities=False option)
            sage: R.<x> = Zmod(8)[]
            sage: (x + 1).roots()
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding with multiplicities for this polynomial not implemented (try the multiplicities=False option)
            sage: R.<x> = Zmod(12)[]
            sage: (x + 1).roots()
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding with multiplicities for this polynomial not implemented (try the multiplicities=False option)

        The zero polynomial has every residue class as a root, but we don't
        support multiplicities even over fields (they would all be infinite).

            sage: R.<x> = Zmod(6)[]
            sage: R.zero().roots()
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding with multiplicities for this polynomial not implemented (try the multiplicities=False option)
            sage: R.zero().roots(multiplicities=False)
            [0, 1, 2, 3, 4, 5]

            sage: R.<x> = Zmod(7)[]
            sage: R.zero().roots()
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding with multiplicities for this polynomial not implemented (try the multiplicities=False option)
            sage: R.zero().roots(multiplicities=False)
            [0, 1, 2, 3, 4, 5, 6]

            sage: R.<x> = Zmod(8)[]
            sage: R.zero().roots()
            Traceback (most recent call last):
            ...
            NotImplementedError: root finding with multiplicities for this polynomial not implemented (try the multiplicities=False option)
            sage: R.zero().roots(multiplicities=False)
            [0, 1, 2, 3, 4, 5, 6, 7]

        This method doesn't support root-finding over rings that aren't Z/nZ:

            sage: R.<x> = Zmod(120)[]
            sage: f = x^2 - 1
            sage: f.base_ring()._roots_univariate_polynomial(f, ring=RR, multiplicities=False)
            Traceback (most recent call last):
            ...
            NotImplementedError

        Sage allows us to coerce polynomials from one modulus to another,
        and that makes the following defined:

            sage: R.<x> = Zmod(100)[]
            sage: (x^2 - 1).roots(Zmod(99), multiplicities=False) == (x^2 - 1).change_ring(Zmod(99)).roots(multiplicities=False)
            True
        """

        # This function only supports roots in an IntegerModRing
        if ring is not self and ring is not None:
            raise NotImplementedError

        deg = f.degree()

        if multiplicities:
            if deg < 0 or not self.is_field():
                raise NotImplementedError(
                    "root finding with multiplicities for this polynomial not"
                    " implemented (try the multiplicities=False option)"
                )
            # Roots of non-zero polynomial over finite fields by factorization
            return f._roots_from_factorization(f.factor(), multiplicities)

        # Zero polynomial is a base case
        if deg < 0:
            # All residue classes are roots of the zero polynomial
            return [*map(self, range(self.cardinality()))]

        # Finite fields are a base case
        if self.is_field():
            return f._roots_from_factorization(f.factor(), False)

        # Otherwise, find roots modulo each prime power
        fac = self.factored_order()
        prime_power_roots = []
        for p, e in fac:
            Zpe = Zmod(p**e)
            fpe = f.change_ring(Zpe)
            fpe_prime = fpe.derivative()
            fp = fpe.change_ring(Zmod(p))

            mod_p_roots = fp.roots(multiplicities=False)

            this_prime_power = []
            for root in mod_p_roots:
                this_prime_power.extend(
                    self._lift_residue_field_root(p, e, fpe, fpe_prime, root)
                )
            prime_power_roots.append(this_prime_power)

        # Combine using Chinese Remainder Theorem
        ppwr_basis = CRT_basis([p**e for p, e in fac])
        result = []
        for res in cartesian_product_iterator(prime_power_roots):
            root = self.zero()
            for c, x in zip(ppwr_basis, res):
                root += c*x.lift()
            result.append(root)
        return result

    #######################################################
    # Suppose for interfaces
    #######################################################
    def _gap_init_(self):
        """
        EXAMPLES::

            sage: R = Integers(12345678900)
            sage: R
            Ring of integers modulo 12345678900
            sage: gap(R)  # indirect doctest                                            # needs sage.libs.gap
            (Integers mod 12345678900)
        """
        return 'ZmodnZ({})'.format(self.order())

    def _magma_init_(self, magma):
        """
        EXAMPLES::

            sage: R = Integers(12345678900)
            sage: R
            Ring of integers modulo 12345678900
            sage: magma(R)  # indirect doctest, optional - magma
            Residue class ring of integers modulo 12345678900
        """
        return 'Integers({})'.format(self.order())

    def degree(self):
        """
        Return 1.

        EXAMPLES::

            sage: R = Integers(12345678900)
            sage: R.degree()
            1
        """
        return integer.Integer(1)


# Register unpickling methods for backward compatibility.

from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.rings.integer_mod_ring', 'IntegerModRing_generic', IntegerModRing_generic)


def crt(v):
    """
    INPUT:

    - ``v`` -- (list) a lift of elements of ``rings.IntegerMod(n)``, for
      various coprime moduli ``n``

    EXAMPLES::

        sage: from sage.rings.finite_rings.integer_mod_ring import crt
        sage: crt([mod(3, 8), mod(1,19), mod(7, 15)])
        1027
    """
    if len(v) == 0:
        return IntegerModRing(1).one()
    x = v[0]
    for i in range(1, len(v)):
        x = x.crt(v[i])
    return x
