r"""
Class groups of number fields

An element of a class group is stored as a pair consisting of both an explicit
ideal in that ideal class, and a list of exponents giving that ideal class in
terms of the generators of the parent class group. These can be accessed with
the ``ideal()`` and ``exponents()`` methods respectively.

EXAMPLES::

    sage: x = polygen(ZZ, 'x')
    sage: K.<a> = NumberField(x^2 + 23)
    sage: I = K.class_group().gen(); I
    Fractional ideal class (2, 1/2*a - 1/2)
    sage: I.ideal()
    Fractional ideal (2, 1/2*a - 1/2)
    sage: I.exponents()
    (1,)

    sage: I.ideal() * I.ideal()
    Fractional ideal (4, 1/2*a + 3/2)
    sage: (I.ideal() * I.ideal()).reduce_equiv()
    Fractional ideal (2, 1/2*a + 1/2)
    sage: J = I * I; J    # class group multiplication is automatically reduced
    Fractional ideal class (2, 1/2*a + 1/2)
    sage: J.ideal()
    Fractional ideal (2, 1/2*a + 1/2)
    sage: J.exponents()
    (2,)

    sage: I * I.ideal()   # ideal classes coerce to their representative ideal
    Fractional ideal (4, 1/2*a + 3/2)

    sage: K.fractional_ideal([2, 1/2*a + 1/2])
    Fractional ideal (2, 1/2*a + 1/2)
    sage: K.fractional_ideal([2, 1/2*a + 1/2]).is_principal()
    False
    sage: K.fractional_ideal([2, 1/2*a + 1/2])^3
    Fractional ideal (1/2*a - 3/2)
"""

from sage.groups.abelian_gps.values import AbelianGroupWithValues_class, AbelianGroupWithValuesElement
from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement
from sage.structure.element import MonoidElement
from sage.rings.integer_ring import ZZ


class FractionalIdealClass(AbelianGroupWithValuesElement):
    r"""
    A fractional ideal class in a number field.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: G = NumberField(x^2 + 23,'a').class_group(); G
        Class group of order 3 with structure C3 of
         Number Field in a with defining polynomial x^2 + 23
        sage: I = G.0; I
        Fractional ideal class (2, 1/2*a - 1/2)
        sage: I.ideal()
        Fractional ideal (2, 1/2*a - 1/2)

        sage: K.<w> = QuadraticField(-23)
        sage: OK = K.ring_of_integers()
        sage: C = OK.class_group()
        sage: P2a, P2b = [P for P, e in (2*K).factor()]
        sage: c = C(P2a); c
        Fractional ideal class (2, 1/2*w - 1/2)
        sage: c.gens()
        (2, 1/2*w - 1/2)
    """
    def __init__(self, parent, element, ideal=None):
        """
        Return the ideal class of this fractional ideal.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 23,'a'); G = K.class_group()
            sage: G(K.ideal(13, a + 4))
            Fractional ideal class (13, 1/2*a + 17/2)
        """
        if element is None:
            element = parent._ideal_log(ideal)
        AbelianGroupWithValuesElement.__init__(self, parent, element, ideal)

    def _repr_(self):
        r"""
        Return string representation of this fractional ideal class.

         EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 23,'a'); G = K.class_group()
            sage: G(K.ideal(13, a + 4))._repr_()
            'Fractional ideal class (13, 1/2*a + 17/2)'
            sage: G(K.ideal(59, a+6))._repr_()
            'Trivial principal fractional ideal class'
        """
        if self.is_principal():
            return 'Trivial principal fractional ideal class'
        return 'Fractional ideal class %s' % self._value._repr_short()

    def _mul_(self, other):
        r"""
        Multiplication of two (S-)ideal classes.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: G = NumberField(x^2 + 23,'a').class_group(); G
            Class group of order 3 with structure C3 of
             Number Field in a with defining polynomial x^2 + 23
            sage: I = G.0; I
            Fractional ideal class (2, 1/2*a - 1/2)
            sage: I*I # indirect doctest
            Fractional ideal class (2, 1/2*a + 1/2)
            sage: I*I*I # indirect doctest
            Trivial principal fractional ideal class

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: G = K.ideal(3,a+1)
            sage: CS(G)*CS(G)
            Trivial S-ideal class
        """
        m = AbelianGroupElement._mul_(self, other)
        m._value = (self.ideal() * other.ideal()).reduce_equiv()
        return m

    def _div_(self, other):
        r"""
        Division of two ideal classes.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: G = NumberField(x^2 + 23,'a').class_group(); G
            Class group of order 3 with structure C3 of
             Number Field in a with defining polynomial x^2 + 23
            sage: I = G.0; I
            Fractional ideal class (2, 1/2*a - 1/2)
            sage: I*I # indirect doctest
            Fractional ideal class (2, 1/2*a + 1/2)
            sage: I*I*I # indirect doctest
            Trivial principal fractional ideal class
        """
        m = AbelianGroupElement._div_(self, other)
        m._value = (self.ideal() / other.ideal()).reduce_equiv()
        return m

    def __pow__(self, n):
        r"""
        Raise this element to the power n.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 3*x + 8)
            sage: C = K.class_group()
            sage: c = C(2, a)
            sage: c^2
            Fractional ideal class (4, a)
            sage: c^3
            Trivial principal fractional ideal class
            sage: c^1000
            Fractional ideal class (2, a)
            sage: (c^2)^2
            Fractional ideal class (2, a)
        """
        # We use MonoidElement's __pow__ routine, since that does
        # repeated squaring, and hence the ideal gets reduced as
        # we go along; actually computing self._value ** n would
        # be disastrous.
        n = n % self.order()
        return MonoidElement.__pow__(self, n)

    def inverse(self):
        r"""
        Return the multiplicative inverse of this ideal class.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 - 3*x + 8); G = K.class_group()
            sage: G(2, a).inverse()
            Fractional ideal class (2, a^2 + 2*a - 1)
            sage: ~G(2, a)
            Fractional ideal class (2, a^2 + 2*a - 1)
        """
        m = AbelianGroupElement.__invert__(self)
        m._value = (~self.ideal()).reduce_equiv()
        return m

    __invert__ = inverse

    def is_principal(self):
        r"""
        Return ``True`` iff this ideal class is the trivial (principal) class.

        EXAMPLES::

            sage: K.<w> = QuadraticField(-23)
            sage: OK = K.ring_of_integers()
            sage: C = OK.class_group()
            sage: P2a, P2b = [P for P, e in (2*K).factor()]
            sage: c = C(P2a)
            sage: c.is_principal()
            False
            sage: (c^2).is_principal()
            False
            sage: (c^3).is_principal()
            True
        """
        return self.is_one()

    def reduce(self):
        r"""
        Return representative for this ideal class that has been
        reduced using PARI's :pari:`idealred`.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: k.<a> = NumberField(x^2 + 20072); G = k.class_group(); G
            Class group of order 76 with structure C38 x C2 of
             Number Field in a with defining polynomial x^2 + 20072
            sage: I = (G.0)^11; I
            Fractional ideal class (33, 1/2*a + 8)
            sage: J = G(I.ideal()^5); J
            Fractional ideal class (39135393, 1/2*a + 13654253)
            sage: J.reduce()
            Fractional ideal class (73, 1/2*a + 47)
            sage: J == I^5
            True
        """
        return self.parent()(self.ideal().reduce_equiv())

    def ideal(self):
        r"""
        Return a representative ideal in this ideal class.

        EXAMPLES::

            sage: K.<w> = QuadraticField(-23)
            sage: OK = K.ring_of_integers()
            sage: C = OK.class_group()
            sage: P2a, P2b = [P for P, e in (2*K).factor()]
            sage: c = C(P2a); c
            Fractional ideal class (2, 1/2*w - 1/2)
            sage: c.ideal()
            Fractional ideal (2, 1/2*w - 1/2)
        """
        return self.value()

    def representative_prime(self, norm_bound=1000):
        r"""
        Return a prime ideal in this ideal class.

        INPUT:

        - ``norm_bound`` -- (positive integer) upper bound on the norm of
          primes tested

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 31)
            sage: K.class_number()
            3
            sage: Cl = K.class_group()
            sage: [c.representative_prime() for c in Cl]
            [Fractional ideal (3),
             Fractional ideal (2, 1/2*a + 1/2),
             Fractional ideal (2, 1/2*a - 1/2)]

            sage: K.<a> = NumberField(x^2 + 223)
            sage: K.class_number()
            7
            sage: Cl = K.class_group()
            sage: [c.representative_prime() for c in Cl]
            [Fractional ideal (3),
             Fractional ideal (2, 1/2*a + 1/2),
             Fractional ideal (17, 1/2*a + 7/2),
             Fractional ideal (7, 1/2*a - 1/2),
             Fractional ideal (7, 1/2*a + 1/2),
             Fractional ideal (17, 1/2*a + 27/2),
             Fractional ideal (2, 1/2*a - 1/2)]
        """
        if self.value().is_prime():
            return self.value()
        c = self.reduce()
        if c.value().is_prime():
            return c.value()
        # otherwise we just search:
        Cl = self.parent()
        K = Cl.number_field()
        from sage.rings.real_mpfr import RR
        for P in K.primes_of_bounded_norm_iter(RR(norm_bound)):
            if Cl(P) == c:
                return P
        raise RuntimeError("No prime of norm less than %s found in class %s" % (norm_bound, c))

    def gens(self):
        r"""
        Return generators for a representative ideal in this
        (`S`-)ideal class.

        EXAMPLES::

            sage: K.<w> = QuadraticField(-23)
            sage: OK = K.ring_of_integers()
            sage: C = OK.class_group()
            sage: P2a, P2b = [P for P, e in (2*K).factor()]
            sage: c = C(P2a); c
            Fractional ideal class (2, 1/2*w - 1/2)
            sage: c.gens()
            (2, 1/2*w - 1/2)
        """
        return self.ideal().gens()


class SFractionalIdealClass(FractionalIdealClass):
    r"""
    An `S`-fractional ideal class in a number field for a tuple `S` of primes.

    EXAMPLES::

        sage: K.<a> = QuadraticField(-14)
        sage: I = K.ideal(2, a)
        sage: S = (I,)
        sage: CS = K.S_class_group(S)
        sage: J = K.ideal(7, a)
        sage: G = K.ideal(3, a + 1)
        sage: CS(I)
        Trivial S-ideal class
        sage: CS(J)
        Trivial S-ideal class
        sage: CS(G)
        Fractional S-ideal class (3, a + 1)

    ::

        sage: K.<a> = QuadraticField(-14)
        sage: I = K.ideal(2, a)
        sage: S = (I,)
        sage: CS = K.S_class_group(S)
        sage: J = K.ideal(7, a)
        sage: G = K.ideal(3, a + 1)
        sage: CS(I).ideal()
        Fractional ideal (2, a)
        sage: CS(J).ideal()
        Fractional ideal (7, a)
        sage: CS(G).ideal()
        Fractional ideal (3, a + 1)

    ::

        sage: K.<a> = QuadraticField(-14)
        sage: I = K.ideal(2, a)
        sage: S = (I,)
        sage: CS = K.S_class_group(S)
        sage: G = K.ideal(3, a + 1)
        sage: CS(G).inverse()
        Fractional S-ideal class (3, a + 2)

    TESTS::

        sage: K.<a> = QuadraticField(-14)
        sage: I = K.ideal(2,a)
        sage: S = (I,)
        sage: CS = K.S_class_group(S)
        sage: J = K.ideal(7,a)
        sage: G = K.ideal(3,a+1)
        sage: CS(I).order()
        1
        sage: CS(J).order()
        1
        sage: CS(G).order()
        2
    """

    def _repr_(self):
        r"""
        Return a string representation of the `S`-ideal class of this fractional ideal.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: J = K.ideal(3, a + 2)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: CS(J)
            Fractional S-ideal class (3, a + 2)
            sage: CS(J^2)
            Trivial S-ideal class
        """
        if self.is_trivial():
            return 'Trivial S-ideal class'
        return 'Fractional S-ideal class %s' % self._value._repr_short()


class ClassGroup(AbelianGroupWithValues_class):
    r"""
    The class group of a number field.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^2 + 23)
        sage: G = K.class_group(); G
        Class group of order 3 with structure C3 of
         Number Field in a with defining polynomial x^2 + 23
        sage: G.category()
        Category of finite enumerated commutative groups

    Note the distinction between abstract generators, their ideal, and
    exponents::

        sage: C = NumberField(x^2 + 120071, 'a').class_group(); C
        Class group of order 500 with structure C250 x C2
        of Number Field in a with defining polynomial x^2 + 120071
        sage: c = C.gen(0)
        sage: c  # random
        Fractional ideal class (5, 1/2*a + 3/2)
        sage: c.ideal()  # random
        Fractional ideal (5, 1/2*a + 3/2)
        sage: c.ideal() is c.value()   # alias
        True
        sage: c.exponents()
        (1, 0)
    """
    Element = FractionalIdealClass

    def __init__(self, gens_orders, names, number_field, gens, proof=True):
        r"""
        Create a class group.

        TESTS::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 23)
            sage: G = K.class_group()
            sage: TestSuite(G).run()
        """
        AbelianGroupWithValues_class.__init__(self, gens_orders, names, gens,
                                              values_group=number_field.ideal_monoid())
        self._proof_flag = proof
        self._number_field = number_field

    def _element_constructor_(self, *args, **kwds):
        r"""
        Create an element of this class group from the given data. This may be:
        an ideal class in this number field; an ideal class in a subfield; or
        anything from which an ideal in this number field can be constructed.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<b> = NumberField(x^2 + 389)
            sage: C = K.class_group()
            sage: C(K.ideal(b)) # indirect doctest
            Trivial principal fractional ideal class
            sage: C(K.ideal(59049, b + 35312)) # indirect doctest
            Fractional ideal class (59049, b + 35312)
            sage: C((59049, b + 35312)) # indirect doctest
            Fractional ideal class (59049, b + 35312)
            sage: C(59049, b + 35312) # indirect doctest
            Fractional ideal class (59049, b + 35312)

            sage: K.<a> = QuadraticField(-23)
            sage: L.<b> = K.extension(x^2 - 2)
            sage: CK = K.class_group()
            sage: CL = L.class_group()
            sage: [CL(I).exponents() for I in CK]
            [(0,), (2,), (4,)]
        """
        if isinstance(args[0], FractionalIdealClass):
            return self.element_class(self, None, self._number_field.ideal(args[0].ideal()))
        else:
            I = self._number_field.ideal(*args, **kwds)
            if I.is_zero():
                raise TypeError("The zero ideal is not a fractional ideal")
            return self.element_class(self, None, I)

    def _ideal_log(self, ideal):
        """
        Compute the exponents from the ``ideal``.

        Used by the element constructor if necessary.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 + 23, 'a')
            sage: G = K.class_group()
            sage: g = G.an_element()
            sage: G._ideal_log(g.ideal())
            (1,)
            sage: g.exponents()
            (1,)
        """
        return tuple(ZZ(order) for order in ideal.ideal_class_log(proof=self._proof_flag))

    def gens_ideals(self):
        r"""
        Return generating ideals for the (`S`-)class group.

        This is an alias for :meth:`gens_values`.

        OUTPUT: a tuple of ideals, one for each abstract Abelian group generator

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^4 + 23)
            sage: K.class_group().gens_ideals()   # random gens (platform dependent)
            (Fractional ideal (2, 1/4*a^3 - 1/4*a^2 + 1/4*a - 1/4),)

            sage: C = NumberField(x^2 + x + 23899, 'a').class_group(); C
            Class group of order 68 with structure C34 x C2 of Number Field
            in a with defining polynomial x^2 + x + 23899
            sage: C.gens()
            (Fractional ideal class (7, a + 5), Fractional ideal class (5, a + 3))
            sage: C.gens_ideals()
            (Fractional ideal (7, a + 5), Fractional ideal (5, a + 3))
        """
        return self.gens_values()

    def __iter__(self):
        r"""
        Return an iterator of all ideal classes in this class group.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^4 + 23)
            sage: G = K.class_group()
            sage: G
            Class group of order 3 with structure C3 of Number Field
            in a with defining polynomial x^4 + 23
            sage: list(G)
            [Trivial principal fractional ideal class,
             Fractional ideal class (2, 1/4*a^3 - 1/4*a^2 + 1/4*a - 1/4),
             Fractional ideal class (2, 1/2*a^2 + 1/2)]
            sage: G.list()
            (Trivial principal fractional ideal class,
             Fractional ideal class (2, 1/4*a^3 - 1/4*a^2 + 1/4*a - 1/4),
             Fractional ideal class (2, 1/2*a^2 + 1/2))

        TESTS::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: G = K.class_group()
            sage: G
            Class group of order 1 of Number Field in a with defining polynomial x^2 + 1
            sage: list(G)
            [Trivial principal fractional ideal class]
            sage: G.list()
            (Trivial principal fractional ideal class,)
        """
        return self._iter_inner(self.one(), 0)

    def _iter_inner(self, i0, k):
        r"""
        Yield all elements of the coset `i0 * \{h in H_k\}`, where
        `H_k` is the subgroup of ``self`` generated by ``self.gens()[k:]``.

        Each new element provided costs exactly one group operation, and is
        not necessarily reduced.

        EXAMPLES::

            sage: x = ZZ['x'].gen()
            sage: K.<v> = NumberField(x^4 + 90*x^2 + 45)
            sage: OK = K.maximal_order()
            sage: G = OK.class_group()
            sage: iter = G._iter_inner(G.gen(0)^2,1)
            sage: all(next(iter) in G for _ in range(4))
            True
        """
        if k == self.ngens():
            yield i0
            return
        gk = self.gen(k)
        for _ in range(self._gens_orders[k]):
            yield from self._iter_inner(i0, k + 1)
            i0 = i0 * gk
        return

    def _repr_(self):
        r"""
        Return string representation of ``self``.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: C = NumberField(x^2 + 23, 'a').class_group()
            sage: C._repr_()
            'Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23'
        """
        s = 'Class group of order %s ' % self.order()
        if self.order() > 1:
            s += 'with structure %s ' % self._group_notation(self.gens_orders())
        s += 'of %s' % self.number_field()
        return s

    def number_field(self):
        r"""
        Return the number field that this (`S`-)class group is attached to.

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: C = NumberField(x^2 + 23, 'w').class_group(); C
            Class group of order 3 with structure C3 of
             Number Field in w with defining polynomial x^2 + 23
            sage: C.number_field()
            Number Field in w with defining polynomial x^2 + 23

            sage: K.<a> = QuadraticField(-14)
            sage: CS = K.S_class_group(K.primes_above(2))
            sage: CS.number_field()
            Number Field in a with defining polynomial x^2 + 14 with a = 3.741657386773942?*I
        """
        return self._number_field


class SClassGroup(ClassGroup):
    r"""
    The `S`-class group of a number field.

    EXAMPLES::

        sage: K.<a> = QuadraticField(-14)
        sage: S = K.primes_above(2)
        sage: K.S_class_group(S).gens()   # random gens (platform dependent)
        (Fractional S-ideal class (3, a + 2),)

        sage: K.<a> = QuadraticField(-974)
        sage: CS = K.S_class_group(K.primes_above(2)); CS
        S-class group of order 18 with structure C6 x C3 of
         Number Field in a with defining polynomial x^2 + 974 with a = 31.20897306865447?*I
        sage: CS.gen(0) # random
        Fractional S-ideal class (3, a + 2)
        sage: CS.gen(1) # random
        Fractional S-ideal class (31, a + 24)
    """
    Element = SFractionalIdealClass

    def __init__(self, gens_orders, names, number_field, gens, S, proof=True):
        r"""
        Create an `S`-class group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: K.S_class_group(S)
            S-class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 + 14 with a = 3.741657386773942?*I
            sage: K.<a> = QuadraticField(-105)
            sage: K.S_class_group([K.ideal(13, a + 8)])
            S-class group of order 4 with structure C2 x C2 of Number Field in a with defining polynomial x^2 + 105 with a = 10.24695076595960?*I
        """
        AbelianGroupWithValues_class.__init__(self, gens_orders, names, gens,
                                              values_group=number_field.ideal_monoid())
        self._proof_flag = proof
        self._number_field = number_field
        self._S = S

    def S(self):
        r"""
        Return the set (or rather tuple) of primes used to define this class group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2, a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S);CS
            S-class group of order 2 with structure C2 of
             Number Field in a with defining polynomial x^2 + 14 with a = 3.741657386773942?*I
            sage: T = tuple()
            sage: CT = K.S_class_group(T);CT
            S-class group of order 4 with structure C4 of
             Number Field in a with defining polynomial x^2 + 14 with a = 3.741657386773942?*I
            sage: CS.S()
            (Fractional ideal (2, a),)
            sage: CT.S()
            ()
        """
        return self._S

    def _ideal_log(self, ideal):
        """
        Compute the exponents from the ``ideal``.

        Used by the element constructor if necessary.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: s = CS.an_element()
            sage: CS._ideal_log(s.ideal())
            (1,)
            sage: s.exponents()
            (1,)
        """
        return tuple(ZZ(order) for order in ideal.S_ideal_class_log(self.S()))

    def _element_constructor_(self, *args, **kwds):
        r"""
        Create an element of this class group from the given data.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: I = K.ideal(2,a)
            sage: S = (I,)
            sage: CS = K.S_class_group(S)
            sage: J = K.ideal(7,a)
            sage: G = K.ideal(3,a+1)
            sage: CS(I)
            Trivial S-ideal class
            sage: CS(J)
            Trivial S-ideal class
            sage: CS(G)
            Fractional S-ideal class (3, a + 1)
        """
        if isinstance(args[0], FractionalIdealClass):
            return self.element_class(self, None, args[0].ideal())
        else:
            I = self.number_field().ideal(*args, **kwds)
            if I.is_zero():
                raise TypeError("The zero ideal is not a fractional ideal")
            return self.element_class(self, None, I)

    def _repr_(self):
        r"""
        Return string representation of this S-class group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-14)
            sage: CS = K.S_class_group(K.primes_above(2))
            sage: CS._repr_()
            'S-class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 + 14 with a = 3.741657386773942?*I'
        """
        s = 'S-class group of order %s ' % self.order()
        if self.order() > 1:
            s += 'with structure %s ' % self._group_notation(self.gens_orders())
        s += 'of %s' % self.number_field()
        return s
