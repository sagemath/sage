r"""
Field `\QQ` of Rational Numbers

The class :class:`RationalField` represents the field `\QQ` of (arbitrary
precision) rational numbers. Each rational number is an instance of the class
:class:`Rational`.

Interactively, an instance of :class:`RationalField` is available as ``QQ``::

    sage: QQ
    Rational Field

Values of various types can be converted to rational numbers by using the
:meth:`__call__` method of :class:`RationalField` (that is, by treating ``QQ`` as a
function).

::

    sage: RealField(9).pi()                                                             # needs sage.rings.real_mpfr
    3.1
    sage: QQ(RealField(9).pi())                                                         # needs sage.rings.real_interval_field sage.rings.real_mpfr
    22/7
    sage: QQ(RealField().pi())                                                          # needs sage.rings.real_interval_field sage.rings.real_mpfr
    245850922/78256779
    sage: QQ(35)
    35
    sage: QQ('12/347')
    12/347
    sage: QQ(exp(pi*I))                                                                 # needs sage.symbolic
    -1
    sage: x = polygen(ZZ)
    sage: QQ((3*x)/(4*x))
    3/4

TESTS::

    sage: Q = RationalField()
    sage: Q == loads(dumps(Q))
    True
    sage: RationalField() is RationalField()
    True
    sage: Q in Fields().Infinite()
    True

AUTHORS:

- Niles Johnson (2010-08): :issue:`3893`: ``random_element()`` should pass on
  ``*args`` and ``**kwds``.

- Travis Scrimshaw (2012-10-18): Added additional docstrings for full coverage.
  Removed duplicates of ``discriminant()`` and ``signature()``.

- Anna Haensch (2018-03): Added function ``quadratic_defect()``
"""

from sage.rings.integer import Integer
from sage.rings.rational import Rational

ZZ = None

import sage.rings.number_field.number_field_base as number_field_base
from sage.misc.fast_methods import Singleton
from sage.misc.superseded import deprecated_function_alias
from sage.structure.parent import Parent
from sage.structure.sequence import Sequence


class RationalField(Singleton, number_field_base.NumberField):
    r"""
    The class ``RationalField`` represents the field `\QQ` of rational numbers.

    EXAMPLES::

        sage: a = 901824309821093821093812093810928309183091832091
        sage: b = QQ(a); b
        901824309821093821093812093810928309183091832091
        sage: QQ(b)
        901824309821093821093812093810928309183091832091
        sage: QQ(int(93820984323))
        93820984323
        sage: QQ(ZZ(901824309821093821093812093810928309183091832091))
        901824309821093821093812093810928309183091832091
        sage: QQ('-930482/9320842317')
        -930482/9320842317
        sage: QQ((-930482, 9320842317))
        -930482/9320842317
        sage: QQ([9320842317])
        9320842317
        sage: QQ(pari(39029384023840928309482842098430284398243982394))                 # needs sage.libs.pari
        39029384023840928309482842098430284398243982394
        sage: QQ('sage')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 'sage' to a rational

    Conversion from the reals to the rationals is done by default using
    continued fractions.

    ::

        sage: QQ(RR(3929329/32))
        3929329/32
        sage: QQ(-RR(3929329/32))
        -3929329/32
        sage: QQ(RR(1/7)) - 1/7                                                         # needs sage.rings.real_mpfr
        0

    If you specify the optional second argument ``base``, then the string
    representation of the float is used.

    ::

        sage: # needs sage.rings.real_mpfr
        sage: QQ(23.2, 2)
        6530219459687219/281474976710656
        sage: 6530219459687219.0/281474976710656
        23.20000000000000
        sage: a = 23.2; a
        23.2000000000000
        sage: QQ(a, 10)
        116/5

    Here's a nice example involving elliptic curves::

        sage: # needs sage.rings.real_mpfr sage.schemes
        sage: E = EllipticCurve('11a')
        sage: L = E.lseries().at1(300)[0]; L
        0.2538418608559106843377589233...
        sage: O = E.period_lattice().omega(); O
        1.26920930427955
        sage: t = L/O; t
        0.200000000000000
        sage: QQ(RealField(45)(t))
        1/5
    """
    _zero_element = Rational(0)
    _one_element = Rational(1)

    def __new__(cls):
        """
        This method actually is not needed for using :class:`RationalField`.
        But it is used to unpickle some very old pickles.

        TESTS::

            sage: RationalField() in Fields() # indirect doctest
            True
        """
        try:
            return QQ
        except BaseException:
            from sage.rings.number_field.number_field_base import NumberField
            return NumberField.__new__(cls)

    def __init__(self):
        r"""
        We create the rational numbers `\QQ`, and call a few functions::

            sage: Q = RationalField(); Q
            Rational Field
            sage: Q.characteristic()
            0
            sage: Q.is_field()
            True
            sage: Q.category()
            Join of Category of number fields
             and Category of quotient fields
             and Category of metric spaces
            sage: Q.zeta()
            -1

        We next illustrate arithmetic in `\QQ`.

        ::

            sage: Q('49/7')
            7
            sage: type(Q('49/7'))
            <class 'sage.rings.rational.Rational'>
            sage: a = Q('19/374'); a
            19/374
            sage: b = Q('17/371'); b
            17/371
            sage: a + b
            13407/138754
            sage: b + a
            13407/138754
            sage: a * b
            19/8162
            sage: b * a
            19/8162
            sage: a - b
            691/138754
            sage: b - a
            -691/138754
            sage: a / b
            7049/6358
            sage: b / a
            6358/7049
            sage: b < a
            True
            sage: a < b
            False

        Next finally illustrate arithmetic with automatic coercion. The
        types that coerce into the rational field include ``str, int,
        long, Integer``.

        ::

            sage: a + Q('17/371')
            13407/138754
            sage: a * 374
            19
            sage: 374 * a
            19
            sage: a/19
            1/374
            sage: a + 1
            393/374

        TESTS::

            sage: TestSuite(QQ).run()
            sage: QQ.variable_name()
            'x'
            sage: QQ.variable_names()
            ('x',)
            sage: QQ._element_constructor_((2, 3))
            2/3

            sage: QQ.is_finite()
            False

            sage: QQ.is_field()
            True
        """
        from sage.categories.basic import QuotientFields
        from sage.categories.number_fields import NumberFields
        Parent.__init__(self, base=self,
                        category=[QuotientFields().Metric(),
                                  NumberFields()])
        self._assign_names(('x',), normalize=False)  # ?????
        self._populate_coercion_lists_(init_no_parent=True)

    _element_constructor_ = Rational

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: QQ # indirect doctest
            Rational Field
        """
        return "Rational Field"

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: QQ._repr_option('element_is_atomic')
            True
        """
        if key == 'element_is_atomic':
            return True
        return super()._repr_option(key)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(QQ) # indirect doctest
            \Bold{Q}
        """
        return r"\Bold{Q}"

    def __reduce__(self):
        r"""
        Used for pickling `\QQ`.

        EXAMPLES::

           sage: loads(dumps(QQ)) is QQ
           True
        """
        return RationalField, tuple()

    def __len__(self):
        """
        Return the number of elements in ``self``.

        Since this does not have a size, this throws a :exc:`TypeError`.

        EXAMPLES::

            sage: len(QQ)
            Traceback (most recent call last):
            ...
            TypeError: len() of unsized object
        """
        raise TypeError('len() of unsized object')

    def construction(self):
        r"""
        Return a pair ``(functor, parent)`` such that ``functor(parent)``
        returns ``self``.

        This is the construction of `\QQ` as the fraction field of `\ZZ`.

        EXAMPLES::

            sage: QQ.construction()
            (FractionField, Integer Ring)
        """
        from sage.categories.pushout import FractionField

        from . import integer_ring
        return FractionField(), integer_ring.ZZ

    def completion(self, p, prec, extras={}):
        r"""
        Return the completion of `\QQ` at `p`.

        EXAMPLES::

            sage: QQ.completion(infinity, 53)                                           # needs sage.rings.real_mpfr
            Real Field with 53 bits of precision
            sage: QQ.completion(5, 15, {'print_mode': 'bars'})                          # needs sage.rings.padics
            5-adic Field with capped relative precision 15
        """
        from sage.rings.infinity import Infinity
        if p == Infinity:
            from sage.rings.real_field import create_RealField
            return create_RealField(prec, **extras)
        else:
            from sage.rings.padics.factory import Qp
            return Qp(p, prec, **extras)

    def _coerce_map_from_(self, S):
        """
        Return a coerce map from ``S``.

        EXAMPLES::

            sage: f = QQ.coerce_map_from(ZZ); f # indirect doctest
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
            sage: f(3)
            3
            sage: f(3^99) - 3^99
            0
            sage: f = QQ.coerce_map_from(int); f # indirect doctest
            Native morphism:
              From: Set of Python objects of class 'int'
              To:   Rational Field
            sage: f(44)
            44

        ::

            sage: L = Localization(ZZ, (3,5))
            sage: 1/45 in L  # indirect doctest
            True
            sage: 1/43 in L  # indirect doctest
            False
        """
        global ZZ
        from . import rational
        if ZZ is None:
            from . import integer_ring
            ZZ = integer_ring.ZZ
        if S is ZZ:
            return rational.Z_to_Q()
        elif S is int:
            return rational.int_to_Q()
        elif ZZ.has_coerce_map_from(S):
            return rational.Z_to_Q() * ZZ._internal_coerce_map_from(S)
        from sage.rings.localization import Localization
        if isinstance(S, Localization):
            if S.fraction_field() is self:
                from sage.structure.coerce_maps import CallableConvertMap
                return CallableConvertMap(S, self, lambda x: x._value, parent_as_first_arg=False)

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        Check to see if the map into ``codomain`` determined by ``im_gens`` is
        a valid homomorphism.

        EXAMPLES::

            sage: QQ._is_valid_homomorphism_(ZZ, [1])
            False
            sage: QQ._is_valid_homomorphism_(QQ, [1])
            True
            sage: QQ._is_valid_homomorphism_(RR, [1])
            True
            sage: QQ._is_valid_homomorphism_(RR, [2])
            False
        """
        try:
            return im_gens[0] == codomain.coerce(self.gen(0))
        except TypeError:
            return False

    def __iter__(self):
        r"""
        Create an iterator that generates the rational numbers without
        repetition, in order of the height.

        See also :meth:`range_by_height()`.

        EXAMPLES:

        The first 17 rational numbers, ordered by height::

            sage: import itertools
            sage: lst = [a for a in itertools.islice(Rationals(), 17r)]
            sage: lst
            [0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3, -3, 2/3, -2/3, 3/2, -3/2, 1/4, -1/4]
            sage: [a.height() for a in lst]
            [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4]
        """
        yield self(0)
        yield self(1)
        yield self(-1)
        height = Integer(1)
        while True:
            height = height + 1
            for other in range(1, height):
                if height.gcd(other) == 1:
                    yield self(other/height)
                    yield self(-other/height)
                    yield self(height/other)
                    yield self(-height/other)

    def __truediv__(self, I):
        """
        Form the quotient by an integral ideal.

        EXAMPLES::

            sage: QQ / ZZ                                                               # needs sage.modules
            Q/Z
        """
        from sage.groups.additive_abelian.qmodnz import QmodnZ
        from sage.rings.ideal import Ideal_generic
        if I is ZZ:
            return QmodnZ(1)
        elif isinstance(I, Ideal_generic) and I.base_ring() is ZZ:
            return QmodnZ(I.gen())
        else:
            return super().__truediv__(I)

    def range_by_height(self, start, end=None):
        r"""
        Range function for rational numbers, ordered by height.

        Returns a Python generator for the list of rational numbers with
        heights in ``range(start, end)``. Follows the same
        convention as Python :func:`range`, type ``range?`` for details.

        See also :meth:`__iter__`.

        EXAMPLES:

        All rational numbers with height strictly less than 4::

            sage: list(QQ.range_by_height(4))
            [0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3, -3, 2/3, -2/3, 3/2, -3/2]
            sage: [a.height() for a in QQ.range_by_height(4)]
            [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3]

        All rational numbers with height 2::

            sage: list(QQ.range_by_height(2, 3))
            [1/2, -1/2, 2, -2]

        Nonsensical integer arguments will return an empty generator::

            sage: list(QQ.range_by_height(3, 3))
            []
            sage: list(QQ.range_by_height(10, 1))
            []

        There are no rational numbers with height `\leq 0`::

            sage: list(QQ.range_by_height(-10, 1))
            []
        """
        if end is None:
            end = start
            start = 1
        start = max(start, 1)
        for height in ZZ.range(start, end):
            if height == 1:
                yield self(0)
                yield self(1)
                yield self(-1)
            for other in ZZ.range(1, height):
                if height.gcd(other) == 1:
                    yield self(other/height)
                    yield self(-other/height)
                    yield self(height/other)
                    yield self(-height/other)

    def primes_of_bounded_norm_iter(self, B):
        r"""
        Iterator yielding all primes less than or equal to `B`.

        INPUT:

        - ``B`` -- positive integer; upper bound on the primes generated

        OUTPUT: an iterator over all integer primes less than or equal to `B`

        .. NOTE::

            This function exists for compatibility with the related number
            field method, though it returns prime integers, not ideals.

        EXAMPLES::

            sage: it = QQ.primes_of_bounded_norm_iter(10)
            sage: list(it)                                                              # needs sage.libs.pari
            [2, 3, 5, 7]
            sage: list(QQ.primes_of_bounded_norm_iter(1))
            []
        """
        try:
            B = ZZ(B.ceil())
        except (TypeError, AttributeError):
            raise TypeError("%s is not valid bound on prime ideals" % B)

        if B < 2:
            return

        from sage.arith.misc import primes
        yield from primes(B + 1)

    def discriminant(self):
        """
        Return the discriminant of the field of rational numbers, which is 1.

        EXAMPLES::

            sage: QQ.discriminant()
            1
        """
        return Integer(1)

    def absolute_discriminant(self):
        """
        Return the absolute discriminant, which is 1.

        EXAMPLES::

            sage: QQ.absolute_discriminant()
            1
        """
        return self.discriminant()

    def relative_discriminant(self):
        """
        Return the relative discriminant, which is 1.

        EXAMPLES::

            sage: QQ.relative_discriminant()
            1
        """
        return self.discriminant()

    def class_number(self):
        """
        Return the class number of the field of rational numbers, which is 1.

        EXAMPLES::

            sage: QQ.class_number()
            1
        """
        return Integer(1)

    def signature(self):
        r"""
        Return the signature of the rational field, which is `(1,0)`, since
        there are 1 real and no complex embeddings.

        EXAMPLES::

            sage: QQ.signature()
            (1, 0)
        """
        return (Integer(1), Integer(0))

    def embeddings(self, K):
        r"""
        Return the list containing the unique embedding of `\QQ` into `K`,
        if it exists, and an empty list otherwise.

        EXAMPLES::

            sage: QQ.embeddings(QQ)
            [Identity endomorphism of Rational Field]
            sage: QQ.embeddings(CyclotomicField(5))                                     # needs sage.rings.number_field
            [Coercion map:
               From: Rational Field
               To:   Cyclotomic Field of order 5 and degree 4]

        The field `K` must have characteristic `0` for an embedding of `\QQ`
        to exist::

            sage: QQ.embeddings(GF(3))
            []
        """
        if K.characteristic() == 0:
            v = [self.hom(K)]
        else:
            v = []
        return Sequence(v, check=False, universe=self.Hom(K))

    def automorphisms(self):
        r"""
        Return all Galois automorphisms of ``self``.

        OUTPUT: a sequence containing just the identity morphism

        EXAMPLES::

            sage: QQ.automorphisms()
            [Ring endomorphism of Rational Field
               Defn: 1 |--> 1]
        """
        return Sequence([self.hom(1, self)], cr=True, immutable=False,
                        check=False)

    def places(self, all_complex=False, prec=None):
        r"""
        Return the collection of all infinite places of ``self``, which
        in this case is just the embedding of ``self`` into `\RR`.

        By default, this returns homomorphisms into ``RR``.  If
        ``prec`` is not None, we simply return homomorphisms into
        ``RealField(prec)`` (or ``RDF`` if ``prec=53``).

        There is an optional flag ``all_complex``, which defaults to
        ``False``.  If ``all_complex`` is ``True``, then the real embeddings
        are returned as embeddings into the corresponding complex
        field.

        For consistency with non-trivial number fields.

        EXAMPLES::

            sage: QQ.places()                                                           # needs sage.rings.real_mpfr
            [Ring morphism:
              From: Rational Field
              To:   Real Field with 53 bits of precision
              Defn: 1 |--> 1.00000000000000]
            sage: QQ.places(prec=53)
            [Ring morphism:
              From: Rational Field
              To:   Real Double Field
              Defn: 1 |--> 1.0]
            sage: QQ.places(prec=200, all_complex=True)                                 # needs sage.rings.real_mpfr
            [Ring morphism:
              From: Rational Field
              To:   Complex Field with 200 bits of precision
              Defn: 1 |--> 1.0000000000000000000000000000000000000000000000000000000000]
        """
        from sage.rings.infinity import Infinity
        if prec is None:
            if all_complex:
                from sage.rings.cc import CC as domain
            else:
                from sage.rings.real_mpfr import RR as domain
        elif prec == 53:
            if all_complex:
                from sage.rings.complex_double import CDF as domain
            else:
                from sage.rings.real_double import RDF as domain
        elif prec == Infinity:
            if all_complex:
                from sage.rings.qqbar import QQbar as domain
            else:
                from sage.rings.qqbar import AA as domain
        else:
            if all_complex:
                from sage.rings.complex_mpfr import ComplexField
                domain = ComplexField(prec)
            else:
                from sage.rings.real_mpfr import RealField
                domain = RealField(prec)
        return [self.hom([domain(1)])]

    def complex_embedding(self, prec=53):
        """
        Return embedding of the rational numbers into the complex numbers.

        EXAMPLES::

            sage: QQ.complex_embedding()                                                # needs sage.rings.real_mpfr
            Ring morphism:
              From: Rational Field
              To:   Complex Field with 53 bits of precision
              Defn: 1 |--> 1.00000000000000
            sage: QQ.complex_embedding(20)                                              # needs sage.rings.real_mpfr
            Ring morphism:
              From: Rational Field
              To:   Complex Field with 20 bits of precision
              Defn: 1 |--> 1.0000
        """
        from . import complex_mpfr
        CC = complex_mpfr.ComplexField(prec)
        return self.hom([CC(1)])

    def residue_field(self, p, check=True):
        r"""
        Return the residue field of `\QQ` at the prime `p`, for
        consistency with other number fields.

        INPUT:

        - ``p`` -- prime integer

        - ``check`` -- (default: ``True``) if ``True``, check the primality of
           `p`, else do not

        OUTPUT: the residue field at this prime

        EXAMPLES::

            sage: QQ.residue_field(5)
            Residue field of Integers modulo 5
            sage: QQ.residue_field(next_prime(10^9))                                    # needs sage.rings.finite_rings
            Residue field of Integers modulo 1000000007
        """
        from sage.rings.finite_rings.residue_field import ResidueField
        return ResidueField(ZZ.ideal(p), check=check)

    def hilbert_symbol_negative_at_S(self, S, b, check=True):
        r"""
        Return an integer that has a negative Hilbert symbol with respect
        to a given rational number and a given set of primes (or places).

        The function is algorithm 3.4.1 in [Kir2016]_. It finds an integer `a`
        that has negative Hilbert symbol with respect to a given rational number
        exactly at a given set of primes (or places).

        INPUT:

        - ``S`` -- list of rational primes, the infinite place as real
          embedding of `\QQ` or as `-1`
        - ``b`` -- a nonzero rational number which is a non-square locally
          at every prime in ``S``
        - ``check`` -- boolean (default: ``True``); perform additional checks on
          input and confirm the output

        OUTPUT:

        - An integer `a` that has negative Hilbert symbol `(a,b)_p` for
          every place `p` in `S` and no other place.

        EXAMPLES::

            sage: QQ.hilbert_symbol_negative_at_S([-1,5,3,2,7,11,13,23], -10/7)         # needs sage.rings.padics
            -9867
            sage: QQ.hilbert_symbol_negative_at_S([3, 5, QQ.places()[0], 11], -15)      # needs sage.rings.padics
            -33
            sage: QQ.hilbert_symbol_negative_at_S([3, 5], 2)                            # needs sage.rings.padics
            15

        TESTS::

            sage: QQ.hilbert_symbol_negative_at_S(5/2, -2)                              # needs sage.modules
            Traceback (most recent call last):
            ...
            TypeError: first argument must be a list or integer

        ::

            sage: QQ.hilbert_symbol_negative_at_S([1, 3], 0)                            # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: second argument must be nonzero

        ::

            sage: QQ.hilbert_symbol_negative_at_S([-1, 3, 5], 2)                        # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: list should be of even cardinality

        ::

            sage: QQ.hilbert_symbol_negative_at_S([1, 3], 2)                            # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: all entries in list must be prime or -1 for
            infinite place

        ::

            sage: QQ.hilbert_symbol_negative_at_S([5, 7], 2)                            # needs sage.libs.pari sage.modules sage.rings.padics
            Traceback (most recent call last):
            ...
            ValueError: second argument must be a nonsquare with
            respect to every finite prime in the list

        ::

            sage: QQ.hilbert_symbol_negative_at_S([1, 3], sqrt(2))                      # needs sage.libs.pari sage.modules sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: second argument must be a rational number

        ::

            sage: QQ.hilbert_symbol_negative_at_S([-1, 3], 2)                           # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: if the infinite place is in the list, the second
            argument must be negative

        AUTHORS:

        - Simon Brandhorst, Juanita Duque, Anna Haensch, Manami Roy, Sandi Rudzinski (10-24-2017)
        """
        from sage.arith.misc import hilbert_symbol, is_prime
        from sage.matrix.constructor import matrix
        from sage.modules.free_module import VectorSpace
        from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
        from sage.rings.padics.factory import Qp
        from sage.sets.primes import Primes

        # input checks
        if type(S) is not list:
            raise TypeError("first argument must be a list or integer")
        # -1 is used for the infinite place
        infty = -1
        for i in range(len(S)):
            if S[i] == self.places()[0]:
                S[i] = -1
        if b not in self:
            raise TypeError("second argument must be a rational number")
        b = self(b)
        if b == 0:
            raise ValueError("second argument must be nonzero")
        if len(S) % 2:
            raise ValueError("list should be of even cardinality")
        for p in S:
            if p != infty:
                if check and not is_prime(p):
                    raise ValueError("all entries in list must be prime"
                                     " or -1 for infinite place")
                R = Qp(p)
                if R(b).is_square():
                    raise ValueError("second argument must be a nonsquare with"
                                     " respect to every finite prime in the list")
            elif b > 0:
                raise ValueError("if the infinite place is in the list, "
                                 "the second argument must be negative")
        # L is the list of primes that we need to consider, b must have
        # nonzero valuation for each prime in L, this is the set S'
        # in Kirschmer's algorithm
        L = []
        L = [p[0] for p in b.factor() if p[0] not in S]
        # We must also consider 2 to be in L
        if 2 not in L and 2 not in S:
            L.append(2)
        # This adds the infinite place to L
        if b < 0 and infty not in S:
            L.append(infty)

        P = S + L
        # This constructs the vector v in the algorithm. This is the vector
        # that we are searching for. It represents the case when the Hilbert
        # symbol is negative for all primes in S and positive
        # at all primes in S'
        V = VectorSpace(GF(2), len(P))
        v = V([1]*len(S) + [0]*len(L))

        # Compute the map phi of Hilbert symbols at all the primes
        # in S and S'
        # For technical reasons, a Hilbert symbol of -1 is
        # respresented as 1 and a Hilbert symbol of 1
        # is represented as 0
        def phi(x):
            v = [(1-hilbert_symbol(x, b, p))//2 for p in P]
            return V(v)

        M = matrix(GF(2), [phi(p) for p in P+[-1]])
        # We search through all the primes
        for q in Primes():
            # Only look at this prime if it is not in our list
            if q in P:
                continue

            # The algorithm terminates when the vector v is in the
            # subspace of V generated by the image of the phi map
            # on the set of generators
            w = phi(q)
            W = M.stack(matrix(w))
            if v in W.row_space():
                break
        Pq = P + [-1] + [q]
        l = W.solve_left(v)
        a = self.prod([Pq[i]**ZZ(l[i]) for i in range(l.degree())])
        if check:
            assert phi(a) == v, "oops"
        return a

    def gens(self) -> tuple:
        r"""
        Return a tuple of generators of `\QQ`, which is only ``(1,)``.

        EXAMPLES::

            sage: QQ.gens()
            (1,)
        """
        return (self(1), )

    def gen(self, n=0):
        r"""
        Return the n-th generator of `\QQ`.

        There is only the 0-th generator, which is 1.

        EXAMPLES::

            sage: QQ.gen()
            1
        """
        if n == 0:
            return self(1)
        else:
            raise IndexError("n must be 0")

    def degree(self):
        r"""
        Return the degree of `\QQ`, which is 1.

        EXAMPLES::

            sage: QQ.degree()
            1
        """
        return Integer(1)

    def absolute_degree(self):
        r"""
        Return the absolute degree of `\QQ`, which is 1.

        EXAMPLES::

            sage: QQ.absolute_degree()
            1
        """
        return Integer(1)

    def ngens(self):
        r"""
        Return the number of generators of `\QQ`, which is 1.

        EXAMPLES::

            sage: QQ.ngens()
            1
        """
        return Integer(1)

    def is_absolute(self):
        r"""
        `\QQ` is an absolute extension of `\QQ`.

        EXAMPLES::

            sage: QQ.is_absolute()
            True
        """
        return True

    def is_prime_field(self):
        r"""
        Return ``True`` since `\QQ` is a prime field.

        EXAMPLES::

            sage: QQ.is_prime_field()
            True
        """
        return True

    def characteristic(self):
        r"""
        Return 0 since the rational field has characteristic 0.

        EXAMPLES::

            sage: c = QQ.characteristic(); c
            0
            sage: parent(c)
            Integer Ring
        """
        return Integer(0)

    def maximal_order(self):
        r"""
        Return the maximal order of the rational numbers, i.e., the ring
        `\ZZ` of integers.

        EXAMPLES::

            sage: QQ.maximal_order()
            Integer Ring
            sage: QQ.ring_of_integers ()
            Integer Ring
        """
        from .integer_ring import ZZ
        return ZZ

    def number_field(self):
        r"""
        Return the number field associated to `\QQ`. Since `\QQ` is a number
        field, this just returns `\QQ` again.

        EXAMPLES::

            sage: QQ.number_field() is QQ
            True
        """
        return self

    def power_basis(self):
        r"""
        Return a power basis for this number field over its base field.

        The power basis is always ``[1]`` for the rational field. This method
        is defined to make the rational field behave more like a number
        field.

        EXAMPLES::

            sage: QQ.power_basis()
            [1]
        """
        return [self.gen()]

    def extension(self, poly, names, **kwds):
        r"""
        Create a field extension of `\QQ`.

        EXAMPLES:

        We make a single absolute extension::

            sage: x = polygen(QQ, 'x')
            sage: K.<a> = QQ.extension(x^3 + 5); K                                      # needs sage.rings.number_field
            Number Field in a with defining polynomial x^3 + 5

        We make an extension generated by roots of two polynomials::

            sage: K.<a,b> = QQ.extension([x^3 + 5, x^2 + 3]); K                         # needs sage.rings.number_field
            Number Field in a with defining polynomial x^3 + 5 over its base field
            sage: b^2                                                                   # needs sage.rings.number_field
            -3
            sage: a^3                                                                   # needs sage.rings.number_field
            -5
        """
        from sage.rings.number_field.number_field import NumberField
        return NumberField(poly, names=names, **kwds)

    def algebraic_closure(self):
        r"""
        Return the algebraic closure of ``self`` (which is `\QQbar`).

        EXAMPLES::

            sage: QQ.algebraic_closure()                                                # needs sage.rings.number_field
            Algebraic Field
        """
        from sage.rings.qqbar import QQbar
        return QQbar

    def order(self):
        r"""
        Return the order of `\QQ`, which is `\infty`.

        EXAMPLES::

            sage: QQ.order()
            +Infinity
        """
        from sage.rings.infinity import Infinity
        return Infinity

    def polynomial(self):
        r"""
        Return a defining polynomial of `\QQ`, as for other number fields.

        This is also aliased to :meth:`defining_polynomial`
        and :meth:`absolute_polynomial`.

        EXAMPLES::

            sage: QQ.polynomial()
            x
        """
        from sage.rings.polynomial.polynomial_ring import polygen
        return polygen(self)

    defining_polynomial = polynomial
    absolute_polynomial = polynomial

    def _an_element_(self):
        r"""
        Return an element of `\QQ`.

        EXAMPLES::

            sage: QQ.an_element() # indirect doctest
            1/2
        """
        return Rational((1, 2))

    def some_elements(self):
        r"""
        Return some elements of `\QQ`.

        See :func:`TestSuite` for a typical use case.

        OUTPUT: an iterator over 100 elements of `\QQ`

        EXAMPLES::

            sage: tuple(QQ.some_elements())
            (1/2, -1/2, 2, -2,
             0, 1, -1, 42,
             2/3, -2/3, 3/2, -3/2,
             4/5, -4/5, 5/4, -5/4,
             6/7, -6/7, 7/6, -7/6,
             8/9, -8/9, 9/8, -9/8,
             10/11, -10/11, 11/10, -11/10,
             12/13, -12/13, 13/12, -13/12,
             14/15, -14/15, 15/14, -15/14,
             16/17, -16/17, 17/16, -17/16,
             18/19, -18/19, 19/18, -19/18,
             20/441, -20/441, 441/20, -441/20,
             22/529, -22/529, 529/22, -529/22,
             24/625, -24/625, 625/24, -625/24,
             ...)
        """
        yield self.an_element()
        yield -self.an_element()
        yield 1/self.an_element()
        yield -1/self.an_element()
        yield self(0)
        yield self(1)
        yield self(-1)
        yield self(42)
        for n in range(1, 24):
            a = 2*n
            b = (2*n + 1)**(n//10 + 1)
            yield Rational((a, b))
            yield Rational((-a, b))
            yield Rational((b, a))
            yield Rational((-b, a))

    def random_element(self, num_bound=None, den_bound=None, *args, **kwds):
        r"""
        Return a random element of `\QQ`.

        Elements are constructed by randomly choosing integers
        for the numerator and denominator, not necessarily coprime.

        INPUT:

        - ``num_bound`` -- positive integer, specifying a bound
          on the absolute value of the numerator.
          If absent, no bound is enforced.

        - ``den_bound`` -- positive integer, specifying a bound
          on the value of the denominator.
          If absent, the bound for the numerator will be reused.

        Any extra positional or keyword arguments are passed through to
        :meth:`sage.rings.integer_ring.IntegerRing_class.random_element`.

        EXAMPLES::

            sage: QQ.random_element().parent() is QQ
            True
            sage: while QQ.random_element() != 0:
            ....:     pass
            sage: while QQ.random_element() != -1/2:
            ....:     pass

        In the following example, the resulting numbers range from
        -5/1 to 5/1 (both inclusive),
        while the smallest possible positive value is 1/10::

            sage: q = QQ.random_element(5, 10)
            sage: -5/1 <= q <= 5/1
            True
            sage: q.denominator() <= 10
            True
            sage: q.numerator() <= 5
            True

        Extra positional or keyword arguments are passed through::

            sage: QQ.random_element(distribution='1/n').parent() is QQ
            True
            sage: QQ.random_element(distribution='1/n').parent() is QQ
            True
        """
        global ZZ
        if ZZ is None:
            from . import integer_ring
            ZZ = integer_ring.ZZ
        if num_bound is None:
            num = ZZ.random_element(*args, **kwds)
            den = ZZ.random_element(*args, **kwds)
            while den == 0:
                den = ZZ.random_element(*args, **kwds)
            return self((num, den))
        else:
            if num_bound == 0:
                num_bound = 2
            if den_bound is None:
                den_bound = num_bound
                if den_bound < 1:
                    den_bound = 2
            num = ZZ.random_element(-num_bound, num_bound+1, *args, **kwds)
            den = ZZ.random_element(1, den_bound+1, *args, **kwds)
            while den == 0:
                den = ZZ.random_element(1, den_bound+1, *args, **kwds)
            return self((num, den))

    def zeta(self, n=2):
        """
        Return a root of unity in ``self``.

        INPUT:

        - ``n`` -- integer (default: 2); order of the root of unity

        EXAMPLES::

            sage: QQ.zeta()
            -1
            sage: QQ.zeta(2)
            -1
            sage: QQ.zeta(1)
            1
            sage: QQ.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: no n-th root of unity in rational field
        """
        if n == 1:
            return Rational(1)
        elif n == 2:
            return Rational(-1)
        else:
            raise ValueError("no n-th root of unity in rational field")

    def selmer_generators(self, S, m, proof=True, orders=False):
        r"""
        Return generators of the group `\QQ(S,m)`.

        INPUT:

        - ``S`` -- set of primes

        - ``m`` -- positive integer

        - ``proof`` -- ignored

        - ``orders`` -- (default: ``False``) if ``True``, output two lists, the
          generators and their orders

        OUTPUT:

        A list of generators of `\QQ(S,m)` (and, optionally, their
        orders in `\QQ^\times/(\QQ^\times)^m`).  This is the subgroup
        of `\QQ^\times/(\QQ^\times)^m` consisting of elements `a` such
        that the valuation of `a` is divisible by `m` at all primes
        not in `S`.  It is equal to the group of `S`-units modulo
        `m`-th powers.  The group `\QQ(S,m)` contains the subgroup of
        those `a` such that `\QQ(\sqrt[m]{a})/\QQ` is unramified at
        all primes of `\QQ` outside of `S`, but may contain it
        properly when not all primes dividing `m` are in `S`.

        .. SEEALSO::

            :meth:`RationalField.selmer_space`, which gives additional
            output when `m=p` is prime: as well as generators, it gives an
            abstract vector space over `\GF{p}` isomorphic to `\QQ(S,p)`
            and maps implementing the isomorphism between this space and
            `\QQ(S,p)` as a subgroup of `\QQ^*/(\QQ^*)^p`.

        EXAMPLES::

            sage: QQ.selmer_generators((), 2)
            [-1]
            sage: QQ.selmer_generators((3,), 2)
            [-1, 3]
            sage: QQ.selmer_generators((5,), 2)
            [-1, 5]

        The previous examples show that the group generated by the
        output may be strictly larger than the 'true' Selmer group of
        elements giving extensions unramified outside `S`.

        When `m` is even, `-1` is a generator of order `2`::

            sage: QQ.selmer_generators((2,3,5,7,), 2, orders=True)
            ([-1, 2, 3, 5, 7], [2, 2, 2, 2, 2])
            sage: QQ.selmer_generators((2,3,5,7,), 3, orders=True)
            ([2, 3, 5, 7], [3, 3, 3, 3])
        """
        gens = list(S)
        ords = [ZZ(m)] * len(S)
        if m % 2 == 0:
            gens = [ZZ(-1)] + gens
            ords = [ZZ(2)] + ords
        if orders:
            return gens, ords
        else:
            return gens

    # For backwards compatibility:
    selmer_group = deprecated_function_alias(31345, selmer_generators)

    def selmer_group_iterator(self, S, m, proof=True):
        r"""
        Return an iterator through elements of the finite group `\QQ(S,m)`.

        INPUT:

        - ``S`` -- set of primes

        - ``m`` -- positive integer

        - ``proof`` -- ignored

        OUTPUT:

        An iterator yielding the distinct elements of `\QQ(S,m)`.  See
        the docstring for :meth:`selmer_generators` for more information.

        EXAMPLES::

            sage: list(QQ.selmer_group_iterator((), 2))
            [1, -1]
            sage: list(QQ.selmer_group_iterator((2,), 2))
            [1, 2, -1, -2]
            sage: list(QQ.selmer_group_iterator((2,3), 2))
            [1, 3, 2, 6, -1, -3, -2, -6]
            sage: list(QQ.selmer_group_iterator((5,), 2))
            [1, 5, -1, -5]
        """
        KSgens, ords = self.selmer_generators(S=S, m=m, proof=proof, orders=True)
        one = self.one()
        from itertools import product

        from sage.misc.misc_c import prod
        for ev in product(*[range(o) for o in ords]):
            yield prod((p**e for p, e in zip(KSgens, ev)), one)

    def selmer_space(self, S, p, proof=None):
        r"""
        Compute the group `\QQ(S,p)` as a vector space with maps to and from `\QQ^*`.

        INPUT:

        - ``S`` -- list of prime numbers

        - ``p`` -- a prime number

        OUTPUT:

        (tuple) ``QSp``, ``QSp_gens``, ``from_QSp``, ``to_QSp`` where

        - ``QSp`` is an abstract vector space over `\GF{p}` isomorphic to `\QQ(S,p)`;

        - ``QSp_gens`` is a list of elements of `\QQ^*` generating `\QQ(S,p)`;

        - ``from_QSp`` is a function from ``QSp`` to `\QQ^*`
          implementing the isomorphism from the abstract `\QQ(S,p)` to
          `\QQ(S,p)` as a subgroup of `\QQ^*/(\QQ^*)^p`;

        - ``to_QSP`` is a partial function from `\QQ^*` to ``QSp``,
          defined on elements `a` whose image in `\QQ^*/(\QQ^*)^p` lies in
          `\QQ(S,p)`, mapping them via the inverse isomorphism to the
          abstract vector space ``QSp``.

        The group `\QQ(S,p)` is the finite subgroup of
        `\QQ^*/(\QQ^*)^p` consisting of elements whose valuation at
        all primes not in `S` is a multiple of `p`.  It contains the
        subgroup of those `a\in \QQ^*` such that
        `\QQ(\sqrt[p]{a})/\QQ` is unramified at all primes of `\QQ`
        outside of `S`, but may contain it properly when `p` is not in `S`.

        EXAMPLES:

        When `S` is empty, `\QQ(S,p)` is only nontrivial for `p=2`::

            sage: QS2, QS2gens, fromQS2, toQS2 = QQ.selmer_space([], 2)                 # needs sage.rings.number_field
            sage: QS2                                                                   # needs sage.rings.number_field
            Vector space of dimension 1 over Finite Field of size 2
            sage: QS2gens                                                               # needs sage.rings.number_field
            [-1]

            sage: all(QQ.selmer_space([], p)[0].dimension() == 0                        # needs sage.libs.pari sage.rings.number_field
            ....:     for p in primes(3, 10))
            True

        In general there is one generator for each `p\in S`, and an
        additional generator of `-1` when `p=2`::

            sage: # needs sage.modules sage.rings.number_field
            sage: QS2, QS2gens, fromQS2, toQS2 = QQ.selmer_space([5,7], 2)
            sage: QS2
            Vector space of dimension 3 over Finite Field of size 2
            sage: QS2gens
            [5, 7, -1]
            sage: toQS2(-7)
            (0, 1, 1)
            sage: fromQS2((0,1,1))
            -7

        The map ``fromQS2`` is only well-defined modulo `p`-th powers
        (in this case, modulo squares)::

            sage: toQS2(-5/7)                                                           # needs sage.modules sage.rings.number_field
            (1, 1, 1)
            sage: fromQS2((1,1,1))                                                      # needs sage.modules sage.rings.number_field
            -35
            sage: ((-5/7)/(-35)).is_square()
            True

        The map ``toQS2`` is not defined on all of `\QQ^*`, only on
        those numbers which are squares away from `5` and `7`::

            sage: toQS2(210)                                                            # needs sage.modules sage.rings.number_field
            Traceback (most recent call last):
            ...
            ValueError: argument 210 should have valuations divisible by 2
            at all primes in [5, 7]
        """
        from sage.rings.number_field.selmer_group import pSelmerGroup
        return pSelmerGroup(self, S, p)

    def quadratic_defect(self, a, p, check=True):
        r"""
        Return the valuation of the quadratic defect of `a` at `p`.

        INPUT:

        - ``a`` -- an element of ``self``
        - ``p`` -- a prime ideal or a prime number
        - ``check`` -- (default: ``True``) check if `p` is prime

        REFERENCE:

        [Kir2016]_

        EXAMPLES::

            sage: QQ.quadratic_defect(0, 7)
            +Infinity
            sage: QQ.quadratic_defect(5, 7)
            0
            sage: QQ.quadratic_defect(5, 2)
            2
            sage: QQ.quadratic_defect(5, 5)
            1
        """
        from sage.arith.misc import legendre_symbol
        from sage.rings.infinity import Infinity
        if a not in self:
            raise TypeError(str(a) + " must be an element of " + str(self))
        if p.parent() == ZZ.ideal_monoid():
            p = p.gen()
        if check and not p.is_prime():
            raise ValueError(str(p) + " must be prime")
        if a.is_zero():
            return Infinity
        v, u = self(a).val_unit(p)
        if v % 2 == 1:
            return v
        if p != 2:
            if legendre_symbol(u, p) == 1:
                return Infinity
            else:
                return v
        if p == 2:
            if u % 8 == 1:
                return Infinity
            if u % 8 == 5:
                return v + 2
            if u % 8 in [3, 7]:
                return v + 1

    #################################
    #  Coercions to interfaces
    #################################
    def _gap_init_(self):
        r"""
        Return the GAP representation of `\QQ`.

        EXAMPLES::

            sage: gap(QQ)  # indirect doctest                                           # needs sage.libs.gap
            Rationals
        """
        return 'Rationals'

    def _magma_init_(self, magma):
        r"""
        Return the magma representation of `\QQ`.

        EXAMPLES::

            sage: magma(QQ)      # optional - magma # indirect doctest
            Rational Field

        TESTS:

        See :issue:`5521`::

            sage: loads(dumps(QQ)) == QQ  # optional - magma
            True
        """
        return 'RationalField()'

    def _macaulay2_init_(self, macaulay2=None):
        r"""
        Return the macaulay2 representation of `\QQ`.

        EXAMPLES::

            sage: macaulay2(QQ)   # optional - macaulay2 # indirect doctest
            QQ
        """
        return "QQ"

    def _axiom_init_(self):
        r"""
        Return the axiom/fricas representation of `\QQ`.

        EXAMPLES::

           sage: axiom(QQ)    #optional - axiom # indirect doctest
           Fraction Integer
           sage: fricas(QQ)   #optional - fricas # indirect doctest
           Fraction(Integer)
        """
        return 'Fraction Integer'

    _fricas_init_ = _axiom_init_

    def _polymake_init_(self):
        r"""
        Return the polymake representation of `\QQ`.

        EXAMPLES::

            sage: polymake(QQ)    #optional - jupymake # indirect doctest
            Rational
        """
        return '"Rational"'

    def _sympy_(self):
        r"""
        Return the SymPy set ``Rationals``.

        EXAMPLES::

            sage: QQ._sympy_()                                                          # needs sympy
            Rationals
        """
        from sage.interfaces.sympy import sympy_init
        from sympy import Rationals
        sympy_init()
        return Rationals

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(QQ, verify=True)
            # Verified
            QQ
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: QQ._sage_input_(SageInputBuilder(), False)
            {atomic:QQ}
        """
        return sib.name('QQ')

    def _factor_univariate_polynomial(self, f):
        """
        Factor the univariate polynomial ``f``.

        INPUT:

        - ``f`` -- a univariate polynomial defined over the rationals

        OUTPUT:

        - A factorization of ``f`` over the rationals into a unit and monic
          irreducible factors

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.factor`.

            This method calls PARI to compute the factorization.

        TESTS::

            sage: # needs sage.libs.pari
            sage: R.<x> = QQ[]
            sage: QQ._factor_univariate_polynomial(x)
            x
            sage: QQ._factor_univariate_polynomial(2*x)
            (2) * x
            sage: QQ._factor_univariate_polynomial((x^2 - 1/4)^4)
            (x - 1/2)^4 * (x + 1/2)^4
            sage: QQ._factor_univariate_polynomial((2*x + 1) * (3*x^2 - 5)^2)
            (18) * (x + 1/2) * (x^2 - 5/3)^2
            sage: f = prod((k^2*x^k + k)^(k-1) for k in primes(10))
            sage: QQ._factor_univariate_polynomial(f)
            (1751787911376562500) * (x^2 + 1/2) * (x^3 + 1/3)^2 * (x^5 + 1/5)^4 * (x^7 + 1/7)^6
            sage: QQ._factor_univariate_polynomial(10*x^5 - 1)
            (10) * (x^5 - 1/10)
            sage: QQ._factor_univariate_polynomial(10*x^5 - 10)
            (10) * (x - 1) * (x^4 + x^3 + x^2 + x + 1)
        """
        G = f._pari_with_name().factor()

        # normalize the leading coefficients
        P = f.parent()
        F = [(P(g).monic(), int(e)) for g, e in zip(*G)]

        from sage.structure.factorization import Factorization
        return Factorization(F, f.leading_coefficient())

    def valuation(self, p):
        r"""
        Return the discrete valuation with uniformizer ``p``.

        EXAMPLES::

            sage: v = QQ.valuation(3); v                                                # needs sage.rings.padics
            3-adic valuation
            sage: v(1/3)                                                                # needs sage.rings.padics
            -1

        .. SEEALSO::

            :meth:`NumberField_generic.valuation() <sage.rings.number_field.number_field.NumberField_generic.valuation>`,
            :meth:`IntegerRing_class.valuation() <sage.rings.integer_ring.IntegerRing_class.valuation>`
        """
        from sage.rings.padics.padic_valuation import pAdicValuation
        return pAdicValuation(self, p)


QQ = RationalField()
Q = QQ


def is_RationalField(x) -> bool:
    """
    Check to see if ``x`` is the rational field.

    EXAMPLES::

        sage: from sage.rings.rational_field import is_RationalField as is_RF
        sage: is_RF(QQ)
        doctest:warning...
        DeprecationWarning: The function is_RationalField is deprecated;
        use 'isinstance(..., RationalField)' instead.
        See https://github.com/sagemath/sage/issues/38128 for details.
        True
        sage: is_RF(ZZ)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38128,
                "The function is_RationalField is deprecated; "
                "use 'isinstance(..., RationalField)' instead.")
    return isinstance(x, RationalField)


def frac(n, d):
    """
    Return the fraction ``n/d``.

    EXAMPLES::

        sage: from sage.rings.rational_field import frac
        sage: frac(1,2)
        1/2
    """
    return Rational(n) / Rational(d)
