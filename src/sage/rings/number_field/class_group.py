r"""
Class groups of number fields

Sage can compute class groups, ray class groups, and `S`-class groups of number
fields, and does so by wrapping the functionality from the PARI C-library. Some
of what can be computed includes the group structure, representative ideals, the
class of a given ideal, generators of the group, and products of elements. This
file also implements moduli of number fields.

AUTHORS:

- ??? (???): Initial version
- Robert Harron (2016-08-15): implemented ray class groups and moduli

EXAMPLES:

Computations with the class group of a quadratic field::

    sage: x = polygen(ZZ, 'x')
    sage: K.<a> = NumberField(x^2 + 23)
    sage: H = K.class_group(); H
    Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23
    sage: I = H.gen(); I
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

Computations with an `S`-class group of a quadratic field::

    sage: K.<a> = QuadraticField(40)
    sage: CS = K.S_class_group(K.primes_above(31)); CS
    S-class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 - 40
    sage: CS.gens()   # random gens (platform dependent)
    (Fractional S-ideal class (3, 1/2*a + 2),)
    sage: CS(2)
    Trivial S-ideal class
    sage: c1 = CS(K.ideal([6, a+2])); c1
    Fractional S-ideal class (6, a + 2)
    sage: c2 = CS(K.ideal([6, a+4])); c2
    Fractional S-ideal class (6, a + 4)
    sage: c1 == c2
    True
    sage: c1.order()
    2

Computations with a ray class group of a quadratic field::

    sage: F = QuadraticField(40)
    sage: m = F.ideal(3).modulus([0, 1]); m
    (Fractional ideal (3)) * infinity_0 * infinity_1
    sage: R = F.ray_class_group(m); R
    Ray class group of order 8 with structure C4 x C2 of Number Field in a with defining polynomial x^2 - 40 of modulus (Fractional ideal (3)) * infinity_0 * infinity_1

Unlike for class groups and `S`-class groups, ray class group elements
do not carry around a representative ideal (for reasons of efficiency).
Nevertheless, one can be demanded. The returned ideal should be somewhat
'small.'

::

    sage: R.gens()
    (c0, c1)
    sage: R.gens_ideals()
    (Fractional ideal (430, 1/2*a + 200), Fractional ideal (-3/2*a + 2))
    sage: c = R.gens()[0]^3 * R.gens()[1]; c
    c0^3*c1
    sage: c.ideal()
    Fractional ideal (10, a)
    sage: c = R(F.ideal(2)); c
    c0^2
    sage: R.gen(0).ideal()^2
    Fractional ideal (30*a - 470)
    sage: R(R.gen(0).ideal()^2).ideal()
    Fractional ideal (2)

Narrow class groups are implemented via ray class groups::

    sage: F.<a> = QuadraticField(3)
    sage: F.class_group()
    Class group of order 1 of Number Field in a with defining polynomial x^2 - 3
    sage: Hn = F.narrow_class_group(); Hn
    Narrow class group of order 2 with structure C2 of Number Field in a with defining polynomial x^2 - 3
    sage: Hn.gens()
    (c,)
    sage: Hn.gens_ideals()
    (Fractional ideal (a + 1),)
"""

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.groups.abelian_gps.values import AbelianGroupWithValues_class, AbelianGroupWithValuesElement
from sage.groups.abelian_gps.abelian_group import AbelianGroup_class
from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement
from sage.structure.element import MonoidElement
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.libs.pari.all import pari


def _integer_n_tuple_L1_iterator(n):
    if n == 1:
        i = 1
        while True:
            yield i
            yield -i
            i += 1
    else:
        from sage.combinat.partition import OrderedPartitions
        from sage.combinat.subset import Subsets
        N = 1
        sign_options_dict = {}
        Subsets_of_n = []
        while True:
            for k in range(1, n + 1):
                Ps = OrderedPartitions(N, k)
                for P in Ps:
                    try:
                        Ss = Subsets_of_n[k - 1]
                    except IndexError:
                        Ss = Subsets(range(n), k)
                        Subsets_of_n.append(Ss)
                    for S in Ss:
                        i = [0] * n
                        for j in range(k):
                            i[S[j]] = P[j]
                        yield i
                        try:
                            sign_options = sign_options_dict[S]
                        except KeyError:
                            sign_options = Subsets(S)[1:]
                            sign_options_dict[S] = sign_options
                        for signs in sign_options:
                            ii = list(i)
                            for index in signs:
                                ii[index] = - ii[index]
                            yield ii
            N += 1

class Modulus(SageObject):
    def __init__(self, finite, infinite=None, check=True):
        r"""
        Create a modulus of a number field.

        INPUT:

        - ``finite`` -- a non-zero fractional ideal in a number field.
        - ``infinite`` -- a list of indices corresponding to real places
          of the number field, sorted.
        - ``check`` (default: True) -- If ``True``, run a few checks on the input.
        """
        self._finite = finite
        if infinite is None:
            self._infinite = ()
        else:
            self._infinite = tuple(ZZ(i) for i in infinite)
        K = self._finite.number_field()
        self._number_field = K
        if check:
            #insert various checks here
            if self._finite == 0:
                raise ValueError("Finite component of a modulus must be non-zero.")
            sgn = K.signature()[0]
            for i in self._infinite:
                if i < 0 or i >= sgn:
                    raise ValueError("Infinite component of a modulus must be a list non-negative integers less than the number of real places of K")
        return

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: K.<a> = NumberField(x^2-5)
            sage: m = K.modulus(K.ideal(31), [0,1]); m
            (Fractional ideal (31)) * infinity_0 * infinity_1
        """
        if len(self._infinite) == 0:
            return str(self._finite)
        str_inf = ''
        for i in self._infinite:
            str_inf += ' * infinity_%s'%(i)
        return '(' + str(self._finite) + ')' + str_inf

    def __eq__(self, other):
        return self._number_field == other._number_field and self._finite == other._finite and self._infinite == other._infinite

    def __mul__(self, other):
        r"""
        Multiply two moduli.

        This multiplies the two finite parts and performs an exclusive or on the real places.

        EXAMPLES::

            sage: K = NumberField(x^3 - 2, 'a')
            sage: m1 = K.modulus(K.ideal(2), [0])
            sage: m2 = K.modulus(K.ideal(3), [0])
            sage: m1 * m2
            Fractional ideal (6)

        A higher degree totally real field::

            sage: K = NumberField(x^5 - x^4 - 4*x^3 + 3*x^2 + 3*x - 1, 'a')
            sage: m1 = K.modulus(K.ideal(5), [2, 3])
            sage: m2 = K.modulus(K.ideal(25), [0, 1, 3, 4])
            sage: m1 * m2
            (Fractional ideal (125)) * infinity_0 * infinity_1 * infinity_2 * infinity_4
            sage: _ == m2 * m1
            True
        """
        inf = tuple(set(self.infinite_part()).symmetric_difference(other.infinite_part()))
        return Modulus(self.finite_part() * other.finite_part(), inf, check=False)

    def lcm(self, other):
        inf = tuple(set(self.infinite_part()).union(other.infinite_part()))
        #Pe_out = []
        self_fact_P, self_fact_e = zip(*self.finite_part().factor())
        self_fact_P = list(self_fact_P)
        self_fact_e = list(self_fact_e)
        #self_facts = self.finite_part().factor()
        #of = other.finite_part()
        other_facts = other.finite_part().factor()
        mf = self._number_field.ideal_monoid().one()
        for P, e in other_facts:
            try:
                i = self_fact_P.index(P)
            except ValueError:
                #Pe_out.append([P, e])
                mf *= P**e
                continue
            #Pe_out.append([P, max(e, self_fact_e[i])])
            mf *= P**max(e, self_fact_e[i])
            del self_fact_P[i]
            del self_fact_e[i]
        for i in range(len(self_fact_P)):
            mf *= self_fact_P[i]**self_fact_e[i]
        return Modulus(mf, inf, check=False)

    def divides(self, other):
        if not set(self.infinite_part()).issubset(other.infinite_part()):
            return False
        return self.finite_part().divides(other.finite_part())

    def number_field(self):
        return self._number_field

    def finite_part(self):
        return self._finite

    def infinite_part(self):
        return self._infinite

    def finite_factors(self):
        try:
            return self._finite_factors
        except AttributeError:
            self._finite_factors = self.finite_part().factor()
            return self._finite_factors

    #def _pari_finite_factors(self):
    #    """
    #    Return
    #    """
    #    return self._number_field.pari_nf().idealfactor(self._finite)

    def equivalent_coprime_ideal_multiplier(self, I, other):
        r"""
        Given ``I`` coprime to this modulus `m`, return a number field element `\beta`
        such that `\beta I` is coprime to the modulus ``other`` and equivalent to
        ``I`` `\mathrm{mod}^\ast m`; in particular, `\beta` will be `1 \mathrm{mod}^\ast m`.

        EXAMPLES:

        An example with two prime factors difference between this modulus and ``other``.

        ::

            sage: F.<a> = QuadraticField(5)
            sage: m_small = F.modulus(3/2*a - 1/2, [0, 1])
            sage: m_big = F.modulus(2*a - 30, [0, 1])
            sage: m_small.equivalent_coprime_ideal_multiplier(F.ideal(6), m_big)
            109/54
        """
        F = self._number_field
        other_Ps = [P for P, _ in other.finite_factors() if self._finite.valuation(P) == 0]
        if len(other_Ps) == 0: #If prime factors of other is a subset of prime factors of this modulus
            return F.one()
        alpha = I.idealcoprime(other._finite)
        if self.number_is_one_mod_star(alpha):
            return alpha
        nf = F.pari_nf()
        other_Ps_facts = [nf.idealfactor(P) for P in other_Ps]
        other_Ps_fact_mat = pari(other_Ps_facts).Col().matconcat()
        self_fact_mat = nf.idealfactor(self.finite_part())
        x = pari([self_fact_mat, other_Ps_fact_mat]).Col().matconcat()
        conditions = [~alpha] * self_fact_mat.nrows() + [1] * len(other_Ps)
        gamma = F(nf.idealchinese(x, conditions))
        beta = alpha * gamma
        if self.number_is_one_mod_star(beta):
            return beta
        from sage.misc.misc_c import prod
        beta_fixed = F.modulus(self._finite * prod(other_Ps), self._infinite).fix_signs(beta) #should be able to do this more efficiently
        return beta_fixed

    def equivalent_ideal_coprime_to_other(self, I, other):
        r"""
        Given ``I`` coprime to this modulus `m`, return an ideal `J` such that `J` is coprime
        to the modulus ``other`` and equivalent to ``I`` `\mathrm{mod}^\ast m`.

        This is useful for lowering the level of a non-primitive Hecke character.

        INPUT:

        - ``I`` -- an ideal relatively prime to this modulus (not checked).
        - ``other`` -- some other modulus.

        OUTPUT:

        an ideal coprime to ``other`` and equivalent to ``I`` in the ray class
        group modulo this modulus.
        """
        return self.equivalent_coprime_ideal_multiplier(I, other) * I

    def number_is_one_mod_star(self, a):
        K = self.number_field()
        am1 = K(a - 1)
        for P, e in self.finite_factors():
            if am1.valuation(P) < e:
                return False
        inf_places = K.places()
        for i in self.infinite_part():
            if inf_places[i](a) <= 0:
                return False
        return True

    def fix_signs(self, a):
        r"""
        Given ``a`` in ``self.number_field()``, find `b` congruent to ``a`` `mod^\ast` ``self.finite_part()``
        such that `b` is positive at the infinite places dividing ``self``.
        """
        if self.is_finite() or a == 0:
            return a
        places = self.number_field().places()
        positive = []
        negative = []
        for i in self.infinite_part():
            if places[i](a) > 0:
                positive.append(i)
            else:
                negative.append(i)
        if not negative:
            return a
        t = self.get_one_mod_star_finite_with_fixed_signs(positive, negative)
        return t * a

    def get_one_mod_star_finite_with_fixed_signs(self, positive, negative):
        if len(negative) == 0:
            return self.number_field().one()
        negative = tuple(negative)
        try:
            return self._one_mod_star[negative]
        except AttributeError:
            self._one_mod_star = {}
        except KeyError:
            pass
        try:
            beta_is, Ainv = self._beta_is_Ainv
        except AttributeError:
            beta_is, Ainv = self._find_beta_is_Ainv()
        d = len(self.infinite_part())
        v = [0] * d
        for i in negative:
            v[i] = 1
        v = (GF(2)**d)(v)
        w = Ainv * v
        t = self.number_field().one()
        for i in range(d):
            if w[i] != 0:
                t *= (1 + beta_is[i])
        self._one_mod_star[negative] = t
        return t

    def _find_beta_is_Ainv(self):
        r"""
        Step 2 of Algorithm 4.2.20 of Cohen's Advanced...
        """
        from sage.matrix.special import column_matrix
        gammas = self.finite_part().basis()
        k = len(self.infinite_part())
        beta_is = []
        Acols = []
        V = GF(2)**k
        it = _integer_n_tuple_L1_iterator(k)
        while len(beta_is) < k:
            e = next(it)
            beta = sum(e[i] * gammas[i] for i in range(k))
            sbeta = V(self._signs(beta))
            Acols_new = Acols + [sbeta]
            A = column_matrix(GF(2), Acols_new)
            if A.rank() == len(Acols_new):
                Acols = Acols_new
                beta_is.append(beta)
        self._beta_is_Ainv = (beta_is, ~A)
        return self._beta_is_Ainv

    def _signs(self, b):
        if b == 0:
            raise ValueError("Non-zero input required.")
        sigmas = self._number_field.real_places()
        return [ZZ.one() if sigmas[i](b).sign() == -1 else ZZ.zero() for i in self.infinite_part()]

    def is_finite(self):
        return len(self._infinite) == 0

    def is_infinite(self):
        return self._finite.is_one()

    def __pari__(self):
        """
        Return the corresponding pari modulus.

        Note that this function performs the conversion between the ordering
        of the real places of the number field in Sage and the ordering of the
        underlying pari nf object.

        EXAMPLES:

        An example where the places in Sage and pari are in a different order.

        ::

            sage: x = polygen(QQ)
            sage: f = x^4 - x^3 - 3*x^2 + x + 1
            sage: F.<a> = NumberField(f(1-2*x))
            sage: F.modulus(1, [2, 3]).__pari__()[1]
            [0, 1, 1, 0]
        """
        inf_mod = [0] * self._number_field.signature()[0]
        conversion = self._number_field._pari_real_places_to_sage()
        for i in self._infinite:
            pari_index = conversion.index(i)
            inf_mod[pari_index] = 1
        return pari([self._finite, inf_mod])

    def __hash__(self):
        return hash((self._finite, self._infinite))

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
        sage: P2a, P2b = [P for P,e in (2*OK).factor()]
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
            sage: P2a, P2b = [P for P,e in (2*OK).factor()]
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
        if self.ideal().is_prime():
            return self.ideal()
        c = self.reduce()
        if c.ideal().is_prime():
            return c.ideal()
        # otherwise we just search:
        Cl = self.parent()
        K = Cl.number_field()
        from sage.rings.real_mpfr import RR
        for P in K.primes_of_bounded_norm_iter(RR(norm_bound)):
            if Cl(P) == c:
                return P
        raise RuntimeError("No prime of norm less than %s found in class %s" % (norm_bound, c))

    def gens(self) -> tuple:
        r"""
        Return generators for a representative ideal in this
        (`S`-)ideal class.

        EXAMPLES::

            sage: K.<w> = QuadraticField(-23)
            sage: OK = K.ring_of_integers()
            sage: C = OK.class_group()
            sage: P2a, P2b = [P for P,e in (2*OK).factor()]
            sage: c = C(P2a); c
            Fractional ideal class (2, 1/2*w - 1/2)
            sage: c.gens()
            (2, 1/2*w - 1/2)
        """
        return self.ideal().gens()


class RayClassGroupElement(AbelianGroupElement):
    #@@def __init__(self, parent, element, ideal=None):
        #@@if element is None:
        #@@    if not parent.modulus().finite_part().is_coprime(ideal):
        #@@       raise ValueError("Ideal is not coprime to the modulus.")
        #@@    element = parent._ideal_log(ideal)
        #Should treat the else case for coprime-ness as well since the code can coerce from different moduli
        #@@FractionalIdealClass.__init__(self, parent, element, ideal)

    #def _repr_(self):
    #    return AbelianGroupWithValuesElement._repr_(self)

    #Should be able to get rid of the operations if make the reduce function a method of the parent
    #def _mul_(self, other):
    #    m = AbelianGroupElement._mul_(self, other)
    #    nf = self.parent()._number_field.pari_nf()
    #    m._value = nf.idealred(nf.idealmul(self.value(), other.value()))
    #    return m

    #def _div_(self, other):
    #    m = AbelianGroupElement._div_(self, other)
    #    nf = self.parent()._number_field.pari_nf()
    #    m._value = nf.idealred(nf.idealdiv(self.value(), other.value()))
    #    return m

    #def inverse(self):
    #    m = AbelianGroupElement.inverse(self)
    #    nf = self.parent()._number_field.pari_nf()
    #    m._value = nf.idealred(nf.idealinv(self.value()))
    #    return m

    #__invert__ = inverse

    def ideal(self, reduce=True):
        """
        Return an ideal representing this ray class.

        If ``reduce`` is ``True`` (by default) the returned ideal is
        reduced to 'small' (this can be slow on large inputs).

        INPUT:

        - ``reduce`` -- (default: ``True``) determine whether or not
          to output a 'small' representative.

        OUTPUT:

        An ideal representing this ray class. If ``reduce`` is ``True``,
        the ideal returned is made 'small' by the ideal's
        ``reduce_equiv`` function (and ``reduce_equiv`` is used at
        each step of computing this representative). Otherwise, the
        output is just the appropriate product of the powers of the
        generators of the ray class group.

        EXAMPLES:

        Over a real quadratic field field of class number 1::

            sage: F.<a> = NumberField(x^2 - 5)
            sage: m = F.ideal(11).modulus([0, 1])
            sage: R = F.ray_class_group(m)
            sage: c0, c1 = R.gens()
            sage: c = c0^4*c1; c.ideal()
            Fractional ideal (-a)
            sage: c.ideal(False)
            Fractional ideal (-6242265*a + 1268055)

        Over a real quadratic field of class number 2::

            sage: F = QuadraticField(40)
            sage: R = F.ray_class_group(F.prime_above(13).modulus([0, 1]))
            sage: for c in R:
            ....:     if R(c.ideal()) != c:
            ....:         print("Bug!")
        """
        #R = self.parent()
        #exps = self.exponents()
        #gens = R.gens_values()
        #i = 0
        #while exps[i] == 0:
        #    i += 1
        #    if i == L:
        #        return R.one()
        #nf = R.number_field().pari_nf()
        #I = nf.idealpow(gens[i], exps[i], flag=1)
        #i += 1
        #while i < L:
        #    e = exps[i]
        #    g = gens[i]
        #    i += 1
        #    if e != 0:
        #        I = nf.idealmul(I, nf.idealpow(g, e, flag=1), flag=1)
        #bnr = R.pari_bnr()
        #return R.number_field().ideal(nf.idealmul(I[0], nf.nffactorback(I[1])))
        R = self.parent()
        #m = R.modulus()
        exps = self.exponents()
        gens = R.gens_ideals()
        L = len(exps)
        #Speed this up later using binary powering
        i = 0
        while exps[i] == 0:
            i += 1
            if i == L:
                return R.one()
        I = (gens[i]**exps[i])
        if reduce:
            #I = I.reduce_equiv(m)
            I = R.ideal_reduce(I)
        i += 1
        while i < L:
            e = exps[i]
            g = gens[i]
            i += 1
            if e != 0:
                I = (I * (g**e))
                if reduce:
                    #I = I.reduce_equiv(m)
                    I = R.ideal_reduce(I)
        return I

    #def reduce(self):
    #    nf = self.parent()._number_field.pari_nf()
    #    return RayClassGroupElement(self.parent(), self.exponents(), nf.idealred(self.value()))


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

class RayClassGroup(AbelianGroup_class):
    Element = RayClassGroupElement

    def __init__(self, gens_orders, names, modulus, gens, bnr, proof=True):
        r"""
        ``gens`` -- a tuple of pari extended ideals
        """
        #AbelianGroupWithValues_class.__init__(self, gens_orders, names, gens,
        #                                      values_group=modulus.number_field().ideal_monoid())
        AbelianGroup_class.__init__(self, gens_orders, names)
        self._gens = gens
        self._proof_flag = proof #TODO: use this, in _ideal_log?
        self._modulus = modulus
        self._number_field = modulus.number_field()
        self._bnr = bnr
        self._is_narrow = (len(modulus.infinite_part()) == self._number_field.signature()[0]) and modulus.finite_part().is_one()

    def _element_constructor_(self, *args, **kwds):
        try:
            L = len(args[0])
        except TypeError:
            L = -1
        if L == self.ngens():
            return self.element_class(self, args[0])
        if isinstance(args[0], RayClassGroupElement):
            c = args[0]
            if c.parent() is self:
                return self.element_class(self, c.exponents())
            nf = self._number_field.pari_nf()
            bnr = self._bnr
            old_exps = c.exponents()
            old_gens = c.parent().gens_values()
            L = len(old_exps)
            if L ==0:
                return self.one()
            i = 0
            while old_exps[i] == 0:
                i += 1
                if i == L:
                    return self.one()
            I = nf.idealpow(old_gens[i], old_exps[i], flag=1)
            i += 1
            while i < L:
                e = old_exps[i]
                g = old_gens[i]
                i += 1
                if e != 0:
                    I = nf.idealmul(I, nf.idealpow(g, e, flag=1), flag=1)
            exps = tuple(ZZ(c) for c in bnr.bnrisprincipal(nf.idealmul(I[0], nf.nffactorback(I[1])), flag = 0))
            return self.element_class(self, exps)
        else:
            I = self._number_field.ideal(*args, **kwds)
            if I.is_zero():
                raise TypeError("The zero ideal is not a fractional ideal")
            exps = self._ideal_log(I)
            return self.element_class(self, exps)

    def _repr_(self):
        if self._is_narrow:
            s0 = 'Narrow '
        else:
            s0 = 'Ray '
        s = 'class group of order %s '%self.order()
        if self.order() > 1:
            s += 'with structure %s '%self._group_notation(self.gens_orders())
        s += 'of %s'%(self._number_field)
        if self._is_narrow:
            return s0 + s
        s += ' of modulus %s'%(self._modulus)
        return s0 + s

    def ray_class_field(self, subgroup=None, names=None, algorithm='stark'):
        r"""
        Two different algorithms are possible: pari's :pari:`bnrstark` and
        :pari:`rnfkummer`. The first one uses the Stark conjecture and only
        deals with totally real extensions of a totally real base
        field. The second one uses Kummer theory and only deals with
        extensions of prime degree.

        INPUT:

        - algorithm -- (default: ``stark``) if the value is ``stark``,
          then pari's :pari:`bnrstark` function is tried first, and if that
          fails, :pari:`rnfkummer` will be attempted. If the value is
          ``kummer``, then pari's :pari:`rnfkummer` is tried first, with
          :pari:`bnrstark` as a backup. Using ``stark_only`` or ``kummer_only``
          will just raise an exception if the first attempt fails.

        OUTPUT:

        The class field corresponding to the given subgroup, or the
        ray class field if ``subgroup`` is ``None``, as a relative
        number field.

        EXAMPLES:

        Class fields of `\QQ(\sqrt{3})`::

            sage: F.<a> = QuadraticField(3)
            sage: m = F.ideal(7).modulus()
            sage: R = F.ray_class_group(m)
            sage: R.ray_class_field(names='b')
            Number Field in b with defining polynomial x^6 + a*x^5 - 4*x^4 - 4*a*x^3 + 2*x^2 + 2*a*x - 1 over its base field
            sage: S = R.subgroup([R.gen()^2])
            sage: R.ray_class_field(S, names='b')
            Number Field in b with defining polynomial x^2 - a*x - 1 over its base field
            sage: m = F.modulus(20)
            sage: R = F.ray_class_group(m)
            sage: S = R.subgroup([R.gens()[0]^2, R.gens()[1]])
            sage: R.ray_class_field(S, names='b')
            Number Field in b with defining polynomial x^2 + (a - 1)*x + 2*a - 4 over its base field

        An example where :pari:`bnrstark` fails, but :pari:`rnfkummer` saves the day::

            sage: F.<a> = NumberField(x^8 - 12*x^6 + 36*x^4 - 36*x^2 + 9)
            sage: m = F.ideal(2).modulus()
            sage: R = F.ray_class_group(m)
            sage: set_verbose(1)
            sage: K = R.ray_class_field(names='b'); K
            verbose 1 (...: class_group.py, ray_class_field) bnrstark failed; trying rnfkummer.
            Number Field in b with defining polynomial x^2 + (1/3*a^6 - 10/3*a^4 + 5*a^2)*x + 1/3*a^6 - 1/3*a^5 - 11/3*a^4 + 3*a^3 + 8*a^2 - 4*a - 5 over its base field
            sage: set_verbose(0)
        """
        if subgroup is not None:
            try:
                test_subgrp = (subgroup.ambient_group() is self)
            except AttributeError:
                subgroup = self.subgroup(subgroup)
                test_subgrp = (subgroup.ambient_group() is self)
            if not test_subgrp:
                raise ValueError("subgroup does not define a subgroup of this ray class group.")
            gens_coords = [h.exponents() for h in subgroup.gens()]
            if len(gens_coords) == 0:
                subgroup = None
                subgroup_mat = None
            else:
                from sage.matrix.special import column_matrix
                subgroup_mat = column_matrix(gens_coords)
        else:
            subgroup_mat = None

        from cypari2.handle_error import PariError
        from sage.misc.verbose import verbose

        bnr = self._bnr
        if algorithm == 'stark_only':
            if len(self._modulus.infinite_part()) > 0 or not self._number_field.is_totally_real():
                raise NotImplementedError("Stark's conjecture algorithm only implemented for totally real extensions of a totally real base field.")
            f = bnr.bnrstark(subgroup=subgroup_mat)
        elif algorithm == 'kummer_only':
            if (subgroup is None and not self.order().is_prime()) or (subgroup is not None and not self.order().divide_knowing_divisible_by(subgroup.order()).is_prime()):
                raise NotImplementedError("Kummer theory algorithm only implemented extensions of prime degree.")
            f = bnr.rnfkummer(subgp=subgroup_mat)
        elif algorithm == 'stark':
            if len(self._modulus.infinite_part()) > 0 or not self._number_field.is_totally_real():
                if (subgroup is None and not self.order().is_prime()) or (subgroup is not None and not self.order().divide_knowing_divisible_by(subgroup.order()).is_prime()):
                    raise NotImplementedError("Ray class fields only implemented for totally real extensions of totally real base fields, or for extensions of prime degree.")
                f = bnr.rnfkummer(subgp=subgroup_mat)
            else:
                try:
                    f = bnr.bnrstark(subgroup=subgroup_mat)
                except PariError:
                    if (subgroup is None and self.order().is_prime()) or (subgroup is not None and self.order().divide_knowing_divisible_by(subgroup.order()).is_prime()):
                        verbose("bnrstark failed; trying rnfkummer.")
                        f = bnr.rnfkummer(subgp=subgroup_mat)
                    else:
                        raise
        elif algorithm == 'kummer':
            if (subgroup is None and self.order().is_prime()) or (subgroup is not None and self.order().divide_knowing_divisible_by(subgroup.order()).is_prime()):
                f = bnr.rnfkummer(subgp=subgroup_mat)
            else:
                f = bnr.bnrstark(subgroup=subgroup_mat)
        else:
            raise ValueError("Value of algorithm must be one of \'stark\', \'stark_only\', \'kummer\', or \'kummer_only\'.")
        if f.type() == 't_VEC':
            raise NotImplementedError("bnrstark returned a list of polynomials. Dealing with this has not been implemented.")
        F = self._number_field
        nf = F.pari_nf()
        f = nf.rnfpolredbest(f)
        d = f.poldegree()
        cs = [F(f.polcoef(i)) for i in range(d + 1)]
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        f = PolynomialRing(F, 'x')(cs)
        return F.extension(f, names=names)

    def _ray_class_field_stark(self, subgroup=None, names=None):
        pass

    def _ray_class_field_kummer(self, subgroup=None, names=None):
        pass

    def _ideal_log(self, ideal):
        return tuple(ZZ(c) for c in self._bnr.bnrisprincipal(ideal, flag=0))

    def ideal_reduce(self, ideal):
        from cypari2.handle_error import PariError
        ideal = pari(ideal)
        try:
            pari_ideal = self._bnr.idealmoddivisor(ideal)
        except PariError as err:
            if err.errtext().find('not coprime') != -1:
                raise ValueError('Ideal in ideal_reduce is not coprime to the modulus of this ray class group.')
            else:
                raise err
        return self._number_field.ideal(pari_ideal)

    def gens_values(self):
        return self._gens

    @cached_method
    def gens_ideals(self):
        return tuple(self._number_field._pari_extended_ideal_to_sage(v) for v in self.gens_values())

    def modulus(self):
        return self._modulus

    def number_field(self):
        return self._number_field

    def pari_bnr(self):
        return self._bnr

    def pari_gens(self):
        return self._bnr[4][2]


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
