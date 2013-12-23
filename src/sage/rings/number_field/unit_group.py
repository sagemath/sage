r"""
Unit and S-unit groups of Number Fields

EXAMPLES::

    sage: x = polygen(QQ)
    sage: K.<a> = NumberField(x^4-8*x^2+36)
    sage: UK = UnitGroup(K); UK
    Unit group with structure C4 x Z of Number Field in a with defining polynomial x^4 - 8*x^2 + 36

The first generator is a primitive root of unity in the field::

    sage: UK.gens()
    (u0, u1)
    sage: UK.gens_values()  # random
    [-1/12*a^3 + 1/6*a, 1/24*a^3 + 1/4*a^2 - 1/12*a - 1]
    sage: UK.gen(0).value()
    -1/12*a^3 + 1/6*a

    sage: UK.gen(0)
    u0
    sage: UK.gen(0) + K.one()   # coerce abstract generator into number field
    -1/12*a^3 + 1/6*a + 1

    sage: [u.multiplicative_order() for u in UK.gens()]
    [4, +Infinity]
    sage: UK.rank()
    1
    sage: UK.ngens()
    2

Units in the field can be converted into elements of the unit group represented
as elements of an abstract multiplicative group::

    sage: UK(1)
    1
    sage: UK(-1)
    u0^2
    sage: [UK(u) for u in (x^4-1).roots(K,multiplicities=False)]
    [1, u0^2, u0^3, u0]

    sage: UK.fundamental_units() # random
    [1/24*a^3 + 1/4*a^2 - 1/12*a - 1]
    sage: torsion_gen = UK.torsion_generator();  torsion_gen
    u0
    sage: torsion_gen.value()
    -1/12*a^3 + 1/6*a
    sage: UK.zeta_order()
    4
    sage: UK.roots_of_unity()
    [-1/12*a^3 + 1/6*a, -1, 1/12*a^3 - 1/6*a, 1]

Exp and log functions provide maps between units as field elements and exponent
vectors with respect to the generators::

    sage: u = UK.exp([13,10]); u # random
    -41/8*a^3 - 55/4*a^2 + 41/4*a + 55
    sage: UK.log(u)
    (1, 10)
    sage: u = UK.fundamental_units()[0]
    sage: [UK.log(u^k) == (0,k) for k in range(10)]
    [True, True, True, True, True, True, True, True, True, True]
    sage: all([UK.log(u^k) == (0,k) for k in range(10)])
    True

    sage: K.<a> = NumberField(x^5-2,'a')
    sage: UK = UnitGroup(K)
    sage: UK.rank()
    2
    sage: UK.fundamental_units()
    [a^3 + a^2 - 1, a - 1]

S-unit groups may be constructed, where S is a set of primes::

    sage: K.<a> = NumberField(x^6+2)
    sage: S = K.ideal(3).prime_factors(); S
    [Fractional ideal (3, a + 1), Fractional ideal (3, a - 1)]
    sage: SUK = UnitGroup(K,S=tuple(S)); SUK
    S-unit group with structure C2 x Z x Z x Z x Z of Number Field in a with defining polynomial x^6 + 2 with S = (Fractional ideal (3, a + 1), Fractional ideal (3, a - 1))
    sage: SUK.primes()
    (Fractional ideal (3, a + 1), Fractional ideal (3, a - 1))
    sage: SUK.rank()
    4
    sage: SUK.gens_values()
    [-1, a^2 + 1, a^5 + a^4 - a^2 - a - 1, a + 1, -a + 1]
    sage: u = 9*prod(SUK.gens_values()); u
    -18*a^5 - 18*a^4 - 18*a^3 - 9*a^2 + 9*a + 27
    sage: SUK.log(u)
    (1, 3, 1, 7, 7)
    sage: u == SUK.exp((1,3,1,7,7))
    True

A relative number field example::

    sage: L.<a, b> = NumberField([x^2 + x + 1, x^4 + 1])
    sage: UL = L.unit_group(); UL
    Unit group with structure C24 x Z x Z x Z of Number Field in a with defining polynomial x^2 + x + 1 over its base field
    sage: UL.gens_values() # random
    [-b^3*a - b^3, -b^3*a + b, (-b^3 - b^2 - b)*a - b - 1, (-b^3 - 1)*a - b^2 + b - 1]
    sage: UL.zeta_order()
    24
    sage: UL.roots_of_unity()
    [b*a, -b^2*a - b^2, b^3, -a, b*a + b, -b^2, -b^3*a, -a - 1, b, b^2*a, -b^3*a - b^3, -1, -b*a, b^2*a + b^2, -b^3, a, -b*a - b, b^2, b^3*a, a + 1, -b, -b^2*a, b^3*a + b^3, 1]

A relative extension example, which worked thanks to the code review by F.W.Clarke::

    sage: PQ.<X> = QQ[]
    sage: F.<a, b> = NumberField([X^2 - 2, X^2 - 3])
    sage: PF.<Y> = F[]
    sage: K.<c> = F.extension(Y^2 - (1 + a)*(a + b)*a*b)
    sage: K.unit_group()
    Unit group with structure C2 x Z x Z x Z x Z x Z x Z x Z of Number Field in c with defining polynomial Y^2 + (-2*b - 3)*a - 2*b - 6 over its base field

TESTS::

    sage: UK == loads(dumps(UK))
    True
    sage: UL == loads(dumps(UL))
    True

AUTHOR:

- John Cremona
"""
#*****************************************************************************
#       Copyright (C) 2009 William Stein, John Cremona
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.abelian_gps.values import AbelianGroupWithValues_class
from sage.structure.sequence import Sequence
from sage.structure.proof.proof import get_flag
from sage.libs.pari.all import pari
from sage.misc.misc import prod
from sage.rings.integer_ring import ZZ

class UnitGroup(AbelianGroupWithValues_class):
    """
    The unit group or an S-unit group of a number field.

    TESTS::

        sage: x = polygen(QQ)
        sage: K.<a> = NumberField(x^4 + 23)
        sage: UK = K.unit_group()
        sage: u = UK.an_element();  u
        u0*u1
        sage: u.value()
        -1/4*a^3 + 7/4*a^2 - 17/4*a + 19/4

        sage: x = polygen(QQ)
        sage: K.<a> = NumberField(x^4 + 23)
        sage: K.unit_group().gens_values() # random
        [-1, 1/4*a^3 - 7/4*a^2 + 17/4*a - 19/4]

        sage: x = polygen(QQ)
        sage: U = NumberField(x^2 + x + 23899, 'a').unit_group(); U
        Unit group with structure C2 of Number Field in a with defining polynomial x^2 + x + 23899
        sage: U.ngens()
        1

        sage: K.<z> = CyclotomicField(13)
        sage: UK = K.unit_group()
        sage: UK.ngens()
        6
        sage: UK.gen(0) # random
        -z^11
        sage: UK.gen(1) # random
        z^5 + z^3
        sage: UK.gen(2) # random
        z^6 + z^5
        sage: UK.gen(3) # random
        z^9 + z^7 + z^5
        sage: UK.gen(4) # random
        z^9 + z^5 + z^4 + 1
        sage: UK.gen(5) # random
        z^5 + z

    An S-unit group::

        sage: SUK = UnitGroup(K,S=21); SUK
        S-unit group with structure C26 x Z x Z x Z x Z x Z x Z x Z x Z x Z x Z of Cyclotomic Field of order 13 and degree 12 with S = (Fractional ideal (3, z^3 - z - 1), Fractional ideal (3, z^3 + z^2 + z - 1), Fractional ideal (3, z^3 + z^2 - 1), Fractional ideal (3, z^3 - z^2 - z - 1), Fractional ideal (7))
        sage: SUK.rank()
        10
        sage: SUK.zeta_order()
        26
        sage: SUK.log(21*z)
        (6, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
    """
    # This structure is not a parent in the usual sense. The
    # "elements" are NumberFieldElement_absolute. Instead, they should
    # derive from AbelianGroupElement and coerce into
    # NumberFieldElement_absolute.

    def __init__(self, number_field, proof=True, S=None):
        """
        Create a unit group of a number field.

        INPUT:

        - ``number_field`` - a number field
        - ``proof`` - boolean (default True): proof flag
        - ``S`` - tuple of prime ideals, or an ideal, or a single
          ideal or element from which an ideal can be constructed, in
          which case the support is used.  If None, the global unit
          group is constructed; otherwise, the S-unit group is
          constructed.

        The proof flag is passed to pari via the ``pari_bnf()`` function
        which computes the unit group.  See the documentation for the
        number_field module.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2-38)
            sage: UK = K.unit_group(); UK
            Unit group with structure C2 x Z of Number Field in a with defining polynomial x^2 - 38
            sage: UK.gens()
            (u0, u1)
            sage: UK.gens_values()
            [-1, 6*a - 37]

            sage: K.<a> = QuadraticField(-3)
            sage: UK = K.unit_group(); UK
            Unit group with structure C6 of Number Field in a with defining polynomial x^2 + 3
            sage: UK.gens()
            (u,)
            sage: UK.gens_values()
            [-1/2*a + 1/2]

            sage: K.<z> = CyclotomicField(13)
            sage: UK = K.unit_group(); UK
            Unit group with structure C26 x Z x Z x Z x Z x Z of Cyclotomic Field of order 13 and degree 12
            sage: UK.gens()
            (u0, u1, u2, u3, u4, u5)
            sage: UK.gens_values() # random
            [-z^11, z^5 + z^3, z^6 + z^5, z^9 + z^7 + z^5, z^9 + z^5 + z^4 + 1, z^5 + z]
            sage: SUK = UnitGroup(K,S=2); SUK
            S-unit group with structure C26 x Z x Z x Z x Z x Z x Z of Cyclotomic Field of order 13 and degree 12 with S = (Fractional ideal (2),)

            """
        proof = get_flag(proof, "number_field")
        K = number_field
        pK = K.pari_bnf(proof)
        self.__number_field = K
        self.__pari_number_field = pK

        # process the parameter S:
        if not S:
            S = self.__S = ()
        else:
            if type(S)==list:
                S = tuple(S)
            if not type(S)==tuple:
                try:
                    S = tuple(K.ideal(S).prime_factors())
                except (NameError, TypeError, ValueError):
                    raise ValueError("Cannot make a set of primes from %s"%(S,))
            else:
                try:
                    S = tuple(K.ideal(P) for P in S)
                except (NameError, TypeError, ValueError):
                    raise ValueError("Cannot make a set of primes from %s"%(S,))
                if not all([P.is_prime() for P in S]):
                    raise ValueError("Not all elements of %s are prime ideals"%(S,))
            self.__S = S
            self.__pS = pS = [P.pari_prime() for P in S]

        # compute the fundamental units via pari:
        fu = [K(u) for u in pK.bnfunit()]
        self.__nfu = len(fu)

        # compute the additional S-unit generators:
        if S:
            self.__S_unit_data = pK.bnfsunit(pS)
            su = [K(u) for u in self.__S_unit_data[0]]
        else:
            su = []
        self.__nsu = len(su)
        self.__rank = self.__nfu + self.__nsu

        # compute a torsion generator and pick the 'simplest' one:
        n, z = pK.nfrootsof1()
        n = ZZ(n)
        self.__ntu = n
        z = K(z)

        # If we replaced z by another torsion generator we would need
        # to allow for this in the dlog function!  So we do not.

        # Store the actual generators (torsion first):
        gens = [z] + fu + su
        values = Sequence(gens, immutable=True, universe=self, check=False)
        # Construct the abtract group:
        gens_orders = tuple([ZZ(n)]+[ZZ(0)]*(self.__rank))
        AbelianGroupWithValues_class.__init__(self, gens_orders, 'u', values, number_field)

    def _element_constructor_(self, u):
        """
        Returns the abstract group element corresponding to the unit u.

        INPUT:

        - ``u`` -- Any object from which an element of the unit group's number
          field `K` may be constructed; an error is raised if an element of `K`
          cannot be constructed from u, or if the element constructed is not a
          unit.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2-38)
            sage: UK = UnitGroup(K)
            sage: UK(1)
            1
            sage: UK(-1)
            u0
            sage: UK.gens()
            (u0, u1)
            sage: UK.gens_values()
            [-1, 6*a - 37]
            sage: UK.ngens()
            2
            sage: [UK(u) for u in UK.gens()]
            [u0, u1]
            sage: [UK(u).exponents() for u in UK.gens()]
            [(1, 0), (0, 1)]
            sage: UK(a)
            Traceback (most recent call last):
            ...
            ValueError: a is not a unit
        """
        K = self.__number_field
        pK = self.__pari_number_field
        try:
            u = K(u)
        except TypeError:
            raise ValueError, "%s is not an element of %s"%(u,K)
        if self.__S:
            m = pK.bnfissunit(self.__S_unit_data, pari(u)).mattranspose()
            if m.ncols()==0:
                raise ValueError, "%s is not an S-unit"%u
        else:
            if not u.is_integral() or u.norm().abs() != 1:
                raise ValueError, "%s is not a unit"%u
            m = pK.bnfisunit(pari(u)).mattranspose()

        # convert column matrix to a list:
        m = [ZZ(m[0,i].python()) for i in range(m.ncols())]

        # NB pari puts the torsion after the fundamental units, before
        # the extra S-units but we have the torsion first:
        m = [m[self.__nfu]] + m[:self.__nfu] + m[self.__nfu+1:]

        return self.element_class(self, m)

    def rank(self):
        """
        Return the rank of the unit group.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(13)
            sage: UnitGroup(K).rank()
            5
            sage: SUK = UnitGroup(K,S=2); SUK.rank()
            6
        """
        return self.ngens()-1

    def _repr_(self):
        """
        Return string representation of this unit group.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: U = UnitGroup(NumberField(x^3 - 2, 'a'))
            sage: U
            Unit group with structure C2 x Z of Number Field in a with defining polynomial x^3 - 2
            sage: U._repr_()
            'Unit group with structure C2 x Z of Number Field in a with defining polynomial x^3 - 2'
            sage: UnitGroup(NumberField(x^3 - 2, 'a'),S=2)
            S-unit group with structure C2 x Z x Z of Number Field in a with defining polynomial x^3 - 2 with S = (Fractional ideal (a),)
        """
        if self.__S:
            return 'S-unit group with structure %s of %s with S = %s'%(
                self._group_notation(self.gens_orders()),
                self.number_field(),
                self.primes())
        return 'Unit group with structure %s of %s'%(
            self._group_notation(self.gens_orders()),
            self.number_field())

    def fundamental_units(self):
        """
        Return generators for the free part of the unit group, as a list.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^4 + 23)
            sage: U = UnitGroup(K)
            sage: U.fundamental_units()  # random
            [1/4*a^3 - 7/4*a^2 + 17/4*a - 19/4]
        """
        return self.gens_values()[1:]

    def roots_of_unity(self):
        """
        Return all the roots of unity in this unit group, primitive or not.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<b> = NumberField(x^2+1)
            sage: U = UnitGroup(K)
            sage: zs = U.roots_of_unity(); zs
            [-b, -1, b, 1]
            sage: [ z**U.zeta_order() for z in zs ]
            [1, 1, 1, 1]
        """
        z = self.gen(0).value()
        n = self.__ntu
        return [ z**k for k in range(1, n+1) ]

    def torsion_generator(self):
        """
        Return a generator for the torsion part of the unit group.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^4 - x^2 + 4)
            sage: U = UnitGroup(K)
            sage: U.torsion_generator()
            u0
            sage: U.torsion_generator().value() # random
            -1/4*a^3 - 1/4*a + 1/2
        """
        return self.gen(0)

    def zeta_order(self):
        """
        Returns the order of the torsion part of the unit group.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^4 - x^2 + 4)
            sage: U = UnitGroup(K)
            sage: U.zeta_order()
            6
        """
        return self.__ntu

    def zeta(self, n=2, all=False):
        """
        Return one, or a list of all, primitive n-th root of unity in this unit group.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<z> = NumberField(x^2 + 3)
            sage: U = UnitGroup(K)
            sage: U.zeta(1)
            1
            sage: U.zeta(2)
            -1
            sage: U.zeta(2, all=True)
            [-1]
            sage: U.zeta(3)
            -1/2*z - 1/2
            sage: U.zeta(3, all=True)
            [-1/2*z - 1/2, 1/2*z - 1/2]
            sage: U.zeta(4)
            Traceback (most recent call last):
            ...
            ValueError: n (=4) does not divide order of generator

            sage: r.<x> = QQ[]
            sage: K.<b> = NumberField(x^2+1)
            sage: U = UnitGroup(K)
            sage: U.zeta(4)
            -b
            sage: U.zeta(4,all=True)
            [-b, b]
            sage: U.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: n (=3) does not divide order of generator
            sage: U.zeta(3,all=True)
            []

        """
        N = self.__ntu
        K = self.number_field()
        n = ZZ(n)
        if n <= 0:
            raise ValueError, "n (=%s) must be positive"%n
        if n == 1:
            if all:
                return [K(1)]
            else:
                return K(1)
        elif n == 2:
            if all:
                return [K(-1)]
            else:
                return K(-1)
        if n.divides(N):
            z = self.torsion_generator().value() ** (N//n)
            if all:
                return [z**i for i in n.coprime_integers(n)]
            else:
                return z
        else:
            if all:
                return []
            else:
                raise ValueError, "n (=%s) does not divide order of generator"%n

    def number_field(self):
        """
        Return the number field associated with this unit group.

        EXAMPLES::

            sage: U = UnitGroup(QuadraticField(-23, 'w')); U
            Unit group with structure C2 of Number Field in w with defining polynomial x^2 + 23
            sage: U.number_field()
            Number Field in w with defining polynomial x^2 + 23
        """
        return self.__number_field


    def primes(self):
        """
        Return the (possibly empty) list of primes associated with this S-unit group.

        EXAMPLES::

            sage: K.<a> = QuadraticField(-23)
            sage: S = tuple(K.ideal(3).prime_factors()); S
            (Fractional ideal (3, 1/2*a - 1/2), Fractional ideal (3, 1/2*a + 1/2))
            sage: U = UnitGroup(K,S=tuple(S)); U
            S-unit group with structure C2 x Z x Z of Number Field in a with defining polynomial x^2 + 23 with S = (Fractional ideal (3, 1/2*a - 1/2), Fractional ideal (3, 1/2*a + 1/2))
            sage: U.primes() == S
            True
        """
        return self.__S


    def log(self, u):
        r"""
        Return the exponents of the unit ``u`` with respect to group generators.

        INPUT:

        - ``u`` -- Any object from which an element of the unit group's number
          field `K` may be constructed; an error is raised if an element of `K`
          cannot be constructed from u, or if the element constructed is not a
          unit.

        OUTPUT: a list of integers giving the exponents of ``u`` with
        respect to the unit group's basis.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<z> = CyclotomicField(13)
            sage: UK = UnitGroup(K)
            sage: [UK.log(u) for u in UK.gens()]
            [(1, 0, 0, 0, 0, 0),
             (0, 1, 0, 0, 0, 0),
             (0, 0, 1, 0, 0, 0),
             (0, 0, 0, 1, 0, 0),
             (0, 0, 0, 0, 1, 0),
             (0, 0, 0, 0, 0, 1)]
            sage: vec = [65,6,7,8,9,10]
            sage: unit = UK.exp(vec); unit  # random
            -253576*z^11 + 7003*z^10 - 395532*z^9 - 35275*z^8 - 500326*z^7 - 35275*z^6 - 395532*z^5 + 7003*z^4 - 253576*z^3 - 59925*z - 59925
            sage: UK.log(unit)
            (13, 6, 7, 8, 9, 10)

        An S-unit example::

           sage: SUK = UnitGroup(K,S=2)
           sage: v = (3,1,4,1,5,9,2)
           sage: u = SUK.exp(v); u
           -997204*z^11 - 2419728*z^10 - 413812*z^9 - 413812*z^8 - 2419728*z^7 - 997204*z^6 - 2129888*z^4 - 1616524*z^3 + 149364*z^2 - 1616524*z - 2129888
           sage: SUK.log(u)
           (3, 1, 4, 1, 5, 9, 2)
           sage: SUK.log(u) == v
           True
        """
        return self(u).exponents()

    def exp(self, exponents):
        r"""
        Return unit with given exponents with respect to group generators.

        INPUT:

        - ``u`` -- Any object from which an element of the unit
          group's number field `K` may be constructed; an error is
          raised if an element of `K` cannot be constructed from u, or
          if the element constructed is not a unit.

        OUTPUT: a list of integers giving the exponents of ``u`` with
        respect to the unit group's basis.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: K.<z> = CyclotomicField(13)
            sage: UK = UnitGroup(K)
            sage: [UK.log(u) for u in UK.gens()]
            [(1, 0, 0, 0, 0, 0),
             (0, 1, 0, 0, 0, 0),
             (0, 0, 1, 0, 0, 0),
             (0, 0, 0, 1, 0, 0),
             (0, 0, 0, 0, 1, 0),
             (0, 0, 0, 0, 0, 1)]
            sage: vec = [65,6,7,8,9,10]
            sage: unit = UK.exp(vec)
            sage: UK.log(unit)
            (13, 6, 7, 8, 9, 10)
            sage: UK.exp(UK.log(u)) == u.value()
            True

        An S-unit example::

           sage: SUK = UnitGroup(K,S=2)
           sage: v = (3,1,4,1,5,9,2)
           sage: u = SUK.exp(v); u
           -997204*z^11 - 2419728*z^10 - 413812*z^9 - 413812*z^8 - 2419728*z^7 - 997204*z^6 - 2129888*z^4 - 1616524*z^3 + 149364*z^2 - 1616524*z - 2129888
           sage: SUK.log(u)
           (3, 1, 4, 1, 5, 9, 2)
           sage: SUK.log(u) == v
           True
        """
        return prod([u**e for u,e in zip(self.gens_values(),exponents)], self.number_field().one_element())


