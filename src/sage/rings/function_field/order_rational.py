r"""
Orders of function fields: rational
"""

# ****************************************************************************
#       Copyright (C) 2010      William Stein <wstein@gmail.com>
#                     2011      Maarten Derickx <m.derickx.student@gmail.com>
#                     2011      Julian Rueth <julian.rueth@gmail.com>
#                     2017-2020 Kwankyu Lee
#                     2019      Brent Baccala
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

import sage.rings.abc

from sage.arith.functions import lcm
from sage.arith.misc import GCD as gcd
from sage.categories.euclidean_domains import EuclideanDomains
from sage.categories.principal_ideal_domains import PrincipalIdealDomains
from sage.rings.number_field.number_field_base import NumberField

from .ideal import FunctionFieldIdeal
from .ideal_rational import FunctionFieldIdeal_rational, FunctionFieldIdealInfinite_rational
from .order import FunctionFieldMaximalOrder, FunctionFieldMaximalOrderInfinite


class FunctionFieldMaximalOrder_rational(FunctionFieldMaximalOrder):
    """
    Maximal orders of rational function fields.

    INPUT:

    - ``field`` -- a function field

    EXAMPLES::

        sage: K.<t> = FunctionField(GF(19)); K
        Rational function field in t over Finite Field of size 19
        sage: R = K.maximal_order(); R
        Maximal order of Rational function field in t over Finite Field of size 19
    """
    def __init__(self, field):
        """
        Initialize.

        TESTS::

            sage: K.<t> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: TestSuite(O).run(skip='_test_gcd_vs_xgcd')
        """
        FunctionFieldMaximalOrder.__init__(self, field, ideal_class=FunctionFieldIdeal_rational,
                                           category=EuclideanDomains())

        self._populate_coercion_lists_(coerce_list=[field._ring])

        self._ring = field._ring
        self._gen = self(self._ring.gen())
        self._basis = (self.one(),)

    def _element_constructor_(self, f):
        """
        Make ``f`` a function field element of this order.

        EXAMPLES::

            sage: K.<y> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: O._element_constructor_(y)
            y
            sage: O._element_constructor_(1/y)
            Traceback (most recent call last):
            ...
            TypeError: 1/y is not an element of Maximal order of Rational function field in y over Rational Field
        """
        F = self.function_field()
        try:
            f = F(f)
        except TypeError:
            raise TypeError("unable to convert to an element of {}".format(F))

        if f.denominator() not in self.function_field().constant_base_field():
            raise TypeError("%r is not an element of %r" % (f,self))

        return f

    def ideal_with_gens_over_base(self, gens):
        """
        Return the fractional ideal with generators ``gens``.

        INPUT:

        - ``gens`` -- elements of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)                                                # needs sage.rings.function_field
            sage: O = L.equation_order()                                                            # needs sage.rings.function_field
            sage: O.ideal_with_gens_over_base([x^3 + 1, -y])                                        # needs sage.rings.function_field
            Ideal (x^3 + 1, -y) of Order in Function field in y defined by y^2 - x^3 - 1
        """
        return self.ideal(gens)

    def _residue_field(self, ideal, name=None):
        """
        Return a field isomorphic to the residue field at the prime ideal.

        The residue field is by definition `k[x]/q` where `q` is the irreducible
        polynomial generating the prime ideal and `k` is the constant base field.

        INPUT:

        - ``ideal`` -- prime ideal of the order

        - ``name`` -- string; name of the generator of the residue field

        OUTPUT:

        - a field isomorphic to the residue field

        - a morphism from the field to `k[x]` via the residue field

        - a morphism from `k[x]` to the field via the residue field

        EXAMPLES::

            sage: # needs sage.modules sage.rings.finite_rings
            sage: F.<x> = FunctionField(GF(2))
            sage: O = F.maximal_order()
            sage: I = O.ideal(x^2 + x + 1)
            sage: R, fr_R, to_R = O._residue_field(I)
            sage: R
            Finite Field in z2 of size 2^2
            sage: [to_R(fr_R(e)) == e for e in R]
            [True, True, True, True]
            sage: [to_R(fr_R(e)).parent() is R for e in R]
            [True, True, True, True]
            sage: e1, e2 = fr_R(R.random_element()), fr_R(R.random_element())
            sage: to_R(e1 * e2) == to_R(e1) * to_R(e2)
            True
            sage: to_R(e1 + e2) == to_R(e1) + to_R(e2)
            True
            sage: to_R(e1).parent() is R
            True
            sage: to_R(e2).parent() is R
            True

            sage: # needs sage.modules sage.rings.finite_rings
            sage: F.<x> = FunctionField(GF(2))
            sage: O = F.maximal_order()
            sage: I = O.ideal(x + 1)
            sage: R, fr_R, to_R = O._residue_field(I)
            sage: R
            Finite Field of size 2
            sage: [to_R(fr_R(e)) == e for e in R]
            [True, True]
            sage: [to_R(fr_R(e)).parent() is R for e in R]
            [True, True]
            sage: e1, e2 = fr_R(R.random_element()), fr_R(R.random_element())
            sage: to_R(e1 * e2) == to_R(e1) * to_R(e2)
            True
            sage: to_R(e1 + e2) == to_R(e1) + to_R(e2)
            True
            sage: to_R(e1).parent() is R
            True
            sage: to_R(e2).parent() is R
            True

            sage: # needs sage.modules sage.rings.number_field
            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x^2 + x + 1)
            sage: R, fr_R, to_R = O._residue_field(I)
            sage: R
            Number Field in a with defining polynomial x^2 + x + 1
            sage: e1, e2 = fr_R(R.random_element()), fr_R(R.random_element())
            sage: to_R(e1 * e2) == to_R(e1) * to_R(e2)
            True
            sage: to_R(e1 + e2) == to_R(e1) + to_R(e2)
            True
            sage: to_R(e1).parent() is R
            True
            sage: to_R(e2).parent() is R
            True

            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x + 1)
            sage: R, fr_R, to_R = O._residue_field(I)
            sage: R
            Rational Field
            sage: e1, e2 = fr_R(R.random_element()), fr_R(R.random_element())
            sage: to_R(e1 * e2) == to_R(e1) * to_R(e2)
            True
            sage: to_R(e1 + e2) == to_R(e1) + to_R(e2)
            True
            sage: to_R(e1).parent() is R
            True
            sage: to_R(e2).parent() is R
            True
        """
        F = self.function_field()
        K = F.constant_base_field()

        q = ideal.gen().element().numerator()

        if F.is_global():
            R, _from_R, _to_R = self._residue_field_global(q, name=name)
        elif isinstance(K, (NumberField, sage.rings.abc.AlgebraicField)):
            if name is None:
                name = 'a'
            if q.degree() == 1:
                R = K
                _from_R = lambda e: e
                _to_R = lambda e: R(e % q)
            else:
                R = K.extension(q, names=name)
                _from_R = lambda e: self._ring(list(R(e)))
                _to_R = lambda e: (e % q)(R.gen(0))
        else:
            raise NotImplementedError

        def from_R(e):
            return F(_from_R(e))

        def to_R(f):
            return _to_R(f.numerator())

        return R, from_R, to_R

    def _residue_field_global(self, q, name=None):
        """
        Return a finite field isomorphic to the residue field at q.

        This method assumes a global rational function field, that is,
        the constant base field is a finite field.

        INPUT:

        - ``q`` -- irreducible polynomial

        - ``name`` -- string; name of the generator of the extension field

        OUTPUT:

        - a finite field

        - a function that outputs a polynomial lifting a finite field element

        - a function that outputs a finite field element for a polynomial

        The residue field is by definition `k[x]/q` where `k` is the base field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(4)
            sage: F.<x> = FunctionField(k)
            sage: O = F.maximal_order()
            sage: O._ring
            Univariate Polynomial Ring in x over Finite Field in a of size 2^2
            sage: f = x^3 + x + 1
            sage: _f = f.numerator()
            sage: _f.is_irreducible()
            True
            sage: K, fr_K, to_K = O._residue_field_global(_f)                                       # needs sage.modules
            sage: K                                                                                 # needs sage.modules
            Finite Field in z6 of size 2^6
            sage: all(to_K(fr_K(e)) == e for e in K)                                                # needs sage.modules
            True

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(2)
            sage: F.<x> = FunctionField(k)
            sage: O = F.maximal_order()
            sage: O._ring                                                                           # needs sage.libs.ntl
            Univariate Polynomial Ring in x over Finite Field of size 2 (using GF2X)
            sage: f = x^3 + x + 1
            sage: _f = f.numerator()
            sage: _f.is_irreducible()
            True
            sage: K, fr_K, to_K = O._residue_field_global(_f)                                       # needs sage.modules
            sage: all(to_K(fr_K(e)) == e for e in K)                                                # needs sage.modules
            True
        """
        # polynomial ring over the base field
        R = self._ring

        # base field of extension degree r over the prime field
        k = R.base_ring()
        a = k.gen()
        r = k.degree()

        # extend the base field to a field of degree r*s over the
        # prime field
        s = q.degree()
        K,sigma = k.extension(s, map=True, name=name)

        # find a root beta in K satisfying the irreducible q
        S = K['X']
        beta = S([sigma(c) for c in q.list()]).roots()[0][0]

        # V is a vector space over the prime subfield of k of degree r*s
        V = K.vector_space(map=False)

        w = K.one()
        beta_pow = []
        for i in range(s):
            beta_pow.append(w)
            w *= beta

        w = K.one()
        sigma_a = sigma(a)
        sigma_a_pow = []
        for i in range(r):
            sigma_a_pow.append(w)
            w *= sigma_a

        basis = [V(sap * bp) for bp in beta_pow for sap in sigma_a_pow]
        W = V.span_of_basis(basis)

        def to_K(f):
            coeffs = (f % q).list()
            return sum((sigma(c) * beta_pow[i] for i, c in enumerate(coeffs)), K.zero())

        if r == 1: # take care of the prime field case
            def fr_K(g):
                co = W.coordinates(V(g), check=False)
                return R([k(co[j]) for j in range(s)])
        else:
            def fr_K(g):
                co = W.coordinates(V(g), check=False)
                return R([k(co[i:i+r]) for i in range(0, r*s, r)])

        return K, fr_K, to_K

    def basis(self):
        """
        Return the basis (=1) of the order as a module over the polynomial ring.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(19))
            sage: O = K.maximal_order()
            sage: O.basis()
            (1,)
        """
        return self._basis

    def gen(self, n=0):
        """
        Return the ``n``-th generator of the order. Since there is only one generator ``n`` must be 0.

        EXAMPLES::

            sage: O = FunctionField(QQ,'y').maximal_order()
            sage: O.gen()
            y
            sage: O.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: there is only one generator
        """
        if n != 0:
            raise IndexError("there is only one generator")
        return self._gen

    def ngens(self):
        """
        Return 1 the number of generators of the order.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().ngens()
            1
        """
        return 1

    def ideal(self, *gens):
        """
        Return the fractional ideal generated by ``gens``.

        INPUT:

        - ``gens`` -- elements of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: O.ideal(x)
            Ideal (x) of Maximal order of Rational function field in x over Rational Field
            sage: O.ideal([x, 1/x]) == O.ideal(x, 1/x)  # multiple generators may be given as a list
            True
            sage: O.ideal(x^3 + 1, x^3 + 6)
            Ideal (1) of Maximal order of Rational function field in x over Rational Field
            sage: I = O.ideal((x^2+1)*(x^3+1), (x^3+6)*(x^2+1)); I
            Ideal (x^2 + 1) of Maximal order of Rational function field in x over Rational Field
            sage: O.ideal(I)
            Ideal (x^2 + 1) of Maximal order of Rational function field in x over Rational Field
        """
        if len(gens) == 1:
            gens = gens[0]
            if not isinstance(gens, (list, tuple)):
                if isinstance(gens, FunctionFieldIdeal):
                    gens = gens.gens()
                else:
                    gens = (gens,)
        K = self.function_field()
        gens = [K(e) for e in gens if e != 0]
        if len(gens) == 0:
            gen = K(0)
        else:
            d = lcm([c.denominator() for c in gens]).monic()
            g = gcd([(d*c).numerator() for c in gens]).monic()
            gen = K(g/d)

        return self.ideal_monoid().element_class(self, gen)


class FunctionFieldMaximalOrderInfinite_rational(FunctionFieldMaximalOrderInfinite):
    """
    Maximal infinite orders of rational function fields.

    INPUT:

    - ``field`` -- a rational function field

    EXAMPLES::

        sage: K.<t> = FunctionField(GF(19)); K
        Rational function field in t over Finite Field of size 19
        sage: R = K.maximal_order_infinite(); R
        Maximal infinite order of Rational function field in t over Finite Field of size 19
    """
    def __init__(self, field, category=None):
        """
        Initialize.

        TESTS::

            sage: K.<t> = FunctionField(GF(19))
            sage: O = K.maximal_order_infinite()
            sage: TestSuite(O).run(skip='_test_gcd_vs_xgcd')
        """
        FunctionFieldMaximalOrderInfinite.__init__(self, field, ideal_class=FunctionFieldIdealInfinite_rational,
                                                   category=PrincipalIdealDomains().or_subcategory(category))
        self._populate_coercion_lists_(coerce_list=[field.constant_base_field()])

    def _element_constructor_(self, f):
        """
        Make ``f`` an element of this order.

        EXAMPLES::

            sage: K.<y> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: O._element_constructor_(y)
            y
            sage: O._element_constructor_(1/y)
            Traceback (most recent call last):
            ...
            TypeError: 1/y is not an element of Maximal order of Rational function field in y over Rational Field
        """
        F = self.function_field()
        try:
            f = F(f)
        except TypeError:
            raise TypeError("unable to convert to an element of {}".format(F))

        if f.denominator().degree() < f.numerator().degree():
            raise TypeError("{} is not an element of {}".format(f, self))

        return f

    def basis(self):
        """
        Return the basis (=1) of the order as a module over the polynomial ring.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(19))
            sage: O = K.maximal_order()
            sage: O.basis()
            (1,)
        """
        return 1/self.function_field().gen()

    def gen(self, n=0):
        """
        Return the `n`-th generator of ``self``. Since there is only one
        generator `n` must be `0`.

        EXAMPLES::

            sage: O = FunctionField(QQ,'y').maximal_order()
            sage: O.gen()
            y
            sage: O.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: there is only one generator
        """
        if n != 0:
            raise IndexError("there is only one generator")
        return self._gen

    def ngens(self):
        """
        Return 1 the number of generators of the order.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().ngens()
            1
        """
        return 1

    def prime_ideal(self):
        """
        Return the unique prime ideal of the order.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(19))
            sage: O = K.maximal_order_infinite()
            sage: O.prime_ideal()
            Ideal (1/t) of Maximal infinite order of Rational function field in t
            over Finite Field of size 19
        """
        return self.ideal( 1/self.function_field().gen() )

    def ideal(self, *gens):
        """
        Return the fractional ideal generated by ``gens``.

        INPUT:

        - ``gens`` -- elements of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order_infinite()
            sage: O.ideal(x)
            Ideal (x) of Maximal infinite order of Rational function field in x over Rational Field
            sage: O.ideal([x, 1/x]) == O.ideal(x ,1/x)  # multiple generators may be given as a list
            True
            sage: O.ideal(x^3 + 1, x^3 + 6)
            Ideal (x^3) of Maximal infinite order of Rational function field in x over Rational Field
            sage: I = O.ideal((x^2+1)*(x^3+1), (x^3+6)*(x^2+1)); I
            Ideal (x^5) of Maximal infinite order of Rational function field in x over Rational Field
            sage: O.ideal(I)
            Ideal (x^5) of Maximal infinite order of Rational function field in x over Rational Field
        """
        if len(gens) == 1:
            gens = gens[0]
            if not isinstance(gens, (list, tuple)):
                if isinstance(gens, FunctionFieldIdeal):
                    gens = gens.gens()
                else:
                    gens = (gens,)
        K = self.function_field()
        gens = [K(g) for g in gens]
        try:
            d = max(g.numerator().degree() - g.denominator().degree() for g in gens if g != 0)
            gen = K.gen() ** d
        except ValueError: # all gens are zero
            gen = K(0)

        return self.ideal_monoid().element_class(self, gen)
