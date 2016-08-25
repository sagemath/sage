"""
Finite residue fields

We can take the residue field of maximal ideals in the ring of integers
of number fields. We can also take the residue field of irreducible
polynomials over `GF(p)`.

EXAMPLES::

    sage: K.<a> = NumberField(x^3-7)
    sage: P = K.ideal(29).factor()[0][0]
    sage: k = K.residue_field(P)
    sage: k
    Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)
    sage: k.order()
    841

We reduce mod a prime for which the ring of integers is not
monogenic (i.e., 2 is an essential discriminant divisor)::

    sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
    sage: F = K.factor(2); F
    (Fractional ideal (1/2*a^2 - 1/2*a + 1)) * (Fractional ideal (-a^2 + 2*a - 3)) * (Fractional ideal (-3/2*a^2 + 5/2*a - 4))
    sage: F[0][0].residue_field()
    Residue field of Fractional ideal (1/2*a^2 - 1/2*a + 1)
    sage: F[1][0].residue_field()
    Residue field of Fractional ideal (-a^2 + 2*a - 3)
    sage: F[2][0].residue_field()
    Residue field of Fractional ideal (-3/2*a^2 + 5/2*a - 4)

We can also form residue fields from `\ZZ`::

    sage: ZZ.residue_field(17)
    Residue field of Integers modulo 17

And for polynomial rings over finite fields::

    sage: R.<t> = GF(5)[]
    sage: I = R.ideal(t^2 + 2)
    sage: k = ResidueField(I); k
    Residue field in tbar of Principal ideal (t^2 + 2) of Univariate Polynomial Ring in t over Finite Field of size 5

AUTHORS:

- David Roe (2007-10-3): initial version
- William Stein (2007-12): bug fixes
- John Cremona (2008-9): extend reduction maps to the whole valuation ring
  add support for residue fields of ZZ
- David Roe (2009-12): added support for `GF(p)(t)` and moved to new coercion
  framework.

TESTS::

    sage: K.<z> = CyclotomicField(7)
    sage: P = K.factor(17)[0][0]
    sage: ff = K.residue_field(P)
    sage: loads(dumps(ff)) is ff
    True
    sage: a = ff(z)
    sage: parent(a*a)
    Residue field in zbar of Fractional ideal (17)
    sage: TestSuite(ff).run()

Verify that :trac:`15192` has been resolved::

    sage: a.is_unit()
    True

    sage: R.<t> = GF(11)[]; P = R.ideal(t^3 + t + 4)
    sage: ff.<a> = ResidueField(P)
    sage: a == ff(t)
    True
    sage: parent(a*a)
    Residue field in a of Principal ideal (t^3 + t + 4) of Univariate Polynomial Ring in t over Finite Field of size 11

Verify that :trac:`7475` is fixed::

    sage: K = ZZ.residue_field(2)
    sage: loads(dumps(K)) is K
    True

Reducing a curve modulo a prime::

    sage: K.<s> = NumberField(x^2+23)
    sage: OK = K.ring_of_integers()
    sage: E = EllipticCurve([0,0,0,K(1),K(5)])
    sage: pp = K.factor(13)[0][0]
    sage: Fpp = OK.residue_field(pp)
    sage: E.base_extend(Fpp)
    Elliptic Curve defined by y^2  = x^3 + x + 5 over Residue field of Fractional ideal (13, 1/2*s + 9/2)

    sage: R.<t> = GF(11)[]
    sage: P = R.ideal(t^3 + t + 4)
    sage: ff.<a> = R.residue_field(P)
    sage: E = EllipticCurve([0,0,0,R(1),R(t)])
    sage: E.base_extend(ff)
    Elliptic Curve defined by y^2 = x^3 + x + a over Residue field in a of Principal ideal (t^3 + t + 4) of Univariate Polynomial Ring in t over Finite Field of size 11

Calculating Groebner bases over various residue fields.
First over a small non-prime field::

    sage: F1.<u> = NumberField(x^6 + 6*x^5 + 124*x^4 + 452*x^3 + 4336*x^2 + 8200*x + 42316)
    sage: reduct_id = F1.factor(47)[0][0]
    sage: Rf = F1.residue_field(reduct_id)
    sage: type(Rf)
    <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_pari_ffelt_with_category'>
    sage: Rf.cardinality().factor()
    47^3
    sage: R.<X, Y> = PolynomialRing(Rf)
    sage: ubar = Rf(u)
    sage: I = ideal([ubar*X + Y]); I
    Ideal ((ubar)*X + Y) of Multivariate Polynomial Ring in X, Y over Residue field in ubar of Fractional ideal (47, 517/55860*u^5 + 235/3724*u^4 + 9829/13965*u^3 + 54106/13965*u^2 + 64517/27930*u + 755696/13965)
    sage: I.groebner_basis()
    [X + (-19*ubar^2 - 5*ubar - 17)*Y]

And now over a large prime field::

    sage: x = ZZ['x'].0
    sage: F1.<u> = NumberField(x^2 + 6*x + 324)
    sage: reduct_id = F1.prime_above(next_prime(2^42))
    sage: Rf = F1.residue_field(reduct_id)
    sage: type(Rf)
    <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_prime_modn_with_category'>
    sage: Rf.cardinality().factor()
    4398046511119
    sage: S.<X, Y, Z> = PolynomialRing(Rf, order='lex')
    sage: I = ideal([2*X - Y^2, Y + Z])
    sage: I.groebner_basis()
    verbose 0 (...: multi_polynomial_ideal.py, groebner_basis) Warning: falling back to very slow toy implementation.
    [X + 2199023255559*Z^2, Y + Z]
    sage: S.<X, Y, Z> = PolynomialRing(Rf, order='deglex')
    sage: I = ideal([2*X - Y^2, Y + Z])
    sage: I.groebner_basis()
    verbose 0 (...: multi_polynomial_ideal.py, groebner_basis) Warning: falling back to very slow toy implementation.
    [Z^2 + 4398046511117*X, Y + Z]
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.ring cimport Field
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.categories.homset import Hom
from sage.rings.all import ZZ, QQ, Integers
from sage.rings.finite_rings.finite_field_constructor import zech_log_bound, FiniteField as GF
from sage.rings.finite_rings.finite_field_givaro import FiniteField_givaro
from sage.rings.finite_rings.finite_field_ntl_gf2e import FiniteField_ntl_gf2e
from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn
from sage.rings.finite_rings.finite_field_pari_ffelt import FiniteField_pari_ffelt
from sage.rings.ideal import is_Ideal
from sage.structure.element cimport Element

from sage.rings.number_field.number_field_element import is_NumberFieldElement
from sage.rings.number_field.number_field_ideal import is_NumberFieldIdeal

from sage.modules.free_module_element import FreeModuleElement
from sage.rings.fraction_field import is_FractionField

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.polynomial_element import is_Polynomial

from sage.structure.factory import UniqueFactory
from sage.structure.element cimport parent_c


class ResidueFieldFactory(UniqueFactory):
    """
    A factory that returns the residue class field of a prime ideal `p`
    of the ring of integers of a number field, or of a polynomial ring
    over a finite field.

    INPUT:

        - ``p`` -- a prime ideal of an order in a number field.

        - ``names`` -- the variable name for the finite field created.
          Defaults to the name of the number field variable but with
          bar placed after it.

        - ``check`` -- whether or not to check if `p` is prime.

    OUTPUT:

         - The residue field at the prime `p`.

    EXAMPLES::

        sage: K.<a> = NumberField(x^3-7)
        sage: P = K.ideal(29).factor()[0][0]
        sage: ResidueField(P)
        Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)

    The result is cached::

        sage: ResidueField(P) is ResidueField(P)
        True
        sage: k = K.residue_field(P); k
        Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)
        sage: k.order()
        841

    It also works for polynomial rings::

        sage: R.<t> = GF(31)[]
        sage: P = R.ideal(t^5 + 2*t + 11)
        sage: ResidueField(P)
        Residue field in tbar of Principal ideal (t^5 + 2*t + 11) of Univariate Polynomial Ring in t over Finite Field of size 31

        sage: ResidueField(P) is ResidueField(P)
        True
        sage: k = ResidueField(P); k.order()
        28629151

    An example where the generator of the number field doesn't
    generate the residue class field::

        sage: K.<a> = NumberField(x^3-875)
        sage: P = K.ideal(5).factor()[0][0]; k = K.residue_field(P); k
        Residue field in abar of Fractional ideal (5, 1/25*a^2 - 2/5*a - 1)
        sage: k.polynomial()
        abar^2 + 3*abar + 4
        sage: k.0^3 - 875
        2

    An example where the residue class field is large but of degree 1::

        sage: K.<a> = NumberField(x^3-875); P = K.ideal(2007).factor()[2][0]; k = K.residue_field(P); k
        Residue field of Fractional ideal (223, 1/5*a + 11)
        sage: k(a)
        168
        sage: k(a)^3 - 875
        0

    And for polynomial rings::

        sage: R.<t> = GF(next_prime(2^18))[]
        sage: P = R.ideal(t - 5)
        sage: k = ResidueField(P); k
        Residue field of Principal ideal (t + 262142) of Univariate Polynomial Ring in t over Finite Field of size 262147
        sage: k(t)
        5

    In this example, 2 is an inessential discriminant divisor, so divides
    the index of ``ZZ[a]`` in the maximal order for all ``a``::

        sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8); P = K.ideal(2).factor()[0][0]; P
        Fractional ideal (1/2*a^2 - 1/2*a + 1)
        sage: F = K.residue_field(P); F
        Residue field of Fractional ideal (1/2*a^2 - 1/2*a + 1)
        sage: F(a)
        0
        sage: B = K.maximal_order().basis(); B
        [1, 1/2*a^2 + 1/2*a, a^2]
        sage: F(B[1])
        1
        sage: F(B[2])
        0
        sage: F
        Residue field of Fractional ideal (1/2*a^2 - 1/2*a + 1)
        sage: F.degree()
        1

    TESTS::

        sage: K.<a> = NumberField(polygen(QQ))
        sage: K.residue_field(K.ideal(3))
        Residue field of Fractional ideal (3)
    """
    def create_key_and_extra_args(self, p, names = None, check=True, impl=None, **kwds):
        """
        Return a tuple containing the key (uniquely defining data)
        and any extra arguments.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-7)
            sage: ResidueField(K.ideal(29).factor()[0][0]) # indirect doctest
            Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)
        """
        if check:
            if not is_Ideal(p):
                if isinstance(p, (int, Integer, Rational)):
                    p = ZZ.ideal(p)
                elif is_NumberFieldElement(p):
                    if p.parent().is_field():
                        p = p.parent().ring_of_integers().ideal(p)
                    else:
                        p = p.parent().ideal(p)
                elif is_Polynomial(p):
                    p = p.parent().ideal(p)
                #elif isinstance(p.parent(), FractionField_1poly_field):
                #    p = p.parent().ring_of_integers().ideal(p)
                # will eventually support other function fields here.
                else:
                    raise ValueError("p must be an ideal or element of a number field or function field.")
            if not p.is_prime():
                raise ValueError("p (%s) must be prime" % p)
            if is_PolynomialRing(p.ring()):
                if not p.ring().base_ring().is_finite():
                    raise ValueError("residue fields only supported for polynomial rings over finite fields")
                if not p.ring().base_ring().is_prime_field():
                    # neither of these will work over non-prime fields quite yet.  We should use relative finite field extensions.
                    raise NotImplementedError
            elif not (is_NumberFieldIdeal(p) or p.ring() is ZZ):
                raise NotImplementedError
        if isinstance(names, tuple):
            if len(names) > 0:
                names = str(names[0])
            else:
                names = None
        if names is None and p.ring() is not ZZ:
            names = '%sbar'%(p.ring().fraction_field().variable_name())
        key = (p, names, impl)
        return key, kwds

    def create_object(self, version, key, **kwds):
        """
        Create the object from the key and extra arguments. This is only
        called if the object was not found in the cache.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(29).factor()[0][0]
            sage: ResidueField(P) is ResidueField(P) # indirect doctest
            True
        """
        p, names, impl = key
        pring = p.ring()

        if pring is ZZ:
            return ResidueFiniteField_prime_modn(p, names, p.gen(), None, None, None)
        if is_PolynomialRing(pring):
            K = pring.fraction_field()
            Kbase = pring.base_ring()
            f = p.gen()
            characteristic = Kbase.order()
            if f.degree() == 1 and Kbase.is_prime_field() and (impl is None or impl == 'modn'):
                return ResidueFiniteField_prime_modn(p, None, Kbase.order(), None, None, None)
            else:
                q = characteristic**(f.degree())
                if q < zech_log_bound and (impl is None or impl == 'givaro'):
                    return ResidueFiniteField_givaro(p, q, names, f, None, None, None)
                elif (q % 2 == 0) and (impl is None or impl == 'ntl'):
                    return ResidueFiniteField_ntl_gf2e(q, names, f, "poly", p, None, None, None)
                elif impl is None or impl == 'pari':
                    return ResidueFiniteField_pari_ffelt(p, characteristic, names, f, None, None, None)
                else:
                    raise ValueError("unrecognized finite field type")

        # Should generalize to allowing residue fields of relative extensions to be extensions of finite fields.
        if is_NumberFieldIdeal(p):
            characteristic = p.smallest_integer()
        else: # ideal of a function field
            characteristic = pring.base_ring().characteristic()
        # Once we have function fields, we should probably have an if statement here.
        K = pring.fraction_field()
        #OK = K.maximal_order() # Need to change to p.order inside the __init__s for the residue fields.

        U, to_vs, to_order = p._p_quotient(characteristic)
        k = U.base_ring()
        R = PolynomialRing(k, names)
        n = p.residue_class_degree()
        gen_ok = False
        from sage.matrix.constructor import matrix
        try:
            x = K.gen()
            if not x:
                LL = [to_vs(1).list()] + [to_vs(x**i).list() for i in range(1,n+1)]
                M = matrix(k, n+1, n, LL)
            else:
                M = matrix(k, n+1, n, [to_vs(x**i).list() for i in range(n+1)])

            W = M.transpose().echelon_form()
            if M.rank() == n:
                PB = M.matrix_from_rows(range(n))
                gen_ok = True
                f = R((-W.column(n)).list() + [1])
        except (TypeError, ZeroDivisionError):
            pass
        if not gen_ok:
            bad = True
            for u in U: # using this iterator may not be optimal, we may get a long string of non-generators
                if u:
                    x = to_order(u)
                    M = matrix(k, n+1, n, [to_vs(x**i).list() for i in range(n+1)])
                    W = M.transpose().echelon_form()
                    if W.rank() == n:
                        f = R((-W.column(n)).list() + [1])
                        PB = M.matrix_from_rows(range(n))
                        bad = False
                        break
            assert not bad, "error -- didn't find a generator."
        # The reduction map is just x |--> k(to_vs(x) * (PB**(-1)))
        # The lifting map is just x |--> to_order(x * PB)
        # These are constructed inside the field __init__
        if n == 1:
            return ResidueFiniteField_prime_modn(p, names, p.smallest_integer(), to_vs, to_order, PB)
        else:
            q = characteristic**(f.degree())
            if q < zech_log_bound and (impl is None or impl == 'givaro'):
                return ResidueFiniteField_givaro(p, q, names, f, to_vs, to_order, PB)
            elif (q % 2 == 0) and (impl is None or impl == 'ntl'):
                return ResidueFiniteField_ntl_gf2e(q, names, f, "poly", p, to_vs, to_order, PB)
            elif impl is None or impl == 'pari':
                return ResidueFiniteField_pari_ffelt(p, characteristic, names, f, to_vs, to_order, PB)
            else:
                raise ValueError("unrecognized finite field type")

ResidueField = ResidueFieldFactory("ResidueField")

class ResidueField_generic(Field):
    """
    The class representing a generic residue field.

    EXAMPLES::

        sage: I = QQ[i].factor(2)[0][0]; I
        Fractional ideal (I + 1)
        sage: k = I.residue_field(); k
        Residue field of Fractional ideal (I + 1)
        sage: type(k)
        <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_prime_modn_with_category'>

        sage: R.<t> = GF(29)[]; P = R.ideal(t^2 + 2); k.<a> = ResidueField(P); k
        Residue field in a of Principal ideal (t^2 + 2) of Univariate Polynomial Ring in t over Finite Field of size 29
        sage: type(k)
        <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_givaro_with_category'>
    """
    def __init__(self, p):
        """
        .. WARNING::

            This function does not call up the ``__init__`` chain, since many
            residue fields use multiple inheritance and will be calling
            ``__init__`` via their other superclass.

            If this is not the case, one should call ``Parent.__init__``
            manually for any subclass.

        INPUT:

           - ``p`` -- the prime (ideal) defining this residue field

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-17)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k = K.residue_field(P) # indirect doctest
            sage: F = ZZ.residue_field(17)  # indirect doctest

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 7)
            sage: k.<a> = P.residue_field() # indirect doctest

            sage: k.category()
            Category of finite fields
            sage: F.category()
            Join of Category of finite fields and Category of subquotients of monoids and Category of quotients of semigroups

        TESTS::

            sage: TestSuite(k).run()
            sage: TestSuite(F).run()
        """
        self.p = p
        # Note: we don't call Parent.__init__ since many residue fields use multiple inheritance and will be calling __init__ via their other superclass.

    def ideal(self):
        r"""
        Return the maximal ideal that this residue field is the quotient by.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k = K.residue_field(P) # indirect doctest
            sage: k.ideal() is P
            True
            sage: p = next_prime(2^40); p
            1099511627791
            sage: k = K.residue_field(K.prime_above(p))
            sage: k.ideal().norm() == p
            True

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 7)
            sage: k.<a> = R.residue_field(P)
            sage: k.ideal()
            Principal ideal (t^3 + t^2 + 7) of Univariate Polynomial Ring in t over Finite Field of size 17
        """
        return self.p

    def _element_constructor_(self, x):
        """
        This is called after ``x`` fails to convert into ``self`` as
        abstract finite field (without considering the underlying
        number field).

        So the strategy is to try to convert into the number field,
        and then proceed to the residue field.

        .. NOTE::

            The behaviour of this method was changed in :trac:`8800`.
            Before, an error was raised if there was no coercion. Now,
            a conversion is possible even when there is no coercion.
            This is like for different finite fields.

        EXAMPLES::

            sage: from sage.rings.finite_rings.residue_field import ResidueField_generic
            sage: K.<i> = NumberField(x^2+1)
            sage: P = K.ideal(-3*i-2)
            sage: OK = K.maximal_order()
            sage: F = OK.residue_field(P)
            sage: ResidueField_generic._element_constructor_(F, i)
            8

        With :trac:`8800`, we also have::

            sage: ResidueField_generic._element_constructor_(F, GF(13)(8))
            8

        Here is a test that was temporarily removed, but newly introduced
        in :trac:`8800`::

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 7)
            sage: k.<a> = P.residue_field()
            sage: k(t)
            a
            sage: k(GF(17)(4))
            4
        """
        K = OK = self.p.ring()
        R = parent_c(x)
        if OK.is_field():
            OK = OK.ring_of_integers()
        else:
            K = K.fraction_field()
        if OK.has_coerce_map_from(R):
            x = OK(x)
        elif K.has_coerce_map_from(R):
            x = K(x)
        else:
            try:
                x = K(x)
            except (TypeError, ValueError):
                raise TypeError("cannot coerce %s" % type(x))
        return self(x)

    def _coerce_map_from_(self, R):
        """
        Returns ``True`` if there is a coercion map from ``R`` to ``self``.

        EXAMPLES::

            sage: K.<i> = NumberField(x^2+1)
            sage: P = K.ideal(-3*i-2)
            sage: OK = K.maximal_order()
            sage: F = OK.residue_field(P)
            sage: F.has_coerce_map_from(GF(13)) # indirect doctest
            True

        TESTS:

        Check that :trac:`11319` is fixed::

            sage: GF(13).has_coerce_map_from(F)
            True

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 7)
            sage: k.<a> = P.residue_field()
            sage: k.has_coerce_map_from(Qp(17)) # indirect doctest
            False
        """
        OK = self.p.ring()
        if OK.is_field():
            OK = OK.ring_of_integers()
        return self.base_ring().has_coerce_map_from(R) or OK.has_coerce_map_from(R)

    def __repr__(self):
        """
        Returns a string describing this residue field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k = K.residue_field(P)
            sage: k
            Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)

            sage: F = ZZ.residue_field(17); F
            Residue field of Integers modulo 17

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 7)
            sage: k.<a> = P.residue_field(); k # indirect doctest
            Residue field in a of Principal ideal (t^3 + t^2 + 7) of Univariate Polynomial Ring in t over Finite Field of size 17
        """
        if self.p.ring() is ZZ:
            return "Residue field of Integers modulo %s"%self.p.gen()
        return "Residue field %sof %s"%('in %s '%self.gen() if self.degree() > 1 else '', self.p)

    def lift(self, x):
        """
        Returns a lift of ``x`` to the Order, returning a "polynomial" in the
        generator with coefficients between 0 and `p-1`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k =K.residue_field(P)
            sage: OK = K.maximal_order()
            sage: c = OK(a)
            sage: b = k(a)
            sage: k.lift(13*b + 5)
            13*a + 5
            sage: k.lift(12821*b+918)
            3*a + 19

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 7)
            sage: k.<a> = P.residue_field()
            sage: k.lift(a^2 + 5)
            t^2 + 5
        """
        if hasattr(self.p, "ring"):
            R = self.p.ring()
            if R.is_field():
                R = R.ring_of_integers()
            return R(x)
        else:
            return x.lift()

    def reduction_map(self):
        """
        Return the partially defined reduction map from the number
        field to this residue class field.

        EXAMPLES::

            sage: I = QQ[2^(1/3)].factor(2)[0][0]; I
            Fractional ideal (a)
            sage: k = I.residue_field(); k
            Residue field of Fractional ideal (a)
            sage: pi = k.reduction_map(); pi
            Partially defined reduction map:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Residue field of Fractional ideal (a)
            sage: pi.domain()
            Number Field in a with defining polynomial x^3 - 2
            sage: pi.codomain()
            Residue field of Fractional ideal (a)

            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 32)
            sage: F = K.factor(2)[0][0].residue_field()
            sage: F.reduction_map().domain()
            Number Field in a with defining polynomial x^3 + x^2 - 2*x + 32
            sage: K.<a> = NumberField(x^3 + 128)
            sage: F = K.factor(2)[0][0].residue_field()
            sage: F.reduction_map().codomain()
            Residue field of Fractional ideal (1/4*a)

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 7)
            sage: k.<a> = P.residue_field(); f = k.reduction_map(); f
            Partially defined reduction map:
              From: Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 17
              To:   Residue field in a of Principal ideal (t^3 + t^2 + 7) of Univariate Polynomial Ring in t over Finite Field of size 17
            sage: f(1/t)
            12*a^2 + 12*a
        """
        return self.convert_map_from(self.p.ring().fraction_field())

    def lift_map(self):
        """
        Returns the standard map from this residue field up to the ring of
        integers lifting the canonical projection.

        EXAMPLES::

            sage: I = QQ[3^(1/3)].factor(5)[1][0]; I
            Fractional ideal (-a + 2)
            sage: k = I.residue_field(); k
            Residue field of Fractional ideal (-a + 2)
            sage: f = k.lift_map(); f
            Lifting map:
              From: Residue field of Fractional ideal (-a + 2)
              To:   Maximal Order in Number Field in a with defining polynomial x^3 - 3
            sage: f.domain()
            Residue field of Fractional ideal (-a + 2)
            sage: f.codomain()
            Maximal Order in Number Field in a with defining polynomial x^3 - 3
            sage: f(k.0)
            1

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 7)
            sage: k.<a> = P.residue_field()
            sage: f = k.lift_map(); f
            (map internal to coercion system -- copy before use)
            Lifting map:
              From: Residue field in a of Principal ideal (t^3 + t^2 + 7) of Univariate Polynomial Ring in t over Finite Field of size 17
              To:   Univariate Polynomial Ring in t over Finite Field of size 17
            sage: f(a^2 + 5)
            t^2 + 5
        """
        OK = self.p.ring()
        if OK.is_field():
            OK = OK.ring_of_integers()
        return self._internal_coerce_map_from(OK).section()

    def __cmp__(self, x):
        """
        Compares two residue fields: they are equal iff the primes
        defining them are equal and they have the same variable name.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-11)
            sage: F = K.ideal(37).factor(); F
            (Fractional ideal (37, a + 9)) * (Fractional ideal (37, a + 12)) * (Fractional ideal (2*a - 5))
            sage: k = K.residue_field(F[0][0])
            sage: l = K.residue_field(F[1][0])
            sage: k == l
            False

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 7)
            sage: k.<a> = P.residue_field()
            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 11)
            sage: l.<b> = P.residue_field()
            sage: k == l
            False
            sage: ll.<c> = P.residue_field()
            sage: ll == l
            False
        """
        c = cmp(type(self), type(x))
        if c: return c
        c = cmp(self.p, x.p)
        if c: return c
        c = cmp(self.variable_name(), x.variable_name())
        return c

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + x + 1)
            sage: hash(K.residue_field(K.prime_above(17))) # random
            -6463132282686559142
            sage: hash(K.residue_field(K.prime_above(2^60))) # random
            -6939519969600666586
            sage: R.<t> = GF(13)[]
            sage: hash(R.residue_field(t + 2)) # random
            3521289879659800254
        """
        return 1 + hash(self.ideal())

cdef class ReductionMap(Map):
    """
    A reduction map from a (subset) of a number field or function field to
    this residue class field.

    It will be defined on those elements of the field with non-negative
    valuation at the specified prime.

    EXAMPLES::

        sage: I = QQ[sqrt(17)].factor(5)[0][0]; I
        Fractional ideal (5)
        sage: k = I.residue_field(); k
        Residue field in sqrt17bar of Fractional ideal (5)
        sage: R = k.reduction_map(); R
        Partially defined reduction map:
          From: Number Field in sqrt17 with defining polynomial x^2 - 17
          To:   Residue field in sqrt17bar of Fractional ideal (5)

        sage: R.<t> = GF(next_prime(2^20))[]; P = R.ideal(t^2 + t + 1)
        sage: k = P.residue_field()
        sage: k.reduction_map()
        Partially defined reduction map:
          From: Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 1048583
          To:   Residue field in tbar of Principal ideal (t^2 + t + 1) of Univariate Polynomial Ring in t over Finite Field of size 1048583
    """
    def __init__(self, K, F, to_vs, to_order, PB, PBinv):
        """
        Create a reduction map.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3 + x^2 - 2*x + 8)
            sage: F = K.factor(2)[0][0].residue_field()
            sage: F.reduction_map()
            Partially defined reduction map:
              From: Number Field in a with defining polynomial x^3 + x^2 - 2*x + 8
              To:   Residue field of Fractional ideal (1/2*a^2 - 1/2*a + 1)

            sage: K.<theta_5> = CyclotomicField(5)
            sage: F = K.factor(7)[0][0].residue_field()
            sage: F.reduction_map()
            Partially defined reduction map:
              From: Cyclotomic Field of order 5 and degree 4
              To:   Residue field in theta_5bar of Fractional ideal (7)

            sage: R.<t> = GF(2)[]; P = R.ideal(t^7 + t^6 + t^5 + t^4 + 1)
            sage: k = P.residue_field()
            sage: k.reduction_map()
            Partially defined reduction map:
              From: Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 2 (using NTL)
              To:   Residue field in tbar of Principal ideal (t^7 + t^6 + t^5 + t^4 + 1) of Univariate Polynomial Ring in t over Finite Field of size 2 (using NTL)
            sage: type(k)
            <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_givaro_with_category'>
        """
        self._K = K
        self._F = F   # finite field
        self._to_vs = to_vs
        self._PBinv = PBinv
        self._to_order = to_order # used for lift
        self._PB = PB # used for lift
        from sage.categories.all import SetsWithPartialMaps
        self._repr_type_str = "Partially defined reduction"
        Map.__init__(self, Hom(K, F, SetsWithPartialMaps()))

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: F = K.factor(2)[0][0].residue_field()
            sage: r = F.reduction_map()
            sage: cr = copy(r) # indirect doctest
            sage: cr
            Partially defined reduction map:
              From: Number Field in a with defining polynomial x^2 + 1
              To:   Residue field of Fractional ideal (a + 1)
            sage: cr == r      # todo: comparison not implemented
            True
            sage: r(2 + a) == cr(2 + a)
            True
        """
        _slots['_K'] = self._K
        _slots['_F'] = self._F
        _slots['_to_vs'] = self._to_vs
        _slots['_PBinv'] = self._PBinv
        _slots['_to_order'] = self._to_order
        _slots['_PB'] = self._PB
        _slots['_section'] = self._section
        return Map._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: F = K.factor(2)[0][0].residue_field()
            sage: r = F.reduction_map()
            sage: cr = copy(r) # indirect doctest
            sage: cr
            Partially defined reduction map:
              From: Number Field in a with defining polynomial x^2 + 1
              To:   Residue field of Fractional ideal (a + 1)
            sage: cr == r      # todo: comparison not implemented
            True
            sage: r(2 + a) == cr(2 + a)
            True
        """
        Map._update_slots(self, _slots)
        self._K = _slots['_K']
        self._F = _slots['_F']
        self._to_vs = _slots['_to_vs']
        self._PBinv = _slots['_PBinv']
        self._to_order = _slots['_to_order']
        self._PB = _slots['_PB']
        self._section = _slots['_section']

    cpdef Element _call_(self, x):
        """
        Apply this reduction map to an element that coerces into the global
        field.

        If ``x`` doesn't map because it has negative valuation, then a
        ``ZeroDivisionError`` exception is raised.

        EXAMPLES::

            sage: K.<a> = NumberField(x^2 + 1)
            sage: F = K.factor(2)[0][0].residue_field()
            sage: r = F.reduction_map(); r
            Partially defined reduction map:
              From: Number Field in a with defining polynomial x^2 + 1
              To:   Residue field of Fractional ideal (a + 1)

        We test that calling the function also works after copying::

            sage: r = copy(r)
            sage: r(2 + a) # indirect doctest
            1
            sage: r(a/2)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot reduce field element 1/2*a modulo Fractional ideal (a + 1): it has negative valuation

            sage: R.<t> = GF(2)[]; h = t^5 + t^2 + 1
            sage: k.<a> = R.residue_field(h)
            sage: K = R.fraction_field()
            sage: f = k.convert_map_from(K)
            sage: type(f)
            <type 'sage.rings.finite_rings.residue_field.ReductionMap'>
            sage: f(1/t)
            a^4 + a
            sage: f(1/h)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in finite field.

        An example to show that the issue raised in :trac:`1951`
        has been fixed::

            sage: K.<i> = NumberField(x^2 + 1)
            sage: P1, P2 = [g[0] for g in K.factor(5)]; (P1,P2)
            (Fractional ideal (-i - 2), Fractional ideal (2*i + 1))
            sage: a = 1/(1+2*i)
            sage: F1, F2 = [g.residue_field() for g in [P1,P2]]; (F1,F2)
            (Residue field of Fractional ideal (-i - 2), Residue field of Fractional ideal (2*i + 1))
            sage: a.valuation(P1)
            0
            sage: F1(i/7)
            4
            sage: F1(a)
            3
            sage: a.valuation(P2)
            -1
            sage: F2(a)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot reduce field element -2/5*i + 1/5 modulo Fractional ideal (2*i + 1): it has negative valuation
        """
        # The reduction map is just x |--> F(to_vs(x) * (PB**(-1))) if
        # either x is integral or the denominator of x is coprime to
        # p; otherwise we work harder.
        p = self._F.p

        # Special code for residue fields of Q:
        if self._K is QQ:
            try:
                return FiniteField_prime_modn._element_constructor_(self._F, x)
            except ZeroDivisionError:
                raise ZeroDivisionError("Cannot reduce rational %s modulo %s: it has negative valuation" % (x, p.gen()))
        elif is_FractionField(self._K):
            p = p.gen()
            if p.degree() == 1:
                return self._F((x.numerator() % p)[0] / (x.denominator() % p)[0])
            else:
                return self._F((x.numerator() % p).list()) / self._F((x.denominator() % p).list())

        try:
            return self._F(self._to_vs(x) * self._PBinv)
        except Exception:
            pass

        # Now we do have to work harder...below this point we handle
        # cases which failed before trac 1951 was fixed.
        R = self._K.ring_of_integers()
        dx = R(x.denominator())
        nx = R(dx*x)
        vnx = nx.valuation(p)
        vdx = dx.valuation(p)
        if vnx > vdx:
            return self(0)
        if vnx < vdx:
            raise ZeroDivisionError("Cannot reduce field element %s modulo %s: it has negative valuation" % (x, p))

        a = self._K.uniformizer(p,'negative') ** vnx
        nx /= a
        dx /= a
        # Assertions for debugging!
        # assert nx.valuation(p) == 0 and dx.valuation(p) == 0 and x == nx/dx
        # assert nx.is_integral() and dx.is_integral()
        # print("nx = ",nx,"; dx = ",dx, ": recursing")

        # NB at this point nx and dx are in the ring of integers and
        # both are p-units.  Recursion is now safe, since integral
        # elements will not cause further recursion; and neither
        # self(nx) nor self(dx) will be 0 since nx, dx are p-units.
        return self(nx)/self(dx)

    def section(self):
        """
        Computes a section of the map, namely a map that lifts elements of the
        residue field to elements of the field.

        EXAMPLES::

            sage: K.<a> = NumberField(x^5 - 5*x + 2)
            sage: P = K.ideal(47).factor()[0][0]
            sage: k = K.residue_field(P)
            sage: f = k.convert_map_from(K)
            sage: s = f.section(); s
            Lifting map:
              From: Residue field in abar of Fractional ideal (14*a^4 - 24*a^3 - 26*a^2 + 58*a - 15)
              To:   Number Field in a with defining polynomial x^5 - 5*x + 2
            sage: s(k.gen())
            a
            sage: L.<b> = NumberField(x^5 + 17*x + 1)
            sage: P = L.factor(53)[0][0]
            sage: l = L.residue_field(P)
            sage: g = l.convert_map_from(L)
            sage: s = g.section(); s
            Lifting map:
              From: Residue field in bbar of Fractional ideal (53, b^2 + 23*b + 8)
              To:   Number Field in b with defining polynomial x^5 + 17*x + 1
            sage: s(l.gen()).parent()
            Number Field in b with defining polynomial x^5 + 17*x + 1

            sage: R.<t> = GF(2)[]; h = t^5 + t^2 + 1
            sage: k.<a> = R.residue_field(h)
            sage: K = R.fraction_field()
            sage: f = k.convert_map_from(K)
            sage: f.section()
            Lifting map:
              From: Residue field in a of Principal ideal (t^5 + t^2 + 1) of Univariate Polynomial Ring in t over Finite Field of size 2 (using NTL)
              To:   Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 2 (using NTL)
        """
        if self._section is None:
            self._section = LiftingMap(self, self._to_order, self._PB)
        return self._section


cdef class ResidueFieldHomomorphism_global(RingHomomorphism):
    """
    The class representing a homomorphism from the order of a number
    field or function field to the residue field at a given prime.

    EXAMPLES::

        sage: K.<a> = NumberField(x^3-7)
        sage: P  = K.ideal(29).factor()[0][0]
        sage: k  = K.residue_field(P)
        sage: OK = K.maximal_order()
        sage: abar = k(OK.1); abar
        abar
        sage: (1+abar)^179
        24*abar + 12

        sage: phi = k.coerce_map_from(OK); phi
        Ring morphism:
          From: Maximal Order in Number Field in a with defining polynomial x^3 - 7
          To:   Residue field in abar of Fractional ideal (2*a^2 + 3*a - 10)
        sage: phi in Hom(OK,k)
        True
        sage: phi(OK.1)
        abar

        sage: R.<t> = GF(19)[]; P = R.ideal(t^2 + 5)
        sage: k.<a> = R.residue_field(P)
        sage: f = k.coerce_map_from(R); f
        Ring morphism:
          From: Univariate Polynomial Ring in t over Finite Field of size 19
          To:   Residue field in a of Principal ideal (t^2 + 5) of Univariate Polynomial Ring in t over Finite Field of size 19
    """
    def __init__(self, K, F, to_vs, to_order, PB, PBinv):
        """
        Initialize ``self``.

        INPUT:

        - ``k`` -- The residue field that is the codomain of this morphism

        - ``p`` -- The prime ideal defining this residue field

        - ``im_gen`` -- The image of the generator of the number field

        EXAMPLES:

        We create a residue field homomorphism::

            sage: K.<theta> = CyclotomicField(5)
            sage: P = K.factor(7)[0][0]
            sage: P.residue_class_degree()
            4
            sage: kk.<a> = P.residue_field(); kk
            Residue field in a of Fractional ideal (7)
            sage: phi = kk.coerce_map_from(K.maximal_order()); phi
            Ring morphism:
              From: Maximal Order in Cyclotomic Field of order 5 and degree 4
              To:   Residue field in a of Fractional ideal (7)
            sage: type(phi)
            <type 'sage.rings.finite_rings.residue_field.ResidueFieldHomomorphism_global'>

            sage: R.<t> = GF(2)[]; P = R.ideal(t^7 + t^6 + t^5 + t^4 + 1)
            sage: k = P.residue_field(); f = k.coerce_map_from(R)
            sage: f(t^10)
            tbar^6 + tbar^3 + tbar^2
        """
        self._K = K
        self._F = F   # finite field
        self._to_vs = to_vs
        self._PBinv = PBinv
        self._PB = PB # used for lift
        self._to_order = to_order # used for lift
        self._repr_type_str = "Reduction"
        RingHomomorphism.__init__(self, Hom(K,F))

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-x+8)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k = K.residue_field(P)
            sage: OK = K.maximal_order()
            sage: phi = k.coerce_map_from(OK)
            sage: psi = copy(phi); psi    # indirect doctest
            Ring morphism:
              From: Maximal Order in Number Field in a with defining polynomial x^3 - x + 8
              To:   Residue field in abar of Fractional ideal (29)
            sage: psi == phi   # todo: comparison not implemented
            True
            sage: psi(OK.an_element()) == phi(OK.an_element())
            True
        """
        _slots['_K'] = self._K
        _slots['_F'] = self._F
        _slots['_to_vs'] = self._to_vs
        _slots['_PBinv'] = self._PBinv
        _slots['_to_order'] = self._to_order
        _slots['_PB'] = self._PB
        _slots['_section'] = self._section
        return RingHomomorphism._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-x+8)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k = K.residue_field(P)
            sage: OK = K.maximal_order()
            sage: phi = k.coerce_map_from(OK)
            sage: psi = copy(phi); psi    # indirect doctest
            Ring morphism:
              From: Maximal Order in Number Field in a with defining polynomial x^3 - x + 8
              To:   Residue field in abar of Fractional ideal (29)
            sage: psi == phi   # todo: comparison not implemented
            True
            sage: psi(OK.an_element()) == phi(OK.an_element())
            True
        """
        RingHomomorphism._update_slots(self, _slots)
        self._K = _slots['_K']
        self._F = _slots['_F']
        self._to_vs = _slots['_to_vs']
        self._PBinv = _slots['_PBinv']
        self._to_order = _slots['_to_order']
        self._PB = _slots['_PB']
        self._section = _slots['_section']

    cpdef Element _call_(self, x):
        """
        Applies this morphism to an element.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-x+8)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k =K.residue_field(P)
            sage: OK = K.maximal_order()
            sage: k.coerce_map_from(OK)(OK(a)^7) # indirect doctest
            13*abar^2 + 7*abar + 21

            sage: R.<t> = GF(next_prime(2^18))[]; P = R.ideal(t - 71)
            sage: k = ResidueField(P); f = k.coerce_map_from(R); f
            Ring morphism:
              From: Univariate Polynomial Ring in t over Finite Field of size 262147
              To:   Residue field of Principal ideal (t + 262076) of Univariate Polynomial Ring in t over Finite Field of size 262147
            sage: f(t^2)
            5041
        """
        # The reduction map is just x |--> F(to_vs(x) * (PB**(-1))) if
        # either x is integral or the denominator of x is coprime to
        # p; otherwise we work harder.

        # No special code for residue fields of Z, since we just use the normal reduction map to GF(p)
        if self._K is ZZ:
            return self._F(x)
        if is_PolynomialRing(self._K):
            p = self._F.p.gen()
            if p.degree() == 1:
                return self._F((x % p)[0])
            else:
                return self._F((x % p).list())
        return self._F(self._to_vs(x) * self._PBinv)
        #return self._F(self._to_vs(x.parent().fraction_field()(x)) * self._PBinv)

    def section(self):
        """
        Computes a section of the map, namely a map that lifts elements of
        the residue field to elements of the ring of integers.

        EXAMPLES::

            sage: K.<a> = NumberField(x^5 - 5*x + 2)
            sage: P = K.ideal(47).factor()[0][0]
            sage: k = K.residue_field(P)
            sage: f = k.coerce_map_from(K.ring_of_integers())
            sage: s = f.section(); s
            Lifting map:
              From: Residue field in abar of Fractional ideal (14*a^4 - 24*a^3 - 26*a^2 + 58*a - 15)
              To:   Maximal Order in Number Field in a with defining polynomial x^5 - 5*x + 2
            sage: s(k.gen())
            a
            sage: L.<b> = NumberField(x^5 + 17*x + 1)
            sage: P = L.factor(53)[0][0]
            sage: l = L.residue_field(P)
            sage: g = l.coerce_map_from(L.ring_of_integers())
            sage: s = g.section(); s
            Lifting map:
              From: Residue field in bbar of Fractional ideal (53, b^2 + 23*b + 8)
              To:   Maximal Order in Number Field in b with defining polynomial x^5 + 17*x + 1
            sage: s(l.gen()).parent()
            Maximal Order in Number Field in b with defining polynomial x^5 + 17*x + 1

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 7)
            sage: k.<a> = P.residue_field()
            sage: f = k.coerce_map_from(R)
            sage: f.section()
            (map internal to coercion system -- copy before use)
            Lifting map:
              From: Residue field in a of Principal ideal (t^3 + t^2 + 7) of Univariate Polynomial Ring in t over Finite Field of size 17
              To:   Univariate Polynomial Ring in t over Finite Field of size 17
        """
        if self._section is None:
            self._section = LiftingMap(self, self._to_order, self._PB)
        return self._section

    def lift(self, x):
        """
        Returns a lift of ``x`` to the Order, returning a "polynomial" in
        the generator with coefficients between 0 and `p-1`.

        EXAMPLES::

            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(29).factor()[0][0]
            sage: k = K.residue_field(P)
            sage: OK = K.maximal_order()
            sage: f = k.coerce_map_from(OK)
            sage: c = OK(a)
            sage: b = k(a)
            sage: f.lift(13*b + 5)
            13*a + 5
            sage: f.lift(12821*b+918)
            3*a + 19

            sage: R.<t> = GF(17)[]; P = R.ideal(t^3 + t^2 + 7)
            sage: k.<a> = P.residue_field(); f = k.coerce_map_from(R)
            sage: f.lift(a^2 + 5*a + 1)
            t^2 + 5*t + 1
            sage: f(f.lift(a^2 + 5*a + 1)) == a^2 + 5*a + 1
            True
        """
        if self.domain() is ZZ:
            return x.lift()
        else:
            return self.section()(x)

cdef class LiftingMap(Section):
    """
    Lifting map from residue class field to number field.

    EXAMPLES::

        sage: K.<a> = NumberField(x^3 + 2)
        sage: F = K.factor(5)[0][0].residue_field()
        sage: F.degree()
        2
        sage: L = F.lift_map(); L
        Lifting map:
          From: Residue field in abar of Fractional ideal (a^2 + 2*a - 1)
          To:   Maximal Order in Number Field in a with defining polynomial x^3 + 2
        sage: L(F.0^2)
        3*a + 1
        sage: L(3*a + 1) == F.0^2
        True

        sage: R.<t> = GF(13)[]
        sage: P = R.ideal(8*t^12 + 9*t^11 + 11*t^10 + 2*t^9 + 11*t^8 + 3*t^7 + 12*t^6 + t^4 + 7*t^3 + 5*t^2 + 12*t + 1)
        sage: k.<a> = P.residue_field()
        sage: k.lift_map()
        Lifting map:
          From: Residue field in a of Principal ideal (t^12 + 6*t^11 + 3*t^10 + 10*t^9 + 3*t^8 + 2*t^7 + 8*t^6 + 5*t^4 + 9*t^3 + 12*t^2 + 8*t + 5) of Univariate Polynomial Ring in t over Finite Field of size 13
          To:   Univariate Polynomial Ring in t over Finite Field of size 13
    """
    def __init__(self, reduction, to_order, PB):
        """
        Create a lifting map.

        EXAMPLES::

            sage: K.<theta_5> = CyclotomicField(5)
            sage: F = K.factor(7)[0][0].residue_field()
            sage: F.lift_map()
            Lifting map:
              From: Residue field in theta_5bar of Fractional ideal (7)
              To:   Maximal Order in Cyclotomic Field of order 5 and degree 4

            sage: K.<a> = NumberField(x^5 + 2)
            sage: F = K.factor(7)[0][0].residue_field()
            sage: L = F.lift_map(); L
            Lifting map:
              From: Residue field in abar of Fractional ideal (-2*a^4 + a^3 - 4*a^2 + 2*a - 1)
              To:   Maximal Order in Number Field in a with defining polynomial x^5 + 2
            sage: L.domain()
            Residue field in abar of Fractional ideal (-2*a^4 + a^3 - 4*a^2 + 2*a - 1)

            sage: K.<a> = CyclotomicField(7)
            sage: F = K.factor(5)[0][0].residue_field()
            sage: L = F.lift_map(); L
            Lifting map:
              From: Residue field in abar of Fractional ideal (5)
              To:   Maximal Order in Cyclotomic Field of order 7 and degree 6
            sage: L.codomain()
            Maximal Order in Cyclotomic Field of order 7 and degree 6

            sage: R.<t> = GF(2)[]; h = t^5 + t^2 + 1
            sage: k.<a> = R.residue_field(h)
            sage: K = R.fraction_field()
            sage: L = k.lift_map(); L.codomain()
            Univariate Polynomial Ring in t over Finite Field of size 2 (using NTL)
        """
        self._K = reduction._K
        self._F = reduction._F   # finite field
        self._to_order = to_order
        self._PB = PB
        Section.__init__(self, reduction)

    cdef dict _extra_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(7)
            sage: F = K.factor(5)[0][0].residue_field()
            sage: phi = F.lift_map()
            sage: psi = copy(phi); psi   # indirect doctest
            Lifting map:
              From: Residue field in abar of Fractional ideal (5)
              To:   Maximal Order in Cyclotomic Field of order 7 and degree 6
            sage: psi == phi             # todo: comparison not implemented
            False
            sage: phi(F.0) == psi(F.0)
            True
        """
        _slots['_K'] = self._K
        _slots['_F'] = self._F
        _slots['_to_order'] = self._to_order
        _slots['_PB'] = self._PB
        return Section._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        Helper for copying and pickling.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(7)
            sage: F = K.factor(5)[0][0].residue_field()
            sage: phi = F.lift_map()
            sage: psi = copy(phi); psi   # indirect doctest
            Lifting map:
              From: Residue field in abar of Fractional ideal (5)
              To:   Maximal Order in Cyclotomic Field of order 7 and degree 6
            sage: psi == phi             # todo: comparison not implemented
            False
            sage: phi(F.0) == psi(F.0)
            True
        """
        Section._update_slots(self, _slots)
        self._K = _slots['_K']
        self._F = _slots['_F']
        self._to_order = _slots['_to_order']
        self._PB = _slots['_PB']

    cpdef Element _call_(self, x):
        """
        Lift from this residue class field to the number field.

        EXAMPLES::

            sage: K.<a> = CyclotomicField(7)
            sage: F = K.factor(5)[0][0].residue_field()
            sage: L = F.lift_map(); L
            Lifting map:
              From: Residue field in abar of Fractional ideal (5)
              To:   Maximal Order in Cyclotomic Field of order 7 and degree 6
            sage: L(F.0) # indirect doctest
            a
            sage: F(a)
            abar

            sage: R.<t> = GF(2)[]; h = t^5 + t^2 + 1
            sage: k.<a> = R.residue_field(h)
            sage: K = R.fraction_field()
            sage: f = k.lift_map()
            sage: f(a^2)
            t^2
            sage: f(a^6)
            t^3 + t
        """
        if self._K is QQ or self._K is ZZ:
            return self._K(x.lift())  # x.lift() is in ZZ
        elif is_FractionField(self._K):
            if self._F.p.degree() == 1:
                return self._K(self._K.ring_of_integers()(x))
            else:
                return self._K(self._K.ring_of_integers()(x.polynomial().list()))
        elif is_PolynomialRing(self._K):
            return self._K(x.polynomial().list())
        # Else the lifting map is just x |--> to_order(x * PB)
        x = self._F(x)
        v = x.polynomial().padded_list(self._F.degree())
        ans = self._to_order(self._PB.linear_combination_of_rows(v))
        if ans.parent() is self._K:
            return ans
        else:
            return self._K(ans)

    def _repr_type(self):
        """
        EXAMPLES::

            sage: K.<theta_12> = CyclotomicField(12)
            sage: F.<tmod> = K.factor(7)[0][0].residue_field()
            sage: F.lift_map() #indirect doctest
            Lifting map:
              From: Residue field in tmod of Fractional ideal (-3*theta_12^2 + 1)
              To:   Maximal Order in Cyclotomic Field of order 12 and degree 4
        """
        return "Lifting"

class ResidueFiniteField_prime_modn(ResidueField_generic, FiniteField_prime_modn):
    """
    The class representing residue fields of number fields that have
    prime order.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^3-7)
        sage: P = K.ideal(29).factor()[1][0]
        sage: k = ResidueField(P)
        sage: k
        Residue field of Fractional ideal (a^2 + 2*a + 2)
        sage: k.order()
        29
        sage: OK = K.maximal_order()
        sage: c = OK(a)
        sage: b = k(a)
        sage: k.coerce_map_from(OK)(c)
        16
        sage: k(4)
        4
        sage: k(c + 5)
        21
        sage: b + c
        3

        sage: R.<t> = GF(7)[]; P = R.ideal(2*t + 3)
        sage: k = P.residue_field(); k
        Residue field of Principal ideal (t + 5) of Univariate Polynomial Ring in t over Finite Field of size 7
        sage: k(t^2)
        4
        sage: k.order()
        7
    """
    def __init__(self, p, name, intp, to_vs, to_order, PB):
        """
        Initialize ``self``.

        INPUT:

        - ``p`` -- A prime ideal of a number field

        - ``name`` -- the name of the generator of this extension

        - ``intp`` -- the rational prime that ``p`` lies over

        EXAMPLES::

            sage: K.<i> = QuadraticField(-1)
            sage: kk = ResidueField(K.factor(5)[0][0])
            sage: type(kk)
            <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_prime_modn_with_category'>

            sage: R.<t> = GF(7)[]; P = R.ideal(2*t + 3)
            sage: k = P.residue_field(); type(k)
            <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_prime_modn_with_category'>
        """
        ResidueField_generic.__init__(self, p)
        FiniteField_prime_modn.__init__(self, intp)
        from sage.rings.finite_rings.integer_mod import IntegerMod_to_IntegerMod, Integer_to_IntegerMod, Int_to_IntegerMod
        K = OK = p.ring()
        if OK.is_field():
            OK = OK.ring_of_integers()
        else:
            K = K.fraction_field()
        if PB is None:
            if OK is ZZ:
                # integer case
                coerce_list = [IntegerMod_to_IntegerMod(GF(intp), self), Integer_to_IntegerMod(self), Int_to_IntegerMod(self)]
            else:
                # polynomial ring case.
                coerce_list = [ResidueFieldHomomorphism_global(OK, self, None, None, None, None), OK.base_ring()]
            self._populate_coercion_lists_(coerce_list=coerce_list, convert_list=[ReductionMap(K, self, None, None, None, None)]) # could be special-cased a bit more.
        else:
            PBinv = PB**(-1)
            self._populate_coercion_lists_(coerce_list=[IntegerMod_to_IntegerMod(GF(intp), self), Integer_to_IntegerMod(self), Int_to_IntegerMod(self), ResidueFieldHomomorphism_global(OK, self, to_vs, to_order, PB, PBinv)], \
                                                 convert_list=[ReductionMap(K, self, to_vs, to_order, PB, PBinv)])

    def _element_constructor_(self, x):
        """
        Construct and/or coerce ``x`` into an element of ``self``.

        INPUT:

           - ``x`` -- something to cast in to ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(29).factor()[1][0]
            sage: k = ResidueField(P)
            sage: k
            Residue field of Fractional ideal (a^2 + 2*a + 2)
            sage: OK = K.maximal_order()
            sage: c = OK(a)
            sage: b = k(a); b
            16
            sage: k(2r)
            2
            sage: V = k.vector_space(); v = V([3])
            sage: type(k.convert_map_from(V))
            <type 'sage.structure.coerce_maps.DefaultConvertMap_unique'>
            sage: k(v) # indirect doctest
            3

            sage: R.<t> = GF(2)[]; P = R.ideal(t+1); k.<a> = P.residue_field()
            sage: V = k.vector_space(); v = V([1])
            sage: k(v)
            1
        """
        if isinstance(x, FreeModuleElement) and len(x) == 1:
            x = x[0]
        try:
            return FiniteField_prime_modn._element_constructor_(self, x)
        except TypeError:
            return ResidueField_generic._element_constructor_(self, x)

class ResidueFiniteField_pari_ffelt(ResidueField_generic, FiniteField_pari_ffelt):
    """
    The class representing residue fields of number fields that have non-prime
    order at least `2^16`.

    EXAMPLES::

        sage: K.<a> = NumberField(x^3-7)
        sage: P = K.ideal(923478923).factor()[0][0]
        sage: k = K.residue_field(P)
        sage: k.degree()
        2
        sage: OK = K.maximal_order()
        sage: c = OK(a)
        sage: b = k(c)
        sage: b+c
        2*abar
        sage: b*c
        664346875*abar + 535606347
        sage: k.base_ring()
        Finite Field of size 923478923

        sage: R.<t> = GF(5)[]; P = R.ideal(4*t^12 + 3*t^11 + 4*t^10 + t^9 + t^8 + 3*t^7 + 2*t^6 + 3*t^4 + t^3 + 3*t^2 + 2)
        sage: k.<a> = P.residue_field()
        sage: type(k)
        <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_pari_ffelt_with_category'>
        sage: k(1/t)
        3*a^11 + a^10 + 3*a^9 + 2*a^8 + 2*a^7 + a^6 + 4*a^5 + a^3 + 2*a^2 + a
    """
    def __init__(self, p, characteristic, name, modulus, to_vs, to_order, PB):
        """
        Initialize ``self``.

        EXAMPLES::

        We create a residue field with implementation ``pari_ffelt``::

            sage: K.<a> = NumberField(x^3-7)
            sage: P = K.ideal(923478923).factor()[0][0]
            sage: type(P.residue_field())
            <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_pari_ffelt_with_category'>
        """
        ResidueField_generic.__init__(self, p)
        FiniteField_pari_ffelt.__init__(self, characteristic, modulus, name)
        K = OK = p.ring()
        if OK.is_field():
            OK = OK.ring_of_integers()
        else:
            K = K.fraction_field()
        if PB is None:
            PBinv = None
        else:
            PBinv = PB**(-1)
        self._populate_coercion_lists_(coerce_list=[self.base_ring(), ResidueFieldHomomorphism_global(OK, self, to_vs, to_order, PB, PBinv)], convert_list = [ReductionMap(K, self, to_vs, to_order, PB, PBinv)])

    def _element_constructor_(self, x):
        """
        Coerce ``x`` into ``self``.

        EXAMPLES::

            sage: K.<aa> = NumberField(x^3 - 2)
            sage: P = K.factor(10007)[0][0]
            sage: P.residue_class_degree()
            2
            sage: ff.<alpha> = P.residue_field(); ff
            Residue field in alpha of Fractional ideal (-12*aa^2 + 189*aa - 475)
            sage: type(ff)
            <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_pari_ffelt_with_category'>
            sage: ff(alpha^2 + 1)
            7521*alpha + 4131
            sage: ff(17/3)
            6677
            sage: V = ff.vector_space(); v = V([3,-2])
            sage: type(ff.convert_map_from(V))
            <type 'sage.structure.coerce_maps.DefaultConvertMap_unique'>
            sage: ff(v) # indirect doctest
            10005*alpha + 3

            sage: R.<t> = GF(5)[]; P = R.ideal(4*t^12 + 3*t^11 + 4*t^10 + t^9 + t^8 + 3*t^7 + 2*t^6 + 3*t^4 + t^3 + 3*t^2 + 2)
            sage: k.<a> = P.residue_field()
            sage: V = k.vector_space(); v = V([1,2,3,4,5,6,7,8,9,0,1,2]); k(v) # indirect doctest
            2*a^11 + a^10 + 4*a^8 + 3*a^7 + 2*a^6 + a^5 + 4*a^3 + 3*a^2 + 2*a + 1
        """
        try:
            return FiniteField_pari_ffelt._element_constructor_(self, x)
        except TypeError:
            return ResidueField_generic._element_constructor_(self, x)

class ResidueFiniteField_givaro(ResidueField_generic, FiniteField_givaro):
    """
    The class representing residue fields of number fields that have non-prime
    order strictly less than `2^16`.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^3-7)
        sage: P = K.ideal(29).factor()[0][0]
        sage: k =K.residue_field(P)
        sage: k.degree()
        2
        sage: OK = K.maximal_order()
        sage: c = OK(a)
        sage: b = k(c)
        sage: b*c^2
        7
        sage: b*c
        13*abar + 5

        sage: R.<t> = GF(7)[]; P = R.ideal(t^2 + 4)
        sage: k.<a> = R.residue_field(P); type(k)
        <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_givaro_with_category'>
        sage: k(1/t)
        5*a
    """
    def __init__(self, p, q, name, modulus, to_vs, to_order, PB):
        r"""
        INPUT:

        - ``p`` -- the prime ideal defining this residue field

        - ``q`` -- the order of this residue field (a power of intp)

        - ``name`` -- the name of the generator of this extension

        - ``modulus`` -- the polynomial modulus for this extension

        - ``to_vs`` -- the map from the number field (or function field) to
          the appropriate vector space (over `\QQ` or `F_p(t)`)

        - ``to_order`` -- the map from a lattice in that vector space to the maximal order

        - ``PB`` -- a matrix used in defining the reduction and lifting maps.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: k = K.residue_field(P)

            sage: R.<t> = GF(3)[]; P = R.ideal(t^4 - t^3 + t + 1); k.<a> = P.residue_field(); type(k)
            <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_givaro_with_category'>
            sage: a^5
            a^3 + 2*a^2 + a + 2
        """
        ResidueField_generic.__init__(self, p)
        FiniteField_givaro.__init__(self, q, name, modulus)
        K = OK = p.ring()
        if OK.is_field():
            OK = OK.ring_of_integers()
        else:
            K = K.fraction_field()
        if PB is None:
            PBinv = None
        else:
            PBinv = PB**(-1)
        self._populate_coercion_lists_(coerce_list=[self.base_ring(), ResidueFieldHomomorphism_global(OK, self, to_vs, to_order, PB, PBinv)], convert_list = [ReductionMap(K, self, to_vs, to_order, PB, PBinv)])

    def _element_constructor_(self, x):
        """
        INPUT:

            - ``x`` -- Something to cast into ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: k =K.residue_field(P)
            sage: k(77*a^7+4)
            2*abar + 4
            sage: V = k.vector_space(); v = V([3,-2])
            sage: type(k.convert_map_from(V))
            <type 'sage.structure.coerce_maps.DefaultConvertMap_unique'>
            sage: k(v) # indirect doctest
            59*abar + 3

            sage: R.<t> = GF(3)[]; P = R.ideal(t^4 - t^3 + t + 1); k.<a> = P.residue_field()
            sage: V = k.vector_space(); v = V([0,1,2,3])
            sage: k(v) # indirect doctest
            2*a^2 + a
        """
        try:
            return FiniteField_givaro._element_constructor_(self, x)
        except TypeError:
            return ResidueField_generic._element_constructor_(self, x)

class ResidueFiniteField_ntl_gf2e(ResidueField_generic, FiniteField_ntl_gf2e):
    """
    The class representing residue fields with order a power of 2.

    When the order is less than `2^16`, givaro is used by default instead.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: K.<a> = NumberField(x^3-7)
        sage: P = K.ideal(29).factor()[0][0]
        sage: k =K.residue_field(P)
        sage: k.degree()
        2
        sage: OK = K.maximal_order()
        sage: c = OK(a)
        sage: b = k(c)
        sage: b*c^2
        7
        sage: b*c
        13*abar + 5

        sage: R.<t> = GF(2)[]; P = R.ideal(t^19 + t^5 + t^2 + t + 1)
        sage: k.<a> = R.residue_field(P); type(k)
        <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_ntl_gf2e_with_category'>
        sage: k(1/t)
        a^18 + a^4 + a + 1
        sage: k(1/t)*t
        1
    """
    # we change the order for consistency with FiniteField_ntl_gf2e's __cinit__
    def __init__(self, q, name, modulus, repr, p, to_vs, to_order, PB):
        """
        INPUT:

        - ``p`` -- the prime ideal defining this residue field

        - ``q`` -- the order of this residue field

        - ``name`` -- the name of the generator of this extension

        - ``modulus`` -- the polynomial modulus for this extension

        - ``to_vs`` -- the map from the number field (or function field) to
          the appropriate vector space (over `\QQ` or `F_p(t)`)

        - ``to_order`` -- the map from a lattice in that vector space to the
          maximal order

        - ``PB`` -- a matrix used in defining the reduction and lifting maps

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: k = K.residue_field(P)

            sage: R.<t> = GF(3)[]; P = R.ideal(t^4 - t^3 + t + 1); k.<a> = P.residue_field(); type(k)
            <class 'sage.rings.finite_rings.residue_field.ResidueFiniteField_givaro_with_category'>
            sage: a^5
            a^3 + 2*a^2 + a + 2
        """
        ResidueField_generic.__init__(self, p)
        FiniteField_ntl_gf2e.__init__(self, q, name, modulus, repr)
        K = OK = p.ring()
        if OK.is_field():
            OK = OK.ring_of_integers()
        else:
            K = K.fraction_field()
        if PB is None:
            PBinv = None
        else:
            PBinv = PB**(-1)
        self._populate_coercion_lists_(coerce_list=[self.base_ring(), ResidueFieldHomomorphism_global(OK, self, to_vs, to_order, PB, PBinv)], convert_list = [ReductionMap(K, self, to_vs, to_order, PB, PBinv)])

    def _element_constructor_(self, x):
        """
        INPUT:

        - ``x`` -- Something to cast into ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K.<a> = NumberField(x^4+3*x^2-17)
            sage: P = K.ideal(61).factor()[0][0]
            sage: k =K.residue_field(P)
            sage: k(77*a^7+4)
            2*abar + 4
            sage: V = k.vector_space(); v = V([3,-2])
            sage: type(k.convert_map_from(V))
            <type 'sage.structure.coerce_maps.DefaultConvertMap_unique'>
            sage: k(v) # indirect doctest
            59*abar + 3

            sage: R.<t> = GF(3)[]; P = R.ideal(t^4 - t^3 + t + 1); k.<a> = P.residue_field()
            sage: V = k.vector_space(); v = V([0,1,2,3])
            sage: k(v) # indirect doctest
            2*a^2 + a
        """
        try:
            return FiniteField_ntl_gf2e._element_constructor_(self, x)
        except TypeError:
            return ResidueField_generic._element_constructor_(self, x)


