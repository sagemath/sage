# sage.doctest: needs sage.libs.pari
r"""
Elements of Quotients of Univariate Polynomial Rings

EXAMPLES: We create a quotient of a univariate polynomial ring over
`\ZZ`.

::

    sage: R.<x> = ZZ[]
    sage: S.<a> = R.quotient(x^3 + 3*x - 1)
    sage: 2 * a^3
    -6*a + 2

Next we make a univariate polynomial ring over
`\ZZ[x]/(x^3+3x-1)`.

::

    sage: S1.<y> = S[]

And, we quotient out that by `y^2 + a`.

::

    sage: T.<z> = S1.quotient(y^2 + a)

In the quotient `z^2` is `-a`.

::

    sage: z^2
    -a

And since `a^3 = -3x + 1`, we have::

    sage: z^6
    3*a - 1

::

    sage: R.<x> = PolynomialRing(Integers(9))
    sage: S.<a> = R.quotient(x^4 + 2*x^3 + x + 2)
    sage: a^100
    7*a^3 + 8*a + 7

::

    sage: R.<x> = PolynomialRing(QQ)
    sage: S.<a> = R.quotient(x^3 - 2)
    sage: a
    a
    sage: a^3
    2

For the purposes of comparison in Sage the quotient element
`a^3` is equal to `x^3`. This is because when the
comparison is performed, the right element is coerced into the
parent of the left element, and `x^3` coerces to
`a^3`.

::

    sage: a == x
    True
    sage: a^3 == x^3
    True
    sage: x^3
    x^3
    sage: S(x^3)
    2

AUTHORS:

- William Stein
"""

#*****************************************************************************
#       Copyright (C) 2005, 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import CommutativeRingElement
from sage.structure.richcmp import richcmp
import sage.rings.polynomial.polynomial_singular_interface as polynomial_singular_interface


class PolynomialQuotientRingElement(polynomial_singular_interface.Polynomial_singular_repr, CommutativeRingElement):
    """
    Element of a quotient of a polynomial ring.

    EXAMPLES::

        sage: P.<x> = QQ[]
        sage: Q.<xi> = P.quo([(x^2 + 1)])
        sage: xi^2
        -1
        sage: singular(xi)                                                              # needs sage.libs.singular
        xi
        sage: (singular(xi)*singular(xi)).NF('std(0)')                                  # needs sage.libs.singular
        -1
    """
    def __init__(self, parent, polynomial, check=True):
        """
        Create an element of the quotient of a polynomial ring.

        INPUT:

        - ``parent`` -- a quotient of a polynomial ring

        - ``polynomial`` -- a polynomial

        - ``check`` -- boolean (default: ``True``); whether or not to
          verify that x is a valid element of the polynomial ring and reduced
          (mod the modulus).
        """
        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
        from sage.rings.polynomial.polynomial_element import Polynomial

        CommutativeRingElement.__init__(self, parent)
        if check:
            if not isinstance(parent, PolynomialQuotientRing_generic):
                raise TypeError("parent must be a polynomial quotient ring")

            if not isinstance(polynomial, Polynomial):
                raise TypeError("polynomial must be a polynomial")

            if polynomial not in parent.polynomial_ring():
                raise TypeError("polynomial must be in the polynomial ring of the parent")

        f = parent.modulus()
        if polynomial.degree() >= f.degree() and polynomial.degree() >= 0:
            try:
                polynomial %= f
            except AttributeError:
                A = polynomial
                B = f
                R = A
                P = B.parent()
                Q = P(0)
                X = P.gen()
                while R.degree() >= B.degree():
                    S = P(R.leading_coefficient()/B.leading_coefficient()) * X**(R.degree()-B.degree())
                    Q = Q + S
                    R = R - S*B
                polynomial = R
        self._polynomial = polynomial

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of this element under the morphism defined by
        ``im_gens`` in ``codomain``, where elements of the
        base ring are mapped by ``base_map``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: Zx.<x> = ZZ[]
            sage: K.<i> = NumberField(x^2 + 1)
            sage: cc = K.hom([-i])
            sage: S.<y> = K[]
            sage: Q.<q> = S.quotient(y^2*(y-1)*(y-i))
            sage: T.<t> = S.quotient(y*(y+1))
            sage: phi = Q.hom([t+1], base_map=cc)
            sage: phi(q)
            t + 1
            sage: phi(i*q)
            -i*t - i
        """
        return self._polynomial._im_gens_(codomain, im_gens, base_map=base_map)

    def __hash__(self):
        return hash(self._polynomial)

    def __reduce__(self):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S.<a> = R.quotient(2*x^3 + 3/2*x -1/3)
            sage: 2 * a^3
            -3/2*a + 1/3
            sage: loads(dumps(2*a^3)) == 2*a^3
            True
        """
        return self.__class__, (self.parent(), self._polynomial, False)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S.<a> = R.quotient(3*x^3 + 3/2*x -1/3)
            sage: 3 * a^3 + S.modulus()
            -3/2*a + 1/3
        """
        # We call _repr since _repr_ does not have a name variable.
        # This is very fragile!
        return self._polynomial._repr(self.parent().variable_name())

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S.<a> = R.quotient(3*x^3 + 3/2*x -1/3)
            sage: latex(a*(3 * a^3) + S.modulus())
            -\frac{3}{2} a^{2} + \frac{1}{3} a
        """
        return self._polynomial._latex_(self.parent().variable_name())

    def __pari__(self):
        """
        Pari representation of this quotient element.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: I = R.ideal(x^10)
            sage: S.<xb> = R.quo(I)
            sage: pari(xb)^10
            Mod(0, x^10)
        """
        return self._polynomial.__pari__().Mod(self.parent().modulus())

    ##################################################
    # Arithmetic
    ##################################################

    def _mul_(self, right):
        """
        Return the product of two polynomial ring quotient elements.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3-2)
            sage: (a^2 - 4) * (a+2)
            2*a^2 - 4*a - 6
        """
        R = self.parent()
        prod = self._polynomial * right._polynomial
        return self.__class__(R, prod, check=False)

    def _sub_(self, right):
        """
        Return the difference of two polynomial ring quotient elements.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 - 2)
            sage: (a^2 - 4) - (a+2)
            a^2 - a - 6
            sage: int(1) - a
            -a + 1
        """
        return self.__class__(self.parent(),
                                             self._polynomial - right._polynomial, check=False)

    def _add_(self, right):
        """
        Return the sum of two polynomial ring quotient elements.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 - 2)
            sage: (a^2 - 4) + (a+2)
            a^2 + a - 2
            sage: int(1) + a
            a + 1
        """
        return self.__class__(self.parent(),
                                             self._polynomial + right._polynomial, check=False)

    def _div_(self, right):
        """
        Return the quotient of two polynomial ring quotient elements.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 - 2)
            sage: (a^2 - 4) / (a+2)
            a - 2
        """
        return self * ~right

    def __neg__(self):
        return self.__class__(self.parent(), -self._polynomial)

    def _richcmp_(self, other, op):
        """
        Compare this element with something else, where equality testing
        coerces the object on the right, if possible (and necessary).

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 - 2)
            sage: (a^2 - 4) / (a+2) == a - 2
            True
            sage: a^2 - 4 == a
            False
        """
        return richcmp(self._polynomial, other._polynomial, op)

    def __getitem__(self, n):
        return self._polynomial[n]

    def __int__(self):
        """
        Coerce this element to an int if possible.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 - 2)
            sage: int(S(10))
            10
            sage: int(a)
            Traceback (most recent call last):
            ...
            TypeError: cannot convert nonconstant polynomial
        """
        return int(self._polynomial)

    def is_unit(self):
        """
        Return ``True`` if ``self`` is invertible.

        .. WARNING::

            Only implemented when the base ring is a field.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S.<y> = R.quotient(x^2 + 2*x + 1)
            sage: (2*y).is_unit()
            True
            sage: (y + 1).is_unit()
            False

        TESTS:

        Raise an exception if the base ring is not a field
        (see :issue:`13303`)::

            sage: Z16x.<x> = Integers(16)[]
            sage: S.<y> =  Z16x.quotient(x^2 + x + 1)
            sage: (2*y).is_unit()
            Traceback (most recent call last):
            ...
            NotImplementedError: The base ring (=Ring of integers modulo 16) is not a field

        Check that :issue:`29469` is fixed::

            sage: S(3).is_unit()
            True
        """
        if self._polynomial.is_zero():
            return False
        if self._polynomial.is_one():
            return True
        try:
            if self._polynomial.is_unit():
                return True
        except NotImplementedError:
            pass

        parent = self.parent()
        base = parent.base_ring()
        if not base.is_field():
            raise NotImplementedError("The base ring (=%s) is not a field" % base)
        g = parent.modulus().gcd(self._polynomial)
        return g.degree() == 0

    def __invert__(self):
        """
        Return the inverse of this element.

        .. WARNING::

            Only implemented when the base ring is a field.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S.<y> = R.quotient(x^2 + 2*x + 1)
            sage: (2*y)^(-1)
            -1/2*y - 1

        Raises a :exc:`ZeroDivisionError` if this element is not a unit::

            sage: (y+1)^(-1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: element y + 1 of quotient polynomial ring not invertible

        TESTS:

        An element is not invertible if the base ring is not a field
        (see :issue:`13303`)::

            sage: Z16x.<x> = Integers(16)[]
            sage: S.<y> =  Z16x.quotient(x^2 + x + 1)
            sage: (2*y)^(-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: The base ring (=Ring of integers modulo 16) is not a field

        Check that :issue:`29469` is fixed::

            sage: ~S(3)
            11
        """
        if self._polynomial.is_zero():
            raise ZeroDivisionError("element %s of quotient polynomial ring not invertible" % self)
        if self._polynomial.is_one():
            return self

        parent = self.parent()

        try:
            if self._polynomial.is_unit():
                inv_pol = self._polynomial.inverse_of_unit()
                return parent(inv_pol)
        except (TypeError, NotImplementedError):
            pass

        base = parent.base_ring()
        if not base.is_field():
            raise NotImplementedError("The base ring (=%s) is not a field" % base)
        g, _, a = parent.modulus().xgcd(self._polynomial)
        if g.degree() != 0:
            raise ZeroDivisionError("element %s of quotient polynomial ring not invertible" % self)
        c = g[0]
        return self.__class__(self.parent(), (~c)*a, check=False)

    def field_extension(self, names):
        r"""
        Given a polynomial with base ring a quotient ring, return a
        3-tuple: a number field defined by the same polynomial, a
        homomorphism from its parent to the number field sending the
        generators to one another, and the inverse isomorphism.

        INPUT:

        - ``names`` -- name of generator of output field

        OUTPUT:

        -  field

        -  homomorphism from ``self`` to field

        -  homomorphism from field to ``self``

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<alpha> = R.quotient(x^3 - 2)
            sage: F.<a>, f, g = alpha.field_extension()
            sage: F
            Number Field in a with defining polynomial x^3 - 2
            sage: a = F.gen()
            sage: f(alpha)
            a
            sage: g(a)
            alpha

        Over a finite field, the corresponding field extension is not a
        number field::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = GF(25,'b')['x']
            sage: S.<a> = R.quo(x^3 + 2*x + 1)
            sage: F.<b>, g, h = a.field_extension()
            sage: h(b^2 + 3)
            a^2 + 3
            sage: g(x^2 + 2)
            b^2 + 2

        We do an example involving a relative number field::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ['x']
            sage: K.<a> = NumberField(x^3 - 2)
            sage: S.<X> = K['X']
            sage: Q.<b> = S.quo(X^3 + 2*X + 1)
            sage: F, g, h = b.field_extension('c')

        Another more awkward example::

            sage: # needs sage.rings.number_field
            sage: R.<x> = QQ['x']
            sage: K.<a> = NumberField(x^3 - 2)
            sage: S.<X> = K['X']
            sage: f = (X+a)^3 + 2*(X+a) + 1
            sage: f
            X^3 + 3*a*X^2 + (3*a^2 + 2)*X + 2*a + 3
            sage: Q.<z> = S.quo(f)
            sage: F.<w>, g, h = z.field_extension()
            sage: c = g(z)
            sage: f(c)
            0
            sage: h(g(z))
            z
            sage: g(h(w))
            w

        AUTHORS:

        - Craig Citro (2006-08-06)

        - William Stein (2006-08-06)
        """
        #TODO: is the return order backwards from the magma convention?

##         We do another example over $\ZZ$::
##
##             sage: R.<x> = ZZ['x']
##             sage: S.<a> = R.quo(x^3 - 2)
##             sage: F.<b>, g, h = a.field_extension()
##             sage: h(b^2 + 3)
##             a^2 + 3
##             sage: g(x^2 + 2)
##             a^2 + 2
##
##         Note that the homomorphism is not defined on the entire
##         ''domain''.   (Allowing creation of such functions may be
##         disallowed in a future version of Sage.)::        <----- INDEED!
##
##             sage: h(1/3)
##             Traceback (most recent call last):
##             ...
##             TypeError: Unable to coerce rational (=1/3) to an Integer.
##
##         Note that the parent ring must be an integral domain::
##
##             sage: R.<x> = GF(25,'b')['x']
##             sage: S.<a> = R.quo(x^3 - 2)
##             sage: F, g, h = a.field_extension()
##             Traceback (most recent call last):
##             ...
##             ValueError: polynomial must be irreducible

        R = self.parent()
        x = R.gen()

        F = R.modulus().root_field(names)
        alpha = F.gen()

        f = R.hom([alpha], F, check=False)

        from sage.rings.number_field.number_field_rel import NumberField_relative
        if isinstance(F, NumberField_relative):

            base_map = F.base_field().hom([R.base_ring().gen()])
            g = F.Hom(R)(x, base_map)

        else:
            g = F.hom([x], R, check=False)

        return F, f, g

    def charpoly(self, var):
        """
        The characteristic polynomial of this element, which is by
        definition the characteristic polynomial of right multiplication by
        this element.

        INPUT:

        - ``var`` -- string; the variable name

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quo(x^3 -389*x^2 + 2*x - 5)
            sage: a.charpoly('X')                                                       # needs sage.modules
            X^3 - 389*X^2 + 2*X - 5
        """
        return self.matrix().charpoly(var)

    def fcp(self, var='x'):
        """
        Return the factorization of the characteristic polynomial of this
        element.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 -389*x^2 + 2*x - 5)
            sage: a.fcp('x')                                                            # needs sage.modules
            x^3 - 389*x^2 + 2*x - 5
            sage: S(1).fcp('y')                                                         # needs sage.modules
            (y - 1)^3
        """
        return self.charpoly(var).factor()

    def lift(self):
        """
        Return lift of this polynomial quotient ring element to the unique
        equivalent polynomial of degree less than the modulus.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 - 2)
            sage: b = a^2 - 3
            sage: b
            a^2 - 3
            sage: b.lift()
            x^2 - 3
        """
        return self._polynomial

    def __iter__(self):
        return iter(self.list())

    def list(self, copy=True):
        """
        Return list of the elements of ``self``, of length the same as the
        degree of the quotient polynomial ring.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 + 2*x - 5)
            sage: a^10
            -134*a^2 - 35*a + 300
            sage: (a^10).list()
            [300, -35, -134]
        """
        v = self._polynomial.list(copy=False)
        R = self.parent()
        n = R.degree()
        return v + [R.base_ring()(0)]*(n - len(v))

    def matrix(self):
        """
        The matrix of right multiplication by this element on the power
        basis for the quotient ring.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 + 2*x - 5)
            sage: a.matrix()                                                            # needs sage.modules
            [ 0  1  0]
            [ 0  0  1]
            [ 5 -2  0]
        """
        # Multiply each power of field generator on the right by this
        # element, then return the matrix whose rows are the
        # coefficients of the result.
        try:
            return self.__matrix
        except AttributeError:
            R = self.parent()
            v = []
            x = R.gen()
            a = R(1)
            d = R.degree()
            for _ in range(d):
                v += (a*self).list()
                a *= x
            S = R.base_ring()
            import sage.matrix.matrix_space
            M = sage.matrix.matrix_space.MatrixSpace(S, d)
            self.__matrix = M(v)
            return self.__matrix

    def minpoly(self):
        """
        The minimal polynomial of this element, which is by definition the
        minimal polynomial of the :meth:`matrix` of this element.

        ALGORITHM: Use
        :meth:`~sage.rings.polynomial.polynomial_zz_pex.Polynomial_ZZ_pEX.minpoly_mod`
        if possible, otherwise compute the minimal polynomial of the :meth:`matrix`.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 + 2*x - 5)
            sage: (a + 123).minpoly()                                                   # needs sage.modules
            x^3 - 369*x^2 + 45389*x - 1861118
            sage: (a + 123).matrix().minpoly()                                          # needs sage.modules
            x^3 - 369*x^2 + 45389*x - 1861118

        One useful application of this function is to compute a minimal
        polynomial of a finite-field element over an intermediate extension,
        rather than the absolute minimal polynomial over the prime field::

            sage: # needs sage.rings.finite_rings
            sage: F2.<i> = GF((431,2), modulus=[1,0,1])
            sage: F6.<u> = F2.extension(3)
            sage: (u + 1).minpoly()                                                     # needs sage.modules
            x^6 + 425*x^5 + 19*x^4 + 125*x^3 + 189*x^2 + 239*x + 302
            sage: ext = F6.over(F2)                                                     # needs sage.modules
            sage: ext(u + 1).minpoly()  # indirect doctest                              # needs sage.modules # random
            x^3 + (396*i + 428)*x^2 + (80*i + 39)*x + 9*i + 178

        TESTS:

        We make sure that the previous example works on random examples::

            sage: # long time, needs sage.rings.finite_rings
            sage: p = random_prime(50)
            sage: K.<u> = GF((p, randrange(1,20)))
            sage: L.<v> = K.extension(randrange(2,20))
            sage: LK = L.over(K)
            sage: a = L.random_element()
            sage: poly = LK(a).minpoly()  # indirect doctest
            sage: poly(a)
            0
            sage: abs_deg = a.minpoly().degree()
            sage: poly.degree() == abs_deg // gcd(abs_deg, K.degree())
            True
        """
        poly = self.lift()
        try:
            return poly.minpoly_mod(self.parent().modulus())
        except AttributeError:
            pass
        return self.matrix().minpoly()

    def norm(self):
        """
        The norm of this element, which is the determinant of the matrix of right
        multiplication by this element.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 - 389*x^2 + 2*x - 5)
            sage: a.norm()                                                              # needs sage.modules
            5
        """
        return self.matrix().determinant()

    def trace(self):
        """
        The trace of this element, which is the trace of the matrix of
        right multiplication by this element.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 - 389*x^2 + 2*x - 5)
            sage: a.trace()                                                             # needs sage.modules
            389
        """
        return self.matrix().trace()

    def rational_reconstruction(self, *args, **kwargs):
        r"""
        Compute a rational reconstruction of this polynomial quotient
        ring element to its cover ring.

        This method is a thin convenience wrapper around
        :meth:`Polynomial.rational_reconstruction`.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = GF(65537)[]
            sage: m = (x^11 + 25345*x^10 + 10956*x^9 + 13873*x^8 + 23962*x^7
            ....:      + 17496*x^6 + 30348*x^5 + 7440*x^4 + 65438*x^3 + 7676*x^2
            ....:      + 54266*x + 47805)
            sage: f = (20437*x^10 + 62630*x^9 + 63241*x^8 + 12820*x^7 + 42171*x^6
            ....:      + 63091*x^5 + 15288*x^4 + 32516*x^3 + 2181*x^2 + 45236*x + 2447)
            sage: f_mod_m = R.quotient(m)(f)
            sage: f_mod_m.rational_reconstruction()
            (51388*x^5 + 29141*x^4 + 59341*x^3 + 7034*x^2 + 14152*x + 23746,
             x^5 + 15208*x^4 + 19504*x^3 + 20457*x^2 + 11180*x + 28352)
        """
        m = self.parent().modulus()
        R = m.parent()
        f = R(self._polynomial)
        return f.rational_reconstruction(m, *args, **kwargs)
