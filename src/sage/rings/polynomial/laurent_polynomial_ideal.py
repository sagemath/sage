# sage.doctest: needs sage.libs.singular sage.modules (because all doctests need laurent_polynomial_mpair, Groebner bases)
r"""
Ideals in Laurent polynomial rings.

For `R` a commutative ring, ideals in the Laurent polynomial ring
`R[x_1^{\pm 1}, x_2^{\pm 1}, \ldots, x_n^{\pm 1}]` are implemented as
ideals in the ordinary polynomial ring `R[x_1, \ldots, x_n]` which are
saturated with respect to the ideal `(x_1 \cdots x_n)`.

AUTHORS:

- Kiran S. Kedlaya (2020): initial implementation
"""
# ****************************************************************************
#       Copyright (C) 2020 Kiran S. Kedlaya <kedlaya@ucsd.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.ideal import Ideal_generic
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing_univariate
from sage.structure.richcmp import op_EQ, op_NE, op_LT, op_LE, op_GT, op_GE
from sage.arith.misc import GCD

class LaurentPolynomialIdeal( Ideal_generic ):
    def __init__(self, ring, gens, coerce=True, hint=None):
        r"""
        Create an ideal in a Laurent polynomial ring.

        To compute structural properties of an ideal in the Laurent polynomial ring
        `R[x_1^{\pm},\ldots,x_n^{\pm}]`, we form the corresponding ideal in the
        associated ordinary polynomial ring `R[x_1,\ldots,x_n]` which is saturated
        with respect to the ideal `(x_1 \cdots x_n)`. Since computing the saturation
        can be expensive, we employ some strategies to reduce the need for it.

        - We only create the polynomial ideal as needed.

        - For some operations, we try some superficial tests first. E.g., for
          comparisons, we first look directly at generators.
        - The attribute ``hint`` is a lower bound on the associated polynomial ideal.
          Hints are automatically forwarded by certain creation operations (such as
          sums and base extensions), and can be manually forwarded in other cases.

        INPUT:

        - ``ring`` -- the ring the ideal is defined in
        - ``gens`` -- list of generators for the ideal
        - ``coerce`` -- whether or not to coerce elements into ``ring``
        - ``hint`` -- an ideal in the associated polynomial ring (optional; see above)

        EXAMPLES::

            sage: R.<x,y> = LaurentPolynomialRing(IntegerRing(), 2, order='lex')
            sage: R.ideal([x, y])
            Ideal (x, y) of Multivariate Laurent Polynomial Ring in x, y
            over Integer Ring
            sage: R.<x0,x1> = LaurentPolynomialRing(GF(3), 2)
            sage: R.ideal([x0^2, x1^-3])
            Ideal (x0^2, x1^-3) of Multivariate Laurent Polynomial Ring in x0, x1
            over Finite Field of size 3

            sage: P.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: I = P.ideal([~x + ~y - 1])
            sage: print(I)
            Ideal (-1 + y^-1 + x^-1) of
             Multivariate Laurent Polynomial Ring in x, y over Rational Field
            sage: I.is_zero()
            False
            sage: (x^(-2) + x^(-1)*y^(-1) - x^(-1)) in I
            True

            sage: P.<x,y,z> = LaurentPolynomialRing(QQ, 3)
            sage: I1 = P.ideal([x*y*z + x*y + 2*y^2, x + z])
            sage: I2 = P.ideal([x*y*z + x*y + 2*y^2 + x + z, x + z])
            sage: I1 == I2
            True
            sage: I3 = P.ideal([x*y*z + x*y + 2*y^2 + x + z, x + z, y])
            sage: I1 < I3
            True
            sage: I1.minimal_associated_primes()
            (Ideal (-1/2*z^2 + y - 1/2*z, x + z) of Multivariate
              Laurent Polynomial Ring in x, y, z over Rational Field,)

            sage: K.<z> = CyclotomicField(4)                                            # needs sage.rings.number_field
            sage: J = I1.base_extend(K)                                                 # needs sage.rings.number_field
            sage: J.base_ring()                                                         # needs sage.rings.number_field
            Cyclotomic Field of order 4 and degree 2
        """
        Ideal_generic.__init__(self, ring, gens, coerce=coerce)
        self._poly_ring = ring.polynomial_ring()
        self._poly_ideal = None  # Create only as needed
        self._saturated = False
        if hint is None:
            self._hint = self._poly_ring.zero_ideal()
        else:
            self._hint = hint.change_ring(self._poly_ring)

    def set_hint(self, hint):
        """
        Set the hint of this ideal.

        The hint is an ideal of the associated polynomial ring, which is
        assumed to be contained in the associated ideal. It is used internally
        to speed up computation of the associated ideal in some cases;
        normally the end user will have no need to work with it directly.

        EXAMPLES::

            sage: P.<x,y,z> = LaurentPolynomialRing(QQ, 3)
            sage: I = P.ideal([x^2*y + 3*x*y^2])
            sage: I.hint()
            Ideal (0) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: I.set_hint(P.polynomial_ring().ideal([x + 3*y]))
            sage: I.hint()
            Ideal (x + 3*y) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        self._hint = hint

    def hint(self):
        """
        Return the hint of this ideal.

        The hint is an ideal of the associated polynomial ring, which is
        assumed to be contained in the associated ideal. It is used internally
        to speed up computation of the associated ideal in some cases;
        normally the end user will have no need to work with it directly.

        EXAMPLES::

            sage: P.<x,y,z> = LaurentPolynomialRing(QQ, 3)
            sage: I = P.ideal([x^2*y + 3*x*y^2])
            sage: I.hint()
            Ideal (0) of Multivariate Polynomial Ring in x, y, z over Rational Field
        """
        return self._hint

    # Comparisons, using the associated polynomial ideal.
    def _richcmp_(self, right_r, op):
        r"""
        Comparison of ``self`` and ``right_r``.

        When testing equality, we first check generators to save time.

        When testing ideal containment, we do not saturate the smaller ideal.

        EXAMPLES::

            sage: P.<x,y,z> = LaurentPolynomialRing(QQ, 3)
            sage: I1 = P.ideal([x*y*z+x*y+2*y^2, x+z])
            sage: I2 = P.ideal([x*y*z+x*y+2*y^2+x+z, x+z])
            sage: I1 == I2
            True
            sage: I1 != I2
            False
            sage: I3 = P.ideal([x*y*z+x*y+2*y^2+x+z, x+z, y])
            sage: I1 <= I3
            True
            sage: I3 >= I1
            True
            sage: I1 < I3
            True
            sage: I3 > I1
            True
        """
        if op in (op_EQ, op_NE):
            if set(self.gens()) == set(right_r.gens()): # Early abort
                return (op == op_EQ)
            return ((self.polynomial_ideal() == right_r.polynomial_ideal()) == (op == op_EQ))
        elif op == op_LE:
            if all(f in right_r.gens() for f in self.gens()): # Early abort
                return True
            return self.polynomial_ideal(saturate=False) <= right_r.polynomial_ideal()
        elif op == op_GE:
            return right_r._richcmp_(self, op_LE)
        elif op == op_LT:
            return self._richcmp_(right_r, op_LE) and self._richcmp_(right_r, op_NE)
        elif op == op_GT:
            return right_r._richcmp_(self, op_LE) and right_r._richcmp_(self, op_NE)
        else:
            raise ValueError("invalid comparison")

    def __contains__(self, f):
        """
        Implement containment testing (in) for Laurent polynomial ideals.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: I = P.ideal([x^2*y + 3*x*y^2])
            sage: x + 3*y in I
            True

        This also works in the univariate case::

            sage: P.<x> = LaurentPolynomialRing(QQ)
            sage: I = P.ideal([x^2 + 3*x])
            sage: 1 + 3*x^-1 in I
            True
        """
        if not f or f in self.gens():
            return True
        f = self.ring()(f)
        if isinstance(self.ring(), LaurentPolynomialRing_univariate):
            g = f.__reduce__()[1][1]
        else:
            g = f.__reduce__()[1][0]
        return (g in self.polynomial_ideal())

    def gens_reduced(self):
        """
        Return a reduced system of generators.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ)
            sage: J = P.ideal([x^2 - y^-2, x * y^3 + 2 * y^2+ y])
            sage: J.gens_reduced()
            (x + 6*y + 5, 3*y^2 + 4*y + 1)
        """
        R = self.ring()
        J = self.polynomial_ideal()
        return tuple([R(p) for p in J.gens()])

    # Operations on ideals

    def change_ring(self, R, hint=None):
        """
        Coerce an ideal into a new ring.

        This operation does not forward hints, but a new hint can be
        specified manually.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: I = P.ideal([x + y])
            sage: Q.<x,y,z> = LaurentPolynomialRing(QQ, 3)
            sage: I.change_ring(Q)
            Ideal (x + y) of Multivariate Laurent Polynomial Ring in x, y, z
             over Rational Field
        """
        return R.ideal(self.gens(), hint=hint)

    def base_extend(self, F):
        """
        Form the base extension of ``self`` to the base ring ``F``.

        This operation forwards hints.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: I = P.ideal([x + y])
            sage: K.<z> = CyclotomicField(3)                                            # needs sage.rings.number_field
            sage: I.base_extend(K)                                                      # needs sage.rings.number_field
            Ideal (x + y) of Multivariate Laurent Polynomial Ring in x, y
             over Cyclotomic Field of order 3 and degree 2
        """
        ring = self.ring()
        return self.change_ring(ring.change_ring(F), hint=self._hint)

    def apply_map(self, f, new_ring=None, new_base_ring=None, apply_to_hint=None):
        """
        Return the new ideal obtained by applying ``f`` to each generator.

        By default, this does not forward hints. To do so, set ``apply_to_hint``
        to specify a function to apply to generators of ``hint``.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: I = P.ideal([x + 1, y - 1])
            sage: I.apply_map(lambda z: z + 2)
            Ideal (x + 3, y + 1) of Multivariate Laurent Polynomial Ring in x, y
             over Rational Field
            sage: K.<i> = CyclotomicField(4)                                            # needs sage.rings.number_field
            sage: I.apply_map(lambda z: z + 2, new_base_ring=K)                         # needs sage.rings.number_field
            Ideal (x + 3, y + 1) of Multivariate Laurent Polynomial Ring in x, y
             over Cyclotomic Field of order 4 and degree 2
        """
        ring = self.ring()
        if new_ring is not None:
            R = new_ring
        elif new_base_ring is not None:
            R = ring.change_ring(new_base_ring)
        else:
            R = ring
        if apply_to_hint:
            hint = self._poly_ring.ideal([apply_to_hint(x) for x in self._hint.gens()])
        else:
            hint = None
        return R.ideal(tuple(f(x) for x in self.gens()), hint=hint)

    def apply_coeff_map(self, f, new_base_ring=None, forward_hint=True):
        """
        Apply a function to coefficients.

        This operation forwards hints by default.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K.<z> = CyclotomicField(3)
            sage: P.<x,y> = LaurentPolynomialRing(K, 2)
            sage: I = P.ideal([x + z, y - z])
            sage: h = K.hom([z^2])
            sage: I.apply_coeff_map(h)
            Ideal (x - z - 1, y + z + 1) of Multivariate Laurent Polynomial Ring
             in x, y over Cyclotomic Field of order 3 and degree 2
            sage: K1.<z1> = CyclotomicField(12)
            sage: h1 = K.hom([z1^4])
            sage: I.apply_coeff_map(h1, new_base_ring=K1)
            Ideal (x + z1^2 - 1, y - z1^2 + 1) of Multivariate Laurent Polynomial Ring
             in x, y over Cyclotomic Field of order 12 and degree 4
        """
        ring = self.ring()
        if new_base_ring is None:
            R = ring
        else:
            R = ring.change_ring(new_base_ring)
        if forward_hint:
            apply_to_hint = lambda x,f=f: x.map_coefficients(f)
        else:
            apply_to_hint = None
        return self.apply_map(lambda x,f=f:
                              x.map_coefficients(f, new_base_ring=new_base_ring),
                              new_ring=R, apply_to_hint=apply_to_hint)

    def toric_coordinate_change(self, M, forward_hint=True):
        """
        Compute the toric change of coordinates defined by the integer matrix ``M``.

        This operation forwards hints by default.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: K.<z> = CyclotomicField(3)
            sage: P.<x,y> = LaurentPolynomialRing(K, 2)
            sage: I = P.ideal([x + 1, y - 1])
            sage: M = Matrix([[2,1], [1,-3]])
            sage: I.toric_coordinate_change(M)
            Ideal (x^2*y + 1, -1 + x*y^-3) of Multivariate Laurent Polynomial Ring
             in x, y over Cyclotomic Field of order 3 and degree 2
        """
        if forward_hint:
            R = self.ring()
            apply_to_hint = lambda x, M=M, R=R: R(x).toric_coordinate_change(M).monomial_reduction()[0]
        else:
            apply_to_hint = None
        return self.apply_map(lambda x, M=M: x.toric_coordinate_change(M),
                              apply_to_hint=apply_to_hint)

    def __add__(self, other):
        """
        Return the sum of two ideals in the same ring.

        Currently this operation does not support coercion.

        This operation forwards hints.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: I = P.ideal([x+y])
            sage: J = P.ideal([y+1])
            sage: (I+J).groebner_basis()
            (x - 1, y + 1)
        """
        ring = self.ring()
        ring2 = other.ring()
        if ring != ring2:
            return ValueError("ambient rings are not equal")
        hint = self._hint
        if hint.is_zero():
            hint = other._hint
        elif not other._hint.is_zero():
            hint += other._hint
        return ring.ideal(self.gens() + other.gens(), hint=hint)

    def normalize_gens(self):
        """
        Redefine the ideal with a normalized set of generators.

        For two ideals returned by this function, equality testing can be done
        by comparing generators.

        This operation forwards hints.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: I = P.ideal([~x+y])
            sage: J = P.ideal([y+1])
            sage: I + J
            Ideal (y + x^-1, y + 1) of Multivariate Laurent Polynomial Ring
             in x, y over Rational Field
            sage: (I + J).normalize_gens()
            Ideal (x - 1, y + 1) of Multivariate Laurent Polynomial Ring
             in x, y over Rational Field
        """
        return self.ring().ideal(self.groebner_basis(), hint=self._hint)

    # Structural queries and properties

    def polynomial_ideal(self, saturate=True):
        """
        Return the associated polynomial ideal.

        By default, the ideal is saturated with respect to the product of the
        polynomial ring generators; this is necessary for testing equality and inclusion.
        As saturation can be quite time-consuming, it can be disabled by setting
        ``saturate=False``; however, the result will then depend not just on the original ideal
        but also on the choice of generators.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: I = P.ideal([x^2*y + 3*x*y^2])
            sage: I.polynomial_ideal()
            Ideal (x + 3*y) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: P.<t> = LaurentPolynomialRing(QQ)
            sage: J = P.ideal(t^2 - t^-1)
            sage: J.polynomial_ideal()
            Principal ideal (t^3 - 1) of Univariate Polynomial Ring in t over Rational Field
            sage: J = P.ideal([t^2 - t^-1, t + t^-1])
            sage: J.polynomial_ideal()
            Principal ideal (1) of Univariate Polynomial Ring in t over Rational Field
            sage: J = P.ideal([t^2 - t^-1, t - t^-1])
            sage: J.polynomial_ideal()
            Principal ideal (t - 1) of Univariate Polynomial Ring in t over Rational Field
        """
        P = self.ring()
        Q = self._poly_ring
        if isinstance(self.ring(), LaurentPolynomialRing_univariate):
            a = [Q(p.polynomial_construction()[0]) for p in self.gens()]
            if P.base_ring().is_field():
                a = GCD(a)
            return Q.ideal(a)
        if self._poly_ideal is not None and (self._saturated or not saturate):
            return self._poly_ideal
        gens = self.gens()
        if len(gens) == 0:
            id = Q.ideal([])
            self._poly_ideal = id
            self._hint = id
            self._saturated = True
            return id
        l2 = [f.__reduce__()[1][0] for f in gens]
        hint = self._hint
        l2 += list(hint.groebner_basis())
        I = Q.ideal(l2)
        if not saturate:
            self._poly_ideal = I
            self._hint = I
            return Q.ideal(l2)
        n = P.ngens()
        I = I.saturation(Q.ideal([Q.monomial(*((1,) * n))]))[0]
        self._poly_ideal = I
        self._hint = I
        self._saturated = True
        return I

    def groebner_basis(self, saturate=True):
        """
        Return the reduced Groebner basis for the specified term order.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: I = P.ideal([x + y])
            sage: J = P.ideal([y + 1])
            sage: (I + J).groebner_basis()
            (x - 1, y + 1)
        """
        l = self.polynomial_ideal(saturate=saturate).groebner_basis()
        return tuple(self.ring()(x) for x in l)

    def is_one(self):
        """
        Determine whether or not ``self`` is the unit ideal.

        This requires saturation of the polynomial ideal.

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: I = P.ideal([~x + ~y - 1, x + y])
            sage: I.is_one()
            True
        """
        return self.polynomial_ideal().is_one()

    def is_binomial(self, groebner_basis=False):
        """
        Determine whether every generator of ``self`` is a binomial.

        If ``groebner_basis`` is True, this becomes intrinsic (for a choice of
        term order).

        EXAMPLES::

            sage: P.<x,y> = LaurentPolynomialRing(QQ, 2)
            sage: I = P.ideal([x + y])
            sage: I.is_binomial()
            True
        """
        if groebner_basis:
            l = self.groebner_basis()
        else:
            l = self.gens()
        return all(not f or f.number_of_terms() == 2 for f in l)

    def associated_primes(self):
        """
        Return associated primes of this ideal.

        These are computed from the polynomial ideal, but it is not necessary to
        saturate. Instead, we omit primes containing any polynomial generator.

        EXAMPLES::

            sage: P.<x,y,z> = LaurentPolynomialRing(QQ, 3)
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = P.ideal((p*q^2, y - z^2))
            sage: tuple(sorted(I.associated_primes(), key=str))
            (Ideal (y + 1, z^2 + 1) of
              Multivariate Laurent Polynomial Ring in x, y, z over Rational Field,
             Ideal (z^2 - y, y*z + 2, y^2 + 2*z) of
              Multivariate Laurent Polynomial Ring in x, y, z over Rational Field)
        """
        l = self.polynomial_ideal(saturate=False).associated_primes()
        l2 = [self.ring().ideal(I.gens(), hint=I) for I in l]
        return tuple(I for I in l2 if not I.is_one())

    def minimal_associated_primes(self, saturate=False):
        """
        Return minimal associated primes of this ideal.

        These are computed from the polynomial ideal, but it is not necessary to
        saturate. Instead, we omit primes containing any polynomial generator.

        EXAMPLES::

            sage: P.<x,y,z> = LaurentPolynomialRing(QQ, 3)
            sage: p = z^2 + 1; q = z^3 + 2
            sage: I = P.ideal((p*q^2, y - z^2))
            sage: tuple(sorted(I.minimal_associated_primes(), key=str))
            (Ideal (z^2 + 1, -z^2 + y) of
              Multivariate Laurent Polynomial Ring in x, y, z over Rational Field,
             Ideal (z^3 + 2, -z^2 + y) of
              Multivariate Laurent Polynomial Ring in x, y, z over Rational Field)
        """
        l = self.polynomial_ideal(saturate=saturate).minimal_associated_primes()
        l2 = [self.ring().ideal(I.gens(), hint=I) for I in l]
        return tuple(I for I in l2 if not I.is_one())

    def radical(self):
        """
        Return the radical of this ideal.

        EXAMPLES::

            sage: P.<x,y,z> = LaurentPolynomialRing(QQ, 3)
            sage: I = P.ideal(((x+1)^2, (y+1)^3, ((x+1)*z)^4 + (y+1)^3 + 10*(x+1)^2))
            sage: I.radical()
            Ideal (y + 1, x + 1) of Multivariate Laurent Polynomial Ring in x, y, z
             over Rational Field
        """
        J = self.polynomial_ideal().radical()
        return self.ring().ideal(J.gens())

    def dimension(self):
        """
        Return the dimension of this ideal.

        EXAMPLES::

            sage: P.<x,y,z> = LaurentPolynomialRing(QQ, 3)
            sage: I = P.ideal(((x+1)^2, (y+1)^3, ((x+1)*z)^4 + (y+1)^3 + 10*(x+1)^2))
            sage: I.dimension()
            1
        """
        return self.polynomial_ideal().dimension()
