# sage.doctest: needs sage.rings.function_field
r"""
Ideals of function fields: extension
"""

# ****************************************************************************
#       Copyright (C) 2010      William Stein <wstein@gmail.com>
#                     2011      Maarten Derickx <m.derickx.student@gmail.com>
#                     2017-2021 Kwankyu Lee
#                     2018      Frédéric Chapoton
#                     2019      Brent Baccala
#                     2021      Jonathan Kliem
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

import itertools
from sage.rings.infinity import infinity
from sage.arith.power import generic_power
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.richcmp import richcmp
from sage.matrix.constructor import matrix

from .ideal import FunctionFieldIdeal, FunctionFieldIdealInfinite


class FunctionFieldIdeal_polymod(FunctionFieldIdeal):
    """
    Fractional ideals of algebraic function fields.

    INPUT:

    - ``ring`` -- order in a function field

    - ``hnf`` -- matrix in hermite normal form

    - ``denominator`` -- denominator

    The rows of ``hnf`` is a basis of the ideal, which itself is
    ``denominator`` times the fractional ideal.

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x^3*y - x)
        sage: O = L.maximal_order()
        sage: O.ideal(y)
        Ideal (y) of Maximal order of Function field in y defined by y^2 + x^3*y + x
    """
    def __init__(self, ring, hnf, denominator=1):
        """
        Initialize.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3*y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: TestSuite(I).run()
        """
        FunctionFieldIdeal.__init__(self, ring)

        # the denominator and the entries of the hnf are
        # univariate polynomials.
        self._hnf = hnf
        self._denominator = denominator

        # for prime ideals
        self._relative_degree = None
        self._ramification_index = None
        self._prime_below = None
        self._beta = None

        # beta in matrix form for fast multiplication
        self._beta_matrix = None

        # (p, q) with irreducible polynomial p and q an element of O in vector
        # form, together generating the prime ideal. This data is obtained by
        # Kummer's theorem when this prime ideal is constructed. This is used
        # for fast multiplication with other ideal.
        self._kummer_form = None

        # tuple of at most two gens:
        # the first gen is an element of the base ring of the maximal order
        # the second gen is the vector form of an element of the maximal order
        # if the second gen is zero, the tuple has only the first gen.
        self._gens_two_vecs = None

    def __bool__(self):
        """
        Test if this ideal is zero.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3*y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y); I
            Ideal (y) of Maximal order of Function field in y defined by y^2 + x^3*y + x
            sage: I.is_zero()
            False
            sage: J = 0*I; J
            Zero ideal of Maximal order of Function field in y defined by y^2 + x^3*y + x
            sage: J.is_zero()
            True

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y>=K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y); I
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: I.is_zero()
            False
            sage: J = 0*I; J
            Zero ideal of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: J.is_zero()
            True
        """
        return self._hnf.nrows() != 0

    def __hash__(self):
        """
        Return the hash of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: { I: 2 }[I] == 2
            True

            sage: # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: { I: 2 }[I] == 2
            True
        """
        return hash((self._ring, self._hnf, self._denominator))

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is in this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal([y]); I
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 + 6*x^3 + 6
            sage: x * y in I
            True
            sage: y / x in I
            False
            sage: y^2 - 2 in I
            False

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal([y]); I
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: x * y in I
            True
            sage: y / x in I
            False
            sage: y^2 - 2 in I
            False

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal([y]); I
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 - x^3 - 1
            sage: x * y in I
            True
            sage: y / x in I
            False
            sage: y^2 - 2 in I
            False

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal([y]); I
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: x * y in I
            True
            sage: y / x in I
            False
            sage: y^2 - 2 in I
            False
        """
        from sage.modules.free_module_element import vector

        vec = self.ring().coordinate_vector(self._denominator * x)
        v = []
        for e in vec:
            if e.denominator() != 1:
                return False
            v.append(e.numerator())
        vec = vector(v)
        return vec in self._hnf.image()

    def __invert__(self):
        """
        Return the inverse fractional ideal of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: ~I
            Ideal ((1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 + 6*x^3 + 6
            sage: I^(-1)
            Ideal ((1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 + 6*x^3 + 6
            sage: ~I * I
            Ideal (1) of Maximal order of Function field in y defined by y^2 + 6*x^3 + 6

        ::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: ~I
            Ideal ((x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: I^(-1)
            Ideal ((x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: ~I * I
            Ideal (1) of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x

        ::

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: ~I
            Ideal ((1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 - x^3 - 1
            sage: I^(-1)
            Ideal ((1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 - x^3 - 1
            sage: ~I * I
            Ideal (1) of Maximal order of Function field in y defined by y^2 - x^3 - 1

        ::

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: ~I
            Ideal ((x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: I^(-1)
            Ideal ((x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: ~I * I
            Ideal (1) of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        R = self.ring()
        T = R._codifferent_matrix()
        I = self * R.codifferent()
        J = I._denominator * (I._hnf * T).inverse()
        return R._ideal_from_vectors(J.columns())

    def _richcmp_(self, other, op):
        """
        Compare this ideal with the other ideal with respect to ``op``.

        INPUT:

        - ``other`` -- ideal

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: I == I + I
            True
            sage: I == I * I
            False

        ::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: I == I + I
            True
            sage: I == I * I
            False
            sage: I < I * I
            True
            sage: I > I * I
            False
        """
        return richcmp((self._denominator, self._hnf), (other._denominator, other._hnf), op)

    def _add_(self, other):
        """
        Add with other ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x + y)
            sage: I + J
            Ideal (y) of Maximal order of Function field in y defined by y^2 + x^3*y + x

            sage: # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x + y)
            sage: I + J
            Ideal (1, y) of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        ds = self._denominator
        do = other._denominator
        vecs1 = [do * r for r in self._hnf]
        vecs2 = [ds * r for r in other._hnf]
        return self._ring._ideal_from_vectors_and_denominator(vecs1 + vecs2, ds * do)

    def _mul_(self, other):
        """
        Multiply with other ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x + y)
            sage: I * J
            Ideal (x^4 + x^2 + x, x*y + x^2) of Maximal order
            of Function field in y defined by y^2 + x^3*y + x

            sage: # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x + y)
            sage: I * J
            Ideal ((x + 1)*y + (x^2 + 1)/x) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        O = self._ring
        mul = O._mul_vecs

        if self._kummer_form is not None:
            p, q = self._kummer_form
            vecs = list(p * other._hnf) + [mul(q, v) for v in other._hnf]
        elif other._kummer_form is not None:
            p, q = other._kummer_form
            vecs = list(p * self._hnf) + [mul(q, v) for v in self._hnf]
        elif self._gens_two_vecs is not None:
            if len(self._gens_two_vecs) == 1:
                g1, = self._gens_two_vecs
                vecs = list(g1 * other._hnf)
            else:
                g1, g2 = self._gens_two_vecs
                vecs = list(g1 * other._hnf) + [mul(g2, v) for v in other._hnf]
        elif other._gens_two_vecs is not None:
            if len(other._gens_two_vecs) == 1:
                g1, = other._gens_two_vecs
                vecs = list(g1 * self._hnf)
            else:
                g1, g2 = other._gens_two_vecs
                vecs = list(g1 * self._hnf) + [mul(g2, v) for v in self._hnf]
        else:
            vecs = [mul(r1, r2) for r1 in self._hnf for r2 in other._hnf]

        return O._ideal_from_vectors_and_denominator(vecs, self._denominator * other._denominator)

    def _acted_upon_(self, other, on_left):
        """
        Multiply ``other`` and this ideal on the right.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x + y)
            sage: J = O.ideal(x)
            sage: x * I ==  I * J
            True

            sage: # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x + y)
            sage: J = O.ideal(x)
            sage: x * I ==  I * J
            True
        """
        from sage.modules.free_module_element import vector

        O = self._ring
        mul = O._mul_vecs

        # compute the vector form of other element
        v = O.coordinate_vector(other)
        d = v.denominator()
        vec = vector([(d * c).numerator() for c in v])

        # multiply with the ideal
        vecs = [mul(vec, r) for r in self._hnf]

        return O._ideal_from_vectors_and_denominator(vecs, d * self._denominator)

    def intersect(self, other):
        """
        Intersect this ideal with the other ideal as fractional ideals.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x + y)
            sage: J = O.ideal(x)
            sage: I.intersect(J) == I * J * (I + J)^-1
            True
        """
        from sage.matrix.special import block_matrix
        from .hermite_form_polynomial import reversed_hermite_form

        A = self._hnf
        B = other._hnf

        ds = self.denominator()
        do = other.denominator()
        d = ds.lcm(do)
        if not d.is_one():
            A = (d // ds) * A
            B = (d // do) * B

        MS = A.matrix_space()
        I = MS.identity_matrix()
        O = MS.zero()
        n = A.ncols()

        # intersect the row spaces of A and B
        M = block_matrix([[I, I], [A, O], [O, B]])

        # reversed Hermite form
        U = reversed_hermite_form(M, transformation=True)

        vecs = [U[i][:n] for i in range(n)]

        return self._ring._ideal_from_vectors_and_denominator(vecs, d)

    def hnf(self):
        """
        Return the matrix in hermite normal form representing this ideal.

        See also :meth:`denominator`

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y*(y+1)); I.hnf()
            [x^6 + x^3         0]
            [  x^3 + 1         1]

        ::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y*(y+1)); I.hnf()
            [x^6 + x^3         0]
            [  x^3 + 1         1]
        """
        return self._hnf

    def denominator(self):
        """
        Return the denominator of this fractional ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y/(y+1))
            sage: d = I.denominator(); d
            x^3
            sage: d in O
            True

        ::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y/(y+1))
            sage: d = I.denominator(); d
            x^3
            sage: d in O
            True
        """
        return self._denominator

    @cached_method
    def module(self):
        """
        Return the module, that is the ideal viewed as a module
        over the base maximal order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: F.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x, 1/y)
            sage: I.module()
            Free module of degree 2 and rank 2 over Maximal order
            of Rational function field in x over Finite Field of size 7
            Echelon basis matrix:
            [          1           0]
            [          0 1/(x^3 + 1)]
        """
        F = self.ring().fraction_field()
        V, fr, to = F.vector_space()
        O = F.base_field().maximal_order()
        return V.span([to(g) for g in self.gens_over_base()], base_ring=O)

    @cached_method
    def gens_over_base(self):
        """
        Return the generators of this ideal as a module over the maximal order
        of the base rational function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x + y)
            sage: I.gens_over_base()
            (x^4 + x^2 + x, y + x)

            sage: # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x + y)
            sage: I.gens_over_base()
            (x^3 + 1, y + x)
        """
        gens, d = self._gens_over_base
        return tuple([~d * b for b in gens])

    @lazy_attribute
    def _gens_over_base(self):
        """
        Return the generators of the integral ideal, which is the denominator
        times the fractional ideal, together with the denominator.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3*y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(1/y)
            sage: I._gens_over_base
            ([x, y], x)
        """
        gens = [sum([c1 * c2 for c1, c2 in zip(row, self._ring.basis())])
                for row in self._hnf]
        return gens, self._denominator

    def gens(self):
        """
        Return a set of generators of this ideal.

        This provides whatever set of generators as quickly
        as possible.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x + y)
            sage: I.gens()
            (x^4 + x^2 + x, y + x)

            sage: # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x + y)
            sage: I.gens()
            (x^3 + 1, y + x)
        """
        return self.gens_over_base()

    @cached_method
    def basis_matrix(self):
        """
        Return the matrix of basis vectors of this ideal as a module.

        The basis matrix is by definition the hermite norm form of the ideal
        divided by the denominator.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x, 1/y)
            sage: I.denominator() * I.basis_matrix() == I.hnf()
            True
        """
        d = self.denominator()
        m = (d * self).hnf()
        if d != 1:
            m = ~d * m
        m.set_immutable()
        return m

    def is_integral(self):
        """
        Return ``True`` if this is an integral ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x, 1/y)
            sage: I.is_integral()
            False
            sage: J = I.denominator() * I
            sage: J.is_integral()
            True

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x, 1/y)
            sage: I.is_integral()
            False
            sage: J = I.denominator() * I
            sage: J.is_integral()
            True

            sage: K.<x> = FunctionField(QQ); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x, 1/y)
            sage: I.is_integral()
            False
            sage: J = I.denominator() * I
            sage: J.is_integral()
            True
        """
        return self.denominator() == 1

    def ideal_below(self):
        """
        Return the ideal below this ideal.

        This is defined only for integral ideals.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x, 1/y)
            sage: I.ideal_below()
            Traceback (most recent call last):
            ...
            TypeError: not an integral ideal
            sage: J = I.denominator() * I
            sage: J.ideal_below()
            Ideal (x^3 + x^2 + x) of Maximal order of Rational function field
            in x over Finite Field of size 2

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x, 1/y)
            sage: I.ideal_below()
            Traceback (most recent call last):
            ...
            TypeError: not an integral ideal
            sage: J = I.denominator() * I
            sage: J.ideal_below()
            Ideal (x^3 + x) of Maximal order of Rational function field
            in x over Finite Field of size 2

            sage: K.<x> = FunctionField(QQ); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x, 1/y)
            sage: I.ideal_below()
            Traceback (most recent call last):
            ...
            TypeError: not an integral ideal
            sage: J = I.denominator() * I
            sage: J.ideal_below()
            Ideal (x^3 + x^2 + x) of Maximal order of Rational function field
            in x over Rational Field
        """
        if not self.is_integral():
            raise TypeError("not an integral ideal")

        K = self.ring().fraction_field().base_field().maximal_order()

        # self._hnf is in reversed hermite normal form, that is, lower
        # triangular form. Thus the generator of the ideal below is
        # just the (0,0) entry of the normal form. When self._hnf was in
        # hermite normal form, that is, upper triangular form, then the
        # generator can be computed in the following way:
        #
        #   m = matrix([hnf[0].parent().gen(0)] + list(hnf))
        #   _,T = m.hermite_form(transformation=True)
        #   return T[-1][0]
        #
        # This is certainly less efficient! This is an argument for using
        # reversed hermite normal form for ideal representation.
        l = self._hnf[0][0]

        return K.ideal(l)

    def norm(self):
        """
        Return the norm of this fractional ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: i1 = O.ideal(x)
            sage: i2 = O.ideal(y)
            sage: i3 = i1 * i2
            sage: i3.norm() == i1.norm() * i2.norm()
            True
            sage: i1.norm()
            x^3
            sage: i1.norm() == x ** F.degree()
            True
            sage: i2.norm()
            x^6 + x^4 + x^2
            sage: i2.norm() == y.norm()
            True

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: i1 = O.ideal(x)
            sage: i2 = O.ideal(y)
            sage: i3 = i1 * i2
            sage: i3.norm() == i1.norm() * i2.norm()
            True
            sage: i1.norm()
            x^2
            sage: i1.norm() == x ** L.degree()
            True
            sage: i2.norm()
            (x^2 + 1)/x
            sage: i2.norm() == y.norm()
            True
        """
        n = 1
        for e in self.basis_matrix().diagonal():
            n *= e
        return n

    @cached_method
    def is_prime(self):
        """
        Return ``True`` if this ideal is a prime ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.is_prime() for f,_ in I.factor()]
            [True, True]

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.is_prime() for f,_ in I.factor()]
            [True, True]

            sage: K.<x> = FunctionField(QQ); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.is_prime() for f,_ in I.factor()]
            [True, True]
        """
        factors = self.factor()
        if len(factors) == 1 and factors[0][1] == 1:  # prime!
            prime = factors[0][0]
            assert self == prime
            self._relative_degree = prime._relative_degree
            self._ramification_index = prime._ramification_index
            self._prime_below = prime._prime_below
            self._beta = prime._beta
            self._beta_matrix = prime._beta_matrix
            return True
        else:
            return False

    ###################################################
    # The following methods are only for prime ideals #
    ###################################################

    def valuation(self, ideal):
        """
        Return the valuation of ``ideal`` at this prime ideal.

        INPUT:

        - ``ideal`` -- fractional ideal

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x, (1/(x^3 + x^2 + x))*y^2)
            sage: I.is_prime()
            True
            sage: J = O.ideal(y)
            sage: I.valuation(J)
            2

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.valuation(I) for f,_ in I.factor()]
            [-1, 2]

        The method closely follows Algorithm 4.8.17 of [Coh1993]_.
        """
        from sage.matrix.constructor import matrix

        if ideal.is_zero():
            return infinity

        O = self.ring()
        F = O.fraction_field()
        n = F.degree()

        # beta_matrix is for fast multiplication with beta. For performance,
        # this is computed here rather than when the prime is constructed.
        if self._beta_matrix is None:
            beta = self._beta
            m = []
            for i in range(n):
                mtable_row = O._mtable[i]
                c = sum(beta[j] * mtable_row[j] for j in range(n))
                m.append(c)
            self._beta_matrix = matrix(m)

        B = self._beta_matrix

        # Step 1: compute the valuation of the denominator times the ideal
        #
        # This part is highly optimized as it is critical for
        # overall performance of the function field machinery.
        p = self.prime_below().gen().numerator()
        h = ideal._hnf.list()
        val = min([c.valuation(p) for c in h])
        i = self._ramification_index * val
        while True:
            ppow = p ** val
            h = (matrix(n, [c // ppow for c in h]) * B).list()
            val = min([c.valuation(p) for c in h])
            if val.is_zero():
                break
            i += self._ramification_index * (val - 1) + 1

        # Step 2: compute the valuation of the denominator
        j = self._ramification_index * ideal.denominator().valuation(p)

        # Step 3: return the valuation
        return i - j

    def prime_below(self):
        """
        Return the prime lying below this prime ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.prime_below() for f,_ in I.factor()]
            [Ideal (x) of Maximal order of Rational function field in x
            over Finite Field of size 2, Ideal (x^2 + x + 1) of Maximal order
            of Rational function field in x over Finite Field of size 2]

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.prime_below() for f,_ in I.factor()]
            [Ideal (x) of Maximal order of Rational function field in x over Finite Field of size 2,
             Ideal (x + 1) of Maximal order of Rational function field in x over Finite Field of size 2]

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.prime_below() for f,_ in I.factor()]
            [Ideal (x) of Maximal order of Rational function field in x over Rational Field,
             Ideal (x^2 + x + 1) of Maximal order of Rational function field in x over Rational Field]
        """
        return self._prime_below

    def _factor(self):
        """
        Return the factorization of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()  # indirect doctest
            True
        """
        O = self.ring()
        F = O.fraction_field()
        o = F.base_field().maximal_order()

        # First we collect primes below self
        d = self._denominator
        i = d * self

        factors = []
        primes = set([o.ideal(p) for p, _ in d.factor()] + [p for p, _ in i.ideal_below().factor()])
        for prime in primes:
            qs = [q[0] for q in O.decomposition(prime)]
            for q in qs:
                exp = q.valuation(self)
                if exp != 0:
                    factors.append((q, exp))
        return factors


class FunctionFieldIdeal_global(FunctionFieldIdeal_polymod):
    """
    Fractional ideals of canonical function fields.

    INPUT:

    - ``ring`` -- order in a function field

    - ``hnf`` -- matrix in hermite normal form

    - ``denominator`` -- denominator

    The rows of ``hnf`` is a basis of the ideal, which itself is
    ``denominator`` times the fractional ideal.

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x^3*y - x)
        sage: O = L.maximal_order()
        sage: O.ideal(y)
        Ideal (y) of Maximal order of Function field in y defined by y^2 + x^3*y + x
    """
    def __init__(self, ring, hnf, denominator=1):
        """
        Initialize.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(5)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3*y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: TestSuite(I).run()
        """
        FunctionFieldIdeal_polymod.__init__(self, ring, hnf, denominator)

    def __pow__(self, mod):
        """
        Return ``self`` to the power of ``mod``.

        If a two-generators representation of ``self`` is known, it is used
        to speed up powering.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^7 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x + y)
            sage: S = I / J
            sage: a = S^100
            sage: _ = S.gens_two()
            sage: b = S^100  # faster
            sage: b == I^100 / J^100
            True
            sage: b == a
            True
        """
        if mod > 2 and self._gens_two_vecs is not None:
            O = self._ring
            mul = O._mul_vecs
            R = self._hnf.base_ring()
            n = self._hnf.ncols()

            I = matrix.identity(R, n)

            if len(self._gens_two_vecs) == 1:
                p, = self._gens_two_vecs
                ppow = p**mod
                J = [ppow * v for v in I]
            else:
                p, q = self._gens_two_vecs
                q = sum(e1 * e2 for e1, e2 in zip(O.basis(), q))
                ppow = p**mod
                qpow = O._coordinate_vector(q**mod)
                J = [ppow * v for v in I] + [mul(qpow, v) for v in I]

            return O._ideal_from_vectors_and_denominator(J, self._denominator**mod)

        return generic_power(self, mod)

    def gens(self):
        """
        Return a set of generators of this ideal.

        This provides whatever set of generators as quickly
        as possible.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x + y)
            sage: I.gens()
            (x^4 + x^2 + x, y + x)

            sage: # needs sage.rings.finite_rings
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(x + y)
            sage: I.gens()
            (x^3 + 1, y + x)
        """
        if self._gens_two.is_in_cache():
            return self._gens_two.cache
        else:
            return self.gens_over_base()

    def gens_two(self):
        r"""
        Return two generators of this fractional ideal.

        If the ideal is principal, one generator *may* be returned.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I  # indirect doctest
            Ideal (y) of Maximal order of Function field
            in y defined by y^3 + x^6 + x^4 + x^2
            sage: ~I  # indirect doctest
            Ideal ((1/(x^6 + x^4 + x^2))*y^2) of Maximal order of Function field
            in y defined by y^3 + x^6 + x^4 + x^2

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I  # indirect doctest
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: ~I  # indirect doctest
            Ideal ((x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        d = self.denominator()
        return tuple(e / d for e in self._gens_two())

    @cached_method
    def _gens_two(self):
        r"""
        Return a set of two generators of the integral ideal, that is
        the denominator times this fractional ideal.

        ALGORITHM:

        At most two generators are required to generate ideals in
        Dedekind domains.

        Lemma 4.7.9, algorithm 4.7.10, and exercise 4.29 of [Coh1993]_
        tell us that for an integral ideal `I` in a number field, if
        we pick `a` such that `\gcd(N(I), N(a)/N(I)) = 1`, then `a`
        and `N(I)` generate the ideal.  `N()` is the norm, and this
        result (presumably) generalizes to function fields.

        After computing `N(I)`, we search exhaustively to find `a`.

        .. TODO::

            Always return a single generator for a principal ideal.

            Testing for principality is not trivial.  Algorithm 6.5.10
            of [Coh1993]_ could probably be adapted for function fields.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x^2, x*y, x + y)
            sage: I._gens_two()
            (x, y)

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3))
            sage: _.<Y> = K[]
            sage: L.<y> = K.extension(Y - x)
            sage: y.zeros()[0].prime_ideal()._gens_two()
            (x,)
        """
        O = self.ring()
        F = O.fraction_field()

        if self._kummer_form is not None:  # prime ideal
            _g1, _g2 = self._kummer_form
            g1 = F(_g1)
            g2 = sum([c1 * c2 for c1, c2 in zip(_g2, O.basis())])
            if g2:
                self._gens_two_vecs = (_g1, _g2)
                return (g1, g2)
            else:
                self._gens_two_vecs = (_g1,)
                return (g1,)

        # ---- start to search for two generators

        hnf = self._hnf

        norm = 1
        for e in hnf.diagonal():
            norm *= e

        if norm.is_constant():  # unit ideal
            self._gens_two_vecs = (1,)
            return (F(1),)

        # one generator; see .ideal_below()
        _l = hnf[0][0]
        p = _l.degree()
        l = F(_l)

        if self._hnf == O.ideal(l)._hnf:  # principal ideal
            self._gens_two_vecs = (_l,)
            return (l,)

        R = hnf.base_ring()

        basis = []
        for row in hnf:
            basis.append(sum([c1 * c2 for c1, c2 in zip(row, O.basis())]))

        n = len(basis)
        alpha = None

        def check(alpha):
            alpha_norm = alpha.norm().numerator()  # denominator is 1
            return norm.gcd(alpha_norm // norm) == 1

        # Trial 1: search for alpha among generators
        for alpha in basis:
            if check(alpha):
                self._gens_two_vecs = (_l, O._coordinate_vector(alpha))
                return (l, alpha)

        # Trial 2: exhaustive search for alpha using only polynomials
        # with coefficients 0 or 1
        for d in range(p):
            G = itertools.product(itertools.product([0, 1], repeat=d + 1), repeat=n)
            for g in G:
                alpha = sum([R(c1) * c2 for c1, c2 in zip(g, basis)])
                if check(alpha):
                    self._gens_two_vecs = (_l, O._coordinate_vector(alpha))
                    return (l, alpha)

        # Trial 3: exhaustive search for alpha using all polynomials
        for d in range(p):
            G = itertools.product(R.polynomials(max_degree=d), repeat=n)
            for g in G:
                # discard duplicate cases
                if max(c.degree() for c in g) != d:
                    continue
                for j in range(n):
                    if g[j] != 0:
                        break
                if g[j].leading_coefficient() != 1:
                    continue

                alpha = sum([c1 * c2 for c1, c2 in zip(g, basis)])
                if check(alpha):
                    self._gens_two_vecs = (_l, O._coordinate_vector(alpha))
                    return (l, alpha)

        # should not reach here
        raise ValueError("no two generators found")


class FunctionFieldIdealInfinite_polymod(FunctionFieldIdealInfinite):
    """
    Ideals of the infinite maximal order of an algebraic function field.

    INPUT:

    - ``ring`` -- infinite maximal order of the function field

    - ``ideal`` -- ideal in the inverted function field

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
        sage: F.<y> = K.extension(t^3 + t^2 - x^4)
        sage: Oinf = F.maximal_order_infinite()
        sage: Oinf.ideal(1/y)
        Ideal (1/x^4*y^2) of Maximal infinite order of Function field
        in y defined by y^3 + y^2 + 2*x^4
    """
    def __init__(self, ring, ideal):
        """
        Initialize this ideal.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/y)
            sage: TestSuite(I).run()
        """
        FunctionFieldIdealInfinite.__init__(self, ring)
        self._ideal = ideal

    def __hash__(self):
        """
        Return the hash of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/y)
            sage: d = { I: 1 }

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/y)
            sage: d = { I: 1 }
        """
        return hash((self.ring(), self._ideal))

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is in this ideal.

        INPUT:

        - ``x`` -- element of the function field

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/y)
            sage: 1/y in I
            True
            sage: 1/x in I
            False
            sage: 1/x^2 in I
            True
        """
        F = self.ring().fraction_field()
        iF, from_iF, to_iF = F._inversion_isomorphism()
        return to_iF(x) in self._ideal

    def _add_(self, other):
        """
        Add this ideal with the ``other`` ideal.

        INPUT:

        - ``ideal`` -- ideal

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I + J
            Ideal (1/x) of Maximal infinite order of Function field in y
            defined by y^3 + y^2 + 2*x^4

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I + J
            Ideal (1/x) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
        """
        return FunctionFieldIdealInfinite_polymod(self._ring, self._ideal + other._ideal)

    def _mul_(self, other):
        """
        Multiply this ideal with the ``other`` ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I * J
            Ideal (1/x^7*y^2) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I * J
            Ideal (1/x^4*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
        """
        return FunctionFieldIdealInfinite_polymod(self._ring, self._ideal * other._ideal)

    def __pow__(self, n):
        """
        Raise this ideal to ``n``-th power.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: J = Oinf.ideal(1/x)
            sage: J^3
            Ideal (1/x^3) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4
        """
        return FunctionFieldIdealInfinite_polymod(self._ring, self._ideal ** n)

    def __invert__(self):
        """
        Return the inverted ideal of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: J = Oinf.ideal(y)
            sage: ~J
            Ideal (1/x^4*y^2) of Maximal infinite order
            of Function field in y defined by y^3 + y^2 + 2*x^4
            sage: J * ~J
            Ideal (1) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: J = Oinf.ideal(y)
            sage: ~J
            Ideal (1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: J * ~J
            Ideal (1) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
        """
        return FunctionFieldIdealInfinite_polymod(self._ring, ~ self._ideal)

    def _richcmp_(self, other, op):
        """
        Compare this ideal with the ``other`` ideal with respect to ``op``.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I * J == J * I
            True
            sage: I + J == J
            True
            sage: I + J == I
            False
            sage: (I < J) == (not J < I)
            True

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x^2*1/y)
            sage: J = Oinf.ideal(1/x)
            sage: I * J == J * I
            True
            sage: I + J == J
            True
            sage: I + J == I
            False
            sage: (I < J) == (not J < I)
            True
        """
        return richcmp(self._ideal, other._ideal, op)

    @property
    def _relative_degree(self):
        """
        Return the relative degree of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: [J._relative_degree for J,_ in I.factor()]
            [1]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        return self._ideal._relative_degree

    def gens(self):
        """
        Return a set of generators of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(x + y)
            sage: I.gens()
            (x, y, 1/x^2*y^2)

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(x + y)
            sage: I.gens()
            (x, y)
        """
        F = self.ring().fraction_field()
        iF, from_iF, to_iF = F._inversion_isomorphism()
        return tuple(from_iF(b) for b in self._ideal.gens())

    def gens_two(self):
        """
        Return a set of at most two generators of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(x + y)
            sage: I.gens_two()
            (x, y)

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(x + y)
            sage: I.gens_two()
            (x,)
        """
        F = self.ring().fraction_field()
        iF, from_iF, to_iF = F._inversion_isomorphism()
        return tuple(from_iF(b) for b in self._ideal.gens_two())

    def gens_over_base(self):
        """
        Return a set of generators of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(x + y)
            sage: I.gens_over_base()
            (x, y, 1/x^2*y^2)
        """
        F = self.ring().fraction_field()
        iF, from_iF, to_iF = F._inversion_isomorphism()
        return tuple(from_iF(g) for g in self._ideal.gens_over_base())

    def ideal_below(self):
        """
        Return a set of generators of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = K[]
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/y^2)
            sage: I.ideal_below()
            Ideal (x^3) of Maximal order of Rational function field
            in x over Finite Field in z2 of size 3^2
        """
        return self._ideal.ideal_below()

    def is_prime(self):
        """
        Return ``True`` if this ideal is a prime ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x^3*y^2) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4)^3
            sage: I.is_prime()
            False
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x)^2
            sage: I.is_prime()
            False
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True
        """
        return self._ideal.is_prime()

    @cached_method
    def prime_below(self):
        """
        Return the prime of the base order that underlies this prime ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x^3*y^2) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4)^3
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True
            sage: J.prime_below()
            Ideal (1/x) of Maximal infinite order of Rational function field
            in x over Finite Field in z2 of size 3^2

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x)^2
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True
            sage: J.prime_below()
            Ideal (1/x) of Maximal infinite order of Rational function field in x
            over Finite Field of size 2
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        F = self.ring().fraction_field()
        K = F.base_field()
        return K.maximal_order_infinite().prime_ideal()

    def valuation(self, ideal):
        """
        Return the valuation of ``ideal`` with respect to this prime ideal.

        INPUT:

        - ``ideal`` -- fractional ideal

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(y)
            sage: [f.valuation(I) for f,_ in I.factor()]
            [-1]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        return self._ideal.valuation(self.ring()._to_iF(ideal))

    def _factor(self):
        """
        Return factorization of the ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: Oinf = F.maximal_order_infinite()
            sage: f = 1/x
            sage: I = Oinf.ideal(f)
            sage: I._factor()
            [(Ideal ((1/(x^4 + x^3 + x^2))*y^2 + 1/x^2*y + 1) of Maximal infinite order of Function field in y defined by y^3 + x^6 + x^4 + x^2,
              1),
             (Ideal ((1/(x^4 + x^3 + x^2))*y^2 + 1) of Maximal infinite order of Function field in y defined by y^3 + x^6 + x^4 + x^2,
              1)]
        """
        if self._ideal.is_prime.is_in_cache() and self._ideal.is_prime():
            return [(self, 1)]

        O = self.ring()
        factors = []
        for iprime, exp in O._to_iF(self).factor():
            prime = FunctionFieldIdealInfinite_polymod(O, iprime)
            factors.append((prime, exp))
        return factors
