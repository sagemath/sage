# sage_setup: distribution = sagemath-flint
r"""
Univariate polynomials over number fields

AUTHOR:

- Luis Felipe Tabera Alonso (2014-02): initial version.

EXAMPLES:

Define a polynomial over an absolute number field and perform basic
operations with them::

    sage: x = polygen(ZZ, 'x')
    sage: N.<a> = NumberField(x^2 - 2)
    sage: K.<x> = N[]
    sage: f = x - a
    sage: g = x^3 - 2*a + 1
    sage: f * (x + a)
    x^2 - 2
    sage: f + g
    x^3 + x - 3*a + 1
    sage: g // f
    x^2 + a*x + 2
    sage: g % f
    1
    sage: factor(x^3 - 2*a*x^2 - 2*x + 4*a)
    (x - 2*a) * (x - a) * (x + a)
    sage: gcd(f, x - a)
    x - a

Polynomials are aware of embeddings of the underlying field::

    sage: # needs sage.rings.padics
    sage: x = polygen(ZZ, 'x')
    sage: Q7 = Qp(7)
    sage: r1 = Q7(3 + 7 + 2*7^2 + 6*7^3 + 7^4 + 2*7^5 + 7^6 + 2*7^7 + 4*7^8
    ....:          + 6*7^9 + 6*7^10 + 2*7^11 + 7^12 + 7^13 + 2*7^15 + 7^16 + 7^17
    ....:          + 4*7^18 + 6*7^19)
    sage: N.<b> = NumberField(x^2 - 2, embedding=r1)
    sage: K.<t> = N[]
    sage: f = t^3 - 2*t + 1
    sage: f(r1)
    1 + O(7^20)

We can also construct polynomials over relative number fields::

    sage: # needs sage.symbolic
    sage: N.<i, s2> = QQ[I, sqrt(2)]
    sage: K.<x> = N[]
    sage: f = x - s2
    sage: g = x^3 - 2*i*x^2 + s2*x
    sage: f * (x + s2)
    x^2 - 2
    sage: f + g
    x^3 - 2*I*x^2 + (sqrt2 + 1)*x - sqrt2
    sage: g // f
    x^2 + (-2*I + sqrt2)*x - 2*sqrt2*I + sqrt2 + 2
    sage: g % f
    -4*I + 2*sqrt2 + 2
    sage: factor(i*x^4 - 2*i*x^2 + 9*i)
    (I) * (x - I + sqrt2) * (x + I - sqrt2) * (x - I - sqrt2) * (x + I + sqrt2)
    sage: gcd(f, x - i)
    1
"""

# ****************************************************************************
#       Copyright (C) 2014 Luis Felipe Tabera Alonso <taberalf@unican.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.polynomial.polynomial_element_generic import Polynomial_generic_dense_field
from sage.rings.rational_field import QQ
from sage.structure.element import coerce_binop
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


class Polynomial_absolute_number_field_dense(Polynomial_generic_dense_field):
    """
    Class of dense univariate polynomials over an absolute number field.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        Create a new polynomial in the polynomial ring ``parent``.

        INPUT:

        - ``parent`` -- the polynomial ring in which to construct the
          element

        - ``x`` -- (default: ``None``) an object representing the
          polynomial, e.g. a list of coefficients.  See
          :meth:`sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field.__init__`
          for more details.

        - ``check`` -- boolean (default: ``True``); if ``True``, make sure that
          the coefficients of the polynomial are in the base ring

        - ``is_gen`` -- boolean (default: ``False``); if ``True``, `x` is the
          distinguished generator of the polynomial ring

        - ``construct`` -- boolean (default: ``False``); unused

        EXAMPLES::

            sage: P.<x> = QQ[I][]
            sage: f = P.random_element()
            sage: from sage.rings.polynomial.polynomial_number_field import Polynomial_absolute_number_field_dense
            sage: isinstance(f, Polynomial_absolute_number_field_dense)
            True
            sage: a = P(x)
            sage: a.is_gen()
            True
        """
        Polynomial_generic_dense_field.__init__(self, parent, x, check, is_gen, construct)

    @coerce_binop
    def gcd(self, other):
        """
        Compute the monic gcd of two univariate polynomials using PARI.

        INPUT:

        - ``other`` -- a polynomial with the same parent as ``self``

        OUTPUT: the monic gcd of ``self`` and ``other``

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: N.<a> = NumberField(x^3 - 1/2, 'a')
            sage: R.<r> = N['r']
            sage: f = (5/4*a^2 - 2*a + 4)*r^2 + (5*a^2 - 81/5*a - 17/2)*r + 4/5*a^2 + 24*a + 6
            sage: g = (5/4*a^2 - 2*a + 4)*r^2 + (-11*a^2 + 79/5*a - 7/2)*r - 4/5*a^2 - 24*a - 6
            sage: gcd(f, g**2)
            r - 60808/96625*a^2 - 69936/96625*a - 149212/96625
            sage: R = QQ[I]['x']
            sage: f = R.random_element(2)
            sage: g = f + 1
            sage: h = R.random_element(2).monic()
            sage: f *= h
            sage: g *= h
            sage: gcd(f, g) - h
            0
            sage: f.gcd(g) - h
            0

        TESTS:

        Test for degree one extensions::

            sage: x = polygen(ZZ, 'x')
            sage: N = NumberField(x - 3, 'a')
            sage: a = N.gen()
            sage: R = N['x']
            sage: f = R._random_nonzero_element()
            sage: g1 = R._random_nonzero_element()
            sage: g2 = g1 * R._random_nonzero_element() + 1
            sage: g1 *= f
            sage: g2 *= f
            sage: d = gcd(g1, g2)
            sage: f.monic() - d
            0
            sage: d.parent() is R
            True

        Test for coercion with other rings and force weird variables
        to test PARI behavior::

            sage: r = polygen(ZZ, 'r')
            sage: N = NumberField(r^2 - 2, 'r')
            sage: a = N.gen()
            sage: R = N['r']
            sage: r = R.gen()
            sage: f = N.random_element(4)*r + 1
            sage: g = ZZ['r']([1, 2, 3, 4, 5, 6, 7]); g
            7*r^6 + 6*r^5 + 5*r^4 + 4*r^3 + 3*r^2 + 2*r + 1
            sage: gcd(f, g) == gcd(g, f)
            True
            sage: h = f.gcd(g); h
            1
            sage: h.parent()
            Univariate Polynomial Ring in r over
             Number Field in r with defining polynomial r^2 - 2
            sage: gcd([a*r + 2, r^2 - 2])
            r + r
        """
        if self.is_zero():
            if other.is_zero():
                return self
            else:
                return other.monic()
        elif other.is_zero():
            return self.monic()
        elif self.degree() == 0 or other.degree() == 0:
            return self.parent().one()

        # If the extension is of degree one, use the gcd from QQ[x]
        if self.base_ring().degree().is_one():
            R = self.base_ring()
            a = self.change_ring(QQ)
            b = other.change_ring(QQ)
            g = a.gcd(b)
            return g.change_ring(R)

        h1 = self._pari_with_name('x')
        h2 = other._pari_with_name('x')
        g = h1.gcd(h2)
        return (self.parent()(g)).monic()


class Polynomial_relative_number_field_dense(Polynomial_generic_dense_field):
    """
    Class of dense univariate polynomials over a relative number field.
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        Create a new polynomial in the polynomial ring ``parent``.

        INPUT:

        - ``parent`` -- polynomial ring in which to construct the
          element

        - ``x`` -- (default: ``None``) an object representing the
          polynomial, e.g. a list of coefficients. See
          :meth:`sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_dense_field.__init__`
          for more details.

        - ``check`` -- boolean (default: ``True``); if ``True``, make sure that
          the coefficients of the polynomial are in the base ring

        - ``is_gen`` -- boolean (default: ``False``); if ``True``, ``x`` is the
          distinguished generator of the polynomial ring

        - ``construct`` -- boolean (default: ``False``); unused

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: f = NumberField([x^2 - 2, x^2 - 3], 'a')['x'].random_element()
            sage: from sage.rings.polynomial.polynomial_number_field import Polynomial_relative_number_field_dense
            sage: isinstance(f, Polynomial_relative_number_field_dense)
            True
        """
        Polynomial_generic_dense_field.__init__(self, parent, x, check, is_gen, construct)

    @coerce_binop
    def gcd(self, other):
        """
        Compute the monic gcd of two polynomials.

        Currently, the method checks corner cases in which one of the
        polynomials is zero or a constant. Then, computes an absolute
        extension and performs the computations there.

        INPUT:

        - ``other`` -- a polynomial with the same parent as ``self``

        OUTPUT: the monic gcd of ``self`` and ``other``

        See :meth:`Polynomial_absolute_number_field_dense.gcd` for
        more details.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: N = QQ[sqrt(2), sqrt(3)]
            sage: s2, s3 = N.gens()
            sage: x = polygen(N)
            sage: f = x^4 - 5*x^2 + 6
            sage: g = x^3 + (-2*s2 + s3)*x^2 + (-2*s3*s2 + 2)*x + 2*s3
            sage: gcd(f, g)
            x^2 + (-sqrt2 + sqrt3)*x - sqrt3*sqrt2
            sage: f.gcd(g)
            x^2 + (-sqrt2 + sqrt3)*x - sqrt3*sqrt2

        TESTS::

            sage: x = polygen(ZZ, 'x')
            sage: R = NumberField([x^2 - 2, x^2 - 3], 'a')['x']
            sage: f = R._random_nonzero_element()
            sage: g1 = R.random_element()
            sage: g2 = R.random_element()*g1 + 1
            sage: g1 *= f
            sage: g2 *= f
            sage: f.monic() - g1.gcd(g2)
            0

        Test for degree one extensions::

            sage: R = NumberField([x - 2, x + 1, x - 3], 'a')['x']
            sage: f = R.random_element(2)
            sage: g1 = R.random_element(2)
            sage: g2 = R.random_element(2)*g1 + 1
            sage: g1 *= f
            sage: g2 *= f
            sage: d = gcd(g1, g2)
            sage: d - f.monic()
            0
            sage: d.parent() is R
            True

        Test for hardcoded variables::

            sage: # needs sage.symbolic
            sage: R = N['sqrt2sqrt3']
            sage: x = R.gen()
            sage: f = x^2 - 2
            sage: g1 = x^2 - s3
            sage: g2 = x - s2
            sage: gcd(f, g1)
            1
            sage: gcd(f, g2)
            sqrt2sqrt3 - sqrt2
        """
        if self.is_zero():
            if other.is_zero():
                return self
            else:
                return other.monic()
        elif other.is_zero():
            return self.monic()
        elif self.degree() == 0 or other.degree() == 0:
            return self.parent().one()

        L = self.parent()
        x = L.variable_name()
        N = self.base_ring()
        c = ''.join(map(str,N.variable_names()))
        M = N.absolute_field(c)
        M_to_N, N_to_M = M.structure()
        R = PolynomialRing(M, x)
        first = R([N_to_M(foo) for foo in self.list()])
        second = R([N_to_M(foo) for foo in other.list()])
        result = first.gcd(second)
        result = L([M_to_N(foo) for foo in result.list()])
        # the result is already monic
        return result
