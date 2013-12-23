"""
Univariate Polynomials over domains and fields

AUTHORS:

- William Stein: first version
- Martin Albrecht: Added singular coercion.
- David Harvey: split off polynomial_integer_dense_ntl.pyx (2007-09)
- Robert Bradshaw: split off polynomial_modn_dense_ntl.pyx (2007-09)

TESTS:

We test coercion in a particularly complicated situation::

    sage: W.<w>=QQ['w']
    sage: WZ.<z>=W['z']
    sage: m = matrix(WZ,2,2,[1,z,z,z^2])
    sage: a = m.charpoly()
    sage: R.<x> = WZ[]
    sage: R(a)
    x^2 + (-z^2 - 1)*x
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.polynomial.polynomial_element import Polynomial, Polynomial_generic_dense
from sage.structure.element import IntegralDomainElement, EuclideanDomainElement

from sage.rings.polynomial.polynomial_singular_interface import Polynomial_singular_repr

from sage.libs.pari.all import pari_gen
from sage.structure.element import coerce_binop

from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ


class Polynomial_generic_sparse(Polynomial):
    """
    A generic sparse polynomial.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(PolynomialRing(QQ, 'y'), sparse=True)
        sage: f = x^3 - x + 17
        sage: type(f)
        <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_sparse'>
        sage: loads(f.dumps()) == f
        True

    A more extensive example::

        sage: A.<T> = PolynomialRing(Integers(5),sparse=True) ; f = T^2+1 ; B = A.quo(f)
        sage: C.<s> = PolynomialRing(B)
        sage: C
        Univariate Polynomial Ring in s over Univariate Quotient Polynomial Ring in Tbar over Ring of integers modulo 5 with modulus T^2 + 1
        sage: s + T
        s + Tbar
        sage: (s + T)**2
        s^2 + 2*Tbar*s + 4
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)
        if x is None:
            self.__coeffs = {}
            return
        R = parent.base_ring()
        if isinstance(x, Polynomial):
            if x.parent() == self.parent():
                x = dict(x.dict())
            elif x.parent() == R:
                x = {0:x}
            else:
                w = {}
                for n, c in x.dict().iteritems():
                    w[n] = R(c)
                # The following line has been added in trac ticket #9944.
                # Apparently, the "else" case has never occured before.
                x = w
        elif isinstance(x, list):
            y = {}
            for i in xrange(len(x)):
                if x[i] != 0:
                    y[i] = x[i]
            x = y
        elif isinstance(x, pari_gen):
            y = {}
            for i in range(len(x)):
                y[i] = R(x[i])
            x = y
            check = True
        elif not isinstance(x, dict):
            x = {0:x}   # constant polynomials
        if check:
            self.__coeffs = {}
            for i, z in x.iteritems():
                self.__coeffs[i] = R(z)
        else:
            self.__coeffs = x
        if check:
            self.__normalize()

    def dict(self):
        """
        Return a new copy of the dict of the underlying
        elements of self.

        EXAMPLES::

            sage: R.<w> = PolynomialRing(Integers(8), sparse=True)
            sage: f = 5 + w^1997 - w^10000; f
            7*w^10000 + w^1997 + 5
            sage: d = f.dict(); d
            {0: 5, 10000: 7, 1997: 1}
            sage: d[0] = 10
            sage: f.dict()
            {0: 5, 10000: 7, 1997: 1}
        """
        return dict(self.__coeffs)

    def coefficients(self):
        """
        Return the coefficients of the monomials appearing in self.

        EXAMPLES::

            sage: R.<w> = PolynomialRing(Integers(8), sparse=True)
            sage: f = 5 + w^1997 - w^10000; f
            7*w^10000 + w^1997 + 5
            sage: f.coefficients()
            [5, 1, 7]
        """
        return [c[1] for c in sorted(self.__coeffs.iteritems())]

    def exponents(self):
        """
        Return the exponents of the monomials appearing in self.

        EXAMPLES::

            sage: R.<w> = PolynomialRing(Integers(8), sparse=True)
            sage: f = 5 + w^1997 - w^10000; f
            7*w^10000 + w^1997 + 5
            sage: f.exponents()
            [0, 1997, 10000]
        """
        return [c[0] for c in sorted(self.__coeffs.iteritems())]

    def valuation(self):
        """
        EXAMPLES::

            sage: R.<w> = PolynomialRing(GF(9,'a'), sparse=True)
            sage: f = w^1997 - w^10000
            sage: f.valuation()
            1997
            sage: R(19).valuation()
            0
            sage: R(0).valuation()
            +Infinity
        """
        c = self.__coeffs.keys()
        if len(c) == 0:
            return infinity
        return ZZ(min(self.__coeffs.keys()))

    def _derivative(self, var=None):
        """
        Computes formal derivative of this polynomial with respect to
        the given variable.

        If var is None or is the generator of this ring, the derivative
        is with respect to the generator. Otherwise, _derivative(var) is called
        recursively for each coefficient of this polynomial.

        .. seealso:: :meth:`.derivative`

        EXAMPLES::

            sage: R.<w> = PolynomialRing(ZZ, sparse=True)
            sage: f = R(range(9)); f
            8*w^8 + 7*w^7 + 6*w^6 + 5*w^5 + 4*w^4 + 3*w^3 + 2*w^2 + w
            sage: f._derivative()
            64*w^7 + 49*w^6 + 36*w^5 + 25*w^4 + 16*w^3 + 9*w^2 + 4*w + 1
            sage: f._derivative(w)
            64*w^7 + 49*w^6 + 36*w^5 + 25*w^4 + 16*w^3 + 9*w^2 + 4*w + 1

            sage: R.<x> = PolynomialRing(ZZ)
            sage: S.<y> = PolynomialRing(R, sparse=True)
            sage: f = x^3*y^4
            sage: f._derivative()
            4*x^3*y^3
            sage: f._derivative(y)
            4*x^3*y^3
            sage: f._derivative(x)
            3*x^2*y^4
        """
        P = self.parent()
        if var is not None and var is not P.gen():
            # call _derivative() recursively on coefficients
            return P(dict([(n, c._derivative(var)) \
                                     for (n, c) in self.__coeffs.iteritems()]))

        # compute formal derivative with respect to generator
        d = {}
        for n, c in self.__coeffs.iteritems():
            d[n-1] = n*c
        if d.has_key(-1):
            del d[-1]
        return P(d)

    def _dict_unsafe(self):
        """
        Return unsafe access to the underlying dictionary of coefficients.

        ** DO NOT use this, unless you really really know what you are doing. **

        EXAMPLES:
            sage: R.<w> = PolynomialRing(ZZ, sparse=True)
            sage: f = w^15 - w*3; f
            w^15 - 3*w
            sage: d = f._dict_unsafe(); d
            {1: -3, 15: 1}
            sage: d[1] = 10; f
            w^15 + 10*w
        """
        return self.__coeffs

    def _repr(self, name=None):
        r"""
        EXAMPLES::

            sage: R.<w> = PolynomialRing(CDF, sparse=True)
            sage: f = CDF(1,2) + w^5 - CDF(pi)*w + CDF(e)
            sage: f._repr()
            '1.0*w^5 - 3.14159265359*w + 3.71828182846 + 2.0*I'
            sage: f._repr(name='z')
            '1.0*z^5 - 3.14159265359*z + 3.71828182846 + 2.0*I'

        AUTHOR:

        - David Harvey (2006-08-05), based on Polynomial._repr()
        - Francis Clarke (2008-09-08) improved for 'negative' coefficients
        """
        s = " "
        m = self.degree() + 1
        if name is None:
            name = self.parent().variable_name()
        atomic_repr = self.parent().base_ring()._repr_option('element_is_atomic')
        coeffs = list(self.__coeffs.iteritems())
        coeffs.sort()
        for (n, x) in reversed(coeffs):
            if x != 0:
                if n != m-1:
                    s += " + "
                x = y = repr(x)
                if y.find("-") == 0:
                    y = y[1:]
                if not atomic_repr and n > 0 and (y.find("+") != -1 or y.find("-") != -1):
                    x = "(%s)"%x
                if n > 1:
                    var = "*%s^%s"%(name,n)
                elif n==1:
                    var = "*%s"%name
                else:
                    var = ""
                s += "%s%s"%(x,var)
        s = s.replace(" + -", " - ")
        s = s.replace(" 1*"," ")
        s = s.replace(" -1*", " -")
        if s==" ":
            return "0"
        return s[1:]

    def __normalize(self):
        x = self.__coeffs
        D = [n for n, z in x.iteritems() if not z]
        for n in D:
            del x[n]

    def __getitem__(self,n):
        """
        Return the `n`-th coefficient of this polynomial if `n` is an integer,
        returns the monomials of self of degree in slice `n` if `n` is a slice.

        Negative indexes are allowed and always return 0 (so you can
        view the polynomial as embedding Laurent series).

        EXAMPLES::

            sage: R.<w> = PolynomialRing(RDF, sparse=True)
            sage: e = RDF(e)
            sage: f = sum(e^n*w^n for n in range(4)); f
            20.0855369232*w^3 + 7.38905609893*w^2 + 2.71828182846*w + 1.0
            sage: f[1]
            2.71828182846
            sage: f[5]
            0.0
            sage: f[-1]
            0.0
            sage: R.<x> = PolynomialRing(RealField(19), sparse=True)
            sage: f = (2-3.5*x)^3; f
            -42.875*x^3 + 73.500*x^2 - 42.000*x + 8.0000
            sage: f[1:3]
            73.500*x^2 - 42.000*x
            sage: f[:2]
            -42.000*x + 8.0000
            sage: f[2:]
            -42.875*x^3 + 73.500*x^2
        """
        if isinstance(n, slice):
            start, stop = n.start, n.stop
            if start < 0:
                start = 0
            if stop is None:
                stop = len(self.__coeffs) + 1
            v = {}
            x = self.__coeffs
            for k in x.keys():
                if start <= k and k < stop:
                    v[k] = x[k]
            P = self.parent()
            return P(v)
        else:
            if not self.__coeffs.has_key(n):
                return self.base_ring()(0)
            return self.__coeffs[n]

    def _unsafe_mutate(self, n, value):
        r"""
        Change the coefficient of `x^n` to value.

        ** NEVER USE THIS ** -- unless you really know what you are doing.

        EXAMPLES::

            sage: R.<z> = PolynomialRing(CC, sparse=True)
            sage: f = z^2 + CC.0; f
            1.00000000000000*z^2 + 1.00000000000000*I
            sage: f._unsafe_mutate(0, 10)
            sage: f
            1.00000000000000*z^2 + 10.0000000000000

        Much more nasty:
            sage: z._unsafe_mutate(1, 0)
            sage: z
            0
        """
        n = int(n)
        value = self.base_ring()(value)
        x = self.__coeffs
        if n < 0:
            raise IndexError, "polynomial coefficient index must be nonnegative"
        if value == 0:
            if x.has_key(n):
                del x[n]
        else:
            x[n] = value

    def list(self):
        """
        Return a new copy of the list of the underlying
        elements of self.

        EXAMPLES::

            sage: R.<z> = PolynomialRing(Integers(100), sparse=True)
            sage: f = 13*z^5 + 15*z^2 + 17*z
            sage: f.list()
            [0, 17, 15, 0, 0, 13]
        """
        zero = self.base_ring()(0)
        v = [zero] * (self.degree()+1)
        for n, x in self.__coeffs.iteritems():
            v[n] = x
        return v

    #def _pari_(self, variable=None):
    #    if variable is None:
    #        return self.__pari
    #    else:
    #        return self.__pari.subst('x',variable)

    def degree(self, gen=None):
        """
        Return the degree of this sparse polynomial.

        EXAMPLES::

            sage: R.<z> = PolynomialRing(ZZ, sparse=True)
            sage: f = 13*z^50000 + 15*z^2 + 17*z
            sage: f.degree()
            50000
        """
        v = self.__coeffs.keys()
        if len(v) == 0:
            return -1
        return max(v)

    def _add_(self, right):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(), sparse=True)
            sage: (x^100000 + 2*x^50000) + (4*x^75000 - 2*x^50000 + 3*x)
            x^100000 + 4*x^75000 + 3*x

        AUTHOR:

        - David Harvey (2006-08-05)
        """
        output = dict(self.__coeffs)

        for (index, coeff) in right.__coeffs.iteritems():
            if index in output:
                output[index] += coeff
            else:
                output[index] = coeff

        output = self.parent()(output, check=False)
        output.__normalize()
        return output

    def _neg_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(Integers(), sparse=True)
            sage: a = x^10000000; a
            x^10000000
            sage: -a
            -x^10000000
        """
        output = { }
        for (index, coeff) in self.__coeffs.iteritems():
            output[index] = -coeff
        output = self.parent()(output, check=False)
        return output

    def _mul_(self, right):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (x^100000 - x^50000) * (x^100000 + x^50000)
             x^200000 - x^100000
            sage: (x^100000 - x^50000) * R(0)
             0

        AUTHOR:
        - David Harvey (2006-08-05)
        """
        output = {}

        for (index1, coeff1) in self.__coeffs.iteritems():
            for (index2, coeff2) in right.__coeffs.iteritems():
                product = coeff1 * coeff2
                index = index1 + index2
                if index in output:
                    output[index] += product
                else:
                    output[index] = product

        output = self.parent()(output, check=False)
        output.__normalize()
        return output

    def _rmul_(self, left):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (x^100000 - x^50000) * (x^100000 + x^50000)
            x^200000 - x^100000
            sage: 7 * (x^100000 - x^50000)   # indirect doctest
            7*x^100000 - 7*x^50000

        AUTHOR:

        - Simon King (2011-03-31)
        """
        output = {}

        for (index, coeff) in self.__coeffs.iteritems():
            output[index] = left * coeff

        output = self.parent()(output, check=False)
        output.__normalize()
        return output

    def _lmul_(self, right):
        r"""
        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: (x^100000 - x^50000) * (x^100000 + x^50000)
            x^200000 - x^100000
            sage: (x^100000 - x^50000) * 7   # indirect doctest
            7*x^100000 - 7*x^50000

        AUTHOR:

        - Simon King (2011-03-31)
        """
        output = {}

        for (index, coeff) in self.__coeffs.iteritems():
            output[index] = coeff * right

        output = self.parent()(output, check=False)
        output.__normalize()
        return output

    def shift(self, n):
        r"""
        Returns this polynomial multiplied by the power `x^n`. If `n` is negative,
        terms below `x^n` will be discarded. Does not change this polynomial.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ, sparse=True)
            sage: p = x^100000 + 2*x + 4
            sage: type(p)
            <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_sparse'>
            sage: p.shift(0)
             x^100000 + 2*x + 4
            sage: p.shift(-1)
             x^99999 + 2
            sage: p.shift(-100002)
             0
            sage: p.shift(2)
             x^100002 + 2*x^3 + 4*x^2

        AUTHOR:
        - David Harvey (2006-08-06)
        """
        n = int(n)
        if n == 0:
            return self
        if n > 0:
            output = {}
            for (index, coeff) in self.__coeffs.iteritems():
                output[index + n] = coeff
            return self.parent()(output, check=False)
        if n < 0:
            output = {}
            for (index, coeff) in self.__coeffs.iteritems():
                if index + n >= 0:
                    output[index + n] = coeff
            return self.parent()(output, check=False)


class Polynomial_generic_domain(Polynomial, IntegralDomainElement):
    def __init__(self, parent, is_gen=False, construct=False):
        Polynomial.__init__(self, parent, is_gen=is_gen)

    def is_unit(self):
        r"""
        Return True if this polynomial is a unit.

        *EXERCISE* (Atiyah-McDonald, Ch 1): Let `A[x]` be a polynomial
        ring in one variable.  Then `f=\sum a_i x^i \in A[x]` is a
        unit if and only if `a_0` is a unit and `a_1,\ldots, a_n` are
        nilpotent.

        EXAMPLES::

            sage: R.<z> = PolynomialRing(ZZ, sparse=True)
            sage: (2 + z^3).is_unit()
            False
            sage: f = -1 + 3*z^3; f
            3*z^3 - 1
            sage: f.is_unit()
            False
            sage: R(-3).is_unit()
            False
            sage: R(-1).is_unit()
            True
            sage: R(0).is_unit()
            False
        """
        if self.degree() > 0:
            return False
        return self[0].is_unit()

class Polynomial_generic_field(Polynomial_singular_repr,
                               Polynomial_generic_domain,
                               EuclideanDomainElement):

    @coerce_binop
    def quo_rem(self, other):
        """
        Returns a tuple (quotient, remainder) where
            self = quotient*other + remainder.

        EXAMPLES::

            sage: R.<y> = PolynomialRing(QQ)
            sage: K.<t> = NumberField(y^2 - 2)
            sage: P.<x> = PolynomialRing(K)
            sage: x.quo_rem(K(1))
            (x, 0)
            sage: x.xgcd(K(1))
            (1, 0, 1)
        """
        P = self.parent()
        if other.is_zero():
            raise ZeroDivisionError, "other must be nonzero"

        # This is algorithm 3.1.1 in Cohen GTM 138
        A = self
        B = other
        R = A
        Q = P.zero_element()
        while R.degree() >= B.degree():
            aaa = R.leading_coefficient()/B.leading_coefficient()
            diff_deg=R.degree()-B.degree()
            Q += P(aaa).shift(diff_deg)
            # We know that S*B exactly cancels the leading coefficient of R.
            # Thus, we skip the computation of this leading coefficient.
            # For most exact fields, this doesn't matter much; but for
            # inexact fields, the leading coefficient might not end up
            # exactly equal to zero; and for AA/QQbar, verifying that
            # the coefficient is exactly zero triggers exact computation.
            R = R[:R.degree()] - (aaa*B[:B.degree()]).shift(diff_deg)
        return (Q, R)

    def _gcd(self, other):
        """
        Return the GCD of self and other, as a monic polynomial.
        """
        g = EuclideanDomainElement._gcd(self, other)
        c = g.leading_coefficient()
        if c.is_unit():
            return (1/c)*g
        return g

    def _xgcd(self, other):
        r"""
        Extended gcd of ``self`` and polynomial ``other``.

        INPUT:

        - ``other`` -- a polynomial defined over the same ring as ``self``

        OUTPUT:

        Polynomials ``g``, ``u``, and ``v`` such that ``g = u * self + v * other``.

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: F = (x^2 + 2)*x^3; G = (x^2+2)*(x-3)
            sage: g, u, v = F.xgcd(G)
            sage: g, u, v
            (x^2 + 2, 1/27, -1/27*x^2 - 1/9*x - 1/3)
            sage: u*F + v*G
            x^2 + 2

        ::

            sage: g, u, v = x.xgcd(P(0)); g, u, v
            (x, 1, 0)
            sage: g == u*x + v*P(0)
            True
            sage: g, u, v = P(0).xgcd(x); g, u, v
            (x, 0, 1)
            sage: g == u*P(0) + v*x
            True

        """
        if other.is_zero():
            R = self.parent()
            return self, R.one_element(), R.zero_element()
        # Algorithm 3.2.2 of Cohen, GTM 138
        R = self.parent()
        A = self
        B = other
        U = R.one_element()
        G = A
        V1 = R.zero_element()
        V3 = B
        while not V3.is_zero():
            Q, R = G.quo_rem(V3)
            T = U - V1*Q
            U = V1
            G = V3
            V1 = T
            V3 = R
        V = (G-A*U)//B
        lc = G.leading_coefficient()
        return G/lc, U/lc, V/lc


class Polynomial_generic_sparse_field(Polynomial_generic_sparse, Polynomial_generic_field):
    """
    EXAMPLES::

        sage: R.<x> = PolynomialRing(Frac(RR['t']), sparse=True)
        sage: f = x^3 - x + 17
        sage: type(f)
        <class 'sage.rings.polynomial.polynomial_element_generic.Polynomial_generic_sparse_field'>
        sage: loads(f.dumps()) == f
        True
    """
    def __init__(self, parent, x=None, check=True, is_gen = False, construct=False):
        Polynomial_generic_sparse.__init__(self, parent, x, check, is_gen)


class Polynomial_generic_dense_field(Polynomial_generic_dense, Polynomial_generic_field):
    def __init__(self, parent, x=None, check=True, is_gen = False, construct=False):
        Polynomial_generic_dense.__init__(self, parent, x, check, is_gen)

############################################################################
# XXX:  Ensures that the generic polynomials implemented in SAGE via PARI  #
# until at least until 4.5.0 unpickle correctly as polynomials implemented #
# via FLINT.                                                               #
from sage.structure.sage_object import register_unpickle_override
from sage.rings.polynomial.polynomial_rational_flint import Polynomial_rational_flint

register_unpickle_override( \
    'sage.rings.polynomial.polynomial_element_generic', \
    'Polynomial_rational_dense', Polynomial_rational_flint)
