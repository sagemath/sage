r"""
Algebraic properties of hypergeometric functions.

[Tutorial]

AUTHORS:

- Xavier Caruso, Florian Fürnsinn (2025-10): initial version
"""

from sage.structure.sage_object import SageObject
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.sequence import Sequence

from sage.misc.misc_c import prod
from sage.misc.functional import log
from sage.functions.other import ceil
from sage.arith.functions import lcm
from sage.misc.cachefunc import cached_method

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.padics.factory import QpFP

from sage.categories.pushout import pushout
from sage.categories.finite_fields import FiniteFields

from sage.matrix.constructor import matrix

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing


def mod2(x, d):
    return d - (-x) % d

class Parameters():
    r"""
    Class for parameters of hypergeometric functions.
    """
    def __init__(self, alpha, beta, add_one=True):
        r"""
        Initialize this set of parameters.

        INPUT:

        - ``alpha`` -- list of top parameters

        - ``beta`` -- list of bottom parameters

        - ``add_one`` -- boolean (default: ``True``),
          if ``True``, add an additional one to the bottom
          parameters.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: p
            ([1/4, 1/3, 1/2], [2/5, 3/5, 1])
            sage: type(p)
            <class 'sage.functions.hypergeometric_algebraic.Parameters'>

        By default, parameters are sorted, duplicates are removed and
        a trailing `1` is added to the bottom parameters::

            sage: Parameters([1/2, 1/3, 2/3], [2/3])
            ([1/3, 1/2], [1])

        We can avoid adding the trailing `1` by passing ``add_one=False``)::

            sage: Parameters([1/2, 1/3, 2/3], [2/3], add_one=False)
            ([1/3, 1/2], [])
        """
        try:
            alpha = sorted([QQ(a) for a in alpha if a is not None])
            beta = sorted([QQ(b) for b in beta if b is not None])
        except TypeError:
            raise NotImplementedError("parameters must be rational numbers")
        i = j = 0
        while i < len(alpha) and j < len(beta):
            if alpha[i] == beta[j]:
                del alpha[i]
                del beta[j]
            elif alpha[i] > beta[j]:
                j += 1
            else:
                i += 1
        if add_one:
            beta.append(QQ(1))
        self.alpha = alpha
        self.beta = beta
        if len(alpha) == 0 and len(beta) == 0:
            self.d = 1
            self.bound = 1
        else:
            self.d =  lcm([ a.denominator() for a in alpha ]
                        + [ b.denominator() for b in beta ])
            self.bound = 2 * self.d * max(alpha + beta) + 1

    def __repr__(self):
        r"""
        Return a string representation of these parameters.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: p  # indirect doctest
            ([1/4, 1/3, 1/2], [2/5, 3/5, 1])
        """
        return "(%s, %s)" % (self.alpha, self.beta)

    def is_balanced(self):
        r"""
        Return ``True`` if there are as many top parameters as bottom
        parameters; ``False`` otherwise.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: p
            ([1/4, 1/3, 1/2], [2/5, 3/5, 1])
            sage: p.is_balanced()
            True

        ::

            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5], add_one=False)
            sage: p
            ([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: p.is_balanced()
            False
        """
        return len(self.alpha) == len(self.beta)

    @cached_method
    def christol_sorting(self, c=1):
        r"""
        """
        d = self.d
        def mod2(x):
            return d - (-x) % d
        A = [mod2(d*c*a), -a, 1) for a in self.alpha]
        B = [mod2(d*c*b), -b, -1) for b in self.beta]
        return sorted(A + B)

    def parenthesis_criterion(self, c):
        AB = self.christol_sorting(c)
        parenthesis = 0
        previous_paren = -1
        for _, _, paren in AB:
            parenthesis += paren
            if parenthesis < 0:
                return False
            previous_paren = paren
        return parenthesis >= 0

    def interlacing_criterion(self, c):
        AB = self.christol_sorting(c)
        previous_paren = -1
        for _, paren in AB:
            if paren == previous_paren:
                return False
            previous_paren = paren
        return True

    def remove_positive_integer_differences(self):
        differences = []
        alpha = self.alpha[:]
        beta = self.beta[:]
        for i in range(len(alpha)):
            for j in range(len(beta)):
                diff = alpha[i] - beta[j]
                if diff in ZZ and diff > 0:
                    differences.append((diff, i, j))
        for _, i, j in sorted(differences):
            if alpha[i] is not None and beta[j] is not None:
                alpha[i] = None
                beta[j] = None
        return Parameters(alpha, beta, add_one=False)

    def has_negative_integer_differences(self):
        return any(a - b in ZZ and a < b for a in self.alpha for b in self.beta)

    def shift(self):
        alpha = [a+1 for a in self.alpha]
        beta = [b+1 for b in self.beta]
        return Parameters(alpha, beta)

    def scalar(self):
	r"""
	Return the scalar of the derivative of the hypergeometric function with
	this set of parameters.

	EXAMPLES::

	    sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: p
            ([1/4, 1/3, 1/2], [2/5, 3/5, 1])
            sage: p.scalar()
            25/144
	"""
        return prod(self.alpha) / prod(self.beta)


class HypergeometricAlgebraic(SageObject, metaclass=ClasscallMetaclass):
    @staticmethod
    def __classcall_private__(cls, parameters, x, scalar):
        base = x.parent().base_ring()
        if base in FiniteFields() and base.is_prime_field():
            return HypergeometricAlgebraic_GFp(parameters, x, scalar)
        if base.characteristic() == 0:
            base = pushout(base, QQ)
            if base is QQ:
                return HypergeometricAlgebraic_QQ(parameters, x, scalar)
            return HypergeometricAlgebraic_charzero(parameters, x, scalar)
        else:
            raise NotImplementedError("the base must have characteristic zero or be a finite field")

    def __init__(self, parameters, x, scalar=None):
        self._base = base = x.base_ring()
        S = x.parent()
        if x != S.gen():
            raise NotImplementedError("x must be the variable")  # find a better message
        if scalar is None:
            self._scalar = base.one()
        else:
            self._scalar = base(scalar)
        if self._scalar == 0:
            self._parameters = None
        else:
            self._parameters = parameters
        self._S = S
        self._x = S.gen()
        self._variable_name = S.variable_name()

    def _repr_(self):
        s = "hypergeometric(%s, %s, %s)" % (self.alpha(), self.beta(), self._variable_name)
        if self._scalar != 1:
            s = str(self._scalar) + "*" + s
        return s

    def alpha(self):
        return tuple(self._parameters.alpha)

    def beta(self):
        return tuple(self._parameters.beta)[:-1]

    def denominator(self):
        return self._parameters.d

    def differential_operator(self, var='d'):
        S = PolynomialRing(self._base, self._variable_name)
        x = S.gen()
        D = OrePolynomialRing(S, S.derivation(), names=var)
        if self._scalar == 0:
            return D.one()
        t = x * D.gen()
        A = D.one()
        for a in self.alpha():
            A *= t + S(a)
        B = D.one()
        for b in self.beta():
            B *= t + S(b-1)
        L = t*B - x*A
        return D([ c//x for c in L.list() ])

    def derivative(self):
        parameters = self._parameters.shift()
        scalar = self._base(self._parameters.scalar()) * self._scalar
        return HypergeometricAlgebraic(parameters, self._x, scalar)


class HypergeometricAlgebraic_charzero(HypergeometricAlgebraic):
    def is_defined(self):
        return not any(b in ZZ and b < 0 for b in self._beta)

    def series(self, prec):
        S = PowerSeriesRing(self._base, name=self._variable_name)
        c = self._scalar
        coeffs = [c]
        for i in range(prec):
            for a in self.alpha():
                c *= a + i
            for b in self.beta():
                c /= b + i
            c /= i + 1
            coeffs.append(c)
        return S(coeffs, prec=prec)


class HypergeometricAlgebraic_QQ(HypergeometricAlgebraic_charzero):
    def mod(self, p):
        F = FiniteField(p)
        x = F['x'].gen()
        return HypergeometricAlgebraic(self._parameters, x, self._scalar)

    __mod__ = mod

    def has_good_reduction(self, p):
        h = self.reduce(p)
        return h.is_defined()

    def is_algebraic(self):
        if any(a in ZZ and a <= 0 for a in self.alpha()):
            return True
        if not self._parameters.is_balanced():
            return False
        simplified_parameters = self._parameters.remove_positive_integer_differences()
        if simplified_parameters.has_negative_integer_differences():
            return False
        d = self._parameters.d
        return all(simplified_parameters.interlacing_criterion(c)
                   for c in range(d) if d.gcd(c) == 1)

    def is_globally_bounded(self, include_infinity=True):
        if include_infinity and len(self.alpha()) > len(self.beta()) + 1:
            return False
        d = self.denominator()
        for c in range(d):
            if d.gcd(c) == 1:
                if not self._parameters.parenthesis_criterion(c):
                    return False
        return True


class HypergeometricAlgebraic_GFp(HypergeometricAlgebraic):
    def __init__(self, parameters, x, scalar=None):
        HypergeometricAlgebraic.__init__(self, parameters, x, scalar)
        self._p = self._base.cardinality()

    def is_defined(self):
        p = self._p
        d = self.denominator()
        if d.gcd(p) > 1:
            return False
        u = 1
        if not self._parameters.parenthesis_criterion(u):
            return False
        u = p % d
        while u != 1:
            if not self._parameters.parenthesis_criterion(u):
                return False
            u = p*u % d
        bound = self._parameters.bound
        if bound < p:
            return True
        prec = 1 + p ** ceil(log(self._parameters.bound, p))
        try:
            self.series(prec)
        except ValueError:
            return False
        return True

    def series(self, prec):
        S = PowerSeriesRing(self._base, name=self._variable_name)
        p = self._p
        pprec = max(1, (len(self.beta()) + 1) * ceil(log(prec, p)))
        K = QpFP(p, pprec)
        c = K(self._scalar)
        coeffs = [c]
        for i in range(prec-1):
            for a in self.alpha():
                c *= a + i
            for b in self.beta():
                c /= b + i
            c /= i + 1
            if c.valuation() < 0:
                raise ValueError("denominator appears in the series at the required precision")
            coeffs.append(c)
        return S(coeffs, prec=prec)

    def is_algebraic(self):
        return self.is_defined()

    def p_curvature(self):
        L = self.differential_operator()
        K = L.base_ring().fraction_field()
        S = OrePolynomialRing(K, L.parent().twisting_derivation().extend_to_fraction_field(), names='d')
        L = S(L.list())
        d = S.gen()
        p = self._p
        rows = [ ]
        n = L.degree()
        for i in range(p, p + n):
            Li = d**i % L
            rows.append([Li[j] for j in range(n)])
        return matrix(rows)

    def dwork_relation(self):
        r"""
        Return (P1, g1), ..., (Ps, gs) such that

            h = P1*g1^p + ... + Ps*gs^p
        """
        if not self._parameters.is_balanced():
            raise ValueError("the hypergeometric function must be a pFq with q = p-1")
        if not self.is_defined():
            raise ValueError("this hypergeometric function is not defined")

        # We compute the series expansion up to x^p
        p = self._p
        S = self._S
        cs = self.series(p).list()

        # We compute the relevant exponents
        exponents = sorted([(1-b) % p for b in self._parameters.beta])
        exponents.append(p)
        Ps = []
        for i in range(len(exponents) - 1):
            e = exponents[i]
            if e < len(cs) and cs[e]:
                P = S(cs[e:exponents[i+1]]) << e
                Ps.append((e, P))

        # We compute the hypergeometric series
        Hs = [ ]
        for r, P in Ps:
            alpha = [ ]
            for a in self.alpha():
                ap = (a + (-a) % p) / p  # Dwork map
                ar = prod(a + i for i in range(r))
                if ar % p == 0:
                    alpha.append(ap + 1)
                else:
                    alpha.append(ap)
            beta = [ ]
            for b in self.beta():
                bp = (b + (-b) % p) / p  # Dwork map
                br = prod(b + i for i in range(r))
                if br % p == 0:
                    beta.append(bp + 1)
                else:
                    beta.append(bp)
            parameters = Parameters(alpha, beta)
            Hs.append((P, HypergeometricAlgebraic_GFp(parameters, self._x)))

        return Hs
