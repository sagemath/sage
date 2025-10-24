r"""
Algebraic properties of hypergeometric functions.

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

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing


def mod2(x, d):
    return d - (-x) % d


class Parameters():
    def __init__(self, alpha, beta, is_one_included=False):
        try:
            alpha = sorted([QQ(a) for a in alpha])
            beta = sorted([QQ(b) for b in beta])
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
        if not is_one_included:
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

    def is_balanced(self):
        return len(self.alpha) == len(self.beta)

    @cached_method
    def parenthesis_criterion(self, c):
        d = self.d
        A = [(mod2(d*c*a, d), -a, 1) for a in self.alpha]
        B = [(mod2(d*c*b, d), -b, -1) for b in self.beta]
        AB = sorted(A + B)
        parenthesis = 0
        previous_paren = -1
        for _, _, paren in AB:
            parenthesis += paren
            if parenthesis < 0:
                return False
            previous_paren = paren
        return parenthesis >= 0

    @cached_method
    def interlacing_criterion(self, c):
        d = self.d
        A = [(mod2(d*c*a, d), 1) for a in self.alpha]
        B = [(mod2(d*c*b, d), -1) for b in self.beta]
        AB = sorted(A + B)
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
        alpha = [a for a in alpha if a is not None]
        beta = [b for b in beta if b is not None]
        return Parameters(alpha, beta, is_one_included=True)

    def has_negative_integer_differences(self):
        return any(a - b in ZZ and a < b for a in self.alpha for b in self.beta)

    def shift(self):
        alpha = [a+1 for a in self.alpha]
        beta = [b+1 for b in self.beta]
        return Parameters(alpha, beta)

    def scalar(self):
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
            raise NotImplementedError

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
    def reduce(self, p):
        F = FiniteField(p)
        x = F['x'].gen()
        return HypergeometricAlgebraic(self._parameters, x, self._scalar)

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
    def __init__(self, alpha, beta, x, scalar=None):
        HypergeometricAlgebraic.__init__(self, alpha, beta, x, scalar)
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
        prec = 1 + p ** ceil(log(self._bound(), p))
        try:
            self.series(prec)
        except ValueError:
            return False
        return True

    def series(self, prec):
        S = PowerSeriesRing(self._base, name=self._variable_name)
        p = self._p
        pprec = max(1, (len(self._beta) + 1) * ceil(log(prec, p)))
        K = QpFP(p, pprec)
        c = K(self._scalar)
        coeffs = [c]
        for i in range(prec-1):
            for a in self.alpha():
                c *= a + i
            for b in self._beta():
                c /= b + i
            c /= i + 1
            if c.valuation() < 0:
                raise ValueError("denominator appears in the series at the required precision")
            coeffs.append(c)
        return S(coeffs, prec=prec)

    def is_algebraic(self):
        return self.is_defined()
