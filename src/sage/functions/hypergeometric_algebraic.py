r"""
Algebraic properties of hypergeometric functions.

[Tutorial]

AUTHORS:

- Xavier Caruso, Florian Fürnsinn (2025-10): initial version
"""

# ***************************************************************************
#    Copyright (C) 2025 Xavier Caruso <xavier.caruso@normalesup.org>
#                       Florian Fürnsinn <florian.fuernsinn@univie.ac.at>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

import operator

from sage.misc.latex import latex
from sage.misc.latex import latex_variable_name
from sage.misc.cachefunc import cached_method

from sage.misc.misc_c import prod
from sage.misc.functional import log
from sage.functions.other import ceil
from sage.arith.misc import gcd
from sage.arith.functions import lcm
from sage.matrix.constructor import matrix

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.sequence import Sequence
from sage.structure.category_object import normalize_names

from sage.categories.action import Action
from sage.categories.pushout import pushout
from sage.categories.map import Map
from sage.categories.finite_fields import FiniteFields
from sage.sets.primes import Primes

from sage.matrix.special import companion_matrix
from sage.matrix.special import identity_matrix

from sage.symbolic.ring import SR
from sage.combinat.subset import Subsets
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.padics.factory import QpFP
from sage.rings.number_field.number_field import CyclotomicField

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing


def insert_zeroes(P, n):
    cs = P.list()
    coeffs = n * len(cs) * [0]
    for i in range(len(cs)):
        coeffs[n*i] = cs[i]
    return P.parent()(coeffs)

def kernel(M, repeat=2):
    n = M.nrows()
    m = M.ncols()
    if n > m + 1:
        raise RuntimeError
    if n <= m:
        K = M.base_ring().base_ring()
        for _ in range(repeat):
            a = K.random_element()
            Me = matrix(n, m, [f(a) for f in M.list()])
            if Me.rank() == n:
                return
    for J in Subsets(range(m), n-1):
        MJ = M.matrix_from_columns(J)
        minor = MJ.delete_rows([0]).determinant()
        if minor.is_zero():
            continue
        ker = [minor]
        for i in range(1, n):
            minor = MJ.delete_rows([i]).determinant()
            ker.append((-1)**i * minor)
        g = ker[0].leading_coefficient() * gcd(ker)
        ker = [c//g for c in ker]
        return ker


# Parameters of hypergeometric functions
########################################

class Parameters():
    r"""
    Class for parameters of hypergeometric functions.
    """
    def __init__(self, top, bottom, add_one=True):
        r"""
        Initialize this set of parameters.

        INPUT:

        - ``top`` -- list of top parameters

        - ``bottom`` -- list of bottom parameters

        - ``add_one`` -- boolean (default: ``True``),
          if ``True``, add an additional one to the bottom
          parameters.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: p
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: type(p)
            <class 'sage.functions.hypergeometric_algebraic.Parameters'>

        By default, parameters are sorted, duplicates are removed and
        a trailing `1` is added to the bottom parameters::

            sage: Parameters([1/2, 1/3, 2/3], [2/3])
            ((1/3, 1/2), (1,))

        We can avoid adding the trailing `1` by passing ``add_one=False``)::

            sage: Parameters([1/2, 1/3, 2/3], [2/3], add_one=False)
            ((1/3, 1/2, 1), (1,))
        """
        try:
            top = sorted([QQ(a) for a in top if a is not None])
            bottom = sorted([QQ(b) for b in bottom if b is not None])
        except TypeError:
            raise NotImplementedError("parameters must be rational numbers")
        i = j = 0
        while i < len(top) and j < len(bottom):
            if top[i] == bottom[j]:
                del top[i]
                del bottom[j]
            elif top[i] > bottom[j]:
                j += 1
            else:
                i += 1
        if add_one:
            bottom.append(QQ(1))
        else:
            try:
                i = bottom.index(QQ(1))
                bottom.append(QQ(1))
                del bottom[i]
            except ValueError:
                bottom.append(QQ(1))
                top.append(QQ(1))
                top.sort()
        self.top = tuple(top)
        self.bottom = tuple(bottom)
        if len(top) == 0 and len(bottom) == 0:
            self.d = 1
            self.bound = 1
        else:
            self.d =  lcm([ a.denominator() for a in top ]
                        + [ b.denominator() for b in bottom ])
            self.bound = 2 * self.d * max(top + bottom) + 1

    def __repr__(self):
        r"""
        Return a string representation of these parameters.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: p  # indirect doctest
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
        """
        return "(%s, %s)" % (self.top, self.bottom)

    def __hash__(self):
        return hash((self.top, self.bottom))

    def __eq__(self, other):
        return (isinstance(other, Parameters)
            and self.top == other.top and self.bottom == other.bottom)

    def is_balanced(self):
        r"""
        Return ``True`` if there are as many top parameters as bottom
        parameters; ``False`` otherwise.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: p
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: p.is_balanced()
            True

        ::

            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5], add_one=False)
            sage: p
            ((1/4, 1/3, 1/2, 1), (2/5, 3/5, 1))
            sage: p.is_balanced()
            False
        """
        return len(self.top) == len(self.bottom)

    @cached_method
    def christol_sorting(self, c=1):
        r"""
        """
        d = self.d
        def mod2(x):
            return d - (-x) % d
        A = [(mod2(d*c*a), -a, 1) for a in self.top]
        B = [(mod2(d*c*b), -b, -1) for b in self.bottom]
        return sorted(A + B)

    def parenthesis_criterion(self, c):
        parenthesis = 0
        previous_paren = -1
        for _, _, paren in self.christol_sorting(c):
            parenthesis += paren
            if parenthesis < 0:
                return False
            previous_paren = paren
        return parenthesis >= 0

    def q_christol_sorting(self, q):
        d = self.d
        A = [(1/2 + (-a) % q, 1) for a in self.top]
        B = [(1 + (-b) % q, -1) for b in self.bottom]
        return sorted(A + B)

    def q_parenthesis_criterion(self, q):
        parenthesis = 0
        previous_paren = -1
        for _, paren in self.q_christol_sorting(q):
            parenthesis += paren
            if parenthesis < 0:
                return False
            previous_paren = paren
        return parenthesis >= 0

    def q_interlacing_number(self, q):
        interlacing = 0
        previous_paren = -1
        for _, paren in self.q_christol_sorting(q):
            if paren == -1 and previous_paren == 1:
                interlacing += 1
            previous_paren = paren
        return interlacing

    def interlacing_criterion(self, c):
        r"""
        Return ``True`` if the sorted lists of the decimal parts (where integers
        are assigned 1 instead of 0) of c*a and c*b for a in the top parameters
        and b in the bottom parameters interlace, i.e., the entries in the sorted
        union of the two lists alternate between entries from the first and from
        the second list. Used to determine algebraicity of the hypergeometric
        function with these parameters with the Beukers-Heckman criterion.

        INPUT:

            - ``c`` -- an integer between 1 and ``self.d``, coprime to ``d``.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/3, 2/3], [1/2])
            sage: p
            ((1/3, 2/3), (1/2, 1))
            sage: p.interlacing_criterion(1)
            True
            sage: p.interlacing_criterion(5)
            True

        ::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/8, 3/8, 5/8], [1/4, 1/2])
            sage: p
            ((1/8, 3/8, 5/8), (1/4, 1/2, 1))
            sage: p.interlacing_criterion(1)
            True
            sage: p.interlacing_criterion(3)
            False
        """
        AB = self.christol_sorting(c)
        previous_paren = -1
        for _, _, paren in AB:
            if paren == previous_paren:
                return False
            previous_paren = paren
        return True

    def remove_positive_integer_differences(self):
        r"""
        Return parameters, where pairs consisting of a top parameter
        and a bottom parameter with positive integer differences are
        removed, starting with pairs of minimal positive integer
        difference.
        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([5/2, -1/2, 5/3], [3/2, 1/3])
            sage: p
            ((-1/2, 5/3, 5/2), (1/3, 3/2, 1))
            sage: p.remove_positive_integer_differences()
            ((-1/2, 5/3), (1/3, 1))

        The choice of which pair with integer differences to remove first
        is important::

            sage: p = Parameters([4, 2, 1/2], [1, 3])
            sage: p
            ((1/2, 2, 4), (1, 3, 1))
            sage: p.remove_positive_integer_differences()
            ((1/2,), (1,))
        """
        differences = []
        top = list(self.top)
        bottom = list(self.bottom)
        for i in range(len(top)):
            for j in range(len(bottom)):
                diff = top[i] - bottom[j]
                if diff in ZZ and diff > 0:
                    differences.append((diff, i, j))
        for _, i, j in sorted(differences):
            if top[i] is not None and bottom[j] is not None:
                top[i] = None
                bottom[j] = None
        return Parameters(top, bottom, add_one=False)

    def has_negative_integer_differences(self):
        r"""
        Return ``True`` if there exists a pair of a top parameter and a bottom
        parameter, such that the top one minus the bottom one is a negative integer;
        return ``False`` otherwise.

        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: p
            ((1/4, 1/3, 1/2), (2/5, 3/5, 1))
            sage: p.has_negative_integer_differences()
            False

        ::

            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/2])
            sage: p
            ((1/4, 1/3, 1/2), (2/5, 3/2, 1))
            sage: p.has_negative_integer_differences()
            True
        """
        return any(a - b in ZZ and a < b for a in self.top for b in self.bottom)

    def dwork_image(self, p):
        try:
            top = [(a + (-a) % p) / p for a in self.top]
            bottom = [(b + (-b) % p) / p for b in self.bottom]
        except ZeroDivisionError:
            raise ValueError("denominators of parameters are not coprime to p")
        return Parameters(top, bottom, add_one=False)

    def decimal_part(self):
        top = [1 + a - ceil(a) for a in self.top]
        bottom = [1 + b - ceil(b) for b in self.bottom]
        return Parameters(top, bottom, add_one=False)

    def frobenius_order(self, p):
        param = self.decimal_part()
        iter = param.dwork_image(p)
        i = 1
        while param != iter:
            iter = iter.dwork_image(p)
            i += 1
        return i


# Hypergeometric functions
##########################

class HypergeometricAlgebraic(Element):
    def __init__(self, parent, arg1, arg2=None, scalar=None):
        Element.__init__(self, parent)
        base = parent.base_ring()
        if scalar is None:
            self._scalar = base.one()
        else:
            self._scalar = base(scalar)
        if self._scalar == 0:
            self._parameters = None
        elif isinstance(arg1, HypergeometricAlgebraic):
            self._parameters = arg1._parameters
            self._scalar *= base(arg1._scalar)
        elif isinstance(arg1, Parameters):
            self._parameters = arg1
        else:
            self._parameters = Parameters(arg1, arg2)

    def __hash__(self):
        return hash((self.base_ring(), self._parameters, self._scalar))

    def __eq__(self, other):
        return (isinstance(other, HypergeometricAlgebraic)
            and self.base_ring() is other.base_ring()
            and self._parameters == other._parameters
            and self._scalar == other._scalar)

    def _repr_(self):
        if self._parameters is None:
            return "0"
        scalar = self._scalar
        if scalar == 1:
            s = ""
        elif scalar._is_atomic():
            scalar = str(scalar)
            if scalar == "-1":
                s = "-"
            else:
                s = scalar + "*"
        else:
            s = "(%s)*" % scalar
        s += "hypergeometric(%s, %s, %s)" % (self.top(), self.bottom(), self.parent().variable_name())
        return s

    def _latex_(self):
        if self._parameters is None:
            return "0"
        scalar = self._scalar
        if scalar == 1:
            s = ""
        elif scalar._is_atomic():
            scalar = latex(scalar)
            if scalar == "-1":
                s = "-"
            else:
                s = scalar
        else:
            s = r"\left(%s\right)" % scalar
        top = self.top()
        bottom = self.bottom()
        s += r"\,_{%s} F_{%s} " % (len(top), len(bottom))
        s += r"\left(\begin{matrix} "
        s += ",".join(latex(a) for a in top)
        s += r"\\"
        s += ",".join(latex(b) for b in bottom)
        s += r"\end{matrix}; %s \right)" % self.parent().latex_variable_name()
        return s

    def base_ring(self):
        return self.parent().base_ring()

    def top(self):
        return self._parameters.top

    def bottom(self):
        return self._parameters.bottom[:-1]

    def scalar(self):
        return self._scalar

    def _add_(self, other):
        if self._parameters is None:
            return other
        if isinstance(other, HypergeometricAlgebraic):
            if other._parameters is None:
                return self
            if self._parameters == other._parameters:
                scalar = self._scalar + other._scalar
                return self.parent()(self._parameters, scalar=scalar)
        return SR(self) + SR(other)

    def _neg_(self):
        if self._parameters is None:
            return self
        return self.parent()(self._parameters, scalar=-self._scalar)

    def _sub_(self, other):
        if self._parameters is None:
            return other
        if isinstance(other, HypergeometricAlgebraic):
            if other._parameters is None:
                return self
            if self._parameters == other._parameters:
                scalar = self._scalar - other._scalar
                return self.parent()(self._parameters, scalar=scalar)
        return SR(self) + SR(other)

    def _mul_(self, other):
        return SR(self) * SR(other)

    def _div_(self, other):
        return SR(self) / SR(other)

    def denominator(self):
        return self._parameters.d

    def differential_operator(self, var='d'):
        S = self.parent().polynomial_ring()
        x = S.gen()
        D = OrePolynomialRing(S, S.derivation(), names=var)
        if self._scalar == 0:
            return D.one()
        t = x * D.gen()
        A = D.one()
        for a in self._parameters.top:
            A *= t + S(a)
        B = D.one()
        for b in self._parameters.bottom:
            B *= t + S(b-1)
        L = B - x*A
        return D([c//x for c in L.list()])

    def derivative(self):
        top = [a+1 for a in self.top()]
        bottom = [b+1 for b in self.bottom()]
        scalar = prod(self._parameters.top) / prod(self._parameters.bottom)
        scalar = self.base_ring()(scalar) * self._scalar
        return self.parent()(top, bottom, scalar)


class HypergeometricAlgebraic_charzero(HypergeometricAlgebraic):
    def is_defined(self):
        return not any(b in ZZ and b < 0 for b in self._bottom)

    def series(self, prec):
        S = self.parent().power_series_ring()
        c = self._scalar
        coeffs = [c]
        for i in range(prec):
            for a in self._parameters.top:
                c *= a + i
            for b in self._parameters.bottom:
                c /= b + i
            coeffs.append(c)
        return S(coeffs, prec=prec)


class HypergeometricAlgebraic_QQ(HypergeometricAlgebraic_charzero):
    def mod(self, p):
        H = HypergeometricFunctions(FiniteField(p), self.parent().variable_name())
        return H(self._parameters, None, self._scalar)

    __mod__ = mod

    def has_good_reduction(self, p):
        h = self.reduce(p)
        return h.is_defined()

    def good_reduction_primes(self):
        r"""
        Return

        (modulus, congruence_classes, exceptionnal_primes)

        ALGORITHM:

        We rely on Christol's criterion ([CF2025]_)
        """
        params = self._parameters
        d = params.d

        # We check the parenthesis criterion for c=1
        if not params.parenthesis_criterion(1):
            return d, [], []

        # We check the parenthesis criterion for other c
        # and derive congruence classes with good reduction
        goods = {c: None for c in range(d) if d.gcd(c) == 1}
        goods[1] = True
        for c in goods:
            if goods[c] is not None:
                continue
            cc = c
            goods[c] = True
            while cc != 1:
                if goods[cc] is False or not params.parenthesis_criterion(cc):
                    goods[c] = False
                    break
                cc = (cc * c) % d
            if goods[c]:
                cc = c
                while cc != 1:
                    goods[cc] = True
                    cc = (cc * c) % d

        # We treat exceptionnal primes
        bound = params.bound
        exceptions = []
        for p in Primes():
            if p > bound:
                break
            if d % p == 0 or not goods[p % d]:
                continue
            q = p
            while q <= bound:
                if not self._parameters.q_parenthesis_criterion(q):
                    exceptions.append(p)
                    break
                q *= p

        goods = [c for c, v in goods.items() if v]
        return d, goods, exceptions

    def is_algebraic(self):
        if any(a in ZZ and a <= 0 for a in self.top()):
            return True
        if not self._parameters.is_balanced():
            return False
        simplified_parameters = self._parameters.remove_positive_integer_differences()
        if simplified_parameters.has_negative_integer_differences():
            return False
        d = simplified_parameters.d
        return all(simplified_parameters.interlacing_criterion(c)
                   for c in range(d) if d.gcd(c) == 1)

    def is_globally_bounded(self, include_infinity=True):
        if include_infinity and len(self.top()) > len(self.bottom()) + 1:
            return False
        d = self.denominator()
        for c in range(d):
            if d.gcd(c) == 1:
                if not self._parameters.parenthesis_criterion(c):
                    return False
        return True

    def monodromy(self, x=0, var='z'):
        params = self._parameters
        if not params.is_balanced():
            raise ValueError("hypergeometric equation is not Fuchsian")
        d = params.d
        K = CyclotomicField(d, names=var)
        z = K.gen()
        S = PolynomialRing(K, names='X')
        X = S.gen()
        if x == 0:
            B = prod(X - z**(b*d) for b in params.bottom)
            return companion_matrix(B, format='right').inverse()
        elif x == 1:
            A = prod(X - z**(a*d) for a in params.top)
            B = prod(X - z**(b*d) for b in params.bottom)
            return companion_matrix(A, format='right').inverse() * companion_matrix(B, format='right')
        elif x is infinity:
            A = prod(X - z**(a*d) for a in params.top)
            return companion_matrix(A, format='right')
        else:
            n = len(params.top)
            return identity_matrix(QQ, n)

    def is_maximum_unipotent_monodromy(self):
        # TODO: check this (maybe ask Daniel)
        return all(b in ZZ for b in self.bottom())

    is_mum = is_maximum_unipotent_monodromy


class HypergeometricAlgebraic_GFp(HypergeometricAlgebraic):
    def __init__(self, parent, arg1, arg2=None, scalar=None):
        HypergeometricAlgebraic.__init__(self, parent, arg1, arg2, scalar)
        self._p = self.base_ring().cardinality()

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

    def is_defined_conjectural(self):
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
        q = p
        while q <= bound:
            if not self._parameters.q_parenthesis_criterion(q):
                return False
            q *= p
        return True

    def series(self, prec):
        S = self.parent().power_series_ring()
        p = self._p
        pprec = max(1, (len(self.bottom()) + 1) * ceil(log(prec, p)))
        K = QpFP(p, pprec)
        c = K(self._scalar)
        coeffs = [c]
        for i in range(prec-1):
            for a in self._parameters.top:
                c *= a + i
            for b in self._parameters.bottom:
                c /= b + i
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

    def p_curvature_corank(self):
        return self._parameters.q_interlacing_number(self._p)

    def dwork_relation(self):
        r"""
        Return (P1, h1), ..., (Ps, hs) such that

            self = P1*h1^p + ... + Ps*hs^p
        """
        if not self._parameters.is_balanced():
            raise ValueError("the hypergeometric function must be a pFq with q = p-1")
        if not self.is_defined():
            raise ValueError("this hypergeometric function is not defined")

        H = self.parent()
        p = self._p

        # We compute the series expansion up to x^p
        cs = self.series(p).list()

        # We compute the relevant exponents and associated coefficients
        S = self.parent().polynomial_ring()
        exponents = sorted([(1-b) % p for b in self._parameters.bottom])
        exponents.append(p)
        Ps = []
        for i in range(len(exponents) - 1):
            ei = exponents[i]
            ej = exponents[i+1]
            P = S(cs[ei:ej])
            if P:
                Ps.append((ei, P << ei))

        # We compute the hypergeometric series
        pairs = [ ]
        for r, P in Ps:
            top = [ ]
            for a in self.top():
                ap = (a + (-a) % p) / p  # Dwork map
                ar = prod(a + i for i in range(r))
                if ar % p == 0:
                    top.append(ap + 1)
                else:
                    top.append(ap)
            bottom = [ ]
            for b in self.bottom():
                bp = (b + (-b) % p) / p  # Dwork map
                br = prod(b + i for i in range(r))
                if br % p == 0:
                    bottom.append(bp + 1)
                else:
                    bottom.append(bp)
            pairs.append((P, H(top, bottom)))

        return pairs

    def annihilating_ore_polynomial(self, var='Frob'):
        # QUESTION: does this method actually return the
        # minimal Ore polynomial annihilating self?
        if not self._parameters.is_balanced():
            raise NotImplementedError("the hypergeometric function is not a pFq with q = p-1")

        K = self.base_ring()
        p = self._p
        S = self.parent().polynomial_ring()
        zero = S.zero()
        Frob = S.frobenius_endomorphism()
        Ore = OrePolynomialRing(S, Frob, names=var)

        # We remove the scalar
        if self._scalar == 0:
            return Ore.one()
        self = self.parent()(self._parameters)

        order = self._parameters.frobenius_order(p)
        bound = self.p_curvature_corank()

        rows = [{self: S.one()}]
        # If row is the i-th item of rows, we have:
        #   self = sum_g row[g] * g**(p**i)
        q = 1
        while True:
            row = {}
            previous_row = rows[-1]
            for _ in range(order):
                row = {}
                for g, P in previous_row.items():
                    for Q, h in g.dwork_relation():
                        # here g = sum(Q * h^p)
                        if h in row:
                            row[h] += P * insert_zeroes(Q, q)
                        else:
                            row[h] = P * insert_zeroes(Q, q)
                previous_row = row
                q *= p  # q = p**i
            rows.append(row)

            i = len(rows)
            Mrows = []
            Mqo = 1
            columns = {}
            for j in range(i-1, max(-1, i-2-bound), -1):
                for col in rows[j]:
                    columns[col] = None
            for j in range(i-1, max(-1, i-2-bound), -1):
                Mrow = []
                for col in columns:
                    Mrow.append(insert_zeroes(rows[j].get(col, zero), Mqo))
                Mrows.append(Mrow)
                Mqo *= p ** order
            M = matrix(S, Mrows)

            ker = kernel(M)
            if ker is not None:
                return insert_zeroes(Ore(ker), order)

    def is_lucas(self):
        p = self._p
        if self._parameters.frobenius_order(p) > 1:
            # TODO: check this
            return False
        S = self.parent().polynomial_ring()
        K = S.fraction_field()
        Ore = OrePolynomialRing(K, K.frobenius_endomorphism(), names='F')
        Z = Ore(self.annihilating_ore_polynomial())
        Ap = self.series(p).polynomial()
        F = Ap * Ore.gen() - 1
        return (Z % F).is_zero()


# Parent
########

class HypergeometricToSR(Map):
    def _call_(self, h):
        from sage.functions.hypergeometric import _hypergeometric
        return h.scalar() * _hypergeometric(h.top(), h.bottom(), SR.var(h.parent().variable_name()))

class ScalarMultiplication(Action):
    def _act_(self, scalar, h):
        return h.parent()(h, scalar=scalar)


class HypergeometricFunctions(Parent, UniqueRepresentation):
    def __init__(self, base, name, category=None):
        self._name = normalize_names(1, name)[0]
        self._latex_name = latex_variable_name(self._name)
        char = base.characteristic()
        if base in FiniteFields() and base.is_prime_field():
            self.Element = HypergeometricAlgebraic_GFp
        elif char == 0:
            base = pushout(base, QQ)
            if base is QQ:
                self.Element = HypergeometricAlgebraic_QQ
            else:
                self.Element = HypergeometricAlgebraic_charzero
        else:
            raise NotImplementedError("hypergeometric functions are only implemented over finite field and bases of characteristic zero")
        Parent.__init__(self, base, category=category)
        self.register_action(ScalarMultiplication(base, self, False, operator.mul))
        self.register_action(ScalarMultiplication(base, self, True, operator.mul))
        if char == 0:
            SR.register_coercion(HypergeometricToSR(self.Hom(SR)))

    def _repr_(self):
        return "Hypergeometric functions in %s over %s" % (self._name, self._base)

    def _element_constructor_(self, *args, **kwds):
        return self.element_class(self, *args, **kwds)

    def _coerce_map_from_(self, other):
        if (isinstance(other, HypergeometricFunctions)
        and other.has_coerce_map_from(self)):
            return True

    def _pushout_(self, other):
        if isinstance(other, HypergeometricFunctions) and self._name == other._name:
            base = pushout(self.base_ring(), other.base_ring())
            if base is not None:
                return HypergeometricFunctions(base, self._name)
        if SR.has_coerce_map_from(other):
            return SR

    def base_ring(self):
        return self._base

    def variable_name(self):
        return self._name

    def latex_variable_name(self):
        return self._latex_name

    def polynomial_ring(self):
        return PolynomialRing(self.base_ring(), self._name)

    def power_series_ring(self, default_prec=None):
        return PowerSeriesRing(self.base_ring(), self._name, default_prec=default_prec)
