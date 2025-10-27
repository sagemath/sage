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

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.sequence import Sequence
from sage.structure.category_object import normalize_names

from sage.categories.pushout import pushout
from sage.categories.map import Map
from sage.categories.finite_fields import FiniteFields

from sage.symbolic.ring import SR
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.finite_field_constructor import FiniteField
from sage.rings.padics.factory import QpFP

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing

from sage.misc.misc_c import prod
from sage.misc.functional import log
from sage.functions.other import ceil
from sage.arith.functions import lcm
from sage.misc.cachefunc import cached_method

from sage.matrix.constructor import matrix

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
        self.top = top
        self.bottom = bottom
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
            ([1/4, 1/3, 1/2], [2/5, 3/5, 1])
        """
        return "(%s, %s)" % (self.top, self.bottom)

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
        for _, _, paren in AB:
            if paren == previous_paren:
                return False
            previous_paren = paren
        return True

    def remove_positive_integer_differences(self):
        differences = []
        top = self.top[:]
        bottom = self.bottom[:]
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
        Returns ``True`` if there exists a pair of a top parameter and a bottom
        parameter, such that the top one minus the bottom one is a negative integer.
        Returns ``False`` otherwise.
        EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: p
            ([1/4, 1/3, 1/2], [2/5, 3/5, 1])
            sage: p.has_negative_integer_differences()
            False

        ::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/2])
            sage: p
            ([1/4, 1/3, 1/2], [2/5, 3/2, 1])
            sage: p.has_negative_integer_differences()
            True
        """
        return any(a - b in ZZ and a < b for a in self.top for b in self.bottom)

    def shift(self):
        r"""
        Return the parameters obtained by adding one to each of them. Is used to 
        define the derivative of the hypergeometric function with these parameters.
	EXAMPLES::

            sage: from sage.functions.hypergeometric_algebraic import Parameters
            sage: p = Parameters([1/4, 1/3, 1/2], [2/5, 3/5])
            sage: p
            ([1/4, 1/3, 1/2], [2/5, 3/5, 1])
            sage: p.shift()
            ([5/4, 4/3, 3/2], [7/5, 8/5, 1])
        """
        top = [a+1 for a in self.top]
        bottom = [b+1 for b in self.bottom]
        return Parameters(top, bottom)

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
        return prod(self.top) / prod(self.bottom)


# Hypergeometric functions
##########################

class HypergeometricAlgebraic(Element):
    def __init__(self, parent, parameters, scalar=None):
        Element.__init__(self, parent)
        base = parent.base_ring()
        if scalar is None:
            self._scalar = base.one()
        else:
            self._scalar = base(scalar)
        if self._scalar == 0:
            self._parameters = None
        else:
            self._parameters = parameters

    def _repr_(self):
        s = "hypergeometric(%s, %s, %s)" % (self.top(), self.bottom(), self.parent().variable_name())
        if self._scalar != 1:
            s = str(self._scalar) + "*" + s
        return s

    def base_ring(self):
        return self.parent().base_ring()

    def top(self):
        return tuple(self._parameters.top)

    def bottom(self):
        return tuple(self._parameters.bottom[:-1])

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
        for a in self.top():
            A *= t + S(a)
        B = D.one()
        for b in self.bottom():
            B *= t + S(b-1)
        L = t*B - x*A
        return D([ c//x for c in L.list() ])

    def derivative(self):
        parameters = Parameters(self.top(), self.bottom(), add_one=False).shift()
        scalar = self.base_ring()(self._parameters.scalar()) * self._scalar
        H = HypergeometricFunctions(self.base_ring(), self.parent().variable_name())
        return H(parameters, scalar=scalar) #in the output the tuple of bottom parameters has an extra comma, 


class HypergeometricAlgebraic_charzero(HypergeometricAlgebraic):
    def is_defined(self):
        return not any(b in ZZ and b < 0 for b in self._bottom)

    def series(self, prec):
        S = PowerSeriesRing(self._base, name=self._variable_name)
        c = self._scalar
        coeffs = [c]
        for i in range(prec):
            for a in self.top():
                c *= a + i
            for b in self.bottom():
                c /= b + i
            c /= i + 1
            coeffs.append(c)
        return S(coeffs, prec=prec)


class HypergeometricAlgebraic_QQ(HypergeometricAlgebraic_charzero):
    def mod(self, p):
        H = HypergeometricFunctions(FiniteField(p), self.parent().variable_name())
        return H(self._parameters, self._scalar)

    __mod__ = mod

    def has_good_reduction(self, p):
        h = self.reduce(p)
        return h.is_defined()

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


class HypergeometricAlgebraic_GFp(HypergeometricAlgebraic):
    def __init__(self, parent, parameters, scalar=None):
        HypergeometricAlgebraic.__init__(self, parent, parameters, scalar)
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

    def series(self, prec):
        S = self.parent().polynomial_ring()
        p = self._p
        pprec = max(1, (len(self.bottom()) + 1) * ceil(log(prec, p)))
        K = QpFP(p, pprec)
        c = K(self._scalar)
        coeffs = [c]
        for i in range(prec-1):
            for a in self.top():
                c *= a + i
            for b in self.bottom():
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
        S = self.parent().polynomial_ring()
        cs = self.series(p).list()

        # We compute the relevant exponents
        exponents = sorted([(1-b) % p for b in self._parameters.bottom])
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
            parameters = Parameters(top, bottom)
            Hs.append((P, HypergeometricAlgebraic_GFp(parameters, self._x)))

        return Hs


# Parent
########

class HypergeometricToSR(Map):
    def _call_(self, h):
        from sage.functions.hypergeometric import _hypergeometric
        return _hypergeometric(h.top(), h.bottom(), SR.var(h.parent().variable_name()))

class HypergeometricFunctions(Parent, UniqueRepresentation):
    def __init__(self, base, name, category=None):
        self._name = normalize_names(1, name)[0]
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
        if char == 0:
            SR.register_coercion(HypergeometricToSR(self.Hom(SR)))

    def _repr_(self):
        return "Hypergeometric functions in %s over %s" % (self._name, self._base)

    def _element_constructor_(self, *args, **kwds):
        return self.element_class(self, *args, **kwds)

    def base_ring(self):
        return self._base

    def variable_name(self):
        return self._name

    def polynomial_ring(self):
        return PolynomialRing(self.base_ring(), self._name)

    def power_series_ring(self, default_prec=None):
        return PowerSeriesRing(self.base_ring(), self._name, default_prec=prec)
