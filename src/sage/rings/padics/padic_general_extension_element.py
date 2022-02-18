"""
Elements of general extensions of p-adic rings and fields; the base ring may also be an extension.

These are implemented as proxy elements, backed by an absolute extension.
"""
# ****************************************************************************
#       Copyright (C)      2019 David Roe <roed.math@gmail.com>
#                     2019-2022 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import deepcopy
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.ring_extension_element import RingExtensionWithBasisElement, RingExtensionElement
from sage.rings.ring_extension_conversion import backend_element, backend_parent
from sage.rings.infinity import infinity
from sage.structure.element import coerce_binop

class pAdicGeneralExtensionElement:
    def _front(self, b):
        if b.parent().is_field():
            P = self.parent().fraction_field()
        else:
            P = self.parent().integer_ring()
        return P.element_class(P, b)

    # Now a bunch of trivial functions
    def frobenius(self, arithmetic=True):
        return self._front(backend_element(self).frobenius(arithmetic=arithmetic))

    def __lshift__(self, shift):
        return self._front(backend_element(self) << shift)

    def __rshift__(self, shift):
        return self._front(backend_element(self) >> shift)

    def lift_to_precision(self, absprec=None):
        return self._front(backend_element(self).lift_to_precision(absprec=absprec))

    def expansion(self, n=None, lift_mode='simple', start_val=None):
        """
        Return a series expansion of this element.

        INPUT:

        - ``n`` -- integer (default: ``None``), if present, return only the
          n-th term of the expansion.

        - ``lift_mode`` -- string (default: ``"simple"``), can be one of
          ``"simple"``, ``"smallest"``, or ``"teichmuller"``.

        .. NOTE::

            Unlike the expansion methods of other element classes, we return
            elements of the residue field here, i.e., the output corresponds to
            the mode ``"…-residue"`` in terms of #33340.

        EXAMPLES:

        The entries of the expansion sum to the original element::

            sage: L.<a> = Qp(2).extension(x^2 + 2*x + 4)
            sage: series = a.expansion()
            sage: series
            2-adic expansion of ...
            sage: series = [L(s).lift_to_precision() for s in series]
            sage: sum([s<<(i + a.valuation()) for (i, s) in enumerate(series)]) == a
            True

        Terms can be selected explicitly, and these also sum to the original element::

            sage: series = [a.expansion(i) for i in range(a.precision_absolute() - a.valuation())]
            sage: series = [L(s).lift_to_precision() for s in series]
            sage: sum([s<<(i + a.valuation()) for (i, s) in enumerate(series)]) == a
            True

        """
        E = ExpansionIterable(self, backend_element(self).expansion(n=None, lift_mode=lift_mode, start_val=start_val), lift_mode)

        if n is None:
            return E

        return E[n]

    def _ext_p_list(self, pos):
        return backend_element(self)._ext_p_list(pos)

    def teichmuller_expansion(self, n=None):
        return self.expansion(n, lift_mode="teichmuller")

    def residue(self, absprec=1, field=None, check_prec=None):
        r = backend_element(self).residue(absprec=absprec, field=field, check_prec=check_prec)
        _, from_backend, _ = backend_parent(self.parent().residue_field(), map=True)
        return from_backend(r)

    def unit_part(self):
        return self._front(backend_element(self).unit_part())

    def slice(self, i, j, k=1, lift_mode="simple"):
        # This is weird
        return self._front(backend_element(self).slice(i, j, k=k, lift_mode=lift_mode))

    def __copy__(self):
        return self._front(backend_element(self))

    def __deepcopy__(self):
        return self._front(deepcopy(backend_element(self)))

    def _poly_rep(self):
        return self.polynomial().change_ring(self.parent()._FP_base())

    # This and _div_ could be moved to a FieldExtension class
    def __invert__(self):
        return self._front(~backend_element(self))

    def _div_(self, right):
        return self._front(backend_element(self) / backend_element(right))

    def inverse_of_unit(self):
        return self._front(backend_element(self).inverse_of_unit())

    def _floordiv_(self, right):
        return self._front(backend_element(self) // backend_element(right))

    def __pow__(self, right, dummy=None):
        return self._front(backend_element(self)**right)

    @coerce_binop
    def quo_rem(self, right):
        q, r = backend_element(self).quo_rem(backend_element(right))
        return self._front(q), self._front(r)

    def _mod_(self, right):
        return self._quo_rem(right)[1]

    def add_bigoh(self, absprec):
        return self._front(backend_element(self).add_bigoh(absprec))

    def _is_exact_zero(self):
        return backend_element(self)._is_exact_zero()

    def _is_inexact_zero(self):
        return backend_element(self)._is_inexact_zero()

    def _is_zero_rep(self):
        return backend_element(self)._is_zero_rep()

    def is_zero(self, absprec=None):
        return backend_element(self).is_zero(absprec)

    def __nonzero__(self):
        return backend_element(self).__nonzero__()

    def _richcmp_(self, other, op):
        return backend_element(self)._richcmp_(backend_element(other), op)

    def is_equal_to(self, right, absprec=None):
        right = self.parent().coerce(right)
        return backend_element(self).is_equal_to(backend_element(right), absprec=absprec)

    def is_integral(self):
        return self.valuation() >= 0

    def is_padic_unit(self):
        return self.valuation() == 0

    def is_unit(self):
        return backend_element(self).is_unit()

    def precision_absolute(self):
        return backend_element(self).precision_absolute()

    def precision_relative(self):
        return backend_element(self).precision_relative()

    def valuation(self, p=None):
        return backend_element(self).valuation(p=p)

    def val_unit(self):
        v, u = backend_element(self).val_unit()
        return v, self._front(u)

    def ordp(self, p=None):
        return self.valuation(p) / self.parent().absolute_e()

    def normalized_valuation(self):
        return self.ordp()

    def _min_valuation(self):
        return backend_element(self)._min_valuation()

    def add_bigoh(self, absprec):
        return self._front(backend_element(self).add_bigoh(absprec))

    def _cache_key(self):
        return backend_element(self)._cache_key()

    def __hash__(self):
        raise TypeError("unhashable type: 'sage.rings.padics.extension_element.pAdicGeneralExtensionElement'")

    def _repr_extension(self, **options):
        from .padic_printing import pAdicPrinter
        P = self.parent().fraction_field()
        printer = P._printer
        D = P._printer.dict()
        if D["mode"] in ["series", "digits", "bars"]:
            with pAdicPrinter(P._backend, D):
                return repr(backend_element(self))
        return printer.repr_gen(self._fieldify(), False)

    def _latex_extension(self, **options):
        from .padic_printing import pAdicPrinter
        P = self.parent().fraction_field()
        printer = P._printer
        D = P._printer.dict()
        if D["mode"] in ["series", "digits", "bars"]:
            with pAdicPrinter(P._backend, D):
                return backend_element(self)._latex_()
        return printer.repr_gen(self._fieldify(), True)

    def str(self):
        return repr(self)

    def additive_order(self, prec=None):
        return backend_element(self).additive_order(prec=prec)

    def artin_hasse_exp(self, prec=None, algorithm=None):
        return self._front(backend_element(self).artin_hasse_exp(prec, algorithm))

    def algdep(self, n):
        """
        Note that this is absolute
        """
        return backend_element(self).algdep(n)

    def algebraic_dependency(self, n):
        return self.algdep(n)

    @coerce_binop
    def gcd(self, other):
        return self._front(backend_element(self).gcd(backend_element(other)))

    @coerce_binop
    def xgcd(self, other):
        r, s, t = backend_element(self).xgcd(backend_element(other))
        return self._front(r), self._front(s), self._front(t)

    def euclidean_degree(self):
        return backend_element(self).euclidean_degree()

    def is_square(self):
        return backend_element(self).is_square()

    def is_squarefree(self):
        return backend_element(self).is_squarefree()

    def square_root(self, extend=True, all=False, algorithm=None):
        try:
            x = backend_element(self).square_root(extend=False, all=all, algorithm=algorithm)
        except ValueError:
            if extend:
                raise NotImplementedError("extending using the sqrt function not yet implemented")
            raise
        if all:
            return [self._front(r) for r in x]
        return self._front(x)

    def sqrt(self, extend=True, all=False, algorithm=None):
        return self.square_root(extend=extend, all=all, algorithm=algorithm)

    def nth_root(self, n, all=False):
        x = backend_element(self).nth_root(n, all=all)
        if all:
            return [self._front(r) for r in x]
        return self._front(x)

    def multiplicative_order(self, prec=None):
        return backend_element(self).multiplicative_order(prec)

    def is_prime(self):
        return backend_element(self).is_prime()

    def _rational_(self):
        return backend_element(self)._rational_()

    def log(self, p_branch=None, pi_branch=None, aprec=None, change_frac=False, algorithm=None):
        if p_branch is not None:
            p_branch = backend_element(p_branch)
        if pi_branch is not None:
            pi_branch = backend_element(pi_branch)
        return self._front(backend_element(self).log(
            p_branch=p_branch,
            pi_branch=pi_branch,
            aprec=aprec,
            change_frac=change_frac,
            algorithm=algorithm))

    def exp(self, aprec=None, algorithm=None):
        return self._front(backend_element(self).exp(aprec=aprec, algorithm=algorithm))

    def polylog(self, n, p_branch=0):
        return self._front(backend_element(self).polylog(n, p_branch=p_branch))

    def abs(self, prec=None):
        return backend_element(self).abs(prec)

    def __abs__(self):
        return self.abs()

    def _is_base_elt(self, p):
        return False

class pAdicGeneralRingExtensionElement(pAdicGeneralExtensionElement, RingExtensionElement):
    def _fieldify(self):
        x = backend_element(self)
        return self._front(x.parent().fraction_field()(x))

    def polynomial(self, var='x', base=None):
        return self._fieldify().polynomial(var, None if base is None else base.fraction_field()).change_ring(self.base_ring())

    def vector(self, base=None):
        if base is not None:
            base = base.fraction_field()
        return self._fieldify().vector(base)

    def _vector_(self, reverse=False, base=None):
        if base is not None:
            base = base.fraction_field()
        return self._fieldify()._vector_(reverse, base)

    def matrix(self, base=None):
        if base is not None:
            base = base.fraction_field()
        return self._fieldify().matrix(base)

    def trace(self, base=None):
        if base is not None:
            base = base.fraction_field()
        t = self._fieldify().trace(base)
        return t.parent().integer_ring()(t)

    def norm(self, base=None):
        if base is not None:
            base = base.fraction_field()
        n = self._fieldify().norm(base)
        return n.parent().integer_ring()(n)

    def charpoly(self, var='x', base=None):
        f = self.matrix(base).charpoly(var)
        return f.change_ring(f.base_ring().integer_ring())

    def minpoly(self, var='x', base=None):
        if base is not None:
            base = base.fraction_field()
        f = self._fieldify().minpoly(var, base)
        return f.change_ring(f.base_ring().integer_ring())

    def minimal_polynomial(self, var='x', base=None):
        return self.minpoly(var, base)

class pAdicGeneralFieldExtensionElement(pAdicGeneralExtensionElement, RingExtensionWithBasisElement):
    def _fieldify(self):
        return self

    def polynomial(self, var='x', base=None):
        return RingExtensionWithBasisElement.polynomial(self, var=var, base=base).univariate_polynomial()


# NotImplementedError: coercion/conversion maps
# NotImplementedError: printing


class ExpansionIterable(object):
    def __init__(self, element, backend_expansion, mode):
        self._element = element
        self._backend_expansion = backend_expansion
        self._mode = mode

    def __iter__(self):
        parent = self._element.parent().residue_field()
        for x in self._backend_expansion:
            yield parent(x)

    def __len__(self):
        return len(self._backend_expansion)

    def __getitem__(self, n):
        parent = self._element.parent().residue_field()
        return parent(self._backend_expansion[n])

    def __repr__(self):
        if self._mode == 'simple':
            modestr = ""
        elif self.mode == 'smallest':
            modestr = " (balanced)"
        else:
            modestr = " (teichmuller)"
        p = self._element.parent().prime()
        return f"{p}-adic expansion of {self._element}{modestr}"
