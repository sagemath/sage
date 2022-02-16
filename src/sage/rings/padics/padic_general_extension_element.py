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
from sage.rings.ring_extension_element import RingExtensionWithBasisElement
from sage.rings.ring_extension_conversion import backend_element, backend_parent
from sage.rings.infinity import infinity
from sage.structure.element import coerce_binop


class pAdicGeneralExtensionElement(RingExtensionWithBasisElement):
    def polynomial(self, var='x', base=None):
        return RingExtensionWithBasisElement.polynomial(self, var=var, base=base).univariate_polynomial()

    def _poly_rep(self):
        return self.polynomial().change_ring(self.parent()._FP_base())

    # Now a bunch of trivial functions
    def frobenius(self, arithmetic=True):
        return self.__class__(self.parent(), backend_element(self).frobenius(arithmetic=arithmetic))

    def __lshift__(self, shift):
        return self.__class__(self.parent(), backend_element(self) << shift)

    def __rshift__(self, shift):
        return self.__class__(self.parent(), backend_element(self) >> shift)

    def lift_to_precision(self, absprec=None):
        return self.__class__(self.parent(), backend_element(self).lift_to_precision(absprec=absprec))

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
        pass

    def unit_part(self):
        return self.__class__(self.parent(), backend_element(self).unit_part())

    def _is_base_elt(self, p):
        return self.parent().prime() == p and self.prime_pow.deg == 1

    def __copy__(self):
        return self.__class__(self.parent(), backend_element(self))

    def __deepcopy__(self):
        return self.__class__(self.parent(), deepcopy(backend_element(self)))

    # This and _div_ could be moved to a FieldExtension class
    def __invert__(self):
        return self.__class__(self.parent().fraction_field(), ~backend_element(self))

    def _div_(self, right):
        return self.__class__(self.parent().fraction_field(), backend_element(self) / backend_element(right))

    def _floordiv_(self, right):
        return self.__class__(self.parent(), backend_element(self) // backend_element(right))

    def __pow__(self, right, dummy=None):
        K = self.parent()
        if isinstance(right, (Integer, Rational, int)) and right < 0:
            K = K.fraction_field()
        return self.__class__(K, backend_element(self)**right)

    def _quo_rem(self, right):
        q, r = backend_element(self)._quo_rem(backend_element(right))
        C, R = self.__class__, self.parent()
        return C(R, q), C(R, r)

    def _mod_(self, right):
        return self._quo_rem(right)[1]

    def add_bigoh(self, absprec):
        return self.__class__(self.parent(), backend_element(self).add_bigoh(absprec))

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

    def precision_absolute(self):
        return backend_element(self).precision_absolute()

    def precision_relative(self):
        return backend_element(self).precision_relative()

    def valuation(self, p=None):
        return backend_element(self).valuation(p=p)

    def val_unit(self):
        v, u = backend_element(self).val_unit()
        return v, self.__class__(self.parent(), u)

    def ordp(self, p=None):
        return self.valuation(p) / self.parent().absolute_e()

    def _cache_key(self):
        return backend_element(self)._cache_key()

    def __hash__(self):
        raise TypeError("unhashable type: 'sage.rings.padics.extension_element.pAdicGeneralExtensionElement'")

    def _repr_extension(self, **options):
        from .padic_printing import pAdicPrinter
        P = self.parent()
        printer = P._printer
        D = P._printer.dict()
        if D["mode"] in ["series", "digits", "bars"]:
            with pAdicPrinter(P._backend, D):
                return repr(backend_element(self))
        return printer.repr_gen(self, False)

    def str(self):
        return repr(self)

    def additive_order(self, prec=None):
        return backend_element(self).additive_order(prec=prec)

    def artin_hasse_exp(self, prec=None, algorithm=None):
        return self.__class__(self.parent(), backend_element(self).artin_hasse_exp(prec, algorithm))

    def algdep(self, n):
        """
        Note that this is absolute
        """
        return backend_element(self).algdep(n)

    def algebraic_dependency(self, n):
        return self.algdep(n)

    @coerce_binop
    def gcd(self, other):
        return self.__class__(self.parent(), backend_element(self).gcd(backend_element(other)))

    @coerce_binop
    def xgcd(self, other):
        C = self.__class__
        P = self.parent()
        r, s, t = backend_element(self).xgcd(backend_element(other))
        return C(P, r), C(P, s), C(P, t)

    def is_square(self):
        return backend_element(self).is_square()

    def is_squarefree(self):
        return backend_element(self).is_squarefree()

    def multiplicative_order(self, prec=None):
        return backend_element(self).multiplicative_order(prec)

    def is_prime(self):
        return backend_element(self).is_prime()

    def _rational_(self):
        return backend_element(self)._rational_()

    def log(self, p_branch=None, pi_branch=None, aprec=None, change_frac=False, algorithm=None):
        P = self.parent()
        # FIXME: allow branch in fraction field
        if p_branch is not None:
            p_branch = backend_element(P(p_branch))
        if pi_branch is not None:
            pi_branch = backend_element(P(pi_branch))
        if change_frac:
            P = P.fraction_field()
        return self.__class__(P, backend_element(self).log(
            p_branch=p_branch,
            pi_branch=pi_branch,
            aprec=aprec,
            change_frac=change_frac,
            algorithm=algorithm))

    def exp(self, aprec=None, algorithm=None):
        ans = backend_element(self).exp(aprec=aprec, algorithm=algorithm)
        P = self.parent()
        if ans.parent().is_field():
            return self.__class__(P.fraction_field(), ans)
        else:
            return self.__class__(P, ans)

    def polylog(self, n, p_branch=0):
        ans = backend_element(self).polylog(n, p_branch=p_branch)
        P = self.parent()
        if ans.parent().is_field():
            return self.__class__(P.fraction_field(), ans)
        else:
            return self.__class__(P, ans)

    def abs(self, prec=None):
        return backend_element(self).abs(prec)

    def __abs__(self):
        return self.abs()

    def _is_base_elt(self, p):
        return False

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
