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
from sage.rings.ring_extension_element import RingExtensionElement
from .padic_generic_element import pAdicGenericElement


class pAdicGeneralExtensionElement(RingExtensionElement, pAdicGenericElement):
    def __init__(self, parent, value):
        RingExtensionElement.__init__(self, parent, value)
        pAdicGenericElement.__init__(self, parent)

    # We start with the interesting functions; need to port these to two step extensions
    def polynomial(self, var='x'):
        raise NotImplementedError

    def _poly_rep(self):
        return self.polynomial().change_ring(self.parent()._FP_base())

    # Now a bunch of trivial functions
    def frobenius(self, arithmetic=True):
        return self.__class__(self.parent(), self._element.frobenius(arithmetic=arithmetic))

    def __lshift__(self, shift):
        return self.__class__(self.parent(), self._element << shift)

    def __rshift__(self, shift):
        return self.__class__(self.parent(), self._element >> shift)

    def lift_to_precision(self, absprec=None):
        return self.__class__(self.parent(), self._element.lift_to_precision(absprec=absprec))

    def expansion(self, n=None, lift_mode='simple', start_val=None):
        """
        Write docstring clarifying that this is in terms of the absolute extension
        """
        def wrap(x):
            if lift_mode == 'teichmuller':
                x = self.__class__(self.parent(), x)
            return x

        E = ExpansionIterable(self, self._element.expansion(n, lift_mode, start_val), wrap, lift_mode)
        if n is None:
            return E
        elif self.parent().is_field():
            v = self.valuation()
            return wrap(E[n-v])
        else:
            return wrap(E[n])

    def _ext_p_list(self, pos):
        return self._element._ext_p_list(pos)

    def unit_part(self):
        return self.__class__(self.parent(), self._element.unit_part())

    def _is_base_elt(self, p):
        return self.parent().prime() == p and self.prime_pow.deg == 1

    def residue(self, absprec=1, field=None, check_prec=True):
        raise NotImplementedError

    def __copy__(self):
        return self.__class__(self.parent(), self._element)

    def __deepcopy__(self):
        return self.__class__(self.parent(), deepcopy(self._element))

    # This and _div_ could be moved to a FieldExtension class
    def __invert__(self):
        return self.__class__(self.parent().fraction_field(), ~self._element)

    def _div_(self, right):
        return self.__class__(self.parent().fraction_field(), self._element / right._element)

    def __pow__(self, right, dummy):
        K = self.parent()
        if isinstance(right, (Integer, Rational, int)) and right < 0:
            K = K.fraction_field()
        return self.__class__(K, self._element**right)

    def _quo_rem(self, right):
        q, r = self._element._quo_rem(right._element)
        K = self.parent()
        C = self.__class__
        return C(K, q), C(K, r)

    def add_bigoh(self, absprec):
        return self.__class__(self.parent(), self._element.add_bigoh(absprec))

    def _is_exact_zero(self):
        return self._element._is_exact_zero()

    def _is_inexact_zero(self):
        return self._element._is_inexact_zero()

    def _is_zero_rep(self):
        return self._element._is_zero_rep()

    def is_zero(self, absprec=None):
        return self._element.is_zero(absprec)

    def __nonzero__(self):
        return self._element.__nonzero__()

    def is_equal_to(self, right, absprec=None):
        right = self.parent().coerce(right)
        return self._element.is_equal_to(right._element, absprec=absprec)

    def precision_absolute(self):
        return self._element.precision_absolute()

    def precision_relative(self):
        return self._element.precision_relative()

    def valuation(self):
        return self._element.valuation()

    def val_unit(self):
        v, u = self._element.val_unit()
        return v, self.__class__(self.parent(), u)

    def _cache_key(self):
        return self._element._cache_key()

    def __hash__(self):
        raise TypeError("unhashable type: 'sage.rings.padics.extension_element.pAdicGeneralExtensionElement'")

# NotImplementedError: coercion/conversion maps
# NotImplementedError: printing


class ExpansionIterable(object):
    def __init__(self, elt, exp, wrap, mode):
        self._elt = elt
        self._exp = exp
        self._wrap = wrap
        self._mode = mode

    def __iter__(self):
        for x in self._exp:
            yield self._wrap(x)

    def __len__(self):
        return len(self._exp)

    def __getitem__(self, n):
        return self._wrap(self._exp[n])

    def __repr__(self):
        if self._mode == 'simple':
            modestr = ""
        elif self.mode == 'smallest':
            modestr = " (balanced)"
        else:
            modestr = " (teichmuller)"
        p = self._elt.prime_pow.prime
        return f"{p}-adic expansion of {self._elt}{modestr}"
