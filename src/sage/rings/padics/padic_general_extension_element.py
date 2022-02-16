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
from sage.rings.ring_extension_conversion import backend_element
from .padic_generic_element import pAdicGenericElement
from sage.rings.infinity import infinity


class pAdicGeneralExtensionElement(RingExtensionWithBasisElement, pAdicGenericElement):
    def __init__(self, parent, value, absprec=infinity, relprec=infinity):
        RingExtensionWithBasisElement.__init__(self, parent, value, absprec=absprec, relprec=relprec)
        pAdicGenericElement.__init__(self, parent)

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

    def __pow__(self, right, dummy=None):
        K = self.parent()
        if isinstance(right, (Integer, Rational, int)) and right < 0:
            K = K.fraction_field()
        return self.__class__(K, backend_element(self)**right)

    def _quo_rem(self, right):
        q, r = backend_element(self)._quo_rem(backend_element(right))
        K = self.parent()
        C = self.__class__
        return C(K, q), C(K, r)

    def add_bigoh(self, absprec):
        return self.__class__(self.parent(), backend_element(self).add_bigoh(absprec))

    def _is_exact_zero(self):
        return backend_element(self)._is_exact_zero()

    def _is_inexact_zero(self):
        return backend_element(self)._is_inexact_zero()

    def _is_zero_rep(self):
        return self._backend._is_zero_rep()

    def is_zero(self, absprec=None):
        return backend_element(self).is_zero(absprec)

    def __nonzero__(self):
        return backend_element(self).__nonzero__()

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
